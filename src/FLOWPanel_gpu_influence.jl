#=##############################################################################
# DESCRIPTION
    Task 051 stage 2: brute-force "rectangular" influence seam.

    When armed (FLOWPANEL_GPU_INFLUENCE env var or `set_gpu_influence!`), the
    per-step cross-influence passes of `_steady_aerodynamics!` — pass 1
    (wake sources -> bodies + wake probes + particles) and pass 3 (bodies ->
    bodies + wake probes + particles) — are evaluated by
    `FastMultipole.direct_rectangular!` (host threaded, or CUDA when
    available) instead of `FastMultipole.fmm!`. The panel SOLVE pass and any
    pass this seam does not recognize fall through to the existing fmm! path
    unchanged. Nothing changes unless armed.

    The seam lives in `influence!(::Tuple, ::Tuple, ::FastMultipoleBackend)`
    (src/FLOWPanel_fmm.jl), after the precalc hook: `_gpu_rect_influence!`
    returns `true` iff it fully handled the pass.

    Pass recognition (see FLOWPanel_simulate.jl `_sa_wake_influence!` /
    `_sa_body_influence!`):
      - pass 1: any source is a wake system (FLOWVPM.ParticleField, PanelWake,
        FilamentWrapper); requested outputs are velocity (+ per-target
        gradient), never the scalar potential.
      - pass 3: all sources are AbstractBody AND the caller passed
        `direct_conditioning`. NOTE (051 audit): `_sa_body_influence!` is the
        per-step caller, but other `direct_conditioning`-tagged body-source
        calls also match (FLOWPanel_kutta.jl:960, FLOWPanel_wake.jl:64,
        FLOWPanel_replay.jl:1053) — the seam's capture surface is
        "conditioned body-source influence", not just pass 3.
    The solve pass requests `scalar_potential=true` (Dirichlet) or velocity
    without `direct_conditioning` (Neumann/Krylov operators) and never
    matches.

    Fallback conditions (returns false; fmm! runs as before):
      - not armed, or the request shape does not match a recognized pass;
      - nonzero `extra_outputs` (induced-vorticity accumulation);
      - a FilamentWrapper source NOT wrapping a PanelWake (a PanelWake-wrapped
        FilamentWrapper's active final-row filaments are packed as tag-3 nv=2
        open-segment columns — 051 FilamentWrapper seam extension);
      - a ParticleField source whose kernel is not gaussianerf;
      - a body element set outside {ConstantSource, ConstantDoublet,
        VortexRing, ConstantSource+VortexRing, ConstantSource+ConstantDoublet}
        or a RigidWakeBody with a semi-infinite attached wake;
      - a self-paired HOST ParticleField with SFS enabled (the FMM post hook
        needs fmm! tree outputs for Estr; Estr stays on the FLOWVPM radix
        path).

    Task 052 extensions (device-resident wake):
      - a CuArray-backed ParticleField target evaluates its SELF influence
        (U, J, and E_str when SFS is enabled) through the FLOWVPM radix GPU
        lifecycle (`UJ_fmm_gpu!`), and receives all other sources through
        device-resident `direct_rectangular!` (no per-step particle H2D/D2H);
        its SFS AfterUJ hook runs via `post_evaluate_influence!(...,
        ::Nothing)` in src/FLOWPanel_gpu_wake.jl.
      - a CuArray-backed ParticleField SOURCE is packed as a device 7 x np
        slice (rows 1:7 == X, Gamma, sigma — the RectangularGaussianErfVortex
        layout).
      - probe passes (all targets FastMultipole.ProbeSystem, request
        `gradient=true` == velocity, `hessian` == velocity gradient) are
        recognized so monitor probe evaluations (BoundCirculationMonitor
        Kelvin slices, KuttaJoukowskiForce) never route a device particle
        field into host fmm!.
      - env-gated pass timers (FLOWPANEL_GPU_TIMERS=1).

    Known deviation from the CPU fmm! path: `direct_conditioning`
    (`_self_panel_core_size_conditioning`) flips the source core_size to
    `core_size_panel` only for NEARFIELD self-pair blocks, while farfield
    self pairs use unregularized expansions. Here ALL self-body pairs are
    direct, so the whole self-body influence is evaluated at
    `core_size_panel`. Beyond `radius_inflation` the regularized and
    singular kernels agree to the fmm tolerance, so the difference is within
    the fmm approximation budget.

# AUTHORSHIP
  * Created by  : task 051 stage 2 (agent), Aug 2026
  * License     : GNU Public License
=###############################################################################

################################################################################
# ARMING SWITCH
################################################################################

"Seam state: :unset (read env on first query), :off, :host, :cuda."
const GPU_INFLUENCE = Ref{Symbol}(:unset)

"Diagnostic: number of influence! calls the seam has fully handled."
const GPU_INFLUENCE_HITS = Ref{Int}(0)

"""
    gpu_influence_enabled()

Whether the rectangular direct-influence seam is armed. Defaults to off; armed
by the `FLOWPANEL_GPU_INFLUENCE` environment variable (`1`/`host` = host
threaded, `cuda`/`gpu` = CUDA with host fallback when CUDA is not functional)
or [`set_gpu_influence!`](@ref). Read lazily (not at precompile time).
"""
function gpu_influence_enabled()
    if GPU_INFLUENCE[] === :unset
        val = lowercase(get(ENV, "FLOWPANEL_GPU_INFLUENCE", "0"))
        GPU_INFLUENCE[] = val in ("1", "host", "true") ? :host :
                          val in ("cuda", "gpu") ? :cuda : :off
    end
    return GPU_INFLUENCE[] !== :off
end

"Set the seam mode programmatically (`:off`, `:host`, `:cuda`)."
function set_gpu_influence!(mode::Symbol)
    mode in (:off, :host, :cuda) || throw(ArgumentError(
        "gpu influence mode must be :off, :host, or :cuda (got $(repr(mode)))"))
    GPU_INFLUENCE[] = mode
    return mode
end

# lazy CUDA readiness: FastMultipole owns the CUDA dependency (lazy include of
# translate_batched_cuda.jl, which installs the CuMatrix direct_rectangular!
# methods). nothing = untested.
const _GPU_CUDA_READY = Ref{Union{Nothing,Bool}}(nothing)

function _gpu_cuda_ready()
    if _GPU_CUDA_READY[] === nothing
        ok = try
            FastMultipole.load_cuda_radix_lifecycle!() &&
                isdefined(FastMultipole, :CUDA)
        catch
            false
        end
        _GPU_CUDA_READY[] = ok
        ok || @info "FLOWPanel gpu influence: CUDA not functional " *
            "($(FastMultipole.cuda_radix_status())); using host rectangular path"
    end
    return _GPU_CUDA_READY[]::Bool
end

_gpu_device() = GPU_INFLUENCE[] === :cuda && _gpu_cuda_ready() ? :cuda : :host

# one-shot informational messages so production logs are not flooded
const _GPU_INFO_ONCE = Set{Symbol}()
function _gpu_info_once(key::Symbol, msg::AbstractString)
    key in _GPU_INFO_ONCE && return nothing
    push!(_GPU_INFO_ONCE, key)
    @info msg
    return nothing
end

################################################################################
# WORK BUFFERS (preallocated, reused across steps)
################################################################################

# keyed by (system identity, role); reallocated only when the required size
# changes (particle counts grow as the wake sheds)
const _GPU_RECT_WORK = Dict{Tuple{UInt,Symbol},Matrix{Float64}}()

function _gpu_workmat!(sys, role::Symbol, nrows::Int, ncols::Int)
    key = (objectid(sys), role)
    m = get(_GPU_RECT_WORK, key, nothing)
    if m === nothing || size(m) != (nrows, ncols)
        m = Matrix{Float64}(undef, nrows, ncols)
        _GPU_RECT_WORK[key] = m
    end
    return m::Matrix{Float64}
end

"Drop all cached seam work buffers (e.g. after geometry rebuilds)."
clear_gpu_influence_buffers!() = (empty!(_GPU_RECT_WORK); nothing)

################################################################################
# SOURCE PACKING
################################################################################

# element-type -> (tag, n strengths) per the RectangularPanelInfluence row
# layout (FastMultipole src/direct_rectangular.jl)
_gpu_panel_tag(::Type{ConstantSource}) = (1, 1)
_gpu_panel_tag(::Type{ConstantDoublet}) = (2, 1)
_gpu_panel_tag(::Type{VortexRing}) = (3, 1)
_gpu_panel_tag(::Type{Union{ConstantSource,VortexRing}}) = (4, 2)
_gpu_panel_tag(::Type{Union{ConstantSource,ConstantDoublet}}) = (5, 2)
_gpu_panel_tag(::Type) = nothing

# active filament regularization family as the functor's Int code
# (1=vatistas, 2=compact, 3=gaussian); enum order is Vatistas=0, Compact=1,
# Gaussian=2 (FLOWPanel_elements_fmm.jl FilamentRegularization)
_gpu_filament_reg() = @isdefined(FILAMENT_REGULARIZATION) ?
    Int(FILAMENT_REGULARIZATION[]) + 1 : 1

@inline function _gpu_pack_column!(A, q, tag, nv, verts, s1, s2, koff)
    @inbounds begin
        A[1, q] = tag
        A[2, q] = nv
        for iv in 1:4
            v = iv <= length(verts) ? verts[iv] : verts[end]
            A[3*iv, q] = v[1]
            A[3*iv + 1, q] = v[2]
            A[3*iv + 2, q] = v[3]
        end
        A[15, q] = s1
        A[16, q] = s2
        A[17, q] = koff
    end
    return nothing
end

# does this body carry a finite attached TE wake that direct!/_induced_wake
# would evaluate alongside its panels?
_gpu_has_te_wake(body::AbstractBody) = false
_gpu_has_te_wake(body::RigidWakeBody) =
    !body.semiinfinite_wake && !body.suppress_attached_wake[]

function _gpu_n_panel_columns(body::AbstractBody)
    n = body.ncells
    if _gpu_has_te_wake(body)
        for ic in 1:body.ncells
            body.shedding_full[1, ic] > 0 && (n += 2)
        end
    end
    return n
end

"""
    pack_panels!(srcmat, body; core_size=body.core_size)

Fill the 17-row `RectangularPanelInfluence` source matrix from an
`AbstractBody`'s panels (strengths from the current solution state,
`core_size` from the pass's active offset). RigidWakeBody finite attached
TE wakes are appended as extra ring/doublet tri columns — the same two
triangles `_induced_wake` (FLOWPanel_elements_fmm.jl:1205-1287) evaluates,
with the same `get_strength_doublet` + `wake_strength_shift` strength.
`srcmat` must have `>= 17` rows and `_gpu_n_panel_columns(body)` columns.
"""
function pack_panels!(srcmat::AbstractMatrix, body::AbstractBody{E};
        core_size::Real=body.core_size) where E
    tagnk = _gpu_panel_tag(E)
    tagnk === nothing && throw(ArgumentError(
        "unsupported element set $(E) for the rectangular panel functor"))
    tag, nk = tagnk
    SV = FastMultipole.StaticArrays.SVector{3,Float64}
    @inbounds for ic in 1:body.ncells
        i1, i2, i3 = body.cells[1, ic], body.cells[2, ic], body.cells[3, ic]
        v1 = SV(body.nodes[1, i1], body.nodes[2, i1], body.nodes[3, i1])
        v2 = SV(body.nodes[1, i2], body.nodes[2, i2], body.nodes[3, i2])
        v3 = SV(body.nodes[1, i3], body.nodes[2, i3], body.nodes[3, i3])
        s1 = body.strength[ic, 1]
        s2 = nk == 2 ? body.strength[ic, 2] : 0.0
        _gpu_pack_column!(srcmat, ic, tag, 3, (v1, v2, v3), s1, s2, core_size)
    end
    col = body.ncells
    if _gpu_has_te_wake(body)
        wtag = get_wake_kernel(body) === VortexRing ? 3 : 2   # tri ring or doublet
        @inbounds for ic in 1:body.ncells
            idx_1 = body.shedding_full[1, ic]
            idx_1 > 0 || continue
            nodes_idx = (body.cells[1, ic], body.cells[2, ic], body.cells[3, ic])
            TE1 = nodes_idx[idx_1]
            TE2 = nodes_idx[body.shedding_full[2, ic]]
            i_surf = body.shedding_full[3, ic]
            c1 = body.shedding_full[5, ic]
            c2 = body.shedding_full[6, ic]
            v1 = SV(body.nodes[1, TE1], body.nodes[2, TE1], body.nodes[3, TE1])
            v2 = SV(body.nodes[1, TE2], body.nodes[2, TE2], body.nodes[3, TE2])
            Da = SV(body.Das[i_surf][1, c1], body.Das[i_surf][2, c1], body.Das[i_surf][3, c1])
            Db = SV(body.Das[i_surf][1, c2], body.Das[i_surf][2, c2], body.Das[i_surf][3, c2])
            vw1 = v1 + Da
            vw2 = v2 + Db
            strength = get_strength_doublet(body, ic)
            if body.wake_correction_active[]
                strength += body.wake_strength_shift[ic]
            end
            # triangles (v1, v2, vw1) and (vw1, v2, vw2), matching _induced_wake
            _gpu_pack_column!(srcmat, col += 1, wtag, 3, (v1, v2, vw1),
                strength, 0.0, core_size)
            _gpu_pack_column!(srcmat, col += 1, wtag, 3, (vw1, v2, vw2),
                strength, 0.0, core_size)
        end
    end
    return srcmat
end

"""
    pack_particles!(srcmat, pfield)

Fill the 7-row `RectangularGaussianErfVortex` source matrix (position, Gamma,
sigma) from a `FLOWVPM.ParticleField` (raw particle rows X=1:3, Gamma=4:6,
sigma=7). `srcmat` must have `>= 7` rows and `pfield.np` columns.
"""
function pack_particles!(srcmat::AbstractMatrix, pfield::FLOWVPM.ParticleField)
    p = pfield.particles
    @inbounds for i in 1:pfield.np, r in 1:7
        srcmat[r, i] = p[r, i]
    end
    return srcmat
end

"""
    pack_filaments!(srcmat, filaments)

Fill the 17-row panel source matrix from a `FilamentWrapper{<:PanelWake}`'s
active final-row filaments as OPEN bound-vortex segments (tag 3, nv=2): the
same segment, strength (`_final_filament_strength`), and `core_size` its
`FastMultipole.direct!` (FLOWPanel_wake.jl:2869-2907) and
`source_system_to_buffer!` (:2790) use. `srcmat` must have `>= 17` rows and
`FastMultipole.get_n_bodies(filaments)` columns (nonzero only once
`wake.overflowed[]`).
"""
function pack_filaments!(srcmat::AbstractMatrix, filaments::FilamentWrapper{<:PanelWake})
    SV = FastMultipole.StaticArrays.SVector{3,Float64}
    wake = filaments.system
    n = FastMultipole.get_n_bodies(filaments)
    koff = wake.core_size
    i_row = wake.nwakes[]
    @inbounds for i in 1:n
        i_surf, j = fmm_to_filament_index(filaments, i)
        nod = wake.nodes[i_surf]
        v1 = SV(nod[1, i_row+1, j], nod[2, i_row+1, j], nod[3, i_row+1, j])
        v2 = SV(nod[1, i_row+1, j+1], nod[2, i_row+1, j+1], nod[3, i_row+1, j+1])
        s1 = _final_filament_strength(wake, i_surf, i_row, j)
        _gpu_pack_column!(srcmat, i, 3, 2, (v1, v2), s1, 0.0, koff)
    end
    return srcmat
end

"Fill the 17-row panel source matrix from a `PanelWake`'s quad vortex rings
(same enumeration as its `FastMultipole.source_system_to_buffer!`,
FLOWPanel_wake.jl:325-372; regularization radius = `core_size`)."
function pack_wake_panels!(srcmat::AbstractMatrix, wake::PanelWake)
    SV = FastMultipole.StaticArrays.SVector{3,Float64}
    n = FastMultipole.get_n_bodies(wake)
    koff = wake.core_size
    @inbounds for i in 1:n
        isurf, irow, icol = global_to_matrix_index(wake, i)
        nod = wake.nodes[isurf]
        v1 = SV(nod[1, irow, icol], nod[2, irow, icol], nod[3, irow, icol])
        v2 = SV(nod[1, irow+1, icol], nod[2, irow+1, icol], nod[3, irow+1, icol])
        v3 = SV(nod[1, irow+1, icol+1], nod[2, irow+1, icol+1], nod[3, irow+1, icol+1])
        v4 = SV(nod[1, irow, icol+1], nod[2, irow, icol+1], nod[3, irow, icol+1])
        s1 = wake.strength[isurf][1, irow, icol]
        _gpu_pack_column!(srcmat, i, 3, 4, (v1, v2, v3, v4), s1, 0.0, koff)
    end
    return srcmat
end

# number of packed source columns per system (0 = nothing to evaluate)
_gpu_source_columns(s::FLOWVPM.ParticleField) = s.np
_gpu_source_columns(s::PanelWake) = FastMultipole.get_n_bodies(s)
_gpu_source_columns(s::FilamentWrapper) = FastMultipole.get_n_bodies(s)
_gpu_source_columns(s::AbstractBody) = _gpu_n_panel_columns(s)

# pack (with buffer reuse) and return (srcmat, functor)
function _gpu_pack_source!(s::FLOWVPM.ParticleField, ::Nothing)
    if !(s.particles isa Array)
        # device-backed field: rows 1:7 are already the
        # RectangularGaussianErfVortex layout (X, Gamma, sigma); slice on
        # the device, no transfer (task 052)
        return s.particles[1:7, 1:s.np], FastMultipole.RectangularGaussianErfVortex()
    end
    m = _gpu_workmat!(s, :src, 7, s.np)
    pack_particles!(m, s)
    return m, FastMultipole.RectangularGaussianErfVortex()
end

function _gpu_pack_source!(s::PanelWake, ::Nothing)
    m = _gpu_workmat!(s, :src, 17, _gpu_source_columns(s))
    pack_wake_panels!(m, s)
    return m, FastMultipole.RectangularPanelInfluence(Int32(_gpu_filament_reg()))
end

function _gpu_pack_source!(s::FilamentWrapper{<:PanelWake}, ::Nothing)
    # key the work buffer on the wrapped wake: get_sources builds a FRESH
    # wrapper every call, so keying on the wrapper would defeat reuse and
    # grow the cache each step (role distinct from the wake's own :src)
    m = _gpu_workmat!(s.system, :src_filaments, 17, _gpu_source_columns(s))
    pack_filaments!(m, s)
    return m, FastMultipole.RectangularPanelInfluence(Int32(_gpu_filament_reg()))
end

function _gpu_pack_source!(s::AbstractBody, koff::Union{Nothing,Real})
    role = koff === nothing ? :src : :src_selfkoff
    m = _gpu_workmat!(s, role, 17, _gpu_source_columns(s))
    koff === nothing ? pack_panels!(m, s) : pack_panels!(m, s; core_size=koff)
    return m, FastMultipole.RectangularPanelInfluence(Int32(_gpu_filament_reg()))
end

################################################################################
# ELIGIBILITY
################################################################################

_gpu_all_zero(::Nothing) = true
_gpu_all_zero(x::Integer) = iszero(x)
_gpu_all_zero(x::Union{Tuple,AbstractVector}) = all(iszero, x)

_gpu_is_wake_source(s) =
    s isa FLOWVPM.ParticleField || s isa PanelWake || s isa FilamentWrapper

_gpu_source_supported(s) = false
_gpu_source_supported(s::FLOWVPM.ParticleField) =
    s.kernel === FLOWVPM.kernel_gaussianerf
_gpu_source_supported(s::PanelWake{TK}) where TK = TK === VortexRing
# non-PanelWake wrappers stay eligible only while empty; a PanelWake-wrapped
# wrapper packs its active final-row filaments (051 seam extension)
_gpu_source_supported(s::FilamentWrapper) = FastMultipole.get_n_bodies(s) == 0
_gpu_source_supported(s::FilamentWrapper{<:PanelWake}) = true
function _gpu_source_supported(s::AbstractBody{E}) where E
    _gpu_panel_tag(E) === nothing && return false
    if s isa RigidWakeBody && !s.suppress_attached_wake[]
        s.semiinfinite_wake && return false
        hasmethod(get_wake_kernel, Tuple{typeof(s)}) || return false
        get_wake_kernel(s) in (VortexRing, ConstantDoublet) || return false
    end
    return true
end

_gpu_target_supported(t) = false
_gpu_target_supported(t::AbstractBody) = true
_gpu_target_supported(t::ProbeWrapper{<:PanelWake}) = true
_gpu_target_supported(t::FLOWVPM.ParticleField) = true
_gpu_target_supported(t::FastMultipole.ProbeSystem) = true

################################################################################
# RESULT WRITE-BACK (same storage the fmm!/buffer path accumulates into)
################################################################################

# out rows: 1:3 velocity; 4:12 gradient with out[3+(j-1)*3+i] = du_i/dx_j

function _gpu_add_result!(body::AbstractBody, out, grad::Bool)
    n = FastMultipole.get_n_bodies(body)
    @inbounds for ic in 1:n
        body.velocity[1, ic] += out[1, ic]
        body.velocity[2, ic] += out[2, ic]
        body.velocity[3, ic] += out[3, ic]
        if grad
            for j in 1:3, i in 1:3
                body.velocity_gradient[i, j, ic] += out[3 + (j-1)*3 + i, ic]
            end
        end
    end
    return nothing
end

function _gpu_add_result!(pw::ProbeWrapper{<:PanelWake}, out, grad::Bool)
    wake = pw.system
    n = FastMultipole.get_n_bodies(pw)
    @inbounds for i in 1:n
        isurf, irow, icol = global_to_matrix_index(pw, i)
        wake.velocity[isurf][1, irow, icol] += out[1, i]
        wake.velocity[isurf][2, irow, icol] += out[2, i]
        wake.velocity[isurf][3, irow, icol] += out[3, i]
        # gradient not stored for PanelWake probes (matches
        # buffer_to_target_system!, FLOWPanel_wake.jl:436-459)
    end
    return nothing
end

function _gpu_add_result!(ps::FastMultipole.ProbeSystemStatic, out, grad::Bool)
    SV = FastMultipole.StaticArrays.SVector{3,Float64}
    SM = FastMultipole.StaticArrays.SMatrix{3,3,Float64,9}
    @inbounds for i in 1:FastMultipole.get_n_bodies(ps)
        ps.gradient[i] += SV(out[1, i], out[2, i], out[3, i])
        if grad
            # out[3+(j-1)*3+i0] = du_i0/dx_j == column-major SMatrix order
            ps.hessian[i] += SM(out[4, i], out[5, i], out[6, i],
                                out[7, i], out[8, i], out[9, i],
                                out[10, i], out[11, i], out[12, i])
        end
    end
    return nothing
end

function _gpu_add_result!(ps::FastMultipole.ProbeSystemArray, out, grad::Bool)
    @inbounds for i in 1:FastMultipole.get_n_bodies(ps)
        ps.gradient[1, i] += out[1, i]
        ps.gradient[2, i] += out[2, i]
        ps.gradient[3, i] += out[3, i]
        if grad
            for j in 1:3, i0 in 1:3
                ps.hessian[i0, j, i] += out[3 + (j-1)*3 + i0, i]
            end
        end
    end
    return nothing
end

function _gpu_add_result!(pf::FLOWVPM.ParticleField, out, grad::Bool)
    p = pf.particles
    @inbounds for i in 1:pf.np
        p[FLOWVPM.U_INDEX[1], i] += out[1, i]
        p[FLOWVPM.U_INDEX[2], i] += out[2, i]
        p[FLOWVPM.U_INDEX[3], i] += out[3, i]
        if grad
            for k in 1:9   # J[(j-1)*3+i] = du_i/dx_j, same order both sides
                p[FLOWVPM.J_INDEX[k], i] += out[3 + k, i]
            end
        end
    end
    return nothing
end

################################################################################
# THE SEAM
################################################################################

function _gpu_fill_targets!(tgt::Matrix{Float64}, sys)
    @inbounds for i in 1:size(tgt, 2)
        x = FastMultipole.get_position(sys, i)
        tgt[1, i] = x[1]
        tgt[2, i] = x[2]
        tgt[3, i] = x[3]
    end
    return tgt
end

# one accumulate call, host or device
function _gpu_direct_batch!(out::Matrix{Float64}, tgt::Matrix{Float64},
        packed::Vector{Any}, grad::Bool, device::Symbol)
    if device === :cuda
        CUDAmod = getglobal(FastMultipole, :CUDA)
        out_d = Base.invokelatest(CUDAmod.CuArray, out)
        tgt_d = Base.invokelatest(CUDAmod.CuArray, tgt)
        for (srcmat, kern) in packed
            # host-packed sources are uploaded; device slices (CuArray-backed
            # particle fields, task 052) are used as-is
            src_d = srcmat isa Matrix ?
                Base.invokelatest(CUDAmod.CuArray, srcmat) : srcmat
            Base.invokelatest(FastMultipole.direct_rectangular!, out_d, tgt_d,
                kern, src_d; gradient=grad)
        end
        copyto!(out, Base.invokelatest(Array, out_d))
    else
        for (srcmat, kern) in packed
            FastMultipole.direct_rectangular!(out, tgt, kern, srcmat; gradient=grad)
        end
    end
    return out
end

"""
    _gpu_rect_influence!(targets, sources; kwargs...)

The rectangular-direct seam body. Returns `true` iff the pass was recognized,
eligible, and fully evaluated (results accumulated into the same per-system
storage `fmm!`'s `buffer_to_target_system!` writes); `false` means the caller
must fall through to the existing `fmm!` path. See the file header for the
recognition and fallback rules.
"""
function _gpu_rect_influence!(targets::Tuple, sources::Tuple;
        scalar_potential=false, velocity=false, velocity_gradient=false,
        gradient=false, hessian=false,
        extra_outputs=nothing, direct_conditioning=nothing, kwargs...)
    gpu_influence_enabled() || return false
    scalar_potential === false || return false

    # probe pass (task 052): FastMultipole.ProbeSystem targets request the
    # potential gradient (= velocity) via `gradient`, and the velocity
    # gradient via `hessian` (BoundCirculationMonitor Kelvin slices,
    # KuttaJoukowskiForce, field probes).
    grad_all = gradient === true ||
        (gradient isa Tuple && length(gradient) > 0 && all(g -> g === true, gradient))
    probe_pass = velocity === false && grad_all && length(targets) > 0 &&
        all(t -> t isa FastMultipole.ProbeSystem, targets)
    if probe_pass
        velocity = true
        velocity_gradient = hessian
    end
    velocity === true || return false
    _gpu_all_zero(extra_outputs) || begin
        _gpu_info_once(:extra_outputs, "FLOWPanel gpu influence: pass requests " *
            "extra outputs (induced vorticity); falling back to fmm!")
        return false
    end

    # pass recognition
    pass1 = any(_gpu_is_wake_source, sources)
    pass3 = !pass1 && all(s -> s isa AbstractBody, sources) &&
        direct_conditioning !== nothing
    (pass1 || pass3 || probe_pass) || return false

    # eligibility
    for s in sources
        if !_gpu_source_supported(s)
            _gpu_info_once(:source, "FLOWPanel gpu influence: unsupported source " *
                "$(typeof(s)); falling back to fmm!")
            return false
        end
    end
    for t in targets
        if !_gpu_target_supported(t)
            _gpu_info_once(:target, "FLOWPanel gpu influence: unsupported target " *
                "$(typeof(t)); falling back to fmm!")
            return false
        end
    end
    # a self-paired SFS-enabled HOST particle field needs fmm! outputs for
    # Estr; a DEVICE (CuArray) field computes Estr through the radix
    # lifecycle below instead (task 052)
    for t in targets
        if t isa FLOWVPM.ParticleField && FLOWVPM.isSFSenabled(t.SFS) &&
                any(s -> s === t, sources) && (t.particles isa Array)
            _gpu_info_once(:sfs, "FLOWPanel gpu influence: SFS-enabled particle " *
                "self-influence needs fmm! Estr outputs; falling back to fmm! " *
                "(Estr stays on the FLOWVPM radix path)")
            return false
        end
    end

    device = _gpu_device()

    # a device-backed particle field anywhere in the pass REQUIRES the CUDA
    # seam: falling back to host fmm! (or the host rectangular path) would
    # scalar-index the CuArray. Fail loudly instead of silently degrading.
    if device !== :cuda
        for s in (targets..., sources...)
            if s isa FLOWVPM.ParticleField && !(s.particles isa Array)
                error("FLOWPanel gpu influence: CuArray-backed particle field " *
                    "requires FLOWPANEL_GPU_INFLUENCE=cuda with functional " *
                    "CUDA (mode=$(GPU_INFLUENCE[]), " *
                    "status: $(FastMultipole.cuda_radix_status()))")
            end
        end
    end

    t_start = time()

    for (it, tsys) in enumerate(targets)
        grad = velocity_gradient isa Tuple ? velocity_gradient[it] === true :
            velocity_gradient === true
        n_t = FastMultipole.get_n_bodies(tsys)
        n_t == 0 && continue

        if tsys isa FLOWVPM.ParticleField && !(tsys.particles isa Array)
            # device-resident particle target (task 052): radix self
            # influence + device rectangular for everything else
            _gpu_device_pfield_target!(tsys, sources, grad)
            continue
        end

        tgt = _gpu_workmat!(tsys, :targets, 3, n_t)
        _gpu_fill_targets!(tgt, tsys)
        out = _gpu_workmat!(tsys, grad ? :out12 : :out3, grad ? 12 : 3, n_t)
        fill!(out, 0.0)
        packed = Any[]
        for ssys in sources
            _gpu_source_columns(ssys) == 0 && continue
            # self-pair core_size conditioning (pass 3): evaluate a body's
            # own influence at core_size_panel, everything else at the
            # active (targets) offset — the rectangular analogue of
            # _self_panel_core_size_conditioning
            koff = (direct_conditioning !== nothing && ssys === tsys &&
                    ssys isa AbstractBody) ? ssys.core_size_panel : nothing
            push!(packed, _gpu_pack_source!(ssys, koff))
        end
        isempty(packed) && continue
        detail_label = if pass1
            tsys isa RigidWakeBody ? :wake_to_rotor_panels :
            tsys isa NonLiftingBody ? :wake_to_ground_panels : :wake_to_probes
        else
            tsys isa FLOWVPM.ParticleField && any(s -> s isa NonLiftingBody, sources) ?
                :ground_panels_to_particles :
            tsys isa FLOWVPM.ParticleField ? :rotor_panels_to_particles :
                :panel_cross_targets
        end
        _step_timer_measure(detail_label; nested=true) do
            _gpu_direct_batch!(out, tgt, packed, grad, device)
        end
        _gpu_add_result!(tsys, out, grad)
    end

    if gpu_timers_enabled()
        device === :cuda && _gpu_cuda_sync()
        label = probe_pass ? "influence_probe" :
                pass1 ? "influence_pass1" : "influence_pass3"
        _gpu_timer_log(label, time() - t_start)
    end

    GPU_INFLUENCE_HITS[] += 1
    return true
end

"Synchronize the CUDA device (timer accuracy); requires the lifecycle loaded."
function _gpu_cuda_sync()
    CUDAmod = getglobal(FastMultipole, :CUDA)
    Base.invokelatest(CUDAmod.synchronize)
    return nothing
end

"""
    _gpu_device_pfield_target!(pfield, sources, grad)

Task 052: accumulate a pass's influence into a CuArray-backed particle
field without leaving the device. The SELF pair (particles -> particles)
goes through the FLOWVPM radix GPU lifecycle — `UJ_fmm_gpu!` with
`reset=false` accumulates U and J, and with SFS enabled additionally
accumulates E_str into SFS_INDEX (matching the host `Estr_fmm!` +=
convention; the AfterUJ dynamic procedure runs in the post hook). All other
sources are evaluated by the device `direct_rectangular!` into a device
output buffer and accumulated into the U (and J when `grad`) rows.
"""
function _gpu_device_pfield_target!(pfield::FLOWVPM.ParticleField,
        sources::Tuple, grad::Bool)
    np = pfield.np
    np == 0 && return nothing

    # self influence via the device radix lifecycle
    if any(s -> s === pfield, sources)
        _step_timer_measure(:particle_self_uj; nested=true) do
            FLOWVPM.UJ_fmm_gpu!(pfield; reset=false, reset_sfs=false,
                sfs=FLOWVPM.isSFSenabled(pfield.SFS))
        end
    end

    packed = Any[]
    for ssys in sources
        ssys === pfield && continue
        _gpu_source_columns(ssys) == 0 && continue
        push!(packed, _gpu_pack_source!(ssys, nothing))
    end
    isempty(packed) && return nothing

    CUDAmod = getglobal(FastMultipole, :CUDA)
    P = pfield.particles
    tgt_d = P[1:3, 1:np]                       # device position slice
    nout = grad ? 12 : 3
    out_d = Base.invokelatest(CUDAmod.zeros, eltype(P), nout, np)
    for ((srcmat, kern), ssys) in zip(packed,
            Tuple(s for s in sources if s !== pfield && _gpu_source_columns(s) > 0))
        label = ssys isa NonLiftingBody ?
            :ground_panels_to_particles : :rotor_panels_to_particles
        _step_timer_measure(label; nested=true) do
            src_d = srcmat isa Matrix ?
                Base.invokelatest(CUDAmod.CuArray, srcmat) : srcmat
            Base.invokelatest(FastMultipole.direct_rectangular!, out_d, tgt_d,
                kern, src_d; gradient=grad)
        end
    end
    view(P, FLOWVPM.U_INDEX, 1:np) .+= view(out_d, 1:3, :)
    if grad
        # out rows 4:12 use out[3+(j-1)*3+i] = du_i/dx_j — the same order as
        # J_INDEX (see the host `_gpu_add_result!` above)
        view(P, FLOWVPM.J_INDEX, 1:np) .+= view(out_d, 4:12, :)
    end
    return nothing
end
