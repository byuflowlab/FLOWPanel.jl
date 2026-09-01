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

    Task 052d extension (panels -> particles host tree-FMM, default OFF):
      - when PANEL_INFLUENCE_FMM (or `set_panel_influence_fmm!(:on)`) is armed,
        AbstractBody sources of a ParticleField target — the
        `rotor_panels_to_particles` / `ground_panels_to_particles` legs — are
        evaluated by host `FastMultipole.fmm!` in ONE call (ProbeSystemArray
        target view; U always, J when the pass requests the gradient; never
        the scalar potential) instead of the dense rectangular kernels, timed
        under `:rotor_panels_to_particles_fmm`. Device-resident particle
        fields round-trip positions D2H and results H2D once per call. All
        other targets (bodies, probes) and all non-body sources are untouched.
        Tunables (read live): PANEL_INFLUENCE_FMM_P (default 8),
        PANEL_INFLUENCE_FMM_THETA (0.4), PANEL_INFLUENCE_FMM_LEAF (50).

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
const GPU_INFLUENCE_ROUTE_HITS = Dict{Symbol,Int}()
const GPU_INFLUENCE_ROUTE_FALLBACKS = Dict{Symbol,Int}()
const GPU_INFLUENCE_ROUTE_LOCK = ReentrantLock()

function reset_gpu_influence_routes!()
    lock(GPU_INFLUENCE_ROUTE_LOCK) do
        GPU_INFLUENCE_HITS[] = 0
        empty!(GPU_INFLUENCE_ROUTE_HITS)
        empty!(GPU_INFLUENCE_ROUTE_FALLBACKS)
    end
    return nothing
end

function gpu_influence_route_snapshot(; since=nothing)
    device = gpu_influence_enabled() ? _gpu_device() : :off
    total, cumulative_hits, cumulative_fallbacks = lock(GPU_INFLUENCE_ROUTE_LOCK) do
        (GPU_INFLUENCE_HITS[], copy(GPU_INFLUENCE_ROUTE_HITS),
            copy(GPU_INFLUENCE_ROUTE_FALLBACKS))
    end
    if since === nothing
        hits, fallbacks, total_hits = cumulative_hits, cumulative_fallbacks, total
    else
        hits = Dict(k => v - get(since.cumulative_hits, k, 0)
            for (k, v) in cumulative_hits if v != get(since.cumulative_hits, k, 0))
        fallbacks = Dict(k => v - get(since.cumulative_fallbacks, k, 0)
            for (k, v) in cumulative_fallbacks
            if v != get(since.cumulative_fallbacks, k, 0))
        total_hits = total - since.cumulative_total_hits
    end
    return (; requested_mode=GPU_INFLUENCE[], device, total_hits, hits, fallbacks,
        cumulative_total_hits=total, cumulative_hits, cumulative_fallbacks)
end

function _gpu_route_hit!(route::Symbol, device::Symbol)
    key = Symbol(device, :_, route)
    lock(GPU_INFLUENCE_ROUTE_LOCK) do
        GPU_INFLUENCE_HITS[] += 1
        GPU_INFLUENCE_ROUTE_HITS[key] =
            get(GPU_INFLUENCE_ROUTE_HITS, key, 0) + 1
    end
    return nothing
end

function _gpu_route_fallback!(route::Symbol, reason::AbstractString)
    lock(GPU_INFLUENCE_ROUTE_LOCK) do
        GPU_INFLUENCE_ROUTE_FALLBACKS[route] =
            get(GPU_INFLUENCE_ROUTE_FALLBACKS, route, 0) + 1
    end
    if GPU_INFLUENCE[] === :cuda &&
            !parse(Bool, get(ENV, "GPU_ALLOW_FALLBACK", "true"))
        error("FLOWPanel gpu influence: required route $route cannot use CUDA " *
            "($reason); GPU_ALLOW_FALLBACK=false")
    end
    return false
end

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

# out rows: 1:3 velocity; 4:12 gradient with out[3+(j-1)*3+i] = du_i/dx_j;
# optional scalar potential is row 4 without a gradient or row 13 with one.

function _gpu_add_result!(body::AbstractBody, out, grad::Bool,
        velocity::Bool=true, scalar_potential::Bool=false)
    n = FastMultipole.get_n_bodies(body)
    @inbounds for ic in 1:n
        if velocity
            body.velocity[1, ic] += out[1, ic]
            body.velocity[2, ic] += out[2, ic]
            body.velocity[3, ic] += out[3, ic]
        end
        if velocity && grad
            for j in 1:3, i in 1:3
                body.velocity_gradient[i, j, ic] += out[3 + (j-1)*3 + i, ic]
            end
        end
        if scalar_potential
            body.potential[ic] += out[grad ? 13 : 4, ic]
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
        packed::Vector{Any}, grad::Bool, scalar_potential::Bool,
        device::Symbol)
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
                kern, src_d; gradient=grad, scalar_potential)
        end
        copyto!(out, Base.invokelatest(Array, out_d))
    else
        for (srcmat, kern) in packed
            FastMultipole.direct_rectangular!(out, tgt, kern, srcmat;
                gradient=grad, scalar_potential)
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
        extra_outputs=nothing, direct_conditioning=nothing,
        production_route=nothing, kwargs...)
    gpu_influence_enabled() || return false

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
    want_potential = scalar_potential === true
    want_velocity = velocity === true
    (want_velocity || want_potential) || return false

    # pass recognition
    pass1 = any(_gpu_is_wake_source, sources)
    body_sources = !pass1 && all(s -> s isa AbstractBody, sources)
    pass3 = body_sources && direct_conditioning !== nothing
    block_cross = body_sources && !isempty(sources) &&
        all(t -> t isa AbstractBody &&
            all(s -> s !== t, sources), targets)
    panel_wake_potential = pass1 && want_potential &&
        all(s -> !(s isa FLOWVPM.ParticleField), sources)
    production_call = production_route !== nothing
    (pass1 || pass3 || probe_pass || block_cross || production_call) || return false
    route = production_call ? Symbol(production_route) :
        panel_wake_potential ? :wake_panel_potential :
        block_cross && want_potential ? :block_cross_potential :
        block_cross ? :block_cross_velocity :
        probe_pass ? :probe :
        pass1 ? :wake_velocity : :conditioned_body_velocity

    _gpu_all_zero(extra_outputs) || begin
        _gpu_info_once(:extra_outputs, "FLOWPanel gpu influence: pass requests " *
            "extra outputs (induced vorticity); falling back to fmm!")
        return _gpu_route_fallback!(route, "extra outputs requested")
    end
    want_potential && !all(s -> s isa AbstractBody || s isa PanelWake ||
        s isa FilamentWrapper{<:PanelWake}, sources) &&
        return _gpu_route_fallback!(route,
            "scalar potential requested from a vector-potential source")

    # eligibility
    for s in sources
        if !_gpu_source_supported(s)
            _gpu_info_once(:source, "FLOWPanel gpu influence: unsupported source " *
                "$(typeof(s)); falling back to fmm!")
            return _gpu_route_fallback!(route,
                "unsupported source $(typeof(s))")
        end
    end
    for t in targets
        if !_gpu_target_supported(t)
            _gpu_info_once(:target, "FLOWPanel gpu influence: unsupported target " *
                "$(typeof(t)); falling back to fmm!")
            return _gpu_route_fallback!(route,
                "unsupported target $(typeof(t))")
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
            return _gpu_route_fallback!(route,
                "host SFS particle self-influence requires FMM")
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
            want_potential && return _gpu_route_fallback!(route,
                "device particle targets do not store scalar potential")
            # device-resident particle target (task 052): radix self
            # influence + device rectangular for everything else
            _gpu_device_pfield_target!(tsys, sources, grad)
            continue
        end

        tgt = _gpu_workmat!(tsys, :targets, 3, n_t)
        _gpu_fill_targets!(tgt, tsys)
        nout = (grad ? 12 : 3) + want_potential
        outrole = want_potential ? (grad ? :out13 : :out4) :
            (grad ? :out12 : :out3)
        out = _gpu_workmat!(tsys, outrole, nout, n_t)
        fill!(out, 0.0)
        # 052d Phase 1 (host-target mirror of the device route below): body
        # sources of a ParticleField target go through the host tree-FMM in
        # one call; U-only or U+J, never the scalar potential. Small fixed
        # targets (bodies, probes) stay dense per the I3 split.
        fmm_bodies = (!pass1 && !want_potential &&
                tsys isa FLOWVPM.ParticleField &&
                panel_influence_fmm_enabled()) ?
            Tuple(s for s in sources if s isa AbstractBody &&
                _gpu_source_columns(s) > 0) : ()
        if !isempty(fmm_bodies)
            _step_timer_measure(:rotor_panels_to_particles_fmm; nested=true) do
                _panel_fmm_evaluate!(out, tgt, fmm_bodies, grad)
            end
            _panel_fmm_maybe_dump(out, tgt, fmm_bodies, n_t)
        end
        # 052h reverse leg: pass1 wake ParticleField sources of a BODY target
        # go through the device reverse cross pass (U-only, no potential);
        # on any fallback the sources stay in the dense pack below
        fmm_wakes = ()
        if pass1 && !want_potential && !grad && tsys isa AbstractBody &&
                panel_wake_fmm_enabled()
            cands = Tuple(s for s in sources
                if s isa FLOWVPM.ParticleField && s.np > 0)
            if !isempty(cands)
                wake_ok = false
                _step_timer_measure(:wake_to_panels_fmm; nested=true) do
                    wake_ok = _panel_wake_reverse_device!(out, tgt, tsys, cands)
                end
                wake_ok && (fmm_wakes = cands)
            end
        end
        packed = Any[]
        for ssys in sources
            _gpu_source_columns(ssys) == 0 && continue
            any(b -> b === ssys, fmm_bodies) && continue
            any(b -> b === ssys, fmm_wakes) && continue
            # self-pair core_size conditioning (pass 3): evaluate a body's
            # own influence at core_size_panel, everything else at the
            # active (targets) offset — the rectangular analogue of
            # _self_panel_core_size_conditioning
            koff = (direct_conditioning !== nothing && ssys === tsys &&
                    ssys isa AbstractBody) ? ssys.core_size_panel : nothing
            push!(packed, _gpu_pack_source!(ssys, koff))
        end
        isempty(packed) && isempty(fmm_bodies) && isempty(fmm_wakes) && continue
        detail_label = if pass1
            tsys isa RigidWakeBody ? :wake_to_rotor_panels :
            tsys isa NonLiftingBody ? :wake_to_ground_panels : :wake_to_probes
        else
            tsys isa FLOWVPM.ParticleField && any(s -> s isa NonLiftingBody, sources) ?
                :ground_panels_to_particles :
            tsys isa FLOWVPM.ParticleField ? :rotor_panels_to_particles :
                :panel_cross_targets
        end
        isempty(packed) || _step_timer_measure(detail_label; nested=true) do
            _gpu_direct_batch!(out, tgt, packed, grad, want_potential, device)
        end
        if tsys isa AbstractBody
            _gpu_add_result!(tsys, out, grad, want_velocity, want_potential)
        else
            want_potential && return _gpu_route_fallback!(route,
                "target $(typeof(tsys)) has no scalar-potential storage")
            _gpu_add_result!(tsys, out, grad)
        end
    end

    if gpu_timers_enabled()
        device === :cuda && _gpu_cuda_sync()
        label = probe_pass ? "influence_probe" :
                pass1 ? "influence_pass1" : "influence_pass3"
        _gpu_timer_log(label, time() - t_start)
    end

    _gpu_route_hit!(route, device)
    return true
end

"Synchronize the CUDA device (timer accuracy); requires the lifecycle loaded."
function _gpu_cuda_sync()
    CUDAmod = getglobal(FastMultipole, :CUDA)
    Base.invokelatest(CUDAmod.synchronize)
    return nothing
end

################################################################################
# 052d PHASE 1: HOST-FMM ROUTE FOR PANELS -> PARTICLES
################################################################################

# Default-off routing of the panels->particles legs (`rotor_panels_to_particles`
# / `ground_panels_to_particles`) through the host tree-FMM instead of the
# dense rectangular kernels. AbstractBody sources of a ParticleField target
# are gathered into ONE `FastMultipole.fmm!` call (one target tree over the
# particles, one combined multi-system source forest); every other
# source/target combination is untouched, so small fixed targets (probes,
# other bodies' control points) stay dense per the I3 split. Armed by the
# PANEL_INFLUENCE_FMM env var or `set_panel_influence_fmm!`.

"052d route state: :unset (read env on first query), :off, :on."
const PANEL_FMM = Ref{Symbol}(:unset)

"Whether the 052d panels->particles host-FMM route is armed (default off)."
function panel_influence_fmm_enabled()
    if PANEL_FMM[] === :unset
        val = lowercase(get(ENV, "PANEL_INFLUENCE_FMM", "0"))
        PANEL_FMM[] = val in ("1", "true", "on", "yes") ? :on : :off
    end
    return PANEL_FMM[] === :on
end

"Set the 052d panels->particles host-FMM route programmatically (:off, :on)."
function set_panel_influence_fmm!(mode::Symbol)
    mode in (:off, :on) || throw(ArgumentError(
        "panel influence fmm mode must be :off or :on (got $(repr(mode)))"))
    PANEL_FMM[] = mode
    return mode
end

# tunables, read live so probe jobs can sweep them between calls without a
# restart; defaults follow the 052d review response (p=8 starting point)
_panel_fmm_p() = parse(Int, get(ENV, "PANEL_INFLUENCE_FMM_P", "8"))
_panel_fmm_theta() = parse(Float64, get(ENV, "PANEL_INFLUENCE_FMM_THETA", "0.4"))
_panel_fmm_leaf() = parse(Int, get(ENV, "PANEL_INFLUENCE_FMM_LEAF", "50"))

################################################################################
# 052d STAGE F: DEVICE CROSS-PASS ROUTE (panels -> particles)
################################################################################
#
# When PANEL_INFLUENCE_FMM is armed AND the CUDA radix lifecycle is up, the
# panels->particles leg runs the fully device-resident cross pass
# (FastMultipole src/cross_stencil_cuda.jl: Stage A producers -> Stage B panel
# B2M/M2M with the TE-wake dipole arm -> Stage C M2L -> Stage D L2L/L2B ->
# Stage E block-sparse near field), accumulating velocity into U on device.
# The host-fmm! route (_panel_fmm_evaluate!) remains the fallback for grad
# requests (the cross pass is U-only by ruling), semi-infinite wakes, missing
# CUDA, or PANEL_INFLUENCE_FMM_DEVICE=0.
#
# Parity note (certified host-FMM semantics): NO wake_strength_shift on this
# leg — the host FMM's near field (_induced_wake buffer variant) and far
# field (RigidWakeBody B2M overload) both read the plain buffer dipole
# strength. The dense pack_panels! path's shifted wake columns are a
# pre-existing dense-vs-FMM asymmetry, not replicated here.
#
# Tunables (certified operating point p32g): PANEL_INFLUENCE_FMM_XQ=12,
# _XELL=5, _XP=6, _XRG=0.006. PANEL_INFLUENCE_FMM_XVERIFY=1 additionally runs
# the host-fmm! route each step and prints the device-vs-host relU (D3).

_panel_fmm_device() =
    lowercase(get(ENV, "PANEL_INFLUENCE_FMM_DEVICE", "1")) in
        ("1", "true", "on", "yes")
_panel_fmm_xverify() =
    lowercase(get(ENV, "PANEL_INFLUENCE_FMM_XVERIFY", "0")) in
        ("1", "true", "on", "yes")
_panel_fmm_xq() = parse(Int, get(ENV, "PANEL_INFLUENCE_FMM_XQ", "12"))
_panel_fmm_xell() = parse(Int, get(ENV, "PANEL_INFLUENCE_FMM_XELL", "5"))
_panel_fmm_xp() = parse(Int, get(ENV, "PANEL_INFLUENCE_FMM_XP", "6"))
_panel_fmm_xrg() = parse(Float64, get(ENV, "PANEL_INFLUENCE_FMM_XRG", "0.006"))

# ---- 052h reverse leg (particles -> panels) gates ----
#
# Two-tier like the forward pass: PANEL_WAKE_FMM (default off) arms the
# route; PANEL_WAKE_FMM_DEVICE (default on within) selects the device cross
# pass; PANEL_WAKE_FMM_XVERIFY additionally evaluates the dense reference and
# prints the relU (debug only). The reverse leg shares the forward tunables
# (_cross_config: XQ/XELL/XP/XRG) — one grid family, both directions.

"052h route state: :unset (read env on first query), :off, :on."
const PANEL_WAKE_FMM = Ref{Symbol}(:unset)

"Whether the 052h particles->panels reverse-FMM route is armed (default off)."
function panel_wake_fmm_enabled()
    if PANEL_WAKE_FMM[] === :unset
        val = lowercase(get(ENV, "PANEL_WAKE_FMM", "0"))
        PANEL_WAKE_FMM[] = val in ("1", "true", "on", "yes") ? :on : :off
    end
    return PANEL_WAKE_FMM[] === :on
end

"Set the 052h particles->panels reverse-FMM route programmatically."
function set_panel_wake_fmm!(mode::Symbol)
    mode in (:off, :on) || throw(ArgumentError(
        "panel wake fmm mode must be :off or :on (got $(repr(mode)))"))
    PANEL_WAKE_FMM[] = mode
    return mode
end

_panel_wake_fmm_device() =
    lowercase(get(ENV, "PANEL_WAKE_FMM_DEVICE", "1")) in
        ("1", "true", "on", "yes")
_panel_wake_fmm_xverify() =
    lowercase(get(ENV, "PANEL_WAKE_FMM_XVERIFY", "0")) in
        ("1", "true", "on", "yes")

# All cross-pass tunables parsed at once into an immutable, validated config.
# The cached-entry key includes it (_cross_entry!), so changing any tunable
# mid-run rebuilds the whole entry instead of silently mixing old tables with
# new parameters — and the multipole/local expansion orders can never skew.
const _CROSS_CFG_LOGGED = Ref(false)
function _cross_config()
    q = _panel_fmm_xq()
    ell_x = _panel_fmm_xell()
    P = _panel_fmm_xp()
    R_guard = _panel_fmm_xrg()
    q >= 0 || throw(ArgumentError(
        "PANEL_INFLUENCE_FMM_XQ must be >= 0 (got $q)"))
    2 <= ell_x <= 8 || throw(ArgumentError(
        "PANEL_INFLUENCE_FMM_XELL must be in 2:8 (got $ell_x)"))
    # upper bound = FastMultipole._CROSS_B2M_MAX_P (device B2M scratch size)
    1 <= P <= 8 || throw(ArgumentError(
        "PANEL_INFLUENCE_FMM_XP must be in 1:8 (got $P)"))
    isfinite(R_guard) && R_guard >= 0 || throw(ArgumentError(
        "PANEL_INFLUENCE_FMM_XRG must be finite and >= 0 (got $R_guard)"))
    return (; q, ell_x, P, R_guard)
end

# per-body persistent device state: objectid(body) => mutable holder.
# `x_min`/`h0` are FROZEN at first build (padded union box) so the cached M2L
# probe tables AND the h0-keyed stencil demotion masks stay valid across
# steps; if bodies escape the box, refresh_cross_producers! signals
# `needs_rebuild` and _cross_run_body! reconstructs the WHOLE entry.
mutable struct _CrossPassEntry
    ns::Int
    np::Int
    x_min::Any
    h0::Float64
    ct::Any            # CrossStencilTables
    ctx::Any           # DeviceCrossProducerContext
    xs::Any            # DeviceCrossExpansionState
    ls::Any            # DeviceCrossLocalState
    m2l_tables::Any    # host (ops, class_slot) for the CURRENT h0
    l2l_table::Any     # host L2L array for the CURRENT h0
    srcmat::Matrix{Float64}
    wakemat::Union{Nothing,Matrix{Float64}}
    cent::Matrix{Float64}
    cfg::Any           # validated tunable NamedTuple this entry was built with
    # persistent device mirrors of the per-step host inputs (uploaded with
    # copyto! each step; reallocated only on a full entry rebuild)
    d_srcmat::Any
    d_cent::Any
    d_wakemat::Any
end

const _CROSS_STATE = Dict{UInt,_CrossPassEntry}()

# 052 leak fix: persistent padded device particle-position mirror (3 x np_cap
# Float64) plus its grow-only capacity. See _gpu_panels_to_particles_cross for
# why the pass runs on a padded capacity rather than the live np.
const _CROSS_POS = Ref{Any}(nothing)
const _CROSS_NP_CAP = Ref(0)

"Drop all cached cross-pass device states (e.g. after geometry rebuilds)."
clear_cross_pass_state!() = (empty!(_CROSS_STATE); _CROSS_POS[] = nothing;
    _CROSS_NP_CAP[] = 0; clear_reverse_wake_state!(); nothing)

# Refresh (and grow, ~6% geometric headroom) the padded position mirror: the
# live prefix converts the particle field's storage eltype to Float64 on
# device; the tail wrap-fills with copies of REAL particle positions so the
# padded set occupies exactly the same cells as the live set (identical
# node/route/block lists, no artificial hot node).
function _cross_padded_positions!(Pd, np::Int, CUDAmod)
    if np > _CROSS_NP_CAP[] || _CROSS_POS[] === nothing
        _CROSS_NP_CAP[] = np + max(np >> 4, 1024)
        _CROSS_POS[] = Base.invokelatest(CUDAmod.CuArray{Float64}, undef, 3,
            _CROSS_NP_CAP[])
    end
    np_cap = _CROSS_NP_CAP[]
    pos = _CROSS_POS[]
    view(pos, :, 1:np) .= view(Pd, 1:3, 1:np)
    i = np
    while i < np_cap
        k = min(np, np_cap - i)
        view(pos, :, (i + 1):(i + k)) .= view(pos, :, 1:k)
        i += k
    end
    return pos, np_cap
end

# fill the per-step host-side inputs for one body: 17-row body-panel columns
# (NO appended wake columns — the cross pass carries the wake through the
# 8-row wake matrix instead), panel centroids for keying, and the wake matrix
# (idx1, Da, idx2, Db per panel; -1 idx when the panel does not shed)
function _cross_fill_inputs!(e::_CrossPassEntry, body::AbstractBody{E}) where E
    tag, nk = _gpu_panel_tag(E)
    SV = FastMultipole.StaticArrays.SVector{3,Float64}
    ns = body.ncells
    @inbounds for ic in 1:ns
        i1, i2, i3 = body.cells[1, ic], body.cells[2, ic], body.cells[3, ic]
        v1 = SV(body.nodes[1, i1], body.nodes[2, i1], body.nodes[3, i1])
        v2 = SV(body.nodes[1, i2], body.nodes[2, i2], body.nodes[3, i2])
        v3 = SV(body.nodes[1, i3], body.nodes[2, i3], body.nodes[3, i3])
        s1 = body.strength[ic, 1]
        s2 = nk == 2 ? body.strength[ic, 2] : 0.0
        _gpu_pack_column!(e.srcmat, ic, tag, 3, (v1, v2, v3), s1, s2,
            body.core_size)
        c = (v1 + v2 + v3) / 3
        e.cent[1, ic] = c[1]; e.cent[2, ic] = c[2]; e.cent[3, ic] = c[3]
    end
    if e.wakemat !== nothing
        w = e.wakemat
        fill!(w, -1.0)
        @inbounds for ic in 1:ns
            idx1 = body.shedding_full[1, ic]
            w[1, ic] = idx1
            idx1 > 0 || continue
            i_surf = body.shedding_full[3, ic]
            c1 = body.shedding_full[5, ic]
            c2 = body.shedding_full[6, ic]
            w[2, ic] = body.Das[i_surf][1, c1]
            w[3, ic] = body.Das[i_surf][2, c1]
            w[4, ic] = body.Das[i_surf][3, c1]
            w[5, ic] = body.shedding_full[2, ic]
            w[6, ic] = body.Das[i_surf][1, c2]
            w[7, ic] = body.Das[i_surf][2, c2]
            w[8, ic] = body.Das[i_surf][3, c2]
        end
    end
    return e
end

# 052d relU attribution (diagnostic only, 2026-08-29): with PANEL_FMM_DUMP_DIR
# set and PANEL_FMM_DUMP_NP=np1,np2,... a comma list of particle-count
# thresholds, the host-fmm! leg dumps — at the first step whose np reaches each
# threshold — the cross-pass inputs (packed 17-row srcmat, 8-row wake matrix,
# centroids), the particle positions, and the host-fmm! velocity result, as
# flat column-major Float64 bins (snapshot472 conventions). Consumed by
# MATRIX_OPERATOR_REFACTOR/prototypes/052d_cross_stencil/p39_relU_attribution.jl.
# The device xverify path passes `deviceU` (3 x np device cross-pass result),
# dumped alongside so the attribution needs no prototype re-run: the dump then
# holds the exact per-step state plus BOTH route results from the SAME step of
# the SAME trajectory (no replay divergence).
const _PANEL_FMM_DUMPED = Set{Int}()

# 052d step-5e (diagnostic, 2026-08-29): eval-time state capture inside
# _cross_run_body!, same env gating as _panel_fmm_maybe_dump but a separate
# claimed set. Dumps the entry's OWN packed inputs (what the device actually
# uploads), the frozen box, list counters, and the far-only d_out so the
# production-vs-dump input question is settled bitwise.
const _PANEL_CROSS_EVAL_DUMPED = Set{Int}()
function _panel_cross_dump_hit(np::Int)
    dir = get(ENV, "PANEL_FMM_DUMP_DIR", "")
    isempty(dir) && return 0, dir
    thresholds = [parse(Int, s) for s in
        split(get(ENV, "PANEL_FMM_DUMP_NP", ""), ","; keepempty=false)]
    for t in thresholds
        (np >= t && !(t in _PANEL_CROSS_EVAL_DUMPED)) && return t, dir
    end
    return 0, dir
end
function _panel_fmm_maybe_dump(out, tgt, fmm_bodies, np::Int; deviceU=nothing)
    dir = get(ENV, "PANEL_FMM_DUMP_DIR", "")
    isempty(dir) && return nothing
    thresholds = [parse(Int, s) for s in
        split(get(ENV, "PANEL_FMM_DUMP_NP", ""), ","; keepempty=false)]
    hit = 0
    for t in thresholds
        if np >= t && !(t in _PANEL_FMM_DUMPED)
            hit = t
            break
        end
    end
    hit == 0 && return nothing
    push!(_PANEL_FMM_DUMPED, hit)
    mkpath(dir)
    pre = joinpath(dir, "dump_np$(np)")
    write(pre * "_positions_3xN_f64.bin", Matrix{Float64}(tgt[1:3, 1:np]))
    write(pre * "_hostU_3xN_f64.bin", Matrix{Float64}(out[1:3, 1:np]))
    deviceU === nothing ||
        write(pre * "_deviceU_3xN_f64.bin", Matrix{Float64}(deviceU[1:3, 1:np]))
    open(pre * "_meta.txt", "w") do io
        println(io, "np=", np)
        println(io, "nbodies=", length(fmm_bodies))
        println(io, "cfg=", _cross_config())
    end
    for (i, body) in enumerate(fmm_bodies)
        ns = body.ncells
        has_wake = body isa RigidWakeBody && _gpu_has_te_wake(body)
        e = _CrossPassEntry(ns, np, nothing, 0.0, nothing, nothing, nothing,
            nothing, nothing, nothing, zeros(17, ns),
            has_wake ? zeros(8, ns) : nothing, zeros(3, ns), nothing,
            nothing, nothing, nothing)
        _cross_fill_inputs!(e, body)
        write(pre * "_body$(i)_srcmat_17xS_f64.bin", e.srcmat)
        e.wakemat === nothing ||
            write(pre * "_body$(i)_wakemat_8xS_f64.bin", e.wakemat)
        write(pre * "_body$(i)_cent_3xS_f64.bin", e.cent)
        open(pre * "_meta.txt", "a") do io
            println(io, "body$(i)=", typeof(body), " ns=", ns,
                " has_wake=", has_wake, " core_size=", body.core_size)
        end
    end
    println("panel_fmm_dump np=$np -> $(pre)_*")
    flush(stdout)
    return nothing
end

# padded cubic union root box over panel nodes + current device particle slice
function _cross_root_box(body::AbstractBody, pos_d, CUDAmod)
    lo = zeros(3); hi = zeros(3)
    for a in 1:3
        lo[a] = min(minimum(view(body.nodes, a, :)),
            Base.invokelatest(minimum, view(pos_d, a, :)))
        hi[a] = max(maximum(view(body.nodes, a, :)),
            Base.invokelatest(maximum, view(pos_d, a, :)))
    end
    ctr = (lo .+ hi) ./ 2
    h0 = 0.75 * maximum(hi .- lo) + 1e-9   # 1.5x padding on the half-width
    x_min = FastMultipole.StaticArrays.SVector{3,Float64}(ctr .- h0)
    return x_min, h0
end

function _cross_entry!(body::AbstractBody, np::Int, pos_d, CUDAmod)
    key = objectid(body)
    ns = body.ncells
    cfg = _cross_config()
    e = get(_CROSS_STATE, key, nothing)
    has_wake = body isa RigidWakeBody && _gpu_has_te_wake(body)
    if e === nothing || e.ns != ns || (e.wakemat !== nothing) != has_wake ||
            e.cfg != cfg
        x_min, h0 = _cross_root_box(body, pos_d, CUDAmod)
        ct = FastMultipole.CrossStencilTables(cfg.q, cfg.ell_x, h0, cfg.R_guard)
        e = _CrossPassEntry(ns, 0, x_min, h0, ct, nothing, nothing, nothing,
            nothing, nothing, zeros(17, ns),
            has_wake ? zeros(8, ns) : nothing, zeros(3, ns), cfg,
            nothing, nothing, nothing)
        e.ctx = FastMultipole.device_cross_producer_context(ct, x_min, h0, ns, np)
        e.np = np
        e.xs = FastMultipole.device_cross_expansion_state(e.ctx, cfg.P)
        e.m2l_tables = FastMultipole.cross_m2l_operators(cfg.P, h0, ct)
        e.l2l_table = FastMultipole.cross_l2l_operators(cfg.P, h0, ct.ell_x)
        e.ls = FastMultipole.device_cross_local_state(e.ctx, cfg.P;
            m2l_tables=e.m2l_tables, l2l_table=e.l2l_table)
        # persistent device input mirrors, refreshed with copyto! each step
        e.d_srcmat = Base.invokelatest(CUDAmod.CuArray, e.srcmat)
        e.d_cent = Base.invokelatest(CUDAmod.CuArray, e.cent)
        e.d_wakemat = has_wake ?
            Base.invokelatest(CUDAmod.CuArray, e.wakemat) : nothing
        _CROSS_STATE[key] = e
        if !_CROSS_CFG_LOGGED[]
            _CROSS_CFG_LOGGED[] = true
            println("panel_cross config: q=$(cfg.q) ell_x=$(cfg.ell_x) " *
                "P=$(cfg.P) R_guard=$(cfg.R_guard)")
            flush(stdout)
        end
    elseif e.np != np
        # particle CAPACITY changed (np here is the padded np_cap, which grows
        # only when the live count outgrows it — every ~6% of growth, not every
        # shedding step): rebuild the two-occupancy context and the
        # particle-sized local state; panel-side expansion state, the frozen
        # root box (=> cached operator tables), and the config carry over —
        # cfg equality above guarantees e.xs was built with e.cfg.P, so the
        # multipole/local expansion orders stay consistent
        e.ctx = FastMultipole.device_cross_producer_context(e.ct, e.x_min, e.h0,
            ns, np)
        e.np = np
        e.ls = FastMultipole.device_cross_local_state(e.ctx, e.cfg.P;
            m2l_tables=e.m2l_tables, l2l_table=e.l2l_table)
        # 052 leak fix: the replaced ctx+ls (~10^2 MB of device arrays) lived
        # many steps, so they are old-generation garbage the host GC will
        # never collect on its own (device bytes don't pressure its
        # heuristics). Rebuilds are rare now — collect eagerly.
        GC.gc(true)
    end
    return e
end

# One body's Stages A–E. ALL fallible work (uploads, producer refresh, box
# rebuild, expansion passes, near field) runs BEFORE the single `.+=` into the
# particle field, so a throw anywhere earlier leaves U untouched for this body.
function _cross_run_body!(e::_CrossPassEntry, body::AbstractBody, np::Int,
        pos_d, Pd, xver, CUDAmod)
    _cross_fill_inputs!(e, body)
    copyto!(e.d_srcmat, e.srcmat)
    copyto!(e.d_cent, e.cent)
    e.wakemat === nothing || copyto!(e.d_wakemat, e.wakemat)
    FastMultipole.refresh_cross_producers!(e.ctx, e.d_cent, pos_d)
    if e.ctx.needs_rebuild
        # box escape: the stencil demotion masks and cached operator tables
        # are keyed on h0, so patching x_min/h0 into live state would misroute
        # far/near work. Rebuild the ENTIRE entry around a fresh padded union
        # box (new tables, masks, contexts, operators) and rerun the producers.
        # (rebuild at the PADDED count — pos_d carries np_cap columns)
        delete!(_CROSS_STATE, objectid(body))
        e = _cross_entry!(body, size(pos_d, 2), pos_d, CUDAmod)
        _cross_fill_inputs!(e, body)
        copyto!(e.d_srcmat, e.srcmat)
        copyto!(e.d_cent, e.cent)
        e.wakemat === nothing || copyto!(e.d_wakemat, e.wakemat)
        FastMultipole.refresh_cross_producers!(e.ctx, e.d_cent, pos_d)
        e.ctx.needs_rebuild && error("cross-pass union-box rebuild failed to " *
            "contain all bodies (freshly padded box escaped within one step)")
    end
    hit, dumpdir = _panel_cross_dump_hit(np)
    if hit != 0
        push!(_PANEL_CROSS_EVAL_DUMPED, hit)
        mkpath(dumpdir)
        pre = joinpath(dumpdir, "ateval_np$(np)")
        write(pre * "_srcmat.bin", e.srcmat)
        e.wakemat === nothing || write(pre * "_wakemat.bin", e.wakemat)
        write(pre * "_cent.bin", e.cent)
        write(pre * "_positions.bin", Matrix{Float64}(Array(pos_d))[:, 1:np])
        open(pre * "_meta.txt", "w") do io
            println(io, "np=", np, " ns=", e.ns)
            println(io, "x_min=", e.x_min)
            println(io, "h0=", e.h0)
            println(io, "n_routes=", e.ctx.n_routes, " n_demoted=",
                e.ctx.n_demoted, " n_blocks=", e.ctx.n_blocks)
        end
        # 052d step-5f: dump the producer route/block lists + occupancy so the
        # production near pair set can be bit-compared against the proto shell
        lists = FastMultipole.download_cross_lists(e.ctx)
        open(pre * "_lists_meta.txt", "w") do io
            wdump = (name, arr) -> begin
                A = collect(arr)
                write(pre * "_lists_" * name * ".bin", A)
                println(io, name, " ", eltype(A), " ", join(size(A), "x"))
            end
            for grp in (:routes, :blocks, :panels, :particles)
                sub = getproperty(lists, grp)
                for f in propertynames(sub)
                    v = getproperty(sub, f)
                    v isa AbstractArray ? wdump(string(grp, "_", f), v) :
                        println(io, grp, "_", f, " = ", v)
                end
            end
            println(io, "x_min = ", lists.x_min)
            println(io, "h0 = ", lists.h0)
            println(io, "needs_rebuild = ", lists.needs_rebuild)
        end
        println("panel_cross_eval_dump np=$np -> $(pre)_*")
        flush(stdout)
    end
    FastMultipole.refresh_cross_multipoles!(e.xs, e.ctx, e.d_srcmat, e.d_wakemat)
    # Stage-B skip counter is a hard error: a skipped panel (unhandled tag or
    # nv < 3) would silently lose its far field while Stage E still evaluates
    # its near field — a plausible-but-wrong velocity
    e.xs.n_skipped == 0 || error("cross-pass Stage B skipped " *
        "$(e.xs.n_skipped) panel(s) (unhandled tag or nv < 3); refusing to " *
        "return a partial far field")
    FastMultipole.refresh_cross_locals!(e.ls, e.ctx, e.xs)
    FastMultipole.finish_cross_locals!(e.ls, e.ctx, pos_d)
    if hit != 0
        write(joinpath(dumpdir, "ateval_np$(np)_farU.bin"),
            Matrix{Float64}(Array(e.ls.d_out))[2:4, 1:np])
    end
    wake_tag = (body isa RigidWakeBody &&
        get_wake_kernel(body) === ConstantDoublet) ? 2 : 3
    FastMultipole.apply_cross_near!(e.ls, e.ctx, e.d_srcmat, pos_d;
        d_wake_buffer=e.d_wakemat, wake_tag=wake_tag,
        reg=_gpu_filament_reg(), potential=false)
    # all fallible work done — accumulate into the particle field LAST
    # (live prefix only: columns np+1:np_cap are padding targets)
    view(Pd, FLOWVPM.U_INDEX, 1:np) .+= view(e.ls.d_out, 2:4, 1:np)
    if xver !== nothing
        xver .+= Matrix{Float64}(Array(e.ls.d_out))[2:4, 1:np]
    end
    return nothing
end

"""
    _panel_cross_device!(pfield, bodies, grad) -> Bool

052d Stage F: run the device cross pass for every body in `bodies`,
accumulating velocity into the particle field's U rows on device. Returns
false (host-fmm! fallback) without touching anything when the leg cannot run
on device. With PANEL_INFLUENCE_FMM_XVERIFY=1, additionally evaluates the
host-fmm! reference and prints the device-vs-host relU each call.
"""
function _panel_cross_device!(pfield::FLOWVPM.ParticleField, bodies::Tuple,
        grad::Bool)
    _panel_fmm_device() ||
        return _gpu_route_fallback!(:panels_to_particles_cross, "device route disarmed")
    grad &&
        return _gpu_route_fallback!(:panels_to_particles_cross,
            "cross pass is U-only (grad requested)")
    FastMultipole.load_cuda_radix_lifecycle!() ||
        return _gpu_route_fallback!(:panels_to_particles_cross,
            "CUDA radix lifecycle unavailable")
    for body in bodies
        body isa AbstractBody ||
            return _gpu_route_fallback!(:panels_to_particles_cross,
                "unsupported source $(typeof(body))")
        _gpu_panel_tag(_gpu_element_type(body)) === nothing &&
            return _gpu_route_fallback!(:panels_to_particles_cross,
                "unsupported element set $(typeof(body))")
        body isa RigidWakeBody && body.semiinfinite_wake &&
            return _gpu_route_fallback!(:panels_to_particles_cross,
                "semi-infinite wake not covered by the cross pass")
    end
    # validate tunables up front so a misconfigured env fails loudly here
    # rather than surfacing as a silent per-step host fallback
    _cross_config()
    np = pfield.np
    CUDAmod = getglobal(FastMultipole, :CUDA)
    Pd = pfield.particles
    # cross-pass state is Float64; production particle fields may be Float32.
    # 052 leak fix: np grows every shedding step, and rebuilding the
    # particle-sized device state (producer context + local state, ~10^2 MB)
    # each step leaks it — the replaced buffers survive a full step, get
    # promoted, and no major GC ever runs (device bytes are invisible to the
    # host GC heuristics; job 13508681 died at step 819 this way). Pad the
    # particle count to a grow-only capacity instead: the tail wrap-fills with
    # REAL particle positions, so occupancy (and thus the route/block lists)
    # is identical and the only cost is ~6% duplicate target work; entries
    # rebuild only when np outgrows the capacity.
    pos_d, np_cap = _cross_padded_positions!(Pd, np, CUDAmod)
    xver = _panel_fmm_xverify() ? zeros(3, np) : nothing
    n_accumulated = 0
    pass_error = nothing
    _step_timer_measure(:rotor_panels_to_particles_cross; nested=true) do
        try
            for body in bodies
                e = _cross_entry!(body, np_cap, pos_d, CUDAmod)
                _cross_run_body!(e, body, np, pos_d, Pd, xver, CUDAmod)
                n_accumulated += 1
            end
        catch err
            pass_error = err
        end
    end
    if pass_error !== nothing
        if n_accumulated == 0
            # U untouched: fall back to the host-fmm! route cleanly
            return _gpu_route_fallback!(:panels_to_particles_cross,
                "device cross pass threw before any accumulation: " *
                sprint(showerror, pass_error))
        end
        # some bodies already accumulated: the field is corrupt — refuse to
        # continue or to double-count via a fallback
        error("device cross pass failed after accumulating $(n_accumulated) " *
            "of $(length(bodies)) bodies into U (partial write; the particle " *
            "field must not be trusted): " * sprint(showerror, pass_error))
    end
    if xver !== nothing
        pos_h = Matrix{Float64}(Array(pos_d))[:, 1:np]
        ref = zeros(3, np)
        _panel_fmm_evaluate!(ref, pos_h, bodies, false)
        relU = sqrt(sum(abs2, xver .- ref)) / max(sqrt(sum(abs2, ref)), eps())
        println("panel_cross_xverify np=$np relU=$relU")
        flush(stdout)
        _panel_fmm_maybe_dump(ref, pos_h, bodies, np; deviceU=xver)
    end
    _gpu_route_hit!(:panels_to_particles_cross, :cuda)
    return true
end

################################################################################
# 052h REVERSE LEG: DEVICE CROSS-PASS ROUTE (particles -> panels)
################################################################################
#
# When PANEL_WAKE_FMM is armed AND the CUDA radix lifecycle is up, the pass1
# wake ParticleField -> body-panel legs run the device reverse cross pass
# (FastMultipole cross_stencil_cuda.jl 052h: reverse producers -> vortex B2M +
# LH M2M over particle nodes -> LH M2L -> LH L2L + dual-channel L2B at the
# panel control points -> vortex near blocks), replacing the dense
# :wake_to_rotor_panels / :wake_to_ground_panels legs. U-only (grad requests
# keep those sources dense); probe targets stay dense in this first landing
# (probe-system identity across steps is unresolved — a per-step entry
# rebuild would swamp the win). Everything else about the entry pattern
# mirrors the forward pass: frozen padded union root box, grow-only padded
# particle buffers (padding particles carry wrapped REAL positions but ZERO
# strength, so occupancy stays step-stable while padded sources contribute
# nothing), rebuild-on-escape, fallback before any accumulation.

mutable struct _ReverseWakeEntry
    nt::Int
    np::Int            # padded capacity the ctx was built with
    x_min::Any
    h0::Float64
    ct::Any            # CrossStencilTables
    ctx::Any           # DeviceCrossProducerContext (build_reverse = true)
    rs::Any            # DeviceCrossReverseState
    cfg::Any
    d_tgt::Any         # 3 x nt device mirror of the target positions
end

const _REVERSE_STATE = Dict{UInt,_ReverseWakeEntry}()

# grow-only padded 8-row particle SOURCE buffer (FMM layout: 1:3 position,
# 5:7 strength, 8 sigma); shares the position padding discipline of
# _CROSS_POS but zeroes the padded strengths
const _REV_BUF = Ref{Any}(nothing)
const _REV_BUF_CAP = Ref(0)

"Drop all cached reverse-leg device states (e.g. after geometry rebuilds)."
clear_reverse_wake_state!() = (empty!(_REVERSE_STATE); _REV_BUF[] = nothing;
    _REV_BUF_CAP[] = 0; nothing)

function _reverse_padded_buffer!(Pd, np::Int, np_cap::Int, pos_d, CUDAmod)
    if np_cap > _REV_BUF_CAP[] || _REV_BUF[] === nothing
        _REV_BUF_CAP[] = np_cap
        _REV_BUF[] = Base.invokelatest(CUDAmod.zeros, Float64, 8, np_cap)
    end
    buf = _REV_BUF[]
    # positions: reuse the padded position mirror wholesale (already wrapped)
    view(buf, 1:3, 1:np_cap) .= view(pos_d, :, 1:np_cap)
    # live strengths + sigma; padding stays zero-strength, sigma 1
    view(buf, 5:7, 1:np) .= view(Pd, 4:6, 1:np)
    view(buf, 8:8, 1:np) .= view(Pd, 7:7, 1:np)
    if np < np_cap
        view(buf, 5:7, (np + 1):np_cap) .= 0.0
        view(buf, 8:8, (np + 1):np_cap) .= 1.0
    end
    return buf
end

function _reverse_entry!(body::AbstractBody, tgt::AbstractMatrix{Float64},
        np_cap::Int, pos_d, CUDAmod)
    key = objectid(body)
    nt = size(tgt, 2)
    cfg = _cross_config()
    e = get(_REVERSE_STATE, key, nothing)
    if e === nothing || e.nt != nt || e.cfg != cfg
        x_min, h0 = _cross_root_box(body, pos_d, CUDAmod)
        ct = FastMultipole.CrossStencilTables(cfg.q, cfg.ell_x, h0, cfg.R_guard)
        e = _ReverseWakeEntry(nt, 0, x_min, h0, ct, nothing, nothing, cfg,
            nothing)
        e.ctx = FastMultipole.device_cross_producer_context(ct, x_min, h0,
            nt, np_cap; build_reverse=true)
        e.np = np_cap
        e.rs = FastMultipole.device_cross_reverse_state(e.ctx, cfg.P)
        e.d_tgt = Base.invokelatest(CUDAmod.CuArray, tgt)
        _REVERSE_STATE[key] = e
    elseif e.np != np_cap
        # padded particle capacity outgrown: rebuild the two-occupancy context
        # + reverse state around the SAME frozen box/config (mirrors
        # _cross_entry!; the LH operator tables depend only on h0/P and are
        # rebuilt inside device_cross_reverse_state)
        e.ctx = FastMultipole.device_cross_producer_context(e.ct, e.x_min,
            e.h0, nt, np_cap; build_reverse=true)
        e.np = np_cap
        e.rs = FastMultipole.device_cross_reverse_state(e.ctx, e.cfg.P)
        GC.gc(true)   # 052 leak fix: eagerly collect the replaced device state
    end
    return e
end

# One target body's reverse Stages A-E. ALL fallible work runs before the
# single accumulation into the host scratch.
function _reverse_run!(e::_ReverseWakeEntry, body::AbstractBody,
        tgt::AbstractMatrix{Float64}, np::Int, pos_d, d_pbuf, acc, CUDAmod)
    copyto!(e.d_tgt, tgt)
    FastMultipole.refresh_cross_producers!(e.ctx, e.d_tgt, pos_d)
    if e.ctx.needs_rebuild
        # box escape: same discipline as the forward pass — rebuild the WHOLE
        # entry around a fresh padded union box and rerun the producers
        delete!(_REVERSE_STATE, objectid(body))
        e = _reverse_entry!(body, tgt, size(pos_d, 2), pos_d, CUDAmod)
        copyto!(e.d_tgt, tgt)
        FastMultipole.refresh_cross_producers!(e.ctx, e.d_tgt, pos_d)
        e.ctx.needs_rebuild && error("reverse-leg union-box rebuild failed " *
            "to contain all bodies (freshly padded box escaped within one step)")
    end
    FastMultipole.refresh_cross_reverse_multipoles!(e.rs, e.ctx, d_pbuf)
    FastMultipole.refresh_cross_reverse_locals!(e.rs, e.ctx)
    FastMultipole.finish_cross_reverse_locals!(e.rs, e.ctx, e.d_tgt)
    FastMultipole.apply_cross_reverse_near!(e.rs, e.ctx, d_pbuf, e.d_tgt;
        sigma_row=8, reg=true)
    # all fallible work done — accumulate LAST
    acc .+= Matrix{Float64}(Array(e.rs.d_out))[2:4, :]
    return nothing
end

"""
    _panel_wake_reverse_device!(out, tgt, tsys, wakes) -> Bool

052h reverse leg: evaluate the wake ParticleFields' velocity at `tsys`'s
target points `tgt` (3 x nt) through the device reverse cross pass,
accumulating into `out[1:3, :]`. Returns false (dense fallback, nothing
touched) when the leg cannot run. With PANEL_WAKE_FMM_XVERIFY=1 additionally
evaluates the dense reference and prints the relU.
"""
function _panel_wake_reverse_device!(out::AbstractMatrix{Float64},
        tgt::AbstractMatrix{Float64}, tsys, wakes::Tuple)
    _panel_wake_fmm_device() ||
        return _gpu_route_fallback!(:wake_to_panels_fmm, "device route disarmed")
    tsys isa AbstractBody ||
        return _gpu_route_fallback!(:wake_to_panels_fmm,
            "reverse leg covers body targets only (got $(typeof(tsys)))")
    FastMultipole.load_cuda_radix_lifecycle!() ||
        return _gpu_route_fallback!(:wake_to_panels_fmm,
            "CUDA radix lifecycle unavailable")
    for w in wakes
        w isa FLOWVPM.ParticleField ||
            return _gpu_route_fallback!(:wake_to_panels_fmm,
                "unsupported wake source $(typeof(w))")
        w.kernel === FLOWVPM.kernel_gaussianerf ||
            return _gpu_route_fallback!(:wake_to_panels_fmm,
                "unsupported particle kernel")
        w.particles isa Array &&
            return _gpu_route_fallback!(:wake_to_panels_fmm,
                "host-backed particle field (device leg needs CuArray particles)")
    end
    _cross_config()   # validate tunables loudly up front
    CUDAmod = getglobal(FastMultipole, :CUDA)
    nt = size(tgt, 2)
    acc = zeros(3, nt)
    for w in wakes
        np = w.np
        np == 0 && continue
        Pd = w.particles
        pos_d, np_cap = _cross_padded_positions!(Pd, np, CUDAmod)
        d_pbuf = _reverse_padded_buffer!(Pd, np, np_cap, pos_d, CUDAmod)
        e = _reverse_entry!(tsys, tgt, np_cap, pos_d, CUDAmod)
        _reverse_run!(e, tsys, tgt, np, pos_d, d_pbuf, acc, CUDAmod)
    end
    if _panel_wake_fmm_xverify()
        ref = zeros(3, nt)
        packed = Any[_gpu_pack_source!(w, nothing) for w in wakes if w.np > 0]
        isempty(packed) ||
            _gpu_direct_batch!(ref, tgt, packed, false, false, _gpu_device())
        relU = sqrt(sum(abs2, acc .- ref)) / max(sqrt(sum(abs2, ref)), eps())
        println("panel_wake_xverify nt=$nt relU=$relU")
        flush(stdout)
    end
    out[1:3, :] .+= acc
    _gpu_route_hit!(:wake_to_panels_fmm, :cuda)
    return true
end

_gpu_element_type(::AbstractBody{E}) where E = E

"""
    _panel_fmm_evaluate!(out, positions, bodies, grad)

052d Phase 1 core: evaluate `bodies`' influence at `positions` (3 x n) through
host `FastMultipole.fmm!` and ACCUMULATE velocity into `out[1:3, :]` and, when
`grad`, the velocity gradient into `out[4:12, :]` in `J_INDEX` linear order
(the same column-major SMatrix order the production fmm! path writes through
`FLOWVPM.buffer_to_target_system!`). The target view is a
`FastMultipole.ProbeSystemArray` — positions in, gradient/hessian out — so no
particle storage is touched here; callers own the accumulation into U/J rows.
"""
function _panel_fmm_evaluate!(out::AbstractMatrix{Float64},
        positions::AbstractMatrix{Float64}, bodies::Tuple, grad::Bool)
    n = size(positions, 2)
    n == 0 && return nothing
    probes = FastMultipole.ProbeSystemArray(n)
    probes.position .= positions
    FastMultipole.fmm!((probes,), bodies;
        expansion_order=_panel_fmm_p(),
        multipole_acceptance=_panel_fmm_theta(),
        leaf_size_source=_panel_fmm_leaf(),
        scalar_potential=false, gradient=true, hessian=grad,
        shrink=true)
    @views out[1:3, :] .+= probes.gradient
    if grad
        @views out[4:12, :] .+= reshape(probes.hessian, 9, n)
    end
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

    # 052d Phase 1: panels -> particles through the host tree-FMM (default
    # off). ONE fmm! call covers ALL AbstractBody sources: positions come down
    # once (D2H, ~6 MB at 242k), results go back up once (H2D) and accumulate
    # on device. Non-body sources stay on the dense device path below.
    fmm_bodies = panel_influence_fmm_enabled() ?
        Tuple(s for s in sources if s isa AbstractBody &&
            _gpu_source_columns(s) > 0) : ()
    # 052d Stage F: device cross pass first; host fmm! only as fallback
    if !isempty(fmm_bodies) && _panel_cross_device!(pfield, fmm_bodies, grad)
        fmm_bodies_handled = fmm_bodies
        fmm_bodies = ()
    else
        fmm_bodies_handled = ()
    end
    if !isempty(fmm_bodies)
        out_h = Matrix{Float64}(undef, 0, 0)
        _step_timer_measure(:rotor_panels_to_particles_fmm; nested=true) do
            CUDAmod = getglobal(FastMultipole, :CUDA)
            Pd = pfield.particles
            pos_h = Matrix{Float64}(Array(Pd[1:3, 1:np]))
            nout_f = grad ? 12 : 3
            out_h = zeros(Float64, nout_f, np)
            _panel_fmm_evaluate!(out_h, pos_h, fmm_bodies, grad)
            outf_d = Base.invokelatest(CUDAmod.CuArray,
                Matrix{eltype(Pd)}(out_h))
            view(Pd, FLOWVPM.U_INDEX, 1:np) .+= view(outf_d, 1:3, :)
            grad && (view(Pd, FLOWVPM.J_INDEX, 1:np) .+= view(outf_d, 4:12, :))
        end
        # PANEL_INFLUENCE_FMM_SNAPSHOT=1 (probe diagnostics, NOT production):
        # additionally evaluate the same leg with the dense device kernels
        # into a scratch buffer (NOT accumulated) and print the rel-RMS
        # FMM-vs-dense field error at production shape. Costs one dense leg
        # per call — probe runs only.
        if lowercase(get(ENV, "PANEL_INFLUENCE_FMM_SNAPSHOT", "0")) in
                ("1", "true", "on", "yes")
            CUDAmod = getglobal(FastMultipole, :CUDA)
            Pd = pfield.particles
            nout_f = grad ? 12 : 3
            ref_d = Base.invokelatest(CUDAmod.zeros, eltype(Pd), nout_f, np)
            tgts_d = Pd[1:3, 1:np]
            for ssys in fmm_bodies
                srcmat, kern = _gpu_pack_source!(ssys, nothing)
                src_d = srcmat isa Matrix ?
                    Base.invokelatest(CUDAmod.CuArray, srcmat) : srcmat
                Base.invokelatest(FastMultipole.direct_rectangular!, ref_d,
                    tgts_d, kern, src_d; gradient=grad)
            end
            ref_h = Matrix{Float64}(Array(ref_d))
            relU = sqrt(sum(abs2, view(ref_h, 1:3, :) .- view(out_h, 1:3, :))) /
                max(sqrt(sum(abs2, view(ref_h, 1:3, :))), eps())
            msg = "panel_fmm_snapshot p=$(_panel_fmm_p()) np=$np relU=$relU"
            if grad
                relJ = sqrt(sum(abs2, view(ref_h, 4:12, :) .- view(out_h, 4:12, :))) /
                    max(sqrt(sum(abs2, view(ref_h, 4:12, :))), eps())
                msg *= " relJ=$relJ"
            end
            println(msg); flush(stdout)
        end
    end

    packed = Any[]
    packed_sys = Any[]
    for ssys in sources
        ssys === pfield && continue
        any(b -> b === ssys, fmm_bodies) && continue
        any(b -> b === ssys, fmm_bodies_handled) && continue
        _gpu_source_columns(ssys) == 0 && continue
        push!(packed, _gpu_pack_source!(ssys, nothing))
        push!(packed_sys, ssys)
    end
    isempty(packed) && return nothing

    CUDAmod = getglobal(FastMultipole, :CUDA)
    P = pfield.particles
    tgt_d = P[1:3, 1:np]                       # device position slice
    nout = grad ? 12 : 3
    out_d = Base.invokelatest(CUDAmod.zeros, eltype(P), nout, np)
    for ((srcmat, kern), ssys) in zip(packed, packed_sys)
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
