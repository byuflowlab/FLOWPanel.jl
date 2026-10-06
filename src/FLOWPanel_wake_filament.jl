#=##############################################################################################
# DESCRIPTION
#   FilamentParticleWake (BRAINSTORM 018 r4): a near wake of trailing-only
#   (streamwise) vortex filaments held for a fixed angular extent, then
#   converted to vortex particles. Parallel to `PanelParticleWake`, with a
#   `TrailingFilamentSheet` replacing the inner `PanelWake`.
#
#   The sheet reuses PanelWake's EXACT storage layout (nodes/strength/velocity
#   arrays, nwakes/overflowed refs), so the FIFO shift, probe maps, propagate!,
#   apply_freestream!, .vts VTK format, and warmstart loader are shared with
#   PanelWake by signature widening (see AbstractWakeSheet), not duplicated.
#
# AUTHORSHIP
#   Created   : Oct 2026
=##############################################################################################

#------- TrailingFilamentSheet -------#

"""
    TrailingFilamentSheet(shedding, TF=Float64; nwakerows, core_size=1e-3,
        filament_core_size=nothing, shed_with_induced_velocity=true,
        freestream_convection=false)
    TrailingFilamentSheet(body; nwakerows, kwargs...)

Near-wake sheet of trailing-only (streamwise) vortex filaments
(BRAINSTORM 018 r4). Storage is identical to [`PanelWake`](@ref) —
`nodes[i_surf]` is `3 × (nwakerows+1) × (ncols+1)`, `strength[i_surf]` is
`1 × (nwakerows+1) × ncols` and stores the column circulation `Γ` (NOT the
trailing jump `ΔΓ`) — but the sheet *evaluates* as trailing filaments only:
the streamwise segment at station `j` of row `irow` runs from node
`(irow, j)` to `(irow+1, j)` and carries the spanwise circulation jump
[`_trailing_filament_strength`](@ref), the byte-for-byte
`LegacyEdgeJumpConversion` trailing decomposition. Row 1 is pinned to the
TE+Das line by `update_TE!`; older rows convect with the induced velocity
(or freestream, per the flags shared with `PanelWake`).

`filament_core_size[i_surf][j]` is the regularization core size `σ_j` of every
trailing filament at wake node column `j` (length `ncols + 1`, always
populated; `Float64` regardless of `TF` per the AD core-size rule). The
filament family is the global `FILAMENT_REGULARIZATION[]` (linegauss in
production).

!!! warning "Dropped physics (deliberate modeling choice)"
    The held sheet carries **no spanwise vorticity**: under unsteady loading
    the trailing-only representation is divergence-inconsistent (the spanwise
    ΔΓ/Δt filaments of a panel row are simply absent until conversion — the
    same modeling choice as the production spanwise `NoShed`, applied during
    the hold instead of at conversion). The row-1 TE bound-vortex image of the
    panel sheet also disappears. Under steady loading (ΔΓ = 0 across rows) the
    representation is exact; the steady-ΔΓ unit test is the exactness
    certificate.

Unsupported `PanelWake` machinery (no fields exist): unsteady/final filaments,
convert-at-shed, surface-vorticity handoff, Kutta Route B live rows
(`live_rows`/`live_step_id` are stored only for index-map parity and guarded
to zero).
"""
struct TrailingFilamentSheet{TF} <: AbstractWakeSheet
    nwakes::Array{Int, 0}
    nodes::Vector{Array{TF, 3}}
    strength::Vector{Array{TF, 3}}
    velocity::Vector{Array{TF, 3}}
    freestream::Vector{TF}
    core_size::Float64
    filament_core_size::Vector{Vector{Float64}}
    overflowed::Array{Bool, 0}
    shed_with_induced_velocity::Bool
    freestream_convection::Bool
    # Per-surface closed-chain flag: true when the shedding edge chain wraps on
    # itself (rotor blade root-tip closure), detected from the pinned row-1
    # node line (see `_refresh_wraps!`) and asserted against the outgoing row
    # at conversion time.
    wraps::Vector{Bool}
    # Mirrored PanelWake live-block metadata for index-map parity; always 0/-1
    # (Kutta Route B is unsupported on this sheet type).
    live_rows::Array{Int, 0}
    live_step_id::Array{Int, 0}
end

function TrailingFilamentSheet(shedding::Vector{Matrix{Int}}, TF=Float64;
        nwakerows=nothing, core_size=1e-3, filament_core_size=nothing,
        shed_with_induced_velocity=true, freestream_convection=false
    )
    nwakerows === nothing && throw(ArgumentError(
        "TrailingFilamentSheet requires an explicit nwakerows"))
    nwakerows >= 1 || throw(ArgumentError(
        "TrailingFilamentSheet requires nwakerows >= 1 (got $(nwakerows)); " *
        "there is no convert-at-shed mode for the filament sheet"))

    # per-station core sizes: always populated (scalar fallback = core_size)
    n_station_cols = [size(s, 2) + 1 for s in shedding]
    if filament_core_size === nothing
        filament_core_size = [fill(Float64(core_size), n) for n in n_station_cols]
    else
        filament_core_size = [Float64.(v) for v in filament_core_size]
    end
    length(filament_core_size) == length(shedding) || throw(ArgumentError(
        "filament_core_size needs one σ vector per shedding surface " *
        "(got $(length(filament_core_size)), expected $(length(shedding)))"))
    for (k, sig) in enumerate(filament_core_size)
        length(sig) == n_station_cols[k] || throw(ArgumentError(
            "filament_core_size for surface $(k) has length $(length(sig)), " *
            "expected n_cols + 1 = $(n_station_cols[k])"))
        all(s -> isfinite(s) && s > 0, sig) || throw(ArgumentError(
            "filament_core_size for surface $(k) must be finite and positive"))
    end

    nwakes = Array{Int,0}(undef)
    nwakes[] = 0

    nodes = [zeros(TF, 3, nwakerows+1, size(s, 2)+1) for s in shedding]
    strength = [zeros(TF, 1, nwakerows+1, size(s, 2)) for s in shedding]
    velocity = [zeros(TF, size(n)) for n in nodes]
    freestream = zeros(TF, 3)

    overflowed = Array{Bool,0}(undef)
    overflowed[] = false

    # all-zero nodes degenerate to "wraps"; the flag is refreshed from real
    # geometry at the first update_TE! (before any solve or conversion)
    wraps = [true for _ in shedding]

    live_rows = Array{Int,0}(undef)
    live_rows[] = 0
    live_step_id = Array{Int,0}(undef)
    live_step_id[] = -1

    return TrailingFilamentSheet{TF}(
        nwakes, nodes, strength, velocity, freestream, Float64(core_size),
        filament_core_size, overflowed,
        Bool(shed_with_induced_velocity), Bool(freestream_convection),
        wraps, live_rows, live_step_id,
    )
end

TrailingFilamentSheet(body::AbstractLiftingBody{TK,NK,TF}; kwargs...) where {TK,NK,TF} =
    TrailingFilamentSheet(body.shedding, TF; kwargs...)

get_probes(sheet::TrailingFilamentSheet) = (ProbeWrapper(sheet),)

# The sheet IS its own (complete) FMM source: the trailing decomposition needs
# no final-filament companion system.
get_sources(sheet::TrailingFilamentSheet) = (sheet,)

_n_wake_source_rows(sheet::TrailingFilamentSheet) = sheet.nwakes[] - sheet.live_rows[]

_logical_nwakerows(sheet::TrailingFilamentSheet) = size(sheet.nodes[1], 2) - 1

"Number of trailing-filament stations of shedding surface `i_surf`: one per
wake node column, except a wrapped (closed) chain where column `ncols + 1`
coincides with column 1."
_n_stations(sheet::TrailingFilamentSheet, i_surf) =
    size(sheet.strength[i_surf], 3) + (sheet.wraps[i_surf] ? 0 : 1)

"Geometric closed-chain test on node row `irow` of a wake node array — the
same `norm(r1 - rend) < 5eps()` test `LegacyEdgeJumpConversion` applies to
the outgoing row of a `PanelWake`."
@inline function _wraps_row_test(nodes::Array{<:Any,3}, irow)
    n_node_cols = size(nodes, 3)
    r1 = SVector{3}(nodes[1, irow, 1], nodes[2, irow, 1], nodes[3, irow, 1])
    rend = SVector{3}(nodes[1, irow, n_node_cols], nodes[2, irow, n_node_cols],
        nodes[3, irow, n_node_cols])
    return norm(r1 - rend) < 5 * eps()
end

function _refresh_wraps!(sheet::TrailingFilamentSheet, system)
    for i_surf in eachindex(sheet.nodes)
        sheet.wraps[i_surf] = _wraps_row_test(sheet.nodes[i_surf], 1)
    end
    return nothing
end

"""
    _trailing_filament_strength(sheet, i_surf, irow, j)

Circulation of the trailing (streamwise) filament at station `j` of row
`irow`: the spanwise jump of the stored column circulations, byte-for-byte the
`LegacyEdgeJumpConversion` trailing decomposition. Wrapped chain (`wraps`):
station 1 carries `Γ(1) - Γ(ncols)` and stations run `1:ncols`; open chain:
station 1 carries `Γ(1)`, station `ncols + 1` carries `-Γ(ncols)`.

This single function feeds BOTH the FMM source buffer and the
particle conversion, so the held filaments and the particles they become can
never drift apart.
"""
@inline function _trailing_filament_strength(sheet::TrailingFilamentSheet, i_surf, irow, j)
    strength = sheet.strength[i_surf]
    ncols = size(strength, 3)
    if sheet.wraps[i_surf]
        jm1 = j == 1 ? ncols : j - 1
        return strength[1, irow, j] - strength[1, irow, jm1]
    else
        j == 1 && return strength[1, irow, 1]
        j == ncols + 1 && return -strength[1, irow, ncols]
        return strength[1, irow, j] - strength[1, irow, j-1]
    end
end

#------- FastMultipole source interface -------#
# Buffer layout is IDENTICAL to FilamentWrapper (12 slots: 1:3 midpoint,
# 4 radius, 5 strength, 6:8 v1, 9:11 v2, 12 per-segment core size), so the
# sheet reuses `_direct_filaments!` and `_filament_body_to_multipole!`.

FastMultipole.numtype(::TrailingFilamentSheet{TF}) where {TF} = TF

FastMultipole.data_per_body(::TrailingFilamentSheet) = 12

FastMultipole.strength_dims(::TrailingFilamentSheet) = 1

# This is what makes `_filter_scalar_potential_sources` exclude the sheet from
# scalar-potential evaluations automatically.
FastMultipole.has_vector_potential(::TrailingFilamentSheet) = true

FastMultipole.get_n_bodies(sheet::TrailingFilamentSheet) =
    _n_wake_source_rows(sheet) * sum(_n_stations(sheet, i) for i in eachindex(sheet.strength))

function global_to_matrix_index(sheet::TrailingFilamentSheet, i_body)
    # which shedding surface (source rows only; mirrors PanelWake's live-row
    # exclusion arithmetic even though live_rows is guarded to 0 here)
    nrows = _n_wake_source_rows(sheet)
    isurf = 1
    i_local = i_body
    n = 0
    for i in eachindex(sheet.strength)
        n += _n_stations(sheet, i) * nrows
        if i_body <= n
            break
        end
        isurf += 1
        i_local -= _n_stations(sheet, i) * nrows
    end

    # local index -> (station, row); row-fastest within a station column,
    # matching PanelWake's ordering
    j, irow = divrem(i_local - 1, nrows)
    j += 1
    irow += 1
    irow += sheet.live_rows[] # skip the reserved live block (row 1 is newest)

    return isurf, irow, j
end

function matrix_to_global_index(sheet::TrailingFilamentSheet, isurf, irow, j)
    nrows = _n_wake_source_rows(sheet)
    i_body = (j - 1) * nrows + (irow - sheet.live_rows[])
    for i in 1:(isurf-1)
        i_body += _n_stations(sheet, i) * nrows
    end
    return i_body
end

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, sheet::TrailingFilamentSheet, i_body)
    isurf, irow, j = global_to_matrix_index(sheet, i_body)
    nodes = sheet.nodes[isurf]

    # streamwise segment at station node column j spanning rows irow -> irow+1
    v1 = SVector{3}(nodes[1, irow, j], nodes[2, irow, j], nodes[3, irow, j])
    v2 = SVector{3}(nodes[1, irow+1, j], nodes[2, irow+1, j], nodes[3, irow+1, j])
    sigma = sheet.filament_core_size[isurf][j]

    buffer[1:3, i_buffer] .= 0.5 * (v1 + v2)
    # regularization reach, not just the core itself (FilamentWrapper precedent)
    buffer[4, i_buffer] = 0.5 * norm(v2 - v1) +
        radius_inflation(VortexRing, sigma, fmm_radius_tolerance(sheet))
    buffer[5, i_buffer] = _trailing_filament_strength(sheet, isurf, irow, j)
    buffer[6:8, i_buffer] .= v1
    buffer[9:11, i_buffer] .= v2
    buffer[12, i_buffer] = sigma
end

function FastMultipole.get_position(sheet::TrailingFilamentSheet, i)
    isurf, irow, j = global_to_matrix_index(sheet, i)
    nodes = sheet.nodes[isurf]
    v1 = SVector{3}(nodes[1, irow, j], nodes[2, irow, j], nodes[3, irow, j])
    v2 = SVector{3}(nodes[1, irow+1, j], nodes[2, irow+1, j], nodes[3, irow+1, j])
    return 0.5 * (v1 + v2)
end

function FastMultipole.body_to_multipole!(sheet::TrailingFilamentSheet, multipole_coefficients, buffer::Matrix, center, bodies_index, harmonics, expansion_order)
    _filament_body_to_multipole!(sheet, multipole_coefficients, buffer, center,
        bodies_index, harmonics, expansion_order)
end

# function barrier: family in the type domain inside the loop (BRAINSTORM 025)
function FastMultipole.direct!(target_system, target_index, switch::FastMultipole.DerivativesSwitch, source_system::TrailingFilamentSheet, source_buffer, source_index)
    _direct_filaments!(target_system, target_index, switch, source_system,
        source_buffer, source_index, Val(FILAMENT_REGULARIZATION[]))
end

function FastMultipole.buffer_to_target_system!(target_system::TrailingFilamentSheet, i_target, ::FastMultipole.DerivativesSwitch, target_buffer, i_buffer)
    @warn "A `::TrailingFilamentSheet` should not be used as a target in an FMM call (its probes are a ProbeWrapper)."
end

#------- VTK output -------#

function write_vtk(name, sheet::TrailingFilamentSheet, idx, t; overwrite=false,
        compress::Bool=true, filament_name=nothing)
    # Same .vts multiblock format as PanelWake (DELIBERATE — keeps the
    # warmstart/replay structured-grid loaders reusable across sheet types).
    _parent, _base = splitdir(name)
    subdir = joinpath(_parent, _base)
    mkpath(subdir)
    block_name = joinpath(subdir, _base)

    vtm = WriteVTK.vtk_multiblock(block_name * ".$idx.vtm")
    if sheet.nwakes[] > 0
        _write_sheet_vts(vtm, sheet, block_name, idx; compress)
    end
    WriteVTK.vtk_save(vtm)
    _pvd_append!(name * ".pvd", t, joinpath(_base, _base * ".$idx.vtm"); overwrite)

    # trailing filaments as polylines (ΔΓ + σ_j cell data) — the evaluated
    # representation, for visual verification against the .vts sheet
    _write_trailing_filaments_vtu(
        isnothing(filament_name) ? name * "_filaments" : filament_name,
        sheet, idx, t; overwrite, compress)
end

function _write_trailing_filaments_vtu(name, sheet::TrailingFilamentSheet, idx, t;
        overwrite=false, compress::Bool=true)
    _parent, _base = splitdir(name)
    subdir = joinpath(_parent, _base)
    mkpath(subdir)
    block_name = joinpath(subdir, _base)

    vtm = WriteVTK.vtk_multiblock(block_name * ".$idx.vtm")
    nrows = _n_wake_source_rows(sheet)
    if nrows > 0
        for i_surf in eachindex(sheet.nodes)
            nodes = sheet.nodes[i_surf]
            n_st = _n_stations(sheet, i_surf)
            n_fils = n_st * nrows

            points = zeros(eltype(nodes), 3, 2 * n_fils)
            cells = Vector{WriteVTK.MeshCell{WriteVTK.VTKCellTypes.VTKCellType, Vector{Int}}}(undef, n_fils)
            strengths = zeros(eltype(sheet.strength[i_surf]), n_fils)
            sigmas = zeros(Float64, n_fils)

            k = 0
            for irow in (1 + sheet.live_rows[]):sheet.nwakes[]
                for j in 1:n_st
                    k += 1
                    ip = 2 * (k - 1)
                    points[:, ip + 1] .= view(nodes, :, irow, j)
                    points[:, ip + 2] .= view(nodes, :, irow + 1, j)
                    cells[k] = WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_LINE, [ip + 1, ip + 2])
                    strengths[k] = _trailing_filament_strength(sheet, i_surf, irow, j)
                    sigmas[k] = sheet.filament_core_size[i_surf][j]
                end
            end

            WriteVTK.vtk_grid(vtm, block_name * ".$(i_surf).$(idx).vtu", points, cells; compress) do vtk
                vtk["strength", WriteVTK.VTKCellData()] = strengths
                vtk["sigma", WriteVTK.VTKCellData()] = sigmas
            end
        end
    end
    WriteVTK.vtk_save(vtm)
    _pvd_append!(name * ".pvd", t, joinpath(_base, _base * ".$idx.vtm"); overwrite)
end

#------- FilamentParticleWake -------#

"""
    FilamentParticleWake(body::AbstractLiftingBody; nwakerows, max_particles=10000,
        core_size=1e-3, filament_core_size=nothing,
        method_trailing=DefaultWakeSheddingMethod(), kwargs...)

Free wake whose near field is a [`TrailingFilamentSheet`](@ref) — trailing-only
vortex filaments held for `nwakerows` steps — and whose far field is a FLOWVPM
particle field; the outgoing filament row is converted to particles by
[`_trailing_filament_strength`](@ref)-exact shedding (BRAINSTORM 018 r4).
Parallel to [`PanelParticleWake`](@ref) with the same keyword names where the
concept carries over.

`nwakerows` is REQUIRED (`>= 1`): the hold extent is a physical modeling
choice (θ* ruling), never a default.

σ-coupling rule: when `method_trailing` resolves stations through a
`StationSigmaOverlap` (possibly wrapped in `OmitStations`) and
`filament_core_size` is not given, the held-filament cores default to the
shedding σ vectors — hold-σ == shed-σ by construction. Passing both is an
error unless they are equal.

Unsupported `PanelParticleWake` machinery (rejected by name): `conversion`
(the trailing decomposition IS the conversion; `SurfaceVorticityConversion`
is PanelParticleWake-only), `method_unsteady` (no spanwise filaments exist to
shed — see the dropped-physics warning on `TrailingFilamentSheet`),
`unsteady_filament`/`include_final_filament` (no such filaments), nwakerows=0
convert-at-shed, and Kutta Route B attachment.
"""
struct FilamentParticleWake{TF,TPF,MT,TPM,TNT} <: AbstractParticleWake
    sheet::TrailingFilamentSheet{TF}
    pfield::TPF                           # FLOWVPM.ParticleField object
    method_trailing::MT                   # particle shedding method
    particle_maintenance::TPM             # particle merge/trim policy chain
    particle_core_size::Float64           # NaN uses source body core_size
    pfield_optargs::TNT                   # resolved FLOWVPM optargs (metadata)
end

_wake_sheet(w::FilamentParticleWake) = w.sheet

# Station-σ extraction for the σ-coupling rule (nothing = not station-based)
_station_sigma_vectors(::WakeSheddingMethod) = nothing
_station_sigma_vectors(m::StationSigmaOverlap) = [Float64.(s) for s in m.sigmas]
_station_sigma_vectors(m::OmitStations) = _station_sigma_vectors(m.method)

function FilamentParticleWake(body::AbstractLiftingBody;
        nwakerows=nothing, max_particles=10000,
        core_size=1e-3,
        filament_core_size=nothing,
        method_trailing::WakeSheddingMethod=DefaultWakeSheddingMethod(),
        particle_maintenance=ParticleMaintenance(),
        particle_core_size::Union{Real,Nothing}=nothing,
        viscous=FLOWVPM.Inviscid(),
        SFS=FLOWVPM.SFS_default,
        relaxation=FLOWVPM.relaxation_correctedpedrizzetti,
        expint=false,
        rk3=false,
        arraytype=Matrix,
        pfield_fmm=FLOWVPM.FMM(autotune_reg_error=false),
        pfield=nothing,
        shed_with_induced_velocity=true,
        freestream_convection=false,
        # named rejections (PanelParticleWake-only concepts; a MethodError
        # here would not say why)
        conversion=nothing,
        method_unsteady=nothing,
    )

    conversion === nothing || throw(ArgumentError(
        "FilamentParticleWake does not take a conversion strategy: the " *
        "trailing-jump decomposition IS the conversion " *
        "(SurfaceVorticityConversion/LegacyEdgeJumpConversion are " *
        "PanelParticleWake-only)"))
    method_unsteady === nothing || throw(ArgumentError(
        "FilamentParticleWake has no unsteady (spanwise) filaments to shed: " *
        "method_unsteady is PanelParticleWake-only (the held sheet drops " *
        "spanwise vorticity by construction; see the TrailingFilamentSheet " *
        "docstring)"))
    nwakerows === nothing && throw(ArgumentError(
        "FilamentParticleWake requires an explicit nwakerows >= 1 (the hold " *
        "extent θ* is a physical modeling choice, never a default)"))
    nwakerows >= 1 || throw(ArgumentError(
        "FilamentParticleWake requires nwakerows >= 1 (got $(nwakerows)); " *
        "there is no convert-at-shed mode for the filament wake"))

    # resolve the sentinel default exactly as the legacy conversion would
    trailing = _resolve_line_policy(LegacyEdgeJumpConversion(), method_trailing,
        "method_trailing")

    # validate station-indexed methods against the body's shedding geometry
    # at construction time (fail here, not mid-conversion)
    for i_surf in eachindex(body.shedding)
        _validate_station_method(trailing, i_surf, size(body.shedding[i_surf], 2) + 1)
    end

    # σ-coupling rule: hold-σ defaults to (and must match) the shedding σ
    station_sigmas = _station_sigma_vectors(trailing)
    if filament_core_size === nothing
        filament_core_size = station_sigmas # nothing when not station-based
    elseif station_sigmas !== nothing
        [Float64.(v) for v in filament_core_size] == station_sigmas || throw(ArgumentError(
            "filament_core_size differs from method_trailing's station σ " *
            "vectors; hold-σ must equal shed-σ (omit filament_core_size to " *
            "couple them by construction)"))
    end

    sheet = TrailingFilamentSheet(body; nwakerows, core_size,
        filament_core_size, shed_with_induced_velocity, freestream_convection)
    TF = FastMultipole.numtype(sheet)

    pfield, pfield_optargs = _make_wake_pfield(pfield, max_particles, TF;
        viscous, pfield_fmm, SFS, expint, rk3, relaxation, arraytype)

    maintenance = ParticleMaintenance(particle_maintenance)
    particle_core_size = particle_core_size === nothing ? NaN : Float64(particle_core_size)
    if !isnan(particle_core_size)
        body.core_size_targets = particle_core_size
    end

    return FilamentParticleWake{TF,typeof(pfield),typeof(trailing),typeof(maintenance),typeof(pfield_optargs)}(
        sheet, pfield, trailing, maintenance, particle_core_size, pfield_optargs,
    )
end

#------- Delegation methods -------#

get_probes(w::FilamentParticleWake) = (get_probes(w.sheet)..., w.pfield)
get_sources(w::FilamentParticleWake) = (get_sources(w.sheet)..., w.pfield)

#------- Shedding and conversion -------#

"""
    _convert_to_particles!(wake::FilamentParticleWake, system=nothing)

Convert the outgoing (oldest) filament row into particles: each trailing
filament at station `j` sheds along its own segment with strength
[`_trailing_filament_strength`](@ref) — the SAME function that fills the FMM
source buffer, so the particles carry exactly the circulation the held
filaments induced (handoff is bookkeeping-exact by construction).

Note: `_shed_particles!` elides exact-zero strengths (NaN guard). This is
exact here too — a zero-strength filament induces nothing.
"""
function _convert_to_particles!(wake::FilamentParticleWake, system=nothing)
    sheet = wake.sheet
    sheet.live_rows[] == 0 || error(
        "TrailingFilamentSheet does not support a reserved live row block " *
        "(Kutta Route B / TEAnchoredAttachment)")
    nwakes = sheet.nwakes[]

    for i_surf in eachindex(sheet.nodes)
        nodes = sheet.nodes[i_surf]
        n_cols = size(sheet.strength[i_surf], 3)
        _validate_station_method(wake.method_trailing, i_surf, n_cols + 1)

        # the cached wraps flag must agree with the geometric test on the
        # outgoing row — a mid-run topology change would silently corrupt the
        # station-1 jump, so fail loudly instead
        wraps_now = _wraps_row_test(nodes, nwakes)
        wraps_now == sheet.wraps[i_surf] || error(
            "TrailingFilamentSheet surface $(i_surf): cached wraps flag " *
            "($(sheet.wraps[i_surf])) disagrees with the geometric test on " *
            "the outgoing row ($(wraps_now)) — shedding-edge chain topology " *
            "changed mid-run?")

        for j in 1:_n_stations(sheet, i_surf)
            Γ = _trailing_filament_strength(sheet, i_surf, nwakes, j)
            r1 = SVector{3}(nodes[1, nwakes, j], nodes[2, nwakes, j], nodes[3, nwakes, j])
            r2 = SVector{3}(nodes[1, nwakes+1, j], nodes[2, nwakes+1, j], nodes[3, nwakes+1, j])
            _shed_particles!(wake.pfield, r1, r2, Γ,
                _station_method(wake.method_trailing, i_surf, j, j))
        end
    end
    return nothing
end

function shed_wake!(wake::FilamentParticleWake, system::AbstractBody)
    sheet = wake.sheet

    n_rows = size(sheet.nodes[1], 2)
    buffer_full = sheet.nwakes[] >= n_rows - 1

    if buffer_full
        # convert the outgoing row to particles before the FIFO shift
        _convert_to_particles!(wake, system)
    end

    # shift filament rows (shared AbstractWakeSheet FIFO; the sheet stores
    # column Γ, so the row-1 strength write is the same mu-jump as PanelWake)
    shed_wake!(sheet, system)
end
