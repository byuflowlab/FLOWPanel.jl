"""
    AbstractFreeWake

Abstract supertype for wake models that exchange influence with
[`AbstractBody`](@ref) systems during time stepping.
"""
abstract type AbstractFreeWake end

"""
    ProbeWrapper(system)

Wrapper that exposes wake probe points to the influence backends without
changing the underlying wake storage layout.
"""
struct ProbeWrapper{P}
    system::P
end

"""
    solve!(body, wake, uinf, t=0.0; body_solver=Backslash(body), backend=FastMultipoleBackend(...))

Perform one coupled body-wake solve using the supplied freestream velocity
vector `uinf`.
"""
function solve!(body::AbstractBody, wake::AbstractFreeWake, uinf::AbstractArray, t=0.0;
        body_solver::AbstractSolver=Backslash(body),
        backend=FastMultipoleBackend(;
            expansion_order=10,
            multipole_acceptance=0.4,
            leaf_size=100,
            shrink=true,
            recenter=false,
        )
    )

    # reset potential/velocity
    reset!(wake)
    reset!(body)
    
    # get probes
    wake_probes = get_probes(wake)

    # wake-on-all velocity
    influence!((body, wake_probes), (wake,), backend;
        velocity=true,
        velocity_gradient=(requires_hessian(body), requires_hessian(wake)))

    # freestream
    apply_freestream!(body, uinf)
    apply_freestream!(wake, uinf)

    # kinematics
    kinematic_velocity!(body, frames; skip_top_level=false)

    # solve body with the panel/solve regularization
    _set_core_sizes!((body,), :core_size_panel)
    solve!(body, body_solver; backend)

    # body-on-all influence with the external-target regularization
    _set_core_sizes!((body,), :core_size_targets)
    influence!((body, wake_probes), (body,), backend;
        velocity=true,
        velocity_gradient=(requires_hessian(body), requires_hessian(wake)),
        direct_conditioning=_self_panel_core_size_conditioning())

    return nothing
end

requires_hessian(b::AbstractBody) = b.needs_velocity_gradient[]
requires_hessian(::AbstractFreeWake) = false # default behavior
requires_hessian(pw::ProbeWrapper) = requires_hessian(pw.system)

#--- Panel Wake ---#

"""
    PanelWake(shedding, kernel, TF=Float64; core_size=1e-3, nwakerows=100,
        shed_with_induced_velocity=true, unsteady_filament=true,
        include_final_filament=true, freestream_convection=false)
    PanelWake(body; kernel=get_wake_kernel(body), nwakerows=100,
        shed_with_induced_velocity=true, unsteady_filament=true,
        include_final_filament=true, freestream_convection=false)

Wake model that stores a panelized wake sheet behind one or more shedding-edge
chains. Set `shed_with_induced_velocity=false` to convect the first wake row
with freestream only when forming newly shed panels. Set
`unsteady_filament=false` to make the final-edge filament cancel the current
last wake row instead of representing the shifted-out previous row.
Set `include_final_filament=false` for a strictly finite, panel-only wake whose
sources all expose scalar potential (for example a ConstantDoublet pressure
oracle); wake-length convergence must then be checked explicitly.
Set `freestream_convection=true` to convect *every* wake row with the freestream
only (no rollup), so the sheet stays straight along `U∞`; this reproduces the
geometry of the semi-infinite wake and exists to test consistency between the
two wake representations. It overrides `shed_with_induced_velocity`.

`convert_at_shed` (BRAINSTORM 024) marks the wake as the storage backing of a
`nwakerows = 0` [`PanelParticleWake`](@ref): single-row (N=1) storage whose
just-shed row is converted to particles inside `shed_wake!` itself, after
which `nwakes[]` is reset to zero so no free sheet ever enters a solve. It is
an internal flag — construct it through `PanelParticleWake(...; nwakerows=0)`,
not directly. A standalone `PanelWake` requires `nwakerows >= 1`.
"""
struct PanelWake{TK,NK,TF} <: AbstractFreeWake
    nwakes::Array{Int, 0}
    nodes::Vector{Array{TF, 3}}
    strength::Vector{Array{TF, 3}}
    velocity::Vector{Array{TF, 3}}
    freestream::Vector{TF}
    core_size::Float64
    overflowed::Array{Bool, 0}
    shed_with_induced_velocity::Bool
    unsteady_filament::Bool
    include_final_filament::Bool
    freestream_convection::Bool
    # Convert-at-shed marker (BRAINSTORM 024): true only for the N=1 storage
    # backing a `nwakerows = 0` PanelParticleWake. `nwakes[]` is 0 whenever a
    # solve sees this wake; the final filament then sits on node row 1 (the
    # TE+Das line) carrying the full just-shed strength.
    convert_at_shed::Bool
    # Live-block reservation metadata (BRAINSTORM 015 Route B / TEAnchoredAttachment).
    # live_rows[] newest panel rows are the reserved live block: their strengths
    # are owned by the body-side attachment operator during the coupled solve and
    # they are excluded from old-wake source views (see get_n_bodies and the
    # index maps). 0 on the legacy path — all arithmetic then reduces exactly to
    # the pre-existing behavior. live_step_id[] is the physical-step identifier
    # of the current live block (-1 = none).
    live_rows::Array{Int, 0}
    live_step_id::Array{Int, 0}
    # Panel-to-particle handoff state (BRAINSTORM 016). False until a
    # `SurfaceVorticityConversion` transaction first commits; from then on the
    # retained final panel filament cancels the *current* final active row's
    # downstream edge (see `_final_filament_strength`), because that edge's
    # circulation now lives in the deposited area particles. Legacy edge-jump
    # wakes never set it, so their filament bookkeeping is untouched.
    particle_handoff_active::Array{Bool, 0}
    # Fraction of the panel/particle interface circulation that the conversion
    # already deposited (`alpha` of BRAINSTORM 016's streamwise attribution).
    # 1.0 = the whole upstream face went to particles, so the filament cancels
    # the interface completely; 0.0 = none did, which reproduces the legacy
    # unsteady filament exactly. Only consulted when the handoff is active.
    particle_handoff_weight::Array{Float64, 0}
end

"""
    get_probes(wake::PanelWake)

Return the wake probe systems that should receive wake-induced influence.
"""
function get_probes(wake::PanelWake)
    return (ProbeWrapper(wake),)
end

"""
    get_sources(wake::PanelWake)

Return the wake source systems used by the active influence backend.
"""
function get_sources(wake::PanelWake)
    return wake.include_final_filament ?
        (wake, FilamentWrapper(wake)) : (wake,)
end

function PanelWake(shedding::Vector{Matrix{Int}}, kernel, TF=Float64;
        core_size=1e-3, nwakerows=100, shed_with_induced_velocity=true,
        unsteady_filament=true, include_final_filament=true,
        freestream_convection=false, convert_at_shed=false
    )
    nwakerows >= 1 || throw(ArgumentError(
        "PanelWake requires nwakerows >= 1 (got $(nwakerows)); " *
        "nwakerows = 0 (convert-at-shed) is a PanelParticleWake mode"))

    # nwakes
    nwakes = Array{Int,0}(undef)
    nwakes[] = 0
    
    # nodes
    nodes = [zeros(TF, 3, nwakerows+1, size(s, 2)+1) for s in shedding]
    
    # strength
    dim = kernel_dim(kernel)
    strength = [zeros(TF, dim, nwakerows+1, size(s, 2)) for s in shedding]
    
    # velocity
    velocity = [zeros(TF, size(n)) for n in nodes]

    # latest freestream applied to this wake
    freestream = zeros(TF, 3)
    
    # overflowed
    overflowed = Array{Bool,0}(undef)
    overflowed[] = false

    # live-block reservation metadata (0/-1 = legacy, no live block)
    live_rows = Array{Int,0}(undef)
    live_rows[] = 0
    live_step_id = Array{Int,0}(undef)
    live_step_id[] = -1

    # panel-to-particle handoff state (BRAINSTORM 016); false = legacy behavior
    particle_handoff_active = Array{Bool,0}(undef)
    particle_handoff_active[] = false
    particle_handoff_weight = Array{Float64,0}(undef)
    particle_handoff_weight[] = 1.0

    return PanelWake{kernel, dim, TF}(
        nwakes, nodes, strength, velocity, freestream, core_size, overflowed,
        Bool(shed_with_induced_velocity), Bool(unsteady_filament),
        Bool(include_final_filament), Bool(freestream_convection),
        Bool(convert_at_shed),
        live_rows, live_step_id, particle_handoff_active, particle_handoff_weight,
    )
end

PanelWake(body::AbstractLiftingBody{TK,NK,TF}, kernel=get_wake_kernel(body); nwakerows=100, kwargs...) where {TK,NK,TF} =
    PanelWake(body.shedding, kernel, TF; nwakerows, kwargs...)

function reset!(wake::PanelWake)
    for vel in wake.velocity
        vel .= zero(eltype(vel))
    end
end

function apply_freestream!(wake::PanelWake, uinf)
    wake.freestream .= uinf
    for vel in wake.velocity
        for ns in axes(vel, 3)
            for nc in 1:wake.nwakes[]+1
                vel[:, nc, ns] .+= uinf
            end
        end
    end
end

# FastMultipole compatibility

"Number of panel rows exposed as old-wake sources: the newest `live_rows[]`
rows are the reserved live block owned by the body-side attachment operator
(BRAINSTORM 015 Route B) and are excluded to prevent double counting. 0 live
rows (the legacy default) reduces every source-view expression to the
pre-existing arithmetic."
_n_wake_source_rows(wake::PanelWake) = wake.nwakes[] - wake.live_rows[]

"Logical `nwakerows` of a wake as configured by the user: 0 for a
convert-at-shed wake (BRAINSTORM 024), whose N=1 storage is an implementation
detail, else the storage row count. This is what metadata/manifests record so
reconstruction reproduces the mode."
_logical_nwakerows(wake::PanelWake) =
    wake.convert_at_shed ? 0 : size(wake.nodes[1], 2) - 1

function global_to_matrix_index(wake::PanelWake, i_wake)

    # determine which shedding surface we're on (old-wake source rows only;
    # the reserved live block is excluded from the source view)
    nwakes = _n_wake_source_rows(wake)
    isurf = 1
    i_wake_local = i_wake
    npanels = 0
    for i in eachindex(wake.strength)
        npanels += size(wake.strength[i], 3) * nwakes
        if i_wake <= npanels
            break
        end
        isurf += 1 # advance to the next surface
        i_wake_local -= size(wake.strength[i], 3) * nwakes # adjust local index
    end

    # convert local index to matrix indices
    icol, irow = divrem(i_wake_local - 1, nwakes)
    icol += 1 # adjust for 1-based indexing
    irow += 1 # adjust for 1-based indexing
    irow += wake.live_rows[] # skip the reserved live block (row 1 is newest)

    return isurf, irow, icol
end

function global_to_matrix_index(wake::ProbeWrapper{<:PanelWake}, i_wake)

    # determine which shedding surface we're on
    nrows = wake.system.nwakes[] + 1
    isurf = 1
    i_wake_local = i_wake
    nnodes = 0
    for i in eachindex(wake.system.nodes)
        nnodes += size(wake.system.nodes[i], 3) * nrows
        if i_wake <= nnodes
            break
        end
        isurf += 1 # advance to the next surface
        i_wake_local -= size(wake.system.nodes[i], 3) * nrows # adjust local index
    end

    # convert local index to matrix indices
    icol, irow = divrem(i_wake_local - 1, nrows)
    icol += 1 # adjust for 1-based indexing
    irow += 1 # adjust for 1-based indexing

    return isurf, irow, icol
end

function matrix_to_global_index(wake::PanelWake, isurf, irow, icol)
    # convert matrix indices to local index (inverse of the source view above:
    # absolute row numbering, with the reserved live block excluded)
    nwakes = _n_wake_source_rows(wake)
    i_wake = (icol - 1) * nwakes + (irow - wake.live_rows[])

    # account for previous surfaces
    for i in 1:(isurf-1)
        i_wake += size(wake.strength[i], 3) * nwakes
    end

    return i_wake
end

function matrix_to_global_index(wake::ProbeWrapper{<:PanelWake}, isurf, irow, icol)
    # convert matrix indices to local index
    i_wake = (icol - 1) * (wake.system.nwakes[] + 1) + irow

    # account for previous surfaces
    for i in 1:(isurf-1)
        i_wake += size(wake.system.nodes[i], 3) * (wake.system.nwakes[] + 1)
    end

    return i_wake
end

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, system::PanelWake{TK}, i_body) where TK

    # get surface index of global `i_body` index
    isurf, irow, icol = global_to_matrix_index(system, i_body)

    # get nodes
    v1x, v1y, v1z = view(system.nodes[isurf], :, irow, icol)
    v2x, v2y, v2z = view(system.nodes[isurf], :, irow+1, icol)
    v3x, v3y, v3z = view(system.nodes[isurf], :, irow+1, icol+1)
    v4x, v4y, v4z = view(system.nodes[isurf], :, irow, icol+1)

    # get centroid
    buffer[1, i_buffer] = (v1x + v2x + v3x + v4x) * 0.25
    buffer[2, i_buffer] = (v1y + v2y + v3y + v4y) * 0.25
    buffer[3, i_buffer] = (v1z + v2z + v3z + v4z) * 0.25

    # get radius
    r1x = v3x - v1x
    r1y = v3y - v1y
    r1z = v3z - v1z
    r1 = sqrt(r1x*r1x + r1y*r1y + r1z*r1z) * 0.5
    r2x = v4x - v2x
    r2y = v4y - v2y
    r2z = v4z - v2z
    r2 = sqrt(r2x*r2x + r2y*r2y + r2z*r2z) * 0.5
    # + regularization reach: the wake panel kernel is regularized over
    # `core_size` (see `FastMultipole.direct!` for `PanelWake`), which the
    # multipole expansions cannot represent
    buffer[4, i_buffer] = max(r1, r2) +
        radius_inflation(TK, system.core_size, fmm_radius_tolerance(system))

    # get strength
    s = view(system.strength[isurf], :, irow, icol)
    nk = length(s)
    for i in 1:nk
        buffer[4+i, i_buffer] = s[i]
    end

    # save vertices
    vs = (v1x, v1y, v1z,
          v2x, v2y, v2z,
          v3x, v3y, v3z,
          v4x, v4y, v4z)
    for iv in 1:4
        for id in 1:3
            buffer[4+nk+(iv-1)*3+id, i_buffer] = vs[(iv-1)*3+id]
        end
    end
end

FastMultipole.numtype(system::PanelWake{TK,NK,TF}) where {TK,NK,TF} = TF

function FastMultipole.numtype(system::ProbeWrapper)
    return FastMultipole.numtype(system.system)
end

FastMultipole.data_per_body(system::PanelWake) = 4 + size(system.strength[1], 1) + 12

FastMultipole.has_vector_potential(::PanelWake{TK,NK,TF}) where {TK,NK,TF} = false # TK<:Union{VortexRing, ConstantVortexSheet}

function FastMultipole.get_position(system::PanelWake, i)

    # get surface index of global `i` index
    isurf, irow, icol = global_to_matrix_index(system, i)

    # get nodes
    v1x, v1y, v1z = view(system.nodes[isurf], :, irow, icol)
    v2x, v2y, v2z = view(system.nodes[isurf], :, irow+1, icol)
    v3x, v3y, v3z = view(system.nodes[isurf], :, irow+1, icol+1)
    v4x, v4y, v4z = view(system.nodes[isurf], :, irow, icol+1)

    # get centroid
    cx = (v1x + v2x + v3x + v4x) * 0.25
    cy = (v1y + v2y + v3y + v4y) * 0.25
    cz = (v1z + v2z + v3z + v4z) * 0.25

    return FastMultipole.StaticArrays.SVector{3}(cx, cy, cz)
end

function FastMultipole.get_position(system::ProbeWrapper{<:PanelWake}, i)

    # get surface index of global `i` index
    isurf, irow, icol = global_to_matrix_index(system, i)

    # get nodes
    v1x, v1y, v1z = view(system.system.nodes[isurf], :, irow, icol)

    return FastMultipole.StaticArrays.SVector{3}(v1x, v1y, v1z)
end

FastMultipole.strength_dims(system::PanelWake) = size(system.strength[1], 1)

FastMultipole.get_n_bodies(system::PanelWake) = _n_wake_source_rows(system) * sum(size(s, 3) for s in system.strength)

FastMultipole.get_n_bodies(system::ProbeWrapper{<:PanelWake}) = (system.system.nwakes[]+1) * sum(size(s, 3) for s in system.system.nodes)

FastMultipole.metadata_per_body(system::ProbeWrapper{<:PanelWake}) = 2
FastMultipole.previous_potential_metadata_index(system::ProbeWrapper{<:PanelWake}) = 1
FastMultipole.previous_gradient_metadata_index(system::ProbeWrapper{<:PanelWake}) = 2

function FastMultipole.metadata_to_buffer!(buffer, switch, i_buffer, system::ProbeWrapper{<:PanelWake}, i_body)
    isurf, irow, icol = global_to_matrix_index(system, i_body)
    vx = system.system.velocity[isurf][1, irow, icol]
    vy = system.system.velocity[isurf][2, irow, icol]
    vz = system.system.velocity[isurf][3, irow, icol]
    buffer[FastMultipole.metadata_index(switch, 1), i_buffer] = zero(eltype(buffer))
    buffer[FastMultipole.metadata_index(switch, 2), i_buffer] = sqrt(vx*vx + vy*vy + vz*vz)
    return nothing
end

function FastMultipole.buffer_to_target_system!(target_system::ProbeWrapper{<:PanelWake}, i_target, switch::FastMultipole.DerivativesSwitch{PS,VS,GS,NO,NM}, target_buffer, i_buffer) where {PS,VS,GS,NO,NM}
    # get surface index of global `i_target` index
    isurf, irow, icol = global_to_matrix_index(target_system, i_target)

    # save potential
    # if PS
    #     target_system.potential[isurf][irow, icol] = target_buffer[i_buffer]
    # end

    # save velocity
    if VS
        vx, vy, vz = FastMultipole.get_gradient(target_buffer, switch, i_buffer)
        target_system.system.velocity[isurf][1, irow, icol] += vx
        target_system.system.velocity[isurf][2, irow, icol] += vy
        target_system.system.velocity[isurf][3, irow, icol] += vz
    end

    # save Hessian (not currently used for PanelWake)
    if GS
        # @warn "Hessian output not currently implemented for PanelWake targets"
    end
end

function rotate_to_panel(v1::FastMultipole.SVector{3,TF}, v2, v3) where {TF}
    # explicit cross(v2-v1, v3-v1)
    dx1 = v2[1] - v1[1]; dy1 = v2[2] - v1[2]; dz1 = v2[3] - v1[3]
    dx2 = v3[1] - v1[1]; dy2 = v3[2] - v1[2]; dz2 = v3[3] - v1[3]
    new_z_x = dy1 * dz2 - dz1 * dy2
    new_z_y = dz1 * dx2 - dx1 * dz2
    new_z_z = dx1 * dy2 - dy1 * dx2
    new_z_norm = sqrt(new_z_x*new_z_x + new_z_y*new_z_y + new_z_z*new_z_z)
    new_z_x /= new_z_norm
    new_z_y /= new_z_norm
    new_z_z /= new_z_norm

    # new x axis: (v3 - v1) / norm(v3-v1)
    new_x_x = dx2
    new_x_y = dy2
    new_x_z = dz2
    new_x_norm = sqrt(new_x_x*new_x_x + new_x_y*new_x_y + new_x_z*new_x_z)
    new_x_x /= new_x_norm
    new_x_y /= new_x_norm
    new_x_z /= new_x_norm

    # new y axis: (orthogonal to x and z)
    new_y_x = new_z_y*new_x_z - new_z_z*new_x_y
    new_y_y = new_z_z*new_x_x - new_z_x*new_x_z
    new_y_z = new_z_x*new_x_y - new_z_y*new_x_x
    R = FastMultipole.StaticArrays.SMatrix{3,3,TF,9}(new_x_x, new_x_y, new_x_z, 
                                                     new_y_x, new_y_y, new_y_z, 
                                                     new_z_x, new_z_y, new_z_z)

    return R
end

function induced(target::AbstractVector{TF}, source_system::PanelWake{TK,NK,<:Any}, source_buffer::Matrix, i_source, derivatives_switch=FastMultipole.DerivativesSwitch(false,true,false),
        fam::Val=Val(FILAMENT_REGULARIZATION[])) where {TF,TK,NK}

    # get vertices
    v1 = FastMultipole.SVector{3}(view(source_buffer, 4+NK+1:4+NK+3, i_source))
    v2 = FastMultipole.SVector{3}(view(source_buffer, 4+NK+4:4+NK+6, i_source))
    v3 = FastMultipole.SVector{3}(view(source_buffer, 4+NK+7:4+NK+9, i_source))
    v4 = FastMultipole.SVector{3}(view(source_buffer, 4+NK+10:4+NK+12, i_source))

    #--- first triangle ---#

    R = rotate_to_panel(v1, v2, v3)

    # get control point and strength
    control_point = (v1 + v2 + v3) * 0.3333333333333333
    strength = FastMultipole.StaticArrays.SVector{NK,TF}(view(source_buffer, 5:4+NK, i_source))

    # evaluate influence
    core_size = source_system.core_size
    potential, velocity, velocity_gradient = _induced(target, (v1, v2, v3), control_point, strength, TK, core_size, R, derivatives_switch, fam)

    #--- second triangle ---#

    R = rotate_to_panel(v1, v3, v4)

    # get control point and strength
    control_point = (v1 + v3 + v4) * 0.3333333333333333

    # evaluate influence
    potential2, velocity2, velocity_gradient2 = _induced(target, (v1, v3, v4), control_point, strength, TK, core_size, R, derivatives_switch, fam)

    # return potential, velocity, velocity_gradient
    return potential+potential2, velocity+velocity2, velocity_gradient+velocity_gradient2
end

# function barrier: family in the type domain inside the loop (BRAINSTORM 025)
function FastMultipole.direct!(target_system, target_index, derivatives_switch::FastMultipole.DerivativesSwitch, source_system::PanelWake, source_buffer, source_index)
    _direct_panelwake!(target_system, target_index, derivatives_switch, source_system,
        source_buffer, source_index, Val(FILAMENT_REGULARIZATION[]))
end

function _direct_panelwake!(target_system, target_index, derivatives_switch::FastMultipole.DerivativesSwitch{PS,GS,HS,NO,NM}, source_system::PanelWake, source_buffer, source_index, fam::Val) where {PS,GS,HS,NO,NM}
    TF = eltype(target_system)
    for i_target in target_index # loop over targets
        target = FastMultipole.StaticArrays.SVector{3,TF}(target_system[1, i_target],
                  target_system[2, i_target],
                  target_system[3, i_target])

        phi_out = zero(eltype(target_system))
        U_out = @SVector zeros(eltype(target_system), 3)
        H_out = zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})

        for i_source in source_index # loop over sources
            # evaluate influence due to this source
            phi, U, H = induced(target, source_system, source_buffer, i_source, derivatives_switch, fam)
            phi_out += phi
            U_out += U
            if HS
                H_out += H
            end
        end

        # store results
        if PS
            FastMultipole.set_scalar_potential!(target_system, derivatives_switch, i_target, phi_out)
        end
        if GS
            FastMultipole.set_gradient!(target_system, derivatives_switch, i_target, U_out)
        end
        if HS
            FastMultipole.set_hessian!(target_system, derivatives_switch, i_target, H_out)
        end
        if HS
            # @warn "Hessian output not currently implemented for PanelWake targets"
        end
    end
end

FastMultipole.body_to_multipole!(system::PanelWake{ConstantDoublet, 1, <:Any}, args...) =
    FastMultipole.body_to_multipole_quad!(FastMultipole.Panel{FastMultipole.Dipole}, system, args...)

FastMultipole.body_to_multipole!(system::PanelWake{VortexRing, 1, <:Any}, args...) =
    FastMultipole.body_to_multipole_quad!(FastMultipole.Panel{FastMultipole.Dipole}, system, args...)

function propagate!(wake::PanelWake, dt; step=0, frames=nothing)
    for i_surf in eachindex(wake.nodes)
        if wake.freestream_convection
            # every row convects with the freestream only: the sheet stays
            # straight along U∞ (semi-infinite-wake geometry). The `reshape` is
            # required because this is a 3×R×C view, not the 3×C view used by
            # the row-1-only branch below.
            view(wake.nodes[i_surf], :, 1:wake.nwakes[]+1, :) .+=
                dt .* reshape(wake.freestream, 3, 1, 1)
        elseif wake.shed_with_induced_velocity
            view(wake.velocity[i_surf], :, 1:wake.nwakes[]+1, :) .*= dt # displacements
            view(wake.nodes[i_surf], :, 1:wake.nwakes[]+1, :) .+= view(wake.velocity[i_surf], :, 1:wake.nwakes[]+1, :) # update nodes
            view(wake.velocity[i_surf], :, 1:wake.nwakes[]+1, :) ./= dt # restore velocities
        else
            view(wake.nodes[i_surf], :, 1, :) .+= dt .* wake.freestream
            if wake.nwakes[] >= 1
                view(wake.velocity[i_surf], :, 2:wake.nwakes[]+1, :) .*= dt # displacements
                view(wake.nodes[i_surf], :, 2:wake.nwakes[]+1, :) .+= view(wake.velocity[i_surf], :, 2:wake.nwakes[]+1, :) # update nodes
                view(wake.velocity[i_surf], :, 2:wake.nwakes[]+1, :) ./= dt # restore velocities
            end
        end
    end
end

# Forward declaration (methods defined later in the vortex filament section)
struct FilamentWrapper{TS}
    system::TS
end

function write_vtk(name, wake::PanelWake, idx, t; overwrite=false, compress::Bool=true,
        filament_name=nothing)
    # Route block files to a subdirectory named after the PVD
    _parent, _base = splitdir(name)
    subdir = joinpath(_parent, _base)
    mkpath(subdir)
    block_name = joinpath(subdir, _base)

    vtm = WriteVTK.vtk_multiblock(block_name * ".$idx.vtm")
    if wake.nwakes[] > 0
        for i_surf in eachindex(wake.nodes)
            pts = view(wake.nodes[i_surf], :, 1:wake.nwakes[]+1, :)
            pts_reshaped = reshape(pts, 3, wake.nwakes[]+1, size(wake.nodes[i_surf], 3), 1)
            WriteVTK.vtk_grid(vtm, block_name * ".$(i_surf).$(idx).vts", pts_reshaped; compress) do vtk
                vel = view(wake.velocity[i_surf], :, 1:wake.nwakes[]+1, :)
                vtk["velocity", WriteVTK.VTKPointData()] = reshape(vel, 3, wake.nwakes[]+1, size(wake.nodes[i_surf], 3), 1)

                str = view(wake.strength[i_surf], :, 1:wake.nwakes[], :)
                vtk["strength", WriteVTK.VTKCellData()] = reshape(str, size(str, 1), wake.nwakes[], size(wake.nodes[i_surf], 3)-1, 1)
            end
        end
    elseif wake.convert_at_shed
        # Convert-at-shed (BRAINSTORM 024): the solve-visible sheet is always
        # empty, but the row-1 node line (TE+Das) and its velocities are live
        # state — warmstart/replay need them to replay the end-of-step shed —
        # so write a single-node-row grid with no cells (hence no strengths;
        # the row-1 strength is carried by the metadata terminal_strength).
        for i_surf in eachindex(wake.nodes)
            n_cols = size(wake.nodes[i_surf], 3)
            pts = reshape(view(wake.nodes[i_surf], :, 1:1, :), 3, 1, n_cols, 1)
            WriteVTK.vtk_grid(vtm, block_name * ".$(i_surf).$(idx).vts", pts; compress) do vtk
                vel = view(wake.velocity[i_surf], :, 1:1, :)
                vtk["velocity", WriteVTK.VTKPointData()] = reshape(vel, 3, 1, n_cols, 1)
            end
        end
    end
    WriteVTK.vtk_save(vtm)
    _pvd_append!(name * ".pvd", t, joinpath(_base, _base * ".$idx.vtm"); overwrite)

    # filaments at trailing edge of last panel row
    write_vtk(isnothing(filament_name) ? name * "_filaments" : filament_name,
        FilamentWrapper(wake), idx, t; overwrite, compress)
end

function write_vtk(name, filaments::FilamentWrapper{<:PanelWake}, idx, t; overwrite=false, compress::Bool=true)
    wake = filaments.system
    if !wake.overflowed[]
        # no filaments yet — still write an empty PVD entry
        _parent, _base = splitdir(name)
        subdir = joinpath(_parent, _base)
        mkpath(subdir)
        block_name = joinpath(subdir, _base)
        vtm = WriteVTK.vtk_multiblock(block_name * ".$idx.vtm")
        WriteVTK.vtk_save(vtm)
        _pvd_append!(name * ".pvd", t, joinpath(_base, _base * ".$idx.vtm"); overwrite)
        return
    end

    i_row = wake.nwakes[]

    _parent, _base = splitdir(name)
    subdir = joinpath(_parent, _base)
    mkpath(subdir)
    block_name = joinpath(subdir, _base)

    vtm = WriteVTK.vtk_multiblock(block_name * ".$idx.vtm")

    for i_surf in eachindex(wake.nodes)
        n_fils = size(wake.strength[i_surf], 3)

        # build points array (3, 2*n_fils) and VTK_LINE cells
        points = zeros(eltype(wake.nodes[i_surf]), 3, 2 * n_fils)
        cells = Vector{WriteVTK.MeshCell{WriteVTK.VTKCellTypes.VTKCellType, Vector{Int}}}(undef, n_fils)
        strengths = zeros(eltype(wake.strength[i_surf]), n_fils)

        for j in 1:n_fils
            ip = 2 * (j - 1)
            points[:, ip + 1] .= view(wake.nodes[i_surf], :, i_row + 1, j)
            points[:, ip + 2] .= view(wake.nodes[i_surf], :, i_row + 1, j + 1)
            cells[j] = WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_LINE, [ip + 1, ip + 2])
            strengths[j] = _final_filament_strength(wake, i_surf, i_row, j)
        end

        WriteVTK.vtk_grid(vtm, block_name * ".$(i_surf).$(idx).vtu", points, cells; compress) do vtk
            vtk["strength", WriteVTK.VTKCellData()] = strengths
        end
    end

    WriteVTK.vtk_save(vtm)
    _pvd_append!(name * ".pvd", t, joinpath(_base, _base * ".$idx.vtm"); overwrite)
end

#--- Wake Shedding Methods ---#

"""
    WakeSheddingMethod

Abstract supertype for particle-shedding strategies used by
[`PanelParticleWake`](@ref).
"""
abstract type WakeSheddingMethod end

"""No particle shedding."""
struct NoShed <: WakeSheddingMethod end

"""Gaussian particle shedding with prescribed smoothing width `sigma`."""
struct SigmaPPS{TF} <: WakeSheddingMethod
    sigma::TF
    p_per_step::Int
end

"""Gaussian particle shedding with fixed `sigma` and minimum overlap."""
struct SigmaOverlap{TF} <: WakeSheddingMethod
    sigma::TF
    overlap::TF
end

"""Overlap-based particle shedding with target overlap ratio."""
struct OverlapPPS{TF} <: WakeSheddingMethod
    overlap::TF
    p_per_step::Int
end

function _shed_particles!(pfield, r1, r2, Γ, method::OverlapPPS)
    dist = LA.norm(r2 - r1)
    dist < eps(typeof(dist)) && return nothing
    sigma = dist * method.overlap / method.p_per_step
    return _shed_particles!(pfield, r1, r2, Γ, SigmaPPS(sigma, method.p_per_step))
end

function _shed_particles!(pfield, r1, r2, Γ, method::SigmaOverlap)
    dist = LA.norm(r2 - r1)
    dist < eps(typeof(dist)) && return nothing
    p_per_step = max(1, ceil(Int, method.overlap * dist / method.sigma))
    return _shed_particles!(pfield, r1, r2, Γ, SigmaPPS(method.sigma, p_per_step))
end

function _shed_particles!(pfield, r1, r2, Γ, method::SigmaPPS)
    # A zero-strength particle carries no vorticity but its zero |Γ| divides
    # relaxation (corrected Pedrizzetti) into NaN, which then contaminates the
    # whole field. Exact zeros occur generically on symmetric wings (mid-span
    # trailing difference) and at impulsive starts (zero-strength first row).
    Γ == 0 && return nothing
    sigma = method.sigma
    p_per_step = method.p_per_step
    distance_vector = (r2 - r1) / p_per_step
    Xp = r1 + distance_vector * 0.5
    Γp = Γ * distance_vector
    for _ in 1:p_per_step
        FLOWVPM.add_particle(pfield, Xp, Γp, sigma; circulation=Γ)
        Xp += distance_vector
    end
end

function _shed_particles!(pfield, r1, r2, Γ, ::NoShed)
    return nothing
end

"""
    StationSigmaOverlap(sigmas, overlap)

Station-indexed Gaussian particle shedding (BRAINSTORM 018 Phase 16 chord–σ
co-scaling): `sigmas[i_surf][j]` is the smoothing width at wake node column
`j` of shedding surface `i_surf` (one vector per surface, each of length
`n_cols + 1`). Trailing filaments at node column `j` shed at
`sigmas[i_surf][j]`; the unsteady filament spanning columns `j → j+1` sheds
at their mean. Each filament resolves to a [`SigmaOverlap`](@ref), so
particle count and placement follow the same overlap rule.

Only [`LegacyEdgeJumpConversion`](@ref) resolves stations; calling
`_shed_particles!` with this method directly is an error.
"""
struct StationSigmaOverlap{TF} <: WakeSheddingMethod
    sigmas::Vector{Vector{TF}}
    overlap::TF

    function StationSigmaOverlap(sigmas::Vector{Vector{TF}}, overlap::Real) where {TF}
        isempty(sigmas) && throw(ArgumentError(
            "StationSigmaOverlap needs at least one per-surface sigma vector"))
        for (k, sig) in enumerate(sigmas)
            all(s -> isfinite(s) && s > 0, sig) || throw(ArgumentError(
                "StationSigmaOverlap: sigmas for surface $(k) must be finite " *
                "and positive"))
        end
        overlap > 0 || throw(ArgumentError("StationSigmaOverlap: overlap must be positive"))
        return new{TF}(sigmas, TF(overlap))
    end
end

function _shed_particles!(pfield, r1, r2, Γ, ::StationSigmaOverlap)
    throw(ArgumentError("StationSigmaOverlap must be resolved to a station " *
                        "(via _station_method) before shedding"))
end

# Station resolution for the legacy conversion loop. Identity for every
# span-uniform method (zero behavior change); the station-indexed method
# resolves to a plain SigmaOverlap at the mean of the two node columns the
# filament touches (j1 == j2 for trailing filaments).
@inline _station_method(m::WakeSheddingMethod, i_surf, j1, j2) = m
@inline function _station_method(m::StationSigmaOverlap, i_surf, j1, j2)
    sig = m.sigmas[i_surf]
    return SigmaOverlap(0.5 * (sig[j1] + sig[j2]), m.overlap)
end

# One-shot per-surface validation so a mis-sized sigma vector fails with the
# expected length instead of an opaque BoundsError mid-conversion. No-op for
# span-uniform methods.
_validate_station_method(::WakeSheddingMethod, i_surf, n_node_cols) = nothing
function _validate_station_method(m::StationSigmaOverlap, i_surf, n_node_cols)
    i_surf <= length(m.sigmas) || error(
        "StationSigmaOverlap: no sigma vector for shedding surface $(i_surf) " *
        "(have $(length(m.sigmas)))")
    length(m.sigmas[i_surf]) == n_node_cols || error(
        "StationSigmaOverlap: sigma vector for surface $(i_surf) has length " *
        "$(length(m.sigmas[i_surf])), expected n_cols + 1 = $(n_node_cols)")
    return nothing
end

"""
Sentinel used as the *implementation* default of `PanelParticleWake`'s
`method_trailing` and `method_unsteady` keywords, so the constructor can tell
"caller said nothing" from "caller explicitly asked for the legacy default".

The legacy conversion resolves it to a fresh `OverlapPPS(1.3, 2)`, which is the
unchanged public default. The surface-vorticity conversion has no line policies
to resolve and stores the sentinel itself, so an accidental use is an error
rather than a silently-different shedding rule.
"""
struct DefaultWakeSheddingMethod <: WakeSheddingMethod end

function _shed_particles!(pfield, r1, r2, Γ, ::DefaultWakeSheddingMethod)
    throw(ArgumentError("DefaultWakeSheddingMethod is an unresolved constructor " *
                        "sentinel and must never be used to shed particles"))
end

#--- Panel-to-particle conversion strategies (BRAINSTORM 016) ---#

"""
    AbstractPanelParticleConversion

Abstract supertype for strategies that convert the outgoing row of a
[`PanelWake`](@ref) into vortex particles inside a
[`PanelParticleWake`](@ref).

Select one through the single `conversion` keyword of `PanelParticleWake`.
There is no automatic fallback between strategies.
"""
abstract type AbstractPanelParticleConversion end

"""
    LegacyEdgeJumpConversion()

Convert the outgoing panel row from panel-edge strength *jumps*: a streamwise
trailing filament carrying the spanwise jump, a spanwise unsteady filament
carrying the streamwise (time) jump, and a closing filament at the tip of a
non-wrapping chain. Line policies come from `method_trailing` and
`method_unsteady`.

This is the exact default and is unchanged by BRAINSTORM 016.
"""
struct LegacyEdgeJumpConversion <: AbstractPanelParticleConversion end

"""
    SurfaceVorticityConversion(sigma; overlap=1.3, rank_rtol, geometry_rtol,
                               diagnose_nearfield=false)

Opt-in conversion that reconstructs the *smooth* surface vorticity of the
outgoing sheet and deposits area-weighted particles

```math
\\boldsymbol\\Gamma_p = \\boldsymbol\\kappa(\\boldsymbol x_p)\\,\\Delta A_p,
\\qquad
\\boldsymbol\\kappa = \\boldsymbol n \\times \\nabla_s \\mu
                    = -\\boldsymbol n \\times \\nabla_s \\hat\\mu ,
```

where ``\\hat\\mu = -\\mu`` is the stored wake strength. Internal panel-to-panel
edges are cancelled rather than reproduced as filaments; only the true open
root and tip streamwise closures are retained, sampled with
`SigmaOverlap(sigma, overlap)`. All surface and line sampling targets spacing
``h = \\sigma/\\mathrm{overlap}``, so every particle this wake deposits shares
the fixed smoothing width `sigma`.

`attribution` selects which streamwise face of the outgoing row this conversion
deposits: `:upstream` (the default retained by the Phase 3 static-equivalence
acceptance evidence -- the whole
upstream face, leaving nothing at the panel/particle partition), `:downstream`
(the legacy attribution), or `:split` (half of each, which makes the streamwise
difference centered and is second-order accurate on a graded wake, at the cost
of a retained half-jump filament at the handoff). All three are exactly
circulation conserving.

`rank_rtol` sets the relative singular-value threshold of the two-point
gradient stencil (rank deficiency is diagnostic, not fatal); `geometry_rtol`
rejects a vanishing metric scale. `diagnose_nearfield` additionally records
panel-induced velocity-gradient distributions at every candidate particle,
separately for interior, root/tip, and perimeter classes.

This strategy does not use `method_trailing` or `method_unsteady`; supplying
either alongside it is a configuration error.

The smooth handoff also requires `unsteady_filament=true` and
`include_final_filament=true`. Its startup ledger assumes the unsteady final
filament carries the physical aft face before the first conversion, and its
steady ledger uses that filament to cancel the panel/particle interface.

See BRAINSTORM item 016.
"""
struct SurfaceVorticityConversion{TF} <: AbstractPanelParticleConversion
    sigma::TF
    overlap::TF
    rank_rtol::TF
    geometry_rtol::TF
    attribution::Symbol
    diagnose_nearfield::Bool
end

"Streamwise attribution modes and their upstream-face fraction `alpha`."
const _ATTRIBUTION_ALPHA = (upstream = 1.0, downstream = 0.0, split = 0.5)

function SurfaceVorticityConversion(sigma;
        overlap = 1.3,
        rank_rtol = sqrt(eps(float(typeof(sigma)))),
        geometry_rtol = sqrt(eps(float(typeof(sigma)))),
        attribution::Symbol = :upstream,
        diagnose_nearfield::Bool = false,
    )
    haskey(_ATTRIBUTION_ALPHA, attribution) || throw(ArgumentError(
        "SurfaceVorticityConversion attribution must be one of " *
        "$(keys(_ATTRIBUTION_ALPHA)) (got :$(attribution))"))
    sigma_f, overlap_f, rank_f, geom_f =
        promote(float(sigma), float(overlap), float(rank_rtol), float(geometry_rtol))

    for (name, val) in (("sigma", sigma_f), ("overlap", overlap_f),
                        ("rank_rtol", rank_f), ("geometry_rtol", geom_f))
        isfinite(val) ||
            throw(ArgumentError("SurfaceVorticityConversion $name must be finite (got $val)"))
    end
    sigma_f > 0 ||
        throw(ArgumentError("SurfaceVorticityConversion sigma must be positive (got $sigma_f)"))
    overlap_f > 0 ||
        throw(ArgumentError("SurfaceVorticityConversion overlap must be positive (got $overlap_f)"))
    rank_f >= 0 ||
        throw(ArgumentError("SurfaceVorticityConversion rank_rtol must be nonnegative (got $rank_f)"))
    geom_f > 0 ||
        throw(ArgumentError("SurfaceVorticityConversion geometry_rtol must be positive (got $geom_f)"))

    return SurfaceVorticityConversion{typeof(sigma_f)}(
        sigma_f, overlap_f, rank_f, geom_f, attribution, diagnose_nearfield)
end

"""
Resolve one of the legacy line-policy keywords against the selected conversion
strategy. `field` names the keyword for error messages.
"""
_resolve_line_policy(::LegacyEdgeJumpConversion, ::DefaultWakeSheddingMethod, field) =
    OverlapPPS(1.3, 2)

_resolve_line_policy(::LegacyEdgeJumpConversion, method::WakeSheddingMethod, field) = method

_resolve_line_policy(::SurfaceVorticityConversion, sentinel::DefaultWakeSheddingMethod, field) =
    sentinel

function _resolve_line_policy(::SurfaceVorticityConversion, method::WakeSheddingMethod, field)
    throw(ArgumentError(
        "SurfaceVorticityConversion does not use $field; it samples its root/tip " *
        "closure with SigmaOverlap(sigma, overlap). Remove the $field keyword " *
        "rather than leaving it silently ignored (got $(method))."))
end

"Validate the final-filament source contract required by a conversion strategy."
_validate_conversion_filaments(::LegacyEdgeJumpConversion, unsteady, included) = nothing

function _validate_conversion_filaments(::SurfaceVorticityConversion,
                                        unsteady::Bool, included::Bool)
    unsteady || throw(ArgumentError(
        "SurfaceVorticityConversion requires unsteady_filament=true: " *
        "without the unsteady final filament the startup aft face is absent."))
    included || throw(ArgumentError(
        "SurfaceVorticityConversion requires include_final_filament=true: " *
        "the final filament is the source that cancels the panel/particle interface."))
    return nothing
end

#--- Surface-vorticity reconstruction core (BRAINSTORM 016 Phase 2 sec. 5) ---#

"""
    WakeGeometryError

A wake panel or gradient stencil failed geometric validation: a repeated or
metric-singular vertex, a vanishing or inconsistently oriented Jacobian, a
folded or inverted bilinear panel, or a non-finite coordinate.
"""
struct WakeGeometryError <: Exception
    msg::String
end

Base.showerror(io::IO, e::WakeGeometryError) = print(io, "WakeGeometryError: ", e.msg)

"""
    _deterministic_tangent_basis(n)

Orthonormal tangent basis `(t1, t2)` of the plane with unit normal `n`, chosen
deterministically: project the Cartesian axis *least* aligned with `n` (ties
broken in x, y, z order) and complete with `t2 = n x t1`.

Fixing the basis this way keeps rigid-rotation fixtures reproducible. It is a
reporting convenience only -- the reconstructed physical gradient is
independent of which tangent basis is used.
"""
function _deterministic_tangent_basis(n::SVector{3,TF}) where {TF}
    ax = 1
    best = abs(n[1])
    for i in 2:3
        if abs(n[i]) < best
            best = abs(n[i])
            ax = i
        end
    end
    e = SVector{3,TF}(ax == 1, ax == 2, ax == 3)
    t1 = e - LA.dot(e, n) * n
    nt1 = LA.norm(t1)
    nt1 > 0 || throw(WakeGeometryError("degenerate normal $(n); no tangent basis exists"))
    t1 = t1 / nt1
    return t1, LA.cross(n, t1)
end

"""
Outcome of the two-point surface-gradient reconstruction on one wake panel.

`gradient` is the physical tangential gradient of the *stored* strength
``\\hat\\mu``. `singular_values`, `rank`, `condition`, and
`observable_directions` describe the scaled geometry matrix and are diagnostic:
rank deficiency is supported, not an error.
"""
struct SurfaceGradientResult{TF}
    gradient::SVector{3,TF}
    rank::Int
    singular_values::SVector{2,TF}
    rank_threshold::TF
    condition::TF
    observable_directions::SVector{2,SVector{3,TF}}
    observable::SVector{2,Bool}
    geometry_scale::TF
end

"""
    _reconstruct_surface_gradient(n, d1, dmu1, d2, dmu2, rank_rtol, geometry_rtol,
                                  reference_length)

Reconstruct the tangential gradient of the stored strength on a panel whose
unit normal is `n`, from two centroid-difference equations

```math
\\boldsymbol d_k \\cdot \\nabla_s\\hat\\mu = \\Delta\\hat\\mu_k ,
\\qquad k = 1, 2 ,
```

where `d1` is the streamwise displacement to the upstream neighbour and `d2`
the spanwise displacement (centered where both neighbours exist, one-sided at a
true root or tip).

The `2 x 2` system is formed in the deterministic tangent basis and scaled by
the largest stencil length before its SVD, so `rank_rtol` is a genuine relative
threshold. Rank 2 uses the full solve, rank 1 the minimum-norm pseudoinverse,
and rank 0 returns a zero gradient; in every case the returned `gradient` is
the same physical vector regardless of basis or SVD sign conventions.
"""
function _reconstruct_surface_gradient(n::SVector{3,TF},
        d1::SVector{3,TF}, dmu1::TF,
        d2::SVector{3,TF}, dmu2::TF,
        rank_rtol::TF, geometry_rtol::TF, reference_length::TF) where {TF}

    all(isfinite, n) && all(isfinite, d1) && all(isfinite, d2) &&
        isfinite(dmu1) && isfinite(dmu2) ||
        throw(WakeGeometryError("non-finite gradient stencil"))

    scale = max(LA.norm(d1), LA.norm(d2))
    scale > geometry_rtol * reference_length ||
        throw(WakeGeometryError(
            "gradient stencil has vanishing metric scale $(scale) relative to " *
            "reference length $(reference_length)"))

    t1, t2 = _deterministic_tangent_basis(n)

    # Column-major fill: A = [d1.t1  d1.t2 ; d2.t1  d2.t2], scaled by `scale`
    # so the singular values are dimensionless and rank_rtol is relative.
    A = SMatrix{2,2,TF}(LA.dot(d1, t1) / scale, LA.dot(d2, t1) / scale,
                        LA.dot(d1, t2) / scale, LA.dot(d2, t2) / scale)
    b = SVector{2,TF}(dmu1 / scale, dmu2 / scale)

    F = LA.svd(A)
    s = SVector{2,TF}(F.S[1], F.S[2])
    threshold = rank_rtol * s[1]
    rnk = count(>(threshold), s)

    g = SVector{2,TF}(zero(TF), zero(TF))
    for i in 1:rnk
        # Minimum-norm solution: only the observable singular directions
        # contribute, so rank 1 falls out of the same expression as rank 2.
        g += (LA.dot(view(F.U, :, i), b) / s[i]) * SVector{2,TF}(F.V[1, i], F.V[2, i])
    end

    dirs = SVector{2,SVector{3,TF}}(
        F.V[1, 1] * t1 + F.V[2, 1] * t2,
        F.V[1, 2] * t1 + F.V[2, 2] * t2)
    obs = SVector{2,Bool}(1 <= rnk, 2 <= rnk)
    cond = rnk == 2 ? s[1] / s[2] : TF(Inf)

    return SurfaceGradientResult{TF}(g[1] * t1 + g[2] * t2, rnk, s, threshold,
                                     cond, dirs, obs, scale)
end

"""
Placeholder [`SurfaceGradientResult`](@ref) for a panel with no observable
stencil leg at all (a single-row, single-column wake). Reported as rank 0 rather
than fabricating a direction.
"""
_null_surface_gradient(::Type{TF}) where {TF} = SurfaceGradientResult{TF}(
    zero(SVector{3,TF}), 0, zero(SVector{2,TF}), zero(TF), TF(Inf),
    SVector{2,SVector{3,TF}}(zero(SVector{3,TF}), zero(SVector{3,TF})),
    SVector{2,Bool}(false, false), zero(TF))

"""
    _surface_vorticity(n, grad_muhat)

Surface vorticity from the *stored* strength gradient, using the package sign
convention ``\\hat\\mu = -\\mu``:

```math
\\boldsymbol\\kappa = \\boldsymbol n \\times \\nabla_s\\mu
                    = -\\boldsymbol n \\times \\nabla_s\\hat\\mu .
```
"""
_surface_vorticity(n::SVector{3,TF}, grad_muhat::SVector{3,TF}) where {TF} =
    -LA.cross(n, grad_muhat)

"""
    _bilinear_position(v1, v2, v3, v4, xi, eta)

Point on the bilinear wake panel in the package's contiguous vertex order
`x(0,0)=v1`, `x(1,0)=v2`, `x(1,1)=v3`, `x(0,1)=v4`.
"""
_bilinear_position(v1, v2, v3, v4, xi, eta) =
    (1 - xi) * (1 - eta) * v1 + xi * (1 - eta) * v2 + xi * eta * v3 + (1 - xi) * eta * v4

"""
    _bilinear_derivatives(v1, v2, v3, v4, xi, eta)

Covariant derivatives `(dx/dxi, dx/deta)` of the bilinear panel.
"""
function _bilinear_derivatives(v1, v2, v3, v4, xi, eta)
    dxi = (1 - eta) * (v2 - v1) + eta * (v3 - v4)
    deta = (1 - xi) * (v4 - v1) + xi * (v3 - v2)
    return dxi, deta
end

"""
    _bilinear_normal(v1, v2, v3, v4, xi, eta)

Oriented unit normal and Jacobian magnitude `|dx/dxi x dx/deta|` at a
parametric point. Throws [`WakeGeometryError`](@ref) where the panel is
metric-singular.
"""
function _bilinear_normal(v1, v2, v3, v4, xi, eta, geometry_rtol, reference_length)
    dxi, deta = _bilinear_derivatives(v1, v2, v3, v4, xi, eta)
    c = LA.cross(dxi, deta)
    jac = LA.norm(c)
    isfinite(jac) ||
        throw(WakeGeometryError("non-finite Jacobian at (xi, eta) = ($xi, $eta)"))
    # A bilinear quad may be warped, but its metric must never vanish: that is
    # a repeated edge, a folded panel, or a self-intersection.
    jac > geometry_rtol * reference_length^2 ||
        throw(WakeGeometryError(
            "vanishing Jacobian $(jac) at (xi, eta) = ($xi, $eta) relative to " *
            "reference area $(reference_length^2)"))
    return c / jac, jac
end

"""
    _subdivision_counts(v1, v2, v3, v4, h)

Number of quadrature subcells in each parametric direction so that subcell
edges do not exceed the target spacing `h`. At least one subdivision is used in
both directions, per the Phase 2 resolution floor.
"""
function _subdivision_counts(v1, v2, v3, v4, h)
    l_xi = max(LA.norm(v2 - v1), LA.norm(v3 - v4))
    l_eta = max(LA.norm(v4 - v1), LA.norm(v3 - v2))
    return (max(1, ceil(Int, l_xi / h)), max(1, ceil(Int, l_eta / h)))
end

# Two-point Gauss-Legendre nodes/weights on [-1, 1].
const _GAUSS2_NODE = 1 / sqrt(3)

"""
    _subcell_area(v1, v2, v3, v4, xi0, xi1, eta0, eta1, geometry_rtol, reference_length)

Physical area of one parametric subcell of the bilinear panel, by a `2 x 2`
Gauss rule applied to the norm of the bilinear Jacobian. Exact for a planar
parallelogram and second-order accurate on a warped quad.
"""
function _subcell_area(v1, v2, v3, v4, xi0, xi1, eta0, eta1,
                       geometry_rtol, reference_length)
    dxi = xi1 - xi0
    deta = eta1 - eta0
    half_xi = dxi / 2
    half_eta = deta / 2
    mid_xi = (xi0 + xi1) / 2
    mid_eta = (eta0 + eta1) / 2

    area = zero(eltype(v1))
    for sx in (-_GAUSS2_NODE, _GAUSS2_NODE), sy in (-_GAUSS2_NODE, _GAUSS2_NODE)
        xi = mid_xi + half_xi * sx
        eta = mid_eta + half_eta * sy
        _, jac = _bilinear_normal(v1, v2, v3, v4, xi, eta, geometry_rtol, reference_length)
        area += jac
    end
    return area * half_xi * half_eta
end

"""
    _validate_wake_panel(v1, v2, v3, v4, geometry_rtol)

Validate one bilinear wake panel before it is sampled: finite vertices, a
nondegenerate scale, a metric that never vanishes, and an orientation that is
consistent across the panel (no fold or inversion). Returns the panel's
reference length.
"""
function _validate_wake_panel(v1, v2, v3, v4, geometry_rtol)
    for v in (v1, v2, v3, v4)
        all(isfinite, v) || throw(WakeGeometryError("non-finite wake panel vertex $(v)"))
    end

    reference_length = max(LA.norm(v2 - v1), LA.norm(v3 - v4),
                           LA.norm(v4 - v1), LA.norm(v3 - v2))
    reference_length > 0 ||
        throw(WakeGeometryError("wake panel has zero extent (all vertices coincide)"))

    # Sample the four corners plus the center. A bilinear Jacobian is bilinear
    # in (xi, eta), so a sign flip anywhere on the panel shows up at a corner;
    # the center guards the warped-but-corner-consistent case.
    reference_normal = nothing
    for (xi, eta) in ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.5, 0.5))
        n, _ = _bilinear_normal(v1, v2, v3, v4, xi, eta, geometry_rtol, reference_length)
        if reference_normal === nothing
            reference_normal = n
        elseif LA.dot(n, reference_normal) <= 0
            throw(WakeGeometryError(
                "wake panel is folded or inverted: normal reverses between the " *
                "reference corner and (xi, eta) = ($xi, $eta)"))
        end
    end
    return reference_length
end

#--- Surface-vorticity transaction state (BRAINSTORM 016 Phase 2 sec. 6-8) ---#

"""
    SurfaceParticleClass

Classification of a candidate particle produced by
[`SurfaceVorticityConversion`](@ref).

`InteriorSurfaceParticle` and `RootTipSurfaceParticle` are both *area*
particles ``\\boldsymbol\\Gamma_p = \\boldsymbol\\kappa\\,\\Delta A_p``; they are
distinguished so the root/tip near-field caveat of Phase 1 sec. 5.7 can be
diagnosed without pooling those samples with interior ones.
`PerimeterLineParticle` is the true open-chain root/tip streamwise closure.
"""
@enum SurfaceParticleClass begin
    InteriorSurfaceParticle
    RootTipSurfaceParticle
    PerimeterLineParticle
end

"""
    PanelParticleCapacityError(requested, available)

The preflight of a surface-vorticity conversion needed `requested` particle
slots but the `ParticleField` had only `available`. Thrown *before* any
mutation, so both the particle field and every panel row are bitwise unchanged.
"""
struct PanelParticleCapacityError <: Exception
    requested::Int
    available::Int
end

Base.showerror(io::IO, e::PanelParticleCapacityError) = print(io,
    "PanelParticleCapacityError: surface-vorticity conversion requested ",
    e.requested, " particles but only ", e.available, " slots remain")

"""
    WakeConversionStateError(msg)

The wake was not in a state in which a conversion transaction is defined (for
example the row buffer is not full, or the continuation information a strategy
needs is unavailable).
"""
struct WakeConversionStateError <: Exception
    msg::String
end

Base.showerror(io::IO, e::WakeConversionStateError) =
    print(io, "WakeConversionStateError: ", e.msg)

"""
    WakeContinuationStateError(msg)

A smooth-conversion replay or warm start could not restore an unambiguous
panel/particle handoff state from its per-step metadata.
"""
struct WakeContinuationStateError <: Exception
    msg::String
end

Base.showerror(io::IO, e::WakeContinuationStateError) =
    print(io, "WakeContinuationStateError: ", e.msg)

"""
Reusable staging buffers for one surface-vorticity conversion.

Every candidate particle is materialized here during preflight; nothing touches
the live `ParticleField` until the whole transaction is validated. The buffers
are `empty!`-ed (not reallocated) at the start of each conversion, so after the
first transaction of a run the steady-state path performs no fresh allocation
for candidate staging.
"""
struct SurfaceVorticityWorkspace{TF}
    positions::Vector{SVector{3,TF}}
    strengths::Vector{SVector{3,TF}}
    sigmas::Vector{TF}
    circulations::Vector{TF}
    classes::Vector{SurfaceParticleClass}
    # Per-panel subcell areas. The deposition needs the *summed* panel area
    # before it can set kappa = V / A, so the subcell areas are computed in a
    # first pass and reused in the second; taking the sum (rather than an
    # independently computed area) is what makes sum(kappa * dA) == V exactly.
    areas::Vector{TF}
end

SurfaceVorticityWorkspace{TF}() where {TF} = SurfaceVorticityWorkspace{TF}(
    SVector{3,TF}[], SVector{3,TF}[], TF[], TF[], SurfaceParticleClass[], TF[])

function _reset_workspace!(ws::SurfaceVorticityWorkspace)
    empty!(ws.positions)
    empty!(ws.strengths)
    empty!(ws.sigmas)
    empty!(ws.circulations)
    empty!(ws.classes)
    empty!(ws.areas)
    return ws
end

_n_pending(ws::SurfaceVorticityWorkspace) = length(ws.positions)

function _stage_particle!(ws::SurfaceVorticityWorkspace{TF}, X, Gamma, sigma,
                          circulation, class::SurfaceParticleClass) where {TF}
    push!(ws.positions, X)
    push!(ws.strengths, Gamma)
    push!(ws.sigmas, TF(sigma))
    push!(ws.circulations, TF(circulation))
    push!(ws.classes, class)
    return nothing
end

"""
Per-ghost-panel diagnostic record (Phase 2 sec. 8.1). Immutable: the cumulative
ledger accumulates alongside it and never rewrites a committed record.
"""
struct SurfaceVorticityPanelRecord{TF}
    i_surf::Int
    j::Int
    rank::Int
    singular_values::SVector{2,TF}
    rank_threshold::TF
    condition::TF
    observable_directions::SVector{2,SVector{3,TF}}
    observable::SVector{2,Bool}
    n_xi::Int
    n_eta::Int
    area::TF
    deposited_strength::SVector{3,TF}
    # Deposited surface vorticity (divergence form) and, for comparison only,
    # what the Stage 2 centroid-stencil reconstruction would have given. They
    # agree identically on a uniform mesh and differ by O(row-to-row stretch)
    # otherwise, so `kappa_difference` is a direct grid-nonuniformity measure.
    kappa_conservative::SVector{3,TF}
    kappa_reconstruction::SVector{3,TF}
    kappa_difference::TF
    handoff::SVector{3,TF}
    downstream_face::SVector{3,TF}
    class::SurfaceParticleClass
    n_requested::Int
    n_elided::Int
    geometry_scale::TF
    min_jacobian::TF
end

"""
Immutable summary of panel-wake-induced velocity-gradient magnitudes sampled
at one class of staged conversion particles. `minimum` through `p95` are `NaN`
when `count == 0`.
"""
struct SurfaceVorticityNearFieldSummary{TF}
    count::Int
    minimum::TF
    maximum::TF
    mean::TF
    rms::TF
    median::TF
    p95::TF
end

"""
Optional near-field diagnostic summaries for one conversion. The three fields
separate interior area, root/tip area, and true-perimeter line particles.
"""
struct SurfaceVorticityNearFieldDiagnostics{TF}
    interior::SurfaceVorticityNearFieldSummary{TF}
    roottip::SurfaceVorticityNearFieldSummary{TF}
    perimeter::SurfaceVorticityNearFieldSummary{TF}
end

"""
Per-conversion diagnostic record (Phase 2 sec. 8.2), written only on a
successful commit.
"""
struct SurfaceVorticityConversionRecord{TF}
    ordinal::Int
    step_id::Int
    panels::Vector{SurfaceVorticityPanelRecord{TF}}
    n_interior::Int
    n_roottip::Int
    n_perimeter::Int
    total_area::TF
    total_area_strength::SVector{3,TF}
    root_strength::TF
    tip_strength::TF
    root_owned::Bool
    tip_owned::Bool
    handoff_active_before::Bool
    handoff_active_after::Bool
    attribution::Symbol
    # True on the conversion that deposits the sheet's physical trailing edge
    # (the starting vortex) whole, because no earlier conversion took a share.
    startup_edge_deposited::Bool
    expected_handoff::SVector{3,TF}
    # Exact filament content the transaction must transfer: sum_j H_j plus every
    # streamwise face jump times its edge vector. `residual_*` measures
    # `deposited_total - expected_total`, which is round-off by construction on
    # *any* geometry -- not only on a uniform affine fixture.
    expected_total::SVector{3,TF}
    deposited_total::SVector{3,TF}
    residual_abs::TF
    residual_rel::TF
    n_elided::Int
    capacity_requested::Int
    capacity_available::Int
    nearfield::Union{Nothing,SurfaceVorticityNearFieldDiagnostics{TF}}
end

"""
Live diagnostics for the surface-vorticity conversion: the per-conversion
records plus a cumulative ledger. Legacy wakes carry `nothing` instead.
"""
mutable struct SurfaceVorticityDiagnostics{TF}
    records::Vector{SurfaceVorticityConversionRecord{TF}}
    total_particles::Int
    total_interior::Int
    total_roottip::Int
    total_perimeter::Int
    total_area::TF
    total_area_strength::SVector{3,TF}
    total_root_strength::TF
    total_tip_strength::TF
    total_expected_handoff::SVector{3,TF}
    total_expected::SVector{3,TF}
    total_deposited::SVector{3,TF}
    total_residual::SVector{3,TF}
    total_elided::Int
end

SurfaceVorticityDiagnostics{TF}() where {TF} = SurfaceVorticityDiagnostics{TF}(
    SurfaceVorticityConversionRecord{TF}[], 0, 0, 0, 0,
    zero(TF), zero(SVector{3,TF}), zero(TF), zero(TF),
    zero(SVector{3,TF}), zero(SVector{3,TF}), zero(SVector{3,TF}),
    zero(SVector{3,TF}), 0)

function _accumulate!(diag::SurfaceVorticityDiagnostics{TF},
                      rec::SurfaceVorticityConversionRecord{TF}) where {TF}
    push!(diag.records, rec)
    diag.total_particles += rec.n_interior + rec.n_roottip + rec.n_perimeter
    diag.total_interior += rec.n_interior
    diag.total_roottip += rec.n_roottip
    diag.total_perimeter += rec.n_perimeter
    diag.total_area += rec.total_area
    diag.total_area_strength += rec.total_area_strength
    diag.total_root_strength += rec.root_strength
    diag.total_tip_strength += rec.tip_strength
    diag.total_expected_handoff += rec.expected_handoff
    diag.total_expected += rec.expected_total
    diag.total_deposited += rec.deposited_total
    diag.total_residual += rec.deposited_total - rec.expected_total
    diag.total_elided += rec.n_elided
    return diag
end

#--- Particle Maintenance Policies ---#

abstract type AbstractParticleTrimPolicy end
abstract type AbstractParticleFunctionalPolicy end

struct ParticleMaintenance{TT<:Tuple,TF<:Tuple}
    trim_policies::TT
    functional_policies::TF
end

ParticleMaintenance() = ParticleMaintenance((), ())
ParticleMaintenance(maintenance::ParticleMaintenance) = maintenance
ParticleMaintenance(policy::Union{AbstractParticleTrimPolicy,AbstractParticleFunctionalPolicy}) =
    ParticleMaintenance((policy,))

function ParticleMaintenance(policies::Tuple)
    trim_policies, functional_policies = _split_particle_policies(policies)
    return ParticleMaintenance(trim_policies, functional_policies)
end

function _split_particle_policies(policies::Tuple)
    isempty(policies) && return ((), ())

    policy = first(policies)
    trim_policies, functional_policies = _split_particle_policies(Base.tail(policies))

    if policy isa AbstractParticleTrimPolicy
        return ((policy, trim_policies...), functional_policies)
    elseif policy isa AbstractParticleFunctionalPolicy
        return (trim_policies, (policy, functional_policies...))
    else
        throw(ArgumentError("Invalid particle maintenance policy of type $(typeof(policy))"))
    end
end

struct ParticleMaintenanceContext{TF,TD}
    frames::TF
    step::Int
    dt::TD
end

struct MinGamma{T} <: AbstractParticleTrimPolicy
    threshold::T
end

struct MaxGamma{T} <: AbstractParticleTrimPolicy
    threshold::T
end

struct RelativeMinGamma{T} <: AbstractParticleTrimPolicy
    fraction::T
    function RelativeMinGamma{T}(fraction) where {T}
        fraction >= zero(fraction) || throw(ArgumentError("RelativeMinGamma fraction must be nonnegative"))
        return new{T}(fraction)
    end
end

RelativeMinGamma(fraction::T) where {T} = RelativeMinGamma{T}(fraction)

struct PreparedRelativeMinGamma{T} <: AbstractParticleTrimPolicy
    absolute_threshold::T
end

MinGammaPolicy(threshold; relative=true) =
    relative ? RelativeMinGamma(threshold) : MinGamma(threshold)

struct GlobalBox{T} <: AbstractParticleTrimPolicy
    xmin::SVector{3,T}
    xmax::SVector{3,T}
end

function GlobalBox(xmin, xmax)
    xmin_s = SVector{3}(xmin)
    xmax_s = SVector{3}(xmax)
    T = promote_type(eltype(xmin_s), eltype(xmax_s))
    return GlobalBox(SVector{3,T}(xmin_s), SVector{3,T}(xmax_s))
end

struct GlobalCylinder{T} <: AbstractParticleTrimPolicy
    origin::SVector{3,T}
    extrude::SVector{3,T}
    radius::T
end

function GlobalCylinder(origin, extrude, radius)
    origin_s = SVector{3}(origin)
    extrude_s = SVector{3}(extrude)
    T = promote_type(eltype(origin_s), eltype(extrude_s), typeof(radius))
    origin_t = SVector{3,T}(origin_s)
    extrude_t = SVector{3,T}(extrude_s)
    radius_t = T(radius)

    radius_t >= zero(T) || throw(ArgumentError("GlobalCylinder radius must be nonnegative"))
    sum(abs2, extrude_t) > zero(T) || throw(ArgumentError("GlobalCylinder extrude vector must be nonzero"))

    return GlobalCylinder(origin_t, extrude_t, radius_t)
end

struct FrameBox{T} <: AbstractParticleTrimPolicy
    i_frame::Int
    xmin::SVector{3,T}
    xmax::SVector{3,T}
end

function FrameBox(i_frame::Integer, xmin, xmax)
    xmin_s = SVector{3}(xmin)
    xmax_s = SVector{3}(xmax)
    T = promote_type(eltype(xmin_s), eltype(xmax_s))
    return FrameBox(Int(i_frame), SVector{3,T}(xmin_s), SVector{3,T}(xmax_s))
end

struct PreparedFrameBox{T,TO,TR} <: AbstractParticleTrimPolicy
    origin_global::TO
    Rg2f::TR
    xmin::SVector{3,T}
    xmax::SVector{3,T}
end

struct MergeParticles{TR,TH,TM} <: AbstractParticleFunctionalPolicy
    every::Int
    r::TR
    r_hash::TH
    sigma_relative::Bool
    max_sigma_ratio::TM
    skip_static::Bool
end

MergeParticles(; every, r=0.5, r_hash=-1.0, sigma_relative=true,
    max_sigma_ratio=2.0, skip_static=true) =
    MergeParticles(Int(every), r, r_hash, sigma_relative, max_sigma_ratio,
                   skip_static)

struct SplitParticles{TO} <: AbstractParticleFunctionalPolicy
    every::Int
    opts::TO
    verbose::Bool
end

SplitParticles(opts; every=1, verbose=false) =
    SplitParticles(Int(every), opts, verbose)

prepare_particle_policy(policy, pfield, ::ParticleMaintenanceContext) = policy

function prepare_particle_policy(policy::FrameBox, pfield, ctx::ParticleMaintenanceContext)
    isnothing(ctx.frames) && throw(ArgumentError("FrameBox particle trimming requires frames"))
    transform = frame_global_transform(ctx.frames, policy.i_frame)
    isnothing(transform) && throw(ArgumentError("FrameBox references unknown frame index $(policy.i_frame)"))
    origin_global, R_f2g = transform
    return PreparedFrameBox(origin_global, R_f2g', policy.xmin, policy.xmax)
end

function prepare_particle_policy(policy::RelativeMinGamma, pfield, ::ParticleMaintenanceContext)
    gmax = zero(typeof(policy.fraction))
    for i in 1:pfield.np
        g = _particle_gamma_magnitude(pfield, i)
        g > gmax && (gmax = g)
    end
    return PreparedRelativeMinGamma(policy.fraction * gmax)
end

prepare_particle_policies(policies::Tuple, pfield, ctx::ParticleMaintenanceContext) =
    Tuple(prepare_particle_policy(policy, pfield, ctx) for policy in policies)

function _particle_gamma_magnitude(pfield, i)
    gamma = FLOWVPM.get_Gamma(pfield, i)
    return sqrt(sum(abs2, gamma))
end

keep(policy::MinGamma, pfield, i, ::ParticleMaintenanceContext) =
    _particle_gamma_magnitude(pfield, i) >= policy.threshold

keep(policy::MaxGamma, pfield, i, ::ParticleMaintenanceContext) =
    _particle_gamma_magnitude(pfield, i) <= policy.threshold

keep(policy::PreparedRelativeMinGamma, pfield, i, ::ParticleMaintenanceContext) =
    _particle_gamma_magnitude(pfield, i) >= policy.absolute_threshold

function keep(policy::GlobalBox, pfield, i, ::ParticleMaintenanceContext)
    x = FLOWVPM.get_X(pfield, i)
    return all(policy.xmin .<= x .<= policy.xmax)
end

function keep(policy::GlobalCylinder, pfield, i, ::ParticleMaintenanceContext)
    x = SVector{3}(FLOWVPM.get_X(pfield, i))
    dx = x - policy.origin
    axis_length2 = sum(abs2, policy.extrude)
    axial = dot(dx, policy.extrude) / axis_length2
    zero(axial) <= axial <= one(axial) || return false
    closest = policy.origin + axial * policy.extrude
    return sum(abs2, x - closest) <= policy.radius^2
end

function keep(policy::PreparedFrameBox, pfield, i, ::ParticleMaintenanceContext)
    x = SVector{3}(FLOWVPM.get_X(pfield, i))
    x_local = policy.Rg2f * (x - policy.origin_global)
    return all(policy.xmin .<= x_local .<= policy.xmax)
end

function _keep_particle(policies::Tuple, pfield, i, ctx::ParticleMaintenanceContext)
    for policy in policies
        keep(policy, pfield, i, ctx) || return false
    end
    return true
end

function apply_particle_policy!(policy::MergeParticles, pfield, ctx::ParticleMaintenanceContext)
    if policy.every > 0 && ctx.step > 0 && ctx.step % policy.every == 0
        FLOWVPM.merge_particles!(pfield;
            r_merge=policy.r,
            r_hash=policy.r_hash,
            sigma_relative=policy.sigma_relative,
            max_sigma_ratio=policy.max_sigma_ratio,
            skip_static=policy.skip_static,
        )
    end
    return nothing
end

function apply_particle_policy!(policy::SplitParticles, pfield, ctx::ParticleMaintenanceContext)
    if policy.every > 0 && ctx.step > 0 && ctx.step % policy.every == 0
        FLOWVPM.split_particles!(pfield, policy.opts; dt=ctx.dt, verbose=policy.verbose)
    end
    return nothing
end

function apply_particle_policies!(policies::Tuple, pfield, ctx::ParticleMaintenanceContext)
    for policy in policies
        apply_particle_policy!(policy, pfield, ctx)
    end
    return nothing
end

function apply_particle_maintenance!(pfield, maintenance::ParticleMaintenance, ctx::ParticleMaintenanceContext)
    # Task 052: device-backed particle fields run maintenance on a host
    # mirror (trim policies and merge_particles! are host scalar-index code);
    # see src/FLOWPanel_gpu_wake.jl.
    if pfield isa FLOWVPM.ParticleField && !(pfield.particles isa Array)
        return _apply_particle_maintenance_device!(pfield, maintenance, ctx)
    end
    apply_particle_policies!(maintenance.functional_policies, pfield, ctx)
    prepared_trim_policies = prepare_particle_policies(maintenance.trim_policies, pfield, ctx)

    for i in pfield.np:-1:1
        if !_keep_particle(prepared_trim_policies, pfield, i, ctx)
            FLOWVPM.remove_particle(pfield, i)
        end
    end

    return nothing
end

#------- Vortex Particle Wake -------#

"""
    PanelParticleWake(body; optargs...)

Hybrid wake model that combines a near-body panel wake with vortex particles
shed downstream from the trailing edge. Panel wake options, including
`shed_with_induced_velocity` and `unsteady_filament`, are forwarded to the
internal [`PanelWake`](@ref). [`SurfaceVorticityConversion`](@ref) requires
`unsteady_filament=true` and `include_final_filament=true`; the legacy
conversion retains the full `PanelWake` option set.

`nwakerows = 0` (BRAINSTORM 024) selects **convert-at-shed**: the just-shed
row is converted to particles inside `shed_wake!` itself, so particles appear
at the TE+Das line and no free wake-panel row ever enters a solve. Requires
both final-filament options true (the retained filament on the Das line
carries the sheet/particle interface) and, with the smooth conversion,
`attribution = :downstream` (the only circulation-conserving attribution when
no sheet survives the shed).
"""
#------- relaxation spatial filters -------#

"""
    RelaxationPlaneFilter(point, normal; i_frame=0)

Per-particle predicate for [`FLOWVPM.Relaxation`](@ref) that restricts relaxation to
the half-space on the side of a plane the `normal` points toward. The plane is
defined by `point` and `normal` expressed in reference frame `i_frame` (use `0` for
the global frame). Particle `p` is relaxed when
`dot(get_X(p) - point_world, normal_world) >= 0`; flipping the sign of `normal`
relaxes the opposite side instead.

The world-space plane is cached and refreshed once per timestep via
[`refresh_relaxation_filter!`](@ref) so the per-particle test is a single dot
product. When `i_frame != 0`, the plane tracks the (possibly moving) frame.

Typical use: relax only particles that have propagated a distance `d` downstream of a
rotor plane. With the rotor frame origin at the hub and downstream `= -ẑ`, place the
threshold plane at `-d*ẑ` with `normal = -ẑ`:

    RelaxationPlaneFilter(SVector(0.0,0.0,-d), SVector(0.0,0.0,-1.0); i_frame=i_rotor)
"""
mutable struct RelaxationPlaneFilter{TF}
    point_local::SVector{3,TF}   # point on the plane, in frame i_frame
    normal_local::SVector{3,TF}  # unit normal; points toward the RELAXED side
    i_frame::Int                 # 0 ⇒ plane defined in the global frame
    point_world::SVector{3,TF}   # cached world-space plane point (refreshed per step)
    normal_world::SVector{3,TF}  # cached world-space plane normal
end

function RelaxationPlaneFilter(point, normal; i_frame::Int=0)
    length(point) == 3 || throw(ArgumentError("RelaxationPlaneFilter point must have length 3"))
    length(normal) == 3 || throw(ArgumentError("RelaxationPlaneFilter normal must have length 3"))
    p = SVector{3}(float.(point))
    n = SVector{3}(float.(normal))
    nrm = LA.norm(n)
    isfinite(nrm) && nrm > zero(nrm) ||
        throw(ArgumentError("RelaxationPlaneFilter normal must be nonzero and finite"))
    n = n / nrm
    TF = promote_type(eltype(p), eltype(n))
    p = SVector{3,TF}(p)
    n = SVector{3,TF}(n)
    # initialize the world cache to the local definition (exact when i_frame == 0,
    # and overwritten by refresh_relaxation_filter! for frame-relative planes)
    return RelaxationPlaneFilter{TF}(p, n, i_frame, p, n)
end

(f::RelaxationPlaneFilter)(p) =
    LA.dot(FLOWVPM.get_X(p) .- f.point_world, f.normal_world) >= 0

(f::RelaxationPlaneFilter)(pfield, i::Integer) =
    LA.dot(FLOWVPM.get_X(pfield, Int(i)) .- f.point_world, f.normal_world) >= 0

FLOWVPM._passes_relaxation_filter(f::RelaxationPlaneFilter, pfield, i::Integer) = f(pfield, i)

"""
    refresh_relaxation_filter!(filter, frames)

Update a relaxation filter's cached world-space geometry from the current frame tree.
A no-op for any filter that does not track a frame (the generic fallback), so it is
safe to call unconditionally each timestep.
"""
refresh_relaxation_filter!(::Any, frames) = nothing

function refresh_relaxation_filter!(f::RelaxationPlaneFilter, frames)
    f.i_frame == 0 && return nothing
    frames === nothing && throw(ArgumentError("RelaxationPlaneFilter requires frames"))
    transform = frame_global_transform(frames, f.i_frame)
    isnothing(transform) &&
        throw(ArgumentError("RelaxationPlaneFilter references unknown frame index $(f.i_frame)"))
    o, R = transform
    f.point_world = R * f.point_local + o
    f.normal_world = R * f.normal_local
    return nothing
end

"""
    plane_filtered_relaxation(base::FLOWVPM.Relaxation, point, normal; i_frame=0)

Return a copy of relaxation scheme `base` (same method, `rlxf`, and `nsteps_relax`)
whose per-particle filter is a [`RelaxationPlaneFilter`](@ref). Convenient for
attaching a spatial filter to a stock scheme, e.g.

    plane_filtered_relaxation(FLOWVPM.relaxation_correctedpedrizzetti,
        SVector(0.0,0.0,-d), SVector(0.0,0.0,-1.0); i_frame=i_rotor)
"""
plane_filtered_relaxation(base::FLOWVPM.Relaxation, point, normal; i_frame::Int=0) =
    FLOWVPM.Relaxation(base.relax, base.nsteps_relax, base.rlxf,
                       RelaxationPlaneFilter(point, normal; i_frame))

struct PanelParticleWake{TK,NK,TF,TPF,MT,MU,TPM,TNT,TC,TW,TD} <: AbstractFreeWake
    panel_wake::PanelWake{TK,NK,TF}
    pfield::TPF                           # FLOWVPM.ParticleField object
    method_trailing::MT                             # particle shedding method
    method_unsteady::MU                             # particle shedding method
    particle_maintenance::TPM             # particle merge/trim policy chain
    particle_core_size::Float64        # NaN uses source body core_size (regularization core radius)
    pfield_optargs::TNT                   # resolved FLOWVPM optargs (for reproduction metadata)
    # Panel-to-particle conversion strategy (BRAINSTORM 016). Legacy instances
    # carry `nothing` for the workspace and diagnostics and must not acquire any
    # new per-step allocation.
    conversion::TC
    conversion_workspace::TW
    conversion_diagnostics::TD
    conversion_count::Array{Int,0}        # successful smooth conversions committed
end

function PanelParticleWake(body::AbstractLiftingBody;
        nwakerows=3, max_particles=10000,
        conversion::AbstractPanelParticleConversion=LegacyEdgeJumpConversion(),
        method_trailing::WakeSheddingMethod=DefaultWakeSheddingMethod(),
        method_unsteady::WakeSheddingMethod=DefaultWakeSheddingMethod(),
        particle_maintenance=ParticleMaintenance(),
        particle_core_size::Union{Real,Nothing}=nothing,
        # deprecated alias for the pre-2026-08-22 name
        particle_kerneloffset::Union{Real,Nothing}=nothing,
        viscous=FLOWVPM.Inviscid(),
        SFS=FLOWVPM.SFS_default,
        unsteady_filament=true,
        include_final_filament=true,
        # Make the relaxation scheme explicit so it is recorded in metadata rather
        # than silently inherited from FLOWVPM. The default matches the FLOWVPM
        # ParticleField default (CorrectedPedrizzetti) to preserve prior behavior.
        relaxation=FLOWVPM.relaxation_correctedpedrizzetti,
        # Frozen-gradient geometric local integrator (BRAINSTORM 020 Phase 2R).
        # Default false = stock forward Euler, bit-identical.
        expint=false,
        # Task 052: array container for the particle field (pass CUDA.CuArray
        # for a device-resident wake) and an optional FMM settings override
        # (a device field REQUIRES all-off autotune; see FLOWVPM_fmm_radix.jl).
        arraytype=Matrix,
        pfield_fmm=FLOWVPM.FMM(autotune_reg_error=false),
        # BRAINSTORM/022 Ruling 7: multiple wakes may shed into ONE shared
        # particle field. Pass another wake's `.pfield` here; `nothing`
        # allocates a private field (legacy behavior). Callers that share a
        # field are responsible for sizing `max_particles` on the OWNING
        # allocation — the simulate! loop deduplicates shared fields by
        # identity, so per-step processing happens once regardless of how
        # many wakes reference it.
        pfield=nothing,
        kwargs...)

    # Resolve the legacy line policies against the selected strategy before any
    # allocation, so a configuration error costs nothing.
    trailing = _resolve_line_policy(conversion, method_trailing, "method_trailing")
    unsteady = _resolve_line_policy(conversion, method_unsteady, "method_unsteady")

    use_unsteady_filament = Bool(unsteady_filament)
    use_final_filament = Bool(include_final_filament)
    _validate_conversion_filaments(conversion, use_unsteady_filament,
                                   use_final_filament)

    # BRAINSTORM 024: nwakerows = 0 selects convert-at-shed — the just-shed
    # row is converted to particles inside shed_wake! itself (particles appear
    # at the Das line) and no free sheet ever enters a solve. Storage is the
    # N=1 layout with the convert_at_shed marker set.
    nwakerows >= 0 || throw(ArgumentError(
        "PanelParticleWake nwakerows must be >= 0 (got $(nwakerows))"))
    convert_at_shed = nwakerows == 0
    if convert_at_shed
        # The retained final filament on the TE+Das line is the only carrier
        # of the sheet/particle interface once no sheet survives: it cannot
        # be opted out of at N=0.
        use_unsteady_filament && use_final_filament || throw(ArgumentError(
            "nwakerows = 0 (convert-at-shed) requires unsteady_filament=true " *
            "and include_final_filament=true: with no surviving sheet the " *
            "final filament on the TE+Das line is the sheet/particle " *
            "interface carrier"))
        # With the row converted in the same solve that produced it, the
        # upstream (rigid-row) face jump is identically zero and the per-step
        # unsteady circulation lives on the downstream face; only :downstream
        # attribution deposits it. :upstream would silently drop ALL unsteady
        # circulation and :split would strand half on a filament row that no
        # longer exists.
        conversion isa SurfaceVorticityConversion &&
            conversion.attribution != :downstream && throw(ArgumentError(
                "nwakerows = 0 (convert-at-shed) requires " *
                "SurfaceVorticityConversion attribution=:downstream (got " *
                ":$(conversion.attribution)): it is the only " *
                "circulation-conserving attribution when no sheet survives " *
                "the shed"))
    end

    panel_wake = PanelWake(body; nwakerows=(convert_at_shed ? 1 : nwakerows),
        unsteady_filament=use_unsteady_filament,
        include_final_filament=use_final_filament,
        convert_at_shed, kwargs...)
    TF = FastMultipole.numtype(panel_wake)

    # Create particle field with default settings (disable autotune_reg_error to avoid convergence issues)
    # `integration` must name the scheme this wake actually steps with
    # (`FLOWVPM._euler`, see `step!` below): `viscousdiffusion` branches on
    # `pfield.integration`, and under the FLOWVPM default (`rungekutta3`) a
    # CoreSpreading scheme hits the RK3 branch with zeroed stage weights —
    # no core spreading, no beta resets, silently inviscid.
    if pfield === nothing
        pfield = FLOWVPM.ParticleField(max_particles, TF;
            viscous,
            fmm=pfield_fmm,
            SFS,
            integration=(Bool(expint) ? FLOWVPM.euler_exp : FLOWVPM.euler),
            relaxation,
            arraytype)
    else
        # Shared-field path (Ruling 7): the supplied field must already step
        # with the scheme this wake expects — `propagate!` branches on
        # `pfield.integration`, and a numtype mismatch would silently promote
        # or truncate shed strengths.
        eltype(pfield.particles) == TF || throw(ArgumentError(
            "shared pfield numtype $(eltype(pfield.particles)) does not " *
            "match this wake's numtype $(TF)"))
        expected_integration = Bool(expint) ? FLOWVPM.euler_exp : FLOWVPM.euler
        pfield.integration === expected_integration || throw(ArgumentError(
            "shared pfield integration $(pfield.integration) does not match " *
            "this wake's expint choice ($(expected_integration))"))
    end

    # Capture the resolved FLOWVPM construction options for reproduction metadata.
    # Read back from the live particle field so the recorded values are authoritative
    # (they reflect actual state, including any FLOWVPM-side defaults).
    pfield_optargs = (;
        viscous     = pfield.viscous,
        SFS         = pfield.SFS,
        relaxation  = pfield.relaxation,
        formulation = pfield.formulation,
        kernel      = pfield.kernel,
        integration = pfield.integration,
    )

    # Infer type params from the actual panel_wake
    WTK = typeof(panel_wake).parameters[1]
    WNK = typeof(panel_wake).parameters[2]
    maintenance = ParticleMaintenance(particle_maintenance)
    particle_core_size = _core_size_alias(particle_core_size, particle_kerneloffset,
                                          NaN, :particle_core_size, :particle_kerneloffset)
    if !isnan(particle_core_size)
        body.core_size_targets = Float64(particle_core_size)
    end
    workspace = _make_conversion_workspace(conversion, panel_wake, TF)
    diagnostics = _make_conversion_diagnostics(conversion, TF)
    conversion_count = Array{Int,0}(undef)
    conversion_count[] = 0

    return PanelParticleWake{WTK,WNK,TF,typeof(pfield),typeof(trailing),typeof(unsteady),typeof(maintenance),typeof(pfield_optargs),typeof(conversion),typeof(workspace),typeof(diagnostics)}(
        panel_wake, pfield, trailing, unsteady, maintenance, Float64(particle_core_size), pfield_optargs,
        conversion, workspace, diagnostics, conversion_count
    )
end

# Legacy conversion allocates no workspace and no diagnostics.
_make_conversion_workspace(::LegacyEdgeJumpConversion, panel_wake, ::Type{TF}) where {TF} = nothing
_make_conversion_diagnostics(::LegacyEdgeJumpConversion, ::Type{TF}) where {TF} = nothing

function _make_conversion_workspace(conversion::SurfaceVorticityConversion,
                                    panel_wake, ::Type{TF}) where {TF}
    ws = SurfaceVorticityWorkspace{TF}()

    # Pre-size the staging buffers from the geometry that will actually be
    # converted, so the steady-state transaction never grows them. This is the
    # explicit pre-transaction capacity path of Phase 2 sec. 7.
    h = TF(conversion.sigma) / TF(conversion.overlap)
    n_est = 0
    N = size(panel_wake.nodes[1], 2) - 1   # full-buffer outgoing row index
    for i_surf in eachindex(panel_wake.nodes)
        nodes = panel_wake.nodes[i_surf]
        n_cols = size(nodes, 3) - 1
        for j in 1:n_cols
            v1 = _wake_node(nodes, N, j)
            v2 = _wake_node(nodes, N + 1, j)
            v3 = _wake_node(nodes, N + 1, j + 1)
            v4 = _wake_node(nodes, N, j + 1)
            # A freshly constructed wake has all-zero nodes, in which case this
            # estimate is simply the floor of one subcell per panel; the buffers
            # grow once on the first real conversion and are reused thereafter.
            n_xi, n_eta = _subdivision_counts(v1, v2, v3, v4, h)
            n_est += n_xi * n_eta
        end
        n_est += 2 * max(1, ceil(Int, TF(conversion.overlap)))  # root/tip closure
    end

    for buf in (ws.positions, ws.strengths, ws.sigmas, ws.circulations, ws.classes)
        sizehint!(buf, n_est)
    end
    return ws
end

function _make_conversion_diagnostics(conversion::SurfaceVorticityConversion,
                                      ::Type{TF}) where {TF}
    return SurfaceVorticityDiagnostics{TF}()
end

"""
Run SFS pre-calculations for particle field before evaluating the velocity field.
"""
function pre_evaluate_influence!(pfield::FLOWVPM.ParticleField)
    if !FLOWVPM.isSFSenabled(pfield.SFS)
        pfield.SFS(pfield, FLOWVPM.BeforeUJ())
        return nothing
    end

    _step_timer_measure(:wake_sfs; nested=true) do
        velocities = copy(view(pfield.particles, FLOWVPM.U_INDEX, 1:pfield.np))
        pfield.SFS(pfield, FLOWVPM.BeforeUJ())
        FLOWVPM._reset_particles(pfield)
        pfield.particles[FLOWVPM.U_INDEX, 1:pfield.np] .= velocities
    end
    return nothing
end

function post_evaluate_influence!(pfield::FLOWVPM.ParticleField,
        source::FLOWVPM.ParticleField, backend::FastMultipoleBackend, outputs;
        i_target::Int=1, i_source::Int=1)
    pfield === source || return nothing
    FLOWVPM.isSFSenabled(pfield.SFS) || return nothing

    _step_timer_measure(:wake_sfs; nested=true) do
        _, _, target_tree, source_tree, _, direct_list, _ = outputs
        FLOWVPM.Estr_fmm!(pfield, pfield, target_tree, source_tree, direct_list;
            i_target_system=i_target, i_source_system=i_source)
        pfield.SFS(pfield, FLOWVPM.AfterUJ())
    end
    return nothing
end

function post_evaluate_influence!(pfield::FLOWVPM.ParticleField,
        source::FLOWVPM.ParticleField, backend::DirectBackend, outputs;
        i_target::Int=1, i_source::Int=1)
    pfield === source || return nothing
    FLOWVPM.isSFSenabled(pfield.SFS) || return nothing

    _step_timer_measure(:wake_sfs; nested=true) do
        FLOWVPM.Estr_direct!(pfield)
        pfield.SFS(pfield, FLOWVPM.AfterUJ())
    end
    return nothing
end

#--- Delegation methods ---#

get_probes(w::PanelParticleWake) = (get_probes(w.panel_wake)..., w.pfield)
get_sources(w::PanelParticleWake) = (get_sources(w.panel_wake)..., w.pfield)

function reset!(w::PanelParticleWake)
    # reset panel wake
    reset!(w.panel_wake)

    # reset particle velocity and Jacobian fields (preserve position and strength)
    FLOWVPM._reset_particles(w.pfield)
    
    # reset particle SFS properties
    FLOWVPM._reset_particles_sfs(w.pfield)
end

function apply_freestream!(w::PanelParticleWake, uinf; include_pfield::Bool=true)
    # apply to panel wake
    apply_freestream!(w.panel_wake, uinf)

    # Ruling 7: when several wakes share one particle field, the caller passes
    # include_pfield=false for every wake after the first so the freestream is
    # ADDED to particle velocities exactly once.
    include_pfield || return nothing

    # Add freestream to particle velocities (broadcast on device-backed
    # fields, task 052; the scalar loop is kept for the host path)
    if !(w.pfield.particles isa Array)
        _gpu_apply_freestream_device!(w.pfield, uinf)
        return nothing
    end
    for i in 1:w.pfield.np
        for d in 1:3
            w.pfield.particles[FLOWVPM.U_INDEX[d], i] += uinf[d]
        end
    end
end

update_TE!(w::PanelParticleWake, sys) = update_TE!(w.panel_wake, sys)

function _particle_gamma_direction_stats(pfield::FLOWVPM.ParticleField;
        vertical=(0.0, 0.0, 1.0), before_gamma=nothing)
    np = pfield.np
    np == 0 && return "np=0"

    vertical_norm = LA.norm(vertical)
    zhat = vertical_norm > 0 ? SVector{3}(vertical...) / vertical_norm : SVector(0.0, 0.0, 1.0)
    mean_abs = zeros(Float64, 3)
    mean_dot_vertical = 0.0
    mean_angle_change = 0.0
    n_angle = 0

    for i in 1:np
        gamma = SVector{3}(view(pfield.particles, FLOWVPM.GAMMA_INDEX, i))
        gamma_norm = LA.norm(gamma)
        gamma_norm > 0 || continue
        ghat = gamma / gamma_norm
        mean_abs .+= abs.(ghat)
        mean_dot_vertical += dot(ghat, zhat)

        if !isnothing(before_gamma) && i <= size(before_gamma, 2)
            gamma0 = SVector{3}(view(before_gamma, :, i))
            gamma0_norm = LA.norm(gamma0)
            if gamma0_norm > 0
                ghat0 = gamma0 / gamma0_norm
                mean_angle_change += acos(clamp(dot(ghat0, ghat), -1.0, 1.0))
                n_angle += 1
            end
        end
    end

    mean_abs ./= np
    mean_dot_vertical /= np
    angle = n_angle == 0 ? NaN : mean_angle_change / n_angle
    return "np=$(np) mean_abs_gammahat=($(mean_abs[1]), $(mean_abs[2]), $(mean_abs[3])) mean_gammahat_dot_vertical=$(mean_dot_vertical) mean_angle_change=$(angle)"
end

function propagate!(w::PanelParticleWake, dt; relax=true, step=0, frames=nothing,
        diagnose_particle_gamma::Bool=false, diagnostic_vertical=(0.0, 0.0, 1.0),
        # Ruling 7: when several wakes share one particle field, the caller
        # passes propagate_pfield=false for every wake after the first so the
        # field is convected and maintained exactly once per step. The panel
        # wake is per-body and always propagates.
        propagate_pfield::Bool=true)

    # panel wake
    propagate!(w.panel_wake, dt)

    propagate_pfield || return nothing

    gamma_before = diagnose_particle_gamma && w.pfield.np > 0 ?
        copy(view(w.pfield.particles, FLOWVPM.GAMMA_INDEX, 1:w.pfield.np)) : nothing
    if diagnose_particle_gamma
        println("particle gamma step=$(step) phase=before_euler " *
            _particle_gamma_direction_stats(w.pfield; vertical=diagnostic_vertical))
    end

    # refresh any frame-tracking relaxation filter to the current frame pose
    refresh_relaxation_filter!(w.pfield.relaxation.filter, frames)

    # convect particles (`pfield.integration` is the single source of truth:
    # stock forward Euler unless the wake was built with `expint=true`)
    _step_timer_measure(:wake_convection; nested=true) do
        if w.pfield.integration === FLOWVPM.euler_exp
            FLOWVPM._euler_exp(w.pfield, dt; relax)
        else
            FLOWVPM._euler(w.pfield, dt; relax)
        end
    end

    if diagnose_particle_gamma
        println("particle gamma step=$(step) phase=after_euler relax=$(relax) " *
            _particle_gamma_direction_stats(w.pfield;
                vertical=diagnostic_vertical, before_gamma=gamma_before))
    end

    # particle maintenance
    _step_timer_measure(:wake_maintenance; nested=true) do
        apply_particle_maintenance!(w.pfield, w.particle_maintenance,
            ParticleMaintenanceContext(frames, step, dt))
    end
end

function write_vtk(name, w::PanelParticleWake, idx, t; overwrite=false, compress::Bool=true,
        # Ruling 7: with a shared pfield the caller passes include_pfield=false
        # for every wake after the first so the particle cloud is written once
        # per step instead of once per referencing wake.
        include_pfield::Bool=true)
    # panel wake (includes filaments)
    write_vtk(name, w.panel_wake, idx, t; overwrite, compress)

    include_pfield || return nothing

    # particle wake — route block files to subdirectory
    vpm_path, vpm_name = splitdir(name)
    particles_pvd_name = joinpath(vpm_path, vpm_name * "_particles")
    particles_subdir = joinpath(vpm_path, vpm_name * "_particles")
    mkpath(particles_subdir)
    particles_block = joinpath(particles_subdir, vpm_name * "_particles")

    np = w.pfield.np
    # Task 052: WriteVTK needs host arrays; on a device-backed field this is
    # the refreshed host mirror (D2H of the live prefix), else a no-op alias.
    host_particles = _wake_monitor_host_pfield(w.pfield).particles
    # Task 052c io fix (Ryan-approved 2026-08-26, single-series ruling same
    # day): ONE uncompressed particle series at a user-chosen precision
    # (zlib costs ~10x the raw write for ~8% size on this data, so it is
    # never used). FLOWPANEL_PARTICLE_PRECISION=f64 (default) keeps warm
    # restarts replay-exact; f32 halves write time and disk for
    # visualization-focused runs at the cost of exact continuation (the
    # warm-start loader warns when restoring from Float32).
    cells = [WriteVTK.MeshCell(WriteVTK.PolyData.Verts(), 1:np)]

    _write_particles_vtp(particles_block * ".$idx.vtp", host_particles, np,
        cells, _particle_vtp_eltype(host_particles))

    vtp_relpath = joinpath(vpm_name * "_particles", vpm_name * "_particles.$idx.vtp")
    _pvd_append!(particles_pvd_name * ".pvd", t, vtp_relpath; overwrite)
end

# Element type of the particle series per FLOWPANEL_PARTICLE_PRECISION:
# "f64" (default) = the pfield eltype (replay-exact restarts), "f32" =
# Float32 (half the write time and disk; restarts are not replay-exact).
function _particle_vtp_eltype(host_particles)
    precision = lowercase(get(ENV, "FLOWPANEL_PARTICLE_PRECISION", "f64"))
    precision in ("f32", "f64") || error(
        "FLOWPANEL_PARTICLE_PRECISION must be f32 or f64, got: $precision")
    return precision == "f32" ? Float32 : eltype(host_particles)
end

# One particle-cloud .vtp at the requested element type (see
# FLOWPANEL_PARTICLE_PRECISION above).
# Always uncompressed; conversion is skipped when the data already matches.
function _write_particles_vtp(filename, host_particles, np, cells, ::Type{T}) where T
    _conv(a) = eltype(a) === T ? a : T.(a)
    X = _conv(view(host_particles, FLOWVPM.X_INDEX, 1:np))
    vtp = WriteVTK.vtk_grid(filename, X, cells; compress=false)

    if np > 0
        vtp["gamma", WriteVTK.VTKPointData()] = _conv(view(host_particles, FLOWVPM.GAMMA_INDEX, 1:np))
        vtp["sigma", WriteVTK.VTKPointData()] = _conv(view(host_particles, FLOWVPM.SIGMA_INDEX, 1:np))
        vtp["vol", WriteVTK.VTKPointData()] = _conv(view(host_particles, FLOWVPM.VOL_INDEX, 1:np))
        vtp["circulation", WriteVTK.VTKPointData()] = _conv(view(host_particles, FLOWVPM.CIRCULATION_INDEX, 1:np))
        vtp["velocity", WriteVTK.VTKPointData()] = _conv(view(host_particles, FLOWVPM.U_INDEX, 1:np))
        vtp["vorticity", WriteVTK.VTKPointData()] = _conv(view(host_particles, FLOWVPM.VORTICITY_INDEX, 1:np))
        vtp["C", WriteVTK.VTKPointData()] = _conv(view(host_particles, FLOWVPM.C_INDEX, 1:np))
        vtp["SFS", WriteVTK.VTKPointData()] = _conv(view(host_particles, FLOWVPM.SFS_INDEX, 1:np))
        vtp["velocity_gradient", WriteVTK.VTKPointData()] = reshape(_conv(view(host_particles, FLOWVPM.J_INDEX, 1:np)), 3, 3, np)
    end

    WriteVTK.vtk_save(vtp)
    return nothing
end

requires_hessian(::FLOWVPM.ParticleField) = true

# --- FastMultipole compatibility ---

# Ensure ParticleField has numtype for FastMultipole integration
FastMultipole.numtype(pf::FLOWVPM.ParticleField) = eltype(pf)

# --- Save and convert last row to particles ---

"""
    _convert_to_particles!(wake::PanelParticleWake, system=nothing)

Convert the outgoing panel-wake row into particles using the wake's configured
[`AbstractPanelParticleConversion`](@ref) strategy. Called by `shed_wake!`
before the panel rows are shifted, so the outgoing row is the current final
active row.

`system` is the shedding body. The legacy strategy ignores it;
[`SurfaceVorticityConversion`](@ref) needs it only when the wake holds a single
panel row, where the upstream half of the streamwise stencil is the row
`shed_wake!` is about to create (see the transaction below).
"""
_convert_to_particles!(wake::PanelParticleWake, system=nothing) =
    _convert_to_particles!(wake, wake.conversion, system)

function _convert_to_particles!(wake::PanelParticleWake, ::LegacyEdgeJumpConversion,
                                system)

    nwakes = wake.panel_wake.nwakes[]
    
    for i_surf in eachindex(wake.panel_wake.nodes)
        nodes = wake.panel_wake.nodes[i_surf]
        strength = wake.panel_wake.strength[i_surf]
        n_cols = size(nodes, 3) - 1
        _validate_station_method(wake.method_trailing, i_surf, n_cols + 1)
        _validate_station_method(wake.method_unsteady, i_surf, n_cols + 1)

        # check if this surface wraps on itself
        r1_le = SVector{3}(nodes[1, nwakes, 1], nodes[2, nwakes, 1], nodes[3, nwakes, 1])
        rend_le = SVector{3}(nodes[1, nwakes, n_cols+1], nodes[2, nwakes, n_cols+1], nodes[3, nwakes, n_cols+1])
        wraps = norm(r1_le - rend_le) < 5*eps()
        Γ_last = wraps ? strength[1,nwakes,n_cols] : zero(eltype(strength))

        for icol in 1:n_cols

            # left trailing filament (streamwise direction)
            Γ = strength[1,nwakes,icol]
            r1_le = SVector{3}(nodes[1, nwakes, icol], nodes[2, nwakes, icol], nodes[3, nwakes, icol])
            r1_te = SVector{3}(nodes[1, nwakes+1, icol], nodes[2, nwakes+1, icol], nodes[3, nwakes+1, icol])
            _shed_particles!(wake.pfield, r1_le, r1_te, Γ-Γ_last,
                _station_method(wake.method_trailing, i_surf, icol, icol))

            # unsteady particle
            r2_te = SVector{3}(nodes[1, nwakes+1, icol+1], nodes[2, nwakes+1, icol+1], nodes[3, nwakes+1, icol+1])
            Γ_tm1 = strength[1,nwakes+1,icol]
            _shed_particles!(wake.pfield, r1_te, r2_te, Γ-Γ_tm1,
                _station_method(wake.method_unsteady, i_surf, icol, icol+1))

            Γ_last = Γ
        end

        if !wraps
            # right trailing particles (streamwise direction)
            Γ = strength[1,nwakes,n_cols]
            r1_le = SVector{3}(nodes[1, nwakes, n_cols+1], nodes[2, nwakes, n_cols+1], nodes[3, nwakes, n_cols+1])
            r1_te = SVector{3}(nodes[1, nwakes+1, n_cols+1], nodes[2, nwakes+1, n_cols+1], nodes[3, nwakes+1, n_cols+1])
            _shed_particles!(wake.pfield, r1_le, r1_te, -Γ,
                _station_method(wake.method_trailing, i_surf, n_cols+1, n_cols+1))
        end
    end
end

#--- Surface-vorticity conversion transaction (BRAINSTORM 016 Phase 3) ---#

"Vertex `(irow, icol)` of a wake node array as a static vector."
@inline _wake_node(nodes, irow, icol) = SVector{3}(
    nodes[1, irow, icol], nodes[2, irow, icol], nodes[3, irow, icol])

"Centroid of the wake panel occupying node rows `irow, irow+1` and node
columns `icol, icol+1`."
@inline _wake_panel_centroid(nodes, irow, icol) = 0.25 * (
    _wake_node(nodes, irow, icol) + _wake_node(nodes, irow + 1, icol) +
    _wake_node(nodes, irow + 1, icol + 1) + _wake_node(nodes, irow, icol + 1))

"""
Stage the particles of one line segment into the conversion workspace, using
exactly the `SigmaOverlap(sigma, overlap)` sampling rule (and the exact-zero
elision of `_shed_particles!`), but writing to the preflight buffers instead of
the live particle field.

Returns `(n_staged, n_elided, staged_vector_strength)`.
"""
function _stage_line_particles!(ws::SurfaceVorticityWorkspace{TF}, r1, r2, Γ,
                                sigma, overlap, class) where {TF}
    dist = LA.norm(r2 - r1)
    dist < eps(TF) && return 0, 0, zero(SVector{3,TF})
    p_per_step = max(1, ceil(Int, overlap * dist / sigma))
    # Exact-zero strength is elided, matching the `_shed_particles!` guard that
    # exists because a zero |Γ| divides corrected-Pedrizzetti relaxation.
    Γ == 0 && return 0, p_per_step, zero(SVector{3,TF})

    distance_vector = (r2 - r1) / p_per_step
    Xp = r1 + distance_vector * TF(0.5)
    Γp = Γ * distance_vector
    for _ in 1:p_per_step
        _stage_particle!(ws, Xp, Γp, sigma, Γ, class)
        Xp += distance_vector
    end
    return p_per_step, 0, p_per_step * Γp
end

function _nearfield_summary(values::Vector{TF}) where {TF}
    n = length(values)
    if n == 0
        nan = TF(NaN)
        return SurfaceVorticityNearFieldSummary{TF}(0, nan, nan, nan, nan, nan, nan)
    end
    ordered = sort(values)
    avg = sum(ordered) / n
    rms = sqrt(sum(abs2, ordered) / n)
    mid = isodd(n) ? ordered[(n + 1) ÷ 2] :
          (ordered[n ÷ 2] + ordered[n ÷ 2 + 1]) / 2
    p95 = ordered[clamp(ceil(Int, TF(0.95) * n), 1, n)]
    return SurfaceVorticityNearFieldSummary{TF}(
        n, first(ordered), last(ordered), avg, rms, mid, p95)
end

"""
Evaluate the retained panel-wake sources, including the final filament, at the
staged nonzero candidates. This uses a fresh probe system and `DirectBackend`,
so existing particles and particle-particle interactions cannot enter the
diagnostic. The caller invokes it during preflight, before any insertion.
"""
function _surface_vorticity_nearfield(pw::PanelWake,
        ws::SurfaceVorticityWorkspace{TF}) where {TF}
    n = _n_pending(ws)
    for i in 1:n
        all(isfinite, ws.positions[i]) && all(isfinite, ws.strengths[i]) ||
            throw(WakeGeometryError(
                "non-finite staged conversion candidate $(i) in near-field preflight"))
    end
    probes = FastMultipole.ProbeSystem(n, TF)
    zero_v = zero(SVector{3,TF})
    zero_h = zero(SMatrix{3,3,TF,9})
    for i in 1:n
        probes.position[i] = ws.positions[i]
        probes.scalar_potential[i] = zero(TF)
        probes.gradient[i] = zero_v
        probes.hessian[i] = zero_h
    end
    n > 0 && influence!((probes,), get_sources(pw), DirectBackend();
        scalar_potential=false, gradient=false, hessian=(true,))

    interior = TF[]
    roottip = TF[]
    perimeter = TF[]
    sizehint!(interior, count(==(InteriorSurfaceParticle), ws.classes))
    sizehint!(roottip, count(==(RootTipSurfaceParticle), ws.classes))
    sizehint!(perimeter, count(==(PerimeterLineParticle), ws.classes))
    for i in 1:n
        value = LA.norm(probes.hessian[i])
        isfinite(value) || throw(WakeGeometryError(
            "non-finite panel-induced velocity gradient at staged conversion candidate $(i)"))
        if ws.classes[i] == InteriorSurfaceParticle
            push!(interior, value)
        elseif ws.classes[i] == RootTipSurfaceParticle
            push!(roottip, value)
        else
            push!(perimeter, value)
        end
    end
    return SurfaceVorticityNearFieldDiagnostics{TF}(
        _nearfield_summary(interior), _nearfield_summary(roottip),
        _nearfield_summary(perimeter))
end

"""
    _convert_to_particles!(wake, conversion::SurfaceVorticityConversion, system)

Smooth surface-vorticity conversion of the outgoing wake row (BRAINSTORM 016).

The outgoing ghost row is the current final active row `N`; conversion runs
before `shed_wake!` shifts the rows, so the ghost and its upstream neighbour are
both readable in place (Phase 2 sec. 13.1 Option B). For each ghost panel the
tangential gradient of the stored strength is reconstructed from a two-leg
centroid stencil -- streamwise to the upstream row, spanwise to the neighbouring
columns -- and the surface vorticity
``\\boldsymbol\\kappa = -\\boldsymbol n\\times\\nabla_s\\hat\\mu`` is deposited as
area-weighted particles ``\\boldsymbol\\Gamma_p=\\boldsymbol\\kappa\\,\\Delta A_p``
on a bilinear subdivision at target spacing ``h=\\sigma/\\mathrm{overlap}``.

Internal panel-to-panel edges are *not* reproduced as filaments; only the true
root and tip streamwise closures of an open chain survive, sampled with the same
orientation and signed strength as the legacy path so the two strategies differ
only in how the *area* is represented. A wrapping chain has no closure at all.

The whole call is a preflight/commit transaction: nothing in the particle field
or the panel rows is mutated until every panel has been reconstructed, every
candidate materialized, and the exact particle count checked against the
remaining `ParticleField` capacity.
"""
function _convert_to_particles!(wake::PanelParticleWake{<:Any,NK,TF},
        conversion::SurfaceVorticityConversion, system) where {NK,TF}

    NK == 1 || throw(ArgumentError(
        "SurfaceVorticityConversion supports scalar-strength wake kernels only " *
        "(got a kernel of strength dimension $(NK))"))

    pw = wake.panel_wake
    ws = wake.conversion_workspace
    diag = wake.conversion_diagnostics
    N = pw.nwakes[]

    capacity = size(pw.nodes[1], 2) - 1
    N == capacity || throw(WakeConversionStateError(
        "conversion is defined only on a full row buffer: nwakes[] = $N but " *
        "capacity is $capacity"))

    sigma = TF(conversion.sigma)
    overlap = TF(conversion.overlap)
    rank_rtol = TF(conversion.rank_rtol)
    geometry_rtol = TF(conversion.geometry_rtol)
    h = sigma / overlap

    # Streamwise attribution. `alpha` is the fraction of the outgoing row's
    # *upstream* face this conversion deposits; the complement stays on the
    # panel side, carried by the retained final filament, until the next
    # conversion deposits it as its own downstream face.
    alpha = TF(_ATTRIBUTION_ALPHA[conversion.attribution])
    handoff_before = pw.particle_handoff_active[]
    # Before the first handoff the aft-most face is the sheet's physical
    # trailing boundary -- the starting vortex -- not an interface with
    # particles, so there is no earlier conversion to have taken the `alpha`
    # share and it must be deposited whole. Omitting this deletes the starting
    # vortex outright (a one-time circulation loss).
    beta = handoff_before ? 1 - alpha : one(TF)
    startup_edge_deposited = !handoff_before

    _reset_workspace!(ws)
    panel_records = SurfaceVorticityPanelRecord{TF}[]

    total_area = zero(TF)
    area_strength = zero(SVector{3,TF})
    line_strength = zero(SVector{3,TF})
    expected_handoff = zero(SVector{3,TF})
    expected_total = zero(SVector{3,TF})
    root_scalar = zero(TF)
    tip_scalar = zero(TF)
    root_owned = false
    tip_owned = false
    n_elided = 0
    # L1 accumulation of every signed contribution to the ledger. Using this as
    # the relative-residual denominator keeps the measure meaningful when the
    # ledger itself sums to zero (a wrapping chain, or a constant strength
    # field), where the norm of the total would be pure round-off.
    ledger_scale = zero(TF)

    #--- Steps 1-4: reconstruct, quadrature, materialize, diagnose. No mutation.

    for i_surf in eachindex(pw.nodes)
        nodes = pw.nodes[i_surf]
        strength = pw.strength[i_surf]
        n_cols = size(nodes, 3) - 1

        # Wrap detection uses the legacy test verbatim, so the two strategies
        # classify topology identically (needed by the Stage 4b A/B).
        r1_le = _wake_node(nodes, N, 1)
        rend_le = _wake_node(nodes, N, n_cols + 1)
        wraps = norm(r1_le - rend_le) < 5 * eps()

        # Root/tip closure: the legacy filaments, verbatim. In divergence form
        # these are exactly the two streamwise faces that have no neighbour to
        # share with, so the stored cell value is the right strength and no
        # correction is needed.
        root_gamma = wraps ? zero(TF) : TF(strength[1, N, 1])
        tip_gamma = wraps ? zero(TF) : -TF(strength[1, N, n_cols])

        for j in 1:n_cols
            # Ghost panel in package contiguous order: xi runs streamwise
            # (upstream -> downstream), eta runs spanwise (increasing column).
            v1 = _wake_node(nodes, N, j)
            v2 = _wake_node(nodes, N + 1, j)
            v3 = _wake_node(nodes, N + 1, j + 1)
            v4 = _wake_node(nodes, N, j + 1)
            muhat_G = TF(strength[1, N, j])

            ref = _validate_wake_panel(v1, v2, v3, v4, geometry_rtol)
            n_centroid, _ = _bilinear_normal(v1, v2, v3, v4, TF(0.5), TF(0.5),
                                             geometry_rtol, ref)
            centroid_G = 0.25 * (v1 + v2 + v3 + v4)

            # Upstream strength. Only the *strength* is needed -- never the
            # upstream panel's geometry -- which is what lets `nwakerows == 1`
            # work without inventing anything: there the upstream row is the one
            # `shed_wake!` is about to create, whose nodes do not exist until the
            # next `update_TE!`, but whose strength is the body's shed value.
            if N >= 2
                muhat_A = TF(strength[1, N - 1, j])
            else
                system === nothing && throw(WakeConversionStateError(
                    "a single-row wake needs the shedding body to supply the " *
                    "upstream strength; call _convert_to_particles!(wake, system)"))
                strengthi, strengthj = _get_wakestrength_mu(system, j, i_surf)
                # Same expression (and sign) shed_wake! writes into row 1.
                muhat_A = TF(strengthi - strengthj)
            end

            # Spanwise neighbours. A wrapping chain has no boundary column; an
            # open chain's outer faces are physical and stay as line particles,
            # contributing nothing to the area.
            if wraps
                jm = j == 1 ? n_cols : j - 1
                jp = j == n_cols ? 1 : j + 1
            else
                jm = j - 1
                jp = j + 1
            end
            has_left = wraps || j > 1
            has_right = wraps || j < n_cols

            #--- Divergence form: the panel carries the filament content of the
            #    faces assigned to it, smeared over its own area.
            #
            #    upstream face  -> whole (Phase 1 eq. 8); the upstream panel
            #                      survives this transaction, so its share
            #                      cannot be deposited inside it
            #    downstream face-> nothing; already cancelled by the retained
            #                      final filament once the handoff is active
            #    spanwise faces -> half each, because both neighbours convert in
            #                      the same transaction (the split moves *where*
            #                      the vorticity sits, never how much)
            handoff = (muhat_A - muhat_G) * (v4 - v1)
            # Downstream face, in the legacy unsteady filament's orientation and
            # sign: v2 -> v3 carrying strength[1,N,j] - strength[1,N+1,j].
            downstream = (muhat_G - TF(strength[1, N + 1, j])) * (v3 - v2)
            face_left = has_left ?
                TF(0.5) * (muhat_G - TF(strength[1, N, jm])) * (v2 - v1) :
                zero(SVector{3,TF})
            face_right = has_right ?
                TF(0.5) * (TF(strength[1, N, jp]) - muhat_G) * (v3 - v4) :
                zero(SVector{3,TF})
            V = alpha * handoff + beta * downstream + face_left + face_right

            expected_handoff += handoff
            expected_total += V
            ledger_scale += LA.norm(alpha * handoff) + LA.norm(beta * downstream) +
                            LA.norm(face_left) + LA.norm(face_right)

            # Diagnostic only: what the Stage 2 centroid-stencil reconstruction
            # would have produced. No geometry is invented for it -- at N == 1
            # the streamwise leg is simply reported unobservable.
            d1 = N >= 2 ? _wake_panel_centroid(nodes, N - 1, j) - centroid_G :
                 zero(SVector{3,TF})
            dmu1 = N >= 2 ? muhat_A - muhat_G : zero(TF)
            if n_cols == 1
                d2 = zero(SVector{3,TF})
                dmu2 = zero(TF)
            else
                # Centered in the interior and on a wrapping chain; one-sided at
                # a true root/tip, where `jm`/`jp` above run off the sheet.
                dm = clamp(jm, 1, n_cols)
                dp = clamp(jp, 1, n_cols)
                d2 = _wake_panel_centroid(nodes, N, dp) - _wake_panel_centroid(nodes, N, dm)
                dmu2 = TF(strength[1, N, dp]) - TF(strength[1, N, dm])
            end
            res = (LA.norm(d1) > 0 || LA.norm(d2) > 0) ?
                _reconstruct_surface_gradient(n_centroid, d1, dmu1, d2, dmu2,
                                              rank_rtol, geometry_rtol, ref) :
                _null_surface_gradient(TF)
            kappa_recon = _surface_vorticity(n_centroid, res.gradient)

            n_xi, n_eta = _subdivision_counts(v1, v2, v3, v4, h)
            class = (!wraps && (j == 1 || j == n_cols)) ?
                RootTipSurfaceParticle : InteriorSurfaceParticle

            # Pass 1: subcell areas. kappa needs the summed panel area, and
            # using that sum (not an independent quadrature) is what makes the
            # deposited total exactly V.
            empty!(ws.areas)
            panel_area = zero(TF)
            min_jac = TF(Inf)
            for a in 1:n_xi, b in 1:n_eta
                xi0 = TF(a - 1) / n_xi; xi1 = TF(a) / n_xi
                eta0 = TF(b - 1) / n_eta; eta1 = TF(b) / n_eta
                _, jac = _bilinear_normal(v1, v2, v3, v4, (xi0 + xi1) / 2,
                                          (eta0 + eta1) / 2, geometry_rtol, ref)
                min_jac = min(min_jac, jac)
                dA = _subcell_area(v1, v2, v3, v4, xi0, xi1, eta0, eta1,
                                   geometry_rtol, ref)
                push!(ws.areas, dA)
                panel_area += dA
            end

            kappa = panel_area > 0 ? V / panel_area : zero(SVector{3,TF})

            # Pass 2: deposit. kappa is uniform over the panel, so no projection
            # onto the local normal -- projecting would break the sum.
            panel_strength = zero(SVector{3,TF})
            n_requested = 0
            n_panel_elided = 0
            k = 0
            for a in 1:n_xi, b in 1:n_eta
                k += 1
                dA = ws.areas[k]
                xi = (TF(a) - TF(0.5)) / n_xi
                eta = (TF(b) - TF(0.5)) / n_eta
                Xp = _bilinear_position(v1, v2, v3, v4, xi, eta)
                Gamma_p = kappa * dA

                n_requested += 1
                if iszero(Gamma_p)
                    n_panel_elided += 1
                    continue
                end
                panel_strength += Gamma_p
                _stage_particle!(ws, Xp, Gamma_p, sigma,
                                 LA.norm(Gamma_p) / sqrt(dA), class)
            end

            total_area += panel_area
            area_strength += panel_strength
            n_elided += n_panel_elided

            push!(panel_records, SurfaceVorticityPanelRecord{TF}(
                i_surf, j, res.rank, res.singular_values, res.rank_threshold,
                res.condition, res.observable_directions, res.observable,
                n_xi, n_eta, panel_area, panel_strength,
                kappa, kappa_recon, LA.norm(kappa - kappa_recon), handoff,
                downstream, class,
                n_requested, n_panel_elided, res.geometry_scale, min_jac))
        end

        # True root/tip closure: the legacy root and tip streamwise filaments,
        # unchanged in position, orientation, signed strength, and sampling.
        # Interior spanwise jumps, the wrap seam, and the downstream handoff
        # edge are all absent.
        if !wraps
            root_owned = true
            tip_owned = true
            r1 = _wake_node(nodes, N, 1)
            r2 = _wake_node(nodes, N + 1, 1)
            _, ne, gs = _stage_line_particles!(ws, r1, r2, root_gamma, sigma,
                                               overlap, PerimeterLineParticle)
            n_elided += ne
            line_strength += gs
            expected_total += root_gamma * (r2 - r1)
            ledger_scale += LA.norm(root_gamma * (r2 - r1))
            root_scalar += root_gamma

            r1 = _wake_node(nodes, N, n_cols + 1)
            r2 = _wake_node(nodes, N + 1, n_cols + 1)
            _, ne, gs = _stage_line_particles!(ws, r1, r2, tip_gamma, sigma,
                                               overlap, PerimeterLineParticle)
            n_elided += ne
            line_strength += gs
            expected_total += tip_gamma * (r2 - r1)
            ledger_scale += LA.norm(tip_gamma * (r2 - r1))
            tip_scalar += tip_gamma
        end
    end

    #--- Step 5: exact capacity check. Still no mutation.

    requested = _n_pending(ws)
    available = wake.pfield.maxparticles - wake.pfield.np
    requested <= available || throw(PanelParticleCapacityError(requested, available))
    nearfield = conversion.diagnose_nearfield ?
        _surface_vorticity_nearfield(pw, ws) : nothing

    #--- Step 6: commit. Capacity is guaranteed above, so no insertion can fail
    #    partway and leave a half-appended row observable.

    for k in 1:requested
        FLOWVPM.add_particle(wake.pfield, ws.positions[k], ws.strengths[k],
                             ws.sigmas[k]; circulation=ws.circulations[k])
    end

    pw.particle_handoff_active[] = true
    pw.particle_handoff_weight[] = Float64(alpha)
    wake.conversion_count[] += 1

    deposited_total = area_strength + line_strength
    residual = deposited_total - expected_total
    denom = ledger_scale
    n_interior = count(==(InteriorSurfaceParticle), ws.classes)
    n_roottip = count(==(RootTipSurfaceParticle), ws.classes)
    n_perimeter = count(==(PerimeterLineParticle), ws.classes)

    _accumulate!(diag, SurfaceVorticityConversionRecord{TF}(
        wake.conversion_count[], pw.live_step_id[], panel_records,
        n_interior, n_roottip, n_perimeter, total_area, area_strength,
        root_scalar, tip_scalar, root_owned, tip_owned,
        handoff_before, true, conversion.attribution, startup_edge_deposited,
        expected_handoff, expected_total, deposited_total,
        LA.norm(residual), denom > 0 ? LA.norm(residual) / denom : zero(TF),
        n_elided, requested, available, nearfield))

    return nothing
end



#------- vortex filament (methods for FilamentWrapper, struct defined above) -------#

FastMultipole.numtype(fw::FilamentWrapper) = FastMultipole.numtype(fw.system)

function fmm_to_filament_index(filaments::FilamentWrapper{<:PanelWake}, n)
    n_counter = 0
    system = filaments.system
    for k in eachindex(system.nodes)
        ns = size(system.strength[k], 3)
        if n_counter + ns >= n # found the surface
            i_filament = n - n_counter
            i_surf = k
            return i_surf, i_filament
        end
        n_counter += ns
    end
end

"""
Strength of the retained filament on the downstream edge of the final active
panel row.

Legacy bookkeeping (`unsteady_filament=true`) reaches one row *past* the final
active row, so the filament carries the not-yet-converted unsteady jump.

Once a surface-vorticity handoff has committed (BRAINSTORM 016), the filament
must instead leave on the panel side exactly the fraction of the interface that
the conversion did *not* deposit. With `alpha = particle_handoff_weight[]` the
required strength is

```
-(alpha * strength[i_row] + (1 - alpha) * strength[i_row+1])
```

which reproduces both legacy forms exactly at its endpoints: `alpha = 1` is the
fully cancelled edge (`unsteady_filament=false` sign of Phase 1 eq. 8) and
`alpha = 0` is the legacy unsteady filament.
"""
function _final_filament_strength(wake::PanelWake, i_surf, i_row, j)
    strength = wake.strength[i_surf]
    if wake.particle_handoff_active[]
        alpha = wake.particle_handoff_weight[]
        # alpha == 0 (:downstream attribution) never reads strength[1, i_row, j].
        # Short-circuiting it is bit-identical for i_row >= 1 (0.0 * finite +
        # 1.0 * b == b) and is the only reachable handoff case at i_row == 0,
        # the convert-at-shed steady state (BRAINSTORM 024), where row `i_row`
        # does not exist and the filament on the TE+Das line must cancel the
        # full just-shed row-1 strength.
        alpha == 0 && return -strength[1, i_row + 1, j]
        return -(alpha * strength[1, i_row, j] +
                 (1 - alpha) * strength[1, i_row + 1, j])
    end
    strength_row = wake.unsteady_filament ? i_row + 1 : i_row
    return -strength[1, strength_row, j]
end

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, filaments::FilamentWrapper{<:PanelWake}, i_body)
    
    # vlm index
    wake = filaments.system
    i_surf, j = fmm_to_filament_index(filaments, i_body)

    # which row of panels we're on
    i_row = wake.nwakes[]

    # get strength
    strength = _final_filament_strength(wake, i_surf, i_row, j)

    # nodes
    v1 = SVector{3}(view(wake.nodes[i_surf], :, i_row+1, j))
    v2 = SVector{3}(view(wake.nodes[i_surf], :, i_row+1, j+1))

    # update buffer
    buffer[1:3, i_buffer] .= 0.5 * (v1 + v2)
    # regularization reach, not just the core itself: at clearance = core_size
    # the Vatistas kernel still differs from the singular kernel by ~30%
    buffer[4, i_buffer] = 0.5 * norm(v2 - v1) +
        radius_inflation(VortexRing, wake.core_size, fmm_radius_tolerance(wake))
    buffer[5, i_buffer] = strength
    buffer[6:8,i_buffer] .= v1
    buffer[9:11,i_buffer] .= v2
    buffer[12,i_buffer] = wake.core_size
end

function FastMultipole.data_per_body(wakes::FilamentWrapper{<:PanelWake})
    return 12
end

function FastMultipole.get_position(filaments::FilamentWrapper{<:PanelWake}, i)
    wake = filaments.system
    i_surf, j = fmm_to_filament_index(filaments, i)
    i_row = wake.nwakes[]

    # nodes
    v1 = SVector{3}(view(wake.nodes[i_surf], :, i_row+1, j))
    v2 = SVector{3}(view(wake.nodes[i_surf], :, i_row+1, j+1))

    return 0.5 * (v1 + v2)
end

function FastMultipole.strength_dims(filaments::FilamentWrapper{<:PanelWake})
    return 1
end

FastMultipole.has_vector_potential(filaments::FilamentWrapper{<:PanelWake}) = true

function FastMultipole.get_n_bodies(filaments::FilamentWrapper{<:PanelWake})
    nwakes = 0
    if filaments.system.overflowed[]
        wake = filaments.system
        for str in wake.strength
            ns = size(str, 3)
            nwakes += ns
        end
    end
    return nwakes
end

function FastMultipole.body_to_multipole!(filaments::FilamentWrapper{<:PanelWake}, multipole_coefficients, buffer::Matrix, center, bodies_index, harmonics, expansion_order)
    # loop over bodies
    for i_body in bodies_index
       
        # extract vertices from buffer
        rtl = FastMultipole.get_vertex(buffer, filaments, i_body, 1)
        rtr = FastMultipole.get_vertex(buffer, filaments, i_body, 2)

        # extract strength from buffer
        gamma = FastMultipole.get_strength(buffer, filaments, i_body)[1]

        # top bound vortex
        body_to_multipole_vl!(multipole_coefficients, harmonics, rtl, rtr, center, gamma, expansion_order)
    end
end

# function barrier: family in the type domain inside the loop (BRAINSTORM 025)
function FastMultipole.direct!(target_system, target_index, switch::FastMultipole.DerivativesSwitch, source_system::FilamentWrapper, source_buffer, source_index)
    _direct_filaments!(target_system, target_index, switch, source_system,
        source_buffer, source_index, Val(FILAMENT_REGULARIZATION[]))
end

function _direct_filaments!(target_system, target_index, switch::FastMultipole.DerivativesSwitch{PS,VS,GS,NO,NM}, source_system::FilamentWrapper, source_buffer, source_index, fam::Val) where {PS,VS,GS,NO,NM}
    TF = FastMultipole.numtype(source_system)
    @inbounds for j_target in target_index
        target = FastMultipole.get_position(target_system, j_target)
        v = SVector{3,TF}(0.0, 0.0, 0.0)
        g = zero(SMatrix{3,3,TF,9})
        w = SVector{3,TF}(0.0, 0.0, 0.0)
        @inbounds for i_source in source_index
            v1 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 1)
            v2 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 2)
            gamma = FastMultipole.get_strength(source_buffer, source_system, i_source)[1]
            cs = source_buffer[12, i_source]
            if VS
                v += _bound_vortex_velocity(target-v1, target-v2, true, cs, fam) * gamma
            end
            if GS
                g += _bound_vortex_gradient(v1-target, v2-target, true, cs, fam) * gamma
            end
            if NO == 3
                gf = _bound_vortex_gradient(v1-target, v2-target, true, cs, fam) * gamma
                w += SVector{3,TF}(gf[3, 2] - gf[2, 3],
                                   gf[1, 3] - gf[3, 1],
                                   gf[2, 1] - gf[1, 2])
            end
        end
        VS && FastMultipole.set_gradient!(target_system, switch, j_target, v)
        GS && FastMultipole.set_hessian!(target_system, switch, j_target, g)
        if NO == 3
            FastMultipole.set_extra_output!(target_system, switch, j_target, 1, w[1])
            FastMultipole.set_extra_output!(target_system, switch, j_target, 2, w[2])
            FastMultipole.set_extra_output!(target_system, switch, j_target, 3, w[3])
        end
    end
end

function FastMultipole.buffer_to_target_system!(target_system::FilamentWrapper, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS,NO,NM}, target_buffer, i_buffer) where {PS,VS,GS,NO,NM}
    @warn "A `::FilamentWrapper` object should not be used as a target in an FMM call."
end
