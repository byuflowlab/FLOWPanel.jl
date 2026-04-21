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
    influence!((body, wake_probes), (wake,), backend; gradient=true, hessian=(requires_hessian(body), requires_hessian(wake)))

    # freestream
    apply_freestream!(body, uinf)
    apply_freestream!(wake, uinf)

    # kinematics
    kinematic_velocity!(body, frames; skip_top_level=false)

    # solve body (shouldn't modify body velocity, but will update body strength)
    solve!(body, body_solver; backend)

    # body-on-all influence
    influence!((body, wake_probes), (body,), backend; gradient=true, hessian=(requires_hessian(body), requires_hessian(wake)))

    return nothing
end

requires_hessian(::AbstractBody) = false # default behavior
requires_hessian(::AbstractFreeWake) = false # default behavior
requires_hessian(pw::ProbeWrapper) = requires_hessian(pw.system)

#--- Panel Wake ---#

"""
    PanelWake(shedding, kernel, TF=Float64; core_size=1e-3, nwakerows=100)
    PanelWake(body; kernel=get_wake_kernel(body), nwakerows=100)

Wake model that stores a panelized wake sheet behind one or more shedding-edge
chains.
"""
struct PanelWake{TK,NK,TF} <: AbstractFreeWake
    nwakes::Array{Int, 0}
    nodes::Vector{Array{TF, 3}}
    strength::Vector{Array{TF, 3}}
    velocity::Vector{Array{TF, 3}}
    core_size::Float64
    overflowed::Array{Bool, 0}
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
    return (wake, FilamentWrapper(wake))
end

function PanelWake(shedding::Vector{Matrix{Int}}, kernel, TF=Float64; 
        core_size=1e-3, nwakerows=100
    )
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
    
    # overflowed
    overflowed = Array{Bool,0}(undef)
    overflowed[] = false

    return PanelWake{kernel, dim, TF}(nwakes, nodes, strength, velocity, core_size, overflowed)
end

PanelWake(body::AbstractLiftingBody{TK,NK,TF}, kernel=get_wake_kernel(body); nwakerows=100) where {TK,NK,TF} =
    PanelWake(body.shedding, kernel, TF; nwakerows)

function reset!(wake::PanelWake)
    for vel in wake.velocity
        vel .= zero(eltype(vel))
    end
end

function apply_freestream!(wake::PanelWake, uinf)
    for vel in wake.velocity
        for ns in axes(vel, 3)
            for nc in 1:wake.nwakes[]+1
                vel[:, nc, ns] .+= uinf
            end
        end
    end
end

# FastMultipole compatibility
function global_to_matrix_index(wake::PanelWake, i_wake)

    # determine which shedding surface we're on
    nwakes = wake.nwakes[]
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
    # convert matrix indices to local index
    i_wake = (icol - 1) * wake.nwakes[] + irow

    # account for previous surfaces
    for i in 1:(isurf-1)
        i_wake += size(wake.strength[i], 3) * wake.nwakes[]
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

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, system::PanelWake, i_body)

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
    buffer[4, i_buffer] = max(r1, r2)

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

FastMultipole.get_n_bodies(system::PanelWake) = system.nwakes[] * sum(size(s, 3) for s in system.strength)

FastMultipole.get_n_bodies(system::ProbeWrapper{<:PanelWake}) = (system.system.nwakes[]+1) * sum(size(s, 3) for s in system.system.nodes)

function FastMultipole.buffer_to_target_system!(target_system::ProbeWrapper{<:PanelWake}, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    # get surface index of global `i_target` index
    isurf, irow, icol = global_to_matrix_index(target_system, i_target)

    # save potential
    # if PS
    #     target_system.potential[isurf][irow, icol] = target_buffer[i_buffer]
    # end

    # save velocity
    if VS
        target_system.system.velocity[isurf][1, irow, icol] += target_buffer[5, i_buffer]
        target_system.system.velocity[isurf][2, irow, icol] += target_buffer[6, i_buffer]
        target_system.system.velocity[isurf][3, irow, icol] += target_buffer[7, i_buffer]
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

function induced(target::AbstractVector{TF}, source_system::PanelWake{TK,NK,<:Any}, source_buffer::Matrix, i_source, derivatives_switch=FastMultipole.DerivativesSwitch(false,true,false)) where {TF,TK,NK}

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
    kerneloffset = source_system.core_size
    potential, velocity, velocity_gradient = _induced(target, (v1, v2, v3), control_point, strength, TK, kerneloffset, R, derivatives_switch)

    #--- second triangle ---#

    R = rotate_to_panel(v1, v3, v4)

    # get control point and strength
    control_point = (v1 + v3 + v4) * 0.3333333333333333

    # evaluate influence
    potential2, velocity2, velocity_gradient2 = _induced(target, (v1, v3, v4), control_point, strength, TK, kerneloffset, R, derivatives_switch)

    # return potential, velocity, velocity_gradient
    return potential+potential2, velocity+velocity2, velocity_gradient+velocity_gradient2
end

function FastMultipole.direct!(target_system, target_index, derivatives_switch::FastMultipole.DerivativesSwitch{PS,GS,HS}, source_system::PanelWake, source_buffer, source_index) where {PS,GS,HS}
    TF = eltype(target_system)
    for i_target in target_index # loop over targets
        target = FastMultipole.StaticArrays.SVector{3,TF}(target_system[1, i_target],
                  target_system[2, i_target],
                  target_system[3, i_target])
        
        phi_out = zero(eltype(target_system))
        U_out = @SVector zeros(eltype(target_system), 3)

        for i_source in source_index # loop over sources
            # evaluate influence due to this source
            phi, U, H = induced(target, source_system, source_buffer, i_source, derivatives_switch)
            phi_out += phi
            U_out += U
        end

        # store results
        if PS
            target_system[4, i_target] += phi_out
        end
        if GS
            target_system[5, i_target] += U_out[1]
            target_system[6, i_target] += U_out[2]
            target_system[7, i_target] += U_out[3]
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

function propagate!(wake::PanelWake, dt; step=0)
    for i_surf in eachindex(wake.nodes)
        view(wake.velocity[i_surf], :, 1:wake.nwakes[]+1, :) .*= dt # displacements
        view(wake.nodes[i_surf], :, 1:wake.nwakes[]+1, :) .+= view(wake.velocity[i_surf], :, 1:wake.nwakes[]+1, :) # update nodes
        view(wake.velocity[i_surf], :, 1:wake.nwakes[]+1, :) ./= dt # restore velocities
    end
end

# Forward declaration (methods defined later in the vortex filament section)
struct FilamentWrapper{TS}
    system::TS
end

function write_vtk(name, wake::PanelWake, idx, t; overwrite=false)
    # Route block files to a subdirectory named after the PVD
    _parent, _base = splitdir(name)
    subdir = joinpath(_parent, _base)
    mkpath(subdir)
    block_name = joinpath(subdir, _base)

    WriteVTK.paraview_collection(name; append=!overwrite) do pvd
        vtm = WriteVTK.vtk_multiblock(block_name * ".$idx.vtm")
        if wake.nwakes[] > 0
            for i_surf in eachindex(wake.nodes)
                pts = view(wake.nodes[i_surf], :, 1:wake.nwakes[]+1, :)
                pts_reshaped = reshape(pts, 3, wake.nwakes[]+1, size(wake.nodes[i_surf], 3), 1)
                WriteVTK.vtk_grid(vtm, block_name * ".$(i_surf).$(idx).vts", pts_reshaped) do vtk
                    vel = view(wake.velocity[i_surf], :, 1:wake.nwakes[]+1, :)
                    vtk["velocity", WriteVTK.VTKPointData()] = reshape(vel, 3, wake.nwakes[]+1, size(wake.nodes[i_surf], 3), 1)

                    str = view(wake.strength[i_surf], :, 1:wake.nwakes[], :)
                    vtk["strength", WriteVTK.VTKCellData()] = reshape(str, size(str, 1), wake.nwakes[], size(wake.nodes[i_surf], 3)-1, 1)
                end
            end
        end
        pvd[t] = vtm
    end

    # filaments at trailing edge of last panel row
    write_vtk(name * "_filaments", FilamentWrapper(wake), idx, t; overwrite)
end

function write_vtk(name, filaments::FilamentWrapper{<:PanelWake}, idx, t; overwrite=false)
    wake = filaments.system
    if !wake.overflowed[]
        # no filaments yet — still write an empty PVD entry
        _parent, _base = splitdir(name)
        subdir = joinpath(_parent, _base)
        mkpath(subdir)
        block_name = joinpath(subdir, _base)
        WriteVTK.paraview_collection(name; append=!overwrite) do pvd
            vtm = WriteVTK.vtk_multiblock(block_name * ".$idx.vtm")
            pvd[t] = vtm
        end
        return
    end

    i_row = wake.nwakes[]

    _parent, _base = splitdir(name)
    subdir = joinpath(_parent, _base)
    mkpath(subdir)
    block_name = joinpath(subdir, _base)

    WriteVTK.paraview_collection(name; append=!overwrite) do pvd
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
                strengths[j] = -wake.strength[i_surf][1, i_row + 1, j]
            end

            WriteVTK.vtk_grid(vtm, block_name * ".$(i_surf).$(idx).vtu", points, cells) do vtk
                vtk["strength", WriteVTK.VTKCellData()] = strengths
            end
        end

        pvd[t] = vtm
    end
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

function _shed_particles!(pfield, r1, r2, Γ, method::SigmaPPS)
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

#------- Vortex Particle Wake -------#

"""
    PanelParticleWake(body; optargs...)

Hybrid wake model that combines a near-body panel wake with vortex particles
shed downstream from the trailing edge.
"""
struct PanelParticleWake{TK,NK,TF,TPF,MT,MU} <: AbstractFreeWake
    panel_wake::PanelWake{TK,NK,TF}
    pfield::TPF                           # FLOWVPM.ParticleField object
    method_trailing::MT                             # particle shedding method
    method_unsteady::MU                             # particle shedding method
    gamma_trim_threshold::TF              # remove particles with |Gamma| below this value
    # Particle merge settings
    merge_every::Int                      # merge every N steps (0 = disabled)
    merge_r::Float64                      # r_merge kwarg
    merge_sigma_relative::Bool            # sigma_relative kwarg
    merge_max_sigma_ratio::Float64        # max_sigma_ratio kwarg
    merge_skip_static::Bool               # skip_static kwarg
    check_neighboring_cells::Bool         # check_neighboring_cells kwarg
    merge_verbose::Bool
end

function PanelParticleWake(body::AbstractLiftingBody;
        nwakerows=3, max_particles=10000,
        method_trailing::WakeSheddingMethod=OverlapPPS(1.3, 2),
        method_unsteady::WakeSheddingMethod=OverlapPPS(1.3, 2),
        gamma_trim_threshold=1e-8,
        merge_every=0,
        merge_r=0.5,
        merge_sigma_relative=true,
        merge_max_sigma_ratio=2.0,
        merge_skip_static=true,
        check_neighboring_cells=true,
        merge_verbose=false,
        kwargs...)

    panel_wake = PanelWake(body; nwakerows, kwargs...)
    TF = FastMultipole.numtype(panel_wake)

    # Create particle field with default settings (disable autotune_reg_error to avoid convergence issues)
    pfield = FLOWVPM.ParticleField(max_particles, TF;
        fmm=FLOWVPM.FMM(autotune_reg_error=false))

    # Infer type params from the actual panel_wake
    WTK = typeof(panel_wake).parameters[1]
    WNK = typeof(panel_wake).parameters[2]
    gamma_trim_threshold_TF = convert(TF, gamma_trim_threshold)
    return PanelParticleWake{WTK,WNK,TF,typeof(pfield),typeof(method_trailing),typeof(method_unsteady)}(
        panel_wake, pfield, method_trailing, method_unsteady, gamma_trim_threshold_TF,
        merge_every, merge_r, merge_sigma_relative, merge_max_sigma_ratio, merge_skip_static, check_neighboring_cells, merge_verbose
    )
end

"""
Run SFS pre-calculations for particle field before evaluating the velocity field.
"""
function pre_evaluate_influence!(pfield::FLOWVPM.ParticleField)
    pfield.SFS(pfield, FLOWVPM.BeforeUJ())
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

function apply_freestream!(w::PanelParticleWake, uinf)
    # apply to panel wake
    apply_freestream!(w.panel_wake, uinf)

    # Add freestream to particle velocities
    for i in 1:w.pfield.np
        for d in 1:3
            w.pfield.particles[FLOWVPM.U_INDEX[d], i] += uinf[d]
        end
    end
end

update_TE!(w::PanelParticleWake, sys) = update_TE!(w.panel_wake, sys)

function trim_weak_particles!(w::PanelParticleWake)
    threshold = w.gamma_trim_threshold
    threshold <= zero(threshold) && return nothing

    for i in w.pfield.np:-1:1
        gamma = FLOWVPM.get_Gamma(w.pfield, i)
        gamma_mag = sqrt(sum(abs2, gamma))
        if gamma_mag < threshold
            FLOWVPM.remove_particle(w.pfield, i)
        end
    end

    return nothing
end

function propagate!(w::PanelParticleWake, dt; relax=true, step=0)

    # panel wake
    propagate!(w.panel_wake, dt)

    # convect particles
    FLOWVPM._euler(w.pfield, dt; relax)

    # particle merging
    if w.merge_every > 0 && step > 0 && step % w.merge_every == 0
        FLOWVPM.merge_particles!(w.pfield;
            r_merge=w.merge_r,
            sigma_relative=w.merge_sigma_relative,
            max_sigma_ratio=w.merge_max_sigma_ratio,
            skip_static=w.merge_skip_static,
            check_neighboring_cells=w.check_neighboring_cells,
        )
    end

    trim_weak_particles!(w)
end

function write_vtk(name, w::PanelParticleWake, idx, t; overwrite=false)
    # panel wake (includes filaments)
    write_vtk(name, w.panel_wake, idx, t; overwrite)

    # particle wake — route block files to subdirectory
    vpm_path, vpm_name = splitdir(name)
    particles_pvd_name = joinpath(vpm_path, vpm_name * "_particles")
    particles_subdir = joinpath(vpm_path, vpm_name * "_particles")
    mkpath(particles_subdir)
    particles_block = joinpath(particles_subdir, vpm_name * "_particles")

    np = w.pfield.np
    X = view(w.pfield.particles, FLOWVPM.X_INDEX, 1:np)
    cells = [WriteVTK.MeshCell(WriteVTK.PolyData.Verts(), 1:np)]

    vtp_filename = particles_block * ".$idx.vtp"
    vtp = WriteVTK.vtk_grid(vtp_filename, X, cells)

    if np > 0
        vtp["gamma", WriteVTK.VTKPointData()] = view(w.pfield.particles, FLOWVPM.GAMMA_INDEX, 1:np)
        vtp["sigma", WriteVTK.VTKPointData()] = view(w.pfield.particles, FLOWVPM.SIGMA_INDEX, 1:np)
        vtp["vol", WriteVTK.VTKPointData()] = view(w.pfield.particles, FLOWVPM.VOL_INDEX, 1:np)
        vtp["circulation", WriteVTK.VTKPointData()] = view(w.pfield.particles, FLOWVPM.CIRCULATION_INDEX, 1:np)
        vtp["velocity", WriteVTK.VTKPointData()] = view(w.pfield.particles, FLOWVPM.U_INDEX, 1:np)
        vtp["vorticity", WriteVTK.VTKPointData()] = view(w.pfield.particles, FLOWVPM.VORTICITY_INDEX, 1:np)
        vtp["velocity_gradient", WriteVTK.VTKPointData()] = reshape(view(w.pfield.particles, FLOWVPM.J_INDEX, 1:np), 3, 3, np)
    end

    pvd = WriteVTK.paraview_collection(particles_pvd_name; append=!overwrite)
    pvd[t] = vtp

    WriteVTK.vtk_save(pvd)
end

requires_hessian(::FLOWVPM.ParticleField) = true

# --- FastMultipole compatibility ---

# Ensure ParticleField has numtype for FastMultipole integration
FastMultipole.numtype(pf::FLOWVPM.ParticleField) = eltype(pf)

# --- Save and convert last row to particles ---

function _convert_to_particles!(wake::PanelParticleWake)

    nwakes = wake.panel_wake.nwakes[]
    
    for i_surf in eachindex(wake.panel_wake.nodes)
        nodes = wake.panel_wake.nodes[i_surf]
        strength = wake.panel_wake.strength[i_surf]
        n_cols = size(nodes, 3) - 1

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
            _shed_particles!(wake.pfield, r1_le, r1_te, Γ-Γ_last, wake.method_trailing)

            # unsteady particle
            r2_te = SVector{3}(nodes[1, nwakes+1, icol+1], nodes[2, nwakes+1, icol+1], nodes[3, nwakes+1, icol+1])
            Γ_tm1 = strength[1,nwakes+1,icol]
            _shed_particles!(wake.pfield, r1_te, r2_te, Γ-Γ_tm1, wake.method_unsteady)

            Γ_last = Γ
        end

        if !wraps
            # right trailing particles (streamwise direction)
            Γ = strength[1,nwakes,n_cols]
            r1_le = SVector{3}(nodes[1, nwakes, n_cols+1], nodes[2, nwakes, n_cols+1], nodes[3, nwakes, n_cols+1])
            r1_te = SVector{3}(nodes[1, nwakes+1, n_cols+1], nodes[2, nwakes+1, n_cols+1], nodes[3, nwakes+1, n_cols+1])
            _shed_particles!(wake.pfield, r1_le, r1_te, -Γ, wake.method_trailing)
        end
    end
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

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, filaments::FilamentWrapper{<:PanelWake}, i_body)
    
    # vlm index
    wake = filaments.system
    i_surf, j = fmm_to_filament_index(filaments, i_body)

    # which row of panels we're on
    i_row = wake.nwakes[]

    # get strength
    strength = -wake.strength[i_surf][1, i_row+1, j]
    # if i_row > 0
    #     strength += wake.strength[i_surf][1, i_row, j]
    # end

    # nodes
    v1 = SVector{3}(view(wake.nodes[i_surf], :, i_row+1, j))
    v2 = SVector{3}(view(wake.nodes[i_surf], :, i_row+1, j+1))

    # update buffer
    buffer[1:3, i_buffer] .= 0.5 * (v1 + v2)
    buffer[4, i_buffer] = 0.5 * norm(v2 - v1) + wake.core_size
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

function FastMultipole.direct!(target_system, target_index, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, source_system::FilamentWrapper, source_buffer, source_index) where {PS,VS,GS}
    TF = FastMultipole.numtype(source_system)
    @inbounds for j_target in target_index
        target = FastMultipole.get_position(target_system, j_target)
        v = SVector{3,TF}(0.0, 0.0, 0.0)
        g = zero(SMatrix{3,3,TF,9})
        @inbounds for i_source in source_index
            v1 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 1)
            v2 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 2)
            gamma = FastMultipole.get_strength(source_buffer, source_system, i_source)[1]
            cs = source_buffer[12, i_source]
            if VS
                v += _bound_vortex_velocity(target-v1, target-v2, true, cs) * gamma
            end
            if GS
                g += _bound_vortex_gradient(target-v1, target-v2, true, cs) * gamma
            end
        end
        VS && FastMultipole.set_gradient!(target_system, j_target, v)
        GS && FastMultipole.set_hessian!(target_system, j_target, g)
    end
end

function FastMultipole.buffer_to_target_system!(target_system::FilamentWrapper, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    @warn "A `::FilamentWrapper` object should not be used as a target in an FMM call."
end
