abstract type AbstractFreeWake end

"""
    solve!(body::AbstractBody{TB}, wake::AbstractFreeWake{TW}, Uinf::Function, t=0.0; solver::AbstractSolver, backend=FastMultipoleBackend)

Solve for the body panel strengths and wake velocities at time `t` given a freestream velocity function `Uinf(t)`. The `solver` keyword argument specifies the linear solver to use for solving the body panel strengths, and the `backend` keyword argument specifies the N-body backend to use for evaluating influence of body and wake on each other.

**Modified Arguments**
- `body`: the body panel strength and velocity fields are modified in-place
- `wake`: the wake velocity field is modified in-place

**Additional Arguments**
- `Uinf::Function`: a function that accepts time `t` and returns the freestream velocity vector at that time
- `t::Real`: the time at which to evaluate the freestream velocity and solve for the body and wake
- `solver::AbstractSolver`: the linear solver to use for solving the body panel strengths (default: `BackslashDirichlet(body)`)
- `backend::AbstractBackend`: the N-body backend to use for evaluating influence of body and wake on each other (default: `FastMultipoleBackend` with expansion order 10, multipole acceptance 0.4, and leaf size 100)

"""
function solve!(body::AbstractBody{TB}, wake::AbstractFreeWake, uinf::AbstractArray, t=0.0;
        body_solver::AbstractSolver=BackslashDirichlet(body), 
        backend=FastMultipoleBackend(;
            expansion_order=10,
            multipole_acceptance=0.4,
            leaf_size=100,
            shrink=true,
            recenter=false,
        )
    ) where {TB}
    
    # get probes
    body_probes = get_probes(body)

    # wake-on-all velocity
    evaluate_influence!((body_probes, wake), (wake,), backend; gradient=true, hessian=(requires_hessian(body), requires_hessian(wake)))

    # freestream
    eachcol(body_probes.gradient) .+= Ref(uinf)
    for vel in wake.velocity
        for ns in axes(vel, 3)
            for nc in axes(vel, 2)
                vel[:, nc, ns] .+= uinf
            end
        end
    end

    # solve body
    check_nans(body)
    @show true in body_probes.gradient
    @show true in body_probes.scalar_potential
    solve2!(body, body_probes.gradient, body_solver; backend)

    # body-on-all influence
    evaluate_influence!((body_probes, wake_probes), (body,), backend; gradient=true, hessian=(requires_hessian(body), requires_hessian(wake)))

    return nothing
end

function check_nans(solver::FGSSolver)
    if any(isnan.(solver.fgs.extra_farfield))
        error("NaN detected in farfield extra points")
    end
    if any(isnan.(solver.fgs.source_tree.buffers[1]))
        error("NaN detected in FMM source buffer")
    end
    if any(isnan.(solver.fgs.target_tree.buffers[1]))
        error("NaN detected in FMM target buffer")
    end
    if any(isnan.(solver.fgs.strengths))
        error("NaN detected in fgs strengths")
    end
    if any(isnan.(solver.fgs.old_influence_storage))
        error("NaN detected in fgs old influence storage")
    end
    if any(isnan.(solver.fgs.extra_right_hand_side))    
        error("NaN detected in fgs extra right hand side")
    end
end

function check_nans(body::AbstractBody)
    if any(isnan.(body.potential))
        error("NaN detected in body potential")
    end
    if any(isnan.(body.velocity))
        error("NaN detected in body velocity")
    end
    if true in isnan.(body.strength)
        error("NaN detected in body strength")
    end
    for das in body.Das
        if true in isnan.(das)
            error("NaN detected in Das")
        end
    end
    for dbs in body.Dbs
        if true in isnan.(dbs)
            error("NaN detected in Dbs")
        end
    end
    if true in isnan.(body.grid._nodes)
        error("NaN detected in body grid nodes")
    end

end

requires_hessian(::AbstractBody) = false # default behavior
requires_hessian(::AbstractFreeWake) = false # default behavior

function get_probes(body::AbstractBody{TK,NK,TF}) where {TK,NK,TF}
    hessian = Array{TF, 3}(undef, 0, 0, 0)  # unused
    normals = _calc_normals(body)
    CPs = _calc_controlpoints(body, normals; off=-1e-10)
    @assert !requires_hessian(body) "`get_probes` must be overloaded to support Hessian output for body type $(typeof(body))"
    body_probes = FastMultipole.ProbeSystemArray(CPs, body.potential, body.velocity, hessian)
    return body_probes
end

#--- Panel Wake ---#

struct PanelWake{TK,NK,TF} <: AbstractFreeWake
    nwakes::Array{Int, 0}
    nodes::Vector{Array{TF, 3}}
    strength::Vector{Array{TF, 3}}
    velocity::Vector{Array{TF, 3}}
    kerneloffset::Float64
end

function get_probes(wake::PanelWake)
    return wake
end

function PanelWake(shedding::Vector{Matrix{Int}}, kernel, TF=Float64; 
        kerneloffset=1e-10, nwakerows=100
    )
    # nwakes
    nwakes = Array{Int,0}(undef)
    nwakes[] = 0

    # nodes
    nodes = [zeros(TF, 3, nwakerows+1, size(s, 2)+1) for s in shedding]

    # strength
    dim = kernel_dim(kernel)
    strength = [zeros(TF, dim, nwakerows, size(s, 2)) for s in shedding]

    # velocity
    velocity = [zeros(TF, size(n)) for n in nodes]

    return PanelWake{kernel, dim, TF}(nwakes, nodes, strength, velocity, kerneloffset)
end

PanelWake(body::AbstractLiftingBody{TK,NK,TF}, kernel=ConstantDoublet; nwakerows=100) where {TK,NK,TF} = 
    PanelWake(body.shedding, kernel, TF; kerneloffset=body.kerneloffset, nwakerows)

# FastMultipole compatibility
function global_to_matrix_index(wake::PanelWake, i_wake)

    # determine which shedding surface we're on
    nwakes = wake.nwakes[]
    isurf = 1
    i_wake_local = i_wake
    npanels = 0
    for i in eachindex(wake.strength)
        npanels += size(wake.strength[i], 2) * nwakes
        if i_wake <= npanels
            break
        end
        isurf += 1 # advance to the next surface
        i_wake_local -= size(wake.strength[i], 2) * nwakes # adjust local index
    end

    # convert local index to matrix indices
    icol, irow = divrem(i_wake_local - 1, size(wake.strength[isurf], 2))
    icol += 1 # adjust for 1-based indexing
    irow += 1 # adjust for 1-based indexing

    return isurf, irow, icol
end

function matrix_to_global_index(wake::PanelWake, isurf, irow, icol)
    # convert matrix indices to local index
    i_wake = (icol - 1) * size(wake.strength[isurf], 2) + irow

    # account for previous surfaces
    for isurf in 1:(isurf-1)
        i_wake += size(wake.strength[isurf], 2) * wake.nwakes[]
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

FastMultipole.data_per_body(system::PanelWake) = 4 + size(system.strength[1], 1) + 12

FastMultipole.has_vector_potential(::PanelWake{TK,NK,TF}) where {TK,NK,TF} = TK<:Union{VortexRing, ConstantVortexSheet}

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

FastMultipole.strength_dims(system::PanelWake) = size(system.strength[1], 1)

FastMultipole.get_n_bodies(system::PanelWake) = system.nwakes[] * sum(size(s, 2) for s in system.strength)

function FastMultipole.buffer_to_target_system!(target_system::PanelWake, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    # get surface index of global `i_target` index
    isurf, irow, icol = global_to_matrix_index(target_system, i_target)

    # save potential
    # if PS
    #     target_system.potential[isurf][irow, icol] = target_buffer[i_buffer]
    # end

    # save velocity
    if VS
        target_system.velocity[isurf][1, irow, icol] = target_buffer[5, i_buffer]
        target_system.velocity[isurf][2, irow, icol] = target_buffer[6, i_buffer]
        target_system.velocity[isurf][3, irow, icol] = target_buffer[7, i_buffer]
    end

    # save Hessian (not currently used for PanelWake)
    if GS
        @warn "Hessian output not currently implemented for PanelWake targets"
    end
end

function rotate_to_panel(v1, v2, v3)
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
    kerneloffset = source_system.kerneloffset
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
            phi, U, _ = induced(target, source_system, source_buffer, i_source, derivatives_switch)
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
            @warn "Hessian output not currently implemented for PanelWake targets"
        end
    end
end

FastMultipole.body_to_multipole!(system::PanelWake{ConstantDoublet, 1, <:Any}, args...) =
    FastMultipole.body_to_multipole_quad!(FastMultipole.Panel{FastMultipole.Dipole}, system, args...)

FastMultipole.body_to_multipole!(system::PanelWake{VortexRing, 1, <:Any}, args...) = @error "Multipole expansion not currently implemented for VortexRing kernel"

function propagate!(wake::PanelWake, dt)
    for i_surf in eachindex(wake.nodes)
        view(wake.velocity[i_surf], 1:wake.nwakes[], :) .*= dt # displacements
        view(wake.nodes[i_surf], 1:wake.nwakes[], :) .+= view(wake.velocity[i_surf], 1:wake.nwakes[], :) # update nodes
        view(wake.velocity[i_surf], 1:wake.nwakes[], :) ./= dt # restore velocities
    end
end

function write_vtk(name, wake::PanelWake, idx, t)
    WriteVTK.paraview_collection(name) do pvd
        vtm = WriteVTK.vtk_multiblock(name * "$idx")
        for i_surf in eachindex(wake.nodes)
            WriteVTK.vtk_grid(vtm, name * "$idx", reshape(view(wake.nodes[i_surf], :, 1:wake.nwakes[], :), 3, wake.nwakes[], size(wake.nodes[i_surf], 2), 1)) do vtk
                vtk["velocity"] = wake.velocity[i_surf]
                vtk["strength"] = wake.strength[i_surf]
            end
        end
        pvd[t] = vtm
    end
end