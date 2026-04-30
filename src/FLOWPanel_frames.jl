#=##############################################################################
# DESCRIPTION
    Definition of reference frame objects for propagation of kinematics.

# AUTHORSHIP
  * Created by  : Ryan Anderson
  * Email       : Ry.M.Anderson@gmail.com
  * Date        : Feb 2026
  * License     : MIT License
=###############################################################################

"""
    ReferenceFrame

Rigid-body kinematic frame used to propagate body motion through a frame tree.
"""
struct ReferenceFrame{TF}
	x::FastMultipole.SVector{3,TF}             # origin in parent frame
	v::FastMultipole.SVector{3,TF}             # velocity in parent frame
	ω_axis::FastMultipole.SVector{3,TF}        # axis of rotation in parent frame
	ω::TF                        # angular velocity in parent frame
	R::FastMultipole.SMatrix{3,3,TF,9}         # basis vectors expressed in parent frame
    Rp2g::FastMultipole.SMatrix{3,3,TF,9}      # parent frame basis vectors expressed in global frame
	name::String                 # name of this frame
	parent_index::Int            # index of parent frame
	child_index::Vector{Int}  # child reference frames
	dependent_index::Vector{Int} # VortexLattice surfaces
end

"""
    ForwardRightDown

Marker type for a forward-right-down coordinate convention.
"""
struct ForwardRightDown end

"""
    BackRightUp

Marker type for a back-right-up coordinate convention.
"""
struct BackRightUp end

"""
    propagate_kinematics!(systems, frames, dt)

Advance all dependent bodies by one timestep `dt` according to the frame tree.
Accepts either a single body or a tuple of bodies.
"""
function propagate_kinematics!(system::Union{AbstractBody}, frames::Vector{<:ReferenceFrame}, dt::Real)

    # translation vector from parent to global frame
    dx_parent_to_global = FastMultipole.SVector{3}(0.0, 0.0, 0.0)

    # global rotation matrix from parent to global frame
    R_parent_to_global = FastMultipole.SMatrix{3,3,Float64}(1.0,0,0,0,1.0,0,0,0,1.0)

    # begin recursion
    propagate_kinematics!(system, 1, frames, dx_parent_to_global, R_parent_to_global, dt)
end

function propagate_kinematics!(system::AbstractBody, i_frame::Int, frames::Vector{<:ReferenceFrame}, dx_parent_to_global, R_parent_to_global::FastMultipole.SMatrix, dt::Real)
    # get frame
    frame = frames[i_frame]

    # origin in global frame
    origin_global = R_parent_to_global * frame.x + dx_parent_to_global # global origin vector

    # differential translation
    dx = frame.v * dt  # translation in parent frame
    dx_global = R_parent_to_global * dx  # global translation vector

    # differential rotation
    dω = frame.ω * dt  # angular displacement in parent frame
    Rω = Rodrigues(frame.ω_axis, dω) # rotation matrix in parent frame
    Rω_global = Rodrigues(R_parent_to_global * frame.ω_axis, dω) # global frame

    # rotate and translate dependent surfaces
    for i in frame.dependent_index
        body = system
        rotate_translate!(body, origin_global, Rω_global, dx_global)
        rotate_Das!(body, Rω_global)
    end

    # Update the frame
    x_new = frame.x + dx
    R_new = Rω * frame.R
    frames[i_frame] = ReferenceFrame(x_new, frame.v, frame.ω_axis, frame.ω, R_new, R_parent_to_global, frame.name, frame.parent_index, frame.child_index, frame.dependent_index)

    # new dx_parent_to_global
    dx_parent_to_global = origin_global + dx_global

    # new R_parent_to_global
    R_parent_to_global = R_parent_to_global * R_new

    # Recursively propagate to child frames
    for i in frame.child_index
        propagate_kinematics!(system, i, frames, dx_parent_to_global, R_parent_to_global, dt)
    end
end

#--- tuple overloads ---#

function propagate_kinematics!(systems::Tuple, frames::Vector{<:ReferenceFrame}, dt::Real)

    # translation vector from parent to global frame
    dx_parent_to_global = FastMultipole.SVector{3}(0.0, 0.0, 0.0)

    # global rotation matrix from parent to global frame
    R_parent_to_global = FastMultipole.SMatrix{3,3,Float64}(1.0,0,0,0,1.0,0,0,0,1.0)

    # begin recursion
    propagate_kinematics!(systems, 1, frames, dx_parent_to_global, R_parent_to_global, dt)
end

function propagate_kinematics!(systems::Tuple, i_frame::Int, frames::Vector{<:ReferenceFrame}, dx_parent_to_global, R_parent_to_global::FastMultipole.SMatrix, dt::Real)
    # get frame
    frame = frames[i_frame]

    # origin in global frame
    origin_global = R_parent_to_global * frame.x + dx_parent_to_global # global origin vector

    # differential translation
    dx = frame.v * dt  # translation in parent frame
    dx_global = R_parent_to_global * dx  # global translation vector

    # differential rotation
    dω = frame.ω * dt  # angular displacement in parent frame
    Rω = Rodrigues(frame.ω_axis, dω) # rotation matrix in parent frame
    Rω_global = Rodrigues(R_parent_to_global * frame.ω_axis, dω) # global frame

    # rotate and translate dependent surfaces
    for i in frame.dependent_index
        body = systems[i]
        rotate_translate!(body, origin_global, Rω_global, dx_global)
        rotate_Das!(body, Rω_global)
    end

    # Update the frame
    x_new = frame.x + dx
    R_new = Rω * frame.R
    frames[i_frame] = ReferenceFrame(x_new, frame.v, frame.ω_axis, frame.ω, R_new, R_parent_to_global, frame.name, frame.parent_index, frame.child_index, frame.dependent_index)

    # new dx_parent_to_global
    dx_parent_to_global = origin_global + dx_global

    # new R_parent_to_global
    R_parent_to_global = R_parent_to_global * R_new

    # Recursively propagate to child frames
    for i in frame.child_index
        propagate_kinematics!(systems, i, frames, dx_parent_to_global, R_parent_to_global, dt)
    end
end

function rotate_translate!(body::AbstractBody, origin, Rω, dx)
    # get grid nodes
    nodes = body.nodes
    nnodes = size(nodes, 2)

    # rotate/translate
    for i in 1:nnodes
        # relative to origin
        for k in 1:3
            nodes[k,i] -= origin[k]
        end

        # rotate
        nodes[:,i] .= Rω * FastMultipole.SVector{3}(nodes[1,i], nodes[2,i], nodes[3,i])

        # translate and shift origin back
        for k in 1:3
            nodes[k,i] += dx[k] + origin[k]
        end
    end
end

function rotate_Das!(body::AbstractLiftingBody, Rω)
    for das in body.Das
        for j in axes(das, 2)
            das[:, j] .= Rω * FastMultipole.SVector{3}(das[1,j], das[2,j], das[3,j])
        end
    end
end
rotate_Das!(body::AbstractBody, Rω) = nothing

"""
    Rodrigues(axis, angle)

Return the rotation matrix associated with rotating by `angle` about `axis`.
"""
function Rodrigues(axis, angle::TF) where TF
    s, c = sincos(-angle)
    t = 1 - c
    x, y, z = axis

    return FastMultipole.SMatrix{3,3,TF,9}(
        t*x*x + c,     t*x*y - s*z, t*x*z + s*y,
        t*x*y + s*z,   t*y*y + c,   t*y*z - s*x,
        t*x*z - s*y,   t*y*z + s*x, t*z*z + c
    )
end

"""
    inverse_Rodrigues(R)

Return the axis-angle rotation vector associated with the rotation matrix `R`.
"""
function inverse_Rodrigues(R::FastMultipole.SMatrix{3,3,TF,9}) where TF
    # inverse Rodrigues rotation matrix
    θ = acos((trace(R) - 1) / 2)
    if θ == 0.0
        return FastMultipole.SVector{3,TF}(0.0, 0.0, 0.0)
    end
    s = sin(θ)
    x = (R[3,2] - R[2,3]) / (2 * s)
    y = (R[1,3] - R[3,1]) / (2 * s)
    z = (R[2,1] - R[1,2]) / (2 * s)
    return FastMultipole.SVector{3,TF}(x, y, z) * θ
end

"""
    ReferenceFrame(system; origin=system.O, v=..., ω_axis=..., ω=..., R=..., Rp2g=..., name="vehicle", child_index=Int[], dependent_index=[1])

Construct the root frame list for a body or body collection.
"""
function ReferenceFrame(system::AbstractBody; 
        # vvv all in global frame vvv
        origin = system.O,
        v = zero(FastMultipole.SVector{3,eltype(system.O)}),
        ω_axis = FastMultipole.SVector{3,eltype(system.O)}(0.0, 1.0, 0.0),
        ω = zero(eltype(system.O)),
        R = FastMultipole.SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        Rp2g = FastMultipole.SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0), # rotation from parent frame to global frame
        name = "vehicle",  # name of this frame
        child_index = Int[],  # indices of child frames
        dependent_index = [1]  # indices of dependent bodies
        # ^^^ all in global frame ^^^
    )
    TF = eltype(system.O)
    parent_index = -1  # no parent frame
    vehicle_frame = ReferenceFrame{TF}(
        origin, v, ω_axis, ω, R, Rp2g,
        name, parent_index, child_index, dependent_index
    )
    frames = [vehicle_frame]
    return frames
end

function update_dependent_indices!(frames::Vector{ReferenceFrame{TF}}, parent_index::Int, surface_indices::Vector{Int}) where TF
    # update dependent indices for the parent frame
    for i in surface_indices
        if !(i in frames[parent_index].dependent_index)
            push!(frames[parent_index].dependent_index, i)
        end
    end
    
    # recursively update parent frames
    grandparent_index = frames[parent_index].parent_index
    if grandparent_index != -1
        update_dependent_indices!(frames, grandparent_index, surface_indices)
    end
end

#------- kinematic velocity helpers -------#

function _kinematic_velocity_te!(body::AbstractLiftingBody, v_global, ω_global, origin_global)
    for ishedding in eachindex(body.Das)
        Vte = body.velocity_te[ishedding]
        shedding = body.shedding[ishedding]

        # nib nodes (columns 1..nshed)
        for j in axes(shedding, 2)
            i_panel = shedding[1, j]
            idx_nib = shedding[3, j]
            node_idx = body.cells[idx_nib, i_panel]
            te = FastMultipole.SVector{3}(body.nodes[1, node_idx], body.nodes[2, node_idx], body.nodes[3, node_idx])
            dv = v_global + cross(ω_global, (te - origin_global))
            Vte[1, j] -= dv[1]
            Vte[2, j] -= dv[2]
            Vte[3, j] -= dv[3]
        end

        # final nia node (column nshed+1)
        if size(shedding, 2) > 0
            i_panel = shedding[1, end]
            idx_nia = shedding[2, end]
            node_idx = body.cells[idx_nia, i_panel]
            te = FastMultipole.SVector{3}(body.nodes[1, node_idx], body.nodes[2, node_idx], body.nodes[3, node_idx])
            dv = v_global + cross(ω_global, (te - origin_global))
            Vte[1, end] -= dv[1]
            Vte[2, end] -= dv[2]
            Vte[3, end] -= dv[3]
        end
    end
end

_kinematic_velocity_te!(::AbstractBody, v_global, ω_global, origin_global) = nothing

#------- kinematic velocity -------#

kinematic_velocity!(system::AbstractBody, frames::AbstractVector{<:ReferenceFrame}) = kinematic_velocity!((system,), frames)

function kinematic_velocity!(systems::Tuple, frames::AbstractVector{ReferenceFrame{TF}}) where TF

    # begin recursion
    kinematic_velocity!(systems, frames, 1, zero(FastMultipole.SVector{3,TF}), FastMultipole.SMatrix{3,3,TF,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))

end

function kinematic_velocity!(systems::Tuple, frames::AbstractVector{ReferenceFrame{TF}}, i_frame::Int, dx_parent_to_global, R_parent_to_global) where TF
    # get frame
    frame = frames[i_frame]

    # this frame's origin in global frame
    origin_global = R_parent_to_global * frame.x + dx_parent_to_global # global origin vector

    # this frame's velocity in global frame
    v_global = R_parent_to_global * frame.v # global velocity vector
    ω_global = R_parent_to_global * frame.ω_axis * frame.ω # global angular velocity
    
    # update the kinematic velocity of each dependent surface
    for isurf in frame.dependent_index
        
        # unpack containers
        body = systems[isurf]
        Vcp = body.velocity
        CPs = body.controlpoints

        # velocity at the control points
        for i in axes(Vcp, 2)

            # extract control point
            cp = FastMultipole.SVector{3}(CPs[1, i], CPs[2, i], CPs[3, i])
            
            # vcp is actually total evaluated velocity, so subtracting kinematic velocity means rigid body motion opposes the freestream
            dv = v_global + cross(ω_global, (cp - origin_global))
            Vcp[:, i] .-= dv
        end

        # trailing edges (only for lifting bodies)
        _kinematic_velocity_te!(body, v_global, ω_global, origin_global)
    end

    # new dx_parent_to_global
    dx_parent_to_global = origin_global

    # new R_parent_to_global
    R_parent_to_global = R_parent_to_global * frame.R
    
    # propagate to child frames
    for i in frame.child_index
        kinematic_velocity!(systems, frames, i, dx_parent_to_global, R_parent_to_global)
    end
end


#------- constructors -------#

"""
    add_frame!(frames, name, parent_index, origin, surface_indices; optargs...)
    add_frame!(frames, name, parent_name, origin, surface_indices; optargs...)

Add a child frame to an existing frame tree and register the affected body
indices as dependents of that frame and its ancestors.
"""
function add_frame!(frames::Vector{ReferenceFrame{TF}}, name::String, parent_index::Int, origin, surface_indices::Vector{Int};
    v = zero(FastMultipole.SVector{3,TF}),  # velocity in parent frame
    ω_axis = FastMultipole.SVector{3,TF}(0.0, 1.0, 0.0),  # axis of rotation in parent frame
    ω = zero(TF),  # angular velocity in parent frame
    R = FastMultipole.SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),  # basis vectors expressed in parent frame
    Rp2g = FastMultipole.SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)  # rotation from parent frame to global frame
) where TF

    # create new frame
    new_frame = ReferenceFrame{TF}(
        origin, v, ω_axis, ω, R, Rp2g,
        name, parent_index, Int[], surface_indices
    )
    push!(frames, new_frame)

    # inform parents
    push!(frames[parent_index].child_index, length(frames))
    update_dependent_indices!(frames, parent_index, surface_indices)
end

function add_frame!(frames::Vector{<:ReferenceFrame}, name::String, parent_name::String, origin, surface_indices::Vector{Int}; optargs...)
    for (i, frame) in enumerate(frames)
        if frame.name == parent_name
            return add_frame!(frames, name, i, origin, surface_indices; optargs...)
        end
    end
end

#------- frame queries -------#

"""
    frame_global_transform(frames, i_frame)

Return `(origin_global, R_frame_to_global)` for the frame at index `i_frame`.

`origin_global` is the frame origin expressed in global coordinates and
`R_frame_to_global` is the rotation matrix that maps vectors in the frame's
local coordinate system to the global coordinate system.
"""
function frame_global_transform(frames::AbstractVector{ReferenceFrame{TF}}, i_target::Int) where TF
    dx = zero(FastMultipole.SVector{3,TF})
    Rp2g = FastMultipole.SMatrix{3,3,TF,9}(1,0,0, 0,1,0, 0,0,1)
    return _frame_global_transform(frames, 1, i_target, dx, Rp2g)
end

function _frame_global_transform(frames::AbstractVector{ReferenceFrame{TF}}, i_frame::Int, i_target::Int, dx_parent_to_global, R_parent_to_global) where TF
    frame = frames[i_frame]

    # origin of this frame in global coordinates
    origin_global = R_parent_to_global * frame.x + dx_parent_to_global

    # rotation from this frame's local axes to global
    R_frame_to_global = R_parent_to_global * frame.R

    # found the target
    if i_frame == i_target
        return (origin_global, R_frame_to_global)
    end

    # recurse into children
    for i_child in frame.child_index
        result = _frame_global_transform(frames, i_child, i_target, origin_global, R_frame_to_global)
        if !isnothing(result)
            return result
        end
    end

    return nothing
end
