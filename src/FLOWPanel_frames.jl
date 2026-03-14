#=##############################################################################
# DESCRIPTION
    Definition of reference frame objects for propagation of kinematics.

# AUTHORSHIP
  * Created by  : Ryan Anderson
  * Email       : Ry.M.Anderson@gmail.com
  * Date        : Feb 2026
  * License     : MIT License
=###############################################################################

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

struct ForwardRightDown end

struct BackRightUp end

function propagate_kinematics!(system::Union{AbstractBody}, frames::Vector{<:ReferenceFrame}, dt::Real)
    
    # translation vector from parent to global frame
    dx_parent_to_global = FastMultipole.SVector{3}(0.0, 0.0, 0.0)

    # global rotation matrix from parent to global frame
    R_parent_to_global = FastMultipole.SMatrix{3,3,Float64}(1.0,0,0,0,1.0,0,0,0,1.0)

    # begin recursion
    propagate_kinematics!(system, 1, frames, dx_parent_to_global, R_parent_to_global, dt)
end

function propagate_kinematics!(system::Union{AbstractBody}, i_frame::Int, frames::Vector{<:ReferenceFrame}, dx_parent_to_global, R_parent_to_global::FastMultipole.SMatrix, dt::Real)
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
        body = system isa MultiBody ? get_body(system, i) : system
        rotate_translate!(body, origin_global, Rω_global, dx_global)
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

function rotate_translate!(body::AbstractBody, origin, Rω, dx)
    # get grid nodes
    nodes = body.grid._nodes
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

function ReferenceFrame(system::Union{AbstractBody}; 
        # vvv all in global frame vvv
        origin = system.O,
        v = zero(FastMultipole.SVector{3,eltype(system.O)}),
        ω_axis = FastMultipole.SVector{3,eltype(system.O)}(0.0, 1.0, 0.0),
        ω = zero(eltype(system.O)),
        R = FastMultipole.SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        Rp2g = FastMultipole.SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0), # rotation from parent frame to global frame
        name = "vehicle",  # name of this frame
        child_index = Int[],  # indices of child frames
        dependent_index = system isa MultiBody ? collect(1:system.nbodies) : [1]  # indices of dependent bodies
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

#------- kinematic velocity -------#

function kinematic_velocity!(Vcp::AbstractMatrix, CPs::AbstractMatrix, system::Union{AbstractBody}, frames::AbstractVector{ReferenceFrame{TF}}; skip_top_level=false) where TF

    # capture the top level frame if requested
    if skip_top_level
        # begin recursion (skipping the top level frame, which is captured by system.fs)
        for i in frames[1].child_index
            frame = frames[1]
            kinematic_velocity!(Vcp, CPs, system, frames, i, frame.x, frame.R)
        end
    else
        kinematic_velocity!(Vcp, CPs, system, frames, 1, zero(FastMultipole.SVector{3,TF}), FastMultipole.SMatrix{3,3,TF,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    end

end

function kinematic_velocity!(Vcp::AbstractMatrix, CPs::AbstractMatrix, system::Union{AbstractBody}, frames::AbstractVector{<:ReferenceFrame}, i_frame::Int, dx_parent_to_global, R_parent_to_global)
    # get frame
    frame = frames[i_frame]

    # this frame's origin in global frame
    origin_global = R_parent_to_global * frame.x + dx_parent_to_global # global origin vector

    # this frame's velocity in global frame
    v_global = R_parent_to_global * frame.v # global velocity vector
    ω_global = R_parent_to_global * frame.ω_axis * frame.ω # global angular velocity
    
    # update the kinematic velocity of each dependent surface
    ncells_cum = 0
    for isurf in frame.dependent_index
        
        # unpack containers
        body = system isa MultiBody ? get_body(system, isurf) : system

        nc = body.ncells

        # velocity at the control points, we subtract this from the matrix because FLOWPanel expects rigid body EOM formulation: (U_inf + U_ind - U_kin) \cdot n = 0 
        for i in 1:nc
            icell = ncells_cum + i
            cp = FastMultipole.SVector{3}(CPs[1, icell], CPs[2, icell], CPs[3, icell])
            
            # vcp is actually total evaluated velocity, so subtracting kinematic velocity means rigid body motion opposes the freestream
            dv = v_global + cross(ω_global, (cp - origin_global))
            Vcp[1, icell] -= dv[1]
            Vcp[2, icell] -= dv[2]
            Vcp[3, icell] -= dv[3]
        end

        ncells_cum += nc

    end

    # new dx_parent_to_global
    dx_parent_to_global = origin_global

    # new R_parent_to_global
    R_parent_to_global = R_parent_to_global * frame.R
    
    # propagate to child frames
    for i in frame.child_index
        kinematic_velocity!(Vcp, CPs, system, frames, i, dx_parent_to_global, R_parent_to_global)
    end
end


# function Freestream(frame::ReferenceFrame{TF}, ref::Reference, vinf_ext) where TF
#     # equivalent freestream about ref.r
#     dx = ref.r - frame.x  # vector from frame origin to reference point
#     ω = frame.ω_axis * frame.ω
#     V = vinf_ext - frame.v - ω × dx  # freestream velocity in South-East-Up convention

#     # create freestream object
#     Omega = frame.ω_axis * frame.ω
#     fs = velocity_to_freestream(V, -Omega)

#     return fs
# end

#------- constructors -------#

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
