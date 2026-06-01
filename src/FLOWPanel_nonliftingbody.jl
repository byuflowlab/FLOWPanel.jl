#=##############################################################################
# DESCRIPTION
    Definition of non-lifting paneled body types (implementations of
    AbstractBody).

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Jun 2018
  * License     : MIT License
=###############################################################################


################################################################################
# NON-LIFTING BODY TYPE
################################################################################

"""
    NonLiftingBody{E, N}
    NonLiftingBody{E}(mesh; DBC=false, optargs...)
    NonLiftingBody{E}(nodes, cells; DBC=false, optargs...)

Concrete body type for non-lifting surfaces discretized with source, doublet,
or vortex-ring panels. Constructors accept a `VSPGeom` mesh or raw node/cell
arrays.
"""
mutable struct NonLiftingBody{E, N, TF, DBC} <: AbstractBody{E, N, TF, DBC}

    # User inputs
    nodes::Matrix{TF}                         # 3xnnodes matrix where nodes[:, i] is the position of the i-th node
    
    # Properties
    vtk_cells::Vector{WriteVTK.MeshCell{WriteVTK.VTKCellTypes.VTKCellType, Vector{Int64}}}      # Vector of WriteVTK cells
    neighbor::Matrix{Int}                     # 3xncells matrix where neighbor[i, j] is the linear index of the cell neighboring the i-th edge of the j-th cell (or 0 if it's a boundary)
    nnodes::Int                               # Number of nodes
    ncells::Int                               # Number of cells
    cells::Matrix{Int}                        # Cell connectivity (each column is a cell)
    Oaxis::Array{TF,2}                  # Coordinate system of original grid
    O::Array{TF,1}                      # Position of CS of original grid

    # Internal variables
    strength::Array{TF, 2}              # strength[i,j] is the stength of the i-th panel with the j-th element type
    potential::Array{TF,1}              # Potential at control points
    velocity::Array{TF,2}               # Apparent fluid velocity at control points (body frame)
    velocity_gradient::Array{TF,3}      # 3x3xncells velocity gradient at control points; only populated when needs_velocity_gradient[]
    induced_vorticity::Matrix{TF}       # 3xncells vorticity at control points; bound surface vorticity initialized in simulate!, then populated by extra_outputs=3
    velocity_kinematic::Matrix{TF}      # Rigid-body kinematic velocity at control points (inertial frame)
    angular_velocity::Vector{TF}        # Net angular velocity (global frame), sum over ancestor frames; populated by kinematic_velocity!
    controlpoints::Matrix{TF}           # 3xncells control points
    normals::Matrix{TF}                 # 3xncells panel normals
    cp_outer::Bool                      # Side of the surface the control point limit is taken from (true = exterior, false = interior)
    kerneloffset::Float64               # Active kernel offset to avoid singularities
    kerneloffset_panel::Float64         # Kernel offset for panel solves/interactions
    kerneloffset_targets::Float64       # Kernel offset for panel influence on targets
    kernelcutoff::Float64               # Kernel cutoff to avoid singularities
    characteristiclength::Function      # Characteristic length of each panel
    watertight::Bool                     # Whether the body is watertight or not
    inside_offset::Float64               # Offset to compute inside control points
    needs_velocity_gradient::Base.RefValue{Bool}  # Set by simulate! when a monitor requires ∇u
end

function NonLiftingBody{E, N, TF, DBC}(
                nodes::Matrix{TF}, cells::Matrix{Int};
                vtk_cells=nothing,
                neighbor=nothing,
                nnodes=size(nodes, 2), ncells=size(cells, 2),
                Oaxis=Array{TF,2}(1.0I, 3, 3), O=zeros(TF,3),
                strength=zeros(size(cells, 2), N),
                potential=zeros(size(cells, 2)),
                velocity=zeros(3, size(cells, 2)),
                velocity_gradient=zeros(TF, 3, 3, size(cells, 2)),
                induced_vorticity=zeros(TF, 3, size(cells, 2)),
                velocity_kinematic=zeros(TF, 3, size(cells, 2)),
                angular_velocity=zeros(TF, 3),
                controlpoints=zeros(3, size(cells, 2)),
                normals=zeros(3, size(cells, 2)),
                cp_outer::Bool=!DBC,
                kerneloffset=1e-8,
                kerneloffset_panel=kerneloffset,
                kerneloffset_targets=kerneloffset,
                kernelcutoff=1e-14,
                characteristiclength=characteristiclength_sqrtarea,
                watertight=false,
                inside_offset=1e-6,
                needs_velocity_gradient=Ref(false),
                ensure_winding::Bool=true,
                flip_normals::Bool=false
              ) where {E, N, TF, DBC}
    if ensure_winding 
        ensure_consistent_winding!(cells, nodes; watertight, flip_normals)
    end
    neighbor = neighbor === nothing || ensure_winding ? calc_neighbors(cells) : neighbor
    vtk_cells = vtk_cells === nothing || ensure_winding ?
                [WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, cells[:, i]) for i in 1:size(cells, 2)] :
                vtk_cells
    ncells = size(cells, 2)

    return NonLiftingBody{E, N, TF, DBC}(
                nodes, vtk_cells, neighbor,
                nnodes, ncells, cells,
                Oaxis, O,
                strength,
                potential,
                velocity,
                velocity_gradient,
                induced_vorticity,
                velocity_kinematic,
                angular_velocity,
                controlpoints,
                normals,
                cp_outer,
                kerneloffset_panel,
                Float64(kerneloffset_panel),
                Float64(kerneloffset_targets),
                kernelcutoff,
                characteristiclength,
                watertight,
                inside_offset,
                needs_velocity_gradient
              )
end

function NonLiftingBody{E, N, TF}(
                nodes::Matrix{TF}, cells::Matrix{Int};
                DBC::Bool=false,
                optargs...
              ) where {E, N, TF}
    return NonLiftingBody{E, N, TF, DBC}(nodes, cells; optargs...)
end

function NonLiftingBody{E, N, TF, DBC}(
                mesh::VSPGeom.TriMesh;
                optargs...
              ) where {E, N, TF, DBC}
    nodes, cells = trimesh2cells(mesh)
    return NonLiftingBody{E, N, TF, DBC}(TF.(nodes), cells; optargs...)
end

function NonLiftingBody{E, N, TF}(
                mesh::VSPGeom.TriMesh;
                DBC::Bool=false,
                optargs...
              ) where {E, N, TF}
    return NonLiftingBody{E, N, TF, DBC}(mesh; optargs...)
end

_count(::Type{ConstantSource}) = 1
_count(::Type{ConstantDoublet}) = 1
_count(::Type{VortexRing}) = 1
_count(::Type{Union{ConstantSource, ConstantDoublet}}) = 2
_count(::Type{Union{ConstantSource, VortexRing}}) = 2
_count(::Type{<:Any}) = error("Unsupported kernel type for NonLiftingBody.")

function (NonLiftingBody{E})(mesh::VSPGeom.TriMesh; DBC::Bool=false, optargs...) where {E}
    nodes, _ = trimesh2cells(mesh)
    return NonLiftingBody{E, _count(E), eltype(nodes), DBC}(mesh; optargs...)
end

function (NonLiftingBody{E, N})(mesh::VSPGeom.TriMesh; DBC::Bool=false, optargs...) where {E, N}
    nodes, _ = trimesh2cells(mesh)
    return NonLiftingBody{E, N, eltype(nodes), DBC}(mesh; optargs...)
end

function (NonLiftingBody{E})(nodes::Matrix{TF}, cells::Matrix{Int}; DBC::Bool=false, optargs...) where {E, TF}
    return NonLiftingBody{E, _count(E), TF, DBC}(nodes, cells; optargs...)
end

function (NonLiftingBody{E, N})(nodes::Matrix{TF}, cells::Matrix{Int}; DBC::Bool=false, optargs...) where {E, N, TF}
    return NonLiftingBody{E, N, TF, DBC}(nodes, cells; optargs...)
end

"""
    save(body::NonLiftingBody, args...; optargs...)

Write a non-lifting body and its stored solution fields using FLOWPanel's VTK
export path.
"""
function save(body::NonLiftingBody, args...; optargs...)
    return save_base(body, args...; optargs...)
end

solved_field_name(::NonLiftingBody{ConstantSource, 1}) = "sigma"
solved_field_name(::NonLiftingBody{ConstantDoublet, 1}) = "mu"
solved_field_name(::NonLiftingBody{VortexRing, 1}) = "gamma"
solved_field_name(::NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2}) = "mu"

#### END OF NON-LIFTING BODY  ##################################################


################################################################################
# CONSTANT-SOURCE SOLVER
################################################################################

_get_Gdims(self::NonLiftingBody{ConstantSource, 1}) = (self.ncells, self.ncells)



################################################################################
# CONSTANT-DOUBLET SOLVER
################################################################################

_get_Gdims(self::NonLiftingBody{ConstantDoublet, 1}) = (self.ncells, self.ncells)



################################################################################
# CONSTANT-SOURCE+DOUBLET SOLVER
################################################################################

_get_Gdims(self::NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2}) = (self.ncells, self.ncells)




################################################################################
# FASTMULTIPOLE BACKEND SUPPORT
################################################################################

FastMultipole.has_vector_potential(::AbstractBody{ConstantSource, 1}) = false

FastMultipole.has_vector_potential(::AbstractBody{ConstantDoublet, 1}) = false

FastMultipole.has_vector_potential(::AbstractBody{Union{ConstantSource, ConstantDoublet}, 2, <:Any}) = false

FastMultipole.has_vector_potential(::AbstractBody{Union{ConstantSource, VortexRing}, 2, <:Any}) = false

FastMultipole.body_to_multipole!(system::AbstractBody{ConstantSource, 1, <:Any}, args...) =
    FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.Source}, system, args...)

FastMultipole.body_to_multipole!(system::AbstractBody{ConstantDoublet, 1, <:Any}, args...) =
    FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.Dipole}, system, args...; scale_strength=1.0)

FastMultipole.body_to_multipole!(system::AbstractBody{Union{ConstantSource,ConstantDoublet}, 2, <:Any}, args...) =
    FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.SourceDipole}, system, args...; scale_strength=FastMultipole.StaticArrays.SVector(1.0, 1.0))

FastMultipole.body_to_multipole!(system::AbstractBody{Union{ConstantSource,VortexRing}, 2, <:Any}, args...) =
    FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.SourceDipole}, system, args...; scale_strength=FastMultipole.StaticArrays.SVector(1.0, 1.0))

function FastMultipole.strength_to_value(strength, source_system::NonLiftingBody{<:Any, 1, <:Any})
    return strength[1]
end

##### END OF FASTMULTIPOLE BACKEND SUPPORT #####################################

function FastMultipole.value_to_strength!(source_buffer, ::NonLiftingBody, i_body, value, rlx)
    prev_value = source_buffer[5, i_body]
    source_buffer[5, i_body] = rlx*value + (1.0 - rlx)*prev_value
end

function FastMultipole.value_to_strength!(source_buffer, ::NonLiftingBody, i_body, value)
    source_buffer[5, i_body] = value
end

function FastMultipole.buffer_to_target_system!(target_system::NonLiftingBody, i_target, switch::FastMultipole.DerivativesSwitch{PS,VS,GS,NO,NM}, target_buffer, i_buffer) where {PS,VS,GS,NO,NM}
    if PS
        phi = FastMultipole.get_scalar_potential(target_buffer, switch, i_buffer)
        target_system.potential[i_target] += phi
    end

    if VS
        vx, vy, vz = FastMultipole.get_gradient(target_buffer, switch, i_buffer)
        target_system.velocity[1, i_target] += vx
        target_system.velocity[2, i_target] += vy
        target_system.velocity[3, i_target] += vz
    end

    # Accumulate the 3x3 velocity gradient from rows 8..16 (column-major
    # SMatrix layout per FastMultipole.get_hessian).
    if GS
        H = FastMultipole.get_hessian(target_buffer, switch, i_buffer)
        @inbounds for j in 1:3, i in 1:3
            target_system.velocity_gradient[i, j, i_target] += H[i, j]
        end
    end

    if NO == 3
        @inbounds for j in 1:3
            target_system.induced_vorticity[j, i_target] +=
                FastMultipole.get_extra_output(target_buffer, switch, i_buffer, j)
        end
    end
end
################################################################################
# COMMON FUNCTIONS
################################################################################

##### END OF COMMON FUNTIONS ###################################################

"""
    FlatGround(center, normal, radius; panel_length=radius/5, bodyoptargs...)

Generate a flat circular ground plane tessellated with equilateral triangles of
uniform size `panel_length`. Returns a `NonLiftingBody{ConstantSource}` suitable
for use with `FlatGroundSolver`.

# Arguments
- `center`: 3-element vector, center position of the ground plane.
- `normal`: 3-element vector, outward normal direction of the ground plane.
- `radius`: scalar, radius of the circular ground plane.

# Keyword Arguments
- `panel_length`: side length of each equilateral triangle (default `radius/5`).
- `bodyoptargs...`: forwarded to the `NonLiftingBody{ConstantSource}` constructor.
"""
function FlatGround(center, normal, radius; panel_length=radius/5, bodyoptargs...)

    TF = promote_type(eltype(center), eltype(normal), typeof(float(radius)), typeof(float(panel_length)))

    # --- Step 1: Generate nodes on an equilateral lattice in the XY plane ---
    h = TF(panel_length)
    e1 = TF[h, 0]
    e2 = TF[h/2, h * sqrt(TF(3))/2]

    N = ceil(Int, radius / h) + 1

    # Map from (i,j) grid index to node index
    node_map = Dict{Tuple{Int,Int}, Int}()
    node_list = Vector{TF}[]

    for j in -N:N, i in -N:N
        x = i * e1[1] + j * e2[1]
        y = i * e1[2] + j * e2[2]
        if x^2 + y^2 <= radius^2
            push!(node_list, TF[x, y, 0])
            node_map[(i, j)] = length(node_list)
        end
    end

    nnodes = length(node_list)

    # --- Step 2: Build equilateral triangles ---
    cell_list = Vector{Int}[]

    for j in -N:N-1, i in -N:N-1
        # Lower triangle: (i,j), (i+1,j), (i,j+1)
        if haskey(node_map, (i,j)) && haskey(node_map, (i+1,j)) && haskey(node_map, (i,j+1))
            push!(cell_list, [node_map[(i,j)], node_map[(i+1,j)], node_map[(i,j+1)]])
        end
        # Upper triangle: (i+1,j), (i+1,j+1), (i,j+1)
        if haskey(node_map, (i+1,j)) && haskey(node_map, (i+1,j+1)) && haskey(node_map, (i,j+1))
            push!(cell_list, [node_map[(i+1,j)], node_map[(i+1,j+1)], node_map[(i,j+1)]])
        end
    end

    ncells = length(cell_list)

    # Assemble into matrices
    nodes_2d = zeros(TF, 3, nnodes)
    for k in 1:nnodes
        nodes_2d[:, k] .= node_list[k]
    end

    cells = zeros(Int, 3, ncells)
    for k in 1:ncells
        cells[:, k] .= cell_list[k]
    end

    # --- Step 3: Rotate and translate to match desired center and normal ---
    n_hat = TF.(normal) / sqrt(sum(TF.(normal).^2))
    from = TF[0, 0, 1]

    axis = cross(from, n_hat)
    sin_angle = sqrt(sum(axis.^2))
    cos_angle = dot(from, n_hat)

    if sin_angle < 1e-12
        # Parallel or anti-parallel
        R = cos_angle > 0 ? Matrix{TF}(1.0I, 3, 3) : TF[-1 0 0; 0 -1 0; 0 0 1]
    else
        axis = axis / sin_angle
        # Rodrigues' rotation formula
        K = TF[0 -axis[3] axis[2]; axis[3] 0 -axis[1]; -axis[2] axis[1] 0]
        R = Matrix{TF}(1.0I, 3, 3) + sin_angle * K + (1 - cos_angle) * (K * K)
    end

    c = TF.(center)
    nodes = R * nodes_2d .+ c

    # --- Step 4: Ensure consistent winding (normals align with desired normal) ---
    v1 = nodes[:, cells[1, 1]]
    v2 = nodes[:, cells[2, 1]]
    v3 = nodes[:, cells[3, 1]]
    computed_normal = cross(v2 - v1, v3 - v1)

    if dot(computed_normal, n_hat) < 0
        # Swap vertex order to flip normals
        cells[[1, 2], :] .= cells[[2, 1], :]
    end

    # --- Step 5: Construct and return body ---
    return NonLiftingBody{ConstantSource}(nodes, cells; bodyoptargs...)
end
