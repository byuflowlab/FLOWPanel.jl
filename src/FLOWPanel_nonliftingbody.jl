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
  `NonLiftingBody{E::AbstractElement, N}(grid::gt.GridTriangleSurface)`

Non-lifting body that is solved using a combination of N panel elements.
`grid` is the grid surface (paneled geometry).

  **Properties**
  * `nnodes::Int`                       : Number of nodes
  * `ncells::Int`                       : Number of cells
  * `fields::Vector{String}`            : Available fields (solutions)
  * `Oaxis::Matrix`                     : Coordinate system of body w.r.t. global
  * `O::Vector`                         : Origin of body w.r.t. global

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

    # Fields
    Cp::Vector{TF}
    F::Matrix{TF}

    # Internal variables
    strength::Array{TF, 2}              # strength[i,j] is the stength of the i-th panel with the j-th element type
    potential::Array{TF,1}              # Potential at control points
    velocity::Array{TF,2}               # Velocity at control points
    controlpoints::Matrix{TF}           # 3xncells control points
    normals::Matrix{TF}                 # 3xncells panel normals
    CPoffset::Float64                   # Control point offset in normal direction
    kerneloffset::Float64               # Kernel offset to avoid singularities
    kernelcutoff::Float64               # Kernel cutoff to avoid singularities
    characteristiclength::Function      # Characteristic length of each panel
    watertight::Bool                     # Whether the body is watertight or not
    inside_offset::Float64               # Offset to compute inside control points
end

function NonLiftingBody{E, N, TF, DBC}(
                nodes::Matrix{TF}, cells::Matrix{Int};
                vtk_cells::Vector{<:WriteVTK.MeshCell}=[WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, cells[:, i]) for i in 1:size(cells, 2)],
                neighbor::Matrix{Int}=zeros(Int, 3, size(cells, 2)),
                nnodes=size(nodes, 2), ncells=size(cells, 2),
                Oaxis=Array{TF,2}(1.0I, 3, 3), O=zeros(TF,3),
                Cp=zeros(TF, size(cells, 2)),
                F=zeros(TF, 3, size(cells, 2)),
                strength=zeros(size(cells, 2), N),
                potential=zeros(size(cells, 2)),
                velocity=zeros(3, size(cells, 2)),
                controlpoints=zeros(3, size(cells, 2)),
                normals=zeros(3, size(cells, 2)),
                CPoffset=1e-14,
                kerneloffset=1e-8,
                kernelcutoff=1e-14,
                characteristiclength=characteristiclength_unitary,
                watertight=false,
                inside_offset=1e-6
              ) where {E, N, TF, DBC}

    return NonLiftingBody{E, N, TF, DBC}(
                nodes, vtk_cells, neighbor,
                nnodes, ncells, cells,
                Oaxis, O,
                Cp, F,
                strength,
                potential,
                velocity,
                controlpoints,
                normals,
                CPoffset,
                kerneloffset,
                kernelcutoff,
                characteristiclength,
                watertight,
                inside_offset
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
                grid::gt.GridTriangleSurface;
                optargs...
              ) where {E, N, TF, DBC}
    
    nodes = grid._nodes
    cells = grid2cells(grid)
    
    # Extract neighbor info from grid
    neighbor = zeros(Int, 3, grid.ncells)
    
    ndivscellsc = Tuple(collect( 1:(d != 0 ? d : 1) for d in grid._ndivscells))
    linc = LinearIndices(ndivscellsc)
    
    for ci in 1:grid.ncells
        for ni in 1:3                   # Iterate over neighbors
            ncoor = gt.neighbor(grid, ni, ci; preserveEdge=true)
            if ncoor[1] != 0
                nlin = linc[ncoor...]
                neighbor[ni, ci] = nlin
            end
        end
    end

    vtk_cells = [WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, cells[:, i]) for i in 1:grid.ncells]

    # check if mesh is watertight
    watertight_guess = false
    if typeof(grid.orggrid) <: gt.Meshes.Mesh
        mesh = grid.orggrid
        watertight_guess = gt.isclosed(mesh)
    end
    
    return NonLiftingBody{E, N, TF, DBC}(
                nodes, cells;
                vtk_cells=vtk_cells, neighbor=neighbor, watertight=watertight_guess, optargs...
              )
end

function NonLiftingBody{E, N, TF}(
                grid::gt.GridTriangleSurface;
                DBC::Bool=false,
                optargs...
              ) where {E, N, TF}
    return NonLiftingBody{E, N, TF, DBC}(grid; optargs...)
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

function (NonLiftingBody{E})(grid::gt.GridTriangleSurface; DBC::Bool=false, optargs...) where {E}
    return NonLiftingBody{E, _count(E), eltype(grid._nodes), DBC}(grid; optargs...)
end

function (NonLiftingBody{E, N})(grid::gt.GridTriangleSurface; DBC::Bool=false, optargs...) where {E, N}
    return NonLiftingBody{E, N, eltype(grid._nodes), DBC}(grid; optargs...)
end

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

function _G_U!(self::AbstractBody{<:Any,NK,TF}, kernel, G, CPs, normals; strength_index=kernel==ConstantDoublet && NK>1 ? 2 : 1, kerneloffset=self.kerneloffset, optargs...) where {NK,TF}
    N = self.ncells
    M = size(CPs, 2)

    if size(G, 1)!=M || size(G, 2)!=N
        error("Matrix G with invalid dimensions;"*
              " got $(size(G)), expected ($M, $N).")
    end

    # Build geometric matrix
    derivatives_switch = FastMultipole.DerivativesSwitch(false,true,false) # only velocity
    
    # store old strength and set to unit
    old_strength = copy(self.strength)
    self.strength .= zero(eltype(self.strength))
    if strength_index > 0
        self.strength[:, strength_index] .= 1.0
    else
        for i in 1:NK
            self.strength[:, i] .= 1.0
        end
    end

    Threads.@threads for i_source in 1:N
        for i_target in 1:M
            # get target
            tx, ty, tz = CPs[1, i_target], CPs[2, i_target], CPs[3, i_target]
            target = FastMultipole.StaticArrays.SVector{3,TF}(tx, ty, tz)

            # compute influence
            _, u, _ = induced(target, self, i_source, derivatives_switch; kerneloffset=kerneloffset)

            # update G
            G[i_target, i_source] = u[1] * normals[1, i_target] + u[2] * normals[2, i_target] + u[3] * normals[3, i_target]
        end
    end

    # restore strength
    self.strength .= old_strength
end

function _G_phi!(self::AbstractBody{<:Any,NK,TF}, kernel, G, CPs; strength_index=kernel==ConstantDoublet || kernel==VortexRing && NK>1 ? 2 : 1, kerneloffset=self.kerneloffset, optargs...) where {NK,TF}
    N = self.ncells
    M = size(CPs, 2)

    if size(G, 1)!=M || size(G, 2)!=N
        error("Matrix G with invalid dimensions;"*
              " got $(size(G)), expected ($M, $N).")
    end

    # Build geometric matrix
    derivatives_switch = FastMultipole.DerivativesSwitch(true,false,false) # only potential
    
    # store old strength and set to unit
    old_strength = copy(self.strength)
    self.strength .= zero(eltype(self.strength))
    self.strength[:, strength_index] .= 1.0

    Threads.@threads for i_source in 1:N
        for i_target in 1:M
            # get target
            tx, ty, tz = CPs[1, i_target], CPs[2, i_target], CPs[3, i_target]
            target = FastMultipole.StaticArrays.SVector{3,TF}(tx, ty, tz)

            # compute influence
            phi, _, _ = induced(target, self, i_source, derivatives_switch; kerneloffset=kerneloffset)

            # update G
            G[i_target, i_source] = phi
        end
    end

    # restore strength
    self.strength .= old_strength
end

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

function FastMultipole.buffer_to_target_system!(target_system::NonLiftingBody, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    if PS
        phi = target_buffer[4, i_buffer]
        target_system.potential[i_target] += phi
    end

    if VS
        vx, vy, vz = target_buffer[5, i_buffer], target_buffer[6, i_buffer], target_buffer[7, i_buffer]
        target_system.velocity[1, i_target] += vx
        target_system.velocity[2, i_target] += vy
        target_system.velocity[3, i_target] += vz
    end
end
################################################################################
# COMMON FUNCTIONS
################################################################################
"""
  `generate_loft(bodytype::Type{B}, args...; bodyoptargs=(), optargs...)  where
{B<:NonLiftingBody}`

Generates a lofted non-lifting body. See documentation of
`GeometricTools.generate_loft` for a description of the arguments of this
function.
"""
function generate_loft(bodytype::Type{B}, args...; bodyoptargs=(),
                        dimsplit::Int64=2, optargs...) where {B<:NonLiftingBody}
    # Loft the surface geometry
    grid = gt.generate_loft(args...; optargs...)

    # Split the quadrialateral panels into triangles
    # dimsplit = 2              # Dimension along which to split
    triang_grid = gt.GridTriangleSurface(grid, dimsplit)

    return bodytype(triang_grid; bodyoptargs...)
end

"""
  `generate_revolution_(bodytype::Type{B}, args...; bodyoptargs=(), optargs...)
where {B<:NonLiftingBody}`

Generates a non-lifting body of a body of revolution. See documentation of
`GeometricTools.surface_revolution` for a description of the arguments of this
function.
"""
function generate_revolution(bodytype::Type{B}, args...; bodyoptargs=(),
                             dimsplit::Int64=2, loop_dim::Int64=2,
                             optargs...)  where {B<:NonLiftingBody}
    # Revolve the geometry
    grid = gt.surface_revolution(args...; loop_dim=loop_dim, optargs...)

    # Split the quadrialateral panels into triangles
    # dimsplit = 2              # Dimension along which to split
    triang_grid = gt.GridTriangleSurface(grid, dimsplit)

    return bodytype(triang_grid; bodyoptargs...)
end

##### END OF COMMON FUNTIONS ###################################################


function generate_revolution_liftbody(bodytype::Type{B}, args...; bodyoptargs=(),
                                                  gridprocessing=nothing,
                                                  dimsplit::Int=1,
                                                  # loop_dim::Int=1,
                                                  loop_dim::Int=2,
                                                  axis_angle=270,
                                                  optargs...) where {B<:NonLiftingBody}
    # Revolves the geometry
    grid = gt.surface_revolution(args...; loop_dim=loop_dim,
                                            axis_angle=axis_angle, optargs...)

    # Intermediate processing of grid: rotate to align centerline with x-axis
    if gridprocessing==nothing
        Oaxis = gt.rotation_matrix2(0, 0, 90)
        O = zeros(3)
        gt.lintransform!(grid, Oaxis, O)

    # User-defined intermediate processing of grid
    else
        gridprocessing(grid)
    end

    # Splits the quadrialateral panels into triangles
    # dimsplit = 2              # Dimension along which to split
    triang_grid = gt.GridTriangleSurface(grid, dimsplit)

    # construct body
    return bodytype(triang_grid; bodyoptargs...)
end
