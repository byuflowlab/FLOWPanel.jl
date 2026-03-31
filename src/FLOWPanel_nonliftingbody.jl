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
struct NonLiftingBody{E, N, TF} <: AbstractBody{E, N, TF}

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
    Uinf::Matrix{TF}
    U::Matrix{TF}
    phi::Vector{TF}
    Cp::Vector{TF}
    Cps::Vector{TF}
    Gamma::Vector{TF}
    F::Matrix{TF}
    solved::Bool

    # Internal variables
    strength::Array{TF, 2}              # strength[i,j] is the stength of the i-th panel with the j-th element type
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

function NonLiftingBody{E, N, TF}(
                nodes::Matrix{TF}, cells::Matrix{Int};
                vtk_cells::Vector{<:WriteVTK.MeshCell}=[WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, cells[:, i]) for i in 1:size(cells, 2)],
                neighbor::Matrix{Int}=zeros(Int, 3, size(cells, 2)),
                nnodes=size(nodes, 2), ncells=size(cells, 2),
                Oaxis=Array{TF,2}(1.0I, 3, 3), O=zeros(TF,3),
                Uinf=zeros(TF, 0, 0),
                U=zeros(TF, 0, 0),
                phi=zeros(TF, 0),
                Cp=zeros(TF, 0),
                Cps=zeros(TF, size(cells, 2)),
                Gamma=zeros(TF, 0),
                F=zeros(TF, 0, 0),
                solved=false,
                strength=zeros(size(cells, 2), N),
                velocity=zeros(3, size(cells, 2)),
                CPoffset=1e-14,
                kerneloffset=1e-8,
                kernelcutoff=1e-14,
                characteristiclength=characteristiclength_unitary,
                check_mesh=true, watertight=false,
                inside_offset=1e-6
              ) where {E, N, TF}
    # check if mesh is watertight
    if check_mesh
        # Need to implement an equivalent check for standalone nodes/cells here, or simply skip.
        # Alternatively, we just depend on the watertight flag explicitly given to the constructor.
    end

    # If the user did not provide half-edges via neighbor matrix, compute them:
    # We will assume that they are either pre-computed via the other constructor or zero unless calculated.
    # Note: `geometrictools` has functions for finding neighbors, but here we require a manual pass or rely on the grid-based constructor extracting it first.
    # If a pure nodes/cells constructor is used, we'll assume the user computes their own neighbors and passes them,
    # or that a boundary defaults to zero.

    return NonLiftingBody{E, N, TF}(
                nodes, vtk_cells, neighbor,
                nnodes, ncells, cells,
                Oaxis, O,
                Uinf, U, phi, Cp, Cps, Gamma, F, solved,
                strength,
                velocity,
                CPoffset,
                kerneloffset,
                kernelcutoff,
                characteristiclength,
                watertight,
                inside_offset
              )
end

function NonLiftingBody{E, N, TF}(
                grid::gt.GridTriangleSurface;
                optargs...
              ) where {E, N, TF}
    
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
    
    return NonLiftingBody{E, N, TF}(
                nodes, cells;
                vtk_cells=vtk_cells, neighbor=neighbor, watertight=watertight_guess, optargs...
              )
end

function (NonLiftingBody{E})(grid::gt.GridTriangleSurface; optargs...) where {E}
    return NonLiftingBody{E, _count(E), eltype(grid._nodes)}(grid; optargs...)
end

function (NonLiftingBody{E})(nodes::Matrix{TF}, cells::Matrix{Int}; optargs...) where {E, TF}
    return NonLiftingBody{E, _count(E), TF}(nodes, cells; optargs...)
end

function save(body::NonLiftingBody, args...; optargs...)
    return save_base(body, args...; optargs...)
end

calc_elprescribe(::NonLiftingBody{ConstantSource, 1}) = Tuple{Int,Float64}[]
calc_elprescribe(body::NonLiftingBody{VortexRing, 1}) = body.watertight ? [(1, 0.0)] : Tuple{Int,Float64}[]
calc_elprescribe(body::NonLiftingBody{ConstantDoublet, 1}) = body.watertight ? [(1, 0.0)] : Tuple{Int,Float64}[]

solved_field_name(::NonLiftingBody{ConstantSource, 1}) = "sigma"
solved_field_name(::NonLiftingBody{ConstantDoublet, 1}) = "mu"
solved_field_name(::NonLiftingBody{VortexRing, 1}) = "gamma"
solved_field_name(::NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2}) = "mu"

#### END OF NON-LIFTING BODY  ##################################################


################################################################################
# CONSTANT-SOURCE SOLVER
################################################################################

function _G_U!(self::AbstractBody{<:Any,NK,TF}, kernel, G, CPs, normals, backend::AbstractBackend=DirectBackend(); strength_index=kernel==ConstantDoublet && NK>1 ? 2 : 1, kerneloffset=self.kerneloffset, optargs...) where {NK,TF}
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

function _G_phi!(self::AbstractBody{<:Any,NK,TF}, kernel, G, CPs, backend::AbstractBackend=DirectBackend(); strength_index=kernel==ConstantDoublet && NK>1 ? 2 : 1, kerneloffset=self.kerneloffset, optargs...) where {NK,TF}
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
            phi, _, _ = induced(target, self, i_source, derivatives_switch; kerneloffset=kerneloffset)

            # update G
            G[i_target, i_source] = phi
        end
    end

    # restore strength
    self.strength .= old_strength
end

function _G_phi!(bodies::Tuple, kernel, G, CPs, backend::AbstractBackend=DirectBackend(); 
                 strength_index=nothing, kerneloffset=nothing, optargs...)
    # Track starting column index in G for each body
    col_start = 1

    for body in bodies
        N = body.ncells
        col_end = col_start + N - 1

        # slice G for this body's columns
        G_body = view(G, :, col_start:col_end)

        # Determine strength_index for this body if not passed
        si = strength_index
        if si === nothing
            si = kernel == ConstantDoublet && body.NK > 1 ? 2 : 1
        end

        # Determine kerneloffset for this body if not passed
        ko = kerneloffset
        if ko === nothing
            ko = getproperty(body, :kerneloffset)
        end

        # Call the single-body version on the slice of G
        _G_phi!(body, kernel, G_body, CPs, backend; 
                strength_index=si, kerneloffset=ko, optargs...)

        # Update column start for the next body
        col_start = col_end + 1
    end
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

##### END OF FASTMULTIPOLE BACKEND SUPPORT #####################################

################################################################################
# ABSTRACT SOLVER INTERFACE
################################################################################

function solve2!(self::NonLiftingBody{TK,1,TFG}, Uinfs::Array{TFS, 2}, solver::AbstractMatrixfulSolver{false};
        backend=DirectBackend(),
        update_G::Bool=false,   # Whether to update the influence matrix G
        strength_name=get_strength_name(self),  # Name of the strength field to solve for
        optargs...              # Additional optional arguments to _G_U!
    ) where {TK, TFG, TFS<:Real}

    # formatting assertions
    @assert size(self.strength, 2) == 1 "AbstractBody{<:Any, 1} expected to have single strength per element; got size $(size(self.strength, 2))."

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)
    CPs_inside = _calc_controlpoints(self, normals; off=-1e-10)
    
    # update influence matrix (if requested)
    if update_G

        # Update geometric matrix (left-hand-side influence matrix)
        _G_U!(self, TK, solver.G, CPs, normals, backend; optargs...)
    end

    # generate RHS
    TF = promote_type(TFG, TFS)
    RHS = zeros(TF, self.ncells)
    calc_bc_noflowthrough!(RHS, Uinfs, normals)

    # Solve system of equations
    solution = zeros(TF, self.ncells)
    solve_matrix!(solution, solver.G, RHS, solver)
    
    # set solved flag
    _solvedflag(self, true)
    
    # Assign solution to body element strengths
    # _assign_elementstrengths!(self, solution)
    self.strength .= solution

    # verify that source strength is equal to -∂ϕ/∂n on the boundary
    us_outside = zeros(TF, 3, self.ncells)
    _Uind!(self, CPs, us_outside, backend; optargs...)
    us_inside = zeros(TF, 3, self.ncells)
    _Uind!(self, CPs_inside, us_inside, backend; optargs...)
    # add_field(self, "us_inside", "vector", collect(eachcol(us_inside)), "cell")
    # add_field(self, "delta_u_normal", "scalar", vec(sum((us_outside - us_inside) .* normals, dims=1)), "cell")

    # save solution fields
    # add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    # add_field(self, strength_name, "scalar", view(self.strength, :, 1), "cell")
    # add_field(self, "normals", "vector", collect(eachcol(normals)), "cell")

    return nothing
end

function FastMultipole.value_to_strength!(source_buffer, ::NonLiftingBody, i_body, value, rlx)
    prev_value = source_buffer[5, i_body]
    source_buffer[5, i_body] = rlx*value + (1.0 - rlx)*prev_value
end

function FastMultipole.value_to_strength!(source_buffer, ::NonLiftingBody, i_body, value)
    source_buffer[5, i_body] = value
end

function FastMultipole.target_influence_to_buffer!(target_buffer, i_buffer, derivatives_switch, target_system::NonLiftingBody, i_target)
    vx, vy, vz = target_system.velocity[1, i_target], target_system.velocity[2, i_target], target_system.velocity[3, i_target]
    target_buffer[5, i_buffer] = vx
    target_buffer[6, i_buffer] = vy
    target_buffer[7, i_buffer] = vz
end

function FastMultipole.buffer_to_target_system!(target_system::NonLiftingBody, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    vx, vy, vz = target_buffer[5, i_buffer], target_buffer[6, i_buffer], target_buffer[7, i_buffer]
    target_system.velocity[1, i_target] = vx
    target_system.velocity[2, i_target] = vy
    target_system.velocity[3, i_target] = vz
end

##### END OF ABSTRACT SOLVER INTERFACE ##########################################


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
