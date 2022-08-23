#=##############################################################################
# DESCRIPTION
    Non-lifting paneled body types definition.
# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jun 2018
  * License   : MIT License
=###############################################################################


################################################################################
# NON-LIFTING BODY TYPE
################################################################################

"""
  `NonLiftingBody(grid::gt.GridTriangleSurface)`

Non-lifting paneled body that is solved using a constant source distribution.
`grid` is the grid surface (paneled geometry).

  **Properties**
  * `nnodes::Int64`                     : Number of nodes
  * `ncells::Int64`                     : Number of cells
  * `fields::Array{String, 1}`          : Available fields (solutions)
  * `Oaxis::Array{T,2} where {T<:Real}` : Coordinate system of original grid
  * `O::Array{T,1} where {T<:Real}`     : Position of CS of original grid

"""
struct NonLiftingBody{E, N} <: AbstractBody{E, N}

    # User inputs
    grid::gt.GridTriangleSurface              # Paneled geometry

    # Properties
    nnodes::Int64                             # Number of nodes
    ncells::Int64                             # Number of cells
    fields::Array{String, 1}                  # Available fields (solutions)
    Oaxis::Array{<:Number,2}                  # Coordinate system of original grid
    O::Array{<:Number,1}                      # Position of CS of original grid

    # Internal variables
    strength::Array{<:Number, 2}             # strength[i,j] is the stength of the i-th panel with the j-th element type
    CPoffset::Float64                        # Control point offset in normal direction

    NonLiftingBody{E, N}(
                    grid;
                    nnodes=grid.nnodes, ncells=grid.ncells,
                    fields=Array{String,1}(),
                    Oaxis=Array(1.0I, 3, 3), O=zeros(3),
                    strength=zeros(grid.ncells, N),
                    CPoffset=0.005
                  ) where {E, N} = new(
                    grid,
                    nnodes, ncells,
                    fields,
                    Oaxis, O,
                    strength,
                    CPoffset
                  )
end


function (NonLiftingBody{E})(args...; optargs...) where {E}
    return NonLiftingBody{E, _count(E)}(args...; optargs...)
end

function _solve(self::NonLiftingBody{ConstantSource, 1},
                                            normals, CPs, G, Vinfs)

    if size(Vinfs, 2) != self.ncells
        error("Invalid Vinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Vinfs))")
    end

    # Define right-hand side
    lambda = [-dot(Vinf, normal) for (Vinf, normal) in
                                        zip(eachcol(Vinfs), eachcol(normals))]

    # Solve the system of equations
    sigma = G\lambda

    add_field(self, "Vinf", "vector", eachcol(Vinfs), "cell")
    add_field(self, "sigma", "scalar", sigma, "cell")
    _solvedflag(self, true)
end

function solve(self::NonLiftingBody{ConstantSource, 1},
                Vinfs::Arr1) where {T1, Arr1<:AbstractArray{T1, 2}}

    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals; off=self.CPoffset)


    # G = zeros(self.ncells, self.ncells)

    return _solve(self, normals, CPs, G, Vinfs)
end


##### INTERNAL FUNCTIONS  ######################################################
# """
# Returns the velocity induced by the body on the targets `targets`. It adds the
# velocity at the i-th target to out[i].
# """
# function _Uind(self::NonLiftingBody, targets::Array{Array{T1,1},1},
#                           out::Array{Array{T2,1},1}) where{T1<:RType, T2<:RType}
#   # Iterates over panels
#   for i in 1:self.ncells
#     # Velocity of i-th  panel on every target
#     U_constant_source(
#                     gt.get_cellnodes(self.grid, i),    # Nodes in i-th panel
#                     get_fieldval(self, "sigma", i; _check=false),  # Strength
#                     targets,                           # Targets
#                     out;                               # Outputs
#                   )
#   end
# end


"""
Computes the geometric matrix (left-hand side matrix of the system of equation)
and stores it under `G`.

**ARGUMENTS**
  * `G::Array{T,2}`                     : Pre-allocated output memory.
  * `CPs::Array{T,2}`                   : Control points.
  * `normals::Array{T,2}`               : Normal associated to every CP.
"""
function _G_U!(self::NonLiftingBody{ConstantSource, 1},
                    G::Arr1, CPs::Arr2, normals::Arr3;
                    optargs...
               ) where{ T1, Arr1<:AbstractArray{T1, 2},
                        T2, Arr2<:AbstractArray{T2, 2},
                        T3, Arr3<:AbstractArray{T3, 2}}

    N = self.ncells

    if size(G, 1)!=size(G, 2) || size(G, 1)!=N
        error("Matrix G with invalid dimension;"*
              " got $(size(G)), expected ($N, $N).")
    end

    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Build geometric matrix
    for (pj, Gslice) in enumerate(eachcol(G))

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pj, lin, ndivscells, cin)

        U_constant_source(
                          self.grid.orggrid.nodes,          # All nodes
                          panel,                             # Index of nodes that make this panel
                          1.0,                               # Unitary strength
                          CPs,                               # Targets
                          # view(G, :, pj);                  # Velocity of j-th panel on every CP
                          Gslice;
                          dot_with=normals,                  # Normal of every CP
                          optargs...
                         )
    end
end
