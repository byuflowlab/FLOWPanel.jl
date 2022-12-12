#=##############################################################################
# DESCRIPTION
    Lifting paneled body types definition.
# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Sep 2018
  * License     : MIT License
=###############################################################################

################################################################################
# ABSTRACT LIFTING BODY BODY TYPE
################################################################################
"""
  Implementations of AbstractLiftingBody are expected to have the following
  fields:
  * `U::Array{Int64,1}`                 : Indices of all panels along the upper
                                          side of the trailing edge.
  * `L::Array{Int64,1}`                 : Indices of all panels along the lower
                                          side of the trailing edge.
  * `ncellsTE::Int64`                   : Number of cells along trailing edge
  * `nnodesTE::Int64`                   : Number of nodes along trailing edge

  NOTE: U and L are assumed to have the same number of points.

  in addition to the same properties and functions expected from an AbstractBody
  implementation. The following functions also need to be implemented.

  ```julia

  # Impose boundary conditions to solve for element strengths
  function solve(self::AbstractLiftingBody, Uinfs::Array{<:Real, 2},
                                            D::Array{<:Real, 2}, args...)
    .
    .
    .
  end

  # Outputs a vtk file with the wake
  function _savewake(self::AbstractLiftingBody, filename::String;
                                  len::Real=1.0, upper::Bool=true, optargs...)
    .
    .
    .
  end
  ```
"""
abstract type AbstractLiftingBody{E, N} <: AbstractBody{E, N} end

"""
    `solve(body::AbstractBody, Uinfs::Array{<:Real, 2})`

Impose boundary conditions to solve for element strengths. `Uinds[:, i]` is the
velocity at the i-th control point used in the boundary condition.

`Das[:, i]` is the unitary vector pointing in the direction that the
semi-infinite wake is shed from the first node of the i-th shedding edge, while
`Dbs[:, i]` is for the second node.
NOTE: These directions are expected to point from the node out to infinite.
"""
function solve(self::AbstractLiftingBody, Uinfs::AbstractMatrix,
               Das::AbstractMatrix, Dbs::AbstractMatrix)
    error("solve(...) for body type $(typeof(self)) has not been implemented yet!")
end

##### COMMON FUNCTIONS  ########################################################
"""
  `generate_loft_liftbody(bodytype::Type{<:AbstractLiftingBody}, args...; optargs...)`
Generates a lofted lifting body of type `bodytype`. See documentation of
`GeometricTools.generate_loft` for a description of the arguments of this
function.
"""
function generate_loft_liftbody(bodytype::Type{B}, args...;
                                bodyoptargs=(), dimsplit::Int=2, optargs...
                               ) where {B<:AbstractLiftingBody}
    # Lofts the surface geometry
    grid = gt.generate_loft(args...; optargs...)

    # Splits the quadrialateral panels into triangles
    # dimsplit = 2              # Dimension along which to split
    triang_grid = gt.GridTriangleSurface(grid, dimsplit)

    ndivs = gt.get_ndivscells(triang_grid)              # Cells in each dimension
    U = [ Base._sub2ind(ndivs, ndivs[1]-1, i) for i in 1:ndivs[2] ] # Upper LE cells
    L = [ Base._sub2ind(ndivs, 2, i) for i in 1:ndivs[2] ]          # Lower LE cells

    nedges = length(U)
    shedding = zeros(Int, 6, nedges)
    for (ei, (u, l)) in enumerate(zip(U, L))
        shedding[1, ei] = u
        shedding[2, ei] = 3
        shedding[3, ei] = 2

        shedding[4, ei] = l
        shedding[5, ei] = 3
        shedding[6, ei] = 2
    end

    return bodytype(triang_grid, shedding; bodyoptargs...)
end

"""
  `generate_revolution_liftbody(bodytype::Type{<:AbstractLiftingBody}, args...; optargs...)`
Generates a lifting body type `bodytype` of a body of revolution. See
documentation of `GeometricTools.surface_revolution` for a description of the
arguments of this function.
"""
function generate_revolution_liftbody(bodytype::Type{B}, args...;
                                                  bodyoptargs=(),
                                                  dimsplit::Int=2,
                                                  loop_dim::Int=2, optargs...
                                      ) where {B<:AbstractLiftingBody}
    # Revolves the geometry
    grid = gt.surface_revolution(args...; loop_dim=loop_dim, optargs...)

    # Splits the quadrialateral panels into triangles
    # dimsplit = 2              # Dimension along which to split
    triang_grid = gt.GridTriangleSurface(grid, dimsplit)

    ndivs = gt.get_ndivscells(triang_grid)                 # Cells in each dimension
    U = [ Base._sub2ind(ndivs, ndivs[1]-1, i) for i in 1:ndivs[2] ] # Upper LE cells
    L = [ Base._sub2ind(ndivs, 2, i) for i in 1:ndivs[2] ]          # Lower LE cells

    nedges = length(U)
    shedding = zeros(Int, 6, nedges)
    for (ei, (u, l)) in enumerate(zip(U, L))
        shedding[1, ei] = u
        shedding[2, ei] = 3
        shedding[3, ei] = 2

        shedding[4, ei] = l
        shedding[5, ei] = 3
        shedding[6, ei] = 2
    end

    return bodytype(triang_grid, shedding; bodyoptargs...)
end

##### COMMON INTERNAL FUNCTIONS  ###############################################
"""
Checks correction definition of trailing edge
"""
function _checkTE(grid, shedding::Array{Int, 2}; tol=1e1*eps())

    nodes = grid.orggrid.nodes

    # Correct number of inputs
    if size(shedding, 1) != 6
        return false
    end

    tricoor, quadcoor, lin, ndivscells, cin = gt.generate_getcellt_args(grid)

    # Check that node position along edge of each side are coincident
    for (pi, nia, nib, pj, nja, njb) in eachcol(shedding)

        if pj != -1

            # Convert node indices from panel-local to global
            pia = gt.get_cell_t(tricoor, quadcoor, grid, pi, nia, lin, ndivscells, cin)
            pib = gt.get_cell_t(tricoor, quadcoor, grid, pi, nib, lin, ndivscells, cin)
            pja = gt.get_cell_t(tricoor, quadcoor, grid, pj, nia, lin, ndivscells, cin)
            pjb = gt.get_cell_t(tricoor, quadcoor, grid, pj, nib, lin, ndivscells, cin)

            for i in 1:3
                if abs(nodes[i, pia] - nodes[i, pjb]) > tol
                    println("rabbit2")
                    return false
                elseif abs(nodes[i, pib] - nodes[i, pja]) > tol
                    println("rabbit3")
                    return false
                end
            end

        end

    end

    return true
end
##### END OF ABSTRACT LIFTING BODY #############################################
