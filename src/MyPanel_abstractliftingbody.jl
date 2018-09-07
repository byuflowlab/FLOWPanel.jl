#=##############################################################################
# DESCRIPTION
    Lifting paneled body types definition.
# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Sep 2018
  * License   : AGPL-3.0
=###############################################################################

################################################################################
# ABSTRACT LIFTING BODY BODY TYPE
################################################################################
"""
  Implementations of AbstractBody are expected to have the following fields:
  * `U::Array{Int64,1}`                 : Indices of all panels along the upper
                                          side of the trailing edge.
  * `L::Array{Int64,1}`                 : Indices of all panels along the lower
                                          side of the trailing edge.
  * `ncellsTE::Int64`                   : Number of cells along trailing edge
  * `nnodesTE::Int64`                   : Number of nodes along trailing edge

  NOTE: U and L are assumed to have the same number of points.

  in addition to the same properties and function expected from an AbstractBody
  implementation with the following functions

  ```julia

  # Solve for distributions of potential flow
  function solve(self::BodyType, Vinfs::Array{Array{T,1},1},
                              D::Array{Array{T,1},1}, args...) where {T<:Real}
    .
    .
    .
  end

  # Outputs a vtk file with the wake
  function _savewake(self::BodyType, filename::String;
                                  len::Rtype=1.0, upper::Bool=true, optargs...)
    .
    .
    .
  end
  ```
"""
abstract type AbstractLiftingBody <: AbstractBody end


# Includes all implementations of AbstractLiftingBody
for header_name in ["liftingbody"]
  include("MyPanel_"*header_name*".jl")
end


# Declares implementations of AbstractLiftingBody
const LBodyTypes = Union{RigidWakeBody}


##### COMMON FUNCTIONS  ########################################################
"""
  `get_TE(self::LBodyTypes, i::Int64; upper::Bool=true)`

Returns the i-th trailing edge point.
"""
function get_TE(self::LBodyTypes, i::Int64; upper::Bool=true)
  if i>self.nnodesTE
    error("Invalid index $i; maximum is $(self.nnodesTE).")
  end

  if upper
    if i==self.nnodesTE
      return gt.get_node(self.grid, gt.get_cell(self.grid, self.U[i-1])[2])
    else
      return gt.get_node(self.grid, gt.get_cell(self.grid, self.U[i])[1])
    end
  else
    if i==self.nnodesTE
      return gt.get_node(self.grid, gt.get_cell(self.grid, self.L[i-1])[1])
    else
      return gt.get_node(self.grid, gt.get_cell(self.grid, self.L[i])[2])
    end
  end
end

"""
  `generate_loft_liftbody(bodytype::Type{LBodyTypes}, args...; optargs...)`
Generates a lofted lifting body of type `bodytype`. See documentation of
`GeometricTools.generate_loft` for a description of the arguments of this
function.
"""
function generate_loft_liftbody(bodytype::Type{LBodyTypes}, args...;
                                                  dimsplit::Int64=2, optargs...)
  # Lofts the surface geometry
  grid = gt.generate_loft(args...; optargs...)

  # Splits the quadrialateral panels into triangles
  # dimsplit = 2              # Dimension along which to split
  triang_grid = gt.GridTriangleSurface(grid, dimsplit)

  ndivs = gt.get_ndivscells(triang_grid)                 # Cells in each dimension
  U = [ sub2ind(ndivs, 1, i) for i in 1:ndivs[2] ]           # Upper LE cells
  L = [ sub2ind(ndivs, ndivs[1], i) for i in 1:ndivs[2] ]    # Lower LE cells

  return bodytype(triang_grid, U, L)
end

"""
  `generate_revolution_liftbody(bodytype::Type{LBodyTypes}, args...; optargs...)`
Generates a lifting body type `bodytype` of a body of revolution. See
documentation of `GeometricTools.surface_revolution` for a description of the
arguments of this function.
"""
function generate_revolution_liftbody(bodytype::Type{LBodyTypes}, args...;
                                                  dimsplit::Int64=2,
                                                  loop_dim::Int64=2, optargs...)
  # Revolves the geometry
  grid = gt.surface_revolution(args...; loop_dim=loop_dim, optargs...)

  # Splits the quadrialateral panels into triangles
  # dimsplit = 2              # Dimension along which to split
  triang_grid = gt.GridTriangleSurface(grid, dimsplit)

  ndivs = gt.get_ndivscells(triang_grid)                 # Cells in each dimension
  U = [ sub2ind(ndivs, 1, i) for i in 1:ndivs[2] ]           # Upper LE cells
  L = [ sub2ind(ndivs, ndivs[1], i) for i in 1:ndivs[2] ]    # Lower LE cells

  return bodytype(triang_grid, U, L)
end

##### COMMON INTERNAL FUNCTIONS  ###############################################
"""
Returns the geometric matrix corresponding to only the vortex rings along a
lifting body impossing the Kutta condition at the trailing edge.
NOTE: This matrix doesn't include the wake.
"""
function _calc_G_lifting_vortexring(grid::gt.GridTriangleSurface,
                                    U::Array{Int64,1}, L::Array{Int64,1})
  return PanelSolver.G_lifting_vortexring(
                  grid.orggrid.nodes,                                  # Nodes
                  [gt.get_cell(grid, i) for i in 1:grid.ncells],       # Panels
                  U, L,                                                # TE
                  [_get_controlpoint(grid, i) for i in 1:grid.ncells], # CPs
                  [gt.get_normal(grid, i) for i in 1:grid.ncells];     # Normals
                )
end

"""
Checks correction definition of trailing edge
"""
function _checkTE(U::Array{Int64,1}, L::Array{Int64,1})
  # Same number of cells
  if size(U)!=size(L)
    return false
  end

  return true
end
##### END OF ABSTRACT LIFTING BODY #############################################
