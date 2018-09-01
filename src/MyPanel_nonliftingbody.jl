#=##############################################################################
# DESCRIPTION
    Non-lifting paneled body types definition.
# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jun 2018
  * License   : AGPL-3.0
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
immutable NonLiftingBody <: AbstractBody

  # User inputs
  grid::gt.GridTriangleSurface              # Paneled geometry

  # Properties
  nnodes::Int64                             # Number of nodes
  ncells::Int64                             # Number of cells
  fields::Array{String, 1}                  # Available fields (solutions)
  Oaxis::Array{T,2} where {T<:Real}         # Coordinate system of original grid
  O::Array{T,1} where {T<:Real}             # Position of CS of original grid

  # Internal variables
  _G::Array{T,2} where {T<:Real}            # Geometric solution matrix

  NonLiftingBody( grid,
                  nnodes=grid.nnodes, ncells=grid.ncells,
                    fields=Array{String,1}(),
                    Oaxis=eye(3), O=zeros(3),
                  _G=_calc_Gsource(grid)
         ) = new( grid,
                  nnodes, ncells,
                    fields,
                    Oaxis, O,
                  _G
         )
end


function solve(self::NonLiftingBody, Vinfs::Array{Array{T,1},1}) where {T<:Real}
  if size(Vinfs,1) != self.ncells
    error("Invalid Vinfs; expected size $(self.ncells), got $(size(Vinfs,1))")
  end

  lambda = [-dot(Vinfs[i], get_normal(self, i)) for i in 1:self.ncells]
  sigma = self._G\lambda

  add_field(self, "Vinf", Vinfs)
  add_field(self, "sigma", sigma)
  _solvedflag(self, true)
end




##### INTERNAL FUNCTIONS  ######################################################
function _calc_Gsource(grid::gt.GridTriangleSurface)
  return PanelSolver.G_constant_source(
                  grid.orggrid.nodes,                                  # Nodes
                  [gt.get_cell(grid, i) for i in 1:grid.ncells],       # Panels
                  [_get_controlpoint(grid, i) for i in 1:grid.ncells], # CPs
                  [gt.get_normal(grid, i) for i in 1:grid.ncells],     # Normals
                )
end

"""
Returns the velocity induced by the body on the targets `targets`. It adds the
velocity at the i-th target to out[i].
"""
function _Vind(self::NonLiftingBody, targets::Array{Array{T,1},1},
                                      out::Array{Array{T,1},1}) where{T<:Real}
  # Iterates over panels
  for i in 1:self.ncells
    # Velocity of i-th  panel on every target
    PanelSolver.Vconstant_source(
                    gt.get_cellnodes(self.grid, i),    # Nodes in i-th panel
                    get_fieldval(self, "sigma", i; _check=false),  # Strength
                    targets,                           # Targets
                    out;                               # Outputs
                  )
  end
end
##### END OF NON-LIFTING BODY ##################################################



################################################################################
# NON-LIFTING VORTEX-RING BODY TYPE
################################################################################
"""
  `NonLiftingBodyVRing(grid::gt.GridTriangleSurface)`

Non-lifting paneled body that is solved using vortex ring panels.
`grid` is the grid surface (paneled geometry).

  **Properties**
  * `nnodes::Int64`                     : Number of nodes
  * `ncells::Int64`                     : Number of cells
  * `fields::Array{String, 1}`          : Available fields (solutions)
  * `Oaxis::Array{T,2} where {T<:Real}` : Coordinate system of original grid
  * `O::Array{T,1} where {T<:Real}`     : Position of CS of original grid

"""
immutable NonLiftingBodyVRing <: AbstractBody

  # User inputs
  grid::gt.GridTriangleSurface              # Paneled geometry

  # Properties
  nnodes::Int64                             # Number of nodes
  ncells::Int64                             # Number of cells
  fields::Array{String, 1}                  # Available fields (solutions)
  Oaxis::Array{T,2} where {T<:Real}         # Coordinate system of original grid
  O::Array{T,1} where {T<:Real}             # Position of CS of original grid

  # Internal variables
  _G::Array{T,2} where {T<:Real}            # Geometric solution matrix

  NonLiftingBodyVRing( grid,
                  nnodes=grid.nnodes, ncells=grid.ncells,
                    fields=Array{String,1}(),
                    Oaxis=eye(3), O=zeros(3),
                  _G=_calc_Gvring(grid)
         ) = new( grid,
                  nnodes, ncells,
                    fields,
                    Oaxis, O,
                  _G
         )
end


function solve(self::NonLiftingBodyVRing, Vinfs::Array{Array{T,1},1}
                                                              ) where {T<:Real}
  if size(Vinfs,1) != self.ncells
    error("Invalid Vinfs; expected size $(self.ncells), got $(size(Vinfs,1))")
  end

  lambda = [-dot(Vinfs[i], get_normal(self, i)) for i in 1:self.ncells]
  Gamma = self._G\lambda

  add_field(self, "Vinf", Vinfs)
  add_field(self, "Gamma", Gamma)
  _solvedflag(self, true)
end



##### INTERNAL FUNCTIONS  ######################################################
function _calc_Gvring(grid::gt.GridTriangleSurface)
  return PanelSolver.G_vortexring(
                  grid.orggrid.nodes,                                  # Nodes
                  [gt.get_cell(grid, i) for i in 1:grid.ncells],       # Panels
                  [_get_controlpoint(grid, i) for i in 1:grid.ncells], # CPs
                  [gt.get_normal(grid, i) for i in 1:grid.ncells],     # Normals
                )
end

"""
Returns the velocity induced by the body on the targets `targets`. It adds the
velocity at the i-th target to out[i].
"""
function _Vind(self::NonLiftingBodyVRing, targets::Array{Array{T,1},1},
                                      out::Array{Array{T,1},1}) where{T<:Real}
  # Iterates over panels
  for i in 1:self.ncells
    # Velocity of i-th  panel on every target
    PanelSolver.Vvortexring(
                    gt.get_cellnodes(self.grid, i),    # Nodes in i-th panel
                    get_fieldval(self, "Gamma", i; _check=false),  # Strength
                    targets,                           # Targets
                    out;                               # Outputs
                  )
  end
end
##### END OF NON-LIFTING VORTEX-RING BODY ######################################





################################################################################
# COMMON FUNCTIONS
################################################################################
"""
  `generate_loft_nonliftbody(args...; optargs...)`
Generates a lofted non-lifting body. See documentation of
`GeometricTools.generate_loft` for a description of the arguments of this
function.
"""
function generate_loft_nonliftbody(args...; vortexring=false, dimsplit::Int64=2,
                                                                    optargs...)
  # Lofts the surface geometry
  grid = gt.generate_loft(args...; optargs...)

  # Splits the quadrialateral panels into triangles
  # dimsplit = 2              # Dimension along which to split
  triang_grid = gt.GridTriangleSurface(grid, dimsplit)

  if vortexring
    return NonLiftingBodyVRing(triang_grid)
  else
    return NonLiftingBody(triang_grid)
  end
end

"""
  `generate_revolution_nonliftbody(args...; optargs...)`
Generates a non-lifting body of a body of revolution. See documentation of
`GeometricTools.surface_revolution` for a description of the arguments of this
function.
"""
function generate_revolution_nonliftbody(args...; vortexring=false,
                                                  dimsplit::Int64=2,
                                                  loop_dim::Int64=2, optargs...)
  # Revolves the geometry
  grid = gt.surface_revolution(args...; loop_dim=loop_dim, optargs...)

  # Splits the quadrialateral panels into triangles
  # dimsplit = 2              # Dimension along which to split
  triang_grid = gt.GridTriangleSurface(grid, dimsplit)

  if vortexring
    return NonLiftingBodyVRing(triang_grid)
  else
    return NonLiftingBody(triang_grid)
  end
end
##### END OF COMMON FUNTIONS ###################################################
