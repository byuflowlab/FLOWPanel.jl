#=##############################################################################
# DESCRIPTION
    Non-lifting paneled body type definition.
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
                  _G=_calc_G(grid)
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

  add_field(self, "sigma", sigma)
end


"""
  `get_field(self::NonLiftingBody, field_name::String)`

Returns the requested field.
"""
function get_field(self::NonLiftingBody, field_name::String)
  if !(field_name in self.fields)
    error("Field $field_name not found! Available fields: $(self.fields)")
  end

  return self.grid.field[field_name]
end

"""
  `get_field(self::NonLiftingBody, field_name::String, i::Int64)`

Returns the requested field value. Give it `_check=false` to skip checking logic
for faster operation.
"""
function get_field_val(self::NonLiftingBody, field_name::String, i::Int64;
                      _check::Bool=true)
  if _check
    if i<1
      error("Invalid index $i.")
    elseif FIELDS[field_name]["entry_type"]=="node" && i>self.nnodes
      error("Invalid index $i. Maximum is $(self.nnodes).")
    elseif FIELDS[field_name]["entry_type"]=="cell" && i>self.ncells
      error("Invalid index $i. Maximum is $(self.ncells).")
    end
  end

  return self.grid.field[field_name]["field_data"][i]
end


"""
  `add_field(self::NonLiftingBody, field_name::String, field_data)`

Adds a new field to the body.
"""
function add_field(self::NonLiftingBody, field_name::String, field_data)

  # ERROR CASES
  if !(field_name in keys(FIELDS))
    error("Invalid field $field_name. Implemeted fields are: $(keys(FIELDS))")
  end

  gt.add_field(self.grid, field_name, FIELDS[field_name]["field_type"],
                field_data, FIELDS[field_name]["entry_type"])

  if !(field_name in self.fields)
    push!(self.fields, field_name)
  end

  nothing
end


"""
  `generate_loft_nonliftbody(args...; optargs...)`
Generates a lofted non-lifting body. See documentation of
`GeometricTools.generate_loft` for a description of the arguments of this
function.
"""
function generate_loft_nonliftbody(args...; dimsplit::Int64=2, optargs...)
  # Lofts the surface geometry
  grid = gt.generate_loft(args...; optargs...)

  # Splits the quadrialateral panels into triangles
  # dimsplit = 2              # Dimension along which to split
  triang_grid = gt.GridTriangleSurface(grid, dimsplit)

  return NonLiftingBody(triang_grid)
end

"""
  `generate_revolution_nonliftbody(args...; optargs...)`
Generates a non-lifting body of a body of revolution. See documentation of
`GeometricTools.surface_revolution` for a description of the arguments of this
function.
"""
function generate_revolution_nonliftbody(args...; dimsplit::Int64=2,
                                                loop_dim::Int64=2, optargs...)
  # Revolves the geometry
  grid = gt.surface_revolution(args...; loop_dim=loop_dim, optargs...)

  # Splits the quadrialateral panels into triangles
  # dimsplit = 2              # Dimension along which to split
  triang_grid = gt.GridTriangleSurface(grid, dimsplit)

  return NonLiftingBody(triang_grid)
end

##### INTERNAL FUNCTIONS  ######################################################
function _calc_G(grid::gt.GridTriangleSurface)
  return PanelSolver.G_constant_source(
                  grid.orggrid.nodes,                                  # Nodes
                  [gt.get_cell(grid, i) for i in 1:grid.ncells],       # Panels
                  [_get_controlpoint(grid, i) for i in 1:grid.ncells], # CPs
                  [gt.get_normal(grid, i) for i in 1:grid.ncells],     # Normals
                )
end

##### END OF NON-LIFTING BODY ##################################################
