#=##############################################################################
# DESCRIPTION
    Abstract body type definition.
# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jun 2018
  * License   : AGPL-3.0
=###############################################################################


################################################################################
# ABSTRACT BODY TYPE
################################################################################
"""
  Implementations of AbstractBody are expected to have the following fields.
  * `grid::gt.GridTriangleSurface `     : Paneled geometry
  * `nnodes::Int64`                     : Number of nodes
  * `ncells::Int64`                     : Number of cells
  * `fields::Array{String, 1}`          : Available fields (solutions)
  * `Oaxis::Array{T,2} where {T<:Real}` : Coordinate system of original grid
  * `O::Array{T,1} where {T<:Real}`     : Position of CS of original grid

  and the following functions

  ```julia

  # Solve for distributions of potential flow
  function solve(self::BodyType, Vinfs::Array{Array{T,1},1}, args...
                                                              ) where {T<:Real}
    .
    .
    .
  end

  # Returns the velocity induced by the body on the targets `targets`. It adds
  # the velocity at the i-th target to out[i].
  function _Vind(self::BodyTypes, targets::Array{Array{T,1},1},
                    out::Array{Array{T,1},1}, args...; optargs...) where{T<:Real}
    .
    .
    .
  end
  ```
"""
abstract type AbstractBody end

# Includes all implementations of AbstractGrid
for header_name in ["nonliftingbody"]
  include("MyPanel_"*header_name*".jl")
end

# Declares implementations of AbstractGrid
const BodyTypes = Union{NonLiftingBody, NonLiftingBodyVRing}




##### COMMON FUNCTIONS  ########################################################

"""
  `Vind(self::BodyTypes, targets::Array{Array{T,1},1},
                  out::Array{Array{T,1},1}, args...; optargs...) where{T<:Real}`

Returns the velocity induced by the body on the targets `targets`. It adds the
velocity at the i-th target to out[i].
"""
function Vind(self::BodyTypes, targets::Array{Array{T,1},1},
                  out::Array{Array{T,1},1}, args...; optargs...) where{T<:Real}

  # ERROR CASES
  if check_solved(self)==false
    error("Body hasn't been solved yet. Please call `solve()` function first.")
  end

  _Vind(self, targets, out, args...; optargs...)
end

"""
  `save(body::BodyTypes, filename::String; opt_args...)`

Outputs a vtk file of this body. See GeometricTools.save(grid, ...) for a
descrition of optional arguments `opt_args...`.
"""
function save(body::BodyTypes, filename::String; out_cellindex::Bool=false,
                                                 out_nodeindex::Bool=false,
                                                 out_controlpoints::Bool=false,
                                                                    opt_args...)
  # Add special fields
  if out_cellindex
    gt.add_field(body.grid, "cellindex", "scalar",
                    [i for i in 1:body.ncells], "cell")
  end
  if out_nodeindex
    gt.add_field(body.grid, "nodeindex", "scalar",
                    [i for i in 1:body.nnodes], "node")
  end

  # Outputs control points
  if out_controlpoints
    save_controlpoints(body, filename; opt_args...)
  end

  gt.save(body.grid, filename; opt_args...)
end

"""
  `save_controlpoints(body::BodyTypes, filename::String;
suffix::String="_cp", opt_args...)`

Outputs a vtk file with the control points of the body along with associated
normals and surface velocities.
"""
function save_controlpoints(body::BodyTypes, filename::String;
                                              suffix::String="_cp", opt_args...)
    # Control points
    CPs = [get_controlpoint(body, i) for i in 1:body.ncells]

    # Normals
    data = [
            Dict( "field_name"  => "normal",
                  "field_type"  => "vector",
                  "field_data"  => [get_normal(body, i) for i in 1:size(CPs,1)])
            ]

    # Surface velocity
    if check_solved(body)
      Vsurf = [get_fieldval(body, "Vinf", i) for i in 1:size(CPs,1)]
      Vind(body, CPs, Vsurf)

      push!(data,
              Dict( "field_name"  => "V",
                    "field_type"  => "vector",
                    "field_data"  => Vsurf)
      )
    end

    # Generates vtk
    gt.generateVTK(filename*suffix, CPs; point_data=data, opt_args...)
end

"""
  `get_controlpoint(body::BodyTypes, i::Int64 or coor::Array{Int64,1})`

Returns the control point on the i-th panel.
"""
function get_controlpoint(body::BodyTypes, args...)
  return _get_controlpoint(body.grid, args...)
end

"""
  `get_ndivscells(body::BodyTypes)`

Returns a tuple with the number of cells in each parametric dimension
"""
get_ndivscells(body::BodyTypes) = body.grid._ndivscells

"""
  `get_ndivsnodes(body::BodyTypes)`

Returns a tuple with the number of nodes in each parametric dimension
"""
get_ndivsnodes(body::BodyTypes) = body.grid._ndivsnodes

"""
  `get_unitvectors(body::BodyTypes, i::Int64 or coor::Array{Int64,1})`

Returns the unit vectors t,n,o of the i-th panel, with t the tanget vector,
n normal, and o oblique.
"""
get_unitvectors(body::BodyTypes, args...) = gt.get_unitvectors(body.grid, args...)

"""
  `get_normal(body::BodyTypes, i::Int64 or coor::Array{Int64,1})`

Returns the normal vector the i-th panel.
"""
get_normal(body::BodyTypes, args...) = gt.get_normal(body.grid, args...)

"""
  `get_field(self::BodyTypes, field_name::String)`

Returns the requested field.
"""
function get_field(self::BodyTypes, field_name::String)
  if !(field_name in self.fields)
    error("Field $field_name not found! Available fields: $(self.fields)")
  end

  return self.grid.field[field_name]
end

"""
  `get_field(self::BodyTypes, field_name::String, i::Int64)`

Returns the requested field value. Give it `_check=false` to skip checking logic
for faster operation.
"""
function get_fieldval(self::BodyTypes, field_name::String, i::Int64;
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
  `add_field(self::BodyTypes, field_name::String, field_data)`

Adds a new field to the body. It overwrites the field if it already existed.
"""
function add_field(self::BodyTypes, field_name::String, field_data)

  # ERROR CASES
  if !(field_name in keys(FIELDS))
    error("Invalid field $field_name. Implemented fields are: $(keys(FIELDS))")
  end

  # Adds field to grid
  gt.add_field(self.grid, field_name, FIELDS[field_name]["field_type"],
                field_data, FIELDS[field_name]["entry_type"]; raise_warn=false)

  # Registers the field
  if !(field_name in self.fields)
    push!(self.fields, field_name)
  end

  nothing
end

"""
  `check_field(self::BodyTypes, field_name::String)`

Returns `true` of the body has the field `field_name`. Returns false otherwise.
"""
check_field(self::BodyTypes, field_name::String) = field_name in self.fields


"""
  `check_solved(self::BodyTypes)`

Returns `true` of the body has been solved. Returns false otherwise.
"""
function check_solved(self::BodyTypes)
  if check_field(self, "solved")
    return get_fieldval(self, "solved", 1)
  else
    return false
  end
end


"""
  `rotate(body::BodyTypes, roll::Real, pitch::Real, yaw::Real;
translation::Array{T, 1}=zeros(3), reset_fields::Bool=true
) where{T<:Real}`

Rotates and translates the body by the given axial angles.

NOTE: Naming follows aircraft convention, with
* roll:   rotation about x-axis.
* pitch:  rotation about y-axis.
* yaw:    rotation about z-axis.
"""
function rotate(body::BodyTypes, roll::Real, pitch::Real, yaw::Real;
                  translation::Array{T, 1}=zeros(3), reset_fields::Bool=true
                ) where{T<:Real}

  M = gt.rotation_matrix2(roll, pitch, yaw)
  gt.lintransform!(body.grid, M, translation; reset_fields=reset_fields)

  body.Oaxis[:,:] = M*body.Oaxis
  body.O[:] += translation

  nothing
end

##### COMMON INTERNAL FUNCTIONS  ###############################################
function _get_controlpoint(grid::gt.GridTriangleSurface, args...)
  cellnodes = gt.get_cellnodes(grid, args...)
  len = maximum([norm(cellnodes[1]-v) for v in cellnodes[2:end]])
  return mean(cellnodes) + 0.005*len*gt.get_normal(grid, args...)
end

function _solvedflag(self::BodyTypes, val::Bool)
  add_field(self, "solved", [val])
end

##### END OF ABSTRACT BODY #####################################################
