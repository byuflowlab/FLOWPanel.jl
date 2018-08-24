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

  and the following functions

  ```julia

  # Solve for distributions of potential flow
  function solve(self::BodyTypes, Vinfs::Array{Array{T,1},1}, args...
                                                              ) where {T<:Real}
    .
    .
    .
  end


  # Returns the requested field.
  function get_field(self::BodyTypes, field_name::String)
    .
    .
    .
  end

  # Returns the requested field value.
  function get_field_val(self::BodyTypes, field_name::String, i::Int64;
                        _check::Bool=true)
    .
    .
    .
  end

  # Adds a new field to the body.
  function add_field(self::BodyTypes, field_name::String, field_data)
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
const BodyTypes = Union{NonLiftingBody}



##### COMMON FUNCTIONS  ########################################################
"""
  `save(body::BodyTypes, filename::String; opt_args...)`

Outputs a vtk file of this body. See GeometricTools.save(grid, ...) for a
descrition of optional arguments `opt_args...`.
"""
function save(body::BodyTypes, filename::String; opt_args...)
  gt.save(body.grid, filename; opt_args...)
end

"""
  `get_controlpoint(body::BodyTypes, i::Int64 or coor::Array{Int64,1})`

Returns the control point on the i-th panel.
"""
function get_controlpoint(body::BodyTypes, args...)
  return _get_controlpoint(body.grid, args...)
end

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
end

##### COMMON INTERNAL FUNCTIONS  ###############################################
function _get_controlpoint(grid::gt.GridTriangleSurface, args...)
  cellnodes = gt.get_cellnodes(grid, args...)
  len = maximum([norm(cellnodes[1]-v) for v in cellnodes[2:end]])
  return mean(cellnodes) + 0.005*len*gt.get_normal(grid, args...)
end

##### END OF ABSTRACT BODY #####################################################
