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
  * ` `             : .

  and the following functions
  ```julia
  ```
"""
abstract type AbstractBody end


for header_name in ["nonliftingbody"]
  include("MyPanel_"*header_name*".jl")
end

# Implementations of AbstractGrid
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

##### COMMON INTERNAL FUNCTIONS  ###############################################
function _get_controlpoint(grid::gt.GridTriangleSurface, args...)
  return mean(gt.get_cellnodes(grid, args...))+0.01*gt.get_normal(grid, args...)
end

##### END OF ABSTRACT BODY #####################################################
