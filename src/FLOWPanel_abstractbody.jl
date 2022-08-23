#=##############################################################################
# DESCRIPTION
    Abstract body type definition.
# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jun 2018
  * License   : MIT License
=###############################################################################


################################################################################
# ABSTRACT BODY TYPE
################################################################################
"""
Abstract type `AbstractBody{N, E<:AbstractElement}` where `N` is the number of
element types of this body and `E` is an Union of containing the `N` element
types.

Implementations of AbstractBody are expected to have the following fields.
* `grid::gt.GridTriangleSurface `     : Paneled geometry
* `nnodes::Int64`                     : Number of nodes
* `ncells::Int64`                     : Number of cells
* `fields::Array{String, 1}`          : Available fields (solutions)
* `Oaxis::Array{T,2} where {T<:Real}` : Coordinate system of original grid
* `O::Array{T,1} where {T<:Real}`     : Position of CS of original grid
* `strength::Array{<:Real,2}`         : Strength of each element of each type
* `CPoffset::Real`                    : Control point offset in normal direction

and the following functions

```julia

    # Solve for distributions of potential flow
    function solve(self::AbstractBody, Vinfs::Array{Array{T,1},1}, args...
    ) where {T<:Real}
        .
        .
        .
    end

    # Returns the velocity induced by the body on the targets `targets`. It adds
    # the velocity at the i-th target to out[i].
    function _Vind(self::AbstractBody, targets::Array{Array{T,1},1},
    out::Array{Array{T,1},1}, args...; optargs...) where{T<:Real}
        .
        .
        .
    end
```
"""
abstract type AbstractBody{E<:AbstractElement, N} end



##### COMMON FUNCTIONS  ########################################################

# """
#   `Uind(self::AbstractBody, targets::Array{Array{T,1},1},
#                   out::Array{Array{T,1},1}, args...; optargs...) where{T<:Real}`
#
# Returns the velocity induced by the body on the targets `targets`. It adds the
# velocity at the i-th target to out[i].
# """
# function Uind(self::AbstractBody, targets::Arr1,
#                   out::Arr2, args...; optargs...
#                                                   ) where{T1, Arr1<:AbstractArray{T1,2},
#                                                           T2, Arr2<:AbstractArray{T2}}
#
#   # ERROR CASES
#   if check_solved(self)==false
#     error("Body hasn't been solved yet. Please call `solve()` function first.")
#   end
#
#   _Uind(self, targets, out, args...; optargs...)
# end

"""
  `save(body::AbstractBody, filename::String; opt_args...)`

Outputs a vtk file of this body. See GeometricTools.save(grid, ...) for a
description of optional arguments `opt_args...`.
"""
function save(body::AbstractBody, filename::String; out_cellindex::Bool=false,
                                                 out_cellindexdim::Array{Int64,1}=Int64[],
                                                 out_nodeindex::Bool=false,
                                                 out_controlpoints::Bool=false,
                                                 out_wake::Bool=true,
                                                 _len::RType=1.0,
                                                 _upper::Bool=true,
                                                 opt_args...)

  str = ""

  # Add special fields
  if out_cellindex
    gt.add_field(body.grid, "cellindex", "scalar",
                    [i for i in 1:body.ncells], "cell")
  end
  if out_nodeindex
    gt.add_field(body.grid, "nodeindex", "scalar",
                    [i for i in 1:body.nnodes], "node")
  end
  for dim in out_cellindexdim
    ndivs = gt.get_ndivscells(body.grid)[1:2]
    data = [ Base._ind2sub(ndivs, i)[dim] for i in 1:body.ncells]
    gt.add_field(body.grid, "cellindexdim$(dim)", "scalar", data, "cell")
  end

  # Outputs control points
  if out_controlpoints
    str *= save_controlpoints(body, filename; opt_args...)
  end

  # # Outputs wake
  # if out_wake
  #   # Case that body is not a LiftingBody
  #   try
  #      body::LBodyTypes
  #      if body.nnodesTE-1 != 0
  #          str *= _savewake(body, filename; len=_len, upper=_upper, opt_args...)
  #      end
  #    catch e
  #      if isa(e, TypeError)
  #        nothing
  #      else
  #        rethrow(e)
  #      end
  #    end
  # end

  # Saves body
  return str*gt.save(body.grid, filename; opt_args...)

end

"""
  `save_controlpoints(body::AbstractBody, filename::String;
suffix::String="_cp", opt_args...)`

Outputs a vtk file with the control points of the body along with associated
normals and surface velocities.
"""
function save_controlpoints(body::AbstractBody, filename::String;
                                              suffix::String="_cp", opt_args...)

    # Normals
    normals = _calc_normals(body)

    # Control points
    CPs = _calc_controlpoints(body, normals; off=body.CPoffset)
    CPs = collect(eachcol(CPs))
    normals = collect(eachcol(normals))

    # Normals
    data = [
            Dict( "field_name"  => "normal",
                  "field_type"  => "vector",
                  "field_data"  => normals)
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
    return gt.generateVTK(filename*suffix, CPs; point_data=data, opt_args...)
end

"""
  `get_controlpoint(body::AbstractBody, i::Int64 or coor::Array{Int64,1})`

Returns the control point on the i-th panel.
"""
function get_controlpoint(body::AbstractBody, args...)
  return _get_controlpoint(body.grid, args...)
end

"""
  `get_ndivscells(body::AbstractBody)`

Returns a tuple with the number of cells in each parametric dimension
"""
get_ndivscells(body::AbstractBody) = body.grid._ndivscells

"""
  `get_ndivsnodes(body::AbstractBody)`

Returns a tuple with the number of nodes in each parametric dimension
"""
get_ndivsnodes(body::AbstractBody) = body.grid._ndivsnodes

"""
  `get_unitvectors(body::AbstractBody, i::Int64 or coor::Array{Int64,1})`

Returns the unit vectors t,n,o of the i-th panel, with t the tanget vector,
n normal, and o oblique.
"""
get_unitvectors(body::AbstractBody, args...) = gt.get_unitvectors(body.grid, args...)

"""
  `get_normal(body::AbstractBody, i::Int64 or coor::Array{Int64,1})`

Returns the normal vector the i-th panel.
"""
get_normal(body::AbstractBody, args...) = gt.get_normal(body.grid, args...)

"""
  `get_field(self::AbstractBody, field_name::String)`

Returns the requested field.
"""
function get_field(self::AbstractBody, field_name::String)
  if !(field_name in self.fields)
    error("Field $field_name not found! Available fields: $(self.fields)")
  end

  return self.grid.field[field_name]
end

"""
  `get_field(self::AbstractBody, field_name::String, i::Int64)`

Returns the requested field value. Give it `_check=false` to skip checking logic
for faster operation.
"""
function get_fieldval(self::AbstractBody, field_name::String, i::Int64;
                      _check::Bool=true)
  if _check
    if i<1
      error("Invalid index $i.")
    end
  end

  return gt.get_fieldval(self.grid, field_name, i)
end

"""
  `get_fieldval(self::AbstractBody, field_name, coor)`

Returns the value of node of coordinates `coor` (1-indexed) in the field
'field_name'.
"""
function get_fieldval(self::AbstractBody, field_name::String, coor::Array{Int64,1})
  return gt.get_fieldval(self.grid, field_name, coor)
end


"""
  `add_field(self::AbstractBody, field_name::String, field_type::String,
                    field_data, entry_type::String)`

Adds a new field to the body. It overwrites the field if it already existed.
`field_type` is either `"scalar"` or `"vector"`. `entry_type` is one out of
`["node", "cell", "system"]`.
"""
function add_field(self::AbstractBody, field_name::String, field_type::String,
                    field_data, entry_type::String; raise_warn=false)

  # Adds field to grid
  gt.add_field(self.grid, field_name, field_type,
                field_data, entry_type; raise_warn=raise_warn)

  # Registers the field
  if !(field_name in self.fields)
    push!(self.fields, field_name)
  end

  nothing
end

"""
  `check_field(self::AbstractBody, field_name::String)`

Returns `true` of the body has the field `field_name`. Returns false otherwise.
"""
check_field(self::AbstractBody, field_name::String) = field_name in self.fields


"""
  `check_solved(self::AbstractBody)`

Returns `true` of the body has been solved. Returns false otherwise.
"""
function check_solved(self::AbstractBody)
  if check_field(self, "solved")
    return get_fieldval(self, "solved", 1)
  else
    return false
  end
end


"""
  `rotate(body::AbstractBody, roll::Real, pitch::Real, yaw::Real;
translation::Array{T, 1}=zeros(3), reset_fields::Bool=true
) where{T<:Real}`

Rotates and translates the body by the given axial angles.

NOTE: Naming follows aircraft convention, with
* roll:   rotation about x-axis.
* pitch:  rotation about y-axis.
* yaw:    rotation about z-axis.
"""
function rotate(body::AbstractBody, roll::Number, pitch::Number, yaw::Number;
                  translation::Array{<:Number, 1}=zeros(3),
                  reset_fields::Bool=true
                )

  M = gt.rotation_matrix2(roll, pitch, yaw)
  gt.lintransform!(body.grid, M, translation; reset_fields=reset_fields)

  body.Oaxis[:,:] = M*body.Oaxis
  body.O[:] += translation

  nothing
end

##### COMMON INTERNAL FUNCTIONS  ###############################################
function _get_controlpoint(grid::gt.GridTriangleSurface, args...)
  return _get_controlpoint(gt.get_cellnodes(grid, args...),
                                                  gt.get_normal(grid, args...))
end

function _get_controlpoint(cellnodes::Array{Array{T1,1},1}, normal::Array{T2,1};
                                    off::Real=0.005) where{T1<:RType, T2<:RType}
  len = maximum([norm(cellnodes[1]-v) for v in cellnodes[2:end]])
  return mean(cellnodes) + off*len*normal
end


function _calc_controlpoints!(grid::gt.GridTriangleSurface,
                                controlpoints, normals; off::Real=0.005)

    lin = LinearIndices(grid._ndivsnodes)
    ndivscells = vcat(grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    for pi in 1:grid.ncells

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                                grid, pi, lin, ndivscells, cin)

        controlpoints[:, pi] .*= 0

        # Average point between nodes
        for ni in panel
            controlpoints[1, pi] += grid.orggrid.nodes[1, ni]
            controlpoints[2, pi] += grid.orggrid.nodes[2, ni]
            controlpoints[3, pi] += grid.orggrid.nodes[3, ni]
        end
        controlpoints[:, pi] /= length(panel)

        # Characteristic length
        min1 = grid.orggrid.nodes[1, first(panel)]
        min2 = grid.orggrid.nodes[2, first(panel)]
        min3 = grid.orggrid.nodes[3, first(panel)]
        max1, max2, max3 = min1, min2, min3

        for ni in panel
            if grid.orggrid.nodes[1, ni] <= min1; min1 = grid.orggrid.nodes[1, ni]; end
            if grid.orggrid.nodes[2, ni] <= min2; min2 = grid.orggrid.nodes[2, ni]; end
            if grid.orggrid.nodes[3, ni] <= min3; min3 = grid.orggrid.nodes[3, ni]; end
            if grid.orggrid.nodes[1, ni] >= max1; max1 = grid.orggrid.nodes[1, ni]; end
            if grid.orggrid.nodes[2, ni] >= max2; max2 = grid.orggrid.nodes[2, ni]; end
            if grid.orggrid.nodes[3, ni] >= max3; max3 = grid.orggrid.nodes[3, ni]; end
        end

        # Offset the controlpoint in the normal direction
        controlpoints[1, pi] += off*(max1-min1)*normals[1, pi]
        controlpoints[2, pi] += off*(max2-min2)*normals[2, pi]
        controlpoints[3, pi] += off*(max3-min3)*normals[3, pi]
    end

end

function _calc_controlpoints(grid::gt.GridTriangleSurface, args...; optargs...)
    controlpoints = zeros(3, grid.ncells)
    _calc_controlpoints!(grid, controlpoints, args...; optargs...)
    return controlpoints
end
_calc_controlpoints(self::AbstractBody, args...; optargs...) = _calc_controlpoints(self.grid, args...; optargs...)

function _calc_normals!(grid::gt.GridTriangleSurface, normals)

    lin = LinearIndices(grid._ndivsnodes)
    ndivscells = vcat(grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    for pi in 1:grid.ncells
        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                                grid, pi, lin, ndivscells, cin)
        normals[1, pi] = gt._calc_n1(grid.orggrid.nodes, panel)
        normals[2, pi] = gt._calc_n2(grid.orggrid.nodes, panel)
        normals[3, pi] = gt._calc_n3(grid.orggrid.nodes, panel)
    end
end
function _calc_normals(grid::gt.GridTriangleSurface)
    normals = zeros(3, grid.ncells)
    _calc_normals!(grid, normals)
    return normals
end
_calc_normals(self::AbstractBody) = _calc_normals(self.grid)

function _solvedflag(self::AbstractBody, val::Bool)
  add_field(self, "solved", "scalar", [val], "system")
end

"Count the number of types in an Union type"
function _count(type::Type)
    if type isa Union
        return _count(type.a) + _count(type.b)
    else
        return 1
    end
end
##### END OF ABSTRACT BODY #####################################################
