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
* `characteristiclength::Function`    : Function for computing the characteristic
                                        length of each panel used to offset each
                                        control point

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

    # Returns the potential induced by the body on the targets `targets`. It
    # adds the potential at the i-th target to out[i].
    function _phi(self::AbstractBody, targets::Array{Array{T,1},1},
    out::Array{Array{T,1},1}, args...; optargs...) where{T<:Real}
        .
        .
        .
    end
```
"""
abstract type AbstractBody{E<:AbstractElement, N} end



##### COMMON FUNCTIONS  ########################################################

"""
  `Uind!(self::AbstractBody, targets, out, args...; optargs...)

Returns the velocity induced by the body on the targets `targets`, which is a
3xn matrix. It adds the velocity at the i-th target to `out[:, i]`.
"""
function Uind!(self::AbstractBody, targets, out, args...; optargs...)

    # ERROR CASES
    if check_solved(self)==false
        error("Body hasn't been solved yet."*
              " Please call `solve()` function first.")
    end

    _Uind!(self, targets, out, args...; optargs...)
end

"""
  `phi!(self::AbstractBody, targets, out, args...; optargs...)

Returns the potential induced by the body on the targets `targets`. It adds the
potential at the i-th target to `out[:, i]`.
"""
function phi!(self::AbstractBody, targets, out, args...; optargs...)

    # ERROR CASES
    if check_solved(self)==false
        error("Body hasn't been solved yet."*
              " Please call `solve()` function first.")
    end

    _phi!(self, targets, out, args...; optargs...)
end

"""
  `save(body::AbstractBody, filename::String; optargs...)`

Outputs a vtk file of this body. See GeometricTools.save(grid, ...) for a
description of optional arguments `optargs...`.
"""
function save(body::AbstractBody, filename::String; out_cellindex::Bool=false,
                                                 out_cellindexdim::Array{Int64,1}=Int64[],
                                                 out_nodeindex::Bool=false,
                                                 out_controlpoints::Bool=false,
                                                 out_wake::Bool=true,
                                                 debug::Bool=false,
                                                 _len::Number=1.0,
                                                 _upper::Bool=true,
                                                 optargs...)

    str = ""

    # Add special fields
    if out_cellindex || debug
        gt.add_field(body.grid, "cellindex", "scalar",
                            [i for i in 1:body.ncells], "cell")
    end

    if out_nodeindex || debug
        gt.add_field(body.grid, "nodeindex", "scalar",
                            [i for i in 1:body.nnodes], "node")
    end

    _out_cellindexdim = debug && length(out_cellindexdim)==0 ? [1, 2] : out_cellindexdim
    for dim in _out_cellindexdim
        ndivs = gt.get_ndivscells(body.grid)[1:2]
        data = [ Base._ind2sub(ndivs, i)[dim] for i in 1:body.ncells]
        gt.add_field(body.grid, "cellindexdim$(dim)", "scalar", data, "cell")
    end

    # Outputs control points
    if out_controlpoints || debug
        str *= save_controlpoints(body, filename; debug=debug, optargs...)
    end

  # # Outputs wake
  # if out_wake || debug
  #   # Case that body is not a LiftingBody
  #   try
  #      body::LBodyTypes
  #      if body.nnodesTE-1 != 0
  #          str *= _savewake(body, filename; len=_len, upper=_upper, optargs...)
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
    return str*gt.save(body.grid, filename; format="vtk", optargs...)

end

"""
  `save_controlpoints(body::AbstractBody, filename::String;
suffix::String="_cp", optargs...)`

Outputs a vtk file with the control points of the body along with associated
normals and surface velocities.
"""
function save_controlpoints(body::AbstractBody, filename::String; debug=false,
                                              suffix::String="_cp", optargs...)

    # Normals
    normals = _calc_normals(body)

    # Control points
    CPs = _calc_controlpoints(body, normals)

    # Normals
    # normals = collect(eachcol(normals))
    normals = eachcol(normals)
    data = [
            Dict( "field_name"  => "normal",
                  "field_type"  => "vector",
                  "field_data"  => normals)
            ]

    if debug
        # Tangent unitary vectors
        tangents = _calc_tangents(body)
        tangents = eachcol(tangents)
        push!(data,
                Dict( "field_name"  => "tangent",
                      "field_type"  => "vector",
                      "field_data"  => tangents))

        # Oblique unitary vectors
        obliques = _calc_obliques(body)
        obliques = eachcol(obliques)
        push!(data,
                Dict( "field_name"  => "oblique",
                      "field_type"  => "vector",
                      "field_data"  => obliques))

    end

    if check_solved(body) && debug
        # Surface velocity
        Usurf = hcat(collect(get_fieldval(body, "Uinf", i) for i in 1:body.ncells)...)
        Uind!(body, CPs, Usurf)


        Usurf = collect(eachcol(Usurf))
        push!(data,
                    Dict( "field_name"  => "U",
                          "field_type"  => "vector",
                          "field_data"  => Usurf)
             )

        # Save surface velocity as a field
        add_field(body, "U", "vector", Usurf, "cell")


        # Surface potential
        phis = zeros(body.ncells)
        phi!(body, CPs, phis)

        push!(data,
                  Dict( "field_name"  => "phi",
                        "field_type"  => "scalar",
                        "field_data"  => phis)
           )

        # Save surface potential as a field
        add_field(body, "phi", "scalar", phis, "cell")
    end

    CPs = collect(eachcol(CPs))

    # Generates vtk
    return gt.generateVTK(filename*suffix, CPs; point_data=data, optargs...)
end

# """
#   `get_controlpoint(body::AbstractBody, i::Int64 or coor::Array{Int64,1})`
#
# Returns the control point on the i-th panel.
# """
# function get_controlpoint(body::AbstractBody, args...)
#   return _get_controlpoint(body.grid, args...)
# end

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
# function _get_controlpoint(grid::gt.GridTriangleSurface, args...)
#   return _get_controlpoint(gt.get_cellnodes(grid, args...),
#                                                   gt.get_normal(grid, args...))
# end
#
# function _get_controlpoint(cellnodes::Array{Array{T1,1},1}, normal::Array{T2,1};
#                                     off::Real=0.005) where{T1<:RType, T2<:RType}
#   len = maximum([norm(cellnodes[1]-v) for v in cellnodes[2:end]])
#   return mean(cellnodes) + off*len*normal
# end

function characteristiclength_bbox(nodes, panel)
    # Characteristic length: Diagonal of bounding box
    min1 = nodes[1, first(panel)]
    min2 = nodes[2, first(panel)]
    min3 = nodes[3, first(panel)]
    max1, max2, max3 = min1, min2, min3

    for ni in panel
        if nodes[1, ni] <= min1; min1 = nodes[1, ni]; end
        if nodes[2, ni] <= min2; min2 = nodes[2, ni]; end
        if nodes[3, ni] <= min3; min3 = nodes[3, ni]; end
        if nodes[1, ni] >= max1; max1 = nodes[1, ni]; end
        if nodes[2, ni] >= max2; max2 = nodes[2, ni]; end
        if nodes[3, ni] >= max3; max3 = nodes[3, ni]; end
    end

    l = sqrt((max1-min1)^2 + (max2-min2)^2 + (max3-min3)^2)

    return l
end

function characteristiclength_maxdist(nodes, panel)

    # Characteristic length: Maximum node distance
    l = 0
    n0 = first(panel)
    for ni in panel
        this_l = (nodes[1, n0] - nodes[1, ni])^2
        this_l += (nodes[2, n0] - nodes[2, ni])^2
        this_l += (nodes[3, n0] - nodes[3, ni])^2
        this_l = sqrt(this_l)

        if this_l > l
            l = this_l
        end
    end

    return l
end

characteristiclength_sqrtarea(nodes, panel) = sqrt(gt._get_area(nodes, panel))

function _calc_controlpoints!(grid::gt.GridTriangleSurface,
                                controlpoints, normals; off::Real=0.005,
                                characteristiclength::Function=characteristiclength_sqrtarea)

    lin = LinearIndices(grid._ndivsnodes)
    ndivscells = vcat(grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    nodes = grid.orggrid.nodes

    for pi in 1:grid.ncells

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                                grid, pi, lin, ndivscells, cin)

        controlpoints[:, pi] .= 0

        # Control point: Average point between nodes
        for ni in panel
            controlpoints[1, pi] += nodes[1, ni]
            controlpoints[2, pi] += nodes[2, ni]
            controlpoints[3, pi] += nodes[3, ni]
        end
        controlpoints[:, pi] /= length(panel)

        # Control point: Convert average point into centroid
        t1, t2, t3 = gt._calc_t1(nodes, panel), gt._calc_t2(nodes, panel), gt._calc_t3(nodes, panel) # Tangent unit vector
        o1, o2, o3 = gt._calc_o1(nodes, panel), gt._calc_o2(nodes, panel), gt._calc_o3(nodes, panel) # Oblique unit vector
        c1, c2, c3 = controlpoints[1, pi], controlpoints[2, pi], controlpoints[3, pi]                # Center of panel

        x = controlpoints[1, pi]*0
        y = controlpoints[1, pi]*0
        for ni in panel
            x += (nodes[1, ni]-c1)*t1
            x += (nodes[2, ni]-c2)*t2
            x += (nodes[3, ni]-c3)*t3
            y += (nodes[1, ni]-c1)*o1
            y += (nodes[2, ni]-c2)*o2
            y += (nodes[3, ni]-c3)*o3
        end
        x /= length(panel)
        y /= length(panel)

        controlpoints[1, pi] += x*t1 + y*o1
        controlpoints[2, pi] += x*t2 + y*o2
        controlpoints[3, pi] += x*t3 + y*o3

        l = characteristiclength(nodes, panel)

        # Offset the controlpoint in the normal direction
        controlpoints[1, pi] += off*l*normals[1, pi]
        controlpoints[2, pi] += off*l*normals[2, pi]
        controlpoints[3, pi] += off*l*normals[3, pi]
    end

end
function _calc_controlpoints(grid::gt.GridTriangleSurface, args...; optargs...)
    controlpoints = zeros(3, grid.ncells)
    _calc_controlpoints!(grid, controlpoints, args...; optargs...)
    return controlpoints
end
function _calc_controlpoints(self::AbstractBody, args...; optargs...)
    return _calc_controlpoints(self.grid, args...; off=self.CPoffset,
                                characteristiclength=self.characteristiclength,
                                optargs...)
end

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


function _calc_tangents!(grid::gt.GridTriangleSurface, tangents)

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
        tangents[1, pi] = gt._calc_t1(grid.orggrid.nodes, panel)
        tangents[2, pi] = gt._calc_t2(grid.orggrid.nodes, panel)
        tangents[3, pi] = gt._calc_t3(grid.orggrid.nodes, panel)
    end
end
function _calc_tangents(grid::gt.GridTriangleSurface)
    tangents = zeros(3, grid.ncells)
    _calc_tangents!(grid, tangents)
    return tangents
end
_calc_tangents(self::AbstractBody) = _calc_tangents(self.grid)


function _calc_obliques!(grid::gt.GridTriangleSurface, obliques)

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
        obliques[1, pi] = gt._calc_o1(grid.orggrid.nodes, panel)
        obliques[2, pi] = gt._calc_o2(grid.orggrid.nodes, panel)
        obliques[3, pi] = gt._calc_o3(grid.orggrid.nodes, panel)
    end
end
function _calc_obliques(grid::gt.GridTriangleSurface)
    obliques = zeros(3, grid.ncells)
    _calc_obliques!(grid, obliques)
    return obliques
end
_calc_obliques(self::AbstractBody) = _calc_obliques(self.grid)


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
