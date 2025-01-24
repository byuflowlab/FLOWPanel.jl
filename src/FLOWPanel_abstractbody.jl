#=##############################################################################
# DESCRIPTION
    Abstract body type definition.
# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Jun 2018
  * License     : MIT License
=###############################################################################


################################################################################
# ABSTRACT BODY TYPE
################################################################################
"""
Abstract type `AbstractBody{N, E<:AbstractElement}` where `N` is the number of
element types in this body and `E` is an Union containing the `N` element
types.

Implementations of AbstractBody are expected to have the following fields
* `grid::GeometricTools.GridTriangleSurface `     : Paneled geometry
* `nnodes::Int`                       : Number of nodes
* `ncells::Int`                       : Number of cells
* `fields::Vector{String}`            : Available fields (solutions)
* `Oaxis::Matrix`                     : Coordinate system of body w.r.t global (3x3 matrix)
* `O::Vector`                         : Origin of body w.r.t. global (3-dim vector)
* `strength::Matrix`                  : Strength of each element of each type (ncells x N matrix)
* `CPoffset::Real`                    : Control point offset in normal direction
* `characteristiclength::Function`    : Function for computing the characteristic
                                        length of each panel used to offset each
                                        control point
* `kerneloffset::Real`                : Kernel offset to avoid singularities
* `kernelcutoff::Real`                : Kernel cutoff to avoid singularities

and the following functions

```julia

    # Imposes boundary conditions to solve for element strengths.
    function solve(self::AbstractBody, Uinfs::Array{<:Real, 2}, args...)
        .
        .
        .
    end

    # Outputs the body as a vtk file
    function save(self::AbstractBody, args...; optargs...)
        .
        .
        .
    end

    # Returns the dimensions of the system of equations solved by `solve(...)`
    # as a tuple `(m, n)`, where `m` is the number of equations and `n` is the
    # number of unknowns (usually `m==n`).
    function _get_Gdims(self::AbstractBody)
        .
        .
        .
    end

    # Returns the velocity induced by the body on the targets `targets`. It adds
    # the velocity at the i-th target to out[:, i].
    function _Uind(self::AbstractBody, targets::Array{<:Real, 2},
                    out::Array{<:Real, 2}, args...; optargs...)
        .
        .
        .
    end

    # Returns the potential induced by the body on the targets `targets`. It
    # adds the potential at the i-th target to out[i].
    function _phi(self::AbstractBody, targets::Array{<:Real, 2},
                    out::Array{<:Real, 1}, args...; optargs...)
        .
        .
        .
    end
```
"""
abstract type AbstractBody{E<:AbstractElement, N} end

"""
    `solve(body::AbstractBody, Uinfs::Array{<:Real, 2})`

Impose boundary conditions to solve for element strengths. `Uinds[:, i]` is the
velocity at the i-th control point used in the boundary condition.
"""
function solve(self::AbstractBody, Uinfs::AbstractArray{<:Number, 2})
    error("solve(...) for body type $(typeof(self)) has not been implemented yet!")
end

##### COMMON FUNCTIONS  ########################################################

"""
    Uind!(self::AbstractBody, targets, out, args...; optargs...)

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
    phi!(self::AbstractBody, targets, out, args...; optargs...)

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
    save_base(body::AbstractBody, filename::String; optargs...)

Outputs a vtk file of this body. See
`GeometricTools.save(::GridExtentions, ::String)` for a description of optional
arguments `optargs...`.
"""
function save_base(body::AbstractBody, filename::String; out_cellindex::Bool=false,
                                                 out_cellindexdim::Array{Int64,1}=Int64[],
                                                 out_nodeindex::Bool=false,
                                                 out_controlpoints::Bool=false,
                                                 debug::Bool=false,
                                                 suffix::String="",
                                                 optargs...)

    str = ""

    # Add special fields
    if out_cellindex || debug
        gt.add_field(body.grid, "cellindex", "scalar",
                        [i for i in 1:body.ncells], "cell"; raise_warn=false)
    end

    if out_nodeindex || debug
        gt.add_field(body.grid, "nodeindex", "scalar",
                        [i for i in 1:body.nnodes], "node"; raise_warn=false)
    end

    _out_cellindexdim = debug && length(out_cellindexdim)==0 ? body.grid.dims==2 ? [1] : [1, 2] : out_cellindexdim
    for dim in _out_cellindexdim
        ndivs = gt.get_ndivscells(body.grid)[1:2]
        data = [ Base._ind2sub(ndivs, i)[dim] for i in 1:body.ncells]
        gt.add_field(body.grid, "cellindexdim$(dim)", "scalar", data, "cell";
                                                            raise_warn=false)
    end

    # Outputs control points
    if out_controlpoints || debug
        str *= save_controlpoints(body, filename; debug=debug,
                                                optargs..., suffix="_cp")
    end

    # Saves body
    str *= gt.save(body.grid, filename*suffix; format="vtk", optargs...)

    # Return path to files
    return str

end

"""
    save_controlpoints(body::AbstractBody, filename::String;
                                suffix::String="_cp", optargs...)

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
        if check_field(body, "U")==false
            add_field(body, "U", "vector", Usurf, "cell")
        end


        # Surface potential
        phis = zeros(body.ncells)
        phi!(body, CPs, phis)

        push!(data,
                  Dict( "field_name"  => "phi",
                        "field_type"  => "scalar",
                        "field_data"  => phis)
           )

        # Save surface potential as a field
        if check_field(body, "phi")==false
            add_field(body, "phi", "scalar", phis, "cell")
        end
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
    get_ndivscells(body::AbstractBody)

Returns a tuple with the number of cells in each parametric dimension
"""
get_ndivscells(body::AbstractBody) = deepcopy(body.grid._ndivscells)

"""
    get_ndivsnodes(body::AbstractBody)

Returns a tuple with the number of nodes in each parametric dimension
"""
get_ndivsnodes(body::AbstractBody) = deepcopy(body.grid._ndivsnodes)


"""
    get_cart2lin_cells(self::AbstractBody)

Returns a `LinearIndices` that converts the coordinates (or "Cartesian index")
of a cell to its linear index.

!!! tip "Example"
```julia
    coordinates = (i, j)                # (i, j) coordinates of an arbitrary cell

    cart2lin = get_cart2lin_cells(body)

    index = cart2lin(coordinates...)    # Linear index of the cell
```
"""
function get_cart2lin_cells(self::AbstractBody)
    # Remove any quasi-dimensions
    ndivscells = [i for i in get_ndivscells(self) if i!=0]

    # Return linear indexing
    return LinearIndices(Tuple(ndivscells))
end

"""
    get_cart2lin_nodes(self::AbstractBody)

Returns a `LinearIndices` that converts the coordinates (or "Cartesian index")
of a node to its linear index.

!!! tip "Example"
```julia
    coordinates = (i, j)                # (i, j) coordinates of an arbitrary node

    cart2lin = get_cart2lin_nodes(body)

    index = cart2lin(coordinates...)    # Linear index of the node
```
"""
function get_cart2lin_nodes(self::AbstractBody)
    # Remove any quasi-dimensions
    ndivsnodes = [i for i in get_ndivsnodes(self) if i!=0]

    # Return linear indexing
    return LinearIndices(Tuple(ndivsnodes))
end

get_strength(self::AbstractBody, i) = view(self.strength, i, :)
get_nstrengths(self::AbstractBody) = size(self.strength, 1)
set_strength(self::AbstractBody, i::Int, val) = (self.strength[i, :] .= val)

"""
    get_unitvectors(body::AbstractBody, i::Int64 or coor::Array{Int64,1})

Returns the unit vectors `t`,`n`,`o` of the i-th panel, with `t` the tanget
vector, `n` normal, and `o` oblique.
"""
get_unitvectors(body::AbstractBody, args...) = gt.get_unitvectors(body.grid, args...)

"""
    get_normal(body::AbstractBody, i::Int64 or coor::Array{Int64,1})

Returns the normal vector the i-th panel.
"""
get_normal(body::AbstractBody, args...) = gt.get_normal(body.grid, args...)

"""
    get_field(self::AbstractBody, field_name::String)

Returns the requested field.
"""
function get_field(self::AbstractBody, field_name::String)
    if !(field_name in self.fields)
        error("Field $field_name not found! Available fields: $(self.fields)")
    end

    return self.grid.field[field_name]
end

"""
    get_field(self::AbstractBody, field_name::String, i::Int)

Returns the requested field value of the i-th cell or node (depending of the
field type). Give it `_check=false` to skip checking logic for faster processing.
"""
function get_fieldval(self::AbstractBody, field_name::String, i::Int;
                      _check::Bool=true)
    if _check
        if i<1
            error("Invalid index $i.")
        end
    end

    return gt.get_fieldval(self.grid, field_name, i)
end

"""
    get_fieldval(self::AbstractBody, field_name::String, coor::Array{Int,1})

Returns the requested field value of the cell or node (depending of the field
type) of coordinates `coor` (1-indexed).
"""
function get_fieldval(self::AbstractBody, field_name::String, coor::Array{Int,1})
    return gt.get_fieldval(self.grid, field_name, coor)
end


"""
    add_field(self::AbstractBody, field_name::String, field_type::String,
                                                field_data, entry_type::String)

Adds a new field `field_name` to the body (overwriting the field if it already
existed).

**Expected arguments**

* `field_type=="scalar"`, then `field_data` is a vector of length `n`.
* `field_type=="vector"`, then `field_data` is an array of 3-dim vectors of
    length `n`.

* `entry_type=="node"`, then `n=length(field_data)` is the number of nodes in
    the body and `field_data[i]` is the field value at the i-th node.
* `entry_type=="cell"`, then `n=length(field_data)` is the number of cells in
    the body and `field_data[i]` is the field value at the i-th cell.
* `entry_type=="system"`, then `n=length(field_data)` is any arbritrary number,
    and `field_data` is a field for the whole body as a system without any
    data structure.
"""
function add_field(self::AbstractBody, field_name::String, field_type::String,
                    field_data, entry_type::String;
                    raise_warn=false, collectfield=true)

    # Add field to grid
    gt.add_field(self.grid, field_name, field_type,
                    collectfield ? collect(field_data) : field_data, entry_type;
                    raise_warn=raise_warn)

    # Register the field
    if !(field_name in self.fields)
        push!(self.fields, field_name)
    end

    nothing
end

"""
    addfields(body::AbstractBody,
                sourcefieldname::String, targetfieldname::String)

Adds field `sourcefieldname` to field `targetfieldname`.
"""
function addfields(body::AbstractBody,
                    sourcefieldname::String, targetfieldname::String)

    srcfield = get_field(body, sourcefieldname)["field_data"]
    trgfield = get_field(body, targetfieldname)["field_data"]

    if get_field(body, sourcefieldname)["field_type"]=="vector"

        for (Fsrc, Ftrg) in zip(srcfield, trgfield)
            Ftrg .+= Fsrc
        end

    else

        trgfield += srcfield

    end

end

"""
    remove_field(self::AbstractBody, field_name)

Removes field from body.
"""
function remove_field(self::AbstractBody, field_name)
    if check_field(self, field_name)

        i = findfirst(name->name==field_name, self.fields)

        splice!(self.fields, i)

        gt.remove_field(self.grid, field_name)
    end
end

"""
    check_field(self::AbstractBody, field_name::String)

Returns `true` of the body has the field `field_name`. Returns false otherwise.
"""
check_field(self::AbstractBody, field_name::String) = field_name in self.fields


"""
    check_solved(self::AbstractBody)

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
    `rotate!(body::AbstractBody, roll::Number, pitch::Number, yaw::Number;
                translation::Vector=zeros(3), reset_fields::Bool=true)`

Rotates the body by the given axial angles. Nomenclature follows the aircraft
convention, with
* roll:   rotation about local x-axis
* pitch:  rotation about local y-axis
* yaw:    rotation about local z-axis

Use keyword argument `translation` to also translate the body.
"""
function rotate!(body::AbstractBody, roll::Number, pitch::Number, yaw::Number;
                  translation::AbstractVector=zeros(3),
                  reset_fields::Bool=true
                )

    # Generate rotation matrix
    M = gt.rotation_matrix2(-roll, -pitch, -yaw)

    return rotatetranslate!(body, M, translation; reset_fields=reset_fields)
end

"""
    rotatetranslate!(body::AbstractBody, M::Matrix, T::Vector; optargs...)

Rotate and translate `body` by rotation matrix `M` and translation vector `T`.
"""
function rotatetranslate!(body::AbstractBody,
                            M::AbstractMatrix, T::AbstractVector; optargs...)

    @assert length(T)==3 ""*
        "Invalid translation vector $(T): it is not three-dimensional!"
    @assert ndims(M)==2 && size(M, 1)==3 && size(M, 2)==3 ""*
        "Invalid"

    # Bring back to the origin
    body.O .*= -1
    gt.lintransform!(body.grid, Im, body.O; optargs...)
    body.O .*= -1

    # Add translation to previous position
    body.O .+= T

    # Rotate and translate to new position
    gt.lintransform!(body.grid, M, body.O; optargs...)

    # Update coordinate system
    body.Oaxis .= M*body.Oaxis

    return nothing
end

"""
    set_coordinatesystem(body::AbstractBody, O::Vector, Oaxis::Matrix)

Redefines the local coordinate system of the body, where `O` is the new origin
and `Oaxis` is the matrix of new unit vectors.
"""
function set_coordinatesystem(body::AbstractBody,
                                O::AbstractVector, Oaxis::AbstractMatrix;
                                optargs...)

    # Undo its current coordinate system
    rotatetranslate!(body, inv(body.Oaxis), -body.O; optargs...)

    # Set new coordinate system
    rotatetranslate!(body, Oaxis, O; optargs...)

    return nothing
end


##### COMMON INTERNAL FUNCTIONS  ###############################################
"""
    characteristiclength_unitary(args...) -> 1
"""
characteristiclength_unitary(args...) = 1

"""
    characteristiclength_bbox(nodes::Matrix, panel::Vector{Int})

Returns the characteristic length of a panel calculated as the diagonal of
the minimum bounding box.
"""
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

"""
    characteristiclength_maxdist(nodes::Matrix, panel::Vector{Int})

Returns the characteristic length of a panel calculated as the maximum distance
between nodes.
"""
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

"""
    characteristiclength_sqrtarea(nodes::Matrix, panel::Vector{Int})

Returns the characteristic length of a panel calculated as the square-root of
its area.
"""
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

    nodes = grid._nodes

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
function _calc_controlpoints!(grid::gt.GridSurface2D,
                                controlpoints, normals; off::Real=0.005,
                                characteristiclength::Function=characteristiclength_sqrtarea)

    out, coor, lin, ndivscells, cin = gt.generate_getcellt_args!(grid)

    nodes = grid._nodes

    for pi in 1:grid.ncells

        coor[1] = pi
        panel = gt.get_cell_t!(out, grid, coor, lin, ndivscells)

        controlpoints[:, pi] .= 0

        # Control point: Average point between nodes
        for ni in panel
            controlpoints[1, pi] += nodes[1, ni]
            controlpoints[2, pi] += nodes[2, ni]
        end
        controlpoints[:, pi] /= length(panel)

        l = characteristiclength(nodes, panel)

        # Offset the controlpoint in the normal direction
        controlpoints[1, pi] += off*l*normals[1, pi]
        controlpoints[2, pi] += off*l*normals[2, pi]
        controlpoints[3, pi] += off*l*normals[3, pi]
    end

end
function _calc_controlpoints(grid::gt.AbstractGrid, args...; optargs...)
    controlpoints = zeros(3, grid.ncells)
    _calc_controlpoints!(grid, controlpoints, args...; optargs...)
    return controlpoints
end
function _calc_controlpoints!(self::AbstractBody, args...; optargs...)
    return _calc_controlpoints!(self.grid, args...;
                                off=self.CPoffset,
                                characteristiclength=self.characteristiclength,
                                optargs...)
end
function _calc_controlpoints(self::AbstractBody, args...; optargs...)
    return _calc_controlpoints(self.grid, args...;
                                off=self.CPoffset,
                                characteristiclength=self.characteristiclength,
                                optargs...)
end

"""
    calc_controlpoints!(body::AbstractBody, controlpoints::Matrix, normals::Matrix)

Calculates the control point of every cell in `body` and stores them in the 3xN
matrix `controlpoints`. It uses `body.CPoffset`, `body.charateristiclength`, and
`normals` to offset the control points off the surface in the normal direction.

**Output:** `controlpoints[:, i]` is the control point of the i-th cell (linearly
indexed).

!!! tip
    Use `normals = calc_normals(body)` to calculate the normals.
"""
const calc_controlpoints! = _calc_controlpoints!

"""
    calc_controlpoints(body::AbstractBody)

Calculates the control point of every cell in `body` returning a 3xN matrix.

See `calc_controlpoints!` documentation for more details.
"""
const calc_controlpoints = _calc_controlpoints


function _calc_normals(grid::gt.AbstractGrid)
    normals = zeros(3, grid.ncells)
    _calc_normals!(grid, normals)
    return normals
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
        normals[1, pi] = gt._calc_n1(grid._nodes, panel)
        normals[2, pi] = gt._calc_n2(grid._nodes, panel)
        normals[3, pi] = gt._calc_n3(grid._nodes, panel)
    end
end
function _calc_normals!(grid::gt.GridSurface2D, normals)

    out, coor, lin, ndivscells, cin = gt.generate_getcellt_args!(grid)

    for pi in 1:grid.ncells
        coor[1] = pi
        panel = gt.get_cell_t!(out, grid, coor, lin, ndivscells)
        normals[1, pi] = gt._calc_n1_2D(grid._nodes, panel)
        normals[2, pi] = gt._calc_n2_2D(grid._nodes, panel)
    end
end

function _calc_normals!(self::AbstractBody, normals; flipbyCPoffset=false)
    _calc_normals!(self.grid, normals)
    if flipbyCPoffset
        normals .*= sign(self.CPoffset) != 0 ? sign(self.CPoffset) : 1
    end

    return normals
end
function _calc_normals(self::AbstractBody; flipbyCPoffset=false)
    normals = _calc_normals(self.grid)
    if flipbyCPoffset
        normals .*= sign(self.CPoffset) != 0 ? sign(self.CPoffset) : 1
    end

    return normals
end

"""
    calc_normals!(body::AbstractBody, normals::Matrix)

Calculates the normal vector of every cell in `body` and stores them in the 3xN
matrix `normals`.

**Output:** `normals[:, i]` is the normal vector of the i-th cell (linearly
indexed).

!!! tip "Tip: Cartesian to linear indices"

    Normals can be accessed through their (i, j) coordinates (or "Cartesian
    indices") as follows:

    ```julia
        coordinates = (i, j)

        ndivscells = get_ndivscells(body)
        lin = LinearIndices(Tuple(ndivscells))
    ```
"""
const calc_normals! = _calc_normals!

"""
    calc_normals(self::AbstractBody)

Calculates the normal vector of every cell in `grid` returning a 3xN matrix.

See `calc_normals!` documentation for more details.
"""
const calc_normals = _calc_normals


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
        tangents[1, pi] = gt._calc_t1(grid._nodes, panel)
        tangents[2, pi] = gt._calc_t2(grid._nodes, panel)
        tangents[3, pi] = gt._calc_t3(grid._nodes, panel)
    end
end
function _calc_tangents!(grid::gt.GridSurface2D, tangents)

    out, coor, lin, ndivscells, cin = gt.generate_getcellt_args!(grid)

    for pi in 1:grid.ncells
        coor[1] = pi
        panel = gt.get_cell_t!(out, grid, coor, lin, ndivscells)
        tangents[1, pi] = gt._calc_t1_2D(grid._nodes, panel)
        tangents[2, pi] = gt._calc_t2_2D(grid._nodes, panel)
    end
end
function _calc_tangents(grid::gt.AbstractGrid)
    tangents = zeros(3, grid.ncells)
    _calc_tangents!(grid, tangents)
    return tangents
end
_calc_tangents!(self::AbstractBody, tangents) = _calc_tangents!(self.grid, tangents)
_calc_tangents(self::AbstractBody) = _calc_tangents(self.grid)

"""
    calc_tangents!(body::AbstractBody, tangents::Matrix)

Calculates the tangent vector of every cell in `body` and stores them in the 3xN
matrix `tangents`.

**Output:** `tangents[:, i]` is the tangent vector of the i-th cell (linearly
indexed).
"""
const calc_tangents! = _calc_tangents!

"""
    calc_tangents(self::AbstractBody)

Calculates the tangent vector of every cell in `grid` returning a 3xN matrix.

See `calc_tangents!` documentation for more details.
"""
const calc_tangents = _calc_tangents


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
        obliques[1, pi] = gt._calc_o1(grid._nodes, panel)
        obliques[2, pi] = gt._calc_o2(grid._nodes, panel)
        obliques[3, pi] = gt._calc_o3(grid._nodes, panel)
    end
end
function _calc_obliques!(grid::gt.GridSurface2D, obliques)
    obliques .= 0
    obliques[3, :] .= 1
end
function _calc_obliques(grid::gt.AbstractGrid)
    obliques = zeros(3, grid.ncells)
    _calc_obliques!(grid, obliques)
    return obliques
end
_calc_obliques!(self::AbstractBody, obliques) = _calc_obliques!(self.grid, obliques)
_calc_obliques(self::AbstractBody) = _calc_obliques(self.grid)

"""
    calc_obliques!(body::AbstractBody, obliques::Matrix)

Calculates the oblique vector of every cell in `body` and stores them in the 3xN
matrix `obliques`.

**Output:** `obliques[:, i]` is the oblique vector of the i-th cell (linearly
indexed).
"""
const calc_obliques! = _calc_obliques!

"""
    calc_obliques(self::AbstractBody)

Calculates the oblique vector of every cell in `grid` returning a 3xN matrix.

See `calc_obliques!` documentation for more details.
"""
const calc_obliques = _calc_obliques


function _calc_areas!(grid::gt.GridTriangleSurface, areas)

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
        areas[pi] = gt._get_area(grid._nodes, panel)
    end
end
function _calc_areas(grid::gt.GridTriangleSurface)
    areas = zeros(grid.ncells)
    _calc_areas!(grid, areas)
    return areas
end
_calc_areas!(self::AbstractBody, areas) = _calc_areas!(self.grid, areas)
_calc_areas(self::AbstractBody) = _calc_areas(self.grid)

"""
    calc_areas!(body::AbstractBody, areas::Matrix)

Calculates the area of every cell in `body` and stores them in the array `areas`.

**Output:** `areas[i]` is the area of the i-th cell (linearly indexed).
"""
const calc_areas! = _calc_areas!

"""
    neighbor(body, ni::Int, ci::Int) -> ncoor
    neighbor(body, ni::Int, ccoor) -> ncoor

Returns the Cartesian coordinates `ncoor` of the `ni`-th neighbor of the cell
of linear indexing `ci` or coordinates `ccoor`.

```@example
# Calculate all normals
normals = pnl.calc_normals(body)

# Identify the second neighbor of the 10th cell
ncoor = pnl.neighbor(body, 2, 10)

# Convert Cartesian coordinates to linear indexing
ndivscells = Tuple(collect( 1:(d != 0 ? d : 1) for d in body.grid._ndivscells))
lin = LinearIndices(ndivscells)
ni = lin[ncoor...]

# Fetch the normal of such neighbor
normal = normals[:, ni]
```
"""
function neighbor(body::AbstractBody, args...; optargs...)
    gt.neighbor(body.grid, args...; optargs...)
end

function isedge(body::AbstractBody, args...)
    gt.isedge(body.grid, args...)
end



"""
    calc_areas(self::AbstractBody)

Calculates the area of every cell in `grid`, returning an array with all areas.

See `calc_areas!` documentation for more details.
"""
const calc_areas = _calc_areas

function _solvedflag(self::AbstractBody, val::Bool)
    # Remove all existing fields
    for field in Iterators.reverse(self.fields)
        remove_field(self, field)
    end

    # Add solved flag
    add_field(self, "solved", "scalar", [val], "system")
end

"Count the number of types in an Union type"
function _count(type::Type)
    if type isa Union
        return _count(type.a) + _count(type.b)
    elseif type isa Core.TypeofBottom
        return 0
    else
        return 1
    end
end
##### END OF ABSTRACT BODY #####################################################
