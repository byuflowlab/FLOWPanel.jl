#=##############################################################################
# DESCRIPTION
    Utilities.
# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Sep 2018
  * License     : MIT License
=###############################################################################

dot(A, B) = sum(a*b for (a,b) in zip(A, B))
norm(A) = sqrt(mapreduce(x->x^2, +, A))
function cross!(out, A, B)
    out[1] = A[2]*B[3] - A[3]*B[2]
    out[2] = A[3]*B[1] - A[1]*B[3]
    out[3] = A[1]*B[2] - A[2]*B[1]
    return out
end
function cross(A::AbstractVector{T1}, B::AbstractVector{T2}) where {T1, T2}
    out = zeros(promote_type(T1, T2), 3)
    return cross!(out, A, B)
end
function cross(A::SVector{3,T1}, B::SVector{3,T2}) where {T1, T2}
    T = promote_type(T1, T2)
    return SVector{3,T}(A[2]*B[3] - A[3]*B[2],
                         A[3]*B[1] - A[1]*B[3],
                         A[1]*B[2] - A[2]*B[1])
end
mean(xs) = sum(xs)/length(xs)

"""
`direction(; alpha=0, beta=0)`

Return the unit vector of a freestream at `alpha` angle of attack (
rotation about -y direction) and `beta` sideslip angle (rotation about +z 
direction).

Both angles must be in degrees.
"""
function direction(; alpha::Number=0, beta::Number=0, gamma::Number=0)
    roll = gamma
    pitch = -alpha
    yaw = beta
    a, b, g = yaw*pi/180, pitch*pi/180, roll*pi/180
    Rz = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1]
    Ry = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)]
    Rx = [1 0 0; 0 cos(g) -sin(g); 0 sin(g) cos(g)]
    return (Rz * Ry * Rx)[:, 1]
end

"""
`simplewing(b, ar, tr, twist_root, twist_tip, lambda, gamma;
bodytype=RigidWakeBody,
span_NDIVS="automatic", rfl_NDIVS="automatic",
airfoil_root="naca6412.dat", airfoil_tip="naca6412.dat",
airfoil_path=def_rfl_path)`

Generate a symmetric single-section wing and return it as the requested body
type.
"""
function simplewing(b::Number, ar::Number, tr::Number, twist_root::Number,
                      twist_tip::Number, lambda::Number, gamma::Number;
                      bodytype::Type{<:AbstractBody}=RigidWakeBody{Union{ConstantSource,ConstantDoublet}},
                      span_NDIVS::ndivstype=nothing,
                      rfl_NDIVS::ndivstype=nothing,
                      airfoil_root::String="naca6412.dat",
                      airfoil_tip::String="naca6412.dat",
                      airfoil_path::String=def_rfl_path,
                      spl_s::Real=0.0000001,
                      rflspl_s::Real=0.00000001,
                      verify_spline::Bool=true,
                      verify_rflspline::Bool=true,
                      b_low=-1.0, b_up=1.0,
                      opt_args...
                      )
    error("simplewing has moved to examples/helper_functions.jl because it depends on GeometricTools. Include that helper file and call simplewing from the example namespace.")
end

direction(dir; X0=zeros(3)) = X -> dot(X - X0, dir)
loop(; Oaxis=Matrix(1.0I, 3, 3), X0=zeros(3)) = X -> -atand(dot(X - X0, Oaxis[:, 2]), dot(X - X0, Oaxis[:, 3]))
nojunction(X) = Inf

noprocessmsh(name, msh, options) = msh

"""
    meshes2nodes_cells(mesh::Meshes.SimpleMesh) -> (nodes, cells)

Extract node and cell matrices from a `Meshes.SimpleMesh`.

Returns `nodes` (3×N Float64 matrix of vertex coordinates) and `cells`
(3×M Int matrix of triangle connectivity, 1-based vertex indices).
"""
function meshes2nodes_cells(mesh::Meshes.SimpleMesh)
    nnodes = length(mesh.vertices)
    ncells = length(mesh.topology.connec)
    nodes = zeros(Float64, 3, nnodes)
    cells = zeros(Int, 3, ncells)
    for (i, v) in enumerate(mesh.vertices)
        nodes[:, i] .= v.coords
    end
    for (i, c) in enumerate(mesh.topology.connec)
        cells[:, i] .= c.indices
    end
    return nodes, cells
end

"""
    read_gmsh(path) -> Meshes.SimpleMesh

Read an ASCII Gmsh 4.1 mesh containing first-order line and/or triangle
elements. This intentionally small reader covers FLOWPanel's tracked surface
and trailing-edge meshes without pulling GeoIO's geospatial dependency stack
into the package environment.

Node tags may be sparse or unordered. Gmsh point elements (type 15) are
ignored. If both line and triangle blocks are present, triangles are returned;
line connectivity is returned only for a line-only mesh. Other element types
are rejected explicitly.
"""
mutable struct _GmshTokenCursor
    data::String
    cursor::Int
    stop::Int
    path::String
end

function _gmsh_section_cursor(data::String, path::AbstractString, name::String)
    section_start = findnext("\$" * name, data, 1)
    section_start === nothing && error("Missing or malformed \$$(name) section in $(path)")
    section_end = findnext("\$End" * name, data, last(section_start) + 1)
    section_end === nothing && error("Missing or malformed \$$(name) section in $(path)")
    return _GmshTokenCursor(data, last(section_start) + 1,
                            first(section_end) - 1, String(path))
end

@inline function _gmsh_skip_whitespace!(tokens::_GmshTokenCursor)
    i = tokens.cursor
    while i <= tokens.stop && codeunit(tokens.data, i) <= 0x20
        i += 1
    end
    tokens.cursor = i
    return i
end

@inline function _gmsh_next_token!(tokens::_GmshTokenCursor)
    start = _gmsh_skip_whitespace!(tokens)
    start <= tokens.stop || error("Unexpected end of Gmsh section in $(tokens.path)")
    i = start
    while i <= tokens.stop && codeunit(tokens.data, i) > 0x20
        i += 1
    end
    tokens.cursor = i
    return SubString(tokens.data, start, i - 1)
end

@inline function _gmsh_next_int!(tokens::_GmshTokenCursor)
    i = _gmsh_skip_whitespace!(tokens)
    i <= tokens.stop || error("Unexpected end of Gmsh section in $(tokens.path)")
    negative = codeunit(tokens.data, i) == UInt8('-')
    i += negative
    value = 0
    found_digit = false
    while i <= tokens.stop
        byte = codeunit(tokens.data, i)
        UInt8('0') <= byte <= UInt8('9') || break
        value = 10value + Int(byte - UInt8('0'))
        found_digit = true
        i += 1
    end
    found_digit || error("Malformed integer in Gmsh file $(tokens.path)")
    tokens.cursor = i
    return negative ? -value : value
end

@inline _gmsh_next_float!(tokens::_GmshTokenCursor) =
    parse(Float64, _gmsh_next_token!(tokens))

function read_gmsh(path::AbstractString)
    data = read(path, String)

    format_tokens = try
        _gmsh_section_cursor(data, path, "MeshFormat")
    catch err
        err isa ErrorException && occursin("Missing or malformed", err.msg) &&
            error("Not a Gmsh file: $(path)")
        rethrow()
    end
    version = _gmsh_next_token!(format_tokens)
    version == "4.1" ||
        error("FLOWPanel.read_gmsh supports Gmsh 4.1; $(path) uses $(version)")
    _gmsh_next_int!(format_tokens) == 0 ||
        error("FLOWPanel.read_gmsh supports ASCII Gmsh files only: $(path)")
    _gmsh_next_int!(format_tokens) # data size

    node_tokens = _gmsh_section_cursor(data, path, "Nodes")
    nblocks = _gmsh_next_int!(node_tokens)
    nnodes = _gmsh_next_int!(node_tokens)
    minimum_tag = _gmsh_next_int!(node_tokens)
    maximum_tag = _gmsh_next_int!(node_tokens)
    coordinates = Vector{NTuple{3,Float64}}(undef, nnodes)

    # A dense lookup is substantially faster for the common Gmsh case while
    # the dictionary fallback avoids allocating according to a huge sparse tag.
    dense_tags = minimum_tag >= 0 && maximum_tag >= minimum_tag &&
                 maximum_tag - minimum_tag <= 8nnodes
    dense_tag_to_index = dense_tags ? zeros(Int, maximum_tag - minimum_tag + 1) : Int[]
    sparse_tag_to_index = dense_tags ? Dict{Int,Int}() : sizehint!(Dict{Int,Int}(), nnodes)

    next_index = 1
    for _ in 1:nblocks
        entity_dim = _gmsh_next_int!(node_tokens)
        _gmsh_next_int!(node_tokens) # entity tag
        parametric = _gmsh_next_int!(node_tokens)
        block_nodes = _gmsh_next_int!(node_tokens)
        tags = Vector{Int}(undef, block_nodes)
        for i in eachindex(tags)
            tags[i] = _gmsh_next_int!(node_tokens)
        end
        for tag in tags
            x = _gmsh_next_float!(node_tokens)
            y = _gmsh_next_float!(node_tokens)
            z = _gmsh_next_float!(node_tokens)
            coordinates[next_index] = (x, y, z)
            for _ in 1:(parametric != 0 ? entity_dim : 0)
                _gmsh_next_float!(node_tokens)
            end
            if dense_tags
                dense_tag_to_index[tag - minimum_tag + 1] = next_index
            else
                sparse_tag_to_index[tag] = next_index
            end
            next_index += 1
        end
    end
    next_index == nnodes + 1 || error("Gmsh node count mismatch in $(path)")

    @inline function node_index(tag)
        if dense_tags
            lookup_index = tag - minimum_tag + 1
            index = checkbounds(Bool, dense_tag_to_index, lookup_index) ?
                    @inbounds(dense_tag_to_index[lookup_index]) : 0
        else
            index = get(sparse_tag_to_index, tag, 0)
        end
        index != 0 || error("Element references unknown node tag $(tag) in $(path)")
        return index
    end

    element_tokens = _gmsh_section_cursor(data, path, "Elements")
    element_blocks = _gmsh_next_int!(element_tokens)
    _gmsh_next_int!(element_tokens) # total element count (includes point elements)
    _gmsh_next_int!(element_tokens) # minimum element tag
    _gmsh_next_int!(element_tokens) # maximum element tag
    line_connectivity = Meshes.Connectivity{Meshes.Segment,2}[]
    triangle_connectivity = Meshes.Connectivity{Meshes.Triangle,3}[]
    for _ in 1:element_blocks
        _gmsh_next_int!(element_tokens) # entity dimension
        _gmsh_next_int!(element_tokens) # entity tag
        element_type = _gmsh_next_int!(element_tokens)
        block_elements = _gmsh_next_int!(element_tokens)
        element_type in (1, 2, 15) ||
            error("Unsupported Gmsh element type $(element_type) in $(path)")
        if element_type == 1
            sizehint!(line_connectivity, length(line_connectivity) + block_elements)
            for _ in 1:block_elements
                _gmsh_next_int!(element_tokens) # element tag
                i = node_index(_gmsh_next_int!(element_tokens))
                j = node_index(_gmsh_next_int!(element_tokens))
                push!(line_connectivity, Meshes.connect((i, j)))
            end
        elseif element_type == 2
            sizehint!(triangle_connectivity, length(triangle_connectivity) + block_elements)
            for _ in 1:block_elements
                _gmsh_next_int!(element_tokens) # element tag
                i = node_index(_gmsh_next_int!(element_tokens))
                j = node_index(_gmsh_next_int!(element_tokens))
                k = node_index(_gmsh_next_int!(element_tokens))
                push!(triangle_connectivity, Meshes.connect((i, j, k)))
            end
        else
            for _ in 1:block_elements
                _gmsh_next_int!(element_tokens) # element tag
                _gmsh_next_int!(element_tokens) # node tag
            end
        end
    end
    connectivity = isempty(triangle_connectivity) ? line_connectivity : triangle_connectivity
    isempty(connectivity) && error("No line or triangle elements found in $(path)")
    return Meshes.SimpleMesh(coordinates, connectivity)
end


function distancetoline(line::Matrix; symmetry=nothing)
    X0 = view(line, :, 1)
    X1 = view(line, :, size(line, 2))

    if isnothing(symmetry)
        return distancetoline(X0, X1)(base)
    else
        X0sym = X0 - 2*dot(X0, symmetry)*symmetry
        X1sym = X1 - 2*dot(X1, symmetry)*symmetry
        fun1 = distancetoline(X0, X1)
        fun2 = distancetoline(X0sym, X1sym)

        return (args...; optargs...) -> min(fun1(args...; optargs...), fun2(args...; optargs...))
    end

end

function distancetoline_symmetric(symmetry)
    return (args...; optargs...) -> distancetoline(args...; optargs..., symmetry=symmetry)
end

function distancetoline(X0::AbstractVector, X1::AbstractVector)

    # Calculate unit direction of line and length
    dir = X1 - X0
    len = norm(dir)
    dir /= len

    function calc_distancetoline(X)

        X2 = X - X0

        # Projection to line
        proj = dot(X2, dir)
        X3 = proj * dir

        # Component normal to line
        Xn = X2 - X3

        # Distance to line
        distance = norm(Xn)

        # Add projection to line if point is beyond the line
        if proj < 0
            distance += abs(proj)
        elseif proj > len
            distance += proj - len
        end

        return distance
    end

    return calc_distancetoline
end


"""
    filter_splitsurfaces(connectivity, vertices, controlsurfaces; offset=zeros(3), invrotation=one(Meshes.QuatRotation), scaling=1.0)

Predicate used during mesh preprocessing to decide whether an element should be
split along a control-surface boundary.
"""
function filter_splitsurfaces(connectivity, vertices, controlsurfaces;
                                        offset=zeros(3),
                                        invrotation::Meshes.Rotation=one(Meshes.QuatRotation),
                                        scaling=1.0)

    # Fetch and untransform the vertices
    Xs = [invrotation * (vertices[vi].coords / scaling) - offset for vi in connectivity.indices]

    # Calculate centroid
    X = mean(Xs)

    # Identify whether this element is intercepted by the edge of a control surface
    for (csi, controlsurface) in enumerate(controlsurfaces)    # Iterate over control surfaces

        (; hinge, side1, side2, boundingbox, tol) = controlsurface

        # Identify whether the centroid is aft the hinge (positive side)
        hingecrit = dot(X - hinge.center, hinge.normal) > tol

        # Identify whether the vertices are split by side 1
        side1crit = length(unique( dot(X - side1.center, side1.normal) > tol for X in Xs )) != 1

        # Identify whether the vertices are split by side 2
        side2crit = length(unique( dot(X - side2.center, side2.normal) > tol for X in Xs )) != 1

        # Identify whether the centroid is inside the bounding box
        boxcrit = all( boundingbox.lower[i] <= X[i] <= boundingbox.upper[i] for i in 1:length(X) )

        # Check for mirrored conditions if requested
        mirror = (controlsurface.mirror==:symmetric || controlsurface.mirror==:antisymmetric)

        # Check whether the mirroring conditions are satisfied
        if !(hingecrit && (side1crit || side2crit) && boxcrit) && mirror

            # Bring vertices from the other side of the symmetry plane to this side
            Xs = [X .* Meshes.Vec(1, -1, 1) for X in Xs]
            X = mean(Xs)

            # Re-evaluation criteria
            hingecrit = dot(X - hinge.center, hinge.normal) > tol
            side1crit = length(unique( dot(X - side1.center, side1.normal) > tol for X in Xs )) != 1
            side2crit = length(unique( dot(X - side2.center, side2.normal) > tol for X in Xs )) != 1
            boxcrit = all( boundingbox.lower[i] <= X[i] <= boundingbox.upper[i] for i in 1:length(X) )

        end

        # Return if identified that element is intercepted
        if hingecrit && (side1crit || side2crit) && boxcrit
            return false
        end

    end

    return true
end


"""
    `find_i(body, xtarget::Number, gdim::Int, xdim::Int; xdir=nothing)`

Find the row or column of cells in structured grid `body` that is the
closest to `xtarget` in the `xdim` spatial dimension. Use `gdim=1` to obtain
a row, and `gdim=2` to obtain a column.

Alternatively, use an arbitrary direction `xdir` in place of `xdim`, if given.

Returns `(itarget, pos, errmin, lin)` where `itarget` is the index of the best
candidate row/column, `pos` is the position of this row/column projected to
`xdim` or `xdir`, and `errmin` is the error between `pos` and `xtarget`. `lin`
is the `LinearIndices` for the user to iterate over the row/column as
`lin[itarget, j]` if `gdim==1`, or `lin[j, itarget]` if `gdim==2`.

"""
function find_i(body::Union{NonLiftingBody, AbstractLiftingBody}, controlpoints,
                xtarget::Number, gdim::Int, xdim::Int; xdir=nothing,
                filter=(args...)->true
               )

    # Define direction on which point will be projected to
    _xdir = xdir!=nothing ? xdir : (i==xdim for i in 1:3)
    _xdir ./= norm(_xdir)

    # Define dimension on which cells will be iterated
    gdimite =   gdim==1 ? 2 :
                gdim==2 ? 1 :
                error("Invalid dimension $gdim; expected 1 or 2.")

    gdims = get_ndivscells(body)            # Grid dimensions
    lin = get_linearindex(body, gdims)      # Linear indexing

    pos = Inf                               # Position of closest candidate
    errmin = Inf                            # Error of closest candidate
    itarget = -1                            # index of closest candidate
    Xmean = zeros(3)                        # Average position along row/column

    for i in 1:gdims[gdim]                  # Iterate over row/column

        # Calculate the average control point over this row/column
        Xmean .= 0
        nmean = 0
        for j in 1:gdims[gdimite] # Iterate over panels in this row/column
            indices = gdim==1 ? (i, j, 1) : (j, i, 1)
            linind = lin[indices...]
            point = view(controlpoints, :, linind)
            if filter(point, linind)
                Xmean .+= point
                nmean += 1
            end
        end

        if nmean != 0
            Xmean ./= nmean

            # Determine proximity of average control point to target position
            err = abs(dot(Xmean, _xdir) - xtarget)

            # Compare this candidate
            if err <= errmin
                pos = dot(Xmean, _xdir)
                errmin = err
                itarget = i
            end
        end
    end

    if itarget==-1
        error("Logic error: no slice found! (itarget=$(itarget))")
    end

    return itarget, pos, errmin, lin
end

"""
    get_linearindex(body::Union{NonLiftingBody, AbstractLiftingBody})

Return the LinearIndex of the grid of `body`.
"""
function get_linearindex(body::Union{NonLiftingBody, AbstractLiftingBody})
    gdims = get_ndivscells(body)            # Grid dimensions
    return get_linearindex(body, gdims), gdims
end

function get_linearindex(body::Union{NonLiftingBody, AbstractLiftingBody}, gdims)
    ndivscells = Tuple(n + 1*(n==0) for n in gdims) # n=0 -> n=1 for quasi-dimensions
    lin = LinearIndices(ndivscells)         # Linear indexing
    return lin
end

"""
    decompose!(out::Matrix, V::Matrix, ihat::Vector, jhat::Vector, khat::Vector)

Project each column of `V` onto the orthonormal bases `ihat`, `jhat`, `khat`.
The projection is calculated
"""
function decompose!(out::AbstractMatrix, V::AbstractMatrix,
                    ihat::AbstractVector, jhat::AbstractVector,
                    khat::AbstractVector)

    # Error cases
    @assert size(out, 1)==size(V, 1) && size(out, 2)==size(V, 2) ""*
        "Invalid `out` matrix. Expected size $(size(V)); got $(size(out))."
    @assert abs(norm(ihat) - 1) <= 2*eps() ""*
        "ihat=$(ihat) is not a unitary vector"
    @assert abs(norm(jhat) - 1) <= 2*eps() ""*
        "jhat=$(jhat) is not a unitary vector"
    @assert abs(norm(khat) - 1) <= 2*eps() ""*
        "khat=$(khat) is not a unitary vector"

    # Project each column into the ihat, jhat, khat bases
    for (o, v) in zip(eachcol(out), eachcol(V))
        o[1] = dot(v, ihat)
        o[2] = dot(v, jhat)
        o[3] = dot(v, khat)
    end

    return out
end

"""
    decompose!(out, V, ihat, jhat)

Similar to `decompose!(out, V, ihat, jhat, khat)`, but automatically calculates
`khat` from `ihat` and `jhat`.
"""
decompose!(out, V, ihat, jhat) = decompose!(out, V, ihat, jhat, cross(ihat, jhat))

"""
    decompose(V, ihat, jhat)

Similar to `decompose!(out, V, ihat, jhat)` but without calculating the
projection in-place.
"""
function decompose(V::AbstractMatrix{T1},
                    ihat::AbstractVector{T2}, jhat::AbstractVector{T3}
                    ) where {T1, T2, T3}
    return decompose!(similar(V, promote_type(T1, T2, T3)), V, ihat, jhat)
end


"""
    slicefield(body::AbstractBody, fieldname::String,
                position::Number, direction::Vector, row::Bool)

Return a slice of the field `fieldname` in `body` corresponding to the row or
column ("row" is the first dimension of the grid, "column" is the second
dimension) that is the closest to `position` calculated as the projection
of the average cell position in the direction `direction`.

**Example:** For a wing with its span aligned along the y-axis, the pressure
along a slice of the wing at the spanwise position y=0.5 is obtained as
`slicefield(wing, "Cp", 0.5, [0, 1, 0], false)`.
"""
function slicefield(body::AbstractBody, fieldname::String, args...; optargs...)

    normals = calc_normals(body)
    controlpoints = calc_controlpoints(body, normals)

    return slicefield(body, controlpoints, fieldname, args...; optargs...)
end

"""
    slicefield(body::AbstractBody, controlpoints::Matrix,
                    fieldname::String,
                    position::Number, direction::Vector, row::Bool)

Same thing, but with the option of providing the control points as to save
avoid memory allocation.
"""
function slicefield(body::AbstractBody, controlpoints::Arr,
                    fieldname::String,
                    position::Number, direction::Vector, row::Bool;
                    reduce=true, filter=(args...)->true
                    ) where {Arr<:AbstractArray{<:Number,2}}

    # Fetch field
    field = getfield(body, Symbol(fieldname))

    # Find index of row or column slicing the field
    gdim = row ? 1 : 2                          # Dimension to slice
    islice, pos, errmin, lin = find_i(body, controlpoints, position, gdim, -1;
                                                xdir=direction, filter=filter)

    # Slice field
    ncell = size(controlpoints, 2)
    indices = collect(row==1 ? lin[islice, j] : lin[j, islice] for j in 1:ncell)

    slice = field[indices]

    # Points along the slice
    slicepoints = controlpoints[:, indices]

    # Reduce the implicit double-column of the triangular grid into one column
    if reduce
        slice = [ (slice[2*(i-1) + 1]+slice[2*(i-1) + 2])/2 for i in 1:Int(length(slice)/2)]
        slicepoints = [ (slicepoints[i, 2*(j-1) + 1]+slicepoints[i, 2*(j-1) + 2])/2 for i in 1:size(slicepoints,1), j in 1:Int(size(slicepoints, 2)/2)]
    end

    return slicepoints, slice
end

function slice_scalarfield(body::AbstractBody, fieldname::Symbol, dim::Int, coord::Number, tol::Number; filter=(args...)->true)
    # find the index of all panels within tolerance
    controlpoints = body.controlpoints
    target = fill(coord, size(controlpoints, 2)) .- view(controlpoints, dim, :)
    target .*= target
    target .= sqrt.(target)

    # Find which panels are within tolerance
    idx = (target .<= tol) .&& filter.(eachcol(controlpoints), eachindex(target))
    field_values = getfield(body, fieldname)[idx]

    if sum(idx) == 0
        error("No panels found within tolerance $(tol) of coordinate $(coord) in dimension $(dim).")
    end

    return controlpoints[:, idx], field_values
end

"""
    `calcfield_winding(body::Union{NonLiftingBody, AbstractLiftingBody})`

Calculate the winding number of a body with a Meshes.jl mesh object. It adds
the winding number of each cell (associated to its corresponding control point)
as the field `"winding"`.
"""
function calcfield_winding(body::Union{NonLiftingBody, AbstractLiftingBody};
                            fieldname="winding", addfield=true)

    # Precalculations
    normals = _calc_normals(body)
    controlpoints = _calc_controlpoints(body, normals)

    vertices = [Meshes.Point(nodes...) for nodes in eachcol(body.nodes)]
    triangles = [Meshes.connect(Tuple(body.cells[:, i])) for i in axes(body.cells, 2)]
    mesh = Meshes.SimpleMesh(vertices, triangles)

    # Calculate winding number associated to control points
    windings = calcfield_winding(mesh, controlpoints)

    # Save field in body
    if addfield
        add_field(body, fieldname, "scalar", windings, "cell")
    end

    return windings
end

function calcfield_winding(mesh::Meshes.Mesh, controlpoints::AbstractMatrix)

    points = ( Meshes.Point(Meshes.Vec(p...)) for p in eachcol(controlpoints) )
    windings = Meshes.winding(points, mesh)

    return windings
end


"""
    `calc_minmax_winding(body::Union{NonLiftingBody, AbstractLiftingBody}) ->
(minw, maxw)`

Calculates the winding number of each cell (associated to its corresponding
control point) in a body with a Meshes.jl mesh object, and returns both minimum
and maximum values.
"""
function calc_minmax_winding(args...; optargs...)
    windings = calcfield_winding(args...; optargs...)

    return minimum(windings), maximum(windings)
end
