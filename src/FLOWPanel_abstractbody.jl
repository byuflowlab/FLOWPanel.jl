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
Abstract type `AbstractBody{N, E<:AbstractElement, TF<:Number}` where `N` is the number of
element types in this body and `E` is an Union containing the `N` element
types. `TF` is the floating point type used in this body.

Implementations of AbstractBody are expected to have the following fields
* `nodes::Matrix{TF}`                 : 3xnnodes matrix where `nodes[:, i]` is the position of the i-th node
* `vtk_cells::Vector{<:WriteVTK.MeshCell}` : Vector of WriteVTK cells
* `neighbor::Matrix{Int}`             : 3xncells matrix where `neighbor[i, j]` is the linear index of the cell neighboring the i-th edge of the j-th cell (or 0 if it's a boundary)
* `nnodes::Int`                       : Number of nodes
* `ncells::Int`                       : Number of cells
* `cells::Matrix{Int}`                : Cell connectivity (each column is a cell)
* `fields::Vector{String}`            : Available fields (solutions)
* `Oaxis::Matrix`                     : Coordinate system of body w.r.t global (3x3 matrix)
* `O::Vector`                         : Origin of body w.r.t. global (3-dim vector)
* `strength::Matrix`                  : Strength of each element of each type (ncells x N matrix)
* `velocity::Matrix{TF}`              : 3xncells apparent fluid velocity at control points (body frame)
* `velocity_kinematic::Matrix{TF}`    : 3xncells rigid-body kinematic velocity at control points (inertial frame)
* `potential::Vector{TF}`             : Total scalar potential at control points
* `CPoffset::Real`                    : Control point offset in normal direction
* `characteristiclength::Function`    : Function for computing the characteristic
                                        length of each panel used to offset each
                                        control point
* `kerneloffset::Real`                : Kernel offset to avoid singularities
* `kernelcutoff::Real`                : Kernel cutoff to avoid singularities
* `watertight::Bool`                  : Whether the body is watertight or not
* `Cps::Vector{TF}`                   : Pressure coefficient at each cell

and the following functions

```julia

    # Imposes boundary conditions to solve for element strengths.
    function solve(self::AbstractBody, Uinfs::Array{<:Real, 2}, solver::AbstractSolver, args...)
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
    # NOTE: `backend` is the N-body backend to be used, allowing support for
    #       different N-body methods (FMM, direct, etc).
    function _Uind(self::AbstractBody, targets::Array{<:Real, 2},
                    out::Array{<:Real, 2}, backend::AbstractBackend, args...; optargs...)
        .
        .
        .
    end

    # Returns the potential induced by the body on the targets `targets`. It
    # adds the potential at the i-th target to out[i].
    # NOTE: `backend` is the N-body backend to be used, allowing support for
    #       different N-body methods (FMM, direct, etc).
    function _phi(self::AbstractBody, targets::Array{<:Real, 2},
                    out::Array{<:Real, 1}, backend::AbstractBackend, args...; optargs...)
        .
        .
        .
    end
```
"""
abstract type AbstractBody{E<:AbstractElement, N, TF, DBC} end

function reset!(body::AbstractBody)
    body.velocity .= 0.0
    body.velocity_kinematic .= 0.0
    body.potential .= 0.0
    body.P .= 0.0
    body.F .= 0.0
    extra_reset!(body)
    return nothing
end

"""
    reset!(body::AbstractBody)

Clear per-step solution fields stored on `body` and delegate subtype-specific
reset work through `extra_reset!`.
"""
has_dirichlet_bc(::AbstractBody{<:Any, <:Any, <:Any, DBC}) where DBC = DBC

extra_reset!(body::AbstractBody) = nothing

"""
    solve(body::AbstractBody, Uinfs)

Solve a body boundary-value problem using the prescribed control-point
velocities `Uinfs`.
"""
function solve(self::AbstractBody, Uinfs::AbstractArray{<:Number, 2})
    error("solve(...) for body type $(typeof(self)) has not been implemented yet!")
end

"""
    trimesh2cells(mesh::VSPGeom.TriMesh)

Converts a `VSPGeom.TriMesh` to the `nodes` and `cells` matrices expected by
FLOWPanel body constructors.
"""
function trimesh2cells(mesh::VSPGeom.TriMesh)
    mesh_local = deepcopy(mesh)
    VSPGeom.setZeroBased!(mesh_local; value=false)

    nodes = reduce(hcat, mesh_local.points)
    cells = Int.(reduce(hcat, mesh_local.cells))

    return nodes, cells
end

"""
    _calc_edge_to_cells(cells::Matrix{Int})

Build a dictionary mapping each edge (as a sorted node-index pair) to the list
of `(cell_index, local_edge_index)` entries that share that edge.

Edge convention for a triangle with vertices `(n1, n2, n3) = cells[:, ci]`:
- Edge 1: `(n1, n2)`
- Edge 2: `(n2, n3)`
- Edge 3: `(n3, n1)`
"""
function _calc_edge_to_cells(cells::Matrix{Int})

    edge_to_cells = Dict{Tuple{Int,Int}, Vector{Tuple{Int,Int}}}()

    ncells = size(cells, 2)
    for ci in 1:ncells
        n1, n2, n3 = cells[1, ci], cells[2, ci], cells[3, ci]
        for (ei, (a, b)) in enumerate(((n1, n2), (n2, n3), (n3, n1)))
            key = a < b ? (a, b) : (b, a)
            list = get!(Vector{Tuple{Int,Int}}, edge_to_cells, key)
            push!(list, (ci, ei))
        end
    end

    return edge_to_cells
end

"""
    iswatertight(nodes, cells; return_open_cells=false) -> (watertight, open_cells)

Determine triangle-mesh watertightness from connectivity alone. A mesh is
watertight iff every undirected edge is referenced by exactly two cells.

When `return_open_cells=true`, the second return value contains the sorted,
unique 1-based cell indices incident to any edge whose incidence count is not
equal to two. Otherwise it is `Int[]`.
"""
function iswatertight(nodes::AbstractMatrix, cells::AbstractMatrix{<:Integer};
                      return_open_cells::Bool=false)
    edge_to_cells = _calc_edge_to_cells(Matrix{Int}(cells))
    watertight = true
    open_cells = Int[]

    for refs in values(edge_to_cells)
        if length(refs) != 2
            watertight = false
            if return_open_cells
                append!(open_cells, first.(refs))
            end
        end
    end

    if return_open_cells
        sort!(unique!(open_cells))
    end

    return watertight, open_cells
end

iswatertight(body::AbstractBody; return_open_cells::Bool=false) =
    iswatertight(body.nodes, body.cells; return_open_cells=return_open_cells)

"""
    calc_neighbors(cells::Matrix{Int})

Compute the `3 × ncells` neighbor matrix from triangle connectivity alone
(no GeometricTools grid required).

`neighbor[i, j]` is the linear index of the cell sharing edge `i` of cell `j`,
or `0` if that edge is on the boundary.  The edge convention matches
[`_calc_edge_to_cells`](@ref).
"""
function calc_neighbors(cells::Matrix{Int})

    ncells = size(cells, 2)
    neighbor = zeros(Int, 3, ncells)
    edge_to_cells = _calc_edge_to_cells(cells)

    for ci in 1:ncells
        n1, n2, n3 = cells[1, ci], cells[2, ci], cells[3, ci]
        for (ei, (a, b)) in enumerate(((n1, n2), (n2, n3), (n3, n1)))
            key = a < b ? (a, b) : (b, a)
            for (other_ci, _) in edge_to_cells[key]
                if other_ci != ci
                    neighbor[ei, ci] = other_ci
                    break
                end
            end
        end
    end

    return neighbor
end

function _shared_edge_direction(cells::AbstractMatrix{Int}, ci::Int, a::Int, b::Int)
    n1, n2, n3 = cells[1, ci], cells[2, ci], cells[3, ci]

    if (n1 == a && n2 == b) || (n2 == a && n3 == b) || (n3 == a && n1 == b)
        return 1
    elseif (n1 == b && n2 == a) || (n2 == b && n3 == a) || (n3 == b && n1 == a)
        return -1
    end

    error("Cell $ci does not contain shared edge ($a, $b).")
end

_flip_triangle!(cells::AbstractMatrix{Int}, ci::Int) = ((cells[2, ci], cells[3, ci]) = (cells[3, ci], cells[2, ci]); cells)

function _component_orientation_score(nodes::AbstractMatrix, cells::AbstractMatrix, component)
    node_ids = sort!(collect(Set(vec(cells[:, component]))))
    component_centroid = vec(sum(view(nodes, :, node_ids); dims=2)) ./ length(node_ids)
    score = zero(promote_type(eltype(nodes), Float64))

    for ci in component
        i1, i2, i3 = cells[1, ci], cells[2, ci], cells[3, ci]
        x1 = nodes[1, i1]; y1 = nodes[2, i1]; z1 = nodes[3, i1]
        x2 = nodes[1, i2]; y2 = nodes[2, i2]; z2 = nodes[3, i2]
        x3 = nodes[1, i3]; y3 = nodes[2, i3]; z3 = nodes[3, i3]

        nx = (y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1)
        ny = (z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1)
        nz = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)

        cx = (x1 + x2 + x3) / 3
        cy = (y1 + y2 + y3) / 3
        cz = (z1 + z2 + z3) / 3

        score += nx * (cx - component_centroid[1]) +
                 ny * (cy - component_centroid[2]) +
                 nz * (cz - component_centroid[3])
    end

    return score
end

function _winding_adjacency(cells::Matrix{Int})
    edge_to_cells = _calc_edge_to_cells(cells)
    adjacency = [Tuple{Int,Int,Int}[] for _ in 1:size(cells, 2)]

    for ((a, b), refs) in edge_to_cells
        if length(refs) >= 2
            c1 = refs[1][1]
            for k in 2:length(refs)
                c2 = refs[k][1]
                push!(adjacency[c1], (c2, a, b))
                push!(adjacency[c2], (c1, a, b))
            end
        end
    end

    return adjacency, edge_to_cells
end

"""
    ensure_consistent_winding!(cells, nodes; watertight=false, flip_normals=false)

Normalize triangle winding in `cells` so each connected component has a
consistent orientation. When `watertight=true`, each component is additionally
flipped so its normals point outward relative to the component centroid.

For open components, the first cell encountered in each component is preserved
as the seed orientation unless `flip_normals=true`.
"""
function ensure_consistent_winding!(cells::Matrix{Int}, nodes::AbstractMatrix;
                                    watertight::Bool=false,
                                    flip_normals::Bool=false)
    adjacency, _ = _winding_adjacency(cells)
    visited = falses(size(cells, 2))
    queue = Int[]

    for start in 1:size(cells, 2)
        visited[start] && continue

        empty!(queue)
        push!(queue, start)
        visited[start] = true

        component = Int[start]
        qidx = 1

        while qidx <= length(queue)
            ci = queue[qidx]
            qidx += 1

            for (cj, a, b) in adjacency[ci]
                if !visited[cj]
                    if _shared_edge_direction(cells, ci, a, b) == _shared_edge_direction(cells, cj, a, b)
                        _flip_triangle!(cells, cj)
                    end

                    visited[cj] = true
                    push!(queue, cj)
                    push!(component, cj)
                end
            end
        end

        if watertight && _component_orientation_score(nodes, cells, component) < 0
            for ci in component
                _flip_triangle!(cells, ci)
            end
        end

        if flip_normals
            for ci in component
                _flip_triangle!(cells, ci)
            end
        end
    end

    return cells
end

"""
    ensure_consistent_winding(nodes, cells; watertight=false, flip_normals=false)

Return a copy of `cells` with triangle winding normalized component-wise.
"""
function ensure_consistent_winding(nodes::AbstractMatrix, cells::Matrix{Int};
                                   watertight::Bool=false,
                                   flip_normals::Bool=false)
    return ensure_consistent_winding!(copy(cells), nodes;
                                      watertight=watertight,
                                      flip_normals=flip_normals)
end

##### COMMON FUNCTIONS  ########################################################


strength_names(::AbstractBody{ConstantSource, <:Any}) = ("sigma",)
strength_names(::AbstractBody{ConstantDoublet, <:Any}) = ("mu",)
strength_names(::AbstractBody{VortexRing, <:Any}) = ("gamma",)
strength_names(::AbstractBody{Union{ConstantSource, ConstantDoublet}, <:Any}) = ("sigma", "mu")
strength_names(::AbstractBody{Union{ConstantSource, VortexRing}, <:Any}) = ("sigma", "gamma")

"""
    write_vtk(name, body::AbstractBody, idx, t; overwrite=false)

Write the body mesh and its solution fields to VTK files at timestep `idx` with
simulation time `t`, following the same pattern as `write_vtk` for `PanelWake`.

Generates:
- `<name>.pvd`            — Paraview collection (append or create)
- `<name>_<idx>.vtm`      — VTK multiblock
- `<name>_<idx>_body.vtu` — unstructured triangular surface mesh with data

All fields stored on the body (strength, `Uinf`, `U`, `phi`, `P`, `Cps`,
`Gamma`, `F`) are written when non-empty.  Body-type-specific fields are added
via the internal `_write_vtk_body_fields!` hook so that subtypes (e.g.
`RigidWakeBody`) can contribute additional data without duplicating the common
code.

**Arguments**
- `name`      : Base filename (no extension)
- `body`      : Body instance to write
- `idx`       : Integer timestep index (embedded in filenames)
- `t`         : Simulation time (used as the PVD collection key)
- `overwrite` : Start a fresh PVD file when `true`; append when `false` (default)
"""
function write_vtk(name::String, body::AbstractBody, idx::Int=0, t::Real=0.0;
                   overwrite::Bool=false)

    # Route block files to a subdirectory named after the PVD
    _parent, _base = splitdir(name)
    subdir = joinpath(_parent, _base)
    mkpath(subdir)
    block_name = joinpath(subdir, _base)

    files = WriteVTK.paraview_collection(name; append=!overwrite) do pvd
        vtm = WriteVTK.vtk_multiblock(block_name * ".$idx.vtm")

        WriteVTK.vtk_grid(vtm, block_name * ".$(idx).vtu", body.nodes, body.vtk_cells) do vtk

            # --- Common solution fields ---

            # normals
            vtk["normals", VTKCellData()] = body.normals

            # Velocity potential  (ncells,)
            vtk["potential", VTKCellData()] = body.potential

            # Surface velocity  (3 × ncells)
            vtk["velocity", VTKCellData()] = body.velocity

            # Gauge pressure (ncells,)
            vtk["gauge pressure", VTKCellData()] = body.P

            # Distributed forces  (3 × ncells)
            vtk["F", VTKCellData()] = body.F

            # add strength fields
            for (i,name) in enumerate(strength_names(body))
                vtk[name, VTKCellData()] = view(body.strength, :, i)
            end

            # Body-type-specific fields (overload _write_vtk_body_fields! for subtypes)
            _write_vtk_body_fields!(vtk, body)
        end
        
        _write_vtk_other_fields!(vtm, block_name, body, idx)

        pvd[t] = vtm
    end

    return join(files, ", ")
end

# Default hook — no extra fields for generic AbstractBody
_write_vtk_body_fields!(vtk, ::AbstractBody) = nothing
_write_vtk_other_fields!(vtm, name, body::AbstractBody, idx) = nothing

function _vtk_stem(filename::AbstractString; path=nothing, num=nothing)
    stem = isnothing(num) ? filename : filename * ".$num"
    return isnothing(path) ? stem : joinpath(path, stem)
end

function _write_vtk_points_or_lines(filename::AbstractString, points;
                                    cells=nothing,
                                    point_data=(),
                                    cell_data=(),
                                    path=nothing,
                                    num=nothing,
                                    override_cell_type=nothing,
                                    optargs...)
    pts = points isa AbstractMatrix ? points : reduce(hcat, collect(points))
    vtk_cells = if isnothing(cells)
        [WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_VERTEX, [i]) for i in axes(pts, 2)]
    else
        celltype = isnothing(override_cell_type) ? WriteVTK.VTKCellTypes.VTK_POLY_LINE :
                   override_cell_type == 4 ? WriteVTK.VTKCellTypes.VTK_POLY_LINE :
                   override_cell_type
        [WriteVTK.MeshCell(celltype, collect(c) .+ 1) for c in cells]
    end

    stem = _vtk_stem(filename; path, num)
    saved = WriteVTK.vtk_grid(stem * ".vtu", pts, vtk_cells) do vtk
        for data in point_data
            vtk[data["field_name"], WriteVTK.VTKPointData()] = data["field_data"]
        end
        for data in cell_data
            vtk[data["field_name"], WriteVTK.VTKCellData()] = data["field_data"]
        end
    end

    files = saved isa AbstractVector{<:AbstractString} ? saved : WriteVTK.vtk_save(saved)
    return join(files, ", ")
end

"""
    get_ndivscells(body::AbstractBody)

Not supported for unstructured body variants.
"""
get_ndivscells(body::AbstractBody) = error("Not supported.")

"""
    get_ndivsnodes(body::AbstractBody)

Not supported for unstructured body variants.
"""
get_ndivsnodes(body::AbstractBody) = error("Not supported.")


"""
    get_cart2lin_cells(self::AbstractBody)

Not supported for unstructured body variants.
"""
function get_cart2lin_cells(self::AbstractBody)
    error("Not supported.")
end

"""
    get_cart2lin_nodes(self::AbstractBody)

Not supported for unstructured body variants.
"""
function get_cart2lin_nodes(self::AbstractBody)
    error("Not supported.")
end

get_strength(self::AbstractBody, i) = view(self.strength, i, :)
get_nstrengths(self::AbstractBody) = size(self.strength, 1)
set_strength(self::AbstractBody, i::Int, val) = (self.strength[i, :] .= val)

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
    M = rotation_matrix2(-roll, -pitch, -yaw)

    return rotatetranslate!(body, M, translation; reset_fields=reset_fields)
end

function rotation_matrix2(roll::Number, pitch::Number, yaw::Number)
    r = roll*pi/180
    p = pitch*pi/180
    y = yaw*pi/180
    Rx = [1 0 0; 0 cos(r) -sin(r); 0 sin(r) cos(r)]
    Ry = [cos(p) 0 sin(p); 0 1 0; -sin(p) 0 cos(p)]
    Rz = [cos(y) -sin(y) 0; sin(y) cos(y) 0; 0 0 1]
    return Rz * Ry * Rx
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

    # Translate back to origin
    for i in 1:body.nnodes
        body.nodes[1, i] -= body.O[1]
        body.nodes[2, i] -= body.O[2]
        body.nodes[3, i] -= body.O[3]
    end

    # Add translation to previous position
    body.O .+= T

    # Rotate and translate to new position
    for i in 1:body.nnodes
        v1 = M[1,1]*body.nodes[1,i] + M[1,2]*body.nodes[2,i] + M[1,3]*body.nodes[3,i]
        v2 = M[2,1]*body.nodes[1,i] + M[2,2]*body.nodes[2,i] + M[2,3]*body.nodes[3,i]
        v3 = M[3,1]*body.nodes[1,i] + M[3,2]*body.nodes[2,i] + M[3,3]*body.nodes[3,i]
        
        body.nodes[1, i] = v1 + body.O[1]
        body.nodes[2, i] = v2 + body.O[2]
        body.nodes[3, i] = v3 + body.O[3]
    end

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
function characteristiclength_sqrtarea(nodes, panel)
    i1, i2, i3 = panel[1], panel[2], panel[3]
    e1x = nodes[1, i2] - nodes[1, i1]
    e1y = nodes[2, i2] - nodes[2, i1]
    e1z = nodes[3, i2] - nodes[3, i1]
    e2x = nodes[1, i3] - nodes[1, i1]
    e2y = nodes[2, i3] - nodes[2, i1]
    e2z = nodes[3, i3] - nodes[3, i1]
    nx = e1y * e2z - e1z * e2y
    ny = e1z * e2x - e1x * e2z
    nz = e1x * e2y - e1y * e2x
    return sqrt(sqrt(nx*nx + ny*ny + nz*nz) * 0.5)
end

function calc_controlpoints!(nodes::AbstractMatrix, cells::AbstractMatrix,
                                controlpoints, normals; off::Real=0.005,
                                characteristiclength::Function=characteristiclength_sqrtarea)

    for pi in axes(cells, 2)
        i1, i2, i3 = cells[1, pi], cells[2, pi], cells[3, pi]

        # Centroid of triangle (average of vertices; equals centroid for triangles)
        controlpoints[1, pi] = (nodes[1, i1] + nodes[1, i2] + nodes[1, i3]) * 0.3333333333333333
        controlpoints[2, pi] = (nodes[2, i1] + nodes[2, i2] + nodes[2, i3]) * 0.3333333333333333
        controlpoints[3, pi] = (nodes[3, i1] + nodes[3, i2] + nodes[3, i3]) * 0.3333333333333333

        l = characteristiclength(nodes, view(cells, :, pi))

        # Offset the controlpoint in the normal direction
        controlpoints[1, pi] += off * l * normals[1, pi]
        controlpoints[2, pi] += off * l * normals[2, pi]
        controlpoints[3, pi] += off * l * normals[3, pi]
    end

    return controlpoints
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
function calc_controlpoints!(self::AbstractBody, normals=self.normals; 
        off=self.CPoffset, characteristiclength=self.characteristiclength)
    return calc_controlpoints!(self.nodes, self.cells, self.controlpoints, normals;
                                off, characteristiclength)
end

function calc_controlpoints(nodes::AbstractMatrix, cells::AbstractMatrix, normals;
        off::Real=0.005,
        characteristiclength::Function=characteristiclength_sqrtarea)
    controlpoints = zeros(promote_type(eltype(nodes), eltype(normals)), 3, size(cells, 2))
    calc_controlpoints!(nodes, cells, controlpoints, normals; off, characteristiclength)
    return controlpoints
end

function calc_controlpoints(self::AbstractBody, normals=self.normals;
        off=self.CPoffset,
        characteristiclength=self.characteristiclength)
    return calc_controlpoints(self.nodes, self.cells, normals; off, characteristiclength)
end

const _calc_controlpoints = calc_controlpoints

"""
    find_panels(body::AbstractBody, dim::Int, coord::Real; tol=1e-6)

Returns the indices of all panels (cells) whose control point is within `tol` of 
`coord` in the `dim` dimension.
"""
function find_panels(body::AbstractBody, dim::Int, coord::Number; tol=1e-6)
    CPs = body.controlpoints
    return findall(cp -> abs(cp[dim] - coord) < tol, eachcol(CPs))
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
function calc_normals!(nodes::AbstractMatrix, cells::AbstractMatrix, normals)
    for pi in axes(cells, 2)
        i1, i2, i3 = cells[1, pi], cells[2, pi], cells[3, pi]

        # Edge vectors from vertex 1
        e1x = nodes[1, i2] - nodes[1, i1]
        e1y = nodes[2, i2] - nodes[2, i1]
        e1z = nodes[3, i2] - nodes[3, i1]

        e2x = nodes[1, i3] - nodes[1, i1]
        e2y = nodes[2, i3] - nodes[2, i1]
        e2z = nodes[3, i3] - nodes[3, i1]

        # Cross product e1 × e2
        nx = e1y * e2z - e1z * e2y
        ny = e1z * e2x - e1x * e2z
        nz = e1x * e2y - e1y * e2x

        # Normalize
        len = sqrt(nx*nx + ny*ny + nz*nz)
        normals[1, pi] = nx / len
        normals[2, pi] = ny / len
        normals[3, pi] = nz / len
    end
end

function calc_normals!(self::AbstractBody, normals=self.normals; flipbyCPoffset=false)
    calc_normals!(self.nodes, self.cells, normals)
    if flipbyCPoffset
        normals .*= sign(self.CPoffset) != 0 ? sign(self.CPoffset) : 1
    end

    return normals
end

function calc_normals(nodes::AbstractMatrix, cells::AbstractMatrix)
    normals = zeros(eltype(nodes), 3, size(cells, 2))
    calc_normals!(nodes, cells, normals)
    return normals
end

calc_normals(self::AbstractBody) = calc_normals(self.nodes, self.cells)

const _calc_normals = calc_normals

function _calc_panel_tangent!(out, nodes::AbstractMatrix, panel)
    i1, i2 = panel[1], panel[2]
    tx = nodes[1, i2] - nodes[1, i1]
    ty = nodes[2, i2] - nodes[2, i1]
    tz = nodes[3, i2] - nodes[3, i1]
    len = sqrt(tx*tx + ty*ty + tz*tz)
    out[1] = tx / len
    out[2] = ty / len
    out[3] = tz / len
    return out
end

function _calc_panel_oblique!(out, nodes::AbstractMatrix, panel)
    tangent = (nodes[1, panel[2]] - nodes[1, panel[1]],
               nodes[2, panel[2]] - nodes[2, panel[1]],
               nodes[3, panel[2]] - nodes[3, panel[1]])
    edge2 = (nodes[1, panel[3]] - nodes[1, panel[1]],
             nodes[2, panel[3]] - nodes[2, panel[1]],
             nodes[3, panel[3]] - nodes[3, panel[1]])
    nx = tangent[2]*edge2[3] - tangent[3]*edge2[2]
    ny = tangent[3]*edge2[1] - tangent[1]*edge2[3]
    nz = tangent[1]*edge2[2] - tangent[2]*edge2[1]
    ox = ny*tangent[3] - nz*tangent[2]
    oy = nz*tangent[1] - nx*tangent[3]
    oz = nx*tangent[2] - ny*tangent[1]
    len = sqrt(ox*ox + oy*oy + oz*oz)
    out[1] = ox / len
    out[2] = oy / len
    out[3] = oz / len
    return out
end

function _calc_tangents!(nodes::AbstractMatrix, cells::AbstractMatrix, tangents)
    for pi in axes(cells, 2)
        panel = cells[:, pi]
        _calc_panel_tangent!(view(tangents, :, pi), nodes, panel)
    end
end
function _calc_tangents(nodes::AbstractMatrix, cells::AbstractMatrix)
    tangents = zeros(3, size(cells, 2))
    _calc_tangents!(nodes, cells, tangents)
    return tangents
end
_calc_tangents!(self::AbstractBody, tangents) = _calc_tangents!(self.nodes, self.cells, tangents)
_calc_tangents(self::AbstractBody) = _calc_tangents(self.nodes, self.cells)

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


function _calc_obliques!(nodes::AbstractMatrix, cells::AbstractMatrix, obliques)
    for pi in 1:size(cells, 2)
        panel = cells[:, pi]
        _calc_panel_oblique!(view(obliques, :, pi), nodes, panel)
    end
end
function _calc_obliques(nodes::AbstractMatrix, cells::AbstractMatrix)
    obliques = zeros(3, size(cells, 2))
    _calc_obliques!(nodes, cells, obliques)
    return obliques
end
_calc_obliques!(self::AbstractBody, obliques) = _calc_obliques!(self.nodes, self.cells, obliques)
_calc_obliques(self::AbstractBody) = _calc_obliques(self.nodes, self.cells)

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


"""
    calc_areas!(body::AbstractBody, areas)

Calculates the area of every cell in `body` and stores them in the vector `areas`.

**Output:** `areas[i]` is the area of the i-th cell (linearly indexed).
"""
function calc_areas!(nodes::AbstractMatrix, cells::AbstractMatrix, areas)
    for pi in 1:size(cells, 2)
        i1 = cells[1, pi]
        i2 = cells[2, pi]
        i3 = cells[3, pi]

        e1x = nodes[1, i2] - nodes[1, i1]
        e1y = nodes[2, i2] - nodes[2, i1]
        e1z = nodes[3, i2] - nodes[3, i1]

        e2x = nodes[1, i3] - nodes[1, i1]
        e2y = nodes[2, i3] - nodes[2, i1]
        e2z = nodes[3, i3] - nodes[3, i1]

        cx = e1y * e2z - e1z * e2y
        cy = e1z * e2x - e1x * e2z
        cz = e1x * e2y - e1y * e2x

        areas[pi] = 0.5 * sqrt(cx*cx + cy*cy + cz*cz)
    end
end
calc_areas!(self::AbstractBody, areas) = calc_areas!(self.nodes, self.cells, areas)

"""
    calc_areas(body::AbstractBody)

Calculates the area of every cell in `body` returning a vector of length N.

**Output:** `areas[i]` is the area of the i-th cell (linearly indexed).
"""
function calc_areas(nodes::AbstractMatrix, cells::AbstractMatrix)
    areas = zeros(size(cells, 2))
    calc_areas!(nodes, cells, areas)
    return areas
end
calc_areas(self::AbstractBody) = calc_areas(self.nodes, self.cells)

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
```"""
function neighbor(body::AbstractBody, ni, ci)
    return body.neighbor[ni, ci]
end

function get_cell(system::AbstractBody, i::Int)
    return system.cells[1,i], system.cells[2,i], system.cells[3,i]
end

function apply_freestream!(body::AbstractBody, uinf)
    eachcol(body.velocity) .+= Ref(uinf)
    for i in axes(body.controlpoints, 2)
        body.potential[i] += uinf[1] * body.controlpoints[1, i] +
                             uinf[2] * body.controlpoints[2, i] +
                             uinf[3] * body.controlpoints[3, i]
    end
    extra_apply_freestream!(body, uinf)
end

function apply_freestream!(bodies::Tuple, uinf)
    for body in bodies
        apply_freestream!(body, uinf)
    end
end

extra_apply_freestream!(body::AbstractBody, uinf) = nothing

#------- FastMultipole interface functions -------#

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, system::AbstractBody, i_body)
    
    # vertex indices for this panel
    i1, i2, i3 = get_cell(system, i_body)
    
    # extract vertices
    nodes = system.nodes
    v1x = nodes[1, i1]
    v1y = nodes[2, i1]
    v1z = nodes[3, i1]
    v2x = nodes[1, i2]
    v2y = nodes[2, i2]
    v2z = nodes[3, i2]
    v3x = nodes[1, i3]
    v3y = nodes[2, i3]
    v3z = nodes[3, i3]

    # normal
    dx1 = v2x - v1x
    dy1 = v2y - v1y
    dz1 = v2z - v1z
    dx2 = v3x - v1x
    dy2 = v3y - v1y
    dz2 = v3z - v1z
    # normal_x = dy1*dz2 - dz1*dy2
    # normal_y = dz1*dx2 - dx1*dz2
    # normal_z = dx1*dy2 - dy1*dx2
    # norm_inv = system.CPoffset / sqrt(normal_x * normal_x + normal_y * normal_y + normal_z * normal_z)
    # normal_x *= norm_inv
    # normal_y *= norm_inv
    # normal_z *= norm_inv

    # centroid
    cx = (v1x + v2x + v3x) * 0.3333333333333333 # + normal_x
    cy = (v1y + v2y + v3y) * 0.3333333333333333 # + normal_y
    cz = (v1z + v2z + v3z) * 0.3333333333333333 # + normal_z

    # get radius
    r = zero(eltype(buffer))
    dx = v1x - cx
    dy = v1y - cy
    dz = v1z - cz
    r = max(r, sqrt(dx*dx + dy*dy + dz*dz))
    dx = v2x - cx
    dy = v2y - cy
    dz = v2z - cz
    r = max(r, sqrt(dx*dx + dy*dy + dz*dz))
    dx = v3x - cx
    dy = v3y - cy
    dz = v3z - cz
    r = max(r, sqrt(dx*dx + dy*dy + dz*dz))

    # update buffer
    buffer[1, i_buffer] = cx
    buffer[2, i_buffer] = cy
    buffer[3, i_buffer] = cz
    buffer[4, i_buffer] = r

    # strength
    ns = size(system.strength, 2)
    for i in 1:ns
        buffer[4+i, i_buffer] = system.strength[i_body,i]
    end

    # save vertices
    buffer[4+ns+1, i_buffer] = v1x
    buffer[4+ns+2, i_buffer] = v1y
    buffer[4+ns+3, i_buffer] = v1z
    buffer[4+ns+4, i_buffer] = v2x
    buffer[4+ns+5, i_buffer] = v2y
    buffer[4+ns+6, i_buffer] = v2z
    buffer[4+ns+7, i_buffer] = v3x
    buffer[4+ns+8, i_buffer] = v3y
    buffer[4+ns+9, i_buffer] = v3z

    # additional fields added by dispatch
    additional_source_system_to_buffer!(buffer, i_buffer, system, i_body, 4+ns+9)

end

function additional_source_system_to_buffer!(buffer, i_buffer, system, i_body, i_start)
    # By default, no additional fields are added to the buffer. This function is
    # meant to be extended by dispatch for specific body types.
    return nothing
end

FastMultipole.numtype(system::AbstractBody) = eltype(system.strength)

FastMultipole.data_per_body(system::AbstractBody) = 4 + size(system.strength, 2) + 9 + additional_data_per_body(system)

additional_data_per_body(system::AbstractBody) = 0 # by default, no additional data is added to the buffer

function FastMultipole.get_position(system::AbstractBody, i)
    return system.controlpoints[1, i], system.controlpoints[2, i], system.controlpoints[3, i]
end

FastMultipole.strength_dims(system::AbstractBody) = size(system.strength, 2)

FastMultipole.get_n_bodies(system::AbstractBody) = system.ncells

function FastMultipole.buffer_to_target_system!(target_system::AbstractBody, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    throw("an <:AbstractBody cannot be used as a target system in FastMultipole calculations")
end

function FastMultipole.target_influence_to_buffer!(target_buffer, i_buffer, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_system::AbstractBody, i_target) where {PS,VS,GS}
    if PS
        target_buffer[4, i_buffer] = target_system.potential[i_target]
    end

    if VS
        vx, vy, vz = target_system.velocity[1, i_target], target_system.velocity[2, i_target], target_system.velocity[3, i_target]
        target_buffer[5, i_buffer] = vx
        target_buffer[6, i_buffer] = vy
        target_buffer[7, i_buffer] = vz
    end
end

function FastMultipole.direct!(target_system, target_index, derivatives_switch::FastMultipole.DerivativesSwitch{PS,GS,HS}, source_system::AbstractBody, source_buffer, source_index) where {PS,GS,HS}
    TF = eltype(target_system)
    for i_target in target_index # loop over targets
        target = FastMultipole.StaticArrays.SVector{3,TF}(target_system[1, i_target],
                  target_system[2, i_target],
                  target_system[3, i_target])
        
        phi_out = zero(eltype(target_system))
        U_out = @SVector zeros(eltype(target_system), 3)

        for i_source in source_index # loop over sources
            # evaluate influence due to this source
            phi, U, _ = induced(target, source_system, source_buffer, i_source, derivatives_switch; kerneloffset=source_system.kerneloffset)
            phi_out += phi
            U_out += U
        end

        # store results
        if PS
            target_system[4, i_target] += phi_out
        end
        if GS
            target_system[5, i_target] += U_out[1]
            target_system[6, i_target] += U_out[2]
            target_system[7, i_target] += U_out[3]
        end
    end
end

function FastMultipole.buffer_to_system_strength!(system::AbstractBody{<:Any,1,<:Any}, i_body, source_buffer, i_buffer)
    system.strength[i_body, 1] = source_buffer[5, i_buffer]
end

function FastMultipole.influence!(influence, target_buffer, source_system::AbstractBody, source_buffer)
    for i in 1:size(target_buffer, 2)
        v = FastMultipole.get_gradient(target_buffer, i)
        n = FastMultipole.get_normal(source_buffer, source_system, i)
        influence[i] = dot(v, n)
    end
end

##### END OF ABSTRACT BODY #####################################################
