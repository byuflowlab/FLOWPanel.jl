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
abstract type AbstractBody{E<:AbstractElement, N, TF} end

function reset!(body::AbstractBody)
    body.velocity .= 0.0
    for vte in body.velocity_te
        vte .= 0.0
    end
    body.potential .= 0.0
    body.Cp .= 0.0
    body.F .= 0.0
    return nothing
end

"""
    `solve(body::AbstractBody, Uinfs::Array{<:Real, 2})`

Impose boundary conditions to solve for element strengths. `Uinds[:, i]` is the
velocity at the i-th control point used in the boundary condition.
"""
function solve(self::AbstractBody, Uinfs::AbstractArray{<:Number, 2})
    error("solve(...) for body type $(typeof(self)) has not been implemented yet!")
end

"""
    `grid2cells(grid::GeometricTools.GridTriangleSurface)`

Converts the cells in a `GeometricTools` continuous grid to a `Matrix{Int}` of size `3 x ncells`.
"""
function grid2cells(grid::gt.GridTriangleSurface)
    cells = zeros(Int, 3, grid.ncells)
    for i in 1:grid.ncells
        cells[:, i] .= gt.get_cell(grid, i)
    end
    return cells
end

##### COMMON FUNCTIONS  ########################################################

"""
    Uind!(self::AbstractBody, targets, out, args...; optargs...)

Returns the velocity induced by the body on the targets `targets`, which is a
3xn matrix. It adds the velocity at the i-th target to `out[:, i]`.
"""
function Uind!(self::AbstractBody, targets, out, backend::AbstractBackend, args...; optargs...)

    # ERROR CASES
    if check_solved(self)==false
        error("Body hasn't been solved yet."*
              " Please call `solve()` function first.")
    end

    _Uind!(self, targets, out, backend::AbstractBackend, args...; optargs...)
end

"""
    phi!(self::AbstractBody, targets, out, args...; optargs...)

Returns the potential induced by the body on the targets `targets`. It adds the
potential at the i-th target to `out[:, i]`.
"""
function phi!(self::AbstractBody, targets, out, backend::AbstractBackend, args...; optargs...)

    # ERROR CASES
    if check_solved(self)==false
        error("Body hasn't been solved yet."*
              " Please call `solve()` function first.")
    end

    _phi!(self, targets, out, backend::AbstractBackend, args...; optargs...)
end

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

All fields stored on the body (strength, `Uinf`, `U`, `phi`, `Cp`, `Cps`,
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
function write_vtk(name::String, body::AbstractBody, idx::Int, t::Real;
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

            # Pressure coefficient (ncells,)
            vtk["pressure coefficient", VTKCellData()] = body.Cp

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

"""
    check_solved(self::AbstractBody)

Returns `true` if the body has been solved. Returns false otherwise.
"""
check_solved(self::AbstractBody) = self.solved

"""
    _solvedflag(self::AbstractBody, val::Bool)

Sets the `solved` flag of the body.
"""
_solvedflag(self::AbstractBody, val::Bool) = self.solved = val




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

function _calc_tangents!(nodes::AbstractMatrix, cells::AbstractMatrix, tangents)
    for pi in axes(cells, 2)
        panel = cells[:, pi]
        tangents[1, pi] = gt._calc_t1(nodes, panel)
        tangents[2, pi] = gt._calc_t2(nodes, panel)
        tangents[3, pi] = gt._calc_t3(nodes, panel)
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
        obliques[1, pi] = gt._calc_o1(nodes, panel)
        obliques[2, pi] = gt._calc_o2(nodes, panel)
        obliques[3, pi] = gt._calc_o3(nodes, panel)
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
    for i in eachindex(body.velocity_te)
        eachcol(body.velocity_te[i]) .+= Ref(uinf)
    end
end


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

    # vertex indices for this panel
    i1, i2, i3 = get_cell(system, i)
    
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

    # get normal vector
    dx1 = v2x - v1x
    dy1 = v2y - v1y
    dz1 = v2z - v1z
    dx2 = v3x - v1x
    dy2 = v3y - v1y
    dz2 = v3z - v1z
    normal_x = dy1*dz2 - dz1*dy2
    normal_y = dz1*dx2 - dx1*dz2
    normal_z = dx1*dy2 - dy1*dx2
    norm_inv = system.CPoffset / sqrt(normal_x * normal_x + normal_y * normal_y + normal_z * normal_z)
    normal_x *= norm_inv
    normal_y *= norm_inv
    normal_z *= norm_inv

    # get centroid
    cx = (v1x + v2x + v3x) * 0.3333333333333333 + normal_x
    cy = (v1y + v2y + v3y) * 0.3333333333333333 + normal_y
    cz = (v1z + v2z + v3z) * 0.3333333333333333 + normal_z

    return cx, cy, cz
end

FastMultipole.strength_dims(system::AbstractBody) = size(system.strength, 2)

FastMultipole.get_n_bodies(system::AbstractBody) = system.ncells

function FastMultipole.buffer_to_target_system!(target_system::AbstractBody, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    throw("an <:AbstractBody cannot be used as a target system in FastMultipole calculations")
end

# function FastMultipole.direct!(target_system, target_index, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, source_system::AbstractBody, source_buffer, source_index) where {PS,VS,GS}
#     throw("FastMultipole.direct! is not implemented for `<:AbstractBody` systems of type $(typeof(source_system))")
# end

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

# function FastMultipole.strength_to_value(strength, source_system::AbstractBody{<:Any, 1, <:Any})
#     return strength[1]
# end

function FastMultipole.influence!(influence, target_buffer, source_system::AbstractBody, source_buffer)
    for i in 1:size(target_buffer, 2)
        v = FastMultipole.get_gradient(target_buffer, i)
        n = FastMultipole.get_normal(source_buffer, source_system, i)
        influence[i] = dot(v, n)
    end
end

##### END OF ABSTRACT BODY #####################################################