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
mean(xs) = sum(xs)/length(xs)

"""
`simplewing(b, ar, tr, twist_root, twist_tip, lambda, gamma;
bodytype=RigidWakeBody,
span_NDIVS="automatic", rfl_NDIVS="automatic",
airfoil_root="naca6412.dat", airfoil_tip="naca6412.dat",
airfoil_path=def_rfl_path)`

Generates a symmetric single-section wing.

**ARGUMENTS**
  * `b::Real`         : Span.
  * `ar::Real`        : Aspect ratio defined as b/c_tip.
  * `tr::Real`        : Taper ratio defined as c_tip/c_root.
  * `twist_root::Real`: (deg) twist of the root.
  * `twist_tip::Real` : (deg) twist of the tip.
  * `lambda::Real`    : (deg) sweep.
  * `gamma::Real`     : (deg) dihedral.

**OPTIONAL ARGUMENTS**
  * `bodytype::Type{LBodyTypes}`: Type of lifting body to generate.
  * `span_NDIVS::ndivstype`     : Spanwise divisions.
  * `rfl_NDIVS::ndivstype`    : Chordwise divisions.
  * `airfoil_root::String`      : File to root airfoil contour.
  * `airfoil_tip::String`       : File to tip airfoil contour.
  * `airfoil_path::String`      : Path to airfoil files.

NOTE: See gt.multidscretize for a description of arguments of type `ndivstype`.
NOTE2: In the current implementation, sweep and dihedral are done about the LE.
"""
function simplewing(b::Number, ar::Number, tr::Number, twist_root::Number,
                      twist_tip::Number, lambda::Number, gamma::Number;
                      bodytype::Type{<:AbstractLiftingBody}=RigidWakeBody{Union{ConstantSource,ConstantDoublet}},
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

    # ----------------- GEOMETRY DESCRIPTION -----------------------------------
    c_tip = b/ar                        # Tip chord
    c_root = c_tip/tr                   # Root chord
    semispan = b/2                      # (m) semi-span length

    y_tip = b/2
    x_tip = y_tip*tan(lambda*pi/180)
    z_tip = y_tip*tan(gamma*pi/180)

    chords = [0.00 c_root/semispan;     # (semi-span position, chord c/semib)
              1.00 c_tip/semispan]

    twists = [0.0 twist_root;           # (semi-span position, twist (deg))
              1.0 twist_tip]

    x_pos = [0.00 0;                    # (semi-span position, LE x-position x/semib)
              1.00 x_tip/semispan]

    z_pos = [0.00 0;                    # (semi-span position, LE z-position x/semib)
              1.00 z_tip/semispan]

    airfoil_files = [(0.0, airfoil_root), # (semi-span position, airfoil file)
                     (1.0, airfoil_tip)]

    # ----------------- DISCRETIZATION 0000-------------------------------------

    # Defines divisions
    if span_NDIVS==nothing
        b_NDIVS = [(1.0, 35, 20.0, true)]  # Span cell sections
    else
        b_NDIVS = span_NDIVS
    end

    if rfl_NDIVS==nothing
        urfl_NDIVS = [(0.25, 7,   10.0, false),   # Cells on upper side of airfoils
                      (0.50,  5,    1.0, true),
                      (0.25,  6, 1/10.0, false)]
    else
        urfl_NDIVS = rfl_NDIVS
    end

    lrfl_NDIVS = urfl_NDIVS             # Cells on lower side of airfoils


    # ----------------- LOFTING PARAMETERS -------------------------------------
    # b_low = -1.0                        # Lower bound of span lofting
    # b_up = 1.0                          # Upper bound of span lofting
    symmetric = true                    # Lofting symmetric about b=0
    spl_k = 1                           # Spline order of distributions along span
    # spl_s = 0.0000001                 # Spline smoothing of distribution along span
    # rflspl_s = 0.00000001             # Spline smoothing of airfoil cross sections.
    # verify_spline = false             # Plots the splined distributions
    # verify_rflspline = true           # Plots the splined airfoil cross sections

    return generate_loft_liftbody(bodytype, airfoil_files, airfoil_path,
                                        urfl_NDIVS, lrfl_NDIVS,
                                        semispan, b_low, b_up, b_NDIVS,
                                        chords, twists, x_pos, z_pos;
                                        loop_dim=1, dimsplit=1,
                                        symmetric=symmetric,
                                        spl_k=spl_k, spl_s=spl_s,
                                        verify_spline=verify_spline,
                                        verify_rflspline=verify_rflspline,
                                        rflspl_s=rflspl_s,
                                        opt_args...
                                    )
end


"""
    `generate_multibody(bodytype::Type{<:AbstractLiftingBody},
                                meshfiles::AbstractVector{Tuple{String, String, Bool}},
                                trailingedges::AbstractVector{Tuple{String, Function, Function, Bool}},
                                reader::Function;
                                read_path=pwd(),
                                offset=zeros(3),
                                rotation::Meshes.Rotation=one(Meshes.QuatRotation),
                                scaling=1.0,
                                tolerance=1e-6,
                                lengthscale=nothing,
                                panel_omit_shedding=Dict(),
                                verbose=true,
                                v_lvl=0
                                )`

Reads a collection of mesh files and stitches them together into a MultiBody.

Returns `(multibody, elsprescribe)` where `multibody` is the MultiBody that
was generated and `elsprescribe` is the array of elements to prescribe in
least-squares solver (`solve(multibody, ...; elprescribe=elsprescribe, ...)`).

# Arguments
* `bodytype`:               Converts all meshes into this body type
* `meshfiles`:              Array with all mesh files to read, defined as
                                `meshfiles[i] = (name, meshfile, flip)` where
                                `name` is the name body name to give to this
                                mesh, `meshfile` is the file name, and
                                `flip=true` indicates that the mesh has normals
                                pointing inside and control points need to be
                                flipped to end up outside of the geometry.
* `trailingedges`:          Array with trailing edges that can be identified in
                                `meshfiles`, defined as `trailingedges[i] =
                                (trailingedge_meshfile, sorting_function,
                                junction_criterion, closed)` where
                                `trailingedge_meshfile` is a mesh file with the
                                line along the trailing edge,
                                `sorting_function = X::Vector -> val::Number` is
                                a function used for sorting the mesh points that
                                are identified to be along this trailing edge,
                                `junction_criterion = X::Vector ->
                                val::Number` is a function used for identifying
                                any wing-body junction along the trailing edge
                                that needs to be ignored (meaning, the point at
                                that junction is not flagged as a TE point), and
                                `closed=true` indicates that the trailing edge
                                is a closed contour (like in the case of a duct
                                , otherwise use `closed=false`)
* `reader`:                 A function for reading mesh files and return a
                                Meshes.jl object.
* `read_path`:              Path where to read files from.
* `offset`, `rotation`, `scaling`: Translates, rotates, and scales each mesh by
                                these inputs.
* `tolerance`:              Tolerance used to identify trailing edge points.
                                Any mesh points read from `meshfiles` that are
                                between a radius `tolerance*lengthscale` from
                                any points read from `trailingedges` will be
                                flagged as trailing edge nodes, so make sure
                                that the lines in `trailingedges` are finely
                                discretized and so you can use a tight (small
                                value) tolerance.

`sorting_function` in `trailingedges[i] = (trailingedge_meshfile,
sorting_function, junction_criterion, closed)` receives a point `X` and needs to
return a scalar value.
This value is used to sort TE points and to identify contiguous TE cells.
For a wing with its span direction aligned along the y-axis, we recommend
defining this function as
```julia
import LinearAlgebra: dot

Xcg = [0, 0, 0]                 # Aircraft center of gravity for reference
spandirection = [0, 1, 0]       # Span direction can be any arbitrary direction

sorting_function(X) = dot(X - Xcg, spandirection)
```
For the trailing edge of a circular duct (or any closed cross section) with
its centerline aligned along the x-axis, we recommend defining this function as
```julia
Xc = [0, 0, 0]                  # Section centroid

sorting_function(X) = -atand(X[2] - Xc[2], X[3] - Xc[3])
```
`junction_criterion` also receives a point `X` and needs to return a scalar
value. A negative or zero value means that this point is at a wing-body junction
and should not be considered a trailing edge point.
For example, for an aircraft where the wing meets the fuselage at y=±0.3, this
function should be defined as follows
```julia
y_junction = 0.3

junction_criterion(X) = abs(X[2]) - (y_junction + 1e-6) # Offset y by a small
                                                        # value to account for
                                                        # discretization error
```
If there are no junctions, simply define this as `junction_criterion(X) = Inf`

`reader` must be a user-defined function to read the mesh files and return a
Meshes.jl object. We recommend the following function:
```julia
import GeoIO

function reader(file)
    msh = GeoIO.load(file)
    return msh.geometry
end
```
"""
function generate_multibody(bodytype::Type{<:AbstractLiftingBody},
                            meshfiles::AbstractVector{Tuple{String, String, Bool}},
                            trailingedges::AbstractVector{Tuple{String, F1, F2, Bool}},
                            reader::Function;
                            read_path=pwd(),
                            offset=zeros(3),
                            rotation::Meshes.Rotation=one(Meshes.QuatRotation),
                            scaling=1.0,
                            tolerance=1e-6,
                            lengthscale=nothing,
                            rediscretize_trailingedge=true,
                            panel_omit_shedding=nothing,
                            debug=false,
                            verbose=true,
                            v_lvl=0
                            ) where {F1<:Union{Function, Any}, F2<:Union{Function, Any}}

    @assert 0 < tolerance <= 1 ""*
        "Received invalid tolerance $(tolerance); it needs to be between 0 and 1."

    @assert !isnothing(lengthscale) && lengthscale > 0 ""*
        "Received invalid length scale $(lengthscale)"

    bodies = bodytype[]
    names = String[]
    elsprescribe = Tuple{Int64, Float64}[]

    # panel_omit_shedding = Dict()
    # Dictionary flagging particle to omit shedding at the junctions between wing
    # and fuselage. This has the form
    # `panel_omit_shedding[body][ei] = (shed incoming, shed outgoing, shed unsteady)`

    # Generate each body
    for (name, meshfile, flip) in meshfiles

        if verbose; println("\t"^(v_lvl)*"Reading $(name) ($(meshfile))"); end

        # Read Gmsh mesh
        msh = reader(joinpath(read_path, meshfile))

        # Transform the original mesh: Translate, rotate, and scale
        msh = msh |> Meshes.Translate(offset...) |> Meshes.Rotate(rotation) |> Meshes.Scale(scaling)

        # Wrap Meshes object into a Grid object from GeometricTools
        grid = gt.GridTriangleSurface(msh)

        this_panel_omit_shedding = Dict()

        # All shedding points (trailing edge) are agglomerated here
        sheddings = [noshedding]

        nremoved = 0

        # Read all trailing edges
        for (trailingedgefile, sortingfunction, junctioncriterion, closed) in trailingedges

            # Read Gmsh line of trailing edge
            TEmsh = reader(joinpath(read_path, trailingedgefile))

            # Apply the same transformations of the mesh to the trailing edge
            TEmsh = TEmsh |> Meshes.Translate(offset...) |> Meshes.Rotate(rotation) |> Meshes.Scale(scaling)

            # Convert TE Meshes object into a matrix of points used to identify the trailing edge
            trailingedge = gt.vertices2nodes(TEmsh.vertices)

            # Sort TE points from "left" to "right" according to span direction
            trailingedge = sortslices(trailingedge; dims=2, by=sortingfunction)

            # Rediscretize the TE line with enough resolution to meet tolerance
            if rediscretize_trailingedge
                ndiscretize = ceil(Int, 1/(0.5*tolerance))
                trailingedge = gt.rediscretize_line(trailingedge, ndiscretize;
                                                    parameterization=sortingfunction)
            end

            # Filter out any points that are close to junctions
            toremove = filter( i -> junctioncriterion(view(trailingedge, :, i)) <= 0.0, 1:size(trailingedge, 2) )
            removedtrailingedge = trailingedge[:, toremove]

            tokeep = filter( i -> junctioncriterion(view(trailingedge, :, i)) > 0.0, 1:size(trailingedge, 2) )
            trailingedge = trailingedge[:, tokeep]

            # Generate TE shedding matrix
            shedding = calc_shedding(grid, trailingedge;
                                            tolerance=lengthscale*tolerance,
                                            periodic=closed)

            # Generate TE shedding matrix of nodes that got filtered out
            removedshedding = calc_shedding(grid, removedtrailingedge;
                                            tolerance=lengthscale*tolerance,
                                            periodic=false)

            nremoved += size(removedshedding, 2)

            # Flag nodes closest to the wing-body junction to not shed particles
            if !isnothing(panel_omit_shedding) && size(shedding, 2)!=0

                nodes = grid._nodes
                prevnshedding = sum(size.(sheddings, 2); init=0)

                tricoor, quadcoor, lin, ndivscells, cin = gt.generate_getcellt_args(grid)

                # Find shedding nodes that are next to the junctions
                for (rpi, rnia, rnib) in eachcol(removedshedding) # Iterate over junction-removed cells

                    # Convert node indices of removed cell from panel-local to global
                    rpia = gt.get_cell_t(tricoor, quadcoor, grid, rpi, rnia, lin, ndivscells, cin)
                    rpib = gt.get_cell_t(tricoor, quadcoor, grid, rpi, rnib, lin, ndivscells, cin)

                    if debug
                        println("\t"^(v_lvl+2), "rpi=$rpi\trpia=$rpia\trpib=$rpib")
                    end

                    for (ei, (pi, nia, nib)) in enumerate(eachcol(shedding)) # Iterate over shedding cells

                        # Convert node indices from panel-local to global
                        pia = gt.get_cell_t(tricoor, quadcoor, grid, pi, nia, lin, ndivscells, cin)
                        pib = gt.get_cell_t(tricoor, quadcoor, grid, pi, nib, lin, ndivscells, cin)

                        # Check if incoming node is also a node in the removed cell
                        omit_a = ( pia==rpia || pia==rpib
                                    # || norm(nodes[:, pia] - nodes[:, rpia]) <= 1e2*eps()
                                    # || norm(nodes[:, pia] - nodes[:, rpib]) <= 1e2*eps()
                                    )

                        # Check if outgoing node is also a node in the removed cell
                        omit_b = ( pib==rpia || pib==rpib
                                    # || norm(nodes[:, pib] - nodes[:, rpia]) <= 1e2*eps()
                                    # || norm(nodes[:, pib] - nodes[:, rpib]) <= 1e2*eps()
                                    )

                        # Case that the cell shared a node with a junction-removed cell
                        if omit_a || omit_b

                            # Flags for FLOWUnsteady to omit shedding or not
                            omitflags = (omit_a, omit_b, false)

                            if debug
                                println("\t"^(v_lvl+3), "ei=$ei", "\t", "pi=$pi", "\t", omitflags, "\t", "pia=$pia", "\t", "pib=$pib")
                            end

                            # Save the cell and node to omit
                            this_panel_omit_shedding[ei+prevnshedding] = omitflags

                        end

                    end

                end

            end

            push!(sheddings, shedding)

        end

        # Combine all TE shedding matrices into one
        shedding = hcat(sheddings...)

        # Generate paneled body
        body = bodytype(grid, shedding; CPoffset=(-1)^flip * 1e-14)

        # Check if watertight
        watertight = gt.isclosed(msh)

        # Function for calculating elements to prescribe in least-square solver
        # (none if open mesh, arbitrary if watertight)
        elprescribe = watertight ? calc_elprescribe(body) : Tuple{Int, Float64}[]

        # Offset the element number by the running total of panels in the multibody
        for (eli, val) in elprescribe
            neweli = eli + sum(getproperty.(bodies, :ncells); init=0)
            push!(elsprescribe, (neweli, val))
        end

        # Store TE nodes to avoid shedding
        if !isnothing(panel_omit_shedding)
            panel_omit_shedding[body] = this_panel_omit_shedding
        end

        # Store all the bodies of the multibody
        push!(bodies, body)
        push!(names, name)

        if verbose; println("\t"^(v_lvl+1)*"Is mesh wateright?\t$(watertight)"); end
        if verbose; println("\t"^(v_lvl+1)*"Number of panels:\t$(body.ncells)"); end
        if verbose; println("\t"^(v_lvl+1)*"Number of sheddings:\t$(body.nsheddings) ($(nremoved) removed)"); end
        if verbose; println("\t"^(v_lvl+1)*"Number of shed omits:\t$(length(this_panel_omit_shedding))"); end

    end

    # Build multibody
    multibody = MultiBody(bodies, names)

    if verbose; println("\t"^(v_lvl)*"Total number of panels:\t$(multibody.ncells)"); end

    return multibody, elsprescribe

end

direction(dir; X0=zeros(3)) = X -> dot(X - X0, dir)
loop(; Oaxis=Matrix(1.0I, 3, 3), X0=zeros(3)) = X -> -atand(dot(X - X0, Oaxis[:, 2]), dot(X - X0, Oaxis[:, 3]))
nojunction(X) = Inf


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

> **NOTE:** `body` cannot be a MultiBody.
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
    field = get_field(body, fieldname)["field_data"]

    # Find index of row or column slicing the field
    gdim = row ? 1 : 2                          # Dimension to slice
    islice, pos, errmin, lin = find_i(body, controlpoints, position, gdim, -1;
                                                xdir=direction, filter=filter)

    # Slice field
    ncell = get_ndivscells(body)[row ? 2 : 1]   # Number of cells in the slice
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

    # Calculate winding number associated to control points
    windings = calcfield_winding(body.grid.orggrid, controlpoints)

    # Save field in body
    if addfield
        add_field(body, fieldname, "scalar", windings, "cell")
    end

    return windings
end

function calcfield_winding(mesh::gt.Meshes.Mesh, controlpoints::AbstractMatrix)

    points = ( gt.Meshes.Point(gt.Meshes.Vec(p...)) for p in eachcol(controlpoints) )
    windings = gt.Meshes.winding(points, mesh)

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
