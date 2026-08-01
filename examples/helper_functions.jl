import FLOWPanel as pnl
import GeometricTools as gt

function gt.Grid(P_min, P_max, NDIVS, loop_dim::Integer)
    return gt.Grid(P_min, P_max, NDIVS; loop_dim=loop_dim)
end

function grid2cells(grid)
    cells = zeros(Int, 3, grid.ncells)
    for i in 1:grid.ncells
        cells[:, i] .= gt.get_cell(grid, i)
    end
    return cells
end

function (bodytype::Type{B})(grid::gt.GridTriangleSurface; kwargs...) where {B<:pnl.NonLiftingBody}
    nodes, cells = pnl.meshes2nodes_cells(grid.orggrid)
    watertight, _ = pnl.iswatertight(nodes, cells)
    return bodytype(nodes, cells; watertight, kwargs...)
end

function (bodytype::Type{B})(grid::gt.GridTriangleSurface, shedding; kwargs...) where {B<:pnl.AbstractLiftingBody}
    nodes, cells = pnl.meshes2nodes_cells(grid.orggrid)
    watertight, _ = pnl.iswatertight(nodes, cells)
    return bodytype(nodes, cells, shedding; watertight, kwargs...)
end

function _surface_revolution_compat(profile::AbstractMatrix{T}, thetaNDIVS::TN;
                                    loop_dim::Integer=0,
                                    axis_angle::Number=0,
                                    low_a::Number=0,
                                    up_a::Number=360,
                                    save_path=nothing,
                                    paraview::Bool=true,
                                    file_name::AbstractString="myrev",
                                    ) where {T<:Real, TN}
    try
        return gt.surface_revolution(profile, thetaNDIVS;
                                     loop_dim=loop_dim,
                                     axis_angle=axis_angle,
                                     low_a=low_a,
                                     up_a=up_a,
                                     save_path=save_path,
                                     paraview=paraview,
                                     file_name=file_name)
    catch e
        if !(e isa MethodError && e.f === gt.Grid)
            rethrow()
        end
    end

    if size(profile, 2) != 2
        error("Invalid point dimensions in `profile`. Expected 2 dimensions, got $(size(profile,2))")
    elseif profile[1, :] == profile[end, :] && loop_dim == 0
        @warn("Received a closed contour but parametric grid wasn't declared to loop.")
    end

    NDIVS = if TN <: Number
        [size(profile, 1) - 1, thetaNDIVS, 0]
    else
        [[(1.0, size(profile, 1) - 1, 1.0, false)],
         thetaNDIVS,
         [(1.0, 0, 1.0, false)]]
    end

    P_min = [0, low_a, 0]
    P_max = [1, up_a, 0]
    grid = gt.Grid(P_min, P_max, NDIVS; loop_dim=loop_dim)

    if axis_angle != 0
        M = gt.rotation_matrix(0, 0, -axis_angle)
        M2D = M[2:3, 2:3]
        M3D = collect(M)'
        points = collect(hcat([M2D * profile[i, :] for i in 1:size(profile, 1)]...))'
    else
        M3D = I
        points = profile
    end

    function my_space_transform(X, ind)
        p_ind = ind[1]
        angle = X[2]
        point = [0, points[p_ind, 1], points[p_ind, 2]]
        return M3D * gt.axis_rotation(Float64[0, 0, 1], Float64(angle)) * point
    end

    gt.transform3!(grid, my_space_transform)

    if save_path != nothing
        gt.save(grid, file_name; path=save_path)
        if paraview
            run(`paraview --data=$(joinpath(save_path, file_name)).vtk`)
        end
    end

    return grid
end

function generate_loft(bodytype::Type{B}, args...; bodyoptargs=(),
                       dimsplit::Int64=2, optargs...) where {B<:pnl.NonLiftingBody}
    grid = gt.generate_loft(args...; optargs...)
    triang_grid = gt.GridTriangleSurface(grid, dimsplit)
    nodes = triang_grid._nodes
    cells = grid2cells(triang_grid)
    watertight, _ = pnl.iswatertight(nodes, cells)
    return bodytype(nodes, cells; watertight, bodyoptargs...)
end

function generate_revolution(bodytype::Type{B}, args...; bodyoptargs=(),
                             dimsplit::Int64=2, loop_dim::Int64=2,
                             gridprocessing=nothing,
                             optargs...) where {B<:pnl.NonLiftingBody}
    grid = _surface_revolution_compat(args...; loop_dim=loop_dim, optargs...)
    if gridprocessing !== nothing
        gridprocessing(grid)
    end
    triang_grid = gt.GridTriangleSurface(grid, dimsplit)
    nodes = triang_grid._nodes
    cells = grid2cells(triang_grid)
    watertight, _ = pnl.iswatertight(nodes, cells)
    return bodytype(nodes, cells; watertight, bodyoptargs...)
end

function generate_loft_liftbody(bodytype::Type{B}, args...;
                                bodyoptargs=(), dimsplit::Int=2,
                                overwrite_shedding=nothing,
                                optargs...) where {B<:pnl.AbstractBody}
    grid = gt.generate_loft(args...; optargs...)
    triang_grid = gt.GridTriangleSurface(grid, dimsplit)
    nodes = triang_grid._nodes
    cells = grid2cells(triang_grid)

    ndivs = gt.get_ndivscells(triang_grid)
    U = [Base._sub2ind(ndivs, ndivs[1] - 1, i) for i in 1:ndivs[2]]
    L = [Base._sub2ind(ndivs, 2, i) for i in 1:ndivs[2]]

    shedding = if isnothing(overwrite_shedding)
        out = zeros(Int, 6, length(U))
        for (ei, (u, l)) in enumerate(zip(U, L))
            out[:, ei] .= (u, 3, 2, l, 3, 2)
        end
        out
    else
        overwrite_shedding
    end

    watertight, _ = pnl.iswatertight(nodes, cells)
    if bodytype <: pnl.AbstractLiftingBody
        return bodytype(nodes, cells, shedding; watertight, bodyoptargs...)
    elseif bodytype <: pnl.NonLiftingBody
        return bodytype(nodes, cells; watertight, bodyoptargs...)
    else
        error("Body type $(bodytype) is not a lifting or non-lifting body.")
    end
end

function generate_revolution_liftbody(bodytype::Type{B}, args...;
                                      bodyoptargs=(),
                                      gridprocessing=nothing,
                                      dimsplit::Int=1,
                                      loop_dim::Int=2,
                                      axis_angle=270,
                                      overwrite_shedding=nothing,
                                      closed_contour=true,
                                      optargs...) where {B<:pnl.AbstractBody}
    grid = _surface_revolution_compat(args...; loop_dim=loop_dim,
                                      axis_angle=axis_angle, optargs...)

    if gridprocessing == nothing
        Oaxis = gt.rotation_matrix2(0, 0, 90)
        O = zeros(3)
        gt.lintransform!(grid, Oaxis, O)
    else
        gridprocessing(grid)
    end

    triang_grid = gt.GridTriangleSurface(grid, dimsplit)
    nodes = triang_grid._nodes
    cells = grid2cells(triang_grid)

    shedding = overwrite_shedding
    if isnothing(shedding)
        ndivs = gt.get_ndivscells(triang_grid)
        U = [Base._sub2ind(ndivs, ndivs[1] - 1, i) for i in 1:ndivs[2]]
        L = [Base._sub2ind(ndivs, 2, i) for i in 1:ndivs[2]]

        shedding = zeros(Int, 6, length(U))
        for (ei, (u, l)) in enumerate(zip(U, L))
            shedding[:, ei] .= (u, 3, 2, closed_contour ? l : -1, 3, 2)
        end
    end

    if bodytype <: pnl.AbstractLiftingBody
        return bodytype(nodes, cells, [shedding]; bodyoptargs...)
    elseif bodytype <: pnl.NonLiftingBody
        watertight, _ = pnl.iswatertight(nodes, cells)
        return bodytype(nodes, cells; watertight, bodyoptargs...)
    else
        error("Body type $(bodytype) is not a lifting or non-lifting body.")
    end
end

"""
    open_boundary_loops(cells) -> Vector{Vector{Int}}

Chain the open (singly-referenced) edges of a triangle mesh into node loops.
Each returned loop lists its nodes in the direction the *surface* cells traverse
them, so a cap face closing edge `(a, b)` must traverse `b -> a` to stay wound
consistently with the surface. Returns an empty vector for a watertight mesh.

This is the connectivity half of `pnl.iswatertight`, kept here because the
library returns open *cells* rather than the ordered boundary chains a cap
needs.
"""
function open_boundary_loops(cells::AbstractMatrix{<:Integer})
    incidence = Dict{Tuple{Int, Int}, Int}()
    for ci in axes(cells, 2)
        n1, n2, n3 = cells[1, ci], cells[2, ci], cells[3, ci]
        for (a, b) in ((n1, n2), (n2, n3), (n3, n1))
            key = a < b ? (a, b) : (b, a)
            incidence[key] = get(incidence, key, 0) + 1
        end
    end

    successor = Dict{Int, Int}()
    for ci in axes(cells, 2)
        n1, n2, n3 = cells[1, ci], cells[2, ci], cells[3, ci]
        for (a, b) in ((n1, n2), (n2, n3), (n3, n1))
            key = a < b ? (a, b) : (b, a)
            incidence[key] == 1 || continue
            haskey(successor, a) && error("open boundary is not a simple loop: " *
                "node $a starts two boundary edges (non-manifold boundary, or " *
                "inconsistent surface winding)")
            successor[a] = b
        end
    end

    loops = Vector{Vector{Int}}()
    visited = Set{Int}()
    for start in sort!(collect(keys(successor)))
        start in visited && continue
        loop = Int[]
        node = start
        while !(node in visited)
            push!(visited, node)
            push!(loop, node)
            haskey(successor, node) || error("open boundary chain dead-ends at " *
                "node $node; the boundary is not a closed loop")
            node = successor[node]
        end
        node == start || error("open boundary loops are not disjoint (chain from " *
            "node $start re-entered an earlier loop at node $node)")
        push!(loops, loop)
    end
    return loops
end

"""
    add_flat_tip_caps(nodes, cells) -> (nodes, cells, cap_nodes, cap_cells)

Close every open boundary loop of a triangle mesh with a flat centroid-fan cap:
one new node at the mean of each loop's nodes, fanned to every loop edge. New
nodes and cells are *appended*, so existing node and cell indices — and hence
any shedding matrix already expressed in them — are unchanged.

Centroid-fan is chosen because it inherits the contour's symmetry: if the loop
nodes are invariant (as a set) under a reflection, so is the fan, since the
centroid maps to the centroid and each fan face to the fan face of the image
edge. A single quad/ear triangulation would not.

`cap_nodes` and `cap_cells` are the appended indices, for symmetry and
watertightness checks. Returns copies; the inputs are not mutated.

Motivation: the interior-Dirichlet Green identity assumes a closed surface, so
an open-tip lofted wing is outside the formulation's assumptions.
"""
function add_flat_tip_caps(nodes::AbstractMatrix, cells::AbstractMatrix{<:Integer})
    loops = open_boundary_loops(cells)
    if isempty(loops)
        return copy(nodes), Matrix{Int}(cells), Int[], Int[]
    end

    n_node0 = size(nodes, 2)
    n_cell0 = size(cells, 2)
    out_nodes = Matrix{eltype(nodes)}(undef, 3, n_node0 + length(loops))
    out_cells = Matrix{Int}(undef, 3, n_cell0 + sum(length, loops))
    out_nodes[:, 1:n_node0] .= nodes
    out_cells[:, 1:n_cell0] .= cells

    cap_nodes = Int[]
    cap_cells = Int[]
    ni, ci = n_node0, n_cell0
    for loop in loops
        ni += 1
        for d in 1:3
            out_nodes[d, ni] = sum(nodes[d, n] for n in loop) / length(loop)
        end
        push!(cap_nodes, ni)
        for k in eachindex(loop)
            a = loop[k]
            b = loop[k == length(loop) ? firstindex(loop) : k + 1]
            ci += 1
            # Reverse the surface's traversal of (a, b) so the cap face is wound
            # consistently with the cell it closes against.
            out_cells[:, ci] .= (b, a, ni)
            push!(cap_cells, ci)
        end
    end
    return out_nodes, out_cells, cap_nodes, cap_cells
end

function simplewing(b::Number, ar::Number, tr::Number, twist_root::Number,
                    twist_tip::Number, lambda::Number, gamma::Number;
                    bodytype::Type{<:pnl.AbstractBody}=pnl.RigidWakeBody{Union{pnl.ConstantSource,pnl.ConstantDoublet}},
                    span_NDIVS=nothing,
                    rfl_NDIVS=nothing,
                    airfoil_root::String="naca6412.dat",
                    airfoil_tip::String="naca6412.dat",
                    airfoil_path::String=pnl.def_rfl_path,
                    spl_s::Real=0.0000001,
                    rflspl_s::Real=0.00000001,
                    verify_spline::Bool=true,
                    verify_rflspline::Bool=true,
                    b_low=-1.0, b_up=1.0,
                    opt_args...)
    c_tip = b/ar
    c_root = c_tip/tr
    semispan = b/2

    y_tip = b/2
    x_tip = y_tip*tan(lambda*pi/180)
    z_tip = y_tip*tan(gamma*pi/180)

    chords = [0.00 c_root/semispan; 1.00 c_tip/semispan]
    twists = [0.0 twist_root; 1.0 twist_tip]
    x_pos = [0.00 0; 1.00 x_tip/semispan]
    z_pos = [0.00 0; 1.00 z_tip/semispan]
    airfoil_files = [(0.0, airfoil_root), (1.0, airfoil_tip)]

    b_NDIVS = span_NDIVS == nothing ? [(1.0, 35, 20.0, true)] : span_NDIVS
    urfl_NDIVS = rfl_NDIVS == nothing ?
                 [(0.25, 7, 10.0, false), (0.50, 5, 1.0, true), (0.25, 6, 0.1, false)] :
                 rfl_NDIVS
    lrfl_NDIVS = urfl_NDIVS

    return generate_loft_liftbody(bodytype, airfoil_files, airfoil_path,
                                  urfl_NDIVS, lrfl_NDIVS,
                                  semispan, b_low, b_up, b_NDIVS,
                                  chords, twists, x_pos, z_pos;
                                  loop_dim=1, dimsplit=1,
                                  symmetric=true,
                                  spl_k=1, spl_s=spl_s,
                                  verify_spline=verify_spline,
                                  verify_rflspline=verify_rflspline,
                                  rflspl_s=rflspl_s,
                                  opt_args...)
end

function simplewing_mirrored(b::Number, ar::Number, tr::Number, twist_root::Number,
                             twist_tip::Number, lambda::Number, gamma::Number;
                             bodytype::Type{<:pnl.AbstractBody}=pnl.RigidWakeBody{Union{pnl.ConstantSource,pnl.ConstantDoublet}},
                             span_NDIVS=nothing,
                             rfl_NDIVS=nothing,
                             airfoil_root::String="naca6412.dat",
                             airfoil_tip::String="naca6412.dat",
                             airfoil_path::String=pnl.def_rfl_path,
                             spl_s::Real=0.0000001,
                             rflspl_s::Real=0.00000001,
                             verify_spline::Bool=true,
                             verify_rflspline::Bool=true,
                             mirror_tol::Real=100eps(Float64),
                             bodyoptargs=(;),
                             caps::Symbol=:none,
                             opt_args...)
    c_tip = b/ar
    c_root = c_tip/tr
    semispan = b/2

    y_tip = b/2
    x_tip = y_tip*tan(lambda*pi/180)
    z_tip = y_tip*tan(gamma*pi/180)

    chords = [0.00 c_root/semispan; 1.00 c_tip/semispan]
    twists = [0.0 twist_root; 1.0 twist_tip]
    x_pos = [0.00 0; 1.00 x_tip/semispan]
    z_pos = [0.00 0; 1.00 z_tip/semispan]
    airfoil_files = [(0.0, airfoil_root), (1.0, airfoil_tip)]

    b_NDIVS = span_NDIVS == nothing ? [(1.0, 35, 20.0, true)] : span_NDIVS
    urfl_NDIVS = rfl_NDIVS == nothing ?
                 [(0.25, 7, 10.0, false), (0.50, 5, 1.0, true), (0.25, 6, 0.1, false)] :
                 rfl_NDIVS
    lrfl_NDIVS = urfl_NDIVS

    grid = gt.generate_loft(airfoil_files, airfoil_path,
                            urfl_NDIVS, lrfl_NDIVS,
                            semispan, 0.0, 1.0, b_NDIVS,
                            chords, twists, x_pos, z_pos;
                            loop_dim=1,
                            symmetric=false,
                            spl_k=1, spl_s=spl_s,
                            verify_spline=verify_spline,
                            verify_rflspline=verify_rflspline,
                            rflspl_s=rflspl_s,
                            opt_args...)
    triang_grid = gt.GridTriangleSurface(grid, 1)
    half_nodes = triang_grid._nodes
    half_cells = grid2cells(triang_grid)

    ndivs = gt.get_ndivscells(triang_grid)
    U = [Base._sub2ind(ndivs, ndivs[1] - 1, i) for i in 1:ndivs[2]]
    L = [Base._sub2ind(ndivs, 2, i) for i in 1:ndivs[2]]
    half_shedding = zeros(Int, 6, length(U))
    for (ei, (u, l)) in enumerate(zip(U, L))
        half_shedding[:, ei] .= (u, 3, 2, l, 3, 2)
    end

    mirror_index = Vector{Int}(undef, size(half_nodes, 2))
    nodes = copy(half_nodes)
    for ni in axes(half_nodes, 2)
        if abs(half_nodes[2, ni]) <= mirror_tol
            mirror_index[ni] = ni
        else
            nodes = hcat(nodes, [half_nodes[1, ni], -half_nodes[2, ni], half_nodes[3, ni]])
            mirror_index[ni] = size(nodes, 2)
        end
    end

    half_centers_y = [sum(half_nodes[2, half_cells[:, ci]]) / 3 for ci in axes(half_cells, 2)]
    pos_order = sort(collect(axes(half_cells, 2)); by=ci -> half_centers_y[ci])
    neg_order = sort(collect(axes(half_cells, 2)); by=ci -> -half_centers_y[ci])

    cells = Matrix{Int}(undef, 3, 2 * size(half_cells, 2))
    out_ci = 0
    for ci in neg_order
        out_ci += 1
        cells[:, out_ci] .= reverse(mirror_index[half_cells[:, ci]])
    end
    for ci in pos_order
        out_ci += 1
        cells[:, out_ci] .= half_cells[:, ci]
    end

    if caps === :flat
        # Close the two open tips before anything reads the connectivity: the
        # interior-Dirichlet Green identity assumes a closed surface. Caps are
        # appended, so `half_cells`-derived TE node indices stay valid.
        nodes, cells, _, _ = add_flat_tip_caps(nodes, cells)
    elseif caps !== :none
        throw(ArgumentError("caps must be :none or :flat; got $caps"))
    end

    watertight, _ = pnl.iswatertight(nodes, cells)
    if bodytype <: pnl.AbstractLiftingBody
        te_nodes = Int[]
        for col in eachcol(half_shedding)
            pi, nia, nib = col[1], col[2], col[3]
            push!(te_nodes, half_cells[nia, pi])
            push!(te_nodes, half_cells[nib, pi])
        end
        full_te_nodes = unique(vcat(mirror_index[te_nodes], te_nodes))
        sort!(full_te_nodes; by=ni -> nodes[2, ni])
        shedding = pnl.calc_shedding(nodes, cells, full_te_nodes, zeros(eltype(nodes), 3, 0))
        final_bodyoptargs = merge((ensure_winding=false,), bodyoptargs)
        return bodytype(nodes, cells, [shedding]; watertight, final_bodyoptargs...)
    elseif bodytype <: pnl.NonLiftingBody
        final_bodyoptargs = merge((ensure_winding=false,), bodyoptargs)
        return bodytype(nodes, cells; watertight, final_bodyoptargs...)
    else
        error("Body type $(bodytype) is not a lifting or non-lifting body.")
    end
end
