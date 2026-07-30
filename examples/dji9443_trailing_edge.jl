import FLOWPanel as pnl
import GeoIO
using LinearAlgebra: dot, norm

"""
    find_dji9443_trailing_edge_indices(msh_file; watertight)

Find the two DJI 9443 trailing-edge seed triples in `msh_file`.  Each returned
triple is `[outermost, second_outermost, innermost]`, using Julia's 1-based node
indices.  The negative-`y` blade is returned first, followed by the positive-`y`
blade.

The detector assumes the repository's DJI 9443 orientation: `x` is axial, `y`
is radial, and the trailing side has matching `y` and `z` signs.
"""
function find_dji9443_trailing_edge_indices(msh_file; watertight)
    isfile(msh_file) || error("DJI 9443 mesh file does not exist: $(msh_file)")

    mesh = GeoIO.load(msh_file).geometry
    nodes, cells = pnl.meshes2nodes_cells(mesh)
    return _find_dji9443_trailing_edge_indices(nodes, cells; watertight,
                                                mesh_label=string(msh_file))
end

function _find_dji9443_trailing_edge_indices(nodes::AbstractMatrix,
                                              cells::AbstractMatrix{<:Integer};
                                              watertight::Bool,
                                              mesh_label::AbstractString="mesh")
    size(nodes, 1) == 3 ||
        error("DJI 9443 trailing-edge detection requires 3-D nodes; " *
              "$(mesh_label) has size $(size(nodes))")
    size(cells, 1) == 3 ||
        error("DJI 9443 trailing-edge detection requires triangular cells; " *
              "$(mesh_label) has $(size(cells, 1)) nodes per cell")

    actual_watertight, open_cells =
        pnl.iswatertight(nodes, cells; return_open_cells=true)
    if actual_watertight != watertight
        detail = actual_watertight ? "no boundary edges" :
                 "$(length(open_cells)) cells incident to boundary/nonmanifold edges"
        error("Topology mismatch for $(mesh_label): watertight=$(watertight) was " *
              "requested, but the mesh is $(actual_watertight ? "watertight" : "not watertight") " *
              "($(detail)).")
    end

    final_cells = pnl.ensure_consistent_winding(
        nodes, Matrix{Int}(cells); watertight=actual_watertight)
    normals = zeros(promote_type(eltype(nodes), Float64), 3, size(final_cells, 2))
    pnl.calc_normals!(nodes, final_cells, normals)
    edge_to_cells = pnl._calc_edge_to_cells(final_cells)

    radius = maximum(abs, view(nodes, 2, :))
    isfinite(radius) && radius > 0 ||
        error("Cannot determine the DJI 9443 radial scale from $(mesh_label)")

    # The physical trailing edge is a sharp, nearly radial two-sided ridge.
    # Dimensionless direction and radius tests keep detection invariant under
    # uniform scaling.
    normal_jump_tol = 0.2
    radial_alignment_tol = 0.7
    candidate_edges = Tuple{Int,Int}[]
    for (edge, adjacent) in edge_to_cells
        length(adjacent) == 2 || continue
        delta = view(nodes, :, edge[2]) .- view(nodes, :, edge[1])
        edge_length = norm(delta)
        edge_length > 0 || continue
        abs(delta[2]) / edge_length >= radial_alignment_tol || continue

        panel_1, panel_2 = adjacent[1][1], adjacent[2][1]
        normal_jump = 1 - dot(view(normals, :, panel_1),
                              view(normals, :, panel_2))
        normal_jump >= normal_jump_tol || continue
        push!(candidate_edges, edge)
    end

    isempty(candidate_edges) &&
        error("No sharp, predominantly radial two-sided edges were found in $(mesh_label)")

    node_neighbors = Dict{Int,Vector{Int}}()
    for (a, b) in candidate_edges
        push!(get!(node_neighbors, a, Int[]), b)
        push!(get!(node_neighbors, b, Int[]), a)
    end

    components = Vector{Vector{Int}}()
    visited = Set{Int}()
    for start in keys(node_neighbors)
        start in visited && continue
        component = Int[]
        stack = Int[start]
        push!(visited, start)
        while !isempty(stack)
            node = pop!(stack)
            push!(component, node)
            for neighbor in node_neighbors[node]
                if !(neighbor in visited)
                    push!(visited, neighbor)
                    push!(stack, neighbor)
                end
            end
        end
        push!(components, component)
    end

    blade_candidates = Dict(-1 => Vector{Vector{Int}}(),
                             1 => Vector{Vector{Int}}())
    for component in components
        endpoints = filter(node -> length(node_neighbors[node]) == 1, component)
        length(endpoints) == 2 || continue
        all(node -> length(node_neighbors[node]) <= 2, component) || continue

        outer = component[argmax(abs.(nodes[2, component]))]
        inner = component[argmin(abs.(nodes[2, component]))]
        outer_radius = abs(nodes[2, outer])
        inner_radius = abs(nodes[2, inner])
        outer_radius >= 0.9radius || continue
        inner_radius <= 0.2radius || continue
        outer in endpoints && inner in endpoints || continue

        # For this fixed rotor orientation, y*z > 0 is the trailing side;
        # y*z < 0 is the similarly sharp leading edge.
        nodes[2, outer] * nodes[3, outer] > 0 || continue
        blade_sign = nodes[2, outer] < 0 ? -1 : 1
        push!(blade_candidates[blade_sign], component)
    end

    triples = Vector{Vector{Int}}()
    for blade_sign in (-1, 1)
        candidates = blade_candidates[blade_sign]
        side = blade_sign < 0 ? "negative-y" : "positive-y"
        isempty(candidates) &&
            error("Missing $(side) DJI 9443 trailing-edge chain in $(mesh_label)")
        length(candidates) == 1 ||
            error("Ambiguous $(side) DJI 9443 trailing edge in $(mesh_label): " *
                  "found $(length(candidates)) complete sharp radial chains")

        component = only(candidates)
        outer = component[argmax(abs.(nodes[2, component]))]
        inner = component[argmin(abs.(nodes[2, component]))]
        neighbors = node_neighbors[outer]
        length(neighbors) == 1 ||
            error("Ambiguous outer tip at node $(outer) on the $(side) blade")
        second_outer = only(neighbors)

        trace = try
            pnl.trace_trailing_edge(
                nodes, final_cells, outer, second_outer;
                end_node=inner, normal_jump_tol, max_turn_angle=pi / 3)
        catch err
            error("Failed to validate the $(side) DJI 9443 trailing-edge chain " *
                  "in $(mesh_label): $(sprint(showerror, err))")
        end
        expected_edge_count = length(component) - 1
        if trace.nodes[end] != inner ||
           length(trace.edges) != expected_edge_count ||
           Set(trace.nodes) != Set(component)
            error("Incomplete $(side) DJI 9443 trailing-edge trace in $(mesh_label): " *
                  "expected $(expected_edge_count) edges ending at node $(inner), " *
                  "got $(length(trace.edges)) edges ending at node $(trace.nodes[end])")
        end

        push!(triples, [outer, second_outer, inner])
    end

    return triples[1], triples[2]
end
