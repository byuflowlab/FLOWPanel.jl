"""
    open_boundary_loops(cells) -> Vector{Vector{Int}}

Chain the open (singly-referenced) edges of a triangle mesh into node loops.
Each returned loop follows the surface-cell traversal direction.
"""
function open_boundary_loops(cells::AbstractMatrix{<:Integer})
    incidence = Dict{Tuple{Int,Int},Int}()
    for ci in axes(cells, 2)
        n1, n2, n3 = cells[1, ci], cells[2, ci], cells[3, ci]
        for (a, b) in ((n1, n2), (n2, n3), (n3, n1))
            key = a < b ? (a, b) : (b, a)
            incidence[key] = get(incidence, key, 0) + 1
        end
    end

    successor = Dict{Int,Int}()
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

Close every open boundary loop with a flat centroid-fan cap. New nodes and
cells are appended, preserving all existing indices.
"""
function add_flat_tip_caps(nodes::AbstractMatrix,
                           cells::AbstractMatrix{<:Integer})
    loops = open_boundary_loops(cells)
    isempty(loops) && return copy(nodes), Matrix{Int}(cells), Int[], Int[]

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
            out_cells[:, ci] .= (b, a, ni)
            push!(cap_cells, ci)
        end
    end
    return out_nodes, out_cells, cap_nodes, cap_cells
end
