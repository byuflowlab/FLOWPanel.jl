#=##############################################################################
# DESCRIPTION
    Definition of methods for post-processing solver results.

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez and Ryan Anderson
  * Email       : Edo.AlvarezR@gmail.com; rymanderson@gmail.com
  * Date        : Oct 2022
  * License     : MIT License
=###############################################################################


################################################################################
# VELOCITY FIELDS
################################################################################

const _GRAD_MU_OPTION_KEYS = Set{Symbol}((
    :basis,
    :quad_normal_dot_min,
    :quad_grow,
    :quad_grow_max_depth,
    :quad_grow_stop,
    :quad_grow_cond_max,
    :tri_robust,
    :tri_robust_ar_threshold,
    :tri_robust_max_depth,
    :tri_robust_target_healthy,
))
const _NORMALIZED_GRAD_MU_OPTION_KEYS = Set{Symbol}((
    :basis,
    :quad_normal_dot_min,
    :quad_grow,
    :quad_grow_max_depth,
    :quad_grow_stop,
    :quad_grow_cond_max,
    :tri_robust,
    :tri_robust_ar_threshold,
    :tri_robust_max_depth,
    :tri_robust_target_healthy,
))

_grad_mu_get(opts::NamedTuple, key::Symbol, default) =
    haskey(opts, key) ? getproperty(opts, key) : default

function _grad_mu_options_namedtuple(grad_mu_options)
    grad_mu_options isa NamedTuple ||
        throw(ArgumentError("grad_mu_options must be a NamedTuple, got $(typeof(grad_mu_options))."))
    for key in keys(grad_mu_options)
        key in _GRAD_MU_OPTION_KEYS ||
            throw(ArgumentError("Unknown grad_mu_options key :$(key)."))
    end
    return grad_mu_options
end

function _validate_grad_mu_basis(basis)
    basis in (:tri, :quad) ||
        throw(ArgumentError("Unknown grad_mu basis $(basis); use :tri or :quad."))
    return basis
end

function _normalize_grad_mu_options(grad_mu_options=(;);
        default_basis::Symbol=:tri,
        default_mode::Union{Nothing,Symbol}=nothing)
    opts = _grad_mu_options_namedtuple(grad_mu_options)

    if default_mode !== nothing
        default_basis = default_mode === :tri ? :tri : :quad
    end

    basis_opt = _grad_mu_get(opts, :basis, default_basis)
    basis_opt = _validate_grad_mu_basis(basis_opt)

    normal_dot_min = Float64(_grad_mu_get(opts, :quad_normal_dot_min, cos(pi / 4)))

    tri_robust_max_depth = Int(_grad_mu_get(opts, :tri_robust_max_depth, 4))
    tri_robust_target_healthy = Int(_grad_mu_get(opts, :tri_robust_target_healthy, 6))
    quad_grow_max_depth = Int(_grad_mu_get(opts, :quad_grow_max_depth, 2))
    stop = _grad_mu_get(opts, :quad_grow_stop, :cond)
    stop in (:cond, :depth) ||
        throw(ArgumentError("Unknown quad_grow_stop=$(stop); use :cond or :depth."))
    quad_grow_cond_max = Float64(_grad_mu_get(opts, :quad_grow_cond_max, 1e3))
    quad_grow = Bool(_grad_mu_get(opts, :quad_grow, true))
    tri_robust = Bool(_grad_mu_get(opts, :tri_robust, false))
    tri_robust_ar_threshold = Float64(_grad_mu_get(opts, :tri_robust_ar_threshold, 10.0))

    normalized_input = Set(Symbol.(keys(opts))) == _NORMALIZED_GRAD_MU_OPTION_KEYS
    if basis_opt === :quad && !normalized_input
        for key in (:tri_robust, :tri_robust_ar_threshold, :tri_robust_max_depth, :tri_robust_target_healthy)
            haskey(opts, key) && throw(ArgumentError("grad_mu_options.$(key) is only valid with basis=:tri."))
        end
    end

    tri_robust_max_depth >= 0 || throw(ArgumentError("tri_robust_max_depth must be nonnegative, got $(tri_robust_max_depth)."))
    tri_robust_target_healthy >= 0 || throw(ArgumentError("tri_robust_target_healthy must be nonnegative, got $(tri_robust_target_healthy)."))
    quad_grow_max_depth >= 0 || throw(ArgumentError("quad_grow_max_depth must be nonnegative, got $(quad_grow_max_depth)."))
    quad_grow_cond_max > 0 || throw(ArgumentError("quad_grow_cond_max must be positive, got $(quad_grow_cond_max)."))
    tri_robust_ar_threshold > 0 || throw(ArgumentError("tri_robust_ar_threshold must be positive, got $(tri_robust_ar_threshold)."))

    return (;
        basis=basis_opt,
        quad_normal_dot_min=normal_dot_min,
        quad_grow=quad_grow,
        quad_grow_max_depth=quad_grow_max_depth,
        quad_grow_stop=stop,
        quad_grow_cond_max=quad_grow_cond_max,
        tri_robust=tri_robust,
        tri_robust_ar_threshold=tri_robust_ar_threshold,
        tri_robust_max_depth=tri_robust_max_depth,
        tri_robust_target_healthy=tri_robust_target_healthy,
    )
end

"""
    calcfield_U!(bodies; backend=DirectBackend(), reset=true, convolve_panels=true, doublet_gradient=true)
    calcfield_U!(body; optargs...)

Compute and store the control-point velocity field for one or more bodies.
"""
function calcfield_U!(bodies::Tuple, uinf=SVector{3,Float64}(0.0, 0.0, 0.0), wakes=Tuple(nothing for _ in bodies);
        backend::AbstractBackend=DirectBackend(),
        reset=true,
        convolve_panels=true,
        doublet_gradient=true,
        grad_mu_options=(;),
    )
    normalized_grad_mu_options = _normalize_grad_mu_options(grad_mu_options;
        default_basis=:quad)

    # reset velocity
    if reset
        for body in bodies
            body.velocity .= zero(eltype(body.velocity))
        end
    end

    # recalculate normals/control points on the target body
    for body in bodies
        calc_normals!(body)
        calc_controlpoints!(body)
    end

    # add wake-induced velocity at control points
    wakes = Tuple(w for w in wakes if !isnothing(w))
    length(wakes) > 0 && influence!(bodies, wakes, backend; scalar_potential=false, velocity=true)

    # Add induced velocity at each control point
    convolve_panels && influence!(bodies, bodies, backend; scalar_potential=false, velocity=true)

    # add doublet gradient (if applicable)
    for body in bodies
        if has_grad_mu(body) && doublet_gradient
            compute_mu_gradient!(body.velocity, body.controlpoints, body.normals,
                body.cells,
                body.neighbor,
                view(body.strength, :, get_Gammai(body)),
                view(body.shedding_full, 1:2, :); scale=0.5,
                nodes=body.nodes,
                grad_mu_options=normalized_grad_mu_options)
            
            # alternatively, comment out the above function and:
            # body.velocity .*= 2.0 # (only works if self-induced velocity is the only one applied)
            # influence!(body, body, backend; scalar_potential=false, velocity=true)
        

            # alternatively, comment out the above and:
            # body.velocity .= zero(eltype(body.velocity))
            # compute_mu_gradient!(body.velocity, body.controlpoints, body.normals,
            #     body.cells,
            #     body.neighbor,
            #     view(body.strength, :, get_Gammai(body)),
            #     view(body.shedding_full, 1:2, :), 
            #     scale=1.0)
        end
    end

    # apply freestream at control points
    for body in bodies
        apply_freestream!(body, uinf)
    end

    return nothing
end

calcfield_U!(targetbody, args...; optargs...) = calcfield_U!((targetbody,), args...; optargs...)

################################################################################
# GRADIENT COMPUTATION
################################################################################

"""
    panel_aspect_ratio_mask(nodes, cells; threshold=10.0)

Flag panels whose longest-edge / shortest-edge ratio exceeds `threshold`.
Such panels are degenerate enough that the 1-ring LS stencil in
`compute_mu_gradient!` cannot reliably recover the in-plane gradient. Set
`grad_mu_options=(; tri_robust=true, tri_robust_ar_threshold=...)` to trigger
the BFS-gathered healthy-only fallback.
"""
function panel_aspect_ratio_mask(nodes::AbstractMatrix{<:Real},
                                 cells::AbstractMatrix{Int};
                                 threshold::Real=10.0)
    _, M = size(cells)
    mask = falses(M)
    @inbounds for i in 1:M
        ns1, ns2, ns3 = cells[1, i], cells[2, i], cells[3, i]
        e1 = sqrt((nodes[1,ns2]-nodes[1,ns1])^2 + (nodes[2,ns2]-nodes[2,ns1])^2 + (nodes[3,ns2]-nodes[3,ns1])^2)
        e2 = sqrt((nodes[1,ns3]-nodes[1,ns2])^2 + (nodes[2,ns3]-nodes[2,ns2])^2 + (nodes[3,ns3]-nodes[3,ns2])^2)
        e3 = sqrt((nodes[1,ns1]-nodes[1,ns3])^2 + (nodes[2,ns1]-nodes[2,ns3])^2 + (nodes[3,ns1]-nodes[3,ns3])^2)
        emin = min(e1, e2, e3); emax = max(e1, e2, e3)
        mask[i] = emin > 0 && emax / emin > threshold
    end
    return mask
end

"""
    compute_mu_gradient!(grad_mu, controlpoints, normals, cells, neighbors, mu, te_info;
                         scale=0.5,
                         nodes=nothing, grad_mu_options=(;))

Compute the surface gradient of doublet or vortex-ring strength using a local
least-squares stencil, with one-sided handling at trailing-edge panels.

When `grad_mu_options.tri_robust=true`, panels with aspect ratio above
`tri_robust_ar_threshold` skip the 1-ring stencil and instead inherit a gradient
computed from a BFS-gathered set of healthy panels reached by walking the
`neighbors` graph and *stepping through* flagged panels. This makes the
reconstruction robust to clusters of degenerate panels (e.g. cap-corner slivers
on capped wing meshes). TE exclusion mirrors the pass-1 behavior: when the bad
panel is itself a TE panel, the cross-TE neighbor is rejected as the first BFS
step. The bad-panel fallback uses `tri_robust_target_healthy` and
`tri_robust_max_depth`.

Use `grad_mu_options=(; basis=...)` to select the reconstruction. `basis=:tri`
uses the per-triangle least-squares stencil. `basis=:quad` reconstructs directly
from paired-quad/singleton agglomerate-averaged strength differences and may grow
its agglomerate stencil when `quad_grow=true`, using `quad_grow_stop`,
`quad_grow_cond_max`, and `quad_grow_max_depth`.
"""
function compute_mu_gradient!(grad_mu,
                            controlpoints::AbstractMatrix{Float64},
                            normals::AbstractMatrix{Float64},
                            cells::AbstractMatrix{Int},
                            neighbors::AbstractMatrix{Int},
                            mu::AbstractVector{Float64},
                            te_info::AbstractMatrix{Int};
                            scale=0.5,
                            nodes::Union{Nothing,AbstractMatrix{Float64}}=nothing,
                            grad_mu_options=(;),
                        )
    opts = _normalize_grad_mu_options(grad_mu_options;
        default_basis=isnothing(nodes) ? :tri : :quad)
    bad_panel_mask = nothing
    if opts.tri_robust
        isnothing(nodes) && throw(ArgumentError("grad_mu_options.tri_robust=true requires nodes for aspect-ratio masking."))
        mask = panel_aspect_ratio_mask(nodes, cells; threshold=opts.tri_robust_ar_threshold)
        any(mask) && (bad_panel_mask = mask)
    end
    return _compute_mu_gradient_masked!(grad_mu, controlpoints, normals, cells,
        neighbors, mu, te_info; scale, bad_panel_mask, nodes,
        grad_mu_options=opts)
end

function _bfs_traverse!(on_discover, neighbors_of, starts;
        max_depth::Int,
        visited::AbstractVector{Bool},
        frontier::Vector{Int},
        next_frontier::Vector{Int},
        initial_depth::Int=0,
        previsited=(),
        stop=() -> false,
        edge_allowed=(p, q, depth) -> true)
    fill!(visited, false)
    empty!(frontier)
    empty!(next_frontier)
    @inbounds for p in previsited
        p > 0 && (visited[p] = true)
    end
    @inbounds for p in starts
        p > 0 || continue
        visited[p] = true
        push!(frontier, p)
    end

    depth = initial_depth
    while !isempty(frontier) && depth < max_depth
        stop() && break
        next_depth = depth + 1
        empty!(next_frontier)
        for p in frontier
            for q in neighbors_of(p)
                q > 0 || continue
                visited[q] && continue
                edge_allowed(p, q, next_depth) || continue
                visited[q] = true
                on_discover(q, next_depth)
                push!(next_frontier, q)
            end
        end
        isempty(next_frontier) && break
        frontier, next_frontier = next_frontier, frontier
        depth = next_depth
    end
    return nothing
end

function _compute_mu_gradient_masked!(grad_mu,
                            controlpoints::AbstractMatrix{Float64},
                            normals::AbstractMatrix{Float64},
                            cells::AbstractMatrix{Int},
                            neighbors::AbstractMatrix{Int},
                            mu::AbstractVector{Float64},
                            te_info::AbstractMatrix{Int};
                            scale=0.5,
                            bad_panel_mask::Union{Nothing,AbstractVector{Bool}}=nothing,
                            nodes::Union{Nothing,AbstractMatrix{Float64}}=nothing,
                            grad_mu_options=(;),
                        )
    opts = _normalize_grad_mu_options(grad_mu_options;
        default_basis=isnothing(nodes) ? :tri : :quad)
    _, M = size(cells)
    local_grad = zeros(Float64, 3, M)

    if opts.basis === :quad
        isnothing(nodes) && throw(ArgumentError("grad_mu_options.basis=:quad requires nodes."))
        bad_panel_mask === nothing ||
            throw(ArgumentError("bad_panel_mask is only valid with grad_mu_options.basis=:tri."))
        _quad_mu_diff_gradient!(local_grad, nodes, cells, neighbors, normals,
            te_info, mu; scale=scale, normal_dot_min=opts.quad_normal_dot_min,
            grow=opts.quad_grow, grow_stop=opts.quad_grow_stop,
            cond_max=opts.quad_grow_cond_max, grow_max_depth=opts.quad_grow_max_depth)
        grad_mu .+= local_grad
        return grad_mu
    end

    # Pre-allocate array for stencil to avoid allocations inside the loop
    stencil = Int[]
    sizehint!(stencil, 10)

    is_bad = bad_panel_mask

    for i in 1:M
        if is_bad !== nothing && is_bad[i]
            continue  # handled in the second (BFS-gather) pass below
        end
        empty!(stencil)

        # Check if current panel is at the Trailing Edge
        is_te = te_info[1, i] > 0 && te_info[2, i] > 0

        if is_te
            # Global vertex IDs of the TE edge
            te_v1 = cells[te_info[1, i], i]
            te_v2 = cells[te_info[2, i], i]
        else
            te_v1 = -1
            te_v2 = -1
        end

        # 1. Gather immediate valid neighbors
        for k in 1:3
            n_idx = neighbors[k, i]
            if n_idx <= 0
                continue
            end

            if is_te
                # Check if this neighbor shares BOTH trailing edge vertices
                has_v1 = (cells[1, n_idx] == te_v1) || (cells[2, n_idx] == te_v1) || (cells[3, n_idx] == te_v1)
                has_v2 = (cells[1, n_idx] == te_v2) || (cells[2, n_idx] == te_v2) || (cells[3, n_idx] == te_v2)

                # If it shares the TE edge, exclude it from the stencil (one-sided bias)
                if has_v1 && has_v2
                    continue
                end
            end
            push!(stencil, n_idx)
        end

        # 2. Expand stencil upstream if it's a TE panel to ensure robust Least Squares
        if is_te
            n_current = length(stencil)
            for s in 1:n_current
                s_idx = stencil[s]
                
                # Look at neighbors of the interior neighbors
                for k in 1:3
                    nn_idx = neighbors[k, s_idx]
                    if nn_idx > 0 && nn_idx != i && !(nn_idx in stencil)
                        push!(stencil, nn_idx)
                    end
                end
            end
        end

        # 3. Build normal equations for Least Squares: (A^T A) * grad = A^T b
        ATA = zeros(Float64, 3, 3)
        ATb = zeros(Float64, 3)
        mean_sq_dist = 0.0

        for j in stencil
            dx = controlpoints[1, j] - controlpoints[1, i]
            dy = controlpoints[2, j] - controlpoints[2, i]
            dz = controlpoints[3, j] - controlpoints[3, i]
            dmu = mu[j] - mu[i]

            ATA[1, 1] += dx * dx; ATA[1, 2] += dx * dy; ATA[1, 3] += dx * dz
            ATA[2, 1] += dy * dx; ATA[2, 2] += dy * dy; ATA[2, 3] += dy * dz
            ATA[3, 1] += dz * dx; ATA[3, 2] += dz * dy; ATA[3, 3] += dz * dz

            ATb[1] += dx * dmu
            ATb[2] += dy * dmu
            ATb[3] += dz * dmu

            mean_sq_dist += dx^2 + dy^2 + dz^2
        end

        mean_sq_dist = length(stencil) > 0 ? (mean_sq_dist / length(stencil)) : 1.0

        # 4. Constrain gradient to the panel surface (Penalty method)
        # This acts as an orthogonal constraint mapping it cleanly to the 3D plane
        nx = normals[1, i]
        ny = normals[2, i]
        nz = normals[3, i]
        
        # Scale penalty by geometry size (mean squared dist) for matrix conditioning
        penalty = 1e4 * mean_sq_dist
        ATA[1, 1] += penalty * nx * nx; ATA[1, 2] += penalty * nx * ny; ATA[1, 3] += penalty * nx * nz
        ATA[2, 1] += penalty * ny * nx; ATA[2, 2] += penalty * ny * ny; ATA[2, 3] += penalty * ny * nz
        ATA[3, 1] += penalty * nz * nx; ATA[3, 2] += penalty * nz * ny; ATA[3, 3] += penalty * nz * nz

        # 5. Add minor Tikhonov regularization in case of an exactly coplanar/rank-deficient system
        reg = 1e-10 * mean_sq_dist
        ATA[1, 1] += reg
        ATA[2, 2] += reg
        ATA[3, 3] += reg

        # 6. Solve the 3x3 local system
        g = ATA \ ATb

        local_grad[1, i] -= g[1] * scale
        local_grad[2, i] -= g[2] * scale
        local_grad[3, i] -= g[3] * scale
    end

    # ---------- Second pass: BFS-gathered healthy-only stencil for bad panels ----------
    if is_bad !== nothing
        visited = falses(M)
        frontier = Int[]
        next_frontier = Int[]
        healthy = Int[]
        sizehint!(frontier, 32)
        sizehint!(next_frontier, 32)
        sizehint!(healthy, 32)

        for i in 1:M
            is_bad[i] || continue

            empty!(healthy)

            is_te_i = te_info[1, i] > 0 && te_info[2, i] > 0
            te_v1 = is_te_i ? cells[te_info[1, i], i] : -1
            te_v2 = is_te_i ? cells[te_info[2, i], i] : -1

            neighbors_of = p -> view(neighbors, :, p)
            edge_allowed = function (p, q, depth)
                # TE exclusion mirrors the pass-1 one-sided rule: only when
                # the bad panel itself is a TE panel and q sits across its TE edge.
                if p == i && is_te_i
                    has_v1 = (cells[1, q] == te_v1) || (cells[2, q] == te_v1) || (cells[3, q] == te_v1)
                    has_v2 = (cells[1, q] == te_v2) || (cells[2, q] == te_v2) || (cells[3, q] == te_v2)
                    return !(has_v1 && has_v2)
                end
                return true
            end
            on_discover = function (q, depth)
                is_bad[q] || push!(healthy, q)
                return nothing
            end
            _bfs_traverse!(on_discover, neighbors_of, (i,);
                max_depth=opts.tri_robust_max_depth,
                visited, frontier, next_frontier,
                stop=() -> length(healthy) >= opts.tri_robust_target_healthy,
                edge_allowed)

            if length(healthy) < 3
                # Best-effort: leave grad_mu unchanged on this panel. Caller
                # already zero-initialized via calcfield_U!'s reset path. This
                # only triggers for fully-enclosed degenerate regions, which
                # indicate an unusable mesh.
                continue
            end

            ATA = zeros(Float64, 3, 3)
            ATb = zeros(Float64, 3)
            mean_sq_dist = 0.0
            for j in healthy
                dx = controlpoints[1, j] - controlpoints[1, i]
                dy = controlpoints[2, j] - controlpoints[2, i]
                dz = controlpoints[3, j] - controlpoints[3, i]
                dmu = mu[j] - mu[i]
                ATA[1,1] += dx*dx; ATA[1,2] += dx*dy; ATA[1,3] += dx*dz
                ATA[2,1] += dy*dx; ATA[2,2] += dy*dy; ATA[2,3] += dy*dz
                ATA[3,1] += dz*dx; ATA[3,2] += dz*dy; ATA[3,3] += dz*dz
                ATb[1] += dx*dmu; ATb[2] += dy*dmu; ATb[3] += dz*dmu
                mean_sq_dist += dx*dx + dy*dy + dz*dz
            end
            mean_sq_dist /= length(healthy)

            nx = normals[1, i]; ny = normals[2, i]; nz = normals[3, i]
            penalty = 1e4 * mean_sq_dist
            ATA[1,1] += penalty*nx*nx; ATA[1,2] += penalty*nx*ny; ATA[1,3] += penalty*nx*nz
            ATA[2,1] += penalty*ny*nx; ATA[2,2] += penalty*ny*ny; ATA[2,3] += penalty*ny*nz
            ATA[3,1] += penalty*nz*nx; ATA[3,2] += penalty*nz*ny; ATA[3,3] += penalty*nz*nz

            reg = 1e-10 * mean_sq_dist
            ATA[1,1] += reg; ATA[2,2] += reg; ATA[3,3] += reg

            g = ATA \ ATb

            local_grad[1, i] -= g[1] * scale
            local_grad[2, i] -= g[2] * scale
            local_grad[3, i] -= g[3] * scale
        end
    end

    grad_mu .+= local_grad

    return grad_mu
end

function _edge_vertices(cells::AbstractMatrix{Int}, edge_i::Int, cell_i::Int)
    edge_i == 1 && return cells[1, cell_i], cells[2, cell_i]
    edge_i == 2 && return cells[2, cell_i], cells[3, cell_i]
    return cells[3, cell_i], cells[1, cell_i]
end

function _shares_te_edge(cells::AbstractMatrix{Int},
                         te_info::AbstractMatrix{Int},
                         i::Int,
                         j::Int)
    te_info[1, i] > 0 && te_info[2, i] > 0 || return false
    te_v1 = cells[te_info[1, i], i]
    te_v2 = cells[te_info[2, i], i]
    has_v1 = (cells[1, j] == te_v1) || (cells[2, j] == te_v1) || (cells[3, j] == te_v1)
    has_v2 = (cells[1, j] == te_v2) || (cells[2, j] == te_v2) || (cells[3, j] == te_v2)
    return has_v1 && has_v2
end

# Normalized signed offset (~sin of half-angle) of node p from the line through
# apex_i with direction (dx,dy,dz), measured in the tangent plane with normal n̂.
@inline function _signed_offset(nodes::AbstractMatrix{Float64}, apex_i::Int, p::Int,
                                dx::Float64, dy::Float64, dz::Float64, dn::Float64,
                                nx::Float64, ny::Float64, nz::Float64)
    px = nodes[1,p]-nodes[1,apex_i]; py = nodes[2,p]-nodes[2,apex_i]; pz = nodes[3,p]-nodes[3,apex_i]
    pn = sqrt(px*px + py*py + pz*pz)
    pn > 0 || return 0.0
    cx = dy*pz - dz*py; cy = dz*px - dx*pz; cz = dx*py - dy*px
    return (cx*nx + cy*ny + cz*nz) / (dn*pn)
end

# Mean |cos| of the four corner angles of the quad with ordered corners p1,p2,p3,p4.
# 0 for any rectangle of any aspect ratio (right angles); ~0.71 for a 45°/135°
# parallelogram (a domino). Aspect-invariant, so it keeps thin TE rectangles.
@inline function _quad_rightangle_cost(nodes::AbstractMatrix{Float64},
                                       p1::Int, p2::Int, p3::Int, p4::Int)
    e1x = nodes[1,p2]-nodes[1,p1]; e1y = nodes[2,p2]-nodes[2,p1]; e1z = nodes[3,p2]-nodes[3,p1]
    e2x = nodes[1,p3]-nodes[1,p2]; e2y = nodes[2,p3]-nodes[2,p2]; e2z = nodes[3,p3]-nodes[3,p2]
    e3x = nodes[1,p4]-nodes[1,p3]; e3y = nodes[2,p4]-nodes[2,p3]; e3z = nodes[3,p4]-nodes[3,p3]
    e4x = nodes[1,p1]-nodes[1,p4]; e4y = nodes[2,p1]-nodes[2,p4]; e4z = nodes[3,p1]-nodes[3,p4]
    l1 = sqrt(e1x^2+e1y^2+e1z^2); l2 = sqrt(e2x^2+e2y^2+e2z^2)
    l3 = sqrt(e3x^2+e3y^2+e3z^2); l4 = sqrt(e4x^2+e4y^2+e4z^2)
    cdot(ax,ay,az,la,bx,by,bz,lb) = (la>0 && lb>0) ? abs(ax*bx+ay*by+az*bz)/(la*lb) : 1.0
    return (cdot(e1x,e1y,e1z,l1,e2x,e2y,e2z,l2) + cdot(e2x,e2y,e2z,l2,e3x,e3y,e3z,l3) +
            cdot(e3x,e3y,e3z,l3,e4x,e4y,e4z,l4) + cdot(e4x,e4y,e4z,l4,e1x,e1y,e1z,l1)) / 4
end

"""
    _quad_best_candidate(nodes, cells, neighbors, normals, te_info, i, taken; normal_dot_min)
        -> (best_j, best_cost, runnerup_cost, ncand)

For unpaired triangle `i`, scan its non-barrier, straddling, unpaired edge-neighbors
and return the lowest right-angle-cost candidate, the runner-up cost, and the number
of candidates. Helper for `_quad_pairing_propagate`.
"""
@inline function _quad_best_candidate(nodes, cells, neighbors, normals, te_info,
                                      i::Int, taken::AbstractVector{Bool};
                                      normal_dot_min::Float64)
    ci1, ci2, ci3 = cells[1, i], cells[2, i], cells[3, i]
    nx, ny, nz = normals[1, i], normals[2, i], normals[3, i]
    bj = 0; bc = Inf; rc = Inf; n = 0
    @inbounds for edge_i in 1:3
        j = neighbors[edge_i, i]
        (j > 0 && !taken[j]) || continue
        _shares_te_edge(cells, te_info, i, j) && continue
        nx*normals[1,j] + ny*normals[2,j] + nz*normals[3,j] >= normal_dot_min || continue
        u, w = _edge_vertices(cells, edge_i, i)
        apex_i = (ci1 != u && ci1 != w) ? ci1 : ((ci2 != u && ci2 != w) ? ci2 : ci3)
        aj = cells[1, j]; (aj == u || aj == w) && (aj = cells[2, j]); (aj == u || aj == w) && (aj = cells[3, j])
        dx = nodes[1,aj]-nodes[1,apex_i]; dy = nodes[2,aj]-nodes[2,apex_i]; dz = nodes[3,aj]-nodes[3,apex_i]
        dn = sqrt(dx*dx + dy*dy + dz*dz); dn > 0 || continue
        su = _signed_offset(nodes, apex_i, u, dx, dy, dz, dn, nx, ny, nz)
        sw = _signed_offset(nodes, apex_i, w, dx, dy, dz, dn, nx, ny, nz)
        su*sw < 0 || continue
        c = _quad_rightangle_cost(nodes, apex_i, u, aj, w)
        n += 1
        if c < bc; rc = bc; bc = c; bj = j; elseif c < rc; rc = c; end
    end
    return (bj, bc, rc, n)
end

"""
    _quad_pairing_propagate(nodes, cells, neighbors, normals, te_info; normal_dot_min, margin, max_rounds) -> partner

Pairing by forced-move (unit) propagation with a confidence/adjacency fallback,
designed to defeat the stable-domino failure of plain mutual-best matching. On a
swept mesh the logical TE cell and an adjacent-cell "domino" are geometrically
congruent and mutually prefer each other, so mutual-best locks onto the domino. But
boundary/TE cells typically have just ONE admissible partner (others are boundary or
trailing-edge edges), so they are *forced*; claiming a forced pair removes the domino
option for the neighbor and cascades inward. Each outer round: (1) run forced moves
(exactly one candidate) to a fixpoint, then (2) accept mutual-best pairs that are
either *confident* (cost margin over the runner-up `> margin`) or *adjacent to an
already-paired cell* (propagation). Genuinely ambiguous / unstructured triangles are
left unpaired for the gradient's grow step (`quad_grow`).
"""
function _quad_pairing_propagate(nodes::AbstractMatrix{Float64},
                                 cells::AbstractMatrix{Int},
                                 neighbors::AbstractMatrix{Int},
                                 normals::AbstractMatrix{Float64},
                                 te_info::AbstractMatrix{Int};
                                 normal_dot_min::Float64=cos(pi / 4),
                                 margin::Float64=0.1,
                                 max_rounds::Int=10000)
    M = size(cells, 2)
    partner = zeros(Int, M)
    taken = falses(M)
    bc(i) = _quad_best_candidate(nodes, cells, neighbors, normals, te_info, i, taken;
                                 normal_dot_min=normal_dot_min)
    has_taken_neighbor(i) = (@inbounds (neighbors[1,i] > 0 && taken[neighbors[1,i]]) ||
                                       (neighbors[2,i] > 0 && taken[neighbors[2,i]]) ||
                                       (neighbors[3,i] > 0 && taken[neighbors[3,i]]))
    pair!(i, j) = (partner[i] = j; partner[j] = i; taken[i] = true; taken[j] = true)

    for _ in 1:max_rounds
        # (1) forced-move cascade to a fixpoint
        while true
            forced = false
            @inbounds for i in 1:M
                taken[i] && continue
                (j, _, _, n) = bc(i)
                if n == 1 && !taken[j]
                    pair!(i, j); forced = true
                end
            end
            forced || break
        end
        # (2) one round of confident / propagation mutual-best pairs
        changed = false
        @inbounds for i in 1:M
            taken[i] && continue
            (j, bcost, rcost, n) = bc(i)
            (j > 0 && !taken[j]) || continue
            (jb, _, _, _) = bc(j)
            jb == i || continue                                  # mutual best
            confident = (n == 1) || (rcost - bcost > margin)
            (confident || has_taken_neighbor(i) || has_taken_neighbor(j)) || continue
            pair!(i, j); changed = true
        end
        changed || break
    end
    return partner
end

# 2x2 symmetric eigen-condition for the LS normal matrix [a b; b c] (MᵀM of the
# tangential design matrix). Returns the 2-norm condition number σmax/σmin
# (= sqrt(λmax/λmin)); Inf when rank-deficient.
@inline function _sym2x2_cond(a::Float64, b::Float64, c::Float64)
    tr = a + c
    d = sqrt(max(((a - c) / 2)^2 + b*b, 0.0))
    lmax = tr/2 + d
    lmin = tr/2 - d
    (lmin <= 0 || lmax <= 0) && return Inf
    return sqrt(lmax / lmin)
end

# Tangential LS gradient of `mu_agg` at agglomerate `ag` from the centroid
# offsets of the agglomerates in `stencil`. Returns (g3::NTuple{3,Float64},
# cond::Float64); cond=Inf signals an unusable (rank-deficient / <2-neighbor)
# stencil and the caller should grow the stencil or report a reconstruction
# failure.
function _agg_ls_gradient(ag::Int, stencil::AbstractVector{Int},
                          aggC::AbstractMatrix{Float64},
                          aggN::AbstractMatrix{Float64},
                          aggMu::AbstractVector{Float64})
    m = length(stencil)
    m >= 2 || return ((0.0, 0.0, 0.0), Inf)

    nx, ny, nz = aggN[1, ag], aggN[2, ag], aggN[3, ag]
    # tangent basis from the first neighbor offset projected onto the tangent plane
    d1x = aggC[1, stencil[1]] - aggC[1, ag]
    d1y = aggC[2, stencil[1]] - aggC[2, ag]
    d1z = aggC[3, stencil[1]] - aggC[3, ag]
    dn = d1x*nx + d1y*ny + d1z*nz
    t1x = d1x - dn*nx; t1y = d1y - dn*ny; t1z = d1z - dn*nz
    nrm = sqrt(t1x^2 + t1y^2 + t1z^2)
    nrm > 0 || return ((0.0, 0.0, 0.0), Inf)
    t1x /= nrm; t1y /= nrm; t1z /= nrm
    t2x = ny*t1z - nz*t1y; t2y = nz*t1x - nx*t1z; t2z = nx*t1y - ny*t1x

    a11 = 0.0; a12 = 0.0; a22 = 0.0; b1 = 0.0; b2 = 0.0
    @inbounds for nb in stencil
        dx = aggC[1, nb] - aggC[1, ag]
        dy = aggC[2, nb] - aggC[2, ag]
        dz = aggC[3, nb] - aggC[3, ag]
        p = dx*t1x + dy*t1y + dz*t1z          # offset in t1
        q = dx*t2x + dy*t2y + dz*t2z          # offset in t2
        dmu = aggMu[nb] - aggMu[ag]
        a11 += p*p; a12 += p*q; a22 += q*q
        b1 += p*dmu; b2 += q*dmu
    end

    cond = _sym2x2_cond(a11, a12, a22)
    det = a11*a22 - a12*a12
    det != 0 || return ((0.0, 0.0, 0.0), Inf)
    g1 = (a22*b1 - a12*b2) / det              # ∂μ/∂t1
    g2 = (a11*b2 - a12*b1) / det              # ∂μ/∂t2
    gx = g1*t1x + g2*t2x
    gy = g1*t1y + g2*t2y
    gz = g1*t1z + g2*t2z
    return ((gx, gy, gz), cond)
end

"""
    _quad_mu_diff_gradient!(local_grad, nodes, cells, neighbors, normals, te_info, mu; ...)

The proven quad-pitch reconstruction: agglomerate paired triangles into logical
quads (`_quad_pairing_propagate`) plus singleton triangles, average μ over each
agglomerate, and difference that agglomerate-averaged μ over the
agglomerate-centroid stencil with a tangential least-squares fit. Writes
`local_grad` (which holds `-scale·∇μ`) on the member triangles of every
agglomerate. Agglomerates with too few/ill-posed neighbors throw an
`ArgumentError` after permitted growth instead of falling back to a triangle
gradient — except *isolated* agglomerates (empty stencil: every edge neighbor
rejected by the TE barrier or fold test, as for degenerate sliver panels),
which receive a zero gradient with a warning, since no growth depth can ever
produce a stencil for them.

When `grow=true`, agglomerates are grown by BFS on the agglomerate graph.
`grow_stop=:cond` grows only ill-posed agglomerates (LS condition number
`> cond_max`, or singletons) until the condition number drops to `cond_max`,
capped at `grow_max_depth` rings; `grow_stop=:depth` grows every agglomerate a
fixed `grow_max_depth` rings. Edges across the trailing edge or folds sharper
than `normal_dot_min` are never crossed.
"""
function _quad_mu_diff_gradient!(local_grad::AbstractMatrix{Float64},
                                 nodes::AbstractMatrix{Float64},
                                 cells::AbstractMatrix{Int},
                                 neighbors::AbstractMatrix{Int},
                                 normals::AbstractMatrix{Float64},
                                 te_info::AbstractMatrix{Int},
                                 mu::AbstractVector{Float64};
                                 scale::Float64,
                                 normal_dot_min::Float64,
                                 grow::Bool,
                                 grow_stop::Symbol,
                                 cond_max::Float64,
                                 grow_max_depth::Int)
    M = size(cells, 2)
    partner = _quad_pairing_propagate(nodes, cells, neighbors, normals, te_info;
                                      normal_dot_min=normal_dot_min)

    # Canonical agglomerate id per triangle: lower member of a pair, else self.
    canon = Vector{Int}(undef, M)
    @inbounds for i in 1:M
        canon[i] = partner[i] == 0 ? i : min(i, partner[i])
    end

    # Per-agglomerate area-weighted geometry / strength, indexed by canonical id
    # (non-canonical entries stay zero and are skipped via aggA > 0).
    aggA = zeros(Float64, M)
    aggC = zeros(Float64, 3, M)
    aggN = zeros(Float64, 3, M)
    aggMu = zeros(Float64, M)
    @inbounds for i in 1:M
        a, b, c = cells[1, i], cells[2, i], cells[3, i]
        e1x = nodes[1,b]-nodes[1,a]; e1y = nodes[2,b]-nodes[2,a]; e1z = nodes[3,b]-nodes[3,a]
        e2x = nodes[1,c]-nodes[1,a]; e2y = nodes[2,c]-nodes[2,a]; e2z = nodes[3,c]-nodes[3,a]
        cx = e1y*e2z - e1z*e2y; cy = e1z*e2x - e1x*e2z; cz = e1x*e2y - e1y*e2x
        w = 0.5*sqrt(cx*cx + cy*cy + cz*cz)
        tcx = (nodes[1,a]+nodes[1,b]+nodes[1,c])/3
        tcy = (nodes[2,a]+nodes[2,b]+nodes[2,c])/3
        tcz = (nodes[3,a]+nodes[3,b]+nodes[3,c])/3
        ag = canon[i]
        aggA[ag] += w
        aggC[1,ag] += w*tcx; aggC[2,ag] += w*tcy; aggC[3,ag] += w*tcz
        aggN[1,ag] += w*normals[1,i]; aggN[2,ag] += w*normals[2,i]; aggN[3,ag] += w*normals[3,i]
        aggMu[ag] += w*mu[i]
    end
    @inbounds for ag in 1:M
        aggA[ag] > 0 || continue
        aggC[1,ag] /= aggA[ag]; aggC[2,ag] /= aggA[ag]; aggC[3,ag] /= aggA[ag]
        aggMu[ag] /= aggA[ag]
        nn = sqrt(aggN[1,ag]^2 + aggN[2,ag]^2 + aggN[3,ag]^2)
        if nn > 0
            aggN[1,ag] /= nn; aggN[2,ag] /= nn; aggN[3,ag] /= nn
        end
    end

    # Agglomerate adjacency: a,b adjacent iff a member triangle of a shares a
    # non-barrier edge with a member of b (barrier = TE edge or fold < dot_min).
    aggnbr = [Int[] for _ in 1:M]
    @inbounds for i in 1:M
        ag = canon[i]
        for k in 1:3
            q = neighbors[k, i]
            q > 0 || continue
            agq = canon[q]
            agq == ag && continue
            _shares_te_edge(cells, te_info, i, q) && continue
            ndot = normals[1,i]*normals[1,q] + normals[2,i]*normals[2,q] + normals[3,i]*normals[3,q]
            ndot >= normal_dot_min || continue
            agq in aggnbr[ag] || push!(aggnbr[ag], agq)
        end
    end

    # BFS frontier buffers reused across agglomerates (grow step)
    visited = falses(M)
    frontier = Int[]; next_frontier = Int[]; stencil = Int[]

    @inbounds for ag in 1:M
        aggA[ag] > 0 || continue              # canonical agglomerates only

        empty!(stencil)
        append!(stencil, aggnbr[ag])

        (g3, cnd) = _agg_ls_gradient(ag, stencil, aggC, aggN, aggMu)

        if grow
            need_grow = grow_stop === :depth ? true : (cnd > cond_max)
            if need_grow
                on_discover = function (agq, depth)
                    push!(stencil, agq)
                    (g3, cnd) = _agg_ls_gradient(ag, stencil, aggC, aggN, aggMu)
                    return nothing
                end
                _bfs_traverse!(on_discover, p -> aggnbr[p], stencil;
                    max_depth=grow_max_depth,
                    visited, frontier, next_frontier,
                    initial_depth=1,
                    previsited=(ag,),
                    stop=() -> grow_stop === :cond && cnd <= cond_max)
            end
        end

        gx, gy, gz = g3
        if !(isfinite(cnd) && isfinite(gx) && isfinite(gy) && isfinite(gz))
            m2 = partner[ag]
            members = m2 == 0 ? string(ag) : string(ag, ",", m2)
            if isempty(stencil)
                # Isolated agglomerate: every edge neighbor was rejected by the
                # TE barrier or the fold test (degenerate sliver panels whose
                # normal disagrees with all neighbors). No stencil can exist at
                # any growth depth, so use a zero surface gradient rather than
                # abort: such panels carry negligible area by construction.
                @warn "basis=:quad found an isolated agglomerate; using zero μ gradient there" agglomerate=ag members=members maxlog=8
                gx = 0.0; gy = 0.0; gz = 0.0
            else
                throw(ArgumentError("basis=:quad failed to reconstruct a finite μ gradient for agglomerate $(ag) (members=$(members), stencil_size=$(length(stencil)), cond=$(cnd), quad_grow=$(grow), quad_grow_stop=$(grow_stop), quad_grow_max_depth=$(grow_max_depth))."))
            end
        end

        # Scatter agglomerate ∇μ to member triangles, projected onto each tangent
        # plane, matching pass-1 sign convention (local_grad = -scale·∇μ).
        m1 = ag
        m2 = partner[ag]                      # 0 if singleton
        for t in (m1, m2)
            t == 0 && continue
            nd = gx*normals[1,t] + gy*normals[2,t] + gz*normals[3,t]
            local_grad[1,t] = -scale*(gx - nd*normals[1,t])
            local_grad[2,t] = -scale*(gy - nd*normals[2,t])
            local_grad[3,t] = -scale*(gz - nd*normals[3,t])
        end
    end

    return local_grad
end

"""
    compute_surface_velocity_gradient!(grad_u, u, controlpoints, normals,
                                       cells, neighbors, te_info; optargs...)

Compute a panel-centered surface gradient of a vector velocity field `u`.
The output layout is `grad_u[k, l, i] = ∂u_k/∂x_l` at panel `i`, with each
row constrained to the local tangent plane. The stencil, trailing-edge
isolation, and optional bad-panel fallback match `compute_mu_gradient!`.
"""
function compute_surface_velocity_gradient!(grad_u::AbstractArray{Float64,3},
                                            u::AbstractMatrix{Float64},
                                            controlpoints::AbstractMatrix{Float64},
                                            normals::AbstractMatrix{Float64},
                                            cells::AbstractMatrix{Int},
                                            neighbors::AbstractMatrix{Int},
                                            te_info::AbstractMatrix{Int};
                                            bad_panel_mask::Union{Nothing,AbstractVector{Bool}}=nothing,
                                            nodes::Union{Nothing,AbstractMatrix{Float64}}=nothing,
                                            grad_mu_options=(;))
    size(u, 1) == 3 || throw(ArgumentError("u must be a 3 × ncells matrix."))
    size(grad_u, 1) == 3 && size(grad_u, 2) == 3 && size(grad_u, 3) == size(u, 2) ||
        throw(ArgumentError("grad_u must have size 3 × 3 × ncells."))
    opts = _normalize_grad_mu_options(grad_mu_options;
        default_basis=isnothing(nodes) ? :tri : :quad)
    if bad_panel_mask === nothing && opts.tri_robust
        isnothing(nodes) && throw(ArgumentError(
            "grad_mu_options.tri_robust=true requires nodes for aspect-ratio masking."))
        mask = panel_aspect_ratio_mask(nodes, cells;
            threshold=opts.tri_robust_ar_threshold)
        any(mask) && (bad_panel_mask = mask)
    end

    fill!(grad_u, 0.0)
    scratch = zeros(Float64, 3, size(u, 2))
    @inbounds for k in 1:3
        fill!(scratch, 0.0)
        # compute_mu_gradient! stores -∇μ by convention. Use scale=-1 to recover
        # the actual surface gradient of the velocity component.
        _compute_mu_gradient_masked!(scratch, controlpoints, normals, cells, neighbors,
            view(u, k, :), te_info;
            scale=-1.0,
            bad_panel_mask=bad_panel_mask,
            nodes=nodes,
            grad_mu_options=opts)
        for i in axes(u, 2)
            grad_u[k, 1, i] = scratch[1, i]
            grad_u[k, 2, i] = scratch[2, i]
            grad_u[k, 3, i] = scratch[3, i]
        end
    end

    return grad_u
end

################################################################################
# PRESSURE FIELDS
################################################################################
"""
    calcfield_P!(out, body, Us, Uinf, rho, phi_dot; correct_kuttacondition=false, clip=nothing)

Compute and store the dimensional gauge pressure field using the Bernoulli
equation:  ``P = \\frac{1}{2} \\rho (U_\\infty^2 - U^2) - \\rho \\frac{\\partial \\phi}{\\partial t}``

`phi_dot` is a per-panel ``\\partial \\phi / \\partial t`` vector; pass
`nothing` to ignore the unsteady term. The field is calculated in-place and
added to `out` (hence, make sure that `out` starts with all zeroes).
Set `correct_kuttacondition=true` to opt into heuristic averaging of pressures
on paired trailing-edge panels.
"""
function calcfield_P!(out::Arr1,
                       body::Union{NonLiftingBody, AbstractLiftingBody},
                       Us::Arr2, Uinf::Number, rho::Number,
                       phi_dot::Union{Nothing, AbstractVector};
                       correct_kuttacondition=false,
                       clip::Union{Nothing, Function}=nothing,
                       ) where {Arr1<:AbstractArray{<:Number,1},
                                Arr2<:AbstractArray{<:Number,2}}

    half_rho = rho / 2
    Uinf2 = Uinf^2

    # Steady Bernoulli: P = 0.5*rho*(Uinf^2 - U^2)
    for (i, U) in enumerate(eachcol(Us))
        out[i] += half_rho * (Uinf2 - norm(U)^2)
    end

    if !isnothing(phi_dot)
        # Unsteady Bernoulli term: -rho * ∂φ/∂t
        for i in eachindex(out)
            out[i] -= rho * phi_dot[i]
        end
    end

    # Kutta-condition correction bringing the pressure on both sides of the TE
    # to be equal (average between upper and lower)
    if correct_kuttacondition && typeof(body) <: AbstractLiftingBody

        # Iterate over TE panels
        for shedding in body.shedding
            for (pi, nia, nib, pj, nja, njb) in eachcol(shedding)
                if pj != -1
                    ave = (out[pi] + out[pj]) / 2
                    out[pi] = ave
                    out[pj] = ave
                end
            end
        end

    end

    # Clip values if requested
    if !isnothing(clip)
        for (i, P) in enumerate(out)
            out[i] = clip(P)
        end
    end

    return out
end

"""
    calcfield_Cp(Ps, uinf, rho)

Compute and return the pressure coefficient field from the dimensional gauge
pressure vector `Ps`.
"""
function calcfield_Cp(Ps::AbstractVector, uinf::Number, rho::Number)
    qinf = 0.5 * rho * uinf^2
    return Ps ./ qinf
end

################################################################################
# FORCE FIELDS
################################################################################
"""
    calcfield_F!(out, body, areas, normals, Ps; correct_kuttacondition=true)

Compute and store distributed surface forces from the current gauge pressure
field:  ``F = -P \\cdot A \\cdot \\hat{n}``
"""
function calcfield_F!(out::Arr0, body::AbstractBody,
                         areas::Arr1, normals::Arr2, Ps::Arr3;
                         correct_kuttacondition=true,
                         ) where {   Arr0<:AbstractArray{<:Number,2},
                                     Arr1<:AbstractArray{<:Number,1},
                                     Arr2<:AbstractArray{<:Number,2},
                                     Arr3<:AbstractArray{<:Number,1}}

    # Error cases
    @assert size(out, 1)==3 && size(out, 2)==body.ncells ""*
        "Invalid `out` matrix."*
        " Expected size $((3, body.ncells)); got $(size(out))."
    @assert length(areas)==body.ncells ""*
        "Invalid `areas` vector."*
        " Expected length $(body.ncells); got $(length(areas))."
    @assert size(normals, 1)==3 && size(normals, 2)==body.ncells ""*
        "Invalid `normals` matrix."*
        " Expected size $((3, body.ncells)); got $(size(normals))."
    @assert length(Ps)==body.ncells ""*
        "Invalid `Ps` vector."*
        " Expected length $(body.ncells); got $(length(Ps))."

    for (i, (P, area, normal)) in enumerate(zip(Ps, areas, eachcol(normals)))
        val = -P * area
        out[1, i] += val*normal[1]
        out[2, i] += val*normal[2]
        out[3, i] += val*normal[3]
    end

    # Kutta-condition correction bringing the pressure on both sides of the TE
    # to be equal (average between upper and lower)
    # NOTE: This overwrites any previous force value instead of accumulating it
    if correct_kuttacondition && typeof(body) <: AbstractLiftingBody

        # Iterate over TE panels
        for shedding in body.shedding
            for (pi, nia, nib, pj, nja, njb) in eachcol(shedding)

                if pj != -1
                    aveP = (Ps[pi] + Ps[pj]) / 2

                    out[1, pi] = -aveP * areas[pi] * normals[1, pi]
                    out[2, pi] = -aveP * areas[pi] * normals[2, pi]
                    out[3, pi] = -aveP * areas[pi] * normals[3, pi]
                    out[1, pj] = -aveP * areas[pj] * normals[1, pj]
                    out[2, pj] = -aveP * areas[pj] * normals[2, pj]
                    out[3, pj] = -aveP * areas[pj] * normals[3, pj]

                end
            end
        end
    end

    return out
end

"""
    calcfield_sectionalforce!(outf, outpos, body, controlpoints, Fs; dimspan=2, dimchord=1, spandirection=[0, 1, 0], ...)
    calcfield_sectionalforce!(outf, outpos, body; optargs...)
    calcfield_sectionalforce(body; optargs...)

Integrate distributed panel forces into sectional loads along a chosen spanwise
direction.
"""
function calcfield_sectionalforce!(outf::Arr0, outpos::Arr1,
                                    body::Union{NonLiftingBody, AbstractLiftingBody},
                                    controlpoints::Arr2, Fs::Arr3;
                                    dimspan=2, dimchord=1,
                                    spandirection=[0, 1, 0],
                                    fieldname="sectionalforce", addfield=true
                                    ) where {   Arr0<:AbstractArray{<:Number,2},
                                                Arr1<:AbstractArray{<:Number,1},
                                                Arr2<:AbstractArray{<:Number,2},
                                                Arr3<:AbstractArray{<:Number,2}}



    lin, gdims = get_linearindex(body)      # LinearIndex and grid dimensions

    # Error cases
    @assert size(outf, 1)==3 && size(outf, 2)==gdims[dimspan] ""*
        "Invalid `outf` matrix."*
        " Expected size $((3, gdims[dimspan])); got $(size(outf))."
    @assert length(outpos)==gdims[dimspan] ""*
        "Invalid `outpos` matrix."*
        " Expected length $(gdims[dimspan]); got $(length(outpos))."
    @assert size(controlpoints, 1)==3 && size(controlpoints, 2)==body.ncells ""*
        "Invalid `controlpoints` matrix."*
        " Expected size $((3, body.ncells)); got $(size(controlpoints))."
    @assert size(Fs, 1)==3 && size(Fs, 2)==body.ncells ""*
        "Invalid `Fs` matrix."*
        " Expected size $((3, body.ncells)); got $(size(Fs))."

    # Pre-allocate memory
    coor = ones(Int, 3)                     # Cartesian coordinates (indices)
    lincoors = zeros(Int, gdims[dimchord])  # Linear coordinate (index)
    outf .= 0

    # Integrate force in the chordwise direction along the span
    for j in 1:gdims[dimspan] # Iterate over span

        for i in 1:gdims[dimchord] # Iterate over chord

            coor[dimchord] = i
            coor[dimspan] = j
            lincoors[i] = lin[coor...]

            # Add force to this section
            outf[1, j] += Fs[1, lincoors[i]]
            outf[2, j] += Fs[2, lincoors[i]]
            outf[3, j] += Fs[3, lincoors[i]]

        end

        # Calculate span position of this section
        spanpos = mean(dot(spandirection, Xcp)
                        for Xcp in eachcol(view(controlpoints, :, lincoors)))
        outpos[j] = spanpos

    end

    # Convert force to be per unit span
    for j in 1:gdims[dimspan] # Iterate over span
        deltapos =  j==1 ?              outpos[j+1]-outpos[j] :
                    j==length(outpos) ? outpos[j]-outpos[j-1] :
                                        (outpos[j+1]-outpos[j-1])/2

        outf[:, j] /= abs(deltapos)
    end

    # Save field in body
    # if addfield
    #     add_field(body, fieldname, "vector", eachcol(outf), "system")
    #     add_field(body, fieldname*"-pos", "vector", eachcol(outpos), "system")
    # end

    return outf, outpos
end

"""
    calcfield_sectionalforce(args...; optargs...)

Similar to [`calcfield_sectionalforce!`](@ref) but without in-place calculation
(`outFs` nor `outpos` are needed).
"""
function calcfield_sectionalforce(body::Union{NonLiftingBody, AbstractLiftingBody}, args...;
                                                        dimspan=2, optargs...)

    lin, gdims = get_linearindex(body)      # LinearIndex and grid dimensions

    outFs = zeros(3, gdims[dimspan])
    outpos = zeros(gdims[dimspan])

    return calcfield_sectionalforce!(outFs, outpos, body, args...;
                                                    dimspan=dimspan, optargs...)
end

"""
    calcfield_Ftot!(out, body, Fs; fieldname="Ftot")
    calcfield_Ftot(body, args...; optargs...)

Integrate the stored distributed force field into a total force vector.
"""
function calcfield_Ftot!(out::AbstractVector, body::AbstractBody,
                            Fs::AbstractMatrix; fieldname="Ftot", addfield=true)

    # Error case
    @assert length(out)==3 ""*
        "Invalid `out` vector. Expected length $(3); got $(length(out))."

    for i in 1:3
        out[i] += sum(view(Fs, i, :))
    end

    # Save field in body
    # if addfield
    #     add_field(body, fieldname, "vector", out, "system")
    # end

    return out
end

"""
    calcfield_Ftot(body, args...; optargs...) = calcfield_Ftot!(zeros(3), body, args...; optargs...)

Similar to [`calcfield_Ftot!`](@ref) but without in-place calculation (`out` is
not needed).
"""
calcfield_Ftot(body, args...; optargs...) = calcfield_Ftot!(zeros(3), body, args...; optargs...)

"""
    calcfield_LDS!(out, body, Fs, Lhat, Dhat, Shat; addfield=true)
    calcfield_LDS!(out, body, Lhat, Dhat, Shat; optargs...)
    calcfield_LDS!(out, body, Lhat, Dhat; optargs...)
    calcfield_LDS(body, args...; optargs...)

Resolve the total force into lift, drag, and side-force components.
"""
function calcfield_LDS!(out::AbstractMatrix, body::AbstractBody,
                        Fs::AbstractMatrix,
                        Lhat::AbstractVector, Dhat::AbstractVector,
                        Shat::AbstractVector;
                        addfield=true)
    # Error case
    @assert size(out, 1)==3 && size(out, 2)==3 ""*
        "Invalid `out` matrix. Expected size $((3, 3)); got $(size(out))."
    @assert abs(norm(Lhat) - 1) <= 2*eps() ""*
        "Lhat=$(Lhat) is not a unitary vector"
    @assert abs(norm(Dhat) - 1) <= 2*eps() ""*
        "Dhat=$(Dhat) is not a unitary vector"
    @assert abs(norm(Shat) - 1) <= 2*eps() ""*
        "Shat=$(Shat) is not a unitary vector"

    # Calculate Ftot (integrated force)
    for i in 1:3
        out[i, 3] += sum(view(Fs, i, :))
    end

    # Project Ftot in each direction
    out[:, 1] = Lhat
    out[:, 1] *= dot(view(out, :, 3), Lhat)
    out[:, 2] = Dhat
    out[:, 2] *= dot(view(out, :, 3), Dhat)
    aux = dot(view(out, :, 3), Shat)
    out[:, 3] = Shat
    out[:, 3] *= aux

    return out
end

"""
    calcfield_LDS(body, args...; optargs...) = calcfield_LDS!(zeros(3, 3), body, args...; optargs...)

Similar to [`calcfield_LDS!`](@ref) but without in-place calculation (`out` is
not needed).
"""
calcfield_LDS(body, args...; optargs...) = calcfield_LDS!(zeros(3, 3), body, args...; optargs...)









################################################################################
# MOMENT FIELDS
################################################################################
"""
    calcfield_Mtot!(out, body, Xac, controlpoints, Fs; fieldname="Mtot", addfield=true)
    calcfield_Mtot(body, args...; optargs...)

Integrate the stored distributed force field into a total moment about a
reference point.
"""
function calcfield_Mtot!(out::AbstractVector, body::AbstractBody,
                            Xac::AbstractVector, controlpoints::AbstractMatrix,
                            Fs::AbstractMatrix;
                            fieldname="Mtot", addfield=true)
    # Error case
    @assert length(out)==3 ""*
        "Invalid `out` vector. Expected length 3; got $(length(out))."
    @assert length(Xac)==3 ""*
        "Invalid `Xac` vector. Expected length 3; got $(length(Xac))."
    @assert size(controlpoints, 1)==3 && size(controlpoints, 2)==body.ncells ""*
        "Invalid `controlpoints` matrix."*
        " Expected size $((3, body.ncells)); got $(size(controlpoints))."
    @assert size(Fs, 1)==3 && size(Fs, 2)==body.ncells ""*
        "Invalid `Fs` matrix."*
        " Expected size $((3, body.ncells)); got $(size(Fs))."

    # Calculate Mtot (integrated moment)
    for (X, F) in zip(eachcol(controlpoints), eachcol(Fs))
        out[1] += (X[2] - Xac[2])*F[3] - (X[3] - Xac[3])*F[2]
        out[2] += (X[3] - Xac[3])*F[1] - (X[1] - Xac[1])*F[3]
        out[3] += (X[1] - Xac[1])*F[2] - (X[2] - Xac[2])*F[1]
    end

    # Save field in body
    # if addfield
    #     add_field(body, fieldname, "vector", out, "system")
    # end

    return out
end

"""
    calcfield_Mtot(body, args...; optargs...) = calcfield_Mtot!(zeros(3), body, args...; optargs...)

Similar to [`calcfield_Mtot!`](@ref) but without in-place calculation (`out` is
not needed).
"""
calcfield_Mtot(body, args...; optargs...) = calcfield_Mtot!(zeros(3), body, args...; optargs...)



"""
    calcfield_lmn!(out, body, Xac, controlpoints, Fs, lhat, mhat, nhat; addfield=true)
    calcfield_lmn!(out, body, Xac, lhat, mhat; optargs...)
    calcfield_lmn(body, args...; optargs...)

Resolve the total moment into user-supplied rolling, pitching, and yaw axes.
"""
function calcfield_lmn!(out::AbstractMatrix, body::AbstractBody,
                        Xac::AbstractVector, controlpoints::AbstractMatrix,
                        Fs::AbstractMatrix,
                        lhat::AbstractVector, mhat::AbstractVector,
                        nhat::AbstractVector;
                        addfield=true)
    # Error case
    @assert size(out, 1)==3 && size(out, 2)==3 ""*
        "Invalid `out` matrix. Expected size $((3, 3)); got $(size(out))."
    @assert length(Xac)==3 ""*
        "Invalid `Xac` vector. Expected length 3; got $(length(Xac))."
    @assert size(controlpoints, 1)==3 && size(controlpoints, 2)==body.ncells ""*
        "Invalid `controlpoints` matrix."*
        " Expected size $((3, body.ncells)); got $(size(controlpoints))."
    @assert size(Fs, 1)==3 && size(Fs, 2)==body.ncells ""*
        "Invalid `Fs` matrix."*
        " Expected size $((3, body.ncells)); got $(size(Fs))."
    @assert abs(norm(lhat) - 1) <= 2*eps() ""*
        "lhat=$(lhat) is not a unitary vector"
    @assert abs(norm(mhat) - 1) <= 2*eps() ""*
        "mhat=$(mhat) is not a unitary vector"
    @assert abs(norm(nhat) - 1) <= 2*eps() ""*
        "nhat=$(nhat) is not a unitary vector"

    # Calculate Mtot (integrated moment)
    for (X, F) in zip(eachcol(controlpoints), eachcol(Fs))
        out[1, 3] += (X[2] - Xac[2])*F[3] - (X[3] - Xac[3])*F[2]
        out[2, 3] += (X[3] - Xac[3])*F[1] - (X[1] - Xac[1])*F[3]
        out[3, 3] += (X[1] - Xac[1])*F[2] - (X[2] - Xac[2])*F[1]
    end

    # Project Mtot in each direction
    out[:, 1] = lhat
    out[:, 1] *= dot(view(out, :, 3), lhat)
    out[:, 2] = mhat
    out[:, 2] *= dot(view(out, :, 3), mhat)
    aux = dot(view(out, :, 3), nhat)
    out[:, 3] = nhat
    out[:, 3] *= aux

    # Save field in body
    # if addfield
    #     add_field(body, "Mroll", "vector", view(out, :, 1), "system")
    #     add_field(body, "Mpitch", "vector", view(out, :, 2), "system")
    #     add_field(body, "Myaw", "vector", view(out, :, 3), "system")
    # end

    return out
end


"""
    calcfield_lmn(body, args...; optargs...) = calcfield_lmn!(zeros(3, 3), body, args...; optargs...)

Similar to [`calcfield_lmn!`](@ref) but without in-place calculation (`out` is
not needed).
"""
calcfield_lmn(body, args...; optargs...) = calcfield_lmn!(zeros(3, 3), body, args...; optargs...)
