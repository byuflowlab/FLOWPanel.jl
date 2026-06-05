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
        gradient_robust=false,
        gradient_ar_threshold=10.0,
        gradient_quad_consistent=true,
    )

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
            mask = nothing
            if gradient_robust
                m = panel_aspect_ratio_mask(body.nodes, body.cells;
                                            threshold=gradient_ar_threshold)
                any(m) && (mask = m)
            end
            compute_mu_gradient!(body.velocity, body.controlpoints, body.normals,
                body.cells,
                body.neighbor,
                view(body.strength, :, get_Gammai(body)),
                view(body.shedding_full, 1:2, :), scale=0.5,
                bad_panel_mask=mask,
                nodes=gradient_quad_consistent ? body.nodes : nothing)
            
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
`compute_mu_gradient!` cannot reliably recover the in-plane gradient. Pass
the returned mask as `bad_panel_mask` to `compute_mu_gradient!` to trigger
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
                         scale=0.5, bad_panel_mask=nothing,
                         bfs_target_healthy=6, bfs_max_depth=4,
                         nodes=nothing, quad_consistent=!isnothing(nodes))

Compute the surface gradient of doublet or vortex-ring strength using a local
least-squares stencil, with one-sided handling at trailing-edge panels.

When `bad_panel_mask` is provided, panels flagged `true` skip the 1-ring
stencil and instead inherit a gradient computed from a BFS-gathered set of
healthy (unmasked) panels reached by walking the `neighbors` graph and
*stepping through* flagged panels. This makes the reconstruction robust to
clusters of degenerate panels (e.g. cap-corner slivers on capped wing meshes).
TE exclusion mirrors the pass-1 behavior: when the bad panel is itself a TE
panel, the cross-TE neighbor is rejected as the first BFS step.

When `nodes` is supplied and `quad_consistent=true`, mutually paired
triangles across their longest shared edge are treated as the two halves of a
logical quadrilateral patch: their reconstructed contributions are averaged
and projected back onto each triangle tangent plane. This removes alternating
triangle-stencil artifacts on structured lofted lifting surfaces while leaving
unpaired or strongly folded neighbors unchanged.
"""
function compute_mu_gradient!(grad_mu,
                            controlpoints::AbstractMatrix{Float64},
                            normals::AbstractMatrix{Float64},
                            cells::AbstractMatrix{Int},
                            neighbors::AbstractMatrix{Int},
                            mu::AbstractVector{Float64},
                            te_info::AbstractMatrix{Int};
                            scale=0.5,
                            bad_panel_mask::Union{Nothing,AbstractVector{Bool}}=nothing,
                            bfs_target_healthy::Int=6,
                            bfs_max_depth::Int=4,
                            nodes::Union{Nothing,AbstractMatrix{Float64}}=nothing,
                            quad_consistent::Bool=!isnothing(nodes),
                            quad_normal_dot_min::Float64=cos(pi / 4),
                        )

    _, M = size(cells)
    local_grad = zeros(Float64, 3, M)

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

            fill!(visited, false)
            empty!(frontier); empty!(next_frontier); empty!(healthy)
            visited[i] = true
            push!(frontier, i)

            is_te_i = te_info[1, i] > 0 && te_info[2, i] > 0
            te_v1 = is_te_i ? cells[te_info[1, i], i] : -1
            te_v2 = is_te_i ? cells[te_info[2, i], i] : -1

            depth = 0
            while !isempty(frontier) && length(healthy) < bfs_target_healthy && depth < bfs_max_depth
                depth += 1
                for p in frontier
                    for k in 1:3
                        q = neighbors[k, p]
                        q <= 0 && continue
                        visited[q] && continue
                        # TE exclusion mirrors the pass-1 one-sided rule: only when
                        # the bad panel itself is a TE panel and q sits across its TE edge.
                        if p == i && is_te_i
                            has_v1 = (cells[1, q] == te_v1) || (cells[2, q] == te_v1) || (cells[3, q] == te_v1)
                            has_v2 = (cells[1, q] == te_v2) || (cells[2, q] == te_v2) || (cells[3, q] == te_v2)
                            if has_v1 && has_v2
                                continue
                            end
                        end
                        visited[q] = true
                        if is_bad[q]
                            push!(next_frontier, q)  # walk through, do not use in fit
                        else
                            push!(healthy, q)
                            push!(next_frontier, q)  # keep walking for multi-ring coverage
                        end
                    end
                end
                frontier, next_frontier = next_frontier, frontier
                empty!(next_frontier)
            end

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

    if quad_consistent && nodes !== nothing
        _enforce_quad_consistent_gradient!(local_grad, nodes, cells, neighbors,
            normals, te_info; normal_dot_min=quad_normal_dot_min)
    end

    grad_mu .+= local_grad

    return grad_mu
end

function _edge_vertices(cells::AbstractMatrix{Int}, edge_i::Int, cell_i::Int)
    edge_i == 1 && return cells[1, cell_i], cells[2, cell_i]
    edge_i == 2 && return cells[2, cell_i], cells[3, cell_i]
    return cells[3, cell_i], cells[1, cell_i]
end

function _edge_length2(nodes::AbstractMatrix{Float64},
                       cells::AbstractMatrix{Int},
                       edge_i::Int,
                       cell_i::Int)
    a, b = _edge_vertices(cells, edge_i, cell_i)
    dx = nodes[1, a] - nodes[1, b]
    dy = nodes[2, a] - nodes[2, b]
    dz = nodes[3, a] - nodes[3, b]
    return dx*dx + dy*dy + dz*dz
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

function _quad_candidate_neighbor(nodes::AbstractMatrix{Float64},
                                  cells::AbstractMatrix{Int},
                                  neighbors::AbstractMatrix{Int},
                                  normals::AbstractMatrix{Float64},
                                  te_info::AbstractMatrix{Int},
                                  i::Int;
                                  normal_dot_min::Float64)
    best = 0
    best_l2 = -Inf
    @inbounds for edge_i in 1:3
        j = neighbors[edge_i, i]
        j > 0 || continue
        _shares_te_edge(cells, te_info, i, j) && continue
        ndot = normals[1, i]*normals[1, j] + normals[2, i]*normals[2, j] + normals[3, i]*normals[3, j]
        ndot >= normal_dot_min || continue
        l2 = _edge_length2(nodes, cells, edge_i, i)
        if l2 > best_l2
            best_l2 = l2
            best = j
        end
    end
    return best
end

function _enforce_quad_consistent_gradient!(grad::AbstractMatrix{Float64},
                                            nodes::AbstractMatrix{Float64},
                                            cells::AbstractMatrix{Int},
                                            neighbors::AbstractMatrix{Int},
                                            normals::AbstractMatrix{Float64},
                                            te_info::AbstractMatrix{Int};
                                            normal_dot_min::Float64=cos(pi / 4))
    M = size(cells, 2)
    paired = falses(M)
    @inbounds for i in 1:M
        paired[i] && continue
        j = _quad_candidate_neighbor(nodes, cells, neighbors, normals, te_info, i;
                                     normal_dot_min=normal_dot_min)
        j > 0 || continue
        paired[j] && continue
        _quad_candidate_neighbor(nodes, cells, neighbors, normals, te_info, j;
                                 normal_dot_min=normal_dot_min) == i || continue

        gx = 0.5 * (grad[1, i] + grad[1, j])
        gy = 0.5 * (grad[2, i] + grad[2, j])
        gz = 0.5 * (grad[3, i] + grad[3, j])

        ni_dot = gx*normals[1, i] + gy*normals[2, i] + gz*normals[3, i]
        nj_dot = gx*normals[1, j] + gy*normals[2, j] + gz*normals[3, j]
        grad[1, i] = gx - ni_dot * normals[1, i]
        grad[2, i] = gy - ni_dot * normals[2, i]
        grad[3, i] = gz - ni_dot * normals[3, i]
        grad[1, j] = gx - nj_dot * normals[1, j]
        grad[2, j] = gy - nj_dot * normals[2, j]
        grad[3, j] = gz - nj_dot * normals[3, j]

        paired[i] = true
        paired[j] = true
    end

    return grad
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
                                            bfs_target_healthy::Int=6,
                                            bfs_max_depth::Int=4,
                                            nodes::Union{Nothing,AbstractMatrix{Float64}}=nothing,
                                            quad_consistent::Bool=!isnothing(nodes))
    size(u, 1) == 3 || throw(ArgumentError("u must be a 3 × ncells matrix."))
    size(grad_u, 1) == 3 && size(grad_u, 2) == 3 && size(grad_u, 3) == size(u, 2) ||
        throw(ArgumentError("grad_u must have size 3 × 3 × ncells."))

    fill!(grad_u, 0.0)
    scratch = zeros(Float64, 3, size(u, 2))
    @inbounds for k in 1:3
        fill!(scratch, 0.0)
        # compute_mu_gradient! stores -∇μ by convention. Use scale=-1 to recover
        # the actual surface gradient of the velocity component.
        compute_mu_gradient!(scratch, controlpoints, normals, cells, neighbors,
            view(u, k, :), te_info;
            scale=-1.0,
            bad_panel_mask=bad_panel_mask,
            bfs_target_healthy=bfs_target_healthy,
            bfs_max_depth=bfs_max_depth,
            nodes=nodes,
            quad_consistent=quad_consistent)
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
    calcfield_P!(out, body, Us, Uinf, rho, phi_dot; correct_kuttacondition=true, clip=nothing)

Compute and store the dimensional gauge pressure field using the Bernoulli
equation:  ``P = \\frac{1}{2} \\rho (U_\\infty^2 - U^2) - \\rho \\frac{\\partial \\phi}{\\partial t}``

`phi_dot` is a per-panel ``\\partial \\phi / \\partial t`` vector; pass
`nothing` to ignore the unsteady term. The field is calculated in-place and
added to `out` (hence, make sure that `out` starts with all zeroes).
"""
function calcfield_P!(out::Arr1,
                       body::Union{NonLiftingBody, AbstractLiftingBody},
                       Us::Arr2, Uinf::Number, rho::Number,
                       phi_dot::Union{Nothing, AbstractVector};
                       correct_kuttacondition=true,
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
