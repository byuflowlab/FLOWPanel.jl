"""
    ParticleBodyOverlapPolicy(; action=:off, core_ratio=1, every=1)

Pre-wake-influence safety policy for vortex particles near panel bodies.
`action` is `:off`, `:warn`, or `:error`; an overlap is reported when a
particle-to-triangle distance is no greater than `core_ratio*sigma`.
`every` controls the timestep sampling interval.
"""
struct ParticleBodyOverlapPolicy
    action::Symbol
    core_ratio::Float64
    every::Int
    function ParticleBodyOverlapPolicy(; action::Symbol=:off,
            core_ratio::Real=1, every::Integer=1)
        action in (:off, :warn, :error) || throw(ArgumentError(
            "particle/body overlap action must be :off, :warn, or :error"))
        core_ratio > 0 || throw(ArgumentError(
            "particle/body overlap core_ratio must be positive"))
        every >= 1 || throw(ArgumentError(
            "particle/body overlap every must be at least 1"))
        new(action, Float64(core_ratio), Int(every))
    end
end

struct ParticleBodyOverlapReport
    checked_particles::Int
    checked_fields::Int
    triangle_evaluations::Int
    overlap_count::Int # particle/body pairs, not unique particles
    min_distance::Float64
    min_distance_over_sigma::Float64
    body_index::Int
    wake_index::Int
    particle_index::Int
    particle_position::SVector{3,Float64}
    particle_sigma::Float64
    threshold_ratio::Float64
end

struct ParticleBodyOverlapError <: Exception
    report::ParticleBodyOverlapReport
end

function Base.showerror(io::IO, err::ParticleBodyOverlapError)
    r = err.report
    print(io, "particle/core overlaps a panel body: distance/sigma=",
        r.min_distance_over_sigma, " <= ", r.threshold_ratio,
        " (body=", r.body_index, ", wake=", r.wake_index,
        ", particle=", r.particle_index, ", x=", r.particle_position, ")")
end

struct _OverlapBVHNode
    bmin::SVector{3,Float64}
    bmax::SVector{3,Float64}
    left::Int
    right::Int
    first::Int
    last::Int
end

struct _OverlapBVH
    body::Any
    permutation::Vector{Int}
    nodes::Vector{_OverlapBVHNode}
end

@inline function _triangle_vertices(body, panel)
    c = view(body.cells, :, panel)
    return SVector{3,Float64}(view(body.nodes, :, c[1])),
        SVector{3,Float64}(view(body.nodes, :, c[2])),
        SVector{3,Float64}(view(body.nodes, :, c[3]))
end

@inline function _triangle_bounds(body, panel)
    a, b, c = _triangle_vertices(body, panel)
    return min.(a, min.(b, c)), max.(a, max.(b, c))
end

@inline function _triangle_centroid_component(body, panel, axis)
    a, b, c = _triangle_vertices(body, panel)
    return (a[axis] + b[axis] + c[axis]) / 3
end

function _build_overlap_bvh(body; leaf_size::Int=8)
    permutation = collect(1:body.ncells)
    nodes = _OverlapBVHNode[]
    function build!(first, last)
        bmin = SVector(Inf, Inf, Inf)
        bmax = SVector(-Inf, -Inf, -Inf)
        for k in first:last
            lo, hi = _triangle_bounds(body, permutation[k])
            bmin = min.(bmin, lo)
            bmax = max.(bmax, hi)
        end
        inode = length(nodes) + 1
        push!(nodes, _OverlapBVHNode(bmin, bmax, 0, 0, first, last))
        if last - first + 1 > leaf_size
            extent = bmax - bmin
            axis = argmax(extent)
            sort!(view(permutation, first:last);
                by=p -> _triangle_centroid_component(body, p, axis))
            middle = (first + last) >>> 1
            left = build!(first, middle)
            right = build!(middle + 1, last)
            nodes[inode] = _OverlapBVHNode(bmin, bmax, left, right, 0, -1)
        end
        return inode
    end
    body.ncells > 0 && build!(1, body.ncells)
    return _OverlapBVH(body, permutation, nodes)
end

@inline function _point_aabb_distance2(p, lo, hi)
    d2 = 0.0
    for d in 1:3
        delta = p[d] < lo[d] ? lo[d] - p[d] :
            (p[d] > hi[d] ? p[d] - hi[d] : 0.0)
        d2 += delta * delta
    end
    return d2
end

# Closest point on triangle, from the Voronoi-region construction in
# Christer Ericson, Real-Time Collision Detection, section 5.1.5.
@inline function _point_triangle_distance2(p, a, b, c)
    ab = b - a
    ac = c - a
    ap = p - a
    d1 = LA.dot(ab, ap)
    d2 = LA.dot(ac, ap)
    d1 <= 0 && d2 <= 0 && return LA.dot(ap, ap)
    bp = p - b
    d3 = LA.dot(ab, bp)
    d4 = LA.dot(ac, bp)
    d3 >= 0 && d4 <= d3 && return LA.dot(bp, bp)
    vc = d1*d4 - d3*d2
    if vc <= 0 && d1 >= 0 && d3 <= 0
        v = d1 / (d1 - d3)
        q = a + v*ab
        return sum(abs2, p - q)
    end
    cp = p - c
    d5 = LA.dot(ab, cp)
    d6 = LA.dot(ac, cp)
    d6 >= 0 && d5 <= d6 && return LA.dot(cp, cp)
    vb = d5*d2 - d1*d6
    if vb <= 0 && d2 >= 0 && d6 <= 0
        w = d2 / (d2 - d6)
        q = a + w*ac
        return sum(abs2, p - q)
    end
    va = d3*d6 - d5*d4
    if va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        q = b + w*(c - b)
        return sum(abs2, p - q)
    end
    denom = inv(va + vb + vc)
    q = a + (vb*denom)*ab + (vc*denom)*ac
    return sum(abs2, p - q)
end

function _overlap_query(bvh::_OverlapBVH, p, radius)
    isempty(bvh.nodes) && return (Inf, 0)
    best2 = Inf
    evaluations = 0
    stack = Int[1]
    while !isempty(stack)
        inode = pop!(stack)
        node = bvh.nodes[inode]
        _point_aabb_distance2(p, node.bmin, node.bmax) > best2 && continue
        if node.left == 0
            for k in node.first:node.last
                panel = bvh.permutation[k]
                a, b, c = _triangle_vertices(bvh.body, panel)
                best2 = min(best2, _point_triangle_distance2(p, a, b, c))
                evaluations += 1
            end
        else
            left = bvh.nodes[node.left]
            right = bvh.nodes[node.right]
            dl = _point_aabb_distance2(p, left.bmin, left.bmax)
            dr = _point_aabb_distance2(p, right.bmin, right.bmax)
            # LIFO stack: push the farther node first so the nearer branch
            # tightens best2 before its sibling is reconsidered.
            if dl <= dr
                push!(stack, node.right, node.left)
            else
                push!(stack, node.left, node.right)
            end
        end
    end
    return best2, evaluations
end

"""
    particle_body_overlap(systems, wakes; core_ratio=1)

Run a device-safe overlap scan. Shared particle fields are inspected once;
device fields incur one contiguous live-prefix D2H mirror refresh. The report's
`overlap_count` counts particle/body pairs. Its minimum-distance fields are
true global nearest-triangle values (not threshold-limited); BVH
branch-and-bound avoids materializing a particle-by-panel Cartesian product.
"""
function particle_body_overlap(systems, wakes; core_ratio::Real=1)
    ratio_limit = Float64(core_ratio)
    ratio_limit > 0 || throw(ArgumentError("core_ratio must be positive"))
    systems_tuple = systems isa Tuple ? systems : (systems,)
    wakes_tuple = wakes isa Tuple ? wakes : (wakes,)
    bvhs = [_build_overlap_bvh(body) for body in systems_tuple]
    seen = IdDict{Any,Int}()
    checked_particles = 0
    checked_fields = 0
    triangle_evaluations = 0
    overlap_count = 0
    best_distance = Inf
    best_ratio = Inf
    best_body = best_wake = best_particle = 0
    best_position = SVector(Inf, Inf, Inf)
    best_sigma = NaN
    for (iwake, wake) in enumerate(wakes_tuple)
        wake isa PanelParticleWake || continue
        pfield = wake.pfield
        haskey(seen, pfield) && continue
        seen[pfield] = iwake
        checked_fields += 1
        host = _wake_monitor_host_pfield(pfield)
        checked_particles += host.np
        for iparticle in 1:host.np
            xraw = FLOWVPM.get_X(host, iparticle)
            x = SVector{3,Float64}(xraw)
            sigma = Float64(FLOWVPM.get_sigma(host, iparticle)[1])
            all(isfinite, x) || throw(ArgumentError(
                "particle $iparticle in wake $iwake has nonfinite position $x"))
            isfinite(sigma) && sigma > 0 || throw(ArgumentError(
                "particle $iparticle in wake $iwake has invalid sigma $sigma"))
            radius = ratio_limit * sigma
            for (ibody, bvh) in enumerate(bvhs)
                d2, neval = _overlap_query(bvh, x, radius)
                triangle_evaluations += neval
                isfinite(d2) || continue
                distance = sqrt(max(d2, 0.0))
                distance_ratio = distance / sigma
                if distance_ratio < best_ratio
                    best_distance = distance
                    best_ratio = distance_ratio
                    best_body = ibody
                    best_wake = iwake
                    best_particle = iparticle
                    best_position = x
                    best_sigma = sigma
                end
                distance_ratio <= ratio_limit && (overlap_count += 1)
            end
        end
    end
    return ParticleBodyOverlapReport(checked_particles, checked_fields,
        triangle_evaluations, overlap_count, best_distance, best_ratio,
        best_body, best_wake, best_particle, best_position, best_sigma,
        ratio_limit)
end

function check_particle_body_overlap!(policy::ParticleBodyOverlapPolicy,
        systems, wakes, i_step::Integer)
    policy.action === :off && return nothing
    mod(i_step, policy.every) == 0 || return nothing
    report = particle_body_overlap(systems, wakes;
        core_ratio=policy.core_ratio)
    report.overlap_count == 0 && return report
    if policy.action === :warn
        @warn "particle/core overlaps a panel body" i_step overlap_count=report.overlap_count min_distance=report.min_distance min_distance_over_sigma=report.min_distance_over_sigma body_index=report.body_index wake_index=report.wake_index particle_index=report.particle_index
    else
        throw(ParticleBodyOverlapError(report))
    end
    return report
end
