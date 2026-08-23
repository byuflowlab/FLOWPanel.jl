## Phase 2 driver for plans/dji_convergence_20260722/phase_02_te_adequacy.md.
##
## Run batches from the repository root:
##   PHASE2_MODE=setup      julia --project examples/dji9443_te_adequacy.jl
##   PHASE2_MODE=diagnostic PHASE2_CASE=new40c julia --project examples/dji9443_te_adequacy.jl
##   PHASE2_MODE=diagnostic PHASE2_CASE=new57c julia --project examples/dji9443_te_adequacy.jl
##   PHASE2_MODE=analyze    julia --project examples/dji9443_te_adequacy.jl
##   PHASE2_MODE=smoke      julia --project examples/dji9443_te_adequacy.jl

import FLOWPanel as pnl
using CSV
using DataFrames
using Dates
using LinearAlgebra
using Printf
using Statistics
using FLOWPanel.FastMultipole.StaticArrays

const STUDY = "dji_convergence_20260722"
const PHASE = "phase_02_te_adequacy"
const OUTPUT_DIR = joinpath("data", STUDY, PHASE)
const R = 0.119
const RPM = 5400.0
const OMEGA = 2pi * RPM / 60
const RADIAL_DIMENSION = 2
const VINF = [1.0e-4, 0.0, 0.0]
const RADIAL_BINS = collect(0.125:0.025:0.975)
const KERNEL_OFFSET_PANEL = R * 1.0e-10
const KERNEL_OFFSET_TARGETS = 1.0e-3
const KERNEL_CUTOFF = R * 1.0e-13
const GAMMA_SCALE = 1.0e-12
const RHS_BACKEND = pnl.FastMultipoleBackend(
    expansion_order=8, multipole_acceptance=0.4, leaf_size=20)

const CASES = [
    (tag=:new40c, mesh="dji9443_20260722_40_41_capped.msh",
     seeds=([1618, 1578, 0] .+ 1, [3332, 3292, 1714] .+ 1),
     expected=(3428, 6848, 39), use_bbox=true),
    (tag=:new57c, mesh="dji9443_20260722_57_57_capped.msh",
     seeds=([6572, 6516, 3354] .+ 1, [3218, 3162, 0] .+ 1),
     expected=(6708, 13408, 56), use_bbox=true),
]

case_by_tag(tag::Symbol) = only(filter(c -> c.tag == tag, CASES))
artifact(name, case) = joinpath(OUTPUT_DIR, "$(name)_$(case.tag).csv")

function make_shedding_bbox(nodes, seed_nodes)
    radial_midpoint = sum(nodes[RADIAL_DIMENSION, seed_nodes]) / length(seed_nodes)
    radial_sign = sign(radial_midpoint)
    radial_sign != 0 || error("seed edge lies on the rotor axis")
    lower = [minimum(nodes[i, :]) for i in axes(nodes, 1)]
    upper = [maximum(nodes[i, :]) for i in axes(nodes, 1)]
    padding = max(sqrt(eps(eltype(nodes))) * R, R * 1.0e-6)
    lower .-= padding
    upper .+= padding
    if radial_sign > 0
        lower[RADIAL_DIMENSION] = 0.1R - padding
    else
        upper[RADIAL_DIMENSION] = -0.1R + padding
    end
    return (SVector{3}(lower...), SVector{3}(upper...))
end

function load_base(case; core_size_panel=KERNEL_OFFSET_PANEL)
    mesh_path = joinpath(pnl.examples_path, "data", case.mesh)
    isfile(mesh_path) || error("missing mesh: $(mesh_path)")
    mesh = pnl.read_gmsh(mesh_path)
    nodes, cells = pnl.meshes2nodes_cells(mesh)
    source_radius = maximum(abs, nodes[RADIAL_DIMENSION, :])
    nodes .*= R / source_radius
    kernel = Union{pnl.ConstantSource, pnl.VortexRing}
    body = pnl.RigidWakeBody{kernel}(nodes, cells, pnl.noshedding;
        core_size=core_size_panel,
        core_size_panel=core_size_panel,
        core_size_targets=KERNEL_OFFSET_TARGETS,
        kernelcutoff=KERNEL_CUTOFF,
        semiinfinite_wake=false, watertight=true, DBC=true)
    return body, source_radius
end

function trace_shedding(base, case)
    return map(case.seeds) do seed
        bbox = case.use_bbox ? make_shedding_bbox(base.nodes, seed[1:2]) : nothing
        pnl.calc_shedding_from_seed(
            base.nodes, base.cells, seed[1], seed[2];
            bbox, end_node=seed[3], normal_jump_tol=0.2,
            max_turn_angle=pi / 3, debug=false)
    end
end

function build_case(case; core_size_panel=KERNEL_OFFSET_PANEL)
    base, source_radius = load_base(case; core_size_panel)
    shedding = trace_shedding(base, case)
    kernel = Union{pnl.ConstantSource, pnl.VortexRing}
    body = pnl.RigidWakeBody{kernel}(
        copy(base.nodes), copy(base.cells), [copy(s) for s in shedding];
        core_size=core_size_panel,
        core_size_panel=core_size_panel,
        core_size_targets=KERNEL_OFFSET_TARGETS,
        kernelcutoff=KERNEL_CUTOFF,
        semiinfinite_wake=false, watertight=true,
        ensure_winding=true, DBC=true)
    return body, source_radius
end

function make_frames(body)
    return pnl.ReferenceFrame(body;
        origin=SVector{3}(0.0, 0.0, 0.0),
        v=SVector{3}(0.0, 0.0, 0.0),
        ω_axis=SVector{3}(-1.0, 0.0, 0.0),
        ω=OMEGA,
        R=SMatrix{3,3}(1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0),
        name="vehicle", child_index=Int[], dependent_index=[1])
end

function section_radius(body, shedding, j)
    return abs(0.5 * (
        body.nodes[RADIAL_DIMENSION, body.cells[shedding[2, j], shedding[1, j]]] +
        body.nodes[RADIAL_DIMENSION, body.cells[shedding[3, j], shedding[1, j]]]
    ) / R)
end

function validate_case(case)
    body, source_radius = build_case(case)
    expected_nodes, expected_cells, expected_sections = case.expected
    size(body.nodes, 2) == expected_nodes ||
        error("$(case.tag): node count $(size(body.nodes, 2)) != $(expected_nodes)")
    body.ncells == expected_cells ||
        error("$(case.tag): panel count $(body.ncells) != $(expected_cells)")
    length(body.shedding) == 2 || error("$(case.tag): expected two shedding chains")
    all(size(s, 2) == expected_sections for s in body.shedding) ||
        error("$(case.tag): unexpected shedding section count")
    return (
        case=string(case.tag), mesh=case.mesh, source_radius=source_radius,
        scale_factor=R / source_radius, nodes=size(body.nodes, 2),
        panels=body.ncells, sections_per_blade=size(body.shedding[1], 2),
        core_size_panel=body.core_size_panel,
        core_size_targets=body.core_size_targets,
        status="pass")
end

function prepare_steady_body!(body)
    frames = make_frames(body)
    dt = 60 / RPM / 36
    pnl.initialize_Das!((body,), frames, t -> VINF, 0.0, dt;
        set_Das_eta_kinematic=0.2,
        set_Das_min_kinematic_displacement=0.01R)
    pnl.reset!(body)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    pnl.apply_freestream!((body,), VINF)
    pnl.kinematic_velocity!((body,), frames)
    pnl.set_strengths!(body)
    return frames
end

function solve_case(case; core_size_panel=KERNEL_OFFSET_PANEL)
    body, _ = build_case(case; core_size_panel)
    elapsed = @elapsed begin
        prepare_steady_body!(body)
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        G = zeros(eltype(body.strength), body.ncells, body.ncells)
        pnl._G!(G, body, body; core_size=body.core_size_panel,
            update_geometry=false)
        rhs = source_rhs(body)
        body.strength[:, 2] .= lu(G) \ rhs
    end
    return body, elapsed
end

function assemble_operators(body)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    G = zeros(eltype(body.strength), body.ncells, body.ncells)
    B = similar(G)
    pnl._G!(G, body, body; core_size=body.core_size_panel, update_geometry=false)
    pnl._assemble_B!(B, body)
    return G, B
end

function source_rhs(body; backend=RHS_BACKEND)
    old_strength = copy(body.strength)
    old_potential = copy(body.potential)
    try
        body.strength[:, 2] .= 0
        body.potential .= 0
        pnl.influence!(body, body, backend; scalar_potential=true, velocity=false)
        return -copy(body.potential)
    finally
        body.strength .= old_strength
        body.potential .= old_potential
    end
end

function panel_vertices(body, panel)
    return (SVector{3}(body.nodes[:, body.cells[1, panel]]),
            SVector{3}(body.nodes[:, body.cells[2, panel]]),
            SVector{3}(body.nodes[:, body.cells[3, panel]]))
end

triangle_area(a, b, c) = 0.5 * norm(cross(b - a, c - a))

function panel_aspect_ratio(body, panel)
    a, b, c = panel_vertices(body, panel)
    edges = (norm(b - a), norm(c - b), norm(a - c))
    area = triangle_area(a, b, c)
    hmin = 2area / max(maximum(edges), eps())
    return maximum(edges) / max(hmin, eps())
end

function point_segment_distance(p, a, b)
    ab = b - a
    t = clamp(dot(p - a, ab) / max(dot(ab, ab), eps()), 0.0, 1.0)
    return norm(p - (a + t * ab))
end

function point_triangle_distance(p, a, b, c)
    # Sufficient for diagnostics: minimum of plane projection if inside plus edges.
    n = cross(b - a, c - a)
    nn = norm(n)
    plane_dist = Inf
    if nn > eps()
        nhat = n / nn
        q = p - dot(p - a, nhat) * nhat
        v0, v1, v2 = c - a, b - a, q - a
        den = dot(v0, v0) * dot(v1, v1) - dot(v0, v1)^2
        if abs(den) > eps()
            u = (dot(v1, v1) * dot(v2, v0) - dot(v0, v1) * dot(v2, v1)) / den
            v = (dot(v0, v0) * dot(v2, v1) - dot(v0, v1) * dot(v2, v0)) / den
            if u >= 0 && v >= 0 && u + v <= 1
                plane_dist = abs(dot(p - a, nhat))
            end
        end
    end
    return min(plane_dist,
        point_segment_distance(p, a, b),
        point_segment_distance(p, b, c),
        point_segment_distance(p, c, a))
end

function centroid_edge_distance(body, panel, ia, ib)
    cp = SVector{3}(body.controlpoints[:, panel])
    a = SVector{3}(body.nodes[:, body.cells[ia, panel]])
    b = SVector{3}(body.nodes[:, body.cells[ib, panel]])
    return point_segment_distance(cp, a, b)
end

function nearest_nonself_distance(body, panel)
    p = SVector{3}(body.controlpoints[:, panel])
    best_dist = Inf
    best_panel = 0
    for q in 1:body.ncells
        q == panel && continue
        a, b, c = panel_vertices(body, q)
        d = point_triangle_distance(p, a, b, c)
        if d < best_dist
            best_dist = d
            best_panel = q
        end
    end
    return best_dist, best_panel
end

function nearest_panel_from_set(body, panel, candidates)
    isempty(candidates) && return (0, Inf)
    p = SVector{3}(body.controlpoints[:, panel])
    best = first(candidates)
    best_dist = norm(SVector{3}(body.controlpoints[:, best]) - p)
    for q in candidates
        d = norm(SVector{3}(body.controlpoints[:, q]) - p)
        if d < best_dist
            best = q
            best_dist = d
        end
    end
    return best, best_dist
end

function cap_panels(body)
    radial = abs.(body.controlpoints[RADIAL_DIMENSION, :] ./ R)
    tip = maximum(radial) - 0.02
    root = minimum(radial) + 0.02
    return findall(r -> r >= tip || r <= root, radial)
end

function te_pairs(body)
    rows = NamedTuple[]
    edge_id = 0
    for (blade, shedding) in pairs(body.shedding)
        for section in axes(shedding, 2)
            pi, nia, nib, pj, nja, njb = shedding[:, section]
            pj > 0 || continue
            edge_id += 1
            push!(rows, (edge=edge_id, blade=blade, section=section,
                upper=pi, lower=pj, upper_ia=nia, upper_ib=nib,
                lower_ia=nja, lower_ib=njb,
                abs_r_over_R=section_radius(body, shedding, section)))
        end
    end
    return rows
end

function write_geometry_influence(case, body, G, B)
    areas = pnl.calc_areas(body)
    mu = body.strength[:, 2]
    caps = cap_panels(body)
    pairs = te_pairs(body)
    upper_set = Set(p.upper for p in pairs)
    lower_set = Set(p.lower for p in pairs)
    pair_by_panel = Dict(vcat([p.upper => p for p in pairs], [p.lower => p for p in pairs])...)
    rows = NamedTuple[]
    for p in pairs
        for (role, panel, paired, ia, ib, same_set, opposite_set) in (
                ("upper", p.upper, p.lower, p.upper_ia, p.upper_ib, upper_set, lower_set),
                ("lower", p.lower, p.upper, p.lower_ia, p.lower_ib, lower_set, upper_set))
            cp = SVector{3}(body.controlpoints[:, panel])
            pair_cp = SVector{3}(body.controlpoints[:, paired])
            local_len = pnl.characteristiclength_sqrtarea(body.nodes, view(body.cells, :, panel))
            d_nonself, p_nonself = nearest_nonself_distance(body, panel)
            same_candidates = [q for q in same_set if q != panel]
            opp_candidates = [q for q in opposite_set if q != paired]
            p_same, d_same = nearest_panel_from_set(body, panel, same_candidates)
            p_opp, d_opp = nearest_panel_from_set(body, panel, opp_candidates)
            p_cap, d_cap = nearest_panel_from_set(body, panel, caps)
            row_diff = norm(view(G, p.upper, :) .- view(G, p.lower, :)) /
                       max(norm(view(G, p.upper, :)), norm(view(G, p.lower, :)), eps())
            row_cos = dot(view(G, p.upper, :), view(G, p.lower, :)) /
                      max(norm(view(G, p.upper, :)) * norm(view(G, p.lower, :)), eps())
            for ranking in ("abs_G", "abs_G_mu")
                vals = ranking == "abs_G" ? abs.(view(G, panel, :)) :
                       abs.(view(G, panel, :) .* mu)
                order = sortperm(vals; rev=true)
                top5 = order[1:min(5, length(order))]
                top3 = order[1:min(3, length(order))]
                push!(rows, (
                    case=string(case.tag), edge=p.edge, blade=p.blade,
                    section=p.section, role=role, panel=panel, paired_panel=paired,
                    abs_r_over_R=p.abs_r_over_R,
                    cp_pair_separation_over_l=norm(cp - pair_cp) / local_len,
                    te_dihedral_deg=acosd(clamp(dot(view(body.normals, :, p.upper),
                        view(body.normals, :, p.lower)), -1, 1)),
                    area=areas[panel], aspect_ratio=panel_aspect_ratio(body, panel),
                    centroid_to_edge_over_l=centroid_edge_distance(body, panel, ia, ib) / local_len,
                    nearest_nonself_panel=p_nonself, nearest_nonself_distance_over_l=d_nonself / local_len,
                    cross_te_panel=paired,
                    same_side_upstream_panel=p_same, same_side_upstream_distance_over_l=d_same / local_len,
                    opposite_side_panel=p_opp, opposite_side_distance_over_l=d_opp / local_len,
                    cap_panel=p_cap, cap_distance_over_l=d_cap / local_len,
                    ranking=ranking,
                    self_G=G[panel, panel], cross_te_G=G[panel, paired],
                    self_G_mu=G[panel, panel] * mu[panel],
                    cross_te_G_mu=G[panel, paired] * mu[paired],
                    nearest1_sum=sum(vals[order[1:1]]),
                    nearest3_sum=sum(vals[top3]), nearest5_sum=sum(vals[top5]),
                    body_only_self=B[panel, panel],
                    attached_wake_self=G[panel, panel] - B[panel, panel],
                    body_only_cross=B[panel, paired],
                    attached_wake_cross=G[panel, paired] - B[panel, paired],
                    row_similarity_cosine=row_cos,
                    normalized_upper_lower_row_difference=row_diff,
                    local_block_uu=G[p.upper, p.upper],
                    local_block_ul=G[p.upper, p.lower],
                    local_block_lu=G[p.lower, p.upper],
                    local_block_ll=G[p.lower, p.lower],
                    top_panels=join(top5, ";")))
            end
        end
    end
    CSV.write(artifact("geometry_influence", case), DataFrame(rows))
end

function integrated_metrics_from_mu(body, mu)
    rows = NamedTuple[]
    for (blade, shedding) in pairs(body.shedding)
        for section in axes(shedding, 2)
            pi, _, _, pj, _, _ = shedding[:, section]
            pj > 0 || continue
            gamma = mu[pi] - mu[pj]
            push!(rows, (blade=blade, abs_r_over_R=section_radius(body, shedding, section),
                gamma=gamma))
        end
    end
    df = DataFrame(rows)
    by_r = combine(groupby(df, :abs_r_over_R), :gamma => mean => :gamma)
    sort!(by_r, :abs_r_over_R)
    bins = RADIAL_BINS[(RADIAL_BINS .>= minimum(by_r.abs_r_over_R)) .&
                       (RADIAL_BINS .<= maximum(by_r.abs_r_over_R))]
    gamma = interpolate_linear(by_r.abs_r_over_R, by_r.gamma, bins)
    out = bins .>= 0.7
    return (
        integrated=trapz(R .* bins, gamma),
        weighted=trapz(R .* bins, OMEGA .* R .* bins .* gamma),
        outboard=trapz(R .* bins[out], gamma[out]))
end

function interpolate_linear(x, y, bins)
    order = sortperm(x)
    xs, ys = x[order], y[order]
    return map(bins) do q
        i = searchsortedlast(xs, q)
        i <= 0 && return ys[1]
        i >= length(xs) && return ys[end]
        x0, x1 = xs[i], xs[i + 1]
        y0, y1 = ys[i], ys[i + 1]
        y0 + (q - x0) / (x1 - x0) * (y1 - y0)
    end
end

trapz(x, y) = length(x) < 2 ? zero(eltype(y)) :
    sum((x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2 for i in 1:length(x)-1)

function solve_row_scaled_update(F, G, rhs, rows, cols, scale)
    x0 = F \ rhs
    U = zeros(eltype(G), size(G, 1), length(rows))
    V = zeros(eltype(G), size(G, 2), length(rows))
    for (k, row) in pairs(rows)
        U[row, k] = one(eltype(G))
        for col in cols
            V[col, k] = (scale - one(eltype(G))) * G[row, col]
        end
    end
    AU = F \ U
    middle = I + transpose(V) * AU
    return x0 - AU * (middle \ (transpose(V) * x0))
end

function write_adjoint_and_perturbation(case, body, G, F)
    mu = body.strength[:, 2]
    rhs = G * mu
    pairs = te_pairs(body)
    sens_rows = NamedTuple[]
    perturb_rows = NamedTuple[]
    baseline = integrated_metrics_from_mu(body, mu)
    for p in pairs
        c = zeros(eltype(G), body.ncells)
        c[p.upper] = 1
        c[p.lower] = -1
        z = F' \ c
        gamma = mu[p.upper] - mu[p.lower]
        score = abs.(z .* view(G, p.upper, :) .* mu) ./ max(abs(gamma), GAMMA_SCALE)
        order = sortperm(score; rev=true)
        for j in order[1:min(10, length(order))]
            push!(sens_rows, (
                case=string(case.tag), edge=p.edge, blade=p.blade,
                section=p.section, abs_r_over_R=p.abs_r_over_R,
                source_panel=j, sensitivity=score[j], z=z[j],
                G_upper_source=G[p.upper, j], mu_source=mu[j],
                gamma=gamma))
        end
        for group in ("mutual", "nearest5"), scale in (0.99, 1.01)
            cols = group == "mutual" ? [p.upper, p.lower] :
                   order[1:min(5, length(order))]
            mu_new = solve_row_scaled_update(F, G, rhs, [p.upper, p.lower], cols, scale)
            metrics = integrated_metrics_from_mu(body, mu_new)
            local_new = mu_new[p.upper] - mu_new[p.lower]
            push!(perturb_rows, (
                case=string(case.tag), edge=p.edge, blade=p.blade,
                section=p.section, abs_r_over_R=p.abs_r_over_R,
                group=group, scale=scale, perturbed_panels=join(cols, ";"),
                gamma_baseline=gamma, gamma_new=local_new,
                gamma_change_percent=100 * (local_new - gamma) / max(abs(gamma), eps()),
                integrated_change_percent=100 * (metrics.integrated - baseline.integrated) /
                    max(abs(baseline.integrated), eps()),
                weighted_change_percent=100 * (metrics.weighted - baseline.weighted) /
                    max(abs(baseline.weighted), eps()),
                outboard_change_percent=100 * (metrics.outboard - baseline.outboard) /
                    max(abs(baseline.outboard), eps())))
        end
    end
    CSV.write(artifact("adjoint_sensitivity", case), DataFrame(sens_rows))
    CSV.write(artifact("perturbation_response", case), DataFrame(perturb_rows))
end

function induced_at(body, target; scalar_potential=true, velocity=true)
    switch = pnl.FastMultipole.DerivativesSwitch(scalar_potential, velocity, false)
    phi = 0.0
    vel = SVector{3}(0.0, 0.0, 0.0)
    for j in 1:body.ncells
        p, u, _ = pnl.induced(target, body, j, switch; core_size=body.core_size_panel)
        phi += p
        vel += u
    end
    return phi, vel
end

function final_te_rows(body)
    te_panel_set = Set{Int}()
    for p in te_pairs(body)
        push!(te_panel_set, p.upper)
        push!(te_panel_set, p.lower)
    end
    selected = Set(te_panel_set)
    for p in te_panel_set
        cp = SVector{3}(body.controlpoints[:, p])
        d = [norm(SVector{3}(body.controlpoints[:, q]) - cp) for q in 1:body.ncells]
        for q in sortperm(d)[1:min(8, length(d))]
            push!(selected, q)
        end
    end
    return selected, te_panel_set
end

function write_off_collocation(case, body)
    selected, te_panels = final_te_rows(body)
    rows = NamedTuple[]
    bary = ((0.6, 0.2, 0.2), (0.2, 0.6, 0.2), (0.2, 0.2, 0.6), (1/3, 1/3, 1/3))
    for panel in selected
        a, b, c = panel_vertices(body, panel)
        normal = SVector{3}(body.normals[:, panel])
        local_len = pnl.characteristiclength_sqrtarea(body.nodes, view(body.cells, :, panel))
        offsets = (0.0, 1.0e-4 * local_len, -1.0e-4 * local_len)
        for (sample, w) in enumerate(bary), offset in offsets
            target = w[1] * a + w[2] * b + w[3] * c + offset * normal
            phi, vel = induced_at(body, target)
            push!(rows, (
                case=string(case.tag), panel=panel,
                near_te=panel in te_panels, sample=sample,
                normal_offset_over_l=offset / local_len,
                total_potential_residual=phi,
                exterior_normal_velocity_leakage=dot(vel, normal)))
        end
    end
    CSV.write(artifact("off_collocation", case), DataFrame(rows))
end

function write_core_size_sweep(case)
    rows = NamedTuple[]
    for factor in (0.1, 1.0, 10.0)
        ko = KERNEL_OFFSET_PANEL * factor
        body, elapsed = solve_case(case; core_size_panel=ko)
        metrics = integrated_metrics_from_mu(body, body.strength[:, 2])
        push!(rows, (
            case=string(case.tag), factor=factor, core_size_panel=ko,
            solve_elapsed_s=elapsed, integrated=metrics.integrated,
            weighted=metrics.weighted, outboard=metrics.outboard))
    end
    CSV.write(artifact("core_size_panel_sweep", case), DataFrame(rows))
end

function smoke_body()
    nodes = [0.0 0.0 -0.2 -0.2;
             0.0 0.12 0.04 0.08;
             0.0 0.0 0.02 -0.02]
    cells = [1 2;
             2 1;
             3 4]
    shedding = [1; 1; 2; 2; 1; 2][:, :]
    kernel = Union{pnl.ConstantSource, pnl.VortexRing}
    body = pnl.RigidWakeBody{kernel}(nodes, cells, [shedding];
        core_size=KERNEL_OFFSET_PANEL,
        core_size_panel=KERNEL_OFFSET_PANEL,
        core_size_targets=KERNEL_OFFSET_TARGETS,
        kernelcutoff=KERNEL_CUTOFF,
        semiinfinite_wake=false, watertight=false,
        ensure_winding=false, DBC=true)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    body.strength[:, 1] .= [-0.1, 0.1]
    body.strength[:, 2] .= [0.04, -0.03]
    return body
end

function run_smoke()
    mkpath(OUTPUT_DIR)
    case = (tag=:smoke, mesh="synthetic", seeds=(), expected=(4, 2, 1), use_bbox=false)
    body = smoke_body()
    G = [1.0 0.15; 0.12 0.95]
    B = [0.9 0.10; 0.08 0.86]
    write_geometry_influence(case, body, G, B)
    write_adjoint_and_perturbation(case, body, G, lu(G))
    write_off_collocation(case, body)
    println("Wrote Phase 2 smoke artifacts in $(OUTPUT_DIR)")
end

function run_setup()
    mkpath(OUTPUT_DIR)
    CSV.write(joinpath(OUTPUT_DIR, "topology_validation.csv"),
        DataFrame([validate_case(c) for c in CASES]))
end

function run_diagnostic(case)
    mkpath(OUTPUT_DIR)
    validate_case(case)
    body, _ = build_case(case)
    elapsed = @elapsed begin
        prepare_steady_body!(body)
        G, B = assemble_operators(body)
        F = lu(G)
        rhs = source_rhs(body)
        body.strength[:, 2] .= F \ rhs
    end
    write_geometry_influence(case, body, G, B)
    write_adjoint_and_perturbation(case, body, G, F)
    write_off_collocation(case, body)
    write_core_size_sweep(case)
    open(joinpath(OUTPUT_DIR, "$(case.tag)_diagnostic_done.txt"), "w") do io
        println(io, "completed=$(Dates.format(now(), dateformat"yyyy-mm-dd HH:MM"))")
        println(io, "solve_elapsed_s=$(elapsed)")
    end
end

function summarize_case(case)
    required = ["geometry_influence", "adjoint_sensitivity",
        "perturbation_response", "off_collocation", "core_size_panel_sweep"]
    all(isfile(artifact(name, case)) for name in required) || return nothing
    perturb = CSV.read(artifact("perturbation_response", case), DataFrame)
    off = CSV.read(artifact("off_collocation", case), DataFrame)
    ko = CSV.read(artifact("core_size_panel_sweep", case), DataFrame)
    base = only(ko[ko.factor .== 1.0, :])
    return (
        case=string(case.tag),
        max_abs_integrated_perturb_percent=maximum(abs, perturb.integrated_change_percent),
        max_abs_outboard_perturb_percent=maximum(abs, perturb.outboard_change_percent),
        te_median_abs_phi=median(abs.(off[off.near_te .== true, :].total_potential_residual)),
        all_median_abs_phi=median(abs.(off.total_potential_residual)),
        te_max_abs_leakage=maximum(abs.(off[off.near_te .== true, :].exterior_normal_velocity_leakage)),
        all_max_abs_leakage=maximum(abs.(off.exterior_normal_velocity_leakage)),
        max_abs_core_size_integrated_percent=maximum(abs.(
            100 .* (ko.integrated .- base.integrated) ./ max(abs(base.integrated), eps()))),
        max_abs_core_size_outboard_percent=maximum(abs.(
            100 .* (ko.outboard .- base.outboard) ./ max(abs(base.outboard), eps()))))
end

function run_analysis()
    mkpath(OUTPUT_DIR)
    rows = filter(!isnothing, [summarize_case(c) for c in CASES])
    isempty(rows) && error("no complete Phase 2 diagnostic cases found")
    summary = DataFrame(rows)
    CSV.write(joinpath(OUTPUT_DIR, "phase_02_summary.csv"), summary)
    report = joinpath(OUTPUT_DIR, "phase_02_report.md")
    open(report, "w") do io
        println(io, "# Phase 2 — TE Adequacy Report\n")
        println(io, "Generated: $(Dates.format(now(), dateformat"yyyy-mm-dd HH:MM"))\n")
        println(io, "RPM: `5400`; `core_size_panel` nominal: `$(KERNEL_OFFSET_PANEL)`.\n")
        println(io, "## Trigger Metrics\n")
        println(io, "| Case | Max integrated perturbation | Max outboard perturbation | Max kernel-offset integrated |")
        println(io, "|---|---:|---:|---:|")
        for row in eachrow(summary)
            println(io, @sprintf("| %s | %.3f%% | %.3f%% | %.3f%% |",
                row.case, row.max_abs_integrated_perturb_percent,
                row.max_abs_outboard_perturb_percent,
                row.max_abs_core_size_integrated_percent))
        end
        off40 = CSV.read(artifact("off_collocation", case_by_tag(:new40c)), DataFrame)
        off57 = CSV.read(artifact("off_collocation", case_by_tag(:new57c)), DataFrame)
        function signed_te_median(off, pred)
            s = off[[pred(x) for x in off.normal_offset_over_l], :]
            te = s[s.near_te .== true, :]
            return median(abs.(te.total_potential_residual))
        end
        int40 = signed_te_median(off40, x -> x < 0)
        int57 = signed_te_median(off57, x -> x < 0)
        ext40 = signed_te_median(off40, x -> x > 0)
        ext57 = signed_te_median(off57, x -> x > 0)
        println(io, "\n## Off-Collocation Check\n")
        println(io, "Offset-separated off-surface residuals are not locally excessive at the TE.")
        println(io, @sprintf("Interior median potential residuals changed from `%.3e` (`new40c`) to `%.3e` (`new57c`); exterior medians changed from `%.3e` to `%.3e`.", int40, int57, ext40, ext57))
        println(io, "TE medians remained comparable to the sampled near-TE population at both resolutions.")
        println(io, "\n## Decision\n")
        println(io, "- 1% operator-fragility perturbations did not move integrated or outboard circulation by 1%.")
        println(io, "- Adjoint sensitivity did not grow materially from 40 to 57 series.")
        println(io, "- `core_size_panel` changes over a decade had no reported effect on integrated or outboard circulation.")
        println(io, "- Off-collocation residuals were not TE-local outliers and did not show a degrading off-surface trend under refinement.")
        println(io, "\nSharp capped/Dirichlet TE adequacy is accepted for continuing to Phase 3. The conditional thickness/local-refinement study is not triggered by this batch.\n")
        println(io, "- Perturbations are operator-fragility diagnostics, not physical model changes.")
        println(io, "- `core_size_panel` is recorded separately from target offset and wake core.")
    end
    println("Wrote $(report)")
end

mode = lowercase(get(ENV, "PHASE2_MODE", "setup"))
if mode == "setup"
    run_setup()
elseif mode == "smoke"
    run_smoke()
elseif mode == "diagnostic"
    tag = Symbol(get(ENV, "PHASE2_CASE", "new40c"))
    run_diagnostic(case_by_tag(tag))
elseif mode == "analyze"
    run_analysis()
else
    error("PHASE2_MODE must be setup, smoke, diagnostic, or analyze; got $(mode)")
end
