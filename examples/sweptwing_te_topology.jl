# Sweptwing TE / wake topology + sign diagnostic.
# Builds the same full-span swept wing two ways, solves both, dumps per-edge
# TE tables and consistency checks. Env: SWEPTWING_TE_N_RFL, SWEPTWING_TE_N_SPAN.

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
import LinearAlgebra
import CSV
import DataFrames: DataFrame, nrow
const LA = LinearAlgebra

n_rfl  = parse(Int, get(ENV, "SWEPTWING_TE_N_RFL",  "8"))
n_span = parse(Int, get(ENV, "SWEPTWING_TE_N_SPAN", "40"))
println("Resolution: n_rfl=$n_rfl, n_span=$n_span ($(24*n_rfl*n_span) panels)")

AOA      = 4.2
magVinf  = 30.0
Vinf     = magVinf * [cos(AOA*pi/180), 0, sin(AOA*pi/180)]
rho      = 1.225
b        = 98 * 0.0254
ar       = 5.0
tr       = 1.0
twist_root = 0; twist_tip = 0
lambda   = 45
gamma_di = 0
airfoil  = "airfoil-rae101.csv"
airfoil_path = joinpath(pnl.examples_path, "data")
NDIVS_rfl = [(0.25, n_rfl, 10.0, false), (0.50, n_rfl, 1.0, true), (0.25, n_rfl, 1/10.0, false)]
NDIVS_span = [(1.0, n_span, 1.0, true)]
bodytype = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}
Vinf_hat = Vinf ./ magVinf

# ---- copy of the local builder from sweptwing.jl --------------------------
function simplewing_mirrored_from_negative(b, ar, tr, twist_root, twist_tip, lambda, gamma;
        bodytype, bodyoptargs=(;), airfoil_root, airfoil_tip, airfoil_path,
        rfl_NDIVS, span_NDIVS, delim=",", mirror_tol=100eps(Float64))
    half = simplewing(b, ar, tr, twist_root, twist_tip, lambda, gamma;
        bodytype=bodytype, bodyoptargs=bodyoptargs,
        airfoil_root=airfoil_root, airfoil_tip=airfoil_tip, airfoil_path=airfoil_path,
        rfl_NDIVS=rfl_NDIVS, delim=delim, span_NDIVS=span_NDIVS,
        b_low=-1.0, b_up=0.0)
    half_nodes = half.nodes; half_cells = half.cells
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
    neg_order = sort(collect(axes(half_cells, 2)); by=ci -> half_centers_y[ci])
    pos_order = sort(collect(axes(half_cells, 2)); by=ci -> -half_centers_y[ci])
    cells = Matrix{Int}(undef, 3, 2 * size(half_cells, 2))
    out_ci = 0
    for ci in neg_order
        out_ci += 1; cells[:, out_ci] .= half_cells[:, ci]
    end
    for ci in pos_order
        out_ci += 1; cells[:, out_ci] .= reverse(mirror_index[half_cells[:, ci]])
    end
    te_nodes = Int[]
    for col in eachcol(half.shedding[1])
        pi, nia, nib = col[1], col[2], col[3]
        push!(te_nodes, half_cells[nia, pi]); push!(te_nodes, half_cells[nib, pi])
    end
    full_te_nodes = unique(vcat(te_nodes, mirror_index[te_nodes]))
    sort!(full_te_nodes; by=ni -> nodes[2, ni])
    shedding = pnl.calc_shedding(nodes, cells, full_te_nodes, zeros(eltype(nodes), 3, 0))
    watertight, _ = pnl.iswatertight(nodes, cells)
    final_bodyoptargs = merge((ensure_winding=false,), bodyoptargs)
    return bodytype(nodes, cells, [shedding]; watertight, final_bodyoptargs...)
end

# ---- build both bodies ----------------------------------------------------
println("\nBuilding body_pos (simplewing_mirrored)...")
body_pos = simplewing_mirrored(b, ar, tr, twist_root, twist_tip, lambda, gamma_di;
    bodytype=bodytype, airfoil_root=airfoil, airfoil_tip=airfoil,
    airfoil_path=airfoil_path, rfl_NDIVS=NDIVS_rfl, delim=",", span_NDIVS=NDIVS_span)
println("Building body_neg (simplewing_mirrored_from_negative)...")
body_neg = simplewing_mirrored_from_negative(b, ar, tr, twist_root, twist_tip, lambda, gamma_di;
    bodytype=bodytype, airfoil_root=airfoil, airfoil_tip=airfoil,
    airfoil_path=airfoil_path, rfl_NDIVS=NDIVS_rfl, delim=",", span_NDIVS=NDIVS_span)

wake_direction = reshape(Vinf_hat, :, 1)
for body in (body_pos, body_neg)
    for i in eachindex(body.Das)
        body.Das[i] .= repeat(wake_direction, 1, size(body.Das[i], 2))
    end
end

# ---- solve both for γ -----------------------------------------------------
backend = pnl.FastMultipoleBackend()
Sref = b^2 / ar; c_ref = b / ar
normalization = pnl.WingNormalization(rho, Sref, c_ref)
Dhat = Vinf ./ LA.norm(Vinf); Shat = [0.0, 1.0, 0.0]
Lhat = LA.cross(Dhat, Shat)

function solve!(body, tag)
    println("Solving $tag...")
    frames = pnl.ReferenceFrame(body)
    pmon = pnl.PressureBernoulli(rho)
    fmon = pnl.ForceMonitor(1, 1; i_frame=-1, normalization=normalization,
        correct_kuttacondition=false, verbose=false)
    @time pnl.steady!(body, frames, Vinf;
        body_solvers=pnl.Backslash(body), backend=backend,
        monitors=(pmon, fmon), path=nothing, verbose=false)
    CL = LA.dot(fmon.force[:, 1], Lhat)
    println("  $tag CL = $(round(CL, digits=4))")
    return CL
end

CL_pos = solve!(body_pos, "body_pos")
CL_neg = solve!(body_neg, "body_neg")

# ---- TE table builder -----------------------------------------------------
function te_table(body)
    sh = body.shedding[1]
    nedges = size(sh, 2)
    Gi = pnl.get_Gammai(body)
    gamma = view(body.strength, :, Gi)
    rows = NamedTuple[]
    for i in 1:nedges
        pi_up = sh[1, i]; nia = sh[2, i]; nib = sh[3, i]
        pj_lo = sh[4, i]; nic = sh[5, i]; nid = sh[6, i]
        node_a = body.cells[nia, pi_up]
        node_b = body.cells[nib, pi_up]
        pos_a = body.nodes[:, node_a]
        pos_b = body.nodes[:, node_b]
        y_mid = 0.5*(pos_a[2] + pos_b[2])
        n_up = body.normals[:, pi_up]
        n_lo = pj_lo > 0 ? body.normals[:, pj_lo] : fill(NaN, 3)
        # Das columns (use shedding_full upper-row mapping)
        Da_col_a = body.shedding_full[5, pi_up]   # for vertex nia (= node_a)
        Da_col_b = body.shedding_full[6, pi_up]   # for vertex nib (= node_b)
        Da_a = body.Das[1][:, Da_col_a]
        Da_b = body.Das[1][:, Da_col_b]
        edge = pos_b .- pos_a
        n_wake_raw = LA.cross(edge, Da_a)
        n_wake = n_wake_raw ./ max(LA.norm(n_wake_raw), eps())
        g_up = gamma[pi_up]
        g_lo = pj_lo > 0 ? gamma[pj_lo] : NaN
        push!(rows, (
            i=i, pi_up=pi_up, pj_lo=pj_lo, y_mid=y_mid,
            nia=nia, nib=nib, nic=nic, nid=nid,
            node_a=node_a, node_b=node_b,
            xa=pos_a[1], ya=pos_a[2], za=pos_a[3],
            xb=pos_b[1], yb=pos_b[2], zb=pos_b[3],
            n_up_x=n_up[1], n_up_y=n_up[2], n_up_z=n_up[3],
            n_lo_x=n_lo[1], n_lo_y=n_lo[2], n_lo_z=n_lo[3],
            Da_col_a=Da_col_a, Da_col_b=Da_col_b,
            Da_a_x=Da_a[1], Da_a_y=Da_a[2], Da_a_z=Da_a[3],
            Da_b_x=Da_b[1], Da_b_y=Da_b[2], Da_b_z=Da_b[3],
            Da_a_dot_Vinf = LA.dot(Da_a, Vinf_hat),
            Da_b_dot_Vinf = LA.dot(Da_b, Vinf_hat),
            n_wake_dot_jump = LA.dot(n_wake, n_up .- n_lo),
            gamma_up=g_up, gamma_lo=g_lo, dGamma=g_up - g_lo,
        ))
    end
    return DataFrame(rows)
end

println("\n#===== building TE tables =====#")
te_pos = te_table(body_pos)
te_neg = te_table(body_neg)
out_dir = joinpath(@__DIR__, "..", "data", "sweptwing_te_diagnostics")
mkpath(out_dir)
CSV.write(joinpath(out_dir, "te_pos_$(n_rfl)x$(n_span).csv"), te_pos)
CSV.write(joinpath(out_dir, "te_neg_$(n_rfl)x$(n_span).csv"), te_neg)
println("Wrote CSVs to $out_dir")

# ---- per-body consistency checks ------------------------------------------
function report_body(label, body, te)
    println("\n#===== $label =====#")
    wt, _ = pnl.iswatertight(body.nodes, body.cells)
    println("  watertight: $wt;  ncells = $(body.ncells)")
    n_pos_y = count(c -> body.controlpoints[2, c] > 0, 1:body.ncells)
    n_neg_y = count(c -> body.controlpoints[2, c] < 0, 1:body.ncells)
    println("  panels with cp_y>0 / cp_y<0 : $n_pos_y / $n_neg_y")

    # winding consistency via dot(normal, controlpoint - centroid)
    centroid = vec(sum(body.controlpoints, dims=2) ./ body.ncells)
    outward = 0; inward = 0
    for p in 1:body.ncells
        r = body.controlpoints[:, p] .- centroid
        d = LA.dot(body.normals[:, p], r)
        d > 0 ? (outward += 1) : (inward += 1)
    end
    println("  normals outward / inward (centroid test): $outward / $inward")

    # shedding-edge sign checks
    n_up_neg_z = count(<(0), te.n_up_z)
    n_lo_pos_z = count(>(0), te.n_lo_z)
    n_Da_back  = count(<(0), te.Da_a_dot_Vinf) + count(<(0), te.Da_b_dot_Vinf)
    println("  shedding edges: $(nrow(te))")
    println("  edges with n_up·ẑ < 0 (upper panel pointing down): $n_up_neg_z")
    println("  edges with n_lo·ẑ > 0 (lower panel pointing up):   $n_lo_pos_z")
    println("  Das columns with Da·V̂∞ < 0 (wake convects upstream): $n_Da_back")
    n_wake_pos = count(>(0), te.n_wake_dot_jump)
    n_wake_neg = count(<(0), te.n_wake_dot_jump)
    println("  edges with n_wake·(n_up-n_lo) > 0 / < 0 : $n_wake_pos / $n_wake_neg")

    # dGamma mirror symmetry
    # Find pairs by y_mid sign
    pos_mask = te.y_mid .> 0
    neg_mask = te.y_mid .< 0
    pos_rows = te[pos_mask, :]
    neg_rows = te[neg_mask, :]
    println("  TE edges +y / −y : $(nrow(pos_rows)) / $(nrow(neg_rows))")
    if nrow(pos_rows) == nrow(neg_rows) > 0
        # match by sorted y magnitude
        psort = sortperm(pos_rows.y_mid)
        nsort = sortperm(neg_rows.y_mid; rev=true)
        y_max_diff = maximum(abs.(pos_rows.y_mid[psort] .+ neg_rows.y_mid[nsort]))
        dG_diff = pos_rows.dGamma[psort] .- neg_rows.dGamma[nsort]
        rel = abs.(dG_diff) ./ max.(abs.(pos_rows.dGamma[psort]), 1e-12)
        println("  max |Δy_mid mirror|: $(round(y_max_diff, sigdigits=4))")
        println("  max |dGamma+ − dGamma−|: $(round(maximum(abs, dG_diff), sigdigits=4))")
        println("  max relative dGamma asym:  $(round(maximum(rel), sigdigits=4))")
    end

    # γ mirror symmetry on all panels
    cps = body.controlpoints
    Gi = pnl.get_Gammai(body)
    gamma = view(body.strength, :, Gi)
    flipped = copy(cps); flipped[2, :] .*= -1
    matched = zeros(Int, body.ncells)
    matched_dist = zeros(body.ncells)
    for p in 1:body.ncells
        # find nearest controlpoint to flipped[:, p]
        best_d = Inf; best_q = 0
        for q in 1:body.ncells
            d = (cps[1,q]-flipped[1,p])^2 + (cps[2,q]-flipped[2,p])^2 + (cps[3,q]-flipped[3,p])^2
            if d < best_d; best_d = d; best_q = q; end
        end
        matched[p] = best_q
        matched_dist[p] = sqrt(best_d)
    end
    g_asym = abs.(gamma .- gamma[matched])
    g_rel  = g_asym ./ max.(abs.(gamma), 1e-12)
    println("  max controlpoint mirror-match distance: $(round(maximum(matched_dist), sigdigits=4))")
    println("  max |γ_p − γ_p'|:        $(round(maximum(g_asym), sigdigits=4))")
    println("  max relative γ asym:      $(round(maximum(g_rel),  sigdigits=4))")
    return (; matched, gamma=collect(gamma))
end

info_pos = report_body("body_pos", body_pos, te_pos)
info_neg = report_body("body_neg", body_neg, te_neg)

# ---- cross-builder check on +y half --------------------------------------
println("\n#===== cross-builder check on +y half =====#")
function match_te_by_position(te_a, te_b)
    # For each row in te_a with y_mid>0, find row in te_b with closest y_mid (same sign).
    a_pos = filter(r -> r.y_mid > 0, eachrow(te_a)) |> collect
    b_pos = filter(r -> r.y_mid > 0, eachrow(te_b)) |> collect
    matches = []
    for ra in a_pos
        best = (Inf, nothing)
        for rb in b_pos
            d = abs(ra.y_mid - rb.y_mid) + 0.1*sqrt((ra.xa-rb.xa)^2 + (ra.za-rb.za)^2)
            d < best[1] && (best = (d, rb))
        end
        push!(matches, (a=ra, b=best[2], dist=best[1]))
    end
    return matches
end
matches = match_te_by_position(te_pos, te_neg)
dG_diff = [abs(m.a.dGamma - m.b.dGamma) for m in matches]
nup_z_diff = [abs(m.a.n_up_z - m.b.n_up_z) for m in matches]
Da_z_diff = [abs(m.a.Da_a_dot_Vinf - m.b.Da_a_dot_Vinf) for m in matches]
n_wake_signs_a = [sign(m.a.n_wake_dot_jump) for m in matches]
n_wake_signs_b = [sign(m.b.n_wake_dot_jump) for m in matches]
n_wake_disagree = count(((sa, sb),) -> sa != sb, zip(n_wake_signs_a, n_wake_signs_b))

println("matched $(length(matches)) TE edges on +y half")
println("  max |dGamma_pos_builder − dGamma_neg_builder|: $(round(maximum(dG_diff), sigdigits=4))")
println("  max |n_up·ẑ diff|:                $(round(maximum(nup_z_diff), sigdigits=4))")
println("  max |Da_a·V̂∞ diff|:               $(round(maximum(Da_z_diff), sigdigits=4))")
println("  edges where sign(n_wake·jump) differs: $n_wake_disagree / $(length(matches))")

println("\n#===== integrated CLs =====#")
println("CL body_pos = $(round(CL_pos, digits=4)),  CL body_neg = $(round(CL_neg, digits=4))")
println("gap = $(round(abs(CL_pos - CL_neg), digits=4))")

# ---- cell-winding cross-builder comparison --------------------------------
println("\n#===== cell-winding cross-builder check =====#")
function panel_signature(body, p)
    # invariant of triangle: sorted node positions (round to defeat tiny float diff)
    a = body.nodes[:, body.cells[1, p]]
    b = body.nodes[:, body.cells[2, p]]
    c = body.nodes[:, body.cells[3, p]]
    pts = sort([round.(a; digits=6), round.(b; digits=6), round.(c; digits=6)])
    return (pts[1], pts[2], pts[3])
end
function compare_meshes(body_pos, body_neg)
    println("  body_pos node bbox y: [$(round(minimum(body_pos.nodes[2,:]),digits=4)), $(round(maximum(body_pos.nodes[2,:]),digits=4))]")
    println("  body_neg node bbox y: [$(round(minimum(body_neg.nodes[2,:]),digits=4)), $(round(maximum(body_neg.nodes[2,:]),digits=4))]")
    println("  body_pos node count: $(size(body_pos.nodes,2)); body_neg: $(size(body_neg.nodes,2))")
    # Match each body_neg panel center to nearest body_pos panel center
    cp_pos = body_pos.controlpoints; cp_neg = body_neg.controlpoints
    Gi_pos = pnl.get_Gammai(body_pos); Gi_neg = pnl.get_Gammai(body_neg)
    g_pos = view(body_pos.strength, :, Gi_pos)
    g_neg = view(body_neg.strength, :, Gi_neg)
    max_d = 0.0; matched_close = 0; dG_max = 0.0; dG_argmax_y = 0.0
    for q in 1:body_neg.ncells
        # brute-force nearest neighbour in pos
        best_d = Inf; best_p = 0
        for p in 1:body_pos.ncells
            d = (cp_pos[1,p]-cp_neg[1,q])^2 + (cp_pos[2,p]-cp_neg[2,q])^2 + (cp_pos[3,p]-cp_neg[3,q])^2
            if d < best_d; best_d = d; best_p = p; end
        end
        d = sqrt(best_d)
        max_d = max(max_d, d)
        d < 1e-3 && (matched_close += 1)
        if abs(g_pos[best_p] - g_neg[q]) > dG_max
            dG_max = abs(g_pos[best_p] - g_neg[q])
            dG_argmax_y = cp_neg[2, q]
        end
    end
    println("  max nearest-controlpoint dist (neg→pos): $(round(max_d, sigdigits=4))")
    println("  panels with nearest dist < 1e-3 (matched): $matched_close / $(body_neg.ncells)")
    println("  max |γ_pos − γ_neg| at nearest-panel match: $(round(dG_max, sigdigits=4)) at y=$(round(dG_argmax_y,digits=4))")

    # Node-set comparison: for each node in body_neg, find nearest node in body_pos
    max_node_d = 0.0; n_node_close = 0
    for q in 1:size(body_neg.nodes, 2)
        best_d = Inf
        for p in 1:size(body_pos.nodes, 2)
            d = (body_pos.nodes[1,p]-body_neg.nodes[1,q])^2 +
                (body_pos.nodes[2,p]-body_neg.nodes[2,q])^2 +
                (body_pos.nodes[3,p]-body_neg.nodes[3,q])^2
            d < best_d && (best_d = d)
        end
        d = sqrt(best_d)
        max_node_d = max(max_node_d, d)
        d < 1e-8 && (n_node_close += 1)
    end
    println("  max nearest-node dist (neg→pos): $(round(max_node_d, sigdigits=4))")
    println("  nodes matched within 1e-8: $n_node_close / $(size(body_neg.nodes, 2))")
end
compare_meshes(body_pos, body_neg)

function compare_winding(body_pos, body_neg)
    sig2p_pos = Dict{Any, Int}()
    for p in 1:body_pos.ncells
        sig2p_pos[panel_signature(body_pos, p)] = p
    end
    flipped = 0; matched_pairs = 0; missing_neg = 0
    Gi_pos = pnl.get_Gammai(body_pos); Gi_neg = pnl.get_Gammai(body_neg)
    g_pos = view(body_pos.strength, :, Gi_pos)
    g_neg = view(body_neg.strength, :, Gi_neg)
    dgamma_at_same_panel = Float64[]
    pos_normal_sign_flips = 0
    for q in 1:body_neg.ncells
        sig = panel_signature(body_neg, q)
        p = get(sig2p_pos, sig, 0)
        if p == 0
            missing_neg += 1; continue
        end
        matched_pairs += 1
        cp = body_pos.cells[:, p]; cq = body_neg.cells[:, q]
        cyclic(a, b) = any(circshift(a, k) == b for k in 0:2)
        if cyclic(cp, cq)
            # same winding
        elseif cyclic(cp, reverse(cq))
            flipped += 1
        end
        # also compare normals directly (independent of node-id renumbering)
        if LA.dot(body_pos.normals[:, p], body_neg.normals[:, q]) < 0
            pos_normal_sign_flips += 1
        end
        push!(dgamma_at_same_panel, g_pos[p] - g_neg[q])
    end
    println("  panels matched by sig: $matched_pairs / $(body_neg.ncells), unmatched: $missing_neg")
    println("  panels with reversed cell winding (sig matched): $flipped")
    println("  panels with normal·normal < 0 between builders:  $pos_normal_sign_flips")
    if !isempty(dgamma_at_same_panel)
        println("  max |γ_pos − γ_neg| at same physical panel: $(round(maximum(abs, dgamma_at_same_panel), sigdigits=4))")
    end
end
compare_winding(body_pos, body_neg)

# ---- non-watertight diagnosis ---------------------------------------------
println("\n#===== watertight diagnosis =====#")
function edge_count(body)
    counts = Dict{Tuple{Int,Int}, Int}()
    for p in 1:body.ncells
        v = body.cells[:, p]
        for (a, b) in ((v[1], v[2]), (v[2], v[3]), (v[3], v[1]))
            key = a < b ? (a, b) : (b, a)
            counts[key] = get(counts, key, 0) + 1
        end
    end
    return counts
end
for (label, body) in (("body_pos", body_pos), ("body_neg", body_neg))
    ec = edge_count(body)
    n1 = count(v -> v == 1, values(ec))
    n2 = count(v -> v == 2, values(ec))
    n3p = count(v -> v >= 3, values(ec))
    println("  $label edge sharing: 1-shared=$n1, 2-shared=$n2, 3+-shared=$n3p (total edges=$(length(ec)))")
end

# ---- dGamma vs y plot ----------------------------------------------------
using GeometricTools: PyPlot as plt
fig = plt.figure(figsize=(7, 4))
ax = fig.subplots(1, 1)
ax.plot(te_pos.y_mid .* 2 ./ b, te_pos.dGamma, "o-", label="body_pos (CL=$(round(CL_pos,digits=3)))", linewidth=1.5)
ax.plot(te_neg.y_mid .* 2 ./ b, te_neg.dGamma, "s--", label="body_neg (CL=$(round(CL_neg,digits=3)))", linewidth=1.5)
ax.set_xlabel("2y/b")
ax.set_ylabel("TE jump dGamma = γ_up − γ_lo")
ax.set_title("TE circulation jump, both mirrored discretizations (n_rfl=$n_rfl, n_span=$n_span)")
ax.legend(frameon=false)
ax.grid(true, alpha=0.3)
fig.tight_layout()
fig.savefig(joinpath(@__DIR__, "..", "sweptwing_dGamma_te.png"), dpi=150)
println("Saved TE dGamma plot to sweptwing_dGamma_te.png")

# ---- diagonal-choice verification ----------------------------------------
# Rebuild a body that uses body_neg's nodes but body_pos's CELL TRIANGULATION,
# via a node-index permutation. If the CL matches body_pos, the discrepancy is
# entirely from the quad-diagonal choice.
println("\n#===== diagonal-choice verification =====#")
function build_node_perm(nodes_from, nodes_to; tol=1e-10)
    n = size(nodes_from, 2)
    @assert size(nodes_to, 2) == n
    perm = zeros(Int, n)
    for i in 1:n
        for j in 1:n
            d = sqrt((nodes_to[1, j] - nodes_from[1, i])^2 +
                     (nodes_to[2, j] - nodes_from[2, i])^2 +
                     (nodes_to[3, j] - nodes_from[3, i])^2)
            if d < tol
                perm[i] = j; break
            end
        end
        perm[i] == 0 && error("no match for node $i")
    end
    return perm
end

perm = build_node_perm(body_pos.nodes, body_neg.nodes)
println("Built node permutation pos→neg (covers all $(length(perm)) nodes)")

# Translate body_pos's cells into body_neg's node indexing
cells_hybrid = similar(body_pos.cells)
for p in axes(body_pos.cells, 2)
    for k in 1:3
        cells_hybrid[k, p] = perm[body_pos.cells[k, p]]
    end
end

# Translate body_pos's TE node list into body_neg's indexing, recompute shedding
te_pos_nodes = Int[]
for col in eachcol(body_pos.shedding[1])
    pi_up, nia, nib = col[1], col[2], col[3]
    push!(te_pos_nodes, body_pos.cells[nia, pi_up])
    push!(te_pos_nodes, body_pos.cells[nib, pi_up])
end
te_hybrid_nodes = unique([perm[n] for n in te_pos_nodes])
sort!(te_hybrid_nodes; by=ni -> body_neg.nodes[2, ni])
shedding_hybrid = pnl.calc_shedding(body_neg.nodes, cells_hybrid, te_hybrid_nodes,
                                    zeros(eltype(body_neg.nodes), 3, 0))

watertight_h, _ = pnl.iswatertight(body_neg.nodes, cells_hybrid)
println("Hybrid body (pos cells on neg nodes): watertight=$watertight_h")

body_hybrid = bodytype(body_neg.nodes, cells_hybrid, [shedding_hybrid];
                       watertight=watertight_h, ensure_winding=false)
for i in eachindex(body_hybrid.Das)
    body_hybrid.Das[i] .= repeat(wake_direction, 1, size(body_hybrid.Das[i], 2))
end

CL_hybrid = solve!(body_hybrid, "body_hybrid (pos diagonals)")

println("\n  CL body_pos    = $(round(CL_pos, digits=4))")
println("  CL body_neg    = $(round(CL_neg, digits=4))")
println("  CL body_hybrid = $(round(CL_hybrid, digits=4))  (= body_pos if hypothesis is right)")
println("  body_pos - body_hybrid = $(round(CL_pos - CL_hybrid, sigdigits=4))")
println("  body_neg - body_hybrid = $(round(CL_neg - CL_hybrid, sigdigits=4))")

# ---- symmetric hybrid: body_pos nodes + body_neg cells -------------------
println("\n#===== symmetric hybrid (body_pos nodes + body_neg cells) =====#")
perm_n2p = build_node_perm(body_neg.nodes, body_pos.nodes)
cells_hybrid_neg = similar(body_neg.cells)
for p in axes(body_neg.cells, 2)
    for k in 1:3
        cells_hybrid_neg[k, p] = perm_n2p[body_neg.cells[k, p]]
    end
end
te_neg_nodes = Int[]
for col in eachcol(body_neg.shedding[1])
    pi_up, nia, nib = col[1], col[2], col[3]
    push!(te_neg_nodes, body_neg.cells[nia, pi_up])
    push!(te_neg_nodes, body_neg.cells[nib, pi_up])
end
te_hybrid_neg_nodes = unique([perm_n2p[n] for n in te_neg_nodes])
sort!(te_hybrid_neg_nodes; by=ni -> body_pos.nodes[2, ni])
shedding_hybrid_neg = pnl.calc_shedding(body_pos.nodes, cells_hybrid_neg, te_hybrid_neg_nodes,
                                        zeros(eltype(body_pos.nodes), 3, 0))
watertight_hn, _ = pnl.iswatertight(body_pos.nodes, cells_hybrid_neg)
body_hybrid_neg = bodytype(body_pos.nodes, cells_hybrid_neg, [shedding_hybrid_neg];
                           watertight=watertight_hn, ensure_winding=false)
for i in eachindex(body_hybrid_neg.Das)
    body_hybrid_neg.Das[i] .= repeat(wake_direction, 1, size(body_hybrid_neg.Das[i], 2))
end
CL_hybrid_neg = solve!(body_hybrid_neg, "body_hybrid_neg (neg diagonals)")

println("\n  CL body_neg        = $(round(CL_neg, digits=4))")
println("  CL body_hybrid_neg = $(round(CL_hybrid_neg, digits=4))  (expect ≈ body_neg)")
println("  body_neg - body_hybrid_neg = $(round(CL_neg - CL_hybrid_neg, sigdigits=4))")
