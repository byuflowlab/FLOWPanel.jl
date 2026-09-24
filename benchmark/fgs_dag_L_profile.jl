#=##############################################################################
BRAINSTORM 021 L-shortening gate-0: weighted critical-path analysis of the
:dagteam lower DAG (pure graph analysis — NO solves, NO timing).

Bounds the prize of edge-level partial pulls (slate item #1,
fgs_lshortening_reset_prompt_20260924.md) before any executor work:
predicted sweep speedup at the L-bound = L_node / L_edge.

Model (byte-weighted, the campaign's established currency — Stage 1/2 showed
the sweep is bandwidth/critical-path bound, and `prio` already uses
nof(i)*ptot[i]):
  edge GEMV  e_ij   = szTM * nof(i) * nof(j)      (stream L_ij block)
  node GEMV  g_i    = szTM * nof(i) * ptot[i]      (aggregated pull)
  leaf LU    lu_i   = szTS * nof(i)^2              (two triangular solves)
  reduction  red_i  = szTS * nof(i) * |preds(i)|   (sum partials, edge mode)

Schedules (infinite processors — pure L):
  node   : finish(i) = max_j finish(j)                        + g_i + lu_i
  edge   : finish(i) = max_j (finish(j) + e_ij)               + red_i + lu_i
  cutoff : preds split at e_ij >= THETA; big edges overlap as partials,
           small edges stay one aggregated GEMV after all small preds:
           finish(i) = max( max_big (finish_j + e_ij),
                            max_small finish_j + Σ_small e_ij ) + red_i + lu_i
  lufloor: finish(i) = max_j finish(j) + lu_i   (no pull split can beat the
           LU dependency chain; reported as the asymptote)

Usage (local, ≤4 threads per house rules; SKIP_B=1 makes the fixture
geometry-only — this script never solves):
  RUNG=R4 SKIP_B=1 THREADING_MODE=multi EXPECT_JULIA_THREADS=4 \
    julia --project=benchmark -t 4 benchmark/fgs_dag_L_profile.jl [outdir]
=###############################################################################

import TOML

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

outdir = isempty(ARGS) ? pwd() : ARGS[1]
mkpath(outdir)

# Champion R4 knobs (retained_r4_champion.toml); tolerance/max_iterations are
# irrelevant here (no solve) but kept for constructor parity.
champ = TOML.parsefile(joinpath(@__DIR__, "retained_r4_champion.toml"))
@assert champ["rung"] == rung "champion TOML rung $(champ["rung"]) != RUNG $rung"

println("Constructing dagteam solver (P=$(champ["P"]), MAC=$(champ["MAC"]), " *
        "leaf=$(champ["leaf"]), precision=$(champ["dagteam_precision"]))...")
t0 = time_ns()
solver = pnl.FGSSolver(rotor; expansion_order=champ["P"],
    multipole_acceptance=champ["MAC"], leaf_size=champ["leaf"],
    inner_iterations=champ["inner"], max_iterations=champ["max_iterations"],
    tolerance=champ["tolerance"], rlx=champ["rlx"], shrink=true, recenter=false,
    reverse_pass=false, cache_leaf_lu=true, sweep_order=:dagteam,
    dagteam_precision=Symbol(champ["dagteam_precision"]),
    verbose=false, project_solution=false, solution_history_length=0)
println("  solver constructed in $(round((time_ns()-t0)/1e9, digits=1)) s")

plan = solver.fgs.dagteam
plan === nothing && error("dagteam plan missing")

n = length(plan.preds)
nof(i) = plan.offset[i+1] - plan.offset[i]
szTM = sizeof(eltype(eltype(plan.Lmat)))   # coefficient bytes
szTS = sizeof(eltype(plan.x))              # sweep-state bytes

g   = [szTM * nof(i) * plan.ptot[i] for i in 1:n]
lu  = [szTS * nof(i)^2 for i in 1:n]
red = [szTS * nof(i) * length(plan.preds[i]) for i in 1:n]
edge(i, j) = szTM * nof(i) * nof(j)

# ---- totals ----
W_node = sum(g) + sum(lu)
W_edge = W_node + sum(red)
W_upper = sum(szTM * plan.mup[j] * nof(j) for j in 1:n)   # filler work (context)
n_edges = sum(length(plan.preds[i]) for i in 1:n)

# ---- node-level L (with backtrace for path composition) ----
fin_n = zeros(Float64, n); pred_of = zeros(Int, n)
for i in 1:n
    a, ai = 0.0, 0
    for j in plan.preds[i]
        fin_n[j] > a && (a = fin_n[j]; ai = j)
    end
    fin_n[i] = a + g[i] + lu[i]; pred_of[i] = ai
end
L_node, iend = findmax(fin_n)
path = Int[]; let i = iend
    while i != 0; push!(path, i); i = pred_of[i]; end
end
reverse!(path)
path_g, path_lu = sum(g[i] for i in path), sum(lu[i] for i in path)

# ---- LU-chain floor ----
fin_f = zeros(Float64, n)
for i in 1:n
    a = 0.0
    for j in plan.preds[i]; fin_f[j] > a && (a = fin_f[j]); end
    fin_f[i] = a + lu[i]
end
L_lufloor = maximum(fin_f)

# ---- edge-level L, with cutoff sweep (THETA=0 is the pure edge variant) ----
KB = 1024
thetas = [0, 4KB, 16KB, 64KB, 256KB, 1024KB]
L_theta = Float64[]
for th in thetas
    fin = zeros(Float64, n)
    for i in 1:n
        t_big, t_small, s_small = 0.0, 0.0, 0.0
        for j in plan.preds[i]
            e = Float64(edge(i, j))
            if e >= th
                v = fin[j] + e; v > t_big && (t_big = v)
            else
                fin[j] > t_small && (t_small = fin[j]); s_small += e
            end
        end
        fin[i] = max(t_big, t_small + s_small) + red[i] + lu[i]
    end
    push!(L_theta, maximum(fin))
end
L_edge = L_theta[1]

# ---- selective split: edge-split only the top-X% leaves by prio ----
# Design trade for the executor: full edge split multiplies tasks/sweep ~46x
# (per-task queue overhead is then first-order at ~16KB/edge); splitting only
# high-prio (critical-path-heavy) leaves keeps the task count low. Reports L
# and tasks/sweep per split fraction. Selection by plan.prio (byte-weighted
# downstream path), the same key the executor already has.
using Statistics: quantile
fracs = [0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.0]
sel_rows = Tuple{Float64,Float64,Int}[]
for f in fracs
    th_p = quantile(plan.prio, 1 - f)
    S = [plan.prio[i] >= th_p for i in 1:n]
    fin = zeros(Float64, n)
    for i in 1:n
        if S[i]
            t = 0.0
            for j in plan.preds[i]
                v = fin[j] + edge(i, j); v > t && (t = v)
            end
            fin[i] = t + red[i] + lu[i]
        else
            a = 0.0
            for j in plan.preds[i]; fin[j] > a && (a = fin[j]); end
            fin[i] = a + g[i] + lu[i]
        end
    end
    ntasks = n + sum(length(plan.preds[i]) for i in 1:n if S[i]; init=0)
    push!(sel_rows, (f, maximum(fin), ntasks))
end

# ---- unit-depth profile (node DAG) ----
lev = zeros(Int, n)
for i in 1:n
    m = 0
    for j in plan.preds[i]; lev[j] > m && (m = lev[j]); end
    lev[i] = m + 1
end
depth = maximum(lev)
width = [count(==(k), lev) for k in 1:depth]

# ---- report ----
mb(x) = round(x / 1e6, digits=3)
println("\n## DAG structure (R4 champion plan)")
println("leaves=$n  edges=$n_edges  avg preds=$(round(n_edges/n, digits=1))  " *
        "unit depth=$depth  mean/max level width=$(round(sum(width)/depth, digits=1))/$(maximum(width))")
println("nof: min/med/max = $(minimum(nof.(1:n)))/$(round(Int, sum(nof.(1:n))/n))/$(maximum(nof.(1:n)))  " *
        "szTM=$szTM szTS=$szTS")
println("W_lower+LU=$(mb(W_node)) MB/sweep  (+red: $(mb(W_edge)))  W_upper(filler)=$(mb(W_upper)) MB")

println("\n## Critical paths (MB streamed)")
rows = [("node (today)", L_node), ("edge (theta=0)", L_edge)]
append!(rows, [("edge theta=$(div(t,KB))KB", L) for (t, L) in zip(thetas[2:end], L_theta[2:end])])
push!(rows, ("LU-chain floor", L_lufloor))
println("| variant | L (MB) | W/L | L_node/L (speedup bound) |")
println("|---|---|---|---|")
for (name, L) in rows
    println("| $name | $(mb(L)) | $(round(W_node/L, digits=1)) | $(round(L_node/L, digits=2)) |")
end
println("\nnode critical path: $(length(path)) leaves; GEMV $(mb(path_g)) MB vs LU $(mb(path_lu)) MB")

println("\n## Selective split (top-X% prio leaves edge-split; node tasks otherwise)")
println("| split frac | L (MB) | L_node/L | tasks/sweep |")
println("|---|---|---|---|")
for (f, L, nt) in sel_rows
    println("| $(round(Int, 100f))% | $(mb(L)) | $(round(L_node/L, digits=2)) | $nt |")
end

println("\n## Plan-C bounds T_j = max(W/j, L), normalized to node@j64")
Tn(j) = max(W_node/j, L_node)
Te(j) = max(W_edge/j, L_edge)
println("| j | node | edge | predicted speedup |")
println("|---|---|---|---|")
for j in (16, 32, 64, 128)
    println("| $j | $(mb(Tn(j))) | $(mb(Te(j))) | $(round(Tn(j)/Te(j), digits=2)) |")
end

# ---- CSV ----
open(joinpath(outdir, "fgs_dag_L_profile_$(rung).csv"), "w") do io
    println(io, "variant,theta_bytes,L_bytes,W_bytes")
    println(io, "node,,$L_node,$W_node")
    for (t, L) in zip(thetas, L_theta)
        println(io, "edge,$t,$L,$W_edge")
    end
    println(io, "lufloor,,$L_lufloor,$(sum(lu))")
end
open(joinpath(outdir, "fgs_dag_L_widths_$(rung).csv"), "w") do io
    println(io, "level,width")
    for k in 1:depth; println(io, "$k,$(width[k])"); end
end
open(joinpath(outdir, "fgs_dag_L_selective_$(rung).csv"), "w") do io
    println(io, "split_frac,L_bytes,tasks_per_sweep")
    for (f, L, nt) in sel_rows; println(io, "$f,$L,$nt"); end
end
# graph dump for offline iteration (no fixture rebuild): leaf sizes + edges
open(joinpath(outdir, "fgs_dag_L_graph_$(rung).csv"), "w") do io
    println(io, "leaf,nof,ptot,prio")
    for i in 1:n; println(io, "$i,$(nof(i)),$(plan.ptot[i]),$(plan.prio[i])"); end
end
open(joinpath(outdir, "fgs_dag_L_edges_$(rung).csv"), "w") do io
    println(io, "target,source")
    for i in 1:n, j in plan.preds[i]; println(io, "$i,$j"); end
end
println("\nCSVs written to $outdir")
