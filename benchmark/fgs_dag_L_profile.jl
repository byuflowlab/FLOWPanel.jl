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
  RUNG=R4 SKIP_B=1 THREADING_MODE=multi EXPECT_JULIA_THREADS=4 BENCH_BLAS_THREADS=1 \
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

# ---- 033 A-T2: cooperative split-graph critical path (COOP_SPLIT=1) ----
# Each leaf's lower-GEMV cost g_i is divided among w cooperative subtasks (LU
# bytes NOT divided); with infinite processors the leaf's elapsed pull becomes
# g_i/w plus elapsed overhead c (swept, expressed in µs-equivalent bytes).
# Equal parallel subtasks each incur c: elapsed overhead is c, summed worker
# overhead is w*c (charged separately in Ww below). Thus c maps to A-T3's
# whole-team elapsed h_w, not to an assumed serialized per-worker cost.
# Split policies:
#   all       : every leaf with g_i > 0 splits
#   selective : split only when it shortens elapsed pull (g_i/w + c < g_i,
#               i.e. g_i > c/(1-1/w) — the measured-cost rule, A-T3)
# Byte<->time conversion from gate-0's own calibration: L_node = 289.5 MB/sweep
# against the measured ~2.0 s sweep floor over 81 sweeps at R4 j64
# => ~11.7 GB/s critical-path streaming rate, so 1 µs ≈ 11.72 KB.
if get(ENV, "COOP_SPLIT", "0") == "1"
    BYTES_PER_US = 289.5e6 * 81 / 2.0 * 1e-6      # bytes per µs ≈ 11.72e3
    dL = abs(L_node - 289.462e6) / 289.462e6
    println("\n## 033 A-T2 cooperative split graph")
    println("unsplit L_node regression vs gate-0 289.462 MB: " *
            "$(mb(L_node)) MB (rel. diff $(round(100dL, digits=3))%)")
    dL < 0.01 || error("gate-0 unsplit L_node not reproduced within 1% — stop per reset prompt")
    println("byte<->time: 1 µs = $(round(BYTES_PER_US/1e3, digits=2)) KB " *
            "(gate-0: 289.5 MB x 81 sweeps / 2.0 s)")
    ovh_grid_us = [0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0]
    coop_rows = NamedTuple[]
    gel = zeros(Float64, n)
    for w in (2, 4), policy in (:all, :selective), c_us in ovh_grid_us
        c = c_us * BYTES_PER_US
        nsplit = 0
        for i in 1:n
            split = g[i] > 0 && (policy == :all || g[i] / w + c < g[i])
            split && (nsplit += 1)
            gel[i] = split ? g[i] / w + c : Float64(g[i])
        end
        fin = zeros(Float64, n); pof = zeros(Int, n)
        for i in 1:n
            a, ai = 0.0, 0
            for j in plan.preds[i]
                fin[j] > a && (a = fin[j]; ai = j)
            end
            fin[i] = a + gel[i] + lu[i]; pof[i] = ai
        end
        Lw, iw = findmax(fin)
        p_g, p_lu, npath = 0.0, 0.0, 0
        let i = iw
            while i != 0
                p_g += gel[i]; p_lu += lu[i]; npath += 1; i = pof[i]
            end
        end
        Ww = W_node + nsplit * w * c                  # work incl. coordination
        Tw64 = max(Ww / 64, Lw)
        push!(coop_rows, (; w, policy, c_us, Lw, sweep_bound=L_node/Lw, nsplit,
            ntasks=n + nsplit*(w - 1), p_g, p_lu, npath, Ww, Tw64,
            j64_bound=Tn(64)/Tw64))
    end
    println("\n| w | policy | ovh (µs) | L_w (MB) | L_node/L_w | n_split | tasks/sweep | path GEMV (MB) | path LU (MB) | LU share | T64 bound |")
    println("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in coop_rows
        println("| $(r.w) | $(r.policy) | $(r.c_us) | $(mb(r.Lw)) | " *
            "$(round(r.sweep_bound, digits=2)) | $(r.nsplit) | $(r.ntasks) | " *
            "$(mb(r.p_g)) | $(mb(r.p_lu)) | $(round(r.p_lu/r.Lw, digits=2)) | " *
            "$(round(r.j64_bound, digits=2)) |")
    end
    open(joinpath(outdir, "fgs_dag_L_coopsplit_$(rung).csv"), "w") do io
        println(io, "w,policy,ovh_us,L_bytes,sweep_speedup_bound,n_split," *
            "tasks_per_sweep,path_gemv_bytes,path_lu_bytes,path_leaves," *
            "W_bytes,T64_bytes,j64_speedup_bound")
        for r in coop_rows
            println(io, "$(r.w),$(r.policy),$(r.c_us),$(r.Lw),$(r.sweep_bound)," *
                "$(r.nsplit),$(r.ntasks),$(r.p_g),$(r.p_lu),$(r.npath)," *
                "$(r.Ww),$(r.Tw64),$(r.j64_bound)")
        end
    end
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
