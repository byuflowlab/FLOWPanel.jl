#!/usr/bin/env julia
# 021 FGS scalability diagnostic — Stage 1 analysis (plan C,
# fgs_scalability_diagnostic_plan_20260921c.md; shipped with the data per the
# plan's "ship the analysis script with the data" clause).
#
# Usage:
#   julia --project benchmark/fgs_stage1_analysis.jl <run_dir> [out_dir]
#
# <run_dir> is the launcher's $COLD_DATA_ROOT/fgs-stage1-<jobid> directory (or a
# local copy). Reads every  <stage>/results/stage1_solves.csv  under it plus the
# per-stage STATUS_* files; writes merged CSVs and a markdown report to
# [out_dir] (default <run_dir>/analysis).
#
# Conclusions structure (binding, from plan C):
#   * The 16->32 PLATEAU and the 64-thread REGRESSION are SEPARATE conclusions,
#     each gated on its own reproduction check against the Stage-0 baseline.
#   * If an effect fails to reproduce, this script reports the comparison
#     deltas and STOPS attributing that effect (no phase decomposition for it).
#   * Phase differences use matched-sample MEANS (additive decomposition);
#     headline solve times use per-process MEDIANS with dispersion.
#   * Instrumented rows are used only after the overhead gate: instrumented vs
#     uninstrumented ladder medians within 5% at every rung AND unchanged
#     scaling shape (same sign of successive-rung ratios).
#   * Instrumented and uninstrumented rows are never pooled (diag_* == -1
#     sentinels mark uninstrumented rows).

using Statistics, Printf

const RUN = length(ARGS) >= 1 ? abspath(ARGS[1]) : error("usage: fgs_stage1_analysis.jl <run_dir> [out_dir]")
const OUT = length(ARGS) >= 2 ? abspath(ARGS[2]) : joinpath(RUN, "analysis")
isdir(RUN) || error("run dir not found: $RUN")
mkpath(OUT)

# ---- minimal CSV reader (cold_csv writes plain comma-separated scalars) ------
function read_csv(path)
    lines = readlines(path)
    isempty(lines) && return NamedTuple[]
    names = Symbol.(split(lines[1], ','))
    rows = NamedTuple[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        vals = map(split(line, ',')) do s
            v = tryparse(Int, s); v !== nothing && return v
            v = tryparse(Float64, s); v !== nothing && return v
            s == "true" ? true : s == "false" ? false : String(s)
        end
        push!(rows, NamedTuple{Tuple(names)}(Tuple(vals)))
    end
    rows
end

# ---- gather ------------------------------------------------------------------
# stage dir name encodes the launcher's variant; j comes from the row-level
# julia thread count recorded in stage1_summary.toml
using TOML
stages = NamedTuple[]
for name in sort(readdir(RUN))
    dir = joinpath(RUN, name)
    isdir(dir) || continue
    solves = joinpath(dir, "results", "stage1_solves.csv")
    isfile(solves) || continue
    status = strip(read(joinpath(RUN, "STATUS_$name"), String) |> String)
    summaryf = joinpath(dir, "results", "stage1_summary.toml")
    summary = isfile(summaryf) ? TOML.parsefile(summaryf) : Dict{String,Any}()
    j = get(summary, "julia_threads", -1)
    rows = read_csv(solves)
    for r in rows
        push!(stages, (; stage=name, status, j, r...))
    end
end
isempty(stages) && error("no stage1_solves.csv found under $RUN")

# merged long-form table
function write_csv(path, rows)
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(names), ','))
        for row in rows
            println(io, join((string(getproperty(row, k)) for k in names), ','))
        end
    end
end
write_csv(joinpath(OUT, "stage1_all_rows.csv"), stages)

trials(rows) = filter(r -> r.phase == "trial", rows)
sel(; kw...) = filter(stages) do r
    all(getproperty(r, k) == v for (k, v) in pairs(kw))
end

# ---- per-process medians + dispersion ---------------------------------------
group_keys = [:stage, :status, :j, :label, :arm, :dagteam_workers, :instrumented, :block]
groups = Dict{Any,Vector{Any}}()
for r in trials(stages)
    key = Tuple(getproperty(r, k) for k in group_keys)
    push!(get!(groups, key, Any[]), r)
end
proc = NamedTuple[]
for (key, rs) in sort(collect(groups); by=first)
    t = [r.solve_seconds for r in rs]
    push!(proc, (; NamedTuple{Tuple(group_keys)}(key)...,
        n=length(t), median_s=median(t), min_s=minimum(t), max_s=maximum(t),
        spread_pct=100 * (maximum(t) - minimum(t)) / median(t),
        iterations=join(unique(r.iterations for r in rs), ';')))
end
write_csv(joinpath(OUT, "stage1_process_medians.csv"), proc)

report = IOBuffer()
println(report, "# 021 FGS Stage-1 analysis — $(basename(RUN))\n")
failed = unique((r.stage, r.status) for r in stages if r.status != "ok")
if !isempty(failed)
    println(report, "**FAILED/partial stages (excluded from conclusions):** ",
        join(first.(failed), ", "), "\n")
end
okproc = filter(p -> p.status == "ok", proc)

# ---- ladder medians (uninstrumented fixed-work arm) --------------------------
ladder = filter(p -> p.label == "ladder" && !p.instrumented, okproc)
med_by_j = Dict(j => median([p.median_s for p in ladder if p.j == j])
                for j in unique(p.j for p in ladder))
js = sort(collect(keys(med_by_j)))
println(report, "## Fixed-work ladder (uninstrumented, medians of block medians)\n")
println(report, "| j | median s | speedup vs j=1 | ratio vs prev rung |")
println(report, "|---|---------|----------------|--------------------|")
for (i, j) in enumerate(js)
    s1 = haskey(med_by_j, 1) ? med_by_j[1] / med_by_j[j] : NaN
    rr = i == 1 ? NaN : med_by_j[js[i-1]] / med_by_j[j]
    @printf(report, "| %d | %.3f | %.2f | %.2f |\n", j, med_by_j[j], s1, rr)
end

# ---- reproduction gates ------------------------------------------------------
# plateau: 32 fails to deliver over 16; regression: 64 slower than 32
plateau_repro = haskey(med_by_j, 16) && haskey(med_by_j, 32) &&
                med_by_j[16] / med_by_j[32] < 1.5   # <1.5x from a 2x thread doubling
regression_repro = haskey(med_by_j, 32) && haskey(med_by_j, 64) &&
                   med_by_j[64] > med_by_j[32]
println(report, "\n## Reproduction gates (SEPARATE conclusions)\n")
println(report, "- 16→32 plateau reproduced: **$plateau_repro**",
    haskey(med_by_j, 16) && haskey(med_by_j, 32) ?
    @sprintf(" (ratio %.2fx of the 2.00x ideal)", med_by_j[16] / med_by_j[32]) : " (rungs missing)")
println(report, "- 64-thread regression reproduced: **$regression_repro**",
    haskey(med_by_j, 32) && haskey(med_by_j, 64) ?
    @sprintf(" (64/32 time ratio %.2f)", med_by_j[64] / med_by_j[32]) : " (rungs missing)")
println(report, "\nA non-reproduced effect gets NO decomposition below — report the deltas and stop.")

# ---- instrumentation overhead gate ------------------------------------------
diagladder = filter(p -> p.label == "ladderdiag" && p.instrumented, okproc)
dmed_by_j = Dict(j => median([p.median_s for p in diagladder if p.j == j])
                 for j in unique(p.j for p in diagladder))
println(report, "\n## Instrumented-vs-uninstrumented overhead gate (≤5% each rung, unchanged shape)\n")
println(report, "| j | uninstr s | instr s | overhead % |")
println(report, "|---|-----------|---------|------------|")
overhead_ok = true
for j in js
    global overhead_ok
    haskey(dmed_by_j, j) || (overhead_ok = false; continue)
    ov = 100 * (dmed_by_j[j] - med_by_j[j]) / med_by_j[j]
    abs(ov) <= 5 || (overhead_ok = false)
    @printf(report, "| %d | %.3f | %.3f | %+.1f |\n", j, med_by_j[j], dmed_by_j[j], ov)
end
shape(m) = [sign(m[js[i+1]] - m[js[i]]) for i in 1:length(js)-1 if haskey(m, js[i]) && haskey(m, js[i+1])]
shape_ok = shape(med_by_j) == shape(dmed_by_j)
println(report, "\nOverhead gate: $(overhead_ok ? "PASS" : "FAIL"); shape unchanged: $(shape_ok ? "PASS" : "FAIL")")
attribute = overhead_ok && shape_ok

# ---- phase decomposition (matched-sample MEANS, instrumented rows only) ------
const PHASES = (:diag_fmm_ns, :diag_influence_mapping_ns, :diag_residual_ns,
    :diag_leaf_solve_ns, :diag_nonself_product_ns, :diag_scatter_ns,
    :diag_remaining_iteration_ns, :diag_final_update_ns, :diag_initialization_ns)
const DAGSPLITS = (:diag_dagteam_spawn_ns, :diag_dagteam_join_ns,
    :diag_dagteam_wait_ns, :diag_dagteam_reduce_ns)

phase_means(j) = begin
    rs = [r for r in trials(stages) if r.status == "ok" && r.label == "ladderdiag" &&
          r.instrumented == true && r.j == j]
    isempty(rs) && return nothing
    Dict(p => mean(getproperty(r, p) for r in rs) / 1e9 for p in (PHASES..., DAGSPLITS..., :diag_total_ns))
end

if attribute
    println(report, "\n## Exclusive-phase reconciliation + doubling shortfall (matched means, s)\n")
    pm = Dict(j => phase_means(j) for j in js)
    hdr = join(string.(PHASES), " | ")
    println(report, "| j | ", hdr, " | phase sum | diag_total | recon % |")
    println(report, "|---|", repeat("---|", length(PHASES) + 3))
    for j in js
        m = pm[j]; m === nothing && continue
        psum = sum(m[p] for p in PHASES)
        recon = 100 * psum / m[:diag_total_ns]
        print(report, "| $j | ", join((@sprintf("%.3f", m[p]) for p in PHASES), " | "))
        @printf(report, " | %.3f | %.3f | %.1f |\n", psum, m[:diag_total_ns], recon)
    end
    for (a, b, gate) in ((16, 32, plateau_repro), (32, 64, regression_repro))
        (haskey(pm, a) && pm[a] !== nothing && haskey(pm, b) && pm[b] !== nothing) || continue
        effect = a == 16 ? "plateau" : "regression"
        if !gate
            println(report, "\n### T_phase($b)−T_phase($a): SKIPPED — $effect not reproduced.")
            continue
        end
        println(report, "\n### T_phase($b) − T_phase($a) and shortfall from doubling (s)\n")
        println(report, "| phase | T($a) | T($b) | Δ = T($b)−T($a) | shortfall = T($b)−T($a)/2 |")
        println(report, "|-------|-------|-------|------------------|----------------------------|")
        for p in (PHASES..., DAGSPLITS...)
            ta, tb = pm[a][p], pm[b][p]
            @printf(report, "| %s | %.3f | %.3f | %+.3f | %+.3f |\n",
                string(p)[6:end], ta, tb, tb - ta, tb - ta / 2)
        end
    end
else
    println(report, "\nPhase decomposition SKIPPED: instrumentation gate failed.")
end

# ---- paired A/Bs (placement, worker cap) -------------------------------------
function paired_section(title, alabel, blabel)
    A = filter(p -> p.label == alabel, okproc)
    B = filter(p -> p.label == blabel, okproc)
    (isempty(A) || isempty(B)) && return
    println(report, "\n## $title\n")
    println(report, "| pair | $alabel s | $blabel s | Δ (B−A) |")
    println(report, "|------|-----------|-----------|---------|")
    deltas = Float64[]
    for pair in sort(unique(p.block for p in A))
        a = [p.median_s for p in A if p.block == pair]
        b = [p.median_s for p in B if p.block == pair]
        (isempty(a) || isempty(b)) && continue
        d = only(b) - only(a)
        push!(deltas, d)
        @printf(report, "| %d | %.3f | %.3f | %+.3f |\n", pair, only(a), only(b), d)
    end
    isempty(deltas) || @printf(report, "\nPaired mean Δ = %+.3f s; all pairs same sign: %s\n",
        mean(deltas), all(>(0), deltas) || all(<(0), deltas))
end
paired_section("Placement A/B @ j=64 (champion socket-0 vs both-socket)",
    "placement-champion", "placement-alt")

for cap in (16, 32)
    capA = filter(p -> p.label == "cap$cap-arm" && p.dagteam_workers == 0, okproc)
    capB = filter(p -> p.label == "cap$cap-arm" && p.dagteam_workers == cap, okproc)
    if isempty(capA) || isempty(capB)
        cap == 32 && println(report, "\n(cap=32 pairs absent — launcher skipped them because cap=16 did not beat baseline, or the stage failed.)")
        continue
    end
    println(report, "\n## Worker-cap A/B @ j=64 (cap=$cap vs all-threads, champion placement)\n")
    println(report, "| pair | workers=0 s | workers=$cap s | Δ |")
    println(report, "|------|-------------|--------------|---|")
    deltas = Float64[]
    for pair in sort(unique(p.block for p in capA))
        a = [p.median_s for p in capA if p.block == pair]
        b = [p.median_s for p in capB if p.block == pair]
        (isempty(a) || isempty(b)) && continue
        push!(deltas, only(b) - only(a))
        @printf(report, "| %d | %.3f | %.3f | %+.3f |\n", pair, only(a), only(b), only(b) - only(a))
    end
    isempty(deltas) || @printf(report, "\nPaired mean Δ = %+.3f s; all pairs same sign: %s\n",
        mean(deltas), all(>(0), deltas) || all(<(0), deltas))
end

# ---- accepted bridge ---------------------------------------------------------
acc = filter(p -> p.label == "accepted", okproc)
if !isempty(acc)
    println(report, "\n## Accepted-accuracy bridge (ties new pins to the Stage-0 verified ladder)\n")
    println(report, "| j | median s | iterations |")
    println(report, "|---|----------|------------|")
    for p in sort(acc; by=p -> p.j)
        @printf(report, "| %d | %.3f | %s |\n", p.j, p.median_s, p.iterations)
    end
    println(report, "\nCaveat: accepted solves run ONE extra fmm+influence+residual vs the")
    println(report, "fixed-work arm (convergence detected at iteration 28's check).")
end

write(joinpath(OUT, "stage1_report.md"), String(take!(report)))
println("Wrote $(joinpath(OUT, "stage1_report.md")) and merged CSVs ($(length(stages)) rows, $(length(proc)) process groups).")
