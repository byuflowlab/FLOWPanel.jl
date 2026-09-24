#!/usr/bin/env julia
# 021 FGS scalability diagnostic — Stage 2 analysis (mechanism stage; shipped
# with the data). Reads a run_r4_fgs_stage2.slurm.sh run dir.
#
# Usage:
#   julia --project benchmark/fgs_stage2_analysis.jl <run_dir> [out_dir]
#
# Conclusions structure (binding):
#   * The discriminator uses the per-worker drain-loop aggregates:
#       - per-task busy time inflating with j at (fixed) task count
#         → memory bandwidth / locality inside the GEMV+LU bodies;
#       - flat per-task busy but idle+lockmgmt worker-time shares growing
#         → scheduling (narrow DAG width and/or queue-lock contention);
#       - a backoff recovery on top isolates the idle lock hammering.
#   * Aggregates are only attributed after the overhead gate: aggdiag vs
#     uninstrumented spin medians within 5% at every rung.
#   * Worker durations overlap: NEVER summed as elapsed time. They are
#     reported as shares of team-time (team_size × nonself_product elapsed).
#   * busy_max/min accumulate per-outer-iteration extrema (sum of per-block
#     extrema), an imbalance indicator, not a single-sweep extremum.
#   * Idle-policy A/Bs compare PAIRED medians (arm order alternates by pair).

using Statistics, Printf, TOML

const RUN = length(ARGS) >= 1 ? abspath(ARGS[1]) : error("usage: fgs_stage2_analysis.jl <run_dir> [out_dir]")
const OUT = length(ARGS) >= 2 ? abspath(ARGS[2]) : joinpath(RUN, "analysis")
isdir(RUN) || error("run dir not found: $RUN")
mkpath(OUT)

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

rows = NamedTuple[]
for name in sort(readdir(RUN))
    dir = joinpath(RUN, name)
    isdir(dir) || continue
    solves = joinpath(dir, "results", "stage1_solves.csv")
    isfile(solves) || continue
    statusf = joinpath(RUN, "STATUS_$name")
    isfile(statusf) || continue   # interrupted stage: no STATUS, exclude
    status = strip(read(statusf, String) |> String)
    summaryf = joinpath(dir, "results", "stage1_summary.toml")
    summary = isfile(summaryf) ? TOML.parsefile(summaryf) : Dict{String,Any}()
    j = get(summary, "julia_threads", -1)
    for r in read_csv(solves)
        push!(rows, (; stage=name, status, j, r...))
    end
end
isempty(rows) && error("no stage1_solves.csv found under $RUN")

function write_csv(path, rs)
    names = propertynames(first(rs))
    open(path, "w") do io
        println(io, join(string.(names), ','))
        for row in rs
            println(io, join((string(getproperty(row, k)) for k in names), ','))
        end
    end
end
write_csv(joinpath(OUT, "stage2_all_rows.csv"), rows)

trials = filter(r -> r.phase == "trial" && r.status == "ok", rows)

# ---- per-process medians ------------------------------------------------------
group_keys = [:stage, :j, :label, :dagteam_workers, :dagteam_idle, :instrumented, :block]
groups = Dict{Any,Vector{Any}}()
for r in trials
    push!(get!(groups, Tuple(getproperty(r, k) for k in group_keys), Any[]), r)
end
proc = NamedTuple[]
for (key, rs) in sort(collect(groups); by=first)
    t = [r.solve_seconds for r in rs]
    push!(proc, (; NamedTuple{Tuple(group_keys)}(key)...,
        n=length(t), median_s=median(t), min_s=minimum(t), max_s=maximum(t)))
end
write_csv(joinpath(OUT, "stage2_process_medians.csv"), proc)

report = IOBuffer()
println(report, "# 021 FGS Stage-2 analysis — $(basename(RUN))\n")
failed = sort(unique(r.stage for r in rows if r.status != "ok"))
isempty(failed) || println(report, "**FAILED/partial stages (excluded):** ", join(failed, ", "), "\n")

# spin, uninstrumented, w=0 medians by j (this node's anchor ladder)
anchor_med(j) = begin
    p = [q.median_s for q in proc if q.j == j && !q.instrumented &&
         q.dagteam_workers == 0 && q.dagteam_idle == "spin"]
    isempty(p) ? NaN : median(p)
end
js = sort(unique(q.j for q in proc))
println(report, "## Spin anchors (uninstrumented, w=0; j64 = section-C spin arms)\n")
println(report, "| j | median s | ratio vs prev rung |")
println(report, "|---|---------|--------------------|")
for (i, j) in enumerate(js)
    rr = i == 1 ? NaN : anchor_med(js[i-1]) / anchor_med(j)
    @printf(report, "| %d | %.3f | %.2f |\n", j, anchor_med(j), rr)
end

# ---- overhead gate ------------------------------------------------------------
agg_med(j, w) = begin
    p = [q.median_s for q in proc if q.j == j && q.instrumented &&
         q.dagteam_workers == w]
    isempty(p) ? NaN : median(p)
end
println(report, "\n## Aggregate-instrumentation overhead gate (≤5% each rung)\n")
println(report, "| config | uninstr s | aggdiag s | overhead % |")
println(report, "|--------|-----------|-----------|------------|")
overhead_ok = true
capanchor = [q.median_s for q in proc if q.j == 64 && !q.instrumented &&
             q.dagteam_workers == 16 && q.dagteam_idle == "spin" && q.label == "anchor-cap"]
gaterows = [("j=16 w=0", anchor_med(16), agg_med(16, 0)),
            ("j=32 w=0", anchor_med(32), agg_med(32, 0)),
            ("j=64 w=0", anchor_med(64), agg_med(64, 0)),
            ("j=64 w=16", isempty(capanchor) ? NaN : median(capanchor), agg_med(64, 16))]
for (tag, a, b) in gaterows
    ov = 100 * (b - a) / a
    (isfinite(ov) && abs(ov) <= 5) || (global overhead_ok = false)
    @printf(report, "| %s | %.3f | %.3f | %+.1f |\n", tag, a, b, ov)
end
println(report, "\nOverhead gate: $(overhead_ok ? "PASS" : "FAIL — aggregates below are NOT attributable")")

# ---- per-worker aggregate decomposition ---------------------------------------
const AGG = (:diag_dagteam_busy_lower_ns, :diag_dagteam_busy_back_ns,
    :diag_dagteam_lockmgmt_ns, :diag_dagteam_idle_ns,
    :diag_dagteam_busy_max_ns, :diag_dagteam_busy_min_ns,
    :diag_dagteam_empty_pops, :diag_dagteam_n_lower, :diag_dagteam_n_back,
    :diag_dagteam_team_size, :diag_nonself_product_ns)
agg_means(j, w) = begin
    rs = [r for r in trials if r.j == j && r.instrumented == true &&
          r.dagteam_workers == w]
    isempty(rs) && return nothing
    Dict(p => mean(getproperty(r, p) for r in rs) for p in AGG)
end
println(report, "\n## Per-worker drain-loop aggregates (aggdiag means per solve)\n")
println(report, "| config | team | n_lower | n_back | per-task lower µs | per-task back µs | busy share | lockmgmt share | idle share | unobs share | imbalance max/min |")
println(report, "|--------|------|---------|--------|-------------------|------------------|------------|----------------|------------|-------------|--------------------|")
decomp = Dict{Tuple{Int,Int},Any}()
for (j, w) in ((16, 0), (32, 0), (64, 0), (64, 16))
    m = agg_means(j, w)
    m === nothing && continue
    decomp[(j, w)] = m
    team = m[:diag_dagteam_team_size]
    teamtime = team * m[:diag_nonself_product_ns]           # worker-time budget
    busy = m[:diag_dagteam_busy_lower_ns] + m[:diag_dagteam_busy_back_ns]
    lockm = m[:diag_dagteam_lockmgmt_ns]
    idle = m[:diag_dagteam_idle_ns]
    unobs = max(0.0, teamtime - busy - lockm - idle)        # spawn/join tails, reduce, epoch spins
    @printf(report, "| j=%d w=%d | %.0f | %.0f | %.0f | %.1f | %.1f | %.3f | %.3f | %.3f | %.3f | %.2f |\n",
        j, w, team, m[:diag_dagteam_n_lower], m[:diag_dagteam_n_back],
        m[:diag_dagteam_busy_lower_ns] / max(m[:diag_dagteam_n_lower], 1) / 1e3,
        m[:diag_dagteam_busy_back_ns] / max(m[:diag_dagteam_n_back], 1) / 1e3,
        busy / teamtime, lockm / teamtime, idle / teamtime, unobs / teamtime,
        m[:diag_dagteam_busy_max_ns] / max(m[:diag_dagteam_busy_min_ns], 1))
end
println(report, """

Reading guide (shares are of team-time = team_size × nonself_product elapsed;
worker durations overlap and are never summed as elapsed):
- per-task lower µs inflating 16→64 at fixed n_lower → bandwidth/locality
  inside the GEMV+LU bodies;
- flat per-task busy, growing idle+lockmgmt shares → scheduling (narrow DAG
  and/or queue-lock contention);
- 'unobs' covers team spawn/join tails, the serial boundary reduce, and
  inter-sweep epoch spins.""")

if haskey(decomp, (16, 0)) && haskey(decomp, (64, 0))
    a, b = decomp[(16, 0)], decomp[(64, 0)]
    pt(m) = m[:diag_dagteam_busy_lower_ns] / max(m[:diag_dagteam_n_lower], 1)
    @printf(report, "\n**Discriminator:** per-task lower busy inflates ×%.2f from j=16 to j=64", pt(b) / pt(a))
    team(m) = m[:diag_dagteam_team_size]
    share(m) = (m[:diag_dagteam_lockmgmt_ns] + m[:diag_dagteam_idle_ns]) /
               (team(m) * m[:diag_nonself_product_ns])
    @printf(report, "; idle+lockmgmt team-time share moves %.3f → %.3f.\n", share(a), share(b))
end

# ---- paired idle-policy A/Bs --------------------------------------------------
for j in (64, 32)
    A = filter(q -> q.label == "idle-spin" && q.j == j, proc)
    B = filter(q -> q.label == "idle-backoff" && q.j == j, proc)
    (isempty(A) || isempty(B)) && continue
    println(report, "\n## Idle-policy A/B @ j=$j w=0 (spin vs backoff, paired medians)\n")
    println(report, "| pair | spin s | backoff s | Δ (backoff−spin) |")
    println(report, "|------|--------|-----------|-------------------|")
    deltas = Float64[]
    for pair in sort(unique(q.block for q in A))
        sa = [q.median_s for q in A if q.block == pair]
        sb = [q.median_s for q in B if q.block == pair]
        (isempty(sa) || isempty(sb)) && continue
        push!(deltas, only(sb) - only(sa))
        @printf(report, "| %d | %.3f | %.3f | %+.3f |\n", pair, only(sa), only(sb), only(sb) - only(sa))
    end
    isempty(deltas) || @printf(report, "\nPaired mean Δ = %+.3f s; all pairs same sign: %s\n",
        mean(deltas), all(>(0), deltas) || all(<(0), deltas))
end

# ---- backoff safety (single blocks; report only, no strong claims) -----------
println(report, "\n## Backoff safety checks (single blocks — deltas are indicative only)\n")
println(report, "| config | spin s | backoff s | Δ % |")
println(report, "|--------|--------|-----------|-----|")
saf16 = [q.median_s for q in proc if q.label == "idle-backoff-safety" && q.j == 16]
isempty(saf16) || @printf(report, "| j=16 w=0 | %.3f | %.3f | %+.1f |\n",
    anchor_med(16), only(saf16), 100 * (only(saf16) - anchor_med(16)) / anchor_med(16))
safcap = [q.median_s for q in proc if q.label == "idle-backoff-cap"]
capspin = isempty(capanchor) ? NaN : median(capanchor)
isempty(safcap) || @printf(report, "| j=64 w=16 | %.3f | %.3f | %+.1f |\n",
    capspin, only(safcap), 100 * (only(safcap) - capspin) / capspin)

verdict = joinpath(RUN, "backoff_verdict.txt")
isfile(verdict) && println(report, "\n`backoff_verdict.txt`: ", strip(read(verdict, String)))

write(joinpath(OUT, "stage2_report.md"), String(take!(report)))
println("Wrote $(joinpath(OUT, "stage2_report.md")) and merged CSVs ($(length(rows)) rows, $(length(proc)) process groups).")
