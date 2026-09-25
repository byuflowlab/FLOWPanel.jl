#=##############################################################################
BRAINSTORM 021 — warm-start R4 head-to-head harvest
(fgs_warmstart_r4_reset_prompt_20260924.md "Harvest deliverable").

Reads one campaign run directory (the launcher's $run: shared unsteady.csv +
per-arm strength snapshots + STATUS_*/COMPLETED_* sentinels) and emits
markdown tables + per-step trace CSVs:

  - per-arm Window A (steps 1-36) and Window B (steps 109-144) summaries:
    t_solve and niter_first mean ± spread (max-min) and median, per-step
    t_project, plus the one-time setup columns (t_setup, t_setup_prime,
    setup_detail) in their own columns — never folded into per-step numbers;
  - convergence-promise audit: count of steps with solved=false, nsolves!=1,
    bcerr_max > bcerr_tol, or bcerr_certified=false, per arm and window;
  - per-step traces (harvest_traces.csv): step, arm, t_solve, niter_first,
    t_project, CT — the transient SHAPE is a deliverable, not noise;
  - cross-arm solution agreement: for each step, rel-L2 delta of each arm's
    solved strengths vs the fgs_cold reference (and ilu_nfcache_cold vs
    fgs_cold explicitly — the known ~2e-3 wake-on fixed-point discrepancy is
    REPORTED, not chased);
  - CT traces per arm.

Usage:
  julia --project benchmark/fgs_r4_warmstart_harvest.jl <run_dir> [out.md]

No conclusions are drawn here beyond the tables: the compete/no-compete
recommendation is written by the harvesting agent in
fgs_warmstart_r4_results_<date>.md, and Ryan rules.
=###############################################################################

using TOML, Statistics, Printf

length(ARGS) >= 1 || error("usage: fgs_r4_warmstart_harvest.jl <run_dir> [out.md]")
rundir = abspath(ARGS[1])
outmd = length(ARGS) >= 2 ? ARGS[2] : joinpath(rundir, "harvest_summary.md")
csvpath = joinpath(rundir, "unsteady.csv")
isfile(csvpath) || error("no unsteady.csv in $rundir")

# --- CSV parse (header-keyed; cells may contain no commas by construction) ----
lines = readlines(csvpath)
cols = Dict(name => i for (i, name) in enumerate(split(lines[1], ',')))
rows = [split(l, ',') for l in lines[2:end]]
getcell(r, name) = String(r[cols[name]])
num(r, name) = (s = getcell(r, name); isempty(s) ? NaN : parse(Float64, s))

arms = unique(getcell(r, "config") * ":" * getcell(r, "warmstart") *
              getcell(r, "warmstart_order") for r in rows)

# arm identity as written by the campaign driver: RUN_NAME carries it, but the
# CSV key is (config, warmstart, warmstart_order)
armkey(r) = (getcell(r, "config"), getcell(r, "warmstart"),
             parse(Int, getcell(r, "warmstart_order")))
armname(k) = k[1] == "fgs" ?
    (k[2] == "cold" ? "fgs_cold" : k[2] == "prev" ? "fgs_prev" : "fgs_proj$(k[3])") :
    (k[2] == "cold" ? "ilu_nfcache_cold" : k[2] == "prev" ? "ilu_nfcache_prev" :
     "ilu_nfcache_proj$(k[3])")

byarm = Dict{Tuple{String,String,Int},Vector{Vector{SubString{String}}}}()
for r in rows
    push!(get!(byarm, armkey(r), Vector{Vector{SubString{String}}}()), r)
end

nt = isempty(rows) ? 36 : parse(Int, getcell(rows[1], "nt"))
restart_of(r) = parse(Int, getcell(r, "restart_step"))
step_of(r) = parse(Int, getcell(r, "step"))
# Global step index: a restarted (winB) leg's CSV rows are locally numbered
# 1..n from the restart point.
gstep(r) = restart_of(r) >= 0 ? restart_of(r) + step_of(r) : step_of(r)

# Window membership under the checkpoint+restart layout (2026-09-24):
#   A — non-restarted rows in steps 1..NT (ckpt legs' first revolution, winA
#       legs, or a full-march arm's first revolution);
#   B — restarted rows (winB legs; the fourth revolution), falling back to a
#       full-march arm's steps 3NT+1..4NT when no restart legs exist for it.
inA(r) = restart_of(r) < 0 && step_of(r) <= nt
function winB_rows(rr)
    restarted = [r for r in rr if restart_of(r) >= 0]
    !isempty(restarted) && return restarted
    [r for r in rr if 3nt < step_of(r) <= 4nt]
end

statline(v) = isempty(v) ? "—" :
    @sprintf("%.4g ± %.2g (med %.4g)", mean(v), maximum(v) - minimum(v), median(v))

io = open(outmd, "w")
println(io, "# Warm-start R4 harvest: $(basename(rundir))\n")
println(io, "Windows (transients INCLUDED, Ryan 2026-09-24): A = steps 1-$nt ",
        "(from the very first step), B = the fourth revolution (steps ",
        "$(3nt + 1)-$(4nt), restarted from the family checkpoint at step $(3nt) ",
        "where restart legs exist).\n")
println(io, "**Reporting note (Ryan 2026-09-24):** solver warm-start histories ",
        "are not serialized in the restart checkpoint, so the first ",
        "(order+1) steps of each restarted WARM leg are effectively cold — ",
        "that history-fill transient is INSIDE Window B's transient-included ",
        "statistics by design (never excluded). Interpret the first few ",
        "Window-B steps of warm arms accordingly; the per-step traces make ",
        "the refill visible.\n")

trace_io = open(joinpath(rundir, "harvest_traces.csv"), "w")
println(trace_io, "arm,global_step,restart_step,t_solve,niter_first,t_project,ct,solved,bcerr_ok")

for wname in ("A", "B")
    println(io, "## Window $wname\n")
    println(io, "| arm | t_solve [s] | niter_first | t_project [s] | ",
            "setup t_setup / t_prime [s] | unconverged | nsolves≠1 | ",
            "bcerr>tol | uncertified |")
    println(io, "|---|---|---|---|---|---|---|---|---|")
    for k in sort(collect(keys(byarm)); by=armname)
        rr = byarm[k]
        inwin = wname == "A" ? [r for r in rr if inA(r)] : winB_rows(rr)
        isempty(inwin) && continue
        ts = [num(r, "t_solve") for r in inwin]
        nf = [num(r, "niter_first") for r in inwin]
        tp = filter(!isnan, [num(r, "t_project") for r in inwin])
        setup = filter(!isempty, [getcell(r, "t_setup") for r in rr])
        prime = filter(!isempty, [getcell(r, "t_setup_prime") for r in rr])
        nbad_solved = count(r -> getcell(r, "solved") != "true", inwin)
        nbad_ns = count(r -> getcell(r, "nsolves") ∉ ("1", "-1"), inwin)
        nbad_bc = count(inwin) do r
            m, t = num(r, "bcerr_max"), num(r, "bcerr_tol")
            !isnan(m) && !isnan(t) && m > t
        end
        nuncert = count(inwin) do r
            c = getcell(r, "bcerr_certified")
            !isempty(c) && c != "true"
        end
        println(io, "| ", armname(k), " | ", statline(ts), " | ", statline(nf),
                " | ", isempty(tp) ? "—" : statline(tp), " | ",
                (isempty(setup) ? "—" : setup[1]), " / ",
                (isempty(prime) ? "—" : prime[1]),
                " | $nbad_solved | $nbad_ns | $nbad_bc | $nuncert |")
    end
    println(io)
end

for k in sort(collect(keys(byarm)); by=armname), r in byarm[k]
    m, t = num(r, "bcerr_max"), num(r, "bcerr_tol")
    println(trace_io, join([armname(k), gstep(r), restart_of(r),
        getcell(r, "t_solve"),
        getcell(r, "niter_first"), getcell(r, "t_project"), getcell(r, "CT"),
        getcell(r, "solved"),
        (isnan(m) || isnan(t)) ? "" : string(m <= t)], ","))
end
close(trace_io)

# --- cross-arm solution agreement (snapshots, leg-aware) ----------------------
# One snapshot file per (arm, leg). Keyed by (arm, window): window A snapshots
# come from non-restarted legs (ckpt/winA/full: columns = global steps 1..n);
# window B from restarted legs (columns = global steps restart_step+1 ..), or
# from a full-march leg's columns 3nt+1:4nt as a fallback. References:
# fgs_cold's matching window.
snapsA = Dict{String,Matrix{Float64}}()
snapsB = Dict{String,Matrix{Float64}}()
for f in filter(f -> endswith(f, "_strength_snapshots.toml"), readdir(rundir))
    meta = TOML.parsefile(joinpath(rundir, f))
    binf = joinpath(rundir, replace(f, ".toml" => ".bin"))
    isfile(binf) || continue
    data = Array{Float64}(undef, meta["ncells"], meta["nsteps"])
    read!(binf, data)
    key = armname((meta["config"], meta["warmstart"], meta["warmstart_order"]))
    rs = get(meta, "restart_step", -1)
    if rs >= 0
        snapsB[key] = data                       # winB leg: rev 4 columns
    else
        haskey(snapsA, key) || (snapsA[key] = data)  # first NT columns = window A
        # full-march layout fallback: its rev-4 columns double as window B
        size(data, 2) >= 4nt && !haskey(snapsB, key) &&
            (snapsB[key] = data[:, (3nt + 1):(4nt)])
    end
end
if (haskey(snapsA, "fgs_cold") || haskey(snapsB, "fgs_cold")) &&
        length(union(keys(snapsA), keys(snapsB))) > 1
    println(io, "## Cross-arm solution agreement (rel-L2 vs fgs_cold, per step)\n")
    println(io, "| arm | winA mean | winA max | winB mean | winB max |")
    println(io, "|---|---|---|---|---|")
    agree_io = open(joinpath(rundir, "harvest_solution_deltas.csv"), "w")
    println(agree_io, "arm,window,local_step,rel_l2_vs_fgs_cold")
    deltas(m, ref, cap) = begin
        ns = min(size(m, 2), size(ref, 2), cap)
        [sqrt(sum(abs2, m[:, i] .- ref[:, i]) /
              max(sum(abs2, ref[:, i]), eps())) for i in 1:ns]
    end
    f4(x) = @sprintf("%.3e", x)
    for k in sort(collect(union(keys(snapsA), keys(snapsB))))
        k == "fgs_cold" && continue
        dA = haskey(snapsA, k) && haskey(snapsA, "fgs_cold") ?
            deltas(snapsA[k], snapsA["fgs_cold"], nt) : Float64[]
        dB = haskey(snapsB, k) && haskey(snapsB, "fgs_cold") ?
            deltas(snapsB[k], snapsB["fgs_cold"], nt) : Float64[]
        foreach(i -> println(agree_io, "$k,A,$i,$(dA[i])"), eachindex(dA))
        foreach(i -> println(agree_io, "$k,B,$i,$(dB[i])"), eachindex(dB))
        println(io, "| $k | ", isempty(dA) ? "—" : f4(mean(dA)), " | ",
                isempty(dA) ? "—" : f4(maximum(dA)), " | ",
                isempty(dB) ? "—" : f4(mean(dB)), " | ",
                isempty(dB) ? "—" : f4(maximum(dB)), " |")
    end
    close(agree_io)
    println(io, "\nKnown context: FGS and Krylov converge to slightly ",
            "different wake-on fixed points (~2e-3, ",
            "rigid_motion_tree_reuse_item.md §5) — reported, not chased. ",
            "Window-B deltas within a solver family share the family ",
            "checkpoint, so they isolate the initial guess; ilu-vs-fgs ",
            "window-B deltas ALSO carry the two checkpoints' divergence.")
else
    println(io, "## Cross-arm solution agreement: snapshots missing — ",
            "cannot compute (SNAPSHOT_STRENGTHS off?)")
end
close(io)
println("harvest written: $outmd")
