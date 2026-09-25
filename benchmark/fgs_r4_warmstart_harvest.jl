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
winA = 1:nt                      # startup revolution, INCLUDING its transient
n_steps = maximum(parse(Int, getcell(r, "step")) for r in rows)
winB = (3nt + 1):min(4nt, n_steps)   # fourth revolution, from its first step

statline(v) = isempty(v) ? "—" :
    @sprintf("%.4g ± %.2g (med %.4g)", mean(v), maximum(v) - minimum(v), median(v))

io = open(outmd, "w")
println(io, "# Warm-start R4 harvest: $(basename(rundir))\n")
println(io, "Windows (transients INCLUDED, Ryan 2026-09-24): A = steps ",
        "$(first(winA))-$(last(winA)), B = steps $(first(winB))-$(last(winB)); ",
        "$n_steps steps total.\n")

trace_io = open(joinpath(rundir, "harvest_traces.csv"), "w")
println(trace_io, "arm,step,t_solve,niter_first,t_project,ct,solved,bcerr_ok")

for win in (("A", winA), ("B", winB))
    wname, wrange = win
    println(io, "## Window $wname (steps $(first(wrange))-$(last(wrange)))\n")
    println(io, "| arm | t_solve [s] | niter_first | t_project [s] | ",
            "setup t_setup / t_prime [s] | unconverged | nsolves≠1 | ",
            "bcerr>tol | uncertified |")
    println(io, "|---|---|---|---|---|---|---|---|---|")
    for k in sort(collect(keys(byarm)); by=armname)
        rr = byarm[k]
        inwin = [r for r in rr if parse(Int, getcell(r, "step")) in wrange]
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
    println(trace_io, join([armname(k), getcell(r, "step"), getcell(r, "t_solve"),
        getcell(r, "niter_first"), getcell(r, "t_project"), getcell(r, "CT"),
        getcell(r, "solved"),
        (isnan(m) || isnan(t)) ? "" : string(m <= t)], ","))
end
close(trace_io)

# --- cross-arm solution agreement (snapshots) ---------------------------------
snaps = Dict{String,Matrix{Float64}}()
for f in filter(f -> endswith(f, "_strength_snapshots.toml"), readdir(rundir))
    meta = TOML.parsefile(joinpath(rundir, f))
    binf = joinpath(rundir, replace(f, ".toml" => ".bin"))
    isfile(binf) || continue
    data = Array{Float64}(undef, meta["ncells"], meta["nsteps"])
    read!(binf, data)
    key = armname((meta["config"], meta["warmstart"], meta["warmstart_order"]))
    snaps[key] = data
end
if haskey(snaps, "fgs_cold") && length(snaps) > 1
    println(io, "## Cross-arm solution agreement (rel-L2 vs fgs_cold, per step)\n")
    println(io, "| arm | winA mean | winA max | winB mean | winB max |")
    println(io, "|---|---|---|---|---|")
    ref = snaps["fgs_cold"]
    agree_io = open(joinpath(rundir, "harvest_solution_deltas.csv"), "w")
    println(agree_io, "arm,step,rel_l2_vs_fgs_cold")
    for (k, m) in sort(collect(snaps); by=first)
        k == "fgs_cold" && continue
        ns = min(size(m, 2), size(ref, 2))
        d = [sqrt(sum(abs2, m[:, i] .- ref[:, i]) /
                  max(sum(abs2, ref[:, i]), eps())) for i in 1:ns]
        foreach(i -> println(agree_io, "$k,$i,$(d[i])"), 1:ns)
        dA = d[intersect(winA, 1:ns)]; dB = d[intersect(winB, 1:ns)]
        f4(x) = @sprintf("%.3e", x)
        println(io, "| $k | ", isempty(dA) ? "—" : f4(mean(dA)), " | ",
                isempty(dA) ? "—" : f4(maximum(dA)), " | ",
                isempty(dB) ? "—" : f4(mean(dB)), " | ",
                isempty(dB) ? "—" : f4(maximum(dB)), " |")
    end
    close(agree_io)
    println(io, "\nKnown context: FGS and Krylov converge to slightly ",
            "different wake-on fixed points (~2e-3, ",
            "rigid_motion_tree_reuse_item.md §5) — reported, not chased.")
else
    println(io, "## Cross-arm solution agreement: snapshots missing — ",
            "cannot compute (SNAPSHOT_STRENGTHS off?)")
end
close(io)
println("harvest written: $outmd")
