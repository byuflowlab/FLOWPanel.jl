#=##############################################################################
p033 R4 thread-scaling harvest: FGS vs ILU-GMRES-nfcache across j=1/8/16/32/64
at fixed (transplanted) champion knobs. Scrapes the Job A (cold) and Job B
(winB warm) run dirs and prints markdown comparison tables.

Usage:
  julia benchmark/p033_ts_harvest.jl <root>

<root> holds the harvested run dirs (rsynced from the data root):
  p033ts-cold-j<J>-<jobid>/fgs-setup/fgs_setup_ab_R4_j<J>.csv
  p033ts-cold-j<J>-<jobid>/ilu/phase2.csv
  p033ts-winb-j<J>-<jobid>/unsteady.csv

Knob provenance is read from the CSVs themselves; every table footnotes the
fixed-across-j transplant and the winB history-fill caveat (first order+1
restarted steps are effectively cold).
=###############################################################################

using Statistics

root = isempty(ARGS) ? pwd() : ARGS[1]
const JS = [1, 8, 16, 32, 64]

readcsv(path) = begin
    isfile(path) || return nothing
    lines = readlines(path)
    isempty(lines) && return nothing
    cols = Dict(String(c) => i for (i, c) in enumerate(split(lines[1], ",")))
    (cols, [split(l, ",") for l in lines[2:end] if !isempty(strip(l))])
end

finddir(pat) = begin
    hits = filter(d -> startswith(d, pat) && isdir(joinpath(root, d)),
                  readdir(root))
    isempty(hits) ? nothing : joinpath(root, sort(hits)[end])  # newest job id
end

fmt(x; d=2) = x === nothing ? "—" : string(round(x; digits=d))

# ---- Job A: cold ------------------------------------------------------------
println("## Job A — cold: one-time setup + isolated cold solve (R4, fixed knobs)\n")
println("| j | FGS ctor serial [s] | FGS ctor threaded [s] | FGS cold solve [s] | " *
        "ILU-nfc setup [s] | ILU factorization [s] | ILU nfcache build [s] | " *
        "ILU-nfc cold solve [s] |")
println("|---|---:|---:|---:|---:|---:|---:|---:|")
for J in JS
    d = finddir("p033ts-cold-j$J-")
    ctor_old = ctor_new = solve_fgs = nothing
    ilu_setup = ilu_fact = ilu_nfc = ilu_solve = nothing
    if d !== nothing
        ab = readcsv(joinpath(d, "fgs-setup", "fgs_setup_ab_R4_j$J.csv"))
        if ab !== nothing
            cols, rows = ab
            vals(arm, phase) = [parse(Float64, r[cols["t_s"]]) for r in rows
                                if r[cols["arm"]] == arm &&
                                   r[cols["phase"]] == phase &&
                                   r[cols["warmup"]] == "0"]
            v = vals("old", "full_ctor"); isempty(v) || (ctor_old = minimum(v))
            v = vals("new", "full_ctor"); isempty(v) || (ctor_new = minimum(v))
            v = vals("new", "cold_solve"); isempty(v) || (solve_fgs = minimum(v))
        end
        p2 = readcsv(joinpath(d, "ilu", "phase2.csv"))
        if p2 !== nothing
            cols, rows = p2
            for r in rows
                r[cols["config"]] == "krylov_ilu_nfcache" || continue
                ilu_setup = parse(Float64, r[cols["t_setup_total"]])
                ilu_fact  = tryparse(Float64, r[cols["t_setup_factorization"]])
                ilu_nfc   = tryparse(Float64, r[cols["nfcache_build_time"]])
                ilu_solve = parse(Float64, r[cols["t_solve_min"]])
            end
        end
    end
    println("| $J | $(fmt(ctor_old)) | $(fmt(ctor_new)) | $(fmt(solve_fgs; d=3)) | " *
            "$(fmt(ilu_setup)) | $(fmt(ilu_fact)) | $(fmt(ilu_nfc)) | " *
            "$(fmt(ilu_solve; d=3)) |")
end
println("\nFGS = retained R4 champion (P8/MAC0.4/leaf100/inner3, dagteam " *
        "f32full, champion tolerance); ILU-nfc apply knobs P12/MAC0.55/leaf48 " *
        "@ 500 GiB. Knobs fixed across j (transplanted j64 set — Ryan " *
        "2026-10-02 ruling, no re-tune).\n")

# ---- Job B: warm (winB) ------------------------------------------------------
println("## Job B — warm rotor hover, winB rev-4 legs (018-ported physics, " *
        "proj2 arms)\n")
println("| j | arm | t_setup [s] | t_prime [s] | median t_solve/step [s] | " *
        "niter_first (step 109) | steps solved |")
println("|---|---|---:|---:|---:|---:|---|")
for J in JS
    d = finddir("p033ts-winb-j$J-")
    u = d === nothing ? nothing : readcsv(joinpath(d, "unsteady.csv"))
    for cfg in ("fgs", "krylov_ilu_nfcache")
        tset = tprime = med = nif = nothing; nsolved = ntot = 0
        if u !== nothing
            cols, rows = u
            lr = [r for r in rows if r[cols["config"]] == cfg]
            if !isempty(lr)
                ts = [parse(Float64, r[cols["t_solve"]]) for r in lr]
                med = median(ts)
                tset = tryparse(Float64, lr[1][cols["t_setup"]])
                tprime = tryparse(Float64, lr[1][cols["t_setup_prime"]])
                nif = tryparse(Float64, lr[1][cols["niter_first"]])
                ntot = length(lr)
                nsolved = count(r -> lowercase(r[cols["solved"]]) == "true", lr)
            end
        end
        arm = cfg == "fgs" ? "fgs_proj2" : "ilu_nfcache_proj2"
        println("| $J | $arm | $(fmt(tset)) | $(fmt(tprime)) | $(fmt(med; d=3)) | " *
                "$(fmt(nif; d=0)) | $nsolved/$ntot |")
    end
end
println("\nCaveat (protocol): warm-start histories are not serialized in the " *
        "checkpoint — the first (order+1) steps of each restarted leg are " *
        "effectively cold and are INCLUDED in the medians above.")
