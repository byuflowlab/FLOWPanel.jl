## Adaptive DJI 9443 direct-potential panel-wake convergence driver.
##
## Orchestrates the wake-row convergence study of plans/20260721_panel_wake.md by
## running examples/rotor_panel_wake_study.jl as an isolated subprocess per case
## (clean checkpoint/resume: each case owns a persistent dir + COMPLETED marker).
##
## Stages (no 80x81):
##   1. Axial 40x40 ladder (constant 4.0 m/s freestream).
##   2. Hover 0.25 m/s 40x40 ladder (ease down from 4.0 m/s).
##   3. Hover zero-flow 40x40 single confirmation at the accepted cap.
##   4. Aggregate report.
##
## Completed valid cases are skipped. Under a wall-time budget the driver stops
## cleanly >=15 min before expiry and resumes at the first incomplete case on
## resubmission.
##
## Env: PANELWAKE_MAX_WALL_SECONDS (soft budget; default from SLURM or 11.25h),
##      THREADS (default Julia thread count).

import TOML

const REPO_ROOT   = normpath(joinpath(@__DIR__, ".."))
const STUDY_JL    = joinpath(@__DIR__, "rotor_panel_wake_study.jl")
const DATA_ROOT   = joinpath(REPO_ROOT, "data", "dji9443_panelwake")
const THREADS     = get(ENV, "THREADS", string(Base.Threads.nthreads()))

# ---- wall-time budget --------------------------------------------------------
const SAFETY_RESERVE_S = 15 * 60
const T_START = time()
function budget_seconds()
    v = tryparse(Float64, get(ENV, "PANELWAKE_MAX_WALL_SECONDS", ""))
    v !== nothing && return v
    # fall back to SLURM time-based default; if unknown assume 11.25h
    return 11.25 * 3600
end
const BUDGET_S = budget_seconds()
elapsed_s() = time() - T_START
remaining_s() = BUDGET_S - elapsed_s() - SAFETY_RESERVE_S

# ---- fixed study parameters --------------------------------------------------
const RPM  = "5400"
const NT   = 36
const START_VINF = "4.0"

# ---- case bookkeeping --------------------------------------------------------
struct Case
    name::String            # persistent dir basename under DATA_ROOT
    run_name::String        # RUN_NAME env -> data/<run_name>
    wake_rows::Int
    terminal_vinf::Float64
    decrease_revs::Int
    settle_revs::Int
end

case_dir(c::Case) = joinpath(REPO_ROOT, "data", c.run_name)

"A case is complete iff its COMPLETED marker exists AND its metadata parses with
a finite CT (never trust a bare VTK/CSV or process exit)."
function is_complete(c::Case)
    dir = case_dir(c)
    isfile(joinpath(dir, "COMPLETED")) || return false
    meta = joinpath(dir, "$(basename(c.run_name))_study_metadata.toml")
    isfile(meta) || return false
    try
        d = TOML.parsefile(meta)
        return get(d, "all_finite", false) === true &&
               isfinite(Float64(get(d, "CT_bernoulli_final_mean", NaN)))
    catch
        return false
    end
end

function read_metadata(c::Case)
    meta = joinpath(case_dir(c), "$(basename(c.run_name))_study_metadata.toml")
    return TOML.parsefile(meta)
end

"Final-revolution mean dCT/d(r/R) on fixed bins with r/R >= 0.1, per blade, from
the spanwise monitor CSV. Returns (rR::Vector, dTdr::Vector) sorted by rR."
function read_loading(c::Case)
    mon = joinpath(case_dir(c), "monitors")
    isdir(mon) || return (Float64[], Float64[])
    csvs = sort(filter(f -> occursin("_spanwise_system1.csv", f), readdir(mon)))
    isempty(csvs) && return (Float64[], Float64[])
    path = joinpath(mon, csvs[1])
    header = String[]; rows = Vector{Vector{String}}()
    open(path) do io
        header = split(strip(readline(io)), ",")
        for ln in eachline(io)
            isempty(strip(ln)) && continue
            push!(rows, String.(split(strip(ln), ",")))
        end
    end
    col(name) = findfirst(==(name), header)
    is, ib, ibc, it = col("step"), col("bin"), col("bin_center"), col("thrust")
    (is === nothing || ib === nothing || ibc === nothing || it === nothing) &&
        return (Float64[], Float64[])
    steps = [parse(Int, r[is]) for r in rows]
    laststep = maximum(steps)
    firststep = max(minimum(steps), laststep - NT + 1)
    R = 0.119
    binvals = Dict{Int,Vector{Float64}}()   # bin -> [bin_center, thrust] accumulation
    bincenters = Dict{Int,Float64}()
    for r in rows
        s = parse(Int, r[is]); (s >= firststep && s <= laststep) || continue
        b = parse(Int, r[ib]); bc = parse(Float64, r[ibc]); th = parse(Float64, r[it])
        push!(get!(binvals, b, Float64[]), th)
        bincenters[b] = bc
    end
    bins = sort(collect(keys(binvals)))
    rR = Float64[]; dTdr = Float64[]
    for b in bins
        rr = bincenters[b] / R
        rr >= 0.1 || continue
        push!(rR, rr)
        push!(dTdr, sum(binvals[b]) / length(binvals[b]))
    end
    order = sortperm(rR)
    return (rR[order], dTdr[order])
end

# ---- subprocess launch -------------------------------------------------------
function run_case(c::Case)
    dir = case_dir(c)
    if is_complete(c)
        println(">> SKIP (already complete): $(c.run_name)")
        return :skipped
    end
    env = copy(ENV)
    env["RUN_NAME"]                 = c.run_name
    env["RHPC_MESH"]                = "40_40"
    env["RPM"]                      = RPM
    env["NT"]                       = string(NT)
    env["WAKE_ROWS"]                = string(c.wake_rows)
    env["SETTLE_REVS"]              = string(c.settle_revs)
    env["START_VINF"]               = START_VINF
    env["TERMINAL_VINF"]            = string(c.terminal_vinf)
    env["FREESTREAM_DECREASE_REVS"] = string(c.decrease_revs)
    env["FORMULATION"]              = "direct"
    env["SAVE_VTK"]                 = "true"
    println("\n" * "="^78)
    println(">> RUN: $(c.run_name)  rows=$(c.wake_rows) Vterm=$(c.terminal_vinf) " *
            "decrease=$(c.decrease_revs) settle=$(c.settle_revs)")
    println("   budget remaining ~$(round(remaining_s()/60, digits=1)) min")
    println("="^78)
    t0 = time()
    cmd = Cmd(`julia --project=$(REPO_ROOT) -t $(THREADS) $(STUDY_JL)`;
              env=env, dir=REPO_ROOT)
    ok = success(pipeline(cmd; stdout=stdout, stderr=stderr))
    dt = time() - t0
    println(">> case wall time: $(round(dt/60, digits=2)) min  (exit ok=$(ok))")
    if !ok || !is_complete(c)
        println(">> case did NOT complete cleanly: $(c.run_name)")
        return :failed
    end
    return :done
end

# ---- convergence metrics -----------------------------------------------------
CT_change(a, b)      = abs(b - a) / max(abs(b), eps())
function loading_change(rRa, La, rRb, Lb)
    # metric on the finer (b) bins; require matching fixed bins
    (length(La) == length(Lb) && length(Lb) > 0) || return Inf
    return sqrt(sum((Lb .- La).^2)) / max(sqrt(sum(Lb.^2)), eps())
end

const CT_TOL = 0.02
const LOAD_TOL = 0.02

"Run a wake-row ladder with the two-consecutive-passing-refinement rule. Returns
(accepted_cap::Union{Int,Nothing}, results::Vector)."
function run_ladder(stage_name, make_case, rows_ladder)
    println("\n########## LADDER: $(stage_name) ##########")
    results = Tuple{Int,Any,Any,Any}[]   # (rows, meta, rR, L)
    passes = 0                            # consecutive passing refinements
    accepted = nothing
    for (i, rows) in enumerate(rows_ladder)
        if remaining_s() <= 0
            println(">> WALL BUDGET reached before rows=$(rows); stopping ladder cleanly.")
            return (accepted, results)
        end
        c = make_case(rows)
        status = run_case(c)
        if status == :failed
            println(">> ladder run failed at rows=$(rows); stopping.")
            return (accepted, results)
        end
        meta = read_metadata(c)
        rR, L = read_loading(c)
        stable = get(meta, "stable", false) === true
        CT = Float64(get(meta, "CT_bernoulli_final_mean", NaN))
        println(">> rows=$(rows)  CT=$(round(CT,sigdigits=6))  stable=$(stable)")
        push!(results, (rows, meta, rR, L))
        if length(results) >= 2
            (pr, pmeta, prR, pL) = results[end-1]
            ctc = CT_change(Float64(pmeta["CT_bernoulli_final_mean"]), CT)
            ldc = loading_change(prR, pL, rR, L)
            prev_stable = get(pmeta, "stable", false) === true
            passed = stable && prev_stable && ctc <= CT_TOL && ldc <= LOAD_TOL
            println("   refinement $(pr)->$(rows): dCT=$(round(ctc,sigdigits=3)) " *
                    "dLoad=$(round(ldc,sigdigits=3)) passed=$(passed)")
            passes = passed ? passes + 1 : 0
            if passes >= 2
                accepted = rows
                println(">> LADDER CONVERGED at cap rows=$(rows) " *
                        "(two consecutive passing refinements).")
                return (accepted, results)
            end
        end
    end
    println(">> ladder exhausted without two consecutive passing refinements.")
    return (accepted, results)
end

# =============================================================================
# STAGE 1 — axial 40x40 ladder (constant 4.0 m/s freestream)
# =============================================================================
axial_case(rows) = Case("axial_rows$(lpad(rows,3,'0'))",
    "dji9443_panelwake/axial_rows$(lpad(rows,3,'0'))",
    rows, 4.0, 0, 3)   # START=TERMINAL=4.0, no decrease, 3 settle revs
axial_ladder = [72, 144, 216, 288, 360, 432]
(axial_cap, axial_results) = run_ladder("axial J=0.1867", axial_case, axial_ladder)
# Fall back to the largest completed rung if the ladder did not formally converge.
axial_ref_cap = axial_cap === nothing ?
    (isempty(axial_results) ? nothing : axial_results[end][1]) : axial_cap

# =============================================================================
# STAGE 2 — hover 0.25 m/s ladder (ease down from 4.0 m/s)
# =============================================================================
round36(x) = max(36, Int(round(x / 36)) * 36)
hover_dir(vlabel, rows) = "hover_vinf$(vlabel)_rows$(lpad(rows,3,'0'))"
hover_decrease_revs = parse(Int, get(ENV, "HOVER_DECREASE_REVS", "4"))
hover025_case(rows) = Case(hover_dir("0p25", rows),
    "dji9443_panelwake/$(hover_dir("0p25", rows))",
    rows, 0.25, hover_decrease_revs, 3)

hover025_cap = nothing
hover025_results = []
if axial_ref_cap === nothing
    println("\n>> No axial reference cap available; skipping hover ladder.")
elseif remaining_s() <= 0
    println("\n>> WALL BUDGET reached; skipping hover ladder (resume on resubmission).")
else
    base = axial_ref_cap
    ladder = unique(round36.([0.5base, 1.0base, 1.5base]))
    # Extend up to 432 in steps equal to the ladder's initial spacing, so every
    # compared refinement has the same row increment (a smaller step would make
    # the 2% change gate systematically easier to pass).
    rung_step = length(ladder) >= 2 ? ladder[2] - ladder[1] : 72
    while ladder[end] < 432
        push!(ladder, min(432, ladder[end] + rung_step))
    end
    (hover025_cap, hover025_results) = run_ladder("hover 0.25 m/s", hover025_case, ladder)
end
hover025_ref_cap = hover025_cap === nothing ?
    (isempty(hover025_results) ? nothing : hover025_results[end][1]) : hover025_cap

# =============================================================================
# STAGE 3 — hover zero-flow single confirmation at the accepted cap
# =============================================================================
hover0_result = nothing
if hover025_ref_cap === nothing
    println("\n>> No hover 0.25 cap available; skipping zero-flow confirmation.")
elseif remaining_s() <= 0
    println("\n>> WALL BUDGET reached; skipping zero-flow (resume on resubmission).")
else
    zero_decrease_revs = parse(Int, get(ENV, "ZEROFLOW_DECREASE_REVS",
        string(hover_decrease_revs + 2)))   # gentler ramp approaching 0
    c0 = Case(hover_dir("0", hover025_ref_cap),
        "dji9443_panelwake/$(hover_dir("0", hover025_ref_cap))",
        hover025_ref_cap, 0.0, zero_decrease_revs, 3)
    status = run_case(c0)
    if status != :failed && is_complete(c0)
        hover0_result = read_metadata(c0)
    end
end

# =============================================================================
# STAGE 4 — aggregate report
# =============================================================================
println("\n\n" * "#"^78)
println("# AGGREGATE REPORT — DJI 9443 direct-potential panel-wake convergence")
println("#"^78)
# CCBlade axial references (from plan): ncrit=4 CT 0.0419673; ncrit=9 CT 0.0325479.
const CCBLADE_NCRIT4 = 0.0419673
const CCBLADE_NCRIT9 = 0.0325479

function report_ladder(title, results, cap)
    println("\n$(title):")
    println("  rows |    CT_bern |   CT_kj | ptp_rel | drift_rel | stable")
    for (rows, meta, _, _) in results
        ct  = Float64(get(meta, "CT_bernoulli_final_mean", NaN))
        ckj = Float64(get(meta, "CT_kj_final_mean", NaN))
        ptp = Float64(get(meta, "CT_bernoulli_final_ptp_rel", NaN))
        drf = Float64(get(meta, "CT_bernoulli_final_drift_rel", NaN))
        stb = get(meta, "stable", false)
        println("  $(lpad(rows,4)) | $(lpad(round(ct,sigdigits=6),10)) | " *
                "$(lpad(round(ckj,sigdigits=5),7)) | $(lpad(round(ptp,sigdigits=3),7)) | " *
                "$(lpad(round(drf,sigdigits=3),9)) | $(stb)")
    end
    println("  accepted cap: $(cap === nothing ? "NOT CONVERGED" : cap)")
end

report_ladder("STAGE 1 — Axial J=0.1867", axial_results, axial_cap)
if !isempty(axial_results)
    ct = Float64(axial_results[end][2]["CT_bernoulli_final_mean"])
    println("  vs CCBlade ncrit=4 ($(CCBLADE_NCRIT4)): CT/CCB = $(round(ct/CCBLADE_NCRIT4,digits=3))")
    println("  vs CCBlade ncrit=9 ($(CCBLADE_NCRIT9)): CT/CCB = $(round(ct/CCBLADE_NCRIT9,digits=3))")
end

report_ladder("STAGE 2 — Hover 0.25 m/s (J~0.0117)", hover025_results, hover025_cap)

println("\nSTAGE 3 — Hover zero-flow confirmation:")
if hover0_result === nothing
    println("  (not run or did not complete)")
else
    ct = Float64(hover0_result["CT_bernoulli_final_mean"])
    println("  cap rows=$(hover0_result["wake_rows_cap"])  CT=$(round(ct,sigdigits=6))  " *
            "stable=$(hover0_result["stable"])")
    println("  references: steady rigid wake ~0.0505, particle wake ~0.062, " *
            "viscous BEM ~0.068, experiment ~0.072")
end

println("\nNOTE: numerical wake-row convergence (the accepted caps above) is " *
        "distinct from physical agreement with the references. Report them " *
        "separately.")
println("\nAggregate elapsed: $(round(elapsed_s()/60, digits=1)) min " *
        "(budget $(round(BUDGET_S/60, digits=1)) min).")
println("Convergence driver finished.")
