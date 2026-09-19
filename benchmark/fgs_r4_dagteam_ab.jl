#!/usr/bin/env julia
# 021 R4 dagteam split-executor A/B campaign driver (v23).
# Baseline arm = the ACCEPTED champion order (colored @ j16, v21); candidate
# arm = the dagteam split dual-layout executor (gate 2d) at DAGTEAM_PRECISION
# (default f32full; f32conv / f64 are the spec's fallback rungs).
# AB_MODE=calibrate  — derive BOTH twins of the retained lexicographic config
#   and calibrate each tolerance by the standard staircase (coloring changes
#   accumulation ordering; dagteam changes it further and — at reduced
#   precision — perturbs the operator, so neither tolerance carries). Every
#   accepted solve passes the independent evaluator at BC rel-L2 <= 1e-6:
#   this stage IS the spec's numerical gate for the selected precision rung.
# AB_MODE=trials     — uninstrumented performance trials, alternating batches
#   colored/dagteam; rankings use total time to accepted accuracy only.
# AB_MODE=activity   — one stage-activity-instrumented solve per order
#   (budget attribution only; never a performance trial). No perf counters.
include(joinpath(@__DIR__, "fgs_r4_counters.jl"))
using Statistics

# /proc does not exist off Linux; the local pre-submit gate executes every mode
# on macOS, where the observer then records spans with no per-thread ticks.
Sys.islinux() || @eval activity_snapshot() = Dict{Int,Tuple{Int,Int}}()

dagteam_precision_env() = get(ENV, "DAGTEAM_PRECISION", "f32full")

# Pure check: the A/B pair differs only in sweep order, calibrated tolerance,
# and the dagteam-only `dagteam_precision` field (the colored baseline must
# carry neither dagteam_precision nor chunks).
function ab_check_pair(colored, dagteam)
    colored["kind"] == "fgs" && colored["sweep_order"] == "colored" ||
        error("Expected the colored baseline FGS configuration")
    haskey(colored, "dagteam_precision") && error("Colored configuration must not carry dagteam_precision")
    haskey(colored, "chunks") && error("Colored configuration must not carry chunks")
    dagteam["sweep_order"] == "dagteam" || error("Expected a dagteam configuration")
    haskey(dagteam, "chunks") && error("Dagteam configuration must not carry chunks")
    Set(keys(dagteam)) == union(Set(keys(colored)), Set(["dagteam_precision"])) ||
        error("A/B configs must share fields up to the dagteam-only dagteam_precision field")
    for key in keys(colored)
        key in ("sweep_order", "tolerance") && continue
        colored[key] == dagteam[key] || error("A/B configs differ in $key")
    end
    colored["tolerance"] > 0 || error("Colored configuration requires calibrated tolerance")
    dagteam["tolerance"] > 0 || error("Dagteam configuration requires calibrated tolerance")
    dagteam["dagteam_precision"] in ("f64", "f32conv", "f32full") ||
        error("Dagteam configuration requires dagteam_precision in f64/f32conv/f32full")
    return nothing
end

function ab_prepare(c, label, records; crosscheck=true)
    reset_cold!()
    solver = cold_make(c)
    compile_row, _ = cold_trial(solver)
    cold_require(compile_row.eligible && compile_row.authoritative_evaluator == "certified_fmm",
        "$label compilation solve failed")
    warm, reference = cold_trial(solver; crosscheck)
    cold_require(warm.eligible && warm.authoritative_evaluator == "certified_fmm",
        "$label warmup failed")
    push!(records, (; sweep_order=label, phase="compile", compile_row...))
    push!(records, (; sweep_order=label, phase="warmup", warm...))
    return solver, reference
end

function ab_calibrate(lex, out)
    colored = copy(lex)
    colored["sweep_order"] = "colored"
    colored["tolerance"] = 0.0
    cold_calibrate!(colored, out)
    dagteam = copy(lex)
    dagteam["sweep_order"] = "dagteam"
    dagteam["dagteam_precision"] = dagteam_precision_env()
    dagteam["tolerance"] = 0.0
    cold_calibrate!(dagteam, out)
    ab_check_pair(colored, dagteam)
    cold_check_config(colored; selected=true)
    cold_check_config(dagteam; selected=true)
    cold_write_toml(joinpath(out, "colored_selected.toml"), colored)
    cold_write_toml(joinpath(out, "dagteam_selected.toml"), dagteam)
    records = NamedTuple[]
    summary = Dict{String,Any}(
        "lexicographic_tolerance" => lex["tolerance"],
        "dagteam_precision" => dagteam["dagteam_precision"])
    for (label, c) in (("colored", colored), ("dagteam", dagteam))
        solver, reference = ab_prepare(c, label, records)
        confirm, x = cold_trial(solver)
        agreement = norm(x-reference)/max(norm(reference), eps())
        cold_require(confirm.eligible && confirm.authoritative_evaluator == "certified_fmm" &&
            agreement <= 1e-8, "$label calibration confirmation failed the repeat gate")
        push!(records, (; sweep_order=label, phase="confirmation", confirm...))
        summary["$(label)_tolerance"] = c["tolerance"]
        summary["$(label)_iterations"] = confirm.iterations
        summary["$(label)_confirmation_repeat_delta"] = agreement
    end
    cold_csv(joinpath(out, "calibration_validation.csv"), records)
    cold_write_toml(joinpath(out, "calibration_summary.toml"), summary)
    return nothing
end

function ab_trials(colored, dagteam, out)
    warmups = NamedTuple[]
    arms = Tuple{String,Any,Vector{Float64}}[]
    for (label, c) in (("colored", colored), ("dagteam", dagteam))
        solver, reference = ab_prepare(c, label, warmups)
        push!(arms, (label, solver, reference))
        cold_csv(joinpath(out, "warmups.csv"), warmups)
    end
    # Informational: both orders are separately certified to COLD_TARGET; their
    # solutions need only agree at that level, not to the repeat gate.
    reference_delta = norm(arms[1][3]-arms[2][3])/max(norm(arms[1][3]), eps())
    reps = parse(Int, get(ENV, "COLD_AB_REPS", "10"))
    reps >= 10 || error("COLD_AB_REPS must be at least 10")
    batches = parse(Int, get(ENV, "COLD_AB_BATCHES", "4"))
    batches >= 2 || error("COLD_AB_BATCHES must be at least 2")
    rows = NamedTuple[]
    for batch in 1:2batches
        label, solver, reference = arms[isodd(batch) ? 1 : 2]
        for trial in 1:reps
            row, x = cold_trial(solver)
            agreement = norm(x-reference)/max(norm(reference), eps())
            cold_require(row.eligible && row.authoritative_evaluator == "certified_fmm" &&
                agreement <= 1e-8, "Performance trial failed correctness/agreement gate")
            push!(rows, (; batch, trial, sweep_order=label, instrumented=false,
                relative_solution_delta=agreement, row...))
            cold_csv(joinpath(out, "ab_trials.csv"), rows)
        end
    end
    summary = Dict{String,Any}("cross_order_solution_rel_l2" => reference_delta,
        "trials_per_order" => batches*reps,
        "dagteam_precision" => dagteam["dagteam_precision"])
    for (label, _, _) in arms
        times = [r.solve_seconds for r in rows if r.sweep_order == label]
        iters = unique(r.iterations for r in rows if r.sweep_order == label)
        summary["$(label)_median_solve_seconds"] = median(times)
        summary["$(label)_min_solve_seconds"] = minimum(times)
        summary["$(label)_iterations"] = iters
    end
    summary["dagteam_speedup_median"] =
        summary["colored_median_solve_seconds"] / summary["dagteam_median_solve_seconds"]
    cold_write_toml(joinpath(out, "ab_summary.toml"), summary)
    return nothing
end

function ab_activity(colored, dagteam, out)
    warmups = NamedTuple[]
    validations = NamedTuple[]
    for (label, c) in (("colored", colored), ("dagteam", dagteam))
        solver, reference = ab_prepare(c, label, warmups)
        cold_csv(joinpath(out, "warmups.csv"), warmups)
        activity = NamedTuple[]
        observer = activity_observer(activity)
        observer(:probe, :start); observer(:probe, :stop); empty!(activity)
        cold_assert_threads()
        cold_reset!(solver)
        start = time_ns()
        pnl._solve!(rotor, solver; stage_observer=observer)
        elapsed = (time_ns()-start)*1e-9
        cold_assert_threads()
        x = copy(rotor.strength[:,2])
        e = cold_validate(x)
        agreement = norm(x-reference)/max(norm(reference), eps())
        cold_require(all(isfinite, x) && solver.solved && e.accepted &&
            e.authoritative_evaluator == "certified_fmm" && agreement <= 1e-8,
            "$label activity solve failed correctness gates")
        push!(validations, (; sweep_order=label, diagnostic_seconds=elapsed,
            iterations=solver.niter, relative_solution_delta=agreement,
            solved=solver.solved, e...))
        cold_csv(joinpath(out, "stage_thread_activity_$(label).csv"), activity)
        cold_csv(joinpath(out, "activity_validation.csv"), validations)
    end
    cold_write_toml(joinpath(out, "activity_scope.toml"), Dict(
        "performance_trial" => false,
        "scope" => "one stage-instrumented prepared solve per sweep order; reset, setup and validation excluded",
        "activity_scope" => "coarse stage boundaries; nearfield_update groups leaf/product/scatter and updates",
        "activity_limitations" => "proc snapshots are sequential and quantized to CLK_TCK; short stages may resolve to zero; incomplete endpoints have cpu_ticks=-1"))
    return nothing
end

function ab_main()
    configs, out, stage = cold_initialize!()
    Base.invokelatest(run_dagteam_ab, configs, out, stage)
end

function run_dagteam_ab(configs, out, stage)
    mode = get(ENV, "AB_MODE", "")
    mode in ("calibrate", "trials", "activity") || error("AB_MODE must be calibrate, trials, or activity")
    length(configs) == 1 || error("Dagteam A/B requires exactly one retained configuration")
    lex = copy(only(configs))
    lex["kind"] == "fgs" && lex["sweep_order"] == "lexicographic" ||
        error("Expected the retained lexicographic FGS configuration")
    cold_calibrate_selected!(lex, out; selected=true)
    cold_write_toml(joinpath(out, "config.toml"), lex)
    if mode == "calibrate"
        ab_calibrate(lex, out)
    else
        colored_file = get(ENV, "COLORED_CONFIG", "")
        isabspath(colored_file) && isfile(colored_file) ||
            error("COLORED_CONFIG must be an existing absolute file")
        dagteam_file = get(ENV, "DAGTEAM_CONFIG", "")
        isabspath(dagteam_file) && isfile(dagteam_file) ||
            error("DAGTEAM_CONFIG must be an existing absolute file")
        colored = cold_check_config(TOML.parsefile(colored_file); selected=true)
        dagteam = cold_check_config(TOML.parsefile(dagteam_file); selected=true)
        ab_check_pair(colored, dagteam)
        cp(colored_file, joinpath(out, "input_colored_config.toml"))
        cp(dagteam_file, joinpath(out, "input_dagteam_config.toml"))
        cold_write_toml(joinpath(out, "config_colored.toml"), colored)
        cold_write_toml(joinpath(out, "config_dagteam.toml"), dagteam)
        mode == "trials" ? ab_trials(colored, dagteam, out) : ab_activity(colored, dagteam, out)
    end
    cold_write_toml(joinpath(out, "status.toml"), Dict("status"=>"completed"))
    return nothing
end

abspath(PROGRAM_FILE) == (@__FILE__) && ab_main()
