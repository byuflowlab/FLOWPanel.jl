#!/usr/bin/env julia
# 021 R4 colored-sweep A/B campaign driver (v21).
# AB_MODE=calibrate  — derive the colored twin of the retained lexicographic
#   config and calibrate its tolerance by the standard staircase (coloring
#   changes accumulation ordering, so the retained tolerance does not carry).
# AB_MODE=trials     — uninstrumented performance trials, alternating batches
#   lexicographic/colored; rankings use total time to accepted accuracy only.
# AB_MODE=activity   — one stage-activity-instrumented solve per order
#   (budget attribution only; never a performance trial). No perf counters.
include(joinpath(@__DIR__, "fgs_r4_counters.jl"))
using Statistics

# /proc does not exist off Linux; the local pre-submit gate executes every mode
# on macOS, where the observer then records spans with no per-thread ticks.
Sys.islinux() || @eval activity_snapshot() = Dict{Int,Tuple{Int,Int}}()

# Pure check: the A/B pair differs only in sweep order and calibrated tolerance.
function ab_check_pair(lex, colored)
    lex["kind"] == "fgs" && lex["sweep_order"] == "lexicographic" ||
        error("Expected the retained lexicographic FGS configuration")
    colored["sweep_order"] == "colored" || error("Expected a colored configuration")
    Set(keys(lex)) == Set(keys(colored)) || error("A/B configs must share fields")
    for key in keys(lex)
        key in ("sweep_order", "tolerance") && continue
        lex[key] == colored[key] || error("A/B configs differ in $key")
    end
    colored["tolerance"] > 0 || error("Colored configuration requires calibrated tolerance")
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
    ab_check_pair(lex, colored)
    cold_check_config(colored; selected=true)
    cold_write_toml(joinpath(out, "colored_selected.toml"), colored)
    records = NamedTuple[]
    solver, reference = ab_prepare(colored, "colored", records)
    confirm, x = cold_trial(solver)
    agreement = norm(x-reference)/max(norm(reference), eps())
    cold_require(confirm.eligible && confirm.authoritative_evaluator == "certified_fmm" &&
        agreement <= 1e-8, "Colored calibration confirmation failed the repeat gate")
    push!(records, (; sweep_order="colored", phase="confirmation", confirm...))
    cold_csv(joinpath(out, "calibration_validation.csv"), records)
    cold_write_toml(joinpath(out, "calibration_summary.toml"), Dict(
        "colored_tolerance" => colored["tolerance"],
        "lexicographic_tolerance" => lex["tolerance"],
        "colored_iterations" => confirm.iterations,
        "confirmation_repeat_delta" => agreement))
    return nothing
end

function ab_trials(lex, colored, out)
    warmups = NamedTuple[]
    arms = Tuple{String,Any,Vector{Float64}}[]
    for (label, c) in (("lexicographic", lex), ("colored", colored))
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
        "trials_per_order" => batches*reps)
    for (label, _, _) in arms
        times = [r.solve_seconds for r in rows if r.sweep_order == label]
        iters = unique(r.iterations for r in rows if r.sweep_order == label)
        summary["$(label)_median_solve_seconds"] = median(times)
        summary["$(label)_min_solve_seconds"] = minimum(times)
        summary["$(label)_iterations"] = iters
    end
    cold_write_toml(joinpath(out, "ab_summary.toml"), summary)
    return nothing
end

function ab_activity(lex, colored, out)
    warmups = NamedTuple[]
    validations = NamedTuple[]
    for (label, c) in (("lexicographic", lex), ("colored", colored))
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
    Base.invokelatest(run_colored_ab, configs, out, stage)
end

function run_colored_ab(configs, out, stage)
    mode = get(ENV, "AB_MODE", "")
    mode in ("calibrate", "trials", "activity") || error("AB_MODE must be calibrate, trials, or activity")
    length(configs) == 1 || error("Colored A/B requires exactly one retained configuration")
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
        colored = cold_check_config(TOML.parsefile(colored_file); selected=true)
        ab_check_pair(lex, colored)
        cp(colored_file, joinpath(out, "input_colored_config.toml"))
        cold_write_toml(joinpath(out, "config_colored.toml"), colored)
        mode == "trials" ? ab_trials(lex, colored, out) : ab_activity(lex, colored, out)
    end
    cold_write_toml(joinpath(out, "status.toml"), Dict("status"=>"completed"))
    return nothing
end

abspath(PROGRAM_FILE) == (@__FILE__) && ab_main()
