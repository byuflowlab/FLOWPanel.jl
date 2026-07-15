#=##############################################################################
# Adaptive convergence driver for examples/pitching_wing.jl.
#
# The driver is intentionally resumable.  Run directories are stable, completed
# cases are recorded in TOML/CSV, and cycle extensions use pitching_wing.jl's
# warm-start path instead of recomputing earlier cycles.
=###############################################################################

include(joinpath(@__DIR__, "pitching_wing.jl"))

import LinearAlgebra: norm
import Printf: @printf, @sprintf
import TOML

const PITCHCONV_FORMAT_VERSION = 1
const PITCHCONV_DEFAULT_ROOT = joinpath("data", "pitching_wing_convergence")

_pitchconv_env_float(name, default) = parse(Float64, get(ENV, name, string(default)))
_pitchconv_env_int(name, default) = parse(Int, get(ENV, name, string(default)))
_pitchconv_env_bool(name, default) = lowercase(get(ENV, name, string(default))) in
    ("1", "true", "yes", "on")

function _pitchconv_env_list(name, default)
    raw = strip(get(ENV, name, join(default, ",")))
    isempty(raw) && return Float64[]
    return parse.(Float64, strip.(split(raw, ',')))
end

function _pitchconv_relative_change(a::Real, b::Real)
    return abs(float(a) - float(b)) / max(abs(float(a)), abs(float(b)), eps(Float64))
end

function _pitchconv_harmonic_fit(time, values, frequency)
    t0 = (first(time) + last(time)) / 2
    scale = max(last(time) - first(time), eps(Float64))
    tau = (time .- t0) ./ scale
    phase = 2pi * frequency .* time
    design = hcat(ones(length(time)), tau, sin.(phase), cos.(phase))
    coefficients = design \ values
    residual = values - design * coefficients
    return (; sse=sum(abs2, residual), coefficients)
end

"""Fit the fundamental amplitude and frequency within one pitching cycle."""
function _pitchconv_cycle_metrics(time, values, period, cycle;
        prescribed_frequency=inv(period), frequency_band=(0.75, 1.25), grid_points=161)
    length(time) == length(values) || throw(ArgumentError("time and values must have equal length"))
    lo = (cycle - 1) * period
    hi = cycle * period
    tolerance = 100eps(Float64) * max(abs(hi), 1.0)
    indices = findall(t -> t >= lo - tolerance && t < hi - tolerance, time)
    length(indices) >= 8 || throw(ArgumentError(
        "cycle $(cycle) contains only $(length(indices)) samples; at least 8 are required"))
    t = collect(float.(time[indices]))
    y = collect(float.(values[indices]))
    all(isfinite, t) && all(isfinite, y) || throw(ArgumentError("cycle data must be finite"))

    fmin, fmax = prescribed_frequency .* frequency_band
    frequencies = collect(range(fmin, fmax, length=grid_points))
    errors = map(f -> _pitchconv_harmonic_fit(t, y, f).sse, frequencies)
    best = argmin(errors)
    left = frequencies[max(1, best - 1)]
    right = frequencies[min(length(frequencies), best + 1)]

    # Golden-section refinement of the variable-frequency least-squares fit.
    phi = (sqrt(5.0) - 1) / 2
    x1 = right - phi * (right - left)
    x2 = left + phi * (right - left)
    e1 = _pitchconv_harmonic_fit(t, y, x1).sse
    e2 = _pitchconv_harmonic_fit(t, y, x2).sse
    for _ in 1:48
        if e1 <= e2
            right, x2, e2 = x2, x1, e1
            x1 = right - phi * (right - left)
            e1 = _pitchconv_harmonic_fit(t, y, x1).sse
        else
            left, x1, e1 = x1, x2, e2
            x2 = left + phi * (right - left)
            e2 = _pitchconv_harmonic_fit(t, y, x2).sse
        end
    end
    frequency = (left + right) / 2
    fit = _pitchconv_harmonic_fit(t, y, frequency)
    amplitude = hypot(fit.coefficients[3], fit.coefficients[4])
    rms_residual = sqrt(fit.sse / length(t))
    return (; cycle=Int(cycle), amplitude, frequency, rms_residual,
        nsamples=length(t), mean=fit.coefficients[1])
end

function _pitchconv_read_history(path)
    lines = readlines(path)
    length(lines) >= 2 || throw(ArgumentError("pitching-wing history is empty: $(path)"))
    header = split(strip(lines[1]), ',')
    itime = findfirst(==("time"), header)
    iCL = findfirst(==("CL"), header)
    isnothing(itime) && throw(ArgumentError("history has no time column: $(path)"))
    isnothing(iCL) && throw(ArgumentError("history has no CL column: $(path)"))
    time = Float64[]
    CL = Float64[]
    for (irow, line) in enumerate(lines[2:end])
        iline = irow + 1
        isempty(strip(line)) && continue
        fields = split(strip(line), ',')
        length(fields) == length(header) || throw(ArgumentError("malformed history row $(iline)"))
        push!(time, parse(Float64, fields[itime]))
        push!(CL, parse(Float64, fields[iCL]))
    end
    return (; time, CL)
end

function _pitchconv_analyze_history(path, period, n_cycles; amplitude_tolerance, frequency_tolerance)
    history = _pitchconv_read_history(path)
    n_cycles >= 2 || throw(ArgumentError("at least two cycles are required for convergence analysis"))
    previous = _pitchconv_cycle_metrics(history.time, history.CL, period, n_cycles - 1)
    last = _pitchconv_cycle_metrics(history.time, history.CL, period, n_cycles)
    amplitude_change = _pitchconv_relative_change(last.amplitude, previous.amplitude)
    frequency_change = _pitchconv_relative_change(last.frequency, previous.frequency)
    return (; previous, last, amplitude_change, frequency_change,
        converged=amplitude_change <= amplitude_tolerance && frequency_change <= frequency_tolerance)
end

function _pitchconv_default_state()
    return Dict{String,Any}(
        "format_version" => PITCHCONV_FORMAT_VERSION,
        "status" => "ready",
        "stage" => "wake",
        "wake_index" => 1,
        "dt_index" => 1,
        "active_cycles" => 3,
        "accepted_wake" => -1.0,
        "accepted_c_per_dt" => -1.0,
        "accepted_cycles" => -1,
        "reference_amplitude" => -1.0,
        "reference_frequency" => -1.0,
        "reference_path" => "",
        "running_requested_cycles" => -1,
        "running_prior_cycles" => -1,
        "records" => Any[],
    )
end

function _pitchconv_state_path(root)
    return joinpath(root, "pitching_wing_convergence.state.toml")
end

function _pitchconv_load_state(root)
    path = _pitchconv_state_path(root)
    state = isfile(path) ? TOML.parsefile(path) : _pitchconv_default_state()
    get(state, "format_version", 0) == PITCHCONV_FORMAT_VERSION || error(
        "Unsupported pitching-wing convergence state in $(path)")
    return state
end

function _pitchconv_save_state(root, state)
    mkpath(root)
    path = _pitchconv_state_path(root)
    temporary = path * ".tmp"
    open(temporary, "w") do io
        TOML.print(io, state; sorted=true)
    end
    mv(temporary, path; force=true)
    return path
end

function _pitchconv_case_name(wake, c_per_dt)
    encode(x) = replace(@sprintf("%.8g", x), "." => "p", "-" => "m")
    return "wake_$(encode(wake))__c_per_dt_$(encode(c_per_dt))"
end

_pitchconv_case_path(root, wake, c_per_dt) = joinpath(root, "cases",
    _pitchconv_case_name(wake, c_per_dt))

function _pitchconv_restart_files(path)
    names = (
        PITCHING_WING_NAME * ".metadata.toml",
        PITCHING_WING_NAME * "_body1.pvd",
        PITCHING_WING_NAME * "_wake1.pvd",
        PITCHING_WING_NAME * "_unsteady.csv",
        PITCHING_WING_CONFIG_FILE,
        PITCHING_WING_PRESSURE_STATE_FILE,
    )
    return [joinpath(path, name) for name in names]
end

function _pitchconv_checkpoint!(path)
    checkpoint = joinpath(path, ".restart_checkpoint")
    isdir(checkpoint) && rm(checkpoint; recursive=true, force=true)
    mkpath(checkpoint)
    for source in _pitchconv_restart_files(path)
        isfile(source) && cp(source, joinpath(checkpoint, basename(source)); force=true)
    end
    return checkpoint
end

function _pitchconv_restore_checkpoint!(path)
    checkpoint = joinpath(path, ".restart_checkpoint")
    isdir(checkpoint) || return false
    for source in readdir(checkpoint; join=true)
        cp(source, joinpath(path, basename(source)); force=true)
    end
    return true
end

function _pitchconv_clear_checkpoint!(path)
    checkpoint = joinpath(path, ".restart_checkpoint")
    isdir(checkpoint) && rm(checkpoint; recursive=true, force=true)
    return nothing
end

function _pitchconv_prune_case!(path)
    isempty(path) && return
    isdir(path) || return
    for suffix in ("_body1", "_wake1")
        target = joinpath(path, PITCHING_WING_NAME * suffix)
        isdir(target) && rm(target; recursive=true, force=true)
    end
    for file in (
            PITCHING_WING_NAME * ".metadata.toml",
            PITCHING_WING_NAME * "_body1.pvd",
            PITCHING_WING_NAME * "_wake1.pvd",
            PITCHING_WING_PRESSURE_STATE_FILE)
        rm(joinpath(path, file); force=true)
    end
    _pitchconv_clear_checkpoint!(path)
    open(joinpath(path, "VTK_PRUNED.txt"), "w") do io
        println(io, "Large restart/VTK artifacts were removed after this case was rejected.")
        println(io, "The unsteady CSV and configuration remain available for analysis.")
    end
    return nothing
end

function _pitchconv_write_summary(root, records)
    path = joinpath(root, "pitching_wing_convergence_summary.csv")
    columns = ["stage", "wake_length_spans", "c_per_dt", "n_cycles", "cycles_added",
        "previous_amplitude", "last_amplitude", "cycle_amplitude_change",
        "previous_frequency_hz", "last_frequency_hz", "cycle_frequency_change",
        "cycle_converged", "reference_amplitude_change", "reference_frequency_change",
        "reference_converged", "elapsed_seconds", "case_path"]
    open(path, "w") do io
        println(io, join(columns, ','))
        for record in records
            println(io, join((get(record, column, "") for column in columns), ','))
        end
    end
    return path
end

function _pitchconv_record!(root, state, record)
    push!(state["records"], record)
    _pitchconv_write_summary(root, state["records"])
    _pitchconv_save_state(root, state)
    return record
end

function _pitchconv_seconds_remaining()
    raw = get(ENV, "SLURM_JOB_END_TIME", "")
    isempty(raw) && return Inf
    return parse(Int, raw) - time()
end

function _pitchconv_estimated_seconds(state, wake, c_per_dt, cycles_added)
    estimates = Float64[]
    for record in state["records"]
        elapsed = float(get(record, "elapsed_seconds", 0.0))
        added = Int(get(record, "cycles_added", 0))
        elapsed > 0 && added > 0 || continue
        dt_scale = float(record["c_per_dt"]) / c_per_dt
        wake_scale = sqrt(max(1.0, wake / float(record["wake_length_spans"])))
        push!(estimates, elapsed / added * cycles_added * dt_scale * wake_scale)
    end
    return isempty(estimates) ? nothing : maximum(estimates)
end

function _pitchconv_has_time(state, wake, c_per_dt, cycles_added; safety_factor, reserve_seconds)
    estimate = _pitchconv_estimated_seconds(state, wake, c_per_dt, cycles_added)
    isnothing(estimate) && return (true, estimate, _pitchconv_seconds_remaining())
    remaining = _pitchconv_seconds_remaining()
    return (remaining >= safety_factor * estimate + reserve_seconds, estimate, remaining)
end

function _pitchconv_run_case!(root, state, wake, c_per_dt, n_cycles, prior_cycles, options)
    path = _pitchconv_case_path(root, wake, c_per_dt)
    history_path = joinpath(path, PITCHING_WING_NAME * "_unsteady.csv")
    extending = prior_cycles > 0
    cycles_added = extending ? n_cycles - prior_cycles : n_cycles

    enough, estimate, remaining = _pitchconv_has_time(state, wake, c_per_dt, cycles_added;
        safety_factor=options.safety_factor, reserve_seconds=options.reserve_seconds)
    if !enough
        state["status"] = "needs_resubmit"
        _pitchconv_save_state(root, state)
        @printf("Not starting the next case: estimated %.1f h, %.1f h remain.\n",
            estimate / 3600, remaining / 3600)
        println("Resubmit the same Slurm script; completed work will be reused.")
        return nothing
    end

    if extending
        _pitchconv_checkpoint!(path)
    elseif isdir(path)
        rm(path; recursive=true, force=true)
    end
    state["status"] = "running"
    state["running_requested_cycles"] = n_cycles
    state["running_prior_cycles"] = prior_cycles
    _pitchconv_save_state(root, state)

    println("\nRunning wake_length_spans=$(wake), c_per_dt=$(c_per_dt), n_cycles=$(n_cycles)")
    started = time()
    result = run_pitching_wing(;
        path,
        wake_length_spans=wake,
        c_per_dt,
        n_cycles,
        include_static_polar=false,
        plot_convergence=false,
        save_vtk=true,
        check_existing=extending,
        restart_input=IOBuffer(extending ? "y\n" : ""),
    )
    elapsed = time() - started
    isfile(history_path) || error("Simulation completed without writing $(history_path)")
    _pitchconv_clear_checkpoint!(path)
    state["status"] = "ready"
    state["running_requested_cycles"] = -1
    state["running_prior_cycles"] = -1
    _pitchconv_save_state(root, state)
    return (; path, history_path, elapsed, cycles_added, period=result.setup.period)
end

function _pitchconv_recover_interrupted!(root, state)
    state["status"] == "running" || return state
    stage = state["stage"]
    wake = stage == "wake" ? state["wake_values"][state["wake_index"]] : state["accepted_wake"]
    c_per_dt = stage == "wake" ? state["baseline_c_per_dt"] : state["dt_values"][state["dt_index"]]
    path = _pitchconv_case_path(root, wake, c_per_dt)
    prior = state["running_prior_cycles"]
    if prior > 0
        _pitchconv_restore_checkpoint!(path) || error(
            "Interrupted extension has no restart checkpoint in $(path)")
        println("Restored the completed $(prior)-cycle checkpoint after an interrupted extension.")
    else
        isdir(path) && rm(path; recursive=true, force=true)
        println("Removed an incomplete initial case; it will be rerun.")
    end
    state["status"] = "ready"
    state["running_requested_cycles"] = -1
    state["running_prior_cycles"] = -1
    _pitchconv_save_state(root, state)
    return state
end

function _pitchconv_initialize_options!(state, options)
    stored = Dict(
        "wake_values" => options.wake_values,
        "dt_values" => options.dt_values,
        "baseline_c_per_dt" => options.baseline_c_per_dt,
        "initial_cycles" => options.initial_cycles,
        "max_cycles" => options.max_cycles,
        "amplitude_tolerance" => options.amplitude_tolerance,
        "frequency_tolerance" => options.frequency_tolerance,
    )
    if !haskey(state, "wake_values")
        merge!(state, stored)
    else
        mismatches = [key for (key, value) in stored if state[key] != value]
        isempty(mismatches) || error("Existing convergence state uses different settings: " * join(mismatches, ", "))
    end
    return state
end

function _pitchconv_comparison(last, state, options)
    has_reference = state["reference_amplitude"] >= 0
    if !has_reference
        return (; has_reference=false, amplitude_change=-1.0, frequency_change=-1.0,
            converged=false)
    end
    amplitude_change = _pitchconv_relative_change(last.amplitude, state["reference_amplitude"])
    frequency_change = _pitchconv_relative_change(last.frequency, state["reference_frequency"])
    return (; has_reference=true, amplitude_change, frequency_change,
        converged=amplitude_change <= options.amplitude_tolerance &&
            frequency_change <= options.frequency_tolerance)
end

function _pitchconv_final_report(root, state)
    path = joinpath(root, "pitching_wing_convergence_report.txt")
    open(path, "w") do io
        println(io, "FLOWPanel pitching-wing convergence")
        println(io, "status = ", state["status"])
        println(io, "wake_length_spans = ", state["accepted_wake"])
        println(io, "n_cycles = ", state["accepted_cycles"])
        println(io, "c_per_dt = ", state["accepted_c_per_dt"])
        println(io, "amplitude_relative_tolerance = ", state["amplitude_tolerance"])
        println(io, "frequency_relative_tolerance = ", state["frequency_tolerance"])
    end
    return path
end

function run_pitching_wing_convergence(;
        root=get(ENV, "PITCHCONV_ROOT", PITCHCONV_DEFAULT_ROOT),
        wake_values=_pitchconv_env_list("PITCHCONV_WAKE_VALUES", [1.0, 2.0, 4.0, 8.0, 16.0]),
        dt_values=_pitchconv_env_list("PITCHCONV_DT_VALUES", [0.25, 0.125, 0.0625, 0.03125]),
        baseline_c_per_dt=_pitchconv_env_float("PITCHCONV_BASELINE_C_PER_DT", 0.5),
        initial_cycles=_pitchconv_env_int("PITCHCONV_INITIAL_CYCLES", 3),
        max_cycles=_pitchconv_env_int("PITCHCONV_MAX_CYCLES", 10),
        amplitude_tolerance=_pitchconv_env_float("PITCHCONV_AMPLITUDE_TOL", 0.01),
        frequency_tolerance=_pitchconv_env_float("PITCHCONV_FREQUENCY_TOL", 0.005),
        safety_factor=_pitchconv_env_float("PITCHCONV_TIME_SAFETY", 1.5),
        reserve_seconds=_pitchconv_env_float("PITCHCONV_TIME_RESERVE_SECONDS", 900.0),
        dry_run=_pitchconv_env_bool("PITCHCONV_DRY_RUN", false))
    isempty(wake_values) && throw(ArgumentError("wake_values cannot be empty"))
    isempty(dt_values) && throw(ArgumentError("dt_values cannot be empty"))
    initial_cycles >= 2 || throw(ArgumentError("initial_cycles must be at least 2"))
    max_cycles >= initial_cycles || throw(ArgumentError("max_cycles must be >= initial_cycles"))
    options = (; wake_values=collect(float.(wake_values)), dt_values=collect(float.(dt_values)),
        baseline_c_per_dt=float(baseline_c_per_dt), initial_cycles, max_cycles,
        amplitude_tolerance=float(amplitude_tolerance), frequency_tolerance=float(frequency_tolerance),
        safety_factor=float(safety_factor), reserve_seconds=float(reserve_seconds))

    if dry_run
        println("Pitching-wing convergence dry run")
        println("  wake stage: c_per_dt=$(baseline_c_per_dt), wake_length_spans=$(wake_values), cycles=$(initial_cycles):$(max_cycles)")
        println("  timestep stage: c_per_dt=$(dt_values), cycles up to $(max_cycles)")
        println("  tolerances: amplitude=$(amplitude_tolerance), frequency=$(frequency_tolerance)")
        return nothing
    end

    state = _pitchconv_load_state(root)
    _pitchconv_initialize_options!(state, options)
    state["status"] == "needs_resubmit" && (state["status"] = "ready")
    _pitchconv_recover_interrupted!(root, state)
    _pitchconv_save_state(root, state)

    while state["stage"] in ("wake", "timestep")
        stage = state["stage"]
        wake = stage == "wake" ? options.wake_values[state["wake_index"]] : state["accepted_wake"]
        c_per_dt = stage == "wake" ? options.baseline_c_per_dt : options.dt_values[state["dt_index"]]
        n_cycles = state["active_cycles"]
        path = _pitchconv_case_path(root, wake, c_per_dt)
        prior_cycles = isfile(joinpath(path, PITCHING_WING_NAME * "_unsteady.csv")) ? n_cycles - 1 : 0
        run = _pitchconv_run_case!(root, state, wake, c_per_dt, n_cycles, prior_cycles, options)
        isnothing(run) && return state

        analysis = _pitchconv_analyze_history(run.history_path, run.period, n_cycles;
            amplitude_tolerance, frequency_tolerance)
        comparison = analysis.converged ? _pitchconv_comparison(analysis.last, state, options) :
            (; has_reference=false, amplitude_change=-1.0, frequency_change=-1.0, converged=false)
        record = Dict{String,Any}(
            "stage" => stage,
            "wake_length_spans" => wake,
            "c_per_dt" => c_per_dt,
            "n_cycles" => n_cycles,
            "cycles_added" => run.cycles_added,
            "previous_amplitude" => analysis.previous.amplitude,
            "last_amplitude" => analysis.last.amplitude,
            "cycle_amplitude_change" => analysis.amplitude_change,
            "previous_frequency_hz" => analysis.previous.frequency,
            "last_frequency_hz" => analysis.last.frequency,
            "cycle_frequency_change" => analysis.frequency_change,
            "cycle_converged" => analysis.converged,
            "reference_amplitude_change" => comparison.amplitude_change,
            "reference_frequency_change" => comparison.frequency_change,
            "reference_converged" => comparison.converged,
            "elapsed_seconds" => run.elapsed,
            "case_path" => run.path,
        )
        _pitchconv_record!(root, state, record)
        @printf("Cycle comparison: amplitude %.3g%%, frequency %.3g%% => %s\n",
            100analysis.amplitude_change, 100analysis.frequency_change,
            analysis.converged ? "converged" : "extend")

        if !analysis.converged
            if n_cycles >= options.max_cycles
                state["stage"] = "failed"
                state["status"] = "inconclusive_cycle_convergence"
                _pitchconv_save_state(root, state)
                break
            end
            state["active_cycles"] = n_cycles + 1
            _pitchconv_save_state(root, state)
            continue
        end

        if !comparison.has_reference
            state["reference_amplitude"] = analysis.last.amplitude
            state["reference_frequency"] = analysis.last.frequency
            state["reference_path"] = run.path
            state["wake_index"] += 1
            state["active_cycles"] = options.initial_cycles
            if state["wake_index"] > length(options.wake_values)
                state["stage"] = "failed"
                state["status"] = "inconclusive_wake_convergence"
                _pitchconv_save_state(root, state)
                break
            end
            _pitchconv_save_state(root, state)
            continue
        end

        @printf("Refinement comparison: amplitude %.3g%%, frequency %.3g%% => %s\n",
            100comparison.amplitude_change, 100comparison.frequency_change,
            comparison.converged ? "converged" : "refine")
        _pitchconv_prune_case!(state["reference_path"])

        if stage == "wake" && comparison.converged
            state["accepted_wake"] = wake
            state["accepted_cycles"] = n_cycles
            state["reference_amplitude"] = analysis.last.amplitude
            state["reference_frequency"] = analysis.last.frequency
            state["reference_path"] = run.path
            state["stage"] = "timestep"
            state["dt_index"] = 1
            state["active_cycles"] = n_cycles
        elseif stage == "timestep" && comparison.converged
            state["accepted_c_per_dt"] = c_per_dt
            state["accepted_cycles"] = n_cycles
            state["stage"] = "complete"
            state["status"] = "converged"
        else
            state["reference_amplitude"] = analysis.last.amplitude
            state["reference_frequency"] = analysis.last.frequency
            state["reference_path"] = run.path
            if stage == "wake"
                state["wake_index"] += 1
                state["active_cycles"] = options.initial_cycles
                if state["wake_index"] > length(options.wake_values)
                    state["stage"] = "failed"
                    state["status"] = "inconclusive_wake_convergence"
                end
            else
                state["dt_index"] += 1
                state["active_cycles"] = n_cycles
                if state["dt_index"] > length(options.dt_values)
                    state["stage"] = "failed"
                    state["status"] = "inconclusive_timestep_convergence"
                end
            end
        end
        _pitchconv_save_state(root, state)
    end

    report = _pitchconv_final_report(root, state)
    println("\nPitching-wing convergence status: $(state["status"])")
    println("Report: $(report)")
    println("Summary: $(_pitchconv_write_summary(root, state["records"]))")
    return state
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_pitching_wing_convergence()
end
