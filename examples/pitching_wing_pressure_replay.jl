# Replay the saved finite-panel-wake pitching-wing case through complete
# unsteady Bernoulli and corrected-edge PressureLaplace.

if !haskey(ENV, "OMP_NUM_THREADS") && haskey(ENV, "FLOWPANEL_OMP_NUM_THREADS")
    ENV["OMP_NUM_THREADS"] = ENV["FLOWPANEL_OMP_NUM_THREADS"]
end

include(joinpath(@__DIR__, "pitching_wing_pressure_comparison.jl"))

const PITCHING_WING_REPLAY_METHODS = (
    (:bernoulli, "Complete unsteady Bernoulli"),
    (:edge_difference, "Corrected-edge PressureLaplace"),
)

function _pitching_wing_replay_steps(spec::AbstractString)
    value = strip(spec)
    lowercase(value) == "all" && return :all
    fields = split(value, ':')
    if length(fields) == 1
        return parse(Int, only(fields))
    elseif length(fields) == 2
        return parse(Int, fields[1]):parse(Int, fields[2])
    end
    throw(ArgumentError("PRESSURE_REPLAY_STEPS must be all, one step, or a contiguous start:stop range; got $(repr(spec))."))
end

function _validate_pitching_wing_replay_steps(steps)
    steps === :all && return steps
    selected = steps isa Integer ? [Int(steps)] : collect(Int, steps)
    isempty(selected) && throw(ArgumentError("Pressure replay requires at least one saved step."))
    all(diff(selected) .== 1) || throw(ArgumentError(
        "Unsteady pressure replay steps must be contiguous and chronological; got $(selected)."))
    return steps
end

function _pitching_wing_replay_case(input_path, run_name)
    metadata = pnl.TOML.parsefile(joinpath(input_path, run_name * ".metadata.toml"))
    config = pnl.TOML.parsefile(joinpath(input_path, run_name * ".config.toml"))
    force_meta = only(m for m in metadata["monitor"] if get(m, "type", "") == "ForceMonitor")
    normalization = force_meta["normalization"]
    normalization["type"] == "WingNormalization" || throw(ArgumentError(
        "Saved pitching-wing force normalization is $(normalization["type"]), expected WingNormalization."))
    frequency = Float64(config["kinematics"]["frequency_hz"])
    return (;
        rho=Float64(normalization["rho"]),
        Sref=Float64(normalization["Sref"]),
        Lref=Float64(normalization["Lref"]),
        period=inv(frequency),
        alpha_mean_deg=Float64(config["kinematics"]["alpha_mean_deg"]),
        alpha_amp_deg=Float64(config["kinematics"]["alpha_amp_deg"]),
        omega=2pi * frequency,
    )
end

function run_pitching_wing_pressure_replay(;
        input_path=get(ENV, "PRESSURE_REPLAY_INPUT_PATH",
            joinpath("data", "pitching_wing_convergence", "cases", "wake_2__c_per_dt_0p25")),
        output_path=get(ENV, "PRESSURE_REPLAY_OUTPUT_PATH",
            joinpath("data", "pitching_wing_pressure_replay")),
        run_name=get(ENV, "PRESSURE_REPLAY_RUN_NAME", "pitching_wing"),
        steps=_pitching_wing_replay_steps(get(ENV, "PRESSURE_REPLAY_STEPS", "0:988")),
        plot=_env_bool("PRESSURE_REPLAY_PLOT", true),
        pressure_itmax=nothing,
        pressure_itmax_per_panel=_env_float("PRESSURE_ITMAX_PER_PANEL", 2.0),
        backend=pnl.FastMultipoleBackend(
            expansion_order=_env_int("FMM_EXPANSION_ORDER", 8),
            multipole_acceptance=_env_float("FMM_ACCEPTANCE", 0.4),
            leaf_size=_env_int("FMM_LEAF_SIZE", 40)),
        verbose=_env_bool("PRESSURE_REPLAY_VERBOSE", true),
        fail_on_nonconvergence::Bool=true)
    _validate_pitching_wing_replay_steps(steps)
    case = _pitching_wing_replay_case(input_path, run_name)
    normalization = pnl.WingNormalization(case.rho, case.Sref, case.Lref)
    report = Ref{Any}()

    monitor_factory = function(systems, wakes, frames, t_range)
        length(wakes) == 1 && wakes[1] isa pnl.PanelWake || throw(ArgumentError(
            "Pressure replay requires one reconstructed PanelWake; got $(typeof(wakes))."))
        wakes[1].include_final_filament && throw(ArgumentError(
            "Pressure replay requires include_final_filament=false for complete unsteady Bernoulli."))
        nt = length(t_range)
        itmax = isnothing(pressure_itmax) ?
            max(1000, ceil(Int, pressure_itmax_per_panel * systems[1].ncells)) :
            Int(pressure_itmax)
        bernoulli = pnl.PressureBernoulli(case.rho;
            unsteady=true, allow_partial=false, correct_kuttacondition=false,
            backend, file=false)
        edge = pnl.PressureLaplace(systems, case.rho;
            unsteady=true, gradient_mode=:edge_difference,
            acceleration_form=:material_derivative, reference_panel=1,
            itmax, verbose=false, file=false)
        pressures = Dict{Symbol,Any}(:bernoulli => bernoulli, :edge_difference => edge)
        forces = Dict{Symbol,Any}()
        recorders = Dict{Symbol,PressureComparisonRecorder}()
        reference_history = Vector{Vector{Float64}}()
        monitors = Any[]
        for (method, _) in PITCHING_WING_REPLAY_METHODS
            pressure = pressures[method]
            force = pnl.ForceMonitor(nt, 1; normalization, i_frame=-1,
                correct_kuttacondition=false, verbose=false, file=false)
            recorder = PressureComparisonRecorder(method, nt; reference_history,
                laplace=method == :bernoulli ? nothing : pressure)
            forces[method] = force.force
            recorders[method] = recorder
            append!(monitors, (pressure, force, recorder))
        end
        report[] = (; pressures, forces, recorders, itmax)
        return Tuple(monitors)
    end

    println("Pitching-wing panel-wake pressure replay")
    @printf("  input = %s\n  output = %s\n", input_path, output_path)
    @printf("  Julia threads = %d, OMP_NUM_THREADS = %s\n",
        Threads.nthreads(), get(ENV, "OMP_NUM_THREADS", "unset"))
    println("  requested steps = $(steps)")
    elapsed = @elapsed result = pnl.replay(input_path, run_name;
        monitor_factory, recompute=(:auto,), steps, backend,
        backend_wake=backend, backend_system=backend, verbose)

    state = report[]
    recorders = state.recorders
    pressure_metrics = Dict(method => (;
        l2=recorders[method].pressure_l2_relative,
        linf=recorders[method].pressure_linf)
        for (method, _) in PITCHING_WING_REPLAY_METHODS[2:end])
    diagnostics = Dict(method => (;
        converged=recorders[method].converged,
        iterations=recorders[method].iterations,
        absolute_residual=recorders[method].absolute_residual,
        relative_residual=recorders[method].relative_residual,
        rhs_l2=recorders[method].rhs_l2,
        gradient_l2=recorders[method].gradient_l2)
        for (method, _) in PITCHING_WING_REPLAY_METHODS[2:end])
    alpha_history = [case.alpha_mean_deg + case.alpha_amp_deg * sin(case.omega * t)
        for t in result.t_range]

    mkpath(output_path)
    csv_path = _comparison_csv(joinpath(output_path, "pitching_wing_pressure_replay.csv"),
        result.t_range, case.period, alpha_history, state.forces, pressure_metrics,
        diagnostics; methods=PITCHING_WING_REPLAY_METHODS)
    rows = _comparison_summary(result.t_range, case.period, state.forces,
        pressure_metrics, diagnostics; skip_first_cycle=true,
        methods=PITCHING_WING_REPLAY_METHODS)
    summary_path = _write_comparison_summary(
        joinpath(output_path, "pitching_wing_pressure_replay_summary.csv"), rows)
    plot_paths = plot ? _plot_comparison(output_path, result.t_range, case.period,
        state.forces, pressure_metrics; methods=PITCHING_WING_REPLAY_METHODS,
        stem="pitching_wing_pressure_replay") : nothing

    edge_diagnostics = diagnostics[:edge_difference]
    failed = findall(!, edge_diagnostics.converged)
    @printf("  completed %d samples in %.3f s; PressureLaplace itmax = %d\n",
        length(result.steps), elapsed, state.itmax)
    println("  history CSV: $(csv_path)")
    println("  summary CSV: $(summary_path)")
    plot && println("  plots: $(plot_paths.loads), $(plot_paths.errors)")
    if !isempty(failed)
        saved_steps = result.steps[failed]
        message = "PressureLaplace failed to converge at $(length(failed)) replay sample(s), saved steps $(saved_steps). Diagnostics were written to $(csv_path)."
        if fail_on_nonconvergence
            error(message)
        else
            @warn message
        end
    end

    return (; result, pressures=state.pressures, forces=state.forces,
        pressure_metrics, diagnostics, rows, csv_path, summary_path, plot_paths,
        elapsed, pressure_itmax=state.itmax)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_pitching_wing_pressure_replay()
end
