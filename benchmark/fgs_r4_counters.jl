#!/usr/bin/env julia
# Separate diagnostics only: these observations are never performance trials.
include(joinpath(@__DIR__, "fgs_cold_common.jl"))

function activity_snapshot()
    result = Dict{Int,Tuple{Int,Int}}()
    for path in readdir("/proc/self/task"; join=true)
        tid = tryparse(Int, basename(path))
        tid === nothing && continue
        # comm is parenthesized and can contain spaces or parentheses.
        stat = try
            read(joinpath(path, "stat"), String)
        catch err
            isdir(path) && rethrow(err)
            continue
        end
        fields = split(stat[findlast(')', stat)+2:end])
        result[tid] = (parse(Int, fields[12]) + parse(Int, fields[13]),
                       parse(Int, fields[37]))
    end
    return result
end

function activity_observer(rows)
    active = Ref{Any}(nothing)
    sequence = Ref(0)
    return function(stage, event)
        if event == :start
            active[] === nothing || error("Nested activity stage")
            sequence[] += 1
            before = activity_snapshot()
            active[] = (stage, time_ns(), before)
        elseif event == :stop
            stop = time_ns()
            active[] !== nothing && active[][1] == stage || error("Unmatched activity stage")
            _, start, before = active[]
            after = activity_snapshot()
            for tid in sort!(collect(union(keys(before), keys(after))))
                complete = haskey(before, tid) && haskey(after, tid)
                b = get(before, tid, (0, -1)); a = get(after, tid, (0, -1))
                push!(rows, (; sequence=sequence[], stage=string(stage), tid,
                    span_seconds=(stop-start)*1e-9,
                    cpu_ticks=complete ? a[1]-b[1] : -1,
                    complete_endpoints=complete, cpu_before=b[2], cpu_after=a[2]))
            end
            active[] = nothing
        else
            error("Unknown activity event")
        end
        nothing
    end
end

function counter_command(control, acknowledgement, command)
    println(control, command)
    flush(control)
    readline(acknowledgement) == "ack" || error("perf did not acknowledge $command")
    nothing
end

function counters_main()
    configs, out, stage = cold_initialize!()
    Base.invokelatest(run_counters, configs, out, stage)
end

function run_counters(configs, out, stage)
    length(configs) == 1 || error("Expected one retained configuration")
    c = copy(only(configs))
    c["kind"] == "fgs" && c["sweep_order"] == "lexicographic" || error("Unexpected configuration")
    cold_calibrate_selected!(c, out; selected=true)
    cold_write_toml(joinpath(out, "config.toml"), c)
    reset_cold!()
    solver = cold_make(c)
    compile_row, _ = cold_trial(solver)
    cold_require(compile_row.eligible && compile_row.authoritative_evaluator == "certified_fmm",
        "Counter compilation solve failed")

    validations = NamedTuple[]
    reference = nothing
    reference_history = Float64[]
    activity = NamedTuple[]
    # Compile observer I/O and perf handshake before enabling workload counters.
    observer = activity_observer(activity)
    observer(:probe, :start); observer(:probe, :stop); empty!(activity)
    open(ENV["COLD_PERF_CONTROL"], "r+") do control
        open(ENV["COLD_PERF_ACK"], "r+") do ack
            counter_command(control, ack, "disable")
            for mode in (:baseline, :activity, :counters)
                cold_assert_threads()
                cold_reset!(solver)
                residuals = Float64[]
                callback = (_, residual) -> push!(residuals, residual)
                observe = mode == :activity ? observer : nothing
                mode == :counters && counter_command(control, ack, "enable")
                start = time_ns()
                try
                    pnl._solve!(rotor, solver; callback, stage_observer=observe)
                finally
                    # Disable before copying the solution, certification, or I/O.
                    mode == :counters && counter_command(control, ack, "disable")
                end
                elapsed = (time_ns()-start)*1e-9
                cold_assert_threads()
                x = copy(rotor.strength[:,2])
                e = cold_validate(x; crosscheck=mode == :baseline)
                if mode == :baseline
                    reference = x
                    reference_history = residuals
                end
                delta = norm(x-reference)/max(norm(reference), eps())
                same_history = residuals == reference_history
                retained = Base.summarysize((rotor, solver))
                peak_rss = Sys.maxrss()
                cold_require(all(isfinite, x) && solver.solved && e.accepted &&
                    e.authoritative_evaluator == "certified_fmm" && delta <= 1e-8 && same_history &&
                    retained <= cold_memory_limit && peak_rss <= cold_memory_limit,
                    "$mode failed diagnostic correctness/history gates")
                push!(validations, (; mode=string(mode), diagnostic_seconds=elapsed,
                    solved=solver.solved, finite_solution=all(isfinite, x),
                    retained_bytes=retained, process_peak_rss_bytes=peak_rss,
                    iterations=solver.niter, relative_solution_delta=delta,
                    identical_history=same_history, e...))
                cold_csv(joinpath(out, "counter_validation.csv"), validations)
                cold_csv(joinpath(out, "stage_thread_activity.csv"), activity)
            end
        end
    end
    cold_write_toml(joinpath(out, "counter_scope.toml"), Dict(
        "scope" => "one zero-reset prepared solve; reset, setup and validation excluded",
        "performance_trial" => false,
        "perf_boundary_overhead" => "enable acknowledgement return and disable command handling included",
        "activity_scope" => "separate solve with coarse stage boundaries; nearfield_update groups leaf/product/scatter and updates",
        "activity_limitations" => "proc snapshots are sequential and quantized to CLK_TCK; short stages may resolve to zero; incomplete endpoints have cpu_ticks=-1",
        "bandwidth_saturation" => "not established by generic cache counters"))
    cold_write_toml(joinpath(out, "status.toml"), Dict("status"=>"completed"))
end

abspath(PROGRAM_FILE) == (@__FILE__) && counters_main()
