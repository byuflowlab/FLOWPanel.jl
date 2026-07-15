using Test
import FLOWPanel as pnl

if !isdefined(@__MODULE__, :run_pitching_wing_convergence)
    include(joinpath(pnl.examples_path, "pitching_wing_convergence.jl"))
end

@testset "Pitching-wing convergence signal analysis" begin
    frequency = 4.01
    period = inv(frequency)
    dt = period / 160
    time = collect(0.0:dt:3period)
    amplitude = 0.27
    CL = @. 0.31 + 0.025 * time + amplitude * sin(2pi * frequency * time + 0.42)

    metrics = _pitchconv_cycle_metrics(time, CL, period, 3)
    @test metrics.amplitude ≈ amplitude rtol=2e-4
    @test metrics.frequency ≈ frequency rtol=2e-4
    @test metrics.nsamples >= 159
    @test metrics.rms_residual < 1e-10

    @test _pitchconv_relative_change(1.01, 1.0) ≈ 0.01 / 1.01
    @test _pitchconv_relative_change(0.0, 0.0) == 0.0

    path = tempname() * ".csv"
    open(path, "w") do io
        println(io, "time,t_over_T,alpha_deg,CL")
        for i in eachindex(time)
            println(io, "$(time[i]),$(time[i] / period),0,$(CL[i])")
        end
    end
    analysis = _pitchconv_analyze_history(path, period, 3;
        amplitude_tolerance=0.01, frequency_tolerance=0.005)
    @test analysis.converged
    @test analysis.amplitude_change < 1e-6
    @test analysis.frequency_change < 1e-6
end

@testset "Pitching-wing convergence state and checkpoints" begin
    root = mktempdir()
    state = _pitchconv_default_state()
    options = (;
        wake_values=[1.0, 2.0],
        dt_values=[0.25, 0.125],
        baseline_c_per_dt=0.5,
        initial_cycles=3,
        max_cycles=10,
        amplitude_tolerance=0.01,
        frequency_tolerance=0.005,
        safety_factor=1.5,
        reserve_seconds=900.0,
    )
    _pitchconv_initialize_options!(state, options)
    state_path = _pitchconv_save_state(root, state)
    @test isfile(state_path)
    loaded = _pitchconv_load_state(root)
    @test loaded["wake_values"] == options.wake_values
    @test loaded["dt_values"] == options.dt_values
    @test_throws ErrorException _pitchconv_initialize_options!(loaded,
        merge(options, (; max_cycles=9)))

    case_path = _pitchconv_case_path(root, 2.0, 0.5)
    mkpath(case_path)
    for file in _pitchconv_restart_files(case_path)
        write(file, "original")
    end
    checkpoint = _pitchconv_checkpoint!(case_path)
    @test isdir(checkpoint)
    write(first(_pitchconv_restart_files(case_path)), "damaged")
    @test _pitchconv_restore_checkpoint!(case_path)
    @test read(first(_pitchconv_restart_files(case_path)), String) == "original"
    _pitchconv_clear_checkpoint!(case_path)
    @test !isdir(checkpoint)

    loaded["records"] = Any[Dict{String,Any}(
        "elapsed_seconds" => 120.0,
        "cycles_added" => 3,
        "c_per_dt" => 0.5,
        "wake_length_spans" => 1.0,
    )]
    estimate = _pitchconv_estimated_seconds(loaded, 4.0, 0.25, 1)
    @test estimate ≈ 160.0
end

