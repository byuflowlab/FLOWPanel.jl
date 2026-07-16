using Test

if !isdefined(Main, :run_pitching_wing_pressure_comparison)
    include(joinpath(@__DIR__, "..", "examples", "pitching_wing_pressure_comparison.jl"))
end

@testset "pitching-wing pressure comparison helpers" begin
    t = [0.0, 0.1, 0.2]
    forces = Dict{Symbol,Any}(:bernoulli =>
        [0.0 0.0 0.0; 0.0 0.0 0.0; 1.0 2.0 3.0])
    for (method, _) in PRESSURE_COMPARISON_METHODS[2:end]
        forces[method] = [0.0 0.0 0.0; 0.0 0.0 0.0; 2.0 100.0 6.0]
    end
    metrics = Dict(method => (; l2=[0.1, 0.2, 0.3], linf=[0.2, 0.3, 0.4])
                   for (method, _) in PRESSURE_COMPARISON_METHODS[2:end])
    diagnostics = Dict(method => (;
        converged=[true, false, true],
        iterations=[4, 9, 5],
        absolute_residual=[1.0e-8, 2.0, 3.0e-8],
        relative_residual=[1.0e-9, 0.2, 3.0e-9],
        rhs_l2=[10.0, 20.0, 30.0],
        gradient_l2=[100.0, 200.0, 300.0],
    ) for (method, _) in PRESSURE_COMPARISON_METHODS[2:end])
    rows = _comparison_summary(t, 0.1, forces, metrics, diagnostics;
        skip_first_cycle=false)
    @test length(rows) == 4
    @test rows[1].method == :bernoulli
    @test rows[1].CL_mean ≈ 2.0
    for row in rows[2:end]
        @test row.converged_samples == 2
        @test row.total_samples == 3
        @test row.CL_mean ≈ 4.0
        @test row.CL_peak_to_peak ≈ 4.0
        @test row.relative_CL_L2 ≈ 1.0
        @test row.pressure_L2_mean ≈ 0.2
        @test row.pressure_Linf_max ≈ 0.4
    end

    failed_diagnostics = Dict(method => merge(diagnostics[method],
        (; converged=falses(length(t))))
        for (method, _) in PRESSURE_COMPARISON_METHODS[2:end])
    failed_rows = _comparison_summary(t, 0.1, forces, metrics, failed_diagnostics;
        skip_first_cycle=false)
    for row in failed_rows[2:end]
        @test row.converged_samples == 0
        @test row.total_samples == 3
        @test isnan(row.CL_mean)
        @test isnan(row.CL_peak_to_peak)
        @test isnan(row.relative_CL_L2)
        @test isnan(row.pressure_L2_mean)
        @test isnan(row.pressure_Linf_max)
    end

    mktempdir() do dir
        csv = _comparison_csv(joinpath(dir, "history.csv"), t, 0.1,
            [1.0, 2.0, 3.0], forces, metrics, diagnostics)
        summary = _write_comparison_summary(joinpath(dir, "summary.csv"), rows)
        @test countlines(csv) == 4
        @test countlines(summary) == 5
        lines = readlines(csv)
        header = split(lines[1], ',')
        @test header[1:3] == ["time", "t_over_T", "alpha_deg"]
        for (method, _) in PRESSURE_COMPARISON_METHODS[2:end]
            for field in ("cg_solved", "absolute_residual", "relative_residual",
                          "rhs_l2", "gradient_l2")
                @test "$(field)_$(method)" in header
            end
        end

        failed_sample = split(lines[3], ',')
        @test parse(Float64, failed_sample[findfirst(==("CL_edge_difference"), header)]) == 100.0
        @test failed_sample[findfirst(==("cg_solved_edge_difference"), header)] == "false"
        @test parse(Float64,
            failed_sample[findfirst(==("absolute_residual_edge_difference"), header)]) == 2.0
    end

    @testset "two-method replay report" begin
        methods = PRESSURE_COMPARISON_METHODS[1:2]
        rows = _comparison_summary(t, 0.1, forces, metrics, diagnostics;
            skip_first_cycle=false, methods)
        @test getfield.(rows, :method) == [:bernoulli, :edge_difference]

        mktempdir() do dir
            csv = _comparison_csv(joinpath(dir, "replay.csv"), t, 0.1,
                [1.0, 2.0, 3.0], forces, metrics, diagnostics; methods)
            summary = _write_comparison_summary(joinpath(dir, "summary.csv"), rows)
            header = split(first(readlines(csv)), ',')
            @test "CL_bernoulli" in header
            @test "CL_edge_difference" in header
            @test "cg_iterations_edge_difference" in header
            @test !any(occursin("corrected_hessian"), header)
            @test countlines(summary) == 3
        end
    end
end
