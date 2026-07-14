using Test

if !isdefined(Main, :run_pitching_wing_pressure_comparison)
    include(joinpath(@__DIR__, "..", "examples", "pitching_wing_pressure_comparison.jl"))
end

@testset "pitching-wing pressure comparison helpers" begin
    t = [0.0, 0.1, 0.2]
    forces = Dict{Symbol,Any}(method => [0.0 0.0 0.0; 0.0 0.0 0.0; 1.0 1.1 1.2]
                              for (method, _) in PRESSURE_COMPARISON_METHODS)
    metrics = Dict(method => (; l2=[0.1, 0.2, 0.3], linf=[0.2, 0.3, 0.4])
                   for (method, _) in PRESSURE_COMPARISON_METHODS[2:end])
    rows = _comparison_summary(t, 0.1, forces, metrics; skip_first_cycle=false)
    @test length(rows) == 4
    @test rows[1].method == :bernoulli
    @test rows[1].CL_mean ≈ 1.1
    @test rows[2].pressure_L2_mean ≈ 0.2

    mktempdir() do dir
        csv = _comparison_csv(joinpath(dir, "history.csv"), t, 0.1,
            [1.0, 2.0, 3.0], forces, metrics)
        summary = _write_comparison_summary(joinpath(dir, "summary.csv"), rows)
        @test countlines(csv) == 4
        @test countlines(summary) == 5
        @test startswith(readline(csv), "time,t_over_T,alpha_deg")
    end
end
