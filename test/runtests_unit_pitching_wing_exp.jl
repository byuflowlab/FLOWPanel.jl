using Test

if !isdefined(@__MODULE__, :PitchingWingExpCurve)
    include(joinpath(@__DIR__, "..", "data", "pitching_wing_exp", "load.jl"))
end

@testset "Pitching-wing experimental lift curves" begin
    fixture_dir = joinpath(@__DIR__, "..", "data", "pitching_wing_exp", "data")

    @test pitching_wing_exp_eta_token(0.25) == "0p250"
    @test pitching_wing_exp_filename(:dynamic, :tripped, 0.9) ==
        "dynamic_tripped_eta_0p900.csv"
    @test pitching_wing_exp_filename(:dynamic, :tripped, 0.9; quantity=:cm) ==
        "dynamic_tripped_eta_0p900_cm.csv"
    @test pitching_wing_exp_cd_filename(:dynamic, :tripped) ==
        "dynamic_tripped_eta_all_cd.csv"
    for invalid_tripping in (:notrip, "notrip")
        @test_throws ArgumentError pitching_wing_exp_filename(
            :dynamic, invalid_tripping, 0.9)
        @test_throws ArgumentError pitching_wing_exp_cd_filename(
            :dynamic, invalid_tripping)
        @test_throws ArgumentError load_pitching_wing_exp_curve(
            :dynamic, invalid_tripping, 0.25; data_dir=fixture_dir)
        @test_throws ArgumentError load_pitching_wing_exp_cd(;
            data_dir=fixture_dir, regime=:dynamic, tripping=invalid_tripping)
    end

    static = load_pitching_wing_exp_curve(:static, :untripped, 0.25; data_dir=fixture_dir)
    @test !isempty(static.aoa_deg)
    @test length(static.aoa_deg) == length(static.values)
    @test issorted(static.aoa_deg)
    @test !pitching_wing_exp_has_hysteresis(static)
    @test static.source == joinpath(fixture_dir, "static_untripped_eta_0p250.csv")

    dynamic_cl = load_pitching_wing_exp_curve(:dynamic, :untripped, 0.25; data_dir=fixture_dir)
    @test !isempty(dynamic_cl.aoa_deg)
    @test length(dynamic_cl.aoa_deg) == length(dynamic_cl.values)
    @test dynamic_cl.source ==
        joinpath(fixture_dir, "dynamic_untripped_eta_0p250.csv")
    @test pitching_wing_exp_has_hysteresis(dynamic_cl)
    @test length(pitching_wing_exp_sweeps(dynamic_cl)) > 1

    dynamic_cm = load_pitching_wing_exp_curve(:dynamic, :tripped, 0.25;
        data_dir=fixture_dir, quantity=:cm)
    @test length(dynamic_cm.aoa_deg) == length(dynamic_cm.values)
    @test pitching_wing_exp_has_hysteresis(dynamic_cm)
    @test any(!sweep.increasing for sweep in pitching_wing_exp_sweeps(dynamic_cm))

    dynamic_cl_09 = load_pitching_wing_exp_curve(:dynamic, :tripped, 0.9;
        data_dir=fixture_dir)
    dynamic_cm_09 = load_pitching_wing_exp_curve(:dynamic, :tripped, 0.9;
        data_dir=fixture_dir, quantity=:cm)
    @test !pitching_wing_exp_has_hysteresis(dynamic_cl_09)
    @test pitching_wing_exp_has_hysteresis(dynamic_cm_09)

    dynamic_cd = load_pitching_wing_exp_cd(data_dir=fixture_dir)
    @test !isempty(dynamic_cd.aoa_deg)
    @test length(dynamic_cd.aoa_deg) == length(dynamic_cd.values)
    @test !pitching_wing_exp_has_hysteresis(dynamic_cd)
    @test dynamic_cd.source == joinpath(fixture_dir, "dynamic_tripped_eta_all_cd.csv")
    dynamic_untripped_cd = load_pitching_wing_exp_cd(; data_dir=fixture_dir,
        regime=:dynamic, tripping=:untripped)
    @test dynamic_untripped_cd.source ==
        joinpath(fixture_dir, "dynamic_untripped_eta_all_cd.csv")
    @test pitching_wing_exp_has_hysteresis(dynamic_untripped_cd)

    exp = load_pitching_wing_exp(data_dir=fixture_dir)
    @test haskey(exp.static.untripped, 0.25)
    @test haskey(exp.static.tripped, 0.25)
    @test haskey(exp.dynamic.untripped, 0.25)
    @test haskey(exp.dynamic.tripped, 0.25)
    @test haskey(exp.dynamic.untripped, 0.9)
    @test haskey(exp.dynamic.tripped, 0.9)
    @test exp.static.untripped[0.25].cl.source == static.source
    @test exp.static.untripped.cd_all !== nothing
    @test exp.static.tripped.cd_all !== nothing
    @test exp.dynamic.untripped.cd_all !== nothing
    @test exp.dynamic.tripped.cd_all !== nothing
    @test length(exp.dynamic.tripped.cd_all.aoa_deg) == length(dynamic_cd.aoa_deg)
    @test length(exp.dynamic.untripped.cd_all.aoa_deg) == length(dynamic_untripped_cd.aoa_deg)
    @test exp.dynamic.untripped.cd_all.source == dynamic_untripped_cd.source
    @test length(exp.dynamic.untripped[0.25].cl.values) == length(dynamic_cl.values)
    @test length(exp.dynamic.tripped[0.25].cm.values) == length(dynamic_cm.values)
    @test exp.dynamic.untripped[0.25].cl.source == dynamic_cl.source
    @test length(exp.dynamic.tripped[0.9].cl.values) == length(dynamic_cl_09.values)
    @test length(exp.dynamic.tripped[0.9].cm.values) == length(dynamic_cm_09.values)
    @test pitching_wing_exp_has_hysteresis(exp.dynamic.tripped[0.9].cm)

    try
        import PythonPlot as plt
        plots = plot_pitching_wing_exp(exp; quantities=(:cl, :cm, :cd), show=false)
        @test haskey(pairs(plots), :cl)
        @test haskey(pairs(plots), :cm)
        @test haskey(pairs(plots), :cd)
        plt.close()
        plt.close()
        plt.close()
    catch err
        @test_skip "PythonPlot unavailable for pitching-wing plot smoke test: $(err)"
    end
end
