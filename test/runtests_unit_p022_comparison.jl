using Test
using TOML

include(joinpath(@__DIR__, "..", "scripts", "p022_compare_vts_hybrid.jl"))

@testset "p022 independent VTS-vs-hybrid comparison artifact" begin
    mktempdir() do dir
        vdir = joinpath(dir, "vts")
        hdir = joinpath(dir, "hybrid")
        mkpath.(joinpath.((vdir, hdir), "monitors"))

        common = Dict{String,Any}(key => 1 for key in MATCH_KEYS)
        common["ground_enable"] = true
        common["nrotors"] = 1
        common["rotor_directions"] = [1]
        common["mesh_key"] = "test"
        common["mesh_file"] = "test.msh"
        common["particle_shedding"] = "overlap_pps"
        common["conversion"] = "legacy"
        common["ground_particle_policy"] = "none"
        common["vpm_arraytype"] = "array"
        common["gpu_influence_mode"] = "host"
        common["particle_body_overlap_action"] = "error"

        vts = copy(common)
        hybrid = copy(common)
        vts["formulation"] = "velocity"
        hybrid["formulation"] = "hybrid"
        vts["CT_window_mean_r1"] = 0.10
        hybrid["CT_window_mean_r1"] = 0.12
        vts["CQ_window_mean_r1"] = 0.01
        hybrid["CQ_window_mean_r1"] = 0.013

        vpath = joinpath(vdir, "case_metadata.toml")
        hpath = joinpath(hdir, "case_metadata.toml")
        open(io -> TOML.print(io, vts), vpath, "w")
        open(io -> TOML.print(io, hybrid), hpath, "w")
        for run_dir in (vdir, hdir)
            open(joinpath(run_dir, "monitors",
                    "run_monitor01_bound_circulation_system1.csv"), "w") do io
                println(io, "step,time,blade,section,r_over_R,circulation_te,circulation_slice")
                println(io, "0,0,1,1,0.5,2.0,2.0")
                println(io, "1,1,1,1,0.5,3.0,3.0")
            end
        end

        output = joinpath(dir, "comparison.toml")
        result = compare_vts_hybrid(vpath, hpath, output)
        @test result["matched_settings_verified"]
        @test result["delta_CT_window_mean_r1"] ≈ 0.02
        @test result["delta_CQ_window_mean_r1"] ≈ 0.003
        @test only(values(result["vts_circulation_final_mean"])) == 3.0
        @test TOML.parsefile(output)["artifact"] == "p022_vts_vs_hybrid"

        hybrid["NT"] = 2
        open(io -> TOML.print(io, hybrid), hpath, "w")
        @test_throws ErrorException compare_vts_hybrid(vpath, hpath, output)
    end
end
