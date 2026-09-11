# Run sequentially on HPC at -t1 and -t4, BLAS=1 before process startup.
using Test
include(joinpath(@__DIR__, "..", "benchmark", "common.jl"))
include(joinpath(@__DIR__, "..", "benchmark", "fgs_cold_common.jl"))

@testset "BLAS control persists after actual work" begin
    facts = assert_and_banner(IOBuffer())
    expected = parse(Int, ENV["BENCH_BLAS_THREADS"])
    @test facts.blas_threads == expected
    for _ in 1:3
        a = ones(128,128)
        @test sum(a*a) == 128^3
        Threads.@threads for i in 1:Threads.nthreads()
            b = ones(64,64)
            @assert sum(b*b) == 64^3
        end
        @test LinearAlgebra.BLAS.get_num_threads() == expected
    end
    for bad in ("0", "-1", "1.5", "abc", "")
        withenv("BENCH_BLAS_THREADS" => bad) do
            @test_throws Exception assert_and_banner(IOBuffer())
        end
    end
end

@testset "Configuration contracts and selected immutability" begin
    for rung in ("R1", "R2", "R3")
        configs = cold_configs(rung, ["fgs", "krylov_ilu"], "screen")
        @test length(configs) == 25
        @test length(unique(cold_id.(configs))) == 25
        @test count(c -> get(c, "diagnostic", false), configs) == 1
        for c in configs
            @test cold_check_config(c) === c
            c["kind"] == "fgs" && (c["tolerance"] = 1e-10)
        end
        mktempdir() do dir
            path = joinpath(dir, "configs.toml")
            cold_write_toml(path, Dict("configs" => configs))
            original = read(path)
            for stage in ("baseline", "screen", "verify")
                loaded = cold_configs(rung, ["fgs", "krylov_ilu"], stage; file=path)
                @test loaded == configs
                for c in loaded
                    before = deepcopy(c)
                    # No fixture globals exist: accidentally calling calibration fails.
                    cold_calibrate_selected!(c, dir; selected=true)
                    @test c == before
                end
                @test read(path) == original
                @test readdir(dir) == ["configs.toml"]
            end
        end
    end
    for kind in ("fgs", "krylov_ilu")
        c = cold_seed("R1", kind)
        for key in keys(c)
            bad = copy(c); delete!(bad, key)
            @test_throws Exception cold_check_config(bad)
            bad = copy(c); bad[key] = nothing
            @test_throws Exception cold_check_config(bad)
        end
        for (key,value) in ("unknown"=>1, "P"=>true, "P"=>1.5, "P"=>0,
                "leaf"=>-1, "MAC"=>NaN, "MAC"=>Inf, "MAC"=>1.1,
                "diagnostic"=>1, "kind"=>"bad", "rung"=>"R4")
            bad = copy(c); bad[key] = value
            @test_throws Exception cold_check_config(bad)
        end
    end
    @test_throws Exception cold_check_config(cold_seed("R1","fgs"); selected=true)
    @test_throws Exception cold_configs("R1", ["fgs"], "verify")
    for (kind,key,value) in (("fgs","sweep_order","bad"), ("fgs","rlx",2.0),
            ("fgs","tolerance",-1.0), ("krylov_ilu","rtol",0.0),
            ("krylov_ilu","atol",-1.0), ("krylov_ilu","ilu_MAC",0.0),
            ("krylov_ilu","memory",0), ("krylov_ilu","diagonal_shift",-1.0))
        c = cold_seed("R1",kind); c[key] = value
        @test_throws Exception cold_check_config(c)
    end
end

@testset "Invalid input has no filesystem effects" begin
    mktempdir() do dir
        out, fixture = joinpath(dir,"output"), joinpath(dir,"fixture")
        withenv("STAGE"=>"baseline", "RUNG"=>"R1", "OUTDIR"=>out,
                "BENCH_CASE_ROOT"=>fixture, "CACHE_B"=>"0", "SKIP_B"=>"0",
                "FLOWPANEL_FILAMENT_REG"=>"linegauss", "CONFIGS"=>"fgs",
                "CONFIG_FILE"=>"", "MEMORY_GIB"=>"500") do
            @test cold_preflight().memory == 500*2.0^30
            for (key,value) in ("STAGE"=>"unknown", "RUNG"=>"R4",
                    "OUTDIR"=>"relative", "BENCH_CASE_ROOT"=>"relative",
                    "BENCH_CASE_ROOT"=>out, "BENCH_CASE_ROOT"=>joinpath(out,"nested"),
                    "CACHE_B"=>"1", "SKIP_B"=>"1", "CONFIGS"=>"unknown",
                    "MEMORY_GIB"=>"NaN", "MEMORY_GIB"=>"Inf", "MEMORY_GIB"=>"0",
                    "MEMORY_GIB"=>"501", "MEMORY_GIB"=>"abc",
                    "EXPECT_JULIA_THREADS"=>"0", "BENCH_BLAS_THREADS"=>"0",
                    "OMP_NUM_THREADS"=>"0", "THREADING_MODE"=>"bad",
                    "KNOBS_MODE"=>"../../outside", "PER_RUNG_DIR"=>"bad", "K_REPS"=>"bad",
                    "CONFIG_FILE"=>joinpath(dir,"missing.toml"))
                withenv(key=>value) do
                    @test_throws Exception cold_initialize!()
                    @test isempty(readdir(dir))
                end
            end
            @test_throws Exception cold_initialize!(; profile=true)
            @test isempty(readdir(dir))
            path = joinpath(dir,"config.toml")
            for doc in (Dict("configs"=>[]), Dict("configs"=>1), Dict("unknown"=>1),
                    Dict("configs"=>[cold_seed("R1","fgs")]),
                    Dict("configs"=>[cold_seed("R1","krylov_ilu")], "unknown"=>1))
                cold_write_toml(path,doc)
                withenv("CONFIG_FILE"=>path) do
                    @test_throws Exception cold_initialize!()
                    @test readdir(dir) == ["config.toml"]
                end
            end
            mkdir(fixture); write(joinpath(fixture,"keep"),"preserve")
            @test_throws Exception cold_initialize!()
            @test !ispath(out)
            @test read(joinpath(fixture,"keep"),String) == "preserve"
        end
    end
end

@testset "Evaluator acceptance and failure propagation" begin
    @test cold_acceptance(1e-6,true,NaN,NaN,"R1").accepted
    @test !cold_acceptance(1.01e-6,true,NaN,NaN,"R1").accepted
    @test !cold_acceptance(1e-8,false,NaN,NaN,"R1").accepted
    @test !cold_acceptance(NaN,true,NaN,NaN,"R1").accepted
    @test cold_acceptance(1e-7,true,1e-7,1e-7,"R1").accepted
    @test !cold_acceptance(1e-7,true,1e-7,1.01e-7,"R1").accepted
    @test !cold_acceptance(2e-6,true,1e-8,1e-8,"R1").accepted
    for rung in ("R1","R2")
        e = cold_acceptance(2e-6,false,1e-7,2e-6,rung)
        @test e.accepted && e.authoritative_evaluator == "direct_fallback"
        @test e.authoritative_rel_l2 == 1e-7
    end
    @test !cold_acceptance(1e-7,false,1e-7,1e-8,"R3").accepted
    @test_throws ErrorException cold_require(false,"profile verification failed")
    @test cold_require(true,"valid") === nothing
end

module ColdFailureHarness
import LinearAlgebra
include(joinpath(@__DIR__, "..", "benchmark", "fgs_cold_common.jl"))
rung = "R1"
cold_selected_file = "/selected.toml"
cold_profile(c, dir; selected=true) = error("injected profile verification failure")
cold_benchmark(c, dir; selected=false) = error("injected selected solve failure")
end

@testset "Failed requested execution cannot report completion" begin
    for profile in (false,true), stage in ("baseline","verify")
        mktempdir() do dir
            c = cold_seed("R1","krylov_ilu"); before = deepcopy(c)
            @test_throws ErrorException ColdFailureHarness.cold_run([c],dir,stage; profile)
            @test c == before
            statuses = [joinpath(root,f) for (root,_,files) in walkdir(dir) for f in files if f == "status.toml"]
            @test length(statuses) == 1
            @test TOML.parsefile(only(statuses))["status"] == "failed"
            @test !isfile(joinpath(dir,"selected.toml"))
        end
    end
end
