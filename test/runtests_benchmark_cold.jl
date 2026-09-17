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
    for rung in ("R1", "R2", "R3", "R4")
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
                "diagnostic"=>1, "kind"=>"bad", "rung"=>"R5")
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

@testset "SCREEN_SET explicit roster" begin
    withenv("SCREEN_SET" => "inner:1,2,3,5") do
        configs = cold_configs("R2", ["fgs"], "screen")
        @test length(configs) == 5
        @test sort([c["inner"] for c in configs]) == [1, 2, 3, 5, 10]
        @test all(c -> c["inner"] isa Int && !(c["inner"] isa Bool), configs)
        @test all(c -> c["P"] == 8 && c["MAC"] == 0.4 && c["leaf"] == 100, configs)
        @test length(unique(cold_id.(configs))) == 5
        @test_throws Exception cold_configs("R2", ["fgs"], "baseline")
        @test_throws Exception cold_configs("R2", ["fgs"], "verify")
        # ILU seed has no "inner" field: unknown axis fails loudly.
        @test_throws Exception cold_configs("R2", ["krylov_ilu"], "screen")
        mktempdir() do dir
            path = joinpath(dir, "configs.toml")
            c = cold_seed("R2", "fgs"); c["tolerance"] = 1e-10
            cold_write_toml(path, Dict("configs" => [c]))
            @test_throws Exception cold_configs("R2", ["fgs"], "screen"; file=path)
        end
    end
    withenv("SCREEN_SET" => "MAC:0.3,0.5;leaf:25,50,200") do
        configs = cold_configs("R2", ["fgs"], "screen")
        @test length(configs) == 6   # seed + two MAC neighbors + three leaves
        @test all(c -> c["leaf"] isa Int && c["MAC"] isa Float64, configs)
    end
    for bad in ("inner", "inner:", "bogus:1,2", "inner:0", "inner:abc",
            "MAC:1.5", "inner:1;inner", "inner:1.5")
        withenv("SCREEN_SET" => bad) do
            @test_throws Exception cold_configs("R2", ["fgs"], "screen")
        end
    end
end

@testset "Screen around saved winners recalibrates every neighbor" begin
    mktempdir() do dir
        path = joinpath(dir, "bases.toml")
        bases = [merge(cold_seed("R2", "fgs"), Dict{String,Any}("inner"=>n, "tolerance"=>1e-9)) for n in (3,5)]
        cold_write_toml(path, Dict("configs"=>bases))
        original = read(path)
        withenv("SCREEN_BASE_FILE"=>path, "SCREEN_SET"=>"leaf:25,50,200") do
            configs = cold_configs("R2", ["fgs"], "screen")
            @test Set((c["inner"], c["leaf"]) for c in configs) == Set(Iterators.product((3,5), (25,50,100,200)))
            @test length(configs) == 8
            @test all(c -> c["tolerance"] == 0.0 && c["P"] == 8 && c["MAC"] == 0.4, configs)
            @test read(path) == original
            @test_throws Exception cold_configs("R3", ["fgs"], "screen")
            @test_throws Exception cold_configs("R2", ["fgs"], "baseline")
            @test_throws Exception cold_configs("R2", ["fgs"], "screen"; file=path)
            withenv("SCREEN_SET"=>"") do
                @test_throws Exception cold_configs("R2", ["fgs"], "screen")
            end
            withenv("SCREEN_BASE_FILE"=>joinpath(dir,"missing.toml")) do
                @test_throws Exception cold_configs("R2", ["fgs"], "screen")
            end
        end
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
            withenv("RUNG"=>"R4") do
                @test cold_preflight().configs[1]["rung"] == "R4"
            end
            for (key,value) in ("STAGE"=>"unknown", "RUNG"=>"R5",
                    "OUTDIR"=>"relative", "BENCH_CASE_ROOT"=>"relative",
                    "BENCH_CASE_ROOT"=>out, "BENCH_CASE_ROOT"=>joinpath(out,"nested"),
                    "CACHE_B"=>"1", "SKIP_B"=>"1", "CONFIGS"=>"unknown",
                    "MEMORY_GIB"=>"NaN", "MEMORY_GIB"=>"Inf", "MEMORY_GIB"=>"0",
                    "MEMORY_GIB"=>"501", "MEMORY_GIB"=>"abc",
                    "EXPECT_JULIA_THREADS"=>"0", "BENCH_BLAS_THREADS"=>"0",
                    "OMP_NUM_THREADS"=>"0", "THREADING_MODE"=>"bad",
                    "KNOBS_MODE"=>"../../outside", "PER_RUNG_DIR"=>"bad", "K_REPS"=>"bad",
                    "COLD_PREPARED_ONLY"=>"2", "COLD_PREPARED_ONLY"=>"bad",
                    "COLD_PROFILE_REPS"=>"0", "COLD_PROFILE_REPS"=>"bad",
                    "COLD_MIN_REPS"=>"0", "COLD_MIN_REPS"=>"bad",
                    "SCREEN_SET"=>"inner:1,2",
                    "SCREEN_BASE_FILE"=>joinpath(dir,"missing.toml"),
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
    @test cold_acceptance(1e-7,true,NaN,NaN,"R4").accepted
    @test !cold_acceptance(1e-7,false,NaN,NaN,"R4").accepted
    @test cold_acceptance(2e-6,false,1e-7,2e-6,"R4").authoritative_evaluator == "direct_fallback"
    @test cold_acceptance(2e-6,false,1e-7,2e-6,"R4").accepted
    @test !cold_acceptance(2e-6,false,2e-6,0.0,"R4").accepted
    @test !cold_acceptance(1e-7,true,1e-7,1.01e-7,"R4").accepted
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

module ColdScreenFailureHarness
import LinearAlgebra
include(joinpath(@__DIR__, "..", "benchmark", "fgs_cold_common.jl"))
rung = "R1"
cold_selected_file = ""
# One roster point fails; the seed (inner=10) succeeds with a stub summary.
cold_benchmark(c, dir; selected=false) = c["inner"] != 10 ?
    error("injected screen candidate failure") :
    [(; mode="prepared", minimum_seconds=1.0, median_seconds=1.0,
        maximum_seconds=1.0, spread_seconds=0.0, repetitions=5, eligible=true)]
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

@testset "Screen records failed candidates and continues" begin
    seed = cold_seed("R1", "fgs")
    bad = copy(seed); bad["inner"] = 3
    # Partial failure: the failed roster point is recorded, the survivor wins.
    mktempdir() do dir
        @test ColdScreenFailureHarness.cold_run([bad, seed], dir, "screen") === nothing
        statuses = Dict(joinpath(root,f) => TOML.parsefile(joinpath(root,f))["status"]
            for (root,_,files) in walkdir(dir) for f in files if f == "status.toml")
        @test sort(collect(values(statuses))) == ["completed", "failed"]
        selected = TOML.parsefile(joinpath(dir,"selected.toml"))["configs"]
        @test length(selected) == 1 && selected[1]["inner"] == 10
    end
    # Total failure: every candidate recorded, no selection is reported.
    mktempdir() do dir
        bad2 = copy(seed); bad2["inner"] = 5
        @test_throws ErrorException ColdScreenFailureHarness.cold_run([bad, bad2], dir, "screen")
        statuses = [joinpath(root,f) for (root,_,files) in walkdir(dir) for f in files if f == "status.toml"]
        @test length(statuses) == 2
        @test all(TOML.parsefile(s)["status"] == "failed" for s in statuses)
        @test !isfile(joinpath(dir,"selected.toml"))
    end
end
