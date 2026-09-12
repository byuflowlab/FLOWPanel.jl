# Harness contracts only: no rotor fixture, matrix assembly, or solver execution.
# Run with -t1 and -t4 to cover both legacy BLAS defaults.
using Test
include(joinpath(@__DIR__, "..", "benchmark", "common.jl"))
include(joinpath(@__DIR__, "..", "benchmark", "fgs_cold_common.jl"))

@testset "Benchmark BLAS control" begin
    original = LinearAlgebra.BLAS.get_num_threads()
    try
        for mode in ("single", "multi")
            withenv("THREADING_MODE" => mode,
                    "EXPECT_JULIA_THREADS" => string(Threads.nthreads())) do
                withenv("BENCH_BLAS_THREADS" => nothing) do
                    facts = assert_and_banner(IOBuffer())
                    @test facts.blas_threads == (mode == "single" ? 1 : Threads.nthreads())
                end
                for count in (1, 2)
                    withenv("BENCH_BLAS_THREADS" => string(count)) do
                        facts = assert_and_banner(IOBuffer())
                        @test facts.blas_threads == count
                        @test occursin("blas_threads   = $count", facts.text)
                    end
                end
                for invalid in ("0", "-1", "1.5", "abc", "")
                    withenv("BENCH_BLAS_THREADS" => invalid) do
                        @test_throws Exception assert_and_banner(IOBuffer())
                    end
                end
            end
        end
    finally
        LinearAlgebra.BLAS.set_num_threads(original)
    end
end

@testset "Cold benchmark configuration contracts" begin
    for rung in ("R1", "R2", "R3")
        baseline = cold_configs(rung, ["fgs", "krylov_ilu"], "baseline")
        @test length(baseline) == 2
        @test all(c -> c["rung"] == rung, baseline)
        screen = cold_configs(rung, ["fgs", "krylov_ilu"], "screen")
        @test length(screen) == 25
        @test length(unique(cold_id.(screen))) == length(screen)
        @test count(c -> get(c, "diagnostic", false), screen) == 1
        @test all(c -> c["kind"] == "fgs", cold_configs(rung, ["fgs"], "screen"))
        mktempdir() do dir
            path = joinpath(dir, "configs.toml")
            cold_write_toml(path, Dict("configs" => screen))
            loaded = cold_configs(rung, ["fgs", "krylov_ilu"], "verify"; file=path)
            @test loaded == screen
            @test cold_id.(loaded) == cold_id.(screen)
            @test_throws Exception cold_configs(rung == "R1" ? "R2" : "R1",
                ["fgs"], "verify"; file=path)
            cold_write_toml(path, baseline[1])
            @test cold_configs(rung, ["fgs"], "verify"; file=path) == baseline[1:1]
            @test isempty(cold_configs(rung, ["krylov_ilu"], "verify"; file=path))
        end
    end
    @test_throws Exception cold_configs("R1", ["fgs"], "verify")
    @test_throws Exception cold_seed("R1", "unknown")
    @test_throws Exception cold_seed("R4", "fgs")

    mktempdir() do dir
        out = joinpath(dir, "output")
        fixture = joinpath(dir, "fixture")
        withenv("STAGE" => "baseline", "RUNG" => "R1", "OUTDIR" => out,
                "BENCH_CASE_ROOT" => fixture, "CACHE_B" => "0", "SKIP_B" => "0",
                "FLOWPANEL_FILAMENT_REG" => "linegauss", "CONFIGS" => "fgs",
                "CONFIG_FILE" => "") do
            for (key, value) in ("STAGE" => "unknown", "RUNG" => "R4",
                    "OUTDIR" => "relative", "BENCH_CASE_ROOT" => "relative",
                    "CACHE_B" => "1", "SKIP_B" => "1",
                    "FLOWPANEL_FILAMENT_REG" => "unknown", "CONFIGS" => "unknown")
                withenv(key => value) do
                    @test_throws Exception cold_initialize!()
                    @test !ispath(out)
                    @test !ispath(fixture)
                end
            end
        end
    end
end
