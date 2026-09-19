#!/usr/bin/env julia
using Test
driver = joinpath(@__DIR__, "..", "benchmark", "fgs_r4_dagteam_ab.jl")
nodes(x) = x isa Expr ? vcat([x], mapreduce(nodes, vcat, x.args; init=Any[])) : Any[x]
findfun(ast, name) = only(filter(x -> x isa Expr && x.head === :function &&
    x.args[1] isa Expr && x.args[1].head === :call && x.args[1].args[1] === name, nodes(ast)))
ast = Meta.parseall(read(driver, String))
@test isempty(filter(x -> x isa Expr && x.head in (:error, :incomplete), nodes(ast)))

module EntryProbe
    run_dagteam_ab(configs, out, stage) = (dynamic_reset!(); :reached)
    function cold_initialize!()
        @eval dynamic_reset!() = nothing
        (nothing, nothing, nothing)
    end
end
Core.eval(EntryProbe, findfun(ast, :ab_main))
@test EntryProbe.ab_main() === :reached

module PairProbe end
Core.eval(PairProbe, findfun(ast, :ab_check_pair))
colored = Dict{String,Any}("kind" => "fgs", "sweep_order" => "colored",
    "tolerance" => 3.5e-7, "leaf" => 100, "P" => 8)
dagteam = copy(colored); dagteam["sweep_order"] = "dagteam"
dagteam["tolerance"] = 1.1e-7; dagteam["dagteam_precision"] = "f32full"
@test PairProbe.ab_check_pair(colored, dagteam) === nothing
for precision in ("f64", "f32conv")
    ok = copy(dagteam); ok["dagteam_precision"] = precision
    @test PairProbe.ab_check_pair(colored, ok) === nothing
end
bad = copy(dagteam); bad["leaf"] = 50
@test_throws ErrorException PairProbe.ab_check_pair(colored, bad)
bad = copy(dagteam); bad["tolerance"] = 0.0
@test_throws ErrorException PairProbe.ab_check_pair(colored, bad)
badc = copy(colored); badc["tolerance"] = 0.0
@test_throws ErrorException PairProbe.ab_check_pair(badc, dagteam)
bad = copy(dagteam); delete!(bad, "P")
@test_throws ErrorException PairProbe.ab_check_pair(colored, bad)
# dagteam_precision cases: missing, invalid, or on the colored side
bad = copy(dagteam); delete!(bad, "dagteam_precision")
@test_throws ErrorException PairProbe.ab_check_pair(colored, bad)
bad = copy(dagteam); bad["dagteam_precision"] = "f16"
@test_throws ErrorException PairProbe.ab_check_pair(colored, bad)
badc = copy(colored); badc["dagteam_precision"] = "f32full"
@test_throws ErrorException PairProbe.ab_check_pair(badc, dagteam)
# chunks must not appear on either side
bad = copy(dagteam); bad["chunks"] = 64
@test_throws ErrorException PairProbe.ab_check_pair(colored, bad)
badc = copy(colored); badc["chunks"] = 64
@test_throws ErrorException PairProbe.ab_check_pair(badc, dagteam)
@test_throws ErrorException PairProbe.ab_check_pair(dagteam, dagteam)
@test_throws ErrorException PairProbe.ab_check_pair(colored, colored)

# The alternating-batch trial loop must interleave the two orders evenly.
trials = findfun(ast, :ab_trials)
@test occursin("isodd(batch)", string(trials))
@test occursin("2batches", string(trials))
# Performance trials stay uninstrumented; perf never appears in this campaign.
@test !occursin("COLD_PERF", read(driver, String))
for launcher in ("run_r4_dagteam_ab.slurm.sh", "run_r4_dagteam_numgate.slurm.sh")
    path = joinpath(dirname(driver), launcher)
    @test !occursin(r"\bperf\b", read(path, String))
    run(`bash -n $path`)
end
println("PASS dagteam A/B driver parse, dynamic initialization, pair check, alternation, and launchers")
