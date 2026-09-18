#!/usr/bin/env julia
using Test
driver = joinpath(@__DIR__, "..", "benchmark", "fgs_r4_chunked_ab.jl")
nodes(x) = x isa Expr ? vcat([x], mapreduce(nodes, vcat, x.args; init=Any[])) : Any[x]
findfun(ast, name) = only(filter(x -> x isa Expr && x.head === :function &&
    x.args[1] isa Expr && x.args[1].head === :call && x.args[1].args[1] === name, nodes(ast)))
ast = Meta.parseall(read(driver, String))
@test isempty(filter(x -> x isa Expr && x.head in (:error, :incomplete), nodes(ast)))

module EntryProbe
    run_chunked_ab(configs, out, stage) = (dynamic_reset!(); :reached)
    function cold_initialize!()
        @eval dynamic_reset!() = nothing
        (nothing, nothing, nothing)
    end
end
Core.eval(EntryProbe, findfun(ast, :ab_main))
@test EntryProbe.ab_main() === :reached

module PairProbe end
Core.eval(PairProbe, findfun(ast, :ab_check_pair))
lex = Dict{String,Any}("kind" => "fgs", "sweep_order" => "lexicographic",
    "tolerance" => 3.5e-7, "leaf" => 100, "P" => 8)
chunked = copy(lex); chunked["sweep_order"] = "chunked"
chunked["tolerance"] = 1.1e-7; chunked["chunks"] = 64
@test PairProbe.ab_check_pair(lex, chunked) === nothing
bad = copy(chunked); bad["leaf"] = 50
@test_throws ErrorException PairProbe.ab_check_pair(lex, bad)
bad = copy(chunked); bad["tolerance"] = 0.0
@test_throws ErrorException PairProbe.ab_check_pair(lex, bad)
bad = copy(chunked); delete!(bad, "P")
@test_throws ErrorException PairProbe.ab_check_pair(lex, bad)
# chunks mismatch cases: missing, non-positive, non-integer, or on the lex side
bad = copy(chunked); delete!(bad, "chunks")
@test_throws ErrorException PairProbe.ab_check_pair(lex, bad)
bad = copy(chunked); bad["chunks"] = 0
@test_throws ErrorException PairProbe.ab_check_pair(lex, bad)
bad = copy(chunked); bad["chunks"] = 64.0
@test_throws ErrorException PairProbe.ab_check_pair(lex, bad)
badlex = copy(lex); badlex["chunks"] = 64
@test_throws ErrorException PairProbe.ab_check_pair(badlex, chunked)
@test_throws ErrorException PairProbe.ab_check_pair(chunked, chunked)
@test_throws ErrorException PairProbe.ab_check_pair(lex, lex)

# The alternating-batch trial loop must interleave the two orders evenly.
trials = findfun(ast, :ab_trials)
@test occursin("isodd(batch)", string(trials))
@test occursin("2batches", string(trials))
# Performance trials stay uninstrumented; perf never appears in this campaign.
@test !occursin("COLD_PERF", read(driver, String))
launcher = joinpath(dirname(driver), "run_r4_chunked_ab.slurm.sh")
@test !occursin(r"\bperf\b", read(launcher, String))
run(`bash -n $launcher`)
println("PASS chunked A/B driver parse, dynamic initialization, pair check, alternation, and launcher")
