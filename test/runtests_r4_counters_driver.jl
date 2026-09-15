#!/usr/bin/env julia
using Test
driver = joinpath(@__DIR__, "..", "benchmark", "fgs_r4_counters.jl")
nodes(x) = x isa Expr ? vcat([x], mapreduce(nodes, vcat, x.args; init=Any[])) : Any[x]
findfun(ast, name) = only(filter(x -> x isa Expr && x.head === :function &&
    x.args[1] isa Expr && x.args[1].head === :call && x.args[1].args[1] === name, nodes(ast)))
ast = Meta.parseall(read(driver, String))
@test isempty(filter(x -> x isa Expr && x.head in (:error, :incomplete), nodes(ast)))

module EntryProbe
    run_counters(configs, out, stage) = (dynamic_reset!(); :reached)
    function cold_initialize!()
        @eval dynamic_reset!() = nothing
        (nothing, nothing, nothing)
    end
end
Core.eval(EntryProbe, findfun(ast, :counters_main))
@test EntryProbe.counters_main() === :reached

module ObserverProbe
    ticks = Ref(0)
    activity_snapshot() = Dict(1 => (ticks[], 2), 2 => (0, 3))
end
Core.eval(ObserverProbe, findfun(ast, :activity_observer))
rows = NamedTuple[]
observe = ObserverProbe.activity_observer(rows)
observe(:fmm, :start)
ObserverProbe.ticks[] = 7
observe(:fmm, :stop)
@test [r.cpu_ticks for r in rows] == [7, 0]
@test all(r -> r.stage == "fmm" && r.complete_endpoints && r.span_seconds >= 0, rows)
@test_throws ErrorException observe(:fmm, :stop)
observe(:residual, :start)
@test_throws ErrorException observe(:fmm, :start)
observe(:residual, :stop)

module ProtocolProbe end
Core.eval(ProtocolProbe, findfun(ast, :counter_command))
control = IOBuffer()
ProtocolProbe.counter_command(control, IOBuffer("ack\n"), "enable")
@test String(take!(control)) == "enable\n"
@test_throws ErrorException ProtocolProbe.counter_command(control, IOBuffer("bad\n"), "disable")

if Sys.islinux()
    Core.eval(Main, findfun(ast, :activity_snapshot))
    snapshot = activity_snapshot()
    @test !isempty(snapshot)
    @test all(v -> v[1] >= 0 && v[2] >= 0, values(snapshot))
end
run(`bash -n $(joinpath(dirname(driver), "run_r4_counters.slurm.sh"))`)
println("PASS counter driver parse, dynamic initialization, stage pairing, protocol, and launcher")
