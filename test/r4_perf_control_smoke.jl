#!/usr/bin/env julia
"""Exercise perf's installed FIFO control protocol from Julia before fixture work."""
using Test

length(ARGS) == 2 || error("usage: r4_perf_control_smoke.jl CONTROL_FIFO ACK_FIFO")
driver = joinpath(@__DIR__, "..", "benchmark", "fgs_r4_counters.jl")
nodes(x) = x isa Expr ? vcat([x], mapreduce(nodes, vcat, x.args; init=Any[])) : Any[x]
findfun(ast, name) = only(filter(x -> x isa Expr && x.head === :function &&
    x.args[1] isa Expr && x.args[1].head === :call && x.args[1].args[1] === name, nodes(ast)))
findstruct(ast, name) = only(filter(x -> x isa Expr && x.head === :struct && x.args[2] === name, nodes(ast)))

module PerfSmokeProbe end
ast = Meta.parseall(read(driver, String))
Core.eval(PerfSmokeProbe, findstruct(ast, :CounterPollFD))
Core.eval(PerfSmokeProbe, findfun(ast, :counter_command))

Base.@noinline function smoke_work(n)
    total = 0
    for i in 1:n
        total += i * i
    end
    total
end

open(ARGS[1], "r+") do control
    open(ARGS[2], "r+") do acknowledgement
        PerfSmokeProbe.counter_command(control, acknowledgement, "disable")
        PerfSmokeProbe.counter_command(control, acknowledgement, "enable")
        result = smoke_work(parse(Int, get(ENV, "PERF_SMOKE_N", "1000000")))
        PerfSmokeProbe.counter_command(control, acknowledgement, "disable")
        @test result > 0
    end
end
println("PASS Julia perf enable/disable acknowledgements")
