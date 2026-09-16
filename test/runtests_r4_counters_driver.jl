#!/usr/bin/env julia
using Test
driver = joinpath(@__DIR__, "..", "benchmark", "fgs_r4_counters.jl")
nodes(x) = x isa Expr ? vcat([x], mapreduce(nodes, vcat, x.args; init=Any[])) : Any[x]
findfun(ast, name) = only(filter(x -> x isa Expr && x.head === :function &&
    x.args[1] isa Expr && x.args[1].head === :call && x.args[1].args[1] === name, nodes(ast)))
findstruct(ast, name) = only(filter(x -> x isa Expr && x.head === :struct && x.args[2] === name, nodes(ast)))
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
Core.eval(ProtocolProbe, findstruct(ast, :CounterPollFD))
Core.eval(ProtocolProbe, findfun(ast, :counter_command))
control = IOBuffer()
ProtocolProbe.counter_command(control, IOBuffer("ack\n"), "enable")
@test String(take!(control)) == "enable\n"
@test_throws ErrorException ProtocolProbe.counter_command(control, IOBuffer("bad\n"), "disable")

if Sys.isunix()
    fifo_dir = mktempdir()
    control_fifo = joinpath(fifo_dir, "control.fifo")
    acknowledgement_fifo = joinpath(fifo_dir, "acknowledgement.fifo")
    run(`mkfifo $control_fifo $acknowledgement_fifo`)
    function fifo_protocol(command, response; timeout_seconds=0.2)
        open(control_fifo, "r+") do fifo_control
            open(acknowledgement_fifo, "r+") do fifo_acknowledgement
                response === nothing || begin
                    open(acknowledgement_fifo, "r+") do peer
                        write(peer, response)
                        flush(peer)
                    end
                end
                result = try
                    ProtocolProbe.counter_command(fifo_control, fifo_acknowledgement, command;
                        timeout_seconds)
                catch err
                    err
                end
                result
            end
        end
    end
    @test fifo_protocol("enable", "ack\n") === nothing
    missing_elapsed = @elapsed missing_result = fifo_protocol("disable", nothing; timeout_seconds=0.1)
    @test missing_result isa ErrorException
    @test missing_elapsed < 1.0
    @test fifo_protocol("disable", "nack\n") isa ErrorException
    partial_elapsed = @elapsed partial_result = fifo_protocol("disable", "a"; timeout_seconds=0.1)
    @test partial_result isa ErrorException
    @test partial_elapsed < 1.0
end

if Sys.islinux()
    Core.eval(Main, findfun(ast, :activity_snapshot))
    snapshot = activity_snapshot()
    @test !isempty(snapshot)
    @test all(v -> v[1] >= 0 && v[2] >= 0, values(snapshot))
end
run(`bash -n $(joinpath(dirname(driver), "run_r4_counters.slurm.sh"))`)
println("PASS counter driver parse, dynamic initialization, stage pairing, protocol, and launcher")
