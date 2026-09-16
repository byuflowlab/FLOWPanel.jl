#!/usr/bin/env julia

driver = length(ARGS) >= 1 ? abspath(ARGS[1]) : joinpath(@__DIR__, "..", "benchmark", "fgs_r4_diagnostics.jl")
launcher = joinpath(dirname(driver), "run_r4_diagnostics.slurm.sh")
nodes(x) = x isa Expr ? vcat([x], mapreduce(nodes, vcat, x.args; init=Any[])) : Any[x]
findfun(ast, name) = only(filter(x -> x isa Expr && x.head === :function &&
    ((x.args[1] === name) || (x.args[1] isa Expr && x.args[1].head === :call && x.args[1].args[1] === name)), nodes(ast)))
ast = Meta.parseall(read(driver, String))
bad = filter(x -> x isa Expr && x.head in (:error, :incomplete), nodes(ast))
isempty(bad) || error("driver AST contains parse error/incomplete nodes")
main_ast, dep_ast = findfun(ast, :main), findfun(ast, :dependency_census)
println("PASS parsed_ast nodes=$(length(nodes(ast)))")

module NewActual
    run_diagnostics(configs, out, stage) = (reset_cold!(); :reached)
    function cold_initialize!()
        @eval begin
            reset_cold!() = :reset
        end
        return ([Dict{String,Any}()], :out, :stage)
    end
end
Core.eval(NewActual, main_ast)
NewActual.main() === :reached || error("actual main did not reach dynamic helper")

bad_body = deepcopy(main_ast.args[2])
for i in eachindex(bad_body.args)
    s = bad_body.args[i]
    if s isa Expr && s.head === :call && occursin("invokelatest", string(s.args[1]))
        bad_body.args[i] = Expr(:call, s.args[2], s.args[3:end]...)
    end
end
bad_main = Expr(:function, main_ast.args[1], bad_body)
module OldActual
    run_diagnostics(configs, out, stage) = (reset_cold!(); :reached)
    function cold_initialize!()
        @eval begin
            reset_cold!() = :reset
        end
        return ([Dict{String,Any}()], :out, :stage)
    end
end
Core.eval(OldActual, bad_main)
old_failed = false
try
    OldActual.main()
catch err
    global old_failed = err isa MethodError && occursin("too new", sprint(showerror, err))
end
old_failed || error("removing invokelatest did not reproduce world-age failure")
println("PASS actual_main new_invokelatest=true removed_invokelatest_fails=true")

module CensusProbe end
Core.eval(CensusProbe, dep_ast)
for entries in ([(2, 1)], [(4, 1), (4, 1)], [(4, 1), (1, 3), (5, 2)])
    # Branch 4 spans leaves 2 and 3; branch 5 is empty. Include repeated
    # entries, asymmetric edges, and reverse edges in the independent oracle.
    fgs = (targets_by_branch=[1:2, 3:5, 6:7, 3:7, 1:0],
           source_tree=(leaf_index=1:3,),
           index_map=[findall(e -> e[2] == i, entries) for i in 1:3],
           direct_list=entries)
    writes, adjacency, colors = CensusProbe.dependency_census(fgs)
    expected = [Set{Int}() for _ in 1:3]
    for (target, source) in entries, dependent in 1:3
        source == dependent && continue
        if !isempty(intersect(fgs.targets_by_branch[target], fgs.targets_by_branch[dependent]))
            push!(expected[source], dependent)
        end
    end
    expected_adjacency = [union(expected[i], Set(j for j in 1:3 if i in expected[j])) for i in 1:3]
    writes == expected || error("directed write/read overlap mismatch")
    adjacency == expected_adjacency || error("symmetrized conflict mismatch")
    all(colors[i] != colors[j] for i in eachindex(adjacency) for j in adjacency[i]) || error("improper coloring")
end
println("PASS dependency_census asymmetric/nonleaf/duplicate/empty cases match brute-force overlap")

isfile(launcher) && run(`bash -n $launcher`)
println("PASS bash_n launcher=$(isfile(launcher))")
