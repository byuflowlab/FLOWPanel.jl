# Dependency-free syntax gate; execute only on an HPC compute node.
function check_syntax(ex, path)
    ex isa Expr || return
    ex.head in (:error, :incomplete) && error("Syntax error in $path: $ex")
    foreach(arg -> check_syntax(arg, path), ex.args)
end

files = ["common.jl", "phase1_case.jl", "fgs_cold_common.jl",
    "cold_precompile.jl", "cold_parse.jl", "rotor_hover_solver_cold.jl",
    "fgs_r4_diagnostics.jl",
    "rotor_hover_solver_cold_smoke.jl", "rotor_hover_solver_phase2_profile.jl",
    "../test/runtests_benchmark_cold.jl"]
for file in files
    path = normpath(joinpath(@__DIR__, file))
    check_syntax(Meta.parseall(read(path, String); filename=path), path)
    println("PARSE_OK ", file)
end
