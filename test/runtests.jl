verbose = true

include("runtests_elements.jl")
include("runtests_semiinfiniteelements.jl")
include("runtests_grid.jl")
include("runtests_solvers.jl")
include("runtests_solvers2.jl")

if VERSION >= v"1.9"
    include("runtests_meshes.jl")
    include("runtests_meshes2.jl")
end
