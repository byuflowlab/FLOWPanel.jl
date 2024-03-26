using Pkg
this_dir = @__DIR__
Pkg.activate(normpath(joinpath(this_dir,"..")))

using FLOWPanel
using FLOWPanel.StaticArrays
using FLOWPanel.WriteVTK
using FLOWPanel.FastMultipole
using Random, Test, ForwardDiff
using LinearAlgebra: norm, cross, dot
Random.seed!(123)

include("auxiliary_functions.jl")
# include("kernel_test.jl")
# include("grid_test.jl")
include("fmm_test.jl")
# include("solver_test.jl")