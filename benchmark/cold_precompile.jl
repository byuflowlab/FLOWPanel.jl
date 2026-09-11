# Separate compute-node preparation; never part of timing or test processes.
using Pkg
Pkg.instantiate(; allow_autoprecomp=false)
Pkg.precompile()
include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "fgs_cold_common.jl"))
assert_and_banner()
cold_packages()
