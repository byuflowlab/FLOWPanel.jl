using Test
import FLOWPanel as pnl
import GeometricTools as gt
using WriteVTK

include("test_helpers.jl")

include("runtests_unit_fmm.jl")
include("runtests_unit_body.jl")
include("runtests_unit_solver.jl")
include("runtests_unit_liftingbody.jl")
include("runtests_unit_wake.jl")
include("runtests_unit_simulate.jl")
include("runtests_unit_warmstart.jl")
include("runtests_analytical.jl")
