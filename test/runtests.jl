using Test
import FLOWPanel as pnl
using WriteVTK

include("test_helpers.jl")

include("runtests_unit_fmm.jl")
include("runtests_unit_kernel_gradient.jl")
include("runtests_unit_regularization.jl")
include("runtests_unit_body.jl")
include("runtests_unit_solver.jl")
include("runtests_unit_instrumentation.jl")
include("runtests_unit_solver_history.jl")
include("runtests_unit_liftingbody.jl")
include("runtests_unit_wake.jl")
include("runtests_unit_postprocess.jl")
include("runtests_unit_simulate.jl")
include("runtests_unit_kutta.jl")
include("runtests_unit_kutta_routeb.jl")
include("runtests_unit_added_mass.jl")
# include("dirichlet_potential_test.jl")
include("runtests_unit_warmstart.jl")
include("runtests_unit_replay.jl")
include("runtests_example_pitching_wing.jl")
include("runtests_example_pitching_wing_convergence.jl")
include("runtests_example_pitching_wing_pressure_comparison.jl")
include("runtests_example_suddenly_started_wing.jl")
include("runtests_example_dji9443_trailing_edge.jl")
include("runtests_analytical.jl")
