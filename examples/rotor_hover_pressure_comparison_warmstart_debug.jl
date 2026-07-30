## Warmstart + SFS-blow-up tracer for rotor_hover_pressure_comparison.
##
## Re-runs the final few timesteps of the canonical run forward in time
## (NOT a monitor-only replay) so that print-statement debugging can be
## inserted into FLOWPanel/FLOWVPM and the SFS / J / Gamma growth that
## causes the simulation to diverge can be localized.
##
## Reuses the geometry, frames, body, wake, and monitors from
## `rotor_hover_pressure_comparison.jl` by including that file with
## RHPC_SETUP_ONLY=1 to skip its simulate! call.
##
## Env knobs (read here, in addition to the canonical script's own):
##   RESTART_STEP      Step to warmstart from. Default = total_saved - 10.
##   DEBUG_NAME        VTK name suffix to write under (default
##                     `<run_name>_warmstart_debug`).
##   DEBUG_PATH        Output dir override.
##   SFS_DEBUG=1       Turns on SFSBlowupTracer prints.
##   SFS_DEBUG_VERBOSE / SFS_DEBUG_FOCUS_INDEX / SFS_DEBUG_*  (see tracer)
##
## Example invocation:
##   SFS_DEBUG=1 SFS_DEBUG_VERBOSE=1 RESTART_STEP=40 \
##     julia --project examples/rotor_hover_pressure_comparison_warmstart_debug.jl

ENV["RHPC_SETUP_ONLY"] = "1"
import FLOWPanel as pnl
include(joinpath(@__DIR__, "rotor_hover_pressure_comparison.jl"))
include(joinpath(@__DIR__, "..", "debug", "sfs_blowup_tracer.jl"))

# After the include, the following names from the canonical script are in
# scope: run_name, save_path, systems, wakes, frames, maneuver!, Uinf,
# t_range, body_solvers, backend, backend_wake, monitors,
# set_Das_min_kinematic_displacement.
#
# Note: `monitors` is a Tuple and `audit_monitors` requires the tracer to
# come AFTER any monitor that provides :P / :F it might depend on. The
# tracer provides/requires nothing, so it is safe to append.

restart_path = save_path
restart_name = run_name
debug_name   = get(ENV, "DEBUG_NAME", run_name * "_warmstart_debug")
debug_path   = get(ENV, "DEBUG_PATH", joinpath("data", debug_name))

# Default restart_step: ten steps before the final saved step.
saved_pvd = joinpath(restart_path, restart_name * "_body1.pvd")
_, saved_idxs = pnl._read_pvd_steps(saved_pvd)
default_restart = saved_idxs[max(1, length(saved_idxs) - 10)]
restart_step = parse(Int, get(ENV, "RESTART_STEP", string(default_restart)))

println("\n=== Warmstart debug run ===")
println("  restart_path : $(restart_path)")
println("  restart_name : $(restart_name)")
println("  restart_step : $(restart_step) (saved steps: $(first(saved_idxs))..$(last(saved_idxs)))")
println("  debug_path   : $(debug_path)")
println("  debug_name   : $(debug_name)")
println("  t_range len  : $(length(t_range))")
println("  SFS_DEBUG    : $(get(ENV, "SFS_DEBUG", "0"))")

# --- Diagnostic: force wake-as-source FMM to behave like direct summation ---
# Gated by WAKE_DIRECT=1. Overrides both the panel-side backend_wake (used in
# _steady_aerodynamics! for wake_sources -> targets) and the VPM-internal
# pfield.fmm (used during RK3 integration for pfield-on-pfield UJ + SFS
# pre-calc). The panel-on-panel `backend` is intentionally NOT touched.
wake_direct      = parse(Bool, get(ENV, "WAKE_DIRECT", "false"))
wake_direct_leaf = parse(Int,  get(ENV, "WAKE_DIRECT_LEAF", "100000"))
if wake_direct
    println("\n*** WAKE_DIRECT=1: forcing wake-as-source FMM to leaf_size=$(wake_direct_leaf) ***")

    backend_wake = pnl.FastMultipoleBackend(;
        expansion_order      = backend_wake.expansion_order,
        multipole_acceptance = backend_wake.multipole_acceptance,
        leaf_size            = wake_direct_leaf,
    )

    pfield = wakes[1].pfield
    old_fmm = pfield.fmm
    pfield.fmm = pnl.FLOWVPM.FMM(;
        p                      = old_fmm.p,
        ncrit                  = wake_direct_leaf,
        min_ncrit              = wake_direct_leaf,
        theta                  = old_fmm.theta,
        shrink_recenter        = old_fmm.shrink_recenter,
        relative_tolerance     = old_fmm.relative_tolerance,
        absolute_tolerance     = old_fmm.absolute_tolerance,
        autotune_p             = false,
        autotune_ncrit         = false,
        autotune_reg_error     = old_fmm.autotune_reg_error,
        default_rho_over_sigma = old_fmm.default_rho_over_sigma,
    )

    println("  backend_wake.leaf_size       = $(backend_wake.leaf_size)")
    println("  pfield.fmm.ncrit / min_ncrit = $(pfield.fmm.ncrit) / $(pfield.fmm.min_ncrit)")
    println("  pfield.fmm.autotune_ncrit    = $(pfield.fmm.autotune_ncrit)")
end

tracer = SFSBlowupTracer()
monitors_with_tracer = (monitors..., tracer)

@time pnl.simulate_warmstart!(systems, wakes, frames, maneuver!, Uinf, t_range;
    name=debug_name,
    path=debug_path,
    restart_path=restart_path,
    restart_name=restart_name,
    restart_step=restart_step,
    body_solvers=body_solvers,
    backend=backend,
    backend_wake=backend_wake,
    monitors=monitors_with_tracer,
    set_Das_eta_kinematic=NaN,
    set_Das_min_kinematic_displacement=set_Das_min_kinematic_displacement,
    verbose=true,
)

println("\n=== Tracer summary ===")
println("  enabled      : $(tracer.enabled)")
println("  seed_step    : $(tracer.seed_step)")
println("  seed_index   : $(tracer.seed_index)")
println("  last |SFS| med = $(tracer.last_sfs_median)")
println("  last |J|   med = $(tracer.last_J_median)")
println("  last |Γ|   med = $(tracer.last_gamma_median)")
