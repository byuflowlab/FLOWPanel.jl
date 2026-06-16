# Rotor Hover CT: Prescribed Helical Wake

## Context

Prior audit has already checked geometry, RPM, density, CT convention, radius/diameter interpretation, and handedness/sign. The active question here is whether a lower-cost prescribed wake can recover the missing rotor loading without the full particle-wake expense.

## Hypothesis

A finite panel helical wake may provide enough downstream vorticity and induced velocity structure to move CT toward the VPM/experiment neighborhood while staying cheaper and more controlled than a fully rolled-up particle wake.

## Proposed Model

- Emit a finite helical wake from each TE strip.
- Keep wake strength constant with wake age for each strip.
- Let strength vary radially based on bound circulation or the TE jump from the solved blade state.
- Iterate helix pitch from azimuthally averaged induced velocity/inflow velocity rather than prescribing a single fixed pitch.
- Keep the wake finite at first so the cost and failure modes remain inspectable.

## Diagnostics

- CT convergence under pitch iteration.
- Radial loading stability across iterations.
- Inflow angle and induced velocity compared with VPM or BEM diagnostics.
- Sensitivity to wake length, azimuthal resolution, and radial TE-strip resolution.

## Acceptance

This track is promising if CT moves toward the VPM/experiment neighborhood without full particle cost and radial loading stabilizes as wake pitch is iterated. If CT remains low or radial loading is unstable under modest wake refinement, the prescribed wake is likely missing essential physics or coupling.

## Caveats

This is a diagnostic model, not a substitute for validating the free-wake path. Keep it scoped as a controlled wake-physics probe unless it consistently explains the discrepancy.

## Progress Log - 2026-06-16

Implemented `examples/rotor_hover_prescribed_helical_wake.jl` as a Bagai/Leishman-style diagnostic free-wake relaxation case. The helical wake is now only the initial condition; subsequent iterations solve the rotor against the current `PanelWake`, relax the TE-jump wake strengths, and relax wake nodes using the solved wake velocity.

Important implementation details:

- `PanelWake(rotor; unsteady_filament=false)` is used so the finite wake has a steady far-edge closure.
- Added controls: `WAKE_RELAX`, `STRENGTH_RELAX`, `WAKE_PSEUDO_DT_FACTOR`, `MAX_NODE_STEP_R`, `GAMMA_TOL`, and `NODE_TOL_R`.
- Added `CT_laplace_lamb` using `PressureLaplace((system,), rho; acceleration_form=:lamb_vector)` followed by its own `ForceMonitor`.
- CSV output now includes Bernoulli CT, Laplace Lamb-vector CT, surface-vorticity CT, Kutta-Joukowski CT, relaxed/candidate strength stats, relative strength change, wake-node motion, wake-row inflow, convergence flag, and stopping reason.
- Fixed a real relaxation issue: the first geometry update added `wake.velocity * pseudo_dt` to every wake node each outer iteration, which has no steady fixed point because a steady hover wake has nonzero convection velocity. The current update relaxes each downstream row toward an age-marched target from the upstream row.
- Fixed a reduced-row test issue: wake pseudo-age `dt` now uses `WAKE_ROWS_PER_REV` rather than the rotor kinematic hard-coded 36 steps/rev.

Runs performed:

1. Smoke, patched short wake:
   `RUN_NAME=/private/tmp/bl_freewake_smoke SAVE_VTK=false MAX_ITER=2 WAKE_REVS=0.25 WAKE_ROWS_PER_REV=12 julia --project examples/rotor_hover_prescribed_helical_wake.jl`

   Result at iteration 2:
   - `CT_bernoulli = 0.06282`
   - `CT_vorticity = 0.07317`
   - `CT_kj = 0.06104`
   - Not converged; only a smoke test.
   - CSV: `/private/tmp/bl_freewake_smoke/iteration_table.csv`

2. Damped short wake with Laplace Lamb-vector, before wake-row `dt` fix:
   `RUN_NAME=/private/tmp/bl_freewake_lamb_damped SAVE_VTK=false MAX_ITER=40 WAKE_REVS=0.25 WAKE_ROWS_PER_REV=12 WAKE_RELAX=0.08 STRENGTH_RELAX=0.10 MAX_NODE_STEP_R=0.015 GAMMA_TOL=5e-3 NODE_TOL_R=2e-3 julia --project examples/rotor_hover_prescribed_helical_wake.jl`

   Result at iteration 40:
   - `CT_bernoulli = 0.21668`
   - `CT_laplace_lamb = 0.20586`
   - `CT_vorticity = 0.23321`
   - `CT_kj = 0.15925`
   - `delta_gamma_rel = 2.209e-02`
   - `max_node_step_R = 0.0012`
   - Not converged; monotone CT/strength growth.
   - CSV: `/private/tmp/bl_freewake_lamb_damped/iteration_table.csv`

3. Damped short wake with Laplace Lamb-vector, after wake-row `dt` fix:
   `RUN_NAME=/private/tmp/bl_freewake_lamb_damped_dtfix SAVE_VTK=false MAX_ITER=40 WAKE_REVS=0.25 WAKE_ROWS_PER_REV=12 WAKE_RELAX=0.08 STRENGTH_RELAX=0.10 MAX_NODE_STEP_R=0.015 GAMMA_TOL=5e-3 NODE_TOL_R=2e-3 julia --project examples/rotor_hover_prescribed_helical_wake.jl`

   Result at iteration 40:
   - `CT_bernoulli = 0.21784`
   - `CT_laplace_lamb = 0.20853`
   - `CT_vorticity = 0.23445`
   - `CT_kj = 0.16084`
   - `delta_gamma_rel = 2.219e-02`
   - `max_node_step_R = 0.0012`
   - Not converged; monotone CT/strength growth persists.
   - CSV: `/private/tmp/bl_freewake_lamb_damped_dtfix/iteration_table.csv`

Verification:

- `julia --project -e 'include("test/runtests_unit_wake.jl")'` passed.
- Parse check of `examples/rotor_hover_prescribed_helical_wake.jl` passed.

Current conclusion:

The diagnostic now predicts CT in the target neighborhood during early iterations, but the relaxed wake does not converge there. For the short 0.25-rev wake tested here, the coupled fixed-point iteration continues increasing circulation and CT well beyond the VPM/BEM/experiment neighborhood. This looks less like a stale reset bug after the node-update fix, and more like the current finite panel-wake relaxation lacks a stabilizing constraint, loss model, or physically appropriate steady free-wake formulation for hover.
