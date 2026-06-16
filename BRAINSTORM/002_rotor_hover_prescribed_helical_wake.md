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

## Root-cause diagnosis - 2026-06-16 (conclusive)

Built `examples/rotor_hover_helical_wake_response.jl`: a frozen-geometry
strength-response sweep that isolates the strength loop from the node-geometry loop.
Freeze the iteration-1 helical wake geometry, set the imposed (uniform-per-column) wake
strength `Gw = alpha * J0` for `alpha in {0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0}` (J0 = the
no-wake TE jump), do ONE body solve per alpha, and record the solved TE jump `J`, CT, and
the signed wake-induced axial velocity at the blade control points.

Result (DJI9443 40_40, 5400 RPM, 1.25-rev wake;
`data/rotor_hover_helical_wake_response/response_sweep.csv`):

| alpha | imposed Gw | solved J | J - Gw | CT(Bern) | wake u_axial |
| --- | --- | --- | --- | --- | --- |
| 0.00 | 0.0      | 0.16997 | 0.16997 | 0.05048 | 0.0    |
| 0.50 | 0.08499 | 0.25454 | 0.16955 | 0.07708 | -2.63  |
| 1.00 | 0.16997 | 0.33911 | 0.16914 | 0.10298 | -5.25  |
| 2.00 | 0.33995 | 0.50824 | 0.16830 | 0.15265 | -10.51 |

Least-squares fit: **`J(Gw) = 0.9951 * Gw + 0.16997`**, i.e. **slope ~= +1, offset ~= +J0**.
`J - Gw` is essentially constant (~0.1697) and the line **never crosses the diagonal
`J = Gw`** -> the relaxation `Gw <- relax*J + (1-relax)*Gw` has **no finite fixed point and
must diverge** (per-step growth `~= relax * offset`, matching the original run's constant
`mean_abs_gamma` increment). CT rises monotonically with imposed strength and the wake
induces axial velocity in `-x` (opposite the slipstream, i.e. upwash-like at the disk),
which *increases* loading -> positive feedback.

**Root cause:** the probe imposes the wake strength externally
(`_steady_aerodynamics!(... update_trailing_edges=false)`) and solves the body freely for
mu, then under-relaxes `Gw` toward the solved jump. This never enforces the Kutta
consistency `Gw = J`. The free solve instead pins the *physical* TE vortex `J - Gw` at the
fixed no-wake value `J0` (set by the kinematic BC), so `J = Gw + J0` for any `Gw` -> the
target stays a constant step ahead forever. This is a non-contractive scheme defect, **not**
missing hover wake physics, wake length, or a loss model. The standard coupled path
converges (steady rigid wake CT 0.0505) because it only sheds the latest jump into the
newest, convecting wake row (`shed_wake!`, `src/FLOWPanel_simulate.jl:843`) rather than
scaling the whole uniform sheet with the current jump.

**Implication for the track:** a usable prescribed/relaxed panel wake here must either
(a) enforce `Gw = J` self-consistently inside the solve (couple the near wake to the Kutta
condition) instead of under-relaxing an externally held sheet, or (b) fix the wake-induced
disk velocity sign (the finite uniform helical sheet currently gives upwash, not downwash).
As written, the relaxation cannot converge by construction. **(See the follow-up section
below: enforcing `Gw = J` alone is NOT sufficient.)**

## Follow-up: would enforcing Kutta consistency fix it? - 2026-06-16

Question considered: solve for the body panel strengths and the *implied* wake strengths
together so the Kutta condition `Gw = J` holds (as in a normal coupled body solve), with
that strength held constant along each helical column. Does that fix convergence?

**Answer: it removes the iteration divergence but converges to a non-physical, hugely
amplified CT. It does not fix the underlying defect.** The convergence pathology and the
physical wrongness are the same issue (open-loop gain `~= +1`), so Kutta consistency cannot
escape it.

Reasoning from the measured affine response `J(Gw) = 0.995*Gw + J0`:

- Two notions of "self-consistent" share the same fixed point:
  - *Picard (what the script does):* hold `Gw` external, free-solve mu, set `Gw = J`,
    repeat. Fixed point needs `J(Gw) = Gw`; but `Gw + J0 = Gw` is impossible -> diverges.
  - *Monolithic (the proposed fix):* fold `Gw = K*mu` into the system,
    `(A_body + B_wake*K)*mu = -rhs`, one linear solve.
- The monolithic solve gives `(I - G)*Gw = J0`, where `G` is the open-loop gain matrix.
  Along the probed (uniform-scaling) mode the sweep measured `G ~= +0.995`, so
  `Gw* = J0 / (1 - 0.995) ~= 200*J0`. The coupled operator `A_body + B_wake*K = M(I + M^-1 N)`
  is **near-singular** (the offending eigenvalue of `M^-1 N` sits at `~= -1`), and its
  "solution" is a ~200x amplified circulation. Extrapolating the sweep's
  `dCT/dGw_mean ~= 0.30` gives `CT ~ O(10)` -- nonsense. So enforcing `Gw = J` converts
  "diverges over iterations" into "converges to nonsense."

**Why gain `~= +1` (the deeper finding):** the clean form `J = Gw + J0` is *doublet-sheet
gauge continuity*, not induced-velocity physics. Pinning the wake's first row exactly on
the TE vertices places a doublet sheet of strength `Gw` immediately behind the TE; the
body's free solve simply rides its usual jump `J0` on top of that pedestal, leaving the
physical bound vortex `J - Gw = J0` unchanged. A correct hover wake should induce *downwash*
at the disk (negative-feedback gain `< 0`, reducing CT toward the BEM/experiment value),
but this finite uniform helical sheet induces *upwash* (`u_axial < 0`, loading-reinforcing).

### Recommendations (fix direction, not yet implemented)

1. **Bring the gain below 1 / restore the correct sign before any coupling scheme matters.**
   This is the actual lever; without it, neither Picard relaxation nor a monolithic solve
   yields a sensible CT.
2. **Do not pin the wake's first row on the TE vertices.** The at-TE placement is what
   produces the gauge-shift / near-unity self-influence. Start the near wake a small step
   downstream with proper Kutta handling (mirror the standard `shed_wake!` near-wake
   treatment), so the wake contributes induced velocity rather than a doublet pedestal.
3. **Verify helix winding / TE-jump sign** so the sheet induces downwash (`u_axial > 0`,
   `+x`) at the disk, not upwash.
4. **Reconsider the "constant strength along the whole helix" assumption.** It makes the
   entire sheet scale 1:1 with the current jump, which is precisely what drives the
   whole-wake gain to `~= 1`. The standard convecting-history wake avoids this because only
   the newest row carries the latest jump.

### Suggested verification before implementing a fix

- Confirm the ~200x amplification empirically: probe a second (non-uniform) strength
  direction to estimate more of `G`'s spectrum, or run a relax-to-consistency / monolithic
  solve on the frozen geometry and watch the strength blow up toward `~200*J0`.
- A/B the near-wake placement: re-run the strength-response sweep with the first wake row
  offset one step downstream of the TE and check whether the slope drops below 1 and the
  sign of `u_axial` flips to `+x`. If so, that isolates the at-TE pinning as the artifact.
