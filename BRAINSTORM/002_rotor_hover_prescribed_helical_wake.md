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

## Near-wake offset / sign A-B results - 2026-06-17

Implemented script-level controls in the diagnostic examples:

- `NEAR_WAKE_MODE=te|das_offset` in `examples/rotor_hover_helical_wake_response.jl`
  and `examples/rotor_hover_prescribed_helical_wake.jl`.
- `NEAR_WAKE_MODE=das_offset` seeds the first helical wake row at `TE + body.Das`,
  matching the production near-wake attachment used by `update_TE!`, rather than exactly
  on the TE vertices.
- `HELIX_WINDING_SIGN` and `WAKE_AXIAL_SIGN` env overrides were added to test whether
  the finite helical sheet's sign/orientation was the dominant issue.

Runs completed:

1. Baseline at-TE placement:
   `RUN_NAME=/private/tmp/helix_response_te SAVE_VTK=false NEAR_WAKE_MODE=te julia --project examples/rotor_hover_helical_wake_response.jl`

   Result:
   - `J(Gw)` fit: slope `0.9951`, offset `1.69972e-01`.
   - Reproduces the old non-contractive failure exactly.
   - CSV: `/private/tmp/helix_response_te/response_sweep.csv`

2. `Das` near-wake offset:
   `RUN_NAME=/private/tmp/helix_response_das SAVE_VTK=false NEAR_WAKE_MODE=das_offset julia --project examples/rotor_hover_helical_wake_response.jl`

   Result:
   - `J(Gw)` fit: slope `0.7918`, offset `1.69972e-01`.
   - The offset weakens the gauge/pedestal behavior, but does **not** produce a useful
     self-consistent solution.
   - The formal fixed point is around `Gw/J0 = 1 / (1 - 0.7918) ~= 4.8`, far beyond the
     externally imposed alpha range that gives plausible CT.
   - Representative rows:

| alpha | imposed Gw | solved J | J - Gw | CT(Bern) | wake u_axial |
| --- | --- | --- | --- | --- | --- |
| 0.25 | 0.04249 | 0.20362 | 0.16113 | 0.06002 | -1.05 |
| 0.50 | 0.08499 | 0.23727 | 0.15228 | 0.06941 | -2.10 |
| 0.75 | 0.12748 | 0.27091 | 0.14343 | 0.07867 | -3.15 |
| 1.00 | 0.16997 | 0.30456 | 0.13459 | 0.08778 | -4.20 |
| 2.00 | 0.33995 | 0.43915 | 0.09920 | 0.12282 | -8.40 |

   Interpretation: alpha `0.25-0.75` can make CT look plausible, but only as an
   externally tuned wake-strength scale. It is not Kutta/self-consistent because solved
   `J` remains much larger than imposed `Gw` in that range.

3. Axial sign flip with offset:
   `RUN_NAME=/private/tmp/helix_response_das_axneg SAVE_VTK=false NEAR_WAKE_MODE=das_offset WAKE_AXIAL_SIGN=-1 ALPHA_LIST=0.0,0.5,1.0,2.0 julia --project examples/rotor_hover_helical_wake_response.jl`

   Result:
   - `J(Gw)` fit: slope `0.8321`, offset `1.69972e-01`.
   - Signed wake axial velocity stayed negative and loading-reinforcing.
   - Not a fix.
   - CSV: `/private/tmp/helix_response_das_axneg/response_sweep.csv`

Current conclusion:

- Moving the first row off the TE helps numerically but does not solve item 002.
- Smaller alpha is not a convergence fix. It is a prescribed wake-strength calibration:
  alpha `~0.5` gives CT near the target, but the coupled iteration would continue driving
  strength upward because the solved TE jump still exceeds the imposed wake strength.
- Do **not** mark item 002 technically complete yet. No converged CT prediction has been
  demonstrated.

Recommended next step:

- Stop tuning alpha/relaxation until the theory is clearer.
- Talk through or review literature for prescribed/free-wake hover formulations that
  specify how near-wake geometry, Kutta enforcement, wake-strength history, and wake
  contraction are coupled without creating a doublet gauge pedestal.
- Specific concepts to check in the literature/theory pass:
  - Bagai/Leishman-style free-wake relaxation: what variables are relaxed, what is held
    fixed, and where the first wake row starts relative to the TE.
  - Rotorcraft prescribed helical wake models with wake contraction: whether circulation is
    held constant along an entire column or convected as a history of prior shed strengths.
  - Panel-method wake Kutta treatments for open lifting surfaces: how the wake doublet
    sheet is coupled to body doublet unknowns without adding a near-singular gauge mode.
  - Whether the wake should be represented as a vortex filament/sheet history rather than
    a constant-strength doublet-panel sheet for this hover diagnostic.
  - Whether a calibrated finite prescribed wake should be explicitly treated as a
    reduced-order model with fitted alpha, rather than advertised as a converged
    self-consistent wake solve.

## One-outer inner-convergence sweep - 2026-06-18

Purpose: hold shed strengths fixed through the inner pitch-geometry loop and allow
exactly one outer strength update, to test whether the wake geometry itself converges
before doing more strength/geometry coupling.

Current corrected defaults used:

- `ITERATION_MODE=nested_pitch`
- `HELIX_WINDING_SIGN=-1`
- `WAKE_AXIAL_SIGN=+1`
- Near wake pinned at `TE + Das` (`NEAR_WAKE_MODE=das_offset` forced by
  `nested_pitch`).
- Fountain hard-stop is based on anchor-relative geometry only: a downstream wake node
  moving upstream of its own local near-wake anchor, or a wake panel crossing that local
  anchor plane. Global raw `x < 0` or negative `u_x` alone is only a warning metric,
  because the two blades and rotated wake columns can produce false positives relative to
  a global `x=0` plane.

Initial visual-check workflow:

`RUN_NAME=/private/tmp/bl_inner_initial SAVE_VTK=false INITIAL_ONLY=true ITERATION_MODE=nested_pitch WAKE_REVS=0.25 WAKE_ROWS_PER_REV=12 julia --project examples/rotor_hover_prescribed_helical_wake.jl`

Result:

- Initial row-2 anchor distance: mean/min/max all `4.1888e-02 R`.
- Upstream row-2 nodes: `0/72`.
- Initial body/wake VTK snapshots written under `/private/tmp/bl_inner_initial/`.

One-outer sweep commands:

1. `RUN_NAME=/private/tmp/bl_inner_006 SAVE_VTK=false ITERATION_MODE=nested_pitch MAX_OUTER_ITER=1 MAX_INNER_ITER=6 WAKE_REVS=0.25 WAKE_ROWS_PER_REV=12 julia --project examples/rotor_hover_prescribed_helical_wake.jl`
2. `RUN_NAME=/private/tmp/bl_inner_012 SAVE_VTK=false ITERATION_MODE=nested_pitch MAX_OUTER_ITER=1 MAX_INNER_ITER=12 WAKE_REVS=0.25 WAKE_ROWS_PER_REV=12 julia --project examples/rotor_hover_prescribed_helical_wake.jl`
3. `RUN_NAME=/private/tmp/bl_inner_024 SAVE_VTK=false ITERATION_MODE=nested_pitch MAX_OUTER_ITER=1 MAX_INNER_ITER=24 WAKE_REVS=0.25 WAKE_ROWS_PER_REV=12 julia --project examples/rotor_hover_prescribed_helical_wake.jl`

Results:

| max inner | stop | inner step | CT(Bern) | CT(Lamb) | max node step/R | min anchor distance/R | anchor crossing | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 6 | `inner-max` | post | 0.05628 | 0.04918 | 0.0045 | 0.02556 | no | geometry not converged |
| 12 | `inner-max` | post | 0.05620 | 0.04920 | 0.0045 | 0.00876 | no | geometry not converged; anchor margin shrinking |
| 24 | `fountain_flow` | 16 | 0.05048 | 0.04427 | 0.0045 | -0.00284 | yes | hard stop before strength update |

Important details:

- The inner convergence criterion is `max_node_step_R < NODE_TOL_R`; with the default
  `NODE_TOL_R=2e-3`, the inner loop never converged. The max node step stayed pinned near
  `0.0045 R`.
- Strengths were held fixed at zero through the inner loops. The `MAX_INNER_ITER=6` and
  `12` runs then performed the single permitted outer update, giving
  `mean_abs_gamma=1.416e-02` from candidate `mean_abs_gamma=1.700e-01`.
- The `MAX_INNER_ITER=24` run found the real "look there first" state before any outer
  strength update: at inner step 16, `2` downstream nodes were upstream of their local
  anchors and `4` wake panels crossed an anchor plane. Diagnostic VTK was written under
  `/private/tmp/bl_inner_024/`:
  `bl_inner_024_iter16_fountain_outer1_inner16_body1*`,
  `bl_inner_024_iter16_fountain_outer1_inner16_wake*`, and
  `bl_inner_024_iter16_fountain_outer1_inner16_wake_filaments*`.
- Negative raw `u_x` and rows with mean signed upstream velocity appeared in all runs
  (`6` downstream rows), but this is not by itself the hard-stop criterion. The 24-inner
  run is the first one that produced anchor-relative geometry crossing.

Conclusion:

Do not run 48-inner or additional outer coupling yet. With corrected signs and `TE + Das`
near-wake placement, the geometry-only pitch relaxation does not converge while strengths
are fixed; it walks the wake toward the local anchor plane and crosses it by inner step 16
on the short 0.25-rev, 12-row/rev wake. The next investigation should inspect the
`/private/tmp/bl_inner_024` fountain state and the geometry update target/capping logic,
not tune outer strength coupling.

## Seeded-strength short-wake follow-up - 2026-06-18

User observation: the ParaView wake strength for `/private/tmp/bl_inner_024` was zero.
That was not a visualization bug. In `nested_pitch`, the first outer iteration holds wake
strength fixed during the inner geometry loop, and the previous diagnostic started from
zero strength; the hard-stop happened before the one permitted outer strength update.

Script controls added:

- `FOUNTAIN_TOL_R`: anchor-relative upstream tolerance before declaring fountain
  crossing. This lets the wake rise slightly upstream of the local anchor plane without
  stopping immediately. The follow-up used `FOUNTAIN_TOL_R=0.05`.
- `INITIAL_STRENGTH_MODE=zero|candidate_latest|candidate_constant`.
  `candidate_constant` first solves the rotor with the initialized wake geometry, reads
  the solved TE-jump candidate strengths, and fills all active wake rows with those
  strengths before the nested pitch loop.

Run:

`RUN_NAME=/private/tmp/bl_inner_024_seed_constant_tol05 SAVE_VTK=false FOUNTAIN_TOL_R=0.05 INITIAL_STRENGTH_MODE=candidate_constant ITERATION_MODE=nested_pitch MAX_OUTER_ITER=1 MAX_INNER_ITER=24 WAKE_REVS=0.25 WAKE_ROWS_PER_REV=12 julia --project examples/rotor_hover_prescribed_helical_wake.jl`

Result:

- Initial strength seeding was nonzero: candidate mean `1.700e-01`, wake mean
  `1.700e-01`, wake max `2.295e-01`.
- The previous zero-strength fountain stop at inner 16 disappeared under the seeded wake
  plus `0.05 R` fountain tolerance.
- Inner-loop CT was in the target range and nearly flat but still rising slowly:
  `CT(Bern) 0.07537` at inner 1 -> `0.07552` at inner 24.
- Laplace Lamb CT rose from `0.06097` to `0.06148`; row inflow ratio rose from `0.0251`
  to `0.0265`.
- `max_node_step_R` stayed pinned at `0.0045`, above `NODE_TOL_R=0.002`, so the geometry
  did **not** converge.
- After the single outer strength update, final row:
  `CT(Bern)=0.07827`, `CT(Lamb)=0.06371`, `CT(vort)=0.09083`, `CT(KJ)=0.05548`,
  `delta_gamma_rel=0.1367`, `mean_abs_gamma=0.1777`, candidate mean `0.2629`.
- CSV: `/private/tmp/bl_inner_024_seed_constant_tol05/iteration_table.csv`.

Conclusion:

Starting with a nonzero wake is materially better than the zero-strength inner loop: it
provides immediate induced velocity, removes the obvious fountain behavior in this short
test, and gives a plausible CT range. It is still not a converged simulation because the
geometry residual does not contract below tolerance and the one allowed strength update
still changes strengths by `13.7%`.

Next staged test after this pause:

Use the same seeded-strength setup and `FOUNTAIN_TOL_R=0.05`, but increase wake extent
before adding more outer iterations. A reasonable next run is:

`RUN_NAME=/private/tmp/bl_inner_long_seed_constant_tol05 SAVE_VTK=false FOUNTAIN_TOL_R=0.05 INITIAL_STRENGTH_MODE=candidate_constant ITERATION_MODE=nested_pitch MAX_OUTER_ITER=1 MAX_INNER_ITER=24 WAKE_REVS=0.75 WAKE_ROWS_PER_REV=24 julia --project examples/rotor_hover_prescribed_helical_wake.jl`

Primary checks for the longer wake:

- Does `max_node_step_R` begin to contract rather than staying pinned at `0.0045`?
- Does CT settle instead of drifting upward?
- Does row inflow remain physically downstream without relying only on a large fountain
  tolerance?
- If one outer update remains stable, then try `MAX_OUTER_ITER=2` with the same longer
  wake before changing relaxation constants.

## Axial advance-ratio pitch seed - 2026-06-18

Question: does a nonzero axial advance ratio make the initialized helical wake closer to a
streamline and easier to converge, while keeping the wake axisymmetric? This branch tests
axial inflow only, not edgewise forward flight.

Implementation:

- Added `AXIAL_ADVANCE_RATIO`, default `0.0`, to
  `examples/rotor_hover_prescribed_helical_wake.jl`.
- The imposed freestream is now axial:
  `Vinf = axial_wake_sign * AXIAL_ADVANCE_RATIO * tip_speed` along the rotor axial
  dimension.
- The initial helix pitch speed is seeded from freestream convection plus induced inflow:
  `(AXIAL_ADVANCE_RATIO + INITIAL_INFLOW) * tip_speed`.
- Startup output now reports axial advance ratio, freestream axial speed, induced seed
  ratio, and total initial helix pitch ratio.
- Derivation added in `docs/src/prescribed_helical_wake.md` and linked into
  `docs/make.jl`.

For the first diagnostic value `AXIAL_ADVANCE_RATIO=0.10` with default
`INITIAL_INFLOW=0.08`, the total initial helix pitch ratio is `0.18`. With
`WAKE_ROWS_PER_REV=12`, the row-2 anchor distance starts at
`2π * 0.18 / 12 = 9.4248e-02 R`, matching the example startup diagnostic.

Smoke run:

`RUN_NAME=/private/tmp/fp_axial_smoke AXIAL_ADVANCE_RATIO=0.10 MAX_OUTER_ITER=1 MAX_INNER_ITER=1 WAKE_REVS=0.25 WAKE_ROWS_PER_REV=12 SAVE_VTK=false julia --project examples/rotor_hover_prescribed_helical_wake.jl`

Result:

- Startup diagnostics were correct: axial advance ratio `0.1000`, freestream axial speed
  `6.7293 m/s` (`0.1000 tip`), induced seed `0.0800`, total pitch ratio `0.1800`.
- Initial row-2 anchor distance: mean/min/max all `9.4248e-02 R`; upstream nodes `0/72`.
- No fountain crossing in the single inner iteration.
- Single-iteration residual metrics remained poor:
  `max_target_residual_R=0.5225`, `capped_node_fraction=1.000`,
  `max_applied_step_R=0.0045`.
- Final post-update row: `CT(Bern)=0.03044`, `CT(Lamb)=0.02440`,
  `CT(vort)=0.03703`, `CT(KJ)=0.03197`, `delta_gamma_rel=0.1336`.
- CSVs: `/private/tmp/fp_axial_smoke/iteration_table.csv` and
  `/private/tmp/fp_axial_smoke/row_residual_table.csv`.

Cheap convergence probe started:

`RUN_NAME=/private/tmp/fp_axial_probe AXIAL_ADVANCE_RATIO=0.10 MAX_OUTER_ITER=1 MAX_INNER_ITER=24 WAKE_REVS=0.25 WAKE_ROWS_PER_REV=12 julia --project examples/rotor_hover_prescribed_helical_wake.jl`

Partial findings before pause request:

- Through inner iteration 14, no fountain crossing occurred.
- `max_target_residual_R` contracted only slightly, from `0.5225` at inner 1 to `0.5153`
  at inner 14.
- `capped_node_fraction` stayed pinned at `1.000`, so every relaxed node was still
  step-limited.
- `max_applied_step_R` stayed pinned at `0.0045`.
- CT stayed nearly flat and low for this axial case: `CT(Bern)≈0.0293`,
  `CT(Lamb)≈0.0230`, `CT(KJ)≈0.0310`.

Current read:

The axial pitch seed fixes the initial geometric pitch accounting and avoids immediate
fountain behavior, but it has not yet shown a material convergence improvement. The
longer-wake diagnostic should wait until the 24-inner cheap probe reaches its natural end
and the final `row_residual_table.csv` confirms whether any rows drop below full capping.

## Reduced-row seeded probes and debug VTK mode - 2026-06-18

The latest reduced-row probes tested whether making the short seeded wake even smaller
helps the capped Picard geometry update converge. It does not.

Commands recorded:

1. One wake row:
   `RUN_NAME=/private/tmp/bl_rows1_seed_constant_tol05 SAVE_VTK=false FOUNTAIN_TOL_R=0.05 INITIAL_STRENGTH_MODE=candidate_constant ITERATION_MODE=nested_pitch MAX_OUTER_ITER=1 MAX_INNER_ITER=24 WAKE_REVS=0.08333333333333333 WAKE_ROWS_PER_REV=12 julia --project examples/rotor_hover_prescribed_helical_wake.jl`

   Final row:
   - `CT(Bern)=0.07990`
   - `CT(Lamb)=0.06808`
   - `max_target_residual_R=0.4168`
   - `capped_node_fraction=1.000`
   - `delta_gamma_rel=0.1398`
   - Stop reason: `inner-max`
   - CSV: `/private/tmp/bl_rows1_seed_constant_tol05/iteration_table.csv`

2. Two wake rows:
   `RUN_NAME=/private/tmp/bl_rows2_seed_constant_tol05 SAVE_VTK=false FOUNTAIN_TOL_R=0.05 INITIAL_STRENGTH_MODE=candidate_constant ITERATION_MODE=nested_pitch MAX_OUTER_ITER=1 MAX_INNER_ITER=24 WAKE_REVS=0.16666666666666666 WAKE_ROWS_PER_REV=12 julia --project examples/rotor_hover_prescribed_helical_wake.jl`

   Final row:
   - `CT(Bern)=0.07865`
   - `CT(Lamb)=0.06507`
   - `max_target_residual_R=0.5108`
   - `capped_node_fraction=1.000`
   - `delta_gamma_rel=0.1379`
   - Stop reason: `inner-max`
   - CSV: `/private/tmp/bl_rows2_seed_constant_tol05/iteration_table.csv`

Conclusion:

Reducing to one or two wake rows does not fix convergence. Seeded short wakes can produce
plausible CT, but the current Picard geometry update remains cap-dominated (`100%` of
relaxed nodes step-limited) and the strength update is still large after the inner loop.
Do not run longer wakes or more outer iterations until the geometry update is changed.

Diagnostic harness update:

- Added `DEBUG_VTK_EVERY_ITER=true|false` to
  `examples/rotor_hover_prescribed_helical_wake.jl`, default `false`.
- When enabled in `nested_pitch`, the script writes one body/wake ParaView time series per
  outer iteration. Open `<run_name>_body1_iter0001.pvd`,
  `<run_name>_wake_iter0001.pvd`, and `<run_name>_wake_filaments_iter0001.pvd` for outer
  iteration 1. Timestep `0` is the original state at the beginning of that outer iteration,
  before any inner wake-shape update; timesteps `1..N` are the fixed-strength inner
  wake-shape states. The final solved state after the strength update is written
  separately without an iteration suffix as
  `<run_name>_body1.pvd`, `<run_name>_wake.pvd`, and `<run_name>_wake_filaments.pvd`.
  Coupled-relax debug output uses `<run_name>_body1_iter.pvd` /
  `<run_name>_wake_iter.pvd`. Fountain/divergence dumps retain their special suffixes.
- `SAVE_VTK` remains the lower-volume output mode. `DEBUG_VTK_EVERY_ITER=true` is an
  explicit high-volume diagnostic override and prints a startup warning.

Smoke/parse check:

`RUN_NAME=/private/tmp/bl_debug_smoke SAVE_VTK=false DEBUG_VTK_EVERY_ITER=true MAX_OUTER_ITER=1 MAX_INNER_ITER=2 WAKE_REVS=0.08333333333333333 WAKE_ROWS_PER_REV=12 julia --project examples/rotor_hover_prescribed_helical_wake.jl`

Expected output includes per-iteration body/wake VTK files for inner 1, inner 2, and the
post-strength update under `/private/tmp/bl_debug_smoke/`.

## Multi-phase wake-solver implementation plan - 2026-06-18

Goal: replace the currently capped Picard-style wake geometry relaxation with a staged,
testable wake-shape solver path. The first implementation target is still the diagnostic
example, not a general public FLOWPanel API.

### Phase 0: lock the residual and safeguards

- [ ] Replace the nested-pitch inner loop with frozen-field streamline integration from
  fixed TE+Das anchors.
- [ ] Add frame-relative effective velocity using the frozen ReferenceFrame state and
  document the transform.
- [ ] Verify the new inner loop leaves wake/body strengths unchanged during geometry
  iterations.
- [ ] Run short hover and axial smoke checks before revisiting outer strength coupling.
- [ ] Define the exact residual to solve first: VSPAero-style arclength alignment,
  age-marched trapezoidal residual, or both behind a solver option.
- [ ] Decide which quantities are frozen inside the inner solve: initially `q`, `Γ_w`,
  and `Δs_e`.
- [ ] Keep row-1 anchoring, fountain checks, anchor-crossing checks, step caps, and
  residual-norm rejection as hard globalization constraints.
- [ ] Add diagnostics needed by all later phases: residual norm, accepted step norm,
  rejected line-search count, GMRES iterations, and per-row residual summaries.

### Phase 1: one-panel-at-a-time continuation

- [ ] Implement a continuation mode that starts with a single wake panel row/age block,
  solves it, then adds the next panel row/block and resolves.
- [ ] Preserve the already-solved upstream rows as the initial condition for each larger
  wake problem.
- [ ] Record whether residual contraction improves as wake extent grows compared with
  solving the full wake at once.
- [ ] Use this continuation path as the default first test for any new nonlinear solver,
  because it should reduce the initial residual and improve convergence robustness.

### Phase 2: abstract the wake solver interface

- [ ] Factor the current Picard relaxation into a solver backend with a shared call shape
  and shared return diagnostics.
- [ ] Add selectable wake solver modes, initially `:picard`, `:continuation_picard`, and
  later `:newton_gmres`.
- [ ] Keep the existing Picard path behaviorally unchanged as the default.
- [ ] Ensure CSV/logging code consumes solver diagnostics generically instead of depending
  on Picard-only field names.

### Phase 3a: derive and test the analytic Jv kernel math

- [ ] Derive the analytic Jacobian-vector product for the exact Vatistas finite-core
  bound-vortex kernel used by `_bound_vortex_velocity`, not only the simplified theory
  kernel.
- [ ] Include derivatives with respect to target point and both segment endpoints.
- [ ] Verify the segment-level derivative against finite differencing over representative
  endpoint/target configurations, including near-core but non-singular cases.
- [ ] Document the final derivation in `docs/src/prescribed_helical_wake.md`.

### Phase 3b: implement residual and analytic Jv

- [ ] Implement wake-coordinate packing/unpacking for movable wake nodes only.
- [ ] Implement residual evaluation for the chosen first residual while preserving pinned
  near-wake rows.
- [ ] Implement analytic `Jv` for the residual, including normalized-velocity projection
  and wake-segment endpoint perturbations.
- [ ] Verify full residual-level `Jv` against finite differencing and, where practical,
  AD on small synthetic wakes.
- [ ] Confirm analytic `Jv` matches the same velocity conventions as current `PanelWake`
  panel and final-filament influence.

### Phase 4: non-preconditioned matrix-free GMRES

- [ ] Wrap `Jv!` in `LinearOperators.LinearOperator`.
- [ ] Solve `J δX = -r` with `Krylov.jl` GMRES and no preconditioner.
- [ ] Add damped Newton updates with line search/trust-region step caps.
- [ ] Verify the nonlinear residual can be driven down arbitrarily on small synthetic
  wakes where the residual has a known consistent solution.
- [ ] Verify residual drops on the rotor diagnostic before using CT as a success metric.

### Phase 5: preconditioner selection and implementation

- [ ] Evaluate preconditioner options after the unpreconditioned path works:
  block-Jacobi by wake row, local edge/segment block preconditioner, arclength-only
  diagonal scaling, or approximate lower-triangular age-marching preconditioner.
- [ ] Choose one preconditioner based on observed GMRES iteration counts and user
  preference.
- [ ] Implement only the selected preconditioner first.
- [ ] Verify it reduces GMRES iterations without changing the accepted nonlinear solution.

### Phase 6: end-to-end solver testing

- [ ] Run the one-panel-at-a-time continuation with `:newton_gmres` on short wakes.
- [ ] Run the seeded-strength short-wake case and compare residual contraction, cap
  fraction, CT drift, and anchor/fountain behavior against Picard.
- [ ] Increase wake extent only after short-wake residual contraction is demonstrated.
- [ ] Test at least one axial advance-ratio case and one hover case.
- [ ] Treat success as residual contraction plus stable accepted updates first; CT in the
  VPM/BEM/experiment range is a secondary validation, not proof of solver correctness.

### Missing risks/checks to keep explicit

- [ ] Globalization is as important as GMRES: line search, step caps, and rejection checks
  must be implemented before judging convergence.
- [ ] The analytic Jv must differentiate the code kernel exactly; a theory/code mismatch
  would make GMRES misleading.
- [ ] Strength coupling remains a separate outer-loop problem. Do not mix `δΓ_w` into the
  first Newton unknown vector until geometry-only convergence is proven.
- [ ] The current wake-sign/upwash issue must remain visible in diagnostics; a better
  nonlinear solver should not hide a physically wrong residual.
- [ ] Add docs whenever equations or solver assumptions change, so item 002 remains
  reproducible.
