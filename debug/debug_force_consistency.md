# Rotor Hover Force Consistency Debug Log

Date started: 2026-05-28

Context: `examples/rotor_hover_pressure_comparison_replay.jl` currently reports
the following rough behavior for the rotor hover replay:

- `bernoulli` and `laplace_md` are similar and consistent with each other, but
  both underpredict thrust by about 3-4x.
- `laplace_lamb` is much smaller, roughly 20x low.
- `vortex_force` is closer than the pressure monitors, but still has roughly
  50% error.
- `kutta_joukowski` is about 3-4x too large.

Important expectation: Bernoulli pressure is not expected to be reliable with a
vortex-element wake because it depends on a scalar potential description of the
velocity field. The non-Bernoulli force paths should be brought into closer
agreement first.

## Replay and Instrumentation

Primary reproducer:

```bash
julia --project examples/rotor_hover_pressure_comparison_replay.jl
```

Useful replay toggles already present:

- `RUN_STAGE1_LOG=true` writes CSV scale/pressure/force diagnostics to
  `debug/logs/rotor_hover_stage1_*.csv`.
- `RUN_STAGE2_LAPLACE=true` enables all four
  `gradient_mode x acceleration_form` variants:
  `laplace_raw_md`, `laplace_surface_md`, `laplace_raw_lamb`,
  `laplace_surface_lamb`.
- `RECOMPUTE=all` recomputes velocity, potential, velocity gradient, and
  induced vorticity from saved state.
- `STEPS=...` can isolate one saved step before running expensive sweeps.

Initial logging goal: for each candidate fix, append the command, selected
steps, mean/final axial `CT`, and any CSV summary that explains whether the
change moved pressure magnitude, force integration, or normalization.

## Next Steps

Procedure:

- Keep this section as the active handoff point for the investigation.
- After finishing any phase, update this section with exactly the recommended
  next step before ending the session.
- If the user has requested a specific next step, record that user-requested
  next step here even if another investigation path looks more likely.
- When the next step changes, move stale detail into the dated run log below
  instead of leaving multiple competing recommendations here.
- Keep the entry concrete enough that a fresh context can resume immediately:
  include the target method, file/function, command or diagnostic to run, and
  the expected decision criterion.

Current next step:

After changing the Lamb kinetic-energy jump to use tangent relative velocity,
the two-step replay gives `laplace_lamb CT ≈ -0.01658` versus `laplace_md CT ≈
-0.02279`. Next, broaden the replay to more saved steps, with KJ still disabled
by default, and compare `laplace_md`, `laplace_lamb`, `bernoulli`, and
`vortex_force`. Decision criterion: if Lamb remains consistently closer to MD
and no sign flip returns, keep this as the Lamb frame fix and move on to the
remaining common underprediction / pressure integration issue.

## Method Survey

### Bernoulli + ForceMonitor

Code path:

- `PressureBernoulli(...)` populates `body.P`.
- `ForceMonitor(...)` calls `calcfield_F!(body)` and integrates `body.F`.
- `calcfield_F!` applies `F_i = -P_i A_i n_i`, optionally averaging pressure
  across Kutta-condition panel pairs.

Status:

- Similar to `laplace_md`, but 3-4x low.
- This agreement is useful as a scale check, but not a correctness target for
  vortex-wake cases.

Likely issues to keep in mind:

- Bernoulli relies on a scalar potential; vortex-ring wakes do not provide the
  full scalar-potential description needed for this pressure recovery.
- If Bernoulli and Laplace-MD agree, a shared downstream issue is still possible:
  pressure-to-force sign/frame/normalization, body surface normals, or rotor
  `CT` reference.

### Laplace Material Derivative + ForceMonitor

Code path:

- `PressureLaplace(...; acceleration_form=:material_derivative)` builds a sparse
  surface pressure Poisson problem.
- Current default `gradient_mode=:raw_hessian` requests body Hessians only for
  this path.
- `_pressure_rhs_from_edge_material_derivative!` uses a two-point edge form:
  midpoint tangential relative velocity dotted with the velocity jump across the
  panel-pair centerline, plus optional unsteady `velocity_dot`.
- `ForceMonitor` then applies `F_i = -P_i A_i n_i`.

Status:

- Similar to Bernoulli, 3-4x low.

Likely investigation points:

- Check whether the edge RHS is using the intended velocity frame. It currently
  uses `body.velocity` for the convective edge velocity/jump, while
  `_pressure_material_acceleration!` separately forms `u_inertial =
  tangent(body.velocity) + body.velocity_kinematic`. The edge RHS does not use
  that inertial field except through optional `velocity_dot`.
- Compare `_pressure_rhs_from_edge_material_derivative!` to the older
  `_pressure_rhs_from_acceleration!` using `m.acceleration`. If the acceleration
  scratch has the right scale but the edge shortcut does not, the issue is in
  the discrete edge form rather than the pressure solve.
- Run `RUN_STAGE2_LAPLACE=true` and compare `laplace_raw_md` against
  `laplace_surface_md`. A large difference would implicate the raw Hessian or
  missing postprocessed surface-gradient contribution.
- Verify pressure integration separately by injecting a known pressure field and
  checking `calcfield_F!` against direct `-P A n` sums in the rotor frame.

### Laplace Lamb Vector + ForceMonitor

Code path:

- `PressureLaplace(...; acceleration_form=:lamb_vector)` uses
  `_pressure_rhs_from_lamb_vector!`.
- It forms an edge flux from optional `velocity_dot`, the kinetic-energy jump
  in `u_inertial`, and `body.induced_vorticity x u_rel`.
- `monitor_requires_induced_vorticity(m::PressureLaplace)` currently returns
  true for every Lamb-vector variant.

Status:

- Much too small, roughly 20x low.

Likely investigation points:

- The Lamb path likely has a derivative-source mismatch. The cleanup plan in
  `debug/pressure_laplace_cleanup_plan_20260526.md` says raw-Hessian Lamb should
  derive vorticity from the curl of `body.velocity_gradient`, while
  surface-velocity Lamb should use `body.induced_vorticity`. Current code always
  uses `body.induced_vorticity`, independent of `gradient_mode`.
- `body.induced_vorticity` may only contain selected nearfield/extra-output
  vorticity, not the bound surface contribution needed by the Lamb-vector
  identity at panel control points. If so, the vorticity term will be too small.
- Check sign and ordering of the Lamb term. The code computes `omega x u_rel`
  and adds it to the kinetic jump. Confirm this matches the pressure-gradient
  convention used by the assembled Laplacian and `F = -P A n`.
- Compare `laplace_raw_lamb` and `laplace_surface_lamb` under
  `RUN_STAGE2_LAPLACE=true`. If both remain tiny, first inspect
  `body.induced_vorticity` scale in the Stage 1 log or add a vorticity summary.

### SurfaceVorticityForce (`vortex_force`)

Code path:

- `SurfaceVorticityForce` reconstructs `grad_mu` from the panel strength field.
- It computes `kappa = n x grad_mu`.
- It integrates `dF = rho * (V x kappa) * dS` at panel control points using
  `body.velocity` as stored when the monitor runs.
- Optional Kutta correction averages force vectors on paired trailing-edge
  panels.

Status:

- Closer than pressure monitors, but still around 50% error.

Likely investigation points:

- Confirm the sign convention. The docstring says `kappa = -n x grad_s(mu)` but
  the implementation uses `cross(n, grad)`. The surrounding text says this is
  intentional for FLOWPanel's stored strength convention, but this needs a
  one-panel or flat-wing sign test.
- Confirm the velocity frame. This monitor uses `body.velocity`; KJ reconstructs
  `induced + uinf - kinematic` at edge probes. These should be compared at panel
  centers before attributing the remaining 50% to the force formula.
- Check whether `compute_mu_gradient!` is robust on the rotor mesh. Enable
  `gradient_robust=true` and sweep `gradient_ar_threshold`; high aspect-ratio
  root/tip/TE panels may dominate the integrated force.
- Inspect whether the Kutta correction should average panel force vectors,
  pressure-like quantities, or reconstructed sheet strength instead.

### KuttaJoukowskiForce

Code path:

- Builds three edge probes per triangular panel at edge midpoints.
- Evaluates induced velocity from all bodies and wakes, manually adds `uinf`,
  subtracts edge kinematic velocity, and sums
  `F_body = -rho * gamma_panel * (Delta_s x V_edge)` over every panel edge.
- Internal shared edges are therefore included from both adjacent panels with
  opposite edge orientation if mesh orientation is consistent.

Status:

- About 3-4x too large.

Likely investigation points:

- Edge probes lie exactly on vortex-ring filaments. This makes the local
  velocity evaluation sensitive to singular/self-induced treatment and
  `kerneloffset_panel`/`kerneloffset_targets`. First test excluding self/body
  contributions or offsetting probes slightly off the filament.
- Verify internal-edge cancellation. For adjacent panels, the effective edge
  circulation should reduce to the strength jump, not a double-counted absolute
  panel strength. A diagnostic over shared edges should report
  `gamma_left - gamma_right` and compare to the current two-panel sum.
- Check whether bound and wake trailing-edge filaments are both being counted.
  If the wake already carries the shed circulation, summing all body panel edges
  may overcount trailing-edge contributions.
- Sweep `KJ_FMM_EXPANSION_ORDER`, `KJ_FMM_ACCEPTANCE`, and a direct backend for
  one step. If the 3-4x factor persists, the issue is formula/convention rather
  than FMM tolerance.

## Cross-Cutting Checks

1. Normalization:
   Verify `RotorNormalization(rho, 2R, 1)` uses the same sign, diameter, density,
   RPM, and axial dimension as the experimental `CT = -0.072`. Stage 1 already
   logs `ct_force_reference`; compare it to `rho * n^2 * D^4`.

2. Frame convention:
   All reported comparisons use `i_frame=1` and `axial_dimension=1`. Confirm
   the rotor-frame axial axis is the same one used by the experiment and by the
   sign of thrust.

3. Pressure-to-force integration:
   Since Bernoulli and Laplace-MD agree, test `calcfield_F!` independently with
   a synthetic pressure distribution on the saved rotor mesh. This separates
   pressure recovery errors from integration/normalization errors.

4. Velocity consistency:
   Add or inspect diagnostics for `body.velocity`,
   `body.velocity + body.velocity_kinematic`, and KJ edge reconstructed
   velocities. The force monitors should not silently mix inertial, body-relative,
   and rotor-frame velocities.

5. Direct vs FMM:
   Run one saved step with direct backends for force diagnostics where practical.
   This is especially important for KJ probes on or near vortex filaments.

## Suggested Investigation Order

1. Validate normalization, frame axis, and `calcfield_F!` using synthetic
   pressure. This can explain common factors shared by pressure methods.
2. Split `laplace_md` into edge-form RHS versus acceleration-scratch RHS, and
   compare raw-Hessian versus surface-velocity modes.
3. Fix/clarify Lamb-vector vorticity sourcing: raw mode from Hessian curl,
   surface mode from `body.induced_vorticity` plus any required bound-surface
   contribution.
4. Compare `SurfaceVorticityForce` and `KuttaJoukowskiForce` on shared
   panel-center velocity conventions and simple sign tests.
5. Audit KJ edge summation for self-induced velocity, internal-edge effective
   circulation, and trailing-edge/wake double counting.

## Run Log

Append entries below as tests are run.

### 2026-05-28 Initial survey

- Read `CLAUDE.md`, `examples/rotor_hover_pressure_comparison_replay.jl`,
  `src/FLOWPanel_simulate_monitors.jl`, and `src/FLOWPanel_postprocess.jl`.
- No code behavior changed.
- Strongest current suspects:
  - Lamb-vector path uses `body.induced_vorticity` for all modes and may miss
    the bound/raw-Hessian vorticity needed for the identity.
  - Material-derivative pressure may be using a body-relative edge velocity
    where the intended pressure equation needs an inertial/surface velocity.
  - KJ edge probes are exactly on vortex filaments and may overcount shared or
    trailing-edge circulation.
  - Surface-vorticity force needs a sign/frame/gradient robustness check before
    treating its 50% error as physical.

### 2026-05-28 Laplace RHS diagnostic

- Added replay-only `RUN_LAPLACE_RHS_DIAG=true` instrumentation in
  `examples/rotor_hover_pressure_comparison_replay.jl`.
- Diagnostic logs two consecutive replay samples so the finite-difference
  `velocity_dot` buffer has a warmup step; use the second logged row for
  conclusions.
- Command:

```bash
env RUN_BERNOULLI=false RUN_LAPLACE_MD=false RUN_LAPLACE_LAMB=false RUN_VORTEX_FORCE=false RUN_KJ=false RUN_LAPLACE_RHS_DIAG=true STEPS=35,36 RECOMPUTE=all STAGE1_LOG_PATH=debug/logs/rotor_hover_rhs_diag_20260528_steps35_36.csv julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

- Step 36 diagnostic summary from
  `debug/logs/rotor_hover_rhs_diag_20260528_steps35_36.csv`:
  - `lamb_induced_vorticity` RMS / `md_edge` RMS: `0.94194`.
  - `lamb_raw_hessian_curl` RMS / `md_edge` RMS: `0.94192`.
  - Correlation between Lamb RHS and MD edge RHS: about `0.962` for both
    induced-vorticity and raw-Hessian-curl variants.
  - Raw-Hessian material-acceleration scratch is pathological here:
    convective RMS / MD edge RMS is `38.19`.
  - Surface-velocity material-acceleration scratch is close in scale:
    convective RMS / MD edge RMS is `1.14`.
- Interpretation:
  - The old hypothesis that Lamb is tiny because `body.induced_vorticity` is
    much smaller than raw-Hessian curl is not supported by this replay sample.
  - The Lamb-vector edge term is dominated by the kinetic-energy jump; the
    vorticity cross term is tiny by comparison.
  - The raw Hessian acceleration scratch should not be used directly on this
    rotor replay without addressing near-singular/surface-gradient behavior.
- End-to-end pressure monitor check on the same two steps:

```bash
env RUN_BERNOULLI=false RUN_LAPLACE_MD=true RUN_LAPLACE_LAMB=true RUN_VORTEX_FORCE=false RUN_KJ=false RUN_STAGE1_LOG=true STEPS=35,36 RECOMPUTE=all STAGE1_LOG_PATH=debug/logs/rotor_hover_md_lamb_20260528_steps35_36.csv julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

- Step 36 force results:
  - `laplace_md`: final `CT = -0.022794666`.
  - `laplace_lamb`: final `CT = 0.0045635053`.
  - `laplace_lamb` pressure RMS is still comparable to MD (`1.245 q_rotor`
    versus `1.475 q_rotor`), but force cancellation is much stronger
    (`20.03` versus `7.73`) and the axial force changes sign.

### 2026-05-28 Laplace pressure comparison

- Added replay-only `RUN_LAPLACE_PRESSURE_COMPARE=true` instrumentation in
  `examples/rotor_hover_pressure_comparison_replay.jl`.
- Command:

```bash
env RUN_BERNOULLI=false RUN_LAPLACE_MD=false RUN_LAPLACE_LAMB=false RUN_VORTEX_FORCE=false RUN_KJ=false RUN_LAPLACE_PRESSURE_COMPARE=true STEPS=35,36 RECOMPUTE=all STAGE1_LOG_PATH=debug/logs/rotor_hover_pressure_compare_20260528_steps35_36.csv julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

- Step 36 summary from
  `debug/logs/rotor_hover_pressure_compare_20260528_steps35_36.csv`:
  - Pressure correlation: `corr(P_lamb, P_md) = 0.922839`.
  - Best affine fit: `P_lamb ≈ 0.781015 P_md + 100.686 Pa`.
  - Affine residual RMS: `1435.31 Pa = 0.479089 q_rotor`.
  - `RMS(P_lamb - P_md) / RMS(P_md) = 0.393512`.
  - MD pressure force: `CT = -0.0227947`.
  - Lamb pressure force: `CT = 0.00456351`.
  - Affine component force: `CT = -0.017803`.
  - Non-affine residual force: `CT = +0.0223665`.
- Interpretation:
  - Gauge/constant offset is not the main cause. The affine/scaled MD-like part
    still gives negative thrust close to scaled MD.
  - The force reversal is caused by the non-affine residual in the Lamb pressure
    field; that residual is only about `0.48 q_rotor` RMS but has the axial
    loading pattern needed to cancel the scaled MD pressure force.
  - This points back to the assembled Lamb RHS shape/sign/frame, not pressure
    integration alone.

### 2026-05-29 Lamb kinetic-frame split

- Extended `RUN_LAPLACE_PRESSURE_COMPARE=true` to log assembled RHS comparisons
  and component solves from the same sparse Laplacian.
- Step 36 RHS comparison:
  - `corr(b_lamb, b_md) = 0.714662`.
  - Best affine residual in `b_lamb` relative to `b_md`: RMS `11983.1 Pa/m`,
    about `0.54778` of `RMS(b_md)`.
  - `L * (P_lamb - P_md)` matches `b_lamb - b_md` to numerical precision, so
    the bad shape is in the assembled Lamb RHS, not solve/gauge reuse.
- Step 36 component solves:
  - MD edge component: `CT = -0.0227946`.
  - Original Lamb kinetic jump (`u_inertial`): `CT = +0.00441694`.
  - Lamb vorticity cross term: `CT = +0.000146615`.
  - Alternate relative/tangent kinetic jump (`body.velocity` tangent): `CT =
    -0.0167266`.
- Interpretation:
  - The vorticity term is not responsible for the sign flip; it is two orders of
    magnitude smaller in force than the kinetic term.
  - The sign/cancellation problem follows the kinetic-energy velocity frame.
    Using `u_inertial` for `∇(|u|²/2)` while using relative `body.velocity` in
    `ω × u` and in the MD edge form is internally inconsistent.

### 2026-05-29 Lamb relative-velocity fix

- Changed `_pressure_rhs_from_lamb_vector!` in
  `src/FLOWPanel_simulate_monitors.jl` so the kinetic-energy jump uses the same
  tangent-projected relative velocity used by the Lamb cross term.
- Temporarily made `RUN_KJ` default to `false` in
  `examples/rotor_hover_pressure_comparison_replay.jl`; enable it explicitly
  with `RUN_KJ=true` when investigating KJ.
- Updated `test/runtests_unit_postprocess.jl` Lamb RHS projection tests for the
  relative-velocity kinetic term.
- Verification:

```bash
julia --project=. -e 'include("test/runtests_unit_postprocess.jl")'
```

  passed: `178/178`.

- Two-step replay after the fix:

```bash
env RUN_BERNOULLI=false RUN_LAPLACE_MD=true RUN_LAPLACE_LAMB=true RUN_VORTEX_FORCE=false RUN_STAGE1_LOG=true STEPS=35,36 RECOMPUTE=all STAGE1_LOG_PATH=debug/logs/rotor_hover_md_lamb_relative_lamb_20260529_steps35_36.csv julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

  Step 36:
  - `laplace_md CT = -0.022794666`.
  - `laplace_lamb CT = -0.016580195` (previously `+0.0045635053`).

- Pressure comparison after the fix:

```bash
env RUN_BERNOULLI=false RUN_LAPLACE_MD=false RUN_LAPLACE_LAMB=false RUN_VORTEX_FORCE=false RUN_LAPLACE_PRESSURE_COMPARE=true STEPS=35,36 RECOMPUTE=all STAGE1_LOG_PATH=debug/logs/rotor_hover_pressure_compare_relative_lamb_20260529_steps35_36.csv julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

  Step 36:
  - `corr(P_lamb, P_md) = 0.97073` (was `0.922839`).
  - `P_lamb ≈ 0.858735 P_md - 420.82 Pa`.
  - Affine residual RMS: `0.312248 q_rotor` (was `0.479089 q_rotor`).
  - Non-affine residual force: `CT = +0.00299439` (was `+0.0223665`).

### 2026-05-29 `check_replay` NaN/Inf replay failure

- Reproduced with:

```bash
env RUN_BERNOULLI=true RUN_LAPLACE_MD=true RUN_LAPLACE_LAMB=true RUN_VORTEX_FORCE=true STEPS=35,36 RECOMPUTE=auto julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

- The fresh `data/rotor_hover_pressure_comparison_check_replay/*.metadata.toml`
  file contained only one `[[step]]` table, for `i_step = 39`; steps 35 and 36
  were missing.
- Replay therefore failed to recover the rotating frame for steps 35 and 36.
  Before the guard added today, it silently fell back to a default frame with
  `omega = 0`, so `RotorNormalization` divided by zero and all force
  coefficients became infinite.
- Root cause: `_append_metadata_step_toml` appended a separately printed TOML
  document for each step. The parsed unified metadata retained only the last
  step table in this fresh run directory.
- Fix:
  - `_append_metadata_step_toml` now parses the unified metadata, updates or
    appends the requested step in the in-memory `step` array, and rewrites the
    file.
  - `_frames_from_metadata` now errors if unified metadata declares frames but
    the requested step frame state is absent, rather than silently constructing
    a zero-rotation default frame.
- Existing broken `check_replay` output cannot recover steps 35 and 36 from
  that metadata file; regenerate the run after this fix before replaying those
  steps.
- The local `check_replay` directory was repaired manually by reconstructing
  the single rigid frame history from `dt`, constant `omega`, `omega_axis`, and
  zero translation. A backup was written at
  `rotor_hover_pressure_comparison_check_replay.metadata.toml.pre_frame_repair`.
- After repair, the same `STEPS=35,36 RECOMPUTE=auto` replay completed:
  - Step 36 `bernoulli CT = -0.031310211`.
  - Step 36 `laplace_md CT = -0.034248868`.
  - Step 36 `laplace_lamb CT = -0.030918479`.
  - Step 36 `vortex_force CT = -0.074739832`.
