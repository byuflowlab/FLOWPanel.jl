# Rotor Hover Debug - 2026-05-26

Purpose: stage the next rotor-hover thrust-coefficient discrepancy experiments from
May 26, 2026. This note is a planning document only; no experiments were run as part
of creating it.

Operational rule: after each experiment stage, append the command, output location,
key values, and interpretation to this file before moving to the next stage.

## Observed CT Split

- `KuttaJoukowskiForce ≈ 0.03`
- `PressureLaplace(lamb) ≈ 0.003`
- `PressureBernoulli` and `PressureLaplace(material_derivative) ≈ 5e-5`
- Experiment target `≈ 0.07`

## Working Hypothesis

Pressure reconstruction is losing rotor-relative velocity, acceleration, and
surface-vorticity content while the circulation/Kutta-Joukowski path is closer to
the expected thrust scale.

## Experiment Stages

### Experiment Checklist

- [ ] Stage 1: Scale and frame sanity log
- [x] Stage 2: `PressureLaplace` variant comparison
- [x] Stage 3: Rotor-relative Bernoulli diagnostic
- [ ] Stage 4: Kutta-Joukowski edge-circulation check
- [ ] Stage 5: Independent kernel-offset sweep

### Stage 1: Scale and Frame Sanity Log

Log scale and reference-frame sanity values for:

- [x] panel strength ranges and norms
- [x] absolute, inertial, body-frame, and rotor-relative velocity scales
- [ ] body kinematics and acceleration terms
- [x] angular velocity magnitude and direction
- [x] pressure ranges from each pressure model
- [x] integrated force and normalized thrust coefficient

Compare the recovered pressure scale against `0.5ρ(ΩR)^2`. In hover, any diagnostic
that effectively normalizes against only `norm(uinf)` is expected to hide the rotor
speed scale.

#### Stage 1 Execution Plan

1. Add a lightweight debug logger to the rotor-hover pressure comparison replay
   workflow, gated behind a local flag so normal replay output stays unchanged.
2. Replay only the final few saved timesteps from
   `examples/rotor_hover_pressure_comparison_replay.jl` to avoid rerunning the
   full rotor-hover simulation.
3. Log one compact row per stored timestep with frame state, `Ω`, `R`, `ΩR`,
   `0.5ρ(ΩR)^2`, and the force normalization used for CT.
4. For each body at the same timestep, log panel-strength extrema, mean absolute
   strength, and representative quantiles to catch one-panel outliers without dumping
   full arrays.
5. Log velocity-scale summaries separately for freestream, kinematic velocity,
   induced velocity, total panel velocity, and rotor-relative velocity. Include min,
   max, mean, and RMS norms.
6. Log pressure summaries for `PressureBernoulli`, `PressureLaplace(lamb)`, and
   `PressureLaplace(material_derivative)` using the same panels and normalization.
7. Log integrated dimensional force and normalized CT for each pressure path plus
   `KuttaJoukowskiForce`.
8. Compare each pressure range against `0.5ρ(ΩR)^2` and flag paths whose pressure
   scale is closer to `0.5ρ|uinf|^2` than the rotor-relative scale.
9. Save the output under `debug/logs/` with the date and experiment label, then append
   the key CT and pressure-scale ratios back to this note.

#### Stage 1 Implementation Progress

- [x] Added an opt-in replay logger to
  `examples/rotor_hover_pressure_comparison_replay.jl` behind `RUN_STAGE1_LOG=true`.
- [x] Logged frame scale, velocity scale, body strength summaries, pressure summaries,
  dimensional axial force, and monitor CT for replayed samples.
- [x] Verified the logger on the final saved timestep with:

```bash
env RUN_STAGE1_LOG=true NSTEPS_REPLAY=1 RUN_KJ=false STAGE1_LOG_PATH=debug/logs/rotor_hover_stage1_smoke.csv julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

Smoke output was written to `debug/logs/rotor_hover_stage1_smoke.csv`.

- [x] Run the full Stage 1 replay on the final few saved timesteps with
  `RUN_KJ=true` so Kutta-Joukowski CT is included in the CSV.
- [x] Append the final Stage 1 pressure-scale ratios and CT comparison to this note.

Final-four replay command:

```bash
env RUN_STAGE1_LOG=true NSTEPS_REPLAY=4 RUN_KJ=true STAGE1_LOG_PATH=debug/logs/rotor_hover_stage1_20260526_final4.csv julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

Final-four replay output was written to
`debug/logs/rotor_hover_stage1_20260526_final4.csv` with 112 data rows.

Key values:

- `0.5ρ(ΩR)^2 = 3424.204964 Pa`
- rotor-relative velocity RMS `≈ 43.116 m/s`, max `≈ 74.769 m/s`
- `PressureBernoulli` CT: first replayed sample `-0.07839`, final sample
  `-7.75e-6`
- `PressureLaplace(material_derivative)` CT: `≈ -5.01e-5`
- `PressureLaplace(lamb)` CT: `≈ -0.00305`
- `KuttaJoukowskiForce` CT: `≈ -0.0306` to `-0.0315`
- pressure RMS over `0.5ρ(ΩR)^2`: Bernoulli `≈ 0.46`, Laplace material
  derivative `≈ 0.42`, Laplace Lamb `≈ 0.066`

Interpretation: pressure magnitudes can reach rotor dynamic-pressure scale in the
replayed fields, but integrated CT remains strongly split by reconstruction path.
The Lamb-vector pressure field has much smaller RMS pressure scale than Bernoulli
or material-derivative pressure, while its integrated CT is still two orders above
the material-derivative force and one order below Kutta-Joukowski.

Reasonableness assessment:

- The velocity and frame scales look reasonable for the replayed hover state:
  `ΩR ≈ 74.77 m/s`, rotor-relative velocity RMS `≈ 43.1 m/s`, and max velocity
  `≈ ΩR`.
- The pressure magnitudes are not simply missing the rotor dynamic-pressure scale.
  Bernoulli and material-derivative Laplace both have RMS pressure of order
  `0.4-0.5 × 0.5ρ(ΩR)^2`.
- The unreasonable part is the force outcome. Material-derivative Laplace has a
  plausible pressure RMS scale but integrates to only `CT ≈ 5e-5`, which points to
  cancellation, sign/orientation structure, or pressure placement over the blade
  surface rather than a simple pressure-amplitude loss.
- Lamb-vector Laplace has a much smaller pressure RMS (`≈ 0.066 q_rotor`) but a
  larger integrated CT (`≈ 0.003`), so the RHS/distribution shape appears more
  important than raw pressure amplitude.
- The first replayed Bernoulli sample (`CT ≈ -0.078`) is suspicious because later
  Bernoulli samples collapse to `O(1e-5)`. Treat it as likely affected by replay
  history initialization for the unsteady `phi_dot` term until checked.

Recommended next diagnostic: inspect panelwise pressure-force correlation and
axial-force cancellation for each pressure path, especially
`PressureLaplace(material_derivative)`.

Panelwise pressure-force diagnostic:

Added panelwise force-cancellation rows to the Stage 1 replay logger. For each
pressure path, the logger now computes per-panel axial CT from
`F_panel = -P_panel A_panel n_panel`, rotated into the rotor frame, then records
positive axial sum, negative axial sum, absolute axial sum, net sum, cancellation
ratio, and pressure/force correlations.

Command:

```bash
env RUN_STAGE1_LOG=true NSTEPS_REPLAY=4 RUN_KJ=true STAGE1_LOG_PATH=debug/logs/rotor_hover_stage1_panelwise_20260526_final4.csv julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

Output: `debug/logs/rotor_hover_stage1_panelwise_20260526_final4.csv` with 208
data rows.

Key final-sample values:

- `PressureBernoulli`: `CT+ ≈ 0.10748`, `CT- ≈ -0.10749`,
  net `CT ≈ -7.75e-6`, cancellation ratio `≈ 2.77e4`
- `PressureLaplace(material_derivative)`: `CT+ ≈ 0.09483`,
  `CT- ≈ -0.09488`, net `CT ≈ -5.01e-5`, cancellation ratio `≈ 3.79e3`
- `PressureLaplace(lamb)`: `CT+ ≈ 0.00454`, `CT- ≈ -0.00759`,
  net `CT ≈ -0.00305`, cancellation ratio `≈ 3.98`
- `KuttaJoukowskiForce`: monitor CT `≈ -0.03145`
- Absolute pressure to absolute panel-CT correlation is moderate for all pressure
  paths: Bernoulli `≈ 0.559`, material-derivative Laplace `≈ 0.570`, Lamb Laplace
  `≈ 0.524`
- Signed pressure to signed panel-CT correlation is near zero for the two nearly
  cancelling paths: Bernoulli `≈ 2.3e-4`, material-derivative Laplace
  `≈ -0.0029`; Lamb Laplace is still weak but larger at `≈ 0.060`

Interpretation: this confirms the Stage 1 suspicion. The material-derivative
Laplace path is not producing a tiny force because pressure is too small; it is
producing nearly equal positive and negative axial panel-force contributions. The
same cancellation pattern also explains the later Bernoulli samples. Lamb-vector
Laplace has much smaller panelwise axial force content overall and far less
cancellation, but it is still an order of magnitude below Kutta-Joukowski.

Next useful check: map the largest positive and negative `CT_panel` contributors
back to panel position, radius, side of blade, and pressure sign. That will show
whether cancellation is physical upper/lower-side cancellation, a pressure sign
pattern issue, or a frame/normal orientation issue.

### Stage 2: PressureLaplace Variant Comparison

Compare the `PressureLaplace` variants side by side:

- [x] `gradient_mode=:surface_velocity`
- [x] `gradient_mode=:raw_hessian`
- [x] Lamb-vector RHS form
- [x] material-derivative RHS form

The immediate goal is to identify whether the CT loss appears when the RHS is built,
when the velocity gradient is reconstructed, or when the pressure field is integrated.

#### Stage 2 Results

Added an opt-in Stage 2 replay mode behind `RUN_STAGE2_LAPLACE=true` in
`examples/rotor_hover_pressure_comparison_replay.jl`. The mode runs four
`PressureLaplace` variants in the same replay pass:

- `laplace_raw_md`: `gradient_mode=:raw_hessian`,
  `acceleration_form=:material_derivative`
- `laplace_surface_md`: `gradient_mode=:surface_velocity`,
  `acceleration_form=:material_derivative`
- `laplace_raw_lamb`: `gradient_mode=:raw_hessian`,
  `acceleration_form=:lamb_vector`
- `laplace_surface_lamb`: `gradient_mode=:surface_velocity`,
  `acceleration_form=:lamb_vector`

Smoke command:

```bash
env RUN_STAGE2_LAPLACE=true RUN_STAGE1_LOG=true RUN_BERNOULLI=false RUN_LAPLACE_MD=false RUN_LAPLACE_LAMB=false RUN_KJ=false NSTEPS_REPLAY=1 STAGE1_LOG_PATH=debug/logs/rotor_hover_stage2_smoke.csv julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

Smoke output was written to `debug/logs/rotor_hover_stage2_smoke.csv`.

Final-four replay command:

```bash
env RUN_STAGE2_LAPLACE=true RUN_STAGE1_LOG=true RUN_BERNOULLI=false RUN_LAPLACE_MD=false RUN_LAPLACE_LAMB=false RUN_KJ=true NSTEPS_REPLAY=4 STAGE1_LOG_PATH=debug/logs/rotor_hover_stage2_20260526_final4.csv julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

Output: `debug/logs/rotor_hover_stage2_20260526_final4.csv`.

Final-four CT summary:

- `laplace_raw_md`: mean/final `CT = -5.0069342e-5`
- `laplace_surface_md`: mean/final `CT = -5.0069342e-5`
- `laplace_raw_lamb`: mean/final `CT = -0.0030493039`
- `laplace_surface_lamb`: mean/final `CT = -0.0030493039`
- `KuttaJoukowskiForce`: mean `CT = -0.030876734`, final
  `CT = -0.031452966`

Final-sample pressure and cancellation values:

- Material-derivative variants: pressure RMS `≈ 0.4220 q_rotor`,
  `CT+ ≈ 0.094831`, `CT- ≈ -0.094881`, net
  `CT ≈ -5.0069e-5`, cancellation ratio `≈ 3789`
- Lamb-vector variants: pressure RMS `≈ 0.06582 q_rotor`,
  `CT+ ≈ 0.004544`, `CT- ≈ -0.007594`, net
  `CT ≈ -0.003049`, cancellation ratio `≈ 3.98`

Interpretation: in the replayed hover state, `gradient_mode` does not affect the
pressure field or integrated CT for either RHS form. The split appears at the RHS
choice and/or the force distribution it produces: material-derivative pressure has
rotor-scale RMS but nearly perfect axial cancellation, while Lamb-vector pressure
has much lower RMS and less cancellation but remains about one order of magnitude
below Kutta-Joukowski. This points away from the raw-Hessian velocity-gradient
reconstruction as the source of the CT loss in replay.

Follow-up cleanup plan: `debug/pressure_laplace_cleanup_plan_20260526.md`.

### Stage 3: Rotor-Relative Bernoulli Diagnostic

Run a Bernoulli diagnostic using local rotor-relative speed rather than only
`norm(uinf)`.

This should isolate whether the hover pressure scale is being lost because the
Bernoulli path is using a near-zero freestream scale where the physical pressure
variation should be set by `ΩR` and local induced velocity.

#### Stage 3 Results

Added two opt-in replay pressure paths behind `RUN_STAGE3_BERNOULLI=true`:

- `bernoulli_rotor_relative`: steady Bernoulli using replayed rotor-relative
  panel velocity, equivalent to `P = -0.5ρ|u_rel|^2` with no unsteady `phi_dot`.
- `bernoulli_local_omega`: local blade-speed reference
  `P = 0.5ρ(|u_kinematic|^2 - |u_rel|^2)`.

Command:

```bash
env RUN_STAGE3_BERNOULLI=true RUN_STAGE1_LOG=true NSTEPS_REPLAY=4 RUN_KJ=true STAGE1_LOG_PATH=debug/logs/rotor_hover_stage3_20260526_final4.csv julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

Output: `debug/logs/rotor_hover_stage3_20260526_final4.csv` with 312 data rows.

Final-four CT summary:

- `PressureBernoulli` with unsteady replay `phi_dot`: mean `CT ≈ -0.0196`,
  final `CT ≈ -7.75e-6`
- `bernoulli_rotor_relative`: `CT ≈ -2.78e-6` for all four samples
- `bernoulli_local_omega`: `CT = 0` for all four samples
- `PressureLaplace(material_derivative)`: `CT ≈ -5.01e-5`
- `PressureLaplace(lamb)`: `CT ≈ -0.00305`
- `KuttaJoukowskiForce`: mean `CT ≈ -0.0309`, final `CT ≈ -0.0315`

Scale and cancellation details:

- `bernoulli_rotor_relative` pressure RMS remains rotor-scale:
  `≈ 0.463 × 0.5ρ(ΩR)^2`, matching the ordinary Bernoulli pressure magnitude.
- Despite rotor-scale pressure, `bernoulli_rotor_relative` has
  `CT+ ≈ 0.10749`, `CT- ≈ -0.10749`, net `CT ≈ -2.78e-6`, and cancellation
  ratio `≈ 7.7e4`.
- `bernoulli_local_omega` is identically zero in replay because
  `|u_rel| == |u_kinematic|` panelwise for these saved states. This matches the
  Stage 1 observation that `body.velocity + body.velocity_kinematic` has zero
  RMS in replay.

Interpretation: Stage 3 rules out a simple "Bernoulli used `norm(uinf)` instead
of rotor-relative speed" explanation for the low hover CT in the replayed state.
The steady rotor-relative Bernoulli pressure has the expected rotor dynamic
pressure scale, but its axial force cancels almost exactly. The local-Ω reference
diagnostic shows that the replayed body velocity is dominated by rigid-body
kinematics with no residual inertial velocity, so this replay path cannot expose
induced/local inflow corrections through `|u_kinematic|^2 - |u_rel|^2`.

Next useful check remains geometric localization of the cancelling panels: identify
the largest positive and negative `CT_panel` contributors by radius, blade side,
normal direction, and pressure sign.

### Stage 4: Kutta-Joukowski Edge-Circulation Check

Compare the current Kutta-Joukowski panel-edge summation against a unique-edge jump
circulation calculation:

```text
γ_edge = μ_left - μ_right
```

This should distinguish duplicated panel-edge accounting from the physical
circulation jump carried by each shared edge.

### Stage 5: Independent Kernel-Offset Sweep

Sweep `kerneloffset_targets` independently from `kerneloffset_panel`.

The goal is to separate near-surface stabilization from thrust suppression. If target
offset dominates the loss, the pressure and force diagnostics should move with
`kerneloffset_targets` even when the panel/source regularization is held fixed.
