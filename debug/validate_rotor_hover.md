# Validate `rotor_hover_pressure_comparison.jl`

Goal: validate `examples/rotor_hover_pressure_comparison.jl`, starting from the saved run
`data/rotor_hover_pressure_cacceleration_form=:lamb_vectoromparison_nt36_2`.

Working hypothesis: this run predicted plausible rotor forces with the
`KuttaJoukowskiForce` monitor, while the pressure-integrated monitors did not yet
match that level. Use this file to keep the run settings and validation checks in
one place.

## Key Parameters

| Parameter | Current example default | Checked saved run value | Source / notes |
| --- | ---: | ---: | --- |
| Run name | `rotor_hover_pressure_comparison_nt72_sfs` | `rotor_hover_pressure_cacceleration_form=:lamb_vectoromparison_nt36_2` | Current script vs. `*.replay.toml` / `*.frames.toml`. Saved name appears to include an accidental string interpolation artifact for `acceleration_form=:lamb_vector`. |
| Geometry file | `examples/data/dji9443_new_40_40.msh` | not stored in replay metadata | Inferred from current example only unless the original command/log is found. |
| Radius `R` | `0.119` m | `0.119` m inferred | Replay has `kerneloffset_panel=1.19e-11 = R*1e-10` and `kernelcutoff=1.19e-14 = R*1e-13`. |
| Density `rho` | `1.071778` kg/m^3 | not stored in replay metadata; replay helper default is `1.225` | Original run value needs confirmation from shell history or stdout log. Existing replay diagnostics used `rho=1.225` unless overridden. |
| RPM | `6000` | `6000` inferred | `omega=628.3185307179587 rad/s` in frames, equal to `2*pi*6000/60`. |
| `omega_axis` | `[-1, 0, 0]` for DJI9443 | `[-1, 0, 0]` | `*.frames.toml`. |
| Axial dimension | `1` for DJI9443 | `1` inferred | Current example and replay helper default. |
| `magVinf` | `0.0001` m/s | effectively `0` in replay diagnostics | `q_uinf=0` in existing replay CSV; saved metadata does not store `Vinf`. |
| `AOA` | `0.0` deg | not stored | Current example only. |
| `nt` | env `NT`, default `72` | `36` inferred | `dt=0.0002777777777777778 s`; at 6000 RPM one rev is `0.01 s`, so `nt=0.01/dt=36`. |
| `nrevs` | env `NREVS`, default `10` | `10` inferred | 360 saved frame steps / 36 steps per rev. |
| `n_steps` | `round(nt*nrevs) + spinup_steps` | `360` | Counted `[[step]]` entries in `*.frames.toml`. |
| Time span | `0:dt:...` | `t=0.0` through `0.09972222222222223` | `*.frames.toml` and body `.pvd`. |
| Spinup revs | env default `0.0` | no spinup evident | Constant `omega` from first saved frame. |
| `CPoffset` | `R*1e-8 = 1.19e-9` | `-1.19e-9` | `*.replay.toml` stores negative `cp_offset`; sign may reflect replay/body convention. |
| `kerneloffset_panel` | `R*1e-10 = 1.19e-11` | `1.19e-11` | `*.replay.toml`. |
| `kerneloffset_targets` | env default `0.003` | `0.003` | `*.replay.toml`. |
| `kernelcutoff` | `R*1e-13 = 1.19e-14` | `1.19e-14` | `*.replay.toml`. |
| Kernel | `Union{ConstantSource,VortexRing}` | `strength_names = ["sigma", "gamma"]` | `*.replay.toml`; `dbc=true` saved. |
| `DBC` | `true` for current kernel expression | `true` | `*.replay.toml`. |
| Shedding root cutoff | `shedding_r_over_R = 0.1` | not stored directly | Saved shedding has two shedding lines with 35 entries each. Need root midpoint recomputation if this parameter matters. |
| Wake type | `PanelParticleWake` | `PanelParticleWake` | `*.replay.toml`. |
| `nwakerows` | `1` | `1` | `*.replay.toml`. |
| `max_particles` | `500000` | `500000` | `*.replay.toml`. |
| Wake core size | env default `kerneloffset_targets` | `0.003` | `*.replay.toml`. |
| Particle kernel offset | `kerneloffset_targets` | `0.003` | `*.replay.toml`. |
| `shed_with_induced_velocity` | `false` | `false` | `*.replay.toml`. |
| `unsteady_filament` | `false` | `false` | `*.replay.toml`. |
| Viscous model | `FLOWVPM.CoreSpreading(wake_nu, wake_core_size, FLOWVPM.zeta_fmm; beta=wake_core_beta)` | not stored in saved replay metadata | Current example only. Do not assume for saved runs unless an original script/stdout captures it. |
| `wake_nu` | `1.85508e-5 / rho = 1.7308435142352243e-5` with current `rho=1.071778` | not stored in saved replay metadata | Current example default, overridable by `WAKE_NU`. |
| `wake_core_beta` | env `WAKE_CORE_BETA`, default `1.5` | not stored in saved replay metadata | Current example only. |
| Viscous kernel | `FLOWVPM.zeta_fmm` | not stored in saved replay metadata | Current example only. |
| SFS model | `FLOWVPM.SFS_Cd_twolevel_nobackscatter` | not stored in saved replay metadata | Current example only. The `nt72_sfs` run name suggests SFS was intended, but the saved files do not encode it. |
| Trailing wake particle method | `SigmaOverlap(R*0.05, 4.0)` | not stored in saved replay metadata | Current example only. With `R=0.119`, first argument is `0.00595`. |
| Unsteady wake particle method | `NoShed()` | not stored in saved replay metadata | Current example only. |
| Particle maintenance | `MergeParticles(every=merge_particles ? 1 : 0, r=0.02R, r_hash=0.02R, sigma_relative=false)` | not stored in saved replay metadata | Current example defaults: `merge_particles=true`, `r=r_hash=0.00238`. |
| Pressure acceleration form | current monitor includes `:lamb_vector` when `lamb_only=true` | run name indicates `:lamb_vector` | Name evidence only; replay can recompute variants independent of the original monitor set. |
| `lamb_only` | env default `true` | likely `true` | Current example default; original command/log needed for certainty. |
| `RUN_KJ` | env default `false` | original unknown; replay diagnostics ran `RUN_KJ=true` | Saved VTK does not encode monitors. Existing debug CSV includes KJ from replay. |
| Main FMM backend | order `8`, MAC `0.4`, leaf `20` | replay diagnostics used same defaults | Current example and replay helper defaults. |
| KJ FMM backend | current example order `3`, MAC `0.4`, leaf `1000` | replay helper default order `4`, MAC `0.4`, leaf `50` unless overridden | Important mismatch: KJ validation should record which script is used. |
| Normalization | `RotorNormalization(rho, 2R, 1)` | same in replay diagnostics | Current example and replay helper. |

## Checked Force Signals

Current replay output from:

```bash
julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

uses the default saved steps `357:359` and reports:

| Monitor / path | CT axial range | Notes |
| --- | ---: | --- |
| `kutta_joukowski` | `-0.06903` to `-0.07046` | Stable and close to the replay target `CT_EXPERIMENT=-0.072`. |
| `laplace_lamb` | about `-0.00300` | About 23x smaller than KJ in this replay check. |
| `laplace_md` | about `-5.0e-5` | Much smaller than KJ. |
| `bernoulli` | first selected step `-0.08204`, then near `-1e-4` | Indicates replay/Bernoulli state handling needs separate scrutiny before using it as a baseline. |

Note: older existing CSVs under `debug/logs/rotor_hover_stage*_20260526_final4.csv`
show `kutta_joukowski` around `-0.031`. Those rows were used in the first draft of
this note, but they are stale relative to the current replay script/worktree.

### nt72_sfs replay notes

Replaying `data/rotor_hover_pressure_comparison_nt72_sfs` with stage-1 logging
enabled gives a stable KJ coefficient near `-0.05` on the last three saved
samples:

| Sample | CT axial |
| --- | ---: |
| step 662 | `-0.050964309` |
| step 663 | `-0.050475825` |
| step 664 | `-0.049400557` |
| mean | `-0.05028023` |

On the same samples, the pressure-derived monitors do not match KJ:

| Monitor | CT axial range / mean |
| --- | ---: |
| `bernoulli` | `-0.13667324` to `0.00031150275` |
| `laplace_md` | about `-5.0069342e-5` |
| `laplace_lamb` | about `-0.002994` |

The saved nt72_sfs run also ends earlier than the nt36 run:

| Run | Saved steps | Final time |
| --- | ---: | ---: |
| `rotor_hover_pressure_comparison_nt72_sfs` | `665` | `0.09222222222222222` s |
| `rotor_hover_pressure_comparison_nt36` | `720` | `0.09986111111111111` s |

That makes the two runs a non-like-for-like late-time comparison, and it is a
plausible contributor to the pressure-force discrepancy.

## Saved Wake Metadata Coverage

The saved replay files for these rotor-hover runs all contain the same wake fields:
`type`, `nwakerows`, `max_particles`, `core_size`, `particle_kerneloffset`,
`shed_with_induced_velocity`, and `unsteady_filament`.

| Run | Saved wake fields checked | Missing from saved files |
| --- | --- | --- |
| `rotor_hover_pressure_comparison_nt36` | `PanelParticleWake`, `nwakerows=1`, `max_particles=500000`, `core_size=0.003`, `particle_kerneloffset=0.003`, `shed_with_induced_velocity=false`, `unsteady_filament=false` | viscous model, `wake_nu`, viscous kernel, `wake_core_beta`, SFS model, trailing method, unsteady method, merge settings |
| `rotor_hover_pressure_comparison_nt72_sfs` | `PanelParticleWake`, `nwakerows=1`, `max_particles=500000`, `core_size=0.003`, `particle_kerneloffset=0.003`, `shed_with_induced_velocity=false`, `unsteady_filament=false` | viscous model, `wake_nu`, viscous kernel, `wake_core_beta`, SFS model, trailing method, unsteady method, merge settings |
| `rotor_hover_pressure_cacceleration_form=:lamb_vectoromparison_nt36_2` | `PanelParticleWake`, `nwakerows=1`, `max_particles=500000`, `core_size=0.003`, `particle_kerneloffset=0.003`, `shed_with_induced_velocity=false`, `unsteady_filament=false` | viscous model, `wake_nu`, viscous kernel, `wake_core_beta`, SFS model, trailing method, unsteady method, merge settings |

Implication: for validation, SFS and viscous settings must be tracked from the
generating script or run command going forward. Current saved files are sufficient
to replay the wake state, but they are not sufficient to prove which SFS/viscous
model generated that state.

New runs should write the unified `*.metadata.toml` manifest, which now records
the wake constructor provenance above together with frame history and backend
settings. Legacy `*.replay.toml` + `*.frames.toml` saved runs remain replayable,
but they still lack the SFS/viscous provenance listed in the table.

## SFS / Panel-Particle Coupling Risk

Investigation date: 2026-06-01.

FLOWVPM's SFS model is not just a diagnostic field. In
`../FLOWVPM.jl/src/FLOWVPM_timeintegration.jl`, the particle update subtracts
`C * SFS * sigma^3 / zeta0` from the circulation/stretching update. For the
dynamic two-level model used by `FLOWVPM.SFS_Cd_twolevel_nobackscatter`,
`BeforeUJ` temporarily scales particle `sigma`, recomputes particle-field U/J and
SFS, stores test-filter stretching/SFS in `M`, and `AfterUJ` compares that state
against the current domain-filter U/J and SFS to calculate `C`.

The hybrid wake path in this repo can feed panel-induced velocity gradients into
that `AfterUJ` stage:

| Code path | Observation | Consequence |
| --- | --- | --- |
| `simulate!` wake-on-target pass | `influence!(targets, wake_sources; precalc=true, velocity_gradient=Tuple(requires_hessian(sys) for sys in targets))` calls `pre_evaluate_influence!` on the particle field, then accumulates wake panel, filament, and particle influence onto particles. | `precalc=true` runs `pfield.SFS(pfield, BeforeUJ())`, but the subsequent coupled influence pass writes U/J from all wake sources onto particles. |
| `simulate!` body-on-target pass | `influence!(targets, systems_tuple; velocity_gradient=Tuple(requires_hessian(sys) for sys in targets))` includes the particle field in `targets`. | Because `requires_hessian(::FLOWVPM.ParticleField) = true`, body/panel-on-particle Hessians are accumulated into particle `J`. |
| `PanelParticleWake.propagate!` | Calls `FLOWVPM._euler(w.pfield, dt; relax)`, not `FLOWVPM.euler(...; custom_UJ=...)`. | `_euler` immediately calls `pfield.SFS(pfield, AfterUJ())` and then updates particles using the already accumulated coupled U/J. There is no particle-only U/J recomputation at this point. |
| FLOWVPM dynamic SFS internals | `dynamicprocedure_pseudo3level_beforeUJ` calls `pfield.UJ(pfield; sfs=true, reset=true, reset_sfs=true)`, which is particle-only. `dynamicprocedure_pseudo3level_afterUJ` uses the current `J` and `SFS` fields. | The dynamic comparison can be asymmetric: test-filter terms are particle-only, while domain-filter terms may include panel-on-particle gradients. |

Assessment: yes, the present coupling can induce spurious SFS predictions when an
SFS model is enabled for `PanelParticleWake`. The most likely mechanism is not
panel velocity convecting particles; that is physically expected in a coupled
panel/VPM wake. The concern is panel-on-particle velocity gradient entering the
VPM stretching tensor `J` and therefore the dynamic SFS coefficient/correction,
even though the FLOWVPM SFS derivation and implementation assume a filtered
particle vorticity field. Near the blade or near the panel wake handoff, panel
gradients are singular/regularization-sensitive and can look like unresolved
strain to the dynamic model. This can alter particle circulation, not just
particle positions.

Current worktree note: `examples/rotor_hover_pressure_comparison.jl` currently
has both `CoreSpreading(...)` and `SFS=pnl.FLOWVPM.SFS_Cd_twolevel_nobackscatter`
commented out in the `PanelParticleWake` constructor, so the current default
fresh run does not exercise this SFS failure mode. The legacy `nt72_sfs` naming
still cannot prove SFS was active because legacy replay metadata does not store
the SFS constructor.

Recommended validation checks:

| Check | Purpose |
| --- | --- |
| Run matched cases with `SFS=FLOWVPM.SFS_default` and `SFS_Cd_twolevel_nobackscatter`, all else fixed. | Quantify whether SFS is changing thrust/torque or particle circulation beyond expected wake diffusion. |
| Log particle `C`, `SFS`, and `J` norms before body-on-particle influence, after body-on-particle influence, and after `pfield.SFS(..., AfterUJ())`. | Directly confirm whether panel gradients dominate dynamic SFS coefficients near the blade/wake handoff. |
| Try `CoreSpreading(...)` plus `SFS_Cd_twolevel_nobackscatter` while forcing particle-target `velocity_gradient=false` for panel-on-particle influence, and likely for all non-particle-on-particle interactions. Keep panel velocity on particles. | Test whether the VPM viscous/SFS model can be used for particle wake evolution without feeding bound-panel or panel-wake strain into the SFS/stretching tensor. |
| Repeat with particle `velocity_gradient=false` for body-on-particle targets while keeping panel velocity on particles. | Minimal version of the gradient-isolation test; separates convective coupling from panel-gradient/stretching/SFS coupling. |
| If SFS is needed with panel coupling, consider a custom propagation path that recomputes particle-only U/J for SFS/stretching or strips panel contributions from `J` before `AfterUJ`. | Avoid treating bound-panel strain as unresolved particle-scale turbulence. |

## Validation Checklist

| Check | Status | Command / artifact | Notes |
| --- | --- | --- | --- |
| Confirm saved metadata loads | done | `data/.../*.metadata.toml` or legacy `data/.../*.replay.toml` + `data/.../*.frames.toml` | New saves should use the unified manifest; old saves remain readable. |
| Confirm saved step count and `nt` | done | count `[[step]]` entries; inspect `dt` | 360 steps, `nt=36`, 10 revs inferred. |
| Recover exact original environment | open | shell history/stdout if available | Need `rho`, `RUN_KJ`, `LAMB_ONLY`, and any FMM/KJ overrides from the original run. |
| Compare current example defaults to saved run | partial | table above | Current default run name and `nt=72` differ from saved `nt=36` run. |
| Replay KJ force on saved run | done | `julia --project=. examples/rotor_hover_pressure_comparison_replay.jl` | Current KJ CT around `-0.070`; stale debug CSVs show older `-0.031` values. |
| Replay pressure-force monitors on saved run | done | `julia --project=. examples/rotor_hover_pressure_comparison_replay.jl` | Pressure-derived CT does not match KJ in current replay checks. |
| Validate KJ backend sensitivity | open | rerun replay varying `KJ_FMM_*` or direct backend if feasible | Current example and replay helper use different KJ backend defaults. |
| Validate `rho` sensitivity | open | rerun replay with `RHO=1.071778` and `RHO=1.225` | Rotor CT normalization and force both depend on density; CT may cancel partly but monitor internals still should be checked. |
| Preserve wake model provenance | done | unified metadata writer | New saves store SFS/viscous model identities and parameters in `*.metadata.toml`; legacy saves do not. |
| Validate force sign convention | open | compare axial dimension and frame/global signs | KJ is negative in axial dimension 1 for this saved run. |
| Validate experiment target | open | document source for `CT_EXPERIMENT=-0.072` | Needed before calling the KJ force "accurate" quantitatively. |

## Reproduction Commands

Narrow replay of the saved run with KJ and stage logging:

```bash
env RUN_STAGE1_LOG=true NSTEPS_REPLAY=4 RUN_KJ=true STAGE1_LOG_PATH=debug/logs/rotor_hover_stage1_validate.csv julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

Replay only PressureLaplace variants plus KJ:

```bash
env RUN_STAGE2_LAPLACE=true RUN_STAGE1_LOG=true RUN_BERNOULLI=false RUN_LAPLACE_MD=false RUN_LAPLACE_LAMB=false RUN_KJ=true NSTEPS_REPLAY=4 STAGE1_LOG_PATH=debug/logs/rotor_hover_stage2_validate.csv julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
```

Fresh run of the current example configured to match the saved run cadence:

```bash
env NT=36 NREVS=10 RUN_KJ=true LAMB_ONLY=true SAVE_VTK=true julia --project=. examples/rotor_hover_pressure_comparison.jl
```

## Notes

- The saved run metadata is enough to reconstruct geometry state, wake state, frame
  motion, and body strengths, but not every original environment variable.
- The current example default `run_name` and `nt` no longer match the saved
  `nt36` run, so validation should explicitly pin `NT=36` until the comparison is
  closed.
- The string `rotor_hover_pressure_cacceleration_form=:lamb_vectoromparison_nt36_2`
  should be preserved in commands for replay, because it is the actual directory
  and metadata name.

## Live-run Stability Sweep

Goal: starting from the current confirmed-stable script defaults of
`examples/rotor_hover_pressure_comparison.jl`, lower `KERNELOFFSET` one step
at a time and identify the smallest value that still produces a stable
final-rev CT. Stability target at this stage is **a stable converged value
(small step-to-step drift over the last rev)**, not a specific magnitude;
matching `CT_EXPERIMENT=-0.072` is *not* required yet. Once a candidate value
is found, verify it is in the converged regime by checking that the adjacent
larger ladder value gives essentially the same converged CT (target: a few
percent on the final-rev mean of the lamb-vector CT).

Working baseline (script defaults, confirmed stable by the user):
`KERNELOFFSET=1e-2`, `ramp_nrev=5`, `magVinf_start=0`, no SFS, no viscous
`CoreSpreading`, `NT=36`, `NREVS=10`, `RUN_KJ=false`, `LAMB_ONLY=false`,
`SAVE_VTK=true`, `MERGE_PARTICLES=true`, all other env vars at their script
defaults. Run under `JULIA_NUM_THREADS=3`. Output is saved to
`data/rotor_hover_pressure_comparison/`.

Launch command:

```bash
env JULIA_NUM_THREADS=3 KERNELOFFSET=<value> \
    julia --project=. examples/rotor_hover_pressure_comparison.jl
```

Planned runs (advance only if the previous row is stable):

| ID | `KERNELOFFSET` | Status | Final-rev CT_axial (lamb / matderiv / bernoulli) | Notes |
| --- | ---: | --- | --- | --- |
| B0 | `1e-2` | done — stable | TBD (record from existing baseline if logged) | Current script default, user-confirmed stable |
| R1 | `5e-3` | done — stable | -0.0188 / -0.189 / -0.300 | Final-rev (steps 325–360) mean; periodic oscillation, no drift. Log: `debug/logs/rotor_hover_R1_kerneloffset_5e-3.log` |
| R2 | `2e-3` | pending | Run only if R1 stable | |
| R3 | `1e-3` | pending | Run only if R2 stable | User's stretch target |

Decision rule per row: if the run completes all `n_steps=360` steps and the
final-rev CT has small step-to-step drift, mark stable and advance to the next
smaller ladder value; if it kills early, NaNs, or visibly diverges, stop the
sweep and keep the immediately-larger value as the answer.
