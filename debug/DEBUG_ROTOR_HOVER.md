# Rotor Hover Debug Log

Purpose: track the investigation into why `examples/rotor_hover.jl` underpredicts thrust by roughly a factor of 2-5. Expected converged thrust coefficient is approximately `0.10-0.14`. The current script only runs three stored timesteps, so these runs are for sensitivity triage, not convergence validation.

## Pre-existing Avenues To Investigate

- [ ] **CPoffset = R·1e-6 is sub-micron.** Interior control point sits at machine-epsilon distance from the panel. `solver.jl:229` flips sign for DBC, `abstractbody.jl:707–710` applies `off · characteristic_length · n` with no minimum-distance guard. Try CPoffset = 1e-3 to 1e-2 and see if Ct moves.
- [ ] **Wake-induced velocity is not in σ.** Implementation uses standard Morino: σ comes from `body.velocity` (freestream + kinematic ω×r only — `solver.jl:1100–1108`, `frames.jl:337–338`); wake influence is added to the doublet RHS as Φ_wake. This diverges from the user's stated formulation ("compute wake and freestream-induced velocity on the body, choose the source panel strengths accordingly"). Decide which formulation is intended and verify the code matches.
- [ ] **kerneloffset = R·1e-3 may smear peak blade suction.** Regularization length ~0.12 mm could be flattening Cp on thin blade sections. Sweep down to R·1e-5 or R·1e-6 and check sensitivity.
- [ ] **Wake under-resolution at 10 revs.** `nwakerows=1`, `p_per_step=2`, aggressive `MergeParticles(every=1, r=R·0.02)`: every step, particles within 2% R get merged. In hover the wake stays close to the disk for many revs; over-merging dissipates the tip-vortex circulation that drives induced inflow. Try `every=10`, smaller `r`, or disable merging and re-check.
- [ ] **Matrixful solver with rotating body — confirm RHS rebuilt each step.** `Backslash` factors the doublet influence matrix once at construction (body-fixed, OK). But Φ_σ on the RHS depends on σ which depends on ω×r at the rotated panel locations. Verify σ and Φ_σ are recomputed every step inside `simulate!`, not cached from t=0.
- [ ] **`set_Das_eta_kinematic=NaN` after init at 0.2.** `initialize_Das!` ran once at t=t₀ with η=0.2; during the run the kinematic refresh is disabled. If Das carries any rotation-dependent quantity that needs refreshing, it's frozen at t₀ orientation. Confirm what Das holds and whether it needs per-step updates for a rotor.
- [ ] **`p_correct_kuttacondition=false`** — resolved per user: linear Kutta is enforced in the matrix; this flag controls a separate iterative TE-loading correction. No action.

## Running Log

### Code-path investigation started: Das, dphidt, force timing

Plan for this pass:
- Inspect `initialize_Das!` and how `Das` is accumulated or refreshed during `simulate!`.
- Inspect the order of `solve!`, `calcfield_U!`, `dphidt`, pressure/force calculation, monitor calls, wake propagation, and kinematic propagation.
- Check whether the first-entry drop can be tied to `dphidt` history, frame update timing, or stale `Das`.

Findings will be appended as they are confirmed.

Confirmed code path:
- `simulate!` stores `-potential_old`, resets body/wake fields, applies freestream and kinematic velocity, snaps the wake TE row, applies wake-induced velocity, solves the body, applies body-on-target influence, computes `dphidt`, computes pressure/force, calls monitors, then propagates wake and body kinematics.
- `reset!(body)` clears `velocity`, `potential`, `P`, and `F`, so pressure/force accumulation across steps is not the source of the thrust drop.

No-`dphidt` diagnostic:

```text
Default:
[0.10403718845370113, 0.05608864677491295, 0.05371702245833515]

Override FLOWPanel._get_dphidt(::AbstractLiftingBody) = nothing:
[0.021337371170963692, 0.03186616843833708, 0.037584986767728236]
```

Interpretation: the unsteady Bernoulli term is adding thrust, especially at the initial step. The sustained underprediction is therefore not caused by an overly large adverse `dphidt`; if anything, the steady pressure part is even lower.

Potential issue found in Dirichlet wake coupling:
- `examples/rotor_hover.jl` uses `kernel = Union{ConstantSource,VortexRing}` and `DBC=true`.
- In the single-body Dirichlet `solve!`, `set_strengths!(body)` sets source strength from the current external velocity, then self source potential is computed, then `_solve!` solves the vortex-ring/doublet column from `solver.rhs .= -self.potential`.
- The main `simulate!` wake-on-target call currently requests `scalar_potential=false, gradient=true`.
- Therefore wake-induced velocity reaches the body before the solve, but wake scalar potential does not appear to be added to `body.potential` for the Dirichlet RHS.

Next diagnostic: temporarily request wake scalar potential on Dirichlet body targets during the wake-on-target influence call, then restore the source.

Wake scalar-potential diagnostic:

```text
Temporary change in src/FLOWPanel_simulate.jl:
scalar_potential = Tuple(sys isa AbstractBody && has_dirichlet_bc(sys) for sys in targets)

Result:
[0.10403718845370113, 0.057085052695284084, 0.054210946421464304]

Warnings:
Scalar potential was requested for a source system inducing a vector potential.
```

Interpretation: this only slightly raises steps 1-2 and produces vector-potential warnings, so the simple "include wake scalar potential" variant is not a clean explanation or fix. The source line was restored after the diagnostic.

Pure VortexRing/Neumann diagnostic on r2 mesh:

```text
kernel = pnl.VortexRing
DBC = false
msh_file = phantom_3_rebuild_r2.msh
watertight = true

Result:
[-7.830256169990838e15, -2.866818796402604e17, 1.178924124421047e18]
```

Interpretation: reject this run. The solver warned that a watertight `RigidWakeBody` with Neumann formulation is rank-deficient, as expected. Rerun with the non-watertight r4 mesh.

Pure VortexRing/Neumann diagnostic on r4 non-watertight mesh:

```text
kernel = pnl.VortexRing
DBC = false
msh_file = phantom_3_rebuild_r4.msh
te_indices_1 = [7, 952, 4] .+ 1
te_indices_2 = [3, 478, 0] .+ 1
watertight = false

Rotor: 5836 nodes, 11620 panels, 604 shedding edges
Result:
[0.10523930581029761, 0.054182084591703046, 0.05114289198202814]
```

Warnings:
- `fmm! called but either sources or targets are empty` at step 0.
- Scalar potential requested for a vector-potential source system on each step.

Interpretation: switching to the non-watertight mesh makes pure VortexRing/Neumann finite, but it does not recover the low sustained thrust. The first entry remains near the expected band; later entries remain near `0.05`.

Scalar-potential warning check:
- Cause: the post-solve system-on-target influence call in `simulate!` requested `scalar_potential=true` unconditionally for all targets.
- Pure `VortexRing` bodies define `FastMultipole.has_vector_potential(::AbstractBody{VortexRing,1}) = true`, so FastMultipole warns when scalar potential is requested from that source.
- The wake-on-target call was not the source of this warning; it already requests `scalar_potential=false`.

Guarded diagnostic:

```julia
scalar_potential = Tuple(sys isa AbstractBody && has_dirichlet_bc(sys) for sys in targets)
```

Pure VortexRing/Neumann r4 result with guarded scalar potential:

```text
[0.01902690748489106, 0.028369762760932164, 0.03348726539227742]
```

Default source+vortex Dirichlet r2 result with guarded scalar potential:

```text
[0.10403718845370113, 0.05608864677491295, 0.05371702245833515]
```

Interpretation: the guard removes the vector-potential scalar warning and is behavior-neutral for the default Dirichlet rotor case. For pure VortexRing/Neumann, it removes scalar-potential history and therefore removes most of the unsteady `dphidt` thrust contribution; this suggests the earlier pure-Neumann r4 result was partly relying on an invalid scalar-potential request.

Correction after FastMultipole interface update:
- User clarified that vortex rings are equivalent to doublets here, so scalar potential is valid; the guarded scalar-potential change was reverted.
- Rerunning pure `VortexRing`/Neumann r4 first exposed a local interface typo in `FastMultipole.body_to_multipole!(::AbstractBody{VortexRing,...})`: the method argument is `source_buffer`, but three `get_vertex` calls still referenced `buffer`.
- Patched those local references to `source_buffer` and reran.

Pure VortexRing/Neumann r4 after FastMultipole warning fix and local `source_buffer` typo fix:

```text
[0.10533669360782248, 0.054343245820128346, 0.051329910359658684]
```

Interpretation: the warning is now gone and the finite r4 Neumann result remains essentially unchanged from the earlier warning-producing run. It still does not recover sustained thrust.

### Follow-up sweep started: kerneloffset, Das initialization, short history

Plan for this pass:
- Sweep `kerneloffset` below baseline (`R*1e-5`, `R*1e-6`) and slightly above baseline (`R*3e-3`).
- Sweep `init_Das_eta_kinematic` at intermediate values (`0.5`, plus repeat/compare `1.0` if needed).
- Run short-history cases with 6 and 10 stored timesteps to see whether thrust stays near `0.05-0.06` or recovers.

Results will be appended below as the runs complete.

Kerneloffset follow-up results:

```text
kerneloffset = R * 3e-3:  [0.10397877907420751, 0.055978363366258575, 0.0535971747335778]
kerneloffset = R * 1e-5:  [0.10402965779761197, 0.056097448198295166, 0.053716774397490834]
kerneloffset = R * 1e-6:  [0.10402965779768465, 0.05609744838243152, 0.053716774479392396]
```

Interpretation: the baseline `R*1e-3` is on a broad flat part of the response. Lowering `kerneloffset` by 100-1000x does not recover thrust, and `R*3e-3` is nearly unchanged. The earlier `R*1e-2` case remains a useful upper bound where excessive smoothing starts to noticeably reduce thrust.

Das initialization follow-up results:

```text
init_Das_eta_kinematic = 0.5:  [0.13589416523718345, 0.06375710937169904, 0.058963778891512046]
init_Das_eta_kinematic = 1.0:  [0.16690276842179008, 0.06438820077846615, 0.06028794129254286]
```

Interpretation: increasing kinematic Das initialization strongly affects the first entry and modestly raises later entries, but the sustained values remain too low.

Short-history follow-up results:

```text
t_range = [1:6], defaults:
[0.10403718845370113, 0.05608864677491295, 0.05371702245833515, 0.05195069345956466, 0.05079937817563612, 0.05022658414714601]

t_range = [1:10], defaults:
[0.10403718845370113, 0.05608864677491295, 0.05371702245833515, 0.05195069345956466, 0.05079937817563612, 0.05022658414714601, 0.050089516882473324, 0.05028745225915549, 0.05074260562857607, 0.05133163974624757]
```

Interpretation: the low thrust is not just a 3-step artifact. It drops after the initial entry, bottoms near `0.050`, and only slowly recovers by 10 stored entries.

Kutta-Joukowski monitor cross-check:

```text
Pressure ForceMonitor y:
[0.10403718845370113, 0.05608864677491295, 0.05371702245833515]

KuttaJoukowskiForce y:
[-0.05665952079110518, -0.06039807714618753, -0.0601846445077248]
```

Interpretation: the Kutta-Joukowski monitor gives comparable sustained magnitude, with opposite y sign convention. The low sustained magnitude is therefore not obviously isolated to pressure integration.

### Verification after sensitivity-parameter grouping

After grouping additional knobs under the `Sensitivity parameters` heading in `examples/rotor_hover.jl`, reran the default script:

```text
Thrust Coefficient: [0.10403718845370113, 0.05608864677491295, 0.05371702245833515]
```

Interpretation: grouping the knobs was behavior-neutral for the current defaults.

### Initial sweep: r2 mesh, default formulation

Source setup:
- Mesh: `phantom_3_rebuild_r2.msh`
- Geometry-specific values left unchanged, including rotor radius `R`
- Thrust output recorded from `monitors[1].force[2,:]`
- Baseline warning observed on first step: `fmm! called but either sources or targets are empty; foregoing calculation`

Baseline:

```text
Thrust Coefficient: [0.10403718845370113, 0.05608864677491295, 0.05371702245833515]
```

Interpretation: first stored entry is inside the expected converged band, but entries 2-3 settle near `0.054-0.056`, which matches the underprediction symptom.

### Timestep sensitivity

Changed `nt` only, which changes `dt = 60 / RPM / nt`.

```text
nt = 18  (dt x2):  [0.07814513399212637, 0.05633241396255565, 0.0552442397956673]
nt = 72  (dt /2):  [0.16116816429707354, 0.05687558787723563, 0.0528207320514813]
```

Interpretation: timestep strongly affects the first stored force entry, but steps 1-2 remain close to baseline. This points toward startup/history/dphidt sensitivity, but does not by itself explain the sustained low value.

### Control-point and kernel regularization sensitivity

```text
CPoffset = R * 1e-4:       [0.10403830245625367, 0.05609018322247158, 0.05371876409813345]
CPoffset = R * 1e-8:       [0.10403717731374273, 0.056088631410526606, 0.05371700504205877]
kerneloffset = R * 1e-2:   [0.10008770597769986, 0.05039004637367297, 0.047270891946884094]
kerneloffset = R * 1e-4:   [0.1040296640196869, 0.056096727677079024, 0.05371646214005884]
kernelcutoff = R * 1e-10:  [0.10403718845370113, 0.05608864677491295, 0.05371702245833515]
kernelcutoff = R * 1e-16:  [0.10403718845370113, 0.05608864677491295, 0.05371702245833515]
```

Interpretation: `CPoffset` and `kernelcutoff` are effectively inactive for this symptom. Increasing `kerneloffset` to `R*1e-2` makes thrust lower, so excessive regularization can worsen the underprediction; decreasing it from the baseline does not recover thrust.

### Wake and force-processing sensitivity

```text
set_Das_eta_kinematic = 1.0:   [0.16690276842179008, 0.06438820077846615, 0.06028794129254286]
p_correct_kuttacondition = true: [0.09990178683259075, 0.05566839141091632, 0.054394514961894194]
p_per_step = 1:                [0.10403718845370113, 0.05608864677491295, 0.05409220516703139]
p_per_step = 4:                [0.10403718845370113, 0.05608864677491295, 0.053607577070894134]
overlap = 1.0:                 [0.10403718845370113, 0.05608864677491295, 0.05360910297940864]
overlap = 4.0:                 [0.10403718845370113, 0.05608864677491295, 0.054094719949877884]
```

Failure case:

```text
set_Das_eta_kinematic = 0.0
Result: failed after reaching step 1 with FastMultipole BoundsError while building the tree.
```

Interpretation: `set_Das_eta_kinematic` has the largest sustained effect found so far, raising later entries to about `0.060-0.064`, but still does not reach `0.10-0.14`. Kutta pressure correction and early wake particle insertion parameters are not major levers in the first three steps.

## Current Conclusions

- Highest-priority path: investigate startup/history handling around `initialize_Das!`, `dphidt`, and the force monitor, because the first entry can be in range while later entries collapse.
- Lower-priority knobs for this symptom: `CPoffset`, `kernelcutoff`, mild decreases in `kerneloffset`, `p_per_step`, and `overlap`.
- Excessive `kerneloffset` regularization is harmful but does not explain the baseline underprediction.
- The `set_Das_eta_kinematic = 0.0` failure is a separate robustness issue worth tracking if zero initialization is expected to work.

## Open Search Avenues

- Compare pressure-integral `ForceMonitor` with the commented `KuttaJoukowskiForce` monitor.
- Inspect `src/FLOWPanel_simulate.jl` around `dphidt` storage/computation and the order of force calculation, wake propagation, and kinematic propagation.
- Run a longer-but-still-short startup history check, for example 6-10 timesteps, to see whether `Ct` remains near `0.05-0.06` or climbs.
- Isolate mesh refinement after the hyperparameter pass:
  - `phantom_3_rebuild_r3.msh`: refined, watertight
  - `phantom_3_rebuild_r4.msh`: refined, not watertight for Neumann boundary conditions
- Keep mesh runs in a separate section because they change TE indices and should not be mixed with pure simulation hyperparameter sensitivity.

### Latent bug encountered: `source_buffer` undef in `body_to_multipole!`

While preparing a `μ`-evolution diagnostic, both the diagnostic script and the unmodified `examples/rotor_hover.jl` crashed with:

```
UndefVarError: `source_buffer` not defined in `FLOWPanel`
  body_to_multipole!  (FLOWPanel_liftingbody.jl:666)
```

Cause: the `RigidWakeBody{Union{ConstantSource,ConstantDoublet,VortexRing},2,...}` overload of `body_to_multipole!` was renamed `source_buffer → buffer` in its parameter list (line 658), but lines 666–668 still referenced `source_buffer`. From an uncommitted in-progress edit. Patched locally by renaming the three references to `buffer` (the sibling overload at 749 still uses `source_buffer` and remains correct). No effect on physics — it was preventing the run.

### μ-evolution diagnostic (10 steps, after fix)

Added a custom `mu_monitor` that logs `(maximum(abs, μ), mean(abs, μ))` per step alongside the `ForceMonitor`. Default settings, otherwise unmodified. Diagnostic lives in `examples/rotor_hover_mudiag.jl`; results write to `rotor_hover_mudiag/`.

```text
Ct (i_step=0..9): 0.10404, 0.05609, 0.05372, 0.05195, 0.05080, 0.05023, 0.05009, 0.05029, 0.05074, 0.05133
|μ|_max:          0.2054,  0.2815,  0.3226,  0.3532,  0.3777,  0.3952,  0.4081,  0.4181,  0.4264,  0.4339
|μ|_mean:         0.0819,  0.1049,  0.1203,  0.1309,  0.1383,  0.1438,  0.1481,  0.1517,  0.1551,  0.1585
```

Interpretation: between steps 1 and 9 the bound circulation `|μ|_mean` grows by ~50% while `Ct` is essentially flat (0.050–0.056). Lift per unit `μ` drops by ~⅓ over the same window. Combined with the prior cross-checks:

- Pressure ≈ Kutta–Joukowski (different formulas, similar magnitude) → force assembly is consistent with circulation; not a pressure-formula bug.
- DBC ≈ Neumann → not a BC-formulation bug.
- `|μ|` rising as expected → bound circulation is responding to wake; not a "Γ collapses" bug.

Leading hypothesis: **wake-induced inflow at the blade is too strong**, so the lift vector tilts further aft as the wake develops. The in-plane (torque) component grows, the out-of-plane (thrust) component stays flat. Consistent with the `Das` sweep (bigger `Das` → wake further from TE → less induced inflow → higher `Ct`).

### Refreshed open avenues

- [x] **Wake-induced velocity at body is over-predicted.** Confirmed. Probe at body CPs over 10 steps:

  ```text
  step  n·V_mean  |V|_mean  |V|_max   ωR_tip
    1    2.15      2.27      18.7    74.8
    5    5.80      6.16      49.0    74.8
    9    7.24      7.72      58.6    74.8
  ```

  Three observations:
  - `mean n·V_wake ≈ 7.2 m/s` at step 9 is already higher than the momentum-theory inflow `v_h ≈ 4.3 m/s` for the produced thrust (T≈2 N), and still climbing.
  - `max |V_wake| ≈ 58.6 m/s` is **78% of tip speed** — somewhere on the blade a wake source is sitting essentially on top of a body CP.
  - `max / mean ≈ 7.6` — the over-induction is concentrated in a small set of panels, not uniformly large.

  This is the bug. Bound circulation Γ is responding correctly to a fluid velocity that has been corrupted by an over-strong wake influence, so lift per unit Γ tilts aft and Ct collapses.

- [x] **Localize within the wake: panels (first row) vs particles.** The over-induction is **entirely in the panel wake** (one freshly-shed row), not in the particles:

  ```text
  step | panel n·V_mean   |V|_max  || particle n·V_mean   |V|_max
    1  |        2.15      18.7    ||         0.00          0.00
    5  |        7.13      52.6    ||         1.36          4.0
    9  |        8.21      62.0    ||         1.18          3.5
  ```

  Particles are well-behaved. The single-row panel wake (sat at `TE + Das` with `Das ≈ 4 mm` for `set_Das_eta_kinematic=0.2`) sees `|V|_max` grow **3.3×** over 9 steps while bound circulation only grows **1.5×**. So strength growth alone doesn't explain it — there is also a geometric drift bringing the wake panel closer to one (or a few) body CPs.

- [x] **Distance from worst body CP to nearest wake-panel vertex.** Logged each step:

  ```text
  step 1: idx=2150  |V|=18.7  d=1.71 mm
  step 5: idx=2150  |V|=52.6  d=1.71 mm
  step 9: idx=2142  |V|=62.0  d=1.69 mm
  ```

  The distance is **constant** across the run — same body CP, same wake vertex. Geometry is not drifting. The 3.3× growth in induced velocity is purely strength-driven (TE doublet jump grows with bound circulation). The body CP sits ~1.7 mm from a wake vertex throughout.

  Geometrically: blade-thickness/2 at this section (~35% radius) is ~1.7 mm. Das = −ω×r·dt·η is purely tangential (zero y-component since ω_axis is y), so the freshly-shed wake panel sits at the same y as the TE node — i.e., on the blade's mean line. Lower-surface CPs are blade-half-thickness *below* that mean line. So the wake panel is *structurally* ~half-blade-thickness from the lower-surface CPs near the TE.

- [x] **Wake `core_size = 1 mm` < blade-half-thickness.** With wake regularization radius 1 mm and structural distance 1.7 mm, the wake-panel kernel evaluated at these body CPs is essentially unregularized, and the doublet is near-singular. This corrupts `body.velocity` at those CPs → σ → the doublet solve → bound circulation. **This is the regularization mismatch driving the under-prediction.**

- [ ] **Test fix: increase wake `core_size` to ≥ blade thickness.** If `core_size = 5e-3` (5 mm) or so, the body CPs at 1.7 mm sit well inside the regularization core, the kernel is smooth, and induced velocity at those CPs should drop substantially. If Ct rises toward 0.10, this is the fix. If Ct doesn't move, the regularization story is incomplete and we need a different angle.
- [ ] **Lift-vector tilt vs. Ct decomposition.** Once wake-induced inflow magnitude is in hand, decide whether the magnitude is physical (then the steady answer really is near `0.054` and the `0.10–0.14` expectation needs revisiting) or non-physical (then trace the wake source).

## Sensitivity Parameters Now Exposed In Script

The following defaults are grouped under the `Sensitivity parameters` heading in `examples/rotor_hover.jl` for future one-at-a-time sweeps:

```julia
CPoffset
kerneloffset
kernelcutoff
p_per_step
overlap
merge_r_factor
merge_r_hash_factor
init_Das_eta_kinematic
p_correct_kuttacondition_flag
```

### Wake `core_size` sweep execution started

First implementation attempt exposed a constructor forwarding bug before any physics result:

```text
WAKE_CORE_SIZE=1e-3 julia --project=. examples/rotor_hover_mudiag.jl
ERROR: no method matching PanelWake(...; nwakerows::Int64, core_size::Float64)
```

Cause: `PanelParticleWake(body; kwargs...)` forwards extra keywords to `PanelWake(body; ...)`, but the `PanelWake(body; ...)` convenience method only accepted `nwakerows` and did not forward `core_size` to the lower-level `PanelWake(shedding, kernel; core_size, nwakerows)` constructor. Patch next: forward keyword arguments through `PanelWake(body; nwakerows, kwargs...)`.

After forwarding `kwargs...` through `PanelWake(body; ...)`, the default 1 mm case is behavior-neutral:

```text
wake_core_size = 1e-3:
Ct = [0.10403718845370113, 0.05608864677491295, 0.05371702245833515, 0.05195069345956466, 0.05079937817563612, 0.05022658414714601, 0.050089516882473324, 0.05028745225915549, 0.05074260562857607, 0.05133163974624757]

panel wake at step 9: n·V_mean = 8.213 m/s, |V|_mean = 8.673 m/s, |V|_max = 62.037 m/s
particle wake at step 9: n·V_mean = 1.177 m/s, |V|_mean = 1.362 m/s, |V|_max = 3.497 m/s
worst CP distance at step 9: 1.6886 mm
```

Interpretation: the core-size plumbing did not change the default result; this is the baseline for the sweep.

2 mm wake core result:

```text
wake_core_size = 2e-3:
Ct = [0.10403718845370113, 0.05364284714235649, 0.05045288923872743, 0.04828441496623624, 0.04695173481087879, 0.046277527696970924, 0.0461202820246109, 0.046256571470534864, 0.04663011839868124, 0.04717810798622627]

panel wake at step 9: n·V_mean = 6.847 m/s, |V|_mean = 7.203 m/s, |V|_max = 33.371 m/s
particle wake at step 9: n·V_mean = 0.946 m/s, |V|_mean = 1.094 m/s, |V|_max = 2.501 m/s
worst CP distance at step 9: 1.9974 mm
```

Interpretation: increasing wake `core_size` from 1 mm to 2 mm cuts the worst panel-wake velocity roughly in half, but `Ct` decreases rather than recovers. This weakens the simple "larger core fixes thrust" hypothesis; continue to 5 mm because that was the proposed value.

5 mm wake core result:

```text
wake_core_size = 5e-3:
Ct = [0.10403718845370113, 0.039817194833772704, 0.03692926893193312, 0.03465543284331525, 0.033687447922625555, 0.0332523001443505, 0.03322945792601928, 0.033437479998539583, 0.03379008566373044, 0.034238214421292976]

panel wake at step 9: n·V_mean = 3.696 m/s, |V|_mean = 3.868 m/s, |V|_max = 9.432 m/s
particle wake at step 9: n·V_mean = 0.533 m/s, |V|_mean = 0.642 m/s, |V|_max = 1.836 m/s
worst CP distance at step 9: 4.7219 mm
```

Interpretation: 5 mm strongly suppresses the local panel-wake velocity but drives sustained `Ct` even lower (~0.034). The proposed larger `core_size` is not the fix. Next direction, per user suggestion: sweep below the 1 mm default.

0.5 mm wake core result:

```text
wake_core_size = 5e-4:
Ct = [0.10403718845370113, 0.0563620407872271, 0.054120061822274794, 0.052409311443605895, 0.05131053917997011, 0.05077935697181704, 0.05066834374875518, 0.050890950883371006, 0.0513590646827884, 0.05196687507410509]

panel wake at step 9: n·V_mean = 8.463 m/s, |V|_mean = 8.947 m/s, |V|_max = 67.201 m/s
particle wake at step 9: n·V_mean = 1.221 m/s, |V|_mean = 1.417 m/s, |V|_max = 3.663 m/s
worst CP distance at step 9: 1.6886 mm
```

Interpretation: reducing wake `core_size` below default slightly raises sustained `Ct`, but only by about `+0.0006` at step 9 while increasing the local peak velocity. This is not enough to explain the factor-of-2 discrepancy. Continue with 0.2 mm to check monotonicity.

0.2 mm wake core result:

```text
wake_core_size = 2e-4:
Ct = [0.10403718845370113, 0.05638260643818949, 0.05414435071554057, 0.0524241822420329, 0.05134340148184559, 0.05080163447794444, 0.05069830230951361, 0.05092425556997373, 0.051412403505964085, 0.05202397337848889]

panel wake at step 9: n·V_mean = 8.483 m/s, |V|_mean = 8.970 m/s, |V|_max = 67.540 m/s
particle wake at step 9: n·V_mean = 1.223 m/s, |V|_mean = 1.418 m/s, |V|_max = 3.679 m/s
worst CP distance at step 9: 1.6886 mm
```

Interpretation: smaller wake cores saturate quickly. Going from 1.0 mm to 0.2 mm raises step-9 `Ct` only from `0.05133` to `0.05202`; the discrepancy is not controlled by wake-panel core size alone.
