# Plan (2026-07-02): CT knob sweep — body strain, surface smoothing radius, relaxation convergence

Stand-alone: a clear-context agent should be able to execute this on a workstation without
prior conversation. All paths relative to the repo root
`/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (or its workstation clone). **Thread
budget: 48 — every simulation subprocess must run with `-t 48` AND BLAS threads 48; the
driver sets `OPENBLAS_NUM_THREADS`/`OMP_NUM_THREADS`/`JULIA_NUM_THREADS` from its `THREADS`
env (default 48).**

## 0. OPERATING RULES — token budget (read first, non-negotiable)

The user has a 5-hour usage limit that this agent MUST NOT exhaust. Simulations cost
wall-clock, not tokens — the only token spend is the agent's own turns. Rules:
Avoid spending tokens during peak token hours of 7am-1pm MST (background runs are fine), and avoid having jobs running at all during the 24-hour period starting at 12am MST on Sundays.

1. Launch scenario batches **in the background** and wait with **long intervals (20–30
   min)** between checks; never poll frequently.
2. Each check is a **minimal probe**: `tail -5` of the scenario log or a glance at
   `data/ct_knob_sweep_summary.csv`. **NEVER read a full scenario log** — each contains a
   ~500-line per-step CT table plus per-step solver output (tens of thousands of lines).
   Parse the CSVs, grep the logs.
3. Substantive analysis/writing only at scenario-batch boundaries — a few real turns per
   multi-hour batch. Between boundaries, do nothing.
4. Keep a rough running self-estimate of spend; target **well under half the budget per
   5-hour window**. When uncertain, idle longer — a late result costs nothing, a blown
   budget kills the session. If the user supplies a real meter reading, recalibrate.
5. When proposing anything token-heavy to the user, include an estimated fraction of the
   5-hour limit (project CLAUDE.md rule).

## 1. Context and goal

Rotor-hover CT under-predicts experiment (~0.072); the settled baseline run
`data/rotor_hover_pressure_comparison` (RHPC_MESH=40_40, NT=48 ⇒ dt=1/3600 s, RPM=6000,
corrected Pedrizzetti rlxf=0.3, 360 steps = 10 revs, freestream fully withdrawn by rev 9)
reads **CT_bernoulli = 0.06238 (final-rev mean; std 3.8e-4, peak-to-peak 1.2e-3)**;
CT_laplace_lamb = 0.06213. Items 003/004 route the shortfall to wake modeling. Three knobs
are tested as **warm-started perturbations**: each scenario resumes the baseline at step
350 (rev 9.72, settled hover) and continues to 14 revs (~154 steps) with ONE knob changed.

Knob rationale (established earlier on this branch, recorded in
`BRAINSTORM/008_particle_overlap_vorticity_evolution.md`, 2026-07-02 sections):

1. **Surface-induced strain on particles.** The run advects particles with the total
   velocity but stretches them with the particle-only gradient
   (`BODY_HESSIAN_TO_PARTICLES=false`, `PANEL_WAKE_HESSIAN_TO_PARTICLES=false`). An
   offline audit (`examples/particle_body_strain_audit.jl`) measured the omitted
   stretching at **3.2× the retained term (energy norm) in the youngest wake band**,
   negligible beyond ~0.11 axial. Restoring it is physical; instability is the suspected
   failure mode (the wake-row half was originally ablated against spurious 1/r² edge
   fields — expect that scenario to be the noisier one).
2. **Surface→particle smoothing radius** `KERNELOFFSET_TARGETS` (baseline 1e-3).
   Confirmed as the correct knob: `FLOWPanel_abstractbody.jl:38` ("panel influence on
   external targets"), activated in `_steady_aerodynamics!` right before the
   body→particles pass. `KERNELOFFSET_PANEL` is the panel-solve offset — do NOT sweep it.
   CRITICAL: the example defaults `WAKE_CORE_SIZE` to `KERNELOFFSET_TARGETS`, so the sweep
   pins `WAKE_CORE_SIZE=1e-3`; otherwise shed-particle cores would be co-swept. Prior
   experience: shrinking raised CT until instability.
3. **Split regularization for the body→particle gradient** (`BODY_GRADIENT_KERNELOFFSET`,
   new src knob, this session). With `BODY_HESSIAN_TO_PARTICLES=true`, the body→particle
   velocity GRADIENT is evaluated at this larger kernel offset while the velocity itself
   keeps `KERNELOFFSET_TARGETS`. Not strictly physical — it smooths the |∇U| "bumpiness"
   that piecewise-constant doublet panels imprint on near-blade particles, and is the
   designated stabilizer if plain `bodyhess` diverges. Implemented as a two-pass influence
   in `_steady_aerodynamics!` (velocity pass at the sharp offset; gradient-only pass at
   the large offset with the co-computed velocity discarded via snapshot/restore, since
   backends may not support hessian-only evaluation). NaN disables; only active when
   `BODY_HESSIAN_TO_PARTICLES=true`. Costs one extra body→particles FMM pass per step.
   Recorded in each run's metadata.toml under solver_options.
4. **Particle merging (`MERGE_PARTICLES`, baseline `true`).** The baseline merges
   particles every step (`MergeParticles` monitor, `r=0.02R`, `r_hash=0.02R`), which
   coarse-grains the wake and diffuses vorticity — a candidate CT suppressor. Turning it
   off (`MERGE_PARTICLES=false`, scenario `merge_off`) keeps every shed particle for the
   ~154-step continuation; if CT rises, merge-induced diffusion is part of the shortfall.
   Cost/risk: particle count grows unbounded over the continuation (per-step walltime and
   memory climb; watch for slowdown, not a hard crash), and un-merged near-duplicate shed
   pairs are structural, not overlap error (see
   `project_particle_shed_pairs_overlap`) — do not read a rising np as instability by
   itself. A finite continuation is short enough that unbounded growth is tolerable.
5. **Pedrizzetti relaxation, treated as artificial dissipation.** rlxf is a per-step
   blend; the physically meaningful parameter is the rate λ = rlxf/dt (baseline
   λ = 0.3·3600 = 1080 s⁻¹). Convergent, tuning-free procedure: run the halving ladder
   rlxf ∈ {0.3 (=control), 0.15, 0.075} at fixed dt, verify monotone CT(λ) over the
   stable range, Richardson-extrapolate CT(λ→0) (the driver does this automatically when
   ≥3 ladder points are stable). Independently, `RELAX_FILTER_DOWNSTREAM_R=0.5` confines
   relaxation to ≥0.5R downstream — if CT plateaus as the unrelaxed near-rotor band grows,
   that is a second, extrapolation-free route to a relaxation-free CT. rlxf→0 directly is
   known to destabilize (item 005); instability at small rlxf is an expected outcome that
   bounds the ladder, not a failure.

## 2. What is already implemented (this session — do not rewrite)

- `src/FLOWPanel_warmstart.jl`: `simulate_warmstart!` now forwards `optargs...` to
  `simulate!` (so `body_hessian_to_particles` etc. take effect in continuations).
- `src/FLOWPanel_warmstart.jl`, two robustness fixes required by the baseline's state:
  the baseline's PVD and metadata.toml manifests are **stale** (truncated to steps 0–4 by
  a later 5-step debug run into the same directory) while all per-step VTK files (0–359)
  are intact. (a) An explicit `restart_step` missing from the PVD is now accepted if its
  body VTU exists (warns about the stale manifest). (b) When the frame-state manifest
  lacks the step, frames are **reconstructed by kinematic replay** — `maneuver!` at
  `t_range[i+1]` + `propagate_kinematics!` with `dt = t_range[i+2]-t_range[i+1]` for steps
  0..restart_step-1 plus the restart step's `maneuver!`, mirroring `simulate!`'s loop
  exactly; body-node motion during replay is harmless (nodes are overwritten from disk
  right after). The continuation-fidelity gate (§4) implicitly validates this
  reconstruction: a frame-angle error would misalign body vs wake and show up immediately
  as a CT excursion at the first continued steps.
- `examples/rotor_hover_pressure_comparison.jl`: `RESTART_STEP` / `RESTART_PATH` /
  `RESTART_NAME` ENV branch calling `pnl.simulate_warmstart!` with the example's full
  kwarg set. Construction-time knobs (RELAX_RLXF, KERNELOFFSET_*,
  RELAX_FILTER_DOWNSTREAM_R, WAKE_CORE_SIZE) take effect via ENV at include time;
  simulate!-kwarg knobs via the forwarding. Also: explicit
  `LinearAlgebra.BLAS.set_num_threads` at startup (from BLAS_NUM_THREADS /
  OMP_NUM_THREADS, printed to the log) and the `BODY_GRADIENT_KERNELOFFSET` ENV knob.
- `src/FLOWPanel_simulate.jl`: new `body_gradient_kerneloffset::Float64=NaN` kwarg
  threaded through `simulate!` → `_steady_aerodynamics!` (and the steady runner), with the
  two-pass body→particles influence in the `body_on_wake` branch and metadata recording
  under solver_options.
- `examples/rotor_hover_ct_knob_sweep.jl`: the driver. Scenario table, one isolated
  subprocess per scenario (crash = recorded outcome), `rm -rf` of each scenario dir
  before its run, logs at `data/ct_knob_<name>.log`, CT parsing (continuation rows only —
  monitor histories are zero before the restart step), np-collapse detection from
  particle-file sizes, merged summary at `data/ct_knob_sweep_summary.csv` (safe to run
  scenarios incrementally across invocations), ΔCT-vs-control table, automatic λ→0
  Richardson extrapolation over the stable rlxf ladder.

Scenarios (driver `SCENARIO_DEFS`): `control` (pure continuation), `bodyhess`,
`bodyhess_gradoff` (BODY_HESSIAN_TO_PARTICLES=true + BODY_GRADIENT_KERNELOFFSET=4e-3 —
knob 3), `wakerowhess`, `koff_5e-4`, `koff_2p5e-4` (both with WAKE_CORE_SIZE pinned),
`rlxf_0p15`, `rlxf_0p075`, `relaxfilter_0p5R`, `merge_off` (MERGE_PARTICLES=false —
knob 4).

**Thread control (verified):** the driver exports `OPENBLAS_NUM_THREADS`,
`OMP_NUM_THREADS`, `BLAS_NUM_THREADS`, and `JULIA_NUM_THREADS` = THREADS to each
subprocess, and the example additionally calls `LinearAlgebra.BLAS.set_num_threads` from
`BLAS_NUM_THREADS`/`OMP_NUM_THREADS` at startup and prints "BLAS threads: N" to the log —
so BLAS threading is enforced regardless of BLAS vendor or library load order. Check that
line in any scenario log to confirm 48/48.

## 3. Prerequisites on the workstation

- Repo clone including `data/rotor_hover_pressure_comparison/` with body/wake/particle
  VTKs for steps 0–359, the `*_body1.pvd`, and `*.metadata.toml` (warmstart reads these).
- `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
- Disk: each scenario writes ~150 steps of VTK (≈ 40% of the baseline dir's size);
  8 scenarios ≈ 3× the baseline dir. Scenario dirs are overwritten in place on reruns.

## 4. Execution sequence

```bash
# 1. Smoke gate (~11-step continuation; also measures per-step walltime)
SMOKE=true SCENARIOS=control THREADS=48 julia --project=. examples/rotor_hover_ct_knob_sweep.jl

# 2. Continuation-fidelity gate + primary suspect
SCENARIOS=control,bodyhess THREADS=48 julia --project=. examples/rotor_hover_ct_knob_sweep.jl

# 3. Remaining scenarios (incremental; summary CSV merges across invocations)
SCENARIOS=bodyhess_gradoff,wakerowhess,koff_5e-4,koff_2p5e-4 THREADS=48 julia --project=. examples/rotor_hover_ct_knob_sweep.jl
SCENARIOS=rlxf_0p15,rlxf_0p075,relaxfilter_0p5R THREADS=48 julia --project=. examples/rotor_hover_ct_knob_sweep.jl
SCENARIOS=merge_off THREADS=48 julia --project=. examples/rotor_hover_ct_knob_sweep.jl
```

**Gates, in order (stop and report on failure):**

1. **Smoke — PASSED locally 2026-07-02** (8 threads): warmstart branch ran, the missing
   frame-state manifest correctly triggered kinematic replay, 11 continuation steps
   completed, CSV parsed. Measured: 9.9 min total of which ~3–4 min is JIT/include/replay
   startup ⇒ ~35–40 s/step at 8 threads; at 48 threads expect ~10–20 s/step ⇒ **~40–60 min
   per full scenario (154 steps) + ~5 min startup**. Re-run the smoke on the workstation
   anyway (environment check) and record its per-step time.
2. **Continuation fidelity — PASSED at smoke scale**: smoke CT_bernoulli = 0.06223 vs
   baseline settled 0.06238 (Δ well inside the 3.8e-4 ripple). Confirm again on the full
   `control` scenario: final-rev CT_bernoulli within **0.06238 ± ~0.001**, std comparable
   to 3.8e-4. Known artifact: CT_laplace_lamb reads low for the first continued steps (the
   PressureLaplace monitor's ∂u/∂t finite difference needs history — smoke read 0.0593 vs
   baseline 0.0621); over a 4-rev scenario this washes out — judge fidelity on Bernoulli.
   If control drifts materially, the warmstart replay is wrong; debug before interpreting
   ANY knob (check the scenario ENV reproduces the baseline construction: RHPC_MESH=40_40,
   NT=48, RPM=6000).
3. Then knobs, comparing each scenario's final-rev CT (and std/ptp) against control's
   from the same sweep, not against the baseline.

## 5. Interpretation and reporting

- **Stability is data**: a crashed/NaN/np-collapsed scenario is recorded with its failure
  step (from the log and the last written VTK index) — for `bodyhess`/`wakerowhess`/`koff_*`
  that is the expected primary result if the user's instability suspicion holds. Distinguish
  "crashed immediately (steps 1–5)" from "ran revs then diverged" — the latter still yields
  a usable CT window before divergence; report it with the window stated.
- ΔCT direction and size per knob vs control, with the control ripple (std) as the
  significance yardstick: |ΔCT| < 2× control std is noise, do not over-interpret.
- Relaxation: report CT(λ) at λ = {270, 540, 1080} s⁻¹, the fitted order p, and CT(λ→0)
  with an honest error estimate (difference between 2- and 3-point extrapolations). If the
  ladder is non-monotone, report raw values and no extrapolation. Compare the
  `relaxfilter_0p5R` CT against the extrapolated value — agreement is strong evidence for
  a relaxation-free CT figure.
- Deliverable: `data/ct_knob_sweep_summary.csv` (written by the driver) plus a dated
  results note at `data/ct_knob_sweep_notes.md` — findings per knob, failure steps,
  extrapolation, and a recommendation of which knob(s) merit a full cold run. The user
  decides which BRAINSTORM item absorbs the results (004/005/006 lineage).

## 6. Autonomous extensions — you have license to follow the results

The user explicitly authorizes trying new things as results come in; wall-clock compute is
cheap (subprocesses), so bias toward running the informative next scenario rather than
waiting. Sanctioned without asking:

- **New ENV scenarios** in `SCENARIO_DEFS` (pure warmstart perturbations), guided by
  results. Prepared contingencies:
  - `bodyhess` unstable → `bodyhess_gradoff` is the designated rescue (already in the
    table): same physics restored but with the gradient regularized at 4e-3. If
    `bodyhess_gradoff` is ALSO unstable, raise `BODY_GRADIENT_KERNELOFFSET` (1e-2); if it
    is stable, bisect downward (2e-3) toward the least-smoothed stable point. This family
    is the single most interesting follow-up — the strain physics is real (measured 3.2×
    the retained term in the young wake) and gradient-only smoothing is the targeted
    stabilizer (it leaves the advecting velocity untouched, unlike raising
    KERNELOFFSET_TARGETS wholesale).
  - `bodyhess` stable and ΔCT > 2× control std → add `bodyhess+wakerowhess` combined, and
    a longer continuation (SETTLE_REVS=8) to confirm the new CT plateau.
  - `koff_5e-4` stable but `koff_2p5e-4` unstable → bisect (3.5e-4); if a CT-vs-offset
    trend emerges, add one larger point (2e-3) to bracket the curve.
  - rlxf ladder non-monotone or 0.075 unstable → insert 0.2 and/or 0.1 to resolve the
    stable range before extrapolating.
  - `relaxfilter_0p5R` stable → sweep d ∈ {0.25, 1.0} to test the CT plateau.
  - Knob interactions that a single-knob result makes compelling (state the reasoning in
    the notes file).
- **Read-only analysis** of scenario outputs, including reusing the existing offline
  diagnostics on new states (`examples/particle_overlap_age_diag.jl`,
  `examples/particle_body_strain_audit.jl` — e.g. run the strain audit on a bodyhess
  output to verify the term is now consistent), and small plotting/summary scripts.
- Editing the **driver's scenario table and reporting** as needed.

**Second-tier / creative knobs (explore after the primary sweep above).** All three
already exist as ENV knobs in `examples/rotor_hover_pressure_comparison.jl` — no src work
needed to try them:

- **Downstream cull plane** `TRUNCATION_DEPTH_R` (default 4). The example installs a
  `GlobalCylinder([-0.5R,0,0], [TRUNCATION_DEPTH_R·R,0,0], 1.5R)` particle-trim policy
  (`examples/rotor_hover_pressure_comparison.jl:298`) that removes every particle past a plane
  `TRUNCATION_DEPTH_R` radii downstream of the rotor. This is a per-step trim applied at
  construction, so it **takes effect in warm-start continuations** — a legitimate sweep
  scenario. Deeper truncation retains more wake and, like `merge_off`, tests whether
  premature wake removal is suppressing CT (but at fixed radial extent 1.5R). Add e.g.
  `trunc_6R`/`trunc_8R` scenarios; watch np/walltime as with `merge_off`.
- **Run length / revolutions** `SETTLE_REVS` (total = ramp+hold+withdraw+settle;
  `NREVS` floor). Already used to set continuation length (§0 base_env sets 5); longer
  settling tests whether the reported CT is a true plateau. Warm-start-compatible.
- **Spinup and freestream schedule** `SPINUP_REVS`, `SPINUP_START_FRACTION`,
  `FREESTREAM_{RAMP,HOLD,WITHDRAW}_REVS`. These shape the **pre-restart** trajectory
  (steps 0–350) and are therefore baked into the baseline being resumed — they **cannot be
  tested by warm-start perturbation** and require a fresh cold run (gated on user approval,
  below). Item 005 is the lineage for these.

Requires checking with the user first: edits to `src/` or to the example's physics, cold
(non-warmstart) full runs, deleting/overwriting the baseline `rotor_hover_pressure_comparison`
data, and anything estimated to consume a nontrivial token fraction (state the estimate).
Every autonomous addition gets one line of justification in `data/ct_knob_sweep_notes.md`.

## 7. Usage-aware pacing

See §0 — the token-budget operating rules there are binding for the whole session.

## 8. Verification summary

- Smoke gate (driver `SMOKE=true`) passes; per-step walltime recorded.
- Control continuation reproduces baseline settled CT (0.06238 ± ripple).
- Band partition of every scenario's CT metrics uses continuation rows only (rev > 9.72).
- Summary CSV merges correctly across incremental invocations.
- On any unexpected driver error: the scenario log at `data/ct_knob_<name>.log` contains
  the subprocess stdout+stderr; the example's own printed config line ("Particle
  diagnostics: ...") confirms which knobs the subprocess actually received.
