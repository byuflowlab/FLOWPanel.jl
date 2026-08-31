# Phase 2e — Unsteady Hover CT Convergence (velocity vs Green reconstruction)

## Purpose and required context

Obtain a **converged unsteady hover CT** for the DJI 9443 on the Phase 2d
production mesh, using a `PanelParticleWake` with a **single doublet panel wake
row** (`nwakerows=1`), in two solve formulations:

1. **velocity** — the production default `VelocityThroughSources()`
2. **green** — `GreenReconstruction(gauge=:area_mean, recompute_interval=1,
   green_solver=nothing)`

and converge the two parameters most likely to move CT: the **timestep**
(`NT`, steps per revolution) and the **wake truncation depth**
(`TRUNCATION_DEPTH_R`, where downstream particles are deleted).

This phase is inserted after Phase 2d at Ryan's direction (2026-07-25) and
**front-runs elements of Phase 4 (Green gate) and Phase 5 (CT convergence)**.
Mitigation for the un-run Phase 4 gate: preflight `test/formulation_test.jl`
(its Stage 6 checks Green consistency) and always report Green CT
**side-by-side with the matched velocity run** — a Green result is never
reported alone.

Read the repository instructions, the dashboard, the top **Current snapshot** of
the [Phase 2e log](../../logs/dji_convergence_20260722/phase_02e_unsteady_ct.md),
and `agent_policies/HPC.md`, `agent_policies/MONITORS.md`,
`agent_policies/WORKFLOW.md`.

Do not begin the next phase until Ryan explicitly approves the Phase 2e decision.

### Binding study directives

- ≤6 active study Slurm jobs at any time (cap raised from 3 by Ryan,
  2026-07-25); check `squeue` before submitting more.
- ≤4 retained full ParaView/VTK histories across this study.
- Slurm submission: Ryan granted explicit permission on 2026-07-25. Cluster
  access is `ssh orc` (the raw `ssh.rc.byu.edu` host requires 2FA); wrap remote
  `squeue`/`sbatch` in `bash -lc "..."`. Still verify md5s and the 3-job cap
  before every submission.
- **Shedding invariant:** never derive shedding from raw mesh cells. Build a
  `noshedding` body first, run `calc_shedding_from_seed` on *its* `.nodes` /
  `.cells`, then rebuild with that shedding. The driver already does this —
  preserve the pattern.
- Local jobs use at most 6 threads.
- Preserve unrelated local changes and untracked files.
- Convergence and agreement with experiment are **separate claims**.

## Configuration (fixed for every case)

| Item | Value |
|---|---|
| Mesh | `examples/data/dji9443_20260725_45_185_capped_captess4.msh` (flat root + round CapUMinTess=4 tip, 45 span, n_airfoil=185; Phase 2d production recipe) |
| RPM | 5400 (study directive; matches `CT_exp`) |
| rho, R | 1.179 kg/m³, 0.119 m (mesh rescaled from ~0.0597 m) |
| Body | `RigidWakeBody{Union{ConstantSource,VortexRing}}`, `DBC=true`, `watertight=true`, `semiinfinite_wake=false` |
| Body solver | `Backslash` (dense direct) |
| Wake | `PanelParticleWake(nwakerows=1)`, `OverlapPPS(3.0, 2)` trailing, `NoShed()` unsteady, `CoreSpreading` viscous, merge + `GlobalCylinder` truncation |
| Monitors | `BERNOULLI_ONLY=true`: `PressureBernoulli(rho; unsteady=false)` → `ForceMonitor(...; RotorNormalization(rho, 2R, 1))` → `BoundCirculationMonitor` |
| Startup | validated staged recipe: `SPINUP_REVS=1.5 SPINUP_START_FRACTION=0.4 MAGVINF_PEAK=5.0 MAGVINF_END=0.0 FREESTREAM_RAMP_REVS=1.0 FREESTREAM_HOLD_REVS=1.5 FREESTREAM_WITHDRAW_REVS=4 SETTLE_REVS=12` (≈20 revs) |
| HPC | 1 node, 64 tasks/threads, **64 GB**, `--time=24:00:00` (48 h for NT=72) |

**Memory deviation from the original plan (24 GB → 64 GB), measured not
guessed.** The production mesh is triangulated: **36,752 panels**, not the
~16.6k quads the plan estimated. Dense-operator footprint at that size:

| Operator | Size |
|---|---:|
| `Backslash` G (`lu!` in place) | 10.1 GB |
| `GreenReconstruction(:area_mean)` build peak (B + bordered K) | 20.1 GB |
| … retained after build (K/LU only) | 10.1 GB |
| **velocity case, steady** | **~10 GB** |
| **green case, peak** | **~30 GB** |

24 GB would OOM every Green case. `--mem=64G` covers both with headroom for the
particle wake and FMM trees; a velocity-only case may be submitted with
`sbatch --mem=24G`.

**Caveat to carry into the report:** `PressureBernoulli(unsteady=false)` on a
moving body was flagged in a past regression (inertial-KE term). It is retained
here because it is exactly what produced the accepted 6000-RPM `CT ≈ 0.062`
plateau, so the comparison stays like-for-like. Note it explicitly.

Do **not** add `PressureLaplace` to routine runs (no formulation has passed the
full unsteady acceptance gate). `BERNOULLI_ONLY=true` already complies.

## Driver

`examples/rotor_hover_pressure_comparison.jl` — one env-driven driver. Phase 2e
additions:

- `RHPC_MESH` now accepts `45_185_ct4` (and honors the env — a previous
  hard-coded override of the 40_40 mesh was removed).
- Trailing-edge seed nodes are detected automatically with
  `find_dji9443_trailing_edge_indices` (`examples/dji9443_trailing_edge.jl`) for
  any mesh other than the legacy hard-coded ones.
- `RHPC_FORMULATION` ∈ {`velocity` (default), `green`} selects the
  `formulation=` kwarg to `simulate!`.
- The metadata TOML records mesh key, mesh file, and formulation.

Launcher: `examples/run_dji9443_hover_ct_hpc.slurm.sh` (case loop; one fresh
Julia process per case).

## Run matrix

**Batch 1 (3 jobs):**

| Case tag | `RHPC_FORMULATION` | `NT` | `TRUNCATION_DEPTH_R` |
|---|---|---:|---:|
| `p2e_vel_nt36_d4`   | velocity | 36 | 4 |
| `p2e_green_nt36_d4` | green    | 36 | 4 |
| `p2e_vel_nt72_d4`   | velocity | 72 | 4 |

For `NT=72`, per-step relaxation is rescaled to preserve the physical rate
(Phase 5 prescription; prior NT=72 recipe: `examples/run_006_nt72_rlxf_sweep.sh`).

**Batch 2 (adaptive, chosen after inspecting Batch 1):**

- `p2e_vel_nt36_d6` (or `d8`) — truncation-depth probe.
- Green confirmation of whichever knob moved CT by ≥1% (matched perturbation).
- Reserve slot: further NT refinement, warm-start settle extension
  (`RESTART_STEP`/`RESTART_NAME`/`RESTART_PATH` → `simulate_warmstart!`), or a
  `KrylovSolver` Green cost check.

## Convergence criterion

**Per-run convergence** (Ryan: "CT settles to a small enough amplitude with a
mean that changes little over 5 revolutions"), over the final **5 complete
self-sustained-hover revolutions** (freestream fully withdrawn):

- (a) the 5 per-rev CT means all lie within **±0.5%** of their average;
- (b) within-rev peak-to-peak CT ≤ **2%** of the mean;
- (c) all state finite (no NaN/Inf in CT, no particle blow-up).

If (a)–(c) fail, extend the settle window by warm start (+5 revs at a time, up
to ~17–20 settle revs) before declaring non-convergence.

**Parameter convergence:** CT mean shift ≤ **1%** between a rung and its
refinement (`NT` 36→72; depth 4→6), evaluated **separately per formulation**.

## Decision rule and deliverables

Phase 2e closes when, for each formulation:

1. a per-run-converged case exists at the finest rung reached, and
2. the NT and truncation-depth refinements each move CT by ≤1% (or the
   non-convergence is quantified and its cause named).

Deliverables under `data/dji_convergence_20260722/phase_02e_unsteady_ct/`:

- per-case CT-vs-revolution CSVs and metadata TOMLs,
- a per-rev block-statistics table (mean, peak-to-peak, drift) per case,
- plateau plots (velocity vs green, per knob),
- `phase_02e_report.md`: CT velocity vs Green per knob; per-run and parameter
  convergence status; comparison against `CT_exp = 0.072` at 5400 RPM, the
  steady Dirichlet baselines (≈0.0505–0.0515, semi-infinite rigid wake), and
  the 6000-RPM unsteady point (≈0.0617–0.0622) — stating clearly that
  convergence and experimental agreement are separate claims; Green cost and
  residual notes; VTK retention ≤4.

## Key files

- Driver: `examples/rotor_hover_pressure_comparison.jl`
- Launcher: `examples/run_dji9443_hover_ct_hpc.slurm.sh`
- TE detection: `examples/dji9443_trailing_edge.jl`
- Formulations: `src/FLOWPanel_formulation.jl`; kwarg in `src/FLOWPanel_simulate.jl`
- Formulation tests: `test/formulation_test.jl`
- Oscillation analysis: `examples/analyze_stable_wake_oscillation.jl`
- Prior stable-wake recipes: `examples/run_rotor_hover_stable_wake.sh`,
  `examples/run_006_nt72_rlxf_sweep.sh`
- Log: `logs/dji_convergence_20260722/phase_02e_unsteady_ct.md`

## Verification

- [x] `test/formulation_test.jl` passes locally (all 10 stages, 2026-07-25).
- [x] Local smoke: both formulations run end-to-end and give finite CT on the
      **40_40** mesh; `RHPC_MESH` env respected; per-rev CSV + metadata written.
      The ct4 mesh itself cannot be smoked locally (10 GB dense G, 30 GB for
      Green, vs 17 GB of laptop RAM), so it was verified through a
      construction-only preflight: auto TE seeds, 36,752 cells, correct body
      type, **0 unpaired shedding edges**, shedding roots at r/R = ∓0.111.
- [x] Cluster md5 checks (mesh + every modified file) recorded in the log
      before submission — all seven match.
- [ ] Each HPC case judged against the 5-rev criterion; knob shifts ≤1%.
- [ ] Every report claim backed by a CSV/log path; ledgers current.

## Checkboxes

- [x] Phase docs + dashboard entry
- [x] Driver edits (mesh env, auto TE seeds, `RHPC_FORMULATION`, metadata)
- [x] Slurm launcher
- [x] Local preflight (formulation tests + both-formulation smoke)
- [x] Deploy + md5 verification + submission (Batch 1: 12894164/12894165/12894166)
- [ ] Batch 1 analysis
- [ ] Batch 2 (adaptive)
- [ ] Decision report + dashboard/log update; stop for Ryan's approval
