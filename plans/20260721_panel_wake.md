# DJI 9443 Direct-Potential Panel-Wake Convergence and HPC Plan

## Context (why this work)

FLOWPanel's steady rigid-wake panel method predicts DJI 9443 hover CT ≈ 0.0505, well below the experimental
≈ 0.072 and viscous-BEM ≈ 0.068; prior investigations route the shortfall to **wake modeling**, and a
corrected particle-wake Bernoulli replay landed near CT 0.05068. This plan implements a new **direct-potential
finite panel-wake** solve formulation (`DirectWakePotential`, the production version of the Task-3 diagnostic)
and runs an HPC study to test whether a *length-converged* scalar-potential panel wake moves panel-method CT
toward the BEM/experimental values — in an axial validation case (vs CCBlade) and in the target hover case.
The deliverable is a defensible statement separating solver correctness, wake-row (length) convergence, and
physical agreement.

This document is self-contained: a fresh agent should be able to execute it without additional codebase
exploration. Read `CLAUDE.md` and the routed policies first (see below), then this plan, then
`plans/20260721_panel_wake_progress.md` (current status + remaining-work checklist).

## STATUS (2026-07-21) — read before doing anything

All local code deliverables are **DONE and verified**. Do not re-implement, re-test, or re-smoke them.

- **`DirectWakePotential` formulation** — implemented in `src/FLOWPanel_formulation.jl`, exported in
  `src/FLOWPanel.jl`. `test/formulation_test.jl` Stages 8/9/10 all pass (Task-3 single-shot equivalence
  rel residual 8.6e-16; all invalid-config rejections; ConstantDoublet↔VortexRing q_f agrees 1.7e-13,
  γ 3.4e-14). Note: the line numbers quoted in the "Solve-formulation subsystem" section below describe
  the file *before* this implementation; treat them as orientation, not current offsets.
- **`examples/rotor_panel_wake_study.jl`** — written, complete (monitors, derived NREVS, CT-vs-rev CSV,
  metadata TOML, COMPLETED marker, final-rev stability stats).
- **Local smoke test** — PASSED. Output in `data/dji9443_panelwake_smoke/`: COMPLETED marker,
  `all_finite=true`, 12 finite steps, VTK + metadata + all monitor CSVs. `stable=false` there is expected
  (1-rev smoke config cannot satisfy the final-rev gates). The plan's "direct/manual potential diagnostic
  agrees" smoke clause is discharged by test Stage 8, not by an in-smoke check. Wall time: 204 s for 12
  tiny steps locally (NT=12, 3 rows, 4 threads) — a lower-bound data point for `--time` estimation.
- **`examples/run_dji9443_panel_wake_convergence.jl`** — written, complete (subprocess-per-case,
  completion-marker checkpointing, all three stages, time-budget awareness, aggregate report).
- **`scripts/dji9443_panel_wake_convergence.slurm.sh`** — written per the launcher template below.

**Remaining work:** HPC deploy → record + submit sbatch → monitor/resubmit → reports + reviewed
conclusion. The step-by-step remaining checklist lives in `plans/20260721_panel_wake_progress.md`.

**Keep the progress log current:** as you complete (or partially complete) each checklist item, update
`plans/20260721_panel_wake_progress.md` immediately — move the item to DONE with its concrete results
(commands run verbatim, Slurm job IDs, per-case wall times, output locations, any deviations from this
plan and why), so the next fresh-context agent can resume without re-deriving anything. That log is the
single source of truth for execution state; this plan is the source of truth for method.

## Mandatory Policy Reads (before any work)

- `CLAUDE.md` (repo routing hub + critical invariants).
- `agent_policies/WORKFLOW.md` (subsystem routing, editing boundaries).
- `agent_policies/TESTING.md` (verification matrix, smoke tests).
- `agent_policies/MONITORS.md` (monitor dependency contracts, PressureLaplace constraints).
- `agent_policies/HPC.md` (Slurm launcher pattern, inputs/checkpoints, submission boundary).

## Repository Orientation (concrete facts — so no exploration is needed)

### Solve-formulation subsystem — `src/FLOWPanel_formulation.jl` (898 lines)

Existing formulations and the exact pattern a new one follows:

- `abstract type AbstractSolveFormulation` — line 25.
- `VelocityThroughSources` (default, empty singleton) — struct line 34; hooks are no-ops.
- `GreenReconstruction{TS}` — struct line 66, kwarg constructor line 71 (validates `gauge`); state struct
  `GreenReconstructionState{TF,TG}` line 241. **This is the closest production-hook precedent** for the new
  formulation: it builds a dense finite-body (bordered) LU once and reuses it.
- `TraceCorrected{TS}` — struct line 112, constructor line 122; state `TraceCorrectedState` line 251.
  **Not the model to copy** — it is Task 5 (velocity trace correction) and is pinned to `DirectBackend`.

The hooks a new formulation must implement (all in this file):

1. **Constructor + kwarg validation** (cf. lines 71, 122).
2. **`initialize_formulation(f, systems_tuple, wakes_tuple, body_solvers, backend_solve, backend_system)`** —
   docstring lines 299–308; dispatch examples: `VelocityThroughSources` returns `nothing` (line 309),
   `GreenReconstruction` (line 312), `TraceCorrected` (line 324). Called **once** before the time loop.
   Shared validation helper **`_validate_formulation_common`** (line 273) enforces: a single Dirichlet
   `RigidWakeBody` with `NK==2`, `semiinfinite_wake==false`, a `Backslash` body solver, and a non-empty wake.
   `TraceCorrected` additionally guards `backend_system isa DirectBackend` at lines 328–331 — **the new
   formulation must NOT add this guard** (it runs on FMM).
3. **`formulation_prewake!(f, state, systems_tuple)`** — line 663; a generic fallback (lines 664–665)
   snapshots pre-wake body velocity into `state.u_prewake`, so the new formulation can rely on it.
4. **`solve_formulation!(f, state, systems, systems_tuple, wakes_tuple, body_solvers; backend_solve,
   backend_wake, i_step)`** — dispatch examples: `VelocityThroughSources` (line 674),
   `GreenReconstruction` (line 681), `TraceCorrected` (line 714). This is the per-step solve.

### How `simulate!` wires formulations — `src/FLOWPanel_simulate.jl`

- `formulation` kwarg (default `VelocityThroughSources()`): outer `simulate!` line 695, inner stepping
  function line 332; companion `formulation_state=nothing` line 333.
- `initialize_formulation(...)` called once before the time loop at **line 739**.
- `formulation_prewake!(...)` at **line 363**.
- `solve_formulation!(...)` at **line 454** (args passed at 787–788).
- `_scalar_potential_sources(w::PanelWake) = get_sources(w)` at **line 269** (particle variant 270,
  aggregation 276) — this is how `simulate!` already treats a `PanelWake` as a scalar-potential source.
- A new formulation is wired in automatically through these kwargs; **no changes to the invocation sites are
  needed**.

### PanelWake API — `src/FLOWPanel_wake.jl`

- `struct PanelWake{TK,NK,TF} <: AbstractFreeWake` — line 92; field `include_final_filament::Bool` at line
  102 (**defaults `true`**; kwarg constructor lines 124–127; also `core_size`, `nwakerows`).
- `get_sources(wake)` (lines 119–122) returns `(wake, FilamentWrapper(wake))` when `include_final_filament`
  is true, else `(wake,)`. The trailing semi-infinite filament is vector-potential-only → the new formulation
  must **reject `include_final_filament=true`**.
- **Scalar potential IS supported via FMM**: `FastMultipole.direct!` for `PanelWake` at line 444 computes
  per-target potential and calls `FastMultipole.set_scalar_potential!(...)` at line 467 when the `PS` switch
  is set. There is no separate `scalarpotential` function; potential is delivered through this `direct!` path.
- Replay round-trips the flag: `src/FLOWPanel_replay.jl` writes `include_final_filament` at line 309,
  reconstructs at line 467.

### Winding-safe RigidWakeBody construction (CRITICAL — silent failure otherwise)

Per the `RigidWakeBody` docstring in `src/FLOWPanel_liftingbody.jl` and `CLAUDE.md`: with
`ensure_winding=true` (default) the constructor re-winds `cells` in place. Compute shedding from the
**constructed** body's `.nodes`/`.cells`, never the raw mesh, or the wake attaches at wrong edges with **no
error** (observed rotor-hover CT collapse ~3.6×). Procedure: (1) build a no-shedding `RigidWakeBody`;
(2) run `calc_shedding_from_seed` on *its* rewound `.nodes`/`.cells`; (3) rebuild with the derived shedding.
Working example: `examples/rotor_hover_convergence.jl`.

### Example templates to copy from (do NOT modify these)

- `examples/rotor_hover.jl` — `PanelWake(rotor; nwakerows=..., ...)` at line 123 (closest PanelWake+rotor
  template).
- `examples/rotor_hover_prescribed_helical_wake.jl` — fills a `PanelWake(rotor; nwakerows=..., ...)` at line
  148 (Bagai/Leishman-style relaxation reference).
- `examples/rotor_hover_pressure_comparison.jl` — commented `PanelWake` at line 273, active
  `PanelParticleWake` at 286; monitors + pressure recovery reference.
- `examples/rotor_axial_j0187_panel.jl` — thin axial wrapper that `include`s
  `rotor_hover_pressure_comparison.jl` (line 26) and overrides freestream to DJI 9443 J=0.1867. Do not modify;
  the new driver is standalone.
- Related axial post-processing: `rotor_axial_j0187_ccblade.jl`, `rotor_axial_j0187_replay.jl`,
  `rotor_axial_j0187_bernoulli_replay.jl`, `rotor_axial_j0187_loading_comparison.jl`.

### Task-3 diagnostic (the formulation's oracle) — `debug/dirichlet_solve/`

- Index: `debug/dirichlet_solve/dirichlet_solve.md`; item file: `debug/dirichlet_solve/task3.md`.
- Theory gate: `docs/wake_solve_schemes.md`, section "Discrete diagnostic construction: Task 3 direct
  fixed-wake potential".
- The diagnostic solves `G_Δ μ_E = −S σ₀ − q_f` with `σ₀ = −U∞·n`; `q_f` evaluated directly from the
  prescribed constant-doublet `PanelWake` at exact body centroids; **wake velocity NOT converted to body
  sources**; **no mean removed** from `q_f`; the finite-body LU applied **manually** because single-body
  `solve!` clears the preassembled potential. Augmented-oracle cross-check: `(G_Δ + P·C) μ* = −S σ₀`.
- Reproduce commands (produce the oracle CSVs/lifts):
  ```
  JULIA_NUM_THREADS=2 julia --project=. debug/dirichlet_solve/dirichlet_solve.jl task3-flat
  JULIA_NUM_THREADS=2 julia --project=. debug/dirichlet_solve/dirichlet_solve.jl task3-march
  ```
  Documented residual gate ≤10⁻¹⁰ on `DirectBackend`; frozen flat cases recover the semi-infinite baseline
  lift within ≈8×10⁻⁵ relative.

### Data and mesh assets (all confirmed present locally)

- Meshes under `examples/data/`: `dji9443_new_40_40.msh` (the study mesh), plus `dji9443_56_57.msh`,
  `dji9443_80_81.msh`, `dji9443_40_40.msh`, `dji9443_50_57.msh` (unused by this plan).
- Steady mesh-convergence reference: `data/rhc_sweep_baseline/convergence_history.csv` (+ per-mesh body dirs
  and `.pvd`).
- CCBlade axial checkpoint: `data/rotor_axial_j0187_ccblade/` (contains `analysis/`, `bernoulli_replay/`,
  polar + CT-vs-J CSV/PNG). Treat the two tagged sectional CSVs, the polar CSV, and the validation report as
  the complete XFOIL/CCBlade checkpoint; rebuild only with `FORCE_CCBLADE=1`.

## Execution Boundary — session-scoped sbatch authorization

- HPC checkout is fixed at `/home/rander39/projects/FLOWPanel.jl`.
- **`agent_policies/HPC.md` (Submission Boundary, lines 110–114) normally forbids agents from running
  `sbatch`/`srun`/`salloc`.** For **this session only**, the user authorizes the agent to run `sbatch` itself
  (via `ssh orc '... sbatch ...'`). This does **not** change the repository policy. The agent still prepares,
  deploys, syntax-checks, and **records the exact submission command** before submitting, and reuses that
  same command verbatim for any resubmission.

## Established Study Parameters

- **Mesh — 40×40 only.** Mesh independence is taken as established by the prior *steady* study (40×40 CT
  0.0504805 vs 80×81 CT 0.0506065 → 0.25%; 56×57 the non-monotonic outlier), recorded in
  `data/rhc_sweep_baseline/convergence_history.csv`. This panel-wake study varies *wake length*, orthogonal
  to mesh resolution, so re-running 80×81 inside the far more expensive unsteady loop is omitted deliberately.
- 40×40 geometry: `examples/data/dji9443_new_40_40.msh`; trailing-edge seeds **`1614,1574,45`** and
  **`3324,3284,1755`** (zero-based ParaView IDs — convert to Julia 1-based indexing).
- Retain default/quad-growing surface-gradient behavior; do **not** use the pathological `quad_nogrow`
  variant.
- Operating geometry: R=0.119 m, ρ=1.179 kg/m³, rotor axis `[-1,0,0]`, axial dim 1, radial dim 2.
  `kerneloffset_panel = R*1e-10`, `kerneloffset_targets = 1e-3`, `kernelcutoff = R*1e-13`.
- CCBlade axial targets at 4 m/s, 5,400 RPM (J=0.18674136): ncrit=4 → CT 0.0419673; ncrit=9 → CT 0.0325479.
- The worktree is dirty with unrelated untracked Dirichlet work. Preserve all unrelated changes; never clean
  or reset the checkout.
- A `PanelWake` is a rolling row buffer; when full its oldest row shifts out. `include_final_filament=false`
  gives a strictly finite, scalar-potential-only wake, hence the explicit wake-length convergence ladder.

## Public Solver/Formulation Change — `DirectWakePotential` (productionize Task 3)

Add and export `DirectWakePotential(; recompute_interval=1)` as a subtype of `AbstractSolveFormulation`.
It is the production version of the Task-3 direct fixed-wake potential diagnostic — **not** `TraceCorrected`.

**Governing equation (each step, on the marched wake):**
```
G_Δ μ_E = −S σ₀ − q_f ,      σ₀ = −U∞ · n
```
`q_f` = scalar potential evaluated directly from the constant-doublet `PanelWake` panels at exact body
centroids (via the FMM `direct!`/`PS` path). Wake velocity is **not** converted into body sources; **no mean
is removed** from `q_f`; the finite-body LU is applied **manually** (single-body `solve!` clears preassembled
potential — see `task3.md` lines 16–19).

Implementation (following the `GreenReconstruction` production-hook pattern):

- Type + docstring alongside the other formulations (`src/FLOWPanel_formulation.jl`, lines ~25–160).
- A state struct holding reusable buffers (cf. `GreenReconstructionState` line 241) for: pre-wake body
  velocity; non-wake source strength σ₀; wake-induced source contribution; direct wake potential `q_f`
  (= `q_wake`); source potential `Sσ₀`; linear-system RHS; last-recompute step.
- `initialize_formulation(::DirectWakePotential, ...)` (cf. lines 309–351): call `_validate_formulation_common`
  (line 273); build/hold the body-only dense LU factorization **once** here (reused every step, like
  `GreenReconstruction`'s bordered LU); validate that every active wake source exposes scalar potential.
- `formulation_prewake!(::DirectWakePotential, state, systems_tuple)` (cf. line 663): the pre-wake body
  velocity snapshot (may reuse the generic fallback that fills `state.u_prewake`).
- `solve_formulation!(::DirectWakePotential, ...)` (cf. line 674): the per-step Task-3 **single-shot** solve
  (NOT the diagnostic's frozen-geometry strength-projection iteration — that is a test-only verification).

Per step:
1. Snapshot body velocity immediately before wake influence (`formulation_prewake!`).
2. Let the standard wake pass add wake velocity to the body.
3. Split body source strength into non-wake σ₀ and wake-induced parts.
4. Evaluate only the free panel wake's scalar potential `q_f` at body control points via `backend_wake`.
5. Evaluate `Sσ₀` via the body solve backend.
6. Solve `G μ = −Sσ₀ − q_f` with the already-factored body matrix.
7. Store σ₀ (not the wake-induced contribution) in the body source column.
8. Restore any temporarily-used body-potential workspace.
9. Clear incompatible affine wake corrections.

**Validation (initialize_formulation must reject clearly):**
- Require exactly one `RigidWakeBody` with Dirichlet BCs, source+doublet storage, `semiinfinite_wake=false`,
  a `Backslash` body solver, and ≥1 free wake (via `_validate_formulation_common`).
- Require all active wake sources to expose scalar potential.
- Reject `PanelWake` with `include_final_filament=true`.
- Reject particle wakes / mixed sources that cannot supply a complete scalar potential.
- Reject Neumann, multi-body, semi-infinite, and non-`Backslash` configurations.

**Backend — FMM (per user direction).** Production runs use the FMM backend (order-4 wake, order-8 body).
Do **not** add the `DirectBackend` guard that `TraceCorrected` uses (Task 3 needs no trace estimation, only a
scalar potential the FMM `PanelWake` path already supplies). FMM approximation error in `q_f` is accepted and
**not** a gating concern. The finite-body `G_Δ` LU is body-only and backend-independent.

**Other requirements:**
- `recompute_interval=1` mandatory for all production runs; larger intervals reserved for later sensitivity
  work only.
- Use source + `ConstantDoublet` for the study body and wake so every wake source has scalar potential. Add a
  regression showing equivalence with the existing source + `VortexRing` representation on a small frozen
  finite-wake case **before** changing the production kernel.
- Do **not** alter defaults/behavior of `VelocityThroughSources`, `GreenReconstruction`, or `TraceCorrected`.

Implementation locations: `src/FLOWPanel_formulation.jl` (type/state/validation/init/solve hooks);
`src/FLOWPanel.jl` (export next to the other formulations); `test/formulation_test.jl` (regression coverage).

## Reusable DJI Panel-Wake Driver — `examples/rotor_panel_wake_study.jl`

Create as a standalone driver (do not modify `examples/rotor_axial_j0187_panel.jl`). Build the body with the
winding-safe procedure above.

Configuration interface (env vars):
- `RUN_NAME`; `RHPC_MESH=40_40` (only supported value); `RPM` (default 5400); `NT` (default 36); `WAKE_ROWS`;
  `SETTLE_REVS` (default 3); `START_VINF` (nonzero-J freestream the run begins at, default 4.0);
  `TERMINAL_VINF` (freestream held at the end); `FREESTREAM_DECREASE_REVS` (revolutions over which freestream
  is linearly decreased `START_VINF`→`TERMINAL_VINF`; `0` = constant freestream, as in the axial study);
  `FORMULATION=direct|velocity|green`; `SAVE_VTK` (default true); plus existing FMM tuning env vars.
- **No RPM spinup, no freestream up-ramp/hold/withdraw.** RPM is constant; the only schedule is a single
  monotone freestream *decrease*.
- **`NREVS` is derived, not fixed:** `NREVS = ceil(WAKE_ROWS/NT) + FREESTREAM_DECREASE_REVS + SETTLE_REVS`.
  (Allow an `NREVS` override for debugging only.) This fixes a latent bug in the prior fixed-`NREVS=12`
  scheme, where top wake-row rungs (e.g. 432 rows = 12 revs just to fill) had zero settled revolutions to
  measure; it also makes short rungs cheap.
- **Schedule ordering (as implemented — deliberate):** the freestream ramp starts at **t=0** (linear
  `START_VINF`→`TERMINAL_VINF` over `FREESTREAM_DECREASE_REVS` revs, then hold); there is *no* leading
  hold at `START_VINF`. With the derived `NREVS`, the run therefore spends
  `ceil(WAKE_ROWS/NT) + SETTLE_REVS` revolutions at the terminal freestream *after* the ramp ends, so the
  final measurement revolution's entire rolling buffer was shed at the terminal condition and has fully
  flushed the ramp-era rows. A fill-first ordering (fill at `START_VINF`, then ramp, then settle) would
  instead leave only `SETTLE_REVS` terminal revolutions — for caps above `SETTLE_REVS·NT` rows the measured
  wake would be mostly axial/ramp-era, contaminating the cap-vs-CT ladder with a cap-dependent bias. Do not
  "fix" the driver to fill first.

Wake and solver setup:
- `PanelWake(rotor; nwakerows=WAKE_ROWS, core_size=1e-3, include_final_filament=false)`.
- Keep `shed_with_induced_velocity=true` unless the smoke test exposes a startup failure.
- `Backslash(rotor)` body solver.
- Body solve backend: existing tuned FMM (expansion order 8, acceptance 0.4, leaf size 20).
- Wake evaluation backend: FMM expansion order 4, acceptance 0.4, leaf size 50.
- `DirectWakePotential(recompute_interval=1)` for production (`FORMULATION=direct`).

Monitors (see `agent_policies/MONITORS.md` for dependency contracts):
- `PressureBernoulli(rho; unsteady=true, allow_partial=false)`.
- Pressure-based `ForceMonitor` with `RotorNormalization`.
- `KuttaJoukowskiForce` — independent total-load cross-check.
- `BoundCirculationMonitor`.
- Spanwise pressure-force loading on fixed radial bins (in-sim or via replay).
- **No `PressureLaplace`** in production convergence runs (adds Hessian/CG cost; not needed here).

Outputs per run: VTK (body, wake, frame state); monitor CSVs; CT vs revolution; final-revolution CT mean,
relative peak-to-peak, relative linear drift; final-revolution mean dCT/d(r/R); formulation name+options;
wake-row cap and active row count; mesh/RPM/NT/`NREVS`/freestream schedule/backend/wall-time; finite-value
checks for body & wake strengths, wake nodes, velocities, pressure, forces; a **completion marker** written
only after all validation and reporting stages succeed.

## Convergence Definitions

A single run is **stable** only if: it completes its (derived) timesteps; all body/wake state and monitored
values stay finite; its final complete revolution has relative CT peak-to-peak ≤5% and relative linear CT
drift ≤2.5%; and its spanwise-loading output is complete (computed in-sim by the driver's spanwise monitor;
a replay pass is optional and NOT required per rung).

For adjacent wake-row caps a, b:
- `CT change = abs(CT_b − CT_a) / max(abs(CT_b), eps())`
- `loading change = norm(L_b − L_a) / max(norm(L_b), eps())`, where `L` is the final-revolution mean
  pressure-derived dCT/d(r/R) on identical fixed radial bins over r/R ≥ 0.1.

Wake length is accepted when **two consecutive** cap refinements both satisfy CT change ≤2% and loading
change ≤2%, with all three involved runs passing the stability gate. CCBlade agreement is reported separately
against ncrit=4 and ncrit=9; it is **never** substituted for numerical convergence.

## Local Verification and Smoke Test (run from repo root)

Focused tests:
```
julia --project=. test/formulation_test.jl
julia --project=. -e 'include("test/runtests_unit_wake.jl")'
julia --project=. -e 'include("test/runtests_unit_liftingbody.jl")'
julia --project=. -e 'include("test/runtests_unit_simulate.jl")'
julia --project=. -e 'include("test/runtests_unit_postprocess.jl")'
```

Required new test scenarios (in `test/formulation_test.jl`):
- **Task-3 oracle equivalence:** run the production `DirectWakePotential` on the same frozen finite-wake state
  the Task-3 diagnostic uses and reproduce its single-shot result to ≤10⁻¹⁰ residual on `DirectBackend`
  (oracle via `debug/dirichlet_solve/dirichlet_solve.jl task3-flat` / `task3-march`; augmented check
  `(G_Δ+P·C)μ* = −Sσ₀`).
- Formulation correctness is proven on `DirectBackend`; production uses FMM. A **non-blocking** sanity check
  may record the FMM-vs-direct `q_f` difference, but FMM error is not a gating criterion.
- `ConstantDoublet` vs `VortexRing` frozen-wake velocity/circulation agree on the small test geometry.
- Invalid configs fail clearly: `include_final_filament=true`, particle wake, Neumann, multi-body,
  semi-infinite, non-`Backslash`.
- Recompute-interval behavior is deterministic.
- Formulation state does not corrupt body potential, source strengths, geometry, wake geometry/strengths.
- Unsteady Bernoulli receives a complete scalar-potential history.

Smoke command (6 steps with a 3-row buffer deliberately exercises wake overflow):
```
RUN_NAME=dji9443_panelwake_smoke RHPC_MESH=40_40 RPM=5400 NT=12 WAKE_ROWS=3 \
SETTLE_REVS=0 START_VINF=4.0 TERMINAL_VINF=4.0 FREESTREAM_DECREASE_REVS=0 \
FORMULATION=direct SAVE_VTK=true \
julia --project=. -t 4 examples/rotor_panel_wake_study.jl
```
(With those settings `NREVS = ceil(3/12)+0+0 = 1` → ~12 steps, still overflowing the 3-row buffer.) Smoke
passes only if tests pass, the run completes, all fields are finite, and VTK+metadata exist; the
direct/manual potential-agreement requirement is discharged by test Stage 8 (Task-3 equivalence), not by an
in-smoke check.

**Status: PASSED** — tests all green and the smoke run completed (`data/dji9443_panelwake_smoke/`:
COMPLETED marker, `all_finite=true`, 12 finite steps, all monitor CSVs; `stable=false` is expected at 1
rev). Do not re-run unless the driver or formulation changes.

## Axial J=0.1867 Study (lightweight formulation validation)

Fixed: 40×40, 5,400 RPM, NT=36, constant axial freestream 4.0 m/s (`START_VINF=TERMINAL_VINF=4.0`,
`FREESTREAM_DECREASE_REVS=0`), `DirectWakePotential`, VTK on. `NREVS` derived per rung.

Row ladder `72, 144, 216, 288, 360, 432` (= 2,4,6,8,10,12 revs of retained wake at 36 rows/rev).

Execution: run 72, 144, 216; evaluate stability + both convergence metrics; continue to 288, 360, 432 only
until two consecutive refinements pass. Axial wake advects downstream and is expected to converge in the first
2–3 rungs — the ladder *demonstrates* convergence, it is not meant to grind. If 432 still fails to produce two
consecutive passing refinements, report wake-length nonconvergence (do not cherry-pick an intermediate result).

Reports: CT histories + final-cycle stats; direct-potential vs Kutta–Joukowski CT; spanwise loading vs CCBlade
ncrit=4 and ncrit=9; integrated CT bias vs each CCBlade model; wake-row convergence table + plots.

## Hover Study (the target case)

Fixed: 40×40, 5,400 RPM (constant), NT=36, `DirectWakePotential`, VTK on. `NREVS` derived per rung as
`ceil(WAKE_ROWS/NT) + FREESTREAM_DECREASE_REVS + SETTLE_REVS`.

**Startup — ease into hover from the nonzero-J operating point.** Each hover run begins at *exactly* the
nonzero-J axial condition (`START_VINF=4.0` m/s, same 5,400 RPM as the axial study) and the freestream is
**gradually and monotonically decreased** from 4.0 m/s to `TERMINAL_VINF` over `FREESTREAM_DECREASE_REVS`
revolutions **starting at t=0**, then held for the rest of the run (see "Schedule ordering" above: the
post-ramp hold spans `ceil(WAKE_ROWS/NT)+SETTLE_REVS` revs, so the measured buffer is entirely
terminal-shed). That is the entire schedule: constant RPM, one freestream ramp-down — no RPM spin-up, no
freestream up-ramp/hold/withdraw. Shedding the early wake under nonzero axial advection and easing the
freestream down lets the wake transition smoothly to the hover-dominated, self-induced-downwash regime.
Each hover run is self-contained (no cross-run warm-start or extra state-loading path needed).

**Ladder once, confirm once:**
1. **Primary ladder on the 0.25 m/s stabilized case** (`TERMINAL_VINF=0.25`, J≈0.0117). Start the wake-row
   ladder at half / equal / 1.5× the axial converged cap, rounded to multiples of 36; apply the same
   two-consecutive-refinement 2% CT/loading rule; if needed, extend up to 432 in steps equal to the
   ladder's initial spacing (0.5× the axial cap) — **not** finer 36-row steps, which would shrink the
   refinement increment and make the 2% change gate systematically easier to pass than the rungs it is
   compared against. This case is the
   ladder target because a small terminal freestream damps the intrinsic non-damping hover oscillation seen in
   prior work.
2. **Zero-flow (0.0 m/s) confirmation — single run at the accepted cap** (`TERMINAL_VINF=0.0`). Same recipe:
   begin at 4.0 m/s and ramp the freestream all the way to 0, letting wake-induced downwash sustain the
   structure once the freestream is gone. Run **only** at the wake-row cap already accepted on the 0.25 m/s
   ladder (no second full ladder). If it holds the stability gate at 0 freestream, report zero-flow as the
   headline hover result. If it destabilizes as the freestream approaches 0, report the 0.25 m/s result as the
   practical hover case (explicitly labeled J≈0.0117, not exact hover) and document the freestream value at
   which it destabilizes. Preserve both runs as a freestream-sensitivity comparison.

**Freestream-decrease rate is a tunable that likely needs adjustment.** How gradually the freestream is
relaxed (`FREESTREAM_DECREASE_REVS`) controls whether the wake stays stable through the transition — too fast
a ramp-down can shock the wake into the non-damping oscillation, especially approaching 0. Start with ~3–4
revs; if a hover run destabilizes during or shortly after the ramp-down, **increase `FREESTREAM_DECREASE_REVS`**
(slower relaxation) and re-run before declaring the case physically nonconvergent. The zero-flow case may need
a gentler ramp than 0.25 m/s. Record the decrease-rate used for each accepted result (`NREVS` scales with it
automatically).

Hover report: terminal freestream + advance ratio; final-cycle CT/loading convergence; direct-potential vs
Kutta–Joukowski CT; comparison with references (experiment ≈0.072, viscous BEM ≈0.068, particle wake ≈0.062,
steady rigid wake ≈0.0505); explicit separation of numerical wake-row convergence from physical agreement.

## Adaptive Driver and Checkpointing — `examples/run_dji9443_panel_wake_convergence.jl`

Stage order (no 80×81 stages):
1. Axial 40×40 ladder (light).
2. Hover 0.25 m/s 40×40 ladder.
3. Hover zero-flow 40×40 single confirmation at the accepted cap.
4. Aggregate report.

Responsibilities:
- Before each case, check for a complete validation report + completion marker; skip valid completed cases.
- Never treat partial VTK/CSV or a bare process exit as completion.
- Record per-case elapsed time; estimate whether the next case fits before the Slurm safety reserve; stop
  cleanly ≥15 min before allocation expiry; on resubmission resume at the first incomplete case.
- Preserve existing CCBlade checkpoints under `data/rotor_axial_j0187_ccblade/`.
- Stable persistent directories, overwriting incomplete attempts in place (no timestamped dirs), e.g.:
  `data/dji9443_panelwake/axial_rows072/`, `data/dji9443_panelwake/hover_vinf0p25_rows216/`,
  `data/dji9443_panelwake/hover_vinf0_rows216/`.

## HPC Launcher and Deployment — `scripts/dji9443_panel_wake_convergence.slurm.sh`

Single-node, 64-thread, per the `agent_policies/HPC.md` launcher pattern:
```bash
#!/usr/bin/env bash
#SBATCH --job-name=fp-dji9443-panelwake
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=<RIGHT-SIZED, see below>   # e.g. 12:00:00; set from measured timing, not a blanket 48h
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail
THREADS=64
EXPECTED_REPO=/home/rander39/projects/FLOWPanel.jl
[[ "$PWD" == "$EXPECTED_REPO" ]] || { echo "ERROR: submit from $EXPECTED_REPO; current dir is $PWD" >&2; exit 2; }
export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
       BLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MPLBACKEND=Agg
echo "DJI 9443 direct-potential panel-wake convergence"
echo "  repo:$PWD threads:$THREADS host:$(hostname) job:${SLURM_JOB_ID:-none}"
julia --project=. -t "$THREADS" examples/run_dji9443_panel_wake_convergence.jl
```

Launcher rules: no derived `REPO_ROOT`; no `cd` in the batch script; no `srun` (single threaded Julia
process); assume the cluster's default account/QOS and Julia-on-PATH; submit from exactly
`/home/rander39/projects/FLOWPanel.jl`; logs at repo top level (Slurm opens them before the script runs).

**Right-size the request for scheduling priority (per user).** Ask for **`--mem=32G`** (body-only LU + a
bounded wake row buffer → not memory-heavy), and set **`--time` from measured timing, not a blanket 48 h** —
smaller, more accurate requests get higher queue priority. Known timing data point: the local smoke took
204 s for 12 steps (NT=12, 3 wake rows, 4 threads) — a lower bound only, since production steps carry NT=36
and up to 432 active rows; re-estimate from the first axial rung's reported wall time on the cluster.
Procedure: from the local smoke plus the first
hover-ladder rung, get per-step wall time and extrapolate the full adaptive run (sum over the axial-light
ladder + the 0.25 ladder + the zero-flow confirmation, using each rung's derived `NREVS`); set `--time` to
that estimate + ~25% margin (round to a clean value like `12:00:00`). Under-estimating is safe because the
adaptive driver checkpoint-resumes — the identical resubmission continues from the first incomplete case.
Re-tune `--time` (and `--mem` only if a run ever pressures 32G) after the first real HPC run reports per-case
times.

Deployment:
1. Compare local implementation files with the cluster checkout; abort if any target file has unrelated
   remote modifications.
2. Copy only new/modified study files into `/home/rander39/projects/FLOWPanel.jl` — never a deleting or
   whole-tree sync.
3. Confirm the 40×40 mesh (`examples/data/dji9443_new_40_40.msh`) and the CCBlade checkpoint CSVs exist on the
   cluster.
4. From the cluster checkout run: Julia parse/include checks; the focused formulation test;
   `bash -n scripts/dji9443_panel_wake_convergence.slurm.sh`; and a one-step setup-only smoke if login-node
   policy permits (else rely on the completed local smoke).
5. **Record and then submit** (session-scoped authorization):
   ```
   ssh orc 'cd /home/rander39/projects/FLOWPanel.jl && sbatch scripts/dji9443_panel_wake_convergence.slurm.sh'
   ```
   If the allocation ends before completion, resubmit the identical command; the adaptive driver resumes from
   validated case checkpoints.

## Final Deliverables

- This plan (`plans/20260721_panel_wake.md`).
- Production `DirectWakePotential` formulation + tests.
- Reusable DJI panel-wake driver (`examples/rotor_panel_wake_study.jl`).
- Adaptive convergence driver (`examples/run_dji9443_panel_wake_convergence.jl`).
- 64-thread Slurm launcher pinned to `/home/rander39/projects/FLOWPanel.jl`, `--mem=32G`, right-sized `--time`.
- Local smoke-test report.
- Axial + hover convergence reports with CT/loading gates.
- CCBlade comparison plots/tables.
- Recorded `sbatch` command + submission confirmation.
- Final conclusion reviewed against methods and raw results, distinguishing: solver/formulation correctness;
  wake-row convergence; mesh independence (cited from prior steady study); stability under terminal freestream;
  physical agreement with CCBlade and hover references.

## Verification (of this study, end-to-end)

1. `julia --project=. test/formulation_test.jl` passes locally (formulation correctness + invalid-config
   rejections + Task-3 oracle equivalence).
2. The smoke command above completes with finite fields, VTK+metadata, and an agreeing direct/manual potential
   diagnostic.
3. On the cluster: parse/include checks + focused test + `bash -n` on the launcher pass; the 40×40 mesh and
   CCBlade CSVs are present.
4. After the batch job: each stage has a completion marker; axial CT converges and relates to CCBlade
   ncrit 4/9; hover CT converges in wake-row length and is compared (not equated) to the reference ladder
   0.0505 → 0.072.
