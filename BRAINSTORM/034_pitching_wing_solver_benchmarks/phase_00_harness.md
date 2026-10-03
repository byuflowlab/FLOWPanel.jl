# 034 Phase 0 — harness adaptation + availability smoke

**Date:** 2026-10-02 · **Status:** deliverables complete, pending clear-context
review + Ryan approval · **Environment:** local Mac, smoke only (≤4 threads),
NOTHING here is publishable. Data of record: `ledger.md` §2026-10-02; CSVs under
`data/phase0_lhs_rebuild/` and `data/phase0_avail_smoke/` (banner beside each).

## Deliverable 1 — solver injection

`solver_factory` kwarg added to `prepare_pitching_wing` and
`run_pitching_wing_static_polar` (`examples/pitching_wing.jl`), default
`pnl.Backslash`:

- unsteady site (was hardcoded `pnl.Backslash(wing)` at `:1008`): now
  `solver = solver_factory(wing)`; consumed unchanged as `body_solvers=(solver,)`
  by `simulate!`/`simulate_warmstart!`;
- static-polar site (was `body_solvers=pnl.Backslash(body)` at `:867`): now
  `solver_factory(body)` per alpha (factory, because the polar builds a fresh
  body per alpha);
- `_run_pitching_wing` forwards `sim.solver_factory` to the static polar.

Default path bit-identical by construction (`pnl.Backslash` as the factory IS
the old call). Verified: `test/runtests_example_pitching_wing.jl` green,
1450/1450 pass (includes the smoke config `n_cycles=0.01, n_span=1,
n_airfoil=21, n_endcap=5`, DirectBackend).

## Deliverable 2 — provisional mesh ladder

4 rungs, 1920 / 6688 / 14336 / 30168 cells, knobs proportional to the example
default; table + formula verification in `ledger.md`. PROVISIONAL — freezes at
the end of Phase 1.

## Deliverable 3 — LHS-rebuild measurement (headline question)

**Answer: the LHS is NEVER rebuilt or refactored during the pitching
time-march; reuse is exact to roundoff.** Measured by
`benchmark/p034_phase0_lhs_rebuild.jl` on R1 (1920 cells, 20 steps):

- `solver.G` byte-identical to its construction state at every step and after
  the final kinematics (max deviation 0.0); `Glu` never reassigned;
- rotation-invariance measured directly: raw G re-assembled from the end-of-run
  rotated geometry matches the t=0 raw G to 1.2e-13 max-abs / 2.9e-14
  rel-Frobenius — so factor-once dense LU stays EXACT under the maneuver;
- counterfactual per-step rebuild on R1 at 1 thread: 1.259 s assembly +
  0.039 s LU (min-of-5) — the cost simulate! avoids each step.

Mechanism (file:line evidence in `ledger.md`): construction-time assembly+LU,
per-step RHS-only `_solve!`, `update_G` never passed; per-step
`transform_body_solvers!` mirrors the rigid delta into persistent FMM state
(no-op for Backslash, `transform_plan!` for persistent-plan Krylov,
tree transform for FGS). The newest wake row couples through the RHS only.
Consequence for the campaign: the Phase 3 crossover question is genuinely
"factor-once LU amortization vs iterative per-step cost" — there is no hidden
per-step refactorization cost on the dense arm.

## Deliverable 4 — Dirichlet confirmation

CONFIRMED (type parameter DBC=true on both static and unsteady bodies; existing
test asserts `has_dirichlet_bc`; `bc_error!` re-asserts at runtime). The
`decision_rules.md` Dirichlet-only caveat is discharged; no amendment needed.
Note the formulation nuance: under pinned `VelocityThroughSources` the wake
enters the BC through the control-point velocity (σ side), so the body-only
`bc_error!` pass measures exactly the equation the solver solved — same
structure 021 Phase 3 relied on (`rotor_hover_solver_unsteady.jl:539-601`).

## Deliverable 5 — availability smoke (4 arms, rung R1)

Driver `benchmark/p034_avail_smoke.jl` (reuses `benchmark/common.jl`: banner,
`bc_error!`, `solver_state_bytes`); arms `backslash`, `krylov_gmres`,
`krylov_ilu_nfcache` (021 production-winner construction), `fgs`
(dagteam+backoff defaults; p=4/mac=0.5/leaf=50/inner=2 and
tol_abs=1e-6·rms_b(t=0) recorded as SEEDS). One process, fresh body per arm
(cold = zero-initial-guess per Ryan 2026-09-23). 10 steps each, per-step
certified `bc_error!` against the fixed t=0 RHS scale (normalization constant
only; thresholds are Phase 1's job). **Result: 4/4 arms completed, every BC
pass certified (certified = the FMM measurement pass met its requested error
bound, `error_success` — NOT "BC ≤ target", which Phase 1 owns); CL identity
agrees to ~5 digits across arms.** The fgs seed
knobs carry an apply-accuracy BC floor (~3–6e-6 rel, growing with wake rows) —
the known 021 floor phenomenon, owned by Phase 1 tuning. Table in `ledger.md`;
CSVs in `data/phase0_avail_smoke/{summary,steps}.csv`.

## Environment caveat

Local macOS OpenBLAS 0.3.31 defeats `assert_and_banner`'s single-mode BLAS pin
(getter reports 8 after any GEMM; timing shows no 1→4 scaling). Phase 0 smokes
launch with `BENCH_BLAS_THREADS=8` so the assert passes and the banner records
the truth. `benchmark/common.jl` untouched; HPC runs keep the strict pin.

## Files touched

- `examples/pitching_wing.jl` — `solver_factory` kwarg (deliverable 1)
- `benchmark/p034_phase0_lhs_rebuild.jl` — NEW (deliverable 3)
- `benchmark/p034_avail_smoke.jl` — NEW (deliverable 5)
- `BRAINSTORM/034_pitching_wing_solver_benchmarks/{ledger,log}.md`, this file,
  control doc `## Current status` + RESET BRIEF

## Exit state

All five deliverables complete; clear-context subagent review: [see log.md].
Commits/pushes/notebook/INDEX ticks Ryan-gated — none performed.
