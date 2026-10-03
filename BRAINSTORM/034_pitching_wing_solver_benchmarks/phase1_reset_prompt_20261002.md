# 034 Phase 1 reset prompt — written 2026-10-02 at Phase 0 exit

You are continuing BRAINSTORM item 034 (pitching-wing solver benchmarks).
**Phase 0 is COMPLETE and reviewed** (clear-context PASS, 0 blockers,
2026-10-02). Everything is **uncommitted** on branch `fastmultipole` unless
Ryan has committed since — check `git status` / `git log` first and do not
assume either way.

## Read first (in this order, nothing else up front)

1. `BRAINSTORM/034_pitching_wing_solver_benchmarks.md` — control doc: RESET
   BRIEF, standing rulings 1–11 (binding), Phase gates table.
2. `BRAINSTORM/034_pitching_wing_solver_benchmarks/decision_rules.md` — binding
   metrics/protocols (certified BC, cold protocol, min-of-k, iteration counting,
   CSV schema, threading modes).
3. `BRAINSTORM/034_pitching_wing_solver_benchmarks/phase_00_harness.md` and the
   `## 2026-10-02` section of `.../ledger.md` — what Phase 0 built and measured.
4. `agent_policies/WORKFLOW.md` + `agent_policies/TESTING.md` before touching
   repo code; `agent_policies/HPC.md` before any HPC work.
Do NOT read 021's item files inline; use `brainstorm-scout` if 021 context is
needed beyond what decision_rules.md carries.

## State you inherit (Phase 0 outputs)

- `solver_factory` kwarg in `examples/pitching_wing.jl` (default `pnl.Backslash`,
  bit-identical; suite 1450/1450 green). Inject arms via
  `prepare_pitching_wing(; solver_factory=body -> ...)`.
- Provisional ladder (freezes at END of Phase 1):
  R1 (7/89/5)=1920, R2 (13/161/9)=6688, R3 (19/233/13)=14336,
  R4 (27/337/19)=30168 cells.
- LHS question ANSWERED BY MEASUREMENT: G is never rebuilt/refactored during the
  march; operator exactly rigid-motion-invariant (2.9e-14 rel-Frobenius).
  Cite `data/phase0_lhs_rebuild/summary.csv` for the headline; do not re-derive.
- Dirichlet CONFIRMED — `bc_error!` (`benchmark/common.jl`) applies as-is.
- Availability 4/4 PASS on R1 (`data/phase0_avail_smoke/`). Drivers to build on:
  `benchmark/p034_avail_smoke.jl` (arm factory + SmokeFormulation per-step
  bc_error! wrapper honoring the velocity-snapshot entry contract),
  `benchmark/p034_phase0_lhs_rebuild.jl`.

## Phase 1 contract (gate row in the control doc)

All four arms (`backslash`, `krylov_gmres`, `krylov_ilu_nfcache`, `fgs`) agree
on **certified BC error ≤ 1e-6** and on force histories (CL hysteresis,
reported never thresholded) on **2 rungs**; then FREEZE the ladder and per-rung
solver/FMM settings (record in the ledger as the Phase 1 freeze). Expect
per-rung FMM-knob and FGS τ tuning (021 values are seeds, recorded in `notes`).

## Flags carried from Phase 0 (do not rediscover)

1. **fgs seed knobs cannot certify 1e-6**: p=4/mac=0.5/leaf=50/inner=2 carries
   an apply-accuracy BC floor ~3–6e-6 rel that GROWS with wake rows
   (steps.csv, 3.1e-6 → 6.2e-6 over 10 steps). Phase 1 must retune (021
   ruling: use the τ=1e-6-tuned config for the solver role).
2. **Local macOS BLAS caveat**: the single-mode pin is unsatisfiable locally
   (getter reports 8 after any GEMM, no observable scaling). Local smokes
   launch with `THREADING_MODE=single EXPECT_JULIA_THREADS=1
   BENCH_BLAS_THREADS=8 julia --project -t 1 ...`. `benchmark/common.jl` is
   UNTOUCHED — keep it that way; HPC/Linux keeps the strict pin.
3. 033's threaded repack in `krylov_ilu_nfcache`/fgs is inert at `-t 1` — it
   has NOT been exercised in this item yet; multi-thread mode first exercises
   it (on HPC).
4. Per-step timing convention: the maneuver callback fires once per `t_range`
   point including t=0 (nsteps instrumented, nsteps−1 timesteps marched).

## Hard constraints (unchanged)

- Local ≤4 threads, smoke only — nothing local is publishable. Any run >20 min
  goes to HPC under campaign-worktree rules (annotated tag
  `campaign/p034-<slug>-YYYYMMDD`, pins in the ledger, outputs to the
  consolidated data root). Read `agent_policies/HPC.md` and memory on the orc
  precompile race before submitting.
- Wake pinned `:panel`; formulation pinned `VelocityThroughSources`; no
  multi-body; no particle wake. Cold = zero-initial-guess (arms may share a
  process); threading modes never mixed in one comparison.
- Judge from harness-written CSVs with a banner beside each; 021 runs.csv
  schema for static benchmarks (`validate_runs_csv`).
- Ledger = data of record (dated, append-only); log.md newest-first narrative;
  create `phase_01_consistency.md` for the phase record; update the control
  doc's `## Current status` + RESET BRIEF in place.
- Phase gate = clear-context subagent review before the box is ticked.
- **Ryan-gated, offer don't do:** commits, pushes, notebook entries, INDEX
  checkbox ticks, and the START of Phase 1 itself if he has not yet approved
  Phase 0 — if `git log` shows Phase 0 still uncommitted and no approval is in
  the control doc's approvals/decision log, ASK before running Phase 1 work.

## Also pending from Phase 0 (offer to Ryan, don't do unprompted)

- Commit of Phase 0 (suggested split: code/harness vs BRAINSTORM records).
- Notebook entry for Phase 0 (LHS measurement + smoke tables; ask how much).
- Phase 0 INDEX/approval ticks.
