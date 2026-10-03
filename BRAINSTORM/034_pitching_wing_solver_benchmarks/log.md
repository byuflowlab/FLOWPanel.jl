# 034 log — session narrative (newest first)

## 2026-10-02 — Phase 0 session: harness adaptation + availability smoke

- Deliverable 1: `solver_factory` kwarg threaded to both hardcoded Backslash
  sites in `examples/pitching_wing.jl`; default path verified bit-identical via
  the full example test suite (1450/1450 pass).
- Deliverable 2: provisional 4-rung ladder defined (1920/6688/14336/30168
  cells), proportional to the example default knobs; table in ledger.
- Deliverable 3 (headline): LHS-rebuild question ANSWERED BY MEASUREMENT —
  Backslash G is byte-identical across all steps (no rebuild, no refactor), and
  the operator re-assembled from the rotated end-of-run geometry matches t=0 to
  2.9e-14 rel-Frobenius: factor-once LU is exact under the pitching maneuver.
  Counterfactual rebuild cost on R1: ~1.30 s/step (1 thread). Driver
  `benchmark/p034_phase0_lhs_rebuild.jl`.
- Deliverable 4: Dirichlet CONFIRMED (DBC=true type parameter, both bodies);
  `bc_error!` applies as-is, no decision_rules amendment.
- Deliverable 5: 4-arm availability smoke on R1 (`benchmark/p034_avail_smoke.jl`,
  ~10 steps each, per-step certified bc_error!). Results in ledger +
  `data/phase0_avail_smoke/`.
- Environment gotcha (recorded in ledger): local macOS OpenBLAS 0.3.31 cannot
  satisfy the strict single-mode BLAS pin (getter reports 8 after any GEMM, no
  observable 1→4 scaling); smokes launch with BENCH_BLAS_THREADS=8, banner
  records the truth, `common.jl` untouched.
- Phase record: `phase_00_harness.md`. Clear-context subagent review: PASS,
  0 blockers, 4 NOTEs (per-step-denominator doc + "certified" wording fixed in
  place; fgs-seed-cannot-certify-1e-6 and threaded-repack-inert-at-t1 carried
  to Phase 1). Reviewer independently re-ran the suite (1450/1450) and
  re-derived ladder counts. Reported to Ryan; Phase 1 NOT started (gated).

## 2026-10-02 — item created (STAGED)

- Ryan directive: new BRAINSTORM item replicating 021 on the pitching-wing case.
- Session rulings via Q&A: scope = core cost + warmstart (no multi-body, no full
  7-rung tuned ladder); wake pinned `:panel`.
- Exploration established: solver hardcoded `pnl.Backslash` at
  `examples/pitching_wing.jl:1008` (unsteady) and `:867` (static polar), no `solver`
  kwarg on `prepare_pitching_wing`; defaults 6688 cells × ~495 steps; smoke config in
  `test/runtests_example_pitching_wing.jl:28-47`.
- Control doc, decision_rules.md (adapted from 021), and this directory created.
  No code changes, no runs. Next = Phase 0; handoff prompt =
  `phase0_reset_prompt_20261002.md` (this dir).
