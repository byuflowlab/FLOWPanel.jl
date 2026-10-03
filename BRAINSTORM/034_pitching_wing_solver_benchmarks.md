# 034 — Pitching-wing solver benchmarks (021 replication on a second case)

**Opened:** 2026-10-02 (Ryan directive: replicate 021 on a different case — the pitching wing under `examples/`)
**Status:** STAGED — plan only, nothing run.
**Item-level approvals:** Technical [ ] · Clear-context [ ] · User-after-AI-discussion [ ]
**Session narrative:** `034_pitching_wing_solver_benchmarks/log.md` · **Data of record:** `.../ledger.md` · **Binding metrics:** `.../decision_rules.md`

## RESET BRIEF

- **What this is:** replicate BRAINSTORM 021's solver-benchmark methodology (Backslash,
  GMRES, GMRES+ILU-nfcache, FGS dagteam+backoff) on the unsteady pitching wing
  (`examples/pitching_wing.jl`), measuring cold setup-vs-per-step cost and warmstart
  behavior over a mesh ladder, ending in a cumulative time-marching crossover statement.
- **Scope rulings (Ryan 2026-10-02):** core cost + warmstart only — NO multi-body phase,
  NO full 7-rung per-rung-FMM-tuned ladder; wake model pinned to `:panel` (`PanelWake`).
- **Status:** Phase 0 deliverables COMPLETE 2026-10-02 (local smoke, uncommitted):
  `solver_factory` kwarg in `examples/pitching_wing.jl` (default bit-identical,
  suite 1450/1450), provisional ladder 1920/6688/14336/30168, LHS question
  ANSWERED (reused, never rebuilt; rotation-invariant to 2.9e-14), Dirichlet
  CONFIRMED, 4-arm availability smoke on R1. Records:
  `034_.../phase_00_harness.md`, ledger §2026-10-02, `data/phase0_*/`.
  Clear-context review PASS 2026-10-02 (0 blockers). Ryan APPROVED Phase 0 +
  Phase 1 go-ahead 2026-10-02 (committed; notebook entry still pending).
  Next-agent entry prompt =
  `034_pitching_wing_solver_benchmarks/phase1_reset_prompt_20261002.md`.
- **Headline scientific angle:** MEASURED (Phase 0): `simulate!` never
  rebuilds/refactors G under the pitching maneuver — Backslash G byte-identical
  across steps, operator exactly rigid-motion-invariant (2.9e-14 rel-Frobenius
  after accumulated pitch); persistent Krylov plans / FGS trees are transformed,
  not rebuilt (`transform_body_solvers!`). The Phase 3 crossover is cleanly
  "factor-once LU amortization vs iterative per-step cost".
- **Key dependencies:** `benchmark/common.jl` (`bc_error!` — Dirichlet-only; Phase 0
  CONFIRMED the wing body qualifies), 033's threaded-repack ILU + dagteam champion
  configs (unexercised at `-t 1` so far), `examples/pitching_wing.jl` (solver
  injection DONE: `solver_factory` kwarg, default bit-identical).
- **Standing conventions inherited from 021:** cold = zero-initial-guess (Ryan
  2026-09-23); threading modes never mixed; judge from CSVs not stdout; campaign runs
  from git-tag-pinned worktrees on HPC (`campaign/p034-<slug>-YYYYMMDD`); ≤4 threads
  local = smoke only.
- **Relation to mission:** solver-infrastructure item like 021/033 — tangential to the
  rotor-hover CT mission statement.

## Current status

Phase 0 deliverables complete 2026-10-02 (local smoke; everything uncommitted).
All five Phase 0 gate items done: solver injection (suite green), provisional
ladder (ledger table), LHS-rebuild measured (reuse, exact to roundoff),
Dirichlet confirmed, 4/4 arms smoked on R1 with certified per-step BC passes.
Clear-context review PASS 2026-10-02 (0 blockers, 4 NOTEs — 2 fixed in place,
2 carried to Phase 1: fgs seed knobs cannot certify 1e-6 without retuning;
threaded repack unexercised at -t 1). Ryan APPROVED Phase 0 and gave the
Phase 1 go-ahead 2026-10-02; Phase 0 committed (code/harness + records split).
Phase 1 (consistency/calibration on 2 rungs) IN PROGRESS. Notebook entry still
pending. Next-agent entry prompt =
`034_pitching_wing_solver_benchmarks/phase1_reset_prompt_20261002.md`.
Local-env caveat: macOS BLAS pin unsatisfiable — smokes use
BENCH_BLAS_THREADS=8 (see ledger).

## Objective and scope

Replicate 021's publishable solver-benchmark methodology on a structurally different
case: a rigid NACA 0015 wing pitching sinusoidally in an unsteady time-march with a
growing panel wake. 021 answered "which solver wins on a rotor-hover solve"; 034 asks
how those answers transfer when

1. the body undergoes rigid motion (LHS plausibly constant or cheaply transformable —
   measure, don't assume; the newest wake row couples to the unknowns, so "constant" is
   a hypothesis, not a given);
2. the RHS grows each step with accumulating `PanelWake` rows;
3. the natural production pattern is a long time-march (~495 steps over 3 pitch
   cycles) where setup amortization and warmstart quality vary over the cycle.

Headline deliverable: a cumulative wall-clock vs step-count crossover table/figure over
a full 3-cycle run per champion arm, analogous to 021's "no break-even horizon" result,
explicitly locating where dense factor-once LU sits relative to the iterative champions.

**Out of scope:** multi-body (021 Phase 4), particle wake (`:particle` — only as a
later robustness check if Ryan approves), full 7-rung ladder with per-rung FMM-knob
descent campaigns (a 3–4-rung ladder with lighter per-rung tuning instead), non-default
formulations (Backslash-only, so pinned off per fairness).

## Case definition

- Case: `examples/pitching_wing.jl` — NACA 0015 rectangular wing, sinusoidal pitch
  about quarter chord, α(t) = 3.94° ± 1.99°, f = 4.01 Hz; static α-polar pass then
  unsteady run. Maneuver enters as frame angular rate (`pitching_maneuver!`, `:971`).
- Mesh: in-script `pitching_wing_mesh` (`:67`); knobs `n_span`, `n_airfoil`,
  `n_endcap`; cells = `2*n_sec*n_span + 4*(n_chord-2)*(n_endcap-1)` (`:93-95`);
  defaults (13/161/9) → 6688 cells. Ladder (Phase 0) spans ~2k → ~30k panels over
  3–4 rungs via these knobs; frozen after Phase 1.
- Time marching: `dt = c/U * c_per_dt`, `c_per_dt = 0.5`, `n_cycles = 3` → ~495 steps
  (~165/cycle). Warmstart/restart machinery exists (`simulate_warmstart!`, `:1159`).
- Wake: **pinned `:panel`** → `pnl.PanelWake`, ConstantDoublet rows,
  `das_chord_fraction = 0.05`, `wake_length_spans = 2.0` (`:954-968`, `:995-999`).
- Solver today: hardcoded `pnl.Backslash` — unsteady `:1008` (consumed as
  `body_solvers=(solver,)` at `:1160`/`:1166`), static polar `:867`.
- Formulation: pinned default `VelocityThroughSources` (non-default formulations are
  Backslash-only). Phase 0 confirms the body is Dirichlet so `bc_error!` applies
  as-is; a Neumann body would need `true_residual!` and a re-derived tolerance.
- Monitors (kept for identity reporting, not benchmark timing):
  `PressureBernoulli(unsteady=true)`, `ForceMonitor` + `WingNormalization`,
  `SpanwiseLoadingMonitor`, `SectionLiftHistoryMonitor` (`:1011-1028`). CL hysteresis
  loop is the physical identity signal (034's analog of 021's CT column) — reported,
  never thresholded.
- Smoke config: `test/runtests_example_pitching_wing.jl:28-47`
  (`n_cycles=0.01, n_span=1, n_airfoil=21, n_endcap=5`, DirectBackend).

## Solver matrix

| arm | constructor | notes |
| --- | --- | --- |
| `backslash` | `pnl.Backslash` (`src/FLOWPanel_solver.jl:442`) | dense LU; the incumbent; `niter = -1` by convention |
| `krylov_gmres` | `pnl.KrylovSolver` (`:1013`) | unpreconditioned GMRES baseline |
| `krylov_ilu_nfcache` | KrylovSolver + near-field ILU (direct-interaction-list pattern) | 021 production iterative winner; 033 threaded repack + `pattern_time` bracket |
| `fgs` | `pnl.FGSSolver` (`:1547`) | 021/033 champion config: dagteam+backoff, f32-full, NUMA-interleave; placement/j re-tuned for this case |

Per-rung FMM knobs and FGS τ-config are re-tuned for this case (seeded from 021 values
as hypotheses, never carried as measurements — 021 decision-rules seeding principle).
Usage patterns to copy: `examples/sweptwing_solverbenchmark.jl:114`,
`examples/suddenly_started_wing.jl:513`.

## Standing rulings (binding on every phase)

1. **Scope (Ryan 2026-10-02):** core cost + warmstart; no multi-body; no full ladder.
2. **Wake pin (Ryan 2026-10-02):** `:panel` for every certified comparison.
3. **Cold = zero-initial-guess solve** (Ryan 2026-09-23, carried from 021/033), not
   fresh-process; harness resets `body.strength` via untimed `setup!` before every
   timed/history/alloc rep; min-of-k (adaptive k per 021 timing protocol) after warmup.
4. **Certification:** static-RHS benchmarks (Phases 1–2) use certified BC error
   (`bc_error!`, `benchmark/common.jl`) ≤ 1e-6 with `error_success=true`; unsteady
   time-marching (Phase 3) uses the per-step BC-satisfaction guard — each arm checked
   against the tolerance IT promised (021 ruling 2026-08-25), never the frozen static
   metric. Never a solver's internal residual.
5. **Separation of costs:** setup (assembly / factorization / tree / preconditioner)
   vs per-step (RHS + solve + warmstart projection) always timed separately;
   `t_step_net` convention for any timestep-share claim.
6. **Threading:** two modes, never mixed in one comparison — (`-t 1` + explicit
   `BLAS.set_num_threads(1)`) and full-thread with BLAS count recorded; banner asserts
   and logs both. Local runs ≤ 4 threads = dev/smoke only, never published.
7. **Fairness pins:** identical mesh per comparison; formulation
   `VelocityThroughSources`; identical kernel/offset; FMM knobs tuned per rung per
   solver family and recorded in the CSV.
8. **Judge from harness-written CSVs, never stdout**; banner written beside each CSV;
   provenance columns per 021 schema (`commit, fm_commit, vpm_commit, julia_version,
   filament_reg, worktrees, solver_settings, backend_settings`).
9. **Iteration counting:** `niter` units are solver-specific (Krylov iterations vs FGS
   sweeps vs FGMRES outer) and never cross-solver comparable; wall-time-to-target is
   the only cross-solver currency; warmstart headline reads `step_niter_first`.
10. **Campaign reproducibility:** published numbers from git-tag-pinned worktrees
    (`campaign/p034-<slug>-YYYYMMDD`) on exclusive HPC nodes; outputs to the
    consolidated data root; md5/clean-tree verification before submission.
11. **Review gates:** each phase passes a clear-context subagent review before its
    box is ticked; Ryan's approvals are separate and may lag.

## Phase gates

| phase | deliverable | file | status | review |
| --- | --- | --- | --- | --- |
| 0 — Harness adaptation + availability | `solver` kwarg threaded to `pitching_wing.jl:1008`/`:867` (default-path bit-identical, existing smoke test passes); 3–4-rung ladder defined with computed cell counts (~2k→~30k); **measured answer to "is G rebuilt/refactored per step under pitching?"** (instrument `simulate!`; account for the newest-wake-row coupling); Dirichlet confirmation for `bc_error!`; all four arms smoke on the smallest rung | `034_.../phase_00_harness.md` | DELIVERABLES COMPLETE 2026-10-02 | PASS 2026-10-02 (clear-context, 0 blockers; Ryan approval separate) |
| 1 — Consistency/calibration | all arms agree on certified BC error and force histories (CL hysteresis) on 2 rungs; ladder + per-rung solver settings FROZEN | `034_.../phase_01_consistency.md` | NOT STARTED | [ ] |
| 2 — Cold cost | setup vs per-step cold cost across the ladder, both threading modes, min-of-k; fitted t ∝ N^p exponents; memory per 021 ruling 8 | `034_.../phase_02_cold_cost.md` | NOT STARTED | [ ] |
| 3 — Warmstart + crossover | warmstart matrix (cold / prev / extrap) at Phase-2 frozen settings; full 3-cycle time-marched runs per arm; cumulative wall-clock vs step crossover incl. factor-once LU; warmstart quality vs pitch phase | `034_.../phase_03_warmstart_crossover.md` | NOT STARTED | [ ] |

## Acceptance criteria (item level)

- Every published static-solve point carries certified BC ≤ 1e-6 (`error_success`);
  every published unsteady step passes the per-step BC-satisfaction guard.
- Cold setup/per-step matrices complete at ≥ 2 rungs in both threading modes.
- Cumulative-cost crossover table over ≥ 495 steps (3 cycles) for all four arms with
  full provenance pins, and an explicit statement of where (or whether) each iterative
  arm overtakes factor-once dense LU.
- The Phase-0 LHS-rebuild question answered with a measurement, cited in the headline.
- All published numbers from tagged-worktree HPC runs.

## Relationship to other items

- **021** (closed at promotion): methodology source — rulings, `bc_error!`, champion
  configs, CSV schema, timing protocol; also `rigid_motion_tree_reuse_item.md`
  (`transform_tree!` under rigid motion) which anticipates 034's LHS-reuse angle.
- **033** (live): FGS setup-cost threading (dagteam repack) and ILU attribution —
  034 inherits whichever defaults Ryan promotes there.
- **Pitching-wing prior art:** `examples/pitching_wing_convergence.jl` (env-knob
  pattern), `code_audit/scripts/task0–task2b` (frozen-state reuse of
  `build_pitching_wing_body`), `test/formulation_test.jl:30-43`.

## Decision log

- **2026-10-02 (Ryan, via session Q&A):** scope = core cost + warmstart (no
  multi-body, no full 7-rung tuned ladder); wake model pinned `:panel`.
- **2026-10-02:** item created; phases 0–3 defined; solver matrix fixed at
  backslash / krylov_gmres / krylov_ilu_nfcache / fgs-champion.
- **2026-10-02 (Ryan):** Phase 0 APPROVED; commit authorized (code/harness +
  BRAINSTORM records split); Phase 1 go-ahead GIVEN. Notebook entry still
  pending (Ryan to specify depth).

## Logging provision

Data of record (tables, certified CSV pointers, frozen settings) accumulate in
`034_pitching_wing_solver_benchmarks/ledger.md` in dated sections. Session-by-session
narrative goes to `.../log.md`, newest first. Binding metric/threshold definitions in
`.../decision_rules.md`, changeable only via this file's Decision log. Do not append
session narrative to this file; `## Current status` and `## RESET BRIEF` are edited in
place, never appended to.
