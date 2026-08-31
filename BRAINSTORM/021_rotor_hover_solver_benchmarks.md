# 021 — Rotor-hover solver benchmark campaign (availability → consistency → cost → warmstart → multi-body)

**Opened:** 2026-08-06 (Ryan). **Status:** ACTIVE — Phase 0 technically complete (W1–W6), gate review in progress; Phases 1–4 not started.

**Item-level approvals:** Technical [ ]; clear-context [ ]; user [ ]

**Session history:** [`021_rotor_hover_solver_benchmarks/log.md`](021_rotor_hover_solver_benchmarks/log.md) (newest first)

## Current status

> **STALE BELOW — last revised 2026-08-13.** Fresh agent: start from
> [`phase_20_context_reset_prompt.md`](021_rotor_hover_solver_benchmarks/phase_20_context_reset_prompt.md)
> (2026-08-25), which chains to `phase_18` for cluster state. As of 2026-08-22
> the campaign is in Phase 2, not "Phase 1 starting". Read
> `021_rotor_hover_solver_benchmarks/log.md` (newest first) for true state; the
> narrative below is kept for Phase-0 provenance only. Live summary: Phase 1
> ladder frozen and tuned through R6 (R7 tune + R6 table failed, causes
> diagnosed 2026-08-22); Phase 2 active since 08-17 with the near-field cache,
> rigid-motion tree reuse and the FGS unsteady staleness fix landed; the
> Phase-2b HPC campaign was 0/8 and awaits relaunch. Note also that
> `fgs_wake_plateau_handoff_prompt.md` is SUPERSEDED — the plateau it stages a
> hunt for was root-caused the same day as H3 (lifecycle ordering), not
> H1/H2/H4.

**Phase 0 TECHNICALLY COMPLETE 2026-08-13 (W1–W6).** W6 (sparse near-field ILU) was
implemented directly by Ryan (Barba direct-list pattern + `ILUZero.jl`; the staged
`phase_00_ilu_plan.md` is superseded), then session-reviewed (correct; ctor deduplicated,
assembly threaded, test-gated) and smoke-run: **acceptance PASS both threading modes** —
ILU-GMRES 290→15 iterations, 223.4→11.5 s single / 55.2→2.86 s multi vs plain GMRES,
setup ~1 s, break-even <0.02 solves; displaces FGMRES+FGS as best iterative smoke config.
Evidence: phase_00 Log 2026-08-13 tables + `benchmark/results/smoke/{single,multi}/`
(CSVs of record = the 08-13 metadata-remediation rerun; earlier numbers persist in the
phase_00 md tables). The first clear-context gate review verified all substance but
DECLINED on evidence-record consistency (CSV `warmstart` mislabel + stale doc text) —
all remediated same day, smoke rerun both modes (concurrent, multi @4 threads, per
Ryan; caveats in the tables); the **second clear-context review APPROVED the Phase 0
gate 2026-08-13** (all remediations verified, tables ≡ CSVs cell-for-cell, tests green;
Ryan's own re-approval checkbox pending his next review). **Next action: Phase 1 in
progress per Ryan 2026-08-13** (`phase_01_consistency.md`: frozen mesh ladder, FMM
residual floor, matched-residual calibration). Standing context: rulings 11 (Krylov honesty;
theory in `theory_nearfield_ilu.md`) and 12 (warmstart = headline deliverable); FGS
thread-degradation (reconfirmed: multi 500 iters this run vs 267 prior — order-dependent)
and FMM-tree-per-apply findings pre-registered for Phase 2.

## Objective and scope

Compare the cost and convergence behavior of FLOWPanel's solvers on the DJI9443 rotor-hover case
across a ladder of increasing panel counts, separating one-time setup cost from per-timestep cost,
and quantify what warmstarting and multi-body coupling do to each. **Publication target**: this
campaign is one story of a paper; a companion static case (wing with rigid wake) is handled
elsewhere and is out of scope here.

- Drivers live in `examples/` and `benchmark/` (new drivers go in `benchmark/`, results in
  `benchmark/results/` and `data/`).
- `src/` changes are allowed only where a phase explicitly authorizes them (Phase 0 unblocking +
  instrumentation, Phase 2b optimizations) and must be test-gated.
- Published numbers run on HPC (see standing ruling 5); local runs (≤6 threads) are
  development/smoke only.

## Solver matrix (verified 2026-08-06, `src/FLOWPanel_solver.jl` unless noted)

The paper's comparison roster is rulings 1a–1f below. What exists today:

| Solver | Status | Notes |
| --- | --- | --- |
| `Backslash` (:213) | exists | Per-body. `G` assembly (`_G!` :137) + `lu!` in **constructor** (one-time); per-step = RHS + `ldiv!`. Neumann/Dirichlet chosen by the body's `DBC` type param (dispatch at :688/:710) — there are no separate structs. No options. |
| `BackslashCoupled` (:969) | exists | Monolithic multi-body; block `G` + `lu!` built inside `solve!` under `update_G` kwarg (:1141–1160), cached in mutable `Glu`. `t_build`/`t_solve` `@elapsed` timers exist but the `return` is commented out (:1175). |
| Tuple-of-solvers block Gauss–Seidel (:851) | exists | The "iterative" multi-body path; what `simulate!` uses (`body_solvers::Tuple`). Couples bodies **through velocity only** (:876) — Neumann-correct, Dirichlet-incorrect. Per-body `Backslash` LUs stay frame-invariant for rigid bodies; only cross-influence is recomputed per outer iteration. Converges on `max |Δstrength|`. |
| `KrylovSolver` (:296) | exists | Matrix-free; `method::Symbol` → any Krylov.jl method, only `:gmres` exercised. Options: `itmax=20, atol=1e-6, rtol=1e-6, preconditioner_cell_size` (>0 builds `FastMultipole.JacobiPreconditioner` once at construction). Struct is **immutable** → no per-timestep preconditioner swap without refactor. Fresh workspace allocated every solve (:430); **no `x0` ever passed**. |
| `KrylovCoupled` (:488) | exists | Monolithic matrix-free multi-body; **no preconditioner field at all** (:596). |
| `FGSSolver` (:626) | exists | Wraps `FastMultipole.FastGaussSeidel` (FMM far-field pass + per-leaf dense block sweeps; `FastMultipole/src/solve.jl:797–1004`). Tree + leaf matrices built once at construction, but leaf blocks are **re-factorized every sweep** (`mat \ rhs` — no stored LUs; caching deferred to Phase 2b per 2026-08-11 D3). Holds the **only existing warmstart**: `solution_history_length` / `project_solution` / `project_solution_order` = rolling history + finite-difference (Taylor-style) extrapolation (`project_solution!` :758, `save_solution!` :775); off by default. |
| Krylov + FGS preconditioner | **missing** | No code path passes FGS as `M`/`N`. Needs `:fgmres` (flexible GMRES; in Krylov.jl, unexercised here), a preconditioner-apply wrapper around fixed FGS sweeps, and the `KrylovSolver` mutability refactor. |

**Blockers and gaps (Phase 0's work items):**

1. ~~Kutta closure hard-errors on Krylov/FGS~~ — RESOLVED BY SCOPE (2026-08-11): the error
   (`src/FLOWPanel_kutta.jl:495-503`) only gates the opt-in 015 pressure-continuity closure;
   the rotor pipeline uses the legacy linear Kutta condition, which works with every solver.
   W1 dropped per Ryan; the error stays.
2. Non-default formulations (`GreenReconstruction` etc., `src/FLOWPanel_formulation.jl:806-900`)
   call `solver.Glu` directly → Backslash-only. Hence ruling 2 pins the default formulation.
3. **No convergence-history capture anywhere**: Krylov `stats` is `@show`n and discarded
   (:436, :474); the FGS residual lives in `solver.fgs.residual_vector` but is never read;
   block-GS prints max-Δstrength only under `verbose=true`.
4. No Krylov warmstart: fresh workspace + zero `x0` every solve; the
   `krylov_solve!(ws, A, b, x0)` form is unused. `set_strengths!` (:1083) zeroes `body.strength`
   before solves — any warmstart insertion must come after it.
5. `BenchmarkTools` appears nowhere in the repo; existing timing is ad-hoc `@elapsed`.

**Reusable assets:**

- `examples/rotor_hover.jl` — canonical unsteady case; solver swap point stubbed (commented
  `FGSSolver` block at :128–140, `Backslash` at :141). Mesh hardcoded; 53 `dji9443_2026*` meshes in
  `examples/data/` (`<date>_<nchord>_<nspan>_{capped,uncapped}`) supply the ladder.
- `examples/rotor_hover_convergence.jl` — ENV-parameterized single `steady!`; Phase 1 harness
  skeleton (mind the shedding-from-constructed-cells invariant in CLAUDE.md).
- `examples/sweptwing_solverbenchmark.jl` — the only prior solver benchmark; superseded by this item.
- `examples/two_ducts.jl` — multi-body comparison pattern (BackslashCoupled vs per-body FGS vs
  tuple) with tangency-error reporting.
- `examples/simple_wing.jl` — NACA 0015 rectangular wing (Phase 4 second body, per Ryan).
- Rotor monitors already produce CT three ways (PressureBernoulli / PressureLaplace /
  KuttaJoukowski) — Phase 1 cross-checks.
- `src/FLOWPanel_metadata.jl:153` `_solver_metadata_dict` — extend for run metadata.

## Standing rulings (binding on every phase)

1. **Solver roster** (the paper's comparison set):
   (a) Backslash-coupled (`BackslashCoupled`); (b) Backslash-iterative (tuple block-GS with
   per-body `Backslash`); (c) Krylov GMRES, unpreconditioned; (d) Krylov GMRES + block-Jacobi
   preconditioner; (e) `FGSSolver`; (f) FGMRES with FGS as preconditioner (implemented in
   Phase 0). Single-body phases (0–3) use the single-body forms; the (a)/(b) distinction
   activates in Phase 4.
2. **Fairness pins**: identical mesh, `formulation = VelocityThroughSources` (default; the
   alternatives are Backslash-only), identical kernel/offset settings; backend recorded per
   solver with FMM knobs pinned and recorded. Tolerance settings come from Phase 1 calibration
   and are then frozen.
3. **Residual metric**: RMS and max of the **true residual** $r = G x - b$ (unpreconditioned,
   physical units), computed by one shared evaluator applied identically to every solver's
   output — never a solver's internal metric (FGS reports MSE; Krylov a preconditioned norm).
   Matrix-free solvers get $r$ via one direct influence evaluation at the converged state.
4. **Convergence history**: iterative solvers record per-iteration residual **and a
   per-iteration wall-clock timestamp**; plots are residual vs. time on shared axes. The
   equal-time-per-iteration assumption is *validated, not assumed* (GMRES iteration cost grows
   with basis size; an FGS iteration contains an FMM pass). If per-iteration times are flat
   within ~10%, later plots may use iteration count × mean time.
5. **Timing protocol**: published numbers from a single dedicated HPC node, exclusive
   allocation; min-of-k (k ≥ 5) after warmup for the isolated solve. Record hardware, commit
   hash, Julia version, and thread counts per run. Local runs (≤6 threads) are never published.
6. **Threading (Ryan 2026-08-06)**: every benchmark runs in BOTH modes — (i) single-threaded:
   Julia 1 thread AND `BLAS.set_num_threads(1)` set explicitly (BLAS defaults to multithreaded
   even under `julia -t 1`, which would silently advantage the dense-LU solvers); (ii)
   multithreaded: 64 threads (`julia -t 64`), BLAS thread count recorded and held consistent
   across solvers. Never mix modes within a comparison. The harness asserts
   `Threads.nthreads()` / `BLAS.get_num_threads()` at startup, writes them to the log banner,
   and every CSV row records both.
7. **Cost split**: every benchmark separates one-time setup (matrix assembly, LU/factorization,
   tree build, preconditioner build) from per-timestep cost (RHS assembly, solve/`ldiv!`,
   warmstart projection). Primary per-step metric = isolated solve on frozen RHS; secondarily,
   one full-`simulate!` timestep per config gives the solve's share of a real step
   (Ryan 2026-08-06: both, solve-first).
8. **Memory** is tracked alongside time (matrixful $O(N^2)$ vs matrix-free $O(N)$ storage is a
   publishable axis; report the crossover N if observed).
9. **Mesh ladder**: ≥4 rungs of increasing panel count from `examples/data/dji9443_*`, chosen in
   Phase 1 setup and then **frozen for all phases**. Fit and report empirical scaling exponents
   $t \propto N^p$ per solver per cost component.
10. Judge runs from harness-written CSVs, not stdout (018 lesson). All `src/` changes
    test-gated. Verify knobs from the log banner.
11. **Krylov honesty and accuracy (Ryan 2026-08-12)**:
    (a) every preconditioned Krylov config is routed as a *right* preconditioner (`N=`)
    so the monitored/stopping residual is the true one — left-routed configs stop on a
    preconditioned norm and can report "converged" at a much worse true residual (observed:
    left-Jacobi stopped 20× worse than plain GMRES at the same nominal tolerance);
    (b) production tolerances come from Phase 1's matched-true-residual calibration, never
    from smoke-scale settings (the Phase 0 smoke ran `rtol=1e-4` for local runtime only —
    GMRES's 2.6e-2 rms residual was the requested tolerance, not a solver ceiling);
    (c) the FMM apply sets a floor on the achievable *direct-evaluated* true residual —
    Phase 1 measures that floor per rung (residual of a Backslash solution evaluated through
    the FMM operator) and sets residual targets above it, tuning `expansion_order` if the
    floor sits above the target;
    (d) preconditioner exploration beyond block Jacobi is a Phase 2b deliverable, with
    **sparse near-field ILU as the priority candidate** — design basis in
    [`theory_nearfield_ilu.md`](021_rotor_hover_solver_benchmarks/theory_nearfield_ilu.md).
12. **Warmstart benefit is a headline deliverable, not a side study (Ryan 2026-08-12)**:
    for every iterative solver the campaign reports the measured per-step savings from
    warmstarting (iterations-to-target and time-to-target vs the cold baseline, plus the
    break-even step count that amortizes any setup overhead), in the wake-developed regime.
    Phase 2 runs strictly cold (clean baselines); Phase 3 owns the warmstart matrix and
    reuses Phase 2's frozen settings so cold-vs-warm is a controlled comparison.

## Phase gates

| Phase | Deliverable | Status | Approval |
| --- | --- | --- | --- |
| [0](021_rotor_hover_solver_benchmarks/phase_00_availability.md) | All roster configs runnable in the rotor pipeline + instrumentation layer + **sparse near-field ILU preconditioner (W6, added 2026-08-12)** | TECHNICALLY COMPLETE 2026-08-13 (W1–W6; W6 smoke PASS both modes) | clear-context APPROVED 2026-08-13 (2nd review, after remediating the 1st's 6 findings); Ryan re-approval [ ] (Ryan 2026-08-13: proceed to Phase 1 meanwhile) |
| [1](021_rotor_hover_solver_benchmarks/phase_01_consistency.md) | **DESCOPED 2026-08-24 (Ryan): solve the ladder up to R7 and verify the solvers agree on the SOLUTION.** Tuning and per-solver settings moved to Phase 2 | IN PROGRESS (ladder FROZEN R1–R7 8,016→419,276 panels; driver = `rotor_hover_solver_phase1_agreement.jl`, ladder extended to R7, optional `NFCACHE`; **BLOCKER: a dense `backslash_ref` is 335 GiB at R6 / 1.31 TiB at R7 — the R6/R7 reference config needs Ryan**) | [ ] |
| [2](021_rotor_hover_solver_benchmarks/phase_02_single_step_benchmarks.md) | Setup vs per-step cost benchmarks, bottleneck attribution, tuning + authorized source speedups, re-benchmark. **Now OWNS tuning** | IN PROGRESS since 2026-08-17. 2026-08-24: caching wired into tuning AND benchmarks per Ryan's ruling; objective is now a REAL SOLVE (`tune_fmm_perturb cost=`) for R1–R5; **memory budget is a swept AXIS** (`rotor_hover_solver_phase2_tune.jl`), cache build reported in its own column, never in `t_solve`. Nothing has run on the cluster yet; R6–R7 objective needs Ryan | [ ] |
| [3](021_rotor_hover_solver_benchmarks/phase_03_warmstart.md) | Warmstart matrix (none / previous / Taylor-extrapolated) × (Krylov, FGS) in unsteady hover. **Validity guard = per-step BC satisfaction** (Ryan 2026-08-25) | NOT STARTED (gated on Phase 2 Step A; approval unticked). Instrumentation READY: FGS per-step iteration capture CLOSED (the `phase_03_launch_prompt.md` "pending callback plumbing" blocker was superseded 2026-08-23, corrected in place 2026-08-25). Guard rebuilt 2026-08-25 after the 08-23 extrap "divergence" was traced to the TEST, not the solver — the old bar read a residual tolerance as a solution tolerance. Ryan then ruled the whole solution-comparison approach out: measure the BC instead. Verified live at R1/fgs (cold+prev): all steps certified, max\|φ\|/tol 1.007 / 0.994, ΔCT 1.3e-5 %. Costs ~62% of `t_solve` in wall time (excluded from `t_solve`; Ryan: acceptable) | [ ] |
| [4](021_rotor_hover_solver_benchmarks/phase_04_blown_wing.md) | Rotor + NACA 0015 wing multi-body benchmarks; coupled-vs-iterative refactorization asymmetry | NOT STARTED | [ ] |

Completing a phase does not authorize the next phase.

## Acceptance criteria (item level)

1. Every roster config (1a–1f) runs the rotor pipeline and emits schema-conforming CSVs.
2. All solvers demonstrably agree on the no-wake solution across the frozen ladder at
   calibrated matched-residual settings (Phase 1 thresholds in `decision_rules.md`), or the
   disagreement is diagnosed and fixed.
3. For every (solver, mesh rung, threading mode): setup and per-step costs reported
   separately, with RMS/max true residual, iteration counts, convergence-vs-time histories
   (iterative solvers), memory, and fitted scaling exponents.
4. **Warmstart effect quantified per solver — headline deliverable (ruling 12)**: per-step
   iteration and wall-time savings vs the cold baseline at matched residual targets, plus
   break-even step counts, in the wake-developed regime (expected ordering — Backslash:
   null; Krylov: modest; FGS: large — is a hypothesis to test, not assume).
5. Multi-body case reproduces the protocol and quantifies the coupled-vs-iterative
   per-step refactorization asymmetry.
6. Everything needed to regenerate the paper's tables/figures lives in `benchmark/`,
   `data/`, and this item's `figures/` (TikZ + CSV per the figure policy).

## Relationship to other items

- **018** (`018_dji9443_hover_convergence_campaign.md`): supplies the DJI9443 physical setup,
  mesh prescriptions, and operating point; 021 borrows its meshes and ops discipline but
  optimizes *cost*, not physics — no 018 ruling is relaxed here.
- Companion paper section: static wing with rigid wake, benchmarked elsewhere; 021's metric
  definitions (rulings 3–8) should be kept compatible so the two stories share axes.
- `examples/sweptwing_solverbenchmark.jl` is the prior art this item supersedes.

## Decision log

- 2026-08-06 — Ryan: benchmark scope = both isolated solve (primary) and full-timestep share;
  published numbers HPC-only; threading = both 1-thread (BLAS pinned to 1) and 64-thread modes;
  Phase 4 second body = NACA 0015 rectangular wing from `examples/simple_wing.jl`.
- 2026-08-06 — Agent proposals accepted into rulings: shared true-residual evaluator (3),
  per-iteration timestamps instead of assumed equal iteration time (4), memory as a tracked
  metric (8), scaling-exponent fits (9), pre-registered Phase 4 refactorization asymmetry.
- 2026-08-11 — Ryan: **W1 dropped** — campaign uses the legacy linear Kutta condition only;
  upgrade to the 015 pressure-continuity closure only if accuracy/convergence demand it (the
  kutta.jl gate stays; archived design in log). **D2**: FGS convergence history via a
  non-breaking `callback` kwarg added to FastMultipole `solve!` (cross-repo edit to the dev
  checkout). **D3**: FastGaussSeidel per-sweep leaf re-factorization discovered; leaf-LU
  caching deferred to Phase 2b, caveat noted for Phase 1 calibration.
- 2026-08-12 — Ryan: **Phase 0 gate APPROVED** (post-remediation of the 5 clear-context
  review findings; corrected evidence in phase_00 Log 2026-08-12). **RESCINDED later the
  same day**: Phase 0 reopened to add **W6 — sparse near-field ILU preconditioner,
  developed now** (theory doc → implementation design → implementation/testing; roster
  gains config (g) GMRES+ILU). Handoff plan for a fresh agent:
  [`phase_00_ilu_plan.md`](021_rotor_hover_solver_benchmarks/phase_00_ilu_plan.md).
- 2026-08-13 — W6 delivered: Ryan implemented the ILU preconditioner directly
  (2026-08-12 evening), resolving the two staged design questions as **Barba direct-list
  pattern** + **`ILUZero.jl` dependency**; roster **ruling 1 grows to config (g)
  GMRES+ILU**. Session review confirmed correctness, deduplicated the ctor onto
  `_ilu_direct_pattern`, and threaded the assembly (test-gated). Smoke acceptance PASS
  both modes — ILU 290→15 iters, ~19× wall time vs plain GMRES, setup ~1 s (evidence:
  phase_00 Log 2026-08-13). Ryan (2026-08-12): autonomous continuation authorized,
  incl. subagent clear-context gate review, then Phase 1.
- 2026-08-12 — Ryan: GMRES accuracy/preconditioning directives adopted as **ruling 11**
  (right-preconditioning for honesty; Phase-1-calibrated tolerances; FMM-floor measurement
  and expansion-order tuning; preconditioner exploration with **sparse near-field ILU as the
  priority candidate**, theory doc written). Warmstart benefit elevated to a **headline
  deliverable** as **ruling 12** (cold Phase 2 baselines, Phase 3 warmstart matrix at frozen
  settings, break-even metric added).

- 2026-08-25 — Ryan: **Phase 3 validity guard = BC satisfaction, not solution
  agreement.** The 08-23 report of the extrapolated arm as DIVERGED was a test
  defect: the guard compared `sum(abs, strength)` between arms against a hard
  `1e-8`, which (a) read a *residual* tolerance as a *solution* tolerance —
  implicitly $\|A^{-1}\|=1$, where Phase 1 measures $\kappa_{\rm eff}\approx12$ at R1,
  putting the observed 2.4e-7 about three decades *inside* the solver's own
  ambiguity — and (b) used an L1 aggregate that cancels equal-and-opposite
  changes. An intermediate relative-L2-between-arms design was built and then
  **discarded on Ryan's ruling**: *"there is no need to measure the divergence of
  the solution itself... it's easier to show that we're satisfying the boundary
  conditions."* Each arm is now checked per step against the tolerance **it**
  promised, via one `bc_error!` pass (021 ruling 3) — for this Dirichlet body the
  control-point potential *is* the BC residual, so no reference solution, no
  dense $G$ and no conditioning argument is needed at any rung. Accuracy is
  *requested and certified* (absolute FMM error $=0.1\times$ the arm's tolerance)
  rather than guessed via expansion order. **Identity is reported, not asserted**
  — CT is the identity signal, per-step `n_particles` is shown for gross
  divergence, and `n_particles_end` is no longer an assertion: the wake is
  chaotic and small differences compound. Ryan accepted the ~1.4× wall-time
  inflation ("run it twice, or use the replay function"). Full rules in
  `decision_rules.md` → *BC satisfaction guard*.
- 2026-08-25 — Sentinel cleanup: `phase1_agreement.jl` and `phase1_solvetime.jl`
  were the last two drivers hard-coding `niter = -1` for non-Krylov solvers; both
  now mirror `unsteady.jl`. Takes effect on future re-runs only (R1–R7 were
  spooled under the old driver). Amendment recorded that the `niter` **column is
  not homogeneous** — Krylov iterations, FGS sweeps and outer FGMRES iterations
  share it.

## Logging provision

**Do not append session narrative to this file.** Session-by-session notes go to
[`021_rotor_hover_solver_benchmarks/log.md`](021_rotor_hover_solver_benchmarks/log.md); results go
to the phase files and [`ledger.md`](021_rotor_hover_solver_benchmarks/ledger.md). This file is
edited **in place**, never appended to. Only four things change here: the header block, the
phase-gate statuses, the decision log, and the `## Current status` section above.
