# 021 decision rules — metrics, thresholds, protocols

Binding definitions for every phase. Change only via the control doc's decision log.

## Residual (ruling 3)

- True residual $r = G x - b$ in the solve's native units, **unpreconditioned**, computed by
  the shared evaluator (Phase 0 W4) for every solver identically.
- Report both $\mathrm{RMS}(r) = \|r\|_2 / \sqrt{n}$ and $\|r\|_\infty$.
- Matrix-free solvers: one direct influence evaluation at the converged state.
- **FMM floor (ruling 11c)**: per rung, record the direct-evaluated residual of the
  `Backslash` solution passed through the FMM operator; matched-residual targets must sit
  above this floor (tune `expansion_order` if not, and record the cost).
- **Honest stopping (ruling 11a)**: preconditioned Krylov configs use right preconditioning
  (`N=`) only, so the solver's stopping criterion is the true residual.

## Primary metric: certified BC error (Ryan rulings 2026-08-15; supersedes the
## 2026-08-14 PROPOSED knob/tolerance freeze, struck per plan Stage 4)

- **Metric**: boundary-condition error of the solved strengths — for the
  Dirichlet campaign case, the interior perturbation-potential residual at the
  control points, relative L2:
  $\mathrm{BC} = \mathrm{RMS}(\varphi_\sigma + G_\mu x)\,/\,\mathrm{RMS}(b)$
  (identically $\|Gx-b\|$-based, but evaluated in ONE influence pass with both
  strength columns loaded — reference-free, no dense G at any rung).
  **Target: BC ≤ 1e-6 per rung.** Frozen BC definition: perturbation gauge
  (φ∞ excluded), σ from the apparent velocity, panel kerneloffset.
- **Evaluator** (`bc_error!`, `benchmark/common.jl`; validated vs direct at
  R1–R3, `bcerror.csv` — discrepancy ≤1% at every resolvable point, evaluator
  floor ~4e-9): FMM with `PowerAbsolutePotential(0.1 × 1e-6 × rms_b)` and
  dynamic-P cap 20. **Certification = the returned `error_success` flag**
  (every M2L met its bound at P ≤ cap); a row with `error_success=false` is
  uncertified — raise the P cap and re-evaluate.
- **Per-run guard**: every benchmark solve records its certified BC error as
  a standard post-solve column. This guard replaces bound-certification
  restrictions on operating points (the 2026-08-14 MAC ≥ 0.6 rejection is
  WITHDRAWN): any tuner point whose runs verify is admissible.

## Per-rung tuning procedures (NO campaign-wide knob freeze — ruling 2)

FMM knobs are tuned per rung per solver family for optimal performance.
Objective (carried judgment, plan ruling 5): per-solve wall time subject to
certified BC ≤ target; setup timed separately (production reuses solver
objects; Phase 3 owns amortization).

1. **Shared Krylov apply knobs** (one set per rung for gmres/jacobi/ilu/
   fgmres outer applies — ruling 3): `tune_fmm` at
   `PowerAbsolutePotential(0.1·1e-6·rms_b)`, max_expansion_order 20, MACs
   0.3–0.7, tuned with a solved strength state loaded. R1–R3 values (tune.csv
   2026-08-14): R1 p17/MAC .5/leaf 21, R2 p17/.5/12, R3 p18/.5/18.
2. **FGS knobs** (purpose-built τ-ladder, `..._fgstune.jl` — tune_fmm cannot
   model leaf-as-GS-block-size or inner_iterations): coordinate descent over
   {p, MAC, leaf, inner} on cost-to-BC-1e-6; one instrumented cold solve per
   candidate yields the full (wall time, certified BC) staircase per outer
   iteration; per-τ winners (τ ∈ 1e-1…1e-6) selected over the whole pool
   (fgstune_selected.csv). Every selected point re-verifies BC ≤ τ on a
   fresh cold solve (fgstune_verify.csv).
3. **Shared FGS knob set, solver ↔ preconditioner (Ryan ruling 2026-08-17,
   supersedes the end-to-end knob-set selection)**: per rung, the shared
   FGS knob set = the **τ=1e-6-tuned config** from the ladder. Rationale:
   coarser-τ-tuned configs have an apply-accuracy BC floor that can sit
   ABOVE 1e-6 (measured R2: the τ=1e-4 config's internal residual falls to
   5.6e-10 while true BC error plateaus at 3.6e-6) — they can never serve
   the solver role; tuning at 1e-6 rules this out by construction
   (`stage3_winner`'s crossing check remains as a tripwire only). The
   preconditioner role wraps the same config with FIXED sweeps (fixed
   counts + zero seed keep the apply linear — tolerance-stopped applies are
   banned); the SWEEP COUNT is still selected by end-to-end FGMRES wall
   time to BC 1e-6 over that config's own staircase crossings
   (fgsprecond.csv).

## FGS reporting at the 1e-6 target (ruling 4 — dual-row rule)

FGS stops on its internal absolute max-abs residual; per-iteration decades
mean it can overshoot the BC target. If the converged row's BC error
overshoots 1e-6 by more than a half decade, report BOTH rows:
`row_kind=target_1e-6` (cost at the target) and `row_kind=last_above_1e-6`
(cost at the last iterate above 1e-6, its BC error snapped to the nearest
half decade). Implemented in `..._table.jl` (solvetable.csv).

## Iteration counting (Phase 3 amendment, 2026-08-22)

Phase 3's headline metric is iterations-to-target per step, so the count must
mean the same thing for every solver.

- **Definition: `niter` = iterative work units actually performed** — Krylov
  iterations for `KrylovSolver` (`ws.stats.niter`, unchanged), Gauss-Seidel
  **sweeps** for `FGSSolver`. `Backslash` records `niter = -1`: a direct solve
  has no iteration count, and that null is itself a reported Phase-3 result.
- **FGS off-by-one.** Upstream FastMultipole keeps its count as a loop local and
  exposes it only through the per-iteration `callback`, which fires at the TOP
  of iteration $k$ on a residual reflecting $k-1$ completed sweeps, *before* the
  tolerance check — hence upstream's own "converged after `iteration-1`
  iterations" (`FastMultipole/src/solve.jl:1249`). FLOWPanel therefore records
  $$\texttt{niter} = \begin{cases} k - 1 & \text{converged at callback iteration } k \\ k & \text{tolerance never met}\end{cases}$$
  (the non-converged branch counts the sweep that still runs after the last
  failed check). Exposed as `FGSSolver.niter` / `.solved`, alongside
  `KrylovSolver.niter` / `.solved`, so the harness reads one uniform accessor.
- **Consequence for comparability:** an FGS sweep and a Krylov iteration are NOT
  the same unit of work. Iteration counts are comparable *within* a solver
  across warmstart types — which is exactly what Phase 3 measures — and
  wall-time-to-target remains the only cross-solver cost currency.
- **`niter` is the LAST inner solve of the step; `niter_first` is the metric**
  (amended 2026-08-23). `solver.niter` / `.solved` are kernel-level: they
  describe whichever `_solve!` ran most recently. When a step contains several
  per-body solves, sampling them after the step reports the *last* one, which
  has already been warmed by its predecessors within the same step. Phase 3
  therefore reads **`step_niter_first(solver)`** — the first per-body solve of
  the step, the one that actually consumes the step-to-step initial guess — as
  the warmstart headline metric, and keeps `niter` only as a legacy diagnostic
  for continuity with Phase 1/2 rows.
- **The `niter` COLUMN is not homogeneous (amendment 2026-08-25).** Three
  different units share it, and a reader comparing them as "cost" will be
  wrong:
  1. `krylov_*` rows — Krylov iterations;
  2. `fgs` rows — Gauss-Seidel **sweeps**;
  3. `fgmres_fgs` rows — the **outer FGMRES** iteration count. Its
     `FGSPreconditioner` holds its own `FGSSolver` with `sweeps` FIXED and
     `tolerance=0` (`src/FLOWPanel_solver.jl:1772-1787`), so the inner sweep
     count is a configured constant and is surfaced nowhere.

  Sentinel cleanup, same date: `phase1_agreement.jl:439` and
  `phase1_solvetime.jl:175` were the last two drivers hard-coding `niter = -1`
  for non-Krylov solvers; both now mirror `unsteady.jl:241`. Every other driver
  in `benchmark/` already read `solver.niter` unconditionally. **The fix reaches
  only future re-runs** — R1–R7 were spooled by Slurm under the old driver and
  stay uniformly `-1` for FGS.

- **`nsolves` and step-level `solved`.** `step_nsolves(solver)` reports how many
  per-body solves the completed step took, and `step_solved(solver)` is the AND
  of those solves' inner convergence flags (it is *not* a multibody outer
  convergence flag — that stays with the `ConvergenceHistory` kwarg and its
  verbose warning). Both return `-1` / `true` for solvers that carry no step
  statistics (`Backslash` and friends); `-1` means **unavailable**, never a
  passed check.
- **Block-GS outer counts** (spec: "record outer counts too") remain a documented
  **non-deliverable for the rotor case** — but NOT because the loop is
  degenerate. **Correction (2026-08-23):** the earlier claim that a 1-tuple's
  outer loop was degenerate was wrong. `sources` is empty for a lone body, so the
  block-GS map has no coupling update to apply, yet the loop still re-solved the
  same system on every outer iteration until `max_delta` fell below tolerance —
  a standing repeated-solve tax, and (with warm-start history advancing per
  `_solve!`) a mechanism that could keep `max_delta` from ever settling. The
  one-body case now **breaks after the first block solve by construction**, so a
  1-tuple records exactly one outer entry and there is nothing to count.
  Multibody outer counts are deferred to Phase 4 (blown wing), for which
  `solve!(bodies::Tuple, solvers::Tuple)` accepts `solver_optargs` to reach the
  inner per-body solve.

## BC satisfaction guard (Phase 3, Ryan's ruling 2026-08-25)

Supersedes the `sum(abs, strength)` vs `1e-8` check, the same-day `1e-8 -> 1e-7`
ruling (withdrawn before implementation), and an intermediate relative-L2
between-arms design (built, then discarded in favour of this).

**Ryan's ruling.** "There is no need to measure the divergence of the solution
itself between warm and cold solves; rather, the residual should be measured...
it's easier to show that we're satisfying the boundary conditions."

- **Why the old guard was invalid.** The solver's acceptance test is
  $\max_i |b - Ax|_i \le \texttt{tol\_abs}$ — absolute and independent of the
  initial guess (`FastMultipole/src/solve.jl:1222`; `residual!` returns max-abs
  despite the `mse` label). Cold and warm call the identical
  `FastMultipole.solve!`; the only difference is `project_solution!`. They land
  at **different points of the same residual ball**, of solution-space diameter
  $\|A^{-1}\| \cdot 2\,\texttt{tol\_abs}$. The old bound assumed
  $\|A^{-1}\| = 1$ — it read a *residual* tolerance as a *solution* tolerance.
  Phase 1 measures $\kappa_{\rm eff} \approx 12$ at R1, so the expected spread is
  near $2\times10^{-4}$ and the observed $2.4\times10^{-7}$ sat ~3 decades
  *inside* the solver's own ambiguity. **extrap was never diverging.** Second
  defect: `sum(abs, x)` is one aggregate and cancels equal-and-opposite changes.
- **What replaced it.** Each arm is checked against the tolerance **it**
  promised, per step, by one `bc_error!` pass (`benchmark/common.jl`, 021 ruling
  3). For this Dirichlet body the control-point potential IS the BC residual
  $\varphi_\sigma + G_\mu x$, so no reference solution, no dense $G$, and no
  conditioning argument is needed at any rung. Recorded per step:
  `bcerr_{max,min,q1,med,q3,rms,tol,eps,certified}`, `t_bcerr`.
- **Accuracy by construction, not by guess.** `bc_error!` *requests* an absolute
  FMM error target and reports whether every M2L pair certified it. The driver
  requests `BCERR_SAFETY` (default 0.1) $\times$ the arm's own stopping
  tolerance, making the instrument ~10x sharper than the acceptance level;
  `bcerr_certified=false` says loudly when it was not.
- **Pass band is $1 + \varepsilon/\texttt{tol}$, not 1.** The arm stopped on its
  own leaf-local estimate at its tuned expansion order; `bc_error!` re-evaluates
  the same quantity independently and more sharply, so the two legitimately
  differ by about the error requested of the measurement. Ratios inside that
  band report as **marginal** rather than OK, so drift toward the bar stays
  visible. Above it: **BC VIOLATED**.
- **Identity is REPORTED, never asserted.** The wake is chaotic — tiny
  differences compound over many steps, so demanding identical particle counts
  or trajectories would flag correct runs. **CT** is the physical identity
  signal and carries no threshold; `n_particles` is emitted every step so gross
  divergence is visible. `n_particles_end` is no longer an assertion.
- **Cost, and what it does and does not contaminate.** The pass measured ~62% of
  `t_solve` at R1 (cold-start smoke), inflating wall time ~1.4x. **Ryan ruled the
  extra wall time acceptable** (run twice, or recover it through replay), so the
  measurement runs INLINE in the same run — `BCERR_EVERY` defaults to 1 (every
  step) and exists only as an escape hatch. Skipped steps write EMPTY cells,
  never a fabricated value or verdict.
  - `t_solve` — the headline warmstart metric — is **clean**: its timer closes
    before the BC pass begins.
  - `t_step_total` is **not** clean: the pass lies inside the
    `maneuver!`-to-`maneuver!` window. Since full-timestep share is a campaign
    deliverable (Ryan 2026-08-06), the driver emits **`t_step_net =
    t_step_total - t_bcerr`** and reports its end-of-run solve share against net
    wall, quoting the raw figure and the measurement total alongside. Measured
    R1 impact: per-step share 35% raw vs 45% net — a ~10-point error if the raw
    column is used unwittingly. **Use `t_step_net` for any timestep-share claim.**
- **Verified on a live R1/fgs run (cold + prev, 5 steps each) 2026-08-25:** all
  steps certified; `max|φ|/tol` = 1.007 (cold) and 0.994 (prev), i.e. both arms
  sit right at the bar they promised; $\Delta$CT between arms = $1.3\times10^{-5}$ %.
- **Neumann bodies are not covered.** `bc_error!` implements the Dirichlet
  metric and errors on a Neumann body; the driver asserts Dirichlet up front. A
  Neumann rotor would need `true_residual!` and a re-derived tolerance.

## Solve-step lifecycle and warm-start history (2026-08-23)

Warm-start history must hold **consecutive timestep solutions**, so it advances
on exactly one event: a **completed top-level solve step**. Three events are now
distinguished in `src/FLOWPanel_solver.jl`:

1. a raw linear-solver kernel call (`_solve!`, `_krylov_launch!`) — updates only
   `solver.niter` / `.solved`;
2. one per-body solve inside a step — `note_step_solve!`;
3. a completed top-level step — `save_step_solution!` (history) then
   `publish_step_stats!` (statistics), together as `finalize_step_solution!`.

Step boundaries live at the orchestration sites only: the two public single-body
`solve!` methods, the tuple `solve!` (which brackets the whole outer loop with
one begin/finalize per solver and calls the single-body methods with the internal
`finalize_step=false`), and `TraceCorrected`'s direct `_solve!` bypass (finalized
after `set_wake_correction!`, not inside the solve's `try`). Raw `_solve!` callers
— the Phase 1/2 benchmarks, the solver-history tests, and `FGSPreconditioner`'s
inner FGS solve — stay deliberately outside this accounting.

Failure semantics: a step that throws leaves the previously published statistics
and the history untouched; the next `begin_step_solution!` discards the partial
transient state. A step that *completes* but whose inner solver hit its iteration
limit still saves (unchanged from before) — the history means "one completed step
state", not "a converged state"; the published `solved` flag is what the driver
rejects on.

**Memory (ruling 8 impact):** `SolveStepStats` is 48 bytes
(`Base.summarysize`), plus an 8-byte reference in each `FGSSolver` /
`KrylovSolver` — 56 bytes per solver, negligible against `solver_state_bytes`
(MB–GB scale) but recorded here because solver memory is a reported metric.

## Phase 3 unsteady CSV schema (2026-08-22; widened 2026-08-25)

`unsteady.csv` uses a bespoke per-step schema, NOT `RUNS_CSV_COLUMNS`, and is
not covered by `validate_runs_csv`. Phase 3 extends it with:
`niter_first, nsolves, solved, warmstart, warmstart_order, restart_step,
skip_steps, strength_checksum`, alongside the existing per-step columns.
The 2026-08-25 BC-satisfaction guard added:
`bcerr_max, bcerr_min, bcerr_q1, bcerr_med, bcerr_q3, bcerr_rms, bcerr_tol,
bcerr_eps, bcerr_certified, t_bcerr, t_step_net, n_particles, CT`.

- `bcerr_*` (2026-08-25) are order statistics of the per-panel Dirichlet BC
  residual from one `bc_error!` pass per step. `bcerr_tol` is **that arm's own**
  effective absolute acceptance level (FGS: `solver.tolerance`; Krylov:
  `atol + rtol*||b||`; Backslash: NaN), so `bcerr_max <= bcerr_tol` is the arm's
  promise checked independently. `bcerr_eps` is the absolute FMM error
  *requested* of the measurement (`BCERR_SAFETY * bcerr_tol`) and
  `bcerr_certified` whether every M2L pair proved it. `t_bcerr` is the
  measurement's own wall time, excluded from `t_solve` but NOT from
  `t_step_total`; `t_step_net = t_step_total - t_bcerr` is the corrected
  full-timestep figure and is what any timestep-share claim must use.
  `BCERR_EVERY > 1`
  subsamples: skipped steps write EMPTY cells, never a fabricated value.
- `n_particles` and `CT` (2026-08-25) are per-step identity **reporting**, not
  assertions — see the BC satisfaction guard section for why the chaotic wake
  makes exact identity the wrong demand. `CT` needs `RUN_MONITORS=true`, which
  the driver defaults on under `PHASE=phase3*` (Bernoulli-only); Phase 2 rows
  are unaffected.

- The driver refuses to append to an `unsteady.csv` whose first line is not the
  current header, so pre-Phase-3 (narrower) rows can never be mixed under one
  header; `benchmark/phase3_analysis.jl` likewise errors on inputs missing
  `niter_first` / `nsolves`.

- `warmstart` ∈ `cold | prev | extrap`; `warmstart_order` is the effective
  polynomial order (0 for `cold` and `prev`).
- `skip_steps` records how many leading steps the analysis must drop from
  wake-developed averages — the warmstart history needs filling after a restart
  (order $p$ needs $p+1$ solves). It is recorded, never applied destructively.
- **Accuracy guard.** A certified BC error is deliberately NOT recorded per
  unsteady step: the frozen BC metric above is defined against a static
  single-step RHS, and the unsteady RHS moves every step, so reusing that column
  would misreport it. The guard is instead (a) the per-step `solved` flag, which
  catches a warm guess that lets a step stop short of tolerance, and (b) a
  trajectory-identity check. **Amended 2026-08-25 (Ryan):** that check is
  replaced by the per-step BC-satisfaction check — each arm against the
  tolerance it promised — see *BC satisfaction guard* above for why the
  checksum form was invalid and why solution-space comparison was dropped.
  `n_particles_end` is no longer an assertion: the wake is chaotic and small
  differences compound. (It remains the discriminator that resolved the
  filament-family question at 328 vs 329 particles, where the difference was
  structural rather than accumulated.)

## GENERATION PIN — LineGauss generation (Ryan's ruling 2026-09-05)

This is **generation break #2** for this campaign; #1 is recorded in `ledger.md`
(2026-08-24, new FMM + isolated checkout). Rows from different generations are kept
for the record but **must never be plotted on the same axes**.

**The pinned triple.** Every campaign run from 2026-09-05 — Phase 2, Phase 3, and any
Phase 1 re-run — uses these three commits, as **detached git worktrees**, and no
others:

| repo | branch of origin | pinned commit | tag |
| --- | --- | --- | --- |
| FLOWPanel.jl | `fastmultipole` | the commit tagged below | `021-lg-gen2` |
| FastMultipole | `flowpanel-20260817` | `0ce3ba60` | `021-lg-gen2` |
| FLOWVPM.jl | `flowpanel` | `a627dd9` | `021-lg-gen2` |

**Resolve the pin by TAG, not by the SHAs above.** All three repos carry the
annotated tag `021-lg-gen2`, and `git worktree add --detach 021-lg-gen2` in each is
the whole setup step. The tag exists because this document cannot name its own
commit: FLOWPanel's pinned tree is the one *containing this file*, so a literal SHA
here would always be the parent commit and would silently pin a tree without its own
pin record. The FastMultipole and FLOWVPM SHAs are listed because they are fixed
independently of this file; FLOWPanel's is not. FLOWPanel's code content is identical
to `63c7f10` — the tagged commit adds only this document and `agent_policies/HPC.md`.

- **All three trees were clean when pinned** (verified 2026-09-05: zero modified
  tracked files in each), so each SHA fully determines its tree. No `-dirty` markers,
  unlike every prior generation of this campaign.
- **One triple for ALL phases.** Mixing SHAs across phases would reintroduce exactly
  the confound found on 2026-09-05 between the Gaussian R1–R5 tuning (FastMultipole
  `4c0f1b8`) and any LineGauss rung: LineGauss did not exist in FastMultipole until
  14 commits after that pin, so a family change and an FMM-build change would have
  moved together and the cost-scaling exponent would have been uninterpretable.
- **`Manifest.toml` is deliberately NOT committed.** The campaign environment's
  Manifest must dev-point at the campaign worktrees, so a committed Manifest would
  bake in one machine's absolute paths and make the pin misleading. The triple of
  SHAs is the pin; the Manifest is generated per worktree via `Pkg.develop`.
- **Provenance is now recorded for all three repos.** `benchmark/common.jl` emits
  `commit`, `fm_commit`, **`vpm_commit`** (new — FLOWVPM was unrecorded for the entire
  campaign before this date, so every earlier row has an indeterminate FLOWVPM
  version) and a **`worktrees`** column of `tag=<dir>@<branch-or-DETACHED>`. All three
  resolve via `pkgdir(loaded_module)`, NOT a fixed sibling path — under a worktree the
  sibling path would silently describe the wrong checkout. A row that reads
  `@DETACHED` is reproducible; one naming a live branch is not, and now says so.
- General policy for campaign worktrees lives in `~/.claude/CLAUDE.md`
  ("Campaign Reproducibility") and cluster mechanics in `agent_policies/HPC.md`;
  neither is restated here.

## Filament regularization — LineGauss (Ryan's ruling 2026-09-05)

- **`LineGaussRegularization` is the campaign default from 2026-09-05**, following the
  codebase default flip of 2026-08-29 (task 052d). All nine 021 launchers had pinned
  `gaussian` since 2026-08-22 — a pin that predates the LineGauss ruling, not a
  considered rejection of it. Now flipped to `linegauss`; `FLOWPANEL_FILAMENT_REG=gaussian`
  reproduces pre-2026-09-05 rows.
- **Phase 2 is RELAUNCHED under LineGauss at every rung R1–R7.** The 2026-08-26..29
  Gaussian descents (R1–R5) are superseded as cost numbers and retained as the SEED
  source and as history. Rationale: LineGauss changes the per-edge kernel cost (4 erf +
  1 exp vs Gaussian's 1 expm1) AND `radius_inflation` (≈6.25σ vs 5.90σ at tol 1e-6),
  so direct-list sizes and therefore the cost optimum move. A ladder half-Gaussian and
  half-LineGauss would confound the scaling exponent, which is Phase 2's deliverable.
- **Phase 1's harvested table is NOT invalidated**: `phase1_case.jl` is a frozen
  single-step solve and the family was measured inert there (2026-08-22). It is
  nevertheless now disclosed as Gaussian-provenance in the table header.

### Descent seeding (2026-09-05)

- The tuner gained `TUNE_SEED` / `TUNE_SEED_B0` (`"P:MAC:LEAF"`), read in
  `rotor_hover_solver_phase2_tune.jl`. Unset = the previous `TUNED`/`REF_START`
  behaviour, unchanged. The seed is recorded in the `notes` column of every row, so a
  row can never hide where its descent started.
- **Two seeds are required because the optimum splits hard by cache regime**, measured
  in the Gaussian descents: uncached (budget 0) converges to **leaf ≈ 6** (R4 6, R5 6;
  MAC 0.60–0.65), cached (budget > 0) to **leaf ≈ 32** (R4 32/32, R5 32/32/32; MAC
  0.50–0.60). Those attractors are 5× apart and no single seed serves both. The
  driver's fallback `REF_START = (17, 0.5, 21)` is wrong for both at high rungs.
- `expansion_order` shows **no trend** and scatters 10–17 even within one rung across
  budgets — the cost surface is flat in P, so the P seed matters little.
- **Seeding does not weaken any guarantee.** The descent still searches and the
  BC ≤ 1e-6 certification gate still decides admissibility, so a stale or wrong seed
  costs optimality, never correctness. Seeds carried over from Gaussian are a
  *hypothesis* about where LineGauss's optimum sits, not a measurement of it.
- **R3 budget 0 (Gaussian) carries `tune_timed_out=true`** — best-so-far, not an
  optimum (§4.1 trap). Excluded from the seed fit and never publishable as a tuned
  result.

### Trace provenance now includes the filament family (2026-09-05)

`load_trace`'s hard guard was `(rung, mem_budget_gib, tune_reps, tune_abandon_factor,
hardware_tag)` with a soft warning on `fm_commit` — **the filament family was absent**.
Since a memoized `t` is a timing and the family changes the objective that timing
measures, a Gaussian trace would have replayed silently into a LineGauss descent and
corrupted it with no error raised. `filament_reg` is now a 13th trace column inside the
HARD guard, and the schema change makes every 12-field Gaussian trace fail loudly at
the header check ("move it aside before resuming") rather than being partially read.

**Relaunch precondition:** the existing R1–R5 `tune_phase2.csv` rows must be moved
aside before resubmission. Row-level resume keys on `(rung, mem_budget_gib)` and would
otherwise skip every already-present rung, producing no LineGauss rows at all.

## Force recovery — Bernoulli-only (Ryan's ruling 2026-09-05)

- **For the remainder of the item, CT is recovered by steady `PressureBernoulli` +
  Force only.** `PressureLaplace` and `KuttaJoukowski` are retired from the forward
  measurement path. This aligns the whole item with what Phase 3 already did
  (`RUN_MONITORS=true`, Bernoulli-only under `PHASE=phase3*`); Phase 2 rows carry no
  CT and are unaffected.
- **Consequence — retired knob**: the `PressureLaplace` CG `itmax` per-rung tuning
  entry (hit its 1000-iteration cap at R3; see `phase_01_hpc_procedure.md`) no longer
  applies to any forward run. Nothing else in the per-rung tuning procedure changes.
- **The already-harvested Phase 1 Table 4 keeps all three columns** — it is history,
  and the `CT_kj` / `CT_laplace` columns are what localised the R7 defect. Footnote it
  as superseded going forward rather than restating it.
- **The R7 `CT_laplace` collapse is DEMOTED to non-blocking** (was Phase 1 next-action
  #3). It no longer gates any 021 deliverable. It remains a real `PressureLaplace`
  defect on a 419,276-panel mesh and is logged as a standalone defect in `ledger.md`,
  not dropped.
- **Redundancy trade, accepted knowingly (Ryan)**: three recoveries are what caught a
  common-mode post-processing defect that on a single recovery would have read as a
  physical trend. On one recovery that check is gone. Bernoulli also carries the
  moving-body inertial-KE caveat — sound for the cross-config deltas Phase 1 used it
  for, doing more work in Phase 3's within-config time series. **If a CT anomaly needs
  debugging later, pull `KuttaJoukowski` and/or `PressureLaplace` back in for that
  diagnosis** — they stay available, they are just not in the standing path.
- **Bernoulli's own R7 behaviour, for the record**: `CT_bernoulli` does NOT diverge at
  R7 (0.053203 vs 0.053391 at R6, −0.35%) — no sign flip, no collapse, and the flattest
  of the three recoveries across the ladder (R1→R7 total spread 0.85%, vs ~4.5% for KJ).
  But that R6→R7 step is ~17× the R3→R6 steps (~0.02–0.04%) after three flat rungs, and
  it cannot be a config-set artifact of R7's dropped `krylov_jacobi` / split provenance
  (configs agree on CT to ~1e-3%, three decades too tight). Treat it as a genuine
  mesh-rung effect worth a note, not a clean point.

## Historical agreement stage (COMPLETE 2026-08-14; thresholds retired)

The reference-based strength-agreement thresholds (relL2 ≤ 1e-5 provisional;
3e-5 recalibration proposal) are SUPERSEDED by the certified BC metric
(ruling 1) — retained here as history: agreement.csv R1–R3 × 7 configs
measured worst |dCT| ≈ 8e-4 % (120× inside the 0.1% cross-check threshold,
which REMAINS in force for any future physical cross-check), and every
config's true residual honestly met rel 1e-6. Any certified-BC failure at
matched settings = bug; fix before Phase 2.

## Timing protocol (rulings 5–7; amended by Ryan 2026-08-18)

- Published numbers: single dedicated HPC node, exclusive allocation
  (enforced: `--constraint=zen3 --exclusive`, m12-class 2×64-core AMD Zen3
  512 GB; `hardware_tag` records the node).
- Isolated solve: ADAPTIVE min-of-k (amends the fixed k ≥ 5): one warmup
  solve (excluded, also selects k) then k timed reps — k=5 below 60 s,
  k=3 below 10 min, k=2 above; `k_reps` column records the k used.
- Cost-ceiling drop-out (extends the ledger's dense-config rule):
  `krylov_gmres` drops from SINGLE-thread mode at R6–R7 (projected
  ~12 h/solve at R7); dense configs drop at R6–R7 both modes. Dropped rows
  can be run later if they matter.
- The per-row direct O(N²) residual is RETIRED from the table protocol
  (Ryan 2026-08-18): redundant with the certified BC column (≤1% agreement,
  bcerror.csv R1–R3); `DIRECT_RESID=1` re-enables for spot checks.
- Setup components timed separately: assembly, factorization, tree build, preconditioner.
- Full-timestep share: one `simulate!` step per config, reported as solve-time / step-time.
- Threading modes (ruling 6): (1 Julia thread + `BLAS.set_num_threads(1)`) and (64 Julia
  threads + recorded consistent BLAS count). Harness asserts and logs both at startup; never
  mix modes in one comparison.

## Memory (ruling 8)

- Report per config: resident solver-state size via `solver_state_bytes(solver)`
  (**finalized 2026-08-12**, `benchmark/common.jl`): `Base.summarysize(solver)` minus the
  summarysize of every referenced `AbstractBody` (found recursively through nested
  operators/preconditioners/tuples) — a comparable object boundary across solver types.
  In: matrices, factorizations, FMM trees/buffers, Krylov workspaces, scratch, histories.
  Out: the body. `summarysize` deduplicates shared references (aliased `G`/`Glu.factors`
  counted once; each body subtracted exactly once).
- Plus allocation during one solve (`@allocated`).

## CSV schema (finalized in Phase 0 W4, provenance columns added 2026-08-12;
## `filament_reg` added 2026-08-22)

`runs.csv` — one row per solve:
`run_id, phase, solver_config, mesh_file, n_panels, threading_mode, julia_threads,
blas_threads, t_assembly, t_factorize, t_tree, t_precond, t_rhs, t_solve_min, k_reps,
iterations, rms_residual, max_residual, mem_state_bytes, alloc_solve_bytes, commit,
fm_commit, julia_version, hardware_tag, filament_reg, solver_settings, backend_settings, notes`

- **`filament_reg` (SCHEMA AMENDMENT, Ryan 2026-08-22)**: the filament regularization
  family in force (`FLOWPanel.FILAMENT_REGULARIZATION[]`), recorded and never modified by
  the harness. Added when the campaign moved from Vatistas to Gaussian and it emerged that
  NOTHING recorded the family — leaving several 2026-08-20/21 wake-carrying arms
  permanently indeterminate (BRAINSTORM/025 landed selectability mid-day on 08-20, so the
  date column cannot resolve them, and every arm shares one `-dirty` commit SHA). Emitted
  by `assert_and_banner()` in `benchmark/common.jl` so all 021 drivers inherit it, and
  required non-empty by `validate_runs_csv`. Rows written before this date carry no
  family column; treat their family as indeterminate unless the run predates 2026-08-20
  (hardcoded Vatistas, not selectable).
- `commit`/`fm_commit`: FLOWPanel and FastMultipole (dev checkout) SHAs with `-dirty`
  markers; `solver_settings`/`backend_settings`: flattened `k=v;…` from the metadata
  dicts (tolerances, itmax, preconditioner knobs, FMM knobs). The startup banner is also
  written to `banner.txt` next to `runs.csv` — provenance never lives only in stdout.
- **Cold-solve protocol (2026-08-12, review finding 1)**: solver-visible state
  (`body.strength`) is reset to zero, untimed, before every timed/history/allocation
  solve; FGS runs carry a tripwire assert that the first cold-solve residual is far above
  tolerance. Warm behavior is measured only where designed (Phase 3).

`history_<run_id>.csv` — per-iteration: `iter, t_wall, residual_internal, residual_true?`
(true residual per iteration only where cheap; internal metric labeled as such).

## Plot specs

- Convergence: residual (log) vs wall-clock time, all iterative solvers on shared axes, one
  panel per (mesh rung × threading mode).
- Cost: log–log time vs N per cost component, with fitted exponents annotated; memory vs N
  alongside.
- Figures follow the global figure policy: standalone TikZ `.tex` + same-named CSV data dir
  under `021_rotor_hover_solver_benchmarks/figures/`, `pdflatex`-compilable.

## Run judging

Judge every run from its CSVs, not stdout; verify knobs from the log banner (018 lesson).
