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

## Phase 3 unsteady CSV schema (2026-08-22)

`unsteady.csv` uses a bespoke per-step schema, NOT `RUNS_CSV_COLUMNS`, and is
not covered by `validate_runs_csv`. Phase 3 extends it with:
`niter_first, nsolves, solved, warmstart, warmstart_order, restart_step,
skip_steps, strength_checksum`, alongside the existing per-step columns.

- `niter_first` (added 2026-08-23) is the iteration count of the FIRST inner
  solve of the step and **is the warmstart comparison metric**. Legacy `niter`
  is the LAST inner solve and is retained as a diagnostic only.
- `nsolves` (added 2026-08-23) is the per-body solve count of the step. For the
  one-body rotor runs it must be `1` on every wake-developed row; anything else
  means the timestep re-solved and the row set is an invalid warmstart
  comparison. `-1` means unavailable (direct/generic solver), not a pass. The
  check is deliberately *not* generalized to multibody, where genuine block-GS
  outer counts vary per timestep.
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
  trajectory-identity check — cold and warm arms of the same config solve to the
  same tolerance and so must agree on `n_particles_end` and
  `strength_checksum`. Divergence means the guess changed the answer rather than
  just reaching it faster. (This is the discriminator that resolved the
  filament-family question at 328 vs 329 particles.)

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
