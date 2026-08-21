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

## CSV schema (finalized in Phase 0 W4, provenance columns added 2026-08-12)

`runs.csv` — one row per solve:
`run_id, phase, solver_config, mesh_file, n_panels, threading_mode, julia_threads,
blas_threads, t_assembly, t_factorize, t_tree, t_precond, t_rhs, t_solve_min, k_reps,
iterations, rms_residual, max_residual, mem_state_bytes, alloc_solve_bytes, commit,
fm_commit, julia_version, hardware_tag, solver_settings, backend_settings, notes`

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
