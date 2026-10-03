# 034 decision rules — metrics, thresholds, protocols

Binding definitions for every phase. Change only via the control doc's decision log.
Adapted from `BRAINSTORM/021_rotor_hover_solver_benchmarks/decision_rules.md`
(2026-10-02); where a rule is copied verbatim the 021 file remains the provenance
record for *why* it exists.

## Primary metric: certified BC error (static-RHS benchmarks, Phases 1–2)

- Metric: relative-L2 Dirichlet BC residual of the solved strengths,
  $\mathrm{BC} = \mathrm{RMS}(\varphi_\sigma + G_\mu x)/\mathrm{RMS}(b)$, evaluated in
  one FMM influence pass by `bc_error!` (`benchmark/common.jl`). **Target BC ≤ 1e-6
  per rung; certification = the returned `error_success` flag** (every M2L met its
  bound at P ≤ cap). Uncertified rows are re-evaluated at a higher P cap, never
  published as-is.
- **Dirichlet-only caveat:** `bc_error!` errors on a Neumann body. Phase 0 must
  confirm the pitching-wing body is Dirichlet under the pinned
  `VelocityThroughSources` formulation before this metric is adopted; otherwise a
  `true_residual!` path and a re-derived tolerance are required (and this file is
  amended via the decision log).

## Unsteady guard (Phase 3)

- A certified static BC error is NOT recorded per unsteady step (the RHS moves every
  step). Instead, per step, each arm is checked against **the tolerance it promised**
  by one `bc_error!` pass: columns `bcerr_{max,min,q1,med,q3,rms,tol,eps,certified}`,
  `t_bcerr`. Pass band is $1 + \varepsilon/\mathrm{tol}$ (marginal above 1.0, VIOLATED
  above the band) — 021 ruling 2026-08-25.
- `t_solve` is clean (timer closes before the BC pass); `t_step_total` is not —
  **use `t_step_net = t_step_total − t_bcerr` for any timestep-share claim.**
- Physical identity is REPORTED, never asserted: the CL hysteresis loop (and per-step
  CL/CM) is 034's identity signal, analogous to 021's CT column. The growing panel
  wake is deterministic (no chaotic particle field), so closer agreement is expected
  than in 021 — but thresholds are still set only after Phase 1 measurements.

## Cold-solve protocol

- **Cold = zero-initial-guess solve** (Ryan 2026-09-23), never fresh-process.
- Solver-visible state (`body.strength`) reset to zero, untimed, before every
  timed/history/allocation solve; FGS runs carry the tripwire assert that the first
  cold residual sits far above tolerance.
- Adaptive min-of-k (021 amendment 2026-08-18): one excluded warmup selects k —
  k=5 below 60 s, k=3 below 10 min, k=2 above; `k_reps` recorded.
- Warm behavior is measured only where designed (Phase 3).

## Iteration counting

- `niter` = iterative work units actually performed: Krylov iterations for
  `KrylovSolver`, Gauss–Seidel **sweeps** for `FGSSolver` (with the upstream
  callback off-by-one handled per 021: converged at callback iteration $k$ →
  $\texttt{niter}=k-1$); `Backslash` records `niter = -1` (no iteration count — a
  reported null, and `-1` always means unavailable, never a passed check).
- The `niter` column is NOT homogeneous across solver families and is never a
  cross-solver cost comparison; **wall-time-to-target is the only cross-solver
  currency.**
- Phase 3 headline reads **`step_niter_first(solver)`** (the per-body solve that
  actually consumes the step-to-step initial guess), with `niter` kept as a legacy
  diagnostic. `step_nsolves` / `step_solved` per 021 semantics.
- `warmstart ∈ cold | prev | extrap`; `warmstart_order` = effective polynomial order
  (0 for cold/prev); `skip_steps` records leading steps the analysis must drop while
  the warmstart history fills (recorded, never applied destructively).

## Timing protocol

- Published numbers: single dedicated HPC node, exclusive allocation, `hardware_tag`
  recorded. Local (≤4-thread) runs are smoke only.
- Setup components timed separately: assembly, factorization, tree build,
  preconditioner (ILU pattern+factor split per 033's `pattern_time` bracket).
- Threading modes: (1 Julia thread + `BLAS.set_num_threads(1)`) and (full threads +
  recorded consistent BLAS count). Harness asserts and logs both at startup; never
  mix modes in one comparison.
- Cost-ceiling drop-outs allowed (e.g. `krylov_gmres` single-thread at the top rung)
  — logged in the ledger, never silently.

## Memory

- Per config: `solver_state_bytes(solver)` (`benchmark/common.jl` — summarysize minus
  every referenced `AbstractBody`), plus `@allocated` during one solve.

## FMM-knob and FGS tuning

- No campaign-wide knob freeze: FMM knobs tuned per rung per solver family, objective
  = per-solve wall time subject to certified BC ≤ target; setup timed separately.
- 021 tuned values are **seeds (hypotheses), never measurements** for this case; the
  seed is recorded in the `notes` column of every row. FGS uses the τ=1e-6-tuned
  config for the solver role (021 ruling 2026-08-17 rationale: coarser-τ configs can
  carry an apply-accuracy BC floor above 1e-6).

## CSV schema

- Static benchmarks: 021 `runs.csv` schema —
  `run_id, phase, solver_config, mesh_file, n_panels, threading_mode, julia_threads,
  blas_threads, t_assembly, t_factorize, t_tree, t_precond, t_rhs, t_solve_min,
  k_reps, iterations, rms_residual, max_residual, mem_state_bytes,
  alloc_solve_bytes, commit, fm_commit, vpm_commit, julia_version, hardware_tag,
  filament_reg, worktrees, solver_settings, backend_settings, notes` — validated by
  `validate_runs_csv`; banner written to `banner.txt` beside every CSV.
- Unsteady (Phase 3): bespoke per-step schema per 021 Phase 3 —
  `niter_first, nsolves, solved, warmstart, warmstart_order, restart_step,
  skip_steps`, the `bcerr_*`/`t_bcerr`/`t_step_net` guard columns, and the identity
  columns (CL, CM, n_wake_rows) replacing 021's (CT, n_particles). The driver
  refuses to append under a mismatched header.
- Per-iteration history: `history_<run_id>.csv` — `iter, t_wall, residual_internal,
  residual_true?` (internal metric labeled as such).

## Run judging

Judge every run from its CSVs, not stdout; verify knobs from the log banner.
