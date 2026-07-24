# Phase 4 — Green Gate

## Purpose and required context

Validate Green reconstruction mechanically and physically before using it as a
production sensitivity in the convergence study.

Read the repository instructions, the dashboard, the top snapshot in the
[Phase 4 log](../../logs/dji_convergence_20260722/phase_04_green_gate.md), and:

- `agent_policies/WORKFLOW.md`
- `agent_policies/TESTING.md`

Do not begin until Phase 3 approval and the verified driver interface appear in
the Phase 4 log.

## Regression and compatibility gate

Review implementation assumptions and run:

- `test/formulation_test.jl`;
- solver, lifting-body, simulate, postprocess, and warm-start suites selected
  according to `agent_policies/TESTING.md`.

Confirm:

- Green reconstruction is restricted to a single capped Dirichlet
  source/doublet `RigidWakeBody` with `semiinfinite_wake=false` and a
  `Backslash` body solver;
- it reconstructs the potential trace from sampled wake velocity and therefore
  supports particle and rolled-up free wakes;
- direct scalar wake potential is used only by `DirectWakePotential` with a
  compatible finite `PanelWake`;
- body potential, source strength, velocity snapshots, and wake-correction
  state are restored correctly;
- failures and solver non-convergence are reported, never silently accepted.

Identify the root cause before modifying formulation code for any failure.

## Rotor frozen-wake oracle

On the new 40-series capped mesh, construct a short finite `PanelWake`. On an
identical frozen body and wake state, compare:

- directly evaluated wake potential `q_f`;
- Green-reconstructed `q` after area-mean gauge matching;
- gauge-invariant Kutta trace `C*q`;
- spanwise circulation from `VelocityThroughSources`,
  `DirectWakePotential`, `GreenReconstruction`, and
  `TraceCorrected(estimator=:green)`.

Preserve enough compact data to reproduce the comparison without retaining a
full unsteady history.

## Solver-route gates

Use dense area-mean bordered LU as the reference. Validate:

- relative Green-system residual at most `1e-8` for dense/direct and `1e-5`
  for matrix-free FMM;
- area-weighted gauge defect at most `1e-8`;
- solver-route `C*q` and circulation agreement at most `1e-3`;
- no non-finite results or unreported Krylov/Picard non-convergence.

Compare matrix-free GMRES and relaxed Picard/FGS with the dense reference.
Prefer GMRES for production if it passes, avoiding a second dense
factorization on the 57-series mesh. FGS remains diagnostic unless it passes
the same gates.

## Accuracy decision

Against direct finite-wake potential, target:

- gauge-matched `q` error no larger than 5%;
- Kutta-trace and resulting circulation error no larger than 1%.

If either target fails, classify Green reconstruction as sensitivity-only, run
at most one short unsteady comparison, and do not describe it as a physically
improved solution.

Use `recompute_interval=1` for accuracy cases. Test interval 2 only after
interval 1 passes and only as a cost/lag sensitivity.

## Deliverables and exit gate

- Focused test results and compatibility audit.
- Frozen-wake direct-versus-reconstructed dataset.
- Dense, GMRES, and FGS residual/gauge/trace/circulation comparison.
- Explicit classification: production-eligible or sensitivity-only.
- Updated Phase 4 log and dashboard.

Review methods, results, and the classification for consistency, report them,
then stop for explicit Phase 5 approval.

