# Phase 4 Log — Green Gate

Plan: [Phase 4 — Green Gate](../../plans/dji_convergence_20260722/phase_04_green_gate.md)

Dashboard: [DJI convergence progress](../../plans/dji_convergence_20260722/dji_convergence_20260722.md)

## Current snapshot

Status: **Blocked on Phase 3 approval**

- Green implementation already exists in `src/FLOWPanel_formulation.jl`.
- `test/formulation_test.jl` is an untracked, substantial validation script;
  preserve it.
- No rotor-specific frozen-wake accuracy gate has been run.
- The current mechanical tests report discrete reconstructed-versus-direct
  `q` error rather than gating it, so physical eligibility remains undecided.

Phase 3 handoff: _pending_

## Established implementation facts

`GreenReconstruction` solves:

```text
(I-B) q = S*sigma
G*muE = -S*sigma0 - q
```

Supported Green controls:

- `gauge=:area_mean`, or dense-only `gauge=:lsq`;
- `recompute_interval >= 1`;
- `green_solver=nothing` for dense bordered LU;
- `KrylovSolver` for matrix-free bordered solve;
- `FGSSolver` for relaxed Picard iteration.

`TraceCorrected(estimator=:green)` uses the reconstructed `C*q` through an
affine Kutta correction. Green routes use sampled wake velocity and support
particle/rolled-up free wakes. `DirectWakePotential` is the scalar-potential
oracle and requires a compatible finite `PanelWake`; particle and mixed wakes
cannot supply its complete scalar potential.

Existing `test/formulation_test.jl` covers:

- Green residual and gauge checks;
- area-mean versus least-squares Kutta-trace agreement;
- dense, Krylov, and FGS route comparisons;
- end-to-end Green versus trace-corrected identities;
- direct-potential oracle checks;
- ConstantDoublet versus VortexRing potential/circulation checks.

## Working records

### Regression tests

| Command | Result | Failure/root cause |
|---|---|---|

### Frozen-wake and solver-route results

| Route/formulation | Residual | Gauge defect | q error | Cq error | Circulation error | Status |
|---|---:|---:|---:|---:|---:|---|

### Classification and next-phase handoff

Green classification: _pending_

Accepted production route/settings: _pending_

## Dated entries

Append exact commands, changed files and why, compact result artifacts,
failures/root causes, and classification decisions here. Keep the snapshot
above current after every material batch.
