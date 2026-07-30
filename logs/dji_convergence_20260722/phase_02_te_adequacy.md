# Phase 2 Log — TE Adequacy

Plan: [Phase 2 — TE Adequacy](../../plans/dji_convergence_20260722/phase_02_te_adequacy.md)

Dashboard: [DJI convergence progress](../../plans/dji_convergence_20260722/dji_convergence_20260722.md)

## Current snapshot

Status: **Phase 2 complete — awaiting Phase 3 approval**

- Ryan approved Phase 1 completion and authorized Phase 2 work.
- Added `examples/dji9443_te_adequacy.jl` for Phase 2 setup, smoke,
  full diagnostic, and analysis batches.
- Real DJI setup validation passes for `new40c` and `new57c`; artifact:
  `data/dji_convergence_20260722/phase_02_te_adequacy/topology_validation.csv`.
- Synthetic smoke validation and an independent small-matrix low-rank update
  check pass. Smoke artifacts validate the driver only; they are not study
  evidence.
- Full `new40c` and `new57c` diagnostics completed with geometry/influence,
  adjoint sensitivity, perturbation response, off-collocation, and
  `kerneloffset_panel` sweep artifacts.
- Maximum integrated/outboard perturbation responses were `0.084%`/`0.232%`
  for `new40c` and `0.084%`/`0.214%` for `new57c`.
- `kerneloffset_panel` decade sweep changed integrated and outboard
  circulation by `0.0%` at reported precision for both cases.
- Offset-separated off-collocation TE residuals were comparable to the sampled
  near-TE population. Interior median potential residual decreased from
  `1.142e-4` to `8.228e-5`; exterior median changed from `8.689e-2` to
  `8.638e-2`.
- Trigger decision: no 1% perturbation trigger, no material sensitivity-growth
  trigger, no off-collocation trigger, no panel-kernel-offset trigger.
- Sharp capped/Dirichlet TE adequacy is accepted for continuing to Phase 3.
  Conditional thickness/local-refinement is not triggered by this batch.
- Keep `kerneloffset_panel`, later `kerneloffset_targets`, and wake core as
  distinct variables in all records and conclusions.

Phase 1 handoff: capped/Dirichlet new-mesh refinement changed integrated
circulation by only `+0.255%`, but uncapped/Neumann circulation was
approximately `5–7%` higher than capped/Dirichlet. Ryan suspects the
capped/Dirichlet formulation may be responsible; Phase 1 evidence does not
separate cap topology from boundary-condition formulation. The old 40-series
cold HPC control remains scheduled for Phase 5 because the outboard re-fit
change exceeded 1%.

## Working records

### Diagnostic implementation and validation

| Date | Command or file | Purpose/result |
|---|---|---|
| 2026-07-23 | `examples/dji9443_te_adequacy.jl` | Added Phase 2 diagnostic driver with setup, smoke, diagnostic, and analyze modes. |
| 2026-07-23 | `PHASE2_MODE=smoke julia --project examples/dji9443_te_adequacy.jl` | Passed synthetic paired-TE writer validation; wrote smoke CSV artifacts only. |
| 2026-07-23 | `PHASE2_MODE=setup julia --project examples/dji9443_te_adequacy.jl` | Passed real DJI topology validation for `new40c` and `new57c`; wrote `topology_validation.csv`. |
| 2026-07-23 | `PHASE2_MODE=diagnostic PHASE2_CASE=new40c julia --project examples/dji9443_te_adequacy.jl` | Interrupted during dense local solve/assembly after several minutes; partial `geometry_influence_new40c.csv` exists and must be regenerated before use. |
| 2026-07-23 | low-rank update explicit check | Small-matrix Sherman-Morrison-Woodbury update matched explicit perturbed solve exactly (`max abs diff = 0.0`). |
| 2026-07-23 | `PHASE2_MODE=diagnostic PHASE2_CASE=new40c julia --project examples/dji9443_te_adequacy.jl` | Completed full diagnostics after optimization; wrote all `new40c` artifacts. |
| 2026-07-23 | `PHASE2_MODE=diagnostic PHASE2_CASE=new57c julia --project examples/dji9443_te_adequacy.jl` | Initial run interrupted during slow direct RHS in kernel-offset sweep after main CSVs; driver optimized to use FMM RHS and avoid `B` assembly in offset sweep. |
| 2026-07-23 | `PHASE2_MODE=diagnostic PHASE2_CASE=new40c julia --project examples/dji9443_te_adequacy.jl` | Reran `new40c` after FMM-RHS change for consistency; completed all artifacts. |
| 2026-07-23 | `PHASE2_MODE=diagnostic PHASE2_CASE=new57c julia --project examples/dji9443_te_adequacy.jl` | Completed all `new57c` artifacts; baseline dense solve block reported `135.0 s`. |
| 2026-07-23 | `PHASE2_MODE=analyze julia --project examples/dji9443_te_adequacy.jl` | Wrote `phase_02_summary.csv` and final `phase_02_report.md`. |

### Per-edge summary artifacts

| Mesh | Artifact | Geometry/influence | Adjoint | Off-collocation | Offset |
|---|---|---|---|---|---|
| synthetic smoke | `data/dji_convergence_20260722/phase_02_te_adequacy/*_smoke.csv` | Driver validation only | Driver validation only | Driver validation only | Not run |
| `new40c` | `data/dji_convergence_20260722/phase_02_te_adequacy/*_new40c.csv` | Complete | Complete | Complete | Complete |
| `new57c` | `data/dji_convergence_20260722/phase_02_te_adequacy/*_new57c.csv` | Complete | Complete | Complete | Complete |

### Trigger assessment

| Trigger | Result | Thickness study required? |
|---|---|---|
| 1% influence perturbation changes circulation by at least 1% | No: max integrated/outboard response `0.084%`/`0.232%` | No |
| Sensitivity grows from 40 to 57 series | No material growth: max sensitivity `17.58` to `19.07`; perturbation response slightly decreased | No |
| TE residuals are excessive or do not refine | No: off-surface TE medians comparable to sampled population; interior decreased, exterior essentially unchanged | No |
| Panel kernel offset changes circulation by at least 1% | No: `0.0%` at reported precision for integrated and outboard metrics | No |

### Adequacy decision and next-phase handoff

Sharp capped/Dirichlet TE adequacy is accepted for proceeding to Phase 3.
Conditional thickness/local-refinement is not triggered. Continue to preserve
the Phase 1 hypothesis as context: capped/Dirichlet versus uncapped/Neumann
still differs by approximately `5–7%`, but Phase 2 did not identify sharp-TE
operator fragility as the cause.

## Dated entries

Append exact commands, changed files, artifact paths, failures/root causes,
quantitative results, and decisions here. Keep the snapshot above current after
every material batch.

### 2026-07-23 — Driver implementation and local validation

- Implemented `examples/dji9443_te_adequacy.jl`.
- The driver reuses the Phase 1 safe construction pattern: load mesh, build a
  no-shedding capped/Dirichlet `RigidWakeBody`, trace shedding from the
  constructed body's rewound cells, then rebuild with shedding.
- Real setup validation passed for `new40c` and `new57c`.
- Synthetic smoke validation passed for the geometry/influence, adjoint,
  perturbation, and off-collocation artifact writers.
- A full `new40c` local diagnostic attempt was manually interrupted after the
  runtime became unsuitable for quick implementation validation. Root cause of
  the interruption was runtime/scale, not a reported code error; the stack was
  in BLAS during dense solve/assembly work.

### 2026-07-23 — Full diagnostics and adequacy decision

- Replaced repeated dense perturbation re-solves with a
  Sherman-Morrison-Woodbury low-rank update around one LU factorization of
  `G`. An independent small-matrix check matched explicit `Gp \ rhs` exactly.
- Changed source-RHS assembly to use the same `FastMultipoleBackend` pattern as
  Phase 1 steady solves; the initial direct-RHS path was too slow in the
  `new57c` kernel-offset sweep.
- Completed full diagnostics for `new40c` and `new57c`.
- Summary metrics:
  - `new40c`: max integrated/outboard perturbation `0.084%`/`0.232%`;
    max panel-offset integrated/outboard change `0.0%`/`0.0%`.
  - `new57c`: max integrated/outboard perturbation `0.084%`/`0.214%`;
    max panel-offset integrated/outboard change `0.0%`/`0.0%`.
- Offset-separated off-collocation check:
  - interior TE median potential residual `1.142e-4` (`new40c`) to
    `8.228e-5` (`new57c`);
  - exterior TE median potential residual `8.689e-2` to `8.638e-2`;
  - TE medians were comparable to all sampled final-row/near-TE medians at
    both resolutions.
- Conclusion after reviewing methods, results, and trigger logic: the sharp
  capped/Dirichlet TE is numerically adequate for continuing. No conditional
  thickness/local-refinement study is required by this batch.
