# Phase 1 Log — Circulation Audit

Plan: [Phase 1 — Circulation Audit](../../plans/dji_convergence_20260722/phase_01_circulation_audit.md)

Dashboard: [DJI convergence progress](../../plans/dji_convergence_20260722/dji_convergence_20260722.md)

## Current snapshot

Status: **Phase 1 complete — awaiting approval**

- The dedicated local driver is
  `examples/dji9443_circulation_audit.jl`.
- The worktree was already dirty. Preserve all unrelated changes and untracked
  files; inspect `git status --short` before editing.
- All six meshes passed setup-only validation with two symmetric, monotone
  shedding chains and the expected node/panel/section counts.
- Old meshes use explicit root endpoints without the new-mesh `|r/R| >= 0.1`
  boxes; their established root endpoints lie inside that cutoff. New meshes
  use the side-specific boxes plus explicit endpoints.
- All six steady cases are complete. Monitor and independent TE jumps agree
  exactly; blade symmetry errors are below `6e-9` relative.
- Re-fit changes on the fixed grid are `+0.534%` integrated / `+4.078%`
  outboard at 40-series and `+0.517%` integrated / `+4.107%` outboard at
  57-series.
- The outboard change exceeds 1% at both resolutions, so the old 40-series
  cold HPC control is required in Phase 5. The stricter Phase 1 direction gate
  is not met because the full integrated increases are below 1%.
- Secondary uncapped/Neumann changes relative to the new capped cases are
  `+6.519%` integrated at 40-series and `+4.853%` at 57-series. These are
  topology/formulation sensitivities, not airfoil-only comparisons.
- Directly against the old meshes, new uncapped/Neumann integrated circulation
  is `+7.089%` at 40-series and `+5.394%` at 57-series.
- Refining the new capped/Dirichlet mesh changes integrated circulation by
  only `+0.255%`; refining the new uncapped/Neumann mesh changes it by
  `-1.314%`, just outside the study's 1% numerical-convergence criterion.
- Ryan's interpretation is that the repeatable approximately 5–7%
  uncapped/Neumann excess raises suspicion that the capped/Dirichlet
  formulation may be wrong. This is a follow-up hypothesis, not established
  by the confounded topology/formulation comparison.
- Reviewed artifacts are under
  `data/dji_convergence_20260722/phase_01_circulation_audit/`.
- Stop here. Do not begin Phase 2 without Ryan's explicit approval.
- Required Phase 1 decision: whether the absolute re-fit circulation change
  reaches 1% and therefore schedules an old 40-series HPC control in Phase 5.

## Established inputs

Verified new-mesh traces use `normal_jump_tol=0.2`,
`max_turn_angle=pi/3`, and a side-specific `|r/R| >= 0.1` box:

| Mesh | Nodes | Panels | One-based blade seed/end triples | Edges/blade |
|---|---:|---:|---|---:|
| `dji9443_20260722_40_41_capped.msh` | 3428 | 6848 | `(1619,1579,1)`, `(3333,3293,1715)` | 39 |
| `dji9443_20260722_40_41_uncapped.msh` | 3200 | 6240 | `(1562,1522,1)`, `(3162,3122,1601)` | 39 |
| `dji9443_20260722_57_57_capped.msh` | 6708 | 13408 | `(3219,3163,1)`, `(6573,6517,3355)` | 56 |
| `dji9443_20260722_57_57_uncapped.msh` | 6384 | 12544 | `(3138,3082,1)`, `(6330,6274,3193)` | 56 |

Without explicit `end_node`, capped traces continued around the root cap and
returned after 79 or 112 edges. That shedding is invalid.

The new coordinates are dimensional at approximately 0.0597 m radius and must
be rescaled consistently to `R=0.119` m.

Old comparison meshes:

- `examples/data/dji9443_new_40_40.msh`: 3646 nodes, 7288 panels.
- `examples/data/dji9443_56_57.msh`: 6956 nodes, 13908 panels.

Established steady 5400 RPM results from
`data/rhc_sweep_baseline/convergence_history.csv`:

| Mesh | CT |
|---|---:|
| old 40-series | `0.05048051147416609` |
| old 56/57-series | `0.051489958562886455` |
| old 80/81-series | `0.05060653180935334` |

These `CT` values are context only; Phase 1's primary result is circulation.

## Working records

### Files and commands

| Date | Command or file | Purpose/result |
|---|---|---|
| 2026-07-22 | `examples/dji9443_circulation_audit.jl` | Added reproducible setup, solve-batch, and analysis driver. |
| 2026-07-22 | `PHASE1_MODE=setup julia --project examples/dji9443_circulation_audit.jl` | All six meshes passed topology validation; artifact: `data/dji_convergence_20260722/phase_01_circulation_audit/topology_validation.csv`. |
| 2026-07-22 | `include("test/test_helpers.jl"); include("test/runtests_unit_liftingbody.jl")` | 40/40 tests passed. |
| 2026-07-22 | `include("test/runtests_unit_postprocess.jl")` | 405/405 tests passed. |
| 2026-07-22 | `PHASE1_MODE=primary julia --project examples/dji9443_circulation_audit.jl` | Four primary steady cases complete; exact monitor/direct TE agreement. |
| 2026-07-22 | `src/FLOWPanel_simulate_monitors.jl` | Corrected `BoundCirculationMonitor` station midpoint lookup from shedding-local slots to global node IDs; final postprocess suite 406/406 passed. |
| 2026-07-22 | `PHASE1_MODE=secondary julia --project examples/dji9443_circulation_audit.jl` | Two secondary uncapped/Neumann cases complete; exact monitor/direct agreement. |
| 2026-07-22 | `PHASE1_MODE=analyze julia --project examples/dji9443_circulation_audit.jl` | Wrote fixed-bin CSVs, metrics, comparisons, and reviewed Markdown report. |

### Case results

| Case | Extraction agreement | Symmetry | Integrated change | Outboard change | Decision |
|---|---:|---:|---:|---:|---|
| old40 → new40c | `0.0` relative | `<8e-11` | `+0.534%` | `+4.078%` | Phase 5 old-40 control required |
| old57 → new57c | `0.0` relative | `<6e-9` | `+0.517%` | `+4.107%` | Confirms resolution-consistent outboard increase |
| new40c → new40u | `0.0` relative | `<3e-14` | `+6.519%` | `+6.395%` | Topology/formulation sensitivity only |
| new57c → new57u | `0.0` relative | `<6e-14` | `+4.853%` | `+4.454%` | Topology/formulation sensitivity only |

### Decisions and next-phase handoff

Fixed comparison bins are `|r/R| = 0.125:0.025:0.975`.

- The re-fit increases integrated circulation by about `0.52%` and outboard
  circulation by about `4.1%` at both resolutions.
- The Phase 1 support threshold (>1% for both integrated and outboard metrics)
  is not met.
- The absolute outboard change exceeds 1%, scheduling the matched old
  40-series cold HPC control in Phase 5.
- Uncapped/Neumann circulation is materially higher than capped/Dirichlet
  circulation, so those secondary cases cannot isolate the airfoil re-fit.
- Mesh refinement has little effect on the new capped result (`+0.255%`) but
  changes the new uncapped result by `-1.314%`.
- Relative to the resolution-matched old meshes, the new uncapped results are
  `+7.089%` (40-series) and `+5.394%` (57-series).
- Ryan suspects the capped/Dirichlet formulation is responsible for the
  approximately 5–7% formulation gap. Keep this as an explicit diagnostic
  hypothesis for the formulation audit; current evidence does not separate
  cap topology from boundary-condition formulation.
- Phase 1 is complete. No Phase 2 work is authorized.

## Dated entries

Append material actions, failures/root causes, quantitative results, and
artifact paths here. Keep the snapshot above current after every batch.

### 2026-07-22 — Setup/topology batch

- Confirmed the six case meshes and formulations specified by the Phase 1
  plan.
- A validation attempt applying the new-mesh radial box to the old 40-series
  mesh stopped at node 214 before its explicit endpoint 46. Root cause: that
  established endpoint is inside `|r/R| = 0.1`. Explicit unboxed old-mesh
  traces and boxed new-mesh traces all pass.
- Validation ranges span approximately `0.0733–0.9853` (old 40),
  `0.1116–0.9856` (new capped 40), `0.0549–0.9890` (old 57), and
  `0.1081–0.9891` (new capped 57), supporting the non-extrapolated fixed grid.
- The first narrow lifting-body test invocation omitted `test_helpers.jl` and
  errored on undefined fixtures. Rerunning with the suite's required helper
  produced 40/40 passes; this was a command/setup issue, not a code failure.

### 2026-07-22 — Primary steady-solve batch

- Initial inspection found constant monitor `r/R` and zero slice circulation.
  Root cause of the station error: shedding rows 2/3 are panel-local node
  slots, but `BoundCirculationMonitor` treated them as global node IDs.
  Correcting the lookup restored the expected radial ranges; a multi-section
  regression was added and the postprocess suite passes 406/406.
- The slice estimator remains zero on these station-aligned rotor meshes
  (`100%` TE/slice discrepancy). Its strict edge-crossing construction
  degenerates when section planes coincide with spanwise mesh vertices.
  A proposed half-open crossing rule failed the existing analytic oracle and
  was not retained. Phase 1 decisions therefore use the exact-agreement TE
  monitor/direct observables; the slice limitation is reported explicitly.
- Primary fixed-grid metrics:
  - old40: `∫Γdr=0.01763757156`, outboard `0.00596576376`;
  - new40c: `∫Γdr=0.01773180521`, outboard `0.00620903396`;
  - old57: `∫Γdr=0.01768571512`, outboard `0.00599212859`;
  - new57c: `∫Γdr=0.01777706344`, outboard `0.00623819565`.

### 2026-07-22 — Secondary and final analysis batch

- new40u: `∫Γdr=0.01888781531`, outboard `0.00660612919`;
  relative to new40c this is `+6.519%` integrated and `+6.395%`
  outboard.
- new57u: `∫Γdr=0.01863971417`, outboard `0.00651603397`;
  relative to new57c this is `+4.853%` integrated and `+4.454%`
  outboard.
- Reviewed raw extraction values, fixed-bin interpolation, trapezoidal
  integrals, comparison percentages, and the generated decision text for
  consistency. The same conclusion holds at both mesh resolutions.
- Final report:
  `data/dji_convergence_20260722/phase_01_circulation_audit/phase_01_report.md`.

### 2026-07-22 — Ryan review observations

- New-mesh refinement comparison on the same fixed radial grid:
  - capped/Dirichlet, new40c → new57c: `+0.255%` integrated;
  - uncapped/Neumann, new40u → new57u: `-1.314%` integrated.
- Direct old-to-new-uncapped comparisons:
  - old40 → new40u: `+7.089%` integrated;
  - old57 → new57u: `+5.394%` integrated.
- Ryan noted that the approximately 5–7% higher uncapped/Neumann circulation
  makes the capped/Dirichlet formulation suspect. The observation is
  resolution-repeatable, but the present cases change both cap topology and
  boundary-condition formulation, so they do not yet identify which change
  causes the discrepancy.
