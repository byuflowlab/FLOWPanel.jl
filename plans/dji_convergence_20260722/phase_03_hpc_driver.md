# Phase 3 — HPC Driver

## Purpose and required context

Build and locally verify the parameterized hover driver and guarded Slurm
interface needed by the later adaptive study.

Read the repository instructions, the dashboard, the top snapshot in the
[Phase 3 log](../../logs/dji_convergence_20260722/phase_03_hpc_driver.md), and:

- `agent_policies/WORKFLOW.md`
- `agent_policies/TESTING.md`
- `agent_policies/MONITORS.md`
- `agent_policies/HPC.md`

Do not begin until Phase 2b approval and its formulation/topology attribution
handoff appear in the Phase 3 log, unless Ryan explicitly cancels Phase 2b.

## Mesh and body contract

Support these four new meshes:

| Mesh | Nodes | Panels | Formulation | One-based blade seed/end triples | Edges/blade |
|---|---:|---:|---|---|---:|
| `dji9443_20260722_40_41_capped.msh` | 3428 | 6848 | Dirichlet source+vortex | `(1619,1579,1)`, `(3333,3293,1715)` | 39 |
| `dji9443_20260722_40_41_uncapped.msh` | 3200 | 6240 | Neumann vortex only | `(1562,1522,1)`, `(3162,3122,1601)` | 39 |
| `dji9443_20260722_57_57_capped.msh` | 6708 | 13408 | Dirichlet source+vortex | `(3219,3163,1)`, `(6573,6517,3355)` | 56 |
| `dji9443_20260722_57_57_uncapped.msh` | 6384 | 12544 | Neumann vortex only | `(3138,3082,1)`, `(6330,6274,3193)` | 56 |

- Capped: `Union{ConstantSource,VortexRing}`, `watertight=true`.
- Uncapped: `VortexRing`, `watertight=false`.
- Build the no-shedding body first, derive shedding from its constructed
  nodes/cells, then rebuild.
- Always use explicit capped-mesh root endpoints; otherwise a trace can
  continue around the cap.
- Rescale the approximately 0.0597 m mesh coordinates consistently to
  `R=0.119` m.

## Driver interface

Expose independent environment-driven controls for:

- mesh and solve formulation;
- timestep;
- wake `core_size`;
- `kerneloffset_targets`;
- mu-gradient reconstruction;
- Green `gauge`, `recompute_interval`, and solver route.

Keep these controls in the study driver; no public FLOWPanel API change is
expected.

Formulation choices:

- `velocity`: production `VelocityThroughSources`;
- `green`: `GreenReconstruction`;
- `direct`: diagnostic finite-panel-wake `DirectWakePotential`;
- `trace_green`: diagnostic `TraceCorrected(estimator=:green)`.

Green solver routes are dense bordered LU, matrix-free GMRES through
`KrylovSolver`, and relaxed Picard through `FGSSolver`. Green and Green trace
reconstruction consume sampled wake velocity and can operate with particle or
rolled-up wakes. Only `DirectWakePotential` requires a compatible finite
`PanelWake` with scalar-potential sources.

Record:

- Bernoulli and Kutta-Joukowski `CT`;
- bound-circulation histories and one-revolution block statistics;
- formulation diagnostics, Green residuals, and iteration counts;
- maximum resident memory, restart provenance, wall time, and every numerical
  setting.

PressureLaplace is excluded from routine convergence runs; use it only if force
recovery becomes suspect.

## Slurm and deployment

Add a generic launcher with:

- `--ntasks=64`, `--mem=32G`, `--time=12:00:00`;
- `JULIA_NUM_THREADS`, `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`,
  `BLAS_NUM_THREADS`, and `MKL_NUM_THREADS` all set to 64;
- Julia invoked with `--project=. -t 64`;
- explicit output/error paths and `set -euo pipefail`.

Before submission, query Slurm and refuse to submit if three study jobs are
already active. Submission itself belongs to Phase 5; Phase 3 only implements
and verifies the interface.

Deploy only the driver, launcher, four meshes, and directly required changed
source files to `/home/rander39/projects/FLOWPanel.jl`. Verify checksums and
preserve unrelated remote changes. Recheck cluster `PATH`, Julia, and Slurm
rather than assuming the earlier environment remains valid.

## Deliverables and exit gate

- Parameterized single-case driver with compact machine-readable summaries.
- Guarded 64-thread/32-GB Slurm launcher.
- Setup-only topology checks and short capped/uncapped local smoke tests.
- Focused tests selected from repository policy.
- Exact deployment manifest and verified local/remote checksums.
- Updated Phase 3 log and dashboard.

Report implementation and local verification, then stop for explicit Phase 4
approval. Do not submit production jobs in this phase.
