# Suddenly Started Wing / Wagner Validation Handoff

Last updated: 2026-07-21 (backend isolation checkpoint)

## Objective

Validate FLOWPanel's unsteady solver against the Wagner function using the setup
in `/Users/ryan/Dropbox/research/projects/VortexLattice.jl/test/suddenly_started_wing_uvlm.jl`:

- suddenly started wing at `alpha = 5 deg`
- `AR = 100`, `c = Uinf = 1`
- VortexLattice reference mesh: 5 chordwise by 48 spanwise panels
- time samples `t Uinf/c = 0:1/8:7` (the reference file's `1/16` header comment
  is stale; its executable timestep is `1/8`)
- compare `CL(t)/CLsteady` with Jones' approximation to Wagner's function,
  `1 - 0.165 exp(-0.09 t*) - 0.335 exp(-0.6 t*)`

The requested first FLOWPanel configuration is an open-tip, extruded NACA 0012
surface with no endcaps, Neumann boundary conditions, and vortex-ring panels.

## Implemented files

- `examples/suddenly_started_wing.jl`
  - Generates a closed-TE cosine-spaced NACA 0012 contour and extrudes it over
    the span without endcaps.
  - Builds `RigidWakeBody{VortexRing,1,Float64,false}` with `watertight=false`.
  - Correctly derives shedding from the already-constructed/rewound base body.
  - Runs a semi-infinite-wake steady case for `CLsteady` and a free-wake
    unsteady case for `CL(t)`.
  - Writes histories, convergence summaries, VTK output, and Wagner/convergence
    plots.
  - Supports `single`, `coarse`, and resumable `convergence` modes through
    `SSW_*` environment variables.
  - The current worktree additionally exposes `SSW_FMM_ORDER`,
    `SSW_FMM_THETA`, and `SSW_FMM_LEAF` for backend sensitivity diagnostics.
- `test/runtests_example_suddenly_started_wing.jl`
  - Covers NACA geometry, open tips, Neumann vortex-ring formulation, TE
    shedding, timestep/inflow, Wagner values, wake allocation, monitor
    contracts, and a one-step unsteady smoke test.
- `test/runtests.jl`
  - Includes the new example test near the other unsteady-wing tests.
- `scripts/suddenly_started_wing.slurm.sh`
  - Single-node, 16-thread, 64-GB, 24-hour resumable convergence launcher.
  - The current worktree adds `module load julia`; confirm that this matches the
    target cluster before committing it.

The example and its dedicated test are already tracked and currently clean in
Git. The `test/runtests.jl` include and launcher's `module load julia` are
uncommitted. There are many unrelated dirty/untracked files in the repository;
do not clean, revert, or include them accidentally.

## Important startup treatment

FLOWPanel solves at `t=0`, unlike the VortexLattice history which reports one
value per completed interval. The implemented startup does the following:

1. Uses genuinely zero freestream for the `t=0` body solve so unsteady
   Bernoulli pressure history is initialized consistently.
2. A final monitor initializes the unattached wake row's velocity to the full
   freestream at step zero.
3. Wake propagation then creates a nondegenerate first shed panel; all later
   body solves use the full suddenly applied freestream.
4. The reported FLOWPanel curve drops the `t=0` sample and uses `dt:dt:t_end`,
   matching VortexLattice.

Without step 2, the first wake panel is degenerate. Monitor tuple order matters:
`(PressureBernoulli, ForceMonitor, initialize_wake_convection!)`.

## Backend and solver settings

The FMM backend is explicitly:

```julia
pnl.FastMultipoleBackend(
    expansion_order=10,
    multipole_acceptance=0.4,
    leaf_size=100,
)
```

For body meshes at or below `SSW_BACKSLASH_MAX_PANELS`, the body solve is dense
`Backslash`/LU. The example default cutoff is 10,000 panels; the Slurm launcher
raises it to 20,000. Above the cutoff it uses matrix-free FMM-backed GMRES with
`atol=1e-9`, `rtol=1e-9`, `itmax=1000`, and no preconditioner.

Consequently, every completed case so far (maximum 3,840 body panels) used a
dense direct body solve. FMM was still used for wake/body induced interactions,
wake rollup, and Bernoulli scalar-potential evaluations. Do not describe these
completed runs as FMM linear solves.

## Validation already performed

The focused example tests passed: 51/51 after adding coverage for configurable
FMM order, acceptance threshold, and leaf size.

The relevant lifting-body, wake, postprocess, and simulation test suites also
passed in the previous session (reported totals included 40 lifting-body, 82
wake, and 405 postprocess assertions). The launcher passed `bash -n`, and the
task's edits passed `git diff --check`.

Useful focused command:

```bash
julia --project=. -e 'include("test/runtests_example_suddenly_started_wing.jl")'
```

## Current numerical results

Canonical completed histories have backend-qualified directory names under
`data/suddenly_started_wing/`. The current spatial series uses requested
`n_airfoil=21` (20 actual contour nodes/40 triangular panels per span strip),
`dt*=0.125`, and the FMM interaction backend:

| `n_span` | panels | Wagner RMS | Wagner max | `CLsteady` |
|---:|---:|---:|---:|---:|
| 12 | 480 | 0.31659 | 1.28717 | 136618.46 |
| 24 | 960 | 1.18552 | 5.29389 | 3270.90 |
| 48 | 1920 | 0.23178 | 0.53911 | 113.116 |
| 96 | 3840 | 0.20367 | 0.57085 | 6.98186 |

The `48 -> 96` normalized-curve change is:

- maximum absolute: 0.24457
- maximum relative: 0.24667
- relative L2: 0.071762

This is not converged to the example's 2% criterion. The very large and highly
mesh-dependent steady lift coefficients are also physically suspect (a 2-D
thin-airfoil estimate at 5 degrees is roughly 0.55). The likely first issue to
investigate is the extremely stretched surface panels caused by combining
`AR=100` with only 12--96 span strips. Do not infer validation success from the
normalized curves alone.

A coarse direct-interaction cross-check also differs materially from FMM:
maximum absolute curve difference 1.5952, maximum relative 0.75374, relative
L2 0.27197. Since the body solve is dense in both cases, this discrepancy is in
the induced-interaction/pressure/wake backend path and deserves investigation.

## Backend isolation result

A shortened `n_span=12`, `t_end*=2.25`, no-VTK diagnostic separated the
interaction and pressure paths:

- Direct and FMM agree closely initially, but at `t*=2.25` Direct gives
  `CL/CLsteady=0.96385360` and standard FMM gives `1.02195642` (absolute
  difference `0.05810283`, relative-L2 history difference `0.0152128`).
- `multipole_acceptance=0` is bit-for-bit identical to standard FMM, so the
  multipole far-field approximation is not the cause in this diagnostic.
- `leaf_size=10000`, which prevents tree subdivision for the coarse system,
  matches Direct bit-for-bit. The difference therefore enters with FMM tree
  partition/direct-list evaluation and is then amplified by wake evolution.
- FMM dynamics with Direct `PressureBernoulli` reproduces the FMM curve, while
  Direct dynamics with FMM pressure reproduces Direct. Pressure-potential
  evaluation is not the source.
- Splitting simulation backends identifies bound-body influence on wake nodes
  (`backend_system`) as dominant: Direct `backend_wake` plus FMM
  `backend_system` ends at `1.02168490`, while FMM `backend_wake` plus Direct
  `backend_system` ends at `0.96260637`.

The robust conclusion is not yet an FMM kernel defect: the extremely
underresolved free wake amplifies tiny tree/summation perturbations into
different trajectories. Span refinement and a targeted `backend_system` leaf
sensitivity are needed before deciding whether a backend implementation issue
remains on credible geometry.

Generated summaries/plots (ignored by `data/*`) are:

- `data/suddenly_started_wing/convergence.csv`
- `data/suddenly_started_wing/wagner_comparison.png`
- `data/suddenly_started_wing/convergence.png`

There are older preliminary directories without `_direct` or `_fmm` suffixes.
Treat the backend-qualified directories as canonical; do not let the old paths
confuse resume checks or comparisons.

## Next work

1. Run a targeted refined-geometry `backend_system`/leaf-size sensitivity
   before accepting the canonical FMM curve. The coarse diagnostic above shows
   that merely setting `multipole_acceptance=0` is insufficient; tree
   subdivision, not M2L acceptance, controls the observed trajectory change.
2. Run/inspect the next span refinement (`n_span=192`, 7,680 panels) and keep
   refining until the curve stabilizes or the formulation problem is clear.
   This is a long run and belongs on HPC under `agent_policies/HPC.md`.
3. Once span convergence is credible, refine the airfoil contour (`21, 41, 81`)
   and then halve the timestep. The driver already implements this sequence.
4. Reassess whether the thick closed NACA body should be compared directly to
   the VortexLattice flat/camber-surface model or only to Wagner. Document that
   modeling distinction in the final conclusions.
5. Only after physical and backend convergence, update the plots and report
   whether FLOWPanel reproduces the Wagner response.

The agent may prepare and syntax-check the Slurm script but must not submit it.
The user submits from repository root with:

```bash
sbatch scripts/suddenly_started_wing.slurm.sh
```

The launcher defaults to `SSW_RESUME=true`, so it reuses canonical completed
histories and continues from the first missing case. Before submission, verify
the site's Julia module name and whether the 20,000-panel dense cutoff is
appropriate for the available memory.

## Restart checklist

1. Read `CLAUDE.md`, then `agent_policies/WORKFLOW.md`, `TESTING.md`,
   `MONITORS.md`, and `HPC.md`.
2. Run `git status --short`; preserve unrelated changes.
3. Read this file and `examples/suddenly_started_wing.jl`.
4. Inspect `data/suddenly_started_wing/convergence.csv` and the two plots.
5. Confirm focused tests still pass before changing startup, pressure recovery,
   force normalization, or backend settings.
