# Phase 1 — Circulation Audit

## Purpose and required context

Determine whether the re-fit DJI 9443 geometry changes bound circulation before
committing to free-wake HPC controls.

Read the repository `AGENTS.md`, both routed `CLAUDE.md` files, and:

- `agent_policies/WORKFLOW.md`
- `agent_policies/MONITORS.md`
- `agent_policies/TESTING.md` when selecting validation

Also read the dashboard and the top snapshot in the
[Phase 1 log](../../logs/dji_convergence_20260722/phase_01_circulation_audit.md).
Do not read later phase documents.

## Matched cases

Run locally without HPC or a free wake. Use one steady rotating solve at
5400 RPM, `R=0.119` m, negligible freestream, normal Kutta/shedding
construction, and otherwise identical numerical settings.

Primary capped/Dirichlet airfoil-only comparisons:

- old `examples/data/dji9443_new_40_40.msh` versus new
  `examples/data/dji9443_20260722_40_41_capped.msh`;
- old `examples/data/dji9443_56_57.msh` versus new
  `examples/data/dji9443_20260722_57_57_capped.msh`.

Record the new uncapped/Neumann meshes as secondary comparisons. Do not call
them airfoil-only comparisons because topology and formulation also change.
The new mesh coordinates are dimensional at about 0.0597 m radius; rescale all
meshes consistently to `R=0.119` m.

Use these zero-based seed values exactly as Julia expressions. Supply explicit
`end_node` values, and trace only after the no-shedding constructor rewinds the
body:

```julia
# dji9443_20260722_40_41_capped.msh
te_indices_1 = [1618, 1578, 0] .+ 1
te_indices_2 = [3332, 3292, 1714] .+ 1

# dji9443_20260722_40_41_uncapped.msh
te_indices_1 = [3161, 3121, 1600] .+ 1
te_indices_2 = [1561, 1521, 0] .+ 1

# dji9443_20260722_57_57_capped.msh
te_indices_1 = [6572, 6516, 3354] .+ 1
te_indices_2 = [3218, 3162, 0] .+ 1

# dji9443_20260722_57_57_uncapped.msh
te_indices_1 = [6329, 6273, 3192] .+ 1
te_indices_2 = [3137, 3081, 0] .+ 1
```

## Measurements and decision

Use `BoundCirculationMonitor` and direct trailing-edge strength jumps to save:

- `Gamma_TE(r/R)` and slice-derived circulation for both blades;
- blade-to-blade symmetry;
- new/old ratios interpolated to fixed radial bins;
- `integral(Gamma dr)` and `integral(Omega*r*Gamma dr)`;
- outboard circulation over `r/R >= 0.7`.

One `steady!` solve per mesh is primary. Add a two-step no-shed consistency run
only if monitor and direct extraction differ by more than `1e-6` relative or
initialization changes the answer.

A consistent increase greater than 1% in integrated and outboard circulation
at both resolutions supports the re-fit direction, but does not validate
thrust. A re-fit change of at least 1% in either direction schedules a matched
old/new 40-series cold HPC control in Phase 5. A sub-1% change is negligible
for that decision.

## Deliverables and exit gate

- Setup/topology validation for every old/new comparison mesh.
- Fixed-bin circulation CSVs and a compact comparison summary.
- A reviewed conclusion covering extraction consistency, symmetry, integrated
  and outboard changes, and whether the old 40-series HPC control is required.
- Updated Phase 1 log and dashboard checkboxes.

Stop after reporting the results. Do not begin Phase 2 until Ryan explicitly
approves it.

