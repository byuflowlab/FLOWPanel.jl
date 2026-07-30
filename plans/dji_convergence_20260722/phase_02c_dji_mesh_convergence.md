# Phase 2c — DJI 9443 Chordwise Mesh Convergence

## Purpose and required context

Verify the Phase 2b attribution on the **actual DJI 9443 rotor mesh**. Phase 2b
attributed the Dirichlet–Neumann bound-circulation gap to **chordwise (n_airfoil)
section under-resolution** on a simple-wing oracle: open+Neumann converges by ~60
chordwise panels/section, capped+Dirichlet (Morino) climbs up to meet it, and the two
agree to 0.11% at 120 panels/section. Extraction (loop-circ = Γ_TE to 3e-5),
solve-formulation/Kutta trace (0.74%), and thinness (thinner → smaller gap) were
excluded. **Ryan approved Phase 2b on 2026-07-24.**

This phase runs a chordwise-refinement ladder on the DJI blade and asks whether the DJI
gap converges like the oracle. It is inserted between Phase 2b and Phase 3.

Read the repository instructions, the dashboard, the top **Current snapshot** of the
[Phase 2c log](../../logs/dji_convergence_20260722/phase_02c_dji_mesh_convergence.md),
and `agent_policies/WORKFLOW.md`, `agent_policies/TESTING.md`.

Do not begin Phase 3 until Ryan explicitly approves the Phase 2c decision.

## Mesh inventory

Ryan supplied a chordwise-refinement ladder of 6 meshes (fixed 30 spanwise sections,
n_airfoil ∈ {81, 97, 121}, triangles, Gmsh `.msh`) in `examples/data/`:

| file | nodes | cells | topology → formulation |
|---|---:|---:|---|
| `dji9443_20260723_30_81_uncapped.msh`  | 4800 | 9280  | open → Neumann |
| `dji9443_20260723_30_81_capped.msh`    | 5268 | 10528 | watertight → Dirichlet |
| `dji9443_20260723_30_97_uncapped.msh`  | 5760 | 11136 | open → Neumann |
| `dji9443_20260723_30_97_capped.msh`    | 6324 | 12640 | watertight → Dirichlet |
| `dji9443_20260723_30_121_uncapped.msh` | 7200 | 13920 | open → Neumann |
| `dji9443_20260723_30_121_capped.msh`   | 7904 | 15800 | watertight → Dirichlet |

## Driver

`examples/dji9443_mesh_convergence.jl` (standalone; assembled from the Phase 1
`dji9443_circulation_audit.jl` pieces). Constants: `STUDY`, `PHASE`, `RPM = 5400.0`;
output under `data/$STUDY/$PHASE/` (+ `raw/`). Cases `dji{81,97,121}{c,u}`; `c` =
watertight/Dirichlet (`Union{ConstantSource,VortexRing}`, DBC=true), `u` = open/Neumann
(`VortexRing`, DBC=false). Kernel offsets, wake, backend, and solve setup identical to
Phase 1 so results are comparable.

- **Seeds**: from `find_dji9443_trailing_edge_indices` (not hardcoded), with
  `make_shedding_bbox` per blade; two-stage build (`noshedding` first, trace on the
  constructed body's nodes/cells, rebuild with `ensure_winding=true`).
- **Bins**: the fixed `|r/R| = 0.125:0.025:0.975` grid, trimmed to the station support
  shared by all 6 cases so `interpolate_linear` never extrapolates and every rung uses
  matched bins.
- **Modes** (env `PHASE2C_MODE`): `setup` (validate all 6, record counts/station
  ranges/bins → `topology_validation.csv`), `solve` (+`PHASE2C_CASE=<tag>`), `all` (six
  solves, coarse→fine), `analyze`.

```bash
PHASE2C_MODE=setup   julia --project examples/dji9443_mesh_convergence.jl
PHASE2C_MODE=all     julia --project examples/dji9443_mesh_convergence.jl
PHASE2C_MODE=analyze julia --project examples/dji9443_mesh_convergence.jl
```

## Decision rule (mirrors Phase 2b criteria)

- **Converged / attribution confirmed**: gap shrinks monotonically with n_airfoil and is
  ≤1% at 121 (note the stricter 0.25% too), Dirichlet per-rung ∫Γ movement decaying —
  expected oracle pattern is Neumann ~flat, Dirichlet climbing to it.
- **Finer mesh desirable**: at n_airfoil=121 the gap is still >1% or Dirichlet per-rung
  change hasn't decayed below ~0.25–0.5%. Report the trend and recommend resolutions to
  generate (e.g. 30_145, 30_161).
- **Attribution challenged**: gap flat or non-monotonic under chordwise refinement —
  report to Ryan before touching Phase 3/5.

Phase 1 reference gaps (uncapped/Neumann above capped/Dirichlet): +6.519% at
n_airfoil=41, +4.853% at n_airfoil=57 (different 40/56 spanwise family — orientation
only). Sanity: the gap at n_airfoil=81 should sit below Phase 1's 4.85% if the
attribution holds.

## Deliverables

- `data/$STUDY/$PHASE/topology_validation.csv` — per-mesh counts, station ranges, bins.
- `data/$STUDY/$PHASE/raw/<tag>_circulation.csv` — per blade/section circulation.
- `data/$STUDY/$PHASE/fixed_bin_circulation.csv`, `case_metrics.csv`.
- `data/$STUDY/$PHASE/phase_02c_report.md` — convergence table (∫Γ, Ωr-weighted,
  outboard) with per-rung Dir–Neu gap, Phase 1 reference rows, and the decision.
- Decision written into the Phase 2c log snapshot and the dashboard Progress line; stop
  for Ryan's approval before Phase 3.

## Checkboxes

- [ ] Setup and topology validation (all 6 meshes)
- [ ] Matched steady solves (3 rungs × 2 formulations)
- [ ] Convergence table and decision report
