# Phase 1 Stage 6 — R4–R7 HPC procedure (SPEC)

**Status:** SPEC drafted 2026-08-15 (autonomous session, per
`phase_01_tuning_plan.md` Stage 6). Execution is a fresh-agent/HPC task; read
`agent_policies/HPC.md` FIRST (allocation, threading, submission boundary;
`ssh orc` needs a live ControlMaster socket, 2FA otherwise). Cluster repo:
`/home/rander39/projects/FLOWPanel.jl` (Manifest ready).

## What runs per rung (R4–R7; same scripts as R1–R3)

Ladder (ledger.md, capped freeze): R4 `65_209` 58,192 / R5 `89_289` 108,240 /
R6 `125_409` 212,108 / R7 `177_577` 419,276 panels
(`dji9443_20260813_<r>_capped_captess4.msh`). Extend `LADDER` and `TUNED` in
`benchmark/phase1_case.jl` per rung as knobs land.

Per rung, in order (each its own Slurm job or job step; judge from CSVs):

1. **Case + frozen b** (`phase1_case.jl`): the O(N²) direct source-potential
   assembly is the dominant fixed cost — projected from R3 (~60 s @ 4T):
   R7 ≈ 210× R3 work; budget hours at low thread counts, minutes at 64T.
   Persist `b`/`rms_b` per rung if jobs are split (cf. `agreement_xref_R3.bin`
   pattern).
2. **tune_fmm** → the rung's shared Krylov apply knobs (ruling 3; no dense
   reference needed). Same call as `rotor_hover_solver_phase1.jl` TUNE=1:
   `PowerAbsolutePotential(0.1·1e-6·rms_b)`, `max_expansion_order=20`,
   MACs 0.3–0.7. Requires a solved strength state — use an FGS solve at
   production knobs first (tune with THAT x loaded), or backslash where dense
   G fits (R4–R5 only). Record to `tune.csv`. Caveat (Phase 1 log): tuner
   degenerates if requested tol < true floor — sanity-check the tuned point's
   measured rel_rms.
3. **BC-evaluator spot-validation at R4 ONLY** (`..._bcerror.jl`): the direct
   path is O(N²) per point — feasible at R4 (~9× R3), pointless-expensive at
   R6–R7. R1–R3 + R4 validation (≤10% agreement) is the certification basis;
   R5–R7 rely on the certified bound (`error_success` must be true on every
   row; any false ⇒ raise the P cap above 20 and re-run that evaluation).
4. **FGS τ-ladder tuner** (`..._fgstune.jl`): reference-free. Watch
   `MAXIT=40` — if a rung needs more outer iterations to reach 1e-6, raise
   it. Descent grids may need extending upward in leaf (GS block size) as N
   grows; add larger LEAF_SET values if the R1–R3 optimum sits at the grid
   edge.
5. **Preconditioner τ-selection** (`..._fgsprecond.jl`) → the rung's shared
   FGS knob set.
6. **Solver × rung table** (`..._table.jl`): all configs at tuned knobs,
   k_reps ≥ 5 (ruling 5 — published numbers are HPC min-of-k), both
   threading modes (ruling 6: 1T+BLAS1, and 64T with recorded BLAS count;
   never mix modes in one comparison). `backslash_ldiv`/dense configs run
   R4–R5 and DROP at R6–R7 (dense G: R5 ≈ 94 GiB fits a large node; R6 ≈
   360 GiB does not — ledger cost-ceiling rule; drop rows per config).

## Per-rung tuning knobs that grow with N (treat as part of tuning)

- **ILU pattern**: the Barba direct-list (leaf 10 / MAC 1.0) density grows
  with N on this mesh — R3 already needed `max_pattern_entries = 2048N`.
  Treat `max_pattern_entries` (and, if the guard still trips, pattern
  leaf/MAC) as per-rung ILU tuning; record the setting in the CSV notes.
- **PressureLaplace CG itmax**: hit its 1000 cap at R3 (any agreement/monitor
  stage on HPC must raise it and record the setting).
- **FGS `max_iterations` / tuner MAXIT** and **LEAF_SET upper end**: see 4.

## Operational cautions (bitten before — 018/Phase-1 lessons)

- **FastMultipole fixes are UNCOMMITTED** in the local `../FastMultipole`
  dev checkout (fm `5adde3b-dirty`: two tune_fmm fixes + the FGS shared-tree
  replay ctor). The cluster CANNOT get them via git pull: rsync the working
  tree (excluding `scripts/20250404_prediction_accuracy.jl`, which is dirty
  from unrelated experimentation and must not travel or be committed), then
  **md5-compare `src/` cluster-vs-local** before any run (the degenerate
  r_over_R incident came from exactly this class of mismatch). Same check
  for FLOWPanel `src/` (also dirty-uncommitted).
- Meshes: confirm the four capped `.msh` files exist on the cluster under
  `examples/data/` (they are large; rsync, don't regenerate — OpenVSP
  chordwise rounding made regeneration non-reproducible, Phase 1 log
  2026-08-13).
- Verify knobs from the banner (`banner.txt` lands next to each CSV); judge
  every run from its CSVs; exit 0 ≠ health.
- Radius fix must be ACTIVE (`FMM_RADIUS_TOL[] = 1e-6` default) — recorded
  per row in the `radius_tol` column; a row showing `Inf` is misconfigured.
- Sliver-kernel fix (relative tan_term tolerance) is in the FLOWPanel
  working tree — the md5 check above covers it; without it every floor and
  BC number is invalid at these meshes.

## Deliverables back from HPC

Per rung: `tune.csv`, `fgstune_{staircase,selected}.csv`, `fgstune_verify.csv`,
`fgsprecond.csv`, `solvetable.csv` (+ `banner.txt`) under
`benchmark/results/phase1/<mode>/` — merge into the local results tree under
a `hpc/` subdirectory or hardware-tagged rows (schema already carries
`hardware_tag` via banner in runs.csv-style files; solvetable carries
threading_mode + commit columns). Then Phase 1 closes with the solver × rung
tables R1–R7 and the decision_rules-frozen procedures.
