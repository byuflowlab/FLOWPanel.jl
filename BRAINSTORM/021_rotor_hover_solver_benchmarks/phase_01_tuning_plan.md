# Phase 1 — BC-error metric + per-rung/per-solver tuning plan (STAGED)

**Status:** STAGED 2026-08-15 (Ryan-approved to stage, NOT yet executed).
Written for a fresh agent; high-level by design. Companion context:
`phase_01_consistency.md` Log (all 2026-08-14/15 entries), `decision_rules.md`
PROPOSED block (to be superseded by Stage 4), CSVs of record under
`benchmark/results/phase1/multi/` (fgs_check.csv, agreement.csv, solvetime.csv
= production-knob BASELINE numbers).

## Ryan rulings (2026-08-15, this plan's basis)

1. **Primary metric = boundary-condition error**: induced potential at
   control points vs the BC (for the Dirichlet campaign case, the interior
   perturbation-potential residual), relative L2, **target 1e-6**, evaluated
   with FMM using the error-tolerance feature (no O(N²) direct needed) —
   works identically at R4–R7. Replaces reference-based strength thresholds
   (the relL2 ≤ 3e-5 proposal is moot).
2. **NO production-knob freeze.** Tune FMM parameters per rung per solver
   for optimal performance.
3. **Shared knobs**: one tuned apply-knob set per rung shared by all Krylov
   configs (gmres / jacobi / ilu / fgmres); one FGS knob set shared between
   FGS-as-solver and FGS-as-preconditioner, selected by the τ-ladder
   procedure (Stage 3).
4. **FGS 1e-6 landing**: wait for tuning; if FGS still overshoots 1e-6
   badly, report BOTH rows — cost at the 1e-6 target AND at the last iterate
   above 1e-6, its error snapped to the nearest half decade.
5. (Carried judgment, Ryan may override) Tuning objective = **per-solve wall
   time subject to certified BC ≤ target, setup timed separately** (production
   reuses solver objects across steps; Phase 3 owns amortization).
6. (Consequence) The earlier MAC ≥ 0.6 bound-certification rejection is
   WITHDRAWN: the per-run certified BC-error verification is the guard at
   every rung, so any tuner point that verifies is admissible.

## Stages

**Stage 1 — Certified BC-error evaluator.** Single FMM pass with BOTH
strength columns loaded (fixed σ + solved μ): the control-point potential IS
the Dirichlet BC residual (zero interior perturbation). Dynamic
`error_tolerance` (PowerAbsolutePotential at 0.1 × target × rms_b) so
evaluation error inflates the metric ≤10% at threshold; re-certifies against
actual solved strengths each call. Normalize by RMS(φ_σ) once per rung.
Validate once against the direct evaluator at R1–R3 (dense machinery already
in place), then it becomes a standard post-solve column in every benchmark
CSV. Frozen shared BC definition (perturbation gauge, φ∞ EXCLUDED, σ from
apparent velocity, panel kerneloffset) — this preserves detection of
convention drift (cf. the coupled-φ∞ gauge finding, phase_01 log).

**Stage 2 — Unified FGS tuner.** Coordinate-descent (or small staged grid)
over {expansion_order, multipole_threshold, leaf_size, inner_iterations},
objective per ruling 5, target = a given BC tolerance τ. tune_fmm cannot do
this (no model of leaf-as-block-size or inner_iterations; measured picking
anti-optimal leaf). Run the ladder τ ∈ {1e-1 … 1e-6} (half-decade or decade
steps) per rung; **the solver role is simply the τ = 1e-6 rung**. Use
per-outer-iteration history (existing `_solve_history!` + BC evaluator) so
one instrumented solve yields the whole (cost, error) staircase.

**Stage 3 — Preconditioner selection via the τ-ladder.** For each τ rung:
translate the tuned config into a FIXED sweep/inner count achieving ≈τ
(fixed counts keep each apply a linear map from a zero seed — the clean
FGMRES regime; do NOT use tolerance-stopped applies), wrap as
FGSPreconditioner, measure end-to-end preconditioned-GMRES wall time to
BC 1e-6 (outer apply = the rung's tuned Krylov apply knobs). Best τ wins;
that single FGS knob set is then shared solver/preconditioner per ruling 3.

**Stage 4 — decision_rules.md rewrite.** Strike the 2026-08-14 knob-freeze
PROPOSED block; replace with: BC-metric definition + target; the two
per-rung tuning procedures (tune_fmm for the shared Krylov apply; the FGS
τ-ladder); shared-knob policy; dual-row FGS reporting rule; certified
verification as the per-run guard.

**Stage 5 — Regenerate the solver × rung table (R1–R3).** All configs at
tuned knobs: setup, per-solve time, iterations, achieved certified BC error;
FGS dual-row where applicable; keep the production-knob solvetime.csv
numbers as the "untuned baseline" column/table to show what tuning buys.
Also fill the small pending item: `backslash_ldiv` rows R1–R3
(config already added to `benchmark/rotor_hover_solver_phase1_solvetime.jl`,
never yet run — the last solvetime pass launched before the edit).

**Stage 6 — HPC (R4–R7).** Same procedures per rung (tune_fmm needs no
dense reference; FGS tuner needs none either; BC evaluator is
reference-free). PressureLaplace CG itmax must be raised (hit 1000-cap at
R3). ILU direct-list pattern density grows with N on this mesh (needed
max_pattern_entries=2048N at R3) — treat ILU pattern knobs as part of its
per-rung tuning.

## Validation gates

- Stage 1: FMM-certified BC error vs direct BC error agree (≤10% relative
  discrepancy at the 1e-6 scale) on R1–R3 across ≥3 solver configs.
- Stage 2/3: every selected point re-verifies BC ≤ target via the certified
  evaluator on a fresh cold solve.
- Stage 5: table numbers reproduce within run-to-run noise (k reps per
  ruling 5 where affordable; local numbers remain indicative — published
  numbers are HPC, ruling 5).

## Operational constraints (carry-forward)

Local ≤4 threads total; long background jobs on Ryan's Mac get killed
(detach + chunk per rung/config; 16 GB RAM — never two dense Gs live);
judge from CSVs not stdout; FastMultipole fixes remain uncommitted in
`../FastMultipole` (don't commit; don't touch
scripts/20250404_prediction_accuracy.jl); never write the lab notebook
without Ryan's approval.
