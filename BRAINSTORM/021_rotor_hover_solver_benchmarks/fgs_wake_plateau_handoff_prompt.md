# Handoff prompt — diagnose + fix the wake-on FGS↔Krylov ~2e-3 fixed-point plateau

**Status: STAGED 2026-08-21 — launching this prompt is Ryan's go.**
Successor to `rigid_motion_tree_reuse_item.md` (COMPLETE — read its
"Execution results" first). Copy everything between the rules into a fresh
agent.

---

Diagnose and FIX the wake-on FGS-vs-Krylov solution discrepancy found while
closing BRAINSTORM 021's rigid-motion tree-reuse item. Required policy reads
before touching src: `agent_policies/WORKFLOW.md`, `agent_policies/TESTING.md`.
Background docs (read in this order, no broad re-exploration needed):
`BRAINSTORM/021_rotor_hover_solver_benchmarks/rigid_motion_tree_reuse_item.md`
("Execution results"), `phase_02_single_step_benchmarks.md` Log 2026-08-20
tree-reuse entry (its "Plateau root cause" block is the problem statement).

## The problem (all measured, R1 = 8016-panel DJI rotor, RHPC unsteady)

Per-step relative L2 divergence of the solved Dirichlet μ column between the
production `FGSSolver` solve and a fresh `KrylovSolver` (:gmres, rtol 1e-9,
trees rebuilt per apply) re-solving the IDENTICAL state:

- step 0 (no free wake particles): 1.7e-7 agreement — 2.9e-8 when the FGS
  tolerance is tightened 100× (so step 0 is tolerance-limited and the
  machinery agrees); μ norms 6.5–12 everywhere (non-vacuous).
- steps ≥1 (free particles exist): ~2.0e-3 → 2.2e-3, FLAT in rotation angle.
- Tolerance ×0.01 leaves steps 1–2 numerically IDENTICAL to 4 digits
  (2.007e-3 / 2.073e-3) ⇒ the FGS CONVERGES — to a fixed point ~2e-3 from
  the Krylov solution. It is a different-linear-system problem, not an
  iteration/tolerance problem.
- Identical in stale-tree and transformed-tree arms at step 1 (2.014e-3 vs
  2.007e-3) ⇒ NOT tree staleness (that bug is fixed and separately verified).

Data: `BRAINSTORM/021_rotor_hover_solver_benchmarks/rigid_motion_tree_reuse/
fgs_staleness.csv`, arms `fgs_vs_fresh_krylov`,
`fgs_transformed_vs_fresh_krylov`, `fgs_transformed_tol0.01_vs_fresh_krylov`.

Already CLEARED — do not re-litigate: FGS stopping tolerance (refuted by the
×0.01 arm); tree staleness (fixed; FM transform_solver! + simulate! hook);
`extra_farfield` plumbing (FLOWPanel passes
`extra_farfield=any(has_semiinfinite_wake.(bodies))` at
`src/FLOWPanel_solver.jl:906`); FGS warm-start projection (initial guess
only); vacuous-norm concerns.

Key discriminating fact to explain: the attached wake (Das strips, Kutta
coupling) exists at step 0 too and the paths AGREE there; the discrepancy
switches on exactly when free wake particles first exist (their influence
enters through `body.velocity` → σ) and the staged-startup freestream ramps.

## Hypotheses (ranked; H-numbers for the log)

- H1 σ/RHS asymmetry: FLOWPanel Dirichlet Krylov RHS = `-body.potential`
  after the outer `solve!` zeroes potential and assembles the interior
  source potential (`src/FLOWPanel_solver.jl:238-257`, `_solve!` at `:577`).
  The FGS path instead builds its external RHS inside FastMultipole
  `solve!` (`FastMultipole/src/solve.jl` ~1080–1110:
  `target_influence_to_buffer!` → `influence!` → `right_hand_side .=`)
  plus FLOWPanel `_solve!(body, ::FGSSolver)` (`FLOWPanel_solver.jl` ~953,
  prior_sigma save/restore, `phi_ext`/`Uext` storage). Any difference in
  what each treats as "external" (source potential, wake correction fields,
  Kutta terms) lands directly in the fixed point.
- H2 operator asymmetry on the SAME μ: Krylov apply = `_apply_dirichlet_G!`
  (includes `_induced_wake` attached-wake term; V0-verified linear). FGS
  operator = dense self/nonself blocks (probed through `direct!` at
  construction) + its own far field. If the attached-wake term is present
  in one and not (or differently) in the other for some pair class, the
  systems differ. Must explain the step-0 agreement gate (e.g. a term
  proportional to something that is zero until particles/shed rows exist —
  wake correction columns, shed-row coupling, Das kinematic refresh).
- H3 kerneloffset state: `_set_kerneloffsets!` sequencing may leave a
  different active offset for the FGS solve than for the probe's reference
  solve (the probe calls `pnl.solve!(body, ref; backend=backend_solve)`
  immediately after the production solve inside a formulation wrapper —
  see the driver). Offset differences change the regularized near field.
- H4 wake-correction fields on the body (columns some path includes as
  external and the other bakes into the solve) — check
  `clear_wake_correction!` users and what `VelocityThroughSources` leaves
  on the body at solve time.

## Recommended attack (in order)

1. **Arbitrate with ground truth**: extend the discriminator driver
   (`benchmark/rigid_motion_fgs_staleness.jl` — knobs: `RUNG`, `N_STEPS`,
   `DIAG=1` relaxes premise guards for short runs, `FGS_TOL_SCALE`, `ARM`
   labels CSV rows) with a THIRD same-state re-solve per step:
   `pnl.Backslash` dense LU (exact; ~0.5 GB at R1, fine locally). One
   `DIAG=1 N_STEPS=2` run answers WHO is wrong: FGS 2e-3 off, Krylov 2e-3
   off, or both. Restore strengths after each probe solve exactly as the
   existing wrapper does (full `body.strength` snapshot — `solve!` mutates
   nothing else; verified).
2. **Bisect system components on the frozen step-1 state**: inside the same
   wrapper, with identical body state, (a) compare RHS vectors — FLOWPanel's
   assembled `-body.potential` vs the FGS `self_matrices.rhs` after its
   external-influence assembly (mind orderings: FGS is leaf/tree-ordered;
   use `sorted_index_2_unsorted_index` or compare through applied results);
   (b) compare operator action — `pnl._apply_dirichlet_G!(body, x, ...)` vs
   the FGS matrices+far-field applied to the same x (a single fixed-sweep
   `solve!` with tolerance=-1 and known input can serve as an apply probe).
   Whichever half disagrees localizes H1 vs H2.
3. **Minimal reproducer**: once the term is identified, reproduce on the
   `make_dirichlet_diamond_body` fixture (`test/test_helpers.jl`; premise
   guard: `nsheddings > 0`) + a few hand-placed particles feeding
   `body.velocity` — no simulate! needed. That becomes the regression test.
4. **Fix test-gated**: FastMultipole changes = commits on
   `flowpanel-20260817` (stage ONLY files you touched; NEVER
   `scripts/20250404_prediction_accuracy.jl`; Ryan lands unrelated commits
   on this branch daily — commit on top, never rebase/amend others').
   FLOWPanel src stays UNCOMMITTED like the rest of its dirty tree.
5. **Verify**: FM full suite (`julia --project=test -t 4 test/runtests.jl`)
   + FLOWPanel `runtests_unit_solver.jl` / `runtests_unit_fmm.jl` /
   `runtests_unit_simulate.jl` at 4T; then the discriminator
   `RUNG=R1 N_STEPS=8` fixed arm should drop from ~2e-3 to the ~1e-7
   tolerance class at ALL angles — that is the acceptance criterion.

## Constraints + gotchas

- Local runs ≤4 threads TOTAL; a full 8-step discriminator run ≈ 25–30 min;
  `DIAG=1 N_STEPS=2` ≈ 8 min. Nothing needs HPC. rsync NOTHING to the
  cluster (Phase-2b chains 13242659–66 + Phase-1 remnants draining); do not
  touch `benchmark/rotor_hover_solver_phase1*` or `benchmark/results/`.
- Julia exit 0 ≠ health (LoadError can exit 0); never pipe verdict-bearing
  output through `tail`/`grep` without keeping the failure lines; silent
  log end = error, rerun foreground.
- Every new test needs a non-vacuous premise guard (campaign convention).
- FLOWVPM.jl is Ryan's active merge-work area — if precompile fails with
  conflict markers there, STOP and report, don't touch it.
- Judgment calls: make them, log them, flag them — don't block. Log as you
  go: `phase_02_single_step_benchmarks.md` Log (dated entry) + `log.md`
  (newest-first) + memory `project_021_solver_benchmarks.md` at the end.
  No lab-notebook writes — offer an entry instead.

Deliverable: root cause named with the H-number and the measurement that
pinned it; fix landed (FM committed / FLOWPanel uncommitted) with a
regression test; discriminator fixed-arm table flat at ~1e-7 class; suites
green; docs + memory updated.

---
