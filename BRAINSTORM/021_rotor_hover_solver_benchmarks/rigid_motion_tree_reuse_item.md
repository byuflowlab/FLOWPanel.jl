# Rigid-motion tree/cache reuse (transform_plan!) — staged item

**Status: COMPLETE 2026-08-20 (executed on Ryan's go; A/B collapse measured
— staleness eliminated, flat 2.1e-3 to 80° vs 21% stale; residual wake-on
FGS fixed-point discrepancy ~2e-3 flagged for Phase 3 (tolerance
hypothesis refuted). See "Execution results" below + phase_02 Log.)** Extends the deferred
near-field-cache commit 4 (`persistent_plan`); scoped during the
2026-08-19 feedback dialogue (Ryan: "seems like it would be simple work to
rotate/translate all the tree cell centers and still reuse it"). Design the
machinery at the TREE level so it serves `FmmPlan` AND `FastGaussSeidel`
(Ryan directive — the same parts should be reusable).

## Idea

Under RIGID body motion (rotation R + translation t), everything the FMM
data structures encode about *relative* geometry is invariant:

- leaf/body assignments: unchanged (same bodies, same subdivision);
- branch radii: rotation-invariant;
- interaction lists (m2l + direct): built from center distances, radii, and
  MAC — **exactly** invariant;
- expansions: recomputed from buffers every solve anyway, so once branch
  CENTERS are transformed (`center → R·center + t`), accuracy is restored.

So a `transform_tree!(tree, R, t)` (rotate/translate branch centers; buffers
refresh from the systems as they already do) + thin wrappers
`transform_plan!(plan, R, t)` / FGS equivalent makes the whole
tree/list/cache bundle reusable across timesteps of a rigidly moving body.
Caveat to verify: any consumer of the stored axis-aligned box half-dimensions
(error bounds, shrink bookkeeping) — radius-based MAC and fixed-p production
applies do not touch them, but audit before trusting adaptive-P paths.

## Payoffs

1. **Cross-timestep NearfieldInfluenceCache persistence for the rotor** (the
   original plan's deferred headline; "frame tricks" no longer out of scope):
   for the Dirichlet scalar-potential operator the cached blocks are
   **exactly invariant** under rigid motion (scalar kernel of relative
   distances); velocity outputs transform per-block as `G → R·G` (cheap
   post-rotation or body-frame evaluation). Build cost then amortizes over
   the RUN, not one solve — this changes the cached-knob economics (the
   per-solve-rebuild objection to big-leaf knobs evaporates) and the
   `*_nfcache` benchmark configs should be re-examined once landed.
2. **Fixes the FGS unsteady staleness bug (below).**
3. Removes per-solve tree/list builds for unsteady Krylov (`cache_tree`
   would persist via the plan across steps).

## FGS UNSTEADY STALENESS FLAG (correctness — verify FIRST when executing)

`FastGaussSeidel` builds trees + interaction lists + dense self/nonself
influence matrices ONCE at construction (`FastMultipole/src/solve.jl:618-665`);
`solve!` refreshes buffer positions/strengths per call (`:1100`) but never
rebuilds trees/lists/matrices. The unsteady driver
(`benchmark/rotor_hover_solver_unsteady.jl:156`) constructs the config solver
once and `simulate!` reuses it every step while `maneuver!` rotates the
rotor. Consequences:

- near-field dense matrices: **exactly valid** under rigid rotation
  (rotation-invariant scalar kernel) — not the problem;
- far field: expansions are formed from CURRENT positions about **t=0 branch
  centers**; as rotation accumulates, bodies leave their branches and the
  MAC that justified the lists no longer holds — error grows with rotation
  angle, silently. KrylovSolver is immune (tree rebuilt per apply/solve).
- Phase-1 campaign (steady single solves): unaffected.
- **Until verified, no Phase-3 unsteady FGS row is trusted.**

Discriminator (cheap, first task on execution): short R1 unsteady run, few
steps, FGS vs krylov_gmres per-step strength divergence vs rotation angle
(expect growth with angle if the bug is real; bitwise-class agreement at
step 1).

Fix via this item: rotate FGS tree centers per step (matrices stay valid);
interim mitigation if needed before this lands: rebuild the FGSSolver per
step in unsteady runs (costly: full matrix rebuild each step).

## Sketch (execution order)

1. FGS staleness discriminator (above) — establishes stakes.
2. `transform_tree!(tree, R, t)` in FastMultipole + invariance audit of box
   consumers; unit test: fmm! on a rotated system with a transformed tree ==
   fresh-tree fmm! (rtol 1e-12), incl. panels with vertices.
3. `transform_plan!` + cache validity: Dirichlet exactness test under
   rotation (cache reused, tree transformed) == fresh build; velocity-output
   handling (R·G) or explicit refusal for gradient switches in v1.
4. FGS wrapper (transform its trees; nonself/self matrices untouched) +
   unsteady A/B vs Krylov.
5. FLOWPanel surface: per-step `transform!` hook from the kinematics
   (rigid frames already carry R, t per step) + `persistent_plan` opt-in
   (the original commit-4 contract, now with motion support).
6. Re-examine `*_nfcache` knob economics with run-amortized builds.

## Execution results (2026-08-20)

1. **Staleness discriminator (sketch step 1): BUG CONFIRMED.** New driver
   `benchmark/rigid_motion_fgs_staleness.jl` (R1, 8 steps, NT=36 ⇒ 10°/step,
   pure-FGS trajectory; after every production FGS solve the SAME state is
   re-solved by a fresh KrylovSolver gmres rtol 1e-9 — trees rebuilt per
   apply — isolating solver error from wake feedback). Per-step relative L2
   divergence of μ: 1.74e-7 at 0° (solver-tolerance class, as pre-registered),
   then monotonic growth: 2.0e-3 @10°, 8.9e-3 @40°, 4.3e-2 @60°, **2.1e-1 @
   80°** (~10⁶× growth). Data:
   `rigid_motion_tree_reuse/fgs_staleness.csv` (arm `fgs_vs_fresh_krylov`).
   **No pre-fix Phase-3 unsteady FGS row is trustworthy.**
2. **FastMultipole machinery landed** (branch flowpanel-20260817):
   - `transform_tree!` (eea944d): centers → R·c+t, box → |R|·box (tight AABB
     of the rotated box; enclosure preserved, error bounds conservative);
     proper-rotation guard. NOTE Ryan's concurrent commit 645cc96 swept the
     src half of this in with his autotune work; eea944d carries the tests.
     Box-consumer audit: post-build, only the adaptive-P/error-bound paths
     read `box` (via `minimum_distance`) — |R|·box keeps them conservative;
     the active MAC is radius-based; m2l operators are formed per call from
     current centers, so index-pair lists stay valid.
   - `transform_plan!` (087bf4a): + target-buffer position refresh (frozen at
     plan build — planned fmm! only zeroes output rows); nearfield cache
     persists EXACTLY for scalar outputs across accumulated motions (tested
     rtol 1e-12, two-step composition); gradient/hessian/extra + cache ⇒
     loud v1 refusal (per-block G→R·G not implemented).
   - `transform_solver!` for FGS (d714544): same tree machinery + target
     refresh + `transformed` flag; dense matrices untouched (scalar rows
     exactly invariant); gradient solves on a transformed solver refuse.
     FM-level A/B (cold fixed-iteration trajectory): transformed replay
     <1e-9 after 63°+translation, stale >1e-3 (>1000×).
   - Test-system findings: a fresh tree on rotated positions is NOT a valid
     reference (axis-aligned subdivision degenerates, measured 160→0 m2l
     pairs — equivariance vs the original run is the right 1e-12 test); the
     vortex-filament test system's velocity carries an origin-dependent
     Lamb-Helmholtz gauge term (t_z/2-class shifts under pure translation) —
     replaced by a self-contained triangular source-panel system for
     vertex-path coverage; latent `size(x,2)` bug fixed in both filament
     generators.
3. **FLOWPanel surface** (uncommitted, in the working tree):
   `propagate_kinematics!` returns per-body rigid (R, t) deltas; `simulate!`
   forwards them via `transform_body_solvers!` → `transform_solver_geometry!`
   (KrylovSolver persistent_plan → transform_plan!; FGSSolver →
   transform_solver!; others no-op — Backslash's dense operators are
   rotation-invariant); `KrylovSolver(persistent_plan=true)` keeps the plan
   (+ nearfield cache) across solves (deferred commit 4, now with motion
   support; kerneloffset restored before every solve keeps it sound);
   metadata `persistent_plan`; unit tests incl. Dirichlet
   shedding-panel cache exactness under rigid motion (rtol 1e-12).
4. **Knob-economics flag (sketch step 6, docs only):** with run-amortized
   builds (persistent_plan), the per-solve-rebuild objection to big-leaf
   `*_nfcache` knobs evaporates — build cost amortizes over the RUN. The
   Phase-2b `*_nfcache` configs should be RE-EXAMINED on HPC once this lands
   in the benchmark drivers. NOT rerun now (chains 13242659–66 in flight;
   rsync freeze).

5. **A/B collapse (production path, measured 2026-08-20 late):** fixed arm
   flat 2.0–2.2e-3 from 10° through 80° (stale: 21% at 80°; ~100×) —
   angle-dependence eliminated. The ~2e-3 plateau is a pre-existing
   FGS-vs-Krylov discrepancy that appears with the first wake-carrying step
   in BOTH arms (step 1: 2.014e-3 stale vs 2.007e-3 fixed).
   Stopping-tolerance hypothesis TESTED AND REFUTED (tolerance ×0.01 arm:
   steps 1–2 numerically unchanged, step 0 tightened to 2.9e-8): the FGS
   converges to a fixed point ~2e-3 from the Krylov solution once free wake
   particles exist — tolerance- and angle-independent, extra_farfield
   consistent; open Phase-3 item to locate the wake-on system asymmetry.
   Also landed: FM test 149bd9f4 "cells exactly follow the bodies" (direct
   drift/containment invariant; documents that under recenter=false `box`
   was never a strict enclosure about `center` even at build — 276/1500
   bodies outside, worst 0.027 — the transform preserves that character).
   Data: all arms in `rigid_motion_tree_reuse/fgs_staleness.csv`.

## Evidence pointers

- Dialogue + rulings: `phase_02_single_step_benchmarks.md` Log 2026-08-19
  entries; `log.md` same date.
- Cache invariance basis: `nearfield_matrix_cache_plan.md` §2 (linearity),
  §3 integration point 5 (the originally-deferred persistence contract).
- Landed prerequisites: FastMultipole 72d3f3d/02f071c/0ef4e83/204188a/1ec0af9.
