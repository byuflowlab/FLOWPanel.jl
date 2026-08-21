# Rigid-motion tree/cache reuse (transform_plan!) — staged item

**Status: EXECUTING 2026-08-20 (Ryan's go via launch prompt).** Extends the deferred
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

## Evidence pointers

- Dialogue + rulings: `phase_02_single_step_benchmarks.md` Log 2026-08-19
  entries; `log.md` same date.
- Cache invariance basis: `nearfield_matrix_cache_plan.md` §2 (linearity),
  §3 integration point 5 (the originally-deferred persistence contract).
- Landed prerequisites: FastMultipole 72d3f3d/02f071c/0ef4e83/204188a/1ec0af9.
