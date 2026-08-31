# Revised plan — BRAINSTORM 021 Phase 3 solution-history / block-GS fix

## Outcome

Make warm-start history advance exactly once per completed top-level solve step,
keep that history frozen throughout block Gauss-Seidel, remove the redundant
second solve for a one-body tuple, and report the first inner-solver work plus
the number of per-body solves performed in the step.

This keeps the original plan's Option 2, but tightens the lifecycle around
three distinct events:

1. a raw linear-solver kernel call (`_solve!`);
2. one per-body solve within block GS;
3. a completed top-level step, where history and published statistics commit.

Only event 3 advances warm-start history. Events 1 and 2 may update transient
step statistics only when they are inside an explicitly begun step.

## Established defect

`save_solution!` (FGS) and the Krylov `x_prev` / `x_history` update currently
run inside every `_solve!`. `solve!(bodies::Tuple, solvers::Tuple)` invokes a
body solve once per block-GS outer iteration, so history contains intra-step
repeats rather than consecutive timestep solutions. Extrapolation can then
alternate around the answer and prevent the outer `max_delta` test from
converging. The driver's current `niter` is also the last outer solve's work,
not the first solve's step-to-step warm-start benefit.

The previously recorded cold/prev/extrap smoke timings and the apparent prev
6-to-0 iteration reduction are therefore void and must be re-measured.

## Design decisions

### Keep Option 2

Move history persistence out of `_solve!` and commit it at a real top-level
step boundary. Keep projection inside `_solve!`: Dirichlet setup clears the
solved strength column before `_solve!`, so the projection must occur after
that setup on every outer iteration. Freezing the saved history makes every
outer iteration start from the same historical guess.

Do not use a body-level `solved` flag or a timestep-indexed history slot. Both
add mutable coordination state whose correctness depends on all callers
maintaining it.

### Preserve completed-iterate behavior

The current code saves the final iterate even when an iterative solver reaches
its limit. Preserve that behavior to avoid an unrelated semantic change: the
new history means one **completed top-level step state**, not necessarily a
mathematically converged state. Publish the convergence flags and let the
driver reject invalid comparisons. A later policy change may choose not to
warm-start from failed steps, but that is out of scope.

### Raw `_solve!` is not a step boundary

Do not put `_note_solve!` in FGS `_solve!` or `_krylov_launch!`. The repository
has many deliberate raw callers in Phase 1/2 benchmarks and solver-history
tests, and `FGSPreconditioner` calls its private FGS solver once per
preconditioner application. Counting there would grow transient accumulators
indefinitely and could pollute a later published step.

Instead, note successful per-body solves at their orchestration sites:

- the two public single-body `solve!` methods;
- tuple `solve!`, through those single-body methods;
- the direct TraceCorrected bypass.

Low-level `_solve!` and `_solve_history!` calls retain their kernel-level
meaning and do not publish step statistics or advance warm-start history.

## Implementation

### 1. Add explicitly scoped per-step statistics

In `src/FLOWPanel_solver.jl`, add:

```julia
mutable struct SolveStepStats
    acc_nsolves::Int
    acc_niter_first::Int
    acc_solved_all::Bool
    nsolves::Int
    niter_first::Int
    solved_all::Bool
end

SolveStepStats() = SolveStepStats(0, -1, true, 0, -1, true)
```

Add a `stats::SolveStepStats` field to `FGSSolver` and `KrylovSolver`, updating
both positional constructor calls. Record the small `summarysize` increase in
`decision_rules.md`, since `solver_state_bytes` is a campaign metric.

Define these lifecycle helpers:

```julia
begin_step_solution!(::AbstractSolver) = nothing
note_step_solve!(::AbstractSolver) = nothing
publish_step_stats!(::AbstractSolver) = nothing

step_nsolves(::AbstractSolver) = -1
step_niter_first(::AbstractSolver) = -1
step_solved(::AbstractSolver) = true
```

For FGS and Krylov:

- `begin_step_solution!` resets the transient fields to `(0, -1, true)`;
- `note_step_solve!` increments `acc_nsolves`, records `.niter` only when the
  count becomes one, and ANDs `.solved` into `acc_solved_all`;
- `publish_step_stats!` copies transient values to the published fields but
  does not own history persistence;
- accessors read the published values.

Beginning a step is load-bearing. Reset-on-publish alone is unsafe after an
exception or an intentionally unfinalized internal solve. On failure, leave
the previously published statistics unchanged; the next real step begins by
discarding any partial transient state.

`step_solved` means all *inner solver calls* reported convergence. It is not a
multibody outer-convergence flag. Phase 3 is a one-body campaign, where the
outer loop exits by construction after the first solve.

### 2. Separate history commit from statistics publication

Add:

```julia
save_step_solution!(body, ::AbstractSolver) = nothing
finalize_step_solution!(body, solver) = begin
    save_step_solution!(body, solver)
    publish_step_stats!(solver)
    nothing
end
```

- FGS `save_step_solution!` delegates to existing `save_solution!`.
- Krylov `save_step_solution!` copies the final workspace solution into
  `x_prev`, sets `have_x_prev`, and shifts/writes `x_history` once.
- Backslash and other solvers use the no-op defaults.

Delete FGS `save_solution!` from `_solve!`. Delete Krylov's `x_prev` and
`x_history` persistence from `_krylov_launch!`. Do not add accounting to either
raw kernel.

The distinction between `save_step_solution!` and `publish_step_stats!` keeps
failure behavior explicit and makes each part independently testable.

### 3. Make top-level step boundaries explicit

Add an internal orchestration keyword to both single-body public methods:

```julia
solve!(body, solver; finalize_step::Bool=true, ...)
```

When `finalize_step=true`:

1. call `begin_step_solution!(solver)` before setup;
2. call `_solve!`;
3. after `_solve!` succeeds, call `note_step_solve!(solver)`;
4. restore temporary Dirichlet state;
5. call `finalize_step_solution!(body, solver)` only after the public solve has
   completed successfully.

When `finalize_step=false`, do not begin or finalize, but do note a successful
per-body solve. Document this as an internal orchestration mode: its caller
must bracket calls with begin/finalize. Ensure `solver_optargs` cannot also
supply `finalize_step`, so duplicate keyword behavior cannot silently change
the lifecycle.

For tuple `solve!(bodies, solvers)`:

1. call `begin_step_solution!` once for every solver before the outer loop;
2. invoke each single-body solve with `finalize_step=false`;
3. restore body velocities after the loop;
4. finalize each body/solver exactly once only after tuple orchestration
   completes normally.

If a body solve throws, do not save history or publish partial statistics. The
next top-level call resets transient counters at begin-step.

For `TraceCorrected` in `src/FLOWPanel_formulation.jl`:

1. begin the step before its direct `_solve!` path;
2. note the solve after `_solve!` succeeds;
3. restore `body.potential`;
4. apply `set_wake_correction!`;
5. only then finalize.

Finalizing inside the existing `try` immediately after `_solve!` is too early:
the formulation has not completed until the physical wake correction is
installed.

GreenReconstruction and DirectWakePotential use Backslash `ldiv!` directly and
need no history hook. Raw Phase 1/2 benchmark calls remain intentionally
unpublished. The FGS preconditioner's inner raw `_solve!` remains completely
outside outer-solver step accounting.

### 4. Short-circuit a one-body tuple

After recording the first block-GS `max_delta`, but before the tolerance test,
exit when `N == 1`:

```julia
if N == 1
    converged = true
    verbose && println("  Single body: one block solve completes the outer loop")
    break
end
```

With no other body, `sources` is empty and the block-GS map has no coupling
update to iterate. The `ConvergenceHistory` contract changes intentionally: a
one-body tuple records exactly one outer entry. Multibody behavior is unchanged.

Do not add a production switch solely to test the old redundant path. Verify
the optimization by comparing the one-tuple result to an equivalent standalone
single-body solve on a cloned body/solver and by a focused frozen-history test
showing two manually repeated raw solves return bit-identical strengths.

### 5. Correct driver instrumentation and make the CSV migration safe

In `benchmark/rotor_hover_solver_unsteady.jl`:

- retain `niter` as the last inner solve for continuity with Phase 1/2;
- add `niter_first` and `nsolves` vectors to `TimedFormulation`;
- read them through `step_niter_first` and `step_nsolves` after the formulation
  returns;
- keep `solved` sourced from `step_solved` and document that it is the AND of
  inner convergence flags;
- make `niter_first` the warm-start headline metric;
- for iterative one-body runs, require `nsolves == 1` on every wake-developed
  row; warn loudly and label the comparison invalid otherwise;
- treat `nsolves == -1` as unavailable for generic/direct solvers, not as a
  successful check;
- do not generalize “constant nsolves” into a multibody validity rule, because
  genuine block-GS outer counts may vary by timestep.

Before appending, compare the exact existing first line with the new expected
header. Refuse to append when it differs, with an error directing the caller to
a fresh/versioned output path. Otherwise old-width and new-width rows could be
mixed under one header and silently misparsed.

In `benchmark/phase3_analysis.jl`:

- require the new `niter_first` and `nsolves` columns for Phase 3 inputs;
- store `nsolves` as an actual `Arm` field, not only in a rendered guard string;
- use `niter_first` for reduction tables, order sweeps, and TikZ data;
- incorporate the one-body `nsolves == 1` and all-inner-solved checks into the
  guard;
- keep `niter` available only as a diagnostic legacy field.

### 6. Correct the campaign record

Update `BRAINSTORM/021_rotor_hover_solver_benchmarks/decision_rules.md`:

- replace the incorrect claim that the one-body block-GS loop was degenerate;
- state the new history/step-boundary and one-body single-solve contracts;
- add `niter_first` and `nsolves` to the Phase 3 schema;
- state that legacy `niter` is the last inner solve and not the comparison
  metric;
- record the measured `SolveStepStats` memory delta.

Update `BRAINSTORM/021_rotor_hover_solver_benchmarks/log.md` and every place the
smoke table is cited:

- mark the prev 6-to-0 reduction and cold 4.378 s / prev 3.872 s / extrap
  81.13 s comparison void;
- explain both artifacts: last-solve sampling and redundant/runaway outer
  solves;
- record that the standing one-body repeated-solve tax is resolved here;
- replace numbers only after the post-fix smoke succeeds.

Do not touch unrelated dirty-worktree changes.

## Verification

### Narrow unit coverage

Run at most four local threads and extend the suites owning these contracts:

- `test/runtests_unit_solver.jl` for public single-body FGS/Krylov lifecycle;
- `test/runtests_unit_solver_history.jl` for block-GS outer-history semantics;
- `test/runtests_unit_fgs_history.jl` for FGS rolling history;
- the relevant formulation test for TraceCorrected's direct bypass.

Required assertions:

1. Repeated standalone FGS and Krylov solves each publish `nsolves == 1`, and
   history advances exactly once per call.
2. `finalize_step=false` followed by a normal standalone solve cannot pollute
   the next published statistics; the normal solve begins by resetting partial
   transient state.
3. A deliberately thrown inner solve neither advances history nor overwrites
   the prior published statistics; the following successful step reports clean
   counts.
4. A one-body tuple with extrapolation advances FGS/Krylov history once,
   publishes `nsolves == 1`, and records exactly one block-GS history entry.
5. A two-body tuple advances each solver history once; each solver's `nsolves`
   equals the block-GS history length; `niter_first` is captured from the first
   per-body solve and `.niter` remains the last.
6. Krylov `x_prev` and `x_history` remain unchanged between outer iterations
   and commit only at tuple completion.
7. A one-tuple solution is bit-identical to the matching standalone solve on a
   cloned body/solver. Separately, two manually repeated solves with frozen
   history produce bit-identical strengths for representative FGS and Krylov
   Dirichlet paths. Do not claim a disabled short-circuit via
   `max_outer_iterations`; that knob cannot disable an unconditional break.
8. TraceCorrected advances history/stats once after the full formulation
   succeeds. The generic exception-path tests above establish that an
   interrupted lifecycle cannot contaminate the next published step; avoid
   adding a production failure-injection seam solely for this test.
9. Existing Backslash and two-body `ConvergenceHistory` behavior remains
   unchanged except for the documented one-tuple length.

Run narrow suites first, then the full suite. The five known
`runtests_unit_kernel_gradient.jl` failures from the unrelated uncommitted
`_onplane_snap` AD-dual mismatch remain out of scope; record them rather than
changing that kernel.

### Defect reproduction

Run the six-step R1/FGS extrapolation diagnostic with four threads:

```bash
RUNG=R1 CONFIG=fgs N_STEPS=6 WARMSTART=extrap WARMSTART_ORDER=1 \
SOLVER_VERBOSE=1 PHASE=phase3_diag EXPECT_JULIA_THREADS=4 THREADING_MODE=multi \
julia --project -t 4 benchmark/rotor_hover_solver_unsteady.jl
```

Expected on every step: one per-body solve, no alternating projected/solved
pair, and `nsolves == 1`.

Then re-run the same six-step R1/FGS smoke for cold, prev, and extrap. Require
matching particle counts and strength checksums to about `1e-7` relative, all
inner solves converged, and `nsolves == 1`. Replace the voided timing table only
with these post-fix results.

Delete only the known generated diagnostic directories when finished:
`benchmark/results/phase3_smoke/` and `benchmark/results/phase3_diag/`.

## Out of scope

- HPC submission or approval-state changes;
- the `_onplane_snap` AD-dual kernel-gradient defect;
- the pending family-attribution relabel for in-flight jobs;
- committing the dirty worktree;
- changing the policy for whether failed completed steps may seed a later
  warm start;
- adding general multibody outer-convergence publication beyond the existing
  `ConvergenceHistory` and verbose warning.
