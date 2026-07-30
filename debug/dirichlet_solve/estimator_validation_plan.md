# GreenReconstruction and TraceCorrected validation handoff

## Objective and approval criterion

Test the finite-wake Dirichlet formulations `GreenReconstruction` and
`TraceCorrected` against the documented semi-infinite and old finite-wake
results for the case in this directory. The main objective is to identify the
method that closes the largest fraction of the old finite-wake lift-coefficient
gap without explicitly evaluating the free wake's scalar-potential influence
at body control points.

Rank methods by a **two-tier** closure metric. Report both fractions per row and
rank primarily on the second.

```math
f_{\mathrm{semiinf}} = 1 -
\frac{|C_{L,\mathrm{method}}-C_{L,\mathrm{semiinf}}|}
     {|C_{L,\mathrm{old}}-C_{L,\mathrm{semiinf}}|},
\qquad
f_{\mathrm{oracle}} = 1 -
\frac{|C_{L,\mathrm{method}}-C_{L,\mathrm{Task3}}|}
     {|C_{L,\mathrm{old}}-C_{L,\mathrm{Task3}}|}.
```

`f_{\mathrm{semiinf}}` is closure toward the Task 1 semi-infinite baseline (keep
as context). `f_{\mathrm{oracle}}` is closure toward the **matching Task 3
direct-potential result** (single-shot vs Task 3 single-shot; iterated vs Task 3
iterated). The latter is the primary ranking quantity because it isolates
*estimator* error — the only thing Tasks 4/5 vary — from the fixed-wake
geometry/strength error that every frozen method shares. On the rolled-wake row
the Task 3 iterated oracle is already −0.50 % (3.94°) / −7.48 % (45°) short of
semi-infinite, so `f_{\mathrm{semiinf}} → 1` is unreachable by *any* frozen
method; do not read a low `f_{\mathrm{semiinf}}` there as estimator failure.

Denominator mapping (each frozen state pairs with its own old case):

- frozen-1 flat `Da=0.5c` ↔ old flat `Da=0.5c` (`C_L=0.2742139824`);
- frozen-2 flat `Da=0.05c` ↔ old flat `Da=0.05c` (`C_L=0.2715662450`);
- frozen-3 rolled `Da=0.05c` ↔ old settled rolled (`C_L=0.265312060869`).

The rolled `Da=0.05c` row is the **primary ranking row**: it has the largest,
best-conditioned old gap (−3.44 %). The flat `Da=0.5c` row's old gap is only
−0.20 % (denominator ≈ `5.5e-4` CL), so its `f_{\mathrm{closed}}` is
numerically noisy — treat flat rows as convergence/consistency checks, not
ranking rows. Because both metrics use absolute values, a method that overshoots
past the target scores like one short of it and the fraction can exceed 1; always
report the **signed** post-correction difference alongside and treat a fraction
`> 1` as overshoot, not extra closure. "Robust gap-closure" (used below and in
the march-selection criterion) means a consistent sign and magnitude of
`f_{\mathrm{oracle}}` across **both angles and both frozen wake states**, not a
single best number.

Do not impose an arbitrary lift-error acceptance threshold. A result is usable
only if its formulation-specific residuals and immutable-state checks pass.
The user explicitly selected:

- both `TraceCorrected` estimator modes;
- the full comparison matrix at both 3.94 and 45 degrees;
- frozen-state comparisons first, followed by a self-consistent march of the
  best method;
- the dense bordered-LU Green solve as the authoritative full-resolution
  accuracy reference.

## Mandatory repository workflow

Before work, read the repository-root `CLAUDE.md`, then
`agent_policies/WORKFLOW.md`, `agent_policies/TESTING.md`, and
`agent_policies/MONITORS.md`. The root `AGENTS.md` requires this routing.

This directory imposes an additional staged workflow in
`dirichlet_solve.md`:

1. Task 4 (`GreenReconstruction`) is the first unchecked task. Complete its
   theory, implementation, results, and `task4.md`, then stop for explicit user
   approval.
2. Do not begin Task 5 until Task 4 is approved.
3. Complete Task 5 (`TraceCorrected`), write `task5.md`, and stop for explicit
   approval.
4. Only after both approvals, run the best eligible formulation as an
   independently marched finite-wake follow-up. Add a separately identified
   follow-up task to the index before doing that work.

Never check a `User approved` box without explicit user instruction. Task item
files contain only final formulation, successful commands, artifacts,
measurements, conclusions, and theory corrections—not chronological failure
logs. Update `dirichlet_solve.md` with status and terse takeaways only.

The worktree is already dirty and contains substantial user changes, including
the formulation implementation and its tests. Inspect `git status --short` and
preserve all unrelated changes. In particular, do not assume the currently
untracked `src/FLOWPanel_formulation.jl`, `test/formulation_test.jl`,
`docs/wake_solve_schemes.md`, or this debug directory can be replaced or
cleaned.

## Exact physical case and established results

The user referred to `debug/dirichlet_solver/dirichlet_solver.md`, but the
actual path is `debug/dirichlet_solve/dirichlet_solve.md`.

Common geometry and flow:

- watertight rectangular NACA 0015 wing;
- chord `c = 1 ft = 0.3048 m`, span `b = 4 ft`, aspect ratio 4;
- `n_airfoil=161`, `n_span=13`, `n_endcap=9`, rounded end caps;
- 6,688 body panels and 13 paired shedding edges;
- `Uinf = 330.2 ft/s = 100.64496 m/s`, `rho = 1.225 kg/m^3`;
- primary angle of attack 3.94 degrees; stress case 45 degrees;
- `dt = 0.5c/Uinf`, constant incidence, no pitching;
- finite attached transition length `|Da|=0.05c` for the physical march;
- outputs under `data/dirichlet_solve/`, with 45-degree outputs under
  `data/dirichlet_solve/alpha45/`;
- compact existing comparison in each output directory's `comparison.csv`.

Lift is computed with steady `PressureBernoulli(rho;
correct_kuttacondition=false)` followed by `ForceMonitor` using
`WingNormalization(rho,Sref,c)`. Project the normalized force onto the lift
direction returned by `_lift_drag_span_directions(Uinf)`. Use `DirectBackend`
for authoritative pressure, physical field, and exterior-probe evaluation.

Established 3.94-degree results:

| Formulation/case | CL | Difference from semi-infinite |
|---|---:|---:|
| Semi-infinite attached wake | 0.2747643938323718 | reference |
| Old flat, `Da=0.5c`, total length 64c | 0.2742139824 | -0.2003% |
| Old flat, `Da=0.05c`, total length 63.55c | 0.2715662450 | -1.1640% |
| Old rolled march, settled step 80 (`tU/c=40`) | 0.265312060869 | -3.44016% |
| Direct-potential rolled, single shot | 0.2684046804 | -2.3146% |
| Direct-potential rolled, fixed-geometry iteration | 0.2733869039 | -0.5013% |

At 45 degrees, the old settled rolled wake is 9.32464% below its
semi-infinite baseline and the iterated direct-potential result is 7.4836%
below it. Read the exact stored values from the existing comparison CSV rather
than rounding these summaries in new output.

The semi-infinite reference residual at 3.94 degrees is about `2.48e-15`.
Task 2 established flat-wake length convergence and rolled-wake settling.
Task 3 is the explicit free-wake scalar-potential oracle.

## Existing formulation interfaces and meaning

Public exports are wired in `src/FLOWPanel.jl`. Implementations are in
`src/FLOWPanel_formulation.jl`; the simulation hook is in
`src/FLOWPanel_simulate.jl`; the affine wake-strength channel is in
`src/FLOWPanel_liftingbody.jl`; operator guards are in
`src/FLOWPanel_solver.jl`.

Use these existing constructors:

```julia
pnl.GreenReconstruction(
    gauge=:area_mean,
    recompute_interval=1,
    green_solver=nothing,
)

pnl.TraceCorrected(
    estimator=:green,              # or :line_integral
    gauge=:area_mean,
    quadrature=:trapezoid,         # or :simpson
    interior_path=:straight,       # or :deformed
    path_depth=1.0,
    recompute_interval=1,
    relaxation=1.0,
    green_solver=nothing,
)
```

Definitions:

- `VelocityThroughSources()` is the old/default route. It samples wake
  velocity at the body, converts it to body source strengths, and solves the
  existing finite attached-wake system.
- `GreenReconstruction` samples wake velocity, constructs wake-only source
  coefficients `sigma`, evaluates the *body source-panel* potential `S*sigma`,
  solves `(I-B)q=S*sigma` with a gauge constraint, and solves
  `G*muE=-S*sigma0-q`. It reconstructs a potential trace but does not directly
  evaluate the free wake's scalar potential at the body.
- `TraceCorrected(estimator=:green)` uses the same reconstructed `q`, sets
  `c=C*q`, solves `G*mutilde=-S*(sigma0+sigma)+W*c`, and uses physical
  circulation `gamma=C*mutilde-c`.
- `TraceCorrected(estimator=:line_integral)` obtains `c` solely by integrating
  wake-only velocity between paired trailing-edge control points. Trapezoid
  reuses endpoint velocity; Simpson adds batched midpoint probes; the deformed
  path routes through an upstream interior waypoint.
- No current `TraceCorrected` setting directly evaluates free-wake scalar
  potential. Task 3 supplies that oracle. Do not add a new public estimator
  merely to duplicate Task 3 unless later requested.

All alternative formulations currently require one finite-wake Dirichlet
source+doublet `RigidWakeBody`, paired shedding for gauge-independent traces,
a `Backslash` body solver, a wake system, and rigid body/`Da` geometry.
`TraceCorrected` additionally requires `backend_system isa DirectBackend`
because its affine attached-wake correction does not yet propagate through FMM
source buffers. `backend_wake` and `backend_solve` may remain FMM, but all
authoritative comparisons should use a Direct physical-system path. Add an old
`VelocityThroughSources`+Direct control so backend changes are not mistaken for
formulation effects.

For the 6,688-panel body, dense Green initialization holds approximately:

- `B`: 357.8 MB;
- bordered Green matrix: 357.9 MB;
- existing body `Backslash` LU: approximately 357.8 MB.

Peak major-matrix storage is therefore roughly 1.1 GB before temporaries. Run
full-resolution dense cases sequentially, drop state between cases, and allow
garbage collection before constructing the next Green state. Do not use
`gauge=:lsq` at full resolution; its dense QR is more expensive. Do not use
lagging or filtering in accuracy comparisons.

Distinguish immutable operator/factorization reuse from mutable solution state.
Because AoA is applied by rotating `Uinf` (`_uinf_from_alpha`) with a fixed body
mesh, the body-geometry-only operators are identical across **all** frozen wake
states and **both** angles and may be built once and reused:

- `B`/`(I−B)` and the bordered Green matrix plus its LU factorization
  (`gauge=:area_mean`) are body-mesh-only;
- `S` (body source-panel potential operator) is body-mesh-only — only the RHS
  vectors `sigma`, `sigma0` change per case.

The finite-wake body `Backslash` LU *does* depend on the attached-`Da` geometry,
so it is reusable only within a fixed (angle, `Da`) group, not globally. The
"no leak across cases" rule (see the harness section) governs mutable
solution/trace state only — `q`, `c`, `last_recompute`, and Krylov warm-start
must stay fresh per case; geometry-keyed factorizations may be shared. Reusing
them cuts the dominant dense-assembly/factorization cost several-fold without
weakening the immutable-state audit. Pin BLAS threads
(`LinearAlgebra.BLAS.set_num_threads`) to a fixed value so the recorded wall
times are comparable across cases.

## Required experiment matrix

For both 3.94 and 45 degrees, evaluate each formulation on these existing
finite-wake states:

1. terminal flat sequence with `Da=0.5c` and 64c total extent;
2. terminal flat sequence with `Da=0.05c` and 63.55c total extent;
3. selected settled rolled wake with `Da=0.05c` (normally step 80, but use the
   step selected by the existing Task 2 settling logic).

Reference rows in each output directory:

- Task 1 semi-infinite attached wake;
- Task 2 old `VelocityThroughSources` using its original backend;
- an additional old `VelocityThroughSources` Direct-system control;
- Task 3 direct-potential single-shot and iterated oracle.

Task 4 formulation:

- `GreenReconstruction(gauge=:area_mean, green_solver=nothing)`.

Task 5 formulations:

- `TraceCorrected(estimator=:green, gauge=:area_mean)`;
- line integral, straight path, trapezoid;
- line integral, straight path, Simpson;
- line integral, deformed path, Simpson, `path_depth=1.0`.

The straight/Simpson/deformed variants are part of the full matrix at both
angles because the user explicitly requested both angles in full. Treat their
spread as a quadrature/path-dependence diagnostic. The principal strictly
velocity-only candidate is the line-integral formulation.

Execute the matrix cheap-to-expensive so an interrupted run still yields a
usable ranking: all methods' 3.94-degree single-shot cases first (fast screen),
then the 3.94-degree fixed-geometry iterations, then the 45-degree cases, with
the rolled `Da=0.05c` state last within each angle. Make the resumable selector
granularity match this order (per method x angle x wake-state x route) so a
resume never reruns completed cells.

If the line-integral formulation is a serious final candidate, add a short
`path_depth` sensitivity check (`0.5, 1.0, 2.0`) on the primary rolled 3.94-degree
row only. It is cheap and directly feeds the quadrature/path-dependence
diagnostic; the deformed path is the near-singular-chord mitigation, so its
result may be `path_depth`-sensitive.

For every frozen case produce:

- a single-shot solve using the Task 2 active wake strengths;
- a fixed-geometry strength projection using the same scheme as Task 3.

The fixed-geometry iteration updates active free-wake strengths toward the
effective transition circulation while body and wake geometry remain fixed.
Use the existing Task 3 criteria: maximum active-row strength defect and
relative lift change both at or below `1e-8` for three consecutive iterations,
maximum 200 iterations, with the existing adaptive relaxation behavior.

A method is **eligible** iff it passes all three: (1) its formulation-specific
residual gate and finite-field checks; (2) the immutable-state audit; and (3)
the velocity-only constraint — it never requests `scalar_potential=true` from the
free-wake sources at body control points. Task 3 and any direct free-wake
potential route are ineligible; `GreenReconstruction` and `TraceCorrected`
(including `:green`) are eligible because their body source-panel potential and
`W*c` evaluations are explicitly permitted (see the Theory section).

After Task 4 and Task 5 are separately approved, select the eligible method
with the largest robust gap-closure fraction (`f_oracle`, per the ranking
definition above). March it independently from step
zero at both angles so its corrected shedding creates its own wake history.
Use the existing Task 2 settling rule: over the latest `10c/Uinf`, both
peak-to-peak lift divided by mean magnitude and least-squares linear change
divided by mean magnitude must not exceed 0.5%; assess first at `40c/Uinf` and
extend in `20c/Uinf` blocks through `80c/Uinf` if necessary.

## Harness implementation details

Extend `debug/dirichlet_solve/dirichlet_solve.jl`; reuse its Task 1-3 geometry,
state construction, snapshot, monitor, probe, CSV/TOML, VTK, and convergence
helpers rather than creating another case driver.

The existing frozen helper `_frozen_solve` assumes the old formulation and
must not be reused unchanged. For each copied body/wake/formulation:

1. Recalculate normals/control points as the existing harness does.
2. Construct or reuse the body `Backslash` solver only when it is valid for
   that exact body geometry/operator.
3. Call

   ```julia
   state = pnl.initialize_formulation(
       formulation, (body,), (wake,), solver,
       backend_solve, backend_system,
   )
   ```

4. Pass `formulation`, `formulation_state=state`, and the correct `i_step` into
   `pnl._steady_aerodynamics!`.
5. Use a fresh formulation state for every independent frozen case. Cached
   `q`, `c`, `last_recompute`, and Krylov objects must not leak across cases.
6. Evaluate pressure, lift, and exterior probes only after the formulation has
   restored physical mode.

For `TraceCorrected`, never compute circulation as `C*mu`. Obtain the effective
physical value via `_get_wakestrength_Gamma` for each shedding edge, or exactly
equivalently as `C*mutilde-state.c`. The first shed row, transition/free-row
mismatch, VTK circulation, pressure, and exterior fields must all see this
corrected circulation.

Snapshot coverage must include:

- body nodes, cells, and shedding maps;
- wake nodes, connectivity/options, active row count, active and inactive
  strength storage as appropriate for the current inner solve;
- body `Da` arrays;
- `wake_strength_correction`, `wake_strength_shift`, and
  `wake_correction_active`;
- source and doublet strengths, explicitly allowing only the arrays named by
  the current solve/outer iteration contract to change.

Use formulation-specific residuals:

- Green reconstruction:
  `(I-B)q + lambda*a - S*sigma`, the gauge row `dot(a,q)`, and
  `G*muE + S*sigma0 + q`.
- Trace corrected:
  `G*mutilde + S*(sigma0+sigma) - W*c` and
  `gamma - (C*mutilde-c)`.

Require relative residuals no larger than `1e-10`, finite fields, and exact
immutable-state checks. For paired edges, compare area-mean and direct-trace
gauges only after removing the direct trace's area-weighted mean.

Diagnostic comparisons to record:

- `CL`, signed and relative difference from Task 1;
- gap-closure fraction relative to the matching old finite-wake case;
- difference from Task 3 single-shot and iterated results;
- norms of `q`, `c=Cq`, and effective `gamma`;
- Green `q` versus Task 3 direct wake trace after gauge matching;
- each trace estimator's `c` versus Task 3's direct trailing-edge
  potential-difference oracle;
- trapezoid-versus-Simpson and straight-versus-deformed `c` spreads;
- transition/free-row mismatch using effective `gamma`;
- relative exterior velocity difference at the existing 16 probes versus
  Tasks 1, 2, and 3;
- initialization, solve/iteration, and total wall times;
- backend and all formulation settings.

Extend `_COMPARISON_COLUMNS` without dropping or silently reinterpreting old
columns. Use formulation+case+route identifiers that cannot overwrite Task
1-3 rows. Write per-task convergence/iteration CSVs, invariant TOML, metadata
TOML, and terminal body/wake VTK for single-shot and converged iterated results.

## Theory and conclusions

Before Task 4 implementation, compare the existing derivation in
`docs/wake_solve_schemes.md` against the current code. Correct the theory first
if signs, gauges, operator definitions, or physical circulation differ. Do the
same for Task 5 before its implementation work.

Sign conventions are high risk: FLOWPanel mesh normals use the GeometricTools
right-hand convention while panel kernels include explicit Hess/Smith and
Katz/Plotkin sign adjustments. Do not change a sign merely to improve `CL`.

When reporting a conclusion, perform the repository-required robustness pass:

1. state the evidence-based ranking;
2. recheck methods, stored results, residuals, frozen-state invariants,
   backend controls, and conclusion for consistency;
3. revise the conclusion if the audit changes it;
4. explicitly distinguish estimator error, rolled-wake geometry/strength
   feedback, backend effects, and quadrature/path sensitivity.

The winning method must not directly request `scalar_potential=true` from the
free-wake sources at body control points. Body source-panel potential
evaluation for the Green identity and attached-transition `W*c` evaluation do
not violate this requirement. Task 3 direct potential is diagnostic-only and
cannot be selected as the final method.

## Verification commands and current test state

Run the narrow mechanical test before physical cases:

```bash
julia --project test/formulation_test.jl
```

It was run successfully on 2026-07-20. All seven stages passed. Reported small
case diagnostics included:

- trapezoid `c` versus direct TE-potential oracle: relative L2 `3.47e-5`;
- Simpson: `3.33e-10`;
- deformed Simpson: `9.07e-10`;
- Green `q` versus direct `q`: relative L2 `0.0748`;
- `TraceCorrected(:green)` versus `GreenReconstruction` circulation gap:
  `3.15e-13`;
- direct/FMM Krylov Green trace agreement: `1.18e-11`;
- FGS versus LU Green trace agreement: `3.44e-15`.

The test took about five minutes, mostly in matrix-free route checks. It is a
standalone diagnostic and is not currently included by `test/runtests.jl`.
Because it is the only coverage of `src/FLOWPanel_formulation.jl`, run it (or a
fast subset) explicitly whenever a theory correction touches that source — the
narrow maintained suites below do not exercise the formulation code.

After any production-source changes, run the narrow maintained suites first:

```bash
julia --project -e 'include("test/runtests_unit_solver.jl")'
julia --project -e 'include("test/runtests_unit_liftingbody.jl")'
julia --project -e 'include("test/runtests_unit_wake.jl")'
julia --project -e 'include("test/runtests_unit_simulate.jl")'
julia --project -e 'include("test/runtests_unit_postprocess.jl")'
```

Then run the broad suite if production code changed:

```bash
julia --project -e 'include("test/runtests.jl")'
```

Use the debug script's CLI pattern for physical runs. Existing selectors cover
Tasks 1-3; add explicit Task 4 and Task 5 selectors plus narrow per-case
selectors so interrupted expensive work can resume without rerunning all
cases. Preserve the existing `--alpha-deg` and `--output-dir` behavior.

## Completion deliverables

- Correct theory in `docs/wake_solve_schemes.md` if required.
- Extended `debug/dirichlet_solve/dirichlet_solve.jl` with resumable Task 4/5
  selectors and formulation-aware diagnostics.
- Final `task4.md`, followed by user approval.
- Final `task5.md`, followed by user approval.
- Updated index takeaways and approval state only when explicitly authorized.
- Comparison/convergence/invariant/metadata artifacts at both angles.
- A separately staged self-consistent march of the winning eligible method at
  both angles.
- Final audited conclusion ranking all methods by gap closure and explaining
  whether remaining error comes primarily from trace estimation, discrete
  Green representation, fixed rolled-wake geometry, or coupled wake evolution.
