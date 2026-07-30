# Phase 2b — Formulation/Topology Attribution

## Purpose and required context

Determine why Dirichlet solves produce lower integrated bound circulation than
Neumann solves in the Phase 1 DJI 9443 steady cases.

Phase 2 accepted the sharp capped/Dirichlet trailing edge for the tested
local-failure mode. This phase does not repeat that adequacy study. Instead,
it separates formulation, topology, and Kutta/wake-coupling effects using a
simple controllable oracle first, then applies only the necessary attribution
tests to the DJI meshes.

Read the repository instructions, the dashboard, the top snapshot in the
[Phase 2b log](../../logs/dji_convergence_20260722/phase_02b_formulation_attribution.md),
and:

- `agent_policies/WORKFLOW.md`
- `agent_policies/TESTING.md`

Do not begin unless Ryan explicitly approves Phase 2b. Do not begin Phase 3
until this phase is complete or Ryan explicitly cancels it.

## Evidence entering this phase

- Phase 1 found uncapped/Neumann integrated circulation higher than
  capped/Dirichlet by `+6.519%` at 40-series and `+4.853%` at 57-series.
- Phase 1 also recorded that this comparison changes both topology and
  formulation, so it cannot identify the cause.
- Phase 2 found no sharp-TE operator-fragility trigger: the largest tested
  integrated perturbation response was `0.084%`, the largest outboard response
  was `0.232%`, `kerneloffset_panel` was inert at reported precision, and
  off-collocation residuals were not TE-local outliers.
- Phase 1 refinement asymmetry: capped/Dirichlet integrated circulation moved
  only `+0.255%` from 40- to 57-series while uncapped/Neumann moved `-1.314%`
  *toward* Dirichlet. "Neumann under-resolved" is therefore a live hypothesis
  and the refinement trend is a primary diagnostic, not a fallback.
- The repository already investigated a "Morino formulation velocity-only
  discrepancy" (commit `be46062`): `debug/dirichlet_solve/` (task1–task5,
  `green_identity_analysis.md`, `unregularized_velocity_plan.md`; CL
  differences up to `-3.4%` between wake routes) and
  `debug/debug_bc_comparison.jl`, a working Dirichlet-vs-Neumann harness that
  includes the drop-one-panel trick for a full-rank closed Neumann solve.
  Phase 2b must intake and reuse this evidence and code before building new
  diagnostics.
- Solver naming fact: there is no `BackslashNeumann`/`BackslashDirichlet`.
  One `Backslash` solver (`src/FLOWPanel_solver.jl`) dispatches on the body's
  `DBC` type parameter. Neumann = `VortexRing` kernel, `DBC=false`,
  `watertight=false`, solved strength column 1, RHS `-U·n`. Dirichlet =
  `Union{ConstantSource,VortexRing}` (or `ConstantDoublet`), `DBC=true`,
  `watertight=true`, solved strength column 2, RHS `-self.potential`
  (Morino interior-potential nullification). The affine Kutta-correction
  channel `γ = μ_upper − μ_lower − c` exists via `set_wake_correction!`
  (`src/FLOWPanel_liftingbody.jl:1628`) and is activated by the
  `TraceCorrected` solve formulation (`src/FLOWPanel_formulation.jl`).

## Phase ordering

Lead with the two near-free steps (prior-evidence intake and the Gate A0
formulation sweep), then the simplified oracle, and answer the
field-circulation question early: do Dirichlet and Neumann produce different
exterior velocity circulation, or only different extracted
`Gamma_TE = mu_upper - mu_lower`? Use those answers to gate all expensive
diagnostics.

1. Prior-evidence intake (read-only, no solves): read the
   `debug/dirichlet_solve/` conclusions (`task5.md`,
   `green_identity_analysis.md`, `unregularized_velocity_plan.md`,
   `dirichlet_solve.md`) and `debug/debug_bc_comparison.jl`. Record in the
   Phase 2b log which mechanisms that work already implicates or excludes,
   and which of its tooling (drop-one-panel Neumann build, Green identity
   residual, source/doublet velocity split) Phase 2b will reuse.
2. Gate A0 formulation sweep on an unchanged capped/Dirichlet body: solve the
   same body and mesh with `formulation = VelocityThroughSources` (current
   baseline), `TraceCorrected` (both `:green` and `:line_integral`
   estimators), and `GreenReconstruction`, then compare extracted `∫Γ dr` on
   the fixed station bins. `TraceCorrected` applies the affine correction
   `γ = μ_upper − μ_lower − c`, which shifts circulation broadly across the
   span — the same spatial shape as the observed gap (integrated and outboard
   percentages nearly equal in Phase 1). If `TraceCorrected` moves Dirichlet
   circulation toward Neumann by order `5%`, attribute the gap to
   wake-potential/Kutta-trace handling and short-circuit the later branches
   down to confirmation-only. Run this first on the pitching-wing oracle body
   (cheap), then on DJI `new40c` if local runtime permits. Constraints:
   `TraceCorrected` requires `DirectBackend`; `steady!` does **not** currently
   accept `formulation=` (only `simulate!` does,
   `src/FLOWPanel_simulate.jl:332`), so Gate A0 needs either a small
   `steady!` wiring change or an equivalent one-step frozen-wake `simulate!`
   invocation — treat that wiring as a driver task.
3. Implement a Phase 2b driver, likely
   `examples/dji9443_formulation_attribution.jl`, with modes:
   `smoke`, `formsweep`, `oracle`, `thickness`, `analyze`, and optional
   `dji_bridge`.
4. Run the simple generated mesh oracle as a default three-point refinement
   ladder (coarse/medium/fine). These are small dense solves; the refinement
   trend of the Dirichlet-Neumann gap is the single most diagnostic curve and
   must not be gated. Refinement gating applies only to the expensive DJI
   bridge cases.
5. Apply Gate A: exterior velocity equivalence, loop-integral circulation, and
   Kutta map consistency on shared probe/station geometry.
6. Branch targeted follow-ups:
   - representation/extraction branch: Dirichlet mix/nonuniqueness diagnostics;
   - physical/formulation branch: wake/Kutta decomposition and
     forcing-versus-operator attribution;
   - TE/closure branch: gated thickness and trailing-edge closure screen.
7. Expand the thickness/refinement matrix only if the screen identifies
   thickness, TE closure, or discretization dependence.
8. DJI 2x2 bridge only if the oracle does not already isolate the mechanism,
   or if DJI-specific confirmation is needed before Phase 3.

Do not perform the geometry-only aerodynamic audit from the earlier discussion
in this phase.

## Driver requirements

Build the Phase 2b driver as a reusable diagnostic layer around the existing
oracle, reusing the Phase 2 machinery in `examples/dji9443_te_adequacy.jl`
rather than creating parallel one-off diagnostics. Also reuse
`debug/debug_bc_comparison.jl` (Neumann body construction including
drop-one-panel, source/doublet velocity split) and `test/formulation_test.jl`
(which already builds the pitching-wing oracle body via
`build_pitching_wing_body` and exercises `_get_wakestrength_Gamma`,
Green-identity, and Kutta-map patterns).

Core driver requirements (needed for smoke/oracle/Gate A0/Gate A):

- full and body-only operator assembly, including the Phase 2 `G`/`B` split
  where applicable;
- fixed-bin circulation extraction and common-grid interpolation;
- exterior velocity probe and loop-integral circulation extraction on shared
  probe geometry, including the loop-integral self-check defined under
  Common observables;
- formulation-sweep support for Gate A0 (`VelocityThroughSources`,
  `TraceCorrected`, `GreenReconstruction` on an unchanged body);
- off-collocation probe generation and residual CSV writers;
- Kutta map consistency comparing `Gamma_TE`, `_get_wakestrength_Gamma`, and
  loop-integral circulation on identical shedding edges.

Conditional tooling — implement only when the corresponding branch is
triggered, not up front: row norm, row similarity, adjoint sensitivity,
low-rank perturbation response, and exterior-invisible mix-mode tools (the
Phase 2 adjoint/SMW helpers in `dji9443_te_adequacy.jl` are the starting
point).

Required driver modes:

- `smoke`: tiny synthetic case that exercises CSV writers, fixed-grid
  interpolation, shared probe generation, loop-integral machinery,
  wake-strength extraction, basic residual calculations, decomposition algebra,
  and analysis inputs without requiring DJI assets or large dense solves;
- `formsweep`: Gate A0 formulation sweep (`VelocityThroughSources`,
  `TraceCorrected` with `:green` and `:line_integral`, `GreenReconstruction`)
  on an unchanged capped/Dirichlet body with `DirectBackend`, reporting the
  fixed-bin circulation for each formulation;
- `oracle`: generated open/capped oracle cases at the requested refinement,
  including the mandatory Gate A field-circulation outputs and the
  capped+Neumann drop-one-panel tiebreaker;
- `thickness`: gated global-thickness and trailing-edge closure screen using
  the same generated oracle, station grid, probe families, and attribution
  outputs selected by Gate A;
- `analyze`: aggregate oracle CSVs into tables and a short report with the
  attribution decision;
- `dji_bridge`: optional integration mode for DJI bridge cases only after the
  oracle identifies the likely mechanism or DJI confirmation is explicitly
  needed.

Write Phase 2b artifacts under `data/dji_convergence_20260722/phase_02b_formulation_attribution/`.
All CSVs must carry enough metadata columns to identify refinement, topology,
formulation, gauge convention, station grid, extraction convention, probe
family, operator/probe domain, and whether the row/operator came from the full
or body-only system.

## Common observables

All Dirichlet and Neumann comparisons must use one fixed station grid and one
extraction convention. Do not compare only common circulation totals. Report the
following matched observables for every oracle case and every DJI bridge case
that is run:

- TE circulation jump, `Gamma_TE(r/R)` or finite-wing spanwise equivalent;
- integrated `∫Γ dr`;
- thrust-weighted `∫Ω r Γ dr` for rotor-like cases, or the documented
  finite-wing weighted analogue if the oracle is wing-like;
- outboard integral using the same outboard station mask in all formulations;
- blade/side symmetry for rotor cases, or finite-wing side symmetry for
  wing-like cases.

Stationwise differences must also be carried through to the decomposition: the
analysis should report which stations contribute to the integrated
Dirichlet-Neumann gap, not just the final integral difference.

Treat `Gamma_TE = mu_upper - mu_lower` as a formulation-dependent extracted
quantity until the diagnostics confirm that it matches the velocity field.
Where feasible, also report velocity-based circulation from shared exterior
loop integrals around selected stations and compare it with `Gamma_TE`,
attached-wake strength, and the Neumann circulation unknown.

Loop-integral procedure (mandatory definition, so Gate A is trusted on the
first pass): at each selected station, integrate `∮ V·dl` around a closed
circuit lying in the spanwise-station plane, enclosing the airfoil section,
and crossing the attached wake sheet exactly once so the integral measures
bound circulation at that station. Keep all probe points offset from both the
body surface and the wake sheet (respecting `kerneloffset_targets`), and use
identical circuits for every formulation. Self-check before any
cross-formulation comparison is drawn: each loop integral must be evaluated at
two loop radii and two quadrature densities, and the four values must agree to
a stated tolerance (default `0.1%` relative); loops failing the self-check are
excluded and reported. Note that off-body velocity evaluation goes through the
induced-kernel path (`influence!`/`_Uind`) and never touches the on-surface
`grad_mu` reconstruction, so Gate A probes are clean of the known
triangulation-sensitive surface-velocity issue — do not spend effort on
`grad_mu` options for field probes.

## Simple generated mesh oracle

Generate an easy-to-refine lifting surface with controlled open and capped
variants. Prefer a rectangular or mildly tapered finite wing/rotor-like blade
with simple sections, fixed angle of attack or equivalent rotating steady
setup, and a trailing edge that can be opened or capped without changing the
outer surface coordinates away from the closure strip.

Use the procedural pitching-wing mesh generator from
`examples/pitching_wing.jl` as the preferred oracle starting point. Adapt it
for Phase 2b rather than duplicating a separate generator:

- split the current construction so the shared extruded lifting-surface mesh
  can be returned before cap cells are appended;
- add a topology option, for example `caps=true/false`, where `caps=false`
  returns the open lifting surface with `watertight=false` and `caps=true`
  appends the existing rounded finite-wing caps with `watertight=true`;
- keep the open and capped variants on identical lifting-surface coordinates
  and differ only by appended cap nodes/cells;
- keep the existing refinement controls (`n_span`, `n_airfoil`, `n_endcap`) and
  use them to define the three rungs of the default refinement ladder.

Make the NACA 4-series section configurable. The current NACA 0015 behavior can
remain the default, but the oracle generator should accept an airfoil option
such as `naca="0015"` or an explicit thickness/camber tuple so Phase 2b can
repeat the same attribution checks at different thicknesses. At minimum, expose
the symmetric thickness parameter needed for NACA 00xx sections; if cambered
4-series support is added, record camber, camber location, and thickness in all
oracle metadata.

Generate cases from the same outer lifting-surface coordinates and run matched
steady, no-free-wake solves at every rung of the default three-point
refinement ladder:

- open + Neumann vortex-only, `watertight=false` (production-like Neumann);
- capped + Dirichlet source/vortex, `watertight=true` (production-like
  Dirichlet);
- capped + Neumann vortex-only with drop-one-panel regularization, the
  designated topology-versus-formulation tiebreaker. Drop one mid-mesh panel
  (away from caps and TE) so the pure-VortexRing closed-surface system is full
  rank, following `debug/debug_bc_comparison.jl:27–34`. Record the dropped
  panel, rank/conditioning, and sensitivity of the observables to the choice
  of dropped panel (one alternate drop location at the coarse rung).

Do **not** run open + Dirichlet by default: the Morino interior-potential
boundary condition has no interior region on an open surface, so the case is
ill-posed and its output cannot be interpreted. Only revisit it if a specific
branch produces a concrete, written justification.

Confound-resolution logic: if open+Neumann and capped+Dirichlet exterior
fields differ, that alone does not implicate the formulation — the caps
physically change the flow. The tiebreaker resolves it: if the
capped+Neumann(drop-one-panel) field matches capped+Dirichlet, the gap is
topology/closure-dominant; if it matches open+Neumann, the gap is
formulation-dominant. The primary evidence is the matched observable set
above, not agreement or disagreement in circulation totals alone.

Measure:

- all common observables from the fixed station grid;
- mesh-refinement trend of the Dirichlet-Neumann gap across the default
  three-point ladder, reported per formulation (Phase 1 showed Neumann moving
  toward converged Dirichlet under refinement, so the per-formulation trend is
  itself attribution evidence);
- solver residuals, condition estimates where available, and blade/side
  symmetry if applicable;
- exterior velocity agreement and loop-integral circulation agreement.

The refinement ladder is not gated: oracle solves are small and dense, and the
trend is a primary diagnostic. Gating applies to broad thickness sweeps and
DJI bridge cases, which remain conditional on Gate A and the ladder result.

Decision from the oracle:

- If the Dirichlet-Neumann gap appears on the simple matched geometry and
  persists in exterior velocity and loop-integral circulation, treat a
  physical formulation or Kutta/wake-coupling difference as the leading cause.
- If the gap disappears on matched simple geometry, treat DJI topology,
  closure, or mesh-specific coupling as the leading cause and continue to the
  DJI bridge.
- If the gap depends strongly on refinement, classify the formulation
  comparison as numerically unresolved and identify the refinement dimension.
- If exterior velocity and loop-integral circulation agree while `Gamma_TE`
  differs, treat the gap as a Dirichlet representation/extraction issue and
  identify the weakly observable mix direction.

## Gated Thickness/TE Closure Screen

Run thickness and trailing-edge closure cases only when Gate A or the oracle
refinement trend points to TE/closure localization, global thickness
sensitivity, or unresolved discretization. The first pass is a cheap screen,
not a full sweep.

Use the same pitching-wing oracle generator, station grid, extraction
convention, residual families, and attribution outputs selected by Gate A.
The initial screen must include:

- a thin and baseline symmetric NACA 00xx section, for example `0006` and
  `0012` or `0015`, at one refinement;
- sharp or near-zero TE and one blunt/finite-TE variant at fixed global
  thickness;
- the same Gate A exterior velocity, loop-integral circulation, and Kutta map
  outputs as the base oracle.

Expand to a fuller sweep, for example `0006`, `0012`, and `0018`, multiple TE
closure variants, and two refinements for the thinnest and baseline cases only
if the screen shows material dependence on global thickness, local TE
thickness, or refinement.

Separate these two axes in the mesh generator and metadata. Do not infer a
trailing-edge pathology from a global thickness change alone. For the
trailing-edge closure sweep, preserve the suction and pressure surfaces away
from the closure neighborhood as much as the procedural section definition
allows.

Report the following diagnostics for every case that is actually run:

- Dirichlet-Neumann TE circulation jump and integrated gap versus global
  `t/c`;
- the same gap versus local TE thickness normalized by nearby panel length;
- exterior velocity agreement and loop-integral circulation agreement versus
  `t/c` and local TE thickness;
- Kutta map consistency versus `t/c` and local TE thickness;
- stationwise contribution to the integrated gap.

Add branch-specific diagnostics only when needed:

- no-penetration leakage when exterior velocity differs materially;
- Dirichlet potential-trace or Green identity residuals only for closed/capped
  cases and only as consistency checks;
- row norms, row similarities, adjoint response, or low-rank response when the
  physical/formulation branch needs operator attribution;
- exterior-invisible mix-mode strength when the representation/extraction
  branch needs weak-observability evidence.

Interpretation rules:

- If the gap follows global `t/c` at fixed TE closure, classify the mechanism as
  global section-thickness/operator-conditioning sensitive.
- If the gap follows TE closure at fixed global `t/c`, classify the mechanism as
  local trailing-edge/closure/Kutta sensitive.
- If the gap collapses primarily with refinement, classify it as numerically
  unresolved and identify the controlling refinement dimension.
- If thick trailing-edge cases improve Dirichlet-Neumann circulation agreement,
  report whether the improvement comes from forcing, operator rows, wake/Kutta
  coupling, or closure topology using the forcing-versus-operator diagnostics.
- If thick trailing-edge cases only improve `Gamma_TE` agreement while exterior
  velocity and loop-integral circulation were already matched, classify the
  improvement as extraction/mix conditioning rather than a physical circulation
  change.

## Matched-observable residual audit

Replace the vague cross-residual audit with mandatory Gate A residuals and
conditional follow-up residuals.
For each oracle solve, evaluate the solved state on shared off-surface probe
points and shared station bins. Use identical probe locations across
formulations wherever the same outer surface exists.

Mandatory Gate A residuals:

1. Exterior velocity residual:
   evaluate Dirichlet and Neumann solved states on shared exterior probe clouds
   near the surface, near the TE/wake, and downstream. Report vector velocity
   error by station and region, plus loop-integral circulation around selected
   stations.
2. Kutta trace/wake-strength residual:
   report the residual at each shedding edge, including wake strength, bound
   TE jump, and the mismatch between the solved wake/Kutta trace and the
   velocity-based loop circulation.

Conditional residuals:

1. No-penetration leakage:
   evaluate each solved state on shared exterior and interior off-surface probe
   points, then report normal-velocity leakage by station and by region
   (`TE`, `mid_chord`, `root`, `tip`, and `closure_adjacent` where applicable).
2. Dirichlet potential-trace residual:
   use a fixed gauge convention and body-only `B` assembly only on closed/capped
   cases. Do not compare raw Dirichlet and Neumann potential traces as primary
   evidence.

Compare residual distributions against the local difference in circulation.
The useful result is spatial attribution: whether the lower Dirichlet
circulation is associated with global boundary-condition mismatch, trailing-edge
rows, root/tip closure, or wake coupling.

## Dirichlet mix/nonuniqueness diagnostics

Add explicit tests for the fact that a closed-surface Dirichlet source/doublet
representation can have many source/doublet mixes that reproduce nearly the
same exterior velocity field. Run these diagnostics only after Gate A shows
that exterior velocity and loop-integral circulation agree while
`Gamma_TE = mu_upper - mu_lower` differs materially.

Gauge scoping fact: a global constant gauge mode cancels exactly in
`mu_upper - mu_lower`, so the constant nullspace of `(I-B)` cannot by itself
explain a `Gamma_TE` gap. Verify the constant-mode cancellation once as an
algebra check, then restrict the mode search to non-constant weakly observable
directions. Any candidate mix mode must be reported with its projection onto
the constant mode removed. Note also that the `TraceCorrected` affine
correction `c` is precisely a per-edge shift of `Gamma_TE` that the plain
solve omits — Gate A0 tests it directly and its result should be folded into
this branch's interpretation.

First run the following compact sequence on the relevant closed/capped oracle
cases:

1. Exterior-invisible circulation modes:
   assemble a probe operator mapping Dirichlet source/doublet perturbations to
   exterior velocity. Use SVD or constrained least squares to find perturbation
   directions with small exterior velocity penalty and large `C*mu` or
   integrated `Gamma_TE` response.
2. Green identity closure:
   compute body-only `S*sigma` and `(I-B)*mu` from the Dirichlet solved state,
   excluding attached-wake columns. Apply the common gauge and report the
   stationwise residual, especially near TE and closure-adjacent panels. Treat
   this as a closed-surface consistency check, not as cross-formulation primary
   evidence.
3. Source-only and doublet-only velocity audit:
   split the Dirichlet exterior tangential velocity and loop-integral
   circulation into source and doublet/vortex contributions.

Then run targeted follow-ups when the compact sequence indicates a likely
representation/mix issue but not enough evidence to classify it:

- Mix perturbation sensitivity:
  add small exterior-invisible perturbations to the solved Dirichlet state and
  recompute `Gamma_TE`, integrated circulation, exterior velocity residual,
  interior potential residual, and Kutta wake strength.
- Neumann-to-Dirichlet reconstruction:
  sample the Neumann exterior velocity field on Dirichlet centroids, compute
  Dirichlet source strengths from the normal velocity, solve/reconstruct a
  compatible Dirichlet trace, and compare the resulting `C*mu` with the
  Neumann circulation.
- Dirichlet-to-Neumann boundary replay:
  use the Dirichlet exterior velocity on Neumann control points as the Neumann
  RHS and solve the Neumann system. If Neumann recovers high circulation from
  Dirichlet velocity, classify the original gap as extraction/mix sensitive.
- Body-only versus wake-coupled mix modes:
  repeat the exterior-invisible mode search with body-only `B` and full
  `G = B + W*C`. If the weak mode appears only in the full operator, attribute
  it to attached-wake/Kutta coupling.
- Region-restricted mix modes:
  constrain perturbations to TE-adjacent, mid-chord, root/tip, and
  closure-adjacent panels separately, then report which region can move
  integrated circulation with the smallest exterior velocity penalty.

Interpretation rules:

- If exterior velocity and loop-integral circulation match but
  `mu_upper - mu_lower` differs, classify the gap as a Dirichlet
  representation/extraction issue.
- If exterior velocity differs in the same stations as `Gamma_TE`, classify the
  gap as a physical/formulation difference and continue with wake/Kutta and
  forcing-versus-operator attribution.
- If a near-null perturbation can move integrated `Gamma_TE` by order `5%`
  while keeping exterior velocity residual small, report that the Dirichlet
  extracted circulation is weakly observable and quantify the mode.
- If the perturbation cost grows materially under refinement, classify the
  effect as discretization-sensitive; if it remains small, carry it forward as
  a formulation/extraction limitation.

## Wake/Kutta contribution decomposition

Run this branch when Gate A shows a material exterior velocity or
loop-circulation difference. Decompose the solve into body-only and
attached-wake/Kutta parts using the same conceptual split as Phase 2, adapted to
each formulation:

- full solve circulation;
- body-only solve with the attached wake suppressed;
- attached-wake delta, `full - body_only`;
- local TE-row contribution;
- root/tip or closure-adjacent contribution;
- stationwise contribution to the integrated Dirichlet-Neumann gap.
- body-only and wake-coupled exterior-invisible mix-mode strength where the
  Dirichlet mix diagnostics indicate weak observability.

If feasible, repeat the decomposition on the 40-series DJI capped/Dirichlet and
uncapped/Neumann cases only after the oracle has identified the most likely
terms to inspect.

## Forcing-versus-operator attribution

Run this branch after the wake/Kutta decomposition if the phase still needs to
distinguish changes caused by different RHS/source terms from changes caused by
operator rows/columns.

Report:

- RHS/source-term magnitudes and spatial distributions for each generated case;
- full and body-only operator row norms near the TE, mid-chord, root, and tip;
- row similarities between matched open/capped and Dirichlet/Neumann rows in
  the same regions;
- adjoint or low-rank response estimates identifying which rows or columns can
  move integrated circulation by the observed `~5%` DJI gap;
- source-only versus doublet-only exterior velocity contributions near the TE
  when the loop-integral diagnostic shows a `Gamma_TE`/velocity-circulation
  mismatch.

Defer Green-system RHS projection error to Phase 4 unless Phase 2b has entered
the representation/reconstruction branch and the compact mix diagnostics remain
inconclusive.

The analysis should explicitly say whether the observed gap is plausible from
forcing changes, operator changes, wake/Kutta coupling, closure topology, or a
mixed mechanism.

## Conditional DJI 2x2 bridge

Run this only if the oracle and residual/decomposition checks do not provide a
clear enough attribution for the Phase 3 handoff.

Construct diagnostic DJI bridge cases that change one factor at a time as
cleanly as possible:

- current capped + Dirichlet baseline;
- current uncapped + Neumann high-circulation case;
- capped-derived open + Neumann, preserving outer lifting-surface coordinates
  as much as possible;
- uncapped-derived capped + Dirichlet, using a minimal closure strip when clean
  meshing permits.

These are attribution cases, not production candidates. Preserve the Phase 1
safe shedding pattern: build the no-shedding body first, derive shedding from
the constructed rewound cells, then rebuild with shedding.

Measure the same fixed-grid circulation metrics used in Phase 1, plus Gate A
outputs (`exterior velocity equivalence`, `loop-integral circulation`, and
`Kutta map consistency`). Add Green identity closure or mix diagnostics to the
bridge only if the oracle classified the issue as representation/extraction or
left that branch unresolved. The bridge answer should identify whether the
`~5%` DJI gap follows cap topology, boundary-condition formulation, Kutta/wake
coupling, Dirichlet representation/extraction, or remains unresolved.

## Test plan

Step 0 (no code, no solves): perform the prior-evidence intake and write its
summary into the Phase 2b log before any driver work.

Step 1 (Gate A0): run the formulation sweep on the pitching-wing
capped/Dirichlet body:

```bash
PHASE2B_MODE=formsweep julia --project examples/dji9443_formulation_attribution.jl
```

with `formulation ∈ {VelocityThroughSources, TraceCorrected(:green),
TraceCorrected(:line_integral), GreenReconstruction}` on an unchanged body and
`DirectBackend`. Extend to DJI `new40c` only if the oracle-body result is
material and local runtime permits.

Then run the new driver smoke mode:

```bash
PHASE2B_MODE=smoke julia --project examples/dji9443_formulation_attribution.jl
```

This must validate CSV writers, common-grid interpolation, residual
calculations, exterior velocity probes, loop-integral circulation,
wake-strength extraction, and decomposition algebra on a tiny synthetic case.
Do not require SVD/null-mode, reconstruction replay, or full thickness-sweep
machinery in smoke.

Then run the generated-oracle Gate A ladder locally with direct dense solves:

```bash
for r in coarse medium fine; do
  PHASE2B_MODE=oracle PHASE2B_REFINEMENT=$r julia --project examples/dji9443_formulation_attribution.jl
done
PHASE2B_MODE=analyze julia --project examples/dji9443_formulation_attribution.jl
```

Run thickness cases only after Gate A points to TE/closure, thickness, or
refinement sensitivity:

```bash
PHASE2B_MODE=thickness PHASE2B_REFINEMENT=<tag> PHASE2B_THICKNESS_SCREEN=true julia --project examples/dji9443_formulation_attribution.jl
PHASE2B_MODE=analyze julia --project examples/dji9443_formulation_attribution.jl
```

For any solver-path change or helper extraction added under `test/`, run:

```bash
julia --project -e 'include("test/runtests_unit_solver.jl")'
julia --project -e 'include("test/runtests_unit_postprocess.jl")'
```

Treat DJI bridge runs as conditional integration checks only after the oracle
identifies the likely mechanism.

When local runtime becomes limiting, prioritize Phase 2b diagnostics in this
order:

1. prior-evidence intake (free);
2. Gate A0 formulation sweep on the oracle body;
3. exterior velocity equivalence and loop-integral circulation (with the
   loop self-check);
4. Kutta map consistency;
5. the capped+Neumann drop-one-panel tiebreaker and the remaining rungs of
   the refinement ladder;
6. Gate A0 on DJI `new40c`;
7. branch-specific wake/Kutta or mix diagnostics;
8. gated thickness/TE closure screen;
9. conditional DJI bridge.

## Decision and exit gate

The phase succeeds if it gives a reviewed attribution for the lower integrated
Dirichlet circulation:

- formulation-dominant;
- topology/closure-dominant;
- attached-wake/Kutta-coupling dominant;
- Dirichlet representation/extraction dominant;
- mesh/refinement unresolved;
- or mixed, with quantified contributions.

If no attribution is defensible, carry the Dirichlet-Neumann gap into Phase 3
and Phase 5 as an unresolved formulation factor with explicit required
comparisons.

Deliverables:

- prior-evidence intake summary in the Phase 2b log (mechanisms implicated or
  excluded by `debug/dirichlet_solve/` and `debug/debug_bc_comparison.jl`);
- Gate A0 formulation-sweep table (circulation per formulation, oracle body
  and any DJI case run) and its attribution consequence;
- Phase 2b driver with `smoke`, `formsweep`, `oracle`, `thickness`,
  `analyze`, and optional `dji_bridge` modes;
- generated-oracle setup description and refinement-ladder table;
- Gate A field-circulation table (with loop-integral self-check results) and
  decision, including the capped+Neumann tiebreaker reading;
- thickness/TE closure screen table and interpretation if run;
- fixed-grid circulation metrics for all oracle cases;
- matched-observable residual summary for mandatory residual families and any
  conditional residual families that were run;
- Dirichlet mix/nonuniqueness diagnostic summary if the
  representation/extraction branch is triggered;
- wake/Kutta decomposition summary if the physical/formulation branch is
  triggered;
- forcing-versus-operator attribution summary if needed after wake/Kutta
  decomposition;
- conditional DJI bridge metrics if run;
- reviewed conclusion and handoff decision for Phase 3.

Update this phase log and the dashboard, then stop for explicit Phase 3
approval.
