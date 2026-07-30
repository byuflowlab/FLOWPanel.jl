# Task 4 follow-up — unregularized wake-velocity Green-identity test

## Objective and decision rule

Determine whether the finite-wake Green-identity residual is materially caused
by using a regularized doublet-panel **velocity** together with an unregularized
analytic doublet-panel **potential**. Repeat the existing frozen-state and mesh-
refinement residual diagnostics with the velocity core removed, changing no
geometry, strengths, potential evaluation, body Green operator, or solver.

This is an additive Task 4 diagnostic. It must not change the production default
kernel, remarch a wake with `core_size=0`, or replace the established Task 4 lift
results. Interpret the paired result as follows:

- residuals substantially smaller with `core_size=0`: velocity regularization /
  velocity-potential inconsistency is an important contributor;
- residuals nearly unchanged: centroid sampling and constant-panel
  trace/normal-flux representation are the leading explanation at these grids;
- mixed changes: report dependence on wake proximity, angle, and resolution;
  do not reduce it to a single cause.

Do not infer causality from convergence slopes alone. The decisive quantities
are paired residual ratios at identical mesh and frozen wake state.

## Mandatory repository context

Before implementing, read repository-root `CLAUDE.md`, then
`agent_policies/WORKFLOW.md` and `agent_policies/TESTING.md`. The worktree is
dirty and much of `debug/dirichlet_solve/` and the formulation implementation is
currently untracked; preserve all unrelated work and do not clean or replace
files wholesale.

Relevant files:

- `debug/dirichlet_solve/dirichlet_solve.jl`: all existing Task 4 and Green-
  identity diagnostic code; make the diagnostic implementation here.
- `debug/dirichlet_solve/task4.md`: established Task 4 formulation/results;
  append the completed follow-up and measurements here.
- `debug/dirichlet_solve/green_identity_analysis.md`: current residual and
  conditioning interpretation; revise the near-field explanation after results.
- `debug/dirichlet_solve/dirichlet_solve.md`: parent checklist/index; add only a
  terse follow-up takeaway and do not change any `User approved` checkbox.
- `src/FLOWPanel_elements_fmm.jl`: inspect only unless a diagnostic-only change
  proves impossible. Production kernel behavior must not be changed by this task.

No production source edit should be necessary.

## Existing mathematical and code behavior

The diagnostic evaluates

```math
r_h = S_h\sigma_h - (I-B_h)q_{\mathrm{direct},h},
\qquad
\sigma_h = -u_{f,h}\cdot n.
```

Current direct doublet-panel behavior in `src/FLOWPanel_elements_fmm.jl` is the
reason for this controlled comparison:

- `compute_source_dipole(..., ConstantDoublet, ...)` computes potential as
  `-strength[1] * tan_term`; `reg_term` does not enter the potential.
- The doublet velocity uses
  `val4 = (ri+rip1)/(ri*rip1*rho + reg_term)`.
- `_induced` computes `reg_term = regularize(minimum_distance(...), core_radius)`.
- `regularize(distance, core_size)` returns zero when `distance >= core_size` and
  `(distance-core_size)^2` otherwise. Passing `core_size=0` therefore selects the
  singular/unregularized analytic velocity formula everywhere without changing
  the potential formula.
- `PanelWake` passes its immutable `core_size` to `_induced`; its default is
  `1e-3`, which is also the value explicitly used by the frozen-state loader.

The existing residual implementation is in these symbols in
`debug/dirichlet_solve/dirichlet_solve.jl`:

- `_green_id_wake_opts(state)`: returns `core_size=0.001` plus final-filament
  options for saved flat/rolled wakes.
- `_green_id_load_wake(...)`: constructs a `PanelWake`, updates its trailing
  edge, then loads saved VTK nodes and strengths.
- `_green_identity_metrics(...)`: removes the area mean from `q_direct`, builds
  `Ssigma` with `pnl._source_potential!`, applies the stored bordered LU forward
  to obtain `(I-B)q_direct`, computes the normalized residual and reconstructed-
  versus-direct trace gap, and verifies the bordered solve.
- `_green_identity_panel_case(...)`: evaluates wake scalar potential with
  `_task3_accumulated_potential`, then wake velocity with `pnl.influence!` and
  forms `sigma` through `pnl._split_sigma!`. It snapshots and audits frozen body
  and wake state.
- `_green_identity_semiinf_core(...)`: obtains attached-sheet potential and
  velocity as the attached-wake-on minus attached-wake-off body influence.
- `run_green_identity_residual(...)`: production-resolution frozen-state test.
- `_GREEN_REFINEMENT_GRIDS`, `_refinement_finite_body`,
  `_refinement_semiinf_solve`, `_refinement_row`, and
  `run_green_identity_refinement(...)`: deterministic refinement ladder.

The existing refinement states are:

- flat free wake with attached transition `Da=0.5c`;
- flat free wake with attached transition `Da=0.05c`;
- semi-infinite attached sheet, which touches the trailing edge and is therefore
  an interpretive rather than clean disjoint-support Green-identity check.

The grids are `(n_airfoil,n_span,n_endcap)` = `(81,7,5)`, `(121,10,7)`,
`(161,13,9)`, and `(201,16,11)`. The last grid (~10,360 panels) has an established
mesh/operator defect: its constant-mode defect is about `6.6e-3`, so retain its
raw paired measurements but exclude it from slope fits and causal conclusions.

The existing production-resolution frozen routes are all combinations of:

- wake evolution/source: shared Task 2 single-shot seed, Task 4
  GreenReconstruction-iterated, and Task 3 oracle-iterated;
- state: `flat`, `flat_das005`, and rolled `march`;
- angle: 3.94 degrees under `data/dirichlet_solve/` and 45 degrees under
  `data/dirichlet_solve/alpha45/`.

Use `DirectBackend` throughout this comparison. FMM approximation is not in
scope, and the body-only `B`, area gauge, and bordered LU do not change between
the two velocity modes.

## Implementation design

### 1. Parameterize diagnostic wake construction

Change only debug helpers, preserving current defaults:

```julia
_green_id_wake_opts(state::Symbol; core_size=0.001)
_green_id_load_wake(body, dir, stem, idx, state::Symbol;
                    nwakerows=200, core_size=0.001)
_flat_wake(body, gamma, length_c, c;
           free_row_length_c=0.5, core_size=1e-3)
```

Forward `core_size` into `pnl.PanelWake`. Existing callers without the keyword
must retain their current behavior exactly.

Do not mutate `wake.core_size`: `PanelWake` is immutable. Construct two wakes
from the same saved or deterministic state, one with `0.001` and one with `0.0`.
Loading/rebuilding must leave their nodes, active row count, strengths, and all
non-core wake options identical.

### 2. Return the sampled flux needed for paired diagnostics

Extend `_green_identity_panel_case` and `_green_identity_semiinf_core` so their
returned named tuple also contains a copy of the `sigma` vector used by
`_green_identity_metrics`. Preserve every existing scalar field and existing
call-site behavior.

Also make the raw, area-mean-removed direct trace available to the paired helper
or compute it once outside both velocity evaluations. The preferred structure is
to split panel evaluation into:

1. one `q_direct_raw = _task3_accumulated_potential(body, (wake,))` evaluation;
2. a small helper that samples velocity/forms `sigma` for the selected wake;
3. `_green_identity_metrics` using the shared `q_direct_raw`.

This makes it impossible for the paired potential comparison to differ merely
because it was evaluated twice. If keeping the current function structure is
less invasive, evaluate both and enforce the equality gate below.

### 3. Frozen production-resolution paired probe

Add `run_green_identity_regularization(; path=DIRICHLET_DATA_PATH,
alpha_deg=DEFAULT_ALPHA_DEG)`. Reuse the same template body, `green` bordered
factorization, areas, route table, state table, and saved VTK indices as
`run_green_identity_residual`.

For each route/state, construct and evaluate:

- regularized wake: `core_size=0.001`;
- unregularized wake: `core_size=0.0`.

Use a **wide paired row**, not two unrelated rows. Write
`green_identity_regularization.csv` in the requested output directory with at
least these columns:

```text
angle,wake_state,route,case,n_panels,
regularized_core_size,unregularized_core_size,
regularized_Ssigma_l2,unregularized_Ssigma_l2,
regularized_residual_rel,unregularized_residual_rel,
residual_ratio_unreg_over_reg,residual_percent_change,
regularized_q_recon_vs_direct_rel,unregularized_q_recon_vs_direct_rel,
trace_gap_ratio_unreg_over_reg,
sigma_relative_change,q_direct_relative_change,
regularized_compat_defect_rel,unregularized_compat_defect_rel,
regularized_recon_solve_rel,unregularized_recon_solve_rel
```

Definitions:

```math
\text{residual ratio}=r_{0}/r_{0.001},
\qquad
\text{percent change}=100(r_0/r_{0.001}-1),
```

so negative percent change means improvement. Define `sigma_relative_change`
as `norm(sigma_unreg-sigma_reg)/max(norm(sigma_reg),eps())` and use the analogous
relative L2 definition for `q_direct_relative_change`.

Print a compact table containing case, route, regularized residual,
unregularized residual, ratio, and sigma change.

### 4. Refinement paired probe

Add `run_green_identity_regularization_refinement(; path=DIRICHLET_DATA_PATH,
alpha_degs=(DEFAULT_ALPHA_DEG,45.0),max_panels=Inf)` rather than changing the
schema of `green_identity_refinement.csv`.

At each grid and angle:

- solve the semi-infinite body once and reuse its `gamma` for both deterministic
  flat wake copies;
- build paired flat wakes with core sizes `0.001` and `0.0` on otherwise
  identical finite bodies;
- reuse the same `green`, areas, `sigma_min`, and `q_direct` per paired state;
- append one wide paired row per state.

Write `green_identity_regularization_refinement.csv`. Include grid identifiers,
`n_panels`, `h_proxy`, angle, state, wake length, `mesh_valid_for_fit`, all paired
metrics listed above, and the unchanged constant-mode proxy `sigma_min`.

For the semi-infinite attached-sheet comparison, there is no `PanelWake.core_size`.
The attached panel uses the body's active `kerneloffset`. Implement a keyword on
`_green_identity_semiinf_core`, for example `velocity_core_size=nothing`:

- compute `q_direct_raw` with the original active kernel setting;
- for regularized velocity, preserve the original active setting;
- for unregularized velocity only, temporarily set `body.kerneloffset=0.0` around
  both attached-wake-on and attached-wake-off velocity calls;
- restore `body.kerneloffset`, `suppress_attached_wake`, potential, velocity, and
  strengths in `finally` blocks.

Record the original semi-infinite active offset as the regularized core size.
This offset is normally much smaller than the free-wake `0.001`, so a negligible
semi-infinite change is expected and should not be treated as a failed test.

Compute log-log slopes versus `h=1/sqrt(N)` separately for regularized and
unregularized residuals and trace gaps using only the first three valid grids.
Also print the ratio of paired residuals at each resolution; a changing ratio is
more informative about core effects than either slope alone.

### 5. CLI selectors and backward compatibility

Add selectors to `main`:

```text
green-regularization / green_regularization
green-regularization-refinement / green_regularization_refinement
```

The first honors `--alpha-deg` and `--output-dir`, matching the current
production-resolution diagnostic convention. The second evaluates both angles
and honors the function-level `max_panels` when called directly for smoke tests.

Do not alter the output schema or numerical path of the existing
`green-identity`, `green-refinement`, or `green-conditioning` tasks.

## Required gates

Apply these gates before interpreting results:

1. **Only the velocity core changes.** Paired bodies/wakes must have identical
   nodes, cells, `Das`, wake nodes, wake strengths, active row counts, and wake
   options other than `core_size`.
2. **Potential invariance.** For panel-wake pairs,
   `q_direct_relative_change <= 1e-13`. The analytic potential does not use
   `reg_term`; failure means the experiment changed more than intended.
3. **Regularized baseline reproduction.** If the corresponding established
   `green_identity_residual.csv` or `green_identity_refinement.csv` exists,
   compare matched regularized rows with `rtol=1e-10, atol=1e-13`. If an artifact
   is absent, warn and continue; do not fabricate a baseline.
4. **Solve integrity.** Preserve the existing bordered residual gate
   `recon_solve_rel <= 1e-8`; report actual values, expected near machine epsilon.
5. **Finite fields.** All sampled velocities, sigma vectors, potentials, and
   reported metrics must be finite. An exact singular encounter with
   `core_size=0` is a failed/unsafe diagnostic case and must be reported rather
   than silently regularized or skipped.
6. **Immutable-state audit.** Keep `_state_snapshot` /
   `_assert_frozen_state` around each individual evaluation. Add an explicit
   cross-pair equality audit that ignores only `core_size`.
7. **Defective fine mesh.** Mark the `(201,16,11)` rows invalid for fits; never
   use them to claim improvement or degradation trends.

## Verification sequence

Run from the repository root.

First perform a reduced-grid smoke test without writing production source:

```bash
julia --project -t4 -e '
include("debug/dirichlet_solve/dirichlet_solve.jl")
run_green_identity_regularization_refinement(
    path="data/dirichlet_solve", alpha_degs=(3.94,), max_panels=2000)
'
```

Then run all production-resolution frozen states at both angles:

```bash
julia --project -t4 debug/dirichlet_solve/dirichlet_solve.jl \
  --task green-regularization --alpha-deg 3.94 \
  --output-dir data/dirichlet_solve

julia --project -t4 debug/dirichlet_solve/dirichlet_solve.jl \
  --task green-regularization --alpha-deg 45 \
  --output-dir data/dirichlet_solve/alpha45
```

Run the full refinement comparison:

```bash
julia --project -t4 debug/dirichlet_solve/dirichlet_solve.jl \
  --task green-regularization-refinement \
  --output-dir data/dirichlet_solve
```

Because this task should modify only the debug driver and Markdown results, the
maintained production unit suites are not mandatory. If implementation requires
any edit under `src/`, stop and reassess the design; if such an edit is truly
necessary, run at minimum:

```bash
julia --project -e 'include("test/runtests_unit_fmm.jl")'
julia --project test/formulation_test.jl
```

## Documentation and final robustness review

After successful runs:

1. Add a clearly labeled “Unregularized velocity-kernel follow-up” section to
   `task4.md` with commands, artifact paths, paired tables, gates, and conclusion.
2. Revise `green_identity_analysis.md` lines currently describing “near-field
   wake-to-body quadrature.” The code uses analytic constant-panel influence
   formulas; distinguish:
   - centroid-sampled constant normal-flux representation error;
   - near-singular numerical sensitivity;
   - velocity/potential inconsistency introduced by the velocity-only core.
3. Update the parent `dirichlet_solve.md` takeaway tersely. Do not alter approval
   state or start Task 5.
4. Audit methods, data, and conclusion for consistency as required by
   `CLAUDE.md`: verify paired inputs, regularized reproduction, potential
   invariance, exclusion of the defective grid, and whether the numerical result
   actually supports the stated cause. Revise the conclusion if any gate fails.

## Explicit scope assumptions

- “Same tests” means the existing Green-identity production-resolution frozen-
  state and mesh-refinement residual tests. It does not mean remarching Task 4
  wakes or recomputing lift with an unregularized production kernel.
- “Remove regularization” means exactly `core_size=0.0` for the diagnostic wake
  velocity evaluation. Do not introduce a core-size sweep unless later requested.
- Scalar-potential evaluation remains the existing analytic expression. This is
  specifically a test of whether making velocity consistent with that singular
  potential improves `Ssigma-(I-B)q_direct`.
- Conditioning analysis need not be rerun: `I-B`, its constant null mode, the
  area gauge, and its singular spectrum are unchanged by wake velocity core size.
