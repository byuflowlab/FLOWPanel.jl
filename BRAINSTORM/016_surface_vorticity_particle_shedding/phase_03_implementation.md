# Phase 3 — Implementation

**Status:** IN PROGRESS (started 2026-08-01)

**Prerequisite:** approved Phase 2 architecture — satisfied 2026-08-01

**Phase approval:** [ ]

Implement the architecture **as amended by Phase 2 §13**. In particular the
ghost is the pre-shift final active row read in place (Option B): `PanelWake`
storage, `capacity` derivation, and every source/probe/filament/VTK view are
unchanged, and only the new row 1 is staged, only when `nwakerows == 1`.
Conversion metadata belongs in `src/FLOWPanel_replay.jl`, and conversion
decoding must not reuse the `NoShed()` unknown-tag fallback.

## Staged execution

| Stage | Content | Status |
| --- | --- | --- |
| 0 | Gate bookkeeping; legacy golden-reference characterization test | **DONE** |
| 1 | Strategy axis and dispatch; no behavior change | **DONE** |
| 2 | Reconstruction core (tangent basis, rank-aware gradient, bilinear quadrature) | **DONE** |
| 3 | Transaction, deposition, root/tip closure, handoff | next |
| 4 | Persistence: metadata, replay, warm start | pending |
| 4b | Static sheet-vs-particle equivalence evidence (no timestepping) | pending |
| 5 | Diagnostics and close-out | pending |

Stage 4b is user-requested acceptance evidence: with no time marching, convert
a prescribed sheet once and compare the deposited particles against the sheet
they replace — induced velocity on an off-body probe shell versus standoff,
far-field circulation and impulse moments, legacy-versus-smooth near-field
difference on the same sheet, and internal-edge cancellation across two
successive conversions.

## Prior-phase handoff

Phase 1 supplies the approved theory only after its gate is passed. Phase 2
must then provide the exact typed strategy API, state ownership,
perimeter/handoff ledger, transactional conversion protocol, diagnostics,
failure behavior, metadata/replay/warm-start representation, and file/test
map.

Read the hub and the approved
[Phase 2 architecture](phase_02_architecture.md) before beginning. Do not
infer unresolved architecture directly from the theory.

## Scope

Implement only the approved architecture around the existing
`PanelParticleWake` conversion seam. The likely source touchpoints are
`PanelWake` storage/source indexing, `PanelParticleWake`,
`_convert_to_particles!`, the current conversion-before-row-shift call in
`shed_wake!`, metadata/replay reconstruction, and warm-start state.

The existing edge-jump strategy remains the exact constructor default and
must retain current `method_trailing`/`method_unsteady` behavior. The
surface-vorticity strategy adds rank-aware gradient reconstruction,
physical-area quadrature, true-perimeter closure, handoff bookkeeping,
capacity preflight, and diagnostics. Its conversion uses a staged post-shift
view: the old final active row becomes a non-source ghost only after the new
active row needed by the stencil is constructed.

## Deliverables

- Approved concretely typed conversion strategies and constructor option.
- Explicit active-row capacity, always-allocated ghost geometry, and source
  views that ignore the ghost.
- Gradient, bilinear-quadrature, perimeter, and handoff implementation.
- A preflight/commit conversion transaction: a partially deposited outgoing
  row is never observable.
- Versioned metadata plus backward-compatible replay and warm-start behavior.
- Deterministic validation, diagnostics, and failure paths.
- Focused unit tests for the new strategy and byte-for-behavior legacy
  regression where practical.
- A progress log with exact verification commands and results.

## Required focused tests

- Constructor default and explicit legacy selection are equivalent.
- Existing legacy locations, strengths, counts, line policies, metadata,
  replay, and warm starts remain unchanged.
- Constant and affine planar reconstructions, rotated sign covariance,
  resolution floor, warped/nonuniform geometry, and a full-rank one-sided
  `nwakerows == 1` ghost stencil.
- Active source counts for `nwakerows=1,2,3`; the ghost never contributes
  panel influence.
- No gradient dependence on `Das`, body strengths, virtual TE panels, or an
  outside-zero value.
- Internal-edge cancellation, true perimeter retention, and single-owned
  panel/particle handoff.
- Root/tip versus interior panel-induced velocity-gradient diagnostics.
- Insufficient capacity, invalid geometry, and invalid configuration leave
  both the particle field and panel-row state unchanged.
- New metadata round-trips and old manifests select the legacy default.

## Acceptance gate and progress log

Phase 3 is complete only when the scoped diff matches the approved
architecture, focused tests pass, legacy behavior remains unchanged, and the
user explicitly approves the implementation. Phase 4 remains blocked until
that approval.

- **2026-08-01 — Stages 1 and 2 complete.**

  *Stage 1 — strategy axis, no behavior change.* Added
  `AbstractPanelParticleConversion`, `LegacyEdgeJumpConversion` (the exact
  default), and `SurfaceVorticityConversion{TF}` with §2.1 validation, next to
  the `WakeSheddingMethod` hierarchy in `src/FLOWPanel_wake.jl`. Added the
  `DefaultWakeSheddingMethod` sentinel so the constructor can distinguish
  "caller said nothing" from "caller explicitly asked for the legacy default":
  `_resolve_line_policy` turns it into a fresh `OverlapPPS(1.3, 2)` for the
  legacy strategy, keeps an explicitly supplied policy, and raises
  `ArgumentError` if either line policy is supplied alongside the smooth
  strategy (§2.2). The sentinel itself throws if it ever reaches
  `_shed_particles!`. `PanelParticleWake` gained `conversion`,
  `conversion_workspace`, `conversion_diagnostics`, and `conversion_count`
  (three new type parameters); legacy instances store `nothing` for workspace
  and diagnostics. `_convert_to_particles!(wake)` now dispatches on
  `wake.conversion`, and the legacy method is the previous body verbatim.

  *Stage 2 — reconstruction core.* Pure, simulation-free functions per §5.2–§5.4:
  `_deterministic_tangent_basis`, `SurfaceGradientResult`,
  `_reconstruct_surface_gradient` (scaled `2x2` SVD with the rank-2/1/0 policy;
  the minimum-norm branch falls out of the same truncated-SVD sum),
  `_surface_vorticity` (the `kappa = -n x grad(muhat)` sign), the bilinear panel
  family (`_bilinear_position`, `_bilinear_derivatives`, `_bilinear_normal`),
  `_subdivision_counts`, `_subcell_area` (`2x2` Gauss on the Jacobian norm),
  `_validate_wake_panel`, and the typed `WakeGeometryError`.

  Verification (`julia --project -e 'include("test/runtests_unit_wake.jl")'`):
  **197/197 pass** in `Free Wakes` — 33 legacy golden reference (still bit-exact
  after the dispatch refactor, which is the point of Stage 0), 23 strategy axis,
  51 reconstruction core. The reconstruction tests cover constant and affine
  planar fields (streamwise, spanwise, combined), the exact package sign and its
  antisymmetry, rigid-rotation covariance of both the gradient and the derived
  surface vorticity, rank-1 minimum-norm and rank-0 behavior with observable
  flags, rejection of a vanishing metric scale and non-finite stencils, exact
  area partitioning and quadrature convergence on a warped panel, the
  subdivision floor and stretched-panel counts, and rejection of coincident,
  collapsed, bow-tie, and non-finite panels.

  Full-suite gate: `test/runtests.jl` passes lines 7–20 (all unit files) with no
  failures. **Pre-existing, unrelated blocker:**
  `test/runtests_unit_pitching_wing_exp.jl:4` aborts on a missing data asset
  (`data/pitching_wing_exp/load.jl`; the directory is absent from the repo and
  is not a pending deletion), which stops `runtests.jl` at line 21. Lines 22–27
  were therefore run separately and all pass. Note also that `julia` exits 0
  here even on such a load error, so the gate must be judged from log content,
  not exit status.

  *Design notes verified for Stage 3 (Option B mechanics).* Per-step ordering is
  `update_TE!` (`FLOWPanel_simulate.jl:460`, placing wake node row 1 at
  `body TE + Das`) → solve → `propagate!` → `propagate_kinematics!` →
  `shed_wake!` (:1125-1127). At conversion time the outgoing ghost is active row
  `N`: node rows `N`, `N+1` and strength row `N`. For `N >= 2` its streamwise
  upstream neighbour is existing active row `N-1`, already in storage — nothing
  is staged. For `N == 1` the upstream neighbour is the row the shift is about
  to create; `shed_wake!(::PanelWake, ...)` (:1190-1212) does **not** build new
  row-1 nodes, it copies row 1 into row 2 and writes only
  `strength[1,1,:] = _get_wakestrength_mu(...)`. So the staged upstream sample
  is: position = centroid of node row 1 (the current `TE + Das` line, i.e. the
  ghost's own upstream edge), strength = `_get_wakestrength_mu(system, ...)`.
  That is a genuine one-sided two-panel stencil with two distinct strengths at
  two distinct positions, and it introduces no new `Das` dependence beyond the
  geometry that already defines the ghost. **Consequence:** `system` must be
  threaded into the conversion, so Stage 3 should give the smooth strategy the
  entry point `_convert_to_particles!(wake, conversion, system)` and pass
  `system` from `shed_wake!(::PanelParticleWake, ...)` (:1224-1236); the legacy
  method can ignore it.

- **2026-08-01 — Stage 0 complete.** Phase 2 approval recorded across the hub,
  the phase document, and `BRAINSTORM/INDEX.md`, with the §13 amendment
  (ghost-storage Option B; the two §10/§9.1 corrections). Added
  `make_conversion_fixture` to `test/test_helpers.jl` — a deterministic
  `PanelParticleWake` whose panel-wake nodes and strengths are written directly
  so `_convert_to_particles!` can be exercised **with no time stepping** — and
  captured its output as
  `test/data/legacy_wake_conversion_reference.jl` for `wraps ∈ {false, true}`
  × `nwakerows ∈ {1, 2, 3}`. The new testset "Legacy edge-jump conversion
  golden reference (BRAINSTORM 016)" in `test/runtests_unit_wake.jl` compares
  positions, vector strengths, smoothing widths, the vestigial scalar
  circulation, and particle counts by **exact equality**, and asserts that the
  wrapping chain deposits fewer particles than the open chain (no root/tip
  closure). This closes the pre-existing coverage gap in which nothing pinned
  `_convert_to_particles!` directly. No source code has been modified yet.

  Command and result:
  `julia --project -e 'include("test/runtests_unit_wake.jl")'` →
  **123/123 pass** (`Free Wakes`), of which 33 are the new reference testset.

- **2026-07-29 — Phase seeded.** Blocked pending explicit Phase 2 approval;
  no implementation has started.
