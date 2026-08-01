# Phase 4 — Testing, Verification, and Validation

**Status:** BLOCKED ON EXPLICIT PHASE 3 APPROVAL

**Prerequisite:** approved Phase 3 implementation

**Phase approval:** [ ]

## Prior-phase handoff

Phase 1 defines the mathematical identities and acceptance scenarios. Phase 2
must define the approved API, diagnostics, and state contracts. Phase 3 must
provide the implementation and focused mechanical-test record. Read all
approved prior phases before beginning.

## 1. Mathematical verification

Use deterministic synthetic sheets independent of aerodynamic validation:

- Constant \(\mu\): zero distributed particles, with circulation represented
  only by the retained true perimeter closure.
- Affine \(\mu\) on planar and rigidly rotated sheets: exact reconstructed
  gradient, correct package sign, and exact quadrature-integrated particle
  strength.
- Nonuniform planar and warped quad grids: measured convergence under panel
  refinement and particle-spacing refinement.
- Telescoping/discrete-Stokes ledger: internal panel boundaries cancel;
  physical perimeter and panel/particle handoff circulation are each counted
  exactly once.
- Resolution floor: for a nonzero full-rank field, one particle lies at the
  bilinear panel centroid when both subdivision counts are one.
- `nwakerows == 1`: exact one-sided streamwise reconstruction from the newly
  shed active row and outgoing non-source ghost, with no artificial
  outside-zero or `Das` sample.
- Cancelled-edge and retained-filament handoff alternatives: signed
  circulation/flux balance for one row, multiple spanwise panels, and
  repeated conversions.

Record tolerances, observed orders, signed vector balances, geometry scales,
and rank thresholds.

## 2. Mechanical regression

Prove that the new strategy axis is opt-in:

- compare omitted strategy, explicit legacy strategy, and the pre-change
  reference;
- require unchanged particle locations, strengths, counts,
  `method_trailing`/`method_unsteady` behavior, metadata, replay, and warm
  starts;
- round-trip smooth-strategy metadata and continuation state;
- verify active source counts remain \(N\) despite \(N+1\) stored panel
  strengths and \(N+2\) node rows;
- verify the ghost never contributes panel influence;
- exercise validation and failure paths; and
- force insufficient capacity and verify that neither particles nor panel
  rows are partially committed.

Run the narrow wake, simulate, replay, and warm-start suites selected by the
approved Phase 3 test map before broader regression.

## 3. Aerodynamic validation

Compare legacy edge-jump and smooth surface-vorticity strategies with
identical geometry, wake, particle, solver, timestep, maintenance, viscous,
SFS, and monitor settings. Only the conversion strategy may differ.

### Case order

1. **Suddenly started wing.** Establish signs, startup circulation,
   lift history, induced velocity, and convergence in a compact transient.
2. **Rotor hover.** Evaluate the production-relevant helical wake only after
   the wing case and all mathematical/mechanical gates pass.

### Quantities and refinement axes

For both cases compare:

- induced velocity at fixed probes and relevant surfaces;
- integrated circulation and impulse;
- particle count and spatial distribution;
- numerical stability and conservation drift; and
- \(C_L\) for the wing or \(C_T\) for the rotor.

Refine panel-row count, spanwise panels, timestep, \(\sigma\), overlap, and
root/tip line-particle spacing. Report cost and particle count alongside
accuracy so apparent smoothing caused only by coarser effective resolution
is visible. Separately report panel-induced velocity-gradient norms at new
root/tip and interior particle locations. Non-finite values, instability, or
nonconvergent growth trigger reconsideration of a conservative diffuse
root/tip closure. Negative results, lack of convergence, or sensitivity
merely transferred to another parameter are valid conclusions.

## Acceptance gate and progress log

Phase 4 and Item 016 are complete only after mathematical verification,
mechanical regression, and both ordered aerodynamic validations satisfy the
approved criteria and the user explicitly approves the conclusions. Only then
may the item-level columns in `BRAINSTORM/INDEX.md` be checked.

- **2026-07-29 — Phase seeded.** Blocked pending explicit Phase 3 approval;
  no tests or simulations have started.
