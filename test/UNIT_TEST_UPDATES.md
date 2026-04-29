# Unit Test Update Candidates

This file tracks unit tests that are currently weak or only partially meaningful,
plus concrete upgrades that would make them verify behavior rather than just
construction, non-throwing execution, or diagnostic output.

## Completed / Removed

### `test/runtests_unit_liftingbody.jl`

- `RigidWakeBody / seeded shedding backslashcoupled`
  - Removed from the lifting-body unit tests because the same seeded `RigidWakeBody`
    `BackslashCoupled` solve is already covered in `test/runtests_unit_solver.jl`
    with a post-solve Dirichlet residual assertion.

## High Priority

### `test/runtests_unit_liftingbody.jl`

- `RigidWakeBody / construction`
  - Current issue: mostly checks dimensions and flags.
  - Improve by asserting construction invariants:
    - `shedding_full` maps seeded trailing-edge panels to the expected cell/edge pairs;
    - `Das` has unit direction columns after initialization or a known default;
    - constructed grid bodies preserve cells/nodes from `GridTriangleSurface`;
    - normals are unit length and consistently oriented for fixed/flipped cases.

- `RigidWakeBody / shedding edge consistency`
  - Current issue: checks only that some trailing-edge panels exist and one row is positive.
  - Improve by asserting each populated `shedding_full` column corresponds to a real shared edge in `body.cells`, and that non-TE panels have sentinel/default values.

### `test/runtests_unit_solver.jl`

- `Solvers / Backslash construction`
  - Current issue: only checks constructors return a `Backslash`.
  - Improve by asserting assembled matrices and solver work buffers are correct for each formulation:
    - matrix size and finite entries;
    - LU factorization can solve a known RHS;
    - Neumann body has exterior control points and Dirichlet body has interior control points.

- `Solvers / FGSSolver construction`
  - Current issue: checks only default option fields.
  - Improve by asserting per-body storage matches the bodies:
    - `Uext`, `phi_ext`, and optional history buffers have expected shapes and element types;
    - constructing with multiple bodies preserves one storage block per body;
    - Dirichlet/Neumann body `CPoffset` values are restored after construction.

- `Solvers / numtype`
  - Current issue: only checks two direct type cases.
  - Improve by adding mixed node/strength precision cases and solver construction checks so `numtype` is validated as the precision source for solver buffers.

### `test/runtests_unit_fmm.jl`

- `FMM Backends / backend construction`
  - Current issue: checks only field defaults and one custom kwarg.
  - Improve by asserting custom kwargs are all honored together, and by adding a small direct-vs-FMM evaluation with non-default backend settings to prove the constructed backend is usable.

## Medium Priority

### `test/runtests_unit_body.jl`

- `Abstract Body Geometry And NonLiftingBody / calc_areas!`
  - Current issue: the body-level path checks only output length.
  - Improve by asserting body-level areas match `calc_areas(body.nodes, body.cells)` and known values for `NODES_2TRI`.

- `Abstract Body Geometry And NonLiftingBody / rotate! reset! get_cell grid2cells`
  - Current issue: reset checks only a subset of mutable fields.
  - Improve by filling `strength`, `velocity`, `potential`, `P`, `F`, and any cached geometry/state that `reset!` promises to clear, then assert all are restored.
  - Improve `grid2cells` by checking exact cell connectivity for `make_basic_triangle_surface()` rather than only bounds and shape.

- `Abstract Body Geometry And NonLiftingBody / NonLiftingBody construction`
  - Current issue: mostly dimension and flag checks.
  - Improve by asserting:
    - initialized arrays are zeroed;
    - `neighbor` connectivity matches expected neighboring cells for `NODES_2TRI`;
    - `vtk_cells` encodes the same connectivity as `cells`;
    - `has_dirichlet_bc` changes the control-point offset orientation used during solves.

### `test/runtests_unit_fgs_history.jl`

- `FGSSolver history & projection / construction with history kwargs`
  - Current issue: checks option plumbing and shapes, but not that invalid history/projection configurations are rejected or normalized.
  - Improve by adding validation tests for negative history length, projection enabled with zero history length, and projection order larger than available saved states.

- `FGSSolver history & projection / allocation-free hot path`
  - Current issue: useful but brittle as a unit assertion on exact allocation count.
  - Improve by keeping this as a performance guard but pairing it with semantic checks that projected/saved values remain correct after the allocation measurement.

### `test/runtests_unit_wake.jl`

- `Free Wakes / ParticleMaintenance separates mixed policies`
  - Current issue: checks type partitioning but not execution order or combined behavior.
  - Improve by applying a mixed maintenance policy to a wake containing particles that should be trimmed and merged, then assert final particle count, positions, and circulation.

- `Free Wakes / PanelParticleWake particle maintenance merges particles`
  - Current issue: checks only that count becomes one.
  - Improve by asserting the merged particle conserves circulation and has the expected position/core-size behavior.

### `test/runtests_unit_warmstart.jl`

- `simulate_warmstart! consistency`
  - Current issue: strong equality checks after restart, but no negative/control case.
  - Improve by adding assertions that restart metadata selects the intended restart step and that changing `restart_step` changes the initialized state as expected.

## Lower Priority / Cleanup

- Remove diagnostic `println` calls from unit tests unless they are behind a debug flag. Current examples include solver and lifting-body solve tests.
- Split mesh-backed solver checks that write VTK files from pure unit tests, or write outputs to a temporary directory so running tests does not mutate tracked `test/check_mesh_body_*` artifacts.
- Rename duplicate or vague testset names, especially repeated `Backslash Dirichlet construction and solve`, `backslash`, and `backslashcoupled`, so failures identify the behavior being tested.
- Review commented-out blocks in `test/runtests_unit_solver.jl`; remove stale comments or convert still-relevant cases into active tests.

## Suggested Iteration Order

1. Strengthen `Backslash construction` and `FGSSolver construction` with storage/matrix/control-point invariants.
2. Make mesh-backed solver tests write to a temp directory or disable VTK output in unit runs.
3. Strengthen `NonLiftingBody construction`, `calc_areas!`, and `grid2cells` with exact expected values.
4. Expand wake maintenance tests to verify combined policy behavior and merge conservation.
