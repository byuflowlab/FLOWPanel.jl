# FLOWPanel.jl Documentation — To-Do List

## New Pages to Add

### Quick Start page

**Goal:** A single copy-paste example that gets a new user to a working simulation with minimal friction.

- [ ] Create `docs/src/quickstart.md` (hand-written, not auto-generated)
- [ ] Add to `make.jl` nav as the second item after "Intro"
- [ ] Contents:
  - Install: `] add FLOWPanel`
  - Minimal working example: a swept wing at one angle of attack
    - Define geometry with `gt.surface_loft` (or point to the existing swept wing example)
    - Create a `RigidWakeBody{VortexRing}` body
    - Call `solve!` with `BackslashNeumann()`
    - Compute and print CL
  - Expected output (CL value)
  - "Next steps" links to the full tutorials

---

### Solvers page

**Goal:** A reference page that documents all available solver types, their parameters, and when to use each one.

- [ ] Create `docs/src/solvers.md` (hand-written)
- [ ] Add to `make.jl` nav (after Quick Start or as a top-level tab)
- [ ] Document each solver with: description, constructor syntax, key parameters, when to use, example snippet

#### Solvers to document:

**`BackslashNeumann()`**
- Direct matrix solve using Julia's `\` operator
- Neumann (no-flow-through) boundary condition
- Default solver; good for small to medium problems (~10,000 panels)
- Syntax: `solve!(body, Uinf, BackslashNeumann())`

**`BackslashNeumann(; leastsquares=true)`** (or `BackslashNeumann{true}`)
- Same as above but assembles an overdetermined system and solves via least squares
- Required for watertight (closed) meshes where the system is rank-deficient
- Syntax: note the `leastsquares` flag

**`BackslashDirichlet()`**
- Direct matrix solve with interior Dirichlet boundary condition formulation
- Alternative to Neumann; useful for certain body types
- Document when Dirichlet vs Neumann is preferred

**`KrylovSolver`**
- Matrix-free GMRES via `Krylov.jl`
- Does not assemble the full influence matrix; applies it as a linear operator
- Designed to work with `FastMultipoleBackend` for large problems
- Key parameters: `solver` (Krylov.jl solver type), `restart`, `atol`, `rtol`, `verbose`
- Syntax example

**`FGSSolver`**
- Block fixed-point / Gauss-Seidel solver for coupled multi-body problems
- Iteratively solves each body while treating others as fixed, then iterates to convergence
- Key parameters: `max_iterations`, `tolerance`, `rlx` (relaxation factor)
- Used with `solve!(bodies::Tuple, FGSSolver(...))`
- Document convergence behavior and when relaxation is needed

#### Also document N-body backends:

**`DirectBackend()`**
- O(N²) direct panel-to-panel summation
- Default; exact to machine precision

**`FastMultipoleBackend(; expansion_order, multipole_acceptance, leaf_size)`**
- O(N log N) fast multipole method via `FastMultipole.jl`
- Required for large problems (>~5,000 panels)
- Document key tuning parameters and accuracy tradeoff
- Syntax: pass as `backend` kwarg to solve functions

---

### Geometry Input Options page

**Goal:** A unified reference for all supported mesh input methods, since the current docs scatter this across multiple tutorial pages.

- [ ] Create `docs/src/geometry/meshinput.md`
- [ ] Add to `make.jl` under "Geometry Engine" (either in Grid Generation or as a new top-level subsection)
- [ ] Contents:

**Method 1: Structured parametric mesh via GeometricTools**
- `gt.surface_loft`, `gt.surface_revolution`, etc.
- Then `gt.GridTriangleSurface(grid, dimsplit)`
- Currently documented in Geometry Engine pages — link there

**Method 2: `.msh` file (Gmsh)**
- Install Gmsh, export an ASCII Gmsh 4.1 `.msh` file
- Import via `FLOWPanel.read_gmsh("mesh.msh")` → `Meshes.jl` object
- Convert to FLOWPanel body
- Minimal code snippet

**Method 3: `.stl` file**
- Import via `GeoIO.load("mesh.stl")`
- Convert to FLOWPanel body
- Note: `.stl` has no cell connectivity info; Gmsh is preferred for lifting bodies

**Method 4: Matrix-based input**
- Directly supply `nodes` (3×N matrix of coordinates) and `cells` (connectivity matrix)
- Show constructor: `NonLiftingBody(elements, nodes, cells)` or equivalent
- Useful for programmatic mesh generation without GeometricTools

**Method 5: OpenVSP**
- Export from OpenVSP as DegenGeom
- Process with VSPGeom utilities
- Link to Cessna example

---

## Updates to Existing Pages

### `src/index.md` (auto-generated from README)

- [ ] Remove or update the "No fast multipole method" bullet under Current Limitations — FMM is now available via `FastMultipoleBackend`
- [ ] Update limitations text to reflect current state (FMM exists but may not yet be fully merged to master)
- [ ] Update tested Julia version if needed (currently says v1.10)

---

### Solver Benchmark Tutorial (`examples/sweptwing-solver.md`)

- [ ] The current page documents `solve_backslash!`, `solve_ludiv!`, `solve_gmres!` (old function-based API)
- [ ] Update `examples/sweptwing_solverbenchmark.jl` to use new solver-object API:
  - `BackslashNeumann()` instead of `solve_backslash!`
  - `KrylovSolver(...)` instead of `solve_gmres!`
- [ ] Add `FGSSolver` example (multi-body benchmark)
- [ ] Regenerate the page after updating the example

---

### Duct Least Squares Tutorial (`examples/duct-leastsquares.md`)

- [ ] Currently documents `solve_leastsquares!` (old API)
- [ ] Update `examples/duct_leastsquares.jl` to use `BackslashNeumann(leastsquares=true)` (or equivalent new syntax)
- [ ] Regenerate the page

---

### Blended Wing Body — Unstructured Meshing (`examples/blendedwingbody-gmsh.md`)

- [x] Verify the Gmsh 4.1 import workflow with `FLOWPanel.read_gmsh` / `Meshes.jl`
- [ ] Update any deprecated API calls

---

## Maintenance Tasks

### Re-run all examples and verify correctness

All 23 tutorial examples need to be re-run with the current codebase to ensure:
- They run without errors
- CL/CD/Cm output values match the tables in the docs
- Plots match the images in `src/assets/images/`

Priority order (most likely to be broken first):
- [ ] `examples/sweptwing.jl` — core validation case
- [ ] `examples/sweptwing_aoasweep.jl`
- [ ] `examples/sweptwing_solverbenchmark.jl` — uses old solver API
- [ ] `examples/centerbody.jl`
- [ ] `examples/duct.jl`
- [ ] `examples/duct_leastsquares.jl` — uses old solver API
- [ ] `examples/blendedwingbody_aero.jl` — uses Meshes.jl / `FLOWPanel.read_gmsh`
- [ ] `examples/blendedwingbody_gpucpu.jl` — GPU path
- [ ] `examples/cessna_aero.jl` — complex multi-body
- [ ] `examples/ll_weber.jl`
- [ ] `examples/ll_weber_aoasweep.jl`
- [ ] `examples/ll_stabilityderivatives.jl`
- [ ] `examples/ll_a50k27.jl`

### Update deprecated API calls throughout examples

- [ ] Audit all `examples/*.jl` files for calls to old solver functions (`solve_backslash!`, `solve_ludiv!`, `solve_gmres!`, `solve_leastsquares!`)
- [ ] Replace with new solver-object API (`BackslashNeumann()`, `KrylovSolver(...)`, etc.)
- [ ] Check for any other API changes (body constructors, field names, post-processing functions)

### Regenerate all example docs

After updating `examples/*.jl` files:
- [ ] Run `julia --project docs/src/generate_examples.jl` to regenerate all `src/examples/*.md` files
- [ ] Run `julia --project docs/make.jl` to build and verify the full site
- [ ] Check for Documenter.jl warnings/errors

### Populate API Reference

- [ ] Uncomment the API Reference section in `make.jl` (lines 84–87)
- [ ] Create `src/api.md`, `src/api-elements.md`, `src/api-abstractbody.md`
- [ ] Use `@docs` blocks to auto-document public API functions (solvers, body constructors, post-processing)
- [ ] Include the new solver types and backend types

### Verify images are current

- [ ] Check that all images referenced in docs exist in `src/assets/images/`
- [ ] Re-generate any plots that are outdated or no longer match example output

### Update GPU/CPU example

- [ ] Verify `examples/blendedwingbody_gpucpu.jl` still runs (CUDA dependency, optional import)
- [ ] Update timing numbers if hardware has changed
