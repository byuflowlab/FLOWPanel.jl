# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

FLOWPanel.jl is a 3D panel method solver for low-speed (inviscid, incompressible) aerodynamics, written in Julia. It supports non-lifting bodies (source/doublet panels), lifting bodies (vortex-ring panels with rigid wake), and unsteady simulation. Computation is accelerated via the `FastMultipole` package (FMM) or direct N-body evaluation.

## Commands

```bash
# Run all tests from shell
julia --project -e 'include("test/runtests.jl")'

# Run a specific test file
julia --project -e 'include("test/runtests_unit_fmm.jl")'
julia --project -e 'include("test/runtests_unit_body.jl")'
julia --project -e 'include("test/runtests_unit_solver.jl")'
julia --project -e 'include("test/runtests_unit_liftingbody.jl")'
julia --project -e 'include("test/runtests_analytical.jl")'

# Run an example
julia --project examples/sphere.jl
```

From the Julia REPL (with `--project`):
```julia
# Run all tests via Pkg
] test

# Run a specific test file
include("test/runtests_analytical.jl")

# Run an example
include("examples/sphere.jl")
```

## Architecture

### Module Load Order (`FLOWPanel.jl`)
Files are included in this order (reflects the dependency chain):
1. `elements` — Panel element types and their direct velocity/potential kernels
2. `fmm` — N-body backend abstraction (`AbstractBackend`, `DirectBackend`, `FastMultipoleBackend`)
3. `abstractbody` — `AbstractBody{E,N,TF,DBC}` interface definition
4. `nonliftingbody` — `NonLiftingBody{E,N,TF,DBC}` concrete type
5. `abstractliftingbody` — `AbstractLiftingBody` interface
6. `liftingbody` — `RigidWakeBody{E,N,TF,DBC}` concrete type
7. `solver` — Solver types: `BackslashNeumann`, `BackslashDirichlet`, `KrylovSolver`, `FGSSolver`
8. `elements_fmm` — FMM integration of panel elements (overloads `_Uind!`, `_phi!`)
9. `frames` — `ReferenceFrame` kinematics for unsteady simulation
10. `liftingline` — Lifting-line coupling utilities
11. `utils` — Geometry helpers (`simplewing`, math utilities)
12. `postprocess` — `calcfield_U!`, `calcfield_P!`, `calcfield_F!`, etc.
13. `wake` — `AbstractFreeWake`, `PanelWake`, and `solve!` for body+wake systems
14. `simulate` — `simulate!` time-marching loop
15. `monitor` — (conditionally loaded) PyPlot-based monitors

### Key Types

**Element types** (`FLOWPanel_elements.jl`):
- `ConstantSource`, `ConstantDoublet`, `VortexRing`, `ConstantVortexSheet`, `UniformVortexSheet`
- Bodies are parameterized as `Body{E,N,TF,DBC}` where `E` is a `Union` of element types, `N` is the number of element types, `TF` is the float type, and `DBC` is a bool flag for Dirichlet boundary conditions.

**Body hierarchy**:
- `AbstractBody{E,N,TF,DBC}` — base interface; holds `nodes`, `cells`, `strength`, `fields`, etc.
- `NonLiftingBody{E,N,TF,DBC}` — for closed non-lifting surfaces
- `AbstractLiftingBody{E,N,TF,DBC}` → `RigidWakeBody{E,N,TF,DBC}` — includes `shedding` edges for rigid wake

**Solver hierarchy** (`FLOWPanel_solver.jl`):
- `AbstractSolver`
  - `AbstractMatrixfulSolver{LS}` — assembles full system matrix; `LS` flag enables least-squares formulation
    - `BackslashNeumann` — direct matrix solve (Neumann formulation)
    - `BackslashDirichlet` — direct matrix solve (Dirichlet/interior-Dirichlet formulation)
  - `AbstractMatrixFreeSolver` — applies system operator without storing the full matrix
    - `KrylovSolver` — matrix-free GMRES via `Krylov.jl`, works with `FastMultipoleBackend`
    - `FGSSolver` — matrix-free Gauss-Seidel-style solver for a single body; key params: `max_iterations`, `tolerance`, `rlx` (relaxation)

**Multi-body solving**: use `solve!(bodies::Tuple, solvers::Tuple)`, including one `FGSSolver(body)` per body when using FGS in coupled solves.

**N-body backends** (`FLOWPanel_fmm.jl`):
- `DirectBackend` — O(N²) direct summation via `FastMultipole.direct!`
- `FastMultipoleBackend` — O(N log N) FMM via `FastMultipole.fmm!`; key params: `expansion_order`, `multipole_acceptance`, `leaf_size`

### Key Conventions

- `gt` is the alias for `GeometricTools` (used everywhere for mesh generation/manipulation).
- Panel influence functions follow Hess & Smith / Katz & Plotkin sign conventions; note GeometricTools defines normals with the right-hand rule, which is opposite to K&P — the element kernels account for this with explicit sign flips.
- `strength[:, j]` stores the strength of the j-th element type for all panels (ncells × N matrix).
- Control points are offset from panel centers along the normal by `CPoffset * characteristiclength(panel)`.
- `Uinf` passed to `solve` is a 3×ncells matrix (freestream at each control point).
- The `fields` vector in each body tracks which solution fields have been computed (used by post-processing).

### Unsteady Simulation Flow

`simulate!(body, frames, maneuver!, Vinf, t_range)` creates a `PanelWake` and calls the time-marching `simulate!` loop. Each step:
1. `propagate_kinematics!` — updates body node positions from `ReferenceFrame` tree
2. `solve!` — evaluates wake influence, solves body BCs, evaluates body influence on wake
3. Wake rollup / panel convection

### GeometricTools Integration

Meshes are built using `GeometricTools.GridTriangleSurface`. Grids are created parametrically (e.g., as `gt.Grid`), transformed (e.g., `gt.transform!`, `gt.lintransform!`), then split into triangles (`gt.GridTriangleSurface(grid, dimsplit)`). Unstructured meshes can also be imported via `GeoIO`/`Meshes.jl`.

### Differentiation Support

The solver supports forward-mode AD (ForwardDiff) and reverse-mode AD (ReverseDiff) through `ImplicitAD`. The `solve_ludiv!` function is overloaded in `ImplicitAD` to efficiently differentiate through the linear solve.
