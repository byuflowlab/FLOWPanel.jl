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
15. `simulate_monitors` — `AbstractMonitor` types called each timestep by `simulate!`
16. `monitor` — (conditionally loaded) PyPlot-based visualization monitors

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
- **Compute `RigidWakeBody` shedding from the *constructed* cells, not the raw mesh.** With `ensure_winding=true` (default) the constructor calls `ensure_consistent_winding!`, which flips/reorders node indices within `cells` in place (and flips whole components when `watertight=true`). Since `shedding` references panels by cell-local node ordering, shedding computed from raw `cells` (e.g. from `meshes2nodes_cells`) becomes inconsistent with the re-wound cells the solver uses — this does **not** error, but the wake attaches at the wrong edges and the body sheds almost no circulation (observed: rotor-hover CT silently collapsed ~3.6×, 0.0505→0.014). Correct procedure: build a `noshedding` body first, run `calc_shedding_from_seed` on *its* `.nodes`/`.cells`, then rebuild with the shedding (see `examples/rotor_hover_convergence.jl`).

**Monitor hierarchy** (`FLOWPanel_simulate_monitors.jl`):

Monitors are callables `(systems, wakes, frames, uinf, i_step, dt) -> nothing` invoked each timestep by `simulate!`. They declare data contracts via two traits:
- `monitor_provides(m)` — tuple of symbols (`:P`, `:F`) written by this monitor
- `monitor_requires(m)` — tuple of symbols that must be written by an *earlier* monitor

`audit_monitors(monitors)` validates the ordering at the start of `simulate!` and throws `ArgumentError` on the first unmet dependency.

Concrete monitors:
- `PressureBernoulli(rho; unsteady, correct_kuttacondition, clip)` — populates `body.P` via steady or unsteady Bernoulli; provides `:P`
- `PressureLaplace(bodies, rho; atol, rtol, itmax, preconditioner, reference_panel, reference_pressure, rebuild_every_step, verbose)` — populates `body.P` by solving a sparse panel-centered surface pressure Poisson equation (CG from Krylov.jl); provides `:P`. Must be constructed with the actual body objects for preallocation. Receives `dt` from `simulate!` at runtime; do **not** pass `dt` at construction.
- `ForceMonitor(nt, i_system; i_frame, normalization, correct_kuttacondition, verbose)` — populates `body.F`, integrates force/moment, stores histories in `.force` and `.moment` (3×nt); requires `:P`, provides `:F`
- `KuttaJoukowskiForce(body, nt, i_system; rho, backend, normalization, verbose)` — independent Kutta–Joukowski cross-check; evaluates `ρ Σ γ (Δs × V)` at edge midpoints via a `FastMultipole.ProbeSystem`; stores history in `.force` (3×nt)

Normalization callables `(CF, CM, systems, frames, uinf) -> (CF_norm, CM_norm)`:
- `WingNormalization(rho, Sref, Lref)` — divides by `0.5 ρ |U∞|² Sref` (and `… Lref` for moments)
- `NoNormalization()` — pass-through, returns dimensional values
- `RotorNormalization(rho, D, i_frame)` — divides by `ρ n² D⁴` / `ρ n² D⁵`, reads `n` (rev/s) from `frames[i_frame].ω / 2π`
- `RotorNormalization2(rho, D, i_frame)` — divides by disk-area dynamic pressure; reads `ω` (rad/s) from `frames[i_frame].ω`

`PressureLaplace` internals:
- Sparse FV surface Laplacian `L` assembled from shared-edge weights `w_ij = ℓ_ij / d_ij`; gauge-fixed by pinning one reference panel to `reference_pressure`
- RHS built from edge-integrated tangential material acceleration: `b_i -= ρ ℓ_ij (a_t,j - a_t,i)·ê_ij`
- Unsteady term `∂u/∂t` approximated by finite difference of successive monitor calls (stored in `velocity_dot` as negative previous velocity)
- Velocity gradient `∇u` needed for convective term obtained analytically: `∇u_induced` is the FastMultipole Hessian populated into `body.velocity_gradient` during the per-step `influence!` calls (gated by `monitor_requires_body_hessian(::PressureLaplace) = true` flipping `body.needs_velocity_gradient[]` in `simulate!`); the kinematic part `[Ω]_×` is reconstructed from `body.angular_velocity` accumulated in `kinematic_velocity!`
- Body count checked at call time against `length(m.b)`; no identity (`objectid`) check — caller is trusted to provide compatible bodies
- `L`, preconditioner, and CG workspace are reused by default; set `rebuild_every_step=true` for deforming geometry that needs a fresh pressure Laplacian each monitor call
- Preconditioners: `JacobiPressurePreconditioner` (default, O(N)), `NoPressurePreconditioner`; `IncompleteCholeskyPressurePreconditioner` and `AMGPressurePreconditioner` reserved but not implemented

### Unsteady Simulation Flow

`simulate!(body, frames, maneuver!, Vinf, t_range)` creates a `PanelWake` and calls the time-marching `simulate!` loop. Each step:
1. `propagate_kinematics!` — updates body node positions from `ReferenceFrame` tree
2. `solve!` — evaluates wake influence, solves body BCs, evaluates body influence on wake
3. Wake rollup / panel convection
4. Each monitor in the `monitors` tuple is called in order

### GeometricTools Integration

Meshes are built using `GeometricTools.GridTriangleSurface`. Grids are created parametrically (e.g., as `gt.Grid`), transformed (e.g., `gt.transform!`, `gt.lintransform!`), then split into triangles (`gt.GridTriangleSurface(grid, dimsplit)`). Unstructured meshes can also be imported via `GeoIO`/`Meshes.jl`.

### Differentiation Support

The solver supports forward-mode AD (ForwardDiff) and reverse-mode AD (ReverseDiff) through `ImplicitAD`. The `solve_ludiv!` function is overloaded in `ImplicitAD` to efficiently differentiate through the linear solve.

## Running Simulations Iteratively

When you launch a simulation from an example (especially diagnostic / repro runs of an in-progress investigation), default to **leaving VTK output on** — the I/O cost is small relative to the value of having ParaView-ready state to inspect after the fact. Do not set `SAVE_VTK=false` unless the user asks for a no-output run.

To avoid filling the disk across repeated iterations, **write each new run over the previous run's directory** rather than creating a new one per attempt. Two acceptable patterns:

- Reuse the example's default `save_path` (typically `data/<run_name>/`) and let `simulate!` overwrite per-step files in place. Before launching, `rm -rf` the previous run's directory so stale steps past the current run's length don't linger.
- If side-by-side comparison is needed, pick one persistent sibling dir per scenario (e.g. `data/<run_name>_nocouple/`) and overwrite it on each rerun of that scenario — don't suffix with timestamps or attempt numbers.

When the user explicitly wants to preserve a previous run for comparison, ask before overwriting and offer to `mv` the old directory aside.

## Response Style

When you respond to prompts, if you ever need to ask me questions where I decide between possible options (not including asking permissions to perform shell commands), and you suspect any of the option will be token-heavy, also include a brief estimate of the token count required as a fraction of my 5-hour limit.
