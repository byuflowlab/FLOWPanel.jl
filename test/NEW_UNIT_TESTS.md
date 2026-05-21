# FLOWPanel.jl — Proposed Unit Test Suite

This document specifies a comprehensive set of unit tests for FLOWPanel.jl, organized by source module. Each test includes setup, assertions, and the property being verified. Priority levels: **P0** (critical), **P1** (important), **P2** (nice-to-have).

Lifting-line functionality is excluded.

---

## 1. Common Test Fixtures

These minimal meshes are reused across many test sections.

```julia
using Test
using LinearAlgebra: norm, dot, I, cross
import FLOWPanel as pnl
import GeometricTools as gt

# --- Single triangle in the XY plane (right triangle, legs=1) ---
const NODES_1TRI = Float64[0 1 0;    # x
                            0 0 1;    # y
                            0 0 0]    # z   (3×3)
const CELLS_1TRI = reshape(Int[1, 2, 3], 3, 1)  # 3×1

# --- Two triangles forming a unit square in the XZ plane ---
const NODES_2TRI = Float64[0 1 1 0;  # x
                            0 0 0 0;  # y
                            0 0 1 1]  # z   (3×4)
const CELLS_2TRI = Int[1 1;
                        2 3;
                        3 4]          # 3×2

# --- Octahedron (6 vertices, 8 faces) — closed watertight surface ---
const NODES_OCT = Float64[ 1 -1  0  0  0  0;
                            0  0  1 -1  0  0;
                            0  0  0  0  1 -1]  # 3×6
const CELLS_OCT = Int[1 1 2 2 1 1 2 2;
                       3 5 3 5 4 6 4 6;
                       5 3 6 6 3 5 3 3]  # 3×8  (verify consistent winding)
```

---

## 2. FMM Backends  (P0)

**Source:** `src/FLOWPanel_fmm.jl`, `src/FLOWPanel_elements_fmm.jl`

### 2.1 Backend construction

| Property | Test |
|----------|------|
| DirectBackend | `DirectBackend()` constructs without error |
| FMM defaults | `FastMultipoleBackend()` has `expansion_order==5`, `multipole_acceptance==0.5`, `leaf_size==10` |
| FMM custom | `FastMultipoleBackend(expansion_order=10)` stores `expansion_order==10` |

### 2.2 Direct vs FMM equivalence — NonLiftingBody

**Setup:** ~20-panel NonLiftingBody{ConstantSource} (octahedron or small sphere). Set unit source strengths. Evaluate induced velocity at 5 external target points.

| Property | Test |
|----------|------|
| Velocity agreement | Direct and FMM (expansion_order=10) results agree within `atol=1e-3` |
| Potential agreement | Same test for scalar potential |

### 2.3 Direct vs FMM equivalence — RigidWakeBody

**Setup:** Minimal 4-panel RigidWakeBody{VortexRing} flat plate with shedding edges. Unit vortex ring strengths.

| Property | Test |
|----------|------|
| Velocity agreement | Direct and FMM backends agree within `atol=1e-3` |

---

## 3. Abstract Body Geometry  (P0)

**Source:** `src/FLOWPanel_abstractbody.jl`

### 3.1 `calc_normals!`

| Property | Test |
|----------|------|
| XY-plane triangle | Normal is `[0, 0, ±1]` (sign depends on winding). Assert unit length. |
| Reversed winding | Flip cell vertex order → normal flips sign |
| Multiple panels | Two-triangle mesh: each panel has correct normal |

### 3.2 `calc_areas!`

| Property | Test |
|----------|------|
| Unit right triangle | `area == 0.5` |
| Scaled triangle | Legs 3 and 4 → `area == 6.0` |
| Equilateral side 2 | `area ≈ sqrt(3)` |
| Allocating version | `calc_areas(nodes, cells)` returns vector of correct length |

### 3.3 `calc_controlpoints!`

| Property | Test |
|----------|------|
| Zero offset | Control point == centroid == mean of 3 vertices |
| Positive offset | CP displaced along +normal direction from centroid |
| Offset magnitude | Displacement = `off * characteristiclength * norm(normal_unnormalized) / norm(normal)` — verify numerically |
| `characteristiclength_unitary` | Offset magnitude = `off` exactly |

### 3.4 Characteristic length functions

| Property | Test |
|----------|------|
| `characteristiclength_unitary` | Always returns `1` regardless of triangle shape |
| `characteristiclength_bbox` | For unit right triangle in XY plane: bounding-box diagonal = `sqrt(2)` |
| `characteristiclength_maxdist` | For unit right triangle: max dist from node 1 to others = `sqrt(2)` |
| `characteristiclength_sqrtarea` | Returns `sqrt(sqrt(0.5))` for unit right triangle (area=0.5) |

### 3.5 `calc_tangents!` and `calc_obliques!`

| Property | Test |
|----------|------|
| Orthogonality | `|dot(tangent, normal)| < 1e-12` and `|dot(oblique, normal)| < 1e-12` and `|dot(tangent, oblique)| < 1e-12` |
| Unit length | `norm(tangent) ≈ 1` and `norm(oblique) ≈ 1` |
| Consistent across calls | Calling twice gives same result |

### 3.6 `rotate!`

**Setup:** NonLiftingBody{ConstantSource} from `NODES_2TRI`/`CELLS_2TRI`.

| Property | Test |
|----------|------|
| 90° yaw | After `rotate!(body, 0, 0, 90)`, x-coords become y-coords (with sign flip). |
| `Oaxis` update | `body.Oaxis` matches expected rotation matrix |
| Translation | `rotate!(body, 0, 0, 0; translation=[1,0,0])` shifts all nodes by `[1,0,0]` |
| Composition | Two sequential 90° yaw rotations ≈ single 180° yaw rotation (within `1e-10`) |
| Field reset | After rotate with `reset_fields=true`, velocity/potential/Cp are zeroed |

### 3.7 `reset!`

**Setup:** Body with manually set nonzero velocity, potential, Cp, F.

| Property | Test |
|----------|------|
| All fields zeroed | `all(body.velocity .== 0)`, `all(body.potential .== 0)`, `all(body.Cp .== 0)`, `all(body.F .== 0)` |

### 3.8 `get_cell`

| Property | Test |
|----------|------|
| Correct indices | `get_cell(body, 1)` returns the 3 node indices of cell 1 matching `cells[:, 1]` |

### 3.9 `grid2cells`

**Setup:** Create a small `gt.GridTriangleSurface`.

| Property | Test |
|----------|------|
| Correct shape | Output is `(3, ncells)` |
| Valid indices | All indices in `[1, nnodes]` |

---

## 4. NonLiftingBody Construction  (P0)

**Source:** `src/FLOWPanel_nonliftingbody.jl`

### 4.1 From nodes and cells

| Property | Test |
|----------|------|
| Dimensions | `body.nnodes == 4`, `body.ncells == 2` for `NODES_2TRI`/`CELLS_2TRI` |
| Strength shape | `size(body.strength) == (2, 1)` for single element type |
| Default CPoffset | `body.CPoffset == 1e-14` |
| Default watertight | `body.watertight == false` |
| VTK cells | `length(body.vtk_cells) == 2` |

### 4.2 Element type parameterization

| Property | Test |
|----------|------|
| ConstantSource | `N == 1`, strength `(ncells, 1)` |
| ConstantDoublet | `N == 1`, strength `(ncells, 1)` |
| Source+Doublet union | `N == 2`, strength `(ncells, 2)` |

### 4.3 Dirichlet BC flag

| Property | Test |
|----------|------|
| DBC=true | `pnl.has_dirichlet_bc(body) == true` |
| DBC=false (default) | `pnl.has_dirichlet_bc(body) == false` |

### 4.4 From GridTriangleSurface

**Setup:** Create a small GT grid (e.g., `gt.Grid` → `gt.GridTriangleSurface`).

| Property | Test |
|----------|------|
| Neighbor populated | `body.neighbor` has correct shape and nonzero entries |
| Node count | Matches grid node count |
| Cell count | Matches grid cell count |

### 4.5 `generate_loft` smoke test

| Property | Test |
|----------|------|
| Runs without error | Generates a body from a simple cross-section loft |
| Valid geometry | `body.ncells > 0`, normals are unit length |

### 4.6 `generate_revolution` smoke test

| Property | Test |
|----------|------|
| Runs without error | Generates a body of revolution from a profile |
| Closed surface | If full 360° revolution, body should be watertight |

---

## 5. Solvers  (P0)

**Source:** `src/FLOWPanel_solver.jl`

### 5.1 `Backslash` dispatch

| Property | Test |
|----------|------|
| NonLiftingBody → Neumann | `typeof(Backslash(nlb))` is `BackslashNeumann` |
| RigidWakeBody → Dirichlet | `typeof(Backslash(rwb))` is `BackslashDirichlet` |

### 5.2 `BackslashNeumann` construction

**Setup:** ~8-panel NonLiftingBody{ConstantSource} (octahedron).

| Property | Test |
|----------|------|
| G matrix shape | `size(solver.G) == (ncells, ncells)` |
| G nonzero diagonal | `all(diag(solver.G) .!= 0)` |
| LU computed | `solver.Glu !== nothing` |
| RHS shape | `length(solver.rhs) == ncells` |

### 5.3 `BackslashNeumann` solve

**Setup:** Octahedron body, freestream `Uinf = [1, 0, 0]`.

| Property | Test |
|----------|------|
| Nonzero strengths | `any(body.strength .!= 0)` after solve |
| BC satisfaction | At each CP: `|dot(U_induced + U_inf, normal)| < 0.1` (coarse mesh tolerance) |

### 5.4 `BackslashDirichlet` construction

**Setup:** Small NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2}(nodes, cells; DBC=true).

| Property | Test |
|----------|------|
| Interior CPs | Control points are offset inward (negative of normal direction from centroid) |
| G matrix shape | `size(solver.G) == (ncells, ncells)` |

### 5.5 `BackslashDirichlet` solve

| Property | Test |
|----------|------|
| Nonzero doublet strengths | After solve, doublet component of `body.strength` is nonzero |
| Source strengths set | Source component = `-dot(Uinf, normal)` for each panel |

### 5.6 `KrylovSolver` construction

| Property | Test |
|----------|------|
| Default params | `solver.method == :gmres`, `solver.itmax == 20` |
| Custom params | `KrylovSolver(body; atol=1e-8, rtol=1e-8)` stores those values |
| Backend stored | `solver.backend` is the supplied backend |

### 5.7 `FGSSolver` construction

| Property | Test |
|----------|------|
| Single body | `FGSSolver(body)` succeeds |
| Multi-body tuple | `FGSSolver((body1, body2))` is unsupported |
| Default params | `solver.max_iterations == 100`, `solver.tolerance == 1e-6` |

### 5.8 Multi-body `solve!`

**Setup:** Two small NonLiftingBody objects offset in space. Backslash solvers for each.

| Property | Test |
|----------|------|
| Both solved | Both bodies have nonzero strengths |
| Runs without error | Completes without exceptions |

### 5.9 `calc_elprescribe`

| Property | Test |
|----------|------|
| ConstantSource NLB | Returns `[]` (empty) |
| VortexRing watertight NLB | Returns `[(1, 0.0)]` |
| VortexRing non-watertight NLB | Returns `[]` |

### 5.10 `numtype`

| Property | Test |
|----------|------|
| Float64 body | `pnl.numtype(body) == Float64` |
| Float32 body | `pnl.numtype(body) == Float32` |

---

## 6. RigidWakeBody  (P0)

**Source:** `src/FLOWPanel_liftingbody.jl`, `src/FLOWPanel_abstractliftingbody.jl`

### 6.1 Construction

**Setup:** 8-panel flat plate mesh with trailing edge shedding edges defined.

| Property | Test |
|----------|------|
| Dimensions | Correct `nnodes`, `ncells` |
| Shedding populated | `length(body.shedding) > 0` |
| Das shape | Each entry has shape `(3, nshed+1)` |
| Strength shape | `(ncells, N)` |

### 6.2 Shedding edge consistency

| Property | Test |
|----------|------|
| Upper/lower pairs | Each shedding surface has matching upper/lower panel counts |
| `shedding_full` indices | Nonzero entries correspond to TE panels |

### 6.3 VortexRing solve

**Setup:** Small flat plate RigidWakeBody{VortexRing} at angle of attack with Das pointing downstream.

| Property | Test |
|----------|------|
| Nonzero circulation | `any(body.strength .!= 0)` |
| Prescribed element | With `elprescribe=[(1, 0.0)]`, `body.strength[1, 1] ≈ 0.0` |

### 6.4 `generate_loft_liftbody` smoke test

**Setup:** Use `pnl.simplewing(...)` or equivalent to create a lifting body via loft.

| Property | Test |
|----------|------|
| Valid body | `body.ncells > 0` |
| Shedding defined | `length(body.shedding) > 0` |
| TE nodes coincident | Upper and lower TE node positions match within `1e-10` |

### 6.5 `extra_reset!`

| Property | Test |
|----------|------|
| TE velocity zeroed | After `reset!(body)`, all `body.velocity_te` entries are zero |

---

## 7. Postprocessing  (P1)

**Source:** `src/FLOWPanel_postprocess.jl`

### 7.1 `calcfield_Cp!`

**Setup:** Body with manually set velocity field and known Uref.

| Property | Test |
|----------|------|
| Stagnation | Velocity = 0 → `Cp = 1.0` |
| Freestream | Velocity magnitude = Uref → `Cp = 0.0` |
| Accelerated | Velocity = 2×Uref → `Cp = 1 - 4 = -3.0` |
| Unsteady term | With `phi_dot[i] = d`, adds `-2d / Uref²` to Cp |

### 7.2 `calcfield_F!`

**Setup:** Body with known Cp, areas, normals.

| Property | Test |
|----------|------|
| Single panel force | `Cp=1, area=2, rho=1, Uinf=10, normal=[0,0,1]` → `F = -Cp * 0.5 * rho * |Uinf|² * area * normal = [0, 0, -100]` |
| Force direction | Force is along `-Cp * normal` |

### 7.3 `calcfield_Ftot!`

| Property | Test |
|----------|------|
| Sum of forces | Total force equals column-wise sum of per-panel forces |
| Zero forces | All-zero per-panel forces → zero total force |

### 7.4 `calcfield_Mtot!`

**Setup:** Single panel at `r = [1, 0, 0]`, force `F = [0, 0, 1]`, reference point at origin.

| Property | Test |
|----------|------|
| Cross product | Moment = `cross([1,0,0], [0,0,1]) = [0, -1, 0]` |
| Shifted reference | Moment about `[1, 0, 0]` is zero (force applied at reference point) |

### 7.5 `calcfield_LDS!`

**Setup:** Total force `F = [1, 0, 0]`, with `Lhat = [0, 0, 1]`, `Dhat = [1, 0, 0]`.

| Property | Test |
|----------|------|
| Decomposition | Lift = 0, Drag = 1, Side = 0 |
| Orthogonal force | `F = [0, 0, 1]` → Lift = 1, Drag = 0 |

### 7.6 `compute_mu_gradient!`

**Setup:** Small mesh with linearly varying doublet strength `mu[i] = x_centroid[i]`.

| Property | Test |
|----------|------|
| Gradient recovery | Interior panels: recovered gradient ≈ `[1, 0, 0]` (within 20% for coarse mesh) |
| Surface constraint | Normal component of gradient near zero |

### 7.7 `calcfield_U!`

**Setup:** Solved NonLiftingBody with known source strengths.

| Property | Test |
|----------|------|
| Velocity populated | After call, `body.velocity` has nonzero entries |
| "U" field added | `"U" in body.fields` |

---

## 8. Reference Frames  (P1)

**Source:** `src/FLOWPanel_frames.jl`

### 8.1 `Rodrigues` rotation matrix

| Property | Test |
|----------|------|
| Identity | `Rodrigues([0,0,1], 0.0) ≈ I(3)` |
| 90° about z | Maps `[1,0,0]` to `[0,1,0]` or `[0,-1,0]` (verify sign convention) |
| 360° | `Rodrigues([0,0,1], 2π) ≈ I(3)` |
| Composition | `Rodrigues(ax, a) * Rodrigues(ax, b) ≈ Rodrigues(ax, a+b)` for same axis |

### 8.2 `inverse_Rodrigues`

| Property | Test |
|----------|------|
| Round-trip | `inverse_Rodrigues(Rodrigues(axis, angle))` recovers `(axis, angle)` within tolerance |
| Identity matrix | Returns angle ≈ 0 |

### 8.3 `ReferenceFrame` construction

| Property | Test |
|----------|------|
| Root frame | `parent_index == -1`, default identity rotation |
| Custom origin | `frame.x == supplied origin` |
| Custom velocity | `frame.v == supplied velocity` |

### 8.4 `add_frame!`

| Property | Test |
|----------|------|
| Frame count | `length(frames) == 2` after adding one child |
| Parent linkage | Parent's `child_index` contains child index |
| Dependent propagation | Parent's `dependent_index` includes child's surface indices |

### 8.5 `propagate_kinematics!` — pure translation

**Setup:** Body at origin, frame velocity `v = [1, 0, 0]`, `dt = 0.1`.

| Property | Test |
|----------|------|
| Node displacement | All nodes translated by `[0.1, 0, 0]` |
| Frame origin update | `frame.x ≈ [0.1, 0, 0]` |

### 8.6 `propagate_kinematics!` — pure rotation

**Setup:** Body centered at origin, frame with `ω_axis = [0,0,1]`, `ω = π/2`, `dt = 1.0`.

| Property | Test |
|----------|------|
| Node rotation | Nodes are rotated 90° about z-axis |
| Frame R update | `frame.R` reflects the 90° rotation |

### 8.7 `kinematic_velocity!`

**Setup:** Body with frame moving at constant velocity `[10, 0, 0]`.

| Property | Test |
|----------|------|
| Velocity subtraction | Each CP velocity includes `[-10, 0, 0]` component (opposing rigid body motion) |

---

## 9. Utilities  (P1)

**Source:** `src/FLOWPanel_utils.jl`

### 9.1 Vector operations

| Property | Test |
|----------|------|
| `dot` | `pnl.dot([1,2,3], [4,5,6]) == 32` |
| `norm` | `pnl.norm([3,4,0]) == 5.0` |
| `cross` | `pnl.cross([1,0,0], [0,1,0]) == [0,0,1]` |
| `cross!` | Writes result into pre-allocated array, matches `cross` result |

### 9.2 `decompose!`

**Setup:** Identity basis `ihat=[1,0,0], jhat=[0,1,0], khat=[0,0,1]`.

| Property | Test |
|----------|------|
| Identity basis | `decompose!([1,2,3], I_hat, J_hat, K_hat)` → `[1, 2, 3]` |
| Rotated basis | Projections onto rotated basis give correct components |

### 9.3 `simplewing`

| Property | Test |
|----------|------|
| Returns valid body | `body.ncells > 0` |
| Correct span | Wing half-span matches `b/2` parameter |
| Shedding edges | Trailing edge shedding is defined |

### 9.4 `slicefield`

**Setup:** Solved body with a populated scalar field.

| Property | Test |
|----------|------|
| Correct extraction | Returns values only for panels matching the slice criteria |
| Empty result | Slice at non-existent location returns empty |

---

## 10. Wake Model  (P1)

**Source:** `src/FLOWPanel_wake.jl`

### 10.1 `PanelWake` construction

**Setup:** Small RigidWakeBody with shedding.

| Property | Test |
|----------|------|
| Zero initial rows | `wake.nwakes[] == 0` |
| Correct node shape | `size(wake.nodes[i]) == (3, nwakerows+1, nshed+1)` per surface |
| Correct strength shape | `size(wake.strength[i]) == (kernel_dim, nwakerows, nshed)` |

### 10.2 `reset!`

| Property | Test |
|----------|------|
| All velocities zero | `all(v .== 0 for v in wake.velocity)` |

### 10.3 `apply_freestream!`

| Property | Test |
|----------|------|
| Uniform addition | After `apply_freestream!(wake, [10, 0, 0])`, all populated velocity entries include `[10, 0, 0]` |

### 10.4 `propagate!`

**Setup:** Wake with 1 shed row, known velocity at each node.

| Property | Test |
|----------|------|
| Node displacement | `nodes_after[:, i, j] ≈ nodes_before[:, i, j] + velocity[:, i, j] * dt` |

### 10.5 `global_to_matrix_index` and `matrix_to_global_index`

**Setup:** Wake with known structure (2 shedding surfaces, 3 rows, 4 shed edges).

| Property | Test |
|----------|------|
| Round-trip | `matrix_to_global_index(wake, global_to_matrix_index(wake, i)...) == i` for all valid `i` |
| Boundary values | First and last global indices map correctly |

### 10.6 `update_TE!`

**Setup:** Body with shedding, wake.

| Property | Test |
|----------|------|
| First row matches TE | After `update_TE!(wake, body)`, first wake row nodes match body's trailing edge node positions |

### 10.7 `shed_wake!`

**Setup:** Body with solved strengths, wake with 1 existing row.

| Property | Test |
|----------|------|
| Row count incremented | `wake.nwakes[]` increases by 1 |
| Strength set | New row strength equals difference of upper/lower TE panel strengths |

---

## 11. Simulation  (P2)

**Source:** `src/FLOWPanel_simulate.jl`

### 11.1 Single time step smoke test

**Setup:** Small RigidWakeBody + PanelWake, constant freestream `Uinf(t) = [1,0,0]`, identity maneuver.

| Property | Test |
|----------|------|
| Runs without error | `simulate!` completes for `t_range = [0.0, 0.1, 0.2]` |
| Strength populated | Body has nonzero strength after simulation |
| Wake grows | `wake.nwakes[] >= 1` |

### 11.2 Steady-state convergence

**Setup:** Same setup, 10 time steps.

| Property | Test |
|----------|------|
| Force convergence | Per-step total force changes decrease over time |
| Finite values | All Cp, F values are finite (no NaN or Inf) |

### 11.3 Monitor callback

**Setup:** Custom monitor function that records forces at each step.

| Property | Test |
|----------|------|
| Called each step | Number of recorded entries == number of time steps |
| Finite forces | All recorded forces are finite |

---

## 12. Analytical Validation  (P0)

**Source:** Full solver pipeline

### 12.1 Sphere in uniform flow — Neumann (ConstantSource)

**Setup:** Sphere of radius `R=1` generated via `generate_revolution` with ~200 panels. `Uinf = [1, 0, 0]`.

| Property | Test |
|----------|------|
| Stagnation Cp | At forward stagnation point: `Cp ≈ 1.0` within 5% |
| Equator Cp | At θ=90°: `Cp ≈ -1.25` within 15% (coarse mesh) |
| Zero drag | `|Ftot[1]| / (0.5 * ρ * |Uinf|² * π * R²) < 0.05` (d'Alembert's paradox) |
| Zero lift | `|Ftot[3]| / (0.5 * ρ * |Uinf|² * π * R²) < 0.05` |
| Symmetry | Cp distribution is symmetric about the flow axis |

### 12.2 Sphere in uniform flow — Dirichlet (Source+Doublet)

**Setup:** Same sphere, `NonLiftingBody{Union{ConstantSource, ConstantDoublet}}` with `DBC=true`, solved with `BackslashDirichlet`.

| Property | Test |
|----------|------|
| Stagnation Cp | `Cp ≈ 1.0` within 5% |
| Zero drag | `|CD| < 0.05` |
| Agreement with Neumann | Cp distributions agree within 10% |

### 12.3 Sphere with FGSSolver

**Setup:** Same sphere, solved with `FGSSolver`.

| Property | Test |
|----------|------|
| Cp agreement | Within 5% of Backslash solution |
| Convergence | Solver converges within max_iterations |

### 12.4 Flat plate lifting body at angle of attack

**Setup:** Thin flat plate RigidWakeBody{VortexRing} at α=5°, aspect ratio 6. Use `simplewing`.

| Property | Test |
|----------|------|
| Nonzero lift | `CL > 0` |
| CL slope | `dCL/dα ≈ 2π * AR/(AR+2)` within 20% (Prandtl correction for finite wing) |
| Small drag | `|CD| << CL` |

---

## 13. FastMultipole Interface  (P1)

**Source:** `src/FLOWPanel_abstractbody.jl` (FastMultipole overloads), `src/FLOWPanel_nonliftingbody.jl`

### 13.1 `source_system_to_buffer!`

**Setup:** NonLiftingBody with known nodes and strengths.

| Property | Test |
|----------|------|
| Centroid in buffer | Buffer centroid matches mean of 3 cell vertices |
| Radius in buffer | Buffer radius = max distance from centroid to any vertex |
| Strength in buffer | Buffer strength matches `body.strength[i, 1]` |

### 13.2 `get_n_bodies` and `get_position`

| Property | Test |
|----------|------|
| Body count | `FastMultipole.get_n_bodies(body) == body.ncells` |
| Position | `get_position(body, i)` returns control point of cell i |

### 13.3 `has_vector_potential`

| Property | Test |
|----------|------|
| All NonLiftingBody variants | Returns `false` |

### 13.4 `buffer_to_target_system!`

**Setup:** Body with known initial potential and velocity. Buffer with known influence values.

| Property | Test |
|----------|------|
| Potential accumulated | `body.potential[i]` increased by buffer's scalar potential |
| Velocity accumulated | `body.velocity[:, i]` increased by buffer's gradient |

---

## 14. Edge Cases  (P2)

**Source:** Multiple files

### 14.1 Degenerate triangle

**Setup:** Triangle with 3 collinear nodes.

| Property | Test |
|----------|------|
| Zero area | `calc_areas!` returns 0.0 without error |
| Normal handling | `calc_normals!` returns finite values (zero vector or NaN — document which) |

### 14.2 Single-panel body

| Property | Test |
|----------|------|
| Construction | 1-cell body constructs without error |
| Normals | `calc_normals!` produces correct unit normal |
| Areas | `calc_areas!` produces correct area |

### 14.3 Very large coordinates

**Setup:** Triangle with vertices at `(1e6, 1e6, 0)`, `(1e6+1, 1e6, 0)`, `(1e6, 1e6+1, 0)`.

| Property | Test |
|----------|------|
| Area accuracy | `area ≈ 0.5` (no catastrophic cancellation) |
| Normal accuracy | Normal is `[0, 0, 1]` |

### 14.4 Float32 body

**Setup:** NonLiftingBody with `TF=Float32`.

| Property | Test |
|----------|------|
| Construction | Succeeds without error |
| Strength type | `eltype(body.strength) == Float32` |
| Solve | `BackslashNeumann` solve produces finite results |

### 14.5 Empty shedding RigidWakeBody

**Setup:** RigidWakeBody with `noshedding`.

| Property | Test |
|----------|------|
| Construction | Succeeds without error |
| No shedding surfaces | `length(body.shedding) == 0` or shedding matrix has 0 columns |

---

## Recommended File Organization

```
test/
├── runtests.jl                      # existing top-level runner
├── runtests_unit_linearsolver.jl    # Section 1
├── runtests_unit_fmm.jl            # Sections 2, 13
├── runtests_unit_body.jl           # Sections 3, 4
├── runtests_unit_solver.jl         # Section 5
├── runtests_unit_liftingbody.jl    # Section 6
├── runtests_unit_postprocess.jl    # Section 7
├── runtests_unit_frames.jl         # Section 8
├── runtests_unit_utils.jl          # Section 9
├── runtests_unit_wake.jl           # Sections 10, 11
├── runtests_analytical.jl          # Section 12
└── runtests_unit_edgecases.jl      # Section 14
```

Each file should begin with:
```julia
using Test
import FLOWPanel as pnl
import GeometricTools as gt
```

and wrap all tests in a top-level `@testset`:
```julia
@testset verbose=true "Linear Solvers" begin
    # ...
end
```
