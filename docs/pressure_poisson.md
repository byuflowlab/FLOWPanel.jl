# Pressure Poisson Solve From the Euler Equation

This note derives a sparse panel-centered pressure solve intended for FLOWPanel
post-processing. The first implementation target is a symmetric finite-volume
surface Laplacian on panel adjacency, solved with conjugate gradients from
Krylov.jl and a Jacobi preconditioner.

## Euler Pressure Equation

For inviscid incompressible flow with constant density, the Euler momentum
equation is

```math
\rho
\left(
    \frac{\partial \mathbf{u}}{\partial t}
    + \mathbf{u} \cdot \nabla \mathbf{u}
\right)
=
-\nabla p ,
```

with incompressibility constraint

```math
\nabla \cdot \mathbf{u} = 0 .
```

Equivalently, the pressure gradient is

```math
\nabla p
=
-\rho
\left(
    \frac{\partial \mathbf{u}}{\partial t}
    + \mathbf{u} \cdot \nabla \mathbf{u}
\right).
```

Taking the divergence gives the pressure Poisson equation

```math
\nabla^2 p
=
-\rho
\nabla \cdot
\left(
    \frac{\partial \mathbf{u}}{\partial t}
    + \mathbf{u} \cdot \nabla \mathbf{u}
\right).
```

For exactly incompressible flow,
`\nabla \cdot \partial \mathbf{u}/\partial t =
\partial(\nabla \cdot \mathbf{u})/\partial t = 0`, so the continuous equation
can also be written as

```math
\boxed{
\nabla^2 p
=
-\rho
\nabla \cdot
\left(
    \mathbf{u} \cdot \nabla \mathbf{u}
\right)
}
```

or, in index notation,

```math
\nabla^2 p
=
-\rho
\frac{\partial u_j}{\partial x_i}
\frac{\partial u_i}{\partial x_j}.
```

In a discrete unsteady panel calculation, the time-derivative contribution may
not vanish exactly because the discrete velocity field is not exactly
divergence-free. The implementation should therefore allow the right-hand side
to include either only the convective source or the full material-acceleration
source,

```math
s =
-\rho \nabla \cdot \mathbf{a},
\qquad
\mathbf{a}
=
\frac{\partial \mathbf{u}}{\partial t}
+ \mathbf{u} \cdot \nabla \mathbf{u}.
```

The pressure is defined only up to an additive constant. A gauge condition is
required before the discrete system is nonsingular.

## Velocity and Acceleration at Moving Control Points

The Euler equation must be evaluated with the inertial-frame fluid velocity
`\mathbf{u}(\mathbf{x}, t)`. In FLOWPanel, `body.velocity` is the apparent fluid
velocity at the moving control points after the rigid-body control-point
velocity has been subtracted by `kinematic_velocity!`. The subtracted
rigid-body velocity is stored in `body.velocity_kinematic`. Therefore the
inertial-frame fluid velocity at the current control points is

```math
\mathbf{u}_i^n
=
\texttt{body.velocity}_{:,i}^n
+
\texttt{body.velocity\_kinematic}_{:,i}^n .
```

The implemented surface solve stores a tangent-projected version in the
rolling finite-difference buffer,

```math
\mathbf{u}_{\mathrm{buf},i}^n
=
\left(I - \mathbf{n}_i\mathbf{n}_i^T\right)
\texttt{body.velocity}_{:,i}^n
+
\texttt{body.velocity\_kinematic}_{:,i}^n .
```

This removes solver-tolerance normal leakage before time differencing. In the
continuous impermeable limit it equals the reconstructed inertial velocity on
the surface.

Let the current control point move with known grid/body velocity

```math
\mathbf{w}_i^n
=
\texttt{body.velocity\_kinematic}_{:,i}^n .
```

A finite difference at the same panel index does not approximate the Eulerian
partial derivative directly. It approximates the derivative following the
moving control point:

```math
\frac{\mathbf{u}_i^n - \mathbf{u}_i^{n-1}}{\Delta t}
\approx
\left.
\frac{d}{dt}\mathbf{u}(\mathbf{x}_i(t), t)
\right|_{\mathrm{grid}}
=
\frac{\partial \mathbf{u}}{\partial t}
+
\mathbf{w}_i^n \cdot \nabla \mathbf{u}.
```

The acceleration required by Euler is the fluid material acceleration,

```math
\mathbf{a}_i^n
=
\frac{\partial \mathbf{u}}{\partial t}
+
\mathbf{u}_i^n \cdot \nabla \mathbf{u}.
```

Combining the two expressions gives the appropriate discrete acceleration from
values sampled at moving body control points:

```math
\boxed{
\mathbf{a}_i^n
=
\frac{\mathbf{u}_i^n - \mathbf{u}_i^{n-1}}{\Delta t}
+
\left[
\left(\mathbf{u}_i^n - \mathbf{w}_i^n\right)
\cdot \nabla
\right]\mathbf{u}^n
}
```

Thus `PressureLaplace` should finite-difference the reconstructed surface
inertial velocity `tangent(body.velocity) + body.velocity_kinematic`, not
`body.velocity` alone.
The convective correction should use the velocity relative to the moving
control-point grid, `u - w`, while the velocity gradients should be gradients
of the inertial velocity `u`. On an impermeable body surface, `u - w` is the
tangential slip velocity, so the surface form uses

```math
\nabla_s p
=
-\rho \mathbf{P}_s \mathbf{a},
\qquad
\mathbf{P}_s = \mathbf{I} - \mathbf{n}\mathbf{n}^T,
```

and the panel-centered right-hand side should approximate

```math
\Delta_s p
=
-\rho \nabla_s \cdot \left(\mathbf{P}_s \mathbf{a}\right).
```

The sparse operator used below is assembled in the symmetric positive
semidefinite `-\Delta_s` convention, so the solved form is equivalently

```math
-\Delta_s p
=
\rho \nabla_s \cdot \left(\mathbf{P}_s \mathbf{a}\right).
```

## Panel-Centered Surface Approximation

The desired field is one pressure value per panel. Let panel `i` have pressure
`p_i`, area `A_i`, control point `x_i`, outward unit normal `n_i`, and neighbor
set `N(i)` from `body.neighbor`. The pressure unknown vector is

```math
\mathbf{p} =
\begin{bmatrix}
p_1 & p_2 & \cdots & p_N
\end{bmatrix}^T ,
```

where `N = body.ncells`.

The sparse solve approximates a surface Poisson equation over the panel mesh,

```math
-\Delta_s p = s ,
```

at the panel control points. This is not a replacement for the full volumetric
Poisson problem around the body. It is a surface reconstruction that uses the
Euler pressure equation to obtain a pressure field consistent with the
available panel velocity and acceleration data.

For each panel pair `(i, j)` sharing an edge, define

```math
\mathbf{r}_{ij} = \mathbf{x}_j - \mathbf{x}_i,
\qquad
d_{ij} = \lVert \mathbf{r}_{ij} \rVert,
```

and let `ell_ij` be the shared-edge length. The panel-centered finite-volume
metric uses the outward surface co-normal from panel `i` toward panel `j`.
With shared-edge tangent `\hat{\mathbf{t}}_{ij}` and averaged unit normal

```math
\bar{\mathbf{n}}_{ij}
=
\frac{\mathbf{n}_i + \mathbf{n}_j}
     {\lVert \mathbf{n}_i + \mathbf{n}_j \rVert},
\qquad
\boldsymbol{\nu}_{ij}
=
\operatorname{orient}_{i\to j}
\frac{\hat{\mathbf{t}}_{ij} \times \bar{\mathbf{n}}_{ij}}
     {\lVert \hat{\mathbf{t}}_{ij} \times \bar{\mathbf{n}}_{ij} \rVert},
```

where the sign is chosen so
`\boldsymbol{\nu}_{ij}\cdot\mathbf{r}_{ij} > 0`, the symmetric two-point
finite-volume weight is

```math
w_{ij}
=
\ell_{ij}
\frac{\boldsymbol{\nu}_{ij}\cdot\mathbf{r}_{ij}}
     {\lVert\mathbf{r}_{ij}\rVert^2}.
```

This is not the standard vertex cotangent Laplace-Beltrami operator. Pressure
unknowns in FLOWPanel remain panel-centered, so the operator uses the
panel-center dual distance projected onto the shared-edge co-normal. The
important algebraic property is symmetry:

```math
w_{ij} = w_{ji} \ge 0 .
```

The integrated finite-volume equation on panel `i` is

```math
\sum_{j \in N(i)}
w_{ij} (p_i - p_j)
=
q_i .
```

This left-hand side is the panel-integrated `-\Delta_s p` approximation. Its
sign convention is important: the RHS must approximate
`\rho \nabla_s \cdot \mathbf{a}_t`, not
`-\rho \nabla_s \cdot \mathbf{a}_t`, where
`\mathbf{a}_t = \mathbf{P}_s \mathbf{a}`.

This produces the sparse linear system

```math
L \mathbf{p} = \mathbf{b},
```

with

```math
b_i = q_i ,
```

and matrix entries

```math
L_{ij} =
\begin{cases}
 \sum_{k \in N(i)} w_{ik}, & i = j, \\
 -w_{ij},                 & j \in N(i), \\
 0,                       & \text{otherwise}.
\end{cases}
```

This matrix is symmetric positive semidefinite. The constant vector is in the
nullspace because only pressure differences appear:

```math
L \mathbf{1} = \mathbf{0}.
```

After applying a pressure gauge, the operator becomes positive definite on the
remaining unknowns and is suitable for conjugate gradients.

## Gauge Condition

The simplest gauge is to pin one reference panel:

```math
p_r = p_\mathrm{ref}.
```

For a default gauge pressure solve, `p_ref = 0` is sufficient. In the assembled
matrix, this can be imposed by replacing row and column `r` with the identity
constraint,

```math
L_{rr} = 1,
\qquad
L_{rj} = L_{jr} = 0 \quad (j \ne r),
\qquad
b_r = p_\mathrm{ref}.
```

Symmetry must be preserved when pinning the matrix. If only the row is replaced,
the matrix becomes nonsymmetric and CG is no longer mathematically appropriate.

A zero-mean gauge,

```math
\sum_i A_i p_i = 0,
```

is often more physically neutral, but it requires either a constrained solve or
a projection onto the mean-free subspace. Pinning one panel is the recommended
first implementation because it keeps the sparse system directly SPD and easy
to use with `Krylov.cg`.

## Sparse Matrix Formation

The topology is already available from `body.neighbor`, where each triangular
panel has up to three edge neighbors. Matrix formation can therefore be done in
one pass over panels and local edges:

1. For each panel `i`, loop over its three local edges.
2. Read neighbor `j = body.neighbor[e, i]`.
3. Skip `j <= 0` boundary edges.
4. Compute the shared edge length `ell_ij` from `body.cells` and `body.nodes`.
5. Compute the averaged normal and edge co-normal `nu_ij`.
6. Form `w_ij = ell_ij * dot(nu_ij, x_j - x_i) / norm(x_j - x_i)^2`.
7. Add `+w_ij` to `L[i, i]` and `-w_ij` to `L[i, j]`.

If every interior edge is visited from both adjacent panels, the weights must be
computed identically from both sides to preserve symmetry. A safer assembly is
edge-based: visit each undirected shared edge once, compute `w_ij`, and insert
the four symmetric contributions

```math
L_{ii} \mathrel{+}= w_{ij},
\qquad
L_{jj} \mathrel{+}= w_{ij},
\qquad
L_{ij} \mathrel{-}= w_{ij},
\qquad
L_{ji} \mathrel{-}= w_{ij}.
```

The edge-based form is preferable because symmetry is exact up to floating-point
roundoff and the nonzero pattern is independent of panel visitation order.

For repeated solves on a fixed mesh, the sparse structure and weights should be
cached. If the body moves rigidly without changing panel shape or connectivity,
the weights are unchanged because edge lengths and control-point distances are
unchanged. If the mesh deforms, the weights must be updated.

## Right-Hand Side Formation

The implemented monitor builds the right-hand side as a conservative
edge-integrated source rather than explicitly computing a panel-centered
divergence and multiplying by panel area. The acceleration is computed from the
available panel velocity data,

```math
\mathbf{a}_i
=
\frac{\partial \mathbf{u}_i}{\partial t}
+
(\mathbf{u}_{t,i} \cdot \nabla_s)\mathbf{u}_i ,
```

where the unsteady term is obtained by finite differencing the current and
previous monitor-call velocities. This acceleration is then projected into the
local tangent plane,

```math
\mathbf{a}_{t,i}
=
\left( I - \mathbf{n}_i \mathbf{n}_i^T \right)
\mathbf{a}_i .
```

This RHS is a finite-volume divergence, not a pointwise divergence evaluated at
the panel center. For panel `i`,

```math
\int_{A_i} \nabla_s \cdot \mathbf{a}_t \, dA
=
\oint_{\partial A_i} \mathbf{a}_t \cdot \boldsymbol{\nu}_i \, d\ell
\approx
\sum_{e \in \partial i}
\ell_e \, \mathbf{a}_{t,e} \cdot \boldsymbol{\nu}_{i,e},
```

where `\boldsymbol{\nu}_{i,e}` is the outward surface co-normal across edge
`e`: tangent to the surface, perpendicular to the shared edge, and pointing out
of panel `i`. The implemented metric uses the same averaged-normal co-normal
`\boldsymbol{\nu}_{ij}` as the sparse operator, rather than approximating the
co-normal with the center-to-center direction. This removes the orthogonal-dual
mesh assumption from the original TPFA form while keeping the pressure field
panel-centered.

The edge contribution is

```math
f_{ij}
=
\ell_{ij}
\left[
\frac{\mathbf{a}_{t,ij,i} + \mathbf{a}_{t,ij,j}}{2}
\right]
\cdot
\boldsymbol{\nu}_{ij}.
```

Here each side is projected with the averaged edge normal,

```math
\mathbf{a}_{t,ij,i}
=
\left(I - \bar{\mathbf{n}}_{ij}\bar{\mathbf{n}}_{ij}^T\right)\mathbf{a}_i,
\qquad
\mathbf{a}_{t,ij,j}
=
\left(I - \bar{\mathbf{n}}_{ij}\bar{\mathbf{n}}_{ij}^T\right)\mathbf{a}_j.
```

The RHS is accumulated with equal and opposite updates on the two adjacent
panels,

```math
b_i \mathrel{+}= \rho f_{ij},
\qquad
b_j \mathrel{-}= \rho f_{ij}.
```

This accumulates an edge-integrated approximation to
`\rho \nabla_s \cdot \mathbf{a}_t` directly into `b`, matching the `-\Delta_s`
sign convention of `L`. The midpoint acceleration is the face value of
`\mathbf{a}_t`; using the difference
`\mathbf{a}_{t,j} - \mathbf{a}_{t,i}` would instead apply a graph derivative to
the acceleration and would not represent the finite-volume flux of
`\mathbf{a}_t`. Panel areas are not used in this v1 RHS. This keeps the
implementation simple and conservative across shared edges. Future
implementations can add explicit dual-cell areas or a vertex-based cotangent
formulation, but that would be a separate pressure discretization rather than a
drop-in replacement for the current panel-centered unknowns.

The first implementation should separate operator assembly from RHS assembly:
the Laplacian depends only on geometry, while the RHS changes with the flow
state.

## Analytic Velocity Gradient via Multipole Hessian

### Motivation

The convective term `(\mathbf{u}_{\text{rel}} \cdot \nabla)\mathbf{u}_{\text{inertial}}`
requires the spatial Jacobian `\nabla\mathbf{u}` at every panel control point.
The v1 monitor estimated `\nabla\mathbf{u}` by a weighted least-squares fit of
the velocity over each panel's 1-ring (with a 2-ring fallback at trailing
edges and boundary panels). On Delaunay-ish triangulations of smooth surfaces
this works well, but on the anisotropic triangulations that arise from rotor
blades, swept lifting surfaces, or any geometry with high local aspect
ratio, the least-squares fit becomes rank-deficient *in the tangent plane*:
the three face-neighbors of a chordwise-elongated panel give nearly parallel
displacement vectors, so the spanwise component of the gradient is
controlled by the Tikhonov noise floor rather than by data. Pressure on a
blade varies predominantly chordwise while panels are aligned chordwise —
the worst-case orientation. The dominant LS error then enters
`(\mathbf{u} \cdot \nabla)\mathbf{u}`, contaminating the surface Poisson
RHS at exactly the locations where pressure resolution matters most.

### Decomposition

The fluid velocity at a moving panel control point splits into induced and
kinematic parts:

```math
\mathbf{u}_{\text{inertial}}(\mathbf{x})
=
\mathbf{u}_{\text{induced}}(\mathbf{x})
+
\mathbf{u}_{\text{kinematic}}(\mathbf{x}),
\qquad
\mathbf{u}_{\text{kinematic}}(\mathbf{x})
=
\mathbf{U}_O + \boldsymbol{\Omega} \times (\mathbf{x} - \mathbf{x}_O).
```

The kinematic part is rigid-body, so its spatial Jacobian is the constant
skew-symmetric tensor of the body's angular velocity:

```math
\nabla \mathbf{u}_{\text{kinematic}}
=
[\boldsymbol{\Omega}]_\times
=
\begin{pmatrix}
0 & -\Omega_z & \Omega_y \\
\Omega_z & 0 & -\Omega_x \\
-\Omega_y & \Omega_x & 0
\end{pmatrix}.
```

The induced part `\mathbf{u}_{\text{induced}} = \nabla\phi + \nabla\times\mathbf{A}`
is the contribution of every source, doublet, vortex-ring, and vortex
filament in the simulation (bodies and wakes). Its spatial Jacobian is
exactly the same quantity that FastMultipole computes as a natural by-product
of the multipole/local expansion to one higher derivative order — the
"velocity gradient" or "Hessian" output in the FastMultipole compatibility
interface. By telling FastMultipole to populate this field at body
control points, we obtain `\nabla \mathbf{u}_{\text{induced}}` analytically
and at no extra FMM tree-build cost.

```math
\boxed{\;
\nabla \mathbf{u}_{\text{inertial}}
=
\underbrace{\nabla \mathbf{u}_{\text{induced}}}_{\text{FastMultipole Hessian}}
+
\underbrace{[\boldsymbol{\Omega}]_\times}_{\text{from frame } \omega}
\;}
```

Both terms are analytic, mesh-independent, and bounded uniformly in panel
aspect ratio.

### Integration into the simulation

A monitor that needs the body's velocity gradient declares its requirement
via the trait

```julia
monitor_requires_body_hessian(::PressureLaplace) = true
```

Before the time loop, `simulate!` walks the registered monitors and sets

```julia
sys.needs_velocity_gradient[] = true
```

on every body in the simulation when any monitor requires `\nabla\mathbf{u}`.
The dispatch

```julia
requires_hessian(b::AbstractBody) = b.needs_velocity_gradient[]
```

then routes the per-target Hessian flag through the existing tuple-keyed
`influence!` calls,

```julia
influence!(targets, sources, backend;
           gradient=true,
           hessian = Tuple(requires_hessian(sys) for sys in targets))
```

so the wake-on-body and body-on-body FastMultipole passes each accumulate
the analytic `\nabla\mathbf{u}_{\text{induced}}` into the per-body
`body.velocity_gradient::Array{Float64,3}` (shape `3 × 3 × n_{\text{cells}}`).
Bodies that don't need the field pay no extra cost because the type-parameter
`HS` is false for them, and the kernel branches are eliminated at compile time.

In parallel, `kinematic_velocity!` accumulates the body's net angular
velocity `\boldsymbol{\Omega}` (summed over every ancestor frame in the
kinematic chain, in global coordinates) into `body.angular_velocity`. The
monitor reads both at each step:

```julia
G_ind = body.velocity_gradient            # 3×3×ncells
Ω     = body.angular_velocity             # length-3
# convective term per panel:
acc[:,i] += (u_t · ∇) u_inertial
         =  G_ind[:,:,i] * u_t  +  Ω × u_t
```

The cross product encodes `[\boldsymbol{\Omega}]_\times \mathbf{u}_t` without
materialising the skew matrix.

### Per-source contributions

Each source kernel must implement the `velocity_gradient` (`HS=true`) branch
of its `_induced` / `compute_*` method. The kernels currently covered are:

- **ConstantSource panel** — closed-form Hess–Smith Hessian of the panel
  scalar potential (`compute_source_dipole, ::ConstantSource`, `GS` branch).
- **ConstantDoublet panel** — same family, with the doublet integrand.
- **`Union{ConstantSource, ConstantDoublet}`** — superposition of the two.
- **VortexRing filament** — the Biot–Savart segment Jacobian, evaluated by
  `_bound_vortex_gradient`. This function uses ForwardDiff dual numbers
  seeded with the spatial-position partials to differentiate
  `_bound_vortex_velocity` directly, which guarantees that the gradient
  inherits whatever core regularisation the velocity uses (currently
  Vatistas n=2) — a mismatch here would otherwise introduce O(1) errors
  near filament endpoints.
- **`Union{ConstantSource, VortexRing}`** — superposition of source and
  vortex-ring contributions.

Vortex-sheet kernels (`ConstantVortexSheet`, `UniformVortexSheet`) do not
yet have a `velocity_gradient` branch and are not supported by the new
`PressureLaplace` path. The monitor will throw at construction if it is
attached to a body whose kernel mix lacks the analytic gradient.

### Computational cost

Asking FastMultipole to populate the velocity gradient raises the
per-influence-pass cost by roughly 10–20% at modest expansion orders: the
local expansion is differentiated once more and the per-source direct path
writes nine extra buffer entries. There is no additional FMM tree build and
no additional `influence!` pass — the gradient piggybacks on the same two
calls (wake-on-body, body-on-body) that compute body velocity. The least-
squares gradient and its scratch storage (`ATA`, `ATb`, `stencil`,
per-component buffers) are removed entirely.

### Failure modes

`\nabla\mathbf{u}_{\text{induced}}` is finite everywhere away from vortex
filaments and source-panel singularities. Near a filament — for example a
control point lying within `core_size` of a vortex-ring edge of an adjacent
panel — the Vatistas-regularised Biot–Savart Jacobian degrades to the
regularised value (smooth and bounded) rather than diverging as `1/r^2`.
Because the gradient is derived from the regularised velocity by automatic
differentiation, the regularisation in `\nabla\mathbf{u}` matches the
regularisation in `\mathbf{u}` exactly — the analytic `(\mathbf{u}\cdot\nabla)\mathbf{u}`
operator remains the spatial derivative of `\mathbf{u}` for all targets.

## Storage And Complexity

For a triangular closed surface with `N` panels:

- each panel has at most three off-diagonal neighbors;
- each row has at most four nonzeros, including the diagonal;
- the sparse matrix has approximately `nnz = 4N`.

Using `SparseMatrixCSC{Float64,Int64}`, storage is approximately

```math
8\,\mathrm{bytes} \cdot \mathrm{nnz}
+ 8\,\mathrm{bytes} \cdot \mathrm{nnz}
+ 8\,\mathrm{bytes} \cdot (N + 1),
```

for values, row indices, and column pointers. With `nnz = 4N`, this is roughly

```math
72N \ \mathrm{bytes},
```

not including solver work vectors. CG needs a small fixed number of vectors,
so total solver memory remains linear in `N`.

Assembly is `O(N)` for triangular meshes. Each sparse matrix-vector product is

```math
O(\mathrm{nnz}) \approx O(N),
```

and is usually memory-bandwidth dominated. The solve cost is

```math
O(k\,\mathrm{nnz}) \approx O(kN),
```

where `k` is the number of CG iterations.

## Solver And Preconditioner

The default solver should be conjugate gradients from Krylov.jl:

```julia
workspace = Krylov.krylov_workspace(Val(:cg), L, b)
Krylov.krylov_solve!(workspace, L, b; M, ldiv=true, atol, rtol, itmax)
p .= workspace.x
```

This choice relies on the matrix being symmetric positive definite after gauge
fixing. Any implementation change that breaks symmetry, such as row-only
pinning or one-sided least-squares correction, should switch the solver to a
nonsymmetric method such as GMRES.

The default preconditioner should be Jacobi scaling:

```math
M^{-1}_{ii} = \frac{1}{L_{ii}} .
```

Jacobi has `O(N)` setup, `O(N)` storage, and `O(N)` application cost. It is
cheap, deterministic, dependency-free, and compatible with CG. It is not as
strong as incomplete Cholesky or algebraic multigrid, but those methods add
implementation and dependency cost. Jacobi is the appropriate first
preconditioner for a sparse surface Laplacian.

## Performance Notes

For performance-sensitive use:

- cache the sparse matrix and Jacobi inverse diagonal when geometry is fixed;
- update only the RHS when only the flow state changes;
- assemble from undirected edges to avoid duplicate work and preserve symmetry;
- preallocate triplet arrays or reusable CSC buffers when rebuilding is needed;
- keep the pressure solve panel-centered so `body.P`, force integration, and VTK
  cell output do not require projection;
- avoid normal-equation formulations because they square the condition number;
- monitor CG iteration count, since total cost is dominated by repeated sparse
  matrix-vector products.

The intended future post-processing interface is a solve that writes the
resulting pressure directly into `body.P`, accepts velocity or acceleration data
for RHS construction, and optionally reuses cached operator and preconditioner
state across timesteps.
