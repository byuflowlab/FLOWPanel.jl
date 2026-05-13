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

Thus `PressureLaplace` should finite-difference the reconstructed inertial
velocity `body.velocity + body.velocity_kinematic`, not `body.velocity` alone.
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

and let `ell_ij` be the shared-edge length. A simple symmetric finite-volume
weight is

```math
w_{ij} =
\frac{\ell_{ij}}{d_{ij}} .
```

More elaborate metric factors can be introduced later, but the important
property for the first implementation is symmetry:

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
5. Compute `d_ij` from `body.controlpoints`.
6. Form `w_ij = ell_ij / d_ij`.
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

For each shared edge between panels `i` and `j`, define the center-to-center
unit direction

```math
\hat{\mathbf{e}}_{ij}
=
\frac{\mathbf{x}_j - \mathbf{x}_i}
     {\lVert \mathbf{x}_j - \mathbf{x}_i \rVert}.
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
of panel `i`. The implemented v1 metric approximates that co-normal direction
with the center-to-center direction `\hat{\mathbf{e}}_{ij}`. This is exact only
for an orthogonal panel/dual mesh, but it is consistent with the same
center-to-center metric used by the two-point Laplacian. A more geometrically
faithful RHS would use explicit edge co-normals and dual-cell areas.

The edge contribution is

```math
f_{ij}
=
\ell_{ij}
\left[
\frac{\mathbf{a}_{t,i} + \mathbf{a}_{t,j}}{2}
\right]
\cdot
\hat{\mathbf{e}}_{ij}.
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
implementation simple and conservative across shared edges, but it is less
geometrically faithful than a full finite-volume surface divergence with
dual-face co-normal metrics. Future implementations can replace this RHS with
an area-weighted divergence without changing the sparse pressure operator or
monitor interface.

The first implementation should separate operator assembly from RHS assembly:
the Laplacian depends only on geometry, while the RHS changes with the flow
state.

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
