# Paired-Triangle Surface-Gradient Reconstruction

This note documents the quad-consistent paired-triangle reconstruction used by
FLOWPanel when recovering panel-centered surface gradients from triangular
surface meshes. The method is used most visibly by the doublet or vortex-ring
half-jump contribution in `calcfield_U!`, and by the surface velocity-gradient
path used by `PressureLaplace`.

The primary implementation entry point is
`compute_mu_gradient!(grad_mu, controlpoints, normals, cells, neighbors, mu,
te_info; scale, bad_panel_mask, nodes, quad_consistent,
quad_normal_dot_min)`. The user-facing velocity path is
`calcfield_U!(body_or_bodies, uinf, wakes; gradient_quad_consistent=true)`,
which passes `body.nodes` so that the paired-triangle reconciliation can run.
The vector-field analogue is
`compute_surface_velocity_gradient!(grad_u, u, controlpoints, normals, cells,
neighbors, te_info; nodes, quad_consistent)`.

## Motivation

Many lifting surfaces in FLOWPanel are generated as structured quadrilateral
lofts and then split into triangles. On such meshes, the two triangle
orientations inside each logical quadrilateral can produce alternating local
least-squares stencils. Even when the panel strength `gamma` or doublet
strength `mu` is smooth, the recovered surface gradient may oscillate from one
triangle half to the other. Because the exterior surface velocity contains a
doublet or vortex-ring jump proportional to this gradient, the artifact can
contaminate reconstructed velocity and pressure fields.

The paired-triangle method treats mutually compatible triangle halves as a
single logical quadrilateral only after each triangle has computed its own
local contribution. The method therefore preserves the existing one-ring
least-squares reconstruction, trailing-edge isolation, and bad-panel fallback,
while removing the structured split-direction bias where the mesh topology
strongly indicates a quad pair.

## Local Least-Squares Reconstruction

For a panel `i` with control point `x_i`, unit normal `n_i`, and scalar sample
`mu_i`, FLOWPanel first builds a local stencil `S_i` from neighboring panels.
For each accepted neighbor `j`, define

```math
\mathbf{r}_{ij} = \mathbf{x}_j - \mathbf{x}_i,
\qquad
\Delta\mu_{ij} = \mu_j - \mu_i .
```

The unconstrained least-squares model is

```math
\mathbf{r}_{ij}\cdot \mathbf{g}_i \approx \Delta\mu_{ij},
\qquad j \in S_i,
```

where `g_i` approximates the surface gradient of `mu` at panel `i`. Since the
gradient should lie in the tangent plane, the implementation adds a normal
penalty to enforce

```math
\mathbf{n}_i\cdot \mathbf{g}_i = 0 .
```

The normal equations have the form

```math
\left(
    \sum_{j\in S_i}\mathbf{r}_{ij}\mathbf{r}_{ij}^T
    + \lambda_i \mathbf{n}_i\mathbf{n}_i^T
    + \epsilon_i \mathbf{I}
\right)\mathbf{g}_i
=
\sum_{j\in S_i}\mathbf{r}_{ij}\Delta\mu_{ij}.
```

Here `lambda_i` is scaled by the mean squared stencil distance for conditioning,
and `epsilon_i` is a small Tikhonov regularization used for nearly rank-deficient
local systems.

The sign convention is important. `compute_mu_gradient!` does not accumulate
`grad(mu)` directly. It accumulates

```math
\texttt{grad\_mu}_{:,i} \mathrel{+}= -\,\texttt{scale}\,\nabla_s\mu_i .
```

For the exterior doublet or vortex-ring half-jump in `calcfield_U!`, FLOWPanel
uses `scale=0.5`, so the added velocity contribution is
`-0.5 grad_s(mu)` in the package convention. When
`compute_surface_velocity_gradient!` needs the actual surface gradient of each
velocity component, it calls `compute_mu_gradient!` with `scale=-1.0`.

## Trailing-Edge Isolation

The local stencil must not mix the two sides of a lifting surface at a trailing
edge. Panels marked by `te_info` identify the two local vertices that form the
trailing-edge edge. When a candidate neighbor shares both of those vertices, it
lies across the trailing edge and is rejected from the stencil. This keeps upper
and lower surface gradients one-sided at the edge where the vortex sheet or
doublet strength can be discontinuous.

For trailing-edge panels, FLOWPanel may expand the remaining one-sided stencil
through accepted interior neighbors to improve least-squares conditioning. The
cross-TE rejection remains tied to the marked trailing-edge edge, so the
expanded stencil does not deliberately bridge upper and lower surfaces.

## Pair Detection

Quad-consistent reconciliation runs only when `nodes` is supplied and
`quad_consistent=true`. The public velocity path enables this by default through
`calcfield_U!(...; gradient_quad_consistent=true)`.

For each triangle `i`, FLOWPanel selects at most one candidate triangle `j`.
The selection rules are:

- candidates come only from the mesh `neighbors` table;
- candidates across a marked trailing-edge edge are rejected;
- candidates must satisfy the normal-alignment threshold
  `n_i dot n_j >= quad_normal_dot_min`, with default
  `quad_normal_dot_min = cos(pi / 4)`;
- among the remaining candidates, the selected neighbor is the one across the
  longest shared edge;
- the pair is accepted only if the selection is mutual: `i` selects `j` and
  `j` selects `i`.

The longest-shared-edge rule identifies the diagonal split of a structured quad
because the two triangle halves share the quad diagonal, which is typically the
longest edge in each triangle. The mutual-pair rule prevents one triangle from
being reconciled with a neighbor that is a better quad candidate for a different
panel. The normal threshold prevents strongly folded or unrelated triangles from
being averaged merely because they share a long edge.

## Paired Reconstruction

The paired method does not change the base least-squares solve. FLOWPanel first
computes each panel's contribution independently, including sign and scale:

```math
\mathbf{q}_i = -\,\texttt{scale}\,\nabla_s\mu_i .
```

For an accepted pair `(i, j)`, it then forms the shared contribution

```math
\bar{\mathbf{q}}_{ij}
=
\frac{1}{2}\left(\mathbf{q}_i + \mathbf{q}_j\right).
```

Because the two triangles may not be exactly coplanar, the averaged vector is
projected back onto each triangle's own tangent plane:

```math
\mathbf{P}_i = \mathbf{I} - \mathbf{n}_i\mathbf{n}_i^T,
\qquad
\mathbf{P}_j = \mathbf{I} - \mathbf{n}_j\mathbf{n}_j^T,
```

```math
\mathbf{q}_i^\star = \mathbf{P}_i\bar{\mathbf{q}}_{ij},
\qquad
\mathbf{q}_j^\star = \mathbf{P}_j\bar{\mathbf{q}}_{ij}.
```

The projected values replace the independent local contributions for the two
triangles. Unpaired panels, rejected trailing-edge candidates, and candidates
that fail the normal-alignment or mutual-selection tests retain their original
least-squares contributions.

After this reconciliation pass, `compute_mu_gradient!` accumulates the local
buffer into the caller's output:

```math
\texttt{grad\_mu} \mathrel{+}= \mathbf{q}^\star .
```

This additive behavior is deliberate. `calcfield_U!` may have already added
wake-induced velocity, body-panel induced velocity, or other contributions
before the doublet-gradient half-jump is applied.

## Pressure-Laplace Path

`PressureLaplace` reconstructs a surface velocity gradient through
`compute_surface_velocity_gradient!`. That routine applies the same scalar
surface-gradient machinery to each velocity component, so trailing-edge
isolation, bad-panel handling, and optional quad-consistent pairing follow the
same rules as `compute_mu_gradient!`.

The output convention is

```math
\texttt{grad\_u}[k,\ell,i] =
\frac{\partial u_k}{\partial x_\ell}
\quad \text{at panel } i,
```

with each row constrained to the local tangent plane. Supplying `nodes` and
leaving `quad_consistent=true` enables the paired-triangle reconciliation for
these component gradients as well, reducing split-direction artifacts before
the pressure Poisson right-hand side is assembled.

## Implementation Notes

The paired-triangle helpers are internal implementation details. Users should
prefer the documented paths:

- call `calcfield_U!(...; gradient_quad_consistent=true)` for the exterior
  surface-velocity reconstruction used in post-processing;
- call `compute_mu_gradient!` directly only when constructing custom scalar
  surface-gradient or half-jump workflows;
- call `compute_surface_velocity_gradient!` for panel-centered gradients of a
  vector velocity field.

The method assumes that `neighbors` reflects true shared-edge triangle
adjacency and that `nodes` and `cells` use the same mesh indexing. It is
therefore a reconciliation of an already valid triangular surface mesh, not a
quad remeshing step.
