# Dirichlet Potential Green Identity

This note records the identity checked by `test/dirichlet_potential_test.jl`.

## Continuous Identity

Let

```math
g(\mathbf{x}, \mathbf{y}) = -\frac{1}{4\pi r},
\qquad
r = \lVert \mathbf{x} - \mathbf{y} \rVert .
```

For a harmonic potential `phi` evaluated on a smooth boundary from the exterior
side, Green's identity gives

```math
\frac{1}{2}\phi(\mathbf{x})
- \int_{\partial\Omega}
    \phi(\mathbf{y})
    \frac{\partial}{\partial n_\mathbf{y}}
    \left(-\frac{1}{4\pi r}\right)
    dS_\mathbf{y}
=
-\int_{\partial\Omega}
    \left(-\frac{1}{4\pi r}\right)
    \frac{\partial\phi}{\partial n}(\mathbf{y})
    dS_\mathbf{y}.
```

The left integral is the doublet influence. The right integral is the source
potential induced by the boundary-normal derivative of `phi`.

## Discrete Form

FLOWPanel's Dirichlet body stores source and doublet panel strengths in the two
columns of `body.strength`. In the test, the target potential `phi` is the
wake-induced potential at the body control points.

The doublet influence matrix `G` is assembled using `pnl._G!` with
`update_geometry=false`, so the caller is responsible for choosing the side of
the boundary where the control points are placed. If no wake panels are coupled
to the body-panel unknowns, the discrete identity can be written as

```math
-G\,\phi = \phi_\mathrm{source}.
```

For a `RigidWakeBody`, however, `G` can also include the influence of wake
panels whose strengths are tied to the trailing-edge body-panel strengths. The
more precise form is

```math
\left(\frac{1}{2}I - K_b - K_w C\right)\phi = \phi_\mathrm{source}.
```

Here `K_b` is the body-panel doublet principal-value influence, `K_w` is the
wake-panel doublet influence at the body control points, and `C` maps the body
doublet unknowns to wake-panel strengths through the trailing-edge/Kutta
coupling. In the implementation, `_G!` assembles `K_b + K_w C` together, with
the limiting body self term included in the diagonal.

## Exterior vs Interior Collocation

The off-diagonal entries of `G` are ordinary panel-to-control-point influence
terms, including any wake-panel influence coupled to a trailing-edge column.
Those entries change only by the small geometric displacement between an
exterior and interior control point. The important difference is the body-panel
diagonal self term, which represents the limiting solid angle of the doublet
kernel.

### Exterior control points

For a smooth flat body panel, the doublet kernel `∂g/∂n` with
`g = −1/(4πr)` integrates to `−½` at the panel's own exterior control point.
Thus the body self part of `G[i,i]` is `−½` on the exterior side, so `−G`
already contains the explicit `+½I` jump term:

```math
(-G)[i,i] = +\tfrac{1}{2} = \tfrac{1}{2}I[i,i] - 0 = ({\tfrac{1}{2}I - G_\text{zeroed}})[i,i].
```

If there are no coupled wake panels, exterior `−G` is equivalent to zeroing the
diagonal and adding `½I` explicitly.

```math
-G_\mathrm{exterior} \approx \tfrac{1}{2}I - G_{\mathrm{exterior},\mathrm{offdiag}}.
```

This is why the exterior test can use the assembled matrix directly:

```math
\phi = (-G_\mathrm{exterior})^{-1}\phi_\mathrm{source}.
```

When wake panels are present, the diagonal of the assembled total matrix may not
be exactly `−½`, because a trailing-edge column can also include same-column wake
influence. That extra wake contribution is not part of the jump term and should
not be discarded. This is why the exterior test uses the assembled matrix
directly:

```math
A_\mathrm{exterior} = -G_\mathrm{exterior}.
```

### Interior control points

On the interior side, `_G!` sees the opposite limiting solid angle and the
body-panel self term is approximately `+½`:

```math
G_\mathrm{interior}[i,i] \approx +\tfrac{1}{2}.
```

If we used raw interior `−G`, the diagonal of the Green matrix would become
`−½`, which has the wrong sign for the identity being checked:

```math
(-G_\mathrm{interior})[i,i] \approx -\tfrac{1}{2}.
```

If there were no wake coupling, the fix would be to replace the diagonal with
`−½` before forming `−G`. With wake coupling, that is too aggressive: it also
removes any same-column wake influence already assembled into `G`.

The correct interior transform is to subtract only the body jump difference
between the two limiting sides. Since the body self term changes from `+½`
inside to the desired exterior-equivalent `−½`, subtract one from the diagonal:

```math
G_\mathrm{interior}[i,i] \leftarrow G_\mathrm{interior}[i,i] - 1,
\qquad
A_\mathrm{interior} = -\left(G_\mathrm{interior} - I\right)
= I - G_\mathrm{interior}.
```

Then

```math
A_\mathrm{interior}
= \frac{1}{2}I - K_{b,\mathrm{offdiag}} - K_w C,
```

up to the small displacement of the interior control points. This keeps the
required `+½I` body jump while preserving wake-panel influence on trailing-edge
columns.

Here `phi_source` is not just `U_wake dot n`; it is the scalar potential induced
by source panels whose strengths are set from the wake-induced normal velocity:

```math
\sigma_i = -\mathbf{u}_{wake,i}\cdot\mathbf{n}_i.
```

This is exactly what `pnl.set_strengths!(body)` does for a Dirichlet
source-doublet body after `body.velocity` has been set to the wake-induced
velocity.

## Test Procedure

The test verifies the identity as follows:

1. Generate a `RigidWakeBody` and a finite `PanelWake`.
2. Directly evaluate the wake-induced potential at exterior body control points.
3. Directly evaluate the wake-induced velocity at those same points.
4. Assemble the doublet matrix `G` at the chosen control-point side.
5. Convert wake-induced velocity to source strengths with `set_strengths!`.
6. Evaluate the source-induced potential `phi_source` from those strengths.
7. Solve `−G φ = phi_source` and compare against the direct wake potential.

The test performs this twice:

- Exterior control points use the unmodified assembled diagonal, since
  `G[i,i] ≈ −½`.
- Interior control points replace the assembled diagonal with `−½` before
  forming `−G`, since raw interior assembly gives `G[i,i] ≈ +½`.

Because normal-derivative data determines potential only up to an additive
constant, the test compares both potentials after subtracting their means.
