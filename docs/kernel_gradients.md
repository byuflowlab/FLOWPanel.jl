# Kernel Gradient Derivations

FLOWPanel stores velocity gradients as

```math
G_{ij} = \frac{\partial u_i}{\partial x_j},
```

where `x` is the target position. The FMM element kernels return potential
`phi`, velocity `u = grad(phi)`, and velocity gradient `G = grad(u)`.

## Constant Source and Constant Doublet Panels

The constant source and doublet panel kernels are evaluated in a panel-local
frame, accumulated edge by edge, and rotated back to the inertial frame. If
`R` maps panel-local vectors to inertial vectors and `G'` is the panel-local
velocity gradient, the returned gradient is

```math
G = R G' R^T.
```

The local formulas in `compute_source_dipole` are the analytic derivatives of
the Hess-Smith edge integrals used by the velocity branches. For the source
kernel, the gradient is the Hessian of the source potential. For the doublet
kernel, the gradient is the derivative of the equivalent vortex-ring edge
velocity. Combined kernels add the source and doublet contributions linearly.

## Finite Vortex Segment

For a finite vortex segment with endpoints `xa` and `xb`, define target-relative
vectors using the implementation convention

```math
r_1 = x_a - x,\qquad r_2 = x_b - x,\qquad s = r_1 - r_2.
```

Then

```math
c = r_1 \times r_2,\qquad
A = c \cdot c,\qquad
B = s \cdot s,\qquad
q = s \cdot \left(\frac{r_1}{|r_1|} - \frac{r_2}{|r_2|}\right).
```

The segment velocity is

```math
u = \frac{1}{4\pi}\frac{c q}{D},
```

with

```math
D = A
```

for the singular kernel and

```math
D = \sqrt{A^2 + r_c^4 B^2}
```

for the Vatistas `n = 2` finite-core kernel.

Because both `r1` and `r2` vary as `dr_i = -dx`, the segment vector `s` is
constant with respect to the target. The cross-product derivative is

```math
dc = -s \times dx.
```

Using

```math
P_i = I - \hat{r}_i \hat{r}_i^T,
```

the scalar derivative is

```math
dq =
\left(
-\frac{P_1 s}{|r_1|}
+\frac{P_2 s}{|r_2|}
\right)\cdot dx.
```

Also,

```math
dA = 2(s \times c)\cdot dx.
```

For the finite-core denominator,

```math
dD = \frac{A}{D} dA,
```

and for the singular denominator, `dD = dA`. Therefore

```math
du =
\frac{1}{4\pi}
\left[
\frac{q}{D}dc
+ c\left(\frac{dq}{D} - \frac{q\,dD}{D^2}\right)
\right].
```

This is assembled directly as the `3 x 3` Jacobian in
`_bound_vortex_gradient`.

## Semi-Infinite Vortex Ray

For a semi-infinite vortex beginning at `p` and extending in unit direction
`d`, define

```math
y = x - p,\qquad
a = y\cdot d,\qquad
h = y - a d,\qquad
H = h\cdot h,\qquad
n = d\times h.
```

With the same Vatistas `n = 2` core radius `r_c`,

```math
E = \sqrt{H^2 + r_c^4}.
```

The semi-infinite velocity branch can be written as

```math
u = -\frac{\Gamma}{4\pi} f n,\qquad
f = \frac{1 + \chi a / |y|}{E},
```

where `chi` is one when the finite bound section from `p` to the projection of
`x` onto the ray is active under the implementation cutoffs, and zero
otherwise.

The derivatives are

```math
dn = d\times dx,
```

```math
d(E^{-1}) = -\frac{2H}{E^3}h\cdot dx,
```

and, when `chi = 1`,

```math
d(a/|y|)
=
\left(\frac{d}{|y|} - \frac{a y}{|y|^3}\right)\cdot dx.
```

Thus

```math
du =
-\frac{\Gamma}{4\pi}
\left[
f\,dn + n\,df
\right].
```

The semi-infinite doublet wake adds the two ray gradients and the finite bound
segment gradient using the same signs as the velocity implementation.
