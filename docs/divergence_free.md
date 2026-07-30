# Regularized Vortex Kernel Divergence

This note checks whether two regularized vortex kernels induce divergence-free
velocity fields on the support of the kernel:

- Gaussian-regularized vortex particles
- Vatistas `n = 2` finite-core vortex segments used by vortex-ring panels

Let a vortex particle at `x_p` carry a constant vector strength `Gamma`. At a
target point `x`, define

```math
r = x - x_p,\qquad \rho = |r|,\qquad \eta = \rho/\sigma .
```

A radially regularized vortex-particle Biot-Savart kernel can be written as

```math
u(r) = \Gamma \times r\, F(\rho),
```

where `F` is a scalar radial factor. For the singular Biot-Savart law,

```math
F(\rho) = \frac{1}{4\pi \rho^3}.
```

For a Gaussian vortex blob, `F` replaces the singular factor by a smooth radial
function. One common form is

```math
F(\rho)
=
\frac{1}{4\pi\rho^3}
\left[
    \operatorname{erf}\left(\frac{\eta}{\sqrt{2}}\right)
    - \sqrt{\frac{2}{\pi}}\eta\exp\left(-\frac{\eta^2}{2}\right)
\right].
```

The exact normalization convention is not important for the divergence check;
the only needed property is that `F` is a scalar function of `rho`.

## Direct Divergence Check

Using index notation with the Levi-Civita symbol,

```math
u_i = \epsilon_{ijk}\Gamma_j r_k F(\rho).
```

Then

```math
\frac{\partial u_i}{\partial x_i}
=
\epsilon_{ijk}\Gamma_j
\frac{\partial}{\partial x_i}\left(r_k F(\rho)\right).
```

Since `r = x - x_p`,

```math
\frac{\partial}{\partial x_i}\left(r_k F(\rho)\right)
= \delta_{ik}F(\rho) +
  r_k F'(\rho)\frac{r_i}{\rho}.
```

---

Since x_p is fixed and r = x - x_p, each component is

```math
r_k = x_k - x_{p,k}.
```

So differentiating r_k F(\rho) with respect to x_i uses the product rule:

```math
\frac{\partial}{\partial x_i}\left(r_k F(\rho)\right)
=
\frac{\partial r_k}{\partial x_i}F(\rho)
+
r_k\frac{\partial F(\rho)}{\partial x_i}.
```

The first derivative is

```math
\frac{\partial r_k}{\partial x_i}
=
\frac{\partial}{\partial x_i}(x_k - x_{p,k})
=
\delta_{ik},
```

because changing x_i changes r_k only when i = k.

For the second derivative, use the chain rule. Since F depends on x_i only through

```math
\rho = |r| = \sqrt{r_1^2 + r_2^2 + r_3^2},
```

we have

```math
\frac{\partial F(\rho)}{\partial x_i}
=
F'(\rho)\frac{\partial \rho}{\partial x_i}.
```

Now

```math
\frac{\partial \rho}{\partial x_i}
=
\frac{\partial}{\partial x_i}(r_j r_j)^{1/2}
=
\frac{1}{2}(r_j r_j)^{-1/2}2r_j\frac{\partial r_j}{\partial x_i}
=
\frac{r_i}{\rho}.
```

Therefore

```math
\frac{\partial F(\rho)}{\partial x_i}
=
F'(\rho)\frac{r_i}{\rho}.
```

Substitute both pieces back into the product rule:

```math
\frac{\partial}{\partial x_i}\left(r_k F(\rho)\right)
=
\delta_{ik}F(\rho)
+
r_k F'(\rho)\frac{r_i}{\rho}.
```

So the two terms are: one from differentiating r_k, and one from differentiating the radial
scalar factor F(\rho).

---



Therefore

```math
\nabla\cdot u
=
\epsilon_{ijk}\Gamma_j\delta_{ik}F(\rho)
+
\epsilon_{ijk}\Gamma_j r_k r_i \frac{F'(\rho)}{\rho}.
```

The first term vanishes because `epsilon_{ijk} delta_{ik} = 0`. The second term
also vanishes because `r_i r_k` is symmetric in `i,k`, while `epsilon_{ijk}` is
antisymmetric in those same indices. Thus

```math
\boxed{\nabla\cdot u = 0}
```

wherever `F` is differentiable.

## Behavior at the Particle Center

The singular kernel is undefined at `rho = 0`, but the Gaussian-regularized
kernel has a finite limit. Expanding the Gaussian form about the origin gives

```math
F(\rho)
=
\frac{1}{6\sqrt{2}\pi^{3/2}\sigma^3}
+ O(\rho^2),
```

so

```math
u(r)
=
\frac{1}{6\sqrt{2}\pi^{3/2}\sigma^3}\Gamma\times r
+ O(\rho^3).
```

This is a smooth local solid-body-rotation field, and its divergence is also
zero at the particle center.

## Support and Cutoff Notes

The mathematical Gaussian has infinite support, so the result applies over all
of `R^3`. If an implementation truncates interactions at a finite radius by
multiplying the velocity by a radial cutoff `c(rho)`, the kernel remains
pointwise divergence free anywhere the cutoff is differentiable because
`c(rho)F(rho)` is still only a radial scalar factor.

For a sharp radial cutoff, the distributional boundary term is also zero for a
spherical cutoff centered on the particle: the velocity is tangent to each
cutoff sphere,

```math
u\cdot\hat{r} = (\Gamma\times r)\cdot\hat{r}\,F(\rho) = 0,
```

so no normal flux is created at the cutoff surface. Non-radial clipping or
interpolation, however, should be checked separately because it can introduce a
nonzero boundary flux.

## Vatistas Finite-Core Vortex-Ring Panels

FLOWPanel's vortex-ring panels are assembled from finite vortex segments. For a
segment with endpoints `x_a` and `x_b`, the implementation uses target-relative
vectors

```math
r_1 = x_a - x,\qquad r_2 = x_b - x,\qquad s = r_1 - r_2 = x_a - x_b.
```

The segment vector `s` is constant with respect to the target. Define

```math
c = r_1\times r_2,\qquad
A = c\cdot c,\qquad
B = s\cdot s,
```

and

```math
q = s\cdot\left(\frac{r_1}{|r_1|} - \frac{r_2}{|r_2|}\right).
```

The finite-core segment velocity is

```math
u = \frac{1}{4\pi}\frac{c q}{D},
```

with the Vatistas `n = 2` denominator

```math
D = \sqrt{A^2 + r_c^4 B^2}.
```

The singular segment uses `D = A`; the finite-core model replaces the singular
`1/A` factor by the smooth `1/D` factor.

### Differential Form

Because `dr_1 = dr_2 = -dx`, the segment vector `s` is fixed and

```math
dc = -s\times dx.
```

The differential of the velocity can be written as

```math
du =
\frac{1}{4\pi}
\left[
    \frac{q}{D}dc
    + c\left(\frac{dq}{D} - \frac{q\,dD}{D^2}\right)
\right].
```

Taking the divergence means taking the trace of the Jacobian mapping `dx` to
`du`. The first term contributes zero trace because the map
`dx -> -s \times dx` is skew-symmetric:

```math
\operatorname{tr}\left(dx \mapsto -s\times dx\right) = 0.
```

For the remaining term, write

```math
dq = a_q\cdot dx,\qquad dD = a_D\cdot dx.
```

Then the second part of the Jacobian is an outer product,

```math
c\left(\frac{a_q}{D} - \frac{q a_D}{D^2}\right)^T,
```

whose trace is

```math
c\cdot\left(\frac{a_q}{D} - \frac{q a_D}{D^2}\right).
```

It remains to show that `c` is orthogonal to both `a_q` and `a_D`.

For `q`, let

```math
P_i = I - \hat{r}_i\hat{r}_i^T.
```

Then

```math
a_q =
-\frac{P_1s}{|r_1|}
+\frac{P_2s}{|r_2|}.
```

Each vector `P_i s` is a linear combination of `s` and `r_i`, hence lies in the
plane spanned by `r_1` and `r_2`. Since `c = r_1 \times r_2` is normal to that
plane,

```math
c\cdot a_q = 0.
```

For the finite-core denominator,

```math
dD =
\frac{A}{D}\,dA,\qquad
dA = 2(s\times c)\cdot dx,
```

so

```math
a_D = \frac{2A}{D}(s\times c).
```

This is also orthogonal to `c`, so

```math
c\cdot a_D = 0.
```

Therefore both trace contributions vanish, and

```math
\boxed{\nabla\cdot u = 0}
```

for every target where the segment formula is differentiable. A vortex-ring
panel is a sum of such segment velocities, so its finite-core induced velocity
is also divergence free away from the segment endpoints and any explicit
implementation cutoffs.

## Conclusion

Both regularizations considered here preserve solenoidal induced velocity
fields on their differentiable support. The Gaussian particle kernel does this
by preserving the `Gamma x r` tangential structure with a radial scalar factor.
The Vatistas finite-core vortex segment does this because its differentiated
terms are either skew-symmetric or orthogonal to the segment-plane normal
`r_1 \times r_2`; summing the segments into a vortex ring preserves zero
divergence.
