# 025 Phase 1 — Regularized Filament Kernel Derivations

**Date:** 2026-08-20. All formulae are derived from the shipped implementation
(`src/FLOWPanel_elements_fmm.jl`, `_bound_vortex_velocity` /
`_bound_vortex_gradient`), so they bind to the code rather than a textbook
variant. Families: Vatistas n=2 (legacy, as built), compact-support (selectable),
Gaussian (DEFAULT — Ryan's matched-core-size ruling 2026-08-20 superseded the
Phase 0 matched-Δr decision; see the phase_00 addendum).

## Code quantities

For one filament side of a vortex-ring panel with unit circulation, the code
forms (target $\mathbf{x}$; vertex positions $\mathbf{v}_i$, $\mathbf{v}_{i+1}$):

$$
\mathbf{r}_1 = \mathbf{v}_i - \mathbf{x}, \quad
\mathbf{r}_2 = \mathbf{v}_{i+1} - \mathbf{x}, \quad
\mathbf{c} = \mathbf{r}_1 \times \mathbf{r}_2, \quad
\mathbf{s} = \mathbf{r}_1 - \mathbf{r}_2,
$$

$$
A = |\mathbf{c}|^2, \qquad B = |\mathbf{s}|^2, \qquad
q = \mathbf{s} \cdot \left(\hat{\mathbf{r}}_1 - \hat{\mathbf{r}}_2\right).
$$

$\mathbf{s} = \mathbf{v}_i - \mathbf{v}_{i+1}$ is independent of the target, so
$B$ is a constant under target differentiation — the property the existing
gradient code already exploits. The perpendicular distance from the target to
the filament *line* is

$$
h = \frac{|\mathbf{r}_1 \times \mathbf{r}_2|}{|\mathbf{s}|} = \sqrt{A/B}.
$$

## Velocity

All families share the numerator and differ only in the scalar denominator
$D$:

$$
\mathbf{u} = \frac{\mathbf{c}\, q}{4\pi\, D}.
$$

- **Singular:** $D_{\rm sing} = A = B h^2$.
- **Vatistas n=2 (as built):** $D_{\rm Vat} = \sqrt{A^2 + r_c^4 B^2}
  = B\sqrt{h^4 + r_c^4}$ — the code comment's $1/h^2 \to 1/\sqrt{h^4+r_c^4}$.
- **Compact-support (selectable):**

$$
D_{\rm cpt} = A + \delta(h)\, B = B\left(h^2 + \delta(h)\right), \qquad
\delta(h) = \begin{cases} (h - r_c)^2, & h < r_c \\ 0, & h \ge r_c \end{cases}
$$

  — the transplant of `regularize` (`:954`) into the filament transverse
  profile. Exactly singular for $h \ge r_c$. At $h \to 0$:
  $D \to B r_c^2$ and $\mathbf{c} \to 0$ linearly, so $\mathbf{u} \to 0$ on
  the filament line (same endpoint guard as today applies).
- **Gaussian:** $\displaystyle D_{\rm gau} = \frac{A}{g(h)}, \qquad
  g(h) = 1 - e^{-h^2/(2 r_c^2)}$, i.e.
  $\mathbf{u} = \mathbf{c}\,q\,g(h)/(4\pi A)$ (Lamb–Oseen transverse profile;
  $r_c \equiv \sigma$). At $h \to 0$: $g/A \to 1/(2 r_c^2 B)$, finite.

## Velocity gradient

The code computes $\partial \mathbf{u}/\partial \mathbf{x}$ as

$$
\frac{\partial \mathbf{u}}{\partial \mathbf{x}} =
\frac{1}{4\pi}\left[\frac{q}{D}\,\frac{\partial \mathbf{c}}{\partial \mathbf{x}}
+ \mathbf{c}\, (\nabla f)^{T}\right], \qquad f = \frac{q}{D}, \qquad
\nabla f = \frac{\nabla q}{D} - \frac{q\, \nabla D}{D^2},
$$

with the target-derivative identities already in the code
($d\mathbf{r}_i/d\mathbf{x} = -I$, $\mathbf{s}$, $B$ constant):

$$
\frac{\partial \mathbf{c}}{\partial \mathbf{x}} = -[\mathbf{s}]_\times, \qquad
\nabla A = 2\, \mathbf{s} \times \mathbf{c}, \qquad
\nabla q = -\frac{(I - \hat{\mathbf{r}}_1 \hat{\mathbf{r}}_1^T)\,\mathbf{s}}{|\mathbf{r}_1|}
          + \frac{(I - \hat{\mathbf{r}}_2 \hat{\mathbf{r}}_2^T)\,\mathbf{s}}{|\mathbf{r}_2|}.
$$

(Correction 2026-08-20, review finding 2: with $\mathbf{s} = \mathbf{r}_1 -
\mathbf{r}_2$ and $d\mathbf{r}_i/d\mathbf{x} = -I$, the $j$-th column of
$\partial\mathbf{c}/\partial\mathbf{x}$ is $\mathbf{e}_j \times \mathbf{s} =
-[\mathbf{s}]_\times \mathbf{e}_j$, i.e. $\partial\mathbf{c}/\partial
\mathbf{x} = -[\mathbf{s}]_\times$; an earlier revision of this doc dropped
the sign. The IMPLEMENTATION was already correct — its `dc_dx` literal
encodes the signed matrix, and the finite-difference unit tests bind it.)

Every family reduces to $\nabla D = \kappa\, \nabla A$ for a scalar $\kappa$
(using $\nabla h = \nabla A / (2 h B)$, from $h^2 = A/B$ with $B$ constant):

- **Singular:** $D = A$, $\kappa = 1$.
- **Vatistas (as built):** $D = \sqrt{A^2 + r_c^4 B^2}$,
  $\nabla D = A \nabla A / D \Rightarrow \kappa = A/D$ — exactly the code's
  `dD_coeff = (A / D) * (2 * cross(s, c))`.
- **Compact ($h < r_c$):** $\nabla D = \nabla A + B\,\delta'(h)\,\nabla h$
  with $\delta'(h) = 2(h - r_c)$:

$$
\nabla D = \nabla A \left(1 + \frac{\delta'(h)}{2h}\right)
         = \left(2 - \frac{r_c}{h}\right) \nabla A
\quad\Rightarrow\quad \kappa = 2 - \frac{r_c}{h},
$$

  and $\kappa = 1$ for $h \ge r_c$. Continuity: $\kappa \to 1$ as
  $h \to r_c^-$, so the gradient is continuous at the support boundary
  ($C^1$ kernel). As $h \to 0$, $\kappa \sim -r_c/h$ diverges but
  $\nabla A = 2\mathbf{s}\times\mathbf{c} \to 0$ linearly in $h$, so
  $\nabla D$ stays finite; $\mathbf{u}$ and its gradient remain bounded (peak
  $2.55\,\Gamma/(2\pi\Delta r^2)$, Phase 0).
- **Gaussian:** with $D = A/g(h)$ and
  $g'(h) = \dfrac{h}{r_c^2}\, e^{-h^2/(2 r_c^2)}$:

$$
\nabla D = \frac{\nabla A}{g} - \frac{A\, g'(h)\, \nabla h}{g^2}
         = \frac{\nabla A}{g}\left(1 - \frac{h\, g'(h)}{2\, g}\right)
\quad\Rightarrow\quad
\kappa = \frac{1}{g}\left(1 - \frac{h g'(h)}{2 g}\right),
$$

  using $A \nabla h = h^2 B \cdot \nabla A/(2hB) = (h/2)\nabla A$. As
  $h \to 0$: $g \sim h^2/(2r_c^2)$, $h g' / (2g) \to 1$, and
  $\kappa$ stays finite after cancellation ($\kappa \to$ finite limit; the
  implementation evaluates the guarded form and the unit tests verify the
  small-$h$ limit against finite differences).

So the implementation only replaces the pair $(D, \kappa)$ per family;
$\partial\mathbf{c}/\partial\mathbf{x}$, $\nabla q$, and the assembly line are
untouched.

## Source/doublet panel kernels (`regularize` path), as shipped

(Added 2026-08-20, review finding 3a — the gate requires the source/doublet
derivations, not only the filament family.) Both kernels are evaluated per
panel side $i \to i{+}1$ in the panel-local frame (rotation $R$; target at
$(x, y, z)$ with $z$ along the panel normal; side vertices $(x_i, y_i)$,
$(x_{i+1}, y_{i+1})$; $d_x = x_{i+1} - x_i$, $d_y = y_{i+1} - y_i$,
$d_s = \sqrt{d_x^2 + d_y^2}$; $r_i, r_{i+1}$ = target-to-vertex distances),
then summed over sides and scaled by $-1/4\pi$ and rotated back
(`_induced` at src/FLOWPanel_elements_fmm.jl:710; side terms in
`compute_source_dipole`). The shipped Hess–Smith side primitives:

$$
L = \ln\frac{r_i + r_{i+1} - d_s}{r_i + r_{i+1} + d_s}, \qquad
T = \operatorname{atan2}\!\big(d_x z\,(a_1 - a_2),\; z^2 d_x^2 + a_1 a_2\big),
$$

with $a_k = (d_y e_k - h_k d_x)/r_k$ the code's `arg1/arg2` combinations
($e_k, h_k$ per `recurse_source_dipole`), plus guards forcing the principal
value $T = 0$ on the self-pair and on panel-side extensions.

**Constant source** (strength $\sigma$; before the $-1/4\pi$ prefactor):

$$
\phi_\sigma = \sigma\left[\frac{(x - x_i)d_y - (y - y_i)d_x}{d_s}\,L
  + z\,T\right], \qquad
\mathbf{u}_\sigma = \sigma\left(\frac{d_y}{d_s}L,\; -\frac{d_x}{d_s}L,\;
  T\right),
$$

and the gradient entries $\phi_{xx}, \phi_{xy}, \ldots$ as coded (lines
498–508), built from $\rho = r_i r_{i+1} + (x-x_i)(x-x_{i+1}) +
(y-y_i)(y-y_{i+1}) + z^2$ and $\lambda = (x-x_i)(y-y_{i+1}) -
(x-x_{i+1})(y-y_i)$. **The `reg_term` does not enter the source branch**:
its potential, velocity, and gradient are the singular Hess–Smith forms
(with the $L$/$T$ guard limits only).

**Constant doublet** (strength $\mu$): $\phi_\mu = -\mu\, T$ (unregularized),
and the velocity

$$
\mathbf{u}_\mu = -\mu\,\big(z\, d_y,\; -z\, d_x,\; \lambda\big)\,
\frac{r_i + r_{i+1}}{r_i\, r_{i+1}\, \rho + \delta_e},
\qquad
\delta_e = \texttt{regularize}(d_e, r_c) =
\begin{cases} (d_e - r_c)^2, & d_e < r_c\\ 0, & d_e \ge r_c\end{cases}
$$

where $d_e$ = `minimum_distance` from the target to the side segment. **This
denominator shift is the ONLY place `regularize` enters the source/doublet
path**: the doublet potential and the doublet gradient (`psi_**`, lines
571–586) are evaluated unregularized — the shipped gradient code carries
`# + reg_term` comments showing the regularization deliberately disabled
there. Beyond $d_e \ge r_c$, $\delta_e = 0$ and every quantity is exactly
the singular kernel, which is what justifies the tol-independent
`radius_inflation = kerneloffset` rule for these kernels.

## Scalar potential of the VortexRing element — KNOWN LIMITATION

(Rewritten 2026-08-20, review finding 3b; the previous text claimed a
velocity/potential family unification that held only under the superseded
compact default.) The VortexRing potential branch evaluates the equivalent
constant-doublet panel (`_induced(..., ConstantDoublet, core_size, ...)`),
whose only regularization is the compact-support $\delta_e$ shift above.
Under the GAUSSIAN default the VortexRing **velocity** is
Gaussian-regularized while its **scalar potential** remains
compact-regularized: inside the core region ($h \lesssim r_c$ of an edge)
$\nabla\phi \ne \mathbf{u}$ for this element. Outside both cores the two
agree with the singular kernel to their respective contracts, so the
mismatch is confined to the regularized neighborhood. Unifying them would
require a Gaussian-regularized doublet-potential derivation — an open Ryan
decision, recorded in the item's Deferred section.

## Radius inflation per family

`radius_inflation(VortexRing, r_c, tol)` = distance beyond which the
regularized kernel matches the singular kernel within relative tolerance
`tol` — for BOTH the velocity and the gradient (revised 2026-08-20, review
finding 1: the rule must be gradient-aware).

For the transverse profile $u(h)$, write $z = h^2/(2 r_c^2)$.

- **Gaussian (default, gradient-aware):** velocity rel. error is
  $1 - g = e^{-z}$; the radial-derivative rel. error follows from

$$
\frac{du_{\rm gau}}{dh} = \frac{d}{dh}\left[\frac{1 - e^{-z}}{h}\right]
= -\frac{1}{h^2}\left[1 - e^{-z}(1 + 2z)\right],
$$

  so the gradient rel. error vs $-1/h^2$ is $e^{-z}(1 + 2z) > e^{-z}$. The
  velocity-derived radius $r_c\sqrt{2\ln(1/{\rm tol})}$ would leave gradient
  error ${\rm tol}\,(1 + 2\ln(1/{\rm tol}))$ — 28.6 tol at $10^{-6}$. The
  shipped rule solves

$$ e^{-z_*}(1 + 2 z_*) = {\rm tol}, \qquad
   \Delta r = r_c \sqrt{2 z_*}, $$

  via the fixed point $z \leftarrow \ln((1+2z)/{\rm tol})$ from
  $z_0 = \ln(1/{\rm tol})$ (contraction rate $2/(1+2z) \approx 0.06$; 5
  iterations). Values: $\Delta r/r_c = 4.99,\ 5.47,\ 5.90$ at
  tol $= 10^{-4}, 10^{-5}, 10^{-6}$.

- **Compact:** exactly singular beyond the support — velocity AND gradient —
  so

$$ \Delta r = r_c \quad (\text{tol-independent}). $$

- **Vatistas (as built, legacy-pinned):** velocity rel. error
  $\approx \tfrac12 (r_c/h)^4 \le$ tol at

$$ \Delta r = r_c\, (2/{\rm tol})^{1/4}. $$

  This shipped rule is deliberately UNCHANGED (legacy-exact reproduction is
  pinned by tests). Honestly stated: from
  $u_{\rm Vat} \approx (1/h)(1 - \tfrac12 x)$ with $x = (r_c/h)^4$,
  $du/dh = -(1/h^2)(1 - \tfrac52 x)$, so the gradient rel. error coefficient
  is $\tfrac52 (r_c/h)^4$ and at the shipped radius ($x = {\rm tol}/2$) the
  gradient error is $\le \tfrac52 \cdot \tfrac{\rm tol}{2} = 1.25\,{\rm
  tol}$, absorbed by the multipole-acceptance clearance margin — same
  argument as the original docstring note.
