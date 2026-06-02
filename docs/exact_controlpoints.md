# Exact Control Points

**Status:** v2 — derivation notes for the kernel self limits.

This note proposes replacing the `CPoffset` regularization currently used to avoid singular self-influence at panel control points with an *exact on-surface* evaluation scheme: evaluate kernels at the true control point on the panel, and substitute a closed-form limit value when the field point is within `ε` of a source-panel control point.

## 1. Motivation and Scope

Panel methods regularize the self-influence singularity in different ways.
FLOWPanel currently displaces every control point off the surface by
`CPoffset * characteristiclength(panel)` along the panel normal, choosing the
sign of the displacement from the body's boundary-condition type parameter
`DBC` (exterior offset for Neumann, interior for Dirichlet). The offset
succeeds in keeping kernels finite, but it couples three things that are
conceptually independent: the geometric definition of *where the control point
lives*, the analytical limit of the influence operator *on the surface*, and
the numerical regularization used by the FMM near field. This note proposes
decoupling them by placing control points exactly on the surface and
substituting closed-form self-limits (§3) wherever a kernel is evaluated
within `ε` of its own control point.

Costs of the current `CPoffset` regularization:

- Formulation-dependent sign flip (positive offset for Neumann, negative for Dirichlet — see `src/FLOWPanel_solver.jl:149`).
- Sensitivity to the offset value (`1e-14` for `NonLiftingBody`, `1e-6` for `RigidWakeBody`).
- FMM evaluation requires a separate `kerneloffset` regularization in `src/FLOWPanel_elements_fmm.jl`.

What this proposal means:

- "Exact control points" — control points sit exactly on the panel surface (centroid, no offset); when a field point is within `ε` of a source-panel control point, the kernel returns its known closed-form limit instead of being numerically evaluated.
- Goal — decouple control-point geometry from the solver formulation and make the diagonal contribution to the influence matrix an explicit per-kernel statement. FLOWPanel uses a fixed self-pair convention: velocity returns the exterior surface limit, while potential returns the interior surface limit.

## 2. Notation and Conventions

This section is a self-contained glossary for the rest of the document.

- **Panel index `j`, panel surface \(S_j\)** — a flat triangular panel with
  vertices \((\mathbf{p}_a, \mathbf{p}_b, \mathbf{p}_c)\), area \(A_j\),
  centroid \(\mathbf{x}_j=(\mathbf{p}_a+\mathbf{p}_b+\mathbf{p}_c)/3\), and
  outward unit normal \(\mathbf{n}_j\) (GeometricTools right-hand-rule
  convention).
- **Characteristic length \(\ell_j\)** — by default
  `characteristiclength_sqrtarea(panel) =` \(\sqrt{A_j}\), defined in
  `src/FLOWPanel_abstractbody.jl`. Used to non-dimensionalize the self-panel
  tolerance \(\varepsilon\).
- **Control point \(\mathbf{x}_{cp,j}\)** — in the current build,
  \(\mathbf{x}_{cp,j}=\mathbf{x}_j+(\mathrm{CPoffset}\cdot\ell_j)\,\mathbf{n}_j\)
  with the sign of `CPoffset` chosen by formulation
  (`src/FLOWPanel_solver.jl:149`). Under the proposed scheme,
  \(\mathbf{x}_{cp,j}\equiv\mathbf{x}_j\) — no offset.
- **Panel frame \((\mathbf{t}_j,\mathbf{o}_j,\mathbf{n}_j)\)** — orthonormal
  triad on the panel. \(\mathbf{t}_j,\mathbf{o}_j\) are tangent;
  \(\mathbf{n}_j\) is outward.
- **Hess–Smith normal \(\mathbf{n}_{HS}\)** — the panel-local \(+z\) direction
  in Hess–Smith derivations, equal to \(-\mathbf{n}_j\). All element kernels
  in `src/FLOWPanel_elements.jl` already absorb this sign; this document uses
  \(\mathbf{n}_j\) throughout for results and only references
  \(\mathbf{n}_{HS}\) when comparing to textbook formulas.
- **Side convention \(\pm\)** — \(+\) denotes the exterior limit
  \(\mathbf{x}\to\mathbf{x}_j+\eta\mathbf{n}_j,\ \eta\downarrow0\); \(-\)
  denotes the interior limit. Jumps are always written as
  \((\cdot)^+-(\cdot)^-\).
- **FLOWPanel self-pair convention** — exact self-pair velocity always uses
  the exterior limit, because that is the physical aerodynamic surface
  velocity. Exact self-pair potential always uses the interior limit, because
  that is the useful boundary-integral value for Dirichlet-style solves. For
  a self panel, exterior doublet or vortex-ring potential is recovered as
  `interior_potential - doublet_or_vortex_strength`.
- **Kernel strengths** — \(\sigma_j\) (source), \(\mu_j\) (doublet),
  \(\Gamma_j\) (vortex-ring circulation),
  \(\boldsymbol{\gamma}_j=\gamma_t\mathbf{t}_j+\gamma_o\mathbf{o}_j\)
  (constant vortex sheet), \(\boldsymbol{\gamma}(\mathbf{y})\) (uniformly
  varying vortex sheet).
- **Principal value (`p.v.`) and finite part (`f.p.`)** — `p.v.` denotes the
  symmetric-shrinkage limit of a Cauchy-singular surface integral; `f.p.`
  denotes the Hadamard finite part used for hypersingular velocity-gradient
  integrals.
- **Self-panel detection** — criterion
  \(\|\mathbf{x}-\mathbf{x}_{cp,j}\|<\varepsilon\). A field point that
  satisfies this criterion for panel \(j\) is treated as "at" that panel's
  control point and receives the §3 self-limit instead of a numerical
  kernel evaluation. The choice of \(\varepsilon\) is discussed in §8.
- **Tie-breaking** — when \(\mathbf{x}\) satisfies the criterion for more than
  one panel (e.g. a shared-edge midpoint), the kernel returns the average of
  the qualifying self-limits weighted by \(1/\|\mathbf{x}-\mathbf{x}_{cp,j}\|^2\)
  (equal weights when distances vanish). See §8 for caveats.
- **`fields` vector** — the per-body bookkeeping array of computed quantities
  is unaffected by this proposal; only the evaluation path for self pairs
  changes.

## 3. Limit Values at the Source Control Point

For each element type defined in `src/FLOWPanel_elements.jl:25–54`, document the on-surface limits at the source panel's own control point of:

- induced potential `φ`,
- induced velocity `u = ∇φ`,
- induced velocity gradient `∇u` (used by `PressureLaplace`, where applicable).

### 3.1 `ConstantSource`

Let a flat source panel \(S_j\) carry constant strength \(\sigma_j\), centroid
\(\mathbf{x}_j\), and GeometricTools outward unit normal \(\mathbf{n}_j\). Define
the scalar single-layer potential with FLOWPanel's sign convention

```math
\phi_\sigma(\mathbf{x})
=
-\frac{\sigma_j}{4\pi}
\int_{S_j}\frac{1}{\|\mathbf{x}-\mathbf{y}\|}\,\mathrm{d}S_y .
```

At an interior point of the panel, including the centroid, this integral is
weakly singular but finite. Therefore the exact self-potential is the ordinary
surface integral

```math
\boxed{
\phi_{\sigma,\mathrm{self}}
=
-\frac{\sigma_j}{4\pi}
\int_{S_j}\frac{1}{\|\mathbf{x}_j-\mathbf{y}\|}\,\mathrm{d}S_y
}
```

or, in the Hess-Smith edge notation already used in
`docs/src/elements/constantsource.md`,

```math
\phi_{\sigma,\mathrm{self}}
=
-\frac{\sigma_j}{4\pi}
\sum_{(a,b)\in\partial S_j}
R_{ab}(\mathbf{x}_j) Q_{ab}(\mathbf{x}_j),
\qquad z=0 .
```

The velocity is the gradient of the single-layer potential. Its tangential
components are ordinary principal-value edge sums at \(z=0\):

```math
\mathbf{u}_{\sigma,\parallel,\mathrm{self}}
=
\operatorname{f.p.}\nabla_\parallel \phi_\sigma(\mathbf{x}_j)
=
\frac{\sigma_j}{4\pi}
\operatorname{f.p.}
\int_{S_j}
\frac{\mathbf{x}_{j,\parallel}-\mathbf{y}_\parallel}
     {\|\mathbf{x}_j-\mathbf{y}\|^3}
\,\mathrm{d}S_y .
```

The normal component is discontinuous across the sheet. With \(\mathbf{n}_j\)
the GeometricTools exterior normal and with FLOWPanel's source sign convention,
the one-sided limits are

```math
\boxed{
\lim_{\eta\downarrow0}
\mathbf{u}_\sigma(\mathbf{x}_j+\eta\mathbf{n}_j)
=
\mathbf{u}_{\sigma,\parallel,\mathrm{self}}
+\frac{\sigma_j}{2}\mathbf{n}_j
}
```

and

```math
\lim_{\eta\downarrow0}
\mathbf{u}_\sigma(\mathbf{x}_j-\eta\mathbf{n}_j)
=
\mathbf{u}_{\sigma,\parallel,\mathrm{self}}
-\frac{\sigma_j}{2}\mathbf{n}_j .
```

Equivalently, in Hess-Smith local coordinates where
\(\hat{\mathbf{n}}_\mathrm{HS}=-\mathbf{n}_j\), the exterior FLOWPanel limit is
the \(z\to -0\) limit and the local normal component is \(-\sigma_j/2\), matching
the note in `docs/src/elements/constantsource.md`.

The velocity gradient has no universal scalar diagonal value analogous to
\(\pm\sigma_j/2\). Away from panel edges the Cauchy finite part is

```math
\boxed{
\nabla\mathbf{u}_{\sigma,\mathrm{self}}^\pm
=
\frac{\sigma_j}{4\pi}
\operatorname{f.p.}
\int_{S_j}
\left(
\frac{\mathbf{I}}{r^3}
-3\frac{\mathbf{r}\mathbf{r}^{T}}{r^5}
\right)\,\mathrm{d}S_y
}
```

with \(\mathbf{r}=\mathbf{x}_j-\mathbf{y}\) and \(r=\|\mathbf{r}\|\). The
\(\pm\) side only affects the distributional derivative of the velocity jump;
the finite part used by an exact control-point implementation is geometry
dependent and must be evaluated by the same edge formulas or automatic
differentiation path as the non-self Hessian, then overridden only where a
pointwise jump term is mathematically defined.

### 3.2 `ConstantDoublet`

Let a flat doublet panel carry constant strength \(\mu_j\). FLOWPanel evaluates
the doublet potential from the source-panel normal velocity identity described
in `docs/src/elements/constantdoublet.md`:

```math
\phi_\mu
=
-\mu_j\frac{\mathbf{n}_\mathrm{HS}\cdot\mathbf{u}_{\sigma=1}}{1}.
```

With the current FMM kernel convention, the source exterior normal velocity is
\(+1/2\) for unit source strength, so the one-sided doublet-potential limits at
the source centroid are

```math
\boxed{
\lim_{\eta\downarrow0}
\phi_\mu(\mathbf{x}_j+\eta\mathbf{n}_j)
=
-\frac{\mu_j}{2}
}
```

and

```math
\lim_{\eta\downarrow0}
\phi_\mu(\mathbf{x}_j-\eta\mathbf{n}_j)
=
+\frac{\mu_j}{2}.
```

The arithmetic average, or principal-value on-surface value, is therefore zero:

```math
\operatorname{p.v.}\phi_{\mu,\mathrm{self}}(\mathbf{x}_j)=0 .
```

This distinction matters for exact control points. A Dirichlet solve that wants
the exterior boundary integral operator must not silently use the principal
value and lose the jump; it must add the \(\frac12 I\) term explicitly if the
self short-circuit returns the principal value.

For constant \(\mu_j\), the ordinary velocity field is equivalent to the
Biot-Savart velocity of a vortex ring on the panel boundary:

```math
\boxed{
\mathbf{u}_{\mu,\mathrm{self}}
=
\frac{\mu_j}{4\pi}
\sum_{(a,b)\in\partial S_j}
\frac{\mathbf{r}_a\times\mathbf{r}_b}
     {\|\mathbf{r}_a\times\mathbf{r}_b\|^2}
\,
\mathbf{r}_{ab}\cdot
\left(
\frac{\mathbf{r}_a}{\|\mathbf{r}_a\|}
-
\frac{\mathbf{r}_b}{\|\mathbf{r}_b\|}
\right)
}
```

where \(\mathbf{r}_a=\mathbf{x}_j-\mathbf{p}_a\),
\(\mathbf{r}_b=\mathbf{x}_j-\mathbf{p}_b\), and
\(\mathbf{r}_{ab}=\mathbf{p}_b-\mathbf{p}_a\). This value is finite for a
centroid strictly inside the panel and not on an edge. There is no additional
\(\pm\frac12\) pointwise velocity jump for a constant doublet sheet; the jump is
in the scalar potential.

The velocity gradient is the derivative of the same boundary-vortex expression:

```math
\boxed{
\nabla\mathbf{u}_{\mu,\mathrm{self}}
=
\mu_j
\sum_{(a,b)\in\partial S_j}
\nabla_{\mathbf{x}}
\mathbf{K}_{ab}(\mathbf{x})
\bigg|_{\mathbf{x}=\mathbf{x}_j}
}
```

where \(\mathbf{K}_{ab}\) is the unit-strength finite-vortex-segment kernel in
the previous equation. This is a geometry-dependent finite value at the
centroid. In implementation terms it is the same analytic segment Jacobian used
for `VortexRing` / doublet-induced gradients, evaluated at the control point.

### 3.3 `VortexRing`

`VortexRing` uses the same finite boundary-vortex velocity as a constant
doublet panel. The scalar-potential path in the FMM implementation is a
constant-doublet surrogate; the physical vortex-ring kernel is a vector-potential
kernel, so the robust self statements are for velocity and velocity gradient.

For ring strength \(\Gamma_j\),

```math
\boxed{
\mathbf{u}_{\Gamma,\mathrm{self}}
=
\frac{\Gamma_j}{4\pi}
\sum_{(a,b)\in\partial S_j}
\frac{\mathbf{r}_a\times\mathbf{r}_b}
     {\|\mathbf{r}_a\times\mathbf{r}_b\|^2}
\,
\mathbf{r}_{ab}\cdot
\left(
\frac{\mathbf{r}_a}{\|\mathbf{r}_a\|}
-
\frac{\mathbf{r}_b}{\|\mathbf{r}_b\|}
\right)
}
```

with the same definitions as in §3.2. At the centroid this is finite for a
non-degenerate panel. For a planar closed ring it points normal to the panel
only when the panel geometry has the corresponding symmetry; in general its
tangential components need not vanish.

The velocity gradient is

```math
\boxed{
\nabla\mathbf{u}_{\Gamma,\mathrm{self}}
=
\Gamma_j
\sum_{(a,b)\in\partial S_j}
\nabla_{\mathbf{x}}
\mathbf{K}_{ab}(\mathbf{x})
\bigg|_{\mathbf{x}=\mathbf{x}_j}.
}
```

The scalar potential used when `scalar_potential=true` follows the
`ConstantDoublet` limits in §3.2 with \(\mu_j=\Gamma_j\):

```math
\phi_{\Gamma,\mathrm{self}}^+=-\frac{\Gamma_j}{2},
\qquad
\phi_{\Gamma,\mathrm{self}}^-=+\frac{\Gamma_j}{2},
\qquad
\operatorname{p.v.}\phi_{\Gamma,\mathrm{self}}=0 .
```

Thus an exact-control-point implementation should not set vortex-ring
self-velocity to zero merely because the target is the source centroid. Zero is
only the limiting value for a target on the filament itself or a deliberate
finite-core regularization choice.

### 3.4 `ConstantVortexSheet`

Let the constant sheet strength be the panel-tangent vector
\(\boldsymbol{\gamma}_j=\gamma_t\mathbf{t}_j+\gamma_o\mathbf{o}_j\), where
\((\mathbf{t}_j,\mathbf{o}_j,\mathbf{n}_j)\) is the GeometricTools panel frame.
There is no scalar velocity potential for a general vortex sheet, so
`φ_self` is not part of this kernel's physical contract:

```math
\boxed{\phi_{\gamma,\mathrm{self}}\ \text{is not defined for a general vortex sheet}.}
```

The Biot-Savart velocity has a tangential jump across the sheet. Denote the
Cauchy principal-value velocity induced by the finite panel at its centroid as

```math
\mathbf{u}_{\gamma,\mathrm{PV}}(\mathbf{x}_j)
=
\operatorname{p.v.}
\frac{1}{4\pi}
\int_{S_j}
\boldsymbol{\gamma}_j
\times
\frac{\mathbf{x}_j-\mathbf{y}}
     {\|\mathbf{x}_j-\mathbf{y}\|^3}
\,\mathrm{d}S_y .
```

The one-sided limits are

```math
\boxed{
\lim_{\eta\downarrow0}
\mathbf{u}_\gamma(\mathbf{x}_j+\eta\mathbf{n}_j)
=
\mathbf{u}_{\gamma,\mathrm{PV}}(\mathbf{x}_j)
+\frac12\,\boldsymbol{\gamma}_j\times\mathbf{n}_j
}
```

and

```math
\lim_{\eta\downarrow0}
\mathbf{u}_\gamma(\mathbf{x}_j-\eta\mathbf{n}_j)
=
\mathbf{u}_{\gamma,\mathrm{PV}}(\mathbf{x}_j)
-\frac12\,\boldsymbol{\gamma}_j\times\mathbf{n}_j .
```

These satisfy the standard vortex-sheet jump condition
\(\mathbf{n}_j\times(\mathbf{u}^+-\mathbf{u}^-)=\boldsymbol{\gamma}_j\).

The velocity gradient again has no universal diagonal constant. The finite part
is obtained by differentiating the Biot-Savart surface integral:

```math
\boxed{
\nabla\mathbf{u}_{\gamma,\mathrm{self}}^\pm
=
\operatorname{f.p.}
\frac{1}{4\pi}
\int_{S_j}
\nabla_{\mathbf{x}}
\left[
\boldsymbol{\gamma}_j
\times
\frac{\mathbf{x}-\mathbf{y}}
     {\|\mathbf{x}-\mathbf{y}\|^3}
\right]_{\mathbf{x}=\mathbf{x}_j}
\,\mathrm{d}S_y
}
```

with the distributional derivative of the velocity jump excluded from the
pointwise finite part.

### 3.5 `UniformVortexSheet`

For a linearly varying sheet, write the tangent strength field over the panel as
\(\boldsymbol{\gamma}(\mathbf{y})\). At the centroid let
\(\boldsymbol{\gamma}_c=\boldsymbol{\gamma}(\mathbf{x}_j)\). As with
`ConstantVortexSheet`, a general scalar velocity potential is not defined:

```math
\boxed{\phi_{\gamma,\mathrm{self}}\ \text{is not defined for a general vortex sheet}.}
```

The local jump depends only on the strength value at the evaluation point, while
the principal-value part depends on the full variation over the panel:

```math
\mathbf{u}_{\gamma,\mathrm{PV}}(\mathbf{x}_j)
=
\operatorname{p.v.}
\frac{1}{4\pi}
\int_{S_j}
\boldsymbol{\gamma}(\mathbf{y})
\times
\frac{\mathbf{x}_j-\mathbf{y}}
     {\|\mathbf{x}_j-\mathbf{y}\|^3}
\,\mathrm{d}S_y .
```

Therefore

```math
\boxed{
\mathbf{u}_{\gamma,\mathrm{self}}^\pm
=
\mathbf{u}_{\gamma,\mathrm{PV}}(\mathbf{x}_j)
\pm
\frac12\,\boldsymbol{\gamma}_c\times\mathbf{n}_j
}
```

where \(+\) denotes the exterior side
\(\mathbf{x}_j+\eta\mathbf{n}_j\). The finite-part velocity gradient is

```math
\boxed{
\nabla\mathbf{u}_{\gamma,\mathrm{self}}^\pm
=
\operatorname{f.p.}
\frac{1}{4\pi}
\int_{S_j}
\nabla_{\mathbf{x}}
\left[
\boldsymbol{\gamma}(\mathbf{y})
\times
\frac{\mathbf{x}-\mathbf{y}}
     {\|\mathbf{x}-\mathbf{y}\|^3}
\right]_{\mathbf{x}=\mathbf{x}_j}
\,\mathrm{d}S_y
}
```

again excluding the distributional derivative of the tangential velocity jump.
The only universal self term is the
\(\pm\frac12\boldsymbol{\gamma}_c\times\mathbf{n}_j\) jump; all other terms are
finite, geometry-dependent, and must be evaluated by the panel quadrature or
closed-form sheet formula.

## 4. Solver Impact

### 4.1 Neumann Formulation (`BackslashNeumann`)

- Current build: `_solve!(body, ::Backslash)` for `DBC=false`
  (`src/FLOWPanel_solver.jl:701–721`) assembles \(G[i,j]=\mathbf{u}_{ij}\cdot\mathbf{n}_i\)
  at offset control points (`_G!` line 179). The diagonal is whatever
  numerical value `induced` returns at the offset point.
- Under exact evaluation, the off-diagonal entries are unchanged to
  \(\mathcal{O}(\mathrm{CPoffset})\): control points moved by
  \(10^{-14}\,\ell_j\) (`NonLiftingBody`) or \(10^{-6}\,\ell_j\)
  (`RigidWakeBody`) onto the surface shifts well-resolved off-diagonal
  kernels only at that order.
- The diagonal \(G[i,i]\) is replaced by the exterior \(+\) self-velocity
  limit from §3 projected onto \(\mathbf{n}_i\), evaluated at unit strength:
  - `ConstantSource`: \(G[i,i]=\mathbf{u}_{\sigma,\parallel,\mathrm{self}}(S_i)\cdot\mathbf{n}_i+\tfrac12\) (§3.1, the \(+\tfrac12\) is the source-sheet normal jump).
  - `ConstantDoublet` / `VortexRing`: \(G[i,i]=\mathbf{u}_{\mu,\mathrm{self}}(S_i)\cdot\mathbf{n}_i\) (§3.2/§3.3). No jump contribution, but a finite geometry-dependent value — it must be evaluated by the boundary-vortex segment sum at the centroid, not zeroed.
  - `ConstantVortexSheet` / `UniformVortexSheet`:
    \(G[i,i]=(\mathbf{u}_{\gamma,\mathrm{PV}}+\tfrac12\boldsymbol{\gamma}_i\times\mathbf{n}_i)\cdot\mathbf{n}_i\) (§3.4/§3.5). The normal projection of the jump \(\boldsymbol{\gamma}_i\times\mathbf{n}_i\) vanishes, so the diagonal reduces to the normal projection of the principal-value velocity.
- Net effect: the `DBC=false` branch of `_set_formulation_geometry!`
  (`src/FLOWPanel_solver.jl:225-234`) becomes a no-op modulo the
  unconditional control-point/normal recomputation.
- Exact self-pair evaluation no longer needs a formulation-side flag. The
  constant-source self-velocity jump is always the exterior
  \(+\tfrac12\sigma\mathbf{n}\) contribution, while doublet and vortex-ring
  self-potential always use the interior \(+\tfrac12\mu\) or
  \(+\tfrac12\Gamma\) contribution. If an exterior self-potential is needed,
  subtract the self panel's doublet or vortex strength from the returned
  interior value.

### 4.2 Dirichlet Formulation (`BackslashDirichlet`)

- Current build: `_solve!` for `DBC=true` mixed-element bodies
  (`src/FLOWPanel_solver.jl:723–739`); `_G!` evaluates `phi` at the
  interior-offset control point (`src/FLOWPanel_solver.jl:177`). The
  \(\tfrac12 I\) doublet diagonal that appears in textbook Dirichlet
  formulations is currently produced *numerically* by the interior offset:
  the interior potential limit of a constant doublet sheet is \(+\mu_j/2\)
  (§3.2).
- Under exact evaluation, the on-surface principal-value doublet potential
  is zero (§3.2), so the \(\tfrac12 I\) term must be added back explicitly:
  \(G[i,i]\mathrel{+}=\tfrac12\) for the doublet column at unit strength.
  For the mixed element bodies actually supported by this `_solve!`
  (\(\mathrm{Union}\{\text{ConstantSource},\text{ConstantDoublet}\}\) and
  \(\mathrm{Union}\{\text{ConstantSource},\text{VortexRing}\}\), `NK=2`), the
  diagonal splits across columns: the source column contributes the finite
  geometry-dependent self-potential \(\phi_{\sigma,\mathrm{self}}(S_i)\)
  (§3.1), and the doublet / vortex-ring column contributes \(+\tfrac12\).
- The interior-potential source term assembled by `influence!` in the
  `DBC=true` `solve!` wrapper (`src/FLOWPanel_solver.jl:250`) is unchanged:
  it depends on *external* influences populating `body.potential`, not on
  the diagonal.
- Cross-reference `docs/dirichlet_potential_theory.md`. The derivation
  there documents the textbook origin of \(\tfrac12 I\); this scheme moves it
  from an emergent numerical artifact of the interior offset to an
  explicit matrix-assembly step.

### 4.3 Matrix-Free Solvers

- `KrylovSolver` (Krylov.jl GMRES, matrix-free) applies the influence
  operator through `influence!` rather than a stored matrix. The
  self-contribution enters whenever `influence!(body, body, …)` reaches a
  self pair. The short-circuit therefore must live inside
  `induced(target, source_system, i_source; kerneloffset)`
  (`src/FLOWPanel_elements_fmm.jl:152`), not at the matrix-assembly level.
  The Dirichlet \(\tfrac12 I\) diagonal that §4.2 adds explicitly must also
  be applied here, as an operator add-on after the body\(\to\)body apply
  (an in-place \(+\tfrac12\,x_i\) on each doublet/vortex-ring component
  before returning).
- `FGSSolver` wraps `FastMultipole.FastGaussSeidel`
  (`src/FLOWPanel_solver.jl:688`). Its sweep computes a per-panel diagonal
  coefficient from the FMM near field. Under the proposed scheme that
  per-panel diagonal coefficient must be the §3 self-limit (projected onto
  \(\mathbf{n}_i\) for Neumann, plus \(\tfrac12\) for Dirichlet
  doublet/vortex-ring columns), not the value emitted by `induced` at an
  offset control point. The precise hook into `FastGaussSeidel` for
  injecting an analytical diagonal will be confirmed at implementation
  time.

### 4.4 Removal of the Formulation-Dependent `CPoffset` Sign Flip

- With diagonals supplied directly from §3 limits, the same set of
  on-surface control points serves both formulations. The two sign-flip
  statements `target_system.CPoffset = abs(CPoffset_old) * (DBC ? -1 : 1)`
  (`_G!`, `src/FLOWPanel_solver.jl:149`) and
  `body.CPoffset = abs(CPoffset_old) * (DBC ? -1 : 1)`
  (`_set_formulation_geometry!`, line 229) are deleted, together with the
  restoration sites at lines 152, 254, 277, 691 (and the equivalent
  restoration in any future `KrylovSolver` / `FGSSolver` constructor sites
  that still carry it).
- `_set_formulation_geometry!` collapses to a `calc_normals!` /
  `calc_controlpoints!` pair (when `update_cps_normals=true`) with no
  state mutation to restore in `finally`. The `try`/`finally` blocks in
  the two `solve!` wrappers at `src/FLOWPanel_solver.jl:236–281` lose the
  `CPoffset` restoration line and reduce to a plain call.

## 5. FMM Path

- The natural site for the self short-circuit is the entry point in
  `src/FLOWPanel_elements_fmm.jl`: the two-argument
  `induced(target, source_system, i_source, …; kerneloffset)` overloads at
  lines 133 and 152 (plus the `VortexRing` non-rotated overloads starting
  at line 173). After computing the source-panel control point
  (line 158: `control_point = (v1 + v2 + v3) / 3`), compare
  \(\|\mathbf{target}-\mathbf{control\_point}\|\) to \(\varepsilon\); if
  within tolerance, return the §3 self-limit (potential, velocity, and
  optional gradient) instead of dispatching to `_induced`.
- `FastMultipole.direct!` reaches `induced` through the same path, so
  handling the self pair inside `induced` covers both `DirectBackend` and
  the FMM near-field passes uniformly. No change is needed at the
  `FastMultipole.fmm!` interface.
- `kerneloffset` (default \(10^{-3}\), body-level overrides
  `kerneloffset_panel`, `kerneloffset_targets`) remains the regularization
  for *non-self* near-singular evaluations — adjacent panels with very
  small gap, or vortex filaments grazing a probe point. In particular,
  wake objects (`AbstractFreeWake`, `PanelWake`, semi-infinite trailing
  filaments) can drift arbitrarily close to body panels during
  rollup/convection and land in singular regions that the §3 self-limits
  do not cover. `kerneloffset` is **not** retired by this proposal; it
  stays in place exactly as today for every non-self pair.

## 6. Removal of `CPoffset` from the Body Types

- No replacement side-convention field is required for self pairs. Exact
  self-pair velocity and potential use the fixed convention stated in §2.
- Fields to retire:
  - `NonLiftingBody.CPoffset` (`src/FLOWPanel_nonliftingbody.jl:51`, default `1e-14` at line 77).
  - `RigidWakeBody.CPoffset` (`src/FLOWPanel_liftingbody.jl:63`, default `1e-6` at line 99).
  Both are `Float64` fields; deleting them changes the struct layout, so
  all downstream code that reads `body.CPoffset` must be migrated in the
  same patch.
- `calc_controlpoints!` (`src/FLOWPanel_abstractbody.jl:701–721`) loses the
  normal-offset block (current lines 715–718) and reduces to the centroid
  average. The `off` keyword argument becomes unused for this purpose and
  can be removed in v2; the `characteristiclength` keyword is retained
  because it is consumed elsewhere (force normalization, the
  \(\varepsilon\) scaling proposed in §8).
- External callers that pass `CPoffset` to body constructors need a
  deprecation path. Grep the wider workspace for `CPoffset=` and
  `body.CPoffset` to enumerate call sites:
  - In-repo callers: the solver-internal restoration sites already listed
    in §4.4; test helpers and examples that set `CPoffset` for
    reproducibility (to be re-enumerated at implementation time).
- v2 deprecation: keep `CPoffset` as an *accepted-but-ignored* keyword for
  one minor release, emitting `@warn maxlog=1` from each body
  constructor; remove in the next minor release.

## 7. Validation Strategy

Concrete test recipes for the eventual implementation. All paths are
relative to the repository root.

1. **Diagonal closed-form check** — in
   `test/runtests_unit_kernel_gradient.jl` (which already exercises
   direct-vs-analytic kernel derivatives), add a
   `@testset "exact CP self-limits"` that, for each element type, builds a
   single flat triangle, calls `induced(centroid, body, 1, …)` with
   `kerneloffset=0` and the new self short-circuit active, and asserts
   equality with:
   - exterior \(+\tfrac12\,\mathbf{n}\) normal-velocity jump for
     `ConstantSource`, independent of `DBC`,
   - interior \(+\tfrac12\) potential for `ConstantDoublet` and `VortexRing`,
     independent of `DBC`; exterior self-potential is recovered by subtracting
     the self panel's doublet or vortex strength,
   - \(\tfrac12\,\boldsymbol{\gamma}\times\mathbf{n}\) velocity jump for
     `ConstantVortexSheet` and `UniformVortexSheet`.
2. **Sphere potential-flow** — extend the existing sphere test in
   `test/runtests_analytical.jl` so the same mesh is solved with both
   `BackslashNeumann` and `BackslashDirichlet`, with and without the new
   exact-CP path, and the surface \(C_p\) is compared to the analytical
   \(1-(\tfrac32\sin\theta)^2\) at the same tolerance as the existing test.
3. **Dirichlet identity** — `test/dirichlet_potential_test.jl` already
   verifies the interior-potential boundary-integral identity. Add a
   parametrization on `exact_controlpoints = true | false` and assert the
   identity holds in both modes to within the existing tolerance.
4. **Diagonal-vs-offset cross-check** — in `test/runtests_unit_solver.jl`,
   add a test that constructs `G` once via the legacy offset path and once
   via the exact path on the same mesh, and asserts
   \(\|G_\mathrm{exact}-G_\mathrm{offset}\|/\|G_\mathrm{offset}\|<\mathcal{O}(\mathrm{CPoffset})\)
   for refined meshes (sphere and a generic wing). Exercise both
   `DBC=true` and `DBC=false`.
5. **FMM near-field** — in `test/runtests_unit_fmm.jl`, add a single
   triangle plus a probe point at the centroid and verify that `induced`
   with `kerneloffset=1e-3` returns the §3 limit (within `kerneloffset`)
   when the short-circuit is enabled. Also probe an off-centroid point to
   verify no behavior change elsewhere.

## 8. Open Questions

Recommended answers:

- **Choice of \(\varepsilon\).** Recommend
  \(\varepsilon=\varepsilon_\mathrm{rel}\cdot\)`characteristiclength_sqrtarea(panel)`
  with \(\varepsilon_\mathrm{rel}=10^{-12}\). Rationale: \(\varepsilon\) must
  (1) scale with panel size so the criterion is mesh-invariant,
  (2) be tight enough that the normal jump is not double-counted at
  neighbouring control points, and (3) be loose enough to absorb
  floating-point round-off in \(\|\mathbf{x}-\mathbf{x}_{cp,j}\|\). A
  relative tolerance of \(10^{-12}\) sits well below typical inter-panel
  spacings on production meshes (\(\sim10^{-2}\,\ell_j\)) and well above
  `eps(Float64)`.
- **Does `VortexRing` self-velocity need the same explicit treatment as
  `ConstantDoublet`?** Partly. The existing edge-segment evaluation in
  `_induced` already returns the correct *velocity* finite part at the
  centroid (a closed vortex ring has no \(\pm\) velocity jump; cf. §3.3),
  so the velocity short-circuit reduces to "use the current path".
  However, the scalar-potential branch of `_induced` (used when
  `scalar_potential=true`) inherits the `ConstantDoublet` jump and still
  requires the explicit \(\pm\Gamma/2\) self value. Implementation: keep
  the edge-segment velocity unchanged; add an explicit branch for the
  potential.

Left open:

- **Shared-edge tie-breaking.** When a field point lies near the shared
  edge of two panels, both \(\|\mathbf{x}-\mathbf{x}_{cp,1}\|\) and
  \(\|\mathbf{x}-\mathbf{x}_{cp,2}\|\) may fall below \(\varepsilon\) for
  very anisotropic meshes. The inverse-square weighting proposed in §2 is
  a placeholder; a principled choice requires numerical experiments on
  shared-edge probe configurations and is tracked here.
- **Deforming geometry in unsteady simulation.** `simulate!` recomputes
  control points each step via `propagate_kinematics!`. The §3 self-limits
  depend only on the current panel geometry, so they remain valid in
  principle. Open: confirm that the per-step recomputation ordering inside
  `simulate!` (`src/FLOWPanel_simulate.jl`) does not cache a stale
  per-panel \(\varepsilon\) (or stale self-limit) between steps.
