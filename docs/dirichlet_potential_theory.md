# Dirichlet potential identity at exact panel centroids

## Scope

This note states the Green identity used by the velocity-through-sources wake
formulation and maps it to FLOWPanel's current Dirichlet operator. The current code
collocates exactly at panel centroids and inserts a selected one-sided self limit. There
is no exterior/interior control-point displacement or runtime `cp_outer` switch.

The companion note, `docs/wake_solve_schemes.md`, applies the identity to the finite
attached-wake Kutta coupling and to comparisons with a semi-infinite wake.

## Symbols and signs

Let \(\Omega\) be the closed body interior, \(\partial\Omega\) its consistently
oriented boundary, and \(\boldsymbol n\) the outward unit normal. Let \(\psi\) be a
single-valued harmonic potential in \(\Omega\), with

\[
q=\left.\psi\right|_{\partial\Omega},
\qquad
\sigma=-\frac{\partial\psi}{\partial n}.
\]

The operators are:

- \(\mathcal S\): the body single-layer potential or its interior trace;
- \(\mathcal D\): the body double-layer potential away from the boundary;
- \(\mathcal B\): the complete interior trace of \(\mathcal D\), including the
  positive half jump used by FLOWPanel;
- \(\mathcal W\): attached-wake potential influence at body control points;
- \(\mathcal C\): Kutta map from body doublet strengths to attached-wake strengths;
- \(\boldsymbol\mu\): body doublet-strength vector.

Thus the body-only boundary operator is \(\mathcal B\), while the lifting-body solve
operator is

\[
\mathcal G=\mathcal B+\mathcal W\mathcal C.
\]

The attached wake is not part of the closed-boundary Green identity. It is included in
\(\mathcal G\) only because its strength is coupled to the body unknowns.

## Continuous derivation

Suppose the singular support that generates \(\psi\) is disjoint from the closed body
interior and boundary. Green's representation in the interior is

\[
\psi(\boldsymbol x)
=\mathcal S\sigma(\boldsymbol x)+\mathcal Dq(\boldsymbol x),
\qquad
\boldsymbol x\in\Omega.
\]

Take the limit from the body interior. With \(\mathcal B\) defined to include the
chosen interior jump,

\[
q=\mathcal S\sigma+\mathcal Bq.
\]

Rearranging gives the boundary identity

\[
\boxed{\mathcal S\sigma=(I-\mathcal B)q.}
\]

The two body layers reconstruct \(\psi\) in the interior. In the exterior component
containing infinity, they reconstruct zero because the singular support responsible
for \(\psi\) lies outside the Green surface:

\[
\mathcal S_{\mathrm e}\sigma+\mathcal D_{\mathrm e}q=0.
\]

This does not say that the single layer is the negative of the wake potential. The
double layer carrying the boundary trace is indispensable.

The exact statement requires:

- a closed boundary with consistent outward orientation;
- singular support disjoint from the body interior and boundary;
- a single-valued harmonic branch of \(\psi\) throughout the interior, with a
  well-defined trace and normal derivative;
- matching kernel signs, one-sided limits, and gauge conventions.

If a vortex-sheet potential needs a branch cut, that cut must remain outside the body
interior and boundary.

## FLOWPanel's exact-centroid operator

`calc_controlpoints!` in `src/FLOWPanel_abstractbody.jl` places every triangular-panel
control point at its centroid. The normals argument retained by that function is for
backwards compatibility and does not offset the point.

At a target equal to its source-panel centroid, the raw doublet kernel supplies its
principal value. `_self_limit` in `src/FLOWPanel_elements_fmm.jl` replaces the
self-potential with the selected interior limit. For a body doublet strength
\(\mu_i\), the self contribution is

\[
\frac{1}{2}\mu_i.
\]

`_G!` in `src/FLOWPanel_solver.jl` consumes that side-aware value directly. It neither
moves the target nor adds another diagonal jump. In operator form, its body-only
Dirichlet block approximates \(\mathcal B\). For a lifting body, `_induced_wake` adds
the same-column attached-wake influence, so the assembled block approximates
\(\mathcal G\).

This convention supersedes descriptions based on separate exterior and interior
control-point locations. The exterior self-potential, when needed for an operator
comparison, is obtained by subtracting the full doublet jump from the stored interior
self value; it is not obtained by moving the centroid.

## Discrete identity and its limits

Let the vector \(\boldsymbol q\) contain samples or coefficients for the trace, and
let \(\boldsymbol\sigma\) contain source-panel coefficients for its negative outward
normal derivative. With no attached-wake coupling in the closed-surface operator, the
discrete target identity is

\[
\boxed{
\mathcal S_h\boldsymbol\sigma
\approx(I-\mathcal B_h)\boldsymbol q.
}
\]

The approximation becomes an exact algebraic equality only if all of the following
hold:

- the trace and normal flux are representable in the chosen panel spaces;
- their coefficients are obtained consistently rather than mixing projection and
  point sampling;
- the same geometry, quadrature, kernel signs, self limits, and gauge are used in both
  terms; and
- evaluation errors, including near-singular wake influence, are absent.

Constant-strength panels filled from centroid samples generally do not meet the first
two conditions for a spatially varying wake field. The continuum Green identity is
exact, while its point-collocated constant-panel realization is normally a consistency
check that should converge under appropriate refinement.

The attached-wake term must not be folded into the body Green operator when checking
this identity. If the assembled lifting-body matrix is used without separating it,
then

\[
(I-\mathcal G)\boldsymbol q
=(I-\mathcal B)\boldsymbol q-\mathcal W\mathcal C\boldsymbol q,
\]

and the last term tests Kutta coupling rather than the closed-boundary identity.

## Gauge and the constant mode

Neumann data determines a harmonic potential only up to an additive constant. The
trace may therefore be replaced by

\[
\boldsymbol q^{\prime}=\boldsymbol q+a\boldsymbol 1.
\]

For a closed doublet shell, the constant mode is exterior-invisible. In the exact
operator,

\[
(I-\mathcal B)\boldsymbol 1=\boldsymbol 0.
\]

Numerically, the corresponding mode may be only approximately null because of
geometry and quadrature error. A comparison should impose one common gauge, for
example zero mean or a pinned reference value, rather than compare raw potentials.

The Kutta map uses differences of upper and lower trailing-edge values and therefore
annihilates an additive constant:

\[
\mathcal C\boldsymbol 1=\boldsymbol 0.
\]

Consequently, gauge selection does not change attached-wake strength. A nonconstant
difference between paired trailing-edge trace samples is different: it produces
\(\mathcal C\boldsymbol q\) and must be retained in the finite-wake derivation.

## Solver-path implications

In a Dirichlet solve, `set_strengths!` maps the velocity at each centroid to a source
strength. Production wake marching in `_steady_aerodynamics!` requests free-wake
velocity and does not request scalar potential.

The single-body `solve!` method saves `body.potential`, clears it before assembling the
body source potential, and restores it after solving. It therefore does not consume a
previously assembled external potential. In contrast, the coupled `KrylovCoupled` and
`BackslashCoupled` solvers save and re-add external potential before forming their
Dirichlet right-hand sides. This distinction matters to an experimental explicit-
potential solve, but not to the production wake path because that path assembles no
free-wake potential in the first place.

## Updated validation design

The most direct closed-surface diagnostic uses a body without Kutta-coupled wake
columns:

1. Construct a closed, consistently oriented source/doublet body.
2. Place a known wake singularity system wholly outside the body.
3. At every exact body centroid, evaluate its potential trace and velocity using the
   same kernel and regularization intended for the comparison.
4. Form source strengths from the negative outward normal velocity.
5. Evaluate the body single-layer potential from those source strengths.
6. Assemble the body-only exact-centroid Dirichlet operator, retaining its positive
   half self-potential limit.
7. Apply the same gauge to both sides and check
   \[
   \mathcal S_h\boldsymbol\sigma
   \approx(I-\mathcal B_h)\boldsymbol q.
   \]
8. Repeat under body-panel refinement and increasing wake/body separation to separate
   representation error from near-singular evaluation error.

Additional checks should be kept distinct:

- Adding attached-wake columns and measuring \(\mathcal C\boldsymbol q\) diagnoses
  Kutta-trace coupling.
- Comparing finite attached and first-row edges at off-body probes diagnoses handoff
  geometry, strength, orientation, and regularization.
- Comparing the complete marched composite with a semi-infinite wake at common probes
  diagnoses wake-shape and far-closure mismatch.

`test/dirichlet_potential_test.jl` does not currently implement this procedure. It is
an ad hoc script outside the maintained test matrix, refers to the removed `cp_outer`
constructor keyword and field, and presently fails before its final Green-identity
assertion. No passing regression should be claimed from that file until it is rewritten
for the exact-centroid operator and admitted to a maintained test suite.
