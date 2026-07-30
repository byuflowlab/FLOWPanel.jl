# Unsteady Pressure Monitors: ALE Derivation and Discrete Design

Date: 2026-07-11

## Scope and notation

This note defines the mathematical contract for revising `PressureBernoulli`
and `PressureLaplace`. It distinguishes quantities that are easy to conflate in
a moving-body panel calculation:

- `u`: inertial fluid velocity at the body control point;
- `w`: inertial velocity of the moving control-point grid/body;
- `q = u - w`: fluid velocity relative to that grid;
- `n`: the current outward panel normal;
- `P = I - n n^T`: projection onto the current tangent plane;
- `G = grad(u)`, stored with `G[k,l] = partial(u_k)/partial(x_l)`.

In FLOWPanel, after `kinematic_velocity!`, `body.velocity` is `q` and
`body.velocity_kinematic` is `w`. The reconstructed surface trace used for
time history is

```math
u_s = w + Pq.
```

For an exactly satisfied impermeability condition, `q dot n = 0` and `u_s=u`.
The projection removes only numerical normal leakage before it is amplified by
a temporal difference. It is not a claim that the inertial velocity itself is
purely tangential.

## Why scalar-potential Bernoulli is incomplete with particle wakes

For globally irrotational flow, `u = grad(phi)` and Euler's equation integrates
to

```math
p + rho partial_t(phi) + (rho/2)|u|^2 = C(t).
```

This scalar-potential Bernoulli relation is not a complete pressure equation in
a rotational region. A vortex-particle field is represented by a vector
potential/curl contribution and generally has nonzero vorticity. There is no
single-valued scalar `phi` whose gradient reproduces the particle-induced
velocity throughout that field. Excluding particle sources from `phi` while
retaining their velocity in `|u|^2` therefore produces a useful partial
diagnostic, not the full unsteady Bernoulli pressure of the coupled flow.

`PressureBernoulli(unsteady=true)` must consequently:

1. evaluate scalar potential locally to the monitor, only from body and panel-
   wake sources for which the FastMultipole interface supplies scalar
   potential, plus the freestream potential;
2. never request scalar potential from vortex particles or other vector-
   potential-only sources;
3. reject vector-potential-capable body sources, whose retained exterior scalar
   trace cannot be separated reliably;
4. reject vector-potential-only wake sources by default;
5. only with `allow_partial=true`, emit one prominent warning per monitor
   instance and continue with a partial diagnostic.

The warning is required even though filtering the sources is numerically safe:
silently presenting the partial result as a complete Bernoulli pressure would
be physically misleading. The opted-in output is a mixed partial diagnostic:
particle velocity remains in inertial kinetic energy while particle potential
history is absent. Its velocity is subtracted from `w dot grad(phi)` so the ALE
contraction uses only retained scalar flow. Detect exclusion from the complete
body and wake source lists, because a prefiltered scalar-source collector has
already discarded capability information.
Steady Bernoulli does not use history and need not inspect wake capabilities.

## Moving-point derivative for Bernoulli

Potential is sampled at a control point `x_p(t)` moving with velocity `w`:

```math
D_g phi := d phi(x_p(t),t)/dt = partial_t phi + w dot grad(phi).
```

Hence the Eulerian derivative required by scalar Bernoulli is

```math
partial_t phi = D_g phi - w dot u.
```

Here ALE means Arbitrary Lagrangian--Eulerian. In an opted-in partial diagnostic,
`u` in this contraction is the retained scalar-potential velocity, while the
kinetic-energy term keeps total inertial exterior velocity.

The monitor must finite-difference the panel-following scalar history and then
subtract `w dot u`. It must also use reconstructed inertial `u_s`, rather than
stored relative `q`, in Bernoulli's kinetic-energy term. Combining inertial
`partial_t(phi)` with `|q|^2/2` is frame-inconsistent. On the first sample the
monitor initializes history and reports zero temporal contribution. With no
earlier state, inventing `phi_old=0` creates an `O(phi/dt)` startup spike with no
mathematical meaning. A Galilean-translation test must verify both terms.

## ALE material acceleration

Euler's equation in an inertial frame is

```math
rho [partial_t u + (u dot grad)u] = -grad(p).
```

The derivative of the inertial velocity sampled by the moving panel is

```math
D_g u = partial_t u + (w dot grad)u.
```

Eliminating `partial_t u` gives the ALE identity

```math
boxed{a = D_g u + [(u-w) dot grad]u = D_g u + G q.}
```

At an impermeable surface, only the tangent relative velocity transports the
surface trace. Numerically the monitor therefore uses

```math
q_t = P q,
qquad
a_t = P [D_g u_s + G_corrected q_t].
```

The history is a history of the inertial surface velocity `u_s`, not of `q`.
Differencing `q` would mix physical fluid acceleration with changes in rigid
grid velocity, while using `q` rather than `q_t` in the convective direction
would amplify boundary-condition residuals.

## Surface projection: what is and is not projected

The surface Euler relation is

```math
grad_s(p) = -rho P a.
```

Projection belongs on the relative transport direction and on the final total
acceleration. A blanket replacement

```math
G -> P G P
```

is inappropriate for the corrected-Hessian formulation:

- right multiplication by `P` is already effected by using `q_t`;
- left multiplication prematurely discards the normal component of each
  acceleration contribution before temporal and convective terms are added;
- `P` varies over a curved surface, so differentiating the projected field
  requires curvature/product-rule terms that a pointwise `PGP` does not add;
- the exterior surface trace includes a jump velocity whose tangent derivative
  must be added explicitly, not hidden by projecting an off-surface Hessian.

The implementation forms the unprojected vector sum first and applies `P` once
to the final acceleration. The surface-velocity alternative instead
reconstructs directional derivatives of the full vector surface trace; those
derivatives already have tangent input directions, but their vector outputs
are again projected only after forming total acceleration.

## Temporal discretization and restart behavior

Let `h_n=t_n-t_(n-1)` and `h_(n-1)=t_(n-1)-t_(n-2)`. After history has been
initialized, the first available derivative uses backward Euler,

```math
D_g f^n = (f^n-f^(n-1))/h_n + O(h_n).
```

Once two previous samples exist, variable-step BDF2 is

```math
D_g f^n approx
  (2h_n+h_(n-1))/(h_n(h_n+h_(n-1))) f^n
- (h_n+h_(n-1))/(h_n h_(n-1)) f^(n-1)
+ h_n/(h_(n-1)(h_n+h_(n-1))) f^(n-2).
```

Equivalently, for `r=h_n/h_(n-1)`,

```math
D_g f^n approx (1/h_n)[
  (1+2r)/(1+r) f^n - (1+r)f^(n-1) + r^2/(1+r)f^(n-2)].
```

For equal steps this reduces to `(3f^n-4f^(n-1)+f^(n-2))/(2h)`.
Both `PressureBernoulli` and unsteady `PressureLaplace` use this history policy.
The first-ever call returns zero derivative, the second uses backward Euler,
and subsequent calls use variable-step BDF2.

History validity must be explicit. It is invalidated when body count, panel
count, or a nonconsecutive step sequence changes, and when a monitor is reused
for a fresh simulation beginning at or before its last observed step. An
invalid history is seeded from the current sample and yields zero derivative.
Warm-starting without serialized monitor history follows the same safe policy:
the restart sample seeds the buffers, the next sample uses backward Euler, and
only then is BDF2 enabled. This loses at most one step of temporal order but
never creates a false derivative. If monitor histories are later serialized,
they may be restored only together with their sample times/step identity.

In the current `simulate!` interface, the `dt` passed while sampling `t_n` is
normally the forward interval `t_(n+1)-t_n`. The derivative at the next call
must therefore use the interval saved by the preceding call as `h_n`; it must
not use the new call's forward `dt` on a nonuniform grid. Two saved forward
intervals supply `h_n` and `h_(n-1)` for BDF2.

All runtime steps must be positive and both BDF2 intervals must be finite and
positive. A repeated/nonmonotone step sequence invalidates and reseeds history
rather than applying invalid coefficients.

## FastMultipole velocity gradients and rigid motion

The analytic kernel output is the inertial-coordinate Jacobian of induced
velocity at the target:

```math
G_kernel[k,l] = partial(u_induced,k)/partial(x_l).
```

Direct and FMM backends must populate the same layout and sign. The rigid grid
velocity is

```math
w(x) = V_O + Omega cross (x-x_O),
```

and therefore

```math
grad(w) = [Omega]_cross,
```

where `[Omega]_cross v = Omega cross v`. This identity is useful for auditing
frames, but **it must not be added to `body.velocity_gradient`**. The current
`influence!` calls accumulate `grad(u_induced)` directly, while
`kinematic_velocity!` subtracts `w` only from velocity values and never
subtracts `grad(w)` from the gradient buffer. With zero freestream gradient,
the stored buffer is already the inertial kernel-resolved gradient. Adding
`[Omega]_cross` would double-count grid rotation. Uniform inertial flow sampled
on a rotating grid is decisive: physical `G=0` and `a=0` even though `q` varies.

The Hessian/velocity-gradient requirement must be configuration-specialized.
A corrected-Hessian monitor requests it; edge-difference and reconstructed-
surface modes do not. This preserves the existing compile-time derivative
switches in the FastMultipole kernels and avoids both storage and kernel work
for modes that cannot consume a Hessian.

Piecewise filament/ray kernels differentiate the active analytic branch, not
the derivative of a cutoff switch. Their gradient is undefined at a cutoff
transition; validation points must stay away from those transitions.

## Exterior surface jump and its gradient

For the doublet/vortex-sheet representation used by lifting bodies, the
off-surface principal-value kernel velocity is not by itself the exterior
surface trace. In terms of stored strength `mu_code`, FLOWPanel postprocessing
adds

```math
u_jump = -(1/2) grad_s(mu_code)
```

This is the concrete `compute_mu_gradient!(...; scale=0.5)` convention: the
helper stores `-scale*grad_s(mu_code)`. A textbook `mu` with the opposite sign
must be translated before using this equation.

The potential history itself also requires the exterior surface limit. The
exact-control-point body evaluation returns FLOWPanel's canonical interior
limit, so every target panel of a body with `has_grad_mu(body)` stores

```math
phi^+ = phi^- - mu_code.
```

The local strength is subtracted exactly once. Source-only bodies, other
bodies' off-surface contributions, and wake potentials are unchanged.

Trailing-edge panel-pressure averaging is an optional heuristic. It is disabled
by default and enabled explicitly with `correct_kuttacondition=true`.

The raw kernel Hessian cannot differentiate a term that was added only after
the influence evaluation. Thus the corrected gradient is

```math
boxed{G_corrected = G_kernel + grad_s(u_jump).}
```

`grad_s(u_jump)` means the tangent directional derivative of the three-
component jump-velocity field. It can be obtained by first reconstructing the
same `u_jump` used by velocity postprocessing and then applying
`compute_surface_velocity_gradient!` to that vector field, with trailing-edge
isolation and the same `grad_mu` basis/options. It is not another factor of
one-half: the half jump is already present in `u_jump`.

Source-only bodies have zero jump correction. At a trailing edge the two
surface traces must not be blended across the sheet; the established shedding
metadata must constrain both reconstructions. This correction should be
tested by prescribing a manufactured `mu` with a nonconstant surface gradient,
not merely a linear `mu` whose jump gradient vanishes.

## Pressure-Laplace discretization

On the surface,

```math
-Delta_s p = rho div_s(a_t).
```

The existing symmetric finite-volume operator uses one shared-edge weight per
interior edge. Let `r_ij=x_j-x_i`, shared-edge length `ell`, and `nu_ij` be the
averaged-surface co-normal oriented from panel `i` to `j`. The two-point weight
is

```math
w_ij = ell (nu_ij dot r_ij)/|r_ij|^2.
```

For a panel-centered acceleration, the conservative pair flux is

```math
F_ij = rho w_ij [ (a_t,i+a_t,j)/2 dot r_ij ],
```

added to row `i` and subtracted from row `j`. This is conservative, but it
equals `rho ell a dot nu_ij` exactly only for orthogonal center-to-center
geometry (or equivalent reconstruction assumptions). Skew-mesh validation
must compare with a known divergence, not only another monitor mode. This is
the corrected-Hessian default: construct ALE acceleration at panels, project
the final vector, then apply the pair flux. The gauge row remains pinned.

## Alternative acceleration formulations

### Corrected Hessian (recommended default)

```math
a_t = P { BDF(u_s) +
          [G_kernel + grad_s(u_jump)] q_t }.
```

This retains analytic spatial derivatives for kernel-resolved velocity, adds
the missing exterior-jump derivative, and performs one final projection. Its
limitations include near-singular Hessian sensitivity and the fact that the
jump correction alone does not prove a regularized off-surface Hessian has the
correct principal-value exterior limit; refinement must verify that limit.

### Corrected edge difference

The naive contraction `q dot (u_s,j-u_s,i)` approximates `q^T G r`, while the
required scalar is `(Gq) dot r = r^T Gq`. They differ for a nonsymmetric
gradient. With repository layout `G-G^T=[omega]_cross`,

```math
r^T Gq = q^T Gr + r dot (omega cross q).
```

The rotationally corrected edge diagnostic is

```math
a_ij dot r_ij approx
  [(D_g u_s,i + D_g u_s,j)/2] dot r_ij
+ q_t,ij dot (u_s,j-u_s,i)
+ (omega_ij cross q_t,ij) dot r_ij.
```

Here `q_t,ij` and `omega_ij` are consistent edge averages. This is conservative
and Hessian-free but needs vorticity of the same exterior trace. Without the
last term it is valid only for locally symmetric/irrotational `G`. Surface-
sheet vorticity conventions make it a cross-check, not a default. A
nonsymmetric manufactured `G` must distinguish the corrected and naive forms.
Accuracy still degrades on anisotropic/non-orthogonal meshes.

### Surface-velocity reconstruction

Reconstruct `grad_s(u_s)` from neighboring panel values, then use

```math
a_t = P [BDF(u_s) + grad_s(u_s) q_t].
```

This automatically includes rigid motion and the exterior jump because both
are already in `u_s`. It is Hessian-free and is a valuable cross-check, but its
least-squares/paired-triangle reconstruction is mesh-sensitive, especially at
high aspect ratio and at stencil boundaries/trailing edges.

### Lamb-vector diagnostic

In an inertial Eulerian frame,

```math
(u dot grad)u = grad(|u|^2/2) + omega cross u.
```

Substituting an ALE panel-following derivative requires undoing the grid
transport term, and replacing `u` by the tangent relative velocity inside only
parts of the identity is not algebraically equivalent. Surface-sheet
vorticity also has a distributional jump whose discrete use is convention-
sensitive. Until a complete ALE surface derivation identifies consistent
kinetic-energy, grid-transport, and sheet-vorticity terms, Lamb-vector mode is
retained only as a diagnostic and is not a default or validation oracle.

## Type-specialized monitor requirements

Monitor configuration should encode whether each expensive data path is used:

- steady Bernoulli: pressure only;
- unsteady Bernoulli: monitor-local scalar potential and scalar history;
- corrected Hessian: velocity history when unsteady, kernel Hessian, and jump
  reconstruction only for bodies that possess the doublet surface jump;
- corrected edge difference: velocity history when unsteady, no Hessian, and
  vorticity for rotational consistency;
- surface velocity: velocity history when unsteady, no Hessian/vorticity;
- Lamb diagnostic: velocity history when unsteady and only its selected
  vorticity source (plus Hessian solely for a Hessian-curl diagnostic).

Traits used by `simulate!` may remain the public query surface, but their
answers must dispatch on specialized monitor configuration rather than force
unused potential, Hessian, vorticity, or history paths through runtime branches.

## Defaults, limitations, and validation gates

Recommended defaults are:

- `PressureBernoulli`: steady unless explicitly requested; unsteady mode warns
  rejects vector-only wake sources by default, or warns once and continues only
  when `allow_partial=true`;
- `PressureLaplace`: corrected-Hessian acceleration, conservative edge
  divergence, unsteady term opt-in as today, BDF history with safe startup;
- corrected edge difference and surface-velocity reconstruction: supported
  alternatives and cross-checks;
- Lamb vector: diagnostic only.

The following evidence is required before treating the revision as complete:

1. manufactured ALE translation with prescribed `u(x,t)` distinguishes `u`,
   `w`, and `q` and recovers `D_g u + Gq`;
2. uniform inertial velocity sampled on a rotating grid verifies that no
   `[Omega]_cross` is added to the already-inertial gradient buffer;
3. a curved/skewed surface case verifies that final `P a`, rather than blanket
   `PGP`, removes normal acceleration without erasing required tangent terms;
4. a nonconstant doublet-strength field verifies the exterior half-jump
   velocity gradient and a trailing-edge case verifies stencil isolation;
5. equal-step refinement demonstrates second-order BDF2 convergence after
   startup, and a nonuniform sequence verifies the variable-step coefficients;
6. first call, monitor reuse, topology mismatch, nonconsecutive steps, and warm
   restart produce zero seeding derivatives rather than spikes;
7. a particle-wake Bernoulli case throws by default, emits exactly one warning
   per opted-in monitor, and continues with finite partial pressure; a
   vector-capable body throws even in partial mode;
8. direct and FMM kernel gradients agree for supported elements;
9. corrected Hessian, corrected edge, and surface reconstruction are compared
   on isotropic and anisotropic meshes, at trailing edges, and in the pitching-
   wing regression. Agreement is expected under refinement, not necessarily
   pointwise equality on coarse meshes.

Validation cases must be discriminating. Uniform inertial velocity on a
rotating grid exposes a spurious kinematic-gradient addition; linear `mu`
cannot expose a missing jump Hessian; a planar surface cannot expose the
projection/curvature distinction; uniform timesteps alone cannot validate
variable-step BDF2; and a panel-only wake cannot exercise the particle warning.
A Galilean translation must distinguish inertial from relative Bernoulli
kinetic energy, and nonsymmetric `G` must distinguish `r^T Gq` from `q^T Gr`.
