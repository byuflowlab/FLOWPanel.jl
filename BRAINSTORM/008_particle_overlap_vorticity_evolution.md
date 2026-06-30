# Particle-Overlap Vorticity Evolution

## Context

Particle-field wakes often use Pedrizzetti or corrected-Pedrizzetti relaxation to
stabilize the vortex-particle field. In the current rotor-hover tracks this has
been treated as a numerical damping or alignment mechanism, but it may also be
compensating for a missing overlap-coupled particle-strength evolution equation.

The hypothesis here is that the usual per-particle stretching update implicitly
assumes that particle filter supports do not overlap strongly. When particles do
overlap, the governing vorticity equation constrains the filtered field as a
whole, not each particle coefficient independently. That suggests solving a
small/large overlap system for `DΓ/Dt` may be more consistent than relaxing the
field after the fact.

## Starting Point

The proposed filtered vorticity representation is:

```math
\vec{\omega}(\vec{x}, t) =
\sum_i \vec{\Gamma}_i(t)\,
\zeta_{\sigma_i}(\vec{x} - \vec{x}_i(t)).
```

For inviscid, incompressible flow with viscosity, baroclinic torque, and
compressibility terms neglected, the vorticity transport equation is:

```math
\frac{D\vec{\omega}}{Dt} =
(\vec{\omega}\cdot\nabla)\vec{u}.
```

These equations are the right starting point for this investigation, with two
notes:

- "particle i" is intended, not "partial i".
- The second equation is the inviscid incompressible form. If viscous core
  spreading or SFS terms are retained, they must appear as additional right-hand
  side terms or as explicit `\dot{\sigma}_i`/model terms in the derivation.

## Theory Task: Generic Kernel System

Let

```math
\vec{\omega}(\vec{x},t)
= \sum_i \vec{\Gamma}_i(t)\zeta_i(\vec{x},t),
\qquad
\zeta_i(\vec{x},t)=\zeta_{\sigma_i(t)}(\vec{x}-\vec{x}_i(t)).
```

The flow material derivative at a fixed field point `\vec{x}` is

```math
\frac{D}{Dt} = \frac{\partial}{\partial t}
              + \vec{u}(\vec{x},t)\cdot\nabla.
```

First apply the product rule:

```math
\frac{D\vec{\omega}}{Dt}
= \sum_i
  \left(
    \zeta_i\frac{d\vec{\Gamma}_i}{dt}
    + \vec{\Gamma}_i\frac{D\zeta_i}{Dt}
  \right).
```

The kernel changes because the field point moves through the flow, the particle
center moves with velocity `d\vec{x}_i/dt`, and the kernel width may change:

```math
\frac{D\zeta_i}{Dt}
=
\frac{\partial\zeta_i}{\partial t}
+ \vec{u}(\vec{x},t)\cdot\nabla\zeta_i
=
-\frac{d\vec{x}_i}{dt}\cdot\nabla\zeta_i
+ \dot{\sigma}_i\frac{\partial\zeta_i}{\partial\sigma_i}
+ \vec{u}(\vec{x},t)\cdot\nabla\zeta_i.
```

Combining the two gradient terms gives:

```math
\frac{D\vec{\omega}}{Dt}
=
\sum_i \zeta_i \frac{d\vec{\Gamma}_i}{dt}
+ \sum_i \vec{\Gamma}_i
  \left(\vec{u}(\vec{x}) - \frac{d\vec{x}_i}{dt}\right)\cdot\nabla\zeta_i
+ \sum_i \vec{\Gamma}_i \dot{\sigma}_i
  \frac{\partial \zeta_i}{\partial \sigma_i}.
```

For Lagrangian particles, `d\vec{x}_i/dt = \vec{u}_i = \vec{u}(\vec{x}_i)`, so the
overlap-aware strength equation should satisfy:

```math
\sum_i \zeta_i(\vec{x}) \frac{d\vec{\Gamma}_i}{dt}
=
(\vec{\omega}(\vec{x})\cdot\nabla)\vec{u}(\vec{x})
- \sum_i \vec{\Gamma}_i
  (\vec{u}(\vec{x}) - \vec{u}_i)\cdot\nabla\zeta_i
- \sum_i \vec{\Gamma}_i \dot{\sigma}_i
  \frac{\partial \zeta_i}{\partial \sigma_i}.
```

The concrete research task is to project or collocate this equation into a linear
system for all particle strength derivatives:

```math
M \dot{\Gamma} = b,
```

where `M` is a scalar kernel-overlap matrix and each row of `b` is a 3-vector.
Candidate discretizations:

- Collocation at particle centers: `M_{ki} = \zeta_i(\vec{x}_k)`.
- Galerkin or least-squares projection with particle/filter test functions.
- Local sparse solves over neighbor lists rather than a global dense solve.

The right-hand side should use the already available particle velocity and
velocity-gradient fields where possible. With the common Jacobian convention
`J_{ab} = \partial u_a / \partial x_b`, the stretching term is `J \omega`.

The first task is purely theoretical and should not depend on FLOWVPM internals:
pose the linear system formally for a generic compact or rapidly decaying kernel
`\zeta`. This should specify the unknown vector, the test/collocation space, the
matrix entries, the right-hand side, and the assumptions needed for solvability
or approximation. It should also make clear which terms are exact consequences
of the filtered vorticity ansatz and which terms are modeling choices, such as
kernel-width evolution, local truncation, mass lumping, or regularization.

## Consistency Checks

- **No-overlap limit:** if only particle `i` contributes at `\vec{x}_i`, and
  `\dot{\sigma}_i = 0`, then the equation should reduce to the standard
  isolated-particle stretching update:

  ```math
  \frac{d\vec{\Gamma}_i}{dt} = J_i \vec{\Gamma}_i.
  ```

- **Equal-velocity limit:** if `\vec{u}(\vec{x}) \approx \vec{u}_i` over each
  particle support, the convective-basis correction should vanish, leaving only
  overlap-coupled stretching and any sigma-growth terms.
- **Core-spreading compatibility:** if viscous core spreading changes `\sigma_i`,
  the `\dot{\sigma}_i \partial\zeta_i/\partial\sigma_i` term must either be
  included or explicitly shown to be handled by the existing viscous model.
- **Conditioning:** high overlap, particle merging, and very different `\sigma_i`
  may make `M` ill-conditioned. The item should identify whether regularization,
  local solves, or mass-lumped approximations are needed.

## Relation to Existing Wake Work

This item is a theory/prototype follow-up to the particle-wake stability work in
items 005 and 006. In item 005, turning relaxation off worsened the stable-wake
baseline, which suggests relaxation is doing useful work. This item asks whether
that useful work can be replaced or reduced by enforcing the filtered vorticity
transport equation with particle overlap included.

After the generic-kernel theory task, the next phase can specialize the proposed
system to the FLOWVPM dependency. The likely implementation landing zone, if the
derivation survives, is the particle propagation path around `PanelParticleWake`
and the FLOWVPM Euler update, where particle velocity, velocity gradient,
relaxation, viscous spreading, SFS, and maintenance are already coordinated.

## Proposed Work

1. Formal theory task: pose `M \dot{\Gamma} = b` for a generic kernel `\zeta`,
   independent of FLOWVPM. Include fixed-`\sigma` first, then the optional
   `\dot{\sigma}` term, and state the test/collocation choice and solvability
   assumptions.
2. Specialize the formal system to the kernels, particle state fields, and
   update hooks exposed by FLOWVPM/FLOWPanel. Decide how existing core spreading,
   SFS, relaxation, and maintenance steps enter the right-hand side or remain
   separate model terms.
3. Prototype an offline diagnostic on saved particle states:
   - assemble local overlap matrices,
   - compare branch-aware isolated updates
     `\mathcal{S}_\mathrm{FLOWVPM}(J_i,\Gamma_i)` against the overlap-solved
     `\dot{Γ}_i`, while separately retaining `J_i\Gamma_i` as the physical
     vorticity-equation operator,
   - measure residual in the field equation before and after the solve.
4. If promising, add an opt-in wake update mode and compare against
   Pedrizzetti/corrected-Pedrizzetti relaxation on the rotor-hover stable-wake
   cases.

## Acceptance

This item is technically complete when it produces either:

- a defensible overlap-coupled `DΓ/Dt` equation with a prototype showing improved
  vorticity-equation residuals or reduced reliance on relaxation, or
- a clear negative result explaining why the overlap solve is ill-posed,
  redundant with existing FLOWVPM machinery, or too expensive for practical wake
  runs.

Any implementation follow-up must report effects on rotor-hover CT mean, CT
ripple, particle-count/conditioning diagnostics, and comparison to the existing
relaxation-on baseline.

## Caveats

Do not treat this as a validated replacement for Pedrizzetti relaxation until it
has been tested on saved particle fields and then in live wake propagation. The
derivation should also keep sign and Jacobian conventions explicit, since a
transpose error in the stretching term would look like a plausible but wrong
stabilization scheme.

## 2026-06-23 Initial Theory Phase

### Notation and assumptions

Use the filtered representation

```math
\omega(x,t) = \sum_{i=1}^{N_p} \Gamma_i(t)\zeta_i(x,t),
\qquad
\zeta_i(x,t) = \zeta_{\sigma_i(t)}(x-X_i(t)),
```

where each `\Gamma_i` is a 3-vector coefficient and each `\zeta_i` is a scalar
kernel. Assume the kernel is smooth enough to evaluate `\nabla\zeta_i`, compact
or rapidly decaying enough to truncate neighbor sums, and normalized
consistently with the meaning of `\Gamma_i`. For the first pass, take particle
centers as Lagrangian markers:

```math
\dot{X}_i = u_i = u(X_i).
```

With the Jacobian convention

```math
J_{ab} = \frac{\partial u_a}{\partial x_b},
```

the inviscid incompressible stretching term is

```math
(\omega\cdot\nabla)u = J\omega.
```

The modeled incompressible vorticity equation used for the expanded RHS is:

```math
\frac{D\omega}{Dt}
= J\omega + \nu\nabla^2\omega + q_\mathrm{SFS}.
```

For the purely inviscid fixed-core limit, set `\nu=0`, `q_\mathrm{SFS}=0`, and
`\dot{\sigma}_i=0`. Viscosity, SFS, core spreading, or filtering corrections are
model terms layered on top of the inviscid overlap algebra. They must be
included consistently: fixed-core viscosity enters through `\nu\nabla^2\omega`,
while core-spreading viscosity enters through the explicit
`\dot{\sigma}_i\partial_{\sigma_i}\zeta_i` basis-evolution term.

### Strong filtered equation

Start from the modeled vorticity equation and the filtered ansatz:

```math
\frac{D\omega}{Dt}
= J\omega + \nu\nabla^2\omega + q_\mathrm{SFS},
\qquad
\omega(x,t)=\sum_i \Gamma_i(t)\zeta_i(x,t).
```

Here `D/Dt = \partial_t + u(x,t)\cdot\nabla` is the material derivative of the
field, not the time derivative following particle `i`. Applying it to the
filtered ansatz gives

```math
\frac{D\omega}{Dt}
= \sum_i
  \left(
      \zeta_i\dot{\Gamma}_i
    + \Gamma_i\frac{D\zeta_i}{Dt}
  \right).
```

For `\zeta_i(x,t)=\zeta_{\sigma_i(t)}(x-X_i(t))`, the chain rule gives

```math
\frac{\partial\zeta_i}{\partial t}
=
-\dot{X}_i\cdot\nabla\zeta_i
+ \dot{\sigma}_i\frac{\partial\zeta_i}{\partial\sigma_i}.
```

Therefore

```math
\frac{D\zeta_i}{Dt}
=
\frac{\partial\zeta_i}{\partial t}
+ u(x,t)\cdot\nabla\zeta_i
=
\left(u(x,t)-\dot{X}_i\right)\cdot\nabla\zeta_i
+ \dot{\sigma}_i\frac{\partial\zeta_i}{\partial\sigma_i}.
```

With Lagrangian particle centers, `\dot{X}_i=u_i=u(X_i,t)`, so

```math
\frac{D\omega}{Dt}
=
\sum_i \zeta_i(x)\dot{\Gamma}_i
+ \sum_i \Gamma_i\left(u(x)-u_i\right)\cdot\nabla\zeta_i(x)
+ \sum_i \Gamma_i\dot{\sigma}_i
    \frac{\partial\zeta_i}{\partial\sigma_i}(x).
```

Equating this expression with the modeled RHS and moving the basis-evolution
terms to the right gives

```math
\sum_i \zeta_i(x)\dot{\Gamma}_i
=
J(x)\omega(x)
+ \nu\nabla^2\omega(x)
+ q_\mathrm{SFS}(x)
- \sum_i \Gamma_i\left(u(x)-u_i\right)\cdot\nabla\zeta_i(x)
- \sum_i \Gamma_i\dot{\sigma}_i
    \frac{\partial\zeta_i}{\partial\sigma_i}(x).
```

This is the exact algebraic consequence of the filtered ansatz, the particle
center kinematics, and the chosen kernel-width evolution. It is not yet a
closed discrete method, because the equation is a vector field identity while
the unknowns are the finite set of particle coefficients `\dot{\Gamma}_i`.

Define

```math
r(x;\dot{\Gamma}) =
\sum_i \zeta_i(x)\dot{\Gamma}_i - f(x),
```

where

```math
f(x) =
J(x)\omega(x)
+ \nu\nabla^2\omega(x)
+ q_\mathrm{SFS}(x)
- \sum_i \Gamma_i\left(u(x)-u_i\right)\cdot\nabla\zeta_i(x)
- \sum_i \Gamma_i\dot{\sigma}_i
    \frac{\partial\zeta_i}{\partial\sigma_i}(x).
```

The overlap-aware strength update is a projection of `r=0` into a finite test
space.

### Viscous and SFS RHS terms

The explicit viscous term depends on how viscosity is modeled. If particle
cores are fixed, then

```math
\nu\nabla^2\omega(x)
=
\nu\sum_i \Gamma_i\nabla^2\zeta_i(x),
```

and this term should be projected directly into `b`.

If viscosity is represented by core spreading, the diffusion contribution is
instead carried by the basis-evolution term. For normalized Gaussian kernels,

```math
G_\sigma(x)
= (2\pi\sigma^2)^{-3/2}
   \exp\left(-\frac{|x|^2}{2\sigma^2}\right),
\qquad
\frac{\partial G_\sigma}{\partial\sigma}
= \sigma\nabla^2G_\sigma.
```

The heat equation identity `\partial_t G_\sigma=\nu\nabla^2G_\sigma` is
therefore reproduced by

```math
\dot{\sigma}_i = \frac{\nu_{\mathrm{eff},i}}{\sigma_i}.
```

With that choice, the explicit diffusion RHS and the
`-\Gamma_i\dot{\sigma}_i\partial_{\sigma_i}\zeta_i` term cancel exactly for the
part of diffusion represented by core spreading. Do not include both
`+\nu\nabla^2\omega` and `\dot{\sigma}` as independent corrections unless the
viscous model is known to leave a residual diffusion term. In FLOWVPM terms,
`CoreSpreading` should first be interpreted as a `\dot{\sigma}_i` model, not as
an additional standalone RHS source.

Represent SFS as a particle strength-rate source `S_i` when the SFS model
provides one:

```math
q_\mathrm{SFS}(x) \approx \sum_i S_i\zeta_i(x).
```

Then its contribution to a collocation RHS is

```math
b^\mathrm{SFS}_k = \sum_i \zeta_i(x_k)S_i = (MS)_k.
```

For the particle-kernel Galerkin projection, the same basis representation gives

```math
b^\mathrm{SFS}_k
= \sum_i \left(\int \zeta_k(x)\zeta_i(x)\,dx\right) S_i
= (MS)_k.
```

Thus SFS is naturally handled by the same overlap operator if the stored SFS
quantity is a particle coefficient-rate source. This must be verified against
FLOWVPM's SFS convention before using the saved `SFS` field as `S_i`; until then
it should be treated as a candidate source term, not a proven RHS contribution.

### Collocation system

Collocate at points `x_k`, usually particle centers. The scalar overlap matrix
and 3-vector right-hand side are

```math
M_{ki} = \zeta_i(x_k),
\qquad
b_k = f(x_k).
```

Then solve

```math
\sum_i M_{ki}\dot{\Gamma}_i = b_k,
\qquad k=1,\ldots,N_p.
```

Equivalently, for components,

```math
M\dot{\Gamma}^{(a)} = b^{(a)}, \qquad a=1,2,3.
```

This is the simplest saved-state diagnostic because it needs only particle
positions, strengths, kernel widths, velocities, vorticity, and velocity
gradients. The downside is that the matrix is generally nonsymmetric and its
conditioning depends strongly on point placement and overlap. Center
collocation also drops most of the convective-basis correction for the self
particle because `u(x_i)-u_i=0`, while neighbor corrections remain.

### Galerkin or least-squares system

A more defensible field projection uses test functions `\phi_k`. Multiplying by
`\phi_k` and integrating gives

```math
M_{ki} = \int_{\Omega} \phi_k(x)\zeta_i(x)\,dx,
\qquad
b_k = \int_{\Omega}\phi_k(x) f(x)\,dx.
```

The natural particle-kernel Galerkin choice is `\phi_k=\zeta_k`, yielding the
Gram matrix

```math
M_{ki} = \int_{\Omega}\zeta_k(x)\zeta_i(x)\,dx.
```

For positive kernels with linearly independent translated/scaled basis
functions, this matrix is symmetric positive definite. This form directly
minimizes the field residual in the kernel-weighted `L^2` sense and is less
arbitrary than center collocation. Its cost is evaluating overlap integrals and
right-hand-side quadrature. For identical Gaussian kernels these integrals may
be analytic; for generic `\zeta` the practical path is neighbor-truncated
quadrature over local support.

The least-squares view is also useful. With quadrature points `x_q` and weights
`w_q`, define

```math
A_{qi} = \sqrt{w_q}\zeta_i(x_q),
\qquad
y_q = \sqrt{w_q}f(x_q).
```

Then solve

```math
\min_{\dot{\Gamma}}\|A\dot{\Gamma}-y\|_2^2,
```

again independently for the three vector components. The normal equations give
the Galerkin Gram matrix when the quadrature and tests are chosen consistently.

### No-overlap and local-velocity checks

For center collocation in the inviscid, no-SFS, fixed-`\sigma` limit, suppose
row `k=i` receives only the self contribution. Then

```math
M_{ii} = \zeta_i(X_i),
\qquad
\omega(X_i) = \Gamma_i\zeta_i(X_i),
\qquad
f(X_i) = J_i\Gamma_i\zeta_i(X_i).
```

The solve gives

```math
\dot{\Gamma}_i = J_i\Gamma_i,
```

so the isolated-particle stretching update is recovered. A Galerkin projection
has the same limit if disjoint supports make `M` diagonal and `J` is effectively
constant over the particle support:

```math
b_i \approx \left(\int \zeta_i^2\,dx\right)J_i\Gamma_i.
```

If `u(x)\approx u_i` across each particle support, the convective-basis term
vanishes. The remaining coupling is then purely the projection of `J\omega`
onto overlapping basis functions, plus whatever model terms are deliberately
included through viscosity, SFS, or `\dot{\sigma}`.

If an SFS particle source `S_i` is included and only the self particle
contributes at `X_i`, then

```math
f(X_i) =
\zeta_i(X_i)\left(J_i\Gamma_i + S_i\right),
```

so the no-overlap collocation update becomes

```math
\dot{\Gamma}_i = J_i\Gamma_i + S_i.
```

For Gaussian core spreading with `\dot{\sigma}_i=\nu_{\mathrm{eff},i}/\sigma_i`,
the fixed-core diffusion RHS and the basis-evolution term cancel for the
represented diffusion. The no-overlap update therefore remains
`J_i\Gamma_i` plus any non-viscous modeled source. This is the main
double-counting check for any later full-RHS diagnostic.

### Solvability and regularization

The overlap system is well posed only to the extent that the active particle
kernels form a numerically independent basis in the chosen test space. Practical
failure modes are:

- too many nearly coincident particles with similar `\sigma`, producing nearly
  duplicate columns;
- very large overlap, which smooths the field so strongly that many coefficient
  distributions represent nearly the same `\omega`;
- very different `\sigma_i`, which can create poor row and column scaling;
- local truncation that cuts important kernel tails and changes the effective
  matrix.

The first prototype should therefore report `cond(M)` or a cheap local
condition proxy, row sums, diagonal dominance, and the relative size of the
regularization. A defensible regularized form is

```math
(M + \lambda D)\dot{\Gamma}^{(a)} = b^{(a)}
```

for collocation, or

```math
(A^T A + \lambda D)\dot{\Gamma}^{(a)} = A^T y^{(a)}
```

for least squares. `D` should be a scale-aware diagonal, for example
`D_{ii}=M_{ii}` or `D_{ii}=(A^TA)_{ii}`. Regularization is a numerical modeling
choice and must be reported separately from the inviscid field equation.

### Local sparse solve option

A global dense solve is probably inappropriate for live wake updates. The
theory still supports local solves if each particle `k` forms a neighbor set
`N(k)` from kernel support or a cutoff radius. For each neighborhood:

```math
M^{(k)}_{pq}=\zeta_q(x_p),
\qquad p,q\in N(k),
```

solve for neighborhood `\dot{\Gamma}` and retain the center particle's value.
Overlapping neighborhoods can be blended by partition-of-unity weights or by
accepting the owner-center value. This is an approximation to the projected
field equation, not an exact global solve, so the offline diagnostic should
measure the actual residual reduction before any live implementation.

### Projected matrix-free solve path

The preferred practical path is a matrix-free Krylov solve, avoiding explicit
formation of the overlap matrix. For a trial coefficient-rate vector `v`, the
matrix-vector product is a scalar kernel sum applied independently to the three
vector components:

```math
(Mv)_k = \sum_i M_{ki}v_i.
```

For center collocation,

```math
M_{ki} = \zeta_{\sigma_i}(X_k-X_i).
```

For Gaussian-kernel Galerkin with normalized Gaussian filters,

```math
M_{ki}
= \int G_{\sigma_k}(x-X_k)G_{\sigma_i}(x-X_i)\,dx
= G_{\sqrt{\sigma_i^2+\sigma_k^2}}(X_k-X_i).
```

Thus the Galerkin mat-vec is also a Gaussian kernel sum, but with a pairwise
effective width. With fixed `\sigma`, this reduces to the ordinary Gaussian
sum with width `\sqrt{2}\sigma`; with slowly varying `\sigma`, particles can be
binned by width or evaluated with a variable-bandwidth near-neighbor kernel.

Two acceleration paths should be compared:

- **Fast Gauss Transform / FMM-style mat-vec:** evaluate the Gaussian kernel
  sums without forming `M`. This is most useful if kernel tails matter enough
  that hard truncation changes the residual or the thrust history.
- **Nearest-neighbor truncated mat-vec:** build a cell-list/hash-grid neighbor
  set and evaluate only pairs inside `c\sigma_\mathrm{eff}`. This is likely
  simpler and can be `O(N)` if the number of meaningful neighbors per particle
  remains bounded by wake maintenance, merging, and core-size choices.

Galerkin is the preferred solve formulation because the Gaussian-overlap matrix
is symmetric positive definite when the active kernels are numerically
independent, making matrix-free CG appropriate. The initial guess should be the
traditional single-particle FLOWVPM estimate already used by the Euler update,
using the same velocity-gradient convention and `pfield.transposed` branch as
FLOWVPM's active `_euler` path:

```math
\dot{\Gamma}^{(0)}_i = \mathcal{S}_\mathrm{FLOWVPM}(J_i,\Gamma_i).
```

This turns the overlap solve into a correction to the existing update rather
than a wholesale replacement. It also gives a useful convergence diagnostic:
if the initial Galerkin residual

```math
\frac{\|M_G\dot{\Gamma}^{(0)}-b_G\|}
     {\|b_G\|+\epsilon}
```

is already small, the overlap system adds little; if CG converges in a few
iterations and materially reduces this same Galerkin residual, overlap is
probably an important missing term. Center-collocation residuals may still be
reported as pointwise diagnostics, but they are not the residual minimized by a
Galerkin solve and should not be used as its acceptance metric.

Collocation remains useful as a cheap diagnostic mat-vec because it uses the
same point-sampled data as the saved particle state. However, the collocation
matrix is generally nonsymmetric, so plain CG is not guaranteed valid unless the
kernel/discretization is symmetrized or the normal equations are solved. If the
collocation path is kept, use it first to test residual scales and conditioning;
then switch to Galerkin for the CG-based production candidate.

### Immediate monitor tasks

Before attempting a live overlap solve, add a diagnostic monitor that evaluates
the overlap matrix-vector product on the particle wake during time marching.
This monitor should start with the center-collocation product because it maps
directly onto saved particle-center data:

```math
(M v)_k = \sum_i \zeta_{\sigma_i}(X_k-X_i)v_i.
```

The first useful vectors are the physical single-particle stretching estimate
and the traditional single-particle FLOWVPM strength-rate estimate. The physical
operator follows the vorticity equation and the stated Jacobian convention:

```math
v^\mathrm{phys}_i = J_i\Gamma_i.
```

The FLOWVPM-consistency operator must be computed exactly as the active FLOWVPM
Euler branch would compute its stretching term:

```math
v^\mathrm{fvpm}_i = \mathcal{S}_\mathrm{FLOWVPM}(J_i,\Gamma_i).
```

For FLOWVPM's default `pfield.transposed=true`, this is the transposed branch in
the stored `velocity_gradient` layout, not the physical `J_i\Gamma_i` operator.
These two operators answer different questions and must not be mixed in the same
residual.

The monitor should not compare `Mv` directly to `v`, because `v` is a
coefficient-rate vector while `Mv` is a sampled/projected field-rate vector. The
first valid reduced diagnostics are field-space residuals of isolated
coefficient-rate estimates against reduced overlap RHS values. Report the
physical-vorticity residual

```math
b^{\mathrm{red},\mathrm{phys}}_k = J_k\omega(X_k),
\qquad
r^{\mathrm{red},\mathrm{phys}}
= Mv^\mathrm{phys}-b^{\mathrm{red},\mathrm{phys}},
\qquad
\frac{\|r^{\mathrm{red},\mathrm{phys}}\|}
     {\|b^{\mathrm{red},\mathrm{phys}}\|+\epsilon},
```

and the FLOWVPM-consistency residual

```math
b^{\mathrm{red},\mathrm{fvpm}}_k
= \mathcal{S}_\mathrm{FLOWVPM}(J_k,\omega(X_k)),
\qquad
r^{\mathrm{red},\mathrm{fvpm}}
= Mv^\mathrm{fvpm}-b^{\mathrm{red},\mathrm{fvpm}},
\qquad
\frac{\|r^{\mathrm{red},\mathrm{fvpm}}\|}
     {\|b^{\mathrm{red},\mathrm{fvpm}}\|+\epsilon}.
```

Here both reduced RHS choices intentionally use the local-velocity/equal-basis
approximation. Neither is enough by itself for the go/no-go decision because both
omit the convective-basis correction. The monitor must also estimate that term
by sampling the same neighbor set and truncation radius used by the collocation
mat-vec:

```math
c_k =
\sum_i \Gamma_i
  \left(u(X_k)-u_i\right)\cdot\nabla\zeta_i(X_k),
\qquad
b^{\mathrm{sampled},\mathrm{phys}}_k = J_k\omega(X_k) - c_k,
\qquad
b^{\mathrm{sampled},\mathrm{fvpm}}_k
= \mathcal{S}_\mathrm{FLOWVPM}(J_k,\omega(X_k)) - c_k,
```

and report

```math
r^{\mathrm{sampled},\mathrm{phys}}
= Mv^\mathrm{phys}-b^{\mathrm{sampled},\mathrm{phys}},
\qquad
r^{\mathrm{sampled},\mathrm{fvpm}}
= Mv^\mathrm{fvpm}-b^{\mathrm{sampled},\mathrm{fvpm}},
\qquad
\frac{\|r^{\mathrm{sampled},\mathrm{phys}}\|}
     {\|b^{\mathrm{sampled},\mathrm{phys}}\|+\epsilon},
\qquad
\frac{\|r^{\mathrm{sampled},\mathrm{fvpm}}\|}
     {\|b^{\mathrm{sampled},\mathrm{fvpm}}\|+\epsilon}.
```

The difference between the reduced and sampled residuals is itself a diagnostic
for whether local velocity variation across overlapping kernels is important.
The difference between the physical and FLOWVPM-consistency residuals is a
separate diagnostic for how much FLOWVPM's transposed/reformulated update differs
from the physical vorticity-equation operator. Do not interpret that difference
as particle-overlap error.
A row-normalized or mass-lumped quick check may also be useful, but it must be
labeled as an approximation; otherwise the no-overlap limit would incorrectly
report a normalization discrepancy because
`M_{ii}v_i = \zeta_i(X_i)v_i`, not `v_i`. Store per-step aggregate diagnostics
such as absolute and relative residual norms, particle count, neighbor count,
and timing metadata. If storage cost is acceptable, also write a per-particle
residual scalar so the discrepancy can be localized in ParaView. This diagnostic
does not change the particle update; its purpose is to answer whether the
overlap equation would materially alter the current FLOWVPM behavior on a
chosen case.

Implementation tasks:

1. Add a monitor for `PanelParticleWake` that reads particle `gamma`,
   `sigma`, `velocity_gradient`, and positions after influence evaluation has
   populated the particle velocity gradient.
2. Implement the collocation mat-vec efficiently. Prototype at least a direct
   reference path and one scalable path, likely a cell-list/hash-grid
   nearest-neighbor truncation. Add FGT/FMM-style Gaussian summation only if
   truncation error is too large or neighbor counts grow too high.
3. Test mat-vec accuracy against the direct reference on controlled particle
   clouds, including fixed and varying `sigma`, near-coincident particles, and
   wake-like nonuniform spacing. Report truncation radius, error, neighbor
   count, and timing.
4. During the selected rotor-hover case, record per-step physical and
   FLOWVPM-consistency residuals for both reduced and sampled RHS choices, not
   `\|Mv-v\|`. Report reduced-vs-sampled differences so the omitted
   convective-basis term cannot be mistaken for an overlap-only effect, and
   physical-vs-FLOWVPM differences so the active transposed/reformulated update
   cannot be mistaken for overlap error.
5. Use the sampled residuals, not the reduced residuals alone, as the go/no-go
   indicators for spending time on the full Galerkin CG solve. The physical
   residual tests the vorticity equation; the FLOWVPM-consistency residual tests
   whether an overlap solve would materially change the current live update.
6. If the collocation monitor shows a material discrepancy at acceptable cost,
   extend the mat-vec machinery to the Gaussian-overlap Galerkin product and
   use the same monitor outputs before adding a solver.
7. Only after the reduced and sampled residuals are understood, add a full-RHS
   monitor mode that includes candidate SFS sources and sigma-evolution
   consistency terms. This mode must report which viscosity and sigma-evolution
   mechanisms are active and must not count fixed-core diffusion and core
   spreading as independent RHS sources.

### FLOWVPM/FLOWPanel mapping for the next phase

The generic system maps cleanly onto saved `PanelParticleWake` state:

- `X_i`: particle point coordinates in the VTP file;
- `\Gamma_i`: `gamma`;
- `\sigma_i`: `sigma`;
- `u_i`: `velocity`;
- `\omega_i`: `vorticity`, useful as a check against
  `\sum_j\Gamma_j\zeta_j(X_i)`;
- `J_i`: `velocity_gradient`; the isolated FLOWVPM estimate must follow the
  active `pfield.transposed` branch in `_euler`, not a hard-coded `J_i\Gamma_i`
  convention. FLOWVPM's default `transposed=true` uses the transposed branch in
  the stored layout, while `transposed=false` uses the classic `J\Gamma`
  branch;
- `S_i`: candidate SFS strength-rate source derived from saved `SFS`, `C`,
  `sigma`, and `zeta0`; raw saved `SFS` is not itself `S_i`. In the ClassicVPM
  Euler update the direct coefficient-rate contribution is
  `-C * SFS * sigma^3 / zeta0`; ReformulatedVPM also includes its
  formulation-dependent `Z` coupling when matching the live update;
- `\dot{\sigma}_i` or `\nu_{\mathrm{eff},i}`: required only for a full RHS
  diagnostic that attempts to include basis-width evolution. This must include
  every active sigma-changing mechanism, not just viscous core spreading:
  continuous `CoreSpreading`, ReformulatedVPM's `Z` update, and discrete
  core-reset/RBF strength recomputation events;
- optional model data: `C`, `vol`, and `circulation`.

The local code path writes and reloads these fields for `PanelParticleWake`
states, and live propagation currently calls `FLOWVPM._euler(pfield, dt;
relax=particle_relax)`. One implementation caveat is that the `ParticleField`
default integration flag is `rungekutta3`, while this FLOWPanel path calls the
private `_euler` routine directly. A full-RHS diagnostic must therefore verify
the actual live `sigma` evolution before assuming `\dot{\sigma}=\nu/\sigma`;
either measure `(\sigma(t+\Delta t)-\sigma(t))/\Delta t` from the wake state used
by the diagnostic or explicitly force and document the intended integration
convention. For ReformulatedVPM, the formulation-dependent `Z` coupling changes
both `\Gamma` and `\sigma`, so matching the live update requires including both
effects consistently.
Therefore the next concrete step should be a monitor and offline diagnostic, not
a live update:

1. Build the collocation mat-vec diagnostic monitor described above.
2. Load a saved particle VTP or restart state, or run the monitor on the chosen
   case once selected.
3. Recompute `\omega(X_k)` from the particle basis and compare it to the saved
   `vorticity` field.
4. Implement a matrix-free overlap mat-vec, first with nearest-neighbor
   truncation and then, if needed, with a Fast Gauss Transform / FMM-style
   Gaussian sum.
5. Form both `b^{\mathrm{red},\mathrm{phys}}_k = J_k\omega(X_k)` and
   `b^{\mathrm{red},\mathrm{fvpm}}_k =
   \mathcal{S}_\mathrm{FLOWVPM}(J_k,\omega(X_k))` first, intentionally invoking
   the local-velocity/equal-basis approximation and omitting SFS,
   `\dot{\sigma}`, viscosity, and convective-basis corrections.
6. Also form the sampled convective-basis estimate
   `c_k=\sum_i\Gamma_i(u(X_k)-u_i)\cdot\nabla\zeta_i(X_k)` over the same
   neighbor set as the collocation mat-vec, and compare against
   both `b^{\mathrm{sampled},\mathrm{phys}}_k=J_k\omega(X_k)-c_k` and
   `b^{\mathrm{sampled},\mathrm{fvpm}}_k=
   \mathcal{S}_\mathrm{FLOWVPM}(J_k,\omega(X_k))-c_k`.
7. Compare the monitored field-space residuals of isolated physical and FLOWVPM
   updates against both RHS choices; use the sampled residuals as the go/no-go
   indicators for whether a solver prototype is justified.
8. If justified, solve the Galerkin system with CG, initialized by the isolated
   FLOWVPM update
   `\dot{\Gamma}^{(0)}_k = \mathcal{S}_\mathrm{FLOWVPM}(J_k,\Gamma_k)`; keep
   center collocation as a cheaper diagnostic path, not the preferred CG
   target.
9. Measure the Galerkin residual
   `\|M_G\dot{\Gamma}-b_G\|/(\|b_G\|+\epsilon)` before and after the overlap
   solve. The point-sampled center-collocation residual may also be reported,
   but only as an auxiliary diagnostic.
10. Add SFS and sigma-evolution terms only in a later full-RHS diagnostic, after
   verifying the saved `SFS` field's sign/scaling and all active `\sigma`
   evolution: viscous core spreading, ReformulatedVPM `Z` coupling, and any
   core-reset/RBF recomputation events.

Only after the fixed-core reduced and sampled residuals are understood should
the diagnostic add SFS terms, sigma-evolution terms, or compare directly to
Pedrizzetti/corrected-Pedrizzetti relaxation.

One FLOWVPM-specific caveat: `CoreSpreading` is not only continuous
`\sigma`-growth. When the core-size growth threshold is reached, FLOWVPM resets
the core sizes and runs an RBF/CG strength recomputation to preserve the
pre-reset vorticity field. That reset solve is already an overlap-style
coefficient correction, so it should be treated as a baseline and possible
reuse target before adding a separate live Galerkin solve.

### Initial theory verdict

The initial theory phase supports a defensible `M\dot{\Gamma}=b` formulation:
the matrix is an overlap/projection matrix for scalar kernels, and the right
hand side is the projected filtered vorticity equation. The result is not yet a
validated replacement for relaxation. Its near-term value is as a residual
diagnostic: if the overlap solve does not materially reduce the saved-state
field residual, or if local matrices are routinely ill-conditioned in the
stable-wake rotor cases, item 008 should pivot to a negative result rather than
live wake changes.

## 2026-06-24 Collocation mat-vec residual diagnostic (offline)

Built the collocation mat-vec residual diagnostic prescribed by the "Immediate
monitor tasks" / "next concrete step" as a **standalone, read-only offline tool**
over saved particle VTP states (no simulation re-run, no package `src/` changes):

- `examples/particle_overlap_residual.jl` — reusable, side-effect-free helper
  module (kernel `ζ_σ` + `∇ζ_σ`, flat CSR cell-list neighbor search, fused
  collocation mat-vec, the physical and FLOWVPM `_euler`-transposed stretching
  operators, basis-vorticity recompute, reduced/sampled RHS, four residuals,
  conditioning proxies). Structured so a future live `AbstractMonitor` can reuse
  the same core.
- `examples/particle_overlap_residual_diag.jl` — driver: loads VTP states, runs
  validation gates, writes `data/<run>/particle_overlap_residual.csv`, prints a
  go/no-go summary.

### Conventions / trust of saved fields (verified, drives the design)

- Step ordering in `src/FLOWPanel_simulate.jl`: solve/influence populates particle
  `J`/`U` → `write_vtk` → `propagate!` (`_euler` strength update + merge). So the
  VTK write **precedes** the strength update — saved `J`/`U` are evaluated at
  exactly the saved `(X,Γ,σ)`: a consistent snapshot, no staleness.
- Run logs confirm `BODY_HESSIAN_TO_PARTICLES=false`,
  `PANEL_WAKE_HESSIAN_TO_PARTICLES=false`, `PARTICLE_HESSIAN_SELF=true` ⇒ saved
  `J` is the **clean particle-only induced gradient with the self/overlap term
  included** — exactly 008's operator, and the *same* `J` the live `_euler`
  consumed (so `S_FLOWVPM(J_i,Γ_i)` from saved `J` reproduces the live update
  exactly).
- Operator layout pinned against `FLOWVPM_timeintegration.jl` `_euler` (lines
  66–76): with column-major reshape `M[a,b]=∂u_a/∂x_b`, the **classic** branch =
  physical `(ω·∇)u`, the **transposed** branch (default `transposed=true`) = the
  live `S_FLOWVPM`. `get_W1=J[6]−J[8]=∂u₃/∂x₂−∂u₂/∂x₃` cross-checks the layout.
  Diagnostic inputs are the trusted primitives `(X,Γ,σ)` + saved `J`; saved
  `velocity` enters the convective term.

### Validation gates (all pass)

1. **No-overlap limit**: isolated particles ⇒ `Γ̇_i → J_iΓ_i` recovered;
   `r_red ≈ 1e-16`.
2. **Mat-vec vs brute-force O(N²)** on a real 4000-particle subsample: rel diff
   `4.1e-6` = pure Gaussian-tail truncation at the `4·σmax` cutoff (negligible vs
   residual signals).
3. **Kernel-normalization calibration**: basis `ω=ΣΓ_jζ_j` vs curl-of-`J`
   (`∇×u_reg=ΣΓ_jζ_j` in regularized VPM) agrees to **0.186** — confirms ζ and the
   `σ³` normalization. NOTE: the saved `vorticity` point-field is **empty (all
   zeros)** in this run, so curl-of-`J` is the calibration reference, not the
   saved field.

### Result on `data/rotor_hover_pressure_comparison` (20 settled states, 340–359)

Highly stable and reproducible across the settled hover window (per-state CSV in
`data/rotor_hover_pressure_comparison/particle_overlap_residual.csv`):

| Quantity | mean (340–359) | note |
| --- | --- | --- |
| basis ω vs curl-of-J | 0.186 | calibration OK |
| neighbors / particle | ~22,500 | heavy overlap (σ_med≈0.020, σ_max≈0.052, np≈39,400) |
| diag-dominance (median) | **0.006** | M strongly NON-diagonally-dominant |
| convective-correction rel size | 0.443 | non-negligible |
| `r_reduced,phys` | 0.315 | |
| `r_reduced,fvpm` | 0.453 | |
| **`r_sampled,phys`** | **0.538** | physical vorticity-equation residual |
| **`r_sampled,fvpm`** | **0.621** | go/no-go: FLOWVPM-consistency residual |

Per-state cost ~13 s (fused single neighbor pass over a flat CSR cell list).

### Interpretation / verdict: GO (materially large residual)

- The live FLOWVPM single-particle strength-rate update sits ~**62%** (relative
  field-space norm), stably, from the overlap-projected filtered-vorticity RHS in
  the settled rotor-hover wake. By the item's own go/no-go rule (the **sampled**
  FLOWVPM-consistency residual), this is materially large ⇒ an overlap solve
  **would** materially change the current live update.
- The gap is not an artifact of omitting the convective-basis term: adding it
  *increases* the residual (reduced→sampled fvpm 0.453→0.621), so the omitted
  convective term does not explain the isolated-update discrepancy.
- The physical-vs-FLOWVPM difference (reduced 0.315 vs 0.453) is the separate
  transposed-reformulation effect the item warned about, **not** overlap error.
- Conditioning caveat: this regime is **heavy overlap** with median diagonal
  dominance 0.006 — collocation `M` is far from diagonally dominant, exactly the
  predicted failure mode. This is why the **next phase must use the Gaussian-
  overlap Galerkin (SPD) Gram matrix with matrix-free CG**, not collocation, as
  the production solve; the collocation residual here is the cheap pointwise
  diagnostic, and it has cleared the go/no-go bar.

### Review note: particle-induced-gradient scope

Clear-context review noted that the saved `velocity_gradient` used here does
**not** include body-panel or panel-wake-row Hessian contributions for the
selected run: `BODY_HESSIAN_TO_PARTICLES=false` and
`PANEL_WAKE_HESSIAN_TO_PARTICLES=false`, while
`PARTICLE_HESSIAN_SELF=true`. This is not an implementation bug for item 008.
The item is intentionally scoped to the particle-induced filtered-vorticity
operator and to whether an overlap solve would materially change the current
particle-only-gradient FLOWVPM strength update.

Consequently, the `physical` residuals in this section should be read as
physical with respect to the particle-induced velocity-gradient operator, not as
a complete hybrid panel+particle vorticity-equation residual. A complete
hybrid-flow residual would also contain a panel/body strain term
`J_panel \omega_particle`, because curl-free potential flow can still stretch or
tilt existing particle vorticity. That extension is out of scope for this item
unless the research question is widened beyond particle-induced overlap.

The convective-basis correction is different: saved particle velocity includes
the velocity that actually advects the particles, including enabled panel/body
velocity contributions. Using that total advecting velocity in
`(u(X_k)-u_i)\cdot\nabla\zeta_i` is consistent with the Lagrangian basis-motion
term, even while `J` is kept particle-induced only. This means the conclusion
that the convective-basis term is non-negligible and increases the sampled
residual remains valid for the implemented particle-update diagnostic; it should
not be reinterpreted as evidence about the omitted panel/body strain term.

### Next phase (justified, not yet done)

Do **not** jump directly to a live wake update. The next step should be a
Gaussian-overlap **Galerkin residual prototype** using the same particle-induced
operator scope as the current diagnostic.

Recommended sequence:

1. Implement a matrix-free Galerkin mat-vec
   `M_G v`, with
   `M^G_{ki}=\int \zeta_k(x)\zeta_i(x)\,dx`. For normalized Gaussian kernels this
   is another Gaussian sum with pairwise width
   `\sqrt{\sigma_k^2+\sigma_i^2}`.
2. Build the Galerkin RHS with the current item scope: particle-induced
   `J`, particle-basis `\omega`, and the convective-basis correction if feasible.
   Keep panel/body strain out of scope for this item, and keep SFS / `\dot{\sigma}`
   terms out until the fixed-core overlap residual is understood.
3. Measure the initial Galerkin residual of the current isolated FLOWVPM update,
   initialized as
   `\dot{\Gamma}^{(0)}_i=\mathcal{S}_\mathrm{FLOWVPM}(J_i,\Gamma_i)`:

   ```math
   \frac{\|M_G\dot{\Gamma}^{(0)}-b_G\|}{\|b_G\|+\epsilon}.
   ```

4. Run matrix-free CG, initialized from `\dot{\Gamma}^{(0)}`, and measure both
   the post-solve Galerkin residual and the correction size:

   ```math
   \frac{\|M_G\dot{\Gamma}^{CG}-b_G\|}{\|b_G\|+\epsilon},
   \qquad
   \frac{\|\dot{\Gamma}^{CG}-\dot{\Gamma}^{(0)}\|}
        {\|\dot{\Gamma}^{(0)}\|+\epsilon}.
   ```

5. Treat this as the next go/no-go gate. If CG substantially reduces the
   Galerkin residual with a coherent, not noise-dominated correction, then a live
   opt-in overlap update is worth prototyping. If the Galerkin residual is small,
   CG is ill-conditioned, or the correction is dominated by localized/noisy
   modes, item 008 should pivot to explaining why the collocation signal does not
   justify a live solve.
6. Only after the Galerkin residual and correction are understood, add SFS /
   `\dot{\sigma}` full-RHS terms and compare against Pedrizzetti or
   corrected-Pedrizzetti relaxation on the stable-wake cases. A second target run
   will be supplied by the user. The diagnostic helper is structured for reuse as
   a live monitor if desired.
