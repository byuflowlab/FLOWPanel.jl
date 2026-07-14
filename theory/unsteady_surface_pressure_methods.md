# Unsteady Surface-Pressure Prediction Methods

Date: 2026-07-14

This note compares pressure-producing methods for unsteady FLOWPanel
simulations. It is a design and theory document only. It does not propose code
changes in this task, and it does not claim that any unvalidated method is a
reliable default.

The derivations reuse the ALE notation and monitor contract in
[`theory/unsteady_pressure_monitors_2026-07-11.md`](unsteady_pressure_monitors_2026-07-11.md),
with supporting details in [`docs/unsteady_bernoulli.md`](../docs/unsteady_bernoulli.md),
[`docs/pressure_poisson.md`](../docs/pressure_poisson.md), and
[`docs/kernel_gradients.md`](../docs/kernel_gradients.md).

## Scope

The eight surface-pressure candidates are:

1. history-based unsteady Bernoulli;
2. `PressureLaplace` with corrected Hessian acceleration;
3. `PressureLaplace` with corrected edge-difference acceleration;
4. `PressureLaplace` with surface-velocity reconstruction;
5. tangent-linear discrete Bernoulli;
6. Kirchhoff-potential decomposition;
7. surface pressure-gradient least squares;
8. thin-volume Euler pressure recovery.

Lamb-vector pressure is not counted as a candidate method here. The current
`PressureLaplace(acceleration_form=:lamb_vector)` path is diagnostic only
because a complete ALE surface derivation with moving grids, exterior doublet
jumps, and wake vorticity has not been closed. Impulse, Kutta-Joukowski, and
virtual-work methods are likewise not pressure methods: they can produce useful
integrated-force checks, but they do not uniquely recover panel pressure.

## Common Notation

Let `x` be an inertial position and `t` time. The body control-point grid moves
with velocity

```math
w(x,t) = V_O(t) + \Omega(t) \times (x - x_O(t))
```

for a rigid frame, or the corresponding prescribed grid velocity for a
deforming surface. The inertial fluid velocity is `u`, the grid-relative
velocity is

```math
q = u - w,
```

and the panel normal is `n`. The tangent projection is

```math
P = I - n n^T.
```

The reconstructed inertial exterior surface trace used by current monitors is

```math
u_s = w + Pq.
```

In FLOWPanel, after `kinematic_velocity!`, `body.velocity` stores `q`,
`body.velocity_kinematic` stores `w`, and `body.velocity_gradient[:,:,i]` stores

```math
G_{kl,i} = \frac{\partial u_k}{\partial x_l}
```

for induced inertial velocity at panel `i`, when a monitor requests Hessians
through `monitor_requires_body_hessian`.

Panel unknowns are the columns of `body.strength`. For
`ConstantDoublet`/`VortexRing` lifting bodies, the relevant scalar potential
jump or circulation column is usually named `"mu"` or `"gamma"` in helper code.
FLOWPanel's sign convention follows its element kernels, which include explicit
sign handling relative to GeometricTools normals. The exterior doublet trace
contains a half-jump contribution; current corrected-Hessian pressure logic
adds a postprocessed exterior `-1/2 grad_s(mu_code)` velocity-jump derivative
through `compute_mu_gradient!` and `compute_surface_velocity_gradient!`.

Let the body solve be written abstractly as

```math
R(\mu, b, z_w, t) = A(b, z_w, t) \mu - r(b, z_w, t) = 0,
```

where `b` denotes body geometry and kinematics, `z_w` wake state, and `r`
contains prescribed flow, source terms, Dirichlet/Neumann data, and wake
influence. A pressure gauge is always required. Gauge choices used here are
either a reference Bernoulli constant `C(t)` or a pinned panel pressure
`p_k = p_ref`.

Pressure force integration is the usual panel sum

```math
F = -\sum_i p_i n_i A_i,
\qquad
M = -\sum_i (x_i - x_ref) \times p_i n_i A_i,
```

with the same sign convention as `ForceMonitor`, which reports force on the
body.

## Method 1: History-Based Unsteady Bernoulli

### Mathematical Derivation

For scalar-potential, incompressible, inviscid flow,

```math
u = \nabla \phi,
```

Euler's equation integrates to

```math
p + \rho \partial_t \phi + \frac{\rho}{2}|u|^2 = C(t).
```

At moving panel control points,

```math
\frac{d}{dt}\phi(x_i(t),t)
= \partial_t \phi_i + w_i \cdot \nabla \phi_i,
```

so

```math
\partial_t \phi_i = D_g \phi_i - w_i \cdot u_i.
```

FLOWPanel evaluates pressure relative to a freestream gauge:

```math
p_i =
\frac{\rho}{2}(|U_\infty|^2 - |u_{s,i}|^2)
- \rho \partial_t \phi_i.
```

Assumptions: a single-valued scalar potential exists for every velocity
contribution entering `|u_s|^2`; density is constant; viscous stresses are
ignored. Particle wakes or vortex-only sources break the complete scalar
potential assumption.

### Practical FLOWPanel Implementation

Current behavior is `PressureBernoulli(rho; unsteady=true)`. The monitor owns
pressure arrays, evaluates scalar potential at body control points with a
`FastMultipole.ProbeSystem`, adds freestream potential `uinf dot x`, and
finite-differences panel-following potential history. Startup is safe:
first sample seeds history with zero temporal contribution, the second sample
uses backward Euler, and later samples use variable-step BDF2.

Required fields are `body.controlpoints`, `body.velocity`,
`body.velocity_kinematic`, `body.normals`, scalar-potential-capable body/wake
sources, `uinf`, `i_step`, and runtime `dt`. Monitor integration is simple:
`PressureBernoulli` provides `:P`; `ForceMonitor(source=:pressure)` requires
`:P`. Direct and FMM backends are both supported through `influence!`, but only
sources with scalar potential are included in the potential history. When
vector-potential-only wake sources are excluded, the monitor warns and
continues with a mixed partial diagnostic.

Pseudocode:

```text
for each body:
    u_s = body.velocity_kinematic + tangent(body.velocity)
    phi = scalar_potential(body + scalar wake sources, controlpoints)
          + dot(uinf, controlpoint)
    Dg_phi = zero/BE/BDF2 panel-following derivative(phi)
    phi_t = Dg_phi - dot(w, u_s)
    p = 0.5*rho*(|uinf|^2 - |u_s|^2) - rho*phi_t
```

Cost is one scalar-potential influence evaluation per pressure step plus
`O(N)` history and pressure work. Direct evaluation is `O(N^2)` in the number
of source/target panels; FMM is approximately `O(N log N)` with the current
backend.

### Advantages, Disadvantages, and Failure Modes

Advantages: it is the exact scalar-potential pressure relation when the wake
and body are potential-flow sources with scalar potential; it produces panel
pressure directly; it is easy to compare against force integration and current
examples.

Disadvantages: it is history-sensitive, has one-step startup order loss, and
can show phase/amplitude error when potential history is underresolved. It is
not complete for particle wakes or vector-potential-only contributions. It is
also sensitive to gauge consistency and exterior doublet trace handling.

Failure modes include false startup spikes if history is not seeded, mixed
frame kinetic energy if `q` is used instead of `u_s`, missing wake-potential
history for particle wakes, and panel-index history invalidation after topology
or panel-count changes.

## Method 2: Pressure Laplace With Corrected Hessian

### Mathematical Derivation

Euler's equation gives

```math
\nabla_s p = -\rho P a,
\qquad
a = \partial_t u + (u \cdot \nabla)u.
```

At moving control points,

```math
D_g u_s = \partial_t u_s + (w \cdot \nabla)u_s,
```

so the ALE acceleration is

```math
a = D_g u_s + G q_t,
\qquad q_t = Pq.
```

The surface Poisson form is

```math
-\Delta_s p =
\rho \nabla_s \cdot(Pa),
```

with a pinned reference pressure. The corrected-Hessian mode uses analytic
`G = grad(u)` from kernel Hessians, plus an explicit surface derivative of the
exterior doublet jump. Projection is applied once to the final acceleration,
not as `PGP`.

### Practical FLOWPanel Implementation

Current behavior is
`PressureLaplace(bodies, rho; unsteady, gradient_mode=:corrected_hessian,
acceleration_form=:material_derivative)`. Construction preallocates sparse
surface-Laplace data, one RHS and pressure vector per body, CG workspaces, and
scratch arrays. The monitor requires Hessian accumulation through
`monitor_requires_body_hessian(::PressureLaplace)`.

Required fields include `body.velocity`, `body.velocity_kinematic`,
`body.velocity_gradient`, `body.normals`, `body.controlpoints`, `body.nodes`,
`body.cells`, `body.neighbor`, doublet/circulation strength if available,
and shared-edge topology. The sparse operator is built from interior shared
edges using the conormal weight returned by `_pressure_edge_conormal_weight`.
`reference_panel` is pinned to `reference_pressure`.

Direct/FMM implication: direct and FMM must agree on `velocity_gradient` layout
and sign. FMM support depends on Hessian-capable source kernels and current
`FastMultipoleBackend` accuracy; direct support is the verification baseline.
Cost is one velocity-gradient influence path during the normal solve plus a
sparse CG solve. Storage is `O(N + E)`.

Startup behavior follows the same velocity-history policy as the ALE note:
steady mode has no temporal derivative; unsteady mode seeds with zero, then BE,
then BDF2.

### Advantages, Disadvantages, and Failure Modes

Advantages: it avoids scalar-potential history and can, in principle, handle
rotational velocity contributions if their velocity gradient is available. It
produces a full pressure distribution and exposes direct/FMM gradient defects.

Disadvantages: it is sensitive to Hessian signs, tensor layout, self limits,
surface jump correction, and nonsymmetric velocity gradients. It cannot repair
an inaccurate acceleration field. Conditioning depends on surface mesh quality
and the pinned gauge.

Failure modes include missing exterior jump derivative, incorrectly adding a
rigid `[Omega]_x` term to `body.velocity_gradient`, blanket `PGP` projection,
FMM Hessian/direct mismatch, nonsymmetric-gradient contraction mistakes, and
poor CG convergence on nonorthogonal or degenerate panel adjacency.

## Method 3: Pressure Laplace With Corrected Edge Difference

### Mathematical Derivation

This mode keeps the same surface Poisson equation but avoids panel-centered
Hessians for the convective term. Across an edge between panels `i` and `j`,
let

```math
r_{ij}=x_j-x_i,\qquad
q_e=\frac{1}{2}(P_i q_i + P_j q_j).
```

The directional convective increment must respect nonsymmetric gradients:

```math
r^T G q = q^T G r + r \cdot (\omega \times q),
```

where `omega` is exterior volumetric vorticity. The edge flux combines a
midpoint unsteady increment and the corrected directional convective increment.
The pressure gauge is again a pinned panel.

### Practical FLOWPanel Implementation

Current behavior is
`PressureLaplace(...; gradient_mode=:edge_difference,
acceleration_form=:material_derivative)`. It is the warned compatibility
fallback when `gradient_mode` is omitted, but the warning explicitly says no
formulation has passed the complete unsteady gate.

Required fields are `u_s` history, tangent-projected `body.velocity`, shared
edges, and `body.induced_vorticity` with bound-sheet kappa removed for the
edge identity. Monitor traits request induced vorticity for this mode through
`monitor_requires_induced_vorticity`. It is Hessian-free, which makes it
attractive for FMM/particle paths where Hessian support is incomplete or
expensive.

Pseudocode:

```text
for each shared edge i-j:
    udot_edge_dot_r = midpoint(velocity_dot) dot (x_j - x_i)
    q_edge = midpoint(tangent(body.velocity))
    du = u_s[j] - u_s[i]
    omega = midpoint(exterior_vorticity_without_bound_kappa)
    conv = dot(q_edge, du) + dot(cross(omega, q_edge), x_j - x_i)
    flux = rho * w_ij * (udot_edge_dot_r + conv)
    b[i] += flux
    b[j] -= flux
pin reference pressure and solve L p = b
```

Cost is `O(E)` RHS assembly plus CG. It avoids a Hessian influence request but
does require vorticity data when rotational wake components are present.

### Advantages, Disadvantages, and Failure Modes

Advantages: no analytic Hessian is needed; the flux is naturally paired with
the same two-point finite-volume pressure operator; direct/FMM differences are
reduced for derivative information.

Disadvantages: edge directional differences are lower-order and mesh-sensitive;
the vorticity correction depends on correct exterior-side vorticity and bound
surface kappa subtraction; it may smear sharp pressure features.

Failure modes include using bound-sheet kappa as exterior volumetric curl,
omitting the nonsymmetric correction, projecting edge differences in a way that
breaks the identity, and treating the fallback as validated production pressure.

## Method 4: Pressure Laplace With Surface-Velocity Reconstruction

### Mathematical Derivation

This mode reconstructs the surface gradient of the full inertial trace:

```math
G_s = \nabla_s u_s.
```

The ALE acceleration is approximated by

```math
a = D_g u_s + G_s q_t,
\qquad
\nabla_s p = -\rho P a.
```

The gauge and finite-volume surface Poisson solve match the other
`PressureLaplace` modes.

### Practical FLOWPanel Implementation

Current behavior is
`PressureLaplace(...; gradient_mode=:surface_velocity)`. It uses
`compute_surface_velocity_gradient!`, whose output layout is
`grad_u[k,l,i] = partial u_k / partial x_l`. The stencil follows
`compute_mu_gradient!`, with panel aspect-ratio masking and trailing-edge
isolation. This mode does not require Hessians but does require reliable
surface reconstruction of vector components.

Required fields are `u_s`, `body.controlpoints`, `body.normals`, `body.cells`,
`body.neighbor`, and trailing-edge information when available. Direct/FMM
implications are weak for derivative data because derivatives are reconstructed
after the influence solve, but the underlying velocity values still depend on
the selected backend.

Cost is roughly three scalar surface-gradient reconstructions plus CG:
`O(stencil work + E + CG)`, with `O(N)` scratch storage.

### Advantages, Disadvantages, and Failure Modes

Advantages: it is backend-agnostic with respect to Hessians; it uses the actual
postprocessed exterior surface velocity trace; it is useful as an independent
diagnostic against analytic Hessians.

Disadvantages: surface differentiation of noisy panel velocities can dominate
the pressure RHS; high-aspect-ratio panels and trailing-edge discontinuities
need special handling; curvature terms are only represented through the chosen
surface stencil.

Failure modes include stencil leakage across wake cuts, bad-panel fallback
masking real pressure gradients, over-smoothing nonconstant doublet strength,
and apparent agreement with a wrong method because both share the same velocity
history.

## Method 5: Tangent-Linear Discrete Bernoulli

### Mathematical Derivation

Instead of finite-differencing scalar potential over time, differentiate the
discrete panel solve at the current time. With

```math
R(\mu,t) = A(t)\mu(t) - r(t) = 0,
```

the tangent equation is

```math
\boxed{
A\dot{\mu} = -\left.\frac{dR}{dt}\right|_{\mu}
= \dot{r} - \dot{A}\mu .
}
```

The surface potential evaluation is also differentiated:

```math
\Phi_s(\mu,t) =
\Phi_body(\mu,b(t)) + \Phi_w(z_w(t)) + U_\infty(t)\cdot x_s(t),
```

```math
\dot{\Phi}_s =
\Phi_\mu \dot{\mu}
+ \left.\frac{d\Phi_body}{dt}\right|_\mu
+ \dot{\Phi}_w
+ \dot{U}_\infty\cdot x_s
+ U_\infty\cdot w_s.
```

The Eulerian partial derivative entering Bernoulli is

```math
\partial_t \phi_s = \dot{\Phi}_s - w_s \cdot u_s.
```

For a pure rigid body under common rigid motion, pairwise source-target
geometry in the body frame is invariant. Therefore the rigid body-to-body
operator `A` is invariant under the common rigid transformation. A pitching
wing with fixed body mesh can reuse the existing factorization and
differentiate only source terms, prescribed flow, relative wake geometry,
shedding strengths, and wake/body interactions. This invariance does not hold
for deformation, changing topology, or non-body-fixed wake/source geometry.

### Practical FLOWPanel Implementation

This is a proposed new method, not current behavior. It should be prototyped
at residual level:

```text
solve A*mu = r
form residual R(mu, state)
directionally perturb current state by +/- eps along known time tangent
J_state = (R(mu, state+eps*tangent) - R(mu, state-eps*tangent)) / (2eps)
solve A*mudot = -J_state
differentiate surface potential evaluation by the same centered direction
phi_t = phidot_surface - dot(w, u_s)
p = 0.5*rho*(|uinf|^2 - |u_s|^2) - rho*phi_t
```

Prototype recommendation: use residual-level centered directional differences
first. For verification, implement direct-kernel JVPs source by source. For
production, use a custom frozen-tree tangent FMM for the influence pieces whose
source/target positions and strengths have a known time tangent. Avoid generic
AD through wake convection or the entire N-body algorithm: it is likely to be
slow, fragile across tree construction, and hard to audit.

Required fields include the solved body unknowns, RHS/source construction,
prescribed flow and its time derivative, body/wake geometry tangents, shedding
strength tangents, and scalar-potential evaluation tangents. Solver integration
should sit beside the existing matrixful solve for the first target case:
single `ConstantDoublet` panel-wake pitching wing with fixed body mesh.

Cost is one extra linear solve per time step if `A` is reused, plus residual
JVP and potential JVP work. Direct prototype cost is several direct influence
evaluations per directional difference; production FMM should target one
tangent influence evaluation with the same asymptotic cost as a normal
influence pass.

Startup behavior is excellent: no finite-history startup is needed, because
`partial_t phi` is obtained at the current state.

### Advantages, Disadvantages, and Failure Modes

Advantages: it preserves exact scalar-potential Bernoulli panel by panel,
removes finite-history startup and phase error, can exploit rigid-operator
invariance and existing factorization, and directly tests discrete consistency
of the solve and pressure formula.

Disadvantages: it is restricted to scalar-potential contributions unless a
rotational pressure extension is added; it needs careful differentiation of
wake state, shedding strengths, and prescribed flow; and it introduces another
linear solve or tangent solve per step.

Failure modes include differentiating the wrong frame, treating wake convection
as body-rigid, missing the exterior potential jump, using a stale factorization
when body-relative geometry changes, and tuning finite-difference epsilon until
truncation/roundoff errors masquerade as agreement.

## Method 6: Kirchhoff-Potential Decomposition

### Mathematical Derivation

Classical Kirchhoff potentials decompose rigid-body potential flow into unit
translation and rotation solutions. For a rigid body,

```math
\phi(x,t) =
\phi_\infty
+ \sum_{\alpha=1}^3 V_\alpha \chi_\alpha(x;b)
+ \sum_{\alpha=1}^3 \Omega_\alpha \psi_\alpha(x;b)
+ \phi_c(x,t),
```

where `chi_alpha` and `psi_alpha` are unit translation/rotation potentials and
`phi_c` contains circulatory and wake contributions. The pressure separates
formally into prescribed rigid-motion kinetic terms, acceleration/added-mass
terms, and circulatory/wake terms:

```math
\partial_t \phi =
\dot{V}_\alpha \chi_\alpha
+ V_\alpha \partial_t \chi_\alpha
+ \dot{\Omega}_\alpha \psi_\alpha
+ \Omega_\alpha \partial_t \psi_\alpha
+ \partial_t \phi_c .
```

Moving-frame derivatives matter. Even if the body-frame basis potentials are
time-independent for a rigid body, their inertial trace derivative includes
frame transport and the subtraction `w dot u` when converting sampled
derivatives to Eulerian `partial_t phi`. The circulatory derivative
`partial_t phi_c` remains; Kirchhoff potentials do not remove the need to
differentiate wake/circulation history.

Gauge: each basis potential needs a consistent constant, and the final
pressure gauge must match the Bernoulli constant used for comparison.

### Practical FLOWPanel Implementation

This is a specialized proposed diagnostic and possible optimization. It would
precompute or solve six unit Neumann/Dirichlet problems for the body operator,
then combine them with current rigid velocities and accelerations. It should
reuse the same `A` and same body discretization as the main solve whenever the
body-relative mesh is unchanged.

Required fields are body geometry, frame velocities and accelerations,
basis-potential values at control points, scalar-potential evaluation for
circulatory/wake components, and wake/circulation time derivatives. Direct/FMM
support follows the underlying basis solves and potential evaluation.

Cost is upfront six basis solves per rigid geometry plus per-step linear
combinations, but any circulatory derivative still requires a history or
tangent-linear path.

### Advantages, Disadvantages, and Failure Modes

Advantages: it cleanly isolates added-mass and acceleration pressure for rigid
bodies; it gives a useful diagnostic for pitching/plunging cases; it may reduce
repeated solve work when rigid kinematics change but geometry does not.

Disadvantages: it is not a universal pressure shortcut. Wake/circulatory
pressure remains unresolved without `partial_t phi_c`; deforming geometry
invalidates the fixed basis; gauges across basis potentials can contaminate
pressure.

Failure modes include presenting added-mass-only pressure as complete pressure,
forgetting moving-frame transport terms, mixing body-frame and inertial
velocities in Bernoulli, and applying rigid basis invariance to deforming or
relative multi-body motion.

## Method 7: Surface Pressure-Gradient Least Squares

### Mathematical Derivation

Starting from the surface Euler relation,

```math
\nabla_s p = -\rho P a,
```

define an edge pressure increment predicted by acceleration:

```math
\delta p_{ij} \approx
-\rho\, \bar{a}_{t,ij}\cdot (x_j-x_i).
```

Let `B` be the oriented panel-edge incidence matrix and `W` edge weights. Solve

```math
\min_p \| W^{1/2}(B p - \delta p)\|_2^2
```

with a gauge, for example `p_k=p_ref` or `sum_i A_i p_i=0`. The normal
equations are

```math
B^T W B p = B^T W \delta p,
```

which is a graph pressure-Poisson solve. The residual

```math
r_e = (B p)_e - \delta p_e
```

measures nonintegrability, or discrete curl, of the pressure-gradient data.

### Practical FLOWPanel Implementation

This is a proposed pressure-producing cross-check. It can ingest acceleration
from corrected-Hessian `PressureLaplace` and surface-velocity reconstruction,
then solve a weighted graph least-squares pressure problem with the same gauge
as `PressureLaplace`.

Required fields are panel centers, shared-edge incidence, edge weights,
surface accelerations, and a gauge. It should report both pressure and edge
residual diagnostics. Direct/FMM support is inherited entirely from the
acceleration provider.

Pseudocode:

```text
for each shared edge e=(i,j):
    delta[e] = -rho * dot(0.5*(a_t[i]+a_t[j]), x[j]-x[i])
    B[e,i] = -1
    B[e,j] =  1
solve min ||sqrt(W)*(B*p - delta)|| with gauge
report edge residual and cycle/curl statistics
```

Cost is one sparse graph solve comparable to the surface Laplacian solve, plus
`O(E)` residual diagnostics.

### Advantages, Disadvantages, and Failure Modes

Advantages: it exposes whether the acceleration field is integrable as a
surface pressure; it gives localized residual diagnostics; it is a useful
independent pressure-producing check against the finite-volume RHS assembly.

Disadvantages: it cannot repair an inaccurate acceleration field. Least
squares distributes inconsistent edge data globally, so a smooth pressure can
hide a bad local acceleration source unless residuals are inspected.

Failure modes include choosing weights that suppress meaningful wake-cut
errors, using a gauge inconsistent with Bernoulli comparison, ignoring cycle
residuals, and treating least-squares pressure as more physical solely because
it is smoother.

## Method 8: Thin-Volume Euler Pressure Recovery

### Mathematical Derivation

Define an exterior offset shell around the body:

```math
x_{i,m} = x_i + \eta_m h_i n_i,
\qquad \eta_m > 0,
```

where `h_i` is a local panel length. Evaluate velocity and ALE acceleration at
several normal layers on the exterior side. Then recover pressure in the
narrow band from

```math
\nabla p = -\rho a
```

or a narrow-band Poisson/least-squares form:

```math
\min_p \| \nabla_h p + \rho a_h \|^2,
```

with boundary and gauge constraints. The surface pressure is the exterior
trace as `eta -> 0+`.

### Practical FLOWPanel Implementation

This is a later research path, not a near-term implementation. It requires a
thin exterior volume mesh or layered graph, velocity and acceleration
evaluation at offset targets, wake-sheet side selection, near-singular
regularization, and a boundary/gauge treatment. Offsets should scale with
local panel length, for example a convergence sequence
`eta in {0.25, 0.125, 0.0625}` times `h_i`, while staying on the exterior side
and avoiding crossing wake sheets.

Wake-sheet side selection is critical: pressure should be traced from the body
exterior side that matches the panel surface, not through the doublet jump or
across a freshly shed wake sheet. For panel wakes, side labels can follow
shedding-edge connectivity. For particle wakes, vorticity is volumetric and
the narrow-band acceleration may be more complete, but offset evaluation still
needs core-size and near-singular audits.

Direct/FMM support requires influence evaluation at many probe points, possibly
with velocity gradients or temporal velocity history. Cost is volumetric:
`O(M N)` direct for `M` shell points and `N` sources, or approximately
FMM influence cost plus a 3D sparse solve over shell points. Storage is far
larger than a surface solve.

### Advantages, Disadvantages, and Failure Modes

Advantages: it is the most natural route for rotational or particle-wake
completeness because it does not force all physics onto a surface scalar
potential. It can separate exterior-side pressure from surface jump artifacts
by convergence in offset distance.

Disadvantages: it is expensive, experimentally delicate, and depends on
near-singular kernel evaluation accuracy. Boundary conditions on the artificial
outer shell and cuts around wake sheets can dominate the result.

Failure modes include offsets crossing the wake sheet, evaluating on the wrong
side of a doublet jump, pressure drift from weak gauges, nonconvergent
near-singular influence, volumetric mesh anisotropy, and mistaking core-radius
regularization for physical pressure smoothing.

## Comparison Tables

### Physical Completeness and Wake Compatibility

| Method | Full surface pressure | Scalar panel wake | Particle/vector wake | Main missing physics |
|---|---:|---:|---:|---|
| History Bernoulli | Yes | Yes | Diagnostic only | scalar potential for rotational wake |
| Laplace corrected Hessian | Yes | Intended | Possible if gradients exist | validated exterior acceleration |
| Laplace edge difference | Yes | Intended | Better path with vorticity | accurate edge acceleration |
| Laplace surface velocity | Yes | Diagnostic | Diagnostic | robust surface derivatives |
| Tangent-linear Bernoulli | Yes | Best first target | No | scalar potential for rotational wake |
| Kirchhoff decomposition | Partial unless `phi_cdot` known | Added-mass diagnostic | No | circulatory/wake derivative |
| Pressure-gradient LS | Yes | Inherits acceleration | Inherits acceleration | accurate acceleration field |
| Thin-volume recovery | Yes, experimental | Possible | Best long-term path | validated offset-volume solve |

### Temporal History and Startup

| Method | Needs history | Startup sensitivity | Gauge |
|---|---:|---|---|
| History Bernoulli | Yes | zero/BE/BDF2, phase possible | Bernoulli constant |
| Laplace corrected Hessian | Only if unsteady | zero/BE/BDF2 velocity history | pinned panel |
| Laplace edge difference | Only if unsteady | zero/BE/BDF2 velocity history | pinned panel |
| Laplace surface velocity | Only if unsteady | zero/BE/BDF2 velocity history | pinned panel |
| Tangent-linear Bernoulli | No | low, current-state derivative | Bernoulli constant |
| Kirchhoff decomposition | Maybe | basis stable, `phi_cdot` still needed | basis plus final gauge |
| Pressure-gradient LS | Inherits provider | inherits acceleration provider | pinned or mean-zero |
| Thin-volume recovery | Inherits acceleration | depends on shell history/JVP | volume gauge |

### Spatial Derivatives and Exterior Jumps

| Method | Spatial derivative source | Exterior jump handling | Key audit |
|---|---|---|---|
| History Bernoulli | scalar potential and velocity | scalar-potential trace | `phi_t = Dg_phi - w dot u_s` |
| Corrected Hessian | analytic `body.velocity_gradient` plus jump gradient | explicit half-jump derivative | Hessian layout/sign |
| Edge difference | edge velocity increments plus vorticity correction | through exterior velocity increments | nonsymmetric identity |
| Surface velocity | reconstructed `grad_s u_s` | included in `u_s` stencil | TE isolation |
| Tangent-linear Bernoulli | differentiated potential evaluation | differentiated exterior trace | residual JVP |
| Kirchhoff | unit potentials and moving-frame transport | basis-gauge dependent | missing `phi_cdot` |
| Pressure-gradient LS | supplied acceleration | inherited | curl residual |
| Thin-volume | volumetric gradients/offset data | exterior-side shell | offset convergence |

### Backend, Solve Count, and Cost

| Method | Direct support | FMM support | Extra solves per step | Extra influence/eval count | Storage | Asymptotic note |
|---|---:|---:|---:|---|---|---|
| History Bernoulli | Yes | Yes | 0 | 1 scalar-potential eval | `O(N)` | direct `O(N^2)`, FMM approx `O(N log N)` |
| Corrected Hessian | Yes | Yes if Hessian kernels valid | 1 CG | Hessian requested in solve | `O(N+E)` | sparse CG plus influence |
| Edge difference | Yes | Yes | 1 CG | vorticity data | `O(N+E)` | sparse CG |
| Surface velocity | Yes | Yes | 1 CG | velocity only | `O(N+E)` | surface stencil plus CG |
| Tangent-linear Bernoulli | Prototype yes | production custom tangent FMM | 1 linear solve | residual and potential JVP | `O(N)` plus solver | reuse factorization for rigid `A` |
| Kirchhoff | Yes | Yes | 0 to many basis solves | basis potential evals | basis storage | amortized for fixed rigid body |
| Pressure-gradient LS | Yes | Yes | 1 graph LS | inherited | `O(N+E)` | graph Poisson-like |
| Thin-volume | Yes, costly | Needed for scale | 1 volume solve | many shell probes | `O(M)` | volumetric cost dominates |

### Body and Wake Cases

| Case | Best candidates | Cautions |
|---|---|---|
| Rigid pure `ConstantDoublet` panel-wake pitching wing | tangent-linear Bernoulli, history Bernoulli oracle, pressure-gradient LS check | validate direct/FMM and wake-strength derivatives |
| Deforming body | PressureLaplace modes, pressure-gradient LS | rigid operator invariance no longer applies |
| Multi-body common rigid motion | tangent-linear may reuse block operators only for invariant relative geometry | relative motion changes `A` |
| Panel wake with scalar potential | Bernoulli methods strongest | final filament/scalar-source choices affect history |
| Particle wake or rotational wake | thin-volume research path, diagnostic Laplace modes | scalar Bernoulli incomplete |
| Force-only validation | ForceMonitor, KuttaJoukowskiForce, SurfaceVorticityForce | cannot uniquely recover pressure |

## Recommendation

The first new implementation should be tangent-linear Bernoulli for the pure
`ConstantDoublet` panel-wake pitching wing. It preserves the scalar-potential
Bernoulli relation, eliminates finite-history startup and phase error, exploits
rigid body-to-body operator invariance, and yields pressure panel by panel. The
prototype should use residual-level centered directional differences; direct
kernel JVPs should verify the prototype; a frozen-tree custom tangent FMM is
the production route. Do not rely on generic AD through wake convection or the
entire N-body algorithm.

Use history Bernoulli as the validation oracle for scalar panel wakes. It is
not perfect because it has temporal-history error, but it exercises the same
physical pressure relation and current monitor integration.

Use pressure-gradient least squares as the independent pressure-producing
cross-check, with acceleration supplied by corrected-Hessian and
surface-reconstruction paths. The least-squares residual should be treated as a
first-class diagnostic, not discarded after pressure reconstruction.

Kirchhoff decomposition should be treated as a specialized added-mass
diagnostic and possible rigid-motion optimization. It is not a universal
primary method because the circulatory/wake potential derivative remains.

Thin-volume recovery should be deferred until the panel-wake problem is solved.
It is the more plausible long-term path for rotational or particle-wake
pressure completeness, but near-singular evaluation, wake-side selection, and
volumetric cost make it experimental.

Lamb-vector pressure, impulse, Kutta-Joukowski, and virtual-work methods should
not be presented as validated surface-pressure predictors. Lamb-vector pressure
currently lacks the complete ALE surface derivation. The force-only methods are
valuable integrated checks, but pressure is not uniquely recoverable from their
outputs.

## Manufactured Checks and Acceptance Gates

Discriminating checks for future implementation:

| Method | Manufactured checks |
|---|---|
| History Bernoulli | Galilean translating grid; rotating grid; scalar panel wake with known potential history; startup without `O(phi/dt)` spike |
| Corrected Hessian | nonsymmetric gradient tensor; nonconstant doublet strength; direct/FMM Hessian agreement; exterior half-jump sign |
| Edge difference | nonsymmetric-gradient identity; bound-kappa removal; nonorthogonal panel pairs; deforming wake panels |
| Surface velocity | high-aspect panels; trailing-edge stencil isolation; nonconstant `mu`; comparison with analytic surface gradients |
| Tangent-linear Bernoulli | rigid-operator invariance under common rigid motion; JVP versus centered directional difference; wake-strength differentiation |
| Kirchhoff | pure translation, pure rotation, acceleration-only added mass; moving-frame derivative sign |
| Pressure-gradient LS | nonintegrable edge pressure data; cycle residual localization; gauge invariance of pressure differences |
| Thin-volume | offset-shell convergence; wake-sheet side selection; near-singular direct/FMM agreement |

Acceptance gates before any new method is made a default:

1. surface-pressure convergence under timestep refinement and spatial
   refinement;
2. force amplitude and phase agreement against history Bernoulli for the pure
   scalar panel-wake case;
3. direct/FMM pressure and integrated-force agreement at fixed discretization;
4. JVP-versus-directional-difference agreement for tangent-linear components;
5. pressure-force agreement with `ForceMonitor` conventions and independent
   Kutta-Joukowski or surface-vorticity force checks where applicable;
6. wake-length, wake-strength, and shedding-time sensitivity reports;
7. runtime and memory reporting for direct, FMM, and tangent variants;
8. explicit warnings or errors for particle/vector-potential cases where a
   scalar-potential pressure method is incomplete.

## References

- D. Bernoulli, *Hydrodynamica*, 1738.
- H. Lamb, *Hydrodynamics*, 6th ed., Cambridge University Press, 1932.
- G. K. Batchelor, *An Introduction to Fluid Dynamics*, Cambridge University
  Press, 1967.
- J. L. Hess and A. M. O. Smith, "Calculation of Potential Flow About
  Arbitrary Bodies," *Progress in Aeronautical Sciences*, vol. 8, 1967.
- L. Morino and C. C. Kuo, "Subsonic Potential Aerodynamics for Complex
  Configurations: A General Theory," *AIAA Journal*, vol. 12, no. 2, 1974.
- J. Katz and A. Plotkin, *Low-Speed Aerodynamics*, 2nd ed., Cambridge
  University Press, 2001.
- L. Greengard and V. Rokhlin, "A Fast Algorithm for Particle Simulations,"
  *Journal of Computational Physics*, vol. 73, no. 2, 1987.
- S. C. Brenner and L. R. Scott, *The Mathematical Theory of Finite Element
  Methods*, Springer, 2008. Used here for standard graph/least-squares
  elliptic-discretization background, not as panel-method authority.
