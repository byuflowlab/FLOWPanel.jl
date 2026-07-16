# Unsteady Bernoulli Pressure

> **2026-07-11 revision.** The reviewed implementation contract is
> [`theory/unsteady_pressure_monitors_2026-07-11.md`](../theory/unsteady_pressure_monitors_2026-07-11.md).
> The kinetic term now uses reconstructed inertial surface velocity; the first
> sample seeds a zero derivative; later samples use backward Euler and then
> variable-step BDF2; history stores the exterior potential; and
> particle/vector-potential-only wake sources require explicit opt-in to a
> warned mixed partial diagnostic. First-order formulas below remain as
> derivational background, not the current temporal scheme.

## Steady Mode

With `unsteady=false`, `PressureBernoulli` evaluates the **body-relative steady
loading formulation**

```math
P = \frac{1}{2}\rho\left(U_\infty^2 - \left|\mathbf{u}_{\mathrm{rel},t}\right|^2\right),
```

where `u_rel,t` is the tangential projection of `body.velocity` (which is
body-relative during monitor execution; see "FLOWPanel Variables" below). This
is valid for flows steady in the body frame, such as a rotor at constant
rotation rate. The complete steady relation in a rotating frame also carries a
centrifugal/reference-potential term `½ρ|w|²`, omitted here: the steady
pressure is defined only up to that rotating-frame reference contribution. The
omitted term is symmetric across a blade section and loading-neutral, so
integrated loads match the historically validated rotor results, but the steady
field is not the complete absolute pressure for every rotating body.

The reconstructed **inertial** surface velocity enters the kinetic term only in
`unsteady=true` mode, where the ALE `∂φ/∂t` term compensates for the frame
motion. Using the inertial kinetic energy *without* that term cancels the
first-order blade loading on a rotating body.

> **Regression note.** Between commit `ef1fe1e` (2026-07-14) and this fix, the
> steady mode incorrectly used the inertial surface velocity. Any
> steady-Bernoulli pressure or force computed on a moving body with code from
> that window is invalid (observed: rotor axial-flow CT collapse).

This note records the discrete unsteady term used by `PressureBernoulli`.
The monitor owns this calculation because pressure is a post-processing
quantity, while `simulate!` should only advance the aerodynamic state needed by
the flow solve.

## Arbitrary Lagrangian--Eulerian Moving-Point Derivative

For an inviscid, incompressible potential flow with scalar potential `phi`,

```math
\mathbf{u} = \nabla \phi .
```

The unsteady Bernoulli pressure contribution is

```math
p =
\frac{1}{2}\rho\left(U_\infty^2 - \left|\mathbf{u}_\mathrm{rel}\right|^2\right)
- \rho \frac{\partial \phi}{\partial t}.
```

In the Arbitrary Lagrangian--Eulerian (ALE) description, panel potentials are
sampled at control points that move with the body. If a
control point has inertial velocity `w`, the time derivative seen by the moving
point is

```math
D_g\phi := \frac{d}{dt}\phi(\mathbf{x}_p(t),t)
=
\frac{\partial \phi}{\partial t}
+ \mathbf{w}\cdot\nabla\phi .
```

Therefore the Eulerian partial-time derivative needed by Bernoulli is

```math
\boxed{
\frac{\partial \phi}{\partial t}
=
\frac{d}{dt}\phi(\mathbf{x}_p(t),t)
- \mathbf{w}\cdot\nabla\phi
}
```

and the first-order discrete monitor formula is

```math
\boxed{
\left.\frac{\partial \phi}{\partial t}\right|_i^n
\approx
\frac{\phi_i^n-\phi_i^{n-1}}{\Delta t}
- \mathbf{w}_i^n\cdot\nabla\phi_i^n .
}
```

The finite difference is panel-following: both `phi_i^n` and
`phi_i^{n-1}` are values stored at the same body panel's control point at two
successive body positions.

## FLOWPanel Variables

`kinematic_velocity!` stores the body control-point velocity in
`body.velocity_kinematic` and subtracts it from `body.velocity`. Thus, during
monitor execution,

```math
\mathbf{w}_i = \texttt{body.velocity\_kinematic[:, i]},
```

and the inertial scalar-potential gradient is

```math
\nabla\phi_i
=
\texttt{body.velocity[:, i]}
+
\texttt{body.velocity\_kinematic[:, i]} .
```

The implemented expression is therefore

```math
\boxed{
\phi_{\mathrm{dot},i}
=
\frac{\phi_i^n-\phi_i^{n-1}}{\Delta t}
-
\texttt{velocity\_kinematic}_i \cdot
\left(
    \texttt{velocity}_i + \texttt{velocity\_kinematic}_i
\right)
}
```

This is not equivalent to subtracting
`velocity_kinematic dot velocity` in the current codebase, because
`body.velocity` is body-relative rather than inertial.

## Scalar-Potential Sources

The scalar potential used in this finite difference is evaluated inside
`PressureBernoulli` only when `unsteady=true`. The monitor computes potential at
body control points using `FastMultipole.ProbeSystem`, then adds the freestream
potential `uinf dot x`.

The stored history is the exterior surface potential. Exact-control-point body
kernels return FLOWPanel's canonical interior limit, so bodies with a doublet
or ring strength use ``\phi^+=\phi^- - \mu_\mathrm{code}`` at the target panel.
Source-only bodies, other bodies' off-surface contributions, and wake
potentials are unchanged.

Only sources with a well-defined scalar potential are included. The complete
body and wake source lists are checked with
`FastMultipole.has_vector_potential(source)`;
sources that induce only a vector potential, such as vortex filaments or vortex
particles, are excluded. This avoids requesting scalar potential from a vector
potential source, which can corrupt the FMM scalar-potential result.
Vector-capable bodies are unsupported. Vector-only wake sources throw by
default; `allow_partial=true` opts into a once-warned partial diagnostic. Their
velocity remains in the total inertial kinetic-energy term but is subtracted
from the retained-scalar ``\mathbf{w}\cdot\nabla\phi`` contraction.

Trailing-edge pressure-pair averaging is an optional heuristic, disabled by
default. Set `correct_kuttacondition=true` to enable it.

The normal velocity/gradient influence calls in `simulate!` should not be used
for this bookkeeping. Keeping scalar-potential evaluation monitor-local avoids
extra cost for steady-pressure runs and prevents invalid scalar-potential
requests from being mixed into the main wake-on-all or system-on-all influence
calls.
