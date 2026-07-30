# Prescribed Helical Wake

This note records the wake-shape formulation behind
`examples/rotor_hover_prescribed_helical_wake.jl` and the next solver direction
for item 002. The current implementation uses a helix only as an initial
condition for a finite `PanelWake`; subsequent iterations solve the rotor
against the explicit wake and relax the wake nodes.

## Literature Review: Free-Wake Vortex-Filament Practice

Free-wake vortex-filament rotor methods usually separate the aerodynamic solve
from the wake-geometry solve. A typical iteration is:

1. solve the bound circulation on the lifting surface,
2. shed trailing vorticity from spanwise gradients in bound circulation and shed
   vorticity from time changes in bound circulation, consistent with Kelvin's
   theorem and the Kutta condition,
3. compute induced velocity from the bound and wake vorticity with a
   desingularized Biot-Savart law, and
4. convect, relax, or solve for wake-node positions until the wake kinematics
   are consistent with the velocity field.

The [QBlade LLFVW theory guide](https://docs.qblade.org/src/theory/aerodynamics/lifting_line/lifting_line.html)
is a compact modern description of this workflow. It represents each blade panel
with vortex rings, evaluates induced velocity by Biot-Savart over the bound and
wake vortex elements, desingularizes the core radius, and advances wake nodes
with explicit, trapezoidal, or predictor-corrector velocity integration. The
same free-vortex-wake pattern appears in recent aeroelastic rotor work such as
Rodriguez and Jaworski's strongly coupled floating-wind-turbine framework
([Part 1](https://doi.org/10.1016/j.renene.2019.04.135),
[Part 2](https://doi.org/10.1016/j.renene.2019.04.136)). For rotorcraft,
Bagai and Leishman's "Rotor free-wake modeling using a pseudo-implicit
relaxation algorithm," *Journal of Aircraft*, 1995, is the relevant classic
reference for stabilizing steady rotor wake relaxation without making every
iteration a fully explicit rollup.

VSPAero is a useful production-code reference because it includes both explicit
and implicit steady wake paths. The steady solver loop in
[`VSP_Solver.C`](https://github.com/OpenVSP/OpenVSP/blob/main/src/vsp_aero/Solver/VSP_Solver.C)
iterates `CurrentWakeIteration_` up to `WakeIterations_`, applies
`WakeRelax_`, and can switch to `ImplicitWake_` after
`ImplicitWakeStartIteration_`. Its explicit steady wake update is implemented in
[`Vortex_Trail.C`](https://github.com/OpenVSP/OpenVSP/blob/main/src/vsp_aero/Solver/Vortex_Trail.C):
each wake edge keeps a fixed arclength increment `dS`, forms the local velocity
direction, and relaxes the residual

```math
\mathbf{r}_e =
\mathbf{x}_{i+1} - \mathbf{x}_i - \Delta s_e \hat{\mathbf{u}}_e .
```

For rotor wakes, VSPAero deliberately projects the residual/update onto the
freestream plus thrust-axis downwash subspace. It removes lateral induced
components that would otherwise let a hover rotor wake drift into fountain or
upwash geometries. The same residual is also assembled in
[`VSP_Edge.C`](https://github.com/OpenVSP/OpenVSP/blob/main/src/vsp_aero/Solver/VSP_Edge.C)
for the implicit path. In implicit wake mode, wake-node coordinates are part of
the nonlinear residual equations; the solver computes corrections with its GMRES
infrastructure and preconditioning, then applies relaxed Newton-style coordinate
updates in `UpdateWakeLocations`.

## Mathematical Formulation

Let

```math
X =
\begin{bmatrix}
\mathbf{x}_1 & \mathbf{x}_2 & \cdots & \mathbf{x}_{N_w}
\end{bmatrix}
```

collect the free wake-node coordinates. Let `Γ_w` denote wake-filament
strengths, `q` denote the body unknowns after solving the panel boundary
conditions, and

```math
\mathbf{U}(X, \Gamma_w, q)
```

denote the total velocity operator evaluated at wake nodes or wake edges. This
operator includes freestream, body motion, body-induced velocity, wake-induced
velocity, and any rotor-frame kinematic terms used by the diagnostic.

The FLOWPanel age-marched residual treats neighboring rows as successive wake
ages and asks each node to satisfy a trapezoidal kinematic step:

```math
\mathbf{r}_k(X, \Gamma_w, q)
= \mathbf{X}_k - \mathbf{X}_{k-1}
- \frac{\Delta t}{2}
\left[
\mathbf{U}_{k-1}(X, \Gamma_w, q)
+ \mathbf{U}_{k}(X, \Gamma_w, q)
\right].
```

A VSPAero-style steady arclength residual instead fixes the wake-edge length and
solves only for alignment with the local velocity direction:

```math
\mathbf{r}_e(X, \Gamma_w, q)
= \mathbf{x}_{i+1} - \mathbf{x}_i
- \Delta s_e \hat{\mathbf{U}}_e(X, \Gamma_w, q),
\qquad
\hat{\mathbf{U}}_e =
\frac{\mathbf{U}_e}{\|\mathbf{U}_e\|}.
```

For hover rotors, the operator used in the arclength residual may be projected
before normalization:

```math
\mathbf{U}_e^\parallel =
\mathbf{U}_\infty
+ \min\left[
(\mathbf{U}_e-\mathbf{U}_\infty)\cdot \hat{\mathbf{t}}, 0
\right]\hat{\mathbf{t}},
```

where `\hat{\mathbf{t}}` is the signed rotor thrust or downstream wake axis.
This is a numerical modeling choice, not a general free-wake identity. It is
useful when the immediate problem is to prevent a hover wake from crossing its
anchor plane or growing upstream under local self-induced velocity.

The inner wake-shape solve can be posed as a weighted nonlinear least-squares
problem:

```math
\min_X \frac{1}{2}\|W_X r_X(X, \Gamma_w, q)\|_2^2.
```

Here `r_X` can be the age-marched residual, the arclength residual, or a mixed
residual with separate rows for axial, radial, and azimuthal constraints. If the
diagnostic also needs to measure wake-strength consistency, add

```math
\frac{1}{2}\|W_\Gamma(\Gamma_w - S(q))\|_2^2,
```

where `S(q)` is the Kelvin/Kutta shedding operator implied by the current body
solution. For item 002, geometry residual contraction should be established
before coupling this strength term into the same inner solve.

The current diagnostics map cleanly onto this notation:

- `max_target_residual_R` is an infinity norm of the requested residual,
  normalized by rotor radius.
- `capped_node_fraction` reports the fraction of coordinates or nodes where the
  proposed correction exceeded the step cap, so the actual map no longer
  follows the residual descent direction.
- `max_applied_step_R` is the largest accepted coordinate update, normalized by
  radius.
- `delta_gamma_rel` tracks outer wake-strength history changes and should be
  interpreted separately from wake-geometry residual contraction.
- Fountain and anchor-crossing checks are geometry validity guards on `X`, not
  residual norms. They should remain hard globalization constraints even when the
  inner residual formulation changes.

## Matrix-Free Newton-GMRES Wake Solve

A Newton-GMRES wake solve treats the wake coordinates `X` as the nonlinear
unknowns for the inner wake-shape problem. GMRES does not solve the nonlinear
problem directly. At Newton iteration `n`, it solves only the linearized
correction equation

```math
J(X^n)\delta X = -r(X^n),
\qquad
X^{n+1} = X^n + \alpha \delta X .
```

Here `\alpha` is supplied by a line search, trust-region step cap, or both. The
matrix-free operator needed by GMRES is the residual directional derivative

```math
J(X)v =
\left.\frac{d}{d\epsilon}r(X+\epsilon v)\right|_{\epsilon=0}.
```

For the VSPAero-style arclength residual,

```math
r_e(X) =
x_{i+1} - x_i - \Delta s_e\hat U_e(X),
\qquad
\hat U_e = \frac{U_e}{\|U_e\|},
```

the Jacobian-vector product for a coordinate perturbation `v = \delta X` is

```math
(Jv)_e =
\delta x_{i+1} - \delta x_i
- \delta(\Delta s_e)\hat U_e
- \Delta s_e\delta\hat U_e .
```

The first prototype should hold `\Delta s_e` fixed during the inner wake-shape
solve, which drops the `\delta(\Delta s_e)` term. If the edge velocity norm is
nonzero, the normalized-velocity perturbation is

```math
\delta\hat U_e =
\frac{1}{\|U_e\|}
\left(I-\hat U_e\hat U_e^T\right)\delta U_e .
```

This follows by differentiating the normalized vector

```math
\hat U_e = U_e s^{-1},
\qquad
s = \|U_e\| = (U_e^T U_e)^{1/2}.
```

The norm perturbation is

```math
\delta s
= \frac{U_e^T\delta U_e}{\|U_e\|},
```

so

```math
\delta(s^{-1})
= -\frac{U_e^T\delta U_e}{\|U_e\|^3}.
```

Substituting into `\delta\hat U_e = \delta U_e\,s^{-1}
+ U_e\,\delta(s^{-1})` gives

```math
\delta\hat U_e
= \frac{\delta U_e}{\|U_e\|}
- U_e\frac{U_e^T\delta U_e}{\|U_e\|^3}
= \frac{1}{\|U_e\|}
\left(I-\hat U_e\hat U_e^T\right)\delta U_e .
```

The projection `I-\hat U_e\hat U_e^T` removes the perturbation component
parallel to `U_e`; to first order, only the perpendicular component changes the
velocity direction.

With the body solution and wake strengths frozen inside this inner solve,

```math
\delta U_e =
\nabla_y U_\mathrm{body}(y_e)\delta y_e
+ \sum_m \Gamma_m \delta K(y_e; a_m,b_m).
```

If the edge velocity is evaluated at the midpoint,

```math
y_e = \frac{x_i+x_{i+1}}{2},
\qquad
\delta y_e = \frac{\delta x_i+\delta x_{i+1}}{2}.
```

For a wake segment from endpoint `a` to endpoint `b`, define

```math
r_1 = y-a,\quad r_2 = y-b,\quad r_0=b-a,
\quad c=r_1\times r_2,\quad D=c\cdot c+\epsilon^2,
```

```math
h =
r_0\cdot
\left(
\frac{r_1}{\|r_1\|}
-
\frac{r_2}{\|r_2\|}
\right),
\qquad
K = \frac{1}{4\pi}\frac{hc}{D}.
```

The segment perturbations are

```math
\delta r_1=\delta y-\delta a,\quad
\delta r_2=\delta y-\delta b,\quad
\delta r_0=\delta b-\delta a,
```

```math
\delta c = \delta r_1\times r_2 + r_1\times\delta r_2,
\qquad
\delta D = 2c\cdot\delta c,
```

and

```math
\delta h =
\delta r_0\cdot
\left(
\frac{r_1}{\|r_1\|}
-
\frac{r_2}{\|r_2\|}
\right)
+
r_0\cdot
\left[
\frac{(I-\hat r_1\hat r_1^T)\delta r_1}{\|r_1\|}
-
\frac{(I-\hat r_2\hat r_2^T)\delta r_2}{\|r_2\|}
\right],
```

where `\hat r_i = r_i/\|r_i\|`. The resulting matrix-free segment-kernel
derivative is

```math
\delta K =
\frac{1}{4\pi}
\left[
\frac{\delta h\,c+h\,\delta c}{D}
-
\frac{hc\,\delta D}{D^2}
\right].
```

If wake strengths are later included in the nonlinear unknown vector, each
segment contribution also gains the strength perturbation term

```math
\delta U_{e,m} = \Gamma_m\delta K_m + \delta\Gamma_m K_m .
```

A practical implementation should use the following outer-to-inner structure:

1. Freeze or update `q` and `Γ_w` at the outer wake-strength/body-solve level.
2. Evaluate the current wake residual `r(X)`.
3. Define a matrix-free `Jv(v)` callback from the directional derivative above.
4. Use GMRES to solve the linear correction approximately.
5. Apply step caps and line-search or trust-region damping.
6. Reject updates that worsen the residual norm or violate fountain and anchor
   checks.
7. Repeat until both the residual and the applied-step diagnostics contract.

## Solver Formulations and Recommendation

The current Picard/under-relaxed explicit update is cheap and easy to diagnose:
compute a target wake position from the velocity field, cap the step, relax it,
and repeat. It becomes fragile when every useful node is capped, because the
iteration is then governed by the cap geometry rather than by the residual
model. In that regime, changing `WakeRelax` alone is unlikely to fix item 002's
observed capped residual.

The recommended next formulation is a damped nonlinear least-squares inner
wake-shape solve. Start with one of the residuals above, retain the existing
step caps as trust-region limits, and use a line search that rejects updates
which increase residual norm or violate the fountain/anchor checks. This keeps
the implementation local to the diagnostic harness while making the residual the
object being minimized, rather than an indirect consequence of repeated capped
node targets.

A VSPAero-style implicit wake Newton/GMRES formulation is the higher-effort
production direction. In that approach, wake-node coordinates are added to the
coupled nonlinear residual beside the panel equations, with matrix-free products
and preconditioning for the enlarged system. This is the right long-term path if
FLOWPanel needs robust steady wake relaxation across general lifting geometries,
but it is more invasive than item 002 requires.

Anderson or Aitken acceleration is a lower-effort improvement only after the
explicit map is no longer fully capped. Acceleration helps a contractive fixed
point map. It does not repair a map whose accepted update is dominated by
globalization limits at nearly every node. Keep the outer wake-strength history
updates separate until the geometry residuals visibly contract.

## Frozen-Field Streamline Update

The `nested_pitch` diagnostic can also rebuild the inner wake geometry by
streamline integration through a velocity field frozen at the start of each
inner geometry iteration. The outer iteration first solves the body against the
current wake and pins the latest shed strengths. During the inner loop, body
strengths, wake strengths, the row-1 `TE + Das` anchors, and the wake geometry
used as the source field are held fixed. Candidate downstream rows are then
generated from the anchors with the physical wake-row time step.

For a global query point `x_g`, the frozen fluid velocity is

```math
\mathbf{U}_g(\mathbf{x}_g)
= \mathbf{U}_\infty
+ \mathbf{U}_\mathrm{body}(\mathbf{x}_g)
+ \mathbf{U}_\mathrm{wake}(\mathbf{x}_g).
```

Wake nodes are material points relative to the moving blade frame, so the
streamline update subtracts the rigid material velocity of the active frame
chain:

```math
\mathbf{V}_{\mathrm{frame},g}(\mathbf{x}_g)
= \mathbf{v}_g
+ \boldsymbol{\omega}_g \times (\mathbf{x}_g-\mathbf{o}_g),
```

accumulated recursively over the same `ReferenceFrame` dependencies used by
`kinematic_velocity!`. The global-frame effective velocity is therefore

```math
\dot{\mathbf{x}}_g
= \mathbf{U}_{\mathrm{eff},g}
= \mathbf{U}_g - \mathbf{V}_{\mathrm{frame},g}.
```

Equivalently, in body-frame coordinates `xi`,

```math
\dot{\boldsymbol{\xi}}
= R_{g2f}
\left[
\mathbf{U}_g
- \mathbf{v}_g
- \boldsymbol{\omega}_g \times (\mathbf{x}_g-\mathbf{o}_g)
\right],
```

with the global point reconstructed as

```math
\mathbf{x}_g = \mathbf{o}_g + R_{f2g}\boldsymbol{\xi}.
```

The current diagnostic integrates each row with RK2 midpoint steps:

```math
\mathbf{x}_{k+1}
= \mathbf{x}_k
+ \Delta t_\mathrm{row}\,
\mathbf{U}_{\mathrm{eff},g}
\left(
\mathbf{x}_k
+ \frac{\Delta t_\mathrm{row}}{2}
\mathbf{U}_{\mathrm{eff},g}(\mathbf{x}_k)
\right),
```

using the frozen source geometry for both velocity evaluations. Convergence is
reported as the maximum candidate-node displacement normalized by rotor radius;
fountain and anchor-crossing checks remain hard rejection diagnostics.

## Axial Helix Pitch

Let the rotor angular speed magnitude be

```math
\Omega = |\omega|,
```

with tip speed

```math
V_\mathrm{tip} = \Omega R .
```

For an axisymmetric prescribed helix, a wake node of age `theta` radians is
rotated by `theta` about the rotor axis and convected axially by

```math
\Delta x = s_x U_p \frac{\theta}{\Omega},
```

where `s_x` is the downstream axial sign and `U_p` is the axial pitch speed used
to initialize the helix. The axial advance per revolution is therefore

```math
\Delta x_{2\pi}
= s_x U_p \frac{2\pi}{\Omega}.
```

Normalizing by radius gives the pitch per revolution

```math
\frac{|\Delta x_{2\pi}|}{R}
= 2\pi \frac{U_p}{\Omega R}
= 2\pi \lambda_p,
```

where

```math
\lambda_p = \frac{U_p}{V_\mathrm{tip}} .
```

The code parameterizes the pitch by this nondimensional speed ratio rather than
by dimensional pitch distance. For row spacing `N` rows per revolution, the
initial axial spacing between adjacent wake rows is

```math
\frac{\Delta x_\mathrm{row}}{R}
= \frac{2\pi}{N}\lambda_p .
```

## Hover Seed

In hover, there is no imposed axial freestream in the baseline case, so the
initial helix pitch comes only from the induced inflow seed:

```math
U_p = \lambda_i V_\mathrm{tip},
\qquad
\lambda_p = \lambda_i .
```

In the example this seed is `INITIAL_INFLOW`, defaulting to `0.08`.

## Axial Advance-Ratio Seed

For the axial-flow diagnostic, the freestream is aligned with the same signed
downstream rotor axis used by the wake:

```math
\mathbf{V}_\infty
= s_x \mu_x V_\mathrm{tip}\,\hat{\mathbf{x}},
```

where `mu_x` is the environment-controlled axial advance ratio
`AXIAL_ADVANCE_RATIO`. The initial wake should convect with the imposed axial
freestream plus the induced inflow seed, so

```math
U_p
= (\mu_x + \lambda_i)V_\mathrm{tip}.
```

Thus the total initial helix pitch ratio is

```math
\boxed{\lambda_p = \mu_x + \lambda_i.}
```

With the first diagnostic value `AXIAL_ADVANCE_RATIO=0.10` and the default
`INITIAL_INFLOW=0.08`, the initialized helix uses

```math
\lambda_p = 0.10 + 0.08 = 0.18.
```

This keeps the wake axisymmetric for the axial test. It does not model edgewise
forward flight, where the wake geometry is no longer a simple axisymmetric
helix.

## References

- QBlade documentation, [Lifting Line Free Vortex Wake](https://docs.qblade.org/src/theory/aerodynamics/lifting_line/lifting_line.html).
- S. N. Rodriguez and J. W. Jaworski, "Strongly-coupled aeroelastic free-vortex
  wake framework for floating offshore wind turbine rotors. Part 1: Numerical
  framework and validation," *Renewable Energy*, 2019.
- S. N. Rodriguez and J. W. Jaworski, "Strongly-coupled aeroelastic free-vortex
  wake framework for floating offshore wind turbine rotors. Part 2:
  Application," *Renewable Energy*, 2019.
- A. Bagai and J. G. Leishman, "Rotor free-wake modeling using a
  pseudo-implicit relaxation algorithm," *Journal of Aircraft*, 1995.
- OpenVSP source, VSPAero steady and implicit wake implementation:
  [`VSP_Solver.C`](https://github.com/OpenVSP/OpenVSP/blob/main/src/vsp_aero/Solver/VSP_Solver.C),
  [`Vortex_Trail.C`](https://github.com/OpenVSP/OpenVSP/blob/main/src/vsp_aero/Solver/Vortex_Trail.C),
  and [`VSP_Edge.C`](https://github.com/OpenVSP/OpenVSP/blob/main/src/vsp_aero/Solver/VSP_Edge.C).
