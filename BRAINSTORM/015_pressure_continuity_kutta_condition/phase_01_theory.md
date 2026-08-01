# Phase 1 — Theory Derivation

**Status:** APPROVED — 2026-07-29

**Prerequisite:** none

**Phase approval:** [x]

## Scope and prior-phase handoff

This phase establishes the theory for independently configurable near-wake
attachment and Kutta closure on one finite-wake, paired-edge, Dirichlet
`RigidWakeBody`. It covers an opt-in nonlinear pressure-continuity closure and
a `Das`-free TE-anchored wake, with eventual body-solver support for
`Backslash`, `KrylovSolver`, and a fully constructed `FGSSolver`. Existing
formulations and defaults must remain unchanged.

There is no prior phase. The motivation comes from Item 014: the first
attached-wake length `Das` is the largest discretization sensitivity yet seen
in the rotor-hover study. Floor-free `Das` changes moved settled \(C_T\) by
about 37%, primarily through bound circulation. The existing linear closure
ties each wake strength to the upper-minus-lower body doublet jump. A
pressure-continuity Kutta condition would instead select the circulation for
which the two sides of each paired trailing edge have equal pressure.

This is not the postprocessing option `correct_kuttacondition=true`. That
option merely replaces two already computed centroid pressures by their
average; it does not change the body solution, wake strength, surface
velocity, or loads consistently.

## Deliverables

- Complete discrete theory, sign conventions, solvability, residual scaling,
  state invariants, and pressure-provider contract.
- Independent wake-attachment and Kutta-closure axes with a four-combination
  attribution matrix.
- Route B newborn-wake deposition alternatives B1–B4 and a selected
  experimental order.
- Recommended primary and fallback nonlinear solvers.
- An explicit decision record suitable as the Phase 2 architecture handoff.

## 1. Existing discrete operators and signs

The following definitions match `src/FLOWPanel_formulation.jl`,
`src/FLOWPanel_solver.jl`, `src/FLOWPanel_liftingbody.jl`, and
`docs/wake_solve_schemes.md`.

- \(N\): number of body panels.
- \(M\): number of **paired** shedding edges.
- \(\boldsymbol\mu\in\mathbb R^N\): body doublet strengths. FLOWPanel uses
  \(\mu=\phi_{\rm inside}-\phi_{\rm outside}\), with the interior
  \(+\mu_i/2\) self limit already included in the matrix.
- \(\boldsymbol\sigma=-\boldsymbol u_{\rm known}\cdot\boldsymbol n\):
  prescribed source strengths from freestream, kinematics, and the frozen
  free-wake velocity.
- \(B\in\mathbb R^{N\times N}\): body-only interior doublet-potential
  operator.
- \(C\in\mathbb R^{M\times N}\): shedding-edge map,
  \[
  (C\boldsymbol\mu)_e=\mu_{\mathrm{upper},e}
  -\mu_{\mathrm{lower},e}.
  \]
- \(W(D_a)\in\mathbb R^{N\times M}\): potential at body centroids from
  unit-strength attached transition panels spanning
  \(\mathrm{TE}\rightarrow\mathrm{TE}+D_a\).
- \(G=B+WC\): the currently assembled finite-wake Dirichlet operator.
- \(\boldsymbol b=-S\boldsymbol\sigma\): the prescribed Dirichlet
  right-hand side for the selected wake-to-body formulation.

The exact definition of \(\boldsymbol b\) must remain formulation-specific.
For the current production `VelocityThroughSources` route it contains the
free-wake velocity through \(\boldsymbol\sigma\), not its scalar-potential
trace. Pressure continuity must not silently change that choice.

For the retained-transition/pressure-continuity combination developed first
below, the existing affine correction channel uses one correction \(c_e\) per
shedding edge:

\[
G\boldsymbol\mu(\boldsymbol c)
  =\boldsymbol b+W\boldsymbol c,
\qquad
\boldsymbol\gamma(\boldsymbol c)
  =C\boldsymbol\mu(\boldsymbol c)-\boldsymbol c.
\tag{1}
\]

A positive \(c_e\) subtracts circulation: the upper transition strength is
shifted by \(-c_e/2\), the lower by \(+c_e/2\), and their difference is
\(\gamma_e\). Substituting \(G=B+WC\) into (1) gives

\[
B\boldsymbol\mu(\boldsymbol c)
  +W\boldsymbol\gamma(\boldsymbol c)=\boldsymbol b.
\tag{2}
\]

Equation (2) is the key consistency identity for that combination: every
nonlinear trial must change the body solution and attached-wake strength
together. Section 5 generalizes the same requirement to both attachment
routes and both Kutta closures. Changing \(\gamma\) after a body solve would
violate the same boundary equation.

The downstream marched wake is frozen during the nonlinear trials. Once a
trial is accepted, the usual physical-time shedding operation copies
\(\gamma_e\) into the new first free-wake row. Older wake rows are history and
must not be modified by the Kutta iteration.

## 2. Trial velocity and pressure residual

Let \(D_{\rm e}\) evaluate the exterior body-doublet velocity, including the
same surface-gradient/half-jump reconstruction used by the accepted
simulation, and let \(W_{\rm e}\) evaluate the attached transition-panel
velocity. All prescribed contributions—freestream, kinematics, source panels,
and the old free wake—are collected in \(\boldsymbol u_{\rm known}\). The
complete trial surface velocity is

\[
\boldsymbol u(\boldsymbol c)
=\boldsymbol u_{\rm known}
 D_{\rm e}\boldsymbol\mu(\boldsymbol c)
 W_{\rm e}\boldsymbol\gamma(\boldsymbol c).
\tag{3}
\]

This expression is conceptual, not permission to assemble dense velocity
operators. The implementation may use the existing influence and
`compute_mu_gradient!` paths, but it must reset to the same snapshot before
each trial so velocity does not accumulate.

A trial-safe pressure provider maps the full trial state to one pressure per
body panel:

\[
\boldsymbol p(\boldsymbol c)
=\mathcal P\!\left(
\boldsymbol u(\boldsymbol c),\boldsymbol\mu(\boldsymbol c),
\boldsymbol\gamma(\boldsymbol c);
\mathcal H^n\right),
\tag{4}
\]

where \(\mathcal H^n\) is frozen physical-time history. The paired-edge
residual is

\[
\boxed{
r_e(\boldsymbol c)
=p_{\mathrm{upper},e}(\boldsymbol c)
-p_{\mathrm{lower},e}(\boldsymbol c),\qquad e=1,\ldots,M.
}
\tag{5}
\]

Thus there are \(M\) nonlinear unknowns and \(M\) equations. The pressure
provider returns panel pressures; pairing and subtraction belong to the
formulation, not to the provider.

### Steady Bernoulli specialization

For the first built-in pressure law,

\[
p_i=p_{\rm ref}
+\frac{\rho}{2}\left(U_{\rm ref}^2-\lVert\boldsymbol u_i\rVert^2\right),
\]

so the common reference cancels:

\[
r_e(\boldsymbol c)
=\frac{\rho}{2}
\left(
\lVert\boldsymbol u_{\mathrm{lower},e}(\boldsymbol c)\rVert^2
-\lVert\boldsymbol u_{\mathrm{upper},e}(\boldsymbol c)\rVert^2
\right).
\tag{6}
\]

For frozen geometry and prescribed wake history, (1) and (3) are affine in
\(\boldsymbol c\). Write

\[
\frac{\partial\boldsymbol\mu}{\partial c_j}=G^{-1}W_{:j},
\qquad
\frac{\partial\boldsymbol\gamma}{\partial c_j}
=C G^{-1}W_{:j}-\boldsymbol e_j,
\]

and

\[
\boldsymbol h_{ij}
=\frac{\partial\boldsymbol u_i}{\partial c_j}
=D_{{\rm e},i}G^{-1}W_{:j}
+W_{{\rm e},i}\left(CG^{-1}W_{:j}-\boldsymbol e_j\right).
\tag{7}
\]

The steady residual is therefore quadratic, with exact Jacobian

\[
J_{ej}
=\rho\left(
\boldsymbol u_{\mathrm{lower},e}\cdot\boldsymbol h_{\mathrm{lower},e,j}
-\boldsymbol u_{\mathrm{upper},e}\cdot\boldsymbol h_{\mathrm{upper},e,j}
\right).
\tag{8}
\]

Equation (8) is a verification oracle and a possible later optimization.
Building all response columns costs at least \(M\) body solves, so it is not
the recommended first implementation.

## 3. Gauges, normalization, and solvability

### Potential and pressure gauges

The retained-transition formulation solves the same \(G\) systems as the
existing finite-wake formulation; it introduces no new body-potential gauge.
The correction vector is physical—it changes circulation—and has no arbitrary
additive constant.

A common additive pressure reference cancels exactly in (5). This is distinct
from numerical normalization: a nonlinear solver still needs pressure and
correction scales to compare residuals and step sizes.

If the body-only \(B\) operator is used by a future no-transition
representation, its constant-potential mode must be tested explicitly. An
area-weighted mean-potential constraint is acceptable because
\(C\boldsymbol 1=0\) for fully paired edges, so it fixes representation
without changing circulation.

### Residual scales and convergence norms

The provider contract should return dimensional pressure. Nonlinear
convergence should use a strictly positive pressure scale \(p_s\), supplied
explicitly or derived from a documented characteristic velocity, and a
strictly positive correction scale \(c_s\). In hover, freestream dynamic
pressure is zero and is therefore not a valid universal default.

The primary pressure gate should be

\[
\lVert \boldsymbol r/p_s\rVert_\infty\le\varepsilon_p,
\tag{9}
\]

because every shedding edge must pass. Also report a span/edge-weighted
\(L_2\) norm so mesh-local noise can be distinguished from a global failure.
Require a step gate as well,

\[
\lVert\Delta\boldsymbol c/c_s\rVert_\infty\le\varepsilon_c,
\tag{10}
\]

and require the inner body solve to satisfy its own convergence test. Scaling
the pressure residual makes the correction invariant to a uniform density or
pressure-unit rescaling, subject to nonlinear tolerances.

### Local solvability

After eliminating \(\boldsymbol\mu\), local solvability requires the
\(M\times M\) Jacobian \(J=\partial\boldsymbol r/\partial\boldsymbol c\) to
have full rank near the root. A singular or poorly conditioned \(J\) can arise
from a degenerate edge pair, a vanishing local velocity scale, duplicate
constraints, insufficient wake leverage, or nearly indistinguishable
spanwise modes. These are formulation failures, not reasons to silently
accept a least-squares pressure mismatch.

All first-version shedding edges must be paired. An unpaired edge has no
defined lower centroid pressure and therefore no residual under (5).

### Spanwise coupling

Although residual \(r_e\) samples one upper/lower pair, changing \(c_j\)
changes the global body solution through \(G^{-1}W_{:j}\) and changes the
wake-induced field. The Jacobian is generally dense. Independent per-edge
secant solves are not a consistent algorithm except as a diagonal
preconditioner or diagnostic.

## 4. What condition is actually being enforced?

Equation (5) matches pressures at the centroids of the two panels adjacent to
each stored shedding edge. It is a discrete surrogate for the sharp-edge
continuum condition

\[
\lim_{\boldsymbol x\to\mathrm{TE}^+}p_{\rm upper}(\boldsymbol x)
=\lim_{\boldsymbol x\to\mathrm{TE}^-}p_{\rm lower}(\boldsymbol x).
\]

The two are not identical at finite resolution:

- the centroids are displaced from the geometric trailing edge;
- upper and lower centroid locations, normals, areas, and velocity
  reconstruction stencils differ;
- kernel offsets regularize the sampled velocity;
- constant-strength panels cannot reproduce the exact edge asymptotics; and
- matching one centroid pair does not prove bounded velocity at every point
  approaching the edge.

Consequently, a small nonlinear residual is an algebraic convergence result,
not physical validation. Chordwise/TE mesh refinement and kernel-offset
studies are mandatory in Phase 4.

## 5. Separate wake attachment from Kutta enforcement

The near-wake geometry and the circulation closure are two different
discretization choices. They must not be represented by one inseparable
“formulation” switch.

Let \(W_R\) denote the wake operator selected by a **wake-attachment route**
\(R\), and let \(K\) denote the selected **Kutta closure**. The common body
equation is

\[
B\boldsymbol\mu+W_R\boldsymbol\gamma=\boldsymbol b.
\tag{11}
\]

Two attachment routes and two closures give four meaningful combinations:

| Wake attachment | Kutta closure | Discrete system | Purpose |
| --- | --- | --- | --- |
| A: rigid `Das` transition | Jump, \(\gamma=C\mu\) | \((B+W_A(D_a)C)\mu=b\) | Current baseline |
| A: rigid `Das` transition | Pressure continuity | \(B\mu+W_A(D_a)\gamma=b,\ r(\mu,\gamma)=0\) | Isolate closure sensitivity at fixed geometry |
| B: TE-anchored wake | Jump, \(\gamma=C\mu\) | \((B+W_BC)\mu=b\) | Isolate and potentially eliminate `Das` |
| B: TE-anchored wake | Pressure continuity | \(B\mu+W_B\gamma=b,\ r(\mu,\gamma)=0\) | Fully composed `Das`-free nonlinear option |

This factorial separation is the central Phase 1 revision. Comparing A/jump
against A/pressure changes only Kutta enforcement. Comparing A/jump against
B/jump changes only wake attachment. The fourth combination tests whether the
two changes compose without introducing an interaction.

### Kutta closure J — deterministic doublet jump

The classical discrete rule is

\[
\boxed{\boldsymbol\gamma=C\boldsymbol\mu.}
\tag{12}
\]

The wake strength is not an additional unknown: it is solved simultaneously
and implicitly with the body strengths through \(W_RC\). Solving the body
first and assigning \(\gamma=C\mu\) afterward would omit the new wake panel's
same-step influence and create a one-step-lagged splitting.

This rule is already deterministic in the current code. It follows that item
014's `Das` sensitivity does not demonstrate ambiguity in how circulation is
computed; it demonstrates sensitivity to the geometry through which that
circulation acts.

### Kutta closure P — pressure-selected circulation

Pressure continuity treats \(\boldsymbol\gamma\) as independent while solving
(11) and (5). For Route A, the existing correction coordinates are

\[
\boldsymbol c=C\boldsymbol\mu-\boldsymbol\gamma,
\]

which recover (1). For Route B it may be simpler to iterate directly on
\(\boldsymbol\gamma\), although Phase 2 should choose coordinates that let
both attachment routes use the same nonlinear driver.

At finite resolution, jump closure and centroid pressure closure are
alternative square closures. Imposing both as independent exact equations
would generally overconstrain the discrete problem. Their agreement under TE
mesh refinement is a verification target, not an assumption.

### Route A — retain the rigid `Das` transition

Here \(W_A(D_a)\) is the existing attached panel from
`TE` to `TE + Das`. Advantages:

- the current A/jump combination is the regression baseline;
- the existing affine correction gives a mechanically direct A/pressure
  implementation;
- the shedding map, circulation sign, and downstream output paths already
  exist; and
- A/jump versus A/pressure cleanly measures the Kutta-closure effect.

Its decisive limitation is that the potential and exterior-velocity operators
remain functions of `Das`. Pressure continuity therefore **cannot eliminate
`Das` dependence by construction**. It may reduce the sensitivity, leave it
unchanged, or relocate it.

### Route B — TE-anchored wake with no rigid `Das`

Route B removes the body-owned `TE` to `TE + Das` panel. The wake owns a node
row at the geometric trailing edge, and its first live panel connects that row
to the first downstream wake row. There is no separate offset magnitude.

For a prescribed steady wake, every wake panel that carries the current
steady circulation belongs to \(W_B\). With jump closure,

\[
\left(B+W_{B,\mathrm{steady}}C\right)\boldsymbol\mu
=\boldsymbol b.
\tag{13}
\]

If wake geometry is relaxed, (13) is the inner strength solve and geometry is
updated in an outer wake-shape iteration.

For an unsteady step \(n+1\), previously shed panels retain their accepted
strengths and contribute known influence to \(\boldsymbol b^{n+1}\). The new
panel connects the current TE to the first previously shed row after that row
has been convected to the new time. Only this live panel depends on the current
strength:

\[
\left(B^{n+1}+W_{B,\mathrm{new}}^{n+1}C\right)
\boldsymbol\mu^{n+1}
=\boldsymbol b^{n+1}(\text{old wake}),
\qquad
\boldsymbol\gamma_{\mathrm{new}}^{n+1}
=C\boldsymbol\mu^{n+1}.
\tag{14}
\]

After acceptance, the new row is committed once. Spanwise differences in
\(\gamma_{\mathrm{new}}\) generate trailing vorticity; differences between
successive time rows generate the unsteady shed/start-vortex system. Old wake
strengths are never re-solved or rewritten.

Route B eliminates the *independent rigid* `Das` parameter. It does not
eliminate numerical resolution: the first live panel length is set by actual
convection and timestep, and its geometry must converge under temporal and
wake refinement.

The current `PanelWake` starts at `TE + Das`, so Route B needs:

1. a TE node row owned by the wake;
2. a live first panel included in the same body solve;
3. an explicit unsteady ordering for convection, body motion, solve, and
   commit;
4. a deterministic startup rule before a previously shed downstream row
   exists;
5. no double counting between body-attached and wake-owned panels;
6. a gauge-fixed body-only solve if \(B\) exposes its constant mode;
7. Kelvin/history semantics that commit only the accepted circulation; and
8. matrix/factorization rebuilding whenever the live panel geometry changes.

The startup choices to assess in Phase 2 are: initialize a prior TE row from
known motion, perform a no-live-panel initial solve and then seed the first
row, or start from a supplied steady wake. A zero-length singular panel is not
an acceptable implicit startup convention.

### Open question — deposition of vorticity created during the newest timestep

Katz and Plotkin's low-order unsteady examples place the newest lumped wake
vortex about \(0.2\)–\(0.3\) of one timestep's travel behind the trailing edge;
their archived Program 15 uses
`DXW = 0.3*UT*DT`. The rationale is numerical: vorticity is produced
continuously over the timestep, and a discrete point placed at the apparently
natural midpoint can underrepresent its strongly weighted near-TE influence.
See [Katz–Plotkin Program
15](https://raw.githubusercontent.com/cibinjoseph/KatzPlotkin/main/p15/p15.f)
and the official [unsteady-potential-flow chapter
listing](https://www.cambridge.org/core/books/lowspeed-aerodynamics/077FAF851C4582F1B7593809752C44AE/listing).

This advice must not be silently discarded, but it also does not establish
that FLOWPanel needs its present rigid `Das` panel. The Katz–Plotkin fraction
is a deposition/placement rule for a particular low-order representation of
vorticity born during \([t^n,t^{n+1}]\). FLOWPanel's `Das` changes the length
of a body-owned transition surface, its influence in the body operator, the
free-wake handoff, and FMM source radii. Those are not demonstrably the same
discretization.

The issue is irrelevant to a one-shot steady boundary-value problem: there is
no physical \(\Delta t\), so \(0.25U_\infty\Delta t\) cannot select a physical
steady wake length. It is important in a genuinely unsteady solve, or in a
pseudo-time march toward steady state.

**Open question for Phase 2:** How should Route B represent the vorticity born
during the most recent physical timestep so that its same-step body influence,
Kelvin bookkeeping, persistent wake storage, and temporal accuracy are
consistent—without replacing `Das` by another empirical release fraction?

Four deterministic candidates are retained:

#### B1 — full TE-anchored constant-strength swept panel

Connect the current TE directly to the first old wake row after convection and
give the whole new strip the accepted \(\gamma^{n+1}\). This is the simplest
fully implicit, first-order rule:

\[
\left(B^{n+1}+W_{\mathrm{swept}}^{n+1}C\right)\mu^{n+1}
=b^{n+1}(\text{old wake}).
\tag{15}
\]

Unlike a lumped vortex, the panel kernel distributes influence over the whole
swept strip, including arbitrarily close to the TE; therefore the
Katz–Plotkin point-placement fraction does not transfer to it directly.
Advantages are exact compatibility with constant-strength `PanelWake`
storage, no empirical length fraction, a linear jump-closure solve, and simple
Kelvin accounting. Risks are first-order temporal bias, near-TE conditioning,
and over-applying the final circulation to material that was born earlier in
the step.

This is the lowest-complexity Route B baseline, not a presumption of adequate
temporal accuracy.

#### B2 — persistent constant-strength subpanels within one timestep

After the B1 timestep study, divide the one physical swept strip into
\(N_s\) streamwise constant-strength subpanels. These are quadrature/storage
subdivisions of one physical timestep, not additional physical timesteps:
provider history, FGS history, monitors, and wake-age bookkeeping still
advance once.

Let \(0=\theta_0<\theta_1<\cdots<\theta_{N_s}=1\) follow the swept trajectory
defined by B4 in (20), with subpanel \(j\) centered at
\(\bar\theta_j=(\theta_{j-1}+\theta_j)/2\). The time-interpolated strength is

\[
\gamma_j
=(1-\bar\theta_j)\gamma^n
+\bar\theta_j\gamma^{n+1}.
\tag{16}
\]

If \(W_j\) is the influence of subpanel \(j\), collect the current and
historical parts:

\[
W_{\rm current}=\sum_{j=1}^{N_s}\bar\theta_j W_j,
\qquad
W_{\rm history}=\sum_{j=1}^{N_s}(1-\bar\theta_j)W_j.
\tag{17}
\]

With jump closure, the same-step system remains linear:

\[
\left(B+W_{\rm current}C\right)\mu^{n+1}
=b_{\rm old}-W_{\rm history}\gamma^n.
\tag{18}
\]

After acceptance, all \(N_s\) subpanels and their interpolated strengths are
stored and convected as old wake. They must not be collapsed into one
constant-strength panel, because that would discard the temporal variation
the experiment is intended to preserve.

B2 has two required substeps:

1. **Geometric null gate:** subdivide the B1 strip but assign every subpanel
   the same \(\gamma^{n+1}\). For a geometrically identical constant-strength
   sheet, internal doublet-panel edges must cancel; the body influence and
   loads should match B1 to the numerical floor. This separates topology and
   kernel bugs from time interpolation.
2. **Interpolated-strength experiment:** use (16), then sweep at least
   \(N_s=2,4,8\). Compare against the B1 timestep-refined limit and check
   convergence with \(N_s\) at fixed physical \(\Delta t\).

This uses existing constant-strength kernels at the cost of more wake
elements, FMM buffering, storage, VTK data, and panel-to-particle handoff
work. A genuinely linear/higher-order panel kernel is a **separate last
resort**, considered only if constant-strength subpanels demonstrate the
needed accuracy but are prohibitively expensive. No higher-order kernel is
part of the first Route B implementation.

#### B3 — Katz–Plotkin-style deterministic lump/release rule

Represent the newborn strip by an element located at

\[
\boldsymbol x_{\rm release}
=\boldsymbol x_{\rm TE}^{n+1}
-k\,\Delta\boldsymbol x_{\rm rel},
\qquad k\in[0.2,0.3],
\tag{19}
\]

where \(\Delta\boldsymbol x_{\rm rel}\) is the local relative TE/wake
convection displacement, not literally \(U_\infty\Delta t\). The latter
vanishes in hover and is not the relevant path.

This is deterministic once \(k\) is selected and is a necessary historical
comparator, with \(k=0.25\) the nominal case. It does **not** eliminate a
placement coefficient; it reclassifies it correctly as low-order temporal
quadrature rather than a physical attached-wake length. It must be swept over
at least \(0.2,0.25,0.3\) and refined in \(\Delta t\). B3 is acceptable as a
benchmark or fallback, not the preferred final claim of a parameter-free
Route B.

#### B4 — time-integrator-derived swept-strip deposition operator (preferred)

Derive the newborn-wake influence from the same discrete time integrator used
to move the body and wake. Let

\[
t_\theta=t^n+\theta\Delta t,\qquad 0\le\theta\le1,
\]

and let \(\boldsymbol X_e^{n+1}(\theta)\) be the position at \(t^{n+1}\) of
material shed from edge \(e\) at birth time \(t_\theta\):

\[
\boldsymbol X_e^{n+1}(\theta)
=\boldsymbol X_{\mathrm{TE},e}(t_\theta)
+\int_{t_\theta}^{t^{n+1}}
\boldsymbol u_{\mathrm w}(\boldsymbol X_e(s),s)\,ds.
\tag{20}
\]

Approximate (20) with the same backward-Euler, trapezoidal, Runge–Kutta, or
other stage trajectory used by wake convection. Interpolate the circulation
over the birth interval with a basis consistent with that integrator:

\[
\boldsymbol\gamma(\theta)
=\sum_{q=0}^{p}\ell_q(\theta)
\boldsymbol\gamma^{\,n+1-q}.
\tag{21}
\]

All histories \(\gamma^n,\gamma^{n-1},\ldots\) are frozen; only
\(\gamma^{n+1}\) is a current trial variable. The newborn-strip influence is
the space–birth-time integral

\[
\mathcal W_{\rm birth}[\gamma]
=\int_{\mathrm{span}}\int_0^1
\mathcal K\!\left(\boldsymbol X(s,\theta)\right)
\gamma(s,\theta)\,J(s,\theta)\,d\theta\,ds.
\tag{22}
\]

After quadrature and collection by time level,

\[
\mathcal W_{\rm birth}[\gamma]
=W_0\gamma^{n+1}
+W_1\gamma^n+\cdots+W_p\gamma^{n+1-p}.
\tag{23}
\]

Only \(W_0\gamma^{n+1}\) belongs in the current operator. All older terms are
known and move to the right-hand side. With jump closure, the body system
remains linear:

\[
\left(B^{n+1}+W_0C\right)\mu^{n+1}
=b^{n+1}_{\rm old}
-\sum_{q=1}^{p}W_q\gamma^{n+1-q}.
\tag{24}
\]

With pressure closure, the same \(W_0\) enters each nonlinear trial while all
old terms stay frozen. Trial evaluation never sheds or commits a row.

The essential deterministic choices are then the physical time integrator,
its interpolation/dense-output rule, and a converged numerical quadrature—not
an empirical fraction \(k\). Required invariants are:

1. **constant-circulation reproduction:** if all time levels have the same
   \(\gamma\), (23) equals the influence of one continuous constant-strength
   swept sheet;
2. **partition of unity:** \(\sum_q\ell_q(\theta)=1\), so no circulation is
   created by interpolation;
3. **interface consistency:** the downstream limit of the newborn strip
   matches the newest old-wake strength, and its TE limit matches the accepted
   current strength;
4. **one-time commit:** only the accepted \(\gamma^{n+1}\) and wake row enter
   persistent history;
5. **same-path geometry:** body motion, TE motion, and wake convection use the
   same stage positions as (20);
6. **no endpoint tuning:** near-TE singular behavior is treated by the
   analytic panel integral, singularity subtraction, or converged interior
   quadrature, not by moving a node until loads look favorable; and
7. **startup order:** missing histories use the same documented order
   reduction as the pressure-history and wake-convection integrators.

For a first-order endpoint scheme—including the current explicit-Euler wake
trajectory combined with a same-step implicit strength solve—B4 reduces
naturally to the B1-type constant new strip. For a second-order method, the
first candidate is linear
interpolation between frozen \(\gamma^n\) at the downstream birth boundary and
trial \(\gamma^{n+1}\) at the TE, with geometry from the same second-order
trajectory. Runge–Kutta methods should use their stage geometry and a
consistent dense-output/interpolation rule.

The main implementation question is whether (22)–(24) can be evaluated using
ephemeral quadrature/influence operators while retaining the existing
constant-strength persistent wake representation without losing the claimed
order at the next step. The selected first answer is to avoid collapse and
persist the B2 constant-strength subpanels. An ephemeral collapsed operator is
deferred. A true higher-order panel kernel is justified only if the persistent
subpanel result succeeds numerically but its cost is unacceptable.

### Candidate priority and verification

The user-selected research order is:

1. **B1 first:** one persistent TE-anchored Euler panel per physical timestep.
   Run at least three successively halved timesteps; report observed order and
   a Richardson-extrapolated limit for circulation, pressure jump, and loads.
2. **B2 second:** first pass the equal-strength geometric null gate, then use
   the linear birth-time interpolation (16) with \(N_s=2,4,8\) persistent
   constant-strength subpanels per physical timestep.
3. **Compare B1 and B2:** determine whether subpanels approach the B1
   timestep-refined limit at materially lower physical-timestep cost, and
   whether either changes the inferred steady/periodic result.
4. **B3 comparator:** retain the Katz–Plotkin \(k=0.2,0.25,0.3\) cases as a
   historical low-order comparison.
5. **Higher-order panel kernel last:** consider it only after the subpanel
   experiment establishes a numerical benefit that cannot be obtained at
   acceptable cost with existing kernels.

B4 remains the organizing derivation: B1 is its first-order endpoint
discretization, while B2 is its first persistent birth-time quadrature.
Acceptance requires:

- Wagner/finite-wing indicial response under timestep refinement;
- convergence of circulation, pressure jump, and loads rather than agreement
  at one chosen \(k\);
- comparison with B3 at \(k=0.2,0.25,0.3\);
- exact constant-\(\gamma\) and Kelvin bookkeeping checks;
- zero-angle symmetry and impulsive-start startup checks; and
- rotor-relative displacement tests where \(U_\infty=0\).

Until these tests are complete, Route B is accurately described as
“independent of rigid `Das`,” not as “free of newborn-wake discretization.”

### Phase 1 representation decision

**Require Phase 2 to make wake attachment and Kutta closure independently
configurable if the solver contracts permit it.** The intended verification
matrix is all four combinations above.

The lowest-risk implementation sequence is:

1. preserve A/jump bitwise as the default;
2. implement B/jump as the clean `Das` experiment without nonlinear Kutta
   changes;
3. implement A/pressure using the existing affine correction machinery; and
4. compose B/pressure through the same pressure driver after both independent
   axes pass their mechanical tests.

This sequencing is not permission to implement before Phase 2 approval. It
prevents a favorable B/pressure result from being ambiguously attributed to
either changed wake geometry or changed Kutta enforcement.

## 6. Nonlinear strategy comparison

One residual evaluation means: apply a trial correction, solve the body,
reconstruct the complete surface velocity from a clean snapshot, and evaluate
pressure. The dominant cost is normally the body solve.

| Method | Body solves per outer update | Memory | Spanwise coupling | Assessment |
| --- | ---: | ---: | --- | --- |
| Relaxed fixed point | 1 | \(O(M)\) | Only through the chosen update map | Cheap, but no natural pressure-to-correction map and weak robustness near strong coupling |
| Independent/local secant | 1 after startup | \(O(M)\) | Diagonal approximation only | Useful preconditioner; not acceptable as the main coupled method |
| Anderson acceleration | 1 | \(O(mM)\) | Captures recent low-rank coupling | Good if a contractive fixed-point map is found; no such map is presently justified |
| Full Broyden | 1 | \(O(M^2)\) | Dense global approximation | Best first balance of cost and coupling; needs damping, scaling, and restart safeguards |
| Finite-difference Newton | \(M+1\) per Jacobian rebuild | \(O(M^2)\) | Full | Strong small-\(M\) oracle, too expensive as the default |
| Finite-difference Newton–Krylov | \(1+\) Krylov Jv trials, plus line search | \(O(mM)\) | Full | More robust fallback for larger \(M\); body-solve count can be high and differencing noise must be controlled |

### Primary recommendation: safeguarded global Broyden

Use a scaled, damped Broyden method with:

- the full correction and residual vectors, never independent edge solves;
- dimensionless residual/correction variables;
- a conservative initial inverse-Jacobian scale or diagonal secant
  preconditioner;
- backtracking on the normalized residual norm;
- rejection of non-finite trials;
- minimum-step and maximum-correction safeguards;
- restart after repeated rejected steps or a badly conditioned update; and
- convergence only when pressure, correction-change, and inner-solve gates all
  pass.

After initialization it needs one new body solve per accepted outer update and
learns dense spanwise coupling. The \(O(M^2)\) storage is acceptable for the
first single-body paired-edge target.

### Fallback recommendation: finite-difference Newton–Krylov

Use matrix-free

\[
J(\boldsymbol c)\boldsymbol v
\approx
\frac{\boldsymbol r(\boldsymbol c+\epsilon\boldsymbol v)
-\boldsymbol r(\boldsymbol c)}{\epsilon}
\]

with a span-local/diagonal preconditioner, forcing terms consistent with the
inner body-solve accuracy, and the same globalization and restoration rules.
This is the fallback when Broyden stagnates or repeatedly restarts. A full
finite-difference Jacobian remains a small-case verification oracle, not a
production default.

## 7. Pressure-law-independent trial contract

The nonlinear formulation owns correction application, body solves, velocity
reconstruction, edge pairing, residual formation, convergence, restoration,
and final commit. A pressure provider owns only repeatable pressure evaluation
from the supplied complete trial state.

Required semantics:

1. Repeated evaluation of the same trial state returns the same pressure.
2. `evaluate` does not mutate physical-time history, body/wake history,
   monitor context, or solver history.
3. Provider requirements—velocity, potential, potential history, gradient,
   acceleration, and so on—are declared before trials begin.
4. A separate `commit` operation, if required, occurs exactly once for the
   accepted state.
5. A common pressure gauge may vary only in ways that cancel from paired
   differences; otherwise the provider must expose a consistent gauge.

Steady Bernoulli is the first concrete provider and requires only the complete
trial surface velocity, density, and a reference speed/pressure convention.
Density changes scale every residual uniformly; normalized convergence and
the accepted correction should therefore be density-invariant.

For unsteady Bernoulli, every trial at physical step \(n+1\) must use the same
frozen \(\phi^n,\phi^{n-1},\Delta t_n,\Delta t_{n-1}\) and moving-point data.
Trial evaluations must not shift the BDF history. After nonlinear convergence,
the accepted potential is committed once. If convergence fails, history is
unchanged. `PressureLaplace` and unsteady Bernoulli adapters are architectural
extensions, not Phase 3 deliverables.

## 8. State and failure invariants implied by the theory

Before the first trial, snapshot all state that is prescribed at the physical
timestep: freestream/kinematic velocity, old free-wake influence, external
potential, body strengths, velocity fields, selected attachment geometry,
live-panel or affine correction state, provider history, and body-solver
history.

Each trial must:

1. restore the prescribed snapshot;
2. apply the trial circulation coordinates for the selected closure;
3. solve (11) without letting an affine correction or stale live-panel
   strength contaminate a linear operator product;
4. reconstruct (3) once, with no accumulation;
5. evaluate (4) without history mutation; and
6. form (5).

`FGSSolver` projection/history is physical-time history. Nonlinear trials must
not call `save_solution!`; projection must use the same pre-trial history for
every trial. Exactly one accepted solution may advance that history.

On nonconvergence, an invalid inner solve, non-finite pressure, or rejected
configuration, restore the last accepted state and throw an informative error
by default. A permissive fallback must be explicit; silently retaining the
last nonlinear trial is forbidden.

## 9. Rejected or deferred alternatives

- **Postprocess pressure averaging:** rejected; it does not close the coupled
  body/wake problem.
- **Change the live wake strength after the body solve:** rejected by (11);
  jump closure must include \(W_RC\) in the same solve.
- **Independent edge solves:** rejected because the Jacobian is dense.
- **Match velocity magnitudes directly:** rejected as the public residual;
  equivalent only for steady Bernoulli and not pressure-law independent.
- **Suppress the attached panel without replacing the TE-to-wake sheet:**
  rejected; it leaves a geometric and circulation gap.
- **Treat TE attachment and Kutta closure as one option:** rejected; it makes
  the `Das` and pressure-continuity effects impossible to attribute.
- **Advance unsteady pressure or FGS history on every trial:** rejected; outer
  iterations are not physical timesteps.
- **Exact affine response matrix as the default:** deferred. It is valuable
  for derivative verification and perhaps repeated steady solves, but costs
  \(M\) body solves to construct and is provider-specific beyond steady
  Bernoulli.
- **Claim that pressure continuity removes `Das`:** rejected for Route A.
  Only Route B removes the independent rigid offset.
- **Claim that Route B is resolution-free:** rejected. It replaces `Das` with
  a physically convected first-panel geometry whose timestep/wake-resolution
  convergence must be demonstrated.

## 10. Approved Phase 1 decision record and acceptance gate

The user explicitly approved the following decisions on 2026-07-29:

1. wake attachment and Kutta closure are independent configuration axes;
2. the target matrix contains A/jump, A/pressure, B/jump, and B/pressure,
   staged in that order after the preserved default;
3. jump closure uses the deterministic relation \(\gamma=C\mu\) inside the
   same solve, while Route A pressure coordinates retain
   \(\gamma=C\mu-c\);
4. Route B owns a TE-anchored live panel, freezes old wake strengths, and has
   explicit startup and one-time commit semantics;
5. newborn-wake deposition follows the selected experiment order: B1
   one-row Euler timestep refinement first; B2 persistent constant-strength
   subpanels second, including the equal-strength null gate and
   \(N_s=2,4,8\) interpolation study; B3 remains the Katz–Plotkin comparator;
   a higher-order panel kernel is last resort;
6. the pressure condition matches paired panel-centroid pressures, with
   sharp-edge convergence deferred to verification;
7. safeguarded global Broyden is primary and finite-difference
   Newton–Krylov is fallback;
8. steady Bernoulli is the only built-in Phase 3 pressure provider;
9. all trial and solver histories are frozen and committed once; and
10. attribution uses the four-combination matrix: Route A is never advertised
   as `Das`-free, while Route B is tested for timestep/wake convergence.

**Acceptance gate:** [x] Theory and algorithm decisions explicitly approved;
Phase 2 may begin, but this approval does not authorize implementation.

---

## Phase 1 progress and approval log

- **2026-07-29 — Phase 1 derivation completed by agent.** Source paths reviewed:
  `src/FLOWPanel_formulation.jl`, `src/FLOWPanel_solver.jl`,
  `src/FLOWPanel_liftingbody.jl`, `src/FLOWPanel_elements_fmm.jl`,
  `src/FLOWPanel_postprocess.jl`, `src/FLOWPanel_simulate.jl`,
  `src/FLOWPanel_wake.jl`, `docs/wake_solve_schemes.md`,
  `BRAINSTORM/007_pressure_monitor_reliability.md`, and
  `BRAINSTORM/014_first_wake_row_offset_selection.md`. No source or test code
  changed. Awaiting user-agent discussion and explicit Phase 1 approval.
- **2026-07-29 — Phase 1 revised after user discussion.** Separated wake
  attachment (A: rigid `Das`; B: TE-anchored) from Kutta closure (jump:
  \(\gamma=C\mu\); pressure continuity). Added the four-combination
  attribution matrix, steady and unsteady Route B equations, deterministic
  shedding/Kelvin semantics, startup requirements, and the requirement that
  Phase 2 make both axes independently configurable if mechanically feasible.
  At this point, Phase 1 was still awaiting explicit approval.
- **2026-07-29 — Newborn-wake deposition question logged after user
  discussion.** Recorded the Katz–Plotkin \(0.2\)–\(0.3\) placement as a
  low-order unsteady deposition rule rather than universal support for
  FLOWPanel's rigid `Das`. Added four deterministic Route B candidates:
  full swept panel (B1), higher-order/subdivided strip as last resort (B2),
  Katz–Plotkin release comparator (B3), and the preferred
  time-integrator-derived swept-strip operator (B4). B4 now specifies
  birth-time geometry, circulation interpolation, operator splitting by time
  level, linear jump-closure coupling, pressure-trial freezing, conservation
  invariants, persistent-storage caveat, and verification gates. The question
  remained open for Phase 2; at this point, Phase 1 was still awaiting
  explicit approval.
- **2026-07-29 — Newborn-wake experiment order selected by user.** B1 now
  comes first as one persistent TE-anchored Euler panel per physical timestep,
  with timestep refinement and Richardson reporting. B2 follows with
  persistent constant-strength subpanels: an equal-strength internal-edge
  cancellation/null test, then linearly interpolated strengths from
  \(\gamma^n\) downstream to \(\gamma^{n+1}\) at the TE for
  \(N_s=2,4,8\). The subpanels advance physical histories only once and are
  retained rather than collapsed. A new higher-order panel kernel is deferred
  to last resort. At this point, Phase 1 was still awaiting explicit approval.
- **2026-07-29 — Phase 1 explicitly approved by user.** The ten decisions in
  the approved decision record are accepted as the architecture handoff.
  Phase 2 is ready to begin; no implementation is authorized by this event.
