# Wake coupling in the Dirichlet solve

## Purpose and current verdict

FLOWPanel has two wake representations that can be compared in a nominally steady
calculation:

- a semi-infinite attached wake whose strength is coupled to the unknown body
  doublets; and
- a finite attached wake followed by a time-marched free wake.

The marched calculation transfers the free wake's velocity to the body and therefore
changes the body source strengths. It does not transfer the free wake's scalar
potential. The two routes can be mathematically equivalent, but that conclusion has
three separate parts: a continuum Green representation, its finite-panel realization,
and equality of the composite finite/free wake with the semi-infinite wake.

The derivation below identifies a Kutta-trace term that could contaminate the
trailing-edge circulation at finite resolution. Its dependence on the first attached
wake length is consistent with the observed circulation deficit, but this is a
hypothesis, not a demonstrated cause. The numerical checks in the final section are
needed to distinguish it from discretization and wake-geometry errors.

## Glossary and conventions

The body occupies a closed domain with boundary $\partial\Omega$ and outward unit
normal $\boldsymbol n$. All boundary limits below use FLOWPanel's kernel signs and
the interior potential limit.

- $\psi$: potential produced by the known free wake.
- $q$: body trace of that potential,
  $$
  q=\left.\psi\right|_{\partial\Omega}.
  $$
- $\sigma$: source density associated with the free-wake normal velocity,
  $$
  \sigma=-\frac{\partial\psi}{\partial n}.
  $$
- $\mathcal S$: body single-layer operator, evaluated with its interior boundary
  trace when its output appears in a solve equation.
- $\mathcal B$: body double-layer interior trace operator. In the discrete code it
  includes the explicit positive half self-potential limit.
- $\mathcal W$: potential influence of the finite attached wake at body control
  points.
- $\mathcal C$: Kutta map from body doublet strengths to attached-wake strengths.
- $\boldsymbol\mu$: vector of body doublet strengths.
- $\mathcal G$: finite attached-wake solve operator,
  $$
  \mathcal G=\mathcal B+\mathcal W\mathcal C.
  $$
- $\mathcal D_{\mathrm e}$, $\mathcal S_{\mathrm e}$, and
  $\mathcal W_{\mathrm e}$: the corresponding fields evaluated at exterior
  observation points.
- $\mathcal O$: an exterior observation operator, such as potential or velocity
  evaluation on a selected exterior set.

Symbols in equations denote mathematical operators and fields. Julia identifiers are
used only in the implementation discussion.

## What the current code solves

### Exact-centroid boundary operator

Body control points are the panel centroids; they are not displaced to geometrically
interior points. `calc_controlpoints!` in `src/FLOWPanel_abstractbody.jl` computes the
triangle centroids exactly. For a self pair, `_self_limit` in
`src/FLOWPanel_elements_fmm.jl` supplies the selected one-sided limit. In the
Dirichlet potential operator the body doublet self contribution is
$+\tfrac12$ times the panel strength. `_G!` in `src/FLOWPanel_solver.jl`
therefore assembles the interior-limit matrix directly and does not add a second jump
term.

For a lifting `RigidWakeBody`, the same matrix column also includes the attached wake
panel tied to that column by the trailing-edge map. Thus the finite-wake matrix is
$\mathcal G$, not the body-only operator $\mathcal B$.

### Source strengths and free-wake transfer

For a Dirichlet source/doublet body, `set_strengths!` sets the source strength from the
normal component of the apparent velocity. In mathematical notation,

$$
\sigma_{0}=-\boldsymbol u_{0}\mathbin{\boldsymbol\cdot}\boldsymbol n,
\qquad
\sigma=-\boldsymbol u_{\mathrm f}\mathbin{\boldsymbol\cdot}\boldsymbol n,
$$

where $\boldsymbol u_{0}$ contains the non-free-wake contributions and
$\boldsymbol u_{\mathrm f}=\boldsymbol\nabla\psi$ is the free-wake velocity.

Production time marching evaluates wake influence in `_steady_aerodynamics!` in
`src/FLOWPanel_simulate.jl` with scalar-potential output disabled and velocity output
enabled. Consequently, the free wake contributes $\sigma$, but its trace $q$ is
not assembled into the Dirichlet right-hand side.

There is a solver-path distinction worth keeping explicit:

- A single-body Dirichlet `solve!` saves `body.potential`, clears it, assembles the
  body's source potential, solves, and restores the saved values. Preassembled
  external potential is therefore discarded for that solve.
- The genuinely coupled `KrylovCoupled` and `BackslashCoupled` paths save external
  potential and add it back after assembling body influence, so they can preserve an
  explicitly assembled external potential.
- The fixed-point tuple-of-solvers path transfers cross-body velocity only and invokes
  the single-body solve for each body.

This distinction does not alter the production free-wake route: no branch of
`_steady_aerodynamics!` assembles free-wake scalar potential on the body.

### The two wake systems

Let $\boldsymbol\sigma_{0}$ be the source strengths without the free-wake velocity.
The semi-infinite attached-wake system is

$$
\mathcal G_{\infty}\boldsymbol\mu_{\infty}
=-\mathcal S\boldsymbol\sigma_{0}.
$$

For the finite attached wake plus a known free wake, an explicit-potential formulation
would solve

$$
\mathcal G\boldsymbol\mu_{E}
=-\mathcal S\boldsymbol\sigma_{0}-q.
$$

The production velocity-through-sources formulation instead solves

$$
\mathcal G\boldsymbol\mu_{V}
=-\mathcal S\boldsymbol\sigma_{0}-\mathcal S\sigma.
$$

The first comparison is between the latter two finite-wake systems. Comparison with
the semi-infinite system requires additional wake-representation conditions.

In implementation terms, `semiinfinite_wake=true` makes `_induced_wake` extend the
attached panel from the trailing edge along `Das` to infinity. With
`semiinfinite_wake=false`, `_induced_wake` ends that panel at the `Das` endpoints and
`update_TE!` places the upstream edge of the first free-wake row at the same handoff.
The attached panel always remains coupled to the current unknown through the Kutta
map; only the downstream free wake is treated as a known time-marched system.

## Discrete diagnostic construction: Task 1 semi-infinite baseline

This section is the implementation gate for Task 1 of
`debug/dirichlet_solve/dirichlet_solve.jl`.  It records the exact discrete
system used for the reference calculation; later tasks must add their own
construction before their solver code is implemented.

### Prescribed fields and unknowns

The prescribed body is the closed, watertight AR=4 rectangular NACA 0015 mesh
constructed by `build_pitching_wing_body` in `examples/pitching_wing.jl`, with
`n_airfoil=161`, `n_span=13`, and `n_endcap=9`.  The freestream is
$U_\infty=330.2\ \mathrm{ft/s}$ at $\alpha=3.94^\circ$.  Its direction also
sets every semi-infinite attached-wake ray.  There is no finite attached panel
and no free-wake system in this task.

The prescribed source coefficient on body panel $i$ is

$$
(\sigma_0)_i=-\boldsymbol U_\infty\mathbin{\boldsymbol\cdot}\boldsymbol n_i.
$$

The only unknown is the constant body-doublet coefficient vector
$\boldsymbol\mu_\infty$.  The paired shedding-edge circulation is a derived
quantity,

$$
\boldsymbol\gamma_\infty=\mathcal C\boldsymbol\mu_\infty.
$$

For shedding edge $e$, FLOWPanel's `shedding` table makes row $e$ of
$\mathcal C$ equal to $+1$ at the upper panel and $-1$ at the paired lower
panel.  All rows in this diagnostic must be paired, and therefore
$\mathcal C\boldsymbol 1=0$ exactly.

### Exact linear system and operator assembly

The diagnostic solves

$$
\boxed{
\mathcal G_\infty\boldsymbol\mu_\infty
=-\mathcal S\boldsymbol\sigma_0,
\qquad
\mathcal G_\infty=\mathcal B+\mathcal W_\infty\mathcal C.
}
$$

`Backslash(body)` assembles this matrix by unit-activating the doublet
coefficient of one body panel at a time and evaluating scalar potential at all
body centroids.  The body part is $\mathcal B$.  For a shedding panel,
`_induced_wake` adds the potential of the semi-infinite panel/ray associated
with that same column.  The upper and lower column orientations implement the
two signs in $\mathcal C$, so the assembled wake block is
$\mathcal W_\infty\mathcal C$ without constructing either factor separately.
The right-hand side is obtained by first assigning $\sigma_0$ from the
freestream and then directly evaluating the source-panel potential
$\mathcal S\boldsymbol\sigma_0$ at the same centroids.

The sign and boundary-limit conventions are those defined above and in
`docs/dirichlet_potential_theory.md`: source strength is negative outward
normal velocity, doublet strength is
$\mu=\phi_{\rm inside}-\phi_{\rm outside}$, and the body-doublet self term is
the interior $+\mu_i/2$ limit already returned by the kernel.  No additional
diagonal jump is added by the diagnostic.

No gauge constraint is imposed on this solve.  The physical zero-at-infinity
convention of the source/doublet and semi-infinite-wake kernels fixes the
represented exterior potential, and the directly assembled
$\mathcal G_\infty$ system is solved as written.  The area-weighted trace gauge
needed by the later Green-reconstruction task is not part of this baseline.

The finite-wake correction operator $\mathcal Wc$ is also not present here:
$c$ is a free-wake trace correction, while Task 1 has no free wake.  For this
task, $\mathcal B$, $\mathcal C$, and $\mathcal W_\infty\mathcal C$ are defined
as above, $\mathcal G_\infty$ is their assembled sum, and $\mathcal Wc=0$.

### Invariants, outputs, and limitations

Before accepting the baseline, the diagnostic checks that the mesh is
watertight, every shedding edge has a lower partner, the shedding map indices
are valid, and $\mathcal C\boldsymbol1=0$.  It reports the direct relative
linear residual

$$
r_{\rm rel}=
\frac{\lVert\mathcal G_\infty\boldsymbol\mu_\infty
-\boldsymbol b_\infty\rVert_2}
{\max(\lVert\boldsymbol b_\infty\rVert_2,\epsilon)},
\qquad
\boldsymbol b_\infty=-\mathcal S\boldsymbol\sigma_0,
$$

requires all coefficients and sampled exterior velocities to be finite, and
records $\boldsymbol\gamma_\infty$, steady $C_L$, and a common set of exterior
probe velocities for later relative comparisons.  Task 1 has no transition to
a free wake, so transition/free-wake handoff checks are explicitly not
applicable.

This is a constant-strength, centroid-collocated discretization.  The
semi-infinite wake is a prescribed flat ray construction, not a developed or
rolled-up wake, and the pressure-derived $C_L$ retains the surface-gradient and
force-integration errors of the chosen mesh.  Consequently it is a numerical
reference for Tasks 2--5, not an exact finite-wing solution.

## Discrete diagnostic construction: Task 2 finite velocity-through-sources wake

Task 2 evaluates the finite-wake system already used by FLOWPanel's production
time marcher.  It does not add the known free-wake scalar-potential trace to the
right-hand side.  Instead, the prescribed free wake supplies velocity
$\boldsymbol u_f$ at every body centroid, and the normal component of that
velocity is converted to an additional body-source coefficient,

$$
\sigma_f=-\boldsymbol u_f\mathbin{\boldsymbol\cdot}\boldsymbol n.
$$

With $\mathcal G_\Delta=\mathcal B+\mathcal W_\Delta\mathcal C$, the exact
discrete system tested here is therefore

$$
\boxed{
\mathcal G_\Delta\boldsymbol\mu_V
=-\mathcal S(\boldsymbol\sigma_0+\boldsymbol\sigma_f).
}
$$

The signs are the same as in Task 1: $\sigma=-\boldsymbol u\cdot\boldsymbol n$,
$\mu=\phi_{\rm inside}-\phi_{\rm outside}$, and each row of the paired Kutta
map is upper-panel doublet minus lower-panel doublet.  The body control points
remain exact centroids and the $+\mu_i/2$ interior self limit is already part of
$\mathcal B$.

### Finite attached/free-wake geometry and unknowns

The body mesh, flow, and exterior probes are exactly those of Task 1.  The
semi-infinite ray is replaced by a finite attached transition panel of length
$U_\infty\Delta t=0.5c$, aligned with the Task 1 freestream ray.  Its upstream
edge is the paired trailing edge and its downstream edge is the first stored
free-wake node row.  The transition-panel circulation is not an independent
unknown: it is $\mathcal C\boldsymbol\mu_V$.  Thus its contribution to the
matrix is $\mathcal W_\Delta\mathcal C$.

Every stored free-wake panel row begins at that downstream transition edge.
Its node geometry and circulation are prescribed during a reported solve and
do not enter the unknown vector.  FLOWPanel stores the row with the contiguous
wake orientation used by `shed_wake!`; with that orientation its scalar
constant-doublet coefficient has the same numerical value as the paired Kutta
circulation $\mu_{\rm upper}-\mu_{\rm lower}$.  The diagnostic verifies this
sign at the transition-to-first-row handoff.

For the primary fabricated cases, every free-wake row is also $0.5c$ long and
collinear with the Task 1 ray.  All rows carry the freshly recomputed Task 1 circulation
$\boldsymbol\gamma_\infty=\mathcal C\boldsymbol\mu_\infty$.  Composite lengths
$L/c=1,2,4,8,16,32,64$ therefore contain one attached half-chord panel and
$2L/c-1$ prescribed free rows.  No downstream final filament is included: the
finite panel sheet is closed by its own finite constant-doublet-panel edges.
Moving that downstream closure farther away makes its field at the body tend
toward the semi-infinite construction.  At any finite length, however, the
closure field, the velocity-through-sources replacement, and the fact that the
transition circulation remains coupled while the free rows are prescribed mean
that equality with Task 1 is neither imposed nor expected.

A second fabricated sequence isolates transition-panel length.  It shortens
the attached panel to $|D_a|=0.05c$ but retains the same $0.5c$ prescribed
free-row spacing and the same row counts $1,3,7,15,31,63,127$.  Its free-wake
lengths are therefore $0.5c,1.5c,3.5c,7.5c,15.5c,31.5c,63.5c$, and its actual
trailing-edge-to-closure extents are $0.55c,1.55c,3.55c,7.55c,15.55c,31.55c,
63.55c$.  Reporting both the free-wake length and total composite extent avoids
relabeling this sequence as though the attached panel were still half a chord.
The free-row circulation, immutable-state checks, convergence criterion, and
Task 1 comparisons are otherwise unchanged.

For the production case, the short $|D_a|=0.05c$ transition is followed by the
default induced-velocity `PanelWake`, with a $0.5c/U_\infty$ time step.  Constant
incidence and freestream are used, with no pitching or body motion.  Each step
follows the production order:
update the transition/free-wake edge, transfer free-wake velocity to the body,
solve the finite Dirichlet system, evaluate lift, convect the wake, and shed the
new circulation.  The calculation is first assessed at $40c/U_\infty$ and may
be extended in $20c/U_\infty$ blocks through the fixed maximum
$80c/U_\infty$.

### Numerical gates and immutable reporting solves

The primary fabricated sequence is called length-converged only when the last
two values, at $32c$ and $64c$, satisfy

$$
\frac{|C_L(64c)-C_L(32c)|}
{\max(|C_L(64c)|,|C_L(32c)|,\epsilon)}\le10^{-3}.
$$

The $|D_a|=0.05c$ sequence applies the same relative criterion to its final two
row-count cases, whose actual composite extents are $31.55c$ and $63.55c$.

Production settling is evaluated over the latest $10c/U_\infty$ window.  Both
the peak-to-peak lift divided by the magnitude of the window mean and the
least-squares linear change across the window divided by that same magnitude
must not exceed $0.5\%$.  Failure at the fixed $80c/U_\infty$ horizon is
reported as nonconvergence; no favorable intermediate step is selected.

Before every reported finite-wake solve, the complete master body/wake state is
deep-copied.  Exact snapshots of body nodes and connectivity, free-wake nodes
and strengths, active row count, and wake options are compared after the solve.
Only the copied body's source/doublet arrays may change; the attached
transition circulation is then recomputed as $\mathcal C\boldsymbol\mu_V$.
For each case the diagnostic records the direct linear residual, $C_L$, that
transition circulation, its mismatch with the first prescribed free row, and
the relative $L_2$ difference at Task 1's sixteen exterior velocity probes.
These measurements report signed and relative differences from Task 1 without
using closeness to Task 1 as an acceptance tolerance.

## Discrete diagnostic construction: Task 3 direct fixed-wake potential

Task 3 retains the complete finite body/attached-transition operator from Task
2 but replaces its velocity-through-sources approximation for the prescribed
free wake.  Let $q_f$ be the scalar potential induced directly by the active
constant-doublet panels of the frozen `PanelWake`, evaluated at the exact body
centroids with `DirectBackend`.  The diagnostic system is

$$
\boxed{
\mathcal G_\Delta\boldsymbol\mu_E
=-\mathcal S\boldsymbol\sigma_0-\boldsymbol q_f,
\qquad
\boldsymbol\sigma_0=-\boldsymbol U_\infty\mathbin{\boldsymbol\cdot}\boldsymbol n.
}
$$

No wake-induced velocity is converted into body-source coefficients and no
mean is subtracted from $q_f$.  The finite panel kernel's zero-at-infinity
potential convention fixes the gauge.  The body source potential and $q_f$ are
assembled separately before the existing finite-body LU is applied manually;
the ordinary single-body `solve!` route is deliberately not used because it
clears preassembled external potential.

### Frozen single shots and active sources

The two fabricated flat sequences and the selected settled rolled geometry are
constructed by the deterministic Task 2 paths.  Task 3 verifies the selected
step, body and transition geometry, wake options, nodes, row count, and initial
strength hashes before solving.  Every single-shot solve preserves the
prescribed wake strengths bit-exactly.  The active scalar source tuple is
explicitly `(wake,)`: it contains the `PanelWake` doublet panels only.  The
rolled Task 2 state retains its `include_final_filament=true` storage option,
but it has not overflowed, so the final-filament wrapper has zero active
elements; the diagnostic additionally asserts that no active source is
vector-potential-only.

For each flat length and the settled rolled state, Task 3 records
$\boldsymbol q_f$, $\mathcal C\boldsymbol q_f$, the direct linear residual,
transition circulation, first-row handoff mismatch, steady lift, and relative
velocity differences at the common exterior probes against both Tasks 1 and 2.
The same panel-only source selection is used to reconstruct the surface and
exterior velocities used by those diagnostics.

### Fixed-geometry wake-strength projection

Only the three terminal geometries are iterated.  Geometry, active row count,
wake options, and inactive storage remain fixed.  During inner solve $k$, every
active wake-panel strength $\eta_{r,e}^{(k)}$ is frozen.  Between solves it is
updated toward the latest transition circulation,

$$
\eta_{r,e}^{(k+1)}=(1-\omega_k)\eta_{r,e}^{(k)}
 +\omega_k(\mathcal C\boldsymbol\mu_E^{(k)})_e.
$$

The iteration starts at $\omega=1$.  If the maximum rowwise strength defect
fails to decrease by at least 5% on two consecutive iterations, $\omega$ is
halved, with a floor of $1/16$.  Convergence requires both the maximum rowwise
strength defect and relative lift change to be at most $10^{-8}$ for three
consecutive iterations.  The fixed limit is 200 iterations; reaching it is
reported explicitly as nonconvergence.

### Augmented fixed-point oracle

For a fixed terminal geometry, define the inexpensive strip-potential matrix
$P$: column $e$ is the centroid potential generated when every active panel row
in wake strip $e$ has unit strength and all other strips have zero strength.
At the fixed point, $\boldsymbol q_f=P\mathcal C\boldsymbol\mu_*$, so the direct
oracle is

$$
\boxed{
(\mathcal G_\Delta+P\mathcal C)\boldsymbol\mu_*
=-\mathcal S\boldsymbol\sigma_0.
}
$$

The diagnostic evaluates this system through the already-factored
$\mathcal G_\Delta$ and the 13-column Woodbury reduction, avoiding a second
dense body matrix.  It nevertheless evaluates the displayed augmented residual
directly and requires every inner and oracle relative linear residual to be
finite and no larger than $10^{-10}$.  Agreement of the converged circulation
and lift with this independent fixed-point solve validates the iteration.

## Exact Green representation

Assume for the moment that each relevant body interior is simply connected and that
the free-wake vorticity and singular support are disjoint from its closure. Then the
wake velocity has a single-valued harmonic potential $\psi$ in that interior, and
that interior potential is all the construction requires. No single-valued potential
need exist throughout the exterior wake domain. This condition can therefore hold
for vortex-particle wakes and rolled-up filament wakes as long as their singular or
vorticity-bearing support does not overlap the closed body. Green's representation
using the trace $q$ and the density
$\sigma=-\partial\psi/\partial n$ gives, at an interior point,

$$
\psi=\mathcal S\sigma+\mathcal D q.
$$

Taking the interior boundary limit, and defining $\mathcal B$ as the complete
interior trace of the double layer, gives

$$
q=\mathcal S\sigma+\mathcal Bq.
$$

Therefore

$$
\boxed{\mathcal S\sigma=(I-\mathcal B)q.}
$$

This is the required identity. The single layer alone does not generally reproduce
the negative of an arbitrary harmonic wake potential. The missing body double layer
is essential.

Because the vorticity and singular support that generate $\psi$ lie outside the body,
the body-layer representation that reconstructs $\psi$ inside has zero continuation
into the exterior component containing infinity:

$$
\mathcal S_{\mathrm e}\sigma+\mathcal D_{\mathrm e}q=0.
$$

This exterior-invisibility statement is the geometric content behind the
velocity-through-sources substitution.

## Finite attached-wake equivalence

Subtract the explicit-potential system from the velocity-through-sources system:

$$
\mathcal G(\boldsymbol\mu_{V}-\boldsymbol\mu_{E})
=q-\mathcal S\sigma.
$$

Using the Green identity,

$$
\mathcal G(\boldsymbol\mu_{V}-\boldsymbol\mu_{E})
=\mathcal Bq.
$$

Since

$$
\mathcal Bq=\mathcal Gq-\mathcal W\mathcal Cq,
$$

an invertible solve gives

$$
\boxed{
\boldsymbol\mu_{V}-\boldsymbol\mu_{E}
=q-\mathcal G^{-1}\mathcal W\mathcal Cq.
}
$$

No gauge fixing is needed for this inversion. The interior double-layer trace has
the constant vector as an (approximate, at finite resolution) eigenvector with
eigenvalue one, $\mathcal B\boldsymbol 1\approx\boldsymbol 1$, and the wake block
annihilates constants because $\mathcal C\boldsymbol 1=0$; hence
$\mathcal G\boldsymbol 1\approx\boldsymbol 1$. The constant mode is mapped to
itself, not to zero, so it is not in the null space of $\mathcal G$ and the
finite-wake systems are solvable as written.

If

$$
\mathcal Cq=0,
$$

then the shift is exactly

$$
\boldsymbol\mu_{V}-\boldsymbol\mu_{E}=q.
$$

The attached-wake strength is unchanged, and the additional source and body-doublet
layers satisfy

$$
\mathcal S_{\mathrm e}\sigma+\mathcal D_{\mathrm e}q=0.
$$

Thus the two finite-wake formulations have the same exterior field under this simple,
physically relevant sufficient condition.

It is not an unconditional necessary condition. Let

$$
\delta\boldsymbol\mu
=q-\mathcal G^{-1}\mathcal W\mathcal Cq.
$$

The actual necessary-and-sufficient condition for equality at the selected exterior
observations is

$$
\boxed{
\mathcal O\!\left(
\mathcal S_{\mathrm e}\sigma
+\mathcal D_{\mathrm e}\delta\boldsymbol\mu
+\mathcal W_{\mathrm e}\mathcal C\delta\boldsymbol\mu
\right)=0.
}
$$

Nonzero Kutta-trace terms could, in principle, be canceled by the rest of the complete
correction field. Requiring $\mathcal Cq=0$ removes that complication and is the
clean sufficient condition relevant to a continuous trace at a sharp trailing edge.

## Suggested alternative formulation for `semiinfinite_wake=false`

The finite-wake solve can retain the continuum exterior-field equivalence without
evaluating the free-wake scalar potential at every body panel. The statement is
conditional: at the discrete level, the equivalence is only as accurate as the panel
approximation to the Green identity and the quadrature described below.

Throughout this section $\mathcal G$ and $\mathcal W$ denote the finite
attached-wake operators written $\mathcal G_\Delta$ and $\mathcal W_\Delta$ in the
Task 2 construction; the subscript is dropped because no semi-infinite operator
appears here.

Let there be $N$ body panels and $M$ shedding edges. Then
$q,\boldsymbol\mu_E\in\mathbb R^N$, while the Kutta map
$\mathcal C\in\mathbb R^{M\times N}$ produces an $M$-component attached-wake
strength. Define the trace-shifted body unknown and the corresponding
trailing-edge trace difference by

$$
\widetilde{\boldsymbol\mu}=\boldsymbol\mu_E+q,
\qquad
c=\mathcal Cq.
$$

The purpose of the shift is to replace the omitted full wake-potential trace $q$ by
the source strengths obtained from wake velocity. To see every algebraic step, start
with the explicit-potential finite-wake system,

$$
(\mathcal B+\mathcal W\mathcal C)\boldsymbol\mu_E
=-\mathcal S\boldsymbol\sigma_0-q,
$$

and substitute $\boldsymbol\mu_E=\widetilde{\boldsymbol\mu}-q$:

$$
\mathcal B\widetilde{\boldsymbol\mu}-\mathcal Bq
+\mathcal W\mathcal C\widetilde{\boldsymbol\mu}
-\mathcal W\mathcal Cq
=-\mathcal S\boldsymbol\sigma_0-q.
$$

Since $c=\mathcal Cq$, collecting the terms that multiply the shifted unknown gives

$$
\mathcal B\widetilde{\boldsymbol\mu}
+\mathcal W\mathcal C\widetilde{\boldsymbol\mu}
=-\mathcal S\boldsymbol\sigma_0-q+\mathcal Bq+\mathcal Wc.
$$

The Green identity $\mathcal S\sigma=(I-\mathcal B)q$ implies

$$
-q+\mathcal Bq=-\mathcal S\sigma.
$$

Substitution therefore gives

$$
\boxed{
\mathcal B\widetilde{\boldsymbol\mu}
+\mathcal W\left(\mathcal C\widetilde{\boldsymbol\mu}-c\right)
=-\mathcal S(\boldsymbol\sigma_0+\boldsymbol\sigma).
}
$$

Equivalently, in terms of the existing finite attached-wake operator,

$$
\boxed{
\mathcal G\widetilde{\boldsymbol\mu}
=-\mathcal S(\boldsymbol\sigma_0+\boldsymbol\sigma)+\mathcal Wc.
}
$$

This second boxed form makes the relation to production exact and minimal: the
matrix and the source right-hand side are identical to the production
velocity-through-sources system, so
$\widetilde{\boldsymbol\mu}=\boldsymbol\mu_V+\mathcal G^{-1}\mathcal Wc$. The
entire solve-side change is one added right-hand-side vector $\mathcal Wc$ — the
potential of the $M$ attached panels carrying strengths $c$ evaluated at the $N$
body centroids, a single transient influence evaluation that requires no
persistent second wake system. The remaining change is downstream: every consumer
of attached-wake strength must apply the affine correction derived next.

The shifted system is also gauge invariant, which the explicit-potential system is
not. Its right-hand side depends on the free-wake potential only through
$\boldsymbol\sigma$ (a velocity, unchanged by adding a constant to $\psi$) and
through $c=\mathcal Cq$ (unchanged for paired rows because
$\mathcal C\boldsymbol 1=0$). The explicit-potential right-hand side contains
$-q$ itself and therefore shifts with the potential gauge. For paired shedding
this is a genuine robustness advantage of the velocity-only route.

The first boxed form exposes the physical attached-wake strength:

$$
\boxed{
\boldsymbol\gamma=\mathcal C\widetilde{\boldsymbol\mu}-c.
}
$$

This affine Kutta relation is essential. Indeed,

$$
\mathcal C\widetilde{\boldsymbol\mu}-c
=\mathcal C(\boldsymbol\mu_E+q)-\mathcal Cq
=\mathcal C\boldsymbol\mu_E.
$$

Thus the shift changes the body coefficient vector without changing the physical
attached-wake circulation. Using $\mathcal C\widetilde{\boldsymbol\mu}$ by itself
would introduce the spurious increment $c$ into both the attached wake and the next
shed free-wake row.

### Obtaining $c$ from velocity

Only the shedding-edge-sized vector $c\in\mathbb R^M$, rather than the full
$N$-component body trace $q$, is needed. For a shedding edge with paired lower and
upper panels, FLOWPanel's Kutta orientation gives

$$
c_e=q_{\mathrm{upper}}-q_{\mathrm{lower}}
=\int_{\boldsymbol x_{\mathrm{lower}}}^{\boldsymbol x_{\mathrm{upper}}}
\boldsymbol u_{\mathrm f}\mathbin{\boldsymbol\cdot}d\boldsymbol\ell.
$$

This follows from the fundamental theorem for line integrals. Within the body
interior, the free-wake field must have a single-valued potential satisfying
$\boldsymbol u_{\mathrm f}=\boldsymbol\nabla\psi$. Hence the integral from the
lower control point to the upper control point is
$\psi(\boldsymbol x_{\mathrm{upper}})-\psi(\boldsymbol x_{\mathrm{lower}})$.
The integrand must be the wake-only velocity, not the freestream, body-induced
velocity, or an already accumulated total velocity.

Choose a path contained in the closed body, with its open portion in the body
interior. The straight centroid-to-centroid chord is convenient only after a geometry
check confirms that the whole chord remains inside; this is not automatic for a
nonconvex body. For a valid straight chord, parameterize

$$
\boldsymbol x(s)=\boldsymbol x_{\mathrm{lower}}+s d,
\qquad
d=\boldsymbol x_{\mathrm{upper}}-\boldsymbol x_{\mathrm{lower}},
\qquad 0\le s\le1.
$$

Then

$$
c_e=\int_0^1\boldsymbol u_{\mathrm f}(\boldsymbol x(s))
\mathbin{\boldsymbol\cdot}d\,ds.
$$

The endpoint trapezoid approximation is

$$
\boxed{
c_e\approx\tfrac12\left(
\boldsymbol u_{\mathrm f,\mathrm{upper}}
+\boldsymbol u_{\mathrm f,\mathrm{lower}}
\right)\mathbin{\boldsymbol\cdot}
\left(\boldsymbol x_{\mathrm{upper}}-\boldsymbol x_{\mathrm{lower}}\right).
}
$$

For a smooth wake field along the chord, its absolute quadrature error is
$O(|d|^3)$. The smoothness premise, however, fails exactly where this integral is
needed: the trailing-edge chord lies at panel-scale distance from the attached-wake
edge and the first free-wake row, so $\boldsymbol u_{\mathrm f}$ varies on the
attached-panel and regularization scale and the error constant (the second
derivative of the integrand) is large. Because $\psi$ is harmonic throughout the
body interior, path independence permits a better-conditioned choice: deform the
path forward into the body, for example lower centroid $\to$ an interior point one
panel length or so upstream of the trailing edge $\to$ upper centroid. This
replaces a near-singular quadrature by a smooth one at the cost of one or two
extra batched velocity targets per edge. If refinement on the straight chord is
preferred instead, evaluate the midpoint
$\boldsymbol x_{\mathrm m}=(\boldsymbol x_{\mathrm{lower}}+
\boldsymbol x_{\mathrm{upper}})/2$ and use Simpson's rule,

$$
\boxed{
c_e\approx\tfrac16\left(
\boldsymbol u_{\mathrm f,\mathrm{lower}}
+4\boldsymbol u_{\mathrm f,\mathrm m}
+\boldsymbol u_{\mathrm f,\mathrm{upper}}
\right)\mathbin{\boldsymbol\cdot}d.
}
$$

All edge midpoints should be evaluated in one batched wake-influence call. A curved
or piecewise interior path instead requires quadrature of its own parametrization.

The endpoint rule can reuse existing control-point values only if the wake-only
velocity contribution is retained or captured during the current wake influence
pass. That reuse avoids additional influence evaluations, but it can require
bookkeeping when the production path immediately accumulates wake velocity with
other contributions. Midpoint or higher-order refinement still requires the
corresponding additional, preferably batched, wake evaluation.

This construction has several analytic limitations. A multiply connected body
interior requires circulation compatibility around every noncontractible loop before
a single-valued interior potential, and hence exact path independence, exists.
Regularized vorticity whose support overlaps the interior likewise invalidates exact
path independence, even if the nominal filament centerline remains outside. For the
particle-wake production path this is not hypothetical: measured near-wake particle
overlap in this project is $\sigma/h\approx4.2$, so regularized cores near the
trailing edge very likely intersect the body interior where the integration path
lives. For particle wakes the two-interior-path comparison below should therefore be
treated as a mandatory diagnostic, not an optional one, and the deformed interior
path above directly reduces the exposure. Discrete
trace/flux representation error and near-singular quadrature error are not removed by
this reformulation.

For a paired edge, compare trapezoid and Simpson values to estimate quadrature error.
Also compare two wholly interior paths between the same control points to detect
vorticity leakage or loss of path independence. An exterior comparison is not an
appropriate leakage test because an exterior path can legitimately link the wake.

An unpaired shedding row does not provide a potential difference: its $c_e$ depends
on an absolute potential reference that velocity data alone cannot supply. Setting
$c_e=0$ for such a row is therefore an explicit approximation or gauge convention,
not a consequence of the velocity integral. Exact velocity-only equivalence is
recommended for paired shedding edges; unpaired edges require separately supplied
gauge information.

If the discrete Kutta trace is neglected or explicitly approximated as zero, the
formulation reduces to the currently documented velocity-through-sources solve,

$$
\mathcal G\widetilde{\boldsymbol\mu}
=-\mathcal S(\boldsymbol\sigma_0+\boldsymbol\sigma).
$$

The shifted vector $\widetilde{\boldsymbol\mu}$ is not the explicit-potential
coefficient vector $\boldsymbol\mu_E$. Nevertheless, when it is combined with the
added source layer and the affine attached-wake strength above, it produces the same
exterior body field: the difference is the exterior-invisible combination
$\mathcal S_{\mathrm e}\sigma+\mathcal D_{\mathrm e}q=0$. Recovering
$\boldsymbol\mu_E$ itself would require reconstructing $q$ by integrating the
surface velocity and choosing a potential gauge; exact recovery of that
gauge-dependent vector is unnecessary for exterior velocity, pressure, and
circulation.

The exterior cancellation can be checked directly. Relative to the
explicit-potential representation, the alternative adds the body layers
$\mathcal S_{\mathrm e}\sigma+\mathcal D_{\mathrm e}q$. Its attached-wake change
is zero because

$$
\left(\mathcal C\widetilde{\boldsymbol\mu}-c\right)
-\mathcal C\boldsymbol\mu_E=0.
$$

The Green construction makes the added body layers exterior-invisible, so the total
exterior potential and velocity are unchanged. Equality of pressure follows only
when pressure is reconstructed from this complete, consistently corrected field; a
postprocessor must not interpret $\widetilde{\boldsymbol\mu}$ as
$\boldsymbol\mu_E$ while omitting the affine wake correction.

### Discrete and implementation qualifications

At constant-panel, centroid-collocated resolution, the continuum identity is usually
only approximated:

$$
\mathcal S_h\boldsymbol\sigma_h
\approx(I-\mathcal B_h)q_h.
$$

Consequently, computing $c$ accurately does not by itself make the alternative and
explicit-potential finite-panel solves exactly identical. Trace/flux representation,
self-limit consistency, wake regularization, and near-singular quadrature remain
possible error sources.

Implementation also requires more than adding $\mathcal Wc$ to the right-hand side.
Any operation that obtains an attached-wake or newly shed strength from the stored
body coefficients must use
$\mathcal C\widetilde{\boldsymbol\mu}-c$. Exterior influence evaluation,
free-wake shedding, circulation reporting, and potential/pressure reconstruction
must all follow the same affine convention. Otherwise the solve equation and the
subsequently evaluated physical system represent different wakes.

Computationally, this approach can reuse the free-wake velocities evaluated at body
control points to form $\sigma$ when their wake-only contributions are retained. The
endpoint rule then avoids additional influence evaluations; refined quadrature adds
only trailing-edge-scale, preferably batched evaluations for $c$. Both variants avoid
free-wake scalar-potential evaluation at all body control points and an auxiliary
dense boundary solve.

### Implementation approaches and FMM compatibility

All implementations must distinguish the linear wake coupling used during the solve
from the affine strength of the physical attached wake. The common workflow is

$$
\mathcal G\widetilde{\boldsymbol\mu}
=-\mathcal S(\boldsymbol\sigma_0+\boldsymbol\sigma)+\mathcal Wc,
\qquad
\boldsymbol\gamma=\mathcal C\widetilde{\boldsymbol\mu}-c.
$$

Matrix assembly and Krylov products remain linear: for a trial vector $x$, an FMM
operator application uses the attached-wake strength $\mathcal Cx$. The constant
$-c$ belongs on the right-hand side, not in the operator. After convergence,
physical field evaluation, first-row shedding, circulation output, and pressure
reconstruction must instead use
$\boldsymbol\gamma=\mathcal C\widetilde{\boldsymbol\mu}-c$. Five implementation
approaches are viable.

1. **Recover $\boldsymbol\mu_E$ from the full wake potential.** Evaluate the body
   trace $q$ in one batched scalar-potential FMM call, form $c=\mathcal Cq$, solve
   the shifted system, and then set
   $\boldsymbol\mu_E=\widetilde{\boldsymbol\mu}-q$. Restore the body sources to
   $\boldsymbol\sigma_0$ and use the existing implicit attached wake
   $\mathcal C\boldsymbol\mu_E$. The resulting stored body and wake strengths have
   the same interpretation as in the explicit-potential formulation, so existing
   physical influence and shedding paths can remain unchanged. This option is
   unavailable when the wake backend exposes only vector potential or velocity, and
   it defeats the objective of a velocity-only formulation.

2. **Add an explicit correction-wake system.** Leave the existing body and FMM
   representation at $\mathcal C\widetilde{\boldsymbol\mu}$, and add a second set of
   attached panels with strengths $-c$. Their known potential contributes
   $-\mathcal Wc$ at the body; moving that contribution to the right-hand side gives
   the required $+\mathcal Wc$. The correction system must participate in every
   physical potential and velocity evaluation so that the combined attached strength
   is $\mathcal C\widetilde{\boldsymbol\mu}-c$. Its panels must use geometry,
   orientation, triangulation, kernel offset, and regularization identical to the
   implicitly coupled attached panels. Any mismatch makes the cancellation
   discretely inconsistent.

3. **Store an affine correction in `RigidWakeBody`.** Continue to store body-panel
   coefficients as $\widetilde{\boldsymbol\mu}$, but give attached-wake evaluation an
   independent effective strength
   $\mathcal C\widetilde{\boldsymbol\mu}-c$. Both direct and FMM paths then need a
   strength channel or buffer for the attached panels rather than always deriving
   their strengths from the body coefficient currently being evaluated. They also
   need an explicit linear/operator mode for matrix assembly and Krylov products; in
   that mode the buffer contains $\mathcal Cx$ and suppresses the constant $-c$.
   Physical mode restores the affine correction after the solve.

4. **Obtain exact $c$ from a trailing-edge-only scalar-potential evaluation.**
   $c$ needs $q$ only at the $2M$ paired trailing-edge control points, not the
   full $N$-point body trace of approach 1. Evaluate wake scalar potential in one
   batched call at just those targets, form $c=\mathcal Cq$ there, and keep every
   other body evaluation velocity-only. This eliminates the velocity-line-integral
   quadrature error entirely at $O(M)$ cost. `PanelWake` already exposes the
   scalar-potential output channel through its derivatives switch in
   `src/FLOWPanel_wake.jl` (panel-only wakes, `include_final_filament=false`), so
   this hybrid is implementable today for panel wakes; particle wakes provide no
   scalar-potential kernel and must fall back to the interior line integral.

5. **Create a first-class attached-wake source system.** Decouple attached panels
   from the body-panel storage and store $\boldsymbol\gamma$ explicitly. During
   matrix assembly or a Krylov product, update that storage to $\mathcal Cx$. After
   convergence, update it to
   $\mathcal C\widetilde{\boldsymbol\mu}-c$. Direct and FMM backends must accept the
   body and attached wake together as source systems, so all target evaluations see
   the same physical strengths. This design makes the attached wake an independent
   source of truth rather than a special side effect of body evaluation.

The preferred choice depends on the objective:

- **Easiest velocity-only implementation:** use a separate correction-wake system,
  because the existing body and FMM operator code remains unchanged.
- **Most accurate overall when scalar potential exists:** evaluate $q$ with one
  batched FMM call and recover $\boldsymbol\mu_E$, because this avoids the discrete
  Green-identity and velocity-line-integral errors.
- **Best accuracy-per-cost for panel wakes:** the trailing-edge-only potential
  evaluation of approach 4, because it obtains exact $c$ at $O(M)$ cost while the
  body solve stays velocity-through-sources.
- **Most robust general velocity-only implementation:** use first-class,
  independently stored attached-wake strengths, because the corrected
  $\boldsymbol\gamma$ becomes the single source of truth for field evaluation,
  shedding, output, and pressure.
- **Most novel or publishable implementation:** use the native affine-Kutta/FMM
  formulation with velocity-derived $c$, independently stored attached strengths,
  and demonstrated equivalence among direct, FMM, and explicit-potential results.

Validation must cover both the linear solve and the complete physical
representation. On identical frozen body and wake geometry, require direct/FMM
agreement for potential and velocity at body control points, wake nodes, and exterior
targets. Compare the corrected shifted formulation with the explicit-$q$ baseline.
Changing $c$ while holding the geometry fixed must change only the right-hand side
and the post-solve physical wake strengths; it must not change any column of
$\mathcal G$. Verify that the first shed row uses
$\mathcal C\widetilde{\boldsymbol\mu}-c$, that $c=0$ exactly recovers current
behavior, and that paired-edge results are invariant to an additive potential gauge.
Finally, any correction or independent attached-wake system must be bundled
automatically with its body. Requiring callers to remember a second system invites
silent omission from pressure reconstruction, probes, output, or other physical
influence evaluations.

### Obtaining $c$ with a vortex-particle wake

A particle wake exposes no scalar potential, which rules out approach 1 (full-trace
recovery of $\boldsymbol\mu_E$) and approach 4 (trailing-edge-only potential
evaluation). The remaining routes to $c$ are alternatives for the c-consuming
approaches 2, 3, and 5, and their viability depends strongly on the configuration.

**Green-identity inversion (preferred for rotors and broken-up wakes).** The trace
can be recovered from data the production path already computes: solve the body-only
system

$$
(I-\mathcal B)q=\mathcal S\boldsymbol\sigma,
\qquad
c=\mathcal Cq,
$$

with $\boldsymbol\sigma$ formed from the sampled wake velocities. This is
velocity-only by construction and indifferent to wake age, origin, or structure. It
introduces one gauge degree of freedom — $(I-\mathcal B)$ has the constant vector in
its approximate null space — but that is exactly the harmless constant mode: fix it
numerically (area-weighted zero-mean trace, or a constrained least-squares solve)
and $c=\mathcal Cq$ is invariant to the choice for paired rows. Unpaired rows feel
the gauge here as in every other route. Because the discrete compatibility condition
$\int\sigma\,dS=0$ holds only approximately, solve the near-singular system as a
constrained least-squares problem rather than a raw inversion. Two further
properties are notable. First, if $q$ is *defined* as this solution, then
$-q+\mathcal Bq=-\mathcal S\boldsymbol\sigma$ holds exactly at the discrete level,
so the shifted and explicit-potential formulations agree exactly with each other;
the discrete Green-identity error becomes a shared discretization error rather than
an inconsistency between routes. Second, if regularized wake vorticity overlaps the
body interior — so that no exact interior potential exists and every route is
approximating — this solve projects the velocity data onto the closest field a
harmonic interior potential can represent, instead of returning whichever
path-dependent line integral was computed. Cost: although
$\mathcal B=\mathcal G-\mathcal W\mathcal C$ is a rank-$M$ update of the factored
$\mathcal G$, the solve operator here is $I-\mathcal B=I-\mathcal G+\mathcal W\mathcal C$,
and $I-\mathcal G$ is a full-rank perturbation — Sherman--Morrison--Woodbury against
the $\mathcal G$ factorization does **not** apply. The implementation instead builds
and caches a one-time bordered LU of $(I-\mathcal B)$ augmented with the
area-weighted gauge row (same $O(N^3)$ order as constructing the `Backslash`
factorization, amortized over all steps for rigid geometry). Without a dense
factorization (Krylov/FMM), it is a second, body-only Krylov solve; see the lagging
remark below for amortization.

**Redundant interior paths and local harmonic fits (cheap fallback and diagnostic).**
Integrate $\boldsymbol u_{\mathrm f}$ along two or three distinct interior paths per
edge and combine by least squares, or fit a local curl-free (harmonic polynomial)
field to a small stencil of sampled velocities around the trailing edge and take the
potential difference from the fit. With overlapping particle cores no single path is
correct; the spread across paths is a direct error bar on $c$ and doubles as the
mandatory leakage diagnostic. Cost is nearly free: $M$ edges $\times$ $P$ paths
$\times$ $n$ nodes of extra velocity targets (hundreds of points against the
thousands of body centroids already evaluated), appended to the existing batched
wake-influence call, plus a tiny per-edge least-squares solve.

**Lagged or filtered $c$ (cost and noise management, not accuracy).** Because $c$
enters only the right-hand side and the affine output correction — never the matrix
— it may be recomputed every $k$ steps, warm-started from the previous step, or
low-pass filtered without affecting the conditioning of the solve. Lagging
amortizes cost on the wake-convection timescale; filtering suppresses the
high-frequency jitter of particles convecting past the integration path or stencil.
Neither corrects a systematic bias in the underlying estimator: filtering reduces
variance and faithfully preserves bias.

**Wing-only options: hybrid near-panel wake and shadow sheets.** Retaining the first
free-wake rows as doublet panels (potential available via approach 4's machinery)
before particle conversion, or evaluating the solid-angle potential of an equivalent
doublet sheet reconstructed from shed-particle connectivity, both rest on the
assumption that the wake nearest the trailing edge is the freshest and most coherent
wake. A rotor violates this: returning wake layers pass close beneath the blade,
dominate the trailing-edge potential difference, and are exactly the old,
particle-converted, possibly rolled-up or turbulent portion of the wake — so the
fresh-row panels capture the wrong subset, and the equivalent-sheet geometry is
ill-defined once connectivity tangles (its strengths remain valid by Kelvin, its
spanning surface does not). These two options are appropriate for wing-like
configurations with short, coherent wake horizons, and should not be used for
rotors or broken-up wakes.

## Three different equivalence questions

### 1. Continuum Green representation

The identity $\mathcal S\sigma=(I-\mathcal B)q$ is exact only under these analytic
preconditions:

- $\partial\Omega$ is closed and consistently oriented with the outward normal.
- Each relevant body interior is simply connected, and the free-wake vorticity and
  singular support are disjoint from its closure. More generally, a multiply
  connected interior must satisfy the circulation compatibility conditions that make
  the interior velocity an exact gradient.
- The resulting single-valued interior wake potential is harmonic and has a
  well-defined trace and normal derivative. No global single-valued exterior
  potential is required.
- Kernel signs, the chosen one-sided limits, and the potential gauge are consistent
  between every term.

### 2. Finite-panel velocity/source versus explicit potential

Exact discrete equality additionally requires $q$ and $\sigma$ to be representable
in the selected doublet and source panel spaces, with the same geometry, quadrature,
self limits, signs, and gauge in both routes. Constant panels populated by centroid
samples generally provide only an approximation to a varying trace and normal flux.
The continuum identity does not make that collocation error disappear.

### 3. Finite attached-plus-free wake versus semi-infinite wake

Even exact equality of the two finite-wake boundary formulations does not prove
equality with the semi-infinite model. The attached/free composite must reproduce the
semi-infinite sheet in geometry, circulation distribution, orientation, continuity at
the handoff, and far-wake closure. A time-marched wake can reach a steady shape that is
rolled up or deflected relative to the prescribed semi-infinite sheet; mere steady
state is therefore insufficient.

## Requirements and error mechanisms

The earlier P1--P5 labels are useful as diagnostics, but they are not a single
necessary-and-sufficient list.

### Analytic assumptions

The closed-boundary, disjoint-support, harmonicity, sign, limit, and gauge assumptions
above govern whether the continuum identity applies at all.

### Discretization errors

- **P1 -- trace/flux representation.** Point-sampled constant source panels may not
  represent the free-wake normal flux, especially near the trailing edge where the
  field varies on the scale of the first attached-wake length.
- **P2 -- near-singular evaluation.** A free-wake edge close to a body centroid can
  amplify quadrature, core-size, and regularization differences. The body attached
  panel and free-wake row use distinct regularization parameters, so numerical
  agreement must be measured rather than inferred.

### Wake-representation errors

- **P3 -- handoff geometry and orientation.** The downstream edge of the attached
  panel and the upstream edge of the first free-wake row must coincide, have compatible
  strength, and cancel with opposite orientation in the composite sheet. Source-route
  and potential-route contributions do not cancel term by term on the boundary; the
  relevant comparison is the complete exterior field.
- **P4 -- shape and closure.** Finite length, the starting vortex, far-wake closure,
  wake roll-up, and any mismatch from the prescribed semi-infinite geometry can all
  change the exterior solution. A long wake addresses truncation only, not shape
  equality.

### Kutta-coupling contamination

- **P5 -- Kutta response to the trace shift.** At finite resolution,
  $\mathcal Cq$ can be nonzero because upper and lower trailing-edge centroid samples
  need not have identical free-wake potential. The shift in the solved doublets can
  then alter both the in-matrix attached wake and subsequently shed strength. This
  mechanism naturally suggests sensitivity to the attached-wake length, but the
  derivation only exposes the term; it does not establish that the term causes the
  observed deficit.

When every shedding row is a difference between paired upper and lower panels, an
additive constant in $q$ is harmless because the Kutta map annihilates the constant
mode,

$$
\mathcal C\boldsymbol 1=0.
$$

This statement is not valid for an unpaired shedding row. Such a row requires an
absolute potential reference; assigning its trace correction to zero is an explicit
approximation or gauge convention. The companion note,
`docs/dirichlet_potential_theory.md`, explains the corresponding
constant-potential gauge for the paired case.

## Numerical checks that distinguish the mechanisms

These are proposed diagnostics, not claims about an existing passing regression.

1. **Reproduce the production linear system.** Freeze a settled marched state,
   assemble the finite attached-wake operator and the source right-hand side from the
   sampled free-wake velocity, and check the solver residual. This verifies that the
   manual system matches production; it does not test Green equivalence.

2. **Test the discrete Green identity without Kutta coupling.** On a closed body with
   no attached-wake columns, evaluate $q$ and $\sigma$ from the same free wake and
   check
   $$
   \mathcal S\sigma\approx(I-\mathcal B)q
   $$
   after applying a common gauge. Failure diagnoses trace/flux representation,
   quadrature, self-limit, or sign inconsistency (P1/P2).

3. **Test the finite attached/free handoff off the body.** With identical prescribed
   strengths, compare the potential and velocity of the finite attached panel plus
   free-wake rows with those of the intended continuous composite at exterior probe
   points. A localized mismatch diagnoses handoff geometry, winding, strength, or
   regularization (P3).

4. **Measure the Kutta trace directly.** Evaluate the free-wake potential at paired
   upper and lower trailing-edge centroids, apply $\mathcal C$, and compare its
   spanwise distribution and attached-length dependence with the circulation deficit.
   Agreement would support P5; disagreement would reject it as the dominant mechanism.

5. **Separate wake shape from boundary formulation.** Evaluate the frozen marched
   composite wake and the prescribed semi-infinite wake at the same off-body probes,
   using matched circulation. Differences persisting after the handoff region diagnose
   wake-shape or far-closure mismatch (P4), not the Green substitution.

6. **Cross-substitute complete systems.** Solve the explicit-potential and
   velocity-through-sources finite systems on identical frozen geometry, then evaluate
   the full exterior correction in the boxed observation identity. This is the direct
   necessary-and-sufficient finite-formulation check; comparing raw doublet vectors
   alone is misleading because they differ by a trace-like shift.

`test/dirichlet_potential_test.jl` is an old ad hoc diagnostic outside the maintained
test matrix. It still passes the removed `cp_outer` keyword and currently fails before
reaching its Green-identity assertion. It must not be cited as an active verified
regression. The companion note, `docs/dirichlet_potential_theory.md`, describes an
updated diagnostic design without changing that script.
