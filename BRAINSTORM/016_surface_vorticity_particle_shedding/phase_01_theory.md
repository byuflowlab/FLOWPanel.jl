# Phase 1 — Surface-Vorticity Theory and Derivation

**Status:** APPROVED

**Prerequisite:** none

**Phase approval:** [x]

**Architecture authorization:** Yes. Phase 2 is unblocked.

## 1. Scope

This phase derives an opt-in replacement for the line-particle conversion in
`PanelParticleWake`. The current `_convert_to_particles!` interprets every
piecewise-constant wake panel as a vortex ring: jumps between neighboring
panel strengths become internal line particles, while uncancelled boundary
edges close the sheet. That representation is retained unchanged as the
legacy default.

The proposed method reconstructs the regular, surface-distributed part of the
same vorticity from neighboring panel-centroid strengths. It deposits that
field by area quadrature and separately retains only the physical
wake-sheet perimeter. It does not change the body solve, the panel-wake
strength definition, or particle evolution after deposition. It does require
an explicit outgoing ghost-panel storage row and a revised conversion
transaction. Section 5 selects cancellation of the artificial
panel/particle handoff edge and rejects retention of that edge.

## 2. Physical and stored sign conventions

Let $\mu$ denote physical doublet strength and let $\hat\mu$ denote the value
stored in `PanelWake.strength`. The wake construction reverses the
body-shedding-edge node order to make the wake-node sequence contiguous and
stores

$$
\boxed{\hat\mu=-\mu.}
\tag{1}
$$

The physical surface-vorticity convention is therefore

$$
\boxed{
\boldsymbol\kappa
=\boldsymbol n\times\nabla_s\mu
=-\boldsymbol n\times\nabla_s\hat\mu,
}
\tag{2}
$$

where

$$
\nabla_s\mu
=\left(\boldsymbol I-\boldsymbol n\boldsymbol n^\mathsf T\right)
\nabla\mu .
$$

This distinction is required by the live package orientation. A `PanelWake`
cell is sent to the element kernel in the order

$$
\boldsymbol v_1\longrightarrow\boldsymbol v_2
\longrightarrow\boldsymbol v_3\longrightarrow\boldsymbol v_4
\longrightarrow\boldsymbol v_1,
\tag{3}
$$

with $\boldsymbol v_1,\boldsymbol v_4$ upstream and
$\boldsymbol v_2,\boldsymbol v_3$ downstream. Thus the downstream spanwise
edge is oriented from
$\boldsymbol r_L=\boldsymbol v_2$ to
$\boldsymbol r_R=\boldsymbol v_3$ and carries the stored ring strength
$\hat\mu$. This is also the orientation assumed by the existing
surface-vorticity reconstruction: using stored samples requires the minus
sign in (2). The Phase 4 affine-sheet sign test remains the executable
authority; an orientation change must fail that test rather than be absorbed
into a convention note.

Units provide an independent check:

$$
[\mu]=[\hat\mu]=L^2T^{-1},\qquad
[\nabla_s\mu]=[\nabla_s\hat\mu]=LT^{-1},\qquad
[\boldsymbol\kappa]=LT^{-1}.
$$

Since $[\Delta A]=L^2$,

$$
[\boldsymbol\Gamma_p]
=[\boldsymbol\kappa\Delta A]=L^3T^{-1},
$$

which is FLOWVPM particle strength (vorticity integrated over volume). The
line-vortex circulation supplied to the legacy shedding helper instead has
units $L^2T^{-1}$; multiplying by a line increment produces the same
particle-strength units.

## 3. Distributional edge jumps versus a smooth reconstruction

For piecewise-constant physical $\mu_h$ (equivalently stored
$\hat\mu_h=-\mu_h$), the exact distributional surface gradient
is supported on panel edges. Representing each panel as a constant-strength
vortex ring therefore gives:

- a line contribution proportional to the strength jump on each shared
  internal edge; and
- an uncancelled line contribution on the physical perimeter.

The current conversion samples those line contributions through
`method_trailing` and `method_unsteady`. It is a consistent representation of
the discontinuous field $\mu_h$, but mesh-scale changes in neighboring
panel values are replicated as discrete internal filaments.

The smooth strategy instead assumes that panel-centroid values sample an
underlying continuous field. It reconstructs $\nabla_s\hat\mu$ from stored
values, deposits the regular field (2) across panel area, and does not also
deposit the internal
edge jumps used to estimate that field. This is a change of reconstruction,
not an algebraic rearrangement at finite resolution; the two strategies are
expected to converge toward the same smooth-sheet limit when their
assumptions hold.

## 4. Ghost storage and tangential-gradient reconstruction

### 4.1 Meaning of `nwakerows` and the ghost row

For user-requested `nwakerows=N`, \(N\) remains the number of active
aerodynamic panel rows. Every `PanelWake` is planned to allocate

$$
N+2\ \text{node rows},\qquad N+1\ \text{strength rows}.
\tag{3a}
$$

Node rows \(1,\ldots,N+1\) bound the \(N\) active panels. The additional node
row \(N+2\), together with the already existing shifted-out strength slot
\(N+1\), can describe one outgoing ghost panel. Active source enumeration,
propagation, VTK panel output, and wake-length metadata continue to expose
only panels \(1,\ldots,N\). A pure `PanelWake` does not use the extra node row
for influence; its shifted-out strength may still supply the existing final
filament policy.

The ghost is the panel being converted, not a previously converted passive
marker. At a conversion transition:

1. propagate the \(N\) active rows;
2. stage the shift of old active row \(N\) into ghost row \(N+1\);
3. construct the newly shed active row 1;
4. reconstruct and preflight the ghost conversion; and
5. atomically commit the shifted wake and all ghost particles.

The ghost is never an aerodynamic panel source. Before the transition it
contributed as active row \(N\); immediately after the transition particles
replace it. The ghost slot is overwritten by the next outgoing active panel.

This supplies a real stencil when `nwakerows == 1`: the newly shed active
panel is row 1 and the outgoing ghost is row 2. The two stored node rows in
the current implementation describe only one panel; the additional \(N+2\)
node row is what makes the ghost a second panel with its own centroid.

### 4.2 Structured panel-centroid differences

Index the ghost panels by spanwise column $j$. Let $G=(N+1,j)$ be a
target ghost and \(A=(N,j)\) its immediately upstream active neighbor. The
streamwise observation is always first-order and one-sided:

$$
\boldsymbol d_s=\boldsymbol x_G-\boldsymbol x_A,\qquad
\delta_s\hat\mu=\hat\mu_G-\hat\mu_A .
\tag{4}
$$

Additional active rows do not change the stencil order. Varying
`nwakerows` changes the panel-wake length and handoff location, so an
aerodynamic difference between `nwakerows=1,2,3,...` must not be attributed
to a higher-order derivative.

For the spanwise direction use actual ghost-panel centroids:

- centered interior:
  $$
  \boldsymbol d_\eta
  =\boldsymbol x_{G,j+1}-\boldsymbol x_{G,j-1},\qquad
  \delta_\eta\hat\mu=\hat\mu_{G,j+1}-\hat\mu_{G,j-1};
  $$
- one-sided at a root/tip boundary:
  $$
  \boldsymbol d_\eta
  =\boldsymbol x_{G,j_{\rm in}}-\boldsymbol x_{G,j},\qquad
  \delta_\eta\hat\mu=\hat\mu_{G,j_{\rm in}}-\hat\mu_{G,j}.
  $$

No `Das` transition, body strength, virtual trailing-edge panel, or
fictitious outside-zero panel participates. Those are not local
panel-centroid samples of the outgoing ghost, can belong to another physical
time, and are not present for every attachment formulation.

Project each displacement into the target tangent plane,

$$
\widetilde{\boldsymbol d}_q
=\left(\boldsymbol I-\boldsymbol n_G\boldsymbol n_G^\mathsf T\right)
\boldsymbol d_q .
\tag{5}
$$

Using an orthonormal tangent basis
$\boldsymbol T=[\boldsymbol t_1\ \boldsymbol t_2]$, assemble

$$
A_{q:}=(\boldsymbol T^\mathsf T\widetilde{\boldsymbol d}_q)^\mathsf T,
\qquad b_q=\delta_q\hat\mu ,
$$

and recover tangent coordinates $\boldsymbol g$ from

$$
\boldsymbol g=A^+ \boldsymbol b,\qquad
\nabla_s\hat\mu=\boldsymbol T\boldsymbol g,\qquad
\boldsymbol\kappa=-\boldsymbol n_G\times\boldsymbol T\boldsymbol g .
\tag{6}
$$

The pseudoinverse is evaluated with a documented scale-aware rank tolerance.
With two independent observations, (5) recovers both tangent components and
is exact for affine $\hat\mu$ (and hence affine $\mu$) on a planar sheet. On
nonuniform grids, using
physical centroid displacements rather than index spacing preserves that
property.

### 4.3 Rank-deficient stencils

Rank deficiency is a supported state:

- rank 2: both tangent components are observable;
- rank 1: return the deterministic minimum-norm gradient in the observable
  tangent direction and mark one component unobservable;
- rank 0: return zero reconstructed gradient and mark both components
  unobservable.

The result must include a diagnostic carrying at least panel identity, rank,
available directions, and whether the streamwise component is observable.
The minimum-norm completion is a numerical representation of missing
information; it must not be reported as evidence that the true missing
derivative is zero.

`nwakerows == 1` is not intrinsically rank deficient under the ghost design:
the newly shed active row supplies the streamwise neighbor. Rank loss can
still arise from a single spanwise panel, coincident centroids, folded
geometry, or nearly collinear streamwise/spanwise observations. Such cases
retain the deterministic observable-component policy; they are not repaired
with `Das`, a virtual TE panel, or an outside-zero value.

### 4.4 Warped panels

Each quadrilateral is treated as the bilinear surface

$$
\boldsymbol x(\xi,\eta)
=\sum_{a=1}^{4}N_a(\xi,\eta)\boldsymbol x_a,
\qquad (\xi,\eta)\in[0,1]^2.
$$

Centroid differences determine the locally observable derivatives. At a
quadrature point, bilinear tangents define the metric and normal; the
reconstructed derivatives are mapped into that local tangent plane before
applying (2). Thus particle strengths follow the physical warped surface
rather than a flat projected area.

This is a local, low-order reconstruction. It is exact for affine strength on
a planar panel grid and covariant under rigid rotation. On a warped grid it
is expected to converge under panel refinement, not reproduce an arbitrary
curved-sheet gradient exactly on a coarse mesh. Folded, inverted,
near-zero-area, or metric-singular bilinear panels are validation errors.

## 5. Signed perimeter and handoff ledger

### 5.1 Discrete-Stokes identity in package signs

For a smooth planar patch with positively oriented boundary tangent
$\boldsymbol t$,

$$
\int_S \boldsymbol\kappa\,dA
=\int_S\boldsymbol n\times\nabla_s\mu\,dA
=\oint_{\partial S}\mu\,\boldsymbol t\,ds
=-\oint_{\partial S}\hat\mu\,\boldsymbol t\,ds .
\tag{7}
$$

The last expression is the negative of the stored-strength ring assembly.
This is the expected distributional ledger: the regular vorticity of a
smooth bounded patch and its stored ring-oriented perimeter have equal and
opposite integrated vectors. For constant strength the regular term is
exactly zero, while the stored ring is a closed constant-strength line whose
integrated vector is also exactly zero.

Every edge of an outgoing conversion patch is classified before deposition:

1. an internal converted edge is replaced by reconstructed surface
   vorticity, so no line-jump particles are deposited there;
2. a true root/tip wake-sheet perimeter remains a package-oriented line
   closure; and
3. the active-panel/particle interface is an artificial partition and is
   cancelled, never treated as an outside-zero physical boundary.

### 5.2 One-panel handoff vector

Let $A$ be the final active panel immediately upstream of outgoing ghost
$G$, and let $\boldsymbol r_L,\boldsymbol r_R$ be their shared spanwise-edge
endpoints in package order. The downstream edge of $A$ contributes

$$
\hat\mu_A(\boldsymbol r_R-\boldsymbol r_L),
$$

whereas the upstream edge of $G$ is traversed in the opposite direction and
contributes

$$
-\hat\mu_G(\boldsymbol r_R-\boldsymbol r_L).
$$

Their assembled handoff vector is therefore

$$
\boxed{
\boldsymbol H
=(\hat\mu_A-\hat\mu_G)
(\boldsymbol r_R-\boldsymbol r_L).
}
\tag{8}
$$

This sign is not discretionary. `_final_filament_strength` returns
$-\hat\mu_A$ when `unsteady_filament=false`, on the same
$\boldsymbol r_L\to\boldsymbol r_R$ edge, and therefore cancels the active
panel's $+\hat\mu_A$ downstream edge. The stored-strength difference in the
existing `_convert_to_particles!` uses the same orientation for its legacy
unsteady line particles.

For a planar affine fixture, let $\boldsymbol e_s$ point from the active
centroid to the ghost centroid, separated by $\Delta s$, and let the shared
edge vector be
$\boldsymbol r_R-\boldsymbol r_L=\Delta\eta\,\boldsymbol e_\eta$, with
$\boldsymbol n=\boldsymbol e_s\times\boldsymbol e_\eta$. The streamwise
part of (2) is

$$
\boldsymbol\kappa_s
=-\boldsymbol n\times
\left(
\frac{\hat\mu_G-\hat\mu_A}{\Delta s}\boldsymbol e_s
\right)
=\frac{\hat\mu_A-\hat\mu_G}{\Delta s}\boldsymbol e_\eta .
\tag{9}
$$

Consequently,

$$
\boxed{
\int_G\boldsymbol\kappa_s\,dA
=(\hat\mu_A-\hat\mu_G)
(\boldsymbol r_R-\boldsymbol r_L)
=\boldsymbol H.
}
\tag{10}
$$

The same result follows in a rotated planar frame because it uses only the
oriented tangent basis. If strength is constant,
$\hat\mu_A=\hat\mu_G$, both $\boldsymbol H$ and
$\int_G\boldsymbol\kappa_s\,dA$ are exactly zero.

### 5.3 Multiple spanwise panels

For a row of ghost panels, each internal spanwise edge appears twice in the
stored ring assembly with opposite traversal. Equal strengths cancel;
unequal strengths leave the familiar jump line. Summing across columns
telescopes those internal terms. Under the smooth reconstruction the
spanwise part of
$-\boldsymbol n\times\nabla_s\hat\mu$ replaces those internal jumps by area
quadrature. It is therefore not also deposited as line particles.

Only the two non-wrapping outer edges survive the telescoping sum. They are
the true root/tip wake-sheet perimeter and remain line particles. For a
planar affine field, the signed row ledger is

$$
\boldsymbol\Gamma_{\rm area}
+\boldsymbol\Gamma_{\rm root/tip}
=\sum_j\boldsymbol H_j .
\tag{11}
$$

The spanwise-gradient part of
$\boldsymbol\Gamma_{\rm area}$ is balanced by the root/tip closure, while
its streamwise-gradient part is precisely the sum of (10). Thus constant,
affine-streamwise, affine-spanwise, and combined planar-affine fields all
satisfy the same ledger. A wrapping surface has no root/tip term.

### 5.4 Repeated conversions

Conversion is a transaction on one identified outgoing ghost row. At
transaction $k$:

1. the former final active row is staged as ghost $G_k$ exactly once;
2. the new final active panel's downstream edge is cancelled with the
   `unsteady_filament=false` sign;
3. the complete reconstructed vorticity of $G_k$ and only its true root/tip
   closure are deposited atomically; and
4. $G_k$ ceases to be a panel source and its storage slot may be reused only
   for $G_{k+1}$.

Induction on $k$ then gives one conversion per physical patch. Equation (10)
replaces the assembled active/ghost difference at each transition, internal
spanwise terms telescope by section 5.3, and neither side retains an
artificial row-boundary filament. The only line particles accumulated over
successive conversions lie on the true physical perimeter.

### 5.5 Selected policy — Alternative A, cancelled edge

Alternative A is selected. The smooth strategy cancels the active-panel
downstream edge and deposits the complete reconstructed surface vorticity
plus true root/tip closure. Equations (8)--(11) show that it preserves the
integrated constant and planar-affine circulation ledger.

The replacement is integrated-strength preserving, not pointwise
near-field identical. It replaces a mesh-scale line representation with
finite-core area particles. Agreement is expected away from the handoff, or
under consistent refinement of panel size, target spacing $h$, and particle
core size $\sigma$. Coarse pointwise velocity or velocity-gradient
differences near the handoff are a representation effect, not by themselves
a conservation failure.

For later diagnostics define the signed handoff residual

$$
\boxed{
\boldsymbol R_{H,j}
=\boldsymbol\Gamma_{{\rm area},s,j}-\boldsymbol H_j,
}
\tag{12}
$$

where $\boldsymbol\Gamma_{{\rm area},s,j}$ is the area-quadrature strength
from the reconstructed streamwise derivative. Also report the row closure
residual

$$
\boldsymbol R_{\rm row}
=\boldsymbol\Gamma_{\rm area}
+\boldsymbol\Gamma_{\rm root/tip}
-\sum_j\boldsymbol H_j .
\tag{13}
$$

Both residuals must vanish for constant and planar-affine fixtures and
converge under consistent refinement for warped or non-affine cases.

### 5.6 Rejected policy — Alternative B, retained filament

Alternative B is rejected. Retaining the unsteady handoff filament leaves
the line vector $\boldsymbol H$ in (8), which already represents the affine
streamwise contribution proved in (10). Adding the complete reconstructed
area field would count that contribution twice.

Subtracting the streamwise component from area deposition could repair the
integrated ledger only by violating the defining particle invariant

$$
\boldsymbol\Gamma_p
=\boldsymbol\kappa(\boldsymbol x_p)\Delta A_p
$$

and by retaining the mesh-scale line representation that the smooth strategy
is intended to replace. It would no longer be a complete
surface-vorticity-area conversion. Alternative B therefore does not advance
to Phase 2.

### 5.7 Root/tip near-field caveat

Newly deposited root/tip surface or closure particles can lie close to the
uncancelled root/tip vortex filament of the retained active panel wake. That
panel filament may induce substantially larger velocity and velocity-gradient
magnitudes at boundary particle locations than at panel-interior particle
locations. The concern is cross-representation panel-on-particle influence,
not alignment among the new particles themselves.

Candidate-particle diagnostics must compare root/tip and interior
panel-induced velocity-gradient norms before commit. Non-finite values,
instability, or nonconvergent growth under spanwise, \(\sigma\), overlap, and
line-spacing refinement trigger reconsideration of a conservative diffuse
root/tip closure.

## 6. Surface quadrature and particle strengths

Let the requested particle smoothing width be $\sigma>0$ and target overlap
be $\mathrm{overlap}>0$. Define

$$
h=\frac{\sigma}{\mathrm{overlap}} .
\tag{14}
$$

For a panel, estimate the two physical parametric lengths $L_\xi,L_\eta$
from its bilinear centerlines and choose

$$
n_\xi=\max\!\left(1,\left\lceil L_\xi/h\right\rceil\right),\qquad
n_\eta=\max\!\left(1,\left\lceil L_\eta/h\right\rceil\right).
\tag{15}
$$

This gives a near-uniform two-dimensional target-spacing grid, not a fixed
particle count. Place one candidate particle at every parametric subcell
centroid,

$$
\xi_a=\frac{a-\tfrac12}{n_\xi},\qquad
\eta_b=\frac{b-\tfrac12}{n_\eta},\qquad
\boldsymbol x_p=\boldsymbol x(\xi_a,\eta_b).
\tag{16}
$$

Use the physical bilinear-surface Jacobian

$$
J(\xi,\eta)
=\left\lVert
\frac{\partial\boldsymbol x}{\partial\xi}
\times
\frac{\partial\boldsymbol x}{\partial\eta}
\right\rVert
$$

to define positive subcell area weights

$$
\Delta A_p
=J(\xi_a,\eta_b)\frac{1}{n_\xi n_\eta}.
\tag{17}
$$

If Phase 2 chooses exact geometric subcell areas, those may replace (17), but
they must remain physical rather than projected areas and preserve the same
refinement behavior.

The deposited strength is

$$
\boxed{
\boldsymbol\Gamma_p
=\boldsymbol\kappa(\boldsymbol x_p)\Delta A_p .
}
\tag{18}
$$

Consequently, by construction,

$$
\sum_{p\in P}\boldsymbol\Gamma_p
=\sum_{p\in P}
\boldsymbol\kappa(\boldsymbol x_p)\Delta A_p ,
\tag{19}
$$

the selected quadrature of the panel's reconstructed surface vorticity.
When $n_\xi=n_\eta=1$, the sole candidate is exactly the bilinear panel
centroid. Exact-zero strengths are not inserted into FLOWVPM because its
relaxation path divides by particle-strength magnitude; hence constant
$\mu$ (equivalently constant $\hat\mu$) produces zero distributed particles
rather than stored zero-strength particles.

## 7. Conservation statements and limitations

The method guarantees only what it constructs:

- panelwise particle-strength sums equal the chosen surface quadrature;
- constant strength creates no distributed vorticity and preserves the true
  perimeter closure;
- planar affine fields are reconstructed exactly when the stencil has rank 2;
- the cancelled-edge ledger is exact for constant and planar-affine fields,
  including streamwise and spanwise variation; and
- repeated conversions neither reconvert a ghost patch nor leave an
  artificial row-boundary filament.

It does not by itself guarantee exact kinetic energy, helicity, impulse, or
all higher vorticity moments after replacing panels by regularized particles.
Finite particle core size, overlap, quadrature, rank-deficient gradients,
warped geometry, later particle merging, relaxation, and viscous/SFS models
can change those quantities. Panel-on-particle near-field gradients at the
root/tip and any rank loss must be diagnosed explicitly.

The smooth method is inappropriate for a genuine discontinuity in $\mu$ (or
$\hat\mu$):
smoothing a physical jump would replace a line vortex by a diffuse sheet.
Topology or user configuration must identify any such discontinuity as a
physical boundary and retain its line contribution.

## 8. Phase 2 requirements

After explicit Phase 1 approval, architecture must decide:

- the concretely typed strategy API and legacy default;
- explicit active-row capacity despite the always-allocated ghost geometry;
- the \(N+2\)-node/\(N+1\)-strength storage, source exclusion, ghost staging,
  and newly shed row construction;
- ownership and reuse of gradient and quadrature workspaces;
- the cancelled-edge handoff transaction selected in section 5.5;
- how $\sigma$ and overlap are supplied for two-dimensional surface
  deposition;
- transactional particle-capacity preflight;
- diagnostics, validation, failure semantics, metadata, replay, and warm
  starts;
- how `method_trailing` and `method_unsteady` remain untouched for legacy
  conversion and whether the smooth strategy reuses either policy only for
  retained perimeter particles.

No architecture choice is approved by this theory document.

## 9. Verification requirements

Phase 4 must include:

1. storage/source mapping for `nwakerows=1,2,3`, including a ghost that never
   contributes panel influence;
2. `nwakerows == 1`: exact one-sided affine reconstruction from the newly
   shed active row and outgoing ghost;
3. constant $\mu$/$\hat\mu$: no distributed particles, zero handoff vector,
   and circulation only on the true physical perimeter;
4. affine $\mu$/$\hat\mu$ in both streamwise and spanwise directions on
   planar and rigidly rotated sheets: exact gradient, package sign, and
   integrated strength;
5. nonuniform planar and warped quad grids: convergence under panel,
   spanwise, and particle-spacing refinement;
6. the selected cancelled-edge policy: one panel, multiple spanwise panels,
   constant strength, affine streamwise strength, affine spanwise strength,
   combined affine strength, and repeated-conversion discrete-Stokes ledgers;
7. resolution floor: a nonzero rank-2 case with one particle at the bilinear
   centroid when both subdivision counts are one;
8. no dependence on `Das`, body strength, virtual TE samples, or an
   outside-zero neighbor;
9. root/tip versus interior panel-induced velocity-gradient diagnostics;
10. exact default regression for legacy locations, strengths, counts,
   metadata, replay, and warm starts;
11. legacy-edge-jump versus selected-smooth-strategy A/B validation with
   identical non-strategy settings: suddenly started wing first, then rotor
   hover, comparing induced velocity, integrated circulation/impulse,
   particle count, stability, and $C_L/C_T$, with panel-row, timestep,
   $\sigma$, and overlap refinement.

## 9b. Amendment (2026-08-01) — divergence-form discretization

Approved by Ryan on 2026-08-01, after Phase 3 Stage 3 implementation exposed
that the continuous statements of §5 admit two inequivalent discretizations.

**Alternative A (cancelled edge) still stands.** What changes is how
$\boldsymbol\kappa$ is *discretized* for deposition.

### 9b.1 The two discretizations

Equations (7)–(13) are continuum identities and are unaffected. But eqs. (9)–(11)
silently assume the *uniform-grid* limit: (10) requires the centroid separation
$\Delta s$ to equal the ghost panel's own streamwise extent $h_G$, which holds
only when neighbouring rows have equal extent. Wake rows stretch as they convect,
so in general $\Delta s\neq h_G$ and the two candidate rules diverge:

- **Reconstruction:** $\partial_s\hat\mu=(\hat\mu_A-\hat\mu_G)/\Delta s$,
  from the centroid stencil of Phase 2 §5.2.
- **Divergence form:** $\partial_s\hat\mu=(\hat\mu_A-\hat\mu_G)/h_G$, i.e.
  the face jump smeared over the panel's own area.

Equivalently: the divergence form deposits, on each panel, exactly the filament
content of the faces assigned to it,

$$
\boldsymbol\kappa_j=\frac{\boldsymbol V_j}{A_j},\qquad
\boldsymbol V_j=\boldsymbol H_j
+\tfrac12\!\!\sum_{\text{spanwise faces of }j}\!\!
(\Delta\hat\mu)\,(\text{edge vector}),
$$

with the upstream face taken whole (panel $A$ survives the transaction, so its
share cannot be deposited inside it), the downstream face taken not at all
(already cancelled by the retained filament under Alternative A), and each
spanwise face split evenly between the two panels that convert together — a
split that moves *where* the vorticity sits, never how much.

### 9b.2 Why divergence form is selected

The conversion is then, exactly, "redistribute the stored ring assembly's
filaments: smear the internal partitions into area vorticity, keep the physical
boundaries as lines." Consequences:

1. **Circulation is conserved exactly on any mesh** — stretched, warped,
   non-affine — not only on uniform affine fixtures. Eq. (11) becomes an
   identity rather than a limit.
2. **No neighbour geometry is ever required**, only $\hat\mu_A$. This removes
   the `nwakerows == 1` difficulty entirely: there the upstream row's nodes do
   not exist at conversion time (`shed_wake!` writes only its strength;
   `update_TE!` writes its nodes on the *next* step), but its strength is the
   body's shed value. `nwakerows == 1` is the production rotor configuration.
3. **§5.3 is vindicated as written.** "Only the two non-wrapping outer edges
   survive the telescoping sum" identifies them with the *stored cell values*
   $+\hat\mu_1$ and $-\hat\mu_{n_{\rm cols}}$ — which is correct in this
   scheme, and is the legacy root/tip filament reused verbatim. Under the
   reconstruction rule those same lines would instead have to carry the field
   extrapolated to the sheet edge, or the row ledger fails by one cell width.

### 9b.3 The trade, stated honestly

This is not a free win, and the cost should not be hidden.

Verified numerically (smooth analytic field, smoothly graded planar mesh,
compared against exact $-\boldsymbol n\times\nabla_s\hat\mu$):

| $M{=}K$ | $h$ | err divergence | err reconstruction | difference | conserv. err div. | conserv. err recon. |
| --- | --- | --- | --- | --- | --- | --- |
| 20 | 4.8e-2 | 2.3e-2 | 1.1e-2 | 3.4e-2 | 1.3e-16 | 3.5e-2 |
| 80 | 1.2e-2 | 5.0e-3 | 3.2e-3 | 8.2e-3 | 8.0e-18 | 8.2e-3 |
| 320 | 3.1e-3 | 1.2e-3 | 8.3e-4 | 2.0e-3 | 1.3e-16 | 2.0e-3 |

Both converge to the same value at first order, and their difference vanishes at
first order; on a *uniform* mesh they are identical to $2\times10^{-15}$. The
difference scales as $(r-1)/2$ with $r$ the row-to-row extent ratio (1 % stretch
→ 0.5 %, 5 % → 2.4 %, 10 % → 4.5 %).

At a *fixed* resolution the reconstruction is the more accurate pointwise
estimator — its $\boldsymbol\kappa$ error falls as $r$ grows while the
divergence form's rises. What the divergence form buys is a conservation error
that is **exactly zero at every resolution**, where the reconstruction leaks
circulation at $O(r-1)$ on *every* conversion. A per-step circulation leak biases
thrust systematically and cumulatively; a local $\boldsymbol\kappa$ error does
not accumulate. That, and not a blanket claim of superiority, is the reason for
the selection.

The reconstruction of Phase 2 §5.2 is retained and still evaluated on every
panel, as a diagnostic: the recorded
$\lVert\boldsymbol\kappa_{\rm cons}-\boldsymbol\kappa_{\rm recon}\rVert$
is a direct measure of grid non-uniformity, and is pinned to converge under
refinement.

### 9b.4 Startup edge and streamwise attribution (2026-08-01, second pass)

Ryan's review of the divergence form raised two consequences that the first pass
of §9b missed.

**The aft face is a physical boundary until the first handoff.** §9b.1 assigns
the outgoing panel's downstream face *nothing*, on the grounds that the retained
filament already cancels it. That is true only once `particle_handoff_active` is
set. On the **first** conversion the filament is still in its legacy form, so the
aft face carries an uncancelled net $\hat\mu_G-\hat\mu_D$ — the starting vortex —
and discarding the row without depositing it destroys that circulation outright.
Before the first handoff the aft face is the sheet's true trailing boundary, not
an interface with particles, and is deposited **whole**, exactly like the root
and tip closures of §5.3.

**The upstream/downstream split is a free parameter.** Let $\alpha$ be the
fraction of the upstream face a conversion deposits and $\beta=1-\alpha$ the
downstream fraction. Every inter-row face is then deposited exactly once, in two
instalments by the two conversions that flank it — the older row takes $\alpha$
as its upstream face, the younger takes $\beta$ as its downstream face one step
later. All values of $\alpha$ conserve circulation exactly. The un-deposited
remainder must stay on the panel side, which fixes the retained filament to

$$
\Gamma_{\rm filament}
=-\bigl(\alpha\,\hat\mu_{\rm last}+(1-\alpha)\,\hat\mu_{\rm converted}\bigr),
$$

an expression that reproduces both pre-existing forms at its endpoints:
$\alpha=1$ is the cancelled edge of Alternative A, $\alpha=0$ is the legacy
unsteady filament.

$\alpha=\tfrac12$ makes the streamwise difference *centered*, so the leading
grid-stretching error cancels. Measured against an exact analytic gradient on a
smoothly graded mesh:

| $M$ | $h$ | $\alpha=1$ | $\alpha=0$ | $\alpha=\tfrac12$ |
| --- | --- | --- | --- | --- |
| 20 | 4.8e-2 | 3.1e-2 | 3.0e-2 | 7.9e-4 |
| 80 | 1.2e-2 | 7.3e-3 (1.04) | 7.2e-3 (1.02) | 5.4e-5 (1.96) |
| 320 | 3.1e-3 | 1.8e-3 (1.01) | 1.8e-3 (1.00) | 3.4e-6 (1.99) |

**Second order versus first**, and 40x more accurate at coarse resolution.

This is **not** the Alternative B rejected in §5.6. That rejection was for
retaining the *full* filament while also depositing the *complete* area field, a
genuine double count. A weighted split is consistent; §5 simply never considered
it.

The cost is fidelity at the handoff. $\alpha=1$ leaves *nothing* at the
panel/particle partition, which is §5.5's entire argument: the partition is
artificial and should carry no line. $\alpha=\tfrac12$ leaves a half-jump
filament there — it halves the mesh-scale line artifact rather than removing it,
and does so exactly where freshly deposited particles sit, which is §5.7's
near-field concern. $\alpha=1$ therefore remains the selected default, with
$\alpha=\tfrac12$ and $\alpha=0$ available and the choice deferred to measured
near-field evidence in Phase 3 Stage 4b.

## 10. Decision and progress log

Phase 1 is technically complete and recommends Alternative A, the
cancelled-edge handoff, for architecture. Alternative B is rejected for
double-counting the affine streamwise component or, after subtraction,
violating the complete-area particle invariant. This technical completion
does not approve Phase 1: the approval checkbox remains unchecked and Phase 2
remains blocked pending explicit user approval.

- **2026-07-29 — Phase 1 APPROVED (clear-context review + explicit user
  approval).** A clear-context reviewer independently re-derived equations
  (7)–(11) by hand (including the spanwise-affine root/tip balance in the row
  ledger) and verified all seven code-anchored claims at exact lines: the
  implicit sign flip `strengthi - strengthj` (`src/FLOWPanel_simulate.jl:1205-1211`),
  vertex ordering (`src/FLOWPanel_wake.jl:289-292`), `_final_filament_strength`
  cancellation (`src/FLOWPanel_wake.jl:1337-1340`), legacy handoff `Γ-Γ_tm1`
  matching $\boldsymbol H$ (`src/FLOWPanel_wake.jl:1276-1315`), current
  N+1/N+1 storage with no ghost node row (`src/FLOWPanel_wake.jl:148-152`),
  the zero-Γ guard (`src/FLOWPanel_wake.jl:673-679`), and body
  κ = n×∇μ consistency (`src/FLOWPanel_simulate.jl:152-181`). The Alternative B
  rejection is provably forced by (10). Non-blocking polish noted: cite the
  sign-chain code anchors in §2; the zero-Γ guard is an uncommitted
  working-tree change and must be committed before Phase 3 lands; duplicate
  equation label (3)/(3a). User approved the review; Phase 2 is authorized.

- **2026-08-01 — Divergence-form amendment (§9b), approved by Ryan.**
  Alternative A stands; its discretization becomes divergence form. Eqs. (9)–(11)
  are recognized as the uniform-grid limit of the ledger, the root/tip closure is
  confirmed as the stored cell value (legacy filament reused verbatim), and the
  `nwakerows == 1` stencil difficulty dissolves because no neighbour geometry is
  needed. Convergence of the two rules to the same limit was verified
  numerically before adoption; the accuracy-vs-conservation trade is recorded in
  §9b.3 rather than glossed.

- **2026-08-01 — §9b.4 added (startup edge + attribution), approved by Ryan.**
  The aft face is a physical boundary until the first handoff and must be
  deposited whole (omitting it deleted the starting vortex). The
  upstream/downstream split is a free parameter `alpha`; all values conserve
  exactly, `alpha=1/2` is second-order where `alpha=1` and `alpha=0` are first,
  and the retained filament is fixed by `alpha`. Default stays `alpha=1`
  (Alternative A, nothing at the artificial partition); Stage 4b decides on
  near-field evidence.

- **2026-07-29 — Handoff ledger completed; Phase 1 technically complete.**
  Distinguished physical $\mu$ from stored $\hat\mu=-\mu$, audited the live
  ring and final-filament orientations, derived
  $\boldsymbol H=(\hat\mu_A-\hat\mu_G)
  (\boldsymbol r_R-\boldsymbol r_L)$, selected cancelled-edge Alternative A,
  rejected retained-filament Alternative B, and defined handoff residuals.
  Awaiting explicit user approval; Phase 2 remains blocked.

- **2026-07-29 — Theory technically completed.** Derivation and acceptance
  scenarios recorded. Awaiting explicit user approval; Phase 2 is blocked.
- **2026-07-29 — Theory reopened by user.** Adopted an always-allocated
  outgoing ghost panel so `nwakerows == 1` has a real one-sided stencil.
  Corrected the root/tip near-field concern to panel-filament-induced
  gradients at new particle locations. Reopened the handoff ledger to derive
  cancelled-edge and retained-filament alternatives. Phase 1 is in progress
  and Phase 2 remains blocked.
