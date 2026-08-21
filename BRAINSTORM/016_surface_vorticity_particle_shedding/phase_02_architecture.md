# Phase 2 — Surface-Vorticity Conversion Architecture

**Status:** APPROVED (2026-08-01), as amended in §13

**Prerequisite:** approved Phase 1 theory

**Phase approval:** [x]

**Implementation authorization:** Yes. The user explicitly approved this
architecture on 2026-08-01. Phase 3 is unblocked; Phase 4 remains blocked on
explicit Phase 3 approval. Read §13 before implementing — it carries the
selected storage decision and two corrections to §10.

## 1. Approved-theory handoff

Phase 1 selected the cancelled-edge handoff. The smooth conversion deposits
the complete reconstructed surface vorticity and true root/tip streamwise
closure, deposits no downstream handoff filament, and makes the retained
panel wake's final filament cancel its current downstream edge after the
first conversion. The retained-filament alternative was rejected because it
double-counts the affine streamwise component or requires an incomplete area
field.

This corrects the original Phase 2 scaffold: Phase 1 did select a handoff
ledger before it was approved. The architecture below makes that selected
ledger concrete; it does not reopen the choice.

Item 015's current `live_rows` source exclusion is the integration baseline.
The ghost-row exclusion introduced here is a separate invariant and must not
replace, weaken, or infer the `live_rows` exclusion.

## 2. Public strategy API

### 2.1 Conversion axis

Add one concrete strategy axis:

```julia
abstract type AbstractPanelParticleConversion end

struct LegacyEdgeJumpConversion <: AbstractPanelParticleConversion end

struct SurfaceVorticityConversion{TF} <: AbstractPanelParticleConversion
    sigma::TF
    overlap::TF
    rank_rtol::TF
    geometry_rtol::TF
    diagnose_nearfield::Bool
end

SurfaceVorticityConversion(
    sigma;
    overlap=1.3,
    rank_rtol=sqrt(eps(float(typeof(sigma)))),
    geometry_rtol=sqrt(eps(float(typeof(sigma)))),
    diagnose_nearfield=false,
)
```

The user selects the strategy through exactly one new keyword:

```julia
PanelParticleWake(body; conversion=LegacyEdgeJumpConversion(), kwargs...)
PanelParticleWake(body; conversion=SurfaceVorticityConversion(sigma), kwargs...)
```

`LegacyEdgeJumpConversion()` is the exact default. Omitting `conversion` and
supplying `conversion=LegacyEdgeJumpConversion()` must enter the same existing
conversion dispatch and produce the same locations, strengths, particle
counts, metadata interpretation, and allocations per legacy conversion step.
There is no automatic or warning-based fallback from the smooth strategy to
the legacy strategy.

The smooth constructor promotes `sigma`, `overlap`, `rank_rtol`, and
`geometry_rtol` to one floating type `TF`. All four must be finite;
`sigma > 0`, `overlap > 0`, `rank_rtol >= 0`, and `geometry_rtol > 0`.
The strategy owns fixed `sigma` and `overlap`, and all of its surface and line
sampling uses

$$
h = \frac{\sigma}{\mathrm{overlap}}.
$$

Every smooth particle deposited by one wake instance therefore has the same
fixed smoothing width `sigma`.

### 2.2 Legacy line-policy compatibility

`PanelParticleWake` changes the implementation defaults of `method_trailing`
and `method_unsteady` to private constructor sentinels. The public legacy
defaults remain exactly `OverlapPPS(1.3, 2)`:

```julia
PanelParticleWake(
    body;
    conversion=LegacyEdgeJumpConversion(),
    method_trailing=_default_wake_shedding_method,
    method_unsteady=_default_wake_shedding_method,
    kwargs...,
)
```

For `LegacyEdgeJumpConversion`, each sentinel resolves to a newly constructed
`OverlapPPS(1.3, 2)` exactly as it does today; an explicitly supplied
`WakeSheddingMethod` is stored and used unchanged.

For `SurfaceVorticityConversion`, both keywords must remain at their
sentinels. Explicitly supplying either legacy line policy is a configuration
error, even if its value equals the legacy default. This prevents silently
ignored options. The smooth strategy samples its only line contributions,
the true open-surface root and tip streamwise closures, with
`SigmaOverlap(sigma, overlap)`. It does not use `method_trailing` or
`method_unsteady`.

**Amended 2026-08-03 after Phase 3 review.** Smooth conversion also requires
`unsteady_filament=true` and `include_final_filament=true`. The former makes
the pre-handoff final filament carry the physical startup edge assumed by the
first-conversion ledger; the latter keeps the cancelling/complement filament
in `get_sources` after handoff. Supplying `false` for either option with
`SurfaceVorticityConversion` is a constructor-time `ArgumentError`. Both
options remain unchanged for `LegacyEdgeJumpConversion`.

### 2.3 Supported domain

The Phase 3 implementation must support:

- scalar-strength `VortexRing` and scalar-strength `ConstantDoublet` panel
  wakes;
- open shedding chains and geometrically wrapping shedding chains;
- every explicit active capacity `nwakerows = N` with `N >= 1`; and
- floating-point body/wake scalar types.

The smooth strategy rejects vector- or multi-component panel strengths,
unsupported panel kernels, integral or otherwise non-floating body types,
and ambiguous topology. A geometric wrap is accepted only when the first and
last node columns coincide within the scale-aware `geometry_rtol`; the
closing duplicate column is not a physical root or tip.

A shedding chain must be smooth except at its true open root and tip. A
physical spanwise discontinuity, blade break, or separate sheet component
must be represented as a separate shedding surface. The reconstruction never
turns an internal strength jump into a line vortex and never guesses that an
internal edge is a physical boundary.

## 3. `PanelWake` storage and active views

### 3.1 Explicit capacity

Add an explicit immutable active-row capacity `capacity::Int` to `PanelWake`.
For `nwakerows=N`, allocate per shedding surface:

- `nodes`: `3 × (N+2) × (nspan+1)`;
- `velocity`: the same shape as `nodes`; and
- `strength`: `kernel_dim × (N+1) × nspan`.

`nwakes[]` remains the number of active panel rows and satisfies
`0 <= nwakes[] <= capacity`. It is never inferred from an array dimension.
The ordinary active panel rows are `1:nwakes[]`, and ordinary active probe
node rows are `1:(nwakes[]+1)`.

When the wake is full, row `N+1` of `strength` and node rows `N+1:N+2`
can describe one non-source ghost panel. The overlap at node row `N+1` is
intentional: it is the downstream boundary of active panel row `N` and the
upstream boundary of ghost panel row `N+1`.

### 3.2 Visibility rules

The ghost row is scratch/continuation state, not an aerodynamic panel source.
It is excluded from:

- direct and FMM panel source counts and index maps;
- panel-wake propagation and velocity application;
- `ProbeWrapper` probe counts and index maps;
- `FilamentWrapper` except for the separately defined active final filament;
- VTK panel cells and VTK probe nodes;
- wake-length and `nwakerows` reporting; and
- overflow and active-row bookkeeping.

Only `nwakes[]+1` active probe-node rows are visible even though storage has
`N+2` rows. Plain `PanelWake` instances keep the unused scratch row zeroed.

Item 015 remains orthogonal. First restrict all views to the explicit active
range, then apply `live_rows[]` to exclude its reserved newest active rows
from old-wake sources. A ghost is never made visible by `live_rows == 0`, and
an active live row never becomes a ghost merely because it is source-excluded.

### 3.3 Startup and final-filament state

Add `particle_handoff_active::Array{Bool,0}` to `PanelWake`, initialized to
`false`. While false, startup behavior preserves the existing physical
downstream/startup edge. After the first successful smooth conversion
transaction it becomes `true`. In that state, the final panel filament uses
the negative strength of the current final active row, so it cancels that
row's downstream edge; it must not use the ghost strength or represent a
shifted-out row.

Legacy wakes never activate this flag and preserve the existing
`unsteady_filament` behavior exactly. For a smooth wake, the active handoff
rule owns final-filament selection once the flag becomes true. The flag is
committed atomically with particles and row state, and is persisted per step.

## 4. Concrete ownership

Extend `PanelParticleWake` with concrete type parameters and fields for:

```julia
conversion::TC
conversion_workspace::TW
conversion_diagnostics::TD
conversion_count::Array{Int,0}
```

`conversion_count[]` starts at zero and increments exactly once after each
successful smooth full-wake conversion. It does not increment on startup
steps, preflight failures, or legacy conversions.

Legacy instances store:

```julia
conversion_workspace === nothing
conversion_diagnostics === nothing
```

and dispatch directly to the present `_convert_to_particles!` behavior.
They must not acquire new per-step staging or diagnostic allocations.

Smooth instances store a concretely typed reusable
`SurfaceVorticityWorkspace{TF}`. It owns the maximum-size centroid, normal,
tangent-basis, neighbor, SVD, gradient, subdivision, quadrature, new-row,
shifted-row, ghost, and pending-particle buffers for all configured shedding
surfaces. Buffers grow only in an explicit capacity-management path before a
transaction; conversion does not expose a partly resized workspace.

`SurfaceVorticityDiagnostics{TF}` owns the latest committed conversion record
and a cumulative ledger. A failed preflight may return or attach a typed
failure diagnostic to the thrown exception, but it must not replace the
latest committed record or advance the cumulative ledger.

`PanelWake` owns active capacity, ghost-valid bookkeeping, and
`particle_handoff_active` because those fields govern panel source/filament
views. `PanelParticleWake` owns the conversion strategy, its workspace and
diagnostics, and conversion count because those govern deposition.

## 5. Virtual post-shift reconstruction

### 5.1 Full-wake preflight

Smooth conversion runs only when `nwakes[] == capacity == N`. Before mutating
the panel wake, build a virtual post-shift state in staging buffers:

1. Construct the new row-1 nodes and strengths from the body using the same
   geometry and stored-strength conventions as ordinary `shed_wake!`.
2. Virtually shift current active rows `1:N-1` to rows `2:N`.
3. Treat the current final active panel row `N` as outgoing ghost `G` at
   panel row `N+1`, including its two node rows and scalar strength.
4. Reconstruct and convert `G`.

The ghost's streamwise upstream neighbor in that virtual state is current
active row `N-1` when `N > 1`. When `N == 1`, it is the newly constructed
row 1. Thus the one-row case has a genuine two-panel one-sided stencil and
uses no body strength, `Das`, virtual trailing-edge panel, outside-zero value,
or previously converted marker.

If the stored smooth continuation state says a conversion should be possible
but the required current final row, upstream continuation row, or step state
is missing, conversion throws a continuation-state error. It must not
approximate the stencil or fall back to legacy conversion.

### 5.2 Deterministic tangent reconstruction

For every ghost panel, use its centroid and stored scalar strength
`muhat = -mu`. Define a deterministic orthonormal tangent basis from the
ghost-centroid normal: choose the Cartesian axis least aligned with the
normal (ties in `x`, `y`, `z` order), project and normalize it as `t1`, and
set `t2 = n × t1`. This makes rigid-rotation sign fixtures reproducible
without allowing a basis sign ambiguity to change the physical vector.

Form two geometric difference equations:

- one streamwise equation from the ghost and its virtual upstream neighbor;
- one spanwise equation from periodic centered neighbors for a wrapping
  chain, centered interior neighbors for an open chain, and one-sided
  root/tip neighbors for an open chain.

Project the two centroid-displacement vectors into `(t1,t2)` and solve the
resulting `2 × 2` system for the tangent components of
`grad_s(muhat)`. Scale the geometry matrix by the largest stencil length
before its SVD. Rank uses singular values of that scaled matrix and the
relative threshold `rank_rtol`; `geometry_rtol` separately rejects a
vanishing metric scale.

The SVD policy is:

- rank 2: use the full reconstructed gradient;
- rank 1: use the SVD pseudoinverse for the minimum-norm observable
  component and mark the orthogonal tangent component unobservable;
- rank 0: return a zero gradient and mark both tangent components
  unobservable.

Rank deficiency is supported and diagnostic, not a geometry failure. Singular
values, effective rank, condition (infinite below rank 2), and observable
right-singular directions are recorded per panel. SVD sign choices may affect
only the reported basis coordinates, never the reconstructed physical vector.

### 5.3 Bilinear panel geometry

Use the package's contiguous wake-panel vertex order exactly:

$$
x(0,0)=v_1,\qquad x(1,0)=v_2,\qquad
x(1,1)=v_3,\qquad x(0,1)=v_4.
$$

At each quadrature/particle point, compute the bilinear covariant derivatives
and oriented normal. Project the centroid-reconstructed physical gradient
into that local tangent plane, then use the approved stored-strength sign:

$$
\boldsymbol\kappa(x_p)
= -\boldsymbol n(x_p)\times
   \nabla_s\widehat{\mu}(x_p).
$$

Let `L_xi_max` and `L_eta_max` be the maximum physical lengths of the two
opposite bilinear-panel edges in each parametric direction. Choose

$$
n_\xi=\max\left(1,\left\lceil\frac{L_{\xi,\max}}{h}\right\rceil\right),
\qquad
n_\eta=\max\left(1,\left\lceil\frac{L_{\eta,\max}}{h}\right\rceil\right).
$$

Place one area particle at each bilinearly mapped parametric subcell centroid.
Compute that subcell's physical area with a `2 × 2` Gauss rule applied to the
norm of the bilinear Jacobian over the subcell. Deposit

$$
\boldsymbol\Gamma_p
= \boldsymbol\kappa(x_p)\,\Delta A_p
$$

with smoothing width `sigma`. FLOWVPM's vestigial scalar circulation field is
set to the nonnegative equivalent magnitude

$$
\Gamma_{\mathrm{scalar},p}
=\frac{\lVert\boldsymbol\Gamma_p\rVert}{\sqrt{\Delta A_p}}.
$$

Root/tip line particles retain signed scalar circulation under the existing
line-particle convention; they do not use the area-particle magnitude rule.

Candidates whose vector strength is exactly zero under the package's
zero-strength guard are elided before capacity accounting. Elision is counted
in diagnostics. No tolerance may silently erase a finite nonzero vector
unless a separately approved policy is added in a later phase.

### 5.4 Geometry validation

Preflight rejects non-finite nodes, strengths, centroids, derivatives,
normals, areas, or particle records. It also rejects, using
`geometry_rtol` relative to local panel/stencil scale:

- repeated or metric-singular vertices/edges;
- zero or near-zero Jacobian anywhere sampled by validation and quadrature;
- a Jacobian orientation inconsistent with the panel's reference orientation;
- folded, self-intersecting, or inverted bilinear panels; and
- inconsistent geometric wrapping.

Warped but consistently oriented bilinear quads are supported. Rank
deficiency of the strength-reconstruction stencil remains nonfatal unless the
geometry itself is metric-singular.

## 6. Perimeter and handoff ledger

Each pending candidate has one concrete classification:

```julia
@enum SurfaceParticleClass begin
    InteriorSurfaceParticle
    RootTipSurfaceParticle
    PerimeterLineParticle
end
```

`InteriorSurfaceParticle` and `RootTipSurfaceParticle` distinguish area
particles on interior versus root/tip-adjacent ghost panels.
`PerimeterLineParticle` is reserved for true open-chain root/tip streamwise
closure.

For an open chain, retain exactly two physical streamwise boundaries: root
and tip. Sample each from the ghost panel's upstream to downstream endpoint
with `SigmaOverlap(sigma, overlap)` and the signed circulation required by
the selected Phase 1 perimeter ledger. For a wrapping chain there is no
root/tip closure.

Do not deposit:

- internal spanwise panel-to-panel jumps;
- a line at a geometric wrapping seam;
- an outside-zero completion line;
- any downstream spanwise handoff filament; or
- any separate line representing the surface area integral.

Before the first full conversion, the physical startup edge remains active.
On the first successful transaction, `particle_handoff_active` becomes true.
Thereafter the retained final panel filament cancels the current active
downstream edge while the converted ghost's complete area vorticity and true
root/tip closure reside in particles. The handoff vector and discrete-Stokes
residual are computed before commit for every conversion.

## 7. Atomic conversion transaction

The smooth path is a preflight/commit transaction:

1. Validate the strategy, kernel, topology, current state, and virtual
   post-shift geometry.
2. Reconstruct every panel and compute every subdivision and quadrature area.
3. Materialize every nonzero area and root/tip line particle in pending
   buffers, including all FLOWVPM fields needed by insertion.
4. Compute panel and conversion diagnostics, signed perimeter ownership,
   handoff vector, and absolute/relative Stokes residual.
5. Compare the exact requested particle count with remaining
   `ParticleField` capacity.
6. Only if all checks pass, append all particles, commit shifted active rows,
   the new row, ghost scratch/bookkeeping, handoff state, diagnostics, and
   `conversion_count += 1` exactly once.

No live particle or panel state is mutated during steps 1–5. Particle
capacity includes only candidates remaining after exact zero-strength
elision. Insufficient capacity throws a typed
`PanelParticleCapacityError(requested, available)` carrying the preflight
counts. Both the particle field and all panel state remain bitwise unchanged
by that failure.

Other preflight failures use typed configuration, geometry, topology,
unsupported-kernel, or continuation-state errors with the shedding-surface
and panel/stencil identity where applicable. There is no partial commit, retry
inside the same call, warning-only continuation, or legacy fallback.

Phase 3 must arrange commit ordering so an unexpected insertion failure is
also atomic. It may reserve/validate the full FLOWVPM insertion range first
or provide rollback from copied affected columns, but a partly appended row
must never be observable.

## 8. Diagnostics contract

Use concrete records rather than untyped dictionaries for live diagnostics.
Names may be private, but the following fields and meanings are mandatory.

### 8.1 Per-panel record

Each ghost panel record contains:

- shedding-surface and spanwise panel indices;
- stencil rank;
- both scaled singular values and rank threshold;
- finite rank-2 condition number or `Inf`;
- two physical observable tangent directions and their observable flags;
- `n_xi` and `n_eta`;
- Gauss-integrated physical area;
- summed deposited area-particle vector strength;
- root/tip-surface ownership classification;
- count of requested and zero-strength-elided area particles; and
- geometry scale and minimum validated Jacobian metric.

### 8.2 Per-conversion record

Each committed conversion record contains:

- conversion ordinal and physical step identity;
- per-panel records;
- counts by `SurfaceParticleClass`;
- total physical area and total integrated area-particle vector strength;
- signed root and tip line strengths and ownership;
- panel-final-filament ownership before and after commit;
- the expected and deposited handoff vectors;
- absolute Stokes residual norm;
- relative Stokes residual using a scale-aware denominator that remains
  defined for a zero ledger;
- total zero-strength elisions;
- exact requested and available particle capacity; and
- handoff state before and after commit.

The cumulative diagnostic ledger sums count, area, integrated strength,
perimeter, handoff, and residual information without changing a previous
record.

### 8.3 Optional near-field record

When `diagnose_nearfield=true`, directly evaluate panel-induced velocity
gradient norms at every pending candidate before particle insertion and
record separate summary distributions for:

- interior surface candidates;
- root/tip surface candidates; and
- perimeter line candidates.

This diagnostic is panel-induced only. It must not mix in particle-particle
alignment effects, and root/tip values must not be pooled with interior
values. The optional work is absent, not zero-filled, when the flag is false.

## 9. Metadata, replay, and warm starts

### 9.1 Versioned conversion metadata

Add a versioned conversion block to `PanelParticleWake` metadata:

```text
conversion:
  version: 1
  type: LegacyEdgeJumpConversion | SurfaceVorticityConversion
  sigma: ...             # smooth only
  overlap: ...           # smooth only
  rank_rtol: ...         # smooth only
  geometry_rtol: ...     # smooth only
  diagnose_nearfield: ...# smooth only
```

Missing conversion metadata means `LegacyEdgeJumpConversion()` for backward
compatibility. An unknown version/type or a missing smooth parameter is an
error, not a legacy default.

Record `nwakerows` from `PanelWake.capacity`, never from node or strength
array dimensions. Legacy metadata keeps its current line-policy records.
Smooth metadata must not fabricate `method_trailing` or `method_unsteady`
settings.

### 9.2 Per-step continuation state

Persist, for every output step needed by replay or warm start:

- active `nwakes`, active rows, and existing Item 015 `live_rows` state;
- `particle_handoff_active`;
- `conversion_count`;
- the smooth strategy identity/parameters; and
- sufficient step identity/bookkeeping to prove whether the outgoing row was
  already deposited.

Smooth replay/warm-start with missing required per-step handoff or conversion
state throws a typed continuation-state error. It must not infer zero, repeat
a conversion, or risk duplicate deposition. Legacy files may continue to use
the historical defaults because missing conversion metadata selects legacy.

### 9.3 VTK and reconstruction

VTK continues to contain active panels and active probe nodes only; it does
not serialize the ghost as an aerodynamic cell. Loading allocates the new
`N+2` node/`N+1` strength storage, restores the active rows, and leaves ghost
scratch empty/invalid.

When a later full-wake smooth conversion needs a ghost, reconstruct it by
staging the then-current outgoing final active row during the next virtual
post-shift transaction. Do not require a passive ghost panel in VTK and do
not treat zeroed scratch as a strength sample.

Old files with `nwakerows=N` allocate the new ghost-capable arrays while
preserving exactly `N` as active capacity and preserving the historical wake
length.

## 10. Phase 3 implementation map

Phase 3 is limited to the following seams:

- `src/FLOWPanel_wake.jl`: strategy and sentinel types, constructor
  validation, `PanelWake` capacity/storage/views, concrete workspaces and
  diagnostics, reconstruction/quadrature/perimeter helpers, final-filament
  handoff selection, and conversion dispatch/transaction;
- `src/FLOWPanel_simulate.jl`: construct the staged new row, supply physical
  step identity, and replace conversion-before-shift with the atomic virtual
  post-shift smooth transaction while preserving the exact legacy call/shift
  ordering;
- `src/FLOWPanel_metadata.jl` and `src/FLOWPanel_replay.jl`: versioned
  conversion metadata, explicit active capacity, per-step state, backward
  legacy interpretation, and ghost-capable reconstruction; and
- `src/FLOWPanel_warmstart.jl`: exact restoration and validation of handoff,
  conversion count, active rows, Item 015 live-row state, and smooth
  continuation identity.

Implementation must preserve the uncommitted/current Item 015 integration
baseline. Phase 3 may not alter the mathematical strategy, add a fallback, or
begin aerodynamic validation without a new approved architecture amendment.

## 11. Focused verification map for Phase 3

Place focused mechanical verification in:

- `test/runtests_unit_wake.jl`;
- `test/runtests_unit_simulate.jl`;
- `test/runtests_unit_replay.jl`; and
- `test/runtests_unit_warmstart.jl`.

The implementation gate requires at least:

1. storage, active source/probe counts, ghost exclusion, final filament, and
   independent `live_rows` exclusion for `N=1,2,3`;
2. constant and affine planar strength fixtures, including streamwise,
   spanwise, combined gradients, exact package sign, and rigidly rotated
   covariance;
3. open and geometrically wrapping topology, with root/tip closure only for
   open chains and physical discontinuities represented as separate
   surfaces;
4. nonuniform and warped-grid refinement, physical-area convergence,
   subdivision floors, and rejection of folded/inverted/metric-singular
   panels;
5. rank-2, rank-1 minimum-norm, and rank-0 observable-direction diagnostics;
6. repeated-conversion discrete-Stokes ledgers and exact root/tip/handoff
   ownership;
7. startup physical-edge preservation followed by the first-conversion
   transition to a cancelled active downstream edge;
8. exact zero-strength elision and vestigial scalar-circulation values;
9. transactional configuration, geometry, continuation, and particle-capacity
   failures that leave both particle and panel state unchanged;
10. optional root/tip-versus-interior panel-induced velocity-gradient
    diagnostics;
11. omitted conversion versus explicit `LegacyEdgeJumpConversion()` exact
    regression, including legacy line policies and no new per-step
    allocations; and
12. new/old metadata, replay, and warm-start round trips, including failure
    on missing required smooth continuation state.

No simulations or Phase 4 aerodynamic validation are authorized by this
verification map.

## 14. Amendment (2026-08-01) — deposition is divergence form

Approved by Ryan. See Phase 1 §9b for the theory and the numerical verification.

**§5.2's rank-aware centroid reconstruction is no longer the deposition rule.**
It is retained, and still evaluated on every panel, purely as a diagnostic.
Deposition is

```
kappa_j = V_j / A_j,   Gamma_p = kappa_j * dA_p
V_j = H_j                                     # upstream face, whole
    + 0.5 * (muhat_j   - muhat_{j-1}) * (v2 - v1)   # left  spanwise face
    + 0.5 * (muhat_{j+1} - muhat_j)   * (v3 - v4)   # right spanwise face
```

with `A_j` the *sum of the subcell areas* (not an independent quadrature), which
is what makes `sum(kappa_j * dA_p) == V_j` exactly. `kappa_j` is **not** projected
onto the local normal — projection would break that sum. Outer faces of an open
chain contribute nothing to the area; they remain line particles at the legacy
strengths `+muhat_1` / `-muhat_ncols`, sampled with `SigmaOverlap(sigma, overlap)`
exactly as before.

**Amended 2026-08-01 (see Phase 1 §9b.4).** The upstream/downstream split is a
selectable parameter `attribution` on `SurfaceVorticityConversion`, with
`alpha = 1` (`:upstream`, default), `0` (`:downstream`), or `1/2` (`:split`):

```
V_j = alpha*H_j + beta_eff*D_j + 0.5*(spanwise face jumps)
D_j      = (muhat_G - muhat_D) * (v3 - v2)      # muhat_D = strength[1,N+1,j]
beta_eff = handoff_active_before ? 1 - alpha : 1
```

`beta_eff = 1` before the first handoff because the aft face is then the sheet's
physical trailing boundary (the starting vortex), with no earlier conversion to
have taken the `alpha` share; omitting it deletes that circulation. `muhat_D`
needs no new storage — `strength` already carries `nwakerows+1` rows and legacy
reads the same value as `Gamma_tm1`. `PanelWake` gains
`particle_handoff_weight` (`alpha`), and `_final_filament_strength` becomes
`-(alpha*strength[i_row] + (1-alpha)*strength[i_row+1])` when the handoff is
active, which reproduces both pre-existing forms at its endpoints.

Consequences for the rest of this document:

- **§2.3, §5.1, §5.2** — the smooth strategy needs only the upstream *strength*,
  never the upstream panel's geometry. The `nwakerows == 1` staging discussion is
  moot: nothing is staged, and `system` is threaded in solely to read
  `_get_wakestrength_mu`. This also means the §13.1 Option A/Option B choice is
  now immaterial to correctness — Option B stands, and its fallback trigger was
  never reached.
- **§6** — unchanged in intent; the "internal spanwise jumps are not also
  deposited as line particles" rule is now enforced structurally, since each face
  is either smeared or kept, exactly once.
- **§8.1** — the per-panel record replaces `streamwise_strength` and
  `handoff_residual` (identically zero by construction now) with
  `kappa_conservative`, `kappa_reconstruction`, and `kappa_difference`.
- **§8.2** — the per-conversion record adds `expected_total` (the exact filament
  content to transfer) and `deposited_total`; `residual_abs`/`residual_rel` now
  measure `deposited_total - expected_total`, which is round-off on any geometry.
  It also records `attribution` and `startup_edge_deposited`, and §8.1 adds
  `downstream_face` alongside `handoff`.
- **§3.3** — "before the first full conversion the physical startup edge remains
  active" is now enforced: that edge is deposited whole by the first transaction.
- **Verification** — note that a residual computed against the transaction's own
  `expected_total` is self-referential and cannot detect an omitted face. The
  binding check is external: every panel ring's perimeter vector sums to zero, so
  the wake's field-relevant content reduces to the retained filament, and the
  particles' gain must equal that quantity's loss across a `shed_wake!` call.

## 12. Acceptance gate and progress log

The exact API, ownership, reconstruction, cancelled-edge ledger, transaction,
diagnostics, validation/failure behavior, persistence contract, and Phase 3
file/test map are technically complete. **Approved 2026-08-01**, subject to the
amendments in §13.

## 13. Amendment (2026-08-01) — ghost storage decision and §10 corrections

### 13.1 Ghost storage: two options, one selected

Sections 3.1 and 3.2 were written in *post-shift* indices, in which the outgoing
row sits at panel row `N+1` and therefore appears to need storage beyond the
`N` active rows. Verification against the live source shows that framing is an
indexing convention rather than a storage requirement, because conversion runs
**before** the row shift (`src/FLOWPanel_simulate.jl:1229`). Both options are
recorded here; the choice was made at the start of Phase 3, as the user
directed.

**Option A — persistent ghost row (as originally drafted in §3.1/§3.2).**
Grow `PanelWake.nodes` to `3 × (N+2) × (nspan+1)`, add an explicit immutable
`capacity::Int`, and exclude the resulting ghost row from every consumer:
direct and FMM source counts and index maps, panel-wake propagation and
velocity application, `ProbeWrapper` counts and index maps, `FilamentWrapper`,
VTK panel cells and probe nodes, wake-length/`nwakerows` reporting, and
overflow bookkeeping. Roughly twenty call sites in `src/FLOWPanel_wake.jl`
(:205, :207-231, :233-281, :371-373, :283-327, :339-356, :521-544, :550-628,
:1323-1416) plus `src/FLOWPanel_simulate.jl` and the persistence layer.

**Option B — ghost read in place, pre-shift (SELECTED).** The ghost *is* the
current final active row. At conversion time the wake is full, so outgoing row
`N` already owns node rows `N` and `N+1` and strength row `N`, and its
streamwise upstream neighbour is existing active row `N-1`. Nothing new is
stored, nothing is hidden, and no consumer changes. The single quantity that
does not yet exist is the not-yet-shed new row 1, which the ghost needs as its
upstream neighbour **only when `nwakerows == 1`**; Phase 3 stages that row in
the conversion workspace, built from the body trailing edge and
`_get_wakestrength_mu` exactly as `shed_wake!` does at
`src/FLOWPanel_simulate.jl:1207-1210`.

**Selection: Option B.** It preserves every mathematical commitment of this
architecture — §5.1's virtual post-shift stencil, the one-sided two-panel
stencil at `N == 1`, the prohibition on `Das`/body-strength/virtual-TE/
outside-zero samples, the cancelled-edge ledger, and §9.3's rule that the
ghost has no cross-step lifetime — while removing the highest-risk part of the
diff. Option A's ghost is write-only scratch that §9.3 already forbids
persisting, so paying for it in every index map buys nothing.

Consequently §3.1 (explicit `capacity`, `N+2` allocation) and §3.2 (ghost
visibility exclusions) do **not** apply under Option B. `capacity` continues to
be derived as `size(wake.nodes[1], 2) - 1`, which is what the existing code and
the persistence layer already do. §3.3 (`particle_handoff_active` and the
final-filament rule) is unaffected and still applies. Item 015's `live_rows`
exclusion likewise remains untouched and orthogonal.

**Fallback trigger.** If Phase 3 finds a consumer that genuinely requires the
ghost to survive a step boundary, stop and re-open this choice rather than
partially adopting Option A.

### 13.2 Correction to §10 — wake serializers live in the replay module

§10 assigns conversion metadata to `src/FLOWPanel_metadata.jl`. That file owns
only the top-level manifest assembly, frame/step dictionaries, and TOML I/O.
All wake serialization is in `src/FLOWPanel_replay.jl`:
`_wake_shedding_manifest` (:66), `_wake_manifest_dict` (:284, with the
`PanelParticleWake` branch at :288-303 and `nwakerows` *derived* at :290),
`_construct_wakes_from_manifest` (:450, `PanelParticleWake` at :471-497), and
`_deserialize_wake_shedding` (:505). Phase 3's conversion metadata belongs
there. Note also that `schema_version` (`FLOWPanel_metadata.jl:326`) is written
but never read, so §9.1's versioned conversion block must carry and check its
own `version` key rather than relying on the manifest schema version.

### 13.3 Correction — the existing unknown-tag fallback violates §2.1

`_deserialize_wake_shedding` (`src/FLOWPanel_replay.jl:505`) responds to an
unrecognized tag with `@warn` plus a silent `NoShed()` fallback (:518) — a
tag it does not understand becomes "shed no particles at all", which changes
physics without failing. Conversion-strategy decoding must **not** reuse that
path: an unknown conversion version or type, or a missing smooth parameter, is
a hard error per §9.1. Only a *wholly absent* conversion block selects
`LegacyEdgeJumpConversion()`, and only for backward compatibility.


- **2026-08-01 — Phase 2 APPROVED by the user; amended in §13.** Approval was
  given after a clear-context review that checked the architecture against the
  live source. The review found one substantive simplification and two factual
  corrections, all recorded in §13: the persistent ghost row and its ~20-site
  visibility exclusion (§3.1/§3.2) are unnecessary because conversion runs
  before the row shift, so the ghost is the pre-shift final active row read in
  place (Option B selected, Option A recorded); wake serializers live in
  `FLOWPanel_replay.jl`, not `FLOWPanel_metadata.jl`; and the existing
  `_deserialize_wake_shedding` unknown-tag `NoShed()` fallback must not be
  reused for conversion decoding. Phase 3 is unblocked. Phase 4 remains blocked.

- **2026-07-29 — Phase 2 technically completed; awaiting approval.** Replaced
  the provisional scaffold with the concrete legacy-default/smooth-opt-in
  strategy API, fixed-width surface and root/tip deposition, explicit active
  capacity plus invisible ghost storage, deterministic rank-aware
  reconstruction, bilinear physical-area quadrature, selected
  cancelled-edge handoff, atomic preflight/commit semantics, concrete
  diagnostics, versioned continuation metadata, and exact Phase 3
  implementation/verification maps. No source code was modified and no tests
  or simulations were run. Phase approval remains unchecked and Phase 3 is
  blocked on explicit user approval.
