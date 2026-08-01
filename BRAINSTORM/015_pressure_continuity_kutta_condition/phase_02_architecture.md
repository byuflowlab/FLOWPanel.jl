# Phase 2 — Pressure-Continuity Kutta Architecture

**Status:** APPROVED — 2026-07-29

**Prerequisite:** Phase 1 theory approved on 2026-07-29

**Phase approval:** [x]

**Implementation authorization:** Yes. Phase 3 may begin, subject to the
concrete-typing requirement added at approval.

## 1. Purpose and Phase 1 handoff

This document is the canonical implementation contract for Phase 3. It turns
the approved Phase 1 theory into a typed public API, runtime ownership model,
trial/commit protocol, compatibility boundary, nonlinear-solver contract,
diagnostic schema, and file/test map.

Phase 1 approved two independent configuration axes:

- wake attachment: Route A retains the body-owned rigid
  `TE → TE + Das` transition; Route B uses a TE-anchored live wake panel and
  has no independent rigid `Das`;
- Kutta closure: jump closure enforces \(\gamma=C\mu\) in the coupled body
  solve; pressure closure selects circulation by paired panel-centroid
  pressure continuity.

All four combinations—A/jump, A/pressure, B/jump, and B/pressure—are in the
architecture. Existing behavior is the A/jump default and must branch directly
into the current call sequence. The new runtime is deliberately restricted to
one finite `PanelWake`, one fully paired Dirichlet `RigidWakeBody`, and
`VelocityThroughSources`.

For the governing derivation and the approved experiment order, read
[Phase 1](phase_01_theory.md).

## 2. Public API

### 2.1 Exported strategy types

Phase 3 adds the following exports:

```julia
abstract type AbstractWakeAttachment end
struct RigidTransitionAttachment <: AbstractWakeAttachment end
struct TEAnchoredAttachment <: AbstractWakeAttachment end

abstract type AbstractKuttaClosure end
struct JumpKutta <: AbstractKuttaClosure end
struct PressureContinuityKutta{P,S,F} <: AbstractKuttaClosure
    # provider, primary/fallback strategies, tolerances, and policy
end

abstract type AbstractKuttaPressureProvider end
struct SteadyBernoulliProvider{T} <: AbstractKuttaPressureProvider
    rho::T
end
```

The user-facing constructors are:

```julia
RigidTransitionAttachment()       # Route A and default
TEAnchoredAttachment()            # Route B/B1

JumpKutta()                       # existing γ = Cμ and default
SteadyBernoulliProvider(rho)

PressureContinuityKutta(provider;
    primary=BroydenKutta(),
    fallback=NewtonKrylovKutta(),
    pressure_scale=:auto,
    correction_scale=:auto,
    pressure_tolerance=1e-6,
    correction_tolerance=1e-6,
    on_failure=:error,
    store_diagnostics=true)
```

`BroydenKutta` and `NewtonKrylovKutta` are named nonlinear-strategy
configuration objects used by `PressureContinuityKutta`; callers may construct
them through their `FLOWPanel.`-qualified names. The required new export list
is the attachment, closure, and provider family shown above.

Constructor validation occurs before runtime state is mutated:

- `provider` must implement the pressure-provider interface in section 6;
- explicit pressure and correction scales must be finite and strictly
  positive;
- both tolerances must be finite and strictly positive;
- `on_failure` must be `:error` or `:jump`;
- primary and fallback strategies must satisfy the nonlinear-strategy
  protocol.

### 2.2 Concrete-typing requirement

Every new concrete struct introduced by Phase 3 must be concretely typed after
construction. This applies to public strategies and providers as well as
internal attachment, closure, trial, diagnostics, live-block, solver-control,
metadata, and restart state.

Implementation requirements are:

- parameterize stored providers, nonlinear strategies, numeric types,
  collection element types, and other pluggable state so concrete instances
  have concrete field types;
- do not introduce fields typed as `Any`, an abstract supertype, bare
  `Function`, or an abstract container solely to simplify dispatch;
- use typed callable parameters instead of `Function` fields;
- make optional or variant state an explicit concretely typed parametric
  representation; and
- add tests using `isconcretetype` and field-type inspection for representative
  instances of every new concrete struct.

The abstract API roots (`AbstractWakeAttachment`, `AbstractKuttaClosure`, and
`AbstractKuttaPressureProvider`) remain abstract by design. The requirement
applies to every instantiable subtype and runtime state beneath them.

### 2.3 `simulate!` and `steady!` keywords

Both entry points gain independent typed keywords:

```julia
wake_attachment::AbstractWakeAttachment = RigidTransitionAttachment()
kutta_closure::AbstractKuttaClosure = JumpKutta()
```

These keywords are orthogonal. Neither is inferred from the other, from the
formulation, or from wake options.

The exact legacy pair

```julia
wake_attachment isa RigidTransitionAttachment
kutta_closure isa JumpKutta
```

must branch into the pre-Phase-3 sequence before allocating attachment,
pressure-trial, or nonlinear-solver state. This preserves the current A/jump
operator assembly, wake update order, solver-history behavior, and numerical
results.

`steady!` supports Route A with either jump or pressure closure. It rejects
`TEAnchoredAttachment`, because a steady API cannot infer the convected
downstream row required to define Route B's live panel. Route B steady support
is deferred until there is an explicit prescribed-wake initialization API.

### 2.4 Shared nonlinear coordinate

Both routes and both closures use

\[
\boldsymbol c=C\boldsymbol\mu-\boldsymbol\gamma .
\]

Here \(C\mu\) is the paired upper-minus-lower bound-strength jump and
\(\gamma\) is the corresponding attached/live wake strength. Therefore:

- jump closure is exactly \(\boldsymbol c=\boldsymbol 0\);
- pressure closure solves for \(\boldsymbol c\);
- Route A applies \(c\) through the existing affine attached-transition
  correction;
- Route B applies \(c\) to the explicit live-panel strength.

This coordinate is the common composition seam. There is no separate Route B
pressure unknown and no post-solve reassignment of \(\gamma\).

### 2.5 Diagnostics accessor

The public accessor is:

```julia
kutta_diagnostics(kutta_closure)
```

For `PressureContinuityKutta`, it returns the accepted per-step diagnostics in
physical-step order. When `store_diagnostics=false`, it returns an empty
collection without disabling failure reporting. For `JumpKutta`, it may return
an empty collection; Route B startup diagnostics belong to the pressure
closure that requested them.

## 3. Supported domain and deterministic validation

### 3.1 Phase 3 support boundary

Any non-default attachment or closure is accepted only when all of the
following are true:

1. there is exactly one body;
2. it is one Dirichlet `RigidWakeBody`;
3. every shedding edge has one valid upper/lower panel pair;
4. the wake is finite (`semiinfinite_wake == false`);
5. unsteady operation has exactly one `PanelWake` associated with the body;
6. the formulation is `VelocityThroughSources`;
7. the selected body solver is `Backslash`, `KrylovSolver`, or a fully
   constructed `FGSSolver`;
8. `bound_strength_rlx == 1`.

`steady!` applies the same body, pairing, formulation, relaxation, and solver
checks but allows no external wake and only Route A.

Validation rejects the whole configuration before changing body, wake,
provider, or solver state. Error messages identify the unsupported object and
the required Phase 3 domain.

### 3.2 Explicitly legacy-only configurations

The following remain available only through the exact default A/jump path:

- `PanelParticleWake`;
- coupled or multi-body solves;
- unpaired, partially paired, or ambiguous shedding edges;
- semi-infinite attached wakes;
- Neumann and other non-Dirichlet bodies;
- formulations other than `VelocityThroughSources`;
- solver types outside the table below.

The new runtime must not silently approximate, ignore, or reinterpret these
configurations.

### 3.3 Option conflicts

Route B rejects every `set_Das_*` option rather than ignoring it:

- non-`NaN` `set_Das_eta_kinematic`;
- non-`NaN` `set_Das_eta_freestream`;
- nonzero `set_Das_min_kinematic_displacement`;
- non-default `set_Das_kinematic_arc`;
- `set_Das_refresh=true`.

The rejection occurs even though Route B never reads `Das`; this prevents a
run from claiming to use options that cannot affect its geometry.

Any non-default attachment or closure rejects `bound_strength_rlx != 1`.
Post-solve relaxation changes \(C\mu\) after the coupled closure and therefore
invalidates both jump and pressure consistency.

### 3.4 Compatibility matrix

| Configuration | `Backslash` | `KrylovSolver` | constructed `FGSSolver` |
| --- | ---: | ---: | ---: |
| A/jump | Existing path | Existing path | Existing path |
| A/pressure | Supported | Supported | Supported |
| B/jump | Supported | Supported | Supported |
| B/pressure | Supported | Supported | Supported |

`FGSSolver{Nothing}` is rejected. It carries options but has no executable
Fast Gauss–Seidel object.

| Entry point | A/jump | A/pressure | B/jump | B/pressure |
| --- | ---: | ---: | ---: | ---: |
| `simulate!` | Supported, legacy path | Supported | Supported | Supported |
| `steady!` | Supported, legacy path | Supported | Rejected | Rejected |

## 4. State ownership

The implementation separates persistent physical state from ephemeral
attachment and nonlinear-trial state.

### 4.1 Persistent owners

`PanelWake` remains the sole owner of persistent wake data:

- wake nodes and panel strengths;
- old-wake topology;
- live-block reservation metadata;
- physical-step identifiers for each deposited block;
- current wake count, velocity, overflow, and convection state.

`RigidWakeBody` continues to own:

- surface geometry, strength, potential, and reconstructed fields;
- paired shedding connectivity;
- Route A's rigid `TE → TE + Das` transition;
- existing affine correction storage used by Route A.

The body does not become the owner of persistent Route B wake panels.

### 4.2 Internal runtime states

Phase 3 introduces three separate internal state families:

1. **attachment runtime** — the active route, paired-edge map, geometry
   version, ephemeral body-facing attachment geometry, live-block view, and
   operator-cache version;
2. **closure runtime** — frozen pressure/correction scales, nonlinear strategy
   workspace, best finite trial, diagnostics, and provider state;
3. **trial snapshot** — restorable body fields, attachment state, provider
   state, wake state, and solver history for one physical step.

These states are initialized only for a non-default configuration. They are
not stored in monitor context and they do not transfer ownership of wake
arrays to the body.

The body receives only an ephemeral attachment-geometry cache sufficient for
its direct, Krylov, and FGS operator paths. The cache is invalidated by a
geometry-version change and is never serialized as authoritative wake state.

### 4.3 Route A

Route A continues to use the body-owned panel from `TE` to `TE + Das`.
Attachment geometry is exactly the current geometry. Pressure trials vary the
existing affine correction \(c\), with operator mode excluding the affine
constant and physical mode including it.

The accepted Route A state is committed through the existing body and wake
sequence. Pressure closure does not change who owns `Das` and does not make
Route A `Das`-free.

### 4.4 Route B/B1

Route B defines a single live panel per paired TE edge:

- its upstream endpoint is the current body TE;
- its downstream endpoint is the first convected downstream wake row;
- it never reads the body's user `Das` values;
- its strength is \(\gamma=C\mu-c\) inside the coupled body solve.

The `PanelWake` reserves the live block. A wake-source view exposed to ordinary
old-wake influence excludes that block. The body attachment operator supplies
the live-panel influence during both the coupled solve and accepted exterior
field reconstruction. This exclusion/substitution rule prevents double
counting.

Acceptance overwrites the reserved live-panel strength exactly once. At the
end of the physical interval, topology advancement:

1. reclassifies the accepted live block as old wake;
2. retains its accepted strength without recomputing circulation;
3. reserves the next live slot;
4. assigns the next physical-step identifier.

Topology advancement is not a second commit and does not call the closure.

### 4.5 Live-block generalization and future seams

Wake insertion is generalized around a block containing one or more panels
per physical step. Every panel in a block shares one physical-step identifier;
subpanels are never counted as extra timesteps.

B1 uses a one-panel block. The block interface must also accept:

- a future B2 \(N_s\)-panel block;
- equal-strength B2 data for the internal-edge cancellation/null gate;
- later birth-time interpolation data without changing physical-step
  bookkeeping.

A separate deposition-strategy seam owns node placement and within-block
strength assignment. B1 is the only implemented strategy in Phase 3. The seam
permits B2 and the B3 Katz–Plotkin comparator later without changing closure,
solver, or ownership contracts.

## 5. Physical-step ordering

### 5.1 Route A

The current A/jump call order is unchanged. For A/pressure, each physical step
does:

1. update frames, body geometry, control points, and prescribed velocities;
2. update Route A attachment geometry and its operator version;
3. freeze wake/external influence and all trial state;
4. run the pressure closure;
5. commit the accepted body/attachment/provider/solver state once;
6. continue the existing wake propagation, shedding, monitor, and output
   sequence.

### 5.2 Cold Route B startup

Route B requires a convected downstream row before a live panel exists. A cold
run therefore has an explicit seed rule:

1. at \(t_0\), solve with jump closure and no live panel;
2. save the accepted \(\gamma^0=C\mu^0\);
3. record a deterministic `:startup_jump` diagnostic when pressure closure
   was requested;
4. convect the TE seed through the first physical interval;
5. at \(t_1\), activate the live-panel coupling between the current TE and
   that convected row.

No artificial `Das`, guessed panel, or pressure solve at \(t_0\) is allowed.
The first Route B pressure-continuity solve is at \(t_1\).

A valid warm start with restored Route B topology, live-block state, prior
circulation, and closure metadata resumes directly and does not repeat cold
startup. An incomplete or incompatible Route B restart is rejected rather
than silently reseeded.

### 5.3 Route B accepted-step order

After startup, each Route B physical step follows:

1. restore or construct the current reserved live block;
2. exclude it from the old-wake source view;
3. update the body-facing attachment geometry cache;
4. freeze the step snapshot and nonlinear scales;
5. solve jump directly or execute pressure trials;
6. commit the accepted body, live-panel strength, provider, and solver history
   once;
7. reconstruct and expose the accepted field with the same no-double-counting
   source split;
8. run wake convection/rollup, monitors, and output in their existing relative
   order;
9. advance the accepted live block into old-wake storage and reserve the next
   block, without recomputing or recommitting circulation.

The implementation must document the exact insertion point relative to the
current wake propagation functions, but it may not reorder existing monitor
dependencies.

## 6. Pressure-provider interface

The provider extension points are:

```julia
pressure_requirements(provider)
evaluate_pressure!(pressure, provider, trial_state)
commit_pressure!(provider, accepted_trial)
rollback_pressure!(provider)
```

`pressure_requirements` declares all data that must be reconstructed before a
trial evaluation—for example surface velocity, potential, potential history,
velocity gradient, or acceleration. Initialization validates that the runtime
can supply every requirement.

`trial_state` is a read-only logical view of the complete restored-and-solved
trial. It contains the body, paired-edge map, current kinematic/freestream
data, attachment state, and only the provider-declared completed fields.
Providers must not solve the body, alter the wake, register monitor data, or
advance time history.

The behavioral contract is:

- the same trial state produces the same pressure;
- evaluation mutates only the caller-provided `pressure` workspace and
  provider scratch explicitly designated as trial-local;
- `commit_pressure!` is called once for the accepted trial;
- `rollback_pressure!` returns provider trial state to the frozen pre-step
  state after every rejected or failed trial;
- pressure gauge behavior must cancel in paired differences or be held fixed.

### 6.1 `SteadyBernoulliProvider`

`SteadyBernoulliProvider(rho)`:

- requires only the completed relative exterior surface velocity;
- computes pressure at the paired panel centroids using steady Bernoulli;
- has no time-history state;
- implements commit and rollback as no-ops;
- does not wrap, call, or reuse the mutable `PressureBernoulli` monitor.

The pressure closure owns reconstruction and calls it exactly once per body
solve. The provider consumes that reconstruction; it must not trigger a second
accumulating field evaluation.

Unsteady Bernoulli, `PressureLaplace`, and monitor adapters remain extension
points, not Phase 3 deliverables.

## 7. Trial transaction

### 7.1 Frozen snapshot

Before the initial clean trial, snapshot every value a failed or repeated
trial could mutate:

- body strengths, potential, velocity, gradients/vorticity when present, and
  `fields`;
- prescribed freestream and kinematic velocity;
- Route A affine correction flags/storage or Route B live attachment state;
- all persistent wake nodes, strengths, counts, velocities, overflow state,
  source-view bounds, and live-block metadata that trial code can reach;
- provider state;
- `FGSSolver` solution history, saved-count, and projection base;
- other solver work/history whose mutation would change a repeat trial.

Old-wake influence and prescribed external fields are computed from this clean
state and remain fixed through the nonlinear solve.

### 7.2 One residual evaluation

For a trial correction \(c\):

1. restore the full frozen snapshot;
2. call `rollback_pressure!`;
3. restore attachment and solver history;
4. apply \(c\), giving \(\gamma=C\mu-c\) through the chosen route;
5. solve the coupled linear body problem;
6. capture the actual inner-solver convergence status;
7. if the inner solve is valid, reconstruct the complete exterior relative
   surface velocity once;
8. call `evaluate_pressure!`;
9. form the paired residual
   \(r_e=p_{\mathrm{upper}(e)}-p_{\mathrm{lower}(e)}\);
10. return the residual, norms, inner status, and a restorable trial record.

A non-finite solve, field, pressure, or residual is a rejected trial. It may
participate in backtracking/restart accounting but can never become the best
or accepted trial.

Repeated evaluation of the same \(c\) from the same snapshot must be
deterministic to the tolerances appropriate for the selected body solver.

### 7.3 Acceptance

Convergence selects an accepted trial record; it does not rerun that trial.
Commit then:

1. restores the accepted body and attachment values;
2. writes Route B's reserved live strength once when applicable;
3. calls `commit_pressure!` once;
4. saves the accepted FGS solution once when applicable;
5. appends and serializes diagnostics;
6. releases the trial snapshot;
7. permits physical wake/time advancement.

Rejected trials never call `save_solution!`, never shift provider history, and
never deposit wake circulation.

### 7.4 Failure and explicit jump fallback

`on_failure=:error` is the default. If convergence fails, an inner solve is
invalid, all trials become non-finite, or commit cannot complete:

1. restore the complete pre-step snapshot;
2. call `rollback_pressure!`;
3. leave wake topology and all histories unchanged;
4. throw `KuttaConvergenceError` carrying the final diagnostics and cause.

The last nonlinear trial is never left installed.

`on_failure=:jump` is an explicit physical fallback. It:

1. restores the pre-step snapshot;
2. performs a fresh coupled solve with \(c=0\);
3. requires that fresh inner solve and field reconstruction to succeed;
4. commits that fresh jump solution once;
5. emits a warning;
6. records a fallback diagnostic disposition/status.

It must not commit the best pressure trial and then relabel it as jump. If the
fresh jump solve fails, state is restored and `KuttaConvergenceError` is
thrown.

## 8. Attachment-aware body solvers

All new combinations use one attachment-aware linear operator. The difference
between jump and pressure closure is the chosen \(c\), not a different body
operator.

### 8.1 `Backslash`

- assemble the attachment-aware matrix after the attachment geometry version
  changes;
- factor it once per geometry version;
- reuse that factorization for all nonlinear trials at the step;
- keep the affine correction out of operator assembly and put it in the
  right-hand side/physical state as required;
- expose a successful direct-solve status to the outer convergence gate.

### 8.2 `KrylovSolver`

- use the same attachment-aware matrix-free product as direct assembly;
- apply Route A affine or Route B live-panel terms consistently in operator
  and right-hand-side modes;
- reset each nonlinear trial to its frozen initial state as required for
  determinism;
- report actual Krylov convergence, iteration count, and residual status to
  the outer solver;
- never allow an unconverged inner result to satisfy outer convergence.

### 8.3 `FGSSolver`

`FGSSolver` gains internal trial controls, conceptually separating:

- projection from frozen physical-time history;
- one trial solve;
- accepted-history save.

Every trial projects from the same frozen pre-step history. Trial solves do not
call `save_solution!`. The accepted solution calls it exactly once.
`solution_history`, `solution_history_nsaved`, and the projection base are
restored after every rejected trial.

A constructed `FGSSolver` must use the same body-facing attachment cache in
its direct/FMM influence route. `FGSSolver{Nothing}` is rejected during
configuration validation.

## 9. Nonlinear pressure closure

### 9.1 Frozen scaling

The initial clean finite trial establishes and freezes both scales for the
physical step.

With `pressure_scale=:auto`,

\[
P_s=\tfrac12\rho U_s^2 ,
\]

where \(U_s\) is the reference speed derived from the initial clean trial.
With `correction_scale=:auto`,

\[
C_s=U_s L_s ,
\]

where \(L_s\) is the median characteristic length of all paired TE panels.

The exact reference-speed reduction used for `U_s` must be one deterministic
global reduction over the completed relative surface velocities and be
recorded in diagnostics. It must not change between trials. A zero,
non-finite, or otherwise degenerate automatic pressure or correction scale
throws an `ArgumentError` requesting an explicit positive scale. Explicit
scales are also frozen for the step.

The nonlinear solver operates on

\[
\hat r=r/P_s,\qquad \hat c=c/C_s.
\]

### 9.2 Safeguarded inverse Broyden

The primary method is dense global inverse Broyden on the scaled vectors. Its
defaults are:

- at most 30 outer iterations;
- at most 12 backtracks per proposed update;
- at most two restarts;
- maximum scaled correction step
  \(\lVert\Delta\hat c\rVert_\infty=2\);
- minimum line-search factor \(2^{-20}\).

It uses Armijo backtracking on the scaled residual norm, rejects non-finite
trials, limits the scaled step before line search, and restarts its inverse
Jacobian approximation after repeated rejected steps or an ill-conditioned
update. It retains the best finite restorable trial.

The inverse approximation is dense across all paired edges; independent
edgewise Broyden updates are not permitted.

### 9.3 Finite-difference Newton–Krylov fallback

Recoverable Broyden stagnation, exhausted restarts, or a persistently unusable
inverse update invokes the fallback from the best finite trial.

Defaults are:

- at most 12 Newton steps;
- at most 30 Krylov iterations per Newton step;
- relative directional-difference step
  \(\sqrt{\operatorname{eps}(\mathrm{Float64})}\);
- the same scaled-step limit, backtracking, non-finite rejection, snapshot
  restoration, and accepted-trial rules as Broyden.

Directional Jacobian products are finite differences of complete residual
transactions. They therefore count body solves and may not reuse a mutated
trial state. Failure of the fallback passes through the configured
`on_failure` policy.

### 9.4 Convergence gates and norms

Acceptance requires all three gates simultaneously:

1. pressure gate:
   \(\lVert r\rVert_\infty/P_s \leq\mathtt{pressure\_tolerance}\);
2. correction-step gate:
   \(\lVert\Delta c\rVert_\infty/C_s
   \leq\mathtt{correction\_tolerance}\);
3. the accepted trial's inner body solve reports convergence.

The initial trial may satisfy the correction-step gate with a defined zero
step only if the pressure and inner gates also pass.

For observability, also report the edge-length-weighted pressure norm

\[
\lVert r\rVert_{2,\ell}
=
\sqrt{\frac{\sum_e \ell_e r_e^2}{\sum_e \ell_e}},
\]

both dimensional and divided by \(P_s\). Edge lengths and panel pairing are
frozen for the trial sequence.

## 10. Diagnostics, errors, metadata, and restart

### 10.1 `KuttaDiagnostics`

Store one `KuttaDiagnostics` record per accepted physical step when enabled.
The schema records at least:

- `status`;
- method used at termination;
- outer iteration count;
- total body-solve count;
- backtrack count;
- restart count;
- dimensional and scaled pressure infinity norms;
- dimensional and scaled edge-weighted pressure L2 norms;
- dimensional and scaled accepted correction step;
- frozen `pressure_scale`, `correction_scale`, `U_s`, and `L_s`;
- inner-solver convergence/status and iteration information;
- whether fallback was entered;
- attachment route;
- startup disposition such as `:startup_jump`;
- failure/fallback disposition, including explicit jump fallback.

Statuses must distinguish ordinary pressure convergence, startup jump,
pressure failure followed by committed jump, and thrown failure. The thrown
`KuttaConvergenceError` includes the uncommitted failure record even when
diagnostic storage is disabled.

### 10.2 Metadata

Run metadata serializes:

- concrete attachment and closure names;
- pressure-provider name and serializable parameters;
- nonlinear strategy names and effective safeguards;
- tolerances, failure policy, and frozen-scale policy;
- per-step accepted diagnostics;
- accepted closure state needed to resume;
- Route B physical-step/live-block identifiers and prior accepted
  circulation needed for warm start.

Metadata must describe only committed state. Rejected trial coordinates,
provider scratch, and ephemeral body attachment caches are not restart state.

### 10.3 Warm-start contract

A warm start validates that saved attachment, closure, provider, wake topology,
pairing, and solver-history state are compatible with the requested run.
Accepted state and diagnostics are restored before the first resumed solve.

For Route B, the restart must identify which block is live, which blocks are
old, their physical-step identifiers, and the prior accepted circulation.
When that state is complete, resumption does not run `:startup_jump`.
Missing or inconsistent new-format state throws an informative error; it does
not guess a live block or reuse `Das`.

Legacy A/jump warm starts retain their current behavior and do not require new
closure state.

## 11. Phase 3 file map

| Concern | Primary files | Required seam |
| --- | --- | --- |
| Public types, exports, formulation integration | `src/FLOWPanel.jl`, `src/FLOWPanel_formulation.jl` | strategy definitions/exports, validation, common \(c\) coordinate, attachment-aware products |
| Body attachment geometry | `src/FLOWPanel_liftingbody.jl` | ephemeral cache, pairing checks, Route A affine and Route B live-panel influence |
| Direct, Krylov, and FGS solves | `src/FLOWPanel_solver.jl` | geometry-version cache, actual inner status, frozen FGS trial/history controls |
| Wake storage and topology | `src/FLOWPanel_wake.jl` | live-block reservation, old-wake view, physical-step IDs, B2/deposition seams |
| Physical-step transaction | `src/FLOWPanel_simulate.jl` | typed keywords, legacy fast branch, startup, snapshots, solve/commit/advance order |
| Provider and pressure reconstruction | pressure/postprocess module selected in Phase 3, plus `src/FLOWPanel_postprocess.jl` as needed | provider traits, complete relative surface velocity, no monitor reuse |
| Diagnostics and run manifests | `src/FLOWPanel_metadata.jl` and the selected closure/diagnostics module | serialization of configuration, accepted state, diagnostics |
| Warm start | `src/FLOWPanel_warmstart.jl` | compatibility validation and Route B live-state restoration |

The authoritative include order in `src/FLOWPanel.jl` determines whether new
definitions live in existing files or in a narrowly scoped new module. Phase 3
may choose that mechanical placement, but it may not change the contracts in
this document without returning to Phase 2 approval.

## 12. Phase 3 mechanical test map

Phase 3 adds focused unit coverage for:

1. all four attachment/closure combinations;
2. exact default dispatch into the legacy A/jump path and numerical
   regression;
3. fully paired-edge validation and rejection of unpaired/ambiguous edges;
4. finite-wake, one-body, Dirichlet, formulation, and solver validation;
5. Route B rejection of every `set_Das_*` option;
6. rejection of `bound_strength_rlx != 1` for non-default configurations;
7. Route B cold `:startup_jump`, activation at \(t_1\), and warm-start bypass;
8. live-panel exclusion from old-wake influence and accepted reconstruction,
   demonstrating no double counting;
9. exactly one Route B live-strength write and one topology advancement per
   accepted physical step;
10. rollback after inner failure, non-finite pressure, line-search rejection,
    Broyden failure, fallback failure, and commit failure;
11. fresh `c=0` solve and diagnostic status for `on_failure=:jump`;
12. repeated-trial determinism for body, wake, provider, and solver history;
13. Broyden-to-Newton–Krylov fallback from the best finite trial;
14. simultaneous pressure, correction-step, and inner-solver convergence
    gates;
15. automatic-scale freezing and degenerate-scale rejection;
16. `Backslash`, `KrylovSolver`, and constructed `FGSSolver` agreement;
17. `FGSSolver{Nothing}` rejection and frozen FGS projection/history with one
    accepted save;
18. direct/FMM attachment-operator and accepted-field agreement;
19. steady Route A pressure support and TE-anchored steady rejection;
20. metadata round trip and warm-start restoration of accepted closure,
    diagnostics, live-block, physical-step, and solver-history state;
21. B2 seam tests showing that a multi-panel block has one physical-step ID
    and that equal-strength internal edges cancel without changing the B1
    ownership/commit protocol.
22. concrete-type tests for every new public and internal concrete struct,
    including rejection of `Any`, abstractly typed, bare-`Function`, and
    abstract-container fields.

The nearest repository suites are:

- `test/runtests_unit_solver.jl`;
- `test/runtests_unit_liftingbody.jl`;
- `test/runtests_unit_wake.jl`;
- `test/runtests_unit_simulate.jl`;
- `test/runtests_unit_postprocess.jl`;
- `test/runtests_unit_fgs_history.jl`;
- `test/runtests_unit_warmstart.jl`.

These are mechanical architecture tests. Phase 3 must not claim physical
validation from them.

## 13. Phase 4 gates retained

Phase 4 retains the physical gates approved in Phase 1:

- pressure continuity must demonstrate sharp-edge refinement/convergence;
- Route B must demonstrate timestep and wake-resolution convergence;
- B1 precedes the B2 equal-strength null gate and \(N_s=2,4,8\) persistent
  subpanel study;
- B3 remains the Katz–Plotkin comparator;
- Route A is never described as `Das`-free;
- Route B is described as independent of rigid `Das`, not independent of
  newborn-wake discretization.

## 14. Decision record

The following Phase 2 decisions are technically complete and await explicit
user approval:

1. The public configuration has independent typed
   `wake_attachment` and `kutta_closure` keywords on `simulate!` and
   `steady!`, defaulting to Route A/jump.
2. The exact default pair branches directly into the unchanged legacy path.
3. Both routes use \(c=C\mu-\gamma\); jump closure is exactly \(c=0\).
4. Phase 3 supports only one fully paired finite-wake Dirichlet
   `RigidWakeBody`, one `PanelWake`, and `VelocityThroughSources`.
5. Route A retains the body-owned rigid transition; Route B derives its live
   endpoint from the first convected row and rejects all `set_Das_*` options.
6. `PanelWake` owns persistent live/old wake storage. The body receives only
   an ephemeral attachment cache used consistently by direct, Krylov, FGS,
   and accepted-field reconstruction paths.
7. Route B's live block is excluded from old-wake source views and represented
   by the attachment operator exactly once.
8. Cold Route B startup is a deterministic jump solve at \(t_0\), convection
   through the first interval, and live-panel activation at \(t_1\).
9. Live insertion is block-based with physical-step identifiers, implementing
   B1 now while preserving B2 and deposition-strategy seams.
10. Pressure providers use explicit requirements, evaluate, commit, and
    rollback hooks. `SteadyBernoulliProvider` consumes completed relative
    surface velocity and does not reuse `PressureBernoulli`.
11. Every trial is a full restore/apply/solve/reconstruct/evaluate
    transaction. Only the accepted state advances wake, provider, or FGS
    history.
12. Backslash caches per geometry version, Krylov exposes actual inner status,
    and FGS projects from frozen history and saves once. `FGSSolver{Nothing}`
    is rejected.
13. Automatic pressure and correction scales are frozen from the initial clean
    trial; degenerate automatic scales require explicit positive values.
14. Safeguarded dense inverse Broyden is primary; finite-difference
    Newton–Krylov is the recoverable-stagnation fallback with the limits in
    section 9.
15. Convergence requires pressure, correction-step, and inner-solver gates
    simultaneously, with both infinity and edge-weighted L2 pressure norms
    reported.
16. Default failure restores pre-step state and throws
    `KuttaConvergenceError`; explicit `on_failure=:jump` restores and commits
    a fresh \(c=0\) solve with warning and diagnostic disposition.
17. Accepted state, configuration, diagnostics, Route B topology, and required
    solver/provider history are serialized for deterministic warm starts.
18. `steady!` supports Route A pressure closure only; TE-anchored steady solves
    await an explicit prescribed-wake initialization API.
19. Every new instantiable public or internal struct must be concretely typed;
    the abstract API roots remain abstract by design.

## 15. Acceptance gate and progress log

Phase 2 is approved. Phase 3 may implement this architecture, including the
binding concrete-typing requirement in section 2.2. Phase 3 remains a separate
approval gate and must not expand the approved scope.

### Phase 2 progress and approval log

- **2026-07-29 — Phase seeded.** Phase 1 approved; architecture work was ready
  but had not started.
- **2026-07-29 — Architecture completed by agent.** Recorded the exact typed
  API, four-combination compatibility boundary, ownership and live-block
  model, Route B startup and no-double-counting rules, provider transaction,
  solver hooks, nonlinear safeguards, diagnostics/metadata/restart contract,
  and Phase 3 file/test seams. Status is technically complete and awaiting
  explicit user approval. Phase 3 remains blocked.
- **2026-07-29 — Phase 2 explicitly approved by user, with concrete typing
  added.** Every new instantiable public or internal struct must be concretely
  typed after construction; representative concrete-type and field-type tests
  are required. Phase 3 is authorized and ready to begin.
