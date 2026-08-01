# Phase 3 — Implementation

**Status:** APPROVED — 2026-07-29 (post-remediation, per the user-approved
review-then-remediate plan)

**Scope amendment (user, 2026-07-29):** the first implementation supports only
the `Backslash` body solver. `KrylovSolver`/`FGSSolver` are validation-rejected
with an explicit "deferred" message; their trial hooks (architecture §8.2-8.3)
and test categories 16-17 (beyond the rejection tests) are deferred until the
Backslash implementation is validated.

**Prerequisite:** approved Phase 2 architecture

**Phase approval:** [x]

## Prior-phase handoff

Phase 1 approved independent wake-attachment and Kutta-closure axes, preserved
legacy defaults, Route B experiment order B1 then B2, safeguarded global
Broyden with finite-difference Newton–Krylov fallback, steady Bernoulli as the
initial pressure provider, and frozen-trial/one-time-commit semantics.

Phase 2 must supply the exact approved API, state ownership, provider traits,
solver hooks, diagnostics, failure behavior, compatibility rules,
implementation map, and test seams. Read the hub and the approved
[Phase 2 architecture](phase_02_architecture.md) before starting. Do not infer
missing design decisions directly from Phase 1.

Phase 2 approval also requires every new instantiable public or internal
struct to be concretely typed after construction. Parameterize stored
providers, strategies, callables, numeric types, and containers; do not use
`Any`, abstract field/container types, or bare `Function` fields as
implementation shortcuts. Add representative `isconcretetype` and field-type
tests for every new concrete struct.

## Scope

Implement only the approved Phase 2 design. Preserve all existing defaults and
existing formulations. Add the approved public configuration, state and trial
machinery, pressure-continuity nonlinear closure, Route B B1 representation,
diagnostics, and deterministic failure paths. Add B2 or B3 only if the
approved Phase 2 architecture explicitly places them in Phase 3.

No unapproved redesign, new pressure law, higher-order panel kernel, or
validation conclusion belongs in this phase.

## Governing decisions

- Treat the approved Phase 2 architecture as the implementation contract.
- Preserve existing defaults and keep both configuration axes independent.
- Freeze trial state and commit physical and solver histories exactly once.
- Keep validation experiments and conclusions in Phase 4.

## Deliverables

- Source changes matching the approved architecture and file map.
- Focused mechanical/unit tests for new API, default preservation, state
  freezing, one-time commit, compatibility checks, and failure paths.
- Documentation of changed locations, defaults, and user-facing behavior.
- A local progress log containing the exact verification commands and results.

## Acceptance gate

Phase 3 is complete only when its focused mechanical tests pass, the scoped
diff matches the approved architecture, legacy defaults remain unchanged, and
the user explicitly approves the implementation. Phase 4 remains blocked
until that approval.

## Implementation record (2026-07-29)

### Changed locations

- **`src/FLOWPanel_kutta.jl` (new, included between `formulation` and
  `simulate_monitors`)** — all public strategy types
  (`RigidTransitionAttachment`, `TEAnchoredAttachment`, `JumpKutta`,
  `PressureContinuityKutta{P,S,F}`, `SteadyBernoulliProvider{T}`, unexported
  `BroydenKutta`/`NewtonKrylovKutta`), the pressure-provider interface
  (`pressure_requirements`/`evaluate_pressure!`/`commit_pressure!`/
  `rollback_pressure!` with `KuttaTrialView`), `KuttaDiagnostics{TF}`,
  `KuttaConvergenceError{TF}`, `kutta_diagnostics`, configuration validation
  (`_validate_kutta_configuration`), the runtime
  (`KuttaRuntime`/`KuttaSnapshot`/`KuttaTrialRecord`), the trial transaction
  (`_kutta_trial!`), frozen scales, the safeguarded inverse-Broyden primary and
  FD Newton–Krylov fallback, the step orchestrator (`_kutta_step!`),
  commit/failure paths, Route B geometry
  (`_update_TE_route_b!`, `_kutta_update_route_b_operators!`,
  `_kutta_advance_topology!`), and metadata/warm-start helpers.
- **`src/FLOWPanel.jl`** — include line; exports
  `AbstractWakeAttachment, RigidTransitionAttachment, TEAnchoredAttachment,
  AbstractKuttaClosure, JumpKutta, PressureContinuityKutta,
  AbstractKuttaPressureProvider, SteadyBernoulliProvider, kutta_diagnostics,
  KuttaDiagnostics, KuttaConvergenceError`.
- **`src/FLOWPanel_simulate.jl`** — pure code-motion stage split of
  `_steady_aerodynamics!` into `_sa_collect`,
  `_sa_reset_freestream_kinematic!`, `_sa_wake_influence!`,
  `_sa_body_influence!`, `_sa_half_jump!` (identical statements/order; the
  legacy function calls them; user-approved); `wake_attachment`/`kutta_closure`
  keywords on `simulate!` and `steady!`; early whole-configuration validation
  before any state mutation; the legacy fast branch (`_is_legacy_kutta` —
  the exact default pair allocates no Kutta state and takes the pre-existing
  call sequence verbatim); the in-loop solve branch; the Route B end-of-step
  advancement hook; metadata pass-through.
- **`src/FLOWPanel_wake.jl`** — `PanelWake` gains live-block reservation
  metadata `live_rows::Array{Int,0}`/`live_step_id::Array{Int,0}` (0/-1
  legacy); `_n_wake_source_rows`; the old-wake source view
  (`FastMultipole.get_n_bodies`, `global_to_matrix_index`,
  `matrix_to_global_index`) excludes the reserved live block. With
  `live_rows[]==0` every expression reduces exactly to the pre-existing
  arithmetic.
- **`src/FLOWPanel_metadata.jl`** — optional `kutta` manifest table and
  per-step `step.kutta` records (accepted `c`, accepted circulation, route,
  live-block identifiers, status).
- **`src/FLOWPanel_warmstart.jl`** — `wake_attachment`/`kutta_closure`
  keywords; `_kutta_warmstart_restore!` (configuration compatibility
  validation; reinstalls the committed correction BEFORE the end-of-step
  replay so the replayed `shed_wake!` deposits γ = Cμ − c; restores and then
  advances Route B live-block metadata). Incomplete or mismatched saved state
  throws `ArgumentError`; legacy warm starts are untouched.
- **Tests** — new `test/runtests_unit_kutta.jl` (602 assertions) and
  `test/runtests_unit_kutta_routeb.jl` (59 assertions), registered in
  `test/runtests.jl`; shared fixture `make_dirichlet_diamond_body` in
  `test/test_helpers.jl` (diamond wedge with 2 chordwise panels per side —
  a flat plate with a coplanar wake is degenerate for pressure continuity:
  the correction has exactly zero leverage).

### Defaults and user-facing behavior

- Both new keywords default to the legacy pair
  (`RigidTransitionAttachment()` + `JumpKutta()`), which branches into the
  unchanged pre-Phase-3 sequence before allocating anything; explicit-default
  and omitted-kwargs runs are bitwise identical (tested), and all pre-existing
  unit suites pass unchanged.
- Non-default configurations enforce the §3.1 support boundary with
  `ArgumentError`s before any state mutation: one fully paired finite-wake
  Dirichlet `RigidWakeBody`, one `PanelWake`, `VelocityThroughSources`,
  `Backslash` (amendment), `bound_strength_rlx == 1`; Route B rejects every
  `set_Das_*` option; `steady!` rejects `TEAnchoredAttachment`.
- `kutta_diagnostics(closure)` returns per-accepted-step
  `KuttaDiagnostics` records; `on_failure=:error` restores the pre-step state
  and throws `KuttaConvergenceError`; `on_failure=:jump` restores, commits a
  fresh c = 0 solve (bitwise the legacy trajectory — tested), warns, and
  records the `:jump_fallback` disposition.

### Mechanical design choices within the architecture's delegation (§11)

1. **The correction channel splits linear-operator and affine parts, and is
   backend-complete** (revised 2026-07-29, follow-up session — user request):
   trials solve on the cached factorization with `c` entering only through
   the right-hand side (`body.potential .-= W*c`, W from the existing
   `_assemble_W!`); trial and committed velocity reconstruction runs the
   influence passes with the correction INACTIVE through **any backend** (the
   proportional γ = Cμ part flows through both the index and FMM buffer
   paths), and the affine −c part — sourced only by the M trailing-edge
   transition strips, one wake row of panels — is completed by the exact
   direct add-on `_add_affine_attached_velocity!` (O(M×N_targets); the same
   per-panel ∓c/2 decomposition and triangle construction as the index-path
   `_induced_wake`, verified against the production full-minus-suppressed
   unit-strip influence to 1e-11). The add-on is applied to the body
   centroids in every trial and to the body plus wake probes at commit. No
   FastMultipole buffer or multipole-moment changes; no backend restriction.
2. **Route B live block = `PanelWake` panel row 1 plus an ephemeral
   attachment cache**: node row 1 is anchored at the bare TE and
   `body.Das .= nodes[:,2,·] − TE` is rebuilt every step, making the body's
   existing attached-transition influence (`_induced_wake`, both index and
   buffer paths) exactly the live panel. The live row is excluded from
   old-wake source views (no double counting — verified against brute force);
   `propagate!`/`shed_wake!` are unchanged, and the row-1 deposit after
   advancement is the reserved next slot's warm-start prior (excluded from
   source views, overwritten at the next commit). Cold startup solves under
   `suppress_attached_wake` (genuinely no live panel) and records
   `:startup_jump`.
3. **Broyden initialization** uses one secant probe along the uniform
   direction to set a conservative inverse-Jacobian scale (Phase 1 §6's
   "conservative initial inverse-Jacobian scale"); with it the primary
   converges in ~2 iterations / 4 body solves per step on the test fixture,
   and the identity-initialized fallback path reproduces the same accepted
   correction to 1e-6 (tested).
4. **Gate mechanics refinement (§9.4 interpretation):** at each outer
   iteration, if the current accepted iterate already passes the pressure and
   inner gates and the *proposed* capped step is below the correction
   tolerance, the current iterate is accepted (a further step cannot improve
   within tolerance — without this, an exactly-converged iterate at the
   round-off floor can never satisfy an Armijo decrease and the step gate
   simultaneously). Symmetrically, the fallback accepts an already-passing
   base point as its initial trial with a defined zero step.
5. **Scale freezing** costs one extra clean solve+reconstruction per step
   (the scaled drivers need C_s before their first trial). `U_s` = serial RMS
   of the relative surface speed over all centroids; `L_s` = median
   `sqrt(panel area)` over the paired TE panels.
6. **Warm start** restores the committed correction onto the body (which also
   seeds the resumed nonlinear solve) and the Route B live-block metadata;
   per-step accepted records live in the metadata TOML. The in-memory
   `diagnostics` vector of the freshly constructed closure object is not
   repopulated from disk (the records remain available in the TOML).
7. `steady!` gains no `formulation` keyword (it is implicitly
   `VelocityThroughSources`, the only supported formulation).
8. The `_kutta_gates_pass`/commit machinery treats a non-finite jump or
   startup solve as `KuttaConvergenceError` for `JumpKutta` combinations too
   (deterministic failure, state restored).

### Verification (all run 2026-07-29, single-threaded, from `test/`)

```bash
# new suites
julia --project=.. -e 'using Test, LinearAlgebra; import FLOWPanel as pnl;
  using WriteVTK; include("test_helpers.jl");
  include("runtests_unit_kutta.jl")'            # 602/602 pass
julia --project=.. -e '... include("runtests_unit_kutta_routeb.jl")'  # 59/59 pass

# legacy regression: every unit suite in runtests.jl order
# (fmm, kernel_gradient, body, solver, liftingbody, wake, postprocess,
#  simulate, kutta, kutta_routeb, added_mass, warmstart, replay) — all pass
# example/analytical suites: pitching_wing_exp, pitching_wing,
#  pitching_wing_convergence, pitching_wing_pressure_comparison,
#  dji9443_trailing_edge, analytical — all pass.
# runtests_example_suddenly_started_wing.jl fails PRE-EXISTINGLY and
# unrelatedly: the worktree's uncommitted item-014 change to
# examples/suddenly_started_wing.jl added an even-n_span mesh-symmetry guard,
# and the example test still passes n_span=1 (throw at
# suddenly_started_wing_mesh, examples/suddenly_started_wing.jl:206 — no
# Phase 3 code in the stack). Left to the item-014 line of work.
```

Functional smoke (scratch, not committed): all four combinations converge on
the diamond fixture (pressure residuals to ~1e-13 of scale); warm-start round
trips for A/pressure, B/jump, B/pressure match the uninterrupted runs to
1e-10; configuration-mismatch resumes rejected.

### Known limitations / deferred

- `KrylovSolver`/`FGSSolver` support deferred (user amendment); rejection
  tested. The FGS trial hooks and the stray Krylov `@show` gating go with
  that follow-up.
- **RESOLVED 2026-07-29 by user decision: `TraceCorrected` is deprecated**
  (construction warns once; kept functional for the archived Gate A0 script,
  `test/formulation_test.jl`, and old-run replay; `GreenReconstruction` is
  the supported formulation). Rationale: Gate A0
  (`logs/dji_convergence_20260722/phase_02b_formulation_attribution.md`)
  rejected the affine-c channel as a lever (~0.74%), the `:green` estimator
  is identical to `GreenReconstruction` (Stage 6 identity), and the defect
  below made its velocity channel vacuous. Original flag retained for the
  record:
  during the FMM follow-up we established empirically that `influence!` never
  applies `wake_strength_shift` under ANY backend — even `DirectBackend`
  routes through the FastMultipole source buffers, whose `_induced_wake`
  variant has no shift (the shift lives only in the index-form `induced`,
  used by matrix assembly). `TraceCorrected`'s
  `backend_system isa DirectBackend` requirement
  (`FLOWPanel_formulation.jl:379-382`) therefore does not deliver the affine
  term into `body.velocity` either: its post-solve surface velocities (hence
  Bernoulli pressures/forces) appear to omit the −c attached-strip velocity
  under every backend. Not changed here (out of 015 scope);
  `_add_affine_attached_velocity!` is the ready-made fix.
- B2 multi-panel blocks are not constructible yet; the live-block seam
  (multi-row `live_rows` source-view arithmetic, single physical-step
  identifier, equal-strength internal-edge cancellation) is implemented and
  tested.
- Pre-existing repo gap, not introduced here: `runtests_unit_fgs_history.jl`
  is not registered in `test/runtests.jl`.

## Remediation after clear-context review (2026-07-29)

A three-agent clear-context review (core file vs architecture; integration
diffs vs legacy defaults; test-suite coverage vs claims) found one blocker and
several contract gaps. All were fixed in this session; both kutta suites and
the full legacy unit regression re-run green afterward.

### Blocker fixed — kernel-offset mismatch in the correction channel

The body-centroid affine add-on was evaluated at `kerneloffset_targets` while
the proportional γ = Cμ part of the SAME trailing-edge strips ran at
`kerneloffset_panel` (the self-pair conditioning rule in
`_self_panel_kerneloffset_conditioning`). On the test fixture both offsets are
equal, so no test could see it; on production rotor settings they differ by
~8 orders of magnitude exactly at the residual-defining TE panels, so the
converged closure physics would have been offset-inconsistent. Fixed:
`_kutta_reconstruct_body_velocity!` and the commit path now pass
`kerneloffset_panel` for body-centroid targets (wake probes are cross pairs
and correctly use `kerneloffset_targets`). New regression test
"reconstruction add-on uses kerneloffset_panel" runs the trial reconstruction
with deliberately distinct offsets and pins the add-on to the panel offset.

### Contract gaps fixed

- **§7.4 commit guard:** `_kutta_commit!` is now a guarded transaction — any
  exception during commit restores the complete pre-step snapshot (which also
  rolls back the provider) and throws `KuttaConvergenceError` with the
  previously never-constructed `:commit` cause. Tested with a provider whose
  `commit_pressure!` throws.
- **Degenerate-scale abort:** the `ArgumentError` from
  `_kutta_freeze_scales!` now restores the pre-step snapshot before escaping
  (previously the initial clean trial's mutation was left installed).
- **Newton–Krylov stats:** an unconverged GMRES direction (including one
  poisoned by a NaN FD product) is now rejected via `stats.solved` instead of
  being silently accepted.
- **Warm-start validation ordering:** `simulate_warmstart!` now runs
  `_validate_kutta_configuration` before any state mutation; an unsupported
  resume (e.g. `PanelParticleWake` + non-default config) throws an
  `ArgumentError` up front instead of a raw field error after mutation.
- **Stale-`W` guard:** Route A now rejects `set_Das_refresh=true` (the
  attachment operator and factorization are assembled once per run; the §8.1
  per-geometry-version refresh is not implemented — `geometry_version` exists
  but is only maintained on the Route B path).
- **Replay guard:** `replay` throws on a manifest containing a `kutta` table
  instead of silently replaying a non-default run as legacy.

### Test-suite hardening (review found the machinery could be a silent no-op)

- **Leverage test:** asserts the accepted correction is genuinely nonzero
  (calibrated ≈8.6e-3·C_s ≫ 1e-3·C_s floor) and that A/pressure moves the
  body solution ≈2.2% vs A/jump — a weak-leverage fixture can no longer pass
  everything with c ≡ 0.
- **Independent residual:** paired-centroid Δp is recomputed directly from the
  post-run committed velocity with the provider law and checked against the
  gate (removes the tautology of re-reading the driver's own gate fields);
  the A/jump baseline is confirmed to violate pressure continuity by ≫ the
  tolerance, so the metric has teeth.
- **Rollback at a settled step:** a provider failing only from the third step
  distinguishes "restored" from "never written"; a full
  `_kutta_restore_snapshot!` round-trip test asserts every snapshot field
  against scribbled-on nonzero state.
- **Rejection messages:** the deferral message and the other support-boundary
  rejections are matched by message substring, not just exception type;
  multi-body, wakeless, and Route A `set_Das_refresh` rejections added;
  rejection-before-mutation asserted on a settled body/wake pair.
- **`:jump` fallback:** upgraded from `isapprox(rtol=1e-12)` to exact `==`
  (verified genuinely bitwise vs the legacy trajectory).
- **Early-accept gate edge case (§11.4):** wide-open tolerances confirm the
  already-passing base point is accepted with c = 0 rather than tripping the
  Armijo/step-gate deadlock.

### Known limitations recorded (not fixed here, by scope decision)

- `_kutta_step!` still drops `diagnose_particle_influence` /
  `diagnostic_vertical`, and Route A always calls `update_TE!` (no
  `update_trailing_edges` pass-through).
- Trial reconstruction targets `(body,)` while commit targets
  `(body, wake_probes…)`; under an FMM `backend_system` the trees differ
  slightly (measured effect ≤ the 1e-6 backend-agreement tolerance; the
  independent-residual test bounds it at ≤ 2× the pressure gate on the
  fixture).
- Dead runtime fields `entry`/`r_prev` (and `geometry_version` on Route A);
  per-restart identity allocation; backtrack cap off-by-one
  (`max_backtracks + 1` evaluations); `KuttaDiagnostics` pinned to `Float64`;
  `include_final_filament` not checked by the finite-wake validation;
  diagnostics in the §11.4 early-accept branch report the proposed
  (never-taken) step.
- Route B live-row restore in `_kutta_warmstart_restore!` is immediately
  overwritten by the replay hook (values agree; the `haskey` validation is
  the operative part).

### Review riders routed elsewhere (NOT item 015)

The integration review confirmed the 015 hunks preserve legacy behavior
exactly (byte-identical code-motion; `live_rows[]==0` algebraic collapse;
`nothing`-gated metadata/warmstart), but flagged two co-resident item-014
changes in this worktree that DO alter legacy trajectories:
`set_Das_kinematic_arc=true` (new default; reroutes kinematic `Das` through
`accumulate_Das_arc!` on every rotating-body run) and the zero-Γ
`_shed_particles!` early return (changes `PanelParticleWake` particle
counts). The "legacy defaults remain unchanged" statements in this document
are scoped to the item-015 changes; the 014 riders need their own
attribution/regression under item 014.

## Phase 3 progress and approval log

- **2026-07-29 — Phase seeded.** Blocked pending explicit Phase 2 approval; no
  implementation has started.
- **2026-07-29 — Phase 2 approved with concrete-typing requirement.** Phase 3
  is ready to begin. No implementation has started.
- **2026-07-29 — User amendment: Backslash-only first implementation.**
  `KrylovSolver`/`FGSSolver` support deferred until the Backslash
  implementation is validated; both are rejected with an explicit message.
- **2026-07-29 — Implementation technically complete.** All four
  attachment/closure combinations implemented and mechanically tested
  (Backslash only); legacy defaults preserved bitwise on the default branch;
  new suites 602+59 assertions green; full unit-suite regression green.
  Awaiting explicit user approval; Phase 4 remains blocked.
- **2026-07-29 — FMM follow-up (user request): backend restriction removed
  and a latent reconstruction defect fixed.** The user asked for FMM
  compatibility, suggesting a hybrid where only the single-wake-row affine
  term is evaluated directly. Implemented as
  `_add_affine_attached_velocity!`: influence passes run correction-inactive
  through any backend; the affine −c strip velocity is added directly at all
  reconstruction targets (body centroids per trial; body + wake probes at
  commit). The `backend_system isa DirectBackend` rejection is deleted. **In
  the process the original channel's premise was refuted empirically:
  `influence!` never applies `wake_strength_shift` under any backend
  (DirectBackend also routes through the FastMultipole source buffers), so
  the pre-follow-up reconstruction had silently omitted the affine velocity
  term — the closure converged, but to a discrete condition where c acted on
  pressures only through μ. The add-on fixes this; accepted corrections
  change accordingly.** New tests: add-on vs production full-minus-suppressed
  unit-strip influence (1e-11); full-FMM (all backend slots) A/pressure and
  B/pressure runs agree with all-direct runs to 1e-6. Both kutta suites and
  the legacy wake/simulate/warmstart/replay suites green. `TraceCorrected`'s
  matching (likely vacuous) DirectBackend restriction is flagged under Known
  limitations, unchanged.
- **2026-07-29 — Clear-context review, remediation, and approval.** A
  three-agent independent review found one blocker (kernel-offset mismatch in
  the correction channel), several §7.4 deterministic-failure gaps, a
  warm-start validation-ordering gap, a replay gap, and two high-severity
  test gaps (no leverage assertion; tautological convergence assertions). All
  were remediated the same day (see "Remediation after clear-context review"
  above); the kutta suites (658 + 62 assertions) and the full legacy unit
  regression re-ran green. Per the user-approved plan for this session
  ("review Phase 3; remediate; if all green, treat the gate as satisfied"),
  Phase 3 is APPROVED and Phase 4 is unblocked. Two review riders belong to
  item 014 (`set_Das_kinematic_arc` default; zero-Γ shed skip) and are logged
  there; the legacy-preservation claims here are scoped to the 015 hunks.
