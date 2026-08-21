# Pressure-Continuity Kutta Formulation

**Opened:** 2026-07-29

**Current phase:** Phase 4 — testing, verification, and validation. V0 complete;
V1 was attempted 2026-07-29, blocked at the wing formulation-sanity gate,
reviewed 2026-07-30, and its gate rebuilt as a scripted experiment (see the
Phase 4 document).

**Current status (2026-07-30):** V1 is still BLOCKED at the wing gate, but the gap is now
attributed: the Dirichlet arm's sensitivity to the first-wake-row offset η is substantially
an artifact of the default `VelocityThroughSources` wake→body transfer, not physics —
switching to `GreenReconstruction` removes 93% of the Dirichlet-vs-Neumann differential and
all of its η-dependence. The residual settled-CL floor is **force reconstruction**, measured
over a 3-rung mesh ladder. Open residue: the **circulation half is not converged** (Green's
KJ gap still falling, 5.58→4.42% at 4440 cells, ~3.7× its own steady anchor), so the sub-5%
gate numbers are reported only — the gate record is untouched. Detail:
`015_pressure_continuity_kutta_condition/phase_04_testing_verification_validation.md`.

**Item-level approvals:** Technical [ ]; clear-context [ ]; user [ ]

## Objective and scope

Design, implement, and validate independently configurable near-wake
attachment and Kutta closure for one finite-wake, paired-edge, Dirichlet
`RigidWakeBody`. The target includes an opt-in nonlinear pressure-continuity
closure, a `Das`-free TE-anchored wake, and eventual support for `Backslash`,
`KrylovSolver`, and a fully constructed `FGSSolver`. Existing formulations and
defaults must remain unchanged.

Item 014 found that floor-free `Das` changes moved settled rotor-hover \(C_T\)
by about 37%, primarily through bound circulation. This item separates the
wake representation from the closure used to select circulation so their
effects can be attributed independently.

## Phase gates

| Phase | Deliverable | Status | Approval |
| --- | --- | --- | --- |
| [1](015_pressure_continuity_kutta_condition/phase_01_theory.md) | Theory and algorithm recommendation | **APPROVED — 2026-07-29** | [x] |
| [2](015_pressure_continuity_kutta_condition/phase_02_architecture.md) | High-level code architecture | **APPROVED — 2026-07-29** | [x] |
| [3](015_pressure_continuity_kutta_condition/phase_03_implementation.md) | Implementation | **APPROVED — 2026-07-29** | [x] |
| [4](015_pressure_continuity_kutta_condition/phase_04_testing_verification_validation.md) | Testing, verification, and validation | **IN PROGRESS** | [ ] |

Completing a phase does not authorize the next phase. Each phase document owns
its detailed progress and approval log. The three item-level approval columns
in `BRAINSTORM/INDEX.md` remain unchecked until all four phases are complete.

## Cross-phase invariants

- Wake attachment and Kutta closure are independent configuration axes. The
  attribution matrix is A/jump, A/pressure, B/jump, and B/pressure.
- Route A retains the rigid `TE → TE + Das` transition and is never described
  as `Das`-free. Route B owns a TE-anchored live panel and has no independent
  rigid `Das`.
- Jump closure enforces \(\gamma=C\mu\) in the same body solve. Pressure
  closure matches paired panel-centroid pressures.
- Existing public behavior and defaults remain unchanged.
- Old wake strengths, physical histories, and solver histories remain frozen
  during nonlinear trials and commit exactly once after convergence.
- Route B testing is ordered B1 persistent Euler panels with timestep
  refinement, then the B2 equal-strength cancellation gate and persistent
  \(N_s=2,4,8\) subpanels. Katz–Plotkin is a comparator; a higher-order panel
  kernel is last resort.
- Safeguarded global Broyden is primary; finite-difference Newton–Krylov is
  fallback. Steady Bernoulli is the first built-in pressure provider.
- Route B must demonstrate timestep/wake-resolution convergence; pressure
  continuity must demonstrate sharp-edge convergence.

## Cross-phase approval events

- **2026-07-29 — Phase 1 explicitly approved.** The complete theory,
  four-combination matrix, B1-then-B2 experiment order, nonlinear strategy,
  provider scope, and frozen-trial/one-time-commit contract are accepted.
- **2026-07-29 — Phase 2 opened.** Architecture is ready but not started. No
  implementation is authorized until Phase 2 is explicitly approved.
- **2026-07-29 — Phase 2 technically completed.** The canonical architecture
  now fixes the typed API, compatibility boundary, state ownership,
  trial/commit and failure contracts, nonlinear safeguards, diagnostics,
  metadata/restart behavior, and Phase 3 file/test seams. Phase 2 remains
  unapproved and Phase 3 remains blocked pending explicit user approval.
- **2026-07-29 — Phase 2 explicitly approved with an addition.** The user
  approved the architecture and added a binding requirement that every new
  instantiable public or internal struct be concretely typed. Phase 3 is ready
  but has not started.
- **2026-07-29 — Phase 3 scope amendment by user.** The first implementation
  supports only the `Backslash` body solver; `KrylovSolver`/`FGSSolver` are
  validation-rejected and their support is deferred until the Backslash
  implementation is validated.
- **2026-07-29 — Phase 3 technically completed.** All four attachment/closure
  combinations implemented (Backslash only) with the trial/commit transaction,
  Broyden + Newton–Krylov drivers, `SteadyBernoulliProvider`, Route B live
  block owned by `PanelWake`, diagnostics, metadata, and warm start; legacy
  defaults preserved bitwise on the default branch; 661 new mechanical test
  assertions plus the full unit-suite regression green. Detail and the exact
  verification commands are in the Phase 3 document. Phase 3 awaits explicit
  user approval; Phase 4 remains blocked.
- **2026-07-29 — FMM follow-up (user request).** The pressure closure is now
  backend-complete: reconstruction influence runs correction-inactive through
  any backend and the affine −c term (one wake row of TE strips) is added by
  an exact direct O(M×N) add-on. En route, the original reconstruction was
  found to have silently omitted the affine velocity term (`influence!` never
  applies `wake_strength_shift` under any backend) — fixed by the add-on;
  `TraceCorrected`'s matching DirectBackend restriction is flagged as likely
  vacuous in the Phase 3 document (separate item). Suites green.
- **2026-07-29 — Phase 3 reviewed, remediated, and approved.** A three-agent
  clear-context review found a blocker (the body-centroid affine add-on used
  `kerneloffset_targets` while the proportional strip part ran at
  `kerneloffset_panel` — invisible on the equal-offset fixture, ~8 orders
  apart in production), plus §7.4 failure-contract gaps, a warm-start
  validation-ordering gap, a replay gap, and two high-severity test gaps
  (nothing asserted the closure had leverage; convergence assertions were
  tautological). All remediated the same day; suites re-ran green (658 + 62
  kutta assertions + full legacy unit regression). Gate satisfied per the
  user-approved review-then-remediate plan. Detail in the Phase 3 document's
  remediation section; two riders routed to item 014. Phase 4 begins.
- **2026-07-29 — Phase 4 V0 complete; V1 attempted and blocked.** The full
  legacy unit regression plus the kutta suites (658 + 62) are green and the
  diamond four-combination record is complete, but the pre-registered
  Dirichlet-vs-Neumann wing sanity gate failed (14–26% CL gap) and the SSW and
  swept-wing matrices were correctly not run or interpreted.
- **2026-07-30 — V1 review.** A three-agent clear-context review found the V1
  gate diagnosis unreproducible (one of ~10 quoted numbers had a script behind
  it) and internally inconsistent (Γ agreeing ~3% while KJ CL differed 2.5×),
  the gate itself sampled at an unsettled \(t^*\le 2\) with no strength-level
  comparison, `examples/sweptwing.jl`'s default body silently flipped to
  Dirichlet with no toggle or formulation token, several diamond-table
  normalization blemishes, and doc/log gaps. The V1 localization claim is
  withdrawn, not merely softened.
- **2026-07-30 — V1 gate remediated and rebuilt.** Decisive new fact (user):
  the gate bodies are **open-tip**, so the interior-Dirichlet Green identity —
  which assumes a closed surface — was being applied outside its assumptions,
  and flat caps are reported to match Neumann-no-caps best. Flat centroid-fan
  tip caps are now the Dirichlet configuration (`add_flat_tip_caps`, shared by
  the SSW and swept-wing paths), the swept-wing example is back on its Neumann
  default behind a `SWEPTWING_BODYTYPE` toggle with a formulation token in every
  output path, the diamond metrics are single-normalization, and the gate is a
  scripted six-cell battery (three formulation/cap configurations × two
  `grad_mu` reconstructions) settled over \(t^*\in[7,8]\) with a pre-registered
  criterion. Detail, results, and commands in the Phase 4 document.
