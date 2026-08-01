# Smooth Surface-Vorticity Particle Shedding

**Opened:** 2026-07-29

**Current phase:** Phase 3 — IMPLEMENTATION IN PROGRESS

**Item-level approvals:** Technical [ ]; clear-context [ ]; user [ ]

## Objective and scope

Design, implement, and validate an opt-in `PanelParticleWake` conversion
strategy that reconstructs the smooth surface vorticity of an outgoing panel
sheet and deposits area-weighted vortex particles. The existing conversion
from panel-edge strength jumps remains the exact default. Switching between
the two methods must require only one `PanelParticleWake` constructor option.

This item separates a physical wake-sheet perimeter from internal
panel-to-panel jumps. The smooth strategy retains the former, replaces the
latter by a finite-difference surface-gradient reconstruction, and assigns
each particle

$$
\boldsymbol\Gamma_p
=\boldsymbol\kappa(\boldsymbol x_p)\,\Delta A_p,
\qquad
\boldsymbol\kappa
=\boldsymbol n\times\nabla_s\mu
=-\boldsymbol n\times\nabla_s\hat\mu,
\qquad
\hat\mu=-\mu .
$$

Only staged planning documents are authorized by this item. No source change,
test execution, or aerodynamic simulation is part of the present phase.

## Phase gates

| Phase | Deliverable | Status | Approval |
| --- | --- | --- | --- |
| [1](016_surface_vorticity_particle_shedding/phase_01_theory.md) | Theory and algorithm recommendation | **APPROVED** | [x] |
| [2](016_surface_vorticity_particle_shedding/phase_02_architecture.md) | High-level code architecture | **APPROVED** (as amended in its §13) | [x] |
| [3](016_surface_vorticity_particle_shedding/phase_03_implementation.md) | Implementation | **IN PROGRESS** | [ ] |
| [4](016_surface_vorticity_particle_shedding/phase_04_testing_verification_validation.md) | Testing, verification, and validation | **BLOCKED ON EXPLICIT PHASE 3 APPROVAL** | [ ] |

Completing a phase does not authorize the next phase. Each phase document owns
its detailed progress and approval log. The three item-level approval columns
in `BRAINSTORM/INDEX.md` remain unchecked until all four phases are complete.

## Cross-phase invariants

- `PanelParticleWake` has a concretely typed panel-to-particle conversion
  strategy axis. Legacy edge-jump conversion is the unchanged default;
  surface-vorticity conversion is opt-in through one constructor option.
- For the legacy strategy, particle locations, strengths, counts,
  `method_trailing`/`method_unsteady` behavior, metadata, replay, and warm
  starts remain unchanged.
- Physical doublet strength $\mu$ and stored wake strength
  $\hat\mu=-\mu$ remain distinct. The smooth strategy uses
  $\boldsymbol\kappa=\boldsymbol n\times\nabla_s\mu
  =-\boldsymbol n\times\nabla_s\hat\mu$ with package signs. It never creates
  an outside-zero strength solely to complete a stencil.
- `nwakerows=N` continues to mean \(N\) active aerodynamic panel rows.
  The outgoing row is converted as if it had been shifted into a ghost panel
  \(N+1\), *after* the newly shed row is formed. **Amended 2026-08-01
  (Phase 2 §13.1):** that shift is an indexing convention, not storage.
  Conversion runs before the row shift, so the ghost is the pre-shift final
  active row read in place; `PanelWake` storage is unchanged and no row is
  hidden from any source, probe, filament, or VTK view. Only the not-yet-shed
  new row 1 is staged, in the conversion workspace, and only when \(N=1\).
- The ghost uses its immediately upstream active panel for a first-order
  streamwise difference, including when `nwakerows == 1`. Spanwise
  reconstruction uses centered interior and one-sided root/tip differences.
  No `Das`, body-strength, virtual-TE, or outside-zero sample is admitted.
- True wake-sheet perimeter circulation is retained and internal panel
  boundaries are not replicated. The selected smooth-policy handoff cancels
  the active-panel downstream edge and deposits the complete reconstructed
  surface vorticity plus true root/tip closure. The retained-filament
  alternative is rejected because it double-counts the affine streamwise
  component or requires an incomplete area field.
- Surface quadrature targets spacing
  \(h=\sigma/\mathrm{overlap}\), uses at least one subdivision in both panel
  directions, maps subcell centroids through the bilinear panel, and uses
  physical-area weights.
- Every deposited particle satisfies
  \(\boldsymbol\Gamma_p=\boldsymbol\kappa(\boldsymbol x_p)\Delta A_p\);
  each panel's particles sum to its chosen quadrature of reconstructed
  surface vorticity.
- Particle-capacity failure is transactional: no fraction of an outgoing
  panel row is deposited.
- Mathematical verification, mechanical regression, and aerodynamic
  validation remain separate evidence tracks.
- Newly deposited root/tip particles may lie close to the uncancelled
  root/tip vortex filament of the retained active panel wake. The resulting
  panel-induced velocity gradients must be diagnosed separately from
  particle-particle alignment effects.

## Cross-phase approval events

- **2026-08-01 — Phase 2 APPROVED; Phase 3 unblocked and started.**
  Clear-context review checked the architecture against the live source and
  recorded one simplification plus two corrections as Phase 2 §13: the
  persistent ghost row and its ~20-site visibility exclusion are unnecessary
  (conversion runs before the row shift, so the ghost is the pre-shift final
  active row read in place — Option B selected, Option A recorded verbatim);
  wake serializers live in `src/FLOWPanel_replay.jl`, not
  `src/FLOWPanel_metadata.jl`; and the existing unknown-tag `NoShed()` fallback
  in `_deserialize_wake_shedding` must not be reused for conversion decoding.
  User explicitly approved. Phase 3 began with a golden-reference
  characterization test pinning the legacy conversion bit-exactly before any
  refactor (`test/data/legacy_wake_conversion_reference.jl`,
  `make_conversion_fixture` in `test/test_helpers.jl`). Phase 4 remains blocked.

- **2026-07-29 — Phase 2 technically complete; awaiting approval.** The
  architecture fixes the exact legacy-default/smooth-opt-in API, explicit
  active capacity and invisible ghost storage, rank-aware reconstruction,
  physical-area deposition, true root/tip closure, cancelled-edge handoff,
  transactional capacity/failure semantics, concrete diagnostics, and
  metadata/replay/warm-start contracts. The Phase 2 approval checkbox remains
  unchecked and Phase 3 remains blocked pending explicit user approval.

- **2026-07-29 — Phase 1 APPROVED; Phase 2 unblocked.** Clear-context review
  independently re-derived the handoff/row ledgers (7)–(11) and verified all
  seven code-anchored claims (sign convention, vertex order, final-filament
  cancellation, legacy handoff deposit, N+1/N+1 storage, zero-Γ guard, body
  surface-vorticity orientation) at exact source lines; the Alternative B
  rejection is forced by equation (10). User explicitly approved. Non-blocking
  polish notes recorded in the Phase 1 log (sign-chain code anchors, commit
  the uncommitted zero-Γ guard before Phase 3, duplicate (3)/(3a) label).
  Phase 2 may begin; Phases 3–4 remain blocked on their upstream approvals.

- **2026-07-29 — Phase 1 handoff ledger completed; technically complete.**
  The live package orientation gives
  $\boldsymbol H=(\hat\mu_A-\hat\mu_G)
  (\boldsymbol r_R-\boldsymbol r_L)$. The planar-affine streamwise area
  integral equals this vector, internal spanwise edges telescope, and
  repeated conversions leave no artificial row-boundary filament.
  Cancelled-edge Alternative A is selected; retained-filament Alternative B
  is rejected. Phase 1 approval remains unchecked and Phase 2 remains
  blocked pending explicit user approval.
- **2026-07-29 — Item opened and Phase 1 technically completed.** The sign
  convention, finite-difference reconstruction, rank-deficiency policy,
  perimeter/handoff ledger, surface quadrature, limitations, and acceptance
  scenarios are recorded. Phase 2 remains blocked pending explicit user
  approval of Phase 1.
- **2026-07-29 — Phase 1 reopened after user review.** `PanelWake` will be
  planned with an always-allocated non-source ghost panel so
  `nwakerows == 1` has a genuine two-panel one-sided stencil. The ghost is the
  outgoing panel being converted, not a previously converted passive marker.
  Cancelled-edge and retained-filament handoff formulations must now be
  derived and compared before Phase 1 can be technically complete.
