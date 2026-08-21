# Smooth Surface-Vorticity Particle Shedding

**Opened:** 2026-07-29

**Current phase:** Phase 4 — all planned measurements complete 2026-08-04
(§1/§2 closed; §3 wing PASS, rotor T4 null PASS, rotor smooth-vs-legacy FAIL
on the unconfounded pair); two Ryan-approved attribution discriminators in
flight 2026-08-05 (larger-d/σ and equal-h); clear-context review and user
approval outstanding

**Item-level approvals:** Technical [ ]; clear-context [ ]; user [ ]

**Session history:**
[`016_surface_vorticity_particle_shedding/log.md`](016_surface_vorticity_particle_shedding/log.md)
(newest first) — see the logging provision at the end of this file before
recording anything.

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
| [3](016_surface_vorticity_particle_shedding/phase_03_implementation.md) | Implementation | **APPROVED 2026-08-03** | [x] |
| [4](016_surface_vorticity_particle_shedding/phase_04_testing_verification_validation.md) | Testing, verification, and validation | **MEASUREMENTS COMPLETE 2026-08-04** (all §1–§3 evidence in; awaiting review/approval) | [x] |

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

## Current status

Phase 3 was explicitly approved by the user on 2026-08-03, following the
clear-context re-review recommendation (R1/R2 guarded at construction time,
R3 closed on pre-registered Stage 4c evidence with `:upstream` retained on
evidence, R4 caveated on the Stage 4b data README; suites wake 583/583,
replay 125/125, warm start 26/26, simulate 81/81, full suite clean).

Phase 4 measurements are complete as of 2026-08-04. §1 mathematical
verification and §2 mechanical regression are closed (coverage map + five
gap-closure tests; wake 605/605, simulate 88/88, full suite clean). §3 case 1
(wing A/B, `examples/ssw_conversion_ab.jl`) **passed**: smooth ≡ legacy to
~0.1% CL/Γ/loading, ledger ≤2e−16, delta σ-linked. §3 case 2 (rotor hover):
T4 attribution null **passed** on both observables (ΔCT 0.008%, ε_Γ 0.028% —
attribution is a free parameter of measured size ~6e−6 CT), but the
smooth-vs-legacy equivalence **FAILED on the unconfounded same-code-state
pair** `p018_conv_legacy`/`p018_conv_smooth` (13036850/13036851): ΔC̄T +0.505%
(marginal vs ≤0.5%), ε_Γ max 1.080% vs ≤1% — reproducing the earlier
confounded numbers, so the ~1% one-signed inboard-peaked circulation delta is
a property of the conversion at rotor σ/chord (d/σ-clearance signature),
window-robust, with the 20-rev settling caveat unable to rescue it. The
legacy control reproduces `p018_b0` to −0.026%, clearing the code-state
question. Net: smooth conversion is validated as an exactly-conservative,
stable opt-in; aerodynamic equivalence is scoped to wing-like σ/chord and is
NOT established at rotor resolution (valid negative result per §3
pre-registration; larger-d/σ rotor test is the follow-up lever). Awaiting
clear-context review and user approval of Phase 4.

**In flight (2026-08-05):** the Ryan-approved larger-d/σ discriminator —
smooth arm `p018_conv_smooth_dasc0p82` (job **13051758**) vs legacy partner
`p018_dasc0p82` (job **13050926**), banners verified to differ only in
`conversion:`, doubling d/σ at unchanged σ. Readout: ε_Γ shrinks vs the
1.080% at Das=0.41c ⇒ clearance hypothesis confirmed (the FAIL is
operating-point-specific); persists ⇒ conversion-operator difference. Second
discriminator: equal-h smooth arm `p018_conv_smooth_h2p0` (job **13051763**,
`CONVERSION_OVERLAP=2.0` so smooth samples at legacy's h = σ/2.0 instead of
its σ/1.3 default; pairs with the completed `p018_conv_legacy`) — ε_Γ
collapses ⇒ the FAIL was a sampling-density artifact; unchanged ⇒ h inert.
The two discriminators are orthogonal (conv_overlap constant within the Das
pair; Das constant within the h pair). Also pending: `p018_dasc0p41` (job
**13050925**, 018's run) to fill the legacy 540/720 rows of the
particle-count table in the Phase 4 doc §2026-08-05.

**Next actions for a fresh agent (in order):**

1. Read `016_surface_vorticity_particle_shedding/log.md` FIRST (standing
   rule: it answers "what has been run"), then the Phase 4 doc's 2026-08-04
   (c) and 2026-08-05 sections.
2. Check the three jobs above (`sacct -j 13050926,13051758,13051763` via
   `ssh orc 'bash -lc "..."'`; ControlMaster socket must be live — Ryan runs
   `! ssh orc -fN` if not; MOTD can eat the first output line, guard with a
   leading `echo SENTINEL`). Ops details: 018's `ops_reference.md`.
3. Harvest each finished run to local `data/<case>/` exactly as done for the
   conv pair: `<case>_CT_vs_rev.csv`, `_CT_per_rev.csv`,
   `_case_metadata.toml`, `monitors/*.csv`, `_conversion_diagnostics.csv`
   (smooth arms only), plus the slurm log. Confirm metadata knob-diff before
   trusting any delta.
4. Reduce with matched windows revs 10–20:
   `python3 scripts/p018_analyze.py m1 <case> --revs 10 20` and
   `m2 <A> <B> --revs-a 10 20 --revs-b 10 20`. Pairings: 13051758 vs
   13050926 (Das=0.82c pair); 13051763 vs the already-harvested
   `p018_conv_legacy` (equal-h pair). Compare each ε_Γ against the baseline
   pair's 1.080% (in-band) / 1.344% (unclipped, r/R≈0.18). Γ̄ delta
   figure/CSV pattern: `data/p016_conversion_gamma/gamma_conv_pair.{png,csv}`
   (built with `p018_analyze.py`'s `gamma_profile`/`_interp`).
5. Record verdicts per the pre-registered readouts (Phase 4 doc §2026-08-05),
   update log.md (new dated entry, newest first), this Current status
   (rewrite, don't append), and INDEX row 016's outcome cell.
6. Then: fill the particle-table rows from `p018_dasc0p41` VTK
   (`NumberOfPoints` in the header of
   `data/p018_dasc0p41/p018_dasc0p41_wake1_particles/*.N.vtp`, N=540/719);
   run the still-outstanding clear-context Phase 4 review; and only when the
   discriminator story is complete, add the finished narrative to the lab
   notebook (Ryan's instruction — notebook waits for the complete story; a
   partial entry for 2026-08-04 already exists in
   `~/Dropbox/research/notebooks/journals/20260803.md`).
7. Settling caveat on every number: 20-rev runs, `M1 settled CHECK` (<15
   settled revs). Extension via matched restart segments is a Ryan decision.

Derivation, evidence, and the dated history of how this state was reached are in
the phase documents and in
[`016_surface_vorticity_particle_shedding/log.md`](016_surface_vorticity_particle_shedding/log.md).

## Logging provision

**Do not append session narrative to this file.** Dated session entries,
approval events, review findings, and measurement outcomes go in
[`016_surface_vorticity_particle_shedding/log.md`](016_surface_vorticity_particle_shedding/log.md),
newest first. Per-phase detail stays in that phase's document.

This file is edited **in place**, never appended to. Only four things change
here: the header block, the phase-gate statuses, the cross-phase invariants, and
the `## Current status` section above — which is rewritten to reflect the present
state rather than extended with new dated bullets.
