# Brainstorm Items for Consideration

This file is a **routing table**, not a record. One row per `.md` file in this directory,
giving its subject, approval state, and headline outcome in one line. All detail — data,
derivations, run commands, and the log of what was tried and refuted — lives in the item's
own `.md` file (and its subdirectory, where one exists).

## Mission Statement

These items should lead to validation of FLOWPanel for the rotor hover $C_T$ case.

## Agent Workflow

Agents should begin by reading this file to understand how to navigate the brainstorm items
and their approval process. When an agent is assigned to work on a specific item, they should
read that item's `.md` file for detailed context, proposed path, acceptance criteria, and any
caveats. The agent then carries out the work described in the item, **keeping detailed notes
and data in the item's own `.md` file** — not here.

### RESET BRIEF convention (required for new items)

Every item file **created from now on** must carry a `## RESET BRIEF` section near the top:
≤40 lines, replace-never-append, stating current status, live jobs, standing rulings, and
next actions — everything a fresh agent needs to resume without reading the rest of the file.
Existing items are not required to retrofit one. To catch up on an item, spawn the
`brainstorm-scout` subagent (or read only the RESET BRIEF where one exists) rather than
reading full item files into the main context.

## Progress Summary

Use the approval columns as gates for each brainstorm item:

- **Technical Completion**: Check this only after the investigation or implementation work
  described by the item is complete enough to evaluate against its own acceptance criteria.
- **Clear-Context Approval**: Check this only when an agent with clear context can approve
  the item immediately as satisfying its requirements. Conditional approval is not allowed:
  if approval depends on further fixes, clarification, or follow-up work, leave this
  unchecked and treat the item as still being in the technical completion phase. Approval
  should seek ways to improve the item, considering the following in order: 1) relevance to
  the mission statement; 2) correctness of implementation and conclusions; 3) clarity and
  conciseness of presentation.
- **User Approval After AI Discussion**: Check this only after the user has discussed the
  item with another clear-context AI agent and explicitly approved the item or its conclusion.

### Rules for the Outcome column — keep this file clean

- **Format:** headline result, plus at most one "next" clause. **Hard cap ≈35 words / two
  displayed lines.** State what is *true now*, not how it was learned.
- **Never** put in this file: refuted hypotheses, superseded numbers, dated sub-findings,
  data tables, derivations, run commands, file-or-line-level detail, process notes, or
  "also / additionally" clauses. Job IDs only for the single in-flight set, if any.
- **The log of failures and refuted mechanisms is valuable — and it belongs in the item's own
  `.md` file**, in a dated section. Deleting it from here is not deleting it.
- **When a result is superseded, replace the cell. Never append to it.** The history is
  appended to the item file instead.
- For an open item, the "next" clause should name the item file's own entry-point heading
  (e.g. `→ item "Current status"`), not restate the plan.
- Each item file should carry a `## Current status` (or equivalent) section that its INDEX
  cell can point at. **This file must never be the sole record of any fact.** Before
  shortening a cell, confirm its content exists in the item file.
- *If a cell no longer fits on two lines, the excess belongs in the item file. Rewrite,
  don't extend.*

| File | Subject | Technical Completion | Clear-Context Approval | User Approval After AI Discussion | Outcome |
| --- | --- | --- | --- | --- | --- |
| `001_rotor_hover_low_effort_checks.md` | Low-effort implementation checks for the rotor-hover CT discrepancy: backend agreement, force-method agreement, Kutta-condition sensitivity, and semi-infinite rigid-wake smoke test. | [x] | [ ] | [ ] | Backend and Kutta cleared; the 16.3% force-method spread is real, but wake modeling dominates (semi-infinite wake doubles CT). Routed to 002/004. |
| `002_rotor_hover_prescribed_helical_wake.md` | Diagnostic finite panel helical-wake idea with TE-strip strengths and pitch iteration from azimuthally averaged inflow. | [ ] | [ ] | [ ] | Open. Seeded short wakes give plausible CT but the geometry will not converge — the capped Picard update is the blocker. Next: residual-minimizing wake-shape solver. |
| `003_rotor_hover_ccblade_bem_comparison.md` | CCBlade/BEM cross-check plan to separate missing wake physics from geometry, airfoil, Reynolds-number, or operating-condition assumptions. | [x] | [ ] | [ ] | BEM on identical geometry/airfoil/Re brackets experiment (CT 0.060–0.071 vs 0.072) and stays 18–41% above panel ⇒ shortfall routes to wake/induced-velocity modeling. |
| `004_rotor_hover_wake_decay_audit.md` | Audit plan for wake-bearing runs that pass through `CT ≈ 0.072` and later drift lower, focused on numerical wake-modeling mechanisms and transient behavior. | [ ] | [ ] | [ ] | Not started — audit plan only. |
| `005_rotor_hover_stable_converged_particle_wake.md` | Constructive plan targeting **stable, non-oscillating** particle-wake thrust (a flat history at any value is a win — magnitude is 002/003/004's job) via gradual RPM spin-up, a ramped-then-withdrawn convecting freestream, and downstream particle truncation. | [x] | [x] | [x] | APPROVED. Staged startup kills the impulsive startup ring; the residual is a bounded 2/rev limit cycle at truncation depth ≥4R, and no numerical damper reduces it. |
| `006_rotor_hover_converge_stable_near_reference.md` | Convergence track: drive the particle-wake CT history to be **both** stable **and** near the 0.068–0.072 reference at validated loose truncation (≥4R). Unlike 005 (any flat value), 006 also requires magnitude. | [ ] | [ ] | [ ] | PARTIAL — magnitude met (η=1.0 → CT 0.0713, in band), stability not (2.23% ripple vs 0.5% tol). Mechanism unexplained ⇒ opened 014. |
| `007_pressure_monitor_reliability.md` | Which surface-pressure monitor to trust for rotor-hover cycle-mean CT: code-verified theory analysis plus numerical tests T1–T4. | [x] | [ ] | [ ] | Steady Bernoulli is the most reliable; defensible 10-rev mean CT 0.0685 ± 0.0013 vs experiment 0.072. Laplace Du/Dt is inflated; unsteady Bernoulli is broken for particle wakes. |
| `008_particle_overlap_vorticity_evolution.md` | First pose the overlap-aware particle-strength linear system for generic kernel `ζ`, then specialize it to FLOWVPM/FLOWPanel particle-wake updates. | [ ] | [ ] | [ ] | GO — the live single-particle update sits ~62% from the overlap-projected RHS, so an overlap solve would materially change it. Next: full-RHS/SFS and Pedrizzetti comparison. |
| `009_particle_field_divergence_free.md` | Theory + ForwardDiff numerical check of whether the vortex-particle-induced velocity field is divergence-free, at particle centers and random points for heavy/light/no overlap. | [x] | [x] | [x] | APPROVED. The particle-induced velocity field is divergence-free to round-off at any regularization or overlap (theory + ForwardDiff). |
| `010_particle_vorticity_curl_vs_basis.md` | Theory + single-particle ForwardDiff check comparing particle vorticity from `curl(velocity_gradient)` against direct `ζ`-basis evaluation, with root-cause analysis and self-consistency improvement options. | [x] | [ ] | [ ] | curl-of-J is the projected vorticity with a 2/3 self-term, not a `ζ` artifact. M6 negative — the block overlap operator is ~20× more ill-conditioned ⇒ lever is overlap reduction (012). |
| `011_pedrizzetti_overlap_residual_efficacy.md` | Quantify how much of the 008/M2/M6 overlap-aware strength correction the cheap, already-shipped Pedrizzetti relaxation captures on its own — gating whether an explicit CG/block overlap solve is worth its cost. | [x] | [ ] | [ ] | Pedrizzetti fixes direction only (median 1.1°); the leftover ~98%-magnitude overlap correction is structurally out of its reach ⇒ routes to 012. |
| `012_robust_resolution_overlap_management.md` | Explore smarter, robust particle-resolution management that bounds the overlap ratio `σ/h` while preserving vortical structure — potentially dissolving the conditioning, divergence-discrepancy, and Pedrizzetti-workload problems at their shared source. | [ ] | [ ] | [ ] | Staged, not started. Bound `σ/h` (≈4.2 measured) via overlap-bounded moment-conserving merging. Next: characterize `σ/h` over wake and time. |
| `013_unsteady_wake_validation_cases.md` | Web-search survey of published validation cases for an unsteady inviscid panel method that would show off a shed free wake (VPM or doublet-panel), prioritizing cases with both existing panel-code AND experimental data. | [x] | [ ] | [ ] | Eleven ranked free-wake validation cases, each with experiment plus panel/free-wake data; high-Re picks 11→4→7→6. Documentation only, no implementation. |
| `014_first_wake_row_offset_selection.md` | **OPEN QUESTION:** on what principle should the first-wake-row offset `Das`/η be chosen, so that CT is a converged prediction rather than a tuned number? | [ ] | [ ] | [ ] | `Das` is a physical length on a log plateau (+0.2%/doubling), not η — but the rotor's +37% has no surviving explanation; the sheet/particle cause was spun off to 017/018. |
| `015_pressure_continuity_kutta_condition.md` | Gated theory, architecture, implementation, and validation of independently configurable wake attachment and Kutta closure for paired-edge Dirichlet lifting bodies. | [ ] | [ ] | [ ] | Phase 4 V1 BLOCKED at the wing gate. The Dirichlet η-sensitivity is largely a `VelocityThroughSources` wake-transfer artifact (Green removes 93%); the circulation half is still unconverged. |
| `016_surface_vorticity_particle_shedding.md` | Gated design of an opt-in smooth surface-vorticity conversion for `PanelParticleWake`, retaining physical root/tip closure and the unchanged legacy edge-jump default. | [ ] | [ ] | [ ] | Phase 3 approved. The smooth conversion is exactly conservative and stable, but rotor aerodynamic equivalence FAILS (ε_Γ 1.08%) — scoped to wing-like σ/chord. Next: two in-flight discriminators → item "Current status". |
| `017_sheet_particle_representation_equivalence.md` | Wing sheet-to-particle equivalence map and mechanism attribution under controlled σ, spacing, and handoff distance; rotor transfer deferred. | [x] | [ ] | [ ] | Fixed-overlap legacy shedding is admissible in 14/15 cells (σ/c≤0.6 over d/c=0.25–4, or d/c≥1 at σ/c≈1.2); the sole failure is near-body chord-scale cores. Phase D deferred. |
| `018_dji9443_hover_convergence_campaign.md` | Publishable hover convergence campaign: converge physical `Das`, timestep, particle σ, shedding distance, truncation, mesh span+chord, Green Δ, and merging null into one converged CT + Γ(r/R) with a closed error budget. | [ ] | [ ] | [ ] | LIVE. etadas A/B NULL in CT (Das law redistributes Γ̄ only) ⇒ carrier offset is an interaction (etadas_n4: the N knob is −2.4%, dominant); Phase 16 OPENED 2026-08-14: chord–σ co-scaling Das law (σ_j=s*·c_j, Das_j=λσ_j satisfies all 3 admissibility bands at once) → `phase_16_chord_sigma_coscaling.md`; dt ladder pins at N=1 uniform-D; drift reframed as burst-sampled equilibration → item "Current status and next actions". |
| `019_sigma_selection_procedure.md` | Formalize the general σ-selection procedure (initialize from σ_eq/σ_stab, probe the ΔtZ margin via tripwire, iterate) and demonstrate it reproducibly on rotor hover, with a regime-map parameter sweep showcasing the blow-up. Publishable quality; reuses 018 runs, gap runs only. | [x] | [ ] | [ ] | COMPLETE 2026-08-12 — σ_stab validated at the screen ε-crossing (~2%); k=0 σ₀ FAILED certification (drift, not per-step margin) ⇒ σ*=0.0381R via L1; P4b: in-band failure is tail-driven ignition, bulk static; drift-arrest v1.1 amendment approved into P1; notebook entry written (journals 20260812). → item dated sections. |
| `020_sigma_aware_subgrid_closure.md` | Physically grounded σ-aware (viscous, σ-channel) subgrid closure for FLOWVPM: theory → gap demonstration → implementation → validation, gated by Ryan discussion at each phase. | [ ] | [ ] | [ ] | Phase 2 complete at revised gate: corrected map crosses Euler death healthy, then tail-runaway fails at step 243; both registered Boolean arms fail/not-confirmed. Closure remains unvalidated candidate; STOP before predictive fixtures/Phase 3. → item "Current status". |
| `021_rotor_hover_solver_benchmarks.md` | Publishable solver benchmark campaign on rotor hover (Backslash coupled/iterative, Krylov ± preconditioner, FGS, FGMRES+FGS): setup vs per-timestep cost, warmstarting, multi-body blown wing, over a panel-count ladder. | [ ] | [ ] | [ ] | Phase 0 clear-context approved (W1–W6; ILU smoke winner, 290→15 iters, ~19× vs plain GMRES). Phase 1 (consistency/calibration) executing → item "Current status". |
| `022_rotor_hover_ground_effect.md` | Rotor hover in ground effect at h/R=1 via existing paneled-ground machinery (FlatGround + FlatGroundSolver + Gauss–Seidel multi-body): CT with vs without ground at matched settings, ground-representation convergence, below-ground particle policy. | [ ] | [ ] | [ ] | Fine-rung headline IGE/OGE = 1.061 ± 0.038 (anchor 1.067, experiment 1.078); coarse IGE blew up ⇒ Phase-3 policy arms (damping + floor-at-ground, 13246557/58) in flight → item "Current status and next actions". |
| `023_018_runtime_cost_profiling.md` | Profile a mature-wake 018 production run (per-phase timing + profiler at 36,752 panels / ~181k particles) and rank quantified cost-reduction levers; recommendations only, implementation routed to Ryan/021. | [ ] | [ ] | [ ] | Opened 2026-08-20; Phase A (wall_s mining) executing → item "Current status". |
| `024_nwakerows_zero_convert_at_shed.md` | Add `nwakerows=0` convert-at-shed mode to `PanelParticleWake` (particles placed at the Das line, no surviving free sheet), bit-identical for N≥1; rigid Das/Kutta row untouched. | [ ] | [ ] | [ ] | Opened 2026-08-20 (Ryan; promoted from 018's parked cleanup, expected result-neutral). Implementation delegated; local-only, cluster deploy Ryan-gated → item Log. |
| `025_small_das_free_row_unsteady_wake.md` | Restore fully unsteady wake modeling: rigid Das shrunk to a measured Kutta-only floor, N free wake-panel rows carrying the near wake time-accurately to the 3–4σ particle handoff; stability-screened against the historical free-row blow-ups. | [ ] | [ ] | [ ] | STAGED 2026-08-20 (Ryan: stage, don't implement). Phase 1 = Kutta-floor mini-sweep (0.25c is NOT a validated floor — clearance confound, see 018 ledger amendment); Phase 2 = free-row stability envelope → item plan. |
| `024_farfield_approximation_direct_kernels.md` | Per-pair far-field approximations inside the direct kernels (singular 1/r for particles, point source/doublet for panels) to cut near-field cost; HARD GATE = eligibility census on the mature 023 state before any implementation. | [ ] | [ ] | [ ] | Staged 2026-08-20 (Ryan directive); Phase 0 census defined, not run → item "Current status". |
| `025_kernel_regularization_update.md` | Selectable panel-kernel regularization with a new default chosen by peak-velocity/gradient comparison (compact-support vs Gaussian; Vatistas retained as legacy); theory + implementation/testing phases; 018 driver-compatibility task. | [ ] | [ ] | [ ] | Staged 2026-08-20 (Ryan directive); Phase 0 comparison not run; nothing implemented → item "Current status". |
