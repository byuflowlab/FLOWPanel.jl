# 020 — σ-aware subgrid closure for FLOWVPM (theory → gap demonstration → implementation → validation)

**Opened:** 2026-08-06 (Ryan, out of the 019 remedies discussion). **Status:** Phase 2 COMPLETE AT REVISED GATE — mixed/negative preregistered verdict; STOP before predictive-closure work.

**Item-level approvals:** Technical [ ]; clear-context [ ]; user [ ]

## Current status

**2026-08-12 — Phase-2R decisive screen COMPLETE.** Corrected
frozen-gradient geometric integrator passes FLOWVPM local 28/28, SFS 505/505,
merging 36/36, filament calibration 40/40, FLOWPanel simulate coverage, and
cluster Julia 1.11; the final focused suite is 28/28 after adding a prescribed
time-varying-gradient convergence test. Checksums match across both cluster checkouts. Job
**13154223**, tag `scr_p020r_geom_s020v`, is the registered $\sigma/R=0.02$
rerun; startup banner verified (`expint:true`, viscous, NT=36, D=3.4).
It crossed the Euler control's step-237 failure healthy ($u_{max}=36$ m/s,
$M=0.173$, median $\sigma/\sigma_{shed}=1.037$), then failed at step 243 after
an abrupt tail runaway ($u_{max}=19.6$ km/s, $M=221.8$, min ratio 0.0137).
The explicit non-finite-ratio guard stopped the next step; peak RSS was only
20.6/64 GB and no negative $\sigma$ occurred.

**Registered rescore:** “threshold moves” FAILS because the corrected run did
not survive 8 revs; “fidelity collapse remains” is NOT CONFIRMED by its literal
Boolean because median stayed above 0.85 and $M>0.2$ was not sustained for ten
consecutive recorded steps. The qualitative tail runaway is real but is not
substituted for the preregistration. Conclusion: the corrected integrator
removes a local numerical death mode but is insufficient at field level; the
proposed closure remains an unvalidated candidate, not a demonstrated need.
Phase 2 is complete at the revised gate; explicit approval is required before
the independent predictive fixtures or Phase 3.

**2026-08-12 — independent theory verification REOPENED Phase 2.** The landed
`euler_exp` froze $S$ rather than the velocity gradient and therefore changed
physical stretching (aligned-strain multiplier tends to $5/3$ instead of
$e^{2\Delta tZ}$). The Stage-A/B mechanism verdict is withdrawn. Binding
corrections also scope the diffusion closure to an isolated/localized blob,
withdraw the unconditional stability and angular-impulse claims, and relabel
the 1.9–3.3-decade gap as conditional. Execution entry point:
`phase_02r_remediation_plan.md`. Leg 1/Leg 3 remain supporting evidence;
Phase-2R requires analytic verification plus one corrected $\sigma/R=0.02$
screen before the gate can be reconsidered.

**2026-08-07 — Phase 2 (gap demonstration) COMPLETE; evidence pack at the
gate.** All three legs measured; results of record + composition in
`020_sigma_aware_subgrid_closure/phase_02_evidence_pack.md`; four figures
built (fig_leg1/fig_stage_a/fig_stage_b/fig_leg3, pdflatex-verified).

- **Leg 1:** ceiling stability-set (p-arm PASS p=3.0; Γ_implied band FAIL
  as registered, boundary-localized — recorded); **gap to physics
  1.9–3.3 decades**.
- **Leg 2 Stage A:** P-2 CONFIRMED (zero spurious events under the exact
  map; fidelity grind to the floor persists).
- **Leg 2 Stage B (13065481, COMPLETED 8:10 h):** BOTH registered arms
  fire — zero negative σ ever (min min_sr +0.0396 = the CoreSpreading
  floor to 3 digits), outlives the Euler twin's step-237 death by
  ~1.3 revs, then the field-coupled Γ-runaway destroys the wake by
  trim-out (n 268k→15k, u>1000 m/s @294). **Closure needed for stability
  AND physics.**
- **Leg 3 (13065482 vs s025v):** SFS viscous-null CONFIRMED (onsets 269
  vs 286 = 5.9% ≤ ±15%; interleaved envelopes; same Γ-route death).

**Composition (§7.4 claim carried):** the subfilter viscous σ-channel
model is the missing piece; the exact integrator is its necessary
companion, not its substitute. **STOPPED at the Phase-2 gate — Ryan
review required before any Phase-3 work.**

### Rulings needed from Ryan at the Phase-2 gate (written for a fresh agent, 2026-08-12)

A fresh agent must NOT start Phase 3 until Ryan rules on these three
points. Everything referenced is self-contained in
`phase_02_evidence_pack.md` (results) and `phase_01_theory.md` (§4.6
constants, §7 registration, §8 Phase-3 tests).

1. **Accept the two criterion-as-written failures as recorded
   deviations?** (a) Leg-1's Γ_implied band criterion — registered as
   "$\Gamma_{\mathrm{implied}} \in [0.5,2]\Gamma_v$ for ≥3 consecutive
   rungs", measured 1/6 viscous rungs in band, because the
   filament-ceiling signature localizes AT the ignition boundary
   (igniters 0.6–3.2, survivors 4–7× below) instead of spanning the
   ladder. The pack's stability-set verdict rests on the p-arm (p = 3.0
   ≥ 1.5 PASS) + boundary localization + 019's ε-crossing at σ_stab.
   (b) Stage-B's "settled median σ/σ_shed" readout — unmeasurable
   because the wake self-destructed by trim-out before any settled
   window; the M > ε arm (sustained from step 225) fired the fidelity
   criterion instead. Options: accept both as recorded (pack stands as
   written), or direct a re-scoring/amendment.
2. **Closure constant κ for Phase 3.** The fixed-filter closure
   $\nu_t = \nu + (\kappa\sigma_0)^2|S|$ needs its one constant fixed
   *before implementation* (the item's no-retuning rule). Proposed
   default $\kappa = 1/\sqrt5 \approx 0.447$: strained cores then
   equilibrate at exactly $\sigma_0$ (maintained scale = shed scale;
   numerically verified in the Phase-1 scratchpad). Larger κ = fatter
   maintained cores (more margin, more diffusion); smaller = closer to
   the laminar behavior. Related sub-choice, also open from the Phase-1
   gate list: the $|S|$ estimator — full strain-rate magnitude from the
   velocity gradient (candidate default) vs the already-computed axial
   projection $s = 5Z$ (cheaper, underestimates off-axis strain).
3. **Phase-3 scope and go/no-go.** If the pack is accepted: implement
   the closure per phase_01 §8 (tests T1–T8 at pre-registered
   tolerances; default-off/bit-identical/test-gated; `ViscousScheme`
   subtype with per-particle ν_t is the natural home —
   `FLOWVPM_viscous.jl`, see phase_01 §3d). Sub-questions: (i) does the
   already-landed `euler_exp` integrator (default-off, tested, deployed)
   ship as a Phase-3 component (theory says closure+integrator are
   complements — the pack measured that closure-less integration is not
   sufficient, and Phase-1 §4.5 shows closure-only still Euler-flips
   below ~1.1σ_stab); (ii) does Ryan want the T7 ignition-fixture bound
   and T2 fixed-point tolerance checked against the *composed*
   (integrator+closure) scheme's own discrete fixed point (Phase-1
   numerical finding: splitting parks σ* ~45% above the continuous
   value at ΔtZ ~ 2 — benign, but the test must assert the right
   number).

---

**2026-08-06 (late) — Phase 2 (gap demonstration) LIVE.** Ryan passed the
Phase-1 gate ("work on phase 2") and authorized HPC submission + the Leg-3
fallback pair. Deliverable in progress:
`020_sigma_aware_subgrid_closure/phase_02_evidence_pack.md` (results of
record; thresholds quoted from phase_01 §7).

- **Leg 1 (strain ceiling): COMPLETE.** Campaign E NT=36 grid fully
  harvested (survivors s030v/sstab/s038/s038v; igniters s015/s015v/s020v/
  s025/s025v). p-arm PASS (survivor p = 3.01); Γ_implied-band arm **FAIL
  as registered** (1/6 in band) — the filament signature localizes at the
  boundary (igniters 0.6–3.2, survivors 4–7× below Γ_v); stability-set
  conclusion SUPPORTED with that on record. **Gap = 1.9–3.3 decades**
  (Z_res ≈ 957 s⁻¹ vs Z_phys at r_c = 0.05c–0.01c). Figure fig_leg1 built.
- **Leg 2 Stage A: COMPLETE (recorded deviation** — no per-particle
  trajectories retained; census + trace-re-integration reinterpretation).
  P-2 CONFIRMED: corpse field carries 908 σ-flip / 1545 Γ⊥-flip
  exceedances (healthy controls: zero); exact map: zero spurious events,
  σ still grinds to ~the floor (argmin trace 1.08× ≤ 1.1 criterion PASS;
  ignition-core 1.28× explained — Euler's near-zeroing at ΔtZ≈1 IS the
  artifact). Figure fig_stage_a built.
- **Leg 2 Stage B: IN FLIGHT** — `euler_exp` exponential integrator landed
  (default-off/bit-identical/test-gated; FLOWVPM 26/26 + FLOWPanel 4/4,
  local AND cluster; deployed md5-verified). Job **13065481**
  (`scr_p020_exp_s020v` = ignited s020v clone + WAKE_EXPINT=true).
- **Leg 3: IN FLIGHT** — no SFS-off run existed anywhere ⇒ registered
  fallback: job **13065482** (`scr_p020_s025v_nosfs`) vs the ignited
  `scr_p019_s025v` as on-member. Threshold: time-to-ignition ±15% + min-σ
  overlay.

Next: harvest both (~4–7 h wall; banner check + terminal-state monitor
armed), score against §7.2/§7.3 verbatim, complete the pack's Composition,
then STOP at the Phase-2 gate (Ryan reviews before any Phase-3 work).

**2026-08-06 — Phase 1 (theory) COMPLETE DRAFT, awaiting the gate discussion
with Ryan.** Deliverable: `020_sigma_aware_subgrid_closure/phase_01_theory.md`
— standalone-readable; every §4 map claim numerically verified (19/19 checks).
Headline results:

1. **The derivation produces the closure** (§2.3): retaining the
   $E_{\mathrm{adv}}$ term the original rVPM treatment zeroed, under a
   gradient-diffusion model, collocation + total-vorticity conservation
   *uniquely* yield $\dot\sigma_p = \nu_t/\sigma_p$ at fixed Γ — a
   per-particle turbulent viscosity in the existing CoreSpreading channel; a
   compensating Γ term is *forbidden*, not merely unneeded. Alvarez/Ning
   themselves state no formal stability analysis exists for the rVPM
   (B p. 4) — 018 §2 + this doc's §4 supply it.
2. **Scale selection** (§3 R1/R2, §4): dynamic-σ Smagorinsky is scale-free
   (no-go, proved); physically-constanted ν_eff (Squire/BL, δ≈2.3 at this
   Re_v — NOT "tens of ×") cannot arrest filament thinning (needs
   ~3000ν); the admissible primary candidate is the **fixed-filter form**
   $\nu_t = \nu + (\kappa\sigma_0)^2|S|$ (019 remedy 2 made local), which
   parks every strained core at $\sqrt5\,\kappa\sigma_0$ regardless of
   strain; proposed $\kappa = 1/\sqrt5$ (maintained scale = shed scale).
3. **Closure and stiff integrator are complements, proved** (§4.5): closure
   alone leaves Euler transverse flips (cap > 2/3 below ~1.1σ_stab);
   integrator alone leaves a super-exponential lagged-Z fidelity runaway;
   together — bounded attractor, no thresholds, no catastrophe (four-quadrant
   numeric exercise at σ₀=0.02R confirms, peak ΔtZ = predicted cap 2.41).
4. **Discussion-record verdicts**: item 1 SUPPORTED and sharpened (the gap is
   the zeroed E_adv, derivation-level); item 3 CORRECTED — the
   "Γ←Γ(1−3ΔtZ)" caricature drops the stretching source; the axial
   multiplier is (1+2ΔtZ), the spurious thresholds are transverse (2/3) and
   σ-flip (1), and the intensification itself is physical (Kelvin), so "blow-up
   is discrete not physical" splits into a discrete part (integrator fixes)
   and a physical part (only a closure fixes).
5. Phase-2 protocol §7 (three legs, thresholds fixed, in-flight 019 Campaign
   E runs named as primary data) and Phase-3 consistency tests §8 (T1–T8
   with tolerances) are PRE-REGISTERED per Ryan's scoping decision.

**Next**: gate discussion (κ default; §4.5 simultaneous-strainer caveat; |S|
estimator choice) → on explicit approval, Phase 2 per §7. No code, no runs,
no notebook were touched.

## Objective

Give the VPM a physically grounded channel for *sub-σ* strain: a subgrid
closure that acts on the core-size (and, coupled, the strength) equation so
that unresolved strain-thinning is balanced by modeled subgrid transport —
removing the σ-collapse/ignition mechanism at its root rather than guarding
against it, and making core evolution meaningful at rotor Reynolds numbers
where the current CoreSpreading model represents *laminar* diffusion at
molecular ν. Success means the simulation can capture the strain rates the
physics predicts, instead of the strain rates the stability limit permits.

## Why (grounding, all measured in 018/019 — do not re-derive)

- The live σ update is σ ← σ(1−ΔtZ) with Γ ← Γ(1−3ΔtZ); the composed map in
  y=σ² is linear with regimes 0<ΔtZ<2 (stable, fixed point → σ_eq=√(ν/Z̄))
  and ΔtZ>2 (unconditional ignition through the Γ equation, which viscosity
  never touches). `018_.../sigma_blowup_mechanism.md` — read it first.
- **SFS cannot touch σ**: its term in Z carries prefactor f/(1+3f) = 0 in
  the live formulation — the model has no subgrid channel for core-scale
  dynamics at all. CoreSpreading with molecular ν is the only σ-side
  physics, and it is laminar (Re_Γ = Γ_v/ν ≈ 2×10⁴ here — the real wake is
  not).
- Consequences measured: ignition boundary σ/R ≈ 0.030 at NT=36; boundary
  nearly Δt-invariant (contraction route, Δt^(−3/2) vs Δt^(−1)); mesh
  refinement TIGHTENS the boundary (finer strain sources); five live
  negative-σ observations inviscid. 019 adds the margin/regime-map
  measurement apparatus (M = max ΔtZ readout, `max_dtZ` column) and the
  fixed-point framing σ\* ≥ 2σ_eq(Z̄(σ\*)).
- Crude cousins already ruled on in 018 §6: σ floor alone RULED OUT (Γ
  route ignites anyway); σ-triggered merging via current trigger RULED OUT;
  per-particle ΔtZ guard identified as the mechanism-level *numerical*
  lever. This item is the *physical* version: replace "clip the update"
  with "model the missing physics that would have prevented the update from
  needing clipping."
- The global-constant version (uniform ν_eff inflation so σ_eq becomes the
  subgrid cutoff, `WAKE_NU_FACTOR`) is 019's remedy 2 — cheap, but uniform
  in space/time and a pure modeling assertion. 020 is its local,
  strain-aware, theoretically derived replacement.

## Discussion record (Ryan ⇄ agent, 2026-08-06) — the arguments this item rests on

Recorded here so the theory phase does not re-derive them from scratch;
each is a claim to *verify or refute*, not an assumption.

1. **The SFS channel gap (why the current SFS does not capture core
   turbulence).** The Alvarez–Ning SFS models subfilter *vortex
   stretching* — enstrophy transfer acting on the **Γ equation** of
   resolved elements. Physical core turbulence acts predominantly by
   *radial transport*: turbulent mixing diffuses core vorticity outward,
   setting core growth (measured tip-vortex ν_eff runs tens of × molecular;
   Squire / Bhagwat–Leishman ν_eff = ν(1 + a·Re_Γ)). In VPM variables that
   is a **σ-channel** effect. The channels do not intersect in the current
   code: (a) the σ equation is structurally closed to SFS (prefactor
   f/(1+3f) = 0, 018-CONFIRMED; corpse σ marched to the laminar floor with
   SFS on); (b) deeper — at our resolution the entire core is one particle
   (core radius ≈ σ), so there is no resolved field *inside* the core for a
   stretching model to act on; the core's internal turbulence is entirely
   subfilter and its effect (fatten the core, relax the swirl profile) is
   inexpressible by a Γ-side transfer term at any coefficient.
2. **The principled derivation route.** Redo the rVPM filtered-equation
   derivation keeping the σ-equation subfilter terms the original treatment
   dropped or zeroed, and see what closure falls out — rather than bolting
   on an ad-hoc ν_t. If the derivation naturally produces something
   Smagorinsky-like (ν_t ~ (Cσ)²|S|) or Bhagwat–Leishman-like (ν_eff(Re_Γ)),
   those literature models become validation anchors instead of competing
   options.
3. **The blow-up is discrete, not physical (contractivity observation).**
   Under sustained Z > 0 the continuous local dynamics decays everything:
   σ² ~ e^(−2Zt), Γ ~ e^(−3Zt), hence the velocity scale Γ/σ² ~ e^(−Zt)
   DECAYS. Regime-2 ignition is forward-Euler overshoot at ΔtZ > 2/3. An
   exponential local update (σ ← σe^(−ΔtZ), Γ ← Γe^(−3ΔtZ)) or local
   sub-stepping is pointwise unconditionally stable, preserves σ > 0
   exactly, and removes the ΔtZ ceiling with zero new physics. Open: full-
   map stability with lagged field-coupled Z; composition with
   CoreSpreading. This is Phase 2's stiff-integration test.
4. **Lever taxonomy (once stability-limited, what resolves more strain).**
   Stiff local integration → fixes ignition/ceiling, not sub-σ physics.
   Adaptive per-particle σ (012 machinery) → fixes cost, not physics.
   Δt → fixes shed-scale stability only; contraction route measured
   ~Δt-invariant; ~Δt^(−2) cost. Models (this item; uniform ν_eff cutoff;
   prescribed-core hybrids) → fix the physics of unresolved scales by
   modeling them. **The wall**: the 1/σ² strain runaway persists only while
   the core is under-resolved; below the physical core radius r_c peak
   strain saturates at ~Γ_v/(2πr_c²) — but r_c at Re_Γ ≈ 2×10⁴ is
   maintained by turbulence *inside* the core, so reaching it faithfully
   means DNS of the core. Short of DNS, refinement below the σ_eq band with
   molecular ν computes stable, affordable, *laminar* — wrong — equilibria.
   Every non-model lever terminates at DNS or a model.
5. **ν_eff-as-cutoff (the global-constant stopgap, 019 remedy 2).** Choose
   affordable σ_target, set ν_eff = (σ_target/2)²·Z̄(σ_target) so
   2σ_eq = σ_target: the fidelity fixed point closes by construction,
   strained cores equilibrate AT the resolution scale (bounding |Γ|/σ² and
   ΔtZ — stability as a byproduct), and the modeling claim ("sub-σ
   strain-thinning balances sub-σ turbulent mixing") is more defensible at
   these Re than molecular ν. Cost: uniform in space/time — fattens
   quiescent cores and adds bulk wake diffusion; needs a loads-sensitivity
   sweep (encouraging: ν-delta at L1 measured aerodynamically NULL). 020's
   closure is this idea made local and derived.
6. **Spatial-resolution taxonomy (why no refinement axis escapes).**
   σ is the true resolution length: refining it resolves strain ∝ 1/σ² and
   costs stability one-for-one (the σ_stab trade). h at fixed σ (PPS/
   overlap) improves quadrature of the σ-smoothed field, resolves NO sub-σ
   strain, and inflates the 010/012 conditioning burden. Blade-mesh
   refinement sharpens shed-Γ gradients and TIGHTENS the ignition boundary
   (018 measured: c249 ended min_sr −1.09, s80 −0.37) — a mesh-convergence
   study at fixed σ is silently a walk toward regime 2. Every axis that
   increases resolved strain moves toward ignition; the one that doesn't
   resolves nothing new.
7. **Diagnostics already in place (019).** The margin-curve scaling
   exponent p from M(σ) ∝ σ^(−p): p ≈ 2 ⟺ compact sub-σ strainer
   (unresolved; Γ_implied = 2πσ²M/Δt ≈ const ≈ Γ_v validates the filament
   picture); p ≈ 0 ⟺ ambient-dominated (resolved). This is Phase 2's
   measurement language for "the gap."

## Phase gates (each ends with a Ryan discussion + explicit approval before the next opens)

### Phase 1 — THEORY (discussion + derivation; no code)

Deliverable: a derivation document (phase file in a `020_.../` subdir)
that:

1. Surveys candidate closures with their provenance, at minimum:
   eddy-viscosity core spreading ν_t = ν + (C_s·σ)²|S| (Smagorinsky-type,
   per-particle, entering σ̇ = ν_t/σ); a strain-conditioned σ production
   term derived from the rVPM enstrophy/energy balance (extend the
   Alvarez–Ning filtering derivation so the σ-equation subfilter term is
   nonzero and *derived* — discussion-record item 2); Γ-side transfer
   consistency (does the closure also modify the −3ΔtZΓ term, and what
   conservation law says it must); vortex-method literature on turbulent
   core spreading (Squire/Bhagwat–Leishman ν_eff(Γ/ν) tip-vortex models,
   Winckelmans CSM, CSM+PSE hybrids); **and the stiffness-aware local
   integrator** (discussion-record item 3) as the mandatory non-model
   companion — derive its full-map stability conditions with lagged Z and
   CoreSpreading composition.
2. Derives the closed-loop discrete map (the 018 §2 analysis redone with
   the closure) and states the stability result: for what closure
   constants is regime 2 unreachable (ΔtZ bounded) for all σ, Γ, Δt — the
   design target is *unconditional* freedom from ignition, proved at the
   map level, not demonstrated empirically.
3. States conservation/consistency constraints: total vorticity, linear/
   angular impulse, energy budget sign (closure must be dissipative on
   enstrophy at sub-σ scales), Burgers-equilibrium recovery in the laminar
   limit (ν_t → ν as |S| → 0 must reproduce σ_eq physics), and the
   zero-strain no-op limit (quiescent wake unchanged).
4. Pre-registers the Phase-2 gap-demonstration protocol (runs, metrics,
   thresholds) and the Phase-3 consistency tests with tolerances.

Gate: theory discussed with Ryan; closure form + constants' meaning agreed;
falsifiable predictions written down (e.g. predicted stable-σ floor vs the
019 regime map).

### Phase 2 — GAP DEMONSTRATION (added by Ryan 2026-08-06; real-simulation evidence that the closure is worth building)

Objective: show from real simulations that there is a gap in the current
model — we cannot resolve strain beyond the stability limit at σ_stab
(unless stiff integration lowers that threshold — test it), and the
existing SFS is stretching-only, not viscous — so a subfilter *viscous*
model is what would let the simulation capture the strain rates the physics
predicts. Success gives Ryan a self-contained evidence pack justifying the
implementation investment.

Deliverables:

1. **The strain ceiling, measured.** From the 019 regime map + margin
   curve: show resolved max strain follows the under-resolved scaling
   (p ≈ 2, Γ_implied ≈ Γ_v) all the way to the ignition boundary — i.e.
   today's max resolvable strain is set by *stability* (Z_max ~ 1/(Δt) at
   σ_stab), not by physics. Against it, the physics-predicted target: peak
   strain at a physical core, Z_phys ~ Γ_v/(2πr_c²) with r_c from measured
   tip-vortex core sizes (literature + any DJI9443 data) — quantify the
   gap in decades. Mostly reuses 019 data; no or few new runs.
2. **The stiff-integration test (does the threshold move?).** Stage A
   (offline, cheap, no src change): re-integrate recorded collapse
   trajectories (`data/p018_L2_visc_forensics/death_trajectory.csv`; the
   ufront corpse wake-health rows) under the exact local map
   σ ← σe^(−ΔtZ), Γ ← Γe^(−3ΔtZ) and show ignition disappears pointwise.
   Stage B (in-code prototype, default-off/bit-identical/test-gated, same
   contract as every src change): one screen-class run below the current
   boundary (e.g. σ = 0.02R viscous NT=36) with the exponential local
   update — does it survive with bounded M where the stock update ignites?
   Outcome either way is evidence: threshold moves ⇒ stiff integration is
   a real stability lever and the *remaining* error is the laminar
   equilibria (fidelity gap isolated, closure still needed for physics);
   threshold does not move ⇒ the ceiling is field-coupled, and the closure
   is needed for stability AND physics.
3. **The SFS-null demonstration.** Show the existing SFS does not and
   cannot supply the missing channel: (a) derivation-level — the σ-update
   prefactor f/(1+3f) = 0 (cite 018, restate standalone); (b) measured —
   SFS-on vs SFS-off A/B at a collapsing σ (reuse 018 evidence where it
   exists; one targeted screen pair if not), showing the σ trajectory and
   ignition are unchanged: SFS is stretching-only, viscous-null.
4. **Evidence pack write-up**: the three results composed into a short,
   self-contained document (publishable-quality figures per the repo
   standard) whose closing claim is: *a subfilter viscous (σ-channel) model
   is the missing piece that would let the simulation approach
   physics-predicted strain rates* — with the measured ceiling, the
   integrator result, and the SFS null as its three legs.

Gate: Ryan reviews the evidence pack and explicitly approves proceeding to
implementation (or redirects — e.g. if stiff integration alone closes
enough of the gap for current purposes, implementation may be re-scoped).

### Phase 3 — IMPLEMENTATION (code + consistency-with-theory tests)

- FLOWVPM change (upstream repo), default-off, bit-identical when off,
  test-gated — same contract as every 018/019 src change. FLOWPanel driver
  exposes it by env knob; `_case_metadata.toml` and the banner carry it.
- Consistency tests = the Phase-1 pre-registered set, at minimum:
  single-particle analytic map vs theory (exact, rtol ~1e-10, mirroring the
  019 `max_dtZ` test pattern); Burgers equilibrium with ν_t(|S|) vs derived
  σ\*; conservation checks; the closed-loop stability bound exercised at
  ΔtZ that ignites the stock model (the 019 σ-ladder ignition fixtures are
  the ready-made harness); off-state bit-identity.

Gate: all consistency tests pass at pre-registered tolerances; deviations
from theory recorded, not patched silently.

### Phase 3a — MINIMAL ν_sgs PROTOTYPE (staged 2026-08-29; Ryan approved STAGING only)

Staged from the σ/VPM illustration campaign (plans/sigma_vpm_illustrations_20260827/),
where the Act-II data made the gap concrete: at σ_shed=0.0249R the viscous run
`scr_p019_s025v_rr` ignites at rev 8.11 (vs 7.83 inviscid) — viscosity arrests σ at
the equilibrium but nothing damps Γ. Research question (Ryan): can an SFS closure
with a viscous interaction term allow a SMALLER stable σ_shed?

- Prototype: implement the candidate closure σ̇_p = −h_σ s σ_p + (ν+ν_sgs)/σ_p with
  ν_t = ν + (κσ₀)²|S|, κ=1/√5 (κ per the Phase-1 downgrade: a resolution-cutoff
  modeling choice), in FLOWVPM; default-off env knob, off-state bit-identical
  (Phase-3 contract). Scope: the σ channel only — this is a stability probe, not
  the full Phase-3 implementation+test battery.
- Success criterion (pre-registered): stabilizes the `scr_p019_s025v` screen arm
  past rev 8.1 (its measured ignition on the post-08-24 stack; wake-health CSV in
  FLOWPanel plans/sigma_vpm_illustrations_20260827/wake_health_csvs/). One FLOWVPM
  edit + one ~2 h 64c screen run. Secondary readout: min σ/σ_shed drift rate vs
  the stock −0.154/rev.
- NOT run within the illustration campaign; Phase-2 gate rulings still stand —
  this stages a bounded experiment, it does not open Phase 3.
- Future-work note (recorded, explicitly NOT prototyped): stable small-σ operation
  likely requires particle refinement in two regimes — (a) diffusion regime:
  ~isotropic splitting driven by ν+ν_sgs spreading (anisotropy near walls TBD);
  (b) stretching regime: splitting along the time-averaged vortex-stretching
  direction. Ryan 2026-08-29.


### Phase 4 — VALIDATION (real simulation + comparison)

- Primary harness: re-run the 019 regime map (or its decisive subset) with
  the closure ON — prediction: previously igniting cells become stable with
  bounded M, and the σ-iterate procedure converges at smaller σ\*.
- Physical validation: DJI9443 hover CT + Γ(r/R) at σ below the stock
  boundary vs experiment (0.072) and vs the stock-model converged results
  (018 B0/L1 anchors); wing-class case from the 013/016/017 suite as the
  second flow family. Closure constants must NOT be retuned per case —
  fixed from Phases 1/3, or the item reports the failure.
- Deliverable: validation table + the "what would falsify this closure"
  statement; 019's procedure re-run with the closure to show σ-selection
  under the new model (cross-item deliverable, noted in 019).

Gate: Ryan approval of the validation verdict (positive or negative — a
clean negative with mechanism is an acceptable close-out).

## Relationship to other items

014 (Das/η selection) and 019 (σ selection) treat wake discretization
parameters as *model parameters* pending exactly this item: if 020 lands,
σ-convergence may become a real axis. 012 (overlap management) attacks a
sibling conditioning problem; a σ closure changes its inputs. 018 §6's
trim-before-merge reordering and ΔtZ guard remain worthwhile *independent*
of 020 (symptom-level robustness).

## Acceptance criteria (item level)

1. Phase 1 theory doc standalone-readable, with the discrete-map stability
   proof and pre-registered protocols for Phases 2–3.
2. Phase 2 evidence pack: measured strain ceiling vs physics target,
   stiff-integration threshold test (both stages), SFS viscous-null — each
   from real simulation data, composed into a self-contained document Ryan
   can present as justification.
3. Phase 3 implementation default-off/bit-identical/test-gated, consistency
   tests pass at pre-registered tolerances.
4. Phase 4 validation on ≥2 flow families with fixed constants, including
   the 019 regime-map re-run, closing with an explicit
   supported/refuted/scoped verdict.

## Log

- 2026-08-06 (night) — Phase 2 launched. **Recorded deviations from the §7
  pre-registration:** (1) Stage A had no tracked per-particle trajectories
  to re-integrate (018 retained only aggregates + two single-step
  snapshots) — reinterpreted as pointwise census + aggregate-trace
  re-integration (`scripts/p020_stage_a.py`), P-2 logic intact. (2) Leg-1
  Γ_implied band criterion FAILS as written (1/6 viscous rungs in
  [0.5,2]Γ_v); supplementary pre-onset M used for ignited rungs per 019
  P1's "last pre-ignition reading" — the stability-set verdict rests on
  the p-arm + boundary localization + the 019 ε-crossing at σ_stab.
  (3) Ignition-core trace exact-floor ratio 1.28 vs the 1.1 criterion
  (registered data source passes at 1.08); cause identified: Euler
  spuriously near-zeroes σ at ΔtZ≈1. (4) FLOWVPM vortex-ring/leapfrog test
  files unrunnable (pre-existing TestEnv/Pkg env rot, fails before FLOWVPM
  loads); all other suites pass local+cluster. (5) Cluster
  `FLOWVPM_merging.jl`/`FLOWVPM_splitting.jl` differ from local
  uncommitted WIP — cluster versions deliberately LEFT IN PLACE (they are
  what every Campaign E reference ran; comparability outranks syncing
  WIP). Implementation: `euler_exp` prototype + `WAKE_EXPINT` knob +
  `expint` kwarg landed and deployed (md5-verified); jobs 13065481
  (Stage B) / 13065482 (Leg-3 SFS-off) submitted.
- 2026-08-06 — staged (this file) per Ryan's instruction during the 019
  strain-resolvability discussion; no work started. Entry context: 019's
  remedies ranking placed this as the mechanism-level fix (remedy 3), with
  uniform ν_eff inflation (remedy 2) as its global-constant stopgap.
- 2026-08-06 (later still) — Phase 1 executed (theory doc drafted, reviewed,
  numerically checked; see Current status). Ryan scoping decisions taken
  during planning: targeted-extension derivation depth; Phase-2 protocol
  pre-registered now against the in-flight 019 Campaign E runs. Plan:
  `~/.claude/plans/work-on-020-per-zesty-star.md`.
- 2026-08-06 (later) — Ryan additions: discussion record written into the
  item (SFS Γ-vs-σ channel gap, contractivity/stiff-integrator observation,
  lever taxonomy and the DNS-or-model wall, ν_eff-as-cutoff, spatial-
  resolution taxonomy); new Phase 2 (gap demonstration) inserted between
  theory and implementation — objective: real-simulation evidence that the
  subfilter viscous model is worth building (strain ceiling measured
  against physics target; stiff-integration threshold test; SFS
  viscous-null); phases renumbered (implementation → 3, validation → 4).
- 2026-08-29 — Phase 3a (minimal ν_sgs prototype) STAGED per Ryan during the
  σ/VPM illustration session; staging only, no code, no run. See Phase 3a section.
