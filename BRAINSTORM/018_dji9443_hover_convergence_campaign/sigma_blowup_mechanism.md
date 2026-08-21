# The σ-collapse blow-up: mechanism, σ_eq, and why coarse σ is safe

**STATUS: COMPLETE (2026-08-06 early AM) — Tracks C (code audit), F (corpse
forensics), and H (in-flight measurements) all landed.** Track H delivered an
unplanned bonus: six tripwire-instrumented ignitions from the inviscid screen
batch, including the first live observations of **negative σ** — the exact
inviscid failure mode the code audit predicted.

Forensic tooling + data: `scripts/p018_sigma_forensics.py` +
`scripts/p018_sigma_forensics_strain.jl`; CSVs/figures in
`data/p018_L2_visc_forensics/` (ten corpse VTPs steps 1032–1041 retained,
351 MB). All Track F gates passed: kernel-vs-FD Jacobian 3.8e-10; recomputed
Z vs the VTPs' recorded velocity-gradient field correlation +1.000; the
ledger's "dt·Z ≈ 1e-3" reproduces as raw Γ̂·S·Γ̂ dt·Z = 2.98e-3 on the healthy
719 aged column — ÷5 for the live MM4 convention gives 6.0e-4, same order as
the recorded number (the prior forensics did not state its convention; treat
raw-vs-MM4 as the reconciliation).

Audience: human reader (Ryan; publication-appendix candidate). Every claim
carries a confidence label:
**CONFIRMED** = directly measured/proved, evidence cited ·
**SUPPORTED** = consistent multi-source evidence, not isolated ·
**HYPOTHESIS** = stated, untested.

---

## 1. What the code actually does to σ (CONFIRMED, file:line)

Each euler step, per particle (`FLOWVPM_timeintegration.jl:103–173`, driven by
`FLOWPanel_wake.jl:2007`):

1. position and Γ update (Γ gains a term **−3·Δt·Z·Γ**, `:156–158`);
2. **stretch: σ ← σ·(1 − Δt·Z)** exactly (`:161`), with
   $$Z = \tfrac{1}{5}\,\frac{(\Gamma\cdot\nabla^{T})\mathbf{U}\cdot\Gamma}{|\Gamma|^{2}}$$
   (live formulation f=0, g=1/5: `FLOWVPM.jl:99,134`; transposed scheme
   default). Z is *lagged* — evaluated from the previous velocity solve.
3. relaxation (Γ only — no σ write anywhere in `FLOWVPM_relaxation.jl`);
4. **viscous (when CoreSpreading active): σ ← √(σ² + 2νΔt)**
   (`FLOWVPM_viscous.jl:161`) — the exact integral of σ̇ = ν/σ, applied to the
   **post-stretch** per-particle σ;
5. FLOWPanel maintenance: merge FIRST, spatial trims SECOND
   (`FLOWPanel_wake.jl:1636–1646`).

**There is no floor, clamp, or abs() on σ anywhere in the live CPU path** of
either repo (the only `abs(sigma)` is GPU-only, `FLOWVPM_gpu.jl:168`). Trim
policies test |Γ| and position, never σ. Merging *skips* pairs with σ ≤ 0
(`FLOWVPM_merging.jl:491`). SFS (on in campaign runs via the driver's
`sfs_choice`, `rotor_hover_pressure_comparison.jl:420–509` — an earlier draft
of this audit wrongly said SFS was absent; corrected) modifies only Γ, and its
term in Z carries prefactor f/(1+3f) = 0 — **no SFS path can arrest a σ
collapse** (CONFIRMED).

## 2. The two regimes of the discrete map (CONFIRMED, derivation)

Compose stretch + viscous in y = σ²:

$$y_{n+1} = (1-\Delta t Z)^{2}\, y_n + 2\nu\Delta t$$

— linear, multiplier a = (1−ΔtZ)², offset 2νΔt > 0.

**Regime 1 — ambient strain (0 < ΔtZ < 2): stable, and the discretization
realizes the Burgers equilibrium.** Fixed point
σ\* = √(2ν/[Z(2−ΔtZ)]) → √(ν/Z) as ΔtZ→0. At the measured ambient
Z̄ = 3.2 s⁻¹, Δt = 3.086e−4 s (ΔtZ ≈ 1e−3): σ\* = 2.166 mm vs continuous
σ_eq = √(ν/Z̄) = 2.165 mm. Stability of the continuous balance: perturbing
σ = σ_eq + δ gives δ̇ = −2Zδ < 0 — strain-thinning and viscous-fattening
(which stiffens as 1/σ²) restore the equilibrium from both sides. e-fold time
1/(2Z̄) ≈ 0.16 s ≈ 500 steps. **With CoreSpreading active, σ additionally has
a hard per-step floor √(2νΔt) = 9.6e−5 m** — the viscous update is
positivity-restoring (only σ² enters).

Consequence (CONFIRMED, corollary): a viscous run can never end a step with
σ below 9.6e−5 m. The recorded 2.1e−5 m collapse was measured on the
**inviscid** run — consistent, since `Inviscid()` viscousdiffusion is a no-op
(`FLOWVPM_viscous.jl:58`) and nothing else touches σ.

**Regime 2 — ignition (ΔtZ > 2): unconditionally unstable, viscosity cannot
save it.** |1−ΔtZ| > 1 multiplies σ geometrically per step; simultaneously
the Γ update's −3ΔtZ·Γ term flips and amplifies Γ by ~3ΔtZ per step. At the
observed ignition velocities (u ~ 10³–3.7e4 m/s over σ ~ 2e−5 m scales,
Z ~ u/σ ~ 5e7–2e9 s⁻¹ ⇒ ΔtZ ~ 1.5e4–6e5), a two-step 733 → 37,088 m/s
explosion is exactly what this map produces. **The crossing is now measured
in the viscous corpse** (Track F): thin-core raw dt·Z median 0.27 / p95 > 2.4
at step 1041, with the Γ-side amplification observed directly (§3 item 3) —
CONFIRMED for `p018_L2_visc`.

Inviscid runs also produce **negative σ** for 1 < ΔtZ < 2 with a proved
downstream consequence: the Gaussian-erf kernel is odd in σ, so a negative-σ
particle induces velocity as if Γ → −Γ, and ζ terms carry 1/σ³ < 0
(`FLOWVPM_kernel.jl:51–53`, `FLOWVPM_viscous.jl:507,560`) — and such
particles are immune to merging (`merging.jl:491`).

## 3. The chain of events in a blow-up (narrative, now measured end-to-end)

1. **Sustained thinning (CONFIRMED)**: ambient strain contracts cores,
   σ(t) = σ₀·e^(−∫Z dt). In the inviscid runs the historical record is the
   aged column-edge layer (field-min σ 3.7e−3 → 2.1e−5 m over ~1000 steps,
   halving every ~180 steps). In the viscous corpse the collapsing
   population (863 particles below 0.25·σ_shed at step 1041) is centered on
   the **near-rotor slipstream edge** (median ax/R 0.46, rad/R 0.82,
   spreading ax/R −0.3…3.1) — at fine σ the contraction is already strong in
   young wake, not only the aged column. (Caveat: the corpse field is
   post-onset; the quiescent-phase spatial distribution was not retained.)
2. **Feedback hypothesis REFUTED (CONFIRMED refutation, with caveat)**: the
   proposed loop "thin cores strain each other" fails the direct test —
   recomputing Z at the collapsing subset with all thin sources EXCLUDED
   changes median |Z| by only ×0.905 (median dt·Z actually *rises*).
   **Ambient thick-core strain does the contracting; thin cores are victims,
   not mutual amplifiers.** (Measured at step 1041, a post-onset state ~8×
   hotter than healthy — a quiescent-phase test would need earlier
   snapshots, which no longer exist.)
3. **Ignition = Γ blow-up on floor-pinned cores (CONFIRMED for the viscous
   run)**: over retained steps 1032→1041, min σ was NOT falling — it sat
   pinned at 9.41e−5 m, which is exactly the CoreSpreading per-step floor
   √(2νΔt) = 9.6e−5 m predicted independently by the code audit (§2) — while
   **max |Γ|/σ² grew ×3.5/step for nine steps (1.7e4 → 1.25e9)**; the
   ignition core's |Γ| went 0.24 → 0.69 → 11.25 over the last three retained
   steps, and recorded max|u| jumped ×64 in the final step (4,330 →
   278,531 m/s). Same signature class as inviscid L2 (733 → 37,088 m/s).
   This is regime 2 of §2 running through the Γ equation (−3ΔtZΓ term) at
   the σ floor: thin-core dt·Z (raw) at death was median 0.27 with p95 > 2.4
   — far beyond any explicit-update comfort zone. The force monitor shows
   the run was already wild from step ~1027 (CFx swinging −0.22 → +0.52 →
   −3.9), so ignition preceded every retained snapshot.
4. **Death by symptom (CONFIRMED, arithmetic attached)**: the 278,531 m/s
   particle moves 86 m ≈ **723R in one step**; the merge hash sizes its
   arrays off the cloud bounding box with cell_size = the absolute r_merge =
   3.213e−4 m (`FLOWVPM_merging.jl:438–457`), giving n_cells ~ 10¹¹–10¹² and
   a ~TB `resize!` at `merging.jl:454` — instant OOM at MaxRSS 39.4 GB.
   Ordering trap CONFIRMED in the stack trace: merge (functional policy)
   runs before the GlobalCylinder trim that would have deleted the runaway
   (`FLOWPanel_wake.jl:1636`); the bbox was still bounded at step 1041, so
   the trim was holding until merge ran first on the post-advection state.
5. **Why the viscous rescue died (RESOLVED)**: CoreSpreading did exactly
   what §2 proves it does — it held σ at its floor and prevented σ→0. It
   cannot arrest the Γ-side of regime 2: at floor σ ≈ 1e−4 m, u ~ Γ/σ² still
   diverges as Γ grows, and no term in the stack damps Γ at ΔtZ ≫ 1.
   **"Viscosity wasn't enough" is now precise: the collapse route was
   blocked; the ignition route is in the Γ update, which viscosity never
   touches.** Also CONFIRMED: merging could not have intervened — only
   1/863 collapsing particles (0.12%) sat within the absolute merge radius
   (NN-dist median 3.1e−3 m vs r_merge 3.2e−4 m); even a 2σ_local criterion
   reaches only 2.3%. The collapsing set is structurally unreachable by the
   current merge trigger.

## 4. Why large-enough σ avoids blow-up

Three independent margins, all growing with σ:

- **Strain-resolution margin (CONFIRMED at screen scale, 2026-08-06)**:
  coarser ambient σ smooths the velocity field, capping the resolvable Z.
  Prior: L2's p99 per-step shrink rate ~6× L1's. New (inviscid 8-rev screens,
  harsh pulse-less startup, all tripwire-instrumented): **the ignition
  boundary sits exactly at ~σ/R 0.030** — σ/R 0.0299 (OV 2.4) survives 324
  steps at min σ/σ_shed 0.082 while σ/R 0.0291 (OV 2.0) ignites by step 284;
  σ/R 0.0199 ignites fastest (step 188, min_sr −2.4); σ/R 0.050 and every
  B0-carrier case survive comfortably (min_sr 0.10–0.38). Five of the six
  deaths recorded **negative σ** (min_sr −1.4 … −38, max|u| 13k–198k m/s) —
  the §2 inviscid prediction observed live. Caveat: the screen startup is
  harsher than production; the production inviscid pair at 0.0299R is alive
  past step 400 (contracting, min_sr 0.145), its viscous twin holds 0.381.
- **Peak-velocity headroom (CONFIRMED as scaling)**: |u|_peak ~ |Γ|/σ².
  At fixed shed circulation, doubling σ buys 4× headroom before any
  neighbor's ΔtZ approaches the regime-2 boundary ΔtZ = 2, i.e. before
  Z ≈ 2/Δt = 6.5e3 s⁻¹.
- **Viscous stiffness margin — with an important nuance (2026-08-06)**: the
  restoring term's core-balance share is (σ_eq/σ)² — ~4% at B0, ~21% at L1,
  ~34% at 0.03R. But **share of the core dynamics ≠ share of the forces**:
  `p018_L1_visc`'s settled window shows the viscosity delta at L1 σ is
  **NULL in circulation** (ε_Γ 0.156% vs inviscid L1; ΔCT −0.38%), because
  the strained thin-core population that viscosity governs carries little of
  Γ̄. So viscosity at moderate σ is a *stability* term, not an *aerodynamics*
  term; the b0_visc→L1_visc σ-pair reproduces the inviscid 3-lobe (+1.06%,
  ε_Γ 8.73%) — the σ-axis error mode is viscosity-independent. Viscosity
  becomes load-bearing for dynamics only near σ_eq (L2), where the inviscid
  model has no equilibrium at all.

Empirical stability boundary (CONFIRMED as observations): σ/R 0.0381 (L1)
completed repeatedly incl. viscous 40-rev; σ/R 0.030 = the screen-condition
inviscid boundary (survivor-side by a hair); σ/R 0.0193 (L2) died twice.

**A-priori initializer for the boundary (HYPOTHESIS, single-case check):**
modeling the worst structure as a filament of circulation Γ_v at core scale σ
(Z ~ Γ_v/(2πσ²)) and requiring ΔtZ ≲ 1 gives σ_stab ≈ √(Γ_v·Δt/2π). With
Γ_v ≈ 0.28 m²/s (measured peak bound circulation) and Δt = 3.09e−4 s this
predicts 0.031R vs the measured 0.029–0.030R boundary — right dimensional
group (Δt- and loading-dependent, σ_eq-free), flow-dependent prefactor.
General procedure: initialize σ at max(~2σ_eq, σ_stab), then verify/iterate
with the tripwire — Δσ/σ = −ΔtZ makes the fastest per-step core contraction
a direct readout of ΔtZ_max, and Z_max ~ 1/σ² makes the iteration converge
in a step or two. The "operate at ≥1.7–2σ_eq" numbers are THIS rotor's
instance of that procedure, not a universal rule.
~~The safe-operating recommendation (viscous production at σ ≥ ~1.7×σ_eq,
i.e. ≥0.03R, hedged at 0.035R) is now SUPPORTED; certification runs
`p018_ufront_n1_visc` (13058534, σ 0.03R) and `p018_ufront_s035_visc`
(13058988, σ 0.0349R, D=3.0) are in flight.~~ **Superseded — see the
2026-08-06 (later) correction at the end of this file: 13058534 DIED; the
1.7σ_eq rung is NOT safe, and the boundary lands where σ_stab (0.031R)
predicts, strengthening the initializer over the σ_eq-multiple rule.**

## 5. Confidence table

| claim | label | evidence |
| --- | --- | --- |
| σ update is σ←σ(1−ΔtZ), no floor/clamp anywhere (CPU) | CONFIRMED | `timeintegration.jl:161`; greps both repos |
| Live Z = (1/5)(Γ·∇ᵀU·Γ)/\|Γ\|²; SFS term in Z ≡ 0 | CONFIRMED | `timeintegration.jl:143–151`, `FLOWVPM.jl:99,134` |
| SFS active in campaign runs but cannot write σ | CONFIRMED | driver `:420–509`; `FLOWVPM_subfilterscale.jl:890,904` (net-zero) |
| CoreSpreading σ←√(σ²+2νΔt), post-stretch, floor √(2νΔt) | CONFIRMED | `FLOWVPM_viscous.jl:161`; composition `timeintegration.jl:161→171` |
| Discrete map stable iff ΔtZ<2; fixed point ≈ σ_eq at ambient | CONFIRMED | §2 derivation; numbers at Z̄=3.2 s⁻¹ |
| Ambient sustained Δt·Z ≈ 1e−3 in column-edge layer | CONFIRMED | prior ledger number reproduced (raw 2.98e−3 ÷5 MM4 ≈ 6e−4, same order); `strain_healthy719_gatea.csv` |
| Ignition = regime-2 crossing in the failing runs | CONFIRMED (L2_visc) | Track F death trajectory: thin-core raw dt·Z p95 > 2.4; \|Γ\|/σ² ×3.5/step ×9 steps |
| Ignition route is Γ growth at floor-pinned σ, not σ→0 | CONFIRMED | min σ pinned 9.41e−5 m across steps 1034–1041 while \|Γ\|/σ² → 1.25e9; ignition core \|Γ\| 0.24→11.25 |
| Corpse min σ = CoreSpreading floor √(2νΔt) | CONFIRMED (cross-validation) | measured 9.41e−5 vs predicted 9.6e−5 m — code audit and corpse agree independently |
| Feedback (thin cores strain each other) | **REFUTED** | thin-source exclusion changes med\|Z\| ×0.905 only — ambient thick-core strain dominates. Caveat: post-onset state |
| Collapsing population lives at near-rotor slipstream edge (this σ) | CONFIRMED (post-onset caveat) | med ax/R 0.46, rad/R 0.82; 863 particles < 0.25σ_shed |
| Merge-OOM = bounding-box symptom of one runaway | CONFIRMED (arithmetic) | 278,531 m/s ⇒ 86 m = 723R/step ⇒ n_cells ~ 1e11–1e12 at `merging.jl:454`; stack trace matches |
| Collapsing set unreachable by current merge trigger | CONFIRMED | 1/863 (0.12%) within absolute r_merge; 2.3% within 2σ_local |
| L2_visc was igniting at death (vs independent failure) | CONFIRMED — igniting | force monitor wild from step 1027; trajectory table above |
| Contraction rate grows with σ refinement; ignition boundary ≈ 0.03R (screen conditions, inviscid) | CONFIRMED | screen batch 2026-08-06: 0.0299R survives / 0.0291R ignites / 0.0199R fastest; phase_14 table |
| Negative σ occurs live in inviscid ignitions | CONFIRMED | five screens with min_sigma_ratio −1.4 … −38 in their wake-health CSVs |
| Viscosity delta at 2.2×σ_eq (L1) | CONFIRMED NULL (aero) | L1_visc settled 25–39: ε_Γ 0.156% vs inviscid L1, ΔCT −0.38%; earlier "positive trend" was scatter |
| σ-axis 3-lobe mode is viscosity-independent | CONFIRMED | b0_visc→L1_visc +1.06%, ε_Γ 8.73% ≈ inviscid pair |
| Viscous floor works live in production | CONFIRMED | ufront_n1_visc tripwire min_sr 0.381 vs inviscid twin 0.145 at comparable age (floor ≠ safety: the run later ignited through Γ — see 2026-08-06 later correction) |
| σ ≥ ~1.7×σ_eq viscous is blow-up-safe operating range | **REFUTED at 1.7σ_eq** (2026-08-06 later) | 13058534 (0.0299R ≈ 1.64σ_eq viscous) ignited and OOM'd; safe bracket now σ ≥ σ_stab ≈ 0.031R, carried by 13058988 (0.0349R, healthy in flight) and L1 0.0381R |

## 6. Mitigation implications (no implementation; Ryan decides)

What Tracks C+F now rule in/out (measured, not argued):

- A **σ floor alone is ruled out as the fix** — the viscous run already HAD
  an effective floor (measured pinned at it) and ignited anyway through the
  Γ equation. Any mitigation must address Γ growth at high ΔtZ, not σ.
- **σ-triggered merging via the current trigger is ruled out** — 0.12% of
  the collapsing set is within the absolute merge radius (2.3% even at
  2σ_local). Rescuing these particles by merging would require a different
  operator (absorb-into-nearest regardless of radius), not a retune.
- **Trim-before-merge reordering is now directly indicated** — the corpse
  proves the cylinder trim was holding the bbox until merge ran first on the
  post-advection state. Reordering converts the fatal OOM into a survivable
  deletion + tripwire event; cheap, low-risk, symptom-level.
- **A ΔtZ guard is the mechanism-level lever** — the instability is the
  explicit update at ΔtZ > O(1); a per-particle cap/sub-step (or clipping
  the −3ΔtZΓ increment) attacks the measured cause. Design question: what
  physics does clipping sacrifice, and does it matter at σ ≥ 0.03R where the
  guard would essentially never fire.
- **Operating at σ ≥ ~1.5–2×σ_eq** avoids the regime entirely (empirical
  bracket: L1 σ/R 0.0381 never ignites; L2 0.0193 ignites both inviscid and
  viscous; **2026-08-06 later: the 0.030R rungs are DOWN — both the inviscid
  13057253 and the viscous 13058534 ignited — so the multiple must be read
  as ≥ ~2×σ_eq, or better, as σ ≥ σ_stab; see the dated correction below**). With the uniform-d_front Das law,
  σ = 0.03R is clearance-admissible — the original *need* for σ = 0.02R is
  gone, so the cheapest complete mitigation may be: σ\* = 0.03R + the
  reordering fix + the tripwire, with the Γ-side guard as insurance.
- **Δt refinement is RULED OUT as a mitigation (measured 2026-08-06)**:
  Z_crit = 2/Δt rises, but the viscous floor √(2νΔt) falls, so attainable
  strain at floored cores scales Δt^(−3/2) — outrunning the threshold's
  Δt^(−1) — and empirically the same-σ half-Δt screen (`scr_ufdt_nt72`,
  σ/R 0.0291) ignited at essentially the same *physical time* as the NT=36
  case (245 vs 284 time units). Related subtlety, corpse-confirmed: strained
  cores equilibrate at σ\*(Z) = √(2ν/[Z(2−ΔtZ)]), not at ambient σ_eq —
  the corpse's min σ was the per-step floor, 22× below ambient σ_eq.

## 2026-08-06 (later) correction — certification outcome at 0.03R

Written after the STATUS-COMPLETE stamp above; supersedes the three
annotated spots (§4 closing paragraph, §5 rows "viscous floor works live" /
"σ ≥ 1.7σ_eq safe", §6 operating-range bullet). Source: campaign ledger
entries of 2026-08-06 (~09:50 MDT); recorded here by the 019 kickoff agent.

- **`p018_ufront_n1_visc` (13058534, σ/R 0.0299 ≈ 1.64σ_eq, viscous, N=1,
  D=3.4) DIED**: Slurm FAILED 1:0 after 11:59:15 with an explicit Julia
  `OutOfMemoryError` through `apply_particle_maintenance!` →
  `FLOWVPM_merging.jl:386` → `:454` (the same allocation line as L2_visc —
  the merge-bbox attribution no longer rests on a single stack trace).
  MaxRSS 67.0 GB. Its wake-health CSV (harvested to
  `data/p018_ufront_n1_visc/monitors/`, 494 rows) records the ignition
  directly: final row step 494, max_u 67,072 m/s, max|Γ|/σ² 3.2e8, with
  min_sr still +0.026 — the regime-2 signature (Γ blow-up at floor-held
  positive σ, the exact §3-item-3 route), not an independent memory failure.
- Its inviscid twin `p018_ufront_n1` (13057253) also died (step ~501,
  min_sr −22.9, max_u 33,977 m/s), so at 0.0299R production startup the
  viscous run merely outlasted the inviscid one; **viscosity delayed, did
  not save** — exactly the §2 regime-2 prediction.
- **Consequence for the safe-operating claim**: "σ ≥ ~1.7σ_eq viscous is
  safe" is refuted at its own boundary case. The surviving bracket is
  0.0349R (13058988, healthy in flight at step 381, min_sr 0.351) and
  0.0381R (L1, repeatedly completed incl. 40-rev viscous). Note
  0.0299R < σ_stab = √(Γ_v·Δt/2π) ≈ 0.031R < 0.0349R: the production
  viscous outcome falls on the side σ_stab predicts, while the σ_eq-multiple
  rule (1.7×) fails — **this upgrades σ_stab from "right dimensional group"
  toward the binding scale**, though still single-flow. Formalization and
  the cross-regime demonstration are item 019's scope
  (`BRAINSTORM/019_sigma_selection_procedure.md`).

## Cross-references

Phase/ledger: `phase_04_sigma_ladder.md` (σ ladder history, L2 deaths),
`phase_14_screen.md` (rate-vs-σ screens), ledger §2026-08-03(c) (original
forensics), §2026-08-05 (L2_visc death record). Notebook:
`journals/20260803.md` §20260805 (σ_eq entry). Data:
`data/p018_L2_visc_forensics/` [PENDING F].
