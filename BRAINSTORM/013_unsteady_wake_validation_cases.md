# 013 — Unsteady Inviscid Panel Validation Cases (VPM / free-wake showcase)

## Purpose

Find published validation cases suitable for an **unsteady inviscid panel method**
that would exercise and show off a shed **free wake** — either the VPM (vortex-particle)
wake or a doublet-panel free wake (per user note, a doublet-panel free wake counts;
the reference panel code need not itself use VPM). Bonus value where **both an existing
panel code AND experiment** provide data to match.

These are exploratory suggestions for a demonstrable, quotable unsteady validation. The
brainstorm mission statement is rotor-hover $C_T$; the two strongest cases (Caradonna–Tung,
and by extension any rotor case) tie back to that, but the flapping/heaving cases are chosen
primarily as clean, well-documented wake showcases that FLOWPanel's unsteady path can hit
without dynamic-stall / viscous physics it does not model.

## Selection criteria applied

- **Inviscid-friendly:** attached flow, thrust/lift dominated by circulation + shed wake,
  not by massive separation or dynamic stall (rules out most McCroskey/AGARD-CT deep-stall sets).
- **Wake is the star:** thrust from a reverse Kármán street or a rolled-up tip/trailing
  wake — exactly what a free wake must capture and a rigid wake cannot.
- **Data to match:** experimental force/efficiency AND/OR an existing panel-method curve.
- **Reachable geometry:** FLOWPanel is 3D; 2D cases are run as high-aspect-ratio wings.

---

## Recommended cases (ranked)

### 1. Heathcote, Wang & Gursul (2008) — heaving rectangular wing  ★ top pick
- **What:** rectangular NACA 0012 wing in a water tunnel, **pure heave** (plunge), measuring
  thrust, lift and propulsive efficiency vs **Strouhal number**, with rigid vs spanwise-flexible
  variants; DPIV of the wake.
- **Why it fits:** canonical 3D flapping-wing dataset; thrust comes entirely from the shed
  reverse-Kármán wake — a direct free-wake showcase. Attached, inviscid-friendly at the rigid
  condition. Finite AR ⇒ genuinely 3D wake with tip roll-up (plays to FMM/VPM strengths).
- **Data to match:** **experiment (forces + efficiency + DPIV)** AND multiple **panel/VPM codes**
  have published curves against it — notably it is a standard validation case in the FLOW ecosystem
  (FLOWUnsteady / Alvarez & Ning VPM), so a like-for-like comparison is available. **Bonus box: both.**
- **Ref:** Heathcote, Wang & Gursul, "Effect of spanwise flexibility on flapping wing propulsion,"
  *J. Fluids and Structures* 24(2), 2008. (Rigid-wing subset is the cleanest inviscid target.)

### 2. Anderson, Streitlien, Barrett & Triantafyllou (1998) — heaving+pitching foil (2D)
- **What:** harmonically **heaving AND pitching** foil; force/power measured at Re≈40k, DPIV at
  Re≈1100. Efficiencies up to ~87%; classic **reverse Kármán street** visualizations.
- **Why it fits:** the textbook thrust-producing oscillating foil; the paper itself compares to
  **linear and nonlinear inviscid panel theory** and identifies the parameter range where inviscid
  agreement is good (weak/no LEV) — a ready-made "inviscid should work here" window.
- **Data to match:** **experiment (thrust/efficiency + wake PIV)** AND **inviscid panel theory
  curves in the same paper. Bonus box: both.** Run as a high-AR wing in FLOWPanel.
- **Ref:** Anderson et al., "Oscillating foils of high propulsive efficiency," *JFM* 360, 1998.

### 3. Wagner / Theodorsen analytic benchmarks — impulsively started & oscillating airfoil  (do first)
- **What:** step change in AoA (Wagner indicial lift build-up) and sinusoidally oscillating airfoil
  (Theodorsen lift/moment vs reduced frequency).
- **Why it fits:** the **exact** unsteady-wake sanity check — a correct shed wake must reproduce the
  Wagner lift-growth lag and the Theodorsen amplitude/phase. Cheap, unambiguous, and isolates the
  wake-shedding + convection logic before any experimental case.
- **Data to match:** **exact theory** (no experiment), and it is the standard verification target for
  unsteady doublet/source panel codes (Katz & Plotkin). Not a bonus-box case, but the correct
  **step 0** so the flapping-case discrepancies are attributable, not questioned.
- **Refs:** Wagner (1925); Theodorsen (NACA TR-496, 1935); Katz & Plotkin, *Low-Speed Aerodynamics*, ch. 13.

### 4. Caradonna–Tung hovering rotor — mission tie-in
- **What:** two-blade hovering rotor, blade **surface-pressure distributions** at several radial stations;
  a long-standing rotor-aero benchmark.
- **Why it fits:** directly serves the brainstorm's rotor-hover mission; validated in the literature with
  an **unsteady panel + (viscous) vortex-particle wake**, so the wake model is exactly this project's.
- **Data to match:** **experiment (pressures/thrust)** AND published **panel+VVPM** results. **Bonus box: both.**
  Caveat: hover, so it overlaps existing items 002–006; use it as the bridge from the flapping showcase back
  to $C_T$, not as the primary "unsteady showoff."
- **Ref:** Caradonna & Tung, NASA TM-81232, 1981; e.g. Tan & Wang, "Simulating unsteady aerodynamics of
  helicopter rotor with panel/viscous vortex particle method," *Aerospace Science & Technology*, 2013.

### 5. AGARD 445.6 pitch–plunge flutter wing — stretch / aeroelastic
- **What:** weakened NACA 0012 rectangular wing, classic **flutter-boundary** experiment.
- **Why it fits:** validated with **unsteady source+doublet panel free-wake** methods; the shed wake drives
  the unsteady loads that set the flutter speed — a strong free-wake test with a crisp scalar to match
  (flutter speed index vs Mach).
- **Data to match:** **experiment (flutter boundary)** AND **unsteady panel free-wake** results. **Bonus box: both.**
  Caveat: requires structural/aeroelastic coupling (modal model), i.e. heaviest lift of the five; incompressible
  panel limits it to the low-Mach subset. Park for later unless aeroelasticity is in scope.

---

## Additional 3D cases (added 2026-07-21, per user: need genuinely 3D, incompressible, high-ish Re so viscous effects are small, no stall / no leading-edge-vortex formation; the reference paper may include violating cases so long as ≥1 case meets these)

> Note on the earlier five: #1 (Heathcote heaving wing) and #4 (Caradonna–Tung) are already 3D and
> qualify; #2 (Anderson foil) and #3 (Wagner/Theodorsen) are 2D/analytic and are kept only as the
> low-effort verification rung. The five below are all 3D and screened against the attached / high-Re /
> no-LEV requirement.

> **Re > 200k screen (added 2026-07-21, per user):** clearing Re > 200k are **#11 marine propeller**
> (Re ≈ 5×10⁵–10⁶, *exactly incompressible* — best fit), **#4 Caradonna–Tung** hover (tip Re ≈ 1.9×10⁶;
> use the low-tip-Mach runs, e.g. M_tip≈0.44, so ~incompressible — higher-M_tip runs are the compressible
> violating case), **#7 MEXICO** (chord Re ≈ 3–6×10⁵), and **#6 isolated propeller** *if a full/large-scale
> prop is chosen* (model props can fall to ~10⁵). **Failing** the screen: #1 Heathcote, #2 Anderson,
> #9 gusts (all water-tunnel Re ≈ 10⁴–10⁵). So for Re > 200k the ranked picks are **#11 → #4 → #7 → #6**.

### 6. Isolated propeller in axial flight  ★ best 3D fit
- **What:** isolated propeller (e.g. APC-class or the propeller from propeller–wing free-wake studies)
  in axial inflow; measured **thrust / power vs advance ratio (J)**.
- **Why it fits:** **3D, incompressible, high blade-section Re, attached, no stall/LEV** across the normal
  operating J-range. The **helical trailing wake** is the whole story — a rigid wake mis-predicts induced
  inflow, a free wake (VPM or doublet) captures slipstream contraction/shear. Directly reuses this repo's
  `rotor_axial` machinery, so lowest implementation cost of the 3D set.
- **Data to match:** **experiment (thrust/power vs J)** AND **free-wake panel/lifting-line+VPM** curves
  (published isolated-propeller validations report thrust within ~5% of experiment). **Bonus box: both.**
  Violating cases available in the same references: very high disk loading / static (J→0) where sections stall.
- **Refs:** Gur & Rosen and the AIAA J. Aircraft "Higher-Order Free-Wake Method for Propeller–Wing Systems"
  (10.2514/1.C034720); ScienceDirect filament free-wake for propeller–wing (S1270963823006715).

### 7. MEXICO wind-turbine rotor  ★ strongest experiment + free-wake-panel pairing
- **What:** 3-bladed model wind-turbine rotor in a large wind tunnel; **blade surface pressures** at multiple
  radii and **PIV tip-vortex tracking**, at several tip-speed ratios (TSR).
- **Why it fits:** at design/high TSR the blades run **attached, high Re, incompressible, no stall** — exactly
  the target regime; the trailing/tip vortex system is a textbook 3D free-wake showcase (PIV gives wake
  geometry to match, not just forces).
- **Data to match:** **experiment (pressures + PIV)** AND a **free-wake *unsteady panel* model** published
  head-to-head against it — the single closest analog to what FLOWPanel does. **Bonus box: both.**
  The same paper's **low-TSR (high-wind) cases stall** — the required in-paper violating case.
- **Refs:** "Comparison and validation of BEM and free-wake unsteady panel model with the MEXICO rotor
  experiment" (ResearchGate 40736376); MEXICO/New-MEXICO database (Snel et al.).

### 8. Propeller–wing slipstream interaction
- **What:** tractor propeller ahead of a wing; measured wing/section loads and slipstream survey with the
  prop on vs off, at zero/low incidence.
- **Why it fits:** **3D, incompressible, high Re, attached** at the validation conditions; the propeller's
  **convecting helical wake washing over the wing** is a demanding unsteady free-wake test (rigid wake cannot
  reproduce slipstream deformation/shear over the wing) — a vivid VPM showcase beyond a single lifting surface.
- **Data to match:** **experiment (wing loads / slipstream)** AND **free-wake panel** results (papers note
  forces are over-predicted by the inviscid method but the free wake captures slipstream deformation well).
  **Bonus box: both.** Higher-incidence / high-thrust cases in the same work separate → in-paper violating case.
- **Refs:** "On the use of filament-based free-wake panel methods for preliminary design of propeller–wing
  configurations" (S1270963823006715); over-the-wing distributed-propulsion study (RG 322311321).

### 9. Wing–gust encounter family (three related cases)

A wing/airfoil at fixed low incidence meets a transverse (vertical) gust; the shed starting +
gust-induced vorticity is a clean, circulatory, inviscid-friendly unsteady-wake benchmark at small
gust ratio (no LEV). **Crucial 2D-vs-3D distinction for a Dirichlet 3D panel code:** the two classical
analytic benchmarks (Sears, Küssner) are **strictly 2D airfoil** theories, whereas only the UMD
towing-tank dataset is a genuine **3D finite wing**. Sears/Küssner are therefore *verification rungs*
(cheap analytic checks of the shed-wake logic on a high-AR wing); the 3D validation-with-data target
is #9a. Finite-AR effects (tip vortices, Helmbold lift roll-off) are exactly what a 3D free wake should
capture and 2D theory cannot — a feature, not a nuisance.

**#9a — UMD transverse gust, finite wing  ★ the 3D one (open data)**
- **What:** flat-plate wing translating through a discrete "sine-squared" transverse gust in a water
  towing tank; **lift time-history** vs gust ratio, **across an aspect-ratio sweep (AR≈1/2/4/∞)**.
- **Dimensionality:** **genuinely 3D** — the Jones group explicitly varied AR and corrected 2D theory
  with Helmbold to recover finite-wing lift; tip vortices / spanwise flow are resolved.
- **Data to match:** **experiment (open/downloadable)** + Küssner-w/-AR-correction theory + CFD.
  Small-gust attached cases match; **large gust ratio (up to 80% of freestream) separates** → in-paper
  violating case. **Bonus box: experiment+theory.** Caveat: plate is ~4.3% thick (fine for Dirichlet),
  Re≈40k moderate ⇒ keep to small gust ratios.
- **Refs (accessible):** Towne et al. open flow database (arXiv:2206.11801) — the AR-4 case + fields are
  downloadable; UC-Davis DAAL project page; closed-loop-pitch preprint (arXiv:2302.02902). Journal
  versions (Biler/Badrya/Sedky/Jones, AIAA J. 10.2514/1.J057646, .J059658, .J059678) are paywalled.

**#9b — Sears sinusoidal (convected) gust  — 2D analytic verification rung**
- **What:** airfoil in a streamwise-convecting **sinusoidal** vertical gust; amplitude/phase of lift vs
  reduced frequency = the exact **Sears function** $S(k)$ (the gust analog of Theodorsen).
- **Dimensionality:** **2D airfoil** (3D needs the Atassi generalization / unsteady lifting line — not a
  standard benchmark). Run as a high-AR wing to compare.
- **Data to match:** **exact theory (Sears)** + open-access wind-tunnel gust-generator validation.
  **Best analytic gust check** — unambiguous frequency response, cheap, no experiment strictly required.
- **Refs (accessible):** MDPI *Fluids* gust-generator paper (10.3390/fluids11030071, open); NASA NTRS
  20200000331 (public domain); Sears (1941); Atassi (1984) for the generalization.

**#9c — Küssner sharp-edged gust  — 2D analytic, public-domain**
- **What:** step (sharp-edged) transverse gust; indicial lift build-up = the exact **Küssner function**.
- **Dimensionality:** **2D airfoil**. Classic experiments are the NACA gust-tunnel program, but those
  measured whole-**aircraft-model** (3D) response, not a clean wing section — data is coarse by modern
  standards. Use mainly as the analytic step-gust companion to #9b's sinusoidal-gust check.
- **Data to match:** **exact theory (Küssner)** + public-domain NACA gust-tunnel reports.
- **Refs (accessible):** Küssner (1932); NACA gust-tunnel reports on NTRS/DTIC (public domain);
  sharp-edged-gust correction (10.2514/1.20298).

### 11. Marine propeller open-water (Wageningen B-series / DTMB 4119 / PPTC)  ★ best high-Re + exactly incompressible
- **What:** a marine screw propeller in uniform axial inflow; measured **open-water curves** — thrust
  $K_T$, torque $K_Q$, efficiency $\eta$ vs advance ratio $J$ (the ITTC-standard characterization).
- **Why it fits (Re > 200k target):** blade-section **Re ≈ 5×10⁵–10⁶+**, and being in **water it is exactly
  incompressible** (zero Mach concern — cleaner than any air case); at design $J$ the blades run **attached,
  no stall/LEV**. The **helical trailing wake** governs induced inflow ⇒ a direct free-wake showcase, and it
  reuses the repo's `rotor_axial` machinery.
- **Data to match:** **experiment** — open-water $K_T/K_Q/\eta$ are tabulated/public (**Wageningen B-series
  polynomials**; **DTMB 4119** and the **PPTC/SVA Potsdam** case ship open geometry + data) — AND **free-wake
  panel** results published head-to-head. **Bonus box: both.** In-paper violating cases: low $J$ (heavily
  loaded, sections stall) and cavitating conditions.
- **Refs:** "Marine propellers performance and flow-field prediction by a free-wake panel method"
  (Greco/Salvatore et al., free on Academia 29607623 / RG 268235190); "Open-Water Thrust and Torque
  Predictions of a Ducted Propeller System With a Panel Method" (RG 235577055).
- **Where to get the data (best→fallback):**
  - **PPTC (VP1304) — best, free download portal:** SVA Potsdam
    (https://www.sva-potsdam.de/en/potsdam-propeller-test-case-pptc/) → geometry + open-water report
    (SVA 3752, https://www.sva-potsdam.de/wp-content/uploads/2016/04/SVA_report_3752.pdf) + LDV wake-velocity
    report (SVA 3754, https://www.sva-potsdam.de/wp-content/uploads/2016/03/SVA-report-3754.pdf). Cite SVA.
    The LDV wake data validates the free wake itself, not just forces.
  - **DTMB/DTRC 4119 — classic panel benchmark:** experiment = Jessup (1989) PhD thesis (K_T/K_Q + LDV);
    geometry tabulated in every 4119 paper + OpenProp's built-in deck.
  - **Wageningen B-series — zero-download pre-check:** Oosterveld & van Oossanen (1975) closed-form
    K_T/K_Q(J, P/D, A_E/A_0, Z) polynomials in any marine-hydro text.
- **Existing panel-method results on PPTC (direct precedents to match):**
  - ★ **smp'17, "Simulate the PPTC propeller with a vortex particle–boundary element method"**
    (https://www.marinepropulsors.com/proceedings/2017/WC2-3.pdf, free) — **BEM panel blade + VPM wake on
    VP1304**, i.e. essentially FLOWPanel's architecture; authors report open-water K_T/K_Q "correlate very
    well" with the SVA experiment (exact J-range/%-errors not extracted — PDF text is binary; read the figures).
  - **smp'11 PPTC workshop proceedings** (https://www.marinepropulsors.com/smp/) — multi-code submissions
    incl. panel/BEM entries + case description (II-1_Barkmann.pdf); workshop summary tabulates each method's
    K_T/K_Q vs experiment across J (see where panel methods sit; snippets quote panel ~2–3% at design J).
  - "A systematic comparison between RANS and panel methods for propeller analysis"
    (free, Academia 28602069) — panel vs RANS vs experiment; calibrates expected panel-method accuracy.

### 10. 3D sudden-start / impulsively-started finite wing (3D Wagner)  — verification companion
- **What:** finite rectangular wing given a step in AoA (or impulsive start) from rest; **lift build-up in time**.
- **Why it fits:** **3D, incompressible, attached, no LEV**; extends #3 to a finite-span shed wake so the
  finite-AR lift lag and **trailing-wake roll-up** are exercised — validates the 3D shedding + convection logic
  before the experimental rotor/propeller cases.
- **Data to match:** **numerical reference only** (unsteady VLM / panel free-wake, Katz & Plotkin) — *no
  experiment*, so this is a **code-to-code / theory verification rung**, the 3D analog of #3, not a bonus-box case.
- **Refs:** Katz & Plotkin, *Low-Speed Aerodynamics*, ch. 13–14 (impulsively started finite wing, unsteady VLM).

---

## Suggested sequencing

1. **Wagner + Theodorsen (#3)** then **3D Wagner finite wing (#10)** — verify 2D then 3D shed-wake correctness
   against exact theory / code (cheap, do first; no experiment, pure verification rungs).
2. **Isolated propeller in axial flight (#6)** — best 3D fit (attached/high-Re/no-LEV) and lowest cost via the
   existing `rotor_axial` path; experiment + free-wake curves.
3. **Heathcote heaving wing (#1)** — headline 3D flapping free-wake showcase; experiment + FLOW-ecosystem VPM.
4. **MEXICO rotor (#7)** — strongest experiment + free-wake-*panel* pairing (pressures + PIV), attached at high TSR.
5. **Propeller–wing (#8)** / **transverse gust (#9)** — advanced 3D wake-interaction showcases.
6. **Anderson foil (#2)**, **Caradonna–Tung (#4)**, **AGARD 445.6 (#5)** — 2D check, mission tie-in, aeroelastic stretch.

**If picking one 3D experimental headline case:** #6 (isolated propeller) for lowest effort, or #7 (MEXICO) for the
closest free-wake-panel-vs-experiment literature match.

## Status / log

- **2026-07-21 (created):** Web-search survey; five candidates identified and ranked. Documentation-only.
  Next step (separate agent): pick #3 as the verification smoke test, then stand up #1 as the headline
  unsteady free-wake validation. No implementation yet.
- **2026-07-21 (3D expansion):** Per user, added five genuinely-3D cases (#6–#10) screened for
  incompressible / high-ish Re / attached / no-stall-no-LEV, each (except the #10 verification rung) with an
  in-paper violating case allowed. New 3D leaders: **#6 isolated propeller** (lowest effort, reuses `rotor_axial`)
  and **#7 MEXICO rotor** (closest free-wake-panel-vs-experiment match). Sequencing updated.
- **2026-07-21 (gust family):** Per user, expanded #9 into three gust cases with an explicit **2D-vs-3D**
  distinction and accessible references: **#9a UMD transverse gust** = the only *genuine 3D finite wing*
  (AR≈1/2/4 sweep, **open downloadable data** via Towne et al. arXiv:2206.11801); **#9b Sears sinusoidal**
  and **#9c Küssner sharp-edged** are *2D-airfoil analytic* verification rungs (run as high-AR wings), with
  open/public-domain validation (MDPI Fluids, NASA NTRS, NACA gust-tunnel). Dirichlet note: all are
  thin-but-finite bodies; thickness is a small predictable steady offset (see thickness caveat), not a shape
  distortion of Sears/Küssner/Wagner.
- **2026-07-21 (Re>200k screen):** Per user, screened for Re>200k and added **#11 marine propeller open-water**
  (Re≈5×10⁵–10⁶, *exactly incompressible* in water, open $K_T/K_Q$ data — Wageningen B-series / DTMB 4119 /
  PPTC — + free-wake panel validation; strongest high-Re incompressible 3D fit). Existing #4 (Caradonna–Tung,
  low-M_tip), #7 (MEXICO), #6 (full-scale prop) also clear the screen; #1/#2/#9 do not. High-Re ranked picks:
  **#11 → #4 → #7 → #6**.

### Sources
- Nonlinear unsteady VLM–VPM (rotorcraft), arXiv:2511.11430 — https://arxiv.org/html/2511.11430v2
- Willis et al., unsteady high-order panel method with vortex-particle wakes (MIT) — https://dspace.mit.edu/handle/1721.1/35592
- Tan & Wang, panel/viscous VPM helicopter rotor — https://www.sciencedirect.com/science/article/abs/pii/S1270963813001508
- Anderson et al. 1998, Oscillating foils of high propulsive efficiency (JFM) — https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/oscillating-foils-of-high-propulsive-efficiency/22E5CD028D92318AFC88ED104E55786B
- Heathcote/Gursul spanwise-flexibility flapping-wing propulsion — https://www.sciencedirect.com/science/article/abs/pii/S0889974613000819
- Unsteady panel method for flapping foil — https://www.sciencedirect.com/science/article/abs/pii/S095579970800129X
- Flutter calcs using unsteady source & doublet panel method (AGARD 445.6) — https://arc.aiaa.org/doi/10.2514/1.C037891
- Free-wake panel method, flexible wing flutter/gusts — https://www.sciencedirect.com/science/article/pii/S0889974623001238
- Unsteady lifting-line w/ Wagner function (3D wings), MDPI — https://www.mdpi.com/2226-4310/5/3/92
- Higher-Order Free-Wake Method for Propeller–Wing Systems (J. Aircraft) — https://arc.aiaa.org/doi/10.2514/1.C034720
- Filament-based free-wake panel for propeller–wing (ScienceDirect) — https://www.sciencedirect.com/science/article/pii/S1270963823006715
- Comparison/validation of BEM & free-wake unsteady panel with MEXICO rotor — https://www.researchgate.net/publication/40736376
- Multilevel panel method for wind-turbine rotor flows — https://www.researchgate.net/publication/311665719
- Experimental & computational investigation of transverse gust encounters (Biler/Sedky/Jones) — https://www.researchgate.net/publication/334690498
- Küssner's function in the sharp-edged gust problem — a correction — https://arc.aiaa.org/doi/10.2514/1.20298
- Towne et al., open database for reduced-complexity modeling of fluid flows (UMD gust data) — https://arxiv.org/pdf/2206.11801
- Experimental mitigation of large-amplitude transverse gusts, closed-loop pitch (open) — https://arxiv.org/pdf/2302.02902
- UC-Davis DAAL, wing response to large transverse gusts — https://daal.ucdavis.edu/wing-response-large-transverse-gusts
- Design, testing & numerical modelling of a low-speed wind-tunnel gust generator (MDPI Fluids, open) — https://doi.org/10.3390/fluids11030071
- Simulation & modeling of flow from a wind-tunnel gust generator (NASA NTRS, public) — https://ntrs.nasa.gov/api/citations/20200000331/downloads/20200000331.pdf
- Marine propellers performance & flow-field by a free-wake panel method — https://www.academia.edu/29607623/Marine_propellers_performance_and_flow_field_prediction_by_a_free_wake_panel_method
- Open-water thrust/torque of a ducted propeller with a panel method — https://www.researchgate.net/publication/235577055
- Caradonna & Tung, experimental & analytical study of a model rotor in hover (NASA TM-81232) — https://www.semanticscholar.org/paper/e773d0c6e2744d0475cf109985035246342d2961
