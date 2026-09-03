# Illustrating the σ/VPM stability models (019–020) for a technical audience

## Context

Ryan wants a set of visuals for a technical audience illustrating four things developed in BRAINSTORM 018–020: (a) the viscous core-spreading model, (b) the reformulated discrete VPM core-spreading mechanism (stretch+diffusion map), (c) the SFS model, and (d) the stability criteria (margin M, σ_stab, drift-arrest). Emphasis is on **why/where each model is needed — or not needed**. Decisions already made with Ryan:

- **Medium:** standalone TikZ figures per the house convention (`<name>.tex` + `<name>/*.csv`, pdflatex-compilable) **plus animations (GIFs)**, all reusable across lab notebook, paper, and slides.
- **Test case:** DJI9443 hover rotor OGE (all existing 018/019/020 data) + ParaView/pyvista renders of retained ignition-corpse VTPs; **optionally** IGE/multirotor wake renders from 022 for visual breadth (no new model evidence, no new simulations anywhere).
- **SFS framing:** "gap + candidate, caveated" — lead with the clean structural null; present the closure as well-motivated but unproven; the 1.9–3.3-decade gap shown only with its "conditional scenario estimate" caveat (binding 2026-08-12 audit correction).
- Deliverable of this task: the figures/animations plan recorded as a **journal entry for today (2026-08-27)** — the entry outlines the plan; figure production itself is follow-on work the entry's checklist tracks.

## 2026-08-31 mechanism-audit correction (binding)

The particle-level audit in `forensics_remaining/` and the full rVPM equations
supersede the scalar-circulation shorthand used below.  The live update is

$$
\boldsymbol\Gamma^{n+1}=\boldsymbol\Gamma^n+
\Delta t\left(\mathbf S-3Z\boldsymbol\Gamma-C\boldsymbol\epsilon\right),
\qquad
\mathbf S=(\boldsymbol\Gamma\mathbin\cdot\nabla)\overline{\mathbf u}.
$$

Thus $1-3\Delta tZ$ is only one term in a vector update; it is not the
particle's circulation multiplier, and $\Delta tZ=2/3$ is not an exact
stability boundary of the full coupled map.  Under aligned stretching the
live homogeneous Euler multiplier is $1+2\Delta tZ$; for a general frozen
gradient the update depends on the full $3\times3$ operator.  The corrected
geometric integrator likewise evolves
$\mathbf q=e^{L\Delta t}\boldsymbol\Gamma$, then rescales by
$r=\lVert\mathbf q\rVert/\lVert\boldsymbol\Gamma\rVert$; it is not the scalar
map $e^{-3\Delta tZ}\boldsymbol\Gamma$.

The $y=\sigma^2$ cobweb remains exact for the frozen-$Z$ explicit-Euler
core-size step followed by molecular CoreSpreading.  It is a **sigma-only**
map, not a unified $\sigma/\Gamma$ mechanism.  Consequently F1 panel (c), the
old joint-stability interpretation, and the $\sqrt{3/2}\,\sigma_{stab}$ claim
require a broader full-operator re-derivation and are withdrawn for audience
use.  The empirical $\sigma_{stab}$ scale may remain a same-$\Delta t$
initializer, but not an exact vector-stability boundary.

The completed 40-revolution `scr_p019_s038v_gpu40` rerun also overturns the
provisional Act-III arrest claim: it has a transient supercritical excursion
at step 875, enters persistent propagation at step 996, and later loses most
of the retained wake.  The existing `actIII_fixed.gif` supports only a clean
short window through about nine revolutions.  It does not certify
$\sigma^*/R=0.0381$ or drift arrest.  The production boundary/initializer now
requires a new long-horizon bracket rather than a silent coefficient patch.

All older mechanism bullets, act rows, and asset specifications below are the
historical production plan; where they conflict with this block or
`forensics_remaining/SUMMARY.md`, they are superseded.

## Source-of-truth summary (from BRAINSTORM scouts — verified against files on disk)

### Mechanism (018, feeds figures 1–2)
- Continuous viscous core spreading: σ̇ = ν/σ ⇒ σ(t)² = σ₀² + 2νt. Discrete code update `σ ← √(σ²+2νΔt)` (FLOWVPM_viscous.jl:161), hard floor σ_floor = √(2νΔt).
- Per-step composition: Γ ← Γ(1−3ΔtZ); σ ← σ(1−ΔtZ) with Z = (1/5)(Γ·∇ᵀ)U·Γ/|Γ|² (lagged); relaxation (Γ only); then viscous update; then merge-then-trim.
- **Key reformulation — 1-D map in y=σ²:** y_{n+1} = (1−ΔtZ)² y_n + 2νΔt. Regime 1 (ΔtZ<2): fixed point → σ_eq = √(ν_eff/Z̄) (measured 2.166 mm vs 2.165 mm predicted at ambient Z̄=3.2 s⁻¹; ν_eff=1.433418e−5). Regime 2 (ΔtZ>2): geometric blow-up; simultaneously Γ's (1−3ΔtZ) flips sign and amplifies ~3ΔtZ/step — **ignition is Γ blow-up on a floor-pinned σ** (corpse: σ pinned at 9.41e-5 m for 9 steps while max|Γ|/σ² grew ×3.5/step → velocity 733→37,088 m/s).
- Why viscosity is needed / not needed: at 2.2σ_eq the viscous-vs-inviscid aero difference is NULL (ε_Γ 0.156%, ΔCT −0.38%) — viscosity is a **stability term, not an aerodynamics term** at moderate σ; without it, inviscid runs produce **negative σ** (measured min ratios −1.4…−38; Gaussian-erf kernel odd in σ ⇒ acts like Γ→−Γ, immune to merging).
- Source doc: `BRAINSTORM/018_dji9443_hover_convergence_campaign/sigma_blowup_mechanism.md` (§1–§3 mechanism, §5 confidence table; internal 2026-08-06 corrections already folded in above).

### Stability criteria (019, feeds figures 3–5)
- σ_eq = √(ν_eff/Z̄) (fidelity floor) and σ_stab = √(Γ_v Δt/2π) (explicit-stability scale); initializer σ₀ = max(2σ_eq, σ_stab).
- Margin M = max ΔtZ (direct `max_dtZ` monitor column; level-difference proxy is DEMOTED, under-reads 3–20×). Target M ≤ ε = 0.2.
- Measured ignition boundary (screen, viscous): (0.0299, 0.0312]R; σ_stab=0.0311R lands on it to ~±2%. Survivor fit M ∝ σ⁻³·⁰.
- **Drift-arrest (v1.1):** M≤ε for 25 revs is NOT sufficient — σ₀=0.0356R held M∈[0.004,0.05] while min σ/σ_shed slid 0.80→0.13 (≈−0.065/rev in ln) into rev-27.6 ignition. Certification = M≤ε AND d ln(min σ/σ_shed)/dt → 0. Production boundary (0.0349,0.0381]R ≈ 1.15–1.25× screen boundary.
- Δt refinement RULED OUT (ignition at constant physical revs; attainable strain ∝ Δt^(−3/2) vs threshold ∝ Δt^(−1); confirmed by fid144).
- Failure is **tail-driven, not bulk**: P4b in-band case shows bulk σ/σ_eq medians static (1.017→1.022) while low tail deepens (p1 0.681→0.592, min→0.14σ_eq).
- Final result: σ* = 0.0381R.

### SFS (020, feeds figure 6)
- **Clean structural result (Leg 3):** shipped Alvarez–Ning SFS enters the σ-equation only via prefactor f/(1+3f) ≡ 0 in the live formulation (FLOWVPM_timeintegration.jl:147) — SFS acts on Γ only, structurally cannot supply σ-channel physics. Empirical twin null: SFS-on vs off ignite at step 286 vs 269 (5.9%, inside ±15% null band).
- Candidate closure (present as motivated, unproven): σ̇_p = −h_σ s σ_p + (ν+ν_sgs)/σ_p with ν_t = ν + (κσ₀)²|S|, κ=1/√5 (κ downgraded to a resolution-cutoff modeling choice). R1 no-go: any closure with the dynamic σ_p as mixing length is scale-free.
- Caveats (binding, from `phase_01_theory.md` top-of-file correction block + `phase_02_evidence_pack.md`): 1.9–3.3-decade strain gap is a **conditional scenario estimate**; Leg-2 "closure needed for stability AND physics" WITHDRAWN → corrected geometric integrator kills the discrete death mode but is not field-level sufficient; necessity of closure = inconclusive. Item STOPPED pending Ryan's Phase-3 ruling — illustrations must not imply Phase 3 happened.

## Mechanism principles (annotate figures with these)

- **The σ floor is natural, not a clamp**: the viscous update adds $2\nu\Delta t$ to $\sigma^2$ every step ⇒ no particle ends a step below $\sigma_{floor}=\sqrt{2\nu\Delta t}$; "floor-pinned" = per-step tug-of-war (strain crushes, diffusion re-inflates), measured exactly at $9.41\times10^{-5}$ m for 9 steps.
- **Two explicit maps share one (lagged) strain rate Z**: $\sigma\times(1-\Delta tZ)$, $\Gamma\times(1-3\Delta tZ)$ — Euler truncations of the exact exponentials $e^{-\Delta tZ}$, $e^{-3\Delta tZ}$; the thresholds ΔtZ=2/3, 1, 2 are *Euler stability boundaries*, not physics. Negative σ and the Γ sign-flip are integration artifacts (020-2R geometric integrator: zero −σ, crosses Euler's death point healthy — but still dies 1.3 revs later by exact-map σ→0 strain collapse ⇒ bug fix ≠ stability; integrator and closure are companions).
  - Regime 1 ($\Delta tZ$ small): strain–diffusion equilibrium $\sigma_{eq}=\sqrt{\nu/\bar Z}$; viscosity is a *stability* term, aero-null at ≥2σ_eq ($\varepsilon_\Gamma$ 0.156%).
  - Regime 2 (ignition loop): tiny σ ⇒ gradients ∼$|\Gamma|/\sigma^2$ ⇒ large Z ⇒ once $\Delta tZ>2/3$ the Γ multiplier leaves stability, Γ sign-flips and grows ×~3ΔtZ/step on floor-pinned σ ⇒ Z grows further ⇒ runaway (×3.5/step measured; 733→37,088 m/s in 1–2 steps).
- **Why it's inescapable in the current model — nothing damps Γ**: viscosity acts on σ only; relaxation reorients Γ without shrinking it; shipped SFS σ-prefactor $f/(1+3f)\equiv 0$. Viscosity blocks the σ-collapse route but leaves the Γ-growth route open.
- **What helps**: shed at $\sigma_0=\max(2\sigma_{eq},\sigma_{stab})$ so the whole population keeps $\Delta tZ\le\varepsilon=0.2$; certify with margin M *and* drift arrest (per-step M≤ε alone missed a 25-rev tail slide to the floor).
- **What doesn't help**:
  - Δt refinement — floor drops with Δt, attainable strain ∝$\Delta t^{-3/2}$ vs threshold ∝$\Delta t^{-1}$; ignition at constant physical revs.
  - Viscosity alone (σ route only) and shipped SFS (5.9% twin null).
  - Merging — negative/tiny-σ particles are merge-immune; merge-before-trim ordering hazard.
- **Confidence caveat**: the y=σ² map is exact; the particle-level "who ignites whom" chain is corpse-forensics inference (§5 confidence table), one run family.

## PRIORITY 1 — The three-mechanism story (create this first)

A three-act narrative, each act = one mechanism with (i) a ParaView blow-up animation, (ii) a mechanism diagram, (iii) a fix diagram + ParaView animation of the fixed run. This section outranks everything below; the F1–F7 inventory becomes its supporting material.

| Act | Mechanism (detail in "Mechanism principles") | Blow-up render (data source) | Mechanism figure | Fix figure | Fixed render (data source) |
|-----|-----------------------------------------------|------------------------------|------------------|------------|----------------------------|
| I — Strain has no counterweight | Inviscid: $\dot\sigma=-Z\sigma$ unopposed ⇒ σ collapse, Euler overshoot ⇒ **negative σ** (min ratios −1.4…−38), Γ→−Γ induced fields, merge-immune | Campaign-E inviscid ignition (e.g. σ/R=0.020–0.030 inviscid cell) — **archived on HPC, data hunt** | F2a `fig_viscous_roles_a` (thinning unopposed; no equilibrium without ν; −σ pathology) | F2b `fig_viscous_roles_b`: viscous equilibrium $\sigma_{eq}=\sqrt{\nu/\bar Z}$ (predicted 2.165 mm vs measured 2.166 mm); aero-null caveat (viscosity = stability term) | Matched viscous twin at same σ/R surviving screen — archived, data hunt |
| II — Explicit Euler leaves its stability region | Viscous: σ floor-pinned at $\sqrt{2\nu\Delta t}$, gradients ∼$\|\Gamma\|/\sigma^2$ ⇒ ΔtZ>2/3 ⇒ Γ sign-flips, ×3.5/step, 733→37k m/s | **LOCAL** — `data/p018_L2_visc_forensics/*.vtp` steps 1032–1041 (color by σ/σ_shed and \|Γ\|/σ²) | F1 `fig_sigma_map` (cobweb; ΔtZ=2/3,1,2 as *Euler boundaries* of exact exponentials) | F7 `fig_stage_b_r` (geometric integrator: zero −σ, crosses Euler death healthy — then exact-map collapse at step 243 ⇒ partial fix, motivates closure) | **REQUIRED (Ryan)**: 2R geometric run (job 13154223) render — survives past its Euler brother's step-237 death, *then* blows up at step 243; narrate "integration error fixed, instability mechanism remains" (alludes to the unresolved-strain gap / closure coda). Archived, data hunt; unpack if needed |
| III — Shedding below the stability scale | Even viscous + well-integrated runs ignite when σ₀ < σ_stab: population ΔtZ climbs, tail drifts to floor over revs (0.80→0.13 ln-slide, rev-27.6 ignition despite M≤ε) | **BOTH (Ryan)**: σ/R=0.0249 fast ignition AND σ₀=0.0356R slow drift-death (25 revs of M≤ε, then rev-27.6 ignition) — two renders; archived, data hunt | F3 margin curve + F5 `fig_drift_arrest` (M≤ε yet monotone drift) | F4 regime map with σ_stab=√(Γ_vΔt/2π) landing on measured boundary ±2%; initializer σ₀=max(2σ_eq,σ_stab) | σ*=0.0381R healthy 40-rev viscous production (CT 0.0723±0.0012) — archived, data hunt |

SFS remains the caveated coda (F6), not an act — its fix is unproven (020 stopped pre-Phase-3).

**Render pipeline (shared across all 7 animations — Act III has two blow-ups):** one pvpython/pyvista script, fixed camera on the rotor disk + near wake, particles as σ-scaled glyphs colored by σ/σ_shed (log scale, fixed range across acts) with a \|Γ\|/σ² inset trace; identical frame timing so blow-up vs fixed play side-by-side; → GIF (<10 MB) + kept PNG frames for slides. Prefer rendering **on the cluster** next to the archived VTKs (ship GIFs back) over pulling VTP series locally.

**Data hunt (gates the section):** hpc-monitor/hpc-storage query of `/nobackup/archive` INDEX.tsv + live run dirs for VTP series of: one inviscid ignition + viscous twin (Act I), 13154223 (Act II fixed), 0.0249R AND 0.0356R ignitions (Act III, both required), and the σ*=0.0381R production run (Act III fixed). Retention policy keeps newest 5 steps on /home — enough for a still, not an animation ⇒ archive unpacking likely needed: **confirm with Ryan before unpacking**, and note which acts survive if a series is gone (Act II blow-up is safe: local; Act II fixed and both Act III blow-ups are Ryan-required — if their VTP series were never written or are lost, flag immediately rather than substituting).

## Deliverables

(Supporting inventory — feeds the story acts above; also standalone-usable.)


New assets live in `~/Dropbox/research/notebooks/img/20260827_sigma_vpm_illustrations/` (house convention: `.tex` + same-named CSV dir per figure; GIF frames scripted). Reused 019/020 figures stay where they are; compiled PDFs (gitignored) get referenced/copied for the notebook.

### Figures (TikZ)
| # | Figure | Status | Content |
|---|--------|--------|---------|
| F1 | `fig_sigma_map` | **new, analytic** | y=σ² map cobweb: regime 1 fixed point (σ_eq) vs regime 2 blow-up; annotate ΔtZ=1 (sign flip), ΔtZ=2 (unconditional); inset: Γ multiplier (1−3ΔtZ). The unifying mechanism figure. |
| F2 | `fig_viscous_roles_a` / `_b` | **new, analytic + CSV** | SPLIT into two figures sharing one data dir, so each pairs with its Act-I render: `_a` (broken) = strain thinning unopposed, σ through zero, negative-σ pathology (min ratios −1.4…−38, corpse numbers from `data/p018_L2_visc_forensics/`); `_b` (fixed) = σ²=σ₀²+2νt diffusion and equilibrium σ_eq=√(ν/Z̄) (2.165 mm predicted vs 2.166 measured), annotated "stability term, not aero term" (aero null ε_Γ 0.156% at 2.2σ_eq). |
| F3 | margin curve | **reuse** `BRAINSTORM/019_sigma_selection_procedure/figures/fig_margin_curve.tex` | M vs σ/R log–log, ε, σ_eq, σ_stab, M∝σ⁻³. |
| F4 | regime map | **reuse** `.../fig_regime_map.tex` | outcome over (σ/R, viscosity), revs-to-ignition, screen vs production boundary. |
| F5 | `fig_drift_arrest` | **new, from run monitor CSVs** | M(t) staying ≤ε while ln(min σ/σ_shed) drifts −0.065/rev → rev-27.6 ignition (σ₀=0.0356R run); **REQUIRED** (load-bearing for Act III's drift-death blow-up, no longer optional); optionally pair with reuse of `fig_p4b_hist.tex` (tail vs bulk). **Data hunt needed** — locate that run's WakeHealthMonitor CSV (local harvest or /nobackup/archive); delegate to `harvester`/`hpc-monitor`. |
| F6 | `fig_sfs_channels` | **new schematic + reuse** `BRAINSTORM/020_sigma_aware_subgrid_closure/figures/fig_leg3.tex` | three-channel particle diagram (x/Γ/σ) with shipped SFS→Γ only, σ prefactor ≡0; candidate closure as dashed arrow into σ; beside the SFS-on/off ignition overlay. Optional: `fig_leg1.tex` with the conditional-estimate caveat printed on the figure. |
| F7 | exponential integrator fix | **reuse** `BRAINSTORM/020_sigma_aware_subgrid_closure/figures/fig_stage_b_r.tex` | Euler vs corrected geometric (exact-exponential) integrator health histories (job 13154223): geometric = zero negative σ, crosses Euler's step-237 death healthy (36 vs 16,095 m/s), then still collapses at step 243 via exact-map σ→0 strain runaway. The "−σ is an integrator bug; fixing it removes the discrete death mode but not the instability — integrator and closure are companions" figure. Do NOT reuse `fig_stage_b.tex`/`fig_stage_a.tex` (superseded mis-frozen prototype). |

### Animations (GIF)
| # | Animation | Source | Tool |
|---|-----------|--------|------|
| A1 | cobweb iteration of the y-map, side-by-side regime 1 / regime 2 (matched styling to F1) | analytic | matplotlib frames → gif (Python, stdlib+mpl; ≤4 threads) |
| A2 | ignition corpse: last 10 steps of the wake, particles colored by σ/σ_shed (or \|Γ\|/σ²), tail runaway visible | `data/p018_L2_visc_forensics/p018_L2_visc_wake1_particles.{1032..1041}.vtp` | pyvista or pvpython frames → gif. Healthy-run comparison: only a single healthy-step CSV exists locally (`map_healthy719.csv`); check archive for a healthy VTP series (hpc-monitor/hpc-storage INDEX.tsv), else static healthy frame vs animated death. |
| A3 (optional) | σ/σ_shed histogram tail deepening over the settle window | `BRAINSTORM/019_sigma_selection_procedure/p4b/*.csv` | matplotlib → gif |
| A4 (stretch) | IGE/multirotor wake beauty shot from 022 runs | VTKs likely archived on HPC (VTK retention policy: archive-first) — availability check first; skip if unpacking is heavy | pvpython on cluster or local after pull |

## Execution steps

1. **Journal entry (first, per the request):** propose a `# 20260827` day header + checklist in `~/Dropbox/research/notebooks/journals/20260821.md` (466 lines — no overflow; no 20260827 header exists yet) with subtitle `## σ/VPM model illustrations (018–020)`; **wait for Ryan's approval of header + verbosity before writing.** Entry content: the motivation, the four-topic outline, the figure/animation table above (self-contained: restate key equations and numbers, per notebook policy), caveat register (SFS audit corrections), and the checklist of production steps.
2. **Story data hunt** (delegate: hpc-monitor/hpc-storage) — locate the seven render sources (Act III: both 0.0249R and 0.0356R) per the PRIORITY-1 table; report availability; get Ryan's go-ahead for any archive unpacking.
3. **Render pipeline + Act II blow-up animation** (local VTPs, no dependencies) — build the shared pvpython/pyvista script here, since Act II data is already on disk.
4. Remaining act animations as data lands (cluster-side rendering preferred).
5. Create `~/Dropbox/research/notebooks/img/20260827_sigma_vpm_illustrations/`; build F1, F2a/F2b (analytic; CSVs emitted by a small Python script for any sampled curves), compile with pdflatex, verify against CSVs.
6. Data hunt for F5 (harvester: find the σ₀=0.0356R certification run's monitor CSVs; provenance map in `BRAINSTORM/019_sigma_selection_procedure/figures/provenance_appendix.csv`); build F5.
7. F6 schematic; recompile reused figures (F3, F4, F7=fig_stage_b_r, fig_leg3, fig_p4b_hist) to confirm they still build; copy PDFs/PNG renders for notebook embedding.
8. A1 (analytic gif), A2 (corpse gif — likely subsumed by Act II render from local VTPs; pyvista preferred — check it's installed, else pvpython), A3 if time.
9. A4/IGE stretch: ask hpc-monitor/hpc-storage what 022 VTK series survive and where; only proceed if cheap to render remotely or already local. Confirm with Ryan before any archive unpacking.
10. Offer notebook progress log with rendered figures embedded (`![...](../img/20260827_sigma_vpm_illustrations/...)`) — approval-gated as always.

## Verification
- Every `.tex` compiles standalone with `pdflatex` and reads only from its same-named CSV dir.
- Numbers on figures cross-checked against the source docs cited above (sigma_blowup_mechanism.md §5 confidence table; 019 item file; 020 phase_02_evidence_pack.md) — especially that no withdrawn 020 claim (old Stage A/B, unconditional decade-gap) appears uncaveated.
- GIFs play (open locally), file sizes sane (<10 MB), colorbars/legends labeled.
- Local compute only, ≤4 threads; no new simulations.

## Size estimate
Moderate-to-heavy: the PRIORITY-1 story adds seven ParaView animations, six of which need archived VTP series located (and likely unpacked) on HPC — that data hunt + cluster-side rendering is the dominant cost and is Ryan-gated. TikZ figures and the Act II animation remain cheap local work on existing data.
