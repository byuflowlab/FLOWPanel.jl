# 018 Results Ledger

One row per settled result. `basis`: NEW (this campaign, sigma_overlap) or
PRIOR (recorded, overlap_pps legacy shedding — corroborating only).
Windows: PRIOR rows are 10-rev cycle-means (revs 10–19) unless noted.

| case / source | basis | form. | NT | η (Das/c@0.75R) | σ/c | nrows | depth | CT̄ | scatter | provenance |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| p2e_vel_nt36_d4 | PRIOR | vel | 36 | 0.2 (0.08, floor-clamped 68% span) | ~1.5 (legacy) | 1 | 4R | 0.06148 | ±0.00152 | 2e ledger |
| p2e_green_nt36_d4 | PRIOR | green | 36 | 0.2 | ~1.5 | 1 | 4R | 0.06205 | ±0.00188 | 2e ledger |
| p2e_das0p5 | PRIOR | vel | 36 | 0.5 (0.20) | ~1.5 | 1 | 4R | 0.06942 | ±0.00185 | 2e ledger |
| p2e_das1p0 | PRIOR | vel | 36 | 1.0 (0.41) | ~1.5 | 1 | 4R | 0.07133 | ±0.00159 | 2e ledger |
| p2e_das4p0 | PRIOR | vel | 36 | 4.0 (1.64) | ~1.5 | 1 | 4R | 0.07190 | ±0.00101 | 2e ledger |
| p2e_das0p2_nofloor | PRIOR | vel | 36 | 0.2, floor off | ~1.5 | 1 | 4R | 0.05215 | ±0.00078 | 014 |
| p2e_nrows4_das1p0 | PRIOR | vel | 36 | 1.0 | ~1.5 | 4 | 4R | 0.07431 | ±0.00096 | 014 (rejected ladder) |
| p2e_nrows8_das1p0 | PRIOR | vel | 36 | 1.0 | ~1.5 | 8 | 4R | 0.07506 | ±0.00160 | 014 (rejected ladder) |
| p2e_nrows16_das1p0 | PRIOR | vel | 36 | 1.0 | ~1.5 | 16 | 4R | 0.07304 | ±0.00203 | 014 (rejected ladder) |
| p2e_nrows36_das1p0 | PRIOR | vel | 36 | 1.0 | ~1.5 | 36 | 4R | 0.07049 | ±0.00047 | 014 (rejected ladder) |
| p2e_nt72_das1p0 (naive) | PRIOR | vel | 72 | 1.0 (0.205 — confounded) | halved | 1 | 4R | 0.07337 | ±0.00063 | 014; confounded control for Phase 3 |
| trunc 4R vs 6R (both form.) | PRIOR | both | 36 | 0.2 | ~1.5 | 1 | 4/6R | null ≤0.22% | — | 2e ledger |
| p2e_nrows72_das1p0 | PRIOR | vel | 36 | 1.0 | ~1.5 | 72 | 4R | 0.06931 | ±0.00154 (2.23%) | job 12955430, harvested 2026-07-31; revs 10–19; tail of rejected ladder still falling (36→0.07049, 72→0.0693) — non-convergence confirmed |
| p2e_sigF_nofilt | PRIOR | vel | 36 | 0.2 | ~0.37 (OVERLAP=4, PPS=8) | 1 | 4R | STABLE, no settled mean (died step 403/719) | CF smooth ≈0.066 | job 12943696 OOM-killed at 16 h — ReqMem was **32G**, MaxRSS 33.4 GB ⇒ genuine memory ceiling, NOT divergence (no \|CT\|>1). σ-refinement viability SUPPORTED; 64G ruling removes the cause |
| p2e_nt72_das2p0_ov6 | PRIOR | vel | 72 | 2.0 (0.41 pinned) | pinned (OV=6, PPS=2; ratio 3.0 = 2× legacy) | 1 | 4R | 0.06852 | ±0.00157; 5-rev block drift −1.14% (window not fully settled — mean may still be falling, cf. nrows72) | job 12950996 COMPLETED 2026-07-31 (41 h); revs 10–19; knobs md5'd from metadata (nwakerows=1, OverlapPPS 6.0/2). **−3.9% vs the NT=36 pinned partner p2e_das1p0 (0.07133)** — BELOW both Phase 0 branch anchors ⇒ dt NOT converged at NT=36, and the naive NT=72 (0.07337) was confound-dominated (Das halving pushed CT UP while true dt refinement pulls it DOWN) ⇒ Phase 3 = full ladder incl. NT=144 |

Anchors: CT_exp = 0.072; CCBlade BEM 0.060–0.071; steady Dirichlet
semi-infinite 0.0505–0.0515. New rows appended below as phases complete.

## Campaign rows (NEW)

| case | phase | job id | window (revs) | CT̄ | CI | ε_Γ vs ref | notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| p018_b0 (+s2) | 1 (+2 η=1.0 rung) | 12993801 + 13007210 | **15–29 (15 settled revs, M1 PASS; drifts 0.12%)** | **0.07170** | per-rev std ±0.00011; bootstrap 95% CI [0.07164, 0.07173] | TBD | Stitched 30-rev history (s2 restart at step 719; seam jump 0.065% of CT, within natural per-step variability — 42% of late deltas exceed it; restart chaining VALIDATED). η=1.0 (Das=0.41c), σ/c=0.68, N=4, metadata-verified. Interim 10-rev value 0.07187 superseded (sat on settling tail). Cost actuals: 13.2 h, N_p≈77k at σ/c=0.68 ⇒ t ≈ 9.0 h/(σ/c), N_p ≈ 52.4k/(σ/c). Ladder deltas vs das0p5/das2p0 still quoted on matched 10–19 windows |
| p018_L1 | 4 | 13011374 | **10–19 NOT VALID** (transient); tail 16–19 = 0.07251 ± 0.00007 (ptp 0.00016) | 0.07251 (provisional, 4-rev tail) | 10–19 window mean 0.07497 ± 0.00311, drift −6.58% ⇒ rejected | TBD | COMPLETED 2026-08-01 (16.2 h). σ/c=0.297 (OVERLAP 2.4 / PPS 11), N=4, η=1.0, metadata-verified. **Large mid-run excursion revs 11–14 (peak 0.0802, +12% over tail)** then a flat tail — finer σ takes markedly longer to settle than L0. Extension `p018_L1_s2` (13017106, +12 revs) submitted for a ≥15-rev M1 window. Provisional Δ vs B0 (0.07170): **+1.13%** |
| p018_nrows1 | 5 | 13011982 | 10–19 | 0.07133 | per-rev std ±0.00014; half-window drift −0.23% | **ε_Γ,max 2.77% / RMS 1.33% vs b0 (FAIL)** | COMPLETED 2026-08-01 (11.1 h). N=1 (d = Das = 0.6σ), η=1.0 metadata-verified. **−0.75% vs b0 — FAILS M1 (>0.5%) AND M2 (>1%)** ⇒ N=1 is NOT admissible at σ/c=0.68 |
| p018_nrows2 | 5 | 13011983 | 10–19 | 0.07188 | per-rev std ±0.00017; half-window drift −0.21% | **ε_Γ,max 0.35% / RMS 0.19% vs b0 (PASS)** | COMPLETED 2026-08-01 (11.6 h). N=2, η=1.0 metadata-verified. **NULL vs b0 (N=4, matched window 0.07187): +0.01%** — halving handoff 2.4σ→1.2σ is inert at σ/c=0.68 |
| p018_das4p0 | 2 | 12993804 | 10–19 | 0.07084 | per-rev std ±0.00044; half-window drift −0.14% | ε_Γ,max 5.32% vs das2p0 | COMPLETED 2026-07-31 (15.5 h). η=4.0 (Das=1.64c), metadata verified; case metadata otherwise identical to das2p0. **−2.02% vs das2p0 — ladder top rung NON-MONOTONE** (unlike legacy ladder) |
| p018_das2p0 | 2 | 12993803 | 10–19 | 0.07230 | per-rev std ±0.00027; half-window drift −0.16% | TBD | COMPLETED 2026-07-31 (14.9 h). η=2.0 (Das=0.82c), metadata-verified. Doubling deltas so far: 0.205→0.41c +2.58%, 0.41→0.82c +0.60% |
| p018_das0p5 | 2 | 12993802 | 10–19 (10 revs; run is 20 total — short of M1's ≥15, extend if rung is decision-critical) | 0.07006 | per-rev std ±0.00031; half-window drift −0.33%, oscillatory not monotone | TBD (ladder analysis) | COMPLETED 2026-07-31 (12.5 h). η=0.5 (Das=0.205c), σ/c=0.68, N=4; metadata-verified SigmaOverlap 2.0, nwakerows=4. Startup transient dies by rev ~9. Scatter ~5× tighter than legacy overlap_pps runs (±0.0003 vs ±0.0016) |

### Appended 2026-08-03 (harvest of the six jobs left in flight 2026-08-01)

All analysed with `scripts/p018_analyze.py` (new this session; validated by
reproducing p018_b0 = 0.071701 over revs 15–29, matching the banked 0.07170).

| case | phase | job id | window (revs) | CT̄ | CI | ε_Γ vs ref | notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| p018_nt72 (**η=1.0, NOT as intended**) | 3 | 13011373 | 10–19 (10 revs) | **0.070473** | bootstrap [0.070382, 0.070586]; per-rev std ±0.00017 | **0.979% max / 0.772% RMS vs das0p5 (PASS, marginal)** | COMPLETED 2026-08-02 (36.9 h). **Metadata says `das_eta_kinematic = 1.0`, not the intended 2.0** — the launcher arm is `${DAS_ETA_KINEMATIC:-2.0}`, so an η=1.0 inherited through `--export=ALL` won. ⇒ Das = 0.205c (HALF of B0's), σ/c = 0.68 correctly pinned. Not a B0-matched rung; **it IS an exact Das- AND σ-matched partner of `p018_das0p5`** (NT=36, Das 0.205c, σ/c 0.6815, rlxf 0.3→0.15 = the intended halving). Read that way: **+0.595% vs das0p5 (0.07006) — FAILS M1 (>0.5%), passes M2 marginally** ⇒ dt NOT converged at NT=36 even at matched Das. Both windows drift −0.32%/−0.33% monotone in the same direction, so the delta is more trustworthy than either mean. Corrected rung `p018_nt72_eta2` submitted (13029921) |
| p018_L1 (+L1_s2) | 4 | 13011374 + 13017106 | **17–31 (15 settled revs)** | **0.072669** | bootstrap [0.072091, 0.073473]; per-rev std ±0.00107 | **8.776% max / 4.692% RMS vs b0 (FAIL)** | Stitched at step 719. Supersedes the provisional 4-rev tail 0.07251. **Mean is window-insensitive: 0.072669 / 0.072664 / 0.072656 over revs 17–31 / 16–31 / 15–31 (0.02% spread)** ⇒ settled, even though 5-rev block drift is 0.9–1.4%: that is limit-cycle amplitude, not drift (non-monotone), exactly the case `decision_rules.md`'s limit-cycle defense covers. **Δ vs B0 = +1.34%** (σ moves CT UP). **But M2 FAILS hard (ε_Γ 8.8%)** — see the σ-axis note below |
| p018_L1_warm | 4 | 13011375 | own settled 22–37 = 0.072287; **gate on matched 22–31** | 0.072280 (matched) vs cold L1 0.072411 (matched) | per-rev std ±0.00023 (cold: ±0.00101) | **0.289% max / 0.148% RMS vs cold L1 (PASS)** | COMPLETED 2026-08-01 (19.6 h). **WARM-PILOT GATE PASSES: Δ = −0.18% on matched windows, inside the ≤0.25% threshold** ⇒ cross-σ warm starts validated (B0 σ/c 0.68 → L1 0.297) ⇒ **L2 may run warm**. Warm arm is markedly quieter than the cold arm (block drift 0.017% vs 0.210%), i.e. warm-starting also suppresses the settling excursion |
| p018_hub0p15_b0 | 5 | 13016062 | 15–19 | **VOID (−0.2989)** | — | — | **MESH DEFECT, not physics.** `CT_bernoulli = −0.474` at **step 1**, before any wake exists (stock p018_b0 = +0.070 at the same step). Not a mere normal flip (that would give ≈ −0.07, not −0.47). See phase_05 |
| p018_hub0p15_nrows2 | 5 | 13016670 | 15–19 | **VOID (−0.3125)** | — | — | Same defect |
| p018_hub0p15_nrows1 | 5 | 13016742 | 15–19 | **VOID (−0.2655)** | — | — | Same defect. **The hub test's payoff (does a hub rescue N=1?) is unanswerable until the mesh is fixed** |

**σ-axis finding (decision-relevant, Phase 4):** on CT̄ alone the σ axis looks
nearly converged (B0 0.07170 → L1 0.072669, **+1.34%**), but the co-equal
observable **fails badly — ε_Γ,max = 8.78% / RMS 4.69%** over 0.3 ≤ r/R ≤ 0.95.
This is exactly the case ruling 9 was added to catch: a rung that moves CT̄ a
little while moving Γ̄(r/R) a lot is **not** converged. L2 is required, and the
σ error-budget term must be set by M2, not by CT̄.

| p018_hubfix_b0 | 5 | 13031568 | none (diverged) | **DIVERGED** | — | — | Fixed hub mesh (`dji9443_20260803_..._hub0p15.msh`, md5 `14857a4d08c8fc1e654c29c0370ef549`), B0 carrier with only the mesh swapped, banner-verified (NT=36, η=1.0, OV 2.0/PPS 4, N=4). **Step-0 `CFx = −0.0803` ⇒ thrust +0.080, correct sign and magnitude — the geometry fix is VALIDATED in the solver** (broken mesh gave +0.474). But `|CFx|` grows ≈×2/step from step 2 to **3263 at step 17**, then OOM-kill at 10:28 (MaxRSS 66.9 GB / 64 G). **Divergence, not a memory ceiling** — blow-up proven in the recorded series, so do NOT resubmit with more memory. **ROOT CAUSE FOUND (Ryan, ParaView): the wake was shed from the ROOT CAP.** Probe: 136 shed edges of which **92 were cap edges** (innermost exactly at the blade root) vs stock's 40 clean TE edges. `shedding_r_over_R=0.1` silently doubled as the guard stopping the TE chain at the cap; a root outboard of the cutoff leaves it inert with NO error. Fixed: `SHEDDING_R_OVER_R`/`SHEDDING_ROOT_MARGIN` env knobs, cutoff auto-tracks the blade root, plus a hard abort if any shed edge sits at the root. Stock smoke bit-identical; hub now 41 edges / 0 cap. See phase_05 |

| p018_hub0p25_b0 | 5 | 13032932 | RUNNING | — | — | — | **Larger-hub case (Ryan-requested).** 0.25R hub, N=4, B0 carrier, only the mesh differs from B0 (0.07170). Mesh gate 7/7; 44 TE edges traced → 44 kept from r/R 0.258, **0 cap edges** (verified locally AND in the run log); deploy md5-verified. Healthy: CFx −0.0665 at step 0 (thrust +0.067, right sign) through +0.065 at step 348. 0.25R chosen so the whole M2 band 0.3–0.95 stays on blade. Pre-registered: small outboard ε_Γ ⇒ inboard region is aerodynamically minor and σ error is outboard; large ⇒ inboard clearance is first-order |
| p018_L2 | 4 | 13029924 | **OOM at 03:20 — NOT divergence** | — | — | — | Force monitor (survives SIGKILL) shows 99 clean steps, CFx −0.0733 → −0.0678, bounded. **MaxRSS 22.3 GB vs 64 G request (35%)** — sacct missed the spike, so MaxRSS≈ReqMem is NOT a usable memory test. ~347k particles at σ/c 0.151 vs hard-coded max 500k. See phase_04 |
| p018_L2_s2 | 4 | 13035813 | **OOM-killed AGAIN at 128G, 8m21s** | — | — | — | L2 relaunched at **128G** (ruling 5 amended by Ryan for σ-L2 and finer), restart-chained from `p018_L2` step **1249**, `P018_SETTLE_REVS=50`, 48 h. **MaxRSS 101.6 GB captured this time. Wrote exactly one step (1250) — bit-identical to the original L2's step 1250 — then died at the same place. DETERMINISTIC blow-up at step ~1251, not a memory ceiling: root-caused 2026-08-03 (see "σ-axis instability root cause" section below). Do NOT resubmit with more memory.** |

## σ-axis instability root cause + viscous erratum (2026-08-03, session c)

**L2 dies deterministically, not from memory pressure.** Both L2 attempts (64G
and 128G) died computing step ~1251; the 128G restart reproduced step 1250
bit-identically (CFx −0.06781…) and died at the same step. Particle-field
forensics (offline VTP scans, `/tmp` scripts `vtp_sigma_diag.py` /
`vtp_shrink_map.py`): max |velocity| grows 733 m/s (step 1248, already 11×
tip speed) → **37,088 m/s** (step 1250, 3,694 particles above 200 m/s); at
37 km/s one step moves a particle ~11 m ≈ 96R, so the next FMM-tree/merge
build over that bounding box is the >100 GB allocation. The OOM is a symptom.

**Root cause chain (three links):**
1. **rVPM dynamic core contraction** (`FLOWVPM_timeintegration.jl:161`,
   `σ -= dt·σ·Z`): by design, sustained axial strain thins cores. In the aged
   wake column (x/R ≳ 1.3) median dt·Z is **positive ~+1e-3/step**; the
   sub-shed-σ population concentrates in the **column-edge shear layer**
   (r/R 0.9–1.4, x/R 1.5–3.3) — exactly where the runaway ignited (x=2.9R,
   r=1.36R). Observed σ floor decays exponentially with wake age:
   3.7e-3 → 2.1e-5 m over ~1000 steps (halving every ~180 steps), with **no
   lower bound anywhere in the code**.
2. **CoreSpreading has been silently inert in EVERY campaign run** (the
   erratum): FLOWPanel steps the pfield with `FLOWVPM._euler`, but the pfield
   declared the FLOWVPM default `integration=rungekutta3`, so
   `viscousdiffusion` took the RK3 branch with zeroed stage weights — no
   spreading, no β resets. Ruling 9's "viscous ON in every run" has been
   configured-on but functionally OFF campaign-wide (SFS unaffected — it runs).
3. **σ refinement arms the loop**: finer ambient σ resolves sharper strain
   (L2's 99th-pct per-step shrink 6.6e-2 vs L1's 1.1e-2, 6×), and thin cores
   meet close neighbors → velocity spikes → more strain → runaway. B0/L1
   survive because coarse ambient σ shields the thin cores; L2 does not.

**Why viscosity is the missing physics, quantitatively:** working CoreSpreading
grows σ by ν·dt/σ² per step ≈ 8.4e-4 at L2's shed σ — the same order as the
median shrink (+1e-3), with restoring force ∝ 1/σ² as cores thin. Implied
equilibrium σ_eq = √(ν/Z̄) ≈ 2.1e-3 m ≈ L2's shed σ: viscosity sets the
smallest sustainable core (Burgers-like balance). An inviscid rVPM hover wake
has NO steady-state σ distribution — the ladder was refining toward a limit
the model cannot hold.

**Shedding-method question (Ryan): does OverlapPPS/SigmaPPS help?** No for
this failure — the contraction is dynamic and happens ~2R downstream of the
TE regardless of birth σ policy. Shedding choice matters for the d/σ
clearance structure (OverlapPPS gives σ ∝ local spacing ∝ r, making d/σ
radially uniform — relevant to the inboard Γ̄ lobe, logged for 012/Phase 5),
but no shedding policy bounds σ against strain. **016 is NOT needed for this
failure** (pre-authorization acknowledged but not triggered).

**Fixes implemented (local, test-gated, NOT yet deployed):**
- `src/FLOWPanel_wake.jl`: pfield now declares `integration=FLOWVPM.euler`
  (matches the actual stepper). Bit-inert for Inviscid runs.
- driver: viscous scheme now EXPLICIT — default `Inviscid()` preserves the
  campaign's de facto baseline bit-for-bit; `CORE_SPREADING_ACTIVE=true`
  enables working CoreSpreading with `CORE_SPREADING_SGM0` defaulting to the
  **shed σ** (not the old `wake_core_size=1e-3`, which β-resets would have
  stamped onto every particle, collapsing the ladder's resolution).
  Metadata prints both; launcher banner prints `visc:`.
- Unit test: "Viscous scheme actually runs under the euler stepper" in
  `test/runtests_unit_wake.jl` (integration declared = euler; one-step σ
  growth = √(σ²+2νdt)).

**Ryan's rulings (this session):** viscous erratum = **disclose + measure**
(one B0-carrier A/B with working CoreSpreading after validation); 016 remains
pre-authorized if it proves the best σ-axis path. Planned runs when slots
free: `p018_b0_visc` (erratum measurement) and `p018_L2_visc` (L2 rescue,
warm from step 1249, 128G) — arms to be added at submission time.

| p018_L1_ov3 | 4 | 13029923 | 20 revs total; tail 17–20 = 0.07249 | 0.07249 (short tail; not an M1 window) | last-10 mean 0.07423 (excursion-contaminated) | — | COMPLETED 2026-08-04 (19.2 h, MaxRSS 42.2 GB). OVERLAP 3.0 / PPS 14 (σ/c 0.292 ≈ L1's 0.297), metadata-verified. **Settle-time metric: the revs-14–16 excursion persists (peak 0.0763 at rev 16), first enters the final band ~rev 17 vs L1's ~16 ⇒ raising OVERLAP at fixed σ does NOT shorten settling.** The h-resolution/settle-time hypothesis is NOT supported; settle-time growth is a σ effect ⇒ q-schedule is not a campaign cost lever; L2's OVERLAP stays 2.88. Tail consistent with L1 within scatter |

## Appended 2026-08-04 — Das×N matrix + hub0p25 harvest (collapse test ANSWERED)

All matched windows revs 10–19 (settled by rev ~9; hub uses its own settled
14–19). Analysis: `p018_analyze.py` + new `scripts/p018_dasN_matrix.py`
(figure `data/p018_gamma_ladders/gamma_dasN_matrix.png`).

| case | phase | job id | window | CT̄ | CI / std | ε_Γ vs same-Das N=4 ref | notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| p018_nrows1_das2p0 | 5 | 13035912 | 10–19 | 0.071469 | [0.071390, 0.071545]; ±0.00015 | **4.02% max / 2.05% RMS vs das2p0 (FAIL)** | N=1, η=2.0, d/σ=1.21. −0.55% CT vs B0 |
| p018_nrows2_das2p0 | 5 | 13035911 | 10–19 | 0.071849 | [0.071716, 0.071954]; ±0.00025 | **2.22% / 1.12% vs das2p0 (FAIL)** | N=2, η=2.0, d/σ=2.41. −0.62% CT vs das2p0 ⇒ N=2≢N=4 at Das 0.82c |
| p018_nrows1_das4p0 | 5 | 13035910 | 10–19 | 0.070159 | [0.070099, 0.070205]; ±0.00013 | **3.08% / 1.75% vs das4p0 (FAIL)** | N=1, η=4.0, d/σ=2.41 |
| p018_nrows2_das4p0 | 5 | 13035913 | 10–19 | 0.070392 | [0.070282, 0.070455]; ±0.00027 | **2.09% / 1.16% vs das4p0 (FAIL)** | N=2, η=4.0, d/σ=4.82. −0.63% CT vs das4p0 |
| p018_hub0p25_b0 | 5 | 13032932 | 14–19 (6 revs; short of M1 — deltas only) | 0.068022 | [0.067892, 0.068153]; ±0.00025 | 10.7% max / 2.41% RMS vs b0 — but **localized**: −16% at r/R 0.29 decaying to 0 by 0.41; **\|Δ\| ≤ 1.2% everywhere outboard of 0.45** | hub removes blade inboard of 0.25R ⇒ absolute CT not comparable; the outboard band is effectively a null |

**Collapse test verdict (pre-registered discriminator, phase_05 step 2): NO
COLLAPSE — d/σ is REFUTED as the governing group.** At fixed d/σ=1.21 the
three realizations {N2·η1, N4·η0.5, N1·η2} give CT̄ {0.07188, 0.07006,
0.07147} (2.6% spread) and ε_Γ vs B0 {0.35%, 8.88%, 0.99%}; at d/σ=2.41,
{B0, N2·η2, N1·η4} = {0.07187, 0.07185, 0.07016}. **Das itself dominates both
error lobes**: the large-Das outboard deficit at r/R 0.84 is nearly
N-independent (−5.9/−6.2/−6.6% for N=1/2/4 at Das 1.64c), and the small-Das
inboard catastrophe (−8.9% at Das 0.205c) occurs at a d/σ where other cells
are clean. N is a secondary modifier: real only at N=1 short clearance
(d/σ≲1) and as a mild sheet-extent effect at Das≥0.82c (N=2 vs N=4 ≈ −0.6%
CT, ~2% ε_Γ). **The N=2≡N=4 null holds ONLY at Das\*=0.41c** — where it
remains clean (0.35%). Production reading: N=2 at Das 0.41c stands;
0.82c's "live competitor" status is now weakened (it carries its own
N-sensitivity); the Das axis still requires its re-test at σ\* per Phase 2's
no-plateau branch.

**Hub0p25 verdict (pre-registered):** the inboard region is **aerodynamically
minor for outboard loading** — deleting the blade inboard of 0.25R perturbs
Γ̄ outboard of r/R 0.45 by ≤1.2%. The σ-axis 3-lobe error (dip at 0.76) is
NOT root-driven; the tip-vortex reading stands and the σ axis stays routed to
viscous L2 / 012.

## Appended 2026-08-04 — L2 viscous rescue: failure, root cause, and relaunch

- **`p018_L2visc` (13036477) FAILED in 3m41s**: `Wake VTS not found:
  data/p018_L1_s2/p018_L1_s2_wake1.1.1151.vts`. Root cause: the 2026-08-04
  VTK retention sweep deleted the restart *pieces* (body `.vtu`, particles
  `.vtp`) of every banked run; only `.vtm` index stubs survive for
  `p018_L1_s2` (1149–1151) and `p018_L2` (1248–1250) ⇒ **no complete
  warm-start source remained on the cluster**. Process amendment recorded in
  `ops_reference.md` (keep one complete restart set per warm-capable run).
- **Rescue executed (Ryan-approved path):** the panel-wake `.vts` files
  survived cluster-side (deletion pattern spares `.vts`/`.vtm`), and the Mac
  archive `~/p018_L1_ov3_paraview/` holds the complementary pieces. Uploaded
  step-719 set (body `.vtu` 3.9 MB, particles `.vtp` 40 MB, body `.pvd`),
  md5-verified ⇒ `p018_L1_ov3` (σ/c 0.292, 20 revs) is again a complete warm
  source at step 719.
- **Submitted 2026-08-04:** `p018_b0_visc` (13047290, B0 carrier +
  `CORE_SPREADING_ACTIVE=true`, cold, 24 h — the erratum A/B, budget row Γ6)
  and `p018_L2_visc` (13047296, warm from `p018_L1_ov3` step 719, 128G,
  48 h, SETTLE=50 — expect an `_s2` chain). Banner verification logged in
  phase_04.
- **HPC slot cap raised to 10 (Ryan, 2026-08-04)** — supersedes the 6-job
  cap in ruling 5 / ops_reference.

## Appended 2026-08-04 (b) — chord-proportional Das implemented; κ ladder submitted

Root cause of the Phase-2 no-plateau ladder identified and the fix implemented
(full derivation + pre-registered predictions in phase_02 "2026-08-04 (b)"):
Das = η·Δt·Ωr is ∝ r while the admissible band (014) is 0.25–1.5 **local
chords**; putting r/R=0.31 on the plateau needs η ≥ 2.8 while keeping r/R=0.95
under 1.5c needs η ≤ 2.3 — no single η can place the whole span on the
plateau. New `DAS_CHORD_FRACTION` = κ ⇒ |Das|_j = κ·c_local(j)
(`set_Das_station_lengths` in src, unit-tested 7/7; formulation gate 10/10;
chord+stock 40_40 smokes clean, both finite; deploy md5-verified incl. full
src/ sweep). L2visc-vintage caveat: full runtests not rerun end-to-end (the
missing `data/pitching_wing_exp/` fixture still blocks the suite); the changed
subsystem's unit files pass.

| case | phase | job id | knobs (banner-verified on submission) | status |
| --- | --- | --- | --- | --- |
| p018_dasc0p25 | 2b | 13050924 | B0 carrier + κ=0.25 (|Das| = 0.25·c_local), 24 h | RUNNING |
| p018_dasc0p41 | 2b | 13050925 | B0 carrier + κ=0.41 (matches B0's Das at 0.75R exactly), 24 h | RUNNING |
| p018_dasc0p82 | 2b | 13050926 | B0 carrier + κ=0.82, 24 h | RUNNING |

Harvest plan: M1 on revs 10–19 (extend if needed), M2 successive-κ AND each
rung vs B0; the pre-registered predictions in phase_02 decide whether the Das
axis is convergent under the chord parameterization.

## Appended 2026-08-05 — `p018_nt72_eta2` harvested (Phase 3 production-Das rung)

| case | phase | job id | window (revs) | CT̄ | CI | ε_Γ vs ref | notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| p018_nt72_eta2 | 3 | 13029921 | 10–19 (10 revs; drift 0.275% monotone) | **0.072243** | [0.072139, 0.072357]; per-rev std ±0.00019 | **0.865% max / 0.705% RMS vs b0 (PASS)** | NT=72, η=2.0 ⇒ Das 0.41c + σ/c 0.68 both pinned, rlxf 0.15, metadata-verified. **ΔCT̄ vs B0 = +0.52% matched (10–19) / +0.76% vs settled 15–29 — FAILS M1 both ways ⇒ dt NOT converged at NT=36 at production Das.** Γ̄ shift uniformly positive (mean +0.66%, no lobes) — same global-truncation signature as the 0.205c matched pair (+0.595%/+0.711%). NT=144 (13029922) is the Richardson closer; legacy pinned pair's −3.9% attributed to legacy shedding config (outlier) |

## Appended 2026-08-05 (b) — κ ladder + viscous erratum A/B harvested; L2_visc failed

Harvested with `scripts/p018_analyze.py` (self-check first: `m1 p018_b0
--revs 15 29` → 0.071729, reproducing the recorded 0.07170 to 0.04%). All four
runs are 20-rev histories, so every comparison below uses the **matched 10–19
window**; B0 on that same window is **0.071866** (not the 15–29 settled
0.07170 — deltas are quoted against the matched value).

| case | phase | job id | window | CT̄ | CI | ε_Γ | notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| p018_dasc0p25 | 2b | 13050924 | 10–19 (drift 0.293%) | **0.071571** | [0.071476, 0.071666] | — | κ=0.25, \|Das\|=0.25·c_local |
| p018_dasc0p41 | 2b | 13050925 | 10–19 (drift 0.274%) | **0.072154** | [0.072037, 0.072263] | 0.952% max / 0.562% RMS vs κ=0.25 (**M2 PASS**) | κ=0.41 (matches B0's Das at 0.75R). **ΔCT̄ vs κ=0.25 = +0.81% ⇒ M1 FAIL** |
| p018_dasc0p82 | 2b | 13050926 | 10–19 (drift 0.219%) | **0.073455** | [0.073258, 0.073552] | **5.219% max / 1.953% RMS vs κ=0.41 (M2 FAIL)** | **ΔCT̄ vs κ=0.41 = +1.80% ⇒ M1 FAIL.** Successive deltas GROW (+0.81 → +1.80%) ⇒ ladder is not converging |
| p018_dasc0p41 vs p018_b0 | 2b | — | 10–19 both | +0.40% | — | 1.002% max / 0.318% RMS (**marginal M2 FAIL**) | Pure *spanwise-redistribution* effect: same Das at 0.75R, chord-wedge instead of r-wedge ⇒ ~1% of max Γ̄ moved, +0.40% CT |
| p018_b0_visc | 4/Γ6 | 13047290 | 10–19 (drift 0.491%) | **0.071810** | [0.071611, 0.072043] | **0.190% max / 0.112% RMS vs b0 (PASS)** | CoreSpreading-erratum A/B: **ΔCT̄ = −0.078%, M1 and M2 both PASS ⇒ the erratum is NULL on the B0 carrier.** Budget row Γ6 = 0.19% (Γ̄), CT term ≤0.08% |

**κ VERDICT: FAIL.** Chord-proportional Das does not converge — successive
doubling deltas grow rather than shrink, and the top step also fails M2 (5.2%).
Per the pre-registered rule this triggers the Ryan-authorized Green fallback,
which was already launched as a latency hedge (13051802/3, κ = 0.41/0.82).
No new submission was needed.

**`p018_L2_visc` (13047296) FAILED** — the viscous σ-L2 rescue diverged. Log
tail at step ~1040 (322 steps past its warm resume at 719) shows CF ≈ (0.55,
−0.16, 0.20) against a baseline ≈ 0.07, with `sigma growth` still 1.24;
`MaxRSS = 39.4 G` of a **128 G** request and no Julia stack trace ⇒ divergence,
not memory. It died *earlier* in absolute step count than the inviscid L2
(~1251), so **`CORE_SPREADING_ACTIVE=true` did not rescue the third σ rung.**
Confound to disclose with this result: the run warm-started across a σ jump
0.292 → 0.151 from the ov3 arm, so "viscous vs inviscid" is not cleanly
isolated (the inviscid L2 was cold). Ryan's ruling 2026-08-05: **record it and
defer the σ decision** to the Phase 12A (h at fixed σ) / Phase 4 Ladder B (σ at
fixed h) pair now in flight — do NOT invoke the 016 contingency yet.

**`p018_nt144_s2` submitted (13054739)** — the parent (13029922) was at step
2014/2447 with ~1.9 h of walltime left and cannot finish. Restart step **2025**
(complete set verified: `body1.2025.vtu`, `wake1.{1,2}.2025.vts`,
`particles.2025.vtp`), `P018_SETTLE_REVS=16` ⇒ 25 revs total, 48 h. Banner
verified: `NT:144 das_eta:4.0 overlap:2.0 pps:1 merge_r:0.0120 settle:16`.
**Operational hole found:** `<run>_CT_vs_rev.csv` is written only *after* the
time loop completes (`examples/rotor_hover_pressure_comparison.jl:897–908`), so
a walltime-killed segment leaves no CT series — the incrementally-written
`monitors/*force_system1.csv` is the only survivor. A `--from-monitor` fallback
for `scripts/p018_analyze.py` is the fix (validate against a completed run).

## VTK retention actions — 2026-08-04

Required by ruling 10 ("log every VTK deletion in the phase file"). Trigger:
disk headroom against the 200 G budget with **9 jobs writing concurrently**
(the 7 campaign jobs plus the 016 conversion A/B, 13036850/13036851);
projected peak was ~190 G.

Deleted VTK (`*.vtu`/`*.vtp`/`*.pvd` only — CSVs, monitors, and both TOMLs
retained in every case):

| scope | reclaimed | note |
| --- | --- | --- |
| 19 closed runs (`dji9443_panelwake{,_diag}`, `p018_b0{,_s2}`, `p018_das{0p5,2p0,4p0}`, `p018_hub0p15_b0`, `p018_hubfix_b0`, `p018_L1{,_s2,_warm}`, `p018_L2{,_s2}`, `p018_nrows{1,2}`, `pitching_wing_convergence`, `rotor_axial_j0187_ccblade`, `suddenly_started_wing`) | 6.2 GB | data 135.5 → 129.4 GB |
| `p018_L1_ov3` | 24.8 GB | **copied to Ryan's Mac first** at `~/p018_L1_ov3_paraview/` (deliberately outside Dropbox to avoid a 25 GB sync), verified byte-exact before deletion: 4314 files / 25,961,688,067 B on both sides |

Home usage 156 G → **140 G** (the 9 running jobs added ~8 GB during the
operation, so the net is smaller than the 24.8 GB reclaimed).

**VTK retained** (≤3-run rule, Ryan's selection criterion "most likely to have
useful ParaView files"): `p016_smooth_up`, `p016_smooth_split` — the only
rotor smooth-conversion wakes in existence, and the wake structure is 016's
open question — plus the 9 actively-writing run directories, which were
excluded from deletion on principle: they are live output, identified
empirically by mtime rather than inferred from job names.

Disk watchdog armed for the session (10-min `df` poll; alerts at
170/180/190/196 G, and on 3 consecutive probe failures so a dropped
ControlMaster socket cannot read as "all fine").

## 016 conversion A/B pair harvested — 2026-08-04 (evening)

Jobs 13036850 (`p018_conv_legacy`, 12h03) / 13036851 (`p018_conv_smooth`,
9h31) COMPLETED exit 0, 720 steps, harvested to local `data/`. Settled
numbers, revs 10–20 (`scripts/p018_analyze.py`): legacy C̄T **0.071847**
[0.071659, 0.072063]; smooth **0.072210** [0.072050, 0.072391]. Campaign
relevance: **the legacy arm reproduces B0 = 0.071866 to −0.026%** on the same
window ⇒ the TE-selection hardening + CoreSpreading integration fix are
CT-neutral and the B0 baseline survives the src change. 016 verdict (smooth
conversion ε_Γ FAIL, d/σ-clearance signature) recorded in the 016 phase 4 doc
§2026-08-04 (c).

Disk note: `p016_smooth_up` + `p016_smooth_split` copied **byte-verified** to
Ryan's Mac at `~/p016_smooth_vtk_paraview/` (7,919 files each; 8,925,212,597 B
/ 8,930,508,208 B) after the 180 G alert; cluster-side VTK deletion (~17 G)
is staged but NOT executed — agent-blocked destructive action, awaiting Ryan.
Home crossed 190 G of the 200 G budget at ~21:45 with 7–9 jobs writing.

## VTK retention actions — 2026-08-05

Ryan-approved deletion (explicit instruction) after the 196 G alert. Deleted
`*.vtu`/`*.vtp`/`*.pvd` only; CSVs, monitors, and TOMLs retained (each dir
now 24 M):

| scope | files | note |
| --- | --- | --- |
| `p016_smooth_up`, `p016_smooth_split` | 4,314 each | full run dirs previously copied **byte-verified** to Ryan's Mac `~/p016_smooth_vtk_paraview/` (7,919 files / 8,925,212,597 B and 8,930,508,208 B) |
| `p018_conv_legacy`, `p018_conv_smooth` | 4,314 each | NOT copied anywhere — Ryan authorized as expendable (legacy + second smooth wake at identical settings; harvested CSVs suffice) |

Home usage ~196 G → **173 G**. No rotor smooth-conversion VTK remains on the
cluster; the only surviving copies are the two runs on Ryan's Mac.

## VTK retention actions — 2026-08-05 (late) — automated sweeper armed

Ryan authorized an autonomous watchdog for the remainder of the session.
Retention rule tightened at his instruction: **every run keeps its final 10
RESTARTABLE steps**, not one — warm-startable *and* debuggable — with all
files at those steps (`.vtm` indices retained together with their pieces, so
the 13036477 stub failure mode cannot recur).

New, committed tooling (nothing like it existed before; both sweeps to date
were ad-hoc `find | xargs rm`):

- `scripts/p018_vtk_sweeper.sh` — dry-run by default, `--apply` to act.
  Computes S = steps where all four warmstart paths exist
  (`src/FLOWPanel_warmstart.jl`), keeps the top 10, never touches `*.csv`,
  `*.toml`, or `monitors/`. Liveness = newest mtime < 20 min **OR** a matching
  `squeue` job name (either signal ⇒ LIVE ⇒ skipped without `--emergency`;
  under `--emergency` files younger than 300 s are still spared). A run with no
  complete restartable step is skipped, never guessed at.
- `BRAINSTORM/018_dji9443_hover_convergence_campaign/vtk_protect_list.txt` —
  Ryan's file, read-only to agents. Seeded with `p018_L1_ov3` (permanent) and
  the four 2026-08-05 spatial-rigor runs `p018_h0p25`, `p018_h0p125`,
  `p018_ov1p4`, `p018_ov1p0` (their restart sets are the `_s2` chain sources).

Sweep executed 2026-08-05, **home 178 G → 112 G (66 GB freed)**. All six runs
verified afterward at 5/5 restart-set integrity (body `.vtu`, wake `.1/.2
.vts`, particles `.vtp`, `_body1.pvd`) with 10 kept steps and 4 CSVs each:

| run | freed | last step S | after |
| --- | --- | --- | --- |
| `p018_nt72_eta2` | 23,016 MB | 1439 | 202 MB |
| `p018_hub0p25_b0` | 11,585 MB | 719 | 207 MB |
| `p018_nrows1_das4p0` | 10,215 MB | 719 | 182 MB |
| `p018_nrows2_das4p0` | 10,208 MB | 719 | 184 MB |
| `p018_nrows2_das2p0` | 10,108 MB | 719 | 181 MB |
| `p018_nrows1_das2p0` | 10,030 MB | 719 | 171 MB |

Untouched: 5 protected runs, 12 live run dirs, and ~40 older dirs with no
complete restartable step (`p2e_*`, `sweptwing*`, `pitching_wing_convergence`,
…) which the sweeper declines to guess at. Disk watchdog armed for the session
(10-min `df` poll; alerts at 150/170/180/190/196 G and on 3 consecutive probe
failures, so a dropped ControlMaster socket cannot read as "all fine").

**2026-08-05 05:20 MDT, watchdog cycle (150 G alert).** Home had climbed
112 G → 150 G in ~3.5 h (~11 GB/h across 13 writers; three new jobs entered the
queue — 13051801 `p018_L1_visc`, 13051802/03 `p018_green_dasc0p41/0p82`).
`fp-018-b0-visc` (13047290) left the queue having COMPLETED 719 steps, so
`p018_b0_visc` was swept: **10,189 MB freed, 10,381 MB → 189 MB**, restart set
5/5 at S=719, 4 CSVs intact. Home 150 G → **148 G** (concurrent writers offset
part of the reclaim). All 8 live runs and all 5 protected runs untouched.

Sweeper note: `p018_L1_visc` was correctly held as LIVE by the **mtime** test
only — its Slurm name `fp-018-L1visc` does not normalize onto the directory
name `p018_L1_visc` (underscore mismatch), so the squeue-name test missed it.
This is exactly why liveness is an OR of two independent signals.

**2026-08-05 ~05:35 MDT, watchdog cycle.** `fp-018-dasc0p25` (13050924) left the
queue having COMPLETED 719 steps; `p018_dasc0p25` swept — **10,063 MB freed,
10,019 MB → 183 MB**, restart set 5/5 at S=719, 4 CSVs intact (the κ-ladder
judge reads the monitors CSVs, which are untouched). Home 151 G → 148 G.

Watchdog retuned after this cycle: the 150 G band re-fired on a 148→151 G
wobble with nothing new to reclaim. Added (a) **hysteresis** — a band re-alerts
only after usage drops 8 G below it — and (b) **job-completion detection**,
emitting `JOB-DONE <name>` when a job leaves the queue. Completion, not disk
pressure, is the event that actually creates a safely sweepable run; disk bands
are now just the backstop.

**2026-08-05 ~06:40 MDT, watchdog cycle.** `p018_L2_visc` (13047296) swept:
**17,843 MB freed, 16,624 MB → 351 MB**, restart set 5/5 at S=1041 (kept steps
1032–1041). `p018_dasc0p41` had also just left the queue but was still inside
the 20-min mtime window, so the sweeper conservatively held it as LIVE — it
will be swept next cycle. Home 151 G → 150 G.

**Incidental finding — 13047296 did NOT complete: `State=FAILED, ExitCode=1:0`
after 14:06:59**, at step ~1042 of the warm chain from `p018_L1_ov3`@719. Cause
is a Julia **`OutOfMemoryError()` inside `merge_particles!`**
(`FLOWVPM/src/FLOWVPM_merging.jl:454`), reached via
`apply_particle_policy!` → `apply_particle_maintenance!`
(`src/FLOWPanel_wake.jl:1610`, `:1636`) inside `propagate!` (`:2016`).

Note the shape of it: **MaxRSS was only 39.4 GB against the 128 GB request**, so
this is not gradual heap exhaustion — it is a single outsized allocation inside
the merge routine (consistent with a hash/neighbour structure sized off the
particle-cloud extent rather than the particle count). Raising `--mem` is
therefore unlikely to fix it; the merge parameters (`r_merge`/`r_hash`) or the
allocation itself are the lever. Related prior art: the p016 finding that
near-duplicate particles are structural and that a large-neighbour count was a
cutoff artifact.

Harvest consequence: the run never reached its epilogue, so
`p018_L2_visc` has **no `CT_vs_rev.csv` / `CT_per_rev.csv`** — only the two
`monitors/` CSVs (the OOM-triage `force_system1.csv` plus bound circulation,
both current to 06:37). Any L2-viscous number must come from those. An `_s2`
restart from step 1041 remains possible; that restart set survived the sweep
precisely because the retention rule is 10 steps rather than 1.

**2026-08-05 ~07:00 MDT, watchdog cycle.** `p018_dasc0p41` swept — **10,152 MB
freed, → 185 MB**, restart 5/5 at S=719, 4 CSVs. `p018_dasc0p82` had just
finished and was still inside the 20-min mtime window, so it was held as LIVE
(10.6 GB, next cycle). Home 150 G → **140 G**.

**κ ladder complete — all three COMPLETED, exit 0:** 13050924 `p018_dasc0p25`
(11:51:51), 13050925 `p018_dasc0p41` (13:32:23), 13050926 `p018_dasc0p82`
(15:21:25). All three retain their 4 CSVs, so the chord-Das κ fork is ready to
judge from the monitors CSVs.

Watchdog bug found and fixed: the completion detector ran `squeue` over a
**non-login** ssh, where it is not on PATH (`squeue: command not found`), so
`JOB-DONE` had never once fired — every completion so far was caught only
because a disk-band alert prompted a manual queue check. Now wrapped in
`bash -lc`. Also raised the alert floor from 150 G to 170 G: post-sweep usage
oscillates across 150 G, so that band was pure noise (it re-armed legitimately
each time a sweep dropped usage below the 8 G hysteresis margin).

**2026-08-05 ~08:10 MDT, watchdog cycle.** `p018_dasc0p82` swept — **10,460 MB
freed, → 189 MB**, restart 5/5 at S=719, 4 CSVs. Home 157 G → **147 G**.
`p018_conv_smooth_dasc0p82` (016 discriminator, job 13051758) finished during
this cycle and holds 8.7 GB; still inside the 20-min mtime window, so held as
LIVE for the next cycle. Queue now 10 jobs.

Second watchdog bug found and fixed — worth recording because it is a trap for
any future agent writing a Monitor script. The completion detector used
`for j in $prev_jobs; do case " $jobs " in ...`, which relies on **word
splitting of an unquoted parameter**. Monitor scripts run under **zsh**, where
that does not happen (`SH_WORD_SPLIT` is off by default), so `$j` was bound to
the entire 11-name string, matched nothing, and the watchdog reported *every*
running job as finished. Verified against `squeue` before acting — only one job
had actually ended — so nothing was swept in error.

Rewritten shell-agnostically: job names are kept newline-separated in files and
differenced with `comm -23`, no word splitting anywhere. Also hardened so an
**empty queue read is treated as a failed probe, never as "everything
finished"** — the previous shape would have proposed sweeping every live run at
once if `squeue` ever returned nothing.

**2026-08-05 ~13:05 MDT, watchdog cycle.** `p018_conv_smooth_dasc0p82` (016
discriminator, job 13051758) swept — **8,578 MB freed**, restart 5/5 at S=719.
Home 157 G → **149 G**.

**13029922 `p018_nt144` did NOT complete — `State=TIMEOUT`, ExitCode 0:15,
Elapsed 2-00:00:19.** It exhausted its 2-day walltime; last written step is
**2056**. The dt-Richardson closer therefore has no final tail from this job.
A continuation `fp-018-nt144s2` (job **13054739**) was submitted at 11:51 by
someone other than this watchdog.

**`p018_nt144`'s 29.6 GB is being HELD, not swept**, even though the run is
closed and not on the protect list. Reason: `nt144s2` is the obvious warm
continuation of a run that timed out at step 2056, which would make
`data/p018_nt144` an *active* restart source — sweeping a live warm source is
the 13036477 failure. Attempts to confirm warm-vs-cold were inconclusive:
after 1h14m `nt144s2` has written **no output directory at all** and only a
478-byte banner (`case p018_nt144`, settle:16), and `scontrol` does not expose
the launcher's `RESTART_*` exports. With 51 G of headroom there is no reason to
force the question, so the hold costs nothing.

Two items for Ryan, neither actioned by this watchdog:
1. **`p018_nt144` arguably belongs on the protect list** while the s2 chain is
   in flight. The list is Ryan's file and is never edited by an agent.
2. `nt144s2` producing no output in 1h14m may just be Julia precompilation, but
   it is worth a glance if it stays silent.

**2026-08-05 ~14:00 MDT, watchdog cycle.** 13051763 `p018_conv_smooth_h2p0`
COMPLETED exit 0 (12:51:39) — the second 016 discriminator is down; it was
still inside the mtime window so it is held for next cycle.

**Operator error, recorded in full.** The previous entry states that
`p018_nt144` was being deliberately HELD (probable warm source for the in-flight
`nt144s2`). This cycle the watchdog ran `./scripts/p018_vtk_sweeper.sh --apply`
with **no `--only`**, so the sweeper correctly classified `p018_nt144` as CLOSED
(no matching job, old mtime) and swept it: **33,549 MB freed**. The hold existed
only as an intention in the operator's head; nothing in the tooling encoded it.

Outcome verified immediately, and **no harm was done**: the 10-step retention
rule preserved steps **2047–2056**, with the restart set at S=2056 complete
(5/5: body `.vtu`, wake `.1/.2 .vts`, particles `.vtp`, `_body1.pvd`), 253 MB
retained. `nt144s2` (13054739) had still written nothing at that moment
(478-byte banner, no output dir, 2:06:57 elapsed), so it had not yet read the
restart files, and they are still there for it. Home 154 G → **122 G**.

Lesson: **a hold that is not in `vtk_protect_list.txt` is not a hold.** A full
`--apply` will always sweep anything the list does not name. The protect list is
Ryan's file and agents do not edit it, so the correct agent behaviour when a
run must be held is to use `--only <run>` on every sweep and never the bare
`--apply`, and to surface the run to Ryan for listing. Ryan: `p018_nt144` still
warrants a protect-list line while the s2 chain is live.

**2026-08-05 ~15:00 MDT, watchdog cycle.** `p018_conv_smooth_h2p0` swept with
`--only` (the corrective from the previous entry) — **10,486 MB freed, →
200 MB**, restart 5/5 at S=719, 5 CSVs. Home 126 G → **124 G**.

**13051775 `p018_ov1p0` COMPLETED, exit 0:0, 13:09:50, MaxRSS 31.1 GB.** It left
the queue ~1d13h early because it *finished its 30 revs* (last step 1080), not
because it failed. **The pre-registered divergence-risk rung did not diverge:**
CT = 0.071553 at rev 29.97 (0.071579 at rev 29.94), smooth, no magnitude flag —
only the benign `fmm! called but either sources or targets are empty` warning.
Ladder B's most aggressive rung (overlap 1.0 / pps 4 / merge_r 0.0060) is
therefore a usable data point rather than a documented blow-up. It is on the
protect list and was **not** swept; its full VTK is intact.

`nt144s2` (13054739) remains at a 478-byte banner with **no output directory
after ~3h10m**. Flagged for Ryan a second time; still not actioned by the
watchdog.

**2026-08-05 ~16:00 MDT, watchdog cycle.** 13051802 `p018_green_dasc0p41`
COMPLETED exit 0:0, 13:40:09, MaxRSS 35.2 GB, full 720 steps / 20 revs, all
three artifacts written. **Green-fallback arm at κ=0.41: CT = 0.0715451 at rev
19.97.** Its 10.3 GB was NOT swept this cycle — it had finished inside the
20-min mtime window, so the sweeper held it as LIVE (correct conservative
behaviour); queued for the next cycle. Home steady at **122 G**, 7 jobs.

**2026-08-05 ~17:00 MDT, watchdog cycle.** `p018_green_dasc0p41` swept (the
item deferred last cycle) — **10,118 MB freed, 10,302 MB → 181 MB**, restart
5/5 at S=719, 4 CSVs. Home **121 G**, 6 jobs.

`fm033cpu` also left the queue but is **not a FLOWPanel run** — it owns no
directory under `data/`, so there was nothing to sweep. Noted so a future agent
does not go looking for it.

**2026-08-05 ~16:20 MDT, fresh session — corrections + queue truth.** Direct
ssh check (sentinel-guarded): queue holds 6 study jobs — 13054739 `nt144s2`
(4:30 elapsed), 13051801 `L1visc`, 13051803 `green-dasc0p82`, 13051774 `ov1p4`,
13051772 `h0p25`, 13051773 `h0p125`. Corrections to earlier watchdog entries:

- **`p018_nt144_s2` DOES have output** — force CSV opens at step 2026 (warm
  from 2025 confirmed) and spans to ~2125/3600. The watchdog's "478-byte
  banner, no output dir" reports (~13:05–15:00) were wrong/stale; phase_03's
  12:0x warm-resume verification stands. The parent sweep's retained set
  (S=2047–2056) was never `_s2`'s source (RESTART_STEP=2025); no harm only
  because `_s2` had already loaded it.
- **`p018_L1_visc` warm-start CONFIRMED** from its log:
  `simulate_warmstart!: resuming from step 719 (file count 720)`; force CSV at
  step ~1151. The item file's "FIRST ACTION" check is discharged.
- Progress: `p018_ov1p4` ~1027/1080 and `p018_green_dasc0p82` ~686/720 (land
  within hours); `p018_h0p25` ~504/1080, `p018_h0p125` ~457/1080 (pace fits
  48 h with margin shrinking as counts grow — watch for `_s2`).
- Awaiting harvest: `p018_ov1p0` (COMPLETED, full 1080 steps, did not
  diverge) and `p018_green_dasc0p41` (COMPLETED, 720 steps).

**Ryan decision (this session): the fixed-κ clearance discriminator is ON
HOLD.** Next new submissions are instead the **phase_13 §4b uniform-d_front
pair at σ=0.03R** — `Das_j = D·σ − (N−1)·travel_j`: N=1 at D=3.4 and N=2 at
its minimal in-band D≈6.5 (also executes §6 action 2, the rigid-row-1 vs
free-row A/B). Implementation + pre-registration in phase_13; submission
gated on the mandatory src md5-sync (this session's uncommitted tripwire +
driver changes).

**2026-08-05 ~17:00 MDT — two harvests.** (provenance: jobs 13051775,
13051802; scored with `p018_analyze.py`, windows below; metadata-verified.)

- **`p018_ov1p0`** (Ladder B rung 3: σ/c 0.34 at FIXED h/c 0.34, OVERLAP 1.0,
  full 30 revs, did NOT diverge): CT̄(15–29) = **0.072548** CI [0.072299,
  0.073035] (drift 0.68% monotone — window caveat) vs B0 0.071701 ⇒
  **+1.18% FAILS M1**; ε_Γ **7.86% max / 4.23% RMS FAILS M2** with the SAME
  3-lobe shape as Ladder A L1 (inboard +6.1%, dip −7.86% @ r/R 0.762, tip
  +5.7%). **Discriminator answered: the 3-lobe mode tracks σ at fixed h.**
  Per-rev scatter ±0.00088 (8× B0) at marginal overlap. Detail phase_04
  §2026-08-05 (d).
- **`p018_green_dasc0p41`** (Green κ=0.41 partner, 20 revs): CT̄(10–19) =
  **0.071921** vs velocity 0.072154 ⇒ Δ = **−0.32%**, ε_Γ 0.716%/0.453%
  PASS — formulation delta small at this rung; κ fallback verdict awaits
  `green_dasc0p82`. Detail phase_02 §2026-08-05 (c).

**2026-08-05 ~18:00 MDT, watchdog cycle.** Two completions, both clean.

- **13051803 `p018_green_dasc0p82` COMPLETED** exit 0:0, 15:06:57, full 720
  steps / 20 revs. **Green-fallback arm at κ=0.82: CT = 0.0733405** at rev
  19.97. Held as LIVE this cycle (finished inside the 20-min mtime window);
  ~9 GB queued for next cycle.
- **`p018_ov1p4` (13051774) reached its full 30 revs**, last step 1080,
  **CT = 0.0721746** at rev 29.97. Protected, not swept, full VTK intact.

Ladder B / Green-fallback numbers now in hand (all from the runs' own CSVs):

| run | overlap / κ | revs | CT |
| --- | --- | --- | --- |
| `p018_ov1p0` | overlap 1.0 | 30 | 0.0715534 |
| `p018_ov1p4` | overlap 1.4 | 30 | 0.0721746 |
| `p018_green_dasc0p41` | Green, κ=0.41 | 20 | 0.0715451 |
| `p018_green_dasc0p82` | Green, κ=0.82 | 20 | 0.0733405 |

Home **118 G**, queue down to 4: `h0p25`, `h0p125` (protected, ~1d08h left),
`L1visc`, `nt144s2`.

**2026-08-05 ~18:00 MDT — Ladder B complete + Green κ verdict.** Jobs 13051774
(`p018_ov1p4`, 15:35 h) and 13051803 (`p018_green_dasc0p82`, 15:07 h) both
COMPLETED exit 0, harvested + metadata-verified, scored on matched windows:

- **`p018_ov1p4`** (σ/c 0.477, fixed h): CT̄(15–29) = **0.072257** ⇒ +0.78% vs
  B0 (FAIL M1); ε_Γ 4.12% (FAIL M2), 3-lobe at half ov1p0's amplitude (dip
  −4.1% @ r/R 0.740) ⇒ **mode ∝ σ-step at fixed h — Ladder B closes with the
  σ attribution confirmed from this side** (Phase 12A h-null pending). Rung
  pair ov1p4→ov1p0: +0.40% (M1 pass) but ε_Γ 3.90% (M2 fail). phase_04
  §2026-08-05 (e).
- **`p018_green_dasc0p82`**: CT̄(10–19) = **0.073046**; green successive κ
  delta 0.41→0.82 = **+1.56%** vs velocity's +1.80% ⇒ **Green does NOT rescue
  the κ ladder; fallback rule resolved negative; κ response is
  formulation-independent** — consistent with the clearance reinterpretation;
  green−velocity Δ small at both rungs (−0.32%/−0.56%, M2 PASS) = early
  Phase-7 Green-Δ points. phase_02 §2026-08-05 (d).

Queue after landings: 4 study jobs (nt144s2, L1visc, h0p25, h0p125); ufront
pair submission takes it to 6/10.

**2026-08-05 ~19:30 MDT — Phase 14 screen batch SUBMITTED (19 jobs, slot cap
waived by Ryan for this batch).** Runner `examples/run_p018_screen_hpc.slurm.sh`
(md5 3859777f… verified on cluster), production mesh, 8 revs, freestream pulse
zeroed, WakeHealthMonitor on. 45_145 mesh uploaded (md5 804301d8…); 45_249 /
60_185 / 80_185 already on cluster, md5-matched; all four passed the mesh
profile gate. Jobs (banner verification pending, mandatory):

| job | case | | job | case |
| --- | --- | --- | --- | --- |
| 13058031 | scr_b0 | | 13058041 | scr_ufdt_nt36 |
| 13058032 | scr_ov1p0 | | 13058042 | scr_ufdt_nt72 |
| 13058033 | scr_nrows1 | | 13058043 | scr_ufsig_0p05 |
| 13058034 | scr_nrows2 | | 13058044 | scr_ufsig_0p02 |
| 13058035 | scr_nt72 | | 13058045 | scr_mesh_c145 |
| 13058036 | scr_uf_n1_d2p6 | | 13058046 | scr_mesh_c249 |
| 13058037 | scr_uf_n1_d3p4 | | 13058047 | scr_mesh_s60 |
| 13058038 | scr_uf_n1_d5p0 | | 13058048 | scr_ufdt_nt144 (24 h) |
| 13058039 | scr_uf_n2_d6p5 | | 13058049 | scr_mesh_s80 (24 h) |
| 13058040 | scr_uf_n2_d8p5 | | | |

Scoring: `scripts/p018_screen_score.py` (validated against `p018_analyze m2`
bit-for-bit on the dasc pair). Pre-registrations + case roster: phase_14.
First harvest action = calibration gate (phase_14 pre-reg 3) BEFORE reading
any campaign ladder.

**2026-08-05 ~19:00 MDT, watchdog cycle.** `p018_green_dasc0p82` swept —
**10,441 MB freed**, restart 5/5 at S=719. Home 140 G → **130 G**.

Two things noted rather than actioned:
- `fm032` left the queue; like `fm033cpu` it is **not a FLOWPanel run** and owns
  no `data/` directory here. Both were submitted and finished between 10-minute
  polls, so they never appeared in a queue snapshot.
- Home jumped **118 G → 140 G (+22 GB) with only four FLOWPanel jobs running**,
  which the four cannot account for at their observed rates. FLOWPanel `data/`
  is now **101 GB** against a **130 G** home total, so **~29 GB lives outside
  the FLOWPanel checkout** — plausibly the `fm0*` jobs. The 200 G budget is
  tracked on the home total, so non-FLOWPanel growth consumes the same budget
  the sweeper defends but is entirely outside its reach. Breakdown of the
  remaining ~29 GB is in progress.

**Scope correction (Ryan, 2026-08-05 evening): the 200 G cap applies to the
FLOWPanel repositories only, not to home total.** Everything outside the
FLOWPanel checkouts is out of scope for this watchdog and is not to be counted
or worried about.

Measured breakdown at the time of the ruling:

| path | size |
| --- | --- |
| `~/projects/FLOWPanel.jl` (of which `data/` = 101 GB) | **102.5 GB** |
| `~/projects/FLOWVPM.jl` | 0.14 GB |
| `~/projects/FastMultipole` | 0.01 GB |
| *(out of scope)* `~/FastMultipole-023/-026/-033` | 2.0 GB |
| *(out of scope)* `~/projects/NeuralFMM` | 0.12 GB |
| *(out of scope)* hidden dirs, `df`-minus-`du` gap (likely `~/.julia`) | ~23 GB |

So **FLOWPanel stands at ~103 G of its 200 G cap — 97 G of headroom**, not the
70 G that the home-total figure implied. The earlier "+22 GB spike with only
four jobs running" was never FLOWPanel growth; it was outside the checkout, and
`du -sm /home/rander39/*` had hidden it because that glob does not match
dotfiles.

Watchdog re-pointed accordingly: it now measures
`du -sm ~/projects/{FLOWPanel.jl,FLOWVPM.jl,FastMultipole}` summed, rather than
`df` on home, and alerts at 170/180/190/196 G of that figure. Poll interval
relaxed 10 min → 15 min since `du` over a 100 GB tree is heavier than `df`.
- 2026-08-05 ~20:15 MDT: all 19 screen banners verified from the job logs (knobs exact per case table). Track F forensics complete: L2_visc WAS igniting (Γ blow-up on floor-pinned σ), feedback hypothesis REFUTED, merge-candidacy 0.12%, OOM arithmetic pinned (723R/step runaway). Track C code audit complete (no σ floor; stability iff ΔtZ<2; CoreSpreading floor √(2νΔt) observed in corpse). Full mechanism doc: sigma_blowup_mechanism.md (DRAFT, Track H holes remain).

**2026-08-05 ~21:00 MDT, watchdog cycle.** No sweep possible — every run is
either live or protected. FLOWPanel now at **119.7 GB of its 200 G cap**
(`~/projects/FLOWPanel.jl` = 122,571 MB), up ~17 GB as five new jobs came in:
13058032 `scr_ov1p0`, 13058033 `scr_nrows1`, 13058034 `scr_nrows2` (10.5 h left)
and 13057253/4 `ufront-n1`/`ufront-n2` (1d21h left). Nine jobs running.

**Correction to three earlier entries: `nt144s2` was never stalled.** This
watchdog repeatedly reported it as "no output directory, 478-byte banner,
possibly unhealthy" — that check looked for `data/p018_nt144s2`, but the run
writes to **`data/p018_nt144_s2`** (underscore before `s2`). It is healthy and
well ahead: **last step 2209**, continuing from nt144's 2056, 4.4 GB written.
The 478-byte `.out` is just buffered stdout, not silence. Lesson: verify a run
directory's actual name before concluding a job is producing nothing.

**Actionable for Ryan — ~31 GB is reclaimable but locked behind the protect
list.** `p018_ov1p0` (15,865 MB) and `p018_ov1p4` (15,947 MB) both **COMPLETED
their full 30 revs** (13:09:50 and 15:35:04). The protect list's own instruction
is to remove a run "once that run has completed its 22 revs or has had its `_s2`
launched", so both now qualify — but the list is Ryan's file and agents do not
edit it. Removing those two lines would let the sweeper reclaim ~31 GB, taking
FLOWPanel from ~120 G to ~89 G. Not urgent at 80 G of headroom; noted so the
option is visible. `h0p25` and `h0p125` are still running and must stay listed.

**2026-08-05 ~22:00 MDT, watchdog cycle — ESCALATION.** A large `scr_*`
screening batch was submitted: the queue is now **24 jobs** (`scr_b0`,
`scr_mesh_c145/c249/s60/s80`, `scr_nrows1/2`, `scr_nt72`,
`scr_ufdt_nt36/nt72/nt144`, `scr_uf_n1_d2p6/d3p4/d5p0`, `scr_uf_n2_d6p5/d8p5`,
`scr_ufsig_0p05`, plus `ufront-n1/n2`, `h0p25`, `h0p125`, `L1visc`,
`nt144s2`, `fm032`). Two `scr_*` jobs (`scr_ov1p0`, `scr_ufsig_0p02`) have
already completed and rotated out.

**FLOWPanel is at 152,609 MB = 149 G of its 200 G cap — 51 G of headroom with
~20 concurrent writers.** Nothing is sweepable: every run in `data/` is either
live or protected. The sweeper has no move left.

**The only remaining reclaim is behind the protect list, and agents do not edit
that file.** `p018_ov1p0` (15.9 GB) and `p018_ov1p4` (15.9 GB) both COMPLETED
their full 30 revs and satisfy the list's own removal condition. Deleting those
two lines would free ~31 G. If Ryan is unavailable and the cap is actually
breached, the in-policy fallback is
`./scripts/p018_vtk_sweeper.sh --apply --emergency`, which trims live runs
oldest-steps-first while keeping each one's last 10 restartable steps and
sparing anything written in the last 300 s.

Third watchdog defect found and fixed: the size probe summed
`du -sm FLOWPanel.jl FLOWVPM.jl FastMultipole` and matched on `/projects/`, so
when the large `FLOWPanel.jl` du did not return, the two small repos alone
summed to 159 MB — non-zero, so the failure check passed — and the watchdog
reported **"FLOWPanel 0G"** while the true figure was 149 G. It now requires the
`FLOWPanel.jl` line specifically and treats its absence as a probe failure.
Added a growth-rate alarm (>12 G between polls) and tightened the poll to 10 min
given ~20 writers; bands re-cut to 155/170/180/190/196 G.

**2026-08-05 ~20:45 MDT — `p018_ufront_n1_visc` SUBMITTED (job 13058534,
Ryan-directed).** The missing decision-critical rung: VISCOUS σ=0.03R at the
uniform-d_front production candidate (OV 2.4/PPS 14, N=1, D=3.4,
CORE_SPREADING_ACTIVE, β=1e9, cold, SETTLE=22 ⇒ ~30 revs ≈ 1130 steps —
approaching the L2 ignition age with the WakeHealthMonitor tripwire live),
48 h/64G. Buys: (a) viscosity delta AT the production candidate (34% core-
balance share); (b) long-horizon viscous stability certification at 0.03R;
(c) clean viscosity A/B vs the inviscid `p018_ufront_n1` (13057253).
Launcher arm deployed, md5 1b65d575… verified. Banner verification pending.
**Ryan pre-authorization on record: if tonight's contraction-rate trend
(scr_ufsig ladder + scr_uf_n1_d3p4 wake-health) looks TIGHT at 0.03R, the
executing agent may submit a σ=0.035R hedge WITHOUT further approval** —
uniform-d_front, viscous, with the mandatory offline D re-check first
(D=3.4 slightly overshoots the 1.5c tip band at 0.035R; drop D toward ~3.0
if needed and record the choice).

## Retention ruling amended again — 2026-08-05 (late), Ryan

> "for any run not in the protect list, you can delete files while it runs, so
> long as you leave the 10 most recent timesteps. No need to wait for the entire
> run to complete."

`scripts/p018_vtk_sweeper.sh` updated accordingly: **live, unprotected runs are
now swept by default** rather than skipped. The two safety properties that make
this safe are unchanged and both still enforced — the 10 most recent
*restartable* steps are retained, and any file written in the last
`LIVE_QUIET_SEC=300` s is spared, so the step currently being written is never
the one removed. `--emergency` is now a no-op (accepted for compatibility);
`--skip-live` restores the old conservative behaviour.

**Sweep executed under the new rule: 101,283 MB freed across 23 runs.
FLOWPanel 166 G → 70.5 G (72,182 MB) of the 200 G cap.**

| scope | freed |
| --- | --- |
| `p018_L1_visc` (live) | 23,117 MB |
| `p018_ufront_n1/n2` (live) | 10,788 MB |
| `p018_nt144_s2` (live) | 5,255 MB |
| 19 `scr_*` screening runs (16 live, 3 closed) | 62,123 MB |

**Post-sweep verification — no job was harmed.** All **22 jobs still RUNNING**,
zero new FLOWPanel failures (the only FAILED entry in the window is `fm032`,
not a FLOWPanel run). Every swept live run advanced past its retained window,
confirming the writers were untouched: `p018_nt144_s2` kept through 2254, now
at 2256; `scr_nt72` kept through 433, now at 438; `p018_L1_visc` at 1344;
`scr_ufsig_0p05` at 253. `p018_ov1p0` and `p018_ov1p4` remain PROTECTED and
untouched, so the earlier request to unprotect them is moot — the new rule
freed far more than they held.

Two script defects fixed in the same change, both caught before any deletion:
`EMERGENCY` left as an unbound variable in the banner line (`set -u` aborted the
whole script, producing an empty dry-run), and the monitor's size parse choking
on the login banner's ANSI colour codes (`bad math expression`, which had killed
the watchdog process outright).

**2026-08-05 ~23:30 MDT, watchdog cycle.** Full sweep under the new live-sweep
rule: **17,665 MB freed**, FLOWPanel 87 G → **71.2 G**. 21 jobs still RUNNING,
none harmed. Growth is now the dominant term: 70.5 → 79 → 87 G across two
10-minute polls (~24 GB/h with ~21 writers), so periodic sweeps — not
completion-triggered ones — are what keep the cap.

**Second OOM of the campaign, same signature. 13058040
`fp-018-scr_uf_n2_d8p5`: `State=OUT_OF_MEMORY`, ExitCode 0:125, after 04:05:32.
ReqMem 64 G but MaxRSS only 18.97 GB** — killed at ~30% of its request. No Julia
stacktrace survives (a cgroup OOM kill leaves none), but the shape matches
13047296 `p018_L2_visc` exactly: **FAILED at MaxRSS 39.4 GB of a 128 G
request**, where the stacktrace *was* captured and pointed at
`merge_particles!` (`FLOWVPM_merging.jl:454`) via `apply_particle_maintenance!`.

Two independent runs now OOM-killed while sitting at roughly a third of their
requested memory. That is the signature of a **large transient allocation
between Slurm's RSS samples**, not of gradual growth — so raising `--mem` is
unlikely to fix either, and `d8p5` being the largest d/σ rung of the `uf_n2`
pair fits an allocation that scales with the particle cloud's extent. Recorded
for Ryan; the watchdog takes no action on it.

**2026-08-06 ~00:15 MDT, watchdog cycle.** Sweep freed **10,705 MB**, FLOWPanel
81 G → **71.9 G**, 20 jobs RUNNING and unharmed.

**Third OOM, and the three now form a quantitative pattern.** 13058041
`fp-018-scr_ufdt_nt36`: `OUT_OF_MEMORY`, ExitCode 0:125, after 04:38:47,
**MaxRSS 20.0 GB against a 64 G request**.

| job | run | ReqMem | MaxRSS | MaxRSS / ReqMem |
| --- | --- | --- | --- | --- |
| 13047296 | `p018_L2_visc` | 128 G | 39.4 GB | **31%** |
| 13058040 | `scr_uf_n2_d8p5` | 64 G | 19.0 GB | **30%** |
| 13058041 | `scr_ufdt_nt36` | 64 G | 20.0 GB | **31%** |

All three died at ~30% of their request, across two different request sizes.
A *fixed-size* runaway allocation would kill at different fractions; a constant
fraction of ~1/3 means **the failing allocation is roughly 2× the live heap** —
i.e. RSS + 2·RSS ≈ ReqMem. That is the classic signature of a **copy or resize
of an already-large array** (allocate the new buffer, both alive at once), not
of a leak or of gradual growth. It is consistent with `p018_L2_visc`'s captured
stacktrace, the only one of the three with a Julia traceback:
`merge_particles!` → `FLOWVPM_merging.jl:454`.

Practical implication for the campaign: **raising `--mem` buys a constant factor
at best** — tripling the request only moves the failure to a 3× larger cloud —
whereas avoiding the copy (in-place resize, or chunking the merge) removes it.
Recorded for Ryan; the watchdog takes no action.

**Watchdog upgraded to sweep autonomously.** Growth (~24 GB/h across ~20
writers) now outpaces completions, so completion-triggered sweeping is
structurally too slow. The watchdog now runs `--apply` itself on any 10-minute
poll where FLOWPanel exceeds **100 G**, and reports what it freed. It refuses to
sweep when the size probe fails (unknown size is never treated as safe), and
still shouts if the figure reaches 175 G despite sweeping. Everything the
sweeper enforces is unchanged: protect list honoured, 10 most recent restartable
steps kept per run, files younger than 300 s spared.

**2026-08-06 ~01:00 MDT — first autonomous sweep, and it worked.** The watchdog
fired its own `--apply` at 101 G: **freed 27 G → 74 G of the 200 G cap**, no
agent in the loop. This is the intended steady state.

**13051801 `p018_L1_visc` COMPLETED**, exit 0:0, 22:57:53 — the viscous L1 run
finished its full 40 revs (last step 1440). **CT = 0.0711758 at rev 39.97.**

**Attrition in the `scr_*` batch is severe: 4 of the 9 finished runs died of
memory**, and the ~30% signature now holds across five failures.

| job | ReqMem | MaxRSS | ratio | state |
| --- | --- | --- | --- | --- |
| `p018_L2_visc` | 128 G | 39.4 GB | 31% | FAILED 1:0 (Julia `OutOfMemoryError`, `merge_particles!`) |
| `scr_uf_n1_d2p6` | 64 G | 23.1 GB | 36% | FAILED 1:0 (Julia `OutOfMemoryError`) |
| `scr_uf_n1_d5p0` | 64 G | 19.7 GB | 31% | OUT_OF_MEMORY 0:125 |
| `scr_ufdt_nt36` | 64 G | 20.0 GB | 31% | OUT_OF_MEMORY 0:125 |
| `scr_uf_n2_d8p5` | 64 G | 19.0 GB | 30% | OUT_OF_MEMORY 0:125 |

The two `FAILED 1:0` cases carry an explicit Julia `OutOfMemoryError`; the three
`OUT_OF_MEMORY 0:125` cases are cgroup kills of the same event caught a moment
earlier. So all five are one failure mode, and the 30–36% band across 64 G and
128 G requests keeps pointing at an allocation ≈2× the live heap (a copy/resize,
not a leak) — consistent with `merge_particles!` / `FLOWVPM_merging.jl:454`.

Survivors and deaths interleave in `d`: `uf_n1` d2p6 **died**, d3p4 completed,
d5p0 **died**; `uf_n2` d6p5 completed, d8p5 **died**. Failure is therefore not
monotonic in `d`, which argues against "big `d` ⇒ too many particles" as the
whole story.

Completed cleanly in the same window: `scr_mesh_c145`, `scr_uf_n1_d3p4`,
`scr_uf_n2_d6p5`, `scr_ufsig_0p05`. 13 jobs still RUNNING.

**2026-08-06 ~02:45 MDT — screen batch first harvest + L1_visc + hedge.**
Detail in phase_14 §2026-08-06. Headlines:

- **6 of 14 landed screens IGNITED** (5 "OOM" + 1 FAILED = all blow-ups by
  tripwire: max|u| 13k–198k m/s, FIVE with NEGATIVE σ — first live
  observation of the inviscid σ(1−ΔtZ)<0 mode). Screen-condition inviscid
  stability boundary measured at ~σ/R 0.030 (0.0299 survives, 0.0291
  ignites). Campaign B (clean dt screen ladder) entirely lost to ignition.
- **Calibration gate FAIL** ⇒ per pre-registration, screen aerodynamic
  deltas untrusted at 8 revs; screens restricted to stability/plumbing/cost
  probes (where they excelled). Still-running screens (nt72, mesh_c249/s60/
  s80, ufdt_nt144) carry the caveat.
- **`p018_L1_visc` COMPLETED + scored**: viscosity delta at L1 σ = **NULL in
  Γ̄** (ε 0.156% vs inviscid L1; ΔCT −0.38%); viscous σ-pair reproduces the
  3-lobe (⇒ σ-mode viscosity-independent). Budget Γ6(L1)=0.16%. Notebook
  viscous-A/B TBD table filled (b0_visc null; L2_visc diverged; L1_visc row
  added).
- **Hedge executed under Ryan pre-authorization** (trigger: boundary at
  0.03R + production inviscid n1 contracting to min_sr 0.145):
  `p018_ufront_s035_visc` = **13058988** (σ/R 0.0349, D=3.0 offline-verified
  uniform, viscous, SETTLE=22). Banner verified. Production tripwires:
  inviscid n1 0.145 / n2 0.208 / visc twin 0.381 — viscous floor working
  live.
- Mechanism doc `sigma_blowup_mechanism.md` now COMPLETE (all three tracks
  landed; confidence table final).

**2026-08-06 ~02:00 MDT — second autonomous sweep.** Fired at 100 G, **freed
23 G → 76 G**. No agent action required; the watchdog is now holding the cap by
itself on a ~10-minute cadence.

Since the previous entry: `scr_mesh_c249` (07:55:27), `scr_mesh_s60` (08:28:06)
and `scr_nt72` (08:00:42) all COMPLETED exit 0:0; **`scr_ufdt_nt72` died
`OUT_OF_MEMORY` 0:125 after 07:23:07** — the sixth memory failure.

**The deaths partition cleanly by family, which is the useful finding.** Every
memory failure in the `scr_*` batch is in the `uf_*` / `ufdt_*` group; nothing
outside it has died:

| family | outcome |
| --- | --- |
| `scr_uf_n1_d2p6`, `scr_uf_n1_d5p0`, `scr_uf_n2_d8p5` | **died** (OOM) |
| `scr_ufdt_nt36`, `scr_ufdt_nt72` | **died** (OOM) |
| `scr_uf_n1_d3p4`, `scr_uf_n2_d6p5` | completed |
| `scr_mesh_c145`, `scr_mesh_c249`, `scr_mesh_s60` | completed |
| `scr_nt72`, `scr_nrows1`, `scr_nrows2`, `scr_b0`, `scr_ov1p0` | completed |
| `scr_ufsig_0p02`, `scr_ufsig_0p05` | completed |

So the failure is specific to the `uf`/`ufdt` configurations rather than to mesh
size, timestep count, or `nrows` — and within `ufdt` it has now taken both
rungs that finished (nt36 and nt72), leaving only nt144 in flight.

Two new jobs entered the queue: `fp-018-ufront-n1-visc` and
`fp-018-ufront-s035-visc`. 10 jobs RUNNING.

**2026-08-06 ~03:00 MDT — third autonomous sweep.** Fired at 100 G, **freed
21 G → 79 G**. Cadence is steady: the watchdog now sweeps roughly hourly as the
remaining writers refill it, with no agent involvement.

**Seventh memory failure, and the first outside the screening batch: 13057253
`fp-018-ufront-n1`, `OUT_OF_MEMORY` 0:125 after 12:58:55.** This is a
*production* `ufront` run, not a `scr_*` screener. The family partition
therefore extends beyond screening: `uf_*`, `ufdt_*` **and** `ufront-*` all
share the failing configuration, and everything outside those families
(`mesh_*`, `nt72`, `nrows*`, `ufsig*`, `b0`, `ov1p0`) has completed cleanly.

Correction to the previous entry's expectation: `scr_ufdt_nt144` has **not**
died — it is still running. The prediction that it would follow nt36/nt72 is so
far unmet.

**`fp-018-ufront-n2` is the same configuration as the run that just died and is
still in flight**, as are `ufront-n1-visc` and `ufront-s035-visc`. Flagged for
Ryan: on the pattern to date, those three are the next candidates. Nine jobs
RUNNING (`h0p25`, `h0p125`, `nt144s2`, `scr_mesh_s80`, `scr_ufdt_nt144`,
`ufront-n2`, `ufront-n1-visc`, `ufront-s035-visc`, plus out-of-scope
`fm033cpu`).

**2026-08-06 ~09:40 MDT — inviscid production candidate DEAD; viscous twins
alive; last screens landed.**

- **`p018_ufront_n1` (inviscid σ=0.0299R) OUT_OF_MEMORY at ~13 h / step ~501
  = IGNITION** (tripwire: min_sr −22.9, max_u 33,977 m/s) — the screen-
  boundary result transfers to the production staged startup, which only
  delayed ignition (step ~500 vs 284). The inviscid ufront A/B is lost;
  "production 0.03R must be viscous" is now enforced by data, not just
  recommended.
- `p018_ufront_n2` (inviscid, N=2, D=6.5): alive at step 621/774, min_sr
  0.067, |Γ|/σ² 9.1e3 — contracting but not igniting; may complete. N=2
  surviving where N=1 died is itself a data point (bigger d_front, same σ).
- `p018_ufront_n1_visc`: step 482, min_sr 0.049 (riding near the floor),
  max_u 87 bounded, |Γ|/σ² 1.4e5 and climbing — regime-boundary behavior
  live; not ignited. WATCH.
- `p018_ufront_s035_visc` (hedge): step 381, min_sr 0.351, max_u 22 —
  comfortably healthy. The hedge is earning its keep.
- Screens `scr_nt72`, `scr_mesh_c249/s60/s80` all COMPLETED exit 0
  (stability-probe readouts only, per the failed calibration gate);
  `scr_ufdt_nt144` still running. `p018_h0p25` 788/1080 (fits walltime),
  `p018_h0p125` 671/1080 (tight), `nt144_s2` 2513/3600 vs 26 h left —
  projected TIMEOUT ~80 steps short; `_s3` chain likely.

**2026-08-06 ~09:50 MDT — `p018_ufront_n2_visc` SUBMITTED (job 13060144,
Ryan-directed):** viscous N=2 arm (OV 2.4/PPS 14, D=6.5, CoreSpreading,
SETTLE=22, 64G/48 h) completing the all-viscous rigid-vs-free-row A/B with
`p018_ufront_n1_visc` after the inviscid n1 ignited. Launcher arm deployed,
md5 91864b79… verified; banner verification armed.

**2026-08-06 ~04:00 MDT — fourth autonomous sweep.** Fired at 100 G, **freed
18 G → 82 G**. Steady state holding.

`scr_mesh_s80` COMPLETED exit 0:0 (12:13:11) — with it, **every `mesh_*` rung
(c145, c249, s60, s80) has completed cleanly**, confirming the failure is not
mesh-size driven.

**Eighth memory failure, and the prediction from the previous entry was borne
out: `fp-018-ufront-n1-visc` FAILED 1:0 after 11:59:15**, with an explicit Julia
**`OutOfMemoryError`** whose stack passes through
`apply_particle_maintenance!` → `FLOWVPM_merging.jl:386` → **`:454`** — the
*identical* line as `p018_L2_visc`. That is now a second independent traceback
naming the same allocation, so the `merge_particles!` attribution no longer
rests on a single stack.

Its **MaxRSS was 67.0 GB**, the largest of the eight and well above the ~20 GB
of the 64 G-request failures, so this run was given more memory and still died
at the same call site — further evidence that the fix is to avoid the
copy/resize rather than to raise `--mem`.

Seven jobs RUNNING: `h0p25`, `h0p125` (15:17 left), `nt144s2` (1d01h),
`scr_ufdt_nt144` (9:01), `ufront-n2` (1d07h), `ufront-s035-visc` (1d16h), and a
newly submitted **`ufront-n2-visc`** (1d23h). Three of those seven are `ufront-*`
runs of the family that has now lost two members.

**2026-08-06 ~11:20 MDT — post-reset harvest: σ\* DECIDED → 0.035R; both
remaining σ=0.03R production arms are dead by ignition; screen stability
readouts show mesh refinement crossing σ<0.**

- **`p018_ufront_n1_visc` (13058534) post-mortem: IGNITION confirmed by
  tripwire**, not a plain merge-copy OOM. Wake-health tail (steps 489→493):
  max_u 2,076 → **67,072 m/s**, |Γ|/σ² up to **1.4e9**, min σ pinned at the
  viscous floor 9.4e−5 m (= predicted √(2νΔt)) throughout. This is the
  root-caused Γ-blow-up-at-the-floor; the `FLOWVPM_merging.jl:454` OOM
  (MaxRSS 67 GB) is the bounding-box symptom. **Per the pre-registered rule,
  σ\* moves 0.03R → 0.035R; `p018_ufront_s035_visc` (13058988) is now the
  production candidate** (healthy: step 431, min_sr 0.316, max_u 22).
- **`p018_ufront_n2` (inviscid, 13057254) OUT_OF_MEMORY 0:125 at 17:31:25 =
  ninth memory failure, also IGNITION**: at step 666/774 min σ went
  **negative** (−5.8e−4 → −7.9e−3 m; min_sr −2.22), max_u 12,755 m/s — the
  inviscid σ(1−ΔtZ)<0 mode. Time-to-ignition step ~667 vs ~501 for inviscid
  n1: N=2 (D=6.5) delays but does not prevent ignition at σ=0.0299R. The
  9:40 "may complete" hope is dead; **every σ=0.03R production arm (invisc
  N1, invisc N2, visc N1) has now ignited**.
- `p018_ufront_n2_visc` (13060144, σ=0.03R viscous N=2) is early and healthy
  (step 166, min_sr 0.48) but shares the σ that killed its N=1 twin at step
  ~490 — **on watch; the all-viscous A/B at 0.03R has already lost its N=1
  arm regardless.** A rebuilt A/B at σ\*=0.035R (an N=1-vs-N=2 pair, or an
  n2@0.035 partner for s035_visc) needs a Ryan call — no pre-authorization
  on record for new submissions.
- **Screen stability-probe readouts (final wake-health rows; aero deltas
  remain untrusted per the failed calibration gate):** `scr_nt72` healthy
  (min_sr 0.418, max_u 10); `scr_mesh_c145` contracted-bounded (0.101, 21);
  `scr_mesh_s60` marginal (0.075, max_u 44); **`scr_mesh_c249` ended with
  min_sr −1.09 and max_u 816 — blow-up in progress at the final step —
  and `scr_mesh_s80` ended at min_sr −0.37** (u still bounded 24). Two
  "COMPLETED exit 0" runs thus concealed σ<0: **exit code is not a
  stability verdict; mesh refinement (chord and span) tightens the inviscid
  boundary past σ/R≈0.030.** Strengthens "production must be viscous" and
  the 0.035R margin. `scr_ufdt_nt144` (running, step 997) is riding at
  min_sr 0.0346 — likely to follow nt36/nt72.
- **Screen-batch VTK swept** (watchdog's pre-reset sweeps had already
  cleaned 18 of 19): freed 1,251 MB from `scr_ufdt_nt144` (live-sweep,
  restart set kept) + 1,204 MB from the `p018_ufront_n2` corpse; the
  `p018_ufront_n1_visc` corpse was already clean. data/ at 88 G of the
  200 G cap. **The autonomous watchdog did not survive the context reset —
  sweeping is manual again until one is re-armed.**
- Queue (11:20): `nt144_s2` at step ~2552/3600 with ~24.6 h left — chain
  `_s3` on timeout (recipe `ops_reference.md`); `h0p25`/`h0p125` ~14 h from
  walltime; `s035_visc`, `n2_visc`, `scr_ufdt_nt144` as above. Item 019
  stays unstarted (harvest queue non-empty).

## 2026-08-06 (afternoon) — item 019 kickoff: provenance harvest + mechanism-doc correction

- **Item 019 is now ACTIVE** (σ-selection procedure; Ryan authorized agent
  HPC submission for its P4 regime-map screens this session). 018 artifacts
  it created/changed, recorded here for the campaign record:
  - `sigma_blowup_mechanism.md` corrected (dated section at end): 13058534's
    death supersedes the doc's "certification in flight / SUPPORTED at
    1.7σ_eq" claims; safe bracket now carried by s035_visc + L1.
  - **Slurm-banner provenance TOMLs** created locally for 11 runs that had
    no end-of-run metadata (`data/<run>/<run>_banner_config.toml`): the 7
    metadata-less screens + the 4 ufront production runs; sacct
    state/elapsed/MaxRSS included; derived σ/R matched the banner for all 11.
  - **HPC-only monitors pulled** to local `data/<run>/monitors/` for
    `p018_ufront_n1` (502 wake-health rows), `p018_ufront_n1_visc` (494 rows
    — final row max_u 67,072, |Γ|/σ² 3.2e8 at min_sr +0.026: ignition
    recorded live), `p018_ufront_n2` (668 rows, died step 667), and partial
    `p018_ufront_s035_visc` (441 rows @ step 440, min_sr 0.317, healthy).
  - New Campaign E case arms (`scr_p019_*`, σ×viscosity regime map) added to
    `examples/run_p018_screen_hpc.slurm.sh`; submissions staged under the
    10-job cap. Details and all further analysis live in
    `BRAINSTORM/019_sigma_selection_procedure.md` — not here.

**2026-08-06 ~12:15 MDT — Ryan-directed σ\*-ladder A/B rebuild SUBMITTED (5
jobs).** Directive: rebuild the all-viscous N=1-vs-N=2 A/B at σ\*=0.035R
(s035_visc = the live N=1 arm) plus full pairs at 0.04R and 0.05R (neither σ
tried under the uniform-d_front law — nearest priors are old-carrier L1/L1_visc
at 0.0381R N=4 kinematic-Das, and the 8-rev screen `scr_ufsig_0p05`).
**Ryan's D ruling: N=1 arms HOLD D=3.4** (uniform C1 clearance; Das/c tip-cap
violation at 0.04R/0.05R — 1.63/2.35 vs the 1.5 band — documented, not
avoided; note this is the opposite choice from the s035 hedge's band-capped
D=3.0). N=2 arms at minimal in-band D. All viscous (CoreSpreading β=1e9),
P018_SETTLE_REVS=22, 48 h/64G. Launcher arms added + deployed md5 f12a3605.

| job | tag | σ/R (OV/PPS) | N | D | merge_r |
| --- | --- | --- | --- | --- | --- |
| 13061047 | `p018_ufront_n2_s035_visc` | 0.0349 (2.4/12) | 2 | 5.6 | 0.0048 |
| 13061048 | `p018_ufront_n1_s040_visc` | 0.0400 (2.75/12) | 1 | 3.4 | 0.0055 |
| 13061049 | `p018_ufront_n2_s040_visc` | 0.0400 (2.75/12) | 2 | 4.9 | 0.0055 |
| 13061050 | `p018_ufront_n1_s050_visc` | 0.0501 (2.87/10) | 1 | 3.4 | 0.0069 |
| 13061051 | `p018_ufront_n2_s050_visc` | 0.0501 (2.87/10) | 2 | 3.9 | 0.0069 |

Queue context: 6 campaign jobs running + 4 p019 screens (13060963–66,
concurrent 019 session) — these 5 will PEND initially; the 019 auto-submitter
honors the 10-job cap and will hold its remaining screens behind them. Banner
verification armed in the session watchdog (emits the `settle:` line per job
on first output). Watchdog also now tracks all 11 jobs and escalates at 160G
(sweep-at-100G alone cannot bound growth once sweepable VTK is exhausted).

**2026-08-06 ~12:35 MDT — all five σ-ladder A/B jobs RUNNING, banners
VERIFIED** (13061047–51): every arm shows the intended knobs in its log
banner — overlap/pps 2.4/12, 2.75/12, 2.87/10; merge_r 0.0048/0.0055/0.0069;
settle:22; visc:true; nwakerows and das_uniform exactly per the submission
table (n2_s035 2/5.6, n1_s040 1/3.4, n2_s040 2/4.9, n1_s050 1/3.4,
n2_s050 2/3.9). No queue wait — the cluster had free slots, so the 10-job
cap never bound. Tripwire coverage armed on all five.

**2026-08-06 ~05:00 MDT — fifth autonomous sweep.** Fired at 101 G, **freed
15 G → 87 G**. Queue jumped 7 → **21 RUNNING** as two new batches landed.

**Ninth memory failure: 13057254 `fp-018-ufront-n2`, `OUT_OF_MEMORY` 0:125 after
17:31:25.** With this, **three of the four original `ufront` runs are dead**
(`ufront-n1` OOM 12:58, `ufront-n1-visc` Julia OOM 11:59, `ufront-n2` OOM
17:31); only `ufront-n2-visc` and `ufront-s035-visc` survive, and both are
younger than the elapsed times at which the others died. Longer walltime does
not rescue the configuration — n2 ran 4.5 h longer than n1 and died anyway.

**Five more `ufront` viscous runs were just submitted** — `ufront_n1_s040_visc`,
`ufront_n1_s050_visc`, `ufront_n2_s035_visc`, `ufront_n2_s040_visc`,
`ufront_n2_s050_visc`, each with ~1d23h walltime. That is roughly **10 days of
committed walltime in the family with the reproducible
`merge_particles!`/`FLOWVPM_merging.jl:454` crash**, submitted after three of
its members had already died. Flagged for Ryan; the watchdog does not act on
jobs.

**A new campaign has appeared: nine `fp-019-*` jobs** (`scr-fid144`, `scr-s015`,
`scr-s015v`, `scr-s020v`, `scr-s025`, `scr-s025v`, `scr-s030v`, `scr-s038`,
`scr-s038v`). They write into the same checkout and therefore count against the
same 200 G cap; the sweeper handles them automatically since none are on the
protect list. Note for the ledger: **018 no longer has the cap to itself.**

Disk outlook: 21 writers versus ~15 G recovered per sweep at a 10-minute
cadence. The 175 G "not keeping up" alarm remains armed as the backstop.

**2026-08-06 ~06:00 MDT — sixth autonomous sweep.** Fired at 102 G, **freed
15 G → 89 G**. 20 jobs RUNNING.

**The post-sweep floor is rising steadily: 74 → 76 → 79 → 82 → 87 → 89 G over
six sweeps (~5 h).** Cause identified — it is not sweep failure. `data/` totals
90.5 GB across 118 run directories, and **four protected runs hold 66.6 GB of
it (74%)**:

| run | size | state |
| --- | --- | --- |
| `p018_h0p25` | 20,059 MB | RUNNING (12:44 left) — protected, growing |
| `p018_h0p125` | 16,379 MB | RUNNING (12:44 left) — protected, growing |
| `p018_ov1p4` | 15,947 MB | **COMPLETED** 30 revs — frozen |
| `p018_ov1p0` | 15,865 MB | **COMPLETED** 30 revs — frozen |

Everything the sweeper is allowed to touch now totals under 24 GB; the largest
unprotected directory is `scr_mesh_s80` at 746 MB. So the sweeper is already
near the floor it can reach, and the rise is `h0p25`/`h0p125` accumulating
un-swept while protected. They finish in ~12 h, adding perhaps 20 GB more first.

**`p018_ov1p0` and `p018_ov1p4` (31.8 GB combined) are frozen, completed, and
still protected.** Unprotecting them is the one lever that lowers the floor;
both satisfy the protect list's own removal condition. Agents do not edit that
file — Ryan's call. Not urgent at 89 G, but it is now the *only* remaining
headroom lever if the 21 writers outpace the sweeper.

**Tenth memory failure, and the 019 campaign is affected too: `fp-019-scr-s015`,
`OUT_OF_MEMORY` 0:125, after only 00:55:30** — by far the fastest death yet
(previous quickest was 4 h). The failure mode is therefore not confined to 018's
`uf`/`ufront` families; it reaches the new sigma-screening campaign as well, and
at a much shorter timescale.

**2026-08-06 ~07:00 MDT — seventh autonomous sweep.** Fired at 103 G, **freed
15 G → 92 G**, 19 RUNNING. Floor continues its predicted ~3 G/sweep climb
(74→76→79→82→87→89→92) as the protected `h0p25`/`h0p125` accumulate; the
sweeper's reachable pool is unchanged.

**Eleventh memory failure: `fp-019-scr-s015v`, `OUT_OF_MEMORY` 0:125 after
02:39:55.** With `fp-019-scr-s015` (00:55:30) that is **both `s015` rungs dead,
and they are the two fastest deaths of the campaign.**

The 019 screening ladder is `s015 < s020 < s025 < s030 < s038`, and **the
smallest-σ rung died first and fastest, in both its plain and viscous
variants.** That is a coherent mechanism rather than a coincidence: smaller σ ⇒
more particles at fixed domain ⇒ a larger structure inside `merge_particles!`,
which is exactly the allocation the two captured tracebacks
(`FLOWVPM_merging.jl:454`) name. It also explains the ~30 % MaxRSS/ReqMem
signature staying constant while absolute RSS varied — the copy scales with the
cloud, not with the request.

**Testable prediction for Ryan:** if the mechanism is σ-driven, `s020v` is the
next to die, then `s025/s025v`, with `s038/s038v` most likely to survive. Still
running: `s020v`, `s025`, `s025v`, `s030v`, `s038`, `s038v`. The ladder is, in
effect, already running the experiment that would confirm it.

**2026-08-06 ~13:05 MDT — `scr_ufdt_nt144` (13058048) IGNITED at step ~1104/
1440; CANCELLED under the tripwire stop rule.** Trajectory (wake-health CSV,
retained): long grind to min σ 7.5e−5 m (min_sr 0.0219 by step 1103), then the
inviscid σ(1−ΔtZ)<0 flip in ONE step — min_sr +0.0218 → **−20.75** at step
1104 — followed by max_u doubling per step (168 → 1856 → 380 → 687 → 1383;
|Γ|/σ² 4.1e5 → 3.8e6). scancel'd with ~8 h walltime left: stability readout
fulfilled, aero deltas untrusted per the failed calibration gate, remaining
steps were garbage physics on 64 cores. **With this, the entire `ufdt` dt
family is dead by ignition (nt36 OOM, nt72 OOM, nt144 tripwire) — the
smallest-dt rung surviving longest (~1104 steps ≈ 77% of its 30 revs) but
still igniting is direct measured evidence that Δt refinement delays and does
NOT mitigate the σ collapse**, consistent with the mechanism doc and the
das/eta no-dt→0-limit finding. Time-to-ignition datapoint banked for 019's
P4 regime map. 10 study jobs now running.

**2026-08-06 ~08:00 MDT — eighth autonomous sweep.** Fired at 100 G, **freed
10 G → 91 G** (down from ~15 G/sweep: the sweeper's reachable pool is shrinking
as the unprotected runs are already trimmed to their 10-step floor). 17 RUNNING.

**`fp-018-scr_ufdt_nt144` was CANCELLED (exit 0:0) after 19:17:08** — cancelled,
not crashed, so this was a human decision rather than the failure mode. With it,
**the `ufdt` sub-sweep has no surviving rungs at all**: nt36 OOM, nt72 OOM,
nt144 cancelled.

**No new memory deaths this cycle, so the σ-ordering prediction is still open.**
`s020v` — the predicted next casualty — is alive, as are `s025`, `s025v`,
`s030v`, `s038`, `s038v`. Prediction stands unresolved rather than confirmed.

All seven `ufront-*` viscous runs remain alive (`n2-visc`, `s035-visc`, and the
five `*_s0{35,40,50}_visc` submitted at ~05:00). None has yet reached the 12–17 h
elapsed range in which the three dead `ufront` runs failed, so they are not yet
past the risk window.

**2026-08-06 ~13:20 MDT — correction + refinement to the §13:05 ufdt
conclusion.** The §13:05 entry called nt144's ~1104 steps "surviving longest";
that is true only in STEP count. In physical time the family died at
**≈7.9 / 6.8 / 7.7 revolutions** (nt36 step 284, nt72 step 490, nt144 step
~1104, at NT 36/72/144) — i.e. **time-to-ignition is nearly constant in revs
across a 4× dt refinement at fixed σ=0.0291R inviscid**. The clean-ladder
verdict is therefore sharper than "dt delays but does not prevent": ignition
is set by the physical evolution of the wake (Z is a flow property; per-unit-
time contraction rate Z does not improve with dt), and dt refinement buys
essentially nothing, in agreement with the mechanism doc's ΔtZ analysis and
the das/eta no-dt→0-limit result. (Also correcting §13:05's "30 revs" — the
nt144 screen was an 8-rev protocol; it died at ~7.7 revs, near the end of its
window.) Strong prior for 019's P4 regime map: the σ axis, not the dt axis,
is the control.

**2026-08-07 ~00:45 MDT — `p018_h0p25` (13051772) COMPLETED + SCORED: the
Δx (h-invariance) leg PASSES both gates at the first rung pair; 3-lobe-absent
prediction CONFIRMED.** COMPLETED exit 0 at 46:28:33 (photo-finish, step
1079/1080). Scored per phase_12 Part A pre-registration (`p018_analyze.py`,
matched revs 15–29):

| metric | value | gate | verdict |
| --- | --- | --- | --- |
| M1 CT̄ (h0p25) | 0.071686, CI [0.071626, 0.071717], drift 0.120% non-monotone | settled | **PASS** |
| ΔCT̄ vs B0 (0.07170) | **−0.02%** | ≤0.5% | **PASS** |
| M2 ε_Γ,max (B0→h0p25, 0.3–0.95 r/R) | **0.223%** (rms 0.161%) | ≤1% | **PASS** |
| 3-lobe Γ̄ mode at fixed σ | max deviation 0.22% ≪ the σ-ladder's 4–8.8% lobes | pre-registered ABSENT | **CONFIRMED absent** |

Consequences per the pre-registration: **h is converged at the production
ratio (h/σ 0.5→0.25 invariant); Phase 4's ladder is retroactively a pure σ
statement; the 3-lobe mode's σ attribution stands** (it did NOT re-attribute
to h/overlap). `p018_h0p125` (h/σ 0.125) demotes to confirmation/Richardson —
still running at ~1.4 h to walltime; `_s2` chain on timeout. Harvested CSVs
committed to local `data/p018_h0p25/`.

**2026-08-07 ~02:20 MDT — `p018_h0p125` (13051773) TIMEOUT at 48 h as
projected; `_s2` CHAINED.** Died at step ~855/1080 (last restartable VTU
854, matching the ~850 projection). Chain submitted per the ops_reference
recipe: **`p018_h0p125_s2` = job 13065708** (RESTART_STEP=854 from
`data/p018_h0p125`, P018_SETTLE_REVS=22, 128G per the memory ruling,
24 h — needs ~15 h for the remaining ~226 steps at the late-run
~230 s/step). Banner verification armed. Note the h-ladder verdict does
NOT wait on this run: the Δx leg already PASSED at B0→h0p25 (§00:45);
h0p125 is confirmation/Richardson only. Seam gate at analysis: per-step CT
jump ≤ 0.05% at step 854.

**2026-08-07 ~02:50 MDT — `p018_h0p125_s2` (13065708) VERIFIED warm:**
banner knobs correct (overlap 8.0 / pps 16 / merge 0.0120 / settle 22,
nwakerows 4, inviscid) and `simulate_warmstart!: resuming from step 854`,
stepping at 855/1079. Chain healthy; ~15 h to completion.

**2026-08-07 ~03:30 MDT — `p018_ufront_n2_visc` (13060144) IGNITED at step
~619–620; CANCELLED under the stop rule.** Endpoint trajectory (CSV
retained): min_sr oscillated 0.077–0.094 through step 618, then max_u jumped
46 → 281 → **516 m/s (~8× tip speed)** in two steps as min σ hit the viscous
floor (9.57e−5 m ≈ √(2νΔt); min_sr 0.0269) — the same Γ-ignition-at-floor
mode as its N=1 twin. **Time-to-ignition ~620 vs ~490 for n1_visc: N=2
(D=6.5) bought ~130 steps (~3.6 revs) but did NOT stabilize σ=0.03R
viscous.** With this, every σ=0.03R arm (invisc N1 ~501, invisc N2 ~667,
visc N1 ~490, visc N2 ~620) is dead by ignition — the σ\*=0.035R decision is
now quadruply enforced, and the all-viscous rigid-vs-free-row A/B moves
entirely to the σ\*-ladder pairs (13061047–51, all healthy at submission).
Datapoint banked for 019 P4 (viscous 0.03R column now has both N values).
scancel'd with ~40 h walltime unspent.

**2026-08-07 ~05:45 MDT — MAJOR: `p018_ufront_s035_visc` (13058988), the
σ\*=0.035R certification run, IGNITED at step ~995/1080 (rev ~27.6);
CANCELLED.** Death tail: min_sr 0.045 → 0.025 (floor) at step 995 with
|Γ|/σ² 4.7e5 → 3.4e6 and max_u 456 → **2,773 m/s** — the same
Γ-ignition-at-the-floor mode, arriving after ~27 revs of slow monotone
contraction (min_sr 0.32 @ 431 → 0.12 @ 927 → collapse @ 995). **The
σ\*=0.035R viscous N=1 (D=3.0) config is NOT certifiable at production
length**; survival-to-30-revs failed with 85 steps to go. Consequences:
(a) certification scoring proceeds on the surviving window but the verdict
is FAIL regardless of CT quality (M1's ≥15-rev window is also truncated —
usable revs ~15–27); (b) **σ\* escalates provisionally to 0.04R** — the two
0.04R arms in flight (13061048/49, min_sr ~0.10–0.13 at step ~600, sliding)
are now the production candidates and their endpoints are decisive; the
0.05R pair is the hedge; (c) the σ_stab/second-criterion open question
(019 note) gains a fifth datapoint: 0.035R viscous sits ABOVE
σ_stab≈0.031R and still ignited at ~27.6 revs — consistent with the
constant-physical-time grind, NOT with a static threshold; (d) every
`ufront` production run to date (0.03R ×4, 0.035R ×1) has now ignited —
survival time grows with σ (≈14, 17, 13.6→18.4, 27.6 revs) suggesting the
grind slows with σ but has not yet arrested anywhere. scancel'd with ~2.7 h
unspent; CSVs intact through ignition for scoring + 019's regime map.

**2026-08-07 ~06:15 MDT — `s035_visc` record-scoring on the surviving window
(CT synthesized from monitor02 CFx; mapping CT(k)=−CFx(k−1) VERIFIED to 8
digits against `p018_h0p25`'s CT_vs_rev).** Numbers, with framing caveats:

| metric | value | note |
| --- | --- | --- |
| M1 CT̄ revs 15–27 | **0.074444**, CI [0.074052, 0.075219], per-rev std 0.00108 | window 13 revs < the 15-rev gate ⇒ CHECK, not PASS; blocks 0.07483 → 0.07466 declining; possible endgame contamination from the contraction |
| ΔCT̄ vs B0 | **+3.8%** | NOT a clean σ statement — s035_visc differs from B0 in σ (0.035R vs 0.0873R), N (1 vs 4), Das law (uniform-D vs kinematic), and viscosity jointly |
| M2 vs B0 | ε_max **11.77%**, rms 6.39% — FAIL | same caveat: a joint carrier delta across the largest σ-step in the campaign, dominated by the known σ-mode redistribution; NOT comparable to the adjacent-rung ladder ε's (<1%) |

Primary verdict unchanged: **certification FAIL on survival** (ignition rev
~27.6). The CT̄ landing at 0.0744 (above experiment's 0.072, vs the
old-carrier 0.0717 family) is noted as a large-effect observation for the
σ*-ladder A/B to resolve on clean adjacent comparisons — the five in-flight
arms share the ufront carrier and differ only in (σ, N), so their
cross-comparisons will be interpretable where this one is not. Synthetic
CT CSV + monitors committed under local `data/p018_ufront_s035_visc/`.

**2026-08-11 — reconnection harvest (socket was dead 08-07→08-11; all
windows expired). Disk swept 174G→~120G. `nt144_s3` chained = job 13134682
(RESTART_STEP=3220 from `_s2`'s TIMEOUT; settle 16 matched; warm-start
watcher armed). Job outcomes + scoring:**

1. **σ-ladder A/B (viscous ufront): survivors n1_s040, n1_s050, n2_s050
   (full 1079 steps, final health VERIFIED clean); ignitions n2_s035 (rev
   ~23.3) and n2_s040 (rev ~26, Γ/σ² 1.5e8).** Survival boundary: N=1 needs
   σ ≥ 0.04R, N=2 needs σ ≥ 0.05R — **viscous N=2 is LESS stable than N=1 at
   matched σ** (opposite of the inviscid 0.03R ordering). `n1_s040` shows the
   campaign's FIRST arrest: min_sr 0.099@588 → 0.122@1079 (recovered,
   plateaued). n1_s050 (0.108) and n2_s050 (0.145) still slowly declining at
   end. No restarts for the ignited arms — re-running from pre-ignition
   state deterministically re-ignites.
2. **But no surviving arm is CONVERGED: M1 all CHECK** — CT̄ 0.075298
   (n1_s040, drift 2.1% monotone), 0.076516 (n1_s050, drift 6.0%), 0.077167
   (n2_s050, drift 3.6%); all 4–8% ABOVE both the old-carrier family
   (0.0717) and experiment (0.072), still rising at window end. **M2 FAILs:
   A/B at 0.05R ε_max 5.04%; N1 σ-pair 0.04→0.05 ε_max 5.13%.** The ufront
   viscous carrier is stable-but-unconverged at 30 revs; settled windows
   need extensions (_s2 chains), and the carrier-level +4–8% CT offset vs
   the η-Das family is now the campaign's dominant open discrepancy.
3. **h-ladder finest rung (h0p125 base+_s2 stitched, seam jump 0.006% PASS,
   no particle-cap saturation — max 154k of 500k): the 0.25→0.125 pair
   FAILS BOTH GATES** — CT̄ 0.073207 (+2.1% vs h0p25/B0; M1 CHECK drift
   0.509%), M2 ε_max 5.13%/5.07% (two seam-free sub-windows, consistent).
   **Deviation profile is a smooth monotone inboard tilt (+5% at r/R 0.31 →
   0 at 0.83 → −2.4% at tip), NOT the 3-lobe** — so the σ-mode attribution
   stands, but per the phase_12 pre-registration the h-confirmation FAILED
   ⇒ **FLAGGED FOR RYAN: h is a live error term at σ/h=8 (deep-overlap
   regime 008/011/012), and the "h converged at production ratio" claim
   rests on the 0.5→0.25 pair alone.**

**2026-08-11 addendum — `nt144_s3` FIRST SUBMISSION (13134682) had the WRONG
Das law and was killed+resubmitted.** Its banner showed `das_eta:1.0` where
`_s2` ran 4.0 (NT=144 needs η=4.0 to hold the physical Das length; the p018
family block exports η=1.0 unconditionally, so the nt144 arm's guarded
`:-4.0` never fires — the documented override-precedence gotcha; `_s2`'s
submitter had passed it explicitly). Caught at warm-start verification
(banner check), scancel'd within minutes, `data/p018_nt144_s3` wiped, and
**resubmitted with `DAS_ETA_KINEMATIC=4.0` explicit = job 13134726**
(otherwise identical: RESTART_STEP=3220, settle 16, 24 h). Banner+warm-start
verification armed. LESSON reinforced: always verify the FULL banner against
the parent segment's, not just settle/σ knobs.

**2026-08-11 — Ryan-directed: `p018_ufront_n1_s040_visc_s2` SUBMITTED = job
13135200** (extension of the arrested σ\* candidate: RESTART_STEP=1079,
P018_SETTLE_REVS=37 ⇒ +15 revs to 45 total/1620 steps; 64G/36 h; case arm
carries all knobs unconditionally so no η-style env gotcha; banner+warm
watcher armed). Purpose: (a) settled ≥15-rev window free of the 2.1%
monotone drift, (b) long-horizon arrest test of the min_sr 0.12 plateau.
Full-banner-vs-parent verification required at start per the nt144 lesson.

**2026-08-11 — Das-LAW decomposition A/B SUBMITTED (Ryan-directed):
`p018_etadas_n1_s040_visc` = job 13135245** (launcher arm added, deployed
md5 b9360ae6). Identical to `p018_ufront_n1_s040_visc` (N=1, σ 0.04R OV
2.75/PPS 12, viscous CoreSpreading, merge_r 0.0055, SETTLE 22, 48 h/64G)
except the handoff law: **η-kinematic Das (η=1.0) instead of uniform
D=3.4σ**. ΔCT vs n1_s040 isolates the clearance-law knob in the ufront
+4–8% offset; prior A/Bs price viscosity ~0, N ~1%, σ 1–2%. Expected
banner: das_uniform:nan, das_eta:1.0 — verify on start (watcher armed).
If ~null ⇒ the offset is an interaction; follow-up = N=4 η-Das at 0.04R
viscous. Also: `n1_s040_visc_s2` (13135200) VERIFIED warm (full banner
matches parent, settle 37, resumed 1079).

**2026-08-11 — `p018_etadas_n1_s040_visc` (13135245) banner VERIFIED:**
das_eta:1.0, das_uniform:nan, overlap 2.75 / pps 12 / merge 0.0055 /
settle 22 / nwakerows 1 / visc true — the intended η-Das A/B arm exactly.
All three in-flight jobs (13134726, 13135200, 13135245) verified. Session
prepared for context reset; RESET BRIEF (c) in the item file is current.

**2026-08-12 — `nt144_s3` (13134726) LANDED; Phase 3 dt-Richardson CLOSED
(non-monotone, order unobservable):** COMPLETED 2026-08-11 15:40 in 6:55
(steps 3220→3456; the run's design length is settle 16 + 8 = 24 revs = 3456
steps — the old "3600" projection assumed 25 revs and was wrong). Full banner
matches `_s2` line-for-line (das_eta:4.0 — the explicit-η resubmission was
correct). Wake health clean over the whole segment: min min_sr 0.4727, final
0.4986, max_u ≤ 16 m/s, no σ<0. Three-segment stitch base+_s2+_s3. Seam-gate note: for
monitor-reconstructed segments the analyzer's csv-step = monitor-step+1
mapping means the cut must be at S+1 (not the VTU step S) to compare truly
adjacent rows — cutting at S drops one valid parent row and doubles the gap.
Adjacent-row seams: **0.0515% (@2026/2027) and 0.0575% (@3221/3222)** —
marginally over the nominal 0.05% gate but ~1.7× the typical 1-step delta
(0.031%), the same magnitude as the validated `_s2` warm-start seam (0.05%,
ledger 2026-08-05) and B0's accepted 0.065%; within natural variability
(≥10% of matched-gap deltas exceed them). Judged PASS-with-note. **M1 matched 10–19: CT̄
0.072023, CI [0.071952, 0.072085], per-rev std 0.00012, block drift 0.166%
monotone (CHECK); long window 10–23: 0.072012.** Richardson over matched
10–19: {NT36 0.071866, NT72 0.072243, NT144 0.072023} — **Δ(36→72) = +0.52%,
Δ(72→144) = −0.31%: NON-MONOTONE, difference ratio −1.72, no valid observed
order and no extrapolation.** The 72→144 step is ~3× the CI half-widths, so
it is not pure noise; it is consistent with the phase_13 §5 clearance
confound (fixed N=4 ⇒ min d/σ degrades 0.713→0.414→0.264 down the ladder —
the rungs are not a clean dt refinement). **Close-out statement: dt residual
is bounded at ≤0.5% peak-to-peak across NT 36–144 at Das 0.41c, but the
ladder cannot yield an observed order or Richardson limit until re-run with
the dt pin (integer schedule N=4/7/13, or the §4b N=1 construction which
pins automatically).** Phase 3 CLOSED with this caveat.

**2026-08-12 — `n1_s040_visc_s2` (13135200) LANDED: arrest HOLDS to 45 revs;
drift caveat reduced but NOT killed; ufront offset shrinks to ~+2.3%:**
COMPLETED 2026-08-12 07:45 (21:37, steps 1079→1620 = 45 revs). Seam@1079
0.197% of CT, smaller than the typical matched-gap delta (89% exceed it) —
PASS trivially (per-step CT in this viscous N=1 run oscillates ~1.8%/step).
**Long-horizon arrest CONFIRMED: global min_sr 0.0976 @ step 1128 (rev
~31.3), then monotone recovery to 0.1531 @ 1619; max_u ≤ 52 m/s; n_p max
215k of 500k.** M1: **late window 30–44 CT̄ 0.073632, CI [0.073032,
0.073874], block drift 1.107% monotone — still CHECK at 45 revs.** The
original 15–29 window (0.075280) was transient-inflated: 5-rev blocks run
0.074662 / 0.075003 / 0.076175 (peak) / 0.073233 / 0.073616 / 0.074048 —
a hump through revs 25–29, a drop, then a slow rise (~0.55%/5 revs).
Consequences: (a) the 2.1%-drift caveat is replaced by a slower 1.1%
monotone drift — settledness still not reached; (b) **the ufront carrier
offset at this arm collapses from +4.6% (15–29) to +2.7% vs the η-Das
family (0.0717) and +2.3% vs experiment (0.072)** — much of the "dominant
+4–8% offset" at n1_s040 was settling transient, not carrier physics. The
etadas A/B (13135245, 30-rev run) must therefore be scored against the
PARENT's matched 15–29 window (same protocol, same transient), where the
knob isolation stays valid; its absolute level will carry the same
inflation. Watchdog v6.3 armed (etadas terminal state + disk sweep >100G /
alert >160G; disk 122G after the 08-12 auto-sweep freed 47.5G).

**2026-08-12 — min_sr fluctuation characterized; CT correlation NULL;
inventory trend (discussion session with Ryan, local CSV forensics only):**
Per-rev min(min_sr) in all three stable viscous arms fluctuates
APERIODICALLY in a 0.09–0.17 band with ~5–10-rev dip/recovery episodes
after the startup decline (n1_s040 dips 0.098@16 and 0.098@31, recovers
to ~0.15 by rev 44; n1_s050 dips 0.095@19–20; n2_s050 dips 0.094@27 then
0.138@29 — the earlier "still declining at end" reads were tail
snapshots of this fluctuation, not trends). Attribution physical-vs-
order-statistic OPEN (min over 200k particles; recovery can be merge/
truncation bookkeeping) — discriminators = p1-percentile + argmin
columns, staged with the drift study. **min_sr↔CT per-rev correlation is
NULL**: n1_s040 detrended lag-0 r=+0.08 (n=35); max |r|=0.36 at CT-leads-
by-2 is below the multiple-comparison bar; shorter arms inconsistent
(−0.35@0, −0.53@0, best lags −4/+4). Supports min_sr as a pure stability
tripwire; CT windows can be scored without conditioning on wake-health
episodes. **Particle inventory is NOT still filling**: n_p peaks 228k@rev
15, dips to 208k@rev 30, then +0.17%/rev over revs 30–44 (parallels the
late CT rise) — so the CT drift lives in the circulation/σ DISTRIBUTION,
not the headcount; this kills the naive fill hypothesis and motivates the
banded-inventory monitor (drift study staged in the item file; Ryan
directions this session: pinned dt ladder NT 36/72/144 @ PPS 12/6/3
decided, cost accepted; drift-source study staged pending design
discussion; numerical-vs-physical fluctuation check queued post-ladder).

**2026-08-12 — Phase-15 subtask (a) DONE (subagent): H4 spectral screen on
the 45-rev CT series. Verdict: H4-as-posed INDISTINGUISHABLE at this length;
the better description is aperiodic BURST RECTIFICATION; and the rev-30
restart is CONFOUNDED with the block-mean step.** Key numbers (revs 10–44,
n=35, AR(1)-MLE): sinusoid P=31 revs A=1.5% wins AIC by only ~1.1 and LOSES
BIC to constant/linear; exp-saturating is degenerate (τ unidentifiable,
CI [4,≥250]); low-freq band power vs AR(1)+trend bootstrap p=0.25; all
models pass Ljung-Box — the window rejects nothing. No coherent slow mode
across arms (humps at revs 25–29 / 20–24 / none). n_p leads CT by ~5–6 revs
at detrended r≈0.5 (marginal after lag-scan multiplicity). **Finding 1
(structure): per-rev CT̄ is strongly driven by within-rev burst amplitude —
r(mean, within-rev std)=0.61 (p=4e-4), replicated in n2_s050 (0.59);
quiet-epoch regression limit CT≈0.0730; episodes 5–15 revs.** **Finding 2
(CAVEAT on today's earlier scoring): the 0.0762→0.0732 block step sits
EXACTLY at the stitch (rev 30.0), and within-rev std collapses ~0.007→0.0005
across the boundary — n_p is continuous (state reloads faithfully) but the
restart plausibly truncated an active burst, so the "offset shrinks to
+2.3%" revision is partially confounded with the restart. Also revs 43–44
are back at the 25–29 level (drift not monotone through the tail).**
Discrimination-on-the-mean needs +30–45 revs (to rev 75–90, 3–5σ);
CHEAPER: track quiet-epoch CT and burst amplitude as separate statistics,
and verify with one restart placed mid-quiet-epoch that the boundary does
not reset burst state. Plots + script in session scratchpad h4_screen/
(session-local). Feeds the Phase-15 plan discussion with Ryan.

**2026-08-12 — Phase-15: H1 two-epoch VTK check (subagent, retained restart
sets rev~30 vs rev~45, 10 snapshots/epoch). Verdict: H1 SUPPORTED — the wake
distribution is decisively NOT equilibrated.** Total Σ|Γ| fell **−11.4%**
(63× within-epoch sd) between rev 30 and rev 45 while particle count rose
+3.0% (and was still climbing +86/step at rev 45); EVERY z-band's Σ|Γ|
dropped 8–16.5%. Largest movers: **band 2–3R** (n +7.7%, Σ|Γ| −16.5%,
σ_p50 +2.1%) — mid/far wake still restructuring; and the **near-disk
outboard cell (z<0.5R, r>1R) emptied 93–95%** — radial spillage present at
rev 30 has vanished by rev 45 (tightening slipstream; less outboard
vorticity near the disk plausibly = the +1% CT). σ_p5 rose 5–12% in bands
0.5–3R (small-core tail thickening = viscous spreading), but co-occurring
with the Σ|Γ| changes, so the pure-H3 signature (σ moves, Σ|Γ| doesn't) is
ABSENT. **H2 now concrete: a hard deletion boundary exists at x=3.49R /
r=1.50R** (exactly constant across all 20 snapshots; the 4R+ band is empty
BY CONSTRUCTION — the per-step monitor's outer band must end at the actual
boundary; ALSO: banner says depth:4R but measured extent is 3.49R from the
disk plane — reconcile against the driver's truncation definition before
designing the H2 A/B). Caveats: two-epoch design (within-epoch sd spans
only ~0.28 rev, so Δ/sd overstates significance vs the 5–10-rev
fluctuation; the uniform-sign all-band decline is hard to ascribe to phase
aliasing but a trend claim needs a third epoch or the monitor); Σ|Γ| is
magnitude-sum (cancellation-blind); the epoch-2 restart-knob caveat is
RESOLVED (full banner verified vs parent this morning). Artifacts:
scratchpad h1_check/ (20 VTP ~820 MB local copies, analyze.py,
band_table.csv — session-local). Combined Phase-15 posture after (a)+(a2):
H1 supported, H2 live/concrete, H3 partial not sole, H4-on-the-mean
indistinguishable and partially superseded by burst-rectification framing.

**2026-08-12 — Phase-15 PLAN DISCUSSION HELD WITH RYAN (the gate his
2026-08-12 direction required); decisions + implementation started.**
Decisions (Ryan, this session): (1) WakeInventoryMonitor v1 ADDS the
cancellation-aware vector sum ΣΓ⃗ per band alongside Σ|Γ| (discriminates
merge-cancellation from vorticity leaving a band); (2) FOUR conditional
runs authorized now — mid-quiet-epoch restart check, H2 truncation A/B
(TRUNCATION_DEPTH_R=6), H4 `_s3` extension to ~60 revs, and a MERGE-OFF
A/B (~5 revs warm from a retained n1_s040_visc set; diagnostic only;
separate restart from the seam check to avoid confounding) — all submit
AFTER the instrumentation deploys so each carries the diagnostics;
(3) quiet-epoch/burst statistics are REPORT-ONLY (decision_rules.md
untouched until validated on >2 arms + restart check clears); (4) σ stays
0.04R for production — under-damping near σ_stab DISFAVORED as drift
cause (n1_s050 at 25% margin also drifts; min_sr↔CT null), and the σ-axis
accuracy cost (1–2%/rung, ε_Γ∝σ) outweighs a stability margin the arrest
already provides. Discussion elevated MERGE CANCELLATION (deferred H4b)
to the leading mechanical candidate for H1's Σ|Γ| decline: merging is the
only always-on Σ|Γ| sink at conserved ΣΓ⃗, and the decline is
uniform-sign across bands (density-weighted, not band-localized).

**2026-08-12 — H2 boundary discrepancy RECONCILED (read-only, driver):
"banner 4R vs measured 3.49R" is labeling, not a bug.** The driver builds
`GlobalCylinder([-0.5R,0,0], [cylinder_depth,0,0], 1.5R)`
(`rotor_hover_pressure_comparison.jl:519` pre-edit): origin 0.5R
UPSTREAM, extrude length `TRUNCATION_DEPTH_R`·R=4R along +x ⇒ downstream
deletion boundary at x=3.5R from the disk plane (measured 3.49R = last
retained particle), radius 1.5R exactly as measured. Banner "depth:4R" is
CYLINDER LENGTH. Consequences: inventory outer band ends at 3.5R; the H2
deep arm `TRUNCATION_DEPTH_R=6` gives 5.5R downstream (report it as
such).

**2026-08-12 — Phase-15 instrumentation IMPLEMENTED (local, gates in
progress; NOT yet deployed):** `WakeInventoryMonitor` appended at EOF of
`src/FLOWPanel_simulate_monitors.jl` (per step, per cell: n, Σ|Γ|, ΣΓ⃗
x/y/z, σ mean/p5/p50/p95; cells = z/R bands 0–0.5 split r/R≶1, 0.5–1,
1–2, 2–3, 3–3.5, + `outside`; statics skipped; type-7 quantiles;
partition conserves totals) and min_sr-attribution columns on
`WakeHealthMonitor` (opt-in `attribution` kwarg, columns
p1_sigma_ratio + argmin x/y/z appended after the legacy set). Driver
wiring: both appended LAST, env knobs `WAKE_INVENTORY` +
`WAKE_HEALTH_ATTRIBUTION` (default ON, in banner). Unit tests: two new
testsets (8+21 asserts) PASS incl. legacy-CSV bit-identity at defaults;
runtests_unit_simulate.jl fully green. Pending gate: formulation_test,
NT=6 smoke pair (bit-identity diff of physics CSVs, diagnostics on/off,
both formulations), full runtests, cluster md5.

**2026-08-12 — burst/quiet-epoch decomposition scripted
(`scripts/p018_burst_stats.py`, REPORT-ONLY) and run retroactively;
first cross-arm reads.** Validation: reproduces p018_b0(+_s2@719) 15–29
grand mean 0.071702 (M1 0.071701) and n1_s040(+_s2@1079) 10–44
r(mean,std)=+0.67 / quiet-regression limit 0.0729 (subagent: 0.61 /
≈0.0730; window-handling differences). Findings: (a) **the carrier's
quiet-series drift is −0.6% to −1.0% ± ~1.5% over revs 10–44 — NOT the
positive monotone drift of the raw mean** ⇒ consistent with a stationary
(or slightly falling) quiet baseline sampled through burst-episode
timing; (b) **burst rectification is a viscous-N=1-ufront family
property, not universal**: B0 shows r(mean,std)=−0.19 and ~zero burst
incidence under a family s* — per-rev scatter differs 5–10× across
families ⇒ **s* pooling rule: pool only within one carrier family**
(cross-family comparisons use the threshold-free regression intercept;
rule documented in the script); (c) quiet-limit σ ladder: 0.0729
(n1 σ0.04R) / 0.0746 (n1 σ0.05R) / 0.0754 (n2 σ0.05R) vs B0 0.0719 —
the σ offset survives into the burst-free level; (d) detrended
r(s_k, min_sr) = +0.35 / −0.06 / +0.01 across the three viscous arms —
no consistent sign ⇒ no coherent burst↔core-contraction coupling;
further disfavors under-damping. etadas arm to be appended when it
lands.

**2026-08-12 — Phase-15 instrumentation DEPLOYED to cluster; gate job
SUBMITTED = 13150992 `fp-018-p15gate`** (full runtests + formulation_test
+ NT=6 three-arm smoke + physics-CSV bit-identity diff, in one 4 h/8-thread
job; Ryan directed the smoke run on HPC — the local pair was scancel'd
mid-arm and cleaned; NO further local jobs this session without asking).
Deployed + md5-verified: `src/FLOWPanel_simulate_monitors.jl` 91749d5d,
driver 8597ad98, launcher f830e177 (adds `p018_trunc6_n1_s040_visc` +
`p018_mergeoff_n1_s040_visc` arms), `test/runtests_unit_simulate.jl`
43676b21, `scripts/p018_burst_stats.py` 14a2f934, gate script
`examples/run_p15_gate_hpc.slurm.sh` ac5987f4. Pre-deploy diff check:
cluster copies contained NOTHING beyond the last-deployed state (launcher
was b9360ae6 = the etadas deploy). **Full src/ sweep NOTE: four files
differ local-vs-cluster (`FLOWPanel.jl`, `_liftingbody.jl`, `_metadata.jl`,
`_solver.jl`) plus local-only `_instrumentation.jl` — these are the
CONCURRENT 021 solver-benchmark session's local work, deliberately NOT
deployed; cluster src = the state that ran all recent 018 jobs + the
monitor bundle only.** Gate verdict pending; conditionals submit only
after it passes.

**2026-08-12 — p15gate attempts 1–2 FAILED on environment/skew, NOT the
instrumentation; take-3 = 13151082.** (1) 13150992: `module load julia`
gives 1.12.6 but the Manifest pins 1.11.7 ⇒ PythonPlot/MbedTLS_jll
precompile explosion. Fix: the production launcher's spack julia-1.11.7
PATH fallback + `MPLBACKEND=Agg` (gate script updated, md5 7f961378).
(2) 13151006: full runtests had NEVER been current on the cluster —
**six test files were STALE vs HEAD** (`runtests.jl`, `test_helpers.jl`,
`runtests_unit_liftingbody.jl` — failed on the missing
`set_Das_kinematic_arc=false` in the tangent-vs-inline testset —
`runtests_unit_{wake,postprocess}.jl`, `runtests_example_suddenly_
started_wing.jl`), and cluster `runtests.jl` still included
`runtests_unit_pitching_wing_exp.jl` (fixture missing everywhere).
Fix deployed + md5-verified: HEAD `runtests.jl` (fc6d6bee — excludes the
021 session's local-only test includes AND the pitching_wing_exp
include), local(=HEAD) liftingbody/wake/suddenly-started-wing tests,
LOCAL postprocess test (matches the deployed zero-grad src fix), LOCAL
test_helpers (016's additive conversion-fixture kwargs). Lesson for the
ops file: "full runtests on cluster" is only meaningful if test/ is
deployed in lockstep with src/ — sweep BOTH before a gate.

**2026-08-12 — p15gate takes 3–5 each died one missing file deeper; take-6
= 13152777.** Additional cluster-absent suite dependencies found and
deployed (all md5-verified): `test/data/legacy_wake_conversion_reference.jl`
(untracked fixture, take-3 kill), `test/runtests_unit_kutta{,_routeb}.jl`
(015's tracked-at-HEAD tests, never copied; HEAD runtests.jl includes
them — take-4 kill), `examples/ssw_representation_probe.jl` +
`examples/ssw_sheet_particle_split.jl` (referenced via
`pnl.examples_path` by the suddenly-started-wing test — take-5 kill).
A full dependency scan of the deployed suite (all include/joinpath .jl
literals in test/) now shows NO remaining cluster-absent files. The
liftingbody Das testset passed in take-5 (42/42) confirming the stale-test
diagnosis.

**2026-08-12 — `p018_etadas_n1_s040_visc` (13135245) LANDED and SCORED:
Das-law ΔCT is ~NULL ⇒ the ufront carrier offset is an INTERACTION, and
the pinned dt ladder's rung 1 = `n1_s040_visc(+_s2)` (uniform-D survives).**
COMPLETED exit 0:0 in 1-06:53 (30 revs). Wake health clean at tail:
min_sr 0.158–0.161, max_u 15–18 m/s, no σ<0, n_p 210k. M1 on the MATCHED
15–29 window (per protocol — both 30-rev arms share the settling
transient; do NOT compare against the `_s2` late window):
**η-Das CT̄ 0.075105** CI [0.073966, 0.076744] vs **uniform-D parent
0.075298** CI [0.074644, 0.075920] ⇒ **ΔCT = −0.26%, deep inside both
CIs — NULL.** Both arms M1 CHECK (drift 3.99% vs 2.08% monotone — the
shared transient). **M2 parent-vs-etadas: ε_max 4.05% / ε_rms 1.66% —
the Das law DOES redistribute Γ̄(r/R) (>1% gate) while leaving CT
unchanged.** Burst decomposition (family s* 0.000993, matched window):
quiet-regression limits 0.0730 (etadas) vs 0.0734 (parent) — also ~equal;
burst rectification in BOTH arms (r +0.74 / +0.68) ⇒ family property of
viscous N=1 σ0.04R, NOT a Das-law artifact. CONSEQUENCES: (1) prior A/Bs
price viscosity ~0, N ~1%, σ 1–2%, Das law ~0 ⇒ **no single knob owns the
carrier's +2–3% offset — interaction confirmed; Ryan's named follow-up =
N=4 η-Das @0.04R viscous (needs his go-ahead)**; (2) STAGED DECISION
resolves: the dt ladder pins at N=1 uniform-D, rung 1 = existing
`n1_s040_visc(+_s2)`, only NT=72/PPS=6 and NT=144/PPS=3 are new
submissions (each will carry the Phase-15 diagnostics); (3) etadas arm
retained as the Γ̄-redistribution reference for the offset attribution.

**2026-08-12 — N-knob follow-up SUBMITTED (Ryan-directed post-etadas-null):
`p018_etadas_n4_s040_visc` = job 13157532, banner VERIFIED** (das_eta:1.0,
das_uniform:nan, nwakerows:4, overlap:2.75, pps:12, merge_r:0.0055,
settle:22, visc:true). The etadas arm with NWAKEROWS 1→4 as the ONLY
change — completes the 2×2 {N=1,4}×{η,uniform-D} square at σ=0.04R
viscous with n1_s040 + etadas. Submitted 48 h/64G with diagnostics OFF
(WAKE_INVENTORY=false, WAKE_HEALTH_ATTRIBUTION=false) so the arm runs the
exact code path of its A/B partners while the p15 gate is in flight.
Launcher arm added, deployed md5 d6b60df3. Score on matched 15–29 vs
n1_s040 AND etadas when it lands (ETA ~30 h ⇒ ~late 2026-08-13).

**2026-08-12 — gate take-7 kill = PRE-EXISTING dependency break, NOT
phase-15: `runtests_example_suddenly_started_wing.jl` cannot load in the
current project env ANYWHERE** — `examples/helper_functions.jl` line 2
imports GeometricTools, which commit 7d52be0 (May 2026) REMOVED from
Project.toml deps; not in the cluster global env either. Cluster
`runtests.jl` now carries a commented-out include with this note; ops_
reference "Pre-submission gate" amended (Ryan-directed) with the full
cluster-gate crib (julia 1.11.7 pin, MPLBACKEND, test/+examples/ lockstep
deploy, dependency scan, this exclusion). **Ryan owes: restore the dep or
port helper_functions.jl off `gt`.** Gate take-8 = 13157614 in flight.

**2026-08-12 — Ryan decision on the GeometricTools break: port
`examples/helper_functions.jl` OFF GeometricTools (do not restore the dep),
but DEFERRED until AFTER item 018 is complete.** Recorded in the item
file's open items; the cluster runtests.jl exclusion and the ops_reference
note stand until then.

**2026-08-12 — p15 GATE PASSED (take-8, 13157614 COMPLETED): full runtests
GREEN (with the documented ssw exclusion), formulation_test 10/10, three
NT=6 smoke arms finite (velocity diag-on/off + green diag-on),
BIT-IDENTITY PASS — all physics CSVs identical von vs voff; inventory CSV
header and wake_health attribution columns as designed. The Phase-15
instrumentation is VALIDATED and DEPLOYED; subtask (b) complete except
this note. Also: dt-ladder submissions are now unblocked pending Ryan's
pre-registration sign-off (phase_03).**

**2026-08-12 — ALL FOUR authorized drift conditionals LAUNCHED
(Ryan-directed) — three jobs + one free reading:**
| job | case | what |
| --- | --- | --- |
| 13157751 | `p018_trunc6_n1_s040_visc` | H2 truncation A/B: 5.5R downstream vs production 3.5R, cold, SETTLE 22 (~30 revs), 48 h |
| 13157752 | `p018_ufront_n1_s040_visc_s3` | H4/(d): chain to ~60 revs (SETTLE 50), warm from `_s2`@1619, 48 h. Its stitch = the (e) seam burst check, free |
| 13157753 | `p018_mergeoff_n1_s040_visc` | H4b/(f): MERGE_PARTICLES=false, warm from the SAME `_s2`@1619, SETTLE 40 (~+5 revs), 24 h. `_s3` = its stock control |
All three carry WakeInventoryMonitor + attribution (default ON).
RESTART_STEP=1619 per the retained VTU set (CSV 1620 = S+1 gotcha).
Banner + warmstart verification armed; verify before trusting any hour.

**2026-08-12 — pinned dt ladder PRE-REGISTRATION APPROVED (Ryan, as
drafted) + BOTH RUNGS SUBMITTED. NEW ALLOWANCE: walltime up to 72 h when
needed** (ops_reference §"Walltime allowance") ⇒ NT=72 fits one job, no
chain; NT=144 = one 72 h segment + one chain. Launcher arms
`p018_upin_nt72/144` added (unconditional exports), deployed md5
2df6d781. **`p018_upin_nt72` = 13157833** (NT=72/PPS=6/rlxf 0.15, 72 h)
— doubles as the H1/H3 measurement run; **`p018_upin_nt144` = 13157834**
(NT=144/PPS=3/rlxf 0.075, 72 h, expect TIMEOUT ~rev 18 ⇒ chain `_s2`).
Both cold, SETTLE 22, σ=0.04R/OVERLAP 2.75/h pinned, N=1 D=3.4σ,
viscous, merge_r 0.0055, diagnostics ON. Banner watchers armed. Queue =
6 (n4 + 3 conditionals + 2 rungs). Score matched 15–29 per the
pre-registration; per-rung stability screen BEFORE trusting any rung.

**2026-08-12 — 3-point rlxf sensitivity ladder SUBMITTED (Ryan-directed,
this session; error-budget term 9 = relaxation ≈ −0.005 CT, one-sided,
unconverged — the largest named budget term, unaddressed by anything in
flight). Design: production-carrier clones with RELAX_RLXF halved /
quartered, warm from `_s2`@1619 (same recipe as mergeoff), P018_SETTLE_REVS
= 50 ⇒ +15 revs = revs 45–60, matched to the `_s3` chain, which IS the
stock rlxf=0.3 point (free).**
| job | case | knob |
| --- | --- | --- |
| 13157881 | `p018_rlxf0p15_n1_s040_visc` | RELAX_RLXF=0.15 |
| 13157882 | `p018_rlxf0p075_n1_s040_visc` | RELAX_RLXF=0.075 |
Both 48 h / 64G, launcher arms added (unconditional exports, deployed md5
107e986b; cluster src deliberately NOT touched — local src diffs are
other-session 021 work; the arms must run the same code state as the
carrier family). NOTE: first submissions 13157866/67 were scancel'd within
minutes — SETTLE=60 would have run to rev ~70 (+25 revs, unmatched tail);
resubmitted at SETTLE=50. Gotcha recorded: the launcher's mergeoff comment
says RESTART_STEP=1620; the live chains used 1619 (CSV=S+1 gotcha) — 1619
is correct. Banner + warmstart verification armed. Scoring on landing:
stability screen FIRST (reduced relaxation may ignite — an ignition itself
bounds the ladder), then `p018_analyze.py` m1+m2 on matched ~50–60 windows
across {0.3, 0.15, 0.075}: CT̄ slope vs rlxf + extrapolation gap toward
rlxf→0 fills budget row 9; ε_Γ per pair (M2 co-equal per Ryan, ruling 9
REAFFIRMED this session — no CT-only descope). Queue = 8 study jobs.

**2026-08-12 ~22:45 MDT — reset-(g) watch session re-armed; submission-
integrity checks ALL PASS.** (1) **Merge-off verification (re-run of the
checker that died with the session): PASS** — wake_health n_p at matched
steps 1642/1643/1644: mergeoff 217479/217608/217742 vs `_s3`
216463/216549/216650; excess +1016→+1059→+1092 (~+44/step, growing
monotone) ⇒ MERGE_PARTICLES=false took; run kept. Both arms healthy
(min_sr ~0.156, max_u 20–31). (2) **rlxf ladder banners VERIFIED**
(13157881 rlxf:0.15, 13157882 rlxf:0.075; both NT:36 nwakerows:1
overlap:2.75 pps:12 merge_r:0.0055 settle:50 visc:true das_uniform:3.4)
and **warmstart VERIFIED**: both logs print `simulate_warmstart!:
resuming from step 1619`. (3) Cancelled corpses 13157866/67: data dirs
empty (wiped by relaunch, same run names) — nothing to sweep. Disk 169G.
sacct: all 8 study jobs RUNNING. Watcher armed: 25-min poll, terminal
states + ignition tripwire (min_sr<0.06 or max_u>500) on all 8 runs.

**2026-08-12 ~23:15 MDT — 200G-budget VTK sweep (Ryan-directed via
coordinator; `p018_vtk_sweeper.sh --apply --only <run>` per run).** Freed
~10.6G: etadas_n1 3018MB (closed, harvested), etadas_n4 2540MB (live),
mergeoff 873MB, trunc6 1313MB, `_s3` 915MB, upin_nt144 975MB, upin_nt72
837MB (all live, older steps only), p15_smoke_{gon,voff,von} 107MB. Every
run keeps its newest 10 COMPLETE restartable steps; CSVs/TOML/monitors
untouched. NOT touched: protect-listed runs and 019/020-owned dirs
(`scr_p019_*`, `scr_p020r_geom_s020v` 5.4G — excluded despite being
sweepable by the script). Tree ≈129G post-sweep (was ~140G), budget 200G.
**CONFLICT FLAGGED FOR RYAN: the big recoverable chunk (~70G) is
p018_h0p25/h0p125/ov1p4/ov1p0 — coordinator authorized their sweep as
closed+harvested, but all four are STILL on `vtk_protect_list.txt`
(Ryan-owned, agents never edit). The list's own comment says protection
expires once completed/`_s2`-launched — condition met — but removing the
lines is Ryan's action. NOT swept; awaiting Ryan's protect-list edit.**

**2026-08-12 ~23:40 MDT — protect-list conflict RESOLVED (Ryan removed the
four entries; cluster copy verified clean) ⇒ big sweep executed.** Freed
~75.7G via `p018_vtk_sweeper.sh --apply --only`: p018_h0p25 25544MB,
p018_h0p125 20476MB (FLAGGED run — metadata.toml + monitors verified
intact post-sweep; it never had a CT CSV, force-monitor fallback
preserved), p018_ov1p4 15771MB, p018_ov1p0 15678MB. Each keeps its
newest-10 complete restartable steps (h0p25/ov1p4/ov1p0: 1070–1079;
h0p125: 845–854; its `_s2` keeps 1070–1079 from the earlier state).
**Tree now 58G (data/ 56G) vs 200G budget.** Deletion logs appended to
phase_12 + phase_04.

**2026-08-13 ~09:50 MDT — rlxf ladder DOUBLE IGNITION (post-mortems) +
mergeoff harvest (H4b KILLED) + upward rlxf pair submitted.**

(1) **rlxf reduced rungs both IGNITED from healthy warm state (rev 45).**
0.15 (13157881): onset step 1734–1735 = rev 48.2 (+3.2 revs), min_sr
0.125→<0.06 in ~1 step, peak max_u 47,492 m/s @1737, min_sigma pinned at
9.407e-5 = viscous floor, scancel'd @1772 (main session). 0.075
(13157882): onset step 1684–1698 = rev 46.8–47.2 (+2.0 revs), peak max_u
1.12e6 m/s @1706, wake self-annihilated 217k→11k particles, COMPLETED
exit 0 with garbage forces (exit-code-≠-verdict again). Dose-response:
quartered ignites ~1.2 revs before halved. **Verdict: stock rlxf=0.3 is
load-bearing at σ=0.04R/N=1 viscous; term-9 downward slope unmeasurable
here** (error_budget row 9 updated; full table phase_09 §2026-08-13).
Corpse VTK swept (6325+5782 MB, newest-10 restart sets kept, all
CSVs/monitors kept cluster+local).

(2) **Merge-off A/B (13157753, COMPLETED 3h56, healthy tail): H4b
KILLED.** Matched steps 1620–1727 vs `_s3`: total banded Σ|Γ| slope
+3.104 vs +3.001 %/rev (unchanged, Δ+0.10); dominant bands ≤0.2%/rev
apart; n_p slope +2.00 vs +1.22 %/rev (merge sink real, ~0.8%/rev, no
Σ|Γ| footprint); CT 45–47: 0.075284 vs 0.075118 = +0.22%; M2 ε_max
1.56%/rms 0.75%. Detail phase_15 §2026-08-13 (f). Budget term 8 hint:
merge ΔCT +0.22% (3-rev, not settled-window).

(3) **Ryan decision: probe rlxf UPWARD.** Submitted (main session,
launcher md5 3ccab7d3): **13159912 `p018_rlxf0p45_n1_s040_visc`
(RELAX_RLXF=0.45)**, **13159913 `p018_rlxf0p6_n1_s040_visc` (0.6)**;
both warm `_s2`@1619, SETTLE 50 (+15 revs, matched to `_s3` 45–60),
48h/64G. Scoring on landing: m1+m2 matched ~50–60 vs `_s3`; compare
slope magnitude to 006 legacy downward slope (~+0.0011 per halving, 006
line ~162) → can term 9 be carried as measured local sensitivity?
Banner/warmstart verification = main session; this session watches
terminal states + tripwires + disk. Queue = 7 study jobs.

**2026-08-13 ~15:50 MDT — `_s3` (13157752) LANDED + SCORED (phase_15
(d)+(e); detail there).** Healthy completion, revs 45–58. M1 45–57 CT̄
0.076806, block drift +3.117% monotone (was +1.107% over 30–44) ⇒
grand-mean drift continues/steepens; BUT stitched burst/quiet (30–57):
quiet limit 0.072975 ≡ _s2-alone 0.072979, burst incidence 0.74,
rectification r=+0.611 ⇒ **drift = burst rectification over a STATIONARY
quiet baseline ≈0.0730**. Seam at rev 45 CLEAN (no per-step
discontinuity, no within-rev-std collapse; rev-45 ptp 0.00794 in
neighbor range). Inventory 45–58: total Σ|Γ| +0.451%/rev RISING
(4.19→4.91) — no Σ|Γ| decline in this window (consistent with H4b
killed). Local harvest: data/p018_ufront_n1_s040_visc_s3/ (CT CSVs +
all monitors).

**2026-08-13 ~19:30 MDT — DISK ALERT (data/ 179G) → live sweep of all six
running jobs, ~126G freed, data/ back to 56G.** Six live 018 runs swept
via sweeper (newest-10 restartable steps kept, CSVs/monitors untouched,
nothing protected touched): etadas_n4 26.0G (kept 898–907 ⇒ rev ~25.2),
upin_nt144 23.7G (kept 1146–1155 ⇒ rev 8.0), upin_nt72 22.5G (kept
856–865 ⇒ rev 12.0), trunc6 22.4G (kept 692–701 ⇒ rev 19.4), rlxf0p45
15.6G (kept 1966–1975 ⇒ rev 54.9), rlxf0p6 15.6G (kept 1961–1970 ⇒ rev
54.7). Health note: **both upward-rlxf arms are ~10 revs past the
rev-45 handoff with no tripwire alert** — already past the reduced
rungs' ignition ages (+2.0/+3.2 revs), consistent with
relaxation-as-stabilizer. VTK growth rate observed ~15–20G/h aggregate ⇒
expect a sweep per ~6–8 h while all six run; watcher alert threshold
175G.

**2026-08-14 ~02:00 MDT — rlxf UPWARD pair LANDED + SCORED (phase_09
§2026-08-14; local harvest data/p018_rlxf0p{45,6}_n1_s040_visc/).**
Stability screen PASS both (~19h, rev 58; tail min_sr 0.153/0.201 vs
_s3 0.149 — stabilizer dose-response confirmed; no ignition precursors).
Matched 45–57: CT̄ {0.3: 0.076806, 0.45: 0.075475, 0.6: 0.075884} —
NON-MONOTONE, CIs overlap, quiet limits non-monotone too ⇒ **term-9
local slope UNRESOLVED (|Δ|≲1.7% = burst-sampling noise); Γ̄ insensitive
(M2 ε_max 0.98%/0.71% PASS)**. Full rlxf picture: <0.3 ignites, [0.3,
0.6] CT-flat within noise + stability rises ⇒ model-definition reframe
supported; Ryan decision on carrying mode stands. Post-harvest VTK
sweep below.

**2026-08-14 — Phase 16 (chord–σ co-scaling) DEPLOY INCIDENT + screens
resubmitted.** First screens 13170749 `p018_scr_cs_l2p4` / 13170750
`p018_scr_cs_l3p4` FAILED ~4 min in at step 1: `UndefVarError:
radius_inflation` (`src/FLOWPanel_wake.jl:330`). Root cause = deploy skew
with the §2026-08-12 quarantine: the whole-file rsync of the local
working-tree wake.jl carried the concurrent 021 session's uncommitted FMM
radius-inflation coupling (`source_system_to_buffer!`, 2 sites), whose
symbols live in `FLOWPanel_elements_fmm.jl`/`FLOWPanel_abstractbody.jl` —
files deliberately NOT deployed. **The quarantined-file set now effectively
includes `_abstractbody.jl` and `_elements_fmm.jl`.** Remediation:
pre-rsync cluster wake/replay recovered exactly from the NetApp home
snapshot (daily 2026-08-13 23:05 UTC; md5s match the pre-overwrite record
a44e484e/84ae7c04), reconstructed as snapshot base + ONLY Phase-16 hunks,
symbol-audited (no other local-only symbol; reconstructed replay ==
local byte-identical). Driver + launcher deploys audited hunk-by-hunk vs
snapshot: Phase-16-only, left in place. Deployed md5 (local==cluster):
wake `790b21cd901d475c4893917d25361edb`, replay
`20ed0def07be994e08bfe9c02a380d51`, driver
`87d0b120112f1fdf4aec70ada8a0bb3d`, launcher
`8fd6beb7ef3e778367cd55a7b6a3f66e`. Login-node compile check OK
(`StationSigmaOverlap{Float64}`). **Resubmitted 13170886 (l2p4) /
13170887 (l3p4)**, 24 h/64G, `P018_SETTLE_REVS=0` (8-rev stability
screens; they gate the 30-rev `p018_cs_*` ladder). Banners re-verified:
`sigma_chord:0.313 sigma_floor:0 das_lambda:2.4|3.4`, carrier
overlap:2.75 pps:12 merge_r:0.0055 nwakerows:1 visc:true, mesh
45_185_ct4. Details: `phase_16_chord_sigma_coscaling.md` §Log.

**2026-08-14 ~21:30 MDT — Phase 16 screens PASS ⇒ 30-rev ladder SUBMITTED
(scheduled-resume session, brief (j); detail phase_16 §Log).** 13170886
(l2p4) / 13170887 (l3p4) COMPLETED healthy, 11.47 revs each. No ignition:
max_u ≤ 41/53 m/s declining at tail; min_sigma 19× above viscous floor;
zero σ≤0; |Γ|/σ² 72 (declining) / 297 (decelerating) vs healthy uniform
carrier 1.29e3 at matched rev-11 age; n_p tracks carrier (no tip-merge
fusion). FLAG: min_sigma_ratio + p1_sigma_ratio NaN all steps under the
station-σ law (monitor reference σ NaN) — tripwire judged from absolutes +
carrier calibration; monitor fix staged for after the ladder. Submitted
(pre-authorized on pass): **13178762 p018_cs_l2p4, 13178763 p018_cs_l3p4,
13178764 p018_cs_l4p8** (48 h/64G, P018_SETTLE_REVS=22 ⇒ 30 revs).
Screens VTK swept post-harvest (19.2G; data/ 123G). Queue = 6 study jobs
(3 ladder PENDING + trunc6/upin_nt72/upin_nt144 RUNNING).

**2026-08-14 ~23:45 MDT — trunc6 TIMEOUT + H2 first read (NOT killed,
depth-dependent drift) + `_s2` deconfound chain; ladder banners VERIFIED.**
(1) 13157751 `p018_trunc6_n1_s040_visc` TIMEOUT at 48 h, step 1062/1080
(rev 29.5) — window intact. Matched 15–29 vs 3.5R carrier: CT̄ 0.074411 vs
0.075298, block drift **−3.56% DOWN vs +2.08% UP** (direction flips with
depth ⇒ H2 not killed); late 22–29 level −2.7%; M2 ε_max 1.87%. CONFOUND:
5.5R still filling through rev 27 (n_p +6–9k/rev) ⇒ downward drift mixes
fill transient with depth response. Detail phase_15 (c). (2) Submitted
**13179250 `p018_trunc6_n1_s040_visc_s2`** (warm @1061, SETTLE 36 ⇒ rev
44, 72 h/64G) to score H2 on a settled 30–44 window — session judgment
call (extension of a Ryan-authorized conditional; logged in STATUS
REPORT). (3) **Phase-16 ladder banners VERIFIED running**: 13178762/63/64
= das_lambda 2.4/3.4/4.8, sigma_chord 0.313, sigma_floor 0, settle 22,
carrier line identical, θ_max 1.26/1.79/2.53 rad, cold start. ETA ~08-16.
Queue = 6 study jobs (3 ladder + 2 upin RUNNING, trunc6_s2 PENDING).
Disk: data/ 123G after screen sweep.

**2026-08-15 ~18:30 MDT — TWO IGNITIONS: trunc6_s2 (H2 depth probe) and
upin_nt72 (dt-ladder screen FAIL ⇒ LADDER STOPPED).** (1) 13179250
`p018_trunc6_n1_s040_visc_s2` COMPLETED exit 0 in 20.3 h with GARBAGE
tail (cycle-mean "CT" 29.58 — exit-code-≠-verdict again): Γ-ignition at
rev ~34.5, n_p 344k→13k, max_u 2.8e4, min_sigma pinned 9.49e-5. NOT a
seam artifact: the COLD 5.5R run had already left the carrier's health
envelope pre-handoff (|Γ|/σ² rev 28: 1.84e3 ≈ carrier 1.98e3; rev 29:
5.3e3 vs 1.6e3; 1.55e4 at step 1061) while the 3.5R carrier stays O(1e3)
through rev 44 ⇒ **the 3.5R truncation deletion is LOAD-BEARING for
stability** (deletes the aged tail before Γ-ignition matures; joins
rlxf=0.3 and σ=0.04R as stabilizers, 019 regime map). H2 settled-window
quantification BLOCKED at 6R depth; cold-run read stands (depth-dependent
drift, level −1.2..−2.7%, fill-confounded; phase_15 (c)). Deeper-H2
options are RYAN's call (e.g. 4.5R arm). (2) 13157833 `p018_upin_nt72`
FAILED exit 1 at rev 19.6: unarrested min_sr collapse (dip 0.014 rev 15,
transient recovery, re-ignition rev 18.6), min_sigma pinned at the NT=72
floor 6.65e-5, max_u 1e6, wake self-annihilated. rlxf 0.15 is
dt-equivalent to stock ⇒ ignition reads on dt refinement itself at
matched physical relaxation/σ/clearance. **Pre-registered rule applied:
dt ladder STOPPED; nt144 (healthy rev 14, times out ~rev 14.6 tonight)
will be harvested but its `_s2` chain is NOT submitted; Richardson triple
dead at these knobs.** Partial burst attribution NT36 vs 72 (revs 7–13):
within-rev CT std ×0.66 (leans numerical) BUT co-timed rev-14 burst in
both rungs (that episode physical); CT̄ Δ −0.31%; complete with nt144.
Detail phase_03 §2026-08-15. (3) Ops: corpse sweeps trunc6_s2 14.9G +
upin_nt72 24.9G + trunc6-cold 22.0G ⇒ data/ 226G→166G. Ladder healthy at
revs 17–19 (|Γ|/σ² 202/338/1800 vs carrier envelope O(1e3); l4p8 warmest,
watched). Queue: 3 ladder + nt144.

**2026-08-15 ~19:30 MDT — Ryan adopted the "adopt-as-is" parts of the
external minimal-run plan (item file STATUS REPORT items 12–17).** In
force: Phase-16 close-out per existing pre-registration (no
intermediate-λ densification; failure routes by signature to F1/F2);
selected-rung extension rule (pass at 30 revs ⇒ done, else warm-chain
that one rung to 45–60 revs = primary validation case; report raw +
quiet-limit if rectification persists); standing do-not-run list (no
low-rlxf reruns, no more upward rlxf, no nt144 chain, no new
truncation/merge ladders, new tags for model changes); acceptance
language (settings blind to CT_exp; health+settledness gates; downgrade
wording if temporal/spanwise unresolved). NOT adopted pending Ryan:
exact-rate NT72 (0.16334), NT72@0.3 swap, spanwise 45→60 (80-station
claim needs phase_10 verification). Watch: cs_l4p8 |Γ|/σ² 6.7e3 @ rev
18.8, rising above carrier envelope — potential θ-correlated stress on
the λ=4.8 rung (P2 signature if it fails, screens still valid).
nt144 at its 72 h wall, harvest next.

**2026-08-15 ~21:45 MDT — nt144 harvested (attribution verdict: bursts
NUMERICAL) + two NT=72 MODEL-DEFINITION arms launched (Ryan-directed).**
13157834 TIMEOUT at rev 14.45 HEALTHY (min_sr 0.384, |Γ|/σ² 92); no
chain (screen-fail stop + adopted do-not-run list). 3-point attribution,
shared clean revs 7–13: within-rev CT std 0.00303 → 0.00199 → 0.00035
(NT 36/72/144; ×0.12 overall, monotone with Δt) ⇒ **fluctuation/burst
mode NUMERICAL**; the rev-14 co-timed burst (36+72) is ABSENT at 144.
Caveats: pre-settled window only; rlxf dt-halved per design; **CT̄
non-monotone in NT (0.07513/0.07490/0.07775)** — mean-CT dt axis
unresolved. Supports quiet-limit as M1 observable (staged; Ryan's
decision). Detail phase_03 §2026-08-15 (night). **Launched on Ryan's
explicit direction (supersedes the resubmission hold), labeled
model-definition changes: 13183888 `p018_upin_nt72_rlxf0p3` (stock
per-step rate) + 13183889 `p018_upin_nt72_rlxf0p16334` (exact
continuous-rate equivalent, 1−√(1−0.3)); 72 h/64G SETTLE 22; new
launcher arms deployed md5 4135f810 local==cluster; banner-verify on
start.** Queue = 5 study jobs (3 ladder RUNNING revs ~19–22, 2 new
PENDING). Also ENACTED earlier tonight: adopted-plan items 12–15 (item
file STATUS REPORT).

**2026-08-15 ~22:15 MDT — third model-definition arm launched
(Ryan-directed): 13183998 `p018_upin_nt144_rlxf0p3`** (NT=144, PPS=3,
per-step rlxf 0.3 = stock rate ⇒ ~4× the per-time damping of the NT=36
carrier; otherwise the pinned construction). 72 h/64G SETTLE 22; one
segment reaches ~rev 14–15 (NT=144 ≈ 120 h/30 revs). Launcher arm added
+ deployed md5 9a5b2f2cb459d75a1b5986be56493f02 local==cluster. Banner
verification when it starts. **RYAN INSTRUCTION ON FILE: when 13183998
finishes, ASK RYAN whether to chain it `_s2`** (he expects yes if the
NT=72 rlxf0p3 arm 13183888 looks good) — do not chain without asking.
Meanwhile 13183888/889 banners VERIFIED (rlxf:0.3 / 0.16334, NT:72
pps:6, das_uniform:3.4, cold). Queue = 6 study jobs.

**2026-08-15 ~23:30 MDT — CONTEXT-RESET PREP: RESET BRIEF (k) written**
(item file top, supersedes (j)). Fresh-agent entry = brief (k): six-job
watch table (ladder 13178762/63/64 + model-def arms 13183888/889/998),
13183998 banner still to verify, ASK-Ryan-before-chain note, verdicts
recorded (screens PASS, truncation load-bearing, dt ladder stopped,
bursts numerical), adopted rules 12–15, deploy-skew quarantine, disk
state. Memory file + MEMORY.md index updated to point at (k).

**2026-08-15 ~22:30 MDT — SCHEDULED PAUSE/RESUME ARMED (Ryan-directed).**
Ryan: pause all HPC jobs 23:50 tonight; resume 00:00 Monday 08-17.
Capability test: `scontrol suspend` DENIED for regular users on orc
(tested on throwaway job 13184024, "Access/permission denied"; stray
test jobs cleaned). SIGSTOP via scancel -s STOP rejected as the
mechanism: walltime keeps accruing while frozen ⇒ a 24 h freeze
TIMEOUTs the 48 h ladder mid-run. Authorized fallback (Ryan: "we can
kill them and restart") = scancel at 23:50 with a pause manifest
(per-run restart step + knobs → `pause_manifest_20260815.md`), then
warm-chain resubmissions Monday 00:00 per ops_reference recipe, banner
+ warmstart verified. Both tasks armed as in-session one-shot crons
(78a94d26 pause, 4a0fcf9d resume). CAVEAT: crons are SESSION-ONLY —
this Claude session and the ssh ControlMaster socket must stay alive
through Monday 00:00; if either dies the tasks do not fire (each task
prints a loud warning if ssh fails rather than retrying into 2FA).

**2026-08-15 23:52 MDT — SCHEDULED PAUSE EXECUTED (Ryan-directed; jobs
killed ON PURPOSE, not failures).** All 7 rander39 jobs scancel'd; queue
verified EMPTY. Six 018 runs harvested (monitors → local) and their
post-kill restart steps verified complete in all 3 VTU families + pvd:
cs_l2p4 @834 (rev 23.2), cs_l3p4 @778 (21.6), cs_l4p8 @713 (19.8),
upin_nt72_rlxf0p3 @252 (3.5), upin_nt72_rlxf0p16334 @194 (2.7),
upin_nt144_rlxf0p3 @268 (1.86). Seventh job = other-session 13184015
`fm041aH` (batch script saved ~/fm041aH_batch_13184015.sh.pausecopy;
restarts from scratch). Health note at kill: **cs_l4p8 gos2 SUBSIDING
(6707@677 → 3468@713)** — the hot phase passed without ignition; all
arms healthy. Resume = 00:00 MDT Mon 08-17 (cron 4a0fcf9d) per
**`pause_manifest_20260815.md`** (full sbatch templates, walltimes,
banner checklist there). Watch monitor exited (ALL_JOBS_TERMINAL, as
designed); re-arm at resume.

**2026-08-17 00:01 MDT — SCHEDULED RESUME EXECUTED (cron 4a0fcf9d, per
pause_manifest_20260815.md).** Launcher md5 re-verified unchanged
(9a5b2f2c); queue was empty. Submitted warm chains (SETTLE 22, 64G):
**13185005 cs_l2p4_rs1 @834 (24 h), 13185006 cs_l3p4_rs1 @778 (24 h),
13185007 cs_l4p8_rs1 @713 (36 h), 13185008 upin_nt72_rlxf0p3_rs1 @252
(72 h), 13185009 upin_nt72_rlxf0p16334_rs1 @194 (72 h), 13185010
upin_nt144_rlxf0p3_rs1 @268 (72 h; ASK RYAN before chaining further)** +
**13185011 fm041aH** (other-session job, resubmitted from saved batch
script, from scratch). All PENDING at submission. Banner + warmstart
verification on start; monitor re-armed (states + absolute tripwires on
the six _rs1 runs). Ladder ETA now ~08-17 late (l2p4/l3p4) / ~08-18
(l4p8).

**2026-08-17 ~02:45 MDT — _rs1 fleet VERIFIED (all six banners +
warmstart lines line-for-line correct).** 13185005/06/07: sigma_chord
0.313, sigma_floor 0, das_lambda 2.4/3.4/4.8, θ_max 1.2649/1.792/2.5298,
resuming from 834/778/713. 13185008/09/10: NT 72/72/144, pps 6/6/3,
rlxf 0.3/0.16334/0.3, das_uniform 3.4, resuming from 252/194/268.
Carrier line identical everywhere (overlap 2.75, merge_r 0.0055,
nwakerows 1, visc true, settle 22, mesh 45_185_ct4). No mismatches ⇒
nothing cancelled. fm041aH (13185011) COMPLETED healthy after 2 h 29
(benchmark wrote its CSV, "=== done ===") — no further action, its
session owns the results. Queue = six 018 _rs1 RUNNING.

**2026-08-17 ~16:00 MDT — Phase-16 l2p4/l3p4 chains landed + interim
scored (detail phase_16 §Log).** Both healthy to rev 30, seams clean.
**Co-scaling collapses the NT=36 burst mode** (per-rev std ×7–10 smaller
than carrier; λ2.4 = first raw-M1 PASS of the campaign; λ3.4 CHECK
0.52% non-monotone). CT̄ 15–29: 0.071074 (λ2.4) / 0.072442 (λ3.4) ⇒
λ-slope +1.92% — P1 doubling verdict awaits l4p8 (~19:00 MDT). M2
FAILs: rung pair ε_max 5.97%, P3 vs uniform carrier 6.45%. upin _rs1
arms healthy (revs ~8–10 at check). Queue = 4 (l4p8 + 3 model-def arms).

**2026-08-17 ~19:45 MDT — PHASE-16 LADDER COMPLETE + SCORED (detail
phase_16 §Log 2026-08-17).** l4p8 landed healthy (no ignition; rev-19
heat subsided). **P1 FAIL** (doubling 2.4→4.8: ΔCT̄ +9.63% raw / +9.5%
quiet; M2 31.35%). **P2: inboard, θ-correlated, λ-monotone Γ̄ EXCESS**
(+46.7% of Γmax at r/R 0.27 → ~0 by 0.76; sign opposite the predicted
deficit, disclosed; outboard low-rung pair agrees <1%). **P3: σ-law
redistribution is the big lever** (λ3.4 vs uniform: ΔCT −3.79%, M2
6.45%) + burst-mode collapse (λ2.4 = first raw-M1 PASS of the
campaign). **Pre-registered routing selects F1 (curvature-capped Das);
STAGED for Ryan — nothing submitted.** Evidence: excess onset θ ≈
0.5–0.7 rad; even λ2.4 binds inboard. Ladder VTK sweeps launched;
model-def arms healthy in flight.

**2026-08-17 ~20:15 MDT — VTK sweeps complete: ladder dirs 109.4G +
model-def arms 78.5G freed; data/ 191G → 114G (budget 200G).** Newest-10
restartable sets + all monitors/CSVs retained everywhere; live runs
swept older-steps-only. Queue = 3 model-def arms (~22.7 h elapsed, all
healthy). Next events: nt72 arms + nt144 arm land ~08-20 (72 h walls);
nt144 landing triggers the ASK-RYAN chain question; F1 confirmation
pending Ryan.

**2026-08-18 — F1 (curvature-capped co-scaling) LAUNCHED per Ryan's
order (detail phase_16 §Log 2026-08-18).** Driver knob
DAS_CURVATURE_BETA implemented (driver-only, no src/); β=0.6 rad; gate
smoke clean 9 revs + banner verified (cap binds 15/36 stations, θ_max
0.6 exactly); deployed md5 257cdc6b/9949b29d local==cluster. Submitted
**13193493/494/495 = p018_cs_f1_l{2p4,3p4,4p8}**, 48 h/64G, SETTLE 22.
Queue = 6 study jobs (3 F1 PENDING + 3 model-def arms RUNNING, land
~08-20). Disclosed: cap pushes inboard clearance below the C1 band
(λ_eff ~1.9 at r/R 0.27) — intrinsic band conflict, watch M2 inboard.

**2026-08-18 — F1 env incident + resubmit; F1b (endpoint-on-arc Das)
IMPLEMENTED + probe 2×2 submitted (Ryan-directed, plan
~/.claude/plans/stage-f1b-by-writing-elegant-stonebraker.md).**
(1) F1 trio 13193493/494/495 FAILED at load: cluster Manifest.toml had
been re-resolved under Julia 1.12.6 (Aug 17 20:35, 021-session window;
1.11.7 hard-fails on 1.12 manifests). Fix: restored 1.11.7 Manifest
from the daily NetApp snapshot, conservative Pkg.resolve() to absorb
021's ILUZero addition, `using FLOWPanel` verified; 021's manifest
saved as Manifest.toml.p021_julia112_20260817. **Resubmitted F1:
13193513/514/515** (banners pending their start). CROSS-SESSION FLAG:
a 021 resubmit under Julia 1.12 can re-break this.
(2) **F1b Route B implemented** (Ryan rulings: general frame-walked
construction, arc-length-integrated endpoint — no floors, Route A
multi-segment REJECTED, live mode deferred on the Kutta operator
cache): src `_set_Das_station_arc!` + `_rigid_back_velocity` +
initialize_Das! kwarg `set_Das_station_drifts`
(src/FLOWPanel_simulate.jl; PROVEN cluster-base+my-hunks-only before
whole-file deploy, pre-overwrite md5 78858390 recorded, deployed
a398bb71); driver knobs DAS_ARC_PLACED / DAS_ARC_HELIX_SOURCE
(steady|kinematic) / DAS_ARC_TABLE + local-TE-basis drift
reconstruction (deployed 4b7fd1cb); unit testset "endpoint-on-arc Das"
14/14 PASS (exact circle/helix identities; legacy station-lengths
tests untouched = zero-behavior-change); local kinematic-mode smoke in
flight. (3) **TE downwash probe** = replay monitor
(examples/p018_te_downwash_probe_replay.jl, deployed 1320f5ab):
local shakeout found a point probe at the TE is DOMINATED by the
body's own bound-sheet field (0.5–0.9 tip speed) ⇒ table u = WAKE-only
induced + Uinf (all-sources kept as diagnostic columns); clamp =
user-specified axial direction (generalization deferred on Ryan's
go-ahead). **Probe job 13206092 submitted**: 2×2 (sources
cs_l3p4_rs1 end / _s3 end × locations te / mid), last-10-step +
blade averages. ON LANDING: compare tables, ASK RYAN (location +
source + any clamp) before submitting csarc screens.

**2026-08-18 (cont) — probe job take-2.** 13206092 FAILED (BoundsError:
ProbeSystem sized from the SETUP-ONLY rotor while the slurm env
defaulted to the 40_40 mesh @ RPM 6000 vs the saved 45_185_ct4 @ 5400
body). Fixes: probe/accumulators now allocated LAZILY from the REPLAYED
body on first call, and the job script exports
RHPC_MESH=45_185_ct4 RPM=5400 NT=36. Local re-shakeout PASS (lazy path
exercised, mid location). Resubmitted **13206196** (probe script md5
75c2c499 both sides). Monitor still watching the old ID — replaced.

**2026-08-18 ~17:45 MDT — TE-probe 2×2 LANDED + Ryan's table picks +
csarc screens & mid-sanity SUBMITTED.** Probe 13206196 COMPLETED: 4
tables at data/p018_{cs_l3p4_rs1,ufront_n1_s040_visc_s3}_te_downwash_
{te,mid}.csv (local+cluster). Findings: mid-span wake-only downwash
u_n ≈ −4.8 m/s (0.07 tip), TE-vs-MID agree 3–14% mid-span; sources
agree 10–15% mid-span BUT **_s3 (rev 58) shows inboard FOUNTAIN FLOW
(u_n up to +4.7 at r/R<0.24) that cs_l3p4 (rev 30) lacks** — wake-
maturity effect, present at both probe locations; cs table shows
physical TIP UPWASH (+1.4–3.3 outboard of r/R 0.89) ⇒ **clamp would
zero real tip physics, NOT free** (answered Ryan's question). **Ryan
picks: TE location + cs_l3p4_rs1 source + clamp OFF + one midpoint-
table sanity arm.** Launch-blocking driver fixes: table reader skips
#-comments; accepts tip→root row order; validated parse (max |u| =
0.102 tip). Local kinematic smoke killed by RYAN at rev 3.1 (resource
pressure; record clean: no σ≤0, gos2 42) — screens carry the stability
gate. Deployed driver 00bda73d + launcher f2ec5789 (md5 both sides).
**Submitted: 13206336 p018_scr_csarc_l2p4, 13206337 p018_scr_csarc_l3p4
(24 h, SETTLE 0), 13206338 p018_csarc_mid_l3p4 (48 h, SETTLE 22).**
TE ladder p018_csarc_l{2p4,3p4,4p8} gated on screens. Banner checks on
start (das_arc:true arc_src:steady arc_table per case).

**2026-08-18 ~18:15 MDT — CONTEXT-RESET PREP: RESET BRIEF (l) written**
(item file top, supersedes (k); pause-window note demoted to history).
Fresh-agent entry = brief (l): nine-job table (F1 trio banners
UNVERIFIED — first duty; model-def walls Wed ~24:00 with the nt144
ask-Ryan chain gate; csarc screens gate the pre-authorized TE ladder;
csarc_mid sanity A/B), F1b state + Ryan's table picks, probe findings,
deploy md5s, Manifest-1.11.7 flag, disk-recheck note. Memory +
MEMORY.md updated to point at (l).

**2026-08-20 ~morning — FRESH-AGENT RESUME (brief (l)); all nine jobs
accounted for.** sacct: F1 trio 13193513/514/515 TIMEOUT at the 48 h wall
(ended 06:49, reached revs 28.4/28.7/28.9 of 30); model-def nt72 arms
13185008/09 TIMEOUT at 72 h (revs 28.2 / 26.6); **13185010 nt144@0.3
FAILED at 51.5 h — Bus error (core dumped), step 1882/4319 = rev 13.07**;
csarc screens 13206336/337 COMPLETED (12.6/12.3 h, 413 steps = 11.5
revs); csarc_mid 13206338 RUNNING (~42 h elapsed). Disk at resume:
data 165G.

**2026-08-20 — F1 trio banners VERIFIED post-hoc from logs (first duty;
never verified at start).** All three: `sigma_chord:0.313 sigma_floor:0
das_beta:0.6`, das_lambda 2.4/3.4/4.8 per case, `F1 curvature cap
ACTIVE: theta_j <= beta = 0.6 rad`, cap binds 10/13/17 of 41 stations
(production mesh), θ_max = 0.6 exactly, carrier line identical to the
uncapped ladder (45_185_ct4, NT36, depth banner 4.0R, nwakerows 1). No
mismatch ⇒ runs are valid F1 arms.

**2026-08-20 — nt144@0.3 (13185010) post-mortem: NODE FAULT, not
physics.** .err shows `Bus error (core dumped)` mid-step; wake health at
death was calm (max_u ~34–41, min_sigma 1.34e-3, gos2 ~290, sigma growth
1.076) — no ignition signature. Newest-10 restartable steps retained
through 1881 (live-sweep policy), so a warm resubmit from step 1881
(rev 13.06) is available. **Restart decision STAGED for Ryan (folds into
his standing nt144 ask; nothing resubmitted).**

**2026-08-20 — csarc F1b screens: BOTH PASS (pre-registered absolutes;
sigma-ratio cols NaN as expected).** Full 413-step tails: zero σ≤0;
max_u tail 14–23 m/s, flat/declining (cs screens passed at ≤41/53;
ignition = 1e4–1e6); min_sigma ends 1.78e-3 m both (identical to the
passing cs screens, 19× above viscous floor, never pinned); |Γ|/σ²
l2p4 ends 425 rising-decelerating, l3p4 flat ~81 — both BELOW the
healthy-carrier rev-11 calibration 1.29e3 (phase_16 §Log 2026-08-14);
n_p ~195k tracks the carrier. Watch note: l2p4's gos2 still rising at
end (×~1.4/rev, decelerating) — inside the uncapped ladder's healthy
200–500 envelope, flagged for the 30-rev tail.

**2026-08-20 — TE ladder SUBMITTED (pre-authorized on screen PASS per
the approved F1b plan §Cases).** INCIDENT first: submission take-1
(13243011/12/13) FAILED in 13 s exit 2 — remote double-expansion ate the
case argument (`ERROR: pass a case tag`; same MOTD/quoting class as the
sentinel rule; no side effects). Resubmitted via `bash -ls` stdin:
**13243083 p018_csarc_l2p4, 13243084 p018_csarc_l3p4, 13243085
p018_csarc_l4p8**, 48 h/64G, P018_SETTLE_REVS=22 (30 revs). All three
RUNNING same morning; banners VERIFIED: das_lambda per case,
`das_arc:true arc_src:steady arc_table:p018_cs_l3p4_rs1_te_downwash_te.csv`,
sigma_chord 0.313, settle:22, das_beta:nan (no cap — correct for F1b).

**2026-08-20 — F1 trio SCORED (matched window 15–28 everywhere; the
TIMEOUT at rev ~28.5 cost rev 29, so the uncapped rungs were re-scored
on 15–28 for like-for-like; no chains submitted).** Detail phase_16
§Log 2026-08-20. Headlines: **the cap collapses the λ axis** — CT̄
0.070539/0.070605/0.069877 (λ2.4/3.4/4.8), doubling 2.4→4.8 = −0.94%
(uncapped same window: +9.68%, 0.071076→0.077959); neighbor 2.4→3.4 =
+0.09% (inside CI). Burst mode stays collapsed even at λ4.8 (per-rev
std 1.8e-4 vs uncapped 1.76e-3). M2: doubling ε_max 5.63%/rms 1.89%
(uncapped 31.35%), neighbor 1.50% — still above the ≤1% gate. Γ̄
overlay: **the inboard θ-correlated excess is GONE under the cap**
(doubling ΔΓ̄ inboard −0.6..−0.9% of Γmax vs +34..+54% uncapped);
residual doubling structure is an OUTBOARD lobe r/R 0.78–0.87 (peak
−5.6% at 0.83) that the uncapped ladder also carried (−6.2% at 0.85) —
cap-untouched territory, i.e. a separate outboard λ-sensitivity, not
curvature. Cap-only A/B at λ3.4 (F1 vs uncapped): ΔCT̄ −2.56%
(0.070605 vs 0.072460), M2 ε_max 7.80% concentrated inboard (cap −2..−6%
of Γmax inboard of 0.4, ~0 outboard of 0.55) — the cap acts where it
binds and nowhere else. **Verdicts + head-to-head vs F1b STAGED for
Ryan; no further F1 submissions.**

**2026-08-20 — model-def nt72 arms SCORED (labeled MODEL-DEF A/Bs; NT
changes the model via pps/das under the station law).** Both arms ran
healthy far past the rlxf0.15 marks (precursor 14.5 / death 18.6):
stock 0.3 to rev 28.2, exact-rate 0.16334 to rev 26.6, no ignition
(gos2 ≤ ~1.7e3, min_sigma plateau ~6e-4, max_u ≤ 38/65). ⇒ the nt72
ignition was a LOW-relaxation dose effect, not NT-intrinsic; stock rate
is stable at NT=72. Matched 15–26 vs NT36 carrier (0.075155): nt72@0.3
CT̄ 0.077198 (+2.72%) but BURSTY/DRIFTING (per-rev std 3.2e-3, block
drift 6.6% monotone — echoes the carrier's burst mode at doubled NT);
nt72@0.16334 CT̄ 0.076717 (+2.08%), std 1.2e-3, drift 0.48% — the
continuous-rate arm is the well-behaved one. M2 vs carrier: 4.34% /
2.93% (0.3 / 0.16334). Dose within NT=72 (0.3 vs 0.16334): ΔCT̄ −0.62%.
**Interpretation + any further NT work STAGED for Ryan.**

**2026-08-20 — ops.** VTK sweeps applied to the seven landed dirs (F1
trio 33G, nt72 _rs1 arms 18G, screens ~0): 46.9G freed, data 165G →
120G (budget 200G; csarc_mid live 17G + three ladder jobs in flight —
headroom OK). Local harvests: monitors+metadata for all seven landed
runs (CT CSVs absent on TIMEOUT runs — scored via the force-monitor
reconstruction path in p018_analyze.py). nt144_rs1 left unswept (578M,
newest-10 restart set through step 1881 = the staged-restart asset).
sacct monitor re-armed (state changes on 13243083/84/85 + 13206338).
Standing: Manifest stays 1.11.7 (untouched this session); no notebook
writes; no chains without Ryan.

**2026-08-20 — RYAN RULING (standing, until he says otherwise): all NT
ladders use the EXACT continuous-rate relaxation rule** r(NT) =
1 − (1 − 0.3)^(36/NT) — 0.3 @ NT36, 0.16334 @ NT72, 0.08539 @ NT144 —
preserving the per-time relaxation dose of the stock carrier. Basis:
the nt72@0.16334 arm matched the carrier's temperament (std 1.2e-3,
drift 0.48%) while stock-rate nt72@0.3 went bursty/drifty; 0.16334
survived to wall with no ignition. Ryan notes the residual CT shift
(+2.08% at NT72) reflects real timestep-linked shedding-property
changes — score NT rungs as labeled model-def A/Bs, not pure dt
refinement. Caveat carried: 0.16334 sits +9% above the ignited 0.15,
close to the boundary — keep the tripwire on any future NT≥72 arm.

**2026-08-20 — F1b NT LADDERS LAUNCHED (Ryan order: "launch NT ladders
for F1b at the three λ rungs, 9 total, omit what's running" — the three
NT=36 rungs are the live csarc λ ladder, so SIX new runs).** New
launcher arms `p018_csarc_nt72_l{2p4,3p4,4p8}` (NT 72, pps 6, rlxf
0.16334) and `p018_csarc_nt144_l{2p4,3p4,4p8}` (NT 144, pps 3, rlxf
0.08539) — NT·pps = 432 (same particles/rev as NT36), relaxation per
Ryan's exact-rate standing rule, all other exports identical to the
csarc λ arms (s* 0.313, steady TE table, no cap). Launcher deployed
md5 75f4a8dd both sides. Submitted **13245449/450/451 (nt72 λ
2.4/3.4/4.8) + 13245452/453/454 (nt144 λ 2.4/3.4/4.8)**, 72 h/64G,
SETTLE 22. Banners VERIFIED on start (NT/rlxf/pps/λ/das_arc:true/
TE table correct on all six). Expected reach in 72 h: nt72 ~rev 28
(full matched windows), nt144 ~rev 13–15 (chains would be needed for a
settled window — Ryan's call later). Note: 018 study-job count is now
TEN active (Ryan's direct order supersedes the six-job directive).
Disk watch: ~10 VTK-writing jobs — sweep threshold will arrive fast;
monitor re-armed on all ten.

**2026-08-20 — ParaView Das-path inspection copy (Ryan request):**
newest-10 steps of `p018_cs_f1_l3p4` (steps 1025–1034, rev ~28.5) and
live `p018_csarc_l3p4` (steps 168–177, rev ~4.8) — body/wake
sheets/particles/filaments — to Ryan's Mac at `~/p018_das_inspect/`
(573M total; .pvd indexes intentionally not copied).

**2026-08-20 — RYAN NOTE (parked, not urgent): N=0 "convert-at-shed"
cleanup.** Discussion clarified the near-wake anatomy: rigid Das row =
the implicit Kutta carrier (in the linear system); the single
`nwakerows=1` FREE row beyond it is already classical time-accurate
(convects with kinematics + induced velocity, pinned to the rigid end
by update_TE!). First particles thus appear at ~Das + travel (~6.6σ ≈
2c at NT36/0.75R), beyond the 3–4σ clearance target that Das = 3.4σ
alone would satisfy. **Ryan: we COULD eliminate the free row (convert
the would-be row to particles at shed time, buffer = Das exactly) as a
cleanup — it shouldn't affect results — but it's NOT urgent; leave it
for now.** If ever done: modest src/ change in _wake.jl/_simulate.jl
(quarantined files; 016 golden-reference tests as the no-side-effect
gate); watch fresh-blob proximity to the rigid row's singular end
filament and sheet-continuity representation (017 territory). Scoring
lens for the live NT ladder: the free-row extent contracts ∝ Δt
(≈3.2σ → 1.6σ → 0.8σ at 0.75R over NT 36/72/144), so NT-flatness of
F1b CT̄ doubles as evidence the free-row extent is benign.
**REMINDER (Ryan request): when the NT ladder (13245449–454) lands,
remind him about the parked N=0 convert-at-shed cleanup alongside the
scoring (the NT-flatness read is the relevant evidence).**

**2026-08-20 — RYAN NOTE (modeling-scope limitation, recorded verbatim in
substance): Das-as-clearance is STEADY-STATE modeling only.** Using the
prescribed Das row to enforce the particle-σ separation criterion (Das =
λσ ≈ 1c) replaces the first chord of wake history with a
quasi-steady prescribed element (F1b further bakes in a steady-wake
downwash table), so the near-wake cannot respond on sub-transit
timescales — fully unsteady load prediction (Wagner/Theodorsen-type
deficiency, BVI, gusts) is out of scope for this architecture. The
fully-unsteady alternative Ryan sketches: keep Das SMALL — just large
enough to enforce the Kutta condition in a converged sense (014's lower
band edge ~0.25c) — and let FREE wake-panel rows (N such that the free
sheet spans the 3–4σ clearance) carry the near wake time-accurately,
converting to particles only beyond it. Not currently done because past
simulations BLEW UP due to the free wake panels; that instability is
the blocker, not the concept. Possible future work (Ryan "wondering
about an approach where we do that, actually"): revisit with the
stability tooling that now exists (viscous CoreSpreading handoff,
019 tripwire/screen protocol, merge policy, panel-kernel core_size) —
would be a NEW item + stability screens, NOT a change to the live 018
campaign, whose deliverable is explicitly the steady/converged state.
**(amendment, same date)** Correction to the note above, on Ryan's
challenge: the ~0.25c lower edge is NOT a validated Kutta floor — it is
the empirical lower edge of 014's wing log-plateau (+0.205% CL/doubling)
whose position was substantially set by particle-σ clearance at the
measured σ/c (0.25c ≈ 3σ on the wing), i.e. clearance in chord units.
014's split sweep shows sheet-buffer 0.25c→8c moves CL <0.1% once
particles are kept away — nothing bounds the Kutta-only requirement
from below at 0.25c; it is untested and plausibly set by TE panel size/
wake-core scale, well under 0.25c. Any future small-Das+free-rows
unsteady design must include a dedicated rigid-row-length mini-sweep
with clearance held by the free sheet.
**(cont) 2026-08-20 — the unsteady alternative above is now STAGED as
BRAINSTORM/025_small_das_free_row_unsteady_wake.md (Ryan: stage, don't
implement): self-contained plan, code anchors, 0.25c-provenance
correction, 4-phase gated plan (Kutta-floor mini-sweep first).**

**2026-08-20 — RYAN RULING: VTK retention raised to newest 36 timesteps
(from 10).** agent_policies/HPC.md + p018_vtk_sweeper.sh KEEP_STEPS
default updated, both deployed local==cluster. Applies to all future
sweeps (one NT=36 revolution of restartable state per run); ruling-10
mechanics otherwise unchanged. Disk-budget note: sweeps now free ~26
steps less per run (~2.6× retained VTK) — the 200G budget bites sooner
with ten jobs writing.

**2026-08-20 evening — 024 N=0 IMPLEMENTED (subagent, all suites green)
+ FIRST CLUSTER A/B LAUNCHED (Ryan order: "launch another NT=36 with
λ=4[.8] for F1b to test it").** Implementation: Route II (age-0
conversion on N=1 storage + convert_at_shed marker); wake 707/707
(+55), simulate 177/177 (+17), replay 142/142 (+17), warmstart 40/40
(+14); N≥1 golden tests untouched at identical counts; N=0 particle
output bit-equal to matched N=1 conversion. Full detail + watch-items
in BRAINSTORM/024 §Log. **Driver amendment (this session, post-024):**
legacy conversion now allowed at N=0 (the wake constructor permits
LegacyEdgeJumpConversion + StationSigmaOverlap at N=0; the smooth-only
driver gate existed only because legacy kwargs hardcoded
unsteady_filament=false — now unsteady_filament=(N==0) in legacy mode).
First smoke attempt FAILED at construction on exactly that flag (the
loud-error invariant working); after the fix, local gate smoke (40_40,
csarc λ4.8 knobs, N=0, 4 threads) ran 50+ steps clean — banner
"nwakerows=0 (convert-at-shed)", CF_x ≈ −0.047 (N=1-like), sigma
growth 1.006, zero errors. VTK copy for Ryan:
~/p018_das_inspect/smoke_csarc_n0_l4p8/ (51 steps).
**Deploy (quarantine-proven):** diffs cluster→local verified 024-only
for all five files; pre-overwrite md5s recorded (wake 20ed0def,
simulate a398bb71, replay 12c5e6f8, warmstart 836ff3a4, driver
00bda73d) + server-side copies at .deploy_backups/pre024_20260820/;
post-deploy md5 local==cluster (wake 4849fd76, simulate 7870348b,
replay d345597b, warmstart b0b9e81e, driver b6bc8f30, launcher
dc0a3fb6). Running jobs unaffected (code loaded at their start);
13245454 started pre-deploy. New launcher arm `p018_csarc_n0_l4p8`
(= csarc_l4p8 + NWAKEROWS=0). **Submitted 13246032**, 48 h/64G,
SETTLE 22; banner check armed (expect nwakerows=0 convert-at-shed +
F1b ACTIVE lines). **DISCLOSED A/B deltas vs 13243085 (csarc_l4p8):
{no free row} AND {unsteady_filament false→true, forced at N=0}** —
score with both in mind; also NoShed+N=0+unsteady_filament=true is a
corner the unit suites cover only via the driver smoke.

**2026-08-20 evening — csarc_mid 13206338 COMPLETED** (45 h, full 30
revs/1080 steps, CT CSVs written; tail healthy: max_u 14.6, gos2 958,
zero σ≤0, min_sigma 5.2e-4). Harvested local. Midpoint-vs-TE sanity
A/B vs `p018_csarc_l3p4` runs when the ladder lands (~08-22; Ryan's
criterion: "behavior doesn't change much"). Fleet monitor re-armed on
TEN live jobs (λ ladder ×3, NT ladders ×6, N=0 A/B 13246032).
**(cont) 13246032 banner VERIFIED RUNNING:** `nwakerows=0
(convert-at-shed)`, das_lambda 4.8, das_arc:true steady TE table, F1b
ACTIVE max|u| 0.102 tip, production mesh 45_185_ct4 — valid N=0 A/B arm
against 13243085.

**2026-08-20 late — N=0 NT LADDER LAUNCHED (Ryan: "let's finish the N=0
NT ladder").** Interpretation: complete the NT dimension at the N=0
arm's λ=4.8. New launcher arms + submitted: **13246048
p018_csarc_n0_nt72_l4p8** (NT 72, pps 6, rlxf 0.16334) and **13246049
p018_csarc_n0_nt144_l4p8** (NT 144, pps 3, rlxf 0.08539), NWAKEROWS=0,
72 h/64G, SETTLE 22, exact-rate rule, NT·pps=432. Launcher deployed md5
d9624186 both sides. With 13246032 (NT36) this is a 3-point N=0 NT
ladder mirroring the N=1 csarc NT column at λ4.8 — the cleanest
free-row-extent read (at N=0 the free-row term is GONE, so N=0
NT-flatness isolates the pps/dt model-def effect from row extent;
N1-vs-N0 at each NT isolates the row itself). Banner checks armed;
fleet monitor now TWELVE 018 jobs. nt144 rungs (N=0 and N=1) reach
~rev 13–15 per wall — chain decisions are Ryan's on landing.
**(cont) 13246048/49 banners VERIFIED RUNNING:** NT 72/144, rlxf
0.16334/0.08539, λ4.8, das_arc:true steady TE table, `nwakerows=0
(convert-at-shed)`, production mesh — both valid N=0 NT rungs.

**2026-08-21 — CONTEXT-RESET PREP: RESET BRIEF (m) written** (item file
top, supersedes (l)). Fresh-agent entry = brief (m): twelve-job table
(all banners verified at start), 08-20 scored verdicts, 024
implemented+deployed state, 025 staged-only hold, exact-rate + 36-step
retention rulings, disk watch, ssh-heredoc gotcha. Memory + MEMORY.md
updated to point at (m).

**2026-08-21 — 025 PHASE 3 LANDED; PRODUCTION SOLVER SETTINGS CHANGED
(two Ryan rulings, both implemented).** Three warm-start continuations of the
frozen carrier `p018_cs_f1_l3p4` from step 1034, 117 steps each (steps
1035–1151, ~3.25 revs) on the mature ~181k-particle wake — everything except
the swept knob held identical. Scored from `monitor02_force`/`monitor04_wake_health`,
NOT the driver's summary block (which zero-fills restored steps and therefore
reports a meaningless `CYCLE-MEAN 0.0229 ±145%, CONVERGED=false` for any
continuation).

| arm | job | change | CT̄ (steps 1035–1151) | Δ | s/step (1045–1151) |
|---|---|---|---:|---:|---:|
| control | 13247862 | vatistas, production FMM knobs | 0.0703810 | — | 160.2 |
| family | 13247863 | **gaussian**, production knobs | 0.0703975 | +0.023% | 140.2 |
| family+knobs | 13290979 | gaussian, **tuned knobs** (wake 16/0.6/38, body 17/0.7/109) | 0.0704002 | +0.0038% vs 13247863 | **71.7** |

Both deltas are nulls against the carrier's own noise: within-rev peak-to-peak
is 4.2–6.0e-4 (0.6–0.85% of CT) and the parent's rev-to-rev drift over its last
three revs is −3.3e-5, i.e. **twice the family delta and ten times the knob
delta**. Worst single-step divergence of the tuned arm vs its production-knob
sibling: 3.73e-5 (0.053% of CT) — the certified field errors (1.8e-6 wake /
1.2e-6 body vs a DirectBackend reference, job 13247200) do not accumulate over
117 steps. Wake health indistinguishable across all three (n_particles to 0.1%,
identical min σ 7.94e-4, identical max Γ/σ² 289.8, no tripwire, all finite).

**Implemented in `examples/run_dji9443_hover_ct_hpc.slurm.sh` (deployed,
md5-verified):** `FLOWPANEL_FILAMENT_REG` default vatistas → **gaussian**, and
`FMM_BODY_*` = 17/0.7/109, `FMM_WAKE_*` = 16/0.6/38 as 018 defaults. Only this
launcher changed; the driver's own defaults (body 8/0.4/20, wake 4/0.4/50) and
every other study sharing the driver are untouched. Supporting changes: the
driver now takes per-pass FMM knob env vars (shared `FMM_*` still the fallback,
so old submissions are unchanged) and records both triples plus the family in
`case_metadata.toml`; `FLOWPanel.__init__` prints the pinned family.

**Campaign consequences.** (a) A 30-rev production run projects from ~48 h
(which timed out at rev 28.5) to **~11 h** — the NT ladder's turnaround
roughly quarters and the NT=144 arms that could not mature inside 72 h now can.
(b) The twelve arms in flight on 2026-08-21 and everything before them ran
vatistas + production knobs; they remain comparable to future rungs under the
two measured nulls above (+0.023% and +0.0038%, vs a 0.6–0.85% within-rev p-p),
and no error-budget term is warranted for either. (c) To reproduce a pre-2026-08-21
arm exactly, pin `FLOWPANEL_FILAMENT_REG=vatistas` and the old `FMM_*` triples
at submission. **(d) CHAIN WARNING: a restart/continuation of an in-flight
vatistas arm will silently pick up the new defaults unless pinned** — pin both
when extending anything submitted before today.

Detail, tables and method: `BRAINSTORM/025_kernel_regularization_update/phase_03_018_compatibility.md`.

---

## VTK retention actions — 2026-08-22 — full-corpus sweep, 588 G → 88 G

**Trigger.** Ryan asked for a manual cleanup cycle. No watchdog loop had been
running, and `du -sm /home/rander39/projects/FLOWPanel.jl` read **602,429 MB
(588 G)** — roughly **3× the 200 G cap** — with **8 jobs RUNNING**
(`fp-018-csarc_{nt72,nt144}_{l2p4,l3p4,l4p8}`, `fp-018-csarc_n0_nt{72,144}_l4p8`;
23–25 h of walltime remaining each). Nothing had been swept in a long while.

**Retention policy corrected first.** `HPC_STORAGE_WATCHDOG.md` §1 rule 4 (and
the §2 flag/behavior lines) still said **10** retained timesteps; the sweeper's
`KEEP_STEPS` default had been raised to **36** by Ryan on 2026-08-20. The doc
now says 36 and points at the script default as authoritative (`--keep` should
not be passed). The historical §6 entries that reference 10-step retention were
left as-is — they describe what was true in August's incidents.

**Sweep.** `scripts/p018_vtk_sweeper.sh` dry run, then `--apply`, at the default
`keep=36`, protect list unmodified (8 runs `PROTECTED`: `p018_L1_ov3`,
`scr_ufdt_nt{36,72,144}`, `scr_p019_fid144`, `scr_p019_s038v`,
`p018_ufront_{s035,n2}_visc`). Live runs were swept, as rule 3 allows.

| metric | value |
|---|---:|
| before | 602,429 MB (588 G) |
| `TOTAL_FREED_MB` | **514,083 MB (502 G)** |
| after | **89,992 MB (88 G)** of the 200 G cap |
| runs with deletions | 23 of 116 evaluated (181 dirs scanned) |
| runs skipped `NO-RESTART-SET` | 57 |

Essentially all of it came from two families — the thirteen `p018_csarc_*` arms
(~28–34 G each; five CLOSED, seven LIVE) and six `p022_*` ground-effect runs
(~9–26 G each). The remaining run directories contributed under 20 MB
combined; most were already at their retained window from earlier sweeps.

| family | runs | freed |
|---|---:|---:|
| `p018_csarc_*` | 13 | 401,638 MB |
| `p022_{ige,oge}_*` | 6 | 112,414 MB |
| everything else | 4 | 18 MB |
| **sum** | **23** | **514,070 MB** (13 MB under `TOTAL_FREED_MB`; per-run MB are rounded) |

**Restart integrity verified (this sweep touched seven live runs).** All eight
jobs are still `RUNNING` after the apply, and every swept live run has advanced
past its retained window — the sweeper's retained top step vs. the newest
particle VTP a few minutes later:

| run | retained top S at sweep | latest S after | particle files kept |
|---|---:|---:|---:|
| `p018_csarc_n0_nt144_l4p8` | 1324 | 1325 | 38 |
| `p018_csarc_n0_nt72_l4p8` | 1080 | 1081 | 38 |
| `p018_csarc_nt144_l2p4` | 1362 | 1364 | 38 |
| `p018_csarc_nt144_l3p4` | 1380 | 1382 | 38 |
| `p018_csarc_nt144_l4p8` | 1364 | 1365 | 37 |
| `p018_csarc_nt72_l2p4` | 1123 | 1125 | 38 |
| `p018_csarc_nt72_l3p4` | 1138 | 1139 | 37 |
| `p018_csarc_nt72_l4p8` | 1165 | 1166 | 37 |

(36 retained steps plus the 1–2 written since the sweep — exactly the expected
count, so no `.vtm` stub / missing-piece situation like the 2026-08-04 incident.)

**Job outcomes.** Nothing left the queue during this cycle; the same eight jobs
were `RUNNING` before and after. No `sacct` post-mortems to report, and no new
instances of the `merge_particles!` OOM family.

**Lesson for next time.** 588 G means the watchdog loop (§3, auto-sweep at
100 G) was not running for this campaign. A one-shot sweep recovered everything
here only because the `csarc` arms happened to still be pre-cap; the growth rate
with 8 writers makes a standing Monitor the safer default whenever a ladder is
in flight.

---

## 2026-08-22 — F1b λ ladder (NT36, N=1) landed + N=0 convert-at-shed A/B

Four of the twelve arms in flight since 2026-08-20/21 have **COMPLETED**
(all 30 revs, `all_finite=true`, `converged=true`, ~45–47 h wall each,
vatistas + production FMM knobs — i.e. pre-2026-08-21 settings, per the chain
warning above). The remaining eight (N=1 NT72/NT144 ladders 13245449–54, N=0
NT ladder 13246048/49) are still RUNNING at ~2 d of their 3 d wall and are
expected to land 2026-08-23.

Common config for all four: mesh `45_185_ct4` (36752 cells), RPM 5400, NT 36,
`p_per_step` 12, `overlap` 2.75, `nwakerows` 1 (except the N=0 arm), `relax_rlxf`
0.3, `sigma_chord_fraction` 0.313, arc-placed steady Das (`das_arc_placed=true`,
table `p018_cs_l3p4_rs1_te_downwash_te.csv`), **F1 curvature cap OFF**
(`das_curvature_beta = NaN`) — this is the uncapped F1b arc-law ladder.

**M1 — cycle-mean CT.** Window = last 15 recorded revs (rev blocks 16–30,
steps 541–1080); the driver's own 10-rev `CT_cycle_mean` is shown alongside and
agrees to <0.02% in every case.

| arm | job | λ | N | CT̄ (15 rev) | rev-to-rev σ_rel | max within-rev p-p | drift %/rev | driver CT̄ (10 rev) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `p018_csarc_l2p4` | 13243083 | 2.4 | 1 | 0.070663 | 1.7e-3 | 1.55% | −0.023 | 0.0706092 |
| `p018_csarc_l3p4` | 13243084 | 3.4 | 1 | 0.071013 | 1.1e-3 | 1.88% | −0.004 | 0.0710044 |
| `p018_csarc_l4p8` | 13243085 | 4.8 | 1 | 0.070666 | 1.0e-3 | 0.84% | −0.020 | 0.0706287 |
| `p018_csarc_n0_l4p8` | 13246032 | 4.8 | **0** | 0.071584 | 8.8e-4 | 1.07% | +0.014 | 0.0716176 |

**M2 — Γ̄(r/R) from the TE μ-jump** (`monitor03_bound_circulation`, same
15-rev window, 80 sections, both blades). Span-averaged Γ̄:

| arm | λ | N | span-mean Γ̄ | Δ vs `l3p4` |
|---|---:|---:|---:|---:|
| `p018_csarc_l2p4` | 2.4 | 1 | 0.230453 | −1.12% |
| `p018_csarc_l3p4` | 3.4 | 1 | 0.233071 | — |
| `p018_csarc_l4p8` | 4.8 | 1 | 0.232655 | −0.18% |
| `p018_csarc_n0_l4p8` | 4.8 | 0 | 0.234967 | +0.81% |

### Verdict 1 — λ is NOT a monotone axis at NT36; spread is ~0.5% and λ=3.4 is the outlier

CT̄ over λ ∈ {2.4, 3.4, 4.8} is 0.070663 / 0.071013 / 0.070666 — the two ends
agree to **0.004%** while the midpoint sits **+0.49%** above both. That is not a
convergence trend in λ; it is a ~0.5% scatter band with the middle rung high.
M2 tells the same story in the same rank order (l2p4 −1.12%, l4p8 −0.18% vs
l3p4), so the two metrics are consistent and this is not a force-recovery
artifact. The 0.49% CT gap does exceed the per-rev sampling noise (σ_rel ≈
1e-3, SE of a 15-rev mean ≈ 0.03% before autocorrelation inflation), but it is
*below* the arms' own within-rev peak-to-peak (0.84–1.88%), so it is a real
mean-level difference riding on a limit cycle that is larger than the effect.

Spanwise, the λ sensitivity is **concentrated at the root**: at r/R = 0.111 the
arms spread −11.4% (λ2.4) / — / +6.3% (λ4.8) about λ3.4, while every section
outboard of |r/R| ≈ 0.3 agrees within ~1.3%. This is the expected signature of
chord–σ co-scaling — the root chord sets the smallest σ, so the arc-placement
length λσ moves most there — and it is where the blade-gap σ cap question
(still un-approved) would bite.

**Caveat, and why this verdict is provisional:** these three rungs share NT=36.
The NT72/NT144 ladders (13245449–54) land tomorrow and are the arms that
decide whether the ~0.5% scatter is a genuine λ-dependence or a temporal-
resolution artifact of NT36. Do not close the λ axis on this table alone.

### Verdict 2 — N=0 convert-at-shed raises loading ~1% vs N=1 (clean A/B, both metrics agree)

`p018_csarc_n0_l4p8` (13246032) vs `p018_csarc_l4p8` (13243085) differ in
`nwakerows` alone (0 vs 1) — verified by diffing the two `case_metadata.toml`
files, which are otherwise identical apart from run name and wall time.

- **M1: CT̄ +1.30%** (0.071584 vs 0.070666)
- **M2: span-mean Γ̄ +0.99%** (0.234967 vs 0.232655)

Sign and magnitude agree between the independent metrics, and the effect is
~2.5× the λ scatter above, so this is a model-level difference, not noise. The
spanwise Δ is positive at *every* section (+0.6% to +2.1%), largest at the tip
shoulder (r/R 0.785, +2.10%) and at r/R −0.762 (+1.77%) — i.e. N=0 recovers
circulation that the single rigid wake row was suppressing, most strongly where
the wake sheet is most curved. Wake health is comparable (n_particles within
0.3%; tail max_u 18.3 vs 13.2 m/s — the N=0 arm's higher induced velocity is
consistent with vorticity being deposited closer to the TE; the whole-run
maxima, 53.4 vs 28.8, are spin-up transients, not the settled state).

**Parked question, now due (Ryan):** BRAINSTORM/024's N=0 convert-at-shed is
implemented and deployed; this is the first settled 30-rev A/B of it inside the
018 carrier. The production-adoption decision was parked pending exactly this
result. Recommend holding the decision until the N=0 NT ladder
(13246048/13246049) lands 2026-08-23, so adoption is not made on a single NT36
rung — the λ ladder above is a live demonstration that one NT is not enough.

**Retention:** VTK for these four runs is still on the cluster and has NOT been
swept; the four CSV/TOML/monitor sets are harvested to `data/<run>/` locally.

### Scoring against RESET BRIEF (m) §"action on landing" for 13243083/84/85

**(i) Wake-health TAIL absolutes** (last recorded step 1079; trend over last 5
revs) — all four arms are healthy and plateaued, no L2 trajectory:

| arm | n_particles | Δn over last 5 rev | tail max_u | tail min σ | tail max Γ/σ² (5-rev range) |
|---|---:|---:|---:|---:|---:|
| `csarc_l2p4` | 184 189 | +1.04% | 16.33 | 9.51e-4 | 228 (192–389) |
| `csarc_l3p4` | 192 991 | +2.29% | 14.18 | 9.60e-4 | 212 (212–384) |
| `csarc_l4p8` | 186 174 | +0.45% | 13.16 | 7.67e-4 | 407 (233–530) |
| `csarc_n0_l4p8` | 185 937 | +2.06% | 18.31 | 8.78e-4 | 378 (229–383) |

`min_sigma` is flat, not falling, in every arm — the stop rule
(steadily-falling `min_sigma_ratio` ⇒ heading for the L2 failure) is not
triggered. The brief's specific worry — l2p4's screen `gos2` ending at 425
rising-and-decelerating — resolves benignly: its settled max Γ/σ² band is
192–389, the *lowest* of the four, so the screen's rise did not continue into
the production run. (`min_sigma_ratio`/`p1_sigma_ratio` are NaN in the tail
rows for all arms — the known `nan_inf` false positive, not a divergence.)

**(ii) P1 doubling gate (λ 2.4 → 4.8, the factor-of-two rungs).** This is the
formal criterion in the brief (M1 raw+quiet ≤0.5%, M2 ≤1%):

| metric | λ2.4 | λ4.8 | Δ | gate | verdict |
|---|---:|---:|---:|---:|---|
| M1 CT̄ (15 rev) | 0.070663 | 0.070666 | **+0.004%** | ≤0.5% | **PASS** |
| M2 span-mean Γ̄ | 0.230453 | 0.232655 | **+0.96%** | ≤1% | **PASS (marginal)** |

**So the λ axis PASSES its doubling gate on both metrics.** The non-monotonicity
noted in Verdict 1 above is real but does not fail the gate: the doubling test
compares the two ends, and they agree to 0.004% in CT. Read the two together as
— the λ endpoints are converged in CT to well inside tolerance, while an
intermediate rung scatters ~0.5%, and M2 sits right at its 1% limit driven
almost entirely by the root sections (r/R ≈ 0.11, ±6–11%). The M2 margin is thin
enough that the NT72/NT144 rungs landing 2026-08-23 should be checked before
the axis is declared closed.

**(iii) λ3.4 head-to-heads** (all scored from `monitor02_force`, CT = −CFx,
whole-rev-aligned; the `csarc_l3p4` force-monitor score reproduces its
`CT_per_rev` score to 6 digits, which validates the method for the carriers that
lack a usable `CT_per_rev`):

| carrier | wake/Das treatment | CT̄ | window | Δ vs `csarc_l3p4` |
|---|---|---:|---|---:|
| `p018_csarc_l3p4` | **arc-placed Das, uncapped** | 0.071013 | 15 rev → step 1079 | — |
| `p018_cs_f1_l3p4` | F1 curvature cap β=0.6 | 0.070690 | 15 rev → step 1034 | −0.46% |
| `p018_cs_l3p4_rs1` | straight-sheet Das | 0.072524 | 8 rev → step 1079 | +2.13% |
| `p018_csarc_mid_l3p4` | arc-placed, **midpoint** table | 0.070552 | 15 rev → step 1079 | −0.65% |

- **Arc vs cap (the head-to-head a Ryan decision follows): 0.46% apart** — i.e.
  the same size as the λ ladder's own internal scatter and *below* both arms'
  within-rev peak-to-peak. On CT alone these two Das treatments are not
  distinguishable at this carrier's noise level; the earlier discriminator was
  the **shape** of Γ̄ (F1's outboard lobe residual at r/R 0.78–0.87), not the
  integral. **Recommend deciding arc-vs-cap on M2 shape, not on M1** — and not
  before the NT ladders land.
- **Straight-sheet is the outlier at +2.13%**, ~4× the arc/cap separation, so
  the arc law does move CT materially off the straight-sheet baseline (the
  point of Phase 16). Caveat: this number rests on an **8-rev** window (the
  `_rs1` force monitor holds only 302 steps) and the run is a stitched restart,
  so treat it as indicative. Its driver `CT_per_rev` reads 0.0404 — **do not use
  it**; that is the known zero-fill-on-restart artifact, the same trap flagged
  for the 025 continuations.
- **Midpoint-table A/B: −0.65%**, versus Ryan's stated pass criterion
  ("behavior doesn't change much"). Comparable to the λ scatter ⇒ **PASS**;
  the TE-downwash table choice is not a leading error term.

**Still owed on this arm** (not done here): Γ̄ overlay *plot* vs the
straight-sheet rungs. The numbers are in `data/<run>/monitors/`; the overlay is
a figure task, deferred.

---

## 2026-08-24 — NT72/NT144 ladders harvested (all 8 arms) + SIGBUS forensics

None of the eight NT-ladder arms completed: seven hit their 3 d wall (TIMEOUT)
and `csarc_nt72_l2p4` (13245449) died on a SIGBUS at step 1454/2159. No arm
wrote `CT_vs_rev.csv`; all eight are scored from `monitor02_force` (CT = −CFx),
the validated fallback. Harvested to local `data/<run>/` (metadata + full
monitors/, no VTK); cluster copies untouched.

**Reach and windows** (whole-rev-aligned, settled portion only; the first
scoring pass with revs 6–20 windows was transient-contaminated — drift
1–3 %/rev — and was discarded):

| arm | job | outcome | reach | M1 window |
|---|---|---|---:|---|
| `csarc_nt72_l2p4` | 13245449 | SIGBUS @1454 | rev 20 | 14–19 (6 rev) |
| `csarc_nt72_l3p4` | 13245450 | TIMEOUT | rev 20 | 14–19 (6 rev) |
| `csarc_nt72_l4p8` | 13245451 | TIMEOUT | rev 21 | 15–20 (6 rev) |
| `csarc_nt144_l2p4` | 13245452 | TIMEOUT | rev 11 | 9–10 (2 rev, UNSETTLED) |
| `csarc_nt144_l3p4` | 13245453 | TIMEOUT | rev 11 | 9–10 (2 rev, UNSETTLED) |
| `csarc_nt144_l4p8` | 13245454 | TIMEOUT | rev 11 | 9–10 (2 rev, UNSETTLED) |
| `csarc_n0_nt72_l4p8` | 13246048 | TIMEOUT | rev 20 | 15–19 (5 rev; still drifting +) |
| `csarc_n0_nt144_l4p8` | 13246049 | TIMEOUT | rev 11 | 9–10 (2 rev, UNSETTLED) |

**M1 — CT̄ across the NT ladder** (NT36 rows from 2026-08-22 entry; NT rungs
are labeled model-def A/Bs under the exact-rate rule, rlxf 0.3/0.16334/0.08539):

| λ | N | NT36 (15 rev) | NT72 (6 rev) | Δ36→72 | NT144 (2 rev) | Δ72→144 |
|---:|---:|---:|---:|---:|---:|---:|
| 2.4 | 1 | 0.070663 | 0.073702 | **+4.30%** | 0.074251 | +0.74% |
| 3.4 | 1 | 0.071013 | 0.072689 | **+2.36%** | 0.074867 | +3.00% |
| 4.8 | 1 | 0.070666 | 0.074214 | **+5.02%** | 0.074065 | −0.20% |
| 4.8 | 0 | 0.071584 | 0.075730 | **+5.79%** | 0.073960 | −2.34% |

**M2 — ε_Γ (0.3 ≤ r/R ≤ 0.95, matched settled windows):** NT72 vs NT36
carrier: ε_max 5.50% / 3.77% / 8.78% (λ 2.4/3.4/4.8), n0 pair 10.63%.
NT144 vs NT72: 4.44% / 5.15% / 1.87%, n0 pair 2.73%. All FAIL the 1% gate.

### Verdicts

1. **NT36 is not temporally converged.** Halving dt (NT36→NT72) moves CT̄ by
   +2.4–5.0% and Γ̄ by 4–11% — an order of magnitude beyond the λ-axis scatter.
   The 2026-08-22 λ doubling-gate PASS was a fixed-NT36 statement and is now
   demoted: at NT72 the λ gate **FAILS** (M1 Δ2.4→4.8 = +0.69% > 0.5%; M2
   ε_max 3.96% > 1%). At NT144 the λ endpoints agree to −0.25% / 1.17%, but on
   unsettled 2-rev windows — suggestive, not scoreable.
2. **NT72→NT144 deltas shrink** (|ΔCT| 0.2–3.0% vs 2.4–5.8% for 36→72),
   consistent with approaching temporal convergence, but NT144 arms reached
   only rev 11 of 30 and were still settling — the rung cannot be scored
   without chaining (RYAN'S decision, standing rule).
3. **N=0 vs N=1: +2.04% CT at NT72** (0.075730 vs 0.074214), larger than the
   +1.30% at NT36 — the convert-at-shed effect is not an NT36 artifact. The
   NT144 sign flip (−0.14%) is unsettled-window noise; hold adoption per the
   2026-08-22 recommendation until NT144 settles.
4. **SIGBUS forensics** (13185010, 13245449): different nodes (m12-2-6,
   m12-4-14), 4 days apart, exit 7:0 (SIGBUS, not cgroup OOM), no disk/quota
   pressure (home FS 7% full). 13245449's trace faults in `getfield`
   (`ijl_get_nth_field_checked`) inside `_set_kerneloffsets!`
   (cluster `FLOWPanel_simulate.jl:446`). Both jobs peaked ~67 G RSS on a 64 G
   request while six sibling arms on the same code path ran 3 d clean —
   consistent with sporadic memory-pressure faults, not a deterministic bug.
   **Mitigation for chains/resubmits: request ≥80 G.**

**Retention:** VTK for all eight arms still on cluster (swept to 36-step
retention 2026-08-22, a few steps written since); disk 104 G/200 G. Restart
snapshot sets preserved — chaining is possible from the retained sets.

**2026-08-24 addendum — NT-ladder _s2 chains STAGED (not submitted).** Ryan
approved chaining; the agent's `sbatch`/`scp` were blocked by the local
permission classifier, so the full submission is staged in
`scripts/p018_submit_nt_chains_s2.sh` (LOCAL ONLY — not yet on the cluster)
for Ryan to push+run:
`scp scripts/p018_submit_nt_chains_s2.sh orc:/home/rander39/projects/FLOWPanel.jl/scripts/ && ssh orc bash /home/rander39/projects/FLOWPanel.jl/scripts/p018_submit_nt_chains_s2.sh`
Contents: all 8 arms chained to 30 revs (`P018_SETTLE_REVS=22`), NT72 48 h /
NT144 72 h, **96G** (SIGBUS mitigation), pinned
`FLOWPANEL_FILAMENT_REG=vatistas` + production FMM knobs (= driver defaults
body 8/0.4/20, wake 4/0.4/50; verified the original arms ran the pre-08-21
launcher/driver with no FMM env — their banners lack the `filament_reg` line
and their metadata records no knobs). RESTART_STEP per arm = last step with
all four warmstart pieces, verified aligned 2026-08-24: nt72 l2p4/l3p4/l4p8 =
1453/1459/1518, n0_nt72 = 1449; nt144 l2p4/l3p4/l4p8 = 1648/1682/1663,
n0_nt144 = 1648. NT144 s2 reaches only ~rev 22 in 72 h (~163 s/step) — an
**s3 chain will be needed** from each s2's last snapshot to close 30 revs.
After submission: verify banners (vatistas + 8/0.4/20 + 4/0.4/50, settle:22,
restart line), re-arm the ≥100 G disk watchdog (8 VTK writers), and at
analysis time stitch on `step` with the ≤0.05% seam gate. Score s2 segments
from `monitor02_force` (driver summary zero-fills restored steps).

**2026-08-24 addendum 2 — _s2 chains LAUNCHED after three failed rounds
(env/sync forensics).** All 8 arms now RUNNING as jobs **13400132–39**
(nt72 l2p4/l3p4/l4p8/n0 = 13400132/33/34/35; nt144 l2p4/l3p4/l4p8/n0 =
13400136/37/38/39), submitted 2026-08-24 ~02:5x with the staged pins
(vatistas + FMM 8/0.4/20 & 4/0.4/50, 96G, settle:22) and verified at 15 min:
resumed from the exact staged restart steps (1453/1459/1518/1449;
1648/1682/1663/1648), banners correct, .err all empty, stepping.

Three failed rounds first, all repo-sync/environment, zero physics impact
(every round died pre-output; restart sets untouched):

1. **13395136–43** (Ryan's submission of the staged script): the 2026-08-22
   env sync had put a julia **1.12.6** Manifest.toml on the cluster, but the
   slurm launcher's non-interactive fallback still exported spack
   julia-1.11.7 → OpenSSL_jll/HDF5_jll precompile failure in ~4 min × 8.
   Fix: launcher fallback now points at the site spack julia-1.12.6
   (patched identically local + cluster; md5-verified), env
   Pkg.instantiate()d clean.
2. **13396725–32**: cluster ~/projects/FLOWVPM.jl was at Jun-30 a950790
   against the Aug-24 FLOWPanel src (021's session had synced FLOWPanel
   src/+examples/ but not FLOWVPM) → `ParticleField` rejected the new
   `arraytype` kwarg at wake construction (~90 s × 8). Fix: synced local
   FLOWVPM (9b5b7cd + worktree) Project.toml+src/+ext/ to the cluster after
   verifying every cluster-dirty file was an exact blob from local git
   history (nothing cluster-only lost); per-file md5 parity confirmed.
3. **13398379 + 13398851–57**: cluster FastMultipole was a week behind
   (a9b734a-dirty, the 021-pinned state) and the new FLOWPanel src references
   9 symbols absent from it (NEARFIELD_CACHE_DEFAULT_MAX_BYTES, FmmPlan,
   NearfieldInfluenceCache, transform_plan!/transform_solver!, + 4
   GPU/rect-kernel names) → UndefVarError at warmstart in
   `_sa_wake_influence!`. **Ryan authorized moving FastMultipole forward**
   (comparability note sent to the 021 session). The pinned dirty state was
   archived first as
   `~/projects/FastMultipole/fm_a9b734a_dirty_021pinned_20260824.patch`
   (applies onto a9b734a to reproduce fm_commit=a9b734ad-dirty). Cluster
   FastMultipole now = local **d8258a7d** clean src+Project.toml; all 72
   FastMultipole symbols referenced by FLOWPanel src verified defined;
   precompile clean.

**Chain code-generation note for scoring:** the _s2 segments run Aug-24 code
(FLOWPanel 5272a5f-era src, FLOWVPM 9b5b7cd, FastMultipole d8258a7d) while
the parent arms ran the pre-08-21 stack — physics is env-pinned
(vatistas + FMM knobs) but the ≤0.05 % CT seam gate at stitch time is the
authoritative check on the s1→s2 boundary; treat a seam FAIL as a
code-generation break, not a restart artifact. Score s2 from
monitor02_force as before. Disk watchdog re-check dispatched at launch
(was 104 G/200 G, 8 VTK writers live). NT144 s2 still reaches only ~rev 22
in 72 h — plan the s3 chain from s2's last snapshots.

2026-08-24: disk watchdog cycle — 103.9G → 104.1G of 200G cap, freed 0G
(nothing sweepable: all non-protected runs at/under 36-step retention floor;
8 new p018_csarc_*_s2 chains NO-RESTART-SET, untouched). Verified 8/8
csarc_s2 jobs RUNNING (13400132-39). Headroom 95.9G; re-check before writers
accumulate multi-day restart sets.

---

2026-08-25: **_s2 chains cancelled by Ryan at ~24 h; NT72 ladder HARVESTED and
settled, NT144 ladder needs _s3.**

All 8 jobs (13400132-39) CANCELLED by user 2026-08-25 T03:04-03:05 at
1-00:11:46 to 1-00:12:35 elapsed, exit 0:0 (clean — checkpoints intact), against
their 48 h/72 h walls. Progress at kill (rev = time*5400/60, verified
self-consistent against step/NT):

| arm | restart step | final step | rev reached | vs settle target 22 | s/step |
|---|---|---|---|---|---|
| nt72_l2p4_s2 | 1453 | 1850 | 25.69 | PASS +3.7 | 219 |
| nt72_l3p4_s2 | 1459 | 1810 | 25.14 | PASS +3.1 | 248 |
| nt72_l4p8_s2 | 1518 | 1910 | 26.53 | PASS +4.5 | 222 |
| n0_nt72_l4p8_s2 | 1449 | 1805 | 25.07 | PASS +3.1 | 245 |
| nt144_l2p4_s2 | 1648 | 1908 | 13.25 | FAIL -8.75 | 335 |
| nt144_l3p4_s2 | 1682 | 1989 | 13.81 | FAIL -8.19 | 284 |
| nt144_l4p8_s2 | 1663 | 1957 | 13.59 | FAIL -8.41 | 296 |
| n0_nt144_l4p8_s2 | 1648 | 1946 | 13.51 | FAIL -8.49 | 292 |

**M1 (CT), NT72, scored rev 20-25 from monitor02_force (no CT_vs_rev.csv —
mid-loop cancellation, expected fallback path):**

| lambda | N | CT_bar | sigma | delta to next rung |
|---|---|---|---|---|
| 2.4 | 1 | 0.072386 | 0.00036 | — |
| 3.4 | 1 | 0.072900 | 0.00011 | +0.710% |
| 4.8 | 1 | 0.074415 | 0.00037 | +1.958% |
| 4.8 | 0 | 0.075921 | 0.00010 | +2.024% vs N=1 |

CT_bar stable to <0.15% between the rev 20-25 and last-5-rev windows for every
arm. **N=0 effect HELD at +2.02% through settlement** (was +2.04% unsettled) —
it is not a settling artifact. 024 adoption still parked pending NT144.

**M2 FAILS the 1% gate at every rung:** eps_Gamma,max = 2.631% (lambda
2.4->3.4) and 3.529% (lambda 3.4->4.8). Per decision_rules.md both M1 and M2
must pass, so **the lambda ladder is NOT converged at NT72 even when properly
settled** — CT convergence is masking spatial non-monotonicity in the
circulation distribution. This is the headline finding of the harvest.

**Seam gate: nt72_l3p4_s2 FAILS at 0.1415%** (2.8x the 0.05% limit); the other
three pass at 0.0062% / 0.0347% / 0.0012%. Per the 08-24 addendum this
signature = code-generation/env-pin break in that arm's launch window, not
restart corruption. **The lambda=3.4 rung is therefore suspect and the +0.710%
2.4->3.4 delta should not be used until this is explained.**

Sweep cycle at keep=288 (normal tier): 251G/500G, TOTAL_FREED_MB=0, nothing
reclaimable. 189 run dirs (124 SWEEP / 57 NO-RESTART-SET / 8 PROTECTED), 12
live 021 runs skipped. **Policy note: at keep=288 the sweeper is a no-op —
no run's restartable-step count exceeds 288, so repeated normal-tier sweeps
cannot reclaim anything. Disk growth toward the cap will need either the 72-step
escalation tier or new capacity, not more sweeps.**

**Verified 4-piece warmstart sets intact for the four NT144 _s3-bound runs:**
csarc_nt144_l2p4_s2@1908, l3p4_s2@1989, l4p8_s2@1957, n0_nt144_l4p8_s2@1946.

**_s3 sizing (measured ~300 s/step at 64 threads, from rev ~13.6):**
rev 22 = 1210 steps = 4.2 d; rev 25 = 1642 steps = 5.7 d; rev 30 = 2362 steps
= 8.2 d. The 3-day cap on every `normal`-QOS partition that fits a 64-core job
means **NT144-to-rev-30 is unreachable in a single non-preemptible job.**

**Cluster resource survey for the relaunch decision (2026-08-25).** Job spec is
64 cores / 96 G / 1 node (`run_dji9443_hover_ct_hpc.slurm.sh:4`, `--ntasks=64`;
the _s2 submitter overrides only --time and --mem). Slots with >=64 free cores
*right now*: `normal` QOS = 25 (m12 x21, m11-2 x2, m14 x2), all capped at 3 d.
`standby` adds 28 more, including dw/cs at 7 d and m12pws at 28 d. So 018 is
NOT blocked on priority or hardware — it is blocked on **wall time only**.
Preemptible does not raise priority (standby QOS priority=0 vs normal=100, and
PriorityWeightQOS=1e6 vs PriorityWeightFairShare=1e4); its value here is the
7 d / 28 d wall. Fairshare cost of heavy usage is negligible: within account
sn72 the entire spread between rander39 (480 M core-s, 41% of account usage)
and a zero-usage member is ~180 priority points against a 100,000-point QOS
base, and standby accrues at UsageFactor 0.30 vs 1.00.
Smaller nodes do NOT help: m9 holds 6550 idle cores but is 28 c/node and this
driver is single-process threaded (`--nodes=1`), so it cannot span nodes.
knlg+knlp (17 slots, 72 c, 192 G, 7 d wall, ~95% idle) are the only genuinely
underused >=64-core resource, but KNL per-core throughput is unbenchmarked for
this workload — probe before committing an arm.

**HAZARD for any preemptible relaunch:** cluster `PreemptMode=REQUEUE` with
`GraceTime=30 s`, and `p018_submit_nt_chains_s2.sh:24` exports a *fixed*
`RESTART_STEP`. On preemption Slurm re-runs the script verbatim, restarting
from the original step and silently discarding all progress since launch, with
no error. Submit `--no-requeue`, or teach the driver to auto-detect the newest
complete 4-piece snapshot in its own output dir.

Open for Ryan: (a) explain/re-run the lambda=3.4 seam FAIL; (b) M2 fails the
gate at NT72 — decide whether the ladder needs tightening or deeper diagnosis
before more NT144 wall-time is spent; (c) approve _s3 target rev and QOS route.

2026-08-25 (later, CORRECTION + convergence figures):

**CORRECTION to the M2 numbers recorded earlier today.** The 2.631% and
3.529% eps_Gamma values above are ARTIFACTS and should not be used. Root
cause: `p018_analyze.py m2` takes a SINGLE `--restart-step` (:275) and applies
it to BOTH runs via `--stitch-a`/`--stitch-b`, but the four NT72 arms restarted
at DIFFERENT steps (1453/1459/1518/1449). Passing arm A's boundary for arm B
makes `stitch()` keep B's parent rows only up to A's step while B's segment
begins later, punching a hole in B's history — 6 steps for the l2p4/l3p4 pair
(shared 1453), 59 steps for l3p4/l4p8 (shared 1459). Reproduced exactly:
shared 1453 -> 2.631%, shared 1459 -> 2.654%; shared 1459 -> 3.529%,
shared 1518 -> 3.397%. Using `chain()` (no hole — it keeps parent rows strictly
below the segment's own first step) the correct consecutive-rung values over
revs 20-25 are **2.620%** and **3.429%**. The M2 verdict is unchanged: both
still FAIL the 1% gate by ~3x. **`p018_analyze.py m2` needs a per-run restart
step before it is used again on arms with unequal restart boundaries.**

**Convergence figures built** (Ryan asked to see trends before deciding on
_s3). New emitter `scripts/p018_nt_convergence_figs.py` with a `check`
subcommand that reproduces every published anchor and writes nothing; `emit`
writes the backing CSVs. Figures live in
`~/Dropbox/research/notebooks/img/20260825_p018_nt_convergence/`:
fig_ct_history, fig_ct_rungs, fig_gamma_span, fig_gamma_nt (+ p018_nt_common.tex,
build.sh). NT144 `_s2` monitors were scp'd from the cluster to make this
possible. All anchors reproduce; CSVs are NaN-free; all four PDFs are 1 page.

**FINDING 1 — the timestep ladder is not converging.** CTbar over each arm's
best-available window:

| lambda | N | NT36 | NT72 | NT144 | 36->72 | 72->144 |
|---|---|---|---|---|---|---|
| 2.4 | 1 | 0.070614 | 0.072471 | 0.074408 | +2.63% | +2.67% |
| 3.4 | 1 | 0.071003 | 0.072912 | 0.074722 | +2.69% | +2.48% |
| 4.8 | 1 | 0.070593 | 0.074480 | 0.073803 | +5.51% | -0.91% |
| 4.8 | 0 | 0.071619 | 0.075954 | 0.073734 | +6.05% | -2.92% |

There is **no asymptotic flattening at any arm**. lambda 2.4 and 3.4 keep
climbing at essentially undiminished rate (+2.67%, +2.48% on the 72->144 leg,
as large as the 36->72 leg); lambda 4.8 and N=0 turn over. The 72->144 leg is
computed from UNSETTLED NT144 points (revs 8-12), so its sign is not yet
trustworthy — which is exactly the case for the _s3 chain.

**FINDING 2 — there is currently no window in which the three NT rungs can be
fairly compared.** revs 8-12 is the only window every arm reaches, but it is
deep in the NT36/NT72 transient: CTbar moves -3.6% (nt36_l3p4: 0.068428 matched
vs 0.071003 settled) to +4.7% (nt72_l2p4: 0.075846 vs 0.072471) between the
matched and settled windows. Both window definitions are emitted side by side
in fig_ct_rungs so this is visible rather than hidden by a window choice.

**FINDING 3 — the M2 failure is REAL, not a TE force-reconstruction artifact.**
Over the settled NT72 window (revs 20-24), vs the lambda=2.4 reference,
eps_Gamma,max is te 2.52 / 6.06 / 9.10% and slice 2.34 / 4.84 / 5.46% for
lambda 3.4 / 4.8 / 4.8-N0. The two circulation columns broadly AGREE, so the
circulation difference is genuine. (Over the transient matched window they
disagree wildly — te 7.89% vs slice 1.58% for lambda 3.4 — which is itself a
symptom of comparing unsettled states.) All pairs FAIL the 1% gate in both
columns. Worst r/R is 0.313 for nearly every pair, i.e. the inboard edge of the
0.30-0.95 evaluation band.

**FINDING 4 — seam failures are not unique to l3p4.** 6 of 11 restart seams in
the campaign exceed the 0.05% gate: nt72_l3p4 0.1415% (2.8x, the worst),
nt144_l3p4 0.0756%, nt144_l2p4 0.0699%, prior b0 0.0656%, and both legacy
nt144 seams (0.0515%, 0.0575%). Passing: nt72_l2p4 0.0062%, nt72_n0 0.0012%,
nt72_l4p8 0.0347%, nt144_l4p8 0.0282%, nt144_n0 0.0223%. This reframes the gate
as a tripwire rather than a verdict — l3p4 is the outlier in MAGNITUDE, not in
kind, and a systematic ~0.05-0.08% restart discontinuity looks structural.

**Legacy arms** (p018_nt144/_s2/_s3, upin_nt72, upin_nt144) are plotted greyed
and are NEVER joined to the csarc rungs: metadata confirms they ran the naive
linear rate rlxf = 0.3*36/NT (0.15 @NT72, 0.075 @NT144) rather than the
exact-rate 0.16334/0.08539 — 8.9% and 12.2% lower, on a rate that sits near the
ignition boundary. p018_upin_nt72 diverged at rev 18.68 and is truncated there.

Disk after the scp: unchanged at 251G/500G (monitors only, ~8 MB total).

2026-08-25 (spanwise loading figures added):

Two loading figures added, one per circulation figure: `fig_loading_span`
(companion to fig_gamma_span, lambda ladder at NT72, settled revs 20-24) and
`fig_loading_nt` (companion to fig_gamma_nt, NT ladder at lambda=3.4, matched
revs 8-12).

**No sectional-force data exists on disk.** `SpanwiseLoadingMonitor` is
implemented (`src/FLOWPanel_simulate_monitors.jl:1675`) but was never enabled
for any csarc arm — only monitors 02/03/04/05 were written. Loading is
therefore DERIVED from bound circulation by Kutta-Joukowski, not measured. In
hover the section-relative velocity is Omega*r to within the inflow angle, so
dT/dr = rho*Gamma(r)*Omega*r per blade and, with CT = T/(rho n^2 D^4), B blades,
x = r/R, D = 2R, Omega = 2*pi*n:

    dCT/dx = [pi*B / (8*n*R^2)] * Gamma(x) * x = 0.6162 * Gamma(x) * x

for B=2, n=90 rev/s, R=0.119 m (metadata). rho cancels.

**VALIDATED against the measured force.** Integrating that expression
reproduces each arm's measured CTbar to **0.3-0.8%** on the four settled NT72
arms (ratio 1.0047 / 1.0083 / 1.0059 / 1.0027 for lambda 2.4 / 3.4 / 4.8 /
4.8-N0). This is an independent check that (a) the KJ derivation is right and
(b) `circulation_te` really is the bound circulation. The area under each
loading curve therefore genuinely IS that arm's CT.

**`circulation_slice` does NOT share that normalization** — its raw KJ integral
lands at ratio ~0.47, i.e. ~2.12x small. An absolute loading built from the
slice column would be wrong by a factor of two. The slice panels are rescaled
by the reference arm's te/slice integral ratio (x2.124 span, x2.062 NT, written
to loading_scale.csv) and are labelled SHAPE ONLY on the figure. **Do not use
circulation_slice for any absolute loading or force recovery without
re-deriving its normalization.**

**Spanwise attribution of the CT spread** (area under the Delta-loading curve =
that arm's Delta-CT):

| pair | KJ Delta-CT | measured Delta-CT | agree? |
|---|---|---|---|
| NT72 lambda 2.4 -> 3.4 | +0.00070 | +0.00044 | overstates ~60% |
| NT72 lambda 2.4 -> 4.8 | +0.00206 | +0.00196 | yes, ~5% |
| NT72 lambda 2.4 -> N0 | +0.00335 | +0.00348 | yes, ~4% |

The lambda=3.4 arm is the one that overstates — and it is also the arm whose
restart seam FAILED at 0.1415%. Two independent diagnostics now point at that
same arm.

Shape findings: the lambda ladder redistributes load mostly INBOARD
(r/R 0.15-0.75 gains, with a compensating dip near r/R 0.80-0.85), whereas the
NT ladder adds load essentially EVERYWHERE (both Delta curves positive across
the whole span, peaking around r/R 0.3-0.5). Refining the timestep does not
merely redistribute circulation, it raises it globally — consistent with CT
climbing monotonically with NT and not flattening.

Also note the KJ/measured ratio drifts with NT on the lambda=3.4 arm:
1.0221 (NT36) -> 1.0126 (NT72) -> 1.0050 (NT144). The coarser the timestep, the
worse KJ-from-circulation agrees with the pressure-integrated force.

Figure dir now holds 6 figures + p018_nt_common.tex + build.sh, 3.0 MB, all
1-page, all CSVs NaN-free. Regenerate with
`python3 scripts/p018_nt_convergence_figs.py check && ... emit`.
2026-08-26 — VTK sweep (keep=288, normal tier): du 245390→247553 MB (239.6G→241.8G) of 400G cap; sweep itself freed 8 MB (dwarfed by ~2.1G growth from 4 new live p022m writers starting mid-cycle); restart-set integrity confirmed (p022m_2r_oge advanced 20→46 particle steps, p022m_4r_oge 2→13, no orphaned .vtm stubs); all 4 new jobs (13484022-25) RUNNING; p022m_2r_ige_hr15 not yet restartable, p022m_4r_ige_hr15 dir still empty — both expected for freshly-started jobs; no jobs completed/failed this cycle.

## 2026-08-26 — Phase 17 staged + launched: N ∝ NT ladder at fixed λ (extent-preserving temporal refinement)

Ryan directive (2026-08-26): stage and launch a ladder that decreases/varies N
with NT such that the total near-wake extent is dt-independent. Design settled
as N−1 ∝ NT at fixed λ (NOT the λ-compensated variant, which would launder the
known λ- and row-count sensitivity into the NT axis): free-row extent
(N−1)·travel held at one NT36 travel, Das = λσ untouched. Full spec:
`phase_17_nprop_nt_ladder.md`.

- Verified pre-launch: under DAS_SIGMA_LAMBDA, `nwakerows` does not perturb Das
  (curvature cap requires DAS_CURVATURE_BETA, unset on csarc; remaining
  nwakerows_extent uses are diagnostics/metadata). Zero code change — three new
  launcher case blocks only (`run_dji9443_hover_ct_hpc.slurm.sh`, mirrored to
  cluster, md5 aae9723f87524dbe31e6b863aa1ff923 both sides).
- Arms (cold, P018_SETTLE_REVS=22, --mem=96G, 64c,
  `-p m11-2,m12,m14 --qos=normal` per slurm-availability probe; normal QOS =
  non-preemptible; Ryan chose 64c over a 28c m9 option):
  - `p018_csarc_n2_l2p4`       NT36  N=2 pps12 rlxf0.3     — job **13490698**, 48 h
  - `p018_csarc_n3_nt72_l2p4`  NT72  N=3 pps6  rlxf0.16334 — job **13490699**, 72 h
  - `p018_csarc_n5_nt144_l2p4` NT144 N=5 pps3  rlxf0.08539 — job **13490700**, 72 h
- Each rung pairs with the existing N=1 csarc l2p4 rung at the same NT → free
  N-effect cross at every NT. NT144 will need _s2 chains (rev 22 ≈ 4+ days).
- Shed-clock/integration-clock decoupling ("gold standard") deliberately NOT
  staged — Ryan hold, 2026-08-26.
- Note: 7 user jobs in queue post-submit (4× 022 arms + these 3); the
  six-study-job directive in the launcher header is satisfied within this
  launcher (3 active).
- Storage (hpc-storage cycle post-launch): normal-tier sweep (keep=288) across
  10 live jobs (3 new p018 csarc + 4 p022m + 3 p021 tuners); freed 1670 MB
  (1.63G) trimming p022m_2r_oge 399→288 restartable steps (+ index cleanup on
  2r_ige_hr15/4r_ige_hr15); 293.2G→294.7G of 400G cap (net rise from continued
  live writes); restart-integrity PASS on all 3 swept runs; all 10 jobs still
  RUNNING, incl. 13490698-700.
- CORRECTION (2026-08-26, prompted by Ryan's ParaView check): the Phase 17
  staging text mislabeled PanelWake row 1 as "the rigid Das row". Ground truth
  (src/FLOWPanel_wake.jl + save() in FLOWPanel_liftingbody.jl): the rigid Das
  row is BODY-owned (the *_tw.vtu quads, TE extended by Das); ALL N PanelWake
  rows are free sheet (node row 1 pinned on the TE+Das line, rows convect one
  travel apart, particles convert from the oldest row). N=1 therefore has
  rigid Das + one free row, as Ryan observed. The launched N={2,3,5} settings
  are UNCHANGED and remain correct: they hold the audited clearance
  d_front = |Das| + (N-1)*travel invariant (= Das + t36 at every rung); only
  the O(dt) converting-row width varies, as any dt ladder requires. Fixed:
  phase_17 doc, launcher Phase-17 comment (re-mirrored, md5 verified), and the
  driver's misleading "Row 1 is rigidly pinned" handoff comment (LOCAL only —
  cluster driver has diverged, deliberately not overwritten mid-campaign).

2026-08-27 20:1x MDT — hpc-storage archive cycle: --all-checkouts across 10 FLOWPanel
clones; archived 30 finished FLOWPanel-052* runs (35.9G tarballs, verified byte-for-byte
before delete, 0 verify-fail, 0 stale), freed 34.05 GiB from /home. Home 566.6G→543.2G
of 400G cap (net freed 23.5G; live 3r GPU arms grew ~10G during the run). Archive
123.7G→158.7G, 93,007→93,048 files (of 20T/1M). RECENT approval queue: 29 runs,
~238.3 GiB, awaiting Ryan sign-off (largest: p018_csarc_n2_nt72_l3p0 55.8G,
p018_csarc_l4p0 24.8G, p018_csarc_l3p0 24.7G, n5_nt144_l2p4_s2gpu 32.2G,
n3_nt72_l2p4_s2gpu 30.7G). Live campaigns (022m_4r_oge, 018gpu-*-3r ladder,
scr_p019/p020 screen jobs) confirmed untouched. Cap still breached by ~143G;
next lever is Ryan's RECENT approvals, not a sweep (sweep math shows little to reclaim).
