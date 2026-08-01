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
