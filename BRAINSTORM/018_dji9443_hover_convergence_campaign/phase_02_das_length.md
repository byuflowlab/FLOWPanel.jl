# Phase 2 — Physical Das Length

**Objective:** find the smallest physical first-wake-row offset that still
converges. Per 014's rule and the campaign reconciliation: "converges" ≡ on
the log plateau (successive-doubling |ΔCT̄| ≤ 0.5%, consistent with the
~0.2%/doubling wing law); "smallest" = the plateau's lower edge. Below-plateau
points are the amplified regime (the +36.8% cliff) and are non-converged by
definition. Expected Das\* = 0.41c (η=1.0 frozen).

## Cases (carrier: NT=36, σ/c=0.68, N=4, velocity; B0 supplies η=1.0)

| tag | η | Das/c(0.75R) | floor clamp (span frac) | time |
| --- | --- | --- | --- | --- |
| `p018_das0p5` | 0.5 | 0.205 | ~20% (report measured) | 24 h |
| `p018_b0` (done, P1) | 1.0 | 0.41 | ~4% | — |
| `p018_das2p0` | 2.0 | 0.82 | ~0% | 24 h |
| `p018_das4p0` | 4.0 | 1.64 | 0% | 24 h |
| `p018_das0p5_floorhalf` (optional) | 0.5, DAS_MIN_DISPLACEMENT_R=0.005 | 0.205 | halved | 24 h |

Run the floorhalf case only if η=0.5 is within 1% of η=1.0 (i.e. the smallest
rung is a plateau candidate and floor contamination must be bounded).
Extend any case by restart if its mean hasn't settled per M1.

## Mandatory side job — FMM |Das|-radius discriminator (017 T4)

Every prior rotor η point used FMM, whose panel radius is inflated by |Das|
(`src/FLOWPanel_liftingbody.jl:996`). Before publishing any Das claim: run
`DirectBackend` vs `FastMultipoleBackend` on `RHPC_MESH=40_40`, NT=6, ~6 revs,
at η ∈ {0.5, 4.0}, locally (≤6 threads) or one 1-h cluster job. Pass = the
Direct-vs-FMM CT delta is η-independent (≤0.1% variation). Fail ⇒ the η ladder
carries an FMM artifact — rerun the two extreme ladder points on
DirectBackend-compatible settings before closing the phase.

## Decision

Das\* = smallest ladder value with |CT̄(Das) − CT̄(2·Das)| ≤ 0.5% and
ε_Γ ≤ 1%. Report the plateau slope (%/doubling) → error-budget term 2. Legacy
η ladder (0.06148/0.06942/0.07133/0.07190) is corroborating shape only.
If NO plateau exists at this σ, do not iterate here — proceed with η=1.0
provisionally and let Phase 4's σ refinement re-test the ladder's top two
points at σ=L1 (Das sensitivity may shrink with σ; record as a finding).

## Exit criteria

Das\* selected with plateau evidence + discriminator pass; ladder + Γ overlays
in the ledger; floor-clamp fractions recorded per case.

## Log

- 2026-07-31 — **T4 FMM discriminator PASSED.** 40_40, NT=6, velocity,
  SFS_OFF=true both arms, matched physical Das via η ∈ {0.083, 0.667}
  (= Das/c 0.20 / 1.64 at NT=6). Tail-mean CT: Das=0.20c fmm 0.056120 vs
  direct 0.056122 (−0.0033%); Das=1.64c fmm 0.056807 vs direct 0.056809
  (−0.0019%). Backend delta η-variation 0.0014% ≪ 0.1% ⇒ the FMM
  |Das|-panel-radius coupling does NOT contaminate the Das axis. CSVs kept at
  `data/p018_t4_discriminator/`. Two execution caveats discovered: (1) η must
  be scaled by 36/NT to hold physical Das in coarse-NT smokes (first attempt
  at η=4, NT=6 ⇒ Das≈9.8c diverged to CF≈430 — physically expected, discarded);
  (2) the SFS *direct* code path has a small-N bug
  (`FLOWVPM_subfilterscale_models.jl:54` `Estr_direct_multithreaded` builds a
  zero-step thread range when np < nthreads) — discriminator ran SFS-off in
  BOTH arms; production (FMM path) unaffected. Driver gained
  `RHPC_BACKEND=fmm|direct` (default fmm, no behavior change), synced +
  md5-verified to cluster.
- 2026-07-31 — Ladder jobs submitted: 12993802 das0p5, 12993803 das2p0,
  12993804 das4p0 (η=1.0 supplied by p018_b0, 12993801).
- 2026-07-31 — **First rung harvested: `p018_das0p5` (12993802) COMPLETED in
  12.5 h.** CT̄(revs 10–19) = **0.07006**, per-rev std ±0.00031, half-window
  drift −0.33% (oscillatory, no monotone trend). Metadata-verified SigmaOverlap
  overlap=2.0, nwakerows=4. Startup transient fully dead by rev ~9; per-rev
  peak-to-peak collapses to ~0.0003 in the settled window — the campaign
  config's scatter is ~5× tighter than every legacy overlap_pps run. Window is
  10 settled revs (run is 20 total): short of M1's ≥15 — extend by restart if
  the plateau-edge decision turns out to hinge on this rung. Ladder analysis
  deferred until b0/das2p0/das4p0 land. CSVs + circulation monitor at
  `data/p018_das0p5/`.
- 2026-07-31 — **`p018_das2p0` (12993803) COMPLETED in 14.9 h.** CT̄(revs
  10–19) = **0.07230**, per-rev std ±0.00027, half-window drift −0.16% (best
  settled of the three so far). η=2.0 metadata-verified. Running ladder:
  Das 0.205c → 0.07006, 0.41c → 0.07187 (+2.58%), 0.82c → 0.07230 (+0.60%).
  The 0.41→0.82c doubling delta sits just above the 0.5% convergence
  threshold and ~3× the 014 wing log-law slope (~0.2%/doubling) — plateau
  judgment (and whether Das\* = 0.41c or 0.82c) waits on das4p0's 0.82→1.64c
  delta, the M2 ε_Γ comparison, and b0's extended (s2) window since b0's mean
  was still drifting −0.38%/half-window. CSVs + circulation monitor at
  `data/p018_das2p0/`.
- 2026-07-31 — **`p018_das4p0` (12993804) COMPLETED (15.5 h) — LADDER ANALYSIS;
  PHASE 2 RESOLVED VIA THE PRE-REGISTERED NO-PLATEAU BRANCH.**

  M1 ladder (CT̄ revs 10–19, per-rev std):

  | Das/c | case | CT̄ | std | doubling Δ |
  | --- | --- | --- | --- | --- |
  | 0.205 | das0p5 | 0.07006 | ±0.00031 | — |
  | 0.41 | b0 | 0.07187 | ±0.00034 | +2.58% |
  | 0.82 | das2p0 | 0.07230 | ±0.00027 | +0.60% |
  | 1.64 | das4p0 | 0.07084 | ±0.00044 | **−2.02%** |

  M2 (ε_Γ from `circulation_te`, blades averaged, same window, 0.3≤r/R≤0.95,
  normalized by max|Γ̄|): das0p5→b0 max 8.97% / RMS 3.25%; b0→das2p0 3.11% /
  1.52%; das2p0→das4p0 5.32% / 2.78%. **No successive pair passes M1
  (|Δ|≤0.5%) AND M2 (ε_Γ,max≤1%).** Top rung is non-monotone (−2.02%),
  unlike the monotone-saturating legacy ladder; das4p0's case metadata is
  otherwise identical to das2p0's and its window is well settled (drift
  −0.14%), so the drop is a real configuration effect. The peak-Γ̄ station
  also migrates inboard with Das (r/R 0.740 → 0.695 → 0.583 → 0.583).

  **Decision (per the pre-registered branch "if NO plateau exists at this σ"):
  Das\* = 0.41c (η=1.0) PROVISIONAL; do not iterate here; Phase 4's σ=L1
  re-tests the ladder's top two points.** Plateau slope / error-budget term 2
  deferred to that re-test (no defensible slope exists at σ/c=0.68).

  Findings & caveats recorded:
  1. **Handoff-distance confound:** at fixed N=4 the sheet extent d≈N·Das spans
     1.2σ → 9.6σ (in σ=0.68c units) across the ladder — the two lowest rungs
     violate the d≳4σ criterion ruling 4 assumed. The Das axis at fixed N is
     inherently coupled to handoff adequacy ⇒ **Phase 5's N=8 spot-check is
     now decision-critical**, and should run at Das\*=0.41c where d only = 2.4σ
     (N=8 ⇒ 4.9σ).
  2. Floor clamps (measured from warm-start metadata |Das| vs 0.01R):
     das0p5 17.1% of stations, b0 2.4%, top rungs ≲2% by monotonicity —
     matches the phase table's predictions; das0p5's rung carries real floor
     contamination (but it is below-plateau anyway).
  3. `circulation_slice` is all-zero in every campaign monitor CSV — the
     slice cross-check estimator never ran (needs a driver flag). M2 rests on
     `circulation_te` alone; enable the slice estimator before Phase 9 quotes
     final Γ̄(r/R).
  4. b0's own window drift (−0.38%) is the largest of the four rungs; its s2
     extension (13007210, running) refines the η=1.0 number before Phase 3
     comparisons.

  T4 discriminator already PASSED (above). Phase 2 exit criteria met under the
  no-plateau branch. **Submitted on resolution (2026-07-31): 13011373
  `p018_nt72` (Phase 3, η=2.0=2η\*, 48 h), 13011374 `p018_L1` (Phase 4 cold,
  36 h — single-segment deviation from the 2×24 h plan, justified by the
  recalibrated ~30 h estimate), 13011375 `p018_L1_warm` (warm pilot from b0
  step 719, SETTLE=30 ⇒ ~8 flush revs + window, 36 h).**
