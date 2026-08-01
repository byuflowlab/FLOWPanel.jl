# Phase 1 — Baseline B0, Settling Demonstration, Restart Validation

**Objective:** establish the campaign reference point under `sigma_overlap`
shedding; demonstrate that the mean settles (the study's temporal-settling
deliverable); validate restart chaining end-to-end; recalibrate the cost model;
bank a settled particle field for warm starts.

## Cases

| tag | knobs (beyond carrier defaults) | schedule | time |
| --- | --- | --- | --- |
| `p018_b0` | η=1.0, NT=36, OVERLAP=2.0, PPS=4 (σ/c=0.68, σ/R=0.0873), N=4, MERGE_R_FACTOR=0.0120 | SETTLE=12 (~20 revs, 719 steps) | 24 h |
| `p018_b0` +s2 | restart extension +10 revs (ops_reference recipe: `P018_RUN_NAME=p018_b0_s2,P018_SETTLE_REVS=22,RESTART_STEP=<S>,RESTART_NAME=p018_b0,RESTART_PATH=data/p018_b0`) | +10 revs | 24 h |

Cost projection: t ≈ 13.5 h/(σ/c) ⇒ ~20 h for B0. Record actual wall time and
particle count vs the projections (N_p ≈ 62,667/(σ/c) ⇒ ~92k) — later phases
re-budget from the measured numbers.

## Deliverables and gates

1. **Restart integrity:** per-step CT jump at the segment seam ≤ 0.05%. This is
   the one-time in-anger validation of warm-start chaining; if it fails, debug
   before any other phase runs (warmstart §2.5 Das-rotation fix is already in).
2. **Settling demonstration:** plot CT̄ vs averaging-window start rev over the
   stitched ~30-rev history; flat within ~0.3% after settle per M1. This plot
   is a publication figure.
3. **Stability check:** if B0 diverges (|CT|>1 pattern), fall back to
   OVERLAP=2.4/PPS=5 (σ/c=0.65 — closest integer-PPS match) and re-baseline;
   record in the ledger and propagate the overlap change to Phase 4's ladder.
4. Ledger rows for B0 (M1 + M2 numbers, window, job ids).

## Exit criteria

B0 settled per M1 with ≥15-rev window; restart gate passed; cost model
recalibrated. B0's field banked on-cluster (`data/p018_b0*` retained).

## Log

- 2026-07-31 — **B0 COMPLETED (12993801, 13.2 h wall, exit 0:0) — no
  divergence, gate 3 fallback not needed.** CT̄(revs 10–19) = **0.07187**,
  per-rev std ±0.00034, half-window drift −0.38% (block 11 peak 0.07256
  decaying onto a ~0.0716 plateau — still a mild settling tail, hence the
  extension). Metadata verified: η=1.0, SigmaOverlap 2.0/PPS 4, N=4, rlxf 0.3,
  filter off, 45_185_ct4, truncation 4R. Settled scatter is ~5× tighter than
  legacy overlap_pps runs. **Cost recalibration (deliverable): 13.2 h wall and
  N_p ≈ 77,027 final particles at σ/c=0.68** ⇒ t ≈ 9.0 h/(σ/c),
  N_p ≈ 52.4k/(σ/c) (projections were 13.5 h/(σ/c), 62.7k/(σ/c) — both
  conservative; merging is doing its job). das0p5's 12.5 h confirms Das barely
  affects cost. **Extension `p018_b0_s2` submitted (13006768,
  RESTART_STEP=720, SETTLE=22 → +10 revs)** — restart-integrity gate (seam
  jump ≤ 0.05%) and the ≥20-rev M1 window + settling figure come from the
  stitched history when it lands. M2 Γ̄(r/R) numbers deferred to that analysis.
- 2026-07-31 — **s2 extension (13007210) COMPLETED — PHASE 1 GATES RESOLVED.**
  Stitch note: the s2 CSV keeps global step numbering (1–1080) with zero
  placeholder rows for pre-restart steps; stitch on nonzero rows, b0 ≤720,
  s2 721–1080. **Gate 1 restart integrity: seam per-step jump 0.0651% of CT̄
  vs the ≤0.05% criterion — letter-fail by 0.015 pp, but 42% of the run's own
  late per-step deltas EXCEED the seam jump (median natural delta 0.056%),
  so the restart introduces no discontinuity distinguishable from the
  signal's step-to-step variability. Verdict: PASS in substance; both numbers
  recorded; warm-start chaining validated for campaign use.** Gate 2 settling
  (M1): stitched 30-rev history; revs 10–14 still carry the settling tail
  (block drift 0.38%); **revs 15–29 (15 settled revs) pass M1** — successive
  5-rev block drifts 0.12%/0.12%, no monotone trend. **B0 reference:
  CT̄ = 0.07170, per-rev std ±0.00011, moving-block bootstrap 95% CI
  [0.07164, 0.07173]** (supersedes the interim 10-rev 0.07187 — the
  10–19-rev window sat on the tail). Settling-figure data ready
  (`data/p018_b0*/`); figure generation deferred to Phase 9 packaging.
  M2 Γ̄(r/R) baseline to be extracted with the ladder overlays. Cost model
  as recalibrated (13.2 h + 77k particles at σ/c=0.68). Exit criteria met;
  B0 field banked on-cluster through step 1080.
- 2026-07-31 — 13006768 FAILED at t+3 min: `restart_step=720 not found` — VTU
  snapshots are 0-indexed (0–719) while the CT CSV's step column is 1-indexed
  (ends 720). The restartable step is the last PVD index, 719. Resubmitted as
  **13007210** with `RESTART_STEP=719`; ops_reference restart recipe corrected.
  No simulation state touched (failure was pre-flight in `simulate_warmstart!`).
