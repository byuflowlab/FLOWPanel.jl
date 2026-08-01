# Phase 3 — Timestep at Fixed Physical Das AND σ

**Objective:** temporal convergence with both known confounds pinned. The
naive NT=72 run (0.07337, +2.9%) halved Das, σ, and rlxf at once — it is the
*confounded control* (a publication figure on why pinning matters), not a
ladder point.

## Pinning recipe (the crux)

At fixed OVERLAP, σ ∝ OVERLAP/(NT·PPS): **halve `P_PER_STEP` as NT doubles**
(the `SigmaOverlap` policy computes the true per-step count itself). Scale η ∝
NT to fix Das. Halve `RELAX_RLXF` per NT doubling (preserves the physical
relaxation rate — same convention as 2e). Never `DAS_REFRESH`.

## Cases (carrier: Das\*, σ from the Phase-2 carrier σ/c=0.68, N=4)

| tag | NT | η | OVERLAP | PPS | rlxf | time | run if |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `p018_b0` (done) | 36 | η\* | 2.0 | 4 | 0.3 | — | — |
| `p018_nt72` | 72 | 2·η\* | 2.0 | 2 | 0.15 | 48 h (or 2×24 h segments) | always (unless 12950996 ≈ 0.0713 AND Ryan accepts the legacy-σ pin as sufficient — then this is the single confirmation) |
| `p018_nt144` | 144 | 4·η\* | 2.0 | 1 | 0.075 | 2–4 × 48 h segments (restart-chained) | only if |CT̄(72) − CT̄(36)| > 0.5% |

Launcher defaults assume Das\* = 0.41c (η\*=1.0); override
`DAS_ETA_KINEMATIC` at submission if Phase 2 chose differently.
Harvested 12950996 (pinned at legacy σ) corroborates whichever way it landed
(see Phase 0 branch table).

## Decision

NT\* = smallest NT with |CT̄(2NT) − CT̄(NT)| ≤ 0.5% and ε_Γ ≤ 1%. With three
rungs, Richardson-extrapolate and report the observed order. dt residual →
error-budget term 3. Expected NT\* = 36.

## Exit criteria

NT\* selected; confounded-control comparison figure data recorded; ledger
updated.

## Log

- 2026-07-31 — Phase 0 branch resolved: `p2e_nt72_das2p0_ov6` (12950996) gave
  CT̄(10–19) = 0.06852 ± 0.00157 at legacy-σ pinning — **−3.9% below the NT=36
  partner (0.07133)**, below both anchors. So: `p018_nt72` runs unconditionally,
  and NT=144 will almost certainly be needed (|Δ| ≫ 0.5% at legacy σ; verify on
  the campaign pair before committing NT=144 segments). dt refinement direction
  is DOWN — the naive NT=72's +2.9% was the Das confound, opposite in sign to
  the true dt effect. Caveat carried: the 12950996 window had −1.14% 5-rev
  drift (not fully settled), so the −3.9% may be partly settling; the campaign
  `p018_nt72` should run long enough to pass M1 before the Richardson step.
  Submission still gated on Phase 2's Das\* (ladder jobs in flight).
