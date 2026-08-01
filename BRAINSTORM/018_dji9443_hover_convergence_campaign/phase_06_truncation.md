# Phase 6 — Truncation-Plane Confirmation

**Objective:** confirm the recorded 4R-vs-6R null (|ΔCT| ≤ 0.22%, both
formulations, legacy σ) at the campaign's converged settings, so the
truncation claim doesn't rest on a different discretization.

## Case (needs Phase 4's σ\*)

| tag | knobs | time |
| --- | --- | --- |
| `p018_trunc6` | final settings (Das\*, NT\*, σ\*, N=4) + TRUNCATION_DEPTH_R=6; pass final-settings env at submission (ops_reference) | 24–48 h (σ\*-dependent; warm start from the L1/final field if the Phase-4 pilot validated) |

## Decision

|ΔCT̄| ≤ 0.3% vs the matched 4R run ⇒ 4R stands; delta + prior null →
error-budget term 6. Larger ⇒ run 8R before concluding (truncation was
depth-sensitive at 40_40 below 4R; the deletion-front signature is
period-scaling ripple — check the per-rev p-p too, not just the mean).

Disclose: truncation cylinder radius 1.5R is hard-coded in the driver —
fixed, not converged (error-budget term 12).

## Exit criteria

Delta recorded on matched windows; ledger updated.

## Log

(append dated entries here)
