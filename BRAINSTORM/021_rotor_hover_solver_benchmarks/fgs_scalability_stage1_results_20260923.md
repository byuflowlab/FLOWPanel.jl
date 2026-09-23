# 021 FGS scalability Stage 1 — results (2026-09-23)

Job **13858983** (`p021-fgs-stage1`, m12-3-29, zen3 exclusive 128c/500G)
COMPLETED 2026-09-23 08:30 MDT after 12:46 elapsed, `failed_count=0`,
**45/45 stages ok** (ladder 12, ladderdiag 12, accepted 3, placement 6,
cap16 6, cap32 6 — `cap16_verdict.txt: cap16_helps=yes` triggered the cap32
extension). Run dir:
`/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/fgs-stage1-13858983`
(now contains `analysis/` with `stage1_report.md`, `stage1_all_rows.csv`
315 rows, `stage1_process_medians.csv` 45 groups). Deploy tree
`/home/rander39/campaigns/p021-fgs-stage1-20260922/` (ARCHIVER_SKIP); pins
in `fgs_scalability_stage1_provenance_20260922.md`.

Analysis script `benchmark/fgs_stage1_analysis.jl` had its first real
execution this session; two fixes landed (this commit):

1. **Soft-scope bug**: `overhead_ok = false` inside the top-level `for`
   created a *local*, so the instrumentation-overhead gate could never
   FAIL. Fixed with `global overhead_ok`. (Gate passes on the real data
   anyway — see below — so no conclusion changes.)
2. **cap32 table was never rendered** — the script only printed a message
   when cap32 pairs were *absent*. Cap sections now loop over cap∈{16,32}
   with paired-mean summaries.

## Headline: both effects reproduced; both localize to `nonself_product`

Fixed-work ladder (uninstrumented, medians of block medians):

| j | median s | speedup vs j=1 | ratio vs prev rung |
|---|---------|----------------|--------------------|
| 1 | 31.754 | 1.00 | — |
| 16 | 4.425 | 7.18 | 7.18 |
| 32 | 4.295 | 7.39 | 1.03 |
| 64 | 5.791 | 5.48 | 0.74 |

- **16→32 plateau REPRODUCED** (1.03× of the 2.00× ideal).
- **64-thread regression REPRODUCED** (64/32 time ratio 1.35).
- Instrumentation overhead gate **PASS** (|overhead| ≤ 1.2% at every rung,
  shape unchanged) → phase decomposition is admissible.

Phase decomposition (matched means, instrumented rows only; recon 100.0%
at every rung — exclusive phases fully account for diag_total):

| phase | T(16) | T(32) | T(64) |
|-------|-------|-------|-------|
| fmm | 2.027 | 1.317 | 1.013 |
| nonself_product | 2.006 | 2.651 | 4.511 |
| initialization | 0.229 | 0.230 | 0.228 |
| residual | 0.127 | 0.125 | 0.130 |
| everything else | ≤0.02 each | | |

FMM keeps scaling through 64 threads. **`nonself_product` (the dagteam
near-field sweep product) anti-scales**: +0.645 s from 16→32 (eating the
FMM gain −0.710 s → net plateau) and +1.860 s from 32→64 (vs FMM's −0.303 s
→ net regression). dagteam spawn/join/wait/reduce splits are all ≤0.025 s —
the loss is inside the sweep work itself, not in team orchestration.

## Placement A/B @ j=64 (3 pairs, alternating arm order)

Champion (socket-0) 5.98–6.15 s vs both-socket 9.40–11.07 s; paired mean
Δ = **+4.43 s**, all pairs same sign. Champion placement confirmed decisively.

## Worker-cap A/B @ j=64 (champion placement; paired medians only)

| cap | workers=0 (baseline) s | capped s | paired mean Δ | same sign |
|-----|------------------------|----------|----------------|-----------|
| 16 | 5.695 / 6.168 / 5.773 | 3.392 / 3.410 / 3.437 | **−2.466 s** | yes |
| 32 | 5.729 / 6.469 / 6.047 | 3.860 / 3.934 / 3.947 | **−2.168 s** | yes |

Capping the dagteam near-field sweep at 16 workers while j=64 gives
**3.41 s median — the fastest configuration in the whole study**, beating
the uncapped j=64 baseline by ~2.5 s, the j=32 ladder (4.30 s), and the
j=16 ladder (4.43 s): FMM keeps all 64 threads while the sweep stays at
its efficient width. cap16 also beats cap32 (3.41 vs 3.91 s) in every pair.

Binding caveat (plan C trap): the cap win localizes the 64-thread loss to
**sweep participation** but does NOT distinguish DAG width vs queue
contention vs memory traffic; worker cap is a participation test, not a
topology test.

## Accepted-accuracy bridge

j=16/32/64 accepted solves: 4.518 / 4.439 / 5.805 s, all 27 iterations —
tracks the fixed-work ladder (accepted runs one extra
fmm+influence+residual), tying the new pins to the Stage-0 verified ladder.

## Status

Stage 1 babysit+harvest COMPLETE. Both plan-C gates passed, decomposition
delivered, single-lever answer: **the near-field sweep's parallel width is
the whole story** — plateau and regression share the mechanism, and a
16-worker cap at j=64 is the current best operating point. Stage 2/3
remain **Ryan-gated**; no further HPC submissions from this thread. The
run dir becomes archive-eligible for the storage flow (deploy tree stays
ARCHIVER_SKIP).
