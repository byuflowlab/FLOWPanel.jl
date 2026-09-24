# 021 FGS scalability Stage 2 — results (2026-09-24)

Job **13875511** (`p021-fgs-stage2`, m12-3-17, zen3 exclusive, submitted
2026-09-23): backfilled immediately (the 15:41-next-day estimate was
conservative) and completed **all 17 primary stages ok** plus
`backoff_verdict.txt: backoff_helps=yes` in ~2 h 50 m. It was canceled at
2:51:34 elapsed — mid-first-pair of the *conditional* j=32 extension — during
a harness discussion (Ryan's "cold" correction; see below); no primary data
was lost. **Resume job 13878514** (2 h wall, `RESUME_FROM_JOB_ID=13875511`)
submitted to finish the 4 conditional j=32 idle pairs; ok stages skip.
Analysis: `benchmark/fgs_stage2_analysis.jl` (two first-contact fixes:
skip interrupted stages with no STATUS file; loop-variable shadowing).
Run dir: `data/p021-cold-20260910/fgs-stage2-13875511` (+ `analysis/`).

## Headline: the mechanism is idle-worker lock hammering, not memory bandwidth

**Bounded backoff at j=64 is the new best configuration: 3.24–3.28 s** —
beating uncapped spin (5.9–6.2 s; paired Δ **−2.86 s**, all 3 pairs same
sign), Stage 1's cap16 champion (3.41 s), and every lower thread count
(j=16/32 ≈ 4.4 s). Safety: backoff does not regress at j=16 (−5.5%, single
block) or at j=64 w=16 (+0.1%).

| j (spin anchors, w=0) | median s | | idle A/B @ j=64 | s |
|---|---|---|---|---|
| 16 | 4.399 | | spin (pair medians) | 5.91–6.25 |
| 32 | 4.387 | | **backoff** | **3.24–3.28** |
| 64 | 6.197 | | cap16 spin anchor | 3.533 |

Per-worker aggregates (aggdiag arms):

| config | per-task lower µs | busy share | lockmgmt share | idle share |
|--------|-------------------|------------|----------------|------------|
| j=16 w=0 | 66.1 | 0.253 | 0.018 | 0.637 |
| j=32 w=0 | 89.8 | 0.128 | 0.047 | 0.748 |
| j=64 w=0 | 161.0 | 0.060 | 0.064 | 0.823 |
| j=64 w=16 | **67.3** | 0.260 | 0.023 | 0.619 |

Three findings, each with independent support:

1. **The per-task busy inflation is contention, not bandwidth.** Per-task
   lower-GEMV time inflates ×2.44 (66→161 µs) from j=16 to uncapped j=64 —
   but at j=64 **w=16** it returns to 67.3 µs, statistically j=16's value,
   with all 64 Julia threads still present and the same memory system. If
   node memory bandwidth were the limit, capping sweep *participation* would
   not restore per-task time. The task bodies are slowed by the co-running
   spinner army (lock-line and cache traffic), not by DRAM saturation.
2. **Idle dominates worker-time everywhere and grows with team size**
   (0.64 → 0.75 → 0.82 of team-time; lockmgmt 0.018 → 0.064). The DAG is
   narrow: even 16 workers are idle ~64% of the sweep. Extra threads
   convert almost entirely into idle spinners — the enabling condition.
3. **The backoff intervention confirms causality.** Removing only the
   idle-loop lock hammering (workers pause on a lock-free queue-length hint)
   recovers −2.86 s at j=64 — more than the entire 64-thread regression
   (Stage 1: +1.5 s vs j=32) plus most of the plateau — while changing no
   arithmetic (solutions bit-identical, verified locally 46/46).

**Overhead-gate caveat (binding):** the aggdiag-vs-anchor gate FAILS at
j=64 w=0 — but in the *negative* direction: instrumentation made runs
**8.5% faster**. Timestamping between pop attempts throttles the idle
hammer rate, i.e. the observer effect itself points at spin-rate-sensitive
contention. Per the rules, the j=64 w=0 aggregate *shares* are therefore
qualitative only (gate passes at 16/32/w16: +2.7/−2.0/+1.7%); the
mechanism conclusion rests on the backoff intervention and the w=16
per-task restoration, which are uninstrumented measurements.

Also collected: spin anchors reproduce Stage 1's shape on a second node
(plateau 1.00×, regression 0.71×); j=64 w=0 spin somewhat slower here
(6.2 vs 5.8 s) — consistent with the qhint-store binary note in the
provenance (same-binary pairs unaffected).

## Status and next steps

- Resume job 13878514 owes the j=32 spin/backoff confirmation pairs.
- **Recommendation (Ryan-gated):** adopt `dagteam_idle=:backoff` as the
  production default for `sweep_order=:dagteam` — it dominates spin at every
  measured operating point, needs no per-machine tuning (unlike cap16), and
  is arithmetic-identical. Worker cap becomes unnecessary at R4.
- Stage 3 (sampled timelines, frozen-input replay) appears unnecessary: the
  scheduling-vs-bandwidth ambiguity is resolved. Residual open question —
  whether the remaining 16→32 plateau (1.00× at spin; check backoff at 32
  from the resume job) is pure DAG-width starvation — is answerable from
  the resume data.

## Harness correction (Ryan 2026-09-23)

"Cold" was always meant as **cold-START solves** (zero initial guess, no
warm start from a previous solution) — which the harness does correctly
(`reset_cold!` zeroes the solution column before every timed solve) — NOT
fresh-process-per-arm. Fresh processes remain necessary only per thread
count (`-t` is a launch flag) and per placement policy (first-touch);
within one (j, placement), arms should share a process. With Stage 2's data
banked, a mid-campaign harness rewrite saves nothing for this stage; the
batched-arm redesign applies to any future FGS benchmark run (runtime
toggles for `dagteam_workers`/`dagteam_idle` + a multi-arm driver).
