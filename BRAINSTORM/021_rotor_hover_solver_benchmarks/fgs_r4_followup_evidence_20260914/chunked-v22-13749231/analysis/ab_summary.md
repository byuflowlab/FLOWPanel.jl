# v22 chunked-sweep A/B — results (job 13749231, 2026-09-18)

Job COMPLETED (10:59 h wall, m12-3-31 zen3 exclusive, 64c, BLAS=1;
2026-09-17 23:21 → 2026-09-18 10:20 UTC). `COMPLETED` sentinel present;
all six stage `status.toml` = completed; harvest 103/103 files
SHA256-verified against remote (`remote-sha256-full.txt`). All cluster
controls PASS including `fgs_chunked_test.jl` 1723/1723.

## Gates (§6 of the plan / diagnostics package) — ALL GREEN

Every one of the 320 trials (40/order/arm, 4 arms): accepted, eligible,
solved, finite, uninstrumented, certified-FMM authoritative, BC rel-L2
≤ 1e-6 (lex 4.780e-7, chunked 7.868e-7), repeat-solution delta **exactly
0 on every trial**. j1 warmups carry direct-vs-FMM equivalence per order
(evaluator delta 5.47e-9 both orders). Iteration count is arm-invariant
per order: lexicographic 27 everywhere, chunked 44 everywhere — the
chunked order's determinism claim is verified, not assumed.
Calibration: chunked tolerance 5.348662427942506e-7 certified in 44
iterations (lex retained 3.479e-7 / 27), confirmation repeat delta
exactly 0.0 — matches the local §7 gate crossing to the last digits
(local 5.348662428097527e-7; ARM↔zen3 libm). Cross-order solution
rel-L2 3.12e-6 (informational; both orders separately certified).

## Ranking — uninstrumented solve_seconds medians (total time to accepted accuracy)

| Arm | Lex median s | Chunked median s | Chunked/Lex | Note |
|---|---:|---:|---:|---|
| j1 | 38.392 | 60.757 | 1.583 | ≈ 44/27 iteration ratio (1.630); per-iteration serial cost ~equal (1.381 vs 1.422 s/it) |
| j4 | 17.328 | 24.335 | 1.404 | |
| j16 | 12.478 | 15.893 | 1.274 | |
| j64 | 10.951 | 14.264 | **1.303** | chunked's best arm; min 13.822 |

(IQRs tight throughout: e.g. j64 chunked [14.115, 14.452], j64 lex
[10.761, 11.214]; j64 lex min 10.627 reproduces the 10.96–11.20
yardstick band's floor.)

**Headline: chunked LOSES at every arm.** Best certified chunked
operating point = j64 median **14.264 s**, vs both yardsticks:

| Yardstick | s | Chunked best vs it |
|---|---:|---:|
| lex @ j64 (v21/v15) | 10.96–11.20 | +27–30% slower |
| colored @ j16 (v21, global best) | **10.116** | **+41% slower** |

Per §5 of `fgs_chunked_hybrid_plan_20260918.md`, the coloring revert is
**NOT executed** — chunked's best certified point does not beat
colored@j16. Coloring KEEPS its production role; this outcome is
reported to Ryan as a finding for a fresh decision.

## Activity attribution (j64 pair; budget only, excluded from rankings)

Per-stage /proc activity (complete endpoints only, CLK_TCK=100):

| Order | Stage | Σ span s | busy CPU s | avg active threads | span/iteration s |
|---|---|---:|---:|---:|---:|
| lex | nearfield_update (27) | 9.33 | 9.32 | 1.00 | 0.346 |
| chunked | nearfield_update (44) | 12.17 | 472.33 | **38.81** | **0.277** |
| lex | fmm (28) | 1.24 | 34.46 | 27.77 | |
| chunked | fmm (45) | 1.83 | 54.86 | 29.94 | |

Unlike colored (42.65 avg active, span per iteration UNCHANGED), chunked
both engages the threads (38.81 avg active) **and shrinks the
per-iteration nearfield span ~20%** (0.346 → 0.277 s/it). The parallel
mechanism works as designed. But the win is nowhere near the Amdahl-style
expectation (~3.5–4.5 s total): per-iteration span shrank 1.25×, not the
~10–20× full-width hope — the serial deferred cross-chunk scatter plus
barrier overhead per inner sweep caps it — and the **iteration count
inflated 27 → 44 (+63%)** from the majority-Jacobi character of 64-way
chunking, more than consuming the per-iteration gain at every arm.

## Interpretation

- Convergence cost dominates: 44 vs 27 iterations is exactly the §1.4
  risk realized. At j1 the entire loss is the iteration ratio; per-
  iteration serial cost is unchanged (scatter refactor is cost-neutral).
- The parallel harvest is real but small: at j64 chunked's per-iteration
  total-solve cost is 0.324 s vs lex 0.406 s (~20% faster per iteration),
  consistent with the nearfield span attribution.
- Net: total-time-to-accepted-accuracy, the binding metric, ranks chunked
  last at every thread count. The v21 result stands: **colored @ j16 =
  10.116 s remains the global best operating point**, lex @ j64 ≈ 10.95 s
  second.
- Possible follow-ups (NEW experiments, Ryan's call per §11): fewer
  chunks (more GS-like → fewer iterations, less parallelism),
  under-relaxation of the Jacobi-like cross-chunk coupling, or the
  deterministic parallel-by-target deferred scatter (§1.3 variant).

Evidence: this directory (`calibrate/`, per-arm `ab_trials.csv` +
`ab_summary.toml`, `warmups.csv`, activity CSVs, control logs; all
SHA256-verified). Pins per
`../v22-deployment/submission-provenance-13749231.md`.
