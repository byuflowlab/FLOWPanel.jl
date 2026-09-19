# v21 colored-sweep A/B — results (job 13738665, 2026-09-17)

Job COMPLETED (10:28 h wall, m12 zen3 exclusive, 64c, BLAS=1). `COMPLETED`
sentinel present; all six stage `status.toml` = completed; harvest
103/103 files SHA256-verified against remote (`remote-sha256-full.txt`).

## Gates (§6 of `fgs_opt_r4_diagnostics_package_20260912.md`) — ALL GREEN

Every one of the 320 trials (40/order/arm, 4 arms): accepted, eligible,
solved, finite, uninstrumented, certified-FMM authoritative, BC rel-L2
≤ 1e-6 (lex 4.780e-7, colored 5.848e-7), repeat-solution delta ≤ 1e-8
(all exactly 0). j1 warmups carry direct-vs-FMM equivalence per order
(evaluator delta 5.47e-9 both orders). Iteration count is arm-invariant
per order: lexicographic 27 everywhere, colored 26 everywhere —
the colored order's determinism claim is verified, not assumed.
Calibration: colored tolerance 5.223e-7 (lex retained 3.479e-7),
confirmation repeat delta 0. Cross-order solution rel-L2 1.53e-6
(informational; both orders separately certified). Matches the local gate
exactly (26 vs 27 iterations, same tolerances).

## Ranking — uninstrumented solve_seconds medians (total time to accepted accuracy)

| Arm | Lex median s | Colored median s | Colored/Lex | Colored effect |
|---|---:|---:|---:|---|
| j1 | 38.059 | 36.766 | 0.966 | = 26/27 iteration ratio (0.963); per-iteration cost unchanged serially |
| j4 | 17.580 | 15.534 | 0.884 | −11.6% |
| j16 | 12.294 | 10.116 | **0.823** | −17.7% |
| j64 | 11.196 | 11.949 | **1.067** | +6.7% — colored LOSES at j64 |

(IQRs tight: e.g. j16 colored [10.079, 10.147]; j64 lex min 10.966
reproduces the v15 10.96 yardstick.)

**Headline: the A/B answer is split.** Colored wins decisively at j4/j16
and loses at j64. The **new global best operating point is colored @ j16 =
10.116 s**, beating the previous champion (lex @ j64, 10.96–11.20 s) by
~8–10% while using a quarter of the cores. Nowhere near the 2.2 s
Amdahl-if-free bound: the ≈6.4k color barriers/solve and the median
16-way color width cap the realizable parallelism.

## Activity attribution (j64 pair; budget only, excluded from rankings)

Per-stage /proc activity (complete endpoints only, CLK_TCK=100):

| Order | Stage | Σ span s | busy CPU s | avg active threads |
|---|---|---:|---:|---:|
| lex | nearfield_update (27) | 10.15 | 10.09 | **0.99** |
| colored | nearfield_update (26) | 11.02 | 469.97 | **42.65** |
| lex | fmm (28) | 1.18 | 34.63 | 29.44 |
| colored | fmm (27) | 1.01 | 33.59 | 33.14 |

Coloring does engage the threads (0.99 → 42.65 avg active) **but the span
does not shrink** (10.15 → 11.02 s): ~470 busy-CPU-seconds buy zero span
reduction at j64. With 79 colors of median size 16, at j64 most of the 64
threads have no leaf in the current color and spin at the barrier
(Julia threads busy-wait), and per-color parallelism is capped at the
color width. This is why j16 — matched to the median color width — is the
sweet spot: enough threads to fill a color, few left idling at barriers.

## Interpretation

- At j1 the entire colored win is the free iteration (26 vs 27) from the
  reordered iterate path; per-iteration serial cost is unchanged.
- j4→j16: genuine parallel harvest of the previously serial chain
  (beyond-iteration per-iteration speedup ~8% at j4, ~15% at j16).
- j64: barrier/spin overhead of 6.4k barriers exceeds the parallel gain;
  colored regresses 6.7% against lex.
- Structural ceiling: color-width distribution (1–27, median 16), not
  memory bandwidth, is the binding constraint at high j.

Evidence: this directory (`rundir/`, per-arm `ab_trials.csv`,
`ab_summary.toml`, `warmups.csv`, calibration outputs, activity CSVs,
control logs; all SHA256-verified). Pins per
`v21-deployment/submission-provenance-13738561.md` (resubmission record).
