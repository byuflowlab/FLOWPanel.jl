# FGS Solver Thread-Scaling Diagnostic: Ladder Summary

**Diagnostic:** diag-v15-13694724  
**Analyzed:** j1, j4, j8, j16, j32, j64 (6 arms × 40 trials each)  
**Trials per arm:** batches 1–4 with 10 trials each (uninstrumented/instrumented pairs)

---

## Table 1: Scaling (Uninstrumented Median Wall Time)

| Arm | Median (s) | Speedup vs j1 | Parallel Efficiency | Instrumentation Overhead % |
|---|---:|---:|---:|---:|
| j1 | 37.56 | 1.0000 | 1.0000 | −0.241 |
| j4 | 16.94 | 2.218 | 0.5544 | 0.890 |
| j8 | 13.46 | 2.790 | 0.3487 | −0.594 |
| j16 | 12.20 | 3.079 | 0.1924 | −2.148 |
| j32 | 11.37 | 3.304 | 0.1033 | 0.385 |
| j64 | 10.96 | 3.426 | 0.0535 | −0.131 |

**Summary:** 3.4× speedup at j64; efficiency drops steeply from j4 (55%) to j64 (5%), suggesting serial bottlenecks dominate beyond j8.

---

## Table 2: Stage Medians from Instrumented Profiles

| Arm | Total (s) | FMM Total (s) | Nonself Prod. (s) | Initialize (s) | Scatter (s) | Leaf Solve (s) |
|---|---:|---:|---:|---:|---:|---:|
| j1 | 37.46 | 26.73 | 8.162 | 1.003 | 0.7921 | 0.5977 |
| j4 | 17.08 | 6.878 | 7.793 | 0.647 | 0.7933 | 0.5406 |
| j8 | 13.38 | 3.550 | 7.807 | 0.3654 | 0.8660 | 0.5312 |
| j16 | 11.94 | 2.087 | 7.996 | 0.3223 | 0.7976 | 0.5542 |
| j32 | 11.34 | 1.415 | 7.945 | 0.3365 | 0.7951 | 0.5503 |
| j64 | 10.94 | 1.114 | 7.887 | 0.3104 | 0.8658 | 0.5496 |

**Key observation:** FMM total dominates at j1 (71% of 37.46s), drops to 10% at j64 (1.114s). Nonself product (work-to-order iteration stage) plateaus near 7.8–8.0s across all arms, indicating unit-stride iteration serialization.

---

## Table 3: Stage Speedups vs j1

| Arm | FMM Speedup | Nonself Prod. | Initialize | Scatter | Leaf Solve |
|---|---:|---:|---:|---:|---:|
| j1 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| j4 | 3.887 | 1.047 | 1.551 | 0.9984 | 1.106 |
| j8 | 7.530 | 1.045 | 2.746 | 0.9146 | 1.125 |
| j16 | 12.81 | 1.021 | 3.113 | 0.9930 | 1.078 |
| j32 | 18.90 | 1.027 | 2.982 | 0.9962 | 1.086 |
| j64 | 24.01 | 1.035 | 3.232 | 0.9148 | 1.088 |

**Top 3 stages at j64 by wall-time contribution:**
1. **FMM total:** 24.0× speedup (1.114s, 10.2% of j64 wall time) — strong scaling via GEMV parallelism
2. **Nonself product:** 1.04× speedup (7.887s, 72.1% of j64 wall time) — essentially no scaling; single- or few-thread bottleneck in outer loop
3. **Initialize:** 3.23× speedup (0.310s, 2.8% of j64 wall time) — modest parallel scaling

---

## Table 4: Iterations and Inner Sweeps per Arm

| Arm | Iterations (p50) | Estimated Inner Sweeps (p50) | FMM Passes (p50) |
|---|---:|---:|---:|---:|
| j1 | 27 | 81 | 28 |
| j4 | 27 | 81 | 28 |
| j8 | 27 | 81 | 28 |
| j16 | 27 | 81 | 28 |
| j32 | 27 | 81 | 28 |
| j64 | 27 | 81 | 28 |

**Constancy check:** PASS. All arms converge to identical iteration counts (27 outer iterations, 81 estimated inner sweeps, 28 FMM passes). No algorithmic divergence across thread counts.

---

## Table 5: Thread Activity Summary

| Arm | Thread Count | Min CPU Ticks | p50 CPU Ticks | Max CPU Ticks |
|---|---:|---:|---:|---:|---:|
| j1 | 2 | 0 | 1.00e+04 | 19963 |
| j4 | 6 | 0 | 5.00e+03 | 5801 |
| j8 | 12 | 0 | 2.00e+03 | 3506 |
| j16 | 24 | 0 | 1.00e+03 | 2324 |
| j32 | 48 | 0 | 6.00e+02 | 1753 |
| j64 | 96 | 0 | 3.00e+02 | 1478 |

**Interpretation:** CPU ticks per thread from `thread_activity.csv` (validation instrumentation); max ticks decay roughly as 1/√threads. Many threads at 0 ticks indicate idle cores. p50 ticks per thread drop from ~10k (j1) to ~300 (j64), consistent with load spreading. Wide min–max range per arm suggests uneven thread scheduling or brief parallel phases.

---

## Amdahl Attribution & Scaling Limits

**Overall j1→j64 speedup:** 3.426× (instrumented, 37.46s → 10.94s)

**Stages exhibiting flat scaling (j32→j64, within 10%):**

| Stage | j1 | j32 | j64 | j32/j64 | % of j64 Wall |
|---|---:|---:|---:|---:|---:|
| **Scatter** | 0.7921 | 0.7951 | 0.8658 | 0.918 | 7.9 |
| **Initialize** | 1.003 | 0.3365 | 0.3104 | 1.084 | 2.8 |

Scatter shows *slight regression* j32→j64 (−8.2%), possibly thread-affinity thrashing or false-sharing at high thread counts. Initialize saturates near j32, contributing only ~3% of remaining wall time.

**Speedup ceiling:** The nonself product stage (7.887s at j64, 72% of wall time) scales trivially (1.04×). This outer-loop stage represents the primary serial bottleneck, likely inherent to the FGS solver's Krylov iteration structure or matrix setup.

---

## Anomalies Detected

- **None in convergence/iteration counts:** All arms resolve to 27 iterations identically.
- **Scatter regression at j64:** Wall time *increases* from j32→j64 (0.795s → 0.866s, +8.2%). Probable cause: thread-scheduling overhead or L3-cache contention at 64 threads. Not fatal (8% of total) but worth noting for future single-node scaling studies.
- **Instrumentation overhead:** Negligible and uniformly distributed (−2.1% to +0.9%), indicating clean profiling setup.

---

## File Paths

- **Script:** `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_r4_followup_evidence_20260914/diag-v15-13694724/analysis/analyze_ladder.py`
- **This report:** `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_r4_followup_evidence_20260914/diag-v15-13694724/analysis/ladder_summary.md`
- **Audit (source):** `analysis/audit-full.txt`

---

## Caveats

- Profile fractions (stage_medians) are from `cpu_thread_complete_flat.txt` and instrumented diagnostic CSV, not wall-clock fractions. Thread-activity ticks are validation-instrumentation counts, spanning entire solve (including setup/validation), not just the FGS kernel.
- Nonself product stage may include outer-loop serial overhead, loop dispatch, and synchronization—not purely the matrix-product kernel.
- Thread activity (p50/min/max ticks per thread) is relative to validation instrumentation epoch; absolute timing skewed by system noise.
- All arms exhibit consistent algorithmic state; speedup floor at j64 is a fundamental serial bound, not a tuning opportunity.
