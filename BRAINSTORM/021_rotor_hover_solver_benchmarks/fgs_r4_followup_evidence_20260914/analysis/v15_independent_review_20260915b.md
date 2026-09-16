# Independent v15 numerical and scaling review — 2026-09-15

All six locally harvested arms of job 13694724 were reviewed programmatically. This is a numerical/internal-consistency review; transfer hashes and executed-source provenance are independently pending with the harvester. Root COMPLETED and six complete arm tables exist. No cluster access or rerun was performed.

## Numerical and structural findings

All 240 timed rows, 12 direct equivalence controls, and six profile validations pass the inspected numerical criteria. Certified FMM BC rel-L2 is 4.77952513e-7 throughout; direct controls give 4.78150591e-7 and direct/FMM disagreement 5.4676859e-9. Reported repeat/control solution differences are zero. Both instrumentation histories have length 28 and recorded exact equality in every arm. The actual histories and solution vectors were not exported, so equality cannot be independently recomputed from these CSVs.

Instrumented counts are **28 outer residual checks, 27 update iterations, 81 inner sweeps, 86,508 leaf visits**. The terminal converged residual check explains why solver.niter=27 differs from outer_count=28. Do not label both as 27 outer calls.

The v15 driver does not explicitly export or assert `all(isfinite, x)`. However, its required finite zero norm-based repeat difference plus finite certified residual strongly supports finite solutions under the executed standard norm semantics: NaN/Inf entries would make the repeat ratio nonfinite and fail `<=1e-8`. This is source-backed evidence, not an independent vector audit. Profile memory values are below the 500-GiB gate; per-trial memory values were checked in the driver but omitted from diagnostic CSVs.

## Scaling and instrumentation robustness

| Julia threads | Uninstrumented median s | Raw instrumented difference % | Difference excluding identified first-instrumented row % |
|---:|---:|---:|---:|
| 1 | 37.55940 | -0.241 | -0.259 |
| 4 | 16.93579 | +0.890 | +0.887 |
| 8 | 13.46425 | -0.594 | -0.642 |
| 16 | 12.19939 | -2.148 | -2.246 |
| 32 | 11.36727 | +0.385 | -0.090 |
| 64 | 10.96265 | -0.131 | -0.231 |

Each arm has four sequential batches U/I/U/I, ten trials each. All raw rows remain retained; exclusion is a clearly labeled sensitivity calculation, not replacement evidence. Broad scaling is substantial: j1→j64 is 3.426× and j4→j64 is 1.545×. The j32→j64 median difference is only 0.405 s (3.56% less wall time), with overlapping ranges and no independent process/allocation replication; do not claim a robust small gain or establish a crossover.

**Warmup artifact:** batch 2 trial 1 in every arm has an outer-solve-minus-internal-timer gap of 0.752–0.781 s. Source warms `diagnostics=nothing`; the instrumented history control has a callback, while measured instrumented solves have no callback. A first-use specialization/JIT cost is therefore a supported hypothesis, not definitively measured attribution. All other gaps are below 0.1 s. This artifact invalidates describing all instrumented rows as fully warmed. Do not rerun the completed ladder merely to remove it; report raw results and sensitivity.

The raw median differences range from −2.148% to +0.890%, not a consistent positive timer tax. Matched adjacent batch differences change sign at j4 (−0.208%, +1.807%) and j64 (−1.575%, +2.620%); j16 remains approximately −1.25%/−2.26%. Noise, drift, allocation/GC, and order confounding prevent interpreting negative values as a speedup or asserting zero instrumentation overhead. Instrumented timings are descriptive stage diagnostics only. Future measured generations should warm the exact no-callback instrumented specialization.

The counter driver baseline uses the same callback type and `stage_observer=nothing` as its later counted solve, so this particular missing-specialization warmup does not automatically transfer to that counter interval.

## Stage budget and placement

Per-trial exclusive-stage sums reconcile to the reported exclusive total, and exclusive+unaccounted reconciles to the internal total, with maximum absolute CSV-rounding discrepancy 9.300001e-8 s. Internal unaccounted work is approximately 0.012–0.013 s. The outer gap, including its first-instrumented artifact, must be reported separately; medians alone conceal it. Independent medians of components need not sum exactly to the median total; per-trial checks are authoritative.

FMM median decreases from 26.735 s at j1 to 6.878 s at j4 and 1.114 s at j64. Products remain roughly 7.8–8.2 s. At j64, median per-trial shares are approximately 72% products, 10% FMM, 8% scatter, 5% leaf solves, and 3% initialization. The old ~1.3% FMM main-task profile share is not its wall fraction. Full numerical per-stage medians and ratio medians are in the stage CSV.

All selected CPUs are on socket 0. j1/4/8/16 use NUMA node 0; j32 spans nodes 0–1; j64 spans 0–3. Memory policy is default/current, with all eight memory nodes allowed. This is a physical-core cpuset, not fixed worker-to-core binding or controlled page placement. First-touch follows setup execution; actual page placement is not measured. Cross-arm memory/NUMA behavior remains a possible contribution. v15 activity spans reset and validation, so it cannot establish stage-specific active-thread counts. Workload cache/DRAM counters are still absent.

## Block/dependency census and next experiment

All three census files are byte-identical across arms. Independently recomputed invariants: 1,068 leaves; 2,862,850,032 Float64 matrix bytes; 5,431,340 scatter entries; 95,390 unique directed dependencies; 48,627 undirected conflicts. Every block satisfies bytes=8*m*n. Independently recomputed directed degrees and symmetrized conflict degrees match the table. All 79 hypothetical colors are conflict-free with exact recorded membership/byte totals. Sizes range 1–27 (median 16); 64-way leaf concurrency is unavailable in a single color. These are independent checks against the exported graph, not a fresh reconstruction from serialized solver geometry.

Actual dimensions: m min/median/max=2,199/4,708.5/19,892; n=1/39/1,450; block bytes=28,136/1,401,636/80,399,352. Percentiles and extremes are in the census CSV; large block imbalance matters for colored scheduling.

The stage budget supports the governing request's **existing colored-sweep experiment**, since the serial product/scatter/leaf sweep dominates j64. Treat it as separately calibrated ordering, not a drop-in numerical-equivalent optimization. First run calibration/acceptance controls at useful thread counts (j4 and j64 are comparable endpoints), then matched uninstrumented accepted timing if controls pass. Record colors, load balance, synchronization/scatter costs, iterations/sweeps, and total solve time. A modest useful-concurrency arm can be motivated by the 27-leaf maximum, but no new broad ladder is warranted. Preserve failed calibration outcomes. Fewer seconds per sweep alone is not success.

Finish stage activity/cache diagnostics before assigning a memory mechanism. Do not claim DRAM saturation or infer speedup from Float32 storage. Slice fixes, Float32 pilots, and new scheduling algorithms remain subsequent isolated experiments. No notebook write or silo deletion is authorized by this review itself; follow the existing evidence/terminal-job cleanup prerequisites.

## Existing audit limitations

The saved audit is useful but PASS alone does not cover explicit finite vectors, raw histories, per-trial memory output, actual counts, census arithmetic/conflicts, or first-instrumented outer gaps. This independent analysis covers exported counts/census/gaps and inspects source-backed finite/history/memory evidence. Its fixed 1e-6-s additive bound is stricter than the observed printed errors but does not substitute for the original audit's decimal-aware tolerances. The original `written_precision` helper ignores scientific-notation exponents; this can over-relax tolerances for negative exponents (or understate them for positive ones), so retain explicit observed residuals rather than relying on that helper generically. Existing direct equivalence controls are properly mandatory and pass; timed/profile NaN direct metrics correctly mean unevaluated.

## Reproduction and artifact index

Run `python3 BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_r4_followup_evidence_20260914/analysis/v15_independent_review.py` from repository root. It writes only these new derived analysis CSVs:

- `v15_independent_scaling.csv`: raw medians/ranges, speedup, sensitivity, first-instrumented outer gap.
- `v15_independent_batches.csv`: all 24 batch medians/ranges and ordering.
- `v15_independent_stage_budget.csv`: exclusive stages, internal total/residual, median per-trial shares.
- `v15_independent_checks.csv`: numerical thresholds, count checks, reconciliation residuals.
- `v15_independent_census.csv`: dimensions/bytes/dependencies/color-size totals and percentiles.

Inputs are `../diag-v15-13694724/j*-b1/results/*.csv`, census TOMLs, and root topology/affinity plus per-arm NUMA records. Source reasoning uses clean local FLOWPanel v15 `benchmark/fgs_r4_diagnostics.jl` and `fgs_cold_common.jl`. Governing task: `../../fgs_opt_r4_diagnostics_handoff_20260912.md`, 2026-09-14 request.
