# 033 A-R1 — Standalone cooperative-leaf-product kernel microbenchmarks (2026-09-26)

**Verdict: column partitioning wins decisively at every measured leaf shape;
row partitioning is rejected for short-wide leaves (it is *slower than
sequential* at w=2 on p25–p90 shapes). Measured team overhead is ~0.1 µs
(barrier) / ~0.4–3 µs (effective, median-class leaves) — a factor ~5–80 below
the 16–32 µs make-or-break range A-T2 swept (only the 1,450-tail h_4 ≈ 24 µs
enters that band, and it still clears its own threshold), putting the
program on the optimistic side of the A-T2 projections on this host.**

## What was run

`benchmark/fgs_coop_kernel_microbench.jl` (new, this session): times the leaf
lower-product y = L·x (the `dagteam_pull!` GEMV) sequentially vs cooperatively
under a persistent spin-barrier team (workers reused across reps; no per-leaf
task spawning — the dagedge anti-pattern), in the two A-T1 layouts:

- **row** — disjoint output row ranges, one BLAS GEMV per worker on a
  row-block view;
- **col** — contiguous column blocks into private partial vectors, then a
  fixed ascending-order reduction by the timing worker (worker 1).

Leaf shapes are real `(nof, ptot)` pairs from the R4 DAG
(`033_atheory_20260926/fgs_dag_L_graph_R4.csv`, 1,068 leaves), selected at nof
quantiles plus the mean-54 leaf and the 1,450 tail; matrices are synthesized
dense Float32 at those shapes (GEMV cost depends on shape, not content),
matching the R4 f32full champion (`szTM=szTS=4`). Per-shape LU `ldiv!` (the
leaf self-solve) is timed for the critical-path model. A zero-work "null"
shape isolates the pure dispatch+barrier round-trip. Timing: `time_ns()` per
product, min/median/mean over 200–20,000 reps after warmup (021 ruling 5/7).

**Environment:** Ryan's local machine (`tmplab-32-117-31.et.byu.edu`), Julia
1.12.5, `-t 4`, BLAS pinned to 1 thread (verified by the common.jl probe;
env-pinned, see "BLAS mystery resolved" below). FLOWPanel `fastmultipole`
@ `c866f6f`-dirty, FastMultipole `flowpanel-20260817` @ `745af760`-dirty.
Local screening run per standing rulings (≤4 threads local); **not** campaign
evidence — HPC re-measurement happens in A-R2.

Evidence: `033_ar1_20260926/fgs_coop_kernel_microbench.csv` (59 rows),
`banner.txt`, `run.log`. Judge from the CSV.

## Results (median µs per product)

| shape | nof×ptot | seq | row w=2 | row w=4 | col w=2 | col w=4 | col w=4 speedup | LU ldiv |
|---|---|---|---|---|---|---|---|---|
| p25 | 22×5836 | 10.4 | 13.1 | 14.0 | 5.6 | **3.3** | 3.19× | 0.08 |
| p50 | 39×5369 | 12.7 | 18.1 | 16.3 | 6.9 | **5.3** | 2.38× | 0.17 |
| mean54 | 54×5889 | 18.4 | 28.3 | 17.3 | 9.6 | **7.5** | 2.44× | 0.21 |
| p75 | 59×5885 | 22.0 | 28.8 | 21.9 | 13.0 | **7.8** | 2.84× | 0.25 |
| p90 | 83×4819 | 24.6 | 33.8 | 29.1 | 12.7 | **8.7** | 2.82× | 0.58 |
| p99 | 410×2920 | 83.1 | 58.2 | 41.3 | 43.5 | **24.9** | 3.34× | 12.3 |
| max | 1450×2878 | 327.9 | 200.5 | 142.8 | 169.2 | **105.6** | 3.10× | 140.9 |

Barrier-only round-trip (null shape): w=2 **0.083 µs**, w=4 **0.125 µs**
(median). Team machinery at w=1 costs nothing measurable (row/col w=1 ≈ seq
within ~8%; where they differ, w=1 is *faster* — the seq median itself carries
that much run-to-run noise).

## Derived overheads and the A-T3 threshold

Effective elapsed team overhead h_w = t_coop − t_seq/w (includes dispatch,
barrier, reduction, load imbalance, and bandwidth contention — the complete
A-T3 elapsed path):

| shape | h_2 (col) | h_4 (col) | A-T3 threshold h_4/(1−1/4) | c_prod (seq) | split wins? |
|---|---|---|---|---|---|
| p25 | 0.40 | 0.66 | 0.9 µs | 10.4 µs | yes |
| p50 | 0.52 | 2.16 | 2.9 µs | 12.7 µs | yes |
| mean54 | 0.42 | 2.94 | 3.9 µs | 18.4 µs | yes |
| p75 | 1.94 | 2.24 | 3.0 µs | 22.0 µs | yes |
| p90 | 0.38 | 2.56 | 3.4 µs | 24.6 µs | yes |
| p99 | 1.96 | 4.10 | 5.5 µs | 83.1 µs | yes |
| max | 5.2 | 23.6 | 31.5 µs | 327.9 µs | yes |

Note: with h_w *defined* as t_coop − t_seq/w, clearing c_prod > h_w/(1−1/w) is
algebraically equivalent to the direct comparison t_coop < t_seq — this table
is A-T3's practical rule ("measure the complete cooperative product against
the sequential product and split only when it is faster") in threshold form,
not an independent screen. The max-leaf h_4 = 23.6 µs is the one measured
value inside A-T2's 16–32 µs band (bandwidth-dominated at that shape); its
327.9 µs product clears the 31.5 µs threshold by 10×.

On this host, **every leaf shape down to p25 clears the A-T3 split threshold
at w=4 in the column layout** — the empirical answer to "which leaves split"
is "essentially all with ptot>0," consistent with A-T3's warning not to
pre-restrict to hot leaves. Mapping to A-T2's swept overhead parameter c:
measured effective h_w ≈ 0.4–3 µs (median-class) sits between A-T2's c=0 and
c=8 µs grid points, i.e. near the 3.21× sweep-ceiling-bound end of the
projection band rather than the 1.26× (c=64 µs) end.

## Findings beyond the headline

1. **Row partitioning rejected for short-wide leaves.** Row w=2 is *slower
   than sequential* for p25–p90 (e.g. mean54: 28.3 vs 18.4 µs). Cause: a
   row-block view is strided (lda = full nof), each worker still traverses
   the full 5,000+-column x, and BLAS falls off its fast path. The plan's
   caution ("do not assume row partitioning wins") was warranted — inverted:
   column wins everywhere, including the tall-thin tail (max: 105.6 vs 142.8).
2. **A-T1 "row = bitwise" does not survive BLAS.** The contract's condition
   (same per-row accumulation order/blocking) is violated by BLAS itself:
   sgemv on a row-block view picks different internal blocking than on the
   full matrix. Bitwise agreement vs sequential is shape- and w-dependent
   (CSV `bitwise` column: row w=2 bitwise at p50/p75/p90 but not p25/mean54).
   Consequence: with BLAS kernels, BOTH layouts must be accepted under the
   column-layout clause of A-T1 (deterministic, exact-arithmetic equivalent,
   time-to-certified-accuracy) — bitwise reproducibility would require
   hand-rolled fixed-order kernels. Determinism (run-to-run identical for
   fixed w/partition) held in every config tested.
3. **Serial LU is the tail-leaf bottleneck at w=4.** For the 1,450 leaf, the
   serial `ldiv!` (140.9 µs) now *exceeds* the w=4 cooperative product
   (105.6 µs) — direct confirmation of A-T2's prediction that LU's
   critical-path share rises (8%→26%) and concrete motivation for Track C
   (exact triangular solver / cooperative trsv).
4. **The A-T2 "8 BLAS threads" mystery is resolved.** On this host, the
   OpenMP-backed OpenBLAS resets a runtime `BLAS.set_num_threads(1)` back to
   8 on first real work; `common.jl`'s probe correctly hard-errors. The A-T2
   session evidently worked around it by setting `BENCH_BLAS_THREADS=8`
   rather than pinning at load time. Fix (used here): set
   `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` in the launch environment.
   This closes the open provenance flag on the A-T2 run — it was a host
   quirk plus a workaround, not an unexplained protocol choice.

## Caveats

- Single host, laptop-class, 4 Julia threads; the w=4 team ran alone on an
  otherwise idle process. In production at j64, ~16 teams of 4 share memory
  bandwidth; h_w and the GEMV throughput itself will differ. These numbers
  gate the **layout choice** (column) and give the **first real overhead
  scale** (µs, not tens of µs); the decision-grade h_w comes from A-R2's
  in-situ prototype on HPC.
- The team kernel closure is dynamically dispatched (`Ref{Any}`), a small
  upper-bias on h_w; a type-stable executor integration can only do better.
- Synthetic matrices at real shapes; admissible for GEMV/ldiv timing (cost is
  shape-driven), inadmissible for anything numerical beyond the determinism/
  tolerance checks performed.
- An earlier version of the harness deadlocked (spin loops without
  `GC.safepoint()`/yield escape — 2h14m hang, killed). Fixed with
  GC-cooperative spinning + rare bounded yields (never on the µs hot path
  once workers own their threads). Lesson recorded for the A-R2 executor:
  **spin waiters inside FastMultipole must call `GC.safepoint()`**
  (the existing dagteam executor already parks via locks, which is safe).

## Consequence for A-R2

Build the prototype cooperative executor with the **column layout** (private
partials + fixed ascending reduction by the publishing worker), team width
w=4 default with w=2 fallback, spin barriers with GC-safepoints, layout/width
as runtime toggles (`dagteam_workers`), and measure h_w in situ on the R4
problem at j-ladder points before any A-T2 re-projection.
