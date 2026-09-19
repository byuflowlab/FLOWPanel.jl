# NUMA-placement investigation — findings (2026-09-18)

Executes `numa_bandwidth_reset_prompt_20260918.md` (fastest path).
Submission of the two diagnostic jobs pre-approved by Ryan 2026-09-18.

## Step 0 — nearfield-cache assembly: threaded or not? (ANSWER: single-threaded)

Checked locally in FastMultipole `flowpanel-20260817` @ `c18e4b46`
(`/Users/ryan/Dropbox/research/projects/FastMultipole`), BEFORE submitting.

The FGS solver's streamed influence cache is `nonself_matrices`, built by the
`FastGaussSeidel` constructor at `src/solve.jl:673` via
`nonself_influence_matrices(...)`:

- **Allocation** — `Matrices(sizes, TF)` (`src/solve.jl:239`) allocates one
  contiguous backing vector through `_calloc_vector` (`src/solve.jl:6-11`),
  which uses `Libc.calloc`. For GB-scale requests glibc serves calloc from
  fresh mmap'd zero pages, so **no page is touched at allocation** — NUMA
  placement is decided entirely by the first write.
- **Fill** — the population loop at `src/solve.jl:256-321`
  (`for (i_target, j_source) in sorted_list ... this_matrix[:,isb] .= this_influence`)
  is a plain sequential loop: no `Threads.@threads`, no `@spawn`. Every
  element of the cache is written by the single constructing thread.
- Same story for `self_influence_matrices` (`src/solve.jl:427-440`,
  `Matrices` at :439): sequential fill.
- The *parallel* assembly path that does exist
  (`NearfieldInfluenceCache`, `src/nearfield_cache.jl:455-478`,
  `@spawn` workers over an atomic chunk pool) is a **different object**,
  not used by the `FastGaussSeidel` constructor — no `NearfieldCacheDonor`
  involvement in the FGS build path.

**Conclusion:** first-touch is single-threaded → under the node's default
first-touch policy all ~2.86 GB of `nonself_matrices.data` lands on the ONE
NUMA node where the constructing thread ran. This is the scenario under
which the placement hypothesis is maximally likely (one memory controller,
~40–50 GB/s local ceiling, 64 threads starved — matches every v22 number).

## Step 1 — standalone dgemv microbenchmark (go/no-go gate)

Design per the handoff: 1,068 Float64 matrices of 650×550
(1068 × 650 × 550 × 8 B = 3.056 GB ≈ the v20 census 2.86 GB; per-block
2.86 MB), per-matrix x/y vectors, BLAS threads = 1, `Threads.@threads
:static`, v22 pinning (first 64 physical cores socket-major → cpubind
nodes 0–3, first-touch default policy), one exclusive zen3 (m12) node.

Arms (median of 20 timed full passes after 5 warm-ups; GB/s = bytes/median):

| arm | threads | first-touch | mempolicy | predicts |
|-----|---------|-------------|-----------|----------|
| a | 1 | single-thread | default | ~29 GB/s (lex analogue) |
| b | 64 static | single-thread | default | ~40–50 GB/s (chunked v22 analogue) |
| c | 64 static | single-thread | `numactl --interleave=0-3` | ~4× (b) |
| d | 64 static | parallel chunk-affine (each thread allocates+fills its partition) | default | ≥ (c) |

Page placement is SHOWN, not inferred: `numastat -p <pid>` +
`/proc/self/numa_maps` captured inside every arm between fill and sweep.

Files: `benchmark/numa_dgemv_bench.jl` + `benchmark/run_numa_dgemv.slurm.sh`
(diagnostic, not an official campaign — no worktree/pin ceremony; pure
Julia + BLAS, touches no repo code). Outputs (CSV + logs only) →
`~/projects/FLOWPanel.jl/data/p021-cold-20260910/numa-bench-<job>/`.

**Gate:** (c) or (d) ≥ ~3× (b) → placement CONFIRMED → step 2 (in-situ j64
interleave pair). < 1.5× → REFUTED → stop and report.

### Results (job 13763831, m12 exclusive, 2026-09-18; ran in ~4 min)

Output: `~/projects/FLOWPanel.jl/data/p021-cold-20260910/numa-bench-13763831/`
(`summary.csv`, per-arm pass CSVs, `numastat_arm{a..d}.txt`,
`numa_maps_arm{a..d}.txt`, topology captures; COMPLETED marker present,
judged by outputs). 3.065 GB streamed per pass, median of 20 passes.

| arm | threads | first-touch / policy | median s/pass | GB/s | vs (b) |
|-----|---------|----------------------|---------------|------|--------|
| a | 1 | serial, default | 0.1044 | **29.4** | — |
| b | 64 | serial, default | 0.0414 | **74.0** | 1.00× |
| c | 64 | serial, `--interleave=0-3` | 0.0199 | **153.8** | 2.08× |
| d | 64 | chunk-affine, default | 0.0186 | **164.4** | 2.22× |

Page placement SHOWN by `numastat -p` (private MB by NUMA node):

- arm b: **2934 of 3151 MB (93%) on node 1** — the single-node pileup the
  hypothesis predicted, directly observed.
- arm c: 751/750/750 MB even across nodes 1–3 (+898 node 0 incl. runtime).
- arm d: 746/700/736 MB spread across nodes 1–3 (+988 node 0).
- arm a: split 1407/1751 across nodes 0–1 (serial process migrated once).

### Verdict: mechanism CONFIRMED, magnitude BELOW the step-2 gate → STOP

- Arm a lands exactly on the bandwidth-math prediction (29.4 vs ~29 GB/s):
  the kernel is bandwidth-bound as modeled.
- Single-thread first-touch demonstrably concentrates the cache on ONE
  NUMA node (arm b numastat), and fixing placement (interleave or
  chunk-affine) roughly **doubles** delivered bandwidth, with chunk-affine
  first-touch the best arm — the *mechanism* is real.
- BUT the *magnitude* is 2.08–2.22×, not the hypothesized ~4×, because
  both ends of the ratio moved: one node under 64-thread load delivers
  74 GB/s (not the assumed 40–50 local ceiling), and the interleaved
  arms at 154–164 GB/s already sit at the practical socket ceiling
  (~160–200 GB/s). There is no 4× to be had within cpubind 0–3.
- Caveat: 74 GB/s from one NPS4 node exceeds a naive 2-channel DDR4
  budget (~51 GB/s), suggesting kernel auto-NUMA balancing may have
  migrated some pages during the timed sweeps (numastat was captured
  before the sweeps). This does not change the decision: any such help
  is equally present in the production configuration, so 2.2× remains
  the achievable improvement over the status quo.
- **Gate outcome:** 2.22× is below the ≥~3× CONFIRMED threshold (and above
  the <1.5× REFUTED threshold). Step 2's submission pre-approval was
  conditioned on ≥~3×, so the in-situ j64 interleave pair was NOT
  submitted. Stopping here per the handoff.

### Implications (for Ryan's decision)

- Best-case in-situ projection at 2.2×: nonself products 0.292 → ~0.13
  s/it, chunked span/iter 0.277 → very roughly 0.15–0.18 s → chunked j64
  ≈ 9.5–11 s at 44 iterations — at best marginal against the colored@j16
  operating point (10.116 s), and the real kernel ran at only ~37 GB/s
  aggregate vs the analogue's 74, so scatter/barrier overheads likely eat
  part of the 2.2× too.
- Placement fixes are therefore NOT the standalone lever the v22
  postmortem hoped for; they are a ~2× bandwidth multiplier that only
  pays if iteration inflation (27→44) is also addressed.
- Levers on the table (all Ryan-gated): **Float32 nearfield storage**
  (already approved; halves streamed bytes, multiplies with placement),
  chunk-affine first-touch assembly in FastMultipole (the code change;
  arm d shows it beats interleave and needs no numactl wrapper),
  iteration reduction for chunked, per-chunk timers/uncore counters.

## Step 2 — in-situ j64 interleave pair

NOT RUN — step 1 gate not met (2.22× < ~3×). See verdict above.
