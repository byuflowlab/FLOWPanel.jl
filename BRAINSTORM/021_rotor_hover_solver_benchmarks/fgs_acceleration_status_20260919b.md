# FGS acceleration — gate-2 rev b harvest & go/no-go (2026-09-19b)

Harvest of Slurm job **13773494** (rev b, corrected socket-0 binding; commit
`b1c6c7af`, provenance `fgs_acceleration_provenance_20260919b.md`). Clean
finish (`== done ==`, empty stderr), all 27 arms produced results. Evidence
copied to `fgs_r4_followup_evidence_20260914/replay-gate2-13773494/`
(slurm .out/.err + `results/` with topology, 26 arm logs, dag.log, numastat
snapshots). Judged by outputs per house rule.

## Verdict (summary)

- **Handoff kill switch: PASS at every team size** — 2.04 / 2.51 / 3.30 /
  4.55 µs avg at t4/t8/t16/t32 (budget < 5.8 µs). The rev-a t32 artifact is
  gone: 86,508 handoffs/solve cost 0.18–0.39 s worst case.
- **Float64 row-parallel alone: NO-GO** — best arm (owner touch, t4) streams
  4.235 s → T ≈ 7.24 s = 1.40×, below the 1.5× minimum at every ladder
  point. Measured 54.8 GB/s vs the 61.9 GB/s the minimum requires.
- **Float32 storage + row-parallel: GO for the 1.5× minimum, with margin**
  — best arm (interleave 0–3, t16) streams 2.843 s → **T ≈ 5.84 s = 1.73×**
  (H budget ≈ 0.90 s before the minimum is at risk); owner t16 gives
  T ≈ 6.10 s = 1.66×. At the top edge of the 4.5–5.9 s planning range.
- **2× design target: not yet supported by measurement** — needs stream
  ≤ 2.06 s (112.7 GB/s F64-equiv); best measured is 81.6. See the F32
  headroom note below for the identified route toward it.
- **Schedule pick: source-major row-parallel confirmed; pull-DAG not
  promoted** (inclusive-comparison rule — details below).
- **Recommendation: proceed to TASK 2 production implementation**
  (rank 1 + rank 2 together, F32 convert-on-load mandatory), pending Ryan's
  go. No further gate-2 reruns needed.

## Measured numbers (true socket 0: cpunodebind 0-3, membind 0-3)

"Useful BW" is F64-equivalent (231.891 GB / 81-sweep stream); 12 timed
sweeps, min reported; 978/1068 leaves parallel (small-bytes threshold
256 KiB). T ≈ 3.0 s remainder + stream + H; baseline colored @ j16 =
10.116 s.

| Arm | t4 | t8 | t16 | t32 |
|---|---:|---:|---:|---:|
| handoff avg µs | 2.04 | 2.51 | 3.30 | 4.55 |
| rowpar F64 serial-touch GB/s | 29.6 | 25.1 | 24.2 | 24.0 |
| rowpar F64 owner GB/s | **54.8** | 46.2 | 47.1 | 49.6 |
| rowpar F64 interleave GB/s | — | — | 48.9 | 53.2 |
| rowpar F32 serial-touch GB/s | 41.5 | 37.4 | 40.3 | 34.3 |
| rowpar F32 owner GB/s | 58.0 | 61.1 | **74.8** | 68.1 |
| rowpar F32 interleave GB/s | — | — | **81.6** | 81.4 |

Serial 1-core references: F64 26.9 GB/s (8.635 s stream), F32-storage 33.4
F64-equiv (6.949 s). Consistent with rev a and the 29.4 GB/s NUMA-diagnostic
single-core arm.

Projected solve times for the candidates:

| Candidate | 81-sweep stream | T (H=0) | Speedup | vs gates |
|---|---:|---:|---:|---|
| F64 rowpar owner t4 | 4.235 s | 7.24 s | 1.40× | fails minimum |
| F32 rowpar owner t16 | 3.102 s | 6.10 s | 1.66× | passes minimum |
| F32 rowpar interleave t16 | 2.843 s | 5.84 s | 1.73× | passes minimum; misses 2× |

Medians sit within 1–10% of mins (worst spread on F32 owner t16,
0.0383→0.0423); the interleave arms are the tightest (0.0351/0.0355).

## Decision rules, itemized

1. **Handoff kill switch: PASS.** Full ladder under the budget with the
   whole socket available; t32's rev-a 26 ms value was pure
   oversubscription, as diagnosed. Matches local M2 (2.52 µs).
2. **First-touch: the comparison is now real, and placement dominates.**
   Serial-touch F64 lands all pages on the coordinator's node and pins the
   team at 24–30 GB/s (single-node ceiling); owner touch reaches 47–55.
   For F32, interleave *beats* owner touch (81.6 vs 74.8): with halved
   bytes, per-worker row tiles more often straddle sub-page boundaries, so
   ownership placement degrades exactly as the spec's sub-page-tile trap
   warned, while round-robin interleave keeps all 8 socket channels evenly
   loaded. **Caveat: every `numastat_<arm>.txt` is empty** (header only,
   PID `(null)` — the snapshot fired after the short-lived process exited
   or captured the wrapper). Placement verdicts therefore rest on the
   bandwidth deltas, which are unambiguous, not on page counts. Fix the
   capture in any future driver rev; not worth a rerun.
3. **Projected solve time: GO at 1.5×, not at 2×.** Numbers above. The
   0.90 s H-budget must absorb: initial nonself setup stream (~1 extra
   sweep, ~0.035 s), the final convergence FMM pass, and any implementation
   overhead beyond the replayed team coordination (which IS included in the
   measured sweeps). Comfortable but not lavish.
4. **Source-major vs pull-DAG: source-major.** The dag arm (priced with
   corrected serial B = 26.9 GB/s, h = 2.04 µs) saturates at its
   byte-weighted span by W=8: 0.0217 s/sweep forward, backward fully
   absorbed as filler, 1.76 s/81-sweep lower stream. But that schedule
   implies ~132 GB/s aggregate (2.863 GB/sweep in 0.0217 s) — 2.4× beyond
   the best *actual* byte rate any real kernel demonstrated on this socket
   (54.8 GB/s, F64 owner t4). Applying the measured rate as a bandwidth
   cap, a perfect DAG schedule is bandwidth-bound at ≈ 0.052 s/sweep —
   the same as the rowpar owner-t4 sweep already measured. The pull-DAG
   offers no inclusive win with margin at equal precision/placement; its
   remaining hypothetical (concurrent leaves spreading traffic across
   more channels) is unproven and costs a scheduler. Keep source-major;
   revisit only if the production kernel stalls below these replay numbers.
5. **Anomaly check (rowpar vs synthetic NUMA bench): explained, flagged.**
   Real kernels plateau at 55 GB/s actual vs the 164.4 GB/s synthetic
   affine ceiling. The replay serializes leaf-by-leaf: mean work per leaf
   is 2.68 MB (168 KB per worker at t16), so each handoff starts a short,
   latency-exposed burst rather than a long stream, and only the nodes
   owning that leaf's pages are active at any instant. The synthetic bench
   streamed large blocks on all 4 nodes concurrently and independently.
   This is a structural property of the lexicographic sequence, not a bug
   — and it is the quantitative reason 2× is out of reach today.

## The identified headroom (route toward the 2× design target)

The F32 arms move only ~41 GB/s of **actual** bytes at their best while F64
arms demonstrate ~55 on the same placement machinery. The F32 kernel is
convert/latency-limited, not bandwidth-limited. If kernel tuning (wider
convert-on-load vectorization, tile sizing, page-aligned tiles to rescue
owner placement) brings F32 actual rates to the demonstrated 55 GB/s, the
stream drops to ~2.12 s → T ≈ 5.1 s ≈ 2×. Second lever: leaf-size retuning
upward (spec rank 4, separately reported) directly attacks the short-burst
limit in item 5. Neither is assumed in the GO numbers above.

## Ready for TASK 2 (Ryan-gated)

Production implementation per the recommendation: persistent adaptive
row-parallel execution of the source-major cache + Float32 nonself storage
with Float64 arithmetic, interleave placement as the initial default at
t16 (owner-with-aligned-tiles as the tested alternative), lex sequence and
RHS semantics unchanged, gate order correctness → numerical (independent
evaluator, BC rel-L2 ≤ 1e-6) → end-to-end interleaved A/B ≥ 1.5× under
full campaign ceremony.

Open items unchanged: OpenBLAS sub-range dgemv bit-identity on the cluster
build (certify via accuracy if it fails); two-system FGS construction
BoundsError at `c18e4b46` (blocking for multi-rotor); Ryan-pending ledger
(origin pushes incl. `p021-fgs-accel-20260918`, 4 notebook entries,
WeakKeyDict/warm-start fix).
