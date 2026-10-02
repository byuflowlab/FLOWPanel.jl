# FGS acceleration — gate-2 harvest & verdict (2026-09-19)

Harvest of Slurm job **13773451** (zen3 m12-3-26, exclusive, qos=test, clean
finish, empty stderr). Provenance: `fgs_acceleration_provenance_20260919.md`.
Evidence copied to
`fgs_r4_followup_evidence_20260914/replay-gate2-13773451/` (slurm .out/.err +
26 result files).

## Verdict (summary)

- **Handoff kill switch: PASS** — 1.81 / 2.18 / 2.92 µs avg at t4/t8/t16
  (budget < 5.8 µs). 86,508 handoffs/solve ≈ 0.16–0.25 s. The t32 point
  (26,028 µs) is an oversubscription artifact, not a scheduler property (see
  below).
- **Bandwidth / projection gate: INCONCLUSIVE — the run is invalidated by a
  binding bug.** Every arm ran under `numactl --cpunodebind=0 --membind=0`,
  which on this NPS4 EPYC 7763 is **one NUMA domain: 16 CPUs and 2 of 8
  socket-0 memory channels**, not socket 0 (nodes 0–3, 64 CPUs). The driver's
  own comment ("socket-0 pinned") shows the intent; the earlier NUMA
  diagnostic that measured 74/164 GB/s used cpubind **nodes 0–3**
  (`numa_placement_findings_20260918.md`). No arm could access the placement
  regime the decision rules assume, so rejecting source-major on these
  numbers would be unsound — as would accepting it.
- **Recommendation: one corrected re-run of the same driver** (same job
  shape, <1 h qos=test) with true socket-0 binding. Needs fresh submission
  approval. No production solver work should start until it reports.

qos=test is **not** the limiter (checked at Ryan's request): sacct shows
AllocCPUS = 128 (the full 2×64-core node; ReqCPUS 64, exclusive) and the
test QOS caps are cpu=256/user, cpu=512 total, no per-node cap. The core
starvation was entirely our own numactl flag.

## Measured numbers (as run — 16-CPU / 2-channel domain)

"Useful BW" is F64-equivalent (231.891 GB / 81-sweep stream time); 12 timed
sweeps, min reported. Model: T ≈ 3.0 s remainder + stream + H.

| Arm | Best ladder point | Useful BW | 81-sweep stream | T (H=0) |
|---|---|---:|---:|---:|
| serial F64, 1 core | — | 27.0 GB/s | 8.588 s | 11.6 s |
| serial F32-storage, 1 core | — | 33.2 GB/s | 6.987 s | 10.0 s |
| rowpar F64 (owner ≈ serial touch) | t4 | 33.1 GB/s | 7.01 s | 10.0 s |
| rowpar F32 serial touch | t16 | 55.8 GB/s | 4.156 s | 7.16 s |
| rowpar F32 owner touch | t4–t16 | 51.8–52.6 GB/s | 4.41–4.49 s | ~7.4 s |
| any t32 arm (incl. interleave controls) | t32 | 0.11 GB/s | ~2,030 s | artifact |

Ladder shape: F64 rowpar is flat/non-scaling (t4 33.1 → t8 28.7 → t16
29.7–30.5 GB/s); F32 arms plateau at ~52–56 GB/s F64-equivalent = ~26–28
GB/s of actual coefficient reads — i.e. both precisions sit at the same
actual byte rate, pinned near the 2-channel domain ceiling with the
non-coefficient traffic (products, strengths, scatter) on top. That is the
signature of a placement-limited run, not a kernel or coordination failure.

Consequences of the bad binding, itemized:

1. **t32 collapse** (~25 s/sweep, all four t32 arms + both interleave
   controls): 33 spin-wait threads (team + coordinator) on 16 CPUs.
   t16 (17 threads on 16 CPUs) is also mildly oversubscribed.
2. **Bandwidth ceiling**: 2 of 8 channels → the 74 GB/s (serial-touch,
   parallel consume) and 164 GB/s (affine) regimes from the NUMA diagnostic
   were physically unreachable.
3. **First-touch comparison nullified**: `--membind=0` forces every page to
   node 0, so serial vs owner touch measured the same placement (logs agree
   to ~1%).
4. **Interleave control nullified**: `--interleave=all` arms ran only at
   t32, so they inherited the oversubscription artifact.

## What the run DID establish

- Handoff cost on zen3 passes the kill switch with 2–3× margin at sane team
  sizes (matches the local M2 2.52 µs); the persistent spin-team machinery
  works on the cluster Julia (1.12) build.
- The F32 convert-on-load kernel delivers its expected ~2× traffic
  reduction (stream 7.0 → 4.16 s at fixed actual byte rate; ratio 1.69 with
  non-coefficient traffic unhalved) with no anomalies.
- Structural checks all exact on the cluster inputs (md5-verified census/
  edges): 1,068 leaves, 2,862,850,032 B/sweep, 48,167 lower edges, unit
  path 279, byte work/span 2.849.
- dag mode ran (priced B=27.0 GB/s, h=1.81 µs): forward makespan saturates
  at the span limit by W=8 → 1.74 s/81-sweep lower stream, backward fully
  absorbed as filler at W≥8; W=4 gives 2.25 s with +0.0044 s/sweep drain.
  **Caution:** the list-schedule sim prices W workers at B each with no
  shared-bandwidth cap (W=8 implies 216 GB/s aggregate), so on this binding
  it is structural evidence only. The source-major vs pull-DAG pick is
  deferred to the corrected run, per the spec's inclusive-comparison rule.

## Projection sanity (why we still expect a GO after the fix)

At the decision thresholds (H=0): minimum 6.744 s needs 61.9 GB/s F64 or
31.0 GB/s actual-F32 (62 F64-equiv); halving needs 112.7 / 56.3. The
corrected socket has measured ceilings 74 GB/s (serial first-touch,
parallel consume) and 164.4 GB/s (affine) from the NUMA diagnostic. F32 +
even the conservative 74 GB/s placement projects T ≈ 3.0 + 115.9/74 + H ≈
4.6 s + H — inside the 4.5–5.9 s planning range. But per the spec these
remain projections; the corrected gate-2 must measure them.

## Proposed corrected re-run (Ryan-gated)

Same driver, one edit pass in `benchmark/fgs_sequence_replay_orc.slurm.sh`:

1. All socket arms: `--cpunodebind=0-3 --membind=0-3` (true socket 0).
2. Interleave controls: `--interleave=0-3 --cpunodebind=0-3`, and run them
   at the best ladder point determined by the rowpar arms (not hardwired
   t32).
3. Keep ladder t∈{4,8,16,32} (fits 64 cores/socket with the coordinator);
   optionally add t48.
4. Add per-arm `numastat -p $$` capture to verify page placement (spec trap:
   never assume placement).
5. Price dag with the corrected serial B and t4 h, and note the sim's
   missing bandwidth cap when reading its output (or add a `--dag-bw-cap`).

Cost: identical job shape (m12 exclusive, qos=test, <1 h). After it
reports: go/no-go + source-major vs pull-DAG pick per the unchanged
decision rules.

## Open items unchanged

- OpenBLAS sub-range dgemv bit-identity on the cluster build (gate-1 trap)
  — still unchecked; certify via accuracy if it fails.
- Two-system FGS construction BoundsError at `c18e4b46` — pre-existing,
  must be resolved or declared before multi-system production use.
- Ryan-pending ledger: origin pushes (incl. `p021-fgs-accel-20260918`),
  4 notebook entries, WeakKeyDict/warm-start fix.
