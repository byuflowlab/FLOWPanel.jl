# Gate-2c/2d benchmarks (full-F32 + real split executor) — provenance (2026-09-19)

Ryan approved 2026-09-19: "I approve hpc jobs to test them both (source-major
+ F32 first, split second without waiting for the first to finish)", plus the
full-F32 question ("If we store everything as F32, would that speed things up
further by avoiding the convert to F64 cost?"). Same lightweight diagnostic
ceremony as gate-2 revs a/b (dev-branch worktree + recorded SHAs; accepted
numbers come later from the end-to-end A/B under full campaign ceremony).

## What these jobs test

- **Gate-2c (job 13773580)**, `benchmark/fgs_replay_gate2c_f32full_orc.slurm.sh`:
  new `--f32-full` precision mode — F32 storage AND state/accumulate/LU, no
  convert-on-load — across serial baselines and the rowpar ladder
  (owner t4–t32, interleave t16/t32), with the rev-b convert-mode champions
  rerun in-job as same-node references. Speed only: full-F32 numerics are a
  different numerical experiment (F32 eps 1.2e-7 vs the 1e-6 BC gate) and
  must separately pass the independent evaluator before production use.
- **Gate-2d (job 13773581)**, `benchmark/fgs_replay_gate2d_dagteam_orc.slurm.sh`:
  new `dagteam` mode — a REAL (not priced) split dual-layout executor:
  target-major lower pulls driven by readiness counters over the 48,167
  directed lower edges (scan-pop priority queue, byte-weighted critical-path
  priority), source-major backward upper products as lower-priority filler
  tasks, target-owned boundary reduction into the next sweep's upper
  accumulator (frozen u^s never overwritten mid-sweep). Ladder t4–t32 F64 +
  t16/t32 at F32conv/F32full, under interleave 0-3 binding (ownership is
  dynamic, owner-touch undefined), with in-job rowpar references at matching
  precision/placement — the spec's inclusive equal-precision comparison.

Local smoke (M2, ≤4 threads, machinery validation only): all modes run;
dagteam structural invariants exact (1,068 lower tasks, 3 roots, 48,167
edges, lower/upper bytes 1,513,294,744 / 1,349,555,288 — matches dag.log).

Driver fix carried in both: numastat capture now polls every 15 s and keeps
the last live snapshot (rev b's single 60 s shot raced process exit; all
rev-b snapshots were empty).

## Jobs

- **13773580** (gate2c) and **13773581** (gate2d), submitted 2026-09-19 from
  `/home/rander39/wt-p021-fgs-gate2`, both `--qos=test`, 1 node exclusive,
  64 cpus-per-task, 100 G, 59 min (same shape as revs a/b). They do NOT run
  concurrently: qos=test enforces MaxJobsPerUser=1 (13773581 pends with
  QOSMaxJobsPerUserLimit and starts automatically when 13773580 finishes).
  Verified at submission: 13773580 RUNNING on m12-2-13 with sane serial
  numbers (F64 27.8 GB/s, F32conv 34.0, F32full 52.2), 13773581 queued.
- Both submitted with `CENSUS`/`EDGES` env overrides (the rev-b trap),
  md5-verified at submission: census `927a9b9988305febf52d448f0127a318`,
  edges `40f7f2d22bc096ae85769777e2e2e020`.

## Code pin

- FastMultipole branch `p021-fgs-accel-20260918`, HEAD **`e904e763`**
  (parent `b1c6c7af` = rev-b HEAD; adds `--f32-full`, `dagteam`, the two
  drivers; generalizes the replay state to a separate storage/state type).
- Synced to orc by pushing ref `p021-fgs-gate2c` to the orc clone via ssh
  (https origin push still Ryan-pending) and ff-merging the worktree;
  verified at `e904e763`, tracked files clean.
- Julia: cluster default via `module load julia` (1.12).

## Outputs (expected)

- `~/wt-p021-fgs-gate2/slurm-p021-fgs-replay-gate2c-13773580.{out,err}` and
  `.../slurm-p021-fgs-replay-gate2d-13773581.{out,err}`
- `~/wt-p021-fgs-gate2/benchmark/replay_gate2c_13773580/` and
  `.../replay_gate2d_13773581/` (topology, per-arm logs, numastat snapshots)

Harvest beside the rev a/b evidence as
`fgs_r4_followup_evidence_20260914/replay-gate2c-13773580/` and
`replay-gate2d-13773581/`.

## Decision rules

- Same T ≈ 3.0 s remainder + stream + H model, baseline colored @ j16 =
  10.116 s, minimum ≤ 6.744 s, design ≤ 5.058 s; 231.891 GB F64 /
  115.945 GB F32 per 81-sweep solve; report useful BW as F64-equivalent.
- Gate-2c question: does f32full beat f32conv at equal placement? If yes and
  it approaches the 2.06 s stream (112.7 GB/s F64-equiv), the 2× design
  target comes into measured reach — but adoption is CONDITIONAL on the
  independent accuracy evaluator at 1e-6, with f32conv as the certified
  fallback and selective per-block F64 retention as the fallback's fallback.
- Gate-2d question (spec rule): does dagteam beat rowpar at equal
  precision/placement with margin enough to justify the split rebuild?
  Compare best-vs-best and same-team-size; check the dagteam sweep against
  the 2.849 work/span structural bound and the priced sim's 1.76 s
  lower-stream claim (which assumed uncapped aggregate bandwidth).
