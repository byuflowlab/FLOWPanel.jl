# Gate-2 sequence replay benchmark rev b — provenance (2026-09-19)

Corrected resubmission of gate 2 after job 13773451 was invalidated by a
numactl binding bug (see `fgs_acceleration_status_20260919.md`). Ryan
approved the resubmission 2026-09-19 ("go ahead"); same lightweight
diagnostic ceremony as rev a (dev-branch worktree + recorded SHAs; final
accepted numbers come from the later end-to-end A/B under full campaign
ceremony).

## What changed vs job 13773451

- All team arms (handoff, rowpar) now bind `--cpunodebind=0-3 --membind=0-3`
  = true socket 0 (NPS4 EPYC 7763: each socket is 4 NUMA nodes × 16 CPUs /
  2 channels; rev a's `--cpunodebind=0 --membind=0` was 1/4 socket).
- Serial 1-core baselines keep single-node `--cpunodebind=0 --membind=0`
  (correct for one thread; matches the 29.4 GB/s reference arm).
- Interleave controls now `--interleave=0-3 --cpunodebind=0-3`, run at t16
  AND t32 (rev a ran `--interleave=all` only at the oversubscribed t32).
- Per-arm `numastat -p` snapshot 60 s into each run
  (`numastat_<arm>.txt`) — verify placement, never assume it.
- dag arm unchanged (priced from this node's serial B and t4 handoff h);
  driver comment now warns the sim has no shared-bandwidth cap (aggregate
  W×B must be sanity-checked against the measured socket ceiling).

## Job

- Slurm job **13773494**, submitted 2026-09-19 from
  `/home/rander39/wt-p021-fgs-gate2`, `--qos=test`, 1 node exclusive,
  64 cpus-per-task, 100 G, 59 min walltime (same shape as rev a; rev a's
  52 min elapsed was dominated by the six ~300 s oversubscribed t32 arms,
  which the fix removes).
- Submitted with `CENSUS`/`EDGES` env overrides pointing at the original
  cluster outputs (paths + md5s below), as rev a was.
- Dead first attempt: job **13773493** FAILED in 1 s at the input check —
  submitted without the `CENSUS`/`EDGES` overrides, and the driver's default
  path (`~/projects/FLOWPanel.jl/BRAINSTORM/.../fgs_r4_followup_evidence_20260914/...`)
  does not exist on orc (evidence commit never pulled there). No results
  produced; no code difference vs 13773494.

## Code pin

- FastMultipole branch `p021-fgs-accel-20260918`, HEAD **`b1c6c7af`**
  (parent `a6492f49` = rev-a HEAD; only
  `benchmark/fgs_sequence_replay_orc.slurm.sh` changed). Pushed to the orc
  clone as ref `p021-fgs-gate2b` (direct ssh push; https origin push still
  Ryan-pending), worktree fast-forwarded and verified at `b1c6c7af`, clean
  apart from rev-a untracked outputs.
- Benchmark script `benchmark/fgs_sequence_replay.jl` unchanged from rev a.
- Julia: cluster default via `module load julia` (1.12).

## Inputs

Identical to rev a, re-verified at this submission:
`~/projects/FLOWPanel.jl/data/p021-cold-20260910/diag-v15-13694724/j64-b1/results/gemv_census.csv`
(md5 `927a9b9988305febf52d448f0127a318`) and `.../dependency_edges.csv`
(md5 `40f7f2d22bc096ae85769777e2e2e020`).

## Outputs (expected)

- `~/wt-p021-fgs-gate2/slurm-p021-fgs-replay-gate2-13773494.{out,err}`
- `~/wt-p021-fgs-gate2/benchmark/replay_gate2_13773494/` (topology, per-arm
  logs + numastat snapshots, dag pricing)

## Decision rules (unchanged from rev a)

- Handoff kill switch: avg < ~5.8 µs (rev a already passed at t4–t16:
  1.81/2.18/2.92 µs; confirm on whole-socket binding).
- Reject candidates whose measured budget cannot reach 6.744 s with margin
  (T ≈ 3.0 + coeff_GB/B + H; 231.891 GB F64 / 115.945 GB F32 per 81-sweep
  solve); design target ≤ 5.058 s.
- Source-major vs pull-DAG: inclusive schedule comparison at equal
  precision/placement.
