# Gate-2 sequence replay benchmark — provenance (2026-09-19)

Diagnostic go/no-go benchmark (NOT a campaign acceptance run; Ryan approved
submission 2026-09-19 and lightweight ceremony — dev-branch worktree + recorded
SHAs — was recommended in `fgs_acceleration_status_20260918.md` and applied.
Final accepted numbers will come from the later end-to-end A/B under full
campaign ceremony.)

## Job

- Slurm job **13773451**, submitted 2026-09-19 ~08:00 from
  `/home/rander39/wt-p021-fgs-gate2` (worktree of the orc
  `~/projects/FastMultipole` clone), `-p m12 --qos=test`, 1 node exclusive,
  64 cpus, 100 G, 59 min walltime.
- ETA probe: m12/test start-now (`slurm_availability.py --cpus 64 --mem-gb
  100 --time 00:59:00 --eta`, 2026-09-19).

## Code pin

- FastMultipole branch `p021-fgs-accel-20260918`, HEAD **`a6492f49`**
  (parent `b0946c36` = gate harnesses; branched off Manifest pin `c18e4b46`).
  Pushed by direct git push to the orc clone (origin push blocked: https
  remote, token is Ryan's); worktree checked out at `a6492f49`, verified.
- Driver: `benchmark/fgs_sequence_replay_orc.slurm.sh`; benchmark:
  `benchmark/fgs_sequence_replay.jl` (modes serial / rowpar / handoff / dag).
- Julia: cluster default via `module load julia` (1.12), stdlib-only script.

## Inputs

R4 census/edges (j64-b1 arm of diag-v15-13694724), read from the original
cluster outputs:

- `~/projects/FLOWPanel.jl/data/p021-cold-20260910/diag-v15-13694724/j64-b1/results/gemv_census.csv`
  — md5 `927a9b9988305febf52d448f0127a318`
- `.../dependency_edges.csv` — md5 `40f7f2d22bc096ae85769777e2e2e020`

Both verified identical to the local evidence copies in
`fgs_r4_followup_evidence_20260914/diag-v15-13694724/j64-b1/results/`
(same md5s, checked 2026-09-19).

## Outputs

- `~/wt-p021-fgs-gate2/slurm-p021-fgs-replay-gate2-13773451.{out,err}`
- `~/wt-p021-fgs-gate2/benchmark/replay_gate2_13773451/` (topology, per-arm
  logs, dag pricing)

## Matrix

serial {F64,F32} ×1 core; handoff ladder t∈{4,8,16,32}; rowpar
{F64,F32}×{serial,owner first-touch}×t∈{4,8,16,32}, all socket-0
`numactl`-pinned; interleave control at t=32; dag mode priced with this
node's measured serial B and t4 handoff h. 12 timed sweeps per arm (min
reported), OPENBLAS_NUM_THREADS=1 throughout.

## Decision rules (from the recommendation)

- Handoff kill switch: avg < ~5.8 µs.
- Reject candidates whose measured budget cannot reach 6.744 s with margin
  (T ≈ 3.0 + coeff_GB/B + H; 231.891 GB F64 / 115.945 GB F32 per 81-sweep
  solve).
- Same benchmark prices the pull-DAG promotion condition (inclusive
  schedule comparison at equal precision/placement).
