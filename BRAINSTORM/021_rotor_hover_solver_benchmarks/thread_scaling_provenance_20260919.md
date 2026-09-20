# Provenance: p021 preliminary R4 thread-scaling study (2026-09-19)

Ryan approved (2026-09-19, this session): commit, sync to HPC, and launch the
R4 thread-scaling campaign, plus the R1–R2 re-runs (recorded separately).
Purpose: before adopting the full machine-class memory ladder
(0/16/32/64/128/500 GiB with per-class thread caps at ~2 GiB/thread, cap 64;
budget 0 uncapped) as general 021 benchmark policy, measure the thread
efficiency of the two production solver families at R4 across
j = 1/8/16/32/64 so that plateaued thread tiers can be pruned from the
general policy. Ryan's decision input: "If one or both plateau at a certain
thread count, that might inform my decision by making some tiers
unnecessary."

## Pins (written before submission; enforced in-process by CAMPAIGN_PINS)

Campaign root: `/home/rander39/campaigns/p021-thread-scaling-20260919/`
(FLOWPanel worktree + `env/` + `pins.toml`). FastMultipole and FLOWVPM reuse
the v23/v1 worktrees unchanged (verified clean at their pinned SHAs before
reuse).

| Package | Worktree | Tag | SHA |
|---|---|---|---|
| FLOWPanel | `<root>/FLOWPanel.jl` | `campaign/p021-thread-scaling-exec-20260919` | `28eef244949c8c936062fa92196d015ec05aeb14` |
| FastMultipole | `/home/rander39/campaigns/p021-r4-dagteam-20260919-v23/FastMultipole` | `campaign/p021-r4-dagteam-fm-20260919-v23` | `f4d6b671b3afc5f7404e0d9b324b09bf77ddde43` |
| FLOWVPM | `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl` | `campaign/p021-cold-exec-20260910-v1` | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` |

Lineage: FLOWPanel exec commit = source tag
`campaign/p021-thread-scaling-source-20260919` (`db4f053`, on branch
`fastmultipole` = v23 wrap `363e29a` + this study's commit) + the standard
data-symlink/site-policy commit. `db4f053` carries: machine-class thread caps
(`machine_max_threads = min(64, GiB/2)`, budget 0 uncapped), tuning identity
extended to (rung, budget, julia_threads) for rows/resume/traces, the
PHASE2_OUTDIR override (tuner + measurement driver), a fix for the
silently-broken row-level resume ("16" vs "16.0" tag mismatch — pre-existing;
standby requeues may have appended duplicate rows to cluster tune_phase2.csv
files, check when harvesting), and `benchmark/run_r4_thread_scaling.slurm.sh`.
Env: `env/{Project,Manifest}.toml` copied from the v23 campaign env with the
FLOWPanel dev path rewritten to this worktree (no dependency changes).
Branch pushed to the orc unified repo as `p021-thread-scaling-20260919`;
origin pushes remain Ryan-pending.

## The job

`benchmark/run_r4_thread_scaling.slurm.sh`, submitted from the worktree top
level with `COLD_PROJECT=<root>/env`, `CAMPAIGN_PINS=<root>/pins.toml`,
`COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`
(consolidated data root; run dirs `thread-scaling-j<J>-<arrayjobid>/`).

Array 0–4 → j = 1/8/16/32/64, one exclusive 128-core zen3 500G node each,
qos=normal, 48 h. Placement = the promoted champion's
(`numactl --interleave=0-3 --cpunodebind=0-3`, BLAS 1) in every measured
stage; `HARDWARE_TAG=orc-m12-zen3-socket0-ilv0-3-blas1` keeps these
rows/traces unmixable with historical both-socket phase-2 data. Per task:

1. **FGS arm** — retained champion knobs (P8/MAC0.4/leaf100/inner3);
   `fgs_r4_dagteam_ab.jl` AB_MODE=calibrate staircase-calibrates the colored
   and dagteam (`DAGTEAM_PRECISION=f32full`) twins AT THIS j (every accepted
   solve evaluator-certified BC rel-L2 ≤ 1e-6 — tolerances never carry across
   thread counts), then AB_MODE=trials (4 batches × 10 per order). The
   colored arm is the fallback-ladder scaling curve for free.
2. **iLU-GMRES arm** — FMM apply knobs re-descended at this j by
   `rotor_hover_solver_phase2_tune.jl` at budgets 0 (uncached endpoint,
   seed 10:0.6:6) and 500 (node, seed 15:0.55:32), 10 h backstop per budget,
   writing to the per-task `PHASE2_OUTDIR`; then measured as
   `krylov_ilu,krylov_ilu_nfcache` by `rotor_hover_solver_phase2.jl` reading
   knobs from the same directory. Retune is REQUIRED, not optional: since
   krylov_ilu's last R4 numbers (2026-08-25, fm `4c0f1b8f`), the FMM apply
   path changed materially (030 block assembly merge, one-division kernel,
   parallel NF cache build + donor retarget — fm `4c0f1b8f..f4d6b671`).

Per-arm verdicts land as `STATUS_*` files; `COMPLETED` means the task ran to
the end (judge by outputs, never sacct). Known risk, accepted: dagteam at
j=1 is uncovered by any gate — a failure there is a scaling finding recorded
by `STATUS_fgs_calibrate=FAILED`, and the iLU arm still runs.

## Submission

Submitted 2026-09-19: **job 13777133** (array _0–_4 → j=1/8/16/32/64), after
a clean `sbatch --test-only` (13777129, immediate start on m12). Local
pre-submit smoke of the new tuner paths ran at R1/j4 on macOS (cap-skip,
landed-skip, PHASE2_OUTDIR all exercised; zero duplicate rows on rerun).

## Decision rules

- Scaling curves = median accepted solve time vs j per family (FGS dagteam,
  FGS colored, krylov_ilu, krylov_ilu_nfcache), each point at its own
  calibrated tolerance/knobs and certified ≤ 1e-6.
- A family "plateaus at j*" if the j > j* points improve by less than ~10%
  per doubling; tiers whose caps exceed every family's j* are candidates for
  pruning from the general ladder policy (Ryan decides).
- No physics null needed: these are solver benchmarks (no production physics
  chain consumes FGS — 2026-09-19 inventory).
