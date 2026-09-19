# Provenance: p021 v23 dagteam numerical gate + end-to-end A/B (2026-09-19d)

Ryan approved launching both the R4 numerical-gate job and the end-to-end A/B
campaign (2026-09-19, this session). Spec gates 3 and 4 of
`fgs_acceleration_recommendation_20260918.md`; measured basis =
`fgs_acceleration_status_20260919c.md` (champion dagteam + F32full +
interleave 0-3 @ t16 → projected 2.10×); implementation state =
`fgs_acceleration_status_20260919d.md` (gate-1 28/28, FLOWPanel plumbing).

## Pins (written before submission; enforced in-process by CAMPAIGN_PINS)

Campaign root: `/home/rander39/campaigns/p021-r4-dagteam-20260919-v23/`
(worktrees + `env/` + `pins.toml`). All tags annotated; worktrees clean;
harness verifies loaded path, HEAD sha, and tag^{commit} per package.

| Package | Worktree | Tag | SHA |
|---|---|---|---|
| FLOWPanel | `<root>/FLOWPanel.jl` | `campaign/p021-r4-dagteam-exec-20260919-v23` | `4b4418f76697bfe2887056486ee94c50d4eb5104` |
| FastMultipole | `<root>/FastMultipole` | `campaign/p021-r4-dagteam-fm-20260919-v23` | `f4d6b671b3afc5f7404e0d9b324b09bf77ddde43` |
| FLOWVPM | `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl` | `campaign/p021-cold-exec-20260910-v1` | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` |

Lineage: FLOWPanel exec commit = source tag
`campaign/p021-r4-dagteam-source-20260919-v23` (`894ba2d`, = `fastmultipole`
branch: plumbing `8aa5511` + v23 campaign assets `894ba2d`) + the standard
data-symlink commit. FastMultipole tag commit = branch
`p021-fgs-accel-20260918` (`c18e4b46` → `b0946c36` → … → `e904e763` →
`29a55bf4` production dagteam → `f4d6b671` gate-1 oracle fix; gate-1 28/28 at
-t4 locally). Refs pushed to the orc clones as `p021-r4-dagteam-v23`
(FLOWPanel) and `p021-fgs-accel-20260918-v23` (FastMultipole); origin pushes
remain Ryan-pending. Env: `env/{Project,Manifest}.toml` copied from the v22
campaign env with the three dev paths rewritten to the v23 worktrees
(no dependency changes in either commit); in-job `cold_precompile.jl`
validates it.

## The two jobs

Both submitted from `<root>/FLOWPanel.jl` (worktree top level), env
`COLD_PROJECT=<root>/env`, `CAMPAIGN_PINS=<root>/pins.toml`,
`COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`
(consolidated data root; run dirs `dagteam-numgate-<jobid>/`,
`dagteam-v23-<jobid>/`).

1. **Numerical gate** (`benchmark/run_r4_dagteam_numgate.slurm.sh`, qos=test,
   ≤1 h, zen3 500G exclusive): AB_MODE=calibrate only — staircase-calibrates
   the colored twin AND the dagteam twin (DAGTEAM_PRECISION=f32full) of the
   retained R4 lex config at j16 under the gate-2d champion placement
   (`numactl --interleave=0-3 --cpunodebind=0-3`). Every accepted solve
   passes the independent evaluator at BC rel-L2 ≤ 1e-6 — a rung that cannot
   certify fails the job (early accuracy verdict). Fallback ladder on
   failure: resubmit with DAGTEAM_PRECISION=f32conv, then f64
   (Ryan-approved ladder from 20260919c; resubmission itself needs no new
   approval per the approved ladder, but a FAILED f32full verdict is
   reported to Ryan first).
2. **End-to-end A/B** (`benchmark/run_r4_dagteam_ab.slurm.sh`, qos=normal,
   10 h, zen3 500G exclusive): controls (driver smoke, parse, precompile,
   benchmark j1/j4, FLOWPanel solver+history units, FastMultipole
   solve/coloring suites + standalone `fgs_dagteam_gate1_test.jl`) →
   self-contained calibrate at j16-interleave → uninstrumented alternating
   trials in three arms: **j16-interleave** (decision arm: champion count +
   champion placement), **j16-native** (taskset physical-core control,
   anchors the colored baseline to its v21 accepted environment),
   **j64-interleave** (robustness) → one activity-instrumented solve per
   order (attribution only). COLD_AB_REPS=10, COLD_AB_BATCHES=4 → 80
   trials/arm-order.

## Decision rules (fixed before results)

- Numerical gate: dagteam f32full staircase must certify (evaluator ≤ 1e-6)
  with a positive calibrated tolerance. FAIL → drop to f32conv rung.
- End-to-end: decision arm = j16-interleave `ab_summary.toml`.
  `dagteam_speedup_median = colored_median / dagteam_median` must be ≥ 1.5
  (acceptance); design target ≥ 2.0. Iteration counts, evaluator
  certification, and repeat-agreement (≤1e-8) are enforced per trial by the
  harness. j16-native colored median is the anchor to v21's accepted
  10.116 s — a large deviation flags an environment shift, judged before any
  promotion claim. Judge by outputs, never sacct status.
- Baseline placement fairness: the decision arm runs BOTH orders under
  interleave; colored's serially-first-touched pages can only benefit from
  interleave relative to its single-node native placement, so the candidate
  is not flattered by the shared environment.

## Local pre-submit gates (this session, M2, ≤4 threads)

- `test/runtests_r4_dagteam_ab_driver.jl`: PASS (parse, dynamic init, pair
  check incl. dagteam_precision/chunks rejections, alternation, `bash -n`
  both launchers).
- Full R4 calibrate (both twins, f32full) at j4 in the scratch env vs the
  FastMultipole worktree: launched detached; result recorded in the status
  file when harvested.

## Submission record

Submitted 2026-09-19 from the campaign worktree top level with
`DAGTEAM_PRECISION=f32full` and the env of §jobs. `sbatch --test-only` at
submission: both start immediately (2026-09-19T12:11:43 on m12-2-11, zen3).
slurm-availability probe: m12 8 idle fitting nodes, `--qos=normal` ETA
immediate.

| Job | Script | JobID | ETA at submission |
|---|---|---|---|
| numerical gate | run_r4_dagteam_numgate.slurm.sh | **13773687** (qos=test) | 12:11 m12-2-11 |
| A/B campaign | run_r4_dagteam_ab.slurm.sh | **13773689** (qos=normal) | 12:11 m12-2-11 |

Storage preflight: /home at 399 G of the 400 G FLOWPanel cap at submission,
with an hpc-storage archive pass RUNNING (≈213 GB of finished p018 runs being
tarred; deletes land per verified run). Submitted under the cap on the
judgment that v23 writes ~1 G of CSV/TOML (no VTK) over hours while the
archiver reclaims two orders of magnitude more in the same window; a
transient ~400.5 G peak is possible before the first archive delete lands.
Storage cycle results recorded in the status file when the pass reports.
