# Provenance: p021 R1–R2 champion-scheme re-run (2026-09-19)

Ryan approved (2026-09-19, this session): re-run the 021 FGS characterization
under the promoted approach — re-tuned per thread count with per-thread saved
parameters, new accuracy recorded, then the new benchmarks — restricted to
**R1–R2** for now (larger rungs wait on the R4 thread-scaling study, job
13777133, which may prune thread tiers).

## What "the new scheme" means here

- **FGS family runs dagteam** via new env plumbing (`FGS_SWEEP_ORDER`,
  `FGS_DAGTEAM_PRECISION` in `benchmark/phase1_case.jl`, splatted into every
  `FGSSolver`/`FGSPreconditioner` ctor in fgstune, fgsprecond, and phase2).
  Precision = **f64**, the safe off-R4 rung (mathematically the
  lexicographic iterate); f32full only after a case's own staircase
  certifies it (resubmit with `FGS_DAGTEAM_PRECISION=f32full` later if
  wanted). Defaults everywhere remain lexicographic/f64 — adoption is per
  launcher, never silent.
- **Per-(rung, j) re-tuning**: fgstune knob descent + tolerance staircase,
  fgsprecond sweep ladder, and the phase-2 apply-knob descent all re-run at
  each thread count; nothing carries across j (021 v23 trap). The staircase
  and phase-2 `bc_rel_l2`/`bc_certified` columns are the recorded new
  accuracy.
- **Full machine-class ladder**: tuner budgets 0/16/32/64/128/500 GiB with
  the 2 GiB/thread caps (a budget whose cap < j is skipped); phase2.jl
  measures the full default CONFIGS table (backslash, krylov family, fgs,
  fgmres_fgs, nfcache variants) and skips budgets with no tuned row.
- **Champion placement** (`numactl --interleave=0-3 --cpunodebind=0-3`,
  socket 0) in every stage; `HARDWARE_TAG=orc-m12-zen3-socket0-ilv0-3`
  keeps rows/traces unmixable with historical both-socket data. **BLAS =
  julia threads** (historical phase-1/2 multi convention, kept so
  backslash_ldiv rows stay comparable; deviation from the R4 champion's
  BLAS-1 is deliberate and recorded per row in blas_threads).
- **Isolation**: per-task `BENCH_CASE_ROOT` + `PHASE2_OUTDIR` under the data
  root — the phase-1 knob CSVs have no sweep_order column and
  `stage3_winner`/`staircase_for` select on rung+knobs only, so per-run
  directories are the correctness boundary; also avoids the 2026-08-18 NFS
  concurrent-append hazard.

Not plumbed (unchanged, still lexicographic): the downstream
`*margin_verify.jl` verification scripts — not part of this pipeline.

## Pins

Own campaign root `/home/rander39/campaigns/p021-r12-champion-20260919/`
(worktree + `env/` + `pins.toml`) — a separate worktree because the
thread-scaling jobs (13777133) are running from the
`p021-thread-scaling-20260919` worktree and a live campaign worktree is
never edited (concurrency rule).

| Package | Worktree | Tag | SHA |
|---|---|---|---|
| FLOWPanel | `<root>/FLOWPanel.jl` | `campaign/p021-r12-champion-exec-20260919` | `f85ac25557c27f5724ab2149fc2d553c06cae599` |
| FastMultipole | `.../p021-r4-dagteam-20260919-v23/FastMultipole` | `campaign/p021-r4-dagteam-fm-20260919-v23` | `f4d6b671b3afc5f7404e0d9b324b09bf77ddde43` |
| FLOWVPM | `.../p021-cold-20260910-v1/FLOWVPM.jl` | `campaign/p021-cold-exec-20260910-v1` | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` |

Source tag `campaign/p021-r12-champion-source-20260919` (`d3ae35e` on
`fastmultipole`); exec = source + data-symlink/site-policy commit. Env copied
from the v23 campaign env with the FLOWPanel dev path rewritten (no
dependency changes). Branch `p021-thread-scaling-20260919` on the orc unified
repo carries both commits; origin pushes remain Ryan-pending.

## The job

`benchmark/run_r12_champion_rerun.slurm.sh`: array 0–13 = (R1, R2) × j
(1/2/4/8/16/32/64) — the full per-class thread ladder (budget 0 uncapped
gets every point; positive budgets pruned by cap inside the tuner). One
exclusive 128-core zen3 500G node per task, qos=normal, 24 h. Stages per
task, sequential, each with its own STATUS_* verdict: fgstune → fgsprecond
(SWEEP_LADDER_1E6=1) → phase2 tuner (MEM_BUDGETS=0:16:32:64:128:500,
TUNE_MAX_SECONDS=14400) → phase2 (MEM_BUDGETS=16:32:64:128:500, full
CONFIGS). K_REPS=3. Run dirs
`r12-champion-<rung>-j<j>-<arrayjobid>/` under
`COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`.

Local pre-submit smoke (macOS, ≤4 threads, R1 @ j2, dagteam f64), full
4-stage chain PASSED: fgstune 16 candidates → τ=1e-6 winner
p8/MAC0.3/leaf150/inner10, verification PASS (bc 2.2e-8); fgsprecond ladder
PASS (niter=1, bc 2.5e-7 MEETS 1e-6); tuner budgets 0+16 certified winners;
phase2 full CONFIGS table incl. additivity check, fgs/fgmres_fgs rows carry
`sweep_order=dagteam;dagteam_precision=f64` in their config strings.

## Submission

Submitted 2026-09-19: **job 13778533** (array _0–_13 → (R1,R2) ×
j=1/2/4/8/16/32/64), after a clean `sbatch --test-only` (13778532, start
estimate 22:16 on m12).

**BLAS-threads ruling (Ryan question, same day):** the pipeline runs
uniformly at BLAS = j within each task (staircase and phase-2 share one BLAS
environment — no cross-BLAS tolerance carry). Measured sensitivity check
(local, R1 @ j4, dagteam f64, champion-style knobs): BLAS 1/2/4
indistinguishable (0.5505/0.5518/0.5519 s min over 5 reps, same iterations) —
the sweep's BLAS calls are per-leaf gemvs on ~leaf-sized blocks, below any
BLAS threading threshold, so the setting is inert for the FGS family while
BLAS = j keeps backslash_ldiv rows honest. Note for the record: the v23 R4
champion was certified and timed at BLAS 1 by harness convention and BLAS > 1
was never swept there; expected null for the same reason (leaf=100 blocks),
zen3 confirmation would be a cheap qos=test AB-trials rider if ever needed.
Tasks _0–_4 started before the estimate; _5–_13 were briefly user-held during
this check and released unchanged.

## Resume rerun 2026-09-21 (24 h-wall timeouts)

Six 13778533 arms hit the 24 h wall mid-p2tune (R1 j1/j2, R2 j1/j2/j4/j8;
fgstune+fgsprecond ok everywhere). The other 8 arms completed. Fix (source
tag `campaign/p021-resume-source-20260921` = `946cec2`; this worktree's exec
= clean pick, tag `campaign/p021-r12-champion-exec-20260921` = `e2b1410`,
pins.toml updated): wall 24→48 h and warm-start resume —
`RESUME_FROM_JOB_ID=13778533` reuses the old run dirs, skips the certified
fgstune/fgsprecond staircases, re-enters p2tune where row-level resume
(rung, budget, julia_threads) skips landed budgets, then runs p2 (stage-ok
skip retained there: p2 appends to phase2.csv). Old logs/pins per run dir in
`logs.before.<new job id>/`. Resubmitted as array 0,1,7,8,9,10 with the same
env + `RESUME_FROM_JOB_ID=13778533` — **job 13829232**, submitted 2026-09-21
after a clean `sbatch --test-only` (13829230, immediate start on m12).
