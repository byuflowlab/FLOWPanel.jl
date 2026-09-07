# Phase 23 — context reset: Phase 3 LAUNCHED (R1–R3), Phase 2 tail running

Written 2026-09-07 (evening). Supersedes `phase_22_phase3_launch_prompt.md` as
the entry point. Binding rules: `decision_rules.md` tail (2026-09-05+); dated
record: `ledger.md`. Ryan approved and I executed the Phase 3 launch on
2026-09-07 — this file records what is live and what remains.

## Pin of record for Phase 3: gen2d

- Tag **`021-lg-gen2d`** (annotated) at FLOWPanel `00390f7`, pushed to the
  cluster repo `~/flowpanel-021/FLOWPanel.jl`. Contains the phase3 knob
  sourcing (`phase2_knobs` + mandatory `KNOBS_BUDGET`), the `p3_warmstart.sh`
  submit guard, and the merge replay-dedup.
- Worktree **`/home/rander39/wt021/FLOWPanel.jl-d`**, branch `021-lg-gen2d-wt`,
  HEAD `e818c63` = tag + data-symlink commit (`data ->
  ~/projects/FLOWPanel.jl/data`). Built with `scripts/prep_campaign_worktree.sh`.
- Manifest copied from worktree-c (Project.toml identical between pins);
  dev-paths verified: `path="."`, `/home/rander39/wt021/FLOWVPM.jl` (`a627dd9`
  = `021-lg-gen2`), `/home/rander39/wt021/FastMultipole` (`0ce3ba6` =
  `021-lg-gen2`). Deps unchanged from gen2c.
- CSVs the driver reads, copied into worktree-d and verified:
  - phase1 FGS tables R1–R3 (tune/fgstune_selected/fgstune_staircase/
    fgsprecond); R1 md5s match the phase_22 record. Gaussian provenance —
    disclose in any published FGS row.
  - phase2 `tune_phase2.csv` R1/R2/R3, each with a **bc_certified budget-0
    row**: R1 17/0.65/6, R2 17/0.65/6 (both from the seed-override relaunches
    13603412/13, which COMPLETED and certified 2026-09-07), R3 16/0.65/6.

## Phase 3 fleet — LIVE (submitted 2026-09-07 from worktree-d)

All jobs `benchmark/slurm/p3_warmstart.sh`, `KNOBS_BUDGET=0`, hardware pinned
in-header (m12/zen3/exclusive/500G), `--time=24:00:00`. Per rung: checkpoint
(CONFIG=backslash, WARMSTARTS=cold, SAVE_VTK=true, RUN_NAME=p3_checkpoint_<R>)
then 6 arms **serialized with `afterany`** (all arms of a rung append to the
same `results/phase3/multi/<rung>/unsteady.csv` — never let them run
concurrently). Arms: RESTART_STEP=-1, RESTART_NAME=p3_checkpoint_<R>,
N_STEPS=72; backslash arm cold-only (null control), the 5 iterative arms run
cold:prev:extrap internally (WARMSTART_ORDER default 1, SKIP_STEPS default 3).

| rung | ckpt | arms (backslash, k_gmres, k_jacobi, k_ilu, fgmres_fgs, fgs) |
| --- | --- | --- |
| R1 | 13603453 | 13603454–59 |
| R2 | 13603519 | 13603520–25 |
| R3 | 13603460 | 13603461–66 |

Sharp edges: unsteady.csv APPENDS (rerun of a guess type duplicates rows);
walltimes are estimates — if an R3 checkpoint or arm hits 24 h, resubmit on
physics2 `--qos=standby` with longer `--time` (the one permitted override).
Judge by outputs, never sacct.

## Phase 2 tail (worktree-c — do NOT edit while its jobs run)

| jobs | what | state 2026-09-07 evening |
| --- | --- | --- |
| 13593015–17 | tune R1–R3 | COMPLETE |
| 13603412/13 | R1/R2 b0 seed-override arms | COMPLETE, certified |
| 13593018–21 | tune R4–R7 | RUNNING (~12 h in; 3/3/5/7 d walltimes) |

- b0 ladder status: R1 ✓ R2 ✓ R3 ✓ R5 ✓ (certified 09-05) R6 ✓ (certified
  09-07 with the OLD seed — the feared gate failure did not occur). Open: R4
  and R7.
- **R4 b0 arm — pre-approved by Ryan**: when 13593018 finishes (one writer per
  rung dir), submit from **worktree-c**:
  `sbatch --job-name=p2lg-tune-R4-b0 --time=8:00:00
  --export=ALL,RUNG=R4,MEM_BUDGETS=0,TUNE_SEED_B0=17:0.65:6
  benchmark/slurm/p2_tune.sh`
- **R7**: check its b0 row / `.err` for "budget 0.0 GiB: tuning FAILED" when
  13593021 lands; if failed, same seed-override treatment (ask Ryan first —
  only R4 is pre-approved).
- When R4–R7 land: full harvest via `benchmark/p021_merge.jl` from a checkout
  containing `00390f7` (worktree-d qualifies); phase_17 DUPLICATE trap stays
  active (last row per resume key wins).

## Disk / storage

**286.6 G / 400 G after the 2026-09-07 evening cycle** (was 412.5 G; per-subdir
du — root du times out, df is a <1%-accurate proxy). hpc-storage cycle with
two Ryan approvals: (a) keep-288 sweep of 7 CLOSED `ARCHIVED-STALE` p018
`_3r` runs (122.7 G gross), then the `_s1p5` run (83.1 G gross; job 13592732
COMPLETED 0:0 before the sweep reached it, so it swept CLOSED — its latest
kept step 2159 IS the run's final step, full 4-path warmstart set verified).
Total 201.0 G gross / 125.8 G net (writers added ~75 G concurrently). Restart
integrity verified on all 8 (288 steps each; monitors/CSVs/.pvd untouched).
Background: archiving was exhausted (zero archivable runs); the stale flags
were proven resumed-chain churn, not interrupted deletes — `--resume-delete`
would fail verification, keep-288 sweep is the right tool; swept runs keep
re-flagging ARCHIVED-STALE by construction (288 > the archiver's 5-step
window) — expected, not an alarm.

**Standing pressure: burn was ~18–30 G/h** during Phase 3 spin-up, queue
winding down (8R/12P at 23:00Z). ~113 G headroom ⇒ another cycle likely
2026-09-08. **Next lever is ARCHIVING, not sweeping**: no sweepable surplus
remains; `_s1p5`'s residual 17.5 G (and other finishing `_3r` arms) become
archivable after 24 h quiet, ~2026-09-08 15:21 local. Remaining stale runs
yield ~3 MB each — ignore.

Hazards: `*_exp_nt` runs misclassified CLOSED by the sweeper while jobs
13603468/9 write them (name mangling — always exclude; consider protect-list
entries if those chains run long); archiver enumerates only 3 checkouts,
skipping wt021/wt018/wt026 worktrees (~6 G, possible alias-dedup bug).
Phase 3 checkpoints RECENT-HOT/protected (R1 step 434, R3 step 199 at sweep
time). Mount has 1.6 T free — cap is policy, not crash risk. Cluster logs:
`~/archiver_dryrun_20260907.log`, `~/sweep_{dryrun,apply}_20260907.log`,
`~/sw_p018_csarc_*.log`, `~/du_watchdog_20260907{,_after}.log`. Full ledger
line drafted in the agent report (this session's transcript), not yet written
to the archiver ledger — the agent's report text has it verbatim.

## Cluster mechanics

`ssh orc` needs a live ControlMaster socket (2FA otherwise). Slurm binaries via
`/apps/slurm/latest/bin/`. `du` on `/home/rander39` root times out — measure
per-subdirectory. No sbatch/scancel beyond the pre-approvals above without
Ryan asking in the moment.

## Open with Ryan

- Ledger entry (not yet written, offered): b0 root cause + seed-override
  relaunches 13603412/13 + FGS CSV copy-forward + gen2d pin + Phase 3 launch
  record (this file has all the facts).
- Notebook entry for the campaign (offer, don't write).
- R7 b0 treatment if its old-seed budget-0 fails.
- Phase 3 harvest plan once arms complete (headline: niter_first cold-vs-warm,
  break-even step count; metrics rulings in decision_rules.md).
