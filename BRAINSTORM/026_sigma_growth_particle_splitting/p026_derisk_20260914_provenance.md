# 026 de-risk wave provenance — campaign/p026-derisk-20260914

Authorization: Ryan 2026-09-14 (three rulings in-session): 3-arm de-risk
first; cold-start (GPU); f_visc enabled at 0.587 (= 4^(1/3)−1, count-matched
tetra4 analog of f_comp = √3−1, picked from AskUserQuestion options);
SIGMA_FLOOR_FRAC=0.1 for the trio (0.25 fallback if the trio fails; keep 0.1
for the remaining 11 arms if it passes).

## Pins (annotated tag `campaign/p026-derisk-20260914`, pushed to origin)

| repo | branch | tag commit |
|---|---|---|
| FLOWVPM.jl | flowpanel | `bf88806` |
| FLOWPanel.jl | fastmultipole | `035f50b` |
| FastMultipole | flowpanel-20260817 | `ac7230a6` |

Worktrees (orc): `~/campaigns/p026-derisk-20260914/{FLOWPanel.jl,FLOWVPM.jl,FastMultipole}`
created from the tag; campaign env `~/campaigns/p026-derisk-20260914/env`
(copy of `~/projects/envs/x86_64` with the three dev-paths repointed at the
campaign worktrees). Worktrees carry no uncommitted state.

## GPU-splitting verification (precondition, 2026-09-14)

- Local: FLOWVPM `runtests_resolution_split.jl` 1055/1055 (incl. t12
  CPU-vs-broadcast parity), FLOWPanel wake/replay suites clean.
- CUDA smoke `scr_p026gpuv_split` (job 13689273, m13h H200, warm-start
  gpu40 s950, 40 steps): **gate PASS** (gpu_gemv=40, cpu_gemv=0, nan=0,
  rc=0); all three mechanisms fired on GPU (viscous/compress/elongate);
  ~7.6 s/step vs CPU twin ~176 s/step (≈23×).
- `scr_p026gpuv_splitmerge` (13689274, MERGE_OVERLAP=3.5): code path ran
  (merge gate active, splits firing); died step ~974 on the euler_exp
  substep-budget guard (dt·|L| ≈ 4.8e3) after an elongate-split storm —
  the ignition continuation blowing up harder under the overlap merge gate.
  Treated as a physics observation (the backstop worked), not a code
  failure; it foreshadows the campaign merge A/B.
- CPU twin `scr_p026cpuv_split` (13689275): COMPLETED all 40 steps.
  **Split-rate parity confirmed**: 34 split lines on both backends; final
  steps CPU elongate 161-199 events/step (480-593 children) vs GPU 172-198
  (511-588 children) — RNG-level agreement. ~176 s/step CPU vs ~7.6 s/step
  GPU (~23x).

## Arms (submit from the campaign FLOWPanel worktree)

Common submission env: `SIGMA_FLOOR_FRAC=0.1 WAKE_SPLIT_VISCOUS=true
WAKE_SPLIT_FRAC_VISCOUS=0.587 WAKE_SPLIT_FRAC_COMPRESS=0.73
WAKE_SPLIT_FRAC_ELONGATE=0.3`; cold start (no RESTART_*); GPU via
`run_p018_screen_gpu052.slurm.sh h200 <case>` with
`P018_REPO_OVERRIDE`/`P018_PROJECT_OVERRIDE` at the campaign worktree/env.

| arm | case | extra env | length |
|---|---|---|---|
| grow cap030 | `scr_p026sp_nt144_cap030` | (linegauss default; `WAKE_SPLIT_SIGMA_MAX=0.030` is in the case def as emission clamp) | `NREVS=20` (through the ~step-2250 cliff region) |
| shrink split | `scr_p026s9_exp_split` | `FLOWPANEL_FILAMENT_REG=vatistas` (exp bracket predates linegauss default) | `NREVS=12` (ignition ~step 210 + margin) |
| merge A/B twin | `scr_p026s9_exp_split_mo35` | `FLOWPANEL_FILAMENT_REG=vatistas` (case def carries `MERGE_OVERLAP=3.5`) | `NREVS=12` |

Acceptance (de-risk): all arms run to completion or die on a *guard* with
interpretable telemetry; cap030 shows no cliff and 6–8 s/step-class cost;
split/skip counters sane; exp pair separates the merge policies; floor 0.1
does not destabilize the healthy phase.

## Data

The dispatcher wipes/moves `data/$RUN_NAME` before a cold run, which kills a
pre-placed symlink, so: runs write into the campaign worktree's `data/<case>/`
during execution (Das arc table symlinked per-file into the worktree), and
each run dir is MOVED to the shared root `~/projects/FLOWPanel.jl/data/`
with a symlink left behind immediately at harvest. Storage/archiver agents
attribute by realpath as usual.

## Submissions (2026-09-14 evening)

| job | arm | pool |
|---|---|---|
| 13691080 | `scr_p026sp_nt144_cap030` (NREVS=20) | m13h H200, 24 h |
| 13691081 | `scr_p026s9_exp_split` (NREVS=12, vatistas) | m13h H200, 8 h |
| 13691082 | `scr_p026s9_exp_split_mo35` (NREVS=12, vatistas) | m13h H200, 8 h |

Submitted from the campaign worktree with the common split/floor env
(§Arms). Banners to be verified at job start (ops rule).

## Extension: cap030 chained restart (2026-09-15, Ryan-approved)

Harvest (see `derisk_harvest_20260915.md`) found all three arms ran the
dispatcher default NREVS=8 (+1 spinup = 9 revs): the dispatcher exports
NREVS unconditionally after sbatch env lands, and the de-risk case defs
carried no NREVS override — the submitted NREVS=20/12 were clobbered.
cap030 therefore stopped at rev 9, short of the historical FMM-adequacy
cliff at rev ~15.3 (steps 2200–2248 @ NT144). Acceptance otherwise PASS →
floor 0.1 kept; wave-2 HELD pending this fix (Ryan 2026-09-15).

Fix: `export NREVS=17` added to the `scr_p026sp_nt144_cap030` case arm
(total 18 revs incl. spinup = 2592 steps; cliff bracket + ~2.4 rev margin)
plus a clobber-warning comment at the dispatcher default. New tag on
FLOWPanel only: **`campaign/p026-derisk-ext-20260915`** (launcher-only
change; FLOWVPM stays `bf88806`, FastMultipole stays `ac7230a6` on the
original tag). Campaign FLOWPanel worktree fast-forwarded to the new tag;
env Manifest dev-paths unchanged (same worktree path).

Run: chained restart, `RESTART_STEP=1295` from the retained state in
`~/projects/FLOWPanel.jl/data/scr_p026sp_nt144_cap030/` (reached through
the worktree's data symlink; restart mode preserves the run dir and appends
to the same VTK series), same common split/floor env as §Arms, m13h H200.
Expected +1297 steps at 10–14 s/step ≈ 4–6 h.
