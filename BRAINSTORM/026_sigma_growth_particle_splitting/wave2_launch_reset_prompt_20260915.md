# 026 reset prompt — wave-2 LAUNCH after slate rulings (2026-09-15 late)

You are picking up BRAINSTORM 026 (resolution-preserving particle splitting)
in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (+ siblings
FLOWVPM.jl, FastMultipole). Read `CLAUDE.md` and the policies it names
(HPC.md before cluster work; `ssh orc` needs a live ControlMaster socket —
ask Ryan to run `! ssh orc echo ok` if 2FA blocks you). Local git branch is
**fastmultipole**. Predecessor context: `wave2_reset_prompt_20260915.md`
(same dir) — state as of the de-risk PASS; this file supersedes its task.

## Session rulings (Ryan, 2026-09-15 evening — all on record)

1. **Linegauss only; NO vatistas** in wave-2. This collapses the §9 matrix:
   the non-lg `ctrl_*`/`exp_*` arms existed solely as the vatistas brackets
   and are DROPPED (cold-start linegauss non-lg ≡ lg twins byte-for-byte).
2. **Keep the three-way mechanism distinction**: `_floor` = floor ONLY
   (no split knobs of any kind), `_split` = splitting ONLY (viscous +
   stretch, NO floor), `_fs` = floor + splitting. The earlier
   floor-on-every-arm reading is DEAD — SIGMA_FLOOR_FRAC goes only on
   `_floor`/`_fs`/cap018.
3. **Vatistas-vs-linegauss comparison**: two runs at the exact de-risk
   parameters under linegauss. The de-risk A-arm config (floor 0.1 +
   stretch + viscous) under lg IS the matrix `explg_fs` arm; a new twin
   case def `scr_p026s9_explg_fs_mo35` (adds MERGE_OVERLAP=3.5) was added
   to the dispatcher as the B arm. Compare against vatistas de-risk pair
   `scr_p026s9_exp_split`/`_mo35` (data in shared root).
4. Floor value stays 0.1 wherever a floor is armed (de-risk ruling).

## Local uncommitted state (edits DONE this session, commit is Ryan-gated)

`examples/run_p018_screen_hpc.slurm.sh`:
- `export NREVS=12` added inside ALL 13 `scr_p026s9_*` case arms (incl.
  already-run exp_split/_mo35, hygiene) — fixes the silent 8+1 clobber
  (dispatcher line ~71 exports NREVS=8 unconditionally by design).
- `export NREVS=17` added to `scr_p026sp_nt144_cap018` (2592 steps at
  NT=144 incl. 1 spinup rev = cliff bracket + ~2.4 rev margin, mirrors
  cap030's revised length; original 20 was pre-revision).
- NEW case def `scr_p026s9_explg_fs_mo35` (clone of explg_fs +
  MERGE_OVERLAP=3.5 + NREVS=12), inserted after explg_fs with comment.
- bash -n clean.

Also uncommitted: `ledger.md`, both wave2 reset prompts, harvest/provenance
appends, 021-arc BRAINSTORM files (bundle into the commit);
`data/rotor_hover_pressure_comparison/rotor_hover_pressure_comparison.metadata.toml`
stays PERPETUALLY uncommitted.

## THE SLATE — 8 arms, all cold-start, all linegauss, m13h H200

| # | case | mechanism | submission env beyond case def | walltime |
|---|---|---|---|---|
| 1 | `scr_p026s9_ctrllg_floor` | floor only | `SIGMA_FLOOR_FRAC=0.1` | 8 h |
| 2 | `scr_p026s9_ctrllg_split` | split only | `WAKE_SPLIT_VISCOUS=true WAKE_SPLIT_FRAC_VISCOUS=0.587 WAKE_SPLIT_FRAC_COMPRESS=0.73 WAKE_SPLIT_FRAC_ELONGATE=0.3` | 8 h |
| 3 | `scr_p026s9_ctrllg_fs` | floor+split | floor env (#1) + split env (#2) | 8 h |
| 4 | `scr_p026s9_explg_floor` | floor only | as #1 | 8 h |
| 5 | `scr_p026s9_explg_split` | split only | as #2 | 8 h |
| 6 | `scr_p026s9_explg_fs` | floor+split = **lg merge-A arm** | as #3 (exact de-risk env under lg) | 8 h |
| 7 | `scr_p026s9_explg_fs_mo35` | fs + overlap merge gate = **lg merge-B arm** | as #3 (MERGE_OVERLAP=3.5 in case def) | 8 h |
| 8 | `scr_p026sp_nt144_cap018` | grow-side cap (σ_max=0.018 clamp in case def) | as #3 (matches cap030 de-risk env; compress 0.73 kept — ruled + passed de-risk) | 24 h |

Notes: `_split` arms (no floor) leave WAKE_SPLIT_SIGMA_MIN defaulting to
NaN (no lower emission clamp) — intentional consequence of "no floor",
driver prints it as "off". The driver HARD-ERRORS on any split fraction
without its mechanism flag and vice versa (rotor_hover_pressure_comparison.jl
~:775-785) — floor arms must carry NO split knobs. WAKE_HEALTH monitors are
in the case defs already. Case defs 2–7 carry WAKE_SPLIT_STRETCH=true and
linegauss; check the banner of EVERY job after submission (mandatory).

Cost: s9 arms ~5.7 s/step early → ~1–1.5 h each (468 steps); cap018 cold
2592 steps, cap030 measured 7.5→18.5 s/step (~8 h total) — tighter cap →
budget 9–12 h. Total ≈ 17–23 GPU-h, wall ~1 day if concurrent.

## OPEN GATES (get Ryan's ruling BEFORE acting)

- (iii) optional telemetry hardening BEFORE launch: per-event merge/floor
  logging + σ census (touches FLOWVPM → new FLOWVPM pin + two-repo tags,
  ~half-day + smoke, delays ~1 day). Unanswered as of reset.
- (iv) approve the bundled commit + annotated tag
  `campaign/p026-wave2-20260915` (FLOWPanel only if FLOWVPM bf88806 /
  FastMultipole ac7230a6 unchanged; tag convention per global CLAUDE.md).
- (v) approve submission of the 8 arms.

## Execution SOP once gated (per provenance p026_derisk_20260914_provenance.md)

1. Commit bundle → tag → push tag+branch.
2. Fast-forward campaign worktree `orc:~/campaigns/p026-derisk-20260914/
   FLOWPanel.jl` to the new tag (env Manifest already points there;
   FLOWVPM/FMM worktrees unchanged).
3. Submit from cwd = campaign worktree: `sbatch --job-name=fp-026-<short>
   --export=ALL,<env from slate table> ~/projects/FLOWPanel.jl/examples/
   run_p018_screen_gpu052.slurm.sh h200 <case>` with
   P018_REPO_OVERRIDE/P018_PROJECT_OVERRIDE at the campaign worktree/env
   (copy the exact pattern from the provenance §Submissions).
4. Verify all 8 banners (nrevs must read 12/17, reg linegauss, floor/split
   knobs per slate; the NaN log-gate greps word NaN — split knobs print
   "off" when unarmed).
5. Record submissions in a new provenance section/appends + `ledger.md`.
6. At harvest: MOVE run dirs to shared root `orc:~/projects/FLOWPanel.jl/
   data/` + symlink (NEVER pre-place symlinks — dispatcher wipes
   data/$RUN_NAME on cold runs). Watch: f_visc=0.587 firing (never fired
   in grow regime; fired late in shrink); monitor04 wake-health CSV does
   NOT append across restarts; sacct state is not evidence.

## Owed / parked (carried)

- Notebook: whole 026 arc unlogged (Ryan "not yet" ×2) — offer at next
  milestone. Earlier arcs also owed: expguard three-arm, Cd transient,
  rlxf derivation, Ladder C forensics.
- scr_p026gpuv_split archive retry via hpc-storage once quiet ≥24 h
  (RECENT-HOT hold); `data/p026_restart_gpu40_s950` left for archiver.
- Vatistas-vs-linegauss merge verdict re-anchor = arms #6/#7 vs de-risk
  pair (flag in harvest).

## Ground rules (carry-over)

Local ≤4 threads; Slurm in non-login shells needs `bash -lc`; MOTD banners
contaminate ssh output — filter. Judge runs by outputs (.out for sim
output). Commits, campaign launches, notebook writes are Ryan-gated. Known
unrelated failure: `runtests_unit_warmstart.jl` first testset. Unrelated
p021 job 13694724 (another session's) may be queued — leave it alone.
