# 026 reset prompt — wave-2 planning after de-risk full PASS (2026-09-15)

You are picking up BRAINSTORM 026 (resolution-preserving particle splitting)
in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (+ siblings
FLOWVPM.jl, FastMultipole). Read `CLAUDE.md` and the policies it names
(HPC.md before cluster work; `ssh orc` needs a live ControlMaster socket —
ask Ryan to run `! ssh orc echo ok` if 2FA blocks you). Local git branch is
**fastmultipole** (not flowpanel).

## State as of 2026-09-15 evening (all Ryan-authorized)

1. **De-risk trio harvested → acceptance fully PASS.** Record =
   `derisk_harvest_20260915.md` (this dir; includes the extension section).
   Highlights: merge A/B at campaign σ shows NO 2.7× retention drop (B/A
   final 1.12 after late crossover; no σ-pump; merge cost free; B fires
   4.7× more viscous splits late); floor 0.1 clean → **Ryan ruled floor
   0.1 kept for the remaining arms**; telemetry gap list (no per-event
   merge/floor logging, no σ mean/max census, no FMM-adequacy metric).

2. **NREVS root cause found & fixed**: the dispatcher
   `examples/run_p018_screen_hpc.slurm.sh` exports NREVS=8 UNCONDITIONALLY
   after sbatch env lands ("case arms override unconditionally" pattern) —
   submission-time NREVS is always clobbered; run length MUST be set inside
   the case arm. The de-risk trio therefore all ran 8+1 revs. Fixed for
   cap030 (`export NREVS=17` in its case arm) in commit `1b59af5`, tag
   **`campaign/p026-derisk-ext-20260915`**, pushed; campaign FLOWPanel
   worktree at `orc:~/campaigns/p026-derisk-20260914/FLOWPanel.jl` is
   checked out on it (FLOWVPM bf88806 / FastMultipole ac7230a6 unchanged).

3. **cap030 extension (job 13694747) COMPLETED: NO CLIFF.** Chained
   restart steps 1296–2591 (rev 18): smooth N-scaling 10→18.5 s/step,
   nothing in the 2200–2248 window, CT cycle-mean 0.07484 ±0.34%,
   Phase-2e CONVERGED=false only on within-rev p-p (0.0326 > 0.02 tol;
   per-rev spread passes). Gotchas learned: monitor04 wake-health CSV does
   NOT append across restarts (CT CSVs do); viscous splits (f_visc=0.587)
   have NEVER fired in the cap030 regime through step 2591 (exp arms do
   fire them late) — keep watching.

4. **Hygiene done**: de-risk + smoke run dirs all in shared root
   `orc:~/projects/FLOWPanel.jl/data/` (campaign worktree keeps symlinks);
   wt026gpu silo REMOVED (worktrees verified clean at pushed pins first);
   scr_p026cpuv_split + scr_p026gpuv_splitmerge archived & verified;
   **scr_p026gpuv_split archive PENDING** (archiver RECENT-HOT hold —
   retry via hpc-storage once quiet ≥24h; its wt026gpu slurm logs are
   inside it at `wt026gpu_logs/`). Ledger = `ledger.md` (this dir).
   `data/p026_restart_gpu40_s950` left for the archiver.

## YOUR TASK — wave-2 plan (present to Ryan BEFORE submitting; launch is gated)

Ryan held wave-2 for the NREVS fix; that fix pattern is now known. Build
the launch plan for the remaining 11 arms (per `campaign_prep_20260912.md`
§2 and §8.4: cap018 + the s020v family minus the two already run — use
brainstorm-scout to pull the arm list/lengths, don't read the file inline):

1. For EACH arm, set the intended run length as `export NREVS=<n>` INSIDE
   its case arm in `examples/run_p018_screen_hpc.slurm.sh` (remember
   n_steps = NT*(NREVS+SPINUP_REVS); screen default SPINUP_REVS=1, NT=36
   unless the case sets NT=144). Audit every wave-2 case def for the same
   missing-NREVS gap.
2. Commit (Ryan-gated — bundle with any pending docs), new tag
   `campaign/p026-wave2-<date>` on FLOWPanel only if FLOWVPM/FMM are
   unchanged; fast-forward the campaign worktree (reuse
   `~/campaigns/p026-derisk-20260914/`, env Manifest already points there).
3. Optional cheap hardening Ryan may want first (from the gap list):
   per-event merge logging + σ census + floor-clamp counter — that would
   touch FLOWVPM → new FLOWVPM pin too. Present as an option with cost.
4. Present plan + cost table (arms × revs × ~6–18 s/step by regime; H200
   m13h; the de-risk trio ran fine there) → Ryan approves → submit with
   the §Arms common env from `p026_derisk_20260914_provenance.md`
   (SIGMA_FLOOR_FRAC=0.1 WAKE_SPLIT_VISCOUS=true
   WAKE_SPLIT_FRAC_VISCOUS=0.587 WAKE_SPLIT_FRAC_COMPRESS=0.73
   WAKE_SPLIT_FRAC_ELONGATE=0.3), wrapper
   `~/projects/FLOWPanel.jl/examples/run_p018_screen_gpu052.slurm.sh h200
   <case>` with P018_REPO_OVERRIDE/P018_PROJECT_OVERRIDE at the campaign
   worktree/env, cwd = campaign worktree. MOVE each run dir to the shared
   root + symlink at harvest (never pre-place symlinks — dispatcher wipes
   data/$RUN_NAME on cold runs).

## Owed / parked

- **Notebook**: whole 026 arc unlogged (Ryan "not yet" ×2, 09-14/15) —
  offer an entry at the next milestone (σ-check verdict, GPU splitting,
  de-risk wave + extension, merge A/B). Also owed from earlier arcs:
  expguard three-arm result, Cd transient, rlxf derivation, Ladder C
  forensics.
- scr_p026gpuv_split archive retry (above).
- f_visc=0.587 never firing in grow-regime arms — flag in wave-2 harvest.
- Uncommitted locally: `data/rotor_hover_pressure_comparison.metadata.toml`
  (stays uncommitted, perpetually rewritten), 021-arc BRAINSTORM files,
  `ledger.md` + this file + harvest-doc/provenance appends (bundle into
  next approved commit).

## Ground rules (carry-over)

Local ≤4 threads; Slurm in non-login shells needs `bash -lc`; MOTD banners
contaminate ssh output — filter. sacct state is not evidence; judge runs by
outputs (.err before .out — though for this harness sim output is in .out).
Commits, campaign launches, notebook writes are Ryan-gated. Known unrelated
failure: `runtests_unit_warmstart.jl` first testset. An unrelated p021 job
(13694724, another session's) may be in the queue — leave it alone.
