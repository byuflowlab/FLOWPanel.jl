# Reset prompt — p033 R4 thread-scaling campaign (2026-10-02)

You are continuing a partially-implemented, Ryan-approved campaign. Read
`BRAINSTORM/033_threadscale_r4_20261002.md` FIRST — it is the campaign doc
(scope, fixed knob table, submission commands) and this prompt only carries
the session state on top of it.

## Approved scope (do not re-litigate)

R4 ONLY, j ∈ {1,8,16,32,64}, latest FGS (threaded setup, FastMultipole
6456c221) vs krylov_ilu_nfcache, FIXED transplanted knobs (NO re-tuning —
Ryan 2026-10-02 ruling, reversing an earlier full-re-tune ruling the same
day). Job A = cold (setup + isolated cold solves), Job B = warm rotor hover
per the 021 winB protocol under 018-ported wake physics (WAKE_ENV_018=1).
Ported physics ⇒ old R4 checkpoints invalid ⇒ fresh family ckpt legs at j64,
then winB legs at every j. Final deliverable: markdown comparison table
across thread counts (Ryan: "Show me the comparison across thread counts").
Approved plan file: ~/.claude/plans/launch-a-new-job-effervescent-squirrel.md.

## State: ALL code changes are written but UNCOMMITTED on `fastmultipole`

(worktree also has pre-existing committed HEAD db18b94 + untracked 018 files
`BRAINSTORM/018_dji9443_hover_convergence_campaign/ntladder_*_20261001.md` —
leave those alone, item 018's.)

New/edited files (all in the live checkout, ready to commit):
1. `benchmark/p033_wake_env_018.jl` — NEW. 018 wake defaults via _setdefault!
   gated on WAKE_ENV_018=1; Das-table sha256 preflight
   (640ba059cf57d6456ded0bf65721326160d90987fa449d1b6cc276d69fe755bf); banner.
   FIX APPLIED after first smoke failure: PARTICLE_SHEDDING=sigma_overlap +
   DAS_ETA_KINEMATIC=1.0 + TRUNCATION_DEPTH_R=4 (SIGMA_CHORD_FRACTION
   hard-requires sigma_overlap; these are p018 dispatcher-level settings).
2. `benchmark/fgs_r4_warmstart_ab.jl` — include of (1) right after
   _setdefault! definition.
3. `benchmark/rotor_hover_solver_unsteady.jl` — FGSSolver ctor gains
   `threaded_setup=fgs_threaded_setup` (env FGS_THREADED_SETUP, default "1");
   recorded in the knobs provenance string. No CSV schema change.
4. `benchmark/fgs_setup_ab.jl` — CERT_SOLVE block now persists cold-solve
   wall times as CSV rows (phase=cold_solve, AB_SOLVE_K passes, default 3).
5. Launchers (new): `benchmark/run_p033_ts_cold_r4.slurm.sh` (array 0-4 over
   j; FGS fgs_setup_ab + ILU phase2.jl at PRE-SEEDED tune_phase2.csv rows the
   launcher writes — budget 0 P15/0.55/21, budget 500 P12/0.55/48),
   `benchmark/run_p033_ts_ckpt_r4.slurm.sh` (array 0-1 family, j64, 108
   steps, VTK, RUN_NAME=p033ts_R4_ckpt_<family>),
   `benchmark/run_p033_ts_winb_r4.slurm.sh` (array 0-4 over j, fgs_proj2 +
   ilu_nfcache_proj2 sequentially, RESTART_NAME/PATH=data/p033ts_R4_ckpt_*,
   submit with --dependency=afterok:<ckpt job>).
6. `benchmark/p033_ts_harvest.jl` — NEW markdown tabulator (cold + warm
   tables; reads p033ts-cold-j*/p033ts-winb-j* run dirs under ARGS[1]).
7. `BRAINSTORM/033_threadscale_r4_20261002.md` — campaign doc (smoke section
   has one <pending> placeholder to fill).

## Smoke status

- fgs_setup_ab smoke (R1 j4): PASS — all certs bitwise (incl. dagteam
  Lmat/Umat), niter 15=15, x rel=0.0, ctor 10.11→2.37 s, cold_solve rows
  land (min 0.279/0.280 s). Evidence:
  /private/tmp/claude-502/-Users-ryan-Dropbox-research-projects-FLOWPanel-jl/30c99ba6-f73b-46bd-bb9a-a86d43e18d1f/scratchpad/smoke_ab/
- warm-start smoke rerun (R1, NT=4, stages fgs_cold:ckpt ilu_nfcache_cold:ckpt
  fgs_proj2:winB ilu_nfcache_proj2:winB, WAKE_ENV_018=1): was RUNNING at
  reset in .../scratchpad/smoke_ws/ (fgs ckpt marching; the first attempt
  failed on the PARTICLE_SHEDDING gate, now fixed). The background process
  may have died with the session — if smoke_ws lacks 4 COMPLETED_* files /
  "SMOKE PASS", RERUN:
    cd ~/Dropbox/research/projects/FLOWPanel.jl && \
    WAKE_ENV_018=1 SMOKE_ROOT=<scratch>/smoke_ws2 \
    STAGES="fgs_cold:ckpt ilu_nfcache_cold:ckpt fgs_proj2:winB ilu_nfcache_proj2:winB" \
    bash benchmark/run_r4_fgs_warmstart_smoke.sh
  Judge: SMOKE PASS + all four STATUS_*=ok + winB legs' unsteady.csv rows all
  solved=true and nsolves==1. Also eyeball the "WAKE PHYSICS (018-ported)"
  banner in a leg log. (≤4 threads locally — hard policy.)

## Remaining steps (in order)

1. Warm smoke to PASS (above). Fill the <pending> smoke line in the campaign
   doc. No src/ changes were made, so the full unit suite is NOT required.
2. Commit everything (one commit on `fastmultipole`), message along
   "033 p033ts: R4 thread-scaling campaign (fixed knobs): 018 wake env port,
   threaded-setup opt-in, cold-solve CSV rows, launchers + harvest".
3. Annotated tag `campaign/p033-threadscale-r4-20261002` in all THREE repos:
   FLOWPanel (the new commit), FastMultipole @ 6456c221 (= current
   flowpanel-20260817 HEAD), FLOWVPM @ eebb984 (= current flowpanel HEAD).
4. Push branch fastmultipole + the tag to the ORC remotes only (GitHub origin
   pushes are Ryan-gated). BLOCKER at reset: `ssh orc` had NO live
   ControlMaster socket (Permission denied) — ask Ryan to open one
   (`! ssh orc` in-session) before any cluster step.
5. On orc: `bash scripts/prep_campaign_worktree.sh campaign/p033-threadscale-r4-20261002
   ~/campaigns/p033-threadscale-r4-20261002`; pinned FastMultipole/FLOWVPM
   worktrees + Manifest dev-paths; pins.toml; verify data symlink, Das table
   (md5 08375291ed2b542ea946d09730e0b629), R4 mesh. Data root:
   ~/projects/FLOWPanel.jl/data/p033_threadscale_r4_20261002/.
6. Run the slurm-availability skill (zen3, 128 cpu, 500G, qos_normal; 48 h
   winB) BEFORE submitting. Precompile-gate job first; submission command
   block is in the campaign doc §Submission. VAR=x sbatch form ONLY (never
   comma --export). julia 1.11.7 (spack path first — 1.12 default breaks).
7. While the ckpt job writes VTK: launch the hpc-storage subagent (policy).
   Monitor via hpc-monitor; judge by STATUS_*/COMPLETED + CSV solved, never
   sacct.
8. Harvest: rsync CSVs local, `julia benchmark/p033_ts_harvest.jl <root>`,
   hand Ryan the comparison tables (label: knobs fixed-across-j transplant;
   winB first order+1 steps effectively cold).

## Gotchas carried

- Session cwd can drift to the FLOWVPM.jl additional working dir — absolute
  paths.
- Local jobs ≤4 threads (global policy).
- The winB launcher preflights data/p033ts_R4_ckpt_{fgs,ilu}; ckpt VTK goes
  through the worktree's data/ symlink (save_path = data/<RUN_NAME>).
- fgs_setup_ab transplant warning at R1 smoke is expected (champion TOML is
  R4); at R4 production it's silent.
- Old 021 R4 checkpoints must NOT be reused (different physics) — names are
  campaign-qualified precisely for this.
- SFS_OFF=true and SIGMA_FLOOR_R=0 are deliberate defaults (flagged to Ryan
  in the approved plan); NT=36 everywhere (RELAX_RLXF=0.3 is NT-dependent).
