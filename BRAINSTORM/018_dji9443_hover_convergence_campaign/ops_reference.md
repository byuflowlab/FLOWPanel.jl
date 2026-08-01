# 018 Ops Reference (read once per session before touching the cluster)

## Cluster access and deployment

- Cluster checkout: `orc:/home/rander39/projects/FLOWPanel.jl` (hard-asserted
  by the launcher). **Ryan must open the ControlMaster socket first**
  (`! ssh orc -fN` in-session; BYU 2FA otherwise) — agents then multiplex
  `ssh`/`scp`/`sbatch`/`squeue` over it, wrapped `bash -lc "..."`.
- Deploy pattern = `scripts/deploy_phase02e.sh` (adapt, don't reinvent): `scp` a
  COPY_FILES list (launcher, driver, any new mesh) to the cluster repo path,
  then **md5-verify every `src/` file against local** — stale cluster `src/`
  silently produced a monitor bug in Phase 2c. Print remote branch/HEAD +
  `squeue` after staging. Nothing about deployment submits jobs.
- `mkdir -p logs/slurm` on the cluster before any submission (Slurm opens log
  paths before the script runs).

## Submission

- One case per job, from the cluster repo top level:
  `sbatch --job-name=fp-018-<tag> examples/run_dji9443_hover_ct_hpc.slurm.sh p018_<tag>`
- All 018 jobs: 64 cores / **64G** / `--time=24:00:00` default; NT=72 cases
  `--time=48:00:00`; NT=144 and σ-L2 segments per their phase files.
- **≤6 active study jobs** — check `squeue` before every submission.
- Final-settings cases (`p018_trunc6`, `p018_green`, `p018_nomerge`,
  `p018_final`) take Das\*/NT\*/σ\* via environment at submission, e.g.
  `sbatch --export=ALL,DAS_ETA_KINEMATIC=1.0,NT=36,OVERLAP=2.4,P_PER_STEP=11,MERGE_R_FACTOR=0.0053 ...`
  (the `p018_*` carrier honors env overrides for every knob).

## Monitoring and harvest

- Status: `ssh orc 'bash -lc "/home/rander39/projects/FLOWPanel.jl/scripts/p2e_status.sh brief"'`
  — parse ONLY the sentinels (`P2E_NJOBS=`, `P2E_JOB=`, `P2E_STEP=`,
  `P2E_CT=`, `P2E_DIVERGED=`). **`P2E_DIVERGED reason=nan_inf` is a false
  positive** (all-NaN `CT_kj` column); only `reason=magnitude` is real.
- Harvest per case: `scp` `data/<run>/<run>_CT_vs_rev.csv`,
  `<run>_CT_per_rev.csv`, `<run>_case_metadata.toml`, and
  `data/<run>/monitors/*_bound_circulation_1.csv` to the matching local path;
  append the settled numbers + provenance (job id, date) to `ledger.md`.
- Retention (ruling 10): launcher wipes `data/$RUN_NAME` on rerun; use
  `RHPC_KEEP_PREV=true` when the previous attempt failed uninspected. VTK
  stays ON during runs, but **after harvest delete a closed run's VTK, keeping
  CSVs/TOML/monitors** — e.g.
  `find data/<run> -name "*.vtu" -o -name "*.vtp" -o -name "*.pvd" | xargs rm`
  — and keep **≤3 runs' VTK at any time**. If a convergence question needs
  wake visualization, ask Ryan to examine the ParaView files before deleting.
  Log every VTK deletion in the phase file.
- **Divergence triage (Ryan-affirmed rule): OOM inside `merge_particles!`
  typically means the simulation blew up** — the diverged cloud wrecks the
  merge spatial hash; memory is the symptom, NOT a ceiling or a merging bug.
  Check the log tail for |CT| > 1 (caveat: stdout is block-buffered, so a
  SIGKILLed log truncates — a sane CT prefix is not evidence of boundedness)
  and MaxRSS vs ReqMem via `sacct`. Genuine memory cap (rare): CT sane AND
  MaxRSS ≈ ReqMem (e.g. `p2e_sigF_nofilt`, 32G request / 33.4 GB RSS). Only
  then is a memory bump justified; otherwise record the divergence and move
  to the phase's branch table.

## Restart chaining (extending a run)

The driver warm-starts when `RESTART_STEP >= 0`
(`examples/rotor_hover_pressure_comparison.jl` ~line 546). Recipe to extend
case `p018_X` by ΔN revs:

1. Determine the last *restartable* step S = the highest index in
   `data/p018_X/p018_X_body1.pvd` (equivalently the last `.vtu` snapshot).
   **Gotcha (hit 2026-07-31, job 13006768 FAILED in 3 min):** the CT CSV's last
   `step` is S+1 — CSV rows are 1-indexed per timestep while VTU snapshots are
   0-indexed — so "CSV last step" is NOT restartable. For a 719-step run:
   CSV ends at 720, restart from 719.
2. Submit the same case tag with:
   `--export=ALL,P018_RUN_NAME=p018_X_s2,P018_SETTLE_REVS=<12+ΔN>,RESTART_STEP=<S>,RESTART_NAME=p018_X,RESTART_PATH=data/p018_X`
   (plus the case's own knob overrides for final-settings cases). The segment
   writes to `data/p018_X_s2/`; the source directory is untouched.
3. Stitch CT/circulation series at analysis time on the `step` column;
   check per-step CT jump at the seam ≤ 0.05% (restart-integrity gate,
   validated once in Phase 1).

Warm start requires construction-compatible mesh/NT/RPM. `Das` rotation across
restarts is handled (warmstart §2.5 fix — do not revert).

## Pre-submission gate (once per code state)

1. `julia --project -t 6 test/formulation_test.jl` passes (10 stages).
2. Coarse smoke: NT=6, ~6 revs, `RHPC_MESH=40_40`, both formulations, locally
   (≤6 threads; ~minutes). Runs clean, finite CT.
3. If `src/` changed: full `test/runtests.jl` + md5 re-verify on cluster.

## Local knob-name crib (driver: examples/rotor_hover_pressure_comparison.jl)

`RPM NT SETTLE_REVS TRUNCATION_DEPTH_R DAS_ETA_KINEMATIC DAS_MIN_DISPLACEMENT_R
DAS_KINEMATIC_ARC DAS_REFRESH NWAKEROWS PARTICLE_SHEDDING OVERLAP P_PER_STEP
MERGE_R_FACTOR MERGE_PARTICLES RHPC_FORMULATION RHPC_MESH RELAX_RLXF
RESTART_STEP RESTART_NAME RESTART_PATH BERNOULLI_ONLY SAVE_VTK`.
Also: `RHPC_BACKEND=fmm|direct` (default fmm; direct only for backend A/B
diagnostics, with `SFS_OFF=true` in BOTH arms — the SFS direct path has a
small-N FLOWVPM bug), `SFS_OFF` (diagnostics only, ruling 9).
Never set: `RELAX_FILTER_DOWNSTREAM_R` (ruling 2), `DAS_REFRESH=true`
(ruling 8). Driver hard-codes `max_particles=500_000` (watch in Phase 8) and
truncation radius 1.5R.

Gotchas learned in anger:
- **`SAVE_VTK=false` suppresses ALL file output** (`save_path=nothing` —
  no CT CSVs, no monitors). Smokes that need numbers must run `SAVE_VTK=true`.
- In coarse-NT smokes, scale η by 36/NT to hold physical Das (η is per-step).
- Nested-ssh MOTD banners can eat the first output line — md5-verify files
  individually when a batch comparison shows a phantom missing first entry.
