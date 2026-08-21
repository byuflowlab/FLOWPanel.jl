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
- **≤20 active study jobs** (Ryan raised the cap 6→10 on 2026-08-04, 10→20 on
  2026-08-06; NOTE: the 019 session's auto-submitter still checks against 10 —
  update or supersede it before relying on it) —
  check `squeue` before every submission.
- Final-settings cases (`p018_trunc6`, `p018_green`, `p018_nomerge`,
  `p018_final`) take Das\*/NT\*/σ\* via environment at submission, e.g.
  `sbatch --export=ALL,DAS_ETA_KINEMATIC=1.0,NT=36,OVERLAP=2.4,P_PER_STEP=11,MERGE_R_FACTOR=0.0053 ...`
  (the `p018_*` carrier honors env overrides for every knob).

## Wake tripwire (added 2026-08-05, phase 13 S0a)

`WakeHealthMonitor` is appended to the driver's monitor tuple by default
(`WAKE_HEALTH=false` disables). It writes
`data/<run>/monitors/<run>_monitor04_wake_health_system1.csv` per step:

| column | read it for |
|:--|:--|
| `n_particles` | cost/cap tracking against the driver's `max_particles` |
| `max_u` | divergence. `p018_L2` went 733 → 37,088 m/s in **two** steps |
| `min_sigma_ratio` | **the early-warning column.** `min sigma / shed sigma`; the rVPM core collapse trends here for hundreds of steps before `max_u` moves |
| `max_gamma_over_sigma2` | velocity-runaway proxy, `\|u\| ~ \|Gamma\|/sigma^2` |
| `wall_s` | per-step wall clock — use it to size time requests instead of guessing |

Appended LAST, so monitor indices 02/03 and every existing CSV filename are
unchanged; verified bit-identical against a `WAKE_HEALTH=false` control.
Freshly shed particles sit at `min_sigma_ratio = 1.0` exactly, so any sustained
decline is contraction, not noise. Per-step wall clock also now prints in the
step line, and stdout is flushed (it was block-buffered, which is why long runs
looked hung).

**Use it as a stop rule, not a post-mortem**: a run whose `min_sigma_ratio` is
falling steadily is heading for the L2 failure, and more memory will not help
(`p018_L2_s2` reproduced the blow-up step bit-identically at 128 G).

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
- **Retention amendment (2026-08-04, after job 13036477 died on a deleted
  restart file): every VTK sweep must preserve one COMPLETE restart set per
  run that could seed a warm start** — `<run>_body1.pvd`, body
  `.<S>.vtu`, wake `.{1,2}.<S>.vts`, particles `.<S>.vtp` at the same last
  restartable step S. `.vtm` stubs without their pieces are worse than
  nothing: they make a run look restartable when it is not. (The warmstart
  loader reads exactly those four path patterns — `src/FLOWPanel_warmstart.jl`;
  filaments are not read.)
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

**Running the gate ON the cluster (learned in anger, 2026-08-12, jobs
13150992–13157614; template `examples/run_p15_gate_hpc.slurm.sh`):**

- **Julia version:** `module load julia` gives 1.12.6, which explodes
  precompile against the 1.11.7-pinned Manifest (PythonPlot/MbedTLS_jll).
  Use the production launcher's spack julia-1.11.7 PATH fallback, and export
  `MPLBACKEND=Agg`.
- **Deploy `test/` and `examples/` in LOCKSTEP with `src/` — an md5 sweep of
  `src/` alone is NOT a sufficient deploy check.** The full suite had never
  been runnable on the cluster: six test files were stale vs HEAD, two
  tracked kutta test files and two `ssw_*` example files had never been
  copied, one untracked fixture (`test/data/legacy_wake_conversion_
  reference.jl`) was local-only, and `examples/suddenly_started_wing.jl` was
  a July-24 copy whose `SSWConfig` predated its test. Each missing file
  costs one ~20-min job iteration to find.
- **Before submitting a suite job, scan its file dependencies:** every
  `include(...)`/`joinpath(...)` `.jl` literal in `test/*.jl` (including
  `pnl.examples_path` references) must exist on the cluster; untracked
  fixtures do not arrive via git.
- **Known exclusion (2026-08-12):** `runtests_example_suddenly_started_wing
  .jl` cannot load in the current project env anywhere —
  `examples/helper_functions.jl` imports GeometricTools, removed from
  Project.toml deps in May 2026 (commit 7d52be0). The cluster `runtests.jl`
  carries a commented-out include with this note. Resolution (restore the
  dep vs port helper_functions.jl off `gt`) is Ryan's call.

## Walltime allowance

Ryan 2026-08-12: submissions may request **up to 72 h** when needed
(supersedes the de-facto 48 h ceiling). Prefer the shortest walltime that
avoids a restart chain.

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

## Post-submission knob verification (MANDATORY — added 2026-08-03)

**Every submission must be verified against what the job itself reports**,
before it burns hours. Within seconds of starting, each run prints a banner to
`logs/slurm/slurm-<name>-<jobid>.out`:

```
  mesh:45_185_ct4 formulation:velocity RPM:5400 NT:72 depth:4R rlxf:0.15
  relax_filter:offR das_eta:2.0 overlap:2.0 pps:2 merge_r:0.0120 settle:12 ...
```

Check `NT`, `das_eta`, `overlap`, `pps`, `merge_r`, `rlxf`, `nwakerows`, `mesh`
against intent. (`_case_metadata.toml` carries the same values but is written
much later.)

**Why this is not optional — the launcher's override precedence is
INCONSISTENT between case arms, and neither direction warns:**

- `p018_nt72` / `p018_nt144` use `export DAS_ETA_KINEMATIC="${VAR:-2.0}"` —
  the **environment wins**. `sbatch --export=ALL` propagates the submitting
  shell, so a stale `DAS_ETA_KINEMATIC=1.0` left over from an earlier
  submission silently halved Das and cost job 13011373 a 36.9 h run at the
  wrong operating point (2026-08-02).
- `p018_L1|p018_L1_warm` use `export OVERLAP=2.4` — **unconditional**, so the
  case arm wins and an `--export` override is silently discarded.

Rules that follow:

1. Verify the banner after every submission. No exceptions.
2. To vary a knob that its case arm sets unconditionally, **add a new case
   arm** (as `p018_L1_ov3` was added) rather than trying to override it.
3. Prefer a clean submitting shell; do not rely on `--export=ALL` inheritance
   to carry intent.

## Trailing-edge selection and the root cap (2026-08-03)

**Cap protection is `end_node`, not the radial cutoff.** Earlier guidance in this
file said the cutoff guarded the cap and should be auto-raised; that was wrong and
is superseded — raising it silently discarded 3 genuine TE edges at the hub.

How it works now in `examples/rotor_hover_pressure_comparison.jl`:

1. `find_dji9443_trailing_edge_indices` returns `[outer, second_outer, inner]`.
   The driver passes `inner` as `end_node`, so `calc_shedding_from_seed` traces
   exactly the detector's validated chain and **errors** if it cannot reach it.
   The bbox is reduced to blade separation (cutoff 0).
2. The root clip is applied **after** tracing, as an explicit modeling choice,
   using the same edge-midpoint criterion the bbox used. `SHEDDING_R_OVER_R`
   (default 0.1) is therefore ONLY a modeling knob now. It is inert on any hub
   variant (blade root outboard of it), so hub meshes shed their full TE.
3. A regression guard aborts if any shed edge runs **circumferentially** (a cap
   edge); it tests edge direction, not radius, because the chain now legitimately
   reaches the blade root.

Why it matters: selecting TE edges by *where they are* fails silently in both
directions — too small a cutoff sheds off the root cap (92 of 136 edges on the
0.15R hub mesh; that run diverged in 17 steps, job 13031568), too large drops
real trailing edge. Selecting by *what they are* (sharp, radial, anchored at both
ends) fails loudly instead.

Stock behavior is unchanged and verified bit-identical (40 edges from r/R 0.111;
5-step CT matches to the last digit).

## OOM triage, corrected 2026-08-03

The rule "OOM inside `merge_particles!` usually means the simulation blew up"
still holds, but **both** of its stated discriminators were wrong in the one case
that mattered (`p018_L2`, job 13029924 — clean run, genuine memory):

1. **Judge |CT| from `data/<run>/monitors/*force_system1.csv`, not the log.**
   That CSV is written incrementally and survives the SIGKILL; the stdout log is
   block-buffered and had **zero** `CF = ` lines. Thrust ≈ +0.068–0.076 over 99
   steps settled it.
2. **`MaxRSS ≈ ReqMem` does NOT establish a memory ceiling.** L2 was OOM-killed
   with `MaxRSS = 22.3 GB` against a **64 G** request (35%). `sacct` samples
   periodically and misses fast allocation spikes (merge hash / FMM rebuild at
   ~350k particles). A low MaxRSS is not evidence against a memory problem.

σ-L2 and finer may be submitted with `--mem=128G` (ruling 5, amended).
