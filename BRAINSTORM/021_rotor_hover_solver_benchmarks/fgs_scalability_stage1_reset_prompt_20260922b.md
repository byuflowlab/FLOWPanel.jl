# Reset prompt: 021 FGS scalability — Stage 1, code done, smoke+submit owed (2026-09-22, supersedes fgs_scalability_stage1_reset_prompt_20260922.md)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

BRAINSTORM 021 is CLOSED AT PROMOTION; the active follow-on is the FGS
scalability diagnostic, plan = `fgs_scalability_diagnostic_plan_20260921c.md`
(plan C). Stage 0 is COMPLETE (`fgs_scalability_stage0_audit_20260922.md` —
read it in full: evidence base, baseline manifest, hypothesis table). Read
`CLAUDE.md` first, and `agent_policies/TESTING.md`/`HPC.md` before the
corresponding work below.

**Stage-1 CODE IS WRITTEN AND PARSE-CHECKED but UNCOMMITTED and UNSMOKED** in
both live checkouts (2026-09-22 session):

### FastMultipole (`~/Dropbox/research/projects/FastMultipole`, branch `flowpanel-20260817`, HEAD f4d6b671, src edits uncommitted)

- `src/solve_dagteam.jl`: (a) `build_dagteam_plan` gained kwarg
  `nworkers::Integer=Threads.nthreads()` — sizes `plan.xg/yb` (team size =
  `length(plan.xg)` everywhere: `dagteam_start_team!` spawns `length(xg)-1`,
  `dagteam_initialize!` partitions over `1:length(xg)`; capped-out workers are
  never spawned, no-deadlock argument in the docstring). (b) `dagteam_sweep!`
  and `dagteam_inner_sweeps!` gained an optional trailing `diagnostics`
  positional (default `nothing`) accumulating coordinator-only timers:
  `:dagteam_spawn_ns`/`:dagteam_join_ns` (team lifecycle, once per outer) and
  `:dagteam_wait_ns` (post-drain laggard wait)/`:dagteam_reduce_ns` (serial
  boundary reduction), once per sweep. All four are SUBSETS of
  `:nonself_product_ns` — existing keys' meaning unchanged.
- `src/solve.jl`: the 4 new keys added to the diagnostics init list
  (solve! ~line 1351); `dagteam_inner_sweeps!` call passes `diagnostics`;
  `FastGaussSeidel` constructor gained `dagteam_workers::Int=0` (0 = all
  threads) forwarded as `nworkers` to `build_dagteam_plan`.

### FLOWPanel (branch `fastmultipole`, uncommitted on top of 4de873c)

- `src/FLOWPanel_solver.jl`: `FGSSolver` gained field + kwarg
  `dagteam_workers::Int=0` (after `dagteam_precision`); forwarded to
  `FastMultipole.FastGaussSeidel` only when `sweep_order===:dagteam` AND
  nonzero (backward compat with pre-kwarg FastMultipole checkouts).
  `FGSPreconditioner` deliberately NOT extended.
- `benchmark/fgs_cold_common.jl`: config key `dagteam_workers` whitelisted
  (dagteam-only, Int ≥ 0) and forwarded by `cold_make`.
- `benchmark/fgs_r4_dagteam_stage1.jl` (NEW): Stage-1 per-process driver.
  One invocation = one fresh process. Env: `DAGTEAM_CONFIG` (calibrated
  dagteam_selected.toml), `STAGE1_ARM=fixed|accepted` (fixed forces
  max_iterations=27/tolerance=0, asserts inner=3 and iterations==27),
  `STAGE1_DIAG=0|1`, `STAGE1_SOLVES` (default 5), `STAGE1_BLOCK`,
  `STAGE1_LABEL`, `DAGTEAM_WORKERS`. Per solve: cold_trial row + `diag_*`
  columns (−1 when uninstrumented) → `stage1_solves.csv`; gates = BC
  certified-accepted + work count + 1e-8 repeat agreement. Also writes
  `numa_pages.csv` (page-placement verification from /proc/self/numa_maps
  after warmup), `residual_history.csv` (one extra recorded solve, validated
  against reference), `stage1_summary.toml`, `status.toml`.
- `benchmark/run_r4_fgs_stage1.slurm.sh` (NEW): single-task exclusive zen3
  launcher (128 cpus, 500G, `--constraint=zen3 --exclusive --qos=normal`,
  48 h), modeled on `run_r4_thread_scaling.slurm.sh`. Needs env
  `COLD_PROJECT`, `CAMPAIGN_PINS`, `COLD_DATA_ROOT` (+optional `CALIB_ROOT`,
  `CALIB_JOB`=13777133, `RESUME_FROM_JOB_ID`). Sequence: parse+precompile →
  fixed ladder j∈{1,16,32,64}×3 blocks (shuffled within block, champion
  placement `--interleave=0-3 --cpunodebind=0-3`, BLAS=1) → matched
  diagnostics ladder (same shape, STAGE1_DIAG=1) → accepted bridge
  16/32/64×1 block → placement A/B @j64 (champion vs `--interleave=0-7
  --cpunodebind=0-7`, 3 alternating pairs) → cap A/B @j64 (workers 16 vs 0,
  3 alternating pairs; cap=32 pairs auto-run only if cap16 paired median
  wins). Per-stage `STATUS_*` files, failures are findings (continue);
  stale-output move-aside on resume; background cpufreq sampler →
  `cpufreq.csv`. Asserts 128 cores + 8 NUMA nodes. Calibrated configs read
  from `$COLD_DATA_ROOT/thread-scaling-j<J>-13777133/fgs-calibrate/results/
  dagteam_selected.toml` (verify these exist for j=1,16,32,64 — launcher
  pre-checks).

`bash -n` and `Meta.parseall` pass on everything; NOTHING has been run.

## NEXT TASKS (in order; HPC submission pre-approved by Ryan for Stage 1 ONLY)

1. **Local smoke, ≤4 threads** (TESTING.md first): narrow FGS/solver checks
   (there is a runtests_unit_solver.jl; `:dagteam` has NO unit test — known
   ledger gap, pre-existing Kutta :jump failure is also known). Then a
   multi-thread progress check of the worker cap: small body, `-t 4`,
   FGSSolver(sweep_order=:dagteam, dagteam_workers=1|2) solves match
   dagteam_workers=0 and lexicographic (dagteam is deterministic at any team
   size). Also smoke the stage1 driver end-to-end on a small rung locally if
   cheap (the AB driver pattern ran R4-only; a full local R4 fixture may be
   too heavy — a 4-thread R1/R2 pass of the DRIVER mechanics needs a
   calibrated dagteam config, so it may be simpler to smoke FGSSolver kwargs
   + diagnostics dict directly in a REPL script and rely on the launcher's
   parse/precompile stages for wiring).
2. **Commit** both repos (FastMultipole src patches; FLOWPanel solver+harness
   +launcher+this file). Do not push without Ryan (standing ledger).
3. **Campaign setup** (global CLAUDE.md rules + HPC.md mechanics): annotated
   tags `campaign/p021-fgs-stage1-20260922` in FLOWPanel + FastMultipole
   (+FLOWVPM pin at its loaded commit — pin the TRIPLE); push branches+tags
   to origin (needed for orc fetch — this push is part of the approved
   campaign flow); on orc `git fetch origin --tags` +
   `scripts/prep_campaign_worktree.sh <tag> <dir>` per repo; campaign env
   (COLD_PROJECT) with Manifest dev-paths at the worktrees; CAMPAIGN_PINS
   toml (packages.{FLOWPanel,FastMultipole,FLOWVPM} path/tag/sha, schema
   consumed by `cold_packages` in fgs_cold_common.jl:272 — worktree
   deployment, annotated-tag check enforced at runtime); provenance file
   `fgs_scalability_stage1_provenance_20260922.md` in BRAINSTORM/021 with
   every pin BEFORE submitting.
4. **Submit**: `ssh orc -fN` if socket cold; slurm-availability skill with
   `--cpus 64 --mem-gb 500` (+`--eta`) to confirm a zen3 node; submit
   `benchmark/run_r4_fgs_stage1.slurm.sh` from the FLOWPanel campaign
   worktree top level with COLD_PROJECT/CAMPAIGN_PINS/COLD_DATA_ROOT
   (=`/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`) set
   via `sbatch --export`. Then monitor via `hpc-monitor` only.
5. **Write the analysis script** (owed by plan C: "ship the analysis script
   with the data") while the job runs: merges stage1_solves.csv across
   processes; process medians + dispersion; exclusive-phase reconciliation
   (phase sums vs diag_total_ns vs solve_seconds); T_phase(32)−T_phase(16),
   T_phase(64)−T_phase(32), and shortfall-from-doubling per phase (matched
   sample MEANS for additive differences, per plan C); instrumented-vs-
   uninstrumented overhead check (≤5% at each rung, unchanged shape);
   paired placement/cap effects. Plateau and regression = SEPARATE
   conclusions. Gate: if either effect fails to reproduce, stop and report
   comparison deltas — no interventions on a non-reproduced effect.

## Owed (carried)

- **R2-j1 (13829232_7, dir `r12-champion-r2-j1-13778533/`)**: was still
  RUNNING ~17 h into 48 h wall on 2026-09-22. Re-check via `hpc-monitor`;
  when `phase2/phase2.csv` lands, harvest via `harvester`, fill R2-j1 column.
- Standing Ryan-gated ledger unchanged (see
  `fgs_scalability_stage1_reset_prompt_20260922.md` §"Standing Ryan-gated
  ledger" — production-adoption ruling, origin pushes, notebook entries,
  WeakKeyDict/warmstart fix, :dagteam unit test, p018 archive approvals).

## Traps (all prior traps bind; new ones first)

- The fixed-work arm's rows have `solved=false`/`eligible=false` BY
  CONSTRUCTION (tolerance=0) — the driver gates on accepted+iterations==27+
  repeat instead. Never filter Stage-1 fixed rows on eligible/solved.
- Accepted solves run ONE extra fmm+influence+residual vs fixed-work
  (convergence detected at iteration 28's check); note in any accepted-vs-
  fixed comparison.
- The worker cap does NOT pin the 16 active workers to specific cores
  (Julia tasks float inside the 64-core cpuset) — the plan's "matching
  nested CPU set" is approximated by the champion cpuset; record as caveat.
- `diag_*` columns are −1 in uninstrumented rows; instrumented and
  uninstrumented rows must never be pooled.
- nfcache citations always carry cache-build (~94 s @ j8) + memory (~8.5 GB
  cache + ~9.8 GB state) caveats. `bc_certified` in phase2.csv is evaluator
  self-certification, NOT the ≤1e-6 test. Run-dir names case-sensitive
  (`R1-j32`). Task logs output-buffered — judge liveness by outputs/CPU,
  never log mtime; judge runs by outputs, never sacct. Resume dirs carry OLD
  job ids (`logs.before.<newid>/`). R4 thread-scaling BLAS=1. `ssh orc`
  needs live ControlMaster socket. Never edit source while a job uses its
  worktree; local runs NEVER >4 threads.

## House rules (binding)

Monitoring `hpc-monitor`; harvesting `harvester`; storage `hpc-storage`
(400 G cap — a VTK-free benchmark job, but archive pressure from other
campaigns may still fire); notebook writes Ryan-gated (offer, don't write);
dated status/provenance files in BRAINSTORM/021. HPC submission approved FOR
STAGE 1 ONLY; Stage 2/3 runs and optional-ledger reruns remain Ryan-gated.
