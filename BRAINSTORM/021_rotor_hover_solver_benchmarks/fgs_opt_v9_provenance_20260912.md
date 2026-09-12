# FGS initialized-CPU optimization: v9 generation provenance (R4 diagnostics)

Prepared 2026-09-12 (America/Boise). Executes the R4 redirect in
`fgs_opt_r4_diagnostics_handoff_20260912.md` under the governing plan
`fgs_initialized_cpu_optimization_plan_20260911.md` (User steering 2026-09-12):
skip R3 and further R2 tuning; R4 control + baselines + inner screen
`{1,2,3,5,10}` + accumulated CPU/allocation profiles on the 58,192-panel
`65_209` mesh. Provisional R4 seed P8/MAC0.4/leaf100/inner3 (R2-derived, NOT
an established R4 optimum).

## v9 harness generation (from the v8 pin, four files changed)

Local worktree `/private/tmp/flowpanel-cold-opt-20260912-v9`, branch
`cold-opt-20260912-v9`. Source tag **`campaign/p021-cold-source-20260912-v9`**
(`721235e86f8db3ecbe6d3c6c00514b0de5172baf`), parent v8 source
`campaign/p021-cold-source-20260911-v8` (`18f97b5`). Branch and both v9 tags
pushed to `byuflowlab/FLOWPanel.jl` 2026-09-12 (exec tag created on ORC,
fetched locally, pushed from the local machine — ORC has no GitHub push auth).

Changes (`benchmark/fgs_cold_common.jl`, `benchmark/run_cold_opt.slurm.sh`,
`benchmark/fgs_cold_README.md`, `test/runtests_benchmark_cold.jl`); no solver
implementation changed. See the handoff for the six-point change list:
R4 rung support (explicit seed, no TUNED lookup), certified evaluator with
reported direct fallback, `SCREEN_BASE_FILE` saved-base screens with
tolerance reset, `COLD_MIN_REPS`, `COLD_OPT_PROFILE_REPS`,
`COLD_OPT_STAGE=attribution`. Local checks: parse of all nine driver/control
files, stdlib-only logic tests, shell syntax, `git diff --check`. Full
controls did NOT run locally (no Krylov in local env); the job runs them at
pinned Julia 1.11.7 j1/b1 and j4/b1 before any measurement.

## Execution pins

| Package | Worktree | Annotated execution tag | SHA |
|---|---|---|---|
| FLOWPanel | `/home/rander39/campaigns/p021-cold-opt-20260912-v9/FLOWPanel.jl` | `campaign/p021-cold-exec-20260912-v9` | `f03ab18a7841e9a38876603a48b1799d83fd1e41` |
| FastMultipole | `/home/rander39/campaigns/p021-cold-20260910-v1/FastMultipole` | `campaign/p021-cold-exec-20260910-v1` | `ef10643a401d6da16e28be87805b67d11bdf1fb5` |
| FLOWVPM | `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl` | `campaign/p021-cold-exec-20260910-v1` | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` |

Exec commit = source v9 + standard data-symlink prep commit
(`scripts/prep_campaign_worktree.sh`, data -> `~/projects/FLOWPanel.jl/data`).
FLOWPanel worktree clean at `f03ab18`; dependency worktrees verified clean at
their pinned SHAs (unchanged from v8/v1). Environment
`/home/rander39/campaigns/p021-cold-opt-20260912-v9/{env,pins.toml}`: the v8
Project/Manifest with the FLOWPanel dev-path repointed at the v9 worktree
(FastMultipole/FLOWVPM dev-paths unchanged). Julia 1.11.7 pinned
(`module load cuda/12.8.1-zkkfiog julia/1.11.7-6bmogfl`).

- `env/Manifest.toml` md5 `b12470a55bc7fea7b3e727d7355bd7c7`
- `env/Project.toml` md5 `1dca07d49de4ea8894eadd4f3df7dbde`

## Launch settings (job 1: R4 all-stage — controls, smoke, baselines, inner screen, profile)

Submitted from the clean v9 FLOWPanel worktree root:

```bash
export COLD_PROJECT=/home/rander39/campaigns/p021-cold-opt-20260912-v9/env
export CAMPAIGN_PINS=/home/rander39/campaigns/p021-cold-opt-20260912-v9/pins.toml
export COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910
export COLD_OPT_RUNG=R4 COLD_OPT_STAGE=all
export COLD_OPT_SCREEN_SET=inner:1,2,5,10   # seed inner=3 always runs; 10 explicit
export COLD_MIN_REPS=10 COLD_OPT_PROFILE_REPS=20
unset CONFIG_FILE SCREEN_SET SCREEN_BASE_FILE COLD_OPT_SCREEN_BASE_FILE
sbatch --time=12:00:00 benchmark/run_cold_opt.slurm.sh
```

Resources: 1 node, 64 cpus-per-task, `--mem=500G`, `--constraint=zen3`,
`--exclusive`, `--qos=normal` (launcher directives); walltime raised from the
6 h default to 12 h because R4 is 3.69x R2 panels (58,192 vs 15,760) and the
completed R2 all-stage job (13653852) took 46m23s — quadratic-scaling stages
could approach 10+ h. Availability probe 2026-09-12: m12 28 fitting idle
nodes, `--test-only` estimated immediate start at 12 h walltime, maxtime
3-00:00:00. Partition left unspecified per site policy; zen3 constraint
routes the job. Launcher pins 64 distinct physical cores via taskset;
BLAS=1 in timed processes; profiling isolated from unprofiled timing.

Storage preflight: `df` home proxy 196,156 MiB used of 2 TiB (~191.6 GiB
against the 400 G FLOWPanel cap); R2 evidence 84 MiB, no VTK in this
campaign's outputs.

## Job record

- Job **13657038**: FAILED at 2m49s in `cold_precompile.jl` — the first
  pins.toml was written with flat `[FLOWPanel]`-style sections, but
  `cold_packages()` requires a top-level `[packages]` table
  (`TOML.parsefile(CAMPAIGN_PINS)["packages"]`). Parse stage passed; no
  measurement ran. Config-only error; pins.toml is untracked, so no new
  generation. Output dir `opt-13657038/` retained (provenance snapshots only).
- Job **13657404** (submitted 2026-09-12, 12 h walltime): resubmission with
  corrected pins.toml (`[packages.FLOWPanel]` / `[packages.FastMultipole]` /
  `[packages.FLOWVPM]`, schema diff-checked against v8). Same launch env as
  above, same worktree/env/tags (unchanged, still clean).
- Output root: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/opt-13657404/`
- Job 13657404 **COMPLETED** 2026-09-12 12:54 (3h53m). All gates PASS; five
  inner candidates accepted; harvested to
  `fgs_opt_evidence_20260912/opt-13657404/` (see harvest_summary.md). Median
  winners inner=3 (10.709 s) and inner=5 (11.005 s) at j64/b1 prepared;
  selected.toml (min-time) agrees on inner=3. R4 profile: dense nonself GEMV
  ~77% of snapshots; FMM passes ~1.3%.

## Launch settings (job 2: R4 leaf-axis screen from saved bases — handoff step 7, stage A)

Config-only rerun on clean v9 (no harness change; same worktree, env, pins,
tags as job 1). Bases file
`/home/rander39/campaigns/p021-cold-opt-20260912-v9/bases_r4_leaf_20260912.toml`:
the two median winners from opt-13657404 (inner=3 tol 3.479128881193055e-7,
inner=5 tol 2.5712316195808637e-7), full config tables copied from their
screen `config.toml`s. Harness resets FGS tolerance to 0 per point
(recalibration); bases file is copied+sha256'd into process provenance.
Roster: 2 bases + leaf {25,50,200} one-factor neighbors each = 8 candidates
(leaf=100 covered by the bases themselves).

```bash
export COLD_PROJECT=/home/rander39/campaigns/p021-cold-opt-20260912-v9/env
export CAMPAIGN_PINS=/home/rander39/campaigns/p021-cold-opt-20260912-v9/pins.toml
export COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910
export COLD_OPT_RUNG=R4 COLD_OPT_STAGE=screen_profile
export CONFIG_FILE=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/opt-13657404/smoke-j4-b1/selected.toml
export COLD_OPT_SCREEN_SET=leaf:25,50,200
export COLD_OPT_SCREEN_BASE_FILE=/home/rander39/campaigns/p021-cold-opt-20260912-v9/bases_r4_leaf_20260912.toml
export COLD_MIN_REPS=10
unset SCREEN_SET SCREEN_BASE_FILE COLD_OPT_PROFILE_REPS
sbatch --time=08:00:00 benchmark/run_cold_opt.slurm.sh
```

`screen_profile` skips parse/precompile/controls/smoke (already validated in
job 1 on this exact worktree/env) and runs prepared j4/j64 baselines from
CONFIG_FILE, the saved-base screen, then the stock 10-rep profile of the
CONFIG_FILE seed. Walltime 8 h: job-1 screen paced ~20 min/candidate × 8
candidates + 2 baselines + profile ≈ 4–5 h with margin (leaf=25/200
calibration cost unknown at R4).

- Job **13660643** submitted 2026-09-12 ~13:25; `--test-only` estimated
  immediate start on m12-1-25 (same node as job 1 — matched hardware).
  Bases file md5 `be538fc013863860e9a181f538076649`; worktree verified clean
  at `f03ab18` immediately before submission. Output root:
  `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/opt-13660643/`
