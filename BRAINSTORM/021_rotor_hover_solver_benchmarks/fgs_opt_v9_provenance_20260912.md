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

## Job 13660643 result (stage A: leaf axis) — COMPLETED

COMPLETED 2026-09-12 17:17 (3h52m52s, m12-1-25, no restarts). All 8 candidates
`status=completed`, eligible, 10 reps; gates PASS (worst authoritative BC
rel-L2 6.68094917e-7 [inner=3/leaf=25]; all certified_fmm authoritative, no
direct fallback; repeat solution delta 0 everywhere). Harvested to
`fgs_opt_evidence_20260912/opt-13660643/`. Median prepared j64/b1 ranking:
i3/l100 11.160 < i5/l100 11.283 < i3/l50 11.576 < i3/l25 11.776 < i5/l50
11.975 < i5/l25 12.390 < i5/l200 12.761 < i3/l200 12.852. **leaf=100 optimal
on the leaf axis for both inners; retained bases = the two stage-A bases
themselves.** Stage-A recalibrated tolerances bit-identical to job-1 values
(deterministic calibration, same node). In-job seed baselines: j4/b1 17.129 s,
j64/b1 10.968 s medians (consistent with job 1).

## Launch settings (jobs 3+4: stage B — P/MAC neighbors, 2-way per-base split)

Config-only reruns on clean v9 (worktree verified clean at `f03ab18` before
submission). Ryan-approved 2-way split: one job per retained base, each
self-anchored (own CONFIG_FILE seed baselines) so ranking comparisons stay
within-node. Bases files (one config each, stage-A calibrated tolerances):

- `/home/rander39/campaigns/p021-cold-opt-20260912-v9/bases_r4_i3l100_20260912.toml`
  md5 `95a946426d85c37eede07bfcb3427765` (inner=3/leaf=100, tol 3.479128881193055e-7)
- `/home/rander39/campaigns/p021-cold-opt-20260912-v9/bases_r4_i5l100_20260912.toml`
  md5 `f0d2d8e26ae7ef0e0fe64d20b585687d` (inner=5/leaf=100, tol 2.5712316195808637e-7)

Env identical to job 2 (COLD_PROJECT/CAMPAIGN_PINS/COLD_DATA_ROOT,
COLD_OPT_RUNG=R4, COLD_OPT_STAGE=screen_profile, CONFIG_FILE=
opt-13657404/smoke-j4-b1/selected.toml, COLD_MIN_REPS=10, COLD_OPT_PROFILE_REPS
unset) except per job: `COLD_OPT_SCREEN_BASE_FILE=<bases file above>` and
`COLD_OPT_SCREEN_SET="P:6,10;MAC:0.3,0.5"` (5 candidates/job: base + 4
one-factor neighbors). Walltime 6 h (2 baselines ~28 min + 5×~20 min screen +
~35 min stock 10-rep profile ≈ 3.2 h + margin). sbatch --test-only first.

- Job **13661797** = stage B base inner=3/leaf=100 (`bases_r4_i3l100_20260912.toml`),
  job **13661798** = stage B base inner=5/leaf=100 (`bases_r4_i5l100_20260912.toml`).
  Both submitted 2026-09-12 ~17:23, 6 h walltime; `--test-only` estimated
  immediate start. Output roots: `.../opt-13661797/`, `.../opt-13661798/`.
- Nodes at start: 13661797 -> m12-1-17, 13661798 -> m12-1-25 (both RUNNING 17:22:58).

## Jobs 13661797/13661798 results (stage B: P/MAC axes) — COMPLETED

13661797 (i3/l100 base, m12-1-17) 2h40m45s; 13661798 (i5/l100 base, m12-1-25)
2h37m38s. Both COMPLETED markers present; harvested to
`fgs_opt_evidence_20260912/opt-1366179{7,8}/`. Per job: base + P10 + MAC0.3
completed (all gates PASS: certified_fmm authoritative, repeat delta 0, worst
BC rel-L2 6.09e-7); **P6 and MAC0.5 FAILED calibration in both jobs** ("FGS
staircase has no certified crossing with a decreasing successor",
fgs_cold_common.jl:442) — looser far field cannot certify 1e-6 at R4; findings,
outputs retained. Median prepared j64/b1, within-node: 13661797 base 11.014 <
P10 11.969 < MAC0.3 14.971 (anchor baselines j4 17.261 / j64 11.212);
13661798 base 11.126 < P10 11.792 < MAC0.3 15.376 (anchors 16.844 / 11.081).
Recalibrated tolerances again bit-identical to prior values on both nodes.
**Finalists = the two bases: i3/l100 (= the seed) and i5/l100. Retained
configuration = seed P8/MAC0.4/leaf100/inner3 → step-4 reprofile NOT needed
(seed already profiled at 20 reps in job 13657404).**

## Launch settings (jobs 5+6: finalist confirmation, 2 independent nodes)

Config-only reruns on clean v9 (worktree verified clean at `f03ab18`).
`confirm_r4_20260912.toml` md5 `800aab6a535e91f5aebc0d2b8efd0e88` at
`/home/rander39/campaigns/p021-cold-opt-20260912-v9/`: configs =
[i3/l100 tol 3.479128881193055e-7, i5/l100 tol 2.5712316195808637e-7].
`cold_run` times both sequentially in one process per baseline stage (>=10
unprofiled trials each, COLD_MIN_REPS=10), j4/b1 then j64/b1; screen stage =
trivial roster reusing stage-A bases file `bases_r4_leaf_20260912.toml`
(md5 `be538fc013863860e9a181f538076649`) with `COLD_OPT_SCREEN_SET=leaf:100`
(neighbor == base value → dedupe → the 2 finalists recalibrated once more);
end profile stage profiles BOTH confirm.toml configs at stock 10 reps.
Env per job (identical; two simultaneous exclusive jobs → distinct nodes):

```bash
export COLD_PROJECT=/home/rander39/campaigns/p021-cold-opt-20260912-v9/env
export CAMPAIGN_PINS=/home/rander39/campaigns/p021-cold-opt-20260912-v9/pins.toml
export COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910
export COLD_OPT_RUNG=R4 COLD_OPT_STAGE=screen_profile
export CONFIG_FILE=/home/rander39/campaigns/p021-cold-opt-20260912-v9/confirm_r4_20260912.toml
export COLD_OPT_SCREEN_BASE_FILE=/home/rander39/campaigns/p021-cold-opt-20260912-v9/bases_r4_leaf_20260912.toml
export COLD_OPT_SCREEN_SET=leaf:100
export COLD_MIN_REPS=10
unset SCREEN_SET SCREEN_BASE_FILE COLD_OPT_PROFILE_REPS
sbatch --time=06:00:00 benchmark/run_cold_opt.slurm.sh   # x2
```

Walltime: baselines ~2x14 min each stage x2 configs + screen 2x~20 min +
profile 2x~25 min ≈ 3 h; 6 h margin.
- Jobs **13663309** (m12-1-17) and **13663310** (m12-1-25) submitted 2026-09-12 20:07, both RUNNING immediately. Output roots: `.../opt-1366330{9,10}/` (opt-13663309, opt-13663310).

## Jobs 13663309/13663310 results (finalist confirmation) — COMPLETED

13663309 (m12-1-17) 2h56m27s; 13663310 (m12-1-25) 2h55m41s. Markers present;
harvested to `fgs_opt_evidence_20260912/opt-1366330{9,10}/` (287 MB each —
two 10-rep profiles per job, one per finalist). All 12 batches (2 nodes x
{j4 baseline, j64 baseline, j64 screen recalibration} x 2 configs) completed;
ALL gates PASS (certified_fmm authoritative everywhere, repeat deltas 0,
accepted/eligible true). Median prepared seconds (config file order i3 then
i5; batches sequential in one process per stage):

| batch | i3/l100 | i5/l100 |
|---|---:|---:|
| 13663309 j4 baseline | 16.897 | 14.653 |
| 13663309 j64 baseline | 11.402 | 11.357 |
| 13663309 j64 screen (recal) | 10.848 | 11.129 |
| 13663310 j4 baseline | 17.066 | 15.501 |
| 13663310 j64 baseline | 11.116 | 11.021 |
| 13663310 j64 screen (recal) | 11.099 | 11.082 |

**Verdict: at j64/b1 the finalists are a statistical tie** (i5 ahead in 3 of
4 paired batches by 0.15–0.9%, i3 ahead once by 2.5%; differences within
batch spread and node variance). **At j4/b1 inner=5 is decisively ~13%
faster on both nodes** (fewer outer iterations dominate at low thread
count). Retained configuration for the code-opt handoff stays the seed
i3/l100 (median winner in every j64 screen ranking: jobs 13657404, 13660643,
13661797); i5/l100 recorded as within-noise equivalent at j64 and preferred
at low thread counts. Both finalists' 10-rep profiles exist in both
confirmation jobs; seed 20-rep profile remains opt-13657404's.

Confirmation harvest audited (confirm_table.csv vs raw summary/convergence
CSVs: match). Final deliverable:
`fgs_opt_r4_diagnostics_package_20260912.md` (configuration conclusion,
measured profile picture, ranked >=5% opportunities, extra measurements).
R4 diagnostics phase COMPLETE 2026-09-12; no solver implementation changed;
no notebook entry written (awaiting Ryan's approval).
