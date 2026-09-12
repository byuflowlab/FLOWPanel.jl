# FGS initialized-CPU optimization: v8 generation provenance and job 13653852

Prepared 2026-09-11 (America/Boise). Executes the first stage of
`fgs_initialized_cpu_optimization_plan_20260911.md`: R2 comparable baseline
(prepared-only), plan step 1 (fixed inner-iteration screen `{1,2,3,5,10}`),
and a longer accumulated R2 CPU profile for attribution. Prior context:
R1 pilot timings (job 13653115) and profiles (job 13653450, both solver
leaves completed and validated) are supporting evidence only.

## v8 harness generation (from the v7 pin, four files changed)

Local development worktree: `/private/tmp/flowpanel-cold-opt-20260911`,
branch `cold-opt-20260911`, commit `18f97b5` from
`campaign/p021-cold-source-20260912-v7` (`7cda26f`). Source tag
**`campaign/p021-cold-source-20260911-v8`** (`18f97b5`); branch and both v8
tags published to `byuflowlab/FLOWPanel.jl` (same pipeline Ryan approved for
the v7 tags).

Changes (`benchmark/fgs_cold_common.jl`, `benchmark/fgs_cold_README.md`,
`benchmark/run_cold_opt.slurm.sh` new, `test/runtests_benchmark_cold.jl`):

- `COLD_PREPARED_ONLY=1`: `cold_benchmark` skips the fresh scope
  (constructor-dominated; the plan excludes construction/warm starts).
  Prepared sampling, calibration, and every gate unchanged.
- `SCREEN_SET="key:v1,v2,...[;key2:...]"`: explicit one-factor screen roster
  off the rung seed, values parsed with the seed field's type (integer/Bool
  axes cannot promote to Float64); requires `STAGE=screen`, conflicts with
  `CONFIG_FILE`, unknown keys/types/ranges fail before filesystem work. The
  seed config always runs; generated candidates still calibrate their own
  stopping tolerance through the existing staircase + accuracy gate.
- Screen stage records a failed generated candidate in its `status.toml` and
  continues (e.g. inner=1 never crossing the gate is a finding, not a job
  kill). Baseline/verify/selected executions still stop at first failure.
  `selected.toml` is now written only when at least one candidate is eligible.
- `COLD_PROFILE_REPS=N`: CPU profile accumulates over N repeated prepared
  solves (resets outside recorded regions; final solution still validated).
  Allocation profile stays single-solve. Provenance records all three
  controls.
- Controls extended: SCREEN_SET roster/typing/rejection tests, new-env
  preflight rejects, and `ColdScreenFailureHarness` partial/total screen
  failure tests. Full control suite passed locally at j1/b1 (Julia 1.12.5,
  logic check only; the job re-runs controls at the pinned 1.11.7).

## Execution pins (job 13653852)

| Package | Worktree | Annotated execution tag | SHA |
|---|---|---|---|
| FLOWPanel | `/home/rander39/campaigns/p021-cold-opt-20260911-v8/FLOWPanel.jl` | `campaign/p021-cold-exec-20260911-v8` | `32fe1c439f701506eee1338e16f3b83e6fec2179` |
| FastMultipole | `/home/rander39/campaigns/p021-cold-20260910-v1/FastMultipole` | `campaign/p021-cold-exec-20260910-v1` | `ef10643a401d6da16e28be87805b67d11bdf1fb5` |
| FLOWVPM | `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl` | `campaign/p021-cold-exec-20260910-v1` | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` |

The exec commit is source v8 + the standard data-symlink prep commit
(`scripts/prep_campaign_worktree.sh`); the site-policy file already matched
the current `/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md`, so no
site-policy commit was needed. Worktree clean at submission; dependency
worktrees verified clean at their pinned SHAs. Dedicated environment
`/home/rander39/campaigns/p021-cold-opt-20260911-v8/{env,pins.toml}` — the
env is the v7 Project/Manifest with the FLOWPanel dev-path repointed at the
v8 worktree (FastMultipole/FLOWVPM paths unchanged).

## Job 13653852 (submitted 2026-09-11 evening, PENDING on m12 at submission)

`sbatch benchmark/run_cold_opt.slurm.sh` from the exec worktree with
`COLD_PROJECT`/`CAMPAIGN_PINS` above and
`COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`.
64 CPU exclusive zen3, 500 G, `--qos=normal`, 6 h (`--eta` estimated start
~2026-09-12 01:04; the fail-fast harness does not survive requeue, so no
standby). Defaults: `RUNG=R2`, `CONFIGS=fgs`,
`COLD_OPT_SCREEN_SET=inner:1,2,3,5`.

Output root: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/opt-13653852/`

Sequence (each its own process, per-process logs in the output root):
parse → precompile → controls j1/b1 + j4/b1 → smoke j4/b1 (R2 seed
P8/MAC0.4/leaf100/inner10, staircase-calibrated) → **prepared-only** seed
baselines j4/b1 and j64/b1 (`STAGE=verify` on the smoke `selected.toml`) →
screen j64/b1 (`STAGE=screen`, inner ∈ {1,2,3,5,10}, per-candidate
recalibration, prepared-only) → 10-rep accumulated CPU profile + allocation
profile of the seed (j64/b1, `COLD_INVESTIGATION=1`).

All existing gates frozen: BC relative L2 ≤ 1e-6 (certified FMM evaluator,
direct crosscheck on R2), evaluator disagreement ≤ 1e-7, repeat agreement
≤ 1e-8. A failing arm never retunes.

## Next actions for a fresh agent

1. Monitor 13653852 via `hpc-monitor` (judge by outputs, not sacct state).
   `COMPLETED` marker + per-process logs under the output root.
2. On completion: harvest `summary.csv`/`trials.csv` per candidate
   (`harvester`), rank inner-iteration roster by prepared total time to
   accepted accuracy including calibration-implied outer iterations
   (`iterations`, `estimated_inner_sweeps`, `estimated_fmm_passes` columns);
   retain the two fastest accepted candidates per plan step 2. Inspect the
   accumulated R2 `cpu_{flat,tree}.txt` for updated attribution (`.jls`
   deserialization only on a compute allocation).
3. Plan step 2 next: `SCREEN_SET="leaf:25,50,200"` (and P/MAC neighbors)
   around the step-1 winners — reuse `run_cold_opt.slurm.sh` with
   `COLD_OPT_SCREEN_SET` (new generation only if code changes).
4. A screen-candidate failure (esp. inner=1) is recorded in its
   `status.toml` — diagnose, don't retune gates.
5. R4 still has no `COLD_SEEDS` entry; add one (seeded from the R3 winner)
   before the plan's R4 stage — that requires a v9 source generation.
