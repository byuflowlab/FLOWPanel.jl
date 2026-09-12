# Resume cold-harness validation and R1 HPC pilot

Ryan paused the previous session. Resume implementation and verification, not merely planning. The authorized scope is **harness validation plus an R1 HPC pilot**, not solver optimization or broader screening. No validated timing comparison exists yet.

## Read first

Read `/Users/ryan/.claude/CLAUDE.md`, repo `CLAUDE.md`, and `agent_policies/{WORKFLOW,TESTING,HPC}.md`. Then read:

1. `BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_cold_validation_pilot_plan_20260910.md` — saved user-intent plan.
2. This handoff — latest state, supersedes older progress status.
3. `BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_cold_validation_review_20260910.md` — development snapshots and initial campaign record; its original pending-job section is historical.

Read current ORC `/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md` and applicable `ai-docs/` policies before cluster work. SSH sandbox failures were resolved with `exec_command(sandbox_permissions="require_escalated")`, not by repairing SSH. Use existing `ssh -o BatchMode=yes -o ConnectTimeout=8 orc`; stop if actual authentication/MFA is unavailable.

**All Julia execution stays on HPC compute nodes**, including parsing, tests and precompilation. Another agent is using the local machine. No local Julia was run this session. Use explicitly lower-cost agents with minimal context for bounded test execution/monitoring/harvesting; primary owns edits and scientific interpretation. Avoid monitoring agents that withhold initial status for a long polling loop. Site scheduler checks must be at least 60 seconds apart; inspect only our jobs.

## Worktrees and commits — preserve everything

The shared Dropbox checkout remains dirty and contains the OLD draft harness. **Do not use its benchmark files as the latest implementation.** It was not cleaned or overwritten. Historical results and unrelated changes remain intact.

Latest local implementation:

- `/tmp/flowpanel-cold-20260910`, branch `cold-harness-20260910`, clean HEAD `3194a2c`.
- `/tmp/fastmultipole-cold-20260910`, branch `cold-pilot-20260910`, source snapshot `720472ee`.
- `/tmp/flowvpm-cold-20260910`, branch `cold-pilot-20260910`, source snapshot `d7402d4`.

FLOWPanel commits:

- `575f2a0`: copies current dirty runtime source into an isolated snapshot (existing formulation/wake/warmstart edits).
- `b7e682e`: cold harness implementation, tests, docs, launchers.
- `cd4d191`: validates inherited fixture controls and adds staged pilot execution.
- `3194a2c`: fixes ORC shell initialization under `set -u`.

Dependency snapshots copy existing runtime source changes: FMM farfield-output hooks; FLOWVPM splitting/timeintegration. This session made no solver API or optimization edits. Unrelated example/test/result changes were not copied into these snapshots.

Published FLOWPanel annotated source tags: `campaign/p021-cold-source-20260910-v1` (`b7e682e`), `...-v2` (`cd4d191`), `...-v3` (`3194a2c`). Dependency annotated source tags `campaign/p021-cold-20260910-v1` are published in each repository.

## What was implemented (UNVALIDATED in Julia)

Read the actual files in `/tmp/flowpanel-cold-20260910`:

- `benchmark/fgs_cold_common.jl`: schema/type/finite/range validation; early paths/memory/thread/mesh checks; immutable selected settings; calibrated seeds; explicit evaluator acceptance; warmup/direct evidence; constructor tuple release; reset checks; independent constructor smoke; thread checks; profile solution verification; fail-fast execution; pinned loaded-package provenance and copied Manifest.
- `benchmark/common.jl`: separate BENCH_BLAS_THREADS and a real GEMM before thread banner.
- `benchmark/phase1_case.jl`: independent fixture root; restore frozen velocity before setting BC sources.
- `benchmark/rotor_hover_solver_cold{,_smoke}.jl`: timing and smoke entry points.
- `benchmark/rotor_hover_solver_phase2_profile.jl`: opt-in cold branch.
- `benchmark/run_cold_process.sh`: separate startup Julia/BLAS/OpenMP controls.
- `benchmark/run_cold_pilot.slurm.sh`: sequential precompile, controls, smoke, timings and profiles; supports `COLD_PILOT_STAGE=all`, `controls_smoke`, or `timing_profiles`.
- `benchmark/cold_precompile.jl`.
- `benchmark/fgs_cold_README.md`.
- `test/runtests_benchmark_cold.jl`: BLAS work, schema/immutability, invalid inputs without filesystem effects, acceptance semantics, injected failure propagation.

Only `bash -n` and `git diff --check` passed. **No Julia parse, control test, smoke, timing or profile has run successfully.** Review the new implementation critically; do not mistake it for validated code.

## Immediate likely parser defect — fix BEFORE another submission

Just before pausing, primary noticed this new command literal in `cold_packages`:

    readchomp(`git -C $path rev-parse $tag^{commit}`)

Julia command literals require unquoted `{}` to be escaped/quoted. Local source inspection (no Julia execution) confirmed `shell_special = "#{}()[]<>|&*?~;"` in `/Applications/Julia-1.11.app/Contents/Resources/julia/share/julia/base/shell.jl`. Treat this as a likely parse blocker; it has **not been fixed or runtime-tested**. Construct the revision string separately (e.g. `revision = tag * "^{commit}"`) and interpolate `$revision`. Inspect for similar command literal defects. Validate syntax on an HPC compute node before expensive precompilation. Do not run local Julia.

Other implementation concerns to review during validation: fixture preflight and inherited controls; memory/thread gates around calibration; persistent evidence from excluded trials; selected tolerance equality and failure propagation; no measured-component claims. Keep inferred work counts labeled estimates (`-1` unavailable), retained object size separate from lifetime RSS, and `begin_step_solution!` outside timed solves.

## Actual cluster state at pause

One job was submitted: **13642618**, from campaign v1. It FAILED BEFORE JULIA because `/etc/profile` reads unset `HISTCONTROL` under `set -u`.

Preserved log:
`/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/slurm-13642618.err`

Exact error: `/etc/profile: line 50: HISTCONTROL: unbound variable`.
No process generation/control logs were produced. No replacement job was submitted. The shell fix (`set +u; source /etc/profile; set -u; module load julia/1.11.7-6bmogfl`) was validated in a lightweight login-shell check using `command -v julia`; it succeeded. That check did not execute Julia.

All preparation commands finished before pause; no intentional background command/monitor remains from this task. Two pre-existing jobs 13593020/13593021 were running from `/home/rander39/wt021/FLOWPanel.jl-c`; they are unrelated and must not be edited/cancelled. Refresh their status only if necessary for concurrency.

### Execution worktrees and tags

v1 root: `/home/rander39/campaigns/p021-cold-20260910-v1`

All v1 execution tags are annotated and published as `campaign/p021-cold-exec-20260910-v1`:

| Package | Execution SHA |
|---|---|
| FLOWPanel | `029b9f3d30bba706dea77418281f0ce4d75a4689` |
| FastMultipole | `ef10643a401d6da16e28be87805b67d11bdf1fb5` |
| FLOWVPM | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` |

v3 preparation ALSO completed, but no v3 job was submitted:

- FLOWPanel worktree: `/home/rander39/campaigns/p021-cold-20260910-v3/FLOWPanel.jl`.
- Actual HEAD: `e3ae9bd89b51d1cb86500b623a2a99195cdffa04`.
- Annotated execution tag: `campaign/p021-cold-exec-20260910-v3`, **currently only on ORC, not yet fetched/published locally**.
- Env and pins: `/home/rander39/campaigns/p021-cold-20260910-v3/{env,pins.toml}`.
- v3 Manifest/pins point FLOWPanel at v3, with unchanged dependency worktrees/tags in v1 above.

The preparation helper creates a data-symlink commit after the input source tag. v3: helper commit `3f46593`, then ORC policy-symlink commit `e3ae9bd`. Always tag the ACTUAL final execution commit. Every loaded dev package must be a clean pinned git worktree.

Cluster GitHub push has no interactive credentials. Do not copy credentials. Fetch execution tags from ORC into the authenticated local repository, then push from local, e.g.:

    git -C /tmp/flowpanel-cold-20260910 fetch orc:/home/rander39/projects/FLOWPanel.jl refs/tags/TAG:refs/tags/TAG
    git -C /tmp/flowpanel-cold-20260910 push origin refs/tags/TAG

Source deployment uses published git tags and cluster git fetch, never scp source then commit. Preparation scripts used this session are `/tmp/prep_cold_remote.sh`, `/tmp/finish_cold_remote.sh`, `/tmp/prep_cold_v3.sh`; logs alongside them. Read before reuse; first scripts deliberately refuse existing generation roots, and the first attempted cluster push failed as explained above.

## Next actions

1. Fix the likely Julia command-literal parse issue in the LOCAL development worktree. Review narrow harness risks, commit, and create a NEW immutable source tag. Preserve v1/v3 artifacts.
2. Deploy a clean FLOWPanel worktree for the corrected pin (unchanged dependency v1 worktrees may be reused), generate dedicated env/pins, annotate actual final execution commit and publish it through local git. Keep current ORC policy symlinks.
3. Refresh availability and validate exact submission with `sbatch --test-only`. Use non-preemptible exclusive Zen3 and 500 GiB. Site policy says leave partition unspecified unless required; the successful prior request used `--qos=normal -C zen3 --exclusive --mem=500G --nodes=1 --ntasks=1 --cpus-per-task=64 --time=04:00:00`, eligible on m12. Idle private Zen3 capacity was standby-only and was rejected for this resumeless driver. Availability/ETAs are stale now.
4. Consider a <=1h `--qos=test` allocation for `COLD_PILOT_STAGE=controls_smoke`, after its exact request passes test-only. It must be non-preemptible. Validate syntax first, then isolated precompile, sequential controls at Julia/BLAS 1/1 and 4/1, then both R1 seed solvers at 4/1. Resolve each failure before continuing.
5. Preserve smoke `selected.toml` (calibrate generated settings once). Submit `COLD_PILOT_STAGE=timing_profiles` with its absolute `CONFIG_FILE`, using smoke runtime to size walltime. Within ONE allocation, run timing processes sequentially at 4/1, 64/1, 64/64, then separate solver profile processes at 64/1. Selected settings/tolerances must never change. All timing arms on same node.
6. Runtime acceptance: convergence, authoritative BC relative L2 <=1e-6, certified FMM/direct disagreement <=1e-7, fresh/prepared agreement <=1e-8, memory limits, requested/observed threads, loaded provenance exclusively pinned worktrees. A thread-arm convergence failure is a failure to investigate, never permission to retune.
7. Harvest all samples/min/median/spread, prepared versus constructor+first solve, validation records and CPU/allocation profiles. Update the review report with final pins/jobs/paths/tests/timings/profile observations/remaining gaps. No broad screening, R2/R3, optimization, scaling plots or full ILU fairness claims.

### Submission environment and output layout

Submit from the pinned FLOWPanel worktree. Export:

- `COLD_PROJECT=/home/rander39/campaigns/<generation>/env`
- `CAMPAIGN_PINS=/home/rander39/campaigns/<generation>/pins.toml`
- `COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`
- `COLD_PILOT_STAGE=controls_smoke` or `timing_profiles`; latter also needs CONFIG_FILE.

Override Slurm output/error to `$COLD_DATA_ROOT/slurm-%j.{out,err}` so logs reside outside source worktrees. Driver claims a unique `pilot-$SLURM_JOB_ID` generation, with distinct OUTDIR/BENCH_CASE_ROOT per Julia process. Timing/profile stage reuses the selected config input, not previous output directories. Shell driver assigns an explicit affinity mask of 64 distinct physical cores.

Do not write a notebook entry without approval. After validated pilot results, offer a concise entry and ask approval/desired detail per global notebook policy. No notebook entry has been written.
