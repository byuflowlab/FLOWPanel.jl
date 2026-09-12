# BRAINSTORM 021: R4 diagnostics handoff

Prepared 2026-09-12. This is the current handoff, superseding the next-action
sequence in `fgs_opt_v8_reset_prompt_20260912.md`.

## User intent and stopping boundary

Ryan explicitly redirected this campaign: skip R3 and remaining R2 tuning,
and go straight to R4 to expose implementation costs at larger problem size.
He accepts carrying R2 settings forward only as provisional starting points
and expects them to need retuning. His priority is code optimization informed
by R4 profiles, not merely finding the best configuration.

This next agent owns diagnostics, configuration tuning, validation, and an
evidence handoff. Run the R4 work below, then stop with validated profile
results and a ranked list of implementation opportunities. Ryan will return
that evidence to the code-optimization agent. Do not begin implementing solver
performance changes in this diagnostics phase. Harness fixes needed to obtain
valid measurements are in scope, with a new clean generation when required.

No R4 job has been submitted. No v9 deployment or execution tag exists yet.
v9 is a clean, locally committed/tagged source handoff, not a validated HPC run.

## Required reads and delegation

Read `~/.claude/CLAUDE.md`, repo `CLAUDE.md`, and `agent_policies/HPC.md`.
Read WORKFLOW/TESTING policies before code or validation work. On ORC read
`/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md` and its applicable
`ai-docs/` references before cluster operations. Follow repository delegation:
`hpc-monitor` for status/output checks, `harvester` for scraping/tabulation,
`test-runner` for controls, and `hpc-storage` for bounded storage preflight.

Governing plan: `fgs_initialized_cpu_optimization_plan_20260911.md`, including
the new **User steering, 2026-09-12** section. Old R2→R3→R4 ordering is
superseded. Original acceptance/reproducibility requirements remain binding.

## Completed R2 evidence — do not rerun R1 or R2

Job **13653852** ran to completion in 46m23s. The output `COMPLETED` marker,
all stage outputs, leaf statuses, and gates were checked; success is not based
only on scheduler accounting. Remote output root:

`/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/opt-13653852/`

Durable local evidence:

`BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_opt_evidence_20260912/opt-13653852/`

Start with `harvest_summary.md` and `screen_rank.csv`. Configuration mappings
were verified against requested/config TOMLs, not directory creation order.
All five candidates passed. Prepared j64/b1 results (five trials each):

| inner | median seconds | max-minus-min seconds | outer iterations | estimated sweeps | estimated FMM passes |
|---:|---:|---:|---:|---:|---:|
| 3 | 1.39697637 | 0.026764171 | 18 | 54 | 19 |
| 5 | 1.44921144 | 0.049635665 | 12 | 60 | 13 |
| 10 | 1.52587803 | 0.017148732 | 7 | 70 | 8 |
| 2 | 1.58072510 | 0.030116398 | 27 | 54 | 28 |
| 1 | 1.99868633 | 0.078553184 | 51 | 51 | 52 |

Worst recorded screen certified BC rel-L2: 7.78557547e-7; evaluator delta:
5.29506755e-9; repeated-solution delta: 0. Seed prepared baseline medians:
j4/b1 1.82464498 s; j64/b1 1.51121875 s. These five-trial screen observations
are not final general speedup claims.

R2 profile text is under
`profile-fgs-j64-b1/R2/fgs_4ea88e39ddb33b21/j64_b1/` in that evidence root.
The accumulated ten-solve main-task profile has 3,857 snapshots: 2,697 at
`compute_nonself_products!`, 2,715 at `gemv!`, 296 at `solve_leaf!`, and 2,712
at BLAS `dgemv_64_`. These overlapping stack counts identify dense nonself
GEMV as the R2 priority; they are not additive wall-time shares. Read the R4
profile independently rather than assuming this attribution carries over.

## Source and dependency pins

Local v9 worktree: `/private/tmp/flowpanel-cold-opt-20260912-v9`

- Branch: `cold-opt-20260912-v9`
- Annotated source tag: `campaign/p021-cold-source-20260912-v9`
- Commit: `721235e86f8db3ecbe6d3c6c00514b0de5172baf`
- Worktree clean at handoff. Branch and tag are **local only, not pushed**.
- Parent: v8 source `campaign/p021-cold-source-20260911-v8` (`18f97b5`).

Never use the Dropbox checkout's stale untracked `benchmark/fgs_cold_*`
drafts for execution or deployment. Do not edit the v8 execution worktree.

Existing v8 deployment/provenance to use as a template:
`fgs_opt_v8_provenance_20260911.md`, and on ORC
`/home/rander39/campaigns/p021-cold-opt-20260911-v8/{env,pins.toml}`.

Reuse the unchanged dependency execution worktrees and annotated tags:

| Package | ORC path | Tag | SHA |
|---|---|---|---|
| FastMultipole | `/home/rander39/campaigns/p021-cold-20260910-v1/FastMultipole` | `campaign/p021-cold-exec-20260910-v1` | `ef10643a401d6da16e28be87805b67d11bdf1fb5` |
| FLOWVPM | `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl` | `campaign/p021-cold-exec-20260910-v1` | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` |

Deploy through the established git workflow: publish the v9 source tag/branch,
fetch on ORC, then use `scripts/prep_campaign_worktree.sh` from the cluster
clone to create a new v9 worktree. Proposed root:
`/home/rander39/campaigns/p021-cold-opt-20260912-v9/`.
The preparation script's data-symlink commit needs its own annotated
`campaign/p021-cold-exec-20260912-v9` tag, published to origin. Check exact
targets and clean state first; never replace an existing worktree or tag.
Create a separate environment with FLOWPanel's Manifest dev-path repointed
to v9 and dependency dev-paths still at the pins above. Record all three
actual execution SHAs, tags, paths, environment hashes, and launch settings
in a new provenance file **before submission**. Do not modify the v8 env.

## What v9 changes

Only four harness/test/docs files changed; no solver implementation changed:
`benchmark/fgs_cold_common.jl`, `benchmark/run_cold_opt.slurm.sh`,
`benchmark/fgs_cold_README.md`, `test/runtests_benchmark_cold.jl`.

1. R4 support uses the existing frozen `65_209` mesh, **58,192 panels**, from
   `phase1_case.jl`. Provisional FGS seed: P8/MAC0.4/leaf100/inner3,
   lexicographic, cached leaf LU, rlx=1, max_iterations=300. It derives from
   the R2 winner; it is not an established R4 optimum. The cold FGS constructor
   uses its explicit configuration, not the shared phase case's `TUNED`
   lookup (which has no hardcoded R4 entry). Do not launch unrelated tuning
   machinery just to silence that unused-table warning.
2. R4 follows the plan's certified evaluator with direct fallback when
   certification is inconclusive. R4 smoke no longer forces direct evaluation
   on every solve. All numerical thresholds and prior R1–R3 behavior remain
   unchanged. Direct-only fallback must be reported as such; NaN disagreement
   when no direct crosscheck was performed is not a measured zero.
3. `SCREEN_BASE_FILE=/absolute/bases.toml` reads one or more saved selected
   configurations. It requires screen stage and explicit `SCREEN_SET`,
   conflicts with `CONFIG_FILE`, includes each base and its one-factor
   neighbors, deduplicates them, and resets FGS tolerance to 0 so every point
   is recalibrated. The file is copied and hashed in process provenance.
   Launcher spelling: `COLD_OPT_SCREEN_BASE_FILE`. v8 could only screen from
   the fixed rung seed, which was insufficient for tuning around winners.
4. `COLD_MIN_REPS=10` raises the minimum unprofiled repetitions to ten.
   Default 1 preserves the old adaptive count. This is distinct from K_REPS.
5. `COLD_OPT_PROFILE_REPS=20` requests twenty accumulated prepared CPU solves
   (launcher default remains ten). Allocation profile is still one solve.
6. `COLD_OPT_STAGE=attribution` runs parse, precompile, j1/j4 controls, smoke,
   j4/b1 and j64/b1 prepared baselines, and profiling, skipping the screen.
   `all`, `controls_smoke`, and `screen_profile` retain their previous roles.

Local checks passed: all nine driver/control files parse; stdlib-only saved
base roster, tolerance reset, input preservation, R4 acceptance and preflight
checks; shell syntax; `git diff --check`. Full controls did **not** execute
locally because the available Julia environment lacks Krylov (initial sandbox
depot pidfile restriction was also observed). This is not an HPC failure and
not evidence that full controls passed. Run the full suite at pinned Julia
1.11.7 j1/b1 and j4/b1 on a compute allocation before measurements. The
launcher already enforces that ordering. Source fixes after this tag need a
new generation; do not move the v9 tag or edit a queued execution tree.

## Diagnostics sequence for the next agent

1. Read policies; confirm existing SSH ControlMaster, deployment pins, mesh,
   environment, output paths, disk headroom, and fresh allocation availability.
   In this Codex environment sandbox SSH initially failed to access the
   control socket; approved `require_escalated` SSH worked. Do not mistake
   sandbox denial for failed authentication or retry into MFA.
2. Establish the R4 control and profile, and screen inner `{1,2,3,5,10}`.
   The existing launcher can aggregate all of this into one job. After v9
   deployment, its relevant environment is:

   ```bash
   export COLD_PROJECT=/home/rander39/campaigns/p021-cold-opt-20260912-v9/env
   export CAMPAIGN_PINS=/home/rander39/campaigns/p021-cold-opt-20260912-v9/pins.toml
   export COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910
   export COLD_OPT_RUNG=R4 COLD_OPT_STAGE=all
   export COLD_OPT_SCREEN_SET=inner:1,2,5,10
   export COLD_MIN_REPS=10 COLD_OPT_PROFILE_REPS=20
   unset CONFIG_FILE SCREEN_SET SCREEN_BASE_FILE COLD_OPT_SCREEN_BASE_FILE
   # Submit benchmark/run_cold_opt.slurm.sh from the clean v9 FLOWPanel root
   # only after the execution pins/provenance and resource choice are recorded.
   ```

   **Include inner=10 explicitly**: v9's R4 seed is inner=3, so the old
   default roster `inner:1,2,3,5` would omit 10. For a baseline/profile first
   without a screen, use `COLD_OPT_STAGE=attribution` instead. The launcher
   profiles the calibrated initial seed, NOT the screen's winner.
3. Use the pinned Julia 1.11.7/module setup. Initial matched hardware is
   64 physical CPUs, BLAS=1, exclusive zen3, 500 GiB ceiling, normal QOS;
   launcher default walltime is six hours. Refresh actual availability and
   assess walltime before submitting. R4 fixture RHS is direct and rebuilt
   per process; construction/calibration/validation are outside solve timing
   but can dominate job runtime. Do not infer runtime from solve medians.
   No preemptible/requeue queue: the launcher is fail-fast, not resumable.
4. Monitor by `COMPLETED`, process logs, per-leaf status, and gate CSVs, not
   sacct alone. Failed inner candidates are findings; retain their outputs.
   Diagnose the first error, never loosen tolerance acceptance gates or
   solver convergence requirements to force a candidate to win.
5. Harvest and rank **median prepared total time to accepted accuracy**,
   not minimum time and not inner sweep count alone. `selected.toml` uses
   minimum time internally; independently choose the two median winners and
   prepare their saved configs. Record iterations, estimated sweeps/FMM
   passes, calibration tolerances, spread, allocations, retained memory/RSS,
   thread/CPU/pin provenance, and all acceptance metrics.
6. Read accumulated CPU flat/tree and allocation profiles on R4. Check sample
   counts, warnings/buffer saturation, thread/task coverage, and final profile
   validation. Do not count overlapping frames as independent percentages.
   Deserialize `.jls` only on a compute allocation. Profile the retained
   inner configuration as well as the initial control if they differ, using
   the existing profile driver with its calibrated CONFIG_FILE in a separate
   clean process/allocation. Keep profiling out of uninstrumented timing.
7. Do a bounded near/far tuning pass where useful to separate configuration
   effects from implementation costs: saved top-two bases plus leaf
   `{25,50,100,200}`, then P/MAC neighbors around retained candidates. Use
   `COLD_OPT_SCREEN_BASE_FILE` with `COLD_OPT_SCREEN_SET`, retaining two by
   median at each stage. Do not run a full Cartesian sweep or presume the
   inherited R2 P/MAC/leaf values are right. Existing colored-vs-lexicographic
   and Julia `{1,4,16,32,64}`/BLAS=1 diagnostics are appropriate if the R4
   profile motivates them; implement no new solver algorithm in this phase.
8. Confirm finalists with >=10 unprofiled trials in matched CPU conditions,
   alternating baseline/candidate batches per the governing plan. The stock
   launcher's initial baselines alone do not provide that alternating
   comparison. Repeat ambiguous results. Reprofile the retained configuration
   so the code-optimization agent sees remaining costs rather than obsolete
   settings. Configuration-only work may reuse clean v9; changed executable
   harness code requires a newly pinned generation.

## Gates and final deliverable

Require solver convergence, finite accepted solutions, BC rel-L2 <=1e-6,
certified evaluation (or the explicit R4 direct fallback), direct/FMM
disagreement <=1e-7 when both certified FMM and direct are evaluated, and
repeat solution agreement <=1e-8. Record which evaluator was authoritative;
`fmm_rel_max` is not evaluator disagreement. Preserve BLAS thread-drift checks,
memory limits, zero-reset semantics, and constructor-free prepared scope.

Save a new durable R4 evidence directory, reproducible source/exec/dependency
pins and commands, configuration/timing tables, text profiles, and a concise
handoff for the code-optimization agent. Rank plausible >=5% total solve
opportunities with supporting R4 stacks/allocations: nonself products and
matrix layout/memory traffic, scatter, leaf solves, color/scheduling overhead,
FMM/residual/scratch costs, or anything newly exposed. Distinguish measured
evidence from hypotheses and identify extra measurements needed.

Review methods/results/conclusions for consistency before reporting. Stop
after handing back the diagnostics; do not claim an R4 code speedup yet.
No notebook entry has been written; offer one only after validated results
and obtain Ryan's approval and desired detail before writing.

Storage preflight at this handoff was bounded: home filesystem `df` proxy
196,156 MiB (~191.6 GiB), not an authoritative global `du` cap measurement;
R2 output 84 MiB and no VTK found. No files were archived or deleted. Refresh
as needed rather than treating this snapshot as permanent capacity.
