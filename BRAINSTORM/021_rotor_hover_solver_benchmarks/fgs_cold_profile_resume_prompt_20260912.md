# Resume R1 CPU profile continuation after timing harvest

Prepared 2026-09-12 02:59 UTC (2026-09-11 evening America/Boise). Ryan requested a context reset while the repaired profiles-only job **13653450** was running. Continue monitoring and profile interpretation; do not repeat the accepted timing arms.

## Read first / scope

Read `/Users/ryan/.claude/CLAUDE.md`, repository `CLAUDE.md`, and `agent_policies/{WORKFLOW,TESTING,HPC}.md`. Read current remote `/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md` and applicable `ai-docs/` before cluster actions. The original intent is in `fgs_cold_validation_pilot_plan_20260910.md`; the previous full handoff is `fgs_cold_profile_resume_prompt_20260911.md`. This file supersedes their pending job/status instructions.

- CPU only. All Julia stays on HPC compute nodes, including parsing, tests and precompilation. Another agent uses the laptop.
- Keep settings frozen. No calibration, optimization, solver changes, R2/R3, broad screening, GPU execution or scaling plots.
- Preserve every existing worktree, environment and result. Never edit the running v7 worktree or deploy the old dirty shared-checkout harness.
- Ryan already approved publishing corrected source and execution tags to **byuflowlab/FLOWPanel.jl**. All v7 tags below are published; do not ask again. Auto-review initially rejected an unverified `origin` push, but accepted after explicit destination verification and evidence of that existing approval. Nothing remains blocked.
- No notebook entry written. Offer only after validated results; obtain approval and desired detail before writing.
- Use bounded lower-cost monitor/harvester agents per repository policy; main agent owns edits and conclusions.

## Verified timing job 13653115: partial pilot, accepted timing arms

Job **13653115 FAILED, exit 1:0**, on **m12-3-5** after **19m58s**. Start/end: 2026-09-12 **02:10:41–02:30:39 UTC** (2026-09-11 20:10:41–20:30:39 America/Boise). It completed all three timing processes, then failed before FGS profiling produced any profile data. ILU profiling never started; no `COMPLETED` marker. Empty Slurm stdout/stderr are expected because each process redirects its log.

All **60 timing samples** passed: 3 Julia/BLAS arms × 2 solvers × 2 scopes × 5 samples. All six leaf status files say completed; configs equal requested configs; selected files match the original hash. Min/median/max/spread arithmetic, repetition policy, accuracy, memory, thread controls and loaded clean pins were harvested and independently reviewed.

| Julia/BLAS | FGS prepared median s | FGS fresh median s | ILU prepared median s | ILU fresh median s |
|---|---:|---:|---:|---:|
| 4/1 | 0.612503448 | 19.0729114 | 11.6191055 | 15.5452754 |
| 64/1 | 0.578412639 | 18.9046045 | 1.36147223 | 4.62458523 |
| 64/64 | 0.875811460 | 19.9931212 | 1.45608379 | 4.64427820 |

Prepared is compiled, constructed/cache-ready isolated zero-reset `_solve!`. Fresh is compiled constructor + first `_solve!`; it excludes process startup, fixtures/RHS, reset and BC diagnostics. Do not compare different scopes or call compile-validation a benchmark sample. At 64/1 FGS prepared is faster, but its constructor median is 18.321465 s versus ILU 3.171055 s, reversing the fresh comparison. This is not a general solver winner claim.

| Solver | Iterations | Certified authoritative BC rel L2 | Max direct BC rel L2 | Max FMM/direct disagreement | Timed retained body+solver bytes |
|---|---:|---:|---:|---:|---:|
| FGS | 6 | 9.53567919e-8 | 9.57141436e-8 | 3.92855040e-9 | 211972080 |
| ILU-GMRES | 7 | 8.57921900e-7 | 8.57157852e-7 | 3.92854753e-9 | 115103934 |

All measured fresh/prepared deltas zero. Independent warmup/convergence crosschecks populated direct/evaluator fields; timed-trial direct fields intentionally NaN. `fmm_rel_max` is infinity-norm residual, NOT FMM/direct disagreement (`evaluator_delta`). Each leaf has one compile validation, one warmup, one fresh warmup, ten trials and seven independently recorded convergence rows. FGS trace has seven callbacks even though timed solver iterations are six; do not reinterpret them as identical counts. Terminal internal residuals are FGS 5.11660132e-9 and ILU 8.53503058e-7; compare authoritative BC separately.

Maximum lifetime RSS across timing validations: 1,958,821,888 bytes, versus 536,870,912,000-byte ceiling. Retained ILU reaches 115,104,222 bytes with convergence history. Allocation bytes, retained memory, process-lifetime RSS and scheduler MaxRSS are different metrics. Initial accounting did not provide MaxRSS.

Timing hardware: AMD EPYC 7763 Zen3, Julia 1.11.7; 64 distinct physical-core affinity; requested 64 CPUs/exclusive, allocated 128 logical processors; 500G, non-preemptible test QOS, one-hour limit, no specified partition/GPU request. BLAS controls/GEMM assertions match 4/1,64/1,64/64.

### Timing evidence already harvested and review updated

- Authoritative ORC: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/pilot-13653115/`.
- Complete local text/CSV/TOML/log snapshot: `/tmp/p021-cold-13653115-evidence/` (no Julia deserialization needed).
- Local scratch reports: `harvest/compact_report.md`, `harvest/detailed_acceptance.md` under that snapshot. The latter has a mislabeled gate-count summary table; prefer original CSVs and the independently checked review for counts. Its numeric timing tables are correct.
- Durable harvested CSVs: `BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_cold_pilot_evidence_20260912/`, with all 60 samples, summaries, excluded compile/warmup, convergence, status/config/provenance tables, and frozen `selected.toml`.
- Scratch Python: `/tmp/p021_harvest_13653115.py`, `/tmp/p021_detail_13653115.py`. They may help reuse structure; inspect before adapting.
- **`fgs_cold_validation_review_20260910.md` already updated**, including full min/median/max/spread tables, memory/scope discussion, failure cause and v7 continuation. Opening status and last paragraph still say profiles running; update after completion. Earlier history deliberately preserved and marked chronological.
- Independent evidence audit found no numeric correction needed. Never shorten “accepted timing evidence from failed/partial pilot” to “pilot passed.”

## Profile failure root cause and repaired pins

Raw failure: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/pilot-13653115/profile-fgs-j64-b1.log`:
`UndefVarError: Allocs not defined in Main` at `rotor_hover_solver_phase2_profile.jl:130`, top-level conditional at line25.

The legacy branch imports `Profile.Allocs` inside the conditional but calls `Allocs.@profile`; Julia expands macros in both branches before the import executes. The syntax-only gate did not exercise macro expansion. The repaired source moves the import before the conditional. A second small change adds `profiles_only` to the existing launcher, skipping controls/smoke/timing but using the identical FGS→ILU profile process loop and final frozen hash verification. Exactly two code files changed; solver, shared cold harness, constructor, timing and validators unchanged. Shell syntax/diff checks and bounded source review passed. Actual profile execution is the HPC validation.

New local development worktree: `/tmp/flowpanel-cold-profile-20260912`, branch `cold-profile-repair-20260912`, clean source HEAD **`7cda26f4d30f6107b78c618422db04536598a596`**; annotated source tag **`campaign/p021-cold-source-20260912-v7`**. Old `/tmp/flowpanel-cold-20260910` (3acf1c4) preserved.

| Package | Execution worktree | Annotated execution tag | SHA |
|---|---|---|---|
| FLOWPanel profiles v7 | `/home/rander39/campaigns/p021-cold-20260912-v7/FLOWPanel.jl` | `campaign/p021-cold-exec-20260912-v7` | `c9e33d29d9dd21b9fbd2424d578c5ce4f0ecd6b6` |
| FLOWPanel timings v6 | `/home/rander39/campaigns/p021-cold-20260911-v6/FLOWPanel.jl` | `campaign/p021-cold-exec-20260911-v6` | `08a0247d173ddd6dcfd8a16b418002de31e3ca5b` |
| FastMultipole both | `/home/rander39/campaigns/p021-cold-20260910-v1/FastMultipole` | `campaign/p021-cold-exec-20260910-v1` | `ef10643a401d6da16e28be87805b67d11bdf1fb5` |
| FLOWVPM both | `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl` | `campaign/p021-cold-exec-20260910-v1` | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` |

V7 dedicated `/home/rander39/campaigns/p021-cold-20260912-v7/{env,pins.toml}`. Manifest contains exactly those three v7/dependency worktree paths. Preparation adds data/site-policy commits; results cite actual execution SHA, not source SHA. Clean trees and expected hash checked before submission. V7 source and execution tags are published. Historical preparation helper `/tmp/prep_cold_profile_v7.sh` is NOT an instruction to recreate anything.

## Frozen configuration — do not retune

`CONFIG_FILE=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/pilot-13644535/smoke-j4-b1/selected.toml`

SHA-256 **`c596a5f895369f04fc1cfd9c8e5b40a3169b2a5d7ff1607ab98c2f3bb9cdf064`**.

- FGS: P6, MAC0.3, leaf50, inner10, max_iterations300, tolerance1.7658196194548e-8, rlx1, lexicographic, cache_leaf_lu=true.
- ILU-GMRES: P17, MAC0.65, leaf6, memory50, itmax500, atol1e-14, rtol1e-6, ilu_leaf10, ilu_MAC1, pattern_entries_per_panel8192, equilibrate=false, diagonal_shift0, cache_tree=true, persistent_plan=true, cache_nearfield=false.
- IDs: FGS409bb3aaf67bbd18; ILUd569863eedc13014.

## Live profiles-only job 13653450

Submitted 2026-09-12 about02:56 UTC after exact test-only acceptance (13653446 is test-only, NOT a second submitted job). Stage **`COLD_PILOT_STAGE=profiles_only`** from v7 launcher, `COLD_PROJECT`/`CAMPAIGN_PINS` v7 paths above, frozen `CONFIG_FILE` above, `COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`.

CPU request unchanged: exclusive Zen3, one node/task,64 requested CPUs,500G, test QOS,1h,no forced partition,no GPU. Scheduler placed it on **m12-3-5** (same node as accepted timings). CUDA12.8.1 module is loaded solely for existing environment extension compilation, allowed by Ryan; Julia1.11.7 pinned. No GPU kernels requested/executed.

Sequence: separate FGS then ILU-GMRES CPU/allocation profiles, each Julia/BLAS64/1. Fail-fast. It does not repeat timings, controls or smoke. Each `unprofiled_trial.csv` is a first solve in that profile process and may include compilation/lazy initialization; do not compare it directly with compiled prepared timing medians. CPU and allocation solves follow it; validation occurs outside recorded regions. Allocation sample rate0.01; raw sampled bytes are not total allocations.

Output and logs:
```
/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/slurm-13653450.out
/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/slurm-13653450.err
/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/pilot-13653450/
  profile-fgs-j64-b1.log
  profile-krylov_ilu-j64-b1.log
  selected.sha256
  COMPLETED
```
Process result leaves: `<process>/R1/<kind>_<config-id>/j64_b1/`, including `status.toml`, config, `unprofiled_trial.csv`, `cpu_{flat,tree}.txt`, `cpu_profile.jls`, `allocations.txt`, `allocation_profile.jls`, `cpu_validation.csv`, `allocation_validation.csv`. No profile artifact was yet harvested at handoff preparation.

## Next actions

1. Check only13653450 (squeue/sacct absolute `/apps/slurm/latest/bin/`), bounded logs/status. Never duplicate pending work. Scheduler periodic checks≥60s apart, back off; no monitor daemon intentionally remains.
2. If failure, diagnose first error before changing anything. Preserve partial outputs, no retuning. Any harness repair needs a new clean annotated campaign generation; do not edit running v7.
3. On completion verify both profile arms actually ran, statuses/COMPLETED and exit0, unprofiled plus CPU/allocation validation, BC≤1e-6, evaluator disagreement≤1e-7, solution agreement≤1e-8, memory, requested/observed64/1, loaded pins and unchanged selected hash/config.
4. Harvest profile text/CSV/TOML/logs in a NEW local snapshot and durable evidence files. Leave raw `.jls` files intact on ORC; any Julia deserialization stays on a compute allocation. Inspect CPU native frames with thread/task grouping and allocation stacks. Distinguish measured evidence from inference; no optimization implementation. CPU profiles include idle/waiting worker stacks, not just useful work; avoid summing inclusive frames as disjoint percentages.
5. Update review opening/current status and final section with profile job/node/pins, validation and measured bottlenecks, interpreting profiles alongside accepted v6 timing scopes. Profiles use a new source/execution generation with only the documented two-file harness delta; say so.
6. Review evidence/conclusions for consistency before reporting. Pilot alone does not establish full ILU fairness, measured component-work attribution, broad winner or scaling claim. Offer notebook entry only after validation and get approval before writing.

## Access / operational notes

SSH `ssh -o BatchMode=yes -o ConnectTimeout=8 orc`, using exec_command require_escalated for network/sandbox when needed. A live session worked throughout this continuation. If actual authentication/MFA fails, stop retries and ask Ryan to reconnect. Do not infer a Slurm outage from PATH, startup or sandbox failures. Remote rg and Python tomllib are unavailable on default login PATH; use bounded POSIX reads/transfers and parse TOML locally using Python. Do not use Julia on login nodes or laptop.

The local shared checkout remains dirty with unrelated changes and the old harness draft. Only the review and new evidence/handoff documents were edited there. All implementation edits are in the isolated new local worktree. All prior v1–v6 worktrees/results remain preserved. No storage cleanup, archive, deletion or notebook write was performed.

### Final reset observation (2026-09-12 about 03:00 UTC)

Final scheduler check: **13653450 RUNNING**, elapsed **3m33s**, node **m12-3-5**, no end time, no `COMPLETED` marker yet; ILU profile process directory exists and its log was still empty. No persistent monitoring loop remains.

FGS leaf already has `status = "completed"`, both validation CSVs, CPU flat/tree text (~857K/~707K), CPU raw profile (~1.4M), allocation raw profile (~3.4M), allocation text (~909K), and unprofiled trial CSV. Allocation header reports `sample_rate=0.01; sampled_bytes=165856`. These were only listed/read in place, not yet harvested or fully validated.

FGS log includes two **“no samples collected” warnings**, but CPU artifacts are substantial and nonempty. Do not assume the entire CPU profile is empty or submit a rerun based on the warning alone: inspect thread/task groups and actual samples first (empty groups may explain warnings; this is an unverified hypothesis). Profile gate values, sample coverage and bottleneck interpretation remain next-agent work.
