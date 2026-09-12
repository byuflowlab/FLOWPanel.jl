# Continue the launched R1 CPU timing/profile pilot

Prepared 2026-09-11 evening America/Boise (2026-09-12 UTC). Ryan requested a context reset **once the full profile job was launched**. That job is now submitted: **13653115**. Resume monitoring/verification and interpretation, not implementation from the old shared-checkout draft.

## Scope and required reads

- Read `/Users/ryan/.claude/CLAUDE.md`, repo `CLAUDE.md`, and `agent_policies/{WORKFLOW,TESTING,HPC}.md` before work. Read current remote `/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md` and applicable `ai-docs/` policies before cluster actions.
- Read this handoff first, then `fgs_cold_validation_pilot_plan_20260910.md` for original user intent. `fgs_cold_resume_status_20260911.md` records the intervening attempts; `fgs_cold_resume_prompt_20260911.md` is the older handoff. This file supersedes their pending/validation status.
- CPU only. CUDA toolkit is loaded **solely for precompiling extensions in the existing environment**; Ryan explicitly allowed that. No GPU request, GPU execution, CUDA speedup exploration, solver optimization, R2/R3 runs, or broad screening is authorized.
- All Julia execution, including parsing/tests/precompilation, stays on HPC compute nodes. Another agent uses the local machine. No local Julia has run.
- Preserve all existing worktrees, environments, and results. Do not edit the executing v6 worktree while its job is queued/running.
- Ryan explicitly approved publishing corrected source and execution tags to `byuflowlab/FLOWPanel.jl` on GitHub. All v4/v5/v6 source and actual execution tags are published. Do not request that permission again.
- No notebook entry has been written. Ask approval and desired detail before writing one, after validated pilot results.

## Current source and immutable execution pins

The shared Dropbox checkout is dirty and contains the **old harness draft**. It is only the home of the BRAINSTORM handoff/review documents. Do not run or deploy its benchmark files.

Latest local implementation: `/tmp/flowpanel-cold-20260910`, branch `cold-harness-20260910`, clean HEAD **`3acf1c4cff71fdee3dfae95e42d10e21012aea7e`**.

| Package | Actual execution worktree | Annotated execution tag | SHA |
|---|---|---|---|
| FLOWPanel | `/home/rander39/campaigns/p021-cold-20260911-v6/FLOWPanel.jl` | `campaign/p021-cold-exec-20260911-v6` | `08a0247d173ddd6dcfd8a16b418002de31e3ca5b` |
| FastMultipole | `/home/rander39/campaigns/p021-cold-20260910-v1/FastMultipole` | `campaign/p021-cold-exec-20260910-v1` | `ef10643a401d6da16e28be87805b67d11bdf1fb5` |
| FLOWVPM | `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl` | `campaign/p021-cold-exec-20260910-v1` | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` |

FLOWPanel source tag: `campaign/p021-cold-source-20260911-v6` (`3acf1c4cff71fdee3dfae95e42d10e21012aea7e`). Preparation adds data and current ORC-policy symlink commits; cite the actual execution SHA above for results.

Dedicated environment and pins:

- `/home/rander39/campaigns/p021-cold-20260911-v6/env`
- `/home/rander39/campaigns/p021-cold-20260911-v6/pins.toml`
- Manifest has exactly three dev paths, all pointing to the pinned worktrees above.
- Clean trees, annotated tag/HEAD equality, and smoke acceptance were rechecked immediately before submitting 13653115.

## Fixes and completed validation

New development commits after the old handoff:

1. `56411b9`: fixed the Julia command-literal parser defect by constructing `revision = tag * "^{commit}"` and interpolating `$revision`; added dependency-free `benchmark/cold_parse.jl` before precompilation.
2. `5355045`: added `cuda/12.8.1-zkkfiog` alongside `julia/1.11.7-6bmogfl` in the launcher. Full environment precompile needs `ptxas` for FLOWVPMCUDAExt even though the pilot is CPU-only. It saves `modules.txt` and `ptxas_path.txt`.
3. `3acf1c4`: changed the candidate-axis collection from an array to a tuple. Julia promoted the numeric vectors in the Krylov axis array, turning integer/Boolean settings into floats; strict configuration validation caught this. Tuples preserve field types. No validator was relaxed.

Historical attempts, all preserved:

- 13642618 (v1): failed before Julia due to unset HISTCONTROL in `/etc/profile`; fixed in earlier `3194a2c`.
- 13644388 (v4): parse passed; precompile failed because CUDA discovery lacked `ptxas`.
- 13644436 (v5): precompile passed; 1-thread BLAS checks passed; configuration tests caught the type-promotion defect.
- **13644535 (v6): COMPLETED, exit 0:0, elapsed 00:07:44 on m12-3-14.** Syntax, precompilation, both control runs, and both R1 seed smoke tests passed.

Each control process passed **530/530** tests at Julia/BLAS **1/1** and **4/1**: BLAS 12, configuration/immutability 424, invalid-input filesystem effects 60, evaluator semantics 14, injected failure propagation 20. The injected error tracebacks are expected tests, not failures.

Smoke evidence root:
`/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/pilot-13644535`

| R1 solver | Trials | Iterations | Certified BC relative L2 | Direct BC relative L2 | FMM/direct disagreement | Max repeat/independent solution delta |
|---|---:|---:|---:|---:|---:|---:|
| FGS | 3/3 accepted | 6 | 9.53567919e-8 | 9.57141435e-8 | 3.92855040e-9 | 0 |
| ILU-GMRES | 3/3 accepted | 7 | 8.57921900e-7 | 8.57157852e-7 | 3.92854753e-9 | 0 |

All smoke rows report solved, certified, independent constructors, zero reset, accepted, and eligible. Max recorded process-lifetime RSS was 1,603,387,392 bytes; retained body+solver bytes were 211,972,080 (FGS) and 115,103,934 (ILU). Scheduler MaxRSS differs due to sampling; do not conflate it with process lifetime RSS or retained size. Provenance records Julia 1.11.7, observed Julia/BLAS 4/1, the clean package pins, and 64 distinct physical cores in the affinity mask.

These are correctness smoke results, **not the final timing comparison**. Smoke included compilation and unmeasured constructor work. A bounded source review found no additional definite mismatch in constructor/reset fields or FGS callback APIs.

## Frozen configuration

Use the existing absolute file, without recalibration:

`/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/pilot-13644535/smoke-j4-b1/selected.toml`

SHA-256: **`c596a5f895369f04fc1cfd9c8e5b40a3169b2a5d7ff1607ab98c2f3bb9cdf064`**.

- FGS: P=6, MAC=0.3, leaf=50, inner=10, max_iterations=300, tolerance=1.7658196194548e-8, rlx=1, lexicographic, cache_leaf_lu=true.
- ILU-GMRES: P=17, MAC=0.65, leaf=6, memory=50, itmax=500, atol=1e-14, rtol=1e-6, ilu_leaf=10, ilu_MAC=1, pattern_entries_per_panel=8192, equilibrate=false, diagonal_shift=0, cache_tree=true, persistent_plan=true, cache_nearfield=false.

## Full CPU timing/profile job: 13653115

Submitted at about **2026-09-12 02:10 UTC** (2026-09-11 20:10 America/Boise), after exact `sbatch --test-only` acceptance. The test-only identifier 13653114 is NOT a separate submitted job.

- Script: v6 `benchmark/run_cold_pilot.slurm.sh`, submitted from the pinned FLOWPanel worktree.
- Stage: `COLD_PILOT_STAGE=timing_profiles`.
- Environment: `COLD_PROJECT`, `CAMPAIGN_PINS` set to v6 paths above; `CONFIG_FILE` is the frozen smoke file above.
- `COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`.
- Request: non-preemptible `--qos=test --time=01:00:00`, exclusive Zen3, 1 node, 1 task, 64 requested CPUs, 500G memory, no partition specified, no GPUs. Exclusive-node scheduling reports 128 allocated processors; the driver explicitly binds to 64 distinct physical cores.
- Walltime sizing: smoke/controls/precompile took 7m44s. Full sequential timings/profiles were estimated comfortably below one hour; profile overhead and high-thread behavior remain unverified. Do not claim a runtime guarantee.
- Order: timing processes at Julia/BLAS **4/1**, **64/1**, **64/64**, then separate **FGS** and **ILU-GMRES CPU/allocation profiling processes at 64/1**. All within this single allocation. Fail-fast: a failing arm prevents later arms.

Logs and output:

```text
/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/slurm-13653115.out
/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/slurm-13653115.err
/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/pilot-13653115/
  timing-j4-b1.log
  timing-j64-b1.log
  timing-j64-b64.log
  profile-fgs-j64-b1.log
  profile-krylov_ilu-j64-b1.log
  selected.sha256
  COMPLETED                 # only if all stages finish successfully
```

Each process has its own same-named output directory and `fixture-<generation>` directory. Timing/profile result leaves are `R1/<kind>_<config-id>/j<Julia>_b<BLAS>/`. Inspect `status.toml` and validation CSVs, not stdout alone.

## Next actions for the new agent

1. Refresh **only job 13653115**, using `/apps/slurm/latest/bin/squeue` and `sacct`; inspect logs and completion/status files. The handoff is intentionally before completion. Do not submit a duplicate because it is pending.
2. If any arm fails, diagnose the first error. A thread-arm convergence failure is a finding to investigate, never permission to retune the frozen configuration. Preserve partial results. Any needed changes require a new clean annotated campaign pin/worktree; never mutate a running one.
3. If complete, verify every requested arm and both profiles actually ran and passed: convergence, authoritative BC relative L2 <=1e-6, certified FMM/direct disagreement <=1e-7, fresh/prepared agreement <=1e-8, memory, observed/requested threads, pinned loaded provenance, unchanged selected hash/config values.
4. Harvest all `trials.csv` and `summary.csv` samples (min/median/max/spread, prepared versus constructor+first solve), warmup/compile/fresh validation, convergence records, and separate unprofiled trial/CPU/allocation profiles. Repetition policy is 5 below 60s, 3 below 600s, otherwise 2. Compare only comparable timing scopes.
5. Review CPU/allocation profiles for measured bottlenecks, distinguish measured evidence from inference, and report findings. Do not implement optimization without further user direction. No full ILU fairness claim, measured component-work attribution, scaling plots, or broad winner claim is supported by this pilot alone.
6. Update `fgs_cold_validation_review_20260910.md` with final outcomes, quantitative tables, pins, job/node/thread/config details, and limitations. Offer a notebook entry after validated results; obtain approval before writing.

## Operational details

- SSH: `ssh -o BatchMode=yes -o ConnectTimeout=8 orc`, via `exec_command(sandbox_permissions="require_escalated")` when needed. After a long pause the session expired; Ryan reauthenticated and it works again. If actual authentication/MFA fails, stop retries and ask Ryan to reconnect.
- Remote login PATH may omit Slurm and `module`; absolute Slurm binaries work. Module commands need `/etc/profile`, which can take tens of seconds. Under `set -u`, disable nounset while sourcing it. Never infer a Slurm outage from PATH/startup/sandbox trouble.
- Scheduler periodic checks must be >=60 seconds apart; use backoff and inspect only our jobs. No monitor daemon is left running by the previous agent.
- System Python on login nodes lacks `tomllib`; use compatible stdlib parsing for tiny known files or a site Python module. `rg` is also absent on that login PATH; use a suitable fallback there. Do not run Julia to work around this on login nodes.
- Prefer bounded lower-cost agents for monitoring/harvesting. A prior monitoring agent hit a usage limit during the session gap; small one-shot checks were then handled inline. Do not wait indefinitely on a failed agent.
- A small local smoke-evidence copy is at `/tmp/p021-cold-smoke-13644535-evidence`; authoritative full results remain on ORC.
- Preparation scripts `/tmp/prep_cold_v{3,4,5,6}.sh` are historical helpers, not instructions to recreate existing generations. Existing-root collisions must not be bypassed.

## Suggested next-agent instruction

Continue the R1 CPU timing/profile pilot using this handoff. Read it before taking action. Job 13653115 has been submitted from the clean v6 pin after syntax, 530/530 controls at both thread settings, and both R1 smoke solvers passed. Monitor and verify the job, harvest timings and CPU/allocation profiles, and update the review. Keep all Julia on HPC; keep settings frozen; preserve every worktree/result; no CUDA execution or optimization implementation. Do not use the shared checkout's old harness draft.

### Final launch observation

At the handoff check, job **13653115 was RUNNING on m12-3-5, elapsed 3m30s**. `pilot-13653115/` exists with `selected.sha256`, module/affinity/pins records, and the first `timing-j4-b1` output/log and fixture directory. Slurm error log had no output. Later timing arms/profiles were not yet verified. The local smoke-evidence copy completed successfully. No intentional monitor loop or background preparation command remains; the Slurm job continues independently.
