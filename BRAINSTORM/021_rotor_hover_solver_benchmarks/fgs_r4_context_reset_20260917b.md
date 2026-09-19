# v21 colored-sweep A/B: context-reset handoff (2026-09-17b)

Supersedes `fgs_r4_context_reset_20260917.md` (its execution phase is DONE:
harness written, gated, tagged, deployed, submitted). **State: job 13738665
RUNNING on the cluster (~12 h wall limit, submitted 2026-09-17 ~01 UTC,
64c zen3 exclusive 500G `--qos=normal`, landed on m12).** Your task: watch it
to a terminal state, harvest, analyze the A/B, report to Ryan. Do NOT
re-implement anything.

## Keep the parent context small

Delegate bounded mechanical work to `.claude/agents/` subagents
(`hpc-monitor`, `harvester`, `hpc-storage`, `test-runner`, `code-scout`),
cheapest reliable model, ≤40-line summaries. Keep inline: physics/numerics
reasoning, conclusions, code edits, job submission, anything needing Ryan.
≥300 s sacct spacing, defensive parsing (banners pollute stdout: `-n`, then
filter for known state words; never `head -1` raw). `TZ=UTC`; `sbatch/sacct`
need `source /etc/profile` over non-interactive ssh. `ssh orc` needs a live
ControlMaster socket with BatchMode — auth failure = STOP (never MFA).

## Required reads (before acting)

`~/.claude/CLAUDE.md`, repo `CLAUDE.md` + `agent_policies/HPC.md`; cluster
work → `/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md`. Then:
1. `fgs_opt_r4_diagnostics_package_20260912.md` **§7 then §6** (measured
   profile + gates — the yardsticks live there).
2. `fgs_r4_followup_evidence_20260914/v21-deployment/submission-provenance-13738561.md`
   (v21 pins, design-as-launched, local gate results, AND the 13738561
   postmortem + resubmission record for 13738665).
3. Tail of `fgs_r4_followup_validation_20260915.md` (chronology through the
   resubmit).

## What v21 is

A/B: does the conflict-free colored leaf-sweep (`sweep_order=:colored`,
79 colors on this tree, median 16-way) beat serial lexicographic on **total
time to accepted accuracy**? Driver `benchmark/fgs_r4_colored_ab.jl`
(AB_MODE=calibrate/trials/activity; no perf anywhere), launcher
`benchmark/run_r4_colored_ab.slurm.sh`: controls → j64 colored-tolerance
calibration (separate calibration; ordering changes the iterate path) →
alternating uninstrumented A/B trials j∈{1,4,16,64}, 8 batches of 10 → 40
trials/order/arm (`ab_trials.csv`, `ab_summary.toml` per arm) → one j64
activity pair (`stage_thread_activity_*.csv`; budget attribution ONLY,
excluded from rankings).

Local gate facts (already established; do not redo): colored calibrates to
tol 5.223e-7 (lex retained 3.479e-7), **26 iterations vs 27 lex**; repeat
delta 0; local j4 signal colored 9.62 s vs lex 11.51 s (M-series — NOT
citable). Yardsticks for the real run: j64 lex uninstrumented median
10.96 s (v15); Amdahl-if-coloring-free ≈ 2.2 s; expect ≈6.4k color
barriers/solve of sync cost eating into that.

## Pins (all verified at deploy; worktrees clean at tags)

| Package | Tag | SHA | Where |
|---|---|---|---|
| FLOWPanel | `campaign/p021-r4-colored-source-20260917-v21` | `5f87a0e7` | `/home/rander39/campaigns/p021-r4-colored-20260917-v21/FLOWPanel.jl` |
| FastMultipole | `campaign/p021-r4-colored-fm-20260917-v21` | `c6185cdd` (merged tip, NOT pre-merge v11) | `.../FastMultipole` |
| FLOWVPM | `campaign/p021-cold-exec-20260910-v1` | `05c658f7` | `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl` |

Env `/home/rander39/campaigns/p021-r4-colored-20260917-v21/env` (julia
1.11.7, dev-paths at the worktrees, `pins.toml` beside it) — **patched
post-13738561 with direct `Meshes` + `StaticArrays`** (env-only fix; see
postmortem). Tags exist in the CLUSTER clones and locally; **NOT on origin**
(Ryan gate on pushing the merged branches — when he approves, push branches
AND both v21 tags).

Run output: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/colored-v21-13738665/`
(subdirs: `parse`, `calibrate/`, `j{1,4,16,64}-trials/`, `j64-activity/`,
`controls-*.log`, `COMPLETED` sentinel on success). Slurm logs:
worktree `logs/slurm/r4-colored-v21-13738665.{out,err}`.

## On completion (judge by outputs, never sacct status alone)

1. Confirm `COMPLETED` sentinel + every `results/status.toml` says completed;
   check each process.log tail.
2. Harvest to `fgs_r4_followup_evidence_20260914/colored-v21-13738665/`
   with `sha256sum` remote-vs-local verification (pattern of
   `counters-v20-13733332/`); include per-arm `provenance.toml`,
   `config*.toml`, `ab_trials.csv`, `ab_summary.toml`, `warmups.csv`,
   calibration outputs, activity CSVs, control logs.
3. Validate gates per arm per order (§6): certified-FMM authoritative, BC
   rel-L2 ≤1e-6, repeat/agreement ≤1e-8, finite, BLAS=1, iterations
   consistent within order; j1 warmups carry the direct-vs-FMM equivalence.
   Check iteration count per order is arm-invariant (colored is claimed
   deterministic; verify, don't assume).
4. Rank by uninstrumented `solve_seconds` medians per arm — total time to
   accepted accuracy. Compare j64 colored vs 10.96 s lex median and the
   2.2 s Amdahl bound; also j4 (coloring may pay earlier at low j).
   From the j64 activity pair: avg active threads over `nearfield_update`
   under :colored vs the 1.00 lex baseline (attribution only).
5. **Report A/B results to Ryan and STOP** — no follow-on experiments, no
   notebook entry without his separate approval (the diagnostics
   ladder-findings notebook entry is ALSO still owed, drafted on request).

## On failure

Harvest evidence FIRST (`colored-v21-13738665-FAILED/` + sha256), root-cause
before any rerun, reproduce locally when possible (local 1.11.8 gate env:
`<scratchpad>/env111` may be gone after reset — rebuild:
juliaup julia-1.11.8, Pkg.develop the three local checkouts, Pkg.add
Meshes StaticArrays). Source fix ⇒ new commit + new annotated tag + fresh
worktree; env-only fix ⇒ patch env, record in provenance, resubmit same
tags (13738561→13738665 is the template). Never edit deployed sources or
move tags. sacct FAILED alone is not proof of failure — check outputs.

## Ryan-pending (do not act without him)

- **Storage**: /home ~654 G of 400 G cap; hpc-storage cycle 2026-09-17
  archived 0 MB — ~369.5 GiB VTK in 15 `RECENT` runs awaits his
  `--include-recent --only` approval (list in the 2026-09-17 hpc-storage
  report / validation-log entry). Re-run `hpc-storage` before any new
  long submission.
- **Origin pushes** of merged branches + v21 tags.
- Notebook entries (v21 results AND the owed diagnostics-ladder entry).

## Memory

`project_021_solver_benchmarks.md` is current through the resubmit
(update it + MEMORY.md hook when the A/B lands).
