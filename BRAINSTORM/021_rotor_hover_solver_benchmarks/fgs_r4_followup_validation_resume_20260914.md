# BRAINSTORM 021: R4 follow-up validation resume handoff (2026-09-14)

## Latest update — 2026-09-15

**Supersedes the historical state below.** Job 13690579 failed before
measurements with a Julia world-age error. Its evidence is harvested and
hash-verified. FLOWPanel v15 fixes the driver entrypoint and adds the missing
dependency census and hardware-counter capability probe.

Replacement job **13694724** was submitted from verified v15 source and was
**PENDING (Priority)** at **2026-09-15T12:43:15Z**. Do not resubmit or modify
its silo source. Read
[`fgs_r4_followup_validation_20260915.md`](fgs_r4_followup_validation_20260915.md)
for the new tag/SHA, exact output path, tests, evidence inventory, and
remaining completion/cleanup conditions. No R4 diagnostics have yet passed.

## Current state — full campaign running

Updated `2026-09-14T23:43:40+00:00`. Ryan said **continue** after the break request. Sol hit
its usage limit and stopped; the parent agent took over. No subagent or
recurring monitor is currently supervising this job. The next session must
resume monitoring; Slurm itself continues independently of this conversation.

**Job 13690579 is RUNNING on m12-3-14**, verified at elapsed 23 seconds.
This is the full v14 campaign, exclusive zen3, 64 requested CPUs, BLAS=1,
500 GiB memory request, 12-hour limit. It starts with the required controls,
then runs j1/4/8/16/32/64. No R4 measurements have yet been validated.

FLOWPanel v14 is now deployed; both complete source manifests passed checksum
verification. The v14 pins were copied explicitly to the silo `pins.toml`.
No executable source changed during the parent continuation.

### Passed short control

Job **13690544 COMPLETED, exit 0:0, elapsed 2m22s**. The exact v14
FastMultipole test invocation passed all nine testsets, including FGS stage
instrumentation, cached LU, threaded repeatability, and colored sweeps. The
fixture/import setup was checked against the package test entry point before
submission. This resolves the missing-import failures from previous attempts.

- Local verified test log and completion marker:
  `fgs_r4_followup_evidence_20260914/fm-v14-control-13690544/`
- Remote control output:
  `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/fm-v14-control-13690544/`
- Durable control wrapper:
  `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/p021-fm-v14-control.sh`
- The harvested control directory also retains both source manifests, all
  package pins, and Project/Manifest environment copies outside the silo.

### Immediate next actions after reset

1. Do not resubmit or edit the running silo. Query `squeue -j 13690579`, or
   `sacct -j 13690579 -X -o JobID,State,Elapsed,ExitCode` if terminal. Respect
   the minimum 60 seconds between periodic scheduler queries.
2. Inspect output root
   `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/diag-v10-13690579/`.
   Scheduler stdout/stderr are currently inside the silo at
   `FLOWPanel.jl/logs/slurm/r4-diag-v10-13690579.out` and `.err`; copy these
   outside the silo before cleanup. Check control logs first, then each
   `j<N>-b1/process.log`, `results/status.toml`, and acceptance CSVs.
3. A successful campaign needs the root `COMPLETED` marker and all required
   per-arm gates. On failure, diagnose the first substantive error; preserve
   failed output and avoid repeated full-campaign retries for runner setup.
4. Remaining work includes separate available hardware counters, evidence
   harvest, instrumentation-overhead/stage-budget analysis, and the detailed
   handoff's final deliverables. Neither submission nor successful unit tests
   establishes R4 accuracy or a performance conclusion.
5. Harvest campaign results and all scheduler logs into the same durable local
   evidence family, audit gates and source/env provenance, then update the
   diagnostics report. Follow the exact silo cleanup condition below only
   after all jobs are terminal and durable evidence is verified.

The older break snapshot below is historical. Its statements that v14 has not
been deployed and no job is running are superseded by this section.

## Historical state at the earlier break

No Slurm job is running. The most recent attempt, job `13689318`, is terminal
`FAILED` (`ExitCode=1:0`, elapsed `00:08:44`) and stopped in the FastMultipole
control command before any R4 timing or profiling measurement. Do not cancel or
resume any of the failed jobs. Their outputs are intentionally retained.

Ryan explicitly authorized creation and rsync deployment of a new isolated ORC
silo, with the binding condition that **only that newly created silo is deleted
after all campaign jobs are terminal and durable evidence outside it is
verified**. The exact silo is:

`/home/rander39/FLOWPanel-p021-r4-diag-v10-silo`

Ryan's authorization overrides the repository's ordinary Git-only deployment
and no-new-silo rules for this campaign. Do not represent the rsync trees as
remote Git worktrees. Do not push the source branches/tags to GitHub. Do not
delete canonical outputs, shared dependency worktrees, existing silos, or the
target of the silo's `data` symlink.

There is no recurring monitor installed because no job remains active.

## Governing scope

Execute the handoff's `2026-09-14 follow-up validation request` through its
instrumentation/diagnostics stopping boundary. Do not rerun completed v9 R4
tuning or finalist confirmations. Retain inner=3, P8/MAC0.4/leaf100,
lexicographic, cached LU, zero reset, constructor-free prepared timing, BLAS=1,
and all certification gates. Required outputs remain: execution/thread/affinity
audit; exclusive stage timers and overhead controls; actual GEMV census;
j1/4/8/16/32/64 scaling on matched exclusive zen3; separate available
bandwidth/cache diagnostics; and a justified conditional follow-up. Do not
bundle a slice fix, mixed precision, or new sweep algorithm.

## Source generations and local worktrees

FastMultipole instrumentation worktree:

- path: `/private/tmp/fastmultipole-p021-r4-diag-v10`
- branch: `p021-r4-diagnostics-v10`
- commit: `87cbc8460b51f24ddf34cc5f41a1d1b6682bf04a`
- annotated tag: `campaign/p021-r4-diag-source-20260914-v10`
- parent campaign pin: `ef10643a401d6da16e28be87805b67d11bdf1fb5`
- content manifest deployed in the silo:
  `p021-r4-diag-fastmultipole-v10.sha256`, SHA256
  `d65252a830f91fe85dacaedf35ccdd0ea86a55eec222376e49906a89582c8fe6`

FLOWPanel harness worktree:

- path: `/private/tmp/flowpanel-cold-opt-20260914-v10`
- branch: `cold-opt-20260914-v10`
- latest local generation: **v14**
- commit: `7c504dfde9ad833c62bc06e3d0cd11e6ce2d3575`
- annotated tag: `campaign/p021-cold-source-20260914-v14`
- parent v9 source: `721235e86f8db3ecbe6d3c6c00514b0de5172baf`
- local v14 content manifest:
  `/private/tmp/p021-r4-diag-flowpanel-v14.sha256`, SHA256
  `80f40cda69575fd7cbe294a96b569f2c77d3a14d448a86a3ec843478287f0a1e`
- local pins file, already updated to v14 but **not yet deployed**:
  `/private/tmp/p021-r4-diag-v11-pins.toml`

The v14 source and pins have not been rsynced to ORC. The silo currently holds
the prior v13 FLOWPanel source/pins. Resume by rsyncing the clean v14 worktree
with `.git` and `data/` excluded, copying the v14 manifest, and copying the local
pins file explicitly to the remote name `pins.toml`. Verify the full manifest
and the canonical data symlink before submission.

Files changed in the executable generations:

- FastMultipole: `src/solve.jl`, `test/solve_test.jl`
- FLOWPanel: `src/FLOWPanel_solver.jl`, `benchmark/fgs_cold_common.jl`,
  `benchmark/fgs_r4_diagnostics.jl`,
  `benchmark/retained_r4_diagnostics.toml`,
  `benchmark/run_r4_diagnostics.slurm.sh`, `benchmark/cold_parse.jl`, and
  `benchmark/fgs_cold_README.md`

## ORC environment and dependencies

Silo paths:

- FLOWPanel source:
  `/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/FLOWPanel.jl`
- FastMultipole source:
  `/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/FastMultipole`
- environment:
  `/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/env`
- pins:
  `/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/pins.toml`
- manifests: silo root
- canonical data link:
  `FLOWPanel.jl/data -> /home/rander39/projects/FLOWPanel.jl/data`

The environment Manifest points to the silo FLOWPanel/FastMultipole trees and
the unchanged pinned FLOWVPM Git worktree:

- FLOWVPM path:
  `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl`
- commit: `05c658f7804ec5f9b68d4cb9826a9f97cfecb373`
- annotated tag: `campaign/p021-cold-exec-20260910-v1`

The isolated Project needed `Meshes` added as a direct dependency for
`runtests_unit_solver.jl`; this package was already pinned in the unchanged
Manifest. Current environment hashes after that correction:

- `env/Project.toml` md5: `cbdc6554b107d457e7aa85ec3ee82d4a`
- `env/Manifest.toml` md5: `b44d8a3c6d5111f52c25004a61fe184c`

Do not resolve or update the environment unless a real package error requires
it. Julia/module pin remains CUDA 12.8.1 plus Julia 1.11.7.

## Completed controls and failed attempts

All failures happened before measurement and are useful provenance:

| Job | Terminal state | Elapsed | First failure | Durable output |
|---:|---|---:|---|---|
| 13688396 | FAILED 1:0 | 2m24s | stale remote `pins.toml` cited the old FLOWPanel manifest hash; compilation passed | `data/p021-cold-20260910/diag-v10-13688396/` |
| 13688430 | FAILED 1:0 | 3m37s | solver unit test could not directly import `Meshes`; benchmark controls passed | `data/p021-cold-20260910/diag-v10-13688430/` |
| 13688474 | FAILED 1:0 | 11m24s | FastMultipole runner included `solve_test.jl` without `gravitational.jl` | `data/p021-cold-20260910/diag-v10-13688474/` |
| 13689318 | FAILED 1:0 | 8m44s | FastMultipole runner loaded the fixture but omitted `using Test` | `data/p021-cold-20260910/diag-v10-13689318/` |

Controls verified in job `13689318` before its final runner error:

- benchmark cold controls j1/b1: PASS
- benchmark cold controls j4/b1: PASS
- FLOWPanel `test/runtests_unit_solver.jl`: PASS, 465/465, 4m09.5s
- FLOWPanel `test/runtests_unit_fgs_history.jl`: PASS
- FastMultipole solve/coloring tests: not executed because `@testset` was
  undefined at parse time

The v14 one-line fix changes the FastMultipole control command to import both
`FastMultipole` and `Test` before including `gravitational.jl`, `solve_test.jl`,
and `fgs_coloring_test.jl`. No solver or instrumentation code changed after the
FastMultipole v10 commit.

## Original deployment sequence (deployment/submission now complete)

1. Re-read the global/repo WORKFLOW, TESTING, HPC, and ORC agent policies.
2. Confirm no job uses the silo and the local v14 worktree is clean at the pin.
3. Rsync FLOWPanel v14 into the silo with `.git` and `data/` excluded. Copy
   `/private/tmp/p021-r4-diag-flowpanel-v14.sha256` to the silo root. Copy
   `/private/tmp/p021-r4-diag-v11-pins.toml` explicitly to remote `pins.toml`.
4. Verify from the silo FLOWPanel root:
   `sha256sum --quiet -c ../p021-r4-diag-flowpanel-v14.sha256`; verify the
   FastMultipole v10 manifest similarly; verify `readlink FLOWPanel.jl/data`.
5. Refresh 64-core/500-GiB/12-hour availability and storage headroom. Prior
   snapshot: m12 had 28 fitting idle nodes and home `df` proxy was 200,083 MiB.
6. Submit `benchmark/run_r4_diagnostics.slurm.sh` from the silo with:

   - `COLD_PROJECT=/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/env`
   - `CAMPAIGN_PINS=/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/pins.toml`
   - `COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`

7. Confirm all controls, especially FastMultipole, before accepting any
   measurement. The launcher then runs j1/4/8/16/32/64 with BLAS=1, two
   ten-trial uninstrumented and two ten-trial instrumented alternating batches,
   direct/FMM/history equivalence controls, stage timers, GEMV census,
   per-thread CPU ticks, and thread/task-complete profiles.
8. Collect separate available DRAM/cache/bandwidth counters scoped as closely
   as ORC permits. If permissions/tools do not expose valid counters, retain
   the failure text and leave bandwidth saturation unresolved.
9. Harvest the terminal job into a new durable local evidence directory under
   `BRAINSTORM/021_rotor_hover_solver_benchmarks/`. Audit every gate and compare
   instrumented versus uninstrumented medians before using stage timings.
   Claim speedups only from accepted uninstrumented solves.
10. Update the diagnostics package/provenance with measured results and a
    justified conditional follow-up. Do not write a notebook entry.
11. Only after all jobs are terminal and outputs/provenance are harvested and
    independently verified outside the silo, delete exactly:
    `/home/rander39/FLOWPanel-p021-r4-diag-v10-silo`. Verify it no longer exists.
    Never follow or delete the canonical `data` symlink target.

## Cleanup gate

The silo must remain while any submitted job uses it. Once the final evidence
is durable outside it, remove only the exact silo path above, then confirm with
an existence check. The failed-job directories and the final result directories
under `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/` are durable
evidence and must not be deleted with the silo.
