# R4 diagnostics: prompt for the next clear-context agent

> **2026-09-15 update:** Read
> [`fgs_r4_followup_validation_20260915.md`](fgs_r4_followup_validation_20260915.md)
> first. Job 13690579 failed with a world-age error before measurements;
> verified replacement v15 job **13694724** was PENDING (Priority) at
> 2026-09-15T12:43:15Z. The historical source generation and job below are
> superseded. Cleanup authorization and numerical requirements still apply.

Continue the already-authorized R4 follow-up validation campaign in
`/Users/ryan/Dropbox/research/projects/FLOWPanel.jl`. Resume execution and
analysis, not just planning. Do not resubmit the current job without checking
its state. Ryan requested a context reset while the campaign runs.

## Read first

1. `/Users/ryan/.claude/CLAUDE.md`, repository `AGENTS.md` and `CLAUDE.md`,
   and the required `agent_policies/WORKFLOW.md`, `TESTING.md`, and `HPC.md`.
   Read current ORC agent policies and applicable references before cluster work.
2. In `BRAINSTORM/021_rotor_hover_solver_benchmarks/`:
   - `fgs_r4_followup_validation_resume_20260914.md`: detailed current state,
     source pins, environment, paths, previous failures, and cleanup contract.
     Its newest section supersedes the historical break/deployment sections.
   - `fgs_opt_r4_diagnostics_handoff_20260912.md`: the **2026-09-14 follow-up
     validation request** specifies required tests and deliverables. The old
     diagnostics assignment below it is historical.
   - `fgs_opt_r4_diagnostics_package_20260912.md`: corrected scientific
     interpretation and retained configuration.

## Last verified state — refresh it

At `2026-09-14T23:43:40Z`, full campaign **13690579** was **RUNNING** on
`m12-3-14` (elapsed 23 seconds). It requests exclusive zen3, 64 CPUs,
500 GiB, 12 hours; BLAS=1. It runs controls before the j1/4/8/16/32/64 ladder.
This is a historical observation, not a guarantee the job is still running.
No R4 measurement results have yet been harvested or accepted.

The isolated 4-core FastMultipole control **13690544 passed** all nine
testsets and completed in 2m22s, exit 0:0. Evidence is already local under
`fgs_r4_followup_evidence_20260914/fm-v14-control-13690544/`, including the
completion marker, test log, source manifests, pins, and environment files.
Previous campaign attempts failed before measurements due to deployment or
test-runner setup issues; their outputs must be retained.

Sol stopped because of a usage limit; the parent agent deployed the already
committed v14 runner fix, verified both full source manifests, passed the
short control, and submitted 13690579. No recurring monitor or active
subagent remains. Resume monitoring yourself or delegate per repository
policy. Do not depend on previous agents surviving the reset.

## Exact paths and source state

- Temporary silo: `/home/rander39/FLOWPanel-p021-r4-diag-v10-silo`
- Campaign output:
  `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/diag-v10-13690579/`
- Scheduler logs inside silo:
  `FLOWPanel.jl/logs/slurm/r4-diag-v10-13690579.out` and `.err`
- Campaign env and provenance: silo `env/` and `pins.toml`.
- FLOWPanel v14: `7c504dfde9ad833c62bc06e3d0cd11e6ce2d3575`, tag
  `campaign/p021-cold-source-20260914-v14`;
  local worktree `/private/tmp/flowpanel-cold-opt-20260914-v10`.
- FastMultipole v10: `87cbc8460b51f24ddf34cc5f41a1d1b6682bf04a`, tag
  `campaign/p021-r4-diag-source-20260914-v10`;
  local worktree `/private/tmp/fastmultipole-p021-r4-diag-v10`.
- FLOWVPM unchanged pin: `05c658f7804ec5f9b68d4cb9826a9f97cfecb373`.

Both changed packages are deployed as verified rsync source trees, not remote
Git worktrees. Exact manifest hashes and environment paths are in the detailed
resume handoff and harvested control evidence. Inspect local Git state before
editing; unrelated changes and untracked campaign evidence are present.

## What to do

1. Check job 13690579 using `squeue` or `sacct`; allow at least 60 seconds
   between periodic scheduler queries. Use the existing `orc` SSH connection.
   Sandbox socket denial previously required approved escalated SSH; do not
   mistake it for authentication failure or retry into MFA.
2. Inspect control logs, per-arm process logs, status files, and numerical
   gates. Require both the root `COMPLETED` marker and valid per-arm results;
   scheduler completion alone is insufficient. Diagnose failures before retries.
3. Complete the handoff's remaining diagnostics, including available hardware
   counters, instrumentation overhead, exclusive stage-time reconciliation,
   block/dependency census, and thread scaling. Counter unavailability must
   be documented; do not infer bandwidth saturation from weak scaling alone.
4. Harvest results and scheduler logs outside the silo, verify numerical gates
   and provenance, then update the report with measured conclusions and a
   justified next experiment. Retained lexicographic GEMVs are serial at
   BLAS=1; main-task sample fractions are not total wall-time fractions.
5. Finish the bounded cleanup below and report the results and remaining
   limitations. No notebook entry without Ryan's separate approval.

Preserve convergence, finite solutions, certified authoritative FMM,
BC rel-L2 ≤1e-6, repeat agreement ≤1e-8, and direct/FMM disagreement ≤1e-7
when both are evaluated. Preserve zero-reset and constructor-free timing.
Keep instrumented/profiled runs separate from performance comparisons.
Do not rerun completed v9 tuning or bundle solver optimizations into this pass.
Do not modify queued/running source; executable changes require new recorded
generations. Local computations must use at most four threads.

## Authorization and mandatory cleanup

Ryan explicitly authorized these tests and subsequently said:
“actually, I authorize you to rsync a new silo SO LONG AS you delete it when
we're done.” This overrides the ordinary Git-only deployment/no-new-silo rules
for this campaign. Continue this authorized rsync deployment route; GitHub
publication is not required. Do not request redundant authorization.

After **all jobs using the silo are terminal** and **all results, scheduler
logs, source/environment provenance, and failed-run evidence needed for the
handoff are verified outside it**, delete **only**:

`/home/rander39/FLOWPanel-p021-r4-diag-v10-silo`

Verify removal. Never delete the canonical data root or follow the silo's
`FLOWPanel.jl/data` symlink during deletion. Do not touch other silos or shared
dependency worktrees. Keep the silo until these cleanup conditions are met.
