# R4 authorized campaign: context-reset handoff

Prepared at Ryan's request, 2026-09-16 02:34 UTC (2026-09-15 MDT).
This supersedes `fgs_r4_context_reset_20260915b.md` for live state. Read that
previous handoff only for unchanged pins, deployment details and policy index.
Ryan requested a reset while work was underway; **no new job was submitted**.

## Mission and context discipline

Continue the authorized campaign through separately pinned counter diagnostics,
measured evidence handoff, and mandatory bounded silo cleanup. Do not resubmit
completed ladder job **13694724**. Keep parent context slim: Luna for mechanical
monitoring/harvesting, Terra for Julia verification, bounded Astra for numerical
review/measurement reasoning/experiment selection. Require compact reports and
artifact indexes. Never load bulk logs/history or full CPU/provenance dumps.
Read `/Users/ryan/.claude/CLAUDE.md`, repository CLAUDE/AGENTS and routed policies;
read current remote BYU ORC instructions and applicable references before work.
No notebook entry without separate approval. Preserve unrelated live changes.

## Completed this turn: v15 harvest and audit

Job 13694724 was already COMPLETED, exit 0:0, elapsed 08:52:49. It was NOT rerun.
All six arms, root COMPLETED, and logs are harvested under:
`fgs_r4_followup_evidence_20260914/diag-v15-13694724/` (relative to this directory).

- `remote-sha256-full.txt`: **131 entries**, exactly 23 root + 96 results files
  (16 per arm) + 12 arm-level files (`process.log`, `numactl_show.txt`).
- `sha256-verification-full.txt`: 131/131 OK. Parent independently recomputed
  all 131 local digests against this remote manifest: zero mismatches.
- Scheduler stdout/stderr were separately checked remote-to-local: both empty,
  SHA256 e3b0c442...2b855. Detailed process logs are in the arms.
- `silo-provenance/`: five additional remote-to-local verified files: pins,
  env Project/Manifest, FLOWPanel v15 and FastMultipole v10 source manifests.
- `harvest_summary-full.md`: corrected compact inventory and provenance.
- `analysis/audit-full.txt`: prepared data audit exit 0; all 240 timed trials,
  all six statuses, profile/equivalence and additive stage checks PASS.
- Tags/SHAs, Julia 1.11.7, BLAS=1, requested arm threads agree across provenance.
  This is not a fresh remote rehash of every source file; distinguish recorded
  content verification from the fresh evidence-file transfer hashes.

| Julia threads | Uninstrumented median s | Instrumented median s |
|---:|---:|---:|
| 1 | 37.5594 | 37.4689 |
| 4 | 16.9358 | 17.0866 |
| 8 | 13.4642 | 13.3843 |
| 16 | 12.1994 | 11.9373 |
| 32 | 11.3673 | 11.4110 |
| 64 | 10.9626 | 10.9483 |

Independent Astra review artifacts are in evidence `analysis/` (NOT the run's
nested analysis). Read its compact review first, then only needed CSV tables.
At j64 nonself products cost about 7.887 s (~72%), FMM about 1.114 s (~10%).
There are 28 outer residual checks, 27 reported update iterations, 81 sweeps.
The j32→64 gain is small and lacks independent process replication.
**Warmup caveat:** every arm's batch 2 trial 1 has outer-minus-internal time
~0.752–0.781 s, consistent with first instrumented/no-callback specialization
compilation. Raw rows are retained; report sensitivity rather than silently
excluding them or claiming all instrumented trials were fully warmed.
The uninstrumented ladder is unaffected. Raw solutions/histories are not exported;
distinguish source-backed gates and exact-history booleans from independent replay.

## Measurement review and v17 correction in progress

`evidence/analysis/measurement_review_20260915b.md` (expand evidence to the root
above) records scope and remaining deliverables. v16 is suitable for generic
whole-solve cache counters and coarse stage CPU activity, NOT measured DRAM
bandwidth or a saturation claim. Perf includes acknowledgement boundary overhead;
its shell execs Julia, and inherited thread scope must be interpreted carefully.
Short /proc tick spans and sequential snapshots limit activity precision.

Review found v16 Julia FIFO acknowledgement could block indefinitely. Python's
smoke was bounded but Julia tests used only IOBuffer. Terra is correcting this in:
`/private/tmp/flowpanel-p021-r4-counters-v17`
branch `campaign/p021-r4-counters-v17-wt`, based on v16 `546037672ac7e41eba5c64cedf6a2e4d6b7394e6`.
**This worktree is currently uncommitted/unpinned. Do not deploy as v16 or claim
v17 is clean/verified until the final verification report and git state are checked.**

Changes: POSIX poll + raw one-byte read loop with monotonic total deadline;
success, missing, malformed and partial ACK tests; actual Julia-under-perf FIFO
smoke before fixture construction; v17 launcher naming. Parent caught partial
ACK blocking in an earlier draft and required the current deadline loop. Julia
unit FIFO writers must use Threads.@spawn (blocking poll prevents @async writers
on that same thread). Terra's final report/log goes in `v17-deployment/`.
Current BYU policy is `/private/tmp/BYU_ORC_AGENTS-current.md`; Terra will add it
and a CLAUDE reference to this new generation before pinning.

FastMultipole remains the clean v11 worktree/tag from the prior handoff; no new
FMM executable changes were requested. FLOWVPM pin also remains unchanged.
Prepared v16 local source/pins and its remote **partial** deployment are unchanged.
No counters-v17 remote directory has been created. No counter job is queued.

## Local deployment preparation / remote access

New `fgs_r4_followup_evidence_20260914/v17-deployment/prepare_deployment.py` builds
NUL-separated tracked regular-file inventories, SHA256 manifests, pins and env
files from CLEAN ANNOTATED worktrees. It has not been run: first commit/tag v17,
review its worktree paths and tag date, then run. It excludes canonical data and
symlinks. It points both packages into a separate remote `counters-v17/` generation.
Use rsync --from0 --files-from, no --delete; verify source hashes remotely before
submission. Reuse the prior v16 environment with the new paths as generated.

SSH existing ControlMaster works; sandbox requires approved `require_escalated`
SSH. Use BatchMode and bounded ConnectTimeout. Stop if authentication expires;
never retry into MFA. Policies were read from current `/apps` copies this turn.

`v17-deployment/availability.csv`: read-only availability at 2026-09-16 02:28 UTC,
64 cores/500 GiB: m11-2 reported two fitting idle nodes; this is NOT a placement
or ETA guarantee. Recheck as needed and use final launcher `sbatch --test-only`.
Keep required Zen3/exclusive/BLAS1 and no unnecessary partition restriction.
A broad home `du -sm` made no progress/output and was interrupted for reset;
`storage-vpm-preflight.txt` is empty, NOT evidence of headroom or VPM validation.
That command's following VPM git checks did NOT run. Tool session 5981 exited 255
after interrupt. No current storage or fresh VPM-clean result was established.

## Next actions, in order

1. Read completed compact numerical review and Terra's v17 verification report.
   Finish any pending narrow checks. Review/pin corrected v17 with a new annotated
   tag; never modify old annotated pins. Local computations ≤4 threads.
2. Generate deployment inventory/env/pins, refresh storage and VPM clean tag/SHA,
   assets and allocation; deploy separate generation under existing authorized
   silo. Verify full remote source hashes, environment/loaded paths and assets.
3. Record provenance first, test allocation, submit corrected counter launcher.
   It must gate expensive fixtures on installed perf + real Julia FIFO smoke and
   FLOWPanel solver/history/FastMultipole controls. Collect j4 and j64 counters
   and activity. No comparisons of counter-run timings as performance trials.
4. Harvest/hash/audit all counter outputs. Check events are supported/counting,
   enabled/running coverage/multiplexing; generic cache misses are NOT DRAM bytes.
5. Reconcile stage budget, placement, dependency/block census and scaling; update
   diagnostics package and provenance. Follow governing 2026-09-14 request for
   conditional separately calibrated colored-order experiment OR concrete
   justified next experiment. Preserve finite/converged, certified FMM BC≤1e-6,
   repeat≤1e-8, evaluated direct/FMM≤1e-7, BLAS1, zero reset, prepared timing.
6. Perform mandatory cleanup only after all conditions below, then final evidence
   handoff. Existing diagnostics package has NOT yet been updated this turn.

## Authorization and mandatory cleanup (unchanged)

Ryan authorized these tests and rsync into a new silo SO LONG AS IT IS DELETED
WHEN DONE. Do not ask again. This exception overrides ordinary no-new-silo/Git-only
deployment rules; clean annotated local pins remain mandatory.
After ALL jobs using ANY generation in the silo are terminal, and ALL results,
scheduler logs, source/environment provenance and failed-run evidence are
hash-verified outside it, delete ONLY:
`/home/rander39/FLOWPanel-p021-r4-diag-v10-silo`
Verify removal. Do not follow data symlinks or delete the canonical data root,
shared dependencies or other silos. This includes unfinished counters-v16 and
any future counters-v17 nested generation. Conditions are NOT yet met; no cleanup
has occurred. Retain earlier failed jobs 13688396/13688430/13688474/13689318/13690579,
passed controls 13690544 and pre-v15 evidence already local.
