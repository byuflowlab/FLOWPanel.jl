# v18 counters/activity submission provenance — job 13712587

Submitted 2026-09-16 ~04:20 UTC (2026-09-15 ~22:20 MDT) from
`/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/counters-v18/FLOWPanel.jl`.
Replaces FAILED job 13711596 (v17).

## Why v18 (postmortem of 13711596)

Job 13711596 FAILED 40 s in (2026-09-16 03:43:30–03:44:10 UTC, node m12-4-11,
exit 1) at the perf-FIFO smoke gate: `test/r4_perf_control_smoke.jl` line 2 had
a triple-quoted string directly above `using Test`, which Julia parses as a
docstring for the `using` statement and rejects (`ERROR: LoadError: cannot
document the following expression: using Test`). No arms ran; perf-smoke.csv
shows all events `<not counted>`.

Why local testing missed it: the local v17 gate ran only
`test/runtests_r4_counters_driver.jl` (driver parse/protocol/launcher tests);
the smoke script's only caller was the Slurm launcher, so it was never executed
anywhere before submission. The v17 claim "Julia perf smoke passed locally" was
inaccurate.

Failed-run evidence harvested and SHA256-verified locally at
`../counters-v17-13711596-FAILED/` (rundir + both scheduler logs +
`remote-sha256-full.txt`).

## v18 fix and verification

Commit `6033dbd` (= v17 `b0b6eec` + two changes):
1. `test/r4_perf_control_smoke.jl`: docstring demoted to a `#` comment.
2. `benchmark/run_r4_counters.slurm.sh`: v17→v18 in job name, log names,
   run-dir prefix. No other changes.

Local verification this time EXECUTED the smoke script: run against a mock
perf FIFO responder (3-command disable/enable/disable sequence with `ack\n`
replies), exit 0, `PASS Julia perf enable/disable acknowledgements`; driver
test suite re-passed. Sibling scripts scanned for the same docstring-before-
`using` pattern: none.

## Pins

| Package | Tag | SHA | Deployment |
|---|---|---|---|
| FLOWPanel | `campaign/p021-r4-counters-source-20260915-v18` | `6033dbdc06df3519e6c623245bcd747dad064555` | rsync into `counters-v18/FLOWPanel.jl` |
| FastMultipole | `campaign/p021-r4-activity-source-20260915-v11` | `adb9967d5b696cf9ab05c557aa2a649e100004dc` | rsync into `counters-v18/FastMultipole` (unchanged from v17) |
| FLOWVPM | `campaign/p021-cold-exec-20260910-v1` | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` | git worktree `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl` (verified clean at pin pre-submit) |

Tags/commits live in the main local repos (worktree
`/private/tmp/flowpanel-p021-r4-counters-v17` now checked out at the v18 tag;
the directory keeps its v17 name).

## Deployment verification (all pre-submit)

- Fresh generation `counters-v18/`, seeded from the hash-verified counters-v17
  via `cp -al` (v17 slurm logs removed from the copy), then delta rsync of the
  git-tracked file lists.
- Remote `sha256sum -c`: FLOWPanel 1,491 files OK; FastMultipole 6,441 files OK.
- Manifest digests match local prepare output:
  `flowpanel.sha256` = `037d000e416028ca8e5c0206d104c3c56397a8fba3a39553043b542e77c7da52`,
  `fastmultipole.sha256` = `93ac52068a66db5aad376e4e03556b50f6626e4fab901d7f1202c3570413341b`
  (FastMultipole digest identical to v17 — source unchanged, pin retained).
- `env/{Project,Manifest}.toml` installed; Manifest dev-paths verified remotely
  to point at the two counters-v18 package paths and the FLOWVPM worktree.
- `pins.toml` at `counters-v18/pins.toml`.
- Data symlink `counters-v18/FLOWPanel.jl/data -> /home/rander39/projects/FLOWPanel.jl/data`;
  `logs/slurm/` created empty.
- Storage preflight: /home/rander39 at 365 G of 400 G cap; filesystem 1.7 T free.

## Job

- Request unchanged from v17: nodes=1, ntasks=1, cpus-per-task=64, mem=500G,
  `--constraint=zen3`, exclusive, qos=normal, time=06:00:00.
- Env: `COLD_PROJECT=$G/env`, `CAMPAIGN_PINS=$G/pins.toml`,
  `COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`
  (G = counters-v18 generation).
- Output: `$COLD_DATA_ROOT/counters-v18-13712587/`; scheduler logs
  `counters-v18/FLOWPanel.jl/logs/slurm/r4-counters-v18-13712587.{out,err}`.
- Gates and caveats unchanged from v17 provenance (no performance comparisons
  from these diagnostic runs; generic cache counters cannot establish DRAM
  bandwidth saturation).

## Silo cleanup obligation (updated)

The silo now holds generations counters-v16 (partial, never used),
counters-v17 (FAILED job, evidence harvested), counters-v18 (live). After
13712587 is terminal and all results/scheduler logs/provenance are
hash-verified locally, delete the whole silo per the standing pre-authorized
cleanup; never follow the data symlink.
