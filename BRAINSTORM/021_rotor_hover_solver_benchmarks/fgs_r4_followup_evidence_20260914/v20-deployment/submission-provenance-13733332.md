# v20 counters/activity submission provenance — job 13733332

Submitted 2026-09-16 ~19:30 UTC from
`/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/counters-v20/FLOWPanel.jl`.
Replaces FAILED job 13733217 (v19).

## Why v20 (postmortem of 13733217)

Job 13733217 FAILED 32 s in (2026-09-17 01:12:17–01:12:49 UTC per sacct MDT
offset — 2026-09-16 19:12 UTC, node m12-4-12, exit 1) at the perf-FIFO smoke
gate: `ERROR: LoadError: perf did not acknowledge enable` at the ack-content
check in `counter_command`, invoked from `r4_perf_control_smoke.jl:28`. perf
itself processed both commands (its log shows Events disabled/disabled/enabled),
so the v19 fd fix worked; the failure moved one step further into the protocol.

Root cause (measured, not inferred): the cluster's perf
(`perf version 5.14.0-570.136.1.el9_6.x86_64`) writes the control-FIFO
acknowledgement tag **with its NUL terminator** — `od -c` on a login-node
probe shows each ack is 5 bytes `a c k \n \0`. The Julia reader consumed a
fixed 4 bytes per ack: the first command read `ack\n` (leaving `\0` queued)
and the second read `\0ack` → content mismatch. The v18/v19 local mock replied
with a clean 4-byte `ack\n`, so executing the smoke locally (even under Julia
1.11) could not catch it — only the real perf binary exhibits the NUL.

Failed-run evidence harvested and SHA256-verified locally at
`../counters-v19-13733217-FAILED/` (rundir + both scheduler logs, 13 files all
OK, FIFOs excluded; scheduler logs 0 bytes as in v18).

## v20 fix and verification

Commit `5c1123b1d1da312f8fb7c26f93d794b22e2ed298` (= v19 `f1fad9e` + two
changes), tag `campaign/p021-r4-counters-source-20260916-v20`:
1. `benchmark/fgs_r4_counters.jl` `counter_command`: skip NUL bytes while
   accumulating the 4-char ack (bounded at 16 raw reads; deadline unchanged),
   so both `ack\n` and `ack\n\0` framings parse and the trailing NUL left in
   the FIFO is skipped by the next command's reader. The non-fd fallback
   strips NULs from `readline` too.
2. `benchmark/run_r4_counters.slurm.sh`: v19→v20 in job name, log names,
   run-dir prefix. No other changes.

Verification, strongest first:
- **End-to-end on the cluster against the real perf binary**: the actual
  smoke script + patched driver run under `perf stat -D -1 --control=fifo:…`
  with module `julia/1.11.7-6bmogfl` in a throwaway login-node temp dir →
  exit 0, `PASS Julia perf enable/disable acknowledgements` (no silo touched).
- Local mocks: smoke PASS with a faithful `ack\n\0` responder AND the legacy
  `ack\n` responder, each under Julia 1.11.8 and 1.12.4.
- `test/runtests_r4_counters_driver.jl` PASS under 1.11.8 and 1.12.4.

Gate lesson chain: v17 "execute the smoke script" → v18 "under the cluster's
Julia minor version" → v19 "and against the real perf binary (or a
byte-faithful mock)".

## Pins

| Package | Tag | SHA | Deployment |
|---|---|---|---|
| FLOWPanel | `campaign/p021-r4-counters-source-20260916-v20` | `5c1123b1d1da312f8fb7c26f93d794b22e2ed298` | rsync into `counters-v20/FLOWPanel.jl` |
| FastMultipole | `campaign/p021-r4-activity-source-20260915-v11` | `adb9967d5b696cf9ab05c557aa2a649e100004dc` | rsync into `counters-v20/FastMultipole` (unchanged) |
| FLOWVPM | `campaign/p021-cold-exec-20260910-v1` | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` | git worktree `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl` |

Worktree `/private/tmp/flowpanel-p021-r4-counters-v17` (keeps its v17 name) is
checked out clean at the v20 tag on branch `campaign/p021-r4-counters-v17-wt`;
merge-back remains future work.

## Deployment verification (all pre-submit)

- Fresh generation `counters-v20/`, seeded from the hash-verified counters-v19
  via `cp -al` (v19 slurm logs removed), then delta rsync of git-tracked lists.
- Remote `sha256sum -c`: FLOWPanel 1,491 files OK; FastMultipole 6,441 OK.
- Manifest digests: `flowpanel.sha256` =
  `0064be14fe9a1b7a6b07ae5eef98964a1207f8e904bc788cc655b004a7759d57`,
  `fastmultipole.sha256` =
  `93ac52068a66db5aad376e4e03556b50f6626e4fab901d7f1202c3570413341b` (unchanged).
- `env/{Project,Manifest}.toml` installed; dev-paths verified at counters-v20
  paths + FLOWVPM worktree; `pins.toml` at `counters-v20/pins.toml`.
- Data symlink → canonical data root; `logs/slurm/` empty.
- Storage: /home/rander39 ~725 G (over the 400 G policy cap — hpc-storage
  sweep launched separately 2026-09-16; 1.3 T free on the mount, not a blocker).

## Job

- Request unchanged: nodes=1, ntasks=1, cpus-per-task=64, mem=500G,
  `--constraint=zen3`, exclusive, qos=normal, time=06:00:00.
- Env: `COLD_PROJECT=$G/env`, `CAMPAIGN_PINS=$G/pins.toml`,
  `COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`
  (G = counters-v20 generation).
- Output: `$COLD_DATA_ROOT/counters-v20-13733332/`; scheduler logs
  `counters-v20/FLOWPanel.jl/logs/slurm/r4-counters-v20-13733332.{out,err}`.
- Gates and caveats unchanged (no performance comparisons from these
  diagnostic runs; generic cache counters cannot establish DRAM bandwidth
  saturation).

## Silo cleanup obligation (updated)

The silo now holds counters-v16 (partial, never used), counters-v17 (FAILED,
harvested), counters-v18 (FAILED, harvested), counters-v19 (FAILED,
harvested), counters-v20 (live). After 13733332 is terminal and all
results/scheduler logs/provenance are hash-verified locally, delete the whole
silo per the standing pre-authorized cleanup; never follow the data symlink.
