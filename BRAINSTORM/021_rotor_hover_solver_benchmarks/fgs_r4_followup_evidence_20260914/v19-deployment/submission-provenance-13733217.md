# v19 counters/activity submission provenance — job 13733217

Submitted 2026-09-16 ~19:15 UTC from
`/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/counters-v19/FLOWPanel.jl`.
Replaces FAILED job 13712587 (v18).

## Why v19 (postmortem of 13712587)

Job 13712587 FAILED 41 s in (2026-09-16 04:13:00–04:13:41 UTC, exit 1) at the
perf-FIFO smoke gate: `ERROR: LoadError: bitcast: argument size does not match
size of target type` at `counter_command` line 67 of
`benchmark/fgs_r4_counters.jl` (`reinterpret(Cint, fd(acknowledgement))`),
invoked from `test/r4_perf_control_smoke.jl:25`. No controls or arms ran;
perf-smoke.csv shows all events `<not counted>`; scheduler logs empty.

Root cause: `fd(::IOStream)` returns a 64-bit `Int` on Julia ≤1.11 (the
cluster module is julia/1.11.7) but a 32-bit `RawFD` on ≥1.12. A
`reinterpret(Cint, ·)` bitcast requires a 32-bit argument, so it throws on the
cluster and passes locally on 1.12. Why the v18 local gate missed it: the v18
lesson was "execute the smoke script locally", which was done — but on Julia
1.12.5, where the bitcast is legal. The failure was reproduced exactly under
local Julia 1.11.8 with a mock perf-FIFO responder before fixing.

Failed-run evidence harvested and SHA256-verified locally at
`../counters-v18-13712587-FAILED/` (rundir + both scheduler logs +
`remote-sha256-full.txt`; 13 files all OK, FIFOs excluded).

## v19 fix and verification

Commit `f1fad9e61b90cb8216f4d7c3aae8ddbbc102dc10` (= v18 `6033dbd` + two
changes), tag `campaign/p021-r4-counters-source-20260916-v19`:
1. `benchmark/fgs_r4_counters.jl` `counter_command`: branch on the type
   returned by `fd` — `RawFD` → `reinterpret(Cint, ·)`, otherwise
   `Cint(raw)`. Version-agnostic across Julia 1.11/1.12.
2. `benchmark/run_r4_counters.slurm.sh`: v18→v19 in job name, log names,
   run-dir prefix. No other changes.

Local verification EXECUTED the smoke script against a mock perf-FIFO
responder (3-command disable/enable/disable, `ack\n` replies) under BOTH
`julia-1.11.8` (cluster line — reproduced the failure pre-fix, PASS post-fix)
and `julia-1.12.4` (PASS), and re-ran
`test/runtests_r4_counters_driver.jl` under both (PASS). New gate lesson:
local pre-submit checks must run under the cluster's Julia minor version.

## Pins

| Package | Tag | SHA | Deployment |
|---|---|---|---|
| FLOWPanel | `campaign/p021-r4-counters-source-20260916-v19` | `f1fad9e61b90cb8216f4d7c3aae8ddbbc102dc10` | rsync into `counters-v19/FLOWPanel.jl` |
| FastMultipole | `campaign/p021-r4-activity-source-20260915-v11` | `adb9967d5b696cf9ab05c557aa2a649e100004dc` | rsync into `counters-v19/FastMultipole` (unchanged from v17/v18) |
| FLOWVPM | `campaign/p021-cold-exec-20260910-v1` | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` | git worktree `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl` |

Tags/commits live in the main local repos; worktree
`/private/tmp/flowpanel-p021-r4-counters-v17` (directory keeps its v17 name)
is checked out clean at the v19 tag on branch
`campaign/p021-r4-counters-v17-wt`. Merge-back to the campaign line remains
future work.

## Deployment verification (all pre-submit)

- Fresh generation `counters-v19/`, seeded from the hash-verified counters-v18
  via `cp -al` (v18 slurm logs removed from the copy), then delta rsync of the
  git-tracked file lists from both local worktrees.
- Remote `sha256sum -c`: FLOWPanel 1,491 files OK; FastMultipole 6,441 files
  OK.
- Manifest digests from local prepare output:
  `flowpanel.sha256` = `e554f00f367fc61e16207f123017b127324956de5dda45d6ec52400e4370c335`,
  `fastmultipole.sha256` = `93ac52068a66db5aad376e4e03556b50f6626e4fab901d7f1202c3570413341b`
  (FastMultipole digest identical to v17/v18 — source unchanged, pin retained).
- `env/{Project,Manifest}.toml` installed; Manifest dev-paths verified remotely
  at the two counters-v19 package paths and the FLOWVPM worktree.
- `pins.toml` at `counters-v19/pins.toml`.
- Data symlink `counters-v19/FLOWPanel.jl/data -> /home/rander39/projects/FLOWPanel.jl/data`;
  `logs/slurm/` empty.
- Storage preflight: /home/rander39 at **725 G** (2.0 T mount, 1.3 T free) —
  over the 400 G policy cap (flagged for an hpc-storage sweep; not a blocker
  for this small-output job).

## Job

- Request unchanged from v17/v18: nodes=1, ntasks=1, cpus-per-task=64,
  mem=500G, `--constraint=zen3`, exclusive, qos=normal, time=06:00:00.
- Env: `COLD_PROJECT=$G/env`, `CAMPAIGN_PINS=$G/pins.toml`,
  `COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`
  (G = counters-v19 generation).
- Output: `$COLD_DATA_ROOT/counters-v19-13733217/`; scheduler logs
  `counters-v19/FLOWPanel.jl/logs/slurm/r4-counters-v19-13733217.{out,err}`.
- Gates and caveats unchanged (no performance comparisons from these
  diagnostic runs; generic cache counters cannot establish DRAM bandwidth
  saturation).

## Silo cleanup obligation (updated)

The silo now holds generations counters-v16 (partial, never used), counters-v17
(FAILED, evidence harvested), counters-v18 (FAILED, evidence harvested),
counters-v19 (live). After 13733217 is terminal and all results/scheduler
logs/provenance are hash-verified locally, delete the whole silo per the
standing pre-authorized cleanup; never follow the data symlink.
