# v17 counters/activity submission provenance — job 13711596

Submitted 2026-09-16 ~02:55 UTC (2026-09-15 evening MDT) from
`/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/counters-v17/FLOWPanel.jl`.

## Pins

| Package | Tag | SHA | Deployment |
|---|---|---|---|
| FLOWPanel | `campaign/p021-r4-counters-source-20260915-v17` | `b0b6eec150183c174f67fee80574c3ca95e77751` | rsync into `counters-v17/FLOWPanel.jl` |
| FastMultipole | `campaign/p021-r4-activity-source-20260915-v11` | `adb9967d5b696cf9ab05c557aa2a649e100004dc` | rsync into `counters-v17/FastMultipole` |
| FLOWVPM | `campaign/p021-cold-exec-20260910-v1` | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` | git worktree (verified clean at pin pre-submit) |

v17 commit = v16 (`5460376`) plus the bounded perf-FIFO acknowledgement fix:
poll(2)+read(2) with one 15 s monotonic deadline per command, Julia perf smoke
(`test/r4_perf_control_smoke.jl`) replacing the Python one in the launcher, and
real-FIFO ack/missing/nack/partial driver tests (all passed locally,
`counter-driver-local-v17.log`). A stray CLAUDE.md edit in the worktree was
discarded, not committed. Tags/commits live in the main local repos (worktrees
`/private/tmp/flowpanel-p021-r4-counters-v17`, `/private/tmp/fastmultipole-p021-r4-activity-v11`).

## Deployment verification (all pre-submit)

- Fresh generation dir `counters-v17/`; the partial `counters-v16/` was not
  reused as a generation (its completed FLOWPanel rsync only seeded v17 via
  `cp -al` before delta rsync).
- Remote `sha256sum -c`: FLOWPanel 1,491 files OK; FastMultipole 6,441 files OK.
- Manifest digests match local prepare output:
  `flowpanel.sha256` = `6ef53cd4cae539c859a61cbbbc7b5ef1ff0d3252f32bf1638bcf577f0ff41282`,
  `fastmultipole.sha256` = `93ac52068a66db5aad376e4e03556b50f6626e4fab901d7f1202c3570413341b`
  (FastMultipole digest identical to v16 — source unchanged, pin retained).
- `env/{Project,Manifest}.toml` installed; Manifest dev-paths point at the two
  counters-v17 package paths (verified remotely) and the FLOWVPM campaign worktree.
- `pins.toml` installed at `counters-v17/pins.toml` (copied into the run dir by
  the launcher as `campaign_pins.toml`).
- Data symlink `counters-v17/FLOWPanel.jl/data -> /home/rander39/projects/FLOWPanel.jl/data`;
  `logs/slurm/` created.
- Storage preflight: /home/rander39 at 353 G of 400 G cap before the ~635 MB
  deployment; filesystem 1.7 T free.
- FLOWVPM remote worktree clean at pinned SHA (verified by read-only monitor).

## Job

- Request: nodes=1, ntasks=1, cpus-per-task=64, mem=500G, `--constraint=zen3`
  (node feature, not a partition), exclusive, qos=normal, time=06:00:00.
- Env: `COLD_PROJECT=$G/env`, `CAMPAIGN_PINS=$G/pins.toml`,
  `COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`.
- Output: `$COLD_DATA_ROOT/counters-v17-13711596/`; scheduler logs
  `counters-v17/FLOWPanel.jl/logs/slurm/r4-counters-v17-13711596.{out,err}`.
- Gates inside the job (in order): perf FIFO smoke under installed perf →
  counter-driver controls → FastMultipole solver/coloring tests → FLOWPanel
  solver/history unit controls → j4 and j64 perf-gated prepared solves
  (events disabled except around one warmed solve per arm). No performance
  comparisons are to be drawn from these diagnostic runs.

## Caveats carried forward

Generic cache counters cannot establish DRAM bandwidth saturation; the coarse
`nearfield_update` activity span combines leaf/product/scatter and updates;
/proc tick resolution and sequential snapshots limit short-stage conclusions.

## Silo cleanup obligation (unchanged)

After all jobs using any generation in
`/home/rander39/FLOWPanel-p021-r4-diag-v10-silo` are terminal and all results,
scheduler logs, and provenance are hash-verified locally, delete the whole silo
(including counters-v16/ and counters-v17/), verify removal, never follow the
data symlink. Job 13711596 must be terminal and harvested first.
