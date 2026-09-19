# v21 colored-sweep A/B — submission provenance (2026-09-17)

Campaign phase directed by Ryan 2026-09-17 (`fgs_r4_context_reset_20260917.md`).
Question: does the conflict-free colored leaf-sweep order parallelize the
~9.5 s serial `nearfield_update` chain (85% of j64 wall) and beat serial
lexicographic on total time to accepted accuracy?

## Pins (recorded pre-submit)

| Package | Tag | SHA | Deployment |
|---|---|---|---|
| FLOWPanel | `campaign/p021-r4-colored-source-20260917-v21` | `5f87a0e7df97f178c98d93019f2d845c942ff97f` | git worktree `/home/rander39/campaigns/p021-r4-colored-20260917-v21/FLOWPanel.jl` |
| FastMultipole | `campaign/p021-r4-colored-fm-20260917-v21` | `c6185cdde2fe9714a6e60f506dc4a7115e28d8cc` | git worktree `/home/rander39/campaigns/p021-r4-colored-20260917-v21/FastMultipole` |
| FLOWVPM | `campaign/p021-cold-exec-20260910-v1` | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` | existing git worktree `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl` |

Pin decisions:
- FLOWPanel v21 = merged `fastmultipole` tip (`dad2ceb` lineage) + one new
  commit adding the v21 harness (driver `benchmark/fgs_r4_colored_ab.jl`,
  launcher `benchmark/run_r4_colored_ab.slurm.sh`, unit test
  `test/runtests_r4_colored_ab_driver.jl`).
- FastMultipole v21 tag sits at the MERGED `flowpanel-20260817` tip
  (`c6185cdd`), not the pre-merge v11 tag (`adb9967d`): the 2216-case
  colored-sweep suite and the 2026-09-17 local validation ran on the merged
  line, so retaining v11 would have pinned code different from what was
  validated. FastMultipole source content is otherwise unchanged since v11.
- Worktrees created with plain `git worktree add -b <tag>-wt` at the tag
  (HEAD == tag commit, clean status — required by `cold_packages`); no
  data-symlink commit, since every v21 path (OUTDIR, BENCH_CASE_ROOT,
  CONFIG_FILE, COLD_DATA_ROOT) is absolute.
- **Tags pushed directly to the cluster clones** (`~/projects/FLOWPanel.jl`,
  `~/projects/FastMultipole`), NOT to origin: the handoff stop condition
  forbids pushing the merged-but-unpushed branches to origin without Ryan.
  Flagged for Ryan: once the branches are approved for pushing, push these
  tags to origin too (`git push origin campaign/p021-r4-colored-source-20260917-v21`,
  FastMultipole analogue) to satisfy the "pinnable by anyone" rule.

Julia env: `/home/rander39/campaigns/p021-r4-colored-20260917-v21/env`
(resolved+instantiated under module `julia/1.11.7-6bmogfl` on login02,
`JULIA_PKG_PRECOMPILE_AUTO=0`; Manifest dev-paths verified at the three
campaign worktrees). Pins file: `.../p021-r4-colored-20260917-v21/pins.toml`.

## Design (as launched)

Launcher `benchmark/run_r4_colored_ab.slurm.sh`: 64 CPU zen3 exclusive 500G
`--qos=normal` 12 h. Steps: controls (new driver unit test, cold_parse,
cold_precompile, benchmark cold j1/j4, FLOWPanel solver+history unit tests,
FMM chain gravitational+solve_test+fgs_coloring_test) → colored tolerance
calibration at j64 (AB_MODE=calibrate; separate calibration because coloring
changes accumulation ordering) → uninstrumented A/B trials at j∈{1,4,16,64}
(AB_MODE=trials; 8 alternating batches of 10 → 40 trials/order/arm;
direct-vs-FMM crosschecked warmups per order per arm) → one j64
activity-instrumented pair (AB_MODE=activity; budget attribution only,
excluded from rankings). perf dropped entirely. Rank ONLY by total time to
accepted accuracy under the §6 gates.

Yardsticks: j64 lexicographic uninstrumented median 10.96 s (v15); Amdahl
bound if coloring were free ≈ 2.2 s (median 16-way over the 9.3 s chain);
79 colors × 81 inner sweeps ≈ 6.4k barriers/solve expected sync cost.

## Pre-submit gate (v17–v19 lesson chain)

Executed locally under Julia **1.11.8** (cluster minor line; no external
binary protocols in v21 — perf dropped):
- `test/runtests_r4_colored_ab_driver.jl` PASS under 1.11.8 AND 1.12.4
  (parse, dynamic init, pair check, alternation, `bash -n` launcher, no-perf
  assertions).
- Driver executed END-TO-END at R4, all three modes, j4/BLAS1, dedicated
  1.11.8 env dev-pathed at the local checkouts:
  - calibrate: colored tolerance staircase → 5.22316207998051e-7
    (lexicographic retained: 3.479128881193055e-7), 26 iterations (vs 27
    lex — the iterate path moved off 27 as §7 predicted), confirmation
    repeat delta 0.0, all gates green, `status = completed`.
  - trials (batches=2 minimum): 40/40 trials certified; local-only signal
    (NOT a campaign measurement): colored median 9.62 s vs lex 11.51 s.
    Cross-order solution rel-L2 1.53e-6 (both orders separately certified
    ≤1e-6 BC rel-L2; informational).
  - activity: both orders solved/certified, repeat delta 0,
    `status = completed`. Per-thread CSVs empty locally by design (macOS has
    no /proc; guard only active off-Linux — cluster path unchanged).
- Unchanged launcher steps (cold_parse/cold_precompile/benchmark-cold/unit
  suites/FMM chain) are cluster-proven from v15/v20 and passed in the
  2026-09-17 local merged-line validation.

## Storage preflight (2026-09-17)

hpc-storage cycle before submission: /home 644→654 G of the 400 G cap;
**0 MB archivable automatically** — every sizeable run is PROTECTED,
RECENT-HOT (live writers), or RECENT (15 runs, ~369.5 GiB VTK, quiet
10–22 h, awaiting Ryan's `--include-recent --only` approval).
STALE_COUNT=0, VERIFY_FAIL_COUNT=0. `data/p021-cold-20260910/` is ~1.3 GB
(NO-VTK) — not the disk driver. v21 writes only CSV/TOML output (~0.3 GB
expected, no VTK), so submission proceeds with the cap breach flagged to
Ryan rather than blocked on it.

## Submission

Submitted from the FLOWPanel campaign worktree top level with
`COLD_PROJECT=$camp/env CAMPAIGN_PINS=$camp/pins.toml
COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`.
`sbatch --test-only` ETA at submission time: start 2026-09-17T17:16 UTC on
m12 (zen3). Job ID recorded below after sbatch.

- Job: **13738561** (submitted 2026-09-17 from login02; ETA at submission ~2026-09-17T17:16 UTC on m12)

## Postmortem: job 13738561 FAILED (env-only) → resubmitted as 13738665

Job 13738561 started ~00:15 UTC 09-17 (backfill, well before the 17:16 ETA)
and died in the `test/runtests_unit_solver.jl` control:
`ArgumentError: Package Meshes not found in current path` (line 4
`import Meshes`). Root cause: the FLOWPanel unit-test controls import
`Meshes` and `StaticArrays` directly, which are dependencies OF FLOWPanel
but not direct dependencies of the campaign env, so `--project=env` cannot
resolve them. The 2026-09-17 local merged-line validation had run these
tests under `--project=.` (FLOWPanel's own project, where they are direct),
which masked the gap; the v21 local gate executed the NEW scripts under the
campaign-style env but relied on that earlier validation for the unchanged
controls. **Gate lesson (extends v17–v19 chain): "execute every script the
launcher invokes" means under the campaign env itself — unchanged scripts
can still fail from env composition.**

Handling per protocol:
- Evidence harvested first: `../colored-v21-13738561-FAILED/` (17 files,
  SHA256-verified against remote, `remote.sha256`).
- Reproduced locally: `import Meshes` under the local 1.11.8 gate env fails
  identically.
- Fix is env-only — **no source change, no new tags, worktrees untouched**
  (still clean at the pinned SHAs): `Pkg.add Meshes StaticArrays` in both
  the local gate env and the cluster campaign env (cluster resolves Meshes
  0.42.2, matching FLOWPanel's compat pin; imports verified under 1.11.7).
- Re-gated locally under the fixed env (Julia 1.11.8):
  `runtests_unit_solver.jl` PASS, `runtests_unit_fgs_history.jl` PASS,
  FMM chain gravitational+solve_test+fgs_coloring_test PASS
  (157 + 2216 colored cases).
- Resubmitted from the same worktree with the same pins: **job 13738665**
  (2026-09-17). Watch re-armed.
