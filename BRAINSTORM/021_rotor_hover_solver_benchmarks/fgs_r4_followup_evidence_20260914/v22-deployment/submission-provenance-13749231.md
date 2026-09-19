# v22 chunked-sweep A/B — submission provenance (2026-09-17)

Campaign phase per Ryan-approved plan `fgs_chunked_hybrid_plan_20260918.md`
(approved 2026-09-17; entry `fgs_r4_context_reset_20260918.md`). Question:
does the chunked hybrid sweep (`sweep_order=:chunked` — Gauss-Seidel within
cost-balanced contiguous leaf chunks, Jacobi across chunks via deferred
cross-chunk scatter, ONE barrier per inner sweep instead of colored's ~6.4k)
parallelize the ~9.3 s serial `nearfield_update` chain (85% of the 11.2 s
j64 wall) and beat BOTH yardsticks — lex@j64 10.96–11.20 s and colored@j16
10.116 s — on total time to accepted accuracy?

## Pins (recorded pre-submit)

| Package | Tag | SHA | Deployment |
|---|---|---|---|
| FLOWPanel | `campaign/p021-r4-chunked-source-20260918-v22` | `67d570f576e8698e82a07a889aeaac72c0c9d8d4` | git worktree `/home/rander39/campaigns/p021-r4-chunked-20260918-v22/FLOWPanel.jl` |
| FastMultipole | `campaign/p021-r4-chunked-fm-20260918-v22` | `c18e4b46116825a201a0cf459b7596242191dcae` | git worktree `/home/rander39/campaigns/p021-r4-chunked-20260918-v22/FastMultipole` |
| FLOWVPM | `campaign/p021-cold-exec-20260910-v1` | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` | existing git worktree `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl` |

Pin decisions:
- FLOWPanel v22 = v21 tip (`5f87a0e`) + one new commit `67d570f` (FGSSolver/
  FGSPreconditioner `chunks` plumbing + metadata, `fgs_cold_common.jl`
  chunked validation + screen axis, v22 harness trio
  `benchmark/fgs_r4_chunked_ab.jl` / `benchmark/run_r4_chunked_ab.slurm.sh` /
  `test/runtests_r4_chunked_ab_driver.jl`, solver-suite chunked plumbing case).
- FastMultipole v22 = v21 tip (`c6185cdd`) + one new commit `c18e4b46`
  (`sweep_order=:chunked` + `chunks::Int=64` kwarg, `build_chunk_map`,
  intra/cross scatter partition, chunked `gs_sweep!` branch,
  `test/fgs_chunked_test.jl` + `test/fgs_chunked_threadcheck.jl`). Coloring
  is still present in the pinned tree (§5 revert is conditional and LAST).
  The unsplit scatter path used by :lexicographic and :colored is unchanged.
- Worktrees created with plain `git worktree add -b <tag>-wt` at the tags
  (HEAD == tag commit, clean status — required by `cold_packages`); no
  data-symlink commit, since every v22 path (OUTDIR, BENCH_CASE_ROOT,
  CONFIG_FILE, COLD_DATA_ROOT) is absolute.
- **Tags pushed directly to the cluster clones** (`~/projects/FLOWPanel.jl`,
  `~/projects/FastMultipole`), NOT to origin: origin pushes of the merged
  branches remain Ryan-pending (v21 precedent). Flagged again: on approval,
  push `campaign/p021-r4-chunked-source-20260918-v22` and the FastMultipole
  analogue to origin to satisfy the "pinnable by anyone" rule.

Julia env: `/home/rander39/campaigns/p021-r4-chunked-20260918-v22/env`
(resolved+instantiated under module `julia/1.11.7-6bmogfl` on login03,
`JULIA_PKG_PRECOMPILE_AUTO=0`; Manifest dev-paths verified at the three
campaign worktrees; **`Meshes` 0.42.2 + `StaticArrays` added as direct deps**
per the 13738561 lesson). Pins file:
`.../p021-r4-chunked-20260918-v22/pins.toml`.

## Design (as launched)

Launcher `benchmark/run_r4_chunked_ab.slurm.sh`: 64 CPU zen3 exclusive 500G
`--qos=normal` 12 h. Steps: controls (new driver unit test, cold_parse,
cold_precompile, benchmark cold j1/j4, FLOWPanel solver+history unit tests,
FMM chain gravitational+solve_test+fgs_coloring_test+**fgs_chunked_test**) →
chunked tolerance calibration at j64 (AB_MODE=calibrate; separate
calibration because the hybrid GS/Jacobi iteration changes accumulation
ordering and iteration character; fixed `chunks=64` makes the chunk map
j-invariant so ONE calibration carries across arms, §4.Q2) →
uninstrumented A/B trials at j∈{1,4,16,64} (AB_MODE=trials; 8 alternating
batches of 10 → 40 trials/order/arm; direct-vs-FMM crosschecked warmups per
order per arm) → one j64 activity-instrumented pair (AB_MODE=activity;
budget attribution only, excluded from rankings — did the chain span shrink
and at what average active width, the number colored failed on). perf
dropped entirely. Rank ONLY by total time to accepted accuracy under the §6
gates. Two-way A/B (chunked vs lexicographic); colored medians cited from
v21 (banked, gate-certified).

Yardsticks: lex@j64 uninstrumented median 10.96–11.20 s; colored@j16
10.116 s (v21). Expectation from the DRAM-bandwidth cap: chain ~1.2–2 s,
j64 wall ~3.5–4.5 s; ~81 barriers/solve (one per inner sweep) at full j-way
width. Divergence risk (majority-Jacobi at 64 chunks × ~17 leaves) is
handled by the separate calibration: a staircase with no certified crossing
is a finding, not a knob-turning license (§11).

## Pre-submit gate (v17–v21 lesson chain)

Executed locally under Julia **1.11.8** in a FRESH campaign-style env
(`Pkg.develop` FLOWPanel.jl + FastMultipole + FLOWVPM.jl local checkouts,
**`Pkg.add Meshes StaticArrays`** — the 13738561 lesson applied from the
start), every launcher-invoked script including unchanged ones:
- `test/runtests_r4_chunked_ab_driver.jl` PASS under 1.11.8 AND 1.12.4
  (parse, dynamic init, pair check incl. chunks-mismatch cases, alternation,
  `bash -n` launcher, no-perf assertions).
- Control chain under the campaign-style env: cold_parse, cold_precompile,
  benchmark cold j1 + j4, `runtests_unit_solver.jl` (incl. new :chunked
  plumbing case), `runtests_unit_fgs_history.jl`, FMM chain
  gravitational+solve_test+fgs_coloring_test (2216) + fgs_chunked_test
  (1723, incl. cross-process -t 1/-t 4 bitwise invariance) — ALL PASS.
- Driver executed END-TO-END at R4, all three modes, j4/BLAS1:
  - calibrate: chunked tolerance staircase → **5.348662428097527e-7**
    (lexicographic retained: 3.479128881193055e-7), **44 iterations** (vs
    27 lex — iteration count moved off 27 as §1.4 anticipated; the ranking
    metric absorbs it), confirmation repeat delta **0.0**, all gates green,
    `status = completed`.
  - trials (batches=2 minimum): 40/40 trials certified; local-only signal
    (NOT a campaign measurement): chunked median 19.26 s vs lex 15.82 s at
    j4 (4 threads / 64 chunks — the campaign target is the j64 arm).
    Cross-order solution rel-L2 3.12e-6 (informational; both orders
    separately certified).
  - activity: both orders solved/certified, `status = completed`. Per-thread
    CSVs empty locally by design (macOS has no /proc; guard only active
    off-Linux — cluster path unchanged).

## Storage preflight (2026-09-17)

hpc-storage cycle before submission: /home **725.9 G → 359.1 G** of the
400 G cap (15 finished runs archived, 272.3 GiB tarballed, 375.6 GiB VTK
freed; newest 5 restartable steps retained per run; STALE_COUNT=0,
VERIFY_FAIL_COUNT=0). **Cap breach RESOLVED — submission proceeds under the
cap** (no flag needed, unlike v21). 3 runs remain RECENT (quiet 17–23 h,
208.0 GiB VTK) awaiting Ryan's `--include-recent --only` approval:
`p018_csarc_l3p0_3r_g25_s2`, `p018_csarc_n2_nt72_l3p0_3r_srlx_g25_s2`,
`p018_csarc_n2_nt72_l3p0_3r_srlx_nv_g25`. v22 writes only CSV/TOML
(~0.3 GB expected, no VTK).

## Submission

Submitted from the FLOWPanel campaign worktree top level with
`COLD_PROJECT=$camp/env CAMPAIGN_PINS=$camp/pins.toml
COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`.
slurm-availability probe (64c/500G/12h): m12 `--qos=normal` ETA immediate.
`sbatch --test-only` at submission: start 2026-09-17T23:21:13 on m12-3-31
(zen3, partition m12).

- Job: **13749231** (submitted 2026-09-17 from login03; ETA at submission
  2026-09-17T23:21 on m12-3-31). Run dir:
  `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/chunked-v22-13749231/`.
