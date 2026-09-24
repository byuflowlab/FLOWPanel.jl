# 021 :dagedge HPC benchmark — provenance (2026-09-24)

Mission: `fgs_dagedge_benchmark_reset_prompt_20260924.md`. Two submissions
PRE-APPROVED by Ryan (2026-09-24): (1) an uninstrumented **performance run**
(dagteam+backoff champion vs dagedge, j ladder, θ probe), then (2) a
diagnostics-instrumented **profile run** (attribution only — never pooled
with performance trials). Anything beyond these two needs a fresh Ryan go.

Design + validation record: `fgs_dagedge_design_20260924.md`. Gate-0:
`fgs_lshortening_gate0_20260924.md` (edge split = 4.25× sweep bound at
θ=4KB, flat in j; expected end-to-end at R4 j64: 3.24 → ~1.7–1.8 s).
Champion baseline: Stage 2 (`fgs_scalability_stage2_results_20260924.md`)
dagteam `dagteam_idle=:backoff` @ j=64 = 3.24 s/solve.

## Campaign pins (annotated tag `campaign/p021-fgs-dagedge-20260924` in all three repos)

| repo | branch | pinned commit | notes |
|---|---|---|---|
| FLOWPanel.jl | `fastmultipole` | the commit carrying this file (tag target; SHA recorded in the deployment section below), on top of `83b8482` (:dagedge plumbing + gate-0 artifacts) | adds the dagedge benchmark harness (below) |
| FastMultipole | `flowpanel-20260817` | `90a60cc3` (:dagedge implementation on top of Stage-2 pin `053c8de7`) | `fgs_dagedge_test.jl` 200/200; dagteam gate-1 28/28 |
| FLOWVPM.jl | `flowpanel` | `8d4a3b4` (unchanged Stage-1/2 pin; new-merge-law default) | untouched; pinned to complete the triple |

## Harness (new in the FLOWPanel pin)

- `benchmark/fgs_cold_common.jl`: `sweep_order="dagedge"` accepted;
  `dagedge_theta` config key (dagedge only, ≥ 0); forwarded by `cold_make`.
- Driver `benchmark/fgs_r4_dagedge.jl`: one process = one (j, placement,
  idle) point running a **batch** of executor arms (`DAGEDGE_ARMS`, e.g.
  `dagteam,dagedge:4096`) — cold = zero-initial-guess solves, arms batched
  per (j, placement) process (Ryan 2026-09-23). Per arm: fresh construction,
  compile + warmup + `DAGEDGE_SOLVES` trials, residual-history replay; every
  accepted solve passes the independent evaluator (`certified_fmm`, f32full
  certified by the evaluator, never bit-compared), runs exactly 27 outer
  iterations (fixed 27×3 work, tolerance=0), and repeats to 1e-8.
  Cross-executor agreement (dagteam vs dagedge, mathematically equal but not
  bitwise) is recorded per arm with a 1e-5 divergence tripwire (R2 f32full
  smoke measured 3.4e-7; the evaluator is the binding accuracy gate). Free schedule metadata
  logged at construction: `plan.ntasks`, big/small/rootfin/back static-task
  counts, `sim_makespan`, `sim_edge_L` (bytes).
- Launcher `benchmark/run_r4_fgs_dagedge.slurm.sh` (m12-class, 128c zen3
  exclusive, 500G, 12 h, champion placement socket 0 interleave 0–3, per-j
  calibrated configs from the verified 13777133 ladder, resume via
  `RESUME_FROM_JOB_ID`):
  - `RUN_MODE=perf` (submission 1): `primary-p{1..3}` — 3 paired A/B blocks
    @ j=64 backoff, dagteam vs dagedge(θ=4KB), in-process arm order
    alternating across blocks; `ladder-j{16,32}` — both executors (j=64 =
    primary); `theta-j64` — dagedge θ ∈ {0, 4096, 16384} in one process
    (θ is a plan-build knob; j/placement stay process-level).
  - `RUN_MODE=profile` (submission 2): `prof-j64-b{1,2}` + `prof-j{16,32}`,
    the same arms instrumented (`DAGEDGE_DIAG=1`; `:dagteam_*` keys, which
    dagedge reuses; lockmgmt identically 0 for dagedge).

Decision thresholds (gate-0 + design doc): sweep ~2.0 → ~0.5 s/solve,
end-to-end 3.24 → ~1.7–1.8 s at j=64. **< ~3× sweep gain ⇒ scheduler
overhead eating the prize** — exactly what the profile run must attribute
(busy/idle/wait/reduce split; per-task busy vs the R2-extrapolated
~0.3 µs/task; busy_max/min imbalance — static-schedule timing skew is the
known risk, recorded fallback = work-stealing deques, next optimization =
NUMA first-touch repack of split Lmat blocks).

## Local verification (2026-09-24, laptop, ≤4 threads, julia 1.12.5)

- `cold_check_config` dagedge validation unit checks: PASS (accept/reject
  matrix incl. θ=0, θ type, dagedge_theta-requires-dagedge).
- 4-thread R2 driver smoke (arms `dagteam,dagedge:4096,dagedge:0`, backoff,
  2 trials/arm): PASS (exit 0, status=completed). All gates green
  (certified_fmm, 27 iterations, 1e-8 repeat, history replay exact);
  schedule metadata populated — dagteam ntasks=655; dagedge θ=4096:
  ntasks=5829 (big 4867 / small 307 / rootfin 1 / back 327, sim_edge_L
  9.69 MB, sim_makespan 55.2 MB); θ=0: ntasks=11278 (big 10623 / small 0),
  sim_edge_L 7.92 MB. Cross-executor delta 3.34e-7/3.38e-7 (f32full
  reduce-order rounding, within the 1e-5 tripwire).
- `:dagedge` unit plumbing (`test/runtests_unit_solver.jl`) and FastMultipole
  `fgs_dagedge_test.jl` 200/200: green at the pinned commits (2026-09-24,
  pre-existing record in `fgs_dagedge_design_20260924.md`).
- `test/runtests_benchmark_cold.jl`: pre-existing environmental failure on
  this laptop (BLAS pin) — not chased per reset prompt.

## Deployment (fill at submission)

- Origin push still deferred (GitHub re-auth owed; `gh auth status` invalid
  token 2026-09-24); `deployment = "rsync"` mode per Ryan's 2026-09-22
  ruling — tagged content shipped via `git archive <tag> | ssh orc tar -x`
  into fresh dirs under `/home/rander39/campaigns/p021-fgs-dagedge-20260924/`,
  sha256 manifests verified on orc, campaign env with Manifest dev-paths at
  the deploy trees.
- [x] Tags created in all three repos (annotated, verified `tag` objects):
      FLOWPanel `a418f8d6a3456cbb65e4d65c168f2e8d00e0be48` (harness commit,
      on top of `83b8482`), FastMultipole
      `90a60cc3576dfcfccfe1db37d1c6eee9e62f2c76`, FLOWVPM
      `8d4a3b4d3012c42fc7d078629234c105b1e570f7`.
- [x] Tagged triple deployed 2026-09-24 via `git archive <tag> | ssh orc
      tar -x` into fresh `/home/rander39/campaigns/p021-fgs-dagedge-20260924/`
      (ARCHIVER_SKIP marked). Manifests generated from the same local export
      and `sha256sum --quiet -c` VERIFIED on orc in all three trees; R4 mesh
      family `dji9443_20260813_*_captess4.msh` confirmed present (24 files).
      | manifest | files | sha256 |
      |---|---|---|
      | MANIFEST.FLOWPanel.jl.sha256 | 3219 | `87910109574c7c1045e113c3092d3b58568b915be49ae12a42c766a41e96198b` |
      | MANIFEST.FastMultipole.sha256 | 6480 | `76ccfda80a408f9c4030e9bee3f0f9f2c19f1e456709d54b8979c40d08c2c3d6` |
      | MANIFEST.FLOWVPM.jl.sha256 | 188 | `876441cb17527d06df61daa849d4b10f125ea83e18f41c21b2e6af77b9f355fa` |
- [x] Campaign env built (julia/1.11.7-6bmogfl; Project.toml from the
      Stage-2 env, `Pkg.develop` on the three deploy trees + instantiate,
      111 deps precompiled); Manifest dev-paths verified to resolve to the
      deploy trees. `pins.toml` written with `deployment = "rsync"` +
      manifests per the table above.
- [x] Availability probed 2026-09-24T05:37Z: m12 access=normal, 0 idle /
      22 mixed / 106 alloc of 136; `sbatch --test-only` (128c/500G/zen3/
      exclusive/12h) estimated start 2026-09-24T11:59 on m12-2-3 —
      conservative; the 12 h wall backfills well.
- [x] **Performance run submitted 2026-09-24 (job 13879622, m12, PENDING at
      submission)** from the deployed FLOWPanel tree top level (`logs/slurm/`
      pre-created): `sbatch -p m12 --export=ALL,RUN_MODE=perf,
      COLD_PROJECT=.../env,CAMPAIGN_PINS=.../pins.toml,
      COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910
      benchmark/run_r4_fgs_dagedge.slurm.sh`. Run dir:
      `data/p021-cold-20260910/fgs-dagedge-perf-13879622`. Resume: resubmit
      with `RESUME_FROM_JOB_ID=13879622` in the same `--export` list;
      STATUS_*=ok stages skip.
- [x] **Profile run submitted 2026-09-24 (job 13879625, m12,
      `--dependency=afterany:13879622`)** — same submit line with
      `RUN_MODE=profile`; runs only after the performance job ends (keeps
      the exclusive-node footprint serial and rankings uncontaminated). Run
      dir: `data/p021-cold-20260910/fgs-dagedge-profile-13879625`.

## Measurement caveats (binding on analysis)

- Stage-1/2 caveats all carry: fixed-work gates, diag −1 sentinels, NEVER
  pool instrumented and uninstrumented rows (Stage 2 measured
  instrumentation shifting timing by 8.5%), plateau vs regression separate.
- Rankings come only from the performance run; the profile run only
  explains them.
- Worker durations overlap — report shares of team-time, never sums as
  elapsed.
- dagedge vs dagteam solutions are mathematically equal, not bitwise; the
  binding accuracy gate is the independent evaluator.
- `sim_makespan`/`sim_edge_L` are build-time byte-cost simulation outputs,
  not measurements — use for schedule sanity, not timing claims.
- θ probe arms share a process with fresh constructions; construction-order
  effects are possible in principle — the primary A/B pairs (alternating
  order) are the ranking-grade comparison.
