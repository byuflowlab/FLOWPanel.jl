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
- [ ] Tags created in all three repos (annotated, verified `tag` objects):
      FLOWPanel `FILL`, FastMultipole `FILL`, FLOWVPM `FILL`.
- [ ] Tagged triple deployed via `git archive <tag> | ssh orc tar -x`;
      manifests verified (`sha256sum --quiet -c`) on orc; R4 mesh present.
- [ ] Campaign env built (julia/1.11.7-6bmogfl; `Pkg.develop` on the three
      deploy trees + instantiate); `pins.toml` written.
- [ ] Availability probed; performance run submitted (job FILL).
- [ ] Profile run submitted after submission 1 (job FILL).

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
