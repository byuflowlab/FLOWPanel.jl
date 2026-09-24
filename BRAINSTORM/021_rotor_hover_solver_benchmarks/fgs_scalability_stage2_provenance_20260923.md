# 021 FGS scalability — Stage 2 provenance (2026-09-23)

Plan: `fgs_scalability_diagnostic_plan_20260921c.md` (plan C) Stage 2, merged
with plan-C Stage 3's per-worker aggregate design per Ryan's 2026-09-23
directive ("prep stage 2 as you suggested. You have my approval to launch when
you're ready. Try to make this test lower-cost than the last one"). Stage 1
results: `fgs_scalability_stage1_results_20260923.md` (job 13858983 — both
effects reproduced, all loss in `nonself_product`; cap16@j64 best at 3.41 s).
HPC submission pre-approved by Ryan for this Stage 2 job (resume path
included).

## What Stage 2 measures

1. **Per-worker drain-loop aggregates** (new, FastMultipole `DagWorkerStats`):
   per worker — lower/back task counts, busy task-body time, busy-path
   lock+pop time, empty-pop idle-streak time. Padded structs, no shared
   counters, no per-spin timing (idle streaks timed at their boundaries only;
   idle-loop lock churn deliberately folded into idle, not lockmgmt).
   Discriminator: per-task busy inflating with j at fixed task count →
   bandwidth/locality; flat busy with growing idle+lockmgmt shares →
   scheduling (narrow DAG / lock contention).
2. **Bounded-backoff idle policy A/B** (`dagteam_idle=:backoff`): starved
   workers pause in a bounded exponential loop gated on a `qhint` atomic
   queue-length hint (maintained under the queue lock) instead of hammering
   the SpinLock. Recovery isolates the idle-lock-hammering contribution.

Cost cuts vs Stage 1 (45 stages, 48 h wall, 12.8 h used): no j=1 rung, no
placement A/B, no cap32, no accepted bridge, 1–2 blocks; 17 stages
(+4 conditional j=32 idle pairs), 12 h wall. Champion placement only.

## Campaign pins (annotated tag `campaign/p021-fgs-stage2-20260923` in all three repos)

| repo | branch | pinned commit | notes |
|---|---|---|---|
| FLOWPanel.jl | `fastmultipole` | the commit carrying this file (tag target; sha recorded in the deployment section below) | Stage-2 commit on top of `a835265` |
| FastMultipole | `flowpanel-20260817` | `053c8de7` on top of `649405cb` (Stage-1 pin) | DagWorkerStats + qhint + `dagteam_idle` |
| FLOWVPM.jl | `flowpanel` | `8d4a3b4` (unchanged Stage-1 pin) | untouched; pinned to complete the triple |

FLOWPanel Stage-2 content: `dagteam_idle` plumbing in
`FGSSolver`/metadata/`fgs_cold_common.jl`, driver `DAGTEAM_IDLE` env + 10 new
`diag_dagteam_*` columns (`benchmark/fgs_r4_dagteam_stage1.jl`, reused as the
Stage-2 driver), launcher `benchmark/run_r4_fgs_stage2.slurm.sh`, analysis
`benchmark/fgs_stage2_analysis.jl` (parse-checked only; first execution
against the harvested CSVs, as with Stage 1).

## Local verification (2026-09-23, laptop, ≤4 threads)

- `test/runtests_unit_solver.jl`: all testsets pass (exit 0).
- `test/runtests_r4_dagteam_ab_driver.jl`: PASS.
- 4-thread Stage-2 smoke (scratchpad `stage2_smoke.jl`): 46/46 + 4/4 —
  solutions **bit-identical** across idle policies (:spin/:backoff) and
  worker caps (0/1/2), dagteam matches lexicographic to 1e-6; the 10 new
  aggregate keys exist and are zero for :lexicographic; self-consistency
  holds (`n_lower == leaf_visit_count`, worker-time observed ≤ team_size ×
  nonself elapsed, coarse dagteam timers ≤ nonself); uninstrumented solves
  leave stats untouched; invalid idle policy rejected; metadata exposes
  `dagteam_idle`; `cold_check_config` accepts spin/backoff on dagteam and
  rejects invalid values and non-dagteam use.
- `test/runtests_benchmark_cold.jl`: pre-existing environmental failure on
  this laptop (BLAS pin reports 8 threads in single mode) — fails identically
  on the unmodified tree (verified via stash); harness targets the cluster
  environment.

## Stage list (launcher `benchmark/run_r4_fgs_stage2.slurm.sh`)

All champion placement (socket 0, interleave 0-3), fixed-work arm
(27×3/tolerance=0), calibrated configs from the verified 13777133 ladder
(j ∈ {16,32,64} only):

- A `aggdiag-*`: diag=1, spin — j16-b1, j32-b1, j64-b1/b2 (w=0);
  `aggdiag-cap-j64-b1/b2` (w=16). [6]
- B `anchor-*`: diag=0, spin — j16-b1, j32-b1 (w=0); `anchor-cap-j64-b1`
  (w=16). j64 w=0 anchors are section C's spin arms. [3]
- C `idle-{spin,backoff}-j64-p{1..3}`: paired A/B, alternating arm order,
  w=0, diag=0. [6]
- D safety: `idle-backoff-j16-b1` (w=0), `idle-backoff-cap-j64-b1` (w=16). [2]
- E conditional (`backoff_verdict.txt`): `idle-{spin,backoff}-j32-p{1,2}`
  only if backoff's paired median beats spin at j=64. [0 or 4]

## Deployment (fill at submission)

- Origin push still deferred (GitHub re-auth owed); `deployment = "rsync"`
  mode again per Ryan's 2026-09-22 ruling — tagged content shipped via
  `git archive <tag> | ssh orc tar -x` into fresh dirs under
  `/home/rander39/campaigns/p021-fgs-stage2-20260923/`, sha256 manifests
  verified on orc, campaign env with Manifest dev-paths at the deploy trees.
- [ ] Tags created in all three repos
- [ ] Tagged triple deployed + manifests verified
- [ ] Campaign env built (COLD_PROJECT), pins.toml written
- [ ] zen3 availability confirmed (slurm-availability, 128c/500G/12h --eta)
- [ ] Submitted (job id, run dir)

## Measurement caveats (binding on analysis)

- Stage-1 caveats all carry (fixed-work gates, diag −1 sentinels, never pool
  instrumented/uninstrumented, plateau vs regression separate).
- Worker durations overlap — never sum them as elapsed time; report shares
  of team-time (team_size × nonself elapsed).
- `busy_max/min_ns` accumulate per-outer-iteration extrema (sums of
  per-block extrema) — imbalance indicators, not single-sweep extrema.
- `lockmgmt_ns` covers busy-path lock+pop only; idle-streak lock churn is in
  `idle_ns` by design (timing every spin would perturb the contention under
  study).
- An empty ready queue has two readings (dependency starvation vs all work
  already running); the aggregates alone do not separate them — that is
  plan-C Stage 3's sampled-timeline territory if it matters.
- The spin baseline now carries one extra atomic store per locked queue
  operation (`qhint` maintenance) vs Stage 1's binary — cross-job absolute
  comparisons should use this job's own anchors; paired A/Bs are same-binary
  by construction.
