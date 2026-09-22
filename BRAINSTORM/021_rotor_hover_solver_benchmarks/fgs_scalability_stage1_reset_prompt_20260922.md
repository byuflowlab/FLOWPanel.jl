# Reset prompt: 021 FGS scalability — Stage 1 execution (2026-09-22, supersedes fgs_scalability_reset_prompt_20260922.md)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

BRAINSTORM 021 is CLOSED AT PROMOTION; the active follow-on is the **FGS
scalability diagnostic**, plan = `fgs_scalability_diagnostic_plan_20260921c.md`
(plan C — self-contained; execute from it). Read first: `CLAUDE.md`, then
`agent_policies/WORKFLOW.md`/`TESTING.md`/`HPC.md` before corresponding work.

**Stage 0 is COMPLETE** — full desk audit in
`fgs_scalability_stage0_audit_20260922.md` (read it in full; it is the
evidence base, baseline manifest, hypothesis table, and exit checklist).
Both of plan C's optional patches were folded in there. Do not redo Stage 0.

Stage-0 highlights (details + file:line pointers in the audit file):

- Baseline VERIFIED: reported dagteam ladder 33.07/5.97/4.50/4.42/6.42 s at
  j=1/8/16/32/64 = dagteam-order medians (40 accepted rows/j) in
  `thread-scaling-j<J>-13777133/fgs-trials/results/ab_trials.csv`, data root
  `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/`. Hardware:
  orc-m12, 2× EPYC 7763, numactl interleave 0–3, physcpubind 0–63, BLAS=1.
- Loss targets: plateau shortfall 2.17 s @ j32 (vs ideal doubling from 16);
  regression 2.00 s absolute (4.21 s vs ideal) @ j64.
- Structural facts: worker team spawn/join EVERY outer iteration (27×/solve,
  inside timing); busy-spin (jl_cpu_pause, never yield) in 3 places; serial
  per-sweep boundary reduction (`dagteam_reduce_u!`); ready-queue pop is a
  linear scan under a SpinLock. H2 (sync/scheduling overhead) leading; H4
  (NUMA/bandwidth) live; H7 demoted (work exactly 27×81 at every j); H5 GC
  branch demoted (median gc=0) but allocations grow ~2.1× with j.
- No worker-cap knob exists (`nw = Threads.nthreads()` hard-coded in
  `build_dagteam_plan`, FastMultipole `src/solve_dagteam.jl:249`); adding one
  is low deadlock-risk (atomic task counter, no fixed-arrival barrier).
- Reusable instrumentation: `FastMultipole.solve!` coarse-phase
  `diagnostics` dict (solve.jl:1348-1357) covers init/fmm/influence/
  residual/leaf-solve/nonself-product/scatter + counts; GAP: spawn/join is
  buried inside `nonself_product_ns`. `stage_observer` hook exists. Harness
  = `benchmark/fgs_r4_dagteam_ab.jl` (AB_MODE=calibrate|trials|activity) +
  `benchmark/fgs_cold_common.jl` per-solve schema (`cold_trial`, :495-513).
  Fixed 27×3 workload = existing kwargs `max_iterations=27,
  inner_iterations=3, tolerance=0.0`.

## NEXT TASK: Stage 1 — reproduce and decompose (HPC submission APPROVED)

**Ryan has pre-approved HPC job submission for this stage** (2026-09-22).
Follow plan C's Stage 1 + measurement contract exactly. Work items:

1. Two small FastMultipole patches (commit before pinning):
   a. Expose `nw` as a kwarg through `build_dagteam_plan` →
      `FastGaussSeidel`/`FGSSolver` (worker-cap knob; capped-out workers are
      never spawned — see audit §2 for why this can't deadlock).
   b. Split team spawn/join (and serial `dagteam_reduce_u!`) out of
      `nonself_product_ns` into new diagnostics keys, keeping existing keys'
      meaning intact.
2. Local smoke at ≤4 threads (TESTING.md narrow FGS/solver checks +
   multi-thread progress check for the worker cap). Local runs NEVER >4
   threads.
3. Campaign setup per house rules: annotated tags
   (`campaign/p021-fgs-stage1-20260922` convention) in BOTH repos
   (FLOWPanel + FastMultipole), worktrees from the tags, Manifest dev-paths
   at the worktrees, pins recorded in a provenance file in BRAINSTORM/021
   BEFORE submitting. Outputs to the consolidated data root.
4. One exclusive zen3-node allocation (orc-m12-class; check availability):
   - Fixed-work ladder (27×3, tolerance=0) at j=1/16/32/64 (+8 only if the
     knee needs it), 3 fresh-process blocks × 5 warmed solves per config,
     thread order randomized within block; accepted-accuracy arm at
     16/32/64 (existing calibration; Stage-0 evidence may substitute for
     screening but not final certification of a winner).
   - Coarse-phase diagnostics in SEPARATE matched runs; verify ≤5% overhead
     and unchanged scaling shape before attributing anything.
   - Within-allocation A/Bs: 64-thread placement A/B (champion interleave
     vs one explicit alternative, fresh first-touch under each) and
     worker-cap A/B (cap=16 @ j=64 on the matching nested CPU set; cap=32
     only if 16 helps).
   - Record the full measurement-contract manifest (topology, affinity,
     page-placement verification, pins, per-solve rows per plan C's schema;
     ship the analysis script with the data).
5. Decompose: T_phase(32)−T_phase(16), T_phase(64)−T_phase(32), and
   shortfall-from-doubling per phase; reconcile phase sums vs totals.
   Plateau and regression are SEPARATE conclusions.

Gate (plan C): if either effect fails to reproduce, stop and report the
comparison deltas — do not proceed to interventions on a non-reproduced
effect. Stage 2 intervention selection follows from the decomposition.

## Owed (small, this session or next)

- **R2-j1 (13829232_7, dir `r12-champion-r2-j1-13778533/`)**: still RUNNING
  at 2026-09-22 (~17 h of 48 h wall; p2 tuning mid-flight). Re-check via
  `hpc-monitor`; when `phase2/phase2.csv` lands, harvest via `harvester`
  and fill the R2 j1 column.
- R1-j32 is DONE and harvested (audit file §5) — the "gap" was a
  case-sensitive dir name (`r12-champion-R1-j32-13778533`, capital R1).
  If the merge script gets fixed, fix its glob.

## Standing Ryan-gated ledger (carried forward, unchanged)

- Production-adoption question at R4 (dagteam vs nfcache, speed-vs-memory;
  what's prunable above 32 threads) — present, let Ryan rule.
- Origin pushes (both repos + v21/v22/v23 + campaign tags); notebook
  entries owed (offer, don't write); WeakKeyDict/warmstart fix
  (`_publish_block_gs_status!` from `7fbd68a`); `:dagteam` unit test
  missing in `test/runtests_unit_solver.jl`; 3 RECENT p018 runs awaiting
  archive approval; hpc-storage archive-pass report.
- Optional reruns: R4 ILU j=1 longer candidate cap; R3+; R1–R2 f32full;
  zen3 BLAS A/B rider.

## Traps (all prior traps still bind)

- nfcache-vs-FGS comparisons must state cache-build exclusion (~94 s @ j8)
  and memory footprint (~8.5 GB cache + ~9.8 GB state) EVERY time.
- `bc_certified` in phase2.csv = the certified-FMM evaluator's
  self-certification (`bc_fmm(x).error_success`), NOT the ≤1e-6 test; R1–R2
  rows can be uncertified-flag yet have rel_l2 < 1e-6 (R1-j32 is exactly
  this). Never discard rows on that flag without checking
  `bc_rel_l2_certified`.
- Run-dir names are case-sensitive and inconsistent (`r2-j1` vs `R1-j32`).
- Task logs are OUTPUT-BUFFERED (can freeze ~10 h while computing) — judge
  liveness by CPU/ps/output files, never log mtime. Judge runs by outputs,
  never sacct.
- Resume run dirs carry OLD job ids; first-pass logs in
  `logs.before.<newid>/` — never double-count.
- BLAS conventions: R4 thread-scaling = BLAS 1; R1–R2 = BLAS j.
- `ssh orc` needs a live ControlMaster socket (2FA otherwise); Slurm CLI
  needs a login shell; filter login banners.
- Do not edit source while a job is using its worktree; each campaign in
  its own tag-pinned worktree, never a shared live checkout.

## House rules (binding)

Monitoring via `hpc-monitor`; harvesting via `harvester`; storage via
`hpc-storage` (400 G cap); notebook writes Ryan-gated (offer, don't write);
dated status/provenance files in BRAINSTORM/021. HPC submission is
approved FOR STAGE 1 ONLY per above — anything beyond (Stage 2/3 runs,
reruns from the optional ledger) remains Ryan-gated.
