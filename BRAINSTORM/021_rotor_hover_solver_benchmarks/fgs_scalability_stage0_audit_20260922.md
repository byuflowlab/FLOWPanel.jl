# Stage 0 desk audit — FGS scalability diagnostic (2026-09-22)

Executes Stage 0 of `fgs_scalability_diagnostic_plan_20260921c.md` (plan C).
Both optional patches folded in: (1) R2-j64 blow-up added to the evidence
table as H2 prior weight; (2) nfcache rows cited by their certified paths
(no re-verification).

## 1. Baseline manifest (raw-record verification)

Reported R4 dagteam ladder (median accepted solve time, s):
33.07 / 5.97 / 4.50 / 4.42 / 6.42 at j = 1/8/16/32/64.

Raw-record verification (harvester, 2026-09-22), source
`<data root>/thread-scaling-j<J>-13777133/fgs-trials/results/ab_trials.csv`
(data root `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/`),
80 rows per j, schema = `fgs_cold_common.jl` cold_trial (30 cols incl.
solve/gc seconds, allocated/retained bytes, peak RSS, iterations, estimated
sweeps/fmm passes, fmm_certified, authoritative evaluator + accepted flag):

Dagteam rows only (sweep_order=dagteam, 40 rows/j, all accepted/eligible/
uninstrumented):

| j | reported | dagteam median | dagteam min | med alloc MB | med gc_s | med iters | med inner sweeps |
|---|---|---|---|---|---|---|---|
| 1 | 33.07 | 33.066 | 33.035 | 60.7 | 0.000 | 27 | 81 |
| 8 | 5.97 | 5.969 | 5.888 | 68.5 | 0.000 | 27 | 81 |
| 16 | 4.50 | 4.502 | 4.455 | 76.9 | 0.000 | 27 | 81 |
| 32 | 4.42 | 4.424 | 4.367 | 93.8 | 0.000 | 27 | 81 |
| 64 | 6.42 | 6.422 | 6.326 | 127.6 | 0.000 | 27 | 81 |

- **Ladder VERIFIED**: the all-row medians above mix dagteam and colored
  trials (40 rows each per j, interleaved in one CSV, all 80 accepted/
  eligible/uninstrumented). Grouped by sweep_order, the dagteam medians are
  33.066 / 5.969 / 4.502 / 4.424 / 6.422 s — matching the reported ladder
  to rounding at every j. Colored medians: 35.58 / 9.07 / 7.26 / 6.91 /
  12.30 (also matches the reported colored column). Same solver variant and
  timing boundary confirmed (single schema, single harness).
- Facts for hypotheses: iterations EXACTLY 27 (81 inner sweeps = 27×3) at
  every j — the production accepted runs already execute the plan-C fixed
  workload, and work counts are thread-invariant (demotes H7 for these
  rows; accepted-arm calibration equivalence still to confirm in Stage 1);
  GC time identically 0.000 in dagteam medians (demotes H5's GC branch as
  a median effect; tail pauses unexcluded); **allocations GROW ~2.1× with
  j** (60.7 → 127.6 MB, j1→j64) — consistent with per-outer-iteration team
  spawn allocating per worker (feeds H2).
- Hardware/placement (all j): node orc-m12, AMD EPYC 7763 64-core ×2 (128
  logical), numactl interleave 0–3, physcpubind 0–63 (socket 0), BLAS=1
  (j1 BLAS unknown). Unknowns: frequency data, page placement verification,
  per-phase times (none recorded).
- j=1 krylov: `ilu/phase2.csv` absent (only tune_phase2.csv) — matches the
  known j1 ILU FAILURE; R4 ILU@j1 stays missing.
- Resume caveat honored: j≥8 first-pass logs in `logs.before.13829231/`.

nfcache context rows (cited, already certified — plan-C patch 2): 4.25 /
3.56 / 2.73 / 2.41 s at j = 8/16/32/64, bc_certified=true @ 1e-6, cold
isolated solves on warm cache, budget-0 knobs P=15/MAC=0.55/leaf=32/32/32/21.
Caveats bind on every citation: excludes ~94 s one-time cache build
(`nfcache_build_time`), ~8.5 GB cache + ~9.8 GB solver state.

Citation paths (bc_certified=true confirmed for all nfcache rows):
`<data root>/thread-scaling-j{8,16,32,64}-13777133/ilu/phase2.csv`.

## 2. Code audit — dagteam executor and harness (code-scout, 2026-09-22)

### Timed solve boundary (plan-C Q: what does the timed solve include?)

- Harness times exactly `pnl._solve!(rotor, solver; diagnostics)` via
  `@timed` (`benchmark/fgs_cold_common.jl:499`); FLOWPanel kernel
  `_solve!(body, ::FGSSolver)` at `src/FLOWPanel_solver.jl:1773-1855` calls
  `FastMultipole.solve!` (FastMultipole `src/solve.jl:1341-1616`).
- INSIDE the timed region: `dagteam_initialize!` (solve.jl:1432), full outer
  loop with per-iteration `fmm!` (1461-1471, buffer resets only — trees are
  built at construction, outside timing), `residual!` (1489), convergence
  check (1506-1511), and `dagteam_inner_sweeps!` (1520-1533).
- **Worker team spawn/join happens EVERY outer iteration** inside
  `dagteam_inner_sweeps!` (`dagteam_start_team!`/`dagteam_stop_team!`,
  solve_dagteam.jl:486-493) — 27 spawn/join cycles per solve, all timed.
- `final_update` second `fmm!` is EXCLUDED (`final_update=false`,
  FLOWPanel_solver.jl:1830).
- Fixed 27×3 workload is directly expressible: `max_iterations=27,
  inner_iterations=3, tolerance=0.0` (recipe already used by
  `FGSPreconditioner`, FLOWPanel_solver.jl:1900-1926, and the A/B harness
  calibration stage, `benchmark/fgs_r4_dagteam_ab.jl:65,70`).

### Worker-cap knob

- **Does not exist**: `build_dagteam_plan` hard-codes
  `nw = Threads.nthreads()` (solve_dagteam.jl:249);
  `dagteam_start_team!` spawns `length(plan.xg)-1` workers (400-406).
- Adding one is low deadlock-risk: completion is an atomic task counter
  (`plan.ndone[] < plan.ntasks`) over a shared queue, not a fixed-arrival
  barrier — fewer workers drain the same queue; capped-out workers are
  simply never spawned (no polling stragglers). Requires exposing `nw` as a
  kwarg. Note: architecture is strictly sequential (team torn down before
  each `fmm!`), so a cap does NOT enable FMM/nearfield concurrency.

### Executor structure (H-relevant mechanics)

- Plan build splits per-leaf nonself operator into lower (`Lmat`,`preds`)
  and upper (`Umat`,`uppers`) triangles; in-degrees `indeg0`, `roots`,
  byte-weighted critical-path priority `prio` (solve_dagteam.jl:96-272).
- Per sweep: coordinator seeds `readyQ` under a `Threads.SpinLock`
  (`qlock`), bumps an epoch to wake workers, and itself drains.
- `dagteam_drain!` (354-377): pop highest-priority ready leaf (LINEAR SCAN
  of `readyQ`, 343-352) or a `backQ` upper-product task; if both empty,
  **busy-spin** (`jl_cpu_pause`), never yield — intentional (comment 38-40).
- Busy-polling in THREE places: idle drain branch (367-368), coordinator
  laggard wait (427-429), inter-sweep worker wait in `dagteam_worker!`
  (392-393).
- `dagteam_reduce_u!` (326-339): **serial** boundary reduction by the
  coordinator, fixed ascending source order, after all tasks finish.

### Reusable instrumentation (plan-C: reuse before building)

- `FastMultipole.solve!` built-in coarse-phase `diagnostics` dict
  (solve.jl:1348-1357): `total_ns, initialization_ns, fmm_ns,
  influence_mapping_ns, residual_ns, leaf_solve_ns, nonself_product_ns,
  scatter_ns, remaining_iteration_ns, final_update_ns, outer_count,
  sweep_count, leaf_visit_count`. Matches the Stage-1 coarse-timer list
  almost exactly. **Gap: team spawn/join is buried inside
  `nonself_product_ns`** (wraps all of `dagteam_inner_sweeps!`) — not
  separately broken out.
- `stage_observer(stage, :start/:stop)` hook exists (used by
  `activity_observer`, `benchmark/fgs_r4_counters.jl:24`).
- Harness: `benchmark/fgs_r4_dagteam_ab.jl` is the dagteam-specific driver
  (env: `AB_MODE=calibrate|trials|activity`, `DAGTEAM_PRECISION` default
  f32full, `COLD_AB_REPS`/`COLD_AB_BATCHES`, config TOML paths).
  `AB_MODE=activity` = existing diagnostic mode (stage-observer solve,
  per-thread activity CSV + validation CSV).
- Per-solve record schema already exists in `cold_trial`
  (`benchmark/fgs_cold_common.jl:495-513`): setup/solve/total seconds,
  allocated bytes, gc seconds, retained bytes, peak RSS, iterations,
  estimated inner sweeps / fmm passes, solved/eligible + validation fields
  (`authoritative_rel_l2`, `accepted`).
- `benchmark/rotor_hover_solver_phase2.jl` (note: benchmark/, not
  examples/) is the solver-family harness (t_solve_min, bc_certified, …,
  resume support) — the thread-scaling data source, but not the
  dagteam-diagnostic driver.

### DAG availability

Materialized data structure, no analyzer needed: `solver.fgs.dagteam ::
DagTeamPlan` (FastMultipole `containers.jl:1131-1160`, field at 1206) with
`preds/nsucc/uppers/indeg0/roots/red/prio`; `ntasks = n_leaves +
count(>0, mup)` (solve_dagteam.jl:266). Level/width summaries = one quick
BFS in a REPL.

### Relaxation / block-Jacobi readiness

`FGSSolver.rlx` exists (FLOWPanel_solver.jl:1485,1504; applied once per
OUTER iteration at solve.jl:1571, not per inner sweep). Previous-sweep
state hooks exist (`strengths_old` solve.jl:1517; `plan.x`/`plan.u`/
`plan.lsum`), but an inner-sweep-granularity Jacobi read path would need
new plumbing in `dagteam_sweep!`/`dagteam_do_lower!`.

## 3. Hypothesis / evidence table

Status: fact = observed in records/code; prior = compatible evidence, no
controlled intervention yet. Nothing is causal from static inspection.

| ID | Candidate explanation | Current evidence (all prior-level) |
|---|---|---|
| H1 | Dependency width / load imbalance | DAG materialized but width unmeasured; `prio` critical path exists. Unresolved. |
| H2 | Queue/scheduling/polling/sync overhead grows with threads | FACTS from code: 27× per-solve team spawn/join inside timing; 3 busy-spin sites (never yield); linear-scan ready-queue pop under a SpinLock. PRIOR from data: R2-j64 blow-up 1.20→10.7 s j32→j64 (~9×), worse on the smaller rung = less work per sync (`champion_adoption_reset_prompt_20260921.md:40,108`); both FGS variants regress at 64 (0.69×/0.56×) while krylov_ilu still scales in the same jobs. **Leading.** |
| H3 | Serial reduction/copies/outer stage | FACT: `dagteam_reduce_u!` is serial per sweep. Magnitude unmeasured. |
| H4 | Memory throughput / NUMA locality | PRIOR: placement is load-bearing (socket membind collapsed dagteam to 1.09×; champion = interleave 0-3). Krylov counter-evidence does NOT rule this out (different kernels/working sets). |
| H5 | Allocations / GC | FACT (dagteam rows): per-solve allocations grow ~2.1× with j (60.7→127.6 MB, j1→j64) but median GC time is identically 0 — allocation growth is real (points at per-iteration team spawn, feeds H2), GC cost negligible in medians; tail pauses unexcluded. |
| H6 | Workers interfere with runtime/OS/GC at full occupancy | PRIOR: j=64 uses all of socket-0 cores, no headroom; spinning workers never yield (couples with H2). |
| H7 | Iteration/calibration changes with threads | DEMOTED for the verified rows: iterations exactly 27 / 81 inner sweeps at every j in the dagteam accepted set — work counts are thread-invariant. Residual-history/calibration equivalence still to record in Stage 1. |
| H8 | SMT siblings / CPU placement semantics | Unresolved — champion CPU-set semantics vs node topology must be resolved before interpreting nominal 64. |
| H9 | Frequency/power/thermal | No data; APERF/MPERF optional per measurement contract. |

## 4. Smallest discriminating experiment (Stage-0 exit)

The Stage-1 baseline can be run with near-zero new engineering:

- Driver: `fgs_r4_dagteam_ab.jl` + `fgs_cold_common.jl` (existing per-solve
  schema, calibrate/trials/activity modes).
- Fixed-work mode: `max_iterations=27, inner_iterations=3, tolerance=0.0`
  (existing kwargs).
- Coarse phases: existing `diagnostics` dict; ONE small patch worth making
  before Stage 1: split team spawn/join (+ optionally serial
  `dagteam_reduce_u!`) out of `nonself_product_ns` so H2's most direct
  signature (27× spawn/join) and H3's (serial reduction) are separately
  visible. Both are wrap-a-timer changes in `solve_dagteam.jl`.
- Highest-value within-allocation A/Bs (already justified by priors):
  64-thread placement A/B (H4), and a worker-cap A/B — requires exposing
  `nw` in `build_dagteam_plan` (small kwarg change, low deadlock risk per
  audit) — cap=16 at j=64 is the single probe that most cleanly separates
  H2/H6 (spinning excess workers) from H1 (no ready work).

### Target quantities (from the verified dagteam medians)

- Plateau: 16→32 gains only 0.078 s (4.502→4.424, 1.7%); shortfall from
  ideal doubling = 4.424 − 4.502/2 = **2.17 s** at j=32.
- Regression: 32→64 LOSES 2.00 s (4.424→6.422, +45%); shortfall from ideal
  = 6.422 − 4.424/2 = **4.21 s** at j=64.
- Anchor: j=1 = 33.07 s → best observed speedup 7.5× at j=32.

### Stage-0 exit checklist (plan-C)

- [x] Baseline manifest: verified above (dagteam medians = reported ladder;
  orc-m12, EPYC 7763 ×2, interleave 0–3, physcpubind 0–63, BLAS=1, 27×3
  work at every j; unknowns recorded: frequency, page placement, per-phase
  times, j1 BLAS threads).
- [x] Reusable measurements/knobs list: §2 (diagnostics dict, stage
  observer, AB_MODE=activity, cold_trial schema, 27×3 kwargs; missing:
  worker-cap knob, spawn/join timer split).
- [x] Smallest discriminating experiment: worker-cap A/B at j=64 (cap=16)
  + placement A/B, on top of the existing-driver fixed-work ladder with
  the diagnostics dict enabled — separates H2/H6 (excess spinning workers)
  from H1 (no ready work) and H4 (placement) with two small code changes
  (expose `nw`; split spawn/join out of `nonself_product_ns`).
- [x] Harness entry points and row locations recorded (§1–2) — no
  rediscovery needed at execution handoff.
- [x] No hypothesis declared causal from static inspection.

## Next step (Ryan-gated)

Stage 1 = HPC submission: prepare the two small FastMultipole patches
(worker-cap kwarg, spawn/join + serial-reduction timer split), tag-pinned
worktrees for FLOWPanel + FastMultipole, launcher for one exclusive zen3
node running the fixed-work ladder (1/16/32/64, +8 only if knee needs it)
with 3 fresh-process blocks × 5 warmed solves, plus within-allocation
64-thread placement A/B and worker-cap A/B. Get Ryan's approval before
submitting.

## 5. Owed-item resolutions (carried from reset prompt)

- **R2-j1 (13829232_7)**: still RUNNING at 2026-09-22 (~17 h elapsed of
  48 h wall; `phase2/tune_phase2.csv` 5 rows, actively written;
  `phase2/phase2.csv` not yet created). Re-check later.
- **R1-j32**: NOT missing — `r12-champion-R1-j32-13778533/phase2/phase2.csv`
  exists with 16 data rows (completed 2026-09-20 07:57; note capital `R1` in
  the dir name — likely why the merge sweep's glob missed it). Harvested
  medians (t_solve_min, s):

  | config | t_solve (s) | rows | bc_rel_l2 |
  |---|---|---|---|
  | backslash_ldiv | 0.0201 | 1 | 5.4e-09 |
  | krylov_ilu_nfcache | 0.167 | 3 | 8.5e-07 |
  | fgmres_fgs_nfcache | 0.568 | 3 | 8.7e-07 |
  | fgs | 0.642 | 1 | 2.8e-07 |
  | fgmres_fgs | 0.880 | 1 | 8.7e-07 |
  | krylov_gmres_nfcache | 1.077 | 3 | 9.9e-07 |
  | krylov_ilu | 2.385 | 1 | 8.3e-07 |
  | krylov_jacobi | 6.557 | 1 | 4.5e-07 |
  | krylov_gmres | 19.73 | 1 | 9.7e-07 |

  (+1 `additivity` diagnostic row.) Caveat: every row has
  `bc_certified=false`. Semantics check (`rotor_hover_solver_phase2.jl:196,
  205`): that column is `bc_fmm(x).error_success` — whether the certified-FMM
  BC evaluator self-certified its error estimate — NOT the ≤1e-6 test. All
  rel-L2 values ARE below 1e-6; the evaluator simply failed to self-certify
  on this small rung. Usable with that caveat. Fits the cross-rung story:
  backslash_ldiv untouchable at R1, nfcache best iterative (0.167 s),
  FGS well behind at high j on small rungs.
