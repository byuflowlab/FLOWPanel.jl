# Reset prompt: 021 FGS critical-path (L) shortening — design + Stage 2 wrap (2026-09-24, supersedes fgs_scalability_stage2_reset_prompt_20260923.md)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Mission (Ryan 2026-09-24)

The FGS scalability mechanism is RESOLVED (Stage 2). Ryan's directive:
**pursue where the factor lives — shortening L, the dagteam sweep's
critical path.** Design/prototype work first; any HPC benchmark run is a new
campaign and Ryan-gated. Read `CLAUDE.md` and `agent_policies/HPC.md` before
HPC work; local runs never >4 threads.

## What Stage 1+2 established (evidence, all verified)

Setting: `FGSSolver(sweep_order=:dagteam)` — FastMultipole
`src/solve_dagteam.jl` executes the near-field Gauss-Seidel sweep as a
pull-DAG over solve leaves (lower-triangle readiness; one aggregated GEMV +
cached leaf-LU per leaf; upper products as filler; serial boundary
reduction). R4 rotor fixture, zen3 128c exclusive, champion placement
(socket 0, interleave 0-3), fixed work 27×3, cold-START solves (zero
initial guess — Ryan's meaning of "cold"; `reset_cold!` does this).

- **Stage 1** (job 13858983, `fgs_scalability_stage1_results_20260923.md`):
  16→32 plateau (1.03× of 2.00) and 64-thread regression (1.35×) both
  reproduced; ALL loss in `nonself_product` (the sweep), which anti-scales
  2.0→2.7→4.5 s while FMM scales cleanly (2.0→1.3→1.0 s). Placement settled
  (both-socket +4.4 s). cap16@j64 = 3.41 s.
- **Stage 2** (job 13875511, `fgs_scalability_stage2_results_20260924.md`):
  mechanism = **idle-worker lock hammering** on the shared queue SpinLock.
  Bandwidth RULED OUT: per-task lower-GEMV 66 µs (j16) → 161 µs (j64
  uncapped) → **67 µs at j64 w=16** with all 64 threads present.
  **`dagteam_idle=:backoff` @ j=64 = 3.24 s, new best** (paired −2.86 s vs
  spin, all pairs same sign; no regression at j16 or under cap; solutions
  bit-identical — backoff = idle workers poll a lock-free `qhint` atomic
  with bounded exponential pauses instead of hammering the lock).
- **The wall that remains:** busy share of team-time is ~0.25 at 16 workers
  and ~0.06 at 64 → **average effective sweep parallelism ≈ 4** at every
  team size and policy. The sweep is critical-path-bound: T_sweep ≈ L
  (byte-weighted longest dependency chain), NOT W/p. No idle policy, worker
  cap, or load balancing can pass L. Plan-C bound: T_j ≥ max(W/j, L).

Pins: tag `campaign/p021-fgs-stage2-20260923` — FLOWPanel `90f7452`,
FastMultipole `053c8de7` (branch `flowpanel-20260817`; DagWorkerStats,
qhint, `dagteam_idle`), FLOWVPM `8d4a3b4`. Later local commits: `d5b6865`
(results+analysis fixes). Origin push still owed after GitHub re-auth.

## Static @threads scheduling — discussed with Ryan, assessed, parked

Ryan proposed cost-estimated static load balancing (`@threads`, one thread
per task when tasks < threads). Assessment (agreed in session 2026-09-24):

- Plain `@threads` over leaves is **incorrect** (breaks GS dependencies).
- Correct variants: (a) level-synchronous wavefront — easy but adds a
  barrier per level and kills cross-level overlap (fatal here: DAG levels
  are single-digit wide at R4); (b) **static list schedule** — per-thread
  task lists + per-task done-flags, cost weights are exact (bytes streamed,
  `n_i × ptot_i`, already used by `prio`), schedule computed once (DAG
  repeats 81×/solve).
- Expected performance ≈ backoff (both hit the L-bound); removes the
  remaining 2–6% lockmgmt and — the one real upside — enables **NUMA
  owner-local first-touch of Lmat/Umat blocks** (locality demonstrably
  matters: +4.4 s both-socket). Estimated single-digit-% gain on top of an
  L-bound. **Not the factor.** Pursue only as "backoff-equivalent + NUMA
  ownership" after L-shortening, measured against backoff as baseline.

## Where the factor lives: shortening L — ranked slate

L ≈ the weighted chain of {wait for last predecessor → aggregate pull GEMV
→ leaf LU solve → publish} along the longest path. With avg parallelism ~4
and sweep elapsed ~2.0 s (j16/cap floor), a 2–4× widening of the DAG is the
prize backoff cannot reach.

1. **Edge-level partial pulls (preserves the exact GS iterate — do first).**
   Today leaf i waits for ALL preds, then one aggregated GEMV. Instead:
   compute partial products y_ij = L_ij·x_j as INDEPENDENT tasks the moment
   pred j publishes; when all of leaf i's partials exist, reduce them **in
   fixed ascending-j order** (determinism preserved — same trick as the
   boundary reduction) and do the leaf LU. Tasks become edges, not nodes →
   DAG width multiplies; L shrinks toward chain of (LU + one edge-GEMV +
   cheap reduction). Costs: more, smaller GEMVs (lower flop efficiency —
   consider a min-bytes threshold: only split pulls above a size cutoff, or
   split only the critical-path leaves using `prio`); per-edge output
   buffers (bounded: Σ n_i over edges = the Lmat row count, one TS vector
   the size of... allocate y-space per edge like `q` does per upper block);
   more queue traffic (mitigated by backoff; or shard the queue). The repack
   already stores L in per-(i,j) column blocks — split points exist
   (`colofs`). Mathematically IDENTICAL iterate to :dagteam/:lexicographic
   — same certification, no recalibration.
2. **Leaf reordering to widen the lower DAG (changes the iterate — valid
   but must recalibrate).** The lower/upper split depends on the leaf
   permutation; lexicographic (Hilbert-ish spatial) order makes chains.
   `sweep_order=:colored` (multicolor GS, 021 Phase 2b) is the extreme:
   maximal width, worst convergence lag. Middle ground: critical-path-aware
   permutations (elimination-tree-style / nested-dissection flavor).
   Score on TIME-TO-ACCEPTED-ACCURACY, never fixed-27-iteration wall.
3. **Ryan's Jacobi-esque chunking — ALREADY IMPLEMENTED as
   `sweep_order=:chunked`** (021 v22; `chunks::Int` kwarg, default 64;
   `build_chunk_map` in solve.jl, `fgs_r4_chunked_ab.jl` benchmark;
   validated by `test/runtests_r4_chunked_ab_driver.jl`): Gauss-Seidel
   within contiguous leaf chunks, Jacobi across chunks via deferred
   cross-chunk scatter — i.e., overlapped couplings use the previous
   iteration's values, exactly Ryan's proposal. Chunks are embarrassingly
   parallel within a sweep → W/p scaling, L ≈ longest single chunk. The
   trade: cross-chunk lag slows convergence → MORE outer iterations.
   Evaluation contract: calibrate tolerance per config (the `accepted` arm
   machinery), compare time-to-accepted, count iterations. Check 021's
   earlier chunked-vs-dagteam A/B history before re-running (why did
   dagteam win? if the answer was "chunked needed more iterations", the
   backoff-fixed executor changes the denominator and chunked may still
   lose; if it was "chunked sweep itself was slow", revisit). Hybrid worth
   noting: dagteam WITH partial pulls inside chunks is compatible.
4. **Sweep-boundary overlap (small, ~0.2 s/solve ceiling):** the serial
   `dagteam_reduce_u!` + laggard wait; per-target incremental reduction
   could release next-sweep roots early. Only worth it after #1.

Discriminating cheap experiment for #1 sizing (local, ≤4 threads or one
short HPC run, Ryan-gated): compute the DAG's weighted critical path L and
width profile directly from an R4 plan (pure graph analysis on
`preds`/`ptot` — no solves needed; `prio[i]` already holds downstream path
weights) for (a) node-level tasks, (b) edge-level tasks, (c) edge-level
with size cutoff. Predicted speedup = L_node / L_edge. Do this FIRST — it
bounds the prize before any executor work.

## Live loose ends

- **Resume job 13878514** (2 h wall, submitted 2026-09-23 ~22:00 MDT,
  `RESUME_FROM_JOB_ID=13875511`): owes the 4 conditional j=32 spin/backoff
  pairs (17 primary stages already ok and analyzed). Babysit via
  `hpc-monitor`; on COMPLETED re-harvest run dir
  `data/p021-cold-20260910/fgs-stage2-13875511` (rsync pattern in
  `fgs_scalability_stage2_results_20260924.md` session; local copy at the
  session scratchpad is gone after reset — re-rsync), re-run
  `julia --project=benchmark -t 1 benchmark/fgs_stage2_analysis.jl <dir>`,
  append the j32 A/B table to the results note. Interest: does backoff at
  32 close the 16→32 plateau (pure DAG-width starvation?) — feeds the
  L-shortening slate.
- **Ryan-gated:** `dagteam_idle=:backoff` as production default
  (recommended); Stage 3 cancellation (recommended — mechanism resolved);
  notebook entry for Stages 1+2 (offer, don't write).
- **Harness batching fix** (Ryan 2026-09-23, see memory
  `feedback-cold-means-cold-start`): "cold" means cold-START solves only;
  batch arms per (j, placement) process in any FUTURE benchmark run
  (runtime-toggleable `dagteam_workers`/`dagteam_idle` + multi-arm driver).
  Build it with (not before) the next benchmark campaign.
- Standing owed: origin push after `gh auth login -h github.com` (branches
  + `campaign/p021-fgs-stage2-20260923` tags, all three repos); 032 ledger
  mirror offer; archiver T5 ruling; stage-1/2 run dirs archive-eligible
  once quiet.

## Traps

- Iterate-changing options (#2, #3) need per-config tolerance calibration
  and time-to-accepted scoring; iterate-preserving (#1, #4) can reuse the
  fixed-work contract directly.
- Fixed-work rows: `solved=false`/`eligible=false` BY CONSTRUCTION; gate on
  certified-accepted + iterations==27 + 1e-8 repeat. Never pool
  instrumented/uninstrumented (`diag_*` = −1 sentinels).
- Aggregate shares at j64-w0 are QUALITATIVE (overhead gate failed −8.5%:
  instrumentation throttles the spin hammer — itself evidence). Backoff and
  cap numbers are uninstrumented and solid.
- dagteam determinism contract: fixed per-row arithmetic order + fixed
  reduction order; any partial-pull design must reduce in fixed
  ascending-source order. Verify bit-identity vs :lexicographic in a
  4-thread smoke (pattern: prior session's `stage2_smoke.jl`, 46 tests).
- `ssh orc` needs a live ControlMaster socket + `bash -lc`; judge runs by
  outputs never sacct; never edit source while a job uses its deployment;
  never touch `/home/rander39/campaigns/*` (ARCHIVER_SKIP).
- `test/runtests_benchmark_cold.jl` fails on the laptop (BLAS pin,
  pre-existing environmental).

## House rules (binding)

Monitoring `hpc-monitor`; harvesting `harvester`; storage `hpc-storage`;
notebook Ryan-gated; dated status/provenance files in BRAINSTORM/021. HPC
submission approval currently covers ONLY resume job 13878514; any
L-shortening benchmark run is a new campaign: worktree/tag pins, provenance
file, and Ryan's explicit go.
