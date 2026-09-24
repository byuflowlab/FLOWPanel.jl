# 021 L-shortening gate-0 RESULTS: R4 dagteam critical-path analysis (2026-09-24)

Discriminating experiment from `fgs_lshortening_reset_prompt_20260924.md`
(slate item #1 sizing): pure graph analysis of the R4 champion dagteam plan —
NO solves, NO timing. Bounds the prize of **edge-level partial pulls** before
any executor work.

- Script: `benchmark/fgs_dag_L_profile.jl` (fixture via `phase1_case.jl` with
  `SKIP_B=1`, geometry-only; champion knobs from
  `retained_r4_champion.toml`: P8/MAC0.4/leaf100/f32full) + offline greedy
  iteration `fgs_lshortening_gate0_20260924/dagL_offline.jl`.
- Run: local laptop, `-t 4` (house rule), FLOWPanel `b86d52c`-dirty,
  FastMultipole `053c8de7`-dirty (branch `flowpanel-20260817`), Julia 1.12.5.
  Graph analysis is machine-independent (plan structure + byte weights).
- Evidence: CSVs + log in `fgs_lshortening_gate0_20260924/`.

## Cost model

Byte-weighted (the campaign's currency; sweep is bandwidth/critical-path
bound), szTM = szTS = 4 (f32full): edge GEMV $e_{ij}=4\,n_i n_j$, node GEMV
$g_i=4\,n_i\,\mathrm{ptot}_i$, leaf LU $4\,n_i^2$, reduction
$4\,n_i|\mathrm{preds}_i|$. Infinite-processor schedules (pure $L$):
node = today's aggregated pull; edge = partials start the moment each pred
publishes, fixed-order reduce + LU after all arrive.

## DAG structure (cross-checks the independent v15 census exactly)

1,068 leaves; 48,167 lower edges (avg 45.1 preds/leaf); unit depth 279;
level width mean/max = 3.8/6; nof min/med/max = 1/54/1450.
W_lower+LU = 809.2 MB/sweep (f32; = v15's 1,513 MB f64 lower bytes / 2 + LU ✓);
upper (filler) 674.8 MB.

## Critical paths (MB streamed per sweep)

| variant | L (MB) | W/L | speedup bound L_node/L |
|---|---|---|---|
| node (today) | 289.5 | 2.8 | 1.00 |
| **edge (θ=0)** | **68.1** | **11.9** | **4.25** |
| edge θ=4KB | 69.2 | 11.7 | 4.19 |
| edge θ=16KB | 105.3 | 7.7 | 2.75 |
| edge θ=64KB | 165.3 | 4.9 | 1.75 |
| edge θ=256KB | 228.0 | 3.5 | 1.27 |
| LU-chain floor | 39.5 | 20.5 | 7.32 |

Node critical path: 261 of 1,068 leaves; composition 265.6 MB GEMV vs
23.9 MB LU — **92% of the chain is pull-GEMV bytes**, exactly what edge
splitting attacks. Byte-weighted W/L = 2.8 (equal-cost census said 3.83;
byte weighting makes the wall slightly worse). Plan-C bound is L-dominated
at every j ≥ 16 in both variants, so **the 4.25× is flat in j** — it neither
needs j64 nor is capped by the plateau.

## The prize does NOT come cheap in task count

| policy | L_node/L | tasks/sweep (today: 1,068) |
|---|---|---|
| split top-1% prio leaves | 1.09 | 1,100 |
| top-10% | 1.34 | 3,259 |
| top-20% | 1.44 | 6,975 |
| top-50% | 2.40 | 20,389 |
| greedy path-driven, converged (815 leaves) | 4.25 | 37,245 |
| θ=4KB size cutoff | 4.19 | 31,605 |
| full split (θ=0) | 4.25 | 49,235 |

- Prio-ranked selective splitting is weak (top-20% → only 1.44×).
- Greedy splitting of successive critical paths converges at **815/1,068
  leaves split** — the critical path is broad, not a thin chain; there is no
  small "hot set".
- Best practical operating point ≈ **θ=4KB: 4.19× at ~31.6k tasks/sweep**
  (~30× today's task count; mean edge task ≈ 16 KB ≈ 1–2 µs of streaming).

## Design implication (feeds the executor spec)

The bound says edge-level partial pulls are worth ~4× on the sweep — but
only if the executor can schedule ~30–50k µs-scale tasks per sweep (×81
sweeps/solve ≈ 2.6–4M tasks/solve) without the scheduler becoming the new
wall. A single SpinLock ready-queue (today's design, 2–6% lockmgmt at 1,068
tasks) will not survive 30× the traffic. This **promotes the parked
static-list-schedule idea from "backoff-equivalent + NUMA" to the natural
edge-task executor**: the DAG repeats 81×/solve, edge costs are exact
(bytes), so per-worker static task lists + per-leaf atomic arrival counters
(no shared queue at all) fit this task granularity, and give NUMA
owner-local first-touch of the split Lmat blocks for free. Alternatives:
sharded/work-stealing deques.

Solve-level expectation at R4 j64 backoff (3.24 s): sweep floor ~2.0 s →
~0.5 s ⇒ solve ~1.7–1.8 s, i.e. **~1.8–1.9× end-to-end**; after that the
FMM (~1.0 s) and the 39.5 MB LU-chain floor set the next walls.

Determinism carries: partials reduced in fixed ascending-source order —
mathematically identical iterate, same certification (no recalibration).

## Caveats

- Infinite-processor L; real executor adds per-task overhead not modeled
  (the reason task count is reported next to L).
- Bytes-as-cost ignores flop-efficiency loss of small GEMVs and cache reuse
  of shared source vectors x_j across a leaf's many out-edges (favorable,
  unmodeled).
- Reduction modeled as one pass over partials; a tree reduction would differ
  slightly. Upper/backward products and the serial boundary reduction are
  outside this L (slate item #4, ~0.2 s ceiling, unchanged).

## Related results banked this session

- Resume 13878514 j=32 pairs harvested (addendum in
  `fgs_scalability_stage2_results_20260924.md`): backoff ladder
  4.155/3.478/3.262 s at j16/32/64 — the 16→32 plateau does NOT close under
  backoff ⇒ pure DAG-width starvation confirmed; Stage 2 fully closed.
- Chunked history check (slate #3): v22 chunked lost on ITERATIONS (44 vs
  27; per-sweep nearfield actually improved) — 14.26 s vs lex 10.95 s.
  Against dagteam+backoff 3.24 s the convergence-lag penalty is even less
  competitive; stays deprioritized.

## Next (Ryan-gated where marked)

1. Draft the edge-partial-pull executor design (static list schedule or
   sharded queues; θ≈4KB aggregation of tiny edges; fixed-order reduction) —
   design doc first, no benchmark.
2. 4-thread local determinism smoke vs :lexicographic once implemented
   (pattern: stage2_smoke.jl, 46 tests).
3. Any HPC benchmark of the new executor = new campaign (worktree/tag pins,
   provenance, Ryan's go). **[Ryan gate]**
4. Standing: backoff-as-default adoption, Stage 3 cancellation, notebook
   entry for Stages 1+2 + this gate-0. **[Ryan gate]**
