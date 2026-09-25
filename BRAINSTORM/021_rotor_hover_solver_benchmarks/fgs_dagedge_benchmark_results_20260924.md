# 021 :dagedge benchmark results — dagedge LOSES (+6.4% @ j64); busy-time inflation, not scheduler idle (2026-09-24)

Campaign: `fgs_dagedge_benchmark_provenance_20260924.md` (binding caveats there).
Jobs: perf **13879622** (m12-4-29, 2026-09-23T23:44 → 2026-09-24T01:27) and
profile **13879625** (m12-3-30, 03:28 → 04:35). Both `COMPLETED` with
`failed_count=0`, all STATUS_* ok, 75 perf + 40 profile solves, zero
failures/outliers, all rows `accepted=true`. Harvested copies:
`data/p021-cold-20260910/dagedge_harvest_20260924/{perf,profile}/`.

## Verdict

**:dagedge does not beat the Stage-2 champion dagteam+backoff — it is
slower at every j.** Gate-0's predicted 4.25× sweep gain (expected ~1.7–1.8
s end-to-end @ j64) did not materialize; the decision threshold ("< ~3×
sweep gain ⇒ overhead ate the prize") trips decisively. The profile run
attributes the loss to **per-task busy-work inflation (4.15× more total
busy time)**, *not* scheduler idle or static-schedule skew — dagedge
actually *improves* both idle and imbalance. The prize was eaten by
per-edge-task fixed cost: 29.5× more task executions at ~8.2 µs each vs
the ~0.3 µs/task the R2 extrapolation assumed (27× miss in the gate-0 cost
model).

Recommendation stack for Ryan (no submissions pre-approved beyond this
campaign): **retain dagteam+backoff as the candidate production default**;
dagedge at θ=4096 is retired as a loss. The recorded fallbacks
(work-stealing deques, NUMA first-touch repack) target idle/imbalance,
which this data shows is *not* the dagedge bottleneck — they would not
rescue it at this granularity. The only dagedge-shaped survivor is the
θ→large limit, which simply converges back to dagteam (see θ probe).

## Q1 — Paired ranking @ j64 (perf run, primary blocks; ranking-grade)

Same-process pairs (same placement/compile state), medians of 5 solves per
executor per block, dagedge θ=4096:

| block | dagteam (s) | dagedge (s) | Δ (s) | dagedge/dagteam |
|---|---|---|---|---|
| primary-p1 | 3.533 | 3.729 | +0.196 | 1.055 |
| primary-p2 (order swapped) | 3.519 | 3.786 | +0.267 | 1.076 |
| primary-p3 | 3.604 | 3.766 | +0.162 | 1.045 |
| **pooled median** | **3.540** | **3.766** | **+0.226** | **1.064** |

Stage-2 champion reference 3.24 s reproduces here as 3.54 s (different
node/day; paired deltas are the ranking-grade signal). Expected ~1.7–1.8 s
did not materialize.

## Q2 — Ladder in j (perf run)

| j | dagteam (s) | dagedge (s) | dagedge penalty |
|---|---|---|---|
| 16 | 4.539 | 5.162 | +13.7% |
| 32 | 3.788 | 3.945 | +4.1% |
| 64 | 3.540 | 3.766 | +6.4% |

dagteam baseline matches the Stage-2 backoff ladder (4.155/3.478/3.262 s)
to within node/day scatter. dagedge is "flat in j" only in the sense that
it tracks dagteam from above — it never crosses.

## Q3 — θ probe (perf run, theta-j64; one shared process — do NOT rank from these)

| θ (B) | plan_ntasks | sim_edge_L (MB/sweep) | dagedge median (s) |
|---|---|---|---|
| 0 | 50300 | 68.1 | 3.829 |
| 4096 | 33640 | 69.1 | 3.767 |
| 16384 | 11246 | 103.5 | 3.540 |

Monotone in task count, not in edge traffic: fewer/bigger tasks win even
at +50% simulated edge bytes. Fixed per-task cost, not memory traffic, is
the controlling term — at θ=16 KB dagedge merely converges to dagteam's
time. Consistent with the attribution in Q4.

## Q4 — Profile attribution (profile run ONLY; never pool with perf — instrumentation shifts timing)

Medians per arm; durations are summed worker-time over one fixed-work
solve (27×3 iterations) and **overlap — never sum as elapsed**. wait/reduce
are small everywhere; `lockmgmt` is **identically 0 for all dagedge rows**
(lock-free path verified).

| arm | executor | solve (s) | busy_lower (ms) | idle (ms) | wait (ms) | reduce (ms) | busy_max/min |
|---|---|---|---|---|---|---|---|
| prof-j16 | dagteam | 4.206 | 4 804 | 17 935 | 0 | 191 | 2.40 |
| prof-j16 | dagedge | 4.722 | 11 118 | 15 745 | 20 | 217 | 1.26 |
| prof-j32 | dagteam | 3.565 | 4 933 | 43 918 | 0 | 187 | 4.17 |
| prof-j32 | dagedge | 3.562 | 13 841 | 29 378 | 30 | 214 | 1.57 |
| prof-j64-b1 | dagteam | 3.244 | 5 054 | 97 292 | 0 | 188 | 8.47 |
| prof-j64-b1 | dagedge | 3.311 | 20 984 | 77 277 | 47 | 207 | 1.97 |
| prof-j64-b2 | dagteam | 3.281 | 5 081 | 97 987 | 0 | 197 | 8.78 |
| prof-j64-b2 | dagedge | 3.382 | 21 860 | 77 871 | 48 | 199 | 1.92 |

Task-level accounting @ j64 (b1):

| | dagteam | dagedge | ratio |
|---|---|---|---|
| plan_ntasks | 2 133 | 33 640 | 15.8× |
| n_lower (task executions/solve) | 86 508 | 2 552 310 | 29.5× |
| total busy_lower | 5.05 s | 20.98 s | 4.15× |
| per-task busy | 58.4 µs | 8.22 µs | 0.141× |

**Reading:** dagedge did exactly what it was designed to do on the
scheduler axis — idle @ j64 drops 97.3 → 77.3 s and busy imbalance
collapses 8.5 → 2.0 — and lost anyway, because each edge task carries
~8.2 µs of busy cost vs the ~0.3 µs/task the R2 extrapolation predicted
(≈27× over model). 29.5× more executions × 0.141× per-task cost = 4.15×
total busy inflation, which swamps the ~20 s of idle it recovered. The
gate-0 4.25× sweep bound assumed the per-task cost would stay near the
measured R2 point; at 4 KB granularity a fixed per-task overhead
(task launch/bookkeeping + partial-pull redundant loads) dominates.
Static-schedule skew — the recorded risk — is refuted as the failure mode.

## Q5 — Sanity

- `sim_edge_L` = 69.1 MB/sweep at θ=4096 @ R4 — matches the gate-0
  expectation of ≈68–69 MB (build-time simulation output, not a
  measurement).
- plan_ntasks stable at 33 640 (θ=4096) across all phases/j — schedule is
  deterministic as designed.
- Cross-executor delta: max **2.24e-7**, far under the 1e-5 tripwire
  (R2 measured 3.4e-7) — dagedge and dagteam compute the same answer.
- reduce time ~190–220 ms everywhere, independent of j, θ, executor.

## Caveats (binding, from provenance)

- Fixed-work gates (27 iterations, certified evaluator, 1e-8 repeat) held;
  all rows accepted.
- Profile solve times carry instrumentation overhead (Stage 2 measured
  8.5%); ranking numbers come from the perf primary pairs only.
- θ-probe arms share one process (construction-order effects possible) —
  used only for the monotonicity read, not ranking.
- diag_* −1 sentinels handled (dagteam wait/lock columns clean).

## Standing Ryan gates (surfaced, not acted on)

1. **backoff as :dagteam production default** — recommendation now
   *strengthened*: dagedge, the last live L-shortening challenger at this
   granularity, lost. Stage-3 cancellation recommendation stands.
2. No new dagedge follow-up submission recommended; if any L-shortening
   thread continues, the data points to reducing per-task fixed cost
   (coarser θ ladder or fusing edge tasks), not scheduling — Ryan's call.
3. Notebook entry (Stages 1+2 + gate-0 + dagedge prototype + this
   campaign): ready to offer.
4. Origin push of branches + tags (`campaign/p021-fgs-stage2-20260923`,
   `campaign/p021-fgs-dagedge-20260924`, all three repos) still pending
   `gh auth login`.
5. Remote run dirs (`fgs-dagedge-{perf,profile}-*` under
   `p021-cold-20260910/`) are small (no VTK) — no archiving urgency.

## Thread-scalability thread PARKED (Ryan 2026-09-24)

Headroom remains — 77 s of summed worker idle at j64 even under dagedge, with
identified levers (NUMA first-touch + small-block repack first, then a
coarser-θ ladder, then per-worker edge fusing) — but it is NOT pursued for
now. Cold-solve verdict: at R4 j≥8, krylov_ilu_nfcache (2.41 s @ j64) beats
the best FGS (dagteam+backoff, 3.24 s @ j64) at every measured j and is the
only config still scaling at 64 threads. dagteam+backoff was adopted as the
FGSSolver default the same day (src/FLOWPanel_solver.jl; FastMultipole
defaults untouched, dagteam_precision stays :f64). The decisive question
moves to the warm-started Phase-3 R4 slice
(`fgs_warmstart_r4_reset_prompt_20260924.md`).
