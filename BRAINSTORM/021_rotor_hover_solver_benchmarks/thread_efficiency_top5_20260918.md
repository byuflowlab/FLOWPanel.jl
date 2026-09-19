# FGS thread-efficiency: top 5 ideas, ranked (2026-09-18)

Synthesis of `thread_efficiency_proposals_20260918.md` (A1–A9, B1–B5) and
`fgs_split_dual_layout_ideas_20260918.md` (split source/target-major cache +
its ranked five). Analysis only. All impact figures are measured against the
colored@j16 = 10.116 s operating point and bounded by the measured ceilings:
29.4 GB/s one core, ~164 GB/s socket, ~330 GB/s two sockets; 8.6 GB streamed
per iteration (3 sweeps × 2.86 GB) at 27 outer iterations.

## Assessment of the split dual layout note

The core design is sound and I consider it the strongest large idea on the
table. Its three load-bearing observations hold up against the code and the
R4 census:

1. A target-major *aggregated* pull (one block per solve leaf, not per direct
   pair) keeps ~1,068 large GEMVs per sweep and reads the lexicographic GS
   strength mixture directly — same mathematics, no chunked-style Jacobi lag,
   no expected iteration inflation.
2. The pure-pull disadvantage (stale RHS → an extra 2.86 GB residual stream,
   or a 15.7% averaged reread penalty) is real, and the triangular split
   ($j<i$ pulled, $j>i$ pushed as deltas) eliminates it exactly: every
   directed interaction read once per sweep, full RHS current at sweep end.
3. The forward half (52.86% of bytes) becomes race-free target-owned reads;
   the backward half (47.14%) is off the current sweep's dependency path.

One weakness, and it matters at j64: **the note's sweep schedule (pull leaf
i → solve i → enqueue backward → next leaf) serializes the forward half into
~1,068 leaf-by-leaf full-team bursts per sweep** — 6.8× more sync points than
colored's 79 colors × 2 phases, whose barrier cost is what made colored lose
at j64. The fix is already latent in the design and is the improvement I
propose below (idea 2): because forward pulls have *no shared writes*, the
write-ordering problem that made barrier-free lexicographic scheduling hard
(proposal A9's per-row ascending-source machinery) disappears entirely, and
the forward sweep pipelines with plain per-leaf dependency counters.

The note's own top-five ranking is close to mine. I differ in three places:
its #2 (persistent adaptive threaded GEMV on the unchanged source-major
layout) should be explicitly subordinated to the split design as the
gate-failure fallback, not a peer; its #3 (both sockets) is strongest *as a
mode of the split layout* (the two halves give a coherence-free socket
decomposition) and folds into idea 2 rather than standing alone; and the
split layout itself belongs in the five, not only in the closing paragraph.

## Top 5, in order

### 1. Float32 cache + consumer-aligned first touch (+ the free exact riders)

Unchanged from both prior docs; still first. Store the nonself operator in
Float32 (strengths, accumulators, self/LU stay Float64 initially); build each
page with the thread that will consume it (measured: chunk-affine first touch
164 GB/s vs 74 serial-touch, and it beats interleave's 154). Halves the
dominant stream and multiplies every scheduling idea below; a complete *dual*
Float32 layout also costs no more memory than today's single Float64 cache,
which makes idea 2's prototype cheap. Ride along the bit-identical
micro-fixes from the first proposals doc: threaded max-abs `residual!`,
threaded per-iteration vector ops (A4), and descending-cost/tail-serial color
scheduling (A5) for whatever life colored mode has left. Numerics change →
standard tolerance staircase + accepted-accuracy certification.
**Impact ~2× alone (10.1 → ~5 s); effort moderate; risk low; unconditional.**

### 2. Pull-DAG split dual layout (improved form of the split design) — the big one

Triangular split exactly as in `fgs_split_dual_layout_ideas_20260918.md`,
with a different forward-half schedule:

- **Forward ($j<i$, 52.86%): dependency-counter pipelined pulls, no
  barriers.** Task graph per sweep: `pull(i)` becomes ready when every
  earlier neighbor $j<i$ of leaf $i$ has been solved this sweep; `solve(i)`
  follows `pull(i)`. Pulls only read strengths (stable after the producing
  solve) and write leaf $i$'s own rows, so there is **no write ordering to
  enforce anywhere** — the property that made A9 expensive. Many leaves'
  pulls run concurrently, each on one or a few workers, saturating bandwidth
  without team-wide handoffs. Deterministic at any thread count because each
  target row's sum is formed in one place with a fixed column order.
- **Backward ($j>i$, 47.14%): work-stealing filler.** After `solve(j)`,
  enqueue the $\Delta x_j$ source-major products; idle workers drain the
  queue between pull dependencies. Applied to the per-target accumulators by
  their owners in ascending source order at the single sweep-boundary sync
  (3 syncs per iteration, total).
- **Two-socket extension (absorbs the note's idea 3):** the halves are
  nearly equal by bytes and touch disjoint matrix data — forward team on
  socket 0, backward on socket 1, meeting only at the boundary through the
  small accumulators and the ~MB-scale strength vector. This is an unusually
  coherence-free dual-socket decomposition (~330 GB/s aggregate) and is the
  *only* proposal on the table that uses socket 1 without paying cross-socket
  barrier latency inside the sweep.

Same lexicographic block-GS iterate as today up to floating-point regrouping
(recalibration needed — but idea 1 already forces that). Supersedes colored,
chunked, the by-target scatter (A3), and the DAG-with-ordering (A9).

Risks, honestly: transposed/target-major GEMV bandwidth on zen3 is unproven
(store $T_i$ with the long source dimension contiguous and it should stream,
but measure it); the critical path through the leaf-dependency chain could be
long in adversarial orderings (mitigating evidence: chunked sustained 38.8
active threads on the same conflict structure, so parallel width exists);
direct-list target branches spanning several solve leaves must be split at
assembly; constructor and transform/reuse paths must carry two
representations. **Gate before committing** (the note's own gate, extended):
an R4-shaped microbenchmark of aggregated transposed pulls vs source-major
pushes at equal precision/placement, plus a toy pipelined-pull scheduler. If
pull bandwidth ≥ push bandwidth, build it.
**With idea 1: products ~0.026 s/it at one socket (~0.013 at two), one sync
per sweep → projected wall ~3–3.5 s (~3×); effort high; risk moderate,
front-loaded into a cheap gate.**

### 3. Joint staircase retune at the new operating point (byte-shedding included)

The note's idea 4 = prior doc's A6, kept at full strength and run *after*
ideas 1–2 change the cost balance: leaf=100 / MAC=0.4 / P=8 / inner=3 were
selected when sweeps ran serially at 29.4 GB/s. Retune all four axes jointly;
in particular test whether higher P buys a looser MAC that moves marginal
direct blocks out of the streamed cache entirely (bytes removed beat bytes
accelerated — this is the only Tier-A entry not bounded by the 164 GB/s
arithmetic). No code change; standard calibration protocol.
**Impact unknown a priori, historically tens of percent; effort one staircase
campaign; risk nil; unconditional.**

### 4. Persistent-team adaptive threaded GEMV on the existing source-major layout — the fallback

The note's idea 2, subordinated: if the idea-2 gate shows target-major pulls
cannot match source-major push bandwidth, keep the current cache and
parallelize *within* each source GEMV with a persistent spinning worker team
(atomic per-leaf handoffs, not fork/join), small blocks serial, large blocks
tiled across workers, pages striped across the participating workers' nodes.
Improve it with the bit-identical parallel-by-target scatter (A3) so the
serial scatter doesn't become the new critical section. Preserves the
incremental RHS semantics almost exactly. The same ~1,068 handoffs/sweep
concern applies as to the note's original split schedule — the persistent
team makes each handoff ~µs instead of a fork/join, which is what makes this
viable where colored's barriers were not. The idea-2 microbenchmark answers
this one for free.
**Impact with idea 1: perhaps 2–2.5× total; effort moderate; risk low;
conditional on the gate.**

### 5. Anderson/FGMRES acceleration of the outer fixed point

After ideas 1–2, per-iteration cost approaches its floor (~0.08–0.1 s/it,
increasingly FMM- and sync-dominated) and the 27 outer iterations × 3 sweeps
= 81 cache streams become the dominant multiplier. Anderson mixing or FGMRES
with one FGS sweep as preconditioner plausibly cuts 27 → 12–18, multiplying
*everything* — kernels, placement, FMM calls, syncs — where the second socket
would shave only the products term (~15% of a post-idea-2 iteration).
Algorithm-changing: new convergence/memory/accuracy tuning, and it is
currently parked — un-park it only once the post-idea-2 profile confirms
iteration count is the binding term.
**Impact up to ~1.5–2× on top of 1+2 (toward ~2 s); effort moderate-high;
risk moderate (restart/window tuning).**

## Cut from the five (and why)

- **Both sockets as a standalone item** — folded into idea 2, where the
  forward/backward split gives it a coherence-free form; standalone (on
  colored or lex) it inherits cross-socket barrier costs and is bounded to
  tens of percent after Float32.
- **A9 dependency-DAG exact-GS with write ordering** — superseded by idea 2,
  which gets the same barrier-free schedule without the per-row ordering
  machinery, at the price of a recalibration idea 1 forces anyway.
- **Chunked rescue (fewer chunks + under-relaxation)** — only worth revisiting
  if both idea 2 and idea 4 fail their gates; the 27→44 inflation is a
  structural majority-Jacobi cost the split layout avoids by construction.
- **THP/huge pages, FMM-sweep overlap, true-reverse SSOR, ACA compression** —
  each ≤~5% or high-effort/low-confidence; keep as opportunistic notes in the
  first proposals doc.

## Recommendation

Run ideas 1 and 2 as one program: idea 1's Float32 + affine-touch work is
needed by every branch, and the dual-Float32 prototype it enables is exactly
the vehicle for idea 2's gate. Decision tree: microbenchmark pull-vs-push
bandwidth (cheap, local-scale first, one m12 job if promising) → pull ≥ push
⇒ build the pull-DAG split (idea 2), else ⇒ persistent-team source-major
(idea 4). Retune (idea 3) after the winner lands; un-park Anderson/FGMRES
(idea 5) if the post-landing profile shows iterations, not bandwidth, as the
binding term. Honest composite ceiling: ~3–3.5 s via 1+2, ~2 s if idea 5
also pays — nothing on this list can beat the bandwidth arithmetic without
cutting bytes (1, 3) or iterations (5), and the projections respect that.
