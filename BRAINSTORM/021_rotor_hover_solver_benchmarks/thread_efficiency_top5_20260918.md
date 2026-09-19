# FGS thread-efficiency: top 5 ideas, ranked (2026-09-18, rev. b)

Synthesis of `thread_efficiency_proposals_20260918.md` (A1–A9, B1–B5),
`fgs_split_dual_layout_ideas_20260918.md` (split source/target-major cache),
and — in this revision — `fgs_acceleration_recommendation_20260918.md` (the
independent code review). Analysis only. Baseline for all impact claims:
colored@j16 = **10.116 s** accepted (26 iterations; lex ran 27). Ceilings:
29.4 GB/s one core, ~164 GB/s socket, ~330 GB/s two sockets; the 81-sweep
solve streams 231.9 GB of Float64 coefficients (115.9 GB in Float32).

**Revision note (rev. b):** after reading the acceleration recommendation I
swapped the order of the two schedule candidates (source-major row
parallelism now precedes the pull-DAG split), corrected the threaded-residual
claim (shared scratch in `residual!` — a threaded leaf loop races without
worker-private scratch, and the 0.06 s remainder is not all residual time),
replaced my ~3–3.5 s projection with that note's explicit traffic budget
(planning range **4.5–5.9 s** for ideas 1+2), and adopted its correctness
caveats (mixed-precision container change, versioned backward accumulators,
warm-start initialization, colored scatter ordering, handoff budget). Its
quantified handoff budget — 1,068 leaf handoffs/sweep × 81 sweeps = 86,508;
10 µs average costs 0.865 s; staying under 0.5 s added coordination requires
< ~5.8 µs — is the single most decision-relevant number added here and gates
ideas 2 and 3 alike.

## Assessment of the split dual layout note (retained from rev. a, amended)

The core design is sound: the aggregated target-major pull preserves the
lexicographic block-GS recurrence

$$
D_i x_i^{s+1}=b_i^{k}-\sum_{j<i}A_{ij}x_j^{s+1}-\sum_{j>i}A_{ij}x_j^{s},
$$

the triangular split (forward $j<i$ pulled at 52.86% of bytes, backward
$j>i$ pushed as deltas at 47.14%) removes pure pull's stale-RHS penalty
exactly, and every directed interaction is read once per sweep. Two
weaknesses temper it: (a) its published schedule serializes the forward half
into ~1,068 leaf-by-leaf full-team bursts per sweep — the same handoff budget
problem as any per-leaf scheme, and pipelining pulls (below) mitigates but
does not remove the lexicographic critical path; (b) a typical target block
exposes only ~39 output rows, so *intra-block* pull parallelism needs
dot-product reductions, whereas a typical source block exposes thousands of
independently writable rows. The pull design's real parallelism argument is
therefore *inter-leaf* (many pulls in flight), which is unproven until the
DAG's weighted longest path and ready-width are measured. Both points move
the split design behind the source-major candidate, not off the list.

## Top 5, in order

### 1. Float32 coefficient storage with Float64 arithmetic + consumer-aligned first touch

Still first: it halves the dominant stream (231.9 → 115.9 GB per solve) and
multiplies every schedule below. Implementation realities from the code
review, adopted: `Matrices{TF}` uses one type for coefficients and
products/RHS, and `FastGaussSeidel` couples self and nonself precision — this
is an explicit container/dispatch change (separate coefficient type from
accumulator type), **not** a cast plus an assumed fast mixed BLAS call.
Convert coefficient tiles on load and accumulate in Float64; do **not**
round strengths to Float32 to unlock `sgemv` — that is a different numerical
experiment. Keep strengths, products, RHS, residual, self matrices, and leaf
LU in Float64 initially. If specific blocks prove precision-sensitive, retain
them in Float64 (traffic reduction $1-f/2$ for converted fraction $f$)
rather than silently relaxing the gate; the independent BC-accuracy evaluator
remains authoritative because the internal residual of the rounded operator
can be small while the true BC residual fails.

Pair with consumer-aligned first touch: pages placed by whoever will read
them. Placement must match the *selected schedule* — a source-affine fill is
wrong for blocks consumed by several row workers (idea 2); page-aware row
tiling is the eventual owner mapping, with interleave as the experimental
control. Verify placement (numastat) rather than assuming it.

Riders (small, kept honest): threaded `residual!` **after fixing the shared
scratch** — every leaf currently reuses `view(residual_vector, 1:length(rhs))`,
so a naive threaded loop races; use worker-private scratch, and treat the
saving as unattributed until profiled (the 0.06 s/it remainder is not all
residual). Descending-cost/tail-serial color scheduling stays as a
colored-mode control experiment only.

**Impact: required bandwidth to halve wall time drops from 112.7 GB/s
(Float64) to 56.3 GB/s (74.4 with H=0.5 s of coordination) — this is why
mixed storage leads. Effort moderate; risk low-moderate (certification);
unconditional.**

### 2. Persistent, adaptive row-parallel execution of the existing source-major cache

Promoted above the split layout (was #4), adopting the code review's
first-candidate design: keep the cache, the lexicographic leaf sequence, the
far-field refresh, and the RHS semantics. Per leaf: coordinator solves the
diagonal block and publishes strengths; a persistent worker team computes
disjoint *row tiles* of that source's tall matrix (contiguous rows within
each column — column-major streaming, never one strided dot per output row).
Each tile accumulates in private scratch, retaining its old product locally
(which can remove the old-product copy on this path), then applies `+= old`,
`-= new` to its owned target rows — disjoint within one active source, so no
races and no coloring. Small blocks run serial, medium on a compact team,
large across more memory controllers; no allocations or task creation inside
the leaf loop.

Why it now leads: it retains large aggregated GEMVs, contiguous source
strengths, one matrix stream per sweep, and a current RHS for the existing
residual; it exposes thousands of parallel rows on typical blocks (vs ~39 on
target-major); and it is the cheapest implementation of the bandwidth story.
Its go/no-go risk is the **handoff budget**: 86,508 leaf handoffs per solve
must average < ~5.8 µs to keep added coordination under 0.5 s. Gate with a
real-shape sequence benchmark (actual block-size distribution and source
order, including row-to-target maps, old/new bookkeeping, handshakes, and
the small-block policy) before any solver rewrite — multithreaded BLAS on
every small product is not an adequate implementation. Row-tiled custom
kernels are mathematically equivalent, not automatically bit-identical
(tiling/SIMD/FMA change rounding): check, else certify accepted accuracy.

**Impact with idea 1, using the review's budget
$T_{32}\approx 3.0 + 115.9/B + H$: planning range 4.5–5.9 s (1.7–2.3×) at
measured useful B = 40–80 GB/s; meeting the strong 5.058 s target needs
~56–74 GB/s depending on H. Effort moderate; risk = measurable handoff/kernel
efficiency, resolved cheaply by the gate.**

### 3. Pull-DAG split dual layout — the promotion path

Demoted from #2 but kept, in the improved (pipelined) form: triangular split
with the forward $j<i$ half executed as dependency-counter tasks —
`pull(i)` ready when leaf $i$'s earlier neighbors are solved, `solve(i)`
after `pull(i)`; pulls read stable strengths and write only their own rows,
so no write-ordering machinery exists anywhere (the property that made the
exact-GS DAG A9 expensive). Backward $j>i$ delta products run as
work-stealing filler and are owner-applied at the sweep boundary. The two
halves are nearly equal by bytes and touch disjoint matrix data, giving the
only coherence-free two-socket decomposition on the table (forward on socket
0, backward on socket 1, meeting through MB-scale accumulators).

Correctness details adopted from the code review: backward products **must
not overwrite the frozen accumulator $u^s=Ux^s$ while any target still needs
it** — use pending buffers or versions; initialize $Ux^0$ for nonzero warm
starts; keep the fixed external RHS separate from the per-iteration
farfield; split direct-list target-branch row ranges at solve-leaf
boundaries before classifying triangles; do not silently change the
reverse-flag behavior. Honest limits, also adopted: the backward stream
leaves the dependency path but still competes for bandwidth (overlap pays
only when controllers would idle — the two-socket form is what makes it
real); per-leaf publication remains ("one sweep boundary" describes the
deferred backward phase, not all synchronization); and inter-leaf pipelining
is hypothesis until the weighted DAG longest path and ready-width are
measured (supporting evidence: chunked sustained 38.8 active threads on the
same conflict structure).

**Promote if and only if:** idea 2's handoff budget fails or its sequence
benchmark stalls below target, AND the same benchmark shows the lower-pull +
batched-upper schedule winning at equal precision and placement, AND the DAG
width measurement supports pipelining. Prototype the Float64 algebra on
small fixtures first (never conflate layout and precision changes in one
experiment); prefer split triangular storage in production — full duplication
has no sustained byte advantage.
**Impact: unlocks the two-socket ceiling and removes shared-RHS coordination;
effort high; risk moderate, front-loaded into gates shared with idea 2.**

### 4. Joint staircase retune at the new operating point (byte-shedding included)

Unchanged in substance: leaf=100 / MAC=0.4 / P=8 / inner=3 were selected for
the serial Float64 cost balance; retune all four axes jointly after the
schedule/precision work lands, testing in particular whether higher P buys a
looser MAC that moves marginal direct blocks out of the streamed cache
(bytes removed beat bytes accelerated — the one lever here not bounded by
the 164 GB/s arithmetic). Two amendments from the review: report it
separately from exact-iterate improvements (it changes the approximate
operator/partition/schedule), and a looser MAC needs independent accuracy
checks — higher P may compensate but is not guaranteed to. The review's
remainder analysis motivates this strongly: after mixed storage, a ~3 s
non-product remainder caps bandwidth-only gains (free products at 80 GB/s
would still leave ~4.45 s → only 1.48× more), so beyond ~2× the wins must
come from less work — this idea and idea 5.
**Impact tens of percent historically; effort one campaign; risk nil;
unconditional.**

### 5. Safeguarded Anderson acceleration of the outer FGS map

Kept fifth, enriched by the review (which independently ranks it first among
algorithm-changing options, above FGMRES): small history (3–8), restart on
ill-conditioned history, accept extrapolated iterates only under a residual
safeguard. The integration cost it flagged is the one that matters: after
mixing, the incremental nonself products and farfield state are stale —
**mixing strengths while retaining stale incremental products is wrong** —
so each acceptance may add an operator pass; compare *total streams and FMM
calls*, not iteration counts. 27 → 12–18 is a hypothesis to be tested, not a
forecast. FGMRES with the existing `FGSPreconditioner`
(`src/FLOWPanel_solver.jl:1873–1972`) is the fallback within this slot: use
the existing path for correctness, but note its apply performs a complete
fixed-iteration FGS solve with save/restore — a cheaper inner apply must be
justified, and a lower Krylov count can conceal more total FGS work.
Un-park only when the post-idea-1/2 profile shows iteration work, not
bandwidth, as the binding term. If a hypothetical 27→15 held at the
optimized per-update cost, 4.5–5.9 s → roughly 2.5–3.3 s *before* added
acceleration costs — the route to ~3–4×, conditional on convergence
evidence.
**Effort moderate-high; risk moderate; conditional on the profile.**

## Cut from the five (and why)

- **Both sockets as a standalone item** — folded into idea 3, its only
  coherence-free form; standalone it is more than a launcher change
  (ownership/pinning must already be correct) and source-level handshakes
  plus remote data can erase the gain.
- **Exact-GS dependency DAG with write ordering (A9)** — superseded by idea
  3's pull form; and per the review, removing barriers does not remove the
  lexicographic critical path — measure before believing.
- **Parallel-by-target scatter (A3) as its own item** — absorbed: idea 2's
  row-tile ownership does the same job within one active source. If applied
  to colored mode, the order to preserve is color-major then ascending
  source *within* the color, and the `+= old` / `-= new` operations must stay
  separate for bit identity.
- **Chunked rescue, THP, FMM-overlap (extra stale farfield), true-reverse
  SSOR, low-rank compression** — each small, high-risk, or dominated: the
  measured FMM cost is too small for overlap to meet any target; compression
  saves only if $r(m+n) \ll mn$ against a ~39-wide small dimension; MAC/P
  retuning (idea 4) is the cheaper byte-shedder. GPU is a separate platform
  project, not a near-term recommendation from these CPU results.

## Recommendation

Run ideas 1 and 2 as one program with the review's gate sequence: (i)
correctness model on small fixtures in Float64 first (nonzero starts,
multiple systems, branch-spanning targets, rigid transforms); (ii) the
real-shape sequence benchmark — serial BLAS vs persistent source rows vs
split lower/upper kernels at equal precision and placement, measuring useful
bandwidth, handoff latency, and verified page placement — which
simultaneously prices idea 2's handoff budget and idea 3's promotion
condition; (iii) independent Float32 certification with recalibrated
tolerance; (iv) end-to-end interleaved A/B against the unchanged champion,
requiring ≥1.5× accepted throughput (≤6.744 s) and designing for ≥2×
(≤5.058 s). Planning range for 1+2: **4.5–5.9 s**. Then use the new profile
to pick among idea 3 (if handoffs bound), idea 4 (if bytes bound), or idea 5
(if iterations bound) — beyond ~2× the arithmetic says the wins come from
less work, not more scheduling machinery. All implementation, submissions,
and campaign ceremony remain Ryan-gated; local checks ≤4 threads.
