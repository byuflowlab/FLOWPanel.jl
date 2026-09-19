# Thread-efficiency proposals for the FGS solve (2026-09-18)

Executes `thread_efficiency_reset_prompt_20260918.md`. Analysis/proposal only —
no code was changed, nothing submitted. Code read directly in FastMultipole
`flowpanel-20260817` @ `c18e4b46` (`src/solve.jl`, `src/containers.jl:1092-1148`,
`src/nearfield_cache.jl:440-489`) and FLOWPanel (`src/FLOWPanel_solver.jl`,
`benchmark/fgs_cold_common.jl`, `benchmark/retained_r4_diagnostics.toml`).

## Frame: what the code actually does, and where the time is

R4 retained config: P=8, MAC=0.4, leaf=100, **inner_iterations=3,
reverse_pass=false**, tolerance 3.479e-7 (`retained_r4_diagnostics.toml`).
Per outer iteration: one threaded `fmm!` farfield pass → influence mapping →
`residual!` (serial) → **3 GS sweeps**. Each sweep streams the full 2.86 GB
nonself cache once (≈8.6 GB/iter) plus the leaf-LU factors.

Structural facts that matter for threading (verified in source):

1. **`:lexicographic` sweeps are 100% serial** (`gs_sweep!`, solve.jl:1277-1302).
   lex@j64's 10.95 s is a *serial* sweep chain (0.346 s/it × 27 ≈ 9.3 s = 85%
   of wall) plus a threaded FMM. The global best colored@j16 (10.116 s) beats
   it by only 8% because the parallel product phase is capped by single-node
   bandwidth (93% of cache pages on one NUMA node) and pays 79 colors × 2
   phases × 81 sweeps of barriers.
2. **Placement is decided by a serial first touch**: `_calloc_vector`
   (solve.jl:6-11) + the sequential fill loop (solve.jl:256-321). Confirmed by
   job 13763831: 74 GB/s delivered vs 164 GB/s chunk-affine ceiling.
3. **The serial scatter preserves a per-target-row ascending-source
   application order** (`scatter_nonself_influence!`, solve.jl:948-984). This
   is the exact property that makes a *bit-identical parallel-by-target*
   scatter possible (proposal A3).
4. **`residual!` returns max-abs** (solve.jl:1803-1839, returns `mae`) — an
   order-independent reduction, so it can be threaded per leaf bit-identically.
5. **A proven parallel assembly pattern already exists in-repo**
   (`nearfield_cache.jl:455-478`: `@spawn` workers, atomic chunk pool,
   per-worker private buffers because `direct!` probing mutates shared target
   buffers) — directly reusable for proposal A1, which needs exactly that
   privatization (the FGS fill loop mutates `target_buffers` via
   `reset!`/`direct!`).

**The ceiling everything is ranked against** (numa findings 2026-09-18): one
core sustains 29.4 GB/s; the socket delivers at most ~164 GB/s regardless of
scheduling → the streamed-products term can never improve more than
164/29.4 ≈ 5.6× per socket. Ideas must therefore do one of: (i) reach that
ceiling (placement + scheduling), (ii) cut bytes (Float32, cache-size retune),
(iii) cut iterations/sweeps, (iv) add memory controllers (socket 1), or
(v) shrink the serial/barrier remainder that the ceiling doesn't cover
(scatter, residual, color barriers, FMM overlap).

Budget at the operating point (per iteration, j64 lex counters): products
0.292 s + scatter 0.032 s + leaf-LU 0.020 s = 0.346 s; everything else (FMM,
influence mapping, residual, vector ops) ≈ 0.06 s/it.

## Tier A — algorithm-preserving (same iterate sequence; A1/A3/A4/A5/A9 bit-identical, A2 numerics-changing but Ryan-approved)

Ranked by (impact ÷ effort), impact measured against colored@j16 = 10.116 s.

### A1. Parallel, sweep-affine first-touch assembly of `nonself_matrices` (and `self_matrices`/LU) — **endorse the staged "chunk-affine" idea, with one refinement**

Parallelize the fill loops (solve.jl:256-321, :450-505) with the
`nearfield_cache.jl` worker pattern (private buffer copies), partitioning **by
source leaf using the same leaf→thread mapping the sweeps will use** (colored:
`@threads` splits each color's ascending list contiguously and identically
every sweep; chunked: `:static` over `chunk_ranges`). That yields *local*
pages at sweep time, not merely spread pages — arm d beat interleave 164 vs
154 GB/s for exactly this reason. Bonus: the constructor itself (currently a
serial multi-GB fill) parallelizes, which shortens solver (re)builds.

- Impact: ×2.2 on the streamed-products term at j64 (74→164 GB/s measured on
  the exact analogue); products 0.292 → ~0.05 s/it *if* the sweep's parallel
  phase can consume it (needs A3/A5 to shed barrier/serial cost, else gains
  saturate around the colored@j16 structure: still ~1.3-1.5× alone).
- Effort: moderate (pattern exists; determinism unaffected — blocks are
  written disjointly, values identical). Risk: low. Bit-identical.

### A2. Float32 nearfield storage (staged, Ryan-approved) — halves bytes

Store `nonself_matrices` (+ leaf-LU factors) as Float32, accumulate products
in Float64 (convert strengths per sweep; `mul!` mixed or manual kernel).
Multiplies with A1: 8.6 GB/iter → 4.3, products → ~0.026 s/it at 164 GB/s.
Iterate changes slightly → tolerance recalibration by the standard staircase;
acceptance gate re-certifies true accuracy. Risk: low-moderate (conditioning
of leaf LU in F32 — keep LU in F64 if the acceptance gate flags it; LU stream
is only ~6% of the chain).

### A3. Parallel-by-target deferred scatter — bit-identical, unlocks j64

Precompute the inverse map target-leaf → list of (source leaf, segment
offset) *in ascending source order* (from `index_map`/`sorted_list`, exactly
the machinery `build_chunk_map` already walks). At each color boundary (or
sweep end), thread over target leaves; each thread applies its row-range's
contributions in ascending source order — the identical per-row floating-point
sequence as the serial loop, hence **bit-identical**, deterministic at any
thread count. Direct win is small (scatter 0.032 s/it), but it removes the
serial section between every pair of color barriers, which is what makes
colored lose at j64 (+6.7% vs lex). Prereq for A1's full value.
Effort: moderate. Risk: low (the ordering argument is checkable in a unit
test against the serial path).

### A4. Thread `residual!` and the per-iteration vector ops — bit-identical

`residual!` streams every self block serially each outer iteration and
returns a max — thread per leaf with a max-reduction (order-independent →
bit-identical). Likewise `strengths_old` delta and `right_hand_side` updates
are serial `broadcast`s over N-body-length vectors. Recovers a slice of the
~0.06 s/it remainder (~0.02-0.03 s/it). Effort: small. Risk: negligible.

### A5. Colored-schedule engineering (barrier-spin reduction) — bit-identical

Three cheap scheduler changes, none touching the iterate:
(a) execute **tail colors serially** — greedy coloring produces many tiny
colors (79 colors, median width 16); a color with < ~2×nthreads·cost leaves
pays more in fork/join+spin than it saves; (b) within a color, schedule
leaves in **descending cost order** (costs known from `sizes`) to cut
straggler spin (writes disjoint → any order valid); (c) keep the leaf→thread
map stable across sweeps and matched to A1's fill affinity. Measured barrier
spin was 10.7 CPU-s per 0.28 s span (chunked); colored@j64's loss is the same
disease. Effort: small. Risk: low.

### A6. Re-tune the staircase (leaf, MAC, P, inner) at the *threaded* operating point — no code change

leaf=100 / MAC=0.4 / inner=3 were selected when sweeps ran serially at
29.4 GB/s. Once A1+A2 make sweeps ~4× cheaper, the optimum moves: the
FMM-vs-nearfield byte split (leaf, MAC) and the sweeps-per-FMM ratio (inner)
should be re-selected on the standard staircase with recalibrated tolerance.
Two directions to test, both byte-cutting or FMM-shifting: smaller leaf /
higher MAC shrinks the streamed cache (moves work into the threads-well,
compute-bound farfield); more inner sweeps per outer iteration amortizes the
FMM if post-fix sweeps become cheap relative to it. This is a parameter
re-tune inside the existing calibration protocol — same algorithm family.
Effort: one staircase campaign. Impact: unknown a priori but historically
the staircase has moved operating points by tens of percent.

### A7. 2 MB pages (THP) for the cache — small, nearly free

2.86 GB streamed through 4 kB pages = ~715k TLB entries/sweep. Check THP mode
on m12 (`/sys/kernel/mm/transparent_hugepage/enabled`); if not `always`, an
`madvise(MADV_HUGEPAGE)` on the calloc'd region (or posix_memalign+madvise in
`_calloc_vector`) typically buys a few percent on pure streams. Note the
interaction with A1: THP granularity (2 MB) is fine-grained enough for
per-leaf placement (blocks ≈ 2.9 MB). Effort: tiny. Impact: ≤5%.

### A8. Use both sockets (cpubind 0-7) — doubles the ceiling, only pays after A1-A5

membind/cpubind 0-7 gives 16 channels ≈ ~330 GB/s ceiling → streaming term up
to ~11× over one core. But cross-socket barriers/scatter latency grow, and at
128 threads the barrier+serial remainder dominates unless A3/A4/A5 are in.
FMM farfield also gains. Sequence this *after* the single-socket fixes; expect
~10-20% additional wall once products are ~0.026 s/it (they'd go to ~0.013,
against a ~0.1 s/it non-product floor). Effort: launcher change only. Risk:
low (placement per A1 must be per-thread-local, which it already is).

### A9. Dependency-DAG (wavefront) exact-GS scheduling — the large algorithm-preserving change

Replace coloring with task-graph execution of the *lexicographic* sweep:
leaf i's solve depends on all scatters into its rows from leaves j<i; scatters
into a row execute in ascending source order (per-row sequential semantics →
**bit-identical to lex**, which is the calibrated 27-iteration reference —
no recalibration, no coloring-induced ordering change). The conflict structure
is already computed (`color_leaves` adjacency); execution via per-leaf atomic
dependency counters + work queue removes *all* global barriers and pipelines
solves, products, and scatters. This strictly dominates colored-mode
semantics; combined with A1/A2 it is the schedule most likely to actually
reach the 164 GB/s ceiling. Effort: high (the one big rewrite in this tier).
Risk: moderate — correctness of the DAG and of per-row ordering enforcement;
mitigated by bit-identity testing against serial lex.

## Tier B — algorithm-changing (authorized second tier)

### B1. Krylov/Anderson acceleration of the outer fixed point (currently parked)

The only lever that multiplies *everything* — 27 outer iterations bound the
streamed bytes, the FMM calls, and the barrier count alike. FGS-preconditioned
FGMRES (already parked) or Anderson mixing on the outer loop plausibly cuts
27 → ~12-18. After A1+A2 push per-iteration cost toward its floor, iteration
count becomes the dominant term; recommend un-parking this *after* the Tier A
bandwidth work lands, since its payoff is multiplicative either way but its
risk (restart tuning, memory) is better spent once per-iteration cost is flat.

### B2. Chunked rescue: fewer chunks + cross-chunk under-relaxation

Chunked's mechanism worked (span 0.277 s/it, 38.8 active threads); it lost on
27→44 iteration inflation from majority-Jacobi coupling. chunks=8-16 (instead
of 64) makes most coupling intra-chunk GS again; optionally under-relax only
the cross-chunk (Jacobi) corrections. If iterations recover to ≤30, chunked +
A1 + A2 projects competitive with colored+A3. Cheap A/B; but A9 supersedes
this direction if attempted — recommend only as fallback if A9 is declined.

### B3. Slim the streamed cache at fixed accuracy (beyond A6's retune)

Low-rank (ACA) compression of the marginal direct-list blocks — pairs that
barely failed MAC=0.4 are nearly-separated and compress well; bytes drop where
rank r ≪ min(m,n). High implementation effort, and A6's cheaper knob (raise
MAC so those pairs move into the FMM) captures most of the same byte cut under
the existing acceptance gate. Rank low; revisit only if A6 shows the
byte-vs-accuracy frontier is binding.

### B4. Overlap the farfield FMM with sweeps via one extra farfield lag

Compute iteration k+1's farfield concurrently with iteration k's sweeps using
one-iteration-stale strengths on dedicated threads. Hides ≤0.06 s/it of FMM
but risks iteration inflation (the same disease that killed chunked) for a
≤6% upside. Rank last.

### B5. Fix the reverse-pass quirk → true symmetric (SSOR-like) sweeps

The "reverse" sweep actually runs forward (flagged in `gs_sweep!` docs).
R4 runs reverse_pass=false so nothing is wrong today, but a *true* backward
sweep enables symmetric GS, which often converges in fewer outer iterations
per byte streamed. Needs recalibration; a cheap convergence experiment can be
run locally at ≤4 threads before any campaign. Rank: opportunistic.

## Projection vs the 10.116 s operating point (27 iterations held fixed except B1)

| stack | products s/it | chain+overhead s/it | projected wall | vs 10.116 s |
|---|---|---|---|---|
| status quo colored@j16 | ~0.12 (≤74 GB/s) | ~0.31 + barriers | 10.116 (measured) | 1.0× |
| A1+A2 @ j16 | ~0.03 | ~0.13 | ~4.5-5.5 s | ~2× |
| A1+A2+A3+A4+A5 @ j64 | ~0.026 | ~0.10-0.12 | ~3.5-4 s | ~2.5-3× |
| + A8 (both sockets) | ~0.013 | ~0.09 | ~3-3.5 s | ~3× |
| + A9 instead of colored | same | less barrier waste | ~3 s | ~3.3× |
| + B1 (27→~15 iter) | — | — | ~2 s | ~5× |

Honesty notes: (i) nothing in Tier A can beat 5.6×/socket on the products
term — the table respects that; (ii) the projections assume the real kernel
reaches the microbenchmark's bandwidth, but v22's real kernel ran at ~37 GB/s
aggregate vs the analogue's 74, so scatter/interleaving overheads may tax
these numbers — treat the wall column as best-case; (iii) A6 and B1 are the
only entries whose upside is not already bounded by the bandwidth arithmetic.

## Recommendation (top 2)

**First: A1+A2 together (parallel sweep-affine first-touch assembly +
Float32), with A3+A4 riding along** — they are one coherent code change in
FastMultipole's solve path, all evidence-backed (164 GB/s measured, bytes
halve by construction), A1/A3/A4 are bit-identical and A2 is already
Ryan-approved; combined honest projection ~2-2.5× (10.1 → ~4-5 s). **Second:
A6, the staircase re-tune at the new operating point** — zero code risk, uses
the existing calibration protocol, and is required anyway because leaf=100 /
inner=3 were optimal for a serial 29 GB/s sweep that will no longer exist. I
would hold A9 (DAG) and B1 (Krylov) as the two big follow-ons and pick
between them based on where the post-A1/A2 profile says the time went
(barriers → A9; iteration count → B1).
