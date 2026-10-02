# Split source/target-major FGS cache and follow-on ideas (2026-09-18)

Companion to `thread_efficiency_proposals_20260918.md`. This note records the
follow-up discussion about storing the FGS near-field operator in both source-
and target-oriented forms. It is analysis only: no solver code was changed and
no jobs were submitted.

## Why the current cache is source-major

`nonself_influence_matrices` sorts the direct list by source leaf and packs one
dense matrix per source leaf. For source leaf $j$,

$$
M_j = \begin{bmatrix} A_{i_1j} \\ A_{i_2j} \\ \vdots \end{bmatrix},
$$

so `compute_nonself_products!` reads the contiguous strength block $x_j$ once
and one tall, narrow GEMV produces its influence on all neighboring targets.
`scatter_nonself_influence!` then applies those target segments to the global
RHS. This organization has three important benefits:

1. R4 performs 1,068 large GEMVs per sweep, rather than one GEMV for each of
   the 95,390 directed leaf interactions.
2. Each GEMV uses a contiguous source vector; the median matrix is about
   $4708\times39$ and 1.40 MB, a favorable non-transposed GEMV shape.
3. Applying every source delta immediately leaves the complete near-field RHS
   current at the end of each sweep, so the next outer residual requires no
   extra matrix pass.

The cost is shared target writes. Products from different source leaves write
overlapping RHS rows, which creates the serial scatter and the coloring/chunking
constraints.

## Why a pure target-major pull is not automatically better

A viable FGS target-major cache would aggregate by *solve leaf*, rather than
copying the existing Krylov near-field cache's one-block-per-direct-entry
layout. For target leaf $i$ it would store

$$
T_i = \begin{bmatrix} A_{ij_1} & A_{ij_2} & \cdots \end{bmatrix},
$$

normally transposed so each target-row dot product reads contiguous data. At
leaf $i$, the global strength vector already contains the lexicographic GS
mixture: leaves $j<i$ have new values and leaves $j>i$ retain old values.
Pulling all neighboring strengths therefore gives the same mathematical
block-GS update.

The R4 census makes this plausible:

| Quantity | Source-major | Aggregated target-major |
|---|---:|---:|
| GEMVs per sweep | 1,068 | about 1,068 |
| Matrix bytes | 2,862,850,032 | unchanged |
| Median block bytes | 1,401,640 | 1,558,288 |
| p90 block bytes | 4,136,700 | 4,472,640 |
| Maximum block bytes | 80,399,400 | 70,620,800 |
| Median small dimension | 39 source strengths | 39 target rows |

Gathering all target-neighbor strengths would write only 47,082,944 bytes per
sweep, 1.64% of the matrix stream. Only 1.95% of matrix bytes belong to target
leaves with fewer than 16 rows, and 10.05% belong to leaves with fewer than 32
rows, so most bytes expose useful intra-leaf parallelism at j16-j32.

The governing disadvantage is final-state maintenance. After target $i$ is
solved, later sources change, making its pulled RHS stale. A full residual pass
would add one 2.86 GB stream per outer iteration, raising three inner-sweep
streams to four. Correcting only interactions with $j>i$ still rereads 47.14%
of the cache after the final inner sweep, a 15.7% matrix-byte penalty averaged
over three sweeps. Pure pull also changes floating-point grouping, needs a
transposed-GEMV kernel, and wants page striping within the one active target
block rather than source-affine placement.

## Preferred form: split dual layout

The best form partitions interactions by their position in the lexicographic
order instead of evaluating two complete copies:

| Interaction | Stored/evaluated as | R4 matrix-byte share |
|---|---|---:|
| $j<i$: earlier source to later target | target-major pull immediately before solving $i$ | 52.86% |
| $j>i$: later source to earlier target | source-major delta push after solving $j$ | 47.14% |

Maintain an upper/backward contribution accumulator for every target leaf.
During a sweep:

1. Before solving leaf $i$, pull the $j<i$ half using current strengths and
   combine it with the maintained $j>i$ contribution and the fixed external +
   far-field RHS.
2. Solve leaf $i$ and form $\Delta x_i$.
3. Enqueue the source-major product from $\Delta x_i$ to targets $k<i$.
   Those targets will not be solved again in this sweep, so this work is off
   the current sweep's dependency path.
4. At the sweep boundary, finish the backward products and apply them by
   target in ascending source order. The updated accumulator is ready for the
   next sweep, and the full RHS is current for the outer residual.

Every directed interaction is read exactly once per sweep. The forward 52.86%
is target-owned and race-free on the critical path; the backward 47.14% is
source-efficient and can run concurrently or as one parallel phase before the
single sweep boundary. This preserves lexicographic block-GS mathematics and
does not introduce chunked mode's cross-chunk Jacobi lag.

### Storage choices

| Representation | R4 matrix storage | Use |
|---|---:|---|
| Current source-major Float64 | 2.86 GB | baseline |
| Two complete Float64 layouts | 5.73 GB | simplest prototype; half of each copy is unused during a split sweep |
| Two complete Float32 layouts | 2.86 GB | prototype with the same matrix footprint as today's Float64 cache |
| Split triangular Float64 layouts | 2.86 GB | production form; each interaction stored once |
| Split triangular Float32 layouts | 1.43 GB | lowest-traffic production form |

The full dual layout is useful for establishing correctness and timing because
it avoids an immediately invasive constructor rewrite. The split layout is the
logical production representation once the schedule is proven.

### Advantages

- Same mathematical lexicographic block-GS iterate; no deliberately stale
  cross-chunk coupling and no expected 27-to-44 iteration inflation.
- One matrix read per interaction per sweep and no extra residual pass.
- Keeps large aggregated GEMVs in both orientations instead of creating
  95,390 small pair GEMVs.
- Removes forward shared-RHS writes: each target pulls and owns its current
  rows.
- Moves nearly half the cache stream off the current sweep's leaf-solve
  dependency path.
- Retains contiguous source strengths and source-affine placement for the
  backward half.
- Allows different placement policies for the two consumers: intra-block
  striping for target-major blocks and source-worker affinity for source-major
  blocks.
- Float32 makes a complete two-layout prototype cost no more matrix memory
  than the present single Float64 cache.

### Disadvantages and risks

- It is not bit-identical to the current incremental source-major path because
  the forward sum is regrouped. It needs the normal tolerance staircase and
  accepted-accuracy certification.
- Direct-list target branches can span multiple solve leaves. Assembly must
  split their row ranges at solve-leaf boundaries before classifying blocks as
  $j<i$ or $j>i$.
- Two matrix builders, maps, and execution kernels replace one simple
  source-major container. Transform/reuse and multi-system paths must carry
  both representations correctly.
- Target-major execution naturally uses a transposed GEMV or long-dot kernel;
  its Zen3 bandwidth is not established by the existing non-transposed GEMV
  measurements.
- Backward products may be dependency-free, but they still consume memory
  bandwidth. Overlap helps only when the critical target pull leaves memory
  controllers idle or the work is divided across sockets.
- Deterministic backward application requires target ownership and ascending
  source order. Concurrent direct writes would race and change summation order.
- A complete dual Float64 prototype raises retained R4 storage by about 2.86
  GB and raises constructor peak memory further. This is comfortable on the
  benchmark node but affects the campaign's smaller memory classes and scales
  poorly at larger rungs.

## Next five ideas, ranked by expected value

The ranking favors general CPU improvements that retain the FGS algorithm;
the fifth item is the strongest algorithm-changing fallback.

### 1. Float32 cache storage plus consumer-aligned first touch

Use Float32 for the nonself operator while retaining Float64 strengths,
accumulators, self matrices, and leaf LU factors initially. Build the
target-major pages with the thread team that consumes target blocks and the
source-major pages with their source-worker owners. Float32 halves the dominant
stream and also makes a complete dual-layout prototype fit in today's Float64
matrix footprint. The numerical operator changes slightly, so the standard
accuracy staircase remains mandatory. Expected value: highest; it directly
multiplies every scheduling improvement and attacks the measured bottleneck.

### 2. Persistent adaptive threaded GEMV on the existing source-major layout

As the lower-risk alternative to the split redesign, retain the current
source-major cache and parallelize within one large source GEMV using a
persistent worker team. Run small blocks serially, split medium blocks across
4-16 workers, and tile the largest blocks more widely; stripe each participating
block's pages across the same workers' NUMA nodes. This preserves the current
incremental RHS and iteration exactly except for parallel reduction details.
It avoids 79 color barriers, although roughly 1,068 leaf-level handoffs per
sweep and small-block efficiency may limit it. A shape-weighted GEMV
microbenchmark can reject or validate it before a solver rewrite.

### 3. Use both sockets after the single-socket schedule is efficient

The second socket roughly doubles the memory-controller ceiling from about
164 GB/s toward 300-330 GB/s. The split layout offers a natural division:
target-pull work can use a compact team while backward source products run on
the other socket, with target-owned application at the boundary. Cross-socket
coherence, remote pages, and the serial leaf-solve chain limit the likely
whole-solve gain to tens of percent after Float32, but this is the remaining
hardware multiplier once one socket is saturated.

### 4. Jointly retune leaf size, MAC, expansion order, and inner sweeps

The retained leaf=100, MAC=0.4, P=8, inner=3 point was selected for the current
serial, Float64 cost balance. After precision and scheduling changes, retune
the axes jointly. In particular, test whether a higher expansion order permits
a looser MAC that moves marginal direct blocks out of the streamed cache while
still satisfying the BC gate; P=10/MAC=0.4 and P=8/MAC=0.5 separately do not
answer that joint question. This retains the solver family and may remove more
bytes than low-level tuning can accelerate.

### 5. Anderson or Krylov acceleration around the FGS sweep

If the optimized sweep remains limited by 81 complete cache streams, reduce
the number of streams. Anderson acceleration of the outer fixed point or
FGMRES with one or a few FGS sweeps as a preconditioner could plausibly reduce
27 outer iterations to roughly 12-18. This changes the algorithm and requires
new convergence, memory, and accepted-accuracy tuning, but it is the strongest
general fallback because it multiplies every kernel and placement improvement.

## Current recommendation

Treat the split dual layout as the leading large algorithm-preserving design,
with a full dual Float32 representation as the simplest experimental form.
Before committing to the full rewrite, compare the actual R4-shaped
non-transposed source GEMVs against aggregated transposed target GEMVs under
equal precision, worker count, and NUMA placement. If target pull cannot beat
or at least match source-product bandwidth, retain source-major storage and
pursue idea 2; if it does, the triangular split removes pure pull's residual
penalty and provides a credible path to one synchronization boundary per
sweep.
