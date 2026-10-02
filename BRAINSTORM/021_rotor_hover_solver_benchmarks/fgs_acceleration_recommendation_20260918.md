# FGS acceleration: recommendation after code review (2026-09-18)

## Recommendation and required result

My first choice is **persistent, adaptive row-parallel execution of the existing source-major cache, with consumer-aligned memory placement**, followed by **Float32 nonself matrix storage with Float64 arithmetic and leaf factors**. Keep the lexicographic leaf sequence, current far-field refresh schedule, and current RHS semantics. Fuse product bookkeeping and scatter where safe. This attacks the dominant cost without changing the block-GS iteration or rebuilding the interaction layout.

The split triangular source/target layout is a sound second candidate. The newer **pull-DAG** proposal improves it enough to justify an early graph/schedule audit alongside the source-kernel experiment, rather than waiting for source-major execution to fail. Its mathematical argument holds; its performance advantage is unproven. Source-major remains my provisional first implementation because it has fewer moving parts. Select between the designs using the complete dependency-constrained schedule, not isolated pull-versus-push bandwidth. Do not begin with a general source-scatter dependency-DAG scheduler.

The acceptance baseline is the best existing accepted R4 result, **colored @ j16 = 10.116 s**, not the slower serial-order reference. “50% speedup” can mean either 1.5× throughput or half the time. Require **≤6.744 s** for the former, and design for the stronger **≤5.058 s** target. A credible planning range for the proposed mixed-storage implementation is **4.5–5.9 s**, conditional on measured kernel bandwidth and bounded synchronization cost. These are projections, not demonstrated results. No implementation or new timing runs were performed for this note.

## Evidence and scope

Read and independently assessed [thread-efficiency proposals](thread_efficiency_proposals_20260918.md), [split-layout ideas](fgs_split_dual_layout_ideas_20260918.md), and, in this revision, [the top-five synthesis](thread_efficiency_top5_20260918.md). Initial code inspection used FLOWPanel HEAD `4e434cb1c2bf94e0fa9beaf4be440f6d1f4a5c97`; the follow-up found HEAD `a6270c8de0683bb9ab66fad11b1ac37048e5cb74`. The Manifest's sibling FastMultipole checkout remains at `c18e4b46116825a201a0cf459b7596242191dcae`; its relevant dependency, residual, and container code was rechecked. Working trees contain unrelated edits; these are inspection identifiers, not clean campaign pins. Only this recommendation document was edited for the revision.

| Evidence | Finding | Qualification |
|---|---|---|
| [v21 accepted A/B](fgs_r4_followup_evidence_20260914/colored-v21-13738665/analysis/ab_summary.md) | Colored j16: 10.116 s, 26 iterations; lex: 27 iterations | Colored and lex use separately calibrated tolerances and different leaf orderings. |
| [v22 accepted A/B](fgs_r4_followup_evidence_20260914/chunked-v22-13749231/analysis/ab_summary.md) | Lex j64: 10.951 s; chunked j64: 14.264 s, 44 iterations | Chunking reduced per-iteration nearfield time but lost overall. |
| [NUMA diagnostic](numa_placement_findings_20260918.md), [driver](../../benchmark/numa_dgemv_bench.jl) | 29.4 GB/s serial; 74.0 with serial first touch and parallel consumption; 164.4 with affine placement | Synthetic 650×550 blocks, not the real highly rectangular distribution. The in-situ placement comparison was not run. |
| [Retained configuration](../../benchmark/retained_r4_diagnostics.toml) | P8, MAC0.4, leaf100, inner3, lex, LU cached | Model below uses 27 completed outer updates / 81 sweeps; do not assign this count to colored. |
| [Saved source census](fgs_r4_followup_evidence_20260914/diag-v15-13694724/j64-b1/results/gemv_census.csv) and [directed edges](fgs_r4_followup_evidence_20260914/diag-v15-13694724/j64-b1/results/dependency_edges.csv) | 1,068 leaves, 2,862,850,032 matrix bytes per sweep; lower/upper shares 52.8597% / 47.1403% | Reconstructed coefficient bytes from leaf-pair sizes exactly match the cache total. The triangular shares are now independently verified; target-kernel performance is not. |

The scope is prepared cold body solves at R4. Constructor improvements must be reported separately. Benefits to an entire unsteady simulation depend on its body-solve fraction and warm-start behavior; they cannot be inferred from these times alone.

## Top five improvements, ranked

Ranked by expected benefit relative to implementation risk, with the cheaper evidence-gathering gates allowed to run before that implementation order. Numbers refer to the quantitative model below, not new measurements.

| Rank | Improvement | Why / condition |
|---|---|---|
| 1 | Persistent source-row workers, matching page placement, fused old/new scatter | First implementation candidate; preserves lexicographic mathematics. Float64 at an effective 80 GB/s projects 5.90 s before added overhead. Reject if leaf handoffs consume the budget. |
| 2 | Float32 nonself coefficients, Float64 arithmetic and self/LU | Multiplies either schedule's benefit. Combined with rank 1, 40–80 useful GB/s projects 5.90–4.45 s before added overhead; must pass independent true-accuracy evaluation. |
| 3 | Split triangular pull-DAG, adaptive teams, deterministic backward reduction | Strongest larger redesign. Audit the actual lower dependency graph early; promote ahead of rank 1 implementation if inclusive scheduling evidence supports it. NUMA-local work queues and a two-socket variant are extensions, not assumed gains. |
| 4 | Joint leaf/MAC/P/inner retuning after the cost balance changes | Removes work or improves its balance; keep accepted accuracy fixed. Changes configuration/operator approximation, so report separately from same-iterate engineering. No defensible speedup estimate yet. |
| 5 | Safeguarded Anderson, or inexpensive FGS-preconditioned FGMRES | **Algorithm-changing; optional separate track.** Best route beyond the bandwidth ceiling if it reduces total work. Charge state reconstruction, residual evaluations, and preconditioner applications; 27→12–18 updates is unverified. |

Ranks 1–3 preserve the mathematical lexicographic iteration for a fixed operator; rank 2 perturbs that operator's coefficients. Rank 4 retains the FGS family but changes its operating point. Rank 5 and the other algorithm-changing alternatives are developed separately below. No individual item is an unconditional 50% win.

## What the implementation implies

The following references are to FastMultipole `src/solve.jl` and `src/containers.jl` in the sibling checkout, except where noted.

1. `gs_sweep!` (`solve.jl:1208–1302`) solves a leaf, evaluates its nonself product, scatters, and only then proceeds in lexicographic mode. The product has many independently writable rows even though leaf solves are ordered. Parallelizing those rows preserves the dependency sequence.
2. `compute_nonself_products!` (`:923`) computes **the full new product** $A_jx_j$, after saving the old product. `scatter_nonself_influence!` (`:948`) performs `rhs += old`, then `rhs -= new`. The current implementation does not multiply $A_j\Delta x_j$. Replacing it with a delta product or `rhs += old - new` changes floating-point behavior.
3. `nonself_influence_matrices` (`:143`, fill at `:256–321`) builds one tall matrix per source leaf, serially. Target branches may contain multiple solve leaves. Do not equate a direct-list entry with one solve leaf when constructing ownership maps.
4. `Matrices{TF}` (`containers.jl:1092`) uses the same type for matrix data and product/RHS storage. `FastGaussSeidel` also couples self and nonself matrix precision. Mixed storage needs an explicit container/type change, not a Float32 conversion followed by an assumed fast mixed-type BLAS call.
5. `residual!` (`solve.jl:1803`) returns max-abs, but every leaf uses `view(residual_vector, 1:length(rhs))`. A threaded leaf loop would race. Use disjoint or worker-private scratch first. Its unused MSE bookkeeping is not the returned convergence criterion.
6. The outer loop (`:1380–1537`) initializes nonself contributions from the actual starting strengths, refreshes far-field influence, checks residual, then performs inner sweeps. There is a final convergence-check FMM pass: v22 records 28 FMM calls for 27 updates. FLOWPanel explicitly sets `final_update=false` (`src/FLOWPanel_solver.jl:1814–1824`), so removing a final direct evaluation is not an available body-solve optimization.

“Algorithm intact” below means the same mathematical block-GS update with the same leaf order and far-field refresh. It does not imply bitwise equality after changing matrix precision or dot-product grouping. I distinguish those changes from new iteration methods below.

## First implementation candidate: source-major row parallelism

### Execution and data placement

Retain the existing cache and LU solves. For each leaf, the coordinator solves its diagonal block and publishes its updated strengths. A persistent team computes disjoint row tiles of that source's tall matrix. Each tile retains its old product, computes its new product, and applies the two separate RHS operations to its owned target rows. Complete all source updates before the next leaf solve.

The direct-list decomposition gives disjoint target segments for one source; establish this explicitly when building the row-to-target map, including branch spans and multiple systems. If a generalized input can contain overlapping segments, assign all occurrences of each target row to one owner and retain their original order. Different sources do overlap, but only one source is active in this first design.

Use contiguous row tiles **within each column**, maintaining column-major streaming. Do not implement one strided scalar dot product per output row. A worker can accumulate a tile in private scratch, then scatter it. This keeps the old product live locally and can remove the full old-product copy on this path. Other sweep modes can retain the old storage until separately adapted. Avoid allocations and task creation inside the leaf loop.

The worker count should depend on block bytes and dimensions: tiny blocks serial, medium blocks on a compact team, large blocks across more memory controllers. Fixed physical worker ownership and page placement must match. A source-affine fill designed for colored mode is **not** the right placement for a block consumed by several row workers. Interleave as a simple experimental control; use page-aware row tiling for the eventual owner mapping. Column boundaries, sub-page tiles, and huge pages can defeat fine-grained ownership, so verify actual page placement rather than assuming it.

Parallel probe-based assembly needs private mutable buffers, as in the separate `nearfield_cache.jl:438–478` builder. An alternative first experiment is a parallel repack into new storage consumed by the worker team; include its peak memory and fresh-solve cost. The steady-state measurement must not hide that cost when reporting constructor-plus-solve performance.

### Why this is my preferred first test

It retains large aggregated products, contiguous source strengths, one matrix stream per sweep, and a current RHS for the existing residual. It exposes thousands of output rows on typical source blocks, rather than roughly 39 outputs on a typical transposed target block. It avoids the 79-color schedule and does not introduce cross-chunk Jacobi lag. Fusing bookkeeping/scatter may also reduce the serial remainder, though the model below does not assume this win.

Its weakness is **1,068 leaf handoffs per sweep**, or 86,508 over 81 sweeps. Just 10 microseconds of extra latency per handoff costs 0.865 s per solve. To keep added coordination below 0.5 s, the average must stay below about 5.8 microseconds. “Persistent” does not make synchronization free. Use a realistic sequence benchmark to reject this design quickly if the handoff budget fails; simply turning on multithreaded BLAS for every small product is not an adequate implementation.

Preserving row arithmetic order is possible in a custom kernel, but reproducing the existing BLAS result bit-for-bit is not automatic: tiling, SIMD, and FMA choices can change rounding. Check it explicitly; otherwise classify the result as mathematically equivalent and certify its accepted accuracy.

### Add mixed storage, retaining sensitive work in Float64

Store only nonself coefficients in Float32 initially; keep strengths, products, RHS, residual, self matrices, and cached LU in Float64. Convert coefficient tiles on load and multiply/accumulate in Float64. Do not round strengths to Float32 merely to get an `sgemv` call: that is a different numerical experiment.

This halves coefficient traffic, but may introduce conversion/compute limits. The required benchmark is the actual mixed kernel, not a Float64 memory-copy rate divided by two. Separate coefficient storage type from accumulator type in the container and dispatch. If some blocks are precision-sensitive, retain those blocks in Float64; the remaining traffic fraction becomes $1-f/2$, where $f$ is the fraction of original coefficient bytes converted to Float32.

Float32 storage perturbs the operator. Float64 accumulation cannot recover discarded coefficients, and the internal residual may be small for the rounded operator while the true BC residual fails. Preserve the independent acceptance evaluator. If precision fails, retain Float64 or introduce selective precision; do not silently relax the accuracy gate. Mixed-precision iterative refinement is a separate algorithm option below.

## A conservative speedup budget

The companion profile gives approximately 0.292 s of products per outer update and 0.032 + 0.020 + 0.06 s for scatter, leaf solves, and other work. This supports an approximate **3.0 s non-product allowance per solve**. It is a planning allowance inferred from the existing profile, not an irreducible measured floor for a new implementation.

For 81 sweeps, the coefficient stream is 231.891 GB in Float64, or 115.945 GB in Float32. Let $B$ be useful coefficient GB/s across the actual sequence, and $H$ be additional coordination/initialization costs beyond that allowance:

$$
T_{64} \approx 3.0 + \frac{231.891}{B} + H,\qquad
T_{32} \approx 3.0 + \frac{115.945}{B} + H.
$$

Use inclusive kernel timings when estimating $B$, or put excluded scheduling time in $H$, but do not count it twice. Initial nonself setup also streams the operator once, and stopping adds the final FMM check; measure those explicitly in the final budget rather than treating 81 as the full solve's exact operation count.

| Candidate | Assumed useful GB/s | Time with $H=0$ | Speedup vs 10.116 s |
|---|---:|---:|---:|
| Float64 row workers | 60 | 6.865 s | 1.47×: fails minimum |
| Float64 row workers | 80 | 5.899 s | 1.71× |
| Float32 storage / Float64 arithmetic | 40 | 5.899 s | 1.71× |
| Float32 storage / Float64 arithmetic | 60 | 4.932 s | 2.05× |
| Float32 storage / Float64 arithmetic | 80 | 4.449 s | 2.27× |

At $H=0$, the minimum requires 61.9 GB/s in Float64 or 31.0 GB/s in mixed storage. Halving wall time requires 112.7 or 56.3 GB/s respectively. With $H=0.5$ s, mixed storage needs about **74.4 GB/s** to halve wall time. This is why the lower-traffic option is valuable and why handoff cost is a go/no-go gate.

Do not multiply a 2.22× placement gain by a 5.6× threading gain: those ratios share the same bandwidth ceiling. Also, the real serial profile's 0.292 s is approximately 8.59 GB / 29.4 GB/s; reaching 164 GB/s would be about 5.6× on that term, not merely 2.22×. The latter ratio compares two already-parallel synthetic arms. Neither establishes whole-solve performance.

## Independent assessment of the split layout

Write the nearfield matrix as $D+L+U$ in lexicographic leaf order, with the far-field term frozen during the inner sweeps. The sweep is

$$
D_i x_i^{s+1}=b_i^{k}
-\sum_{j<i}A_{ij}x_j^{s+1}
-\sum_{j>i}A_{ij}x_j^{s}.
$$

The proposed lower-triangular target pull computes the first sum. An upper accumulator $u^s=Ux^s$ supplies the second. Backward source products can be buffered while the sweep proceeds, then applied at its boundary to obtain $u^{s+1}$. **They must not overwrite the frozen accumulator while any target still needs it**; separate pending buffers or versions make this explicit. At the boundary, the saved lower sums are already $Lx^{s+1}$ because all their sources were solved earlier. Combining them with $u^{s+1}$ yields a complete nearfield RHS without another matrix pass. The mathematical claim is valid.

Implementation must also initialize $Ux^0$ for nonzero warm starts, keep fixed external RHS separate from the current outer farfield, preserve the sign convention, split branch ranges at solve-leaf boundaries, and map all systems consistently. Actual reverse-order sweeps would interchange the two triangles. The current reverse flag repeats forward order; do not silently change that behavior as part of this optimization.

My reservations are about speed and development cost:

- The reported 47.14% backward stream leaves the immediate dependency path but does not disappear. It competes for bandwidth with the forward stream; overlap only helps while resources would otherwise be idle.
- The lower solve chain remains. A threaded pull still needs publication/completion at each leaf. “One sweep boundary” describes the deferred backward phase, not all synchronization in the method.
- Typical target blocks expose about 39 output dots. Parallelizing one block across more workers requires dot-product reductions, changes grouping, and adds coordination. The source layout already exposes thousands of rows. A pull-DAG can instead process multiple ready target blocks concurrently; therefore the 39-row limitation matters most when the ready frontier is narrow.
- A full dual Float32 prototype combines layout and precision changes, making a failure hard to attribute. Establish the Float64 algebra on small fixtures first; then compare same-precision, same-placement real-shape kernels. Avoid allocating two full production-size Float64 layouts merely to prove the recurrence.

Promote the split design if its lower-pull plus batched-upper schedule wins on an inclusive sequence benchmark with sufficient margin to justify the extra implementation cost; source-major need not fail first. Prefer split triangular storage for production; full duplication has no sustained byte advantage. Pure target pull is lower priority because it additionally needs a complete or partial final-RHS refresh.

### Verified improvement from the top-five note: pull-DAG execution

For each solve leaf $i$, define $\mathcal P_i=\{j<i:A_{ij}\text{ is a retained direct interaction}\}$. Wait until all leaves in $\mathcal P_i$ have published their new strengths, then evaluate its aggregated lower pull and diagonal solve. Independent ready leaves may execute out of index order: by induction, each reads exactly the lower new strengths and frozen upper contribution required by the lexicographic recurrence. This preserves the mathematical iterate even when execution order differs. Build these **directed lower edges** from the actual source/target ranges. `color_leaves` (`solve.jl:1126–1160`) instead symmetrizes conflicts for coloring; its adjacency is not the exact pull graph.

This removes the forward shared-RHS scatter and its per-row source-order scheduling burden. It is a meaningful simplification over a source-push DAG, and my original blanket recommendation to defer DAG work was too broad. It does not remove publication/dependency synchronization, the backward reduction, or the lower graph's critical path. A counter must observe completed strength writes before making a target ready. Determinism requires fixed per-row arithmetic/reduction order, not just ownership; a thread-dependent BLAS reduction can still change results.

The note's 38.8 active chunked threads do **not** establish a broad exact-GS frontier. The inspected chunked code (`solve.jl:1244–1274`) defers cross-chunk dependencies; the [v22 activity report](fgs_r4_followup_evidence_20260914/chunked-v22-13749231/analysis/ab_summary.md) measures that different schedule. A simple counterexample is a lower chain $1\to2\to\cdots\to n$: splitting it into Jacobi-coupled chunks creates parallel work, while its exact lower solve has width one. Moreover, CPU activity includes spinning. Neither activity nor the number of colors measures the directed critical path.

**The saved R4 edges permit a stronger check than that counterexample.** The [census generator](../../benchmark/fgs_r4_diagnostics.jl) emits directed `writes_to` edges before symmetrizing its separate coloring graph. Filtering `source_leaf < dependent_leaf` therefore gives the required lower dependencies. Read-only CSV analysis, independently repeated in this review, gives:

| R4 structural quantity | Recomputed result |
|---|---:|
| Directed nonself leaf pairs / lower edges | 95,390 / 48,167 |
| Longest lower path, counting solve leaves | 279 of 1,068 leaves |
| Largest earliest-level cohort | 6 leaves |
| Equal-cost work / longest path | 3.828 |
| Lower coefficient bytes per sweep | 1,513,294,744 |
| Upper coefficient bytes per sweep | 1,349,555,288 |
| Longest path weighted by each target's lower-product bytes | 531,147,048 bytes |
| Lower byte work / byte-weighted longest path | 2.849 |

For reproducibility, read each leaf's `n` from `gemv_census.csv`; assign a directed pair $(j,i)$ the weight $8n_jn_i$. The sum over all edges is exactly 2,862,850,032 bytes. Sum incoming lower-edge weights to obtain each target's pull bytes, then use the recurrence below with those bytes as task weights. The unit-weight recurrence gives 279 levels. No solver construction or numerical benchmark was needed.

The six-leaf cohort is specific to earliest-level scheduling, **not a proof that every possible ready set is at most six**. The stronger product-only bound is the 2.849 work/span ratio: with one worker per pull and uniform seconds per coefficient byte, average lower-product parallelism cannot exceed it, even with unlimited workers and zero scheduler cost. Actual kernel/solve costs require reweighting. This substantially weakens the case for many independent one-worker pulls; it does not rule out socket saturation using larger teams on critical pulls plus backward filler work. The adaptive-team experiment is therefore central, not an optional refinement.

Before choosing this design, replace the structural weights above with measured per-leaf lower-product/solve costs. With one worker per task, compute

$$
C_i=t_i^{\rm pull}+t_i^{\rm solve}+\max_{j\in\mathcal P_i} C_j,
\qquad C=\max_i C_i,
$$

where the empty maximum is zero. Also measure ready-frontier width, byte-weighted work, edge-publication overhead, and the release times/drain tail of backward tasks. A sweep cannot beat its critical path, worker-work bound, or memory-traffic bound. Recompute task costs if narrow frontiers trigger intra-leaf teams. A toy scheduler must replay the actual edges and variable costs; a wide synthetic DAG would not answer the question.

Use a hybrid policy: one/few workers per pull when many leaves are ready, a larger row/dot team when the frontier narrows, and backward work as lower-priority filler. Keep critical pulls from waiting behind long backward blocks. Work stealing should start within each NUMA domain; unrestricted stealing conflicts with consumer-affine first touch. These policies need inclusive timing, not assumptions about free scheduling.

The proposed gate “pull bandwidth ≥ push bandwidth ⇒ build split” is **neither sufficient nor necessary**. A fast pull can lose to a long dependency chain or a backward drain tail; a somewhat slower pull can still win by eliminating full-team handoffs and overlapping independent work. Compare complete sweeps and projected accepted solve time.

### Further improvement: rebuild the upper contribution from full products

For the split design, consider computing $q_j=U_{:,j}x_j^{s+1}$ once after source $j$ finishes, then forming $u^{s+1}=\sum_j q_j$ in fixed source order into a separate next-sweep buffer. This is an alternative to accumulating $U\Delta x$ into the previous upper sum. It reads the upper matrix once per sweep either way, needs no saved old strengths for the backward product, and avoids carrying incremental upper-sum roundoff across sweeps. The current source-major path already computes full new products, so this choice is closer to its kernel semantics.

Keep $u^s$ immutable while current pulls consume it. Source-private product segments and a target-owned ordered reduction can produce $u^{s+1}$; publish it only after all work for the sweep completes. Budget the output buffers and reduction traffic. This does not promise bit identity with the old scatter or eliminate floating-point error, but it offers a simpler state invariant and a useful robustness improvement at the same coefficient-stream count.

### Claims from the top-five note that are not established

| Claim | Verification / correction |
|---|---|
| Float32 + first touch gives ~2× “unconditionally” | Container types still couple matrix/RHS precision; actual mixed-kernel speed and true-operator accuracy remain unmeasured. Retain the numerical and schedule gates. |
| Threaded residual is a free exact rider | Rechecked `residual!`: all leaves reuse the same scratch prefix. Private/disjoint scratch is required before threading; max reduction alone does not make it safe. |
| Split halves give a coherence-free two-socket decomposition | Matrix streams can be disjoint, but updated strengths and per-source readiness must cross sockets during the sweep. End-of-sweep accumulators also cross ownership. Pad/batch publications and measure traffic/latency; the sockets do not meet only at the boundary. |
| ~330 GB/s is a measured two-socket ceiling | The cited NUMA experiment measured up to 164.4 GB/s on one socket; ~330 is an extrapolation. No two-socket result is supplied by that evidence. Nearly equal matrix bytes also do not imply balanced times: lower dependency stalls and different kernels matter. |
| 1,068 handoffs versus 158 color phases proves a 6.8× synchronization penalty | The count ratio is correct; their latency, team size, and useful work differ. Benchmark actual handshakes. For a pull-DAG also count edge notifications and queue operations. |
| Dual Float32 is as cheap as current Float64 | Coefficient capacity is equal for two complete copies. Product buffers, graph metadata, scratch, and construction peak memory are extra. |
| Joint retuning has nil risk | New MAC/P/leaf/inner combinations can fail accuracy or convergence. Calibration controls that risk; it does not guarantee a faster accepted point. |

The projected 3–3.5 s also needs an explicit remainder improvement. Using this note's unchanged 3.0 s allowance, even 164.4 GB/s mixed-storage streaming gives about **3.71 s**, before added overhead. At a hypothetical 330 GB/s it gives about **3.35 s**. A single-socket 3–3.5 s result is possible only if the new schedule also reduces other costs; it does not follow from coefficient bandwidth alone. Keep it as a stretch target rather than replacing the more conservative budget above.

## Disposition of the other existing proposals

| Idea | My assessment |
|---|---|
| Parallel first-touch assembly | Accept, but ownership must match the selected consumer. Placement alone cannot parallelize a lexicographic GEMV. Constructor speed is separate from prepared solve speed. |
| Parallel target-owned scatter | Accept within a dependency-safe phase. For colored mode preserve **color-major, then ascending source within that color** order. Globally sorting all sources or deferring all colors to sweep end does not preserve that iteration. Two separate old/new RHS operations matter for bit identity. |
| Thread residual and vector work | Secondary. First fix shared residual scratch; preserve nonfinite handling. The whole 0.06 s remainder is not residual time, so do not promise 0.02–0.03 s savings without attribution. |
| Colored scheduling improvements | Useful fallback/control: persistent workers, serial tiny colors, cost balancing with stable ownership. Reordering work inside a color is safe only if scatter retains its defined order. Useful by itself, but current results do not establish a 50% gain. |
| Exact-GS dependency DAG | Distinguish two designs: defer the general source-scatter DAG, but audit the simpler split pull-DAG early. The latter removes forward scatter ordering, not the lower critical path or backward reduction. No basis yet for guaranteed socket saturation. |
| Huge pages | Opportunistic after placement. A typical block can be smaller than 2 MB; coarse pages can cross owners. Measure TLB misses and placement rather than budgeting a guaranteed percentage. |
| Second socket | Test after the one-socket kernel. More controllers can help, but source-level handshakes and remote data can erase the gain. It is more than a launcher change unless ownership/pinning are already correct. |
| Joint leaf/MAC/P/inner retuning | High-value follow-on, separately reported from exact-iterate improvements. It changes the approximate operator, partition, or iteration schedule. A looser MAC needs independent accuracy checks; a larger expansion order may compensate but is not guaranteed to do so. |

After mixed storage, a 3 s remainder limits further gains. At 80 GB/s, even free products would improve the projected 4.45 s by only 1.48×. Therefore pursuing much beyond 2× overall should target the remainder, reduce cache size at accepted accuracy, or reduce iteration work—not just add bandwidth.

## Algorithm-changing options, ranked separately

1. **Safeguarded Anderson acceleration of the outer FGS map.** Start with a small history (for example 3–8), restart on ill-conditioned history, and accept extrapolated iterates only under a defined residual safeguard. The important integration cost is rebuilding consistent nonself and farfield state after mixing; mixing strengths while retaining stale incremental products is wrong. This can add an operator pass, so compare total streams and FMM calls, not just iteration count. The claimed 27→12–18 is a hypothesis, not a forecast. I rank this first because it can reuse the optimized sweep without nesting a long solver inside every Krylov step.
2. **FGS-preconditioned FGMRES with deliberately cheap fixed applications.** FLOWPanel already has `FGSPreconditioner` (`src/FLOWPanel_solver.jl:1873–1972`) and a prior benchmark. Its apply zeroes strengths, invokes a complete fixed-iteration FGS solve, and saves/restores body state. The existing phase-1 script selected sweep counts from a convergence ladder; that is not evidence that one inexpensive sweep will suffice. Include preconditioner initialization, farfield work, the outer operator, and orthogonalization. A lower Krylov iteration count alone can conceal more total FGS work. Use the existing path as a correctness starting point, then determine whether a cheaper internal apply is warranted.
3. **Mixed-precision correction against the original operator.** If Float32 stalls the true-accuracy gate, use it as an approximate solver/preconditioner and compute residual corrections with the original operator. This can preserve final accuracy while changing the iteration. Count the accurate residual passes and any retained Float64 storage; refinement is not a free certification layer.
4. **Fewer, coupling-aware GS/Jacobi domains or a real symmetric sweep.** Existing 64-way chunking lost despite faster sweeps. Fewer domains could recover convergence; geometric/graph partitions should minimize strong cut couplings rather than assume contiguous equal leaf counts are best. A true backward pass is a separate algorithm experiment. Both must beat the accepted wall-time target including changed sweep counts.
5. **Compress weak direct blocks.** Low-rank storage could cut bytes beyond Float32, but the small source dimension limits rank savings. For an $m\times n$ block, factors only save coefficient storage if $r(m+n)<mn$; demand a substantial margin for two products and metadata. First try the existing MAC/P controls. Compression introduces operator error and needs a rank/accuracy study.

Do not prioritize an extra stale farfield iteration: the measured FMM cost is too small to meet the target alone, and the iteration penalty is unknown. A GPU port could offer a larger hardware bandwidth budget, but the serial leaf sequence would require persistent device execution and resident solver state to avoid thousands of launches/transfers. It is a separate platform project, not an evidence-backed near-term recommendation from these CPU results.

For perspective, a hypothetical reduction from 27 to 15 updates at the same optimized per-update cost could take a 4.5–5.9 s solve toward roughly 2.5–3.3 s **before added acceleration costs and fixed work**. That is the route to roughly 3–4×, but it is conditional on convergence evidence and must not be sold as an achieved or expected result.

## Decision gates for implementation

1. **Correctness model:** compare complete per-sweep strengths, products, RHS, and residual against lex on small fixtures, including nonzero starts, multiple systems, branch-spanning targets, and supported rigid transforms. Test the source schedule in Float64 before changing precision. Test the split recurrence separately if pursued.
2. **Graph and real-shape performance gate:** audit the actual lower graph early, alongside source-kernel work. Replay its edges and measured variable costs, including backward release/drain, reduction, and NUMA policy. For source-major, replay the actual matrix-size distribution and source order, including row-to-target mapping, old/new bookkeeping, leaf handshakes, and small-block policy. Compare complete schedules at equal precision; record useful bandwidth, latency, allocations, and placement. Do not select split from isolated pull bandwidth. Reject candidates whose measured budget cannot reach 6.744 s with margin.
3. **Numerical gate:** certify the Float32 variant independently, recalibrating its tolerance where needed. Keep authoritative BC relative L2 ≤1e-6, evaluator certification/disagreement checks, memory limits, and repeatability requirements from the [cold harness](../../benchmark/fgs_cold_README.md). Do not apply within-configuration repeatability thresholds as an unqualified cross-algorithm equality test.
4. **End-to-end gate:** interleave uninstrumented repeats of the unchanged champion and candidate on the same node. Report distributions, completed updates, actual sweeps/FMM calls, setup cost, retained/peak memory, and accepted accuracy. Require ≥1.5× accepted solve throughput and aim for ≥2×. A faster isolated kernel does not satisfy this requirement.
5. **Continue beyond the first win:** use the new profile to choose between reduced scatter/residual cost, split layout, two sockets, joint parameter retuning, or acceleration. If bandwidth already saturates, prioritize less work over more scheduling machinery.

Relevant implementation checks include the solver, FMM, and FGS-history suites; add simulation/warm-start coverage for changed persistent state. Any local checks stay at ≤4 threads. Official timing campaigns use the required clean tagged worktrees and pinned dependencies. This document proposes those experiments; it does not initiate them.

The updated recommendation retains source-row execution as the provisional first implementation, but promotes pull-DAG analysis to an early competing experiment. Both preserve the required mathematical dependencies under the stated state-ownership rules. The remaining choice depends on measured handoff costs versus directed critical-path and backward-drain costs; the current evidence does not settle it. Accepted end-to-end accuracy and time remain the deciding criteria.
