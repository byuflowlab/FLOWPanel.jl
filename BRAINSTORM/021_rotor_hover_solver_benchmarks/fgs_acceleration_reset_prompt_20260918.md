# Reset prompt: implement the FGS acceleration recommendation (2026-09-18)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Task

Implement the program specified in
`BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_acceleration_recommendation_20260918.md`
(BRAINSTORM item 021). That document is the authoritative spec — read it in
full before writing any code, along with its companion
`thread_efficiency_top5_20260918.md` (compact execution ordering; where the
two disagree, the recommendation wins). Background proposals, only if needed:
`thread_efficiency_proposals_20260918.md`, `fgs_split_dual_layout_ideas_20260918.md`.

Also read `CLAUDE.md` and `agent_policies/WORKFLOW.md` + `agent_policies/TESTING.md`
before editing, and `agent_policies/HPC.md` before any cluster work.

## Goal and acceptance

Accelerate the prepared cold FGS body solve at R4 (rotor-hover benchmark,
retained config P8 / MAC0.4 / leaf100 / inner3 / lex / LU cached — see
`benchmark/retained_r4_diagnostics.toml`). Baseline = best accepted R4
result: **colored @ j16 = 10.116 s** (26 iterations).

- Minimum acceptance: **≤ 6.744 s** (1.5× throughput).
- Design target: **≤ 5.058 s** (2×).
- Planning range for the recommended implementation: 4.5–5.9 s (projection,
  not a promise).
- Accuracy gate unchanged: authoritative BC relative L2 ≤ 1e-6 with the
  independent evaluator; see `benchmark/fgs_cold_README.md` for the cold
  harness, certification, and repeatability requirements.

## What to build (from the recommendation, in order)

1. **Persistent, adaptive row-parallel execution of the existing source-major
   nonself cache**, with consumer-aligned memory placement. Keep the
   lexicographic leaf sequence, far-field refresh schedule, and current RHS
   semantics (`rhs += old` then `rhs -= new`, two separate ops — never a
   delta product). Coordinator solves each leaf's diagonal block and
   publishes strengths; a persistent worker team computes disjoint row tiles
   (contiguous rows *within each column*; never one strided dot per row) of
   that source's tall matrix, accumulating in private scratch, then scatters
   to owned target rows. Adaptive team size by block bytes: tiny serial,
   medium compact team, large across controllers. No allocations or task
   creation inside the leaf loop.
2. **Float32 nonself coefficient storage with Float64 arithmetic** and
   Float64 self matrices / LU / strengths / products / RHS / residual. This
   is an explicit container/dispatch change — `Matrices{TF}`
   (FastMultipole `src/containers.jl:1092`) couples coefficient and
   product/RHS types, and `FastGaussSeidel` couples self and nonself
   precision. Convert tiles on load, accumulate in Float64. Do NOT round
   strengths to Float32 to unlock `sgemv`. Selective per-block Float64
   retention is the fallback if precision fails (traffic 1−f/2), never a
   relaxed gate.
3. **Early pull-DAG audit in parallel with 1** (analysis + replay scheduler,
   not a solver rewrite): the split triangular design is the competing
   candidate. Structural facts already recomputed and verified (see the
   recommendation's table): 48,167 directed lower edges; longest lower path
   279 of 1,068 leaves; byte-weighted work/span **2.849**; largest
   earliest-level cohort 6. Replace structural weights with measured
   per-leaf pull/solve costs via the recurrence in the recommendation
   before deciding. Promote split only if a complete inclusive schedule
   comparison at equal precision/placement beats source-major with margin.

Deferred / out of scope unless the profile says otherwise: joint
leaf/MAC/P/inner retune (separate campaign, changes operator), safeguarded
Anderson / FGS-preconditioned FGMRES (algorithm track), second socket, huge
pages, low-rank compression, GPU.

## Decision gates (do not skip or reorder)

From the recommendation §"Decision gates", condensed:

1. **Correctness model** on small fixtures in Float64 FIRST: compare
   complete per-sweep strengths, products, RHS, residual against lex —
   including nonzero starts, multiple systems, branch-spanning targets,
   rigid transforms. Precision change comes after the schedule is proven.
2. **Real-shape sequence benchmark** (go/no-go): replay the actual R4
   block-size distribution and source order (census:
   `fgs_r4_followup_evidence_20260914/diag-v15-13694724/j64-b1/results/gemv_census.csv`,
   edges: `.../dependency_edges.csv` — filter `source_leaf < dependent_leaf`
   for lower edges; edge weight 8·n_j·n_i sums to 2,862,850,032 bytes/sweep).
   **Handoff budget is the kill switch**: 1,068 leaf handoffs/sweep × 81
   sweeps = 86,508; average must stay < ~5.8 µs to keep added coordination
   under 0.5 s. Reject candidates whose measured budget cannot reach 6.744 s.
   This same benchmark prices the pull-DAG promotion condition.
3. **Numerical gate**: independent Float32 certification, recalibrated
   tolerance if needed; the internal residual of the rounded operator can be
   small while the true BC residual fails — the external evaluator rules.
4. **End-to-end A/B**: interleaved uninstrumented repeats of unchanged
   champion vs candidate on the same node; report distributions, completed
   updates, sweeps/FMM calls, setup cost, memory, accepted accuracy.
   ≥1.5× accepted solve throughput required.

## Cost model (for planning, not reporting)

T ≈ 3.0 s remainder + coeff_GB/B + H, with 231.891 GB (F64) or 115.945 GB
(F32) per 81-sweep solve. Halving wall time needs 112.7 GB/s (F64) or
56.3 GB/s (F32; 74.4 with H=0.5 s). Measured ceilings: 29.4 GB/s one core,
164.4 GB/s one socket with affine placement (synthetic —
`benchmark/numa_dgemv_bench.jl`, findings in `numa_placement_findings_20260918.md`).
Never multiply the 2.22× placement ratio by a threading ratio — shared
ceiling.

## Code anchors

Implementation lives mostly in the sibling **FastMultipole** checkout
(`../FLOWVPM.jl` is unrelated; FastMultipole path is whatever the FLOWPanel
Manifest dev-path points to — verify, last seen at commit `c18e4b4`):

- `src/solve.jl:1208–1302` `gs_sweep!` (lex leaf loop: solve → product → scatter)
- `src/solve.jl:923` `compute_nonself_products!` (full new product, old saved)
- `src/solve.jl:948` `scatter_nonself_influence!` (`+= old`, `-= new`)
- `src/solve.jl:143`, fill `:256–321` `nonself_influence_matrices` (serial
  build, one tall matrix per source leaf; target branches may span multiple
  solve leaves — a direct-list entry ≠ one solve leaf)
- `src/solve.jl:1803` `residual!` — **all leaves share
  `view(residual_vector, 1:length(rhs))`; threading it races without
  worker-private scratch**
- `src/solve.jl:1126–1160` `color_leaves` — symmetrized adjacency, NOT the
  directed pull graph
- `src/solve.jl:1380–1537` outer loop (init from actual strengths, farfield
  refresh, final convergence FMM pass — FLOWPanel sets `final_update=false`
  at `src/FLOWPanel_solver.jl:1814–1824`, so don't "optimize away" a final
  direct evaluation that isn't there)
- `src/containers.jl:1092` `Matrices{TF}`
- `nearfield_cache.jl:438–478` — pattern for private-buffer parallel assembly
- FLOWPanel `src/FLOWPanel_solver.jl:1873–1972` `FGSPreconditioner`
  (reference only; algorithm track is out of scope)

## Known traps (each cost a prior agent time)

- Row-tiled custom kernels are mathematically equivalent to BLAS, not
  bit-identical (tiling/SIMD/FMA). Check explicitly; else certify accuracy.
- Source-affine first touch is WRONG placement for row-worker consumption;
  placement must match the selected schedule; verify with numastat, don't
  assume. Sub-page tiles / huge pages can defeat ownership.
- Chunked v22 LOST end-to-end (14.264 s, 44 iterations) despite faster
  sweeps; its "38.8 active threads" measures a Jacobi-lagged schedule plus
  spinning — not evidence of exact-GS frontier width.
- Colored and lex use separately calibrated tolerances / orderings — don't
  mix their iteration counts.
- The reverse flag currently repeats forward order — don't silently change it.
- For split (if promoted): backward products must not overwrite the frozen
  upper accumulator u^s while targets still read it; initialize Ux^0 for
  warm starts.

## House rules (binding)

- **Local runs ≤ 4 threads**, always. Full-scale timing only on HPC.
- HPC job **submission stays with Ryan-approved flow**; monitoring/log
  tailing goes through the `hpc-monitor` subagent; HPC policy in
  `agent_policies/HPC.md`. `ssh orc` needs a live ControlMaster socket.
- Official timing campaigns run from **clean tagged git worktrees**
  (`campaign/<item>-<slug>-YYYYMMDD` annotated tags) with Manifest dev-paths
  pointed at the campaign worktrees; pins recorded in a provenance file
  before submission. Both live checkouts (FLOWPanel `fastmultipole` branch
  and FastMultipole) currently carry uncommitted work from other threads —
  **do not commit, revert, or build on top of unrelated dirty state; do your
  development in your own worktree/branch** and ask Ryan how to reconcile if
  a needed file is dirty.
- Tests: pick via `agent_policies/TESTING.md`; relevant suites = solver,
  FMM, FGS-history; a pre-existing Kutta `:jump` failure is known (item 030).
- Notebook writes, campaign submissions, and anything outward-facing are
  Ryan-gated. Sequence-benchmark results, gate verdicts, and the
  source-major vs pull-DAG decision come back to Ryan before the end-to-end
  campaign.
- Log progress as dated status/provenance files in
  `BRAINSTORM/021_rotor_hover_solver_benchmarks/` (existing naming
  convention: `<topic>_status_YYYYMMDD.md`, `<topic>_provenance_YYYYMMDD.md`).

## Suggested first moves

1. Read the recommendation and top-5 note end to end.
2. Locate the FastMultipole checkout via FLOWPanel's Manifest; confirm the
   anchors above still match (HEADs may have moved since `a6270c8` /
   `c18e4b4`).
3. Build the small-fixture correctness harness (gate 1) and the census/edge
   replay benchmark (gate 2) — these are cheap, local, ≤4 threads, and
   decide everything downstream.
4. Report gate-2 numbers (handoff latency, useful bandwidth, DAG replay) to
   Ryan with a go/no-go before touching the production solver path.
