# Thread-efficiency ideation: context-reset handoff (2026-09-18)

Follows the NUMA-placement investigation (commit `8e3e01c`; findings in
`numa_placement_findings_20260918.md` — read it first, it is short and
carries the ceiling numbers everything below cites).

**Task (Ryan 2026-09-18):** propose anything else that could speed up the
FGS solve by using threads more efficiently. Examine the code/algorithm
directly. **Keep the algorithm intact first** — exhaust
algorithm-preserving improvements before proposing algorithmic changes;
algorithmic changes ARE authorized as a second tier, and larger changes
are acceptable after smaller ones are covered if you judge them more
impactful. This is an ANALYSIS/PROPOSAL task: deliverable is a ranked
proposals doc, not code. No job submissions are pre-approved; local
experiments are fine at ≤4 threads (Ryan's global rule); any repo code
change or HPC campaign is Ryan-gated.

## State of play (don't re-derive)

Rotor-hover R4 FGS benchmark, m12 node (dual-socket zen3, NPS4 → 8 NUMA
nodes × 16 cpus; jobs pin the 64 physical cores of socket 0 = NUMA nodes
0–3, first-touch default policy).

- **Operating points (medians, time-to-accepted-accuracy):**
  colored@j16 **10.116 s = global best**; lex@j64 10.95 s; lex@j16
  12.48 s; chunked@j64 14.26 s (loses everywhere — iteration inflation
  27→44 from majority-Jacobi outweighs a real span win).
- **Per-iteration chain (v20 counters, lex, j64):** nonself dgemv
  products 0.292 s + cross-leaf scatter 0.032 s + leaf-LU solves 0.020 s
  = 0.346 s/it, ~85% of j64 wall. Chunked shrank the span to 0.277 s/it
  at 38.8 avg active threads — the parallelization mechanism works, the
  memory doesn't keep up.
- **It is bandwidth-bound:** each inner sweep streams the full 2.86 GB
  nonself influence cache (1,068 leaf-pair block matrices), ×3 inner
  sweeps/iter ≈ 8.6 GB/iter. One core sustains 29.4 GB/s (measured =
  predicted). NUMA microbenchmark (job 13763831): serial first-touch
  piles 93% of pages on ONE NUMA node → 64 threads get only 74 GB/s;
  interleave 154; chunk-affine first-touch 164 GB/s ≈ the practical
  socket ceiling. So **placement fixes give ~2.2× on the products, no
  more; a 64-thread socket cannot exceed ~164 GB/s ÷ 29 GB/s ≈ 5.6× over
  one core on this kernel no matter how threads are scheduled.** Real
  speedups beyond that must come from moving FEWER bytes (Float32 —
  already Ryan-approved, halves bytes; fewer sweeps/iterations; cache
  blocking) or from overlapping/eliminating serial+barrier time.
- **Known thread-efficiency sinks already measured:** barrier spin
  (busy-time 10.7 CPU-s per 0.28 s chunked span = starved stragglers
  spinning); colored@j64 loses to colored@j16 via 79 colors × 81 sweeps
  of barriers with median color width 16; chunked's cross-chunk scatter
  is serial by design; leaf-LU solves are small and serial-ish.

## Code map (FastMultipole, branch `flowpanel-20260817` @ `c18e4b46`, local checkout `/Users/ryan/Dropbox/research/projects/FastMultipole`)

- `src/solve.jl` — the whole FGS: `_calloc_vector` :6, `Matrices` ctor
  :13; sequential `nonself_influence_matrices` fill :256-321 (the
  single-thread first-touch culprit); `self_influence_matrices` ~:427;
  `FastGaussSeidel` ctor ~:640-700; the iteration/sweep loops, coloring
  (`color_leaves`), chunk machinery (`build_chunk_map`, scatter
  partition) further down — grep `sweep_order`, `nearfield_update`.
- `src/containers.jl:1092-1150` — `Matrices`, `LeafLUCache`,
  `FastGaussSeidel` struct (incl. `sweep_order`, `chunk_ranges`,
  `scatter_intra/cross`).
- FLOWPanel side: `src/FLOWPanel_solver.jl` (`FGSSolver` plumbing);
  harness `benchmark/fgs_r4_chunked_ab.jl`.
- The parallel `NearfieldInfluenceCache` build (`src/nearfield_cache.jl`
  :455-478, atomic chunk pool) exists but is NOT the FGS path — its
  pattern may be reusable.

## Seed ideas already on the table (evaluate, don't duplicate blindly)

Ryan-gated but staged: **Float32 nearfield storage** (approved; ~2× bytes)
and **chunk-affine first-touch assembly** (arm d beat interleave 164 vs
154 GB/s, needs no numactl). Follow-ups mentioned but undecided: fewer
chunks / under-relaxation for chunked's iteration inflation,
parallel-by-target deferred scatter, per-chunk timers / uncore counters,
Krylov-accelerated outer loop (parked). Your job is to go beyond these:
e.g. think about sweep/iteration count (why 3 inner sweeps? convergence
vs streamed bytes trade), overlap of farfield/M2L work with nearfield
streaming, using BOTH sockets (membind 0–7 is available; cpubind was 0–3
by choice), leaf-block sizing vs cache/TLB, batching dgemv into dgemm
across leaves sharing a source, software prefetch/nontemporal patterns,
barrier-free scheduling (work-stealing, dependency-DAG sweeps), or
anything else the code supports — but rank honestly against the 164 GB/s
socket ceiling: an idea that doesn't cut bytes, cut iterations, or add
memory controllers cannot beat 5.6×/socket.

## Ground rules / gotchas

- Local runs ≤4 threads; don't submit anything; don't edit repo code —
  scratch experiments go in the scratchpad or `benchmark/` proposals only.
- Evidence to cite: `numa_placement_findings_20260918.md`,
  `fgs_r4_followup_evidence_20260914/chunked-v22-13749231/analysis/ab_summary.md`,
  `fgs_opt_r4_diagnostics_package_20260912.md` §7, item decision log.
- Pre-existing dirty files from other campaigns (018/026, data/,
  examples/) — leave them alone. Known unrelated test failure:
  WeakKeyDict × `WarmstartNoopSolver` (`src/FLOWPanel_solver.jl:2331`) —
  out of scope.
- **Wrap:** proposals doc `thread_efficiency_proposals_<date>.md` in this
  directory, ranked (impact estimate vs the 10.116 s operating point,
  effort, algorithm-preserving vs algorithm-changing, risk), plus a
  one-paragraph recommendation of the top 1–2. Update the 021 item
  Current status + log.md + memory, commit. All implementation and any
  notebook entry (4 already owed) remain Ryan-gated — OFFER, never write.
