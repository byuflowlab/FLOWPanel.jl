# 033 B-R1 + B-R2 — Existing-machinery audit and threaded FGS setup prototype (2026-09-26)

Companion to [`033_bt1_setup_profile_20260926.md`](033_bt1_setup_profile_20260926.md)
(B-T1: `nonself_influence_matrices` = 95.2% of R4 FGS setup, single-threaded).
Evidence: `BRAINSTORM/033_br2_20260926/` (CSVs, banners, logs). Same protocol
and host as B-T1 (MacBook, `OPENBLAS/OMP=1`, `BENCH_BLAS_THREADS=1`,
`julia --project=. -t {4,1}`, banner via `common.jl`, judged from CSVs;
FLOWPanel `c866f6f`-dirty, FastMultipole `745af760`-dirty).

## B-R1(i) — Hook-adoption gap analysis

B-T1 already established the first half: the FGS setup path
(`FastMultipole/src/solve.jl` `nonself_influence_matrices` /
`self_influence_matrices`) predates BRAINSTORM 030 and does **not** use the
`assemble_influence_block!` machinery. The other half — what adopting the
hook would actually compute vs what the probe computes — traced from source:

| aspect | FGS unit-strength probe (`solve.jl`) | 030 analytic hook (`nearfield_cache.jl:270`, FLOWPanel overload `src/FLOWPanel_abstractbody.jl:1354`) |
|---|---|---|
| block rows | ONE `influence!`-projected scalar per target body — FLOWPanel: `dot(U_induced, n̂_target)` (`FLOWPanel_abstractbody.jl:1427`), normals read from the *target system's* source-format buffer | `n_out` RAW kernel-output rows per target: `vec(target_buffer[output_range(switch), target_range])` (φ + U for the PS,GS switch) |
| block columns | ONE column per source body, strength set by `value_to_strength!` — champion `RigidWakeBody{Union{ConstantSource,VortexRing},2}` sets (σ, μ) = (0, 1) (`FLOWPanel_liftingbody.jl:1022`), i.e. **doublet-only** | `strength_dims = 2` columns per body (unit σ AND unit μ) |
| kernel work | one `induced` eval per (target, source-body) pair | `strength_dims`× that — for the NK=2 champion the analytic path evaluates the unit-σ column **only to have it discarded** by the (0,1) pattern: ~2× the kernel evaluations FGS needs |
| equivalence | — | FGS block = row-projection (`influence!`) × hook block × column-combination (`value_to_strength!` weights). Exact under the same linearity assumption FGS already makes, but **not bitwise** (different accumulation/composition order), and needs new composition code |

**Conclusion:** the 030 hook is *reusable in principle* but not
shape-compatible — adoption requires a projection/composition layer, breaks
bitwise equivalence, and (for the production NK=2 body) roughly doubles the
kernel-evaluation count. Meanwhile B-T1 showed the setup gap at HPC j64 is
overwhelmingly a **threading** gap (serial loop on 1 of 64 threads), not a
per-column-overhead gap. The cheapest credible fix is therefore to
**parallelize the existing probe**, not to reroute it through the hook (see
design fork below).

## B-R1(ii) — Unmerged branch audit

| branch | tip | ahead of production | contents | verdict |
|---|---|---|---|---|
| `influence_matrices` | `7dbc1a76` 2024-06-20 "add influence matrices, singlethreaded" | 0 commits (full **ancestor** of `flowpanel-20260817`) | the origin of the current serial probe itself | **DEAD/MERGED** — nothing to harvest |
| `faster-influence-matrices` | `5f95d08f` 2024-07-03 "sort_direct_list in progress" | 1 commit, on a base 665 commits behind production | a 17-line **unfinished stub** `sort_direct_list` in `src/fmm.jl`: empty loop body, no return, never called | **DEAD** — production already has working `sort_by_source`/`sort_by_target` counting sorts; superseded, not rebaseable-useful |

B-G gate: candidate *external* machinery is dead, but the dominant component
is 95% of setup and the fix below cost <200 lines — Track B proceeds
(trivially passes the gate on the "component ≥30%" side).

## B-R2 — Implementation: threaded influence-matrix population (toggle, default off)

Design: parallelize the **existing probe** over disjoint per-source-leaf
matrices. Each worker pulls whole matrices from an atomic counter and probes
on **private copies** of both buffer sets (`reset!`/`direct!` mutate target
buffers; unit strengths are set on the shared source buffers *before* the
copies so workers inherit them). Every column is produced by the identical
`reset!` → `direct!(single column)` → `influence!` → copy sequence as the
serial loop, and matrix writes are disjoint (one matrix per source leaf), so
the result is **bitwise identical at any worker count**.

Files (all uncommitted, on top of the existing dirty state):

- `FastMultipole/src/solve.jl`
  - `nonself_influence_matrices(...; setup_threads::Integer=0)` (line 143);
    `setup_threads=0` (default) runs the legacy serial populate **unchanged**;
    ≥1 branches at line 248 into the new path
  - `_populate_nonself_threaded!` (line 375): source-leaf grouping of the
    source-sorted list, save/set/restore strengths, `@spawn` worker pool with
    atomic counter + per-worker buffer copies
  - `_populate_nonself_group!` (line 429): faithful replay of the serial
    per-entry logic (incl. the multi-system row-cursor semantics) for one
    matrix
  - `self_influence_matrices(...; setup_threads::Integer=0)` (line 577);
    threaded branch at 597; `_populate_self_threaded!` (679),
    `_populate_self_leaf!` (710)
  - `FastGaussSeidel(...; threaded_setup::Bool=false)` (line 854): passes
    `setup_threads = threaded_setup ? Threads.nthreads() : 0` (line 912) to
    both matrix builders
- `FLOWPanel/src/FLOWPanel_solver.jl`
  - `FGSSolver(...; threaded_setup::Bool=false)` (line 1571), forwarded to
    `FastGaussSeidel` only when `true` (line 1608 — checkout-compat elision,
    same pattern as the dagteam kwargs)
- `FLOWPanel/benchmark/fgs_setup_ab.jl` (new): paired A/B harness — both
  arms in one process on the same trees/direct list, element-wise matrix
  certification, optional full-ctor + cold-solve certification

## Equivalence certification

Element-wise comparison of the two arms' assembled matrices in the same
process (`sorted_list`, `sizes`, `data`, `rhs` for both nonself and self):
**bitwise identical** (`max|Δ| = 0.0` exactly) at R1 j4, R1 j1, and R4 j4
(R4 j1 ran the arms in separate processes for the time budget — the path is
thread-count-independent and was certified bitwise at both other
configurations).

Solve certification (R1, `SKIP_B=0`, both j1 and j4): full `pnl.FGSSolver`
built both ways (champion knobs, `:dagteam`/`:f32full`), solver-internal
`nonself_matrices.data`/`self_matrices.data` bitwise identical; one cold
Dirichlet solve each: **niter 15 = 15, solved true = true, solution strength
column bitwise identical (rel Δ = 0.0)**. Since every downstream setup
artifact (leaf LU cache, dagteam plan, F32 repack) is a deterministic
function of these matrices, bitwise-equal matrices ⇒ an identical solver.

## Paired setup A/B (min over timed passes, warmup excluded; 3 timed passes at R1, 1 at R4 — R4 warmup passes agree with the timed pass within 0.7%)

| rung | threads | phase | old [s] | new [s] | speedup |
|---|---|---|---:|---:|---:|
| R4 | j4 | nonself | 135.66 | 37.84 | **3.59×** |
| R4 | j4 | self | 5.06 | 1.55 | 3.26× |
| R4 | j4 | nonself+self | 140.72 | 39.39 | **3.57×** |
| R4 | j1 | nonself | 137.02 | 136.71 | 1.00× |
| R4 | j1 | self | 5.06 | 5.14 | 0.98× |
| R1 | j4 | nonself | 7.34 | 2.08 | 3.53× |
| R1 | j4 | self | 0.48 | 0.14 | 3.4× |
| R1 | j1 | nonself | 7.40 | 7.39 | 1.00× |
| R1 | j1 | self | 0.48 | 0.50 | 0.96× |

j4 efficiency ≈ 0.9 (3.59/4); j1 shows the new path costs **nothing** when
run with one worker (buffer copies are noise). Full-ctor cross-check at R1
j4: 10.0 s → 2.31 s end-to-end (the ctor rows are compile-order-confounded
— old runs first and pays residual compile — so the phase rows are the
evidence of record; at R4 the two phases are 98.8% of the ctor per B-T1).

## Projected R4 HPC setup (proportional argument — CAVEATED)

Warm-start R4 campaign (zen3, j64): FGS setup 315–320 s, of which ~98.8%
(~300–316 s) is the serial nonself+self probe. Local j4 parallel efficiency
is ~0.9; assuming it degrades to 0.4–0.7 at j64 (memory-bandwidth-bound
kernel, shared-buffer copies ×64, per-leaf load imbalance absorbed by the
atomic pool over ~10³ matrices), the probe phase projects to

$$ t_{probe} \approx \frac{300\ \mathrm{s}}{64 \times (0.4\ldots0.7)} \approx 7\ldots12\ \mathrm{s}, $$

giving total FGS setup ≈ **10–20 s** (plus the ~4 s of other phases,
HPC-scaled) vs today's 315–320 s — i.e. the ~210 s FGS-vs-ILU setup gap
would not just close but invert (ILU setup+prime is ~108–110 s). Even a
pessimistic 25% j64 efficiency lands at ~19 s probe / ~25 s setup. This is a
**proportional projection, not a measurement** — the actual j64 number
belongs to B-I1's HPC re-measure.

## Design forks chosen

1. **Threaded probe over 030-hook adoption.** Rationale in B-R1(i): the hook
   needs a projection/composition layer (semantic gap), loses bitwise
   equivalence, and does ~2× the kernel work for the NK=2 production body;
   the profiled gap is a threading gap. The analytic hook remains available
   as a *further* serial-efficiency lever (no per-column `reset!`/dispatch,
   block-level `Val` regularization barrier) if B-I1's j64 measurement shows
   thread scaling saturating — bounded below by ~2 s (137 s / 64) either way.
2. **Kwarg toggle (`threaded_setup`), not env var**, default `false`; the
   legacy serial code path is byte-for-byte the same lines, taken by default.
3. **Whole-matrix work units** (one source leaf per pull) rather than finer
   block/column chunks — simplest disjointness argument; the atomic pool
   over hundreds-to-thousands of matrices absorbs leaf-size imbalance.

## Caveats

- Local MacBook screening at ≤4 threads; j64 scaling and memory-bandwidth
  behavior are projections (see above). Per-worker buffer copies scale as
  workers × (target+source buffer bytes) — at R4 ≈ tens of MB per worker,
  ~1–2 GB transient at j64; fine for zen3 nodes but worth a glance at R7.
- R4 j1 arms ran in separate processes (10-min local command budget); their
  bitwise certification is inherited from R1 j1 / R1 j4 / R4 j4 in-process
  checks plus thread-count-independence of the per-column computation.
- Champion R4 knobs transplanted onto R1 (same convention as B-T1).
- Both checkouts dirty (pre-existing A-R2 + MATRIX_OPERATOR_REFACTOR
  uncommitted state; nothing reverted) — base-revision provenance, not a
  campaign pin. Nothing committed, per task boundary.
- The full FLOWPanel test suite was not run; the default path is untouched
  (kwarg default off, serial lines unchanged) and the R1 solve
  certification exercises the new path end-to-end.
- `threaded_setup=true` at `-t 1` still takes the new path with one worker
  (measured parity) — it does not silently fall back to the legacy loop.
