# 033 B-T1 — FGS setup-cost attribution (2026-09-26)

Local screening profile of `FGSSolver` construction at the R4 champion
operating point (P8 / MAC 0.4 / leaf 100 / `:dagteam` / `:f32full` /
`cache_leaf_lu=true`), decomposed into the actual phases of
`FastMultipole.FastGaussSeidel`. Harness: `benchmark/fgs_setup_profile.jl`
(new; replays the constructor's own call sequence non-invasively — no
library edits, no solver-behavior changes). Evidence:
`BRAINSTORM/033_bt1_20260926/` (CSVs, banners, logs).

**Protocol.** `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 THREADING_MODE=multi
BENCH_BLAS_THREADS=1`, `julia --project=. -t {4,1}`, `SKIP_B=1` (geometry-only
fixture; nothing is solved). Timing = `time_ns` around each phase, min over
k timed passes after 1 compile-warmup pass (R4 j4: k=2 decomposed + 2
end-to-end ctor; R4 j1: k=1+1, adequate for a phase this flat and
multi-minute; R1: k=3 + 3). Judged from the CSVs. Host = Ryan's MacBook
Pro (`mecsrs-MacBook-Pro-188.local`), FLOWPanel `c866f6f`-dirty
(`fastmultipole`), FastMultipole `745af760`-dirty (`flowpanel-20260817`,
uncommitted A-R2 coop-executor changes, default-off). Local screening only —
NOT HPC evidence; the warm-start R4 campaign numbers quoted below are from
the zen3 node at j64.

## Phase table — R4 (58,192 panels), 4 Julia threads

Min over 2 timed passes; alloc = GC bytes during the phase (pass 1).

| phase | what it is | t [s] | % of sum | alloc [GB] | threaded? |
|---|---|---:|---:|---:|---|
| trees | source tree + target-tree replay + topology assert | 0.026 | 0.02% | 0.08 | yes (`tree.jl` `@threads`) |
| lists | interaction lists + canonical double-sorts | 0.003 | 0.00% | 0.01 | no |
| **nonself** | **near-field non-self influence blocks (serial per-source-column unit-strength probe)** | **135.2** | **95.2%** | **4.93** | **NO** |
| bookkeep | self-interaction list, index maps | 0.001 | 0.00% | 0.00 | no |
| self | leaf self-influence blocks (same serial probe) | 5.01 | 3.5% | 0.18 | NO |
| lu | per-leaf Float64 `lu!` cache | 0.18 | 0.13% | 0.11 | NO (serial `map`, BLAS=1) |
| dagplan | dagteam DAG + Lmat/Umat split repack (F64→F32) + F32 leaf LU + scratch | 1.48 | 1.0% | 1.54 | NO |
| phase_sum | | 141.9 | 100% | 6.8 | |
| full ctor (cross-check) | real `pnl.FGSSolver` end-to-end | 142.4 | — | 7.3 | |

Phase sum accounts for 99.7% of the real constructor (min-of-2, per
protocol; residual ≈ 0.5 s: `calc_normals!/controlpoints!`, solver-struct
allocs, replay-loop overhead).
Inside `dagplan`, the standalone Float32 leaf-LU refactorization
(`build_leaf_lu_cache_as`) re-times at 0.12 s — the repack loops dominate
the 1.48 s.

## Phase table — R1 (8,016 panels), 4 threads

Min over 3 timed passes:

| phase | t [s] | % |
|---|---:|---:|
| nonself | 7.35 | 93.0% |
| self | 0.47 | 6.0% |
| dagplan | 0.038 | 0.5% |
| lu | 0.009 | 0.1% |
| trees + lists + bookkeep | 0.004 | 0.05% |
| phase_sum | 7.90 | 100% |
| full ctor | 8.03 | — |

**Rung scaling (R1→R4, 7.26× panels):** nonself 18.4× (≈ n^1.47 at these two
points — near-field block area grows faster than panel count at fixed
leaf/MAC), self 10.6×, dagplan 39×, lu 21×. For Tier-2 arithmetic: setup is
one phase to first order, and it grows superlinearly with rung.

## Thread scaling — 1 vs 4 threads

Min-of-k per phase (R4 j1: k=1 timed pass after warmup; others as above):

| phase | R4 j1 [s] | R4 j4 [s] | j4/j1 speedup | R1 j1 [s] | R1 j4 [s] |
|---|---:|---:|---:|---:|---:|
| nonself | 136.14 | 135.16 | 1.01× | 7.45 | 7.35 |
| self | 5.03 | 5.01 | 1.00× | 0.48 | 0.47 |
| dagplan | 1.62 | 1.48 | 1.09× | 0.18 | 0.038 |
| lu | 0.181 | 0.182 | 1.00× | 0.009 | 0.009 |
| trees | 0.037 | 0.026 | 1.44× | 0.003 | 0.004 |
| phase_sum | 143.0 | 141.9 | 1.01× | 8.15 | 7.90 |
| full ctor | 143.1 | 142.4 | 1.01× | 8.12 | 8.03 |

**Measured: setup is flat in thread count (1.01× total at R4, 4 vs 1
threads)** — every material phase is a serial loop, exactly as the source
reads. Only `trees` shows real scaling and it is 0.02% of setup. (The R1
`dagplan` j1/j4 spread is systematic across passes but sub-0.2 s absolute —
consistent with the ~0.14 s j1 excess R4 also shows — and immaterial.) Track B's lever is therefore wide open: at HPC j64 the dominant
phase currently uses 1 of 64 threads.

## Dominant phase and its mechanism

`nonself_influence_matrices` (`FastMultipole/src/solve.jl:143`) is ~95% of
FGS setup and is **single-threaded**. Mechanism (from source): for every
direct-list block it probes one source body at a time — `reset!` the target
block, `direct!` with a unit strength for that single source column, then
`influence!` and copy the column into the matrix. So each matrix column pays
a full per-column kernel dispatch + buffer-reset + influence-projection
round trip, serially over all near-field columns.
`self_influence_matrices` (`solve.jl:439`, the remaining ~3.5%) is the same
probe pattern over the leaf diagonal blocks.

By contrast, the **Krylov/ILU near-field cache** builds the *same class of
near-field blocks* via `FastMultipole/src/nearfield_cache.jl`: parallel
workers over key chunks (`Threads.@spawn`, atomic chunk queue) calling the
`assemble_influence_block!` hook — which FLOWPanel **opts into analytically**
(`src/FLOWPanel_abstractbody.jl:1354`), writing block entries directly with
no unit-strength probing. That is the BRAINSTORM 030 machinery (merged in
FastMultipole `ac7230a6`); the FGS setup path in `solve.jl` predates it and
**does not use it** (this also answers half of B-R1(i)).

## Gap attribution vs krylov_ilu_nfcache (~210 s at R4 j64)

Warm-start R4 campaign (fgs_warmstart_r4_results_20260925.md, zen3 j64):
FGS setup 315–320 s; ILU setup+prime 108–110 s (nfcache build included) —
gap ≈ 210 s. Applying this profile's proportions to the HPC number:
~300 s of the 315–320 s FGS setup is the serial near-field probe
(nonself+self ≈ 98.8% at R4 locally), while ILU assembles its near-field
blocks threaded + analytic inside its ~82 s setup. **The gap is therefore
attributable almost entirely to one phase: FGS's serial, probe-based
near-field assembly.** Tree/DAG planning, LU factorization, and F32
conversion are collectively ~1.5% and irrelevant to the gap.

Closing arithmetic: adopting the ILU-style assembly (threaded +
`assemble_influence_block!`) for `nonself_influence_matrices` /
`self_influence_matrices` bounds FGS setup at roughly the ILU cache-build
cost plus the ~4 s of other phases (HPC-scaled). Even the *analytic-serial*
intermediate (hook without threading) should cut the probe overhead
substantially; threading at j64 multiplies that. A residual FGS-only cost
that cannot be removed: the dagteam split repack + F32 LU (~1.5 s local,
second copy of the near-field coefficients), which is noise at this scale.

## Recommended Track B lever

**B-R2 candidate = route FGS near-field assembly through the existing
BRAINSTORM 030 block-assembly machinery** (or directly reuse
`nearfield_cache.jl`'s parallel chunked builder) for both
`nonself_influence_matrices` and `self_influence_matrices`:

1. per-block `assemble_influence_block!` (FLOWPanel already opts in,
   analytic, probe-free, bit-comparable per 030's acceptance), with the
   probe as fallback for systems without the overload;
2. parallel over blocks (writes are disjoint by construction — each block
   belongs to one (target-range, source-leaf) key).

This is exactly the B-R1(i) audit's hypothesis; the profile confirms it is
the only phase worth touching. B-G stop-gate check: dominant component is
≥30% of setup by a wide margin (95%), so Track B proceeds.

## Caveats

- Local MacBook screening at ≤4 threads, BLAS=1; HPC magnitudes differ
  (315–320 s vs 143 s here) but the phase *proportions* are structural
  (serial loops scale with the same near-field work everywhere).
- Champion TOML is R4-specific; R1 rows transplant the same knobs
  (P8/MAC0.4/leaf100/f32full), not an R1-tuned champion.
- Both checkouts dirty (uncommitted, default-off A-R2 changes in
  FastMultipole; harness additions here) — base-revision provenance, not a
  campaign pin. `SKIP_B=1`: nothing was solved; this is setup-cost evidence
  only.
- Two rungs only; the n^1.47 scaling exponent is a two-point estimate.
- The B-T1 checklist asks for "R4 (and one larger mesh)"; this local screen
  used R4 + the *smaller* R1 (a larger-than-R4 rung is not feasible on the
  MacBook at these serial-probe costs). The larger-mesh point, if still
  wanted, belongs to the HPC re-measure in B-I1.
- `gc_bytes` measures allocation traffic, not peak RSS.
