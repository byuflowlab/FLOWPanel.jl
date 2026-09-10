# 030 — Generic influence-block assembly (`assemble_influence_block!`)

**Opened:** 2026-09-09 (Ryan directive). **Status:** staged; design frozen, implementation delegated.

**Item-level approvals:** Technical [ ]; clear-context [ ]; user [ ]

## RESET BRIEF

**Goal:** let any FastMultipole system expose a per-pair influence primitive so dense
influence blocks (near-field cache, FGS solver matrices) are ASSEMBLED directly instead
of recovered by unit-strength probing through `direct!` — while `direct!` remains the
only REQUIRED kernel interface (probing stays as the universal fallback).

**Motivation (measured 2026-09-09, local M-series, 4 threads):** for a cheap kernel
(gravitational), the kernel arithmetic is only ~17–18% of the serial near-field cache
build; the other ~82% is probe bookkeeping — per column: zero the target output rows,
one `direct!` call per (block × source column) (~leaf-size targets per call), copy the
column into the packed matrix (5 memory touches/entry), plus a full repeat of the
geometry per strength component when `sd > 1`. Threaded scaling of the probe caps at
~2.5–3× on 4 threads because the residual work is memory ops on shared cache. For the
FLOWPanel panel kernel (~µs/pair) the kernel is >99% of the build and none of this
matters — this item is for cheap kernels and for making the FGS block builders cheaper.

**The pattern already exists in three places (reconnaissance 2026-09-09):**

| builder | pattern | stores |
|---|---|---|
| FLOWPanel `_G!` (`src/FLOWPanel_solver.jl:237`) | direct per-pair assembly via `induced` | solver projection (φ or u·n̂) |
| FGS `self_influence_matrices` (`FastMultipole/src/solve.jl:408`), `nonself_influence_matrices` (`:124`) | unit-strength probe via `direct!` | projection via `influence!` |
| near-field cache `_build_nearfield_cache` (`FastMultipole/src/nearfield_cache.jl`) | unit-strength probe via `direct!` | raw buffer output rows |

FLOWPanel `_G!` proves the target pattern: unit-activate all strengths once, threaded
double loop over (source, target) calling `induced` per pair, write `G[i,j]` directly.
This item promotes that pattern into a generic FastMultipole interface.

## Design (frozen)

New optional interface function in FastMultipole:

```julia
assemble_influence_block!(block, target_buffer, target_range, switch,
                          source_system, source_buffer, source_range)
```

- `block` is the m×n dense block (m = n_out·|targets|, n = sd·|sources|), column layout
  matching the probe's: column j = (i_body − first(source_range))·sd + i_comp; rows =
  vec of (out_range × target_range).
- **Default method = the existing probe loop** (unit columns through `direct!`), so
  systems that only define `direct!` are unaffected. The near-field cache builder calls
  the hook per block; its worker-pool parallel structure, packed `Matrices` storage,
  size/time guards, and bit-identical-at-any-thread-count guarantee are unchanged
  (each block is one work unit either way). Opted-in systems need no per-worker buffer
  copies (they never write the target buffer) — keep the copies for the fallback path.
- **Semantics:** entries are ASSIGNED (not accumulated); the hook must write every
  entry of `block`. The kernel contract it asserts is linearity in the declared
  strength rows — the same contract the probe certifies implicitly.
- **Mandatory exactness test per opt-in:** assembled block vs probed block, rtol 1e-12
  (bitwise NOT required — accumulation order may differ), plus a mixed build (one
  opted-in system, one fallback system) in the same cache.
- **Projection split (decision):** the cache stores raw output rows → the hook returns
  raw rows. The FGS builders store post-`influence!` scalars → they keep probing in
  Phase 1–2 and migrate in Phase 3 by projecting the assembled raw block (or a
  projected hook variant — decide at Phase 3 with measurements).
- calloc-backed `Matrices` (zeros instead of undef): pairs with assignment semantics,
  free zeroing via OS zero pages, first-touch NUMA placement by the owning worker.

**Depends on (uncommitted, session 2026-09-09):** threaded probe build (worker pool +
atomic counter, private buffer copies per worker), `n_threads` plumbing,
`NearfieldCacheDonor`/`retarget_nearfield_cache`, FLOWPanel donor plumbing
(`KrylovSolver → KrylovOperator → _apply_*_G! → influence!`), tune-driver `NF_DONOR`.
All tests green: FastMultipole cache test file (25+8+17+4 pre-existing, 6 parallel,
8 retarget), FLOWPanel `runtests_unit_solver.jl` 465/465. Commit these first.

## Phases

- **Phase 1 (FastMultipole):** `assemble_influence_block!` hook + probe default;
  `_build_nearfield_cache` calls it per block; calloc `Matrices`; opt-in for the
  gravitational test system; exactness + determinism + mixed-system tests; measure
  speedup on the cheap-kernel benchmark (expect ~3–5× serial over the probe, better
  thread scaling).
- **Phase 2 (FLOWPanel):** opt-in overload for `AbstractBody` (both Neumann/Dirichlet
  raw-output forms) built on `induced` — the same per-pair primitive `_direct_body!`
  and `_G!` already use, honoring the `Val(FILAMENT_REGULARIZATION[])` function
  barrier once per block; exactness tests vs probed blocks on bodies WITH shedding
  panels (attached-wake term exercised, cf. the existing cache_nearfield testset);
  full `runtests_unit_solver.jl` green.
- **Phase 3 (Ryan-gated):** migrate FGS `self_influence_matrices` /
  `nonself_influence_matrices` to the hook (projection question above); measure.
- **Phase 4 (staged, not planned):** route block assembly through the 051 GPU
  rectangular-influence seam (`FLOWPANEL_GPU_INFLUENCE`) — the hook is exactly the
  interface a GPU path plugs into.

## Log

- 2026-09-09: opened; design frozen from session reconnaissance (probe cost anatomy
  measured, three existing builders identified, `_G!` pattern chosen for promotion).
  Implementation prompt below.

## Implementation prompt (Phase 1 + 2)

See `BRAINSTORM/030_implementation_prompt_20260909.md`.
