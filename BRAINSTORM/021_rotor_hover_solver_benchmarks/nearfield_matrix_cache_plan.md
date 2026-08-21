# Near-field influence-matrix cache — implementation plan (EXECUTING 2026-08-19, Ryan-directed)

**Status: IMPLEMENTED 2026-08-19 (commits 1–2 = FastMultipole 72d3f3d +
02f071c; commit 3 = FLOWPanel surface in the working tree, uncommitted per
the launch constraint; commit 4 `persistent_plan` DEFERRED pending Ryan's
go). V0 PASSED (no refusal path needed). R1 local smoke: 265× isolated
near-field, 5.5× cold krylov_gmres solve, break-even 8.3 applies — see
`phase_02_single_step_benchmarks.md` Log 2026-08-19 and `ledger.md`.
AMENDED same day per Ryan's feedback round (FM 0ef4e83 + 204188a): the
tune=true refusal (§6) is superseded — cached-path tuning is supported
(provided cache tunes expansion_order only; `tune_fmm
tune_nearfield_cache=true` uses per-trial throwaway caches with
stop-at-cap semantics), and the builder gains an estimate-first
`max_build_time` guard beside `max_bytes`.** Repos: FLOWPanel.jl (this
repo) + FastMultipole dev checkout `../FastMultipole` (branch
`flowpanel-20260817`, HEAD 752c525). All facts below were verified against
code 2026-08-19 — spot-check line numbers, don't re-explore. Priorities
(Ryan): performance > robustness > user-friendliness. Test-gated; separate
commits; nothing rsyncs to the cluster until the Phase 2b A/B.

## 1. Idea and payoff

For frozen relative geometry, every near-field direct interaction is a LINEAR
map from source strengths to target outputs. Cache each direct-list entry's
dense block once, then re-evaluate the whole near-field as packed BLAS
matvecs (~2 flops/entry) instead of analytic kernel calls (~10²–10³
flops/pair for `compute_source_dipole`). Measured motivation: R1 `krylov_ilu`
profile — ~40% of solve samples in `nearfield_multithread!` →
`induced`/`_induced`/`compute_source_dipole`/`rotate_to_panel`
(`benchmark/results/phase2/profile/multi/R1/`). Amortization: within one
Krylov solve (7–59 applies at R1), across FGS sweeps, and — the headline —
across TIMESTEPS whenever a `Tree` (and hence relative positions) persists.

## 2. Verified code facts the design builds on

- **Packed container to REUSE verbatim**: `Matrices{TF}`
  (`FastMultipole/src/containers.jl:279-285`): `data` (all blocks
  concatenated, contiguous), `rhs`, `sizes::Vector{Tuple{Int,Int}}`,
  `matrix_offsets`, `rhs_offsets`. Built by `Matrices(sizes, TF)`
  (`src/solve.jl:3-27`, prefix-sum offsets; `EmptyMatrices` :28-31).
  Accessors `get_matrix_range`/`get_rhs_range`/`get_matrix_vector`
  (`src/solve.jl:32-51`) — `get_matrix_vector` returns a zero-copy `reshape`
  of a contiguous view, so `mul!` hits BLAS. Trap: `get_matrix_range`
  ignores its passed `(m,n)` and re-reads `ms.sizes[k]`.
  `LeafLUCache` + `build_leaf_lu_cache` (`containers.jl:287`,
  `solve.jl:53-67`) is the precedent for bolting a cache object on without
  changing `Matrices`, incl. exposing `build_time`/`bytes`.
- **The FGS builder is NOT reusable as-is** (`self_influence_matrices`
  `src/solve.jl:408-495`, `nonself_influence_matrices` :124-315): it is a
  scalar-strength → scalar-influence map — one column per BODY via
  `set_unit_strength!` (`solve.jl:83-99`, scalar 1.0 through
  `value_to_strength!`) and one row per TARGET via the system's `influence!`
  reduction (FLOWPanel's = `dot(gradient, normal)`,
  `FLOWPanel_abstractbody.jl:1298-1304`). Full `DerivativesSwitch` outputs
  and multi-component strengths are not representable. A NEW builder is
  required (Section 4); only the container and the column-fill *pattern* are
  reused.
- **`InteractionList{TF}` is DEAD CODE** (`containers.jl:233-238`, ctor
  `interaction_list.jl:682-708`, exported `FastMultipole.jl:138`): its ctor
  calls `add_influence_matrix!`, defined nowhere — any invocation throws.
  Prior art only. Do NOT silently repurpose the exported name; flag removal
  to Ryan in the PR/commit message.
- **Evaluation paths to replace** — per direct-list entry
  `(i_target, i_source)` both paths make the single call
  `direct!(target_buffer, target_index, derivatives_switch, source_system,
  source_buffer, source_index)`:
  serial `nearfield_singlethread!` → `nearfield_loop!` →
  `nearfield_direct_list_loop!` (`src/fmm.jl:3-49`); threaded
  `nearfield_multithread!` (`fmm.jl:129-207`) via `make_direct_assignments!`
  (`fmm.jl:66-110`, requires target-sorted list, cuts thread boundaries only
  at target-branch changes → each thread owns whole target branches) and
  `execute_assignment!` (`fmm.jl:112-127`). Call sites inside `fmm!`:
  `fmm.jl:1253` (serial), `:1324` (threaded). Standalone no-tree entry:
  `src/direct.jl:25-45` (`direct!`), `direct_singlethread!` :48-92,
  `direct_multithread!` :94-140 (full index ranges, no lists).
  **`direct!` ACCUMULATES (`+=`) into target output rows** — the cached
  evaluator must accumulate too.
- **Buffer layouts** (`src/derivativesswitch.jl`, contract at
  `src/compatibility.jl:6-22`): target buffer rows = 1:3 position,
  `metadata_range` 4:3+NM, then `scalar_potential_index`/`gradient_range`
  (3)/`hessian_range` (9); `output_range(switch)` = the rows a cached block
  must produce; `reset!` zeroes exactly `output_range`
  (`compatibility.jl:712`). FLOWPanel panels: NM=2
  (`FLOWPanel_abstractbody.jl:1212`). Source buffer: rows 1:3 position, 4
  radius, **5:4+strength_dims strengths**, vertices after
  (`get_strength` `compatibility.jl:314`; FLOWPanel writes strengths at
  `FLOWPanel_abstractbody.jl:1170-1174`, vertices next, `strength_dims` =
  `size(strength,2)` :1208).
- **Linearity (the correctness premise)**: FLOWPanel `induced`
  (`FLOWPanel_elements_fmm.jl:236-259`) reads
  `strength = SVector{NK}(view(source_buffer, 5:4+NK, i_source))`; `_induced`
  and `_self_limit` are linear in it. **VERIFY-FIRST item V0**: confirm
  `_induced_wake` (`:256`, attached-wake term for shedding panels) is ALSO
  linear in the SAME buffer strength rows (expected — the dense Backslash `G`
  exists and is linear in μ) and reads no independent strength-like inputs
  from the additional buffer rows
  (`additional_source_system_to_buffer!`, `FLOWPanel_liftingbody.jl:968`).
  If any system carries independent strength inputs outside rows
  5:4+strength_dims, v1 REFUSES to cache that system (loud error), rather
  than silently mis-caching.
- **`direct_conditioning`** (`src/direct_conditioning.jl`;
  applied at `fmm.jl:26-32`/`:138-144` and `direct.jl:73-87`/`:119-135`)
  mutates source buffers in place around near-field evaluation. v1:
  **refuse caching when `has_direct_conditioning(rules)`** (loud error at
  build AND at evaluation) — conditioning may change non-strength rows.
- **Threading/determinism**: target-sorted direct lists put all of a target
  branch's entries contiguous (`sort_by_target`,
  `interaction_list.jl:601-638`; applied at `fmm.jl:1007`/`:1086`).
  Target-major keying → owner-partitioned parallel accumulate with NO
  serial scatter (deliberately opposite to FGS's source-major
  `nonself_matrices`, whose scatter is serial, `solve.jl:868`). Fixed
  owner-major order ⇒ deterministic at any thread count.
- **Memory-guard precedent**: FLOWPanel `ILUPreconditioner`
  `max_pattern_entries` (default `512 * ncells`,
  `FLOWPanel_solver.jl:1124/:1134/:1161`) — size checked BEFORE allocation,
  surfaced in metadata (`FLOWPanel_metadata.jl:179`).

## 3. Design

New FastMultipole type (name: `NearfieldInfluenceCache{TF}`; do not reuse the
exported dead `InteractionList`):

```julia
struct NearfieldInfluenceCache{TF,TDS}
    matrices::Matrices{TF}       # one block per cached direct-list entry, target-major order
    entries::Vector{SVector{4,Int32}}   # (i_target_branch, i_source_branch, i_target_system, i_source_system)
    derivatives_switches::TDS    # captured at build; per-entry switch index if multi-system
    n_out::Vector{Int}           # output rows per target body, per entry (== length(output_range(switch)))
    n_comp::Vector{Int}          # strength components per source body, per entry (strength_dims)
    strengths_scratch::Vector{TF}   # gather buffer, max block width (per-thread copies at eval)
    target_tree_id::UInt         # objectid of the trees at build (identity guard)
    source_tree_id::UInt
    n_target_bodies::Int
    n_source_bodies::Int
    build_time::Float64
    bytes::Int
end
```

- **Block definition** for entry `(i_target, i_source)` and system pair:
  `m = n_out * n_target_bodies_in_branch`, `n = n_comp * n_source_bodies_in_branch`.
  Row layout: `i_row = i_out + (i_body_local - 1) * n_out` — i.e. the block
  is the column-major flattening of the `(n_out × n_targets)` output
  submatrix, so evaluation can `reshape` the product and do
  `@views target_buffer[output_range, target_cols] .+= reshape(y, n_out, :)`
  (one strided broadcast, no scatter). Column layout mirrors it:
  `i_col = i_comp + (i_body_local - 1) * n_comp`, gathered from source-buffer
  rows `5:4+n_comp`.
- **Ordering**: entries stored in the target-sorted direct-list order
  (`sort_by_target` output) so `make_direct_assignments!`-style owner
  partitioning applies unchanged.
- **`Matrices.rhs` reuse**: sized `sum(m)` as the per-block matvec output
  scratch (same role as in FGS), giving disjoint per-entry scratch for the
  parallel accumulate. NOTE: entries of one TARGET branch may be evaluated by
  the same thread (owner partitioning), so per-entry rhs slices are
  race-free by construction.

### Build (`NearfieldInfluenceCache(target_systems, target_tree, source_systems, source_tree, direct_list, derivatives_switches; max_bytes)`)

Column-by-column unit-strength probing, generalizing the FGS pattern
(`solve.jl:281-291`) to full outputs and multi-component strengths:

1. Size pass: compute all `(m,n)`; `bytes = 8 * (sum(m*n) + sum(m))`; throw a
   descriptive `ArgumentError` if `bytes > max_bytes` BEFORE allocating
   (guard precedent above). Default `max_bytes = 4 GiB` (judgment call —
   Ryan may override).
2. `save_strengths` (`solve.jl:100-112`) on all source buffers; zero all
   strength rows.
3. For each source branch (source-sorted traversal for locality), for each
   source body, for each strength component `k`: write `1.0` DIRECTLY into
   source-buffer row `4+k` of that body (bypass `value_to_strength!` — it
   implements the FGS scalar convention, e.g. FLOWPanel's zeroes the source
   column, `FLOWPanel_liftingbody.jl:1015-1018`, which would silently drop
   the σ column); for each direct-list entry with that source branch:
   `reset!` the entry's target rows, `direct!(target_buffer, target_index,
   switch, source_system, source_buffer, i_body:i_body)`, copy
   `output_range` rows into the block column. Restore strengths.
   Cost: `n_comp ×` (one full near-field pass) — pays for itself after
   ~`n_comp + 1` cached evaluations even before the flops-per-entry gap.
4. Record `objectid` of both trees, body counts, build_time, bytes.

Multi-system pairs: loop system pairs exactly as `nearfield_loop!` does; one
block per (entry × system pair) — fold the pair into `entries`.

### Evaluation (`nearfield_matvec!(target_buffers, cache, source_buffers; multithread)`)

1. Guards (cheap, performance-first): `objectid(tree) == stored id` for both
   trees, body counts match; else a descriptive error naming the contract
   ("rebuild the cache for a new Tree"). NO geometry fingerprinting — Tree
   identity IS the validity key (Ryan's spec); callers own invalidation.
2. Per entry (owner-partitioned threads over target branches, static
   assignment via the `make_direct_assignments!` pattern): gather strengths
   from source-buffer rows into the column-layout scratch,
   `mul!(rhs_slice, block, x)`, accumulate
   `@views buffer[output_range, cols] .+= reshape(rhs_slice, n_out, :)`.
   Accumulation preserves `direct!` semantics (callers still `reset!`
   targets themselves, as today).
3. Determinism: fixed entry order per owner + serial accumulate within an
   owner ⇒ bitwise repeatable at any thread count (test it).

### Integration points (each its own commit)

1. **FastMultipole core** (container + build + `nearfield_matvec!` + tests).
2. **`fmm!` + `FmmPlan`**: the `fmm.jl:1020` method gains
   `nearfield_cache=nothing`; when set, the calls at `fmm.jl:1253`/`:1324`
   become `nearfield_matvec!` (assert no `direct_conditioning`). `FmmPlan`
   (added 2026-08-18, `fmm.jl` after :1018, commit b3cf4ad) gains an optional
   cache slot + `build_nearfield_cache!(plan, systems)` — the natural owner,
   since plan lifetime ≡ tree lifetime. This composes with the Krylov
   `cache_tree` lever: within one solve the near-field becomes matvecs.
3. **Standalone `direct!`**: `direct!(target, source; cache)` degenerate
   no-tree form — one block per system pair over full index ranges (the
   "general, usable by any direct! function" requirement).
4. **FLOWPanel surface**: kwarg on `KrylovSolver`
   (`cache_nearfield::Bool=false`, requires/implies `cache_tree=true`) built
   lazily on first apply alongside the plan; metadata entries mirroring
   `cache_tree` (`FLOWPanel_metadata.jl` KrylovSolver branch). FGSSolver's
   own matrices are OUT of scope (different, scalar semantics).
5. **Timestep-to-timestep reuse** (the headline; second stage): requires
   cross-solve `FmmPlan` persistence — currently the Krylov plan is
   per-solve by design (kerneloffset safety). Extend `cache_tree` with an
   explicit opt-in `persistent_plan=true` whose contract is "caller
   guarantees frozen geometry AND fixed active kerneloffset across solves"
   (document: the Dirichlet solve path flips kerneloffset between solve and
   target passes — `_set_kerneloffsets!`; the plan+cache pair is keyed to
   the PANEL-offset state used inside `_krylov_launch!`, which is restored
   before every solve, so persistence is sound for rigid bodies whose global
   positions don't change; a ROTATING rotor changes global positions every
   step ⇒ new tree ⇒ no reuse — that is correct behavior, not a bug; frame
   tricks are explicitly out of scope).

## 4. Tests (FastMultipole `test/nearfield_cache_test.jl` + FLOWPanel solver units)

1. **Exactness, scalar system**: Gravitational (`test/gravitational.jl`)
   fmm! with cache vs without — outputs `isapprox` rtol 1e-12 (NOT bitwise:
   BLAS matvec sums in a different order than the kernel's per-source
   accumulation). Premise guards: nonempty direct list AND nonempty m2l.
2. **Exactness, panel system (the V0 linearity verification)**: FLOWPanel
   RigidWakeBody Dirichlet R1-class body (construction recipe in
   `benchmark/phase1_case.jl`) — `apply_G!` through a cached backend vs the
   uncached FMM apply, rtol 1e-12; MUST include shedding panels so
   `_induced_wake` is exercised (this test IS the wake-linearity check; if
   it fails, implement the V0 refusal path and document).
3. **Strength-change reuse**: change strengths, re-evaluate, compare vs
   fresh direct — proves only strengths are read at eval.
4. **Tree-identity guard**: evaluating with a rebuilt tree throws.
5. **Conditioning refusal**: build/eval with a `DirectConditioningRule`
   present throws.
6. **Memory cap**: tiny `max_bytes` throws before allocation (premise:
   estimator > cap).
7. **Determinism**: repeated evaluations, fixed inputs, bitwise-identical
   target buffers; run under `-t 1` and `-t 4` (in-process repeatability +
   the serial-owner argument mirrors the colored-sweep test pattern,
   `test/fgs_coloring_test.jl`).
8. **FLOWPanel**: `cache_nearfield` KrylovSolver solve == uncached solve
   (rtol 1e-12), metadata round-trip, DirectBackend combo rejected.

## 5. Performance acceptance (021 Phase 2b A/B, after implementation)

- Isolated: cached near-field apply vs kernel near-field at R1–R3 campaign
  knobs (expect ≥5× on the panel kernel; report measured + build time +
  bytes + break-even apply count).
- End-to-end: `benchmark/rotor_hover_solver_phase2.jl` rows for
  krylov_{gmres,jacobi,ilu} and fgmres_fgs with `cache_nearfield` on/off;
  certified BC unchanged (rtol-level solution shifts are expected — record
  them); ledger before/after rows.
- Memory column: cache bytes vs the dense-G 8N² baseline (the cache is the
  interesting middle point between matrix-free and matrix-ful — a
  publishable axis per ruling 8).

## 6. Judgment calls already made (Ryan may override)

- Refuse-and-error over silent fallback (cap exceeded, conditioning present,
  tree mismatch): performance-first means no hidden kernel-path surprises.
- Tree `objectid` identity as the whole validity key; no fingerprinting.
- `value_to_strength!` bypassed in the builder (direct buffer-row writes) so
  ALL strength components get columns — the FGS convention would drop some.
- v1 excludes FGS-internal matrix replacement and any rotating-frame
  invariance tricks.
- Dead `InteractionList` left in place but flagged for removal (exported
  name; Ryan's call).

## 7. Suggested commit split

1. FastMultipole: `NearfieldInfluenceCache` + builder + `nearfield_matvec!`
   + standalone `direct!` integration + tests 1/3–7.
2. FastMultipole: `fmm!`/`FmmPlan` integration + test (fmm-with-cache ==
   without).
3. FLOWPanel: `cache_nearfield` on KrylovSolver + metadata + tests 2/8.
4. (After Ryan's go on persistence) cross-solve `persistent_plan` +
   timestep-reuse test.
