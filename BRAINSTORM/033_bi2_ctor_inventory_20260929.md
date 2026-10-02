# B-I2 prep: FGS constructor stage inventory (2026-09-29)

Code-trace of the `FastGaussSeidel` constructor
(`FastMultipole/src/solve.jl:846–1039`, live checkout at fm ≥ b4c35f67),
made when staging B-I2. Line numbers are as of 2026-09-29 — re-verify before
editing. Context: after B-I1 threaded the probes, the remaining ctor tail is
~20 of 30.5 s at R4 j64 and ~53 of 85.8 s at R5 (composition UNMEASURED —
profiling it is B-I2 step 1).

| # | Stage | Location | Threaded today? | Notes |
|---|---|---|---|---|
| 1 | Octree build (source + target replay) | `solve.jl:878,880` → `tree.jl:14/202` | Partially | `child_branches_multithread_parents!` (`tree.jl:547`, `@threads` 564/616), `sort_bodies_multithread!` (`tree.jl:1190`, `@threads` 1204/1230/1257), `shrink_recenter_*_multithread!` (`tree.jl:1844/1877`); serial residue = root setup / cumulative-offset bookkeeping |
| 2 | `assert_shared_topology` | `solve.jl:881` (fn 826) | Serial | cheap sanity check |
| 3 | Interaction-list build | `solve.jl:895` → `interaction_list.jl:3` (recursion at 69ff) | **Serial** | recursive dual-tree walk; no threaded variant exists |
| 4 | List sorts | `solve.jl:905–908` | Mixed | `sort_by_source` serial (`interaction_list.jl:640`); `sort_by_target` auto-threads (`:601–604` → `:543`, `@threads :static`). Mirror the threaded twin for `sort_by_source` |
| 5 | Non-self influence probes | `solve.jl:913` (fn 143) | **Threaded (B-I1)** | `setup_threads>=1` branch at 248 → `_populate_nonself_threaded!` (375, `@sync`/`@spawn` 404–405) |
| 6 | `add_self_interactions` | `solve.jl:918` (fn 801) | Serial | plain loop |
| 7 | `index_by_source` | `solve.jl:922` (fn 491) | Serial | plain loop |
| 8 | Self influence probes | `solve.jl:926` (fn 577) | **Threaded (B-I1)** | `setup_threads>=1` at 596 → `_populate_self_threaded!` (679) |
| 9 | Leaf LU cache | `solve.jl:927` → `build_leaf_lu_cache` (72–86) | **Serial** | `map` + `lu!` over independent square leaf blocks — embarrassingly parallel. Only ctor stage with an existing timer (`build_time`, lines 73/83/85) |
| 10 | strengths / `map_by_leaf` / `map_by_branch` | `solve.jl:931–936` (fns 756, 783) | Serial | bookkeeping loops |
| 11 | RHS / influence allocations | `solve.jl:940–951` | Serial | trivial |
| 12 | Residual vector sizing | `solve.jl:953–963` | Serial | trivial |
| 13 | `color_leaves` (`:colored` only) | `solve.jl:967–973` (fn 1388) | Serial | not on production path |
| 14 | `build_chunk_map` (`:chunked` only) | `solve.jl:977–985` (fn 1307) | Serial | not on production path |
| 15 | dagteam plan build (production default) | `solve.jl:989–1006` → `build_dagteam_plan` (`solve_dagteam.jl:131–335`) | **Fully serial** | edge derivation 167–225, split-storage sizing 193–219, **L/U repack incl. F32 conversion 227–266** (second serial pass over all near-field bytes — prime suspect), sweep-precision LU refactorization (`build_leaf_lu_cache_as`, 80–89), critical-path priorities 274–281 |
| 16 | Struct assembly | `solve.jl:1008–1038` | Serial | trivial |

Inherently sequential: only the critical-path priority pass
(`solve_dagteam.jl:274–281`, reverse-topological over ~1k leaves) —
negligible. Everything else costly is per-leaf/per-pair independent work.

Instrumentation gap: apart from stage 9's `build_time`, the ctor has NO
per-stage timers; the `diagnostics[:...]` keys (`dagteam_*_ns` etc.,
`solve.jl:1543–1871`, `solve_dagteam.jl:759–859`) are all solve-time, not
ctor-time. B-I2 step 1 = add `time_ns()` brackets per stage above.

NF-cache non-fact: `NearfieldInfluenceCache` (`nearfield_cache.jl`, threaded
`@spawn` chunk queue at :463) is used only by `autotune.jl`/`direct.jl`/
`fmm.jl` — the FGS ctor does NOT use it; FGS near-field comes solely from
stages 5/8.

ILU side (for B-I3): `ILUPreconditioner` (`FLOWPanel
src/FLOWPanel_solver.jl:2295ff`) stats dict records `tree_time`,
`assembly_time`, `factorization_time`, `total_time`; assembly is threaded
(`@threads` at `:2219`), nfcache build threaded, serial residue =
`_ilu_direct_pattern` tree/lists + `sparse()` CSC assembly + `ILUZero.ilu0`
(`:2348`). Measured split at R4 j64 (warm-start R4 campaign): t_setup
≈ 80–82 s, t_prime ≈ 28 s.
