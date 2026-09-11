# 030 — context-reset prompt after Phases 1+2 COMPLETE (2026-09-10, session 2)

You are picking up BRAINSTORM 030 (generic influence-block assembly) with **Phases
1 and 2 finished, tested, and committed**. Read FIRST, in order:
`BRAINSTORM/030_generic_influence_block_assembly.md` (design of record, frozen —
its Log section now carries both implementation-session entries),
`BRAINSTORM/030_implementation_prompt_20260909.md` (original instructions), then
this file. Work ONLY in the `~/wt030` worktrees (`~/wt030/FastMultipole`,
`~/wt030/FLOWPanel.jl`, both on branch `030-block-assembly`); never touch the live
checkouts under `~/Dropbox/research/projects/` except the single allowed append to
the 030 item's Log at the very end of your session. Local only, max 4 threads,
commit but do not push or merge (merge-back is Ryan-gated).

## State at handoff — everything below is DONE and committed

Branch `030-block-assembly`, FastMultipole: `da1bd13a` → `dd70fb19` (Phase 1) →
`fff72d29` (hook speedup). FLOWPanel: `8dce66c` → `6798197` (Phase 2) → `431ea46`
(session-1 handoff doc). Worktrees clean at handoff.

- **Phase 1** (`dd70fb19`): `assemble_influence_block!` hook + `_probe_influence_block!`
  default + `overrides_block_assembly` which-based trait (compute once per source
  system per build — per-block reflection cost 2.4×) + calloc-backed
  `Matrices(sizes, TF)` + gravitational test-system opt-in +
  `use_block_assembly=false` diagnostic knob on both cache constructors.
- **Hook speedup** (`fff72d29`, Ryan-approved pass): `_assemble_source_keys!` hands
  overrides a plain-`Matrix` via `unsafe_get_block_matrix` (`src/solve.jl`,
  `GC.@preserve`d `unsafe_wrap`; must not escape the builder) instead of the
  `ReshapedArray{SubArray}` from `get_matrix_vector`; grav override computes
  `rinv2`/`rinv` once (one division + one sqrt per pair, `@fastmath`).
- **Phase 2** (`6798197`): `FastMultipole.assemble_influence_block!(... ::AbstractBody ...)`
  in `src/FLOWPanel_abstractbody.jl` (after `_direct_body!`) → generic
  `_assemble_influence_block!` with `Val(FILAMENT_REGULARIZATION[])` hoisted once
  per block; per-source column copied into a single-column scratch, strengths
  unit-activated per component, rows in `output_range(switch)` order, TS/extra
  rows zeroed to match the probe. `PanelWake`/`FilamentWrapper` keep the probe
  (own `direct!`, not `AbstractBody` — by design). New testset
  "assemble_influence_block! opt-in (030 Phase 2)" in
  `test/runtests_unit_solver.jl` (24 tests: diamond Dirichlet WITH shedding,
  sphere Neumann, operator applies, full cached KrylovSolver assembled vs
  probe-swapped on a persistent plan).

**Test totals (verbatim, all green)**: FastMultipole cache test file
25+8+17+4+6+9443+7+8 = **9518/9518** (driver = runtests prelude:
gravitational/vortex/vortex_filament/panels, then nearfield_cache_test.jl, run
`julia --project=$HOME/wt030/FLOWPanel.jl -t 4 <driver>`). FLOWPanel
`julia --project=. -t 4 test/runtests_unit_solver.jl` → **Solvers | 489 489**
(465 baseline + 24 new). `runtests_unit_replay.jl` not required (replay does not
touch the near-field cache path; TESTING.md routing).

**Benchmark of record** (grav n=20k, leaf 40, MAC 0.5, pot+grad, max_bytes 6 GiB,
min-of-3; machine SHARED with another 4-thread julia job throughout — both
sessions' numbers carry the same confound):

| path  | serial (s) | 4T (s) |
|-------|-----------|--------|
| probe | 1.745 | 0.602 |
| hook  | 0.602 | 0.253 |

Serial 2.90×, 4T 2.38×. **Attribution (measured, do not re-derive)**: assembly
proper on pre-touched pages is 0.293 s = 15.8 GB/s vs a 68 GB/s single-thread
warm-fill floor — COMPUTE-bound (~2 ns/pair sqrt+div); the residual ~0.3 s of the
serial build is first-touch page faults on the 4.31 GiB calloc region, paid
identically by both paths. So ~3× is the practical ceiling for this kernel and
the design's 3–5× expectation over-counted the shared fault cost. Ryan reviewed
this attribution 2026-09-10 and the follow-up: the fault cost is operationally
irrelevant — the phase2 tuner EXCLUDES build cost from its objective by
construction (untimed warm-up solve; Ryan's ruling), production reuses the cache
across timesteps via persistent_plan/transform_plan!, and tuner candidates adopt
blocks via `retarget_nearfield_cache`, which ALIASES donor storage (no new
allocation, no faults). The "reuse Matrices allocation across rebuilds" idea was
considered and DROPPED as unmotivated.

## Root causes & gotchas established this session (do not re-derive)

- **The "sphere FmmPlan hang" was NOT an 030 regression**: `make_sphere_source_body`
  (test_helpers.jl) never calls `calc_controlpoints!`, so planning on a raw body
  targets ALL-ZERO controlpoints; with >leaf_size coincident targets the octree
  target subdivision recurses without bound (`tree.jl:520` guard
  `exceeds(...) && (target || child_radius >= max_body_radius)` — the radius stop
  applies to SOURCE trees only, and no depth guard exists). Baseline was green
  because every prior sphere test reached the plan through `solve!`, which
  initializes geometry. Fixed in the testset (initialize before planning);
  **tree.jl deliberately untouched** — the missing depth/degeneracy guard is a
  pre-existing upstream gap flagged to Ryan, unfixed.
- Neumann exactness: assembled vs probed sphere cache measured BITWISE identical
  (worst rel diff 0.0); asserted at rtol 1e-12 per the design (bitwise is not
  promised).
- `unsafe_get_block_matrix` wrappers must stay inside the builder loop's
  `GC.@preserve`; do not let them escape.
- Background julia runs: pipe to a FILE, not `| tail`. `timeout` does not exist
  on this Mac — use a background watchdog (`kill -0` loop) instead.
- Benchmark needs `max_bytes = 6 GiB`; machine has 8 cores/16 GiB, check
  `ps aux | grep julia` for competing agent jobs before timing anything.

## Remaining work — ALL Ryan-gated; do NOT start any without his explicit go

1. **Phase 3**: migrate FGS `self_influence_matrices` / `nonself_influence_matrices`
   (`FastMultipole/src/solve.jl:408` / `:124`) to the hook. Open design question
   (decide WITH measurements, per the design): projected-hook variant vs
   projecting the assembled raw block through `influence!`. Note the FGS builders
   store post-projection scalars, not raw rows. The `unsafe_get_block_matrix`
   wrapper will benefit either route.
2. **Merge-back** of `030-block-assembly` into `fastmultipole` (FastMultipole) /
   the FLOWPanel working branch — Ryan decides when; watch for conflicts with
   052e (another agent works in live FastMultipole).
3. **Upstream robustness** (optional, Ryan's call): depth cap or duplicate-point
   error in `tree.jl` so degenerate geometry fails loudly instead of hanging.
4. Item-level approval checkboxes in the 030 item remain unticked (Ryan's).

## Log discipline

The 030 item's Log in the LIVE checkout
(`~/Dropbox/research/projects/FLOWPanel.jl/BRAINSTORM/030_generic_influence_block_assembly.md`)
already has the session-2 completion entry (2026-09-10). Append-only, that file
only, dated entries under the Log section.
