# Implementation prompt — BRAINSTORM 030 Phase 1+2 (2026-09-09)

You are implementing BRAINSTORM item 030 (generic influence-block assembly) in
FastMultipole, then wiring FLOWPanel to use it. Read
`BRAINSTORM/030_generic_influence_block_assembly.md` FIRST — it is the design of
record (frozen; do not redesign). Then read `agent_policies/WORKFLOW.md` and
`agent_policies/TESTING.md` before touching code. Never use more than 4 threads
locally. Local work only — no cluster deploys, no Slurm.

## Repos and state — WORK ONLY IN THE 030 WORKTREES

- Your working trees are `~/wt030/FLOWPanel.jl` (branch `030-block-assembly`) and
  `~/wt030/FastMultipole` (branch `030-block-assembly`), side by side so
  FLOWPanel's relative Manifest dev-paths (`../FastMultipole`, `../FLOWVPM.jl`)
  resolve inside `~/wt030/`. A `~/wt030/FLOWVPM.jl` worktree exists purely as a
  read-only dependency — never edit it. Run all julia commands with
  `--project=~/wt030/FLOWPanel.jl` from within the worktrees.
- NEVER touch the live checkouts under `~/Dropbox/research/projects/` — other
  agents are working there concurrently (052e in FastMultipole, others in
  FLOWPanel/FLOWVPM). The single exception: appending your dated Log entry to
  `BRAINSTORM/030_generic_influence_block_assembly.md` in the live FLOWPanel
  checkout at the end (append-only, that file only).
- The prerequisite work from 2026-09-09 is COMMITTED and already in your
  worktrees: threaded near-field cache build (worker pool + `Threads.Atomic`
  chunk counter, private buffer copies per worker) and `NearfieldCacheDonor`/
  `retarget_nearfield_cache` in `FastMultipole/src/nearfield_cache.jl` (+ exports,
  + `n_threads` plumbing in `src/fmm.jl`) — FastMultipole commit `da1bd13a`; donor
  plumbing in FLOWPanel (`src/FLOWPanel_solver.jl` KrylovOperator/KrylovSolver,
  `src/FLOWPanel_fmm.jl` influence!, `src/FLOWPanel_instrumentation.jl`
  `_apply_*_G!`), `NF_DONOR` in `benchmark/rotor_hover_solver_phase2_tune.jl`,
  and new testsets (`FastMultipole/test/nearfield_cache_test.jl` parallel +
  retarget; `test/runtests_unit_solver.jl` donor testset). All tests green at
  handoff; re-run the FastMultipole cache test file and the FLOWPanel solver
  unit tests once at the start to confirm your worktree baseline before editing.

## Key code to read before writing anything

- `FastMultipole/src/nearfield_cache.jl` — the whole file. The probe loop is
  `_probe_source_keys!`; the builder is `_build_nearfield_cache` (worker pool,
  guards). Blocks: `entries[k] = (i_target_branch, i_source_branch, i_ts, i_ss)`,
  column j = (i_body − first(source_range))·sd + i_comp, rows =
  vec(target_buffer[out_range, target_range]).
- `FastMultipole/src/solve.jl:3-52` — `Matrices` packed storage,
  `get_matrix_vector`.
- FLOWPanel `src/FLOWPanel_solver.jl:237` (`_G!`) — the assembly pattern being
  promoted: unit-activate strengths, per-pair `induced`, write entries.
- FLOWPanel `src/FLOWPanel_abstractbody.jl:1308` (`direct!` → `_direct_body!`) —
  the per-pair `induced(target, source_system, source_buffer, i_source, switch,
  fam; core_size)` primitive and the `Val(FILAMENT_REGULARIZATION[])` function
  barrier (read it ONCE per block; per-pair reads measured +34–49% regression,
  see the comment there).
- Test conventions: `FastMultipole/test/nearfield_cache_test.jl` (esp. "parallel
  build is bit-identical" and the rtol-1e-12 cached-vs-kernel pattern) and
  FLOWPanel `test/runtests_unit_solver.jl` testset "KrylovSolver cache_nearfield
  (021 Phase 2b)" (shedding-body premise guards, `make_dirichlet_diamond_body`,
  `make_sphere_source_body`).

## Phase 1 — FastMultipole

1. Add `assemble_influence_block!(block, target_buffer, target_range, switch,
   source_system, source_buffer, source_range)` with a DEFAULT method that
   reproduces the current unit-strength probe for that one block (zero output
   rows → `direct!` per source column → copy column). Semantics: every entry of
   `block` is ASSIGNED. Export it.
2. Restructure `_build_nearfield_cache` so the per-block inner work goes through
   the hook. Preserve exactly: worker-pool parallel structure, key grouping,
   strength save/zero/restore around the build, size and `max_build_time` guards,
   bit-identical results at any thread count. The fallback path still needs the
   per-worker private target+source buffer copies (probing writes target output
   rows and pokes strengths, and source branch body ranges OVERLAP across tree
   levels — this raced before and was fixed; do not reintroduce it). An opted-in
   system's blocks must not require the copies; simplest correct structure:
   detect per-source-system whether an override exists (e.g.
   `hasmethod`-based trait computed once per build) and route each block.
3. Switch `Matrices(sizes, TF)` to calloc-backed zeros. Check callers don't rely
   on undef-then-fill semantics (FGS `build_leaf_lu_cache` copies data; fine).
4. Opt-in override for the gravitational TEST system (in `test/gravitational.jl`
   or the test file — it is a test-only system; do not pollute src with it).
5. Tests (extend `test/nearfield_cache_test.jl`):
   - assembled vs probed block data, rtol 1e-12, on a case with a real direct
     list (reuse the existing premise guards);
   - mixed build: two source systems, one opted-in + one fallback, cache matvec
     matches the kernel near-field to rtol 1e-12;
   - parallel bit-identical testset still passes with the hook in place
     (n_threads=1 vs 4) — for the fallback AND for the opted-in system
     separately (opted-in assembly must also be deterministic across thread
     counts: entries assigned independently, no accumulation across tasks);
   - existing 4 testsets in the file unchanged and green.
6. Benchmark (report numbers, don't commit the script): gravitational n=20k,
   leaf 40, MAC 0.5, pot+grad — serial and 4-thread build, hook vs probe.
   Expectation from the design: ~3–5× serial for the opted-in path; flag if not
   met rather than tuning silently.

## Phase 2 — FLOWPanel

1. Opt-in `assemble_influence_block!` overload(s) for `AbstractBody` covering
   both output forms the cache stores (Neumann: gradient rows; Dirichlet:
   scalar-potential rows — follow what `output_range(switch)` selects, do NOT
   hand-roll the projection; the cache stores RAW output rows, not φ/u·n̂
   projections). Build on `induced` with the `Val` family barrier hoisted once
   per block, mirroring `_direct_body!`'s loop structure. Respect
   `strength_dims` for the column layout; note `_direct_body!` accumulates
   φ/U/H per target — the override assigns per (target, source) pair instead.
   Careful: `induced` at self pairs returns the side-aware self limit — keep it
   (the probe path gets the same values through `direct!`, so exactness tests
   will catch any deviation).
2. Exactness tests in `test/runtests_unit_solver.jl` (new testset, mirror the
   existing cache_nearfield one): assembled-vs-probed cache on
   `make_dirichlet_diamond_body(nspan=40)` (shedding panels present — premise
   guards) and `make_sphere_source_body` (Neumann), rtol 1e-12; then a full
   cached KrylovSolver solve must match the probe-built solve to rtol 1e-12
   (NOT bitwise — assembly may sum in a different order, so the blocks can
   differ in last bits; if they happen to be bitwise equal, still assert only
   rtol). The donor-retarget testset must stay green (retarget is orthogonal —
   it rebinds whatever blocks were built).
3. Run: full `FastMultipole/test/nearfield_cache_test.jl` (via the include
   prelude used by its runtests: gravitational/vortex/vortex_filament/panels),
   FLOWPanel `julia --project -t 4 test/runtests_unit_solver.jl` (all 465+ must
   pass), and `julia --project -t 4 test/runtests_unit_replay.jl` if WORKFLOW.md
   routes cache-path changes there (check TESTING.md matrix; when in doubt run
   the narrower test and say which you ran).

## Ground rules

- Do not change `direct!` or any kernel signatures; the hook is additive.
- Do not touch the FGS builders (`solve.jl` self/nonself matrices) — Phase 3 is
  Ryan-gated.
- Do not modify the 021 benchmark drivers beyond what already exists.
- Sign conventions are the #1 regression source in FLOWPanel (see CLAUDE.md);
  the rtol-1e-12 exactness tests are the guard — if they fail, suspect YOUR
  assembly loop, not the probe.
- Commit on the `030-block-assembly` branches in the worktrees (FastMultipole:
  hook + calloc + tests; FLOWPanel: overloads + tests; clean, separate commits)
  but DO NOT push and DO NOT merge into `flowpanel-20260817` or `fastmultipole`
  — merging back is Ryan-gated. List the commits in your report.
- Report at the end: what changed (file:line), test results verbatim totals,
  benchmark table (probe vs hook, serial vs 4T), deviations from this prompt,
  and open questions for Ryan. Update the 030 item's Log section with a dated
  entry (append-only).
