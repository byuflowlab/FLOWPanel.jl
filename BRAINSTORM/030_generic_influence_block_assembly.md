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
- 2026-09-10 (implementation session 1, wt030 worktrees): **Phase 1 COMPLETE** —
  FastMultipole `030-block-assembly` commit `dd70fb19`: hook + probe default +
  `overrides_block_assembly` trait + calloc `Matrices` + gravitational opt-in +
  tests (nearfield_cache_test.jl 9518/9518, incl. mixed grav+vortex build and
  parallel bit-identical for both paths). Benchmark (grav n=20k, leaf 40, MAC 0.5,
  pot+grad): probe 2.06 s / hook 0.75 s serial (2.74×), 0.76 / 0.36 s at 4T
  (2.13×) — below the 3–5× design expectation, flagged not tuned (memory-bound
  4.3 GiB write; gotcha: trait must be computed per SYSTEM, per-block `which()`
  cost 2.4×). **Phase 2 WIP, BLOCKED** — FLOWPanel `030-block-assembly` commit
  `f2e28fe`: `AbstractBody` overload (scratch-column unit strengths over
  buffer-based `induced`, fam barrier per block) verified BITWISE-equal to the
  probe on the tiny Dirichlet diamond (0.87 ms vs 195 ms build), but FmmPlan
  Tree construction now HANGS for sphere bodies (`sort_bodies!` runaway,
  ~1e11 allocs) before any cache code runs — root cause open, Neumann exactness
  and full test suite pending. Entry point for the next agent:
  `BRAINSTORM/030_reset_prompt_20260910.md` (in the wt030 worktree, commit
  `f2e28fe`).
- 2026-09-10 (implementation session 2, wt030 worktrees): **Phase 2 COMPLETE**.
  Hang root cause: NOT the session's edits — `make_sphere_source_body` never calls
  `calc_controlpoints!`, so the WIP testset planned on a body whose control points
  were ALL ZEROS; with >leaf_size coincident targets the octree target subdivision
  recurses without bound (pre-existing FastMultipole gap: `tree.jl` `exceeds` has
  no depth/degeneracy stop for target trees — the radius stop applies to source
  trees only). Baseline was green because every prior sphere test reached the plan
  through `solve!`, which initializes geometry. Fix: initialize the sphere in the
  testset (tree.jl untouched); flagged upstream gap for Ryan. Neumann exactness:
  assembled vs probed sphere cache measured bitwise-identical (worst rel diff 0.0),
  operator applies identical. **Hook speedup pass (Ryan-approved)**: serial
  2.37×→**2.90×**, 4T 2.21×→2.38× (probe 1.745/0.602 s serial, 0.602/0.253 s 4T;
  same benchmark, machine shared with another 4T job throughout). Changes
  (FastMultipole commit `fff72d29`): builder hands overrides a plain-`Matrix`
  `unsafe_wrap` per block instead of `ReshapedArray{SubArray}`; grav override uses
  one division + one sqrt per pair (`@fastmath`). Attribution: assembly proper
  15.8 GB/s serial vs 68 GB/s warm-fill floor — compute-bound (sqrt+div ~2 ns/pair);
  the remaining ~0.3 s residual is first-touch page faults on the 4.31 GiB calloc
  region, paid identically by both paths, so 3× is the practical ceiling for this
  kernel and the 3–5× design expectation over-counted the shared fault cost.
  Tests: FastMultipole cache file 9518/9518; FLOWPanel `runtests_unit_solver.jl`
  **489/489** (465 baseline + 24 new). `runtests_unit_replay.jl` not run — replay
  does not touch the near-field cache path (TESTING.md routing). Final commits:
  FLOWPanel `6798197` (Phase 2) + `431ea46` (handoff doc), FastMultipole
  `fff72d29`, all on `030-block-assembly`, not pushed/merged (Ryan-gated).
  Phase 3 open question: reuse the projected-hook variant for FGS migration.
- 2026-09-10 (Ryan rulings, session 2 follow-up): **Phase 3 projection question
  DECIDED — route A** (project-after-assembly: reuse the exact raw
  `assemble_influence_block!` hook, collapse rows through the `influence!`-owned
  projection via a small adapter; per-worker raw-block scratch is L2-resident, so
  the round-trip is expected to be noise against the compute-bound kernel;
  `influence!` stays the single owner of sign-sensitive projection semantics).
  **New Phase 3b STAGED (Ryan directive): far-pair point-panel approximation
  prototype** — the µs-per-pair panel kernel is the FGS/cache build cost for
  FLOWPanel bodies, and every direct-list pair pays it regardless of separation
  (core_size radius inflation makes much of the interaction direct; recon
  2026-09-10: NO existing point-panel approximation in active kernels — only
  commented legacy code in `FLOWPanel_elements.jl`; FastMultipole's panel
  multipole expansions act at branch level only). Prototype: per-pair switch in
  the 030 assembly hook — target farther than η·(panel characteristic length) →
  point source σA (ConstantSource) / point dipole μA·n̂ (ConstantDoublet,
  VortexRing ≡ constant-doublet panel) instead of the full integral; error
  O((L/r)²), η tuned so the error is commensurate with FMM p=8 truncation
  (~1e-8). Exact mode stays the default and keeps the rtol-1e-12 tests; the
  approximation is opt-in behind a knob and gets its own error-bound and
  solve-level acceptance tests. Details in the reset prompt
  (`030_reset_prompt_20260910c.md`, wt030 worktree).

- 2026-09-10 (implementation session 3, wt030 worktrees): **Phase 3 COMPLETE
  (route A), STOP-rule triggered; Phase 3b prototype COMPLETE, adoption not
  recommended without a Ryan ruling.** Phase 3: FGS `self_influence_matrices`
  / `nonself_influence_matrices` route opted-in systems through
  `assemble_influence_block!` (raw block into a reusable max-size scratch,
  per-body strength-component columns weighted by the unit-value strength
  vector, written into the target buffer's output rows, collapsed through the
  SAME `influence!` call as the probe); fallback systems keep the probe
  verbatim; `use_block_assembly` knob on both builders, `FastGaussSeidel`,
  and `FGSSolver`. Exactness: grav self+nonself and FLOWPanel diamond
  (shedding, self, fixed-iteration solve trajectory) + sphere (multi-leaf,
  nonself) all rtol 1e-12; FastMultipole driver 597435/597435, FLOWPanel
  `runtests_unit_solver.jl` 504/504 at Phase 3. **Measured (machine shared
  with a 4T 052e job, same confound as the benchmark of record): grav FGS
  build hook 0.99–1.02x vs probe; (hook−raw)/hook = 86% — the FGS build is
  dominated by the per-body projection bookkeeping (reset!/influence!/copy)
  that route A retains, so the >~10% stop rule fired: NO projected-hook
  variant built, Ryan decides.** Second finding: for NK=2 panel bodies the
  hook is ~1.9x SLOWER than the probe on the self matrix (route A evaluates
  `induced` once per strength component; the probe's `direct!` does both in
  one pass) — production rotor FGS bodies (VortexRing, NK=1) are at parity;
  scratch round-trip itself is noise for µs kernels (−4%). Phase 3b
  (far-pair point-panel, exact-by-default `FARFIELD_ETA[]=Inf` +
  `with_farfield_eta` + `FGSSolver(farfield_eta=...)`): point kernels
  calibrated against `_induced` (phi = −A/4π(σ/r + µ n̂·r/r³), U = +∇phi,
  GT normal; VortexRing ≡ constant doublet exactly); attached wake always
  exact; HS blocks always exact. Sweep (η ∈ {2,3,4,5,8}): sphere block error
  ≈ 0.14/η² (clean (L/r)² trend), solve-level strength diff 8.5e-4→6.1e-5,
  cache build 12.2x→1.9x, far fraction 97%→50%; per-pair point 76 ns vs
  exact 370–620 ns (~5x). **The η target commensurate with FMM p=8 (~1e-8)
  is UNREACHABLE with a monopole/dipole far field (needs η~4e3); and on the
  motivating all-direct shedding fixture (diamond) only 1.2% of pairs are
  far at η=2, 0% at η≥3 — the <~30% clause fired.** Prototype shipped with
  honest-bound tests (rtol-1e-12 exactness at Inf untouched; sphere η-bound
  + trend + solve-level acceptance; diamond far path with wake). Commits:
  FastMultipole `417489d5` (Phase 3); FLOWPanel `36428aa` (Phase 3),
  3b `bc139c3`; all on `030-block-assembly`, not
  pushed/merged (Ryan-gated). Entry point for details:
  session report + `030_reset_prompt_20260910c.md`.

- 2026-09-11 (Ryan ruling): **Phase 3 (route A FGS migration) and Phase 3b
  (far-pair point-panel dipole prototype) RETIRED** — route A does not help the
  FGS build (0.99-1.02x cheap kernel, 1.9x SLOWER for NK=2 panel bodies; the
  build is projection-bookkeeping-bound) and the dipole far field can neither
  reach useful accuracy (0.14/eta^2) nor find far pairs on the shedding
  fixtures that motivate it. Reverted on the worktree branches: FLOWPanel
  `0c444a4` (reverts `bc139c3`) + `3fdf5af` (reverts `36428aa`); FastMultipole
  `00bd5f5d` (reverts `417489d5`). Prototypes remain recoverable from those
  commits. Phases 1+2 (near-field cache hook, 2.90x serial) stay SHIPPED on
  `030-block-assembly`. Successor exploration: BRAINSTORM item 031
  (`031_quadrupole_panel_farfield.md`) — quadrupole/higher-order panel
  far-field, opened 2026-09-11.

- 2026-09-12 (merge-back session, per `030_mergeback_prompt_20260912.md`):
  **Stage 1 COMPLETE, Stage 2 BLOCKED on preconditions — production untouched.**
  Production identified as: FLOWPanel `fastmultipole` @ `8dce66c`, FastMultipole
  `flowpanel-20260817` @ `da1bd13a`. Both production tips are already ancestors
  of the dev branches (production had NOT moved since the worktrees were cut),
  so both Stage 1 merges returned "Already up to date" — no merge commits, no
  conflicts; Stage 2 would be a fast-forward of each live branch to the dev tips
  (FLOWPanel `00c5fcd`, FastMultipole `00bd5f5d`). Stage 1 validation in wt030
  (Manifest verified: FastMultipole → `../FastMultipole` = wt030 worktree),
  verbatim totals:
  - FastMultipole cache/FGS driver: **597299/597299** (matches pre-merge total).
  - FLOWPanel `runtests_unit_solver.jl`: **Solvers | 489 489**.
  - Full `runtests.jl` (first run against 030 code): every testset green except
    **Kutta closure (BRAINSTORM 015) | 656 2 658** — 2 failures in
    "explicit jump fallback (:jump)" at `test/runtests_unit_kutta.jl:539-540`
    (exact-equality strength checks). Attribution: re-ran that file with both
    worktrees detached at the production pins (`8dce66c`+`da1bd13a`) —
    **identical 656/658 failure**, i.e. PRE-EXISTING in production, not
    introduced by 030. Root cause not investigated (production-side).
  Stage 2 preconditions failed at check time: (1) two live Julia jobs were using
  the live checkouts (052e realsim script; a rotor-hover smoke chain); (2) live
  FastMultipole's uncommitted 052e edits touch `src/FastMultipole.jl` and
  `src/fmm.jl`, which the 030 payload also modifies; (3) live FLOWPanel has an
  untracked `BRAINSTORM/031_quadrupole_panel_farfield.md` byte-identical to the
  dev copy (benign, delete or leave at merge time). Stopped and reported per
  prompt. Stage 2 remains: fast-forward both live branches once jobs finish and
  052e state is committed/stashed by its owner, then smoke (driver +
  `runtests_unit_solver.jl`) in the live checkouts. Pending items (a) tree.jl
  guard, (b) notebook entry, (c) wt030 worktree deletion: untouched, Ryan's call.

## Implementation prompt (Phase 1 + 2)

See `BRAINSTORM/030_implementation_prompt_20260909.md`.
