# 030 — context-reset prompt for Phase 3 (route A) + Phase 3b prototype (2026-09-10)

SUPERSEDES `030_reset_prompt_20260910b.md` (state snapshot there is still accurate;
this file adds Ryan's Phase 3 decisions and the Phase 3b staging). Read FIRST, in
order: `BRAINSTORM/030_generic_influence_block_assembly.md` (design of record —
the Log carries all session entries including the 2026-09-10 rulings),
`030_reset_prompt_20260910b.md` (complete Phases 1+2 state: commits, test totals,
benchmark of record, gotchas — do not re-derive anything listed there), then this
file. Also read `agent_policies/WORKFLOW.md` and `agent_policies/TESTING.md`
before touching code.

Ground rules (unchanged): work ONLY in the `~/wt030` worktrees
(`~/wt030/FastMultipole`, `~/wt030/FLOWPanel.jl`, branch `030-block-assembly`);
never touch the live checkouts under `~/Dropbox/research/projects/` except the
single allowed append to the 030 item's Log at the end. Local only, max 4
threads, commit but do not push or merge. Sign conventions are FLOWPanel's #1
regression source; rtol-1e-12 exactness tests are the guard.

## Ryan authorizations (2026-09-10)

1. **Phase 3 IS authorized, route A** (project-after-assembly). Decision
   rationale, recorded: reuse the exact raw hook; keep `influence!` the single
   owner of projection semantics (avoids duplicating sign-sensitive logic); the
   raw-block scratch is L2-resident (~50 KB at leaf 40) so the round-trip should
   be noise against the compute-bound kernel. Fall back to a projected-hook
   variant ONLY if the cheap-kernel benchmark shows the scratch round-trip
   costing >~10% — measure, don't assume.
2. **Phase 3b IS authorized as a prototype**: far-pair point-panel approximation
   (details below).
3. Merge-back into `fastmultipole`/main FLOWPanel branches remains GATED — do
   not merge. Item approval checkboxes remain Ryan's.

## Phase 3 — FGS migration, route A

The FGS builders `self_influence_matrices` (`FastMultipole/src/solve.jl:408`,
probe loop ~:483-497) and `nonself_influence_matrices` (`:124` region) currently
probe: unit strength → `direct!` one column into the target buffer →
`influence!` collapses to ONE scalar per target (FLOWPanel Neumann:
`dot(gradient, normal)`, `src/FLOWPanel_abstractbody.jl:1427`) → copy the
column. FGS blocks are |targets|×|sources| POST-projection scalars, unlike the
cache's raw `n_out`·|targets| rows.

Migration shape:
- Per block: assemble the RAW block via `assemble_influence_block!` into a
  per-worker scratch (reuse one scratch sized to the largest block, not per-block
  allocs), then collapse rows through the projection into the FGS matrix column
  space. Adapter note: a raw block column reshaped `(n_out, n_targets)` has the
  same layout as the target buffer's output rows, so the projection can be
  applied by either (a) copying the column into the target-buffer output rows
  and calling stock `influence!`, or (b) a thin `influence!`-equivalent that
  reads the reshaped scratch directly — prefer whichever keeps `influence!` the
  single source of truth; if you write (b), it must CALL the same primitives
  (`get_gradient`-equivalent row reads + `get_normal`), not re-derive the math.
- Route each block through `overrides_block_assembly` exactly like the cache
  builder: opted-in systems skip the probe; fallback systems keep the existing
  probe loop verbatim (do not break non-opted-in users).
- Preserve: `Matrices` packed storage and `get_matrix_vector` consumers
  (`build_leaf_lu_cache` LU-factors self matrices — blocks must stay square and
  dense), strength save/restore, deterministic results.
- Tests: FGS matrices assembled vs probed at rtol 1e-12 on the gravitational
  system AND a FLOWPanel body with shedding (the FGSSolver/FGSPreconditioner
  testsets in `test/runtests_unit_solver.jl` are the harness — add an
  assembled-vs-probed testset mirroring the 030 Phase 2 one); existing FGS
  testsets must stay green. Benchmark: FGS matrix build hook-vs-probe on the
  standard grav case AND on a FLOWPanel fixture; report the route-A scratch
  overhead explicitly (the >~10% fallback trigger).

## Phase 3b — far-pair point-panel approximation (PROTOTYPE)

Motivation (Ryan 2026-09-10): for FLOWPanel bodies the ~µs-per-pair panel
integral IS the build cost (>99%), and every direct-list pair pays it regardless
of separation — core_size radius inflation makes much of the interaction direct
(memory: the diamond fixture is all-direct at any MAC). Reconnaissance
2026-09-10: NO existing point-panel approximation in active kernels — `induced`
(`src/FLOWPanel_elements_fmm.jl:243-322`) always runs the full Hess-Smith
integral; only commented-out legacy code in `FLOWPanel_elements.jl`;
FastMultipole's panel→multipole (`bodytomultipole.jl`) acts at BRANCH level
beyond the MAC, never per-pair inside direct blocks.

Physics: at target distance r ≫ panel size L, a ConstantSource panel → point
source of strength σ·A at the panel centroid; ConstantDoublet → point dipole
μ·A·n̂; VortexRing ≡ constant-doublet panel → same dipole form (check the
strength-column convention against `induced`'s VortexRing method, and mind the
GeometricTools normal-orientation sign flips). Leading error O((L/r)²) — a
quadrupole correction is OPTIONAL scope if η must otherwise be large.

Prototype scope (deliberately narrow):
- Implement at the FLOWPanel assembly-hook level ONLY
  (`_assemble_influence_block!` in `src/FLOWPanel_abstractbody.jl`): per pair,
  compute r² to the panel centroid; if r² > (η·L)² use the point kernel, else
  the full `induced`. Precompute per-source-column centroid, area, n̂, and L
  (max edge or sqrt(A)) ONCE per column from the buffer vertices — not per pair.
  Do NOT touch `_direct_body!`/`direct!` or production evaluation paths — the
  prototype changes only what gets ASSEMBLED into cached/FGS blocks.
- Exact mode stays the DEFAULT. Knob: e.g. `farfield_eta::Float64 = Inf` (Inf =
  exact) threaded through the cache/FGS build entry points the same way
  `use_block_assembly` is. The rtol-1e-12 assembled-vs-probed tests run at
  Inf and MUST stay green untouched.
- η selection: sweep η ∈ {2, 3, 4, 5, 8} on the diamond (Dirichlet, shedding)
  and sphere (Neumann) fixtures; report max entrywise relative error vs the
  exact block and the (L/r)² trend. Target: error commensurate with FMM p=8
  truncation (~1e-8) — propose the η that achieves it.
- Acceptance beyond entrywise error: a full cached KrylovSolver solve at the
  proposed η vs exact — report strength-vector rel diff and boundary-condition
  residual delta (assert_boundary_residuals pattern). The operator changes, so
  gate on solve-level effect, not just per-entry error.
- Measure: per-pair time point vs full kernel; end-to-end cache/FGS build
  speedup at the proposed η on the FLOWPanel fixtures; report the fraction of
  pairs that qualified as far at each η (if <~30% qualify at the accurate η,
  say so — the phase may not pay and Ryan decides).
- This is a PROTOTYPE: production adoption (using approximated blocks in real
  campaigns, or extending to `direct!` itself) is a separate Ryan gate.

Suggested order: Phase 3 first (correctness anchored by exactness tests), then
3b on top (its speedup claim needs the migrated FGS builders to matter for FGS).
If Phase 3 route A hits the >10% scratch overhead, STOP and report before
building the projected-hook variant.

## Practical notes

- Test drivers and totals: see `030_reset_prompt_20260910b.md` (FastMultipole
  cache file 9518/9518 via the runtests-prelude driver; FLOWPanel
  `runtests_unit_solver.jl` 489/489). Re-run both at the start to confirm your
  baseline, and after every change.
- Machine: 8 cores/16 GiB; check `ps aux | grep julia` for competing agent jobs
  before benchmarking; no `timeout` binary on this Mac (background watchdog
  pattern); pipe long runs to a FILE.
- The zero-controlpoints FmmPlan hang gotcha and the `unsafe_get_block_matrix`
  GC.@preserve rule are in the b-file — read them before writing tests.
- Commit clean, separate commits per phase on `030-block-assembly` in each repo;
  list them in your report. Append ONE dated Log entry to the live item at the
  end. Final report: file:line changes, verbatim test totals, benchmark tables
  (route-A overhead %, η sweep, build speedups), deviations, open questions.
