# FGS multithreading: determinism + performance investigation (COMPLETE)

**Status:** COMPLETE 2026-08-17. Root cause was a threaded M2L owner overlap
caused by a deterministic but non-owner-major interaction list. FastMultipole
now canonicalizes `(target, source)` before threaded assignment and dependent
matrix/index construction. The full four-way matrix is exact within every
configuration; Julia 4T/BLAS 1T fixed-ten timing changed +0.87%, passing the
5% gate. All 18 R1--R3 post-crossing margins pass certified BC error. Results:
`results/fgs_determinism/summary.md`.

This is FastMultipole work
(Ryan's dev checkout `../FastMultipole`, currently `5adde3b-dirty` with three
uncommitted campaign fixes — coordinate any new edits with that pending
commit; never touch `scripts/20250404_prediction_accuracy.jl`).

## Motivating observations (evidence of record)

- Phase 0 (smoke, pre-registered): identical FGS config, 4 threads — 267
  outer iterations one run, 500 the next; single-threaded runs exactly
  reproducible (204 every time).
- Phase 1 Stage 2 (fgstune_verify.csv R2, 2026-08-15): staircase run crossed
  BC 1e-6 at 7 outer iterations; fresh verification solve in the SAME process
  (new FGSSolver instance, stopping tol = the crossing residual, zero margin)
  needed 15 iterations / ~2× wall time — trajectories genuinely diverge, and
  zero-margin threshold stopping amplifies the divergence into large cost
  scatter near the convergence floor.
- Consequence for the campaign: local FGS timings carry scatter; HPC
  min-of-k bounds it but doesn't remove it; determinism would also make the
  τ-ladder calibration (stopping tolerances) transferable between runs.

## What the code already guarantees (verified by reading, 2026-08-15 —
## saves the executing agent the survey)

All threaded loops use `:static` scheduling with deterministically
constructed assignments:

- `solve.jl` (the FGS outer/inner loops, leaf `mat \ rhs` solves,
  `update_nonself_influence!`, `residual!`): **NO threading at the Julia
  level** — the Gauss–Seidel sweep is serial by leaf in tree order. Only
  BLAS/LAPACK inside `mat \ rhs` is (potentially) multithreaded.
- In-loop far-field `fmm!` (one per outer iteration) passes an **empty
  direct list** — near-field runs through the serial self/nonself matrices,
  so `nearfield_multithread!` (fmm.jl:153, with its commented-out data-race
  debug scaffolding) is NOT in the FGS loop path.
- `assign_m2l!` (fmm.jl:518) and `make_direct_assignments!` (fmm.jl:66) both
  finish the current TARGET before cutting an assignment → no two threads
  write the same target expansion; assignment boundaries are deterministic
  functions of the sorted list.
- Threaded interaction-list build (interaction_list.jl:39) concatenates
  thread-local lists in ASSIGNMENT-INDEX order (interaction_list.jl:56) —
  deterministic given the tree and `start_list`.

NOT yet verified (first targets for the experiments below):

- `assign_m2m!` / the upward- and downward-pass task splits
  (fmm.jl:~440–470, :725, :780): whether a PARENT's child-accumulations can
  split across threads (concurrent `+=` into one parent expansion would be a
  true race, not just reduction-order noise).
- Threaded tree construction (tree.jl, 7 sites): whether `sort_index` /
  branch ordering is bit-reproducible across builds in one process.
- `sort_by_target` stability, and BLAS determinism at 4 threads for the leaf
  `getrf`/`getrs` sizes involved (leaf 100–300).

## Part A — localize the nondeterminism (cheap, do first)

Fixture: the R1 campaign case (phase1_case.jl), shared FGS knob set
(p=6/MAC 0.3/leaf 150/inner 5), 4 threads. Instrument with the existing
`_solve_history!` callback; compare full residual SEQUENCES and iterate
snapshots (the fgstune mid-solve `buffer_to_system_strength!` pattern)
bit-for-bit (`===` on bits via `reinterpret(UInt64, x)`).

- **E1 (in-loop?):** ONE solver object, N repeated cold solves. Bit-identical
  sequences ⇒ the loop (incl. its fmm! passes + BLAS) is deterministic and
  the variation is ctor-time; differing ⇒ in-loop source (go to E4/E5).
- **E2 (ctor-time?):** two FGSSolver instances in one process; diff their
  trees (`sort_index`, `bodies_index`, branch centers/radii), `m2l_list`,
  `direct_list`/`full_direct_list`, and each self-matrix bit-for-bit. Any
  diff localizes to tree/list construction.
- **E3 (thread-count matrix):** {1T, 4T} × {BLAS 1, BLAS 4} × 3 repeats.
  Separates Julia-thread from BLAS-thread effects (Phase 0 says 1T is
  reproducible — confirm BLAS was pinned there).
- **E4 (M2M/L2L race check):** read `assign_m2m!` + downward equivalents for
  parent-boundary splitting; if suspicious, empirically hammer: repeated
  upward passes on fixed strengths, diff expansions bitwise.
- **E5 (BLAS):** repeated `lu(mat)\rhs` on one captured leaf block at 4 BLAS
  threads, diff bitwise.

Deliverable: a table attributing the variation to {ctor, upward/downward,
M2L, BLAS, other} + a log entry. (Estimated: small — a scratchpad script and
a few hours.)

## Part B — determinism fixes (gated on Part A's attribution)

In likely-order of attribution:

1. **IMPLEMENTED — canonicalize interaction ordering at ctor:** stable
   O(list + branches) FastMultipole counting sorts by source then target give
   lexicographic `(target, source)` order before threaded M2L assignments,
   matrices, and index maps. This fixes the demonstrated target-owner overlap.
2. **Serialize-or-partition parent accumulation in M2M/L2L** (if E4 hits):
   make assignments respect parent boundaries (mirror the assign_m2l!
   pattern) — removes the race AND fixes the order.
3. **Pin BLAS threads for leaf solves** (if E5 hits): `BLAS.set_num_threads(1)`
   scoped around the sweep (leaf blocks ≤300² are small enough that BLAS
   threading likely doesn't pay anyway — measure; this may be a performance
   WIN too).
4. **IMPLEMENTED — correctly directed stopping margin:** find the first
   finite residual *after* the BC crossing that is strictly lower than the
   crossing residual, then use their geometric mean. If no such successor
   exists, extend the staircase; R3 tau 1e-6 exercised this guard. Never use
   the unsafe pre-crossing mean or a fallback.

Each fix lands in FastMultipole with its own unit test (mirror the
shared-topology-fix precedent: premise-guarded, non-vacuous), then: rerun
E1–E3 to demonstrate bit-reproducibility, rerun fgstune verify rows R1–R3,
and quantify any speed delta (determinism must not cost >~5% or it needs a
Ryan ruling).

## Part C — performance improvements (concrete candidates, measure-first)

Profile one outer iteration at R2/R3 shared knobs FIRST (far-field fmm! vs
leaf solves vs `update_nonself_influence!` vs residual) — then attack in
measured-cost order. Candidates, biggest expected first:

1. **Cache leaf LU factorizations** (known deferred item, decision D3 /
   FGSPreconditioner docstring): `mat \ rhs` re-factorizes every dense leaf
   block on EVERY sweep and inner iteration — O(L³) each — while the blocks
   are geometry-fixed at ctor. Factor once at ctor (or first solve), reuse
   `getrs` (O(L²)) thereafter. At leaf 150–300 and inner 5–20 this
   plausibly changes the tuner's whole leaf/inner tradeoff — re-run the
   τ-ladder after landing it.
2. **Thread the leaf sweep deterministically**: strict GS is sequential by
   leaf, but a fixed graph coloring of leaves by direct-list adjacency
   (computed at ctor) allows parallel same-color solves with deterministic
   color order. CHANGES the iteration (colored GS ≠ lexicographic GS) —
   convergence must be re-measured; needs Ryan's sign-off as a solver-variant
   decision, not a drop-in.
3. **Batch `update_nonself_influence!`**: currently invoked after every leaf
   solve (O(direct-list touch each time)); accumulate per color/segment or
   restructure to touch each direct-list entry once per sweep.
4. **BLAS-thread audit for leaf sizes** (see B3): possibly faster single-
   threaded at these block sizes; frees Julia threads for anything colored.

Constraint: items 2–3 alter numerical behavior → separate commits, each with
convergence A/B on the R1–R3 ladder (fgstune staircases make this cheap);
item 1 is bit-neutral in exact arithmetic but WILL perturb rounding — treat
as a determinism-relevant change (rerun Part A E1 after).

## Ordering & gates

Part A → (B fixes gated on attribution; B4 immediately) → Part C item 1 →
re-run τ-ladder R1–R3 → C items 2–4 only with Ryan's ruling. Everything
lands test-gated in FastMultipole; campaign CSVs regenerate only where knobs
or numerics changed. Local runs ≤4 threads; HPC for anything >20 min.
