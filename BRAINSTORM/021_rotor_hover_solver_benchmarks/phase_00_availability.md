# Phase 0 — Availability: every roster config runnable + instrumentation

**Status:** TECHNICALLY COMPLETE 2026-08-13 (W1–W6; W6 delivered and smoke-accepted —
see the 2026-08-13 Log entry; the 2026-08-12 corrected tables supersede 2026-08-11's).
**Approval:** [ ] (pending clear-context review + Ryan's re-approval; the 2026-08-12
approval was rescinded to add W6).

### W6 — Sparse near-field ILU preconditioner (added 2026-08-12, Ryan; DELIVERED 2026-08-13)

Develop the ILU preconditioner now, in Phase 0 (pulled forward from Phase 2b), adding
roster config **(g) GMRES + near-field ILU** to the smoke harness. Delivered by Ryan's
direct implementation (Barba direct-list pattern + `ILUZero.jl`); the staged execution
plan [`phase_00_ilu_plan.md`](phase_00_ilu_plan.md) is **superseded** (retained as design
context). Evidence: Log 2026-08-13 below.

## Scope

Make all six roster configurations (control doc ruling 1) runnable end-to-end in the rotor
pipeline, and build the instrumentation layer every later phase depends on. This is the only
phase authorized to touch `src/` before Phase 2b, and every `src/` change must land with tests.

## Work items

### W1 — ~~Unblock the Kutta closure~~ → RESCOPED (Ryan 2026-08-11): Backslash bug fixes

Diagnosis found the `src/FLOWPanel_kutta.jl:495-503` hard error only gates the opt-in 015
pressure-continuity closure — the rotor pipeline uses the legacy linear Kutta condition, which
already works with every roster solver. Ryan dropped the generalization (upgrade only if
accuracy/convergence demand it; archived design in `log.md` 2026-08-11). W1 is now the two
bugs the diagnosis surfaced:

- **`BackslashCoupled` identity-LU**: ctor's dummy `lu!` (`FLOWPanel_solver.jl:989`) is never
  replaced inside `simulate!` (`solve!` default `update_G=false`) — config 1a silently solves
  against the identity. Fix: build on first solve (`built::Bool` guard). Test-gated.
- **`Backslash` stale-`Glu`**: `_solve!(update_G=true)` factorizes into a local LU, leaving
  `solver.Glu` with new factors + old pivots (aliasing corruption) for direct consumers
  (Kutta Route A, GreenReconstruction, DirectWakePotential). Fix: `mutable struct` + write-back
  `solver.Glu = lu!(solver.G)`. Test-gated.

### W2 — FGMRES + FGS preconditioner (roster config 1f)

- Route `method=:fgmres` through the existing `method::Symbol` mechanism (Krylov.jl supports it;
  flexible GMRES tolerates a preconditioner that changes between iterations/solves).
- Preconditioner-apply wrapper: a fixed small number of FGS sweeps (`inner_iterations`-style
  knob) exposed as an `ldiv!`-style operator over the existing `FGSSolver`/`FastGaussSeidel`
  machinery. Reuse the constructed FGS tree/leaf-LUs; do not rebuild per apply.
- Refactor `KrylovSolver` (immutable struct, `src/FLOWPanel_solver.jl:296-333`) so the
  workspace and preconditioner persist across solves and can be swapped per timestep. Persisting
  the workspace also kills the per-solve allocation (:430) — a pre-registered Phase 2b candidate
  that lands here if it falls out naturally.

### W3 — Krylov warmstart plumbing (`x0`)

Support passing an initial guess via Krylov.jl's `krylov_solve!(ws, A, b, x0)` form (or
workspace warm start). Must be inserted **after** `set_strengths!` (:1083) zeroes
`body.strength`, and after the Dirichlet potential save/restore (:236-255). Off by default;
Phase 3 flips it on. Test: warmstarted solve reproduces the cold solution to solver tolerance.

### W4 — Instrumentation layer

- **Shared true-residual evaluator** (ruling 3): RMS and max of $r = Gx - b$, identical code
  path for all solvers; matrix-free solvers evaluate $r$ via one direct influence pass.
- **Convergence-history capture** with per-iteration wall-clock timestamps (ruling 4): Krylov
  `workspace.stats` (currently `@show`n and dropped at :436/:474), FGS
  `solver.fgs.residual_vector` (currently never read) + a timestamped callback if FastMultipole
  allows, block-GS outer-iteration history (currently `verbose`-print only, :892-894).
- Restore/expose `BackslashCoupled` `t_build`/`t_solve` (return commented out at :1175).
- **CSV schema** (documented in `decision_rules.md`): one row per solve with solver config,
  mesh, phase, threading mode (Julia + BLAS thread counts), setup-time components (assembly,
  factorization, tree, preconditioner), solve time, iteration count, RMS/max true residual,
  memory, commit hash, hardware tag; convergence histories in a sidecar CSV keyed by run id.
- Harness startup asserts `Threads.nthreads()` / `BLAS.get_num_threads()` and writes the log
  banner (ruling 6).

### W5 — Smoke harness

`benchmark/` driver that runs all six configs on one small `dji9443` mesh (both threading modes,
local scale is fine for smoke) and emits schema-valid CSVs.

## Exit criteria

1. All six roster configs run a `rotor_hover.jl`-style setup end-to-end on one small mesh.
2. Instrumentation CSVs validate against the schema; residual evaluator agrees with a dense
   hand computation on a tiny case.
3. New/changed `src/` code covered by tests; full test suite green.

## Log

(dated entries appended newest-last)

### 2026-08-11 — Phase 0 implemented; single-mode smoke complete; multi-mode in flight

**Code landed (all test-gated, full suite green — 45 testsets, 0 failures):**

- `src/FLOWPanel_instrumentation.jl` (new): shared operator apply
  `apply_G!`/`_apply_neumann_G!`/`_apply_dirichlet_G!` (the same code path the
  matrix-free solvers use), `true_residual!` (ruling 3; verified vs dense
  `G*x−b` to 1e-12 both BC types), `assemble_rhs!`, `ConvergenceHistory`
  (iter, `time_ns()` stamp, internal residual labeled by metric).
- `KrylovSolver` → mutable with persistent operator/rhs/workspace (kills the
  per-solve workspace+RHS allocation), `method=:fgmres` supported, `memory`
  exposed, preconditioner routing: Jacobi → `M=, ldiv=true` (left, status
  quo); new `FGSPreconditioner` → `N=, ldiv=true` (flexible right). W3
  warmstart (`warmstart::Bool=false`, solver-owned `x_prev`); `record_history`
  + per-iteration timestamps; `niter`/`solved` exposed. `@show stats` removed.
- `FGSPreconditioner` (config 1f): fixed-sweep FGS apply, zero-start (linear
  map — verified by a linearity test to 1e-10), side-effect-free, packs the
  input through the body's BC channel.
- `KrylovCoupled`: persistent workspace, stats captured, warmstart via its
  `x` vector.
- `BackslashCoupled`: **identity-LU availability bug fixed** (builds on first
  solve; previously simulate! solved config 1a against the ctor's dummy
  identity); `t_build`/`t_solve` fields live.
- `Backslash`: **stale-Glu corruption fixed** (mutable + `Glu` write-back under
  `update_G=true`; regression test on the factors/ipiv pair).
- FGS history per D2: non-breaking `callback=nothing` kwarg in FastMultipole
  `solve!` (dev checkout, `src/solve.jl`), fired once per outer iteration with
  the loop's own residual (max-abs despite the MSE label); FLOWPanel
  `_solve_history!(body, ::FGSSolver, hist)`; block-GS tuple `solve!` gained
  `history=` kwarg. FastMultipole's own solve tests pass (3166+20000).
- **New availability fix found by the smoke run**: FastMultipole
  `strength_to_value` had no overload for `RigidWakeBody{<:Any,1}` — FGSSolver
  (and FGSPreconditioner) had *never* worked on the uncapped Neumann rotor
  body type. One-line overload in `src/FLOWPanel_liftingbody.jl`, test-gated
  (FGS ≡ Backslash to 1e-10 on the plate-vortex fixture).
- Metadata: `_solver_metadata_dict` now records Krylov preconditioner
  type+knobs, `memory`, `warmstart`; `KrylovCoupled` branch added.
- Tests: `test/runtests_unit_instrumentation.jl` (35),
  `test/runtests_unit_solver_history.jl` (15), solver suite extended to 122.
- Harness: `benchmark/common.jl` (banner asserts per ruling 6, decision_rules
  CSV schema writer + validator, `min_of_k` via `time_ns()`, `summarysize` /
  `@allocated` memory) and `benchmark/rotor_hover_solver_smoke.jl` (all six
  configs; smoke-scale Krylov `itmax=400, rtol=1e-4` — Phase 1 calibrates the
  real frozen tolerances).

**Single-mode smoke (6240-panel `dji9443_20260722_40_41_uncapped`, 1 Julia
thread, BLAS pinned to 1, local laptop — NOT publishable numbers):** all six
configs ran end-to-end (isolated frozen-RHS solve → history → shared true
residual → 3-step `simulate!`), `runs.csv` schema-validated (6 rows) +
history sidecars in `benchmark/results/smoke/single/`.

| config | t_solve_min [s] | iters | rms true residual | notes |
| --- | --- | --- | --- | --- |
| backslash_coupled | 2.74 | — | 5.8e-12 | t_build 4.10 s (assembly+lu) |
| backslash_iterative | 0.0125 | 1 | 5.9e-12 | assembly 3.39 s + lu 0.95 s; block-GS runs 2 identical inner solves |
| krylov_gmres | 224.9 | 290 | 2.6e-2 | ill-conditioned; 290 iters to the LOOSE smoke rtol=1e-4 |
| krylov_jacobi | 280.9 | 362 | 4.8e-1 | cell_size=R/4 *hurt*: more iters AND worse true residual |
| fgs | 0.047 | 1 | 9.8e-7 | converged in ONE outer iteration |
| fgmres_fgs | 22.4 | 26 | 1.3e-2 | 11× fewer iters and 10× less time than plain GMRES at matched rtol |

Headline smoke findings (to be confirmed at scale in Phases 1–2):
1. The FGS leaf-sweep machinery is decisively superior on this operator:
   `FGSSolver` converges in one outer iteration (residual 9.8e-7) where
   unpreconditioned GMRES needs 290 iterations for rtol=1e-4.
2. Block-Jacobi (R/4 cells) *worsens* the true residual at equal iteration
   count — preconditioner choice is load-bearing, not decorative.
3. Krylov matrix-free applies rebuild the FMM tree every iteration
   (~4.8–6.1 GB allocated per solve) — pre-registered Phase 2b optimization.
4. The two Backslash configs agree at round-off residual but their 3-step CT
   differs (−2.8e-5 vs −4.0e-5, tiny absolute transients) — Phase 1's
   agreement table owns this.
5. `iterations` CSV column initially recorded the post-simulate warmstarted
   solve; fixed to capture the benchmark solve (single-mode CSV predates the
   fix for configs c/d/f — their `iterations` cells read 0/2; histories are
   correct).

**Multi-mode smoke (6 Julia threads, BLAS 4 recorded — local, not
publishable):** all six configs end-to-end, 6 schema-valid rows +
sidecars in `benchmark/results/smoke/multi/`. Iteration counts identical to
single mode (right: same math, different threading); residuals match to
round-off.

| config | t_solve_min [s] (multi) | vs single | notes |
| --- | --- | --- | --- |
| backslash_coupled | 0.646 | 4.2× | t_build 1.71 s (2.4×) |
| backslash_iterative | 0.0119 | 1.05× | ldiv! BLAS-bound, already fast |
| krylov_gmres | 60.6 | 3.7× | FMM apply threads well |
| krylov_jacobi | 76.0 | 3.7× | |
| fgs | 0.0272 | 1.7× | but its 3-step simulate! ran 1070 s vs 145 s single — investigate in Phase 2 |
| fgmres_fgs | 6.83 | 3.3× | |

Known warts (all noted, none blocking):
- The `iterations` cells for configs c/d/f in both modes were initially the
  post-simulate warmstarted counts (script bug, now fixed in
  `rotor_hover_solver_smoke.jl`); the CSVs were patched from the history
  sidecars (290 / 362 / 26), which were always correct.
- `fgs` multi-mode `simulate!` slowdown (1070 s vs 145 s single-thread) is
  unexplained — flagged for Phase 2 bottleneck attribution.
- Smoke Krylov tolerances are LOOSE (`rtol=1e-4`, `itmax=400`) so the run
  finishes locally; the true residuals for c/d/f in the CSVs reflect that,
  not solver quality at production tolerance.

**Exit criteria review:**
1. All six roster configs run a `rotor_hover.jl`-style setup end-to-end on
   one small mesh — **MET** (both threading modes).
2. Instrumentation CSVs validate against the schema; residual evaluator
   agrees with a dense hand computation on a tiny case — **MET**
   (`validate_runs_csv` passes; evaluator ≡ dense `G*x−b` to 1e-12 in
   `runtests_unit_instrumentation.jl`).
3. New/changed `src/` code covered by tests; full test suite green — **MET**
   (45 testsets, 0 failures, incl. 35 instrumentation + 15 history tests and
   the extended solver suite).

Phase 0 is technically complete pending clear-context review + Ryan's gate
approval.

### 2026-08-12 — Review remediation: cold-state rerun REPLACES the tables above

A clear-context review (5 confirmed findings; assessment + fixes in `log.md` 2026-08-12)
invalidated parts of the 2026-08-11 evidence. **The 2026-08-11 tables above are superseded
by the tables below** — kept per append-only policy; do not quote them. What was wrong:
the timing loop never reset `body.strength`, so the standalone-FGS "1 iteration / 0.047 s"
and block-GS "1 outer" numbers were warm-state no-ops (FGS seeds from current strengths and
breaks before sweeping); the Jacobi row predates ruling 11a right-routing. Configs a/c/f
were genuinely cold and reproduce below to ~1% — the cross-run consistency check.

Remediation landed: untimed cold reset before every solve + FGS tripwire; provenance
columns (`fm_commit`, `julia_version`, `solver_settings`, `backend_settings`), `t_rhs`,
`banner.txt`; comparable `solver_state_bytes` (summarysize minus referenced bodies);
FastMultipole callback committed (`5adde3b`) + portable patch under `patches/`.

**Corrected cold single-thread table** (6240 panels, BLAS=1, smoke tolerances
`rtol=1e-4`/`itmax=400` for Krylov, FGS `tolerance=1e-6`; local, not publishable):

| config | t_solve_min [s] | iters | rms true residual | mem_state [MB] |
| --- | --- | --- | --- | --- |
| backslash_coupled | 2.63 | — | 5.8e-12 | 311.8 |
| backslash_iterative | 0.023 | 2 | 5.9e-12 | 311.8 |
| krylov_gmres | 225.2 | 290 | 2.6e-2 | 17.7 |
| krylov_jacobi (right-routed) | 310.8 | 400 (cap) | 4.9e-2 | 22.4 |
| fgs | 35.2 | 204 | 9.8e-7 | 110.6 |
| fgmres_fgs | 22.1 | 26 | 1.3e-2 | 133.2 |

**Corrected multi-thread table** (6 Julia threads, BLAS 4, clean machine — the first multi
rerun was discarded: Ryan had a concurrent job; this one ran quiet):

| config | t_solve_min [s] | vs single | iters | note |
| --- | --- | --- | --- | --- |
| backslash_coupled | 0.61 | 4.3× | — | |
| backslash_iterative | 0.024 | 1.0× | 2 | ldiv! BLAS-bound |
| krylov_gmres | 55.8 | 4.0× | 290 | |
| krylov_jacobi | 77.4 | 4.0× | 400 (cap) | |
| fgs | 77.2 | **0.46×** | **267** | slower AND more iterations threaded |
| fgmres_fgs | 6.40 | 3.4× | 26 | |

**Corrected headline findings:**

1. Cold standalone FGS costs **35.2 s / 204 outer iterations** to reach 9.8e-7 (max-abs) —
   not "1 iteration / 0.047 s". It still beats plain GMRES decisively on time-to-residual
   (~6× faster to a ~4-decade tighter residual), but the e-vs-f ranking is now genuinely
   open: 1f reached only 1.3e-2 in 22.1 s (loose smoke rtol). Phase 1's matched-residual
   calibration settles it.
2. **"Jacobi (R/4) hurts" survives honest re-derivation**: right-routed, it hits the
   400-iteration cap where plain GMRES converges at 290, with a worse true residual
   (4.9e-2 vs 2.6e-2) and more time. The preconditioner exploration (ruling 11d, ILU
   priority) stands motivated.
3. **FGS threads badly — now isolated to the solve, not just simulate!**: multi-thread FGS
   is 2.2× *slower* than single (77.2 vs 35.2 s) and needs more iterations (267 vs 204 —
   threaded leaf sweeps change the Gauss–Seidel update order, degrading convergence toward
   Jacobi-like behavior; residual endpoint shifts in the 3rd digit, consistent with
   order-dependence). Pre-registered for Phase 2 bottleneck attribution; affects 1e AND
   every 1f apply.
4. Cross-run reproducibility: untouched configs (a, c, f) match the 2026-08-11 run to ~1%
   single-mode; residuals bit-identical. The harness is deterministic enough to trust at
   smoke scale.

Exit-criteria review unchanged (all three MET — the criteria concern availability, schema
validity, and tests, all of which the rerun re-confirms with the corrected protocol).
Status returns to: technically complete, pending clear-context review + Ryan's gate
approval.

### 2026-08-13 — W6 delivered: FMM-direct-list ILU(0) is the decisive smoke winner

Ryan implemented W6 directly on 2026-08-12 evening (collapsing the staged plan's Stages
A–C; `phase_00_ilu_plan.md` is superseded as an execution plan, retained as context). The
two open design decisions were resolved by the implementation itself: **pattern source =
FastMultipole Barba direct interaction list** (a dedicated tree build at ctor time — the
staged plan's radius-cutoff fallback was skipped entirely) and **factorization =
`ILUZero.jl`** (new registered dep in `Project.toml`, compat 0.2.0).

What landed (all in the W1–W5 uncommitted stack):

- `src/FLOWPanel_solver.jl` § "FMM-DIRECT-LIST ILU(0) PRECONDITIONER":
  `ILUPreconditioner(body; leaf_size=10, multipole_acceptance=0.4, max_pattern_entries=
  512N, equilibrate=false, diagonal_shift=0.0, keep_matrix=false)`. Sparse near-field
  operator assembled from expanded Barba direct-list panel pairs in tree order via the
  same `induced` kernel path `_G!` uses (side-aware self limits, attached-wake terms,
  operator mode, `kerneloffset_panel`); linear pattern-size guard fires before any sparse
  allocation or kernel evaluation; missing diagonals topped up; ILU(0) via `ILUZero.ilu0`
  with pivot checks; `ldiv!` = two sparse triangular solves with body↔tree permutation
  and optional row equilibration. Right-routed automatically by `_krylov_launch!`'s
  existing else-branch (ruling 11a, no Krylov-side changes).
- `src/FLOWPanel_metadata.jl`: `ILUPreconditioner` branch → `solver_settings` CSV column.
- Tests (`test/runtests_unit_solver.jl`): dense-pattern exactness (single leaf ⇒ ILU(0) ≡
  exact LU ≡ `G\b`, incl. the wake-aware plate-vortex fixture), pattern ≡ `_G!` entries on
  a nontrivially-ordered sphere, guard behavior + state restoration, equilibration
  algebra, ldiv! linearity, Krylov end-to-end vs Backslash, metadata keys, construction
  ladder linearity, direct-list ≡ `fmm!`'s lists.
- Harness: config (g) `krylov_ilu` + end-of-run ILU-vs-GMRES break-even summary.

Session review of the implementation (this agent, 2026-08-12/13): verdict CORRECT —
assembly semantics match `_G!` exactly; permutation/equilibration algebra verified;
ILUZero-internal assumptions (U diagonal stored last per column; 3-arg `ldiv!`; `nnz`)
audited against ILUZero's source. Two fixes applied to `src` (test-gated): (1) ctor now
calls the shared `_ilu_direct_pattern` helper instead of duplicating its ~45 lines of
tree/list/guard logic (helper extended to return tree/list timings); (2) the assembly
kernel loop is now `Threads.@threads` over direct-list pairs with precomputed triplet
offsets, mirroring `_G!` (`induced` is pure). Two crash bugs fixed in the harness's new
uncommitted code: escaped quotes inside `$()` string interpolation (Julia 1.12 parse
error), and an empty-generator `sum` in `run_config!`'s `setup_total` fallback for
configs with no setup timings (`init=0.0`). Unit solver suite + full suite green.

**Smoke results** (7 configs, 6240 panels, smoke tolerances `rtol=1e-4`/`itmax=400`;
full provenance in `benchmark/results/smoke/*/runs.csv`; the harness truncates
`runs.csv` per launch). **CSVs of record = the 2026-08-13 metadata-remediation rerun**
(first gate review's finding 1: `solver_settings` was captured after the `simulate!`
leg's `warmstart=true` flip; the harness now captures it before that leg, and all
Krylov-family rows correctly record `warmstart=false` with the sim-leg flip in `notes`).
Rerun caveats, per Ryan: **the two modes ran concurrently** (5 threads total), so wall
times carry some contention — iteration counts and residuals are deterministic and
reproduced exactly; and **multi mode ran at 4 Julia threads** (prior 2026-08-12 tables:
6 threads, quiet machine), so vs-single ratios are not comparable across reruns. The
first W6 run's quiet-machine numbers (2026-08-12/13: ILU 11.51 s single / 2.86 s multi
@6t; GMRES 223.4 / 55.2) remain quotable from this entry's history and reproduce within
contention noise.

Single-thread (BLAS=1, concurrent with the 4-thread multi run):

| config | t_solve_min [s] | iters | rms true residual | t_precond [s] | mem_state [MB] |
| --- | --- | --- | --- | --- | --- |
| backslash_coupled | 2.59 | — | 5.8e-12 | — | 311.8 |
| backslash_iterative | 0.023 | 2 | 5.9e-12 | — | 311.8 |
| krylov_gmres | 281.4 | 290 | 2.6e-2 | — | 17.7 |
| krylov_jacobi | 303.0 | 400 (cap) | 4.9e-2 | 0.57 | 22.4 |
| fgs | 34.8 | 204 | 9.8e-7 | — | 110.6 |
| fgmres_fgs | 22.0 | 26 | 1.3e-2 | 2.99 | 133.2 |
| **krylov_ilu** | **11.47** | **15** | 1.8e-2 | **1.14** | 29.9 |

Multi-thread (4 Julia threads per Ryan, concurrent with the single run):

| config | t_solve_min [s] | vs single | iters | note |
| --- | --- | --- | --- | --- |
| backslash_coupled | 0.80 | 3.2× | — | |
| backslash_iterative | 0.024 | 1.0× | 2 | ldiv! BLAS-bound |
| krylov_gmres | 80.3 | 3.5× | 290 | |
| krylov_jacobi | 101.1 | 3.0× | 400 (cap) | |
| fgs | 33.8 | 1.0× | 423 | thread-order nondeterminism again (423 here; 500 @6t 08-13; 267 @6t 08-12; single stable at 204) — pre-registered for Phase 2 |
| fgmres_fgs | 8.15 | 2.7× | 26 | |
| **krylov_ilu** | **3.88** | 3.0× | **15** | ILU: nnz 1.15M (184/panel, 2.9% dense), factor nnz = S nnz |

**W6 acceptance: PASS, both modes** (and in both the quiet first run and the
remediation rerun). At matched requested tolerance ILU-GMRES beats plain GMRES on
iterations (290 → 15) and wall time (281.4 → 11.47 s single; 80.3 → 3.88 s multi;
quiet-machine first run: 223.4 → 11.51 / 55.2 → 2.86 @6t), with total setup ~1.0–1.2 s
(trees 0.10 s + interaction lists 0.04 s, assembly ~0.42 s single, ILU(0) factorization
~0.25 s, remainder ~0.4 s untracked ctor overhead — state copies, sparse build,
summarysize) — amortization break-even 0.004 solves (single) / 0.014 (multi).
It also displaces FGMRES+FGS (26 iters / 22.0 s / setup 2.99 s) as the best iterative
config at smoke scale, at a quarter of FGS's memory (29.9 vs 133.2 MB) and ~3× fewer
allocations per solve. True-residual endpoint at the loose smoke rtol is comparable to
plain GMRES (1.8e-2 vs 2.6e-2); matched-residual calibration is Phase 1's job.

Exit criteria: all three still MET, now including config (g) end-to-end with schema-valid
rows and populated `t_precond`/`solver_settings`. Status: W6 technically complete;
Phase 0 pending clear-context review + Ryan's gate re-approval.
