# Phase 2 — Single-step benchmarks: setup vs per-step cost, then speedups

**Status:** IN PROGRESS — leaf-LU cache implementation and local R1 cache-only
A/B complete 2026-08-17; dedicated-node ladder and production retuning pending.
**Approval:** [ ]

## Prior-phase handoff

Phase 1 delivers: frozen mesh ladder, frozen per-solver matched-residual settings.

## Stage 2a — Baseline benchmarks

For every (roster config × ladder rung × threading mode), on HPC per rulings 5–6:

- **One-time setup**, componentized: influence-matrix assembly (`_G!`), `lu!`, FMM tree build,
  leaf-LU factorization, preconditioner build. (`Backslash` does assembly+LU in its
  constructor; time the constructor's components, not just `solve!`.)
- **Per-step cost**: isolated solve on frozen RHS, min-of-k (k ≥ 5) after warmup; RHS-assembly
  cost timed separately. **Strictly cold — warmstart off everywhere** (ruling 12; Phase 3
  owns warm measurements against these baselines).
- **Full-timestep share**: one `simulate!` timestep per config (kinematics + solve + Kutta +
  wake) to report the solve's fraction of a real step.
- Convergence history vs wall-clock time for iterative solvers (shared axes per rung); RMS/max
  true residual for every solve; memory per config; validate the equal-time-per-iteration
  assumption (ruling 4).
- Fit scaling exponents $t \propto N^p$ per solver per cost component; note any
  matrixful/matrix-free crossover in time or memory.

**Bottleneck attribution**: profile the dominant per-step cost of each solver
(`Profile`/allocation tracking on a mid-ladder rung, both threading modes) and rank candidate
levers.

## Stage 2b — Speedups (two tiers, re-benchmark after each)

1. **User-facing tuning first**: Krylov tolerances/`itmax`/restart behavior, preconditioner
   cell size, FGS knobs (`inner_iterations`, `rlx`, `expansion_order`, `leaf_size`,
   `multipole_acceptance`), FMM backend knobs. Tuning must not break Phase 1 matched-residual
   agreement — re-verify agreement at any adopted setting.
2. **Authorized source improvements** (test-gated; each lands separately with a before/after
   ledger row). Pre-registered candidates:
   - Persistent Krylov workspace (kill the per-solve `krylov_workspace` allocation,
     `src/FLOWPanel_solver.jl:430`) — LANDED in Phase 0 W2.
   - Whatever 2a profiling reveals (allocation hotspots, redundant influence passes,
     `set_strengths!` zeroing cost, RHS assembly).
   - **GMRES preconditioner exploration beyond block Jacobi** (Ryan 2026-08-12; motivated by
     the Phase 0 smoke, where Jacobi at `cell_size=R/4` *increased* iterations 290→362 AND
     worsened the true residual, while the FGS preconditioner cut iterations 11×). Candidates,
     roughly cheapest-first:
     - Jacobi done right: `cell_size` sweep (R/4 was arbitrary; try leaf-scale cells), and
       route Jacobi as a *right* preconditioner (`N=`) so the monitored residual stays the
       true one — the left-routing is why the jacobi config stopped at a 20× worse true
       residual than plain GMRES at the "same" tolerance.
     - Diagonal/row equilibration of the operator (panel areas vary root→tip on the rotor
       mesh; a symmetric scaling is nearly free and attacks that part of the conditioning).
     - **Sparse near-field ILU — PRIORITY candidate (Ryan 2026-08-12)**: assemble the
       near-neighbor influence block sparsely (pattern = FMM tree direct lists) and
       ILU(0)/ILUT it — the classic fast-BEM preconditioner. Full design basis, cost model,
       wake subtleties, and validation plan in [`theory_nearfield_ilu.md`](theory_nearfield_ilu.md);
       implement per that doc, right-routed (`N=`) per ruling 11a.
     - FGS-sweep preconditioner tuning (sweeps × inner_iterations × leaf_size) — config 1f's
       apply is 1 fmm! per sweep; the trade is apply cost vs Krylov iterations.
     - Deflation/coarse-space correction only if the above stall (research-tier; likely out
       of scope for the paper).
   - **LANDED 2026-08-17:** FGS leaf-LU caching in FastMultipole (deferred D3)
     benefits configs 1e and 1f; dedicated-node ladder A/B and retuning remain.
   - **Near-field influence-matrix cache — STAGED (Ryan 2026-08-19)**: store
     each direct-list entry's dense strength→output block in one packed
     `Matrices`-style vector (target-major, owner-parallel matvec evaluation),
     keyed to a `Tree`, general over any `direct!` kernel, reusable across
     applies/timesteps while relative geometry is frozen. Motivated by the R1
     profile (panel kernel ≈ 40%+ of a Krylov solve's samples). Full
     self-contained implementation plan for a fresh agent:
     [`nearfield_matrix_cache_plan.md`](nearfield_matrix_cache_plan.md).
   - Krylov matrix-free apply rebuilds the FMM tree every iteration (~5 GB allocated per
     290-iteration smoke solve) — cache the tree across applies within a solve (geometry is
     frozen during a solve by construction).
   - Investigate the FGS multi-thread `simulate!` anomaly (1070 s vs 145 s single-thread on
     the smoke case) before trusting any multi-mode FGS numbers.

## Exit criteria

1. Complete baseline table + convergence-vs-time figures across the ladder, both threading
   modes, in `ledger.md` + `figures/`.
2. Bottleneck attribution written up per solver.
3. Before/after tables for every adopted speedup; final re-benchmark supersedes the baseline
   as the paper's numbers; scaling-exponent fits reported.

## Log

(dated entries appended newest-last)

### 2026-08-17 — FGS leaf-LU cache landed; local R1 A/B verified

FastMultipole now preserves the original contiguous self-influence buffer and,
by default, builds one second contiguous buffer containing checked LU factors.
`cache_leaf_lu=false` retains the historical per-sweep `mat \\ rhs` path.
FLOWPanel propagates the keyword through `FGSSolver` and `FGSPreconditioner`
and records it in metadata. Factor build time and retained numeric/pivot bytes
are exposed on `fgs.leaf_lu_cache`.

The corrected local R1 smoke used the frozen Phase-1 τ=1e-6 knobs, 4 Julia
threads / BLAS 1, cold starts, min-of-5 after warmup. Constructor compilation
was excluded from setup timing. Cached/uncached iteration counts and certified
BC errors were identical. Leaf pass: 0.0129974 -> 0.00046275 s (28.1x);
FGS solve: 1.19617 -> 0.455358 s (2.63x); FGMRES+FGS solve: 3.30494 ->
2.22390 s (1.49x). Factor payload+pivots: 10,656,608 bytes (10.16 MiB), built
in 9.77--11.10 ms. Raw rows and histories:
`benchmark/results/phase2/leaf_lu_cache/multi/`.

This is local smoke evidence, not the published Phase-2 result. R4--R7 still
lack frozen Phase-1 settings, so the dedicated-node R1--R7 cache A/B and the
cache-enabled production retuning remain pending rather than inventing knobs.
The complete FastMultipole and FLOWPanel regression suites pass at 4T/BLAS-1.

### 2026-08-18 — Phase 2 prep landed concurrent with Phase 1 (Ryan authorization; harness + all staged src levers)

Ryan (2026-08-18, in session): do the Phase 2 prep that can run concurrently
with the in-flight Phase 1 HPC chains, scope = harness + scripts + ALL staged
src levers — his selection explicitly covered the Krylov FMM-tree cache AND
`fgs_determinism_performance_plan.md` Part C2/C3 (colored parallel sweeps,
batched `update_nonself_influence!`), i.e. **the sign-off C2/C3 were gated
on**; plus a corrected unsteady campaign driver in `benchmark/` modeled on
RHPC (`examples/rotor_hover.jl` untouched). Isolation held throughout:
nothing rsynced to the cluster; Phase 1 scripts of record untouched except
one pure code-motion extraction; all local runs ≤4T dev/smoke
(hardware_tag-distinguished). Plan file:
`~/.claude/plans/resume-work-on-021-functional-valley.md`.

**Harness (all smoke-validated at R1, multi @4T):**

- `benchmark/phase1_knobs.jl` — knob-selection + adaptive-min-of-k helpers
  extracted verbatim from `..._phase1_table.jl` (which now includes it;
  refactor verified by a `backslash_ldiv` R1 row reproducing the historical
  shape: setup 7.48 s / solve 0.0182 s / BC 4.0e-9 certified).
- `benchmark/rotor_hover_solver_phase2.jl` (Stage 2a driver) — componentized
  setup (dense assembly on a scratch G with LU derived as ctor−assembly to
  avoid double O(N³); Jacobi/ILU/FGSPreconditioner builds split out; ILU
  tree/lists/assembly/factorization from its stats dict; FGS leaf-LU cache
  build time/bytes), adaptive-min-of-k cold solves (ruling 12), history
  sidecars with wall-clock stamps for every iterative config, memory columns
  (`solver_state_bytes` + `@allocated`). All six configs certified at R1
  (`results/phase2/multi/phase2.csv`); the FGS rows required regenerating the
  local R1 `SWEEP_LADDER_1E6` fgsprecond rows for the post-LU-cache-tuned
  config 6/0.3/50/10 (the 08-17 killed re-tune had left `fgstune_selected`
  pointing at a config with no sweep-ladder rows — `stage3_winner`'s loud
  error caught it, as designed). fgs 0.459 s / 7 it / BC 9.5e-8 — consistent
  with the pre-kill LU-cache observation (0.55 s / 7 it).
- `benchmark/rotor_hover_solver_unsteady.jl` — includes RHPC with
  `RHPC_SETUP_ONLY=1` (construction NOT duplicated → the stale-convention bug
  class is structurally excluded), `CONFIG` env selects the roster solver at
  Phase-1 frozen knobs, per-step instrumentation with NO src changes: a
  `TimedFormulation` wrapper times `solve_formulation!` (the body solve
  inside `_steady_aerodynamics!`, incl. the per-step Dirichlet source
  assembly inside `solve!`), full-step wall from `maneuver!` deltas.
  R1 smokes: backslash 5 steps / solve share 36.8%; fgs 40.2%; identical
  249-particle wake both configs. Rows in `unsteady.csv`. Phase 3 reuses this
  driver (ruling 12). LIMITATION flagged: FGS per-step iteration counts need
  callback plumbing through `solve!` (niter=-1 recorded) — a Phase 3
  prerequisite, since iterations-to-target is the warmstart headline metric.
- `benchmark/slurm/p2_table.sh`, `p2_unsteady.sh` — prepared, NOT submitted
  (submission waits on the Phase 1 freeze + Ryan; zen3+exclusive documented).
- `benchmark/rotor_hover_solver_phase2_profile.jl` — Profile + sampled
  Profile.Allocs per config → flat/tree/alloc-site files + summary CSV
  (`results/phase2/profile/`); R1 smoke: krylov_ilu dominated by the
  per-apply `fmm!` (the tree-rebuild lever), fgs by `body_to_multipole`.
- `benchmark/phase2_analysis.jl` — merges flat + per-rung CSVs (latest row
  per key), emits per the figure policy into this item's `figures/`:
  `solver_scaling.tex` + data dir + fitted exponents (dev R1–R3: fgs solve
  ~N^1.15, dense LU setup ~N^2.3), `convergence_<rung>.tex` from history
  sidecars, equal-time-per-iteration validation (ruling 4: all R1 iterative
  configs flat within 5%), memory CSV (ruling 8). Both figures compile under
  pdflatex. Regenerate from HPC CSVs when they land — the pipeline is the
  deliverable, current numbers are dev-only.

**Src levers (local only — NOT rsynced). FastMultipole side committed on
`flowpanel-20260817`: `b3cf4ad` (FmmPlan) + `752c525` (colored sweeps),
Ryan's unrelated dirty script untouched; FLOWPanel src changes left
uncommitted per repo convention (Ryan commits FLOWPanel):**

- **Krylov FMM-tree cache**: FastMultipole gains `FmmPlan` (Cache + both
  trees + sorted lists + switches built once; `fmm!(targets, sources, plan)`
  refreshes source strengths in sorted order, zeros target outputs, and
  dispatches to the prebuilt-tree method — a per-apply strength refresh the
  phase docs' assumption missed: buffers fill only at Tree construction).
  FLOWPanel gains `cache_tree=false` on `KrylovSolver`: per-SOLVE scope —
  plan built on first apply, dropped at `_krylov_launch!` exit — so geometry
  and the active kerneloffset (which feeds buffer radii via
  `radius_inflation`) are frozen for the plan's lifetime by construction;
  ctor rejects non-FMM backends; `cache_tree` recorded in solver metadata.
  Tests: FastMultipole `fmm_plan_test.jl` 11/11 (bitwise plan-vs-fresh across
  strength changes + body-count guard); FLOWPanel solver unit suite 235/235
  incl. the new 11-test set (bitwise cached-vs-uncached solve, allocation
  bound with premise floor, geometry-change invalidation — the kerneloffset
  flip premise was unmeasurable on the sphere fixture (offset inert there),
  so the discriminating input is a node perturbation). `KrylovCoupled`
  deferred (same pattern, optional commit 3).
- **Colored FGS sweeps (`sweep_order=:colored`), Part C3 folded in**:
  `update_nonself_influence!` split into thread-safe per-leaf products +
  serial scatter (pure code motion on the default path);
  `color_leaves` greedy coloring on rhs row-interval conflict (deterministic,
  thread-count-independent); `gs_sweep!` runs per-color parallel
  solve+products then a serial ascending scatter. Part C3 is NOT a standalone
  mode: the per-sweep direct-list touch was already one pass (the plan doc's
  premise was stale), and per-sweep batching would degrade GS→Jacobi — its
  value is the race-free per-color batch inside C2. Tests
  (`fgs_coloring_test.jl`, 2216/2216 @4T): coloring validity (premise-guarded
  non-vacuous), the batching THEOREM (colored sweep bitwise ≡ sequential
  color-major immediate GS), cold fixed-iteration repeatability, iteration
  health. **Measured finding: multicolor GS is a weaker iteration —
  it DIVERGES on the synthetic gravity case (residual ~1e149 @30 sweeps)
  where lexicographic converges.** Adoption is therefore hard-gated on the
  campaign convergence A/B (fgstune staircases R1–R3 + determinism-probe
  matrix rerun) before `:colored` appears in any campaign config; ships
  opt-in, default `:lexicographic` bit-identical. Ryan's scheduling question
  answered by A/B on the real sweep structures (8k/16k bodies, 512 leaves,
  34–37 colors, 4T): default `Threads.@threads` beats `:static` (2.75 vs
  3.06 ms; 7.00 vs 7.25 ms) — default kept; determinism is scatter-order
  based, not schedule-based. Legacy quirk found and preserved: the historical
  `reverse_pass` loop iterated `enumerate(reverse(leaf_index))` but used only
  the counter, so it was never actually reversed — behavior kept bit-identical,
  flagged for Ryan. FLOWPanel plumbing mirrors `cache_leaf_lu`
  (`FGSSolver`/`FGSPreconditioner` kwarg + metadata; unit-tested).

Suites at end of session: FLOWPanel solver unit 235/235; FastMultipole
`fmm_plan_test` 11/11 + `fgs_coloring_test` 2216/2216 standalone; full
FastMultipole suite result recorded in `log.md` (run last, after all edits).
Pending Phase 2 work unchanged: dedicated-node ladder benchmarks (blocked on
Phase 1 freeze), leaf-LU cache HPC A/B, Stage 2b adoption decisions.

### 2026-08-20 — Phase-2b HPC chains LAUNCHED (nfcache campaign)

Ryan's go. **Jobs of record (per-rung afterok chains, zen3 exclusive,
multi/64T)**: R1 13242659→13242660, R2 13242661→13242662,
R3 13242663→13242664, R4 13242665→13242666. Job A per rung =
`p2b_nearfield.sh` (full §5 A/B: isolated near-field + gmres solve
cache on/off, THEN the cached tune with full MAC sweep 0.3–0.8 →
`tune_cached.csv`); job B = `p2_table.sh` with
CONFIGS=krylov_gmres:krylov_ilu:fgmres_fgs + the three `_nfcache` twins
(R4 drops the fgs-family — no Stage-3 winner CSVs for R4). Outputs land in
`benchmark/results/phase2/multi/<rung>/` on the cluster (PER_RUNG_DIR; the
phase2b script was patched to honor it so tune_cached.csv sits where
phase2.jl reads it; TUNE_MACS accepts colons for sbatch --export).

**Ops judgment calls (logged)**: (1) Phase-1 chains still in flight
(R5-single running, R6 tune running/tables pending) ⇒ the campaign ran from
a SEPARATE cluster checkout pair `~/projects/p2b/{FLOWPanel.jl,FastMultipole}`
— main checkout src/results untouched (evidence integrity); code rsynced
from local (FLOWPanel uncommitted tree + FastMultipole @ 1ec0af9), meshes
R1–R4 + `results/phase1/multi/` knob CSVs/bcaches COPIED cluster-side
read-only, Manifest copied from the main checkout (no resolve).
(2) FLOWVPM path dep satisfied by symlink `~/projects/p2b/FLOWVPM.jl` →
main checkout (read-only, untouched by this work). (3) FastMultipole
`src/error.jl` includes `../test/evaluate_multipole.jl` ⇒ FM `test/*.jl`
must ship with src (gotcha for future partial syncs). (4) PythonPlot/
CondaPkg precompile fails in the p2b env (no `.CondaPkg`, 1.1 GB — not
copied); benchmark path never plots ⇒ jobs set
`JULIA_PKG_PRECOMPILE_AUTO=0` (FLOWPanel/FastMultipole/GeoIO verified
LOAD_OK_WITH_CACHE on the login node). (5) R5–R7 excluded until their
Phase-1 knob chains land. Harvest: merge per-rung phase2 CSVs to local,
check bc_certified on every row, then A/B tables + ledger.

### 2026-08-19 (evening) — nfcache benchmark configs + rigid-motion staging (Ryan directives)

Ryan directives, session dialogue: (1) add one near-field-cache variant PER
GMRES-family case — `krylov_gmres_nfcache`, `krylov_ilu_nfcache`,
`fgmres_fgs_nfcache` (Jacobi skipped per his enumeration); (2) re-tune for
the cached path — ONE cached tune per rung SHARED by all three variants
(the build-blind tuner's knobs are case-independent; ruling-3 style);
(3) convention match confirmed: per-solve state is inside `t_solve_min`
("tree rebuild per Krylov apply included"), so the cache build — per-solve
state, since a rotating rotor invalidates it every step — is included in
per-step and reported in row notes, never in setup columns; tuning
disregards one-time solver builds (they were never in the tuner objective).

**Landed (benchmark scripts only, uncommitted)**: phase2b script's
TUNE_CACHED stage now also writes machine-readable
`results/phase2/<mode>/tune_cached.csv` (csv-first, latest row per rung);
`rotor_hover_solver_phase2.jl` gains `tuned_cached_knobs()`/`nfcache_setup()`
(shared backend + one throwaway cache build for the notes), an
`apply_knobs` kwarg on `measure_and_emit!`, the two Krylov nfcache names in
`run_krylov!`, `run_fgmres_fgs_nfcache!` (Stage-3 winner precond knobs
untouched — only the operator apply is cached), and the three names in the
default CONFIGS + dispatch.

**R1 local smoke (multi 4T, cached shared knobs p=14/MAC=0.5/leaf=326 from
tune_cached.csv; cache 810.1 MB / 31.7 s build per solve; all
bc_certified=true, cold min-of-5)**:

| config | t_solve | niter | bc rel-L2 |
| --- | --- | --- | --- |
| krylov_gmres_nfcache | 33.79 s | 59 | 9.67e-7 |
| krylov_ilu_nfcache | 33.05 s | 7 | 8.30e-7 |
| fgmres_fgs_nfcache | 34.12 s | 1 | 5.88e-7 |

**Reading: BUILD-DOMINATED across the board** — ~94% of every cold solve is
the per-solve cache build; with applies nearly free at the cached knobs, the
preconditioner choice is worth <1 s (59 vs 7 vs 1 iterations,
indistinguishable wall). This is the per-solve-rebuild economics made
explicit, and it points squarely at the staged rigid-motion persistence
item as the real lever (build amortized over the run would leave ~2 s
solves). HPC rows + per-rung cached tunes (full TUNE_MACS sweep) wait for
the window.

**STAGED (Ryan's go required): `rigid_motion_tree_reuse_item.md`** —
`transform_tree!`-level rigid-motion reuse serving `FmmPlan` AND
`FastGaussSeidel` (shared machinery per Ryan), unlocking cross-timestep
plan+cache persistence for the rotor (Dirichlet cache blocks exactly
rotation-invariant). **Includes a CORRECTNESS FLAG (flag-only per Ryan,
verify inside that item): FGS unsteady staleness** — `FastGaussSeidel`
trees/lists/matrices are built once at ctor and never refreshed while the
unsteady driver rotates the body ⇒ far-field expansions about t=0 branch
centers, error grows with rotation angle silently; steady campaign
unaffected; NO Phase-3 unsteady FGS row is trusted until the discriminator
runs. (Krylov is immune: trees rebuilt per apply/solve.)

### 2026-08-19 (later) — Ryan feedback round: cached-path autotuning + build-time guard

Ryan reviewed the landed cache and ruled (dialogue in-session): (1) cache
ownership STAYS on `FmmPlan` (Tree-field idea examined and rejected — the
validity key is the tree PAIR + sorted direct list + switches, which is what
the plan bundles); (2) `tune=true` refusal SUPERSEDED — cached-path tuning
approved at full scope; (3) build-time safety must be ESTIMATED up front
(one kernel sample), never discovered by timing the build; (4) cap semantics
during tuning = "stop growth at the cap" for both `max_bytes` and
`max_build_time`; (5) with a PROVIDED cache, `tune=true` DOES return a
leaf_size_source suggestion (from the cached per-interaction timing) — it is
the user's call whether to rebuild trees + cache at the suggested leaf.
(A first implementation restricted provided-cache tune to expansion_order
only, misreading Ryan's ruling as fmm!-level when it concerned `tune_fmm`;
corrected in FM 1ec0af9 same session.) The `tune_fmm` throwaway-cache route
remains the leaf/MAC tuning path that manages rebuilds itself.

**Landed (FastMultipole)**:
- `0ef4e83` — `estimate_nearfield_cache(...)` → `(; bytes, est_build_time,
  n_blocks, total_probe_pairs)`: bytes from the exact builder size pass
  (shared helper, cannot drift); `est_build_time` = warmed-up min-of-3
  single-source kernel sample × probe-pair count. Non-throwing (tuner
  feasibility checks). Both ctors + `build_nearfield_cache!` gain
  `max_build_time::Real=Inf`, checked BEFORE probing.
- `204188a` (+ correction `1ec0af9`) — `tune=true` allowed on the cached
  path. Provided cache: full suggestions returned incl. leaf (see ruling 5).
  `tune_nearfield_cache=true` on `fmm!`/`tune_fmm`
  (+ `nearfield_cache_max_bytes`/`_max_build_time`): throwaway cache per
  tuning call, build cost excluded from `t_direct` AND subtracted from the
  tuner's wall-clock MAC comparisons; over-cap trials are not built —
  `optargs.nearfield_cache_feasible=false` makes `tune_fmm` pin
  `leaf_size_source` to the last cache-feasible trial (direction-agnostic;
  measured at n=800 the model can suggest SMALLER leaves with MORE bytes).
  First-trial-infeasible errors loudly. `nearfield_cached!` now attributes
  cached time per source system ∝ block interaction counts under tune.
  Deviation from the approved plan text (logged): `cache_capped` is returned
  as a THIRD `tune_fmm` return value `(; cache_capped)` instead of inside
  `tuned_params`, because `tuned_params` is documented as splat-able into
  `fmm!` kwargs and an extra field would break that contract; existing
  2-value destructuring is unaffected.

**R1 cached-economics tuner smoke** (TUNE_CACHED=1 section added to
`benchmark/rotor_hover_solver_phase2b_nearfield_cache.jl`, MAC fixed at 0.5,
caps 4 GiB / 60 s build): tuned knobs **p=14, MAC=0.5, leaf=342** vs the
kernel-tuned frozen p=17/MAC=0.5/leaf=21 — the leaf optimum grows ~16× under
cached economics exactly as the `sqrt(t_m2l/t_direct)` model predicts;
tuner-reported apply cost 0.0627 s (excl. build); cache at tuned knobs
862.1 MB, est. build 33.9 s; `cache_capped=false`; tune wall 85.1 s. Row
`tune_cached` appended to `nearfield_cache_ab.csv` (local-smoke). GMRES
re-tune of record still belongs to the HPC window.

**Suites**: FastMultipole full suite exit 0 (54 cache tests incl. 17 new
tuning tests); FLOWPanel solver units 256/256; fmm units 31+18 — all 4T.
The old "tune refused" test was replaced by the P-only provided-cache test.

**Ruling (Ryan, same session): small-n false low-leaf optimum = ACCEPTED
AS-IS, no mitigation.** Measured cached-apply leaf sweeps (min-of-20, 4T):
n=800/p=6 — tuner pick leaf 10 = 4.15 ms vs true optimum leaf 160 = 0.88 ms
(4.7× off, but single-digit ms absolute; the true small-n optimum is
"leaf ≥ 20 empties the m2l list" — pure dense, a regime boundary the
sqrt(t_m2l/t_direct) model structurally can't see). n=8000/p=8 — flat
plateau 35–37 ms across leaf 50–342 streaming at ~56 GB/s ≈ memory
bandwidth, i.e. per-interaction cost is 8 B/bandwidth and leaf-INDEPENDENT
in the memory-bound regime, so the model's assumption is valid exactly at
campaign scale and its pick (342) lands mid-plateau. Real hazard is
OVERSHOOT (leaf 600 = 79.7 ms, 2.3× worse: ~13 leaves starve the
owner-partitioned threads), but the model doesn't move upward and the TUNE
stage's perturbation double-check guards it on benchmark hardware. A
bracket-verify step ({pick, 2×, 4×} measured, keep fastest) was designed
and OFFERED but declined as unnecessary — revisit only if a tuned pick ever
lands off-plateau on HPC.

**p-selection verification (Ryan question, measured R1+R2): the cached
tuner's expansion order is NOT warranted — descend p at the tuned leaf.**
R1 @ leaf=342/MAC=0.5 (m2l only 28 entries, nearly dense): p=4 already
rel-L2 1.8e-7 vs p=20; tuner's p=14 costs +15 ms/apply for nothing. R2 @
tuned leaf=196/MAC=0.5 (real far field: m2l=516, cache 1.6 GB; tune
`capped=true` — the 60 s build-time cap bound the leaf): p=4→2.4e-6,
**p=6→3.1e-7 (meets the 1e-6 target)**, p=8→9.8e-9, vs tuner p=15→1.6e-11
at 77.9 ms — ~4 orders tighter than required at 1.8–2× the p=6–8 apply
cost (38.8–43.4 ms). R3+ verification folded into the HPC retune (R2 local run took
~17 min; R3 would breach the 20-min local rule). Error-request semantics:
the tune calls asked for ABSOLUTE `PowerAbsolutePotential(0.1·1e-6·rms_b)`
(campaign convention). **MAX-error check (Ryan, R2): the algorithm is
nearly VINDICATED** — max pointwise abs error at leaf=196: p=6 = 4.90e-7
(414× OVER the requested 1.18e-9), p=8 = 24× over, p=15 = 5.02e-11 (meets,
24× margin ⇒ minimal compliant p≈12–13, only ~1 decade of bound
conservatism). The 4-decade gap to the campaign metric lives in the REQUEST
construction (0.1 safety + max-abs-per-point vs rel-L2-BC intent +
rms_b-scale epsilon on an O(1) linear-operator field). Fix item restaged
accordingly as `tune_fmm_p_selection_item.md` (epsilon-construction
semantics, Ryan's go required).

**RULING (Ryan, closing the thread): KEEP the tuner's full p.** Since the
worst-case point genuinely needs the picked order under the requested
max-abs criterion, and the cost of honoring it is ~2× on the apply (R2:
77.9 ms at p=15 vs 38.8–43.4 ms at p=6–8), the ~2× is accepted as
insurance for that worst-case point. RETUNE PROTOCOL for the HPC window is
therefore: adopt the cached tuner's knobs AS-IS — (leaf, MAC, p) — no
descending-p step; per-run BC certification remains the guard as always.
An earlier note in this entry prescribing a descending-p sweep is
superseded by this ruling (revisit only via the staged
epsilon-construction item, which stays open for a future session).

### 2026-08-19 — Near-field influence-matrix cache IMPLEMENTED (commits 1–3; commit 4 deferred)

Executed `nearfield_matrix_cache_plan.md` (status flipped to EXECUTING; Ryan's
scheduled launch = the go). All three implementation stages landed; the
plan's commit 4 (`persistent_plan` cross-solve persistence) is explicitly
deferred pending Ryan's go — stopped before it.

**What landed** (FastMultipole `flowpanel-20260817`, on top of 752c525):

- `72d3f3d` — `NearfieldInfluenceCache` (new `src/nearfield_cache.jl`):
  one dense block per (direct-list entry × system pair) in a reused
  `Matrices` container, target-major order; built by unit-strength column
  probing through `direct!` with `value_to_strength!` bypassed (direct
  buffer-row writes so ALL strength components get columns). Evaluator
  `nearfield_matvec!` gathers strengths, runs one BLAS matvec per block into
  disjoint `Matrices.rhs` slices, accumulates into target output rows;
  threading is owner-partitioned with cuts only at target-branch boundaries
  ⇒ deterministic and bitwise thread-count-invariant (tested at 1T and 4T).
  Guards: `max_bytes` (default 4 GiB) checked BEFORE allocation;
  `direct_conditioning` refused at build AND eval; Tree `objectid` is the
  whole validity key (`check_cache_trees`) + body counts at eval. Standalone
  `direct!` gains `nearfield_cache=` with a degenerate treeless form.
  Tests `test/nearfield_cache_test.jl`: exactness vs kernel (rtol 1e-12,
  bitwise-inequality premise proves the BLAS path ran), strength-change
  reuse, tree guard, conditioning refusal, memory cap, determinism.
- `02f071c` — `fmm!` integration: prebuilt-tree `fmm!` gains
  `nearfield_cache=nothing`; `FmmPlan` carries a `RefValue` cache slot filled
  by `build_nearfield_cache!(plan, targets, sources)` and picked up
  automatically by `fmm!(targets, sources, plan)`. `tune=true` is refused on
  the cached path (per-interaction kernel timing unmeasurable). Test:
  plan-with-cache fmm! == plan-without (rtol 1e-12) across strength changes.
- FLOWPanel surface (UNCOMMITTED per the launch constraint — working tree
  only, on top of the 2026-08-18 prep edits): `KrylovSolver(...;
  cache_nearfield=true)` (implies `cache_tree`; DirectBackend rejected),
  plumbed `KrylovOperator` → `_apply_*_G!` → `influence!` where the cache is
  built lazily alongside the fresh per-solve plan; metadata key
  `"cache_nearfield"`. New testset in `runtests_unit_solver.jl` (26 tests).

**V0 (`_induced_wake` linearity): PASSED — refusal path not needed.**
Dirichlet diamond-wedge body WITH shedding panels (nspan=40, finite attached
wake, Das≠0): cached `_apply_dirichlet_G!` matches the kernel FMM apply to
rel max err 8.6e-16; R1-scale (8016-panel rotor, real TE shedding) apply
matches to rel-L2 4.3e-16. The attached-wake term is linear in the buffer
strength rows as expected.

**R1 local smoke A/B** (Phase-1 frozen knobs p=17/mac=0.5/leaf=21, multi 4T,
`benchmark/results/phase2/multi/nearfield_cache_ab.csv`, hardware_tag
local-smoke; publishable numbers wait for the HPC window):

| Measurement | kernel | cached | ratio |
| --- | --- | --- | --- |
| isolated near-field pass (min-of-7, shared plan buffers) | 1.4331 s | 0.005406 s | **265×** |
| cold `krylov_gmres` solve (min-of-3, atol 1e-10/rtol 1e-8) | 130.97 s | 23.72 s | **5.52×** |

Cache build 11.89 s, 272.4 MB (dense-G 8N² baseline at R1 = 514 MB — the
cache is 53% of matrix-ful storage while staying matrix-free in the far
field), break-even **8.3 applies** (vs 83 gmres iterations ⇒ pays off ~10×
over within ONE solve). niter identical (83/83); certified BC rel-L2
8.4747e-9 both, `certified=true` both (agree to 8 digits); solution shift
rel-L2 3.5e-15. Script: `benchmark/rotor_hover_solver_phase2b_nearfield_cache.jl`
(new, phase2 pattern; also uncommitted). Note the uncached leg's 131 s
includes per-apply work the plan already absorbs — this is the
cache_tree-on baseline (both legs `cache_tree=true`), so the 5.5× is the
near-field-cache increment alone.

**Suites at end of session (all green, 4T)**: FastMultipole full suite
(incl. 30 new cache tests; exit 0); FLOWPanel `runtests_unit_solver.jl`
256/256; `runtests_unit_fmm.jl` 31+18.

**Judgment calls (flagged for Ryan, made without blocking)**:
1. "Commits 1–3" vs "commit only in FastMultipole": commits 1–2 landed in
   FastMultipole; the plan's commit 3 (FLOWPanel surface) is implemented but
   left uncommitted in the working tree, consistent with the inherited state.
2. Plan doc erratum: `reset!` zeroes buffer rows 4:end (`reset_outputs!` is
   the exact-`output_range` one); the builder zeroes exactly the probed
   output rows, so target metadata rows are never clobbered.
3. Cache struct stores per-block target/source/output ranges at build so
   EVALUATION needs no trees at all (enables the treeless standalone
   `direct!` form); per-thread gather scratch is preallocated at build
   (eval thread count capped at build-time `Threads.nthreads()`).
4. The plan's test 2 called for an R1-class body in the unit suite; used the
   asset-free 320-panel shedding diamond fixture there (identical code path,
   seconds not minutes) and moved the R1-scale exactness assert into the
   smoke script (`@assert rel_l2 < 1e-12`, measured 4.3e-16).
5. Solve-equality tolerance in the unit test relaxed from the plan's 1e-12
   to rtol 1e-10 (Krylov tolerance amplification of the BLAS-order
   perturbation); measured R1 shift is 3.5e-15 anyway. Sphere-solve tangency
   residual asserted at 1e-6, which is FMM p=8 truncation (measured
   identical, 2.904e-8, on cached AND uncached).
6. `tune=true` + cache → loud `ArgumentError` (refuse-and-error convention).
7. Cached near-field timing is reported in `t_direct` slot 1 only
   (per-source-system attribution is meaningless for the fused matvec).
8. Smoke A/B scope: isolated near-field + `krylov_gmres` end-to-end only;
   the full §5 matrix (jacobi/ilu/fgmres_fgs × on/off) belongs to the HPC
   window where the numbers are publishable anyway.
9. Local run used Julia 1.12.5 (the default `julia`); NO Pkg resolve was
   performed (the 1.11.7-Manifest gotcha from 018 concerns re-resolving —
   Manifest untouched, tests/benchmarks all green).

### 2026-08-19 — Near-field influence-matrix cache STAGED (plan only, per Ryan)

Ryan directive: new Stage 2b lever — cache near-field direct interactions as
packed dense influence blocks (FGS `Matrices` single-vector storage reused;
avoid many small allocations), **general** (any `direct!` kernel),
**cache-able** across timesteps for frozen relative geometry, keyed to a
`Tree`; priorities performance > robustness > user-friendliness. NOT
implemented — staged as a self-contained plan for a fresh clear-context
agent: [`nearfield_matrix_cache_plan.md`](nearfield_matrix_cache_plan.md).
Design facts verified in-session: the FGS matrix builder is scalar→scalar
only (a new full-`DerivativesSwitch`, multi-strength-component builder is
specified); exported `InteractionList` is dead code (ctor calls an undefined
function) — flagged, not reclaimed; target-major block keying gives an
owner-parallel scatter-free evaluator (opposite of FGS's source-major +
serial scatter); `direct_conditioning` and Tree identity define the refusal/
validity contract; V0 verification item = `_induced_wake` linearity in the
buffer strength rows (exercised by the plan's mandatory shedding-panel
exactness test). Composes with the landed `FmmPlan`/`cache_tree` lever;
timestep-to-timestep reuse is stage 2 (needs the cross-solve
`persistent_plan` extension).
