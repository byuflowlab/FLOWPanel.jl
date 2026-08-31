# 021 Session Log

Newest first. Narrative only — results go to the phase files and `ledger.md`.

## Dated entries

### 2026-08-24 (latest) — campaign ladder set to machine classes; HPC launch handoff written

Ryan set the memory ladder in terms of machines people actually have rather
than an arbitrary sweep: *"16G for a laptop, 128 for a workstation, and 500 for
a supercomputer."* Final campaign setting is `MEM_BUDGETS = 0:16:128:500` GiB of
TOTAL memory, where **0 is a forced-uncached sentinel** (cache forbidden,
memory unconstrained) rather than a machine. Budget 0 is REQUIRED and must run
first: it is what tunes the uncached configs, and at R1-R5 every real machine
can afford a cache, so no positive budget ever returns an uncached winner.

Why this beats the geometric ladder it replaced: a flat `4:8:...:500` sweep was
checked against the W0 map first and yields only **13 distinct operating points
out of 56 descents** (R1-R2 saturate at every budget; R7 is uncached at five of
eight). The machine ladder is 4 x 7 = 28 descents and its duplicates ARE the
result. Per-rung meaning: R1-R3 all three classes saturate; R4-R5 laptop
constrained; R6-R7 laptop cannot cache at all. R6/R7 saturate at 158 and 428
GiB, both under the 500 GiB node, so the whole axis is reachable on pinned
hardware — no need for the 1 TB (`cs`) or 2 TB (`m11-2`, zen2) partitions,
which would confound the axis with a hardware change.

**CORRECTION to the same-day floor claim** (detail + table in `ledger.md`): I
had read the W0 map's smallest sampled leaf as a floor and written that cache
bytes "never fall below the small-leaf value". That was extrapolation past the
sampled range. Measured: bytes keep falling below leaf 9 — R6 23.8 -> 19.3 GiB
(-19%), R7 78.2 -> 67.0 GiB (-14%) at leaf 2 — while clearly asymptoting toward
a real MAC-set minimum. The R6/R7 conclusion survives (a 16 GiB laptop still
cannot cache) but the R6 margin is 23%, not the ~50% the quoted 24 GiB implied.
The substantive consequence is a code fix: total memory has an INTERIOR
minimum, because the FMM plan grows as the leaf shrinks while the cache
shrinks, so pricing a floor at any single leaf is wrong in both directions. The
tuner's floor check now SCANS `FLOOR_LEAVES` (2:4:9:21) and minimises the
actual total, reporting cheapest-uncached and cheapest-cached separately.

**Confirmed to Ryan, from code and measurement, that the cache build is NOT in
the tuning objective**: structurally the warm-up solve that builds the cache is
not wrapped in `@elapsed` (only the min-of-reps loop is, with `reset_cold!`
outside the timer), and empirically at R1 the build is 11.7-14.2 s while
`t_warm` is 4.8 s — impossible if the build were included.

**Launched nothing.** Wrote two Slurm wrappers reusing `p1_tune.sh`'s pinned
hardware preamble verbatim — `benchmark/slurm/p2_tune.sh` (one job per rung x
budget) and `benchmark/slurm/p1_agreement.sh` — and a self-contained launch
handoff at `phase_16_hpc_launch_prompt.md`. Everything is implemented and
locally smoke-tested; nothing has ever run on the cluster. `sbatch` remains
Ryan's call.

### 2026-08-24 (latest) — error metric is the BC, not a reference; memory budgets are TOTAL memory

Two Ryan corrections to the same-day work below, both of which simplify things.

**1. Measure error against the boundary condition, not a reference.** Ryan:
*"we should be measuring error based on the boundary condition - not a
reference. That way, it doesn't matter."* Correct, and the instrument was
already in the driver: `true_residual!` evaluates `||Ax-b||` through the
DirectBackend on each config's OWN solution, so `rel_rms_direct = rms_d/rms_b`
is an exact, reference-free correctness metric, identically defined at R1 and
R7.

This **dissolves the R6/R7 blocker** raised in the entry below — a dense
`backslash_ref` would be 335 GiB at R6 and 1.31 TiB at R7, and now nothing
requires it. `rotor_hover_solver_phase1_agreement.jl` changes:
- the BC residual is the primary exit criterion; reference-relative columns
  (`x_relL2_vs_ref`, `dCT_*_pct`) are now OPTIONAL and blank unless
  `REFERENCE=<config>` is nominated. The dense-reference guard is gone, and
  **config order no longer matters** (previously `backslash_ref` had to run
  first or the script errored).
- cross-config CONSISTENCY is now measured against the ENSEMBLE, also
  reference-free: every config persists its solution
  (`agreement_x_<rung>_<config>.bin`) and a summary pass reports each config's
  deviation from the ensemble mean plus the **max pairwise** relative-L2
  difference → `agreement_spread.csv`. Works across separate processes and
  across arbitrary config subsets.
- Backslash is still worth running where it FITS (25 GiB at R4) as an
  independent non-iterative cross-check — just not required and not
  privileged.

**BC metric is evaluated through the FMM, not DirectBackend** (Ryan, same
session): the direct evaluation is O(N^2) — 1.8e11 pair evaluations at R7,
prohibitive. `bc_error!` (already the campaign's BC evaluator, validated by
`rotor_hover_solver_phase1_bcerror.jl`) evaluates it through the FMM under an
EXPLICIT `PowerAbsolutePotential(safety * target_rel * rms_b)` with
`safety=0.1` — an order of magnitude tighter than the 1e-6 criterion under
test — and returns `error_success`, the FMM's own certification that it met
that tolerance. The driver records `bc_certified` per row and warns when it is
false, so an uncertified evaluation is visible rather than silently reported
as an error. `RESIDUAL_BACKEND=direct` restores the exact O(N^2) evaluation,
worth running at R1-R2 as a cross-check that the FMM metric matches the exact
one. The Phase 2 tuner already used `bc_error!` on this contract, so no change
was needed there.

Also hardened while in the file: `agreement.csv` gains a header guard and all
rows now go through `emit_row!`, which asserts the row width matches the
header — the FAILED-path row had been hand-counted blanks, which is exactly
the kind of thing that silently shifts every column when a schema changes.

Caveat recorded in the driver: at matched BC residual the solution difference
is bounded only up to `cond(A)`, so the residual demonstrates CORRECTNESS and
the spread table demonstrates CONSISTENCY. Both are needed; neither needs a
reference.

**2. Memory budgets are TOTAL memory, and rungs that admit no cache get tuned
uncached.** Ryan: *"if memory bandwidths prohibit a cache, then that particular
run will need to run without a cache... we'll need to tune those rungs that
don't admit a cache just so - without a cache"* and *"memory rungs should
correspond to the total memory available - not just that used by the cache."*

Both landed in `rotor_hover_solver_phase2_tune.jl`:
- a budget now covers **body + solver state + FMM plan + cache**, using
  `solver_state_bytes` (the ruling-8 comparable metric, which excludes the body
  by construction — hence the separate term). This puts `backslash_ldiv`'s
  dense 8N² on the SAME axis and turns the question into "which solver is
  fastest on a node with X GB". It also matters mechanically: **the plan grows
  as the leaf shrinks**, partly offsetting the cache's growth, so pricing the
  cache alone would misrank small-leaf candidates.
- **whether a candidate is cached is now an OUTCOME, not an input**: the cache
  is used iff the whole configuration still fits. A rung or budget admitting no
  cache is therefore tuned uncached automatically — no separate code path, no
  separate driver. R6/R7 need exactly this (cache floors 24 and 78 GiB).
- budgets below a rung's floor (even uncached) are reported as the infeasible
  end of the axis; budgets whose winner duplicates the previous budget's are
  flagged `DUPLICATE` so the curve is not read as having independent evidence
  at every point.
- schema gains `mem_total_predicted` / `mem_total_measured`. The budget is
  enforced against the PREDICTED total because the plan is live during a solve
  even when `cache_tree=false` drops it afterwards, so a post-hoc measurement
  would under-report.

`rotor_hover_solver_phase2.jl` follows the tuned row's `cached` flag rather
than assuming a `*_nfcache` config is cached, and records
`mem_total_measured`.

**Ladder top = 500 GiB.** Node inventory (read from the Slurm config, NOT live
`sinfo` — verify before relying on it): the campaign's pinned physics2 nodes
`phys-1-[1-10]` are zen3 / 128 cores / **512 GB**, which is what `--mem=500G`
targets. Larger exists — `cs` (zen3, 1 TB) and `m11-2` (zen2, 2 TB) — but
reaching either means changing hardware, which folds a CPU change into the
memory axis and destroys timing comparability across the ladder. Extending past
500 GiB needs Ryan's go plus an overlapping-budget cross-check on the new
nodes; m11-2 is zen2 and would confound the axis outright.

### 2026-08-24 (latest) — caching wired into Phase 2 tuning AND benchmarks; memory budget becomes an AXIS; Phase 1 descoped

Executes Ryan's caching ruling (`phase_15_caching_and_objective_prompt.md` §1)
plus two decisions taken in the same session.

**Audit of what actually existed** (verified against code, not docs):

| Where | Cached? | Defect found |
| --- | --- | --- |
| Phase 2 tuning (`rotor_hover_solver_phase2b_nearfield_cache.jl:195`) | yes | ran on stock `tune_fmm` — the accuracy-only cost-model objective, i.e. the very defect fixed for the uncached path on 2026-08-24. Build cost was already excluded and recorded, so that half of the ruling was already policy. |
| Phase 2 benchmarks (`rotor_hover_solver_phase2.jl`) | yes | cache build folded INTO `t_solve_min`, violating the ruling's second half. Its justifying comment ("per-SOLVE state, a rotating rotor invalidates it every step") was STALE — `persistent_plan` landed afterwards with rigid-motion support, and the Dirichlet operator + cache blocks are exactly invariant under rigid motion. |
| Phase 1 | no, nothing | by design |

**The finding that reshaped the work.** An uncapped cached tune walks to a
DENSE operator. At R1 the cached optimum sat at leaf 275–342 / 687–862 MB
against a dense 8N² of 514 MB (1.3–1.7× dense) with `cache_capped=false` — the
caps never bound, the cost optimum genuinely wanted it. That is structural: if
assembly is free (and Ryan's ruling makes it free for tuning), matrix-ful always
beats matrix-free wherever memory allows. So `max_bytes` does not guard the
answer, it DETERMINES it. **Ryan's ruling: sweep the memory budget as an axis**,
publishing per-step cost vs memory with `backslash_ldiv` (8N²) and matrix-free
Krylov as the two endpoints of one curve. Budget 0 is how the uncached configs
get tuned — same code path, no second driver.

**W0 — feasibility map** (`benchmark/nfcache_feasibility.jl`, new). Sweeps
`estimate_nearfield_cache` over R1–R7 × a leaf ladder without building or
solving; `bytes` is exact size-pass arithmetic, hardware-independent, so it is
valid locally. Anchor check: leaf 340 → 0.803 GiB reproduces the recorded
862104272 B at leaf 342 exactly. Full map in
`benchmark/results/phase2/multi/nfcache_feasibility.csv`; the headline numbers
are in `ledger.md`. Two results matter:
- cache/dense at the kernel leaf FALLS with N — 0.53 (R1), 0.34 (R2), 0.26
  (R3), 0.157 (R4), 0.114 (R5), 0.082 (R6), 0.068 (R7) — so the dense
  degeneracy is an R1–R2 artifact and matrix-free character genuinely returns
  at scale;
- but the ABSOLUTE floor is brutal. A cache has a floor (shrinking the leaf
  raises block count without shrinking total near-field area): the cheapest
  cache on the ladder is 0.20 GiB at R1 but **78 GiB at R7** and 24 GiB at R6.
  The handoff's "naive R7 ≈ 32 GiB" underestimated by **2.4-2.8×** (78 GiB at
  the smallest leaf, 89 GiB at the kernel leaf 21). A budget below a
  rung's floor admits NO cached configuration at all, and the drivers now
  report that as the infeasible end of the memory axis rather than as a
  descent failure.
- The cache asymptotes at ~2× dense, not 1×: the direct list holds both (i,j)
  and (j,i), so symmetric blocks are stored twice. A dense solver is strictly
  more memory-efficient at the dense end — a publishable point on the memory
  axis.

**W1 — `tune_fmm_perturb` gains `cost=`** (FastMultipole `src/autotune.jl`). A
caller-supplied objective, called as `cost(; expansion_order,
multipole_acceptance, leaf_size_source)` returning a time / `(t, success)` /
`(; t, success)`. The closure owns the ENTIRE cost model, so
`tree_amortization`, `error_tolerance` and `fmm!` kwargs are REFUSED alongside
it rather than silently ignored. Everything else — min-of-reps, abandonment,
memoization, the descent — is unchanged. This is Ryan's option (a) for the
objective blocker, and it also subsumes `tree_amortization` (a real solve
builds one plan and reuses it, which IS the correct amortization) and makes
cached-vs-uncached tuning fall out for free with zero cache plumbing in the
tuner. Verified by a synthetic-surface smoke: descent reaches the analytic
minimum, memoization evaluates nothing twice, infeasible start is rejected
loudly, all three mutually-exclusive knobs throw, abandonment still fires, all
three return shapes work, a bad return type throws.

**W2 — `benchmark/rotor_hover_solver_phase2_tune.jl`** (new). Per rung × memory
budget: estimate bytes first (an over-budget candidate costs one tree build,
not one solve); build a `KrylovSolver` at the candidate knobs with
`cache_nearfield`+`persistent_plan`; run ONE UNTIMED warm-up solve, which is
what lazily builds the plan and its blocks; then min-of-reps timed solves from
a cold SOLUTION (`reset_cold!`, warmstart off, so `niter` is honest) against a
WARM cache. **The build is excluded from the objective by construction, not by
subtraction.** Success is the certified BC verdict on the converged solution —
strictly stronger than the per-apply tolerance the proxy used. Writes
`tune_phase2.csv` (new filename on purpose: the old `tune_cached.csv` holds
accuracy-only rows and must not be mixed).

R1 local smoke, both arms green: budget 0 (uncached) → 113.2 s, niter 59, BC
9.67e-7 certified; budget 2 GiB → warm 4.2–5.6 s at leaf 21, cache 0.254 GiB,
**build 12.14 s**. The build alone is ~2× a warm solve, so the old driver's
folding it into `t_solve` inflated the reported per-step cost roughly 3× —
Ryan's ruling made concrete. (Local timings are indicative ONLY; the 4.219 vs
5.591 s spread between the descent and the re-measure at identical knobs is
the documented 22–39% Mac scatter.)

**W3 — `rotor_hover_solver_phase2.jl` reworked.** `*_nfcache` configs now run
`persistent_plan=true` with the plan+cache built in the SETUP region, so
`t_solve` is the warm per-step cost and the build is its own column. Schema
gains `mem_budget_gib, nfcache_build_time, nfcache_bytes, nfcache_capped`
(header-guarded); resume identity gains the budget, since a cached config now
emits one row per budget. Apply knobs now come from `tune_phase2.csv` (Phase 2
owns tuning), falling back to `TUNED` with a loud warning. Added a D4
additivity check (R1 only, `row_kind="additivity"`): verifies
`t_cold ≈ t_build + t_warm` once so every other rung can rely on it instead of
paying for a second timed pass — Ryan: "a cold solve would merely add the cache
build time, right?" **VERIFIED at R1 (2026-08-24):** `t_cold = 16.059 s` vs
`t_build + t_warm = 11.800 + 4.223 = 16.024 s`, a **0.22% discrepancy**. The
warm test alone therefore recovers both quantities, and no rung needs a second
timed pass. The published warm row reads `t_solve_min = 4.22 s` with
`nfcache_build_time = 12.27 s` in its own column — under the old convention
that same row would have read ~17 s, i.e. the build was inflating the reported
per-step cost by ~2.9x.

**W4 — Phase 1 descoped** (Ryan): its remaining purpose is to solve the ladder
up to R7 and verify that every solver agrees on the SOLUTION. Tuning moved to
Phase 2. `rotor_hover_solver_phase1_agreement.jl` gains the full R1–R7 ladder
(was R1–R3), a knob fallback (accuracy matters here, speed does not), and
`NFCACHE=<GiB>` opportunistic caching that is skipped with a notice when it
does not fit. **Blocker surfaced:** the default reference `backslash_ref`
assembles a DENSE G — 25 GiB at R4, 87 GiB at R5, 335 GiB at R6, **1.31 TiB at
R7** — so an exact reference is impossible above ~R5. Added `REFERENCE=<config>`
plus a guard that refuses a dense reference rather than OOM-ing hours in. Which
config becomes the reference at R6/R7 is **Ryan's call** (see phase_15 §Open).

**W5 — two measurement biases fixed.**
- `rotor_hover_solver_phase2b_nearfield_cache.jl`: the end-to-end A/B had no
  warm-up rep, so the first arm absorbed all the JIT — and the first arm is the
  UNCACHED one, biasing the A/B in favour of the feature under test. The
  previously reported **5.5× cold-solve speedup is therefore an upper bound**.
- `fm_leaf_ab.jl`: corrected diagnosis. Its `bc_rel_l2` column was invalid not
  because of `bc_error!` arguments but because its local `reset!` omitted the
  frozen-potential restore that `reset_cold!` does; the Dirichlet rhs is
  `-rotor.potential` and a solve clobbers it, so every solve after the first
  ran against a stale b. Now uses `reset_cold!`. Timings were unaffected (same
  operator, niter stable to ±1), so the R2 2.18× leaf result stands.

Also: `phase1_case.jl` gains `SKIP_B=1` for geometry-only consumers (the O(N²)
assembly is tens of minutes at R7 and the feasibility map never solves); it
sets `rms_b=NaN` so anything that tries to solve fails loudly.

**Not done / next:** nothing has run on the cluster. Every timing above is
local and therefore indicative only. Next is the R1–R5 real-solve tune on an
exclusive node, then Phase 1 R7 verification, then the Phase 2 tables. R6–R7
objective and the R6/R7 reference config both still need Ryan.

### 2026-08-24 (latest) — the tuner's stall was BIAS, not noise: tree build was charged to every apply

**Revises the "reps is the fix" diagnosis in the entry below**, which is now
half the story. Found while planning the `reps` change; the code fix landed the
same day.

**The defect.** `FastMultipole.tune_fmm_perturb` timed every candidate through
`fmm!(targets, sources, cache; ...)` — the `Cache` path, which **rebuilds both
octrees and the interaction lists inside every timed call** (`src/fmm.jl`,
`fmm!(::Tuple, ::Tuple, ::Cache)`). A steady iterative solve does not pay that
per apply: FLOWPanel builds one `FmmPlan` and reuses its trees/lists across
every Krylov iteration (`src/FLOWPanel_fmm.jl`, the `plan_slot` branch), so one
build is amortized over ~57 applies.

**Why it biases toward large leaves.** Tree and interaction-list construction
get MORE expensive as `leaf_size` shrinks. Charging a full build to every apply
therefore adds a leaf-dependent penalty pointing in exactly the direction of
the observed stall.

**Why it is not just noise.** At R1 the descent moved 68 -> 45 (t=0.444 s), then
measured leaf 30 at 0.450 s and stopped. The leaf sweep implies leaf 30 should
be ~15-20% FASTER than 45. Min-of-2 timing jitter does not plausibly erase a
20% gap; a leaf-dependent additive overhead does. Independent confirmation from
a 20k-body smoke test after the fix: at the same knobs the tree-inclusive
objective reads 0.026 s and the amortized one 0.018 s — the build is ~30% of the
measured cost, and it is precisely the part that grows as leaf shrinks.

**Ryan's ruling: this is a switch, not a replacement.** An unsteady run whose
geometry moves every step DOES rebuild the tree per apply, so the old objective
is correct there (e.g. 023's mature-wake tuning).

**What changed.**

`FastMultipole/src/autotune.jl`, `tune_fmm_perturb` — three additive keywords,
defaults reproduce the old behaviour exactly:

- `tree_amortization::Int=1` — how many applies share one tree/list build.
  Candidate cost = `t_build/tree_amortization + t_apply`. `1` keeps the legacy
  `Cache` path verbatim; `n>1` builds one `FmmPlan` per candidate (timed
  separately, released with a `GC.gc()` after), then min-of-`reps` over applies
  that reuse it. Caller `kwargs` are split between `FmmPlan` and the per-apply
  `fmm!` via a new `FMMPLAN_STRUCTURAL_KWARGS` list in `src/fmm.jl` (`FmmPlan`
  has no catch-all `optargs...`, so an unrecognized key would error).
- `max_seconds=Inf` — wall-clock guard checked between candidates; on expiry it
  returns best-so-far, `@warn`s "max_seconds EXPIRED", and reports
  `timed_out=true`.
- Third return value `info = (timed_out, t_elapsed, n_candidates, t_best)`.
  Existing two-value destructuring still works.

`benchmark/rotor_hover_solver_phase1_tune_hpc.jl`:

- `TUNE_REPS` — default **5 at R1-R4, 2 at R5-R7**, asserted <= 5 (Ryan's cap).
  Not 1: min-of-2 is what discards first-touch/page-fault cost on the first call
  at a new leaf. (`expansion_order` is a runtime `Int`, not a type parameter, so
  there is no per-candidate recompilation to worry about.)
- `TUNE_TREE_AMORTIZATION` — default **50**. Phase 1 is a frozen single-step
  solve that never sheds a wake, and its `KrylovSolver` runs ~50-60 applies
  against one `FmmPlan` (R1 measured ~57 iterations).
- `TUNE_MAX_SECONDS` — per-rung defaults 900/1800/3600/7200/14400/21600/43200 s
  for R1-R7.
- `tune.csv` gains `tune_reps,tree_amortization,tune_timed_out`, plus a header
  guard that refuses to append to a pre-fix schema (the driver APPENDS and
  `phase1_case.jl` reads knobs back by column name, so mixing objectives would
  be silent). Purge a rung's CSVs and re-run; keep `bcache_R*.bin`.

`benchmark/phase1_case.jl`: the hardcoded R1-R3 `TUNED` fallbacks are now marked
as pre-fix accuracy-only picks, not valid campaign knobs.

`benchmark/p023_fmm_tune.jl`: left at `tree_amortization=1` deliberately, with a
comment saying why (mature-wake unsteady step rebuilds the tree every step).

**Ryan's ruling on caching (2026-08-24, latest): cache in BOTH tuning and
benchmarks for Phase 2.** Ignore the cache build cost while tuning (legitimate:
panels-on-self has frozen geometry, so the cache is built once a priori);
MEASURE it for benchmarks but keep it in its own column, never folded into the
per-step solve cost. Only panels-on-self is worth caching — wake->body geometry
moves every step, so those matrices change every step and are never reused
inside a solve (the campaign already freezes `b` via `bcache_R*.bin`, which is
the same fact). Ryan has NOT ruled on whether Phase 1 should also go cached.

Gap this exposes: phase2b's cached tune (`TUNE_CACHED=1`,
`rotor_hover_solver_phase2b_nearfield_cache.jl:195`) already excludes build cost
— but drives it with stock `tune_fmm`, i.e. the accuracy-only objective we just
fixed for the uncached path. `tune_fmm_perturb` has no cache support at all.
Feasibility unmeasured: cache bytes ~ panels x leaf, R4 ~ 4.5 GiB recorded,
naive scaling puts R7 near ~32 GiB (more once caching pushes leaf up), against
a 32 GiB cap and a SERIAL build that `max_build_time` may bind first. Use
`FastMultipole.estimate_nearfield_cache` (no build) to map where it caps.

**BLOCKER: the tuner's objective still disagrees with a real solve.** Under
`tree_amortization=Inf` at R1 (job 13447582, p=17/MAC=0.5, reps=5) the descent
measured leaf 43 = 0.434 s and leaf 29 = 0.435 s — flat — and stopped after one
move, while real solves at the same knobs fall monotonically to leaf 9
(0.272 s/iter at 9 vs 0.458 at 40). The tuner minimizes ONE apply; Phase 1
publishes solve wall clock, and below leaf ~40 the proxy stops tracking.
RULED OUT: `tune=true` counting overhead (<=3%, only at leaf 9) and
per-candidate plan allocation / NUMA first-touch (4 fresh plans scatter the same
as 1 plan re-measured). Residual is ~+-9% against the solve-implied curve —
enough to erase the true ~15% gap between leaf 43 and 29 at `improve_tol=0.02`.
Unresolved whether that is resolution or a real apply-vs-solve divergence.

**INSTRUMENT WARNING.** Ryan's laptop cannot measure any of this: min-of-5 apply
timings scatter **22-39%** at fixed knobs locally, against a 15% effect. All
timing evidence must come from an exclusive cluster node. An earlier local
`tree_amortization` 1/50/Inf comparison logged in this session sat inside that
noise band — the `Inf` code path is verified to RUN correctly, its effect on the
answer is NOT verified.

Handoff for the next agent: `phase_15_caching_and_objective_prompt.md`.

**CONFIRMED: the PUBLISHED ladder carries the same ~2x leaf error
(job 13443627, 2026-08-24).** This was inferred from R1 and never measured. It
is now measured. `benchmark/fm_leaf_ab.jl` at R2 (15,760 panels), fixed
p=17/MAC=0.5, on the current FM:

| leaf | 9 | 21 | 40 | 60 | 71 |
| --- | --- | --- | --- | --- | --- |
| `t_solve_min` [s] | **36.71** | 46.72 | 54.58 | 73.17 | 80.20 |
| niter | 62 | 63 | 62 | 62 | 60 |

Monotone, **2.18x across the range at essentially constant iteration count**.
R2's old campaign pick was leaf **59**, i.e. sitting at ~73 s where ~37 s was
available — the same ~2x error R1 had, in a rung whose numbers were already
published. By the same argument R3 (old pick leaf 58) is affected too.

**Consequence for the published Phase 1 comparison.** The apply knobs are shared
across configs within a rung (ruling 3), so the inflation hits every FMM-based
config at R2/R3 roughly alike — but `backslash_ldiv` runs no FMM in its solve
and is untouched. The published rungs therefore biased the comparison AGAINST
the iterative FMM solvers by up to ~2x. Any conclusion drawn from those rows
about backslash-vs-Krylov crossover needs re-reading once the retuned ladder
lands.

(`fm_leaf_ab.jl`'s `bc_rel_l2` column read ~1.02-1.09 here; that column is known
INVALID — wrong args to `bc_error!` — and only the timings are sound. Accuracy
parity rests on the campaign rows.)

**Ryan's ruling: ignore the tree cost entirely for panel tuning
(2026-08-24, latest).** `TUNE_TREE_AMORTIZATION` default is now **`Inf`**, not
50. The panels-on-panels operator has FROZEN geometry: its tree is built once a
priori and reused across every solve iteration *and* every timestep, so the
build is a one-off that must not influence the choice of apply knobs at all.
The right objective for panel tuning is pure apply cost.

`tree_amortization=1` is retained and is the value to use when tuning a
PARTICLE field — that tree must be rebuilt every timestep, so its build IS paid
per apply and belongs in the objective. A finite `n` prices a build reused over
exactly `n` applies. `Inf` makes `t_build/tree_amortization` exactly `0.0`; the
plan is still built per candidate (the applies need it) but never charged.

Verified on a 20k-body descent at `reps=5`: `t_best` 0.0346 / 0.0287 / 0.0264 s
for `tree_amortization` 1 / 50 / Inf, argmin moving leaf 33 -> 22 -> 22 as the
build charge drops — monotone in the expected direction.

**R1's first validation run (job 13443626, tree_amortization=50) is therefore
superseded.** For the record it seeded at p=14/MAC=0.4/leaf=58 (0.658 s) and
descended 58 -> 39 -> 26 -> 17 to p=14/MAC=0.45/leaf=17 at 0.312 s — the
large-leaf pull is gone, which was the point. It did NOT reach leaf 9, for two
reasons visible in the trace: stock `tune_fmm` seeded a different (p, MAC)
family this time (MAC 0.4/p 14 vs the earlier MAC 0.5/p 17 — its cost model
rated them 0.662 vs 0.679, a 2.6% coin flip that picks the basin the local
descent explores), and within that family the bottom is genuinely flat (leaf 17
= 0.318 s vs leaf 11 = 0.315 s, 0.9%, below `improve_tol`). Leaf numbers are not
comparable across (p, MAC) families.

**Early abandonment (Ryan 2026-08-24, later).** `abandon_factor=1.3` stops a
trial as soon as its running min passes 1.3x the fastest COMPLETE,
error-satisfying candidate so far. Such a point can no longer be accepted —
acceptance needs `t < t0*(1-improve_tol)` and the incumbent `t0` is never below
that best — so the remaining reps are waste. This is what makes `reps=5`
affordable. Under `tree_amortization>1` the applies are skipped entirely when
the amortized build alone already loses. The threshold only tightens, so an
abandoned point stays rejected and its over-estimated time is safe to memoize;
abandoned candidates never tighten the threshold themselves, and one abandoned
at the build stage has NOT had its error tolerance verified (flagged
`abandoned=true` in `history`, counted in `info.n_abandoned`).

Why 1.3 and not tighter: the risk is a candidate whose first rep is unluckily
slow but whose true min is competitive. 1.3 kills a point only when even its
best rep is 30% above the best — well outside jitter. Anything near 1.0 would be
exposed to exactly the sampling noise this campaign is trying to remove.
`TUNE_ABANDON_FACTOR=Inf` disables it.

Measured on a 20k-body descent at `reps=5`: `abandon_factor` Inf / 1.3 / 1.05
all reach the identical argmin `(leaf 22, P 5, MAC 0.75)` and the same
`t_best=0.018 s`, at **9.0 / 7.8 / 6.4 s** wall with 0 / 8 / 20 of 44 candidates
abandoned. The saving is bounded below by one full trial per candidate, so it
grows with `reps` and with the cost spread across the neighborhood — both larger
on the real rungs than on this toy.

**Verification so far.** `FastMultipole` loads; `test/fmm_plan_test.jl` passes
(11/11); a 20k-body gravitational descent runs to convergence under both
`tree_amortization=1` and `=50` and reaches the same argmin there; the
`max_seconds` backstop fires, warns, and reports `timed_out=true`. **Not yet
validated at R1 on the cluster** — that is the next step, and the pass criterion
is a walk 68 -> 45 -> 30 -> 20 -> 13 -> 9 (or lower) with `krylov_gmres/standard`
`t_solve_min` near 15.5 s.

### 2026-08-24 — tuner objective was wrong; cost-based tuning added, R1 re-measuring, R7 tuning

**The apparent FastMultipole regression was not one.** R1 under the new FM came
out 1.5-2.5x slower on every FMM path at identical iteration counts and BC
agreeing to 4 figures, while `backslash_ldiv` (no FMM in its solve) was flat at
-4%. A leaf sweep on the current FM (job 13407693, p=17/mac=0.5) settled it:
leaf 9 = 15.51 s, 21 = 20.65, 40 = 26.10, 60 = 32.64, 71 = 34.95 — monotone,
2.25x across the range at constant niter. **At matched knobs (leaf 9) the two FM
generations agree to +6%** (14.59 -> 15.51 s). The remaining ~140% was the knob.

**Root cause.** Stock `FastMultipole.tune_fmm` optimizes a SINGLE evaluation
against an error tolerance; it never measures iterative-solve cost, so it accepts
a large leaf that is accurate but ~2x more expensive per apply — exactly the
quantity Phase 1 measures.

**Correction to an earlier claim in this log.** I first wrote that the old tuner
picked leaf 9 and the new one picks 60, i.e. that the tuner had regressed. That is
true only at R1. The old campaign's picks were 9 / 59 / 58 / 32 / 25 / 20 across
R1-R6 — the same large-leaf behaviour. R1's 9 was a lucky pick, not a better
tuner. **Implication: parts of the OLD accepted ladder may also be mis-tuned, R2
and R3 most suspiciously.** Inferred from R1's sweep, not measured at R2; the
cheap test is `benchmark/fm_leaf_ab.jl` with `RUNG=R2` against its old leaf 59.

**Fix.** `rotor_hover_solver_phase1_tune_hpc.jl` now runs stock `tune_fmm` for an
accuracy-valid seed, then `FastMultipole.tune_fmm_perturb` (`src/autotune.jl:337`)
to coordinate-descend on BENCHMARKED wall-clock under the same error tolerance.
`PERTURB_TUNE=0` restores the old behaviour (variant `tuned_seed_only`).

**Cost tuner validated at R1 and found INSUFFICIENT (2026-08-24, later).**
`tune_fmm_perturb` seeded from stock `tune_fmm` moved R1 from leaf 68 to leaf 45
(~35 -> ~28 s) but stopped there, ~1.8x short of the ~15.5 s available at leaf 9.
Trace: it probed BOTH directions correctly (leaf 102 upward = 0.943 s rejected,
leaf 45 downward = 0.444 s accepted), then at iter 2 measured leaf 30 at 0.450 s
— 1.4% WORSE than 45 — and terminated. That is a measured regression, not a
below-threshold gain, so `improve_tol` is not the lever; `reps=2` simply cannot
resolve the real difference (the sweep predicts ~0.40 s at leaf 30).
**Fix agreed with Ryan: `reps=2 -> 7`, plus a `TUNE_MAX_SECONDS` wall-clock guard
and possibly rung-scaled reps.** A custom leaf ladder was considered and rejected
— `perturb` already does bidirectional coordinate descent over (p, mac, leaf) and
correctly rejects candidates on accuracy; the only defect is the noise floor.
Note the sweep never sampled below leaf 9, so the true optimum may be lower; FMM
cost is U-shaped in leaf and R1/p=17 sits on the "too large" branch throughout.

**Ryan's ruling: retune EVERY rung, not just R7** — the old knobs came from the
same accuracy-only objective (leaf 9/59/58/32/25/20 across R1-R6), so reusing any
of them is unsafe.

**Campaign state.** All jobs cancelled and all results CSVs purged (the seven
`bcache_R*.bin` kept — FM- and tuning-independent). Next agent picks up from
`phase_14_cost_tuner_relaunch_prompt.md`: adjust the tuner in plan mode, retune
all rungs, then relaunch. Earlier state: all earlier-generation jobs cancelled; all bad-knob CSVs
purged (the seven `bcache_R*.bin` kept — direct N^2 assembly, FM-independent).
Two validation chains launched: R1 tune+tables (13425064/65/67) as the check
against the known-good leaf 9, and R7 tune s1-s3 + tables (13425069-77) since R7
has never been tuned. R2-R6 wait on R1 validating. Ryan's ruling: R7 does not
gate the others.

**Trap recorded.** 21 held jobs had to be cancelled: `afterany` treats a CANCELLED
upstream as satisfied, so releasing them would have run tables with no knobs.
Tune->table chains now use `afterok`. Second trap: the new resume support skips
any row already in the CSV, so stale bad rows are silently inherited — purge a
rung's CSVs before re-running it.

Handoff for the next agent: `phase_14_cost_tuner_relaunch_prompt.md`.

### 2026-08-24 (late) — campaign restarted: new FMM generation, isolated checkout, half-node allocation

The morning's relaunch (13396780-13396817) was cancelled before it ran. Reason:
the source it would have used could not have worked. The 2026-08-24 FLOWPanel
references nine FastMultipole symbols absent from the pinned `a9b734a` —
including `NEARFIELD_CACHE_DEFAULT_MAX_BYTES` at `src/FLOWPanel_solver.jl:504,644`,
i.e. the solver every 021 driver runs. My earlier claim that the mixed generation
was "deliberate and, for 021, the correct state" was WRONG; the `import FLOWPanel`
load test passed only because those symbols appear as default argument values,
evaluated at call time rather than at load. Loading proved nothing about running.

**Provenance was the deeper problem.** The old pin `fm_commit = a9b734ad...-dirty`
was an uncommitted worktree whose `src/autotune.jl` blob (`47beaa59`) exists in no
git history anywhere (checked with `git cat-file` and `--find-object`). When new
FastMultipole content was first placed on the cluster it was laid on top of that
same commit, so two entirely different FMM generations would both have recorded
the identical string `a9b734ad...-dirty` — a validated CSV column turned actively
misleading. Ryan then committed and pushed FastMultipole, FLOWPanel and FLOWVPM,
and the campaign moved to clean clones.

**New home: `~/flowpanel-021/`** — isolated from `~/projects/`, which is shared
with 018/023 and where two agents collided earlier today. Layout is
`FLOWPanel.jl` + sibling `FastMultipole` (symlink to the `.jl` clone; the Manifest
dev path is `../FastMultipole`) + `FLOWVPM.jl`. All three clones zero-dirty:
FLOWPanel `321473f`, FastMultipole `d8258a7d`, FLOWVPM `4494c25`. Two gaps had to
be filled by rsync because they are untracked upstream: `benchmark/` (only 7 of
its files are tracked) and `Manifest.toml`. 47/47 md5-verified.

`Project.toml` was left ALONE. `PythonPlot` in it is a `[weakdeps]` entry, inert
unless the extension triggers — instantiate ran clean and created no `.CondaPkg`.
The 08-22 killer was `GeoIO`, a hard dep, and it is absent from the commit.

**Knobs are re-tuned, not reused**: `tune_fmm` tunes FMM parameters, so its output
is FM-dependent by construction. All of R1-R7 re-tune.

**Half-node allocation (Ryan ruling).** 64 CPUs per job, not an exclusive 128.
Julia runs 64 threads either way, so `--exclusive` wasted half of every node and
made each job wait for a fully free one. At 64, two jobs pack per node: 24 of 32
started immediately instead of queuing ~15 h. Trade-off recorded in `ledger.md`
and in every launcher header — `OverSubscribe=YES:4` means a co-tenant can now
contend for memory bandwidth, which is a real if small contaminant on a
bandwidth-bound FMM/ILU workload.

**Two self-inflicted errors, both caught before any row was written.**
(1) `--ntasks=64 --cpus-per-task=64` requests 4096 CPUs; the correct spec for one
64-thread process is `--nodes=1 --ntasks=1 --cpus-per-task=64`. (2) The runtime
CPU assert used `nproc`, which honours `OMP_NUM_THREADS` — and Slurm presets that
to 1 whenever `--cpus-per-task` is used, so it read 1 CPU while the real affinity
was the full `0-63`. That false failure killed all 32 jobs in ~20 s each. Net
positive: the assert exists precisely to stop a silent 64:1 thread
oversubscription from inflating every timing, and the fix (assert on
`SLURM_CPUS_ON_NODE` for the allocation, `getconf _NPROCESSORS_ONLN` for the node)
is stronger than what it replaced.

**Resume support added**, required because `standby` is preemptible and the
launchers now set `--requeue`. `rotor_hover_solver_phase1_table.jl` and
`..._phase2.jl` skip on `(rung, config, row_kind)`; `..._phase2b_nearfield_cache.jl`
skips on `(rung, row_kind)`, driving its existing `SKIP_AB`/`TUNE_CACHED` gates.
Rows flush as each measurement completes, so a row in the CSV means that
measurement is finished and is never repeated. Parser validated against real
data: recovers R6's 2 partial rows and R1's complete 7 including the `fgs` dual
row.

**Walltimes right-sized** from the measured 08-18/19 timings with ~2x margin
instead of a blanket 7 days — on a priority-0 standby QOS an oversized request
backfills far worse than an accurate one.

**Launched: 32 jobs.** Tune R1-R6 plus R7 as three `afterany` stages; tables
R1-R7 x multi/single with the ledger's cost-ceiling drops (`backslash_ldiv` out
at R6-R7 both modes, `krylov_gmres` additionally out of single at R6-R7); Phase 2a
R1-R5 multi (the `*_nfcache` variants wait on 2b's `tune_cached.csv`); Phase 2b
R1-R4. Tables chain `afterany` on their own rung's tune. First R1 tune row landed
clean under the new FMM: `p=17 MAC=0.5 leaf=71`, `meets_target=true`,
`fm_commit=d8258a7d...` with no `-dirty` suffix.

**Note for whoever compares generations:** R1's `rms(b)` moved from
`0.011505099103413874` to `0.011505099103103312` — agreement to ~10 significant
figures, differing in the 11th. b is a direct N^2 assembly, so this is the new
FLOWPanel source, not the FMM; small, but it confirms the break is real rather
than nominal.

### 2026-08-24 — all ten 08-22 jobs diagnosed (one cause), repaired, 13 relaunched

**Ryan's report was "all HPC jobs failed"; the cause was single and load-time.**
`sacct`: eight `FAILED` exit `1:0` after 52 s - 3 min 17 s; `13306604`/`13306605`
`CANCELLED` at `00:00:00` as `afterok` casualties of `13306603` — chain
casualties, not independent failures, the same confusion that cost a wrong
diagnosis on 08-22. Every `.err` on both checkouts:

```
ERROR: LoadError: ArgumentError: Package GeoIO not found in current path.
  [6] top-level scope @ ~/projects/.../benchmark/phase1_case.jl:27
in expression starting at ~/projects/.../examples/dji9443_trailing_edge.jl:2
```

**Root cause — two individually-correct decisions that were jointly fatal.** On
08-22 GeoIO was dropped from `Project.toml` on both cluster checkouts to close
the `.CondaPkg/lock` EACCES class. Locally `examples/dji9443_trailing_edge.jl`
lost its `import GeoIO` in the same generation (it reads meshes through the
native `pnl.read_gmsh`; see `src/FLOWPanel_utils.jl:112`). But the standing sync
ruling covers `src/` + `benchmark/` only, so `examples/` was never pushed. Both
cluster copies kept `import GeoIO` on line 2 (md5 `0f8dce88…` on both vs local
`d75e8df3…`), and `benchmark/phase1_case.jl:27` includes that file — so all five
drivers died before touching solver code. An md5 sweep showed **45+ files under
`examples/` were a full generation stale**; the TE detector was simply the one
the 021 drivers reach.

**Lesson for the catalogue: `examples/` is a real dependency of every Phase-1
driver.** The sync ruling must be `src/` + `benchmark/` + `examples/`.

**Repair (Ryan's rulings 2026-08-24).** Whole-tree rsync of `src/` +
`benchmark/` + `examples/*.jl` to both checkouts, no `--delete`, excluding
`Project.toml` / `Manifest.toml` / `examples/data/` **and `benchmark/results/`**
— the last of these was nearly missed and would have clobbered the cluster's
authoritative accepted data (561 files, per-rung R1–R6) with the thinner local
tree. Verified **0 md5 mismatches over 236 files on both checkouts**; the only
unmatched entries are two intentional cluster-only files preserved by omitting
`--delete`. `import FLOWPanel` plus the exact failing include now load clean on
both. Not committed (Ryan's call).

**Hardware consistency audit (Ryan's request, mid-session).** Every job —
failed and accepted alike (13206077, 13242659/61/63/65) — already ran
`m12` / `zen3` / `AllocCPUS=128` / exclusive, so the existing ladder *is*
mutually comparable. Only `ReqMem` drifted (64/128/192/256 G), and only because
the flags were typed by hand: `p1_table.sh`'s documented submit line omitted
`--constraint`/`--exclusive` entirely. Now pinned identically into all seven 021
SBATCH headers — `--partition=m12 --constraint=zen3 --exclusive --mem=0` — plus a
runtime assert that exits 1 if `SLURM_CPUS_ON_NODE != 128`, and a stable
`HARDWARE_TAG=orc-m12-zen3-2x64c-512g`. Node spec measured by `scontrol`:
Sockets=2, CoresPerSocket=64, ThreadsPerCore=1 ⇒ CPUTot=128, RealMemory=524288 MB,
MaxMemPerNode=524288. `--mem=0` retires the per-rung memory knob entirely and
clears both the R5 dense-G (~94 GiB) and the 32 GiB near-field cache cap.
Pre-existing rows keep their old tags (`m12-2-30`, `zen3-exclusive`); all three
labels denote the same verified hardware.

**Chain semantics.** `p1_tune.sh` now documents `afterany`, not `afterok`, and
carries a stage-input assert that names the missing CSV and exits 1 — so a
downstream stage produces a diagnosable failure instead of vanishing.

**Append-duplication swept before relaunch.** The failed jobs wrote nothing, but
older partial output did: R6-multi already holds accepted `krylov_gmres` +
`krylov_jacobi` rows (08-20, from 13206077 before it died on the ILU limit), so
the R6 relaunch runs only `krylov_ilu:fgs:fgmres_fgs`. Phase 2b's R3 partial (3
rows, no `tune_cached` — it died on the old 4 GiB/60 s caps) and R4 header-only
file were quarantined to `*_20260820.csv`.

**Relaunched 13 jobs.** Phase 1: `13396780` R6-multi (3 configs), `13396781`
R6-single-a, `13396782` R6-single-jacobi, `13396783/784/785` R7 tune s1→s2→s3
chained `afterany`. Phase 2a (R2–R5 multi, 6 configs — the `*_nfcache` variants
need 2b's `tune_cached.csv`, which R2–R5 do not have yet): `13396805/807/809/811`.
Phase 2b (R2–R4, raised caps): `13396813/815/817`. All validated with
`sbatch --test-only` (128 processors, m12) before submission. The eight
BRAINSTORM-018 jobs (`13396725`–`13396732`) are running and were not disturbed.

**Two things left open by this round:** R6-single is a timeout risk (multi-mode
gmres alone took 1364 s/solve at R6; single-threaded on 212k panels against a
hard 72 h ceiling), and the R7 *table* still has no splitting strategy (~100 h vs
72 h). R7's table is downstream of the tune chain, so it did not block.

**Also corrected this session (local, no HPC impact).** The one-body early break
in the tuple `solve!` made four prose sites wrong; see `ledger.md` "Data
quarantine" for the note fixes and the pre-fix row quarantine, and
`phase_12_relaunch_prompt.md` for the corrected "neither path is altered" claim
(the single-body `solve!` and the FGS/Krylov kernels *were* changed; the
no-measured-change conclusion survives on `assemble_source_potential=false`,
which no driver opts into).

### 2026-08-23 — Phase 3 solution-history / block-GS defect FIXED; smoke table replaced

Implemented the approved plan (`plans/validated-scribbling-heron-reviewed.md`).
Uncommitted, local only; no HPC submission.

**What was wrong.** Warm-start persistence lived inside the raw kernels — FGS
`save_solution!` at the end of `_solve!`, and Krylov's `x_prev` / `x_history`
writes at the end of `_krylov_launch!`. `solve!(bodies::Tuple, solvers::Tuple)`
calls a body solve once per block-GS *outer iteration*, so a single timestep
wrote several history slots. Order-1 extrapolation across a duplicated pair then
alternates around the answer, which can keep the outer `max_delta` test from
settling. Separately, a lone body has no coupling to iterate (`sources` is
empty), yet the outer loop still re-solved the same system every iteration — a
standing repeated-solve tax that the campaign had noted since Phase 0 but never
removed.

**The fix — three events, not one.** `src/FLOWPanel_solver.jl` now distinguishes
(1) a raw kernel call, which only updates `solver.niter`/`.solved`; (2) one
per-body solve inside a step (`note_step_solve!`); (3) a completed top-level step
(`save_step_solution!` then `publish_step_stats!`, together
`finalize_step_solution!`). Only (3) advances warm-start history. Step boundaries
sit at the orchestration sites: the two public single-body `solve!` methods (new
internal `finalize_step` kwarg), the tuple `solve!` (one begin/finalize per
solver bracketing the whole outer loop), and `TraceCorrected`'s direct `_solve!`
bypass (finalized after `set_wake_correction!`). Raw `_solve!` callers — Phase
1/2 benchmarks, solver-history tests, `FGSPreconditioner`'s inner solve — stay
outside the accounting by design. A one-body tuple now breaks after the first
block solve, so it records exactly one outer entry. Failure semantics are
explicit: a step that throws leaves history and published statistics untouched.
`SolveStepStats` costs 48 bytes + an 8-byte reference per solver (ruling 8).

**Instrumentation was also wrong, independently.** The driver sampled
`solver.niter` after the step, which returns the *last* inner solve — already
warmed by its predecessors within the same step. `niter_first` (first inner
solve of the step) and `nsolves` (per-body solves in the step) are new CSV
columns; `niter_first` is the warmstart metric, legacy `niter` is a diagnostic.
For these one-body runs the driver and `phase3_analysis.jl` both require
`nsolves == 1` on every scored row, treat `nsolves == -1` as *unavailable*
rather than as a pass, and refuse pre-Phase-3 (narrower) CSV inputs instead of
mixing row widths under one header.

**Voided.** The 2026-08-22 smoke table (`phase_03_history_fix_handoff.md`) —
the "prev 6 → 0 iterations" reduction and the cold 4.378 s / prev 3.872 s /
extrap 81.13 s comparison — is void and marked so in place. Both artifacts are
now explained: last-solve sampling for the iteration column, redundant (and for
extrap runaway) repeated solves for the wall times. Its particle/checksum
columns are not comparable to anything current either, since Gaussian
regularization landed in between (62d72db).

**Post-fix smoke, R1/FGS/local 4T, `N_STEPS=6`** (7 marched steps;
wake-developed window = steps 4–7; `nsolves == 1` on every scored row in all
three arms):

| arm | niter_first mean ± sd | Δiter vs cold | t_solve mean ± sd [s] | Δt vs cold | particles | checksum |
| --- | --- | --- | --- | --- | --- | --- |
| cold | 6.00 ± 0.00 | — | 2.3050 ± 0.0366 | — | 407 | 114595.9408684665 |
| prev (order 0) | 4.75 ± 0.50 | −20.8% | 2.1967 ± 0.0685 | −4.7% | 407 | 114595.9503456424 |
| extrap (order 1) | 5.00 ± 0.00 | −16.7% | 2.2339 ± 0.0365 | −3.1% | 407 | 114595.9683812670 |

All inner solves converged in all three arms. The extrapolation diagnostic
(`SOLVER_VERBOSE=1`, order 1) shows the projected-minus-solved residual (first
strength, column 2) going +1.09e-3 at the first projected step — the history is
still filling there — then holding **one sign** and decaying monotonically
across the wake-developed window: −8.01e-4, −2.83e-4, −2.06e-4, −1.45e-4. The
pre-fix log's alternating ±7.8e-4 pair is gone.

Note the modest, believable size of the warm-start benefit here: at R1 with a
6-sweep cold solve, one or two sweeps saved is what a good guess buys. The
extravagant pre-fix numbers were measurement artifacts in both directions.

**Open for Ryan — trajectory-identity guard threshold.** The `phase3_analysis.jl`
guard flags a warm arm as DIVERGED above 1e-8 relative checksum difference. Both
warm arms match cold on particle count (407) but differ by 7.9e-8 (prev) and
2.4e-7 (extrap) relative, so both are flagged. The arms all converge to the same
solver tolerance (`tol_abs = 1.77e-8`, absolute, on ~8016 strengths) and the
difference amplifies through the wake over steps, so 1e-8 on a *summed* checksum
looks too tight to be the right test. The threshold is untouched pending your
call; the plan's own acceptance bar was "about 1e-7 relative", which prev meets
and extrap misses by ~2×.

**Tests.** New/extended coverage in `test/runtests_unit_solver.jl` (per-step
lifecycle, 51 assertions), `test/runtests_unit_solver_history.jl` (one-body and
two-body tuple step semantics, `solver_optargs` guard),
`test/runtests_unit_fgs_history.jl` (step-boundary writer), and
`test/formulation_test.jl` (TraceCorrected). One documented limitation: the
formulation validator admits only `Backslash` for `TraceCorrected`, which carries
no step statistics, so that hook is currently inert — the test pins the no-op
contract and the ordering, and says what to extend if the validator widens.

### 2026-08-22 (night) — Gaussian rerun executed: scope was 3 artifacts, not the campaign

Executed the pre-Phase-3 rerun. The scope question ("what changes under
Gaussian?") resolved to **two mechanisms**, and knowing which one is active in a
given driver decides everything.

**Channel A — the filament kernel.** Only `_bound_vortex_velocity` /
`_bound_vortex_gradient` (`src/FLOWPanel_elements_fmm.jl:938-1067`) branch on the
family, reached from vortex-ring panel edges (`:862,867`), panel-wake trailing
filaments (`src/FLOWPanel_wake.jl:2891-2897`), and the semi-infinite gradient
helper (`:1627`) — and **only on velocity/gradient passes**. Potential-only
evaluation routes through `ConstantDoublet` and never reaches a filament.
Vortex particles / FLOWVPM are entirely unaffected.

**Channel B — FMM buffer radii (the non-obvious one).**
`radius_inflation(::Type{VortexRing}, core_size, tol)`
(`src/FLOWPanel_elements_fmm.jl:1130-1142`) reads the same global but dispatches
on the body's ELEMENT TYPE, not on whether a filament is evaluated. It feeds the
source-buffer radius (`src/FLOWPanel_abstractbody.jl:1200`,
`src/FLOWPanel_wake.jl:355,2815`), moving the near/far split. So it can change
timings in a run with no wake at all.

**Measured, geometry only (no solve), both ends of the ladder:**

| offset | Δ(buffer) vat−gau | vs smallest panel radius, R1 | vs smallest, R7 |
| --- | --- | --- | --- |
| `core_size_panel = R·1e-10` (all solve + tuning runs) | 3.773e-10 m | 1.7e-5 | 1.5e-4 |
| `core_size_targets = 1e-3` (`steady!` target pass) | 3.171e-2 m | 1413× | 12321× |

(R1 panel radius min 2.24e-5 max 4.43e-3 m; R7 min 2.57e-6 max 5.66e-4 m.)

**This is why the FMM tuning runs do NOT need redoing** (Ryan's question). Every
tuning artifact — `tune.csv`, the FGS τ-ladder, `fgsprecond.csv`, the p2b
cached-economics tune — runs at the panel offset, where the family moves the
buffer by 0.0017% of the smallest panel at R1 and 0.015% at R7. The optimal
(p, MAC, leaf) is family-invariant. Corroborated by an R1 A/B through
`FastMultipoleBackend` at tuned knobs (so the near/far split was genuinely
exercised): rms_b, ‖x‖, Σx bit-identical and the SAME 118 iterations.
CAVEAT, stated honestly: the panel-offset perturbation GROWS ~9× from R1 to R7
as panels shrink, and 1.5e-4 relative is not machine-zero. A marginal MAC
comparison sitting within 0.015% of threshold could in principle flip at R7.
Bit-identity is measured at R1 only; R7 rests on the scaling argument.

#### What was rerun

| # | Artifact | Result |
| --- | --- | --- |
| 1 | `results/phase1/multi/agreement.csv` | **NOT sensitive — prediction refuted by measurement, NO rerun needed.** See the correction below; the file is restored as the record of evidence. |
| 2 | `results/phase2/multi/unsteady.csv` | Rerun R1 backslash+fgs under Gaussian; now 12 rows, all `filament_reg=GaussianRegularization`. |
| 3 | `rigid_motion_tree_reuse/fgs_staleness.csv`, arm `fgs_geometry_order_fixed_8step` | Rerun under Gaussian. **Reproduced BIT-FOR-BIT** — 171/171 comparable cells identical, divergence 1.740e-7 → 1.886e-7 flat, growth 1.1×. |

Result 3 is the useful one: it **resolves the provenance indeterminacy**. The
arm was already Gaussian (as suspected — the in-repo default became Gaussian
mid-day 2026-08-20 and this arm ran 08-21), so the Phase-3 trust gate
(`phase_02...md:313`, "NO Phase-3 unsteady FGS row is trusted until the
discriminator runs") is certified in the campaign's kernel. No further action.

#### CORRECTION — Channel B is a COST channel, not an ACCURACY channel

Predicted, from the 12321× buffer difference at the targets offset, that
`agreement.csv` would move. **It does not.** Rerunning R1 (all 7 configs) under
Gaussian reproduces `CT_laplace` to 8 decimals against rows dated 2026-08-14 —
which are definitively Vatistas, since the kernel was hardcoded and not
selectable until 08-20:

| config | CT_lap gaussian | CT_lap vatistas (08-14) | Δ% |
| --- | --- | --- | --- |
| backslash_ref | 0.06040613 | 0.06040613 | 0.0000 |
| backslash_coupled | 0.06040565 | 0.06040565 | 0.0000 |
| krylov_gmres | 0.06040636 | 0.06040636 | 0.0000 |
| krylov_jacobi | 0.06040614 | 0.06040614 | 0.0000 |
| krylov_ilu | 0.06040613 | 0.06040613 | 0.0000 |
| fgs | 0.06040613 | 0.06040613 | 0.0000 |
| fgmres_fgs | 0.06040607 | 0.06040607 | 0.0000 |

**Why the prediction was wrong.** `radius_inflation` changes WHICH panels are
classified near vs far — it does not change the VALUE the FMM returns, because
both the near and far paths evaluate the same influence to the same tolerance.
So Channel B moves **cost**, never numbers. Only Channel A (the filament kernel
actually evaluated on a velocity/gradient pass) moves numbers, and that requires
filaments to exist — i.e. a wake.

Revised rule, and the one to carry forward:

- **Channel A → numbers**, wake-carrying runs only. Confirmed by the 328 vs 329
  particle split below.
- **Channel B → timings only**, any `VortexRing`-typed body, wake or not.

Consequence: the rerun scope was even smaller than planned — only the
wake-carrying unsteady artifacts ever needed it. R2/R3 of `agreement` were
cancelled once R1 established the null; `agreement.csv` is restored unchanged as
the record, and the Gaussian R1 confirmation is kept as
`agreement_gaussian_verify_R1_20260822.csv`.

Superseded originals are in `quarantine_vatistas_20260822/` beside each artifact,
with a README; they carry no `filament_reg` column and must not be pooled with
post-08-22 rows.

#### Same-code family A/B (new measurement) — `unsteady_family_ab_20260822.csv`

Old-vs-new could NOT isolate the family: the 08-18 rows predate three days of
code changes (and `nsim` is derived from `maneuver!` deltas, so `N_STEPS=5` now
yields 6 recorded steps vs the old 5). So both families were run on TODAY's code:

| metric | gaussian | vatistas | vat−gau |
| --- | --- | --- | --- |
| Σ t_solve, backslash | 26.432 s | 23.723 s | **−10.25%** |
| Σ t_solve, fgs | 30.333 s | 30.851 s | **+1.71%** |
| particles at final step | 328 | 329 | +1 |

**The robust finding is the particle count: 328 vs 329.** That is deterministic
and load-independent, and it proves the family changes the unsteady trajectory —
so 025's CT null (+0.023%) does not license treating wake-carrying 021 work as
family-invariant.

**The timing columns are NOT trustworthy as stated** and are recorded as
indicative only: they are single runs (the campaign's own ruling 5/7 requires
min-of-k, k≥5, after warmup), the two families ran in separate batches so
machine thermal/load state is confounded with the family, and the signs
DISAGREE between configs (−10% backslash vs +2% fgs) — which is what noise
domination looks like, since a pure near/far-cost effect should move both the
same way. A proper measurement belongs on HPC under the min-of-k protocol.

#### Instrumentation fix — `filament_reg` now recorded (SCHEMA AMENDMENT)

Root cause of the whole indeterminacy: nothing recorded the family. Fixed
centrally in `assert_and_banner()` (`benchmark/common.jl`) so every 021 driver
inherits it — banner text, returned NamedTuple, `RUNS_CSV_COLUMNS`, and the
`validate_runs_csv` required-non-empty list. Pattern lifted from
`benchmark/fm051_pass_parity.jl:922` (previously the ONLY recorder in the tree).
Row-construction sites updated in `rotor_hover_solver_phase1.jl`,
`rotor_hover_solver_smoke.jl`, and the three rerun drivers.
`decision_rules.md` amended (the schema is pinned there).
Verified live: the new unsteady rows carry `filament_reg=GaussianRegularization`.

#### Latent bug found, NOT fixed

`_U_boundvortex` (`src/FLOWPanel_elements_fmm.jl:1559-1615`) is **hardcoded
Vatistas** and ignores the global, while its gradient counterpart
`_U_boundvortex_gradient` (`:1627`) honours it. Only fires on the
`semiinfinite_wake=true` path, which 021 does not use — latent, not live.
Separate fix; flagged here so it is not rediscovered.

## ~~OPEN TODO — relabel the 10 in-flight jobs after they finish~~ — CLOSED 2026-08-24

**Moot.** All ten jobs died at load time before reaching any solver code and
produced **zero** rows (`find benchmark/results -newermt 2026-08-23` empty on
both checkouts). There is nothing to relabel. The relaunched generation
(13396780-13396817) runs Gaussian from the launcher default. Original text kept
below for the record.

### Original text (2026-08-22)

The ten jobs queued on 2026-08-22 (13306549/50/54/56/57, 13306599/600,
13306603/604/605) were submitted BEFORE the launcher default was flipped to
Gaussian. Slurm snapshots the batch script at submit time, so they will run and
self-report as `FLOWPANEL_FILAMENT_REG=vatistas`.

**Ryan's ruling: let them run.** Their numbers are family-invariant by
construction (Phase 1 / 2b are potential-only at the panel kernel offset —
measured bit-identical, see the 2026-08-22 evening entry), so the label is a
provenance artifact, not a correctness problem. Cancelling would have cost queue
position on a cluster at 128/136 nodes for zero numerical change.

**ACTION WHEN THEY FINISH — do not skip:** update the family attribution on
their outputs so the campaign reads as uniformly Gaussian.

- [ ] 13306549 p2b-nearfield-R1
- [ ] 13306550 p2b-nearfield-R2
- [ ] 13306599 p2b-nearfield-R3
- [ ] 13306600 p2b-nearfield-R4
- [ ] 13306554 p1-table-R6-multi
- [ ] 13306556 p1-table-R6-single-a
- [ ] 13306557 p1-table-R6-single-jacobi
- [ ] 13306603 p1-tune-R7-s1
- [ ] 13306604 p1-tune-R7-s2
- [ ] 13306605 p1-tune-R7-s3

For each: relabel the `filament_reg` attribution on the produced rows
(vatistas -> gaussian) and annotate in `ledger.md` that the row is
family-invariant by measurement rather than by the pin it ran under. Note these
runs predate the `filament_reg` CSV column reaching the cluster, so the
attribution lives in `banner.txt` / the launcher snapshot, not in the rows
themselves — check `scontrol show job <id>` output or the job's `slurm-*.out`
banner to confirm what each actually ran under before relabelling.

### 2026-08-22 (evening) — Campaign moved to Gaussian regularization; MEASURED inert for Phase 1/2b

**Ryan's ruling (2026-08-22):** the campaign uses Gaussian filament
regularization from here on; the earlier vatistas-era rungs are to be RERUN
under Gaussian later. **TODO carried forward — see "Gaussian rerun backlog"
below.**

**Measured first, before acting.** Slurm snapshots the batch script at submit
time, so editing a launcher cannot reach an already-queued job — applying the
switch to the 10 in-flight jobs would have meant cancelling and resubmitting all
of them. That is only worth doing if the flag actually changes the answer, so it
was tested rather than assumed. R1, GMRES to rtol 1e-12, tuned knobs, 2 threads,
identical in every respect but `FLOWPANEL_FILAMENT_REG`:

| reg | rms_b | norm(x) | sum(x) | niter |
| --- | --- | --- | --- | --- |
| vatistas | 0.011505099103103303 | 6.437344805815382 | 11.41656294891941 | 118 |
| gaussian | 0.011505099103103303 | 6.437344805815382 | 11.41656294891941 | 118 |

Bit-identical to the last digit, same iteration count. **Root cause: the
regularization is unreachable in this case.** `phase1_case.jl` builds a frozen
single-step Dirichlet solve — RHS is a pure body->body direct source-potential
assembly, `semiinfinite_wake=false`, and there is no time-marching, so no wake
is ever shed and no filament kernel is evaluated. Every one of the 10 queued
jobs (p1_tune, p1_table, p2b_nearfield) includes `phase1_case.jl`, so the flag
is inert for all of them.

**Therefore the 10 queued jobs were NOT cancelled** — resubmitting would have
cost a queue cycle and changed no number. The launcher DEFAULTS were flipped to
`gaussian` in all six `benchmark/slurm/*.sh` (md5-verified identical on both
checkouts), so the switch applies from the next submission onward and, more
importantly, to the wake-carrying `p2_unsteady.sh` arm where it is NOT inert.

**Consequence for the ledger:** the existing Phase-1 rows are unaffected by this
switch — they are filament-free by construction, so vatistas-era and Gaussian-era
Phase-1 rows ARE mutually comparable. The regularization only partitions the
evidence for wake-carrying work (Phase 2 unsteady, and campaigns 018/022/025).

#### Gaussian rerun backlog (opened 2026-08-22)

To revisit once the current queue drains. Nothing here is urgent, and per the
measurement above most of it may be a no-op:

- [ ] Phase-2 unsteady arms (`p2_unsteady.sh`) — the ONE place in 021 where the
      regularization is known to be live. Any pre-08-22 unsteady rows were taken
      under vatistas and need rerunning under Gaussian before being pooled with
      post-08-22 rows.
- [ ] Confirm the inertness result holds at a wake-carrying Phase-2 state (the
      08-21 FGS work used a "first particle-carrying R1 state"); the R1 probe
      above only establishes it for the frozen wake-free case.
- [ ] Phase-1 R1-R6 rungs + `tune.csv` / `bcache_R*.bin`: rerun NOT required on
      regularization grounds (filament-free), but note in the ledger that they
      were produced pre-switch so the provenance is explicit.

### 2026-08-22 (later still) — Walltime audit: R7 tune restructured into stages; p2b R3/R4 raised

Before any of the 8 resubmitted jobs started, an elapsed-vs-N fit on the
COMPLETED rungs showed two of them could not finish in the walltime requested.
The fit validates: it predicts R5-single at 30-57 h (mid 38 h) against an actual
44.6 h, so mid-range numbers are good to roughly +-40%.

| pipeline | measured rungs (h) | exponent in N | R6/R7 extrapolation |
| --- | --- | --- | --- |
| p1-tune | R1 0.29, R3 1.29, R4 4.32, R5 10.46, R6 39.05 | 1.18-1.96 | **R7 87-148 h (mid ~114)** |
| p1-table multi | R1 0.08, R2 0.22, R3 0.36, R4 0.87, R5 2.32 | 0.77-1.58 | R6 4-7 h |
| p1-table single (all cfg) | R1 0.81, R2 3.37, R3 6.38, R4 15.37, R5 44.6 (a+gmres) | 1.06-2.12 | R6 ~61 h across the two split jobs |

`n_panels` doubles per rung (R1 8,016 -> R7 419,276) and the `m12` partition
MaxTime is **3-00:00:00**. So R7's tune at ~114 h does not fit a single job at
ANY allowed walltime, let alone the 48 h first requested. This is structural,
not a sizing slip, and it extends to the rest of the R7 row:
`p1-table-R7-single-a` extrapolates to ~100 h and the ledger already had
`R7-single-jacobi` at ~69 h ~= its limit. **The R7 row needs a splitting
strategy before it is submitted, not a bigger `--time`.**

**Actions taken.** Cancelled 13306551/52 (p2b R3/R4) and 13306553 (p1-tune-R7)
while all were still PENDING, so nothing was wasted.

- p2b R3/R4 resubmitted as 13306599/13306600 at `3-00:00:00`. Reason: raising
  the per-trial cache-build cap 60 s -> 1800 s raised the tuner's worst-case
  cost ~30x across a 6-MAC sweep, but the script's 24 h walltime had not been
  raised to match. R1/R2 left at 24 h. NOTE these are the weakest estimates in
  the campaign — no p2b job has completed at ANY rung, so there is no measured
  reference, only the R1 local smoke.
- `p1_tune.sh` gained a `STAGES` gate (colon/comma separated subset of `1:2:3`,
  **default `1:2:3`** so every other rung's invocation is unchanged). R7
  resubmitted as three chained jobs 13306603 (s1 tune_fmm) -> 13306604
  (s2 fgstune) -> 13306605 (s3 fgsprecond), `afterok`, 72 h and 256G each.
  `afterok` is correct here — unlike the cascade that killed 7 jobs on 08-20,
  this is a genuine data dependency (s2 consumes s1's `tune.csv`). Verified
  `results/phase1/multi/R7/` is empty first: the stage scripts APPEND rather
  than checkpoint-skip, so re-running a completed stage duplicates rows.

**Latent bug found and fixed: `p1_tune_s3.sh` was missing the filament pin.**
It is a pre-08-20 fork of `p1_tune.sh` that never received the BRAINSTORM/025
`FLOWPANEL_FILAMENT_REG=vatistas` export, so any stage-3 recovery run since
08-20 would have silently used the new Gaussian default and produced rows not
comparable with the rest of Phase 1. It was also **cluster-only** — never in the
local repo, and it survived the 08-22 sync only because the `benchmark/` rsync
carried no `--delete`. Now pulled into the repo, pin added, and marked
DEPRECATED in favour of `STAGES=3`; forking the launcher again is what caused
the drift, so future splits go through the gate.

### 2026-08-22 (later) — Fixes synced to both cluster checkouts; all 8 jobs resubmitted

**Sync method (Ryan's ruling):** rsync the WHOLE local `src/` (with `--delete`)
and `benchmark/` to both `~/projects/FLOWPanel.jl` and
`~/projects/p2b/FLOWPanel.jl`, rather than committing or cherry-picking files.
Rationale: the local tree is dirty with mixed 018/022 work and all of
`benchmark/` is untracked, so any partial sync would have left a mixed
Aug-19/Aug-22 generation inside `src/` — the exact failure that already cost
this campaign a run. Whole-tree rsync makes the cluster source byte-for-byte
verifiable against local.

**Provenance (verified post-sync, both checkouts):**

| item | value |
| --- | --- |
| local branch / HEAD | `fastmultipole` @ `d6bf8b6` + uncommitted 021 fixes |
| `src/FLOWPanel_fmm.jl` | `4128f687d3bf1d161e0cdabea640dac9` |
| `src/FLOWPanel_instrumentation.jl` | `a08ac333b20d17da523d6d298024b493` |
| `src/FLOWPanel_solver.jl` | `114428fa7c6aa70528fb2f77a1004dca` |
| `benchmark/rotor_hover_solver_phase1_table.jl` | `3f9d5b8ff6c104e7f16b6091dbbca9df` |
| `benchmark/rotor_hover_solver_phase2b_nearfield_cache.jl` | `58c2e32fd93ca8617bd259ee74f0d8cb` |
| `benchmark/slurm/p2b_nearfield.sh` | `7b9d338ddcbeb352e881c86c93ada974` |
| full `src/*.jl` per-file md5 list | diffs CLEAN vs local on both checkouts |
| `grep -c nearfield_cache_max_bytes src/FLOWPanel_solver.jl` | 7 (was 0) |
| `grep -c 'max_pattern_entries=8192' .../phase1_table.jl` | 1 (was 0) |
| `import FLOWPanel` | `LOAD_OK_MAIN`, `LOAD_OK_P2B` (FLOWPanel recompiled, 240 deps cached) |

The first rsync pass silently transferred nothing and was caught only by the
md5/grep check — a direct vindication of the "verify on the cluster, never
assume" rule. Re-run individually, it landed. Cluster-generated flat-layer CSVs
under `benchmark/results/phase1/multi/` were NOT clobbered (identical size+mtime,
rsync skipped them); the per-rung `R1`–`R6` dirs were untouched and the p2b
R1–R4 knob CSVs + `bcache_R*.bin` prerequisites are intact.

**Failure class 3 (EACCES) is closed by dependency surgery, not by lock cleanup.**
Ryan removed GeoIO and PythonPlot from `Project.toml` on both checkouts on 08-22,
so `import FLOWPanel` no longer runs PythonCall's `__init__` -> `CondaPkg.resolve()`
and there is no `.CondaPkg/lock` to contend for. Deps lists on the two checkouts
are now identical to local's. The new `FLOWPanel_gpu_influence.jl` /
`FLOWPanel_gpu_wake.jl` (task 051) add no dependencies and are inert unless armed.

**Correction to the 08-22 entry above:** 13242660/62/64/66 were NOT
"never submitted". They are in `sacct` as `CANCELLED`, Reason `Dependency`,
00:00:00 elapsed — their `afterok` parent failed and Slurm cancelled the chain.
Same for `p1-table-R6-single-a/-jacobi` (13206078/79) and ALL THREE
`p1-table-R7-*` (13206081/82/83). So the campaign lost 7 jobs to dependency
cascade on top of the 3 that actually failed. The afterok chain was not
malformed; it was doing what afterok does. Resubmission therefore drops the
cross-job chaining wherever the jobs do not share an output directory.

**Resubmitted 2026-08-22, all zen3 + `--exclusive` (128-core / 512 GB nodes):**

| job | id | notes |
| --- | --- | --- |
| p2b-nearfield-R1..R4 | 13306549/50/51/52 | from the p2b checkout, four INDEPENDENT jobs (no afterok), `--mem=256G`, `MODE=multi`, `TUNE_MACS=0.3:...:0.8` |
| p1-tune-R7 | 13306553 | retry of 13206080, `--mem=256G`, `--time=48:00:00` |
| p1-table-R6-multi | 13306554 | retry of 13206077, needs the 8192*N ILU fix; 5 configs (drops `backslash_ldiv` per the R6/R7 drop-out recipe) |
| p1-table-R6-single-a | 13306556 | `krylov_ilu:fgs:fgmres_fgs`, 192G/3-00:00:00 |
| p1-table-R6-single-jacobi | 13306557 | `krylov_jacobi`, 128G/3-00:00:00, `afterany:13306556` |

`--mem` raised to 256G on the p2b rungs (from the script's 64G default) because
the near-field cache cap is now 32 GiB; under `--exclusive` this is free and does
not affect timing comparability. The two R6 `-single` jobs remain serially
chained — they share `results/phase1/single/R6/` and the NFS single-writer rule
applies — but via `afterany`, so an upstream failure no longer silently deletes
downstream work. `p1-table-R7-*` cannot go back in until 13306553 produces R7's
knob CSVs.

Not yet resubmitted / still open: `p1-table-R7-multi`, `p1-table-R7-single-a`,
`p1-table-R7-single-jacobi` (blocked on p1-tune-R7).

### 2026-08-22 — Cluster campaign was 0/8; three independent causes found and two fixed

Live `sacct` contradicted the docs: of the eight Phase-2b jobs treated as in
flight, **four were never submitted** (13242660/62/64/66 absent from Slurm
history — the afterok chain needs checking before relaunch) and four failed.
Phase-1 was also incomplete: `p1-tune-R5/R6` and `p1-table-R5-multi` COMPLETED,
but `p1-tune-R7` (13206080) and `p1-table-R6-multi` (13206077) FAILED. The four
p2b failures were NOT one cause, contrary to the single error text first seen:

1. **Near-field cache caps (R3, R4).** R4 needed 4.47 GiB for 207,344 blocks
   against a 4 GiB cap — a 12% overshoot. R3 tripped the tuner's first-trial
   guard at leaf 58 / MAC 0.3. Root cause of the abort-rather-than-degrade
   behaviour: `tune_fmm` clamps back to `last_feasible_leaf`, but that is
   initialised to `nothing` and `&&` short-circuits, so when the FIRST trial is
   infeasible the `error()` fires and the clamp path below it is dead code
   (`FastMultipole/src/autotune.jl:109-121`). The leaf ladder only ever grows,
   so the message's own advice (start lower) is the only route the code supports.
   Compounding it, the MAC sweep starts at 0.3, its most expensive point — at R1
   the near field is 85% of the full matrix there vs 53% at MAC 0.5.
   FIXED: caps threaded from `KrylovSolver` down to `build_nearfield_cache!`
   (they were previously **unreachable from FLOWPanel** — `FLOWPanel_fmm.jl`
   called the builder with no `max_bytes`, so R4 was hard-blocked); driver
   defaults raised to 32 GiB / 1800 s, per-rung overridable via
   `TUNE_MAX_BYTES` / `TUNE_MAX_BUILD_TIME`.
   NOTE the cache build is **serial** (`nearfield_cache.jl` `_build_nearfield_cache`
   has no `@threads`; the threaded code is `nearfield_matvec!`, the evaluation),
   so the time cap is wall-clock, not a per-thread budget — measured R1 19.5 s at
   4 threads vs 20.0 s at 1. `estimate_nearfield_cache` overshoots actual build
   time by a consistent ~1.65x.
   The old 60 s cap was distorting the result, not just aborting jobs: with caps
   raised the R1 cached-economics tuner picks leaf=275 (0.69 GiB, 57.6 s build)
   rather than the kernel-tuned leaf=21 — the answer sat just past the old cap.

2. **ILU pattern cap (p1-table-R6).** `max_pattern_entries=2048*N`, set at R3,
   undershot R6 by 0.6% (needs 2,059.5 entries/row). Raised to 8192*N across all
   seven benchmark scripts. See the new ILU scaling finding in `ledger.md` +
   `ilu_scaling/` — the pattern grows ~N^1.5, which is a Phase-1 result in its
   own right, not just a cap constant. The guard now reports the total required
   rather than the subtotal at which it tripped.

3. **Cluster filesystem EACCES (p2b R1, R2; p1-tune-R7).** NOT a repo bug.
   `open(".../p2b/FLOWPanel.jl/.CondaPkg/lock") — permission denied` for the p2b
   jobs, and `open("~/.julia/compiled/v1.12/Pkg/*.ji.pidfile") — permission
   denied` for the R7 tune. FLOWPanel hard-depends on PythonPlot
   (`Project.toml:30`), so `import FLOWPanel` always runs PythonCall's `__init__`
   -> `CondaPkg.resolve()` even though the benchmark path never plots;
   `JULIA_PKG_PRECOMPILE_AUTO=0` does not prevent this, as it is module init and
   not precompilation. Ryan is clearing the stale lock state cluster-side
   (2026-08-22 ruling); submit scripts deliberately left untouched. Same depot as
   the in-flight 018 jobs, so worth watching there too.

Verification: 347 unit tests pass (`runtests_unit_solver.jl` +
`runtests_unit_fmm.jl`), including the ILU guard and 021 rigid-motion testsets;
R1 Phase-2b driver smoke passes all three sections and writes both CSVs.

Still open: resubmission (waiting on the EACCES cleanup), the four never-submitted
p2b jobs, and reconciling the item file's `## Current status` / phase-gate table,
which still read "Phase 2 NOT STARTED" and are stale since 08-13.

### 2026-08-21 — Wake-on FGS plateau diagnosed and fixed (H3 lifecycle ordering)

Dense-LU arbitration on the first particle-carrying R1 state identified FGS
as the incorrect path (FGS solution error 2.073e-3 and true relative residual
1.941e-2; fresh Krylov 5.337e-9 and 9.542e-10). The underlying H3 bug was a
one-step control-point lag: `simulate!` transformed persistent solver target
buffers before recalculating moved control points. Reordering geometry refresh
before `transform_body_solvers!` removes the jump. Fixed 0°–80° divergence is
flat at 1.740e-7–1.893e-7; dense wake-on agreement is 1.879e-7. Added the
triaxial panel rigid-motion gate and simulation ordering regression. Details:
`phase_02_single_step_benchmarks.md` 2026-08-21 entry.

### 2026-08-20 — Rigid-motion tree/cache reuse executed; FGS unsteady staleness confirmed + fixed

`rigid_motion_tree_reuse_item.md` executed (Ryan's go via launch prompt).
The pre-registered FGS staleness discriminator CONFIRMED the bug: per-step μ
divergence vs a fresh-Krylov reference grows ~10⁶× with rotation angle
(1.7e-7 at 0° → 21% at 80°, R1, 8 steps) — pre-fix unsteady FGS results are
untrustworthy. FastMultipole gained `transform_tree!`/`transform_plan!`/
`transform_solver!` (commits eea944d/087bf4a/d714544; note Ryan's concurrent
645cc96 swept in the transform_tree! src half) with exact scalar
nearfield-cache persistence under rigid motion and loud v1 refusals for
direction-carrying outputs; FLOWPanel (uncommitted) mirrors per-step rigid
kinematics deltas into persistent solver state inside `simulate!` and adds
`KrylovSolver(persistent_plan=true)` (deferred nearfield-cache commit 4).
All suites green. Details + tables: `phase_02_single_step_benchmarks.md` Log
2026-08-20 (tree-reuse entry) and the item file's "Execution results".

### 2026-08-20 — Phase-2b HPC chains launched

R1–R4 per-rung afterok chains submitted (p2b-nearfield → p2-table-nf:
13242659–13242666, zen3 exclusive, multi 64T) from a SEPARATE cluster
checkout `~/projects/p2b/` (Phase-1 chains still in flight — main checkout
untouched). Full detail + five ops judgment calls (p2b checkout, FLOWVPM
symlink, FM test/-dir ship requirement, JULIA_PKG_PRECOMPILE_AUTO=0, R5–R7
deferred): `phase_02_single_step_benchmarks.md` Log 2026-08-20.

### 2026-08-19 (evening) — nfcache benchmark configs + rigid-motion staging

Per Ryan: three nfcache variants added to the phase2 driver
(gmres/ilu/fgmres_fgs; Jacobi skipped), shared cached-tuned knobs per rung
via new `tune_cached.csv` plumbing; cache build counts inside per-step
(per-solve state, matching the Krylov tree-rebuild convention). R1 smoke:
**all three ≈33–34 s cold, BUILD-DOMINATED (~94% = the 31.7 s per-solve
cache build; preconditioner worth <1 s)** — points at the newly STAGED
`rigid_motion_tree_reuse_item.md` (transform trees under rigid motion,
shared FmmPlan+FGS machinery, cross-timestep cache persistence) as the real
lever. That item also carries a **CORRECTNESS FLAG: FGS unsteady staleness**
(ctor-frozen trees under a rotating body — far field silently degrades with
rotation angle; steady campaign unaffected; verify before any Phase-3
unsteady FGS row). Detail: `phase_02_single_step_benchmarks.md` Log
2026-08-19 "(evening)" entry.

### 2026-08-19 (later) — Cached-path autotuning + build-time guard (Ryan feedback round)

Ryan rulings: FmmPlan cache ownership confirmed (Tree-field rejected);
tune=true refusal superseded. Landed FM 0ef4e83
(`estimate_nearfield_cache` + `max_build_time` guard — estimate from one
kernel sample, never time the build to discover it) + 204188a (tune on the
cached path: provided cache tunes expansion_order ONLY; `tune_fmm
tune_nearfield_cache=true` builds throwaway caches per trial, build cost
excluded, leaf changes stop at the bytes/build-time caps). R1 smoke:
cached-economics knobs **p=14/MAC0.5/leaf=342** vs kernel-tuned
p17/MAC0.5/leaf21 — leaf optimum ~16× larger, cache 862 MB / est build 34 s,
uncapped. All suites green 4T. Detail + rulings:
`phase_02_single_step_benchmarks.md` Log 2026-08-19 "(later)" entry.

### 2026-08-19 — Near-field influence-matrix cache IMPLEMENTED (scheduled autonomous session)

Executed the staged plan: FastMultipole commits 72d3f3d
(`NearfieldInfluenceCache` + builder + deterministic owner-partitioned
`nearfield_matvec!` + standalone `direct!` form + tests) and 02f071c
(`fmm!`/`FmmPlan` integration); FLOWPanel `cache_nearfield` KrylovSolver
surface implemented but UNCOMMITTED per the launch constraint. V0
(`_induced_wake` linearity) PASSED (8.6e-16 on the shedding diamond, 4.3e-16
at R1 scale) — no refusal path needed. R1 local smoke A/B: isolated
near-field **265×**, cold krylov_gmres solve 131→23.7 s (**5.5×**, same 83
iters, certified BC 8.47e-9 both, solution shift 3.5e-15), build 11.9 s /
272 MB / break-even 8.3 applies. All suites green at 4T. Commit 4
(`persistent_plan` cross-solve reuse) deferred pending Ryan. Full detail +
9 flagged judgment calls: `phase_02_single_step_benchmarks.md` Log
2026-08-19 (IMPLEMENTED entry); rows in `ledger.md`; CSV
`benchmark/results/phase2/multi/nearfield_cache_ab.csv`.

### 2026-08-19 — Near-field influence-matrix cache lever STAGED (Ryan)

Plan-only session: new Stage 2b lever (packed dense near-field blocks keyed
to a Tree, general over `direct!`, timestep-reusable; performance-first)
staged as a self-contained implementation plan at
`nearfield_matrix_cache_plan.md`; pointer + design summary in
`phase_02_single_step_benchmarks.md` Log 2026-08-19. No code touched.

### 2026-08-18 (eve) — Phase 2 prep executed concurrent with Phase 1 HPC chains

Ryan authorized concurrent Phase 2 prep at full scope (harness + ALL staged
src levers — the selection constitutes the sign-off Part C2/C3 of the FGS
performance plan were gated on) + a corrected unsteady driver via RHPC
setup-only include. Everything landed and smoke-validated locally; full
detail in `phase_02_single_step_benchmarks.md` Log 2026-08-18. Headlines:
Phase 2a componentized driver + unsteady driver + p2 Slurm launchers (not
submitted) + profiling harness + analysis/TikZ pipeline all working at R1;
`FmmPlan`/`cache_tree` per-solve Krylov tree cache (bitwise-verified,
235/235 solver units); `sweep_order=:colored` FGS sweeps with the batching
theorem proven bitwise (2216/2216) — **but multicolor GS measurably diverges
on a synthetic case where lexicographic converges, so campaign convergence
A/B is a hard adoption gate**; scheduling A/B per Ryan: default
`Threads.@threads` beats `:static` ~4–10%, default kept. Nothing rsynced to
the cluster; Phase 1 evidence untouched (one pure code-motion extraction
from table.jl, smoke-verified). Legacy FGS `reverse_pass` found to be
never-actually-reversed (counter-vs-branch quirk) — preserved bit-identical,
flagged for Ryan. FGS per-step iteration capture in `simulate!` flagged as a
Phase 3 prerequisite. FastMultipole full suite at session end (4T, incl. the
two new test files): exit 0, 73 testset summaries, zero failures/errors.

### 2026-08-17 — Leaf-LU caching implemented and local R1 cache-only A/B complete

Implemented the deferred D3 improvement across FastMultipole and FLOWPanel.
Original self-influence blocks remain untouched for residual matvecs; checked
LU factors live in one distinct contiguous buffer and are reused through
three-argument `ldiv!`. Default is cached; `cache_leaf_lu=false` is the
low-memory/control path. Wrapper keywords and metadata are covered by tests.

FastMultipole's full package suite passes (including 153 cache
invariants/trajectory assertions); FLOWPanel's full regression suite also
passes (focused solver 217/217 and history suites 43/43). Local R1
4T/BLAS-1 min-of-5 at frozen Phase-1 settings preserved iterations and
certified BC errors exactly while improving leaf pass 28.1x, FGS 2.63x, and
FGMRES+FGS 1.49x. The factor cache is 10.16 MiB and built in ~10--11 ms.
Raw evidence is under `benchmark/results/phase2/leaf_lu_cache/multi/`.

No HPC results or cache-enabled production retuning are claimed: Phase 1 has
not frozen R4--R7 settings, so those prerequisites remain pending. No notebook
entry was written.

### 2026-08-17 — FGS nondeterminism attributed to M2L target overlap; fixed and verified

Executed the staged determinism plan on R1 with a true
`{Julia 1T,4T} x {BLAS 1T,4T}` subprocess matrix. The pre-fix 4T/BLAS-1T
replay had exact constructors, tree geometry/topology, lists, packed matrices,
upward/L2L/L2B passes, and representative leaf LU solves, but fixed-input M2L
differed in every repeat. Root cause: FGS fed the raw non-owner-major M2L list
to `assign_m2l!`, whose partitioning requires contiguous target ownership;
separate tasks could `+=` into one target expansion. Canonical stable counting
sorts by source then target now run before matrices/index maps. A premise-guarded
FastMultipole regression requires multiple leaves, nonempty M2L, and exact
repeated cold residual/strength histories.

Post-fix, all four thread coordinates pass five repeated cold solves and three
fresh constructors bit-for-bit. Selected 4T/BLAS-1T fixed-ten minimum changed
0.886844 -> 0.894550 s (+0.87%, gate PASS). FastMultipole counting sort was
3.14x faster than generic lexicographic `sort!` on the 1,760-pair R1 list.
Corrected post-crossing geometric-mean margins pass all 18 R1--R3 checks; R3
tau 1e-6 correctly forced a staircase extension before verification. Evidence
and full attribution table: `results/fgs_determinism/summary.md`. FastMultipole
solve tests 582,335/582,335; FLOWPanel history 15/15; solver units 200/200.
Stopped before LU caching or any performance redesign. No notebook entry.

### 2026-08-13 (later) — First gate review DECLINED on evidence-record consistency; remediated + rerun

The clear-context gate review (subagent, authorized by Ryan) **verified all W6 substance**
(code correctness incl. permutation/equilibration algebra, tests green, tables ≡ CSVs,
acceptance claims recomputed) but **DECLINED** on six findings — chief (Medium): the CSV
of record labeled the cold benchmark solves `warmstart=true`, because `run_config!`
captured `solver_settings` after the `simulate!` leg's `body_solvers_for` mutation
(ruling 10 makes that column authoritative; the timed solves were traced genuinely cold).
Plus: stale control-doc header ("STAGED — not started"), stale phase_00 header/W6 section,
a wrong setup-breakdown claim ("tree+lists ~0.5 s"; actual 0.15 s + ~0.4 s untracked ctor
overhead), and the 08-12 log entry's misordered approve/rescind paragraphs. One
informational note (untested missing-diagonal top-up path — dead under Barba
`self_induced=true`, guarded loudly; no action).

Remediation, same day: harness captures `solver_settings` before the sim leg (flip
recorded in `notes`); all four doc findings fixed; **both smoke modes rerun** (agent
judgment: reviewer deemed rerun optional, but clean CSVs beat an annotation for a
publishable campaign). Per Ryan (back online): the reruns ran **concurrently** (5
threads total) with **multi at 4 threads** — accepted as CSVs of record with caveats
noted in the phase_00 tables; iteration counts/residuals reproduced exactly, wall times
carry contention (e.g. plain GMRES single 281 s vs 223 s quiet). Ryan rulings: proceed
into Phase 1 on the re-review's approval; no notebook entry yet.

**Second gate review: APPROVED** (same day). All six findings verified fixed/no-action;
both 08-13 tables cross-checked cell-for-cell against the new CSVs; acceptance holds in
both modes; unit solver suite green. Two informational notes, no action: the sim-leg
note in `notes` is appended unconditionally (cosmetic; `solver_settings` authoritative
and correct), and "lists 0.04 s" matches multi (single shows 0.058 s, within "~").
Phase 0 gate: clear-context approved; Ryan's re-approval checkbox pending his review.
**Phase 1 begins** per Ryan's 2026-08-13 ruling.

### 2026-08-13 — W6 complete: ILU implemented (Ryan), reviewed, smoke-run PASS both modes

Ryan implemented W6 directly on the evening of 2026-08-12 — `phase_00_ilu_plan.md`'s
Stages A–C collapsed into one sitting; the plan file is marked superseded. Design
decisions made by the implementation: pattern = FastMultipole **Barba direct interaction
list** (dedicated ctor-time tree, not the staged radius-cutoff fallback); factorization =
**`ILUZero.jl`** (new dep). This session (agent): reviewed the implementation for
correctness and performance — verdict CORRECT (assembly ≡ `_G!` semantics; permutation +
equilibration algebra checked; ILUZero internals audited at source). Fixes applied:
ctor deduplicated onto `_ilu_direct_pattern` (helper now returns timings); assembly
threaded over direct-list pairs (mirrors `_G!`; `induced` is pure). Two crash bugs in the
harness's new code fixed at launch time: `\"` inside `$()` interpolation (1.12 parser
rejects) and empty-generator `sum` in the `setup_total` fallback. Unit + full suites
green.

Smoke runs relaunched from scratch (both modes, quiet machine; the harness truncates
`runs.csv` per launch, so the 2026-08-12 corrected-rerun CSVs are replaced on disk —
carried-over configs reproduce to ~1–2%, and the old numbers persist in
`phase_00_availability.md`'s tables). **W6 acceptance PASS**: ILU-GMRES 290→15 iters,
223.4→11.51 s (single) and 55.2→2.86 s (multi) vs plain GMRES; setup ~1.0–1.2 s,
break-even <0.02 solves; displaces FGMRES+FGS as best iterative config at smoke scale.
Full tables + stats in `phase_00_availability.md` Log 2026-08-13. Side observation: multi
FGS hit 500 iterations this run (prior rerun: 267) — thread-order nondeterminism,
consistent with the pre-registered Phase 2 finding.

Phase 0 status: technically complete incl. W6; pending clear-context review + Ryan's
gate re-approval. Ryan authorized (2026-08-12) autonomous continuation incl. a subagent
clear-context gate review, then Phase 1.

### 2026-08-12 — Rulings 11/12; clear-context review findings; remediation

Morning: Ryan's GMRES directives adopted as ruling 11 (right-preconditioning honesty —
landed in code, Jacobi now routed `N=`; Phase-1-calibrated tolerances; FMM-floor
measurement; preconditioner exploration with sparse near-field ILU as priority — theory doc
`theory_nearfield_ilu.md` written) and ruling 12 (warmstart benefit = headline deliverable;
Phase 2 strictly cold, Phase 3 owns the warmstart matrix incl. break-even step counts).

Afternoon: a clear-context review declined Phase 0 approval with 5 findings. Verified
assessment — all confirmed, one scope correction:

1. (High, CONFIRMED) **Warm-state FGS/block-GS evidence**: the smoke timing loop never
   reset `body.strength`; FGS seeds from current strengths and breaks on the residual check
   *before* sweeping, so every timed standalone-FGS rep after warmup was a no-op residual
   check. The "FGS converges in ONE outer iteration / 0.047 s" headline was an artifact and
   had propagated to the phase log, control doc, INDEX, and agent memory. Block-GS's
   1-outer history was the same warm-entry artifact. **Scope correction**: config 1f
   (FGMRES+FGS) was NOT warm-biased — Krylov never seeds from `body.strength` and
   `FGSPreconditioner.ldiv!` zeroes strengths every apply (the linearity contract) — so
   "1f ≫ plain GMRES, 26 vs 290 iterations" survives; 1c/1d timings were also genuinely
   cold per rep.
2. (High, CONFIRMED) Config 1d CSVs predate the ruling-11a right-routing change; the
   left-routed "Jacobi hurts, 362 iters" conclusion must be re-derived.
3. (Medium, CONFIRMED) CSV provenance gaps (`julia_version`, tolerances, FMM knobs,
   `t_rhs`, banner stdout-only).
4. (Medium, CONFIRMED) `Base.summarysize(solver)` had inconsistent object boundaries
   (Krylov holds the body; Backslash doesn't).
5. (Medium, CONFIRMED) FastMultipole callback was an uncommitted dev-checkout edit.

Remediation landed same day: cold-state reset before every timed/history/alloc solve
(`min_of_k` gained an untimed `setup!`; FGS cold tripwire assert); provenance columns
`fm_commit, julia_version, solver_settings, backend_settings` + populated `t_rhs` +
`banner.txt` (schema updated in `decision_rules.md`); memory metric finalized as
`solver_state_bytes` (summarysize minus referenced bodies — comparable boundary);
FastMultipole callback committed as `5adde3b` (dev checkout, local commit) + portable diff
at `patches/fastmultipole_callback.patch`.

**Rerun complete same day** (single mode clean; the first multi rerun was discarded —
Ryan had a concurrent job on the machine — and repeated on a quiet machine). Corrected
tables + headline findings in `phase_00_availability.md` Log 2026-08-12; all four doc
surfaces replaced. Biggest corrections: cold FGS = 35.2 s / 204 iterations (not 1/0.047 s);
"Jacobi hurts" re-confirmed honestly; NEW measured finding — FGS multi-thread solve is
2.2× slower with more iterations (267 vs 204; threaded sweep-order convergence
degradation), now isolated to the solve itself. Untouched configs reproduced to ~1%.
Phase 0 back to technically complete. **Ryan approved the Phase 0 gate later the same day**
(recorded in the control-doc decision log and gate table).

**Evening: that Phase 0 approval was RESCINDED by Ryan** — ILU pulled forward from
Phase 2b into Phase 0 as **W6** (develop now: theory doc → implementation design →
implementation/testing; roster gains config (g) GMRES+ILU). Self-contained execution
plan for a fresh agent staged at `phase_00_ilu_plan.md`; W1–W5 evidence and remediation
stand unchanged. Next action: execute W6 — resolved in the 2026-08-13 entry above.

### 2026-08-11 — Phase 0 started; W1 dropped; D2/D3 ruled

Ryan gave the Phase 0 go-ahead. Planning exploration (2026-08-07) verified the solver
inventory and produced three corrections now folded into the control doc:

1. The rotor pipeline is **not** blocked by the Kutta gate — `rotor_hover.jl` uses the legacy
   linear Kutta condition; the `FLOWPanel_kutta.jl:495-503` hard error only fires for the
   opt-in 015 pressure-continuity closure (Dirichlet paired bodies).
2. `BackslashCoupled` availability bug: ctor installs a dummy identity `lu!`
   (`FLOWPanel_solver.jl:989`) awaiting `update_G=true`, but `solve_formulation!` never passes
   `update_G` (default false) — inside `simulate!`, roster config 1a silently solves against
   the identity. Fixed in Phase 0 (WP-A).
3. `Backslash` stale-`Glu` corruption: ctor's `lu!(G)` aliases `G` as `Glu.factors`; a later
   `_solve!(update_G=true)` factorizes into a *local* LU, leaving `solver.Glu` with new
   factors + old pivots — garbage for direct consumers (Kutta Route A, GreenReconstruction,
   DirectWakePotential). Fixed in Phase 0 (WP-D).

**Ryan rulings (recorded in control-doc decision log):**

- **W1 dropped.** The campaign runs the legacy linear Kutta condition throughout; upgrade to
  the 015 pressure-continuity closure only if accuracy/convergence later demand it. The
  kutta.jl hard error stays. Phase 0 W1 is rescoped to the two bug fixes above.
- **D2 — FGS convergence history via FastMultipole callback**: add a non-breaking
  `callback=nothing` kwarg to `FastGaussSeidel` `solve!` in the dev checkout
  (`../FastMultipole/src/solve.jl`), invoked once per outer iteration with
  (iteration, residual) so FLOWPanel records exact production-loop timestamps (ruling 4).
  Cross-repo edit — noted here per the repo boundary.
- **D3 — leaf-LU caching deferred to Phase 2b.** Discovery: FastGaussSeidel leaf blocks are
  dense and re-factorized every sweep (`mat \ rhs`, FastMultipole `solve.jl:934`) — the
  control doc's "leaf LUs built once" was wrong (corrected). This inflates FGSSolver and every
  FGS-preconditioner apply; Phase 1's 1e/1f calibration numbers must be read with that caveat.

**Archived W1/Kutta design (revivable as its own item if the closure is ever needed):**
option (a) solver-generic inner solve — `Backslash` mutable + `Glu` write-back;
`KuttaRuntime{TS<:Backslash}` widened to `AbstractSolver`; `_kutta_inner_solve!(rt, body)`
replacing the direct `_kutta_ldiv!` (Backslash = ldiv path; Krylov = production `_solve!` with
forced warmstart from `rt.current.mu`, essential at ~10–50 trials/step; FGS = seed
`strength[:,2]`, finite-only status); validator relaxed to accept
Backslash|KrylovSolver|FGSSolver, still rejecting coupled/tuple≠1; Route B for matrix-free =
`_assemble_W!` + version bump (no G rebuild); FGS × TEAnchoredAttachment needs per-step
FastGaussSeidel rebuild (stale leaf matrices) — the contested piece. Correctness test design:
unsteady wing, Krylov/FGS vs Backslash at Phase 1 provisional thresholds.

Implementation plan for this session: `/Users/ryan/.claude/plans/work-on-brainstorm-item-bright-lobster.md`
(work packages A→B→C→E→F→D). Also noted: single-body tuple block-GS always runs a second
identical solve to observe Δ=0 — config 1b timing must record solve counts (harness notes).

**Same-day completion:** all work packages landed and Phase 0 exit criteria met (results and
per-config smoke tables in `phase_00_availability.md` Log, 2026-08-11 entry). Beyond the
planned scope, the smoke run exposed and fixed a third availability bug: FastMultipole's
`strength_to_value` had no overload for `RigidWakeBody{<:Any,1}`, so FGSSolver had never
actually worked on the uncapped Neumann rotor body type (the commented FGSSolver block in
`rotor_hover.jl` was aspirational). Fixed in `src/FLOWPanel_liftingbody.jl`, test-gated.
Cross-repo edits to `../FastMultipole` this session: the D2 `callback` kwarg only
(`src/solve.jl`); FastMultipole solve tests pass. Removed stdout `@show workspace.stats` from
Krylov solves (aligned with judge-from-CSVs). Full FLOWPanel suite green (45 testsets).
Phase 0 gate approval + notebook entry pend Ryan.

### 2026-08-07 — Item opened and staged

Item created from Ryan's 2026-08-06 brief plus a solver-subsystem inventory (control doc
"Solver matrix"). Key discoveries baked into Phase 0: Kutta closure hard-errors on
Krylov/FGS; FGS-as-preconditioner and Krylov `x0` don't exist; no convergence-history capture
anywhere. Ryan's decisions recorded in the control-doc decision log (benchmark scope, HPC,
dual threading modes, NACA 0015 wing). No code touched. Next: Phase 0 on go-ahead.
