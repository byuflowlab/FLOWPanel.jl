# 021 Session Log

Newest first. Narrative only — results go to the phase files and `ledger.md`.

## Dated entries

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
