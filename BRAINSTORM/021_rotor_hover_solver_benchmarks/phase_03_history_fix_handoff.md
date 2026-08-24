# Handoff — BRAINSTORM 021 Phase 3: fix the solution-history/outer-loop interaction

Written 2026-08-22 night at a context reset. **Ryan's instruction for the next
session: enter plan mode and design the fix.** Do not implement before the plan
is approved.

---

## What this session did (do NOT re-derive)

Phase 3 was planned and its harness built. Plan file:
`~/.claude/plans/gentle-mixing-leaf.md` (approved). Then a local smoke found a
real defect, which is now the whole subject of the next session.

**Ryan's Phase-3 scope decisions, already made:** multi-thread only; R1 + R2 to
start; ALL solvers and ALL warmstart types; additional rungs later. Wake-developed
state comes from ONE shared checkpoint per rung that every arm restarts from.
Checkpoint retention needs **no protect-list entry** — the VTK sweeper's standing
newest-36-timestep window always leaves the tail intact and `RESTART_STEP=-1`
resolves to the last PVD entry.

**HPC disk sweep: DONE.** 118.5 G → 87.7 G of the 200 G cap, 30.7 G freed, all 8
live fp-018 jobs confirmed still RUNNING afterward, 10 pending 021 jobs
untouched. A ledger line for `018/ledger.md` was produced and is NOT yet pasted.

### Code landed this session (all uncommitted)

| File | Change |
| --- | --- |
| `src/FLOWPanel_solver.jl` | `_extrapolation_coefficient(order, j)` + materialized `_extrapolation_coefficients(order)`, shared by FGS and Krylov; `project_solution_order = 0` now means plain previous-solution reuse (needs only 1 sample; order ≥ 1 keeps the old 2-sample guard so existing runs are bit-unchanged); `FGSSolver.niter` / `.solved` captured via a chaining callback; `KrylovSolver.warmstart_order` + `x_history` + `x0_scratch` (allocated ONLY at order ≥ 1, to keep `solver_state_bytes` unchanged by default — ruling 8); `solve!(bodies::Tuple, solvers::Tuple)` gained `solver_optargs` |
| `src/FLOWPanel_metadata.jl` | records `warmstart_order` |
| `test/runtests_unit_solver.jl` | 5 new testsets, 48 assertions, all passing |
| `benchmark/rotor_hover_solver_unsteady.jl` | `WARMSTART`, `WARMSTART_ORDER`, `RESTART_*`, `SKIP_STEPS`, `PHASE`, `SOLVER_VERBOSE`; restart via `simulate_warmstart!`; new CSV columns `solved, warmstart, warmstart_order, restart_step, skip_steps, strength_checksum`; per-step `solved` guard; warmstart-suffixed `RUN_NAME`; a guard refusing to write into the checkpoint it restarts from |
| `benchmark/slurm/p3_warmstart.sh` | NEW launcher, `WARMSTARTS` gate (cold:prev:extrap), 16 arms/rung → 6 jobs |
| `benchmark/phase3_analysis.jl` | NEW; reduction table, break-even, order sweep, TikZ figure + data dir |
| `.../decision_rules.md` | amended — **CONTAINS A KNOWN ERROR, see below** |

**Test status:** full suite passes EXCEPT 5 pre-existing failures in
`test/runtests_unit_kernel_gradient.jl`, confirmed by `git stash` to fail
identically without this session's changes. Root cause: the uncommitted
`_onplane_snap` in `src/FLOWPanel_elements_fmm.jl` snaps `target_Rz` to exact
zero for plain floats but passes AD duals through, so at on-plane targets
`velocity` takes the PV branch while `FD.gradient` takes the smooth atan branch
and the two no longer agree. **Not Phase 3's bug**, but it lives in the filament
kernel that Phase 3 is the first 021 phase to actually exercise.

---

## The defect to fix

**Symptom.** FGS with `WARMSTART=extrap` (order ≥ 1): solves per time step go
from 2 to **100**, wall time 4.4 s → ~90 s per step (20×), while the recorded
`niter` reads 0–2 and hides it completely.

**Measured, R1/FGS/local 4T, 6 steps** (wake-developed window = steps 4–6):

> **VOID (2026-08-23) — do not cite these numbers.** They were measured on the
> defective code and through defective instrumentation, so both columns are
> unusable:
>
> * the `niter/step` column sampled `solver.niter` *after* the step, which
>   returns the **last** inner solve of the step. In a step containing many
>   repeated solves that solve is already warm, so the "6 → 0" reduction is a
>   last-solve sampling artifact, not a warm-start benefit. `niter_first` (added
>   in the fix) is the count that measures the step-to-step guess;
> * the wall times (cold 4.378 s / prev 3.872 s / extrap 81.13 s) were inflated
>   by the redundant and, for extrap, runaway repeated solves within each step.
>
> The **standing one-body repeated-solve tax** these numbers exposed is resolved:
> a 1-tuple now breaks out of the block-GS loop after the first block solve, and
> warm-start history advances once per completed step. Post-fix numbers replace
> this table in `log.md` (2026-08-23 entry).

| arm | niter/step | t_solve mean | vs cold | particles | checksum |
| --- | --- | --- | --- | --- | --- |
| cold | 6 (VOID) | 4.378 s (VOID) | — | 328 | 113678.5589 |
| prev (order 0) | 0 (VOID) | 3.872 s (VOID) | −11.6% (VOID) | 328 | 113678.5642 |
| extrap (order 2) | 1 (VOID) | 81.13 s (VOID) | **+1753%** (VOID) | 328 | 113678.5764 |

All three reach the same answer (checksums agree to ~1e-7 relative, identical
particle counts), so this is a cost and measurement-validity defect, not a
wrong-answer defect.

**Root cause, verified.** `save_solution!` fires on every `_solve!` call, but
`solve!(bodies::Tuple, solvers::Tuple)` (`src/FLOWPanel_solver.jl:1625`) calls
`_solve!` repeatedly *within a single time step*. So history slots are intra-step
repeats, not consecutive time steps. Consequences, in order:

1. Order-1 extrapolation computes `2·H₁ − H₂` across a duplicated pair, so the
   guess **alternates** around the answer. Straight from the verbose log:
   ```
   projected  0.0021339633652076324 → solved to 0.002917102160351033  (Δ +0.00078)
   projected  0.0037002409554944340 → solved to 0.002917102970431909  (Δ −0.00078)
   ```
2. Once the guess is good, FGS's tolerance check — which runs at the TOP of the
   loop, before any sweep — accepts the projected vector verbatim: 0 sweeps,
   `actual - projected == [0.0, 0.0, 0.0, 0.0, 0.0]`. **The solve becomes the
   identity map.**
3. The outer loop's fixed-point test therefore compares two *different*
   projections rather than two solves. `max_delta` never falls under
   `outer_tolerance = 1e-8` (`:1677`), so it runs to `max_outer_iterations = 50`
   (`:1676`) every step. Per-step solve counts measured directly from the log:
   2 (cold) vs 100 (extrap).

**This is pre-existing.** `project_solution` / `project_solution_order` predate
this work; this session added order-0 support and mirrored the scheme onto
Krylov. Any current use of FGS `project_solution=true, project_solution_order ≥ 1`
through a tuple solve hits it. **The `KrylovSolver.warmstart_order ≥ 1` path added
this session inherits the identical flaw** — it also pushes history inside
`_krylov_launch!`, once per solve. Fix both.

**Two hypotheses that were tested and are WRONG — do not revisit:** it is not
dynamic expansion order (FGS runs at fixed `p=6`; FastMultipole's solve paths
pass `error_tolerance=nothing`), and it is not the inner solver running to
`max_iterations` (it converges in 0–2 sweeps throughout).

---

## Corrections owed

> **ALL RESOLVED 2026-08-23.** The fix landed; see `log.md`, 2026-08-23 entry,
> and the amended `decision_rules.md`. Kept below as the original statement of
> what was owed.

- **`decision_rules.md` is wrong** where it says the block-GS outer loop is
  "degenerate" for the 1-tuple rotor case and defers outer counts to Phase 4. It
  is not degenerate: it runs ≥ 2 iterations always, and up to 50 under
  extrapolation. Rewrite that subsection.
- **The `prev` 6→0 result overstates the benefit.** The driver reads `niter`
  after `solve_formulation!` returns, i.e. from the LAST outer solve, which is
  seeded with the same step's own answer. The genuine step-to-step measurement is
  the FIRST outer solve, which nothing currently records. Phase 3's headline
  metric is not trustworthy until this is fixed too.
- **Standing ~2× cost, independent of Phase 3:** for a 1-tuple, `sources` is
  empty and nothing changes between outer iterations, so the second solve is
  provably identical to the first (verified: bit-identical MSE ladders). Every
  single-body unsteady step in the whole campaign pays 2× for it. Worth its own
  item — it is not Phase 3's to fix, but Phase 3's numbers sit on top of it.

---

## The three fix options Ryan raised, and the assessment

Ranked on Ryan's stated priorities: **robustness > performance > minimal
invasiveness.**

### Option 2 — move the history update out of `_solve!` — RECOMMENDED

Keep `project_solution!` where it is (it must stay inside `_solve!`, since
`set_strengths!` at `:245` zeroes μ before `_solve!` runs — that constraint is
why it lives there), but move **`save_solution!` only** to a hook that fires once
the step's solve is complete.

- **Robustness — decisively best.** With the history frozen for the duration of a
  step, every outer iteration projects the *identical* guess, so the solve returns
  the identical result and `max_delta` is exactly 0. The outer loop is a fixed
  point again **by construction**, not by a flag being maintained correctly. It
  also makes the history mean what the extrapolation already assumes: one entry
  per timestep, holding converged strengths.
- **Performance — best.** One projection-worth of arithmetic and one save per
  step instead of one per outer iteration, and the runaway disappears.
- **Invasiveness — moderate.** Needs a small generic hook (e.g.
  `finalize_step_solution!(body, solver)`, no-op by default) called after the
  outer loop in `solve!(bodies::Tuple, solvers::Tuple)` and on the direct
  single-body path (`src/FLOWPanel_formulation.jl:803`, and the
  `GreenReconstruction` route at `:899` which calls `_solve!` directly). Then
  delete two lines from `_solve!`. The same hook serves the Krylov history.

### Option 3 — slot selected by timestep index

- **Robustness — middling.** Fixes the pollution, but the guess can still differ
  between outer iterations (slot 1 gets overwritten mid-step), so the outer loop's
  fixed point is not guaranteed the way Option 2 guarantees it. Also depends on an
  externally maintained index; a caller that forgets to bump it fails silently.
- **Performance — fine**, same as Option 2.
- **Invasiveness — low-ish.** `simulate!` already iterates `body_solvers` each step
  at `src/FLOWPanel_simulate.jl:1359` (`transform_body_solvers!`), a natural home
  for a counter bump.

### Option 1 — `solved` flag on `AbstractBody`

- **Robustness — worst**, despite being the most explicit. It adds mutable state
  that must be reset and set correctly on every path that solves — `simulate!`,
  `steady!`, the static Phase-1 drivers, tests, examples. A stale `true` silently
  disables history; a stale `false` double-saves. That is a new silent-failure
  class, and this campaign has repeatedly been bitten by exactly that (the
  shedding-from-raw-cells gotcha, the FGS staleness bug).
- **Performance — fine.**
- **Invasiveness — worst.** A new field on `AbstractBody{E,N,TF,DBC}` touches
  every body type and constructor, plus metadata/replay/warmstart-restore paths,
  and Ryan's own note anticipates needing a higher-level solve dispatch on top.

**Recommendation: Option 2.** It is the only one that restores the invariant the
outer loop's convergence test actually depends on — *within a timestep, the solve
is a deterministic function of body state* — structurally rather than
conditionally. Option 3 is the acceptable fallback if the hook plumbing proves
uglier than expected. Option 1 should be declined.

Design the fix to cover **both** FGS `solution_history` and the Krylov
`x_history`, and to also fix the "record the first outer solve, not the last"
measurement problem, since both are the same handoff.

---

## Housekeeping

- Whole tree uncommitted except `benchmark/` which is entirely **untracked**;
  cluster checkouts are rsync copies with no git anchor. Ryan elected to skip
  committing earlier this session.
- Smoke/diagnostic artifacts to delete when done:
  `benchmark/results/phase3_smoke/`, `benchmark/results/phase3_diag/`.
- Reproduce the defect in ~10 min:
  `RUNG=R1 CONFIG=fgs N_STEPS=6 WARMSTART=extrap WARMSTART_ORDER=1 SOLVER_VERBOSE=1 PHASE=phase3_diag EXPECT_JULIA_THREADS=4 THREADING_MODE=multi julia --project -t 4 benchmark/rotor_hover_solver_unsteady.jl`
  Count solves per step with:
  `awk '/projected first|projection skipped/{n++} /step [0-9]+\/6/{print n" solves"; n=0}'`
- Local runs: **never more than 4 threads.**
- OPEN TODO still carried at the top of `log.md`: relabel the family attribution
  on the 10 in-flight jobs when they finish (~2026-08-25).
- Phase 3 remains `NOT STARTED` with its approval box unticked. Nothing has been
  submitted to HPC. Checkpoint jobs, arm submission and results all still need
  Ryan's explicit go.
