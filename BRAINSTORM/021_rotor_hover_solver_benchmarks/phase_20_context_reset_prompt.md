# BRAINSTORM 021 — Phase 20: context-reset handoff

**Written 2026-08-25 (late) while 7 jobs run. Read this file, then
`phase_18_context_reset_prompt.md` for cluster state and the precompilation
root cause (§2, §3 there are still current and still authoritative). If you need
the campaign's rules, read `decision_rules.md`. You should NOT need `ledger.md`,
`log.md`, or the item file to start. Do not read `phase_01_consistency.md`
(73 KB) or `phase_02_single_step_benchmarks.md` (41 KB) end-to-end — grep them.**

Predecessor: `phase_19_fgs_niter_followups.md`. Its §5 banner is now a record of
two *superseded* designs — read it only for that history, not for instructions.

---

## 0. START HERE — check HPC job progress first

**Ryan's instruction at the reset: the next agent checks on HPC job progress.**
Do that before anything else in this file, then report to him.

Delegate it — do NOT ssh or tail logs inline (per `CLAUDE.md`, that output is
exactly what fills a fresh context):

- **`hpc-monitor`** subagent — job status, log tails, disk usage. Strictly
  read-only; it never submits, cancels, or deletes.
- **`hpc-storage`** subagent — only if disk needs sweeping (see the cap below).

**Connection gotcha:** `ssh orc` needs a **live ControlMaster socket** or it
hits 2FA and hangs. If the subagent reports a stalled/failed connection, that is
the likely cause — ask Ryan to open a session rather than retrying blindly.

**Ask the monitor for:** per-job state (R/PD/CD/F/TO), elapsed vs walltime,
whether each job is *producing output* (not just running — the 08-24 launch sat
12 h in package loading while looking healthy, see `phase_18` §3), newest CSV row
timestamps, and total usage under `/home/rander39`.

### The seven jobs in flight (as launched; verify, do not trust this table)

All `physics2`/`standby`, `--constraint=zen3 --exclusive --mem=500G --requeue`,
64 threads, from `~/flowpanel-021/FLOWPanel.jl` on `orc`.

| job id | name | walltime | started | lands ~ | what it unblocks |
| --- | --- | --- | --- | --- | --- |
| 13477933 | p2-tune-R1 | 24 h | 08-25 22:32 | **08-26 22:30** | Phase 2 Step A; first rung of the Phase 3 gate |
| 13477934 | p2-tune-R2 | 2 d | 08-25 22:32 | 08-27 | Step A |
| 13477935 | p2-tune-R3 | 4 d | 08-25 22:32 | 08-29 | Step A |
| 13477936 | p2-tune-R4 | 6 d | 08-25 22:32 | 08-31 | Step A |
| 13477937 | p2-tune-R5 | 7 d | 08-25 22:32 | 09-01 | Step A; **needs a 2nd cycle** |
| 13469157 | p1-agree-R6 | 2 d | 08-25 10:25 | **08-27 10:25** | Phase 1 R6 agreement |
| 13469158 | p1-agree-R7 | 3 d | 08-25 10:25 | 08-28 | Phase 1 R7 agreement |

**Nothing was queued or pending at reset.**

### What to do with what you find

| finding | action |
| --- | --- |
| all running and producing | report the table; **do nothing else on the cluster** |
| a job produced nothing for hours | suspect the precompilation race — `phase_18` §3 has the verbatim symptom and the fix. Report; **do not cancel** |
| a tuning job hit its wall timer again | that was the 08-25 finding (descents timing out, not converging). Report the evidence; the 5× timer raise was Ryan's call and any further change is too |
| a job finished | harvest per `phase_17_harvest_prompt.md` — its §4 trap list and §7 quarantine list are **still authoritative** before publishing any table |
| disk over 300 G | escalate retention to 72 timesteps; **400 G cap across ALL of `/home/rander39`**, standing retention 288, sweep trigger 150 G. Use `hpc-storage`; only the sweeper script mirrors to the cluster |

**Every `sbatch` / `scancel` / `scontrol` is Ryan's call.** Reporting a dead job
is the deliverable; fixing it without being asked is not.

---

## 1. What this session did, in one paragraph

Closed the FGS per-step iteration-capture thread (it was already implemented; a
stale doc said otherwise, now corrected in place), fixed the last two drivers
hard-coding `niter = -1`, and — the substantive change — **rebuilt the Phase 3
validity guard twice**. The 2026-08-23 report that the extrapolated warm arm had
DIVERGED was traced to the test, not the solver. A relative-L2 replacement was
built, then discarded on Ryan's ruling in favour of measuring the **boundary
condition** directly. That design is implemented, verified on a live R1 run, and
documented. **No cluster interaction of any kind occurred.**

---

## 2. The ruling that governs everything here

> **Ryan, 2026-08-25:** "There is no need to measure the divergence of the
> solution itself between warm and cold solves; rather, the residual should be
> measured... it's easier to show that we're satisfying the boundary
> conditions."

**Why the old guard was invalid** (this reasoning stands and is worth not
re-deriving): the solver's acceptance test is `max|b − Ax| ≤ tol_abs`, absolute
and initial-guess-independent (`FastMultipole/src/solve.jl:1222`; `residual!` at
`:1530-1566` returns max-abs — the `mse` name is vestigial, see its `:1210`
comment). Cold and warm call the *identical* `FastMultipole.solve!`
(`src/FLOWPanel_solver.jl:1694`); the only difference is `project_solution!` at
`:1659`. They land at different points of the **same residual ball**, of
solution-space diameter `‖A⁻¹‖·2·tol_abs`. The old bound assumed `‖A⁻¹‖ = 1` —
it read a *residual* tolerance as a *solution* tolerance. Phase 1 measures the
amplification at ≈12 (R1), so the observed 2.4e-7 sat ~3 decades *inside* the
solver's own ambiguity. **extrap was never diverging.**

Full rules: `decision_rules.md` → **BC satisfaction guard**.

---

## 3. What is implemented and verified — do not rebuild it

| what | where | state |
| --- | --- | --- |
| `niter = -1` sentinels | `phase1_agreement.jl:439`, `phase1_solvetime.jl:175` | FIXED, mirror `unsteady.jl:241`. Future re-runs only — R1–R7 were spooled under the old driver |
| per-step BC check | `rotor_hover_solver_unsteady.jl`, `TimedFormulation` | ONE `bc_error!` pass per step; emits `bcerr_{max,min,q1,med,q3,rms,tol,eps,certified}`, `t_bcerr` |
| per-panel φ exposure | `benchmark/common.jl` `bc_error!` | new `phi_out` kwarg, defaulted — all 10 existing callers verified unaffected |
| identity reporting | same driver | per-step `n_particles` + `CT` (needs `RUN_MONITORS=true`, defaulted on under `PHASE=phase3*`, Bernoulli-only) |
| timing correction | same driver | `t_step_net = t_step_total − t_bcerr`; end-of-run share quoted against net wall (see §5) |
| analysis | `phase3_analysis.jl` | guard + "BC satisfaction and identity" table; old `1e-8` constant deleted |
| docs | `decision_rules.md`, `phase_03_launch_prompt.md:103`, `phase_19` §5, item Phase-3 gate row + decision log | all updated; superseded designs marked in place, never deleted |

**Live verification, R1/fgs, cold + prev, 5 steps each (2026-08-25):**

| arm | max \|φ\|/tol | band | med \|φ\|/tol | tail max/med | certified | CT final | ΔCT % |
| --- | --- | --- | --- | --- | --- | --- | --- |
| cold | 1.007 | 1.1 | 0.889 | 15.65 | 5/5 | 0.1058 | 0.0 |
| prev | 0.9939 | 1.1 | 0.8642 | 19.03 | 5/5 | 0.1058 | **1.3e-5** |

Both arms sit right at the bar they promised — that is what a correct solve
looks like. Solver unit suite 397/397 at the time of the sentinel fix.

### Two design points that will look wrong if you skim them

- **The pass band is `1 + ε/tol`, not `1`.** The arm stops on its own leaf-local
  estimate at its tuned expansion order; `bc_error!` re-evaluates the same
  quantity independently and more sharply, so they legitimately differ by about
  the error *requested* of the measurement. That band is the measurement's own
  **certified** uncertainty, not a fudge factor. Ratios inside it report as
  `marginal`, so drift toward the bar stays visible. Do not widen it, and do not
  narrow it to 1.0 — a bare 1.0 flagged the *cold* arm at 1.007.
- **Accuracy is requested, not guessed.** Ryan asked for ~10× sharper than the
  solve. `bc_error!` requests an absolute FMM error (`BCERR_SAFETY`, default 0.1,
  × the arm's own tolerance) and reports whether every M2L pair certified it —
  strictly better than bumping the expansion order and hoping. Verified
  `bcerr_eps/bcerr_tol = 0.100` exactly, certified on every step.

---

## 4. Traps

- **The rotor is DIRICHLET, not Neumann.** `kernel = Union{ConstantSource,
  VortexRing}` ⇒ `DBC = true` (`rotor_hover_pressure_comparison.jl:283-287`), and
  `solution_column = 2`. This killed an earlier design that captured the RHS
  before the solve: for a Dirichlet body the RHS is only reproducible *inside*
  the solve. `bc_error!` sidesteps it — it reloads σ and the doublets together in
  one pass. The driver asserts Dirichlet up front and errors otherwise.
- **`bc_error!`'s entry contract:** `body.velocity` must hold the apparent
  velocity at the control points. The driver snapshots it *before* the inner
  solve and restores it for the measurement — the unsteady analogue of Phase 1's
  `rotor.velocity .= frozen_velocity`. Do not remove that snapshot.
- **`tol_abs = 1.77e-8` is CORRECT** for the Phase 3 winner knob set
  (`p=6;mac=0.3;leaf=50;inner=10`). A previous session "corrected" it to 9.93e-8
  from `fgstune_margin_verify_R1.csv` — that row is a *different* knob set
  (leaf=150, inner=5). The live run emits `1.7658196265827666e-8`. Retraction is
  recorded in `phase_19` §5.
- **`niter` is not a homogeneous column** — Krylov iterations, FGS sweeps, and
  outer FGMRES iterations all share it. `phase1_table.jl:254` is the outer FGMRES
  count *by design*; do not "fix" it.
- **Do not re-read `phase_19` as a task list.** Its items are done. Its §5 is
  history.

---

## 5. Cost, and the one thing that has NOT been re-measured

The BC pass ran **~62% of `t_solve`** at R1 (`t_bcerr` column). **Ryan ruled the
wall-time cost acceptable** ("run it twice, or use the replay function"), so the
measurement runs **inline in the same run** — no second pass, no replay — and
`BCERR_EVERY` defaults to 1 (every step) as an escape hatch only. Skipped steps
write EMPTY cells, never a fabricated verdict.

**What that inline choice contaminates, and how it is handled:**

| quantity | clean? | why |
| --- | --- | --- |
| `t_solve` (headline warmstart metric) | **yes** | its timer closes before the BC pass begins |
| `t_step_total`, end-of-run solve share | **no** | the pass lies inside the `maneuver!`-to-`maneuver!` window |

Because full-timestep share is a campaign deliverable (Ryan 2026-08-06), the
driver emits **`t_step_net = t_step_total − t_bcerr`** and prints its solve share
against net wall, quoting raw and the measurement total alongside. Measured R1
impact: **35% raw vs 45% net**. **Any timestep-share claim must use
`t_step_net`.** This was missed on the first pass and caught when Ryan asked how
the overhead had been treated.

**Not re-measured:** that 62% is from a *cold-start* R1 smoke with a wake growing
from zero. Phase 3 arms restart from a developed wake, where `t_solve` is larger
and the fraction should fall. Worth confirming from the first real run rather
than assuming.

---

## 6. Recommended next actions, in order

1. **Check HPC job progress and report — see §0.** This is Ryan's explicit
   instruction for the reset. Delegate to `hpc-monitor`; do not ssh inline.
2. **Offer Ryan a notebook entry** for the guard rebuild. Nothing has been
   written to `~/Dropbox/research/notebooks/journals/` this session, and this
   result is hard to reconstruct later: the finding is not "a threshold moved",
   it is that a validity guard was asking a question the solver never answers.
   Ask how much to write — natural scope: the problem, why the old guard was
   invalid (the residual-vs-solution tolerance confusion), Ryan's ruling, the
   design, the R1 verification table. **Offer, never write.**
3. **When Phase 2 Step A lands** (R1 first, ~08-26), the Phase 3 launch needs
   `PHASE=phase3` (which auto-enables the monitors `CT` depends on) and the
   frozen knobs from `results/phase1/multi/*`. On the first real rows, confirm
   `bcerr_certified` is true before trusting any BC verdict, and re-read
   `t_bcerr / t_solve` — the 62% figure is from a cold-start smoke and should
   fall with a developed wake (§5).
4. **Low-priority, deliberately not done:** the `phase_19` items are all
   complete; `benchmark/slurm/p2_unsteady.sh` has still never been submitted; the
   whole tree remains uncommitted (~100 modified, `benchmark/` untracked in full)
   and both cluster checkouts are rsync copies with no git anchor — flagged
   2026-08-22, still true. That last one is the item least well served by another
   context reset, but it is a separate conversation and Ryan's call.

---

## 7. Standing limits (unchanged)

- **No `sbatch` / `scancel` / `scontrol`** without Ryan asking in the moment.
- **No sync to `~/flowpanel-021/FLOWPanel.jl`** while Phase 1 jobs run.
- **No git commit.**
- **Notebook entries and BRAINSTORM item edits: offer, don't write** (the item
  edit this session was explicitly requested).
