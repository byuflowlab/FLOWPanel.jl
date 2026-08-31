# BRAINSTORM 021 — Phase 19: FGS iteration capture — status correction and follow-ups

> **CLOSED 2026-08-25. All five items are done.** Do not read this as a task
> list. Start from `phase_20_context_reset_prompt.md` instead; this file is
> retained for the history in §5 (two superseded guard designs and why).

**Written 2026-08-25. Read this file, then `phase_18_context_reset_prompt.md`
for what is running. You should not need anything else to start.**

## READ THIS FIRST — the task you were probably handed is already done

Ryan asked for "plan the plumbing for FGS per-step iteration capture,
implement, and update the item." **The plumbing is already implemented.** Do
not re-implement it. This file exists because a stale planning document says
otherwise and has now misled at least one agent.

**The stale source:** `phase_03_launch_prompt.md:103` lists as a "Known blocker
to design around" that FGS per-step iteration capture "is not plumbed —
`unsteady.csv` records `niter=-1` for Backslash and FGS ... pending callback
plumbing (Phase 3)". That document was written the night of **2026-08-22** and
was **superseded on 2026-08-23**. It has not been corrected in place.

**Verified current state of the source (2026-08-25):**

| what | where | state |
| --- | --- | --- |
| `FGSSolver.niter` field | `src/FLOWPanel_solver.jl:1383` | EXISTS — "GS sweeps performed by the last solve" |
| `FGSSolver.solved`, `.stats::SolveStepStats` | `:1384-1385` | EXISTS |
| counting callback into `FastMultipole.solve!` | `:1674-1709` | EXISTS, chains a user callback |
| off-by-one correction | `:1709` | `niter = solved ? max(last_iteration[]-1, 0) : last_iteration[]` |
| per-step lifecycle (`SolveStepStats`, `note_step_solve!`, `publish_step_stats!`, `step_niter_first`) | `:53-120`, `:1431-1466` | EXISTS, under section header `PER-STEP SOLVE STATISTICS AND WARM-START LIFECYCLE (BRAINSTORM 021 Phase 3)` |
| unsteady driver reads it | `benchmark/rotor_hover_solver_unsteady.jl:241` | `_niter_of(s) = s isa KrylovSolver \|\| s isa FGSSolver ? s.niter : -1` |
| FGS convergence history sidecar | `:1728-1744` `_solve_history!` | EXISTS, metric `:fgs_maxabs` |
| the phrase "pending callback plumbing" | `src/`, `benchmark/` | **ABSENT** — appears only in the stale doc |

The 2026-08-23 log entry records the working result: `niter_first` of
6.00 ± 0.00 (cold), 4.75 ± 0.50 (prev), 5.00 ± 0.00 (extrap) at R1/FGS, with
`nsolves == 1` enforced on every scored row. It also fixed a real measurement
bug — the driver had been sampling `solver.niter` *after* the step, returning
the last inner solve, already warmed by its predecessors. `niter_first` is the
warmstart metric; legacy `niter` is now a diagnostic. Tests are in
`test/runtests_unit_solver.jl`.

---

## How to work with Ryan — TWO MANDATORY CHECKPOINTS

This is not an autonomous task. Ryan works interactively and expects to be
consulted. Do not run items 1-5 end-to-end and present a finished result.

**Checkpoint A — BEFORE you edit anything.** Do the read-only survey first:
grep the sibling drivers for the `-1` sentinel (item 1), confirm the source
state in the table above still matches, and locate the guard constant in
`phase3_analysis.jl` (item 5). Then report to Ryan:

- exactly which files you propose to touch and the one-line diff for each;
- the sibling-driver survey result, because "fix all of them or none" is a
  judgement call about column semantics that is HIS, not yours;
- anything in the table above that does NOT match what you find — if the
  source has moved, stop and say so rather than adapting silently.

Wait for his go-ahead before the first edit.

**Checkpoint B — BEFORE updating the BRAINSTORM item (item 4).** Show him the
proposed item text and wait. Item edits are durable campaign record. The same
applies with more force to any notebook entry: **offer, never write.**

**Bring back, as an explicit decision request, not a footnote:**

1. ~~**extrap still fails at 1e-7** (2.4e-7, ~2x over). State the number plainly
   and ask whether it is a real finding about order-1 extrapolation or a defect
   in the test. Do NOT loosen the threshold further to make it pass.~~
   **ANSWERED 2026-08-25 — see the superseded banner atop §5. It is neither: the
   TEST is mis-specified (a residual tolerance read as a solution tolerance).
   2.4e-7 sits ~3 decades inside the solver's own solution ambiguity; extrap is
   fine. Do not report it as DIVERGED.**
2. The two structural weaknesses in the trajectory test (summed checksum
   permits cancellation; a fixed bar ignores step count).
3. Anything you found that contradicts this document.

**Hard limits — no exceptions without Ryan asking in the moment:**

- **No `sbatch`, `scancel`, or `scontrol`.** Seven jobs are running (see
  `phase_18_context_reset_prompt.md` §2). Nothing here needs the cluster.
- **No syncing to `~/flowpanel-021/FLOWPanel.jl`.** These edits take effect
  only on a future re-run; pushing them mid-campaign buys nothing and risks
  the running jobs' reproducibility story.
- **No git commit.** The tree is uncommitted by long-standing choice; that is
  a separate conversation.

---

## Your task

Five items. **Items 1 and 5 are code; 2, 3 and 4 are documentation.** Item 5
(the threshold) is independent of 1-4 — do it in either order, but do not let
it absorb the caveat attached to it.

### 1. The one genuine gap — `phase1_agreement.jl` still hard-codes `-1`

`benchmark/rotor_hover_solver_phase1_agreement.jl:439`:

```julia
niter = solver isa pnl.KrylovSolver ? solver.niter : -1
```

This predates the Phase 3 FGS work and was never updated. It is why today's
Phase 1 rows show `fgs` with `niter=-1` while `krylov_ilu` shows 11. Since
`FGSSolver.niter` now exists, this discards a real, useful cost number — FGS's
sweep count for a cold solve is directly comparable in spirit to a Krylov
iteration count, and Phase 1's table is exactly where a reader wants it.

The obvious fix mirrors the unsteady driver:

```julia
niter = solver isa pnl.KrylovSolver || solver isa pnl.FGSSolver ? solver.niter : -1
```

**Before you make it, know these three things:**

- **CORRECTED 2026-08-25 — the sibling survey was run; this paragraph's
  guesses were wrong.** The original text (kept below, struck) named
  `phase1_table.jl:254` and "the phase2 drivers" as candidates. Neither is a
  sentinel site, and the one real sibling was not named. Measured result:

  | file:line | expression | verdict |
  | --- | --- | --- |
  | `rotor_hover_solver_phase1_agreement.jl:439` | `solver isa pnl.KrylovSolver ? solver.niter : -1` | **sentinel — fix** |
  | `rotor_hover_solver_phase1_solvetime.jl:175` | same | **sentinel — fix** (never named in this doc) |
  | `phase1_table.jl:149`, `:254` | `solver.niter` unconditional | not a sentinel. `:254` is the `fgmres_fgs` OUTER count — the item-3 caveat case. Editing it would corrupt that row's semantics. |
  | `phase1_table.jl` `fgs` branch (~`:185`) | `niter = length(hist)` | not a sentinel |
  | `phase2.jl:352/412/453`, `phase2_tune.jl:415`, `phase2b_nearfield_cache.jl:162`, `phase1_fgsprecond.jl:124`, `phase1_fgs_check.jl:246`, `smoke.jl` (×4) | `solver.niter` unconditional | not sentinels — already record FGS sweeps |

  So "fix all or none" is NOT an open question about column semantics. Every
  driver in `benchmark/` except those two already writes `solver.niter`
  unconditionally; the column is already non-homogeneous everywhere (see item
  3). The two sentinel files are the outliers. Fix both; edit nothing else.

  ~~Original text: Grep for `isa pnl.KrylovSolver ? solver.niter : -1` (and
  near variants) across `benchmark/` — `rotor_hover_solver_phase1_table.jl:254`
  and the phase2 drivers are candidates. Fix them consistently or not at all; a
  column that means "sweeps" in one file and "unavailable" in another is worse
  than one that is uniformly unavailable.~~
- **No schema break, but no retroactive fill.** The `niter` column already
  exists, so nothing in the header changes and `p021_merge.jl`'s header-drift
  guard will not fire. But R1–R5 are already written with `-1`, and R6/R7 are
  RUNNING right now under the old code — a mid-flight edit does not reach
  them (Slurm spooled the script at submit; the driver is read at job start).
  So this generation will be uniformly `-1` for FGS regardless. The fix takes
  effect only on a future re-run. Say so wherever you record it.
- **Do not sync this to the cluster while Phase 1 jobs are running** unless
  you have a reason to. There is no benefit until a re-run, and the campaign's
  discipline is to change as little as possible mid-flight.

**Verification:** `agreement.csv` is produced by a full HPC solve, so do not
try to regenerate one locally to test. Instead run the narrowest unit test per
`agent_policies/TESTING.md` covering `FGSSolver.niter`, and confirm by
inspection that the expression matches the unsteady driver's.

### 2. Correct the stale document IN PLACE

Edit `phase_03_launch_prompt.md:103` so it no longer states the blocker is
open. Do **not** delete the paragraph — mark it superseded with a pointer to
this file and to the 2026-08-23 log entry, matching the campaign's convention
of marking artifacts void in place rather than removing them (see how the
2026-08-22 smoke table was voided). A future agent planning Phase 3 will read
that document first, exactly as happened here.

### 3. Record a caveat that affects how the Phase 1 table is READ

`fgmres_fgs` reports a real `niter` (10 at R3) but it is the **outer FGMRES
iteration count**, not inner FGS sweeps. The `FGSPreconditioner` holds its own
`FGSSolver` with `sweeps` FIXED and `tolerance=0`
(`src/FLOWPanel_solver.jl:1772-1787`), so its per-apply sweep count is a
configured constant and is not surfaced anywhere.

This means the `niter` column in the Phase 1 agreement table is **not
homogeneous**: Krylov iterations, FGS sweeps (once item 1 lands), and outer
FGMRES iterations are three different units. Anyone comparing them as "cost"
will be wrong. Put this in `decision_rules.md` next to the existing iteration
counting amendment (`decision_rules.md:81`, "Iteration counting (Phase 3
amendment, 2026-08-22)"), which already governs this exact topic.

### 4. Update the item

Per Ryan's instruction, update the BRAINSTORM item to reflect that the FGS
capture is closed. Keep it to the item's own conventions; do not restate this
file's contents there.

---

## 5. Trajectory-identity threshold — SUPERSEDED 2026-08-25 (later same day)

> **SUPERSEDED. Do not implement the 1e-7 edit described below, and do not
> implement the relative-L2 design that briefly replaced it.** Ryan asked why
> extrap missed by 2.4x; the investigation found the guard itself was
> mis-specified, and Ryan then ruled that the whole approach of measuring
> SOLUTION divergence should be replaced by measuring the BOUNDARY CONDITION.
>
> **WHY THE OLD GUARD WAS WRONG (this part stands).** The solver's acceptance
> test is `max|b - Ax| <= tol_abs` — ABSOLUTE and independent of the initial
> guess (`FastMultipole/src/solve.jl:1222`; `residual!` at `:1530-1566` returns
> `mae`, the `mse` name is vestigial per `:1210`). Cold and warm arms call the
> identical `FastMultipole.solve!` (`src/FLOWPanel_solver.jl:1694`); the only
> difference is `project_solution!` at `:1659`. Both therefore land at different
> points of the SAME residual ball, whose solution-space diameter is
> `||A^-1|| * 2*tol_abs`. The old bound assumed `||A^-1|| = 1` — it read a
> RESIDUAL tolerance as a per-strength SOLUTION tolerance. Phase 1 measures the
> amplification at ~12x (R1), putting the expected spread near 2e-4, so the
> observed 2.4e-7 was ~3 decades INSIDE the solver's own ambiguity.
> **extrap was never diverging.** A second defect: `sum(abs, x)` is one
> aggregate over all panels and cancels equal-and-opposite changes.
>
> **RYAN'S RULING 2026-08-25 — measure the BC, not the solution.** "There is no
> need to measure the divergence of the solution itself between warm and cold
> solves; rather, the residual should be measured... it's easier to show that
> we're satisfying the boundary conditions." This removes the entire
> conditioning argument: instead of comparing arms and justifying the gap, each
> arm is checked against the tolerance IT promised.
>
> **IMPLEMENTED (see `decision_rules.md` -> "BC satisfaction guard"):**
> `bc_error!` (`benchmark/common.jl`, 021 ruling 3) already did exactly this and
> needs only ONE influence pass — for this Dirichlet body the control-point
> potential IS the BC residual `phi_sigma + G_mu x`. The driver now records per
> step: `bcerr_{max,min,q1,med,q3,rms,tol,eps,certified}`, `t_bcerr`,
> `n_particles`, `CT`. Accuracy comes from REQUESTING an absolute FMM error
> `BCERR_SAFETY` (0.1) x the arm's own tolerance and reporting certification —
> Ryan asked for ~10x sharper than the solve; this delivers it by construction
> rather than by guessing an expansion order. Identity is REPORTED, not
> asserted (CT + per-step particle count): the wake is chaotic, so exact
> trajectory identity is not expected.
>
> **RETRACTION — `tol_abs = 1.77e-8` was CORRECT.** An earlier revision of this
> banner claimed the value "could not be reproduced" because
> `fgstune_margin_verify_R1.csv` gives 9.93e-8 at R1/tau=1e-6. That was wrong:
> that row is for a DIFFERENT knob set (p=6, mac=0.3, leaf=150, inner=5) than
> the Phase 3 winner. A live R1/fgs run on 2026-08-25 emitted
> `tol_abs=1.7658196265827666e-8`. The number in the original document was
> right; the check against it was not.
>
> **Verified on a live R1/fgs run (cold + prev, 5 steps each), 2026-08-25:**
> every step certified, `max |phi|/tol` 1.007 (cold) and 0.994 (prev) — i.e.
> both arms sit right at the bar they promised, as they should. `Delta CT` between
> the arms is **1.3e-5 %**. Pass band is `1 + eps/tol` (the measurement's own
> certified uncertainty), not a bare 1.0: the arm stops on its own coarser
> internal estimate, so a sharper independent evaluation legitimately differs by
> about the error requested of it.
>
> Everything below is retained as the record of the superseded reasoning.

**Ryan ruled 2026-08-25: 1e-7 is sufficient.** Implement it in
`phase3_analysis.jl` (the guard that flags a warm arm as DIVERGED). Record the
ruling and the reasoning below wherever the campaign records rulings.

**Exactly what the number is a tolerance ON.** Not wake trajectories, not
particle positions:

- `benchmark/rotor_hover_solver_unsteady.jl:364` —
  `strength_checksum = sum(abs, rotor.strength)` — the L1 norm of the BODY's
  panel-strength vector, written once at the final step (`i == nsim`).
- `benchmark/phase3_analysis.jl:171` —
  `rel = abs(a.checksum - bi.checksum) / abs(bi.checksum)`, warm arm vs cold
  baseline.

So the guard is a relative tolerance on ONE SCALAR: the L1 norm of ~8016 panel
strengths at the end of the run. Wake identity is covered separately and
discretely by `n_particles_end` (407 in all three arms).

**Why 1e-8 was wrong.** Each solve converges to `tol_abs = 1.77e-8` ABSOLUTE
per strength, and the checksum sums 8016 of them to ~114596 (mean |strength|
14.30). The worst case a single step's tolerance ball can produce is every
strength moving the full tolerance in the same direction:
`8016 * 1.77e-8 / 114596 = 1.24e-9` relative. Random signs give
`sqrt(8016) * 1.77e-8 / 114596 = 1.38e-11`. Recomputed from the 2026-08-23
checksums (cold 114595.9408684665, prev 114595.9503456424, extrap
114595.9683812670): **prev 8.27e-8, extrap 2.40e-7** — 67x and 194x that
worst case, consistent with amplification across 7 marched steps rather than a
physics difference. A 1e-8 bar sits below that, so it flags correct runs.

(The 08-23 log entry quotes prev as 7.9e-8; the checksums in the same table
give 8.27e-8. Immaterial to the ruling, but use the recomputed values.)

**CARRY THIS CAVEAT — 1e-7 does not make both arms pass.** prev (7.9e-8)
passes; **extrap (2.4e-7) still fails by ~2x**. So after this change the guard
reports the campaign's own order-1 warm start as DIVERGED. That is EITHER a
real finding about order-1 extrapolation OR evidence the test itself is
mis-specified. **Do not paper over it** by loosening further on your own
authority, and do not report extrap as passing. Surface it to Ryan with the
number.

**Two structural weaknesses in the test, for Ryan when that comes up.** Raise
them; do not act on them unilaterally.

- **One scalar aggregate can under-report.** `sum(abs, ...)` has no SIGN
  cancellation, but it is still a single number over 8016 panels: a change
  that raises `|strength|` on one panel and lowers it by the same amount on
  another leaves the checksum untouched. Phase 1 answers the analogous
  question with relative L2 over the full strength vector (max pairwise
  1.73e-5 at R1, 1.90e-5 at R3 across seven solvers), which is element-wise
  and cannot compensate. Adopting the same metric here would make the two
  phases directly comparable.
- **A fixed threshold ignores step count.** If the difference amplifies through
  the wake, the correct bar scales with marched steps — a constant is too loose
  early and too tight late. The 08-23 data shows the extrapolation residual
  holding one sign and decaying monotonically (-8.01e-4, -2.83e-4, -2.06e-4,
  -1.45e-4) across the wake-developed window, which is the behaviour a
  step-aware bar would want to exploit.

**Verification.** Both arms match cold on particle count (407); keep that as a
separate, independent assertion — it is a discrete check that a relative
tolerance cannot substitute for.

---

## Still Ryan's call — do not decide

- Starting Phase 3 at all: `NOT STARTED`, approval box unticked.
- Every `sbatch`.
- Any notebook entry — offer, don't write.
- Whether extrap's continued failure at 1e-7 is a finding or a test defect.

---

## Context you should NOT re-derive

- Phase 3's numeric baselines come from `phase2.csv` (Step C output), which has
  never been run and is gated on Step A tuning landing — 1 to 7 days out as of
  2026-08-25. Frozen knobs from `results/phase1/multi/*` are fine and
  family-invariant.
- `benchmark/slurm/p2_unsteady.sh` exists but has **never been submitted**, so
  there are no unsteady HPC rows at all.
- Phase 0 delivered the `x0`/`warmstart` plumbing; it landed and is unit-tested.
- Phase 3 is Gaussian-from-birth and carries no filament rerun debt.
- `unsteady_family_ab_20260822.csv` timing columns are noise-confounded and
  **must not be cited**.
- The whole working tree is uncommitted (~100 modified, 34 untracked,
  `benchmark/` untracked in full) and both cluster checkouts are rsync copies
  with no git anchor. Flagged 2026-08-22, still true.
