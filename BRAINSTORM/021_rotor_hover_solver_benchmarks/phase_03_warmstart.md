# Phase 3 — Warmstart: previous-solution and Taylor-extrapolated initial guesses

**Status:** NOT STARTED. **Approval:** [ ]

**Headline deliverable (ruling 12, Ryan 2026-08-12):** how much warmstart helps *each*
solver is one of the campaign's key published results, not a side study. Phase 2's
benchmarks run strictly cold so this phase's cold-vs-warm comparison is controlled: same
frozen Phase 1 settings, same rungs, same residual targets — the only variable is the
initial guess.

## Prior-phase handoff

Phase 0 delivers Krylov `x0` plumbing (`warmstart` field, `x_prev`, positional-x0 launch —
landed and unit-tested, incl. `KrylovCoupled` via its persistent `x`); Phase 2 delivers
tuned cold baseline settings and costs.

## Scope

Measure the per-step savings from warmstarting each iterative solver in the **unsteady**
rotor-hover simulation (wake on), over enough steps to reach the wake-developed regime (steady
periodic iteration counts; startup transient excluded from averages).

Warmstart matrix:

| Initial guess | Backslash / BackslashCoupled | Krylov (GMRES / GMRES+precond / FGMRES+FGS) | KrylovCoupled | FGS |
| --- | --- | --- | --- | --- |
| None (zero) | null control — direct solves are guess-independent; documents that their per-step cost cannot benefit | baseline | baseline | baseline |
| Previous solution | — | via Phase 0 `x0` | via Phase 0 `warmstart=true` (`solver.x` seed) | `solution_history_length=1`-equivalent |
| Taylor/polynomial extrapolation | — | extrapolated `x0` (reuse or mirror `project_solution!` logic) | same | existing `project_solution=true`, sweep `project_solution_order` |

Pre-registered expectation (to test, not assume): Backslash null; Krylov modest; FGS large.
Phase 0 smoke anecdotes to test properly: warmstarted GMRES steps late in a 3-step sim
converged in ~0–2 iterations vs a 290-iteration cold solve — if that holds in the
wake-developed regime it reorders the per-step cost ranking entirely.

**Reported metrics per (warmstart type × solver × rung × threading mode)** — ruling 12:
- iterations-to-target and wall-time-to-target, wake-developed mean ± spread, vs cold;
- **break-even step count**: steps needed for cumulative warm savings to amortize any
  setup/bookkeeping overhead of the warmstart machinery (history storage, projection cost);
- interaction with the block-GS outer loop for tuple solvers (a warmstarted inner solver
  changes outer-iteration convergence — record outer counts too).

Notes:
- FGS already owns the history/extrapolation machinery
  (`project_solution!`/`save_solution!`, `src/FLOWPanel_solver.jl:758-785`); Krylov's
  extrapolated guess should reuse the same finite-difference coefficients for comparability.
- Warmstart insertion must survive `set_strengths!` zeroing (:1083) — verified in Phase 0 W3.
- Convergence criterion stays the Phase 1 matched-residual target: report iterations-to-target
  and time-to-target reductions vs. the zero-guess baseline, per rung, both threading modes.
- Guard against silent accuracy drift: spot-check RMS/max true residual and CT at the accepted
  settings.

## Exit criteria

1. Iteration-count and wall-time reduction table (warmstart type × solver × rung × threading
   mode), wake-developed averages with spread.
2. Recommended production warmstart setting per solver, with the order sweep for the
   extrapolated variants.

## Log

(dated entries appended newest-last)
