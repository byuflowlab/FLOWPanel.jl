# Phase 4 — Blown wing: rotor + NACA 0015 wing multi-body benchmarks

**Status:** NOT STARTED. **Approval:** [ ]

## Prior-phase handoff

Phases 0–3 deliver: instrumented solvers, frozen ladder + settings, tuned costs, production
warmstart settings.

## Scope

Add a second body — the NACA 0015 rectangular wing from `examples/simple_wing.jl` (Ryan
2026-08-06) positioned in the rotor wake (blown-wing configuration; exact placement fixed at
phase start and recorded) — and repeat the Phase 2 benchmark protocol for the multi-body
solve. The rotor rotates relative to the wing, so cross-body influence changes every step.

Multi-body roster:

| Config | Per-step refactorization | Notes |
| --- | --- | --- |
| `BackslashCoupled` | **full block-`G` rebuild + `lu!` every step** — relative motion invalidates the monolithic factorization; LU is now per-timestep cost (timers restored in Phase 0) | run with `update_G=true` per step |
| Tuple block-GS (per-body `Backslash`) | per-body LUs frame-invariant; only cross-influence recomputed per outer iteration | **pre-registered headline asymmetry** vs coupled |
| `KrylovCoupled` | none (matrix-free) | currently has no preconditioner field — either extend (small, authorized) or benchmark unpreconditioned and say so |
| Tuple block-GS with per-body Krylov / FGS | per-body warmstart from Phase 3 | outer-iteration history recorded |
| FGMRES + FGS (multi-body) | as implementable | scope decided at phase start from Phase 0's architecture |

Method caveats to document (paper's methods section):

- The tuple block-GS path couples bodies **through velocity only**
  (`src/FLOWPanel_solver.jl:876`) — correct for Neumann bodies (both bodies here are Neumann,
  `watertight=false` → `VortexRing`), incorrect for Dirichlet. State this explicitly.
- Verify solution agreement between coupled and iterative paths at matched residual
  (Phase 1-style check) before timing; outer tolerance (`outer_tolerance`,
  `max_outer_iterations`) calibrated to the matched-residual target.
- Ladder: rotor rungs from the frozen ladder × a small set of wing panel counts (≥2), so
  scaling in total N is still measurable.

Benchmarks per rulings 5–8: setup vs per-step split (noting the coupled solver's LU now lives
in per-step), isolated solve primary + full-timestep share, both threading modes,
convergence-vs-time for iterative configs, memory, scaling fits.

## Exit criteria

1. New example/driver (rotor + NACA 0015 wing) committed with a smoke test.
2. Agreement check passed across multi-body configs at matched residual.
3. Full benchmark table + figures; the coupled-vs-iterative per-step refactorization
   asymmetry quantified across the ladder.

## Log

(dated entries appended newest-last)
