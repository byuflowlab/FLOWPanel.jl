# Handoff — plan BRAINSTORM 021 Phase 3 (warmstart)

Written 2026-08-22 night, at a context clear. Ryan's instruction: **enter plan
mode and prepare Phase 3.**

## Your task

Prepare an implementation plan for **Phase 3 — warmstart** (`phase_03_warmstart.md`),
the campaign's headline deliverable (ruling 12 / acceptance criterion 4).
Phase 3 is `NOT STARTED` and its approval box is unticked; "completing a phase
does not authorize the next," so **starting the phase needs Ryan's explicit go** —
plan it, do not launch it.

Read `log.md` newest-first for the narrative. Do NOT read the item file's
`## Current status` (stale, carries a banner) and do NOT read
`phase_01_consistency.md` (73 KB) or `phase_02_single_step_benchmarks.md` (41 KB)
end-to-end — grep them.

## Two things to raise with Ryan BEFORE planning Phase 3

Both were flagged this session and neither is done. Ask whether to clear them
first; they are cheap and both are pure risk reduction.

1. **HPC disk is at 117 G / 200 G and climbing** (~98 G earlier the same
   evening, so ~19 G in a few hours) with **8 BRAINSTORM 018 jobs still writing
   VTK**. `agent_policies/HPC.md` says to sweep at ~100 G *because sweeping late
   becomes a race*. The fix is to launch the `hpc-storage` subagent — it is
   restart-set-aware and keeps the newest 36 timesteps. Do NOT hand-roll `rm`
   (a bad sweep killed job 13036477 on 2026-08-04). This session did not run it
   only because of a standing "no subagents unless asked" instruction.
2. **The entire working tree is uncommitted** — 103 modified, 34 untracked,
   including all of `benchmark/` (untracked in full), the 021 near-field-cache
   and ILU fixes, this session's `filament_reg` instrumentation, unrelated
   018/022 work, and new `src/FLOWPanel_gpu_*.jl`. Both cluster checkouts are
   **rsync copies with no git anchor**, so a 5-day campaign is running on source
   that exists in exactly one place. Two things already went wrong from this
   class of drift: `p1_tune_s3.sh` was cluster-only and had silently lost the
   filament pin, and an `rsync` in a `for` loop silently transferred nothing.
   Suggest committing in logical chunks.

## State you are inheriting

**Ten 021 jobs queued, all still PENDING** (cluster at 128/136 nodes):

| job | id | walltime |
| --- | --- | --- |
| p2b-nearfield-R1 / R2 | 13306549 / 50 | 24 h |
| p2b-nearfield-R3 / R4 | 13306599 / 600 | 72 h |
| p1-table-R6-multi | 13306554 | 24 h |
| p1-table-R6-single-a → -jacobi | 13306556 → 13306557 (`afterany`) | 72 h |
| p1-tune-R7-s1 → s2 → s3 | 13306603 → 04 → 05 (`afterok`) | 72 h each |

Estimates: everything except R7 lands ~Tue 2026-08-25; the R7 tune chain ~5 days.
**Do not disturb the 8 fp-018 jobs.**

**OPEN TODO carried at the top of `log.md`:** those ten jobs snapshot-ran with
`FLOWPANEL_FILAMENT_REG=vatistas` (Slurm captures the batch script at submit
time, and the launcher default was flipped after submission). Ryan ruled **let
them run** — their numbers are family-invariant — but their family attribution
must be relabelled when they finish. Checklist is in `log.md`.

## Context you must NOT re-derive

- **The filament-family question is CLOSED.** Two channels: **Channel A** = the
  filament kernel (`_bound_vortex_velocity/_gradient`,
  `src/FLOWPanel_elements_fmm.jl:938-1067`), reached only on velocity/gradient
  passes where filaments exist → moves **numbers**, wake-carrying runs only.
  **Channel B** = `radius_inflation` (`:1130-1142`), dispatches on element type
  so it fires for any `VortexRing` body → moves the near/far split, i.e.
  **cost only, never values**. Measured, not assumed.
- Phase 1 / 2b are family-invariant (R1 A/B bit-identical incl. iteration count).
  **The FMM tuning runs do NOT need redoing** — they run at
  `core_size_panel = R·1e-10`, where the family moves the buffer by 0.0017% of
  the smallest panel at R1 and 0.015% at R7.
- The Gaussian rerun is DONE. `agreement.csv` was predicted sensitive and
  **measured not to be** (CT identical to 8 decimals) — restored unchanged.
  `unsteady.csv` rerun under Gaussian. The staleness certification arm
  `fgs_geometry_order_fixed_8step` reproduced **bit-for-bit**, proving it was
  already Gaussian, so the Phase-3 trust gate (`phase_02...md:313`: "NO Phase-3
  unsteady FGS row is trusted until the discriminator runs") **is satisfied**.
- FGS unsteady staleness + `simulate!` lifecycle ordering: CONFIRMED and FIXED
  2026-08-21. ILU scales ~N^1.5. Both written up; do not re-measure.
- `filament_reg` is now recorded by `assert_and_banner()` (`benchmark/common.jl`)
  and required by `validate_runs_csv`; `decision_rules.md` carries the schema
  amendment.

## What Phase 3 needs — start here

`phase_03_warmstart.md` is the spec (small, ~3.5 KB, read it in full). Key facts
already established:

- **The driver already exists**: `benchmark/rotor_hover_solver_unsteady.jl` says
  in its own header "Phase 3 warmstart runs reuse this file, ruling 12". It
  includes `examples/rotor_hover_pressure_comparison.jl` with `RHPC_SETUP_ONLY=1`,
  so the case construction is NOT duplicated. `benchmark/slurm/p2_unsteady.sh` is
  its launcher — **PREPARED 2026-08-18 but NEVER SUBMITTED**, so no unsteady HPC
  rows exist at all.
- **Phase 3's numeric prerequisites are safe.** Frozen knobs come from
  `results/phase1/multi/*` and the cold single-step baselines from
  `phase2.csv` — all wake-free, hence family-invariant.
- **~~Known blocker to design around~~ — SUPERSEDED 2026-08-23, corrected in
  place 2026-08-25.** This blocker is **CLOSED**: `FGSSolver.niter` /`.solved`
  and the per-step lifecycle (`step_niter_first`, `step_nsolves`,
  `step_solved`) all exist and are unit-tested; `unsteady.jl` records FGS
  counts alongside Krylov. The phrase "pending callback plumbing" no longer
  appears anywhere in `src/` or `benchmark/`. **Do not re-implement it.** See
  `phase_19_fgs_niter_followups.md` and the 2026-08-23 log entry. Only
  `Backslash` still records `niter = -1`, which is a deliberate null result,
  not a gap.

  ~~Original text: FGS per-step iteration capture is not plumbed —
  `unsteady.csv` records `niter=-1` for Backslash and FGS, flagged in the driver
  header as "pending callback plumbing (Phase 3)". Phase 3's headline metric is
  *iterations*-to-target per step, so this must be solved. Krylov already
  reports `solver.niter`.~~
- Phase 0 delivered the `x0` / `warmstart` plumbing (`warmstart` field, `x_prev`,
  positional-x0, `KrylovCoupled` persistent `x`), landed and unit-tested. FGS
  already owns history/extrapolation (`project_solution!` / `save_solution!`,
  `src/FLOWPanel_solver.jl:758-785`); the spec says Krylov's extrapolated guess
  should reuse the same finite-difference coefficients for comparability.
- Phase 3 is **unsteady/wake-on** — the one place in 021 where the filament
  family is live. It will be Gaussian from birth (all six
  `benchmark/slurm/*.sh` now default to gaussian), so it carries no rerun debt.

## Open questions worth putting to Ryan while planning

- **The R7 row does not fit the partition.** `m12` MaxTime is 3-00:00:00;
  `p1-table-R7-single-a` extrapolates to ~100 h and the ledger already has
  `R7-single-jacobi` at ~69 h. The tune was split via the new `STAGES` gate in
  `p1_tune.sh`, but the R7 *table* jobs need their own splitting strategy
  (per-config is the obvious move) before they can be submitted. Not urgent —
  they are blocked on 13306605, ~5 days out.
- Whether Phase 3 runs at R1 only or across rungs, and how many steps are needed
  to reach the "wake-developed regime" the spec requires (startup transient
  excluded from averages).
- A same-code family A/B on the unsteady driver produced one solid number —
  particles 328 (gaussian) vs 329 (vatistas), a real trajectory change — but its
  **timing columns are noise-confounded and must not be cited** (single runs
  against a protocol requiring min-of-k k≥5; families ran in separate batches;
  signs disagree between configs). Artifact:
  `unsteady_family_ab_20260822.csv`. If per-step *cost* under the two families
  matters to Phase 3's baseline, it needs a proper HPC measurement.

## Also flagged, not fixed

`_U_boundvortex` (`src/FLOWPanel_elements_fmm.jl:1559-1615`) is hardcoded
Vatistas and ignores the global, while `_U_boundvortex_gradient` (`:1627`)
honours it. Latent only — it fires on the `semiinfinite_wake=true` path, which
021 does not use.
