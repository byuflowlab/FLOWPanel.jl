# Phase 0 — driver construction + local smoke + GS convergence gate

**Objective:** a working `examples/rotor_hover_ground_effect.jl` (fork of RHPC
@ b251071) whose IGE and OGE paths both run locally, with the block
Gauss–Seidel rotor↔ground coupling **verified to converge** (Ryan ruling 3 —
HARD GATE before any HPC submission).

## Driver design summary

- Fork, not include: RHPC is frozen for live 018; it also hard-codes ρ and R.
- New ENV knobs: `RHO` (1.16), `ROTOR_R` (0.1195), `GROUND_ENABLE` (true),
  `GROUND_H_R` (1.0), `GROUND_RADIUS_R` (4.0), `GROUND_PANEL_LENGTH_R` (0.15),
  `GROUND_PARTICLE_POLICY` (none|cull, default none), `TRUNC_RADIUS_R`
  (3.0 with ground, 1.5 without), `GS_LOG` (true), `GS_VERBOSE`,
  `GS_MAX_OUTER` (50), `GS_TOL` (1e-8).
- Ground = `FlatGround` disc at `x = x_rotor_plane + GROUND_H_R*R` (downstream
  = +axial; x_rotor_plane = mean axial node coordinate, printed in banner),
  normal toward the rotor; solved with `FlatGroundSolver` (O(N) diagonal,
  exact for a planar source sheet); wake entry `nothing`.
- The ground belongs to NO frame: `propagate_kinematics!`/`kinematic_velocity!`
  only touch frame dependents, so it stays static while frame 1 spins the
  rotor. (Do NOT use `add_frame!` as a child of the rotor frame — children
  are pushed into the parent's dependent list and would co-rotate.)
- GS instrumentation: the driver OVERWRITES the one-line
  `solve_formulation!(::VelocityThroughSources, ...)` method to forward
  `history`/`verbose`/`max_outer_iterations`/`outer_tolerance` into the
  EXISTING `solve!` outer loop (`src/FLOWPanel_solver.jl:1454`) — no new
  iteration logic, no src edits. Logs per-solve iteration count + final
  max-strength-delta to `*_gs_convergence.csv`; warns loudly on cap-hit
  (stock behavior is silent unless verbose).
- Ground diagnostics monitor (appended last, CSV `*_ground_diagnostics.csv`):
  per-step RMS/max ground U·n + below-ground particle count and Σ|Γ|.
- Metadata TOML extended with all ground/GS keys.

## Cases

| tag | knobs | where | status |
|---|---|---|---|
| smoke_ige | 40_40, NT=6, NREVS≈3, GROUND_ENABLE=true, defaults otherwise | local ≤4 threads | pending |
| smoke_oge | same, GROUND_ENABLE=false | local ≤4 threads | pending |

## Decision (pre-registered)

PASS Phase 0 iff ALL of:
1. Both smokes run to completion and write CT CSVs + metadata.
2. Banner checks pass: shedding counts sane, zero circumferential edges,
   RPM=6000 / ρ=1.16 / R=0.1195 echoed, ground geometry echoed.
3. Wake convects toward +axial (toward the ground) — confirmed via the
   below-ground census becoming nonzero and/or VTK inspection.
4. **GS gate:** every solve converges (`gs_nonconverged = 0`) with
   `gs_iters_max` well below the cap (≤ ~15 of 50); residual trajectory
   monotone on a spot check.
5. Ground tangency RMS(U·n) small vs tip speed and stable step-to-step.
6. OGE smoke CT trace is sane (finite, right order of magnitude).

FAIL on (4) → contingency chain (item file): diagnose, consider
`BackslashCoupled`; NO HPC submission.

## Exit criteria

Phase 0 PASS recorded in the Log below with the measured GS iteration stats
and tangency scale; then Phase 1/2 submission is unblocked.

## Log

- 2026-08-18 (staging): driver forked and parse-checked; smokes not yet run.
- 2026-08-18 smoke_ige (40_40, NT=6, 24 steps, 4 threads, defaults MAGVINF_PEAK=10):
  COMPLETED. First attempt crashed on wrong `ConvergenceHistory` field names
  (`residuals` → actual fields `iter`/`residual_internal`); fixed and rerun.
  Results:
  - Ground built: 4752 panels, x_rotor_plane/R = −0.0425, ground at
    axial/R = 0.9575 ⇒ downstream = +x CONFIRMED (particles crossed the plane
    late in the run: census 3 @ step 20 → 105 @ step 23; leave-be working).
  - **GS gate PASS:** 24/24 solves converged, 6–7 outer iterations each
    (cap 50), final max-delta ≤ 1e-8 every solve, `gs_nonconverged = 0`.
  - Tangency: RMS(U·n) max 1.14 m/s at the 10 m/s pulse peak (≈11% of onset —
    coarse-smoke scale: σ/R = 0.24 particles sit against the ground), final
    hover value 0.12 m/s = 0.16% of tip speed. Judge production scale on the
    56_57/45_185 runs.
  - Banner echoed RPM 6000 / ρ 1.16 / R 0.1195, shedding 35 edges/blade after
    root clip, zero circumferential edges. Note: uniform-D Das at NT=6 gives
    |Das|/R 0.82, Das/c 3–11 (out of the 014 band — an NT=6 smoke artifact,
    σ is 6× production; not a production concern at NT=36).
  - CT numbers at NT=6 are not meaningful; machinery-only smoke.
- 2026-08-18 smoke_oge (same knobs, GROUND_ENABLE=false): COMPLETED, exit 0.
  Single-body path clean, banner echoes "Ground plane DISABLED (OGE
  baseline)", CT finite and right order (0.068 ± noise at NT=6).
- **2026-08-18 PHASE 0 DECISION: PASS** (all six criteria met; GS gate
  decisive at 6–7 iterations of 50). HPC submission unblocked.
