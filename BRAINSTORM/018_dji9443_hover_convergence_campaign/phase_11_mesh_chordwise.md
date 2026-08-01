# Phase 11 — Mesh Chordwise Refinement

**Objective:** converge CT̄ AND Γ̄(r/R) (TE μ-jump) in chordwise panel count
at fixed spanwise resolution (45) and fixed converged wake settings. The
steady-state chordwise ladder is already converged (Phase 2c: Dir ∫Γ changes
≤0.22%/rung from n=121; 006: production-mesh refinement moved CT less than
cycle scatter) — this phase establishes the same under the UNSTEADY particle
wake, which has never been measured.

## Cases (final settings Das\*/NT\*/σ\*/N=4 via env; velocity; meshes exist)

| tag | mesh | ~panels | time | note |
| --- | --- | --- | --- | --- |
| `p018_chord145` | `45_145_ct4` (known key) | ~28.8k | 24 h | coarse rung |
| final-settings run (exists) | 45_185_ct4 | 36,752 | — | middle rung |
| `p018_chord249` | `dji9443_20260725_45_249_capped_captess4.msh` (RHPC_MESH_FILE) | ~49.5k | 48 h | fine rung; dense G ~20 GB |

All cold starts (no cross-mesh warm start). Disclose alongside results: the
Dirichlet tangency thrust contamination GROWS with chordwise-only refinement
(2d App. G: 0.84% of CT at n=145 → 3.24% at n=249, aspect-ratio driven) — the
n=249 rung's CT delta must be interpreted with that systematic in view, which
is exactly why Phase 10 (spanwise) runs separately.

## Decision

|ΔCT̄| ≤ 0.5% and ε_Γ ≤ 1% between successive rungs. Chord delta →
error-budget term 15. If 145→185 and 185→249 deltas bracket zero within
threshold, the production mesh stands; a monotone drift larger than threshold
at 249 triggers a look at the tangency contamination before any mesh change.

## Exit criteria

Chordwise rung deltas recorded (CT and Γ); budget term filled.

## Log

(append dated entries here)
