# Phase 10 — Mesh Spanwise Refinement

**Objective:** converge CT̄ AND Γ̄(r/R) (TE μ-jump) in spanwise panel count at
fixed chordwise resolution (n_airfoil=185) and fixed converged wake settings.
Spanwise refinement is the natural resolver of the circulation distribution —
it adds Γ(r/R) sampling stations — and 2d Appendix G specifically warns that
chordwise-only refinement grows panel aspect ratio and the Dirichlet tangency
residual: span must be refined alongside chord for a defensible mesh claim.

## Prerequisite: mesh generation (local, Mac-only OpenVSP)

Rungs 60 and 80 spanwise don't exist yet. Generate with the 2d production cap
recipe (flat root + round tip, CapUMinTess=4):

```
scripts/generate_dji9443_mesh.sh 60 185 flat round 4
scripts/generate_dji9443_mesh.sh 80 185 flat round 4
```

Verify panel counts and cap topology (2d prescription); scp the `.msh` files
to the cluster `examples/data/` before submission.

## Cases (final settings Das\*/NT\*/σ\*/N=4 via env; velocity)

| tag | mesh | ~panels | dense G | time | note |
| --- | --- | --- | --- | --- | --- |
| final-settings run (exists from P4/P6–8) | 45_185_ct4 | 36,752 | 10.1 GB | — | span rung 1 |
| `p018_span60` | 60_185 ct4 (RHPC_MESH_FILE) | ~49k | ~19 GB | 48 h | span rung 2 |
| `p018_span80` | 80_185 ct4 (RHPC_MESH_FILE) | ~65k | ~34 GB | 48 h+ | contingent: only if 45→60 delta > threshold; watch memory (G + LU workspace + wake at 64G is tight — submit `--mem=128G` if the preflight estimate says so) |

Warm starts do NOT transfer across meshes (body construction differs) — all
rungs run cold with the full staged startup.

## Decision

|ΔCT̄| ≤ 0.5% and ε_Γ ≤ 1% between successive rungs (Γ̄ compared on the
common r/R grid by interpolation; report both TE-jump and slice estimators).
Span delta → error-budget term 14. Also report the Dirichlet tangency residual
trend (rms U·n) across rungs if cheaply available from the monitors — it
corroborates the 2d App. G aspect-ratio mechanism.

## Exit criteria

Spanwise rung deltas recorded (CT and Γ); budget term filled; mesh files
retained in `examples/data/` and on cluster.

## Log

- 2026-07-31 — Meshes GENERATED and shipped:
  `dji9443_20260731_60_185_capped_captess4.msh` and
  `..._80_185_capped_captess4.msh` (2d cap recipe), md5-verified on cluster.
  Panel-count / cap-topology / TE-detection verification still pending — do it
  before first submission (the RigidWakeBody shedding-from-constructed-cells
  invariant applies to any new mesh).
