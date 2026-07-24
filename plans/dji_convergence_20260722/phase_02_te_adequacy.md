# Phase 2 — TE Adequacy

## Purpose and required context

Decide whether the sharp trailing edge on capped meshes makes the low-order
Dirichlet circulation locally inaccurate or excessively sensitive to nearby
panel influence. Keep these failure modes separate:

- poor conditioning of the full lifting operator;
- an ill-conditioned Kutta observable
  `Gamma_TE = mu_upper - mu_lower`;
- small centroid residuals but poor boundary-condition satisfaction between
  control points.

Read the repository instructions, the dashboard, the top snapshot in the
[Phase 2 log](../../logs/dji_convergence_20260722/phase_02_te_adequacy.md), and:

- `agent_policies/WORKFLOW.md`
- `agent_policies/TESTING.md`

Do not begin unless Phase 1 approval and its required handoff decisions appear
in the Phase 2 log.

## Baseline diagnostic

Run locally with steady, no-free-wake solves. Start on the new 40-series
capped/Dirichlet mesh; add the new 57-series only after the diagnostic works.
Keep the body-panel kernel offset distinct from the surface-to-particle target
offset used in later phases.

For every paired shedding edge, record:

- upper/lower control-point separation normalized by local panel length;
- trailing-edge dihedral angle from the two panel normals;
- panel area, aspect ratio, and centroid distance from the shedding edge;
- minimum control-point-to-triangle distance to every non-self panel;
- identities and roles of the nearest paired cross-TE, same-side upstream,
  opposite-side, and cap panels.

Decompose the Dirichlet lifting operator into body-only `B` and attached-wake
Kutta terms, `G = B + W*C`. Assemble `B` using
`suppress_attached_wake`. For upper and lower TE rows, rank source panels by
both `abs(G[i,j])` and `abs(G[i,j] * mu[j])`, recording:

- self and paired cross-TE influence;
- cumulative nearest one, three, and five contributions;
- body-only and attached-wake contributions separately;
- upper/lower row similarity and normalized row differences;
- the local upper/lower 2-by-2 influence block.

## Observable sensitivity

For each shedding edge `e`, define `c_e` such that
`Gamma_e = c_e' * mu`, solve `G' * z_e = c_e`, and report nearest-entry
sensitivity as:

```text
abs(z_e[i] * G[i,j] * mu[j]) / max(abs(Gamma_e), Gamma_scale)
```

Confirm the adjoint prediction by scaling only the mutual upper/lower entries,
then the nearest three to five TE contributions, by `0.99` and `1.01`.
Re-solve and measure local, integrated, thrust-weighted, and `r/R >= 0.7`
circulation. Label this an operator-fragility diagnostic, not a physical
modeling change.

## Between-centroid and offset checks

On the final two or three chordwise rows, evaluate the solved total field at
multiple barycentric locations per panel and at short normalized
interior/exterior normal offsets safely outside the regularization scale.
Compare total-potential residual and exterior normal-velocity leakage near the
TE with surface-wide median and maximum values. Repeat at both mesh
resolutions and test whether local errors decrease.

Log `kerneloffset_panel` independently and test a small decade around its
nominal value. It must remain much smaller than local panel length and the
upper/lower control-point separation. Never describe this as a
`kerneloffset_targets` or wake-core sensitivity.

## Conditional thickness study

Trigger this subphase if any condition holds:

- a 1% nearest-influence perturbation changes integrated or outboard
  circulation by at least 1%;
- circulation-observable sensitivity grows materially from 40 to 57 series;
- TE off-collocation residuals are much larger than surface-wide values or do
  not decrease under refinement;
- circulation changes by at least 1% under the body-panel offset diagnostic.

If triggered, generate matched watertight capped/Dirichlet geometries targeting
`t_TE/h_TE = 0, 0.05, 0.1, 0.2, 0.4` where clean meshing permits. Preserve
camber, thickness distribution, and paneling outside the final one or two
chordwise rows; close each blunt TE with a narrow strip and retain the same
shedding/wake construction. Also refine only the final two or three chordwise
rows on the sharp geometry.

Measure whether local/integrated circulation, Kutta-Joukowski thrust,
observable sensitivity, and off-collocation residuals approach stable limits
as TE thickness decreases and local resolution increases.

## Decision and exit gate

The sharp TE is adequate only if:

- circulation is stable to the larger of 1% or twice extraction variability;
- off-collocation errors decrease under refinement;
- circulation sensitivity does not grow without bound.

Otherwise report it unresolved and propose a geometry, local-refinement, or
higher-order remedy before proceeding.

Preserve a compact per-edge CSV and summary covering geometry, influences,
adjoint sensitivities, perturbation responses, residuals, and any
thickness/refinement results. Review the methods, results, and conclusion for
consistency; update the log and dashboard; then stop for explicit Phase 3
approval.

