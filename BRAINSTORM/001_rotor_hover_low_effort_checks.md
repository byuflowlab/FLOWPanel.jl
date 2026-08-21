# Rotor Hover CT: Low-Effort Implementation Checks

## Context

Prior audit has already checked geometry, RPM, density, CT convention, radius/diameter interpretation, and handedness/sign. Treat those as background, not active hypotheses. Do not spend more time on the ParaView symmetry/handedness path except as context for interpreting results.

## Hypothesis

The CT discrepancy may come from an implementation choice or force-recovery path that can be exposed without changing the physical model. The useful checks are the ones that either move CT materially or create disagreement between otherwise equivalent evaluations on the identical solved state.

## Checks

- Compare `DirectBackend` and FMM on the same rotor-hover pressure/velocity pass. If they disagree beyond expected numerical tolerance, isolate whether the error is in induced velocity, pressure recovery, or force integration.
- Evaluate Bernoulli pressure, `PressureLaplace` with lamb/convective terms, and `SurfaceVorticityForce` on the same steady solved state. This should distinguish pressure-recovery issues from circulation/velocity-field issues.
- Sweep `correct_kuttacondition` in pressure and force monitors. A material CT shift would point toward TE treatment rather than missing wake physics.
- Run a semi-infinite rigid-wake smoke test using the existing `Das` setup. This is not expected to be final physics, but it can reveal whether a cheap wake representation moves the solution in the right direction.

## Acceptance

This track survives only if at least one check materially changes CT or exposes a force-method/backend disagreement on the same solved state. If all checks agree and CT remains low, demote implementation-error explanations and focus on wake modeling.

## Verification Notes

Prefer the smallest script that reuses the rotor-hover setup and writes comparable CT values for each variant. If this touches code, start with the relevant postprocess/simulate unit tests before any long example run.

## Results — 2026-06-15

Audit run via `examples/rotor_hover_force_method_audit.jl` (dji9443 40_40 mesh, 5400 RPM,
steady rigid wake). Reference band: experiment 0.072, BEM 0.068 (item 003), particle-wake
0.062.

| Check | Result | Verdict |
| --- | --- | --- |
| 1 — backend | `DirectBackend` vs FMM agree to **0.000%** on the same solved state | no backend error |
| 2 — force method | **16.3% spread** on the SAME state | real disagreement, follow-up |
| 3 — Kutta condition | Bernoulli Δ0.4%, surface vorticity Δ0.0% | negligible, not the cause |
| 4 — semi-infinite rigid wake | Bernoulli 0.0505 → **0.110 (+118%)** | materially moves CT |

Check 2 detail — CT by force method on one identical solved state:

| Method | CT |
| --- | --- |
| Bernoulli | 0.0505 |
| `PressureLaplace` | 0.0520 |
| `SurfaceVorticityForce` | 0.0603 |
| Kutta–Joukowski | 0.0506 |

Surface vorticity sits closest to the particle-wake value 0.062.

Check 4 overshoots the 0.068–0.072 reference rather than landing in it.

### Verdict

The track survives on the force-method spread (Check 2), but the dominant lever is clearly
wake modeling (Check 4); backend and Kutta condition are cleared. Focus shifts to wake
physics (items 002/004), consistent with 003 (BEM = 0.068).

### Setup gotcha

Shedding must be derived from the watertight-reoriented `RigidWakeBody.cells`, **not** the raw
mesh cells. Using raw mesh cells attaches the wake at the wrong edges with no error and
collapses CT to ~0.014. (Now documented as a critical invariant in `CLAUDE.md`.)
