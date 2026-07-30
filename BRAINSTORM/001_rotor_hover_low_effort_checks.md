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
