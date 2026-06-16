# Rotor Hover CT: CCBlade/BEM Comparison

## Context

Prior audit has already checked geometry, RPM, density, CT convention, radius/diameter interpretation, and handedness/sign. CCBlade/BEM should be used as a diagnostic cross-check, not as proof that either solver is correct by itself.

## Hypothesis

A BEM comparison can separate missing wake/inflow physics from geometry, airfoil, Reynolds-number, or operating-condition assumptions. Starting at nonzero advance ratio should avoid hover-specific singular behavior and make the comparison easier to interpret before approaching hover.

## Proposed Path

- Build a CCBlade case using the same blade geometry, operating point, density, RPM, and reference definitions.
- Start at a small but nonzero advance ratio.
- March toward hover only after the nonzero-advance case has sensible radial trends.
- Compare panel and BEM outputs in dimensional and nondimensional form.

## Quantities To Compare

- Radial loading.
- Inflow angle.
- Induced velocity.
- Circulation.
- Integrated CT.
- Any local section assumptions that differ from the panel model, especially airfoil polar and Reynolds-number treatment.

## Acceptance

This track is useful if it clearly separates missing-wake physics from geometry/airfoil/Re assumptions. For example, if BEM agrees with VPM/experiment while the panel result remains low under equivalent inputs, focus on wake and induced-velocity modeling. If BEM also predicts low CT under the same assumptions, revisit section data and operating-condition interpretation.

## Caveats

BEM and panel methods do not represent the same physics. Use agreement or disagreement to route the investigation, not as a final validation by itself.
