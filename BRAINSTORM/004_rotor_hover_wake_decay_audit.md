# Rotor Hover CT: Wake-Decay Audit

## Context

Wake-bearing rotor-hover runs can reportedly pass through `CT ≈ 0.072` and later drift toward lower values. Prior audit has already checked geometry, RPM, density, CT convention, radius/diameter interpretation, and handedness/sign, so this track focuses on the time evolution of the wake-bearing solution.

## Hypothesis

The early CT agreement may be a transient, or the later CT decay may be numerical wake modeling. The key is to determine whether particle wake diffusion, merging, SFS, core-size evolution, or body-on-wake coupling is reducing induced loading over time.

## Audit Targets

- Particle wake diffusion settings and their time scale relative to rotor revolutions.
- Particle merging and whether circulation or near-wake structure is being smoothed too aggressively.
- SFS and V-SFS sensitivity.
- Core size growth and its effect on blade-induced velocities.
- Body-on-wake coupling and whether blade influence on the wake remains consistent over long runs.
- CT, circulation, inflow, and radial loading histories around the moment where `CT ≈ 0.072` occurs.

## Diagnostics

- Compare CT histories with diffusion, merging, and SFS variants.
- Track circulation conservation and near-wake strength by wake age.
- Compare early-agreement states against later states using the same force recovery method.
- Probe whether radial loading changes uniformly or decays first in specific blade regions.

## Acceptance

This track is useful if it identifies whether late-time CT decay is numerical wake modeling, a startup transient, or expected physics. If disabling or modifying one wake mechanism preserves CT near the expected neighborhood without introducing obvious artifacts, that mechanism becomes the next focused target.

## Caveats

Do not treat the first crossing of `CT ≈ 0.072` as validation by itself. It must be tied to stable wake and loading behavior over multiple subsequent revolutions.
