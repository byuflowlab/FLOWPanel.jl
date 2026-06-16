# Rotor Hover: Stable (Non-Oscillating) Thrust with Particle Wake

## Goal

**The primary deliverable is a stable — non-oscillating, or nearly so — predicted
thrust history.** Hitting `CT ≈ 0.072` is *not* required for this item to be a
win; a CT that settles to a flat (or only weakly rippling) plateau at *any* value
is the success condition. This track exists to kill the oscillations seen in the
particle-wake CT history, not to chase the reference magnitude. Magnitude
agreement is the concern of 002/003/004; here, oscillation amplitude is the metric.

## Context

Audits 001–004 converge on wake/induced-velocity modeling as the dominant lever
for the rotor-hover CT shortfall (panel ≈ 0.0506 vs experiment ≈ 0.072). Item 004
specifically flags that wake-bearing runs pass through `CT ≈ 0.072` and then
drift — i.e. the thrust history is unsteady/oscillatory rather than settled. This
track is the constructive counterpart to 004: rather than ablating individual wake
mechanisms, assemble a startup/forcing schedule that lets the particle wake
develop smoothly into a self-sustaining hover state so the thrust stops swinging.

## Hypothesis

The thrust oscillation is a symptom of an undeveloped / impulsively-started wake.
Starting impulsively at full RPM in a quiescent (near-zero freestream) domain
slams circulation onto the blades before any wake column exists and leaves the
far wake unconvected, exciting a startup transient that rings rather than damps.
A staged, *smoothly ramped* startup — gradual RPM ramp plus a temporary convecting
freestream that is then withdrawn slowly enough for wake self-induction to take
over — should establish a developed wake column with no sharp excitation, so the
thrust converges to a steady (non-oscillating) plateau and stays there.

## Approach

All required features already exist in
`examples/rotor_hover_pressure_comparison.jl`; this track is about combining them
into one schedule and tuning the timing, not new implementation.

- **Gradual RPM spin-up.** Use the existing `SPINUP_REVS` /
  `SPINUP_START_FRACTION` smoothstep schedule (`spinup_fraction(t)`,
  `set_frame_omega!`, applied in `maneuver!`) to ramp from near-zero to full RPM
  rather than starting impulsively.
- **Ramped convecting freestream, then withdrawn.** The example already ramps
  `magVinf` up via `ramp_magVinf(t)` / `Uinf(t)` over `ramp_nrev`. Extend this to
  a ramp-up-then-ramp-down profile: start at (near-)zero, raise the axial
  freestream high enough to convect the shed wake clear of the rotor and build a
  wake column, then lower it slowly enough that wake-induced (self-induced)
  velocity can sustain the developed structure as the artificial convection is
  removed. The end state should be hover (freestream → near-zero) with a
  self-sustaining wake.
- **Ramp rates must be chosen intelligently — and will require iteration.** The
  RPM and freestream ramp rates/durations are not known a priori; they are the
  central tuning variables of this track and should be set from physical scales
  (e.g. ramp durations measured in rotor revolutions and in wake-convection
  transit times across the rotor disk), then refined empirically. Expect several
  iterations: read the resulting thrust history, see whether the transient damps
  or rings, and adjust ramp duration/shape accordingly. A too-fast ramp re-excites
  the very oscillation we are trying to remove; a too-slow ramp wastes compute.
  Treat finding the rates as an iterative loop, not a one-shot guess.
- **Downstream particle truncation for cost.** Use the existing
  `ParticleMaintenance` / `GlobalCylinder` mechanism (currently
  `GlobalCylinder([-0.5R,0,0], [cylinder_depth,0,0], 1.5R]`, `cylinder_depth = 4R`)
  to truncate particles beyond a chosen number of diameters downstream, keeping
  particle count bounded over long runs. Pick the truncation distance far enough
  that it does not perturb the near-disk induced velocity that sets CT.
- Keep merging (`MergeParticles`), viscous core spreading (`CoreSpreading`), and
  SFS at their current defaults initially; vary only after a stable baseline
  exists (this is where 004's ablations plug in).

## Diagnostics

- **CT-history oscillation amplitude** (Bernoulli, Laplace ∇u, Laplace λ) vs
  revolution — the primary metric. Quantify peak-to-peak ripple over the final
  revolutions after the freestream is fully withdrawn; the goal is to drive this
  toward zero, regardless of the plateau value.
- Bound-circulation and radial-loading histories (`BoundCirculationMonitor`) to
  confirm loading is steady, not drifting or oscillating region-by-region.
- Particle count vs time, to confirm the downstream truncation bounds cost
  without starving the near wake.
- Wake column structure in VTK (developed, coherent helical column at the end
  state) to confirm self-sustainment after freestream withdrawal.

## Acceptance

This track succeeds if there is a startup/forcing schedule under which the
particle-wake **thrust history settles to a stable, non-oscillating (or nearly
non-oscillating) plateau** over multiple revolutions after the convecting
freestream is removed, with bounded particle count and a coherent self-sustaining
wake. **The plateau need not be near 0.072 — a flat history at any value is the
win.** The converged CT value and its distance from the 0.068–0.072 reference is
reported but is secondary (that magnitude question belongs to 002/003/004). If no
schedule produces a sustained, non-oscillating plateau, that itself is an
informative negative tied back to 004's wake-mechanism ablations.

## Caveats

- The freestream withdrawal rate is the delicate knob: too fast and the wake
  column disperses before self-induction can hold it; too slow and the run is
  expensive and the "hover" CT is contaminated by residual freestream. Sweep it.
- A converged plateau is only meaningful at hover (freestream ≈ 0); report the
  residual `magVinf` at the readout so a plateau under nonzero convection is not
  mistaken for hover.
- Downstream truncation distance must be validated as non-perturbing to the
  near-disk inflow before its CT is trusted.
