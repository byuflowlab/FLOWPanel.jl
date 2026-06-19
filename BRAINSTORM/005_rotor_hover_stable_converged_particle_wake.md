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

---

## 2026-06-17 — Implementation + first tuning loop (QUALIFIED RESULT)

### What was built

`examples/rotor_hover_pressure_comparison.jl` gained the staged-startup machinery
this track called for (all ENV-driven, defaults leave prior behavior intact except
the freestream is now a pulse rather than a no-op ramp to `magVinf=0.0001`):

- **Four-phase smoothstep freestream pulse** replacing the old monotonic
  `ramp_magVinf`: ramp-up → hold → withdraw → hover, in rotor revolutions.
  Hooks: `MAGVINF_PEAK`, `MAGVINF_END`, `FREESTREAM_RAMP_REVS`,
  `FREESTREAM_HOLD_REVS`, `FREESTREAM_WITHDRAW_REVS`, `SETTLE_REVS`,
  `MAGVINF_START`. Reuses the existing `smoothstep`; `Uinf(t)` signature unchanged.
- **Truncation depth hook** `TRUNCATION_DEPTH_R` (× R) for the `GlobalCylinder`.
- **Run-length bookkeeping**: total revs = `max(NREVS, ramp+hold+withdraw+settle)`
  so a run never ends mid-withdrawal; phase boundaries printed at startup.
- **`BERNOULLI_ONLY`** monitor set (skips both CG pressure solves + the per-step
  FMM Hessian) for cheaper tuning iterations.
- **Plateau diagnostics** printed at the end: peak-to-peak CT ripple over the
  final `SETTLE_REVS` plus the residual `magVinf` at readout.
- Driver `examples/run_rotor_hover_stable_wake.sh` (withdraw-rate sweep, bash-3.2
  safe), modeled on the 004 driver.

RPM spin-up used the pre-existing `SPINUP_REVS` / `SPINUP_START_FRACTION`.
Freestream sign verified correct: the +x pulse convects the wake downstream (same
axis as the truncation cylinder, opposite the −x thrust); confirmed by particles
being *retained* inside the downstream `GlobalCylinder` (`keep` deletes anything
with axial projection outside `[origin, origin+extrude]`).

### Headline finding — startup ring killed; intrinsic hover oscillation persists

Representative run `data/stable_wake_iter4` (40_40, RPM 6000, NT 15, Bernoulli;
spinup 1.5 rev from 0.4×RPM; pulse peak 5 m/s ≈ v_h, ramp 1 / hold 1 / withdraw
2.5 / settle 3 rev; `TRUNCATION_DEPTH_R=2`). CT(Bernoulli) vs rev:

| phase (rev) | behavior |
| --- | --- |
| spinup+ramp (0–1) | small 0.089 bump at rev 0.2, settles to ~0.047 by rev 1. **No 0.11 overshoot ring** (cf. impulsive baseline). |
| hold, 5 m/s (1–2) | **flat ~0.055**, peak-to-peak ≈ 0.006 |
| withdraw 5→0 (2–4.5) | smooth monotonic rise 0.055 → 0.086 (climb→hover re-loading) |
| settle, magVinf=0 (4.5–8.8) | **non-damping oscillation** 0.086→0.075→0.093→0.082, period ~1.5–2 rev, 2nd peak *higher* than 1st; final-3-rev peak-to-peak **0.0183** |

**Positive:** the staged startup eliminates the impulsive *startup* transient — the
baseline's CT≈0.11 overshoot followed by a ~6-rev ring (0.11→0.067→0.071→0.061→…)
is replaced by a smooth approach to a flat ~0.055 band by rev ~1.3 while the
freestream is held.

**Open question (do NOT overclaim):** withdrawing the freestream did not produce a
flat plateau in this run — a CT oscillation appears in the settle window. An
initial read called it "intrinsic / non-damping / growing," but a closer look at
`stable_wake_iter4` shows a **clean, regular oscillation of period ≈ 3.5 rev**
riding on a **rising mean** (the withdraw relaxation + the wake column still
building) and spanning only **~1.3 cycles** — far too few to decide grow vs steady
vs decay. The run also used aggressive `TRUNCATION_DEPTH_R=2`, which may itself
shape the oscillation. So whether this is an intrinsic wake instability (which
would corroborate 002's gain≈1/no-fixed-point picture and route work to 004's
ablations), a long transient, a truncation artifact, or withdrawal-excited is
**not yet established**.

**Phase 2 plan** (`plans/20260617_brainstorm_005.md`) characterizes it: a long run
for period/amplitude consistency & damping; a truncation-depth control (1/2/4); a
withdraw-gradualness sweep; and then tests the full menu of damping schemes
(viscous core-spreading↑, stronger SFS, kernel overlap/core growth, resolution;
then numerical: relaxation, merging, kernel regularization, body-strength
low-pass) on the best non-oscillating baseline. Driver
`examples/run_rotor_hover_stable_wake.sh` (EXPERIMENT=e1_long|e2_depth|e3_withdraw|e4_damping).

### Caveats on this run

- Magnitude (~0.085 plateau) is **inflated by the aggressive `TRUNCATION_DEPTH_R=2`**
  (deletes far-wake downwash → higher CT); not comparable to the 0.062 baseline.
  The *oscillation* (shape/ripple) is the valid signal, not the level.
- Settle was only 3 rev and partly overlaps the withdraw relaxation; a looser
  truncation (≥4R) + slower withdrawal + longer settle would sharpen whether the
  oscillation amplitude is steady or slowly growing. Not yet run (each full run is
  ~1–1.5 h wall at this fidelity).

### Tuning path (iterations)

- iter1: peak 8 m/s, spinup_start 0.05 → violent negative-CT *windmilling* during
  spinup (strong freestream on a near-static rotor). Rejected: raise
  `SPINUP_START_FRACTION` and/or lower peak.
- iter2/iter3: spinup_start 0.4, peak 5 m/s → clean positive-CT startup; confirmed
  flat ~0.051–0.055 through ramp+hold. Slow 5 m/s convection let particle count
  climb (~140/step, no plateau) → expensive.
- iter4: shorter schedule + `TRUNCATION_DEPTH_R=2` to finish; revealed the
  post-withdrawal oscillation above.

### Next iteration (if pursued)

One confirmation run at `TRUNCATION_DEPTH_R≥4`, `FREESTREAM_WITHDRAW_REVS≥6`,
`SETTLE_REVS≥6` to test (a) whether a slower withdrawal avoids exciting the mode
and (b) the un-truncated plateau magnitude. Expectation from the above: the mode is
intrinsic and re-appears regardless — making 004's ablations the actual lever.
Use `examples/run_rotor_hover_stable_wake.sh` (sweeps the withdraw rate).

---

## 2026-06-19 — Phase 2: characterized the oscillation (RESULT — overcall corrected)

Ran the full Phase 2 campaign (`plans/20260617_brainstorm_005.md`): E1 long run,
E2 truncation-depth gate, E3 withdraw-rate sweep, E4 13-scheme damping menu. New
tooling: `examples/analyze_stable_wake_oscillation.jl` (detrends the rising mean
with a low-order polynomial, peak-picks the residual, reports period consistency +
per-cycle detrended-amplitude trend). E4.8 (`BOUND_STRENGTH_RLX` body bound-
circulation low-pass) implemented in `simulate!` (the only new solver code the plan
required). All runs: 40_40 mesh, RPM 6000, Bernoulli-only, viscous core-spreading +
SFS ON.

**Bottom line: the earlier "intrinsic / non-damping / GROWING" call was wrong on
two of three counts.** The post-withdrawal hover signal is a **bounded, regular,
*steady* limit cycle**, and most of its amplitude in the old `iter4` run was a
**near-disk wake-truncation artifact**, not wake physics.

### E1 — long run (30 settle revs, depth=1): it's a STEADY limit cycle (Q1/Q3)
`data/stable_wake_e1_long_e1_long`. The 1-rev boxcar mean reaches a **flat plateau
~0.0845 by rev ~12** (the rev 0–12 "rise" is a startup/withdraw-reloading hump that
*decays*, not a slow mode). On the genuinely-settled window (rev 12.9–36.9, **48
cycles**): raw ptp 0.0189 ≈ detrended ptp 0.0178 (confirming a flat mean), **period
0.498 ± 0.021 rev** (≈ **2/rev**, i.e. blade passage), **amplitude growth ≈
−0.00003/cycle ⇒ STEADY limit cycle** — *non-damping* (won't vanish by waiting) but
*not growing* (not an instability). The original ~3.5-rev "period" was the decaying
startup hump on a too-short (~1.3-cycle) window, not the settled signal.

### E2 — truncation-depth gate (1/2/4): depth=1 is an ARTIFACT proxy (Q4)
Matched schedule, 12-rev settle window:

| depth | mean CT | detrended ptp | period (rev) |
| --- | --- | --- | --- |
| 1 | 0.0719 | 0.0335 | 0.52 ± 0.03 (clean 2/rev) |
| 2 | 0.0698 | 0.0122 | 0.83 ± 0.82 |
| 4 | 0.0609 | 0.0068 | 2.20 ± 1.59 |

Amplitude shrinks **~5×** from depth 1→4, **period scales with depth** (the plan's
"deletion-front transit time" signature), and mean CT drops toward the more physical
0.061. ⇒ **The large `iter4`/depth=1 ripple is substantially a near-disk wake-
truncation artifact; depth=1 is NOT a trustworthy fast proxy.** All subsequent
characterization was run at **depth ≥ 4R**. Adequate truncation depth is itself the
single biggest, most *physical* "damper."

### E3 — withdraw-rate sweep (depth=4): gradualness is a (modest) lever (Q2)
Final-10-rev window:

| config (depth=4) | mean CT | detrended ptp |
| --- | --- | --- |
| withdraw 2.5 rev | 0.0601 | 0.0091 |
| withdraw 6 rev | 0.0602 | 0.0087 |
| **withdraw 12 rev** | 0.0599 | **0.0042** |
| spinup-only (no pulse) | 0.0849 | 0.0231 |

Slower withdrawal **shrinks** the settle oscillation (~2× from 2.5→12 rev). The
pure-spinup control (no freestream pulse) is **worse** on *both* amplitude (0.023)
and mean (0.085) ⇒ the staged freestream pulse is **beneficial**, it is *not* the
thing that rings. Mean CT is withdraw-rate-invariant (~0.060) ⇒ a robust hover
plateau. **Best non-damping baseline = depth=4 + withdraw_12** (detrended ptp 0.0042).

### E4 — damping menu on the depth=4 + withdraw_12 baseline (the deliverable)
13 schemes, one variable each, final-8-rev window. Baseline detrended ptp **0.00231**
(≈ **3.8 % of CT 0.060** — already near-flat). Ranked by ptp (smaller = more damping);
P = physical, N = numerical, D = diagnostic:

| scheme | type | det ptp | Δ vs base | mean CT | reading |
| --- | --- | --- | --- | --- | --- |
| sfs_off | D | 0.00079 | −66 % | 0.0594 | **UNRELIABLE** (1.9 cyc, slow mode; not a real reduction) |
| corebeta_hi (β 1.5→2) | P | 0.00216 | −6 % | 0.0602 | within noise |
| **baseline** | — | 0.00231 | — | 0.0602 | reference |
| nu_x10 | P | 0.00234 | +1 % | 0.0599 | within noise |
| nu_x3 | P | 0.00244 | +6 % | 0.0601 | within noise |
| merge_aggressive | N | 0.00251 | +9 % | 0.0591 | within noise |
| bound_rlx α=0.5 (E4.8) | N | 0.00262 | +13 % | 0.0598 | no benefit at this baseline |
| pps_hi (p/step 2→4) | P/res | 0.00318 | +38 % | 0.0621 | worse (resolution) |
| relax_off | N | 0.00624 | +170 % | 0.0722 | **worse**; relax-ON stabilizes |
| sfs_strong | P | 0.00741 | +221 % | 0.0603 | **worse** |
| nt_hi (NT 18→36) | P/res | 0.01043 | +352 % | 0.0562 | **worse, period stays clean 0.543 ± 0.029** |
| kernoff_hi | N | 0.01110 | +381 % | 0.0630 | worse |
| overlap_hi (3→4) | P | 0.02844 | +1131 % | 0.0974 | **breaks the wake** (CT way off) |

**Menu conclusions (physical first, as the plan asked):**
1. **No tested damper meaningfully *reduces* the ripple.** The baseline is already so
   flat (~4 % of CT) that the two nominal "reductions" are within measurement noise
   (`corebeta_hi`, −6 %) or an unreliable few-cycle artifact (`sfs_off`). The real
   damping was the **physical configuration (truncation depth ≥4R + slow withdrawal)**
   established in E1–E3 — not any add-on scheme.
2. **The residual 2/rev limit cycle is PHYSICAL, not under-resolution.** Refining
   resolution (`nt_hi` NT=36, `pps_hi`) *sharpens/worsens* it while keeping a clean
   stable period (`nt_hi` 0.543 ± 0.029 rev, the most regular of all 13) — the plan's
   E4.4 diagnostic ⇒ "**damp, don't fix**." But at ~4 % of CT it is small enough to
   live with; there is little to damp.
3. **Several knobs are actively harmful:** `overlap_hi` destroys the wake (mean CT
   0.097), `relax_off` lets CT drift to 0.072 with 2.7× the ripple, `sfs_strong` and
   `kernoff_hi` worsen it. Viscous diffusion (`nu_×3/×10`) and core-β are inert here.
4. **E4.8 `bound_rlx` (α=0.5) gives no benefit** on this baseline (ripple already
   physical + tiny). The body-strength low-pass would only matter on a baseline with a
   genuine feedback-driven oscillation; logged as an honest negative for the last-
   resort numerical knob (code retained behind `BOUND_STRENGTH_RLX`, default off).

### Verdict vs the earlier overcall
- "intrinsic" → **partly true**: a small 2/rev limit cycle is physical and resolution-
  robust, but most of the *observed* `iter4` amplitude was a **depth=2 truncation
  artifact**, removed by going to depth ≥4R.
- "non-damping" → **true** (steady limit cycle; it does not self-decay).
- "growing" → **false** (growth ≈ 0 over 48 cycles).
- Net: at **depth=4 + slow withdrawal** the thrust history *does* settle to a flat
  ~0.060 plateau with only a ~4 %-of-CT physical 2/rev ripple — close to the item's
  "stable, non-oscillating plateau" acceptance, achieved by **physical** levers. The
  remaining ripple is not worth chasing with numerical dampers (all marginal or
  harmful). Magnitude (0.060 vs 0.068–0.072 ref) remains 002/003/004's concern.

### Reproduce
`examples/run_rotor_hover_stable_wake.sh` (`EXPERIMENT=e1_long|e2_depth|e3_withdraw|
e4_damping`); analyze with `examples/analyze_stable_wake_oscillation.jl <run_dir>
[settle_revs] [rpm]`. Data: `data/stable_wake_e{1,2,3,4}_*`. Best baseline:
`TRUNCATION_DEPTH_R=4 FREESTREAM_WITHDRAW_REVS=12 SETTLE_REVS≥10`.
