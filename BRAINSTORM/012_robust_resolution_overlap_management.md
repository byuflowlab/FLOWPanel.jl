# Robust Resolution / Overlap Management

## Context

Items 008 and 010 (M1/M2/M4/M6) converge on a single root cause: the settled
rotor-hover particle wake sits at **heavy overlap**. Backing the overlap ratio out
of the diagnostics (item 011, `eff_overlap_sample`, measured directly from
nearest-neighbor spacing and *local* core size), the **true local overlap is
`\sigma/h\approx4.2`** — about `4\times` the textbook `\sigma/h\approx1\text{–}1.5` —
with `\approx394` particles within `2\sigma_{\mathrm{local}}` and a kernel-weighted
effective overlap `\approx181` (diag-dominance `\approx0.007`), versus an *order-10*
neighborhood at nominal overlap. (The often-quoted `\sim2\times10^4` neighbor count
is **not** the overlap: it is a `4\cdot\sigma_{\max}` cutoff artifact — `\sigma_{\max}`
is the *global* max core size, `\sigma_{\max}/\sigma_{\min}\approx15`, so that cutoff
`\approx1.7\,R` engulfs half the wake. Use `\sigma/h` and the effective count, not
the raw cutoff count.) Three otherwise separate problems all scale with this overlap:

- **Conditioning** (010 M6): the overlap operator (scalar `M` *and* block `K_\sigma`)
  is near-singular (diagonal-dominance `\approx0.006`–`0.010`), so the M2/M6 overlap
  solve is intractable.
- **Vorticity divergence discrepancy** (010): `\eta=\|\omega_\zeta-\omega_J\|/\|\omega_\zeta\|\approx0.22`,
  large precisely because the field is far from a smoothly-overlapped solenoidal one.
- **Pedrizzetti workload** (010 M12, item 011): the relaxation correction angle is
  *entirely* overlap-induced; with no overlap, relaxation is idle.

010's repeated finding is that **existing machinery (stretching, Pedrizzetti, RBF
core reset) works fine at low overlap** — the failures are all a heavy-overlap
regime the wake *drifted into*, not intrinsic to the method or kernel. This item
asks whether the most leverage comes not from a richer kernel or an expensive solve,
but from **smarter, more robust particle-resolution management that bounds the
overlap ratio while preserving vortical structure**. If `\sigma/h` can be held near
nominal, all three problems above soften at the source.

Relation to the mission (rotor-hover `C_T`): this is potentially the highest-leverage
wake-modeling lever, but also the one most **coupled to the validation target** —
over-aggressive merging smears tip-vortex peak vorticity and changes `C_T`. So it is
both the most promising and the most dangerous direction, and must be evaluated with
`C_T` and peak-vorticity preservation as first-class criteria.

## Proposed Path

1. **Characterize the overlap state (prerequisite for everything).** Measure the
   `\sigma/h` distribution across the wake and over time on the
   `rotor_hover_pressure_comparison` states — local spacing `h` from the CSR neighbor
   grid (already built in `examples/particle_overlap_residual.jl`), `\sigma` from the
   saved state. Map where overlap is worst (tip-vortex roll-up vs diffuse wake) and
   how it grows (core spreading vs particle concentration). This is the common
   independent variable for items 008/010/011/012.

2. **Audit the current handling.** Read `FLOWVPM.jl/src/FLOWVPM_merging.jl` and the
   RBF core-reset / core-spreading in `FLOWVPM_viscous.jl`: what triggers merging
   today (proximity, strength, count)? Does any trigger target the overlap ratio
   `\sigma/h`? Quantify what the run's current merging actually removes
   (`data/rotor_hover_pressure_comparison/particles_merging_savings.csv` exists —
   relate it to the overlap map).

3. **Design an overlap-bounded, moment-conserving policy.** Merge/remesh only where
   `\sigma/h` exceeds a target, conserving the moments that define the local
   structure: 0th (circulation, always) and ideally 1st (centroid / impulse, to keep
   tip vortices in place and un-smeared); consider 2nd (core size) for fidelity. Note
   the conceptual tie: **merging is itself a local overlap refit** — the same
   projection problem as the 008/M2 solve, but small and bounded, so it inherits
   their machinery and caveats.

4. **Validate offline, then live.** Offline: confirm the policy reduces `\sigma/h`
   and the overlap residual / `\eta` / conditioning without collapsing local peak
   vorticity. Live: rerun the rotor-hover case and compare the **`C_T` history and
   tip-vortex peak vorticity** before/after — the policy is only acceptable if it
   bounds overlap *and* preserves the validation observable.

## Acceptance Criteria

- A measured `\sigma/h` distribution over the rotor-hover wake and time, with the
  worst-overlap regions identified.
- An audit of the current merging / core-reset triggers and what they remove.
- A proposed overlap-bounded moment-conserving merging policy, with offline evidence
  that it lowers overlap (and thereby conditioning / `\eta` / relaxation workload)
  **without** damaging tip-vortex peak vorticity.
- A live before/after on `C_T` (and peak vorticity) demonstrating the policy does not
  corrupt the validation target — or a clear statement of the trade-off if it does.

## Caveats

- **Highest physics coupling of the improvement tracks.** Aggressive merging diffuses
  vortical structures, lowers peak tip-vortex vorticity, and changes `C_T` — the exact
  quantity items 003/005/006 are trying to validate. Treat `C_T` and peak-vorticity
  preservation as first-class acceptance, not afterthoughts; collides directly with
  005/006's stable-wake work.
- Bigger project than 008/010/011; stage it **behind** item 011 (is the solve even
  needed?) and 010's M6 consumer audit (does the divergence discrepancy even matter?),
  so effort targets the confirmed lever.
- Conserving only circulation (0th moment) is insufficient to protect structure —
  first-moment/impulse conservation is likely needed to avoid smearing; M12's
  Stokes-circulation-vs-volume-moment distinction is relevant to what "conserve" means.
- Do not assume a free lunch: the well-conditioned/low-`\eta` limit is also the limit
  where the overlap correction shrinks toward self-only (010), so reducing overlap may
  make the solve unnecessary rather than tractable — which is itself a fine outcome,
  but should be stated, not assumed away.
