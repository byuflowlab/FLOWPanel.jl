# 026 — Particle splitting to cap sigma growth (design)

Status: DESIGN + implementation sketch (2026-08-27). No implementation yet.
Origin: the 018 NT144 GPU performance cliff; this item is deliberately
standalone (the mechanism is general FLOWVPM physics, not campaign-specific).
Companion band-aid: an absolute σ cap in the rVPM update (see §8), deployed to
unblock 018 NT144 arms while this item is designed/implemented properly.

## 1. Motivation

Diagnosed 2026-08-27 on the 018 campaign's NT144 GPU arm
(`018_dji9443_hover_convergence_campaign/gpu_nt144_cliff_findings_20260827.md`
has the investigation; conclusion below is self-contained):

- The rVPM area evolution `σ ← σ − Δt·σ·Z` (`FLOWVPM_timeintegration.jl:310-312`)
  grows σ without bound for a particle in sustained compression (Z < 0). The
  052c sigma guard caps only the shrink side (`Δt·Z > cap`) and floors σ —
  growth is unguarded.
- Measured: a **single runaway particle** (σ = 0.0389, 17 % above the #2
  particle at 0.0333) growing at +0.036 %/step while its nearest-size peers
  held steady or shrank. Growth attribution across quantiles: per-step
  $\Delta\sigma^2$ spans ~50× from p90 to max — incompatible with viscous core
  spreading (uniform $\Delta\sigma^2 = 2\nu\Delta t \approx 10^{-8}$/step,
  ~100× too small), incompatible with merging (smooth continuous growth of one
  particle). It is the rVPM compression term.
- Consequence: FLOWVPM's radix-FMM geometry rule sizes the grid from the
  single scalar `sigma_max` (`_radix_sigma_outgrown!`,
  `FLOWVPM_fmm_radix.jl:643`). When sigma_max crossed the cached ell=3
  adequacy limit (computed 0.0388975 at the measured box; crossing landed
  between steps 2252→2253, exactly the observed cliff), the cache rebuilt at
  ell=2 → 64 cells for 267k particles → quasi-dense near field → 6.5→52
  s/step, permanent. One oversized particle degraded the entire solver.
- Viscous core spreading is a secondary, slower driver: with the RBF reset
  disabled on the GPU path (zeta is CPU-only; `Critical:1e9`), the *bulk*
  σ distribution drifts up secularly and pushes the population toward any
  fixed threshold over a 30-rev run.

Splitting is the physical fix: a σ-grown particle is an under-resolved one;
splitting restores resolution instead of distorting the dynamics (as clamping
Z or σ would).

Relation to `018.../sigma_blowup_mechanism.md`: same σ-equation, opposite
regime (that doc: Z > 0 collapse + Γ ignition at Δt·Z > 2). The existing
splitting triggers target the shrink side; nothing here may interfere with
them — growth triggers act on σ above a cap, shrink triggers on σ/σ₀ below a
ratio.

## 2. Existing infrastructure (extension, not new build)

`FLOWVPM.jl/src/FLOWVPM_splitting.jl` (2655 lines) already provides:

- `SplitOptions`, trigger types (`SigmaShrinkTrigger` :303-319, `HoldTrigger`
  :328-344, `All/AnyTrigger` :352-394, an `H_chi` trigger ~:271-294),
  `compute_split_direction` with `STRENGTH`/`STREAMLINE`/`STRAIN1` axes
  (:407-426), `split_particles!` (:517-575, maxparticles guard at :562),
  `_do_split!` (:580-). Default split: symmetric 2-child, Γ/2, offset
  ±a·e_split with 2a = κ_split·σ, `preserve_sigma=true` (σ-constant — built
  for the collapse problem, NOT σ-reducing).
- `SplittingState{R}` per-particle side-buffer (`sigma_0`, `H_chi`,
  hold/cooldown counters), kept in lockstep by `add_particle` /
  `remove_particle` swap-with-last (`FLOWVPM_particlefield.jl:367-384,
  :672-676, :710-713`). `accumulate_H_chi!` (splitting.jl:186-198), called
  per accepted step from `nextstep` (particlefield.jl:769), is the template
  for any new per-particle time-averaged quantity.
- Merging conservation rule to invert (`FLOWVPM_merging.jl:112-184`): Γ vector
  summed exactly, vol summed, σ_new = cbrt(Σσ³), X weight-averaged.
- Design heritage: `FLOWVPM.jl/sfs_musings.md` §"Particle subdivision"
  (~:940-1050) — moment constraints, σ-ratio trigger, hold/cooldown.

Gaps: (a) no growth-side trigger; (b) no σ-reducing child geometry; (c) the
stretching vector S=(MM1,MM2,MM3) and scalar Z=MM4 are stack locals discarded
each substep — no persistent (averaged) stretching axis exists; (d)
`split_particles!` and `merge_particles!` are CPU-only (bare scalar indexing;
`add_particle` has a `_add_particle_broadcast!` GPU path but splitting/merging
do not).

## 3. Two mechanisms, treated separately

Growth attribution is separable per particle: the viscous contribution is
analytically known ($2\nu\,t$; RK3 keeps the accumulator in `M[7]`,
`FLOWVPM_viscous.jl:197-198`), so rVPM contribution = total − viscous. Route
each split event to the mechanism that dominates that particle's growth.

### 3a. Mechanism A — viscous-spreading split (isotropic)

Spreading is isotropic ($\sigma^2 \mathrel{+}= 2\nu\Delta t$ for every
particle), so the split is isotropic-on-average:

- **Geometry**: m = 4 children at the vertices of a regular tetrahedron,
  random orientation per event (avoids lattice bias), centroid at the parent
  position, each child Γ_p/4 parallel to the parent Γ. Total Γ and linear
  impulse exact.
- **Sizing**: conserve volume as the inverse of merging's rule:
  $\sigma_c = \sigma_p \cdot 4^{-1/3} \approx 0.63\,\sigma_p$. Vertex radius
  from second-moment matching $\sigma_p^2 \approx \sigma_c^2 + a^2$:
  $a \approx 0.78\,\sigma_p$, giving $a/\sigma_c \approx 1.24$ — children
  remain overlapped and the combined Gaussian support approximates the
  parent's.

### 3b. Mechanism B — rVPM-compression split (anisotropic; Ryan 2026-08-27)

Kinematics: stretching (Z > 0) thins the core while the tube lengthens;
compression (Z < 0) fattens the core while the tube **shortens**. The
fattened element is a bundle of thinner parallel filaments, so the split
re-discretizes the cross-section:

- **Geometry**: 3 children on the plane whose normal is the **time-averaged
  stretching axis** (≈ Γ̂ for a coherent tube), equilateral triangle centered
  on the parent, random in-plane orientation, each child Γ_p/3 **parallel to
  the parent Γ** (i.e. normal to the plane). Total Γ and centroid exact by
  construction. Mass conservation in the sense of Alvarez's σ-equation
  derivation (tube mass under compression): E. J. Alvarez (2022),
  *Reformulated Vortex Particle Method and Meshless Large Eddy Simulation of
  Multirotor Aircraft*, PhD dissertation, BYU — full reference and σ-equation
  re-derivation in `020_sigma_aware_subgrid_closure/phase_01_theory.md`
  (:971-973 and §theory).
- **Direction tracking (new state)**: running/exponentially-averaged unit
  vector of the stretching axis, added to `SplittingState`, updated by an
  `accumulate_*!` hook from `nextstep` and mirrored in add/remove exactly as
  `H_chi` is. First cut: the existing `STRENGTH` axis (Γ̂) with zero new
  state; the averaged axis is the refinement (compare in validation).
- **Sizing vs overlap trade-off** (decide at implementation; work the
  algebra then):
  - mass-per-length rule $3\sigma_c^2 = \sigma_p^2$ →
    $\sigma_c = \sigma_p/\sqrt{3} \approx 0.58\,\sigma_p$; exact transverse
    second-moment matching then forces triangle radius
    $a \approx 1.15\,\sigma_p$ → inter-child spacing ≈ 3.5 σ_c
    (**under-overlapped**, lumpy cross-section);
  - volume rule $\Sigma\sigma_c^3 = \sigma_p^3$ →
    $\sigma_c \approx 0.69\,\sigma_p$; moment matching gives
    $a \approx 1.02\,\sigma_p$ → spacing ≈ 2.6 σ_c (still loose vs the 018
    campaign's overlap ≈ 2.75 convention).

  Resolution: treat a as free — shrink below the moment-matching value to
  meet an overlap floor (spacing/σ_c ≲ 1.5–2), accepting a slightly narrower
  combined support (small, local, diffuses). The (i)/(ii) σ-rule choice is
  open for Ryan.
- Expected behavior note: the child trio co-rotates slowly under mutual
  induction like a physical vortex bundle — not an artifact.

## 4. Triggers

- New `SigmaGrowthTrigger`: σ > σ_cap (absolute), with hold/cooldown counters
  reusing the `HoldTrigger` pattern. Optionally a ratio form σ/σ₀ > C using
  `SplittingState.sigma_0`.
- Mechanism routing per §3 preamble (viscous share from `M[7]`/t_sgm vs
  total growth).
- Children start with fresh `sigma_0` and a cooldown, and must be protected
  from immediate re-merge (see §6).

## 5. Threshold vs FMM adequacy geometry

At the measured 018 NT144 box (L ≈ 0.512): ell=3 admissibility limit 0.0389,
ell=4 limit 0.0194 (limits scale with box size L, which grows as the wake
extends — cap conservatively). Two candidate operating points:

| σ_cap | FMM depth held | expected effect |
|---|---|---|
| ≈ 0.030 | ell=3 (pre-cliff state) | restores 6–8 s/step; ~handful of splits/run |
| ≈ 0.018 | ell=4 admissible | 8× fewer bodies/cell — possibly net FASTER than pre-cliff; splits ~0.5 % of particles |

Both are worth campaign experiment arms once splitting exists.

## 6. System concerns

- **Merge/split ping-pong**: children carry cooldown counters (exists);
  additionally ensure merge skips fresh children or child spacing exceeds the
  merge radius (`merge_r` policy wired at `FLOWPanel_wake.jl:1720`,
  `MERGE_R_FACTOR` at `rotor_hover_pressure_comparison.jl:105`).
- **Capacity**: maxparticles headroom guard exists (splitting.jl:562); at
  σ_cap ≈ 0.030 the split count is negligible.
- **GPU**: splitting is CPU-only today. At σ-cap rarity, a host round-trip
  during wake maintenance (mirroring however merging is invoked for CuArray
  fields — verify that path) is acceptable; a broadcast port is only needed
  if the σ_cap ≈ 0.018 regime is adopted.
- **Defaults**: both mechanisms OFF by default; per-case env knobs in the
  driver following the `MERGE_R_FACTOR` pattern.

## 7. Verification plan

1. Unit: a single-particle split conserves Γ and impulse exactly and matches
   the parent's induced far-field velocity to tolerance; child overlap and
   combined-support residual within spec for both mechanisms.
2. A/B: rerun the 018 NT144 λ=2.4 arm from a pre-cliff snapshot with
   σ_cap = 0.030 — the cliff must not occur, CT̄ within campaign arm-to-arm
   scatter (±0.36 %).
3. Optional: σ_cap = 0.018 / ell=4 arm measuring s/step against the 6–8
   s/step baseline.

## 8. Interim band-aid (deployed separately, 2026-08-27)

Until splitting lands: an absolute σ ceiling in the rVPM σ update (clamping
`new_sig` from above, alongside the existing 052c floor), env-switchable and
off by default. Physically crude — it freezes the runaway particle's σ
instead of re-resolving it — but at the 018 operating point it touches ~8
particles and keeps `sigma_max` inside the ell=3 adequacy limit, restoring
FMM performance. Remove when this item is implemented.
