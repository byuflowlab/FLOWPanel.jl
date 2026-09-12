# 026 — Particle splitting to cap sigma growth (design)

Status: Phase 0 COMPLETE 2026-09-03 (§11 — W1 blocker fixed, W2 accumulators,
W3 merge guard, W4 snapshots located, W5 σ-rule decided, W6 census run; code
on `026-phase0` branches). Design 2026-08-27, revised 2026-09-02 after review
(corrected §3a algebra, M[7] attribution replaced, Phase-0 blockers and
phasing added in §10). Split mechanisms A/B still unimplemented by design.
Origin: the 018 NT144 GPU performance cliff; this item is deliberately
standalone (the mechanism is general FLOWVPM physics, not campaign-specific).
Companion band-aid: an absolute σ cap in the rVPM update (see §8), deployed to
unblock 018 NT144 arms while this item is designed/implemented properly.

## RESET BRIEF

- **What**: particle splitting for FLOWVPM — growth-side triggers (σ cap) +
  two child geometries (isotropic viscous §3a, anisotropic rVPM §3b), plus
  validation of the *existing* shrink-side split on the gpu40/LineGauss
  stretching ignitions (§9). No code written; design only.
- **Why (Ryan 2026-09-02)**: beyond curing the 018 NT144 FMM cliff, the hope
  of the rVPM anisotropic splitting is to conserve overlap of vortical
  structures and give vorticity another dof to deform under high strain →
  improved method stability. Test it paired with a σ floor (floor-only vs
  split-only vs floor+split, §3b purpose + §9.2); the floor only clamps σ
  from below and never gates splitting — trigger must fire on strain
  exposure / on-floor state, not realized σ/σ₀.
- **Status 2026-09-03**: **Phase 0 COMPLETE** (§11). W1 blocker FIXED
  (SplittingState persists in checkpoints; legacy checkpoints reconstruct
  armed-at-ratio-1), W2 Δσ² accumulators implemented + tested, W3 merge
  reconciliation guard landed, W4 snapshots all located (archived tarballs),
  W5 decided rule (i) mass-per-length, W6 census run (all-rVPM on
  p022lg_hr10). Also fixed: `accumulate_H_chi!` was dead in production
  (FLOWPanel bypasses `nextstep`). Branches `026-phase0` in both repos.
- **Status 2026-09-02**: design revised after review. §3a algebra corrected
  (a ≈ 1.35σ_p, child geometry OPEN pending kernel-fit study); M[7]
  attribution replaced by a to-be-built persistent Δσ² accumulator; merge
  does not reconcile split state.
- ~~**BLOCKER W1**~~ FIXED 2026-09-03 (§11): warm start neither restored nor
  initialized `SplittingState` (`FLOWPanel_warmstart.jl:269ff` bypassed
  `add_particle`), and `SigmaShrinkTrigger` refuses `sigma_0 ≤ 0`
  (`FLOWVPM_splitting.jl:309`) — no §9 arm could trigger until fixed.
- **Standing rulings**: no-split discriminator arms run FIRST (020's
  exponential/local-substep update, then `control_no_backscatter_projection`;
  never a positive-Cd floor — clip fires exactly when SFS would amplify).
  Both mechanisms OFF by default; §8 band-aid stays until this lands.
- **Next actions**: Phase 1 no-split discriminator arms (warm-start the
  gpu40/LineGauss ignitions from the archived tarballs in §11-W4: 020 stable
  integrator, safe SFS projection, combination) → Phase 2 shrink-split gate
  with the σ-floor pairing matrix → W7 averaged stretching axis (§10
  Phase-3 prerequisite: sign-invariant EMA director in `SplittingState`,
  persisted via the W1 schema, new `STRETCH_AVG` direction) → Mechanism B →
  Mechanism A (gated; note the W6 local census found zero viscous-dominated
  crossings, §11-W6).
- **Entry points**: §10 phasing; §9 warm-start continuation gate; forensics
  in `018_.../gpu_nt144_cliff_findings_20260827.md` and
  `FastMultipole/MATRIX_OPERATOR_REFACTOR/052-handoff-prompt-2026-08-31v.md`.

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
do not); (e) **warm start does not restore or initialize `SplittingState`**
(`FLOWPanel_warmstart.jl` bulk-writes `pf.particles` directly, bypassing
`add_particle`), so `sigma_0 = 0` on every continued particle and
`SigmaShrinkTrigger` refuses them all (`FLOWVPM_splitting.jl:309`) — the §9
gate cannot fire until this is fixed; (f) merging changes the
representative's σ without reconciling `sigma_0`/hold/cooldown/provenance,
and does not consult split cooldowns at all.

## 3. Two mechanisms, treated separately

Growth attribution is separable per particle: the viscous contribution is
analytically known ($2\nu\,\Delta t$ per accepted step), so rVPM contribution
= total − viscous. Route each split event to the mechanism that dominates
that particle's growth.

**Attribution caveat (2026-09-02): `M[7]` cannot carry this.** Only the RK3
branch of CoreSpreading writes it (`FLOWVPM_viscous.jl:197-198`); `M[7:9]` is
aliased as RBF target-vorticity scratch (`FLOWVPM_viscous.jl:239, 461`); and
under Euler production it accumulates nothing. Attribution needs a new
persistent per-particle pair of accepted-step $\Delta\sigma^2$ accumulators
(viscous, rVPM) in `SplittingState`, with defined split/merge/restart
semantics — see §10 Phase 0.

### 3a. Mechanism A — viscous-spreading split (isotropic)

Spreading is isotropic ($\sigma^2 \mathrel{+}= 2\nu\Delta t$ for every
particle), so the split is isotropic-on-average:

- **Geometry**: m = 4 children at the vertices of a regular tetrahedron,
  pseudo-random orientation per event (avoids lattice bias; drawn
  reproducibly from lineage state, see §4), centroid at the parent
  position, each child Γ_p/4 parallel to the parent Γ. Total Γ and linear
  impulse exact.
- **Sizing (corrected 2026-09-02; child geometry OPEN)**: conserve volume as
  the inverse of merging's rule:
  $\sigma_c = \sigma_p \cdot 4^{-1/3} \approx 0.63\,\sigma_p$. Second-moment
  matching must project the vertex radius per axis (tetrahedron vertices give
  $\overline{vv^T} = (a^2/3)I$):

  $$\sigma_p^2 = \sigma_c^2 + a^2/3 \;\Rightarrow\;
  a = \sqrt{3\left(1-4^{-2/3}\right)}\,\sigma_p \approx 1.35\,\sigma_p,$$

  i.e. nearest-child spacing ≈ 3.5 σ_c — **under-overlapped**, the same
  trade-off as Mechanism B. (An earlier draft matched
  $\sigma_p^2 = \sigma_c^2 + a^2$, dropping the /3 projection, and wrongly
  concluded $a \approx 0.78\,\sigma_p$ with children "remaining overlapped".)
  Shrinking `a` to an overlap floor is NOT a small concession: at spacing
  2σ_c the effective per-axis width is ≈ 0.77 σ_p, at 1.5σ_c ≈ 0.71 σ_p —
  a 23–29 % support contraction. The child geometry is therefore left
  unresolved pending a pointwise/$L_2$ kernel-fit and induced-field
  optimization over child count and placement (larger sets, e.g. 6–8 on a
  shell, on the table).

### 3b. Mechanism B — rVPM-compression split (anisotropic; Ryan 2026-08-27)

> **AMENDMENT (Ryan 2026-09-07, third ruling — two-regime stretch
> mechanism):** the stretch mechanism is sign-dependent, superseding this
> section's routing of *both* regimes to the 3-child triangle. Negative
> stretching (compression: tube shortens/fattens) keeps the triangle
> geometry below (`_split_compress_tri3!`). Positive stretching
> (elongation: tube lengthens/thins — the shrink-side events) now routes
> to a **2-child in-line kernel** (`_split_elongate_pair2!`): children at
> ±b along the averaged stretch axis (b = 0.5 σ_p default, spacing 1.0
> σ_p), each Γ_p/2 ∥ Γ_p, and **σ_c = σ_p** — the cross-section is fine,
> the split re-discretizes the LENGTH. Total Γ, centroid, linear impulse
> exact; angular impulse exact by symmetry. Same mechanism/enable flag
> (`enable_stretch_split`); the two kernels are dispatched by event side
> (grow → compression → tri3, shrink → elongation → pair2). Implemented in
> FLOWVPM `src/FLOWVPM_resolution_split.jl` (Phase 2 commit 2).

**Purpose (Ryan 2026-09-02)**: the hope of the rVPM anisotropic splitting is
to (1) **conserve overlap of vortical structures** and (2) **give vorticity
another degree of freedom to deform under high strain**, hopefully improving
method stability. Stability is a first-class success metric alongside
re-resolution — score the §9 arms on it. This also motivates **pairing
splitting with a σ floor** (052c-style). Roles are complementary, not
alternatives (Ryan 2026-09-02): the floor merely keeps σ from shrinking
below a set value — it must **never prevent or gate splitting** — while
splitting supplies the spatial dof. Design consequence: the trigger must not
be starved by the clamp (a σ/σ₀-ratio trigger stops advancing once σ sits on
the floor), so the shrink/stretch trigger should fire on strain exposure
(accumulated ΔtZ) or on a sitting-on-the-floor state, not only on realized
σ. Test floor-only vs split-only vs floor+split (§9).

Kinematics: stretching (Z > 0) thins the core while the tube lengthens;
compression (Z < 0) fattens the core while the tube **shortens**. The
fattened element is a bundle of thinner parallel filaments, so the split
re-discretizes the cross-section:

- **Geometry**: 3 children on the plane whose normal is the **time-averaged
  stretching axis** (≈ Γ̂ for a coherent tube), equilateral triangle centered
  on the parent, pseudo-random in-plane orientation (reproducible from
  lineage state, see §4), each child Γ_p/3 **parallel to
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
  meet an overlap floor (spacing/σ_c ≲ 1.5–2), accepting a narrower combined
  support (quantify per event; the §3a numbers show this can reach tens of
  percent). The (i)/(ii) σ-rule choice is open for Ryan, and is **gated on
  the Phase-0 `vol`-consumer audit (§10)**: rule (i) drops σ³-implied volume
  42 %/event and breaks merge-inversion (σ = cbrt(Σσ³)). Note SFS evolution
  uses σ³, not `vol`; the actual dynamic `vol` consumers are the
  PSE/core-spreading RBF reconstruction, merging, and restart I/O — the audit
  must decide whether children carry `vol_p/m` (conserving Σvol, breaking
  vol∝σ³) or σ³-consistent vol.
- Expected behavior note: the child trio co-rotates slowly under mutual
  induction like a physical vortex bundle — not an artifact.

## 4. Triggers

- New `SigmaGrowthTrigger`: σ > σ_cap (absolute), with hold/cooldown counters
  reusing the `HoldTrigger` pattern. Optionally a ratio form σ/σ₀ > C using
  `SplittingState.sigma_0`.
- Mechanism routing per §3 preamble, using the new persistent Δσ²
  accumulators (NOT `M[7]` — see the §3 attribution caveat and §10 Phase 0).
- Children start with fresh `sigma_0` and a cooldown, and must be protected
  from immediate re-merge (see §6).
- **Reproducibility / lineage state (new)**: split orientations are
  pseudo-random but must reproduce across warm-start continuations. Hashing
  particle index/step/position is fragile (swap-with-last moves indices,
  positions are FP-sensitive, restarts lose an RNG sequence). Instead add to
  `SplittingState`: a persistent particle-lineage ID, split generation, split
  event ordinal, and the run seed; orientation = f(seed, lineage ID,
  ordinal). All of it persists across warm starts (§10 Phase 0, blocker
  W1).

## 5. Threshold vs FMM adequacy geometry

At the measured 018 NT144 box (L ≈ 0.512): ell=3 admissibility limit 0.0389,
ell=4 limit 0.0194. Limits scale with box size L, which grows as the wake
extends — so a fixed σ_cap *gains* margin over time at fixed depth (a growing
box raises the limit). The limit can still move against the cap when depth,
near-field radius, bounds, or cache state change, so: **log
`sigma_max / current_adequacy_limit` at every radix rebuild and warn above
~0.9**. Two candidate operating points:

| σ_cap | FMM depth held | expected effect |
|---|---|---|
| ≈ 0.030 | ell=3 (pre-cliff state) | restores 6–8 s/step; ~handful of splits/run |
| ≈ 0.018 | ell=4 admissible | 8× fewer bodies/cell — possibly net FASTER than pre-cliff; splits ~0.5 % of particles |

Both are worth campaign experiment arms once splitting exists.

## 6. System concerns

- **Merge/split ping-pong AND state reconciliation**: children carry
  cooldown counters (exists), but the merge routine does not consult them —
  it can recombine fresh children immediately — and it changes the
  representative's σ without reconciling `sigma_0`, exposure/hold/cooldown
  counters, or provenance. The design must specify merge callbacks/state
  rules (what merging does to every `SplittingState` field) and the
  maintenance functional-policy ordering, plus a spacing guard vs the merge
  radius (`merge_r` policy wired at `FLOWPanel_wake.jl:1720`,
  `MERGE_R_FACTOR` at `rotor_hover_pressure_comparison.jl:105`).
- **Capacity**: maxparticles headroom guard exists (splitting.jl:562); at
  σ_cap ≈ 0.030 the split count is negligible. Log candidates *skipped* by
  `max_fraction`/`maxparticles`, not only successful splits.
- **GPU**: splitting is CPU-only today, but the seam already exists —
  `_apply_particle_maintenance_device!` does D2H → host maintenance → H2D
  with explicit side-buffer sync (`FLOWPanel_gpu_wake.jl:109-123`). The
  requirement is a correctness/performance test that `split_particles!`
  works through that seam (not an existence audit); a broadcast port is only
  needed if the σ_cap ≈ 0.018 regime is adopted.
- **Defaults**: both mechanisms OFF by default; per-case env knobs in the
  driver following the `MERGE_R_FACTOR` pattern.

## 7. Verification plan

1. Unit: a single-particle split conserves Γ and linear impulse exactly;
   matches the parent's induced far field to tolerance in **velocity AND
   velocity gradient/strain and the SFS inputs** (ignition is strain-driven —
   velocity alone is insufficient); child overlap and combined-support
   residual within spec for both mechanisms; **quantify angular-impulse
   error** (linear impulse is exact by symmetry, the transverse child
   geometries generally do not conserve angular impulse); measure child
   mutual-induction timescales vs Δt (splitting must not introduce new
   stiffness).
2. Collective: a single vortex ring, force-split every particle once per
   mechanism; verify ring translation speed, circulation, impulse, and
   enstrophy drift vs the no-split control over ~1 convective time. Catches
   collective geometry errors that far-field unit checks miss.
3. A/B: rerun the 018 NT144 λ=2.4 arm from a pre-cliff snapshot with
   σ_cap = 0.030 — the cliff must not occur, CT̄ within campaign arm-to-arm
   scatter (±0.36 %), **s/step back in the 6–8 band**, split/skip telemetry
   recorded.
4. Optional: σ_cap = 0.018 / ell=4 arm measuring s/step against the 6–8
   s/step baseline.

## 8. Interim band-aid (deployed separately, 2026-08-27)

Until splitting lands: an absolute σ ceiling in the rVPM σ update (clamping
`new_sig` from above, alongside the existing 052c floor), env-switchable and
off by default. Physically crude — it freezes the runaway particle's σ
instead of re-resolving it — but at the 018 operating point it touches ~8
particles and keeps `sigma_max` inside the ell=3 adequacy limit, restoring
FMM performance. Remove when this item is implemented.

## 9. 2026-08-31 gpu40 ignition: shrink-side splitting validation case

The completed 40-revolution GPU run `scr_p019_s038v_gpu40` supplies a second,
more diagnostic motivation for splitting than the sigma-growth/FMM cliff that
originated this item. Particle-level forensics found a smooth, localized
stretching runaway rather than an FMM or GPU error:

- Patient zero was particle index 102340 (stable in the saved VTP ordering),
  in the aged outer wake at approximately 1.5R off-axis and 2R downstream.
  Its strength grew from |Gamma| = 4.3e-5 at step 850 to 4.0e-3 at step 993,
  then 0.39 at step 998. Over the same interval sigma contracted from
  1.17e-3 m to 4.9e-4 m, while the persistent local velocity-gradient norm
  rose from roughly 500--1000 1/s to 9e3 1/s.
- The saved Gamma increments agree with the implemented stretching update and
  a constant dt of approximately 2.3--2.4e-4 s. Direct O(N) regularized
  Biot--Savart summation reproduced the saved velocity at patient zero to
  3.2e-4 relative error. The GPU/FMM far field is therefore exonerated at the
  ignition site; the later FMM adequacy failure was a symptom of the already
  corrupted wake.
- A partly antiparallel partner (index 179085) remained about 4.1e-3 m away,
  approximately 8 sigma near ignition, and ignited one step later. At 8 sigma
  the Gaussian particle kernel is already effectively in its singular
  far-core regime, so sigma contraction did not simply "uncover" this
  partner. The supported positive feedback is instead growing Gamma ->
  stronger pair/ambient strain at nearly fixed geometry -> further Gamma
  growth. Sigma contraction is the resolution-loss signal and increases the
  danger of any closer interactions.
- This has a direct physical interpretation for splitting: a material vortex
  element under extension becomes longer, thinner, and generally curved. Its
  increasing vectorial strength must be distributed along that increasing
  length. Keeping it as one isotropic shrinking blob concentrates the growing
  moment at one point after the element has lost overlap with its neighbors.
  Splitting along the strength/stretching direction supplies the missing
  spatial degrees of freedom; merely flooring sigma does not.

### SFS and other damping available in the failed run

The run used `DynamicSFS(Estr_fmm, pseudo3level_beforeUJ,
pseudo3level_positive_afterUJ)` with alpha=0.999, Lagrangian relaxation
`rlxf=0.005`, 0 <= Cd <= 1, `clipping_backscatter`, and no magnitude or
directional controls. Although the implementation is named pseudo-three-level,
alpha=0.999 is the code's effective two-level configuration.

Cd did not converge at patient zero. It repeatedly switched between its upper
bound and zero: approximately 0.995 at steps 850--855, zero at 856--860,
0.990 at 985, zero at 986--987, 0.952 at 988, zero at 989, 0.996/0.987 at
990/991, zero at 992--993, 1.0 at 994, 0.436 at 995, and zero at 996--998.
The backscatter clip prevents an antidissipative SFS contribution by setting
Cd=0, but supplies no replacement forward-scatter dissipation. Thus the only
continuous term capable of directly damping an individual particle's
|Gamma| was absent through much of the ignition window.

Other nominal stabilizers did not provide a Gamma sink: CoreSpreading changed
sigma but not Gamma; `WAKE_CORE_BETA=1e9` made its RBF strength-redistribution
reset unreachable; corrected-Pedrizzetti relaxation (`RELAX_RLXF=0.3`)
preserved |Gamma| while changing direction; and the current absolute merge
radius was about 0.62 mm, far below the 4.1 mm patient/partner separation.
There was no active split/remesh operation. Spatial trimming acted only after
the runaway particles left the retained wake.

### Proposed warm-start continuation gate for shrink-side splitting

(Wording note: this is a *warm-start continuation* — the repo's "replay" mode
runs monitors without evolving the wake. Prerequisite: warm start must
restore/reconstruct `SplittingState`, §10 blocker W1, or no trigger fires.)

**No-split discriminator arms run FIRST** — before any split-geometry work,
continue the same pre-ignition window with, separately and combined:

- (a) the pointwise-stable local σ,Γ update from item 020 (§"contractivity
  observation", `020_sigma_aware_subgrid_closure.md:~252-262`): ignition is
  forward-Euler overshoot at ΔtZ > 2/3; the exponential update
  σ ← σe^(−ΔtZ), Γ ← Γe^(−3ΔtZ) (or local sub-stepping) removes the ΔtZ
  ceiling with zero new physics. Item 020 owns this closure/stability
  question (028 is the separate variable-filter-width consistency question).
- (b) `control_no_backscatter_projection`
  (`FLOWVPM_subfilterscale.jl:536`), which removes only the amplifying SFS
  component. Do NOT instead force a positive Cd floor: `clipping_backscatter`
  fires exactly when the SFS term would amplify |Γ|, so re-enabling any
  positive Cd there amplifies.

If (a) or (b) alone arrests ignition, the lever is integrator/SFS repair and
splitting's role narrows to resolution maintenance — that outcome redirects
this item before geometry work starts.

Then use the ignition as an end-to-end validation case once the split
mechanism is wired through `PanelParticleWake`:

1. Warm-start the unmodified case from a retained pre-ignition state
   (prefer step 950 for sufficient lead time, with step 985 as the short test)
   and confirm that the no-split control reproduces patient-zero growth and
   ignition near steps 995--998.
2. Enable shrink/stretch-triggered two-child splitting along `STRENGTH` first;
   compare `STRAIN1` when a persistent strain-axis history is available. The
   trigger should fire before local overlap is lost, not after the explicit
   update reaches dt*Z = O(1). Per the §3b purpose statement, also run the
   **σ-floor pairing matrix**: floor-only, split-only, and floor+split. The
   floor only clamps σ from below and never gates splitting; in the
   floor+split arm verify the trigger still fires for particles sitting on
   the floor (strain-exposure or on-floor trigger, not realized σ/σ₀ —
   §3b). Score all arms on stability (bounded max|Gamma|/sigma^2, no
   compounding leader) as well as re-resolution.
3. Require exact parent-to-child Gamma-vector and impulse conservation,
   bounded far-field velocity error at the split, restored local overlap, and
   no immediate merge/split ping-pong.
4. Dynamic acceptance: carry the replay through at least step 1070 with
   bounded max|Gamma|/sigma^2 and max|u|, no fixed compounding leader, and no
   FMM sigma-adequacy fallback. Pre-trigger rotor loads and the unsplit bulk
   wake should remain within the replay/control numerical tolerance.
5. Record split count, locations, child ancestry, Cd/clipping state, dt*Z,
   nearest-neighbor distance/sigma, and the resolved versus SFS contributions
   to d|Gamma|/dt. These diagnostics distinguish successful re-resolution
   from a split rule that merely delays ignition.

The same experiment can be repeated on the independent LineGauss run, whose
root-region patient zero ignited around steps 490--516. Passing both events
would be substantially stronger evidence than curing one kernel-specific
trajectory. The forensic assets and provenance are summarized in the task-052
handoff `FastMultipole/MATRIX_OPERATOR_REFACTOR/052-handoff-prompt-2026-08-31v.md`;
the local gpu40 particle windows cover steps 850--998 and the full ignition
window 985--1010, while the LineGauss window covers steps 450--520.

## 10. Phase structure (added 2026-09-02 review)

**Phase 0 — audits and blockers (no split geometry until these close):**

- **W1 (BLOCKER)**: warm-start `SplittingState` semantics — restore or
  defensibly reconstruct `sigma_0`, exposure/hold/cooldown counters,
  provenance, and lineage IDs on warm start (§2 gap (e);
  `FLOWPanel_warmstart.jl:269ff`, `FLOWVPM_splitting.jl:307-319`).
- **W2 (BLOCKER for §4 routing)**: replace `M[7]` attribution with persistent
  accepted-step viscous/rVPM Δσ² accumulators (§3 caveat). First pin the
  production integrator and per-step `M` clearing behavior, then spec the
  accumulator's split/merge/restart semantics.
- **W3**: merge/split state reconciliation rules + policy ordering (§6).
- **W4**: snapshot availability — confirm the §7.3 pre-cliff and §9
  step-950/985 restart states still exist locally or in `/nobackup/archive`
  tarballs (VTK retention keeps only the newest 5 steps).
- **W5**: `vol`-consumer audit (§3b) → the (i)/(ii) σ-rule decision.
- **W6**: census of σ_cap crossings by attribution (viscous- vs
  rVPM-dominated) on an existing mature run — gates Mechanism A (see below).

**Phase 1** — no-split warm-start continuations of the gpu40 (and LineGauss)
ignitions: 020 stable local integrator, safe SFS projection, and their
combination (§9 discriminator arms).

**Phase 2** — validate the already-existing two-child longitudinal shrink
split through the §9 gate, with the full §7.1 diagnostics.

**Phase 3** — Mechanism B, only if growth events remain relevant after
Phases 1–2; §7.2 ring test + §7.3 018 NT144 A/B.

- **W7 (Phase-3 prerequisite, added 2026-09-03)**: persistent time-averaged
  stretching axis for the §3b split direction — closes §2 gap (c) (S and Z
  are stack locals discarded each substep; the only axes
  `compute_split_direction` offers today are *instantaneous* Γ/U/J). Build
  it as the W2 pattern extended to a vector:
  - `stretch_axis` (3 × maxparticles) in `SplittingState`, updated as an
    exponential moving average of the applied stretching direction at the
    same integrator sites that accumulate `dsigma2_rvpm` (euler CPU loop has
    S = MM1–MM3 in hand; euler_exp has `L*G`).
  - **Sign-invariant director averaging**, not a raw-S EMA: raw S cancels
    under oscillating strain, and the §3b geometry needs only the plane
    normal to the axis. Align each sample to the running mean before
    accumulating (EMA of ±S with the sign that gives positive dot with the
    current mean), or a per-particle structure tensor if that proves noisy.
  - Same lockstep/reset semantics as W2 (add/remove swap-with-last; zero on
    split children, merge representative, and RBF reset) and persistence
    through the W1 VTP schema (one more optional `split_*` array) — an EMA
    is exactly the state a warm start would otherwise silently zero.
  - New `STRETCH_AVG` member of `SplitDirection` reading it, falling back to
    `STRAIN1` while the EMA is unconverged (e.g. the first ~1/α steps after
    birth/split, tracked by comparing particle age against the EMA
    timescale or by a norm threshold on the accumulated director).
  - Open: the EMA timescale α (relate to the §4 trigger's exposure window so
    the axis averages over the same history that armed the trigger).

**Phase 4** — Mechanism A, only if the W6 census shows viscous-dominated cap
crossings AND the §3a kernel-fit study lands an acceptable child geometry.

## 11. Phase 0 results (2026-09-03)

All six W-tasks closed. Code on branches `026-phase0` of FLOWPanel
(`fastmultipole` base) and FLOWVPM (`flowpanel` base); worktree session under
`~/Dropbox/research/projects/worktrees/026/`.

**Design-doc amendments discovered during implementation:**

- `SplittingState` has exactly four fields (`sigma_0`, `H_chi`,
  `hold_counter`, `cooldown_counter`) — the "provenance and lineage IDs"
  assumed by §2/§4 **do not exist** in code. They are also unnecessary for
  reproducible split orientations: `compute_split_direction`
  (`FLOWVPM_splitting.jl`) derives the axis from Γ/U/J, all of which the
  checkpoint already persists. `FilamentEdgeGraph` persistence remains a
  non-goal (warm-started runs begin with an empty edge graph; affects
  filament-edge arms only).
- Checkpoints are VTK (`.vtp` per step) + TOML with **no version field**;
  field-presence probing is the established back-compat idiom and is what W1
  uses.
- **`accumulate_H_chi!` was dead in production**: FLOWPanel drives the wake
  via `_euler`/`_euler_exp` directly (`FLOWPanel_wake.jl`), bypassing
  `FLOWVPM.nextstep` which hosts the H_chi hook — the same class of bug as
  W1. Fixed: called explicitly after the integrator step.

**W1 (blocker) — FIXED.** Checkpoints now carry six optional `split_*`
point-data arrays (four state vectors + the two W2 accumulators; counters as
Int32, reals at series precision). Loader: all six present → exact restore;
none present (legacy campaign checkpoints, incl. fp64-sidecar era) →
reconstruction `sigma_0 :=` restored post-SIGMA_CEIL σ (armed, ratio exactly
1, one-shot warn); partial set → typed error. `sigma_0` is never clamped
(creation-time reference; clamping would spuriously arm overgrown
particles). Fallback semantics: shrink accrued pre-restart will not
re-trigger — under-triggers only; faithful trigger continuity requires
post-change checkpoints. Tests: exact round-trip incl. trigger-decision
equality and a split firing on a warm-started field (the blocker
regression); legacy fallback incl. SIGMA_CEIL cross-check; IGE warm-start
suite 113/113. (Pre-existing, unrelated: testsets using the immutable
`WarmstartNoopSolver` error on Julia 1.12 via a `WeakKeyDict` finalizer
issue — fails identically on the untouched checkout.)

**W2 — implemented.** `SplittingState` gains `dsigma2_visc`/`dsigma2_rvpm`
(Δσ², accumulated as the *applied* post-guard delta at each σ-update site:
euler CPU rVPM update, euler_exp geometric contraction, RK3 per-stage, and
the three CoreSpreading branches). Cleared on split children, merge
representative, and CoreSpreading RBF reset (numerical re-projection, not
physics). Device-backed broadcast paths skip accumulation (splitting is
CPU-only; follow-up documented in-code if GPU splitting lands). Persisted
with W1. Invariant tested per integrator × viscous scheme:
$$\sigma^2(t) - \sigma^2(t_0) = \Delta\sigma^2_\mathrm{visc} + \Delta\sigma^2_\mathrm{rVPM}$$
(44 assertions, `FLOWVPM/test/runtests_dsigma2_accumulators.jl`).

**W3 — ruling + guard.** Confirmed gap: `_finalize_merged_particle!` never
touched splitting state — the representative kept a stale `sigma_0` while σ
jumped to `cbrt(Σσ³)`. Guard landed: on merge the representative's splitting
state is reset wholesale (`sigma_0 :=` merged σ; exposure/counters/
accumulators zeroed) — a merged particle is a new entity. Removed members
are handled by `remove_particle`'s existing swap-with-last lockstep.
**Policy ordering ruling**: maintenance order inside a step is the
`ParticleMaintenance` chain order (functional policies run in tuple order,
then trim); when both merge and split policies are active, list
`MergeParticles` before `SplitParticles` so a fresh child cannot be
immediately re-absorbed by its sibling, and rely on `N_cooldown` for the
converse (child re-splitting).

**W4 — snapshots all located (gate PASS).**

| snapshot | where | steps |
|---|---|---|
| gpu40 ~950 | `/nobackup/archive/usr/rander39/FLOWPanel_runs/projects_FLOWPanel.jl/scr_p019_s038v_gpu40__052-h200.tar.zst` | 0–1062 |
| gpu40 ~985 | `…/scr_p019_s038v_gpu40__018-gpu-gh200.tar.zst` | 985–1475 |
| 018 NT144 pre-cliff (cliff @2253) | `…/p018_csarc_n5_nt144_l2p4_s2gpu.tar.zst` | 1333–2278 (use 2200–2248) |

No `.bson` restarts exist; all warm-start material is `.vtp` (consistent
with W1). NT144 pre-cliff states must be extracted from the tarball at
Phase-1 launch.

**W5 — decision: rule (i), mass-per-length, σ_c ≈ 0.58 σ_p.** Full
`vol`-consumer audit: `vol` feeds no vortex dynamics (zero reads in UJ/FMM
kernels, SFS, relaxation, integrators). Only consumers: RBF CG *initial
guess* (convergence, not the converged answer; tolerates vol=0), PSE
*overwrites* vol from σ every step, merge `vol_sum` is pass-through
bookkeeping, rest is I/O. Production FLOWPanel already sheds particles with
vol=0. Volume conservation buys nothing physical; take the
overlap-friendlier rule. Hygiene: any future split writes children's
vol = (4/3)πσ_c³ (matches the edge-refinement convention).

**W6 — census script + first result.**
`scripts/p026_sigma_cap_census.jl` (FLOWPanel): offline pass over a particle
VTP series; NN position-matching across adjacent steps; classifies σ_cap
crossings by realized Δσ² vs the exact viscous budget 2ν·dt. Smoke run,
p022lg_hr10 steps 300–329 (σ_cap = 7.5e-3 m ≈ p99.9, ν = 1.48e-5):
1741 crossings, **all rVPM-dominated, zero viscous-dominated** — p99
Δσ²/step ≈ 800× the viscous budget, max ≈ 24,600×. Known limitation: a
fast-moving runaway particle racks up false crossings via NN mismatch
(observed: the top-10 table is one runaway matched against successive small
neighbors), inflating *counts*; the dominance *split* is robust since
nothing approaches the viscous budget. Local evidence therefore keeps
**Mechanism A gated off**; run the census on a gpu40/018 series before any
final ruling.

## 12. Phase 1 launch record (2026-09-03)

**§11-W4 tarball mapping CORRECTION (found during launch).** The Vatistas
gpu40 and its LineGauss twin share the SAME run name
`scr_p019_s038v_gpu40` in different silos, and the W4 mapping confused
them. Verified from run-residue wake-health monitors:

| archive tarball (`/nobackup/archive/.../projects_FLOWPanel.jl/`) | actual identity | evidence |
|---|---|---|
| `scr_p019_s038v_gpu40__052-h200.tar.zst` (steps 0–1062) | **LineGauss twin** (052-h200 silo, job 13518861 line) | ignition at ~500–530 (max_u 66→2.2e5); post-ignition corpse by 950 (max_u ≈3e3) |
| `scr_p019_s038v_gpu40.crashed1061.todelete.tar.zst` (0–1061) | **Vatistas gpu40, part 1** (018-gpu-gh200 silo) | healthy at 940–990 (max_u ≈17–24), ignition 995–1000 (54.6→213, min_σ collapse) |
| `scr_p019_s038v_gpu40__018-gpu-gh200.tar.zst` (1060–1475) | Vatistas gpu40, part 2 (post-crash chain) | starts at 1060 already post-ignition (max_u 2.4e4) |

A first launch (jobs 13569074–77) warm-started from step 950 of the
`__052-h200` tarball believing it the Vatistas run — i.e. from the LG
post-ignition corpse — and was **scancel'd ~40 min in**; its outputs
(`data/scr_p026ph1_*` dirs + monitor CSVs) were deleted. Silver lining:
the mis-mapping's discovery located the LG pre-ignition state, so the §9
LineGauss repeat launched immediately instead of pending a snapshot hunt.

**Corrected launch — eight arms, both ignition events** (m12
`--qos=normal`, ~29 s/step CPU, ≈1–1.5 h each):

| job | case | event | arm |
|---|---|---|---|
| 13569125 | `scr_p026ph1_ctrl_gpu40` | Vatistas, restart 950 → 1100 | control |
| 13569127 | `scr_p026ph1_exp_gpu40` | " | (a) `WAKE_EXPINT` |
| 13569129 | `scr_p026ph1_proj_gpu40` | " | (b) `SFS_NO_BACKSCATTER_PROJECT` |
| 13569131 | `scr_p026ph1_expproj_gpu40` | " | (a+b) |
| 13569126 | `scr_p026ph1_ctrl_lg` | LineGauss, restart 450 → 600 | control |
| 13569128 | `scr_p026ph1_exp_lg` | " | (a) |
| 13569130 | `scr_p026ph1_proj_lg` | " | (b) |
| 13569132 | `scr_p026ph1_expproj_lg` | " | (a+b) |

- Cases appended to the `run_p018_screen_hpc.slurm.sh` case table (cluster
  + local copies in sync; cluster backup `.bak_p026`). All clone the
  `scr_p019_s038v_gpu40` env (OVERLAP 2.4, PPS 11, MERGE_R 0.00524, N=1,
  DAS_UNIFORM 3.4, viscous CoreSpreading β=1e9, rlxf 0.3) +
  `WAKE_HEALTH_DTZ=true` + `WAKE_HEALTH_ATTRIBUTION=true`. The `_lg`
  cases add `FLOWPANEL_FILAMENT_REG=linegauss` (verified from the LG
  twin's own job banner: linegauss reg, otherwise identical env).
- Warm-start sources extracted from the archive tarballs (archived run
  dirs untouched — avoids ARCHIVED-STALE):
  `data/p026_restart_gpu40_s950/` (Vatistas step 950, from crashed1061
  tarball) and `data/p026_restart_lg_s450/` (LG step 450, from
  `__052-h200`). Submission env: `RESTART_STEP={950|450}
  RESTART_NAME=scr_p019_s038v_gpu40 RESTART_PATH=data/p026_restart_*`.
  Neither tarball holds fp64 particle VTPs; the loader uses fp32 `.vtp`.
- Run lengths: gpu40 `NREVS=29.5833` → total step 1100 (spinup_steps=35 +
  round(36·29.5833)), past the §9 gate step 1070; LG `NREVS=15.6944` →
  total step 600 (ignition window 490–516 + ~85 steps margin).
- **Backend deviation from the source run (deliberate)**: all four arms run
  the CPU host-array pfield (`VPM_ARRAYTYPE=array` driver default) on the
  cluster's `unified-052` stack, because `_euler_exp` has no GPU/broadcast
  path (`FLOWVPM_timeintegration.jl:429`, `Threads.@threads` scalar loop).
  Backend is therefore matched across arms (each differs by exactly one
  knob) and the ctrl arm doubles as the §9-step-1 check that ignition
  (steps ~995–998) reproduces off-device. No 026 code deployed to the
  cluster — Phase 1 needs none (splitting off; `euler_exp`, the projection
  control, and cross-tag warm-start all pre-exist in unified-052); the W2
  Δσ² attribution columns are consequently absent from these arms' logs.
- Provenance for the LG event (from the 052 handoff
  `FastMultipole/MATRIX_OPERATOR_REFACTOR/052-handoff-prompt-2026-08-31v.md`):
  LG patient zero = particle idx 160932, root-vortex region ~0.13R
  off-axis, |Γ| 4e-4 (450) → 2.2e-2 (490) → 37 (517). It is NOT
  `p022lg_hr10` (that is the 022 ground-effect corpse, ignition step 710).
- Scoring next: §9 criteria per event — does each ctrl reproduce its
  patient-zero growth (gpu40: idx 102340, ignition ~995–998; LG: idx
  160932, ~490–516); do (a)/(b)/(a+b) arrest it (bounded max|Γ|/σ², no
  compounding leader) — from
  `monitors/scr_p026ph1_*_monitor04_wake_health*.csv` (max_dtZ,
  attribution columns) and force monitors. Passing BOTH events is the
  strong-evidence outcome §9 names. If (a) or (b) alone arrests ignition,
  the item redirects to integrator/SFS repair before any split geometry.

## 13. Phase 1 harvest and §9 ruling (2026-09-03)

All eight arms ran to full length (gpu40 → step 1100, LG → step 600; no
crashes; completion judged from logs + wake-health CSV coverage, not sacct).
Scored from `monitors/scr_p026ph1_*_monitor04_wake_health_system1.csv` and
`*_CT_vs_rev.csv` (`CT_bernoulli`; the force monitor's CFz column is NOT CT).

**End-window summary** (max over run; min for σ; "u>100" = first step
max_u exceeds 100):

| case | max_u | max γ/σ² | min σ | max dtZ | u>100 |
|---|---|---|---|---|---|
| ctrl_gpu40 | 1.3e6 | 4.1e10 | 9.41e-5 | 3.3e3 | 1011 |
| proj_gpu40 | 3.2e6 | 2.2e10 | 9.41e-5 | 4.7e4 | 1005 |
| exp_gpu40 | 48 | 2.7e4 | 2.8e-4 | 0.314 | never |
| expproj_gpu40 | 38 | 2.8e4 | 3.3e-4 | 0.225 | never |
| ctrl_lg | 2.1e4 | 2.9e7 | 9.41e-5 | 2.7e2 | 504 |
| proj_lg | 3.5e5 | 1.2e9 | 9.41e-5 | 9.7e3 | 487 |
| exp_lg | 35 | 5.2e3 | 3.9e-4 | 0.262 | never |
| expproj_lg | 37 | 3.8e3 | 4.0e-4 | 0.136 | never |

**Findings:**

- **Ctrl reproduction (§9 step 1): PASS, with a timing note.** LG ctrl
  ignites in the reference window (u>100 at 504; ref 490–530). gpu40 ctrl
  ignites ~10–15 steps late (u>100 at 1011 vs ref 995–1000) but with the
  identical signature: min_σ collapse to the 9.41e-5 floor, γ/σ² → 1e10,
  dtZ crossing 2/3 before blowup (1.29 at step 1040). Attributed to
  warm-start perturbation (fp32 VTP restore) + CPU-vs-GPU backend shifting
  chaotic timing; ignition itself is backend-independent.
- **(a) `WAKE_EXPINT` alone ARRESTS BOTH ignitions.** Both exp arms stay
  bounded to end-of-run with max_dtZ ≤ 0.31 < 2/3, no min_σ collapse, no
  compounding leader. (a+b) likewise. Consistent with the §9 mechanism:
  ignition is forward-Euler ΔtZ overshoot, removed by the exponential
  σ/Γ update with zero new physics.
- **(b) SFS no-backscatter projection alone FAILS both events**, igniting
  marginally *earlier* than ctrl (1005 vs 1011; 487 vs 504 — within
  chaotic-timing noise, but clearly no arrest). Backscatter removal is not
  the lever.
- **Pre-trigger loads agree.** Pre-ignition-window mean CT: gpu40 arms
  0.0720–0.0723 (<0.4% spread), LG arms 0.0744–0.0746 (<0.3% spread).

**§9 binding ruling: the redirect branch FIRES.** (a) alone arrests both
independent ignitions → the lever is integrator repair (item 020 owns the
closure/stability question; the integrator's definitive write-up is
`020_sigma_aware_subgrid_closure/phase_02r_integrator.md` — corrected
frozen-gradient map + exact CoreSpreading composition); splitting's role
narrows to resolution maintenance. Phase 2 split-geometry work is ON HOLD pending Ryan's ruling.
Run dirs `data/scr_p026ph1_*` remain on orc unarchived.

## 14. Phase 1b staging (2026-09-04, Ryan-approved ordering)

Follow-on to the §13 ruling, staged in this order:

**Task 1 — GPU/broadcast path for `euler_exp`** (currently CPU-only:
`FLOWVPM_timeintegration.jl:429` is a `Threads.@threads` scalar loop).

- Two sites need device counterparts: (i) the `_euler_exp` frozen-gradient
  update (position, `exp(dt·L)·Γ` + `r^(-3g)`/`r^(-g)` rescale, M[9］Zeff
  stash, SFS Lie split); (ii) the CoreSpreading euler_exp branch
  (`FLOWVPM_viscous.jl:175-194`, the `expm1` blended diffusion using M[9]).
- Follow the existing split pattern: `pfield.particles isa Array` → scalar
  loop, else broadcast (`update_particle_states_broadcast_reformulated!` /
  `_corespreading_euler_broadcast!` are the templates — row-slice views +
  preallocated scratch rows).
- The 3×3 matrix exponential is the only broadcast-unfriendly piece.
  Options: closed-form 3×3 exp via Cayley–Hamilton/Putzer in elementwise
  broadcast form, or a custom CUDA/KernelAbstractions kernel with
  per-thread StaticArrays (likely simpler and matches the CPU code 1:1).
- Keep existing invariants: device paths skip the `dsigma2_*` splitting
  accumulators (precedent + comment at `FLOWVPM_timeintegration.jl:383`);
  `sigma_guard`/`SIGMA_CEIL` stays incompatible with euler_exp; the
  non-finite-ratio DomainError guard must survive on device (or be checked
  post-hoc).
- Tests: CPU-vs-GPU parity on a small random field (positions, Γ, σ to
  rtol ~1e-6 fp32); the 020 Phase-2R suite still green; GPU smoke per
  `agent_policies/HPC.md` (m13h or mgh) confirming `FLOWVPMCUDAExt`.

**Task 2 — RK3 discriminator arms** (prediction ON RECORD: RK3 raises the
linear-stability ceiling only to ΔtZ ≈ 2.51/3 ≈ 0.84 vs forward Euler's
2/3, while the §13 runaway ramps ΔtZ ~0.3 → 1.3 in ~20 steps → expect
DELAY, NOT ARREST. Either outcome is valuable: arrest ⇒ GPU-ready fix
exists today; failure ⇒ falsifies order-of-accuracy as the lever and
justifies Task 1 for production).

- **Wiring is NOT an env knob**: FLOWPanel's `step!` calls
  `_euler`/`_euler_exp` directly with U/J pre-evaluated by its own
  panel+wake orchestration (`FLOWPanel_wake.jl:2226-2237`); RK3 needs
  UJ+SFS re-evaluated at each of 3 stages. Use FLOWVPM's
  `rungekutta3(pfield, dt; custom_UJ)` hook with a closure that re-runs
  the wake's particle UJ evaluation per stage. Design decision to make
  explicitly: freeze the panel-surface influence within the step
  (recommended for the discriminator — panels don't move mid-step and the
  solve happens once) and re-evaluate particle–particle+panel-on-particle
  velocities/gradients per stage; document whatever is chosen. Add a
  `WAKE_INTEGRATOR=rk3` (or similar) knob; keep `WAKE_EXPINT` semantics
  unchanged. Note `viscousdiffusion` already has the RK3 per-stage
  CoreSpreading branch (`aux1/aux2`), so β=1e9 CoreSpreading composes.
- Two arms, CPU host-array, backend-matched to §12:
  `scr_p026ph1_rk3_gpu40` (restart 950 → 1100) and `scr_p026ph1_rk3_lg`
  (restart 450 → 600, `FLOWPANEL_FILAMENT_REG=linegauss`), same restarts
  (`data/p026_restart_gpu40_s950`, `data/p026_restart_lg_s450`,
  `RESTART_NAME=scr_p019_s038v_gpu40`), same env clone + wake-health
  knobs, `-p m12 --qos=normal`. Cost ≈ 3× UJ ⇒ ~85–90 s/step, ~4 h/arm.
- Score identically to §13 (max_u / γ/σ² / min_σ / max_dtZ trajectories,
  u>100 step, pre-trigger CT); compare ignition delay vs ctrl.
- **Deployment caution**: this one DOES require new code on the cluster
  (unlike §12). The `~/projects` trees serve live 018/022 campaigns — use
  a `git worktree` on orc (per HPC.md) with the RK3 branch and point the
  launcher's `*_REPO_OVERRIDE`/`*_PROJECT_OVERRIDE` at it; do not mutate
  the live checkouts.

Still pending Ryan: archiving the eight §12 run dirs (hpc-storage);
notebook entry for Phase 0 + Phase 1 (deferred once, "not yet").

## 15. Phase 1b RK3 discriminator harvest and ruling (2026-09-04)

Arms per §14 Task 2: `scr_p026ph1_rk3_gpu40` (warm start 950, target 1100;
orc job 13582167) and `scr_p026ph1_rk3_lg` (warm start 450, target 600;
job 13582168), run from the `~/wt026` worktrees (FLOWPanel `c9b411d`,
`WAKE_INTEGRATOR=rk3`). **Neither reached its target: both ignited and
died with `ERROR: PARTICLE OVERFLOW` (500k cap) during blowup-driven
shedding** — gpu40 at step 983 (n_particles 178k→227k in one step before
the cap), LG at step 525. The ignition event is fully captured in both
wake-health CSVs, so the truncation does not affect scorability (arrest
is falsified by the ignition itself). Scored from
`monitors/scr_p026ph1_rk3_*_monitor04_wake_health_system1.csv`.

**End-window summary** (§13 format; ctrl rows repeated for reference):

| case | max_u | max γ/σ² | min σ | max dtZ | u>100 | dtZ>2/3 | end |
|---|---|---|---|---|---|---|---|
| rk3_gpu40 | 3.1e6 | 6.3e10 | 4.70e-5 | 9.8e3 | **959** | 957 | overflow @983 |
| rk3_lg | 2.8e6 | 4.1e9 | 4.71e-5 | 2.6e4 | **517** | 497 | overflow @525 |
| ctrl_gpu40 (§13) | 1.3e6 | 4.1e10 | 9.41e-5 | 3.3e3 | 1011 | — | ran to 1100 |
| ctrl_lg (§13) | 2.1e4 | 2.9e7 | 9.41e-5 | 2.7e2 | 504 | — | ran to 600 |

**Findings:**

- **Ruling: NO ARREST — and no systematic delay either.** Ignition
  timing vs ctrl: gpu40 −52 steps (959 vs 1011), LG +13 steps (517 vs
  504). The signature is identical to §13 ctrl in both: min_σ collapse,
  γ/σ² → 1e9–1e10, compounding leader, dtZ through the ceiling. RK3's
  linear-stability ceiling (ΔtZ ≈ 0.84) was crossed at steps 957/497 and
  saved neither arm. This is the §14 prediction's core claim
  (order-of-accuracy is not the lever) confirmed *more strongly than
  predicted* — the expected DELAY did not even materialize.
- **gpu40's early ignition is chaotic-timing scatter, not a wiring
  fault.** rk3_gpu40 departs ctrl essentially immediately after
  warm-start (dtZ 0.38 by step 958 vs ctrl's 0.07 at the same step;
  u>100 nine steps in), while rk3_lg tracked its ctrl twin closely for
  ~65 steps (max_u, γ/σ², dtZ all matching to within noise through step
  ~515) before igniting slightly *later* than ctrl. The LG arm is
  therefore the wiring control: if the RK3 stage re-evaluation were
  broken, it would not reproduce ctrl's trajectory for 65 steps. The
  gpu40 restart state sits near-critical (§13: the runaway ramps
  ΔtZ 0.3→1.3 in ~20 steps; ctrl itself ignited 10–15 steps late vs its
  own reference), so a different integrator's truncation error picks a
  different chaotic realization; timing scatter of ±tens of steps at
  gpu40 is expected. Net: (−52, +13) brackets zero → no delay signal.
- **min σ reached 4.70e-5, below the 9.41e-5 floor seen in every §13
  arm** — the §13 "floor" is evidently the forward-Euler practical bound,
  not a hard clamp; RK3's multi-stage positions let the collapse run
  deeper before death.
- **Pre-ignition loads: CT_bernoulli is UNAVAILABLE for both rk3 arms.**
  The driver writes `*_CT_vs_rev.csv` only after `simulate!` returns
  (`examples/rotor_hover_pressure_comparison.jl:1468`), and both arms
  crashed inside `simulate!`. Fallback check via `monitor02_force` CFz
  window means (same-step windows, rk3 vs ctrl): gpu40 steps 951–958
  5.75e-4 vs 5.45e-4 (+5.6%, acceptable given the 8-step window); LG
  steps 451–500 −1.6e-6 vs +9.8e-5 — both ≈0 (CFz oscillates about zero
  at this normalization), so the LG window mean is an inconclusive load
  check, not a discrepancy. No CT-level pre-trigger comparison is
  possible for this pair.
- **Cost: RK3 measured ~5–10× the euler arms, well above the 3×-UJ
  estimate.** gpu40 ~300 s/step and LG ~150 s/step pre-ignition (vs
  ~29 s/step for the §12 euler arms; §12 predicted ~85–90 s/step). The
  per-stage re-evaluations rebuild FMM trees and panel-on-particle
  influence each stage; post-ignition steps ballooned to ~900 s.

**Binding ruling: the §13 redirect stands and is strengthened.** RK3 is
falsified as a production path (no arrest, no delay, 5–10× cost); the
lever remains the exponential integrator (`WAKE_EXPINT`/euler_exp), now
GPU-capable per Task 1. Splitting's role stays narrowed to resolution
maintenance. Logs:
`~/wt026/FLOWPanel.jl/logs/slurm/slurm-fp-p026ph1-rk3{,-lg}-1358216{7,8}.{out,err}`;
run dirs `data/scr_p026ph1_rk3_{gpu40,lg}` on orc (unarchived).

## 16. Expint-fails hunt: shortlist (2026-09-04, per Ryan's step-4 ask)

Goal: events where the exponential integrator does NOT arrest blowup, to
map the limits of integrator repair and re-motivate Phase 2 splitting.
Sources: 020 evidence pack + a login-node sweep of every
`monitor04_wake_health` CSV under `~/projects/FLOWPanel.jl/data`
(ignition = first max_u>100; "pre-crossing dtZ" = max over rows strictly
before that step).

**Candidate 1 (STRONGEST — the discriminator pair already exists, no new
runs needed).** The 020 Phase-2R σ/R=0.02 viscous screen:
`scr_p019_s020v` (euler ctrl) ignites @213 with pre-crossing dtZ 0.46;
`scr_p020r_geom_s020v` (job 13154223, the CORRECTED frozen-gradient
geometric map — same family as production euler_exp) stays healthy past
the ctrl's death, then a tail-localized contraction ignites @242
(pre-crossing dtZ 0.47 < 2/3; min σ/σ_shed 0.0137, M=222, 19.6 km/s,
non-finite-ratio guard stop @243); the rerun `scr_p020r_geom_s020v_rr`
ignites @210 (pre-crossing dtZ 0.39) — *earlier than ctrl*. Adjudicated
in `020_.../phase_02_evidence_pack.md` ("stiff integration is a real
local numerical remedy, but not a sufficient field-level remedy"): the
regime is under-resolved (Leg-1 pre-onset M≈28), so field-coupled Γ
amplification runs away regardless of integrator — a resolution-loss
failure, the original 026 splitting motivation. Residue: monitor CSVs on
/home (`data/scr_p020_exp_s020v`, `data/scr_p020r_geom_s020v{,_rr}`,
`data/scr_p019_s020v`) + committed figure CSVs
(`020_.../figures/fig_stage_b{,_r}/`). NO VTP checkpoints and no archive
tarballs → any warm-started variant would need a cold rerun (~250 steps,
cheap). Superseded-prototype caveat does not apply to the 2R runs.

**Candidate 2 (categorical, no run needed).** σ-growth/SIGMA_CEIL cliff
events (018 NT144 cliff class): `sigma_guard`/`SIGMA_CEIL` is
*incompatible* with euler_exp by design, so growth-side cliffs sit
outside expint protection by construction. Cite, don't test.

**Candidate 3 (weak until instrumented).** `p022lg_hr10`: ignition
@643 (transient; terminal collapse ~710 with min σ → 9.2e-5), but its
wake-health CSV predates the max_dtZ column, so the not-overshoot claim
can't be scored from residue. Would need offline dtZ reconstruction or a
re-run with current monitors — only worth it if Candidate 1 is rejected.

**Sweep caveats:** ~33 older runs report dtZ≡0 (column absent/unpopulated)
and were excluded from dtZ-based ranking. Low pre-crossing dtZ alone is a
WEAK discriminant — even `scr_p026ph1_ctrl_gpu40` (which expint arrests)
shows pre-crossing dtZ 0.34, because dtZ crosses 2/3 only mid-runaway;
the load-bearing evidence for "expint fails" is the direct 2R pair, not
the dtZ census.

**§16 reruns LAUNCHED (Ryan-approved 2026-09-05):** cold reruns of the
pair on the current stack, from the `~/wt026` worktree, m12/normal, full
VTP series retained for future warm starts. Vatistas pair:
`scr_p026ef_ctrl_s020v` (job 13591760) / `scr_p026ef_exp_s020v`
(13591761). LineGauss twins (same knobs +
`FLOWPANEL_FILAMENT_REG=linegauss`): `scr_p026ef_ctrl_s020v_lg`
(13591762) / `scr_p026ef_exp_s020v_lg` (13591763). **Ryan's ruling on
record (2026-09-05): if the LineGauss pair showcases the failure (expint
arm still blows up), the campaign default filament regularization
changes to linegauss from then on.**

**Proposed next (Ryan-gated, not launched):** adopt the s020v pair as the
expint-fails validation target. Optional tightening: one cold §12-style
pair `ctrl` vs `WAKE_EXPINT=true` on today's production stack (the 2R
run used the 08-12 code; a rerun on current euler_exp + exact
CoreSpreading composition would close the residual version gap), ~250
steps CPU, hours not days. An expint-fails event validated on current
code becomes the Phase-2 splitting motivation case.

## 17. §16 expint-fails rerun harvest and ruling (2026-09-05)

Four cold s020v rerun arms (§16 candidate 1, launch record in
`phase1b_handoff_prompt_3.md`: jobs 13591760–63, m12 CPU, 12 h wall,
`~/wt026` worktree, 323 steps, full VTP retention). All four arms
IGNITED and died; completion judged from logs + wake-health CSVs, not
sacct (sacct briefly showed ctrl COMPLETED while it was still RUNNING —
the PATH/socket gotcha struck again mid-babysit and produced a false
"ctrl stopped at 208" alarm, later retracted).

Scored from `monitors/scr_p026ef_*_monitor04_wake_health_system1.csv`.
No `*_CT_vs_rev.csv` exists (written post-`simulate!`; no arm finished)
— load agreement uses monitor02 **CFx** (thrust axis is x,
`diagnostic_vertical=(1,0,0)`; CFz is a near-zero lateral component and
must not be averaged for load comparisons).

**Summary** (max over run; min for σ; "u>100" = first step max_u > 100;
first dtZ>2/3 coincided with ignition in all four arms):

| case | max_u | max γ/σ² | min σ | max dtZ | u>100 | death (step, cause) |
|---|---|---|---|---|---|---|
| ctrl_s020v (euler, Vatistas) | 7.7e5 | 7.6e9 | 9.41e-5 | 1.3e4 | 225 | 263, PARTICLE OVERFLOW (500k) |
| exp_s020v (expint, Vatistas) | 4.3e9 | 1.9e16 | 1.94e-5 | 1.9e7 | 210 | 211, DomainError NaN (frozen-gradient ratio, `FLOWVPM_timeintegration.jl:486`) |
| ctrl_s020v_lg (euler, LineGauss) | 3.2e6 | 1.2e11 | 9.41e-5 | 1.5e4 | 230 | 260, PARTICLE OVERFLOW (500k) |
| exp_s020v_lg (expint, LineGauss) | 1.5e3 | 6.5e6 | 4.61e-5 | 3.9e1 | 285 | 300, DomainError Inf (same site) |

**Findings:**

- **Expint-FAILS confirmed on the current stack → §16 candidate 1 is
  VALIDATED as the Phase-2 splitting motivation case.** In the
  resolution-loss regime the exponential integrator does not arrest —
  the opposite of the §13 gpu40/LG events. Vatistas expint ignites
  *earlier* than its euler ctrl (210 vs 225) and dies within one step
  (NaN); LineGauss expint delays ignition 55 steps (285 vs 230) but
  still ignites and dies (+15 steps, Inf). Death mode differs by
  integrator: euler arms grind to the 500k particle cap (+38/+30 steps
  post-ignition; the overflow is intra-step — last wake-health rows
  show 190k/27k particles because post-ignition shedding explodes
  within a step); expint arms die at the frozen-gradient ratio
  the moment the field goes non-finite.
- **Ryan's linegauss ruling (2026-09-05) FIRES: exp-lg blew up.** Per
  the ruling the campaign default filament regularization changes to
  linegauss from then on. Proposed implementation (Ryan scope sign-off
  pending): flip the dispatcher default at
  `examples/run_p018_screen_hpc.slurm.sh:43` from
  `${FLOWPANEL_FILAMENT_REG:-vatistas}` to
  `${FLOWPANEL_FILAMENT_REG:-linegauss}` and update the 025 comment
  block above it; per-case `FLOWPANEL_FILAMENT_REG=vatistas` overrides
  remain available for A/Bs. No other Vatistas pin found in the
  dispatcher (the `_lg` cases already pin linegauss explicitly).
- **Pre-trigger loads agree.** Mean CFx over steps 150–200:
  Vatistas −0.0762 (ctrl) vs −0.0757 (exp), −0.72%; LineGauss −0.0768
  (ctrl-lg) vs −0.0752 (exp-lg), −2.15%.
- **Do not score against the original 019/020 verdicts** (dep stack
  moved 08-24); ignition timing here (225/230 for euler ctrls) vs the
  original refs (213/242) confirms timing shifted, signature did not.

**Pacing note:** ~78–108 s/step at steps 135–169 rising with particle
count to ~167–200 s/step near death; euler ctrl reached step 263 in
~5h55m, exp-lg reached 300 in ~5h50m. A 323-step survivor would have
needed ~7–8 h of the 12 h wall.

**Warm-start brackets (full VTP series retained, verified on orc):**
ctrl_s020v 224/226, exp_s020v 209/211, ctrl_s020v_lg 229/231,
exp_s020v_lg 284/286 (steps below/above each arm's first u>100).

Run dirs `data/scr_p026ef_*` remain on orc unarchived (archive after
Ryan reviews; hpc-storage must NOT touch them before then).

## 18. RK3 confound check on the s020v expint-fails event (2026-09-06)

Arm `scr_p026ef_rk3_s020v_lg` (job 13593720, m12 CPU, `~/wt026` at tag
`campaign/p026ef-rk3-20260905`): cold s020v rerun, same knobs as ctrl-lg
plus `WAKE_INTEGRATOR=rk3`. Banner verified live (`WAKE_INTEGRATOR=rk3`,
`filament regularization = LineGaussRegularization (pinned by
FLOWPANEL_FILAMENT_REG)`). Purpose (§17 follow-up): if RK3 arrests the
event, the failure is time-integration accuracy and the Phase-2
splitting motivation weakens; if it ignites, the resolution-loss
reading stands.

**Ruling: NO ARREST — RK3 ignited far *earlier* than every euler/expint
arm.** First max_u>100 at step **74** (ctrl-lg 230, exp-lg 285); first
dtZ>2/3 the same step; died at step 92's shedding with `ERROR: PARTICLE
OVERFLOW` (500k) after ~1.5 h (~50–55 s/step).

| case | max_u | max γ/σ² | min σ | max dtZ | u>100 | death (step, cause) |
|---|---|---|---|---|---|---|
| rk3_s020v_lg | 1.3e6 | 5.0e9 | 4.71e-5 | 5.1e3 | **74** | 92, PARTICLE OVERFLOW |
| ctrl_s020v_lg (§17) | 3.2e6 | 1.2e11 | 9.41e-5 | 1.5e4 | 230 | 260, PARTICLE OVERFLOW |
| exp_s020v_lg (§17) | 1.5e3 | 6.5e6 | 4.61e-5 | 3.9e1 | 285 | 300, DomainError Inf |

- Signature identical to the §17 family: min-σ collapse (again below the
  euler 9.41e-5 practical floor, matching §15's multi-stage-positions
  observation), γ/σ² → 5e9, dtZ through the ceiling, particle-count
  thrashing (15k→107k→26k) in the last five steps before overflow.
- Pre-trigger loads agree: mean monitor02 CFx over steps 40–70 −0.0824
  (rk3) vs −0.0845 (ctrl-lg), +2.4% — comparable to §17's exp-lg −2.15%.
  Only ~70 clean steps exist given the early ignition.
- The very early ignition (74 vs 230) reads as §15's chaotic-timing
  scatter, amplified by the near-critical s020v regime; §15's rk3_lg arm
  already validated the RK3 wiring (tracked its ctrl for ~65 steps), so
  this is not a wiring fault.

**Consequence: the §17 verdict stands and is strengthened —
order-of-accuracy is falsified as the lever for the resolution-loss
event; Phase-2 splitting proceeds as motivated.** Ryan confirmed the
implementation direction the same day (plan drafted; see
`~/.claude/plans/shimmering-wobbling-sunrise.md`). Log:
`~/wt026/FLOWPanel.jl/logs/slurm/slurm-fp-p026ef-rk3-lg-13593720.{out,err}`;
run dir `data/scr_p026ef_rk3_s020v_lg` (on /home, unarchived; no VTP
warm-start value — ignition fully bracketed but the event class is
already covered by the §17 arms).

---

# §19 Design amendment — per-mechanism fractional gating, σ bounds as clamps (Ryan rulings 2026-09-08)

Ryan reviewed the shipped trigger set and replaced it. Four rulings
(AskUserQuestion, 2026-09-08), all implemented the same day on the live
checkouts:

1. **Gate metric**: each mechanism fires on its OWN accumulated growth
   fraction relative to the particle's `sigma_0`, using the per-mechanism
   Δσ² accumulators (`dvisc`, `drvpm`) that previously only routed:
   mechanism k fires when `sqrt(σ₀² + Δσ²_k)/σ₀` leaves
   `[1 − f_elong, 1 + f_k]` on its side. The trigger IS the mechanism —
   the `dvisc ≥ drvpm` routing tie-break is deleted.
2. **Shrink side unified**: accumulators record ATTEMPTED (pre-clamp) Δσ²,
   so floor/ceil-pinned particles keep accruing credit and **triggers still
   fire at the clamps**. This subsumes and deletes the exposure integral
   (`log_stretch_max`) and the floor-pin trigger (both were workarounds for
   realized-σ starvation on the floor).
3. **rVPM accumulator is signed net**: compression (+) and elongation (−)
   cancel; a particle that compresses then relaxes back never splits.
   Net > 0 → tri3, net < 0 → pair2.
4. **σ_max/σ_min are clamps, not triggers**: the integrator sigma_guard
   floor/ceil stays as a permanent partner of splitting (the "band-aid"
   framing and the SIGMA_CEIL-vs-splitting mutual-exclusion error are
   retired — **commit 8 (§8 SIGMA_CEIL removal) is CANCELLED**), and
   emitted children are additionally clamped into `[sigma_min, sigma_max]`
   at emission (vol recomputed from the clamped σ_c).

Anti-refire stays cooldown-free: `_rsplit_reset_slot!` restamps
`sigma_0 = σ_c` and zeroes accumulators on every child and merged
representative, so a still-pinned child re-arms only after fresh attempted
deformation.

Knob surface (old → new):
`ResolutionSplitOpts.sigma_max/sigma_growth_ratio_max/log_stretch_max/sigma_floor`
→ `f_visc/f_comp/f_elong` (NaN-disabled fractions) + `sigma_min/sigma_max`
(NaN-unbounded emission clamps). Driver:
`WAKE_SPLIT_SIGMA_MAX`(trigger)/`_SIGMA_GROWTH_RATIO_MAX`/`_LOG_STRETCH_MAX`/
`_ON_FLOOR` → `WAKE_SPLIT_FRAC_VISCOUS`/`_FRAC_COMPRESS`/`_FRAC_ELONGATE` +
`WAKE_SPLIT_SIGMA_MIN`/`_SIGMA_MAX` (clamps, defaulting to the 052c
sigma_guard floor and `SIGMA_CEIL` so the two clamp layers agree).
`ResolutionSplitState` drops `exposure` (state file: sigma_0, axis, weight,
dvisc, drvpm); particle-VTK persistence is now 5 `rsplit_*` fields, with a
warn-and-zero-accumulators migration path for legacy 6-field saves
(old accumulators recorded applied post-clamp Δσ², not comparable).

Consequences for pending work:
- **Commit 7 arms** (§8.4 cap030/cap018, §9 s020v matrix): still Ryan-gated;
  the cap arms' `WAKE_SPLIT_SIGMA_MAX` is now the emission clamp and the
  operating point moves to `WAKE_SPLIT_FRAC_COMPRESS` (re-derive vs the §5
  adequacy table so splits fire before particles sit long at the clamp);
  shrink arms move from `WAKE_SPLIT_LOG_STRETCH_MAX`+`WAKE_SPLIT_ON_FLOOR`
  to `WAKE_SPLIT_FRAC_ELONGATE` (D4 value must be re-derived in fraction
  space). Dispatcher comments updated in `run_p018_screen_hpc.slurm.sh`.
- **Commit 8 is cancelled** (ruling 4); D5 closed accordingly.
- GPU gap note (§ device paths): the exposure-trigger half of the gap is
  gone with the trigger; the accumulator half remains — device-resident
  integrator twins still do not maintain `dvisc`/`drvpm`, so splitting
  remains host-mirror-only.

# §20 Separation-criteria audit — child spacing/count vs vortex-tube physics (2026-09-09)

Audit of the split kernels' *geometry* (how far apart, how many children)
against the physics that fires the §19 triggers. Elongation first (the
production-relevant mechanism). All facts pinned from the live checkouts.

## 20.1 Conventions and code facts

- **Overlap convention**: $\Phi \equiv \sigma/h$ (spacing $h = \sigma/\Phi$).
  Pinned at `FLOWPanel_wake.jl:730` (`sigma = dist*overlap/p_per_step`) and
  `:2113` (`h = sigma/overlap`). Driver `OVERLAP` default 3.0; the 018
  campaign convention is 2.75. A shed particle is born tiling filament
  length $\ell_0 = \sigma_0/\Phi$ ($\approx 0.36\,\sigma_0$ at $\Phi=2.75$).
- **rVPM σ-law** (`FLOWVPM_timeintegration.jl` MM4; 020
  `phase_01_theory.md` §2.1): $\dot\sigma = -\sigma Z$ with $Z = h_\sigma s$
  and $h_\Gamma = 2h_\sigma$ at $(f,g)=(0,1/5)$, so $|\Gamma|$ grows at rate
  $2Z$ while $\sigma^2$ decays at rate $2Z$. Since circulation is conserved
  along the tube, $|\Gamma| \propto \Gamma_{\rm circ} L \propto L$: **the
  represented segment length grows exactly as $|\Gamma|$, and
  $\sigma^2 L = \mathrm{const}$ holds exactly for the rVPM channel.** Hence

  $$\lambda \equiv \frac{L}{L_0} = \frac{\sigma_0^2}{\sigma_{\rm att}^2},
  \qquad \sigma_{\rm att}^2 = \sigma_0^2 + \texttt{drvpm},$$

  and at the elongation trigger $\lambda^* = (1-f_{\rm elong})^{-2}$
  (e.g. 2.04 at $f_{\rm elong}=0.3$). Caveat: the *material* line stretch is
  $e^{\int s\,dt}$; the represented tube stretches as $e^{2\int Z\,dt} =
  e^{2h_\sigma\int s\,dt}$ — at $h_\sigma = 1/5$ only 2/5 of kinematic
  stretching becomes represented length. $\lambda$ from `drvpm` measures the
  represented tube, which is the right quantity for re-discretization.
- **"Attempted" λ under clamps**: each accumulation step records
  $\Delta\sigma^2$ computed *off the realized (clamped) σ*, so `drvpm` is a
  chained linearization, not the free trajectory. Floor-pinned at
  $\sigma_f$: per-step $\Delta\sigma^2 \approx -2\sigma_f^2 Z\,dt$, giving
  $\lambda_{\rm att} = [1 - 2(\sigma_f/\sigma_0)^2\!\int\!Z\,dt]^{-1}$ vs
  the free $\lambda = e^{2\int Z dt}$ — first-order equal, with
  $\lambda_{\rm att}$ running high for long pinned windows at
  $\sigma_f \approx \sigma_0$ and low for $\sigma_f \ll \sigma_0$. Adequate
  for split sizing (splits fire at small fractions, so windows are short).

## 20.2 Elongation (pair2) audit — VERDICT: fixed spacing is arbitrary

Current kernel: $m=2$ children at $\pm 0.5\,\sigma_p$ (spacing
$s = 1.0\,\sigma_p$), $\sigma_c = \sigma_p$, $\Gamma_c = \Gamma_p/2$,
`circulation` unchanged (correct — series division of a tube keeps
circulation).

- **Tiling mismatch.** The parent represents $\lambda\ell_0$ of tube; $m$
  children should be spaced $s_{\rm phys} = \lambda\ell_0/m$. Unclamped
  ($\sigma_p = (1-f)\sigma_0$, $m=2$):

  $$\frac{s_{\rm current}}{s_{\rm phys}}
  = \frac{1.0\,(1-f)\sigma_0}{\lambda\sigma_0/(2\Phi)}
  = 2\Phi(1-f)^3 .$$

  At $(\Phi,f)=(2.75,0.3)$: **1.89 — children placed ~89 % too far
  apart**; crossover at $f\approx 0.43$; at $f=0.5$ the same constant
  UNDER-covers (0.69). No fixed ratio matches the physics as $f_{\rm elong}$
  varies — confirming the coupling the reset prompt suspected.
- **Overlap.** Child–child overlap $\Phi_{cc} = \sigma_c/s = 1.0$ for the
  defaults vs the 2.75 shedding convention — under-overlapped by 2.75×
  relative to birth discretization (but see §20.3: merging, not shedding,
  sets the wake's effective spacing).

**Matched design (closed form, no new state).** Choose target overlap
$\Phi_t$; child spacing $s_t = \sigma_c/\Phi_t$; child count from tiling
$\lambda\ell_0$ with $\ell_0 = \sigma_0/\Phi_t$:

$$m^\* = \frac{\lambda_{\rm att}\,\ell_0}{s_t}
       = \lambda_{\rm att}\,\frac{\sigma_0}{\sigma_c}
       \;\;\xrightarrow{\ \sigma_c=(1-f)\sigma_0,\ \text{unclamped}\ }\;
       (1-f_{\rm elong})^{-3} = \lambda_{\rm att}^{3/2}.$$

Properties (all exact pre-rounding):

1. **The volume rule emerges**: $m^\*\sigma_c^3 = \sigma_0^3$ — the
   overlap-matched tiling is precisely merging's inverse
   ($\sigma = \sqrt[3]{\Sigma\sigma^3}$), restoring the σ³-consistency the
   W5 hygiene rule gave up.
2. **Self-consistent without per-particle $\ell_0$ state**: a child's
   implied birth tile $\sigma_c/\Phi_t$ equals its true tile
   $\lambda\ell_0/m^\*$ identically — including when floor-pinned
   ($\sigma_c = \sigma_f$ enters both sides). The NEW-state option in the
   reset prompt is NOT needed; $\Phi_t$ as a knob suffices.
3. **Composition is exact**: $m_1 m_2 = \lambda_1\frac{\sigma_0}{\sigma_1}
   \cdot \lambda_2\frac{\sigma_1}{\sigma_2} = \lambda_{\rm tot}
   \frac{\sigma_0}{\sigma_2}$ — repeated small splits produce the same
   child count, σ, and tiled span as one big split. Only integer rounding
   breaks this, boundedly.
4. Numbers at $(f,\Phi_t,\sigma_0)=(0.3, 2.75, \sigma^\*{=}0.0381R)$:
   $\lambda^\* = 2.04$, $m^\* = 2.92 \to 3$ children, spacing
   $s = \lambda\ell_0/3 = 0.247\,\sigma_0 = 0.0094R$.

## 20.3 The binding constraint is MERGING, not the tube physics

Production merge (driver lines 108–110, 856–859): every step,
`sigma_relative=false`, radius $r_m = 0.02R$ ABSOLUTE (= $0.525\,\sigma^*$),
pairing = nearest-within-radius (no Γ-alignment gate; `gamma_align_cos`
defaults off), representative $\sigma = \sqrt[3]{\Sigma\sigma^3}$.

- $r_m = 0.02R$ **exceeds the shed spacing** $\sigma^\*/\Phi = 0.0139R$:
  the wake's effective spacing floor is set by merging, not by `OVERLAP`.
  The overlap "maintained everywhere else" is at most
  $\sigma^\*/r_m \approx 1.9$, less after merge σ-growth.
- Overlap-matched children ($s = 0.0094R \ll r_m$) are **immediate merge
  candidates** (identical σ, mutual nearest). Split leaves σ unchanged;
  merge returns $\sigma = m^{1/3}\sigma_c$. A split→merge cycle is
  therefore a **σ-pump**: ×$m^{1/3}$ per cycle (+26 % at $m=2$, +44 % at
  $m=3$) with both states re-armed fresh each time — positive feedback,
  not mere churn. The §19 anti-refire argument covers accumulators only;
  it cannot prevent this.
- The current $1.0\,\sigma_p$ spacing ($\approx 0.7\sigma_0 = 0.0267R$ at
  $f=0.3$, $\sigma_0=\sigma^*$) clears $r_m$ by only 34 % — the shipped
  constant sits near the minimum merge-safe spacing, apparently by
  accident. Note the margin shrinks as $\sigma_p$ floors: at
  $\sigma_p = \sigma_f < 0.02R/1.0$ the CURRENT kernel is merge-unsafe too.

Resolution options (Ryan): (a) **merge-safe spacing floor** in the kernel:
$s \ge \kappa\,r_m$ (κ ≈ 1.2), reducing $m$ to
$\max(2, \lfloor\lambda\ell_0/(\kappa r_m)\rfloor)$ when the floor binds —
physics-tiling whenever merging permits, graceful degradation to
pair2-at-safe-spacing otherwise (requires passing the merge radius to
`ResolutionSplitOpts`); (b) merge exemption for fresh children (age state —
against the no-cooldown ruling); (c) Γ-alignment gate on merging (would
also stop the wholesale wake coarsening production currently relies on);
(d) shrink $r_m$. Recommendation: (a).

## 20.4 Compression (tri3) — analysis only

Tube picture: `drvpm` > 0, $\lambda < 1$; cross-sectional area grows by
$1/\lambda = (1+f_{\rm comp})^2$ at trigger. Re-discretizing the fattened
bundle into children of birth-sized cores needs
$m \approx (1+f_{\rm comp})^2$ filaments: $m=3$ is matched to
$f_{\rm comp} \approx 0.73$; at the likelier 0.3–0.5, $m=2$ suffices. So
yes — child count should scale with accumulated compression, with exponent
1 in $\sigma_{\rm att}^2/\sigma_0^2$ (2-D cross-section) vs the elongation
$\lambda^{3/2}$ (1-D length at fixed cross-section). The FIXED ring radius,
unlike pair2's spacing, already scales with realized compression through
$\sigma_p$; it under-represents attempted fattening only when
ceiling-pinned ($\sigma_p < \sigma_{\rm att}$ — substituting
$\sigma_{\rm att} = \sqrt{\sigma_0^2 + \texttt{drvpm}}$ for $\sigma_p$ in
the radius would fix that if it ever matters). Ring children at the default
$0.6\sigma_p$ (spacing $1.04\,\sigma_p$, grow-side $\sigma_p \ge
(1{+}f)\sigma_0$) are merge-safe under production numbers ($\ge 0.05R$ at
$f{=}0.3$, $\sigma_0{=}\sigma^*$). D2 (0.6 vs moment-match 1.155) stays
open pending the §3a kernel-fit study. **Hygiene flag**: tri3/tetra4 pass
the parent's `circulation` to every child, but parallel-filament division
should carry `circ`$/m$ (pair2's unchanged `circ` is correct — series
division). `circulation` is diagnostic-only today; fix opportunistically.

## 20.5 Viscous (tetra4) — audit note only

Offset $1.3503\,\sigma_p$ is the per-axis second-moment match (§3a);
spacing $\approx 3.5\,\sigma_c$ under-overlapped as documented. Zero
production events; D3 stays gated on the §3a kernel-fit study. No change.

## 20.6 Proposed elongation kernel (pending Ryan's rulings)

In-line adaptive-$m$: $m = \mathrm{clamp}(\mathrm{round}(
\lambda_{\rm att}\sigma_0/\sigma_c),\ 2,\ m_{\max})$, spacing
$s = \lambda_{\rm att}\ell_0/m$ (exact tiling; $\approx \sigma_c/\Phi_t$),
subject to the §20.3 merge-safe floor; symmetric offsets
$\big(k - \tfrac{m+1}{2}\big)s$, $k=1{:}m$ along the averaged axis;
$\Gamma_c = \Gamma_p/m$, $\sigma_c = \sigma_p$, `circulation` unchanged.
Γ-total, centroid, linear impulse exact for any $m$; angular impulse exact
by symmetry. Knobs per the established pattern: `elongate_overlap`
($\Phi_t$; NaN → legacy fixed-ratio pair2), `elongate_m_max`,
merge-floor passthrough; driver `WAKE_SPLIT_ELONGATE_OVERLAP` /
`_ELONGATE_M_MAX`. Capacity check becomes $m-1$ appended slots; verbose
counters gain a children-emitted tally. Tests: extend t1/t2/t9 to adaptive
$m$; new composition test (two small splits ≡ one big: same $m\sigma^3$
and span); overlap + merge-floor assertions; re-derive t3/t4 pins for
$m>2$ in-line geometry.

## 20.7 Rulings and implementation (Ryan, AskUserQuestion 2026-09-09)

Theory now lives in `splitting_theory.md` (same directory) — a SINGLE living
draft, updated in place with no change history (Ryan's requested format);
this design doc remains the decision/history record. Rulings on the §20
audit, all implemented the same day on the live checkouts (uncommitted, on
top of the still-uncommitted §19 redesign):

1. **Adaptive in-line elongation kernel SHIPPED** (`_split_elongate_line!` +
   `_elongate_plan`): `m = clamp(round(λ_att·σ₀/σ_c), 2, elongate_m_max)`
   children tile the accumulated stretch at target child overlap
   `elongate_overlap` (Φ_t). New `ResolutionSplitOpts` fields
   `elongate_overlap` (NaN → legacy fixed pair2, the FLOWVPM default) and
   `elongate_m_max` (default 4); driver `WAKE_SPLIT_ELONGATE_OVERLAP`
   (defaults to the shedding `OVERLAP`, so the driver runs adaptive by
   default) and `WAKE_SPLIT_ELONGATE_M_MAX`. `split_particles!` return
   gained `n_children_elongate`.
2. **Merge–split interplay → overlap-gated merging as an OFF-by-default
   knob**: `MERGE_OVERLAP=Φ_merge` flips the driver's MergeParticles to
   `sigma_relative=true, r = 1/Φ_merge` (merging already supports the
   criterion natively — no FLOWVPM change). Recommended Φ_merge = 3.5 >
   Φ_t = 2.75 kills the §20.3 σ-pump by construction; production keeps the
   absolute 0.02R radius pending an A/B (the gate merges strictly less, so
   counts/cost rise).
3. **Compression stays m = 3** even though the count-match is
   `m ≈ (1+f_comp)²` (triangle matched to f_comp ≈ 0.73): m = 2 would
   impose artificial transverse anisotropy on an axisymmetric fattening.
   Noted in the kernel docstring + theory doc §5 that higher m with a more
   complicated child shape is the future extension; this effectively pairs
   the mechanism with f_comp ≈ 0.73 when count-matching matters.
4. **circulation bookkeeping fixed**: tri3/tetra4 children now carry
   `circ/3` / `circ/4` (lengthwise bundle division splits the vorticity
   flux); elongation kernels keep `circ` unchanged (crosswise cut).
   Diagnostic-only field, no dynamics change.

Correction to §20.6 as written: the elongation kernel's angular impulse is
NOT exact by ± symmetry (quadratic offset terms don't cancel); it obeys the
t2 bound `|ΔA| ≤ a²|Γ|/3`, `a = (m−1)s/2`, vanishing when the axis ∥ Γ.
Tests: t10 (plan math, saturation, conservation + geometry + angular bound,
exact composition, capacity, legacy fallback, validation) and t11
(circulation shares) added to `runtests_resolution_split.jl` — suite
1011/1011 green.

# §21 Launch-prep rulings and GPU directive (Ryan, 2026-09-11)

Cost models re-derived against the §20.7 adaptive kernel (launch-prep
session, 2026-09-11):

- **Elongation cost is now fraction-independent.** Exact composition +
  volume rule make the full-descent multiplier
  $(\sigma_0/\sigma_{\rm floor})^3$ — ≈64× to the 0.25σ₀ floor regardless
  of $f_{\rm elong}$ (~72× at $f=0.3$ with rounding; plan unclamped for
  $f \le 1-4^{-1/3} \approx 0.37$). Against measured populations: healthy
  field (60/180k below floor) +~4k (+2%); cs0p002-style ignited tail
  (4401) +~277k (~1.5× field), and firing continues past 64× at m_max per
  fire while attempted collapse keeps accruing (capacity `break` is the
  backstop). $f_{\rm elong}$ therefore controls granularity/timing, not
  total cost; the cost levers are `elongate_m_max` (m_max=2 reproduces
  pair2 economics, ~15× at $f=0.3$) and the merge policy.
- **Compression bound unchanged**: σ_c = σ_p/√3 keeps the self-limiting
  bound $f_{\rm comp} < \sqrt3-1 \approx 0.732$; per §20.7 ruling 3, m=3 is
  count-matched at that same fraction, making $f_{\rm comp}=0.73$ the
  principled operating point (σ-neutral per fire, matched cross-section
  re-discretization).

**RULED (Ryan 2026-09-11): settings adopted for the gated arms** —
`WAKE_SPLIT_FRAC_ELONGATE=0.3` (m=3/fire, confirmed firing in the
2026-09-10 driver smoke), `WAKE_SPLIT_FRAC_COMPRESS=0.73`, adaptive
elongation at the default Φ_t (= shedding `OVERLAP`), and a
`MERGE_OVERLAP=3.5`-vs-production-absolute A/B as the first merge
discriminator (exp bracket first).

**RULED (Ryan 2026-09-11): implement splitting and merging on GPU before
launch.** The 09-11 session re-confirmed splitting is silently CPU-only
(device integrator twins skip all `_rsplit_accumulate*`; no fail-fast), and
per-arm GPU capability is wanted per the standing GPU-default preference.
Scoping fact: GPU maintenance already runs merge + the split pass on the
host mirror (`_apply_particle_maintenance_device!`,
`src/FLOWPanel_gpu_wake.jl`); the gap is ONLY device-side accumulation of
the trigger/direction state and its D2H sync. Handoff:
`gpu_split_merge_reset_prompt_20260911.md` (this directory). Campaign
launches remain Ryan-gated.
