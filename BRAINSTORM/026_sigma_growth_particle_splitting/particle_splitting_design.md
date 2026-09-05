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
