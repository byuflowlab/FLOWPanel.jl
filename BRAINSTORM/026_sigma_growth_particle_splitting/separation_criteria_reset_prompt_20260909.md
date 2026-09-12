# 026 reset prompt — audit the split SEPARATION criteria against vortex-tube physics (2026-09-09)

You are picking up BRAINSTORM 026 (resolution-preserving particle splitting)
in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (panel solver) +
`/Users/ryan/Dropbox/research/projects/FLOWVPM.jl` (VPM core). Read
`CLAUDE.md` and the routing policies it names before touching code.

## Where the item stands (2026-09-09, all UNCOMMITTED on both live checkouts)

The fractional-gating redesign (Ryan's four 2026-09-08 rulings) is
implemented, tested, and smoke-verified but **not committed**:

- Triggers are per-mechanism growth fractions vs each particle's `sigma_0`:
  `ResolutionSplitOpts.f_visc / f_comp / f_elong` (NaN-disabled). Mechanism k
  fires when `sqrt(σ₀² + Δσ²_k)/σ₀` crosses `1 + f_k` (grow) or `1 − f_elong`
  (shrink). The trigger IS the mechanism — no routing tie-break.
- `ResolutionSplitState` = `sigma_0, axis, weight, dvisc, drvpm` (the
  `exposure` integral is gone). `dvisc`/`drvpm` accumulate **attempted
  (pre-clamp) Δσ²** since the particle's last split; `drvpm` is **signed
  net** (rVPM compression +, elongation −). Clamp-pinned particles keep
  accruing, so triggers still fire at the sigma_guard floor/ceiling.
- `sigma_min`/`sigma_max` are CLAMPS (emission clamp in
  `_rsplit_emit_child!` + the integrator sigma_guard, now a permanent
  partner of splitting; commit 8 = SIGMA_CEIL removal is CANCELLED).
- Full FLOWVPM `Pkg.test` green; FLOWPanel wake/replay/rsplit-warm-start
  suites green; rotor-hover driver smokes verified (config errors, split
  telemetry, SIGMA_CEIL coexistence). Known unrelated failure:
  `test/runtests_unit_warmstart.jl` FIRST testset dies pre-existing
  (immutable `WarmstartNoopSolver` as WeakKeyDict key in
  `_publish_block_gs_status!`) — not yours.
- History/design record: item doc
  `BRAINSTORM/026_sigma_growth_particle_splitting/particle_splitting_design.md`,
  **§19** (this redesign) on top of §3a/§3b (kernel geometry + algebra) and
  the D2/D2b/D3 open offset-ratio questions. Commit-7 campaign arms remain
  Ryan-gated with fractions to re-derive.

Key files: `FLOWVPM.jl/src/FLOWVPM_resolution_split.jl` (state, opts,
kernels `_split_viscous_tetra4!` / `_split_compress_tri3!` /
`_split_elongate_pair2!`, trigger `_rsplit_check`, main pass
`split_particles!`); accumulation sites in `FLOWVPM_timeintegration.jl`
(rVPM σ update + euler_exp geometric step) and `FLOWVPM_viscous.jl`;
FLOWPanel plumbing in `src/FLOWPanel_wake.jl` (`ResolutionSplit` policy,
VTK), `src/FLOWPanel_warmstart.jl`, driver knobs in
`examples/rotor_hover_pressure_comparison.jl` (`WAKE_SPLIT_*`); tests in
`FLOWVPM.jl/test/runtests_resolution_split.jl`.

## Your task: audit the SEPARATION criteria (child spacing/geometry) against the physics

The triggers now encode *when* to split. This session asks whether the
kernels' *geometry* — how far apart children are placed, and how many —
matches the physics that fired the trigger. Work mechanism by mechanism,
elongation first (it is the production-relevant one).

### 1. Elongation (pair2) — the central question

Current kernel: exactly 2 children at `±0.5·σ_p` along the averaged stretch
axis (spacing `1.0·σ_p`), `σ_c = σ_p`, `Γ_c = Γ_p/2`, regardless of how much
stretching accumulated. Interrogate this against the vortex-tube picture:

- **Does the spacing match the accumulated stretch?** For an incompressible
  tube, `σ²·L = const`, so the attempted net accumulator gives the length
  stretch factor directly: `λ = σ₀² / (σ₀² + drvpm)` (with `drvpm < 0` on
  the shrink side, `λ > 1`; at the trigger, `λ = 1/(1 − f_elong)²`). The
  particle's represented tube segment has grown from `L₀` to `λ·L₀` —
  where `L₀` is the inter-particle spacing the discretization was born with
  (tie this to the campaign overlap convention; verify what OVERLAP=2.75
  means in this codebase — σ/h vs h/σ — before using it). Does a fixed
  `1.0·σ_p` spacing re-discretize that length, or is it an arbitrary
  constant that under- or over-covers depending on `f_elong`?
- **Is child overlap retained?** Children at spacing `s` with cores `σ_c`
  must keep `s/σ_c` within the overlap the campaign maintains everywhere
  else; quantify the post-split overlap for the current defaults and for a
  physically-matched spacing.
- **Ryan's explicit opening (2026-09-09): more than 2 children is
  acceptable, even an adaptive number.** Design (and, if the audit supports
  it, implement) an elongation kernel where the child count `m` and spacing
  are chosen together so that (a) the children tile the stretched length
  `λ·L₀` implied by the accumulated attempted stretch since the last split,
  AND (b) neighboring children retain the target overlap. E.g.
  `m = clamp(ceil(λ·L₀ / s_target), 2, m_max)` with `s_target` set by the
  overlap convention, children in-line along the axis, `Γ_c = Γ_p/m`,
  `σ_c = σ_p` (cross-section untouched), centroid/±-symmetric placement so
  Γ-total, centroid, linear impulse stay exact (angular impulse by
  symmetry). Mind: capacity headroom becomes `m−1` appended slots; the
  verbose counters and tests assume fixed child counts; `f_elong` and the
  spacing rule are now coupled (a larger trigger fraction should produce
  proportionally more/farther children — check the design stays consistent
  as `f_elong` varies).
- Check the interaction with `sigma_min` clamping and with merging (children
  spaced at overlap are exactly what `MergeParticles` likes to re-merge —
  is there a ping-pong risk at the chosen spacing? The anti-refire argument
  only covers the accumulators, not merge–split cycles).

### 2. Compression (tri3) and viscous (tetra4) — same lens, second priority

- tri3: ring radius `0.6·σ_p` (spacing `1.8·σ_c`) was a compromise flagged
  under-overlapped in §3b (D2: moment-match would be 1.155). Under the tube
  picture, compression shortens/fattens the segment — does a FIXED ring
  radius match the accumulated compression `λ < 1`, and should the ring
  radius (or child count) scale with `sqrt(1 + drvpm/σ₀²)` the same way the
  elongation spacing scales with λ?
- tetra4: offset `1.3503·σ_p` is a second-moment match (D3, gated on a §3a
  kernel-fit study; mechanism has zero observed production events — audit
  only, don't invest in implementation).
- Deliverable here can be analysis + recommendation; implement only what the
  elongation redesign makes natural to share.

### 3. Ground rules

- Derive before you code: write the σ²L-conservation algebra for each
  regime, including what "attempted" accumulation means for λ when clamps
  engaged mid-window, and check the limit behaviors (λ→1⁺ should reproduce
  something sane; repeated small splits should compose to the same
  discretization as one big split — test that composition property).
- Conservation invariants are non-negotiable: total Γ, centroid, linear
  impulse exact; angular impulse error bounded (see t1/t2 in
  `runtests_resolution_split.jl` — extend them to adaptive m).
- Far-field equivalence and support-residual pins (t3/t4) must be
  re-derived, not just loosened, if geometry changes.
- Keep the trigger layer untouched unless the audit finds the trigger and
  spacing cannot be made consistent — flag that for Ryan instead of
  changing rulings.
- New knobs follow the existing pattern: NaN/disabled defaults in
  `ResolutionSplitOpts`, `WAKE_SPLIT_*` env plumbing in the driver, VTK/
  warm-start untouched unless per-particle state is added (if you need the
  birth spacing `L₀` per particle, that is a NEW `ResolutionSplitState`
  field → follow the lockstep-hook + VTK + warm-start + legacy-migration
  pattern just established for the exposure removal).
- Anything ambiguous about the overlap convention, `m_max`, or whether to
  ship vs analyze-only: AskUserQuestion Ryan. Campaign launches stay
  Ryan-gated. Local runs ≤ 4 threads. Commit only when Ryan says so —
  note the fractional-gating redesign itself is still uncommitted, so
  coordinate with Ryan on commit sequencing before layering more changes.

### Verification expectations

Unit: extended t1–t5/t8/t9 for the new geometry (including adaptive-m
composition and overlap assertions). Integration: the sane-fraction driver
smoke (`WAKE_SPLIT_STRETCH=true WAKE_SPLIT_FRAC_COMPRESS=0.5
WAKE_SPLIT_FRAC_ELONGATE=0.3 SIGMA_CEIL=1.0 NREVS=0.2 RUN_MONITORS=false
julia --project examples/rotor_hover_pressure_comparison.jl`) must complete
with sane counters. Document findings + any new algebra as a dated §20
appendix in the item design doc (append-only).
