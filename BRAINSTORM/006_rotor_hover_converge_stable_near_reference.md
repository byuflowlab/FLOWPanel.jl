# Rotor Hover: Converge to a Stable AND Near-Reference CT

## Goal

Drive the particle-wake rotor-hover prediction to a CT history that is **both**:

1. **Stable** — non-oscillating (or nearly so) over many revolutions at hover
   (residual `magVinf ≈ 0`), peak-to-peak ripple driven toward zero; and
2. **Near the reference magnitude** — CT in the **0.068–0.072** band, evaluated
   at a **validated loose truncation (≥4R)** so the level is trustworthy.

This is the convergence counterpart to item 005. **005 deliberately accepts a
flat plateau at *any* value** (oscillation amplitude is its only metric); **006
additionally requires the plateau to land near the experimental/BEM reference.**
The two halves are coupled: a stable history at the wrong magnitude (005) and a
brief pass through 0.072 that then drifts (004) are each only half the answer.

Temporal + grid convergence (insensitivity to `NT` / `P_PER_STEP` / mesh) is
**out of scope here** and split to a future item (007), per the user — chased
only once stable+near-reference is in hand.

## Context

Items 001–003 cleared backend, Kutta, geometry, airfoil, and Reynolds number as
causes and routed the panel shortfall (panel CT ≈ 0.0506 vs experiment ≈ 0.072,
BEM 0.060–0.071) to **wake / induced-velocity modeling**. Item 004 (still an
unexecuted audit plan) flags that wake-bearing runs pass through `CT ≈ 0.072` and
then drift. Item 005 built the staged-startup machinery in
`examples/rotor_hover_pressure_comparison.jl` (4-phase smoothstep freestream pulse
+ ENV hooks; downstream `GlobalCylinder` truncation; plateau diagnostics) and
showed the **staged startup kills the impulsive startup ring**, but withdrawing
the freestream to hover re-introduces a CT oscillation whose grow/steady/decay
character is **not yet established**.

## Dependency — gated on 005 Phase 2

**006's specific sweeps are deliberately deferred until 005's Phase-2 plan
(`plans/20260617_brainstorm_005.md`, experiments E1–E4) has been executed.** That
phase is the execution vehicle that 006 consumes; 006 does not duplicate it. 006
*begins* by reading the E1–E4 outcomes and only then selects which levers to pull
toward the near-reference target:

- **E1 (long run):** is the post-withdrawal oscillation growing / a limit cycle /
  a decaying long transient, and what is its period?
- **E2 (truncation depth 1/2/4R):** is the oscillation a deletion-front boundary
  artifact, and is tight depth a trustworthy fast proxy?
- **E3 (withdraw-rate sweep + pure-spinup control):** is the oscillation
  withdrawal-excited (ramp is the lever) or self-sustained by the hover wake
  (damping is the lever)?
- **E4 (damping menu):** ranked detrended-amplitude reduction for each physical
  and numerical damper, applied on the best non-damping baseline.

The verdict from E1–E4 decides 006's lever set: if oscillation is
withdrawal-excited → refine the ramp; if self-sustained → apply the top-ranked
physical damper(s); if a truncation artifact → re-baseline at ≥4R first.

## Candidate levers + expected effects (selection deferred)

Carried over from 005 Phase 2; ranked physical-first. **Which to pull is decided
after E1–E4**, not now.

- **Truncation depth (`TRUNCATION_DEPTH_R`, 1/2/4R).** Dual role: boundary-artifact
  diagnostic *and* magnitude gate. Tight depth inflates CT (deletes far-wake
  downwash) — so the **near-reference magnitude judgment must be made at ≥4R**.
- **Velocity ramp-up-then-down / withdraw rate (`FREESTREAM_WITHDRAW_REVS`, plus a
  pure-spinup `MAGVINF_PEAK≈0` control).** Tests withdrawal-excited vs
  self-sustained. Slower withdraw shrinking the amplitude ⇒ the ramp is the lever.
- **VPM relaxation scheme.** *Only on/off (`PARTICLE_RELAX`) is exposed today.* The
  real "different scheme" knob — std vs **corrected**-Pedrizzetti, `rlxf`,
  `nsteps_relax` — requires a **small new ENV hook** in
  `examples/rotor_hover_pressure_comparison.jl`, wired into how the
  `PanelParticleWake`'s `FLOWVPM.Relaxation` is constructed (relaxation is applied
  at `src/FLOWPanel_wake.jl:1054`, `FLOWVPM._euler(w.pfield, dt; relax)`; the
  schemes are `relax_pedrizzetti` / `relax_correctedpedrizzetti` in
  `FLOWVPM_relaxation.jl`). Classed as **numerical (second-resort)** damping per
  005 — it realigns Ω with local velocity at some physical cost.
- **Physical knobs.** Viscous diffusion (`WAKE_NU` ×3/×10 — most physical, on by
  default at ~1.69e-5/ρ); SFS strength (`SFS_MAXC`/`SFS_RLXF`, on by default);
  kernel overlap / core-beta (`OVERLAP`, `WAKE_CORE_BETA`); merging
  (`MERGE_R_FACTOR`).
- **Resolution (`NT`, `P_PER_STEP`).** Diagnostic as much as cure: refine-removes
  ⇒ under-resolution artifact; refine-sharpens/persists ⇒ a real hover-wake
  instability (then damp, don't "fix"). Belongs to 007's grid-convergence scope
  but is also the cleanest discriminator here.

## Magnitude-recovery angle (006-specific)

Unlike 005, 006 must land *near 0.072*, so the magnitude half leans on the
wake-physics levers 003/004 flagged, not just oscillation damping:

- 003 showed BEM brackets experiment (0.060–0.071) on identical geometry/airfoil/Re
  — lower-ncrit (noisier inflow, physical for a rotor ingesting its own wake)
  raised CT toward 0.072. The panel analogue is a **more energetic, less-dissipated
  near-disk induced field**, which argues *against* over-damping.
- 001's semi-infinite rigid wake **overshot to 0.110** (+118%) — confirming wake
  extent/induced velocity is the dominant magnitude lever, but also that "more
  wake" is easy to overshoot. The near-reference target is a balance: enough
  developed wake column to reach ~0.072, dissipated/relaxed just enough to stop
  ringing — **without** truncating away the far-wake downwash (hence ≥4R) or
  damping away the near-disk inflow that sets the level.

This tension (damp to stabilize vs preserve inflow to hit magnitude) is the core
of 006 and is revisited concretely once the E4 ranked damping menu exists.

## Execution constraints

- **Thread budget: 60 total** (workstation; raised from 005's 36). Partition so
  `CONCURRENCY × THREADS_PER_RUN ≤ 60`. Suggested: E2 (3 cfgs) 3×20, E3 (4 cfgs)
  4×15, E4 (~12 cfgs) 6×10 in two waves; E1 single long run one 20–24-thread job.
- **OMP/BLAS tied to each run's Julia thread count** (so BLAS never exceeds
  budget): each job exports `OMP_NUM_THREADS` / `OPENBLAS_NUM_THREADS` /
  `MKL_NUM_THREADS` / `VECLIB_MAXIMUM_THREADS` = its `julia --threads N`.
  `examples/run_rotor_hover_stable_wake.sh` already derives all of these from
  `THREADS_PER_RUN`; only the budget cap changes 36→60.
- `BERNOULLI_ONLY=true` for tuning; confirm conclusions with the full pressure
  monitors and at ≥4R truncation once. Keep VTK on (CLAUDE.md disk hygiene:
  overwrite each scenario's `data/<run>` in place).
- **This item is documentation-only at creation.** Implementation (execute 005
  Phase 2, then 006's chosen sweeps) is to be carried out by a separate agent on
  the 60-thread workstation *after* 005 Phase 2 completes.

## Acceptance

006 succeeds if there is a configuration under which the particle-wake CT history
**settles to a stable, non-oscillating (or nearly so) plateau** over multiple
post-withdrawal revolutions **AND** that plateau sits in **0.068–0.072** at
**≥4R truncation** with residual `magVinf ≈ 0` and bounded particle count. Report
the converged CT, its peak-to-peak ripple, and its distance from the reference. A
stable-but-off-magnitude or near-magnitude-but-oscillating result is a partial
result that names the remaining lever.

## Caveats

- Magnitude is only meaningful at hover (`magVinf ≈ 0`) and at loose truncation
  (≥4R) — never trust a near-0.072 reading at tight truncation (it is inflated by
  deleted far-wake downwash).
- Over-damping to stop the oscillation can suppress the very near-disk inflow that
  supplies the magnitude — track CT level and ripple together, never ripple alone.
- If no configuration achieves both halves, that negative result routes back to
  004's wake-mechanism ablations and to whether a fundamentally different wake
  treatment (e.g. the 002 prescribed-helical-wake direction) is required.

## 2026-07-07 — CT knob sweep result (executed 2026-07-03/04): relaxation is the dominant CT suppressor

Full execution of `plans/20260702_convergence.md` on the workstation (14 warm-start
scenarios + 1 extension; run table `data/rotor_hover_sweeps.md`, analysis
`data/ct_knob_sweep_notes.md`, machine-readable `data/ct_knob_sweep_summary.csv`).
Baseline: workstation-regenerated `data/rotor_hover_pressure_comparison` (40_40,
NT=36, RPM=6000, rlxf=0.3, 360 steps) — a different chaotic realization than the
local baseline (0.0663 still decaying at rev 10 vs local settled 0.06238), so all
knob effects were judged against the same-sweep control.

**Headline: Pedrizzetti relaxation at the stock setting (rlxf=0.3, λ=rlxf/dt=1080
s⁻¹) suppresses hover CT by ≈5.5e-3.** Control continuation plateaus at
CT_bernoulli 0.06403 ± 2.8e-4; `RELAX_FILTER_DOWNSTREAM_R=0.5` (relaxation only
≥0.5R downstream), extended to 18 revs so it actually settles, equilibrates at
**CT ≈ 0.0695** (per-rev means 0.06839 → 0.06928 → 0.06953 → 0.06948 over revs
14–18). That closes ~70% of the gap to experiment (~0.072) and lands in this
item's 0.068–0.072 acceptance band.

Supporting structure:
- rlxf halving ladder (fixed dt): CT(λ=1080)=0.06403, CT(540)=0.06513,
  CT(270)=0.06576, CT(135)=0.06638 — monotone, ALL stable (rlxf→0 instability did
  not appear down to 0.0375 over 4.3 revs).
- Filter-depth: 0.25R→0.06625, 0.5R→0.06727, 1.0R→0.06757 (4.3-rev reads).
- **Censored-equilibria caveat:** low-relaxation states need ~5–8 hover revs to
  equilibrate; every 4.3-rev continuation number above is a lower bound, the
  4-point ladder is non-geometric (increments 1.10, 0.63, 0.62 e-3), and its λ→0
  Richardson extrapolation (0.0666) is an underestimate. Only the 18-rev extension
  value (0.0695) is a settled number.

**Exonerated as CT suppressors** (all stable, |ΔCT| ≤ ~2σ = 5.7e-4):
surface→particle strain restoration (BODY_HESSIAN_TO_PARTICLES,
PANEL_WAKE_HESSIAN_TO_PARTICLES, and both combined — see 008 closure note),
KERNELOFFSET_TARGETS (5e-4, 2.5e-4), MERGE_PARTICLES=false.

**Infrastructure finding:** warm-start restarts were silently broken when the
metadata frame manifest is intact — `Das` TE shed vectors rotate with the body
(`rotate_Das!`) but are not persisted, so manifest-path restores left them at t=0
orientation, misplacing the first wake row and blowing up the Kutta condition
(Γ_TE +40% at the first continued step, CT diverging ~−0.03/step). Fixed
(UNCOMMITTED as of 2026-07-07) in `src/FLOWPanel_warmstart.jl` §2.5: always
reconstruct kinematics by replay; manifest is now a cross-check. Replay unit tests
89/89; smoke gate passed (warm-start control reproduced baseline tail to 2.3e-4).
The fix is REQUIRED for any warm-start work — do not revert.

**Next step (plan `plans/20260707_rotor_hover.md`):** full cold run with
RELAX_FILTER_DOWNSTREAM_R=0.5, 17 revs, targeting a stable 2-rev plateau near
0.072; sanctioned follow-ups include settled filter-depth points (1.0R, 2.0R),
filter+reduced-rlxf combos, and deeper truncation.
