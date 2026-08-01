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

## Results (2026-07-08): 0.5R-filter cold run — no plateau; stable limit cycle, cycle-mean CT = 0.0680

Cold run `rotor_hover_relaxfilter0p5_cold` (17 revs) + four warm-start extensions
(`_ext`…`_ext4`, +3 revs each, revs 17→29; all warm-starts on the replay-fix
warmstart path). All five runs stable: no NaN, per-rev std ≤ 4.5e-4. Full tables in
`data/rotor_hover_sweeps.md` ("Cold run" section).

**Headline: the 0.5R-filtered configuration does not converge to a plateau — it
settles into a slow, stable limit cycle of CT_bernoulli:**

- post-withdrawal dip 0.0660 (rev ~12) → peak **0.0695** (rev ~18) → trough
  **0.0666** (rev ~23) → second peak **0.0696** (rev ~27) → declining again at
  rev 29. Period ~9–9.5 revs, amplitude ±1.5e-3.
- **Cycle-mean over one full period (revs 18–27): CT = 0.0680** (per-rev means
  0.06920, 0.06823, 0.06710, 0.06698, 0.06665, 0.06686, 0.06800, 0.06910,
  0.06960; within-rev std ≤ 4.5e-4).
- The previously quoted "settled CT ≈ 0.0695" (18-rev extension read) was taken
  near the oscillation *peak* — a 2-rev window anywhere on this cycle biases the
  read by up to ±1.5e-3. **All low-relaxation hover CT numbers must be quoted as
  cycle-means over ≥1 period (~10 revs).** The strict 2-rev plateau gate can
  never fire for this config.
- Laplace cross-check (COLD run only, rev 16–17): matderiv 0.07165, lamb 0.06813
  — bracket Bernoulli; no red flag. Extension Laplace columns ignored
  (warm-start ∂u/∂t corruption).

**Deliverable line:** stable at CT = 0.0680 ± 0.0015 (limit cycle, revs 18–29);
gap to experiment (0.072) ≈ 4.0e-3.

The oscillation is plausibly a physical/numerical hover wake breathing mode
(periodic accumulation and shedding of the starting-vortex-like structure below
the rotor) that stock relaxation was previously damping along with the thrust.

**Follow-ups in flight** (warm-start perturbations from ext4 step 1043 = rev 29,
+5 revs each, judged against the 0.0680 cycle-mean): filter depth 1.0R (running
2026-07-08), then 2.0R, RELAX_RLXF ∈ {0.15, 0.075} @ 0.5R, TRUNCATION_DEPTH_R
∈ {6, 8}. Next-lever recommendation pending those results.

## Results (2026-07-09): follow-up sweep closed — depth is the only lever, saturates at 2.0R; best CT = 0.0686 (cycle-mean)

Six warm-start perturbations from the 0.5R-baseline rev-29 state (step 1043), +5
revs each (plus a +5-rev extension of the winner); all stable, no NaN. Tables in
`data/rotor_hover_sweeps.md` ("Follow-ups" section).

- **Filter depth is the only active lever** and its phase-matched 5-rev means
  were 0.0672 (0.5R) → 0.0683 (1.0R) → 0.0691 (2.0R) → 0.0691 (4.0R): the depth
  curve **saturates at ≈2.0R**. Deeper filters raise within-rev noise 2–3×
  (std up to 1.4e-3) but stay stable.
- **Null knobs** (each within +0.1e-3 of baseline, tight std): RELAX_RLXF=0.15
  in the damped ≥0.5R region (⇒ far-field relaxation *strength* irrelevant;
  0.075 skipped), TRUNCATION_DEPTH_R=6 (⇒ far-wake retention irrelevant;
  8 skipped).
- **Methodology correction:** the 2.0R extension (revs 34→39) shows the same
  ~9-rev limit cycle (amplitude ±1.7e-3), and its **true full-cycle mean is
  CT = 0.0686** (revs 29–39) — the phase-matched +1.9e-3 estimate was inflated
  because the knob change itself excites the cycle. Judge every knob on a
  ≥1-period cycle-mean, never a 5-rev window.

**Deliverable:** baseline (0.5R filter): stable at **CT = 0.0680 ± 0.0015**
(limit cycle, revs 18–27); best relaxation configuration (2.0R filter):
**CT = 0.0686 ± 0.0017** (revs 29–39); **gap to experiment 0.072 ≈ 3.4e-3**.
Relaxation-side knobs are exhausted — the remaining ~5% thrust deficit is not
recoverable by relaxation scheduling.

**Next-lever recommendation (needs user approval):** spatial/azimuthal
resolution — RHPC_MESH refinement beyond 40_40 and/or NT>36 (smaller dt) — since
the deficit persists with an essentially undamped, untruncated wake; secondarily,
a cold run at 2.0R to remove the warm-start-transient caveat on the 0.0686
number (~1 day walltime).

## Results (2026-07-09, later): replay post-processing of the 2.0R cycle — monitor spread brackets the gap; hub is not a factor

Post-processed the saved revs-29–39 chain (no re-simulation; `pnl.replay` +
saved per-panel pressures; only the new Kutta–Joukowski cross-check needed
N-body work). New `select` panel-mask option added to
`ForceMonitor`/`KuttaJoukowskiForce` for hub exclusion (r < 0.1R, 1220/7288
panels). Full table in `data/rotor_hover_sweeps.md` ("Replay post-processing");
plot + CSV in `data/rotor_hover_replay2p0_forces/`. Replayed CT reproduces the
original runs' CSVs to 1e-12 (loader/integration validated).

Cycle means, revs 30–39: **Laplace(Du/Dt) 0.0713 ± 0.0015**, Bernoulli
0.0685 ± 0.0013, Laplace(lamb) 0.0665 ± 0.0012, KJ rev-means 0.060–0.063.
Hub ΔCT ≤ +2e-4 on all pressure monitors (−1.7e-3 on KJ edges).

Two conclusions that update the convergence picture:

1. **The "gap to experiment" is monitor-dependent and the monitors bracket
   0.072.** The material-derivative Laplace pressure lands within 1e-3 of the
   reference; the cross-monitor spread (~4.8e-3) exceeds the Bernoulli-based
   deficit (3.4e-3). The remaining shortfall is therefore within
   pressure-recovery uncertainty, not clearly a wake/relaxation deficiency.
2. **Hub exclusion is a dead end** for closing the gap (≤ +0.35% of CT).

KJ diagnostic: instantaneous signal is unusable (fixed-azimuth ±0.1–0.2 1/rev
spikes from probes near the first wake row); per-rev means are tight and sit
~7e-3 below Bernoulli.

**Deliverable line:** settled 2.0R cycle, revs 30–39, cycle-mean CT ± std:
Laplace(Du/Dt) 0.0713 ± 0.0015 (gap −0.7e-3), Bernoulli 0.0685 ± 0.0013
(gap −3.5e-3), Laplace(lamb) 0.0665 ± 0.0012 (gap −5.5e-3), KJ ~0.061–0.063;
hub contribution ≤ +2e-4 (pressure) / −1.7e-3 (KJ).

## Results (2026-07-09, later still): spanwise loading percentiles for the settled 2.0R cycle

Replayed the same revs-30–39 window (steps 1080:1403, 324 samples) through
per-source (Laplace Du/Dt / Laplace lamb / Bernoulli) spanwise binning of the
saved panel pressures — zero N-body work, ~75 s. Script:
`examples/rotor_hover_spanwise_replay2p0.jl`; outputs (stats CSV + median/
IQR/min–max percentile plot per source, 35 bins/blade) in
`data/rotor_hover_replay2p0_forces/spanwise/`.

- Binned dT/dr integrates back to the forces-replay CT cycle-means at ratio
  1.0000 for all three sources (self-consistency of binner + integration).
- Blade 1 / blade 2 medians overlay indistinguishably — settled axisymmetric
  hover confirmed at the loading-distribution level, not just in total CT.
- All sources peak at r/R ≈ 0.76–0.78; the matderiv-vs-Bernoulli CT excess
  accumulates over the outboard half. Laplace variants have wide min–max
  scatter (and a small negative median dip) inboard of r/R ≈ 0.3 where bins
  are hub-influenced; Bernoulli is tight everywhere.

## Cross-reference: monitor reliability study (2026-07-09)

The three-monitor CT spread quoted throughout this document is resolved in
`BRAINSTORM/007_pressure_monitor_reliability.md`: steady Bernoulli (0.0685) is
the defensible cycle-mean estimate, cross-validated to 1.2e-4 by the ω-free
Laplace form; the matderiv excess (+2.9e-3) rides on normal-residual velocity
content at near-singular hub/tip panels, and the lamb deficit (−1.8e-3) is
entirely the quad-basis bound-κ injection. The −5% gap to experiment 0.072 is
physics/discretization, not monitor choice.

## Results (2026-07-27): production-mesh test — mesh refinement does not close the gap; the 2.0R filter triggers an instability there (unfiltered cycle-means 0.0615 / 0.0621)

Executed as Phase 2e of the DJI convergence study
(`plans/dji_convergence_20260722/`; all tables and run provenance in
`logs/dji_convergence_20260722/phase_02e_unsteady_ct.md`). This carries out the
2026-07-09 next-lever recommendation above — mesh refinement beyond 40_40 — on
the production mesh `examples/data/dji9443_20260725_45_185_capped_captess4.msh`
(36,752 panels vs 7,288 for 40_40), 5400 RPM, `nwakerows=1` particle wake,
truncation 4R, `BERNOULLI_ONLY`. Six cold HPC cases: two unfiltered baselines
completed 719 steps (20 revs); four 2.0R-filtered cases were killed early.

- **Mesh refinement did not raise CT.** Unfiltered 10-rev cycle-means (revs
  10–19): velocity **0.06148 ± 0.00152**, Green **0.06205 ± 0.00188**. Both sit
  *below* the 0.068–0.072 acceptance band and close to the 40_40 unfiltered
  control (0.06403, 2026-07-07 section). A 5× panel-count increase moved CT by
  less than the cycle scatter, so the resolution lever is spent in the **mesh**
  direction.
- **Stock relaxation does not plateau here either.** Revs 10→19 drift
  monotonically downward 0.0628 → 0.0596 (velocity), per-rev spread 2.7–2.8%,
  peak-to-peak 4.4–5.7%. The limit-cycle behaviour of the 2026-07-08 section is
  not confined to low-relaxation configurations.
- **Scope qualification on the 0.0686 deliverable.** All four filtered cases
  died before reaching hover, at steps 421, 420, 219 and 128 of 1007. The
  failure surfaces as an `OutOfMemoryError` thrown inside `merge_particles!`
  (`FLOWVPM_merging.jl:454`); ParaView inspection of the retained
  `p2e_OOM_sigF_das4p0` history confirms the wake had blown up, and job
  12894481 reached CT ≈ 90 before dying the same way. The 2026-07-09 numbers
  remain valid **on 40_40**; they are not reproducible at this mesh size as-is.
- **What that isolates — and what it does not.** The filter is the **trigger**:
  the unfiltered baselines ran the identical schedule, formulation and mesh to
  step 719 while filtered cases failed well short of it, and the *σ-stock*
  filtered cases — particle resolution identical to the healthy baselines —
  failed too, so neither run length nor σ refinement is needed to provoke it.
  **It does not follow that the filter is at fault.** Pedrizzetti relaxation is
  artificial dissipation — measured in the 2026-07-07 section suppressing hover
  CT by ≈5.5e-3 — so disabling it near the disk removes a damper that may have
  been masking a pre-existing fragility: the coarse wake did not need it, the
  finer wake does. The filter may be usable here once that fragility is
  repaired, and it is worth repairing, being the only lever this item found
  that moves CT toward the reference band.
- **Prime suspect for the fragility: wake resolution / overlap conditioning**
  (items 008/011/012). Measured particle `sigma/R` on this mesh is
  **0.022 / 0.145 / 0.379** (min/median/max) — the tip vortex is smeared over a
  third of the rotor radius. See the cross-reference below.
- **The experiment that would separate the two is untested:** an *unfiltered*
  run at refined σ. Every refined-σ case so far also carried the filter, so
  "the filter destabilises the wake" and "an under-resolved wake cannot
  tolerate reduced damping" remain indistinguishable on present evidence.
- **Truncation-depth null reproduced** independently on a new mesh and RPM: 4R
  vs 6R gives |ΔCT| ≤ 0.22% at matched revs 7–10 for both formulations,
  confirming the 2026-07-09 null.
- **First rotor measurement of `GreenReconstruction`**
  (`src/FLOWPanel_formulation.jl`): Green − velocity = **+0.93%** on the 10-rev
  cycle-mean (+0.34% on the last 5) — smaller than the ±2.5–3% cycle scatter, so
  **not resolved as an improvement**, and far short of the ~95%-of-a-10.8%-
  deficit recovery projected for the equivalent scheme A′ in `code_audit/log.md`
  (Task 1b, pitching-wing oracle).
- **`NT>36` is still the live half of the 2026-07-09 recommendation.** A
  cancelled NT=72 partial reached rev 7 at CT 0.0696 against NT=36's 0.0643
  (+8%), *opposite* in sign to 005's `nt_hi` (−6.6%). Still inside freestream
  withdrawal, so a lead rather than a result; partial histories retained for a
  warm start.

**Deliverable:** on the 36,752-panel production mesh at 5400 RPM the unfiltered
particle wake gives **CT = 0.0615 ± 0.0015 (velocity)** and **0.0621 ± 0.0019
(Green)** as 10-rev cycle-means — **~15% below experiment (0.072)** and below
this item's acceptance band. The 2.0R relaxation filter that supplied the best
40_40 result cannot presently be run at this mesh size, because the wake goes
unstable once near-disk relaxation is removed. Mesh refinement is therefore
**not** the missing lever, and the magnitude question is reopened.

**Next-lever recommendation:** treat the near-disk relaxation dependence as the
central question rather than a nuisance — establish whether reduced damping is
survivable at adequate wake resolution (item 012), since the CT gain lives on
that side. Concretely: (i) an unfiltered refined-σ run to separate trigger from
cause; (ii) if stable, re-apply the 2.0R filter on top of it; (iii)
independently close out `NT>36`, untested at this mesh and currently pointing
opposite to 005.

**Follow-up, resolved same day — the divergence reading is confirmed.** The four
filtered cases were resubmitted at 256/384 GB with run length matched to the
baselines (719 steps), as a direct test: if the failure were a memory ceiling
the extra headroom would carry them; if it were divergence it would not. Both
σ-fine cases **failed again at 384 GB** (12915036 after 2 h 41 m, 12915037 after
1 h 41 m), each with `ERROR: LoadError: OutOfMemoryError()` raised inside
`merge_particles!` — the same signature as before, now thrown by Julia itself
rather than as a Slurm OOM-kill because the allocation ceiling was no longer the
binding limit. **A 3× memory increase changed nothing**, which is what
divergence predicts and a memory ceiling does not.

**Evidence caveat for anyone re-reading these logs:** Julia's stdout is
block-buffered to file and a SIGKILL discards the buffer, so a killed run's log
is truncated tens of steps before the failure. A sane `max|CT|` over that prefix
is **not** evidence the run stayed bounded — an earlier reading of these same
jobs wrongly concluded "genuine memory ceiling, not divergence" on exactly that
basis. Magnitude tests are only valid on flushed (completed) logs.

## Cross-reference: particle resolution and overlap (2026-07-27)

The stability limit found above is squarely
`BRAINSTORM/012_robust_resolution_overlap_management.md` territory, and Phase 2e
supplies the first production-mesh measurements that item lacked. Particle
`sigma/R` in the settled wake is **0.022 / 0.145 / 0.379** (min/median/max), so
the tip-vortex core is a third of the rotor radius wide — concrete support for
012's warning that smearing tip vortices changes `C_T`. Mechanism note for
anyone tuning it: with `nwakerows=1` the trailing shed segment spans one
timestep of trailing-edge travel, so `sigma = dist * OVERLAP / P_PER_STEP` is
set by the per-step travel and is **independent of the `Das` first-row offset**;
`MERGE_R_FACTOR` must be scaled with `sigma` or merging undoes any refinement.

There is also a **floor** on the overlap ratio that qualifies 012's "bound
`σ/h`" direction. Refining to `σ_med ≈ 0.024R` by dropping nominal `OVERLAP`
3 → 2 at `P_PER_STEP=8` diverged almost immediately (CT ≈ 90 by step 58), while
`OVERLAP=4` at the same particle spacing survived an order of magnitude longer.
Reducing overlap is therefore not monotonically beneficial: below roughly the
stock ratio the particle field stops representing a smooth vorticity field.
011's "the lever is overlap reduction" should be read as reducing the **true
local** `σ/h ≈ 4.2` toward its nominal value, not as lowering the nominal
setting.

## Staged studies (2026-07-27): next runs, ranked

Staged at Ryan's direction after the production-mesh results above. All assume
the **unfiltered** configuration as the carrier, since it is the only one that
currently completes at this mesh size; the filter is re-introduced only where a
study explicitly tests it. Baseline cost is ~10 h per 719-step case at 64 GB,
so five cases fit the six-job cap comfortably. Driver knobs referenced below are
already wired in `examples/rotor_hover_pressure_comparison.jl` and case-tagged
in `examples/run_dji9443_hover_ct_hpc.slurm.sh`.

- **(A) `Das` length convergence — the requested study, and the largest
  unquantified discretization uncertainty.** Sweep `DAS_ETA_KINEMATIC` ∈
  {0.2, 0.5, 1.0, 2.0, 4.0} unfiltered at NT=36, stock σ, truncation 4R, 719
  steps. `Das = max(η·dt·|V_te|, min_disp)` is the offset of the first shed wake
  row from the trailing edge; the driver has always hard-coded **η = 0.2**,
  which is short — η = 1.0 is the geometrically natural "shed where the TE was
  one step ago". Transient revs 0–4 already show η=4 sitting **12–25% above**
  η=0.2, larger than the relaxation filter's +7% and larger than the oracle-wing
  analogue (+13% over a 4× `das` range, `code_audit/log.md`). **This sweep is
  clean — now measured, not assumed**: with `nwakerows=1` the shed segment spans
  one step's TE travel, so particle σ is set by `dist·OVERLAP/P_PER_STEP` and is
  **independent of `Das`** — η moves one thing only. Checked on the live η=1.0
  and η=4.0 jobs (step 9): a 4× change in η moves median σ/R by +1.2% (0.0685 →
  0.0693) and the outboard 0.6–0.95R median by −0.2%, so the tangential-offset
  curvature residual that could have confounded the sweep (naively +35% at η=4)
  is absent. Detail in `logs/dji_convergence_20260722/phase_02e_unsteady_ct.md`,
  2026-07-27 23:40 entry. *Acceptance:* CT-vs-η flattens, naming a converged η
  and an uncertainty band. *Discipline:* η is a discretization parameter, not a
  tuning knob — the deliverable is where CT stops depending on it, **not** the η
  that lands nearest 0.072. If CT never flattens, that is itself the finding and
  it would make first-wake-row placement the dominant model uncertainty.

- **(B) Filter-depth survivability at production mesh — the cheapest direct
  attack on the reopened magnitude question.** The 2.0R filter is unstable here,
  but on 40_40 the depth curve was nearly flat: 0.0672 (0.5R) → 0.0683 (1.0R) →
  0.0691 (2.0R). A **shallower** filter leaves a much smaller unrelaxed
  near-disk band and may be survivable while still buying most of the CT. Run
  `RELAX_FILTER_DOWNSTREAM_R` ∈ {0.5, 1.0} unfiltered-carrier otherwise stock.
  *Acceptance:* does any filter depth both complete 719 steps and raise the
  10-rev cycle-mean above 0.0615? A stable 0.5R result near 0.068 would restore
  this item's magnitude claim on the production mesh.

- **(C) Unfiltered refined-σ — separates trigger from root cause.** Every
  refined-σ case so far also carried the filter, so "the filter destabilises the
  wake" and "an under-resolved wake cannot tolerate reduced damping" are still
  indistinguishable. Run `OVERLAP=4, P_PER_STEP=8` (σ_med ≈ 0.048R) with **no**
  filter. *Acceptance:* if it completes, the filter is implicated directly and
  (B) becomes the priority; if it also fails, resolution is the root cause and
  the work belongs in item 012.

- **(D) Temporal convergence `NT` ∈ {18, 36, 72} — with a confound to remove
  first.** 005 measured NT 18→36 at **−6.6%**; our cancelled NT=72 partial
  pointed **+8%** the other way. Before spending jobs, note that refining NT at
  fixed η changes **three things at once**: `dt`, the `Das` length (∝ `dt`), and
  particle σ (∝ per-step TE travel ∝ 1/NT). A naive NT ladder is therefore *not*
  a temporal-convergence study and may explain the sign conflict. *Do it
  properly:* hold `Das` and σ fixed in physical units while refining `dt` — i.e.
  scale `DAS_ETA_KINEMATIC` ∝ NT and pin σ via `PARTICLE_SHEDDING=sigma_overlap`
  — then the ladder isolates the timestep. Gate this behind (A), which supplies
  the converged η.

- **(E) `shed_with_induced_velocity` — a free null-test.** The hover driver sets
  it `false` while `examples/rotor_panel_wake_study.jl` sets it `true`, an
  untested inconsistency with no measured CT sensitivity anywhere in the repo.
  Shedding at the local induced velocity is arguably the more defensible choice.
  Env `SHED_WITH_INDUCED_VELOCITY` is wired. One job, unfiltered, stock.

- **(F) Coupling ablations — instability diagnostics, not production configs.**
  Ryan suggested disabling the surface-on-wake velocity gradient. **That is
  already off**: `BODY_HESSIAN_TO_PARTICLES` defaults `false`
  (`examples/rotor_hover_pressure_comparison.jl`), and `simulate!` requests no
  `velocity_gradient` for a `ParticleField` target when it is false; the
  attached-wake-row gradient is off too
  (`PANEL_WAKE_HESSIAN_TO_PARTICLES` false via `WAKEROW_NO_HESSIAN_TO_PARTICLES`).
  Every run in this study already had both disabled, so neither can explain or
  cure the blow-up. The **only live velocity gradient is particle-self**
  (`PARTICLE_HESSIAN_SELF=true`), which is what drives **vortex stretching** —
  the classic VPM blow-up mechanism. Two ablations, both run *with* the 2.0R
  filter so there is an instability to diagnose: `PARTICLE_HESSIAN_SELF=false`
  (no stretching) and `BODY_ON_WAKE=false` (remove body-induced *velocity* on
  the wake). *Acceptance:* if disabling stretching alone keeps the filtered case
  alive to 719 steps, the blow-up is stretching amplification in the unrelaxed
  near-disk band, which points squarely at resolution (012) rather than at the
  filter. Neither ablation is physically admissible for a production CT — they
  are attribution tools only.

**Suggested order:** A and B in parallel (both attack live questions and neither
blocks the other), C alongside them if slots allow since it settles the
filter-vs-resolution attribution; F only if the filtered cases keep failing and
the mechanism is still unclear; then D once A has fixed η; E opportunistically.
**Deliberately not staged:** further relaxation-*strength* tuning, truncation
depth, kernel offsets, merge settings, viscous/SFS knobs — all measured null in
the 2026-07-09 section or reconfirmed null on the production mesh.

## Interim (2026-07-28): study A's ladder has flattened — η=0.2 was under-converged by ~16%

> **SUPERSEDED** by the settled-results section below. This section was written
> from the revs 0–4 transient while the jobs were still running; its "do not
> quote 0.071" caution has since been satisfied by completed 10-rev cycle-means.
> Retained because its transient prediction (+15–18%) is what the settled result
> confirmed, and because the arc-artifact argument still stands.

**Lead, not a result.** Full tables and provenance in
`logs/dji_convergence_20260722/phase_02e_unsteady_ct.md` (2026-07-28 entries).

Matched per-rev CT (revs 0–4, velocity, production mesh, NT=36, 4R, stock σ),
η=0.2 column being the completed unfiltered baseline 12894164:

| η | rev 2 | rev 3 | rev 4 |
|---:|---:|---:|---:|
| 0.2 | 0.04446 | 0.04863 | 0.06025 |
| 1.0 | 0.05119 | 0.05663 | 0.07034 |
| 4.0 | 0.05102 | 0.05661 | 0.07054 |

- **η=1.0 and η=4.0 agree to <0.5%** while **0.2→1.0 is +15–18%**. Study A's
  acceptance criterion (CT-vs-η flattens) is met with **η=1.0 converged**; the
  driver's hard-coded **η=0.2 is the outlier**, and it underlies *every* prior
  006 number.
- **Not a geometric artifact.** The legacy straight-vector `Das` lands the first
  wake row at radius `r√(1+θ²)` — 0.015R too far out at η=1.0 but **0.22R at
  η=4.0**. If that drove CT, η=4 would exceed η=1; it does not. The +16% is a
  genuine first-wake-row-placement effect.
- **Why it matters:** applied to the settled unfiltered baseline
  0.06148 ± 0.00152, a ~16% elevation lands near **0.071** — inside this item's
  0.068–0.072 band, and **without** the near-disk relaxation filter that
  destabilises the production mesh. That would make first-wake-row placement,
  not wake resolution, the missing magnitude lever, and would reframe the
  2026-07-27 "magnitude question REOPENED" conclusion.
- **Caveat that must not be dropped:** revs 0–4 sit inside the staged freestream
  ramp/withdrawal (`SETTLE_REVS=12`), so all three cases are transient. The
  deliverable is the 10-rev cycle-mean at revs 10–19 when 12918582/12918583
  finish. Do not quote 0.071.
- **Next rung:** η=0.5 (not η=2.0) — the entire jump is between 0.2 and 1.0, so
  η=0.5 locates the knee; η≥1 is already flat.

**Code change (unreleased, not on the cluster):** `Das` now follows the swept arc
rather than the tangent, behind `DAS_KINEMATIC_ARC` (default on) — exact rigid
backward path, 25-assertion testset, regressions green. At the production η=0.2
the correction is 6e-4 R, so no existing 006 number changes. The in-flight
straight-vector ladder runs to completion as the A/B reference.

## SETTLED RESULTS (2026-07-28): magnitude half MET at η=1.0 (CT = 0.0713, in band); stability half NOT met

Three Batch-3 jobs completed cleanly (`COMPLETED` exit 0:0; their
`P2E_DIVERGED reason=nan_inf` flags are a false positive that the known-good
baselines also carry). 10-rev cycle-means over revs 10–19, computed from each
case's `monitors/*_monitor02_force_system1.csv`; the extraction reproduces the
logged baseline exactly (0.06148 ± 0.00152), validating the method. Production
mesh (36,752 panels), 5400 RPM, velocity, NT=36, truncation 4R:

| case | `Das` η | filter | CT | std | rel | per-rev p-p |
|---|---:|---:|---:|---:|---:|---:|
| baseline `p2e_vel_nt36_d4` | 0.2 | off | 0.06148 | 0.00152 | 2.47% | 7.79% |
| **`p2e_das1p0`** | **1.0** | **off** | **0.07133** | **0.00159** | **2.23%** | **7.15%** |
| `p2e_das4p0` | 4.0 | off | 0.07190 | 0.00101 | 1.40% | 4.19% |
| `p2e_filt0p5` | 0.2 | 0.5R | 0.06616 | 0.00109 | **1.64%** | **4.34%** |
| `p2e_filt1p0` | 0.2 | 1.0R | 0.06578 | 0.00170 | 2.59% | 7.96% |

**Study A — ANSWERED.** `Das` η is the magnitude lever this item has been
hunting. η=1.0 gives **CT = 0.0713 ± 0.0016**, **+16.0%** over the η=0.2
baseline and **inside the 0.068–0.072 acceptance band**, at 4R truncation,
**unfiltered**, on the production mesh. η=1.0 → 4.0 is only +0.80%, below either
run's own cycle scatter, so the ladder is flat for η ≥ 1 in the settled state.
Converged η = 1.0; converged CT ≈ 0.0713–0.0719; uncertainty ±0.0016 (cycle
scatter, which dominates the η sensitivity above η=1).

**This reverses the 2026-07-27 verdict.** The magnitude gap was **not** wake
resolution and **not** mesh refinement. It was the first-wake-row offset — a
discretization parameter hard-coded at η=0.2 from the beginning and never
varied, inherited by every prior number in this item.

**Study B — ANSWERED, positively.** Both shallow filters completed 719 steps and
raised CT above 0.0615 (0.0662 at 0.5R, 0.0658 at 1.0R). The 2.0R filter's
instability at this mesh is therefore a **depth** effect, not a filter effect.

**006 acceptance status: PARTIAL — magnitude met, stability not.** η=1.0 lands
in band but its per-rev spread is 2.23% and p-p 7.15%, against tolerances of
0.5% and 2%. This is a limit cycle at the right level — precisely the
"near-magnitude-but-oscillating" partial result this item's acceptance
anticipates. The named remaining lever is the **combination**: `filt0p5` has the
lowest ripple measured at this mesh (1.64% / 4.34%) *and* raises CT, so the two
levers are independent and both point the right way. `p2e_das1p0_filt0p5`
(η=1.0 + 0.5R filter, job 12921072) is in flight. Caveat: the filter added
+0.0047 at η=0.2, so if the effects are additive the combination lands near
0.076 — above the band — which would make the filter a magnitude knob to back
off rather than a free stabiliser.

### Do not read 0.0713 as settled physics yet — see BRAINSTORM/014

The **mechanism for the +16% is unexplained**, and three candidate explanations
have been proposed and refuted (near-wake pushed out of influence; particles
seeded inside the attached panel's span; shed cores straddling the TE). The
question "on what principle should η be chosen?" is now its own open item,
**`BRAINSTORM/014_first_wake_row_offset_selection.md`**. Its verdict governs how
this section should be read:

- **If η is determinable on principle**, 0.0713 is a prediction and this item's
  magnitude claim stands.
- **If η is irreducibly a model parameter**, 0.0713 is *tuned*, and ±16% must be
  carried as a model-form uncertainty on every CT reported here.

014 also records a **confound discovered after these runs**: the hard-coded
`min_displacement = 0.01R` floor sets `Das` — not η — over **~68% of the
shedding span at η=0.2** versus 3.6% at η=1.0, so this ladder was not a clean
single-parameter sweep. 014's recommended next run (η=0.2 with the floor
removed) can reframe the attribution above and should be run before further η
rungs.

Also still open and in flight: **NT=72 at η=1.0** (job 12921559) tests whether
0.0713 survives timestep refinement, and **`p2e_das1p0_arc`** (job 12921071)
re-confirms it under the corrected arc `Das` construction — all results above
used the legacy tangent construction.

## Update (2026-07-28, evening): combination completed — overshoots; stability half still unmet by every configuration

`p2e_das1p0_filt0p5` (η=1.0 + 0.5R filter) completed: **CT = 0.07625 ± 0.00140**
(rel 1.84%, p-p 5.38%) — **above** the 0.068–0.072 band. The two levers are
**cleanly additive** (predicted 0.07601 from η=1.0's 0.07133 plus the filter's
+0.00468; measured 0.07625, within 0.3%), confirming the caveat recorded when
the run was launched: the 0.5R filter is a **magnitude knob**, not a free
stabiliser.

| case | η | filter | CT | rel | p-p | in band | stable |
|---|---:|---:|---:|---:|---:|:--:|:--:|
| baseline | 0.2 | off | 0.06148 | 2.47% | 7.79% | no | no |
| — | 0.5 | off | 0.06942 | 2.66% | 7.30% | no | no |
| — | 1.0 | off | 0.07133 | 2.23% | 7.15% | **yes** | no |
| — | 4.0 | off | 0.07190 | **1.40%** | **4.19%** | **yes** | no |
| — | 0.2 | 0.5R | 0.06616 | 1.64% | 4.34% | no | no |
| combination | 1.0 | 0.5R | 0.07625 | 1.84% | 5.38% | no (over) | no |

**006 status unchanged: PARTIAL — magnitude met, stability NOT.** No
configuration tested meets the gates (rel ≤ 0.5%, p-p ≤ 2%); the best is η=4.0
at 1.40% / 4.19%, still ~2.8× and ~2.1× over. The combination was this item's
named remaining lever and it did not deliver the stability half.

**Do not tune to the target.** Since the levers superpose, η≈0.35 + 0.5R filter
would land ~0.070. That is fitting two knobs to a known answer and must not be
reported as a converged result — see BRAINSTORM/014.

**Where the stability half goes next:** the lowest ripple measured at this mesh
comes from **η=4.0 unfiltered** (1.40% / 4.19%), not from damping — so the
remaining lever may be on the η/wake-representation side rather than the damper
side. That is now entangled with 014's verdict, whose evening results show CT
collapsing as `Das → 0` regardless of what is downstream (floor-off and
`nwakerows=2` both), pointing at Kutta enforcement rather than wake extent.

## Ruling 2026-07-28 (Ryan): relaxation policy for production configs

**Rule: production rotor configs use stock relaxation (corrected-Pedrizzetti,
default `rlxf`) with NO near-disk filter.** Basis: at η=1.0 the 0.5R filter
adds +0.005 CT (0.0713 → 0.0763, overshooting band and experiment) while
delivering only marginal ripple improvement (per-rev std 2.2%→1.8%), and the
2.0R filter destabilizes the production mesh outright.

The filter and relaxation-strength knobs are retained **as diagnostics only**:
relaxation is artificial dissipation whose ~−0.005 CT suppression is monotone
in `rlxf` and not converged (2026-07-07 ladder), so it is carried as an
explicit term in the CT error budget alongside the ~±1% `Das` log-law term
(BRAINSTORM/014). Verifying convergence in that term later (rlxf→0 / filtered
limits) is expected to require adequate near-disk wake resolution first
(BRAINSTORM/012) — removing damping at 2.0R currently diverges.

## HPC JOB STATUS (2026-07-29) — supersedes every earlier "in flight" statement in this file

**All earlier "in flight" / "pending" claims above are STALE.** This block is
authoritative. Everything relevant to 006 has completed except one job.

| Job | Case | State | 10-rev cycle-mean CT (revs 10–19) |
|---|---|---|---|
| 12894164 | `p2e_vel_nt36_d4` (η=0.2 baseline) | COMPLETED | 0.06148 ± 0.00152 |
| 12920967 | `p2e_das0p5` | COMPLETED | 0.06942 ± 0.00185 |
| 12918582 | `p2e_das1p0` | COMPLETED | **0.07133 ± 0.00159** |
| 12918583 | `p2e_das4p0` | COMPLETED | 0.07190 ± 0.00101 |
| 12921071 | `p2e_das1p0_arc` | COMPLETED | 0.07108 ± 0.00184 (arc immaterial) |
| 12918580 | `p2e_filt0p5` | COMPLETED | 0.06616 ± 0.00109 |
| 12918581 | `p2e_filt1p0` | COMPLETED | 0.06578 ± 0.00170 |
| 12921072 | `p2e_das1p0_filt0p5` | COMPLETED | 0.07625 ± 0.00140 (**over** the band) |
| 12925945 | `p2e_das0p2_nofloor` | COMPLETED | 0.05215 ± 0.00078 |
| 12925947 | `p2e_nrows2_dassmall` | COMPLETED | 0.03362 ± 0.00046 |
| 12927923 | `p2e_nrows4_das1p0` | COMPLETED | 0.07431 ± 0.00096 (**over** the band) |
| 12927924 | `p2e_das1p0_refresh` | **CANCELLED — INVALID** | `Backslash` factorized `G` with `Das = 0`; Γ collapsed 21.7× |

### ONE JOB STILL RUNNING — check this after a context clear

**`12921559` — `p2e_nt72_das1p0`** (NT=72, η=1.0, `rlxf` 0.15, unfiltered, 4R,
48 h wall). At 21 h it was at step ~1174 of 1438. **This is the only outstanding
run that can change 006's status.**

It asks: **is CT = 0.0713 dt-converged?** Caveat: at fixed η, halving `dt` also
halves `Das` *and* σ, so three things refine together (`Das`/σ = 0.267 is
invariant). Holding ⇒ strong simultaneous-convergence result. Moving ⇒ the run
alone cannot attribute the change.

Harvest (steps/rev = **72** for this case, not 36):
```bash
ssh orc 'bash -lc "cd /home/rander39/projects/FLOWPanel.jl && \
  awk -F, \"NR>1{r=int(\\\$1/72); if(r>=10&&r<=19){s[r]+=-\\\$3;n[r]++}} \
  END{c=0;for(i=10;i<=19;i++) if(n[i]>0){m[i]=s[i]/n[i];tot+=m[i];c++}; mu=tot/c; \
  for(i=10;i<=19;i++) if(n[i]>0){d=m[i]-mu;v+=d*d}; \
  printf \\\"CT=%.5f std=%.5f\\n\\\", mu, sqrt(v/(c-1))}\" \
  data/p2e_nt72_das1p0/monitors/p2e_nt72_das1p0_monitor02_force_system1.csv"'
```
Check the terminal state with `sacct -j 12921559` first — **a running job's prefix
is not evidence**; only terminal states count. Note `scripts/p2e_status.sh` flags
completed runs as `P2E_DIVERGED reason=nan_inf` — that is a **false positive** on
the all-NaN `CT_kj` column; only `reason=magnitude` has ever indicated a real
blow-up.

**006's verdict does not change while this runs:** PARTIAL — magnitude met
(0.0713 in band), stability NOT met by any configuration (best `das4p0` at
1.40% / 4.19% vs gates 0.5% / 2%).
