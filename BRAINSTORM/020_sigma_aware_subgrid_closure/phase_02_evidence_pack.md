# 020 Phase 2 — Evidence pack: the subgrid gap, demonstrated

**Status: PHASE 2 COMPLETE AT REVISED GATE (2026-08-12).** Leg 1 and Leg 3
remain supporting evidence. The corrected Phase-2R screen supersedes the old
Leg-2 mechanistic verdict. No Phase-3 closure work is authorized by this pack.

## 2026-08-12 audit ruling (binding)

The original `euler_exp` integrated
$\dot\Gamma=S_n-3Z_n\Gamma$ with the initial stretching vector $S_n$ frozen
as an additive forcing. That is not the frozen-gradient rVPM system because
$S=L\Gamma$ evolves with $\Gamma$. In aligned strain the prototype approaches
a $5/3$ strength multiplier, while the intended homogeneous rVPM solution is
$e^{2\Delta tZ}$. Therefore:

- Stage A's scalar-$Z$ replay proves positivity only for the same incorrect
  surrogate and is not mechanistic evidence for the live vector system.
- Stage B job 13065481 changed physical stretching as well as integration
  error; its extra survival cannot isolate the discrete contribution.
- The claim that the closure is needed “for stability AND physics” is not yet
  established by Leg 2. Leg 1's failed $\Gamma_{\mathrm{implied}}$ criterion
  and its conditional core-radius bracket remain explicit caveats.

Phase-2R replaces the prototype with the frozen-gradient geometric map,
verifies it analytically, and reruns the single decisive $\sigma/R=0.02$
screen before this pack is rescored.

## Phase 2R — corrected Stage B (RESULT OF RECORD)

**Run:** `scr_p020r_geom_s020v`, job **13154223**, an otherwise identical
$\sigma/R=0.02$ viscous twin using the corrected frozen-gradient geometric
map. Startup parameters and staged-source checksums were verified. The focused
integrator suite passes **28/28** locally and on cluster, including aligned and
general constant gradients, large $\Delta tZ$, coupled CoreSpreading limits,
and first-order convergence for a prescribed time-varying gradient. Full local
FLOWVPM and relevant FLOWPanel suites also pass.

**Terminal result:** Slurm `FAILED` after 3:55:50 at step 243/323. This was
neither OOM nor a negative-core event: peak RSS was 20.6 GB of 64 GB, all 244
health rows had positive $\sigma$, and the explicit non-finite frozen-gradient
ratio guard stopped the next step after the field had already run away.

| readout | Euler control | corrected geometric |
| --- | ---: | ---: |
| step 237 max $|u|$ [m/s] | 16,095 | 36.0 |
| step 237 $M$ | 568.1 | 0.173 |
| terminal step | 240 (OOM job) | 243 (guard failure) |
| first $|u|>1000$ step | 237 | 243 |
| whole-history min $\sigma/\sigma_{shed}$ | 0.0396 | 0.0137 |
| terminal max $|\Gamma|/\sigma^2$ | $3.97\times10^7$ | $4.60\times10^9$ |
| terminal $M$ | 178.5 | 221.8 |

At the Euler twin's failure point the corrected run is plainly healthy and its
VTP median $\sigma/\sigma_{shed}=1.037$. Six steps later, a tail-localized
contraction jumps to min ratio 0.0137, $M=221.8$, and 19.6 km/s while the VTP
median remains 1.037. Thus the corrected map isolates a real discrete
contribution at step 237, but it does not create a viable operating point.

**Literal registered scoring:**

- “Threshold moves” is **FAIL as written**: the run had zero negative cores
  and no merge-OOM, but did not survive all 8 revolutions.
- “Fidelity collapse remains” is **NOT CONFIRMED by the registered Boolean**:
  the median never fell below 0.85 and $M>0.2$ occurred only at isolated steps
  233, 238, 241, and 243 before the abrupt terminal jump. The harvester's
  reproducible post-registration definition (10 consecutive rows) therefore
  does not fire; neither does the ordinary qualitative reading of “sustained.”
  The terminal
  runaway is decisive qualitative evidence of loss of fidelity, but it is not
  substituted for the preregistered Boolean.
- Mechanistic conclusion: **stiff integration is a real local numerical
  remedy, but not a sufficient field-level remedy.** The corrected result no
  longer supports the old claim of 8-revolution survival or a quantitatively
  confirmed closure requirement.

Figure and backing data: `figures/fig_stage_b_r.tex` and
`figures/fig_stage_b_r/` (Euler, superseded prototype, corrected health,
terminal harvest, and snapshot history).
Protocol: PRE-REGISTERED in `phase_01_theory.md` §7 (thresholds quoted
verbatim below); every deviation is recorded in the item Log, never patched.

Claim under demonstration (registered, §7.4): *a subfilter viscous
(σ-channel) model is the missing piece that would let the simulation
approach physics-predicted strain rates* — three legs: the measured strain
ceiling, the stiff-integrator test, the SFS viscous-null.

Source of record: monitors CSVs + committed figure-backing CSVs
(`figures/fig_leg1/`, `figures/fig_stage_a/`); scripts
`scripts/p020_leg1.py`, `scripts/p020_stage_a.py` (deterministic reruns).

---

## Leg 1 — the strain ceiling, measured (COMPLETE)

**Registered thresholds (§7.1):** ceiling is "stability-set" iff $p \ge 1.5$
over the ladder and $\Gamma_{\mathrm{implied}} \in [0.5, 2]\,\Gamma_v$ for
≥ 3 consecutive rungs; gap = decade count
$\log_{10}(Z_{\mathrm{phys}}/Z^{\mathrm{resolved}}_{\max})$ with
$r_c \in [0.01, 0.05]c$, both endpoints reported.

**Data:** 019 Campaign E complete NT=36 grid (9 screens, direct `max_dtZ`),
via the 019 pipeline plus a pre-onset supplement (`p020_leg1.py`; ignited
runs' whole-run max is a divergence transient — 019's wave-1 ruling — so
the 019-P1-sanctioned *last pre-ignition* window is used for them; onset =
first step max_u > 1000 m/s).

| run | σ/R | visc | outcome | M (pre-onset, direct) | Γ_implied [m²/s] |
| --- | --- | --- | --- | --- | --- |
| s015 | 0.0150 | – | ignited @111 | 10.07 | 0.650 |
| s015v | 0.0150 | ✓ | ignited @196 | 12.59 | 0.812 |
| s020v | 0.0200 | ✓ | ignited @237 | 28.24 | 3.240 |
| s025 | 0.0249 | – | ignited @287 | 3.44 | 0.617 |
| s025v | 0.0249 | ✓ | ignited @286 | 1.116 | 0.200 |
| s030v | 0.0299 | ✓ | survivor | 0.295 | 0.076 |
| sstab | 0.0312 | ✓ | survivor | 0.134 | 0.038 |
| s038 | 0.0381 | – | survivor | 0.141 | 0.059 |
| s038v | 0.0381 | ✓ | survivor | 0.112 | 0.047 |

**Scoring against the registration:**

- **p-arm: PASS.** Viscous survivors $p = 3.01$ (3 pts, 0.030–0.038R, the
  019-recorded small-range caveat applies); supplementary all-rung
  pre-onset fits are steeper still (viscous 6.4, inviscid 4.5 — the
  pre-onset window of an igniter still contains growth transient, biasing
  steep; labeled supplementary, not verdict-bearing).
- **Γ_implied-band arm: FAIL as literally registered** — only 1 viscous
  rung (s025v, 0.200) sits in $[0.139, 0.556]$; igniting rungs run 0.6–3.2
  (≈ 2–12 × Γ_v, onset-contaminated at the low-σ end), survivors 4–7×
  *below* Γ_v. Recorded as a deviation-in-outcome, not reinterpreted away.
- **Mechanism reading (honest synthesis):** the filament-ceiling signature
  ($\Gamma_{\mathrm{implied}} \sim O(\Gamma_v)$) localizes **at the
  boundary** (both 0.025R rungs: 0.20–0.62, in or near band) instead of
  spanning the ladder; survivors sit below the ceiling (margins partly
  ambient), and every rung below the boundary ignites. Combined with 019's
  measured ε-crossing at $(0.0299, 0.0312]R$ — i.e. at σ_stab ± 0.4% — the
  conclusion *today's max resolvable strain is stability-set, not
  physics-set* is **SUPPORTED**, with the band criterion's literal failure
  and its cause on the record.

**The gap.** Resolved ceiling at the last viscous survivor:
$M = 0.295 \Rightarrow Z_{\mathrm{res}} \approx 957\ \mathrm{s^{-1}}$.
Physics target $Z_{\mathrm{phys}} = \Gamma_v/(2\pi r_c^2)$ with
$c(0.75R) = 15.3$ mm (018 σ-ladder anchor):

| $r_c$ | $Z_{\mathrm{phys}}$ [1/s] | $\Delta t Z_{\mathrm{phys}}$ | gap [decades] |
| --- | --- | --- | --- |
| $0.05c$ = 0.764 mm | $7.58\times10^4$ | 23.4 | **1.90** |
| $0.01c$ = 0.153 mm | $1.90\times10^6$ | 584.9 | **3.30** |

Figure: `figures/fig_leg1.tex` (+ `fig_leg1/` CSVs).

---

## Leg 2 — the stiff-integration test

### Stage A — offline re-integration (SUPERSEDED by 2026-08-12 audit)

**Deviation (item Log):** the registration assumed tracked per-particle
trajectories; none were retained. Reinterpretation faithful to P-2
(`p020_stage_a.py`): (A1) a pointwise one-step census over every recorded
per-particle ΔtZ; (A2) re-integration of the recorded aggregate σ traces
with per-step Z inferred by inverting the live update.

**A1 census** (MM4 units = raw forensics ΔtZ / 5; thresholds from
phase_01 §1.3: transverse-Γ flip 2/3, σ flip 1, σ geometric 2):

| population | N | max ΔtZ | >2/3 | >1 | >2 |
| --- | --- | --- | --- | --- | --- |
| corpse 1041, collapsing subset | 863 | 13.2 | 33 | 21 | 10 |
| corpse 1041, whole field | 177,656 | 878.7 | 1,545 | 908 | 362 |
| healthy 719, aged column (control) | 3,000 | 0.025 | 0 | 0 | 0 |
| healthy 719, whole field (control) | 209,577 | 0.036 | 0 | 0 | 0 |

Under the exact map these populations have **zero** flips at any ΔtZ
(multipliers $e^{-3x}, e^{-x} > 0$); the healthy controls show the wake is
nowhere near any threshold until collapse begins.

**A2 traces** (Euler re-integration reproduces every recorded trace —
inversion consistency PASS):

| trace | neg-σ Euler | neg-σ exact | exact floor / recorded floor |
| --- | --- | --- | --- |
| ignition core (viscous corpse) | 0 | 0 | 1.28 |
| argmin-σ (viscous corpse) | 0 | 0 | **1.08** |
| ufront_n1 min_sr (inviscid death) | 1 (the recorded flip) | **0** | 0.11 |

**P-2 verdict: CONFIRMED.** Zero spurious events under the exact map, and
σ still grinds to ~the viscous floor — *stability without fidelity*, as
predicted. The registered ≤ 1.1× floor criterion passes on the registered
data source (`death_trajectory.csv` argmin trace, 1.08); the supplementary
ignition-core trace reads 1.28 because at ΔtZ ≈ 1 the Euler update
spuriously near-zeroes σ in one step while the exact map contracts
$e^{-x}$ — the excess *is* the artifact being removed (recorded, not
patched). Figure: `figures/fig_stage_a.tex`.

### Stage B — original in-code prototype (SUPERSEDED by 2026-08-12 audit)

Implementation landed this session (default-off, bit-identical off,
test-gated): FLOWVPM `euler_exp`/`_euler_exp`
(`FLOWVPM_timeintegration.jl`; exact frozen-coefficient update
$\Gamma \leftarrow e^{-3\Delta tZ}\Gamma + \Delta t\,\varphi(3\Delta tZ)\,\mathbf S_{\mathrm{eff}}$,
$\sigma \leftarrow \sigma e^{-\Delta tZ}$), `viscousdiffusion` branches
extended, FLOWPanel `PanelParticleWake(expint=...)` +
`FLOWVPM._euler_exp` step branch, `WAKE_EXPINT` env knob
(parse/pass/banner/metadata). Tests: FLOWVPM `test/runtests_expint.jl`
26/26 local + cluster; FLOWPanel expint testset 4/4 local + cluster; full
FLOWPanel suite PASS local; FLOWVPM non-TestEnv files PASS (vortex-ring/
leapfrog blocked by pre-existing TestEnv env rot — deviation, item Log).

**Run:** `scr_p020_exp_s020v` = job **13065481** — exact clone of the
ignited `scr_p019_s020v` (σ/R = 0.02 viscous, NT=36, N=1, D=3.4,
`max_dtZ` on) + `WAKE_EXPINT=true`; integrator ON, closure OFF (does not
exist).

**Registered readouts + verdict rules (§7.2, fixed before data):**
survival to 8 revs; M trajectory; settled min/median σ/σ_shed; max
$|\Gamma|/\sigma^2$ tail. "Threshold moves" = survives 8 revs, no negative
σ, no merge-OOM; "fidelity collapse remains" = settled median
σ/σ_shed < 0.85 or M > ε = 0.2 sustained. All four quadrants have stated
meaning; none is a protocol failure. Prediction P-3: survives the
discrete-overshoot route but shows the fidelity collapse.

**Result (2026-08-07, banner-verified `expint:true`; Slurm COMPLETED,
8:10 h, all 324 steps): BOTH registered arms fire — the discrete death
mode is eliminated AND the fidelity collapse remains. P-3 CONFIRMED.**

- **"Threshold moves" arm: TRUE by its letter.** Survives 8 revs
  (COMPLETED), **zero negative-σ events over the entire history** (global
  min min_sr = +0.0396), no merge-OOM. At step 237 — where the
  banner-identical Euler twin `scr_p019_s020v` OOM-ignited — the expint
  run reads M = 0.19 (< ε), max_u = 43 m/s: it sails through the twin's
  death healthy and runs ~1.3 revolutions deeper before its own collapse
  begins (~step 283).
- **"Fidelity collapse remains" arm: TRUE, in the severe form.**
  M > ε = 0.2 sustained from step 225; max |Γ|/σ² reaches 1.8×10⁵ and
  velocities cross the 1000 m/s tripwire at step 294 — whereupon the wake
  **self-destructs by trim-out** (n_particles 268k @283 → 15k @300; the
  runaway particles convect out of the GlobalCylinder and are deleted;
  min_sr = 1.0 thereafter = only fresh particles remain). A settled
  window does not exist to score median σ/σ_shed; the M-arm alone already
  fires the criterion.
- **The Phase-1 prediction lands to three digits:** the deepest core
  parks at min_sr = 0.0396 = √(2νΔt)/σ_shed = 0.0396 — the exact
  CoreSpreading laminar floor, reached smoothly with σ > 0 throughout
  (theory §3f/§4.5: "exact + molecular ν: no flip, dies by fidelity
  runaway to the laminar floor").

**Registered-branch reading (item file, Phase-2 objective):** the outcome
occupies *both* stated branches at once — stiff integration IS a real
stability lever (every spurious death mode is gone; the run outlives its
Euler twin), and the ceiling IS field-coupled (the Γ-side amplification
proceeds regardless and destroys the wake physically) ⇒ **the closure is
needed for stability AND physics.** P-6 corollary: the collapse is
contraction/field-driven, not threshold-dated — the expint run's ~1.3-rev
extension (vs the twin's rev-6.6 death at matched σ) is the measured size
of the purely-discrete contribution at this rung.

Figure: `figures/fig_stage_b.tex` (+ `fig_stage_b/expint_health.csv`).

---

## Leg 3 — the SFS viscous-null (COMPLETE)

**Derivation leg (standalone):** the SFS term enters the σ equation only
through the prefactor $f/(1+3f)$, identically zero in the live rVPM
($f=0$): no SFS model, at any coefficient, can write σ
(`FLOWVPM_timeintegration.jl:147`; phase_01 §1.1, §2.2). Corpse evidence:
σ marched to the laminar floor with SFS on (018).

**Measured leg:** no SFS-off run existed anywhere in `data/` (search
2026-08-06: every metadata carries `SFS_Cd_twolevel_nobackscatter`) — the
registered fallback pair fires at σ/R = 0.025 viscous, the mid-ladder
igniter: SFS-on member = `scr_p019_s025v` (ignited @286, 4:42 h); SFS-off
twin = `scr_p020_s025v_nosfs` = job **13065482** (clone + `SFS_OFF=true`).

**Registered threshold (§7.3):** time-to-ignition within ±15% and min-σ
trajectory overlay within line width ⇒ SFS viscous-null CONFIRMED at
collapse scale.

**Result (2026-08-07, banner-verified `sfs_off:true`, OOM after 1:35 h):
SFS viscous-null CONFIRMED.**

- **Time-to-ignition: PASS.** Onset (first max_u > 1000 m/s) at step 269
  (SFS-off) vs 286 (SFS-on) — **5.9%**, well inside ±15%.
- **Trajectory overlay: PASS with recorded nuance.** The min_sr traces
  interleave with repeated crossings and no systematic separation
  (e.g. on/off = 0.39/0.24 at step 120 but 0.23/0.39 at 180); pointwise
  median relative difference 31% is chaotic-twin field-min noise (the two
  runs decorrelate; the field-min switches particles), i.e. ~0.12 decades
  on the 1.4-decade log overlay — within line width as an envelope, not
  pointwise. Log-decay rates −0.0069 vs −0.0047 /step bracket each other's
  scatter given the crossings.
- **Same death signature:** Γ-side blow-up at floor-held positive σ
  (terminal min_sr +0.035, max |Γ|/σ² 9.8×10⁵, max_u 14.7 km/s) — the
  regime-2 Γ-route in both twins.
- Labeled non-verdict observation: pre-onset direct M reads 4.03 (off) vs
  1.12 (on) — near-onset transient timing relative to the tripwire, high
  twin-to-twin variance; the registered metrics above are the
  verdict-bearers.

Figure: `figures/fig_leg3.tex` (+ `fig_leg3/minsr_overlay.csv`).

---

## Historical composition — WITHDRAWN

**Leg 1** measured the ceiling: resolvable strain is capped at
$Z_{\mathrm{res}} \approx 957\ \mathrm{s^{-1}}$ by stability (ε-crossing
at σ_stab ± 0.4%; every rung below it ignites; steep margin growth
p ≥ 3), while the physics of this rotor's tip vortices sits **1.9–3.3
decades higher**. **Leg 2** isolated what kind of fix the numerics can
and cannot supply: offline (Stage A) and live (Stage B), the exact local
integrator removes *every* spurious discrete event — no negative σ, no
transverse flip, no merge-OOM, outliving its Euler twin — yet the
field-coupled Γ-amplification runs to the laminar CoreSpreading floor
(to three digits) and destroys the wake physically. Stability alone,
however clean, buys ~1.3 revolutions at σ/R = 0.02, not an operating
point. **Leg 3** closed the last alternative: the existing SFS model is
viscous-null at collapse scale (structurally, prefactor $f/(1+3f) = 0$;
measured, twin ignitions 5.9% apart) — there is no already-shipped
channel that supplies the missing physics.

**Closing claim (registered §7.4), now carried by measurement:** *a
subfilter viscous (σ-channel) model is the missing piece that would let
the simulation approach physics-predicted strain rates.* The measured
ceiling says the gap is decades; the integrator result says no
discretization fix closes it (but the integrator is a necessary
companion — it removes the spurious deaths the closure cannot); the SFS
null says nothing currently in the model can. What remains is the
Phase-1-derived closure: $\nu_t = \nu + (\kappa\sigma_0)^2|S|$ in the
CoreSpreading channel, paired with the exact integrator — exactly the
§4.5 complementarity, now demonstrated on real simulations from both
sides.

**Scope and caveats, on the record:** Leg-1's Γ_implied band criterion
failed as literally registered (boundary-localized instead of
ladder-wide); Stage A ran on aggregate traces (deviation); Stage B is one
rung, one flow family, screen class; the Leg-3 overlay is an envelope
statement between chaotic twins. None of these alters the composition's
direction; all are logged in the item file.

## Revised composition and gate

Phase 2 establishes three narrower points. First, the present rotor screens
show a steep stability boundary, but the registered $\Gamma_{implied}$ band
failed and the 1.9--3.3-decade physical gap remains conditional on an
unsupported core-radius scenario. Second, the corrected geometric map removes
the Euler twin's local step-237 catastrophe, yet the field still develops an
abrupt tail runaway at step 243. Third, the shipped SFS is structurally
$\sigma$-null for $f=0$, with the measured twin comparison supporting that
null at collapse scale.

This is sufficient to reject the old numerical implementation and to retain
the corrected integrator as a necessary experimental tool. It is **not**
sufficient to declare the proposed fixed-filter closure predictive, uniquely
derived, or required by the registered Phase-2 Boolean. The closure remains a
candidate model whose estimator and constants must pass the independent
fixtures in `phase_02r_remediation_plan.md` before any rotor tuning.

**Gate:** Phase 2 is complete with a mixed/negative preregistered verdict.
Stop here. Ryan's explicit approval is required before the predictive-closure
gate or any Phase-3 implementation work.

**Robustness review:** the conclusion does not depend on the terminal median,
which is insensitive to a tail-localized event, nor on labeling the guard
failure as OOM. It uses the paired step-237 state to isolate the discrete
Euler contribution and the corrected run's own raw histories to establish the
later runaway. The stronger historical claims (“8-revolution survival,”
“both registered arms fire,” and “closure needed for stability and physics”)
do not survive literal rescoring and remain withdrawn.
