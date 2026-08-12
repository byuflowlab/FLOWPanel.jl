# 019 — General σ-selection procedure: formalize, then demonstrate on rotor hover

**Opened:** 2026-08-06 (Ryan). **Status:** COMPLETE 2026-08-12 (P1–P6 done;
notebook entry written journals/20260803.md §20260812, Ryan-approved;
drift-arrest v1.1 amendment Ryan-approved and added to P1).

**Item-level approvals:** Technical [x] (2026-08-12); clear-context [ ]; user [ ]

## Objective

Formalize the core-size selection procedure that fell out of the 018 σ-blow-up
investigation (`018_.../sigma_blowup_mechanism.md`; notebook entry "Core Size
Impact on Viscous VPM Stability") as a **general, flow-agnostic algorithm**,
and demonstrate it end-to-end on the DJI9443 hover case as a **fully
reproducible** exercise: initialize σ from the a-priori estimates, probe,
iterate to the stability margin (or show the initializer already lands stable),
and plot (a) the σ-iterate convergence and (b) the stability margin vs σ.

**Data policy (Ryan): reuse existing HPC runs wherever a grid/probe point is
config-matched; run new jobs only where gaps exist.** Reproducibility is
carried by the *pipeline*, not by re-running: one scripted entry point
ingests every data point from its in-repo CSVs with exact provenance (run
name, job id, config hash from `_case_metadata.toml`), recomputes all derived
quantities, and regenerates every figure. A reader with the repo can rebuild
the entire demonstration from raw monitors CSVs; a reader with the cluster
can also re-run any probe. **Quality bar: polished and publishable** — this
item is written and figured as a paper section (methods-grade procedure
statement, publication-quality figures, no campaign jargon in the final
text).

## The procedure (to be formalized as the item's first deliverable)

Inputs: ν_eff, Δt, a loading estimate Γ_v (strongest circulation the wake
concentrates — for a rotor, ~peak bound circulation from BEM or a cheap
coarse probe), and an ambient-strain estimate Z̄ (see P1).

1. **Initialize** σ₀ = max( c_eq·σ_eq, σ_stab ) with
   $$\sigma_{eq} = \sqrt{\nu_{eff}/\bar Z}, \qquad
     \sigma_{stab} = \sqrt{\Gamma_v\,\Delta t/2\pi},$$
   c_eq ≈ 2 (viscous-fidelity margin). σ_stab derives from ΔtZ ≲ 1 with
   Z ~ Γ_v/(2πσ²); status HYPOTHESIS — single-case validation (predicts
   0.031R vs measured 0.029–0.030R ignition boundary on this rotor).
2. **Probe**: short instrumented run at σ_k (screen-class, ~8 revs,
   WakeHealthMonitor on; viscous whenever σ_k < 3σ_eq). Measure the margin
   $$M_k = \max_{\text{steps, particles}} \Delta t\,Z
         = \max\left(-\Delta\sigma/\sigma\right)_{\text{per step}}$$
   — directly observable because the discrete update is Δσ/σ = −ΔtZ.
3. **Iterate**: pre-registered target M ≤ ε. If M_k ≤ ε: done (σ\* = σ_k).
   Else update with the Z ~ 1/σ² Newton-like step
   $$\sigma_{k+1} = \sigma_k \sqrt{M_k/\varepsilon}$$
   and re-probe. Divergence of a probe (tripwire ignition) is itself a
   measurement: σ_k is inside the boundary; bisect upward.
4. **Certify**: one production-length run at σ\* (staged startup, settled
   window, M1+M2-class scoring) — the procedure selects; certification
   verifies.

Open design points to settle in Phase P1 (they are part of "formalize"):
ε's value (margin target; candidate 0.1–0.3 — justify against the measured
boundary's sharpness); how Z̄ is estimated *a priori* vs from the first probe
(bootstrapping note: σ_eq's input can come from the σ₀ probe itself,
one-cycle lag is acceptable); whether M uses the max over all particles or a
high percentile (single-particle max is noisy; the corpse/screens suggest max
is the right tail to watch but p99.9 may iterate more smoothly — decide and
pre-register).

## P1 — Formalized procedure (PRE-REGISTERED 2026-08-06)

This section is the standalone deliverable: a reader should be able to apply
it to a new flow with no 018 context. Deviations discovered during the
demonstration are recorded in the Log, never silently patched here (P6 rule).

### Algorithm

**Inputs.** Effective viscosity ν_eff; timestep Δt; loading scale Γ_v (the
strongest circulation the wake concentrates — for a rotor, peak bound
circulation from BEM or a cheap coarse probe); optionally an ambient-strain
scale Z̄ (may be bootstrapped, see below).

**Scales.**
$$\sigma_{eq} = \sqrt{\nu_{eff}/\bar Z}, \qquad
  \sigma_{stab} = \sqrt{\Gamma_v\,\Delta t/2\pi}.$$
σ_eq is the Burgers strain–diffusion equilibrium core size: cores below it
have no viscous equilibrium and contract toward the per-step floor
$\sqrt{2\nu\Delta t}$. σ_stab is the explicit-update stability scale: a
filament of circulation Γ_v at core scale σ self-strains at
$Z \sim \Gamma_v/(2\pi\sigma^2)$, and the discrete map $\sigma \leftarrow
\sigma(1-\Delta t Z)$, $\Gamma \leftarrow \Gamma(1 - 3\Delta t Z)$ is
unstable for $\Delta t Z > 2$ (and produces negative σ inviscidly for
$\Delta t Z > 1$), so require $\Delta t Z \lesssim 1$ at $\sigma$.

1. **Initialize** $\sigma_0 = \max(c_{eq}\,\sigma_{eq},\; \sigma_{stab})$
   with $c_{eq} = 2$ (viscous-fidelity margin: keeps the strained-core
   balance share $(\sigma_{eq}/\sigma)^2 \le 25\%$). If no a-priori Z̄ is
   available, take $\sigma_0 = \sigma_{stab}$ and check
   $c_{eq}\sigma_{eq} \le \sigma_0$ after the first probe (one-cycle lag).
2. **Probe** at σ_k: a short instrumented run (screen-class; run viscous
   whenever $\sigma_k < 3\sigma_{eq}$), measuring the stability margin
   $$M_k = \max_{\text{steps},\ \text{particles}} \Delta t\, Z
        \;=\; \max\left(-\Delta\sigma/\sigma\right)_{\text{per step}},$$
   directly observable because the discrete core update is exactly
   $\Delta\sigma/\sigma = -\Delta t Z$.
3. **Iterate** against the pre-registered target $M \le \varepsilon$.
   If $M_k \le \varepsilon$: stop, $\sigma^* = \sigma_k$. Else
   $$\sigma_{k+1} = \sigma_k\sqrt{M_k/\varepsilon}$$
   (Newton-like, from $Z_{\max} \sim 1/\sigma^2$) and re-probe. A diverging
   probe (tripwire ignition) is itself a measurement — σ_k is inside the
   unstable region; take the update from the last pre-ignition M reading, or
   bisect upward if none is clean.
4. **Certify**: one production-length run at σ\*. The procedure *selects*;
   certification *verifies*. **Survival of a probe is never the acceptance
   metric — M is.** (Measured case in point: the 0.0299R screen *survived*
   its 8 revs yet the production twin ignited at step ~501; the probe's
   margin, min σ/σ_shed slid to 0.082 ≈ cumulative 12× contraction, already
   flagged it.) *(v1.1, see amendment below: certification additionally
   requires contraction arrest — per-step M alone is drift-blind.)*

### Pre-registered choices (fixed before the demonstration)

- **ε = 0.2.** Justification: the hard failure thresholds are ΔtZ = 1
  (inviscid, σ sign flip) and ΔtZ = 2 (unconditional); ε = 0.2 leaves 5–10×
  headroom. The measured boundary is sharp in σ (0.0299R survives a screen
  while 0.0291R ignites — ~3% apart) but time-to-ignition shifts strongly
  with startup class and run length (screen step 284 vs production ~501 at
  matched σ), so ε must absorb probe-vs-production duration bias; 5–10×
  headroom with $Z \sim 1/\sigma^2$ costs only $\sqrt{5}$–$\sqrt{10}$ ≈
  2.2–3.2× in σ relative to the cliff edge, an acceptable resolution price.
- **M estimator.** Primary: the direct per-particle readout
  `max_dtZ` = $\max_p \Delta t Z_p$ (WakeHealthMonitor column, added for
  019; same Z as the σ update, MM4 convention). Fallback for reused runs
  that predate the column: the level proxy from `min_sigma_ratio` —
  $\tilde M = \max_{\text{steps}}\left[-\Delta \ln(\text{min\_sr})\right]$
  over downward steps only, valid only while min_sr > 0; it is biased low
  when the field-min switches particles (a fresh particle enters at 1.0)
  and undefined post-ignition. Method is a labeled attribute of every
  reported M. Cross-validation of proxy vs direct on the first new probe
  carrying both is a required deliverable. The **max** is the pre-registered
  statistic (the ignition tail is what kills); p99.9 is reported alongside
  for smoothness commentary, not used in the update rule.
- **Z̄ estimation.** A-priori path: from a cheap coarse probe or an existing
  run at any stable σ (ambient ΔtZ readout); bootstrap path: from the σ₀
  probe itself, accepting a one-cycle lag on the $c_{eq}\sigma_{eq}$ check.
  For this demonstration Z̄ is *recomputed* from probe data by the pipeline,
  not imported from 018 prose (the 018 value 3.2 s⁻¹ serves only as an
  anchor to reproduce).
- **Probe class.** Screen-class (8 revs, harsh pulse-less startup, tripwire
  + `max_dtZ` on) for all *new* probes, uniformly; reused points keep their
  native class as a labeled attribute. Startup class is plotted, never
  pooled silently.

### P1 addendum (2026-08-06, post-registration extension; ε/M/update rule unchanged) — strain resolvability

Raised by Ryan: the resolvable strain is capped by σ itself (coarser cores
smooth the velocity field), so Z̄ measured at σ is Z̄(σ), not the flow's
"true" strain. Three consequences, now part of the procedure's
interpretation (they add diagnostics; they do not alter the pre-registered
ε, M, or update rule):

1. **For stability, resolved strain is the correct input, not a bias.** The
   discrete σ/Γ updates respond to the *resolved* Z — the blow-up is a
   phenomenon of the resolved field — so probing M at the operating σ is
   self-consistent by construction. No correction is needed or wanted.
2. **For the σ_eq arm, coarse probing is conservative, and the honest
   statement is a self-consistency condition.** Z̄(σ) decreasing in σ makes
   σ_eq(σ) = √(ν/Z̄(σ)) increasing in σ; a coarse probe therefore
   over-estimates σ_eq and pushes σ up. The bootstrap should be read as a
   fixed-point condition σ\* ≥ 2σ_eq(Z̄(σ\*)) evaluated with the certifying
   probe's own Z̄. If that fixed-point iteration runs away (Z̄ steepening
   faster than 1/σ² as σ drops), the flow has **no σ that is simultaneously
   stable and viscously faithful at this (Δt, ν)** — a named failure mode
   the procedure should report rather than paper over.
3. **Unresolved strain is detectable from data already collected — the
   M(σ) scaling exponent.** Fit M(σ) ∝ σ^(−p) across the probe ladder:
   p ≈ 2 means the dominant strainer is a compact structure below σ
   (Γ_implied = 2πσ²M/Δt independent of σ — the σ_stab filament regime:
   strain UNresolved, stability governed by σ_stab); p ≈ 0 means
   ambient-dominated (resolved); intermediate p = partial resolution. The
   pipeline will report p and Γ_implied(σ) from the margin curve, and the
   P5 margin figure gains the fitted-slope annotation. Pre-registered
   interpretation: Γ_implied ≈ Γ_v over the ladder validates the σ_stab
   model's mechanism, not just its number.

   Caveat recorded up front: the procedure certifies *stable and
   self-consistent at the chosen resolution* — it does not certify that σ
   resolves all strain-carrying scales of the physical flow. Below σ_stab
   that resolution is unreachable at this Δt regardless (ignition precedes
   convergence), which is why σ is a model parameter here, not a
   convergence axis (018 ruling).

**Remedies when strain is unresolved (p ≈ 2) or the fixed point runs away
(Ryan question, 2026-08-06).** Ranked by soundness for this flow class:

1. *Δt refinement is the nominal lever but empirically suspect here.* In
   the σ_stab model σ_min ∝ √Δt, so resolvable strain ∝ 1/Δt — but 018
   measured the contraction-route boundary to be nearly Δt-invariant
   (attainable strain at floored cores ∝ Δt^(−3/2) vs threshold ∝ Δt^(−1);
   nt36/nt72 ignited at the same physical time; nt144 at 0.0291R degrading).
   Whether a *genuinely stable* window opens between σ_stab(Δt/4) and the
   old boundary is exactly what `fid144` (13061089) measures: if it
   survives, Δt buys real (expensive, ~Δt^(−2) cost) resolution; if it
   ignites, the boundary is contraction-controlled and Δt is NOT the
   adjustment. Even in the best case, resolving the physical tip-vortex
   core (~0.01–0.05c vs our σ ladder 0.12–0.68c) is 1–2 decades away —
   Δt cannot close that.
2. *Raise ν_eff so σ_eq becomes an honest subgrid cutoff (the LES stance) —
   the principled fix for fixed-point runaway.* Choosing ν_eff ≈
   (σ_target/2)²·Z̄(σ_target) places 2σ_eq = σ_target and closes the fixed
   point by construction; core spreading then *is* the subgrid transport
   model rather than laminar physics with molecular ν. The knob exists
   (`WAKE_NU_FACTOR`); the cost is a modeling claim needing its own
   validation (loads-level ν-sensitivity was measured NULL at L1, which is
   encouraging: the wake tolerates ν_eff inflation aerodynamically).
3. *A σ-aware subgrid closure / ΔtZ guard (code change, FLOWVPM).* SFS
   currently cannot touch σ (prefactor exactly 0) — the model has no
   channel for sub-σ cascade except CoreSpreading. 018 §6's mechanism-level
   lever (per-particle ΔtZ cap / clipping the −3ΔtZΓ increment) is the
   crude version; a real closure transferring unresolved strain into core
   growth or Γ decay is the right long-term answer. Out of 019 scope —
   **staged as item 020 (`020_sigma_aware_subgrid_closure.md`), Ryan
   2026-08-06**, with gated theory → implementation → validation phases;
   its Phase-3 validation harness is this item's regime map re-run.
4. *Accept and bound (always do this regardless):* report p and
   Γ_implied(σ), state which structures are modeled-not-resolved, and
   demonstrate QoI insensitivity (or measure the sensitivity, e.g. the
   σ-axis 3-lobe mode ε_Γ 8.7% b0→L1 shows loads are NOT fully insensitive
   — that residual is 018's error-budget line, not a σ-procedure defect).

### P1 amendment v1.1 (2026-08-12, Ryan-approved; post-demonstration — NOT part of the pre-registered v1.0, which stands unmodified above)

**Certification (step 4) requires BOTH criteria:**

1. $M \le \varepsilon$ (per-step margin, as in v1.0), AND
2. **contraction arrest**: $d\,\ln(\min \sigma/\sigma_{shed})/dt \to 0$
   over the settle window (no monotonic cumulative drift).

Basis (P3 verdict, 2026-08-11): the k=0 certification at σ₀ = 0.0356R
held $M \le \varepsilon$ (0.004–0.05) for 25 revolutions while
min σ/σ_shed slid 0.80 → 0.13 (≈ −0.065/rev in ln, never arresting)
into terminal runaway. Per-step M ≤ ε is necessary but not sufficient
at production horizon; the cumulative-contraction drift is the leading
indicator. Screen-class probes *select* (steps 1–3 unchanged); only the
drift criterion *certifies* at long horizon. Measured duration bias this
absorbs: production boundary (0.0349, 0.0381]R ≈ 1.15–1.25× the screen
ε-crossing (0.0299, 0.0312]R.

### Validated scope (to be finalized in P6)

One flow family (small hover rotor, explicit Euler VPM with CoreSpreading).
σ_stab carries a flow-dependent prefactor (filament model); what would
falsify it elsewhere: a flow whose measured ignition boundary departs from
$\sqrt{\Gamma_v \Delta t/2\pi}$ by more than the prefactor band this
demonstration measures, at matched probe protocol.

## P4b — viscous-fidelity discriminator at NT=144 (added 2026-08-06, Ryan)

**Question under test.** The c_eq arm of the initializer asserts: below
~2σ_eq a run can be *stable* (M ≤ ε) yet *physically compromised* — its
core-size field migrates to the local Burgers equilibrium σ\*(Z) (a laminar,
ν_eff-set scale) instead of remaining the user-chosen shed σ. At NT=36 this
is untestable here because σ_stab ≈ 0.031R sits above the band; at NT=144,
σ_stab(Δt/4) = 0.0155R drops below it, opening a stable-but-in-band window.

**The run.** `scr_p019_fid144`: NT=144 (Δt = 7.716e−5 s), OV 2.4 / PPS 5 →
σ/R = 0.02094 = **1.18σ_eq = 0.59·(2σ_eq) = 1.35·σ_stab(NT144)**; viscous
(CoreSpreading β=1e9), D=3.4, N=1, rlxf 0.075, 8 revs = 1152 steps,
`max_dtZ` on, screen class.

**Pre-registered outcomes.**
1. *Stability arm*: survives with M ≤ ε = 0.2 (it sits 1.35× above the Δt/4
   stability scale). If it instead ignites, that is itself a finding: the
   long-horizon boundary does not scale as √Δt (already hinted by
   `scr_ufdt_nt144` at 0.0291R degrading), and σ_stab must be read as a
   per-Δt *shed-scale* criterion, not the full boundary — recorded, not
   patched.
2. *Fidelity arm*: the σ/σ_shed distribution migrates materially toward
   σ\*(Z_local) — pre-registered discriminator: the settled-wake **median**
   σ/σ_shed drops well below the above-band reference's median (threshold:
   in-band median < 0.85 while reference median ≥ 0.95), and the low tail
   parks near σ_eq/σ_shed ≈ 0.85 (ambient) to far below (strained regions)
   rather than at the reference's level.

**Tail-vs-bulk honesty note (from discussion with Ryan 2026-08-06).** How
confident are we that *above* 2σ_eq cores do NOT drop toward σ_eq? Graded,
and tail ≠ bulk. Evidence in hand: at B0 (σ_shed = 4.9σ_eq, viscous) the
tripwire **minimum** settles at min_sr ≈ 0.35 ⇒ even the single worst core
sits at 1.7σ_eq — nothing reaches σ_eq. But at the hedge (σ_shed = 1.96σ_eq,
right at the band edge) min_sr ≈ 0.32–0.35 ⇒ the worst core is at
**0.62σ_eq — the tail already dips below σ_eq at 2σ_eq**. The protection at
coarse σ is not only the (σ_eq/σ)² balance-share argument: coarser σ smooths
the velocity field and caps resolvable Z̄(σ) (018's strain-resolution
margin, measured: L2's p99 shrink rate ~6× L1's), so σ\*(Z̄(σ)) tracks
closer to σ_shed. The *bulk* distribution above the band is bounded below by
min_sr but otherwise unmeasured — hence the scoring below compares full
histograms, and the claim to be evaluated is about medians/quartiles, not
minima.

**Scoring plan.** From retained VTP restart sets (last ~10 steps are kept by
the sweep policy) + monitors: (a) σ/σ_shed histogram of `scr_p019_fid144`
settled wake vs an above-band reference — primary: `scr_p019_s038v`
(σ = 2.14σ_eq, NT=36, viscous, same screen class); secondary:
`p018_ufront_s035_visc` restart set (production class, labeled);
(b) max_u and max|Γ|/σ² monitor tails (artifact-velocity check);
(c) M from `max_dtZ` (stability verdict); (d) Γ(r/R) and CT deltas reported
but NOT verdict-bearing (screens fail the aero calibration gate). VTPs of
fid144's final steps must be added to the protect list at harvest.

## Current status

**2026-08-12 — SIGN-OFF SESSION COMPLETE (see Log): notebook entry
written (journals/20260803.md § 20260812), drift-arrest v1.1 approved
into P1, Technical Completion marked. Still open: clear-context + user
approval checkboxes; item 020 start decision (Ryan).** The 2026-08-11
status below stands as the analysis record.

**2026-08-11 — MEASUREMENT + ANALYSIS COMPLETE; awaiting Ryan.** All 10
Campaign E runs terminal and harvested; P2/P3/P4/P4b verdicts written as
dated sections below; P5 figures authored (4 .tex, compile-verified);
P6 reconciliation written. Nothing is in flight and no jobs are needed.
Open items — ALL gated on Ryan: (1) notebook entry (policy: ask before
writing; proposal offered 2026-08-11: short context + results table +
margin-curve and P4b figures + conclusion bullets); (2) item-level
approval checkboxes; (3) the drift-arrest v1.1 amendment (P6 §
reconciliation item 4); (4) optional checkpoint commit (~150 modified/
untracked files); (5) item 020 start decision. Headline conclusions:
see the "cliff notes" ordering in the 2026-08-11 P6/P3/P4b sections —
σ_stab validated at the screen ε-crossing (~2%); σ₀ failed production
certification via cumulative contraction drift (not per-step margin) ⇒
σ\* = 0.0381R via the bisect-upward branch onto L1's certification;
in-band failure mode is tail-driven ignition (bulk static); proxy
estimator demoted (3–140× low); Δt not a mitigation lever.

## Next actions for a fresh agent (superseding note, 2026-08-11)

The 2026-08-06 handoff list below is COMPLETE — do not re-execute it.
A fresh agent should: read "Current status" above, then wait for Ryan's
decisions on the five open items. If drafting the notebook after
approval: source material = the 2026-08-11 dated sections + P5 captions;
target = journals entry "Core Size Impact on Viscous VPM Stability"
(notebook policy in ~/.claude/CLAUDE.md applies — append-only, propose
header first). To regenerate anything: `python3
scripts/p019_sigma_procedure.py` (add `--allow-missing-provenance`; two
non-019 legacy runs still lack provenance files), then `pdflatex` each
`figures/fig_*.tex` from the figures dir. P4b reduction rerun (only if
VTPs change): `scripts/vtp_sigma_reduce.py` on the cluster (usage in
its docstring; σ_shed values in the P4b scoring section).

## Archived status (2026-08-06 — superseded)

**2026-08-06 — ACTIVE; P4 grid partially in flight.** Campaign E screen jobs
(banner-verified 11:40 MDT): 13060963 `scr_p019_s015v`, 13060964
`scr_p019_s020v`, 13060965 `scr_p019_s025v`, 13060966 `scr_p019_s030v`
(all viscous, CoreSpreading β=1e9, `WAKE_HEALTH_DTZ=true`). Submission queue
(automatic, 10-job cap, checks every 15 min), in order: `fid144` (P4b, 36 h
walltime), then `s015`, `s025`, `s038`, `s038v`. **Clearance ruling (Ryan
2026-08-06): C1 clearance outranks the tip chord-band cap.** Every Campaign
E arm including `s038`/`s038v` now runs uniform d/σ = 3.4 — exactly at the
C1 bound (phase_12/13: admissible d/σ\* 2.6 median / 3.4 p95),
dt-independent by construction. Consequence: the s038 pair exceeds 1.5
local chords at the tip (the cap the hedge respected by dropping to D=3.0);
this is a deliberate, labeled attribute — note it when comparing s038v
against the D=3.0 hedge. The D=3.0 arms were never launched (caught in
queue before submission). `scr_p019_fid144` = job **13061089** (RUNNING,
banner-verified 12:13 MDT). Queue anomaly, needs Ryan's awareness: five
new `fp-018-ufront-*` production jobs (13061047–51) + `fm029p2` appeared
~12:05 from outside this session, putting the queue at 16 > the 10-job cap;
fid144 was submitted anyway on Ryan's explicit launch order. **Cap lifted
by Ryan (2026-08-06 ~12:20): all four remaining screens submitted and
RUNNING, banners verified (all D=3.4)** — 13061166 `s015`, 13061167 `s025`,
13061168 `s038`, 13061169 `s038v`. Full Campaign E grid (9 jobs incl.
fid144) is now in flight simultaneously. Monitor `max_dtZ` column
implemented, test-gated (exact match to FLOWVPM `_euler` −Δσ/σ, rtol 1e-10),
deployed to cluster with full md5 verification.

## Archived handoff (2026-08-06 ~13:00 MDT — COMPLETE, superseded above)

Read this item top to bottom first (P1 pre-registration + addendum, P2/P3
first pass, P4b) — it is the single source of truth. Session-local watchers
and submitters did NOT survive the reset: poll jobs yourself.

1. **Poll the 9 in-flight Campaign E jobs** (`ssh orc`, needs live
   ControlMaster socket; `squeue -u rander39` via `bash -lc`): 13060963
   `scr_p019_s015v`, 13060964 `s020v`, 13060965 `s025v`, 13060966 `s030v`
   (started 11:40), 13061089 `fid144` (36 h wall, started 12:13), 13061166
   `s015`, 13061167 `s025`, 13061168 `s038`, 13061169 `s038v` (started
   ~12:25), **and 13064696 `scr_p019_sstab`** (started 15:53; Ryan-ordered
   probe at exactly σ_stab: OV2.5/PPS14 → σ/R 0.03117 = σ_stab+0.4%,
   viscous, D=3.4 — answers "would σ_stab alone work as initializer" by
   measurement; add it to the pipeline REGISTRY and the margin figure).
   All banners verified correct (all D=3.4, `WAKE_HEALTH_DTZ=true`
   confirmed live in CSV headers). Screens end 2–7 h after start (igniters
   die early — that's data, not failure).
2. **Harvest each finished run**: rsync `data/<run>/monitors/*.csv` +
   `<run>_case_metadata.toml` (completed runs) or build
   `<run>_banner_config.toml` from the slurm log banner (ignited runs —
   schema examples in any `data/scr_uf*/`); record sacct state/elapsed/
   MaxRSS. Protect fid144's final-step VTPs (P4b scoring needs them; add to
   `018_.../vtk_protect_list.txt`).
3. **Extend `scripts/p019_sigma_procedure.py`**: add the 9 new runs to
   REGISTRY (fid144 has NT=144 ⇒ own dt); parse the `max_dtZ` column
   (direct M) alongside the proxy; implement the **proxy-vs-direct
   cross-validation** (required deliverable, P1); implement the **p-exponent
   fit** M(σ) ∝ σ^(−p) + Γ_implied(σ) = 2πσ²M/Δt per the P1 addendum
   (screen-class, per viscosity row). Rerun full pipeline; append results
   here as a dated section.
4. **P3 gate**: when 13058988 (`p018_ufront_s035_visc`) completes, harvest
   it, recompute M(σ₀) with final data, and write the P3 verdict (k=0 PASS
   currently provisional at M=0.0198 proxy).
5. **P4b scoring** (after fid144 ends): outcomes pre-registered in section
   P4b — stability verdict from max_dtZ, σ/σ_shed histograms from VTPs vs
   `scr_p019_s038v` reference. Either branch (survive/ignite) is a finding;
   ignite ⇒ also record in 020's gap-demonstration evidence.
6. **P5 figures**: TikZ per Ryan's figure policy (standalone .tex +
   same-named CSV dir, pdflatex-compilable), data dirs already emitted
   under `019_sigma_selection_procedure/figures/`; regenerate CSVs from the
   pipeline after step 3 before authoring.
7. **P6 last**: notebook update requires Ryan's approval BEFORE writing
   (notebook policy; Ryan 2026-08-06: explicitly deferred — do NOT write
   the notebook yet); also close out INDEX cell + procedure-text
   reconciliation. The σ_stab probe is no longer optional — it is in
   flight (13064696, item 1 above).

Related: item 020 (σ-aware subgrid closure) is STAGED with a discussion
record and 4 gated phases — do not start it without Ryan; its Phase 2 will
reuse this item's regime map.

**2026-08-06 — kickoff state.** P1 pre-registration written (above); WakeHealthMonitor
`max_dtZ` column, Slurm-banner provenance harvest, and the
`scripts/p019_sigma_procedure.py` pipeline in progress; P4 full grid
(σ/R ∈ {0.015, 0.020, 0.025, 0.030, 0.038} × {inviscid, viscous},
new jobs where no config-matched 018 run exists) authorized for submission
this session. P3 certification gated on hedge job 13058988
(`p018_ufront_s035_visc`, healthy in flight). Note for P2: the 018 record
carries both σ_eq ≈ 0.0175R and 0.0182R (ν 1.5e−5 vs 1.4334e−5 m²/s) — the
pipeline recomputes from declared constants and resolves this; do not
hand-pick.

**2026-08-06 ~14:00 MDT (Ryan-directed note, from the 018 session).**
OPEN QUESTION for the demonstration: the clean-dt screen family
`scr_ufdt_nt{36,72,144}` (σ/R 0.0291 inviscid, N=1 D=3.4, PPS 12/6/3)
**all ignited at ≈7–8 physical revolutions** (steps 284 / 490 / ~1104; ledger
§13:05 + §13:20) even though σ_stab = √(Γ_vΔt/2π) shrinks as √Δt — so the
nt72 (σ_stab ≈ 0.022R) and nt144 (≈ 0.0155R) rungs sat **on the stable side
of the σ_stab criterion and blew up anyway**. Investigate why the criterion
fails there — Ryan's suggested lead: **the second stability criterion**
(the σ_eq strain–diffusion side of the initializer / the ΔtZ<2 map bound, as
distinct from the ΔtZ<1 σ-flip) may be the operative one; the constant-revs
result (dσ/dt = −Zσ is Δt-independent, so the flow grinds σ into the runaway
in fixed physical time) is the candidate mechanism to reconcile. This bears
directly on P1's framing: σ_stab is an instantaneous threshold, not a static
safe/unsafe classifier of the shed σ. Evidence preserved: all three corpses
are on the VTK protect list with their final-10-step ignition windows
(nt36 275–284, nt72 481–490, nt144 1097–1106).

## Phases

| phase | deliverable |
| --- | --- |
| P1 | Procedure formalized: algorithm box + measurement definitions (Γ_v, Z̄, M) + pre-registered ε and update rule; **tripwire extension if needed** — WakeHealthMonitor currently logs `min_sigma_ratio` per step (level); add a max-per-step-contraction column (= ΔtZ_max readout) if the level-difference proxy proves too noisy. Any src change: default-off/bit-identical, test-gated |
| P2 | Rotor demonstration — initialize: compute σ_eq, σ_stab from stated inputs (Γ_v from measured peak bound circulation ≈ 0.28 m²/s **and** independently from a BEM estimate, to show the a-priori path); expected σ₀ = max(2·0.0175R, 0.031R) = **0.035R** |
| P3 | Probe + iterate: margin measured at σ₀ (and iterates if M > ε). **Reuse-first**: `p018_ufront_s035_visc` IS the σ₀ = 0.035R viscous probe if its config matches the P1 spec (verify from its metadata); `p018_ufront_n1_visc` covers 0.03R viscous; new probes only where the iteration demands an unprobed σ |
| P4 | **Parameter-space exploration showcasing the blow-up**: a probe grid spanning the unstable regime(s), not just the stable side — σ ladder crossing the boundary (e.g. σ/R ≈ 0.015, 0.020, 0.025, 0.030, 0.038) × {inviscid, viscous}, each screen-class with tripwire. Deliverables per cell: outcome (stable / ignites), time-to-ignition, peak M, and the ignition trajectory (min σ/σ_shed and max\|Γ\|/σ² vs step) for the dying cells — the measured companion to the mechanism flow chart. Optional third axis (one row): NT=72 at a dying σ to *show* the Δt non-mitigation. Several cells already exist as 018 runs (screens, L2, L1) — reuse where the config matches the grid exactly, run fresh where it does not, label provenance on the figure |
| P5 | Plots + write-up: (a) σ-iterate convergence (σ_k and M_k per iteration vs ε band); (b) stability margin M(σ) across all probed σ with the measured ignition boundary and σ_eq/σ_stab marked; (c) the P4 regime map — outcome grid over (σ, viscosity) with time-to-ignition contours/colors and the two initializer scales overlaid; enriched (clearly labeled) with the 018 historical points |
| P6 | Close-out: procedure text finalized against what the demonstration actually required (any deviation from P1's pre-registration is recorded, not silently patched); **update the notebook's "Core Size Impact on Viscous VPM Stability" entry** (journals/20260803.md, # 20260805) with the demonstrated procedure — measured margin curve, regime map, σ_stab validated-or-revised — keeping that entry the single human-readable synthesis (concise, takeaway-focused; a dated pointer back to this item for detail) |

## Reproducibility + polish contract

- One entry point: `scripts/p019_sigma_procedure.py` — declares all constants
  (ν, Δt, R, RPM, Γ_v source), computes the initializers, ingests every grid/
  probe point from in-repo monitors CSVs with per-point provenance (run name,
  job id, `_case_metadata.toml` config check — a mismatched config is a hard
  error, not a warning), scores margins, applies the update rule, emits job
  scripts for any gap runs, and regenerates all figures.
- No numbers imported from conversation history or ledger prose — everything
  recomputed from CSVs or derived from declared constants.
- New gap probes use the phase_14 screen runner pattern unless P1 decides the
  staged startup is required for a faithful margin — decide once,
  pre-register, apply uniformly; startup class is a labeled attribute of
  every point on every figure (screen vs production startup margins are known
  to differ — the boundary moved ~step 284 → ~step 500 between them).
- Figures at publication standard: vector or ≥300 dpi, consistent fonts/
  units, self-contained captions drafted in the item file; the final write-up
  reads as a methods/results section (no p018_* tags in figure text — a
  provenance appendix table maps display labels to runs/jobs).

## Grounding from 018 (context, not inputs)

Known anchors the demonstration should reproduce or be consistent with:
ignition boundary σ/R ≈ 0.030 (inviscid, screen and production startup);
L1 (0.038R) stable 40 revs viscous; L2 (0.019R) ignites viscous; hedge
0.035R viscous healthy in flight; ΔtZ ambient ~1e-3 at B0. Mechanism and
confidence table: `018_.../sigma_blowup_mechanism.md`. If the demonstration
*contradicts* an anchor, that is a finding for the 018 record, not something
to reconcile silently.

## Acceptance criteria

1. The algorithm is stated standalone (a reader needs no 018 context to
   apply it to a new flow).
2. The rotor demonstration runs from the single entry point; every
   verdict-bearing number traces to its own probe outputs.
3. All three plots produced; the margin plot shows σ₀, any iterates, ε, the
   measured boundary, and σ_eq/σ_stab; the regime map covers both sides of
   the boundary including at least two measured ignitions in the unstable
   region (with their trajectories), with per-cell provenance labeled.
4. Explicit statement of the procedure's validated scope (one flow family)
   and what would falsify σ_stab elsewhere.
5. The notebook's "Core Size Impact on Viscous VPM Stability" entry updated
   with the demonstrated results (P6) — it remains the concise human-readable
   synthesis; this item holds the detail.

## P2/P3 — first pipeline pass (2026-08-06, existing data only)

`scripts/p019_sigma_procedure.py` (stdlib+tomllib; subcommands constants /
gamma_v / initializers / ingest / margins / iterate / gaps / figures;
byte-identical on rerun; hard-error config contract active, known-pair
σ-derivation asserts pass).

**P2 initializers (all recomputed, nothing imported from prose):**
Γ_v = 0.27792 m²/s (p018_L1_visc monitor03, settled last-25% window,
per-station median |Γ_te|, peak over stations, max over blades — matches the
~0.28 anchor); ν_eff = 1.433418e−5 m²/s (metadata `wake_nu`, authoritative;
the 018-prose 1.5e−5 printed but unused); dt = 3.0864e−4 s cross-checked.

| scale | m | σ/R |
| --- | --- | --- |
| σ_eq = √(ν_eff/Z̄) | 0.002116 | 0.01779 |
| σ_stab = √(Γ_v dt/2π) | 0.003695 | 0.03105 |
| **σ₀ = max(2σ_eq, σ_stab)** | 0.004233 | **0.03557** |

σ₀ is set by the 2σ_eq arm (2σ_eq = 0.0356 > σ_stab = 0.0311), consistent
with the P2 expectation of ~0.035R.

**P3 iterate table (k=0):** nearest probe to σ₀ is `p018_ufront_s035_visc`
(σ/R 0.03491, production class, in flight): M = 1.98e−2 ≤ ε = 0.2 →
**PASS, σ\* = σ₀ at k=0** — pending (a) the hedge run completing (P3 gate)
and (b) proxy→direct cross-validation, since every M so far is the fallback
proxy (no run carries `max_dtZ` yet; the Campaign E screens will).

**Margin snapshot (proxy method, pre-onset window):** survivors read M
0.0097–0.12 (B0 0.0097, 0.050R 0.034, s035 0.020); ignited runs read
0.022–1.2. Known proxy pathologies confirmed and recorded: (1) cliff-blind —
pre-onset proxy on violently ignited runs is small (d2p6: M 0.036 then
min_sr −38); (2) boundary bias — `scr_uf_n1_d3p4` (0.0299R screen survivor)
reads M 0.12 < ε yet its production twin ignited: **a proxy-based ε pass at
0.0299R is explicitly untrusted; the direct-column cross-validation is
load-bearing**; (3) p99.9 degenerates to max at screen lengths (~200–500
samples) — max is the reported statistic, as pre-registered. Also already
visible in ingested data: NT=72 vs NT=36 at σ=0.0291R ignite at ~the same
physical time (steps 447 vs 277) — the Δt non-mitigation row of the regime
map needs no new run.

Figure data emitted: `019_sigma_selection_procedure/figures/`
(fig_margin_curve/, fig_regime_map/ incl. 9 ignition trajectories,
fig_iterates/, provenance_appendix.csv). TikZ authoring pending final data.

## 2026-08-06 ~16:15 MDT — Campaign E first harvest (3 ignitions) + pipeline extension

**Poll (16:05 MDT).** Still RUNNING: `s025v`, `s030v` (started 11:40),
`s025`, `s038`, `s038v` (~12:25), `fid144` (12:13; 36 h wall), `sstab`
(15:53), and the P3 gate `p018_ufront_s035_visc` 13058988 (13.5 h elapsed).
Finished and harvested (monitors CSVs + full `.metadata.toml` rsynced;
hand-built `<run>_banner_config.toml` provenance from the slurm banners;
slurm logs mirrored to `data/p019_slurm_logs/`):

| run | job | state | elapsed | MaxRSS | onset step | onset rev (of 9 incl. spinup) |
| --- | --- | --- | --- | --- | --- | --- |
| `scr_p019_s015` (inv) | 13061166 | OUT_OF_MEMORY | 0:55:30 | 15.9 G | 97 | 2.7 |
| `scr_p019_s015v` | 13060963 | OUT_OF_MEMORY | 2:39:55 | 67.0 G | 196 | 5.4 |
| `scr_p019_s020v` | 13060964 | OUT_OF_MEMORY | 3:59:12 | 66.5 G | 237 | 6.6 |

All three OOMs are ignitions, not resource failures: max_u crosses the
1000 m/s tripwire at the onset steps above (terminal max_u 1.8e4 / 1.2e4 /
4.4e4 m/s), and the inviscid s015 additionally reaches min_sigma < 0. The
viscous 0.020R ignition corroborates the L2 (0.0193R, production) viscous
ignition at screen class; inviscid 0.015R dies fastest. VTK protect list:
added `scr_p019_fid144` + `scr_p019_s038v` (P4b scoring inputs), pushed to
cluster.

**Pipeline extension (step 3 of the handoff list).** 10 new REGISTRY
entries (9 Campaign E + `sstab`, with job ids; fid144 carries NT=144 and
its own dt); direct `max_dtZ` parsed; BOTH estimators now computed per run
(`M_direct`/`M_proxy` retained side by side); proxy-vs-direct
cross-validation implemented (`xval_proxy_vs_direct.csv` + printed table);
p-exponent fit + Γ_implied implemented (`pfit.csv`, `gamma_implied.csv`);
ignition-trajectory CSVs now carry `max_dtZ`. Idempotency re-verified
(md5-identical repeat runs). All 10 P4 grid cells are now covered — the
gap manifest is empty.

**Proxy-vs-direct cross-validation (first data carrying both):**

| run | σ/R | outcome | M_direct | M_proxy | proxy/direct |
| --- | --- | --- | --- | --- | --- |
| `scr_p019_s015` | 0.01496 | ignited | 1.75 | 0.151 | 0.086 |
| `scr_p019_s015v` | 0.01496 | ignited | 12.6 | 1.38 | 0.110 |
| `scr_p019_s020v` | 0.01995 | ignited | 28.2 | 1.42 | 0.050 |

The proxy under-reads the direct margin by 9–20× on these windows —
consistent with the pre-registered particle-switch bias, and confirming
the direct column is load-bearing. Caveat: all three points are pre-onset
windows of *ignited* runs (divergence transients); the survivor-side
cross-validation — the one that bears on P3's k=0 PASS — lands when
`s025v`/`s030v`/`s038v`/`sstab` complete.

**p-fit: survivor-only ruling (recorded deviation, P6 rule).** The P1
addendum pre-registered "fit M(σ) ∝ σ^(−p) across the probe ladder"
without restricting to survivors. Fitting the two ignited viscous points
gives p = −2.8 (M *increasing* with σ) with Γ_implied 0.81 and 3.24 m²/s
(≫ Γ_v = 0.278): pre-onset M of an ignited run (12.6, 28.2 ≫ the ΔtZ = 2
unconditional threshold) is a divergence transient, not the stationary
margin the scaling law is about. Ruling (recorded here, not silently
patched into P1): the fit uses survivor points only; ignited direct-M
points are reported alongside, unfitted. Fit currently skipped (0 fittable
survivors carry `max_dtZ` yet); it will populate from the in-flight
survivors.

## 2026-08-06 ~18:20 MDT — wave-1 terminal: first survivor with direct M

Watcher fired 18:09; harvested. `scr_p019_s030v` **COMPLETED** all 8 revs
(6:24 h, case_metadata emitted). `scr_p019_s025v` ignited (max_u 5.9e4 m/s
at step 288; Slurm state FAILED = Julia crash one step post-ignition, not
a resource event; banner TOML records this). `scr_p019_s025` ignited
(OOM, min_sigma < 0 by step 289). Still running: `s038`, `s038v`,
`sstab`, `fid144`, P3 gate 13058988.

**Headline: the 0.030R viscous survivor carries M_direct = 0.295 > ε =
0.2 → FAIL by the margin rule despite surviving its 8 revs.** This is the
pre-registered principle ("survival is never the acceptance metric — M
is") now *measured* with the direct column, and it is consistent with the
production twin 13058534 (0.0299R viscous) igniting at step ~501: the
screen's margin already disqualified the point. Its proxy reads 0.079 —
which would have (wrongly) PASSED ε. Survivor-side proxy/direct = 0.27.

**Cross-validation status (6 runs, both estimators):** ignited windows
proxy/direct 0.05–0.11; boundary-adjacent ignitions 0.53–1.38 (s025v/
s025 — proxy noise-dominated near onset); the one survivor 0.27. Ruling
implied for P3: the k=0 PASS (hedge, proxy M = 0.0198) scales to a
direct-equivalent ~0.07–0.10 at the survivor ratio — likely still < ε,
but the verdict must come from direct data: 13058988 predates the
`max_dtZ` column, so the P3 certification will lean on `sstab` (0.0312R,
direct, in flight) plus the hedge's own proxy with the measured
survivor-ratio band quoted as its uncertainty. Note for the P3 write-up.

**Boundary sharpening (viscous row, screen class, direct M):** 0.0249R
ignites (M pre-onset 1.12), 0.0299R survives at M 0.295, σ_stab = 0.0311R,
hedge 0.0349R healthy. The measured viscous ignition boundary sits in
(0.0249, 0.0299]R with the margin still above ε at 0.0299R — the ε = 0.2
acceptance line therefore sits *above* the survival boundary, exactly the
headroom the pre-registration intended. Γ_implied at the survivor: 0.076
m²/s (vs Γ_v = 0.278); p-fit still needs a second survivor (s038v/sstab).

Figure CSVs regenerated; all three figures recompiled clean.

## 2026-08-06 ~22:35 MDT — wave-2 terminal: NT36 screen grid complete; sstab PASSES

`s038` (6:16), `s038v` (6:04), `sstab` (6:28) all **COMPLETED** their full
8 revs with case_metadata; harvested. Still running: `fid144` (784/1296
steps, current max_dtZ ~0.04, healthy) and the P3 gate 13058988 (step
~769, ~20 h). The full NT=36 Campaign E grid is now measured.

**sstab verdict (Ryan's "would σ_stab alone work as initializer" probe):
M_direct = 0.134 ≤ ε = 0.2 → PASS.** σ = σ_stab + 0.4% (0.0312R,
viscous) is stable *with margin headroom* by direct measurement.
Combined with s030v (0.0299R, M_direct = 0.295 > ε): **the ε = 0.2
acceptance boundary falls in (0.0299, 0.0312]R — σ_stab = 0.0311R lands
almost exactly on the measured ε-crossing.** The σ_stab initializer arm
is thereby validated far more tightly than the pre-registration hoped
(the prefactor band the P6 scope statement needs is ~±2%, not the
~2–3× resolution-price band ε was sized for).

**Survivor margin ladder (viscous, direct M):** 0.0299R → 0.295;
0.0312R → 0.134; 0.0381R → 0.112. Inviscid: 0.0381R → 0.141.

**Cross-validation, survivor side (now 4 points):** proxy/direct
0.21–0.36 — the proxy under-reads by ~3–5× uniformly on survivors
(vs 9–20× on deep-ignition windows). Consequence for P3 stands: hedge
proxy 0.0198 × (1/0.21–1/0.36) → direct-equivalent ≈ 0.055–0.094 < ε,
now bracketed by measured ratios on four survivors spanning the relevant
σ range; plus the directly-measured neighbors sstab (0.134) and s038v
(0.112) both PASS. The k=0 PASS is thus effectively certified pending
only the formal 13058988 completion write-up.

**p-fit (P1 addendum), first populated row:** viscous survivors give
p = 3.0 with Γ_implied 0.038–0.076 m²/s (≪ Γ_v = 0.278). Caveats
recorded: 3 points spanning only 0.030–0.038R, and the 0.0299R point sits
near the boundary where M may carry a transient component — p is
dominated by that point. Read against the pre-registered key (p≈2
sub-σ filament / p≈0 ambient): the steep slope + sub-Γ_v Γ_implied
suggests near-boundary margin growth steeper than the pure filament
model, not ambient-resolved behavior; defer interpretation to P5/P6
with the fit's small-range caveat attached. Inviscid row still has one
survivor (s038) — fit skipped.

Figure CSVs regenerated; all three figures recompiled clean.

## 2026-08-11 — final Campaign E harvest: fid144 AND the σ₀ hedge both ignited

Harvested the last two runs (both ended 2026-08-07, unharvested since).

**1. `scr_p019_fid144` (P4b): IGNITED — onset step 865 = 6.0 revs
(OOM at 13:52 h; pre-onset M_direct = 3.23).** This is the pre-registered
*ignite* branch of the P4b stability arm, recorded as written: the
long-horizon boundary does NOT scale as √Δt; σ_stab is a per-Δt
*shed-scale* criterion, not the full boundary; and per the P1-addendum
remedies ranking, Δt refinement is confirmed NOT the adjustment lever
(contraction-controlled boundary) — this goes to 020's gap-demonstration
evidence. Ignition at 6.0 physical revs continues the constant-revs
pattern of the `scr_ufdt` ladder (2.7–8 revs across every σ/Δt probed).
Proxy/direct on its pre-onset window: 0.007 — the proxy is essentially
blind at NT=144 (particle-switch bias worsens with 4× more shed events).
505 VTPs retained incl. final steps 915–917; the fidelity-arm histogram
scoring remains possible but must use a pre-onset window with an
explicit not-settled caveat (it never reached a settled wake).

**2. `p018_ufront_s035_visc` (P3 gate): IGNITED at step 996 = 27.6 revs —
the run was cancelled (by Ryan's UID, step 997/1080) at the exact step
the max_u tripwire fired (2773 m/s at 996, min_sr 0.033).** The full CSV
overturns the partial-harvest picture: proxy M pre-onset = 0.72, not the
0.0198 read at step ~587. **Consequence: the k=0 PASS is OVERTURNED at
certification — σ\* = σ₀ = 0.0356R is NOT certified; the procedure's
certify step did its job.** The 018 anchor "hedge 0.035R viscous healthy
in flight" is contradicted by completion — a finding for the 018 record.
Duration pattern: production 0.0299R died ~14 revs, 0.0349R at 27.6 revs.

**Restart assessment (Ryan's question): no cancelled run merits the
restart function.** Both `13058988` (s035) and `13060144` (n2_visc,
018-scope) were cancelled *mid-ignition* — n2_visc's last rows show
max_u 46→281→516 m/s over 3 steps at step 620 — so restarting would only
complete the blow-up. Their retained 10-step restart sets (987–996 and
611–620) are the ignition windows; both added to the VTK protect list
(pushed to cluster) as corpse evidence. The TIMEOUT runs (h0p125,
nt144s2) are 018 chain business, and the nt144 chain is being actively
managed on the cluster as of today (job 13134726 running) — not touched.

**P3 status after this:** k=0 FAIL at certification. The pre-registered
update rule with the hedge's own (proxy) M would give σ₁ = σ₀·√(M/ε) —
but the M estimate at production length, the screen-vs-production
duration bias (screens at 0.0312–0.0381R all passed ε while σ₀ = 0.0356R
died at production length), and whether ε itself must absorb more
duration bias are now the P3/P6 analysis questions. Dedicated session
needed; do not fold this into a quick verdict.

## 2026-08-11 — P4b scoring (fid144 vs s038v)

Scored per the pre-registered plan. Reduction: `scripts/vtp_sigma_reduce.py`
(stdlib VTP parser; also deployed to cluster) → per-step σ/σ_shed quantiles
+ pooled 60-bin histograms in `019_sigma_selection_procedure/p4b/`.
Windows (chosen at scoring time — the plan pre-registered the comparison,
not the step ranges): fid144 "mid" steps 700–709 (rev 4.9), "late"
840–849 (rev 5.9, ~15 steps pre-onset); reference `scr_p019_s038v`
314–323 (settled, last 10). Artifact-velocity gate passed in all windows
(max_u ≤ 36 m/s, window max_dtZ ≤ 0.11).

**Stability arm: IGNITED** (onset step 865 = 6.0 revs; pre-onset
M_direct = 3.23) — pre-registered outcome branch 2, recorded in the
2026-08-11 harvest section: σ_stab does not scale to NT=144 as a full
long-horizon boundary.

**Fidelity arm: the bulk-migration discriminator does NOT trigger.**
Pre-registered threshold: in-band settled median < 0.85 with reference
median ≥ 0.95. Measured medians: fid144 1.017 (mid) → 1.022 (late);
reference 0.984. The in-band bulk sits *at the shed scale* (marginally
above it — mean viscous growth), statistically indistinguishable from
the above-band reference; nothing migrated toward σ\*(Z). What moves is
only the low tail: fid144 p1 0.68 → 0.59, min 0.163 → 0.121 (= 0.14σ_eq)
while the reference's tail holds (p1 0.66, min 0.22 = 0.47σ_eq).

| quantile | fid144 mid | fid144 late | s038v ref |
| --- | --- | --- | --- |
| median | 1.017 | 1.022 | 0.984 |
| p5 | 0.885 | 0.857 | 0.784 |
| p1 | 0.681 | 0.592 | 0.665 |
| min | 0.163 | 0.121 | 0.221 |

**Interpretation (for P6):** the c_eq arm's feared failure mode in-band —
bulk core-size migration to the laminar equilibrium — is not what
happens at 1.18σ_eq. The observed in-band failure is *tail-driven
ignition*: a thin, deepening low-σ tail ignites while the bulk is
unmoved. The 2σ_eq arm's protective role is therefore about tail/
ignition control (via the Z̄(σ)-resolution capping of coarse cores), not
bulk viscous fidelity. Caveat attached as planned: fid144 never reached
a settled wake (ignited at 6 revs), so "settled median" means the
pre-onset plateau (median was static 1.017→1.022 over a full rev,
supporting the plateau reading). CT proxies (not verdict-bearing,
screens fail the aero gate): mean |CFx| 0.0793 (fid144 pre-onset) vs
0.0772 (s038v settled), +2.7%.

## 2026-08-11 — P3 verdict

**k=0: FAIL. σ\* ≠ σ₀.** The certification run at σ₀
(`p018_ufront_s035_visc`, 0.0349R production viscous — the config-matched
reuse probe) ignited at 27.6 revs. The screen-window readings that
underwrote the provisional PASS were not wrong, they were insufficient:
windowed proxy margin stayed at 0.004–0.05 for 25 revolutions (any
screen-length window PASSES ε), while min σ/σ_shed slid monotonically
0.80 → 0.13 (≈ −0.065/rev in ln, never arresting) into terminal runaway.
**Per-step M ≤ ε is necessary but not sufficient at production horizon;
the cumulative-contraction drift is the leading indicator** (P1 already
cited exactly this signature on the 0.0299R screen; the certification
step existed for this reason and did its job).

**Iterate branch applied (pre-registered).** The diverging-probe rule
says: σ_k is inside the boundary; take the last *clean* pre-ignition M,
else bisect upward. No M reading here is clean in the required sense —
the clean-window M (0.005–0.05) is drift-blind and would absurdly
prescribe σ < σ₀; the full-record M (0.716) is onset-saturated and
prescribes 0.067R. So the applicable branch is **bisect upward** in
(0.0349, 0.0381]R — and the upper endpoint is already
production-certified: `p018_L1_visc`, 0.0381R viscous, 40 revs, CT
0.0723 ± 0.0012. **The demonstration therefore terminates at σ\* =
0.0381R with zero new jobs.** (Optional refinement, not required: one
production-class probe at the interval midpoint 0.0365R, ~24 h, would
tighten σ\* by 4%; a screen there cannot decide it — screens pass ε
throughout this interval.)

**Duration bias, now measured on both sides:** screen-class ε-crossing
(0.0299, 0.0312]R vs production boundary (0.0349, 0.0381]R — the
production requirement is ~1.15–1.25× the screen crossing in σ. This is
the quantitative content of "ε must absorb probe-vs-production duration
bias" (P1), and the honest statement is that ε alone cannot absorb it:
the drift criterion (contraction arrest, d ln min_sr/dt → 0 over the
settle window) is the production-horizon observable. Recorded here as a
P6 procedure-text finding — the pre-registered algorithm is NOT patched
retroactively.

## P5 — figures and captions (2026-08-11)

Four standalone TikZ figures in `019_sigma_selection_procedure/figures/`
(each `<name>.tex` + same-named CSV data dir, pdflatex-verified,
regenerated from the pipeline; build artifacts gitignored). Encoding
conventions shared across figures: color = viscosity (blue inviscid /
vermillion viscous, CVD-validated pair); marker shape = startup class
(circle screen / square production); fill = outcome (filled ignited /
open survivor); edge weight = M estimator (bold direct / thin proxy).
Draft captions (self-contained, per the polish contract):

- **fig_margin_curve** — Stability margin $M=\max\Delta tZ$ vs core size
  $\sigma/R$ (log–log), all probes and production runs, with the failure
  thresholds $\Delta tZ{=}1,2$, the pre-registered target
  $\varepsilon=0.2$, the initializer scales $\sigma_{eq}$,
  $\sigma_{stab}$, $\sigma_0=2\sigma_{eq}$, and the survivor-only fit
  $M\propto\sigma^{-3.0}$ (viscous, NT=36). Ignited runs plot their
  pre-onset $M$ (filled): every point left of $\sigma_{stab}$ exceeds
  the $\Delta tZ$ thresholds, while survivors cluster at or below
  $\varepsilon$ — the measured $\varepsilon$-crossing falls in
  $(0.0299,0.0312]R$, on top of $\sigma_{stab}=0.0311R$.
- **fig_regime_map** — Outcome over $(\sigma/R,\ \text{viscosity})$ with
  revolutions-to-ignition (log y): ignited runs at their onset time,
  survivors in the shaded "no ignition" band. Time-to-ignition rises
  toward the boundary at $\sigma_{stab}$; the production squares (L2 at
  28.8 revs, 0.0299R at 13.5, $\sigma_0$-hedge at 27.6) show the
  longer-horizon boundary sitting above the screen boundary.
- **fig_iterates** — The procedure trace: $\sigma_k/R$ per iteration
  against the initializer scales (top) and $M_k$ against
  $\varepsilon$ (bottom); k=0 at $\sigma_0=0.0356R$ FAILs at
  certification, terminating via the bisect-upward branch at
  $\sigma^*=0.0381R$ (existing production certification).
- **fig_p4b_hist** — $\sigma/\sigma_{\mathrm{shed}}$ distributions
  (log y): in-band NT=144 run ($\sigma=1.18\sigma_{eq}$, pre-onset, two
  windows one rev apart) vs the above-band settled reference
  ($2.14\sigma_{eq}$). The bulk peaks coincide at the shed scale
  (medians 1.02 vs 0.98 — no migration toward the laminar equilibrium);
  the in-band run instead grows a thin low-$\sigma$ tail (min
  $0.14\sigma_{eq}$) that subsequently ignites — the in-band failure
  mode is tail-driven ignition, not bulk fidelity loss.

Display-label → run/job mapping: `figures/provenance_appendix.csv`
(pipeline-emitted).

## P6 — Close-out (2026-08-11)

### Procedure-text reconciliation (P1 pre-registration vs demonstration)

The P1 algorithm and pre-registered choices are unchanged; the
demonstration required four recorded deviations/refinements (all Logged
at the time, none silently patched):

1. **Survivor-only p-fit** (dated section 2026-08-06 ~16:15): ignited
   runs' pre-onset M are divergence transients (1.1–28 ≫ 2); fitting
   them yields unphysical p. Fit restricted to survivors; ignited points
   reported unfitted. Result: p = 3.0 (viscous, 3 points, narrow range
   0.030–0.038R — quoted with that caveat), Γ_implied 0.038–0.076 m²/s
   ≪ Γ_v: steeper than the pure sub-σ filament model (p=2), not
   ambient (p=0); consistent with near-boundary margin growth.
2. **Proxy estimator demoted to unusable for verdicts** (P1 already
   ranked it a fallback): measured proxy/direct = 0.21–0.36 on
   survivors, 0.05–0.11 on deep ignitions, 0.007 at NT=144. Every
   pre-Campaign-E margin is a lower bound only.
3. **Diverging-certification branch interpretation** (P3 verdict): when
   a *certification* (not a probe) diverges, no M reading is "clean" —
   per-step M is drift-blind, full-record M is onset-saturated — so the
   bisect-upward branch is the applicable rule; it terminated on the
   existing certified rung (σ\* = 0.0381R, L1).
4. **Missing production-horizon criterion identified** (the substantive
   finding): per-step M ≤ ε passed for 25 revs while min σ/σ_shed slid
   0.80 → 0.13 without arresting. **Recommended v1.1 addition (for
   Ryan's sign-off, not retro-applied): certification requires BOTH
   M ≤ ε and contraction arrest, d ln(min σ/σ_shed)/dt → 0 over the
   settle window.** Screen-class probes select; only the drift
   criterion certifies at long horizon. Duration bias measured:
   production boundary (0.0349, 0.0381]R ≈ 1.15–1.25× the screen
   ε-crossing (0.0299, 0.0312]R.

### Validated scope (final)

One flow family: small hover rotor (DJI 9443 class), explicit-Euler VPM
with CoreSpreading viscosity, NT = 36–144, screen and production startup
classes. Validated: σ_stab = √(Γ_vΔt/2π) predicts the *screen-class,
same-Δt* ε-crossing to ~2% (0.0311R vs measured (0.0299, 0.0312]R).
NOT validated / bounded: σ_stab does not extrapolate across Δt as a
long-horizon boundary (NT=144 at 1.35·σ_stab(Δt/4) ignited at 6 revs —
constant-revs contraction governs, Δt is not a mitigation lever); the
bulk-fidelity reading of the c_eq·σ_eq arm (in-band failure is
tail-driven ignition, bulk stays at shed scale). Falsifier elsewhere: a
flow whose screen-class ignition boundary departs from √(Γ_vΔt/2π) by
more than the ~±5% prefactor band measured here at matched protocol.

### Acceptance criteria status

1. Standalone algorithm — **met** (P1 + addendum).
2. Single-entry-point reproducibility — **met**
   (`scripts/p019_sigma_procedure.py`; P4b reduction documented in
   `scripts/vtp_sigma_reduce.py` + `p4b/` CSVs; every verdict number
   traces to in-repo CSVs with provenance).
3. Three plots with required annotations — **met** (four authored incl.
   the P4b histogram; margin plot carries σ₀/ε/boundary/scales/slope;
   regime map has both sides of the boundary with 9 measured ignitions
   and per-cell provenance via appendix CSV).
4. Validated-scope statement — **met** (above).
5. Notebook entry — **OPEN: requires Ryan's approval** (explicitly
   deferred by Ryan 2026-08-06; proposal pending).

Item-level approvals (header checkboxes) await Ryan.

## Log

- 2026-08-12 — sign-off session (Ryan, interactive): notebook entry
  written (journals/20260803.md § 20260812, "σ-Selection Procedure
  Demonstrated (019)" — context + verdict table + margin-curve and P4b
  figures + conclusions; figures copied to
  notebooks/img/20260812_019_sigma_procedure/ per image policy);
  drift-arrest v1.1 amendment APPROVED and added to P1 (subsection +
  step-4 pointer, v1.0 untouched); Technical Completion checkbox marked
  (header + INDEX). Clear-context and user approvals remain open.
- 2026-08-11 — P5 polish (fig_p4b_hist authored; p-fit slope annotation
  added to margin curve; captions drafted) and P6 close-out written
  (reconciliation: 4 recorded deviations incl. the drift-arrest v1.1
  recommendation; validated scope finalized; acceptance 1–4 met, 5 =
  notebook pending Ryan). INDEX updated.
- 2026-08-11 — P4b scored (fidelity discriminator NOT triggered: bulk
  static at shed scale, medians 1.02 vs 0.98; failure mode is tail-driven
  ignition) and P3 verdict written (k=0 FAIL; bisect-upward branch →
  σ\* = 0.0381R via the existing L1 certification; duration bias
  measured ~1.2× in σ; drift-arrest identified as the missing
  production-horizon criterion — P6 item, not a retro-patch).
- 2026-08-11 — final harvest (see dated section above): fid144 ignited
  (P4b ignite branch, Δt not a lever, proxy blind at NT144);
  σ₀ hedge ignited at 27.6 revs ⇒ k=0 OVERTURNED at certification;
  no cancelled run merits restart (both died mid-ignition; restart sets
  = ignition windows, protected). Pipeline + figures regenerated.
- 2026-08-06 ~22:35 — wave-2 harvest (see dated section above): NT36 grid
  complete; sstab PASS (M 0.134) ⇒ ε-crossing measured in
  (0.0299, 0.0312]R, σ_stab validated at ~±2%; survivor proxy/direct
  0.21–0.36; viscous p-fit 3.0 (small-range caveat). Remaining: fid144
  (~01:00–03:00), P3 gate (Aug 7 afternoon).
- 2026-08-06 ~18:20 — wave-1 harvest (see dated section above): s030v
  survivor M_direct 0.295 > ε (proxy would have passed — survivor
  proxy/direct 0.27); s025 pair ignited; registry updated; figures
  recompiled. P3-verdict note: 13058988 lacks max_dtZ — certification
  route via sstab + quoted proxy-ratio band.
- 2026-08-06 ~16:45 — first Campaign E harvest session (see dated section
  above): 3 ignition corpses harvested with banner provenance; pipeline
  extended (REGISTRY +10, direct max_dtZ, xval, survivor-only p-fit +
  Γ_implied, plot-ready figure CSVs); survivor-only p-fit ruling recorded
  as a deviation; three P5 TikZ figures authored
  (`figures/fig_{margin_curve,regime_map,iterates}.tex`, pdflatex-verified,
  auto-refresh from pipeline CSVs; build artifacts gitignored); protect
  list += fid144, s038v (pushed to cluster). Remaining: 7 screens + P3
  gate still RUNNING at 16:40 — steps 4 (P3 verdict) and 5 (P4b scoring)
  blocked on them.
- 2026-08-06 (later) — kickoff. P1 pre-registration section written (ε = 0.2,
  M = max ΔtZ direct-or-proxy with labeled method, screen-class probes,
  Z̄ recomputed). 018 mechanism doc corrected in step: 13058534
  (`p018_ufront_n1_visc`, 0.0299R viscous) DIED — refutes "1.7σ_eq safe" at
  its boundary case and lands the viscous production outcome on the side
  σ_stab predicts (0.0299 < 0.031 < 0.0349-healthy). Ryan decisions: full P4
  grid; add `max_dtZ` monitor column now; build pipeline now, gate P3 on
  13058988. HPC submission authorized for this session. Plan:
  `~/.claude/plans/work-on-019-per-recursive-nest.md`.
- 2026-08-06 — staged (this file); no work started. Expected cheap: if the
  σ₀ = 0.035R probe meets ε immediately (likely, given the in-flight hedge's
  health), the demonstration is "initialize → verify stable → plot margin
  curve from probes at 2–3 σ values", ~3 screen-class jobs total.
