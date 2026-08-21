# Phase 2 — Physical Das Length

**Objective:** find the smallest physical first-wake-row offset that still
converges. Per 014's rule and the campaign reconciliation: "converges" ≡ on
the log plateau (successive-doubling |ΔCT̄| ≤ 0.5%, consistent with the
~0.2%/doubling wing law); "smallest" = the plateau's lower edge. Below-plateau
points are the amplified regime (the +36.8% cliff) and are non-converged by
definition. Expected Das\* = 0.41c (η=1.0 frozen).

## Cases (carrier: NT=36, σ/c=0.68, N=4, velocity; B0 supplies η=1.0)

| tag | η | Das/c(0.75R) | floor clamp (span frac) | time |
| --- | --- | --- | --- | --- |
| `p018_das0p5` | 0.5 | 0.205 | ~20% (report measured) | 24 h |
| `p018_b0` (done, P1) | 1.0 | 0.41 | ~4% | — |
| `p018_das2p0` | 2.0 | 0.82 | ~0% | 24 h |
| `p018_das4p0` | 4.0 | 1.64 | 0% | 24 h |
| `p018_das0p5_floorhalf` (optional) | 0.5, DAS_MIN_DISPLACEMENT_R=0.005 | 0.205 | halved | 24 h |

Run the floorhalf case only if η=0.5 is within 1% of η=1.0 (i.e. the smallest
rung is a plateau candidate and floor contamination must be bounded).
Extend any case by restart if its mean hasn't settled per M1.

## Mandatory side job — FMM |Das|-radius discriminator (017 T4)

Every prior rotor η point used FMM, whose panel radius is inflated by |Das|
(`src/FLOWPanel_liftingbody.jl:996`). Before publishing any Das claim: run
`DirectBackend` vs `FastMultipoleBackend` on `RHPC_MESH=40_40`, NT=6, ~6 revs,
at η ∈ {0.5, 4.0}, locally (≤6 threads) or one 1-h cluster job. Pass = the
Direct-vs-FMM CT delta is η-independent (≤0.1% variation). Fail ⇒ the η ladder
carries an FMM artifact — rerun the two extreme ladder points on
DirectBackend-compatible settings before closing the phase.

## Decision

Das\* = smallest ladder value with |CT̄(Das) − CT̄(2·Das)| ≤ 0.5% and
ε_Γ ≤ 1%. Report the plateau slope (%/doubling) → error-budget term 2. Legacy
η ladder (0.06148/0.06942/0.07133/0.07190) is corroborating shape only.
If NO plateau exists at this σ, do not iterate here — proceed with η=1.0
provisionally and let Phase 4's σ refinement re-test the ladder's top two
points at σ=L1 (Das sensitivity may shrink with σ; record as a finding).

## Exit criteria

Das\* selected with plateau evidence + discriminator pass; ladder + Γ overlays
in the ledger; floor-clamp fractions recorded per case.

## Log

- 2026-07-31 — **T4 FMM discriminator PASSED.** 40_40, NT=6, velocity,
  SFS_OFF=true both arms, matched physical Das via η ∈ {0.083, 0.667}
  (= Das/c 0.20 / 1.64 at NT=6). Tail-mean CT: Das=0.20c fmm 0.056120 vs
  direct 0.056122 (−0.0033%); Das=1.64c fmm 0.056807 vs direct 0.056809
  (−0.0019%). Backend delta η-variation 0.0014% ≪ 0.1% ⇒ the FMM
  |Das|-panel-radius coupling does NOT contaminate the Das axis. CSVs kept at
  `data/p018_t4_discriminator/`. Two execution caveats discovered: (1) η must
  be scaled by 36/NT to hold physical Das in coarse-NT smokes (first attempt
  at η=4, NT=6 ⇒ Das≈9.8c diverged to CF≈430 — physically expected, discarded);
  (2) the SFS *direct* code path has a small-N bug
  (`FLOWVPM_subfilterscale_models.jl:54` `Estr_direct_multithreaded` builds a
  zero-step thread range when np < nthreads) — discriminator ran SFS-off in
  BOTH arms; production (FMM path) unaffected. Driver gained
  `RHPC_BACKEND=fmm|direct` (default fmm, no behavior change), synced +
  md5-verified to cluster.
- 2026-07-31 — Ladder jobs submitted: 12993802 das0p5, 12993803 das2p0,
  12993804 das4p0 (η=1.0 supplied by p018_b0, 12993801).
- 2026-07-31 — **First rung harvested: `p018_das0p5` (12993802) COMPLETED in
  12.5 h.** CT̄(revs 10–19) = **0.07006**, per-rev std ±0.00031, half-window
  drift −0.33% (oscillatory, no monotone trend). Metadata-verified SigmaOverlap
  overlap=2.0, nwakerows=4. Startup transient fully dead by rev ~9; per-rev
  peak-to-peak collapses to ~0.0003 in the settled window — the campaign
  config's scatter is ~5× tighter than every legacy overlap_pps run. Window is
  10 settled revs (run is 20 total): short of M1's ≥15 — extend by restart if
  the plateau-edge decision turns out to hinge on this rung. Ladder analysis
  deferred until b0/das2p0/das4p0 land. CSVs + circulation monitor at
  `data/p018_das0p5/`.
- 2026-07-31 — **`p018_das2p0` (12993803) COMPLETED in 14.9 h.** CT̄(revs
  10–19) = **0.07230**, per-rev std ±0.00027, half-window drift −0.16% (best
  settled of the three so far). η=2.0 metadata-verified. Running ladder:
  Das 0.205c → 0.07006, 0.41c → 0.07187 (+2.58%), 0.82c → 0.07230 (+0.60%).
  The 0.41→0.82c doubling delta sits just above the 0.5% convergence
  threshold and ~3× the 014 wing log-law slope (~0.2%/doubling) — plateau
  judgment (and whether Das\* = 0.41c or 0.82c) waits on das4p0's 0.82→1.64c
  delta, the M2 ε_Γ comparison, and b0's extended (s2) window since b0's mean
  was still drifting −0.38%/half-window. CSVs + circulation monitor at
  `data/p018_das2p0/`.
- 2026-07-31 — **`p018_das4p0` (12993804) COMPLETED (15.5 h) — LADDER ANALYSIS;
  PHASE 2 RESOLVED VIA THE PRE-REGISTERED NO-PLATEAU BRANCH.**

  M1 ladder (CT̄ revs 10–19, per-rev std):

  | Das/c | case | CT̄ | std | doubling Δ |
  | --- | --- | --- | --- | --- |
  | 0.205 | das0p5 | 0.07006 | ±0.00031 | — |
  | 0.41 | b0 | 0.07187 | ±0.00034 | +2.58% |
  | 0.82 | das2p0 | 0.07230 | ±0.00027 | +0.60% |
  | 1.64 | das4p0 | 0.07084 | ±0.00044 | **−2.02%** |

  M2 (ε_Γ from `circulation_te`, blades averaged, same window, 0.3≤r/R≤0.95,
  normalized by max|Γ̄|): das0p5→b0 max 8.97% / RMS 3.25%; b0→das2p0 3.11% /
  1.52%; das2p0→das4p0 5.32% / 2.78%. **No successive pair passes M1
  (|Δ|≤0.5%) AND M2 (ε_Γ,max≤1%).** Top rung is non-monotone (−2.02%),
  unlike the monotone-saturating legacy ladder; das4p0's case metadata is
  otherwise identical to das2p0's and its window is well settled (drift
  −0.14%), so the drop is a real configuration effect. The peak-Γ̄ station
  also migrates inboard with Das (r/R 0.740 → 0.695 → 0.583 → 0.583).

  **Decision (per the pre-registered branch "if NO plateau exists at this σ"):
  Das\* = 0.41c (η=1.0) PROVISIONAL; do not iterate here; Phase 4's σ=L1
  re-tests the ladder's top two points.** Plateau slope / error-budget term 2
  deferred to that re-test (no defensible slope exists at σ/c=0.68).

  Findings & caveats recorded:
  1. **Handoff-distance confound:** at fixed N=4 the sheet extent d≈N·Das spans
     1.2σ → 9.6σ (in σ=0.68c units) across the ladder — the two lowest rungs
     violate the d≳4σ criterion ruling 4 assumed. The Das axis at fixed N is
     inherently coupled to handoff adequacy ⇒ **Phase 5's N=8 spot-check is
     now decision-critical**, and should run at Das\*=0.41c where d only = 2.4σ
     (N=8 ⇒ 4.9σ).
  2. Floor clamps (measured from warm-start metadata |Das| vs 0.01R):
     das0p5 17.1% of stations, b0 2.4%, top rungs ≲2% by monotonicity —
     matches the phase table's predictions; das0p5's rung carries real floor
     contamination (but it is below-plateau anyway).
  3. `circulation_slice` is all-zero in every campaign monitor CSV. Originally
     read as "the estimator never ran (needs a driver flag)" — **superseded
     2026-08-03**: the estimator was found to be **telescoping, i.e. a
     theorem-zero** — the loop it integrates closes on itself, so its exact
     value is zero and no driver flag would have changed that. A **Kelvin-loop
     replacement was in progress** as of that date (completion state not
     recorded elsewhere in these files — verify before relying on it). M2
     therefore rests on `circulation_te` alone; do NOT spend effort "enabling"
     the slice estimator. Any instruction elsewhere in this item to report or
     enable the slice estimator as a cross-check is superseded by this note
     (see `decision_rules.md` §M2, `phase_10_mesh_spanwise.md` §Decision).
  4. b0's own window drift (−0.38%) is the largest of the four rungs; its s2
     extension (13007210, running) refines the η=1.0 number before Phase 3
     comparisons.

  T4 discriminator already PASSED (above). Phase 2 exit criteria met under the
  no-plateau branch. **Submitted on resolution (2026-07-31): 13011373
  `p018_nt72` (Phase 3, η=2.0=2η\*, 48 h), 13011374 `p018_L1` (Phase 4 cold,
  36 h — single-segment deviation from the 2×24 h plan, justified by the
  recalibrated ~30 h estimate), 13011375 `p018_L1_warm` (warm pilot from b0
  step 719, SETTLE=30 ⇒ ~8 flush revs + window, 36 h).**

## 2026-08-03 — Γ̄(r/R) explains the non-monotone top rung: TWO opposed mechanisms

Figure: `data/p018_gamma_ladders/gamma_das.png` (`scripts/p018_plot_gamma.py das`).

### Successive-rung M2 (matched revs 10–19)

| step | ε_Γ max | ε_Γ RMS |
| --- | --- | --- |
| 0.205 → 0.41c | **8.963%** | 3.243% |
| 0.41 → 0.82c | **3.112%** | 1.520% |
| 0.82 → 1.64c | **5.320%** | 2.783% |

Non-monotone in Γ̄ just as in CT̄, with the minimum at the 0.41 → 0.82c step.
**No rung pair passes M2**, consistent with Phase 2's no-plateau branch — and
now corroborated on the co-equal observable rather than on CT̄ alone.

### The mechanism: two errors of opposite sign in opposite places

The distribution plot separates what CT̄ had summed into one confusing number:

- **Small Das (0.205c): an INBOARD deficit** — −8.9% at r/R 0.31, decaying
  monotonically to ≈ −1% by r/R 0.6, roughly flat outboard.
- **Large Das (1.64c): an OUTBOARD deficit** — −6.6% at r/R 0.84, while
  *inboard it has saturated*: 0.82c and 1.64c lie on top of each other for
  r/R < 0.5 (both ≈ +3 to +4% above B0).

So as Das grows, the inboard gain **saturates** while the outboard loss
**keeps growing**. Their sum is a CT̄ that rises, peaks, and turns over — which
is exactly the observed ladder (0.07006 → 0.07187 → 0.07230 → 0.07084).
**Phase 2's non-monotone top rung is no longer unexplained.**

Physical reading of each half:

- The inboard deficit at small Das is a **clearance** effect, shared with the
  N=1 failure — see phase_05, where the same inboard signature appears for the
  other small-`d` case. Both reduce d = N·Das.
- The outboard deficit at large Das is the **rigid-sheet** limit: the first
  wake row is a straight frozen sheet, and at 1.64c it extends far enough that
  it misrepresents the rolled-up tip vortex, precisely where wake curvature
  matters most.

These two mechanisms bound Das from below and above respectively, which is a
more useful statement than "no plateau": the ladder is not failing to
converge, it is being squeezed between two distinct errors. **Das\* = 0.41c
remains provisional, but 0.82c is now a live competitor** — it has the smallest
successive ε_Γ of any step and sits between the two mechanisms. Re-test both at
σ = L1 as already planned (Phase 4's deferred re-test (a)); the σ refinement
may move the outboard mechanism substantially since it is tip-vortex-mediated.

**Update 2026-08-04 (Das×N matrix, see phase_05):** the collapse test refutes
d/σ as the governing group — **both lobes are functions of Das itself** (the
outboard deficit at Das=1.64c is nearly N-independent; the inboard deficit at
Das=0.205c occurs at a d/σ where other cells are clean), i.e. the
"clearance" reading of the small-Das rung is corrected to an
attachment-length effect. Also **0.82c's competitor status is weakened**: at
Das=0.82c, N=2 vs N=4 fails M1+M2 (−0.62% CT, ε_Γ 2.2%), whereas at 0.41c the
N=2≡N=4 null is clean (0.35%).

## 2026-08-04 (b) — ROOT CAUSE OF THE NO-PLATEAU LADDER + chord-proportional Das (Ryan-directed)

**The η parameterization is structurally unable to converge Γ̄(r/R).** `Das` is
set per station as η·Δt·|V_te| with V_te = Ωr, so Das ∝ r — but the admissible
band 014 established is in **local-chord units** (≈0.25–1.5c), and the DJI
chord falls outboard. With the actual mesh chords
(`scripts/p018_mesh_profile.py show`), Das/c_local at the ladder's rungs:

| r/R | η=0.5 | η=1 (B0) | η=2 | η=4 |
| --- | --- | --- | --- | --- |
| 0.31 | 0.04c | 0.09c | 0.18c | 0.35c |
| 0.75 | 0.21c | 0.41c | 0.82c | 1.64c |
| 0.95 | 0.33c | 0.66c | 1.32c | 2.64c |

Putting r/R=0.31 on the plateau (≥0.25c) needs η ≥ 2.8; keeping r/R=0.95
under 1.5c needs η ≤ 2.3 — **the admissible windows do not overlap.** The
ladder was sliding a wedge through a fixed band, which reproduces every
Phase-2 observation: inboard deficit at small η, outboard deficit at large η,
inboard "saturation" between 0.82c and 1.64c (r/R 0.31 crosses 0.25c exactly
there), and the non-monotone CT̄ sum. The matrix result above ("Das owns both
lobes") is the same statement: the two lobes are the two edges of the local
band. Note 014's draft rule already said "fraction of local chord / TE
travel" — the implementation froze the ∝r kinematic length and the
local-chord half was never built.

**Implemented (Ryan 2026-08-04): `DAS_CHORD_FRACTION` = κ ⇒ |Das|_j = κ ·
local chord at station j**, direction along the local kinematic TE tangent.
src: `set_Das_station_lengths` kwarg in `initialize_Das!`
(`src/FLOWPanel_simulate.jl`; mutually exclusive with the η path; unit-tested,
7/7). Driver: `station_chords` computes per-station chords from the body
geometry (vertex-based, matching `_kinematic_velocity_te!`'s station
convention: station j = edge j's nib node, station n+1 = last nia); banner
prints `das_chord:`; metadata records `das_chord_fraction`. dt-independent by
construction (no more η∝NT bookkeeping), freeze factor and floor irrelevant.
Launcher arms `p018_dasc0p25|0p41|0p82` (B0 carrier; κ=0.41 matches B0's Das
at 0.75R exactly). Disclosure: rows 2+ keep per-step-travel lengths, so the
row-1/row-2 length mismatch becomes spanwise-varying (the η scheme already
carried a uniform 2.5× mismatch via the freeze factor; 014 established the
partition tiles cleanly at any Das).

**Pre-registered κ-ladder predictions:** (1) at κ=0.41 vs B0 (same Das at
0.75R), the inboard lobe should move — inboard stations get up to ~4.6× more
Das — while the mid-outboard stays close; (2) successive-κ ε_Γ should be
dominated by a uniform (log-law) shift, not the two opposed lobes; (3) if the
κ ladder passes M1+M2 where the η ladder had no passing pair, the Das axis is
convergent under the chord parameterization and the campaign's Das\* becomes
κ\* with the wing-law ~0.2%/doubling model-form term.

**Pre-registered FALLBACK (Ryan authorization 2026-08-04: "if it doesn't
work... you have my permission to try the most likely thing"):** if the κ
ladder fails (no successive pair passes M1+M2, or the two opposed lobes
persist), the next approach — chosen NOW, before the data — is a
**GreenReconstruction Das pair on the rotor**: `RHPC_FORMULATION=green` at
κ=0.41 and κ=0.82 (B0 carrier otherwise; existing knob, Phase-0 `p2e_green`
precedent). Rationale: 015 Battery II measured that swapping only the
wake→body transfer from `VelocityThroughSources` to Green removes **93% of
the wing's Das differential and all of its η-dependence**, i.e. much of the
Das sensitivity is a transfer-scheme artifact. Decision rule: if
|ΔCT̄(κ 0.41→0.82)| and ε_Γ under Green shrink to ≤ the M1/M2 thresholds
where the velocity formulation failed them, Green becomes the production
formulation (this folds Phase 7 into Phase 2) and the Das axis is closed
with the Green pair's residual as its budget term; 016 remains the σ-axis
contingency (ruling 7), second in line for a Das-axis failure.

## 2026-08-05 — κ LADDER HARVESTED: FAIL. Fallback executed (Green pair already in flight)

Jobs 13050924/25/26, all COMPLETED exit 0 (11.9 / 13.5 / 15.4 h), 20-rev
histories. Scored with `scripts/p018_analyze.py` on **matched 10–19 windows**
(self-check `m1 p018_b0 --revs 15 29` → 0.071729 first). Numbers in `ledger.md`
§2026-08-05 (b).

| κ (\|Das\|=κ·c_local) | CT̄ (10–19) | ΔCT̄ vs previous | ε_Γ vs previous |
| --- | --- | --- | --- |
| 0.25 | 0.071571 | — | — |
| 0.41 | 0.072154 | **+0.81% (M1 FAIL)** | 0.952% (M2 PASS) |
| 0.82 | 0.073455 | **+1.80% (M1 FAIL)** | **5.219% (M2 FAIL)** |

**Verdict against the three pre-registered predictions:**

1. *Inboard lobe moves at κ=0.41 vs B0* — **CONFIRMED, and it is small**:
   ΔCT̄ = +0.40%, ε_Γ = 1.002% max / 0.318% RMS. Replacing the r-wedge with the
   chord-wedge at fixed Das(0.75R) redistributes ~1% of max Γ̄ — real, but an
   order below the 8.78% σ-axis redistribution.
2. *Successive ε_Γ dominated by a uniform log-law shift* — **REFUTED** at the
   top step: 5.2% max against 1.95% RMS is lobe-shaped, not uniform.
3. *Ladder passes M1+M2 ⇒ Das axis convergent under κ* — **FAILED.** No
   successive pair passes M1, and the deltas **grow** (+0.81 → +1.80%): this is
   not a ladder approaching a limit. Fixing the local-chord band removed the
   non-monotonicity (the sequence is now monotone) but did not produce a
   plateau.

⇒ **The pre-registered Ryan-authorized fallback applies.** It required no new
submission: the Green κ pair (`p018_green_dasc0p41` 13051802, `p018_green_dasc0p82`
13051803) was launched 2026-08-05 (c) as a latency hedge and is running. Decision
rule unchanged (above): if Green's κ 0.41→0.82 step passes M1+M2 where velocity
failed, Green becomes production and Phase 7 folds into Phase 2.

### Code finding that reframes the Das/σ relationship (record before the Green verdict)

Read while planning the follow-up, from the shedding path itself:

- Trailing particles are deposited along the **streamwise** segment of the
  handed-off wake row, `r1_le → r1_te` (`src/FLOWPanel_wake.jl:2099–2101`), and
  `OverlapPPS` sets `σ = |r2−r1|·overlap/pps` (`:678–683`). So
  **σ_trailing ∝ Ω·r·Δt** — the *same* ∝r wedge that motivated κ. (Unsteady
  particles are the other family, shed along the spanwise TE edge `r1_te→r2_te`,
  σ ≈ spanwise-uniform.)
- Therefore, under the **old η prescription** `Das = η·Δt·|V_te|`,
  **Das/σ_trailing = η·pps/overlap** — constant along span *and* dt-independent.
  B0 (η=1, pps 4, overlap 2) sits at Das/σ = 2, and with N=4 doublet rows the
  body→first-particle clearance is ≈ N·Das ≈ 8σ, i.e. **2–3× above the Phase-12
  C1 admissible clearance d/σ\* ≈ 2.6–3.4**.

So the two admissibility criteria are in direct conflict: **κ·c_local satisfies
014's local-chord band and breaks the σ-clearance constancy; η satisfies the
σ-clearance constancy and breaks the chord band.** This does not by itself
explain the κ failure — C1 says B0 already clears the kernel threshold, and the
Das×N matrix found the Γ̄ lobes N-independent (⇒ independent of total clearance
at fixed Das) — both of which argue *against* clearance being the driver and in
favour of the sheet-length / transfer-scheme reading that the Green fallback
tests. Hence: **no σ-referenced ladder is being launched blind.** Next
discriminator is free and offline (extend `examples/particle_surface_regularization_diag.jl`
to measure per-station σ_shed(r/R) and the actual body→nearest-particle
clearance in local-σ units, accounting for the N rows covering [TE, TE+N·Das],
on B0 and the three κ runs, then correlate with the measured Γ̄ lobes). If that
diagnostic *does* show a clearance violation, the ladder to run is
λ ≡ Das/σ ∈ {2, 3.4, 5} at N=2 — implementable as
`DAS_ETA_KINEMATIC = λ·overlap/pps` at fixed overlap/pps, **no new code**.

## 2026-08-05 (c) — Green κ=0.41 arm harvested: formulation delta small at this rung

`p018_green_dasc0p41` (13051802) COMPLETED (exit 0, 720 steps / 20 revs).
Metadata-verified: `formulation="green"`, κ=0.41 (`das_chord_fraction`),
OVERLAP 2.0 / PPS 4, N=4, settle 12 — exact partner of `p018_dasc0p41`.
Matched 10–19 windows (`p018_analyze.py`):

| | CT̄ (10–19) | 95% CI | drift |
| --- | --- | --- | --- |
| `p018_dasc0p41` (velocity) | 0.072154 | [0.072037, 0.072263] | 0.27% |
| `p018_green_dasc0p41` | **0.071921** | [0.071806, 0.072003] | 0.27% |

ΔCT̄(green − velocity) = **−0.32%**; M2 ε_Γ,max **0.716% / RMS 0.453% — PASS**.
At κ=0.41 the two wake→body transfer schemes agree to within the M1
threshold — an early Phase-7 Green-Δ data point. **No verdict yet on the
pre-registered fallback question**: that requires the κ=0.82 arm
(13051803, ~95% done) to see whether Green's successive κ delta
(0.41→0.82) stays bounded where velocity's grew (+0.99% on top of +0.81%).
Interpretation caveat (phase_13 §(e)) stands: the velocity κ ladder was
collinear with a clearance ladder, so even a growing Green delta would not
cleanly indict chord-Das itself.

## 2026-08-05 (d) — Green κ pair COMPLETE: Green does NOT rescue the κ ladder

`p018_green_dasc0p82` (13051803) COMPLETED (exit 0, 15:07 h, 720 steps).
Metadata-verified (green, κ=0.82, OVERLAP 2.0 / PPS 4, settle 12). Matched
10–19 windows:

| | κ=0.41 | κ=0.82 | successive Δ |
| --- | --- | --- | --- |
| velocity | 0.072154 | 0.073455 | **+1.80%** |
| green | 0.071921 | **0.073046** CI [0.072804, 0.073187] | **+1.56%** |

**Pre-registered fallback decision rule RESOLVED negative:** Green shows the
same growing κ delta (+1.56% vs velocity's +1.80%, both ≫ M1's 0.5%) ⇒ the
κ-ladder behavior is **formulation-independent**, so Green does NOT become
the production formulation off this test (the rule required Green to pass
M1/M2 where velocity failed). The green−velocity delta itself stays small and
roughly constant (−0.32% at κ0.41, −0.56% at κ0.82; M2 0.72%/0.41% PASS both)
— banked as early Phase-7 Green-Δ data points.

**Interpretation:** a transfer-scheme-independent, clearance-correlated κ
response is exactly what the phase_13 §(e) reinterpretation predicts (the κ
ladder was collinear with a d/σ clearance ladder wholly inside the
inadmissible region). The discriminating experiment is now the
uniform-clearance pair `p018_ufront_n1/n2` (phase_13 §7, being submitted),
which moves Das while HOLDING d/σ uniform at an admissible value — not
another κ rung. The Das axis stays OPEN with Das\* = 0.41c provisional.
