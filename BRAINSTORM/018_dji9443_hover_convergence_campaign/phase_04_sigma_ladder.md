# Phase 4 — Particle σ Ladder (h refined faster than σ)

**Objective:** classical VPM convergence in σ under `SigmaOverlap` with
σ/h > 1.3 everywhere and h ~ σ^q, q > 1 (OVERLAP increases down the ladder;
017 theory). This is the axis most likely to need the contingency chain.

## Ladder (carrier: NT=36, Das\*, N=4, velocity)

| tag | OVERLAP | PPS | σ/c(0.75R) | σ/R | h/c | q vs prev | MERGE_R_FACTOR | N_p est. | wall est. |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `p018_b0` (L0, done) | 2.0 | 4 | 0.68 | 0.0873 | 0.34 | — | 0.0120 | ~92k | ~20 h |
| `p018_L1` | 2.4 | 11 | 0.297 | 0.0381 | 0.124 | ~1.22 | 0.0053 | ~210k | ~45 h (2×24 h segments) |
| `p018_L2` | 2.88 | 26 | 0.151 | 0.0193 | 0.052 | ~1.24 | 0.0027 | ~415k | ~90 h cold / ~50 h warm |

Re-budget from Phase 1's measured cost before submitting L2. Driver
`max_particles=500_000`: L2 sits at ~83% of the cap — if the count trends
above ~480k mid-run, stop and raise the constant in
`examples/rotor_hover_pressure_comparison.jl:343` (one-line change; rerun the
pre-submission gate + redeploy).

## Warm-start pilot (gates L2's cost)

Run `p018_L1` twice: cold (restart-chained from its own start) AND
`p018_L1_warm` warm-started from B0's banked 20-rev field
(`RESTART_STEP=<S_B0>,RESTART_NAME=p018_b0,RESTART_PATH=data/p018_b0` + 8
flush revs before the averaging window). Cycle-means agree ≤ 0.25% ⇒ cross-σ
warm starts validated, L2 runs warm from L1. Else L2 runs cold and all later
final-settings runs chain cold (+~40 h budget).

## Decision

Fit CT̄(σ) = CT∞ + A·σ^p over {L0, L1, L2}. Accept if p > 0 and
|CT̄(L1) − CT∞| ≤ 0.75% ⇒ **production σ\* = L1** (0.297c), extrapolation gap
→ error-budget term 4; ε_Γ(L1, L2) ≤ 1%. Harvested 12943696 (σ≈0.37c, legacy
policy) plots as an off-ladder consistency point.

**Contingency chain (in order; record each trigger in the ledger):**
(a) raise OVERLAP one notch across the ladder (e.g. 2.4/2.88/3.46) — the
overlap floor is real (OVERLAP=2 diverged under legacy shedding);
(b) `SHEDDING_R_OVER_R` 0.1 → 0.2/0.3 root-clip test at L1 (cheap hub-cutoff
proxy — near-hub particles barely move relative to the body; 40_40 hub
exclusion measured ≤ +2e-4 CT, so this should be a null unless root particles
drive the instability);
(c) hub-radius mesh variant: add `HUB_R_OVER_R` to
`scripts/generate_dji9443_mesh.sh` (`SetParmVal(prop_id,"RadiusFrac","XSec_0",...)`,
fold into output suffix; Mac-only OpenVSP; keep the 2d cap recipe);
(d) **016 surface-vorticity shedding — PRE-AUTHORIZED** (item ruling 7):
implement per `BRAINSTORM/016_surface_vorticity_particle_shedding/phase_03_implementation.md`
+ its phase_04 V&V, then rerun the ladder with the smooth conversion.

## Exit criteria

σ\* selected with fit + extrapolation term, or contingency outcome recorded;
Γ overlays across rungs in the ledger.

## Log

- 2026-08-01 — **`p018_L1` (13011374) COMPLETED (16.2 h) — σ/c=0.297 rung, but
  NOT M1-settled within its 20 revs.** Knobs metadata-verified (OVERLAP 2.4,
  PPS 11, N=4, η=1.0). History: settled ≈0.0711 at rev 9, then a **large
  excursion over revs 11–14 (peak 0.0802, +12% above the eventual tail)**,
  decaying to a very flat tail — revs 16–19 give **0.07251 ± 0.00007
  (ptp 0.00016)**. The 10–19 window is contaminated (mean 0.07497 ± 0.00311,
  drift −6.58%) and is REJECTED; the 4-rev tail is quoted as provisional only.
  **Finding: settling time grows markedly with σ refinement** — L0 was settled
  by rev ~10, L1 needs ~16. Budget L2 accordingly (its cold-start settling may
  exceed the segment length; the warm pilot matters more than planned).
  Provisional Δ vs B0 (0.07170) = **+1.13%**, i.e. the σ axis moves CT UP —
  opposite in sign to dt refinement. Extension **`p018_L1_s2` (13017106,
  RESTART_STEP=719, SETTLE=24 ⇒ +12 revs, knob overrides re-passed)**
  submitted for a ≥15-rev window before any σ fit is attempted. Do NOT fit
  CT̄(σ)=CT∞+Aσ^p until L1 has a settled mean — a 4-rev tail cannot anchor an
  extrapolation. Warm-pilot comparison (13011375) also waits on this.

  **Ryan hypothesis (2026-08-01): raising OVERLAP (higher q) may itself
  REDUCE the required settle time.** Rationale: at fixed σ, larger OVERLAP
  means smaller inter-particle spacing h, so the discrete particle set
  represents the vorticity field more smoothly — less spurious small-scale
  content to be shed, convected, and dissipated before the wake column
  reaches its attractor. If true, the L1 excursion (revs 11–14) is partly an
  under-resolved-*spacing* transient rather than an intrinsic property of
  finer σ, and the "settling time grows with σ refinement" finding above is
  really "settling time grows when h is under-resolved relative to σ".

  **Consequences if true:**
  (a) L2 should use a *more* aggressive OVERLAP than the planned 2.88 rather
  than economizing, because the extra per-step cost may be repaid in fewer
  settle revs;
  (b) the ladder's q schedule becomes a lever on total campaign cost, not
  just on formal convergence order;
  (c) B0's fast settling (rev ~10 at OVERLAP=2.0) would need explaining — at
  coarse σ the resolvable scales are few, so the transient may be short for a
  different reason; the hypothesis predicts the settle-time *peak* sits at
  intermediate σ with lagging h;
  (d) **σ and Das do not converge independently — raising OVERLAP can
  COMPETE with the Das/clearance convergence** (Ryan). Because
  σ = 2πR·OVERLAP/(NT·PPS), buying overlap by raising OVERLAP *alone* raises
  σ, which shrinks the clearance ratio Das/σ (0.60 at L0, 1.38 at L1) and
  pushes the configuration back toward the regime where Phase 5 measured a
  real penalty (N=1 at d=0.6σ failed by −0.75%/ε_Γ 2.77%). So the overlap
  test must hold σ fixed by raising PPS in step (as the staged case below
  does); otherwise a settle-time win is bought with a clearance loss, and the
  Das\*/N\* selections would have to be re-derived at the new σ. This also
  cuts the other way: the σ ladder is *itself* a clearance sweep at fixed
  physical Das, which is why Phase 4 owns the deferred Das top-two and N=1
  re-tests — the three axes share one dimensionless group.

  **Cheap test (staged, not yet submitted):** rerun L1's σ at a raised
  OVERLAP — `p018_L1_ov3` with `OVERLAP=3.0, P_PER_STEP=14, MERGE_R_FACTOR=0.0053`
  (σ/c ≈ 0.297 held **deliberately fixed** per (d); h/σ 0.42→0.33) — and
  compare the rev at which the history first enters its final ±0.3% band
  against L1's ~16. A shorter settle at equal σ confirms; equal/longer
  refutes. One 24 h job; run it when a slot frees, ideally before committing
  L2's segments.

- 2026-08-01 — **`p018_L1_warm` (13011375) COMPLETED 19:34:17.** Status
  bookkeeping only: the run finished on the cluster and **its output has not
  been retrieved, read, or interpreted** (recorded per Ryan's instruction to log
  completion without touching results). The warm-pilot gate it exists to answer
  — ≤0.25% against the *settled* cold L1 ⇒ L2 may run warm — remains unevaluated
  and is still blocked on its comparator: the cold L1 extension `p018_L1_s2`
  (13017106) was still running at the time of this entry, and without a settled
  cold L1 mean there is nothing to gate against. Do not treat this completion as
  progress on the σ fit.

## 2026-08-03 — L1 stitched and settled; warm-pilot PASSES; M2 fails the σ axis

Analysis via `scripts/p018_analyze.py` (new; validated against the banked
`p018_b0` = 0.07170 over revs 15–29 before any number below was recorded).

### L1 settled mean (supersedes the provisional 4-rev tail)

`p018_L1` (13011374) + `p018_L1_s2` (13017106), stitched on `step` at 719:

| window | CT̄ | per-rev std | 5-rev block drift |
| --- | --- | --- | --- |
| revs 17–31 (15 revs) | **0.072669** | ±0.00107 | 1.17% (non-monotone) |
| revs 16–31 (16 revs) | 0.072664 | ±0.00103 | 0.88% (non-monotone) |
| revs 15–31 (17 revs) | 0.072656 | ±0.00100 | 1.41% (non-monotone) |

The **0.02% spread across three windows** is the settledness evidence, not the
block drift: the drift is limit-cycle amplitude (non-monotone, period ~7–9
revs, consistent with 006), which is precisely the signal
`decision_rules.md`'s limit-cycle defense says must be judged on cycle-means
rather than block tolerances. Adopt **CT̄(L1) = 0.07267, revs 17–31**,
bootstrap CI [0.07209, 0.07347].

**Δ vs B0 (0.07170) = +1.34%** (previous provisional: +1.13%). σ moves CT
**up**, opposite in sign to dt refinement (Phase 3) — the two axes partly
cancel, which is worth stating explicitly in the error budget rather than
letting the cancellation flatter the final number.

### M2 on the σ axis — FAILS, and it is the binding result

ε_Γ (B0 vs L1, matched to each run's settled window, 29 sections on
0.3 ≤ r/R ≤ 0.95): **max 8.776%, RMS 4.692% — FAIL** (threshold 1%).

This is the decision-relevant finding of the session. On CT̄ alone the σ axis
would read as nearly converged (+1.34%, close to the 0.5% gate); on the
co-equal circulation observable it is off by ~9%. Ruling 9 exists for exactly
this case. **Consequences:** (a) L2 is mandatory, not optional; (b) the σ
error-budget term must be set by M2, not CT̄; (c) any temptation to stop the σ
ladder at L1 on cost grounds is refuted.

### Warm-pilot gate — PASSES

`p018_L1_warm` (13011375, 19.6 h; warm-started from B0 at step 719, i.e.
across a σ/c change 0.68 → 0.297). Gate is on **matched** windows (revs 22–31,
per the reporting rules):

| arm | CT̄ (revs 22–31) | per-rev std | block drift |
| --- | --- | --- | --- |
| cold `p018_L1(+s2)` | 0.072411 | ±0.00101 | 0.210% |
| warm `p018_L1_warm` | 0.072280 | ±0.00023 | 0.017% |

**Δ = −0.18%, inside the ≤0.25% threshold ⇒ PASS.** ε_Γ = 0.289% max /
0.148% RMS ⇒ PASS. Cross-σ warm starts are validated; **L2 runs warm.**

Secondary observation worth keeping: the warm arm is ~4× quieter than the cold
arm and shows none of the cold run's revs 11–14 excursion (peak 0.0802). Warm
starting does not merely save wall-clock, it suppresses the settling
excursion — relevant to the `p018_L1_ov3` settle-time question below, and a
reason to prefer warm chaining for the remaining ladder.

(Reporting caveat: on each arm's *own* settled window the values are 0.072669
cold (17–31) and 0.072287 warm (22–37), a −0.53% spread. The gate is defined
on matched windows and the unmatched spread is inside the cold arm's own CI
[0.07209, 0.07347], but it is recorded here so the choice of window is not
silently doing work.)

### Submitted 2026-08-03

| job | case | knobs (verified from the run's log echo) |
| --- | --- | --- |
| 13029923 | `p018_L1_ov3` | OVERLAP 3.0 / PPS 14 ⇒ **σ/c = 0.292** (vs L1's 0.297, −1.8% — σ held only approximately, recorded not glossed), MERGE_R_FACTOR 0.0052 = 0.138·(σ/R), NT=36, η=1.0, N=4, 30 h |
| 13029924 | `p018_L2` | OVERLAP 2.88 / PPS 26 ⇒ σ/c = 0.151, MERGE_R_FACTOR 0.0027, **warm from `p018_L1_s2` step 1151**, `P018_SETTLE_REVS=42` (⇒ ~18 new revs), 48 h |

**A launcher hazard was found and fixed en route.** `p018_L1_ov3` needed its
own case arm: the `p018_L1|p018_L1_warm` arm assigns
`export OVERLAP=2.4; export P_PER_STEP=11` **unconditionally**, so submitting
ov3 as an `--export` override of the `p018_L1` tag would have silently run at
OVERLAP=2.4 — the mirror image of the Phase 3 η leak, in the opposite
direction. The two arms had opposite override precedence and neither warned.
The new `p018_L1_ov3` arm is deployed and md5-verified on the cluster
(`b2a8a2d48028c521662dd88253ef7f2f`).

`p018_L1_ov3` metric (unchanged): the rev at which the history first enters
its final ±0.3% band, vs L1's ~16.

### Still open in Phase 4

- **The σ fit stays blocked until L2 lands** — two points are not a fit of
  CT̄(σ) = CT∞ + Aσ^p. Note also that the fit must be checked against M2, which
  currently fails on the L0→L1 step.
- Deferred re-tests, unchanged: (a) Das top-two rungs at σ=L1 (Phase 2's
  no-plateau branch); (b) N=1 at L1, where d/σ = 1.38 clears the bracketed
  threshold — now additionally framed by Ryan's revised N preference
  (2026-08-03): prefer N=2 over N=1 **if** N=2 permits a significantly smaller
  particles-per-step, since N buys clearance that would otherwise have to be
  bought with a finer (more expensive) σ.

## 2026-08-03 (b) — WHERE the σ-axis ε_Γ lives: it is a redistribution, not a bias

Figures: `data/p018_gamma_ladders/gamma_sigma.png` and
`gamma_sigma_matched.png`, from `scripts/p018_plot_gamma.py` (reproduces
`p018_analyze.py m2`'s 8.776% exactly, so figure and metric agree).

### The shape

The 8.78% is **not** noise and **not** a uniform offset. It is a smooth,
three-lobed redistribution (L1 minus B0, normalized by max_r|Γ̄_B0|):

| region | sign | magnitude |
| --- | --- | --- |
| inboard, 0.30 ≤ r/R ≲ 0.55 | **L1 higher** | up to **+6%** |
| mid-outboard, 0.60 ≲ r/R ≲ 0.85 | **L1 lower** | dip to **−8.8%**, centred r/R ≈ 0.76 |
| tip, 0.88 ≲ r/R ≤ 0.95 | **L1 higher** | +5 to +6% |

### The number that explains the CT̄/Γ̄ paradox

∫Γ̄ dr over the metric band (0.3–0.95): B0 **0.163875**, L1 **0.164067** —
**+0.117%.** Over the full span: +3.21% (that net sits mostly *inboard of the
band*, r/R < 0.3, where L1 is much higher).

So within the band the refinement moves ~9% of peak Γ̄ around **at essentially
constant total**. That is precisely why CT̄ moved only +1.34% while ε_Γ hit
8.78%: **CT̄ is structurally near-blind to this error mode.** Two consequences,
both binding:

1. **The σ error-budget term cannot be taken from CT̄** — CT̄ would report
   ~1%, the distribution says ~9%. Take it from M2.
2. **A σ ladder judged on CT̄ alone would have declared convergence
   prematurely.** Ruling 9 earned its keep here.

### Robustness

Re-running on the *matched* window (revs 22–31, both arms) gives ε_max 8.847%
and a band integral change of −0.094% — same magnitude, same shape, same
conclusion. The finding is not an artifact of the two arms using different
settled windows.

### Mechanism reading (hypothesis, testable at L2)

B0 carries a **flat-topped Γ̄ plateau out to r/R ≈ 0.85**; L1 peaks near
r/R ≈ 0.6 and declines outboard. A flat plateau that far outboard is not what a
resolved tip vortex produces — it is what happens when the tip vortex is
smeared: too little outboard downwash, so the outboard blade stays
over-loaded. Refining σ 0.68c → 0.297c sharpens the tip vortex, increases
outboard induced velocity, and unloads r/R 0.6–0.85 — exactly the observed dip.

If that reading is right, **L2 (σ/c 0.151) should continue to deepen the same
dip rather than introduce a new lobe**, and the band integral should stay near
constant. That is a sharp, falsifiable prediction to check when 13029924 lands,
and it routes the σ axis to BRAINSTORM/012 (σ management / tip-vortex
preservation) rather than to a generic "more resolution" argument.

## 2026-08-03 (c) — `p018_L2` OOM-KILLED: genuine memory, NOT divergence

Job 13029924 (warm from `p018_L1_s2` step 1151, σ/c 0.151, OVERLAP 2.88 / PPS 26)
was **OUT_OF_MEMORY at 03:20:46**. It is the σ ladder's third rung, so this
blocks the CT̄(σ) fit.

**Diagnosis: the rare genuine-memory branch, not a blow-up.** The
BoundCirculation/Force monitor CSV — which is written incrementally and
therefore *survives the SIGKILL* — records 99 clean steps:

| step | 1152 | 1248 | 1249 | 1250 |
| --- | --- | --- | --- | --- |
| CFx | −0.07328 | −0.07602 | −0.07600 | −0.06781 |

Thrust ≈ +0.068–0.076 throughout, bounded and smooth, with **no growth**. The
standing OOM-means-blow-up rule does not apply here.

**Two refinements to the diagnostic rule, both earned the hard way:**

1. **Use the force-monitor CSV, not the stdout log, to test |CT|.** The log had
   **zero** `CF = ` lines because block-buffered stdout was truncated by the
   kill — exactly the failure mode `ops_reference.md` warns about ("absence of
   insane CT in the log is NOT evidence of boundedness"). The monitor CSV closes
   that gap and gives positive evidence either way.
2. **`MaxRSS ≈ ReqMem` is NOT a reliable memory-ceiling test.** Here
   `MaxRSS = 22.3 GB` against a **64 G** request — only 35% — yet the job was
   OOM-killed. `sacct` samples periodically and missed the peak, so a fast
   allocation spike (merge hash / FMM tree rebuild at ~350k particles) is
   invisible to it. A low MaxRSS therefore cannot be used to argue "not a memory
   problem"; it was only the CT history that settled this.

**Why L2 is memory-hungry:** the cost model gives N_p ≈ 52.4k/(σ/c) ⇒ **≈ 347k
particles** at σ/c = 0.151, against the driver's hard-coded
`max_particles = 500_000`. L2 is the first campaign case to approach that.

**Decision needed from Ryan — it collides with a standing ruling.** Ruling 5
fixes memory at **flat 64G**. L2 appears not to fit. Options, cheapest first:

- **Raise this one job to 128G** (breaks ruling 5 for a single case; simplest and
  most likely to work — note `p2e_sigF_nofilt` was the earlier genuine-memory
  case at 32G, and the 64G ruling was itself the fix for that one).
- Reduce the retained particle count at fixed σ: shorten `TRUNCATION_DEPTH_R`
  from 4R for L2 only (changes a converged axis — would need its own null).
- Lower `max_particles` and accept truncation (silently changes physics — not
  acceptable without a null).

**L2 is still the binding item for Phase 4**, and the σ error-budget term must
come from M2 (which failed at 8.78% on L0→L1), so the ladder cannot simply stop
at L1. The 99 recorded steps are usable as a restart source if L2 is relaunched
with more memory (last VTU step is restartable per the usual S = last-VTU rule).

## 2026-08-03 (session c): L2 is NOT a memory problem — deterministic σ-instability, root-caused

The memory framing above is SUPERSEDED. The 128G relaunch (`p018_L2_s2`,
13035813) reproduced step 1250 bit-identically and died at the same step
~1251 with MaxRSS 101.6 GB: the blow-up is deterministic. Root cause (full
forensics + fix inventory in `ledger.md` "σ-axis instability root cause",
2026-08-03 session c):

- rVPM core contraction (σ̇ = −Zσ) under sustained strain in the aged wake
  column thins cores exponentially with NO floor (σ_min 3.7e-3 → 2.1e-5 m
  over ~1000 steps), concentrated in the column-edge shear layer
  (r/R 0.9–1.4, x/R 1.5–3.3);
- the counterbalancing CoreSpreading has been silently inert campaign-wide
  (integration-scheme mismatch, now fixed in src + driver, test-gated);
- finer ambient σ arms the feedback (99th-pct shrink rate 6× L1's), so L2
  ignites where B0/L1 merely accumulate thin cores.

Consequences for this phase: **the σ ladder cannot proceed inviscid** — an
inviscid rVPM hover wake has no steady-state σ distribution, so "converge
σ→0 at fixed Das" was refining toward a limit the model cannot hold. The
implied viscous equilibrium σ_eq = √(ν/Z̄) ≈ 2.1e-3 m sits at L2's shed σ,
i.e. L2 is exactly where viscosity becomes load-bearing. Path (Ryan-ruled
this session): disclose the erratum + measure a B0 viscous A/B
(`p018_b0_visc`), and rescue L2 with `CORE_SPREADING_ACTIVE=true`
(`p018_L2_visc`, warm from step 1249, 128G) once validation passes and slots
free. If viscous L2 still fails (blow-up or M2), 016 remains pre-authorized.

Shedding-method note (Ryan's question): OverlapPPS/SigmaPPS do not address
this failure — contraction is dynamic, ~2R downstream, independent of birth-σ
policy. OverlapPPS's σ ∝ r property is separately relevant to the inboard
d/σ clearance lobe and is logged for Phase 5/012.

## 2026-08-04 — L2 viscous rescue: first attempt lost to the retention sweep; rebuilt and relaunched

**`p018_L2visc` (13036477) FAILED in 3m41s** — `Wake VTS not found:
data/p018_L1_s2/p018_L1_s2_wake1.1.1151.vts` (stderr log; the .out banner was
otherwise correct, incl. `visc:true`). Root cause is operational, not
physics: the 2026-08-04 VTK retention sweep (ledger) deleted the `.vtu`/`.vtp`
restart pieces of every banked run, leaving only `.vtm` index stubs —
**no complete warm-start source survived on the cluster.** Restart needs, per
`src/FLOWPanel_warmstart.jl`: `<run>_body1.pvd`, body `.vtu`, wake
`.{1,2}.<step>.vts`, particles `.vtp` (+ frames from the metadata TOML);
filaments are not read.

**Rescue (Ryan-approved: restore from the Mac archive):** the panel-wake
`.vts` files had survived cluster-side (the sweep's pattern spares
`.vts`/`.vtm`), and `~/p018_L1_ov3_paraview/` holds the complementary pieces.
Uploaded and md5-verified the step-719 set (body `.vtu`, particles `.vtp`,
body `.pvd`) into `data/p018_L1_ov3/` ⇒ `p018_L1_ov3` (σ/c 0.292, 20 revs,
settled tail ≈ 0.0725) is again a complete warm source at step 719.

**Submitted 2026-08-04 (banner-verified):**

| job | case | knobs (from the run's own banner) |
| --- | --- | --- |
| 13047290 | `p018_b0_visc` | B0 carrier exactly (NT=36, η=1.0, OV 2.0/PPS 4, merge_r 0.0120, rlxf 0.3, N=4) + **visc:true**, cold, 24 h — the CoreSpreading erratum A/B (budget row Γ6) |
| 13047296 | `p018_L2_visc` | OV 2.88/PPS 26 ⇒ σ/c 0.151, merge_r 0.0027, η=1.0, N=4, **visc:true**, warm from `p018_L1_ov3` step 719, 128G, SETTLE=50, 48 h — expect an `_s2` chain for the M1 window |

Note the σ jump for the warm start is 0.292 → 0.151 (ratio 1.9, vs the
validated pilot's 0.68 → 0.297 ratio 2.3) and the source is the ov3 arm
(OVERLAP 3.0) rather than L1 — both differences are inside what the
warm-pilot gate covered in spirit but should be disclosed with the result.
L2_visc's decision criteria unchanged: survival past step ~1251-equivalent
wake age, then M2 vs L1 (falsifiable prediction: the r/R≈0.76 dip deepens
along the same 3-lobe shape).

## 2026-08-05 — Phase 12 cross-reference (h-conditionality of this ladder)

This ladder co-refines h faster than σ by design (OVERLAP 2.0 → 2.4 → 2.88,
h ~ σ^q, q ≈ 1.22), so its rung deltas are joint (σ, h) statements. Phase 12
(`phase_12_spatial_rigor.md`) supplies the missing h-invariance leg — an
h-only ladder at fixed B0 σ (h/σ 0.5 → 0.25 → 0.125, submitted 2026-08-05).
The campaign's σ claim is the conjunction (this ladder) ∧ (Phase 12A pass);
Phase 12A carries a pre-registered discriminator for the L1 3-lobe Γ̄ mode
(if it appears at fixed σ, the σ attribution re-attributes to h/overlap).

## 2026-08-05 — Ladder B: σ at FIXED h (Ryan-requested alternative ladder)

The main ladder (now "Ladder A") co-refines h; Ladder B varies σ ALONE at
fixed particle spacing h, completing the one-factor decomposition around B0
(Phase 12A is the third arm: h alone at fixed σ). Mechanics: with
`P_PER_STEP` fixed the `SigmaOverlap` shed count is invariant, so `OVERLAP`
alone scales σ = OVERLAP·h at exactly fixed h = 2πR/(NT·PPS).

**Range limitation, pre-registered:** at fixed h the admissible σ is bounded
below by the overlap condition σ ≳ h (OVERLAP ≥ 1) — below it the particle
field no longer represents a smooth vorticity field. So Ladder B is NOT a
σ→0 convergence ladder; it is the complementary one-factor arm whose rung
deltas, combined with Phase 12A's, test whether Ladder A's joint deltas
separate as Δ(σ,h) ≈ Δσ|h + Δh|σ.

| tag | OVERLAP | PPS | σ/c | h/σ | merge_r | mem/time | status |
|---|---|---|---|---|---|---|---|
| `p018_b0` | 2.0 | 4 | 0.68 | 0.5 | 0.0120 | done | reference 0.07170 |
| `p018_ov1p4` | 1.4 | 4 | 0.477 | 0.71 | 0.0084 | 64G / 48 h | submitted 2026-08-05 |
| `p018_ov1p0` | 1.0 | 4 | 0.34 | 1.0 | 0.0060 | 64G / 48 h | submitted 2026-08-05 |

Cold starts, `P018_SETTLE_REVS=22` (~30 revs) for ≥15-rev M1 windows; cost
~B0 (particle count unchanged — this ladder is cheap). Scoring: M1+M2
successive pairs vs B0, `scripts/p018_analyze.py`, matched windows.

Pre-registered expectations/caveats:

- d/σ is NOT held: σ halving at fixed Das doubles clearance (0.60→1.2 in
  N=1 units at ov1p0) — the same direction as Ladder A, noted not confounded
  away (clearance penalty was only at d/σ ≲ 1, and it improves here).
- σ/chord falls toward the wing-like regime (0.34c at ov1p0) where 014/017
  found the sheet/particle split benign — a Γ̄ move here is informative.
- `p018_ov1p0` (h/σ = 1.0, marginal overlap) is the divergence-risk rung; a
  legacy OVERLAP=2 divergence (12894481) was in the filtered sigC config and
  does not directly predict this one. **If it diverges, that IS the finding
  (the overlap condition binds; ladder is range-limited) — do not chase it
  with dampers.** Shed σ/c 0.34 is far above the L2 σ-collapse regime, so
  the inviscid-blow-up mechanism should not bind; runs are inviscid to match
  B0 exactly.
- Discriminator shared with Phase 12A: does the L1 3-lobe Γ̄ mode track σ
  (should appear here, weaker, ∝ σ-step) or h (should not appear here)?

## 2026-08-05 (c) — viscous L1 submitted (unconfounding the viscous σ pair)

Gap noticed while auditing launchable work: the viscous σ ladder had rungs 1
and 3 in flight (`p018_b0_visc`, `p018_L2_visc`) but NO viscous L1 — so the
critical L1→L2 M2 comparison (pre-registered: "the r/R≈0.76 dip deepens")
would have compared an INVISCID L1 against a viscous L2, mixing physics.
`p018_L1_visc` (job 13051801, Ryan-approved beyond the 10-cap): OV 2.4 /
PPS 11 / merge_r 0.0053, `CORE_SPREADING_ACTIVE=true`, `WAKE_CORE_BETA=1e9`,
warm from `p018_L1_ov3`@719 (same source and −1.8% σ caveat as L2_visc),
SETTLE=32, 64G, 48 h. Readouts: (a) vs `p018_b0_visc` = the same-physics
viscous B0→L1 pair; (b) vs inviscid L1 = the viscosity delta at L1 σ (Γ6);
(c) partner for L2_visc's M2.

## 2026-08-05 (b) — viscous erratum A/B is NULL; `p018_L2_visc` diverged

**`p018_b0_visc` (13047290) COMPLETED (12.4 h).** B0 carrier with
`CORE_SPREADING_ACTIVE=true`, cold, 20 revs. Matched 10–19 window vs B0
(0.071866):

| metric | value | gate |
| --- | --- | --- |
| CT̄ | 0.071810, CI [0.071611, 0.072043] | **ΔCT̄ = −0.078% ⇒ M1 PASS** |
| ε_Γ vs `p018_b0` | 0.190% max / 0.112% RMS | **M2 PASS** |

⇒ **The CoreSpreading erratum is null on the B0 carrier**: campaign runs made
with the silently-inert core-spreading path are not biased at B0's σ. Error
budget: row Γ6 = 0.19% (Γ̄), CT term ≤ 0.08%. This is the disclosure the
erratum needs — the fix changes nothing at production σ. (It says nothing about
the fine-σ rungs, where the instability lives; see below.)

**`p018_L2_visc` (13047296) FAILED after 14.1 h — the viscous σ rescue did not
work.** Warm from `p018_L1_ov3`@719, OV 2.88 / PPS 26 (σ/c 0.151), viscous ON,
128 G. Log tail at step ~1040 (322 steps past resume) shows CF ≈ (0.55, −0.16,
0.20) against a baseline ≈ 0.07 — already diverged — with `sigma growth` still
1.24 and **MaxRSS 39.4 G of 128 G**, no Julia stack trace ⇒ divergence, not
memory (per the corrected OOM triage in `ops_reference.md`, low MaxRSS is not
itself evidence, but the diverged forces are). It died *earlier* in absolute
step count than the inviscid L2 (~1251).

**Confound to disclose with this result:** the run warm-started across a σ jump
0.292 → 0.151 from the ov3 arm while the inviscid L2 was cold, so
"viscous vs inviscid" is not cleanly isolated — what is established is that
`CORE_SPREADING_ACTIVE=true` *as configured here* did not extend the run's life.

**Ryan's ruling 2026-08-05: record and defer.** No σ decision is taken on this;
in particular the **016 contingency is NOT triggered** (016's own rotor A/B
currently fails at rotor σ/chord, phase_04 of that item). The σ axis now waits
on the one-factor decomposition already in flight — Phase 12A (h at fixed σ,
13051772/3) and Ladder B (σ at fixed h, 13051774/5) — which may set the σ budget
term without an L2 rung at all. `p018_L1_visc` (13051801, warm start verified)
remains the viscous middle rung for a same-physics L1 comparison.

## 2026-08-05 (d) — Ladder B rung 3 (`p018_ov1p0`) harvested: the 3-lobe Γ̄ mode tracks σ at FIXED h

Job 13051775 COMPLETED (exit 0, full 1080 steps / 30 revs, MaxRSS 31.1 GB —
the pre-registered divergence-risk rung did NOT diverge; marginal overlap
h/σ=1.0 is usable). Metadata-verified: velocity, OVERLAP 1.0 / PPS 4 ⇒
σ/c 0.34 at fixed h/c 0.34, merge_r 0.0060, N=4, η=1.0, settle 22, cold.
Scored with `p018_analyze.py`, windows 15–29 both arms (B0 stitched):

| | CT̄ (15–29) | 95% CI | per-rev std | block drift |
| --- | --- | --- | --- | --- |
| `p018_b0` (σ/c 0.68) | 0.071701 | [0.071651, 0.071752] | ±0.00011 | 0.12% |
| `p018_ov1p0` (σ/c 0.34) | **0.072548** | [0.072299, 0.073035] | ±0.00088 | 0.68% monotone |

- **M1: ΔCT̄ = +1.18% — FAIL** (window caveat: ov1p0's 15–29 still drifts
  0.68% monotone, so the mean is an upper-ish bound; CI half-width ±0.5%).
- **M2: ε_Γ,max = 7.86% / RMS 4.23% — FAIL**, and the signed profile is the
  **same 3-lobe redistribution as Ladder A's L1**: inboard +6.1% (r/R 0.27),
  dip **−7.86% at r/R 0.762**, tip +5.7% (r/R 0.92).

**Pre-registered discriminator ANSWERED (σ side):** the 3-lobe mode appears
at *fixed particle spacing h*, at nearly the magnitude Ladder A produced with
h co-refined (L1: +1.34%, ε 8.78%, dip −8.8% at r/R 0.76; here σ/c 0.34 vs
L1's 0.297). The mode is a **σ phenomenon, not an h/overlap phenomenon** —
consistent with Phase 4's under-resolved-tip-vortex attribution (→ 012) and
with the h-ladder's complementary prediction (Phase 12A rungs, in flight,
should show it ABSENT); that null remains the closing leg of the argument.
Also note per-rev scatter grew 8× over B0 at h/σ=1.0 (±0.00088 vs ±0.00011) —
marginal overlap is noisy even when stable, worth carrying into the σ*
discussion. `p018_ov1p4` (rung 2) lands shortly and fills the middle point
(σ-linearity check of the mode amplitude).

VTK: ov1p0 remains on the protect list (full VTK intact) until ov1p4 is
scored, in case Ryan wants a ParaView look at the marginal-overlap wake.

## 2026-08-05 (e) — Ladder B COMPLETE (`p018_ov1p4` harvested): mode amplitude ∝ σ-step

Job 13051774 COMPLETED (exit 0, 15:35 h, full 1080 steps). Metadata-verified
(velocity, OVERLAP 1.4 / PPS 4 ⇒ σ/c 0.477 at fixed h, merge_r 0.0084,
settle 22, cold). Windows 15–29 throughout; B0 stitched:

| rung | σ/c | CT̄ (15–29) | Δ vs prev rung | ε_Γ,max vs B0 | dip |
| --- | --- | --- | --- | --- | --- |
| `p018_b0` | 0.68 | 0.071701 | — | — | — |
| `p018_ov1p4` | 0.477 | **0.072257** CI [0.072221, 0.072290], std ±0.000065 | **+0.78% FAIL M1** | **4.12% FAIL M2** | −4.12% @ r/R 0.740 |
| `p018_ov1p0` | 0.34 | 0.072548 | **+0.40% PASS M1** (pair ε_Γ 3.90% FAIL M2) | 7.86% | −7.86% @ r/R 0.762 |

**Ladder B verdict:** CT̄ per-rung deltas shrink (+0.78% → +0.40%, first-order-
like in σ), but **M2 fails every pair** — the 3-lobe Γ̄ mode appears at fixed h
with amplitude ∝ σ-step (inboard +3.8/+6.1, dip −4.1/−7.9, tip +3.4/+5.7 at
σ/c 0.477/0.34), dip stationary at r/R ≈ 0.74–0.76. Combined with Ladder A's
L1 (σ/c 0.297, joint h-refine, ε 8.78%, dip −8.8% @ 0.76): **σ alone
reproduces the Ladder A signature at matched σ — the h co-refinement
contributed little**, subject to the Phase 12A h-ladder null (in flight)
closing the argument from the other side. The σ budget term remains M2-set;
per ruling 9 no σ\* can be declared from CT̄'s shrinking deltas while Γ̄ still
moves 4%/pair. Also noteworthy: per-rev scatter is non-monotone in OVERLAP
(B0 ±0.00011, ov1p4 ±0.000065, ov1p0 ±0.00088) — only the marginal-overlap
rung is noisy.

**VTK deletion log 2026-08-12 ~23:40 MDT (200G budget sweep, protect-list
entries removed by Ryan):** Ladder-B runs swept of VTK — p018_ov1p4
15771MB, p018_ov1p0 15678MB freed; newest-10 restartable steps
(1070–1079) + all CSVs/TOML/monitors retained for both.
