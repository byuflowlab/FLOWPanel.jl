# Phase 17 — N ∝ NT ladder at fixed λ (extent-preserving temporal refinement)

Staged 2026-08-26 (Ryan directive). Status: **LAUNCHED 2026-08-26** — jobs
13490698 (NT36/N2, 48 h) / 13490699 (NT72/N3, 72 h) / 13490700 (NT144/N5,
72 h), cold, P018_SETTLE_REVS=22, 64c/96G, `-p m11-2,m12,m14 --qos=normal`
(non-preemptible). See ledger 2026-08-26 entry.

## Motivation

Wake-row geometry (corrected 2026-08-26; the original staging text wrongly
called PanelWake row 1 "the rigid Das row"): the **rigid Das row is
body-owned** — the `*_tw.vtu` quads, TE extended by Das — and is not one of
the N PanelWake rows. **All N PanelWake rows are free sheet**: node row 1
(upstream edge of the newest row) is pinned on the TE+Das line, node rows
2..N+1 convect freely one step of travel (= |ω×r|·dt) apart, and particles
convert from the *oldest* row, whose upstream edge is the audited clearance
`d_front = |Das| + (N−1)·travel` (driver comment + Phase 12 C1).

The csarc NT ladder (Phase 3, F1b construction) pins Das = λσ (dt-independent)
but holds **N = 1**: one free row spanning [Das, Das+travel], so the converting
row's position is dt-independent (d_front = Das) but its width — and thus where
particle vorticity first appears — shrinks with dt. Meanwhile CT climbs
+2.5–2.7% per NT doubling with no decay (ledger 08-25 FINDING 1), and every
rung is honestly a model-def A/B, not a discretization refinement.

This phase adds the missing rung type: a ladder where **the conversion-front
clearance d_front is dt-independent and only the discretization of the
free-sheet buffer refines with NT**. d_front − Das = (N−1)·travel is held at
one NT36-travel by scaling N−1 ∝ NT; λ is held fixed so Das = λσ is unchanged;
all other csarc pinnings (σ chord-based, NT·PPS = 432, exact-rate rlxf) carry
over unchanged. The only geometric quantity that still varies across rungs is
the O(dt) width of the converting row (total free-sheet extent
N·travel = t₃₆ + travel(NT) → t₃₆ as NT→∞) — a genuinely vanishing
discretization quantity that any dt ladder must shrink, with the deposit
interval converging onto the fixed Das + t₃₆ line from above.

## Construction (all inherited from csarc F1b except N)

| knob | rule | NT36 | NT72 | NT144 |
|---|---|---|---|---|
| N (NWAKEROWS, all free rows) | N−1 ∝ NT, anchor N−1=1 @ NT36 | **2** | **3** | **5** |
| clearance d_front − Das = (N−1)·travel | fixed = travel@NT36 | 1·t₃₆ | 2·(t₃₆/2) | 4·(t₃₆/4) |
| converting-row width = travel | O(dt), shrinks by design | t₃₆ | t₃₆/2 | t₃₆/4 |
| Das | λσ, λ = **2.4** fixed | = | = | = |
| σ | SIGMA_CHORD_FRACTION=0.313 | = | = | = |
| PPS | NT·PPS = 432 | 12 | 6 | 3 |
| rlxf | r(NT)=1−(1−0.3)^(36/NT) | 0.3 | 0.16334 | 0.08539 |
| mesh / FMM / kernel | production 45_185_ct4, Gaussian + tuned FMM | = | = | = |

Why N = {2,3,5} and not {1,2,4} or {2,4,8}: the audited clearance is
`d_front = |Das| + (N−1)·travel` (rotor_hover_pressure_comparison.jl handoff
block) — the upstream edge of the converting (oldest) row. N−1 ∝ NT holds
d_front exactly invariant at every rung. The alternative N ∝ NT ({2,4,8})
would instead pin the downstream edge (N·travel = const) and let d_front
wander by O(dt); the two constructions differ only at O(dt) and converge to
the same limit, but {2,3,5} matches the campaign's Phase-12 clearance metric
and costs fewer rows. {1,2,4} holds neither edge ({0, ½, ¾}·t₃₆).

Why λ = 2.4: the N=1 csarc rungs already exist at λ2.4 for NT {36,72,144}
(p018_csarc_l2p4, _nt72_l2p4, _nt144_l2p4), so each new rung gets a free
N-effect cross (N=1 vs N∝NT at fixed NT) with no extra jobs; λ2.4 is also the
best-behaved λ on the existing ladder.

Verified before launch (2026-08-26): under the DAS_SIGMA_LAMBDA branch,
`nwakerows` does NOT perturb Das — the curvature cap (line 800) requires
DAS_CURVATURE_BETA (unset on csarc), and the other `nwakerows_extent` uses
(lines 820, 828–832) are print/metadata diagnostics (`theta_max_cs` is written
to metadata only, line 1623). **No code change; launcher case blocks only.**

## Arms

| arm | NT | N | cold/chain | walltime |
|---|---|---|---|---|
| `p018_csarc_n2_l2p4` | 36 | 2 | cold, P018_SETTLE_REVS=22 (~30 revs) | 48 h |
| `p018_csarc_n3_nt72_l2p4` | 72 | 3 | cold, P018_SETTLE_REVS=22 | 72 h (may need _s2) |
| `p018_csarc_n5_nt144_l2p4` | 144 | 5 | cold, P018_SETTLE_REVS=22; **will** need _s2/_s3 chains (rev 22 ≈ 4+ days) | 72 h + chains |

Launcher blocks live in `examples/run_dji9443_hover_ct_hpc.slurm.sh` next to
the csarc NT-ladder blocks; identical to the corresponding
`p018_csarc[_ntXX]_l2p4` block except NWAKEROWS. Submission: cold sbatch from
`/home/rander39/projects/FLOWPanel.jl`, `--mem=96G` (SIGBUS at ~67 G RSS on
64 G, per p018_submit_nt_chains_s2.sh), `--export=ALL,P018_SETTLE_REVS=22`.

## Scoring

- Primary: CT̄ over the settled window (chain-method M2 rules per
  decision_rules.md), NT36→72→144 deltas on the N∝NT rungs vs the N=1 rungs.
- If the N∝NT ladder flattens where N=1 does not → the residual NT drift on
  the N=1 ladder is a handoff/free-row discretization artifact, and this
  becomes the production temporal-convergence construction.
- If it still climbs +2.5%/doubling undiminished → remaining suspects are
  integration error proper and the PPS redistribution; the candidate follow-on
  (shed-clock/integration-clock decoupling — shed at fixed Δψ, sub-step the
  dynamics) is deliberately NOT staged yet (Ryan hold, 2026-08-26).
- Secondary: N-effect cross at each NT (new rung vs existing N=1 rung) —
  quantifies the free-row-extent sensitivity as a function of NT.

## Caveats

- Rungs differ from the archived N=1 rungs in row count, so N-cross deltas are
  model-def A/Bs (disclosed); the *ladder itself* (across NT at N∝NT) is the
  rigorous object.
- PPS redistribution (12/6/3) still rides the NT axis — unchanged from csarc,
  deliberately, so the only new axis vs the N=1 ladder is the row construction.
  A broken-rule PPS arm remains a possible Phase-17b addendum.
- NWAKEROWS≥2 adds free wake-panel rows to the solve; per-step cost delta
  expected small vs the particle field but unmeasured — check s/step in the
  first log against the N=1 rung before sizing chains.

## 3R ladder outcome (2026-08-29) — λ converges, NT does not: brainstorm

Result of the 9-arm 3R ladder (SIGMA_CEIL=0.030, TRUNCATION_RADIUS_R=3.0,
fresh from step 0): at fixed NT, CT̄ is nearly flat in λ (NT36:
0.0707→0.0708→0.0710; NT72: 0.0711→0.0718→0.0725; NT144:
0.0732→0.0738→0.0735), but across NT the climb persists
(~+1.5–2%/doubling at λ3.0: 0.07078 → 0.07184 → 0.07383) and the NT144 arms
are the only ones failing the Phase-2e scatter tolerance.

Measured facts undermining "same wake per physical time" (from the
wake-health monitors, λ3.0 arms):

| arm | dnp/step (avg) | dnp/rev | min σ over run |
|---|---|---|---|
| NT36 | 212.8 | 7,660 | 7.13e-4 |
| NT72 | 150.9 | 10,862 | 4.15e-4 |
| NT144 | 100.4 | 14,453 | 4.07e-4 |

So NO — the wake discretization is not NT-invariant: particles per physical
time grows roughly $\propto \sqrt{NT}$ (×1.42, ×1.33 per doubling), and the
minimum core size drops ~40% from NT36 to NT72. The NT axis is therefore
simultaneously a *wake spatial-resolution* axis, not a pure time-step axis.

Candidate mechanisms (brainstorm, unranked):

1. **Shedding geometry tied to the step clock.** Per-step blade travel
   $\Delta\psi = 2\pi/NT$ sets the streamwise spacing of shed particles;
   with sigma_overlap conversion, σ follows spacing, so halving Δt halves
   spacing and shrinks cores → finer, less-diffused wake → higher induced
   loading. CT̄ then converges only as the *wake resolution* converges,
   much slower than the integrator order suggests.
2. **PPS redistribution (12/6/3) riding the NT axis** (already flagged in
   Caveats). The measured dnp/rev is *not* constant, so whatever invariance
   the 12/6/3 rule was meant to enforce, it does not hold — check whether
   the rule targets sheds/rev, particles/shed-row, or overlap ratio, and
   which one actually changes.
3. **Per-step (not per-time) operators applied more often at high NT:**
   rVPM relaxation, divergence filtering, and the σ compression/stretching
   update all act once per step. E.g. relaxation with fixed per-step blend
   factor applies $NT\times$ more often per rev at NT144 — an effective
   NT-dependent dissipation/regularization.
4. **First-shed distance / panel-to-particle handoff** scales with
   $\Delta t \cdot V$: at high NT the first particle sits closer to the
   trailing edge, changing the bound-to-free circulation handoff and the
   blade's own induced velocity (the N∝NT free-row construction compensates
   row *count*, but not necessarily the near-wake σ/spacing).
5. **σ floor/ceil interactions:** SIGMA_CEIL=0.030 is NT-independent, but
   the σ *distribution* shifts with NT (min σ 7.1e-4 vs 4.1e-4); any
   σ-dependent term (Hessian-self, PSE, FMM adequacy) samples a different
   population per arm.
6. **Physical resolution of instabilities:** finer cores + more particles
   at NT144 resolve vortex pairing / wake instabilities that NT36 smears —
   would explain both the CT̄ climb *and* the rev-to-rev scatter failure at
   NT144 (physics, not error; converges only from above, slowly).
7. **Integration error proper** (the pre-registered suspect): but λ-flatness
   at fixed NT and the non-invariant shedding make discretization-of-the-wake
   the stronger prior; a shed-clock/integration-clock decoupling arm (shed at
   fixed Δψ, sub-step dynamics — Ryan hold 2026-08-26) is the clean test.
8. **CT sampling/averaging:** per-rev blocks contain NT samples; azimuthal
   quadrature of the thrust ripple differs across NT (weak — the ±band is
   tiny compared to the climb).

ParaView endpoint sets (first-10 + last-10 steps, λ3.0 arms) downloaded to
`FLOWPanel.jl/data/p018_l3p0_3r_endpoints/{nt36,nt72,nt144}/`. CAVEAT: the
runs keep only the last ~1,500 VTPs (rolling retention), so true steps 0–9
exist only for NT36; NT72 "first" = steps 644–653, NT144 "first" =
steps 2788–2797 (earliest retained). Regenerating true early steps would
need a short (~10-step) rerun per arm.

### Follow-up (2026-08-29): the dnp/rev difference is NOT shed-side — it is removal-side

Determination of why particles/rev differ across NT (λ3.0 arms):

**Gross shedding per rev is nearly NT-invariant, as designed.** Early-wake
np slopes (steps 10 → 10+NT, before any removal matters): NT36 14,662/rev,
NT72 14,707/rev, NT144 17,211/rev. The csarc σ-overlap construction works:
`_shed_particles!(… ::SigmaOverlap)` (FLOWPanel_wake.jl:734) sheds
`p_per_step = max(1, ceil(overlap·dist/σ))` particles at fixed station σ
(σ_j = 0.313·c_j, NT-independent) along each trailing filament, and since
`dist` ∝ 1/NT while the count ∝ dist, particles-per-rev per station is
NT-independent in the continuum. Shedding is trailing-only
(`method_unsteady = NoShed()`, driver line ~650), so no per-step unsteady
contribution. The NT144 +17% excess is the `ceil`/`max(1,·)` quantization:
at NT144 `dist` is 4× smaller, so more stations sit at small
`overlap·dist/σ` where the round-up is a proportionally larger overshoot
(and the `max(1,…)` floor forces ≥ NT sheds/rev/station wherever
`overlap·dist < σ`).

**The mature-wake population difference (229k/325k/433k) is
age-progressive removal, ~2× stronger at NT36.** Evidence:

- np(t) trajectories are identical through ~rev 10, diverging only as the
  wake ages; NT36 peaks at 244k (rev 14), then *shrinks* to ~217k before
  recovering — its mid-run net slope is −4,654/rev while shedding ~14k/rev.
- Axial density in the final frames (0.5R slabs along x): near the rotor
  (x < 0.5R) the arms agree within ~30% (equal shed influx), but the
  NT144/NT36 count ratio grows monotonically with wake age: 1.3 (0–0.5R),
  2.1 (1–1.5R), 2.6 (2.5–3R). Old wake is being eaten at NT36.
- Survivor σ distributions are nearly identical across arms (median
  0.0041–0.0045 m), so whatever removes particles does not leave a
  grown-σ signature in the survivors.

The only removal processes active in the mid-wake are merging
(MERGE_PARTICLES=true, MERGE_R_FACTOR=0.0055, identical across arms) and
3R/4R truncation (barely engaged — final wake extents are ±~2.8R lateral,
3.5R deep). So the working hypothesis is **NT-dependent merging**: identical
merge settings acting on wakes whose *positional micro-structure* differs —
coarser per-step transport at NT36 (4× larger convection error per step,
plus per-step Pedrizzetti relaxation at rlxf 0.3/0.163/0.085 — the
exact-rate rule equates the linear rate but not the nonlinear per-step
alignment) plausibly produces more close pairs in the aged wake. Not yet
proven: FLOWVPM does not log merge counts. A cheap instrumentation
(println(stderr,…) of merges/step — NOT @info, which is suppressed) on a
short restart would settle it, and would also say whether the CT̄ climb with
NT is really "NT36 merged its old wake away" (i.e. the NT ladder is
secretly a *merge-error* ladder).

### Follow-up 2 (2026-08-29): why NT144's gross shed/rev is +17–21% — the max(1,ceil) floor binds the inner blade

Measured mature per-station shed counts (final frames; particles with
pristine station σ and x < 0.02 m = exactly the newest shed row, counts over
2 blades; 41 stations). Per-filament counts k_j = count/2:

- NT36: k_j runs 2 (r/R 0.15) → 21 (tip); total 562/step → 20,232/rev.
- NT72: k_j ≈ ceil(k_j(NT36)/2) at every station; total 297 → 21,384/rev (+6%).
- NT144: k_j = 1 for ALL 18 stations inboard of r/R ≈ 0.47 (floor-bound),
  reaching 6 at the tip; total 170 → 24,480/rev (+21%).

Mechanism: `p_per_step = max(1, ceil(overlap·dist/σ))` with dist ∝ 1/NT.
The count is per-step-quantized; per rev this is NT·max(1, ceil(x_j/NT·36))
for continuum x_j. When the per-step argument is large (NT36: 2–21) the
ceil overshoot is a few %. At NT144 the inner half of the blade has
per-step arguments ≲ 1, where max(1,·) pins the count at 1/step → 144/rev
per filament regardless of the continuum value (≈ 40–140/rev) — a 1×–3.6×
local inflation. Summed over stations that yields the observed +21%
(mature) / +17% (early) excess. NT72's counts (1–11/step) sit in between
(+6%). So the NT·PPS=432 invariance holds well outboard (where the
circulation lives) and breaks progressively inboard as NT grows; the excess
particles are small-x_j inboard particles carrying little Γ, consistent
with the CT̄ effect of NT NOT being primarily shed-side.

(Empirical trick for future reference: newly shed particles are
identifiable in any VTP frame as those whose σ exactly matches the station
σ set — core spreading/rVPM σ evolution re-stamps every older particle —
AND x within one row of the TE; counts come out as clean even numbers
(2 blades × integer), confirming exactly one pristine row per frame.)
