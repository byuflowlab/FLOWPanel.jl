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
