# Phase 3 — Timestep at Fixed Physical Das AND σ

**Objective:** temporal convergence with both known confounds pinned. The
naive NT=72 run (0.07337, +2.9%) halved Das, σ, and rlxf at once — it is the
*confounded control* (a publication figure on why pinning matters), not a
ladder point.

## Pinning recipe (the crux)

At fixed OVERLAP, σ ∝ OVERLAP/(NT·PPS): **halve `P_PER_STEP` as NT doubles**
(the `SigmaOverlap` policy computes the true per-step count itself). Scale η ∝
NT to fix Das. Halve `RELAX_RLXF` per NT doubling (preserves the physical
relaxation rate — same convention as 2e). Never `DAS_REFRESH`.

## Cases (carrier: Das\*, σ from the Phase-2 carrier σ/c=0.68, N=4)

| tag | NT | η | OVERLAP | PPS | rlxf | time | run if |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `p018_b0` (done) | 36 | η\* | 2.0 | 4 | 0.3 | — | — |
| `p018_nt72` | 72 | 2·η\* | 2.0 | 2 | 0.15 | 48 h (or 2×24 h segments) | always (unless 12950996 ≈ 0.0713 AND Ryan accepts the legacy-σ pin as sufficient — then this is the single confirmation) |
| `p018_nt144` | 144 | 4·η\* | 2.0 | 1 | 0.075 | 2–4 × 48 h segments (restart-chained) | only if |CT̄(72) − CT̄(36)| > 0.5% |

Launcher defaults assume Das\* = 0.41c (η\*=1.0); override
`DAS_ETA_KINEMATIC` at submission if Phase 2 chose differently.
Harvested 12950996 (pinned at legacy σ) corroborates whichever way it landed
(see Phase 0 branch table).

## Decision

NT\* = smallest NT with |CT̄(2NT) − CT̄(NT)| ≤ 0.5% and ε_Γ ≤ 1%. With three
rungs, Richardson-extrapolate and report the observed order. dt residual →
error-budget term 3. Expected NT\* = 36.

## Exit criteria

NT\* selected; confounded-control comparison figure data recorded; ledger
updated.

## Log

- 2026-07-31 — Phase 0 branch resolved: `p2e_nt72_das2p0_ov6` (12950996) gave
  CT̄(10–19) = 0.06852 ± 0.00157 at legacy-σ pinning — **−3.9% below the NT=36
  partner (0.07133)**, below both anchors. So: `p018_nt72` runs unconditionally,
  and NT=144 will almost certainly be needed (|Δ| ≫ 0.5% at legacy σ; verify on
  the campaign pair before committing NT=144 segments). dt refinement direction
  is DOWN — the naive NT=72's +2.9% was the Das confound, opposite in sign to
  the true dt effect. Caveat carried: the 12950996 window had −1.14% 5-rev
  drift (not fully settled), so the −3.9% may be partly settling; the campaign
  `p018_nt72` should run long enough to pass M1 before the Richardson step.
  Submission still gated on Phase 2's Das\* (ladder jobs in flight).

## 2026-08-03 — `p018_nt72` harvested; η leak found; corrected ladder relaunched

**The run did not execute the intended configuration.** `p018_nt72`
(job 13011373, COMPLETED 2026-08-02, 36.9 h) reports
`das_eta_kinematic = 1.0` in its `_case_metadata.toml`, not the intended 2.0.

Root cause: the launcher arm is
`p018_nt72) export NT=72; export DAS_ETA_KINEMATIC="${DAS_ETA_KINEMATIC:-2.0}"; ...`
— a *defaulted* assignment. `sbatch --export=ALL` propagates the submitting
shell's environment, so an `DAS_ETA_KINEMATIC=1.0` exported there (the B0
carrier value, likely left over from a previous submission in the same shell)
silently won over the case default. Nothing warned.

Consequence: Das = 0.41·1.0·(36/72) = **0.205c**, half the intended 0.41c.
σ/c = 1.363·(36/72)·2.0/2 = **0.68**, correctly pinned. So dt and Das were
refined *together* — the exact confound Phase 0 diagnosed in the naive NT=72.

**Salvage — it is a valid rung, just a different one.** At Das = 0.205c and
σ/c = 0.68 it is an exact match for `p018_das0p5` (NT=36, η=0.5, PPS=4,
OVERLAP=2.0; rlxf 0.3 → 0.15 is the intended halving-with-dt). Matched
windows, revs 10–19:

| | NT | Das | σ/c | CT̄ | per-rev std |
| --- | --- | --- | --- | --- | --- |
| `p018_das0p5` | 36 | 0.205c | 0.68 | 0.070056 | ±0.00031 |
| `p018_nt72` | 72 | 0.205c | 0.68 | 0.070473 | ±0.00017 |

**ΔCT̄ = +0.595% — FAILS M1 (>0.5%).** M2: ε_Γ,max = 0.979%, RMS 0.772% —
passes, but only just (threshold 1%). Both windows drift −0.32%/−0.33%,
monotone and in the same direction, so the *delta* is more trustworthy than
either mean; neither window reaches the ≥15-rev M1 length.

**Reading:** dt is not converged at NT=36 even with Das held fixed. This is
weaker than the legacy pinned pair's −3.9% but the same sign of conclusion,
and it is now free of the Das confound. **NT=144 is justified.**

Note the sign structure, which is the whole reason the confound mattered:
refining dt at fixed Das moves CT **up** slightly (+0.60%), whereas the legacy
pinned pair moved it **down** 3.9%. These disagree, and the legacy pair used
`overlap_pps` shedding with nwakerows=1 — the corrected rung below is what
settles it.

### Submitted 2026-08-03

| job | case | knobs (verified from the run's own log echo) |
| --- | --- | --- |
| 13029921 | `p018_nt72_eta2` | NT=72, **η=2.0** ⇒ Das 0.41c, OVERLAP 2.0/PPS 2 ⇒ σ/c 0.68, rlxf 0.15, N=4, 48 h |
| 13029922 | `p018_nt144` | NT=144, **η=4.0** ⇒ Das 0.41c, OVERLAP 2.0/PPS 1 ⇒ σ/c 0.68, rlxf 0.075, N=4, `P018_SETTLE_REVS=9` (⇒ ~17 total revs, ~45 h, inside the 48 h wall), 48 h |

`p018_nt144` is sized to fit one wall-clock segment; extend it with a
restart-chained `_s2` segment for the ≥15-rev M1 window (recipe in
`ops_reference.md` — restart from the last **VTU** step, one less than the
CSV's last `step`).

**Process change adopted (see also `ops_reference.md`):** after every
submission, read back the knobs the job itself reports — either the log
banner (`mesh:… NT:… das_eta:… overlap:… pps:…`, printed within seconds) or
`_case_metadata.toml` — and confirm them against intent *before* the job burns
hours. All four jobs submitted 2026-08-03 were verified this way.

### VTK retention (ruling 10), 2026-08-03

Deleted VTK for `p018_nt72` (23 G → 42 M), keeping its CSV/TOML/monitors: it is
harvested, and superseded as a Phase 3 rung by `p018_nt72_eta2`. Its salvaged
matched-Das value (0.070473) is fully preserved in the CSVs and `ledger.md`.

## 2026-08-03 (b) — the dt error is broad and near-uniform, unlike every other axis

Figure: `data/p018_gamma_ladders/gamma_dt.png`.

For the matched-Das pair (`das0p5` NT=36 vs `nt72` NT=72, both Das 0.205c,
σ/c 0.68): ε_Γ max **0.979%**, RMS **0.772%**, mean signed difference
**+0.711%**, worst point at r/R 0.470 — i.e. mid-span, with max only 1.27× the
RMS.

That ratio is the informative part. Every other axis in this campaign produces
a *localized* error (σ: a 3-lobe redistribution peaking at r/R 0.76; Das and N:
inboard deficits at r/R 0.31). The dt error is instead a **broad, nearly
uniform positive shift of the whole distribution** — no lobe structure, no
localization. Band integral ∫Γ̄ dr rises +0.779%, tracking the pointwise mean
almost exactly.

This is what a genuine global time-integration error should look like, and it
supports reading the +0.595% CT̄ delta as real dt truncation rather than a
geometric or clearance artifact. It also means dt, unlike σ, is an axis where
**CT̄ and Γ̄ agree** (+0.595% vs +0.711%) — so the dt budget term can be taken
from either, and the NT=144 rung will be interpretable on CT̄ alone.

## 2026-08-05 — `p018_nt72_eta2` harvested: the production-Das dt rung confirms NT=36 unconverged

Job 13029921 COMPLETED (~28 h). Metadata-verified: NT=72, η=2.0 ⇒ **Das=0.41c
pinned**, OVERLAP 2.0/PPS 2 ⇒ **σ/c=0.68 pinned**, rlxf 0.15, N=4.

| | CT̄ (revs 10–19) | per-rev std | block drift |
| --- | --- | --- | --- |
| `p018_b0` (NT=36) | 0.071866 (matched window; settled 15–29 = 0.07170) | ±0.00034 | — |
| `p018_nt72_eta2` (NT=72) | **0.072243**, CI [0.072139, 0.072357] | ±0.00019 | 0.275% monotone |

**ΔCT̄ = +0.52% matched / +0.76% vs B0's settled mean — FAILS M1 both ways.**
M2 vs B0: ε_max **0.865% / RMS 0.705% — PASS**, signed shift **uniformly
positive** (mean +0.66%, range +0.08…+0.87%, no lobe structure) — the same
broad global-truncation signature as the matched-Das 0.205c pair (+0.595% CT,
+0.711% mean Γ̄). Two Das values now give the same sign, magnitude class, and
spatial structure for the dt effect ⇒ the axis is well-behaved and **NT=144
(13029922, running) is decision-critical for the Richardson step**. Window
caveat: 10 settled revs (short of M1's ≥15) with a 0.275% monotone drift —
extend by restart if the Richardson order comes out implausible; the legacy
pinned pair's −3.9% (overlap_pps, N=1) is now clearly the outlier and its
sign disagreement is attributed to the legacy shedding config, not to dt.

## 2026-08-05 (b) — `p018_nt144` chained (`_s2`, job 13054739)

The parent (13029922) reached step 2014/2447 with ~1.9 h of walltime left at
~82 s/step and could not finish. Chained per `ops_reference.md` §Restart
chaining from **step 2025** (complete restart set verified: `body1.2025.vtu`,
`wake1.{1,2}.2025.vts`, `wake1_particles.2025.vtp`):

```
sbatch --job-name=fp-018-nt144s2 --time=48:00:00 \
  --export=ALL,P018_RUN_NAME=p018_nt144_s2,P018_SETTLE_REVS=16,DAS_ETA_KINEMATIC=4.0,\
RESTART_STEP=2025,RESTART_NAME=p018_nt144,RESTART_PATH=data/p018_nt144 \
  examples/run_dji9443_hover_ct_hpc.slurm.sh p018_nt144
```

`P018_SETTLE_REVS=16` ⇒ `schedule_revs = 2+3+4+16 = 25` revs total (3600 steps),
so both the matched 10–19 window and a ≥15-rev M1 window (10–24) will exist;
~1575 steps ≈ 36 h remain, inside the 48 h request. `DAS_ETA_KINEMATIC=4.0` is
passed explicitly because this case arm uses `${VAR:-4.0}` (environment wins —
the precedence that cost 36.9 h in Phase 3). **Banner verified:**
`NT:144 das_eta:4.0 overlap:2.0 pps:1 merge_r:0.0120 settle:16`.

**Operational hole found and worth fixing before the next walltime kill:**
`<run>_CT_vs_rev.csv` is written only *after* the time loop completes
(`examples/rotor_hover_pressure_comparison.jl:897–908`), so the parent segment
will die leaving no CT series file — only the incrementally-written
`monitors/<run>_monitor02_force_system1.csv` survives. `scripts/p018_analyze.py`
reads `CT_vs_rev.csv` exclusively, so it needs a `--from-monitor` fallback
(reconstruct `step,revolution,CT_*`; `rev = step·dt·RPM/60`; validate by
reproducing a completed run's file, e.g. `p018_dasc0p41`) before the Richardson
stitch can use the parent segment.

Phase 3 close-out on landing: Richardson over {NT 36: 0.071866, NT 72: 0.072243,
NT 144: TBD} on matched 10–19 windows at Das 0.41c.

**Warm-resume verified (2026-08-05, 12:0x MDT):** the `simulate_warmstart!`
line was still inside the block-buffered stdout, so the check was made on the
incrementally-written monitor instead — `p018_nt144_s2`'s force CSV opens at
**step 2026** (a cold run opens at 0/1) with CT = 0.071885 against the parent's
0.071919 at step 2025, a **seam jump of 0.05%**, at the restart-integrity gate.
Both monitors (force + bound circulation) are being written.

**Re-confirmed 2026-08-05 ~16:20 MDT (fresh session, direct ssh):** job
13054739 RUNNING (4:30 elapsed of 48 h), output directory + both monitors
present, force CSV spans steps **2026 → ~2125** of 3600; ~34 h of work remain
against ~43 h of walltime — on track, no `_s3` expected. The ledger watchdog's
repeated "no output directory after 3h10m" observations (~13:05–15:00 MDT)
were **wrong/stale** — this verification stands. Note also the parent sweep
(ledger ~14:00) retained steps 2047–2056, not 2025; that was harmless only
because `_s2` had already read the step-2025 set before the sweep. Richardson
close-out proceeds on landing, with the phase_13 §5 caveat: at fixed N=4 the
NT ladder degrades inboard handoff clearance (min d/σ 0.713 → 0.414 → 0.264),
so the observed order is reported but NT\* stays provisional pending a
clearance-pinned re-run.

## Close-out (2026-08-12)

NT=144 chain landed (13054739 `_s2` + 13134726 `_s3`; banner das_eta:4.0
verified against parent; wake health clean, min_sr ≥ 0.47). Three-segment
stitch seams 0.0515%/0.0575% on adjacent rows — at the nominal gate, within
natural per-step variability (ledger §2026-08-12). Matched 10–19 window:
**CT̄(144) = 0.072023, CI [0.071952, 0.072085]** (long 10–23: 0.072012).

Richardson over matched 10–19 windows at Das 0.41c, fixed N=4:

| NT | CT̄ | Δ vs previous rung |
| --- | --- | --- |
| 36 | 0.071866 | — |
| 72 | 0.072243 | +0.52% |
| 144 | 0.072023 | −0.31% |

**The sequence is NON-MONOTONE (difference ratio −1.72): no observed order,
no Richardson extrapolation.** The 72→144 step is ~3× the bootstrap CI
half-widths, so it is not pure sampling noise; it is consistent with the
phase_13 §5 confound — at fixed N=4 the ladder degrades inboard handoff
clearance (min d/σ 0.713 → 0.414 → 0.264), so the rungs are not a clean dt
refinement. **Phase 3 verdict: dt residual bounded at ≤0.5% peak-to-peak
over NT 36–144; observed order unobtainable without a clearance-pinned
ladder (N = 4/7/13, or the §4b N=1 construction which pins automatically).
NT\* remains provisional per the standing caveat.** Error-budget dt entry
should carry ±0.5% with this caveat.

## PRE-REGISTRATION (2026-08-12): pinned dt ladder — DRAFT FOR RYAN, no submission before his sign-off

Reopens Phase 3's closed caveat with the clearance-pinned construction.
Decisions already made by Ryan 2026-08-12 (item file "STAGED DECISION" +
ledger): schedule NT 36/72/144 at PPS 12/6/3; carrier = §4b N=1 uniform-D
(the etadas A/B came back ~NULL in CT, so uniform-D survives as the
carrier law and rung 1 pins automatically); cost ACCEPTED.

### Construction (what is pinned, and how)

| quantity | law | rung values (NT 36/72/144) |
| --- | --- | --- |
| NT·PPS | = 432 fixed | 432/432/432 |
| σ | 2πR/NT · OVERLAP/PPS, OVERLAP=2.75 | 0.04R at every rung |
| h (streamwise sheet sampling) | ∝ 1/(NT·PPS) | pinned |
| d_front | Das_j = D·σ − (N−1)·travel_j, D=3.4, N=1 | 3.4σ span-uniform, **dt-independent by construction** |
| merge_r | 0.138·σ/R (σ pinned) | 0.0055 at every rung |
| rlxf | halved per NT doubling (standing Phase-3 rule) | 0.3 / 0.15 / 0.075 |
| viscosity | CoreSpreading ON (production physics) | all rungs |
| mesh | production 45_185_ct4 | all rungs |

Contrast with the closed ladder: there min d/σ degraded 0.713→0.414→0.264
down the rungs (fixed N=4, η∝NT); here the handoff clearance is identical
at every rung, so ΔCT across rungs reads dt alone.

### Rungs

| rung | run | status | knobs (launcher arm) |
| --- | --- | --- | --- |
| NT=36 | `p018_ufront_n1_s040_visc` (+`_s2`, 45 revs) | EXISTS (rung 1 by the etadas-null verdict) | OVERLAP 2.75 / PPS 12 / rlxf 0.3 |
| NT=72 | `p018_upin_nt72` | NEW | NT=72, PPS=6, rlxf 0.15, OVERLAP 2.75, merge_r 0.0055, N=1 D=3.4, visc |
| NT=144 | `p018_upin_nt144` | NEW | NT=144, PPS=3, rlxf 0.075, otherwise ↑ |

Both new rungs carry the Phase-15 diagnostics (WakeInventoryMonitor +
attribution, default ON) — **the NT=72 rung doubles as the H1/H3
measurement run**, per the Ryan-approved sequencing (instrumentation gates
the ladder).

### Windows and scoring

- **Primary: matched revs 15–29** (rung 1 has it; both new rungs run
  SETTLE 22 → ~30 revs cold, same protocol as the σ-ladder A/B family).
  `p018_analyze.py` M1 + M2 per `decision_rules.md`, unchanged.
- **Richardson** on the matched-window CT̄ triple; report observed order;
  update budget rows 3/Γ3. Pre-registered expectation: with clearance
  pinned, successive deltas shrink and the sign is consistent (the closed
  ladder's non-monotonicity was the confound's signature).
- **Burst decomposition REPORTED alongside (report-only, per Ryan):**
  quiet-regression limit, family-pooled quiet-rev mean (s* pooled over the
  three rungs — same carrier family), burst incidence, quiet-series drift
  with AR(1) CI (`scripts/p018_burst_stats.py`). NOT a gate; context for
  the M1 CHECK flags the family carries.
- **Post-ladder discriminator (queued, comes free):** compare the CT
  fluctuation's band, episode timescale (~5–10 revs), and per-rev std
  across rungs — a numerical mode moves with Δt, a physical wake mode does
  not. Fold into the min_sr attribution.

### Per-rung stability screens (before trusting any rung)

σ_stab = √(Γ_v Δt/2π) and the viscous floor √(2νΔt) both SHRINK with Δt,
but the 014 screens showed refinement TIGHTENS the inviscid ignition
boundary — do not assume halving Δt buys margin. Per rung: WakeHealthMonitor
tripwire ON (it is, by default); kill + flag on σ<0, unarrested min_sr
collapse, or max_u runaway; record min_sr trajectory vs the NT=36 arrest
shape (dip ~0.098 @ rev ~31, recovery). If a rung ignites, STOP the ladder
and bring the regime map (019) before resubmitting anything.

### Cost (accepted, restated)

Per-rev particle load is NT-independent (PPS ∝ 1/NT), so wall-clock ≈
∝ NT: NT=36 ran ~30 h/30 revs ⇒ NT=72 ≈ 60 h, NT=144 ≈ 120 h. **Ryan
2026-08-12 (at approval): submissions may use up to 72 h walltime when
needed** ⇒ NT=72 fits ONE job (no chain); NT=144 = one 72 h segment + one
`_s2` chain (`ops_reference.md` recipe). 64G flat; NT=144 to be watched
for particle-count cap (500k) though counts should match the NT=36 rung's
~215k.

**APPROVED by Ryan 2026-08-12 as drafted (with the 72 h allowance).**

### Gates before ANY submission

1. p15 instrumentation gate job PASSES (full suite + formulation_test +
   smoke + bit-identity) — take-8 13157614 in flight at drafting time.
2. Launcher arms `p018_upin_nt72` / `p018_upin_nt144` added, deployed,
   md5-verified (NOT yet written — waiting on this pre-registration's
   approval).
3. Banner verification within minutes of each start: NT, pps, overlap
   2.75, das_uniform:3.4, nwakerows:1, merge_r:0.0055, rlxf, visc:true,
   WAKE_INVENTORY=true.
4. **Ryan approves this pre-registration.**

**VTK deletion log 2026-08-12 ~23:15 MDT (200G budget sweep):** dt-ladder
rungs swept of older-step VTK via `p018_vtk_sweeper.sh` while live —
upin_nt72 837MB (kept 128–137), upin_nt144 975MB (kept 171–180); newest-10
restartable sets retained (nt144's `_s2` chain source is safe by
construction), monitors/CSVs untouched.

**VTK deletion log 2026-08-13 ~19:30 MDT:** live sweep (disk alert) —
upin_nt72 22,519MB freed (kept 856–865), upin_nt144 23,740MB (kept
1146–1155); newest-10 restartable sets + monitors/CSVs retained.

## STAGED (Ryan 2026-08-13): burst-attribution discriminator on the upin ladder

On landing of 13157833 (nt72) / 13157834 (nt144), IN ADDITION to the
pre-registered dt-ladder scoring, run the burst attribution comparison
across NT 36/72/144 (rung 1 = the n1_s040_visc chain, incl. `_s3` revs
45–58):

1. min_sr fluctuation band — carrier reference 0.09–0.17;
2. burst-episode timescale — carrier reference ~5–10 revs (per-rev
   within-rev std series, `p018_burst_stats.py` episodes);
3. CT per-rev std magnitude.

Rule: a NUMERICAL mode moves with Δt; a PHYSICAL wake mode does not.
The verdict decides whether bursts belong in the M1 mean — i.e. whether
the carrier settledness claim rests on the raw mean or on the stationary
quiet-limit ≈0.0730 (`_s3` harvest, phase_15 §2026-08-13) — and GATES A
RYAN DECISION on the M1 observable. Record verdict + numbers in
ledger.md and fold into the min_sr attribution (per-rev min_sr–CT
correlation NULL at 45 revs, detrended r≈+0.08). Note nt144 rungs run
per-step rlxf 0.075 (physically equivalent to stock 0.3 by dt-scaling)
— cite the rlxf-ignition asymmetry (phase_09 §2026-08-13) when
interpreting any nt144 burst differences.

## 2026-08-15 — upin_nt72 IGNITED (stability screen FAIL) ⇒ dt ladder STOPPED

Job 13157833 FAILED (exit 1) at step 1411 = rev 19.6 after 2 d 9.6 h.
Wake-health forensics (local harvest `data/p018_upin_nt72/monitors/`):
precursor at rev 14–15 (min_sr dip 0.0967 → **0.0141** with |Γ|/σ² spike
1.48e5 at rev 15, then a TRANSIENT recovery to 0.0975 at rev 16 — an
arrest attempt), re-deepening 0.086/0.072 revs 17–18, then catastrophic
Γ-ignition rev 18.6–19.6: n_p 242k→50k→28k (self-annihilation), max_u
1.04e6 m/s, min_sigma pinned at 6.65e-5 (the NT=72 viscous floor
√(2νΔt)), |Γ|/σ² 4.8e10, force blow-up (CF ~ −5.3e5) and crash.

**Reading vs the NT=36 arrest shape** (dip ~0.098 @ rev ~31, recovered):
at NT=72 the dip comes ~2× earlier in revs, ~7× deeper, and does NOT
arrest. Since rlxf 0.15 at NT=72 is dt-equivalent to stock 0.3 (the
pre-registered halving rule), this ignition reads on **dt refinement
itself at matched physical relaxation, matched σ=0.04R, matched
clearance** — consistent with the 014/019 expectation that refinement
TIGHTENS the ignition boundary (σ_stab and the viscous floor both shrink
with Δt). Per the pre-registered screen rule: **LADDER STOPPED — no
resubmission before the 019 regime map goes to Ryan.** Consequence
applied: nt144 (13157834, healthy at rev 14.0 at check time, |Γ|/σ²=92)
will be allowed to TIMEOUT at its 72 h wall (~rev 14.6) and harvested,
but its pre-authorized `_s2` chain is NOT submitted (the chain assumed a
passing screen).

**Richardson triple: DEAD at these knobs** — no rung ≥ NT=72 can reach
the matched 15–29 window. Salvage below.

### Burst attribution — PARTIAL (NT 36 vs 72, shared clean window revs 7–13)

| metric (revs 7–13) | NT=36 carrier | NT=72 | read |
| --- | --- | --- | --- |
| CT̄ | 0.07513 | 0.07490 (−0.31%) | dt effect on mean small |
| mean within-rev CT std | 0.00303 | 0.00199 (×0.66) | moves with Δt ⇒ leans NUMERICAL |
| min_sr band | 0.187–0.575 | 0.206–0.607 | similar (early-rev band, not the settled 0.09–0.17 reference) |

Caveats: window is early (bursts of the carrier's settled era start
later); at rev 14 BOTH rungs burst simultaneously with similar magnitude
(0.0077 vs 0.0070) — co-timed across Δt ⇒ that episode looks PHYSICAL.
Mixed verdict pending nt144's revs 7–14 (3-point std-vs-Δt trend) —
complete this on its harvest, then record the verdict in the ledger and
fold into the min_sr attribution. The M1-observable decision stays
RYAN'S, now with the added caveat that the discriminator window is
pre-settled-era only.

## 2026-08-15 (night) — nt144 harvested; attribution verdict: NUMERICAL; two NT=72 model-definition arms launched (Ryan-directed)

**nt144 (13157834) TIMEOUT at 72 h wall, rev 14.45, HEALTHY throughout**
(rev 14: min_sr 0.384, |Γ|/σ² 92, max_u 12 — far inside every envelope).
No `_s2` chain per the screen-fail stop (and per the adopted do-not-run
list). Local harvest `data/p018_upin_nt144/monitors/`.

**Burst attribution — 3-point trend (shared clean window revs 7–13):**

| rung | CT̄ (7–13) | mean within-rev CT std |
| --- | --- | --- |
| NT=36 (rlxf 0.3) | 0.07513 | 0.00303 |
| NT=72 (rlxf 0.15) | 0.07490 | 0.00199 (×0.66) |
| NT=144 (rlxf 0.075) | 0.07775 | 0.00035 (×0.18; ×0.12 vs NT36) |

**Verdict: the CT fluctuation/burst mode is NUMERICAL — it moves
strongly and monotonically with Δt.** The rev-14 burst that hit NT=36
and NT=72 co-timed (0.0077/0.0070) is ABSENT at NT=144 (0.00043), so
even that episode is Δt-sensitive. Caveats: (1) window is the
pre-settled era only (revs 7–13; no rung ≥72 reaches the settled era);
(2) per-step rlxf differs across rungs by the pinned dt-halving design —
the two new arms below will separate rlxf-model from Δt; (3) **CT̄ is
NON-MONOTONE in NT (−0.31% then +3.8%)** on this window — the dt axis
remains unresolved for the mean, only the fluctuation attribution is
answered. Implication staged for Ryan's M1 decision: bursts being
numerical rectification supports the QUIET-LIMIT as the physical
observable (raw-mean drift = rectification artifact), subject to
caveat (1). Fold into the min_sr attribution when scoring the new arms.

**Launched (Ryan-directed 2026-08-15, explicitly MODEL-DEFINITION
CHANGES, not dt-convergence rungs; his direction supersedes the
stop-rule's resubmission hold):** new launcher arms + tags (provenance
rule), deployed md5 4135f810c9efa645f7872dde4ec4a263 local==cluster:
- **13183888 `p018_upin_nt72_rlxf0p3`** — NT=72 pinned construction,
  per-step rlxf 0.3 (stock per-step rate ⇒ ~2× the per-time damping of
  the NT=36 carrier).
- **13183889 `p018_upin_nt72_rlxf0p16334`** — same, rlxf = 1−√(1−0.3)
  = 0.16334 (exact continuous-rate equivalent of stock; tests whether
  the 9% dose over the ignited 0.15 changes the outcome).
Both 72 h/64G, SETTLE 22 (30 revs), cold. Banner-verify on start
(NT:72 pps:6 rlxf:0.3|0.16334, overlap 2.75, das_uniform:3.4,
merge_r 0.0055, visc:true). Score: stability first (ignition times vs
nt72's rev 14.5 precursor / 18.6 catastrophe), then matched windows vs
the NT=36 carrier if either survives; label all comparisons
model-definition A/Bs.

- 2026-08-15 ~22:15 — Third model-def arm (Ryan-directed): **13183998
  `p018_upin_nt144_rlxf0p3`** (NT=144/PPS=3, per-step rlxf 0.3 = stock
  rate, pinned construction otherwise; 72 h/64G SETTLE 22, reaches ~rev
  14–15/30). Launcher md5 9a5b2f2c deployed local==cluster. NT=72 arms
  13183888/889 banners VERIFIED (rlxf:0.3/0.16334). **STANDING
  INSTRUCTION (Ryan): when 13183998 finishes, ASK RYAN before any `_s2`
  chain — he expects to approve it if 13183888 (NT=72 @ 0.3) looks
  good. Do not chain unprompted.**
