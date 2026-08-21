# Phase 16 — Chord–σ Co-Scaling (Das as a non-conflicted choice)

**Opened:** 2026-08-14 (Ryan-directed after the fig_das review: "none of these
look great... can you think of a better way to choose Das?"). Status: staged →
implementing.

## Motivation: why every Das law so far failed

Das is doing three jobs at once, each with a different natural length unit:

1. **Kutta/attachment buffer** — the doublet sheet must be long enough to
   carry the shed circulation before handoff. Natural unit: **local chord**
   (014's admissible band ≈ 0.25–1.5·c_local; violating it from below is the
   inboard Γ̄ deficit).
2. **Particle clearance** — the first particles must sit far enough from the
   body that the regularized kernel does not corrupt the BC. Natural unit:
   **σ** (Phase 12 C1 band d/σ ≈ 2.6 median / 3.4 p95).
3. **Frozen near-wake extent** — the rigid straight sheet of length N·Das
   replaces wake that is actually curving. Natural unit: **local wake
   curvature radius** (the outboard deficit at Das = 1.64c, N-independent).

Every one-parameter law tried satisfies at most one job somewhere on the span:

- **η-kinematic** (`Das = η·Δt·Ωr`): holds Das/σ_trailing constant (job 2)
  but slides an r-wedge through the chord band — the admissible η windows for
  r/R = 0.31 and 0.95 do not overlap (phase_02 §2026-08-04 (b)). Non-monotone
  ladder, two opposed Γ̄ lobes, no plateau.
- **κ chord-proportional** (`Das = κ·c_local`): holds Das/c constant (job 1)
  but σ_trailing ∝ Ωr·Δt does not follow chord, so the κ ladder is collinear
  with a clearance ladder (phase_13 §(e)). Monotone but growing deltas
  (+0.81 → +1.80%), formulation-independent (Green pair). FAIL.
- **uniform-D** (`Das_j = D·σ − (N−1)·travel_j`): holds clearance span-uniform
  (job 2) but ignores jobs 1 and 3; carrier carries an unattributed +2–3%
  interaction offset, and the etadas A/B showed the Das law redistributes
  Γ̄(r/R) at fixed CT (ε_max 4.05%).

**The root inconsistency is σ's spanwise law, not Das's.** The chord band and
the clearance band conflict *only because* `OverlapPPS`/travel-referenced σ is
∝ Ωr·Δt (an artifact of tying σ to per-step TE travel — nothing physical picks
that wedge) while chord falls outboard. The ladders do not converge because
there is no small parameter: the frozen-sheet error grows with Das while the
attachment/transfer error grows as Das shrinks.

## The law

Co-scale the whole near-wake discretization with the local physical scale:

$$\sigma_j = s^{*} \, c_j \qquad |\mathrm{Das}|_j = \lambda \, \sigma_j = \lambda s^{*} c_j$$

Then Das/c is span-uniform (job 1), Das/σ is span-uniform (job 2), both are
dt-independent by construction, and the previously incompatible admissibility
windows are satisfied simultaneously by a single dimensionless pair (s*, λ).
This is the panel-method instinct applied to the wake: nobody meshes the blade
with span-uniform absolute panel size; the shedding currently does exactly
that in σ while varying it in Das (or vice versa).

**The λ window exists — the structural point.** At fixed s*, the chord band
0.25 ≤ λs* ≤ 1.5 and the clearance band λ ≥ 2.6 give

$$\lambda \in \left[\max(2.6,\ 0.25/s^{*}),\ 1.5/s^{*}\right],$$

non-empty iff s* ≤ 0.58. At s* = 0.313 the window is **λ ∈ [2.6, 4.8]** —
the η and κ ladders never had a non-empty joint window at any parameter value.

## Parameter choices (pre-registered)

- **s\* = 0.313** ⇒ σ(0.75R) = 0.0400R exactly — matched to the production
  candidate `p018_ufront_n1_s040_visc`, so the λ = 3.4 rung is a *pure
  spanwise-co-scaling A/B* against it (same σ at 0.75R, same d/σ = 3.4, same
  carrier knobs; the only change is σ and Das becoming ∝ c_local).
  Chord anchors (mesh-profile derived, verify at implementation and record
  exact values in the log): c(0.75R) ≈ 0.1277R, c_root ≈ 0.265R,
  c_tip ≈ 0.0723R ⇒ σ ranges ≈ 0.0226R (tip) – 0.083R (root).
- **λ ladder = {2.4, 3.4, 4.8}**: doubling pair 2.4→4.8 for M1's convention;
  3.4 = production-carrier match; 2.4 probes just below the C1 median (its
  marginality is disclosed, and its failure mode is pre-registered below);
  4.8 sits at the chord band's upper edge (Das/c = λs* = 0.75 / 1.06 / 1.50).
- **Carrier**: the viscous production-candidate family otherwise unchanged —
  N = 1, NT = 36, rlxf 0.3 (stock; load-bearing per the rlxf-down closure),
  CoreSpreading, merge_r 0.0055 m, OVERLAP 2.75 for pps derivation,
  SETTLE 22 / 30 revs, same mesh as `p018_ufront_n1_s040_visc`.
- Both shed families (trailing + unsteady) take σ from the station law.

## Curvature-cap diagnostic (job 3, diagnostic-only this phase)

Per station, report the arc subtended by the rigid extent on the local helix
circle: **θ_j = N·|Das_j| / r_j** (radians), table + max, in the banner and
metadata. NOTE the sign flip vs the η law: under chord-scaling the binding
station moves **inboard** (large c, small r — θ can reach ~1 rad at r/R ≈ 0.3
for λ = 4.8) while the tip is easy (θ ≈ 0.1 rad). Γ̄ is small inboard, which
may forgive it — that is exactly what the diagnostic will show. If the λ = 4.8
rung develops an *inboard* deficit that correlates with θ_j, the curvature cap
graduates from diagnostic to active law (fallback F1).

## Stability risk (the #1 threat) and screen

σ_tip ≈ 0.0226R is **below the 0.030–0.035R that ignited every uniform-σ
viscous run** (phase_13/ledger: 0.03R ×4 and 0.035R ×1 all ignited; 0.04R
survives), and 019 established in-band failure is tail-driven ignition at the
tip. Co-scaling makes σ smallest exactly where Γ concentrates. Unknown whether
the uniform-σ threshold transfers (surrounding σ distribution differs), so:

- **Stability screen GATES the ladder**: short runs (≈8 revs, NT = 36) at
  λ = 2.4 and λ = 3.4, floor off, 019 tripwire (min_sr trend, max_u, σ<0,
  |Γ|/σ²) — wake-health TAIL check, never exit-0.
- **Contingency knob `SIGMA_FLOOR_R`** (σ_j = max(s*·c_j, floor·R), default
  off): if the tip ignites, re-screen at floor 0.030R then 0.035R and carry
  the floored stations as a disclosed deviation from pure co-scaling (the
  floor re-introduces σ/c non-uniformity only outboard of where s*·c_j <
  floor·R — record the crossing station).
- Implementation checks required before submission: `CORE_SPREADING_SGM0`
  reference vs per-station σ (scalar sgm0 must not clamp/reset tip cores
  wrongly); absolute merge_r 0.0055 m ≈ 0.046R against tip σ 0.0226R
  (merge radius > 2σ_tip — verify merging cannot fuse the tip filament;
  if in doubt run the λ=3.4 screen with merge on AND off for 3 revs).

## Cases

| tag | λ | purpose | length | gate |
| --- | --- | --- | --- | --- |
| `p018_scr_cs_l2p4` | 2.4 | stability screen, smallest clearance | ~8 revs | — |
| `p018_scr_cs_l3p4` | 3.4 | stability screen, production match | ~8 revs | — |
| `p018_cs_l2p4` | 2.4 | ladder bottom / doubling pair | 30 revs | screens pass |
| `p018_cs_l3p4` | 3.4 | ladder mid + pure co-scaling A/B vs `ufront_n1_s040_visc` | 30 revs | screens pass |
| `p018_cs_l4p8` | 4.8 | ladder top / doubling pair / chord-band edge | 30 revs | screens pass |

All with s* = 0.313, floor off unless the screen forces it, diagnostics per
the current production defaults, banner-verified at submission (knobs from
log banner; judge from monitors CSV; exit 0 ≠ health).

## Pre-registered predictions & decision rule

1. **P1 (the ladder converges):** the doubling step λ 2.4→4.8 passes M1
   (|ΔCT̄| ≤ 0.5%) AND M2 (ε_Γ,max ≤ 1%) on matched 15–29 windows, where
   every previous Das ladder failed both; residual slope ~ the wing log-law
   (≈0.2%/doubling). Score BOTH raw-mean and quiet-limit M1 (phase_15
   machinery) pending Ryan's M1-observable decision.
2. **P2 (the lobes are gone):** no inboard Γ̄ deficit at λ = 2.4 (Das/c = 0.75
   is mid-band everywhere) and no outboard deficit at λ = 4.8 (Das/c = 1.50
   at the edge but span-uniform). A *new* inboard deficit at λ = 4.8
   correlating with θ_j is the curvature mechanism, not the attachment one.
3. **P3 (co-scaling vs uniform):** `p018_cs_l3p4` vs `p018_ufront_n1_s040_visc`
   isolates spanwise σ/Das redistribution at fixed everything-else. Expect a
   Γ̄ redistribution (M2 order 1–4%, cf. the etadas-null's redistribution) and
   a CT shift that reads on the carrier's +2–3% interaction offset.

**Decision:** if P1 passes, the Das axis CLOSES with Das* = (s* = 0.313,
λ = 3.4), budget term 2 = the measured slope, and the campaign's Das
prescription becomes the co-scaling law. If the ladder fails, execute
fallbacks in order (below), and the failure's lobe signature (inboard vs
outboard, θ-correlated or not) selects which.

## Fallbacks (pre-registered, in order)

- **F1 — curvature-capped law** (if failure is θ-correlated / inboard at
  large λ): `|Das|_j = min(λ·s*·c_j, β·r_j)`, β from the measured θ-deficit
  correlation; diagnostic column already in place.
- **F2 — Green transfer under co-scaling** (if failure sits at the small-λ
  edge): 015 Battery II showed Green removes 93% of the wing's Das
  differential; the rotor Green pair "failed" only under the confounded κ
  ladder. One A/B: `RHPC_FORMULATION=green` at λ = 2.4 and 4.8.
- **F3 — representation equivalence** (if Das sensitivity survives F1/F2):
  016/017's line — make the first particle row's induced field/impulse match
  the sheet segment it replaces, so handoff location loses leverage by
  construction ($dC_T/d\ln \mathrm{Das} \to 0$ rather than "small on a
  plateau"). Largest code risk; ruling-7 precedent (016 is already the σ-axis
  contingency).

## Implementation plan (knobs)

- src (`FLOWPanel_wake.jl`): new station-indexed shedding method (per-station
  σ vector + overlap, resolving to the existing `SigmaOverlap` per column in
  `_convert_to_particles!`; existing methods untouched — zero behavior change
  when unused). Unit tests: per-station σ lands on the right particles both
  families; wrap and non-wrap; exact-σ match against manual `SigmaOverlap`
  per station.
- driver (`examples/rotor_hover_pressure_comparison.jl`):
  `SIGMA_CHORD_FRACTION` (s*, NaN = off), `SIGMA_FLOOR_R` (default 0),
  `DAS_SIGMA_LAMBDA` (λ, NaN = off; Das via the existing
  `set_Das_station_lengths` path — no new src for Das). Mutual-exclusion
  errors vs `DAS_ETA_KINEMATIC` / `DAS_CHORD_FRACTION` / `DAS_UNIFORM_DSIGMA`
  and vs plain `PARTICLE_SHEDDING` σ. Banner lines (`sigma_chord:`,
  `das_lambda:`, `sigma_floor:`, θ table) + metadata fields.
- launcher: arms `p018_scr_cs_l{2p4,3p4}`, `p018_cs_l{2p4,3p4,4p8}` with
  unconditional exports (the η-inheritance gotcha, phase_03's nt72 lesson).
- Deploy: rsync + md5 verify src/driver/launcher (cluster-src lesson,
  DJI Phase 2c).

## Post-ladder rule (ADOPTED, Ryan 2026-08-15 via minimal-run plan)

If the winning rung passes settledness at 30 revs, no extension.
Otherwise warm-chain THAT ONE RUNG to 45–60 revs (`_s2` recipe); the
extended run becomes the primary hover-validation case. If burst
rectification persists, report raw mean AND quiet-limit side-by-side.
No intermediate-λ rungs to tune CT; on ladder failure route by the
registered signature to F1/F2 (below) — no densification unless
genuinely ambiguous.

## Exit criteria

Screens dispositioned (pass, or floor selected with crossing station
recorded); ladder run and scored (M1 raw + quiet-limit, M2, Γ̄ overlays into
`fig_das` as a fourth series); P1–P3 verdicts against the pre-registrations;
Das axis closed or a fallback selected with its evidence.

## Log

- 2026-08-14 — Phase opened; doc written from the fig_das review discussion
  (three-jobs/three-bands diagnosis, λ-window derivation, curvature
  diagnostic, fallback ordering per Ryan's direction). Implementation
  starting.
- 2026-08-14 — **IMPLEMENTED locally (not yet deployed/submitted).** Files:
  `src/FLOWPanel_wake.jl` (`StationSigmaOverlap` per-surface/per-station
  method; `_station_method` identity dispatch in `_convert_to_particles!` —
  existing methods bit-identical, pinned by the 016 golden-reference
  testset), `src/FLOWPanel_replay.jl` (manifest round trip),
  `examples/rotor_hover_pressure_comparison.jl` (`SIGMA_CHORD_FRACTION`,
  `SIGMA_FLOOR_R`, `DAS_SIGMA_LAMBDA` knobs; θ table; metadata fields;
  `station_chords`/`station_radii` hoisted above the shedding-method build),
  launcher arms `p018_scr_cs_l{2p4,3p4}` / `p018_cs_l{2p4,3p4,4p8}`
  (unconditional exports; screens = `P018_SETTLE_REVS=0` ⇒ 8 revs).
  Tests: wake 652/652 (new station testset 41, metadata round trip +4),
  simulate all pass (new λ·σ_j identity block), replay 125/125.
  **Chord anchors (production mesh `45_185_ct4`, driver `station_chords`,
  profile-only run):** c(0.751R) = 0.1266R ⇒ σ(0.75R) = **0.0396R at
  s\* = 0.313** — 0.9% BELOW the uniform arm's 0.0400R (disclosed, kept;
  the 0.1277R planning anchor holds to 0.7%; `p018_mesh_profile.py`'s ring
  measure reads 0.1395R — different chord definition, not used by the law).
  σ range 0.0226–0.0829R, floor binds 0/41 stations, Das/c = 1.064 uniform
  at λ = 3.4, θ_max = 1.79 rad (innermost station; θ(0.75R) = 0.18).
  **Flags:** (1) CoreSpreading `sgm0` is INERT in production (β = 1e9 ⇒
  resets unreachable; β_cur trigger is sgm0's only role) — driver now ERRORS
  on sigma-chord + CoreSpreading + β < 1e6, since a reset would stamp the
  uniform sgm0 over the per-station σ. (2) Absolute stock merge
  (`MERGE_R_FACTOR=0.0055`): r_merge/σ = 0.139 at 0.75R (stock ratio) but
  0.243 at the tip and r_merge/h_tip = 0.67 (vs stock 0.38) — tip merging is
  ~1.75× more aggressive relative to σ; watch merge counts in the screens;
  `MergeParticles(sigma_relative=true)` exists as an untouched fallback knob.
  (3) Screens' 8-rev length ends inside the withdraw transient — screens
  gate STABILITY only, never CT.
- 2026-08-14 — **Local smoke PASSED** (40_40 coarse mesh, 2 revs / 89 steps,
  NT=36, s\*=0.313 λ=3.4 N=1, viscous spreading-only, 4 threads): exit 0,
  all-finite, banner + θ table print (40_40: σ/R 0.0226–0.0829, floor 0/36,
  θ_max 1.72), metadata carries `sigma_chord_fraction/sigma_floor_r/`
  `das_sigma_lambda/theta_max`. Per-particle σ verified from the final-step
  VTP (43,343 particles): no NaN, no σ≤0, 40k distinct σ values spanning
  ~5.6× (a uniform-σ run has one shed value + growth) — the spanwise law is
  live in the pfield; exact per-station mapping is pinned particle-by-
  particle by the unit testset. Smoke CT meaningless at 2 revs (disclosed).
  VTK swept; CSVs+metadata kept at `data/smoke_cs_l3p4/`. NOT deployed, NOT
  submitted — awaiting review.
- 2026-08-14 — **Review passed (s\* = 0.313 KEPT as pre-registered, no
  retune; the 0.0396R disclosure stands). DEPLOYED + SCREENS SUBMITTED.**
  Deploy: 4 runtime files rsync'd to the cluster checkout, md5-verified
  local==cluster: wake `12c5e6f8…e25`, replay `20ed0def…d51`, driver
  `87d0b120…b3d`, launcher `8fd6beb7…66e`. (Pre-overwrite check: cluster
  copies were an older local deploy state — the known 2026-08-12 four-file
  divergence; local supersedes, in-flight jobs unaffected.) Submissions
  (24 h/64G, `--export=ALL,P018_SETTLE_REVS=0` ⇒ 8 revs):
  **13170749 `p018_scr_cs_l2p4`**, **13170750 `p018_scr_cs_l3p4`**, both
  RUNNING within a minute. Banners VERIFIED line-for-line: mesh 45_185_ct4,
  NT:36, rlxf:0.3, overlap:2.75 pps:12 merge_r:0.0055 settle:0, nwakerows:1,
  visc:true, das_chord:nan das_uniform:nan,
  `sigma_chord:0.313 sigma_floor:0 das_lambda:2.4|3.4`. Screens gate the
  30-rev ladder (`p018_cs_l{2p4,3p4,4p8}` NOT submitted). Watch duties for a
  later session: 019 tripwire on wake-health CSVs (min_sr trend, max_u,
  σ<0, |Γ|/σ²) + merge counts (tip-merge aggressiveness flag above);
  ignition ⇒ SIGMA_FLOOR_R contingency per this doc.
- 2026-08-14 — **INCIDENT + REMEDIATION: first screens (13170749/750) DIED
  at step 1, `UndefVarError: radius_inflation` (wake.jl:330), ~4 min each.**
  Root cause: rsyncing the LOCAL working-tree `FLOWPanel_wake.jl` carried
  the concurrent 021 session's uncommitted FMM radius-inflation coupling
  (`source_system_to_buffer!`, two sites), which depends on symbols in
  `FLOWPanel_elements_fmm.jl`/`FLOWPanel_abstractbody.jl` that are
  DELIBERATELY not deployed (ledger §2026-08-12 quarantine — note the
  quarantined set now effectively includes `_abstractbody.jl` and
  `_elements_fmm.jl`). Remediation: pre-rsync cluster state recovered
  EXACTLY from the NetApp home snapshot
  (`.snapshot/daily_...2026-08-13_23_05`; md5s match the pre-overwrite
  record: wake `a44e484e…788`, replay `84ae7c04…43f`); reconstructed
  cluster files = snapshot base + ONLY the Phase-16 hunks; a mechanical
  symbol audit of every identifier in the added lines against the cluster
  base found NO other local-only symbol (the reconstructed replay is in
  fact byte-identical to the local file — replay carried no 021 content).
  Deployed driver + launcher audited hunk-by-hunk vs the snapshot: deltas
  exclusively Phase 16 — left in place. Deployed + md5-verified both
  sides: wake `790b21cd901d475c4893917d25361edb`, replay
  `20ed0def07be994e08bfe9c02a380d51`; login-node compile check passes
  (`StationSigmaOverlap{Float64}` constructs; needs `module load julia
  python`, else an unrelated PyCallExt init error). **Resubmitted:
  13170886 `p018_scr_cs_l2p4`, 13170887 `p018_scr_cs_l3p4`** (same
  24 h/64G/`P018_SETTLE_REVS=0`), banners re-VERIFIED (sigma_chord:0.313,
  das_lambda:2.4/3.4, sigma_floor:0, carrier line identical). OPS LESSON:
  never deploy src by whole-file rsync from a worktree carrying another
  session's work — reconstruct target = cluster state + intended hunks;
  md5 the cluster copies BEFORE overwriting (that record is what made the
  snapshot recovery provable).
- 2026-08-14 ~21:30 MDT — **SCREENS LANDED: BOTH PASS ⇒ LADDER SUBMITTED
  (scheduled-resume session, brief (j)).** 13170886/887 COMPLETED
  (~9.5/10 h, 414 steps = 11.47 revs; settle 0 ⇒ 10 required + 1.47
  spinup, longer than the ~8 planned — more tail evidence, kept). Tripwire
  verdict on the full tail, both arms: **no ignition** — max_u ≤ 41 (l2p4)
  / 53 (l3p4) m/s all steps, both DECLINING in the last rev (ignition
  signature is 1e4–1e6); min_sigma ends 1.87e-3 / 1.79e-3 m (0.0156R /
  0.0149R), 19× above the 9.4e-5 viscous floor, never pinned; zero σ≤0.
  **Calibration against the healthy uniform carrier at matched cold age
  (ufront_n1_s040_visc revs 0–12, later survived 60 revs):** carrier
  rev-11 max|Γ|/σ² = 1.29e3, min_sigma 1.09e-3 ⇒ screens sit 4–18× BELOW
  the healthy reference on |Γ|/σ² (l2p4 72 peaked-and-declining; l3p4 297
  rising, growth decelerating ×2.4→×1.15 per rev) and ABOVE it on
  min_sigma; the sub-shed-σ min_sigma drift exists in the uniform carrier
  too (normal, not a station-law artifact). n_p end-of-rev tracks the
  carrier (202.9k / 215.2k vs 213.2k), same saturation shape — **no
  tip-merge fusion signature** (the r_merge/σ_tip=0.243 flag did not
  materialize as an n_p anomaly). **INSTRUMENTATION GAP (flag):**
  `min_sigma_ratio` and `p1_sigma_ratio` NaN at ALL steps in both arms —
  the monitor's reference σ is NaN under StationSigmaOverlap; tripwire
  evaluated via absolutes + calibration instead. Fix candidate (reference
  = min shed station σ) staged for AFTER the ladder — no redeploy tonight.
  Screen CT disclosed, not judged (ends inside withdraw transient): l2p4
  rev-11 0.0651, l3p4 0.0863. **Ladder submitted (pre-authorized):
  13178762 `p018_cs_l2p4`, 13178763 `p018_cs_l3p4`, 13178764
  `p018_cs_l4p8`** — 48 h/64G, `P018_SETTLE_REVS=22` (=30 revs), one case
  per job. Banner verification recorded here once RUNNING. Screens' VTK
  swept post-harvest (19.2G freed; local harvest
  `data/p018_scr_cs_l{2p4,3p4}/` monitors+CT+metadata).
- 2026-08-14 ~23:45 MDT — **Ladder banners VERIFIED RUNNING** (all three
  started ~21:25 MDT): 13178762 `das_lambda:2.4` θ_max 1.2649, 13178763
  `das_lambda:3.4` θ_max 1.792, 13178764 `das_lambda:4.8` θ_max 2.5298;
  all with `sigma_chord:0.313 sigma_floor:0`, settle:22 (1080 steps = 30
  revs), carrier line identical to the screens (NT:36 rlxf:0.3
  overlap:2.75 pps:12 merge_r:0.0055 nwakerows:1 visc:true, mesh
  45_185_ct4, depth banner 4R = 3.5R actual), cold start (no warmstart
  line, correct). No mismatch ⇒ no scancel. ETA ~48 h (2026-08-16).
  Watch note for the next session: in-flight tripwire must use min_sigma
  and max_u ABSOLUTES (sigma-ratio columns are NaN under the station law,
  see flag above); healthy-carrier calibration values in this log.
- 2026-08-17 ~16:00 MDT — **l2p4 + l3p4 chains LANDED healthy through rev
  30; INTERIM scoring (final P1/P2 verdicts await l4p8, ETA ~19:00).**
  Health: gos2 200–500 all the way (carrier envelope O(1e3)); pause
  seams CLEAN (l2p4 0.098% ≈ its typical 0.095% per-step move; l3p4
  0.037% < 0.05% gate). **Headline: co-scaling COLLAPSES the burst mode
  at NT=36** — per-rev CT std 0.000433 (λ2.4) / 0.000301 (λ3.4) vs the
  uniform carrier's 0.00303 on the same window (×7–10 smaller), block
  drift 0.18%/0.52% NON-monotone vs the carrier's +2.08% monotone;
  **λ2.4 M1 raw = PASS — the campaign's first raw-M1 pass** (λ3.4 CHECK
  at 0.52%, non-monotone). Consistent with the burst-mode-NUMERICAL
  verdict: the station-σ law removes most of that mode without touching
  dt. Matched 15–29 (stitched at 834/778): CT̄ λ2.4 = 0.071074
  [0.070937, 0.071352], λ3.4 = 0.072442 [0.072234, 0.072618] ⇒
  **Δ(2.4→3.4) = +1.92%** — the λ slope is NOT yet small, P1's
  2.4→4.8 doubling test rides on l4p8. **M2: FAIL both pairs** —
  λ2.4-vs-λ3.4 ε_max 5.97%/rms 2.40%; P3 A/B (λ3.4 vs uniform carrier)
  ε_max 6.45%/rms 3.81% (larger than the predicted 1–4%
  redistribution). Γ̄ spanwise structure analysis deferred to the full-
  ladder scoring (P2 lobe signatures). Local harvests
  `data/p018_cs_l{2p4,3p4}_rs1/`.
- 2026-08-17 ~19:30 MDT — **LADDER COMPLETE (l4p8 landed healthy, gos2
  1.0–2.7e3 tail, no ignition — the rev-19 heat subsided as recorded).
  FULL PRE-REGISTERED VERDICTS:**
  - **P1: FAIL.** Doubling λ 2.4→4.8 matched 15–29: raw CT̄ 0.071074 →
    0.077920 = **+9.63%** (quiet-limit 0.0717 → 0.0785 = +9.5% — same
    verdict in both observables); M2 ε_max **31.35%**/rms 11.28%. The
    λ axis does not converge as a scalar ladder.
  - **P2: the lobes are NOT gone, and the failure signature is the
    CURVATURE one:** Γ̄(λ4.8)−Γ̄(λ2.4) is an **INBOARD, θ-correlated,
    λ-monotone EXCESS** — +46.7% of Γmax at r/R 0.268, decaying
    monotonically to ~0 at r/R 0.76, small outboard dip (−6.0% at
    0.85); λ3.4 shows the identical shape at ~1/5 amplitude (+8.6%
    inboard). Sign disclosure: pre-registration predicted an inboard
    DEFICIT; the measured sign is EXCESS (the long rigid sheet at
    small r, θ up to 2.53 rad, OVER-carries bound circulation), but
    the location/correlation selector is unambiguous. Outboard
    (r/R > 0.65) the 2.4-vs-3.4 pair agrees to <1% — the tip side is
    already consistent between low rungs.
  - **P3 (co-scaling vs uniform, λ3.4 pure A/B):** ΔCT̄ = −3.79%
    (0.072442 vs 0.075298), M2 ε_max 6.45% — a much larger
    redistribution than the etadas-null's 4.05%, i.e. the σ spanwise
    law (not the Das law) owns the bigger lever, PLUS the burst-mode
    collapse (per-rev std ×7–10 down; λ2.4 = first raw-M1 PASS).
  - **Decision per the pre-registered routing: failure is inboard +
    θ-correlated ⇒ F1 (curvature-capped law |Das|_j = min(λ·s*·c_j,
    β·r_j)) is SELECTED by the doc — STAGED for Ryan's confirmation,
    NOT enacted** (fallback selection is a Ryan decision per brief
    (j)). Evidence for β: the excess onset sits near θ ≈ 0.5–0.7 rad
    (deviation <1% for r/R > 0.74 where θ(λ4.8) < ~0.5; grows
    superlinearly inboard of that). Note the 2.4-vs-3.4 inboard gap
    (6–8%) implies even λ2.4 (θ_inboard ≈ 0.75 rad) is not yet below
    the curvature threshold inboard — F1's cap would bind there too.
  - Quiet-limits: 0.0717 / 0.0724 / 0.0785 (λ 2.4/3.4/4.8) — raw ≈
    quiet at the low rungs (bursts gone), l4p8 bursty again (per-rev
    std 0.0017, 5.6× l2p4) — burstiness itself tracks the curvature
    violation.
- 2026-08-18 — **F1 LAUNCHED (Ryan order 2026-08-17 "launch F1").**
  Implementation: driver knob `DAS_CURVATURE_BETA` (β rad; |Das|_j =
  min(λ·s*·c_j, β·r_j/N) ⇒ θ_j ≤ β; validation: requires
  DAS_SIGMA_LAMBDA, β > 0; banner cap-count line; metadata
  `das_curvature_beta`; θ table/θ_max reflect the cap automatically),
  launcher banner `das_beta:` + arms `p018_cs_f1_l{2p4,3p4,4p8}` with
  **β = 0.6 rad** (measured excess onset θ ≈ 0.5–0.7). Gate smoke
  (coarse 40_40, λ4.8+cap, 4 threads): ran 9.0 revs/324 steps clean
  (zero σ≤0, zero NaN, tail max_u 27.7, |Γ|/σ² 16; over-long because a
  fictitious REQUIRED_REVS knob defaulted the run to full length —
  disclosed) + banner-capture rerun verified **cap ACTIVE, binds 15/36
  stations (40_40), θ_max = 0.6 rad exactly**. Smoke CSVs kept at
  `data/smoke_cs_f1_l4p8/`, VTK removed. No `src/` change (driver-only)
  ⇒ no quarantine exposure; station-σ law unchanged ⇒ no new stability
  screen (stability was σ-driven; the cap only SHORTENS inboard Das).
  Disclosed tradeoff: at capped stations the local clearance
  λ_eff = β·r_j/σ_j drops below the C1 band inboard (e.g. ~1.9 at
  r/R 0.27) — that conflict is intrinsic to adding the curvature band;
  watch for inboard BC-corruption signatures in M2. Deploy: driver md5
  `257cdc6b9fc8ff91e9a1f72764cee25a`, launcher
  `9949b29dd2339809a63162245a765eb7`, local==cluster. **Submitted:
  13193493 `p018_cs_f1_l2p4`, 13193494 `p018_cs_f1_l3p4`, 13193495
  `p018_cs_f1_l4p8`** (48 h/64G, SETTLE 22 ⇒ 30 revs). Banner-verify on
  start (incl. `das_beta:0.6` + cap-count line + θ_max 0.6). Scoring on
  landing: P1-style doubling 2.4→4.8 under the cap (M1 raw+quiet, M2),
  Γ̄ overlay vs the uncapped rungs (did the inboard excess collapse?),
  l3p4-F1 vs uncapped l3p4 (cap-only A/B at the production λ).
- 2026-08-20 — **F1 TRIO LANDED (TIMEOUT at 48 h wall, revs 28.4/28.7/
  28.9 of 30) + SCORED on matched 15–28 (uncapped rungs re-scored on the
  same window; banners post-hoc VERIFIED first — cap ACTIVE, binds
  10/13/17 of 41 production stations, θ_max 0.6).**
  - **M1 under the cap:** CT̄ 0.070539 (λ2.4) / 0.070605 (λ3.4) /
    0.069877 (λ4.8); per-rev std 2.9e-4/2.2e-4/1.8e-4 — the λ4.8
    burstiness is GONE (uncapped 1.76e-3). Doubling 2.4→4.8 = **−0.94%**
    (uncapped same-window +9.68%); neighbor 2.4→3.4 = +0.09%, CIs
    overlap. The λ axis is now flat to ~1% — a ×10 collapse of the P1
    failure.
  - **M2:** doubling ε_max 5.63%/rms 1.89% (uncapped 31.35%/11.28%);
    neighbor pair 1.50%/0.61%. Still above the ≤1% gate.
  - **Γ̄ overlay (the P2 selector):** the inboard θ-correlated excess
    COLLAPSED — doubling ΔΓ̄ is −0.6..−0.9% of Γmax across r/R
    0.11–0.74 (was +34..+54%). The residual is an OUTBOARD lobe, r/R
    0.78–0.87, peak −5.6% at 0.83, nearly identical to the uncapped
    ladder's own outboard dip (−6.2% at 0.85) where the cap does not
    bind — i.e. a distinct outboard λ-sensitivity (tip-roll-up region),
    NOT curvature residue.
  - **Cap-only A/B at λ3.4** (f1_l3p4 vs cs_l3p4, matched 15–28): ΔCT̄
    −2.56%; M2 ε_max 7.80% concentrated inboard of r/R 0.4 (−2..−6% of
    Γmax), ~0 outboard of 0.55 — the cap acts only where it binds.
  - No chains submitted (the lost rev 29 does not block any matched
    comparison). **Next decision (Ryan's): F1-vs-F1b head-to-head once
    the csarc ladder lands.**
- 2026-08-20 — **F1b screens PASS ⇒ TE ladder SUBMITTED (pre-authorized).**
  Screens 13206336/337 full-tail absolutes: zero σ≤0, max_u 14–23
  flat/declining, min_sigma 1.78e-3 (matches the passing cs screens),
  gos2 l2p4 425 rising-decelerating / l3p4 flat 81 — below the 1.29e3
  healthy-carrier calibration. Ladder = **13243083/084/085
  p018_csarc_l{2p4,3p4,4p8}**, 48 h/64G, SETTLE 22, banners VERIFIED
  RUNNING (das_arc:true, steady table p018_cs_l3p4_rs1_te_downwash_te.csv,
  λ per case, no cap). Watch note: csarc_l2p4's rising gos2 tail from the
  screen — judge the 30-rev tail before scoring. (Take-1 submission
  13243011/12/13 died in 13 s: quoting bug ate the case arg; no side
  effects.)
- 2026-08-20 — **F1b NT ladders LAUNCHED (Ryan order; can't wait for the
  λ ladder).** Six new arms csarc_nt{72,144}_l{2p4,3p4,4p8} =
  13245449–454, 72 h/64G, SETTLE 22; NT·pps = 432 held fixed; rlxf by
  the exact-rate rule r(NT)=1−(1−0.3)^(36/NT) (Ryan standing ruling
  2026-08-20). Banners verified on start. Together with the live NT=36
  λ ladder this forms the full 3λ×3NT F1b matrix; NT rungs score as
  labeled model-def A/Bs (dt changes shedding-linked properties — Ryan).
  nt144 rungs reach only ~rev 13–15 per wall; settled windows there
  need chains (Ryan's call when they land).
