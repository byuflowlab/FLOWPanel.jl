# 022: Rotor Hover In Ground Effect (IGE)

**Opened:** 2026-08-18 (Ryan request, staged by agent)
**Current phase:** Phase 0 (driver + smoke) → Phases 1–2 submission
**Item-level approvals:**
- [ ] Technical Completion
- [ ] Clear-Context Approval
- [ ] User Approval After AI Discussion

## Current status and next actions

> **RESET BRIEF (b) — 2026-08-18 (submitted).** Phase 0 PASS (both smokes
> clean; GS gate decisive: 6–7 outer iters of 50, zero non-converged, wake
> direction +x confirmed by below-ground census). Four Phase 1/2 jobs
> SUBMITTED and **banners VERIFIED 2026-08-18** (all four: correct mesh,
> RPM 6000/ρ 1.16/R 0.1195, σ = 0.0400 R, d/σ = 3.4 uniform, clean shedding
> traces; IGE pair shows ground 4752 panels @ h/R 1.0, GS logging on; ige_
> coarse step-0 tangency RMS 1e-3 m/s). Harvest per `decision_rules.md` M1
> after walltime (24 h coarse / 72 h fine); headline =
> CT_IGE/CT_OGE per rung vs momentum-theory anchor ≈1.07. Deploy note in
> `ledger.md`: cluster src deliberately NOT overwritten (local has 021 WIP);
> tuple solve! diff-verified identical.
>
> 2026-08-19: Phase 5 (h/R sweep 0.5/1.5/2.0 vs Ryan-provided experiment +
> CB theory, capped at 2.0) staged; gated on the Phase-2 h/R = 1 harvest.
>
> **RESET BRIEF (d) — 2026-08-20.** Ryan swapped Phases 3↔4 and refined the
> below-ground policy (ruling 6): linear vertical velocity cutoff within R/10
> of the ground + truncation floored at the ground. Driver implements
> `GROUND_DAMP_BAND_R` (propagate! overwrite, ground-ward damping only, with
> n_damped/n_inband diagnostics) and `GROUND_PARTICLE_POLICY=truncate`;
> ground/GS CSVs now written incrementally (walled-run lesson). Phase 3
> executing per `phase_03_particle_policy.md`: A/B "without" arm = cancelled
> 13207682's logged census (steps ≤630); new arms = p022_ige_coarse_damp
> (damp only, census still counting) and p022_ige_coarse_trunc (damp +
> truncate, production candidate), 48 h each. Watch p022_ige_fine near rev
> 17–18 (mesh-independence of the blow-up) — still the most informative
> pending observation.

> **RESET BRIEF (c) — 2026-08-19.** Phase 2 gate NOT met and the coarse route
> to it is gone: **p022_ige_coarse (13207682) blew up at rev 17.6 and was
> cancelled** (CF_x −0.074 → −371.9 in one step; ground max|U·n| 4.9 → 20.2;
> below-ground sum|Γ| 0.38 → 9.86 after a monotone 456 → 6362 census climb).
> Full signature in `ledger.md`, analysis in `phase_02_ige_first_light.md`.
> Suspect = below-ground particle accumulation under the leave-be policy, but
> CT was *stationary* right up to the event, so the contingency chain's "AND CT
> drifts" half was not met — sudden instability, not slow corruption. Mesh
> independence unknown: **watch p022_ige_fine (13207681) around its rev 17–18
> — that is the single most informative thing pending.** The gate now rests on
> the fine pair (13207679/81), both healthy as of 2026-08-19.
>
> Ops: measured ~87 s/step coarse ⇒ the 1007-step schedule needs 25–31 h, so
> the 24 h coarse walltime is undersized (13207680 will wall near rev 25 of 28)
> and no IGE run can be chained (warm-start with ground errors out by design).
> Phase 5 now sizes coarse runs at 48 h. Preliminary matched read over revs
> 7–17 (pre-blow-up): IGE/OGE = **1.045 ± 0.083**, right direction, not
> discriminating. Resolvability analysis in `phase_05_hr_sweep.md` shows the
> per-rev scatter puts h/R = 2.0 (0.4% signal) below the noise floor and 1.5
> (1.0%) at the margin — **Ryan asked to revisit CT smoothing after the runs
> land** (open question recorded there, deliberately not designed yet).
>
> Phase 5 pre-work DONE (gate still closed, nothing submitted):
> `scripts/p022_harvest.py`, `figures/fig_hr_sweep.tex` + CSVs,
> launcher arms `p022_hr05/hr15/hr20`, local h/R 0.5 & 2.0 geometry smokes
> (GS 0 non-converged; below-ground allowance invariant 2.5425R confirmed).

> **RESET BRIEF (a) — 2026-08-18 (staging).** Item staged. Driver
> `examples/rotor_hover_ground_effect.jl` forked from RHPC at commit b251071
> (RHPC itself is frozen for the live 018 campaign — never modify it). Ground
> machinery is entirely pre-existing: `FlatGround` + `FlatGroundSolver` +
> block Gauss–Seidel `solve!(bodies::Tuple, solvers::Tuple)`
> (`src/FLOWPanel_solver.jl:1454`), precedent `examples/flat_ground.jl`.

In-flight jobs:

| job | case | action on landing |
|---|---|---|
| 13207679 | p022_oge_fine (45_185_ct4, OGE, 72 h) | verify banner → harvest M1 → ledger |
| 13207680 | p022_oge_coarse (56_57, OGE, 24 h) | verify banner → harvest M1 → early CT_OGE cross-check |
| 13207681 | p022_ige_fine (45_185_ct4, IGE, 72 h) | verify banner + GS/tangency health → harvest M1 → headline ratio |
| ~~13207682~~ | p022_ige_coarse (56_57, IGE, 24 h) | **CANCELLED 2026-08-19 — blew up at rev 17.6**; no M1 window. See `ledger.md` |

## Objective and scope

Simulate the DJI-9443 rotor hovering **one rotor radius above a flat ground
plane** using existing machinery (paneled ground disc; there is no image/mirror
system anywhere in the stack, verified 2026-08-18), and report **CT with and
without the ground at matched simulation settings** on two mesh rungs.
Converge the ground representation (disc radius, panel length, truncation
radius) and settle the below-ground particle policy. Out of scope: circulation
-preserving particle reassignment (Deferred), image-system implementation
(Deferred), h/R sweeps beyond 2.0 (Phase 5 owns 0.5–2.0, gated on Phase 2;
experiment shows the effect is gone by ≈1.004 at h/R = 2.0).

## Standing rulings (binding on every phase)

1. **(Ryan 2026-08-18)** Operating point: RPM = 6000, ρ = 1.16 kg/m³,
   R = 0.1195 m, rotor plane 1.0 R above the ground. These deliberately
   override 018's 5400 / 1.179 / 0.119 — 018 CT numbers are anchors only; the
   OGE baseline must be re-run at this operating point.
2. **(Ryan 2026-08-18)** Below-ground particles: start with the simplest
   policy — leave them be — and iterate (cull, then interpolation-based
   reassignment) only if the diagnostics justify it.
   **Refined by ruling 6 after the p022_ige_coarse blow-up.**
3. **(Ryan 2026-08-18, plan review)** The rotor↔ground coupling relies on the
   existing block Gauss–Seidel outer loop in
   `solve!(bodies::Tuple, solvers::Tuple)`; **verify it converges before any
   HPC submission**. Existing machinery only — no new iteration logic.
4. **(Ryan 2026-08-18)** OGE and IGE come from the *same driver file* with
   `GROUND_ENABLE` toggled and every other knob identical — that is the
   matched-settings contract.
5. Non-018 operating point aside, simulation knobs start at the 018 production
   carrier values (see Fixed operating point below); deviations are logged in
   `ledger.md` per case.
6. **(Ryan 2026-08-20)** Refined below-ground policy, superseding ruling 2's
   pure leave-be: (1) a **linear vertical velocity cutoff** for particles
   within **R/10** of the ground (ground-ward axial component scaled by
   d/band; `GROUND_DAMP_BAND_R`), with a measured pass-through A/B (with vs
   without) since the cutoff is not foolproof; then (2) **truncation floored
   at the ground** (`GROUND_PARTICLE_POLICY=truncate`, cylinder floor raised
   to the plane). Phases 3 and 4 swapped so this lands before the ground
   convergence ladders.

## Fixed operating point and unit conversions

- RPM = 6000 ⇒ Ω = 628.32 rad/s, sec/rev = 0.01, tip speed ΩR = 75.08 m/s.
- ρ = 1.16 kg/m³; ν from NASA-paper μ: 1.69e-5/ρ = 1.457e-5 m²/s (driver default).
- R = 0.1195 m; ground plane at axial (downstream, +x for dji9443 meshes)
  distance 1.0 R below the rotor plane, normal pointing at the rotor.
- CT convention: negated axial force channel (thrust), `RotorNormalization(rho, 2R, 1)`
  — identical to 018/RHPC.
- Carrier knobs (from 018 production carrier `p018_ufront_n1_s040_visc`):
  NT = 36, OVERLAP = 2.75, P_PER_STEP = 12 ⇒ σ = (2π/36)(2.75/12) R = 0.0400 R;
  Das uniform-D with `DAS_UNIFORM_DSIGMA = 3.4`; NWAKEROWS = 1;
  MERGE_R_FACTOR = 0.0055; RELAX_RLXF = 0.3 (stock, load-bearing);
  CoreSpreading ON with β = 1e9, DynamicSFS ON; formulation = velocity;
  staged startup per the validated Item-005 recipe (SPINUP 1.5 revs @ 0.4
  start fraction, MAGVINF_PEAK = 5, ramp/hold/withdraw = 1/1.5/4 revs);
  BERNOULLI_ONLY=true; settle ≥ 20 revs (IGE recirculation may settle slower
  than OGE — judge per-rev blocks, not a fixed count).
- Ground first-light: disc radius 4 R, panel_length 0.15 R, truncation
  cylinder radius 3 R (OGE keeps 1.5 R), depth 4 R, particle policy `none`.
- Sanity anchor (direction + rough magnitude only): momentum-theory /
  Cheeseman–Bennett IGE at h/R = 1 gives T/T_OGE ≈ 1/(1−(R/4h)²) ≈ 1.07 at
  fixed power.

## Phase gates

| Phase | Deliverable | Status | File |
|---|---|---|---|
| 0 | Driver + local smoke + GS convergence gate | staged | `022_rotor_hover_ground_effect/phase_00_driver_smoke.md` |
| 1 | OGE baseline CT at new operating point (45_185_ct4 + 56_57) | staged | `022_rotor_hover_ground_effect/phase_01_oge_baseline.md` |
| 2 | IGE first light CT + headline IGE/OGE comparison | staged | `022_rotor_hover_ground_effect/phase_02_ige_first_light.md` |
| 3 | Below-ground particle policy: velocity cutoff + truncation-at-ground (order swapped with 4, Ryan 2026-08-20) | executing | `022_rotor_hover_ground_effect/phase_03_particle_policy.md` |
| 4 | Ground convergence ladders (disc radius, panel length, truncation radius) | staged (blocked on 3) | `022_rotor_hover_ground_effect/phase_04_ground_convergence.md` |
| 5 | h/R sweep (0.5, 1.5, 2.0) vs experiment + CB theory | staged (gated on Phase 2 verify) | `022_rotor_hover_ground_effect/phase_05_hr_sweep.md` |

Cross-cutting: `ledger.md` (single running results table),
`decision_rules.md` (M1 metric + acceptance thresholds),
`ops_reference.md` (HPC submission, inherited 018 ops lessons).

## Contingency chain

- GS outer loop stalls or hits the iteration cap → diagnose; fallback is the
  monolithic `BackslashCoupled` solve (exists in `src/FLOWPanel_solver.jl`).
  Do NOT submit HPC jobs until resolved (ruling 3).
- Staged-startup pulse misbehaves against the ground (tangency monitor blows
  up during ramp/hold) → reduce `MAGVINF_PEAK`; log as a ruling candidate,
  never silently change.
- Leave-be particle policy corrupts the solution (below-ground census grows
  without bound AND CT drifts) → Phase 4 cull ablation moves earlier.
- Wall jet hits the truncation radius visibly (VTK) even at 3 R → widen
  before running the Phase 3 ladder.

## Deferred

- Circulation-preserving reassignment of below-ground particles (interpolate
  Γ back onto surviving above-ground particles). Only if Phase 4's cull
  ablation moves CT̄ beyond the convergence tolerance.
- Image-system (mirror) implementation — would be new machinery; out of scope.
- h/R > 2.0 rungs of the sweep — Phase 5 now owns h/R ≤ 2.0; the experiment
  (table in `phase_05_hr_sweep.md`) shows the effect is gone by ≈1.004 at 2.0.
- Warm-start with the ground body (driver errors out deliberately; multi-body
  restart state untested).
