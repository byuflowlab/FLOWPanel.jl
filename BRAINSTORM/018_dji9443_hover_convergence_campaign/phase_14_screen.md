# Phase 14 — Screening ladders: short time-marched runs for parameter search

Opened 2026-08-05 (Ryan). Purpose: the campaign has been paying full production
price (12–48 h/job) to *explore* parameter space. Phase 14 runs short
time-marched screens so that full HPC runs become confirmations. Origin: the
never-executed S1 stage of the 2026-08-05 escalation plan
(`~/.claude/plans/work-on-018-per-squishy-flute.md`); Ryan's amendments:
**time-march every case** (no frozen-state screening), include a **mesh
ladder**, and **run the screens on HPC ignoring slot limits** — which upgrades
the screen to the PRODUCTION mesh (no mesh-transfer caveat; calibration is
purely about run length) and lets each case be one 8-rev job whose per-step
monitors contain every stage window (no warm chaining).

## Pre-registrations (written before any screen data)

1. **Evidence class**: screen results rank parameters, prune ladders, and size
   HPC requests. They can NEVER satisfy M1/M2 (≥15 settled revs, staged
   startup, HPC). Screen CT is not a CT prediction — the freestream pulse is
   zeroed (1-rev spinup only), so absolute levels are transient-contaminated;
   only rung-to-rung deltas at matched windows carry information.
2. **Stage windows**: revs 1–1, 2–3, 4–7 (rev 0 = spinup, never scored).
   A pair verdict is trusted only if the CT delta AND the Γ̄ shape are
   consistent across the last two windows (delta-vs-length stability).
3. **Calibration gate** (decides whether screening works at all): the
   calibration ladder must reproduce the sign, spatial shape, and ordering of
   ≥2 verdicts already measured at full length on the same mesh:
   σ-step (scr_b0 vs scr_ov1p0 → expect CT up + 3-lobe with mid-outboard dip),
   N=1 vs N=2 (→ expect inboard-localized deficit at N=1), dt pair
   (scr_b0 vs scr_nt72 → expect small uniform positive shift). PASS ⇒ screens
   may steer ladder pruning. FAIL ⇒ record which observable does not transfer
   (a run-length result in its own right); screens restrict to
   plumbing/cost/tripwire checks.
4. **σ-collapse horizon**: ignition was at step ~1250; the 8-rev screen
   (~324 steps at NT=36) cannot see it. `scr_ufsig_0p02` reads the
   CONTRACTION RATE (`min_sigma_ratio` slope, WakeHealthMonitor) — a rate
   ranking, not a stability certificate.
5. **Mesh-ladder anchor**: chordwise convergence was already demonstrated at
   full length on HPC (2c/2d ladders to n_airfoil=249 ⇒ production
   45_185_ct4); spanwise is untested (Phase 10). The screen ladder brackets
   the production point on both axes at the ufront carrier. All rungs passed
   `scripts/p018_mesh_profile.py show` (max chord ≈ 1.328 = stock, no fold).

## Machinery

- Runner: `examples/run_p018_screen_hpc.slurm.sh` (one case/job; production
  mesh default; schedule zeroed; NREVS=8; unconditional case-arm exports;
  banner echo; WakeHealthMonitor on by default). 64 cores / 64G / 12 h
  (24 h for `scr_ufdt_nt144` and `scr_mesh_s80`).
- Scorer: `scripts/p018_screen_score.py` (`pair` / `health` / `ladder`) —
  reuses `p018_analyze.py` internals; validated by reproducing `m2` exactly
  (ε 9.882%/3.106% on `dasc0p41`→`dasc0p82` revs 4–7).

## Case roster (submitted 2026-08-05 evening; job IDs in ledger)

| ladder | cases | knobs axis |
| --- | --- | --- |
| Calibration | `scr_b0`, `scr_ov1p0`, `scr_nrows1`, `scr_nrows2`, `scr_nt72` | known HPC verdicts at screen length |
| A ufront map | `scr_uf_n1_d{2p6,3p4,5p0}`, `scr_uf_n2_d{6p5,8p5}` | D × N at σ=0.03R (N=2 floor D≥6.45) |
| B dt (clean) | `scr_ufdt_nt{36,72,144}` | first dt ladder with d_front pinned by construction (σ 0.0291R, PPS 12/6/3, rlxf halving) |
| C σ ladder | `scr_ufsig_0p05`, `scr_ufsig_0p02` (+ `scr_uf_n1_d3p4` as 0.03R rung) | σ at ufront N=1 D=3.4; 0.02R = σ_eq rung |
| D mesh | `scr_mesh_c145`, `scr_mesh_c249`, `scr_mesh_s60`, `scr_mesh_s80` (+ `scr_uf_n1_d3p4` anchor) | chord 145/185/249 at span 45; span 45/60/80 at chord 185 |

19 jobs. Ryan waived the slot cap for this batch (2026-08-05); the 10-cap
resumes for non-screen study jobs.

## Log

### 2026-08-05 — staged
Harness written; scorer validated; mesh gates run; submission + banner
verification recorded in the ledger as they happen.

### 2026-08-06 (early AM) — first harvest: 8 complete, 6 IGNITED; calibration gate FAILS; hedge triggered

**Completions/deaths (sacct + tripwire triage; all six deaths are IGNITIONS,
not memory — recorded max|u| 13k–198k m/s at last write, and FIVE show
NEGATIVE min_sigma_ratio (−1.4 … −38), the inviscid σ(1−ΔtZ)<0 failure mode
predicted by the code audit and never before observed live):**

| case | state | last step | tripwire at last write |
| --- | --- | --- | --- |
| scr_b0 / nrows1 / nrows2 / ov1p0 | COMPLETED | 323 | min_sr 0.35/0.38/0.38/0.18, max_u ≤ 50 |
| scr_uf_n1_d3p4 | COMPLETED | 323 | min_sr 0.082, max_u 27 — survivor at the boundary |
| scr_uf_n2_d6p5 | COMPLETED | 323 | min_sr 0.118 |
| scr_ufsig_0p05 | COMPLETED | 323 | min_sr 0.098 |
| scr_mesh_c145 | COMPLETED | 323 | min_sr 0.101 |
| scr_uf_n1_d2p6 | FAILED | 317 | **min_sr −38.1, max_u 197,951** |
| scr_uf_n1_d5p0 | OOM | 272 | **min_sr −7.6, max_u 13,243** |
| scr_uf_n2_d8p5 | OOM | 242 | min_sr +0.053, max_u 17,813 (igniting) |
| scr_ufdt_nt36 | OOM | 284 | **min_sr −1.4, max_u 17,204** |
| scr_ufdt_nt72 | OOM | 490 | **min_sr −5.7, max_u 32,870** |
| scr_ufsig_0p02 | OOM | 188 | **min_sr −2.4, max_u 23,062** |

Reading: at ~σ/R 0.029–0.030 INVISCID under the screen's harsh startup
(no freestream pulse, 1-rev spinup), the stability margin is razor-thin —
OVERLAP 2.4/PPS 14 (σ/R 0.0299) survives 324 steps while OVERLAP 2.0/PPS 12
(σ/R 0.0291) ignites by step 284; D=3.4 survives while D=2.6 and D=5.0 die at
the same σ. σ/R 0.0199 ignites fastest (step 188). All Campaign B (clean dt
ladder) rungs died ⇒ no dt screen data. **Caveat: screen startup is harsher
than production staged startup — these deaths bound the screen condition, not
production; the production inviscid pair at 0.0299R is alive at step ~410
(min_sr 0.145 and falling) and the viscous twin holds min_sr 0.381.**

**Calibration gate: FAIL (pre-reg 3).** σ-step pair: dip location (r/R 0.76)
and tip sign reproduce at revs 4–7, but dCT is wildly window-unstable
(−0.97% → −15.8% → +5.6%); N-axis: wrong-sign null (nrows2→nrows1 +0.07% vs
settled −0.75% with inboard deficit) AND a false positive (b0→nrows2 −3.4%
growing vs settled +0.01% null). **Consequence per pre-registration: screen
aerodynamic deltas are NOT trustworthy at 8 revs; screens are restricted to
stability/plumbing/cost probes.** The still-running mesh/nt cases will be
recorded but carry this caveat. The stability probe role, meanwhile, paid for
the whole batch in one night (six tripwire-instrumented ignition datasets).

**Hedge executed (Ryan pre-authorization, trigger met):**
`p018_ufront_s035_visc` = job **13058988** (σ/R 0.0349 via OV 2.4/PPS 12,
N=1, **D=3.0** — offline-verified d/σ exactly 3.0 uniform, Das/c 0.395–1.448
in band; D=3.4 would breach the 1.5c tip cap at this σ), viscous, cold,
SETTLE=22, 64G/48 h. Banner verified.

Also harvested: **`p018_L1_visc` COMPLETED** (22:58 h) — settled window
25–39: CT̄ 0.072391 (per-rev std ±0.00106, rippled). **Viscosity delta at
L1 σ (2.2×σ_eq): NULL in circulation** — vs inviscid L1: ε_Γ 0.156%/0.082%,
ΔCT −0.38% (the earlier "trending positive" was per-step scatter). The
viscous σ-pair b0_visc→L1_visc (+1.06%, ε_Γ 8.73%) reproduces the inviscid
3-lobe ⇒ the σ-axis error mode is viscosity-independent. Budget Γ6(L1) =
0.16%.
