# Phase B — HPC per-phase timing at production scale

**Protocol (2026-08-20).** New harness `benchmark/p018_mature_wake_timing.jl`
(successor to the 2026-06 `rotor_hover_pressure_steady_aero_timing.jl`):
restores `p018_cs_f1_l3p4` at its highest retained step (~1034, rev ~28.5,
~181k particles) by mirroring `simulate_warmstart!` sections 2.5–5 —
including the mandatory Das kinematic replay and
`_restore_wake_continuation!`, both missing from the old harness — then runs
5 measured production steps timing each stage helper of
`_steady_aerodynamics!` (`_sa_wake_influence!`, `solve_formulation!`,
`_sa_body_influence!`, …) plus monitors, VTK/metadata output,
`propagate!` and `shed_wake!`, then one profiled step (`Profile.@profile`,
flat report). Wrapper `benchmark/slurm/p023_mature_timing.sh` carries a
verbatim copy of the launcher's `p018_cs_f1_l3p4` env (SETTLE_REVS=22 ⇒
1081-step t_range). Scratch output goes to node tempdir; `data/` is read-only.

**Validation before submission:**
- Local smoke (Mac, 4 threads, 40_40 baseline restart at step 59, NWAKEROWS=4
  to match that run's construction): full pipeline ran; steady-step phase sum
  reconciled with the step wall total to <0.1 s. (The smoke's frame-deviation
  warning is expected — its env intentionally differs from the June baseline;
  the warning's ABSENCE on the production run is part of the sanity gate.)
- Cluster `src/` + `examples/rotor_hover_pressure_comparison.jl` md5-verified
  identical to local; harness + wrapper md5-verified after scp.

**Acceptance gates:**
1. Banner shows production config (mesh 45_185_ct4, NT 36, σ* 0.313, λ 3.4,
   β 0.6, settle 22) and restored particle count ≈ 181k.
2. No "replayed frame state deviates" warning (env matches parent run).
3. Per-phase sums reconcile with step totals to ±10%; step totals consistent
   with the parent run's late-window `wall_s` (~155–200 s, node variance per
   Phase A finding 4).
4. First-step thrust (bench monitor CSV CF_x) within noise of the parent's
   final steps (parent monitor02 at steps 1030–1034).

## Runs

| job | what | status |
|---|---|---|
| 13245638 | 64-thread, 5 steps + profile step, output arm on | **DONE 2026-08-20 15:38** (21 min), all gates PASS |
| 13245763 | sensitivity: FMM_LEAF_SIZE=100 (both backends), 3 steps | submitted 2026-08-20 |
| 13245764 | sensitivity: FMM_ACCEPTANCE=0.5 (both backends), 3 steps | submitted 2026-08-20 |
| 13245765 | thread scaling: 32 threads, 3 steps | submitted 2026-08-20 |

## HEADLINE REFRAME (2026-08-20, nf-split job 13245869)

`fmm!` nearfield/farfield ablation on the frozen mature state shows the wake
particle influence call (production knobs p=4/mac=0.4/leaf=50) costs
**140.8 s with the Dynamic-SFS postcalc (`Estr_fmm!`) and 7.4 s without it**
(nearfield-only 137.1 / farfield-only 4.0 with postcalc; 7.8 / 4.2 without).
So ~95% of `sa_wake_influence` — ≈133 s, ~75% of the whole production step —
is the SFS estimator walking the nearfield interaction lists, NOT the FMM
velocity evaluation (which splits roughly half nearfield / half farfield at
~4 s each). This also explains the MAC sensitivity: raising MAC shrinks the
nearfield lists Estr_fmm! traverses.

The body pass is the mirror image: **45.7 s full, 45.6 s nearfield-only,
0.35 s farfield-only** — 99% nearfield direct at p=8/mac=0.4/leaf=20.
Suspected mechanism: kerneloffset-inflated panel radii (FastMultipole warns
"radius inflation exceeds 10x the panel radius ... most interactions will be
evaluated directly" on the 40_40 smoke; the same inflation applies at
production scale), which forces pairs into the direct lists regardless of
MAC/leaf. If the tune job's body pass cannot pull this down, the lever is the
radius-inflation/kerneloffset interaction (route to Ryan), not the knobs.

Same-node MAC ladder (m12-2-8, steady-step means, thrust continuity 4-digit
in every arm): MAC 0.4 → 176.7 s (repeat arm 13245823; original-node baseline
172 s ⇒ node variance ~3% here), MAC 0.5 → 123.1 s, **MAC 0.6 → 85.6 s
(2.06×)** — monotone through the tested range; achieved-error quantification
rides on the tune job's verification stage. Leaf 100 (both backends) → 268 s
(+52%): leaves are already past the direct-work optimum in the growing
direction.

## FMM tuning stage (Ryan 2026-08-20 directive) — debugging trail

Three findings from the failed tune attempts, before the live run (13246014):

1. **Stale first wake row after warm-start restore** (root cause, localized
   by local bisection): after `warmstart_restore!`'s replayed `shed_wake!`,
   the first panel-wake row is stale/sentinel until the next step's
   `update_TE!` rebuilds it from the body TE + Das. Any influence evaluation
   that includes the panel-wake sources BEFORE `update_TE!` returns a 100%
   non-finite target field (reproduced locally at 40_40: sources, strengths,
   and kerneloffset all clean; reset clean; the panel-wake-source call alone
   NaNs all 7288 body CPs + 9378 particles). Production and the timing
   harness are immune because they always run the full step head. Initial
   "high expansion order is broken" reading was WRONG — p=16 failed for this
   same reason. Standing rule for any future harness: run maneuver! + reset/
   freestream/kinematic + `update_TE!` before any influence evaluation on a
   restored state.
2. This NaN field is what poisoned both tolerance mechanisms: FastMultipole's
   `PowerRelativeGradient` computed a NaN tolerance internally (job 13245865,
   "Error tolerance NaN not reached, eps_mp=5.3e-13"), and the absolute
   tolerance derived from the NaN reference RMS did the same (13245991). With
   NaN tolerance every comparison is false, so tune_fmm silently rejects all
   MACs and returns its defaults with cost=Inf — a misleading failure mode.
3. Fix: the achieved-error reference is now `DirectBackend` (exact, feasible
   at this scale: ~4e10 pair interactions on 64 cores), and the tolerance is
   absolute = rel-target x RMS(direct reference) per the 021 convention.
   Cluster `tune_fmm` (older FastMultipole rev a9b734a) returns a 2-tuple,
   not the local 3-tuple — the driver indexes rather than destructures.

`FastMultipole.tune_fmm_perturb` (new, reusable, appended additively to the
CLUSTER's `FastMultipole/src/autotune.jl` and in the local copy): greedy
one-at-a-time descent (P±1, MAC±0.05 in [0.25,0.85], leaf x/÷1.5) on measured
`fmm!` wall time, min-of-reps, memoized, error-tolerance-guarded, stops when
no neighbor improves ≥2%.

## FMM tuning RESULTS (job 13246037, DONE 2026-08-20 18:22)

Mature state (181,307 particles @ step 1034), 64 threads. Reference =
DirectBackend exact evaluation of the wake-induced field (RMS 5.64 m/s wake
pass / 15.31 m/s body pass); tolerance = rel-target x that RMS; cost = the
full production influence phase (reset + update_TE! + influence incl. SFS
postcalc for the wake pass), min-of-2. Data:
`benchmark/results/p023_fmm_tune.csv` (+ `_history.csv` for every descent
candidate). Zero non-finite columns this run.

| pass | variant | p | MAC | leaf | phase t [s] | rel err | meets |
|---|---|---:|---:|---:|---:|---:|---|
| wake (1e-4) | production | 4 | 0.4 | 50 | 139.2 | 2.1e-4 | **NO** |
| wake | tune_fmm | 15 | 0.6 | 36 | 47.1 | 3.3e-6 | yes |
| wake | **perturbed min** | 16 | 0.6 | 24 | **37.8** | 1.7e-5 | yes |
| wake | p16/mac0.3 diag | 16 | 0.3 | 50 | 252.5 | 8.8e-13 | yes |
| body (1e-5) | production | 8 | 0.4 | 20 | 45.2 | 5.9e-10 | yes |
| body | tune_fmm | 15 | 0.8 | 531 | 64.1 | 2.8e-9 | yes |
| body | **perturbed min** | 15 | 0.8 | 236 | **36.3** | 3.5e-9 | yes |
| body | p16/mac0.3 diag | 16 | 0.3 | 50 | 66.1 | 5.4e-13 | yes |

Findings:

1. **Production wake knobs MISS the 1e-4 relative target** (2.1e-4) while
   costing 3.7x the tuned point: p=16/MAC=0.6/leaf=24 is simultaneously
   3.7x faster (139->38 s) and 12x more accurate. Most of the saving is the
   SFS postcalc shrinking with the MAC-0.6 direct lists (the pure fmm!
   velocity trial at these knobs is ~1.7 s).
2. **The perturbation stage earns its keep**: on the body pass, tune_fmm's
   cost model mispredicted leaf (531 -> 64.1 s, SLOWER than production);
   the measured descent corrected it to leaf=236 -> 36.3 s. On the wake pass
   it shaved a further 20% off the tuned point (47 -> 38 s).
3. **The body pass is error-saturated, cost-floored**: production error
   5.9e-10 has ~4 orders of headroom against the 1e-5 target, yet cost only
   drops 45 -> 36 s (MAC 0.4 -> 0.8!) because kerneloffset-inflated radii
   force near-field direct pairs regardless of MAC. The remaining ~36 s is
   NOT reachable through knobs — the lever is the radius-inflation /
   kerneloffset interaction (design question for Ryan).
4. Projected production step at the perturbed knobs: wake 38 + body 36 +
   solve 16 + I/O+misc ~4 ~= **~94 s vs 172 s (≈1.8x)**, with wake accuracy
   improved from out-of-tolerance to 1.7e-5. A confirming full-step A/B
   (5-step timing at the perturbed knobs + CT check) is the natural next
   measurement; NOT yet run.

## Results — job 13245638 (production config, 64 threads)

**Gates:** banner correct (mesh 45_185_ct4, NT 36, σ* 0.313, λ 3.4, β 0.6,
settle 22 ⇒ 1080-step schedule); restored 181,307 particles; **no**
frame-deviation warning; phase sums reconcile with step totals to <0.1 s
(after step 1's 3.9 s JIT residue); thrust continuity CF_x −0.070110 (parent
step 1034, monitor02) → −0.07012 (bench step 1035) — 4-digit agreement.
Restore is nearly free: kinematic replay 0.1 s, state load 1.5 s (RHPC setup
incl. dense-G factorization 70.8 s).

**Per-phase split** (mean over steady steps 1036–1039; step totals 142–250 s,
mean 172 s — matching the parent run's late-window `wall_s` 155–200 s):

| phase | mean s | % of step |
|---|---:|---:|
| sa_wake_influence (wake+particle sources → body + 181k particle targets) | 110.5 | 64.2 |
| sa_body_influence (36,752-panel body → body + particle targets) | 43.5 | 25.3 |
| sa_solve (dense Backslash re-solve) | 16.1 | 9.3 |
| output_wake_vtk + body_vtk + metadata | 1.8 | 1.1 |
| everything else (monitors, propagate, shed, TE, resets) | <0.3 | <0.2 |

Two structural surprises vs the 2026-06 40_40/39k measurement (which had
wake 87.9 / solve 6.3 / body 5.3): body influence is now a quarter of the
step (particle targets × 8.5× more panels), and `propagate_total` is
~0.06 s — the VPM integrator is trivial because ALL particle velocity
evaluation happens inside the two influence phases.

**Profiler (flat, one full step):**
1. **Thread utilization is 49%** across 64 threads — half the machine idles
   in `poptask`/`wait` during an average sample. Parallel efficiency, not
   arithmetic, is co-limiting.
2. The busy time is dominated by **near-field direct kernels**, not far-field
   expansions: `induced`/`_induced` panel evaluations
   (`FLOWPanel_elements_fmm.jl`) routed through `direct!` +
   `FastMultipole.execute_assignment!` nearfield lists, then FLOWVPM
   particle-particle `direct!` with `Estr_direct` (Dynamic-SFS) attached.
   ⇒ leaf-size / multipole-acceptance retuning at production scale directly
   attacks the dominant cost; hence the sensitivity jobs above.
