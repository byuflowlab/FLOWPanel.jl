# 025 Phase 3 — 018 production compatibility (warm-start A/B on the frozen carrier)

**Status:** COMPLETE 2026-08-21. Two Ryan rulings the same day, both
implemented: 018 production runs **Gaussian**, and at the **tuned FMM knobs**
(wake 16/0.6/38, body 17/0.7/109). Net production step 178 → 71.7 s with CT
unchanged to 0.03%.

## Question

Phases 0–2 chose Gaussian on kernel peaks and certified it on a 5-step
warm-start bench. Two things that bench cannot answer:

1. **Physics compatibility.** Does swapping the filament regularization family
   move CT on the frozen 018 production carrier — i.e. is the campaign's CT
   ladder still comparable across the switch?
2. **Real-pipeline cost.** The bench timed a stripped harness. What does the
   actual RHPC driver (monitors, VTK, merging, relaxation, I/O) cost per step
   under each family?

## Method

Two RHPC-driver warm-start continuations of the production carrier
`p018_cs_f1_l3p4` (mesh `45_185_ct4`, velocity formulation, NT=36, 4R,
rlxf 0.3, N=1, σ*=0.313, λ=3.4, β=0.6, viscous), each resuming from the
parent's step 1034 and marching to step 1151 — **117 steps ≈ 3.25 revs** on an
already-mature ~181k-particle wake.

| arm | job | family pin | node | elapsed | exit |
|---|---|---|---|---:|---|
| control | 13247862 | `FLOWPANEL_FILAMENT_REG=vatistas` | m12-1-3 | 05:18:17 | COMPLETED 0:0 |
| treatment | 13247863 | `FLOWPANEL_FILAMENT_REG=gaussian` | m12-1-7 | 04:40:34 | COMPLETED 0:0 |

Everything else — mesh, schedule, knobs, seed state, `SETTLE_REVS=24`,
production FMM knobs (body 8/0.4/20, wake 4/0.4/50) — identical. Scored from
the monitor CSVs, not the driver's summary block (see *Reporting caveat*).

### Family provenance

Neither log printed the family (fixed since: `FLOWPanel.__init__` now prints
the pinned family, and `case_metadata.toml` records
`filament_regularization`). The arms are nevertheless separable after the
fact, from the FMM radius-inflation warning in `*.err`:

| arm | reported Δr | Δr / r_c | expected rule |
|---|---:|---:|---|
| vatistas | 0.03760603 m | 37.6 | $r_c(2/\mathrm{tol})^{1/4}$, tol $10^{-6}$ |
| gaussian | 0.00589798 m | 5.90 | gradient-aware $e^{-z}(1+2z)=\mathrm{tol}$ |

Both match Phase 1 exactly ⇒ the arms really ran different families.

## Result 1 — CT compatibility (the headline)

$C_T = -C_{F_x}$ from `monitor02_force_system1.csv`, per 36-step revolution:

| steps | parent (vatistas, pre-handoff) | vatistas arm | gaussian arm | gau − vat |
|---|---:|---:|---:|---:|
| 927–962 | 0.0704647 | — | — | — |
| 963–998 | 0.0704548 | — | — | — |
| 999–1034 | 0.0704314 | — | — | — |
| 1035–1070 | — | 0.0703744 | 0.0703765 | +2.1e-6 (+0.003%) |
| 1071–1106 | — | 0.0704469 | 0.0704698 | +2.29e-5 (+0.033%) |
| 1107–1142 | — | 0.0703522 | 0.0703760 | +2.38e-5 (+0.034%) |
| **1035–1151 (all)** | — | **0.0703810** | **0.0703975** | **+1.65e-5 (+0.023%)** |

**Compatibility number: +0.023%** — CT identical to 3.5 significant figures.
Three things make that a null rather than a small effect:

- within-rev peak-to-peak is 4.2–6.0e-4 (0.6–0.85% of CT), so the family delta
  is ~3% of one revolution's own oscillation;
- the *parent's* rev-to-rev drift over its last three revolutions is −3.3e-5 —
  **twice the family delta**, in the same code, same family;
- the handoff is continuous: parent step 1034 $C_{F_x}=-0.0701100$ →
  gaussian step 1035 $-0.0701215$, and both arms continue the parent's mild
  downward drift with no discontinuity.

⇒ The Gaussian default is physics-compatible with the 018 carrier. Campaign
rungs run before and after the switch remain comparable; the family carries no
CT bias worth a term in the error budget.

## Result 2 — cost through the real driver

`wall_s` from `monitor04_wake_health_system1.csv`, mature steps 1045–1151
(n=107; step 1035 is the restore step and logs NaN by construction):

| arm | mean s/step | median | min–max | time-marching wall |
|---|---:|---:|---:|---:|
| vatistas | 160.2 | 152.7 | 112.1–223.3 | 18880.4 s |
| **gaussian** | **140.2** | 133.6 | 93.8–223.3 | 16586.1 s |
| parent, steps 927–1034 (same carrier, older src) | 178.1 | 178.4 | 134.4–222.3 | — |

**Real-pipeline speedup 1.14×** (mean and median; whole-run 1.138×), i.e.
~20 s/step. Cross-check against the Phase 2 bench at production knobs
(Vatistas 183.4 s, Gaussian 141.8 s): **the Gaussian arm reproduces the bench
to 1.1%**. The Vatistas arm is 13% below the bench control and 10% below the
parent's own 178 s — node/state variance, which 023 measured at ±3–30%.

**Separating the two levers.** The 2.4× headline (72.2 s/step) is Gaussian
**plus its own tuned FMM knobs**. At the time of this A/B the launcher shipped
production knobs, so the family alone buys ~1.14×; the knobs are the larger
half and are scored separately in Result 5.

## Result 3 — monitor and tripwire health (steps 1035–1151)

- No tripwire fired in either log; `all_finite = true`.
- `n_particles`: vat 180 885–182 546, gau 180 800–182 525 (Δ < 0.1%).
- `min_sigma`: both 1.0664e-3 → 7.947e-4 — identical trend, shared
  merge/stretch history, no negative-σ events.
- `max_gamma_over_sigma2`: both 298.0 → ~177, monotone-ish **down** — no
  Γ-ignition in either arm.
- **`max_u`: vatistas peaks at 31.24 m/s, gaussian at 24.27 m/s (−22%).**
  This is the Phase 0 matched-$r_c$ peak ranking (0.45 vs 0.71) showing up in a
  live wake — the conditioning argument now has in-situ evidence, not only
  closed-form peaks.

## Result 4 — artifact-chain integrity

Parent force CSV ends at step 1034; both arms start at 1035 — no gap, no
overlap. `restart_step = 1034` in both `case_metadata.toml`. Full artifact sets
present (CT_vs_rev, CT_per_rev, case_metadata, run metadata TOML, body / wake /
filament / particle pvd + VTK, four monitor CSVs). 117 timesteps of VTK
retained per arm, 4.3 G each — the 36-step retention sweeper has not been run
on these two directories (8.6 G recoverable).

## Reporting caveat (applies to every warm-start continuation)

The driver's headline convergence block is **invalid for continuations**:
`CT_vs_rev.csv` zero-fills the 1034 restored steps, so revolution blocks 1–28
read 0.0, block 29 is half-zero, and the 10-rev window reports
`CYCLE-MEAN CT = 0.02288 ± 145%`, `CONVERGED = false`. That is arithmetic on
zeros, not divergence. Score continuations from per-rev blocks 30–32 or
directly from the monitor CSVs, as done above.

## Changes made under this phase

- `examples/run_dji9443_hover_ct_hpc.slurm.sh`: default
  `FLOWPANEL_FILAMENT_REG` **vatistas → gaussian** (Ryan ruling 2026-08-21);
  banner now echoes the family and both FMM knob triples. Pre-2026-08-21 arms
  are reproduced by pinning `vatistas` at submission.
- `examples/rotor_hover_pressure_comparison.jl`: per-pass FMM knobs
  `FMM_BODY_{EXPANSION_ORDER,ACCEPTANCE,LEAF_SIZE}` and `FMM_WAKE_*`, falling
  back to the shared `FMM_*` names (behavior unchanged when only the shared
  names are set); `case_metadata.toml` now records both knob triples and
  `filament_regularization`.
- `src/FLOWPanel.jl`: `__init__` prints the pinned family for log provenance.

## Result 5 — tuned FMM knobs (job 13290979, Ryan follow-on)

Same warm-start continuation, same seed state (step 1034), same 117 steps,
Gaussian family, but at the tuned knobs (**wake 16/0.6/38, body 17/0.7/109**)
instead of production knobs. Family held fixed ⇒ the comparison against the
`reg_gau` arm isolates FMM accuracy alone. COMPLETED 0:0 in 02:26:28.

| steps | gaussian, production knobs | gaussian, tuned knobs | tuned − production |
|---|---:|---:|---:|
| 1035–1070 | 0.0703765 | 0.0703738 | −2.69e-6 (−0.0038%) |
| 1071–1106 | 0.0704698 | 0.0704720 | +2.20e-6 (+0.0031%) |
| 1107–1142 | 0.0703760 | 0.0703801 | +4.14e-6 (+0.0059%) |
| **1035–1151 (all)** | **0.0703975** | **0.0704002** | **+2.70e-6 (+0.0038%)** |

Largest single-step divergence over the 117 steps: 3.73e-5 (0.053% of CT),
mean 7.8e-6. So 117 steps of accumulated FMM error move CT by **+0.0038%** —
an order of magnitude below the family null (+0.023%), and the worst *instant*
is a fifteenth of one revolution's own peak-to-peak. The static certification
(1.8e-6 wake / 1.2e-6 body vs a DirectBackend reference, job 13247200) is
confirmed to hold dynamically: the tuned knobs do not accumulate.

Cost and health, mature steps 1045–1151 (n=107):

| arm | knobs | mean s/step | median | max_u | n_particles | min_σ | max Γ/σ² |
|---|---|---:|---:|---:|---|---:|---:|
| vatistas | production | 160.2 | 152.7 | 31.24 | 180 885–182 546 | 7.947e-4 | 289.8 |
| gaussian | production | 140.2 | 133.6 | 24.27 | 180 800–182 525 | 7.947e-4 | 289.8 |
| **gaussian** | **tuned** | **71.7** | **72.7** | 25.34 | 180 830–182 521 | 7.940e-4 | 289.8 |

**2.23× vs the Vatistas control, 2.48× vs the parent run's own 178.1 s/step**,
and 71.7 s reproduces the Phase 2 bench (72.2 s) to 0.7%. Wake health is
indistinguishable across all three arms — identical particle counts to 0.1%,
identical min σ, identical max Γ/σ², no tripwire, all finite.

### Ruling implemented

Ryan 2026-08-21: **these become the 018 production FMM defaults.**
`examples/run_dji9443_hover_ct_hpc.slurm.sh` now exports
`FMM_BODY_*` = 17/0.7/109 and `FMM_WAKE_*` = 16/0.6/38 (overridable at
submission; arms submitted before 2026-08-21 ran the production knobs
8/0.4/20 and 4/0.4/50). Only the 018 hover-CT launcher is changed — the driver
defaults, and therefore every other study sharing it, are untouched.

Campaign consequence: a 30-rev production run projects from ~48 h (which timed
out at rev 28.5) to **~11 h** at 71.7 s/step — the NT-ladder turnaround roughly
quarters, and the NT=144 arms that could not reach maturity inside 72 h now can.
Cross-comparability across the switch is carried by the two nulls measured here:
family +0.023%, knobs +0.0038%, against a 0.6–0.85% within-rev p-p.
