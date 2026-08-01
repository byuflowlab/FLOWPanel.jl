# Publishable DJI 9443 Hover Convergence Campaign

**Opened:** 2026-07-30

**Current phase:** Phases 0–2 CLOSED; 3/4/5 executing (see status below)

**Item-level approvals:** Technical [ ]; clear-context [ ]; user [ ]

## Current status and next actions (as of 2026-08-01, session handoff)

A fresh agent should read: this file → `ops_reference.md` → `decision_rules.md`
→ the phase file being executed. No repo exploration should be needed.
`ledger.md` holds every harvested number with provenance.

### Where the campaign stands

**Phases 0, 1, 2 are CLOSED.** Headline results:

| quantity | value | source |
| --- | --- | --- |
| **Baseline B0** (η=1.0, σ/c=0.68, N=4, NT=36) | **CT̄ = 0.07170, 95% CI [0.07164, 0.07173]**, revs 15–29, M1 PASS | phase_01 |
| Das ladder {0.205, 0.41, 0.82, 1.64}c | {0.07006, 0.07187, 0.07230, 0.07084} — top rung NON-monotone; **no plateau at this σ** ⇒ **Das\* = 0.41c PROVISIONAL** | phase_02 |
| dt (legacy pinned NT=72) | 0.06852, −3.9% vs NT=36 ⇒ **dt NOT converged; Phase 3 = full ladder incl. NT=144** | phase_00 |
| **N-sensitivity** (ruling 11) | N=1 (d=0.6σ) −0.75% / ε_Γ 2.77% **FAIL**; N=2 (1.2σ) +0.01% / 0.35% **PASS**; ⇒ **adopt N=2 provisionally** | phase_05 |
| σ refinement (L1, σ/c 0.297) | NOT settled in 20 revs — big excursion revs 11–14 (peak 0.0802), flat tail 0.07251 ± 0.00007 ⇒ **provisional +1.13% vs B0** (σ moves CT **UP**, opposite to dt). Extension 13017106 running | phase_04 |
| Campaign scatter | per-rev std ±0.0002–0.0003 (≈5× tighter than legacy `overlap_pps` runs) | ledger |

**Cross-cutting insight (do not lose):** Das, σ, and N are NOT independent
axes — they share the dimensionless clearance group **d/σ** (d = handoff
distance ≈ N·Das). Phase 5 measured the threshold between 0.6σ (fails) and
1.2σ (passes); the σ ladder is therefore itself a clearance sweep at fixed
physical Das, and raising OVERLAP without raising PPS would raise σ and undo
clearance. Any future "converge one axis at a time" plan must account for
this coupling (phase_04 log consequence (d); phase_05 STAGED section).

Cost model (recalibrated from B0's 13.2 h / 77k particles):
**t ≈ 9.0 h/(σ/c)**, **N_p ≈ 52.4k/(σ/c)**.

**Tooling/infrastructure done:** scaffolding + `p018_*` launcher deployed
(md5-verified); T4 FMM discriminator PASSED (Das axis clean of the FMM
|Das|-radius confound); Phase 10 spanwise meshes (60/80 × 185 ct4) shipped;
restart chaining validated end-to-end; **hub-radius mesh variant implemented**
(`scripts/generate_dji9443_mesh.sh --hub-r-over-r X`, ceiling ≈ 0.15R).

### Jobs on the cluster

Check `squeue`/`sacct`; re-open access with Ryan's `! ssh orc -fN` if 2FA
blocks. **No in-session monitors survive a context clear — poll directly.**
Harvest recipe + retention rules: `ops_reference.md`. Respect the 6-job cap.

Six jobs were submitted, verified against `squeue` at 2026-08-01 10:06 MDT.
**As of 2026-08-01, two have finished and four are still RUNNING**, so two
slots are free under the 6-job cap. Elapsed is as of the 10:06 timestamp;
"expect done" assumes ~13 h total for b0-carrier cases, ~16 h for σ=L1 cases,
~26 h for NT=72, all well inside their time limits. Every case's exact knobs
are reproducible from its launcher tag plus the `--export` overrides listed
here.

**Completed, NOT yet harvested — no output has been read or interpreted:**
`13016062` (`p018_hub0p15_b0`, COMPLETED 13:05:13) and `13011375`
(`p018_L1_warm`, COMPLETED 19:34:17). Their "action on completion" entries
below are still outstanding. Note the dependency order: the `L1_warm` gate
cannot be evaluated until the settled cold L1 exists, which waits on
`13017106`.

| job | case (tag + overrides) | phase | elapsed → expect done | action on completion |
| --- | --- | --- | --- | --- |
| 13011373 | `p018_nt72` — NT=72, η=2.0, PPS=2, rlxf 0.15 (launcher-set) | 3 | 17:31 → ~19:00 MDT | harvest; compare to B0 (0.07170) on matched settled windows; if \|Δ\| > 0.5% submit `p018_nt144` as restart-chained 48 h segments |
| 13011375 | `p018_L1_warm` — `RESTART_STEP=719,RESTART_NAME=p018_b0,RESTART_PATH=data/p018_b0,P018_SETTLE_REVS=30` | 4 | **COMPLETED 19:34:17** (not harvested) | warm-pilot gate vs the *settled* cold L1 (needs 13017106 first): ≤0.25% ⇒ L2 may run warm |
| 13016062 | `p018_hub0p15_b0` — `P018_RUN_NAME=p018_hub0p15_b0,RHPC_MESH_FILE=dji9443_20260731_45_185_capped_captess4_hub0p15.msh` | 5 | **COMPLETED 13:05:13** (not harvested) | harvest; hub-vs-stock delta at N=4. NB hub mesh removes blade area inboard of r/R 0.136 ⇒ compare DELTAS and Γ̄(r/R) outboard of 0.2R only, never absolute CT vs experiment |
| 13016670 | `p018_hub0p15_nrows2` — as above + `NWAKEROWS=2` | 5 | 5:34 → ~18:00 MDT | harvest; hub×N cell |
| 13016742 | `p018_hub0p15_nrows1` — as above + `NWAKEROWS=1` | 5 | 4:05 → ~19:00 MDT | **the hub test's payoff**: does the hub mesh rescue N=1 (were near-axis particles the cause of its −0.75% / ε_Γ 2.77% failure)? |
| 13017106 | `p018_L1_s2` — `P018_RUN_NAME=p018_L1_s2,P018_SETTLE_REVS=24,RESTART_STEP=719,RESTART_NAME=p018_L1,RESTART_PATH=data/p018_L1,OVERLAP=2.4,P_PER_STEP=11,MERGE_R_FACTOR=0.0053` | 4 | 0:18 → ~02:00 MDT (Aug 2) | harvest; stitch with `p018_L1` on the `step` column (skip s2's zero placeholder rows) for a ≥15-rev M1 window; **gates the entire σ fit and the warm-pilot comparison** |

All earlier jobs are harvested and recorded in `ledger.md` (12993801–04,
13007210, 13011374, 13011982/83, 12950996; failures 13006768 and 13015838
diagnosed — see phase_01 / phase_05 logs).

**Staged but NOT submitted** (waiting on free slots; full rationale in the
named phase files):
- `p018_L1_ov3` — `OVERLAP=3.0, P_PER_STEP=14, MERGE_R_FACTOR=0.0053` (σ/c
  held at 0.297): tests Ryan's hypothesis that higher overlap shortens settle
  time, with σ deliberately fixed so a settle-time win isn't bought with a
  clearance loss (phase_04 log, consequence (d)). **Run before committing
  L2's segments.**
- Coarse-σ Das×N cross matrix, N ∈ {1,2} × η ∈ {2.0, 4.0} (phase_05 STAGED
  step 2) — collapses all {N, Das} runs onto d/σ to test whether the Phase-2
  non-monotonicity is clearance physics.
- `p018_nrows8` (phase_05, N=8 upper-side sensitivity point).
- Phase 3 `p018_nt144`, Phase 4 `p018_L2`, and all final-settings cases
  (`p018_trunc6/green/nomerge/final`, Phases 10/11 mesh cases) — gated on the
  ladders above.

### Next actions, in priority order

1. **Harvest the six in-flight jobs** as they land (recipe in
   `ops_reference.md`; M1 = cycle-mean over ≥15 settled revs, M2 = ε_Γ from
   `circulation_te`). Append to `ledger.md` + the owning phase file.
2. **Phase 3 decision** from `p018_nt72`: if |Δ| > 0.5% vs B0 (expected),
   submit `p018_nt144` as restart-chained 48 h segments.
3. **Phase 4 σ fit** — gated on `p018_L1_s2` (13017106), NOT on the raw L1
   run: the cold L1 has no valid settled window (see phase_04). Stitch first,
   then fit CT̄(σ)=CT∞+Aσ^p, then L2 (`p018_L2`, warm if the pilot passed).
   **Run `p018_L1_ov3` before committing L2's segments** (overlap-vs-settle-
   time test, phase_04). **Phase 4 also owns two deferred re-tests:** (a) the
   Das top-two rungs at σ=L1 (Phase 2's no-plateau branch), (b) N=1 at L1,
   where d/σ=1.38 clears the bracketed threshold — if N=1 passes there,
   Ryan's minimal-N preference is satisfiable at production σ.
4. **Phase 5 wrap** once the three hub runs land: fill the Das×N×σ
   interaction table (phase_05 STAGED section, step 2's coarse-σ cross matrix
   is still unrun) and settle the clearance criterion (d/σ threshold vs 017's
   (d/σ)^−3.4 prior).
5. Phases 6–8, 10, 11 remain gated on final settings; Phase 9 last.

### Open items needing Ryan

- **`data/pitching_wing_exp/` fixture is MISSING** from disk (untracked, so
  not recoverable from git; absent on cluster/Trash; no session command
  touched it). `test/runtests_unit_pitching_wing_exp.jl` cannot run until it
  is restored via **Dropbox web → Deleted files**. Everything else in the
  suite passes.
- 20 GB ParaView set for the pinned-NT=72 run is at Ryan's
  **`~/p2e_nt72_das2p0_ov6/`** (cluster copy deleted after verified transfer)
  if a visual check of that below-anchor result is wanted.

### Uncommitted source changes this session (all test-gated)

- `src/FLOWPanel_postprocess.jl` — isolated agglomerate (empty stencil) gets a
  zero μ gradient + warning instead of an `ArgumentError` abort; needed for
  the hub mesh's root slivers. Non-isolated failures still throw.
- `~/Dropbox/research/projects/FLOWVPM.jl/src/FLOWVPM_subfilterscale_models.jl`
  — `Estr_direct_multithreaded` guards np=0 (this WAS the documented
  "small-N SFS direct bug"; the T4-era `SFS_OFF` workaround is now unnecessary).
- `test/runtests_unit_postprocess.jl` (new contract),
  `test/runtests_unit_simulate.jl` (missing `cross` import).
- `scripts/generate_dji9443_mesh.sh` — `--hub-r-over-r`.
Both source files are deployed and md5-verified on the cluster (FLOWVPM lives
at `~/projects/FLOWVPM.jl` there). Full suite green except the missing-fixture
test above.

**Standing diagnostic rule (Ryan, re-affirmed 2026-07-31): an OutOfMemoryError
thrown inside `merge_particles!` typically means the simulation blew up** —
the diverged particle cloud wrecks the merge spatial hash; memory is the
symptom, not the cause. Confirm via |CT| > 1 in the log tail (buffered stdout
truncates on SIGKILL — absence of insane CT in the log is NOT evidence of
boundedness); do not resubmit with more memory. The one recorded counterexample
is a run that was genuinely memory-capped: `p2e_sigF_nofilt` (32G request,
MaxRSS 33.4 GB, smooth CF tail) — distinguishable because CT stayed sane and
RSS approached the request.

## Objective and scope

Deliver a publishable convergence study of the DJI 9443 hover rotor (5400 RPM,
$C_{T,\mathrm{exp}} = 0.072$): a high-confidence converged thrust coefficient
and bound-circulation distribution $\Gamma(r/R)$, with **each discretization
axis converged separately** and an explicit error budget. Context: CT currently
spans 0.034–0.110 across plausible discretizations, so no single number is yet
a prediction. This item runs an HPC campaign on existing machinery
(`examples/rotor_hover_pressure_comparison.jl` +
`examples/run_dji9443_hover_ct_hpc.slurm.sh`, `p018_*` case tags) and reuses
recorded Phase-2e/006/014 results as data points.

All operational detail lives in the subdirectory
`018_dji9443_hover_convergence_campaign/`. A clean-context agent should read,
in order: this file → `ops_reference.md` → `decision_rules.md` → the phase file
it is executing. `ledger.md` is the single running results table.

## Standing rulings (binding on every phase)

1. **RPM = 5400 for all convergence jobs** (Ryan 2026-07-30).
2. **Pedrizzetti relaxation FILTER stays OFF in every run** — do not explore it
   in this study; recommending it in the final report is allowed if evidence
   suggests it would help (Ryan 2026-07-30). Stock corrected-Pedrizzetti
   relaxation itself stays ON with stock `rlxf` (Ryan 2026-07-28); relaxation is
   carried as an error-budget term (≈ −0.005 CT, monotone in rlxf, unconverged).
3. **Do not tune to the target**: $C_{T,\mathrm{exp}}$ is compared only after
   the error budget closes; it never selects settings (Ryan 2026-07-28).
4. **No nwakerows ladder** — the legacy-σ ladder was non-monotone and rejected
   (Ryan 2026-07-30, BRAINSTORM/014). Instead: choose `NWAKEROWS` high enough
   that particle-on-body influence isn't stunted by smoothing-radius overlap
   (handoff distance $d \gtrsim 4\sigma$; N=4 here), plus one N=8 spot-check.
5. **Memory: flat 64G, 64 cores per job, ≤6 concurrent** (Ryan 2026-07-30).
6. **Convergence metric: long averaging window; the key requirement is that the
   mean settles** (Ryan 2026-07-30). The old 5-rev gate is retired with a
   limit-cycle defense; see `decision_rules.md`.
7. **BRAINSTORM/016 (surface-vorticity shedding) is PRE-AUTHORIZED** as the
   contingency if particle properties cannot be converged with existing knobs
   (Ryan 2026-07-30). Try existing machinery first.
8. `DAS_REFRESH` stays **false** everywhere (`Backslash` LU caching bug);
   `DAS_KINEMATIC_ARC=false` (tangent) for comparability with the recorded
   ladders; PressureBernoulli(unsteady=false) caveat disclosed in publication.
9. **SFS (DynamicSFS) and viscous (CoreSpreading) particle models stay ON in
   every campaign run** (Ryan 2026-07-31; driver defaults already comply —
   `SFS_OFF` is permitted only in backend-A/B diagnostics like T4, never in a
   ladder rung). **Mesh spanwise and chordwise refinement are convergence axes
   in their own right** (Phases 10/11, respectively at fixed chord and fixed
   span). **The convergence observable is CT̄ AND the spanwise circulation
   distribution Γ̄(r/R), determined from the trailing-edge μ jump**
   (BoundCirculationMonitor `circulation_te`; slice estimator as cross-check)
   — every axis must pass BOTH (Ryan 2026-07-31).
10. **HPC disk retention (Ryan 2026-07-31):** after a run is harvested and
    deemed no longer useful, **delete its VTK files, keeping the CSV
    summaries** (`*_CT_vs_rev.csv`, `*_CT_per_rev.csv`, `*.toml`,
    `monitors/*.csv`). **Keep VTK for at most 3 runs at a time** (prefer the
    actively-running and most decision-relevant). If ParaView inspection of a
    run would help answer a convergence-behavior question, **ask Ryan to
    examine it** rather than hoarding histories.
11. **nwakerows is not a convergence axis (Ryan 2026-07-31, superseding the
    framing of ruling 4):** "convergence of nwakerows doesn't even make sense"
    — the quantities to converge are **particle shedding (σ, spacing) and
    Das**. Direction: try N=1 (and N=2) at campaign settings; treat N as a
    modeling choice whose sensitivity is measured and disclosed, not driven
    to a limit. (Ruling 4's d≳4σ criterion and the N=8 spot-check stand as
    sensitivity probes, not as a convergence claim.)

## Fixed operating point and unit conversions

Mesh `45_185_ct4` (36,752 panels), $R = 0.119$ m, $\rho = 1.179$,
$c(0.75R) \approx 0.128R$; staged startup (SPINUP 1.5 revs @ 0.4 start
fraction, freestream peak 5 m/s, ramp/hold/withdraw 1.0/1.5/4 revs, SETTLE ≥ 12);
velocity formulation unless stated; `PARTICLE_SHEDDING=sigma_overlap`;
truncation 4R (radius 1.5R hard-coded — disclosed, not converged);
`NWAKEROWS=4`. Conversions (freeze factor 0.4 included; tangent Das):

$$\mathrm{Das}/c(0.75R) = 0.41\,\eta\,(36/\mathrm{NT}), \qquad
\sigma/R = \frac{2\pi}{\mathrm{NT}}\cdot\frac{\mathrm{OVERLAP}}{\mathrm{PPS}}, \qquad
\sigma/c(0.75R) = 1.363\,\frac{36}{\mathrm{NT}}\cdot\frac{\mathrm{OVERLAP}}{\mathrm{PPS}}$$

$h = \sigma/\mathrm{OVERLAP}$; classical convergence needs $h$ refined faster
than $\sigma$ ($\sigma^q/h$ const, $q>1$) ⇒ OVERLAP increases down the σ
ladder. `MERGE_R_FACTOR` $= 0.138\,(\sigma/R)$ holds $r_\mathrm{merge}/\sigma$
at the stock ratio. **Pinning σ across NT** (at fixed OVERLAP): halve
`P_PER_STEP` when NT doubles (σ ∝ OVERLAP/(NT·PPS)); the `SigmaOverlap` policy
then sheds the correct per-step count on its own.

## Phase gates

| Phase | Deliverable | Status |
| --- | --- | --- |
| [0](018_dji9443_hover_convergence_campaign/phase_00_harvest.md) | Harvest 3 in-flight cluster jobs; seed ledger; deploy scaffolding | **CLOSED 2026-07-31** — all 3 harvested; nt72match ⇒ Phase 3 full ladder incl. NT=144 |
| [1](018_dji9443_hover_convergence_campaign/phase_01_baseline_settling.md) | Baseline B0 (η=1.0, σ/c=0.68) + settling demo + restart validation + cost recalibration | **CLOSED 2026-07-31** — B0 = 0.07170 CI [0.07164,0.07173] (revs 15–29, M1 PASS); restart chaining validated (seam within natural variability); cost model recalibrated |
| [2](018_dji9443_hover_convergence_campaign/phase_02_das_length.md) | Physical Das ladder {0.205, 0.41, 0.82, 1.64}·c → Das\* = plateau lower edge; FMM/Direct discriminator | **RESOLVED 2026-07-31 (no-plateau branch): Das\* = 0.41c PROVISIONAL** — no pair passes M1+M2; top rung non-monotone (−2.02%); handoff confound flagged ⇒ Phase 5 decision-critical; T4 PASSED |
| [3](018_dji9443_hover_convergence_campaign/phase_03_timestep.md) | NT ladder at pinned physical Das AND σ (36 → 72 [→ 144]) → NT\* | `p018_nt72` RUNNING (13011373); legacy pinned pair says full ladder incl. NT=144 likely |
| [4](018_dji9443_hover_convergence_campaign/phase_04_sigma_ladder.md) | σ ladder 0.68 → 0.297 → 0.151 c with q ≈ 1.2 → production σ + extrapolation term | L1 cold (13011374) + warm pilot (13011375) RUNNING; also re-tests Das top-two at σ=L1 per Phase 2's no-plateau branch |
| [5](018_dji9443_hover_convergence_campaign/phase_05_shedding_distance.md) | N-sensitivity + **staged Das×N×σ interaction study seeking convergent clearance criteria** (ruling 11; Ryan directive 2026-07-31 — see phase file's STAGED section) | nrows1 (13011982) + nrows2 (13011983) RUNNING; then coarse-σ Das×N matrix → σ-refined N=1 → fine-σ Das re-test |
| [6](018_dji9443_hover_convergence_campaign/phase_06_truncation.md) | 6R vs 4R confirmation at final settings | OPEN |
| [7](018_dji9443_hover_convergence_campaign/phase_07_green_delta.md) | GreenReconstruction ΔCT and ΔΓ(r/R) at converged settings | OPEN |
| [8](018_dji9443_hover_convergence_campaign/phase_08_merging.md) | Merging null demonstration + speedup report | OPEN |
| [9](018_dji9443_hover_convergence_campaign/phase_09_final_error_budget.md) | ≥30-settled-rev production run; CT̄ ± CI, Γ̄(r/R) ± CI; closed error budget | OPEN (last) |
| [10](018_dji9443_hover_convergence_campaign/phase_10_mesh_spanwise.md) | Mesh spanwise ladder 45 → 60 [→ 80] at n_airfoil=185, final wake settings (needs local mesh generation) | OPEN — staged 2026-07-31 |
| [11](018_dji9443_hover_convergence_campaign/phase_11_mesh_chordwise.md) | Mesh chordwise ladder 145 → 185 → 249 at span 45, final wake settings (meshes exist) | OPEN — staged 2026-07-31 |

Guessed parameter values are starting points — executing agents may adapt
based on results, recording deviations in the phase file and `ledger.md`.
Dependencies: P0 → all; P2 needs the B0 carrier definition (can run alongside
P1); P3 needs P0+P2; P4 needs P2; P5 needs P1; P6/P7/P8/P10/P11 need P4
(final settings); P9 last, after 10/11 fill their budget terms.
Respect the 6-job cap throughout.

## Contingency chain (Phase 4 failure path, in order)

(a) raise OVERLAP; (b) `SHEDDING_R_OVER_R` root-clip test (existing knob —
cheap proxy for a larger hub: near-hub particles barely move; 40_40 hub
exclusion measured ≤ +2e-4 CT); (c) hub-radius mesh variant — **IMPLEMENTED 2026-07-31**
(`scripts/generate_dji9443_mesh.sh --hub-r-over-r X`; recipe cuts interior
XSecs inside the hub then moves the stock root out; validated ceiling ≈ 0.15R,
larger hubs fold the loft non-watertight; hub0p15 mesh md5-verified on
cluster; 3 comparison runs staged in phase_05);
(d) **016 implementation (pre-authorized)** per its phase_03 architecture doc.

## Error budget (assembled in Phase 9)

Skeleton in `018_dji9443_hover_convergence_campaign/error_budget.md`: Das
model-form (~0.2%/doubling) · dt residual · σ extrapolation · N spot-check
spread · truncation null · Green Δ · merge Δ · relaxation (≈ −0.005 CT,
unconverged) · PressureBernoulli(unsteady=false) disclosure · Dirichlet
tangency contamination (2d App. G, ~0.8–1% of CT at n=185) · radial
truncation 1.5R fixed.

## Decision log

- 2026-07-30 — Item opened; scaffolding written; rulings 1–8 recorded (planning
  session with Ryan). Correction vs the planning doc: pinning σ across NT
  requires **halving** `P_PER_STEP` as NT doubles, not doubling it
  ($\sigma = 2\pi R\,\mathrm{OVERLAP}/(\mathrm{NT}\cdot\mathrm{PPS})$).
- 2026-07-31 — Ruling 9 added (Ryan): mesh spanwise + chordwise refinement
  staged as Phases 10/11; SFS+viscous stay ON in all runs; Γ̄(r/R) from the
  TE μ jump is a co-equal convergence observable with CT̄ on every axis.
  Launcher gained `p018_chord145/chord249/span60/span80`; span meshes need
  local generation (Phase 10 prerequisite).
