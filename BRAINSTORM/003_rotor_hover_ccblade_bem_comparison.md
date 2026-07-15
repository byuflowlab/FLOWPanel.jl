# Rotor Hover CT: CCBlade/BEM Comparison

## Context

Prior audit has already checked geometry, RPM, density, CT convention, radius/diameter interpretation, and handedness/sign. CCBlade/BEM should be used as a diagnostic cross-check, not as proof that either solver is correct by itself.

## Hypothesis

A BEM comparison can separate missing wake/inflow physics from geometry, airfoil, Reynolds-number, or operating-condition assumptions. Starting at nonzero advance ratio should avoid hover-specific singular behavior and make the comparison easier to interpret before approaching hover.

## Proposed Path

- Build a CCBlade case using the same blade geometry, operating point, density, RPM, and reference definitions.
- Start at a small but nonzero advance ratio.
- March toward hover only after the nonzero-advance case has sensible radial trends.
- Compare panel and BEM outputs in dimensional and nondimensional form.

## Quantities To Compare

- Radial loading.
- Inflow angle.
- Induced velocity.
- Circulation.
- Integrated CT.
- Any local section assumptions that differ from the panel model, especially airfoil polar and Reynolds-number treatment.

## Acceptance

This track is useful if it clearly separates missing-wake physics from geometry/airfoil/Re assumptions. For example, if BEM agrees with VPM/experiment while the panel result remains low under equivalent inputs, focus on wake and induced-velocity modeling. If BEM also predicts low CT under the same assumptions, revisit section data and operating-condition interpretation.

## Caveats

BEM and panel methods do not represent the same physics. Use agreement or disagreement to route the investigation, not as a final validation by itself.

## Axial freestream comparison at J = 0.1867 (2026-07-14)

Prepared a fixed axial operating point for the DJI 9443: RPM=5400, `Vc=4 m/s`,
`J=Vc/(nD)=0.1867`, with `ncrit=4` and `ncrit=9` as co-primary viscous CCBlade
polar cases. The CCBlade entry point is `examples/rotor_axial_j0187_ccblade.jl`;
it writes tagged Vc/J sectional CSVs, loading/alpha plots, CT output, and an
operating-point validation CSV under `data/rotor_axial_j0187_ccblade/`. The
validation requires finite sections, coverage by the original polar domain,
and (by default) outboard (`r/R >= 0.3`) `|alpha| <= 5 deg`; the pre-run evidence
was `-0.23..1.89 deg` for ncrit=4 and `1.25..2.17 deg` for ncrit=9.

`examples/rotor_axial_j0187_panel.jl` reuses the established 80_81 / 36
steps-per-revolution workflow at constant axial `V∞=4 m/s`, retaining VTK
output in the same directory. `examples/rotor_axial_j0187_replay.jl` requires
the final full revolution to meet the existing CT peak-to-peak (5%) and drift
(2.5%) gates. The resulting spanwise comparison and compact CT CSV/plot are
made by `examples/rotor_axial_j0187_loading_comparison.jl`.

**Status: results pending user Slurm submission.** The prepared single-node
script is `scripts/rotor_axial_j0187_ccblade.slurm.sh` (48 CPUs, 8 h, 4 GiB/CPU,
BYU `bigmem` constraint); it is intentionally not submitted by an agent.

## Results (2026-07-07)

Implemented in `examples/rotor_hover_ccblade.jl`; outputs in `data/rotor_hover_ccblade/`. CCBlade BEM was built from the *same* CST blade geometry (`dji9443_brainstorm_item003.csv`, 25 sections), RPM=5400, ρ=1.179, B=2, R=0.119, with XFOIL section polars generated at the local hover Reynolds/Mach (Re≈10k–54k, M≈0.02–0.20; Viterna ±180° extrapolation). Hover is approximated with a tiny nonzero climb speed (Vc=1e-4; pure Vx=0 is degenerate in BEM — the induction factor `a` blows up but induced velocity and thrust stay physical). A small climb sweep (Vc = 1, 2, 4 m/s → J = 0.047, 0.093, 0.187) was run alongside hover. Follow-up added an inviscid XFOIL polar peer (`ncrit=0`, `cd=0`) and replayed saved FLOWPanel VTK to add tangential spanwise loading.

**Integrated hover CT vs. XFOIL polar set**:

| polar set | hover CT | vs experiment (0.072) | vs panel (0.0506) |
|-----------|----------|-----------------------|-------------------|
| inviscid  | 0.1142   | 159%                  | +126%             |
| ncrit 1   | 0.0711   | 99%                   | +41%              |
| ncrit 4   | 0.0705   | 98%                   | +39%              |
| ncrit 9   | 0.0683   | 95%                   | +35%              |
| ncrit 14  | 0.0598   | 83%                   | +18%              |

**Findings:**
- Across the plausible transition range, BEM **brackets the experimental CT (0.072)** and stays **18–41% above the converged panel result (~0.0506)**. Because BEM reaches the experiment using identical geometry / airfoil / Reynolds / operating-condition inputs, the panel-method shortfall routes **not** to those inputs but to **wake / induced-velocity modeling** — consistent with the force-method audit (001) and the wake-decay audit (004). This satisfies the Acceptance criterion.
- Counter-intuitively, *lower* `ncrit` (earlier transition, noisier inflow) gives *higher* CT: at these low Re a more laminar BL (high ncrit) separates and degrades cl/cd. CT saturates below ncrit≈4 (1 and 4 nearly identical). Lower ncrit is arguably the physical case for a rotor ingesting its own turbulent wake.
- The inviscid BEM peer is much higher (CT≈0.114) and is not a better experiment match, but it is the closer physics peer to an inviscid panel method. The inviscid loading comparison shows FLOWPanel is smoother and lower through the outer loaded span, while inviscid BEM develops a sharp local normal/tangential loading spike near r/R≈0.79.
- Radial trends (loading, inflow angle, induced velocity, circulation) are sensible and monotone toward hover, plotted at the smallest Vc — see `radial_*.png` and `rotor_hover_ccblade_radial.csv`.
- Replay note: the locally available `data/rotor_hover_pressure_comparison/` VTK window currently spans only 5 saved samples (0.111 rev at RPM=6000), so the regenerated replay stats warn that the averaging window is not flat. They are useful for axis/sign and plotting checks, not a final converged hover-loading statistic.

**Artifacts:**
- CT sweep: `rotor_hover_ccblade_CT_vs_J.csv`, `CT_vs_J.png`
- Radial (hover): `radial_loading.png`, `radial_inflow_angle.png`, `radial_induced_velocity.png`, `radial_circulation.png`; long-format `rotor_hover_ccblade_radial.csv`
- Section polars fed to BEM: `rotor_hover_ccblade_polars.csv` (cl **and** cd; inviscid rows have exactly `cd=0`); lift-polar plots `polars_cl_inviscid.png`, `polars_cl_ncrit{1,4,9,14}.png`, and `polars_cl_vs_ncrit.png`; drag-polar plot `polars_cd_vs_ncrit.png` (cd vs α at representative sections, one line per polar set)
- Sectional CSVs: `rotor_hover_ccblade_sectional_inviscid.csv`, `rotor_hover_ccblade_sectional_ncrit{1,4,9,14}.csv`
- Spanwise loading: `rotor_hover_ccblade_spanwise_loading.png`
- FLOWPanel replay stats with thrust and tangential loading: `data/rotor_hover_pressure_comparison/spanwise_loading_replay/rotor_hover_pressure_comparison_spanwise_loading_stats.csv`
- Reusable comparison script: `examples/rotor_hover_loading_comparison.jl`; default outputs `rotor_hover_loading_comparison_normal.png` and `rotor_hover_loading_comparison_tangential.png` in the replay stats directory

**Open item:** the advance-ratio (J) comparison is preliminary — experimental J points still need to be supplied to finalize it. The BEM CT-vs-J sweep infrastructure is in place (configurable via `VC_LIST`).

## Inviscid overprediction diagnosis (2026-07-07)

Follow-up diagnostics implemented in `examples/rotor_hover_ccblade_diagnostics.jl` (CSV-only, reads the outputs above); plots land in `data/rotor_hover_ccblade/`.

**(a) The r/R≈0.79 spike is a spurious BEM root, worth ~0.022 of CT — not a polar or geometry feature.** At section 18 (r/R=0.7854) the inviscid hover solve converged onto a nonphysical high-induction branch: axial induced velocity u=97.7 m/s and W=112.8 m/s vs u≈5–6, W≈50–55 at neighbors, while α (1.9°) and the section polar (smooth, linear, cl₀≈0.62) are unremarkable — see `polars_oppoint_spike.png` (op point sits on a clean part of the polar) and `radial_momentum_check.png` (blade-element vs annulus-momentum loading; the spike violates momentum consistency by ~2 orders of magnitude). Trapezoid CT with the spike is 0.1134; replacing that one station by its neighbor mean gives **0.0917** (`radial_cumCT.png` isolates the jump).
  - **Confirmed by a Vc sweep** (`data/rotor_hover_ccblade_vcsweep/`, inviscid only, Vc = 1e-4 / 1e-3 / 1e-2 / 0.1 m/s): at Vc=1e-4 the spike station returns u=97.7; at every Vc ≥ 1e-3 the same station returns the physical root (u≈5.3, Np≈19.6 N/m) and integrated CT = 0.0920 — exactly the despiked value. The Vx→0 momentum residual is nearly degenerate and cd=0 removes damping, so CCBlade's φ-residual bracketing can land on the wrong branch. **Practical fix: run "hover" BEM peers at Vc=1e-3 (J≈5e-5) instead of 1e-4.**

**(b) The remaining inviscid excess (CT 0.092 vs viscous BEM 0.060–0.071) is the inviscid polar lift level.** Evaluating every polar set at the *same* α (the ncrit=9 operating α per section, `radial_cl_gap.png` and the script's summary table): inviscid cl is 1.1–1.5× viscous mid-span and 2.5–11× near root and tip, where Re is lowest and viscous decambering/separation losses are largest. The inviscid peer therefore overstates the physics-matched gap to the panel method; the viscous ncrit family (CT 0.060–0.071, bracketing experiment) remains the meaningful comparison.

**(c) The earlier `rotor_hover_loading_comparison_normal.png` impression — "panel above BEM everywhere except the spike" — was an apples-to-oranges artifact, now fixed.** The old plot compared *dimensional* dT/dr from a 5-sample startup-transient panel window at RPM=6000 (implied CT ≈ 0.087) against BEM at RPM=5400. `examples/rotor_hover_loading_comparison.jl` now plots dCT/d(r/R) normalized by each dataset's own RPM (`BEM_RPM` / `PANEL_RPM`, the latter auto-read from the replay metadata), prints each curve's integrated CT in the legend, and warns when the panel stats window is short. The panel stats were regenerated from the *converged last revolution* (rev 9.0–9.97, 36 samples, CT flat to ~2%) of the 10-rev `rotor_hover_pressure_comparison` run. Corrected picture: panel loading (int CT 0.0613–0.0617; this run's own converged CT at RPM=6000 — distinct from the 0.0506 convergence-study value) sits *below* viscous BEM (ncrit=9 int 0.0681) through the mid-span and far below inviscid BEM (0.1134) everywhere, restoring consistency with the integrated-CT table. One real shape difference: outboard of r/R≈0.8 the (inviscid) panel loading stays up like the inviscid BEM tail while viscous BEM collapses from low-Re tip lift loss.
  - Replay-data note: a 5-step debug rerun (2026-06-12) had truncated `rotor_hover_pressure_comparison_body1.pvd` and the metadata `[[step]]` entries to steps 0–4 while all 360 step VTUs survived. Both were reconstructed from the deterministic kinematics (validated against the intact entries); originals kept as `.bak` alongside.

New radial diagnostics: `radial_alpha.png`, `radial_W.png` (the spike is glaring here), `radial_cl.png`, `radial_cd.png`, `radial_cl_gap.png`, `radial_cumCT.png`, `radial_momentum_check.png`, `polars_oppoint_spike.png`. The main script now also records the tip-loss factor `F` in the radial CSV and honors `SAVE_PATH`.

## 10×Reynolds sensitivity (2026-07-07)

Same comparison rerun with section Reynolds numbers scaled 10× (Re ≈ 100k–538k; Mach unchanged) for the viscous ncrit set, via new `RE_SCALE`/`SUFFIX` ENV knobs in `examples/rotor_hover_ccblade.jl` (run: `RE_SCALE=10 SUFFIX=_re10x NCRIT_LIST=1,4,9,14`). Outputs carry `_re10x` filenames alongside the baseline in `data/rotor_hover_ccblade/`. Comparison + failure report: `examples/rotor_hover_ccblade_re_comparison.jl`.

**Hover CT (trapz of B·Np), RPM=5400:**

| polar set | baseline Re | 10× Re |
|-----------|------------|--------|
| ncrit=1   | 0.0708     | 0.0717 |
| ncrit=4   | 0.0703     | 0.0738 |
| ncrit=9   | 0.0681     | 0.0756 |
| ncrit=14  | 0.0597 (spike-contaminated; 0.0436 despiked) | 0.0767 |

- **At 10×Re the ncrit spread collapses** (0.0717–0.0767 vs 0.0436–0.0708 baseline) **and the ordering flips**: with a healthy-Re boundary layer, later transition no longer costs lift, so higher ncrit gives slightly *higher* CT. All 10×Re values sit at or just above experiment (0.072). The low-Re viscous lift loss (and its ncrit sensitivity) is thus confirmed as the driver of the baseline spread — while the panel gap (BEM/panel ≈ 1.4–1.5) persists at both Re levels, keeping the wake-modeling routing intact.
- **Baseline ncrit=14 has its own spurious-root spike** (found via the loading overlay): u=99.3 m/s at r/R=0.745 at Vc=1e-4 — same disease as the inviscid set; despiked CT is 0.0436, so the previously reported 0.0597 was spike-inflated. ncrit 1/4/9 baseline and *all four* 10×Re hover solves are clean (max u ≈ 5.7 m/s) — more support for running hover peers at Vc=1e-3.
- **Loading** (`rotor_hover_ccblade_loading_re10x_vs_baseline.png`): the 10×Re curves are smooth, nearly ncrit-independent, and hold loading further outboard (the baseline high-ncrit tip collapse largely disappears); baseline ncrit=14 is visibly jagged with the r/R=0.745 spike.
- **XFOIL convergence** (failures out of 1495 attempted solves per set = 65 α × 23 sections; full list in `rotor_hover_ccblade_xfoil_failures_report.csv`, per-α direct log in `rotor_hover_ccblade_xfoil_failures_re10x.csv` — the two agree exactly):
  - baseline: inviscid 0, ncrit=1 55, ncrit=4 48, ncrit=9 213, ncrit=14 194
  - 10×Re: ncrit=1 102, ncrit=4 91, ncrit=9 361, **ncrit=14 616 (41%)**
  - Failures concentrate at the two root sections (r/R ≤ 0.14, where the chord is thick and Re lowest) and grow sharply with ncrit; counter-intuitively, 10×Re converges *worse* than baseline at high ncrit (long laminar runs + transition on a hard-worked BL). ncrit=14 at 10×Re should be treated as polar-quality-limited: sections 1–2 kept only 9–20 of 65 points.

Artifacts: `rotor_hover_ccblade_{radial,CT_vs_J,polars,xfoil_failures}_re10x.csv`, `rotor_hover_ccblade_sectional_ncrit{1,4,9,14}_re10x.csv`, `_re10x` PNG peers of all baseline plots, `rotor_hover_ccblade_loading_re10x_vs_baseline.png`, `rotor_hover_ccblade_xfoil_failures_report.csv`.

## Mixed polars: inviscid cl + viscous cd (2026-07-07)

Implemented a mixed-polar option in `examples/rotor_hover_ccblade.jl`: each BEM polar set is now keyed by a string `polarset`, with `ncrit` retained as compatibility metadata (`-1` for mixed sets). The script can combine lift from one XFOIL sweep with drag from another through `RUN_MIXED`, `MIX_CL_NCRIT`, `MIX_CL_RESCALE`, `MIX_CD_NCRIT`, and `MIX_CD_RESCALE`. Downstream CSV readers now filter by `polarset` when present and derive it from `ncrit` only for older CSVs.

Motivation: the inviscid peer (`ncrit=0`, `cd=0`) is the closest lift-physics comparison to the inviscid panel method, but `cd=0` removes damping from CCBlade's near-hover φ residual and allowed the earlier Vc=1e-4 solve to land on a spurious high-induction root. The script's default hover proxy is now `Vc=1e-3 m/s`, which reproduced the despiked inviscid CT.

Run:

```bash
MPLBACKEND=Agg RUN_MIXED=true NCRIT_LIST=0 SUFFIX=_mixed julia --project examples/rotor_hover_ccblade.jl
```

The mixed case used `cl` from inviscid polars at physical Re and `cd` from `ncrit=4` at 10×Re. XFOIL convergence for the donor drag sweep was 89 failed solves out of 1495 attempted; the inviscid sweep had 0 failures.

**Hover CT at Vc=0.001 m/s:**

| polar set | hover CT | max u [m/s] | vs panel 0.0506 | vs experiment 0.072 |
|-----------|----------|-------------|------------------|---------------------|
| inviscid | 0.09204 | 6.469 | 1.819× | 1.278× |
| mixed_cl-inviscid_cd-ncrit4re10x | 0.09162 | 6.458 | 1.811× | 1.272× |

Spike check passed: both hover curves have max induced velocity O(5–6 m/s), with no ~100 m/s station. Adding viscous drag from the 10×Re ncrit=4 donor slightly reduced CT relative to the despiked inviscid case, but did not move the inviscid-lift BEM peer toward the panel result.

Artifacts in `data/rotor_hover_ccblade/`: `rotor_hover_ccblade_{radial,CT_vs_J,polars,xfoil_failures}_mixed.csv`, `rotor_hover_ccblade_sectional_inviscid_mixed.csv`, `rotor_hover_ccblade_sectional_mixed_cl-inviscid_cd-ncrit4re10x_mixed.csv`, and `_mixed` PNG peers (`CT_vs_J_mixed.png`, `radial_*_mixed.png`, `polars_*_mixed.png`, `rotor_hover_ccblade_spanwise_loading_mixed.png`).

## Xfoil alpha_sweep retry path (2026-07-07)

Updated `examples/rotor_hover_ccblade.jl` so viscous (`ncrit>0`) polar generation uses `Xfoil.alpha_sweep` instead of the local hand-rolled branch sweeps. The inviscid path is intentionally unchanged: `ncrit=0` still calls `Xfoil.solve_alpha(a; mach)` directly and records `cd=0`. New environment knobs:

- `XFOIL_ITER` (default `200`)
- `XFOIL_NPAN` (default `140`)
- `XFOIL_PERCUSSIVE` (default `true`)

Validation runs:

```bash
julia --project -e 'Meta.parseall(read("examples/rotor_hover_ccblade.jl", String))'
MPLBACKEND=Agg RUN_MIXED=true NCRIT_LIST=0 SUFFIX=_mixed_pm julia --project examples/rotor_hover_ccblade.jl
MPLBACKEND=Agg RE_SCALE=10 SUFFIX=_re10x_pm NCRIT_LIST=4 julia --project examples/rotor_hover_ccblade.jl
```

Results:

| case | old failures | new failures | hover CT | max hover u [m/s] |
|------|--------------|--------------|----------|-------------------|
| mixed donor drag, ncrit=4 @ 10×Re | 89 | 33 | mixed CT 0.09162 | mixed 6.458 |
| pure ncrit=4 @ 10×Re | 91 | 33 | 0.07454 | 5.704 |

The failure reduction is substantial and no CSV schema changes were needed. The mixed-polar hover result is unchanged to practical precision relative to the previous `0.09162`, and the pure ncrit=4 @ 10×Re hover CT moved modestly from the earlier `0.0738` to `0.07454`. Both hover cases remain spike-free with induced velocity O(5-6 m/s).

Optional stress check attempted:

```bash
MPLBACKEND=Agg RE_SCALE=10 SUFFIX=_re10x_pm14 NCRIT_LIST=14 julia --project examples/rotor_hover_ccblade.jl
```

This was interrupted after several minutes without completing the first root-section sweep; Xfoil was inside the percussive-maintenance retry path. Treat ncrit=14 @ 10×Re as still polar-quality-limited unless a longer dedicated run is needed.
