# Plan: 006 — R/2 downstream-filtered relaxation vs full-wake relaxation

> Self-contained recipe. **First execution step: copy this file to
> `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/plans/006_rlxf_filter.md`**
> (the user wants it persisted there so it survives a context clear), then proceed.

## Context

BRAINSTORM/006 (`006_rotor_hover_converge_stable_near_reference.md`) wants a rotor-hover
particle wake that is **both stable and near the reference CT** (0.068–0.072; experimental
anchor ≈0.072). Item 005 found a tension at NT=36, SFS+viscous on, rlxf=0.3:
- **Full-wake corrected-Pedrizzetti relaxation (rlxf=0.3)** → stable plateau but CT too low
  (mean ≈0.056–0.062).
- **Relaxation off everywhere** → CT jumps into the reference band (≈0.072) but the wake
  destabilizes (≈2.7× the 2/rev ripple).

Hypothesis under test: applying relaxation **only downstream** of the rotor plane (beyond
R/2) — i.e. **not** relaxing the near-rotor zone where blade circulation is strongest —
keeps the far wake stabilized while preserving near-disk strength, recovering CT magnitude
without losing stability. We measure **how much CT changes** (and whether stability holds)
between full-wake relaxation and the R/2-filtered relaxation.

The spatial-filter machinery already exists in the library
(`RelaxationPlaneFilter`, `plane_filtered_relaxation`, per-step `refresh_relaxation_filter!`
in `src/FLOWPanel_wake.jl`; verified by unit tests). The only code change needed is to wire
an env-gated filter into the rotor-hover example, which currently never passes `relaxation=`.

**Decision (confirmed with user): run two matched runs** (full vs R/2-filtered) on the
**full 005 stable-wake schedule**, 40_40 mesh, NT=36, 2 threads.

## Geometry of the filter (rotor-specific)

- Rotor radius `R = 0.119`; rotor reference frame is **index 1** in `frames`.
- `axial_dimension = 1` (**+x is downstream / below the disk**, the direction the wake
  convects). Spin axis is x, i.e. normal to the R/2 plane → frame tracking is a correct
  no-op but harmless.
- R/2 plane: `point = (R/2, 0, 0)`, `normal = (+1, 0, 0)`, `i_frame = 1`. The predicate
  `dot(X - point, normal) ≥ 0` relaxes particles with `x ≥ R/2` (downstream); the
  near-rotor band `0 ≤ x < R/2` is left **unrelaxed** — exactly "no relaxation closest to
  the rotor plane."

## Code change (one file)

`examples/rotor_hover_pressure_comparison.jl` — add an env-gated relaxation scheme and pass
it to the wake. The example already imports `StaticArrays` (line 6) and defines
`axial_dimension`.

1. Insert after the `particle_relax` / `sfs_off` env block (≈ line 211), before
   `sfs_choice`:
   ```julia
   # Item 006: optional spatially-filtered relaxation. When RELAX_FILTER_DOWNSTREAM_R is
   # set, apply Pedrizzetti relaxation only to particles that have propagated at least
   # RELAX_FILTER_DOWNSTREAM_R*R downstream (+axial) of the rotor plane, leaving the
   # near-rotor band unrelaxed. Unset/NaN => unfiltered full-wake relaxation (default).
   relax_filter_downstream_R = parse(Float64, get(ENV, "RELAX_FILTER_DOWNSTREAM_R", "NaN"))
   base_relaxation = pnl.FLOWVPM.relaxation_correctedpedrizzetti
   relaxation_scheme = if isnan(relax_filter_downstream_R)
       base_relaxation
   else
       d = relax_filter_downstream_R * R
       plane_point  = SVector{3,Float64}(ntuple(i -> i == axial_dimension ? d   : 0.0, 3))
       plane_normal = SVector{3,Float64}(ntuple(i -> i == axial_dimension ? 1.0 : 0.0, 3))
       pnl.plane_filtered_relaxation(base_relaxation, plane_point, plane_normal; i_frame=1)
   end
   ```
2. In the `pnl.PanelParticleWake(rotor; …)` call (≈ line 248), add the kwarg:
   ```julia
       SFS=sfs_choice,
       relaxation=relaxation_scheme,
   ```

This records the filter in run metadata automatically (the wake reads `pfield.relaxation`
into `pfield_optargs`; `FLOWPanel_replay.jl` serializes `RelaxationPlaneFilter`), giving
provenance for each run.

## The two runs (full 005 stable-wake schedule, NT=36, 40_40 mesh, 2 threads)

Shared env (the validated 005 non-damping baseline; SFS+viscous are ON by default):
```
RHPC_MESH=40_40 NT=36 RPM=6000 BERNOULLI_ONLY=true SAVE_VTK=true
SPINUP_REVS=1.5 SPINUP_START_FRACTION=0.4
MAGVINF_START=0.0 MAGVINF_PEAK=5.0 MAGVINF_END=0.0
FREESTREAM_RAMP_REVS=1.0 FREESTREAM_HOLD_REVS=1.5 FREESTREAM_WITHDRAW_REVS=12 SETTLE_REVS=12
TRUNCATION_DEPTH_R=4 NREVS=10
```
(`NREVS=10` is intentionally below the schedule total so `required_revs ≈ ramp+hold+
withdraw+settle ≈ 26.5` drives length ⇒ ≈28 rev incl. spin-up, ≈1000+ steps.)

Threads: `julia --project --threads 2`, with `export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2`.

Run them **sequentially** (each uses 2 threads; avoid oversubscription). Overwrite each
run's own directory per CLAUDE.md (`rm -rf` the target dir first).

**Run A — baseline (full-wake relaxation, unfiltered):**
`RUN_NAME=rotor_hover_relax006_full`, `RELAX_FILTER_DOWNSTREAM_R` unset.

**Run B — R/2 downstream filter:**
`RUN_NAME=rotor_hover_relax006_filt_r05`, `RELAX_FILTER_DOWNSTREAM_R=0.5`.

Driver script to create and launch in the background (writes a log per run):
```bash
cd /Users/ryan/Dropbox/research/projects/FLOWPanel.jl
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
COMMON="RHPC_MESH=40_40 NT=36 RPM=6000 BERNOULLI_ONLY=true SAVE_VTK=true \
SPINUP_REVS=1.5 SPINUP_START_FRACTION=0.4 MAGVINF_START=0.0 MAGVINF_PEAK=5.0 MAGVINF_END=0.0 \
FREESTREAM_RAMP_REVS=1.0 FREESTREAM_HOLD_REVS=1.5 FREESTREAM_WITHDRAW_REVS=12 SETTLE_REVS=12 \
TRUNCATION_DEPTH_R=4 NREVS=10"

rm -rf data/rotor_hover_relax006_full data/rotor_hover_relax006_filt_r05

env $COMMON RUN_NAME=rotor_hover_relax006_full \
  julia --project --threads 2 examples/rotor_hover_pressure_comparison.jl \
  > data/006_full.log 2>&1

env $COMMON RUN_NAME=rotor_hover_relax006_filt_r05 RELAX_FILTER_DOWNSTREAM_R=0.5 \
  julia --project --threads 2 examples/rotor_hover_pressure_comparison.jl \
  > data/006_filt_r05.log 2>&1
```
Launch with Bash `run_in_background: true`; monitor the logs for the per-step progress and
the printed **"Item 005 plateau diagnostics"** block (mean + peak-to-peak CT over the settle
window). Expect a few hours per run.

## Reading the result / comparison

For each run:
- CT history: `data/<RUN_NAME>/<RUN_NAME>_CT_vs_rev.csv`, column **`CT_bernoulli`** (already
  positive). With `BERNOULLI_ONLY=true` the other CT columns are NaN.
- Plateau stats: the run prints the "Item 005 plateau diagnostics" block (mean CT and
  detrended peak-to-peak over the final `SETTLE_REVS` window); also re-derive with
  `julia --project examples/analyze_stable_wake_oscillation.jl data/<RUN_NAME> 12 6000`.

Report: **ΔCT** = mean(filtered) − mean(full) over the settle window, plus each run's
peak-to-peak ripple (stability). Reference points from item 005 at this config: full-wake
NT=36 ≈ **0.0562** mean / ripple ≈0.0104; relax-off ≈ **0.0722** / ripple ≈0.0062×(2.7).
The R/2-filtered run is expected to land between these — the experiment quantifies where.

## Critical files

- `examples/rotor_hover_pressure_comparison.jl` — add `relax_filter_downstream_R` env block
  + `relaxation=relaxation_scheme` kwarg (only change).
- Reused library API (no change): `pnl.plane_filtered_relaxation`,
  `RelaxationPlaneFilter`, `refresh_relaxation_filter!` (`src/FLOWPanel_wake.jl`), already
  called each step from `propagate!`.
- Analysis: `examples/analyze_stable_wake_oscillation.jl`.

## Verification

1. **Smoke / wiring check (fast, before the long runs):** run the filtered case with
   `RHPC_SETUP_ONLY=true` (skips `simulate!`) to confirm the example constructs the
   `PanelParticleWake` with the `RelaxationPlaneFilter` and the metadata records it — and
   `DIAGNOSE_PARTICLE_GAMMA=true` for a few steps (override schedule to a tiny `NREVS=1`,
   `FREESTREAM_*`/`SETTLE_REVS` small) to confirm near-rotor particles keep raw Γ while
   downstream particles realign.
2. **Both full runs complete** and each writes `<RUN_NAME>_CT_vs_rev.csv` plus the plateau
   diagnostics block (non-NaN `CT_bernoulli`).
3. **Comparison:** tabulate mean CT and ripple for full vs filtered over the settle window;
   state ΔCT and whether the filtered run stays stable. Inspect particle VTK in ParaView if
   the CT delta is surprising (confirm the R/2 boundary in the unrelaxed-vs-relaxed Γ field).
