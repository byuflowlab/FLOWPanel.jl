# Plan (2026-07-07): rotor hover cold run with downstream-filtered relaxation — stable CT plateau near 0.072

Stand-alone: a clear-context agent can execute this on the workstation without prior
conversation. All paths relative to the repo root. **Thread budget: 36 — every
simulation runs with `julia -t 36` AND BLAS threads 36 (set
`BLAS_NUM_THREADS`/`OMP_NUM_THREADS`/`JULIA_NUM_THREADS=36`; the example calls
`BLAS.set_num_threads` at startup and prints "BLAS threads: N" to the log —
verify 36/36 there).**

## 0. OPERATING RULES — token budget (read first, non-negotiable)

The user has a 5-hour usage limit this agent MUST NOT exhaust. Simulations cost
wall-clock, not tokens. Rules (carried over from `plans/20260702_convergence.md`):

1. Launch runs **in the background**; wait with **long intervals (20–30 min+)**
   between checks; never poll frequently.
2. Each check is a **minimal probe**: `tail -5` of the run log or a parse of the
   run's `*_CT_vs_rev.csv`. **NEVER read a full run log** (hundreds of steps ×
   per-step solver output). Parse CSVs, grep logs.
3. Substantive analysis/writing only at run boundaries. Between boundaries, idle.
4. Keep a rough running self-estimate of spend; target well under half the budget
   per 5-hour window. When uncertain, idle longer.
5. Avoid token spend during peak hours **7am–1pm MST** (background runs are fine);
   **no jobs running at all during the 24 h starting 12am MST on Sundays.**
6. Token-heavy proposals to the user must include an estimated fraction of the
   5-hour limit (project CLAUDE.md rule).

## 1. Context and goal

Rotor-hover CT under-predicts experiment (~0.072). The 2026-07-03/04 knob sweep
(results table `data/rotor_hover_sweeps.md`; analysis `data/ct_knob_sweep_notes.md`;
summary CSV `data/ct_knob_sweep_summary.csv`; dated section in
`BRAINSTORM/006_rotor_hover_converge_stable_near_reference.md`) established:

- **Pedrizzetti relaxation (stock rlxf=0.3) suppresses hover CT by ≈5.5e-3.**
  Control (pure continuation of the 40_40/NT=36/RPM=6000 baseline) plateaus at
  CT_bernoulli **0.06403 ± 2.8e-4**; with `RELAX_FILTER_DOWNSTREAM_R=0.5`
  (Pedrizzetti applied only to particles ≥0.5R downstream of the rotor plane) the
  settled value is **CT ≈ 0.0695** (18-rev extension run
  `ct_knob_relaxfilter_0p5R_ext`, per-rev means 0.06839→0.06928→0.06953→0.06948
  over revs 14–18).
- **Equilibration is slow at low relaxation: ~5–8 hover revs.** Any 4-rev read of
  a low-relaxation state is a censored lower bound (the sweep's rlxf ladder and
  filter-depth points all under-read for this reason; the ladder's λ→0
  extrapolation 0.0666 was an underestimate). Always verify a plateau before
  quoting a level.
- The rlxf halving ladder was monotone and **stable down to rlxf=0.0375**;
  filter depths 0.25R/0.5R/1.0R read 0.06625/0.06727/0.06757 at 4.3 revs
  (censored). Strain restoration (BODY/PANEL_WAKE_HESSIAN_TO_PARTICLES, alone and
  combined), KERNELOFFSET_TARGETS (5e-4, 2.5e-4), and MERGE_PARTICLES=false are
  **null** — exonerated as CT suppressors.
- Judge CT on **CT_bernoulli**. The PressureLaplace monitors are valid in cold
  runs but corrupted in warm-start continuations (∂u/∂t history warm-up).

**Goal:** a full cold run with `RELAX_FILTER_DOWNSTREAM_R=0.5`, ≥16 revs, that is
(a) stable, (b) CT plateaus to a consistent value over ≥2 revolutions, (c) as near
0.072 as achievable. Expected landing ≈0.069–0.070 (confirming the warm-start
prediction); §5 gives sanctioned follow-ups to chase the remaining ~2–3e-3.

## 2. Repo state prerequisites

- `src/FLOWPanel_warmstart.jl` must contain the **Das-rotation fix** (§2.5 of
  `simulate_warmstart!`: comment "ALWAYS reconstruct by replaying the kinematics"
  + unconditional replay loop). UNCOMMITTED as of 2026-07-07 — verify present, do
  NOT revert. Without it, every warm-start restart from a manifest-complete run
  silently misplaces the first wake row (Das TE shed vectors rotate with the body
  via `rotate_Das!` but are not persisted) and CT diverges from the first
  continued step. Details: BRAINSTORM/006 2026-07-07 section.
- `julia --project=. -e 'using Pkg; Pkg.instantiate()'` if fresh clone.
- `examples/rotor_hover_pressure_comparison.jl` is the driver example; all knobs
  are ENV vars read at include time. `examples/rotor_hover_ct_knob_sweep.jl` holds
  the prior 14-scenario warm-start table (reusable pattern, not needed here).
- Disk: the run writes ~600 steps of VTK (comparable to
  `data/rotor_hover_pressure_comparison`, ~a few GB). Overwrite the run dir on
  reruns (`rm -rf` first); do not create timestamped copies. Keep VTK output ON.

## 3. The primary run

```bash
cd <repo-root>
rm -rf data/rotor_hover_relaxfilter0p5_cold
RUN_NAME=rotor_hover_relaxfilter0p5_cold RELAX_FILTER_DOWNSTREAM_R=0.5 \
RHPC_MESH=40_40 NT=36 RPM=6000 SETTLE_REVS=8 \
BLAS_NUM_THREADS=36 OMP_NUM_THREADS=36 JULIA_NUM_THREADS=36 \
nohup julia --project=. -t 36 examples/rotor_hover_pressure_comparison.jl \
  > data/rotor_hover_relaxfilter0p5_cold.log 2>&1 &
```

- `SETTLE_REVS=8` ⇒ schedule 2 (ramp) + 3 (hold) + 4 (withdraw) + 8 (hover)
  = **17 revs = 612 steps** (dt = 1/3600 s). Freestream fully withdrawn by rev 9
  ⇒ 8 hover revs: covers the ~5–8-rev filtered-relaxation transient plus a 2-rev
  plateau window. All other knobs stay at example defaults (stock rlxf=0.3 applies
  beyond 0.5R; `TRUNCATION_DEPTH_R=4`; merging on).
- Walltime ≈ 6–9 h at 36 threads (~20 s/step early, ~50 s/step in developed hover).
- Startup check (one probe after ~15 min): log shows "BLAS threads: 36 (Julia
  threads: 36)", the config line "Total run length: 17.0 revs (612 steps)", and
  stepping has begun. Then idle with long intervals.
- On completion the example writes
  `data/rotor_hover_relaxfilter0p5_cold/rotor_hover_relaxfilter0p5_cold_CT_vs_rev.csv`
  (columns: step, revolution, CT_bernoulli, CT_laplace_matderiv, CT_laplace_lamb,
  CT_kj) — the analysis input.

## 4. Gates / acceptance

1. **Stability**: run completes 612 steps, no NaN/blowup; per-rev CT_bernoulli std
   ≤ ~4e-4 over the final revs.
2. **Plateau (primary goal)**: last two per-rev means agree within ≤5e-4 AND no
   monotone trend >~2e-4/rev across the last 3 revs. If still trending at rev 17:
   **extend by warm-starting the run from its own last step** — same ENV plus
   `RESTART_STEP=611 RESTART_NAME=rotor_hover_relaxfilter0p5_cold
   RESTART_PATH=data/rotor_hover_relaxfilter0p5_cold
   RUN_NAME=rotor_hover_relaxfilter0p5_cold_ext SETTLE_REVS=11` (⇒ +3 revs; ~2.5 h).
   Repeat/lengthen as needed. Do NOT cold-rerun to extend.
3. **Level**: report the settled CT_bernoulli against the 0.068–0.072 band
   (item 006 acceptance). Expected ≈0.069–0.070. Also report the plateau
   CT_laplace values (valid in a cold run) as a cross-check — disagreement
   >~1.5e-3 between Bernoulli and Laplace at the plateau is worth flagging.

## 5. Flexibility — sanctioned follow-ups toward CT ≈ 0.072

The user authorizes iterating without asking; primary goal is a stable, 2-rev-
consistent plateau as near 0.072 as possible. Each addition gets one line of
justification in the results notes. Warm-start perturbations restart from the cold
run's settled last step (pattern in §4.2; construction-compatible = same
RHPC_MESH/NT/RPM; construction-time knobs like RELAX_* take effect via ENV).

- **Filter depth, settled**: RELAX_FILTER_DOWNSTREAM_R ∈ {1.0, 2.0} as warm-start
  continuations (+4–6 revs each so they settle). The censored sweep showed
  1.0R ≈ 0.5R + 0.3e-3; a settled depth curve tests whether removing more
  near-wake relaxation buys the remaining gap. If a depth destabilizes, record the
  failure rev and fall back to the deepest stable depth — stability is data.
- **Filter + weaker far-field relaxation**: RELAX_FILTER_DOWNSTREAM_R=0.5 with
  RELAX_RLXF ∈ {0.15, 0.075} (ladder was stable to 0.0375 unfiltered).
- **Deeper truncation**: TRUNCATION_DEPTH_R ∈ {6, 8} on top of the 0.5R filter
  (null at stock relaxation, but wake retention may matter once the near wake is
  unrelaxed). Watch np/walltime.
- **Longer settle** to confirm any promising new plateau (2-rev consistency is the
  acceptance bar for every quoted number).
- Combinations that a single result makes compelling (state reasoning in notes).

**Requires user approval first**: src edits; mesh/NT refinement (item 007 scope);
deleting/overwriting `data/rotor_hover_pressure_comparison/`; anything estimated
to consume a nontrivial token fraction (state the estimate).

## 6. Reporting

- Extend the run table in `data/rotor_hover_sweeps.md` with each run (config,
  window, per-rev means, settled CT, stability).
- Append a dated results section to
  `BRAINSTORM/006_rotor_hover_converge_stable_near_reference.md`: settled
  CT_bernoulli mean ± std, plateau window, Laplace cross-check, gap to 0.072, and
  which follow-up (if any) closed it further.
- Deliverable sentence: "stable at CT = X ± Y over revs A–B; gap to experiment
  Z" — plus a recommendation for the next lever if the gap remains.
