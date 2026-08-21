# 023: Runtime Cost Profiling of Mature-Wake 018 Production Runs

**Opened:** 2026-08-20 (Ryan request, staged and executed same session)
**Current phase:** Phase A (wall_s mining) → Phase B (HPC per-phase timing)
**Item-level approvals:**
- [ ] Technical Completion
- [ ] Clear-Context Approval
- [ ] User Approval After AI Discussion

## Current status and next actions

> **RESET BRIEF (a) — 2026-08-20 (staging).** Item opened because 018
> production runs (`45_185_ct4` = 36,752 panels, NT=36, ~181k particles at
> maturity) cost ~170–230 s/step on 64 cores and wall out one rev short of
> the 30-rev target in 48 h; NT=144 arms reach only rev ~13–15 in 72 h. No
> runtime-cost analysis exists anywhere in 018 (the `src/FLOWPanel_simulate.jl`
> per-step timer comment explicitly notes early-vs-late step cost was never
> measured). Plan: Phase A mines the per-step `wall_s` already recorded by
> `WakeHealthMonitor` CSVs in local `data/p018_*`; Phase B warm-starts
> `p018_cs_f1_l3p4` from its retained mature restart set (steps 1025–1034,
> rev ~28.5) on HPC and re-runs the existing `benchmark/` per-phase timing
> harnesses at production scale; Phase C ranks cost-reduction levers with
> measured evidence. Prior small-scale verdict to re-test (40_40 mesh, 39k
> particles, `benchmark/results/rotor_hover_pressure_steady_aero_timing_summary.md`):
> `sa_wake_influence` 87.9%, `sa_solve` 6.3%, `sa_body_influence` 5.3%.

> **Update (same day):** Phase A COMPLETE (`phase_A_wall_s_mining.md`) —
> per-step cost is ~linear in particle count (intercept ~10–20 s at 36,752
> panels, slope ~50–100 s per 100k particles, R² up to 0.92); the particle
> census, not mesh or dt, is the dominant lever (6R-truncation run: 342k
> particles at 328 s/step); node-to-node variance ~±30% confounds cross-run
> slope comparisons. Phase B harness written, locally smoked, deployed
> (cluster src md5-verified in lockstep), job SUBMITTED.

> **Update 2 (same day):** Phase B main run 13245638 DONE, all gates PASS
> (thrust continuity to 4 digits, reconciliation <0.1 s, 181,307 particles).
> Production split: **wake influence 64.2% / body influence 25.3% / solve
> 9.3% / I+O 1.1%**; profiler shows 49% thread utilization and near-field
> direct kernels dominating ⇒ three sensitivity jobs submitted (leaf 100,
> mac 0.5, 32 threads). Results in `phase_B_hpc_profile.md`.

> **Update 3 (2026-08-20 evening) — FMM tuning stage (Ryan directive) DONE.**
> All sensitivity arms + nf-split + tune/perturb landed; results consolidated
> in `phase_B_hpc_profile.md`. Headlines: (a) ~75% of a production step is
> the Dynamic-SFS estimator (`Estr_fmm!` postcalc) walking the nearfield
> lists — the wake FMM velocity evaluation itself is ~7 s; (b) production
> wake knobs MISS the 1e-4 relative wake-error target (2.1e-4) while the
> tuned+perturbed point (p16/MAC0.6/leaf24) is 3.7x faster AND 12x more
> accurate; body pass floor ~36 s is kerneloffset-radius-inflation-bound,
> not knob-bound; (c) `FastMultipole.tune_fmm_perturb` (new reusable
> function) corrected a tune_fmm leaf misprediction on the body pass
> (64->36 s). Gotcha for any future restored-state harness: run the full
> step head (maneuver!+reset/freestream/kinematic+update_TE!) before ANY
> influence evaluation — the restored first wake row is stale until
> update_TE! and silently NaNs every target (cost ~3 failed jobs).

In-flight jobs:

| job | case | action on landing |
|---|---|---|
| 13246437 | full-step A/B at tuned knobs (wake 16/0.6/24, body 15/0.8/236), 5 steps | steady mean vs 172/177 s baseline + CT continuity → confirms projected ~94 s step |
| ~~13245638/763/764/765/815/869/13246037~~ | Phase B timing / sensitivity / nf-split / tree / tune | **ALL DONE** — `phase_B_hpc_profile.md` + `benchmark/results/` |

## Objective and scope

Measure where wall-clock goes in a production-configuration 018 timestep with
a **mature** (~181k-particle) wake, and produce a ranked, quantified list of
cost-reduction recommendations (projected s/step and hours/run per lever).
Out of scope: **implementing** any recommendation (each lands as a proposal
routed to Ryan; body-solver levers route to item 021, which owns solver cost
and has already measured ILU-GMRES at 19× vs plain GMRES); modifying the
frozen RHPC driver (`examples/rotor_hover_pressure_comparison.jl`); touching
any live 018 run or its retained restart sets (read-only use only).

## Standing rulings (binding on every phase)

1. Inherited from 018: RHPC driver is frozen — env-var overrides and the
   `benchmark/` drivers are the only sanctioned edit surface.
2. Inherited HPC ops: Julia **1.11.7** (never `module load julia`),
   single-node threaded, `THREADS` = `--ntasks` exported through all five
   `*_NUM_THREADS` vars, banner verified within minutes of submission,
   md5-verify deployed sources against local, disk ≤ 200 G total.
3. Local work ≤ 4 threads; anything > ~20 laptop-minutes goes to HPC.
4. Profiling arms suppress VTK output (timing results go to
   `benchmark/results/`, not `data/`); the retained restart sets are never
   moved, modified, or swept by this item.

## Reference cost datapoints (as of opening)

| config | measured cost | source |
|---|---|---|
| NT=36 production carrier, mature wake | ~170–230 s/step, 64 cores | `data/p018_cs_f1_l3p4/.../monitor04_wake_health` `wall_s` |
| NT=72 rung | ~82 s/step | `018.../phase_03_timestep.md` |
| 022 coarse mesh `56_57` | ~87 s/step | 022 RESET BRIEF (c) |
| f1_l3p4 particle census | 108k @ step 198 → 192k @ 598 → ~181k plateau | monitor04 |
| 40_40 mesh / 39k particles phase split | wake_influence 87.9% / solve 6.3% / body 5.3% | `benchmark/results/rotor_hover_pressure_steady_aero_timing_summary.md` |

## Phase gates

| Phase | Deliverable | Status | File |
|---|---|---|---|
| A | wall_s mining of landed runs: cost-vs-step, cost-vs-n_particles fit, NT-ladder comparison | open | `023_018_runtime_cost_profiling/phase_A_wall_s_mining.md` |
| B | HPC per-phase timing + `Profile`/`Profile.Allocs` at production scale from the `p018_cs_f1_l3p4` mature restart | staged | `023_018_runtime_cost_profiling/phase_B_hpc_profile.md` |
| C | Ranked recommendation table (lever, measured % of step, projected savings, risk, routing) | gated on A+B | this file, dated section |

Phase B acceptance: per-phase sums reconcile against the monitor `wall_s` for
the same steps to ±10%; restart sanity = first post-restart CT within noise of
the parent run's final CT (guards the shedding/warm-start gotchas).

## Contingency chain

- `p018_cs_f1_l3p4` restart set incomplete on cluster → fall back to
  `p018_nt144_rs1` step 1881 (578 M, kept unswept; NT=144, so scale
  conclusions carefully).
- Warm-start drifts or first-step CT off → discard step 1 as warm-up and time
  steps 2+ only; if still off, diagnose before trusting any timing.
- Cluster inaccessible (no ControlMaster socket) → Phase A + prior-art
  analysis stand alone; Phase B jobs staged but unsubmitted.
- Timing harness incompatible with production restart config → thin new
  wrapper `benchmark/p018_mature_wake_timing.jl` over `benchmark/common.jl`
  rather than deep edits to the 2026-06 harnesses.

## Deferred

- Implementing any lever (FMM knob retune, plan/tree reuse, ILU-GMRES swap,
  particle-count policy changes) — proposal-only here.
- Thread-scaling study beyond a 2-step spot check (16/32/64).
- Profiling the NT=72/144 arms (only if the NT=36 split doesn't explain the
  ladder costs via the Phase A fit).
