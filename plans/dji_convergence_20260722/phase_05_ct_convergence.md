# Phase 5 — CT Convergence

## Purpose and required context

Run the adaptive steady and free-wake sequence needed to assess numerical
convergence of `CT` and circulation. Submit only cases justified by completed
evidence.

Read the repository instructions, the dashboard, the top snapshot in the
[Phase 5 log](../../logs/dji_convergence_20260722/phase_05_ct_convergence.md),
and:

- `agent_policies/WORKFLOW.md`
- `agent_policies/TESTING.md`
- `agent_policies/MONITORS.md`
- `agent_policies/HPC.md`

Do not read earlier phase documents. Their accepted geometry, driver,
airfoil-control decision, and Green classification must be summarized in the
Phase 5 log before this phase is approved.

## Preflights and first cold batch

Run steady preflights for all four new mesh/formulation combinations. Compare
normal quad-growing and triangular mu gradients first on the 40-series pair;
extend to 57 series only if the difference reaches 1%. Do not include
quad-with-growth-disabled, which prior histories found pathological and
mesh-dependent.

Initial cold settings:

- `NT=36`;
- wake core `1e-3` m;
- target offset `1e-3` m;
- default quad-growing gradient;
- 5400 RPM;
- production velocity formulation.

Launch no more than three jobs:

1. new 40-series capped/Dirichlet;
2. new 40-series uncapped/Neumann;
3. old 40-series capped control if Phase 1 found an absolute re-fit
   circulation change of at least 1%; otherwise a new 57-series case.

Submit remaining required baselines only as slots open. Query Slurm before
every submission and enforce the three-active-job limit.

## Plateau handling

Use a 2-revolution freestream ramp, 3-revolution hold, 4-revolution withdrawal,
and initially 2 hover revolutions.

A plateau requires:

- final two complete revolution means within 0.5%;
- final-revolution relative standard deviation at most 2%;
- finite body, wake, force, and circulation state.

Warm-start a failed plateau by two revolutions at a time, up to 17 total
revolutions, then label it temporally unresolved.

Mesh, formulation, RPM, and timestep changes require cold runs. Warm starts are
only for construction-compatible cases. Record complete restart provenance.

## Matched effects and perturbations

Measure:

- old/new re-fit effect only if Phase 1 scheduled the HPC control;
- mesh effect within each formulation;
- combined cap/formulation effect at each resolution.

Never describe capped/Dirichlet versus uncapped/Neumann as
formulation-only.

On retained 40-series checkpoints, independently test:

- wake core `1e-3 -> 5e-4` m with target offset fixed;
- target offset `1e-3 -> 5e-4` m with wake core fixed;
- quad-growing versus triangular mu gradient;
- validated Green reconstruction versus matched velocity continuation on the
  capped mesh.

Continue wake-changing perturbations for at least three revolutions so newly
shed particles replace the near wake. Confirm a sub-1% change with one
additional halving. When unstable, bisect between stable and unstable values.

Green rules:

- use it only on capped Dirichlet cases;
- follow the Phase 4 classification and accepted solver route;
- confirm every material result with a matched velocity continuation from the
  same checkpoint, comparing circulation as well as `CT`;
- add a cold confirmation if the warm-started difference remains at least 1%.

## Timestep and final confirmations

Run cold `NT=48` cases for both formulations, scaling per-step relaxation to
preserve its physical rate. Stop timestep refinement if both effects are at
most 1%; otherwise run `NT=72` only for affected formulations.

Finish with cold 57-series capped/Dirichlet and uncapped/Neumann cases using
the finest stable timestep, accepted independent regularizations, and accepted
mu-gradient method. If the last refinement remains above 1%, report that
dimension as unresolved.

## Case records and retention

For every case, log:

- exact settings, formulation, mesh, and restart provenance;
- job state, elapsed time, peak memory, and output/error logs;
- final block `CT` mean/variability and circulation metrics;
- completion/failure status, root cause, and reuse safety;
- whether full VTK is retained and the decision caused by the result.

Preserve no more than four full study histories. Never delete an output needed
for a warm start, an incomplete case lacking a safe summary, or anything
outside this study.

## Exit gate

Produce a compact convergence summary identifying resolved and unresolved
dimensions and the final candidates for Phase 6. Review methods, evidence, and
conclusions; update the log and dashboard; then stop for explicit Phase 6
approval.

