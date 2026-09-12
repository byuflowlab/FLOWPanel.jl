# Cold investigation implementation checkpoint — 2026-09-10

Stopped at Ryan's request to discuss usage/cost. Implementation is incomplete,
uncommitted, and **not campaign-ready**. No HPC campaign was launched. Do not
interpret the smoke timings as performance evidence or the plan as completed.

## Saved implementation

- `benchmark/common.jl`: positive `BENCH_BLAS_THREADS` override; unset retains defaults.
- `benchmark/phase1_case.jl`: optional `BENCH_CASE_ROOT`; restore frozen velocity
  before deriving BC source strengths in `reset_cold!`.
- `benchmark/fgs_cold_common.jl`: draft shared constructors, explicit seed configs,
  bounded screening roster, calibration, cold reset, accuracy checks, prepared/fresh
  timing, TOML provenance/configuration, CSV trials, independent CPU/allocation profiles.
- `benchmark/rotor_hover_solver_cold.jl`: new entry point.
- `benchmark/rotor_hover_solver_phase2_profile.jl`: opt-in `COLD_INVESTIGATION=1` branch.
- `test/runtests_benchmark_cold.jl`: agent-authored control tests; NOT executed.

Six affected/new Julia files passed recursive syntax parsing on 2026-09-10.
The helper loaded and generated 25 R1 screening configurations before the smoke.
Neither of these checks establishes runtime correctness.

## Smoke result and blocker

Command used Julia 4 / requested BLAS 1, R1, both seed solvers, `MEMORY_GIB=12`,
LineGauss, and fresh temporary outputs:

- `/tmp/fgs-cold-smoke-20260909.log`
- `/tmp/fgs-cold-smoke-20260909-results/`
- `/tmp/fgs-cold-smoke-20260909-case/`

Banner reported BLAS 1; after fixture loading, provenance and the generated
configuration directory reported BLAS 8. Root cause remains unidentified.
The run was interrupted with SIGINT; the execution session returned exit 130.
No benchmark remains intentionally running. The last code edit adds assertions
after fixture loading and before/after each timing trial so this drift fails
instead of silently producing accepted timings. It does not fix the cause.

The FGS calibration CSV contains certified FMM BC relative L2 values
5.90100783e-7 at callback 6 and 9.53567916e-8 at callback 7. These are diagnostic
observations only. The trial/profiling/accuracy protocol has not completed
validation; no comparative speedup or winner is established.

## Known unfinished work (address before another substantive run)

1. Identify BLAS reset during fixture loading, resolve it, then run the control
   tests and a fresh R1 smoke at at most four local threads. Avoid concurrent
   test processes while collecting timings.
2. Validate CONFIG_FILE schema/types/ranges before loading the fixture. Move
   MEMORY_GIB validation before filesystem effects and expensive initialization.
3. Verification/profiling must preserve selected tolerances: current Krylov
   calibration can mutate rtol even with CONFIG_FILE. Disable recalibration for
   selected configurations; validate and record failures instead.
4. Review fresh-constructor memory lifetimes: the excluded fresh trial's `setup`
   tuple can retain the previous solver into the next constructor. Clear it.
5. Review accuracy semantics: direct fallback currently retains FMM's false
   certification field while allowing authoritative direct acceptance. Add an
   explicit evaluator/backend label, save the direct crosscheck result (the
   current excluded warmup discards it), and validate evaluator agreement.
6. Finish actual component/work instrumentation, cache/leaf/interaction memory
   reporting, fixed-work experiments, thread-screen launching, and fairness
   accounting. Krylov FMM work is currently -1 (unavailable), not measured.
   Phase timers, worker/color attribution and hardware counters are not implemented.
7. Add usage documentation (driver currently refers to a README not yet created),
   validated configuration selection, report/plots and campaign launch preparation.
   No solver optimizations have been selected or implemented.

## Provenance / campaign boundary

Runtime package loading was verified locally:

| Package | Checkout HEAD | State |
|---|---|---|
| FLOWPanel | 8755eb24fe8608cf645f5328127cab428b3960b5 | dirty |
| FastMultipole | 8c2500678989d5f9173e4b5f60b0af7f1570efe9 | dirty |
| FLOWVPM | 8b0b70da91297f84cf7a48f33d70c4d6f50373bb | dirty |

All resolve to the live sibling development checkouts under
`/Users/ryan/Dropbox/research/projects/`. Unrelated edits were preserved. Before
an official campaign, select/capture the intended baseline, create annotated
tags and clean worktrees for all three packages, and point the campaign Manifest
at those worktrees. Follow the existing plan and repository HPC policy.

## Usage control for resumption

Both support agents hit the usage limit before finishing their assigned checks.
They inherited the parent model; no lower-cost model override was selected.
Ryan explicitly asked about cheaper subagents. On resumption, use explicit
lower-cost model selection for bounded searches/test execution, pass minimal
context rather than full-history forks, and keep reasoning/edits with the main
agent. Do not spawn agents or start additional investigation until Ryan resumes.
