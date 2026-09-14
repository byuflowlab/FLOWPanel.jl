# Frozen cold-solve harness

Run from an isolated FLOWPanel campaign worktree with a dedicated Julia environment
whose Manifest points FLOWPanel, FastMultipole and FLOWVPM at clean, annotated-tagged
worktrees. See repository HPC policy and `BYU_ORC_AGENTS.md` on ORC. All Julia work
for this pilot, including tests and precompilation, runs on a Slurm compute node.

## Commands

Set `COLD_PROJECT` to the dedicated environment, `CAMPAIGN_PINS` to the absolute
TOML pins file, and `RUNG=R1`. Pins have `[packages.FLOWPanel]`,
`[packages.FastMultipole]`, `[packages.FLOWVPM]` tables with `path`, `tag`, `sha`.
Use an absolute, new `OUTDIR` and a separate absolute, empty `BENCH_CASE_ROOT` for
**each process** under the consolidated data root. Existing output generations are
never overwritten. The wrapper sets Julia and BLAS counts independently before
startup, including OpenMP controls needed by OpenMP-backed OpenBLAS.

```bash
# Sequential controls; separate processes from timing/precompilation.
bash benchmark/run_cold_process.sh 1 1 test/runtests_benchmark_cold.jl
bash benchmark/run_cold_process.sh 4 1 test/runtests_benchmark_cold.jl

# Fresh output/fixture paths must be set before each command.
STAGE=baseline bash benchmark/run_cold_process.sh 4 1 benchmark/rotor_hover_solver_cold_smoke.jl
# Preserve smoke/selected.toml and use its absolute path as CONFIG_FILE below.
STAGE=verify bash benchmark/run_cold_process.sh 4 1 benchmark/rotor_hover_solver_cold.jl
STAGE=verify bash benchmark/run_cold_process.sh 64 1 benchmark/rotor_hover_solver_cold.jl
STAGE=verify bash benchmark/run_cold_process.sh 64 64 benchmark/rotor_hover_solver_cold.jl
COLD_INVESTIGATION=1 STAGE=verify bash benchmark/run_cold_process.sh 64 1 benchmark/rotor_hover_solver_phase2_profile.jl
```

`CONFIGS=fgs:krylov_ilu` selects both seeds; use one kind for separate profile
processes. Profiling requires `CONFIG_FILE`. The CPU and allocation recordings
occupy separate solve regions and are validated after recording. CPU reports use
native frames and thread/task grouping. The raw allocation profile retains stacks;
its report lists sampled allocations (sample rate 0.01), not exact attribution.

## Optimization-campaign controls (v8)

The initialized-FGS optimization campaign
(`BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_initialized_cpu_optimization_plan_20260911.md`)
adds three environment controls; all default to the pilot behavior when unset.

- `COLD_PREPARED_ONLY=1` skips the fresh scope entirely (constructor-dominated;
  the plan excludes construction and warm starts). Prepared sampling, gates and
  summaries are unchanged.
- `SCREEN_SET="key:v1,v2,...[;key2:...]"` replaces the default `screen` axes
  with an explicit one-factor-at-a-time roster off the rung seed (e.g.
  `inner:1,2,3,5` for the plan's step-1 inner-iteration screen; the seed value
  always runs). Values parse with the seed field's type; unknown keys, wrong
  types, out-of-range values, non-`screen` stages, and combination with
  `CONFIG_FILE` all fail before any filesystem work. Generated screen
  candidates still calibrate their own stopping tolerance through the
  staircase, and a failed screen candidate is recorded in its `status.toml`
  and skipped rather than stopping the roster; baseline/verify/selected
  executions still stop at the first failure.
- `COLD_PROFILE_REPS=N` accumulates the CPU profile over N repeated prepared
  solves (resets outside recorded regions; final solution still validated).
  The allocation profile remains a single solve.

`benchmark/run_cold_opt.slurm.sh` drives the campaign sequence on one node:
controls → smoke → prepared-only seed baselines (j4/b1, j64/b1) → `SCREEN_SET`
screen (j64/b1) → accumulated seed profile. `COLD_OPT_RUNG` (default `R2`),
`COLD_OPT_SCREEN_SET` (default `inner:1,2,3,5`) and
`COLD_OPT_STAGE ∈ {all, controls_smoke, screen_profile}` parameterize it.

## Calibration and gates

Only generated baseline/screen candidates calibrate. The smoke calibrates once and
checks two independent constructors and repeated zero-guess solves per solver.
Every `CONFIG_FILE` is immutable, even with `STAGE=baseline`; selected FGS tolerances
must be positive. Calibration snapshots and final selected settings are preserved.
A failing thread arm never retunes. `screen` defines a future bounded roster; the
pilot does not run it. All current stages stop at the first failed candidate.

Configuration schema accepts all fields produced by `cold_seed` plus optional
Boolean `diagnostic`. Missing/unknown fields, incorrect types, nonfinite numbers,
invalid ranges/enums, inconsistent cache settings, invalid memory/thread controls,
and output/fixture collisions fail before directories or fixture work. The memory
ceiling is positive and at most 500 GiB. Mesh existence is checked before loading.

Each accepted solve requires solver convergence, authoritative BC relative L2
at most `1e-6`, and the memory gate. A certified FMM residual remains authoritative;
if direct is evaluated, its residual must also pass and evaluator disagreement must
be at most `1e-7`. Uncertified FMM may use explicitly labeled direct acceptance on
R1/R2. FMM certification and residual, direct residual, authoritative evaluator and
residual, and disagreement remain separate columns. Warmups, excluded fresh solves,
convergence recordings, CPU profiles and allocation profiles save validation.
Repeated/fresh solutions must agree with the prepared reference within `1e-8`.
Failures write candidate status and propagate a nonzero exit; launchers stop.

## Measurement boundaries and limitations

Prepared timing measures only `_solve!` with a constructed and exercised solver.
Fresh timing measures constructor plus first zero-guess `_solve!`; previous solver
references and constructor timing tuples are released before the next constructor.
Reset, `begin_step_solution!`, frozen direct RHS assembly, accuracy evaluation,
retained-size traversal, and explicit GC are outside measured regions. Reset clears
solution guesses/history while retaining prepared caches. Solver constructors get
independent cache donors. No public solver API changes are introduced.

After excluded warmups, each mode takes 5 samples below 60 seconds, 3 below 600
seconds, otherwise 2. `trials.csv` records every sample including setup and solve;
`summary.csv` records min, median, max, spread and repetition count. Profile
measurements are not timing samples; `unprofiled_trial.csv` is saved separately.

`retained_bytes` is deduplicated object size for body plus solver.
`process_peak_rss_bytes` is the lifetime process RSS peak, including earlier
initialization/diagnostics, and is not retained object size. Both must fit the
ceiling. Allocation bytes and GC seconds are separate metrics. The Slurm limit
provides the shared process memory ceiling.

`estimated_inner_sweeps` and `estimated_fmm_passes` are inferred counts; `-1`
means unavailable (Krylov). They are **not measured component attribution**.
No phase timers, worker/color accounting, hardware counters, cache-size breakdown,
fixed-work instrumentation, solver optimization or complete ILU fairness audit is
claimed. The pilot compares fixed settings; it cannot identify an optimized winner.

Each process saves the actual loaded package paths, SHAs, tags and status, the
Project/Manifest and its hash, thread controls, affinity, fixture/RHS hashes,
configuration hash, Julia/BLAS versions and memory ceiling. `CAMPAIGN_PINS` enforces
clean worktrees and annotated execution tags before fixture loading. The campaign
preparation helper adds a data-symlink commit after its input tag: create and record
an annotated tag for that actual execution commit before running.

The Slurm driver also supports `COLD_PILOT_STAGE=controls_smoke` for a
non-preemptible test allocation, followed by `COLD_PILOT_STAGE=timing_profiles`
with the preserved absolute `CONFIG_FILE`. All three timing settings and both
profile processes still execute sequentially on the same node in the second job.

`benchmark/run_r4_diagnostics.slurm.sh` is the pinned R4 follow-up diagnostic
launcher. It runs the retained lexicographic configuration at Julia thread
counts 1, 4, 8, 16, 32, and 64 with BLAS=1 on explicit physical-core lists.
Each arm alternates two ten-trial uninstrumented and two ten-trial instrumented
batches, checks solution and convergence-history equivalence, writes exclusive
stage timers and an actual nonself GEMV census, samples per-thread CPU ticks,
and retains a thread/task-complete profile. Prepared timing excludes fixture and
solver construction; validation and reset remain outside each timed region.
The default `all` mode runs every stage in one allocation.
# Saved bases for staged screens (v9)

R4 now uses the frozen 58,192-panel `65_209` mesh. Its starting settings
P8/MAC0.4/leaf100/inner3 are provisional R2-derived settings, not an R4
optimum. R4 uses certified evaluation and direct fallback when certification
is inconclusive, as required by the initialized-CPU plan. Existing thresholds
and R1–R3 evaluation behavior are unchanged.

`COLD_OPT_STAGE=attribution` runs controls, smoke, prepared j4/b1 and j64/b1
baselines, then profiling, without a tuning screen. `COLD_MIN_REPS=10` sets a
minimum number of unprofiled repetitions; its default 1 preserves the old
adaptive count. `COLD_OPT_PROFILE_REPS=20` requests 20 accumulated prepared
CPU solves in the launcher (default 10). Timed and profiled solves remain
separate, with resets and validation outside the timed/profiled regions.

`SCREEN_BASE_FILE=/absolute/bases.toml` accepts the selected-config TOML schema
with one or more accepted configurations. It requires `STAGE=screen` and an
explicit `SCREEN_SET`, and conflicts with `CONFIG_FILE`. The roster contains
each base plus its one-factor neighbors, deduplicated by configuration. FGS
tolerances are reset and independently staircase-calibrated for every point.
The input file is copied and hashed in process provenance. For the campaign
launcher, pass `COLD_OPT_SCREEN_BASE_FILE` with `COLD_OPT_SCREEN_SET`; the base
file applies only to the screen process. v8 had no way to screen around saved
winners: its screen always started from the frozen rung seed.
