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
