# Investigate FGS scaling in item 021

Prepared and revised 2026-09-09. Status: investigation plan; no new profiling or campaign runs performed. Ryan's revised objective: make FGS as fast as practical with reasonable effort, initially using cold-started solves only, and ensure ILU-GMRES has no obvious comparative handicaps.

## Summary

Prioritize cold-started time to a common physical accuracy target on frozen rotor problems. Defer warm-start benchmarks, RHS-sequence replay, and end-to-end unsteady performance campaigns. Use existing fixtures and profiling machinery to keep setup inexpensive.

Here, **cold start means a zero solution guess**, not an uncompiled Julia process. Compile both paths before timing. Report (1) solve time using an already constructed solver and (2) construction/cache/preconditioner setup plus one cold solve. Repeated timing trials must reset the guess, residual state, and convergence history; they must never inherit the preceding solution. Separate eager and lazy setup costs explicitly so neither solver gets free setup in the first-solution comparison.

The current Phase 3 data supports the concern, but primarily for cold solves:

| Mean solve time | R1 | R2 | R3 |
|---|---:|---:|---:|
| FGS, cold | 1.46 s | 3.31 s | 22.23 s |
| ILU-GMRES, cold | 2.61 s | 6.47 s | 12.67 s |

These are historical wake-developed averages of cold-started solves, excluding the first three steps. They motivate the investigation but are not the new benchmark protocol. These measurements concern mesh scaling at fixed 64 Julia/64 BLAS threads; they do not establish thread scaling. Source: `figures/R1.csv`, `figures/R2.csv`, and `figures/R3.csv` beside this plan, using last-wins deduplication, excluding checkpoint rows and steps at or below `skip_steps=3` (70 retained rows per arm).

Thread utilization is a credible hypothesis: FGS's default nearfield sweep and residual evaluation are serial. However, convergence growth, inherited tuning, memory traffic, and BLAS overhead also need attribution. Leaf-LU caching already exists.

## 1. Establish an accurate, reproducible comparison

- Use the existing `benchmark/phase1_case.jl` frozen rotor fixture at R1–R3, with identical geometry, operator, RHS, precision, and zero initial guess for both solvers. Begin with current code captured in a reproducible baseline; do not require rerunning the historical unsteady campaign.
- For official runs, pin the selected baseline in clean campaign worktrees with annotated tags for FLOWPanel, FastMultipole, and FLOWVPM. Preserve LineGauss and record loaded package paths, Julia version, BLAS implementation, hardware, and effective thread counts. Keep before/after results within this generation; preserve unrelated development edits.
- Validate the frozen benchmark's common physical BC evaluator and its `1e-6` acceptance target before declaring a performance winner. Calibrate each solver's internal stopping tolerance to this target; do not assume identical numeric internal tolerances mean identical accuracy. Cross-check R1/R2 using direct evaluation outside the timed region, and require certified evaluation for larger cases.
- Record historical unsteady BC-certification failures as an unresolved limitation of the motivating data. They do not block the frozen benchmark if its own evaluator passes independent checks. Diagnose unsteady RHS/timing/wake issues only if they reproduce in the frozen case; otherwise defer that separate investigation.
- If the frozen fixture does not reproduce the R3 slowdown, report that explicitly. Only then add one captured R3 operator/RHS as a cold-started secondary fixture if an existing snapshot makes this inexpensive; do not launch a trajectory-replay campaign to obtain it.

### ILU-GMRES fairness audit

- Give each solver its best measured Julia/BLAS thread combination within the same CPU allocation and memory ceiling. Also show fixed-resource, fixed-knob comparisons for attribution; fairness does not require identical algorithm-specific knobs.
- Verify right-preconditioned GMRES, convergence status, restart/memory settings, and iteration limits. Check that restart is not causing obvious stagnation or excess iterations; screen the existing restart setting against twice that value within the memory ceiling.
- Verify ILU factor reuse, FMM plan reuse, and available nearfield-cache behavior. Enable applicable existing optimizations for each solver; include their memory and construction costs. Compare under a common total memory ceiling, counting FGS's mandatory dense caches as well as ILU factors and optional caches. Label uncached runs as diagnostic when they handicap a solver.
- Check ILU pattern/drop settings and the FMM apply parameters for obvious poor choices. Use a bounded local search around existing settings, including neighboring leaf/acceptance settings and drop/fill controls where exposed. Tune to actual cold-solve time at the common accuracy target, not only matvec time.
- Apply identical timing boundaries, compilation warmup, repetition policy, diagnostic exclusion, and numerical validation. Report extra tuning effort spent on FGS and disclose any unresolved ILU limitation rather than claiming a fully optimized global optimum.

## 2. Profile where the time goes

Extend `benchmark/rotor_hover_solver_phase2_profile.jl`, which currently profiles one warmed-up cold solve and hides native frames.

- Start with frozen R2 and R3; use R1 for inexpensive correctness checks. Profile both solvers on the same inputs, initially at their existing settings.
- Separate constructor/cache/preconditioner setup, formulation and RHS work, iterative solver time, and BC measurement. Report cold solve time and setup plus one cold solve; exclude independent diagnostic BC measurement from both performance metrics and report its cost separately.
- Within FGS, time farfield stages, leaf solves, nonself matrix products, scatter, and residual checks. Within ILU-GMRES, separate FMM applications, ILU application, and Krylov operations.
- Count outer updates, actual inner sweeps, FMM passes, and formulation subsolves. FGS iteration counts are not interchangeable with GMRES iterations; R3 also changes the configured inner-sweep count.
- Capture CPU profiles grouped by thread/task with native frames enabled, plus separate allocation profiles, actual allocated bytes, GC time, and peak memory. Use [Julia's profiling facilities](https://docs.julialang.org/en/v1.10.2/stdlib/Profile/); sampled allocation bytes must not be labeled total allocation.
- Record leaf-size distributions, interaction counts, cached matrix bytes, per-worker work, and—for colored sweeps—color widths and barrier/scatter time. Collect Linux CPU utilization and hardware counters where available; retain a timers-and-profiles fallback.
- Measure performance without profiling enabled. Report the campaign's warmed minimum-of-k timing alongside median and spread.

## 3. Test threading and tuning independently

`benchmark/common.jl:55` explicitly sets BLAS threads equal to Julia threads, overriding environment-only attempts to separate them.

- Add an explicit benchmark BLAS-thread override, preserving the existing default and recording the effective setting. No public solver API change is initially required.
- On exclusive Zen3 hardware, screen Julia threads `{1,4,16,32,64}` with BLAS fixed at 1. Compare the historical `64/64` setting and test BLAS `{4,16}` at Julia 64. Record affinity and socket placement; distinguish physical cores from SMT.
- Hold solver knobs and captured inputs fixed during attribution. Measure both fixed-work cost and time to the accuracy target, so faster sweeps cannot conceal worse convergence.
- Then tune standalone FGS specifically for cold-solve time. The old harness calls `stage3_winner()`, but that helper actually restricts candidates to the standalone `stage2_selected(1e-6)` settings before selecting the preconditioner's sweep count. Do not claim the standalone settings were optimized solely for FGMRES. The concrete concern is inherited Gaussian-era selections and dependence on historical selection files while Krylov apply settings were retuned for LineGauss.
- Tune expansion order, acceptance, leaf size, and inner iterations using the existing search machinery, with zero guesses on every trial. Start with neighboring values around the existing settings; expand only when the best admissible result lies on the search boundary. Apply the ILU-GMRES fairness audit before publishing the comparison.
- Expand to R4 only after R2/R3 identifies a reproducible cause and promising change. Keep local work at four threads or fewer.

## 4. Optimize according to measured evidence

| Measured bottleneck | Candidate change |
|---|---|
| BLAS overhead or contention | Use fewer BLAS threads; choose the best measured Julia/BLAS combination. |
| Serial residual work | Parallelize independent leaf residuals using private scratch and a deterministic reduction. |
| Serial nearfield sweeps | Evaluate the existing colored sweep option; measure time to convergence and available parallelism. |
| Scatter or memory traffic | Improve locality and buffering; investigate deterministic target-owned accumulation while preserving required dependencies. |
| Farfield work or worker imbalance | Optimize the dominant FMM stage and partition work by measured interaction cost. |
| Increasing outer iterations | Retune expansion/acceptance and inner-versus-outer work for standalone FGS, measuring cold-started time to the common target. |

Do not simply thread the lexicographic Gauss–Seidel loop: updates depend on preceding leaves. Colored sweeps already exist but have historical divergence evidence, so they require real-rotor convergence testing before adoption.

### Effort limit and order

- First complete the benchmark/fairness fixes, thread screen, and bounded parameter tuning. These existing controls take priority over new algorithms.
- Then implement at most two localized code optimizations, selected by measured cost and expected benefit. Candidates include scratch reuse, removing redundant work, threaded residual evaluation, or a small scheduling/locality improvement. Assess the existing colored option before considering new sweep scheduling.
- Stop expanding this initial effort when remaining opportunities require a solver redesign, a new parallel GS algorithm, extensive scatter restructuring, GPU work, or another solver family. Record those as follow-up proposals with estimated benefit and complexity. FGS-preconditioned FGMRES and warm starts are deferred.
- Do not promise a particular speedup or that FGS will beat ILU-GMRES. Success is the fastest validated FGS found within this bounded effort, a credible explanation of the remaining bottleneck, and a comparison without known easy-to-remove ILU handicaps.

## Validation and deliverables

- Run `test/runtests_unit_solver.jl` and `test/runtests_unit_fgs_history.jl` for solver changes; add `test/runtests_unit_fmm.jl` and the dependency's relevant tests for FMM changes, and a short rotor simulation for integration coverage. Verify repeatability at fixed thread settings and numerical agreement across settings.
- Keep changes only when accuracy gates pass and unprofiled cold-solve timing improves reproducibly beyond measured variability. Prefer simple gains even below 10%; do not impose an arbitrary speedup threshold on low-effort improvements. Investigate regressions above 5% at R1/R2 and retain per-size settings when warranted rather than forcing one global configuration. Report setup and memory tradeoffs explicitly.
- Deliver a ranked bottleneck report, phase/thread scaling plots with backing CSVs, reproducible profiling commands, and before/after cold-start results for both solvers. Show solve time, setup plus first solve, achieved BC error, actual work counts, memory, and chosen knobs. Include the completed ILU fairness audit and deferred opportunities. Quantify the maximum plausible benefit from removing each dominant serial fraction.
- Preserve historical results and unrelated development edits. Follow repository campaign policies for future runs; writing this plan does not launch a campaign.

## Fresh-agent implementation handoff

### First actions and code map

Read the policy files listed below, then inspect `git status --short` in each loaded repository. This plan is the task brief; reading the entire item ledger or regenerating Phase 3 figures is unnecessary. Source line numbers below are navigation hints; use the named functions if they move.

| Entry point | What to reuse or change |
|---|---|
| `benchmark/common.jl:35`, `assert_and_banner` | Add optional `BENCH_BLAS_THREADS`; unset preserves existing behavior, set must be a positive integer. Assert and record actual BLAS threads. |
| `benchmark/phase1_case.jl` | Frozen fixture, `LADDER`, `rotor`, `b`, `rms_b`, `reset_cold!`. Include after `common.jl`. Add optional `BENCH_CASE_ROOT` to isolate its results/knob directories; retain old paths when unset. |
| `benchmark/common.jl:372`, `bc_error!` | Shared physical accuracy check; returns `rel_l2`, `rel_max`, `error_success`, `t_eval`, `epsilon_requested`. |
| `benchmark/rotor_hover_solver_phase2_profile.jl` | Extend with an opt-in cold-investigation mode that uses the new shared constructors and output root; preserve its legacy invocation. Keep CPU and allocation profiles separate from timing trials. |
| `benchmark/rotor_hover_solver_phase2.jl` | Reference for separated setup accounting, ILU stats, cached variants, and memory accounting. Do not include this top-level driver from the new driver: it executes the old roster and depends on historical tuning files. |
| `benchmark/phase1_knobs.jl` | Reuse `adaptive_min_of_k` logic and `margin_tol`; avoid automatic `stage3_winner()`/CSV dependency in the new path. Adaptive timing uses 5 repetitions below 60 seconds, 3 below 600 seconds, otherwise 2. Save all samples to compute median/spread. |
| `benchmark/rotor_hover_solver_phase1_fgstune.jl` | Reference for FGS staircase callback alignment and tolerance calibration. Callback iteration 1 is the cold iterate; handle the final iterate if the iteration cap is exhausted. |
| `src/FLOWPanel_solver.jl` | `KrylovSolver` near 1000, `FGSSolver` near 1500, FGS `_solve!` near 1765, `ILUPreconditioner` near 2209. |
| Loaded FastMultipole `src/solve.jl` | `gs_sweep!`, `solve_leaf!`, `compute_nonself_products!`, `scatter_nonself_influence!`, `residual!`, and the outer solve loop are the targeted profiling/optimization sites. |

Create `benchmark/fgs_cold_common.jl` for shared explicit constructors/reset/validation and `benchmark/rotor_hover_solver_cold.jl` for the new timing/screening driver. Share that helper with the profiling driver's opt-in mode so timed and profiled configurations cannot drift. Avoid adding public solver API unless the selected measured optimization requires it.

### Fixture and cold-reset contract

- R1/R2/R3 have 8,016/15,760/28,752 panels. Mesh names in `LADDER` are `dji9443_20260813_23_73_capped_captess4.msh`, `dji9443_20260813_33_105_capped_captess4.msh`, and `dji9443_20260813_45_145_capped_captess4.msh`, under `examples/data/`. The helper also includes `examples/dji9443_trailing_edge.jl`. Check these assets exist in the campaign worktree before running.
- Fixture is Dirichlet with source strength in column 1, solved strength in column 2; RPM 6000, radius 0.119 m. Retain its constructed-body shedding procedure and panel core size. Do not reconstruct the fixture independently.
- `reset_cold!()` restores apparent velocity, BC sources, zero solved strengths, and the frozen source potential. RHS is `b = -potential_frozen`, with `rms_b = norm(b)/sqrt(n)`. Direct source assembly is outside solver timing and reported as common fixture preparation cost.
- Explicitly disable Krylov `warmstart` and FGS `project_solution`, with FGS `solution_history_length=0`. Reset any convergence recorder per trial. Verify the reused-solver trial agrees with a fresh-solver trial and that solver counters describe only the current call.
- Before BC evaluation, copy solved column 2 and restore `rotor.velocity .= frozen_velocity`; the evaluator derives sources from that velocity. Accept only `error_success && rel_l2 <= 1e-6`. Start with evaluator `safety=0.1`, cap 20, MAC 0.5, leaf 20, independently of solver knobs. If uncertified, increase evaluation accuracy or use direct evaluation on small cases; never turn a certification failure into a pass. Record `rel_max` diagnostically.
- `CACHE_B=1` currently trusts `bcache_<rung>.bin` without provenance. Initially use `CACHE_B=0` and share the once-assembled RHS within each process. Do not load an old cache merely because its filename matches. Do not set `SKIP_B=1` for any solve.
- `phase1_case.jl` prints obsolete `TUNED` fallback values when tune files are absent. The new driver must use its explicit configuration, not those fallbacks. Route helper outputs to the new `BENCH_CASE_ROOT` so old CSVs cannot silently select settings.

### Explicit starting configurations and fair cache lifetimes

These are **seeds to validate and retune**, not certified winners for the new baseline:

| Rung | FGS `(P, MAC, leaf, inner)` | Krylov FMM `(P, MAC, leaf)` |
|---|---|---|
| R1 | `(6, 0.3, 50, 10)` | `(17, 0.65, 6)` |
| R2 | `(8, 0.4, 100, 10)` | `(17, 0.65, 6)` |
| R3 | `(6, 0.3, 100, 5)` | `(16, 0.65, 6)` |

- FGS: `max_iterations=300`, `rlx=1.0`, `shrink=true`, `recenter=false`, `reverse_pass=false`, `cache_leaf_lu=true`, `sweep_order=:lexicographic`, `verbose=false`. Calibrate `tolerance` from a new BC staircase and `margin_tol`, rather than reusing historical absolute thresholds.
- ILU-GMRES: `method=:gmres`, `itmax=500`, `atol=1e-14`, initial `rtol=1e-6`, `memory=50`; validate and adjust tolerance against the common metric. Screen memory 50 versus 100. Start ILU with `leaf_size=10`, `multipole_acceptance=1.0`, `max_pattern_entries=8192*rotor.ncells`, `equilibrate=false`, `diagonal_shift=0.0`. The current implementation exposes pattern controls, equilibration, and shift, **not a numerical drop tolerance**; do not invent one. Screen ILU leaf `{5,10,20}` and MAC `{0.8,1.0}` before broadening.
- Krylov defaults `cache_tree=false` and `cache_nearfield=false`; leaving them untouched can be a real handicap. Screen plan reuse with `cache_tree=true`, then the cached variant under the common memory ceiling. `persistent_plan=true` allows geometry-dependent plan reuse across independent zero-guess trials; this is permitted setup reuse, not a solution warm start. FGS likewise retains its tree, matrices, and LU.
- Report fresh-constructor plus first-solve timings directly, including lazy build work. Separately report prepared-state cold solves with permitted persistent plans/caches populated. Clear cache donors for every fresh-setup trial so previous candidates cannot donate uncharged work. Record reuse flags and build times, and never mix these two timing modes in one speedup.
- For the first official R1–R3 screen, use the existing 500 GiB exclusive Zen3 allocation as a shared ceiling, accounting for body, solver, factors, plans, and caches without double-counting shared objects. Do not run the full 16/128/500 GiB campaign now. Local smoke cases must fit available memory and at most four threads.

### Driver controls, artifacts, and execution order

Implement these controls for the new driver: `RUNG`, `CONFIGS=fgs:krylov_ilu`, `EXPECT_JULIA_THREADS`, `THREADING_MODE`, `BENCH_BLAS_THREADS`, `BENCH_CASE_ROOT`, `OUTDIR`, and `STAGE=baseline|screen|verify`. Baseline validates/times the seeds; screen performs the thread/parameter experiments specified above (one Julia process per thread setting); verify reruns selected candidates. Store selected full configurations in TOML and let verification/profiling load that same file through `CONFIG_FILE`. The new profiling mode is selected by `COLD_INVESTIGATION=1` and uses the same controls/configuration file. These interfaces are planned additions, not commands available before implementation.

1. Build the shared helper and baseline driver, isolate outputs, and validate R1 frozen BC against direct evaluation. Confirm zero guesses, compilation exclusion, counters, and fresh/prepared timing boundaries.
2. Record R2/R3 baseline timing and profiles for both solvers. Perform thread screening and the ILU fairness checks before selecting FGS code changes.
3. Run bounded FGS knob tuning, then at most two evidence-backed localized optimizations. Keep separate baseline, tuning, and optimized configurations/results. A failed or capped candidate is a recorded result, never silently omitted.
4. Validate selected settings and changes with the test commands below, then rerun unprofiled R1–R3 comparisons. Add R4 only under the existing evidence gate.
5. Write `cold_scaling_report.md` beside this plan with the measured explanation, fair winner comparison, remaining limitations, and stopping rationale. Update the item ledger with concise results; notebook writing follows its separate approval policy.

Outputs go under an explicit unique directory in the consolidated data root, with a per-rung/per-configuration/per-Julia-BLAS subdirectory and one writer per directory. Save `provenance.toml`, `config.toml`, individual `trials.csv`, summary/component CSVs, convergence sidecars, CPU profile data plus thread-grouped text, and allocation reports. Record fixture hash, all three package tags/full SHAs, actual threads, memory/caching settings, accuracy/status, and timing mode. Never append this new schema to historical `phase2.csv` or reuse resume keys that omit configuration/thread/provenance identity.

Example **after implementing the new controls**, from the selected worktree, with both output paths set explicitly to a fresh local temporary or campaign data directory:

```bash
RUNG=R1 CONFIGS=fgs:krylov_ilu EXPECT_JULIA_THREADS=4 \
THREADING_MODE=multi BENCH_BLAS_THREADS=1 CACHE_B=0 \
BENCH_CASE_ROOT=/absolute/new-case-dir OUTDIR=/absolute/new-results-dir \
STAGE=baseline FLOWPANEL_FILAMENT_REG=linegauss \
julia --project=. --startup-file=no -t 4 benchmark/rotor_hover_solver_cold.jl
```

Local tests, each with no more than four threads:

```bash
julia --project=. --startup-file=no -t 4 -e 'include("test/runtests_unit_solver.jl")'
julia --project=. --startup-file=no -t 4 -e 'include("test/runtests_unit_fgs_history.jl")'
julia --project=. --startup-file=no -t 4 -e 'include("test/runtests_unit_fmm.jl")'
```

Only run the FMM suite when that path changes; add FastMultipole's relevant suite and a short rotor integration check as specified above. For HPC, adapt `benchmark/slurm/p2_table.sh` only as a launcher reference: its historical partition/QOS and precompile behavior are not a current availability result. Follow `HPC.md`, refresh availability, pin Julia deliberately, serialize precompilation before fan-out, and set explicit allocation/thread/affinity controls. No remote job status was checked for this plan.

## Handoff references and cautions

- Read `/Users/ryan/.claude/CLAUDE.md`, repository `CLAUDE.md`, and the task-specific policies before execution. Applicable policies include `WORKFLOW.md`, `TESTING.md`, and `HPC.md`; read `MONITORS.md` before diagnosing BC/monitor recovery behavior.
- Current campaign context: `phase_23_context_reset_prompt.md` and `decision_rules.md` beside this plan. Phase 3 records tag `021-lg-gen2d` at FLOWPanel `00390f7`, actual HEAD `e818c63` including the symlink commit, FastMultipole `0ce3ba6`, and FLOWVPM `a627dd9`. Verify full identities against raw provenance before preparing worktrees. Later handoffs describing R2/R3 FGS as pending are superseded by the locally available CSV rows, not by a fresh remote-status check.
- `benchmark/rotor_hover_solver_unsteady.jl:379` measures `solve_formulation!`, not only the iterative kernel; the subsequent BC measurement is excluded from `t_solve`. Use `t_step_net = t_step_total - t_bcerr` for production timestep comparisons.
- FastMultipole `src/solve.jl`: colored sweep around line 1047, serial default around 1074, residual evaluation around 1530. Colored scatter remains serial. Existing leaf-LU caching is default-on. Resolve the actual loaded dependency path rather than assuming a sibling checkout.
- `phase_02_single_step_benchmarks.md` documents existing LU caching and colored-sweep experiments. `results/fgs_determinism/summary.md` documents the previously fixed M2L ownership race and historical threading measurements; these are leads, not current performance evidence.
- Both FLOWPanel and FastMultipole had unrelated uncommitted changes during this review. Do not attribute current source behavior or future timings to historical campaign pins without checking the differences. Old `figures/solver_scaling/fits.csv` contains development-era isolated-solve fits and must not be used as the current Phase 3 scaling result.
