# Cold harness validation and HPC pilot — execution plan 2026-09-10

User-authorized scope: harness validation plus an R1 HPC pilot. Preserve dirty worktrees and historical results. All Julia execution is on HPC; another agent uses the local machine. No public solver API changes.

## Review and harness fixes
Installed OpenBLAS uses OpenMP (parallel64 == 2); OMP_NUM_THREADS=1 before startup preserved one BLAS thread through GEMMs and Julia-thread work. Do not replace the runtime unless further checks fail.
- Validate schema, types, finite numbers, enums, constructor ranges, memory, thread settings and output/fixture paths before mkdir or fixture load. Accept seed fields plus diagnostic only.
- CONFIG_FILE is immutable selected settings in timing and profiling. Calibrate generated baseline/screen only. Selected FGS requires a positive calibrated tolerance; failures never silently retune.
- Launch separate Julia and BLAS settings. Export OMP_NUM_THREADS, OPENBLAS_NUM_THREADS and applicable controls to requested BLAS count before startup. Keep BENCH_BLAS_THREADS and runtime assertions; GEMM before banner, assertions after load/constructors/solves/profiles.
- Clear every constructor timing tuple including excluded fresh trial. Verify independent cache donors and zero initial guesses/history. Keep begin_step_solution! outside solve timing.
- Persist warmup and direct crosschecks. Record authoritative evaluator/residual, FMM certification/residual, direct residual and disagreement separately. BC relative L2 <=1e-6; certified FMM/direct disagreement <=1e-7. Explicit direct acceptance allowed for uncertified FMM on R1/R2.
- Validate CPU/allocation profile solutions outside recorded regions. Verification/requested-profile failures fail run. Keep unprofiled timings separate.
- Mark inferred work as estimates and unavailable counts explicitly. Retained object bytes and process peak RSS are distinct.
- Document commands, boundaries, calibration, failure rules and instrumentation limits.

## HPC preparation and pilot
1. Isolated development and execution worktrees; inventory loaded FLOWPanel/FastMultipole/FLOWVPM and commit required current code separately, preserving unrelated changes.
2. Annotated campaign tags for all three; deploy through git, clean worktrees, dedicated environment; record Manifest and actual loaded paths/tags/SHAs.
3. prep_campaign_worktree.sh creates a data-symlink commit AFTER its input tag: tag and record the actual execution commit.
4. Refresh ORC availability and test submission. Exclusive Zen3, 500 GiB ceiling, explicit affinity, non-preemptible allocation. No resumeless standby run.
5. Sequential control tests at Julia 1 and 4, then R1 smoke on compute node. Both seeds at Julia/BLAS 4/1, calibrate once, preserve selected config.
6. Exact configs at 64/1 and 64/64 in sequential processes on same node; fixed configuration comparisons only.
7. Separate CPU/allocation profiles for both at 64/1, native frames and thread/task grouping, unique generations in consolidated data root.
8. Stop on correctness/thread/provenance/memory failure; resolve before continuing. No automatic R2/R3, screening, instrumentation expansion or optimization.

## Acceptance and reporting
Extend test/runtests_benchmark_cold.jl: invalid inputs have no filesystem effects; selected settings immutable; BLAS stable after actual work; evaluator semantics; failure propagation. Tests and precompilation separate from timings. Actual R1 repeats must demonstrate zero resets, independent constructors, convergence, BC accuracy, memory and fresh/prepared agreement <=1e-8. Adaptive repetitions: 5 below 60s, 3 below 600s, otherwise 2. Report all samples, min/median/spread, prepared and setup+first solve.
Acceptance requires successful R1 timing/profile artifacts, unchanged tolerances, requested/observed threads, provenance exclusively pinned worktrees. Thread-arm convergence failures are failures, never retuning permission.
Final review: fixes, tests, tags/SHAs, jobs, outputs, timing/profile observations, gaps. Use explicitly lower-cost agents with minimal context for bounded test execution, monitoring and harvesting; primary owns edits and interpretation. Offer notebook entry after pilot and obtain approval before writing. Pilot cannot establish optimized winner or complete ILU fairness audit; instrumentation, broader screening, scaling plots and optimization remain subsequent work.
