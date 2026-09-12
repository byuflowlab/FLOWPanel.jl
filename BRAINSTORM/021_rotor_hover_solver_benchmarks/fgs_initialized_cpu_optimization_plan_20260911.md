# Item 021: optimize initialized, CPU-only FGS solves

## User steering, 2026-09-12

Skip R3 and the remaining R2 tuning; proceed to the frozen R4 case (58,192
panels). R2 settings are provisional starting points, not presumed R4 winners.
The immediate goal is R4 baseline measurements, inner-iteration/configuration
diagnostics, and longer CPU/allocation profiles that expose implementation
costs at larger problem size. A separate diagnostics agent will collect and
validate that evidence, then hand it back to the code-optimization agent.
Do not implement solver performance changes in the diagnostics phase. The
acceptance, reproducibility, and eventual code-validation requirements below
remain binding; the original R2 → R3 → R4 sequence is superseded by this
instruction. See `fgs_opt_r4_diagnostics_handoff_20260912.md` for the handoff.

## Objective and evidence

Optimize time to an accepted solution with an already-initialized FGS solver and zero initial strengths. Exclude construction and warm starts.

The completed R1 profile identifies nonself dense matrix-vector products as the first target: 85 of 118 main-task samples reach BLAS `gemv!`, versus 12 reaching cached leaf solves. These overlapping stack counts establish priority, not precise wall-time shares. Prepared R1 medians are 0.6125 s at Julia/BLAS 4/1, 0.5784 s at 64/1, and 0.8758 s at 64/64. R1 remains supporting evidence; optimization proceeds in workload order R2, then R3, then R4.

Pursue several focused optimization passes, including algorithm changes. Treat roughly 5% total solve improvement as significant; replace the prior trajectory and its two-change limit.

## Establish a comparable baseline

- Start from the profile campaign’s pinned implementation, preserving unrelated changes in the live checkout. Use isolated development worktrees and annotated, clean campaign pins for FLOWPanel, FastMultipole, and FLOWVPM.
- Extend the existing cold benchmark harness with a prepared-only execution path. Reuse solver state, reset strengths and mutable solve state before each trial, and exclude compilation, construction, reset, validation, and profiling from solve timing.
- Use the existing frozen R2, R3, and R4 cases and FGS seeds, in that order. Start with R2 and BLAS=1, promote only validated improvements to R3, then validate and tune the promoted configuration on R4. Retain current lexicographic, cached-LU behavior as the control.
- Collect longer repeated profiles and separate component timings for nonself products, scatter, leaf solves, FMM, residual evaluation, and remaining solve overhead. Record outer iterations, actual inner sweeps, FMM passes, allocations, and matrix/leaf dimensions.
- Obtain attribution on R2 first, then R3 and R4; use R1 only as supporting profile evidence. Keep profiler results separate from uninstrumented performance measurements.

## Optimization sequence

1. **Reduce repeated near-field work first.** Screen fixed inner-iteration counts `{1,2,3,5,10}` at otherwise frozen settings. Rank by total time to accepted accuracy, including any extra outer iterations and FMM passes.
2. **Retune the near/far split around the winners.** Test leaf sizes `{25,50,100,200}`, then neighboring expansion orders and MAC values around each existing seed. Use staged screening, retaining the two fastest accepted candidates at each stage rather than a full Cartesian sweep. Recalibrate stopping tolerance through the existing independent accuracy gate.
3. **Test existing colored sweeps.** Compare against lexicographic sweeps at matched settings, then retune inner iterations for the winner. Measure color sizes, parallel product time, serial scatter, and scheduling overhead. Screen Julia threads `{1,4,16,32,64}` on HPC with BLAS=1; use at most four threads locally.
4. **Improve the remaining dominant implementation cost.**
   - If nonself products still dominate, investigate matrix layout, memory traffic, and per-leaf product overhead. Do not batch dependent lexicographic updates across leaves.
   - If scatter becomes significant, fuse its old/new update traversal and cache interaction offsets. For colored execution, evaluate target-owned parallel scatter with fixed accumulation order.
   - If small colors cause excessive scheduling overhead, add serial execution below a measured work threshold.
   - Consider residual threading, reusable FMM scratch, or other FMM changes only when updated attribution predicts a significant total solve gain.
5. After each retained change, rerun attribution and reorder remaining work by estimated total benefit versus implementation effort. Do not assume higher thread count is better.

Keep existing public solver defaults and signatures unless validated results justify a general default change. Store case-specific tuning in benchmark configurations; keep new scratch and scheduling machinery internal.

## Acceptance and verification

- Preserve the existing BC relative-L2 target of `1e-6`, evaluator certification, and evaluator disagreement limit of `1e-7`. Cross-check R2 and R3 directly where feasible; use certified evaluation for R4, with direct evaluation when certification is inconclusive.
- Require finite solutions, successful convergence, and repeated-run solution agreement within the existing `1e-8` threshold. Algorithm variants may follow different convergence histories; validate their final physical accuracy independently.
- Measure finalists with at least ten unprofiled trials per case in matched CPU environments, alternating baseline and candidate batches. Report median, spread, iterations, allocations, and memory. Repeat ambiguous comparisons.
- Retain broadly enabled changes only when their gains exceed measurement noise and they introduce no unexplained regression above 5% on another rung. Keep workload-specific winners explicitly configured.
- Run FastMultipole’s relevant solver tests plus FLOWPanel solver and FGS-history suites. Add tests for any changed update ordering, scratch isolation, and threaded determinism. Cover scalar-potential and gradient paths when shared machinery changes, then run the nearest relevant integration test.

## Completion and deliverables

Stop when every identified reasonable-effort opportunity with plausible ≥5% total benefit has been implemented and validated, experimentally rejected, or shown to require substantial redesign. Reprofile the final combination to check for newly exposed opportunities.

Update item 021 with reproducible commands, pins, baseline/final timing tables, chosen configurations, and a ranked disposition of remaining opportunities. Claim measured CPU solve speedups only; constructor improvements, warm starts, GPUs, and comparisons with other solver families remain outside this effort.
