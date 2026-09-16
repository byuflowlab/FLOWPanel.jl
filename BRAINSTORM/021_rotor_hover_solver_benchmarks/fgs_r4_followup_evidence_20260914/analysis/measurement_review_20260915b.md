# R4 measurement review — 2026-09-15

Read-only source review of FLOWPanel v16 (`546037672ac7e41eba5c64cedf6a2e4d6b7394e6`) and FastMultipole v11 (`adb9967d5b696cf9ab05c557aa2a649e100004dc`). No cluster measurements were performed for this review.

## Conclusion and measurement scope

The generation is suitable for separate coarse-stage activity and whole-prepared-solve generic cache-counter evidence alongside v15's stage budget and thread ladder. It cannot establish DRAM bandwidth or saturation. No numerical or prepared-scope defect was identified.

- Reset, construction, calibration, certification, and solution copying are outside the enabled counter region. Baseline precedes activity and counters, warming the unobserved callback solve before measurement.
- `perf stat` launches `bash benchmark/run_cold_process.sh`; that script uses `exec julia`. Thus the shell becomes Julia rather than leaving a concurrently active wrapper. Default perf inheritance is essential to include Julia worker threads; do not add `--no-inherit`. This is a workload-process/thread measurement, not a system-wide or socket-wide count. Child activity created by the workload could also be inherited; verify actual scope with a bounded Linux multithread smoke if uncertainty remains.
- PMU events use `:u`, excluding kernel execution; task-clock is a separate software event. Report event semantics, support, running/enabled coverage, multiplexing, and errors from the actual perf output. Acknowledgements and exit status alone do not establish useful counts.
- Counters are enabled before the enable acknowledgement returns and disabled after the disable command is processed. Their scope includes handshake tails around the solve. `diagnostic_seconds` starts after enable acknowledgement and ends after disable acknowledgement, so its interval differs from the counter interval. **Handshake overhead has not been measured**; the source discloses it but supplies no numerical bound. Do not call either interval exact solver-only time or use these diagnostic times for speedup claims. A bounded repeated empty enable/disable control can quantify the boundary cost if precision requires it.
- Coarse observer spans enclose complete FMM calls, initialization, influence mapping, residual, nearfield update, final update, and strength copy without changing sweep order. Nearfield update combines leaf solves, products, scatter, and updates; v15 timers provide the finer wall-time decomposition.
- Sequential `/proc` snapshots and clock-tick quantization limit activity estimates, especially for short spans. Aggregate repeated spans, flag incomplete endpoints, and do not infer exact short-stage utilization. Observer runs are separate from the counted run.
- Explicit physical-core cpusets constrain j4/j64; a cpuset still permits migration within its CPUs. Endpoint CPU IDs are not continuous placement traces. Preserve topology, NUMA, allocation, and first-touch evidence.
- Generic cache references/misses do not measure DRAM bytes. Uncore inventory is not bandwidth measurement. Record unavailable suitable counters and leave saturation unresolved; any saturation claim requires a measured bandwidth reference under matching placement.

## Protocol correction recommended before expensive fixtures

Julia `counter_command` uses unbounded `readline` on an `r+` FIFO. A missing acknowledgement can hang until the Slurm timeout. Python smoke has a 15-second timeout, but Julia controls use `IOBuffer`, not actual Linux FIFOs. Add bounded actual-Julia handshake validation before fixture construction, preferably with negative timeout/error coverage. Any executable correction requires fresh annotated pins and clean worktrees; do not rewrite v16/v11 tags.

## Exact completion requirements

1. Audit root completion, all six v15 arm statuses, numerical gates, source/environment provenance, and retained failed-run reasons. Do not rerun the completed ladder merely to reproduce existing evidence.
2. Deliver raw timing/counter/config CSVs, thread-complete text profiles and sample coverage, package pins/loaded paths, reproducible commands, and control outcomes.
3. Report matched alternating batches, order, medians/spread, and instrumentation overhead. Replicate apparent crossovers or small claimed gains independently. Only matched uninstrumented accepted solves support speedup claims.
4. Reconcile an exclusive wall-time budget with total prepared time; include outer/sweep counts, initialization/final work, and remaining work. Do not add nested timers twice or interpret profile sample shares as wall fractions.
5. Report actual GEMV dimensions, matrix bytes/types, leaves, interactions, scatter sizes, percentiles/extremes, and dependency/conflict census. Preserve solve→product→scatter ordering in any proposed batching.
6. Report per-stage activity and the j1/4/8/16/32/64 BLAS=1 ladder, exact physical CPUs, sockets, NUMA policy, allocation/first-touch conditions, and separate measured facts from serial/per-core/NUMA/aggregate-bandwidth hypotheses.
7. Revise the diagnostics package ranking using measured wall time. Its historical ~1.3% main-task FMM share cannot justify low priority: the already-reviewed j1 budget has about 26.7 seconds FMM and 8.16 seconds products. Final ranking awaits the complete ladder.
8. Preserve inner=3 reference, BLAS=1, zero reset, constructor-free prepared timing, finite converged solutions, authoritative certified FMM, BC rel-L2 ≤1e-6, repeat difference ≤1e-8, and direct/FMM disagreement ≤1e-7 when evaluated. NaN unevaluated metrics are not zero.

If the complete stage budget and dependency/scheduling evidence support it, test existing colored ordering as a separately calibrated configuration: color sizes, synchronization/scatter cost, iterations, and total time to accepted accuracy. Otherwise finish with a concrete justified next experiment. Slice fixes and Float32-storage implementations are subsequent isolated experiments, not prerequisites for completing this diagnostics pass. Stop after the evidence handoff; no notebook entry without separate approval.

## Cleanup inventory and boundary

Preserve hash-verified evidence outside the silo for failed jobs 13688396, 13688430, 13688474, 13689318, and 13690579; passed 13690544 controls; pre-v15 provenance; complete v15 results; and any counter-generation successes/failures. Preserve scheduler/process/control logs, statuses, source manifests and verification records, all package tags/SHAs/loaded paths, Project/Manifest, pins, commands, hardware/placement records, and deployment/test evidence.

Only after every job using every generation is terminal and this evidence is verified outside it, delete `/home/rander39/FLOWPanel-p021-r4-diag-v10-silo` and verify removal. Do not follow data symlinks or delete canonical data, shared dependencies, or other silos. Cleanup was already authorized; prerequisites remain evidence-dependent.

## Artifact index

- `../../fgs_opt_r4_diagnostics_handoff_20260912.md`, 2026-09-14 section: governing requirements.
- `../../fgs_r4_context_reset_20260915b.md`: current pins, deployment state, cleanup authorization.
- `../../fgs_r4_followup_validation_20260915.md`: failed-run preservation and v15 provenance.
- `../../fgs_opt_r4_diagnostics_package_20260912.md`: historical ranking requiring measured revision.
- `/private/tmp/flowpanel-p021-r4-counters-v16/benchmark/{fgs_r4_counters.jl,run_r4_counters.slurm.sh,run_cold_process.sh}`: driver, launcher, exec/inheritance boundary.
- Same worktree, `test/{r4_perf_control_smoke.py,runtests_r4_counters_driver.jl}`: Python FIFO smoke and mocked Julia protocol tests.
- `/private/tmp/fastmultipole-p021-r4-activity-v11/src/solve.jl`, lines 1174–1356: observed solve stages.
