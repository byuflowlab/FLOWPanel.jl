# R4 follow-up validation — 2026-09-15

## Refreshed state

Job **13690579 failed**, exit 1, elapsed 45m09s. Its controls passed, but
the first j1 process stopped before measurements. There is no accepted R4
diagnostic result from that job and no completed thread ladder. The user queue
was empty at **2026-09-15T12:35:43Z**.

The full traceback identifies a Julia world-age error: `cold_initialize!`
dynamically includes `phase1_case.jl`, defining `reset_cold!()`, while the
diagnostics `main()` continues in its earlier world. The zero-argument method
exists; this is not a missing argument. The failure occurred after the
58,192-panel direct source assembly (2049.5 s).

## Evidence preserved before changing the silo

Local evidence root: `fgs_r4_followup_evidence_20260914/`.

- Five failed runs: 13688396, 13688430, 13688474, 13689318, 13690579.
  The bounded harvest contains 91 files, each checked against its remote
  SHA256 digest using an explicit path mapping.
- `pre-v15-silo-provenance/`: 19 scheduler-log, environment, pin, and
  source-manifest files, checked the same way. The ten scheduler logs are
  empty; substantive output is in each run's process/control logs.
- `sha256-normalized-verification.txt` records both successful comparisons.
  Older nonempty `.diff` files compare different path prefixes and are not
  verification results.
- The existing passed control evidence for 13690544 remains intact.
- No VTK output was found in these R4 runs. The excluded `fixture-j1`
  directory is empty. No archive action was necessary.

Home storage snapshot: 186180 MiB used, 1910972 MiB available; headroom under
the 400-GiB campaign cap is 223420 MiB.

## v15 preparation

Isolated worktree: `/private/tmp/flowpanel-p021-r4-diag-v15`.

FLOWPanel annotated tag: `campaign/p021-cold-source-20260915-v15`
(`39ec4e3630bc6c04f0865a1d3feecce7130904fe`). FastMultipole remains at
`campaign/p021-r4-diag-source-20260914-v10`
(`87cbc8460b51f24ddf34cc5f41a1d1b6682bf04a`); FLOWVPM remains at
`campaign/p021-cold-exec-20260910-v1`
(`05c658f7804ec5f9b68d4cb9826a9f97cfecb373`).

The v15 FLOWPanel SHA256 content-manifest digest is
`3b2843a4a449eb54859b997cbb6eff7d159a0e9013a170f1bdfe64e0acc14f9c`.
It covers 1486 tracked regular files, excluding the canonical `data/` tree
and the documentation-image symlink, matching the previous manifest scope.

- Split initialization from `run_diagnostics`, crossing the dynamic-loading
  boundary once with `Base.invokelatest`, following `cold_main`.
  This boundary is outside all solve timings.
- Require `certified_fmm` explicitly during compilation, warmup, timed
  trials, and instrumentation-equivalence controls. The generic acceptance
  helper permits a direct fallback and is insufficient by itself.
- Export actual directed leaf dependencies, symmetrized conflict degrees,
  and hypothetical greedy-color group sizes/bytes. These use the same
  write/read-overlap definition as FastMultipole's coloring code; the
  executed solver remains lexicographic.
- Record compute-node hardware-counter capability separately from timings,
  plus `lscpu` and the CPU-tick frequency. Capability probes are not workload
  measurements. Successful PMU access still requires a separately scoped
  workload-counter run.

Local validation passed with one Julia thread: actual `main` AST entrypoint,
a negative world-age control with the boundary removed, three independent
brute-force dependency-overlap cases (asymmetry, non-leaf spans, duplicate
entries, empty target ranges), complete-driver parse, and `bash -n` launcher.
The test is committed at `test/runtests_r4_diagnostics_driver.jl` and runs
first in the cluster launcher. Local Julia is 1.12; the cluster reruns it
under the campaign's Julia 1.11.7 before expensive fixture assembly.

The exact exclusive Zen3/64-core/500-GiB/12-hour request passed Slurm's
non-allocating start estimate at 2026-09-15T12:33:39Z, targeting m12.
Both remote source manifests passed full verification at
**2026-09-15T12:41:47Z**; the data symlink and FLOWVPM pin matched, and the
user queue was empty. Evidence: `v15-deployment/remote-preflight.txt`.

Submitted **13694724**, with no explicit partition and the recorded resource
request, from the verified silo source. Environment overrides:

```text
COLD_PROJECT=/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/env
CAMPAIGN_PINS=/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/pins.toml
COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910
```

Output: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/diag-v15-13694724/`.
Scheduler logs: silo `FLOWPanel.jl/logs/slurm/r4-diag-v15-13694724.out` and
`.err`. Do not edit this source or resubmit while the job is active. Initial
live state and campaign validation results remain to be recorded.

## Remaining completion gates

Require root `COMPLETED`, all six per-arm completed statuses, numerical gates,
verified source/environment provenance, overhead analysis, exclusive-stage
reconciliation, dependency/block census, and matched uninstrumented scaling.
Thread CPU ticks span reset/validation work as well as the timed solve; they
are not per-stage activity measurements. Profile sample fractions are not
wall-time fractions. Bandwidth saturation remains unproven without valid
workload counters.

Delete only `/home/rander39/FLOWPanel-p021-r4-diag-v10-silo` after all jobs
using it are terminal and required results, logs, provenance, and failed-run
evidence are verified outside it. Preserve the canonical data root and do
not follow the silo's data symlink. No notebook entry has been written.


## Resumed validation — 2026-09-15, 19:59 UTC snapshot

Job **13694724 remains RUNNING** on `m12-2-26`, elapsed 5:41:37.
Slurm timestamps are cluster local MDT: start 08:17:32 MDT is **14:17:32 UTC**;
the 12-hour deadline is **2026-09-16 02:17:32 UTC**. Earlier monitor output
incorrectly labeled local timestamps as UTC; use these conversions.

The j1 arm completed. Its 19 files and 24 immutable root evidence files were
harvested to `fgs_r4_followup_evidence_20260914/diag-v15-13694724/`; all **43
files** matched remote SHA256 digests. Manifest and verification records are
in that directory. The j4 arm has produced its census, but no accepted timing
results were available at this snapshot. The other arms and root completion
marker remain outstanding; this is **partial evidence, not campaign acceptance**.

### Accepted j1 diagnostics

| Quantity | Result |
|---|---:|
| Uninstrumented trials | 20 (two batches of 10) |
| Uninstrumented median / range | 37.5594 s / 37.3886–37.6857 s |
| Instrumented trials | 20 (two batches of 10) |
| Instrumented median / range | 37.4689 s / 37.0734–38.2301 s |
| Instrumented minus uninstrumented median | -0.0905 s (-0.241%) |
| Certified FMM BC relative L2 | 4.77952513e-7 |
| Direct BC relative L2, equivalence controls | 4.78150591e-7 |
| Direct/FMM disagreement, equivalence controls | 5.4676859e-9 |
| Repeat and instrumentation solution differences | 0 |
| Iterations / inner sweeps | 27 / 81 |

All 40 timed rows pass the solved/accepted/certified-FMM and repeatability gates;
the instrumentation controls have identical convergence histories and the
thread-complete profile passes separate validation. Unevaluated direct metrics
are NaN, not measured zero. The negative median timing difference is not a
speedup claim: trial ranges overlap, and independent process replication is
absent. No positive median instrumentation overhead is observed in this arm.

Instrumented median stage times are FMM 26.73485 s, nonself product 8.16198 s,
initialization 1.00337 s, scatter 0.79205 s, leaf solve 0.59771 s, residual
0.18580 s, mapping 0.01468 s, and remaining iteration work 0.01515 s.
The internal total is 37.46101 s, exclusive sum 37.44830 s, and internal
unaccounted time 0.01268 s; median outer-solve minus internal-total time is
0.00788 s. Per-trial additive checks pass after accounting for CSV decimal
rounding (largest discrepancy below 1e-7 s). Medians of individual components
need not add exactly to the median total. These are diagnostic stage budgets,
not uninstrumented speed comparisons or profile sample fractions.

Census: 1,068 nonempty leaves, 357,856,254 Float64 matrix elements
(2,862,850,032 bytes), 95,390 directed dependencies and 48,627 undirected
conflicts. All exported degrees agree with an independent edge scan. The
hypothetical greedy coloring has 79 groups of 1–27 leaves (median 16), with
zero same-color conflicts. Execution remains lexicographic and BLAS=1;
colors describe a possible schedule, not an executed speedup. See
`analysis/j1_summary.md` for block distributions and profile coverage.

The compute-node PMU capability probe succeeded for cycles, instructions,
cache references and cache misses. **No workload-scoped counter evidence yet
exists**, so bandwidth saturation is unresolved. A separate local diagnostic
generation is being prepared in `/private/tmp/flowpanel-p021-r4-counters-v16`
and `/private/tmp/fastmultipole-p021-r4-activity-v11`; it is not yet pinned,
deployed, or submitted. It adds coarse stage CPU observations and gated perf
counters around a warmed prepared solve; it does not change solver ordering.
The running v15 source is untouched. Silo cleanup conditions are not yet met.


## Context reset — 2026-09-16 02:16 UTC

Current continuation state is saved in `fgs_r4_context_reset_20260915b.md`.
The last verified v15 job state is stale (19:59 UTC); refresh it immediately.
The counter generation is pinned and locally tested but only partially
deployed, with no job submitted. SSH authentication failed on its second
source transfer; retries stopped. The silo and all canonical data remain.


### SSH restored; job complete — 2026-09-16 02:19:05 UTC

Ryan re-established SSH. Fresh `TZ=UTC sacct` confirms job **13694724
COMPLETED, exit 0:0**, elapsed **08:52:49**, end **2026-09-15 23:10:21 UTC**.
Root `COMPLETED` and all six arm `status.toml` files exist remotely. Remaining
arms have not yet been harvested or independently audited locally. Do not
resubmit v15; next work is full harvesting, hash/provenance/numerical auditing,
and the separate counter diagnostic. The v16 deployment remains partial;
no counter job was submitted. See `fgs_r4_context_reset_20260915b.md` for the
current continuation instructions and delegation strategy.

## 2026-09-16 ~03:15 UTC update (supersedes the 20260915b reset prompt's live state)

- v15 ladder job 13694724: fully harvested locally (all arms + root + scheduler
  logs), 131/131 SHA256 verified, `analysis/audit-full.txt` PASS on all six
  arms. Ladder tables in `.../diag-v15-13694724/analysis/ladder_summary.md`.
- Measured conclusions written as **§7 of
  `fgs_opt_r4_diagnostics_package_20260912.md`** (ranking revision: serial
  leaf-sweep chain ≈85% of j64 wall; FMM far field scales 24×; colored-sweep
  condition SATISFIED).
- v17 counters generation: FLOWPanel `b0b6eec` (tag
  `campaign/p021-r4-counters-source-20260915-v17`, = v16 + bounded poll(2)
  perf-FIFO acks + Julia perf smoke; local driver tests PASS), FastMultipole
  v11 pin unchanged. Deployed to silo `counters-v17/` (fresh generation; v16
  partial dir only seeded the rsync), remote `sha256sum -c` PASS on both
  packages, env/pins/data-symlink/logs installed, FLOWVPM worktree verified
  clean at pin. Storage preflight 353 G/400 G.
- **Job 13711596 submitted 2026-09-16 ~02:55 UTC** (64 CPU exclusive
  `--constraint=zen3`, 6 h, qos=normal). Provenance:
  `fgs_r4_followup_evidence_20260914/v17-deployment/submission-provenance-13711596.md`.
- Remaining: harvest+verify 13711596 when terminal; then the mandatory
  authorized deletion of `/home/rander39/FLOWPanel-p021-r4-diag-v10-silo`
  (whole silo incl. counters-v16/-v17) once ALL its jobs are terminal and all
  evidence is hash-verified locally; then final provenance touch-up. No
  notebook entry without Ryan.
