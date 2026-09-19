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

## 2026-09-16 ~04:25 UTC — 13711596 FAILED at smoke gate; fixed, resubmitted as 13712587 (v18)

- **Job 13711596 FAILED in 40 s** (03:43:30–03:44:10 UTC, m12-4-11, exit 1):
  `test/r4_perf_control_smoke.jl` line 2 was a triple-quoted string directly
  above `using Test` — Julia treats it as a docstring for the `using` statement
  and errors (`cannot document the following expression`). No arms ran.
- Why local testing missed it: the local v17 gate ran only
  `runtests_r4_counters_driver.jl`; the smoke script's sole caller is the Slurm
  launcher, so it was never executed pre-submit. The v17 "smoke passed locally"
  claim was inaccurate.
- Failed evidence harvested + SHA256-verified:
  `fgs_r4_followup_evidence_20260914/counters-v17-13711596-FAILED/`.
- Fix = commit `6033dbd`, tag `campaign/p021-r4-counters-source-20260915-v18`
  (docstring → comment; launcher names v17→v18). This time the smoke script was
  EXECUTED locally against a mock perf FIFO responder (exit 0, correct 3-command
  ack sequence) and the driver suite re-passed; no sibling files share the
  pattern.
- Fresh generation `counters-v18/` deployed (cp -al seed from verified v17 +
  delta rsync), remote `sha256sum -c` PASS (1,491 + 6,441 files), Manifest
  dev-paths verified at v18 paths, FLOWVPM worktree clean at pin, storage
  365 G/400 G.
- **Job 13712587 submitted 2026-09-16 ~04:20 UTC** (same request: 64 CPU
  exclusive zen3, 6 h, qos=normal). Provenance:
  `fgs_r4_followup_evidence_20260914/v18-deployment/submission-provenance-13712587.md`.
- Silo cleanup obligation now covers counters-v16/-v17/-v18; unchanged
  otherwise. Handoff steps 1–5 of `fgs_r4_context_reset_20260916.md` now target
  13712587.

## 2026-09-16 ~19:20 UTC — v18 (13712587) FAILED at smoke gate; root-caused; v19 job 13733217 submitted

- **Job 13712587 FAILED in 41 s** (04:13:00–04:13:41 UTC, exit 1) at the perf
  FIFO smoke gate: `ERROR: LoadError: bitcast: argument size does not match
  size of target type` in `counter_command` (`reinterpret` at driver line 67,
  called from `r4_perf_control_smoke.jl:25`). perf-smoke.csv all
  `<not counted>`; no controls or arms ran; both scheduler logs empty.
- Failed evidence harvested + SHA256-verified (13 files, all OK; FIFOs
  excluded): `fgs_r4_followup_evidence_20260914/counters-v18-13712587-FAILED/`.
- **Root cause (reproduced locally):** `reinterpret(Cint, fd(acknowledgement))`
  is a 32-bit bitcast. On Julia ≤1.11, `fd(::IOStream)` returns a 64-bit `Int`
  → bitcast fails; on ≥1.12 it returns a 32-bit `RawFD` → passes. The cluster
  module is julia/1.11.7; the local v18 gate ran on 1.12.5, so executing the
  smoke script locally could not catch it. Reproduced exactly under local Julia
  1.11.8 with a mock FIFO responder; same error, same line.
- Fix = commit `f1fad9e61b90cb8216f4d7c3aae8ddbbc102dc10`, tag
  `campaign/p021-r4-counters-source-20260916-v19` (branch on the returned fd
  type: `RawFD` → reinterpret, otherwise `Cint(raw)`; launcher names v18→v19;
  no other changes). Smoke script + mock responder and the driver suite now
  PASS under BOTH Julia 1.11.8 and 1.12.4 locally. Lesson recorded: version-gate
  local checks against the cluster Julia line (1.11.x), not just "execute the
  script".
- Fresh generation `counters-v19/` deployed (cp -al seed from verified v18 +
  delta rsync of git-tracked lists); remote `sha256sum -c` PASS (FLOWPanel
  1,491, FastMultipole 6,441 — FMM digest unchanged from v18); Manifest
  dev-paths verified at v19 paths; data symlink → canonical root; pins at
  `counters-v19/pins.toml`. FastMultipole/FLOWVPM pins unchanged.
- **Storage flag:** /home/rander39 preflight shows 725 G used (2.0 T mount,
  1.3 T free) — well over the 400 G policy cap (v18 preflight was 365 G).
  Not a blocker for this small-output job; flagged for an hpc-storage sweep.
- **Job 13733217 submitted 2026-09-16 ~19:15 UTC** (request unchanged: 64 CPU
  exclusive zen3, 6 h, qos=normal). Provenance:
  `fgs_r4_followup_evidence_20260914/v19-deployment/submission-provenance-13733217.md`.
- Silo cleanup obligation now covers counters-v16/-v17/-v18/-v19. Handoff steps
  1–5 of `fgs_r4_context_reset_20260916b.md` now target 13733217.

## 2026-09-16 ~19:35 UTC — v19 (13733217) FAILED at smoke gate; perf NUL-ack root cause measured; v20 job 13733332 submitted

- **Job 13733217 FAILED in 32 s** (19:12 UTC, m12-4-12, exit 1): `perf did not
  acknowledge enable` at the ack-content check. The v19 fd fix worked (the
  first ack was read); the failure moved one command deeper.
- **Root cause measured on the cluster** (login-node `od -c` probe of the real
  perf 5.14.0-570.136.1.el9 control FIFO): each ack is **5 bytes `ack\n\0`** —
  perf writes the tag with its NUL terminator. A fixed 4-byte reader gets
  `ack\n` on the first command and `\0ack` on the next. The v18/v19 mock
  replied a clean 4-byte `ack\n`, so no local gate — even under Julia 1.11 —
  could catch this; only the real perf binary exhibits it.
- Failed evidence harvested + SHA256-verified (13 files OK):
  `fgs_r4_followup_evidence_20260914/counters-v19-13733217-FAILED/`.
- Fix = commit `5c1123b1d1da312f8fb7c26f93d794b22e2ed298`, tag
  `campaign/p021-r4-counters-source-20260916-v20` (NUL-skipping ack reader,
  bounded at 16 raw reads; NUL-stripping fallback; launcher v19→v20).
- Verification now includes the missing gate: **end-to-end on the cluster
  against real perf** (smoke script + patched driver under
  `perf stat --control` with module julia/1.11.7 in a throwaway temp dir) →
  PASS, exit 0. Plus local mocks in BOTH ack framings under Julia 1.11.8 and
  1.12.4, and the driver suite under both. Gate lesson chain: execute the
  script → under the cluster Julia → against the real perf binary.
- Fresh generation `counters-v20/` deployed (cp -al seed from verified v19 +
  delta rsync); remote `sha256sum -c` PASS (1,491 + 6,441); Manifest dev-paths
  at v20 paths; pins/symlink/logs verified. FMM/FLOWVPM pins unchanged.
- **Job 13733332 submitted 2026-09-16 ~19:30 UTC** (request unchanged).
  Provenance:
  `fgs_r4_followup_evidence_20260914/v20-deployment/submission-provenance-13733332.md`.
- Silo cleanup obligation now covers counters-v16 through -v20. Handoff steps
  1–5 now target 13733332.

## 2026-09-17 — v20 job 13733332 COMPLETED; all gates PASS; §7 updated; silo cleanup executed

- **Job 13733332 COMPLETED** (started 2026-09-17 01:27 sacct-local, elapsed
  44 min 26 s). Smoke PASS; all four controls PASS (2216/2216, 11/11, 3/3 +
  counter-driver PASS); both arms `status=completed`, solved/finite, 27
  iterations, BC rel-L2 4.78e-7, repeat delta 0, identical history, certified
  FMM authoritative, BLAS=1, zero-reset counter scope, no perf multiplexing.
- Evidence harvested + SHA256-verified 46/46:
  `fgs_r4_followup_evidence_20260914/counters-v20-13733332/` with analysis
  tables at `analysis/counters_summary.md` (script `analysis/tabulate.py`).
- **Headline finding:** direct /proc thread-activity proves the
  `nearfield_update` chain executes at avg ≈1.00 active thread at BOTH j4 and
  j64 with a non-scaling ~9.5 s span, while `fmm` runs at 3.9/31.0 active
  threads and scales 6.3× — direct confirmation of §7's serial leaf-sweep
  inference; spans reconcile with the v15 wall budget. Counters: IPC 3.59 (j4)
  / 3.26 (j64), miss ratio 5.8/6.0%, cache-ref volume arm-invariant. DRAM
  bandwidth saturation unresolved (pre-declared generic-counter limitation);
  colored sweep proceeds on the serial-execution evidence.
- §7 of `fgs_opt_r4_diagnostics_package_20260912.md` updated: pending-13711596
  paragraph replaced with measured v20 findings + the v17/v18/v19 failure
  chain (docstring / fd bitcast / NUL-terminated perf ack).
- Silo cleanup (pre-authorized): all four silo jobs terminal
  (13711596, 13712587, 13733217, 13733332), all evidence hash-verified local →
  deleted `/home/rander39/FLOWPanel-p021-r4-diag-v10-silo` entirely
  (counters-v16…v20), canonical data root and campaign worktrees untouched.
- Remaining per handoff step 6: STOP — colored-sweep experiment is a new
  campaign phase requiring Ryan; notebook entry awaits Ryan's approval.

## 2026-09-17: v21 colored-sweep A/B campaign launched (job 13738561)

Ryan directed the colored-sweep experiment (handoff
`fgs_r4_context_reset_20260917.md`). Executed:

- Plumbing verified end-to-end (no harness changes needed): config
  `"sweep_order"` → `cold_check_config` → `cold_make` →
  `FGSSolver(sweep_order=…)` → `FastMultipole.FastGaussSeidel` colored sweep.
- New v21 harness (commit `5f87a0e`, tag
  `campaign/p021-r4-colored-source-20260917-v21`): driver
  `benchmark/fgs_r4_colored_ab.jl` (AB_MODE = calibrate / trials / activity),
  launcher `benchmark/run_r4_colored_ab.slurm.sh` (controls → j64 colored
  calibration → alternating A/B trials j∈{1,4,16,64}, 40 uninstrumented
  trials/order/arm → j64 activity pair; perf dropped entirely), unit test
  `test/runtests_r4_colored_ab_driver.jl`.
- Pre-submit gate (v17–v19 lesson chain): every new script executed
  END-TO-END locally at R4 under Julia 1.11.8 (all three modes; unit test
  also under 1.12.4). Colored calibration: tolerance 5.223e-7 (lex retained
  3.479e-7), **26 iterations vs 27 lex** — iterate path moved as predicted;
  repeat delta 0; all gates green. Local-only signal (not citable): colored
  median 9.62 s vs lex 11.51 s at j4/M-series; cross-order solution rel-L2
  1.53e-6.
- Deployment per no-more-silos: pinned git worktrees at
  `/home/rander39/campaigns/p021-r4-colored-20260917-v21/` (FLOWPanel
  `5f87a0e7…`, FastMultipole `c6185cdd…` = merged tip, new tag
  `campaign/p021-r4-colored-fm-20260917-v21` since the validated line is the
  MERGE, not pre-merge v11; FLOWVPM retains `…cold-exec-20260910-v1`
  worktree). Tags pushed to the CLUSTER clones only (origin push of merged
  branches still awaits Ryan). Env resolved under julia/1.11.7; Manifest
  dev-paths verified. Provenance:
  `fgs_r4_followup_evidence_20260914/v21-deployment/submission-provenance-13738561.md`.
- Storage preflight: hpc-storage cycle archived 0 MB — /home at ~654 G of
  400 G cap, all reclaimable mass (~369.5 GiB VTK, 15 runs) is RECENT and
  awaits Ryan's `--include-recent --only` approval. v21 writes only
  CSV/TOML (~0.3 GB), so submission proceeded with the breach flagged.
- **Job 13738561** submitted 2026-09-17 (64c zen3 exclusive 500G normal
  12 h); `--test-only` ETA 2026-09-17T17:16 UTC on m12. Watch armed
  (30-min sacct spacing, defensive parsing).
- Ledger (storage): du /home 659712→669785 MB during cycle (live writers);
  archive tier unchanged 4.027T/93,345 files; STALE=0, VERIFY_FAIL=0.

Stop conditions honored: A/B results go to Ryan before follow-ons; no
notebook entry without approval; no origin pushes.

## 2026-09-17 (later): 13738561 FAILED (env-only) → resubmitted as 13738665

Job 13738561 started immediately via backfill and died in the
`runtests_unit_solver.jl` control: `Meshes`/`StaticArrays` are FLOWPanel
deps but were not DIRECT deps of the campaign env (`--project=env` cannot
resolve a test script's `import Meshes`). Earlier local validation had run
these tests under `--project=.`, masking the gap. Evidence harvested +
SHA256-verified (`colored-v21-13738561-FAILED/`, 17 files); reproduced
locally; fix env-only (`Pkg.add Meshes StaticArrays` both envs — no source
change, no new tags, worktrees untouched); re-gated locally under 1.11.8
(unit solver + history + FMM chain incl. 2216 colored cases all PASS).
**Resubmitted: job 13738665.** Gate lesson appended to the provenance
postmortem: "execute every launcher-invoked script" means under the
campaign env itself.

## 2026-09-17 (later still): 13738665 COMPLETED — A/B landed, all gates green

Job 13738665 COMPLETED (started 06:56 UTC, 10:28 h wall, m12). Sentinel
present, all six stage `status.toml` = completed. Harvested to
`colored-v21-13738665/` (103/103 files SHA256-verified); full analysis in
its `analysis/ab_summary.md`. §6 gates green on all 320 trials; colored
deterministic at 26 iterations every arm (lex 27); calibration matched the
local gate exactly (tol 5.223e-7 vs lex 3.479e-7).

Uninstrumented solve medians (lex → colored, s): j1 38.06→36.77 (= the
26/27 iteration ratio), j4 17.58→15.53 (−11.6%), j16 12.29→10.12
(−17.7%), j64 11.20→11.95 (+6.7%: colored LOSES at j64; lex j64 min
10.97 reproduces the v15 yardstick). **New best operating point: colored
@ j16 = 10.12 s** (~8% under the old lex@j64 champion on 1/4 the cores);
far from the 2.2 s Amdahl bound. j64 activity pair (attribution only):
colored raises nearfield_update avg active threads 0.99→42.65 but the
span does not shrink (10.15→11.02 s; ≈470 busy-CPU-s of barrier spin) —
the 79-color / median-16-width schedule caps realizable parallelism, and
j16 matches the color width.

Reported to Ryan; STOPPED per handoff (no follow-ons, no notebook entry,
no origin pushes without him).
