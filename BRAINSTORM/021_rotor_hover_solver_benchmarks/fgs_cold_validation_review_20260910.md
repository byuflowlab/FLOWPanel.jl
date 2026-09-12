# Cold harness validation / R1 pilot review — 2026-09-10

Status as of 2026-09-12 UTC: all three R1 timing arms passed (60 samples) in job 13653115; that job failed at profile startup. The harness import is repaired in a new v7 pin, and profiles-only continuation 13653450 is submitted. Detailed timing results and continuation status are below. Earlier sections are retained as chronological history.

## Preserved work and scope
Original FLOWPanel/FastMultipole/FLOWVPM dirty checkouts and historical outputs were left intact. Development worktrees: `/tmp/flowpanel-cold-20260910`, `/tmp/fastmultipole-cold-20260910`, `/tmp/flowvpm-cold-20260910`. Runtime source was copied into dedicated snapshot commits (FLOWPanel `575f2a0`, FMM `720472ee`, FLOWVPM `d7402d4`); no scientific/source optimization edits were made in this task. FLOWPanel source snapshot includes existing formulation/wake/warmstart changes, FMM includes existing farfield output hooks, FLOWVPM includes existing splitting/timeintegration changes. Unrelated examples, tests, notebooks and results were not included.

## Harness implementation
FLOWPanel implementation `b7e682e`, annotated source tag `campaign/p021-cold-source-20260910-v1`. Dependency source tags: `campaign/p021-cold-20260910-v1` in each dependency repository. Published through git; execution tags will include data-symlink/site-policy infrastructure commits.

Changes: pure input/schema/path/memory/thread preflight; immutable selected settings; pre-start BLAS/OpenMP controls and real GEMM thread check; constructor/reset checks; persisted direct/FMM/authoritative residual evidence; fresh tuple release; convergence/profile validation; fail-fast execution; explicit estimated/unavailable work counts; separate retained memory and lifetime RSS; usage docs; compute-only sequential launcher.

Static verification: shell syntax and `git diff --check` passed. Julia has not run locally.

## HPC preparation
Current user jobs 13593020/13593021 use `/home/rander39/wt021/FLOWPanel.jl-c`; that checkout was not edited. Exact pilot resource request passed `sbatch --test-only`: 1 node, 1 task, 64 CPUs, exclusive Zen3, 500 GiB, normal QOS, 4 hours, no explicit partition (site policy). Scheduler selected m12 as eligible; private idle Zen3 capacity was standby-only and excluded. Current ORC site policies were read over authorized SSH; initial sandbox-only connection restriction was resolved through escalation.

Campaign root: `/home/rander39/campaigns/p021-cold-20260910-v1`.
Output root: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910`.

Pending: actual execution pins, job IDs, control/smoke outcomes, fixed-config timings, profiles, final acceptance review. No optimized winner or completed ILU fairness audit is claimed. Notebook entry has not been written; approval remains required.

## Published execution pins and submitted job

All three repositories use annotated `campaign/p021-cold-exec-20260910-v1`:

| Package | Execution SHA |
|---|---|
| FLOWPanel | `029b9f3d30bba706dea77418281f0ce4d75a4689` |
| FastMultipole | `ef10643a401d6da16e28be87805b67d11bdf1fb5` |
| FLOWVPM | `05c658f7804ec5f9b68d4cb9826a9f97cfecb373` |

FLOWPanel's helper made data-symlink commit `361bf22` after input source tag; the
execution pin additionally records the ORC policy symlink. Dependency execution
commits add the policy reference/symlink only. Cluster GitHub push lacked interactive
credentials; execution tags were fetched over existing SSH and published using the
authenticated local git client. No credentials were copied.

Job `13642618` uses `benchmark/run_cold_pilot.slurm.sh` from the pinned FLOWPanel
worktree. Environment: `.../p021-cold-20260910-v1/env`; pins file:
`.../p021-cold-20260910-v1/pins.toml`. Logs are
`/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/slurm-13642618.{out,err}`;
process outputs/logs use `.../p021-cold-20260910/pilot-13642618/`. Requested 4 hours,
normal QOS, Zen3, exclusive, 500G, 64 CPUs; no partition was forced. The final dry
submission check accepted the request on m12. Each Julia process is bound to the
same explicit mask of 64 distinct physical cores inside the allocation.

## Pause checkpoint — 2026-09-11

Ryan requested a break. Job 13642618 failed before Julia (`HISTCONTROL` unset in
`/etc/profile` under `set -u`); no Julia tests or performance artifacts exist.
Launcher shell initialization is fixed and checked without executing Julia.
Latest clean local development HEAD is `3194a2c`; v3 cluster preparation completed
at FLOWPanel execution SHA `e3ae9bd89b51d1cb86500b623a2a99195cdffa04`, but no v3 job
was submitted. Its execution tag remains cluster-only. A likely Julia parse defect
in the unquoted `$tag^{commit}` command literal was noticed immediately before
pause and remains unfixed. See `fgs_cold_resume_prompt_20260911.md` for authoritative
continuation instructions and complete paths/pins. No background work from this
task remains intentionally active, and no notebook entry was written.

## 2026-09-11 continuation: validation passed; CPU timing/profile job launched

This section supersedes earlier pending/unvalidated startup status. The detailed next-agent handoff is `fgs_cold_profile_resume_prompt_20260911.md`.

The documented Julia command-literal defect was fixed in `56411b9`, with a compute-node syntax gate added before precompilation. Runtime validation then exposed two further startup defects: missing CUDA toolkit discovery during extension precompilation (`5355045` loads the pinned site toolkit), and numeric promotion of generated configuration axes (`3acf1c4` preserves their types with tuples). CUDA is used solely for precompilation dependencies; Ryan explicitly restricted the investigation to CPU execution. Strict validation was preserved, and no solver algorithm or optimization changed.

Latest local implementation is clean `/tmp/flowpanel-cold-20260910` at `3acf1c4cff71fdee3dfae95e42d10e21012aea7e`, source tag `campaign/p021-cold-source-20260911-v6`. Actual execution FLOWPanel pin is `campaign/p021-cold-exec-20260911-v6` (`08a0247d173ddd6dcfd8a16b418002de31e3ca5b`), worktree `/home/rander39/campaigns/p021-cold-20260911-v6/FLOWPanel.jl`. Dependencies retain v1 execution tags: FastMultipole `ef10643a401d6da16e28be87805b67d11bdf1fb5`; FLOWVPM `05c658f7804ec5f9b68d4cb9826a9f97cfecb373`. All execution tags are annotated/published; the dedicated Manifest points to the three clean pinned worktrees.

Job **13644535 completed successfully**, exit 0:0, 7m44s on m12-3-14. The dependency-free syntax gate and precompilation passed. Controls passed **530/530 at Julia/BLAS 1/1 and 530/530 at 4/1**. Each total comprises 12 BLAS, 424 configuration/immutability, 60 invalid-input/filesystem, 14 evaluator, and 20 failure-propagation checks.

Both R1 smoke solvers passed all three repeat/independent-constructor trials:

| Solver | Iterations | Certified BC relative L2 | Direct BC relative L2 | FMM/direct disagreement | Maximum solution difference |
|---|---:|---:|---:|---:|---:|
| FGS | 6 | 9.53567919e-8 | 9.57141435e-8 | 3.92855040e-9 | 0 |
| ILU-GMRES | 7 | 8.57921900e-7 | 8.57157852e-7 | 3.92854753e-9 | 0 |

Every smoke row passed convergence, certification, independence, reset, accuracy, memory, and agreement gates. Retained body+solver sizes were 211,972,080 bytes (FGS) and 115,103,934 bytes (ILU); maximum recorded process-lifetime RSS was 1,603,387,392 bytes. These quantities have different scopes and are not interchangeable with scheduler-sampled RSS. All smoke rows and clean execution pins were rechecked before the timing submission. Smoke timings are not presented as the final benchmark comparison.

Authoritative evidence: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/pilot-13644535/`, especially controls logs, `smoke-j4-b1/{fgs,krylov_ilu}/smoke.csv`, statuses, and provenance. Frozen configuration: `smoke-j4-b1/selected.toml`, SHA-256 `c596a5f895369f04fc1cfd9c8e5b40a3169b2a5d7ff1607ab98c2f3bb9cdf064`. FGS tolerance is `1.7658196194548e-8`; ILU-GMRES keeps rtol `1e-6`. No recalibration is permitted in timing/profile arms.

Full CPU timing/profile job **13653115 was submitted at approximately 2026-09-12 02:10 UTC** (2026-09-11 20:10 America/Boise) after exact test-only acceptance. Request: one-hour non-preemptible test QOS, exclusive Zen3, 64 requested CPUs, 500G, no partition specified or GPUs. The driver binds 64 distinct physical cores and runs timing arms 4/1, 64/1, 64/64, then separate FGS and ILU CPU/allocation profiles at 64/1, sequentially within one allocation. The one-hour budget was sized from the completed 7m44s controls/smoke job and observed small R1 solves; profile overhead and high-thread behavior remain unverified.

Logs are `data/p021-cold-20260910/slurm-13653115.{out,err}` and generation `pilot-13653115/` under the cluster FLOWPanel project. **Submission is not timing/profile validation.** Remaining work is to verify all arm/status/validation/config/provenance artifacts, harvest every sample and min/median/spread, inspect both CPU/allocation profiles, and report measured bottlenecks with appropriate limitations. No optimized-winner, measured component-work, scaling, or complete ILU fairness claim exists. All historical worktrees/results were preserved; no notebook entry was written.

Final launch check: 13653115 was RUNNING on m12-3-5 at elapsed 3m30s, with the first 4/1 timing generation present and an empty Slurm error log. Later arms and profiles were not yet verified when Ryan requested the context reset.

## 2026-09-12 UTC: verified timing results and profile recovery

Job **13653115** ran on **m12-3-5**, AMD EPYC 7763 (Zen3), from 2026-09-12 02:10:41 to 02:30:39 UTC (2026-09-11 20:10:41–20:30:39 America/Boise), elapsed **19m58s**, and ended **FAILED, exit 1:0**. All three timing processes completed before the profile startup failure; the allocation as a whole is not a successful pilot. Slurm stdout/stderr are empty because process output is redirected into per-arm logs. `COMPLETED` is absent. Scheduler MaxRSS was unavailable in the initial accounting query; process-lifetime RSS is separately recorded below.

The timing processes used Julia **1.11.7**, requested/observed Julia/BLAS **4/1, 64/1, 64/64**, one exclusive node, 64 requested CPUs (128 allocated logical processors), 500G memory, non-preemptible `qos=test`, one-hour limit, no partition or GPU request. Affinity is 64 distinct physical cores. BLAS environment controls and runtime GEMM checks match each arm. All timing provenance reports clean loaded package worktrees and the v6 execution pins listed above, with Manifest SHA-256 `d8b46abf61005f76044625c0624a7042bdedf008747fc45cfbbb37897b0ba612`.

Frozen selected config remains byte-identical to smoke: SHA-256 **`c596a5f895369f04fc1cfd9c8e5b40a3169b2a5d7ff1607ab98c2f3bb9cdf064`**. FGS config ID `409bb3aaf67bbd18`; ILU-GMRES config ID `d569863eedc13014`. No recalibration or solver changes occurred.

### Timing scope and all summary statistics

Prepared = compiled, already constructed solver/cache, zero-reset isolated `_solve!`. Fresh = compiled constructor plus first `_solve!` on a newly constructed solver. Reset, fixture/RHS preparation, independent BC diagnostics and convergence recording are excluded. Fresh is not a full application or process startup time. Compile-validation, excluded warmup, and fresh warmup are separate records, not benchmark samples. All 12 groups have **5 samples**, as required by their excluded warmup times below 60 s. Recomputed min/median/max/spread agree with `summary.csv`.

| Julia/BLAS | Solver | Scope | n | Min (s) | Median (s) | Max (s) | Spread (s) |
|---|---|---|---:|---:|---:|---:|---:|
| 4/1 | FGS | prepared | 5 | 0.609735 | 0.612503 | 0.652487 | 0.042751 |
| 4/1 | FGS | fresh | 5 | 19.060163 | 19.072911 | 19.125670 | 0.065507 |
| 4/1 | ILU-GMRES | prepared | 5 | 11.608074 | 11.619105 | 11.660354 | 0.052279 |
| 4/1 | ILU-GMRES | fresh | 5 | 15.451471 | 15.545275 | 15.607489 | 0.156018 |
| 64/1 | FGS | prepared | 5 | 0.563589 | 0.578413 | 0.585494 | 0.021904 |
| 64/1 | FGS | fresh | 5 | 18.847121 | 18.904605 | 18.982933 | 0.135812 |
| 64/1 | ILU-GMRES | prepared | 5 | 1.350179 | 1.361472 | 1.499375 | 0.149196 |
| 64/1 | ILU-GMRES | fresh | 5 | 4.573305 | 4.624585 | 4.644882 | 0.071577 |
| 64/64 | FGS | prepared | 5 | 0.778418 | 0.875811 | 0.904556 | 0.126138 |
| 64/64 | FGS | fresh | 5 | 19.837486 | 19.993121 | 21.808331 | 1.970845 |
| 64/64 | ILU-GMRES | prepared | 5 | 1.346065 | 1.456084 | 1.583649 | 0.237584 |
| 64/64 | ILU-GMRES | fresh | 5 | 4.507552 | 4.644278 | 4.829430 | 0.321878 |

All 60 raw samples and excluded warmup/compile records are retained in [the evidence directory](fgs_cold_pilot_evidence_20260912/README.md); authoritative raw CSV/TOML/log files remain on ORC under `data/p021-cold-20260910/pilot-13653115/`.

### Accuracy and memory acceptance

All six solver/thread leaves report successful status; all measured trials are solved, accepted, eligible and repeatable. Each leaf has one compile validation, one warmup, one fresh warmup, ten measured trials, and a separate convergence record/validation. All measured fresh/prepared deltas are zero. Independent direct/evaluator comparisons are performed outside measured regions; timed-trial direct fields are intentionally NaN and are not counted as extra direct checks.

| Solver | Iterations | Authoritative certified BC relative L2 | Maximum direct BC relative L2 | Maximum FMM/direct disagreement | Retained body+solver bytes (timed trials) |
|---|---:|---:|---:|---:|---:|
| FGS | 6 | 9.53567919e-8 | 9.57141436e-8 | 3.92855040e-9 | 211,972,080 |
| ILU-GMRES | 7 | 8.57921900e-7 | 8.57157852e-7 | 3.92854753e-9 | 115,103,934 |

These pass BC relative L2 ≤1e-6, certified evaluator disagreement ≤1e-7 and fresh/prepared agreement ≤1e-8. `fmm_rel_max` is a residual infinity norm and must not be mistaken for evaluator disagreement. Independent convergence traces contain seven rows per solver/thread leaf; terminal internal residuals are 5.11660132e-9 (FGS) and 8.53503058e-7 (ILU), with separate authoritative BC validation. Internal residual definitions differ between solvers.

| Julia/BLAS | Solver | Scope | Median constructor (s) | Median solve (s) | Median allocated bytes | Median GC (s) | Max process-lifetime RSS bytes in samples |
|---|---|---|---:|---:|---:|---:|---:|
| 4/1 | FGS | prepared | 0.000000 | 0.612503 | 7,428,312 | 0.000000 | 1,286,459,392 |
| 4/1 | FGS | fresh | 18.460905 | 0.611184 | 445,012,112 | 0.068809 | 1,795,244,032 |
| 4/1 | ILU-GMRES | prepared | 0.000000 | 11.619105 | 15,909,096 | 0.006205 | 1,795,244,032 |
| 4/1 | ILU-GMRES | fresh | 3.924745 | 11.622434 | 1,764,887,360 | 0.467142 | 1,845,620,736 |
| 64/1 | FGS | prepared | 0.000000 | 0.578413 | 18,141,928 | 0.011497 | 1,325,727,744 |
| 64/1 | FGS | fresh | 18.321465 | 0.576861 | 455,848,488 | 0.089592 | 1,759,100,928 |
| 64/1 | ILU-GMRES | prepared | 0.000000 | 1.361472 | 68,834,552 | 0.020788 | 1,759,100,928 |
| 64/1 | ILU-GMRES | fresh | 3.171055 | 1.448700 | 1,818,220,288 | 0.455987 | 1,768,755,200 |
| 64/64 | FGS | prepared | 0.000000 | 0.875811 | 18,141,896 | 0.000000 | 1,312,411,648 |
| 64/64 | FGS | fresh | 19.052314 | 0.938215 | 455,848,424 | 0.067160 | 1,736,876,032 |
| 64/64 | ILU-GMRES | prepared | 0.000000 | 1.456084 | 68,837,144 | 0.008694 | 1,736,876,032 |
| 64/64 | ILU-GMRES | fresh | 3.179026 | 1.465252 | 1,818,221,536 | 0.457596 | 1,842,569,216 |

Medians of setup and solve need not sum to the median of their per-trial total. Allocation bytes are cumulative allocations, not live memory. Lifetime RSS can include earlier solvers and validation in the same process; it is not per-solver incremental memory. Across all timing validation records, maximum lifetime RSS is 1,958,821,888 bytes; ILU retained size reaches 115,104,222 bytes when convergence history is retained. Both remain below the 536,870,912,000-byte eligibility ceiling.

### Bounded interpretation

At 64/1, FGS prepared median is 0.578413 s versus ILU-GMRES 1.361472 s, while fresh medians are 18.904604 s versus 4.624585 s. FGS's measured constructor median is 18.321465 s versus 3.171055 s for ILU. The scope changes the comparison; this pilot does not establish a general winner.

Increasing Julia threads from 4 to 64 with BLAS fixed at 1 changes prepared medians from 0.612503 to 0.578413 s for FGS, and 11.619105 to 1.361472 s for ILU. Raising BLAS to 64 at Julia 64 changes FGS prepared median to 0.875811 s and ILU to 1.456084 s. These are observations from a single-node, fixed-order pilot, not a scaling campaign or causal attribution to a particular component. Full ILU fairness, component-work instrumentation, broader screening, scaling plots and optimization remain outside scope.

### Profile startup failure and immutable continuation

FGS profile process failed before fixture/provenance/profile artifacts with `UndefVarError: Allocs not defined in Main`, at legacy `Allocs.@profile` line 130 inside the top-level conditional at line 25. Julia expands macros in both branches before executing the branch-local import; syntax-only parsing did not exercise this. Consequently job 13653115 contains **no CPU/allocation profiles for either solver**, and ILU profiling was never started.

Repair source commit **`7cda26f`**, annotated source tag **`campaign/p021-cold-source-20260912-v7`**, moves `import Profile.Allocs` before that conditional. A `profiles_only` launcher stage runs the existing two profile processes without repeating accepted timings. These are the only two changed code files; shared solver, constructor, timing and validation code is unchanged. Shell syntax and diff checks pass; actual profile execution is the compute-node validation. No Julia ran locally or on a login node.

New immutable execution tree: `/home/rander39/campaigns/p021-cold-20260912-v7/FLOWPanel.jl`, annotated execution tag **`campaign/p021-cold-exec-20260912-v7`**, SHA **`c9e33d29d9dd21b9fbd2424d578c5ce4f0ecd6b6`**. Dedicated `env/` and `pins.toml` are alongside it; the same two v1 dependency worktrees/tags/SHAs are retained. Source and execution tags were published to the previously authorized GitHub repository. Clean worktrees, exact three Manifest dev paths, and the frozen config hash were checked before submission.

Exact `sbatch --test-only` request was accepted (test identifier 13653446, not a submitted job). Profiles-only job **13653450** was submitted with the same CPU-only exclusive Zen3 request and 1-hour test QOS; it is running on **m12-3-5**. Outputs: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/pilot-13653450/`. Profile validation, final status and measured profile findings remain pending. All worktrees, environments, previous results and failure logs are preserved. No notebook entry has been written.

### Final reset observation (2026-09-12 about 03:00 UTC)

Final scheduler check: **13653450 RUNNING**, elapsed **3m33s**, node **m12-3-5**, no end time, no `COMPLETED` marker yet; ILU profile process directory exists and its log was still empty. No persistent monitoring loop remains.

FGS leaf already has `status = "completed"`, both validation CSVs, CPU flat/tree text (~857K/~707K), CPU raw profile (~1.4M), allocation raw profile (~3.4M), allocation text (~909K), and unprofiled trial CSV. Allocation header reports `sample_rate=0.01; sampled_bytes=165856`. These were only listed/read in place, not yet harvested or fully validated.

FGS log includes two **“no samples collected” warnings**, but CPU artifacts are substantial and nonempty. Do not assume the entire CPU profile is empty or submit a rerun based on the warning alone: inspect thread/task groups and actual samples first (empty groups may explain warnings; this is an unverified hypothesis). Profile gate values, sample coverage and bottleneck interpretation remain next-agent work.

Ryan requested a context reset at this point. Continue from `fgs_cold_profile_resume_prompt_20260912.md`; no further execution change is authorized by the reset itself.
