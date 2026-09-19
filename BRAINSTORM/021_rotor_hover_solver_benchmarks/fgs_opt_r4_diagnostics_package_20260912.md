# BRAINSTORM 021: R4 FGS diagnostics package (2026-09-12)

Final deliverable of the R4 diagnostics phase per
`fgs_opt_r4_diagnostics_handoff_20260912.md`. Audience: the code-optimization
agent. **No solver implementation was changed in this phase.**

Reviewed 2026-09-14: corrected sweep-threading interpretation, sample-share
limits, finalist comparisons, and optimization prerequisites. The handoff's
**2026-09-14 follow-up validation request** governs the next diagnostics pass.

**This file is the sole entry point for the next agent.** Before acting, also
read `~/.claude/CLAUDE.md`, repo `CLAUDE.md` (+ the policy files it routes to:
WORKFLOW/TESTING, and HPC.md before any cluster work), and on ORC
`/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md`. Everything
campaign-specific is in this directory:

- Pins, environments, launch commands, per-job results:
  `fgs_opt_v9_provenance_20260912.md`. Key pins — FLOWPanel exec `f03ab18`
  (tag `campaign/p021-cold-exec-20260912-v9`, worktree
  `/home/rander39/campaigns/p021-cold-opt-20260912-v9/FLOWPanel.jl`, env in
  `.../env`); FastMultipole `ef10643a` and FLOWVPM `05c658f7` (tag
  `campaign/p021-cold-exec-20260910-v1`, worktrees under
  `/home/rander39/campaigns/p021-cold-20260910-v1/`).
- Durable evidence: `fgs_opt_evidence_20260912/opt-<jobid>/` for jobs 13657404
  (all-stage + 20-rep profile), 13660643 (stage A: leaf axis),
  13661797/13661798 (stage B: P/MAC axes), 13663309/13663310 (finalist
  confirmation, independent nodes). Each has `harvest_summary.md` plus
  `screen_rank.csv` or `confirm_table.csv`; all tables audited against raw
  CSVs. Raw profile text (`cpu_flat.txt`, `cpu_tree.txt`, `allocations.txt`)
  lives under each job's `profile-fgs-j64-b1/R4/fgs_*/j64_b1/`.
- Ground rules that carry over: **any solver/harness code change requires a
  new pinned generation (v10+) built per the worktree workflow in the
  provenance file — never move v9 tags or edit a live/queued worktree.**
  Benchmarks rank by median prepared time to accepted accuracy; never loosen
  the gates in §6. Speedup claims need the matched alternating-batch
  comparison used in the confirmation jobs. No notebook entry without Ryan's
  approval.

Case: R4 = frozen `65_209` mesh, **58,192 panels**, Dirichlet, direct-assembled
RHS; prepared-scope timing (constructor excluded), BLAS=1, 64 pinned physical
cores (j64) on exclusive zen3 (m12), Julia 1.11.7, exec pin `f03ab18` (v9).

## 1. Configuration conclusion (measured)

**The provisional seed survived full bounded tuning: retained R4 configuration
= FGS P=8, MAC=0.4, leaf=100, inner=3, lexicographic, cached leaf LU, rlx=1,
calibrated tolerance 3.479128881193055e-7** (BC rel-L2 4.78e-7 ≤ 1e-6 gate).
Tolerance recalibration is deterministic: bit-identical values across 4 jobs
and 2 nodes.

Landscape (median prepared j64/b1 seconds, within-job comparisons):

| axis | values tried | result |
|---|---|---|
| inner (13657404) | 1,2,3,5,10 | 3 best (10.71); 5 close (11.00); 1 worst (13.05) |
| leaf (13660643) | 25,50,100,200 per inner∈{3,5} | 100 best both (11.16/11.28); 200 worst (~12.8) |
| P (13661797/8) | 6,8,10 | 8 best (11.01/11.13); 10 +8%/+6%; **6 FAILS certification** |
| MAC (13661797/8) | 0.3,0.4,0.5 | 0.4 best; 0.3 +36%/+38%; **0.5 FAILS certification** |

P=6 and MAC=0.5 fail calibration in both jobs ("FGS staircase has no certified
crossing with a decreasing successor") — the looser far field cannot certify
the 1e-6 BC gate at R4. Outer-iteration counts are nearly invariant along P/MAC at
fixed inner (28/28/27 for i3; 18/18/18 for i5): the P=10 and MAC=0.3 penalties
primarily reflect per-iteration cost; i3 has a one-iteration variation.

**Finalist confirmation (13663309 m12-1-17, 13663310 m12-1-25; ≥10 unprofiled
trials per batch, both configs sequential in one matched process per stage;
all 12 batches pass all gates):**

| batch | i3/l100 | i5/l100 |
|---|---:|---:|
| 09 j4 baseline | 16.897 | **14.653** |
| 09 j64 baseline | 11.402 | 11.357 |
| 09 j64 screen | **10.848** | 11.129 |
| 10 j4 baseline | 17.066 | **15.501** |
| 10 j64 baseline | 11.116 | 11.021 |
| 10 j64 screen | 11.099 | 11.082 |

At j64/b1 there is **no consistent winner observed**: matched median gaps
range from 0.15% to 2.59% and change direction between batches (within-batch
spread ~0.6 s). This is not a formal statistical equivalence result. At j4/b1
**inner=5 is decisively faster on both nodes: 13.3% (13663309) and 9.2%
(13663310)** (fewer outer iterations →
less per-outer overhead, which weighs more at low thread counts). i3/l100
stays retained as the established seed: it wins the original screen and job
13663309's confirmation screen, while i5 wins job 13663310's confirmation
screen (11.082 vs 11.099 s). i5/l100 (tolerance
2.5712316195808637e-7) is the recorded low-thread-count preference.

## 2. Measured R4 cost profile (retained config, 20-rep accumulated CPU profile, job 13657404; 6,717 main-task snapshots; corroborated by four 10-rep profiles in 13660643/13663309/13663310)

Reported attribution within main-task snapshots (self and inclusive frames
are distinguished below; these rows are not an additive wall-time budget):

| cost | share | evidence |
|---|---:|---|
| dense nonself products → BLAS dgemv (`dgemv_kernel_4x4` 5170 + 4x2/4x1 variants; `dgemv_n_ZEN` inclusive 5367) | **~77–80%** | cpu_flat.txt |
| `scatter_nonself_influence!` | ~6% | job-1 harvest (tree) |
| `solve_leaf!` (cached LU backsolves) | ~5% | job-1 harvest (tree) |
| `daxpy` (residual/update axpys) | ~2.4% | cpu_flat.txt |
| Base `setindex!`/`getindex`/broadcast machinery | ~3–5% | cpu_flat.txt |
| all `fmm!` passes (far field) | ~1.3% | job-1 harvest |

Structure per solve (i3): 27 outers = 81 inner sweeps + 28 FMM passes. i5:
17–18 outers = 85–90 sweeps. **Total sweep count is roughly conserved between
finalists** → per-sweep dense-product cost dominates and the j64 tie follows.
Dividing the j4 time gaps by the roughly ten fewer outer iterations gives
0.224 and 0.157 s per avoided outer on the two nodes. This is a heuristic,
not a measured stage cost: sweep counts differ, and scatter occurs per leaf
within each sweep, not once per outer. Stage timers are needed to separate
FMM, residual, sweep, and other costs.

Allocations (Allocs.@profile at 1% sampling + unprofiled trial): top site by
far is **`influence!` at FastMultipole `solve.jl:1336`** — array-slice
`getindex` (`similar` → `Array` per call). One prepared solve allocates
**642 MB** total; GC time is only 0.09 s (~0.5% — measured, small); retained
3.2 GB; peak RSS 4.5 GB. Same top site in all five profiles (both configs).

Thread scaling (measured, same config): j4/b1 17.0 s → j64/b1 11.1 s =
**1.53× from 16× more cores**. Idle worker threads show as futex waits;
this does not establish socket memory-bandwidth saturation. Source inspection
at FastMultipole pin `ef10643a` (`src/solve.jl`, `gs_sweep!`) shows that the
retained **lexicographic sweep is a serial leaf loop**. Its per-leaf GEMVs
also run serially with BLAS=1. Only the alternative `:colored` sweep uses
`Threads.@threads` across leaves. Workers may participate in other solve
stages; their activity must be measured separately.

The main-task profile establishes a dense-product hotspot, but does not by
itself establish total solve wall-time fractions. In particular, the ~1.3%
FMM attribution cannot account for worker execution and waits without a
thread-complete profile and stage wall timers.

## 3. Implementation opportunities (provisional priorities)

Measured = supported by the profiles above. Hypothesis = plausible mechanism,
needs the listed measurement before implementation.

1. **Dense nonself product engine (~77–80% of main-task samples, MEASURED).**
   This is the leading hotspot; its wall-time share remains to be measured.
   The matrix blocks are revisited every sweep (81 sweeps/solve); actual
   DRAM traffic depends on cache residency and must be measured.
   Sub-levers, all HYPOTHESES pending the bandwidth measurement in §4:
   a. **Reduced-precision influence blocks (Float32 storage, Float64
      accumulate)** — halves matrix-storage bytes, but conversion,
      accumulation, and kernel efficiency determine actual traffic and time.
      No total-solve speedup estimate is established. Benchmark the mixed-type
      kernel and re-certify the 1e-6 BC gate with recalibrated tolerance.
   b. **Cache-blocking / sweep fusion** — reuse a resident block across
      multiple RHS updates or fuse leaf blocks to raise arithmetic intensity;
      gemv is O(1) flops/byte. Cross-leaf or cross-sweep reuse must preserve
      dependencies: each leaf solve consumes preceding leaves' RHS updates.
      Otherwise this is an algorithm change requiring fresh convergence
      validation, not just a storage optimization.
   c. **Batching leaf gemvs** (grouped gemm-like kernels over leaves sharing
      source blocks) — attacks kernel-launch/loop overhead and enables
      blocking; single RHS per sweep limits naive gemm conversion.
2. **Thread-scaling recovery (MEASURED symptom: 16× cores → 1.53×).**
   The retained sweep is serial by construction. First measure stage wall
   time and active threads; distinguish serial execution from per-core
   memory limits, NUMA effects, and aggregate bandwidth saturation. Then
   assess dependency-preserving parallelism, the existing colored route,
   and block placement. Reduced precision is not the only possible route
   to gains. Fewer allocated cores are a cost hypothesis, not an established
   equal-time alternative.
3. **`scatter_nonself_influence!` (~6%, MEASURED).** Fuse the scatter into the
   gemv output epilogue or reorder targets so scattered writes become
   contiguous while preserving leaf-update order. Measure its wall-time share
   before assigning a total-solve gain bound.
4. **`solve_leaf!` (~5%, MEASURED).** Batch the per-leaf cached-LU triangular
   solves (dtrsv → blocked/batched form) only where leaf dependencies permit.
   Measure its wall-time share before assigning a gain bound.
5. **`influence!` slice temporaries at solve.jl:1336 (MEASURED site; total
   impact unmeasured).** 642 MB allocated per solve but GC is only 0.5%; the
   remaining cost is allocator traffic + bandwidth churn interleaved with
   lever 1. Preallocating views/buffers is cheap, low-risk, and de-noises
   every future profile. After baseline instrumentation, test an isolated
   view/buffer change; report allocations and matched wall time separately.
6. **Convergence-rate levers (sweep-count reduction) — HYPOTHESIS.** Sweeps
   are similar across inner∈{3,5}; that does not bound other algorithmic gains.
   Outers are nearly invariant along the tested P/MAC neighbors.
   Warm starts and relaxation tuning belong to other 021 threads;
   at j4-like thread counts prefer inner=5 (measured 9–13%) as a config win.
7. **FMM far field (~1.3% of main-task samples): lower initial priority.**
   Measure full-stage wall time before dismissing it. P8/MAC0.4 is best among
   the tested certified neighbors, not a proof of a global optimum or absence
   of implementation opportunities.

## 4. Required follow-up measurements

Execute the detailed test sequence and deliverables in the handoff's
**2026-09-14 follow-up validation request**. In order:

1. Confirm executed sweep order, BLAS threads, CPU placement, and per-stage
   thread activity against the pinned source and saved configs.
2. Add stage wall timers, with an instrumentation-equivalence control and a
   measured timing-overhead check. Separate FMM, influence mapping, residual,
   leaf solves, nonself products, scatter, and remaining solve work.
3. Record actual per-source-leaf GEMV dimensions, total matrix bytes, and
   interaction/scatter counts; these GEMVs aggregate target blocks rather
   than necessarily representing one leaf pair each.
4. Run the retained-config thread ladder j∈{1,4,8,16,32,64}, BLAS=1, with
   pinned physical CPUs and recorded socket/NUMA placement. Collect thread
   activity and bandwidth/cache counters in separate instrumented runs.
   Interpret scaling stage by stage, without assuming a bandwidth knee.
5. Use the evidence to select the next isolated experiment. A colored-sweep
   control is conditional on dependency/scheduling costs; a slice fix or
   Float32-storage pilot belongs to a subsequent implementation experiment,
   with new pins and unchanged accuracy gates. Neither is a prerequisite
   for completing this follow-up diagnostics pass.

## 5. Dev-branch delta since the campaign pin (checked 2026-09-14)

The campaign measured FastMultipole at pin `ef10643a`. Its development branch
(`flowpanel-20260817`, local checkout
`~/Dropbox/research/projects/FastMultipole`) has since changed `src/solve.jl`
only in **constructor scope** (BRAINSTORM 030 block-assembly hook, commits
`dd70fb19` + `fff72d29`): calloc-backed `Matrices` storage and
`unsafe_get_block_matrix` (plain-Matrix block wrappers; serial hook-vs-probe
2.90×). The historical profile remains evidence for the campaign pin, and
**the `influence!` slice op is textually unchanged**. Constructor time is
excluded, but changed allocation/first-touch placement can still affect
prepared-solve performance; remeasure before transferring timing claims to HEAD.
The allocating update `sorted_influences[range] .-= view(influences, index)`
sits at pin line 1336, now ~line 1355 at HEAD. Opportunity 5 therefore still
stands on the dev branch. Two facts the implementing agent should use:

- The 030 hook machinery (`assemble_influence_block!`, `_calloc_vector`,
  `unsafe_get_block_matrix`) is shipped and is the natural toolbox for
  buffer-reuse/layout work in this file — prefer extending it over inventing
  parallel plumbing.
- A Phase-3 attempt routing FGS matrix assembly through `influence!`
  projection (`417489d5`) was **reverted** (`00bd5f5d`); 030 Phases 3/3b were
  retired 2026-09-11 because FGS is projection-bound and a dipole far field
  was insufficient. Don't re-walk that route; the reverted commit is the
  reference for what was tried.

## 6. Gates and validity

Every accepted measurement in this package: solver converged, finite accepted
solutions, certified FMM evaluator authoritative (no direct fallback anywhere),
BC rel-L2 ≤ 1e-6 (worst 6.68e-7, stage-A leaf=25), repeat-solution delta = 0,
≥10 unprofiled repetitions, BLAS=1, pinned cores, profiling isolated from
unprofiled timing. Failed candidates (P6, MAC0.5 ×2 jobs) are certification
findings with outputs retained, not gate relaxations. Config↔directory
mappings audited via each run's `config.toml`/`requested_config.toml`.
Known benign quirks: job-root `selected.sha256` cites the CONFIG_FILE it was
given (job-1's selected.toml); the profile process's single unprofiled trial
(17.28 s) exceeds the 10-rep median (11.1 s) — first-solve/cold-cache effect,
use the rep medians for speed claims.

## 7. Measured thread-ladder update (2026-09-15/16, jobs 13694724 + 13733332)

Follow-up tests 1–4 of the 2026-09-14 request are complete (v15 generation,
FLOWPanel `39ec4e36` tag `campaign/p021-cold-source-20260915-v15`, FastMultipole
`87cbc846` tag `campaign/p021-r4-diag-source-20260914-v10`). Job 13694724
(exclusive zen3, 64 CPUs, BLAS=1, pinned physical cores) ran the retained
config at j∈{1,4,8,16,32,64}, 40 trials/arm in alternating
uninstrumented/instrumented batches of 10. Evidence + audit scripts:
`fgs_r4_followup_evidence_20260914/diag-v15-13694724/` (all 131 files SHA256
verified against remote; `analysis/{audit-full.txt,ladder_summary.md}`).
All gates in §6 held on every arm: 27 iterations / 81 inner sweeps / 28 FMM
passes invariant across arms, identical instrumented/uninstrumented histories,
BC rel-L2 4.78e-7, direct/FMM disagreement 5.5e-9 (j1 equivalence controls),
repeat delta 0. Instrumentation overhead −2.1%…+0.9%, inside batch spread —
instrumented stage numbers below are budget attribution, not speed claims.

### Measured wall-time budget (medians)

| Arm | Uninstr. wall s | Speedup | Efficiency | FMM total s | Nonself prod. s | Init s | Scatter s | Leaf solve s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| j1 | 37.56 | 1.00 | 1.00 | 26.73 | 8.16 | 1.00 | 0.79 | 0.60 |
| j4 | 16.94 | 2.22 | 0.55 | 6.88 | 7.79 | 0.65 | 0.79 | 0.54 |
| j16 | 12.20 | 3.08 | 0.19 | 2.09 | 8.00 | 0.32 | 0.80 | 0.55 |
| j64 | 10.96 | 3.43 | 0.05 | 1.11 | 7.89 | 0.31 | 0.87 | 0.55 |

Exclusive stage sums close the budget on every arm (unaccounted ≤0.013 s).

### Ranking revisions (measured wall, superseding §2/§3 sample shares)

- **The §2 "nonself product 77–80%" was a main-task-sample share, not wall.**
  Measured wall at j1: the outer FMM far-field stage dominates (26.73 s, 71%);
  nonself products are 8.16 s (22%). The main-task view undercounted FMM
  worker time exactly as §2's caveat warned.
- **FMM far field scales 24× to j64** and shrinks to 1.11 s (10% of j64 wall).
  Opportunity 7's "lower priority" is thereby confirmed at high thread counts
  and refuted at j1 — but j1 is not a production point.
- **The serial leaf-sweep chain is the measured j64 wall owner:** nonself
  products 7.89 s (72%, 1.04× speedup), scatter 0.87 s (7.9%, flat with a
  −8.2% j32→j64 regression), leaf solve 0.55 s (5%, flat) — together ≈85% of
  j64 wall. This is the measured basis for opportunities 1–4; the whole-solve
  3.43× ceiling is Amdahl on this chain.
- Initialization saturates by j32 (0.31 s, 2.8% at j64); not a lever.

### Census (j1; per-arm census files retained in every arm)

1,068 leaves; 2.86 GB Float64 influence-matrix bytes; 95,390 directed edges;
48,627 undirected conflicts. Hypothetical conflict-free coloring: 79 colors,
sizes 1–27 (median 16), verified no same-color conflicts. Executed sweep
confirmed serial lexicographic in the loaded source (test 1).

### Test-5 condition: SATISFIED — colored sweep is the justified next experiment

The stage budget supports the conditional colored-sweep experiment: the only
non-scaling stage cluster owns ~85% of j64 wall while everything else scales,
and the census provides a concrete 79-color schedule (median 16-way
parallelism per color; sync cost bounded by 79 colors × 81 sweeps ≈ 6.4k
color barriers per solve). It changes accumulation ordering, so
it must be run as a separately calibrated configuration passing all §6 gates,
ranked by total time to accepted accuracy (iteration count may move off 27).

### Counters + stage activity (2026-09-16, job 13733332, v20 generation)

The counters/activity run completed on the fourth attempt: **13711596 (v17)**
died at the smoke gate (docstring-above-`using` parse error), **13712587
(v18)** at `reinterpret(Cint, fd(io))` (`fd(::IOStream)` returns 64-bit `Int`
on the cluster's Julia 1.11.7 vs 32-bit `RawFD` on ≥1.12), and **13733217
(v19)** at ack parsing (cluster perf 5.14 EL9 writes each control-FIFO ack as
5 bytes `ack\n\0`, NUL terminator included — measured by `od`; the fixed
4-byte reader de-synced on the second command). Each failure is harvested at
`fgs_r4_followup_evidence_20260914/counters-v{17,18,19}-*-FAILED/` with its
postmortem in the next version's `v*-deployment/submission-provenance-*.md`.
Job **13733332** (v20, FLOWPanel `5c1123b` tag
`campaign/p021-r4-counters-source-20260916-v20`, FastMultipole `adb9967d` tag
`campaign/p021-r4-activity-source-20260915-v11`) passed every gate: smoke
PASS, all four controls PASS, both arms solved/finite, 27 iterations, BC
rel-L2 4.78e-7, repeat delta 0 with identical history, certified FMM
authoritative, BLAS=1, zero-reset counter scope. Evidence: 46/46 files
SHA256-verified at `fgs_r4_followup_evidence_20260914/counters-v20-13733332/`
(tables: `analysis/counters_summary.md`). Not performance trials
(perf-boundary + instrumentation overhead; baseline diagnostic walls 17.57 s
j4 / 11.92 s j64 vs §7 uninstrumented medians 16.94 / 10.96).

**Direct thread-activity measurement confirms the serial leaf-sweep chain.**
/proc-based per-stage activity (CLK_TCK=100, no incomplete endpoints, summed
over all 27 iterations):

| Arm | Stage | Σ span s | busy CPU s | avg active threads |
|---|---|---:|---:|---:|
| j4 | fmm (28 passes) | 6.90 | 26.98 | 3.91 |
| j4 | nearfield_update (27) | 9.46 | 9.44 | **1.00** |
| j64 | fmm (28 passes) | 1.10 | 34.04 | 30.96 |
| j64 | nearfield_update (27) | 9.59 | 9.63 | **1.00** |

The `nearfield_update` span (which contains the nonself-product / scatter /
leaf-solve chain) executes at almost exactly one active thread at BOTH j4 and
j64, and its span does not scale (9.46→9.59 s) — while `fmm` runs at
3.9/31 active threads and its span shrinks 6.3× (6.90→1.10 s, matching §7's
6.88→1.11 s wall medians). This is the direct activity-level proof of the §7
inference that ≈85% of j64 wall is a serially executed chain; the coarse span
totals reconcile with §7's per-stage sums (chain ≈9.1–9.3 s) within the
declared span-mixing granularity.

**Hardware counters (one warmed, FIFO-gated prepared solve; user-space; no
multiplexing):** j4 — 119.6e9 cycles, 429.9e9 instructions (IPC 3.59), 11.16e9
cache refs, 5.76% miss ratio, task-clock 37.6 s; j64 — 139.9e9 cycles,
456.7e9 instructions (IPC 3.26), 11.17e9 refs, 5.97% misses, task-clock
44.8 s. Cache-reference volume is essentially arm-invariant; the extra j64
cycles/task-clock are parallel overhead, not extra work.

**Bandwidth saturation remains unresolved — recorded as the pre-declared
limitation.** Generic cache counters cannot establish DRAM saturation
(uncore/IMC events are not process-attributable on this node), and the coarse
`nearfield_update` span mixes leaf/product/scatter, so no per-stage IPC can be
attributed. The whole-solve aggregate IPC of 3.3–3.6 argues against the solve
being memory-stall-dominated overall, but says nothing dispositive about the
serial chain's window. Per the follow-up plan, the colored-sweep experiment
proceeds on the serial-execution evidence, which the activity data now
establishes directly.
