# R4 FGS scalability diagnosis — consolidated plan (20260921c)

## Objective and strategy

Explain **both the 16–32-thread plateau and the 64-thread regression** of
`dagteam + f32full`, and identify improvements that reduce certified total
solve time. Use a gated experiment: reuse evidence, reproduce the effect,
locate the lost time with coarse measurements, then change one suspected
cause. Build detailed tracing or replay only when it will separate remaining
hypotheses. A useful workaround and an established mechanism are distinct
outcomes; report either without overstating the evidence.

Standalone synthesis of [plan A](fgs_scalability_diagnostic_plan_20260921.md)
and [plan B, revision 3](fgs_scalability_diagnostic_plan_20260921b.md).
All execution requirements are stated here; the earlier plans are provenance.
The sequence is: audit → reproduce and decompose → cheap controlled changes
→ only the detailed probes still needed → certify the best improvement.
Batch inexpensive comparisons within one allocation where practical; gate
engineering effort and diagnostic runs by the question they can answer.

Reported evidence to verify against raw records in Stage 0:

- The latest dagteam ladder is **33.07 / 5.97 / 4.50 / 4.42 / 6.42 s** at
  **1/8/16/32/64 threads**. See the
  [handoff](champion_adoption_reset_prompt_20260921.md).
- Earlier accepted trials give **4.475 s at 16 versus 6.223 s at 64**;
  interleaving improved the 16-thread result from **5.951 to 4.475 s**.
  See the [acceleration summary](fgs_acceleration_summary_20260919.md).
- The companion reports continued Krylov scaling on the same nodes, worse
  scaling for smaller/colored cases, and no FGS benefit from BLAS threading.
  These prioritize probes; different kernels and working sets mean Krylov
  scaling does **not** rule out FGS-specific memory saturation.
- Older lexicographic phase timings do not diagnose today's dagteam executor.
  Its shared queue, spinning workers, dependencies, and boundary reduction
  are candidate causes, not conclusions.

**Scope:** R4 only: 58,192 panels, P8/MAC0.4/leaf100/inner3, cached leaf LU,
champion precision, BLAS=1. Reuse existing R2/coloring evidence if useful;
no new R2, coloring, Krylov, or broad tuning campaign is required.

## Stage 0 — establish the baseline and prune unnecessary work

Reuse the existing cold/A/B harness and evidence before writing a new driver.
Locate the records associated with jobs 13777133/13829231 and the linked
handoffs; record missing fields as unknown rather than inferring them.
Extract per-run thread count, hardware, placement, settings/calibration,
iterations, actual sweep/FMM counts where available, solve time, allocations,
GC time, and certification. Check that compared rows represent the same
solver variant and timing boundary. Aggregated medians alone cannot recover
process-level dispersion.

Maintain a hypothesis/evidence table, distinguishing observed facts,
compatible explanations, and causes supported by controlled interventions:

| ID | Candidate explanation |
|---|---|
| H1 | Dependency width / load imbalance limits useful parallel work |
| H2 | Queue, scheduling, polling, or synchronization overhead grows with threads |
| H3 | Serial reduction, copies, or another outer stage limits total speedup |
| H4 | Memory throughput, coherence traffic, or NUMA locality limits kernels |
| H5 | Allocations, collection, or GC coordination delay the solve |
| H6 | Compute workers interfere with runtime/OS/GC work at high occupancy |
| H7 | Iteration counts, calibration, or other numerical work change with threads |
| H8 | SMT siblings / CPU placement change the effective compute resources |
| H9 | Frequency, power, or thermal behavior changes the operating point |

Placement is both a potential cause and an intervention on H4/H6/H8.
Resolve champion CPU-set semantics against the diagnostic node topology now;
flag any sibling sharing before interpreting a nominal 64-thread result.
Flat iteration counts demote iteration growth but do not prove identical
work or certified accuracy. Flat allocation bytes do not rule out GC pauses,
safepoint delays, or collections triggered by earlier allocations.

Inspect only the current executor and harness paths needed to answer:

- What does the timed solve include? Are reset, team lifetime, tree updates,
  stopping checks, and FMM calls inside it? Is 27 outer × 3 inner truly the
  intended fixed workload, including initialization/finalization calls?
- Does a near-field worker-cap knob already exist? If added, can excess
  workers be parked outside the queue and barriers without deadlock?
- Which dependency edges, upper tasks, and sweep barriers constrain execution?
  Can coarse phase timers or GC statistics already be reused?

If the actual R4 DAG is readily available, summarize task count, block sizes,
level sizes, and dependencies. **Do not build an expensive DAG analyzer yet.**
Level width is not the dynamic ready-queue width; block size alone is not a
validated task-time model. Defer weighted bounds to Stage 3 unless measured
weights already exist. Do not compare a sweep critical path directly with
whole-solve time.

**Exit:** an exact baseline manifest, a list of reusable measurements/knobs,
and the smallest experiment that can discriminate the leading hypotheses.
Record the exact harness entry point, diagnostic options, and row locations
so the execution handoff does not require rediscovery. No hypothesis is
declared causal from static inspection alone.

## Measurement contract (all stages)

**Platform.** Start on one exclusive zen3 node where the precision is
certified. Resolve the champion's socket/NUMA IDs against actual topology;
do not assume NUMA nodes 0–3 always mean the same socket or that 64 threads
occupy every physical core. Use pinned physical cores with explicit nested
CPU sets, balanced across the champion NUMA domains where possible. Keep
memory policy fixed across the initial ladder. At j=1 document the unavoidable
asymmetry. Record host/CPU model, sockets/NUMA/SMT, allowed CPU and memory
sets, actual worker affinity, GC/runtime thread counts, Julia/BLAS versions,
clock observations, and loaded package pins/paths. Verify page distribution
after setup/warmup, not just the requested `numactl` policy. If the historical
configuration used SMT siblings, preserve a separately labeled reproduction
of that configuration and compare with physical-core placement where
feasible. A comparison that also crosses sockets or memory domains tests a
combined placement change, not SMT alone.

**Frequency.** Record effective frequency during representative timed regions
when accessible with low overhead (e.g. supported APERF/MPERF counters).
Use separate matched diagnostic runs if measurement perturbs timing; missing
counter access does not block the campaign. Raw elapsed time remains the
primary metric. Do not subtract a “DVFS fraction” or require whole-solve
frequency normalization: memory stalls and waiting need not scale inversely
with core frequency. A normalized estimate is optional sensitivity analysis
only after validating that model for the relevant phase. If frequency is
suspected to explain material loss, use a permitted matched frequency control
when available, or leave the causal share uncertain. Lower spinning that also
raises useful-work frequency is a real intervention benefit.

**Fixed work.** Identical geometry, matrices, ordering, initial state,
precision, and numerical settings; exactly **27 outer updates × 3 inner
sweeps**, after confirming the harness semantics. Disable convergence stopping
only in diagnostic mode. Assert and record actual FMM calls and sweeps;
record residual histories and reject nonfinite or mismatched-work runs.
Reset all mutable solver state between solves, including queue/dependency
state and buffers; keep immutable matrices and cached LU prepared.

**Accepted accuracy.** Keep this as a separate arm using the existing
per-thread calibration and independent **BC relative-L2 ≤ 1e-6** certification.
Record the calibration, work counts, and residual history. If calibration
changes workload, call this a production comparison, not pure thread scaling.
A baseline/variant pair uses identical settings unless the intervention
explicitly changes the algorithm; any retuning is a separately labeled result.

**Timing.** Compile and warm up every exercised path. Measure prepared,
zero-reset solves with the same timing boundary as the baseline. Report
construction, warmup, reset, and validation/output separately; if the existing
metric includes reset, retain that metric and also report the decomposition.
Never quietly move required work outside the timed boundary to claim a gain.
Keep normal diagnostic output outside solve timing. GC remains enabled by
default; do not force a collection before every solve unless testing that
separate condition. Preserve the sequence of individual timings/GC events.

**Replication.** Begin with three independent fresh-process blocks per
configuration and five warmed solves per process; extend to ten if needed
for stability or intermittent GC. Randomize thread-count order within each
block. For interventions, run baseline and variant in fresh processes close
in time, alternating/randomizing their order. Analyze process medians and
paired ratios; inner solves are repeated observations, not independent
replicates. Report all process medians and dispersion; retain outliers with
GC/clock metadata rather than dropping inconvenient runs. Three blocks are a
screen, not grounds for a narrow CI. Use a separate, fixed-size confirmation
batch for selected winners (normally at least five fresh baseline/variant
pairs), chosen before seeing its results; report the CI method/sample count.
Do not repeatedly add runs until significance appears. If the budget leaves
an effect unresolved, report its interval and stop that comparison. This supports same-node claims;
a second matched node is only needed to establish portability.

Predeclare a **5% total-time gain** as the default practical screening
threshold. Require confirmation beyond observed variability; a useful phase
improvement below this threshold can still explain a small part of the loss.
Do not rerun clear failures simply to complete a matrix.

## Stage 1 — reproduce and locate the loss cheaply

Prepare one allocation with a baseline block, coarse diagnostics, and a
short ordered list of optional comparisons. Reuse existing harness switches;
implement a small conditional runner only if it is less work than a fixed
short list. A logical stage boundary does not require another submission or
human checkpoint. Estimate startup, construction, memory, and timed solves
from a pilot; reserve enough walltime for the highest-value follow-ups.
Checkpoint each completed configuration. Stop on numerical failure and do
not edit code while the allocation is using it.

1. Run the uninstrumented fixed-work ladder at **16/32/64**, plus **1 thread**
   as the work/speedup anchor. Add 8 only to resolve the knee or confirm the
   historical full ladder. Run the baseline accepted-accuracy arm at
   16/32/64 using the same replication screen. Existing matched, complete
   evidence may replace redundant accepted-accuracy screening, but not final
   certification of a winning change.
2. Collect per-solve allocation and GC statistics where supported. In separate
   matched runs, add only coarse, exclusive wall timers for initialization,
   FMM, influence mapping, residual, near-field sweeps, and final update.
   Check their sum against total elapsed time and expose any remainder.
   Within sweeps, time team lifecycle, boundary reduction, and copies if these
   are naturally coarse regions. No per-poll timers or full traces yet.
3. Confirm the instruments preserve numerical results and the scaling shape.
   Target ≤5% total overhead at each rung, with overhead well below the effect
   being attributed. If this fails, reduce sampling or use separate probes;
   do not subtract a guessed constant overhead.

Quantify the plateau and regression separately with phase differences
`T_phase(32) - T_phase(16)` and `T_phase(64) - T_phase(32)`; include negative
contributions. Also report the shortfall from ideal doubling,
`T_phase(32) - T_phase(16)/2`: a flat phase can explain the plateau even
though its absolute difference is zero. Ideal doubling is a reference, not
an assumed achievable target. Reconcile exclusive phases per solve and use
matched sample means for additive phase differences; separately computed
phase medians need not sum to the median total. For a serial phase, report
the optimistic whole-solve saving from eliminating it; this prevents
optimizing a visibly slow but irrelevant region.

**Within-allocation follow-ups:** after the baseline is valid, include one
64-thread placement A/B by default because prior evidence makes it a strong
candidate, provided setup is already available. Include a worker-cap A/B if
a validated knob already exists. Add 48 threads if the onset remains useful
to locate; add 60 or 62 only if a near-full-occupancy comparison could change
the operating recommendation. A topology flag triggers the feasible SMT/
physical-core comparison before attributing the loss elsewhere. GC settings
are conditional on measured GC symptoms. Co-run tests remain optional and
late: they offer less specific evidence than these direct interventions.
If automatic selection would require substantial new machinery, preselect
this short list from Stage 0; do not build an orchestration framework.

**Gate:** if either effect fails to reproduce, report exactly which effect
is absent. Compare original versus current hardware, placement, pins,
calibration, and timing contract before changing the solver. Failure to
reproduce on a different/exclusive node does **not** establish co-tenancy or
OS interference. Continue only on a reproducible effect; otherwise return a
bounded reproduction report with the smallest missing comparison.

## Stage 2 — targeted, inexpensive interventions

Select from the table based on Stage 1; do not run a Cartesian product.
Reuse any comparisons already completed in Stage 1. Start at **64 threads
with a 16-thread control** and confirm a promising change at 32. Prioritize
an existing worker-cap knob and placement comparison; then build the smallest
change aimed at the phase that actually owns the loss. Each variant must
pass the numerical checks below before its speed is considered.

| Probe and control | What it can establish; next action |
|---|---|
| Cap near-field workers at 16 while retaining 64 for FMM; keep those 16 on the matching nested CPU set. Inactive workers must not poll or join active-team barriers. Account for parking/wakeup costs. | Recovery localizes a loss to sweep participation, including its synchronization/locality effects. It does not alone distinguish narrow DAG width, queue contention, or memory traffic. If useful, try cap=32 only to select the better operating point. |
| If waiting/queue overhead is implicated, replace busy polling with one bounded backoff policy, preserving dependency visibility, progress, and arithmetic ordering. | Lower wait/CPU cost plus recovered wall time supports a waiting-policy cost. No benefit does not rule out lock contention while managing nonempty queues. Check for starvation/deadlock and lost low-thread performance. |
| If placement is suspect, compare champion interleave with one explicit alternative at fixed j and CPU affinity; construct/first-touch the state afresh under each policy. | Verified page-placement change plus recovery establishes an operational locality benefit. Prefer parallel first-touch/owner-local allocation if available; serial first-touch under “local allocation” may concentrate pages and is not a clean locality test. Confirm the winner at 16. |
| If GC pauses/allocations are material, compare one supported GC-thread setting or targeted allocation removal with baseline, holding compute threads fixed. | Reduced GC time and matched solve-time recovery support GC contribution. A bounded GC-disabled solve is diagnostic only: check memory headroom, collect outside timing, restore GC in `finally`, and confirm the eventual fix with GC enabled. |
| If full-core occupancy/runtime interference is suspected, add 48 or 60 threads; optionally 62 to bracket onset, with explicit affinity. | Recovery establishes a useful thread cap, not its mechanism: fewer workers also change contention and topology. Isolating headroom requires fixed compute-worker count/placement and a verified change to auxiliary-thread placement or available spare cores. If that control is infeasible, leave the mechanism unresolved. |
| If coarse timing shows substantial boundary reduction/copy cost, parallelize across disjoint target leaves, retaining ascending source order within each target. | Recovery in the predicted phase and total supports a serial bottleneck. Implement only when its recoverable time justifies the effort. |

A thread cap may cure the 64-thread regression while leaving the plateau
unexplained. Keep those conclusions separate. If fixed-work solves scale
but accepted-accuracy solves do not, follow H7 first: compare calibration,
work counts, and residual histories, then test a matched-setting accepted
solve or one controlled calibration correction. If another outer stage owns
the loss, investigate that stage instead of adding sweep instrumentation.

**Gate:** stop adding diagnostics once measured phase loss and a controlled
intervention explain the relevant effect within uncertainty. Quantify
individual contributions first; then test only combinations likely to help.
If the remaining plateau/loss is material and ambiguous, choose the specific
Stage 3 probe that can resolve it. Do not claim that an intervention's total
gain is a unique causal percentage when it also changes other phases. Report
those interactions, and test combinations only after individual comparisons.
Follow the applicable execution and submission policies when implementing.

## Stage 3 — resolve only the remaining ambiguity

### Dependency availability versus executor overhead

Start with preallocated, padded per-worker aggregates: completed task counts,
useful work, queue management/lock wait, empty-queue wait, and boundary wait.
Avoid shared instrumentation counters and timing every spin. Add a sampled
timeline only if aggregates cannot distinguish lack of ready work from
failure to dispatch it. Capture task type, worker, start/end, and dependency
readiness; label an interval unknown if readiness was not observed. Sample
early/middle/late sweeps and validate observer overhead/scaling again.

Report useful occupancy, imbalance, and queue availability. Never sum
overlapping worker durations as elapsed time. Include both lower solves and
upper products plus barriers in the execution model.

For measured task costs, define total useful work $W$ and weighted longest
DAG path $L$. With fixed task costs and $S_j = W/T_j$, an optimistic
sweep bound at $j$ workers is

$$
T_j \geq \max(W/j, L), \qquad S_j \leq \min(j, W/L).
$$

State which costs/serial sections the model includes. Check **$W$ against
one-worker useful sweep time**, and compare $L$ with the corresponding
multithreaded sweep, not the full 33 s solve. Costs can change with cache,
NUMA, and contention; use representative measured weights and sensitivity
bounds. Static level sizes and a greedy schedule are supporting diagnostics,
not proofs of actual occupancy or optimal scheduling.

An empty ready queue has two different interpretations: dependencies can
starve otherwise idle workers, or all useful tasks can already be running
while excess workers spin. Use active-task counts and completion timing to
separate these before attributing empty-queue time to H1. This also prevents
counting harmless idle worker-time as recoverable elapsed time.

If ready work exists while workers wait, prioritize queue/scheduling changes
(e.g. batched queue access or ownership changes), one at a time. If little
work is ready and time approaches a credible dependency bound, test the
bounded dependency-relaxing pilot below; static width alone is insufficient.

### Kernel throughput versus memory/locality limits

Build **frozen-input replay** only if it will decide between dependency/
executor limitations and kernel throughput. Replay actual lower/upper block
products and cached leaf solves with frozen inputs, independent outputs, and
balanced static scheduling. Include gather/copy costs, the full representative
working set, and enough sweeps to match production reuse. Check outputs are
consumed. Match precision, layout, affinity, page policy, and warmup; document
traffic/order/cache differences. This is an empirical throughput reference,
not a valid FGS solve or a mathematically guaranteed ceiling.

When available, collect phase-scoped memory-controller bandwidth, relevant
cache/coherence events, and CPU time in separate runs. Establish that the
replay itself reaches a stable memory-throughput limit before calling it a
bandwidth ceiling. Approaching replay throughput alone cannot distinguish
memory from compute limits. Generic cache misses, high CPU utilization, or
calculated “useful GB/s” do not establish bandwidth saturation.

An optional co-run probe uses two independently prepared 32-thread processes,
solo and simultaneous on identical disjoint CPU sets, with documented page
placement and verified shared memory-controller domains. Use it only if cheap
and relevant. Slowdown shows shared-resource interference (possibly memory,
cache, power, or frequency); no slowdown on separate domains says little
about one 64-thread process. **It does not replace bandwidth evidence.** If
counters are unavailable, report locality/traffic interventions and leave the
specific bandwidth-saturation claim unresolved where necessary.

## Effort and production tradeoffs

Keep diagnosis and optimization decisions separate. Each further diagnostic
must address a named material ambiguity; each larger implementation should
have a measured upper estimate of recoverable total time and a plausible
path to enough benefit to justify its effort. Finishing the FGS diagnosis
does not require it to beat a different solver.

Plan B reports `krylov_ilu_nfcache` at **4.25/3.56/2.73/2.41 s** for
**8/16/32/64 threads**, excluding cache build and using about **8.5 GB** of
cache. Verify these records before citing them. Use them as context, not a
performance bound or an automatic stop condition. A production comparison
must match hardware, accuracy, timing boundaries, and anticipated reuse;
include setup/cache rebuild amortization and peak memory. No new Krylov
campaign is required. The report may recommend FGS for memory reasons or
limit further optimization effort while still reporting an unresolved
mechanism honestly.

## Algorithm follow-up, only with demonstrated headroom

If dependency relaxation has substantial potential, pilot block-Jacobi while
retaining the FMM outer loop, leaf partition, matrices, cached leaf solves,
and precision. Read previous-sweep strengths and write into a separate
buffer. Screen relaxation **1.0/0.8/0.6/0.4**, with the existing **300-iteration
cap**. Reject divergence/nonfinite results early. Screen at 16/64; extend
only promising settings to 32 and independent confirmation.

Compare **certified total solve time**, including buffering and additional
iterations, against the best validated dagteam configuration. A faster sweep
or higher parallel efficiency that loses in total time is not an improvement.
No ILU–Krylov replacement is in scope. If replay/bounds show little headroom,
retain the best worker/placement policy and avoid an algorithm rewrite.

## Validation, effort limits, and handoff

- Keep instrumentation and experimental variants optional; default solver
  API/behavior remains unchanged. Test semantics-preserving variants against
  baseline residual history and the existing **1e-8 repeat-solution gate**,
  using the established norm/normalization. Independently certify **BC
  relative-L2 ≤ 1e-6**; algorithm-changing variants need not reproduce the
  baseline history. Validate the final recommendation with GC enabled and
  uninstrumented accepted-accuracy timings at 16/32/64.
- Before execution, read `agent_policies/TESTING.md` and select the narrow
  FGS/solver checks plus a local smoke at **at most four threads**. Include
  multi-thread progress checks for worker caps/backoff; local smoke cannot
  establish 64-thread correctness or performance. High-thread confirmation
  stays on HPC.
- Follow current HPC/site submission policies. Campaigns use clean worktrees
  pinned by annotated tags for FLOWPanel and every development dependency,
  including FastMultipole/FLOWVPM when loaded. Verify Manifest paths and
  record actual loaded tags/SHAs before submission. Outputs go to the
  consolidated data root. Keep source fixed while jobs use it.
- Persist one row per solve with configuration/process/replicate IDs, pins,
  thread/worker counts, affinity and placement, numerical settings, work
  counts, elapsed/phase times, allocations/GC, available frequency data, and
  numerical-gate results. Record unavailable fields explicitly. Keep a small
  analysis script with the report so phase accounting and paired effects
  can be reproduced without manually scraping logs.
- Stage 0 should remain a small evidence audit. Budget Stage 1 from measured
  startup/preparation cost plus solves × repetitions; do not assume “one
  short job” or a day without a pilot. Checkpoint completed rows. Spend more
  runs on ambiguous or winning comparisons, not failed branches. Stage 3
  requires a named unresolved question and a probe that could change the
  recommendation.

Deliver one compact report with provenance/raw-row locations; baseline and
accepted-accuracy ladders; exclusive phase times and their contributions to
each loss; paired intervention effects/uncertainty; numerical gates; and a
recommended operating point or code change. Include DAG analysis, timelines,
replay, and counters only if used. Classify each hypothesis as supported,
demoted, or unresolved and explain the basis.

**Success:** both scaling effects have quantitative explanations supported
by controlled comparisons, and an improvement wins on independently
confirmed, certified total time. State how much of each effect is explained
and the uncertainty; exact recovery is not required when a measured resource
or dependency limit accounts for the remaining loss.

**Bounded partial outcome:** a useful workaround, an identified limitation
without a winning fix, or an unresolved effect is still reportable. State
which objective remains unfinished, why the next probe was deferred, and
the smallest experiment that would discriminate the surviving explanations.
A faster alternative solver or an exhausted run budget is not proof that
FGS cannot improve. Do not manufacture a root cause to satisfy the stopping
rule.
