# 033 — FGS scalability and setup-cost program

Opened 2026-09-26 (Ryan directive, promoting
[`021/fgs_cooperative_leaf_scaling_plan_20260926.md`](021_rotor_hover_solver_benchmarks/fgs_cooperative_leaf_scaling_plan_20260926.md)
to a standalone item). Goal: the most performant FGS product we can produce —
high performance first, strong robustness/correctness, lean engineering
effort where possible.

## RESET BRIEF

- **Status:** Track A theory COMPLETE 2026-09-26 (A-T1 + A-T2 + A-T3
  technical-complete and clear-context approved). Trials A-R1 AND A-R2 both
  COMPLETE 2026-09-26 session 2 (technical + clear-context approved): column
  layout selected, prototype coop executor implemented and validated, local
  screening done. A-R3 (HPC j-ladder A/B) not started. No live jobs. Nothing
  here authorizes a campaign; A-R3 requires pinned worktrees (checkouts
  currently dirty with uncommitted A-R2 changes).
- **Session-1 headline (A-T2):** split-graph model reproduces gate-0 exactly
  (289.462 MB unsplit); w=4 selective sweep ceilings 3.21×/2.53×/2.10× at
  elapsed overhead 0/8/16 µs (1 µs ≡ 11.72 KB, gate-0 calibration) ⇒
  projected solve ~1.86/2.03/2.19 s vs 3.24 s (1.74×/1.60×/1.48×) — Tier-1
  >15% margin survives to c≈32 µs; beating cold ILU 2.41 s needs c ≲ 16–32 µs.
  Critical path stays span-bound (breadth never binds at j64); LU share rises
  8%→26% at w=4; unconditional `all` splitting goes sub-1.0× at c=64 µs, so
  A-T3 measured-cost selection is mandatory. Tasks/sweep ≤4.3k (dagedge's
  fatal regime was 31–49k). Evidence: `033_atheory_20260926/`.
- **What this is:** four method tracks to attack (i) FGS solve-phase thread
  scaling (dagteam+backoff gains only 1.07× from 32→64 threads at R4) and
  (ii) FGS setup cost (~210 s more than krylov_ilu_nfcache at R4). Each track
  carries its own Theory → Trial → Implementation sub-items with per-sub-item
  gates: `[ ] technical completion` / `[ ] clear-context review` (INDEX
  semantics; conditional approval not allowed).
- **Tracks:** A cooperative leaf products (priority 1) · B setup/construction
  cost (priority 1, independent) · C exact triangular solver via subdomains
  (SPIKE-like) · D overlapping local GS + coarse correction.
  **Ryan ruling 2026-09-26: C and D are unconditional** — attempted regardless
  of A/B outcomes; their stop gates fire only on genuine infeasibility or
  measured failure, never on "A was good enough." **But both are RYAN-GATED:
  converse with Ryan on the theory and implementation plan before starting any
  C or D sub-item** (see the gate notes under each track header).
- **Acceptance (pinned, two-tier):** Tier 1 per-track = beat dagteam+backoff
  in paired, uninstrumented cold solves by the track's stated threshold.
  Tier 2 item headline = explicit final best-FGS vs krylov_ilu_nfcache
  assessment (cold and cumulative-with-setup), recorded whatever the answer.
- **Standing rulings:** no GMRES/Krylov outer iteration (Ryan); ≤64
  shared-memory threads, distributed memory out of scope; "cold" =
  zero-initial-guess solves, arms batched per (threads, placement) process;
  BLAS single-threaded in experiments; local runs ≤4 threads; production
  defaults unchanged until a candidate wins paired solves; official campaigns
  from tagged worktrees per `agent_policies/HPC.md`.
- **Next actions:** A-R1 technical work DONE 2026-09-26 (column layout wins
  at every shape, row rejected short-wide; local barrier ~0.1 µs, effective
  h_w ≈ 0.4–3 µs median leaves, ~24 µs at the 1,450 tail —
  see [`033_ar1_microbench_results_20260926.md`](033_ar1_microbench_results_20260926.md));
  A-R2 prototype DONE + reviewed 2026-09-26 (coop executor in FastMultipole
  behind `dagteam_coop`, R1 validation 14/14, local R4 j=4 no-gain as
  expected — see [`033_ar2_prototype_20260926.md`](033_ar2_prototype_20260926.md)).
  B-T1 setup profiling DONE + reviewed 2026-09-26: 95.2% of FGS setup is the
  single-threaded `nonself_influence_matrices` per-column probe; the ~210 s
  gap vs ILU is attributable to it; lever = route nonself+self through the
  030 `assemble_influence_block!` hook, parallel over disjoint blocks (see
  [`033_bt1_setup_profile_20260926.md`](033_bt1_setup_profile_20260926.md)).
  B-R1/B-R2 DONE + reviewed 2026-09-26: hook adoption rejected (semantic
  mismatch), fix = threaded existing probe, bitwise-certified, 3.59× at j4,
  projected j64 setup ≈10–20 s (gap vs ILU inverts); B-G gate passes (see
  [`033_br12_machinery_and_fix_20260926.md`](033_br12_machinery_and_fix_20260926.md)).
  Next (both HPC-bound, need commits + pinned worktrees): **B-I1** (land
  threaded setup, j64 re-measure, feed Tier-2 arithmetic) and **A-R3 paired
  j-ladder A/B** (w∈{1,2,4} at fixed total worker budget, uninstrumented
  paired data for the A-G gate).
  Ryan decisions pending: accept back-product cooperation scope extension or
  request model-pure variant before A-R3; nw%teamw==0 rule at A-I1; accept
  B-T1's R4+R1 substitution for "R4 and one larger mesh"; `threaded_setup`
  production default after B-I1; run full test suite before committing. Tracks C
  and D wait on their Ryan conversations (Track C freshly motivated: serial
  LU exceeds the w=4 cooperative product on the 1,450 tail leaf).

## Motivation and provenance

This item executes the recommendations of the 021 scalability review,
[`fgs_cooperative_leaf_scaling_plan_20260926.md`](021_rotor_hover_solver_benchmarks/fgs_cooperative_leaf_scaling_plan_20260926.md)
(read it first — it carries the quantitative record, caveats, and links).
Companion handoffs, also under `021_rotor_hover_solver_benchmarks/`:

- [`fgs_exact_triangular_solver_plan_20260926.md`](021_rotor_hover_solver_benchmarks/fgs_exact_triangular_solver_plan_20260926.md) — Track C build order.
- [`fgs_overlapping_gs_coarse_correction_strategy_20260926.md`](021_rotor_hover_solver_benchmarks/fgs_overlapping_gs_coarse_correction_strategy_20260926.md) — Track D principles.
- [`fgs_overlapping_gs_coarse_implementation_handoff_20260926.md`](021_rotor_hover_solver_benchmarks/fgs_overlapping_gs_coarse_implementation_handoff_20260926.md) — Track D prototype plan.

Evidence base (all in `021_rotor_hover_solver_benchmarks/`):
[Stage 2 results](021_rotor_hover_solver_benchmarks/fgs_scalability_stage2_results_20260924.md)
(dagteam+backoff champion; j64 3.262 s; 32→64 only 1.07×),
[DAG analysis / gate 0](021_rotor_hover_solver_benchmarks/fgs_lshortening_gate0_20260924.md)
(lower GEMVs = 92% of modeled critical-path bytes; work/span 2.8),
[dagedge results](021_rotor_hover_solver_benchmarks/fgs_dagedge_benchmark_results_20260924.md)
(fine-grained edge tasks LOST 3.766 vs 3.540 s — granularity, not concept, was the failure),
[warm-start R4 results](021_rotor_hover_solver_benchmarks/fgs_warmstart_r4_results_20260925.md)
(fgs_proj2 per-step ties ilu_nfcache_proj2, ~210 s more setup ⇒ no crossover).

**Competitive record, stated plainly:** FGS currently trails krylov_ilu_nfcache
at every measured operating point (cold R4 j64: 2.41 vs 3.24 s; warm per-step:
statistical tie with ~210 s extra setup). The no-Krylov constraint is Ryan's;
this item's Tier-2 headline must report the ILU gap honestly either way.

## Acceptance criteria (pinned — do not re-litigate per review)

- **Tier 1 (per-track gate):** the track's candidate beats **dagteam+backoff**
  (current FGS-internal champion, FLOWPanel/FastMultipole defaults) in paired,
  uninstrumented, cold (zero-initial-guess) solves at R4 j64 by the track's
  threshold, arms batched per (threads, placement) process. Thresholds:
  Track A team-of-4 > ~15% wall clock; Tracks C/D thresholds set in their
  theory phases before any prototype benchmark, and must at minimum beat
  dagteam+backoff outright at matched accuracy.
- **Tier 2 (item headline):** on close, an explicit assessment of best-FGS vs
  krylov_ilu_nfcache — cold solve and cumulative-with-setup — recorded as the
  item's INDEX outcome whether favorable or not. Track B's setup savings feed
  this arithmetic directly.
- Ranking metric everywhere: **time to independently certified accuracy**,
  iteration count included; setup cost and memory reported separately.
- Correctness bar: zero and nonzero starts, 1 and multiple inner sweeps,
  repeatability, and mathematical agreement with the existing recurrence (or,
  where FP grouping changes — e.g. column-partition reductions — iterate
  equivalence in exact arithmetic plus accuracy certification).

## Standing rules

- No GMRES/Krylov outer iteration. Target ≤64 shared-memory threads; larger
  problems of interest; distributed memory out of scope.
- "Cold means cold-start" = zero-initial-guess solves, NOT fresh process per
  arm (Ryan 2026-09-23); batch arms per (threads, placement) process.
- BLAS single-threaded in all experiments. Local checks ≤4 threads;
  1/8/16/32/64-thread scaling on HPC only.
- Production defaults unchanged until a candidate wins paired, uninstrumented
  solves. No public API change for experimental prototypes (runtime toggles in
  the style of `dagteam_idle`/`dagteam_workers` instead).
- Official campaigns: pinned annotated-tag worktrees, provenance file before
  submission, outputs to the consolidated data root (`agent_policies/HPC.md`).

## Shared infrastructure (build once, in Track A trial; reused by all tracks)

- **S1 — Certification + paired-benchmark harness.** Accuracy certification
  (independent residual/accuracy check), paired same-process A/B driver,
  time-to-certified-accuracy ranking, setup/memory side channels. Extends the
  Stage-2 benchmark tooling rather than replacing it.
  - [ ] technical completion
  - [ ] clear-context review

## Track A — Cooperative leaf products (priority 1)

Small persistent teams (2–4 workers) cooperate on expensive leaf GEMVs inside
the existing dagteam DAG; GS dependencies retained. Theory phase is thin by
design — resist over-deriving. Source: plan doc §"First recommendation".

### A Theory

- **A-T1 — Correctness contract.** Record officially: iterate equivalence for
  output-row partitions (disjoint writes, bitwise-reproducible); for
  contiguous-column partitions with private partials + fixed-order reduction,
  exact-arithmetic equivalence + deterministic reduction order + accuracy
  certification as the contract. Zero/nonzero starts, sweep counts.
  - [x] technical completion (2026-09-26 — see "A-T1 — Correctness contract" below)
  - [x] clear-context review (2026-09-26 — see review record below)
- **A-T2 — Split-graph critical-path model.** Extend
  `benchmark/fgs_dag_L_profile.jl` (~50 lines) to recompute the critical path
  over the whole split graph (node GEMV cost ÷ w subtasks) with a per-subtask
  overhead parameter swept over a plausible range. Modeled ceilings from
  gate 0 (w=2 ≤1.85×, w=4 ≤3.21× sweep) are conditional estimates, not bounds.
  - [x] technical completion (2026-09-26 — see "A-T2 — Split-graph critical-path model" below)
  - [x] clear-context review (2026-09-26 — see review record below)
- **A-T3 — Split-size threshold rule.** Under ideal product scaling, split
  only when product cost > h_w/(1−1/w), with h_w the cooperative
  implementation's own measured elapsed overhead for the whole team.
  Confirm the complete cooperative product is faster in measurements.
  Do NOT seed h_w from dagedge's 8.2 µs (that
  figure includes useful computation). Which leaves split is an empirical
  outcome; do not restrict to a handful of hot leaves (top-1% splitting barely
  moved the bound).
  - [x] technical completion (2026-09-26 — see "A-T3 — Split-size threshold rule" below)
  - [x] clear-context review (2026-09-26 — see review record below)

### A Trial

- **A-R1 — Standalone kernel microbenchmarks.** Row vs column layouts at
  w=2/4 on representative leaf shapes (median 54-unknown short-wide; 1,450
  tail). These measurements pick or gate the executor variants — do not assume
  row partitioning wins.
  - [x] technical completion (2026-09-26 —
    [`033_ar1_microbench_results_20260926.md`](033_ar1_microbench_results_20260926.md):
    column layout wins everywhere, row rejected short-wide; effective h_w
    ≈ 0.4–3 µs on median leaves, barrier ~0.1 µs)
  - [x] clear-context review (2026-09-26, fresh-context AI audit: code +
    numbers recomputed clean, OpenBLAS reset mechanism reproduced live;
    verdict fixed-minor — 9 wording/consistency fixes, no rerun needed)
- **A-R2 — Prototype cooperative executor.** Persistent teams reused across
  the inner-sweep block (no per-leaf task spawning — the dagedge granularity
  failure is the anti-pattern); one worker does the cached LU solve and
  publishes successors after the product completes; team size a runtime
  toggle (`dagteam_workers`-style) for paired same-process arms; measured-cost
  split selection per A-T3. NUMA placement deferred (Stage 2 showed lock
  contention, not bandwidth).
  - [x] technical completion (2026-09-26 session 2 —
    [`033_ar2_prototype_20260926.md`](033_ar2_prototype_20260926.md):
    column-layout coop executor behind `dagteam_coop`/`teamw[]` toggle,
    default bit-identical to production; R1 validation 14/14 PASS; in-situ
    h_w 0.5–3 µs median / 17–27 µs tail confirms A-R1; local R4 j=4 shows
    no gain, expected — j=4 is work-bound, decision gate is the HPC ladder)
  - [x] clear-context review (2026-09-26, fresh-context AI audit: concurrency
    protocol verified sound, numbers recomputed clean, refactor verified pure
    code motion; verdict fixed-minor — disclosure-only fixes: bitwise
    run-to-run reproducibility requires nw % teamw == 0, held in all
    evidence; back-product cooperation flagged as A-T2 scope extension)
- **A-R3 — Paired j64 A/B vs dagteam+backoff** (1/2/4 workers per leaf,
  fixed total worker budget), measuring product speed, synchronization, LU
  time, and the time-weighted critical path before predicting solve-level gain.
  - [ ] technical completion
  - [ ] clear-context review
- **A-G — Stop gate.** Team-of-4 gain <~15% paired wall clock ⇒ Track A
  implementation does not proceed (C/D continue regardless — Ryan ruling).
  Reviewer verifies the gate was applied to uninstrumented paired data.
  - [ ] technical completion
  - [ ] clear-context review

### A Implementation

- **A-I1 — Integrate winning variant** behind runtime toggles; full
  correctness/certification suite (A-T1 contract) green.
  - [ ] technical completion
  - [ ] clear-context review
- **A-I2 — Scaling study.** 1/8/16/32/64 threads on HPC at R4 plus one larger
  mesh (fixed-problem scaling AND how available parallelism changes with mesh
  size — more leaves ≠ wider DAG). Production-default decision.
  - [ ] technical completion
  - [ ] clear-context review

## Track B — FGS setup/construction cost (priority 1, independent of A)

Ryan 2026-09-26: FGS construction is slow (~210 s more than ILU at R4);
suspected bottleneck is influence-matrix construction; profile before
spending effort anywhere.

### B Theory / profiling

- **B-T1 — Setup-cost attribution.** Profile FGS setup at R4 (and one larger
  mesh): influence-matrix/near-field block assembly vs leaf LU factorizations
  vs tree/DAG planning vs allocation/copies. No optimization work before this
  attribution exists and is reviewed.
  - [x] technical completion (2026-09-26 session 2 —
    [`033_bt1_setup_profile_20260926.md`](033_bt1_setup_profile_20260926.md):
    `nonself_influence_matrices` = 95.2% of R4 setup, single-threaded
    per-column probe; setup flat 1 vs 4 threads; near-field assembly 98.8%;
    ILU's equivalent is threaded+analytic via the 030
    `assemble_influence_block!` hook which the FGS path bypasses ⇒ the ~210 s
    setup gap is attributable to this one phase; B-G gate trivially passes.
    SPEC DEVIATION: profiled R4 + smaller R1 instead of "R4 and one larger
    mesh" — larger-mesh point to ride along with B-I1's HPC re-measure,
    pending Ryan's acceptance)
  - [x] clear-context review (2026-09-26, fresh-context AI audit: all
    mechanism claims verified in source, numbers recomputed from CSVs;
    verdict fixed-minor — doc transcription slips fixed to min-of-k values,
    no rerun; n^1.47 rung exponent is a two-point estimate, first-order only)

### B Trial

- **B-R1 — Existing-machinery audit.** (i) BRAINSTORM 030 generic block
  assembly (`assemble_influence_block!`) is already merged into FastMultipole
  production `flowpanel-20260817` (merge `ac7230a6`) — confirm whether the FGS
  setup path actually uses it, and what adopting it saves on the profiled
  bottleneck. (ii) Audit unmerged FastMultipole branches
  `faster-influence-matrices` and `influence_matrices`: contents, and whether
  they are alive against current production.
  - [x] technical completion (2026-09-26 session 2 —
    [`033_br12_machinery_and_fix_20260926.md`](033_br12_machinery_and_fix_20260926.md):
    (i) 030 hook adoption is the WRONG lever — shape/semantics mismatch (raw
    n_out rows × strength_dims columns vs influence!-projected scalar rows ×
    one value_to_strength! column; (σ,μ)=(0,1) champion ⇒ ~2× kernel work,
    σ column discarded, loses bitwise); (ii) both unmerged branches DEAD:
    `influence_matrices` is an ancestor of production,
    `faster-influence-matrices` is a 17-line never-called stub on a
    665-commit-stale base)
  - [x] clear-context review (2026-09-26, fresh-context AI audit: hook-gap
    claims verified in source, branch verdicts reproduced via git; fixed-minor)
- **B-R2 — Prototype the cheapest credible fix** for the top profiled cost;
  paired setup-time A/B with solve-phase results certified unchanged.
  - [x] technical completion (2026-09-26 session 2 — threaded the EXISTING
    probe over disjoint per-source-leaf matrices (atomic-counter pool,
    worker-private buffer copies), bitwise-identical by construction;
    `threaded_setup` kwarg, default off, legacy path untouched. Paired A/B:
    nonself 135.7→37.8 s (3.59×) at R4 j4, full ctor R1 j4 10.0→2.31 s;
    matrices + R1 cold solve certified bitwise (max|Δ|=0.0, niter 15=15).
    Projected R4 j64 setup ≈10–20 s vs ILU ~110 s — gap INVERTS;
    measurement belongs to B-I1)
  - [x] clear-context review (2026-09-26, fresh-context AI audit: concurrency
    clean — scratch is per-matrix-disjoint rhs slices, buffers worker-private,
    shared read-only source systems match the stock concurrent direct!
    contract; numbers recomputed from CSVs; verdict fixed-minor — doc-only
    fixes incl. disclosure of min-of-1 at R4 (warmups agree within 0.7%))
- **B-G — Stop gate.** Stop if profiling shows no component ≥~30% of setup,
  or if candidate machinery is dead against production and revival costs more
  than the projected saving.
  - [x] technical completion (2026-09-26 — gate PASSES, do not stop:
    dominant component 95.2% ≥ 30%; dead branches irrelevant since the
    chosen fix threads the existing probe rather than reviving them)
  - [x] clear-context review (2026-09-26 — gate inputs verified within the
    B-T1 and B-R1/B-R2 fresh-context reviews; no separate evidence to audit)

### B Implementation

- **B-I1 — Land the winning fix**; re-measure setup at R4 + larger mesh; feed
  the new setup number into the Tier-2 cumulative-crossover arithmetic.
  - [ ] technical completion
  - [ ] clear-context review

## Track C — Exact triangular solver via subdomains/interface (UNCONDITIONAL)

Preserve the GS recurrence with a partitioned lower solve: subdomain
interiors concurrent, coupling through a reduced interface system (SPIKE-like
precedent, Torun et al.). Ryan 2026-09-26: attempt regardless of A/B results.
Build order: [`fgs_exact_triangular_solver_plan_20260926.md`](021_rotor_hover_solver_benchmarks/fgs_exact_triangular_solver_plan_20260926.md).

**RYAN GATE (2026-09-26): have a conversation with Ryan before proceeding on
this track — he wants to be clear on the theory and implementation plan before
any C sub-item begins.**

### C Theory

- **C-T1 — Stage A graph-only feasibility** (laptop-scale, may run
  concurrently with Track A): interface size and fill estimates on the real
  R4 DAG (~45 predecessors/leaf makes a large interface plausible — measure,
  don't presume). Record the partitioning algebra officially.
  - [ ] technical completion
  - [ ] clear-context review
- **C-T2 — Tier-1 threshold + cost model** for the exact solve (setup,
  storage, arithmetic vs shortened sequential path), set before any prototype
  benchmark.
  - [ ] technical completion
  - [ ] clear-context review

### C Trial

- **C-R1 — Tiny algebra reference implementation** verifying exact agreement
  with the sequential lower solve on small systems; then a prototype partition
  at R4 with measured interface solve cost.
  - [ ] technical completion
  - [ ] clear-context review
- **C-G — Stop gate** (infeasibility only, not "A sufficed"): interface
  growth/fill measured prohibitive at Stage A, or prototype interface solve
  dominates the saved path. Reviewer verifies against C-T2's model.
  - [ ] technical completion
  - [ ] clear-context review

### C Implementation

- **C-I1 — Production-quality integration** + certification + HPC scaling
  study, same protocol as A-I2.
  - [ ] technical completion
  - [ ] clear-context review

## Track D — Overlapping local GS + coarse correction (UNCONDITIONAL)

Change the stationary iteration: parallel subdomain updates with halos,
coarse residual correction carrying global error. The real derivation work of
this item lives here — equations recorded officially. Ryan 2026-09-26:
attempt regardless of A/B results. No guaranteed convergence for this
operator; multigrid+multipole has precedent (Caspar) but needs
formulation-specific validation. Guidance:
[strategy](021_rotor_hover_solver_benchmarks/fgs_overlapping_gs_coarse_correction_strategy_20260926.md) ·
[implementation handoff](021_rotor_hover_solver_benchmarks/fgs_overlapping_gs_coarse_implementation_handoff_20260926.md).

**RYAN GATE (2026-09-26): have a conversation with Ryan before proceeding on
this track — he wants to be clear on the theory and implementation plan before
any D sub-item begins.**

### D Theory

- **D-T1 — Formulation.** Full derivation: overlapping local GS with halos,
  weighted reconciliation, full-residual coarse correction; iteration operator
  written out; convergence conditions/heuristics for this operator stated with
  their assumptions. Coarse-space design (what global patterns it represents,
  interface-to-interior sensitivity modes, role of cached local inverses).
  - [ ] technical completion
  - [ ] clear-context review
- **D-T2 — Tier-1 threshold + iteration budget** set before prototyping:
  the candidate must beat dagteam+backoff at matched certified accuracy
  including its iteration count (prior chunked-GS precedent: sweeps shortened
  but iterations 27→44 — the trap this design must beat).
  - [ ] technical completion
  - [ ] clear-context review

### D Trial

- **D-R1 — Reference prototype** per the implementation handoff: operator
  contracts, reference algebra, conditional enrichment, tests. Convergence
  measured on R4 (zero/nonzero starts); time-to-certified-accuracy vs
  dagteam+backoff.
  - [ ] technical completion
  - [ ] clear-context review
- **D-G — Stop gate** (infeasibility only): iteration count blows up and the
  coarse correction (with enrichment per the handoff's stop rules) does not
  recover it, or time-to-accuracy cannot beat dagteam+backoff under D-T2's
  budget.
  - [ ] technical completion
  - [ ] clear-context review

### D Implementation

- **D-I1 — Production-quality integration** + certification + HPC scaling
  study, same protocol as A-I2.
  - [ ] technical completion
  - [ ] clear-context review

## Item close-out

- **Z1 — Tier-2 headline.** Best-FGS configuration (winning tracks combined,
  including B's setup savings) vs krylov_ilu_nfcache: cold solve and
  cumulative-with-setup crossover, stated plainly either way; INDEX outcome
  cell updated; notebook entry offered to Ryan.
  - [ ] technical completion
  - [ ] clear-context review

## Deprioritized (recorded so they are not re-proposed)

- Generic coloring (tested; color boundaries limited scaling — `:colored` in
  `FastMultipole/src/solve.jl`).
- Cross-sweep pipelining (versioned state, only 3 inner sweeps, can't cross
  the outer FMM/residual boundary).
- Parallelizing the upper boundary reduction (~0.2 s/solve — smaller prize).
- NUMA team placement (follow-up only; Stage 2 showed lock contention, not
  bandwidth).
- Distributed memory; GMRES/Krylov outer iteration (out of scope by ruling).
- Another fine-grained edge scheduler (dagedge LOST; granularity was the
  failure mode).

## A-T1 — Correctness contract (official) — 2026-09-26

Contract for cooperative leaf products inside the dagteam executor. Reference
implementation being modified: `FastMultipole/src/solve_dagteam.jl:319–358`
(`dagteam_pull!` gathers predecessor strengths in ascending order and streams
the aggregated lower block; `dagteam_do_lower!` forms
$x_i \leftarrow L_i^{-1}(b_i - (Lx)_i - u_i)$ via the cached LU and publishes
successors under the queue lock). Line range verified 2026-09-26 against
branch `flowpanel-20260817`, commit `745af760`.

**Setting.** Per leaf $i$ per sweep, the update is
$$x_i \leftarrow \mathrm{LU}_i^{-1}\left(b_i - L_i\,\hat{x}_{\mathrm{pred}(i)} - u_i\right),$$
where $\hat{x}_{\mathrm{pred}(i)}$ is the gathered predecessor state at the
moment all lower predecessors have published. Cooperation changes **only** the
internal evaluation of the product $y_i = L_i\,\hat{x}_{\mathrm{pred}(i)}$;
gather order, the $b - y - u$ update, the LU solve, and publication order are
untouched.

**1. Row-partitioned product (disjoint output writes).** The $w$ workers own
disjoint contiguous row ranges of $y_i$; each entry $y_i[k]$ is computed by the
same per-row accumulation as the sequential kernel. Since no two workers write
the same entry and no per-row arithmetic order changes, $y_i$ is
**bitwise identical** to the sequential product, hence every iterate is
bitwise identical to the existing recurrence. Condition for this claim: the
row partition must not alter the per-row dot-product accumulation order
(same kernel loop order/SIMD blocking within each row) — partitioning chooses
*which worker* owns a row, never *how* a row is accumulated.

**2. Column-partitioned product (private partials + fixed-order reduction).**
Workers compute partial products $y_i^{(p)} = L_i[:, c_p]\,\hat{x}[c_p]$ over
disjoint contiguous column blocks into private vectors, then one worker
reduces in **fixed ascending block order**. In exact arithmetic
$\sum_p y_i^{(p)} = L_i\hat{x}$, so the iterates are **mathematically
equivalent** to the sequential recurrence. In floating point the summation
grouping changes (block-wise partial sums vs. one straight pass), so bitwise
match with the sequential kernel is *not* part of the contract. The contract
is instead: (a) **deterministic reduction order** — for fixed $w$ and fixed
partition, results are bitwise run-to-run reproducible; (b) exact-arithmetic
iterate equivalence as above; (c) **accuracy certification** — acceptance is
judged by time to independently certified accuracy (item ranking metric), not
by bit comparison. This matters doubly under `:f32full` (TS = Float32), where
regrouping perturbations are larger.

**3. Both layouts — synchronization contract.** The cached LU solve and
successor publication are performed **once, by one worker, only after the full
product is complete** (team barrier on all subtasks of leaf $i$). Subtasks may
start only after all lower predecessors of leaf $i$ have published (unchanged
wait-for-all semantics), so every subtask reads the same published
$\hat{x}_{\mathrm{pred}(i)}$ snapshot — no torn reads. Publication keeps the
existing lock-release fence semantics.

**4. Starts and sweep counts.** The per-leaf update is a pure function of the
gathered predecessor state; the split changes only the evaluation of
$L_i\hat{x}$. By induction over DAG topological order within a sweep and over
sweeps, equivalence of each leaf update implies equivalence of the entire
iterate sequence — for zero **and** nonzero starts and for any inner-sweep
count. (Bitwise for row partition; exact-arithmetic + certified for column
partition.) No recalibration of the certification harness is needed for the
row layout; the column layout must pass certification per (b)/(c).

## A-T2 — Split-graph critical-path model — 2026-09-26

Whole-split-graph recomputation of the byte-weighted critical path with each
leaf's lower-GEMV cost divided among $w$ cooperative subtasks (LU bytes NOT
divided) and elapsed overhead $c$ added once to each split leaf's parallel
path. In this model, each of the $w$ parallel subtasks incurs $c$, so the
additional total worker cost is $wc$, while additional elapsed cost is $c$.
Tooling:
`benchmark/fgs_dag_L_profile.jl` new `COOP_SPLIT=1` mode (existing gate-0
analysis untouched). Fixture: `phase1_case.jl` `SKIP_B=1`, champion knobs from
`retained_r4_champion.toml` (P8/MAC0.4/leaf100/f32full); run local `-t 4`,
`BENCH_BLAS_THREADS=8` (matches gate-0's banner), FLOWPanel `fastmultipole`,
FastMultipole `745af760` (`flowpanel-20260817`); both checkouts were dirty
according to the banner, so these commits are base revisions, not complete
source pins. Evidence: CSV + log in
`033_atheory_20260926/` (`fgs_dag_L_coopsplit_R4.csv`, `coopsplit_run.log`).

**Provenance correction (2026-09-26 review):** the original log records host
`tmplab-32-117-31.et.byu.edu`, four Julia threads and eight BLAS threads.
Ryan confirmed on 2026-09-26 that this is his local computer and that he does
not recall deliberately selecting eight BLAS threads. The host provenance
question is resolved; the reason for the eight-thread setting is unconfirmed.
Eight BLAS threads conflict with this item's
single-threaded-BLAS rule; matching an older banner is not an exception.
Retain the log unchanged as structural graph evidence: the model uses plan
dimensions, dependencies and scalar byte arithmetic, not measured BLAS speed.
`SKIP_B=1` means no solver validation was performed. The logged 147.0 s
construction time is not admissible setup-performance evidence. Future runs
must explicitly set `BENCH_BLAS_THREADS=1`; compliant timing and clean tagged
campaign provenance remain required for performance acceptance.

**Regression check:** unsplit L_node = 289.462 MB — exact match to gate-0
(0.000% diff), and the $c=0$ ceilings reproduce the plan doc's hand model
(w=2: 156.7 MB, w=4: 90.3 MB).

**Byte↔time conversion (gate-0 calibration):** L_node = 289.5 MB/sweep against
the measured ~2.0 s j64 sweep floor over 81 sweeps/solve ⇒ ~11.7 GB/s
critical-path streaming rate, so **1 µs ≡ 11.72 KB**. Overhead swept
0–64 µs-equivalent. Split policies: `all` (every leaf) and `selective`
(split only when it shortens elapsed pull, $g_i > c/(1-1/w)$ — the A-T3 rule).

**Results (selective policy; sweep-speedup ceiling = L_node/L_w):**

| ovh (µs) | w=2 L (MB) | w=2 ceiling | w=4 L (MB) | w=4 ceiling | w=4 n_split | w=4 tasks/sweep |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 156.7 | 1.85 | 90.3 | 3.21 | 1065 | 4263 |
| 2 | 162.7 | 1.78 | 96.4 | 3.00 | 1028 | 4152 |
| 4 | 168.7 | 1.72 | 102.4 | 2.83 | 985 | 4023 |
| 8 | 180.6 | 1.60 | 114.3 | 2.53 | 939 | 3885 |
| 16 | 203.5 | 1.42 | 137.9 | 2.10 | 808 | 3492 |
| 32 | 236.0 | 1.23 | 179.8 | 1.61 | 536 | 2676 |
| 64 | 263.6 | 1.10 | 229.7 | 1.26 | 204 | 1680 |

(`all` policy is numerically close to `selective` for $c \le 16$ µs but goes
**below 1.0×** at w=4, $c=64$ µs (0.98) and w=2, $c=64$ µs (0.81); the
selective rule floors the model at ≥1.0 and cuts task count — measured-cost
selection is not optional.)

**Where the critical path migrates.** LU share of the critical path rises
from 8% unsplit (23.9/289.5 MB) to **26% at w=4, $c=0$** (66.4 MB GEMV +
23.9 MB LU) — LU becomes material but not dominant; the un-split LU chain
(39.5 MB floor, 7.32× ceiling from gate-0) is the asymptote. Breadth never
binds at j64: $W/64 \approx 12.6$ MB $\ll L_w$ at every point, so the T64
bound equals the sweep ceiling everywhere in this model. Here $W$ is the
unsplit work; the calculation actually includes coordination using
$W_w = W + n_{\mathrm{split}}wc$ and still finds $W_w/64 < L_w$ at every
sampled point. This does not establish freedom from contention in a real
executor. Overhead re-inflates the GEMV leg, so at
large $c$ the path is overhead-dominated GEMV again, not LU.

**Implied solve-level projections** (illustrative: sweep ≈ 2.0 s of the
3.24 s j64 solve scales with the ceiling; non-sweep 1.24 s fixed):
w=4 at $c$ = 0/8/16 µs → solve ≈ 1.86/2.03/2.19 s (**1.74×/1.60×/1.48×**);
w=2 at $c$ = 0/8 µs → 2.32/2.49 s (1.40×/1.30×). The Tier-1 gate (>15%
wall-time reduction against dagteam+backoff for team-of-4, i.e. sweep ceiling
> 1.321) passes at the sampled $c=32$ µs point (1.61× sweep ⇒ 2.48 s,
23.4% less solve wall time) and fails at $c=64$ µs (1.26× sweep ⇒ 2.83 s,
12.8% less solve wall time); beating cold krylov_ilu_nfcache (2.41 s,
i.e. sweep ceiling > 1.71)
requires roughly $c \lesssim 16$–$32$ µs at w=4 (2.19 s at 16 µs, 2.48 s at
32 µs).
Task counts stay modest (≤4.3k/sweep at w=4 vs dagedge's fatal 31–49k), so
the dagedge granularity failure mode has real margin here.

**Caveats (binding):** bytes-as-cost (ignores flop-efficiency loss of short
products and cache reuse); infinite processors (real teams contend);
per-subtask overhead $c$ uncalibrated until A-R1/A-R2 measure it — dagedge's
8.2 µs must not seed it; the reduction cost of the column layout is not
separately modeled (folded into $c$); single mesh (R4), one plan; the 2.0 s
sweep-floor calibration is itself a Stage-2 estimate. These are conditional
estimates, not bounds.

## A-T3 — Split-size threshold rule (official) — 2026-09-26

**Corrected in the 2026-09-26 review:** distinguish elapsed overhead from
summed worker cost. Let $c_{\mathrm{prod}}$ be the measured sequential product
time and $h_w$ the implementation's **own measured elapsed overhead** for
the whole team (dispatch, synchronization and reduction on the elapsed path).
Under ideal product scaling, split only when
$$\frac{c_{\mathrm{prod}}}{w} + h_w < c_{\mathrm{prod}}
\quad\Longleftrightarrow\quad
c_{\mathrm{prod}} > \frac{h_w}{1 - 1/w}.$$
A-T2's parameter $c$ is $h_w$ after byte-to-time conversion. Its assumption
of parallel equal-overhead subtasks gives total worker overhead $wc$, which
belongs in the work budget, not again in the elapsed threshold. The earlier
formula $w c_{\mathrm{ovh}}/(1-1/w)$ applies only if coordination is serialized
so that $h_w = w c_{\mathrm{ovh}}$; that was not established.

The practical rule is to measure the **complete cooperative product elapsed
time** against the sequential product for the same shape and layout, and
split only when it is faster. The ideal formula is a screening model:
short products, layout changes, imbalance and contention can invalidate
the assumed $1/w$ compute scaling. Paired whole-solve acceptance is still
required at the fixed total worker budget.

Binding constraints on applying the rule:

- $h_w$ must be measured from the cooperative prototype itself
  (A-R2), **never seeded from dagedge's 8.2 µs** — that figure is average
  whole-task busy time including useful computation, reductions, atomics, and
  waits, not an isolated overhead.
- If $h_w$ lands anywhere near the dagedge scale, thresholds of
  order tens of µs overlap the median-leaf product size (~35–60 µs), so
  **which leaves split is an empirical outcome** of the measured threshold —
  do not pre-restrict to a handful of hot leaves (gate-0: splitting the top 1%
  of prio-ranked leaves moved the bound by only 1.09×; the critical path is
  broad, converging at 815/1,068 leaves split).

## Clear-context review — 2026-09-26

**APPROVED: A-T1, A-T2 and A-T3**, each against its theory-phase scope.
These address item 033's solver-cost contribution to the rotor-hover
validation program. Review checked correctness, evidence and clarity after
resolving the elapsed-overhead, projection-percentage and provenance issues.

- **A-T1:** the contract matches `dagteam_pull!` / `dagteam_do_lower!`:
  ascending lower predecessors, unchanged RHS update and cached LU, followed
  by successor publication under the queue lock. Row-layout bitwise
  equivalence requires unchanged per-row arithmetic; column-layout
  reproducibility requires deterministic partial products and reduction in
  a fixed arithmetic environment. The unchanged sweep-boundary operation
  extends the leaf-update argument across sweeps and initial states.
  These are implementation requirements, not claims about an existing
  cooperative executor.
- **A-T2:** source inspection confirms topological dynamic programming over
  lower predecessors (`j < i`), product-only splitting, unchanged LU cost,
  and separate elapsed/work overhead accounting. The original log reproduces
  the rounded gate-0 value; all 32 CSV rows passed work-budget, path-component
  and T64 arithmetic checks. Projections remain conditional estimates.
  Ryan confirmed the local host; the BLAS deviation does not invalidate
  dimension/dependency evidence. Constructor timing is excluded. No rerun
  was needed for this structural-model review.
- **A-T3:** the corrected elapsed-time break-even is algebraically consistent
  with A-T2. Measured complete-product selection and later fixed-budget
  whole-solve acceptance cover the limits of ideal scaling. No performance
  conclusion is being inferred from the screening rule alone.

No unresolved correction blocks these three theory deliverables. A-R1/A-R2
measurements and later certification are separate trial/implementation work,
not conditions on this approval. This records the AI clear-context review;
item-wide technical completion, clear-context approval and user approval
remain unchecked in INDEX because the overall program is unfinished.

## A-R1 — kernel microbenchmark results (2026-09-26, session 2)

Full record:
[`033_ar1_microbench_results_20260926.md`](033_ar1_microbench_results_20260926.md);
evidence in `033_ar1_20260926/` (CSV, banner, log); harness =
`benchmark/fgs_coop_kernel_microbench.jl` (new). Local screening run (Ryan's
machine, `-t 4`, BLAS pinned 1 — env-pinned, see below), real R4 leaf shapes,
synthetic Float32 matrices, persistent spin-barrier team. Headlines:

- **Column partitioning wins at every shape** (w=4 speedups 2.4–3.3×);
  **row partitioning rejected** — slower than sequential at w=2 on p25–p90
  short-wide shapes (strided row-block views break BLAS's fast path).
- **Measured overhead is µs-scale on median leaves**: barrier round-trip
  0.083/0.125 µs (w=2/4); effective h_w (incl. reduction + imbalance)
  0.4–3 µs on median-class leaves (~24 µs at the 1,450 tail,
  bandwidth-dominated, still clears its 31.5 µs threshold). Every shape down
  to p25 clears the A-T3
  split threshold at w=4 ⇒ empirically, essentially all leaves with ptot>0
  split. Maps to the optimistic (c≈0–8 µs) end of A-T2's projection band.
- **A-T1 amendment needed:** "row ⇒ bitwise" does not survive BLAS (sgemv
  blocking differs on row-block views; bitwise agreement is shape/w-dependent
  — recorded per config in the CSV). Both layouts must be accepted under the
  column-layout clause (deterministic + exact-arithmetic + certified
  accuracy). Determinism held in every config.
- **Tail-leaf LU now binds:** serial ldiv! on the 1,450 leaf (140.9 µs)
  exceeds its w=4 cooperative product (105.6 µs) — confirms A-T2's LU-share
  prediction and freshly motivates Track C.
- **A-T2 8-BLAS-thread mystery RESOLVED:** this host's OpenMP OpenBLAS resets
  runtime `set_num_threads(1)` on first real work (common.jl's probe
  hard-errors); the A-T2 session evidently worked around via
  `BENCH_BLAS_THREADS=8`. Correct fix: `OPENBLAS_NUM_THREADS=1
  OMP_NUM_THREADS=1` in the launch env (used here). Provenance flag closed.
- Harness lesson for A-R2: spin waiters must call `GC.safepoint()` and keep a
  rare bounded-yield escape (first harness version deadlocked without them).

## Current status

TRACK A THEORY APPROVED 2026-09-26: A-T1/A-T2/A-T3 are technically complete
and their clear-context review gates are checked. Review corrections on
2026-09-26 distinguish elapsed coordination cost from total worker cost,
correct solve wall-time percentages, and qualify the original run provenance.
A-R1 TECHNICALLY COMPLETE 2026-09-26 (session 2): column layout selected,
row rejected short-wide, effective h_w ≈ 0.4–3 µs on median leaves (barrier
~0.1 µs; ~24 µs at the 1,450 tail), A-T3 threshold cleared by all
shapes at w=4; the A-T2 eight-BLAS-thread deviation is now explained (OpenMP
OpenBLAS reset; env pinning is the fix) and the flag is closed. A-R1
clear-context APPROVED 2026-09-26 (fixed-minor).
A-R2 COMPLETE + clear-context APPROVED 2026-09-26 (session 2): column-layout
cooperative executor implemented in FastMultipole (`dagteam_coop` kwarg,
`teamw[]` runtime toggle, default bit-identical to production — refactor
verified pure code motion), R1 validation 14/14 PASS (zero/nonzero starts ×
inner 1,3 × w 1,2,4; identical iteration counts; w=4 bitwise-deterministic),
in-situ h_w on the live R4 plan 0.5–3 µs median / 17–27 µs tail (confirms
A-R1; A-T2 c-band 0–8 µs holds). Local R4 j=4 paired screening: no gain
(0.90–0.98×), expected — j=4 is work-bound (work/span 2.8); the decision
regime is span-bound j≳16. Review disclosures: run-to-run bitwise
reproducibility requires nw % teamw == 0 (held in all evidence; enforce or
re-disclose at A-I1); back-product cooperation is an A-T2 scope extension
(Ryan to accept or request model-pure variant before A-R3). Full record:
033_ar2_prototype_20260926.md + 033_ar2_20260926/. Uncommitted changes in
both checkouts; A-R3 needs commits + pinned worktrees + HPC 1/8/16/32/64
ladder.
B-T1 COMPLETE + clear-context APPROVED 2026-09-26 (session 2): FGS setup is
dominated (95.2% at R4) by the single-threaded per-source-column probe in
`nonself_influence_matrices` (FastMultipole src/solve.jl:143); setup is flat
1 vs 4 threads; near-field assembly (nonself+self) = 98.8%; tree/DAG/LU
phases ≈1.5% combined. The ~210 s setup gap vs krylov_ilu_nfcache is
attributable to this phase — ILU's near-field cache is threaded AND analytic
via the 030 `assemble_influence_block!` hook (FLOWPanel_abstractbody.jl:1354)
which the FGS setup path bypasses. Track B lever: route nonself+self through
the 030 hook, parallel over disjoint blocks; B-G gate trivially passes.
Evidence: 033_bt1_20260926/ (R1/R4 × j1/j4), harness
benchmark/fgs_setup_profile.jl (zero library edits). Spec deviation (R4+R1
instead of R4+larger) pending Ryan; larger-mesh point rides with B-I1.
B-R1/B-R2/B-G COMPLETE + clear-context APPROVED 2026-09-26 (session 2):
030-hook adoption rejected on audit (shape/semantics mismatch, ~2× kernel
work, loses bitwise); both unmerged influence-matrix branches confirmed dead.
The landed prototype instead threads the existing per-column probe over
disjoint per-source-leaf matrices (worker-private buffers, atomic-counter
pool) behind `threaded_setup` (default off, legacy path byte-untouched):
bitwise-certified (matrices max|Δ|=0.0 at R1 j1/j4 + R4 j4; R1 cold solve
niter 15=15, solution bitwise), nonself 135.7→37.8 s (3.59×) at R4 j4,
projected R4 j64 FGS setup ≈10–20 s vs ILU ~110 s — the ~210 s setup gap
inverts, pending B-I1's actual j64 measurement. Evidence: 033_br2_20260926/,
harness benchmark/fgs_setup_ab.jl, doc 033_br12_machinery_and_fix_20260926.md.
Session-2 total: A-R1, A-R2, B-T1, B-R1, B-R2, B-G all double-checked.
Remaining work is HPC-bound (A-R3 j-ladder, B-I1 landing+re-measure — both
need commits + pinned worktrees) or Ryan-gated (C/D conversations; scope and
default decisions listed in RESET BRIEF Next actions). Entry point: RESET
BRIEF.
