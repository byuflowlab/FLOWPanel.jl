# Plan: chunked hybrid sweep (`sweep_order=:chunked`) — v22 A/B campaign

Status: **APPROVED by Ryan 2026-09-17** (planned 2026-09-17, plan-mode
directive in `fgs_r4_context_reset_20260917c.md`; approval covers both
flagged items — the §4.Q1 deferred-scatter mechanics and the §5
conditional-revert sequencing). Nothing is implemented, committed, or
submitted yet; implementation entry point =
`fgs_r4_context_reset_20260918.md`. The implementing agent should need no
exploration beyond the file:line anchors given here; all anchors were
verified 2026-09-17 against FastMultipole `flowpanel-20260817` @ `c6185cdd`
and FLOWPanel `fastmultipole` @ `5f87a0e`.

Required reads for the implementer (before touching anything):
`~/.claude/CLAUDE.md` (Campaign Reproducibility), repo `CLAUDE.md`,
`agent_policies/HPC.md`; on the cluster,
`/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md`. Campaign context:
`fgs_opt_r4_diagnostics_package_20260912.md` §6–§7,
`fgs_r4_followup_evidence_20260914/colored-v21-13738665/analysis/ab_summary.md`,
`fgs_r4_followup_evidence_20260914/v21-deployment/submission-provenance-13738561.md`.

## 0. Objective and measured basis (cite, don't re-derive)

Parallelize the serial FGS `nearfield_update` chain, which is the measured
j64 wall owner: ~9.3 s ≈ 85% of the 11.2 s j64 solve, executing at ~1.00
active thread at both j4 and j64 (§7 + v20 counters). The colored-sweep v21
A/B engaged 42.65 avg threads at j64 but its span did not shrink (79 colors
× 81 sweeps ≈ 6.4k barriers, median color width 16); best colored point is
j16 = 10.116 s. Ryan's chunked hybrid replaces ~6.4k barriers with ~81 (one
per inner sweep) at full j-way width. Chain streams ~2.86 GB influence data
per sweep (~230 GB/solve) → DRAM-bandwidth cap expectation ~1.2–2 s on the
chain, j64 wall ~3.5–4.5 s.

**Yardsticks to beat (uninstrumented medians): lex@j64 10.96–11.20 s AND
colored@j16 10.12 s.** Ranked ONLY by total time to accepted accuracy under
the §6 gates.

## 1. The algorithm (Ryan's design, made concrete)

Partition the leaf sweep sequence (tree/Morton `leaf_index` order) into
`nchunks` contiguous chunks. Within a chunk: Gauss–Seidel (fresh values, the
existing serial loop). Across chunks: Jacobi — cross-chunk influence is
evaluated at end-of-previous-sweep strengths. One barrier per inner sweep.

### 1.1 Where it slots in (FastMultipole, all `src/solve.jl` @ `c6185cdd`)

- Constructor kwarg + validation: `FastGaussSeidel(...; sweep_order=:lexicographic)`
  at `solve.jl:611-619`; colored-only precompute at `solve.jl:725-731`.
- Sweep dispatch: `gs_sweep!` at `solve.jl:1061-1121` — `:lexicographic`
  serial branch `1091-1118`, `:colored` branch `1067-1089`.
- Inner-sweep loop (`inner_iterations`, retained config inner=3):
  `solve.jl:1302-1324`, inside the outer loop `solve.jl:1233`; residual gate
  `1266-1293`; `strengths_old` outer-convergence snapshot `1230/1299/1328-1330`
  (NOT a sweep buffer — leave untouched).
- Strength storage: flat `strengths::Vector{TF}` (`solve.jl:689`,
  `containers.jl:1118`); leaf→range map `strengths_by_leaf` from
  `map_by_leaf` (`solve.jl:521-538`, contiguous per-leaf, tree leaf order).
- Per-leaf work: `solve_leaf!` (`solve.jl:88-100`, cached-LU `ldiv!`),
  `update_nonself_influence!` (`solve.jl:896-900`) =
  `compute_nonself_products!` (`905-925`) + `scatter_nonself_influence!`
  (`930-965`). The compute/scatter split is SHARED infrastructure (used by
  both existing orders) — build on it, don't duplicate it.
- Colored machinery (for §5's revert): `color_leaves` `solve.jl:985-1042`
  (~58 lines, greedy conflict-graph coloring; `overlapping_leaves`
  `999-1004`), colored `gs_sweep!` branch `1067-1089` (~23 lines), struct
  fields `containers.jl:1128-1132` (`sweep_order`, `leaf_colors`,
  `leaves_by_color`). Introduced by single commit `752c5259`. Dedicated test
  `test/fgs_coloring_test.jl` (154 lines) + 2-line `test/runtests.jl`
  registration. (`417489d5` matched a "color" grep but is the reverted 030
  Phase-3 refactor, believed coloring-unrelated — confirm with
  `git show --stat 417489d5` before writing the revert commit.)

### 1.2 Chunk map (pure function of tree + matrices + nchunks)

Precompute at construction when `sweep_order === :chunked`:

- Per-leaf cost estimate from data already on the constructed solver:
  `cost_i = m_i*n_i + n_i^2`, where `(n_i,n_i) = self_matrices.sizes[i]`
  (leaf DOF count) and `(m_i,n_i) = nonself_matrices.sizes[i]` (total
  nonself target rows this leaf influences) — `containers.jl:1092-1098`,
  populated in `self_influence_matrices` / consumed at `solve.jl:908-965`.
  This is the GEMV-flops proxy the handoff asks for; leaf count alone is
  forbidden.
- Contiguous partition of `1:length(leaf_index)` into `nchunks` chunks by
  prefix-sum of `cost`: cut at the smallest index where cumulative cost
  ≥ k/nchunks of the total (deterministic, no RNG, no thread dependence).
  Store `chunk_ranges::Vector{UnitRange{Int}}` (chunks tile the leaf index,
  each nonempty; if `nchunks > n_leaves`, clamp to `n_leaves`).
- Per-source-leaf scatter partition: for each source leaf, split its target
  segments (`index_map[i_leaf]`, built by `index_by_source`, consumed at
  `solve.jl:940`) into *intra-chunk* (target leaf in the same chunk) and
  *cross-chunk* lists, precomputed once. Pure function of (tree, chunk map).
  Determining the owning chunk of a target row uses the same
  binary-search-on-`leaf_starts` trick as `overlapping_leaves`
  (`solve.jl:999-1004`).

New struct fields (`containers.jl`, alongside `1128-1132`):
`chunks::Int`, `chunk_ranges::Vector{UnitRange{Int}}`, and the per-leaf
intra/cross scatter partition (empty for other orders).

### 1.3 The `:chunked` sweep branch

```
function gs_sweep! ... :chunked branch:
  Threads.@threads :static for c in 1:nchunks          # one barrier at end
      for i_leaf in chunk_ranges[c]                    # ascending = GS order
          solve_leaf!(...)                             # solve.jl:88-100
          compute_nonself_products!(...)               # full product, as today
          scatter_intra!(...)                          # ONLY same-chunk targets
      end
  end
  for c in 1:nchunks, i_leaf in chunk_ranges[c]        # serial, ascending
      scatter_cross!(...)                              # remaining targets
  end
```

Why this realizes Ryan's double-buffered Jacobi **without copying the
strength vector** (open question 1, resolved — see §4.Q1): during the
parallel phase a thread writes only (a) its own chunk's strength blocks and
(b) its own chunk's RHS rows. Cross-chunk contributions sit in the per-leaf
nonself product buffers (which already persist per source leaf) and are
applied after the barrier in a fixed ascending order. A target leaf in sweep
k therefore sees cross-chunk influence evaluated at each source's
end-of-sweep-(k−1) strengths — exactly the "snapshot taken at the top of
each sweep" semantics — and every RHS row receives its additions in one
deterministic order. Repeat delta is exactly 0; the result is invariant to
thread count and scheduler for a fixed chunk map, because no shared location
is ever written concurrently and no application order depends on threads.

Implementation notes (binding):

- **Verify the RHS incremental semantics first.** The scatter is
  `sorted_influences[range] .-= view(influences, index)`
  (`solve.jl:930-965`; the allocating slice at ~1355 is 021 opportunity 5 —
  do not fix it in this campaign). Establish whether products are applied as
  deltas or against a rebased RHS across the 3 inner sweeps, and mirror that
  convention exactly in `scatter_intra!`/`scatter_cross!`. Splitting the
  scatter must not change what is subtracted, only *when* and *by whom*.
- `Threads.@threads :static for c in 1:nchunks` — determinism does NOT
  depend on `:static` (chunks are write-disjoint during the parallel phase),
  but `:static` gives reproducible thread placement for the activity stage.
- The serial cross-chunk scatter is expected cheap (total scatter is 0.87 s
  of the j64 wall today, §7). If the activity stage later shows it binding,
  a deterministic parallel-by-target-chunk variant (each thread applies all
  cross contributions into its own chunk's rows, source-ascending) is the
  approved follow-up — do not build it preemptively.
- Preserve the diagnostics/observer hooks (`diagnostics[:leaf_solve_ns]`,
  `:sweep_count`, `:leaf_visit_count` at `solve.jl:1097/1108/1113/1304-1306`,
  `stage_observer`): the v22 activity mode and FLOWPanel's
  `fgs_determinism_probe.jl` consume them.
- `inner_iterations`, `tolerance`, `rlx`, `max_iterations` all live on
  `solve!`'s kwarg list (`solve.jl:1129-1134`) and are order-agnostic —
  nothing to add there. `rlx` applies only at `update_by_leaf!`
  (`solve.jl:1334`), outside the sweep.

### 1.4 Convergence expectation and risk

With 1,068 leaves / 64 chunks ≈ 17 leaves per chunk and ~89 directed edges
per leaf (95,390/1,068, §7 census), most edges are cross-chunk: the
iteration is majority-Jacobi. v21 measured ordering robustness (colored 26
vs lex 27 iterations, arm-invariant), but `fgs_coloring_test.jl:120-129`
documents a synthetic case where the Jacobi-like color-major order
*diverges* while lexicographic converges. This is exactly what the separate
calibration stage exists to catch: if the R4 staircase has no certified
crossing, that is a certification finding (like P=6/MAC=0.5), not a gate
relaxation. Expect the chunked iteration count to move off 27; the ranking
metric (total time to accepted accuracy) absorbs it.

## 2. FastMultipole changes (Phase A)

Work in a fresh local worktree off `flowpanel-20260817` (never edit the live
checkout while anything using it is queued/running; commit onto the branch
when done — new commits only, no history rewrites, no tag moves).

1. Constructor: accept `:chunked` in the `sweep_order` validation
   (`solve.jl:611-619`); add kwarg `chunks::Int=64`; build
   `chunk_ranges` + scatter partitions per §1.2 when chunked
   (`solve.jl:725-731` neighborhood); add struct fields (§1.2).
2. `gs_sweep!` `:chunked` branch per §1.3, plus the
   `scatter_intra!`/`scatter_cross!` split layered on
   `scatter_nonself_influence!` (`solve.jl:930-965`) — refactor so the
   existing orders keep calling the unsplit path unchanged
   (bit-identical lexicographic is a test gate below).
3. New test `test/fgs_chunked_test.jl` (register in `test/runtests.jl`),
   covering at minimum:
   - `nchunks=1` chunked ≡ lexicographic **bitwise** (strengths and
     `self_matrices.rhs`), same fixture style as the coloring test's
     batching theorem;
   - thread-count invariance: fixed chunk map, `JULIA_NUM_THREADS` 1 vs 4 →
     bitwise-identical results;
   - chunk-map validity + purity: chunks tile `leaf_index`, contiguous,
     nonempty, identical across repeated construction;
   - lexicographic bit-identical to pre-change reference on the same
     fixture (regression guard for the scatter refactor);
   - convergence + repeat-solution delta 0 on the gravitational fixture.
4. Run the FMM suite locally: `gravitational.jl`, `solve_test.jl`,
   `fgs_coloring_test.jl` (coloring is still present at submission — see
   §5 sequencing), `fgs_chunked_test.jl`.

## 3. FLOWPanel changes (Phase B)

On branch `fastmultipole` (again: commit before tagging; the campaign
worktree must contain no uncommitted state).

1. `src/FLOWPanel_solver.jl`: `FGSSolver` kwarg `sweep_order` passes through
   at `:1508/:1531/:1537` (and the second constructor `:1908/:1916`) — add
   `chunks::Int=64` pass-through to `FastMultipole.FastGaussSeidel`, store it
   on `FGSSolver`, and update the `:1508` comment. Record `chunks` in both
   metadata emitters (`src/FLOWPanel_metadata.jl:167,229`).
2. `benchmark/fgs_cold_common.jl`:
   - `:60` validation → `("lexicographic","colored","chunked")`, plus: if
     `sweep_order=="chunked"`, require integer `chunks ≥ 1` (default-fill 64
     in the calibrate path rather than in `cold_seed`, so existing lex
     configs are untouched);
   - `cold_make` (`:335-350`) → pass `chunks=get(c,"chunks",64)` and
     `sweep_order=Symbol(c["sweep_order"])` (already present);
   - `:141` screen-axis `"sweep_order" => ["colored"]` → update to
     `["chunked"]` (or drop) so future screens don't propose colored; this
     is independent of the §5 revert decision.
3. v22 A/B harness — near-clone of the v21 trio (commit `5f87a0e`, tag
   `campaign/p021-r4-colored-source-20260917-v21`):
   - `benchmark/fgs_r4_chunked_ab.jl` ← `fgs_r4_colored_ab.jl` (175 lines,
     read in full 2026-09-17). Deltas only: labels/filenames
     colored→chunked; `ab_calibrate` sets `sweep_order="chunked"`,
     `chunks=64`, `tolerance=0.0` then `cold_calibrate!`; `ab_check_pair`
     allows the key-set to differ by exactly `{"chunks"}` (and skips
     `sweep_order`/`tolerance` as today); env var `CHUNKED_CONFIG` replaces
     `COLORED_CONFIG`. Keep everything else — `ab_prepare` gates,
     alternating batches (`isodd(batch)`), per-trial repeat gate ≤1e-8,
     cross-order informational delta, activity mode with the macOS /proc
     guard, `status.toml` sentinel.
   - `benchmark/run_r4_chunked_ab.slurm.sh` ← `run_r4_colored_ab.slurm.sh`
     (103 lines, read in full). Deltas: job name `p021-r4-chunked-v22`, run
     dir `chunked-v22-$SLURM_JOB_ID`, driver/unit-test filenames, FMM
     control chain gains `fgs_chunked_test.jl` (keep `fgs_coloring_test.jl`
     — coloring is still in the pinned tree at submission). Same shape:
     64c zen3 exclusive 500G `--qos=normal` 12 h, controls → calibrate@j64 →
     trials j∈{1,4,16,64} (`COLD_AB_REPS=10`, `COLD_AB_BATCHES=4` → 40
     trials/order/arm, alternating batches) → one j64 activity pair. No
     perf anywhere.
   - `test/runtests_r4_chunked_ab_driver.jl` ← the 44-line v21 unit test:
     parse/AST probes retargeted, pair-check cases extended with a `chunks`
     mismatch case, alternation + no-perf + `bash -n` assertions kept.
4. **Two-way A/B (chunked vs lexicographic), colored cited from v21.**
   v21 two-way took 10.5 h of the 12 h limit; a three-way would not fit and
   adds nothing — colored's medians are banked and gate-certified. (Handoff
   explicitly permits this minimum.)

## 4. The three open design questions — resolutions

**Q1 — Snapshot mechanics: RESOLVED — neither full-vector copy nor
boundary buffering.** The deferred cross-chunk scatter (§1.3) *is* the
double buffer: per-leaf nonself product buffers hold the cross-chunk
contributions across the barrier, and their post-barrier application in
fixed order reproduces top-of-sweep-snapshot Jacobi semantics exactly.
Zero copies, no new arrays beyond the chunk map. Flagged for Ryan because
the mechanics deviate from the literal "copy the strength vector" wording
while preserving its semantics bit-for-bit; if Ryan prefers the literal
snapshot-read formulation (cross-chunk sources read from a copied vector,
target-centric), it is mathematically identical but requires a transpose
edge map and a rebased RHS — strictly more new code for the same result.

**Q2 — Chunk-count decoupling: RESOLVED — option (a), fixed `chunks=64`,
independent of j.** The chunk map and therefore the iterate path are then
j-invariant, so ONE calibration at j64 carries across arms exactly as v21's
did (v21 calibrated once because colored's ordering was j-invariant; fixed
chunks restores that property, which `chunks=j` would destroy). At low j
each thread simply runs several chunks serially; §1.3's construction makes
the result bitwise identical at any thread count. 64 matches the largest
arm and the node's physical cores. No work-stealing (would break the fixed
map). `chunks` is recorded in every selected/config TOML for provenance.
This follows the handoff's stated preference — no deviation to flag.

**Q3 — inner=3 interaction: RESOLVED — snapshot-per-inner-sweep.** The
barrier + cross-chunk application run inside `gs_sweep!`, i.e. once per
inner sweep (~81/solve), so cross-chunk data ages at most one inner sweep —
the faithful Jacobi analogue Ryan named. Per-outer snapshots would triple
staleness to save nothing (the barrier count is already the design's
selling point).

## 5. Coloring keep-or-revert — verdict: **REVERT** (grounds below), executed as the conditional final step

**Algorithmic composition — no value found.**
- The hybrid needs no conflict graph and no color schedule: cross-chunk
  correctness comes from deferred application, not conflict avoidance, so
  `color_leaves`/`overlapping_leaves`/`leaves_by_color` have no role.
- Color-aware chunk boundaries would break the contiguity that keeps edges
  intra-chunk and the cost balance that sets the critical path, for zero
  mathematical benefit — cross-chunk edges are *handled*, not avoided.
- The conflict census that motivated coloring (79 colors, median width 16)
  is banked diagnostic evidence in v20/v21; nothing at runtime needs it.
- The one genuinely shared asset — the `compute_nonself_products!` /
  `scatter_nonself_influence!` split (~60 lines) — predates-in-function and
  is used by `:lexicographic` too; it **stays regardless of the revert**.

**Structural drag from keeping it — real and growing.**
- ~90 colored-specific lines (`color_leaves` ~58, sweep branch ~23, 3
  struct fields) + the 154-line `fgs_coloring_test.jl` + control-suite time
  in every future campaign launcher.
- A second parallel-sweep semantics that every future change to the
  product/scatter path must preserve and re-verify — and the very next
  staged item (Float32 nearfield storage) rewrites exactly that path, so
  the drag is immediate, not hypothetical.
- A third `sweep_order` branch in `gs_sweep!` and in every config
  validation/screen-axis list.

**Performance case for keeping it — expected to evaporate.** Colored's only
win is j16 = 10.12 s. Chunked targets ~3.5–4.5 s at j64 with full width and
81 barriers; if it merely beats 10.12 s anywhere, colored has no production
role left.

**Sequencing (safety):** the revert lands as the final implementation step,
**conditional on the v22 gate: chunked's best certified operating point
beats colored@j16 (10.116 s)**. This keeps the only proven win in the tree
until its replacement is certified, while still making the revert part of
this approved plan per Ryan's directive. If chunked fails to beat 10.12 s,
keep coloring and report — that outcome would itself be a finding worth a
fresh decision.

**Revert contents (new commits only; v21 tags stay, preserving colored
evidence reproducibility):**
- FastMultipole (`flowpanel-20260817`): remove the `:colored` `gs_sweep!`
  branch, `color_leaves` + helpers, the three struct fields
  (`containers.jl:1128-1132` → keep `sweep_order`, drop
  `leaf_colors`/`leaves_by_color`), `:colored` from the kwarg validation,
  `test/fgs_coloring_test.jl` + its `runtests.jl` registration. Verify
  `417489d5` is untangled first (`git show --stat`).
- FLOWPanel (`fastmultipole`): drop `"colored"` from
  `fgs_cold_common.jl:60`; remove the v21 harness trio
  (`benchmark/fgs_r4_colored_ab.jl`, `benchmark/run_r4_colored_ab.slurm.sh`,
  `test/runtests_r4_colored_ab_driver.jl`) — they cannot run without
  `:colored` and remain reachable at the v21 tags; update the
  `FLOWPanel_solver.jl:1508` comment.
- Full test suites of both repos green after the revert commits.

## 6. Gates (binding, verbatim from §6 + v21 practice)

Every accepted measurement: solver converged, finite accepted solutions,
certified-FMM evaluator authoritative (no direct fallback), BC rel-L2 ≤
1e-6, **repeat-solution delta exactly 0**, ≥10 unprofiled reps per batch,
BLAS=1, pinned physical cores, instrumentation isolated from rankings,
iteration count arm-invariant per order (verify, don't assume), j1 warmup
direct-vs-FMM equivalence per order, rankings from uninstrumented trials
only. Chunked is a separately calibrated configuration (new staircase, v21
calibrate pattern); a no-certified-crossing outcome is a finding, not a
failure to suppress. Cross-order solution rel-L2 is informational only.

## 7. Local pre-submit gate (v17–v21 lesson chain, binding)

Execute EVERY script the launcher invokes, locally, under a campaign-style
env — including unchanged ones (the 13738561 lesson: env composition, not
script content, killed the first v21 submission):

1. Env: `juliaup` julia-1.11.8; fresh env; `Pkg.develop` the three local
   checkouts (FLOWPanel.jl, FastMultipole, FLOWVPM.jl); **`Pkg.add`
   `Meshes` and `StaticArrays`** (direct deps of the unit-test controls).
2. Driver unit test under 1.11.8 AND 1.12.4.
3. `cold_parse.jl`, `cold_precompile.jl`, `runtests_benchmark_cold.jl`
   (j1 + j4), `runtests_unit_solver.jl`, `runtests_unit_fgs_history.jl`,
   FMM chain (`gravitational` + `solve_test` + `fgs_coloring_test` +
   `fgs_chunked_test`) — all under the campaign-style env.
4. Driver END-TO-END at R4, all three `AB_MODE`s, j4/BLAS1: calibrate must
   certify a chunked tolerance with confirmation repeat delta 0.0 and
   `status=completed`; trials with `COLD_AB_BATCHES=2` minimum; activity
   (per-thread CSVs empty on macOS by design — guard already in the driver).
5. Record local calibrate outputs (tolerance, iterations) in the provenance
   file; they are the cluster cross-check.

## 8. Campaign deployment (v22)

Per `~/.claude/CLAUDE.md` Campaign Reproducibility + `agent_policies/HPC.md`:

- Commit everything first (both repos; campaign worktrees carry no
  uncommitted state). Annotated tags:
  - FLOWPanel `campaign/p021-r4-chunked-source-20260918-v22`
  - FastMultipole `campaign/p021-r4-chunked-fm-20260918-v22`
  - FLOWVPM unchanged at `campaign/p021-cold-exec-20260910-v1`
    (`05c658f7`, existing worktree
    `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl`).
- Push tags to the **cluster clones** (`~/projects/FLOWPanel.jl`,
  `~/projects/FastMultipole`) only — origin pushes remain Ryan-pending
  (v21 precedent); flag again for origin push on approval.
- Worktrees: `git worktree add -b <tag>-wt` at the tags under
  `/home/rander39/campaigns/p021-r4-chunked-20260918-v22/{FLOWPanel.jl,FastMultipole}`
  (HEAD == tag commit, clean status). No data-symlink commit needed — all
  v22 paths are absolute.
- Env `/home/rander39/campaigns/p021-r4-chunked-20260918-v22/env` under
  module `julia/1.11.7-6bmogfl` on a login node,
  `JULIA_PKG_PRECOMPILE_AUTO=0`: dev-path the three campaign worktrees,
  **plus `Pkg.add Meshes StaticArrays`**; verify Manifest dev-paths; write
  `pins.toml` (tags + SHAs) before submitting.
- Preflight: `hpc-storage` cycle (mandatory — /home was 654 G of the 400 G
  cap at v21; v22 writes only ~0.3 GB CSV/TOML, so proceed with the breach
  flagged to Ryan, per v21 precedent, unless it has worsened);
  `slurm-availability` probe `--cpus 64 --mem-gb 500` + `--eta` before
  choosing the partition; `ssh orc` needs a live ControlMaster socket —
  auth failure = STOP, never trigger MFA.
- Submit from the FLOWPanel worktree top level:
  `COLD_PROJECT=$camp/env CAMPAIGN_PINS=$camp/pins.toml
  COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910
  sbatch benchmark/run_r4_chunked_ab.slurm.sh`. Record job ID + `--test-only`
  ETA in the provenance file
  `fgs_r4_followup_evidence_20260914/v22-deployment/submission-provenance-<job>.md`
  (template: the v21 provenance file — pins table, pin decisions, design as
  launched, pre-submit gate evidence, storage preflight, submission record).

## 9. Watch, harvest, analysis, wrap

- Watch via `hpc-monitor`; sacct spacing ≥300 s; **judge by outputs, never
  exit status**. On failure: harvest evidence FIRST to
  `fgs_r4_followup_evidence_20260914/chunked-v22-<job>-FAILED/`, reproduce
  locally, postmortem into the next version's provenance (v17–v21 pattern).
- Harvest (delegate to `harvester`): full run dir → 
  `fgs_r4_followup_evidence_20260914/chunked-v22-<job>/`, SHA256-verified
  against remote. Analysis `analysis/ab_summary.md` mirroring v21's: gates
  table, per-arm lex/chunked medians + IQRs, comparison against BOTH
  yardsticks (lex@j64 10.96–11.20 s, colored@j16 10.12 s, cited from v21),
  iteration counts, j64 activity attribution (did the chain span shrink and
  at what avg active threads — the number colored failed on).
- Apply the §5 revert rule to the result; execute the revert commits (or
  the keep-and-report path) accordingly; run both repos' full suites.
- Update `../021_rotor_hover_solver_benchmarks.md` `## Current status` +
  decision log, `log.md`, and the memory file
  `project_021_solver_benchmarks.md` (+ MEMORY.md hook). Notebook entry:
  offer, never write without Ryan's approval (two entries already owed:
  v21 A/B + diagnostics ladder).

## 10. Acceptance criteria

1. All §6 gates green on every arm, both orders, including repeat delta
   exactly 0 and a certified chunked tolerance.
2. `fgs_chunked_test.jl` green, including nchunks=1 ≡ lexicographic bitwise
   and thread-count bitwise invariance; lexicographic bit-identical through
   the scatter refactor.
3. A/B verdict vs both yardsticks stated in `ab_summary.md`; keep-or-revert
   rule applied and executed.
4. Provenance file complete (pins, gate evidence, postmortems if any);
   evidence SHA256-verified locally.

## 11. Risks / contingencies

- **Chunked staircase fails to certify** (majority-Jacobi divergence risk,
  §1.4): finding, not failure — harvest, report; fallbacks to discuss with
  Ryan: more chunks→fewer? (fewer chunks = more GS-like, less parallel), or
  under-relaxation, both NEW experiments, not silent knob turns.
- **Iteration count inflates enough to eat the parallel win**: the ranking
  metric handles it; report against the Amdahl-style expectation
  (~3.5–4.5 s) explicitly.
- **Serial cross-chunk scatter becomes the new serial term**: visible in
  the activity pair; the deterministic parallel-by-target variant (§1.3) is
  the pre-approved follow-up.
- **12 h wall**: v21 two-way used 10.5 h; chunked arms should be faster,
  but if the local gate suggests otherwise, drop j1 trials to
  `COLD_AB_REPS=10`×`BATCHES=2` for that arm only — flag any such deviation
  in the provenance file.

## 12. Ryan-pending (unchanged; do not act)

- Origin pushes of merged branches + v21 (and future v22) tags.
- Notebook entries: v21 A/B results + owed diagnostics-ladder entry.
- Storage: ~370 GiB RECENT VTK awaiting `--include-recent --only` approval;
  /home over cap.
- Approval of THIS plan — including the §4.Q1 mechanics note and the §5
  conditional-revert sequencing — before any implementation.
