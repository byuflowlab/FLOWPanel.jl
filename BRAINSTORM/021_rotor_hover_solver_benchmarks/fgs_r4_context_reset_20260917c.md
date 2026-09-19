# Chunked hybrid sweep: context-reset handoff (2026-09-17c)

Supersedes `fgs_r4_context_reset_20260917b.md` (its task is DONE: v21
colored A/B landed and was reported — job 13738665 COMPLETED, all §6 gates
green, evidence + analysis at
`fgs_r4_followup_evidence_20260914/colored-v21-13738665/`).

**Your task (directed by Ryan 2026-09-17, invoked in plan mode): plan the
implementation of the chunked hybrid sweep** — his design, staged as item
next-step (1) in the `## Current status` block of
`../021_rotor_hover_solver_benchmarks.md`. Produce a plan file with enough
context that the implementing agent needs no exploration; save it in this
directory and stop for Ryan's approval.

## The design (Ryan's, verbatim intent)

Partition the FGS leaf sweep sequence into j contiguous chunks, one thread
each. Each thread runs Gauss–Seidel *within* its chunk (fresh values);
*across* chunks the update is Jacobi: cross-chunk leaf strengths are read
from a double-buffered snapshot taken at the top of each sweep. One barrier
per sweep (~81/solve vs colored's ~6.4k). Requirements binding on the plan:

- **Deterministic**: repeat-solution delta must be exactly 0 (§6 gate).
  The double-buffered snapshot is what makes it so — no chaotic relaxation,
  no reading a neighbor mid-write. Same result at any thread count for a
  fixed chunk map, and the chunk map must be a pure function of
  (tree, j) so runs are reproducible.
- **Separately calibrated**: chunking changes accumulation ordering and the
  iterate path, exactly like colored — new tolerance staircase per the v21
  calibrate pattern, ranked ONLY by total time to accepted accuracy.
- **Load balance**: chunk by estimated leaf cost (e.g. member count or
  influence-matrix rows), not leaf count alone; contiguous in the existing
  (Morton/tree) leaf order so most conflict edges stay intra-chunk.
- Expose as a `sweep_order` value (e.g. `:chunked`) alongside
  `:lexicographic`.
- **Coloring keep-or-revert (Ryan 2026-09-17, supersedes the earlier
  "defer until the A/B lands"): decide it in THIS plan.** Assess whether
  the chunked hybrid is more or less effective with the coloring code
  present — both algorithmically (is there any composition worth having,
  e.g. color-aware chunk boundaries or reuse of the conflict census?) and
  structurally (shared sweep-loop paths, dispatch complexity, test-suite
  weight, maintenance surface that would constrain or slow the chunked
  implementation). If the hybrid is more effective without coloring,
  **plan the revert as part of the implementation**: remove
  `sweep_order=:colored`, its scheduling machinery, and its dedicated
  tests (the 2216-case suite / `fgs_coloring_test`) via NEW commits on the
  merged branches — never rewrite history or move tags; the v21 tags
  preserve the colored evidence's reproducibility. If coloring materially
  helps the hybrid, keep it and say why. Either way the plan must state
  the verdict and its grounds explicitly for Ryan's approval. Note the
  current default expectation: the hybrid needs no conflict graph, no
  color schedule, and per-sweep barriers only — so "revert" is the likely
  answer unless the code inspection finds genuine shared value.

## Why we expect it to win (measured basis — cite, don't re-derive)

`fgs_opt_r4_diagnostics_package_20260912.md` §7 + §6: serial
`nearfield_update` chain ≈9.3 s = ~85% of j64 wall (11.2 s), proven
~1.00 active thread. v21 A/B (`colored-v21-13738665/analysis/ab_summary.md`):
colored engaged 42.65 avg threads but span did NOT shrink (barrier spin,
median-16 color width); medians lex→colored: j1 38.06→36.77, j4
17.58→15.53, j16 12.29→**10.12** (best point), j64 11.20→11.95 (loses).
Iteration robustness to ordering: colored 26 vs lex 27, arm-invariant.
Chain streams ~2.86 GB influence data/sweep (~230 GB/solve) → expect a
DRAM-bandwidth cap ~1.2–2 s on the chunked chain, j64 wall ~3.5–4.5 s.
Yardsticks to beat: lex@j64 10.96–11.20 s AND colored@j16 10.12 s.

Staged AFTER this experiment (do not plan it now): Float32 nearfield
storage (item next-step 2). Parked: Krylov outer acceleration.

## Where the code lives

- Sweep implementation: FastMultipole FGS solver, merged
  `flowpanel-20260817` tip `c6185cdd` (local checkout
  `/Users/ryan/Dropbox/research/projects/FastMultipole` — verify branch).
  `sweep_order=:colored` and its 2216-case suite + `fgs_coloring_test`
  entered on this line; the chunked path slots in beside it. Locate the
  sweep loop + strength storage before writing the plan (delegate to
  `code-scout`; do not bulk-read inline).
- A/B harness template: FLOWPanel `benchmark/fgs_r4_colored_ab.jl` +
  `benchmark/run_r4_colored_ab.slurm.sh` + unit test
  `test/runtests_r4_colored_ab_driver.jl` (commit `5f87a0e`, tag
  `campaign/p021-r4-colored-source-20260917-v21`). The v22 chunked A/B
  should be a near-clone: controls → calibrate → uninstrumented trials
  j∈{1,4,16,64} (40/order/arm, alternating batches) → one j64 activity
  pair. Three-way trials (lex/colored/chunked) are acceptable if cheap;
  minimum is chunked vs lex with colored medians cited from v21.

## Gates and process (binding)

- §6 of `fgs_opt_r4_diagnostics_package_20260912.md`: certified-FMM
  authoritative, BC rel-L2 ≤1e-6, repeat delta 0, finite, BLAS=1,
  uninstrumented rankings, iterations arm-invariant per order (verify).
- v17–v21 gate lessons: execute EVERY launcher-invoked script locally
  under a campaign-style env (direct-dep gap killed 13738561 — Meshes/
  StaticArrays must be direct deps of the campaign env); local gate env
  recipe: juliaup julia-1.11.8, Pkg.develop the three local checkouts,
  Pkg.add Meshes StaticArrays.
- Campaign reproducibility policy (`~/.claude/CLAUDE.md`): commit, tag
  (`campaign/p021-r4-chunked-*-20260918-v22` convention), fresh worktrees
  from tags, pins in provenance, env dev-pathed at worktrees, outputs to
  `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/`.
- Never edit deployed sources or move tags; sacct status is unreliable —
  judge by outputs; ≥300 s sacct spacing; `ssh orc` needs a live
  ControlMaster socket (auth failure = STOP, never MFA).

## Required reads (before planning)

`~/.claude/CLAUDE.md`, repo `CLAUDE.md`, `agent_policies/HPC.md`; cluster
work → `/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md`. Then:
1. `fgs_opt_r4_diagnostics_package_20260912.md` §7 then §6.
2. `colored-v21-13738665/analysis/ab_summary.md` (results + interpretation).
3. `v21-deployment/submission-provenance-13738561.md` (pins, design,
   postmortem — the v22 provenance template).
4. Item `## Current status` + 2026-09-17 decision-log entry
   (`../021_rotor_hover_solver_benchmarks.md`).

## Open design questions for the plan (resolve or stage for Ryan)

- Snapshot mechanics: full strength-vector copy per sweep vs per-chunk
  boundary-only buffering (copy cost vs complexity; full copy of the
  strength vector is small next to 2.86 GB of influence data — measure).
- Chunk-count decoupling: chunks == nthreads, or fixed chunk count with
  work-stealing? Fixed chunk map is required for determinism across a
  calibration↔trials pair; simplest deterministic choice = chunks == j,
  recalibrate per arm? NO — v21 calibrated once at j64 and reused across
  arms because ordering was j-invariant. Chunked ordering DEPENDS on j, so
  either (a) fix the chunk count (e.g. 64) independent of j — ordering
  j-invariant, one calibration, threads share chunks at low j — or
  (b) calibrate per arm. Prefer (a); flag for Ryan if the plan deviates.
- Inner-sweep count (`inner=3`) interaction: stale cross-chunk data ages
  over inner sweeps; decide snapshot-per-inner-sweep vs per-outer-iteration
  (per-inner-sweep is the faithful Jacobi analogue of Ryan's design).

## Ryan-pending (unchanged; do not act)

- Origin pushes of merged branches + v21 tags. (The coloring
  keep-or-revert verdict is part of the plan you produce — see the design
  section — but executing a revert still lands inside the Ryan-approved
  plan, not before it.)
- Notebook entries: v21 A/B results AND the owed diagnostics-ladder entry.
- Storage: /home ~654 G of 400 G cap; ~370 GiB RECENT VTK awaits his
  `--include-recent --only` approval. Re-run `hpc-storage` before any new
  long submission.

## Memory

`project_021_solver_benchmarks.md` is current through the v21 report;
update it (+ MEMORY.md hook) when the plan is approved / v22 lands.
