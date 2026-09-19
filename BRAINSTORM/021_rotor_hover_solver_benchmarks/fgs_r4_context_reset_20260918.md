# Chunked hybrid sweep: context-reset handoff (2026-09-18)

Supersedes `fgs_r4_context_reset_20260917c.md` — its task is DONE: the
chunked-hybrid implementation plan was produced and **Ryan APPROVED it
2026-09-17**, including both flagged items (deferred-scatter mechanics and
the conditional coloring-revert sequencing).

**Your task: IMPLEMENT the approved plan**
[`fgs_chunked_hybrid_plan_20260918.md`](fgs_chunked_hybrid_plan_20260918.md)
(this directory). The plan is the complete spec — read it IN FULL before
acting and follow it phase by phase:

- §2 Phase A — FastMultipole `sweep_order=:chunked` (fresh worktree off
  `flowpanel-20260817`; new commits only) + `test/fgs_chunked_test.jl`.
- §3 Phase B — FLOWPanel plumbing (`FGSSolver` pass-through,
  `fgs_cold_common.jl` validation) + the v22 A/B harness trio cloned from
  v21 (`benchmark/fgs_r4_colored_ab.jl`,
  `benchmark/run_r4_colored_ab.slurm.sh`,
  `test/runtests_r4_colored_ab_driver.jl` are the templates — read them).
- §7 local pre-submit gate (v17–v21 lesson chain; campaign-style 1.11.8 env
  with `Pkg.add Meshes StaticArrays`; execute EVERY launcher-invoked
  script; driver end-to-end all three AB_MODEs at R4 j4).
- §8 v22 campaign deployment (tags
  `campaign/p021-r4-chunked-{source,fm}-YYYYMMDD-v22`, worktrees from tags,
  pins.toml, env dev-pathed at worktrees; tags to cluster clones only —
  origin pushes are Ryan-pending).
- §9 watch/harvest/analysis → `chunked-v22-<job>/` evidence +
  `ab_summary.md` vs BOTH yardsticks (lex@j64 10.96–11.20 s, colored@j16
  10.116 s).
- §5 coloring revert — verdict REVERT, executed LAST and only if chunked's
  certified best point beats 10.116 s; if it loses, keep colored and report.

Approved design decisions binding on you (details/rationale in the plan):

- Deferred cross-chunk scatter realizes the double-buffered Jacobi — NO
  strength-vector copy (§1.3/§4.Q1). First implementation act: verify the
  RHS incremental semantics at FastMultipole `src/solve.jl:896-965` and
  mirror them exactly in the scatter split.
- Fixed `chunks=64`, j-invariant chunk map, ONE calibration at j64 (§4.Q2);
  barrier + cross-chunk application per INNER sweep (§4.Q3).
- Cost-balanced contiguous chunks from `self/nonself_matrices.sizes`
  (`cost_i = m_i*n_i + n_i^2`), never leaf count alone.
- Determinism gates: repeat delta exactly 0; nchunks=1 ≡ lexicographic
  bitwise; thread-count bitwise invariance; lexicographic bit-identical
  through the scatter refactor.

All file:line anchors in the plan were verified 2026-09-17 at FastMultipole
`flowpanel-20260817` @ `c6185cdd` and FLOWPanel `fastmultipole` @ `5f87a0e`.
Do not re-explore; DO re-verify anchors if either branch has moved since.

## Required reads (before acting)

`~/.claude/CLAUDE.md` (Campaign Reproducibility), repo `CLAUDE.md` (+
WORKFLOW/TESTING policies for source edits), `agent_policies/HPC.md`;
cluster work → `/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md`. Then
the plan in full, and as background only: `ab_summary.md` in
`fgs_r4_followup_evidence_20260914/colored-v21-13738665/analysis/` and
`v21-deployment/submission-provenance-13738561.md` (the v22 provenance
template).

## Standing cautions

- Never edit deployed sources or a live checkout with jobs queued/running;
  never move tags; new commits only.
- sacct status is unreliable — judge by outputs; ≥300 s sacct spacing;
  `ssh orc` needs a live ControlMaster socket (auth failure = STOP, never
  MFA).
- /home was ~654 G of the 400 G cap at v21 — run an `hpc-storage` cycle
  before submission; v22 writes only ~0.3 GB CSV/TOML, proceed with the
  breach flagged (v21 precedent) unless it has worsened.
- A chunked staircase with no certified crossing is a FINDING (harvest +
  report), not a knob-turning license (§11).

## Ryan-pending (unchanged; do not act)

- Origin pushes of merged branches + v21/v22 tags.
- Notebook entries: v21 A/B results AND the owed diagnostics-ladder entry
  (offer, never write without approval).
- Storage: ~370 GiB RECENT VTK awaiting `--include-recent --only` approval.

## Memory

`project_021_solver_benchmarks.md` is current through plan approval; update
it (+ MEMORY.md hook) when v22 lands or the campaign state changes
materially.
