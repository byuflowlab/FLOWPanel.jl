# Reset prompt: FGS production dagteam — evaluator gate + end-to-end A/B (2026-09-19e)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

You are continuing BRAINSTORM 021 FGS acceleration, PRODUCTION
IMPLEMENTATION (TASK 2). Read first: `CLAUDE.md`,
`agent_policies/WORKFLOW.md` + `TESTING.md` + `HPC.md` (before corresponding
work), and in BRAINSTORM/021:
`fgs_acceleration_recommendation_20260918.md` (authoritative spec),
`fgs_acceleration_status_20260919c.md` (gates 2c+2d measured verdict),
`fgs_acceleration_status_20260919d.md` (gate-1 closure + plumbing — the
state you are inheriting).

**Champion (unchanged):** dagteam split executor + F32full + interleave 0-3
@ t16 → projected T ≈ 4.81 s = 2.10× vs baseline colored@j16 10.116 s.
Fallback ladder: f32full → f32conv (T≈5.03 s) → per-block F64 (t8 T≈6.32 s).
Accuracy gate NEVER relaxed: independent evaluator BC rel-L2 ≤ 1e-6
(`benchmark/fgs_cold_README.md`).

## What is DONE (2026-09-19d session)

- **Gate-1 CLOSED: 28/28 at -t4.** The 2 failing multi-system
  strength-recovery oracles were a fixture-class problem, NOT a code bug:
  the 1/r first-kind potential operator has block-GS spectral radius ≫ 1
  for any multi-leaf partition (union-600 SINGLE-system control diverges
  identically; the only convergent recovery test in the old suite used
  leaf_size=n_bodies = dense direct solve). Replaced with a one-sweep
  fixed-point oracle (all-direct MAC=0.01, exactly 1 sweep via
  max_iterations=1/inner_iterations=1/tolerance=0.0): clean dev ~4e-8,
  in-test negative controls (1% operator corruption) trip at ~1.6e4 for
  both schedules. FastMultipole worktree
  `/private/tmp/fastmultipole-p021-fgs-accel-20260918`, branch
  `p021-fgs-accel-20260918`, commit **`f4d6b671`** (on `29a55bf4`).
- **Regressions green:** fgs_rowpar_gate1 15/15; solve_test.jl FGS portion
  fully passing.
- **FLOWPanel plumbing DONE (live checkout, UNCOMMITTED):**
  `src/FLOWPanel_solver.jl` (FGSSolver `dagteam_precision` field+kwarg,
  forwarded ONLY when sweep_order===:dagteam so the live FastMultipole
  `c18e4b46` keeps working; FGSPreconditioner passes through),
  `src/FLOWPanel_metadata.jl` (recorded in both solver dicts),
  `benchmark/fgs_cold_common.jl` (sweep_order="dagteam" +
  `dagteam_precision` config key, v22-chunks pattern). Smoke-verified in a
  scratch env against the worktree (f64 3e-15 / f32conv 1.6e-9 / f32full
  3.4e-7 vs lex on a sphere) and live-env lex compatibility verified.
  These 3 files need a clean FLOWPanel commit (stage ONLY them — the live
  checkout carries other threads' state).

## Remaining plan (in order; ALL HPC submission Ryan-gated)

1. Commit the FLOWPanel plumbing (3 files above) once Ryan approves the
   staging; propose the commit message referencing 021 + 20260919d status.
2. **Numerical gate** (spec §gates-3): R4 cold run with
   sweep_order=dagteam, dagteam_precision=f32full; independent evaluator
   BC rel-L2 ≤ 1e-6 decides f32full vs fallback rung. F32 eps 1.2e-7 — one
   decade of headroom; the internal residual can flatter the rounded
   operator; recalibrate F32 tolerance per spec if needed. orc worktree
   `/home/rander39/wt-p021-fgs-gate2` is this thread's own; sync it to
   `f4d6b671` + FLOWPanel plumbing before any run (sync dance per
   20260919c prompt §house rules).
3. **End-to-end gate**: interleaved uninstrumented A/B vs unchanged
   colored@j16 champion on one node, ≥1.5× accepted throughput required,
   2× target; full campaign ceremony (annotated tags
   `campaign/<item>-<slug>-YYYYMMDD`, worktrees from tags, Manifest pins,
   provenance BEFORE submission; outputs to the data root).
4. Report to Ryan between gates; stop before any HPC submission.

## Traps (20260919c+d lists still apply)

- The one-sweep fixed-point oracle MUST stay at exactly 1 sweep — GS
  amplification ~5e5/sweep on this fixture blows the fixed point past
  tolerance by sweep 2. Never "strengthen" it to multi-sweep recovery.
- `dagteam.Lmat[1]` is an empty root block — corrupt/inspect ALL blocks.
- Julia 1.12 threadid can exceed nthreads (interactive pool) — never index
  scratch by threadid.
- Spinning dagteam workers must never be alive across the FMM call.
- dagteam vs lex is tolerance-equivalent, NOT bitwise; per-run determinism
  IS bitwise at any thread count.
- Local runs ≤4 threads; worktree tests:
  `julia --project=. -t4 test/fgs_dagteam_gate1_test.jl`.
- The scratch smoke env lives at
  `<scratchpad>/dagteam_env` (session-specific; rebuild if gone:
  Pkg.develop worktree + FLOWVPM.jl + FLOWPanel.jl).

## House rules (binding, unchanged)

HPC submission Ryan-gated; `ssh orc` needs a live ControlMaster socket;
monitoring via `hpc-monitor`; notebook writes Ryan-gated; dated
status/provenance files in BRAINSTORM/021. Ryan-pending ledger: origin
pushes (merged branches, v21/v22 tags, `p021-fgs-accel-20260918` now at
`f4d6b671`), FLOWPanel plumbing commit, 4+1 notebook entries,
WeakKeyDict/warmstart fix.
