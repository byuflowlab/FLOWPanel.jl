# Reset prompt: FGS production dagteam — finish gate-1 multi-system oracle, FLOWPanel plumbing, evaluator + A/B (2026-09-19d)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

You are continuing BRAINSTORM 021 FGS acceleration, now in the PRODUCTION
IMPLEMENTATION phase (TASK 2). Read first: `CLAUDE.md`,
`agent_policies/WORKFLOW.md` + `TESTING.md` + `HPC.md` (before corresponding
work), and in BRAINSTORM/021:
`fgs_acceleration_recommendation_20260918.md` (authoritative spec),
`fgs_acceleration_status_20260919c.md` (gates 2c+2d verdict — the measured
basis for everything below).

**Gates 2c+2d verdict (harvested, Ryan approved "try it" 2026-09-19):**
champion = dagteam split executor + F32full + interleave 0-3 @ t16 →
128.2 GB/s F64-equiv, stream 1.810 s, projected T ≈ 4.81 s = **2.10×** vs
baseline colored@j16 10.116 s (beats the 2× design target). F32full beats
F32conv at every matched placement; dagteam beats rowpar everywhere.
Fallback ladder: f32full → f32conv (dagteam T≈5.03 s) → per-block F64
(dagteam F64 t8 T≈6.32 s = 1.60×). Accuracy gate NEVER relaxed:
independent evaluator BC rel-L2 ≤ 1e-6 (`benchmark/fgs_cold_README.md`).

## What is already implemented (committed)

FastMultipole worktree `/private/tmp/fastmultipole-p021-fgs-accel-20260918`,
branch `p021-fgs-accel-20260918`, commit **`29a55bf4`** (parent `e904e763`):

- `src/solve_dagteam.jl` — production split dual-layout executor:
  `FastGaussSeidel(...; sweep_order=:dagteam, dagteam_precision=:f64|:f32conv|:f32full)`.
  `build_dagteam_plan` repacks the source-major nonself matrices into
  target-major lower / source-major upper split storage from the ACTUAL
  direct-list ranges (asserts: leaf rows tile 1:n; branch rows tile whole
  leaves; no self-containing target branch; no duplicate blocks; full
  byte-tiling verified in the gate test). Readiness-counter pulls
  (one BLAS GEMV per pull, custom convert-on-load kernel for :f32conv),
  backward upper products as filler, target-owned serial boundary
  reduction, byte-weighted critical-path priority. Worker teams are
  spawned per inner-sweep block and stopped before returning (spinning
  workers would starve the threaded FMM between iterations).
  `dagteam_initialize!` primes u^0 = Ux^0 (warm starts) and reproduces the
  lex init rhs; `dagteam_inner_sweeps!` rebuilds the production invariant
  `rhs = ext + ff − Lx − Ux` at iteration boundaries, so residual/callback/
  delta/rlx paths in solve! are untouched. :f32full runs sweep state +
  leaf LU in F32 (shadow vectors, F32 LU cache); residual and outer
  bookkeeping stay F64. NO clamp (that was replay-only).
  The source-major `nonself_matrices` are RETAINED (other sweep orders,
  compatibility) — dagteam holds a second split copy of the coefficients
  (memory 2× coefficients; fine at R4, note for huge meshes).
- `src/containers.jl` — `DagTeamPlan{TM,TS,TF,TLU}`; `FastGaussSeidel`
  gains a `dagteam` field and 5th type param (only `{TF,N}` partial
  parameterizations existed elsewhere; FLOWPanel does not parameterize).
- `src/solve.jl` — constructor + solve! branches; **multi-system fill fix**:
  the nonself fill double-advanced its row cursor once per SOURCE system
  (the pre-existing two-system construction BoundsError). Rows now reuse
  the same segment rows for every source-system column block.
- `test/fgs_dagteam_gate1_test.jl` — gate-1 vs the REAL implementation.

**Gate-1 results (local, -t4): 24/26 pass.** Sweep-level dagteam vs lex
~1e-14 (zero + nonzero starts, 1 and 3 sweeps); rerun bitwise-identical;
cross-thread determinism verified bitwise t1/t3/t4 (separate probe);
end-to-end solve! lex-vs-dagteam 1.7e-11; transformed-solver fixture
3e-15; multi-system dagteam-vs-lex 6.8e-15 (constructor now works);
:f32conv 1.3e-6 / :f32full 1.6e-5 vs F64 lex (sanity, not certification).

## The 2 failing tests — multi-system truth oracle (IN PROGRESS)

The new strength-recovery oracle (solve for known strengths from their
inverted potential, 2×300 gravitational bodies, e4/MAC0.5/leaf40) returns
NaN: the multi-system solve DIVERGES (×1e6+/iteration) for lex AND dagteam
equally. Diagnosis so far (all reproduced with throwaway scripts):

- All-direct control (MAC 0.01): one sweep from the true solution moves
  strengths only 3.8e-8 → nearfield fill/scatter/strength mappings are
  CONSISTENT after the fix. All direct pairs are leaf-leaf (verified;
  `map_by_branch` gives non-leaf branches EMPTY ranges, so non-leaf direct
  targets could never scatter anyway).
- MAC 0.4: per-leaf init residual ≤ 6.6e-7 (truncation level, correct),
  yet one sweep moves leaf 43 (n=7) by 10.8 — amplification ~1e8 along the
  sweep. Leaf self-block condition numbers are modest (max ~1e4), so it is
  NOT leaf-local conditioning alone.
- **Single-system control (300 bodies, same settings) fails at
  construction with `SingularException(1)` in `build_leaf_lu_cache`** —
  the fixture produces a singular (probably 1-body) leaf self block. This
  strongly suggests the divergence is a FIXTURE conditioning problem
  (tiny near-singular leaves → block-GS spectral radius ≫ 1), not a
  multi-system code bug. NOT yet conclusive.

Next steps for this item: (1) run the union-of-600-bodies-as-ONE-system
control (constructor `Gravitational(vcat(sa.bodies, sb.bodies), zeros(34,
600))` after including test/gravitational.jl + the 6 solver-compat
overloads copied at the top of `test/fgs_dagteam_gate1_test.jl`) — if it
also diverges, the oracle fixture is the problem; (2) make the oracle
fixture well-conditioned (larger leaf_size so leaves aren't tiny, e.g.
leaf_size=100, and/or fewer bodies / larger radius_factor), confirm both
schedules recover strengths, and keep the oracle in the test; (3) if the
union control CONVERGES while 2-system diverges, there IS still a
multi-system inconsistency — bisect fill vs scatter with dense-matrix
comparison as in `test/solve_test.jl` "all influence matrices".

## Remaining plan (in order)

1. Finish the multi-system oracle (above); all gate-1 tests green at -t4.
2. Run regression suites: `test/fgs_rowpar_gate1_test.jl` (already green
   13/13 after the fill fix), FGS/solver portions of the FastMultipole
   suite, and note the FLOWPanel side is untouched so far.
3. **FLOWPanel plumbing**: pass `sweep_order=:dagteam` +
   `dagteam_precision` through FLOWPanel's `FGSSolver` (same pattern as
   the v22 `chunks` plumbing, `src/FLOWPanel_solver.jl`) and through
   `benchmark/fgs_cold_common.jl` so the R4 cold harness can select it.
4. **Numerical gate**: independent evaluator at 1e-6 on R4 decides
   :f32full vs fallback rung (F32 eps 1.2e-7 — one decade of headroom;
   internal residual can flatter the rounded operator; F32 tolerance may
   need recalibration per spec §gates-3). HPC submission is Ryan-gated.
5. **End-to-end gate**: interleaved uninstrumented A/B vs the unchanged
   colored@j16 champion on one node, ≥1.5× accepted throughput required,
   2× target; full campaign ceremony (tagged worktrees
   `campaign/<item>-<slug>-YYYYMMDD`, Manifest pins, provenance BEFORE
   submission; run from worktrees, never live checkouts).
6. Report to Ryan between gates; stop before any HPC submission.

## Traps (beyond the 20260919c list, which still applies)

- Julia 1.12: `Threads.threadid()` can EXCEED `Threads.nthreads()`
  (interactive pool) — never index per-thread scratch by threadid; use
  explicit worker indices (already fixed in dagteam_initialize!).
- Spinning dagteam workers must never be alive across the FMM call —
  keep the per-inner-sweep-block spawn/stop pattern.
- The dagteam iterate is mathematically equivalent to lex, NOT bitwise —
  compare at tolerance; per-run determinism IS bitwise at any thread count.
- reverse_pass quirk preserved: dagteam runs extra FORWARD sweeps, like lex.
- The oracle fixtures start AT the true solution (systems carry their
  strengths) — a correct solver's iteration 1 residual is already tiny;
  divergence from there means wrong operator OR unstable GS, distinguish
  via the all-direct control.
- Local runs ≤ 4 threads (house rule); `julia --project=. -t4` in the
  worktree; test standalone: `julia --project=. -t4 test/fgs_dagteam_gate1_test.jl`.

## House rules (binding, unchanged)

HPC submission Ryan-gated; `ssh orc` needs a live ControlMaster socket;
monitoring via `hpc-monitor`; notebook writes Ryan-gated; dated
status/provenance files in BRAINSTORM/021; orc worktree
`/home/rander39/wt-p021-fgs-gate2` is this thread's own; sync dance per
20260919c prompt §house rules. Ryan-pending ledger: origin pushes (merged
branches, v21/v22 tags, `p021-fgs-accel-20260918`), 4 notebook entries,
WeakKeyDict/warmstart fix.

## Suggested first moves

1. Read the two files named at the top; skim `src/solve_dagteam.jl`.
2. Run the gate-1 test at -t4 to confirm the 24/26 baseline.
3. Execute the union-600 single-system control; fix the oracle fixture
   conditioning (or the real bug, if the control exonerates the fixture).
