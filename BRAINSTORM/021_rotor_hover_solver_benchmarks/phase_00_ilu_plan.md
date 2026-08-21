# Phase 0 W6 execution plan — sparse near-field ILU preconditioner

> **SUPERSEDED 2026-08-13.** Ryan implemented W6 directly on 2026-08-12 evening,
> collapsing Stages A–C and resolving the two open design decisions differently from the
> recommendations below: pattern = FastMultipole **Barba direct-list** (not radius
> cutoff), factorization = **`ILUZero.jl`** (as recommended). Delivered, reviewed,
> smoke-run PASS: see `phase_00_availability.md` Log 2026-08-13 and `log.md` 2026-08-13.
> Retained as design context only — do not execute.

**For a fresh agent with clean context.** Written 2026-08-12. Ryan's directive: Phase 0
approval rescinded; develop the ILU preconditioner NOW as Phase 0 W6, in three stages —
(A) theory doc, (B) implementation design, (C) implementation/testing. Read
`BRAINSTORM/INDEX.md` → control doc `021_rotor_hover_solver_benchmarks.md` (rulings,
esp. 11) → this file. The code facts below were verified 2026-08-11/12 — trust them;
spot-check line numbers only if the branch has moved (repo state at writing: branch
`fastmultipole`, everything from Phase 0 W1–W5 uncommitted on top of `b251071`).

## What already exists (do NOT rebuild)

- **Theory doc**: `theory_nearfield_ilu.md` (this dir) — the §1–§7 design basis: split
  $G = S + F$ on the FMM near-field pattern, ILU(0)/ILUT of sparse $S$, right-routed apply,
  cost model, wake subtleties, validation plan. Stage A reviews/finalizes it; it is ~90%
  done.
- **Preconditioner plumbing (no Krylov-side work needed)**: `KrylovSolver` is mutable with
  a `preconditioner` field; `_krylov_launch!` (`src/FLOWPanel_solver.jl`, search
  `_krylov_launch!`) routes `P isa FGSPreconditioner → N=, ldiv=true` and **everything else
  non-nothing → also `N=, ldiv=true`** (right preconditioning, ruling 11a). Any object
  implementing 3-arg `LinearAlgebra.ldiv!(y, P, x)` drops in via
  `KrylovSolver(body; preconditioner=P)`. `FGSPreconditioner` (same file, section
  "FGS PRECONDITIONER") is the pattern to mirror: struct holding body + saved-state
  buffers, ctor, documented `ldiv!` that is linear and side-effect-free.
- **Operator/kernel path**: dense entries come from `_G!` (`src/FLOWPanel_solver.jl:137`,
  kernel + strength column via `_G_kernel_and_strength_index`, side-aware `induced` self
  limit, sign conventions handled); matrix-free equivalents `apply_G!` /
  `_apply_neumann_G!` / `_apply_dirichlet_G!` in `src/FLOWPanel_instrumentation.jl`.
  `true_residual!` (same file) is the shared evaluator.
- **Tests to mirror**: `test/runtests_unit_solver.jl` testsets "FGSPreconditioner linearity
  and side effects" (linearity to 1e-10 + state restoration + single-leaf-exactness
  pattern) and "KrylovSolver fgmres + FGSPreconditioner (config 1f)" (solution vs
  `Backslash` to 1e-6, iteration-count comparison). Fixtures in `test/test_helpers.jl`:
  `make_octa_source_body()` (Neumann sources), `make_dirichlet_diamond_body()` (Dirichlet
  paired), `make_plate_vortex_body()` (Neumann VortexRing RigidWakeBody **with shedding**
  — use it to test the wake-terms question). Dense reference: `dense_G(body)` helper in
  `test/runtests_unit_instrumentation.jl` (assembles via `_G!`).
- **Harness**: `benchmark/rotor_hover_solver_smoke.jl` runs roster configs a–f with a
  cold-state protocol (untimed `rotor.strength .= 0` before every solve via `min_of_k`'s
  `setup!`; FGS tripwire) and writes the `decision_rules.md` CSV schema (provenance
  columns incl. `solver_settings` via `pnl._solver_metadata_dict`). Add config (g) by
  copying the config-(d) block. Rotor smoke case is **Neumann** (`RigidWakeBody{VortexRing,
  1, Float64, false}`, uncapped mesh, 6240 panels).
- **Corrected cold baselines to beat** (single-thread smoke, `rtol=1e-4`, `itmax=400`,
  from `benchmark/results/smoke/single/runs.csv`): plain GMRES 225.2 s / 290 iters to
  rms 2.6e-2; right-routed Jacobi R/4 **hurts** (400-cap, 4.9e-2); FGMRES+FGS 22.1 s /
  26 iters. ILU's success criterion at smoke scale: fewer iterations AND less wall time
  than plain GMRES at the same tolerance, setup (`t_precond`) recorded.
- `import SparseArrays` already exists in `src/FLOWPanel.jl`. `LA` = `LinearAlgebra`
  (module-level `import LinearAlgebra as LA`; also `lu!`, `ldiv!` imported).

## Stage A — Finalize the theory doc (small)

1. Read `theory_nearfield_ilu.md` end-to-end against the code; fix anything stale.
2. Resolve its two open design questions *in the doc* (they gate Stage B):
   - **Pattern source**: tree-lists vs radius-cutoff. Investigate what FastMultipole
     exposes: `FastGaussSeidel` (dev checkout `../FastMultipole`, struct in
     `src/containers.jl` ~:286-304) carries `direct_list`, `full_direct_list`,
     `strengths_by_leaf`, `source_tree` (leaf indices via `source_tree.leaf_index`;
     sort order via `source_tree.sort_index_list`) — a throwaway
     `FastMultipole.FastGaussSeidel((body,); leaf_size=...)` build is one way to get
     leaf-pair lists; `FastMultipole.Tree` directly is another. The fallback — a plain
     $O(N^2)$ radius-cutoff pair scan at construction ($\|c_i-c_j\| \le \alpha\,\bar h$,
     $\bar h$ from mean panel `characteristiclength`) — is completely adequate at smoke
    scale and has no FastMultipole coupling; **recommended default: radius cutoff first,
    tree-lists as a Phase 2b refinement** unless the tree API turns out trivial.
   - **Factorization dependency**: hand-rolled ILU(0) (~50 lines, zero new deps) vs
     `ILUZero.jl` (registered, tiny, provides `ilu0` + `ldiv!`) vs `IncompleteLU.jl`
     (ILUT, threshold knob). Adding a dep touches `Project.toml`/`Manifest.toml` — needs
     Ryan's OK. Recommendation: propose `ILUZero.jl` to Ryan (battle-tested triangular
     solves beat hand-rolled), with hand-rolled ILU(0) as the zero-dep fallback.
3. Log a dated Stage-A entry in this dir's `log.md`; get Ryan's nod on the two decisions
   (one short question, both defaults pre-recommended) before Stage B code.

## Stage B — Implementation design (doc, no src edits)

Write `design_nearfield_ilu.md` (this dir), covering concretely:

1. **Struct + API** (mirror `FGSPreconditioner`):
   ```julia
   struct ILUPreconditioner{TF,TFACT}
       S::SparseArrays.SparseMatrixCSC{TF,Int}   # near-field operator entries
       fact::TFACT                               # ILU(0)/ILUT factorization
       cutoff::Float64                           # pattern radius (physical length)
       nnz_ratio::Float64                        # nnz(S)/N² — record in metadata
   end
   ILUPreconditioner(body; cutoff_factor=..., variant=:ilu0, ...)
   LinearAlgebra.ldiv!(y, P::ILUPreconditioner, x)   # two sparse triangular solves
   ```
   No body reference needed at apply time (unlike FGSPreconditioner) — the apply is pure
   linear algebra; that also keeps `solver_state_bytes` trivially correct.
2. **Assembly**: `_S_nearfield!` — for each pattern pair $(i,j)$ evaluate the same kernel
   `_G!` uses (unit-activate the solved strength column of panel $j$, evaluate `induced`
   at control point $i$, Neumann: dot with normal; Dirichlet: potential with sources
   zeroed). Do NOT call `_G!` per pair (it does whole-matrix work); extract the inner
   kernel evaluation. Reuse `_G_kernel_and_strength_index`. Threads.@threads over columns
   like `_G!`. **Decide and document**: assemble per-pair via `induced(target, body, j,
   switch; kerneloffset)` exactly as `_G!`'s inner loop does (`src/FLOWPanel_solver.jl:175-189`)
   — this automatically includes attached-wake contributions for shedding panels
   (theory §4), because `induced` on the constructed body is wake-aware.
3. **BC channels**: Neumann (`DBC=false`) and Dirichlet (`DBC=true`) both supported,
   dispatch mirroring `_G!`'s `DBC` switch. The rotor smoke case exercises Neumann; unit
   tests cover both.
4. **Knobs**: `cutoff_factor` (pattern radius in mean-panel-lengths; sweep candidates
   {2, 4, 8}), `variant` (`:ilu0` default; `:ilut` + `τ` if IncompleteLU chosen),
   optional row equilibration toggle (theory §3 risk mitigation).
5. **Metadata**: extend `_preconditioner_metadata_dict` (`src/FLOWPanel_metadata.jl`,
   already handles FGSPreconditioner/Jacobi/nothing) with an `ILUPreconditioner` branch
   (variant, cutoff, nnz_ratio) so `solver_settings` in the CSVs is complete.
6. **Failure-mode handling** (theory §6): zero/small-pivot detection at ctor time with a
   clear error naming the mitigation (equilibration / larger cutoff / ILUT).
7. Dated log entry; brief Ryan check-in on the design (defaults pre-recommended so the
   ask is one message).

## Stage C — Implementation + testing

1. **src**: new section in `src/FLOWPanel_solver.jl` after the FGS-preconditioner section
   (or a new `FLOWPanel_ilu.jl` header if it grows past ~200 lines — register in the
   include loop in `src/FLOWPanel.jl` after "instrumentation"). No `_krylov_launch!`
   changes (the else-branch already right-routes any `ldiv!`-capable object).
2. **Tests** (extend `test/runtests_unit_solver.jl`; run
   `julia --project -e 'include("test/runtests_unit_solver.jl")'` first, then the full
   suite):
   - $S$ pattern-restricted entries ≡ `dense_G(body)` entries at pattern pairs, both BC
     types (octa + diamond fixtures), atol 1e-12.
   - With cutoff large enough that the pattern is dense: `ldiv!(y, P, b)` ≡ `G \ b`
     (ILU(0) of a dense-pattern S is exact LU) — the analogue of the FGS single-leaf
     exactness test.
   - `ldiv!` linearity (α x₁ + β x₂, 1e-10) and no body mutation.
   - End-to-end: `KrylovSolver(body; method=:gmres, preconditioner=ILUPreconditioner(body))`
     ≡ `Backslash` strengths to 1e-6 on octa (Neumann) and diamond (Dirichlet); iteration
     count ≤ unpreconditioned GMRES.
   - Wake-aware: same end-to-end on `make_plate_vortex_body()` (shedding attached) —
     guards the theory-§4 claim that `induced`-based assembly includes wake terms.
   - Metadata branch returns the new keys.
3. **Harness**: add config **(g) `krylov_ilu`** to `benchmark/rotor_hover_solver_smoke.jl`
   (copy the config-(d) block; `t_precond` = ILU ctor time; notes record cutoff + nnz
   ratio). Update the roster comment block and, in the control doc, ruling 1 (add (g)) via
   a decision-log line. Run BOTH threading modes (nohup + monitor; **the Bash tool kills
   backgrounded jobs at 10 min — use `nohup ... &` in a `dangerouslyDisableSandbox` call
   and watch the PID**, commands at the top of the smoke file; machine must be otherwise
   quiet — ask Ryan). Judge from `runs.csv`, not stdout.
4. **Acceptance** (smoke scale, single-thread): config (g) runs end-to-end; schema-valid
   rows; at matched smoke tolerance ILU-GMRES beats plain GMRES on iterations AND wall
   time (if it does not, that is a *finding to report*, not to hide — document and
   discuss with Ryan before any tuning safari).
5. **Bookkeeping**: dated entries in this dir's `log.md` + phase_00 Log (append-only;
   newest-first in log.md); update control doc `## Current status` + phase gate row;
   INDEX outcome cell; memory file `project_021_solver_benchmarks.md` (+ MEMORY.md hook
   line). Phase 0 then returns to gate-approval-pending (Ryan re-approves after W6).

## Constraints and gotchas (hard-won — do not rediscover)

- **Cold-state protocol**: FGS-family solvers seed from `body.strength`; the harness's
  untimed reset + tripwire must stay. ILU-GMRES is state-independent but gets the reset
  anyway.
- **Right-routing only** (ruling 11a): never `M=`.
- **Local runs ≤ 6 threads total** (Ryan global policy); published numbers are HPC-only
  (ruling 5) — smoke is availability evidence, not paper data.
- Smoke Krylov tolerances (`rtol=1e-4`, `itmax=400`) are local-runtime conveniences
  (ruling 11b).
- `Backslash` ctor `lu!`s in place — `solver.G` holds LU factors after construction; use
  `dense_G(body)`/`_G!` for reference entries, never `solver.G`.
- The 021 control doc is edited in place (status/gates/decision log only); session
  narrative goes to `log.md`, results to phase files. INDEX outcome cell ≤ ~35 words.
- FLOWPanel W1–W5 work is uncommitted; FastMultipole dev checkout has local commit
  `5adde3b` (callback) — don't touch FastMultipole for W6.
- Notebook entries require Ryan's explicit approval — offer, don't write.

## Suggested effort split

Stage A ~30 min (read + two decisions + one Ryan question). Stage B ~1 h (design doc).
Stage C: src+tests ~2–3 h, smoke runs ~2 h wall (background), bookkeeping ~30 min.
