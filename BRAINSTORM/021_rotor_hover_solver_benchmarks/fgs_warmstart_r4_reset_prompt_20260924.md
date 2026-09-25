# Reset prompt: 021 warm-start R4 head-to-head — FGS(dagteam+backoff) vs krylov_ilu_nfcache (2026-09-24; supersedes fgs_dagedge_babysit_reset_prompt_20260924.md and the champion-adoption redo plan)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

**DRAFT — awaiting Ryan approval. Nothing below is pre-approved; every HPC
submission is Ryan-gated.**

---

## Mission and decision rule (Ryan, 2026-09-24)

The cold-solve verdict is in: at R4 j≥8, `krylov_ilu_nfcache` beats the best
FGS (dagteam+backoff, 3.24 s @ j64) at every measured j (2.41 s @ j64), and
it is the only config still scaling at 64 threads. dagedge (edge-level
pulls) LOST (+6.4%; `fgs_dagedge_benchmark_results_20260924.md`) and further
thread-scalability work is PARKED.

Production solves are not cold: unsteady stepping warm-starts every solve.
**This campaign runs the warm-started portion of Phase 3
(`phase_03_warmstart.md`, ruling 12) as a focused R4 slice**: FGS
(dagteam+backoff) vs krylov_ilu_nfcache, warm-started, wake-on,
wake-developed regime.

**Decision rule:** if warm FGS does not compete with warm
krylov_ilu_nfcache on per-step time-to-target at R4 (pre-sim setup costs
excluded from the comparison, reported separately), the FGS
acceleration thread goes **back to the drawing board** — no further FGS
optimization campaigns; solver-direction discussion with Ryan instead.
"Compete" threshold: Ryan rules on the harvested table; present per-step
medians, spreads, and cumulative-cost crossover curves, don't pre-judge.

Required reads first: `CLAUDE.md`, `agent_policies/WORKFLOW.md`,
`agent_policies/TESTING.md`, `agent_policies/HPC.md`, top-level
`BYU_ORC_AGENTS.md`. Delegation per CLAUDE.md (hpc-monitor / harvester /
hpc-storage / test-runner / code-scout). `ssh orc` needs a live
ControlMaster socket (`ssh orc -fN` if cold; never retry into 2FA).

## Step 0 (in-session, Ryan-approved 2026-09-24, before the campaign)

1. **Adopt dagteam+backoff as the FGS default** in
   `src/FLOWPanel_solver.jl`: both `FGSSolver` constructors
   (~lines 1513/1517 and ~1932): `sweep_order = :lexicographic → :dagteam`,
   `dagteam_idle = :spin → :backoff`. `dagteam_precision` stays `:f64`
   (accuracy knob, explicit opt-in); FastMultipole defaults untouched.
   Update docstrings; add the missing `:dagteam` case to
   `test/runtests_unit_solver.jl` if still absent; verify via test-runner
   (unit solver tests, ≤4 threads local); commit on `fastmultipole`.
2. **Park thread-scalability note**: append to
   `fgs_dagedge_benchmark_results_20260924.md` + ledger: headroom remains
   (77 s summed idle @ j64 even under dagedge; levers = NUMA first-touch +
   small-block repack, coarser-θ ladder, per-worker edge fusing) but is NOT
   pursued for now (Ryan 2026-09-24).

## What exists (no new plumbing expected — verify, don't rebuild)

- **Krylov warm start**: Phase 0 `x0` plumbing (`warmstart` field,
  `x_prev`, positional-x0) landed + unit-tested.
- **nfcache persistence across steps**: `KrylovSolver(persistent_plan=true,
  cache_nearfield=true)` keeps plan + near-field cache across solves with
  rigid-motion support (`transform_solver_geometry!` called by `simulate!`;
  Dirichlet operator exactly invariant under rigid motion; core_size
  restored each solve). See docstring `src/FLOWPanel_solver.jl:~945-970`
  and `rigid_motion_tree_reuse_item.md` (COMPLETE).
- **FGS warm start**: `solution_history` / `project_solution` +
  `project_solution_order` (`src/FLOWPanel_solver.jl:758-785`);
  `transform_solver!` handles rigid motion; dagteam+backoff executor from
  Stage 2. Warm-start insertion survives `set_strengths!` zeroing
  (Phase 0 W3).
- **Design check owed before the driver is written**: where does the ILU
  factorization live across steps for the nfcache config — does it persist
  with the plan (rigid-invariant like the cache) or refactorize per step?
  Locate in src (code-scout), and make the driver record
  `t_precond`/factorization events per step either way.

## Campaign design (Phase-3 slice; provenance file before submission)

- **Fixture**: R4 (58,192-panel DJI rotor, production mesh prescription per
  `project_dji_tipcap_phase2d` memory), unsteady wake-on rotor hover via
  `simulate!`, frozen Phase-1 settings, 36 steps/revolution. Run 4 full
  revolutions (144 steps). **Two scoring windows, BOTH including their
  transients (Ryan 2026-09-24 — the transient response to a changing wake
  is part of what is being tested; no settling/exclusion anywhere):**
  - **Window A — startup:** steps 0–35 (the first revolution, from the very
    first step);
  - **Window B — developed:** the fourth revolution, steps 108–143, counted
    from its first step with no wait.
  Confirm walltime against prior 021 unsteady runs at R4 before staging.
- **Convergence contract**: fixed-accuracy (matched residual target from
  Phase 1), time-to-target and iterations-to-target — NOT the 27-iteration
  fixed-work gate (that was for cold executor ranking). Both solvers
  certified at the same target; per-step true-residual spot checks.
- **Arms (all @ j64, champion placement/numactl, one process per arm,
  arms as separate Slurm tasks):**
  1. FGS dagteam+backoff, cold each step (zero guess) — baseline; "cold"
     per Ryan 2026-09-23 = zero-initial-guess solves, NOT fresh process.
  2. FGS dagteam+backoff, warm: previous solution
     (`solution_history_length=1` equivalent).
  3. FGS dagteam+backoff, warm: `project_solution=true`, order sweep
     {1, 2} (two sub-arms or in-process sequence — record order).
  4. krylov_ilu_nfcache (persistent_plan), cold each step (zero x0).
  5. krylov_ilu_nfcache (persistent_plan), warm: x0 = previous solution.
  6. (optional, Ryan may strike) krylov_ilu_nfcache, warm: extrapolated x0
     reusing `project_solution!` coefficients for comparability.
- **Knobs**: FGS = R4 champion (P8/MAC0.4/leaf100) + f32full ONLY if its
  certification transfers to the warm fixture — otherwise f64 both sides;
  ILU/nfcache = certified budget-0/500 knobs from the 2026-09-22 harvest
  (P=15/MAC=0.55, budget-500 for nfcache). BLAS pinned to 1 for BOTH arms
  and stated in provenance (ILU never BLAS-swept — flag any deviation).
- **Metrics per step** (CSV row per step per arm): step index, t_solve,
  iterations (inner + outer where applicable), residual achieved, any
  in-step precond/factorization time, warm-start projection cost, CT.
  **Pre-sim/setup costs (nfcache build, plan build, initial ILU
  factorization, solver construction) are EXCLUDED from all per-step cost
  comparisons (Ryan 2026-09-24) — reported once per arm in a separate
  setup-cost column of the summary table, not amortized in.** No
  break-even/amortization curves. Summary per window (A and B): per-step
  mean ± spread and median of t_solve and iterations, per arm.
- **Accuracy guard (binding)**: the known ~2e-3 FGS-vs-Krylov wake-on
  fixed-point discrepancy (`rigid_motion_tree_reuse_item.md` §5, open
  Phase-3 item) means the two arms converge to slightly different
  solutions once free wake exists. Record per-step cross-solver solution
  deltas and CT traces; report alongside timing — a speed win at a
  different fixed point is not a clean win. Do NOT chase the discrepancy's
  root cause in this campaign; report it.
- **Harness**: new driver `benchmark/fgs_r4_warmstart_ab.jl` reusing
  `fgs_cold_common.jl` conventions (STATUS_* per arm, COMPLETED sentinel,
  per-solve CSV schema extended with step/warm columns) + launcher
  `benchmark/run_r4_fgs_warmstart.slurm.sh`. Local smoke on R1 (≤4
  threads, few steps, all arms) before any submission.
- **Campaign discipline**: annotated tags `campaign/p021-fgs-warmstart-<date>`
  in FLOWPanel + FastMultipole + FLOWVPM; deploy per no-silos rule
  (worktrees/git-archive under `/home/rander39/campaigns/`, ARCHIVER_SKIP);
  pins in provenance file BEFORE submission; outputs to
  `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/` (or a
  new dated data root if the unsteady fixture writes heavier output);
  monitoring cadence ≥60 s; judge by outputs never sacct. **Submission =
  Ryan's explicit go on the staged sbatch line(s).**

## Job 2 — cold R4 FGS under the new default (Ryan 2026-09-24)

A second, separate job: benchmark the **cold solve** (zero-initial-guess)
on R4 with FGS under the new default (`sweep_order=:dagteam`,
`dagteam_idle=:backoff`) for direct comparison against the certified
iLU-GMRES columns already harvested (2026-09-22 table in
`fgs_scalability_reset_prompt_20260922.md`).

- **Reuse the thread-scaling harness** (`rotor_hover_solver_phase2.jl`
  measurement path / thread-scaling launcher), FGS config only — the
  krylov_ilu and krylov_ilu_nfcache columns are NOT rerun; compare against
  the existing certified rows (same rung, same data root, same BLAS=1
  convention, same placement).
- j ladder matching the existing table: j ∈ {1, 8, 16, 32, 64}, run dirs
  alongside `thread-scaling-j<J>-13777133/` (new job id dirs; do not write
  into the old ones).
- Same certification contract as those rows (staircased tolerance,
  certified BC evaluator, t_solve_min on cold isolated solves); champion
  knobs P8/MAC0.4/leaf100; f32full only if its cold R4 certification
  applies unchanged (it does — cold R4/zen3 is exactly what was certified);
  record precision in provenance either way.
- Expected reading: dagteam+backoff ≈ 3.2–3.3 s @ j64 replacing the stale
  spin-era 6.42 s regression; plateau shifts j16 → j64.
- Deliverable: the 2026-09-22 R4 table re-issued with a
  "FGS-dagteam+backoff (new default)" column, in the same results file as
  the warm harvest or its own dated file.
- Submission Ryan-gated like Job 1 (stage both sbatch lines together for
  one approval).

## Harvest deliverable

`fgs_warmstart_r4_results_<date>.md` in BRAINSTORM/021: per-arm tables for
Window A (steps 0–35, startup revolution) and Window B (fourth revolution,
steps 108–143) — iterations + time-to-target mean±spread/median per step —
with the one-time setup cost in its own column (never folded into per-step
numbers), per-step cost traces vs step index (the transient shape is a
deliverable, not noise), solution-agreement report, and the
**compete/no-compete recommendation** for Ryan's ruling. Notebook entry:
offer, don't write.

## Scratched (Ryan 2026-09-24)

The previously drafted redo plan is scratched EXCEPT the piece revived as
Job 2 above (cold R4 FGS rerun under the new default, vs the existing
certified ILU columns — no ILU reruns, no R1/R2 cold re-anchoring). The
warm comparison (Job 1) is the decisive question.

## Standing Ryan gates (carried, surface don't act)

- Notebook entry for FGS Stages 1+2 + gate-0 + dagedge campaign (offer
  once; verbosity per Ryan).
- Origin pushes: branches + tags `campaign/p021-fgs-stage2-20260923`,
  `campaign/p021-fgs-dagedge-20260924` (three repos) after
  `gh auth login -h github.com`.
- Stage-3 cancellation: confirmed by dagedge outcome.
- 018 NT-ladder jobs (13878882, 13879081–94) run concurrently — disk
  alarms are their VTK; launch hpc-storage, don't touch their queue.

## Traps

- Cold = zero-initial-guess, warm = seeded guess, SAME process either way
  (Ryan 2026-09-23).
- nfcache validity contract: rigid motion only + core_size restore before
  each solve — the production Dirichlet path does this via `simulate!`;
  assert it in the driver rather than assuming.
- Tolerances are per rung AND per environment — re-staircase, never carry.
- f32full certified at R4/zen3 COLD ONLY — do not assume it transfers to
  warm/unsteady without re-certification; default the campaign to f64.
- Task logs are output-buffered (freeze ~hours mid-compute) — judge
  liveness by CPU/outputs, not log mtime.
- Local runs ≤4 threads; laptop `runtests_benchmark_cold.jl` failure is
  pre-existing (BLAS pin) — don't chase.
