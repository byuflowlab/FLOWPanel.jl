# 033 session-5 reset prompt (written 2026-09-28, after B-I1 close)

Copy-paste for a fresh agent:

---

Work on BRAINSTORM item 033 (FGS scalability & setup cost) in
`~/Dropbox/research/projects/FLOWPanel.jl`. Read
`BRAINSTORM/033_fgs_scalability_and_setup_cost.md` RESET BRIEF + Track C
section + Item close-out notes first, then
`BRAINSTORM/033_hpcwave1_provenance_20260926.md` in full (failure history +
standing gotchas). Use the worker-subagent → fresh clear-context-review-
subagent pattern for every sub-item; tick both checkboxes only after a clean
review. Keep main-session context slim (delegate reading/harvesting;
hpc-monitor for job status).

State (2026-09-28, end of session 4):
- **Track A STOPPED at A-G** (ticked): j64 w=4 paired uninstrumented gain
  −18.85% vs ≥+15% gate. Per-leaf coop speedups were real (2–3.7×); the loss
  is outside-leaf sync/scheduling. `033_ar3_results_20260928.md`.
- **Track B COMPLETE through B-I1** (ticked, clean-reviewed): `threaded_setup`
  (fm b4c35f67, tag campaign/p033-hpcwave1-20260926) measured at R4/R5 j64:
  full-ctor 308.2→30.5 s (10.11×) and 789.1→85.8 s (9.20×), bitwise-certified,
  j1 parity 1.002. Tier-2 flip: FGS setup now ~78 s CHEAPER than
  krylov_ilu_nfcache; cumulative@36 dead heat (~434 vs 438 s); ILU keeps the
  cold-solve win (2.41 vs 3.24 s). `033_bi1_results_20260928.md` +
  `033_bi1_rerun_20260928/` (evidence) + `033_bi1_20260928/` (degraded
  wave-1c harvest, superseded but retained).
- No live jobs. New/modified BRAINSTORM files may be uncommitted — check
  `git status` and ask Ryan before committing.

**PRIMARY GOAL THIS SESSION: the Track C conversation with Ryan** (Ryan gate,
his directive 2026-09-28). Prepare, then converse — do NOT start any C
sub-item before he signs off on theory + implementation plan. Prep reading:
item file Track C section (C-T1/C-T2/C-R1/C-G/C-I1);
`BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_exact_triangular_solver_plan_20260926.md`
(build order); fresh motivation from A-R1/A-R3: serial ldiv! on the 1,450
tail leaf (140.9 µs) exceeds its w=4 cooperative product (105.6 µs), and
Track A's failure was outside-leaf sync — an exact partitioned lower solve
attacks the same span bound differently. Bring to the conversation: the
partitioning algebra sketch, expected interface size on the real R4 DAG
(~45 predecessors/leaf — measure, don't presume), cost model shape for C-T2,
how the Tier-1 threshold will be set, and what C-G infeasibility would look
like. Track D has the same Ryan gate if he wants to cover both.

Also on the table (Ryan-owned, raise as he wishes):
- **Ctor-tail note (Ryan 2026-09-28, recorded in Item close-out):** before
  the item completes, consider the single-threaded constructor tail that now
  dominates new setup (~20 of 30.5 s at R4, ~53 of 85.8 s at R5) — attack it
  (possible B-I2) or record as accepted residual in Z1.
- `threaded_setup` production default (B-I1 evidence strong); B-T1 R4+R1
  substitution acceptance; 021 compete decision (fed by B-I1's Tier-2
  arithmetic); push campaign tags to GitHub origin (interactive https creds);
  notebook entry for sessions 2–4 (ask verbosity per topic); empty wave-1
  output dirs on orc (`*-139012??`) deletable at his say-so.

Gotchas (all bitten this campaign, details in provenance doc): gate any
multi-job wave on ONE precompile job (VAST pidfile EACCES race); julia
pinned 1.11.7 in launchers (juliaup drift breaks HDF5_jll under 1.12);
never pass comma-valued env via `--export` lists — use `VAR=x sbatch`;
slurm bins need `export PATH=/apps/slurm/latest/bin:$PATH` in
non-interactive ssh; `ssh orc` alias needs the live ControlMaster socket
(harvester subagents may fail to connect — pull small files inline);
orc worktree HEADs are one benign `data/` symlink commit ahead of the
campaign tags (verified, provenance intact); local runs ≤4 threads +
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`; "cold" = zero-initial-guess,
arms batched per (threads, placement) process.

---
