# 033 session-6 reset prompt (written 2026-09-28, after Track C dropped)

Copy-paste for a fresh agent:

---

Work on BRAINSTORM item 033 (FGS scalability & setup cost) in
`~/Dropbox/research/projects/FLOWPanel.jl`. Read
`BRAINSTORM/033_fgs_scalability_and_setup_cost.md` RESET BRIEF + Track D
section + Item close-out notes first, then
`BRAINSTORM/033_hpcwave1_provenance_20260926.md` in full (failure history +
standing gotchas). Use the worker-subagent → fresh clear-context-review-
subagent pattern for every sub-item; tick both checkboxes only after a clean
review. Keep main-session context slim (delegate reading/harvesting;
hpc-monitor for job status).

State (2026-09-28, end of session 5):
- **Track A STOPPED at A-G** (ticked): j64 w=4 paired uninstrumented gain
  −18.85% vs ≥+15% gate; per-leaf coop speedups real (2–3.7×), loss is
  outside-leaf sync/scheduling. `033_ar3_results_20260928.md`.
- **Track B COMPLETE through B-I1** (ticked, clean-reviewed): `threaded_setup`
  (fm b4c35f67, tag campaign/p033-hpcwave1-20260926): R4/R5 j64 full-ctor
  308.2→30.5 s (10.11×) / 789.1→85.8 s (9.20×), bitwise-certified. Tier-2
  flip: FGS setup ~78 s CHEAPER than krylov_ilu_nfcache; cumulative@36 dead
  heat (~434 vs 438 s); ILU keeps cold-solve win (2.41 vs 3.24 s).
  `033_bi1_results_20260928.md`.
- **Track C DROPPED (Ryan 2026-09-28), NOT formally C-G'd** — gate
  conversation held; prep measurement on the real R4 DAG (from A-T2 exports,
  laptop-only) showed 45.1 mean preds/leaf ⇒ interface 59–91% of unknowns at
  every contiguous K; serial-interface sweep ceiling 0.29–0.75× (guaranteed
  loss, interface work alone exceeds the whole baseline sweep); parallel-
  interface rescue ≤~3× sweep-phase but re-fights the fine-grained
  scheduling battle dagedge and Track A lost, +~1.6 GB, thousands of setup
  column-solves. All C checkboxes stay unticked. Record + evidence:
  `033_trackc_prep_20260928.md` + `033_trackc_prep_20260928/`. Do NOT
  restart Track C without Ryan.
- No live jobs. New/modified BRAINSTORM files may be uncommitted — check
  `git status` and ask Ryan before committing.

**PRIMARY GOAL THIS SESSION: the Track D conversation with Ryan** (same Ryan
gate as C — converse on theory + implementation plan BEFORE starting any D
sub-item; do not touch D-T1 until he signs off). Track D = overlapping local
GS + coarse correction: changes the stationary iteration itself (parallel
subdomain updates with halos + coarse residual correction), so it is NOT
bound by the exact-span argument that killed C — but it has no guaranteed
convergence for this operator. Prep reading: item file Track D section
(D-T1/D-T2/D-R1/D-G/D-I1);
`BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_overlapping_gs_coarse_correction_strategy_20260926.md`
(principles) and `fgs_overlapping_gs_coarse_implementation_handoff_20260926.md`
(prototype plan); `033_trackc_prep_20260928.md` (the density measurement —
45 preds/leaf, interface fractions — directly informs halo sizing and
subdomain coupling strength for D). Bring to the conversation: the iteration
being proposed (operator, halo definition, coarse space), convergence risk
and how D-T2's iteration budget guards the chunked-GS trap (sweeps shortened
but iterations rose 27→44), the Tier-1 threshold (must beat dagteam+backoff
outright at matched accuracy), what D-G failure looks like (iteration
blow-up unrecoverable, or can't win under the D-T2 budget), and expected
cost/memory shape. Reuse the C-prep pattern: measure on the real R4 DAG
exports (`033_atheory_20260926/`) where a claim can be checked cheaply
before the conversation.

Also on the table (Ryan-owned, raise as he wishes):
- Ctor-tail note (Ryan 2026-09-28, in Item close-out): single-threaded
  constructor tail now dominates new setup (~20 of 30.5 s at R4, ~53 of
  85.8 s at R5) — attack (possible B-I2) or record as accepted residual in Z1.
- `threaded_setup` production default (B-I1 evidence strong); B-T1 R4+R1
  substitution acceptance; 021 compete decision (fed by B-I1's Tier-2
  arithmetic); push campaign tags to GitHub origin (interactive https creds);
  notebook entry for sessions 2–5 (ask verbosity per topic); empty wave-1
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
