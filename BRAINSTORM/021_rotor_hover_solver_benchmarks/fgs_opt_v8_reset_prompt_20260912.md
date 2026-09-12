# Resume: FGS initialized-CPU optimization — R2 v8 job 13653852

Prepared 2026-09-12 (just after submission, 2026-09-11 evening America/Boise).
You are continuing BRAINSTORM 021's initialized-FGS CPU optimization campaign.

## Read first

1. `~/.claude/CLAUDE.md` (Campaign Reproducibility section is binding) and repo
   `CLAUDE.md`.
2. `agent_policies/HPC.md` before any cluster action;
   `/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md` on ORC.
3. The governing plan:
   `BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_initialized_cpu_optimization_plan_20260911.md`.
4. Full provenance, pins, and job description:
   `BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_opt_v8_provenance_20260911.md`
   — this file supersedes nothing there; it only adds the reset context below.

## State at reset

- R1 pilot thread CLOSED: timings 13653115 accepted (FGS prepared medians
  0.6125/0.5784/0.8758 s at 4/1, 64/1, 64/64), profiles 13653450 completed and
  validated for both solvers. R1 is supporting evidence only; do not rerun it.
- **Job 13653852 submitted 2026-09-11 evening, PENDING on m12** (normal QOS,
  zen3 exclusive, 64c/500G, 6 h; --eta estimated start ~2026-09-12 01:04).
  First execution of the v8 generation on R2: parse → precompile → controls
  j1+j4 → smoke j4/b1 (R2 seed P8/MAC0.4/leaf100/inner10) → **prepared-only**
  seed baselines j4/b1 + j64/b1 → screen j64/b1 with inner ∈ {1,2,3,5,10}
  (per-candidate staircase recalibration; failed candidates recorded and
  skipped) → 10-rep accumulated CPU profile + allocation profile of the seed.
- Output root:
  `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/opt-13653852/`
  (per-process logs `<generation>.log`, leaves under
  `<process>/R2/fgs_<id>/jN_bM/`, `COMPLETED` marker on success).
- Pins (all clean at submission; details in the provenance file): FLOWPanel
  `campaign/p021-cold-exec-20260911-v8` (`32fe1c4`) at
  `~/campaigns/p021-cold-opt-20260911-v8/FLOWPanel.jl`; FastMultipole +
  FLOWVPM at the v1 exec tags. Env + pins.toml in
  `~/campaigns/p021-cold-opt-20260911-v8/{env,pins.toml}`.
- Local development worktree `/private/tmp/flowpanel-cold-opt-20260911`
  (branch `cold-opt-20260911`, clean, pushed). The Dropbox checkout's
  `benchmark/fgs_cold_*` files are stale drafts — never use or deploy them.
- Memory (`project_021_solver_benchmarks.md`) and MEMORY.md updated. No
  notebook entry written; Ryan wants results validated first, then approval +
  detail level before writing.

## Next actions

1. Check 13653852 via the `hpc-monitor` subagent (absolute
   `/apps/slurm/latest/bin/` paths; sacct FAILED is unreliable — judge by
   outputs/logs/status.toml). If it fails, diagnose the first error in the
   per-process log; preserve partial output; a harness fix requires a new
   clean v9 generation (never edit v8 while its results are interpreted).
2. On completion, verify per leaf: `status.toml` completed, solver converged,
   BC rel-L2 ≤ 1e-6 certified (direct crosscheck active on R2), evaluator
   disagreement ≤ 1e-7, repeat agreement ≤ 1e-8, thread/pin provenance.
   A failed screen candidate (esp. inner=1) is a finding — diagnose, never
   loosen gates.
3. Harvest (use `harvester`): per-candidate `summary.csv`/`trials.csv`,
   calibration tolerances, `iterations`/`estimated_inner_sweeps`/
   `estimated_fmm_passes`. Rank the inner roster by prepared total time to
   accepted accuracy; retain the two fastest accepted candidates. Copy durable
   CSVs into a `BRAINSTORM/021_.../fgs_opt_evidence_*` dir like the pilot did.
4. Read the accumulated R2 `cpu_{flat,tree}.txt` for updated attribution
   (`.jls` deserialization only on a compute allocation). Compare against the
   R1 finding (nonself dense gemv dominant, 85/118 samples).
5. Then plan step 2: retune near/far split around the winners —
   `COLD_OPT_SCREEN_SET="leaf:25,50,200"` (and P/MAC neighbors staged after)
   reusing `run_cold_opt.slurm.sh`; a config-only screen needs no new source
   generation, code changes do (v9). Step 3 (colored sweeps + thread screen
   {1,4,16,32,64} BLAS=1) follows per the plan.
6. Report conclusions with the evidence-review pass (methods/results/
   conclusion consistency) before stating them; offer a notebook entry.

## Gotchas carried forward

- `COLD_SEEDS` has no R4 entry; add one (seed from the R3 winner) in a v9
  before the plan's R4 stage.
- ssh `orc` needs a live ControlMaster socket (2FA otherwise); if auth fails,
  stop and ask Ryan to reconnect.
- All campaign Julia execution on HPC compute nodes at the pinned Julia
  1.11.7 (cuda module loaded for precompile only). Local machine: ≤4 threads,
  logic/parse checks only.
- BLAS thread-drift assertions and the frozen gate values must never be
  relaxed; `fmm_rel_max` is an internal residual, not the evaluator
  disagreement (`evaluator_delta`).
- Never claim a general solver winner from these fixed-setting runs; claims
  are per-rung, prepared-scope, measured medians.
