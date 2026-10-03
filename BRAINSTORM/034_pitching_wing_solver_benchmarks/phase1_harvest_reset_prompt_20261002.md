# 034 Phase 1 harvest reset prompt — written 2026-10-02, campaign in flight

You are continuing BRAINSTORM item 034 (pitching-wing solver benchmarks).
Phase 0 is COMPLETE, reviewed, committed, and Ryan-approved. Phase 1 is IN
PROGRESS: the harness is built and committed, the FGS retune is done, and the
**Phase 1 consistency campaign (8 jobs) was submitted to orc 2026-10-02** and
is likely finished or near-finished by the time you read this. Your job is to
harvest it, judge the gate, and execute the Phase 1 freeze.

## Read first (in this order)

1. `BRAINSTORM/034_pitching_wing_solver_benchmarks.md` — control doc: standing
   rulings 1–11 (binding), decision log (Ryan rulings incl. campaign scope and
   the NO-notebook ruling), Phase gates table.
2. `.../decision_rules.md` — binding metrics (certified BC, gate metric
   definition, arm-promise guard, CSV conventions).
3. `.../ledger.md` — the THREE `## 2026-10-02` Phase 1 sections: FGS retune
   table, smoke/probe results, campaign pins + job table (IDs 13961478,
   13961567–573).
4. `agent_policies/HPC.md` before touching the cluster.
Do NOT read 021 item files inline; `brainstorm-scout` if needed.

## State you inherit

- Commits on `fastmultipole` (local, UNPUSHED): Phase 0 = `1d235d7`+`911a82f`;
  Phase 1 harness = `f686bc2`, `aff5272`, `43f9763` (tagged), + ledger
  commits. **GitHub origin pushes PENDING (gh auth token invalid; Ryan runs
  `gh auth login -h github.com`)** — needed for the b-tags AND the FM/VPM dev
  branches (030-era commits exist on no GitHub remote). Branch pushes
  Ryan-gated.
- **WAVE 1 (jobs 13961478, 13961567–573) FAILED — wrong dep pins** (orc
  unified-052 HEADs lack 030's `assemble_influence_block!`; ledger entry has
  the full root cause). CORRECTED pins = the local smoke banner: tag
  `campaign/p034-phase1-20261002b` on all three repos — FLOWPanel `43f9763`,
  FastMultipole `6456c221`, FLOWVPM `eebb9848`. Worktrees under orc
  `~/campaigns/p034-phase1-20261002/` rebuilt at the b pins; clean trees;
  Manifest dev-paths relative. LESSON: pins come from the loaded-package
  banner, never a clone's HEAD.
- WAVE 2 jobs (m9, --qos=normal, non-exclusive, 1 CPU + 24 G, `-t 1` strict
  single mode, one (rung, arm) per job via `benchmark/slurm/p034_phase1.sh`):
  R1 backslash/ilu/fgs/gmres = 13961625/13961629/13961630/13961631;
  R2 backslash/ilu/fgs/gmres = 13961632/13961633/13961634/13961635. ALL 8
  verified RUNNING 2026-10-02 with FLOWPanel precompiled clean; the gate job's
  banner verified. Walltimes 12–70 h (R2 gmres is the long pole). No VTK
  written (path=nothing). Verify every banner_<arm>.txt: commit **d9e4432**
  (the worktree HEAD = tag + data-symlink commit, whose PARENT is the tagged
  43f9763 — this is expected, not a mismatch), fm_commit 6456c22…, vpm_commit
  eebb984…, fm/vpm DETACHED, threading_mode single, blas_threads 1.
- Outputs: orc `~/projects/FLOWPanel.jl/data/p034_phase1/R{1,2}/` — per arm
  `summary_<arm>.csv`, `steps_<arm>.csv`, `banner_<arm>.txt`. JUDGE FROM THESE
  CSVs, never stdout/logs; verify knobs from each banner (commit must be
  43f9763, worktree DETACHED).
- FGS R1-tuned config (driver default, recorded in ledger): p=6, mac=0.3,
  leaf=150, inner=5, tol_abs=1e-6*rms_b_t0, rlx=1.0, shrink=true,
  dagteam+backoff, f64. Held on the local R2 probe (9.1e-8).
- Local smoke precedent (R1, 17 steps): all four arms meet the gate; CL agrees
  to 7 digits across the <=1e-8 arms; fgs offset 1.1e-5 rel at its 9e-8 BC.

## What to do

1. Check job status via `hpc-monitor` (read-only). Slurm FAILED/NODE_FAIL is
   advisory — judge runs by their CSVs (completed=true rows, nsteps=495).
2. Harvest (use `harvester` for the tables): per rung per arm —
   max bcerr_rel (gate: <= 1e-6 with every pass certified), n_promise_viol,
   CL/CM hysteresis agreement across arms (REPORTED, never thresholded:
   max |ΔCL| per step and over the final cycle, hysteresis loop deltas),
   niter trends, t_solve medians (context only — non-exclusive nodes).
3. Gate verdict: all four arms certify BC <= 1e-6 on BOTH rungs → Phase 1
   PASSES. Then FREEZE in the ledger: the 4-rung ladder
   (1920/6688/14336/30168) and per-rung solver/FMM settings (R1+R2 FGS =
   the tuned config above unless the campaign contradicts it; krylov seeds
   as-is given their ~1e-8/9 BC; note R3/R4 settings inherit the freeze
   prescription and get confirmed when Phase 2 first touches those rungs).
4. Write `.../phase_01_consistency.md` (phase record), update ledger +
   log.md + control doc `## Current status` / RESET BRIEF / gates row.
5. Clear-context subagent review BEFORE the Phase 1 gate box is ticked.
6. If an arm FAILED physically or missed the gate: diagnose from steps CSV
   (bcerr trend vs wake rows, promise violations), retune/resubmit the single
   arm from the SAME worktree/tag (new job IDs in the ledger); do not
   improvise new scope.

## Hard constraints (unchanged)

- Ryan-gated, offer don't do: commits, pushes (branch AND the origin tag
  push), INDEX ticks, Phase 2 start. **NO notebook entries for 034 ever**
  (Ryan ruling 2026-10-02).
- Local = smoke only, <= 4 threads, macOS BLAS caveat (launch local smokes
  with THREADING_MODE=single EXPECT_JULIA_THREADS=1 BENCH_BLAS_THREADS=8,
  `-t 1`); `benchmark/common.jl` stays untouched.
- Wake `:panel`, formulation VelocityThroughSources, cold = zero-initial-guess,
  threading modes never mixed in one comparison.
- 033's threaded repack remains UNEXERCISED in 034 (single-mode campaign);
  first exercised by Phase 2 multi mode — carry the flag.
- Phase 2 preview (do not start without Ryan): cold setup-vs-per-step cost
  across the ladder, BOTH threading modes, min-of-k, exclusive-node timing
  per 021 ruling; needs a new campaign tag.
