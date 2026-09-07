# Phase 22 — context reset: launch Phase 3 (warmstart), R1–R3

Written 2026-09-07. Supersedes `phase_21_context_reset_prompt.md` as the entry
point. Binding rules live in `decision_rules.md` (read its tail, everything
dated 2026-09-05+); the dated record is `ledger.md` (tail: gen2c fleet entry,
2026-09-07). This file is the state + task handoff for starting Phase 3.

## Mission

Prepare and (on Ryan's explicit go) launch the Phase 3 warmstart campaign for
rungs R1–R3, per `phase_03_warmstart.md`: guess type {cold, prev, extrap} ×
solver {backslash (null control, cold only), krylov_gmres, krylov_jacobi,
krylov_ilu, fgmres_fgs, fgs}, multi-thread only. Headline deliverable:
cold-vs-warm iteration/time savings + break-even step count per solver in the
wake-developed regime. **The launch itself is Ryan's call — the Phase 3
approval checkbox in the item file is unticked. Prep everything, then ask.**

## Pin & worktree state (verify, don't trust)

- Pin of record: **gen2c** — FLOWPanel `54ff10c2` (one commit past tag
  `021-lg-gen2c`), FastMultipole `0ce3ba6` (`021-lg-gen2`), FLOWVPM `a627dd9`
  (`021-lg-gen2`). Cluster worktree: `/home/rander39/wt021/FLOWPanel.jl-c`
  (verified a proper git worktree of `~/flowpanel-021/FLOWPanel.jl`, not a
  silo). `filament_reg = LineGaussRegularization` everywhere.
- **Local commit `00390f7` (branch `fastmultipole`) is NOT on the cluster
  worktree** and Phase 3 cannot run without it: it added (a) mandatory
  Phase-2 knob sourcing under `PHASE=phase3*` in
  `benchmark/rotor_hover_solver_unsteady.jl` (`phase2_knobs`, env
  `KNOBS_BUDGET`), (b) the `KNOBS_BUDGET` submit-time guard in
  `benchmark/slurm/p3_warmstart.sh`, (c) resume-replay dedup in
  `benchmark/p021_merge.jl`. Worktree-c has jobs running from it (13593018–21
  + b0 arms) and MUST NOT be edited while they run. Path forward: when ready
  to launch Phase 3, tag a new pin (convention: `campaign/…` or the item's
  `021-lg-gen2d`) containing `00390f7`, build `FLOWPanel.jl-d` via
  `scripts/prep_campaign_worktree.sh`, point the campaign Manifest dev-paths
  at it, and copy forward the results CSVs Phase 3 reads (see below).

## Fleet state (as of 2026-09-07 ~18:00 UTC)

| jobs | what | state |
| --- | --- | --- |
| 13593015–17 | Phase 2 LineGauss tune R1–R3 | COMPLETE, harvest-clean |
| 13593018–21 | same, R4–R7 | RUNNING (3–7 d walltimes) |
| 13603412/13 | R1/R2 budget-0-only arms, `TUNE_SEED_B0=17:0.65:6` | RUNNING |

- b0 root cause (2026-09-07): the default b0 seed (16, 0.65, 6) fails the
  1e-6 BC certification gate at the START point on R1/R2 (bc 1.030e-6 /
  1.076e-6), so `tune_fmm_perturb` refuses to descend and the budget is
  skipped with only a `.err` warning. Seeds are optimality-only (recorded in
  `notes`), so seed-override relaunches are provenance-safe.
- **R4 needs the same b0 arm** once 13593018 finishes (one writer per rung
  dir — do not submit while it runs):
  `sbatch --job-name=p2lg-tune-R4-b0 --time=8:00:00
  --export=ALL,RUNG=R4,MEM_BUDGETS=0,TUNE_SEED_B0=17:0.65:6
  benchmark/slurm/p2_tune.sh` from worktree-c. **R6/R7**: their running jobs
  carry budget 0 with the OLD seed — check their `.err` for
  "budget 0.0 GiB: tuning FAILED" when they land; if failed, same treatment.
- When R4–R7 land: full harvest via `benchmark/p021_merge.jl` (now collapses
  `resumed_from_trace` replay rows; a real clobber still errors). Run it from
  a checkout containing `00390f7`.

## Phase 3 inputs — all in place on worktree-c except the driver commit

1. **Krylov apply knobs**: `phase2_knobs(rung)` reads
   `results/phase2/multi/<rung>/tune_phase2.csv`, selects the bc_certified
   row at `KNOBS_BUDGET` (mandatory env; errors list available budgets).
   The unsteady driver is UNCACHED, so **`KNOBS_BUDGET=0` is the consistent
   choice** — b0 rows exist for R3/R5; R1/R2 landing via 13603412/13.
2. **FGS knobs**: Gaussian-era `tune.csv`, `fgstune_selected.csv`,
   `fgstune_staircase.csv`, `fgsprecond.csv` copied 2026-09-07 from
   `~/projects/FLOWPanel.jl/benchmark/results/phase1/multi/` (flat) into
   worktree-c `results/phase1/multi/R{1,2,3}/`. Justified reg-inert:
   `phase1_case.jl` is a frozen single-step solve, no filament kernel runs,
   LineGauss/Gaussian measured bit-identical 2026-08-22; Ryan chose copy over
   re-run 2026-09-07. **Disclose Gaussian-provenance in any published FGS
   row.** md5s: tune 5c1421df2c7f14e9cee64a68b0b1949f, fgstune_selected
   264909990896d8b443fca178bfc14987, fgstune_staircase
   4db4dfaa688e8ce1ce7804f43274c52b, fgsprecond
   ff9d7a36c646a58dfb8d725431cd7f05. Consumer chain verified
   (stage3_winner→staircase_for→margin_tol): R1 p6/mac0.3/leaf50/inner10/
   sweeps5, R2 p8/0.4/100/10/6, R3 p6/0.3/100/5/13.
   FGS tables exist for R1–R3 ONLY and that is fine: Phase 3 is scoped R1–R3;
   R4+ never had them (agreement driver used hardcoded p10/0.4/150) and needs
   the full stage 1→2→3 chain only if Ryan ever extends Phase 3 upward.
3. **Launcher**: `benchmark/slurm/p3_warmstart.sh` — checkpoint job first
   (CONFIG=backslash runs the 13-rev staged-startup schedule once per rung),
   then each arm restarts from that checkpoint (`RESTART_STEP=-1`).
   `SKIP_STEPS=3` drops post-restart history-refill steps. Chain with
   `afterany`, never `afterok`. Hardware pinned: zen3 exclusive 128-core
   500G; never override.

## Binding rulings for Phase 3 rows (decision_rules.md)

- BC-satisfaction guard, per-arm, against the tolerance the arm promised
  (`bc_error!` every step; `t_step_net = t_step_total - t_bcerr` for any
  timestep-share claim). CT and n_particles reported, never asserted.
- Warmstart metric is `niter_first` (first per-body solve of the step);
  `niter` is non-homogeneous across solver families — never "fix" it.
- Force recovery Bernoulli-only (auto under PHASE=phase3).
- Krylov preconditioning must be right-side (`N=`), never left.
- Timing: adaptive min-of-k, exclusive node, BLAS pinned & recorded.
- `unsteady.csv` APPENDS (no skip-on-resume): rerunning a guess type
  duplicates rows; driver refuses pre-Phase-3 (narrower) headers.

## Cluster mechanics

`ssh orc` needs a live ControlMaster socket (2FA otherwise). Slurm binaries
via `/apps/slurm/latest/bin/` explicitly. Judge jobs by outputs, never by
sacct state. No sbatch/scancel without Ryan asking in the moment. Disk was
418G/400G on 2026-09-07 — being handled in a separate session; check before
adding VTK-heavy runs.

## Open with Ryan

- The Phase 3 launch itself (and the KNOBS_BUDGET=0 recommendation).
- Ledger entry for the b0 root cause + relaunch jobs + FGS CSV copy (offered
  2026-09-07, not yet approved/written — the facts are all in this file).
- Notebook entry for the campaign (offer, don't write).
