# Reset prompt: 021 warm-start R4 campaign — fix the winB restart smoke, redeploy, get Ryan's go (2026-09-24b; supersedes fgs_warmstart_r4_reset_prompt_20260924.md)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Mission

Finish staging the 021 warm-start R4 campaign
(`fgs_warmstart_r4_provenance_20260924.md` is the provenance of record;
the original spec `fgs_warmstart_r4_reset_prompt_20260924.md` is fully
implemented EXCEPT one open smoke defect). In order:

1. **Root-cause and fix the winB restart smoke failure** (details below —
   the ONLY blocker).
2. Re-run the leg smoke to PASS; record it in the provenance.
3. Commit; refresh the campaign tag + orc deployment (steps below).
4. Surface the staged sbatch lines to Ryan and **submit only on his
   explicit go** (decision points listed at the end).
5. After submission: babysit (hpc-monitor), then harvest
   (`benchmark/fgs_r4_warmstart_harvest.jl` + agent-written
   `fgs_warmstart_r4_results_<date>.md` with the compete/no-compete
   recommendation — Ryan rules; see the provenance "Harvest" section).

Required reads first: `CLAUDE.md`, `agent_policies/WORKFLOW.md`,
`agent_policies/TESTING.md`, `agent_policies/HPC.md`, cluster
`BYU_ORC_AGENTS.md`, and the provenance file. `ssh orc` needs a live
ControlMaster socket (`ssh orc -fN` if cold; never retry into 2FA).
Local runs ≤4 threads.

## State (all committed on `fastmultipole`, all local — origin push owed)

- `08667c2` — dagteam+backoff = FGSSolver/FGSPreconditioner default;
  t_project instrumentation. FastMultipole `745af760`
  (branch flowpanel-20260817) fixes the dagteam empty-direct-list
  BoundsError the default flip exposed. Unit solver suite 513/513.
- `f88e35a` — campaign harness Jobs 1+2 (per-arm driver
  `benchmark/fgs_r4_warmstart_ab.jl` over the extended
  `benchmark/rotor_hover_solver_unsteady.jl`; launchers
  `run_r4_fgs_warmstart.slurm.sh` / `run_r4_fgs_cold_newdefault.slurm.sh`;
  harvest; smoke). Includes `KrylovSolver(rtol_rhs=…)`: Krylov.jl's rtol is
  r0-relative so warm starts gained nothing by construction — rtol_rhs
  anchors the stopping target to ‖b‖ (cold-identical; phase3-only in the
  driver). Full-march smoke (7 arms, R1) PASSED; three defects found and
  fixed pre-HPC (zero-rhs Dirichlet priming; r0-relative stopping;
  effective-tol read after potential overwrite — read `solver.rhs`).
- `db445cf` — pins filled; tags `campaign/p021-fgs-warmstart-20260924` cut
  in all three repos (FLOWPanel `f88e35a`, FastMultipole `745af760`,
  FLOWVPM `8d4a3b4`); rsync-mode deployment to
  `/home/rander39/campaigns/p021-fgs-warmstart-20260924/` sha256-VERIFIED;
  campaign env dev-pointed; sbatch lines staged in the provenance.
- `0e38c21` — **checkpoint+restart leg layout (Ryan's second 2026-09-24
  ruling)**: TWO family checkpoints (fgs / ilu). Legs per (arm, leg)
  process via `WSR4_LEG`: `ckpt` = cold arm marches revs 1–3 (3NT steps)
  with SAVE_VTK=true (its own Window A + the family restart source);
  `winA` = warm arms march rev 1; `winB` = ALL 7 arms restart rev 4 from
  the family checkpoint (`RESTART_STEP=3NT` via `simulate_warmstart!`);
  `full` = old single-march (kept). Launcher runs 14 legs sequentially,
  winB gated on family-ckpt STATUS. Harvest is leg-aware (Window A =
  non-restarted steps 1..NT; Window B = restarted rows; global step =
  restart_step + local step) and prints the REQUIRED reporting note:
  solver warm-start histories are not serialized in checkpoints, so each
  restarted warm leg's first (order+1) steps are effectively cold —
  inside Window B's transient-included stats by design (Ryan: fine, note
  when reporting).

## THE BLOCKER: winB restart legs fail the R1/NT=4 smoke

Smoke (`benchmark/run_r4_fgs_warmstart_smoke.sh`, output in
`benchmark/results/wsr4_smoke_legs/`): both `ckpt` legs PASS, both `winA`
legs PASS. The restart legs are broken in BOTH families:

- `fgs_cold_winB`, `fgs_prev_winB`: die at the SECOND restarted step
  (step 14/16) with "block Gauss-Seidel produced a nonfinite physical
  residual at outer iteration 1" (`src/FLOWPanel_solver.jl:~2625`). The
  first restarted step converges but slowly (~35 s vs ~10 s forward).
- `ilu_nfcache_prev_winB`: sentinel reads ok but is PHYSICALLY GARBAGE —
  every restarted step hits itmax=500, CF ~1e7–1e8, checksum 1e7. (Also a
  harness gap: the leg sentinel should not read ok when every step is
  unconverged — consider failing the wrapper when n_unconverged > 0.)

Evidence banked:
- Restored VTK state is ALL FINITE (probe scripts in the session
  scratchpad; checkpoint step-12 particle VTP + body VTU checked field by
  field). Restart resolves correctly ("resuming from step 12 (file count
  13)").
- So the restored state goes bad at/after the first restarted solve, for
  both solver families ⇒ NOT an FGS-tree staleness issue alone.
- A Float32-particle-VTP warning fires on restore ("continuation will not
  be replay-exact") — the ckpt wrote f32 particles. Worth checking
  FLOWPANEL_PARTICLE_PRECISION (f64 default claimed in
  `src/FLOWPanel_warmstart.jl:~250`) — why is the smoke ckpt f32?
- A "core_size 0.001 large relative to panel radius" warning appears in
  both forward and restarted logs — probably not the discriminator.
- **The decisive control was IN FLIGHT at handoff**: a FORWARD full-march
  control at the same fixture (ARM=fgs_cold WSR4_LEG=full RUNG=R1 NT=4
  N_STEPS=16) to discriminate "NT=4 fixture is physically unstable past
  step 12" vs "the restart path is inconsistent". Its output dir:
  `<scratchpad>/debug_fwd_out` (session-scratchpad, likely gone) — just
  RERUN it (one env line, ~8 min; copy the env block from
  `run_r4_fgs_warmstart_smoke.sh` and override ARM/WSR4_LEG/N_STEPS).
  - If the forward control ALSO blows up at step ~14: the smoke fixture
    (NT=4, dt=T/4) is the problem, not the restart. Fix the smoke, e.g.
    restart at 2NT=8 and march to 12 (`RESTART_STEP=8 N_STEPS=12` env
    overrides for the winB stages — the wrapper's `_setdefault!` lets
    explicit env win), or raise NT to 8 and accept a longer smoke.
  - If the forward control is CLEAN through step 16: the restart path has
    a real inconsistency. Suspects, in order: Das/eta kinematic state
    across restart (RHPC uses DAS_ARC_HELIX_SOURCE=steady; check
    set_Das_* replay in `simulate_warmstart!`), the f32 particle
    round-trip, wake attachment/Kutta terminal-strength restore, and
    `body.core_size` restore vs the solver-construction value. Note the
    known ~1092-particle restored field is finite, and 018 production
    restarts at NT≥36 work routinely on HPC — a tiny-NT-only bug is
    plausible but must be PROVEN before trusting the R4 campaign's
    restart legs (the campaign runs NT=36; if the bug is NT-scale only,
    demonstrate a clean restart smoke at a fixture closer to production,
    e.g. NT=12 ckpt=36 steps, before staging).
- Debug harness: `<scratchpad>/debug_winb.jl` pattern — run the wrapper
  env with a top-level try/catch include, then inspect
  rotor/wakes/pfield/solver for nonfinite (all snippets reproducible from
  this file's evidence; the include-in-catch trick gives you the live
  state post-mortem).

## After the fix

1. Rerun the leg smoke to PASS (all 7 stages); also rerun
   `test/runtests_unit_solver.jl` if src changed (test-runner, ≤4
   threads). Record results in the provenance smoke section.
2. Commit. **Refresh the pin**: delete + recreate the annotated tag
   `campaign/p021-fgs-warmstart-20260924` at the new FLOWPanel commit
   (tags are local-only; origin push still owed after `gh auth login`).
   FastMultipole/FLOWVPM pins unchanged unless you touched them.
3. **Refresh the orc deployment** (rsync mode, Ryan 2026-09-22 ruling):
   re-export `git archive <tag>` locally, rsync the FLOWPanel tree diff to
   `/home/rander39/campaigns/p021-fgs-warmstart-20260924/FLOWPanel.jl/`
   (beware: macOS tar folds git-tracked `__MACOSX/._*` files into xattr
   headers — the original deploy needed an rsync repair pass; rsync the
   whole export, don't re-tar), regenerate
   `MANIFEST.FLOWPanel.jl.sha256` (find -type f | sort | xargs shasum -a
   256, run `sha256sum --quiet -c` on orc to VERIFY), update the
   FLOWPanel `sha` + `content_manifest_sha256` in the campaign root's
   `pins.toml`, and update the pins table + deployment section in the
   provenance file. The campaign env
   (`…/env`, dev-pointed at the deploy trees) needs no change for
   src-only edits.
4. Update the staged sbatch lines in the provenance if any env names
   changed. Job 2 (`run_r4_fgs_cold_newdefault.slurm.sh`) is untouched by
   the leg work.
5. Clean up smoke VTK from the repo: `data/fgs_wsr4_R1_ckpt_*` (local
   checkout data dir — smoke artifacts, not tracked).

## Ryan's gates at submission (surface, don't decide)

Staged sbatch lines live in the provenance "Staged submissions" section
(orc login shell, from the deployed FLOWPanel tree). Open decision points:
1. `FGS_TOL_ABS` — staged with the champion cold tolerance
   3.4309419310610173e-7 (found j-independent across the 13777133 ladder);
   per-step certified BC check is the binding gate. Alternative:
   re-staircase first.
2. Walltime — `--time=36:00:00`, estimate-based (NO prior R4 unsteady data
   exists anywhere; verified on-cluster). 648 simulated steps total under
   the leg layout.
3. Optional arm `ilu_nfcache_proj1` — keep or strike.
4. Sequential-legs-on-one-exclusive-node interpretation of "arms as
   separate Slurm tasks" (chosen for constant hardware; dagedge
   precedent) — veto means array-job rework.
Apply knobs staged: KNOBS_P=12 KNOBS_MAC=0.55 KNOBS_LEAF=48 (certified
budget-500 row of the j64 13777133 tune CSV; the original spec's
"P=15/MAC=0.55" was the budget-0 row — already flagged in provenance).

## Standing gates (carried; surface, don't act)

- Notebook entry (FGS Stages 1+2, gate-0, dagedge campaign + verdict,
  default adoption, this campaign): offer once, verbosity per Ryan.
- Origin pushes: branches + tags `campaign/p021-fgs-stage2-20260923`,
  `campaign/p021-fgs-dagedge-20260924`, `campaign/p021-fgs-warmstart-20260924`
  (three repos) after `gh auth login -h github.com`.
- dagedge jobs 13879622 (perf) / 13879625 (profile, afterany) may still be
  queued/running on m12 — this campaign queues behind them; their harvest
  is a SEPARATE task (`fgs_dagedge_benchmark_provenance_20260924.md`).
- 018 NT-ladder jobs run concurrently — disk alarms are their VTK; launch
  hpc-storage, don't touch their queue.

## Traps

- Cold = zero-initial-guess, warm = seeded guess, SAME process (Ryan
  2026-09-23). Pre-sim setup costs excluded from per-step comparisons,
  reported in the setup columns (never amortized).
- The winB sentinel currently reads ok even when every step is
  unconverged — judge legs by the CSV `solved` column and the driver's
  n_unconverged warning, not the sentinel alone (and consider fixing).
- Tolerances are per rung AND environment; f32full is certified cold-R4
  ONLY (campaign FGS runs f64; Job 2 runs f32full where certified).
- Task logs are output-buffered — judge liveness by CPU/outputs.
- Local runs ≤4 threads; laptop `runtests_benchmark_cold.jl` failure is
  pre-existing (BLAS pin) — don't chase.
- Never edit the deployed orc trees in place except via the
  manifest-verified rsync refresh above; never submit from
  `~/projects` live clones (no-silos rule).
