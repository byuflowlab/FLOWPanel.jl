# Provenance: 021 warm-start R4 head-to-head + cold rerun under the new FGS default (2026-09-24)

Campaign spec: `fgs_warmstart_r4_reset_prompt_20260924.md` (Jobs 1 and 2).
**Status: DRAFT — pins unfilled, submission NOT approved. Both sbatch lines
are staged for one Ryan approval; nothing runs until his explicit go.**

## Decision rule (Ryan 2026-09-24)

If warm FGS does not compete with warm krylov_ilu_nfcache on per-step
time-to-target at R4 (pre-sim setup excluded, reported separately), the FGS
acceleration thread goes back to the drawing board. Ryan rules on the
harvested table: per-step medians, spreads, and cumulative-cost crossover
curves are presented, not pre-judged.

## Step 0 (completed in-session, Ryan-approved 2026-09-24)

- FLOWPanel `08667c2` on `fastmultipole`: FGSSolver + FGSPreconditioner
  defaults `sweep_order=:dagteam`, `dagteam_idle=:backoff`;
  `dagteam_precision` stays `:f64`; kernel-level `t_project` timing on both
  solver families (+ `step_t_project`); unit solver suite 508/508.
- FastMultipole `745af760` on `flowpanel-20260817`: dagteam plan BoundsError
  on empty direct lists (single-leaf trees) fixed + regression testset —
  found by the default flip, required for the default to be usable.
- Thread-scalability thread PARKED (ledger + dagedge results note).

## Job 1 — warm-start head-to-head

- **Driver**: `benchmark/fgs_r4_warmstart_ab.jl` (per-arm entry) over
  `benchmark/rotor_hover_solver_unsteady.jl` (the ruling-12 Phase-3 vehicle),
  extended this session with: `krylov_ilu_nfcache` config (persistent_plan +
  cache_nearfield + one-time ILU), setup-cost split + pre-march priming solve
  (state-restored) so the plan/nfcache build lands in `t_setup_prime`, never
  in step 1; per-step `t_project`; per-step solved-strength snapshots
  (`*_strength_snapshots.bin`); explicit-knob env overrides;
  `OUTDIR_OVERRIDE` (data-root discipline).
- **Launcher**: `benchmark/run_r4_fgs_warmstart.slurm.sh` — ONE exclusive
  128-core zen3 node, 14 legs SEQUENTIAL (2 ckpt → 5 winA → 7 winB, winB
  gated on its family checkpoint's STATUS), one fresh Julia process per
  (arm, leg), j=64, champion placement
  (`--interleave=0-3 --cpunodebind=0-3`), BLAS=1 both families (ILU never
  BLAS-swept). *Interpretation note for Ryan: the
  reset prompt says "arms as separate Slurm tasks"; sequential arms on one
  exclusive node was chosen to keep hardware constant across the
  cross-arm comparison (dagedge precedent). Veto if separate array jobs
  (possibly different nodes) are preferred.*
- **Fixture**: R4 (58,192 panels, production mesh prescription), unsteady
  wake-on hover via `simulate!`, RHPC frozen settings, NT=36, Bernoulli CT
  monitor on.
- **Checkpoint + restart layout (Ryan 2026-09-24, second ruling)**: TWO
  restart checkpoints, one per solver family. Legs per (arm, leg) process:
  - `ckpt` (2): each family's COLD arm marches revs 1–3 (108 steps) from
    scratch with `SAVE_VTK=true` — doubles as that cold arm's Window A and
    as the family's shared restart source (VTK to the shared data root via
    the deploy tree's `data` symlink; restartable under standard retention).
  - `winA` (5 warm arms): rev 1 (36 steps) from scratch, VTK off.
  - `winB` (all 7 arms, cold included for symmetric treatment): rev 4
    restarted from the FAMILY checkpoint at step 108 via
    `simulate_warmstart!`. Revs 2–3 are thus simulated once per family, and
    future warm work can restart from the same checkpoints without
    re-marching.
- **Windows (harvest-side, BOTH including transients, SKIP_STEPS=0)**:
  A = steps 1–36 (startup revolution from the very first step, non-restarted
  legs); B = steps 109–144 (the winB restart legs).
  **Reporting note (Ryan 2026-09-24)**: solver warm-start histories are not
  serialized in the checkpoint, so each restarted WARM leg's first (order+1)
  steps are effectively cold — this history-fill transient sits inside
  Window B's transient-included stats by design and is flagged in the
  harvest, never excluded. Within a family, Window B arms share one
  checkpoint, so B isolates the initial guess; across families the two
  checkpoints' wakes differ at the known fixed-point level.
- **Arms (7)**: fgs_cold, fgs_prev, fgs_proj1, fgs_proj2 (order sweep {1,2}),
  ilu_nfcache_cold, ilu_nfcache_prev, ilu_nfcache_proj1 (optional — Ryan may
  strike; uses the SAME shared extrapolation coefficients as FGS).
  Cold = zero-initial-guess solves in the same process (Ryan 2026-09-23).
- **Knobs**: FGS = R4 champion P8/MAC0.4/leaf100/inner3
  (`benchmark/retained_r4_champion.toml`), **f64 both sides** (f32full is
  certified cold-R4 only; not carried to the warm fixture). ILU/nfcache =
  certified 2026-09-22 apply knobs passed explicitly (KNOBS_P/MAC/LEAF, copied at staging from the budget-500 row of the R4
  `tune_phase2.csv` in the 13777133 run dirs); ILU build knobs leaf10/MAC1.0/8192·n (the certified
  configuration); `NFCACHE_MAX_GIB=500`.
- **Convergence contract**: fixed-accuracy. **Stopping-semantics finding
  (2026-09-24, caught by the R1 smoke)**: Krylov.jl's `rtol` is relative to
  the INITIAL residual r0 = b − A·x0, so a warm x0 shrinks the target with
  it — identical iteration counts cold vs warm, deeper solve, zero
  time-to-target gain by construction. Fixed with a new
  `KrylovSolver(rtol_rhs=...)` mode: stop at ‖r‖ ≤ atol + rtol_rhs·‖b‖
  (absolute target anchored to the current rhs; numerically identical for
  cold solves where r0 = b). The driver enables it for PHASE=phase3* only,
  so Phase-2 rows keep their historical semantics. Krylov rtol_rhs=1e-6
  (atol=1e-14);
  FGS stops at an explicit absolute tolerance `FGS_TOL_ABS` — **value to be
  fixed at staging (open item below)**. The binding per-step check is the
  arm's own promise measured independently: one certified `bc_error!` pass
  per step (`bcerr_max <= bcerr_tol`, instrument 10× sharper by
  construction, `bcerr_certified` flags failures loudly).
- **Metrics per step (CSV row)**: step, t, t_step_total/t_step_net, t_solve,
  niter (legacy), niter_first (headline), nsolves (must be 1), t_project,
  solved, bcerr_* order statistics + certification, n_particles, CT.
  Setup (once per arm, separate columns, never amortized): t_setup (solver
  construction incl. ILU build), t_setup_prime (plan + nfcache build via the
  priming solve), setup_detail (ILU tree/lists/assembly/factorization split).
- **ILU persistence design check (owed, answered 2026-09-24)**: the ILU
  factorization is built once at `ILUPreconditioner` construction and never
  refactorized; `transform_solver_geometry!` transforms only the plan
  (+nfcache). For this Dirichlet body the operator is exactly invariant
  under rigid motion, so the factors remain valid by the same argument as
  the nfcache. There are NO per-step factorization events; the build lands
  in `setup_detail.ilu_factorization`. The driver asserts the Dirichlet
  context (`has_dirichlet_bc`) before the per-step BC instrument.
- **Accuracy guard (binding)**: the known ~2e-3 FGS-vs-Krylov wake-on
  fixed-point discrepancy is REPORTED, not chased: per-step solved-strength
  snapshots per arm → harvest computes cross-arm per-step solution deltas;
  CT traces reported alongside timing. A speed win at a different fixed
  point is not a clean win.
- **Local smoke**: `benchmark/run_r4_fgs_warmstart_smoke.sh` — all 7 arms,
  R1, 9 steps, 4 threads, functional knobs: **PASS 2026-09-24** (all
  STATUS/COMPLETED sentinels, CSV schema, snapshots, setup split, harvest
  script exercised on the output). Warm-start machinery demonstrably works:
  FGS niter_first 14 (cold, flat) → 13→7 (prev) / 14→7 (proj1/2, engaging at
  step order+2); Krylov (post-rtol_rhs) 8 flat cold → 7→5 warm; t_project
  ~30–70 µs/step. Three defects found AND fixed by the smoke before any HPC
  time: (1) Dirichlet priming solve had a zero rhs (rhs = −potential, not
  velocity) → 0 iterations, no plan build — primes under a synthetic unit
  potential now; (2) the Krylov r0-relative stopping (rtol_rhs finding,
  above); (3) `_effective_tol` read `body.potential` AFTER the solve
  overwrote it, collapsing the recorded Krylov promise to atol=1e-14 and
  flagging every row — now reads the solver's own untouched `rhs` vector.

## Job 2 — cold R4 FGS under the new default

- **Launcher**: `benchmark/run_r4_fgs_cold_newdefault.slurm.sh`, array over
  j ∈ {1,8,16,32,64}, FGS only (NO ILU reruns), new run dirs
  `fgs-cold-newdefault-j<J>-<jobid>` alongside the 13777133 ladder.
- Reuses the per-j staircase-calibrated dagteam configs from
  `thread-scaling-j<J>-13777133/fgs-calibrate/results/dagteam_selected.toml`
  (champion knobs P8/MAC0.4/leaf100, f32full — cold R4/zen3 is exactly what
  was certified; the launcher asserts the knob set) with
  `dagteam_idle=backoff` injected. Justification for carrying the calibrated
  tolerance: idle policy changes scheduling only, never arithmetic
  (FastMultipole solve_dagteam.jl contract; dagedge campaign precedent,
  cross-executor delta 2.24e-7 vs 1e-5 tripwire), and every trial is still
  independently gated by the certified evaluator (STAGE=verify), so a
  transfer failure is caught, not absorbed.
- Measurement contract identical to the existing rows: cold isolated solves,
  COLD_PREPARED_ONLY=1, t_solve minimum, champion placement, BLAS=1.
- **Deliverable**: the 2026-09-22 R4 table re-issued with an
  "FGS-dagteam+backoff (new default)" column (same results file as the warm
  harvest). Expected ≈3.2–3.3 s @ j64, plateau moving j16 → j64.

## Pins (annotated tags created 2026-09-24, all three repos)

Tag (all repos): `campaign/p021-fgs-warmstart-20260924`

| repo | branch | SHA | content |
|---|---|---|---|
| FLOWPanel.jl | `fastmultipole` | `f88e35afcd76eadd9eee45ac25d03fcaa2506594` | campaign harness on top of `08667c2` (dagteam+backoff default + t_project) |
| FastMultipole | `flowpanel-20260817` | `745af76022ea782146901fd65fc301f8474cbc8c` | dagteam empty-direct-list fix on Stage-2 pin `053c8de7` (+dagedge `90a60cc3`) |
| FLOWVPM.jl | `flowpanel` | `8d4a3b4d3012c42fc7d078629234c105b1e570f7` | unchanged production pin (new-merge-law default) |

## Deployment (completed 2026-09-24; submission still owed to Ryan)

- [x] Tagged triple deployed via `git archive <tag>` + rsync repair of 25
  git-tracked `__MACOSX/._*` files that macOS tar folded into xattr headers,
  into `/home/rander39/campaigns/p021-fgs-warmstart-20260924/`
  (ARCHIVER_SKIP marked). All three trees `sha256sum -c` VERIFIED on orc
  against `MANIFEST.<name>.sha256`.
- [x] `pins.toml` written at the campaign root (deployment="rsync", tag +
  SHA + manifest hash per repo, hashes:
  FLOWPanel `5d1c5fb0…`, FastMultipole `b5d11503…`, FLOWVPM `b4d0fa8d…`).
- [x] Campaign env `…/env` cloned from the dagedge campaign env (same dep
  versions; no Project.toml changes since) and `Pkg.develop`-repointed at
  the three deploy trees — Manifest paths verified.
- [x] `data` symlink → `/home/rander39/projects/FLOWPanel.jl/data` in the
  FLOWPanel deploy tree; Das arc table confirmed reachable;
  `logs/slurm/` pre-created.
- [x] Cluster facts (hpc-monitor 2026-09-24): all five
  `thread-scaling-j{1,8,16,32,64}-13777133` dagteam_selected.toml exist,
  P8/MAC0.4/leaf100/dagteam confirmed, tolerance identical at every j
  (3.4309419310610173e-7 — j-independent on this fixture/env). Certified
  j64 apply-knob rows: **budget-500 = P12/MAC0.55/leaf48**, budget-0 =
  P15/MAC0.55/leaf21 (the reset prompt's "P=15/MAC=0.55" parenthetical
  conflated the budget-0 row; nfcache uses the budget-500 row). No prior
  R4 unsteady walltime data exists anywhere in the repo — Job 1 walltime is
  an estimate (see below). /home usage 79 G, far under cap.
- [x] Availability probed 2026-09-25T04:39Z: m12 access=normal
  (qos normal/test, maxtime 3-00:00:00), 0 idle / 41 mixed / 87 alloc —
  submittable, expect queue wait for an exclusive node (dagedge jobs
  13879622/25 are in the same queue).

## Staged submissions (RYAN'S GO REQUIRED — nothing submitted)

From an orc login shell:

```bash
DEPLOY=/home/rander39/campaigns/p021-fgs-warmstart-20260924
cd $DEPLOY/FLOWPanel.jl

# ---- Job 1: warm-start head-to-head (one exclusive zen3 node, 7 arms) ----
export WSR4_PROJECT=$DEPLOY/env
export CAMPAIGN_PINS=$DEPLOY/pins.toml
export WSR4_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910
export FGS_TOL_ABS=3.4309419310610173e-7   # open item 1: recommended value
export KNOBS_P=12 KNOBS_MAC=0.55 KNOBS_LEAF=48   # certified budget-500 row
export CONTENT_MANIFESTS="$DEPLOY/FLOWPanel.jl:$DEPLOY/MANIFEST.FLOWPanel.jl.sha256 $DEPLOY/FastMultipole:$DEPLOY/MANIFEST.FastMultipole.sha256 $DEPLOY/FLOWVPM.jl:$DEPLOY/MANIFEST.FLOWVPM.jl.sha256"
sbatch -p m12 --export=ALL benchmark/run_r4_fgs_warmstart.slurm.sh

# ---- Job 2: cold R4 FGS under the new default (array j in {1,8,16,32,64}) ----
export COLD_PROJECT=$DEPLOY/env
export CAMPAIGN_PINS=$DEPLOY/pins.toml
export COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910
sbatch -p m12 --export=ALL benchmark/run_r4_fgs_cold_newdefault.slurm.sh
```

Walltime basis (no prior R4 unsteady data exists): R4 j64 cold solve is
2.4–3.3 s; per step add wake evolution, the certified BC pass, and
monitors → est. 15–40 s/step. Checkpoint layout totals 648 simulated steps
(2×108 ckpt + 5×36 winA + 7×36 winB, vs 1008 for full marches) ≈ 2.7–7.2 h
plus ~1–2 h of per-process setup across 14 legs; `--time=36:00:00` carries
large margin and fits m12's 3-day cap. Job 2's 8 h/task covers the j=1
rung's slow verify. Resume paths exist in both launchers
(RESUME_FROM_JOB_ID; landed legs skip by STATUS).

Deployment: origin push still deferred (GitHub re-auth owed), so
`deployment = "rsync"` mode per Ryan's 2026-09-22 ruling, exactly as the
dagedge campaign: `git archive <tag> | ssh orc tar -x` into fresh dirs under
`/home/rander39/campaigns/p021-fgs-warmstart-20260924/` (ARCHIVER_SKIP),
sha256 content manifests generated from the same local export and verified
on orc (Job 1's launcher re-verifies via CONTENT_MANIFESTS at job start; the
cold harness verifies Job 2 via CAMPAIGN_PINS), campaign env with Manifest
dev-paths at the deploy trees, `data` symlink to the shared data root in the
FLOWPanel deploy tree (RHPC resolves `data/...` relative paths; launcher
preflights the Das arc table). Under this mode the CSV `commit` columns read
"unknown" (no .git in an archive export) — provenance rests on the pins +
manifests. Outputs to the consolidated data root
(`/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/` unless the
unsteady fixture warrants a new dated root — decide at staging).

## Open items before the staging gate

1. **FGS_TOL_ABS for the warm fixture** — tolerances are per rung AND per
   environment (never carry). Options for Ryan: (a) carry the champion cold
   tolerance 3.4309419310610173e-7 with the per-step certified BC check as
   the binding gate, (b) re-staircase cold on the campaign node first
   (cheap; Job 2's calibration machinery), then use that value. The per-step
   `bcerr_max <= bcerr_tol` column reports compliance either way.
2. **Walltime**: no prior R4 unsteady data exists (verified on-cluster
   2026-09-24), so `--time=36:00:00` is estimate-based (see the staging
   section's basis).
3. Optional arm ilu_nfcache_proj1 — keep or strike.
4. Sequential-legs-on-one-node interpretation (see Job 1 launcher note).

RULED by Ryan 2026-09-24 (second round): checkpoint+restart layout with TWO
family checkpoints (implemented as the ckpt/winA/winB legs above);
history-fill transient in early Window B is fine, noted in reporting;
per-arm from-scratch Window A kept.

## Smoke + harvest

- Local smoke (R1, 4 threads, 8 steps, all 7 arms): _result recorded here
  when complete._
- Harvest deliverable: `fgs_warmstart_r4_results_<date>.md` — per-arm
  Window A/B tables (iterations + time-to-target mean±spread/median), setup
  cost in its own column, per-step cost traces vs step index (transient
  shape is a deliverable), solution-agreement report, compete/no-compete
  recommendation for Ryan's ruling. Notebook entry: offered, not written.
