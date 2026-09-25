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
  128-core zen3 node, arms SEQUENTIAL, one fresh Julia process per arm,
  j=64, champion placement (`--interleave=0-3 --cpunodebind=0-3`), BLAS=1
  both families (ILU never BLAS-swept). *Interpretation note for Ryan: the
  reset prompt says "arms as separate Slurm tasks"; sequential arms on one
  exclusive node was chosen to keep hardware constant across the
  cross-arm comparison (dagedge precedent). Veto if separate array jobs
  (possibly different nodes) are preferred.*
- **Fixture**: R4 (58,192 panels, production mesh prescription), unsteady
  wake-on hover via `simulate!`, RHPC frozen settings, NT=36, 144 steps
  (4 revolutions), no restart, VTK off, Bernoulli CT monitor on.
- **Windows (harvest-side, BOTH including transients, SKIP_STEPS=0)**:
  A = steps 1–36 (startup revolution from the very first step);
  B = steps 109–144 (fourth revolution, from its first step).
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

## Pins (TO FILL before submission — annotated tags, three repos)

| repo | tag | SHA | worktree |
|---|---|---|---|
| FLOWPanel.jl | `campaign/p021-fgs-warmstart-20260924` | _pending_ | _pending_ |
| FastMultipole | `campaign/p021-fgs-warmstart-20260924` | _pending_ | _pending_ |
| FLOWVPM.jl | `campaign/p021-fgs-warmstart-20260924` | _pending_ | _pending_ |

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
2. **Walltime confirmation** against prior 021 unsteady runs at R4 before
   staging (currently `--time=36:00:00` for 7 sequential arms).
3. Optional arm ilu_nfcache_proj1 — keep or strike.
4. Sequential-arms-on-one-node interpretation (see Job 1 launcher note).

## Smoke + harvest

- Local smoke (R1, 4 threads, 8 steps, all 7 arms): _result recorded here
  when complete._
- Harvest deliverable: `fgs_warmstart_r4_results_<date>.md` — per-arm
  Window A/B tables (iterations + time-to-target mean±spread/median), setup
  cost in its own column, per-step cost traces vs step index (transient
  shape is a deliverable), solution-agreement report, compete/no-compete
  recommendation for Ryan's ruling. Notebook entry: offered, not written.
