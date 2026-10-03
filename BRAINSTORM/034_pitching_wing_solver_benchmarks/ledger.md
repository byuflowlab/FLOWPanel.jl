# 034 ledger — data of record

Append-only, dated sections. Tables, certified CSV pointers, frozen settings, and
provenance pins land here. Nothing yet — item STAGED 2026-10-02, no runs.

## 2026-10-02 — item opened

- Scope rulings (Ryan): core cost + warmstart; wake `:panel`. See control doc
  Decision log.
- No data.

## 2026-10-02 — Phase 0: harness adaptation (local smoke, NOT publishable)

### Provisional mesh ladder (freezes at END of Phase 1)

Knob triples preserve the example default's proportions (13/161/9); cells from
`2*n_sec*n_span + 4*(n_chord-2)*(n_endcap-1)` (`examples/pitching_wing.jl:93-95`),
verified against the mesh formula (`n_sec = 2*max(11,cld(n_airfoil,2))-2`,
`n_chord = cld(n_sec,2)+1`).

| rung | n_span | n_airfoil | n_endcap | cells |
| --- | --- | --- | --- | --- |
| R1 | 7 | 89 | 5 | 1920 |
| R2 (example default) | 13 | 161 | 9 | 6688 |
| R3 | 19 | 233 | 13 | 14336 |
| R4 | 27 | 337 | 19 | 30168 |

### LHS-rebuild measurement (headline question) — ANSWERED: reused, never rebuilt

Driver `benchmark/p034_phase0_lhs_rebuild.jl`; CSVs =
`data/phase0_lhs_rebuild/{summary,per_step}.csv` + `banner.txt`. Rung R1
(1920 cells), 20 steps, panel wake, default FMM backend, `-t 1`.

| measurement | value |
| --- | --- |
| max per-step \|G − G_construction\| over 20 steps + post-final | 0.0 (byte-identical) |
| `solver.Glu` objectid stable (no refactorization) | true |
| rotation invariance: raw G from end-of-run rotated geometry vs t=0, max-abs | 1.20e-13 |
| rotation invariance, rel Frobenius | 2.94e-14 |
| counterfactual per-step rebuild cost (min-of-5): G assembly | 1.259 s |
| counterfactual per-step rebuild cost (min-of-5): LU | 0.0393 s |
| constructor-recorded g_assembly_s / lu_s (first call) | 2.390 / 0.0652 s |
| simulate! total / per step (context incl. compilation) | 107.8 s / 5.67 s |

Code evidence: `Backslash` assembles+LUs G once in its constructor
(`src/FLOWPanel_solver.jl:489-497`, `lu!` aliases G); per-step `_solve!` only
refreshes the RHS unless `update_G=true` (`:1714-1736`), which nothing in the
simulate! path ever passes. Per step `propagate_kinematics!` rigidly rotates
nodes+Das (`src/FLOWPanel_simulate.jl:1495`) and `transform_body_solvers!`
(`:1505`) mirrors the delta into persistent solver state — documented no-op for
Backslash (dense Dirichlet operator rotation-invariant,
`src/FLOWPanel_solver.jl:1248-1249`), `transform_plan!` for persistent-plan
Krylov, `transform_solver!` trees for FGS (021 rigid_motion_tree_reuse item).
The newest-wake-row coupling enters through the RHS (wake influence →
control-point velocity/potential before the solve), not the LHS.

### Dirichlet confirmation — CONFIRMED, bc_error! applies as-is

`build_pitching_wing_body` constructs
`RigidWakeBody{Union{ConstantSource,ConstantDoublet},2,Float64,true}`
(`examples/pitching_wing.jl:242-243`); the 4th parameter is DBC
(`has_dirichlet_bc`, `src/FLOWPanel_abstractbody.jl:156`), identical for the
static (semiinfinite) and unsteady bodies. Asserted by the existing suite
(`test/runtests_example_pitching_wing.jl:271`) and re-asserted at runtime by
`bc_error!` itself (`benchmark/common.jl:382`). No decision_rules amendment
needed.

### Solver injection

`solver_factory` kwarg (default `pnl.Backslash`, bit-identical default path) on
`prepare_pitching_wing` and `run_pitching_wing_static_polar`, threaded to both
hardcoded sites (`:1008` unsteady via the returned `sim.solver`; `:867` static
polar builds `solver_factory(body)` per alpha). Verified: full
`test/runtests_example_pitching_wing.jl` green (1450/1450 pass).

### Local-environment caveat (recorded for every Phase 0 CSV)

macOS OpenBLAS 0.3.31: `BLAS.get_num_threads()` reports 8 after any GEMM
regardless of pin, and a 2000x2000 GEMM times identically at 1 vs 4 threads —
the pin is neither observable nor effective locally. Phase 0 smokes launch with
`BENCH_BLAS_THREADS=8` so `assert_and_banner`'s assert passes and the banner
records the truth; the strict single-mode pin is untouched for HPC/published
runs.

### Availability smoke (021 W1–W6 analog) — 4/4 PASS

Driver `benchmark/p034_avail_smoke.jl`; CSVs =
`data/phase0_avail_smoke/{summary,steps}.csv` + `banner.txt`. Rung R1
(1920 cells), 10 unsteady steps per arm, panel wake, `VelocityThroughSources`,
one process (cold = zero-initial-guess), per-step certified `bc_error!` against
the fixed t=0 scale rms_b_t0 = 1.0682. `-t 1`, BENCH_BLAS_THREADS=8 (see
caveat above). Timings are availability context only, NOT benchmarks
(compilation included in each arm's first step).

| arm | completed | all cert | max bcerr_rel | final CL | t_setup s | mem state MB |
| --- | --- | --- | --- | --- | --- | --- |
| backslash | true | true | 2.42e-9 | 0.2923276 | 3.16 | 29.6 |
| krylov_gmres | true | true | 8.53e-9 | 0.2923276 | 0.038 | 0.98 |
| krylov_ilu_nfcache | true | true | 2.63e-9 | 0.2923276 | 1.57 | 25.0 |
| fgs (seed knobs) | true | true | 6.19e-6 | 0.2923228 | 2.88 | 42.9 |

Notes:
- CL identity agrees across arms to ~5 digits; fgs's offset in the 5th digit is
  consistent with its BC floor.
- fgs seed knobs (p=4/mac=0.5/leaf=50/inner=2, tol_abs=1e-6·rms_b_t0, f64,
  dagteam+backoff defaults) converge in 36–41 sweeps/step but carry an
  apply-accuracy BC floor ~3–6e-6 rel, GROWING with wake rows (3.1e-6 step 0 →
  6.2e-6 step 9) — exactly the floor phenomenon 021's decision rules warn about;
  per-rung FMM/τ tuning is Phase 1's job. Availability unaffected (metric
  certified; value reported, not thresholded).
- krylov arms hit 1e-8-level BC as promised (rtol=1e-8); unpreconditioned gmres
  is ~8× slower per step than ilu_nfcache already at this rung.

## 2026-10-02 — Phase 1: FGS knob retune on R1 (local smoke, NOT publishable)

Ryan approved Phase 0 + gave the Phase 1 go-ahead 2026-10-02 (control doc
decision log); Phase 0 committed as `1d235d7` (code/harness) + `911a82f`
(records).

Driver `benchmark/p034_phase1_fgs_tune.jl`; CSVs =
`data/phase1_fgs_tune/{summary,steps}.csv` + `banner.txt`. Rung R1 (1920
cells), 17-step unsteady march per config (enough to expose the wake-row
growth trend), panel wake, `-t 1`, BENCH_BLAS_THREADS=8 (macOS caveat above).
Grid: Phase 0 seed as control + 021 tau=1e-6 rotor-rung winners
(fgstune_verify.csv: R1 6/0.3/150/5, R2 8/0.4/100/10) + FGSSolver constructor
defaults, with a tol_factor dimension (tol_abs = tolf*1e-6*rms_b_t0).
Gate: max-over-steps bcerr_rel <= 1e-6, every pass certified.

| config (p/mac/leaf/inner/tolf) | max bcerr_rel | med t_solve s | med niter | meets |
| --- | --- | --- | --- | --- |
| 4/0.5/50/2/1.0 (Phase 0 seed, CONTROL) | 7.10e-6 | 1.290 | 40 | NO |
| 7/0.4/10/2/1.0 (ctor defaults) | 1.88e-7 | 1.398 | 39 | yes |
| **6/0.3/150/5/1.0 (021 R1 seed) — WINNER** | **9.34e-8** | **1.265** | **15** | **yes** |
| 8/0.4/100/10/1.0 (021 R2 seed) | 1.49e-7 | 1.261 | 8 | yes |
| 7/0.4/10/2/0.3 | 5.74e-8 | 1.427 | 44 | yes |
| 6/0.3/150/5/0.3 | 3.03e-8 | 1.278 | 17 | yes |
| 7/0.4/10/2/0.1 | 1.90e-8 | 1.461 | 48 | yes |

Winner rationale: per-step time statistically tied with 8/0.4/100/10 (1.265 vs
1.261 s); 6/0.3/150/5 has the larger BC margin (10.7x vs 6.7x) and finer sweep
granularity (inner=5 vs 10 — less overshoot when warmstart shrinks the work in
Phase 3). Per-step bcerr trend FLAT over 17 steps for both finalists (no
wake-row growth; the control's growth phenomenon is absent once apply accuracy
is adequate). R1 FGS setting (PROVISIONAL until the Phase 1 freeze):
p=6, mac=0.3, leaf=150, inner=5, tol_abs=1e-6*rms_b_t0, rlx=1.0, shrink=true,
dagteam+backoff defaults, f64. Per-rung confirmation on R2 pending.

## 2026-10-02 — Phase 1: consistency-driver smoke R1 + sizing probe R2 (local, NOT publishable)

Driver `benchmark/p034_phase1_consistency.jl` (per-step certified bc_error!
gate metric vs fixed t=0 scale, 021 Phase 3 arm-promise absolute stats, CL/CM
identity columns, env rung/arms/cycles/FGS knobs). Local smoke CSVs in session
scratchpad (not retained — campaign CSVs are the record); headline numbers:

R1, 17 steps, all four arms, `-t 1`/BENCH_BLAS_THREADS=8:

| arm | max bcerr_rel | promise viol | final CL | t_sim s |
| --- | --- | --- | --- | --- |
| backslash | 3.37e-9 | 0 | 0.33359571 | 71.1 |
| krylov_gmres | 9.89e-9 | 0 | 0.33359577 | 448.3 |
| krylov_ilu_nfcache | 3.62e-9 | 0 | 0.33359571 | 54.0 |
| fgs (6/0.3/150/5) | 9.34e-8 | 0 | 0.33359203 | 54.9 |

4/4 meet the gate on this smoke; CL agrees to 7 digits across the <=1e-8 arms,
fgs offset 1.1e-5 relative (consistent with its BC level). R2 probe
(backslash+fgs, 5 steps): both meet gate; **R1-tuned FGS knobs HOLD at R2**
(max bcerr_rel 9.09e-8); fgs rel_MAX approaches 1e-6 at R2 (9.8e-7, reported —
gate metric is rel L2 per decision_rules). Sizing at `-t 1`: R2 ~27 s/step
total (fgs solve 11.8 s) → full 3-cycle 495-step march ≈ 3.7 h/arm; gmres
R1 is ~26 s/step → R2 est. 15–28 h. HPC campaign required for the 2-rung
certified record.

## 2026-10-02 — Phase 1 campaign pins (campaign/p034-phase1-20261002)

Scope (Ryan 2026-10-02): rungs R1+R2, all four arms, 3-cycle marches,
single-thread mode, non-exclusive allocations. One (rung, arm) per job via
`benchmark/slurm/p034_phase1.sh`; outputs to the consolidated data root
`~/projects/FLOWPanel.jl/data/p034_phase1/R{1,2}/` via the worktree data
symlink. Judge from the CSVs there.

| repo | tag | commit | worktree |
| --- | --- | --- | --- |
| FLOWPanel.jl | campaign/p034-phase1-20261002 | 43f9763 | ~/campaigns/p034-phase1-20261002/FLOWPanel.jl (HEAD d9e4432 = tag + data-symlink commit; clean) |
| FastMultipole | campaign/p034-phase1-20261002 | 3da58a1a | ~/campaigns/p034-phase1-20261002/FastMultipole |
| FLOWVPM.jl | campaign/p034-phase1-20261002 | d896145 | ~/campaigns/p034-phase1-20261002/FLOWVPM.jl |

Campaign env = the FLOWPanel worktree's own `--project=.`; Manifest dev-paths
are RELATIVE (`../FastMultipole`, `../FLOWVPM.jl`) and resolve to the pinned
worktrees above; tree verified clean after resolve. FLOWVPM pin d896145
satisfies the new-merge-law floor (>= 8d4a3b4). Tag pushed to the orc clone;
**origin (GitHub) tag push PENDING — local gh auth token invalid** (Ryan:
`gh auth login -h github.com`, then `git push origin campaign/p034-phase1-20261002`).

## 2026-10-02 — Phase 1 campaign submitted (8 jobs, m9 --qos=normal, non-exclusive)

Submitted from the campaign worktree (pins above); 1 CPU + 24 G per job,
`-t 1` strict single mode. First job gated the precompile (flock guard +
serialized warm; "precompile stage done" confirmed before fan-out). Walltimes
sized from local `-t 1` estimates x broadwell margin. Outputs →
`~/projects/FLOWPanel.jl/data/p034_phase1/R{1,2}/` (judge from CSVs:
`summary_<arm>.csv`, `steps_<arm>.csv`, `banner_<arm>.txt` per arm).

| job | rung | arm | walltime |
| --- | --- | --- | --- |
| 13961478 | R1 | backslash | 12 h |
| 13961567 | R1 | krylov_ilu_nfcache | 12 h |
| 13961568 | R1 | fgs | 12 h |
| 13961569 | R1 | krylov_gmres | 24 h |
| 13961570 | R2 | backslash | 36 h |
| 13961571 | R2 | krylov_ilu_nfcache | 36 h |
| 13961572 | R2 | fgs | 36 h |
| 13961573 | R2 | krylov_gmres | 70 h |

Slurm exit status is advisory (judge by outputs, not sacct). No VTK written
(driver runs path=nothing) — no storage pressure expected from this wave.
