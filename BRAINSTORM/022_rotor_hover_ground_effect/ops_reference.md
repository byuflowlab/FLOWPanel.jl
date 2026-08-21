# 022 ops reference

Inherits the 018 ops discipline (`018_dji9443_hover_convergence_campaign/ops_reference.md`);
this file records only what 022 needs and what differs.

## Where things run

- Production: BYU cluster `orc`, checkout `/home/rander39/projects/FLOWPanel.jl`.
  `ssh orc` needs a live ControlMaster socket (2FA otherwise) — ask Ryan if dead.
- Local (Mac): coarse smokes only (RHPC_MESH=40_40, NT=6, ≤4 threads), harvest,
  analysis.
- **Manifest stays julia 1.11.7.** Never re-resolve (021's 1.12 re-resolve
  killed the first 018 F1 trio).

## Submission

```bash
mkdir -p logs/slurm
sbatch --job-name=fp-022-oge-fine                     examples/run_rotor_ground_effect_hpc.slurm.sh p022_oge_fine
sbatch --job-name=fp-022-oge-coarse --time=24:00:00   examples/run_rotor_ground_effect_hpc.slurm.sh p022_oge_coarse
sbatch --job-name=fp-022-ige-fine                     examples/run_rotor_ground_effect_hpc.slurm.sh p022_ige_fine
sbatch --job-name=fp-022-ige-coarse --time=24:00:00   examples/run_rotor_ground_effect_hpc.slurm.sh p022_ige_coarse
```

64 cores / 64G / 72 h (fine), 24 h (coarse). Check `squeue` first — respect
the ≤20 concurrent budget and the live 018 jobs.

## Deploy

`scripts/deploy_phase02e.sh` pattern: scp changed files, **md5-verify every
transferred file** (018 lesson: uncommitted local src drift vs cluster
produced degenerate results). 022 transfers:
- `examples/rotor_hover_ground_effect.jl`
- `examples/run_rotor_ground_effect_hpc.slurm.sh`
No `src/` changes in 022 staging (the GS instrumentation is driver-local by
design) — if `src/` ever diverges, STOP and reconcile with the live 018 state
first.

## Post-submission (mandatory before declaring healthy)

1. Wait for the job to enter R and the log to show the driver banner.
2. Verify banner: mesh key/file, RPM=6000 rho=1.16 R=0.1195, shedding trace
   counts + zero circumferential edges, ground line (h/R, disc, panels,
   policy, GS settings), staged-startup schedule. **Exit 0 ≠ health; banner
   verification is mandatory** (018 override-precedence incident).
3. IGE runs: spot-check early `Ground step N:` lines — tangency RMS bounded,
   GS not warning about cap hits.
4. Record job IDs in `ledger.md` and the item's in-flight table.

## Harvest

- `data/<tag>/<tag>_case_metadata.toml` → CT_cycle_mean/std, gs_iters_max,
  gs_nonconverged, ground_tangency_rms_max, below_ground_count_max.
- `<tag>_CT_per_rev.csv` → M1 window per `decision_rules.md`.
- `<tag>_gs_convergence.csv`, `<tag>_ground_diagnostics.csv` → health gates.
- Disk: VTK sweeper policy as in 018 (`scripts/p018_vtk_sweeper.sh` pattern);
  the ground body doubles the per-step VTK file count — watch the ~200 G
  budget.

## Known constraints

- Warm-start with GROUND_ENABLE=true errors out in the driver (untested
  multi-body restart) — extend runs by re-running longer cold, or chain OGE
  only.
- Wipe policy: launcher `rm -rf data/$RUN_NAME` on resubmit;
  `RHPC_KEEP_PREV=true` to preserve a failed attempt for inspection.
- The ground body requires `BERNOULLI_ONLY=true`: `PressureLaplace` is built
  for a single body and throws `ArgumentError: PressureLaplace was constructed
  for 1 bodies, got 2 bodies` on the two-body IGE path. The launcher's shared
  block sets it; hand-rolled smokes must too.

## Measured step rates (2026-08-19)

From the live Phase 1/2 jobs, for sizing future walltimes:

| rung | mesh | s/step (avg) | s/step (recent) | 1007-step schedule |
|---|---|---|---|---|
| coarse | 56_57 | ~87 | ~110 | **25–31 h** |
| fine | 45_185_ct4 | ~130 | ~170–250 | ~40–55 h |

**The 24 h coarse walltime used in Phases 1–2 is undersized** — p022_oge_coarse
(13207680) lands around rev 25 of 28, short of a ≥10-rev settled window at
`SETTLE_REVS=20`. Combined with the no-warm-start constraint above (a walled
IGE run cannot be chained, only re-run cold), first-submission walltime is
load-bearing. Phase 5 sizes its coarse runs at **48 h**.

Also note `squeue`/`scancel`/`sbatch` are not on `PATH` in a non-interactive
`ssh orc` shell (they are shell functions in the login profile). Use
`/apps/slurm/latest/bin/<cmd>` directly for scripted checks.

## Harvest without a finished run

`<tag>_CT_vs_rev.csv` is only written after the time loop completes, so a
walled, cancelled, or diverged run leaves none. `scripts/p022_harvest.py` falls
back to `monitors/<tag>_monitor02_force_system1.csv` (written incrementally):
`CT = -CFx`, `revolution = time*RPM/60`, `step_csv = step_monitor + 1`. Use it,
not `scripts/p018_analyze.py`, which pins RPM 5400 in that same fallback with
no CLI pass-through and would misplace every 022 rev boundary by 11%.
