# 018 GPU — ladders PAUSED 2026-08-30 ~01:15 (context-reset handoff)

Written for the next agent. Read this FIRST. Companions (do not re-derive):
`gpu_nomerge_ladder_status_20260829.md` (3R ladder results, merge-error
hypothesis, recipes/quirks §4, guardrails §5 — ALL still true) and
`phase_17_nprop_nt_ladder.md` last three sections (the NT-climb analysis).

## 1. What is paused and why

Ryan ordered all jobs paused 2026-08-30 ~01:15 for exactly 24 h 2 min.
All five ladder jobs were `scancel`ed (Slurm suspend is not user-level):

| cancelled job | arm | ladder | state at cancel |
|---|---|---|---|
| 13513063 | NT144 λ3.0 | merge-off | RUNNING mgh-1-2, step 808/4319 (~46 min) |
| 13513064 | NT72 λ3.0 | merge-off | RUNNING eng-1-1, step 618/2159 (~45 min) |
| 13513065 | NT36 λ3.0 | merge-off | PENDING eng |
| 13513913 | NT72 λ3.0 rlxf=0.3 | fixed-rlxf (gaussian) | PENDING eng |
| 13513914 | NT144 λ3.0 rlxf=0.3 | fixed-rlxf (gaussian) | PENDING eng |

Both running arms were < 1 h in, so the relaunch is a FRESH restart from
step 0 (no `_s2` chaining; the dispatcher wipe policy clears the partial
`data/p018_*_3r_nm` dirs on resubmit — intended). Before the pause, both
running merge-off arms were verified merge-off (`every = 0` in the
`MergeParticles` block of each run's `.metadata.toml`) and stepping cleanly
(2–6 s/step, no cliff).

## 2. RELAUNCH (the next agent's job, at pause end)

One command; submits all five with the exact original submit lines
(recovered via `sacct ... SubmitLine`), from the correct silos:

    ssh orc 'bash -lc "bash ~/relaunch_018_ladders_20260830.sh"'

Script is already on orc (`~/relaunch_018_ladders_20260830.sh`, syntax-checked).
It ends with an `squeue` so you see the five new job IDs — record them.
After relaunch: re-verify merge-off metadata (`every = 0`) for the three
`_nm` arms once they start, and for the two `rlxf0p3` arms verify the run
banner prints `rlxf:0.3` and `filament_reg:gaussian` (banner is in the
slurm .out, ~20 lines after "DJI 9443 Phase 2e").

## 3. The two ladders (what they test)

- **merge-off λ3.0 ladder** (`*_3r_nm`, MERGE_PARTICLES=false): if CT̄(NT)
  flattens vs the with-merge λ3.0 column (0.0708/0.0718/0.0738), the NT
  climb was merge error (phase_17 hypothesis: NT-dependent merging eats
  NT36's old wake). Also check np vs no-merge estimates (~440–520k @ rev
  30, saturation 520–690k, cap 1.5M) and NT144 late-run step cost (~2×).
- **fixed-rlxf λ3.0 ladder** (`*_3r_rlxf0p3`, gaussian, merge ON): pins
  per-step Pedrizzetti rlxf at 0.3 for every NT instead of the exact-rate
  values (0.16334 @ NT72, 0.08539 @ NT144), testing phase-17 mechanism 3
  (per-step relaxation frequency riding the NT axis). NT36 rung = existing
  completed run 13507289 (CT̄ 0.0707752) — rlxf 0.3 is the NT36 default, so
  only NT72/NT144 have jobs. Readout: within-ladder NT trend vs
  0.0708/0.0718/0.0738.

## 4. Session changes (2026-08-29/30, all uncommitted — ask Ryan before commits)

- `examples/run_dji9443_hover_ct_hpc.slurm.sh`: `RELAX_RLXF_OVERRIDE` hook
  after the case table (the table exports RELAX_RLXF unconditionally, so
  plain `--export` cannot win).
- `examples/run_dji9443_hover_ct_gpu.slurm.sh`: `P018_REPO_OVERRIDE` /
  `P018_PROJECT_OVERRIDE` hooks.
- Both mirrored (`rsync --checksum`) to `~/FLOWPanel-018-gpu-h200` and
  `~/FLOWPanel-018-gpu-gh200` on orc; inert without the env vars.
- **Unused `-lg` LineGauss snapshot silos exist on orc** (built, then Ryan
  switched the rlxf ladder back to gaussian for comparability):
  `~/FLOWPanel-018-lg-h200` (018 layout + local working-tree src/examples ≈
  FLOWPanel-052-h200 tested state), `~/FastMultipole-018-lg-h200`
  (snapshot of FastMultipole-052-h200 — the state that passed the 052d GPU
  LineGauss work; local FastMultipole WIP was deliberately NOT deployed),
  env `~/p018env-lg-h200` (dev-paths repointed; FLOWVPM still the shared
  band-aid tree). Keep them: they are ready if Ryan wants the linegauss
  ladder later (add `FLOWPANEL_FILAMENT_REG=linegauss` +
  `P018_REPO_OVERRIDE=$HOME/FLOWPanel-018-lg-h200` +
  `P018_PROJECT_OVERRIDE=$HOME/p018env-lg-h200`; eng ONLY — no gh200
  LineGauss stack). Note: `FILAMENT_REGULARIZATION[]` applies to VortexRing
  panel-edge filaments AND TE-wake filaments (panels too); wake VPM
  particles stay Gaussian σ regardless.

## 5. Monitoring recipe (post-relaunch)

Background poll script pattern used this session (local scratchpad,
`watch_nomerge_ladder.sh`): every 10 min ssh orc for
`sacct -j <ids> -X -P -n -o JobID,State` + `sinfo -p mgh -N -h -O
NodeList,GresUsed`; exit (waking the agent) on any terminal state, on
first RUNNING of an arm needing banner verification, or on a free mgh-1-1
GH200 while the merge-off NT36 pends (hop recipe:
gpu_nomerge_ladder_status_20260829.md §4 — scancel + resubmit from the
gh200 silo, never `scontrol update`; do NOT hop the rlxf arms' eng jobs
unless silo/arch constraints are rechecked, and never hop anything needing
LineGauss to gh200).

## 6. Guardrails (unchanged, ask Ryan first)

Notebook + ledger entries for the whole arc: approved-pending-draft, still
undrafted. Commits: launcher hooks above + band-aid still uncommitted.
No FLOWVPM/FastMultipole source edits (incl. merge-count instrumentation).
026 splitting on hold. Don't delete 1.5R-era data. ssh orc needs
`bash -lc` (or `bash -ls` + heredoc); `println(stderr,…)` not `@info`;
exit 124 = valid partial → `_s2` chain; sacct authoritative; rsync
`--checksum`. Ryan also wants to brainstorm other NT-problem approaches —
phase_17 candidate list + the two live ladders are the prep.
