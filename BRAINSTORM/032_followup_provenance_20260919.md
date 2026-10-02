# 032 follow-up arms — provenance (2026-09-19)

Ryan-approved 2026-09-19 (AskUserQuestion): launch arms 1+2 of the ranked
slate from `032_reset_prompt_20260919b.md`; feature-A move-vs-delete and
NT144 un-park NOT approved this round; ΔCT baseline declined (accept the
dose-response bound). Parent records: `032_rerun_provenance_20260918.md`
(B15/B20/B-floor, FINAL 032 VERDICT),
`026_sigma_growth_particle_splitting/rerunslate_provenance_20260918.md`
(A1–A5 slate autopsies).

## Arms

| arm | question | case | RUN_NAME_OVERRIDE | mech env | partition |
|---|---|---|---|---|---|
| C-ctrl | does root shed also drive the pump-independent ctrl blow-up (A3 onset 253–288)? | scr_p026s9_ctrllg_fs | scr_p032om15_ctrllg_fs | FLOOR+SPLIT + omission 0.15 (WAKE_EXPINT=false from case def) | eng (backend-match A3) |
| C-oldlaw | was the §22 merge redesign necessary once root shed is gone? (prediction on record: near-null vs B15) | scr_p026s9_explg_fs | scr_p032om15_explg_fs_oldlaw | FLOOR+SPLIT + omission 0.15 | m13h (backend-match A1/B15) |

FLOOR = `SIGMA_FLOOR_FRAC=0.1`; SPLIT = `WAKE_SPLIT_VISCOUS=true
WAKE_SPLIT_FRAC_VISCOUS=0.587 WAKE_SPLIT_FRAC_COMPRESS=0.73
WAKE_SPLIT_FRAC_ELONGATE=0.3`; omission = `PARTICLE_OMIT_ROOT_R_OVER_R=0.15`.
Wall 12 h, launcher `~/projects/FLOWPanel.jl/examples/run_p018_screen_gpu052.slurm.sh h200 <case>`.

## Pins

C-ctrl reuses the p032 pins verbatim (no new tag needed):
FLOWPanel `af92740` (`campaign/p032-rootomit-20260918`), FLOWVPM `8d4a3b4`
(`campaign/p026-rerunslate-20260918`), FastMultipole `ac7230a6`
(`campaign/p026-derisk-20260914`); worktrees/env
`orc:~/campaigns/p032-rootomit-20260918/{FLOWPanel.jl,env}` +
`orc:~/campaigns/p026-derisk-20260914/{FLOWVPM.jl,FastMultipole}`.

C-oldlaw pin delta (annotated tag `campaign/p032-oldlaw-20260919`):

| repo | tag | commit | note |
|---|---|---|---|
| FLOWVPM.jl | `campaign/p032-oldlaw-20260919` | `2b253db` | OLD merge law (pre-`8d4a3b4` §22 redesign) = wave-2 pin `campaign/p026-wave2-20260915` |
| FLOWPanel.jl | `campaign/p032-rootomit-20260918` | `af92740` | unchanged |
| FastMultipole | `campaign/p026-derisk-20260914` | `ac7230a6` | unchanged |

Tag created locally at `2b253db`, pushed to `orc:projects/FLOWVPM.jl` by
name. Worktree `orc:~/campaigns/p032-oldlaw-20260919/FLOWVPM.jl` detached
at the tag, `dirty=0` verified 2026-09-19. Env
`orc:~/campaigns/p032-oldlaw-20260919/env` = copy of the p032 env with the
FLOWVPM dev-path repointed at the old-law worktree (FLOWPanel /
FastMultipole dev-paths untouched). Verified pre-submission:
`Pkg.instantiate()` clean; `pathof(FLOWPanel)` → p032 worktree;
`pathof(FLOWVPM)` → old-law worktree; `isdefined(FLOWPanel, :OmitStations)`
= `true` (login node; CUDA-init warning expected off-GPU).

## Submissions (2026-09-19, cwd = p032 FLOWPanel worktree, wall 12 h)

C-ctrl (eng):
`sbatch --time=12:00:00 -p eng --qos=eng --gres=gpu:h200:1 --export=ALL,P018_REPO_OVERRIDE=…/p032-rootomit-20260918/FLOWPanel.jl,P018_PROJECT_OVERRIDE=…/p032-rootomit-20260918/env,FLOWPANEL_FILAMENT_REG=linegauss,SIGMA_FLOOR_FRAC=0.1,WAKE_SPLIT_VISCOUS=true,WAKE_SPLIT_FRAC_VISCOUS=0.587,WAKE_SPLIT_FRAC_COMPRESS=0.73,WAKE_SPLIT_FRAC_ELONGATE=0.3,RUN_NAME_OVERRIDE=scr_p032om15_ctrllg_fs,PARTICLE_OMIT_ROOT_R_OVER_R=0.15 ~/projects/FLOWPanel.jl/examples/run_p018_screen_gpu052.slurm.sh h200 scr_p026s9_ctrllg_fs`

C-oldlaw (m13h default): same export string but
`P018_PROJECT_OVERRIDE=…/p032-oldlaw-20260919/env`,
`RUN_NAME_OVERRIDE=scr_p032om15_explg_fs_oldlaw`, case
`scr_p026s9_explg_fs`, no partition flags.

| job | arm | submitted |
|---|---|---|
| 13773490 | C-ctrl | 2026-09-19 (eng) |
| 13773491 | C-oldlaw | 2026-09-19 (m13h) |

Run dirs land worktree-local (launcher `rm -rf`s pre-made symlinks) —
MOVE to `orc:~/projects/FLOWPanel.jl/data/` + symlink back at harvest.
Slurm logs `logs/slurm/slurm-fp-052-scr-gpu-<jobid>.{out,err}` in the p032
worktree; tracebacks land in the `.err`.

## Acceptance / read-out

- C-ctrl: completion + converged CT → root shed also drives the ctrl
  channel (omission is a general rescue); in-window blow-up (~253–290)
  → ctrl channel is root-shed-independent (separate investigation stands).
- C-oldlaw: CT within noise of B15 (0.07162 ± 0.038%) and completion →
  merge redesign unnecessary once root shed is omitted (prediction
  holds); death or CT shift → the redesign contributes independently.
- Banner verification owed on start: omission ACTIVE 3/41 stations per
  blade, CONVERSION=legacy, wave-2 checklist (linegauss, FLOOR+SPLIT
  knobs, SIGMA_CEIL=Inf, merge log + sigma telemetry); C-ctrl must show
  WAKE_EXPINT=false; C-oldlaw must show the OLD merge law loaded (FLOWVPM
  path = oldlaw worktree in the banner/metadata).
- Judge by outputs, not sacct.

## Banner verification

| job | arm | verdict |
|---|---|---|
| 13773490 | C-ctrl | **PASS 7/7** (2026-09-19, eng): repo/project = p032 worktree/env; omission ACTIVE 3/41 stations per blade (\|r\|/R < 0.15) both blades; **WAKE_EXPINT=false (ctrl — correct)**; linegauss, split trio 0.587/0.73/0.3, floor 0.1 guard=on, SIGMA_CEIL=Inf, CONVERSION=legacy; 12 revs/468 steps; merge log + sigma telemetry at `data/scr_p032om15_ctrllg_fs/`; .err clean (benign FMM warnings only). Step 37/467 at check. |
| 13773491 | C-oldlaw | **PASS 8/8** (2026-09-19, m13h): repo = p032 FLOWPanel worktree, **project = p032-oldlaw env; FLOWVPM dev-path = oldlaw worktree (`2b253db`) confirmed in banner `project:` line + Manifest**; omission ACTIVE 3/41 per blade; WAKE_EXPINT=true; linegauss, split trio, floor 0.1, SIGMA_CEIL=Inf, CONVERSION=legacy; merge log + telemetry at `data/scr_p032om15_explg_fs_oldlaw/`; .err clean. Step 26/467 at check. |

Verdict window (un-omitted deaths 253–328) expected ~40–70 min after
start; harvest owed on completion (MOVE + symlink, CT readout, merge
counts, wake-health tail).

## Outcomes (2026-09-19, harvested same day)

| job | arm | outcome |
|---|---|---|
| 13773490 | C-ctrl | **COMPLETED 467/467** (09:02→10:13, 1 h 10 m 39 s): NO guard trip — sailed through A3's 253–288 window. CYCLE-MEAN CT 0.0714872 ± 4.77e-5 (±0.0667%, final 2 revs), CONVERGED=true, GATE gpu_gemv=468 cpu_gemv=0 nan_lines=0, .err clean. Omission totals: deleted Σ\|Γ\|·Δl = 0.4855 m³/s over 2796 filaments (3/blade). Final np 361,916; merge events 2,742. Wake-health tail (steps 418–467): min σ 9.97e-4 (≈4.2× floor), max Γ/σ² 107.5 — no ignition tail. Run dir 16 G, moved to `~/projects/FLOWPanel.jl/data/` + symlink. |
| 13773491 | C-oldlaw | **COMPLETED 467/467** (09:02→10:16, 1 h 14 m 14 s): NO guard trip. CYCLE-MEAN CT 0.0716405 ± 9.84e-7 (±0.00137% — tightest plateau in the campaign), CONVERGED=true, GATE clean, .err clean. Omission totals: 0.4862 m³/s over 2796 filaments (matches B15's 0.4863 — shedding consistent). Final np 368,444; merge events 3,088 (vs B15's 2,991 under the NEW law — near-identical). Wake-health tail: min σ 9.60e-4 (≈4× floor), max Γ/σ² 54.0. Run dir 16 G, moved + symlinked. |

### Verdict (both acceptance criteria met on the "prediction holds" side)

**C-ctrl: root shed also drives the ctrl channel — omission is a general
rescue.** The A3 config (WAKE_EXPINT=false, FLOOR+SPLIT) that blew up
physically at 253–288 completes converged under omission 0.15 with the
same healthy-tail signature as B15/B-floor (min σ ~4× floor, max Γ/σ²
~108 vs ignition-class 10⁴–10⁵). The "separate ctrl-channel
investigation" contingency is CLOSED — no root-shed-independent ctrl
mechanism remains in evidence.

**C-oldlaw: the §22 merge redesign is unnecessary once root shed is
omitted — the on-record near-null prediction holds.** CT 0.0716405 falls
inside B15's ±0.038% noise band (Δ = +0.026%), the run completes, and the
OLD law's merge count (3,088) is statistically indistinguishable from the
NEW law's (2,991) — confirming the merge-burden read that the fs-family
merge frenzy was root-shed-driven, not law-driven, and the law choice is
inert in the omission regime. The merge redesign remains justified only
for un-omitted / cap-riding regimes.

Completed-family CT summary under omission 0.15: B15 0.07162, B-floor
0.07165, C-oldlaw 0.07164, C-ctrl 0.07149 (ctrl case def, expint off),
B20 (0.20 clip) 0.07132 — the exp-family cluster spans 0.03% across
merge law and floor/split policy; the ctrl offset (−0.19%) is the
WAKE_EXPINT=false model difference, not noise.
