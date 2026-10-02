# 032 root-shed-omission reruns — provenance (2026-09-18)

Ryan directive 2026-09-18: re-run the 026 A1-class blow-up arms with the
particle-side root omission active (feature B, `OmitStations`) to test
whether deleting the root-shed particles rescues the ~274–328-step
fountain Γ-ignition. Charter: `BRAINSTORM/032_reset_prompt_20260918.md`;
item: `BRAINSTORM/032_omit_shed_locations.md`; 026 baseline autopsies:
`026_sigma_growth_particle_splitting/rerunslate_provenance_20260918.md`.

## Pins (annotated tag `campaign/p032-rootomit-20260918`)

| repo | tag | commit | note |
|---|---|---|---|
| FLOWPanel.jl (`fastmultipole`) | `campaign/p032-rootomit-20260918` | `af92740` | 026 rerun-slate pin `6e52628` + 032 commit (filter_shedding + OmitStations + driver/launcher knobs) |
| FLOWVPM.jl (`flowpanel`) | `campaign/p026-rerunslate-20260918` | `8d4a3b4` | UNCHANGED from 026 slate (§22 merge redesign) |
| FastMultipole | `campaign/p026-derisk-20260914` | `ac7230a6` | UNCHANGED from 026 slate |

Worktrees:

- FLOWPanel: `orc:~/campaigns/p032-rootomit-20260918/FLOWPanel.jl`,
  detached at the tag (`af92740`), `dirty=0` verified 2026-09-18.
- FLOWVPM + FastMultipole: reuse the existing detached campaign worktrees
  `orc:~/campaigns/p026-derisk-20260914/{FLOWVPM.jl,FastMultipole}`
  (read-only, pins unchanged; the running 026 arms use the same trees).

Env: `orc:~/campaigns/p032-rootomit-20260918/env` — copy of the 026 slate
env with the FLOWPanel dev-path repointed at the p032 worktree; FLOWVPM /
FastMultipole dev-paths untouched (still the p026 worktrees). Verified by
`Pkg.instantiate()` + `pathof(FLOWPanel)` → p032 worktree,
`isdefined(FLOWPanel, :OmitStations)` → `true` (fill: see Submissions).

Tag pushed to `orc:projects/FLOWPanel.jl` by name (branches NOT pushed —
divergence ruling still parked with Ryan).

## Feature verification (local, pre-submission)

- Unit suites at `af92740`: liftingbody 56/56, wake 755/755 PASS.
- Feature-B smoke C (knob 0.15, 40_40, NT=36, 5 steps; predecessor's
  smokes C/D had both been killed with the session before post-march —
  rerun 2026-09-18): banner `Particle root-shed omission ACTIVE:
  sheddingN omits 2/36 stations (|r|/R < 0.15)` on both blades;
  post-march totals `Σ|Γ|·Δl = 1.593e-3 m³/s over 12 filaments`;
  `*_case_metadata.toml` records `particle_omit_root_r_over_R = 0.15`,
  `particle_omit_stations_{1,2} = 2`, circulation/filament totals.
- Knob-off smoke D (same recipe, knob unset): zero banners,
  `particle_omit_root_r_over_R = 0.0`, no station/total keys. Knob-off
  bit-identity is additionally unit-tested.
- Smoke recipe gotcha: `NREVS` alone does NOT shorten the run — the
  freestream schedule floors it at 13 revs (`required_revs =
  max(nrevs, schedule_revs)`); short smokes need
  `FREESTREAM_RAMP_REVS=0.15 FREESTREAM_HOLD_REVS=0
  FREESTREAM_WITHDRAW_REVS=0 SETTLE_REVS=0` as well.

## Arms

Base config = 026 A1 exactly (case `scr_p026s9_explg_fs`, FLOOR+SPLIT env
verbatim, linegauss, 12 h wall, m13h H200 — backend-matched to A1;
launcher `~/projects/FLOWPanel.jl/examples/run_p018_screen_gpu052.slurm.sh
h200 <case>`, cwd = p032 FLOWPanel worktree, `P018_REPO_OVERRIDE` /
`P018_PROJECT_OVERRIDE` at the p032 worktree/env). Only delta per arm:

| arm | RUN_NAME_OVERRIDE | delta vs A1 | question it answers |
|---|---|---|---|
| B15 | `scr_p032om15_explg_fs` | `PARTICLE_OMIT_ROOT_R_OVER_R=0.15` | does deleting the innermost root-shed particles (incl. chain-closing filament) rescue A1's step-328 ignition? |
| B20 | `scr_p032om20_explg_fs` | `PARTICLE_OMIT_ROOT_R_OVER_R=0.20` | dose response: wider root clip |

FLOOR = `SIGMA_FLOOR_FRAC=0.1`; SPLIT = `WAKE_SPLIT_VISCOUS=true
WAKE_SPLIT_FRAC_VISCOUS=0.587 WAKE_SPLIT_FRAC_COMPRESS=0.73
WAKE_SPLIT_FRAC_ELONGATE=0.3` (verbatim from the 026 slate).

Optional feature-A discriminator (`SHEDDING_R_OVER_R=0.2`, "move" vs
"delete" the root vortex) NOT submitted — queue-budget call parked with
Ryan. Extension to the A2 (floor-only) config: TRIGGERED by B15's clean
completion — B-floor arm in the Submissions table (knob 0.15, eng,
backend-matched to A2's eng H200 slot).

## Acceptance / read-out

- Success signal: no euler_exp guard trip in the ~274–328 window (A1 died
  328, A2 274, wave-2 twins 274/280); fountain-region Γ concentration
  absent or displaced in the VTK (A1 steps 278–327 VTPs local at
  `~/scr_p026s9r2_explg_fs_last50steps/` for before/after).
- Banner verification mandatory: `Particle root-shed omission ACTIVE`
  with sane masked counts (stock TE starts r/R 0.111 → 0.15 masks 2/36
  stations per blade at 40_40; production mesh counts recorded below),
  `CONVERSION=legacy`, plus the wave-2 checklist (linegauss, FLOOR+SPLIT
  knobs, expint on, SIGMA_CEIL=Inf, merge log + sigma telemetry).
- Report run-end `Particle root-shed omission totals` (deleted Σ|Γ|·Δl)
  alongside ΔCT vs A1 as the modeling cost.
- Judge by outputs, not sacct.

## Submissions (2026-09-18, cwd = p032 FLOWPanel worktree, wall 12 h)

Env verification pre-submission: `Pkg.instantiate()` clean;
`pathof(FLOWPanel)` = p032 worktree; `isdefined(FLOWPanel, :OmitStations)`
= `true` (login node; CUDA-init warning expected off-GPU).

`sbatch --time=12:00:00 --export=ALL,P018_REPO_OVERRIDE=…/p032-rootomit-20260918/FLOWPanel.jl,P018_PROJECT_OVERRIDE=…/p032-rootomit-20260918/env,FLOWPANEL_FILAMENT_REG=linegauss,SIGMA_FLOOR_FRAC=0.1,WAKE_SPLIT_VISCOUS=true,WAKE_SPLIT_FRAC_VISCOUS=0.587,WAKE_SPLIT_FRAC_COMPRESS=0.73,WAKE_SPLIT_FRAC_ELONGATE=0.3,RUN_NAME_OVERRIDE=<run>,PARTICLE_OMIT_ROOT_R_OVER_R=<knob> ~/projects/FLOWPanel.jl/examples/run_p018_screen_gpu052.slurm.sh h200 scr_p026s9_explg_fs`

| job | arm | case | RUN_NAME_OVERRIDE | partition |
|---|---|---|---|---|
| 13771353 | B15 | scr_p026s9_explg_fs | scr_p032om15_explg_fs | m13h (default) |
| 13771354 | B20 | scr_p026s9_explg_fs | scr_p032om20_explg_fs | m13h (default) |
| 13772015 | B-floor (A2 ext.) | scr_p026s9_explg_floor | scr_p032om15_explg_floor | eng (`--qos=eng --gres=gpu:h200:1`, backend-matched to A2; FLOOR only, no split knobs; knob 0.15; submitted 2026-09-18 ~23:50 after B15's clean completion, per the "extend if B rescues" directive; queues behind A4/A5 on eng) |

Run dirs land worktree-local (`data/scr_p032om*_explg_fs/`) — the launcher
`rm -rf`s any pre-existing `data/$RUN_NAME` (symlinks included), so the
026 precedent applies: MOVE to `orc:~/projects/FLOWPanel.jl/data/` +
symlink back at harvest. Slurm logs:
`logs/slurm/slurm-fp-052-scr-gpu-<jobid>.{out,err}` in the p032 worktree
(remember: DomainError tracebacks land in the `.err`).

## Banner verification

| job | arm | verdict |
|---|---|---|
| 13771353 | B15 | **PASS 10/10** (2026-09-18, m13h-1-1): repo/project = p032 worktree/env; 12 revs/468 steps; linegauss (pinned); floor 0.1 guard=on; split trio 0.587/0.73/0.3 (Resolution splitting ACTIVE, clamp=[2.374e-4, off]); WAKE_EXPINT=true; SIGMA_CEIL=Inf; CONVERSION=legacy; merge log + sigma telemetry at `data/scr_p032om15_explg_fs/`; **omission ACTIVE: 3/41 stations per blade masked (\|r\|/R < 0.15)**. Step 70/467 at check, ~4–5 s/step, no errors (FMM constant-redefinition warnings in .err are the known benign load-order noise). |
| 13771354 | B20 | **PASS 10/10** (2026-09-18, m13h-2-2): identical checklist; **omission ACTIVE: 5/41 stations per blade masked (\|r\|/R < 0.2)**; merge log at `data/scr_p032om20_explg_fs/`. Step 70/467 at check, no errors. |
| 13772015 | B-floor | **PASS** (2026-09-19, eng): omission ACTIVE 3/41 stations per blade (\|r\|/R < 0.15); floor 0.1 guard=on; **split OFF as intended** (zero "Resolution splitting ACTIVE" lines — floor-only arm matching A2); SIGMA_CEIL=Inf, expint on, CONVERSION=legacy. Step 35/467 at check, ~4 s/step. A2 died @274 — verdict window in ~1 h. |

## Outcomes

| job | arm | outcome |
|---|---|---|
| 13771353 | B15 | **COMPLETED 467/467** (2026-09-18 23:43, ~71 min wall): NO guard trip — the ~274–328 ignition window passed at a smooth ~9 s/step. CT converged: CYCLE-MEAN 0.0716217 ± 2.73e-5 (±0.038%, final 2 revs), Phase 2e CONVERGED=true, GATE gpu_gemv=468 cpu_gemv=0 nan_lines=0. Omission totals: deleted Σ\|Γ\|·Δl = 0.4863 m³/s over 2796 filaments (3 stations/blade). np ended ~367k — never approached the 500k cap A1 rode from step ~318; late steps ~11.6 s. wake-health tail sane (mean σ 9.8e-4, no floor-clamp column growth). **A1-class blow-up RESCUED by deleting root-shed particles.** ΔCT vs baseline pending (A1 died pre-readout; nearest completed baselines = wave-2/slate arms at harvest). |
| 13771354 | B20 | **COMPLETED 467/467** (2026-09-18 23:47, ~75 min wall): NO guard trip. CYCLE-MEAN CT 0.0713228 ± 1.02e-4 (±0.144%), CONVERGED=true, GATE clean. Omission totals: deleted Σ\|Γ\|·Δl = 0.7219 m³/s over 4660 filaments (5 stations/blade). np ended ~371k. Dose response: knob 0.15→0.20 deletes 48% more circulation for ΔCT of only −0.42% (0.07162→0.07132) — the rescue is not knife-edge in the knob. |

### Fountain-region before/after (2026-09-19)

Quantitative (wake-health extreme tail over A1's death window, steps
315–328): A1 ran min σ pinned AT the floor (2.374e-4) with max Γ/σ²
826→9.1e4; **B15 holds max Γ/σ² flat at ~46–47 and min σ at 1.19e-3 —
5× above the floor. The σ-at-floor, Γ-concentrated close-pair population
that drove the ignition never forms once the root shed is deleted.**
np at step 327: B15 288k vs A1 499k (cap-pinned).

ParaView setup for Ryan: local side-by-side dirs, same step range 278–327 —
`~/scr_p026s9r2_explg_fs_last50steps/` (A1) vs
`~/scr_p032om15_explg_fs_steps278_327/` (B15: 50 particle VTPs + 50 body
VTUs + monitors + case metadata, 2.3 G; open the VTP series directly, the
.pvd references all 468 steps). Offender coordinates to inspect:
x=−0.37R, r=0.34R (A1 step-327 knot).

### Merge-burden / merge-law interaction analysis (2026-09-19, Ryan question)

Total merge events, same case family across policies:

| arm | merge law | omission | merges | fate |
|---|---|---|---|---|
| wave-2 explg_fs | OLD | — | 83,705 | died @295 (overflow) |
| r2 A1 explg_fs | NEW | — | 87,287 | died @328 |
| B15 explg_fs | NEW | 0.15 | **2,991** | COMPLETED |
| B20 explg_fs | NEW | 0.20 | **2,816** | COMPLETED |
| wave-2 explg_floor | OLD | — | 2,819 | died @274 |
| r2 A2 explg_floor | NEW | — | 2,754 | died @274 |
| B-floor explg_floor | NEW | 0.15 | 2,190 | COMPLETED |

Reads: (1) the fs-family merge frenzy (~84–87k events) is driven by the
root-shed/fountain dynamics, NOT the merge law — omission collapses it
~30× to ~3k. (2) With so little merging left, the merge-law choice is
likely nearly inert under omission in this family: the new law's
demonstrated advantages (holding np at the 500k cap, +33-step survival)
acted exactly in the regime omission eliminates (np now peaks ~370k).
(3) The floor arm never had a merge frenzy (2.2–2.8k events under every
policy) — its death was direct Γ-ignition, consistent with zero
merge-redesign effect on its death step (274 twice). (4) Per-event σ
growth cannot be compared from merge_events.csv (no merged-σ column —
known telemetry gap). Definitive old-vs-new-under-omission A/B = one arm
(OLD-law FLOWVPM pin + omission 0.15); prediction: near-null on CT and
survival. Ryan-gated offer.

CT plateau quality of all completed runs in the family, for the
convergence question: B-floor ±0.014%, B15 ±0.038%, cap018 (old law,
σ-ceil) ±0.07%, B20 ±0.144%, cap030 (NT144) ±0.34% — the omission arms
are the tightest plateaus in the campaign to date. Separately, omission's
~30% lower np directly relieves the count-driven FMM σ-adequacy gate
(the killer of every 018 NT72 rung and A1's warning source steps
310–328) — candidate hidden improvement: NT ladders may now be runnable
at exp settings; the parked NT144 rung is the natural probe.

### ΔCT vs available baselines (2026-09-19; no clean twin exists — A1 died pre-readout, which is the point)

| baseline | CT | B15 Δ | B20 Δ | caveat |
|---|---|---|---|---|
| wave-2 cap018 (s9 family, σ-ceil 0.018, no omission) | 0.07395 ± 0.07% | **−3.15%** | −3.55% | nearest same-NT config, but carries a binding-adjacent σ ceiling instead of omission |
| de-risk cap030 chained rung | 0.07484 ± 0.34% | −4.30% | −4.69% | NT144 family — different time resolution, weak comparator |

Bounding argument: the omission's own dose response is −0.42% CT for a
48% increase in deleted circulation (0.15→0.20), so the *marginal*
CT-sensitivity to the clip is small; the ~3% offset vs cap018 conflates
omission cost with the cap arm's σ-ceiling physics and cannot be
attributed cleanly. A clean attribution needs an un-omitted A1-config
survivor, which does not exist — option for Ryan: accept the bound, or
commission a baseline at a survivable config.

| 13772015 | B-floor | **COMPLETED 467/467** (2026-09-19 08:22, ~65 min): NO guard trip — sailed through A2's twice-fatal step 274. CYCLE-MEAN CT 0.071654 ± 1.02e-5 (±0.014%), CONVERGED=true, GATE clean, zero .err errors. Omission totals: deleted Σ\|Γ\|·Δl = 0.4865 m³/s over 2796 filaments (matches B15's 0.4863 — shedding consistent across arms). Final np 318k, min σ 1.10e-3 (≈4.6× floor), max Γ/σ² 143 — extreme tail never approached ignition. Run dir moved to consolidated root + symlink (13 G). |

### FINAL 032 VERDICT (2026-09-19, all three arms in)

**Root-shed omission rescues every tested exp-family s9 class.** A1-config
(fs) at knobs 0.15 AND 0.20, and A2-config (floor-only) at 0.15, all
complete 467/467 converged with CT tightly clustered 0.0713–0.0717,
while every un-omitted variant of these configs has died at steps
274–328 across two campaign waves. Mechanism confirmed three ways:
(1) survival through the kill window; (2) the σ-at-floor Γ-concentrated
close-pair tail never forms (min σ stays 4.6–5× above the floor, max
Γ/σ² ~46–143 vs A1's 9.1e4); (3) np stays 290–370k where un-omitted arms
surged to the 500k cap. Ryan's hypothesis stands: the blade-root shed
circulation — of dubious physicality (root cutout / hub interference) —
was feeding the fountain-region ignition. Cost: explicit, accounted
(~0.49–0.72 m³/s deleted Σ\|Γ\|·Δl; CT clip-sensitivity −0.42% per +48%
deletion). Remaining discriminators are OPTIONAL (feature-A "move vs
delete"; ctrl-channel + omission — both Ryan-gated offers).

**Early verdict (2026-09-18): both omission arms rescue the A1-class
blow-up.** The exp-family s9 fountain Γ-ignition (A1 @328, A2 @274,
wave-2 twins @274/280) does not occur when the innermost root-shed
particles (incl. the chain-closing root filament) are deleted at the
panel→particle handoff; both arms complete 467/467 converged with CT
insensitive to the clip width (0.15 vs 0.20: −0.42%). Discriminators
still open: B-floor (13772015, A2-class floor-only) queued on eng;
feature-A "move vs delete" arm unsubmitted (Ryan queue-budget call).
