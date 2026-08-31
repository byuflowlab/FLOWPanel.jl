# Handoff — σ/VPM illustration campaign (2026-08-27, session 2 entry point)

Read `sigma_vpm_illustrations_plan.md` (same dir) first: it holds the approved plan —
mechanism principles, the PRIORITY-1 three-act story, figure specs F1–F7, animations
A1–A4, SFS caveats. This file is only the live state.

## Done (session 1)

- **Journal**: `# 20260827` entry in `~/Dropbox/research/notebooks/journals/20260821.md`,
  approved and written; all built figures/GIFs already linked in per-act groups; Status
  paragraph current as of jobs-running. NEVER tick its checkboxes (Ryan does).
- **Figures built** in `~/Dropbox/research/notebooks/img/20260827_sigma_vpm_illustrations/`
  (.tex + data dirs, pdflatex-verified, `\pagecolor{white}`, previews 150 dpi opaque):
  F1 `fig_sigma_map`, F2a/b `fig_viscous_roles_a/b`, F5 `fig_drift_arrest` (from real
  `p018_ufront_s035_visc` monitor CSV, local at `fig_drift_arrest/wake_health.csv`;
  measured drift −0.080/rev revs 3–25 vs −0.065/rev on the 019 window — both quoted),
  F6 `fig_sfs_channels`. Reused F3/F4/F7/fig_leg3 re-verified in their BRAINSTORM dirs
  and copied as white 150 dpi PNGs into the img dir.
- **GIFs built**: `a1_sigma_map_cobweb.gif`, `actII_blowup.gif` (local corpse VTPs
  `data/p018_L2_visc_forensics/`, steps 1032–1041, σ_shed=2.3006e-3).
- **Generators** (repo): `scripts/illus_render_particle_series.py` (shared VTP→GIF
  renderer: fixed camera, color=log10(σ/σ_shed), CLIM (−1.15,0.35), white bg),
  `scripts/illus_fig_sigma_map.py`, `scripts/illus_fig_viscous_roles.py`.

## In flight (SUPERSEDED — all four jobs landed 08-28; see Session-3 final state)

- **Jobs 13505318–21** (`fp-il-*`, m12/normal, 64c/64G, ~20 s/step ⇒ ≈2 h/run;
  submitted ~20:45 MDT 08-27): `scr_p019_s025_rr` (Act I blow-up, inviscid 0.0249R,
  ignites ~7.7 revs), `scr_p019_s025v_rr` (Act III fast, viscous 0.0249R, ignites
  ~7.9 revs), `scr_p019_s030v_rr` (Act I fixed, viscous 0.0299R survivor),
  `scr_p020r_geom_s020v_rr` (Act II fixed, expint, dies ~step 243 originally).
  Check `sacct -j 13505318,13505319,13505320,13505321 -X` — the session-1 monitor died
  with the session. First submission (13505310–13) FAILED on a cluster driver/src
  mismatch: driver passed `sigma_guard` kwarg the cluster `src/FLOWPanel_simulate.jl`
  doesn't accept. FIXED by patching the CLUSTER driver only (splat, passes kwarg only
  when guard non-empty); backups: cluster
  `examples/rotor_hover_pressure_comparison.jl.bak_20260827` + scratchpad copy. Do NOT
  sync cluster `src/` (live campaign chains restart from that checkout).
- Launcher `_rr` case arms added in `examples/run_p018_screen_hpc.slurm.sh` (local,
  uncommitted) and synced to cluster (md5-verified).

## Session-3 (2026-08-28) — superseded by Session 4

Session 3 iterated the actI_blowup GIF through the `gamma` colour scheme (opaque
sphere glyphs, colour+opacity from log10(|Gamma|/median)). That scheme is NO
LONGER USED by any GIF; `--color-by gamma` still exists in the renderer but is
dead code for this campaign. Its findings that still hold: plain world-size
gaussian splats read as fog (core overlap sigma/h~4) and were rejected; GIF
assembly must go through `magick -delay`, never imageio.

## Session-4 final state (2026-08-29, entry point for next agent)

**All seven act GIFs are now on ONE scheme (`--color-by abssigma`) and were
re-rendered together.** The renderer is `scripts/illus_render_particle_series.py`
(UNCOMMITTED, as is this whole plans/ dir). The exact commands that produced the
current files are in `plans/sigma_vpm_illustrations_20260827/render_all_acts.sh`
— run that script to reproduce all seven; it is the source of truth, read it
before changing anything.

### The `abssigma` scheme (Ryan-approved, arrived at over ~10 iterations)

Two glyph sets per particle, both coloured by |sigma|:

- **Spheres**, diameter `2|sigma|` (world units, UNCAPPED — `--diam-cap 0`),
  opacity **0.03**. These carry the wake's shape.
- **Arrows**, oriented along `gamma`, length = **|Gamma| at true 1:1 world
  scale** (`--arrow-scale 1.0 --arrow-cap 0`; no normalization, no clipping).
  Opacity 100% up to |Gamma| = R, fading linearly to 1% at |Gamma| >= 4R,
  R = 0.11938 m (DJI 9443). Ryan explicitly likes that arrows are INVISIBLE
  early (|Gamma| sits decades below R) and only appear as ignition runs away.
- **Colour**: `--sigma-log --sigma-clim 0.1:10 --sigma-clim-units shed`, i.e.
  log10 of |sigma|/sigma_shed over two decades centred on the shed size,
  viridis. Ryan called this "perfect" after trying 0.01:100, 1e-4:1 and
  1e-3:10 (all rejected: too flat, or one side of the story clamped away).
- `--yaw 0`, `--font-size 18` (bar text 21 via `--bar-font-size`, default
  font_size+3), `--nt 36 --spinup-steps 36`.
- HUD is title / step+rev / u_max only. The `min sig/sig_shed` line was
  REMOVED at Ryan's request, for all colour modes.

### Renderer gotchas discovered this session (all cost real time — read these)

1. **VTK silently drops Greek glyphs.** A sigma in any text renders as
   *nothing*, no error. Fixed by loading a Unicode TTF via
   `SetFontFamily(4)` + `SetFontFile(...)`; flag `--font-file`, default
   `/System/Library/Fonts/Supplemental/Arial Unicode.ttf`. On a box without
   that font the sigma vanishes again.
2. **Scalar-bar tick labels cannot be padded off the bar.** VTK centres
   annotation text ON the bar edge; `SetTextPad` and trailing-space hacks both
   do nothing. The decade labels are therefore drawn as our OWN right-justified
   text actors at `--bar-label-gap` left of the ramp, with VTK's annotations
   disabled.
3. **`SetVerticalTitleSeparation` (`--bar-title-gap`, 22 px) shrinks the ramp
   from the top by exactly that many pixels.** Our hand-drawn labels must
   follow or they drift off their decades — hence the measured inset in the
   code (5 px bottom, `20 + gap` top). If you change the bar geometry, re-measure.
4. **A single runaway glyph blanks the whole frame.** At a death step one
   particle can reach |Gamma| = 7.2e6, so a true-scale arrow is 7.2e6 m long;
   VTK auto-fits the depth range around it and clips the 0.35 m wake out of the
   frustum. Symptom: pure white frame with a colorbar (hit on
   actII_fixed_geom step 211). Fixed by pinning
   `pl.camera.clipping_range = (0.05*dist, 3.0*dist)`.
5. Arrow length spans ~3 decades of |Gamma|; median-anchored auto-scaling put
   single arrows at 7x the frame span. The percentile anchor
   (`--arrow-ref-pct/-len/-cap`) still exists but is unused now that everything
   runs at true 1:1.

### Current GIF inventory (all 7 rebuilt 2026-08-29, same scheme)

In `~/Dropbox/research/notebooks/img/20260827_sigma_vpm_illustrations/`:

| gif | series | sigma_shed | steps | frames |
|---|---|---|---|---|
| actI_blowup | vtp_s025 | 2.967e-3 | 266:305 | 40 |
| actI_fixed | vtp_s030v | 3.560e-3 | 264:323 | 60 |
| actII_blowup | data/p018_L2_visc_forensics (in repo) | 2.3006e-3 | 1032:1041 | 10 |
| actII_fixed_geom | vtp_geom_s020v | 2.374e-3 | 180:211 | 32 |
| actIII_fast_ignition | vtp_s025v | 2.967e-3 | 285:323 | 39 |
| actIII_drift_death | vtp_ufront_s035 | 4.236e-3 | 987:996 | 10 |
| actIII_fixed | vtp_s038v | 4.531e-3 | 314:323 | 10 |

Windows are now contiguous stride-1: the sparse lead-in frames that used to be
spliced in (actII_fixed_geom 150/170, actIII_fast_ignition 260/280) were dropped
because they jumped ~20 steps in a single frame. actI_fixed's prescribed
rev-6.75 milestone overlay was NOT carried over — that rev came from the old
record and ignition shifted on the rerun stack; ask Ryan for the rev he wants
before re-adding `--milestone`.

Two frames worth knowing: actII_fixed_geom step 211 shows a visually pristine
wake at u_max 4.3e9 m/s (the Act II point made visible), and actIII_drift_death
carries almost no arrows because its death is a slow sigma drift, not a Gamma
runaway — the spheres do the work there.

### OPEN — needs Ryan, do not resolve alone

- **Clipping-range change to actI_blowup.** Gotcha 4's fix is applied to all
  seven, including actI_blowup, which Ryan had ALREADY approved before the fix.
  It only alters the last few frames (giant arrows no longer blanket the view;
  normal frames are bit-identical). Flagged to him twice, never answered.
  The other six need the fix regardless; reverting actI alone is a one-flag job.
- **Caption rewrite for all seven.** THIS IS THE NEXT TASK. Every caption in the
  journal describes a superseded scheme. The actI_blowup caption at
  `journals/20260821.md:527` was rewritten mid-session but describes the
  LINEAR-scale version and is now stale too. Several captions quote a
  `min sigma/sigma_shed` number that is no longer on the HUD, and
  actII_blowup / actII_fixed_geom / actIII_fast_ignition quote step ranges that
  the window changes invalidated. Ryan said "update the caption" for actI only;
  he has NOT yet green-lit the seven-caption pass — offer it, then do all seven
  in ONE edit. Journal rules: append-only for past days, never tick checkboxes.

### Data (all VTPs local; NO cluster access needed to re-render)

6.7 G in the SESSION-3 scratchpad, which is in /private/tmp and may be swept:
`/private/tmp/claude-502/-Users-ryan-Dropbox-research-projects-FLOWPanel-jl/2c6862dd-d8a5-430e-a379-84e453dda661/scratchpad/`
- `vtp_s025/` 266-323 + 240/260, `vtp_s030v/` 264-323, `vtp_geom_s020v/`
  180-211 + 150/170, `vtp_s025v/` 285-323 + 260/280, `vtp_ufront_s035/`
  987-996, `vtp_s038v/` 314-323.
- Wake-health HUD CSVs for the `_rr` runs + s038v sit in that scratchpad root.
- In-repo (safe): `data/p018_L2_visc_forensics/` (actII blow-up VTPs +
  death_trajectory.csv) and `data/p018_ufront_s035_visc/monitors/`.
**If that scratchpad is gone, the VTPs must be re-pulled from the `_rr` run dirs
on /home before any re-render.** `render_all_acts.sh` hardcodes the path in `$S`.

### Background jobs / cluster

Nothing in flight. All four `_rr` render jobs (13505318-21) are DONE and
harvested; 13505321's "FAILED" status is its in-sim death at step 211, not a
launcher problem. `_rr` run dirs still on /home; keep/archive is Ryan's call.

## Pending Ryan decisions carried over (do not act without him)

- **RECENT archive approvals**: /home 543 G vs 400 G cap; 29 quiet runs ~238 GiB
  need his per-checkout sign-off
  (`--include-recent --only RUN,... --root <checkout>`). Ledger line already
  appended to the 018 ledger.
- **Protect list additions** (his file): `p018_ufront_s035_visc`, `scr_p019_s038v`.
- Whether to re-sync the hybrid cluster checkout after the _s2/_3r chains land.

## Gotchas (campaign-level, still current)

- Original run dirs (`scr_p019_s025` etc.) hold the 019/020 monitor CSVs the
  published figures cite — the launcher `rm -rf`s an existing run dir on
  resubmit; that is WHY the reruns use `_rr` tags. Never resubmit under the
  original tags.
- Rerun results are render-only (post-08-24 dep stack; julia 1.12.6 env):
  ignition steps may shift vs the record. Do not score them against 019/020
  verdicts.
- `ssh orc` needs the live ControlMaster socket; `sbatch`/`squeue` need `bash -lc`.
- Journal policy: append-only for past days; ask Ryan before new headers or
  content beyond this campaign's established entry.

## Session-5 (2026-08-29) — regime-based restructure (entry point for next agent)

Ryan approved restructuring the acts to the regime flowchart
(../img/20260805_sigma_blowup/regime_flowchart.png). Approved plan:
`~/.claude/plans/read-plans-sigma-vpm-illustrations-20260-wiggly-hartmanis.md`.
Key data finding driving it (wake-health CSVs, now durably in
`plans/sigma_vpm_illustrations_20260827/wake_health_csvs/`): the viscous twin
s025v at the SAME sigma_shed ignites too (rev 8.11 vs 7.83 inviscid, HUD revs
7.11/6.83 — HUD rev = (step-36)/36, one less than step/36; BEWARE the two
conventions); s030v is NOT a survivor on the rerun stack (ignites rev 8.47) and
is retired; negative sigma in s025 appears 1 step AFTER ignition (post-mortem
artifact, not mechanism).

New act set (render_all_acts.sh rewritten, source of truth):
- actI_blowup — unchanged frames (s025 266:305). Title still says "sigma through
  zero"; proposed retitle pending Ryan (frames otherwise identical).
- actI_fixed — NEW: s025v 266:288, stops before its own ignition (step 292),
  --milestone "6.833:inviscid twin died here". BUILT.
- actI_sidebyside.gif — composite variant (magick +append, frames 266:288),
  BUILT; Ryan to pick milestone version vs composite.
- actII_runaway — s025v 285:323 (ex-actIII_fast_ignition retitled). BUILT.
- actIII_drift_death — unchanged (s035 987:996).
- actIII_fixed_matched / actIII_fixed_steady — pending job 13508659.
- RETIRED: old actI_fixed (s030v), actII_blowup, actII_fixed_geom (F7 coda only).

Long sigma* run: job 13508659 `scr_p019_s038v_gpu40` RUNNING on mgh-1-1 (gh200)
since 09:52 MDT 08-29 — s038v knobs + NREVS=40, GPU stack, banner verified.
Cluster-side edits (gh200 silo ~/FLOWPanel-018-gpu-gh200 ONLY, backups .bak_20260829):
screen dispatcher accepts $P018_REPO cwd + honors P018_JULIA/P018_PROJECT + new
case arm scr_p019_s038v_gpu40; NEW wrapper examples/run_p018_screen_gpu.slurm.sh.
Main checkout untouched. Run dir: ~/FLOWPanel-018-gpu-gh200/data/scr_p019_s038v_gpu40
(~25-30G VTPs expected — mind the disk cap). If it dies before ~rev 35, that is a
REAL finding (stack shifted the s030 boundary); bump sigma with a small safety
factor (0.040-0.042R) and relaunch per Ryan.

019 production hunt result: the 40-rev sigma* run is `p018_L1_visc` (local+cluster
data/p018_L1_visc), CT 0.0716-0.0723 healthy, but NO wake-health monitor and NO
VTPs — stability point only. So 13508659 is the only source of a long min-sigma/M
trace.

New figure: fig_sigma_drift_overlay.tex (+ data dir) in the img dir — min
sigma/sigma_shed vs revs for s025/s025v/s035/s038v, criterion numbers annotated;
PROVISIONAL until the long-run curve is appended. Compiles.

020 item: Phase 3a (minimal nu_sgs prototype, staged only) inserted before
Phase 4 + log line — Ryan approved STAGING only.

OPEN: per-GIF message approval loop with Ryan (success criterion #1); actI_blowup
retitle; actIII_fixed pair after job lands; journal caption pass (ONE edit, ask
first); drift-overlay final; HUD-vs-record rev convention decision.

### Session-5 addendum (later 2026-08-29)

- actI_blowup.gif REVERTED to the approved 266:305 faded-arrow version (Ryan);
  the opaque 266:288 twin lives on as actI_blowup_opq/ and is ONLY the
  composite's left pane. actI_sidebyside.gif = milestone-free, opaque-arrow,
  matched-window (266:288) pair — THE Act-II "fix helps" presentation.
  render_all_acts.sh reflects all of this.
- Journal entry `# 20260827` RESTRUCTURED per Ryan's instruction: per-act
  subheaders (### Act I/II/III), integrator substory after Act III, all
  captions rewritten to the regime story with measured numbers; status
  addendum appended. Retired links (actI_fixed s030v, actII_blowup,
  actIII_fast_ignition) removed from the entry.
- GIF messages: Ryan approved actII_runaway (ex-fast-ignition) and
  actIII_drift_death formatting; actI_blowup approved in reverted form.
  Still pending: actIII_fixed pair from job 13508681 (40-rev sigma* GPU run,
  RUNNING on mgh-1-1 with SCR_GPU_RESERVE_GIB=16 env hook).

## Session-5 FINAL STATE (2026-08-29, context reset — ENTRY POINT for next agent)

Read the approved plan first:
`~/.claude/plans/read-plans-sigma-vpm-illustrations-20260-wiggly-hartmanis.md`
(regime-based three-act restructure; Ryan-approved). Everything below
supersedes earlier sections of this file where they conflict.

### DONE and Ryan-approved
- Act structure = regime flowchart. Journal entry `# 20260827` in
  `~/Dropbox/research/notebooks/journals/20260821.md` RESTRUCTURED (Ryan
  ordered the in-place edit): subheaders `### Act I/II/III`, integrator
  substory after Act III, SFS coda; all captions rewritten with measured
  numbers; retired GIF links removed; status addendum appended.
- GIFs (in `~/Dropbox/research/notebooks/img/20260827_sigma_vpm_illustrations/`),
  all on the abssigma scheme; rebuild = `render_all_acts.sh` (source of truth):
  - actI_blowup.gif — REVERTED to approved 266:305 faded-arrow version.
    Opaque 266:288 twin kept as actI_blowup_opq/ = composite left pane ONLY.
  - actI_sidebyside.gif — Act-II "fix helps" presentation: inviscid|viscous,
    matched steps 266:288, arrows FULLY opaque (--arrow-opacity-min 1.0),
    no milestone text. APPROVED.
  - actI_fixed.gif — standalone right pane (building block, not in journal).
  - actII_runaway.gif — s025v 285:323 (ex actIII_fast_ignition). APPROVED.
  - actIII_drift_death.gif — unchanged. APPROVED.
  - actIII_fixed.gif — 9-rev remnant, PLACEHOLDER in journal until cpu40 lands.
- fig_sigma_drift_overlay.tex + data dir (img dir) — compiles; PROVISIONAL:
  add cpu40 (and optionally gpu40 partial) curves when available.
- 020 item: Phase 3a nu_sgs prototype STAGED (section before Phase 4 + log
  line). Staging only.
- 052 note filed: `../FastMultipole/MATRIX_OPERATOR_REFACTOR/
  052-gpu-memory-growth-longrun-2026-08-29.md` + pointer in
  052d-handoff-prompt-2026-08-28j.md (GPU leak evidence, regression recipe).

### IN FLIGHT — job 13508968 `scr_p019_s038v_cpu40` (THE gating item)
40-rev sigma* rerun, CPU, main checkout /home/rander39/projects/FLOWPanel.jl,
launched ~11:10 MDT 08-29 on m12-2-18, ~19-20 s/step x 1475 steps ~ 9 h
(+12 h wall). Case arm `scr_p019_s038v_cpu40` added to LOCAL
examples/run_p018_screen_hpc.slurm.sh and synced to cluster main checkout
(md5 100ef206466c0e0aa28b27fda7b7f3bf, verified). Run dir
data/scr_p019_s038v_cpu40 will be LARGE (~60-80 G VTPs) — sweep after harvest
per archive-first policy, keep the two render windows.
GPU attempt history: 13508659 (OOM at S-upload, reserve 32 GiB), 13508678
(same — fm052 env bundle clobbers env; fixed via SCR_GPU_RESERVE_GIB hook in
the NEW silo wrapper examples/run_p018_screen_gpu.slurm.sh), 13508681 (ran
819/1475 steps at 6.3 s/step then died on the ~35 MB/step device leak — see
052 note). Salvaged gpu40 monitor CSV (820 rows, 22.7 revs) IS ALREADY in
plans/sigma_vpm_illustrations_20260827/wake_health_csvs/ — key numbers:
min sig ratio 0.491@rev3 -> 0.145@rev14 -> 0.103-0.104 plateau revs 19-22.7
(drift ARREST emerging, local slope ~ -0.003/rev vs s035's -0.08/rev);
peak max_dtZ after rev 10 = 0.365 > eps=0.2, transient, no ignition — supports
Ryan's "certify at the peak M of the converged fluctuation band" framing.
GPU run dir ~/FLOWPanel-018-gpu-gh200/data/scr_p019_s038v_gpu40 (18 G) is
redundant once cpu40 lands — flag for cleanup/archive, do NOT delete alone.

### NEXT STEPS when 13508968 lands (a monitor was armed in-session; after
reset, just check `sacct -j 13508968` / the run dir)
1. Judge by outputs (sacct FAILED is unreliable): harvest
   monitors/scr_p019_s038v_cpu40_monitor04_wake_health_system1.csv ->
   drift rate + plateau + peak M in the converged band. If it IGNITES before
   ~rev 35: REAL finding (stack shifted the boundary; s030v already ignites on
   this stack) — report to Ryan, then bump sigma to 0.040-0.042R ("criterion +
   small safety factor") and relaunch per his standing instruction.
2. Pull VTPs for two windows: matched (steps 987-996, HUD revs 26.42-26.67 —
   same as actIII_drift_death) and steady loop (~revs 30-35, pick an integer
   number of blade passages so the GIF loops cleanly).
3. Render actIII_fixed_matched + actIII_fixed_steady (add stanzas to
   render_all_acts.sh; sigma_shed for s038 family = 4.531e-3, hud CSV = the
   cpu40 monitor). Captions must quote peak M vs eps and measured arrest rate.
4. Replace the journal PLACEHOLDER caption (actIII_fixed) with the pair;
   append cpu40 + gpu40 curves to fig_sigma_drift_overlay (rev axis note:
   figure currently uses step/36; HUD uses step/36 - 1 — journal preamble
   documents the HUD convention).
5. Offer Ryan a status addendum for the journal; NEVER tick checkboxes.

### Open Ryan decisions / loose ends
- actI_blowup in-frame title still says "sigma through zero" (Ryan kept the
  reverted frames; a retitle would re-render — only do it if he asks).
- Overlay figure rev-axis convention (step/36 vs HUD) — pick when finalizing.
- RECENT archive approvals + protect-list additions (carried over, untouched).
- _rr run dirs + gpu40 run dir keep/archive: Ryan's call.

### Gotchas that cost time this session
- Background Bash here is ZSH: unquoted $FLAGS strings do NOT word-split —
  wrap multi-flag renders in a bash script file (scratchpad/build_actI.sh
  pattern) or renders fail with "unrecognized arguments".
- ssh orc needs the live ControlMaster socket; sbatch/squeue/sacct need
  `bash -lc`; MOTD pollutes ssh output — grep defensively in monitors, and
  double-check "job left queue" via sacct before believing it.
- The fm052 GPU env bundle overrides sbatch --export env (that is why
  SCR_GPU_RESERVE_GIB hook exists in the screen GPU wrapper, silo only).
- Journal HUD rev convention: HUD rev = step/36 - 1 (spinup excluded).

## Session-6 FINAL STATE (2026-08-29 ~13:30 MDT, context reset — ENTRY POINT)

Short session. Everything in "Session-5 FINAL STATE" above still stands
(read it in full — plan pointer, DONE list, NEXT STEPS 1-5, gotchas).
Only delta:

- Job 13508968 `scr_p019_s038v_cpu40` (sacct name fp-il-s038v-cpu40) did
  NOT start at ~11:10 as Session-5 assumed — it was found RUNNING on
  2026-08-29 ~13:00 MDT at step 48/1475, ~25-32 s/step (CPU, slower than
  the 19-20 s/step estimate) → lands ~01:00-02:00 MDT 2026-08-30
  (~12 h from step 48). Health at step 48: nominal — sigma_ratio 0.923
  @step 40, sigma_growth ~1.010, forces steady, no ignition signatures.
- A persistent sacct poll monitor (10-min cadence) was armed in the OLD
  session; it dies with the context reset. NEXT AGENT: re-check
  `ssh orc "bash -lc 'sacct -j 13508968 ...'"` first thing and re-arm a
  completion watch if still running.
- No files were changed this session besides this HANDOFF appendix.
  No renders, no journal edits, no cluster changes.

When the job lands: execute Session-5 "NEXT STEPS when 13508968 lands"
items 1-5 verbatim (judge by outputs, not sacct state).

## Session-5 addendum (2026-08-29) — fig_velocity_blowup

New figure in the img dir, linked in the journal's Act II right after
`actI_sidebyside.gif`: `fig_velocity_blowup.tex` (+ data dir, + `_preview.png`),
generated by `scripts/illus_fig_velocity_blowup.py` (UNCOMMITTED, like the rest
of the campaign). Two panels vs `sigma` on a log x axis, single Gaussian-erf
particle of UNIT strength: (a) `max|u| = C_u/sigma^2`, (b)
`max|grad u| = C_g/sigma^3`, `C_u = 1.702951e-2`, `C_g = 2.993119e-2` (computed
in the script by maximizing over rho and theta; a vectorized and an
explicit-tensor path cross-check each other via asserts).

Ryan's styling ruling on this figure (2026-08-29): bare axes only — no grid, no
boundary box, no title, no legend, no in-figure caption, no marker/annotation
lines. Deliberately stripped from an earlier draft. Do not re-add decoration.

Y axis is log (curves span 6 and 9 decades; a linear y collapses them into a
wall at the left edge). The .tex header documents the switch.

**Correctness note (a trap I fell into and backed out of):** do NOT overlay the
death-step field maxima as a "check" on these curves. Measured off the step-1041
VTP directly, max|u| = 2.7853e5 m/s and max||grad u||_F = 5.2406e8 1/s, and BOTH
are attained at a particle with sigma = 2.0543e-3 m — not at the floor — whose
own |Gamma| has already run away to 11.25 (vs 1.2956e-2 at the healthy step
1032). Those are collective-field numbers at a runaway Gamma; the figure is a
single-particle fixed-strength result, so they are not comparable. The mechanism
doc's quoted "733 -> 37,088 m/s" is from a different step than 1041 — read
numbers off the VTPs, not off the prose.

## Status note (2026-08-29, evening) — gpu40 RELAUNCHED with the leak fix (from the 052 agent)

- **Job 13512297** `scr_p019_s038v_gpu40`, RUNNING on **mgh-1-2** (GH200) since
  ~18:5x MDT 08-29. Run dir `~/FLOWPanel-018-gpu-gh200/data/scr_p019_s038v_gpu40`
  (the OLD 18 G dir from 13508681 was auto-preserved by the dispatcher as
  `data/scr_p019_s038v_gpu40.prev` — nothing deleted). Expect ~3-6.3 s/step
  (np-dependent), 1475 steps, well inside the 24 h wall.
- **Leak-fix status: VALIDATED before launch.** Root cause of the 13508681 OOM
  (~35 MB/step): FastMultipole's cached device finalize/SFS output buffers
  reallocate on every particle-count change; the replaced buffers survive one
  step, get promoted, and the host GC never collects them (device bytes don't
  trigger its heuristics). Fixed with grow-only capacity buffers + prefix
  views, patched into the gh200 018 silo (FastMultipole + FLOWVPM,
  `.bak-052leak` backups beside each file). Validation probe job 13511798
  (`scr_p019_s038v_leakprobe`, 10 revs = 396 steps, same physics): pool_reserved
  froze at 15.27 GB from gemv 40→390, free pinned at 21.36 GB, pool_used
  trendless 11.7-14.5 GB — the old run grew steadily over the same window.
  (Probe sacct shows FAILED — that is the wrapper's NaN-line gate tripping on
  the short run's empty CT-plateau printout, `CT ... plateau mean=NaN`; the run
  itself finished 396/396 steps, all-GPU gemvs, all finite = true.)
- **Regularization: gaussian-era filaments** (the 018 silo stack predates the
  LineGauss pin and the device cross pass entirely) — Ryan's explicit choice
  for cpu40 comparability. Do NOT compare it against LineGauss expectations.
- cpu40 (13508968) still RUNNING at ~9 h elapsed on m12-2-18, on track.
- When 13512297 lands, harvest per the "NEXT STEPS when 13508968 lands" recipe
  (monitor CSV `monitors/scr_p019_s038v_gpu40_monitor04_wake_health_system1.csv`,
  VTP windows, renders); judge by outputs, not sacct state.
