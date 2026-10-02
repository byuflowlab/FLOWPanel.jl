# Omission re-open prep — 022 / 018 / 020 arms (2026-09-19)

Ryan-approved 2026-09-19 (AskUserQuestion, this session): PREPARE — do not
submit — omission-enabled reruns for (1) 022 hr10 IGE, (2) 018 NT72 `_3r`,
(3) 020 Phase-2R discriminator. NT144 (018/026) is CONDITIONAL on the 018
NT72 arm succeeding — do not prepare/submit before that verdict. 026
NT144 un-park and feature-A move-vs-delete remain declined. Sweep evidence
base: scout briefs summarized in the session of this date;
032 verdicts in `032_rerun_provenance_20260918.md` +
`032_followup_provenance_20260919.md` (C-ctrl general rescue; C-oldlaw
merge-law inert under omission).

Common delta for all three arms: add `PARTICLE_OMIT_ROOT_R_OVER_R=0.15`
(feature `OmitStations`, FLOWPanel `af92740` = tag
`campaign/p032-rootomit-20260918`). Banner gate on every start: `Particle
root-shed omission ACTIVE` with sane masked counts, plus each arm's
original banner checklist. Judge by outputs, not sacct. Campaign rules:
worktrees from annotated tags, env dev-paths at worktrees, run dirs MOVE
to `orc:~/projects/FLOWPanel.jl/data/` + symlink at harvest.

## Arm P1 — 022 hr10 IGE + omission (rank 1)

Question: does root-shed omission rescue the below-ground Γ-ignition
(13548847 died step 717/1008, rev 19.9, Γ-sum 0.25→618)? Success past
rev ~20 unparks the hr05/15/20 sweep.

Original executed submission (sacct SubmitLine, 2026-09-01, elapsed
1-05:19:26 to death):

```
sbatch --job-name=fp-022lg-hr10 --partition=m12 --qos=normal \
  examples/run_rotor_ground_effect_hpc.slurm.sh p022lg_hr10
```

ran from LIVE checkout `~/projects_unified/FLOWPanel.jl` + env
`~/projects_unified/envs/x86_64`, m12 CPU 64 threads. Banner:
mesh 45_185_ct4, RPM 6000, NT 36, h/R 1.0, disc 4R, panel 0.15R,
policy none, overlap 2.75, pps 12, merge_r 0.0055, das_uniform 3.4,
linegauss, settle 20, 26.5 revs / 1008 steps.

Rerun plan:
- Repo/env: p032 worktree + env (`~/campaigns/p032-rootomit-20260918/`),
  pins verbatim from `032_rerun_provenance_20260918.md`. No new tag
  needed unless pre-flight patches are required (then new tag).
- Backend match: keep m12 CPU / qos normal for a clean A/B vs 13548847
  (GPU-default preference noted, but the discriminator should be
  backend-matched; a GPU repeat can follow on success). Wall: original
  died at 29 h; full 1008 steps ≈ 36–40 h → request 48 h.
- Submission (draft): original line + `--time=48:00:00
  --export=ALL,P022_REPO_OVERRIDE/…,PARTICLE_OMIT_ROOT_R_OVER_R=0.15,
  RUN_NAME_OVERRIDE=p022lg_hr10_om15` — exact override-variable names
  TBD by pre-flight check 2.
- PRE-FLIGHT (blocking, ~30 min local):
  1. Verify `af92740` contains the 022 ground-effect driver + launcher
     state that 13548847 ran (it ran from projects_unified live checkout,
     not a pinned worktree — diff the driver/launcher against af92740).
  2. Verify the 022 driver reads `PARTICLE_OMIT_ROOT_R_OVER_R` (the 032
     commit wired the knob into the p018 screen driver; the 022 driver
     may need the same ~5-line patch → commit + new tag if so).
  3. Verify the launcher supports repo/project overrides and
     RUN_NAME_OVERRIDE (022 launcher predates the override convention).
  4. Feature smoke (5 steps, local or login): banner shows omission
     ACTIVE with ground plane present.
- Acceptance: survives past rev 20 and completes 26.5 revs with
  converged CT_IGE → fountain mechanism confirmed general incl. IGE;
  h/R sweep unparks with omission active. Dies ~step 650–720 again →
  IGE channel root-shed-independent; mechanism hunt reopens (escalate).

## Arm P2 — 018 NT72 `_3r` + omission (rank 2)

Question: does omission (np −30%, σ tail held 4–5× above floor) carry an
NT72 `_3r` rung past the FMM near-set adequacy gate (13592729 sfs3nb died
1550/2159; 13603725 merge2 died 1550; 13605984 cs0p18 died 1876)?

Baseline recipe = 13592729 (cleanest: explicit exports; 13605984 used
bare `--export=ALL` with shell-inherited env). Original executed
submission (sacct SubmitLine, 2026-09-05):

```
sbatch --job-name=fp-018gpu-n2_nt72_l3p0-3r-sfs3nb --partition=eng \
  --qos=eng --constraint=intel --gres=gpu:h200:1 --cpus-per-task=64 \
  --no-requeue --mem=192G --time=14:00:00 \
  --output=.../logs/slurm/slurm-%x-%j.out --error=.../slurm-%x-%j.err \
  --export=ALL,SIGMA_CEIL=0.030,TRUNCATION_RADIUS_R=3.0,\
MAX_PARTICLES=1500000,P018_SETTLE_REVS=22,SFS_THREELEVEL=true,\
P018_REPO_OVERRIDE=/home/rander39/wt018/FLOWPanel-pin-h200,\
P018_PROJECT_OVERRIDE=/home/rander39/p018wtenv-h200,\
P018_RUN_NAME=p018_csarc_n2_nt72_l3p0_3r_sfs3nb \
  examples/run_dji9443_hover_ct_gpu.slurm.sh h200 p018_csarc_n2_nt72_l3p0
```

Banner: mesh 45_185_ct4, NT 72, rlxf 0.16334 (exact-rate rule — KEEP),
sigma_chord 0.313, das_lambda 3.0, arc table
`p018_cs_l3p4_rs1_te_downwash_te.csv`, gaussian filament reg, N=2,
28.5 revs / 2160 steps.

Rerun plan:
- Delta vs original: `+PARTICLE_OMIT_ROOT_R_OVER_R=0.15`,
  `P018_RUN_NAME=p018_csarc_n2_nt72_l3p0_3r_sfs3nb_om15`,
  repo/project overrides → p032 worktree/env.
- Pins: FLOWPanel `af92740`; env = p032 env (FLOWVPM `8d4a3b4`,
  FastMultipole `ac7230a6`). NOTE the original ran the pre-§22 FLOWVPM;
  C-oldlaw showed the merge law is inert under omission, so the new-law
  env is acceptable — record this as a known model-def delta.
- PRE-FLIGHT (blocking):
  1. `git merge-base --is-ancestor f46c3fe af92740` in FLOWPanel — the
     smooth-ladders pin must be an ancestor so the case/driver features
     (das_arc, sigma_chord, SIGMA_CEIL, SFS_THREELEVEL) are all present
     at af92740; if any driver drift occurred between f46c3fe and
     af92740, diff the case defaults and record.
  2. Confirm the arc table CSV exists in the p032 worktree.
  3. SIGMA_CEIL=0.030 interaction: 032 arms ran SIGMA_CEIL=Inf. Keep the
     original 0.030 for backend-matching the adequacy-gate question
     (ceiling was part of that arm's model def) — flag for Ryan if he
     prefers Inf.
- Cost: eng H200, wall 14 h (original died at 2:26 @1550; full 2160
  likely ~4–6 h with np −30%).
- Acceptance: passes step 1550→1876 window and completes 2160 with M1
  settle → NT72 axis reopens; then M2 ε_Γ scoring per decision_rules
  (M1 AND M2 both required). Trips the same adequacy gate → σ-growth
  channel is root-shed-independent at NT72; NT144 stays parked.
- CONDITIONAL FOLLOW-ON (Ryan 2026-09-19): if this arm succeeds,
  offer NT144 rungs (exact-rate r=0.0854) as the next AskUserQuestion —
  not before.

## Arm P3 — 020 Phase-2R discriminator + omission (rank 4)

Question: is the Phase-2R tail-runaway (13154223, scr_p020r_geom_s020v,
died step 243/324) independent of root shed? This is a DISCRIMINATOR, not
an expected rescue (death is tail-localized field-coupled Γ-runaway, min
σ/σ_shed 0.0137, distinct from the root-fountain close-pair signature).
Phase 3 of 020 stays STOPPED regardless of outcome.

Original executed submission (sacct SubmitLine, 2026-08-12, m12 CPU,
elapsed 3:55:50 to death @243):

```
sbatch --job-name=fp-020r-geom examples/run_p018_screen_hpc.slurm.sh \
  scr_p020r_geom_s020v
```

ran from LIVE `~/projects/FLOWPanel.jl`. Banner: mesh 45_185_ct4, NT 36,
rlxf 0.3, σ/R 0.02 (sigma 0.0023736), overlap 2.4, pps 21, merge_r
0.00275, expint true, 8 revs / 324 steps, no freestream pulse.

Rerun plan:
- Delta: `+PARTICLE_OMIT_ROOT_R_OVER_R=0.15`,
  `RUN_NAME_OVERRIDE=scr_p020r_geom_s020v_om15`, repo/project overrides
  → p032 worktree/env (the p018 screen driver already reads the omission
  knob — proven by the B/C arms; CPU launcher variant
  `run_p018_screen_hpc.slurm.sh` shares the driver).
- Backend match: m12 CPU as original. Wall: 8 h (original died @243
  after ~4 h; full 324 steps ≈ 5–6 h).
- PRE-FLIGHT (blocking):
  1. Verify the corrected Phase-2R frozen-gradient/geometric integrator
     and its guards are present in the p032 pins (FLOWPanel af92740 +
     FLOWVPM 8d4a3b4) — the rig was built ~2026-08-12 against the then-
     live checkouts; if the integrator lives in uncommitted or since-
     diverged state, this arm needs its own pin ruling from Ryan.
  2. Confirm the CPU screen launcher supports REPO/PROJECT overrides +
     RUN_NAME_OVERRIDE (the GPU052 variant does).
  3. Check `scr_p020r_geom_s020v` case def still exists at af92740.
- Acceptance: survives step 243 and completes 324 → root shed co-drives
  the 020 testbed death; the closure-validation testbed may be viable
  under omission (report to Ryan; Phase 3 still gated). Dies ~243 with
  the same tail signature → mechanisms independent, 020 unchanged.

## Pre-flight results (2026-09-19, local checks at af92740)

- **P2 (018) is fully wired**: the 018 GPU launcher sources the CPU
  dispatcher, which runs `examples/rotor_hover_pressure_comparison.jl` —
  the driver that reads `PARTICLE_OMIT_ROOT_R_OVER_R` (proven live by the
  032 B/C arms). Arc table present at af92740
  (`data/p018_cs_l3p4_rs1_te_downwash_te.csv`). Local tag
  `campaign/p018-nt-20260905` (82a8e34) IS an ancestor of af92740; the
  scout-reported pin `f46c3fe` / tag `campaign/p018-smooth-ladders-20260905`
  does not exist in the local repo — verify on orc what commit
  `wt018/FLOWPanel-pin-h200` is detached at before claiming model-def
  equivalence (non-blocking: the rerun uses the p032 worktree anyway).
- **P3 (020) driver is wired, launcher needs a ~2-line patch**:
  `run_p018_screen_hpc.slurm.sh` reads the omission knob and
  `RUN_NAME_OVERRIDE`, and carries the `scr_p020r_geom_s020v` case def,
  but hard-codes `EXPECTED_REPO=/home/rander39/projects/FLOWPanel.jl`
  with no `P018_REPO_OVERRIDE`/`P018_PROJECT_OVERRIDE` — patch it to
  accept the overrides like the GPU052 wrapper, or route through the
  GPU052 wrapper (changes backend). Patch → commit → tag bump.
- **P1 (022) needs the most work, as anticipated**:
  `examples/rotor_hover_ground_effect.jl` does NOT read the omission
  knob (0 matches), and `run_rotor_ground_effect_hpc.slurm.sh` has only
  `P022_PROJECT_OVERRIDE` (no repo override; `P022_EXPECTED_REPO` guard
  defaults to the live checkout) and no `RUN_NAME_OVERRIDE`. Required:
  port the ~5-line omission-knob block from
  `rotor_hover_pressure_comparison.jl` into the 022 driver, add
  repo/run-name overrides to the launcher, commit, new tag
  (`campaign/p032-reopen-2202xxxx`), smoke with ground plane + banner.
- Integrator provenance for P3 (whether the Phase-2R frozen-gradient
  integrator survives intact in the p032 pins) is still UNVERIFIED —
  needs an orc-side diff of the 020 rig vs the pinned worktrees.

## Pre-flight results, round 2 (2026-09-19, post-reset session)

- **wt018 pin RESOLVED**: `wt018/FLOWPanel-pin-h200` is detached at
  `f46c3fe` = tag `campaign/p018-smooth-ladders-20260905` (tag exists on
  orc, not local). `f46c3fe` is **NOT an ancestor** of `af92740` — it sits
  on the unified-052/silo-snapshot lineage (merge-base `5615ada`; its
  parent `8f3ca07` is the exact 018-gpu silo snapshot). Model-def check
  done the direct way instead: the `p018_csarc_n2_nt72_l3p0` case-def
  export line is **byte-identical** at `f46c3fe` and `af92740`; driver
  drift between them is +266 lines of default-off feature additions
  (omission, split knobs, telemetry) + the launcher's path repointing.
  Record the rerun as model-def A vs the original with the FLOWVPM
  new-law delta (C-oldlaw showed the law inert under omission) and the
  driver generation noted.
- **P2 arc table**: NOT present in the p032 worktree (`data/` there hides
  tracked fixtures); fixed by symlinking
  `~/projects/FLOWPanel.jl/data/p018_cs_l3p4_rs1_te_downwash_te.csv` into
  `~/campaigns/p032-rootomit-20260918/FLOWPanel.jl/data/` (2026-09-19).
- **P3 integrator provenance RESOLVED**: the Phase-2R frozen-gradient
  integrator IS the `WAKE_EXPINT=true` euler_exp path — present at
  `af92740` (driver reads the knob; rk3/expint mutual-exclusion guards
  intact) and in FLOWVPM `8d4a3b4` (exercised live by 032's
  B15/B-floor/C-oldlaw arms). No separate pin ruling needed.
- **P3 modern-stack controls already exist**: the 026 sec.16 reruns
  (2026-09-05, current stack) reproduce the expint death without
  omission — `scr_p026ef_exp_s020v` (vatistas) blew up at step 211
  (CF ~1e11 at 211), `scr_p026ef_exp_s020v_lg` ended ~step 300. Score
  the om15 discriminator against these, not only against 13154223.
  NOTE: the CPU launcher now defaults linegauss — the P3 rerun must
  export `FLOWPANEL_FILAMENT_REG=vatistas` to backend-match the original
  and the ef control.
- **P3 launcher patch APPLIED (uncommitted)**:
  `run_p018_screen_hpc.slurm.sh` now maps
  `P018_REPO_OVERRIDE`/`P018_PROJECT_OVERRIDE` → `P018_REPO`/`P018_PROJECT`
  (2-line block; the run line already honored `P018_PROJECT`).
- **P1 driver port APPLIED (uncommitted)**: `PARTICLE_OMIT_ROOT_R_OVER_R`
  knob + mask/`pnl.OmitStations` wrap ported from
  `rotor_hover_pressure_comparison.jl` into
  `rotor_hover_ground_effect.jl` (gated to CONVERSION=legacy AND
  NROTORS=1 — station radii are measured on the rotor-1 node set, wrong
  for offset rotors); launcher banner gains `omit_root:`. Feature smoke
  (local, 40_40 mesh, ground plane on, om15): banner
  `Particle root-shed omission ACTIVE: shedding{1,2} omits 2/36 stations`,
  ground `h/R=1.0 ncells=4752`. First smoke attempt died at step 0 on the
  DRIVER-DEFAULT `GS_TOL=1e-8`/`GS_MAX_OUTER=50` (residual 4.2e-8 — the
  documented FMM metric floor, unrelated to omission; the hr10 case def
  sets 1e-5/100); second attempt showed `BERNOULLI_ONLY=true` is required
  with GROUND_ENABLE (PressureLaplace monitor is built for 1 body; the
  launcher always sets it). Final smoke with case-matched GS +
  BERNOULLI_ONLY: **PASS** — omission ACTIVE 2/36 both blades, ground
  h/R=1.0, 3/3 steps clean (~14 s/step local). **P1 smoke gate CLEARED.**
- **P1 launcher needs NO patch**: `P022_EXPECTED_REPO`,
  `P022_PROJECT_OVERRIDE`, and `P022_RUN_NAME` already exist — the rerun
  overrides all three at submission.
- **P1 driver drift vs 13548847 is inert**: diff of the 09-01 state
  (`8413697`) vs `af92740` for the 022 driver = behavior-identical
  `clip_shedding_root` refactor + `sigma_guard` now passed to euler_exp
  (inert: `p022lg_hr10` does not set WAKE_EXPINT) + launcher default
  paths (overridden anyway).
- **Commit/tag needed** (Ryan-gated): one commit bundling the CPU-launcher
  override patch + the 022 omission port, tagged
  `campaign/p032-reopen-20260919`; new orc worktree + env for P3/P1
  (P2 can run from the existing p032 worktree unchanged).

## Suggested execution order

P3 (cheapest, ~6 h CPU) and P2 (one 14 h GPU slot) can go first and in
parallel; P1 is the most expensive (≈2 d CPU) and has the most pre-flight
risk (022 driver knob wiring). All submissions Ryan-gated after
pre-flight results are reported.
