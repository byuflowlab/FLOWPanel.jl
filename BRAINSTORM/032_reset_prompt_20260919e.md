# Reset prompt — merge-law merge finish (Ryan-gated) + omission re-open submissions (2026-09-19e)

You are picking up work in
`/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (branch
**fastmultipole**, HEAD `af92740` = tag `campaign/p032-rootomit-20260918`,
plus UNCOMMITTED edits listed below) + sibling
`/Users/ryan/Dropbox/research/projects/FLOWVPM.jl` (branch `flowpanel`,
HEAD `8d4a3b4`). Read `CLAUDE.md` + the policies it names. Supersedes
`032_reset_prompt_20260919d.md` (Task 1 VERIFIED but blocked on
permissions; Task 2 pre-flight COMPLETE). Read FIRST — do NOT re-derive:

- `BRAINSTORM/032_omission_reopen_prep_20260919.md` — now carries a
  "Pre-flight results, round 2 (2026-09-19, post-reset session)" section
  with every resolved item (wt018 pin = `f46c3fe` on the silo lineage,
  NOT an ancestor of af92740, but P2 case-def line byte-identical; P3
  integrator = WAKE_EXPINT path, present in pins; ready-made P3 controls
  `scr_p026ef_exp_s020v` died step 211 / `_lg` ~300; P1 driver drift
  inert; arc-table symlink created in the p032 worktree).
- `BRAINSTORM/032_followup_provenance_20260919.md` — C-ctrl/C-oldlaw
  verdicts (omission general rescue; merge law inert under omission).

## TASK 1 — merge-law production merge: VERIFIED, finish is Ryan-gated

The FLOWVPM test suite **PASSED on the merge tree** (`d896145`,
`mergetest-unified052`, worktree `orc:/tmp/rander39_vpm_mergetest`): the
pre-reset background run survived the session — final line
`Testing FLOWVPM tests passed`, radix/SFS suites green, CUDA subtests
skipped off-GPU (expected). Queue was empty. The 2026-09-19d session's
attempt to execute the finish steps was **blocked by the permission
classifier** (remote git state changes over ssh). Bring them to Ryan
(they are part of the single AskUserQuestion below); once approved, run
or have Ryan run:

1. `ssh orc 'cd ~/projects/FLOWVPM.jl && git merge --ff-only mergetest-unified052'`
   (live checkout is on `unified-052` @ `3315b22`, tracked-clean).
2. Tag orc's stale `flowpanel` ref (at `eaf257c`):
   `ssh orc 'cd ~/projects/FLOWVPM.jl && git tag -a archive/flowpanel-wip-20260902 eaf257c -m "wip snapshot superseded by d896145 merge"'`,
   then from LOCAL FLOWVPM.jl: `git push --force-with-lease orc flowpanel`.
3. Cleanup — ONLY after the duplicate test exits (see below):
   `ssh orc 'cd ~/projects/FLOWVPM.jl && git worktree remove /tmp/rander39_vpm_mergetest --force; rm -rf /tmp/rander39_vpm_testenv /tmp/rander39_vpm_test_20260919d.log'`.

NOTE: the d-session, assuming the original test had died, launched a
DUPLICATE `Pkg.test` on the orc login node
(`julia --project=/tmp/rander39_vpm_testenv`, log
`/tmp/rander39_vpm_test_20260919d.log`); kill attempts were
classifier-blocked. It is harmless (sandboxed, 4 threads) and will exit
on its own ~1–2 h after 11:22 orc time — do not re-diagnose it, just
gate cleanup step 3 on `pgrep -u rander39 -f vpm_testenv` being empty.

## TASK 2 — omission re-open: pre-flight DONE, everything left is Ryan-gated

UNCOMMITTED edits in the LOCAL live checkout (yours to commit once Ryan
approves; verify with `git diff --stat`):

- `examples/run_p018_screen_hpc.slurm.sh` — maps
  `P018_REPO_OVERRIDE`/`P018_PROJECT_OVERRIDE` → `P018_REPO`/`P018_PROJECT`
  (P3 enabler) + this session's `omit_root:` banner echo was added to the
  022 launcher, see next lines.
- `examples/rotor_hover_ground_effect.jl` — PARTICLE_OMIT_ROOT_R_OVER_R
  knob + `pnl.OmitStations` wrap ported from the p018 driver (gated
  CONVERSION=legacy AND NROTORS=1).
- `examples/run_rotor_ground_effect_hpc.slurm.sh` — banner `omit_root:`.
- `BRAINSTORM/032_omission_reopen_prep_20260919.md` — round-2 results.
- This file.

P1 feature smoke (local, 40_40 mesh, ground on, om15): banner PASS
(`omission ACTIVE: shedding{1,2} omits 2/36 stations`, ground h/R=1.0
ncells=4752). Two setup gotchas hit and understood (documented in the
prep file): driver-default GS_TOL=1e-8 trips the documented FMM metric
floor (case def uses 1e-5/100 — not an omission bug), and
BERNOULLI_ONLY=true is required with GROUND_ENABLE (PressureLaplace is
built for 1 body; the launcher always sets it). Final smoke **PASSED**:
omission ACTIVE 2/36 both blades, ground h/R=1.0, 3/3 steps clean. The
P1 smoke gate is CLEARED — no smoke work remains.

### The single AskUserQuestion to bring Ryan (do this first)

Bundle these decisions:

1. **Task-1 finish** — approve the three command groups above (or Ryan
   runs them himself with `! ssh orc ...`).
2. **Commit + tag** — one commit on fastmultipole bundling the three
   uncommitted example/launcher edits, annotated tag
   `campaign/p032-reopen-20260919`; push branch+tag to the orc FLOWPanel
   remote; on orc create worktree
   `~/campaigns/p032-reopen-20260919/FLOWPanel.jl` detached at the tag +
   env copy of the p032 env with the FLOWPanel dev-path repointed
   (FLOWVPM/FastMultipole dev-paths unchanged = `8d4a3b4`/`ac7230a6`);
   verify `pathof` + `isdefined(FLOWPanel, :OmitStations)`; give the
   worktree a `data/` with the run-dir conventions. P2 needs NONE of
   this — it runs from the existing p032 worktree.
3. **Submissions** (suggested order: P3 + P2 in parallel, P1 after; all
   paths literal `/home/rander39`):

P3 (m12 CPU, backend-matched incl. explicit vatistas; submit from the
reopen worktree):

```
sbatch --job-name=fp-020r-geom-om15 --partition=m12 --qos=normal --time=8:00:00 \
  --export=ALL,P018_REPO_OVERRIDE=/home/rander39/campaigns/p032-reopen-20260919/FLOWPanel.jl,P018_PROJECT_OVERRIDE=/home/rander39/campaigns/p032-reopen-20260919/env,FLOWPANEL_FILAMENT_REG=vatistas,PARTICLE_OMIT_ROOT_R_OVER_R=0.15,RUN_NAME_OVERRIDE=scr_p020r_geom_s020v_om15 \
  examples/run_p018_screen_hpc.slurm.sh scr_p020r_geom_s020v
```

P2 (eng H200; submit from the p032 worktree, no new code):

```
sbatch --job-name=fp-018gpu-n2_nt72_l3p0-3r-sfs3nb-om15 --partition=eng --qos=eng \
  --constraint=intel --gres=gpu:h200:1 --cpus-per-task=64 --no-requeue --mem=192G --time=14:00:00 \
  --output=/home/rander39/campaigns/p032-rootomit-20260918/FLOWPanel.jl/logs/slurm/slurm-%x-%j.out \
  --error=/home/rander39/campaigns/p032-rootomit-20260918/FLOWPanel.jl/logs/slurm/slurm-%x-%j.err \
  --export=ALL,SIGMA_CEIL=0.030,TRUNCATION_RADIUS_R=3.0,MAX_PARTICLES=1500000,P018_SETTLE_REVS=22,SFS_THREELEVEL=true,PARTICLE_OMIT_ROOT_R_OVER_R=0.15,P018_REPO_OVERRIDE=/home/rander39/campaigns/p032-rootomit-20260918/FLOWPanel.jl,P018_PROJECT_OVERRIDE=/home/rander39/campaigns/p032-rootomit-20260918/env,P018_RUN_NAME=p018_csarc_n2_nt72_l3p0_3r_sfs3nb_om15 \
  examples/run_dji9443_hover_ct_gpu.slurm.sh h200 p018_csarc_n2_nt72_l3p0
```

P1 (m12 CPU 48 h; submit from the reopen worktree):

```
sbatch --job-name=fp-022lg-hr10-om15 --partition=m12 --qos=normal --time=48:00:00 \
  --export=ALL,P022_EXPECTED_REPO=/home/rander39/campaigns/p032-reopen-20260919/FLOWPanel.jl,P022_PROJECT_OVERRIDE=/home/rander39/campaigns/p032-reopen-20260919/env,P022_RUN_NAME=p022lg_hr10_om15,PARTICLE_OMIT_ROOT_R_OVER_R=0.15 \
  examples/run_rotor_ground_effect_hpc.slurm.sh p022lg_hr10
```

4. **SIGMA_CEIL for P2** — draft keeps the original 0.030
   (backend-matches the adequacy-gate question); flag that 032 arms ran
   Inf if Ryan prefers.
5. **Storage** — /home at 396.4/400 G. (a) Re-dispatch hpc-storage for
   the two new arms (`scr_p032om15_ctrllg_fs`,
   `scr_p032om15_explg_fs_oldlaw`, ~31 G): they were still ~40 min short
   of the 2 h quiet window at 11:36 MDT — clear after ~12:16 MDT;
   archiver: `cd /home/rander39/projects/FLOWPanel.jl && ./scripts/run_archiver.sh --only scr_p032om15_ctrllg_fs,scr_p032om15_explg_fs_oldlaw --apply`.
   (b) The 212 GB `p018_csarc_*_3r_*` reclaim (quiet ≥24 h, unprotected)
   is still a PENDING Ryan decision — the P2 rerun will add ~10–20 G, so
   surface it in the same ask.

Banner gates on every start: omission ACTIVE with sane masked counts
(expect 3/41 per blade on 45_185_ct4), each arm's original checklist
(P3: expint true, vatistas, s020v knobs; P2: NT72 rlxf 0.16334
sigma_chord 0.313 arc table SFS_THREELEVEL SIGMA_CEIL=0.030; P1:
linegauss, h/R 1.0, GS_TOL 1e-5, omit_root in launcher banner).
Acceptance criteria per arm: see the prep file. NT144 offer ONLY after a
P2 success verdict.

## Owed / parked (carried)

- Ryan-gated commit bundle: 026 ledger line, provenance edits (incl.
  `032_followup_provenance_20260919.md`, the reopen-prep file, the d and
  e reset prompts), 032 item-file Log updates; INDEX.md rows for 032
  outcome refresh + missing `031_quadrupole_panel_farfield.md` row.
- github pushes of branches+tags (orc pushes partly pending via Task 1).
- k=3 cap retune + merged-σ/clamp telemetry gap (A4); 021 silo cleanup;
  `scr_p026gpuv_split` retry once quiet ≥24 h;
  `data/scratch_p032_smoke{A,B,C,D}` disposable. 026 NT144 un-park and
  feature-A move-vs-delete remain DECLINED.
- Sweep verdicts for the record: unchanged from the d prompt (020
  tail-seated field-coupled runaway; 005/006 truncation-seated NO
  re-open; 017 500k-cap-blocked; 018 Phase-16 λ-ladder TIMEOUT+M2).

## Ground rules

Local ≤4 threads. `ssh orc` needs a live ControlMaster socket (`! ssh
orc echo ok` if 2FA blocks); `bash -lc` for Slurm; strip MOTD/ANSI from
ssh output; ssh commands that spawn remote background children hold the
channel open — use `-o ControlMaster=no` for one-shots. All submissions,
commits, remote git state changes, and notebook writes Ryan-gated (the
permission classifier also enforces this — do not fight it). Pre-existing
dirty files (018/026 docs not from this arc, rotor_multi slurm script,
data TOML) are NOT yours. Judge runs by outputs, never sacct.
