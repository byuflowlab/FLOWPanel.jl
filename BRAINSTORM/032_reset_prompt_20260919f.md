# Reset prompt — 032 reopen arms LAUNCHED; babysit + finish deferred items (2026-09-19f)

You are picking up work in
`/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (branch
**fastmultipole**, HEAD `7b49168` = tag `campaign/p032-reopen-20260919`)
+ sibling `/Users/ryan/Dropbox/research/projects/FLOWVPM.jl` (branch
`flowpanel`, `8d4a3b4`). Read `CLAUDE.md` + the policies it names.
Supersedes `032_reset_prompt_20260919e.md` — that prompt's Tasks 1+2
were EXECUTED this session (Ryan approved all via one AskUserQuestion;
he additionally authorized HPC jobs and said **use rsync for now, no
git pushes yet**). Read for context, do NOT re-derive:
`BRAINSTORM/032_omission_reopen_prep_20260919.md` (acceptance criteria
per arm), `BRAINSTORM/032_followup_provenance_20260919.md`.

## DONE this session (2026-09-19 afternoon MDT)

**Task 1 (merge-law finish), partially:**
- orc `~/projects/FLOWVPM.jl` ff-merged to **`d896145`** (test suite had
  PASSED on the merge tree pre-session).
- Tag `archive/flowpanel-wip-20260902` created on orc at `eaf257c`.
- **DEFERRED by Ryan**: `git push --force-with-lease orc flowpanel` from
  local FLOWVPM (all pushes deferred — rsync era).
- **NOT DONE — cleanup** (gate: `ssh orc 'pgrep -u rander39 -f
  vpm_testenv'` empty; duplicate Pkg.test was still alive ~12:10 MDT,
  expected to self-exit; my session monitor died with the session):
  `ssh orc 'cd ~/projects/FLOWVPM.jl && git worktree remove /tmp/rander39_vpm_mergetest --force; rm -rf /tmp/rander39_vpm_testenv /tmp/rander39_vpm_test_20260919d.log'`

**Task 2 (omission re-open), fully launched:**
- Local commit **`7b49168`** (the three launcher/driver edits: p018
  screen `P018_REPO_OVERRIDE`/`P018_PROJECT_OVERRIDE` mapping;
  `rotor_hover_ground_effect.jl` PARTICLE_OMIT_ROOT_R_OVER_R +
  `pnl.OmitStations`; ground launcher `omit_root:` banner), annotated
  tag **`campaign/p032-reopen-20260919`**. NOT pushed anywhere (Ryan).
- Reopen checkout staged on orc **via rsync of a `git archive` tarball**
  (not a git worktree — no push): tree at
  `orc:~/campaigns/p032-reopen-20260919/FLOWPanel.jl` with
  `PROVENANCE_PIN.txt` ("campaign/p032-reopen-20260919 (7b49168)"),
  `data/`, `logs/slurm/`. Env at `.../p032-reopen-20260919/env` =
  byte-copy of the p032 env Project.toml + Manifest.toml with the
  FLOWPanel dev-path sed-repointed. Verified under the JOB Julia (see
  gotcha): pathof → reopen tree, `OmitStations=true`, FLOWVPM →
  `~/campaigns/p026-derisk-20260914/FLOWVPM.jl` @ `8d4a3b4` clean.
- **All three arms RUNNING with verified banners** (omission ACTIVE
  3/41 stations BOTH blades on every arm; run dirs relocated to the
  consolidated root `/home/rander39/projects/FLOWPanel.jl/data/` with
  symlinks back into the worktree data/):

| Arm | Job | Where | Banner verified |
|---|---|---|---|
| P3 `scr_p020r_geom_s020v_om15` | **13774448** (m12, 8 h, ~16:47 start) | reopen wt | vatistas pinned, expint:true, s020v knobs (overlap 2.4, pps 21, merge_r 0.00275, nrevs 8) |
| P2 `p018_csarc_n2_nt72_l3p0_3r_sfs3nb_om15` | **13774449** (eng H200, 14 h, ~16:40 start) | p032 wt (`af92740`) | NT72, rlxf 0.16334, sigma_chord 0.313, arc table, SFS 3-level NB, **SIGMA_CEIL=Inf (Ryan chose Inf over 0.030)** |
| P1 `p022lg_hr10_om15` | **13776680** (m12, 48 h, ~17:2x start; RETRY) | reopen wt | omit_root 0.15, linegauss, h/R=1.0, GS_TOL 1e-5, ground ncells 4752 |

- First P1 attempt **13775310 DIED in precompile** (OpenSSL_jll→HDF5→
  FLOWVPM). ROOT CAUSE + GOTCHA: orc login default is Julia **1.12.7**
  (juliaup); jobs use the spack module **julia/1.11.7-6bmogfl**. An
  initial `Pkg.develop` under 1.12 re-resolved/polluted the reopen
  Manifest. Fixed by re-copying the p032 Manifest + sed path repoint,
  verified with
  `/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin/julia`.
  **Do ALL env ops with that binary; never let 1.12 touch a campaign
  Manifest.** (`import HDF5` directly in the env fails benignly — HDF5
  is transitive via FLOWVPM, not a project dep.)
- Backend note (Ryan asked): P3/P1 on CPU is deliberate backend-matching
  to their CPU controls; P2 on H200 matches its GPU comparator. A GPU
  P1-class pair would need a GPU control rerun — offered, not requested.

**Storage:** hpc-storage agent dispatched for the two APPROVED actions:
(a) archive `scr_p032om15_ctrllg_fs` + `scr_p032om15_explg_fs_oldlaw`
(~31 G, `./scripts/run_archiver.sh --only ... --apply` from
`~/projects/FLOWPanel.jl`), (b) the **212 GB `p018_csarc_*_3r_*`
reclaim (Ryan APPROVED)**. Agent was still mid-measurement (slow du) at
session end — VERIFY completion first thing (`hpc-monitor`/ledger
lines; archiver logs under the usual archive ledger). Unresolved
oddity: `df -h /home/rander39` showed 193 G used on a 2.0 T VAST mount
vs the morning's 396.4/400 quota figure — reconcile before trusting
either.

## NEXT ACTIONS

1. Verify hpc-storage finished both actions; re-dispatch if it died.
2. Babysit the three arms (delegate to `hpc-monitor`): judge by outputs
   (never sacct). Acceptance criteria per arm = prep file
   `032_omission_reopen_prep_20260919.md`. P3 control context: ready-made
   controls died step 211 (`scr_p026ef_exp_s020v`) / ~300 (`_lg`).
   P1 control = `p022lg_hr10` 13548847 physical blow-up (below-ground
   Γ-ignition). On any death: tail log, classify (Γ-ignition vs infra).
3. Task-1 cleanup once `pgrep -f vpm_testenv` is empty (commands above).
4. Harvest owed: C-ctrl **13773490** + C-oldlaw **13773491** (launched
   09-19 morning, see `032_followup_provenance_20260919.md`).
5. NT144 offer ONLY after a P2 success verdict.

## Owed / parked (carried, all Ryan-gated)

- Docs commit bundle: 026 ledger line, provenance edits (followup +
  reopen-prep files, d/e/f reset prompts), 032 item-file Log updates,
  INDEX.md rows (032 refresh + missing `031_quadrupole_panel_farfield.md`).
- ALL pushes (github + orc, incl. the deferred FLOWVPM `flowpanel`
  force-with-lease) — Ryan said "don't worry about pushes yet".
- Notebook entries (×4 from 021 + this arc's launches) — offer, don't write.
- k=3 cap retune + merged-σ/clamp telemetry (A4); 021 silo cleanup;
  `scr_p026gpuv_split` retry once quiet ≥24 h; `data/scratch_p032_smoke*`
  disposable; 026 NT144 un-park + feature-A move-vs-delete DECLINED.
- 021 entry point remains `fgs_acceleration_reset_prompt_20260919e.md`
  (evaluator+A/B Ryan-gated).

## Ground rules

Local ≤4 threads. `ssh orc` needs a live ControlMaster socket (`! ssh
orc echo ok` if 2FA blocks); `-o ControlMaster=no` for one-shots;
`bash -lc` for Slurm; strip MOTD/ANSI. Submissions/commits/remote git
state/notebook writes Ryan-gated (HPC job submission for THIS arc was
pre-authorized and is done; new submissions need a fresh ask).
Pre-existing dirty files (018/026 docs, rotor_multi slurm script,
pressure-comparison TOML, this arc's BRAINSTORM docs) are expected —
commit only with Ryan's approval. Judge runs by outputs, never sacct.
