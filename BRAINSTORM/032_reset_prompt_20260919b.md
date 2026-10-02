# Reset prompt — 026 next-runs recommendation + 032 notebook entry (2026-09-19b)

You are picking up work in
`/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (branch
**fastmultipole**, HEAD `af92740` = tag `campaign/p032-rootomit-20260918`)
+ sibling `/Users/ryan/Dropbox/research/projects/FLOWVPM.jl` (branch
`flowpanel`). Read `CLAUDE.md` + the policies it names (`HPC.md` +
`BYU_ORC_AGENTS.md` before ORC work). Supersedes
`032_reset_prompt_20260919.md` (all its tasks DONE). Read these two
records FIRST — they carry the complete 032/026 state; do NOT re-derive:

- `BRAINSTORM/032_rerun_provenance_20260918.md` — pins, arms, banner/
  outcome tables, FINAL 032 VERDICT, fountain before/after, ΔCT vs
  baselines, merge-burden analysis.
- `BRAINSTORM/026_sigma_growth_particle_splitting/rerunslate_provenance_20260918.md`
  — slate verdict (all 5 arms dead/blown; redesign rescues nothing),
  harvest table.

## State in one paragraph

032 is CLOSED as an A/B: root-shed omission (`OmitStations`, knob
`PARTICLE_OMIT_ROOT_R_OVER_R`) rescues every tested exp-family s9 class.
B15/B20 (A1 config, knobs 0.15/0.20, jobs 13771353/54) and B-floor (A2
config, 0.15, job 13772015) all completed 467/467 CONVERGED, CT
0.0713–0.0717 (tightest plateaus in the campaign: ±0.014–0.144%);
un-omitted variants died @274–328 in two waves. Mechanism: σ-at-floor
Γ-concentrated close pairs never form (min σ ≥4.6× floor, max Γ/σ²
46–143 vs A1's 9.1e4); np peaks ~370k (A1 rode the 500k cap); merge
frenzy collapses ~30× (2.2–3k events vs 84–87k). Deleted circulation
accounted: 0.486/0.722/0.487 m³/s. All run dirs at
`orc:~/projects/FLOWPanel.jl/data/scr_p032om{15,20}_explg_fs`,
`.../scr_p032om15_explg_floor` (+symlinks from the p032 worktree).
Local ParaView staging: `~/scr_p026s9r2_explg_fs_last50steps/` (A1,
WITH root shed, steps 278–327) vs `~/scr_p032om15_explg_fs_steps278_327/`
(B15, omitted, same steps; particles+body VTPs/VTUs, monitors, metadata).

## TASK 1 — recommend the next 026 HPC runs (root-shedding omitted)

Form a ranked recommendation and put it to Ryan via AskUserQuestion
(multiSelect) BEFORE submitting anything. Candidates already identified
(032 provenance "Merge-burden analysis" + verdict sections; add your own
judgment):

- **ctrl-channel discriminator**: A3 config (`scr_p026s9_ctrllg_fs`,
  WAKE_EXPINT=false) + omission 0.15 — does root shed also drive the
  pump-independent ctrl blow-up (onset steps 253–288, np collapsed to
  17k)? Highest information value.
- **old-merge-law + omission A/B**: A1 config + omission 0.15 with the
  OLD FLOWVPM pin (pre-`8d4a3b4`, wave-2 tag `campaign/p026-wave2-20260915`
  = FLOWVPM `2b253db`) — definitive "was the merge redesign necessary
  once root shed is gone"; prediction on record: near-null. NOTE: needs
  its own FLOWVPM worktree+env (different pin!) — new campaign tag/env
  per Campaign Reproducibility rules.
- **NT144 cap-ladder rung un-park**: omission's ~30% lower np relieves
  the count-driven FMM σ-adequacy gate that killed every NT72 rung —
  probe whether NT ladders now run at exp settings. (Case
  `scr_p026sp_nt144_cap030`-family; check case defs for what a rung
  needs; exact-rate NT ladder rule r(NT)=1−(1−0.3)^(36/NT) applies to
  ladders per Ryan 2026-08-20.)
- **feature-A "move vs delete"**: A1 config + `SHEDDING_R_OVER_R=0.2`
  (edge omission relocates the chain-closing filament outboard instead
  of deleting it) — mechanism discriminator, parked queue-budget call.
- **un-omitted survivable baseline for clean ΔCT** — only if Ryan wants
  ΔCT attribution beyond the dose-response bound (−0.42% per +48%
  deletion).

## TASK 2 — launch the approved set (campaign rules)

- New-code arms can REUSE the p032 pins/worktree/env
  (`orc:~/campaigns/p032-rootomit-20260918/{FLOWPanel.jl,env}`, FLOWPanel
  `af92740`, FLOWVPM `8d4a3b4`, FMM `ac7230a6`) — FLOWPanel is unchanged,
  no new tag needed. The OLD-law arm is the exception: pin FLOWVPM
  `2b253db` in a fresh worktree + its own env (tag convention
  `campaign/p032-<slug>-20260919`), record in a provenance file BEFORE
  submitting.
- Submission pattern (exactly what B15/B20/B-floor used — copy from the
  032 provenance Submissions section): cwd = p032 FLOWPanel worktree,
  `sbatch --time=12:00:00 --export=ALL,P018_REPO_OVERRIDE=…,P018_PROJECT_OVERRIDE=…,
  FLOWPANEL_FILAMENT_REG=linegauss,<mech env>,RUN_NAME_OVERRIDE=<run>,
  PARTICLE_OMIT_ROOT_R_OVER_R=0.15
  ~/projects/FLOWPanel.jl/examples/run_p018_screen_gpu052.slurm.sh h200 <case>`;
  m13h default partition (add `-p eng --qos=eng --gres=gpu:h200:1` for
  eng). FLOOR env = `SIGMA_FLOOR_FRAC=0.1`; SPLIT trio =
  `WAKE_SPLIT_VISCOUS=true WAKE_SPLIT_FRAC_VISCOUS=0.587
  WAKE_SPLIT_FRAC_COMPRESS=0.73 WAKE_SPLIT_FRAC_ELONGATE=0.3`; ctrl arm
  = fs env on the ctrllg case (WAKE_EXPINT=false comes from the case def).
  Backend-match any arm that A/Bs a specific predecessor.
- Banner verification owed on start (omission ACTIVE + masked counts +
  wave-2 checklist; `CONVERSION=legacy`). Run dirs land worktree-local
  (launcher `rm -rf`s pre-made symlinks) — MOVE + symlink at harvest.
  Judge by outputs, not sacct; tracebacks in the `.err`.

## TASK 3 — notebook entry summarizing 032 (Ryan-directed, but STILL ask before writing)

Notebook rules in `~/.claude/CLAUDE.md` apply: open the newest
`~/Dropbox/research/notebooks/journals/YYYYMMDD.md`; if no `# 20260919`
header exists, PROPOSE the day header + checklist and WAIT for approval;
never tick checkboxes; append-only.

Ryan's spec for this entry (2026-09-19): **less than 50 words of prose
TOTAL** covering context, methods, analysis, conclusions — tables/figure
captions don't count as prose, so lean on: one compact results table
(arm | knob | fate | CT | deleted Σ|Γ|·Δl; un-omitted twins died 274–328)
and the GIF. Plus **a GIF of a ParaView visual showing with-root-shedding
and omitted side-by-side**:

- Frame data is already staged locally, same steps 278–327:
  `~/scr_p026s9r2_explg_fs_last50steps/` (A1, WITH) vs
  `~/scr_p032om15_explg_fs_steps278_327/` (B15, OMITTED). Particle VTPs
  carry `velocity_gradient` (9-comp) and Gamma; the ignition knot is at
  x=−0.37R, r=0.34R (rotor axis = x, thrust −x; R=0.119 m). Color by
  |Γ|/σ² or |Γ|; fix the SAME color scale + camera on both views
  (side-by-side layout), annotate step number.
- Render via `pvpython` (check `which pvpython` / the ParaView.app
  bundle path on macOS) writing PNG frames per step, then assemble
  (`ffmpeg -f image2 … .gif` or ImageMagick `convert -delay`). ≤4 local
  threads. If pvpython is unavailable, fall back to asking Ryan to
  screen-record, or build the side-by-side from raw VTP parsing +
  matplotlib scatter (VTPs are appended-binary; meshio canNOT read them
  — a manual parser precedent exists from the A1 offender localization).
- Store as `~/Dropbox/research/notebooks/img/20260919_p032_root_omission/
  p032_sidebyside.gif` (+ the render script in the same dir), reference
  by relative path `../img/20260919_p032_root_omission/…`. Offer Ryan
  the draft entry text + GIF for approval BEFORE writing anything into
  the journal file.

## Owed / parked (carried — do not silently drop)

- hpc-storage archiving: 026 slate dirs (40.3 G) + 032 dirs (45 G) once
  Ryan's ParaView pass / the GIF render no longer needs them;
  `scr_p026gpuv_split` retry once quiet ≥24 h.
- Uncommitted by design (Ryan-gated commit bundle): 026 ledger line,
  both provenance-file edit sets, 032 item-file Log updates, this file.
- github pushes of branches+tags; orc branch-divergence ruling.
- k=3 cap retune + merged-σ/clamp-engagement telemetry gap (A4 finding).
- 021 silo cleanup; `031_quadrupole_panel_farfield.md` missing INDEX row.
- Local scratch: `data/scratch_p032_smoke{A,B,C,D}` disposable;
  `~/scr_p032om15_explg_fs_steps278_327/` (2.3 G) keep until GIF done.

## Ground rules

Local ≤4 threads. `ssh orc` needs a live ControlMaster socket (`! ssh orc
echo ok` if 2FA blocks); `bash -lc` for Slurm; MOTD + ANSI escapes
contaminate ssh output (strip before parsing — a `^`-anchored grep broke
once). All submissions beyond Task-2 approvals, all commits, and all
notebook writes are Ryan-gated. Pre-existing dirty files (018/026 docs
not from this arc, rotor_multi slurm script, data TOML) are NOT yours.
