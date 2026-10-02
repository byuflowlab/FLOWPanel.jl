# Reset prompt — BRAINSTORM omission-sweep + 032 follow-up harvest (2026-09-19c)

You are picking up work in
`/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (branch
**fastmultipole**, HEAD `af92740` = tag `campaign/p032-rootomit-20260918`)
+ sibling `/Users/ryan/Dropbox/research/projects/FLOWVPM.jl` (branch
`flowpanel`). Read `CLAUDE.md` + the policies it names (`HPC.md` +
`BYU_ORC_AGENTS.md` before ORC work). Supersedes
`032_reset_prompt_20260919b.md` (all three tasks DONE: follow-up arms
launched, notebook entry written). Read these records FIRST — do NOT
re-derive state:

- `BRAINSTORM/032_followup_provenance_20260919.md` — the two live arms
  (pins, recipes, banner PASS tables, acceptance criteria).
- `BRAINSTORM/032_rerun_provenance_20260918.md` — FINAL 032 VERDICT,
  B15/B20/B-floor outcomes, merge-burden analysis, ΔCT bounds.

## State in one paragraph

032 is CLOSED as an A/B: root-shed omission
(`PARTICLE_OMIT_ROOT_R_OVER_R=0.15`, feature `OmitStations`) rescues
every tested exp-family s9 class (B15/B20/B-floor 467/467 CONVERGED, CT
0.0713–0.0717) while every un-omitted variant died @253–328 across two
waves. Mechanism: the σ-at-floor Γ-concentrated fountain close pairs
never form; np peaks ~370k instead of riding the 500k cap; merge frenzy
collapses ~30×; omission also relieves the count-driven FMM σ-adequacy
gate (the killer of every 018 NT72 rung). Two Ryan-approved follow-up
arms are RUNNING (submitted 2026-09-19 ~09:00, ~71–75 min expected):
**C-ctrl** = job 13773490 (A3 ctrl config + omission 0.15, eng, run
`scr_p032om15_ctrllg_fs`) and **C-oldlaw** = job 13773491 (A1 config +
omission 0.15 on OLD FLOWVPM `2b253db`, m13h, run
`scr_p032om15_explg_fs_oldlaw`, own worktree/env
`orc:~/campaigns/p032-oldlaw-20260919/`, tag
`campaign/p032-oldlaw-20260919`). Banners PASS on both. Declined this
round (do not resubmit without Ryan): feature-A move-vs-delete, NT144
un-park, un-omitted ΔCT baseline (dose-response bound accepted:
−0.42% CT per +48% deletion). Notebook entry for 032 already written
(`journals/20260901.md` § `# 20260919`); GIF done.

## TASK 1 — harvest the two follow-up arms when they finish

Judge by outputs, not sacct; DomainError tracebacks land in the `.err`
(`~/campaigns/p032-rootomit-20260918/FLOWPanel.jl/logs/slurm/slurm-fp-052-scr-gpu-<jobid>.{out,err}`).
Delegate status/tailing to `hpc-monitor`. Per arm:

- Outcome vs acceptance (in the 20260919 provenance): C-ctrl —
  completion+converged CT ⇒ root shed also drives the ctrl channel;
  in-window blow-up ~253–290 ⇒ ctrl channel root-shed-independent
  (A3 blew up physically while sacct said COMPLETED — check CT_per_rev,
  not just exit). C-oldlaw — CT within noise of B15 (0.07162 ± 0.038%)
  ⇒ merge redesign unnecessary under omission (prediction on record:
  near-null); death or CT shift ⇒ redesign contributes independently.
- Harvest numbers: CT cycle-mean ± plateau, GATE line, omission totals
  (deleted Σ|Γ|·Δl — expect ≈0.486 for exp fs; ctrl case may differ),
  final np, merge-event count, wake-health extreme tail (min σ, max
  Γ/σ²). Compare vs the merge-burden table in the 20260918 provenance.
- Run dirs land worktree-local: MOVE to
  `orc:~/projects/FLOWPanel.jl/data/` + symlink back at the worktree
  path (launcher `rm -rf`s pre-made symlinks — 026 precedent).
- Record outcomes in `032_followup_provenance_20260919.md` (Outcomes
  section + verdict paragraph); offer Ryan a 1–2-line notebook
  follow-up appended under the existing `# 20260919` 032 section
  (approval required before writing).

## TASK 2 — BRAINSTORM sweep: what should omission re-open?

Ryan directive: now that root-shed omission is identified as the s9
stability lever, sweep `BRAINSTORM/` and recommend which items should
be **continued or re-opened with omission active**. Method:

- Start from `BRAINSTORM/INDEX.md` for the item list. Fan out cheap
  subagents — `brainstorm-scout` (repo agent, built for this) or
  `Explore` on the CHEAPEST reliable model (haiku) — one per candidate
  item; NEVER read full item files inline (many are 1000+ lines). Ask
  each scout: item goal, why it stalled/parked/closed, whether the
  blocker is (a) an exp-family/fountain blow-up in the ~253–328 band,
  (b) the 500k particle cap, (c) the count-driven FMM σ-adequacy gate,
  or (d) unrelated; and what one omission-enabled run would test.
- Obvious priors to check (from the 032/026 record — verify, don't
  assume): 018 NT72 ladder rungs (all died at the count-driven
  adequacy gate; omission lowers np ~30%); 026 NT144 cap-ladder rungs
  (parked; same gate + Ryan declined un-park 09-19 — note it was
  declined, only re-offer with new justification); 022 ground-effect
  arms (below-ground Γ-ignition — different geometry, root fountain
  may or may not be the driver; IGE fountain flow is STRONGER, check
  the autopsies); 020 σ-closure Phase 3 (parked on Ryan; wake died
  physically); 005 intrinsic hover-wake oscillation; 026 remaining
  discriminators (feature-A move-vs-delete — declined 09-19). Also
  scan for anything citing euler_exp DomainError, overflow @500k, or
  "adequacy" as its death mode.
- Synthesize a RANKED recommendation table (item | blocker class |
  what omission changes | proposed run + est. cost | information
  value) and put it to Ryan via AskUserQuestion (multiSelect). LAUNCH
  NOTHING without approval. Remember exact-rate NT ladder rule
  r(NT)=1−(1−0.3)^(36/NT) for any NT rung; backend-match any A/B;
  campaign worktree rules apply (p032 pins reusable for new-code arms).

## Owed / parked (carried — do not silently drop)

- hpc-storage archiving: 026 slate dirs (40.3 G) + 032 dirs (45 G) +
  the two new arm dirs once harvested; local staging dirs
  `~/scr_p026s9r2_explg_fs_last50steps/`,
  `~/scr_p032om15_explg_fs_steps278_327/` (2.3 G) now GIF-done —
  archivable pending Ryan's OK; `scr_p026gpuv_split` retry once quiet
  ≥24 h; `data/scratch_p032_smoke{A,B,C,D}` disposable.
- Uncommitted by design (Ryan-gated commit bundle): 026 ledger line,
  all provenance-file edit sets (incl. `032_followup_provenance_20260919.md`
  and this file), 032 item-file Log updates.
- github pushes of branches+tags; orc branch-divergence ruling.
- k=3 cap retune + merged-σ/clamp-engagement telemetry gap (A4).
- 021 silo cleanup; `031_quadrupole_panel_farfield.md` missing INDEX row.
- Notebook checklist item under `# 20260919` remains untickable (Ryan
  ticks on approval).

## Ground rules

Local ≤4 threads. `ssh orc` needs a live ControlMaster socket (`! ssh
orc echo ok` if 2FA blocks); `bash -lc` for Slurm; MOTD + ANSI escapes
contaminate ssh output (strip before parsing). All submissions, all
commits, and all notebook writes are Ryan-gated. Pre-existing dirty
files (018/026 docs not from this arc, rotor_multi slurm script, data
TOML) are NOT yours. GIF render script + frames precedent:
`~/Dropbox/research/notebooks/img/20260919_p032_root_omission/render_sidebyside.py`
(pvpython at `/Applications/ParaView-5.12.0.app/Contents/bin/pvpython`).
