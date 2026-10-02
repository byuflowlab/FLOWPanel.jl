# Reset prompt — 018 NT-ladder r3 COMPLETE: 4/4 PASS harvested, prediction confirmed, ladder NOT converged; everything remaining is Ryan-gated (2026-09-28)

You are picking up in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl`
(branch **fastmultipole**). Read `CLAUDE.md` + policies. Authoritative
campaign doc (pins, all arm matrices, harvested §Results and §Round 3
results incl. the campaign-deciding reads):
`BRAINSTORM/018_dji9443_hover_convergence_campaign/ntladder_provenance_20260923.md`
— read it FIRST. Supersedes `ntladder_reset_prompt_20260926b.md`.

## State at reset (verified 2026-09-28)

1. **r1/r2 TERMINAL 8/8 PASS** (notebook entry `journals/20260901.md`
   `# 20260925`, checkbox unticked; 3 TikZ figs in
   `notebooks/img/20260925_018_ntladder/`).
2. **r3 TERMINAL 4/4 PASS, fully harvested + settings-audited**
   (worktree 6e9a815, tag `campaign/p018-ntladder-r3-20260925`, NOT
   pushed). CT (window convention, revs 20–30): S18 0.069629, S36
   0.070816, S144 0.072402, N8-144 0.071396. Key reads (provenance
   §Round 3 results, recomputed inline):
   - **Prediction CONFIRMED**: S144 = 0.072402 < 0.0733;
     S144−NT144(r2) = −1.26% (sign flip vs +1.67%/+0.51% coarse rungs —
     SFS-timescale mechanism confirmed; magnitude asymmetric).
   - **TWO ladders (framing per Ryan 2026-09-29)**: constant-handoff
     (N∝NT, design intent) S18→S36→L72→**N8-144** = +1.70% / −0.23% /
     **+1.05%** (endpoint 0.071396); fixed-N pairs S36→ANCH (N2) =
     +0.84%, L72→S144 (N4) = **+2.47%** (endpoint 0.072402). NEITHER
     closes; policies diverge (endpoints differ by the N-effect). Never
     quote the mixed S18/S36/L72/S144 sequence as "the" ladder.
   - **N-effect**: −1.06% @NT72, −1.39% @NT144 per doubling — grows
     with refinement, so N axis not converged either (cross-partition
     caveat: N8 eng vs S144 m13h).
   - **Settings audit 2026-09-28** (banner-diff L36/L72/NT144r2/S-arms):
     only intended diffs + N covary + σ-ceil off @NT144 (moot: run-max
     σ 0.0174 < 0.030 everywhere) + MAX_PARTICLES (never binds).
     Non-convergence is NOT a settings artifact.
   - **CT harvest convention**: use `case_metadata.toml` `CT_window_mean`
     ± `CT_cycle_std` (window p2p = `CT_ptp_rel`) — NEVER hand-window
     monitor02 (final-rev-only gave 0.072267, wrong). Harvester
     subagents also botched CT arithmetic a third time.
3. No jobs running; watcher exited (all pairs resolved). 021 job
   13899021 is separate and NOT yours.
4. **Storage**: /home was 622 G before r3 VTK growth (S144 particles
   peaked 1.07M, VTK every step both NT144 arms) — next hpc-storage
   apply pass is needed but Ryan-gated.

## NEXT ACTIONS (all Ryan-gated — present, don't do)

1. Propose notebook additions for r3 (extend `# 20260925` tables/figures,
   add sfsr3 series to `notebooks/img/20260925_018_ntladder/` CSV dirs;
   ask Ryan how much to write).
2. Propose a possible r4 (NT288 rung?) — ladder unconverged; Ryan
   decides.
3. Storage apply pass (hpc-storage), commit bundle + tag pushes.

## BLOCKED on Ryan (carried)

1. Task-1 orc cleanup one-liner (hung PID 1628006; command in
   `BRAINSTORM/032_reset_prompt_20260922b.md` §BLOCKED).
2. Docs commit bundle + ALL pushes: ntladder provenance (now incl. all 4
   r3 rows + reads), reset prompts, 018/026/032 item Logs, INDEX rows,
   tag `campaign/p018-ntladder-r3-20260925`.
3. Round-1 partial run dirs `data/*_lg_eng/` in the ntladder wt.
4. Notebook r3 additions.

## Owed / parked (carried, Ryan-gated)

P3 clean 12 h rerun; `simulate_warmstart!` monitor-clobber fix; k=3 cap
retune + merged-σ telemetry; 021 silo cleanup; `scr_p026gpuv_split` retry
once quiet ≥24 h.

## Ground rules (carried verbatim)

Local ≤4 threads; macOS no `timeout`, bash 3.2, shell zsh. `ssh orc` needs
live ControlMaster (else ask Ryan `! ssh orc echo ok`). Plain ssh has no
slurm — wrap `bash -lc`; strip MOTD/ANSI (nested quoting breaks — pipe a
script via `ssh orc 'bash -s' < file` instead). Never `pkill -f` your own
pattern. Judge runs by outputs, never sacct. Delegate monitoring/harvest
to repo subagents; conclusions inline (and re-check subagent arithmetic).
Submissions/commits/remote git/notebook Ryan-gated. Never run threelevel
SFS. orc julia = spack `julia/1.11.7-6bmogfl`. VTP coords: rotor axis X,
radial sqrt(y²+z²). Pre-existing dirty files expected — commit only with
approval. Notebook: today's-entry edits OK, previous days append-only.
