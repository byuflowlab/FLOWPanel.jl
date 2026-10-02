# Reset prompt — 018 NT-ladder campaign: L18/L36 harvested PASS, 4 arms pending, archiver detached (2026-09-24 ~08:00)

You are picking up in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl`
(branch **fastmultipole**). Read `CLAUDE.md` + policies. Authoritative
campaign doc (pins, arm matrix, acceptance, job IDs, round-1 post-mortem,
**Results section now holds the L18/L36 harvest**):
`BRAINSTORM/018_dji9443_hover_convergence_campaign/ntladder_provenance_20260923.md`
— read it FIRST. Supersedes `ntladder_reset_prompt_20260923.md`.

## State at reset (verified ~06:30–08:00 09-24; verify, don't re-derive)

1. **L18 13879081 + L36 13879082 (both m13h winners): COMPLETED, acceptance
   PASS, harvested** — full numbers + banner verification in provenance
   §Results (L18 CT 0.068488±5.65e-5; L36 0.070457±6.86e-5; mon04 bounded;
   gate_rc=0 both). Eng twins watcher-cancelled.
2. **NT144-eng 13878882: Ryan ruled KEEP (2026-09-24)** despite running
   first. Was RUNNING healthy at 68% (step 2943/4320, ~14.7 s/step,
   ETA ~late 09-24). Banner PASS.
3. **G-hi 13879087 (m13h winner): RUNNING**, banner PASS (floor 0.5,
   ceil Inf), step 1335/2160 at ~06:30, ~12.2 s/step. Eng twin cancelled.
4. **L72 13879083/92, ANCH 13879084/91, G-lo2 13879085/90, G-lo1
   13879086/89: both twins still PENDING** (Priority/Resources). Pair
   watcher `ntl2_pair_watcher.sh` alive on orc login (pgrep -af, launch
   PID 2570978); loser cancellation remains the only pre-authorized cancel.
5. **Storage: Ryan approved archiving the 3 finished 032 round-2 runs** in
   `/home/rander39/campaigns/p032-rootomit-20260918/FLOWPanel.jl`
   (`..._omi1_split_mo4_eng`, `..._om15_eng`, `..._omi1_split_eng`,
   ~284 GB would-free). Archiver launched DETACHED on orc ~07:55 09-24
   (nohup, PID 402849 at launch), log
   `/home/rander39/archiver_apply_p032rootomit_20260924.log`. /home was
   661.5 G vs 400 G cap; expect ~377 G after. VERIFY completion + before/
   after on resume; the old "~31 G scr + ~212 GB p018_csarc_3r" reclaim
   targets are STALE (already reclaimed in the 09-05 pass; do not chase).
   hpc-storage subagent CANNOT run --apply (permission gate) — apply runs
   must be launched inline with Ryan's in-session approval.
6. Previous local monitor/session died at session end. Lessons carried:
   ANSI-strip all orc output; `bash -lc` for slurm; debounce ×2; a
   terminal-state check must test the CURRENT read, not the debounce seed.

## NEXT ACTIONS

1. Re-arm local job monitor on 13878882, 13879083–87, 13879089–92
   (ANSI-safe, debounced ×2).
2. Banner-verify each new winner as it starts (per provenance §Acceptance);
   NT144-eng and G-hi already verified.
3. Harvest arms as they finish (delegate to hpc-monitor; ask for CT
   cycle-mean, mon04 max/end, particle peak/final, min_sigma,
   floor_clamp_cum, gate_rc, banner NT/N/rlxf/pps) and append to
   provenance §Results.
4. Verify archiver finished + freed space; relay ledger line to Ryan.
5. Key reads once L72/ANCH/NT144 land: ANCH−A5 (0.07136±5.6e-5) =
   linegauss effect; L72−ANCH = N-effect; guard arms vs ANCH = floor dose
   (incl. radial circulation overlay); CT slope NT36→72→144 (current:
   L18 0.0685 → L36 0.0705, rising toward A5).

## BLOCKED on Ryan (carried)

1. Task-1 orc cleanup one-liner (hung PID 1628006, /tmp vpm worktree;
   command in `BRAINSTORM/032_reset_prompt_20260922b.md` §BLOCKED).
2. Docs commit bundle + ALL pushes: 032 round-2 provenance addendum,
   ntladder provenance (incl. new §Results) + reset prompts, 018/026/032
   item Logs, INDEX rows. Notebook: offer NT-ladder entry when the ladder
   completes.
3. Round-1 partial run dirs `data/*_lg_eng/` in the ntladder wt (
   merge_events.csv only) — clean only with Ryan's OK.

## Owed / parked (carried, Ryan-gated)

P3 clean 12 h rerun; `simulate_warmstart!` monitor-clobber fix; k=3 cap
retune + merged-σ telemetry; 021 silo cleanup; `scr_p026gpuv_split` retry
once quiet ≥24 h. 021 campaign is separate and NOT yours (its jobs
13879622/13879625 must not be touched).

## Ground rules (carried verbatim)

Local ≤4 threads; macOS no `timeout`, bash 3.2, shell zsh. `ssh orc`
needs live ControlMaster (else ask Ryan `! ssh orc echo ok`). Plain ssh
has no slurm — wrap `bash -lc`; strip MOTD/ANSI. Never `pkill -f` your
own pattern. Judge runs by outputs, never sacct. Delegate monitoring/
harvest to repo subagents; conclusions inline. Submissions/commits/remote
git/notebook Ryan-gated (standing exception: cancelling a losing ntl2
twin). Never run threelevel SFS. orc julia = spack `julia/1.11.7-6bmogfl`.
VTP coords: rotor axis X, radial sqrt(y²+z²). Pre-existing dirty files
expected — commit only with approval.
