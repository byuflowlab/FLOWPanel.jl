# Reset prompt — 018 NT-ladder campaign: r2 twins queued, NT144-eng running (2026-09-23 ~23:00)

You are picking up work in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl`
(branch **fastmultipole**). Read `CLAUDE.md` + the policies it names.
Authoritative campaign doc (pins, arm matrix, acceptance, job IDs, failure
post-mortem): `BRAINSTORM/018_dji9443_hover_convergence_campaign/ntladder_provenance_20260923.md`
— read it FIRST; do not re-derive anything in it.

## Context in one paragraph

032 P2-rerun round 2 closed 2026-09-23: A4/A5/A6 all passed; champion =
A5 config (0.12R omission 1/41, split+stretch, MERGE_OVERLAP=4, DynamicSFS
rlxf=0.0025031 backscatter, σ guard floor 0.25→0.00119 m). Ryan then approved
this NT-ladder campaign: switch gaussian→**linegauss**, constant-handoff
ladder at 20° (NT18/N1 exploratory, NT36/N2, NT72/N4), **ANCH** = NT72/N2
champion-geometry anchor (linegauss A/B vs A5 + N-effect link + guard-sweep
base), σ-guard floor sweep at ceiling-off (0.0625/0.125/0.5 vs ANCH=0.25),
and NT144/N4 fine rung (floor-only 0.25, MAX_PARTICLES=3M, 48 h). All from
pinned wt `orc:/home/rander39/campaigns/p018-ntladder-20260923/FLOWPanel.jl`
(tag `campaign/p018-ntladder-20260923`=`3d7fba4`=`af92740`+case rows), env
`.../p018-ntladder-20260923/env`. Open physics question stands: family CT
~0.0714 ≈ experiment vs historical 0.0506 — unexplained; NT144-only rung and
L18 exploratory caveats in provenance.

## State at reset (established ~23:00 09-23 — verify, don't re-derive)

1. **Round-1 submission FAILED wholesale**: untracked asset
   `data/p018_cs_l3p4_rs1_te_downwash_te.csv` missing from the fresh
   worktree → all 7 non-NT144 arms died at init on eng; watcher had
   cancelled the m13h twins. FIXED: CSV copied+md5-verified into `<wt>/data/`.
2. **NT144-eng 13878881→13878882 pair**: 13878882 RUNNING since 22:41:27,
   healthy banner (linegauss, N4, NT144, pps 3), CSV race won. It runs
   FIRST, violating Ryan's "NT144 last" — **Ryan-gated: keep or kill**;
   nobody cancels it without his word. n144m 13878881 watcher-cancelled.
3. **r2 resubmission of all 7 arms** (~22:58, job IDs in provenance table:
   13879081–87 m13h, 13879088–94 eng; run names `..._lg_r2_{m13h,eng}`,
   job names `fp-018gpu-ntl2-*`). Not yet banner-verified — none had
   started at reset.
4. Pair watcher `ntl2_pair_watcher.sh` alive on orc login (launch PID
   2570978; check `pgrep -af ntl2_pair_watcher`), log
   `<wt>/logs/slurm/ntl2_pair_watcher.log`. Loser cancellation is the ONLY
   pre-authorized cancel. Old `ntl_pair_watcher` should have exited.
5. The previous session's local Monitor died with it. Its lessons: strip
   ANSI (`sed $'s/\x1b\[[0-9;]*m//g'`) from ALL orc output before parsing
   (MOTD escape codes fuse onto squeue's first line — cost us 3 false
   alarms), zsh doesn't word-split unquoted vars, macOS bash is 3.2.

## NEXT ACTIONS (in order)

1. Re-arm a local job-state monitor on 13879081–94 + 13878882 (ANSI-safe,
   debounced ×2; template lesson above).
2. Banner-verify each r2 WINNER as it starts (logs
   `<wt>/logs/slurm/slurm-fp-018gpu-ntl2-<arm><m|e>-<jobid>.out`):
   linegauss ("LineGaussRegularization (pinned by ...)"), correct
   NT/NWAKEROWS/rlxf/pps per arm (provenance matrix), guard floor/ceil per
   arm (sweep arms have NO SIGMA_CEIL → prints Inf, guard still on via
   floor), omission 1/41 @0.12R, split+stretch active, mesh 45_185_ct4,
   H200, 64 threads, pinned paths, no ERRORs (32 benign FastMultipole
   constant-redefinition warnings in .err are expected).
3. Babysit vs acceptance (provenance §Acceptance): complete steps
   (540/1080/2160/4320 by NT), finite CT, bounded monitor04 Γ/σ², gate_rc=0.
   Judge by outputs, never sacct. L18 pre-declared exploratory (may fail
   without contaminating the trend).
4. Harvest as arms finish: CT cycle-mean (monitor02 −CFx IS CT, verified),
   monitor04 max/end, particle peak/final, min_sigma/floor_clamp_cum,
   record twin winners. Key reads: ANCH−A5 (0.07136±5.6e-5) = linegauss
   effect; L72−ANCH = N-effect; guard arms vs ANCH = floor dose (CT +
   radial circulation distribution overlay); CT slope NT36→72(→144).
5. Report the NT144-first question to Ryan if unresolved (item 2 above).
6. Failed round-1 partial run dirs `data/*_lg_eng/` (merge_events.csv only)
   — ignore; clean only with Ryan's OK.

## BLOCKED on Ryan (carried — do NOT launder)

1. NT144-eng keep-or-kill (above).
2. hpc-storage re-dispatch (carried from 032 prompts: ~31 G scr archives +
   212 GB old `p018_csarc_*_3r_*` reclaim — EXCLUDE `om15/omi1/srlx_g25`
   AND `_m13h/_eng/_lg` run dirs — live/new runs match the old glob!).
3. Task-1 orc cleanup one-liner (hung PID 1628006, /tmp vpm worktree; command
   in `BRAINSTORM/032_reset_prompt_20260922b.md` §BLOCKED).
4. Docs commit bundle + ALL pushes: 032 round-2 provenance addendum (see
   `032_reset_prompt_20260923.md` §5), THIS campaign's new files
   (ntladder provenance + this prompt), 018/026/032 item Logs, INDEX rows.
   Notebook: 20260923 entry EXISTS (032 ladder summary, checkbox untucked);
   offer an NT-ladder-campaign entry only when results land.

## Owed / parked (carried, Ryan-gated)

P3 clean 12 h rerun (ask before submit); `simulate_warmstart!`
monitor-clobber fix; k=3 cap retune + merged-σ telemetry; 021 silo cleanup;
`scr_p026gpuv_split` retry once quiet ≥24 h. 021 campaign is separate and
NOT yours (entry `fgs_scalability_stage2_reset_prompt_20260923.md`; its jobs
must not be touched).

## Ground rules

Local ≤4 threads; macOS has NO `timeout`, bash is 3.2, default shell zsh.
`ssh orc` needs a live ControlMaster socket — if it hangs ask Ryan to run
`! ssh orc echo ok` (2FA); NEVER retry into a 2FA loop. Plain
`ssh orc '<cmd>'` has NO slurm in PATH — wrap in `bash -lc` (a resubmission
was lost to this at 22:55); strip MOTD/ANSI from every parsed output. Never
`pkill -f` a pattern in your own command line. Judge runs by outputs, never
sacct. Delegate monitoring/harvest/scouting to repo subagents; conclusions
inline. New submissions/commits/remote git/notebook writes Ryan-gated
(standing exception: cancelling a losing twin of the 7 ntl2 pairs). Never
run threelevel SFS. orc login julia is 1.12 — campaign env ops use spack
`julia/1.11.7-6bmogfl`. VTP coords: rotor axis X, radial sqrt(y²+z²).
Pre-existing dirty files (018/026/032 docs, rotor_multi slurm script,
pressure-comparison TOML) expected — commit only with approval.
