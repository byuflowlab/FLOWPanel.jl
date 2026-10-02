# Reset prompt — 018 NT-ladder: 6/8 arms harvested PASS, G-lo2 running, ANCH-eng pending (2026-09-24 ~PM)

You are picking up in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl`
(branch **fastmultipole**). Read `CLAUDE.md` + policies. Authoritative
campaign doc (pins, arm matrix, acceptance, job IDs, round-1 post-mortem,
**§Results now holds 6 harvested arms + harvest-2 notes**):
`BRAINSTORM/018_dji9443_hover_convergence_campaign/ntladder_provenance_20260923.md`
— read it FIRST. Supersedes `ntladder_reset_prompt_20260924.md`.

## State at reset (verified ~PM 09-24; verify, don't re-derive)

1. **Harvested acceptance-PASS (6 arms, rows in provenance §Results):**
   L18 0.068488, L36 0.070457, L72 0.070656 (m13h 13879083),
   G-hi 0.0714012 (m13h 13879087), G-lo1 0.0713961 (eng 13879089),
   NT144 0.073327 (eng 13878882). All gate_rc=0, verbatim banners
   verified (incl. floors 0.25/0.5/0.125/0.25 — an earlier monitor pass
   mis-parsed sigma_floor=0; refuted). mon04 max==end caveat noted in
   provenance (values all « A5 ref 154.8).
2. **Early reads (details in provenance):** ladder CT rising, NT72→144
   +3.8% but ceil-confounded (NT144 ceil off, L-rungs 0.030); floor
   dose-response ~null at NT72/N=2 (G-hi vs G-lo1 Δ0.007% despite 6300×
   clamp difference).
3. **Live:** G-lo2 eng 13879090 RUNNING (banner PASS, floor 0.0625/Inf);
   ANCH eng 13879091 PENDING (Resources). **Ryan manually killed ANCH
   m13h 13879084 (2026-09-24, acknowledged fine)** — record ANCH winner
   as "eng (m13h twin killed manually by Ryan)". Pair watcher
   `ntl2_pair_watcher.sh` alive on orc login (pgrep -af; launch PID
   2570978); loser cancellation remains the only pre-authorized cancel.
4. **Archiver (032 round-2): COMPLETED** — 283.9 GB freed, 0 verify
   fails, DF_AFTER=467G. Still over the 400 G cap (~90 G ≈ live ladder
   growth). Another storage pass needed once the ladder finishes;
   hpc-storage subagent cannot --apply (permission gate) — apply runs
   launch inline with Ryan's in-session approval only.
5. Previous session's local monitor died with the session — re-arm.

## NEXT ACTIONS

1. Re-arm local job monitor on 13879090 + 13879091 (poll squeue via
   `ssh orc 'bash -lc ...'`, ANSI-strip, debounce ×2 on the CURRENT
   read; script pattern: scratchpad `ntl2_local_monitor.sh` from prior
   session, or rewrite).
2. Banner-verify ANCH 13879091 when it starts (expect N=2 NT=72
   rlxf=0.16334 pps=6 floor 0.25/ceil 0.030, linegauss, omission 1/41
   @0.12R, split+stretch, ct4, SFS=DynamicSFS(rlxf=0.0025031, maxC=1.0,
   alpha=0.999, clippings=backscatter, controls=none, nostatic=false)).
3. Harvest G-lo2 + ANCH as they finish (delegate hpc-monitor; ask for
   CT cycle-mean, mon04 max/end, particle peak/final, min_sigma,
   floor_clamp_cum, gate_rc, verbatim banner) → append provenance
   §Results rows matching existing columns.
4. Once ANCH lands, compute the key reads: ANCH−A5 (0.07136±5.6e-5) =
   linegauss effect; L72−ANCH = N-effect (2→4) at NT72; guard arms
   (G-lo2/G-lo1/G-hi/ANCH) = floor dose-response incl. radial
   circulation overlay; full CT-vs-NT slope with the ceil confound
   flagged. Offer Ryan the notebook entry (Ryan-gated).
5. After ladder terminal: propose next storage pass (Ryan-gated apply).

## BLOCKED on Ryan (carried)

1. Task-1 orc cleanup one-liner (hung PID 1628006, /tmp vpm worktree;
   command in `BRAINSTORM/032_reset_prompt_20260922b.md` §BLOCKED).
2. Docs commit bundle + ALL pushes: 032 round-2 provenance addendum,
   ntladder provenance (incl. §Results harvests 1+2) + reset prompts,
   018/026/032 item Logs, INDEX rows. Notebook: offer NT-ladder entry
   when the ladder completes.
3. Round-1 partial run dirs `data/*_lg_eng/` in the ntladder wt
   (merge_events.csv only) — clean only with Ryan's OK.

## Owed / parked (carried, Ryan-gated)

P3 clean 12 h rerun; `simulate_warmstart!` monitor-clobber fix; k=3 cap
retune + merged-σ telemetry; 021 silo cleanup; `scr_p026gpuv_split`
retry once quiet ≥24 h. 021 campaign is separate and NOT yours (jobs
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
