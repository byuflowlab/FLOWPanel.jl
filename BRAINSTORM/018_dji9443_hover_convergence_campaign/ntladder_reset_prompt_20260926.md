# Reset prompt — 018 NT-ladder: r1/r2 COMPLETE (8/8 PASS, notebook+figures done), r3 SUBMITTED (SFS-rlxf fix + N=8) (2026-09-26)

You are picking up in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl`
(branch **fastmultipole**). Read `CLAUDE.md` + policies. Authoritative
campaign doc (pins, arm matrices r1/r2/**r3**, acceptance, job IDs, all 8
harvested §Results rows, key reads, r3 rationale):
`BRAINSTORM/018_dji9443_hover_convergence_campaign/ntladder_provenance_20260923.md`
— read it FIRST. Supersedes `ntladder_reset_prompt_20260924b.md`.

## State at reset (verified ~AM 09-26; verify, don't re-derive)

1. **r1/r2 ladder TERMINAL, 8/8 arms PASS and harvested** (provenance
   §Results). Key reads: linegauss NULL (ANCH−A5 +0.07%, 0.5σ); floor AND
   ceil dose-response flat at NT72/N=2 (spread 0.02% over clamp 0→74.3M);
   pure N-effect (2→4, equal pps 6) −1.06%; CT-vs-NT 0.0685/0.0705/0.0707/
   0.0733 (+2.9/+0.28/+3.8%). mon04 = col 7 `max_gamma_over_sigma2` of
   monitor04_wake_health CSV; G-lo2 peaked 164.1 > A5 ref 154.8 but decayed
   (no ignition). NT144/G-hi/ANCH missed the informational Phase 2e p2p
   criterion (within-rev ripple >2%).
2. **Notebook entry WRITTEN** (`~/Dropbox/research/notebooks/journals/
   20260901.md` under `# 20260925`): 4 ladder tables (incl. shedding/Das/
   rlxf/SFS-rlxf columns + pps clarification), p2p-criterion note, 3 TikZ
   figures (Γ(r/R), KJ spanwise loading, CT-vs-rev) in
   `notebooks/img/20260925_018_ntladder/` (compiled, KJ integrals verified
   0.5–0.7% vs measured CT). Checkbox unticked (Ryan approves).
   Circulation read: NT144's CT excess is distributed mid-span
   (0.15–0.8 r/R), NOT tip; NT36/NT72 collapse.
3. **SFS-rlxf confound found (Ryan 2026-09-25)**: DynamicSFS rlxf = Δt/T;
   r1/r2 froze 0.0025031 (= exact-rate NT72 image of default 0.005@NT36)
   → averaging time T varied 4× down the ladder. L72 already has the
   correct NT72 value (reusable). Full rationale in provenance §Round 3.
4. **r3 SUBMITTED 2026-09-26** (Ryan-approved), 4 arms × m13h/eng twins:
   S18 13898944/45, S36 13898946/47, S144 13898948/49, N8-144 13898950/51.
   SFS_RLXF 0.009975/0.005/0.0012524/0.0012524; N8-144 = NEW case
   `p018_csarc_n8_nt144_l3p0` (N=8, wall 72 h; S144 48 h); NT144 arms
   afterany-held behind short arms. Worktree commit **6e9a815**, tag
   **campaign/p018-ntladder-r3-20260925** (NOT pushed). At reset:
   s18-m13h RUNNING, rest PENDING. Watcher `ntl3_pair_watcher.sh` detached
   on orc login (launch PID 2629745; verify `pgrep -af ntl3_pair_watcher`);
   loser cancellation is the ONLY pre-authorized cancel. Submit script
   `/home/rander39/campaigns/p018-ntladder-20260923/p018_submit_ntl3.sh`.
5. **Storage pass 2026-09-25/26 COMPLETE**: TOTAL_FREED_MB=375,762,
   VERIFY_FAIL_COUNT=0, STALE_COUNT=0; /home 937→622 G. Remaining bulk =
   r2 G-lo2/G-lo1 VTK (~230 G, classified RECENT at pass time, now aged) +
   live r3 growth. Next apply pass is Ryan-gated; note the auto-mode
   classifier blocks the detached `--apply` launch — Ryan flipped to
   manual mode / `dangerouslyDisableSandbox` last time.

## NEXT ACTIONS

1. Babysit r3: verify each arm's banner as it starts (expect linegauss,
   correct N/NT/pps per §Round 3 matrix, **SFS=DynamicSFS(rlxf=<per-arm
   value>, maxC=1.0, alpha=0.999, clippings=backscatter)**, floor 0.25,
   ceil 0.030 for S18/S36 / guard-off for the NT144 arms, omission 1/41
   @0.12R, split+stretch, ct4). Delegate log checks to hpc-monitor.
2. Harvest arms as they finish (same columns as §Results; mon04 = col 7;
   CT from monitor02 −CFx). Reduced-CSV extraction pattern:
   `/home/rander39/ntl_extract.sh` (adapt run names `..._lg_sfsr3_*`).
3. Key r3 reads once landed: S18−L18, S36−L36, S144−NT144(r2) =
   SFS-timescale effect per rung; corrected ladder slope S18→S36→L72→S144;
   N8-144−S144 = N-effect at NT144. Extend the notebook figures/tables
   (Ryan-gated entry addition; figure dirs exist, add sfsr3 series).
4. N8-144 watch: expect ~1.5M particles (cap 3M); if the 72 h wall looks
   tight mid-run, warm-start options exist — ask Ryan.
5. After r3 terminal: propose storage pass (Ryan-gated apply) — r2 RECENT
   runs will have aged and r3 VTK adds ~150+ G.

## BLOCKED on Ryan (carried)

1. Task-1 orc cleanup one-liner (hung PID 1628006, /tmp vpm worktree;
   command in `BRAINSTORM/032_reset_prompt_20260922b.md` §BLOCKED).
2. Docs commit bundle + ALL pushes: ntladder provenance (now incl.
   harvests 1–3 + §Round 3), reset prompts, notebook-adjacent figure
   dirs are in Dropbox (no commit needed), 018/026/032 item Logs, INDEX
   rows, campaign tag push (`campaign/p018-ntladder-r3-20260925`).
3. Round-1 partial run dirs `data/*_lg_eng/` in the ntladder wt
   (merge_events.csv only) — clean only with Ryan's OK.

## Owed / parked (carried, Ryan-gated)

P3 clean 12 h rerun; `simulate_warmstart!` monitor-clobber fix; k=3 cap
retune + merged-σ telemetry; 021 silo cleanup; `scr_p026gpuv_split` retry
once quiet ≥24 h. 021 campaign is separate and NOT yours.

## Ground rules (carried verbatim)

Local ≤4 threads; macOS no `timeout`, bash 3.2, shell zsh. `ssh orc` needs
live ControlMaster (else ask Ryan `! ssh orc echo ok`). Plain ssh has no
slurm — wrap `bash -lc`; strip MOTD/ANSI. Never `pkill -f` your own
pattern. Judge runs by outputs, never sacct. Delegate monitoring/harvest
to repo subagents; conclusions inline. Submissions/commits/remote git/
notebook Ryan-gated (standing exception: cancelling a losing ntl3 twin).
Never run threelevel SFS. orc julia = spack `julia/1.11.7-6bmogfl`. VTP
coords: rotor axis X, radial sqrt(y²+z²). Pre-existing dirty files
expected — commit only with approval. Notebook: today's-entry edits OK,
previous days append-only.
