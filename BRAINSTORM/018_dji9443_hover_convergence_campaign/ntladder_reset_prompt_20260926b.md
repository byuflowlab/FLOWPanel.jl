# Reset prompt — 018 NT-ladder r3: S18/S36 PASS+harvested, S144 + N8-144 in final quarter (~ETA 03:00–04:00 on 09-27), harvest next (2026-09-26 late PM)

You are picking up in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl`
(branch **fastmultipole**). Read `CLAUDE.md` + policies. Authoritative
campaign doc (pins, r1/r2/r3 arm matrices, acceptance, job IDs, all
harvested §Results and §Round 3 results rows, r3 rationale):
`BRAINSTORM/018_dji9443_hover_convergence_campaign/ntladder_provenance_20260923.md`
— read it FIRST. Supersedes `ntladder_reset_prompt_20260926.md`.

## State at reset (verified 22:20 on 09-26; verify, don't re-derive)

1. **r1/r2 TERMINAL, 8/8 PASS, harvested, notebook entry written**
   (`journals/20260901.md` under `# 20260925`, checkbox unticked; 3 TikZ
   figs in `notebooks/img/20260925_018_ntladder/`). CT-vs-NT (frozen SFS):
   0.0685 / 0.0705 / 0.0707 / 0.0733.
2. **r3 (SFS-rlxf correction + N=8 rung), worktree 6e9a815, tag
   `campaign/p018-ntladder-r3-20260925` (NOT pushed):**
   - **S18 13898944 m13h: COMPLETE, PASS, HARVESTED** (provenance §Round 3
     results). **S18−L18 = +1.67%** (0.069629 vs 0.068488).
   - **S36 13898946 m13h: COMPLETE, PASS, HARVESTED.** **S36−L36 =
     +0.51%**; corrected coarse-end ladder 0.069629 → 0.070816 → 0.070656
     (L72), i.e. +1.70% / −0.23%. Monotone dose–response confirmed.
   - **S144 13898948 m13h: RUNNING**, step 3298/4319 at 22:20, ~19.7
     s/step → ETA ~04:00 on 09-27. Node m13h-1-2. Nominal (clean forces,
     sigma growth ~1.15, no NaNs/stalls). Run dir
     `<wt>/data/p018_csarc_n4_nt144_l3p0_3r_srlx_g25nc_omi1_split_mo4_lg_sfsr3_m13h/`.
   - **N8-144 13898951 eng: RUNNING**, step 3295/4319 at 22:20, ~15.6
     s/step → ETA ~02:50 on 09-27. Node eng-1-1. Nominal. Run dir
     `.../p018_csarc_n8_nt144_..._sfsr3_eng/`.
   - Both banners already verified vs the r3 matrix incl. ceil OFF
     (SIGMA_CEIL=Inf) and SFS rlxf 0.0012524. Loser twins
     (13898945/47/49/50) cancelled by watcher; watcher exited "all pairs
     resolved" — **no watcher needed anymore**.
3. Slurm logs:
   `/home/rander39/campaigns/p018-ntladder-20260923/FLOWPanel.jl/logs/slurm/slurm-fp-018gpu-ntl3-{s144m,n8144e}-139889{48,51}.{out,err}`.
   Harvest conventions: CT = monitor02 −CFx cycle mean; mon04 = col 7
   `max_gamma_over_sigma2` of monitor04_wake_health; same columns as
   provenance results tables. **A prior harvester misreported deltas
   twice — recompute all percent deltas yourself.**
4. **Storage**: /home was 622 G after the 09-25/26 pass; r2 G-lo2/G-lo1
   VTK (~230 G) aged, r3 VTK growing. Next apply pass Ryan-gated.
5. orc ssh intermittently slow; direct short `ssh orc 'bash -lc ...'`
   works. Plain ssh has no slurm — wrap `bash -lc`.

## NEXT ACTIONS

1. **Check S144 (13898948) + N8-144 (13898951)** via hpc-monitor — they
   are probably COMPLETE by the time you read this. Judge by outputs
   (gate line at log tail: gpu_gemv count, 0 NaN, rc=0), never sacct.
2. **Harvest each arm**: provenance §Round 3 results row + auxiliary
   reads (p2p, mon04 decay, wall, particles, clamp counts), exactly as
   the S18/S36 rows. Recompute all percent deltas yourself.
3. **Campaign-deciding reads once both land:**
   - **Test the standing prediction S144 < 0.0733** (provenance §Round 3
     results, dose–response paragraph): frozen T at NT144 was 2× too
     SHORT → correction should LOWER CT there (opposite sign vs coarse
     rungs). State whether the mechanism is confirmed.
   - S144−NT144(r2) = SFS-timescale effect at NT144.
   - Corrected ladder slope S18→S36→L72→S144 (L72 = 0.070656 reusable).
   - N8-144−S144 = N-effect at NT144.
4. Then propose notebook additions (Ryan-gated; extend tables/figures,
   add sfsr3 series to `notebooks/img/20260925_018_ntladder/` dirs) and
   a storage pass (Ryan-gated apply).

## BLOCKED on Ryan (carried)

1. Task-1 orc cleanup one-liner (hung PID 1628006, /tmp vpm worktree;
   command in `BRAINSTORM/032_reset_prompt_20260922b.md` §BLOCKED).
2. Docs commit bundle + ALL pushes: ntladder provenance (incl. r3
   S18/S36 rows), reset prompts, 018/026/032 item Logs, INDEX rows,
   campaign tag pushes (`campaign/p018-ntladder-r3-20260925`).
3. Round-1 partial run dirs `data/*_lg_eng/` in the ntladder wt
   (merge_events.csv only) — clean only with Ryan's OK.
4. Notebook r3 additions (entry edits are Ryan-gated).

## Owed / parked (carried, Ryan-gated)

P3 clean 12 h rerun; `simulate_warmstart!` monitor-clobber fix; k=3 cap
retune + merged-σ telemetry; 021 silo cleanup; `scr_p026gpuv_split` retry
once quiet ≥24 h. 021 campaign (job 13899021 on m12) is separate and NOT
yours.

## Ground rules (carried verbatim)

Local ≤4 threads; macOS no `timeout`, bash 3.2, shell zsh. `ssh orc` needs
live ControlMaster (else ask Ryan `! ssh orc echo ok`). Plain ssh has no
slurm — wrap `bash -lc`; strip MOTD/ANSI. Never `pkill -f` your own
pattern. Judge runs by outputs, never sacct. Delegate monitoring/harvest
to repo subagents; conclusions inline (and re-check subagent arithmetic).
Submissions/commits/remote git/notebook Ryan-gated (standing exception:
cancelling a losing ntl3 twin was pre-authorized; all pairs now resolved).
Never run threelevel SFS. orc julia = spack `julia/1.11.7-6bmogfl`. VTP
coords: rotor axis X, radial sqrt(y²+z²). Pre-existing dirty files
expected — commit only with approval. Notebook: today's-entry edits OK,
previous days append-only.
