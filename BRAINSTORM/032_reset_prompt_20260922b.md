# Reset prompt — P2-rerun ROUND 2: A6/A4/A5 twin-queued (m13h+eng), watcher v2 live (2026-09-22 pm)

You are picking up work in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl`
(branch **fastmultipole**) + sibling `/Users/ryan/Dropbox/research/projects/FLOWVPM.jl`
(branch `flowpanel`, `8d4a3b4`). Read `CLAUDE.md` + the policies it names.
Supersedes `032_reset_prompt_20260922.md` (its round-1 arms are all terminal;
outcomes below). Context files, read as needed, do NOT re-derive:
`BRAINSTORM/032_p2rerun_provenance_20260921.md` (arm matrix, pins, acceptance
criteria, submit template — still authoritative for acceptance).

## Round-1 outcomes (2026-09-22, established facts — do NOT re-harvest)

| Arm | Job | Outcome |
|---|---|---|
| A1 ctrl (3nb SFS, om15=0.15R 3/41, guard) | 13842791 eng | **FAIL — true Γ-ignition**: monitor04 bounded ~10⁶ through step 1756, exploded 3.0e7→8.8e8 over 1757–1776, erratic particle counts (156k→23k→138k), died step 1776 at particle overflow. Same P2 ignition delayed ~7 steps — σ guard alone did NOT rescue. |
| A2 base (DynamicSFS rlxf=0.0025031, omi1=0.12R 1/41, guard) | 13843184 m13h | **PASS**: 2160/2160, monitor04 max 451 end 356, CT(−CFx, last 144 steps) 0.0714. |
| A3 (A2+MERGE_OVERLAP=4) | 13843185 m13h | **PASS**: 2160/2160, monitor04 max 513 end 503, CT 0.0713 → merge-gate CT null at NT72. |
| A4/A5 (split arms) | 13843727/28 m13h | **INIT FAIL, zero steps**: driver guard `rotor_hover_pressure_comparison.jl:813` — FRAC_COMPRESS/ELONGATE set without `WAKE_SPLIT_STRETCH=true`. The 026 quartet summary omitted STRETCH (026 launcher case-defs carried it internally). Banners otherwise clean; no run dirs/metadata created. Eng twins 13842794/95 were watcher-cancelled. |

CT CAVEAT: 0.0714/0.0713 ≈ experiment (0.072) but far above historical
panel-converged ~0.0506 — harvester read −CFx from monitor02 raw; verify
normalization before citing. Guard offset +0.39%: compare slopes not levels.

**Ryan rulings (2026-09-22):** (1) confound-resolution arm = **DynamicSFS +
om15** ("A6"); **never run threelevel SFS again**. (2) Resubmit A4/A5 with
`WAKE_SPLIT_STRETCH=true`. (3) Twin all three onto eng in REVERSE order
(a5e first) with **distinct run names** (`_eng` suffix) so double-starts
can't clobber. All three submissions + eng twins DONE this session.

## Round-2 job table (all pinned wt `orc:/home/rander39/campaigns/p032-rootomit-20260918/FLOWPanel.jl`, tag `campaign/p032-rootomit-20260918`=`af92740`, env `.../env`; case `p018_csarc_n2_nt72_l3p0`; guard env SIGMA_FLOOR_FRAC=0.25, SIGMA_CEIL=0.030, TRUNC 3.0R, MAX_PARTICLES=1.5M, SETTLE 22)

| Arm | m13h job | eng twin | Run names (m13h / eng) | Extra env | Wall |
|---|---|---|---|---|---|
| A6 srlx+om15 | 13860852 | 13861315 (a6e) | `..._3r_srlx_g25_om15` / `..._om15_eng` | SFS_RLXF=0.0025031, PARTICLE_OMIT_ROOT_R_OVER_R=0.15 | 16 h |
| A4 retry | 13860853 | 13861314 (a4e) | `..._omi1_split_m13h` / `..._omi1_split_eng` | +WAKE_SPLIT_VISCOUS=true, **WAKE_SPLIT_STRETCH=true**, FRACs 0.587/0.73/0.3 | 24 h |
| A5 retry | 13860854 | 13861313 (a5e) | `..._omi1_split_mo4_m13h` / `..._split_mo4_eng` | A4 env + MERGE_OVERLAP=4 | 24 h |

All 6 PD at reset (~21:05 09-22). m13h estimates (worst-case): a6m 09-23
01:31, a4m2 18:00, a5m2 20:15; eng estimate ~09-23 23:48 (twins are
backfill insurance only). All 40 H200s busy at probe; mgh (2 free GH200)
is ARM — unusable for the x86-pinned env; preempt pools still ruled out
(warmstart clobbers monitor CSVs). Exact submit lines: recover via
`sacct -j <id> --format=SubmitLine%2000 -P`.

**Pair watcher v2** detached on orc login (PID 695377 at launch):
`/home/rander39/campaigns/p032-rootomit-20260918/p2rr_pair_watcher_v2.sh`,
log `<wt>/logs/slurm/p2rr_pair_watcher_v2.log`. Hardcoded to ONLY these 6
IDs; cancels a pair's loser when the other RUNs; both-running tiebreak by
start time; exits when all pairs resolve. Loser cancellation (manual too)
is pre-authorized; nothing else is. Launch pattern if dead:
`ssh -f orc 'bash -lc "setsid nohup <script> >/dev/null 2>&1 </dev/null &"'`.
NOTE: the previous session's local Monitor died with that session — re-arm
your own squeue poller (or delegate to hpc-monitor periodically).

## NEXT ACTIONS (in order)

1. Check job/watcher state (hpc-monitor). Banner-verify each WINNER when
   it starts. Logs `<wt>/logs/slurm/slurm-fp-018gpu-p2rr-a[456]{m2,m,e}-<jobid>.out|.err`
   (job names a6m/a4m2/a5m2/a5e/a4e/a6e). Checklist: guard on (floor_frac
   0.25=0.00119 m, ceil 0.03 m), DynamicSFS rlxf=0.0025031
   clippings=backscatter, mesh 45_185_ct4, NT72, RPM5400, depth 4R, rlxf
   0.16334, H200, 64 threads, pinned paths, no ERRORs. A6: omission 3/41
   at 0.15R. A4/A5: omission 1/41 at 0.12R, WAKE_SPLIT_VISCOUS=true +
   STRETCH=true + FRACs 0.587/0.73/0.3 (must PASS the line-813 guard now);
   MERGE_OVERLAP never prints — check metadata TOML `sigma_relative`
   (A5 true, A4 false). Record which twin won (run-dir name) for harvest.
2. Babysit vs acceptance (032_p2rerun_provenance_20260921.md): survive
   step 1750, complete 2160, finite CT, bounded monitor04. Judge by
   outputs, never sacct. **A6 decision tree**: A6 survives → A1's round-1
   death attributed to threelevel SFS (retired anyway) → om15 viable,
   omission-width lever cleared under champion SFS. A6 ignites → wide
   omission implicated → escalate to Ryan. Split arms: also watch particle
   count vs 1.5M cap + merged-σ behavior.
3. Verify A2/A3 CT normalization (caveat above) before anyone cites 0.0714.
4. **Provenance addendum owed** (Ryan-gated docs): round-1 outcomes table
   above + twin story (round-1 IDs 13843184/85/727/728 + watcher) + the
   WAKE_SPLIT_STRETCH post-mortem + round-2 matrix + A6 arm + Ryan
   rulings + still-owed P3 splice note (13774448→13829223, restart step
   293, monitor history lost).

## BLOCKED on Ryan (do NOT launder via peers)

1. hpc-storage re-dispatch — STILL never run: (a) archive
   `scr_p032om15_ctrllg_fs` + `scr_p032om15_explg_fs_oldlaw` (~31 G);
   (b) 212 GB reclaim of OLD `p018_csarc_*_3r_*` runs — **exclude
   `om15`/`omi1`/`srlx_g25` names AND `_m13h`/`_eng` run dirs** (live runs
   match the old glob!). Reconcile df ~287 G/2.0 T vs 396.4/400 quota.
2. Task-1 orc cleanup one-liner (hung PID 1628006, /tmp vpm worktree):
   `! ssh orc 'kill 1628006; cd ~/projects/FLOWVPM.jl && git worktree remove /tmp/rander39_vpm_mergetest --force; rm -rf /tmp/rander39_vpm_testenv /tmp/rander39_vpm_test_20260919d.log'`

## Owed / parked (carried, Ryan-gated)

- Docs commit bundle: 026 ledger line; provenance edits; 032 item-file Log
  update (reopen verdicts + P3 splice + rerun matrix + round-1 outcomes +
  THIS round-2 session); INDEX.md rows (032 refresh + missing
  `031_quadrupole_panel_farfield.md`); new files incl. this prompt.
- ALL pushes (github + orc, incl. deferred FLOWVPM force-with-lease).
- Notebook entries (×4 from 021 + the 032 reopen/rerun arc) — offer only.
- P3 clean 12 h rerun (citable CT) — pre-authorized in principle, still
  ask before submitting. `simulate_warmstart!` monitor-clobber fix (gated).
- k=3 cap retune + merged-σ/clamp telemetry (026 A4); 021 silo cleanup;
  `scr_p026gpuv_split` retry once quiet ≥24 h; NT144 offer OFF (A1 died;
  revisit only if Ryan re-raises after A6). 021 entry point:
  `fgs_scalability_reset_prompt_20260922.md` (see MEMORY for 021 state).

## Ground rules

Local ≤4 threads; macOS has NO `timeout` (poll-and-kill wrappers).
`ssh orc` needs a live ControlMaster socket — after a reset ask Ryan to
run `! ssh orc echo ok` (2FA) if it hangs; NEVER retry into a 2FA loop.
Strip MOTD/ANSI. Plain `ssh orc '<cmd>'` has NO slurm in PATH — wrap
cluster commands in `bash -lc`. Never `pkill -f` a pattern contained in
your own remote command line. Judge runs by outputs, never sacct.
Delegate monitoring/harvest/scouting to repo subagents (`hpc-monitor`,
`harvester`, `brainstorm-scout`, `code-scout`); keep conclusions inline.
New submissions/commits/remote git state/notebook writes Ryan-gated
(exceptions already granted: round-2 twins are DONE; cancelling a losing
twin of the three pairs is pre-authorized). Never run threelevel SFS
again (Ryan 2026-09-22). orc login julia is 1.12 — never let it touch a
campaign Manifest (spack `julia/1.11.7-6bmogfl` for env ops). VTP
particle coords: rotor axis is X — radial is sqrt(y²+z²). Pre-existing
dirty files (018/026/032 docs, rotor_multi slurm script,
pressure-comparison TOML) are expected — commit only with approval.
