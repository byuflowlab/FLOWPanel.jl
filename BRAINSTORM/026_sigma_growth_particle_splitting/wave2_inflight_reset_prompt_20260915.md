# 026 reset prompt — wave-2 IN FLIGHT, banner checks + harvest (2026-09-15 night)

You are picking up BRAINSTORM 026 (resolution-preserving particle splitting)
in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (+ sibling
FLOWVPM.jl). Read `CLAUDE.md` and the policies it names (HPC.md before
cluster work; `ssh orc` needs a live ControlMaster socket — ask Ryan to run
`! ssh orc echo ok` if 2FA blocks you). Local branch **fastmultipole**
(FLOWVPM local branch **flowpanel**). Predecessor:
`wave2_launch_reset_prompt_20260915.md` (same dir) — its task is DONE;
this file supersedes it. Authoritative record: §"Wave-2 launch" appended to
`p026_derisk_20260914_provenance.md` + the 2026-09-15 wave-2 LAUNCH ledger
line.

## What happened this session (all three gates ruled by Ryan and executed)

1. **Telemetry hardening (gate iii) — DONE, committed, smoked.**
   - FLOWVPM: `merge_particles!` kwargs `event_io`/`event_tag` → one CSV row
     per accepted merge pair `step,np,sigma_i,sigma_j,dist`; new cumulative
     `SIGMA_FLOOR_HITS` atomic counting sigma_guard FLOOR engagements in all
     four integrator paths (`_euler`/`_euler_exp` × scalar/broadcast,
     restricted to live `1:np` prefix — capacity/stale lanes were a real
     bug caught in testing).
   - FLOWPanel: `MergeParticles` policy carries `event_io`; driver opens
     `data/<run>/merge_events.csv` in APPEND mode (restart-safe, unlike
     monitor CSVs) under env `MERGE_EVENT_LOG` (default true) and prints a
     "Merge event log:" + "Sigma telemetry:" banner line; wake-health CSV
     gained `mean_sigma,max_sigma,floor_clamp_cum` columns appended after
     all existing/optional columns. `floor_clamp_cum` is per-process
     cumulative → resets to 0 on restart.
   - NO sigma0 columns in the merge log (no birth-σ storage exists in the
     particle field); classify vs local shed σ at harvest.
   - Verified: FLOWVPM merging 50/50, resolution-split 1055/1055,
     filament-edge-graph 477/477; FLOWPanel unit_wake 730/730, unit_replay
     142/142, unit_simulate PASS after updating its two wake-health schema
     testsets; scalar-vs-broadcast floor-count parity exact; 467-step local
     CPU driver smoke end-to-end (64k merge events, census populated,
     floor counter engaged, zero errors).

2. **Commit + tags (gate iv) — DONE.** Annotated tag
   `campaign/p026-wave2-20260915` on BOTH repos, pushed to orc:
   - FLOWPanel **4d6d2c0** (telemetry commit 212556c + bundle commit
     4d6d2c0: slurm NREVS fixes, new mo35 case def, 026/021 docs incl.
     ~1.2 G of 021 evidence dirs — largest file 77 M, precedented).
   - FLOWVPM **2b253db** (supersedes bf88806 as the campaign pin).
   - FastMultipole unchanged **ac7230a6** (`campaign/p026-derisk-20260914`).
   - **FLAG for Ryan (unresolved):** orc's `fastmultipole` and `flowpanel`
     BRANCHES were NOT pushed — each holds a 2026-09-02 "wip snapshot
     before unified-052 consolidation" commit absent from local history
     (FLOWPanel 20c9da8, FLOWVPM eaf257c). Do NOT force-push without
     Ryan's ruling; tags carry the pins so the campaign is unaffected.
   - Deliberately left uncommitted (other arcs' in-flight files):
     018 gamma_distribution_*/mechanism_tests_*/rlxfscaled_* files,
     `expguard_provenance_20260908.md`,
     `run_rotor_multi_ground_effect_gpu.slurm.sh` (022), and the
     perpetually-uncommitted
     `data/rotor_hover_pressure_comparison/...metadata.toml`.
   - NEW uncommitted since the tag: provenance §Wave-2-launch append,
     ledger line, this file → bundle into the next Ryan-gated commit.

3. **Worktrees + submission (gate v) — DONE.** Campaign worktrees under
   `orc:~/campaigns/p026-derisk-20260914/` fast-forwarded (FLOWPanel.jl →
   4d6d2c0, FLOWVPM.jl → 2b253db, both clean/detached); env Manifest
   dev-paths verified unchanged. Pre-flight confirmed no stale run
   dirs/symlinks for the 8 case names. Submitted 2026-09-15 ~21:45 from
   cwd = campaign FLOWPanel worktree, wrapper
   `~/projects/FLOWPanel.jl/examples/run_p018_screen_gpu052.slurm.sh h200
   <case>` with `P018_REPO_OVERRIDE`/`P018_PROJECT_OVERRIDE` at the
   campaign worktree/env:

   | job | case | mechanism | submission env beyond case def | wall |
   |---|---|---|---|---|
   | 13712311 | scr_p026s9_ctrllg_floor | floor | FLOOR | 8 h |
   | 13712312 | scr_p026s9_ctrllg_split | split | SPLIT | 8 h |
   | 13712313 | scr_p026s9_ctrllg_fs | floor+split | FLOOR+SPLIT | 8 h |
   | 13712314 | scr_p026s9_explg_floor | floor | FLOOR | 8 h |
   | 13712315 | scr_p026s9_explg_split | split | SPLIT | 8 h |
   | 13712316 | scr_p026s9_explg_fs (lg merge-A) | floor+split | FLOOR+SPLIT | 8 h |
   | 13712317 | scr_p026s9_explg_fs_mo35 (lg merge-B) | fs + MERGE_OVERLAP=3.5 in case def | FLOOR+SPLIT | 8 h |
   | 13712318 | scr_p026sp_nt144_cap018 | grow cap σ_max=0.018 | FLOOR+SPLIT | 24 h |

   FLOOR = `SIGMA_FLOOR_FRAC=0.1`; SPLIT = `WAKE_SPLIT_VISCOUS=true
   WAKE_SPLIT_FRAC_VISCOUS=0.587 WAKE_SPLIT_FRAC_COMPRESS=0.73
   WAKE_SPLIT_FRAC_ELONGATE=0.3`. At session end: 11/12/13 RUNNING on
   m13h, 14–18 PENDING (QOSMaxCpuPerUser — they cycle in as arms finish,
   s9 arms ≈1–1.5 h each). Slurm .out files land in the campaign
   FLOWPanel worktree (`slurm-<jobid>.out`).

## YOUR TASK (in order)

1. **Banner verification, MANDATORY, all 8 jobs** (a session-local watcher
   died with the old session — re-check yourself; delegate log reading to
   `hpc-monitor`). For each `slurm-<jobid>.out` in the campaign worktree
   verify: nrevs=12 (cap018: 17) in "Total run length" (468 steps s9 /
   2592 steps cap018 incl. spinup); `FLOWPANEL_FILAMENT_REG` → reg
   linegauss; `SIGMA_FLOOR_FRAC=0.1 (guard=on)` ONLY on floor/fs/cap018
   arms (split-only arms: floor 0.0, WAKE_SPLIT_SIGMA_MIN prints "off" —
   intentional); split knobs per slate (viscous 0.587 / compress 0.73 /
   elongate 0.3 on split/fs/cap018; ALL off on floor-only arms — driver
   hard-errors on mismatch, so a running sim past the banner implies
   consistency); "Merge event log: data/<case>/merge_events.csv" +
   "Sigma telemetry:" lines present. The NaN log-gate greps the word NaN —
   split knobs print "off" when unarmed. Record verdicts in a provenance
   append. If a banner is WRONG: do not cancel without Ryan unless it is
   unambiguously mis-armed (wrong mechanism knobs) — then scancel, fix,
   resubmit, and record.
2. **Watch runs** (hpc-monitor; judge by outputs, sacct state is not
   evidence). Watch for: f_visc=0.587 firing (never fired in grow regime,
   fired late in shrink); floor_clamp_cum engagement on floor arms;
   monitor04 wake-health CSV does NOT append across restarts (but
   merge_events.csv DOES); WakeGeometryError-class deaths.
3. **Harvest when arms finish** (delegate to `harvester`): MOVE run dirs
   from the worktree `data/` to shared root `orc:~/projects/FLOWPanel.jl/
   data/` + symlink back (NEVER pre-place symlinks — dispatcher wipes
   data/$RUN_NAME on cold runs). Key comparisons: (a) three-way mechanism
   attribution floor vs split vs fs, ctrl and exp families; (b)
   vatistas-vs-linegauss merge verdict re-anchor = arms 13712316/17 vs
   de-risk pair `scr_p026s9_exp_split`/`_mo35` (shared root); (c) cap018
   vs cap030 (both NREVS=17-class, cap030 data in shared root) for the
   grow-side cap ladder; (d) NEW telemetry: merge_events.csv event
   rates/σ distributions, floor_clamp_cum engagement curves, mean/max
   sigma census vs the σ-pump hypothesis.
4. Gate any conclusions per the global CLAUDE.md evidence rule; offer (not
   write) a notebook entry — the whole 026 arc is still unlogged (Ryan
   "not yet" ×2), plus older owed arcs (expguard three-arm, Cd transient,
   rlxf derivation, Ladder C forensics).

## Owed / parked (carried)

- orc branch divergence ruling (see FLAG above).
- scr_p026gpuv_split archive retry via hpc-storage once quiet ≥24 h;
  `data/p026_restart_gpu40_s950` left for the archiver.
- 021 silo cleanup owed (see 021 package §7); unrelated p021 job may be
  queued — leave it alone. 018mech/022g jobs in queue are other sessions'.
- Next commit bundle: provenance/ledger appends + this file.

## Ground rules (carry-over)

Local ≤4 threads; Slurm needs `bash -lc` in non-login shells; MOTD
banners contaminate ssh output — filter. Commits, campaign launches,
notebook writes are Ryan-gated. Known unrelated failure:
`runtests_unit_warmstart.jl` first testset. Judge runs by outputs (.out
for sim output).
