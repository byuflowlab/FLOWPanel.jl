# Phase 25 — context reset: re-archives DONE, archiver v2 DEPLOYED (uncommitted), 4 Phase-3 arms left, R6+R7 both hung

Written 2026-09-08 ~19:45Z. Supersedes `phase_24_context_reset_prompt.md` as the
entry point. Binding rules: `decision_rules.md` tail (2026-09-05+); dated record:
`ledger.md`. Pin of record unchanged: gen2d (`021-lg-gen2d` @ FLOWPanel
`00390f7`; worktree-d `/home/rander39/wt021/FLOWPanel.jl-d`).

## 1. Phase 3 fleet (verified 2026-09-08 ~19:30Z, by outputs + squeue)

| rung | status |
| --- | --- |
| R1 | all 6 arms COMPLETE (unchanged) |
| R2 | backslash/gmres/jacobi/ilu complete in CSV (219 rows per Krylov solver, 541 backslash). **fgmres_fgs 13603524 RUNNING step 71/72** (17:00 elapsed, 7:00 walltime left — tight but should finish); fgs 13603525 PENDING behind it |
| R3 | backslash/gmres/jacobi/ilu COMPLETE incl. extrap (13603462/63/64 all PASS, 70 wake-developed rows each; net t_solve gmres 5064.7 s, jacobi 1974.3 s, ilu 1634.7 s). **fgmres_fgs 13603465 RUNNING step 70/72** (11:24 elapsed); fgs 13603466 PENDING behind it |

So 4 arms outstanding: 2 fgmres_fgs finishing (~hours), then the 2 fgs arms
(the most interesting absentee: does the FGS warm-start win grow with rung?).
Walltime contingency unchanged: resubmit on physics2 `--qos=standby` with
longer `--time` (the one permitted override). Watch R2-fgmres_fgs's 7 h margin.

## 2. Phase 2 tail — **R6 AND R7 BOTH HUNG (new finding; surface to Ryan)**

- R4 b0 arm 13603860 **PASS** (completed 2026-09-08 02:07Z; wrote
  `wt021/FLOWPanel.jl-c/benchmark/results/phase2/multi/R4/tune_phase2.csv`,
  1 trace point (17, 0.65, 6) → 143.4 s). R5 complete + b0 certified.
- **R6 13593020**: previously believed running normally — it is NOT. Same
  pathology as R7: Julia heartbeats alive (5-min cadence), log shows header
  only ("021 Phase 2 real-solve tuning"), tune_phase2.csv never created, zero
  progress in 37+ h. 3d10h walltime left.
- **R7 13593021**: heartbeat alive, tune_phase2.csv 0 bytes since 2026-09-05
  15:14, zero progress 37+ h. 5d10h walltime left.
- Both launched Julia but never reached a solve step — likely a common startup
  hang. NO ACTION without Ryan (kill/resubmit not pre-approved).
- **Investigation plan (run read-only steps freely; kill/resubmit needs Ryan)**:
  1. Read both jobs' stderr/stdout fully (slurm-1359302{0,1}.out in the
     submission dir under `wt021/FLOWPanel.jl-c`) — look for the LAST line
     before silence: precompile banner, package load, mesh build, first tune
     point. R4/R5 logs are the healthy reference; diff the startup sequence.
  2. Suspect list, in order: (a) Pkg/precompile lock contention — R6+R7
     launched simultaneously 2026-09-05 sharing worktree-c's Julia depot;
     check `~/.julia/logs/manifest_usage.toml` mtimes and any
     `*.ji.pidfile`/`P.lock` under `~/.julia/compiled` from 09-05 15:1x;
     (b) both tuning at large rungs → first real solve simply enormous
     (but 37+ h with no CSV header row makes this unlikely — R6's CSV was
     never even created); (c) an interactive prompt/2FA-style stall in the
     driver (grep the driver for `readline`/`ask`).
  3. `scontrol show job` for both → node names; `ssh <node> ps -o etime,args
     -u rander39` and `py-spy dump`/`gdb -p` if available to see where Julia
     is parked (read-only, allowed).
  4. If a depot lock is confirmed: the fix is stagger + per-job
     `JULIA_DEPOT_PATH` overlay or pre-warmed precompile; propose resubmit of
     ONE (R7, longer walltime) to Ryan before touching R6.
  5. Heartbeats come from the harness shell loop, not Julia progress — do not
     read them as liveness (add to decision_rules if Ryan agrees).
- When R4–R7 + b0 land: full harvest via `benchmark/p021_merge.jl` from a
  `00390f7` checkout (worktree-d); phase_17 DUPLICATE trap (last row per
  resume key wins).

## 3. Phase 3 results + figures (unchanged from phase_24 §3)

Provisional headline (from legacy `niter` column): warm-start helps only FGS
(R1 mean 5.93→4.82→3.99, −33%); plain Krylov slightly worse warm; CT identity
clean (<4.3e-5). **Real harvest must use the `niter_first` COLUMN** (the
warmstart metric — not first-step niter). Cost columns: use `t_step_net`.
Four TikZ figures `figures/p3_cost_per_step_{cold,prev,extrap}.tex` +
`p3_total_cost_vs_rung.tex`; when the 4 pending arms land, re-run
`figures/p3_cost_figures_extract.py` (edit SRC; scp fresh rung CSVs with
distinct names) and add the `\addplot` lines flagged in each .tex header.
Investigate FGMRES+FGS once-per-rev cost humps before ranking it. Harvest
plan to offer Ryan: headline niter_first cold-vs-warm + break-even step count
per solver (metrics rulings in decision_rules.md) + per-step cost figures.

## 4. Storage saga — RESOLVED items this session (2026-09-08 ~19:20–25Z)

1. **Both approved re-archives are DONE and verified** (run detached on login
   node, log `/home/rander39/rearchive_20260908.log`):
   - `p018_csarc_n2_nt72_l3p0_3r_sv_h2p0`: 3178 files, src 10.2 GB → tar
     7.9 GB, kept [699–703], freed 9.9 GB, verify PASS.
   - `p018_csarc_n2_nt72_l3p0_3r_sv_s1p5`: first attempt refused by the 24 h
     quiet guard (`RECENT quiet=21h`); root-caused the mtimes to the previous
     session's own verify/trim churn (newest real data = step 2159 outputs,
     no queue job) → re-ran with `--include-recent --only` (the sanctioned
     scoped override). 3181 files, src 17.5 GB → tar 13.5 GB, kept
     [2155–2159], freed 17.2 GB, verify PASS.
   - Old tarballs remain versioned aside (`.v20260903`/`.v20260904`).
2. **Disk: 548 G** (df, after re-archives) vs 400 G policy cap. Remaining
   levers: Ryan's decision on the 6 rewound runs (~100+ GiB), live chains
   finishing, and the 3 `p3_checkpoint_R*` dirs (6.9 GiB, deferred until
   R2/R3 arms finish — they showed RECENT quiet=19–23h in today's dry-run).
3. The 6 rewound runs remain FROZEN awaiting Ryan (old tarballs are the only
   copy of their high-step data): n4_nt144_l3p0_3r_sv (4319→209),
   n2_nt72_l3p0_3r_sv (2159→736), n2_nt72_l3p0 (2159→1439), l3p0_3r_sv
   (1011→536), l3p0 (1079→719), n2_nt72_l3p0_3r_sfs3nb (1625→1550).

## 5. Archiver v2 — DEPLOYED to cluster 2026-09-08, dry-run clean, COMMIT PENDING

- Deployed via scp after re-archives finished; md5 `754484a1...` matches local;
  `bash -n` clean on cluster.
- Dry-run `--all-checkouts` results (all expected): alias-dedup collapsed all
  14 checkouts/worktrees onto the one data root; STALE_COUNT=7 = the 6 rewound
  runs + `l3p0_3r_exp_nt` (its writer 13603853 finished cleanly this morning,
  so it is now genuinely CONTINUED-stale); all 7 show the generic ASK-RYAN
  message because pre-v2 tarballs lack `.fp` sidecars — expected, subtypes
  appear only after re-archiving under v2. RECENT_COUNT=12 (live/recent p018
  chains + p3 checkpoints), VERIFY_FAIL=0, no false ALIAS results.
- **TODO next session: commit** — stage ONLY `scripts/run_archiver.sh` (repo
  has many unrelated modified files). Suggested message:
  `archiver: v2 — supersede flow, verify-before-promote, .fp sidecars + STALE subtypes, suffix queue-match, resume-delete liveness guard (deployed to cluster 2026-09-08)`.
- **New wart found**: on the login node `squeue` is not in PATH, so the
  archiver logs "squeue unavailable — liveness rests on mtime alone" and v2's
  queue_match/liveness guards silently degrade to mtime-only. One-line fix:
  fall back to `/apps/slurm/latest/bin/squeue` if `squeue` absent. Suggest to
  Ryan alongside follow-ups #5 (root discovery/UNCOVERED audit) and #6
  (INDEX.tsv supersede-awareness) — none approved yet.

## 6. p018 GPU chains — RESOLVED: chained Slurm segments; ALL STOPPED, 4 arm deaths root-caused

The "detached writer" mystery is resolved: the p018 GPU chains run as
**short chained Slurm segments** (1.5–6 h each, new job ID per segment, same
job name; logs at `wt018/<worktree>/logs/slurm/slurm-<name>-<id>.{out,err}`).
The phase_24 job IDs (13603735/744/853) were segments that ended 09-07; the
09-08 segments were 13605974/983/984/985/986. **As of 2026-09-08 ~19:40Z NO
p018 job is running or queued — every chain is terminal or stalled:**

| chain | last segment | outcome |
| --- | --- | --- |
| n2_nt72_3r_exp_nt | 13605974 | **COMPLETE** step 2159/2159, gate_rc=0, clean (GH200) |
| l3p0_3r_exp_nt | 13603853 | COMPLETE 1080 steps, gate_rc=0 |
| l3p0_3r_cs0p18 | 13605983 | COMPLETE |
| n2_nt72_3r_cs0p18 | 13605984 | **DIED step 1876/2159** — FMM σ-gate (below) |
| n2_nt72_3r_nosfs | 13603744 | **STALLED at step 1458/2159 since 09-07** — FMM σ-gate |
| l3p0_3r_cs0p002_exp_nt | 13605985 | **DIED step 898/1079** — dt·\|L\|=8.6e9 exceeds euler_exp substep budget (Γ-ignition blow-up) |
| l3p0_3r_cs0p34_exp_nt | 13605986 | **DIED step 361/1079** — dt·\|L\|=8.2e4 exceeds euler_exp substep budget |

- The σ-gate deaths are the KNOWN 018 blocker (FastMultipole "regularized
  nearfield near-set adequacy failed": σ_max ≈ 0.19–0.20 grew past the ell=2
  direct-stencil cutoff, "admissible depth ell <= 1"; remedy = BRAINSTORM 026
  particle splitting, Phase 2 in flight). nosfs died at 1458, cs0p18 at 1876
  — same class as the deterministic ~step-1550 gate.
- The euler_exp budget deaths are physical/numerical blow-ups in the new
  const-Cs SFS ladder arms (cs0p002 ≈ no SFS, cs0p34 heavy) — consistent with
  the known Γ-ignition regime map.
- **Surface to Ryan**: whether to rewind+restart the σ-gate arms with reduced
  ell / smaller σ, wait for 026 splitting, or accept truncated series; and
  whether the cs-ladder blow-ups kill those arms or warrant a schedule change.
- Storage note: with all chains stopped, the RECENT classifications in §5's
  dry-run will age into archivable within ~24 h; the completed exp_nt chains
  are large fresh archive candidates once Ryan confirms no immediate restart.

## 7. Open with Ryan

- Decision: the 6 rewound runs (release old high-step data, or `--supersede`
  under v2 which preserves them versioned).
- **NEW: R6+R7 both hung** — approve kill/diagnose/resubmit? (Only R4 b0 was
  pre-approved.)
- Archiver v2 commit approval (+ squeue-PATH fix, follow-ups #5/#6).
- Ledger entry (owed): b0 root cause + seed-override relaunches 13603412/13 +
  FGS CSV copy + gen2d pin + Phase 3 launch + full storage saga (facts:
  phase_23 §Disk, phase_24 §4, this file §4–5).
- Notebook entry for the campaign (offer, don't write).
- Phase 3 harvest plan once the 4 arms complete (§3).

## 8. Mechanics

`ssh orc` needs a live ControlMaster socket (2FA otherwise — ask Ryan to run
`! ssh orc true`). Slurm via `/apps/slurm/latest/bin/`. Cluster local time =
UTC−6 (ls timestamps are local; `date -u` for logs). Root `du` on /home times
out — per-subdir or df. No sbatch/scancel beyond standing approvals without
asking Ryan in the moment. scp of multiple `unsteady.csv`: distinct local
names. Never overwrite `run_archiver.sh` while an archiver process runs.
Re-archive log: `/home/rander39/rearchive_20260908.log`.
