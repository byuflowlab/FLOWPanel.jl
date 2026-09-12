# Phase 24 — context reset: Phase 3 mostly landed, storage saga mid-flight, archiver v2 undeployed

Written 2026-09-08 morning. Supersedes `phase_23_context_reset_prompt.md` as the
entry point. Binding rules: `decision_rules.md` tail (2026-09-05+); dated record:
`ledger.md`. Pin of record unchanged: gen2d (`021-lg-gen2d` @ FLOWPanel
`00390f7`; worktree-d `/home/rander39/wt021/FLOWPanel.jl-d`; details in phase_23).

## 1. Phase 3 fleet (as of 2026-09-08 ~03:00Z; judge by outputs, never sacct)

| rung | ckpt | arms status |
| --- | --- | --- |
| R1 | 13603453 done | **all 6 arms COMPLETE** (unsteady.csv 1637 rows, clean) |
| R2 | 13603519 done | backslash/gmres/jacobi/**ilu** complete (ilu all 3 modes in CSV); fgmres_fgs 13603524 + fgs 13603525 pending |
| R3 | 13603460 done | backslash complete; gmres cold+prev in CSV, extrap mid-run (13603462); jacobi/ilu/fgmres_fgs/fgs 13603463-66 pending |

Arms serialized with afterany per rung (shared unsteady.csv, APPENDS — rerun
duplicates rows; last-wins per (solver,guess,step)). No walltime risk observed
(longest job <5 h vs 24 h). Contingency if one ever times out: resubmit on
physics2 `--qos=standby` with longer `--time` (the one permitted override).

## 2. Phase 2 tail (worktree-c; do NOT edit while jobs run)

- R4 tune 13593018 COMPLETE (p=15/0.6/32). **b0 arm submitted 2026-09-07 as
  13603860** (pre-approved sbatch, 8 h) — check its outcome.
- R5 complete + b0 certified. R6 13593020 RUNNING normally (5 d walltime).
- **R7 13593021 WATCH ITEM**: heartbeat alive but `tune_phase2.csv` 0 bytes
  since submission 09-05; no failure strings in logs. Possibly silently hung
  with a 7 d walltime. Surface to Ryan before touching; if its b0 fails, ask
  Ryan (only R4 was pre-approved).
- When R4–R7 + b0 arms land: full harvest via `benchmark/p021_merge.jl` from a
  checkout with `00390f7` (worktree-d qualifies); phase_17 DUPLICATE trap
  (last row per resume key wins).

## 3. Phase 3 results so far + figures (NEW)

- Quick-harvest headline (PROVISIONAL — see caveat): warm-start helps only FGS
  (mean niter R1 cold 5.93 → prev 4.82 → extrap 3.99, −33%); plain Krylov
  slightly WORSE warm (gmres +5–10% iters with extrap on R1/R2/R3);
  fgmres_fgs flat. CT identity across guess modes clean (max |ΔCT| <4.3e-5).
- **CAVEAT**: that harvest summarized the `niter` column = LAST inner solve
  (legacy diagnostic). The CSV has a dedicated `niter_first` column which is
  THE warmstart metric (see the notes field in any unsteady.csv row). Real
  harvest must use `niter_first`; conclusions may sharpen.
- Per-step cost columns exist for every step: `t_step_total`, `t_solve`,
  `t_bcerr`, `t_step_net` (= total − bcerr; use this for cost trajectories).
- **Four TikZ figures** in `BRAINSTORM/021_rotor_hover_solver_benchmarks/figures/`
  (Ryan-requested, built 2026-09-07): `p3_cost_per_step_{cold,prev,extrap}.tex`
  (t_step_net vs revs-since-restart; color=solver, style=rung) and
  `p3_total_cost_vs_rung.tex` (solid=cold, dashed=extrap). Backing CSVs in
  same-named dirs; regenerate with `figures/p3_cost_figures_extract.py` (edit
  its SRC to point at fresh rung CSVs; scp them with distinct names). When
  pending arms land: re-run extract + add the `\addplot` lines flagged in each
  .tex header comment. pdflatex-verified; artifacts gitignored.
- Observations worth carrying: FGMRES+FGS shows large once-per-rev cost humps
  (5→25 s/step, crosses above plain GMRES at peaks) in R1 cold — investigate
  before ranking it. Total-cost ordering stable: backslash < FGS ≲ ILU <
  Jacobi < FGMRES+FGS < GMRES; warm-vs-cold is second-order vs solver choice.
  Most interesting absentee: R3 fgs (does the FGS win grow with rung?).

## 4. Storage saga (018 runs; state 2026-09-08 morning)

Timeline that got us here (all in this session, 2026-09-07→08):
1. hpc-storage cycle found 12 ARCHIVED-STALE runs (~158 GiB) + nothing else
   archivable. Ryan **released the 288-step warm-start sets** on them.
2. `--resume-delete` pass: BLOCKED correctly — all 11 attempted runs failed
   byte-for-byte re-verification (chains drifted after archiving); the 12th
   (`p018_csarc_l3p0_3r_exp_nt`) was LIVE (job 13603853) despite STALE label.
3. Step comparison (tar vs disk, via INDEX.tsv kept_steps): **6 runs REWOUND**
   (old tarball holds steps disk no longer has): n4_nt144_l3p0_3r_sv
   (4319→209), n2_nt72_l3p0_3r_sv (2159→736), n2_nt72_l3p0 (2159→1439),
   l3p0_3r_sv (1011→536), l3p0 (1079→719), n2_nt72_l3p0_3r_sfs3nb (1625→1550).
   **These 6 are FROZEN awaiting Ryan** — their old tarballs are the only copy
   of the high-step data.
4. Ryan approved version-aside + re-archive for the 5 NON-rewound runs. Agent
   completed 3 (`l3p0_3r_sfs3nb`, `l3p0_3r_sv_h2p0`, `l3p0_3r_sv_s1p5` — new
   verified tarballs 2026-09-08 ~13:2xZ, old preserved as `.v20260905`) then
   was killed by a stall watchdog during a long tar. Locks are CLEAN, no
   .partial found, no archiver process left.
5. **REMAINING (approved, just unfinished): re-archive
   `p018_csarc_n2_nt72_l3p0_3r_sv_h2p0` and `p018_csarc_n2_nt72_l3p0_3r_sv_s1p5`.**
   Their old tarballs are already versioned aside (`.v20260903` / `.v20260904`
   in `/nobackup/archive/usr/rander39/FLOWPanel_runs/projects_FLOWPanel.jl/`),
   so a plain `scripts/run_archiver.sh --root /home/rander39/projects/FLOWPanel.jl
   --only <run> --apply` treats them as fresh. The s1p5 old tar was 82 GB —
   expect a LONG tar; run it detached/nohup on the login node, not through an
   agent watchdog window.
- **Disk: 574 G (df) vs 400 G policy cap**, burned overnight by the 3 live p018
  GPU writers (13603735/13603744/13603853, the `_3r_{exp_nt,nosfs,sfs3nb_m2}`
  chains, RECENT-HOT ~115 GiB and growing). Cap is policy not crash risk
  (mount 1.5 T free). After the 2 pending re-archives (~25–30 GiB yield), the
  levers are: Ryan's decision on the 6 rewound runs (~100+ GiB), and eventually
  the live chains finishing. The 3 `p3_checkpoint_R*` dirs (6.7 GiB) are
  deferred until R2/R3 arms finish.

## 5. Archiver v2 — edited locally, NOT deployed (NEW, Ryan-directed)

Ryan approved robustifications #1–4; implemented in the LOCAL repo
`scripts/run_archiver.sh` (Mac side, uncommitted), `bash -n` clean, suffix
matcher unit-tested. Cluster copy is still the old version (md5 baseline was
identical before edits). Changes:
1. `--supersede --only RUN[,...]`: sanctioned STALE exit — versions old tarball
   to `<tarball>.superseded.<utc>` (+ its `.fp`), falls through to normal
   tar/verify/trim; refuses live/hot runs BEFORE touching the archive; dry-run
   prints SUPERSEDE-PLAN. Scoped like --include-recent.
2. Verify-before-promote: tarball stays `.partial` until verify_tarball passes,
   only then `mv` to canonical. (Old order could leave an unverified file at
   the canonical name after a kill, or clobber a pre-existing tarball.)
3. Fingerprint sidecar `<tarball>.fp` (max archived body step + byte size, from
   the manifest); STALE messages now say `subtype=CONTINUED/REWOUND/DIVERGED/
   INTERRUPTED` with the right remedy named.
4. `queue_match` common-suffix test (≥12 normalized chars; catches the two
   observed false negatives incl. `_3r_exp_nt`); `--resume-delete` gains the
   liveness guard it entirely lacked (queue match or <2 h quiet → refuse,
   exit 9).

TODO: deploy via scp AFTER the two pending re-archives finish (NEVER overwrite
the script while an archiver is running — bash reads scripts incrementally);
dry-run `--all-checkouts` on the cluster to sanity-check classifications (expect
STALE subtypes on the 6 rewound runs: REWOUND requires .fp sidecars which old
tarballs lack, so they'll show the generic message until re-archived — that is
expected); then commit (repo convention: "archiver: ..." message, note deployed
date). Consider follow-ups #5 (root discovery/UNCOVERED audit) and #6
(INDEX.tsv supersede-awareness) — suggested to Ryan, not yet approved.

## 6. Open with Ryan

- Decision: the 6 rewound runs (release old tarballs' high-step data too, or
  `--supersede` after deploying v2, which preserves them versioned?).
- Ledger entry (still owed, offered): b0 root cause + seed-override relaunches
  13603412/13 + FGS CSV copy + gen2d pin + Phase 3 launch + the whole storage
  saga (facts in phase_23 §Disk, this file §4, and agent reports).
- Notebook entry for the campaign (offer, don't write).
- R7 b0 treatment / hung-tune investigation.
- Phase 3 harvest plan once arms complete (headline: niter_first cold-vs-warm
  and break-even step count per solver; metrics rulings in decision_rules.md;
  include per-step cost figures §3).

## 7. Mechanics

`ssh orc` needs a live ControlMaster socket (2FA otherwise — if it dies, ask
Ryan to run `! ssh orc true`). Slurm via `/apps/slurm/latest/bin/`. Root `du`
on /home times out — per-subdir or df. No sbatch/scancel beyond standing
approvals without asking Ryan in the moment. scp of multiple `unsteady.csv`:
use distinct local names (brace-expansion overwrites silently).
