# 018 GPU — ladders RELAUNCHED 2026-08-31 ~01:05 (context-reset handoff)

## 0. STATE OF PLAY as of 2026-08-31 ~16:45 (read this, then the sections below for detail)

Results so far (CT̄ over the 10-rev convergence window, revs 21–30):

| rung | merge-off `_nm` | fixed-rlxf 0.3 | with-merge ref (Aug 29) |
|---|---|---|---|
| NT36 | 0.070811 (done) | 0.0707752 (=13507289) | 0.0707752 |
| NT72 | 0.071934 (done) | 0.071709 (done) | 0.0718444 |
| NT144 | 13518981 RUNNING m13h-1-2, ETA ~23:00+ | 13518980 RUNNING m13h-2-2, ETA ~20:30 (partial revs 15–20: CT ≈ 0.0732–0.0733) | 0.073829 |

Both interventions look like NULL results: the NT climb (~+1.5–2%/doubling)
persists with merging off AND with rlxf pinned at 0.3. NT144 completions
tonight confirm. Readout recipe: `<run>_CT_per_rev.csv` in the run's data
dir, mean of `CT_mean` over `in_convergence_window == true` rows (10 revs);
mid-run, use `monitors/*_monitor02_force_system1.csv`, CT = −CFx averaged
per 144/72/36-step rev block.

Traditional-Pedrizzetti (`RELAX_SCHEME=pedrizzetti`) ladders queued —
see the section below for the full table, twin pairs (m13h↔eng:
773↔903, 774↔904, 775↔905, 776↔906; 777 mgh no twin), cancel-the-loser
policy, and metadata verification recipe.

MONITORING after a context reset: background tasks die with the session.
Plain `run_in_background` watcher loops were being externally killed midday;
the Monitor tool (persistent watch task) proved reliable — re-arm it with
the current JOBS/pending list (poll `sacct -X -P` via `ssh orc bash -lc`
every 900 s; wake on any terminal state or any pending job turning RUNNING).
`watch_018_ladders_20260831.sh` in this dir is the older script equivalent
(update its JOBS/UNVERIFIED before reuse).

Node notes: mgh-1-2 = 3× NODE_FAIL under NT144 load today, avoid + ticket-
worthy (rc.byu.edu/ticket, not yet ticketed). mgh-1-1 rebooted clean
(0 MiB, usable — 13542777 targets it). m13h user cap: QOS `m13h` cpu=192;
eng-1-1 CPU-blocked ~24 h by Ryan's other campaign job 13542825 (see below).
Other active runs under this account (do NOT touch): 13542825
`fp-il-s038v-lgpz-e` (eng), 13518479 `fp-il-s038v-cpu40-r3` (m12), plus any
`fp-il-*`/`fp052*` jobs — different campaigns.

Written for the next agent. Read this FIRST. Supersedes the pause doc
`gpu_paused_ladders_status_20260830.md` (its §3 ladder rationale, §4 session
changes, and §6 guardrails ALL still apply — do not re-derive). Companions:
`gpu_nomerge_ladder_status_20260829.md` (3R results, merge-error hypothesis,
recipes §4, guardrails §5) and `phase_17_nprop_nt_ladder.md` last three
sections (NT-climb analysis, brainstorm prep).

## 1. Current jobs (all five relaunched 2026-08-31 ~01:05 by Ryan's order —
he cut the 24 h pause short)

| job | arm | ladder | state as of 04:30 | verified |
|---|---|---|---|---|
| ~~13518894~~ | NT144 λ3.0 | merge-off `_nm` | **NODE_FAIL 07:06 on mgh-1-2** (2h43m; 2nd mgh-1-2 NODE_FAIL) | was verified |
| ~~13518931~~ | NT144 λ3.0 | merge-off `_nm` | **NODE_FAIL 12:18 on mgh-1-2 (4h44m) — THIRD mgh-1-2 strike (2h12/2h43/4h44). mgh retries stopped per plan; eng 13518932 (est 09-01 15:00) + m13h 13518981 (behind own 192-cpu cap) are the live copies.** | was verified |

12:2x updates: mgh-1-1 shows `reboot^` ("bootaction 538") — admin reboot in
progress, which should clear the 61.6 GB leak; once it returns idle, probe
memory and if clean (<1 GB) resubmit NT144-nm there with `--exclude=mgh-1-2`
(mgh-1-2 is the 3-strike node; mgh-1-1's suspension was leak-conditional).
13518485 (NT72 rlxf eng spare) started 10:04 eng-1-1, VERIFIED banner
`rlxf:0.3 filament_reg:gaussian NT:72` — races m13h copy 13518979.
13518486 (NT144 rlxf eng spare) est start today 20:04.

~13:0x: NT72-rlxf race resolved — m13h copy 13518979 COMPLETED first, eng
twin 13518485 scancel'd per Ryan. The freed slots let both NT144 spares
start: 13518981 (nm `_m13h`) RUNNING m13h-1-2, 13518486 (rlxf) RUNNING
eng-1-1 — BOTH VERIFIED ~13:2x (981: fresh `_nm_m13h` metadata `every = 0`;
486: banner `NT:144 rlxf:0.3 filament_reg:gaussian`). mgh-1-1 moved reboot→`maint` (still
unschedulable; probe memory when it returns to idle before any use).
Remaining queued: 13518932 (nm eng backup) only.

~13:5x: NT144-rlxf race decided by Ryan's order — at equal step 621 the m13h
copy 13518980 was ~10% faster/step (3.37 s vs 3.73 s) and 2123 steps ahead
(2744/4319 vs 621/4319), so eng copy 13518486 scancel'd. 13518980 ETA ~20:30
(last steps ~15 s). Live: 13518980 (rlxf m13h), 13518981 (nm m13h),
13518932 (nm eng backup — started 13:35 on freed eng-1-1, VERIFIED
`every = 0` ~14:1x).

~14:3x per Ryan: nm race pre-killed after all — at equal step 779 the m13h
copy 13518981 was ~3% faster/step (3.93 s vs 4.06 s) and 586 steps ahead, so
eng copy 13518932 scancel'd. NT144-rlxf partial readout (force monitor,
revs 15–20): CT ≈ 0.0732–0.0733, stable — fixed-rlxf does NOT flatten
NT72→NT144 (+2.1% vs rlxf NT72 0.0717087; with-merge NT144 was 0.073829,
so only −0.7%, within with-merge NT144's ±1.31% scatter). Live jobs:
13518980 (rlxf m13h, ETA ~20:30), 13518981 (nm m13h, ETA ~23:00–24:00).

## Traditional-Pedrizzetti ladders (Ryan order 2026-08-31 ~15:30)

Confirmed the runs use `pnl.FLOWVPM.relaxation_correctedpedrizzetti`
(driver `examples/rotor_hover_pressure_comparison.jl:469`). Added a
`RELAX_SCHEME` env hook to the driver (default `correctedpedrizzetti`,
unset behavior unchanged; `pedrizzetti` selects the traditional scheme, no
`/sqrt(b2)` renorm) + a `relax_scheme = "..."` line in the run metadata.
Edited locally, parse-checked, rsync --checksum to BOTH silos (local and
silo drivers were line-identical pre-edit). Dispatcher untouched —
`RELAX_SCHEME` rides `--export=ALL,...`.

Two ladders, merge ON, λ3.0, 3R, settle 22 (mirrors the with-merge column
and today's fixed-rlxf ladder; NT36 rung SHARED since progressive rlxf at
NT36 = 0.3 = the fixed value). Progressive rlxf = case defaults
(NT36 0.3 / NT72 0.16334 / NT144 0.08539).

| job | arm | run name | where | time |
|---|---|---|---|---|
| 13542773 | NT36 tp (shared) | `p018_csarc_l3p0_3r_tp` | m13h | 6h |
| 13542774 | NT72 tp prog | `..._n2_nt72_l3p0_3r_tp` | m13h | 10h |
| 13542775 | NT72 tp rlxf0.3 | `..._n2_nt72_l3p0_3r_tp_rlxf0p3` | m13h | 10h |
| 13542776 | NT144 tp rlxf0.3 | `..._n4_nt144_l3p0_3r_tp_rlxf0p3` | m13h | 18h |
| 13542777 | NT144 tp prog | `..._n4_nt144_l3p0_3r_tp` (gh200 silo) | mgh-1-1 (`--exclude=mgh-1-2`) | 24h |

All PENDING at submit (~16:00): 777 est 17:35 today; m13h ones est
overnight, will pull in as 13518980/981 finish. VERIFY on start via run
metadata: `relax_scheme = "pedrizzetti"` AND expected `relax_rlxf`
(and `merge_particles = true`).

~16:2x eng twins (Ryan: "request eng in addition, cancel one after the
other starts"): 13542903 (NT36 `_e`), 13542904 (NT72 prog `_e`), 13542905
(NT72 fix `_e`), 13542906 (NT144 fix `_e`) — qos=eng, `-C hopper`, run
names suffixed `_e`. Twin pairs: 773↔903, 774↔904, 775↔905, 776↔906
(777 has no twin). On first RUNNING of a pair member: scancel the twin,
then verify. Note eng-1-1 is currently CPU-blocked ~24h by Ryan's OTHER
campaign job 13542825 (`fp-il-s038v-lgpz-e`, 64c) + st392's 16c eng job —
preemption (gstandby only) can free at most 32 of the needed 64 cpus, so
the m13h copies will likely win; eng twins are insurance. A brief
scontrol move of 775/776 to eng (before the twin decision) was reverted.

~17:5x: **13542777 STARTED mgh-1-1 17:36 and VERIFIED** (no twin). Verification
method CORRECTION: the `relax_scheme = "..."` line lands in
`<run>_case_metadata.toml`, which the driver writes only at END of run
(driver line ~1543) — useless for start-time checks. Instead verify via the
big step-metadata file `<run>.metadata.toml` (fresh mtime > job start):
relaxation `type = "FLOWVPM.relax_pedrizzetti"` (default scheme would say
`relax_correctedpedrizzetti`), `rlxf = <expected>` (777 showed 0.08539 ✓),
and MergeParticles `every = 1` (merge ON ✓). Use this same recipe for
773–776/903–906. Also: a leftover watcher task from the prior session was
found alive and stopped; single fresh Monitor rearmed (poll 900 s).

~15:0x: mgh-1-1 back idle after admin reboot; probe shows **0 MiB used** —
the 61.6 GB leak is CLEARED, mgh-1-1 usable again (no ticket needed for it).
mgh-1-2 remains the suspect node: 3× NODE_FAIL under NT144 load today
(2h12/2h43/4h44) while reporting healthy — still ticket-worthy on its own.
No resubmits needed; both NT144 arms run on m13h.
| 13518932 | NT144 λ3.0 | merge-off `_nm` **eng backup** | PENDING eng (est. 09-02 05:00); run name `p018_csarc_n4_nt144_l3p0_3r_nm_eng`, h200 silo | not started |

Ryan's order 2026-08-31 ~08:00: race 13518931 (mgh) vs 13518932 (eng backup,
distinct run name `*_nm_eng` in `~/FLOWPanel-018-gpu-h200`); keep whichever
finishes, cancel the other once one completes cleanly.

Ryan's second order ~08:05: m13h (25/32 H200 free via `--qos=gpu`) starts
almost immediately, so the pending eng arms were double-queued there too —
keep the winners. Same h200 silo, run names suffixed `_m13h`. m13h nodes lack
the `hopper` feature and the slurm script's baked-in `#SBATCH --constraint=arm`
otherwise applies, so m13h submits need `--constraint=intel`. Partition QOS
`m13h` caps each user at cpu=192 (other running m13h jobs 13518861/13518964
count against it).

| job | arm | state ~08:15 |
|---|---|---|
| 13518484 | NT36 nm (eng) | **COMPLETED ~10:07, exit 0:0, 2h13m — CT̄(10-rev window) 0.0708113, half-range ±0.11% (vs with-merge 0.0707752: +0.05%, unchanged)** (was VERIFIED `every = 0`) |

Merge-off ladder readout so far (2 of 3 rungs): NT36 0.070811, NT72 0.071934
— both within 0.15% of the with-merge column, so the NT36→NT72 climb (+1.6%)
persists with merging off. Early evidence AGAINST the merge-error hypothesis;
NT144 nm (13518931 mgh / 13518932 eng / 13518981 m13h) will decide.
| 13518979 | NT72 rlxf `_m13h` | **COMPLETED 12:46, exit 0:0, 4h49m — CT̄(10-rev window) 0.0717087, half-range ±0.20%. vs NT36 rlxf rung 0.0707752: +1.3% → fixed-rlxf does NOT flatten NT36→NT72 either.** Eng twin 13518485 cancelled ~12:5x per Ryan (was at step 1479/2159) — m13h copy won the race. |
| 13518980 | NT144 rlxf `_m13h` | RUNNING m13h-2-2 since 08:03 (trimmed 64→48 cpus for the 192-cpu cap; banner still says threads:64 — possible mild oversubscription, perf-only); VERIFIED banner `rlxf:0.3` + `filament_reg:gaussian` + `NT:144` ~08:2x |
| 13518981 | NT144 nm `_m13h` | PENDING QOSMaxCpuPerUser at 64 cpus (fine — mgh copy 13518931 already running verified) |
| 13518483 | NT72 λ3.0 | merge-off `_nm` | **COMPLETED 07:55, exit 0:0, 5h38m — CT̄(10-rev, revs 21–30) 0.0719342, half-range ±0.22% (vs with-merge 0.0718444: +0.13%, i.e. NT72 rung ~unchanged by merge-off)** | YES: `every = 0` |
| 13518484 | NT36 λ3.0 | merge-off `_nm` | PENDING eng | not started |
| 13518485 | NT72 λ3.0 rlxf=0.3 | fixed-rlxf gaussian | PENDING eng | not started |
| 13518486 | NT144 λ3.0 rlxf=0.3 | fixed-rlxf gaussian | PENDING eng | not started |

Verification recipe (per arm, once RUNNING): merge-off arms → fresh
`data/p018_*_3r_nm/*.metadata.toml` (in the arm's silo) must show `every = 0`
inside the `MergeParticles` functional-policy block (the `every` line sits a
few lines ABOVE `type = "MergeParticles"` — grep -B6). rlxf arms → run banner
in the slurm .out (~20 lines after "DJI 9443 Phase 2e") must print `rlxf:0.3`
and `filament_reg:gaussian`.

## 2. NT144 saga (why it is job #3, and why mgh-1-1 is OFF LIMITS)

- 13518482 (first relaunch) ran 2h12m on mgh-1-2, verified merge-off, then
  **NODE_FAIL** at 03:21 (transient; node returned to idle, no reason set).
- Resubmit 13518892 landed on mgh-1-1 and **FAILED in 4.5 min**: `insufficient
  free GPU memory for resident S: free≈36.1 GB, required≈45.2 GB`.
  Probe (`srun -p mgh -w mgh-1-1 -C arm --gres=gpu:gh200:1 --mem=8G
  --cpus-per-task=2 nvidia-smi`) showed **61,624 MiB used with NO processes**
  — leaked/orphaned GPU memory; Slurm still reports the node idle. Only a
  reboot/admin can reclaim it.
- 13518894 resubmitted with `--exclude=mgh-1-1`, ran on mgh-1-2 — then
  **NODE_FAIL again at 07:06 (2h43m in)**. mgh-1-2 has now killed NT144 twice
  (2h12m, 2h43m), both times returning to idle with no reason set. Probe at
  ~07:50 showed mgh-1-2 GPU memory clean (0 MiB) and mgh-1-1 still leaked
  (61,624 MiB). 13518931 is the 3rd attempt (same line, `--exclude=mgh-1-1`);
  if mgh-1-2 NODE_FAILs a 3rd time, stop retrying mgh — move NT144-nm to eng
  (h200 silo `~/FLOWPanel-018-gpu-h200`, pattern per the rlxf NT144 line in
  `~/relaunch_018_ladders_20260830.sh`) and ticket ORC for BOTH mgh nodes.
- Consequences: (a) do NOT hop anything to mgh-1-1 until its memory is
  confirmed clean (re-probe with the srun above; expect <1 GB used); the
  pause-doc §5 mgh-1-1 hop recipe for NT36 is SUSPENDED. (b) Any future mgh
  resubmit should carry `--exclude=mgh-1-1` until then. (c) Ryan may want to
  ticket ORC (rc.byu.edu/ticket): mgh-1-1 GH200, 61.6 GB used / 0 processes,
  observed 2026-08-31 04:23. Not yet ticketed. NODE_FAIL restarts are FRESH
  (dispatcher wipes the partial `_nm` dir) — no `_s2` chaining for those;
  `_s2` chaining is only for exit 124 (timeout with valid partial).

Verbatim NT144 resubmit line (from `~/FLOWPanel-018-gpu-gh200` on orc; other
arms' lines are in `~/relaunch_018_ladders_20260830.sh`, still on orc):

    sbatch --job-name=fp-018gpu-n4_nt144_l3p0-3r-nm --partition=mgh \
      --gres=gpu:gh200:1 --constraint=arm --qos=normal --no-requeue \
      --exclude=mgh-1-1 --cpus-per-task=72 --mem=192G --time=24:00:00 \
      --output=logs/slurm/slurm-%x-%j.out --error=logs/slurm/slurm-%x-%j.err \
      --export=ALL,SIGMA_CEIL=0.030,TRUNCATION_RADIUS_R=3.0,MAX_PARTICLES=1500000,P018_SETTLE_REVS=22,MERGE_PARTICLES=false,P018_RUN_NAME=p018_csarc_n4_nt144_l3p0_3r_nm \
      examples/run_dji9443_hover_ct_gpu.slurm.sh gh200 p018_csarc_n4_nt144_l3p0

## 3. Monitoring (restart this after a context reset — background tasks die
with the session)

Script: `watch_018_ladders_20260831.sh` (this directory; run locally in the
background, no args). 30-min poll of sacct + mgh sinfo; exits (waking the
agent) on any terminal state or on first RUNNING of a job in `UNVERIFIED`
(currently `13518484 13518485 13518486`). ANSI codes from orc's ssh banner
are stripped before grepping (earlier bug, fixed). mgh-1-1 hop wake is
disabled in-script (see §2). Keep `JOBS`/`UNVERIFIED` in the script current
as jobs start/get replaced. On wake for a new RUNNING arm: verify per §1
recipe (metadata lags job start by ~5–10 min of Julia startup — poll for the
fresh `*.metadata.toml` with `-newermt <job start>`).

Beware stale dirs when verifying: old `p018_*_3r` (with-merge, Aug 28) dirs
sit next to the `_nm` dirs — always check the metadata file's mtime is after
the job's start before trusting it.

## 4. Readout when arms complete (unchanged from pause doc §3)

With-merge λ3.0 column: CT̄ = 0.0708 (NT36) / 0.0718 (NT72) / 0.0738 (NT144).
- merge-off ladder: flattening of CT̄(NT) ⇒ the NT climb was merge error.
  Also check np vs no-merge estimates (~440–520k @ rev 30, saturation
  520–690k, cap 1.5M) and NT144 late-run step cost (~2×).
- fixed-rlxf ladder: NT36 rung = existing completed run 13507289
  (CT̄ 0.0707752); compare NT72/NT144 within-ladder trend.
Then Ryan wants to brainstorm further NT-problem approaches — phase_17
candidate list + these two ladders are the prep.

## 5. Guardrails (unchanged — ask Ryan first)

Notebook + ledger: whole arc approved-pending-draft, still undrafted.
Commits: launcher hooks + band-aid still uncommitted (mirrored to silos).
No FLOWVPM/FastMultipole source edits. 026 on hold. Don't delete 1.5R data.
Keep unused `-lg` LineGauss silos. ssh orc needs `bash -lc`; sacct
authoritative; rsync `--checksum`; exit 124 → `_s2` chain.
