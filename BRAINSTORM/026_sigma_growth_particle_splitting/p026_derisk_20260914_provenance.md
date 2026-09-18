# 026 de-risk wave provenance — campaign/p026-derisk-20260914

Authorization: Ryan 2026-09-14 (three rulings in-session): 3-arm de-risk
first; cold-start (GPU); f_visc enabled at 0.587 (= 4^(1/3)−1, count-matched
tetra4 analog of f_comp = √3−1, picked from AskUserQuestion options);
SIGMA_FLOOR_FRAC=0.1 for the trio (0.25 fallback if the trio fails; keep 0.1
for the remaining 11 arms if it passes).

## Pins (annotated tag `campaign/p026-derisk-20260914`, pushed to origin)

| repo | branch | tag commit |
|---|---|---|
| FLOWVPM.jl | flowpanel | `bf88806` |
| FLOWPanel.jl | fastmultipole | `035f50b` |
| FastMultipole | flowpanel-20260817 | `ac7230a6` |

Worktrees (orc): `~/campaigns/p026-derisk-20260914/{FLOWPanel.jl,FLOWVPM.jl,FastMultipole}`
created from the tag; campaign env `~/campaigns/p026-derisk-20260914/env`
(copy of `~/projects/envs/x86_64` with the three dev-paths repointed at the
campaign worktrees). Worktrees carry no uncommitted state.

## GPU-splitting verification (precondition, 2026-09-14)

- Local: FLOWVPM `runtests_resolution_split.jl` 1055/1055 (incl. t12
  CPU-vs-broadcast parity), FLOWPanel wake/replay suites clean.
- CUDA smoke `scr_p026gpuv_split` (job 13689273, m13h H200, warm-start
  gpu40 s950, 40 steps): **gate PASS** (gpu_gemv=40, cpu_gemv=0, nan=0,
  rc=0); all three mechanisms fired on GPU (viscous/compress/elongate);
  ~7.6 s/step vs CPU twin ~176 s/step (≈23×).
- `scr_p026gpuv_splitmerge` (13689274, MERGE_OVERLAP=3.5): code path ran
  (merge gate active, splits firing); died step ~974 on the euler_exp
  substep-budget guard (dt·|L| ≈ 4.8e3) after an elongate-split storm —
  the ignition continuation blowing up harder under the overlap merge gate.
  Treated as a physics observation (the backstop worked), not a code
  failure; it foreshadows the campaign merge A/B.
- CPU twin `scr_p026cpuv_split` (13689275): COMPLETED all 40 steps.
  **Split-rate parity confirmed**: 34 split lines on both backends; final
  steps CPU elongate 161-199 events/step (480-593 children) vs GPU 172-198
  (511-588 children) — RNG-level agreement. ~176 s/step CPU vs ~7.6 s/step
  GPU (~23x).

## Arms (submit from the campaign FLOWPanel worktree)

Common submission env: `SIGMA_FLOOR_FRAC=0.1 WAKE_SPLIT_VISCOUS=true
WAKE_SPLIT_FRAC_VISCOUS=0.587 WAKE_SPLIT_FRAC_COMPRESS=0.73
WAKE_SPLIT_FRAC_ELONGATE=0.3`; cold start (no RESTART_*); GPU via
`run_p018_screen_gpu052.slurm.sh h200 <case>` with
`P018_REPO_OVERRIDE`/`P018_PROJECT_OVERRIDE` at the campaign worktree/env.

| arm | case | extra env | length |
|---|---|---|---|
| grow cap030 | `scr_p026sp_nt144_cap030` | (linegauss default; `WAKE_SPLIT_SIGMA_MAX=0.030` is in the case def as emission clamp) | `NREVS=20` (through the ~step-2250 cliff region) |
| shrink split | `scr_p026s9_exp_split` | `FLOWPANEL_FILAMENT_REG=vatistas` (exp bracket predates linegauss default) | `NREVS=12` (ignition ~step 210 + margin) |
| merge A/B twin | `scr_p026s9_exp_split_mo35` | `FLOWPANEL_FILAMENT_REG=vatistas` (case def carries `MERGE_OVERLAP=3.5`) | `NREVS=12` |

Acceptance (de-risk): all arms run to completion or die on a *guard* with
interpretable telemetry; cap030 shows no cliff and 6–8 s/step-class cost;
split/skip counters sane; exp pair separates the merge policies; floor 0.1
does not destabilize the healthy phase.

## Data

The dispatcher wipes/moves `data/$RUN_NAME` before a cold run, which kills a
pre-placed symlink, so: runs write into the campaign worktree's `data/<case>/`
during execution (Das arc table symlinked per-file into the worktree), and
each run dir is MOVED to the shared root `~/projects/FLOWPanel.jl/data/`
with a symlink left behind immediately at harvest. Storage/archiver agents
attribute by realpath as usual.

## Submissions (2026-09-14 evening)

| job | arm | pool |
|---|---|---|
| 13691080 | `scr_p026sp_nt144_cap030` (NREVS=20) | m13h H200, 24 h |
| 13691081 | `scr_p026s9_exp_split` (NREVS=12, vatistas) | m13h H200, 8 h |
| 13691082 | `scr_p026s9_exp_split_mo35` (NREVS=12, vatistas) | m13h H200, 8 h |

Submitted from the campaign worktree with the common split/floor env
(§Arms). Banners to be verified at job start (ops rule).

## Extension: cap030 chained restart (2026-09-15, Ryan-approved)

Harvest (see `derisk_harvest_20260915.md`) found all three arms ran the
dispatcher default NREVS=8 (+1 spinup = 9 revs): the dispatcher exports
NREVS unconditionally after sbatch env lands, and the de-risk case defs
carried no NREVS override — the submitted NREVS=20/12 were clobbered.
cap030 therefore stopped at rev 9, short of the historical FMM-adequacy
cliff at rev ~15.3 (steps 2200–2248 @ NT144). Acceptance otherwise PASS →
floor 0.1 kept; wave-2 HELD pending this fix (Ryan 2026-09-15).

Fix: `export NREVS=17` added to the `scr_p026sp_nt144_cap030` case arm
(total 18 revs incl. spinup = 2592 steps; cliff bracket + ~2.4 rev margin)
plus a clobber-warning comment at the dispatcher default. New tag on
FLOWPanel only: **`campaign/p026-derisk-ext-20260915`** (launcher-only
change; FLOWVPM stays `bf88806`, FastMultipole stays `ac7230a6` on the
original tag). Campaign FLOWPanel worktree fast-forwarded to the new tag;
env Manifest dev-paths unchanged (same worktree path).

Run: chained restart, `RESTART_STEP=1295` from the retained state in
`~/projects/FLOWPanel.jl/data/scr_p026sp_nt144_cap030/` (reached through
the worktree's data symlink; restart mode preserves the run dir and appends
to the same VTK series), same common split/floor env as §Arms, m13h H200.
Expected +1297 steps at 10–14 s/step ≈ 4–6 h.

Extension submission: job **13694747** (m13h H200, 2026-09-15), cwd =
campaign FLOWPanel worktree @ `campaign/p026-derisk-ext-20260915` (1b59af5),
`RESTART_STEP=1295` + §Arms common split/floor env, wrapper =
`~/projects/FLOWPanel.jl/examples/run_p018_screen_gpu052.slurm.sh h200
scr_p026sp_nt144_cap030`.

Extension outcome (2026-09-15): job 13694747 COMPLETED, gate_rc=0, steps
1296–2591 (rev 18). **No cliff** — smooth N-scaling 10→18.5 s/step through
the 2200–2248 window; wake health flat; CT cycle-mean 0.07484 ±0.34%.
De-risk acceptance fully PASS. Details: `derisk_harvest_20260915.md`.

## Wave-2 launch (2026-09-15 late — 8-arm linegauss slate)

Rulings on record (Ryan, 2026-09-15 evening): linegauss only (vatistas
arms dropped); three-way mechanism split `_floor` / `_split` / `_fs`;
floor 0.1 wherever armed; telemetry hardening (gate iii) APPROVED and
implemented before launch.

**Pins** (annotated tags, worktrees under `~/campaigns/p026-derisk-20260914/`):

- FLOWPanel `campaign/p026-wave2-20260915` (**4d6d2c0**) — NREVS fixes
  (12 on all s9 arms, 17 on cap018), new `scr_p026s9_explg_fs_mo35` case
  def, wave-2 telemetry (below), 026/021 doc bundle.
- FLOWVPM `campaign/p026-wave2-20260915` (**2b253db**, supersedes bf88806)
  — merge event sink + `SIGMA_FLOOR_HITS` floor-clamp counter.
- FastMultipole unchanged `campaign/p026-derisk-20260914` (**ac7230a6**).

Both worktrees fast-forwarded to the tags; env Manifest dev-paths
unchanged. NOTE: tags pushed to orc; the orc `fastmultipole`/`flowpanel`
BRANCHES were NOT updated (they hold a 2026-09-02 WIP-snapshot commit
absent from local history — flagged to Ryan, tags are the campaign pins).

**Telemetry (gate iii, smoke-verified locally end-to-end):**

- `data/<run>/merge_events.csv` (append mode, restart-safe): one row per
  accepted merge pair, `step,np,sigma_i,sigma_j,dist`. Env
  `MERGE_EVENT_LOG` (default true). No sigma0 columns — no birth-sigma
  storage exists in the particle field; classify vs local shed sigma at
  harvest.
- wake-health CSV: new `mean_sigma,max_sigma,floor_clamp_cum` columns
  appended after all existing/optional columns. `floor_clamp_cum` is the
  cumulative per-process count of sigma_guard floor-clamp engagements
  (FLOWVPM `SIGMA_FLOOR_HITS`, all four euler/euler_exp scalar+broadcast
  paths, live 1:np prefix only); resets to 0 on restart.

**Submissions (2026-09-15, m13h H200, cwd = campaign FLOWPanel worktree,
wrapper `run_p018_screen_gpu052.slurm.sh h200 <case>`, overrides at the
campaign worktree/env):**

Common env: floor arms `SIGMA_FLOOR_FRAC=0.1`; split arms
`WAKE_SPLIT_VISCOUS=true WAKE_SPLIT_FRAC_VISCOUS=0.587
WAKE_SPLIT_FRAC_COMPRESS=0.73 WAKE_SPLIT_FRAC_ELONGATE=0.3`; fs/cap018
arms both sets.

| job | case | mechanism | walltime |
|---|---|---|---|
| 13712311 | `scr_p026s9_ctrllg_floor` | floor only | 8 h |
| 13712312 | `scr_p026s9_ctrllg_split` | split only | 8 h |
| 13712313 | `scr_p026s9_ctrllg_fs` | floor+split | 8 h |
| 13712314 | `scr_p026s9_explg_floor` | floor only | 8 h |
| 13712315 | `scr_p026s9_explg_split` | split only | 8 h |
| 13712316 | `scr_p026s9_explg_fs` | floor+split (lg merge-A) | 8 h |
| 13712317 | `scr_p026s9_explg_fs_mo35` | fs + MERGE_OVERLAP=3.5 (lg merge-B) | 8 h |
| 13712318 | `scr_p026sp_nt144_cap018` | grow-side cap, floor+split env | 24 h |

At submission 13712311/12/13 started immediately; 14–18 PENDING on
QOSMaxCpuPerUser. Banner verification owed on ALL 8 at job start (nrevs
12/17, linegauss, floor/split knobs per slate, merge event log line).

## Wave-2 banner verification + outcomes (2026-09-16 morning)

**Banner verification: PASS on all 8 arms** (hpc-monitor sweep of
`slurm-<jobid>.out` in the campaign worktree, spot-checks inline):
nrevs 12.0/468 steps on the seven s9 arms and 17.0/2592 on cap018;
LineGauss regularization pinned everywhere; `SIGMA_FLOOR_FRAC=0.1
(guard=on)` on 11/13/14/16/17/18 and `0.0 (guard=off)` on the split-only
arms 12/15; split knobs `f_visc=0.587 f_comp=0.73 f_elong=0.3
elongate_overlap=2.4` ACTIVE on 12/13/15/16/17/18 and absent on
floor-only 11/14; "Merge event log:" + "Sigma telemetry:" lines present
on all 8. No mis-arm; nothing cancelled.

`MERGE_OVERLAP=3.5` on 13712317 does not print in the stdout banner and
the settings dump (`*_case_metadata.toml`, driver line ~1784) is written
only AFTER `simulate!` returns, so the early-dead arm has no dump.
Arming confirmed **empirically** from `merge_events.csv`: accepted-pair
max(dist/sigma_min) = 0.2857 = 1/3.5 exactly (hard overlap gate), vs
1.288/1.279 on `explg_fs`/`ctrllg_fs` (absolute-radius gate). VERIFIED.

**Outcomes (all 8 terminal by 2026-09-16 morning; judged by outputs):**

| job | case | end state | last step | max np (merge log) |
|---|---|---|---|---|
| 13712311 | ctrllg_floor | **COMPLETED** 467/467, 3h01 | 467 | 272,234 |
| 13712312 | ctrllg_split | 8 h walltime truncation | 370/467 | 485,793 |
| 13712313 | ctrllg_fs | 8 h walltime truncation | 355/467 | 473,932 |
| 13712314 | explg_floor | DIED: DomainError `dt*|L|`=2133 exceeds euler_exp broadcast substep budget | 274 | 266,180 |
| 13712315 | explg_split | DIED: same DomainError (2229) | 280 | 367,219 |
| 13712316 | explg_fs | DIED: PARTICLE OVERFLOW at 500,000 cap | 295 | 496,691 |
| 13712317 | explg_fs_mo35 | DIED at np=496,445 — .err ends with FMM sigma-adequacy ratio 0.999 warnings then no traceback (abrupt exit, rc=1); overflow-class death at the same cap, traceback lost | 321 | 496,445 |
| 13712318 | cap018 | **COMPLETED** 2591/2591 (17 revs incl. spinup), 8h50 | 2591 | 435,749 |

Notes: zero NaN gate lines, zero WakeGeometryError anywhere. The whole
exp (expint) family runs far hotter in particle count than ctrl at the
same step and died at 59–69%: floor-only and split-only by physical
blow-up (the `dt*|L|` substep-budget guard), both fs arms by exhausting
the 500k particle capacity (mo35's overlap gate merges strictly less, as
predicted, and still reached the cap ~26 steps later than merge-A).
13712317's radix-FMM sigma-adequacy warnings (ratio 0.999, limit
~0.0081) explicitly recommend arming `ResolutionSplitOpts.sigma_max` —
the cap018-style grow-side cap. cap018 itself completed with max np
435,749 and no such deaths. ctrl split arms flirted with the same cap
(486k/474k) and were saved only by the 8 h wall. Slurm sacct states
(FAILED/COMPLETED) recorded but not used as evidence per policy.

Data on /home: ~130 G across the 8 worktree `data/` dirs (cap018 93 G,
s9 arms 5.3–8.8 G each). Harvest/move to shared root next.

## Wave-2 harvest (2026-09-16; data MOVED to shared root + symlinks back)

All 8 run dirs moved `worktree data/` → `~/projects/FLOWPanel.jl/data/`
(symlinks back; cap018 93 G, s9 arms 5.3–8.8 G). Thrust column = CFx
(rotor axis x, CT = −mean CFx) from `*_monitor02_force_system1.csv`; a
first harvester pass mistakenly used CFz and mis-grepped the split
lines — every number below re-verified inline against CFx and the .out.

**Ctrl (non-expint) family blew up mid-run, all three arms** — first
|CFx|>1 at step 259 (fs) / 277 (split) / 294 (floor); ctrllg_floor then
merge-collapsed np 246k→54k (steps 300–350, ~1500 merges/step,
88,914 floor clamps) and "completed" with headline CT 270.8 — unusable
past ~step 250. ctrllg_split σ_max spiked to 0.192. Floor-clamp
engagement on ctrl arms (88k–493k) is a blow-up symptom, not a healthy-
run driver: exp arms logged only 87–6,230 clamps, cap018 zero.

**CT (mean −CFx ± sd), exp family common window steps 200–273:**

| arm | CT | sd |
|---|---|---|
| explg_floor | 0.0767 | 0.0085 |
| explg_split | 0.0766 | 0.0032 |
| explg_fs | 0.0798 | 0.0041 |
| explg_fs_mo35 | 0.0783 | 0.0036 |
| ctrllg_floor (200–290, pre-collapse) | 0.0815 | 0.0069 |

Mechanism attribution (exp): floor-only and split-only agree to 0.1%;
split cuts CT noise 2.6×; fs sits ~+4% above either single mechanism.

**Merge-A/B linegauss re-anchor:** fs (200–295) 0.0787 vs fs_mo35
(200–320) 0.0766 → mo35 −2.6%; the vatistas de-risk pair
(exp_split/_mo35, 200–273) was a null (0.0780 vs 0.0782, +0.2%).
mo35 overlap-gate arming verified (max dist/σ_min = 1/3.5 exactly);
merge rate cut ~3× (279 vs 835 events/step). Caveat: both linegauss fs
arms were racing to the 500k particle cap during the window.

**Cap ladder:** cap018 healthy end-to-end — headline CT 0.07395 ±
0.0703% over final 2 revs (fails only within-rev p-p 0.038 vs tol
0.02); cap030 final-2-rev CT 0.0748 → cap018 −1.1% vs cap030. cap018
final np 435,749, σ_max 0.030, zero floor clamps.

**Viscous splitting fired** in the exp fs arms (max/step: 175 on fs,
452 on mo35) — first observed f_visc firing in this regime.

Open for Ryan: rerun exp s9 arms with σ_max cap and/or larger particle
capacity (adequacy warnings recommend the cap; all four exp arms died
at a ceiling); ctrl-family linegauss instability at steps ~260–290 is a
new unexplained blow-up channel; resubmit ctrl split/fs truncations?

## Merge-σ redesign implementation + smoke (2026-09-17/18)

§22 rulings implemented in FLOWVPM (branch `flowpanel`, UNCOMMITTED,
Ryan-gated):

- **§22.1 second-moment σ** (`FLOWVPM_merging.jl`): σ_new² = ⟨σ²⟩_w +
  (1/3)⟨|xᵢ−x̄|²⟩_w, w=|α|, variance about the actual placement
  (already |α|-weighted; unweighted fallback only when Σ|Γ| ≤ √eps —
  reviewed by Ryan 2026-09-18: fine). Works for n-member clusters.
  Pair gates, strength/position, circulation (uses Σσ), vol (=Σvol,
  Ryan 2026-09-18: leave), merge_events.csv schema all unchanged.
- **§22.2 ledger lineage** (`_rsplit_merge_lineage!`,
  `FLOWVPM_resolution_split.jl`): merge-path reset replaced; σ₀², dvisc,
  drvpm := same |α|-weighted means; separation term → drvpm. **Axis
  lineage (Ryan 2026-09-18: do it)**: axis/weight := sign-aligned
  |α|-weighted means (see §22.2 implementation-review block). Host +
  device dispatch mirrors `_rsplit_reset_slot!`; split-path resets
  untouched.

Tests all green: FLOWVPM merging 78/78 (5 ruled + 1 axis test added),
resolution-split 1058/1058, filament-edge-graph 477/477; FLOWPanel
unit_wake/unit_replay/unit_simulate exit 0. Three old cbrt-law
assertions updated to the new law (merging σ-volume testset,
resolution-split merge-reset testset, filament σ-cube-root testset —
its observed σ matched the new law to 1e-12 pre-update).

**467-step local CPU smoke** (`data/smoke_mergesigma2m_20260917/`,
wave-2 fs env: SIGMA_FLOOR_FRAC=0.1, split 0.587/0.73/0.3 + STRETCH,
MERGE_EVENT_LOG on; NREVS=0.2, 4 threads, exit 0, zero errors; run
PRE-axis-lineage — σ law + ledger only):

| metric | new law | wave-2 baseline |
|---|---|---|
| merge events (467 steps) | 4,763 | 64k (gate-iii telemetry smoke) |
| per-event σ growth vs larger member | median 0.24%, p95 0.52%, max 1.3% | same pairs under cbrt law: median 25.8% |
| max σ / shed σ at end | 1.109 (0.0346 m, drift 0.033→0.035 over steps 116→467) | exp arms died at 2.5–3.4×; cap018 held 0.030 only via cap |
| floor_clamp_cum | 32 | healthy-range |

np peaked ~42k; CT cycle-mean 0.0634 ± 6.4% (smoke window, not a
convergence claim). The σ-pump is gone; §22.3 rerun-slate caps should
be re-derived, not reused at k≈2–3 blindly.
