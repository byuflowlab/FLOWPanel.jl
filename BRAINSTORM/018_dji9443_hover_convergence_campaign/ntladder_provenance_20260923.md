# 018 NT-ladder campaign — provenance (2026-09-23)

Ryan-approved 2026-09-23 (this session): constant-handoff-azimuth NT ladder +
linegauss switch + σ-guard floor sweep + NT144 fine rung, all on the 032
round-2 champion configuration (= A5, `..._omi1_split_mo4_eng`, PASS
2026-09-23: CT 0.07136±5.6e-5, monitor04 max 154.8).

## Pins

| Repo | Pin | Worktree |
|---|---|---|
| FLOWPanel | tag `campaign/p018-ntladder-20260923` (`3d7fba4`, = `af92740` + case-table rows only) | `orc:/home/rander39/campaigns/p018-ntladder-20260923/FLOWPanel.jl` |
| FLOWVPM | `8d4a3b4` (new merge law, production default) | `orc:/home/rander39/campaigns/p026-derisk-20260914/FLOWVPM.jl` |
| FastMultipole | as in p026-derisk env | `orc:/home/rander39/campaigns/p026-derisk-20260914/FastMultipole` |
| Julia env | `orc:/home/rander39/campaigns/p018-ntladder-20260923/env` (copied from p032 env, FLOWPanel dev-path repointed, instantiated with spack julia/1.11.7) | — |

The only commit on top of `af92740` adds launcher case rows
`p018_csarc_{n1_nt18,n2_nt36,n4_nt72,n4_nt144}_l3p0` (clones of the champion
case `p018_csarc_n2_nt72_l3p0` changing NWAKEROWS/NT/P_PER_STEP/RELAX_RLXF).
No Julia source changes vs `af92740`.

## Design rationale (Ryan rulings this session)

- **linegauss**: all runs pin `FLOWPANEL_FILAMENT_REG=linegauss` (exact
  segment-convolved kernel, 052d) replacing the campaign-pinned `gaussian`.
- **Constant handoff azimuth** (N·360/NT fixed): ladder at 20° (NT18/N1,
  NT36/N2, NT72/N4; NT54 skipped per Ryan), champion geometry 10° for
  ANCH (NT72/N2) and NT144/N4. Handoff = age of a panel-wake row at
  panel→particle conversion (one row converts per step once the N-row sheet
  fills; SigmaOverlap shedding keeps filament sampling σ-adequate at coarse NT).
- **Exact-rate relaxation** r(NT)=1−(1−0.3)^(36/NT): 0.51 / 0.3 / 0.16334 /
  0.08539 at NT 18/36/72/144. Rungs score as labeled model-def A/Bs.
- **P_PER_STEP** scaled 1/NT per F1b ladder convention: 24/12/6/3.
- **σ-guard sweep** (Ryan: ceiling off, iterate floor from champion 0.25):
  floors 0.0625/0.125/0.5 vs ANCH (=0.25 point; its ceiling 0.030 kept for the
  exact A5 A/B but verified inert in round 2 — peak max_sigma 0.017–0.018 m).
  Acceptance: no ignition; CT within +0.39%-guard-offset expectations; radial
  circulation distribution overlays ANCH.
- **NT144/N4**: floor-only guard at 0.25, ceiling off, MAX_PARTICLES=3M
  (A5 ended at 757k still growing; NT144 sheds 2× rows/rev), 48 h wall.
- **L18 pre-declared exploratory**: allowed to fail / fall off trend without
  contaminating the NT36→72→144 convergence claim (rlxf 0.51 beyond tested
  range; 27-step spin-up; 18 rows/rev vs tip helix).

## Arm matrix (8 arms × m13h/eng twins = 16 jobs)

Common env (all arms): TRUNCATION_RADIUS_R=3.0, P018_SETTLE_REVS=22,
SFS_RLXF=0.0025031 (DynamicSFS backscatter), PARTICLE_OMIT_ROOT_R_OVER_R=0.12,
MERGE_OVERLAP=4, WAKE_SPLIT_VISCOUS=true, WAKE_SPLIT_STRETCH=true,
WAKE_SPLIT_FRAC_{VISCOUS,COMPRESS,ELONGATE}=0.587/0.73/0.3,
FLOWPANEL_FILAMENT_REG=linegauss, MAX_PARTICLES=1.5M (NT144: 3M),
P018_REPO_OVERRIDE/P018_PROJECT_OVERRIDE = campaign worktree/env above.
Case carries mesh 45_185_ct4, RPM 5400, depth 4R, OVERLAP 2.75,
MERGE_R_FACTOR 0.0055, SIGMA_CHORD_FRACTION 0.313, DAS λ=3.0 arc/steady,
30-rev horizon (28.5 + 1.5 spin-up).

| Arm | Case | N | NT | rlxf | pps | FLOOR_FRAC | CEIL | steps | wall | run name (…_{m13h,eng}) |
|---|---|---|---|---|---|---|---|---|---|---|
| L18 | p018_csarc_n1_nt18_l3p0 | 1 | 18 | 0.51 | 24 | 0.25 | 0.030 | 540 | 6 h | p018_csarc_n1_nt18_l3p0_3r_srlx_g25_omi1_split_mo4_lg |
| L36 | p018_csarc_n2_nt36_l3p0 | 2 | 36 | 0.3 | 12 | 0.25 | 0.030 | 1080 | 8 h | p018_csarc_n2_nt36_l3p0_3r_srlx_g25_omi1_split_mo4_lg |
| L72 | p018_csarc_n4_nt72_l3p0 | 4 | 72 | 0.16334 | 6 | 0.25 | 0.030 | 2160 | 12 h | p018_csarc_n4_nt72_l3p0_3r_srlx_g25_omi1_split_mo4_lg |
| ANCH | p018_csarc_n2_nt72_l3p0 | 2 | 72 | 0.16334 | 6 | 0.25 | 0.030 | 2160 | 12 h | p018_csarc_n2_nt72_l3p0_3r_srlx_g25_omi1_split_mo4_lg |
| G-lo2 | p018_csarc_n2_nt72_l3p0 | 2 | 72 | 0.16334 | 6 | 0.0625 | off | 2160 | 12 h | p018_csarc_n2_nt72_l3p0_3r_srlx_gf0625nc_omi1_split_mo4_lg |
| G-lo1 | p018_csarc_n2_nt72_l3p0 | 2 | 72 | 0.16334 | 6 | 0.125 | off | 2160 | 12 h | p018_csarc_n2_nt72_l3p0_3r_srlx_gf125nc_omi1_split_mo4_lg |
| G-hi | p018_csarc_n2_nt72_l3p0 | 2 | 72 | 0.16334 | 6 | 0.5 | off | 2160 | 12 h | p018_csarc_n2_nt72_l3p0_3r_srlx_gf50nc_omi1_split_mo4_lg |
| NT144 | p018_csarc_n4_nt144_l3p0 | 4 | 144 | 0.08539 | 3 | 0.25 | off | 4320 | 48 h | p018_csarc_n4_nt144_l3p0_3r_srlx_g25nc_omi1_split_mo4_lg |

Twin scheme (round-2 pattern): every arm submitted to m13h (qos gpu) and eng
(qos eng), both `--constraint=intel --gres=gpu:h200:1 --cpus-per-task=64
--mem=192G --no-requeue`, distinct `_m13h`/`_eng` run names; pair watcher
cancels the loser when the other twin RUNs (loser cancellation pre-authorized).
m13h submitted in table order, eng in reverse; **NT144 twins submitted last on
both partitions with `--dependency=afterany:<all 14 other jobids>`** so NT144
runs only after every other arm is terminal (Ryan: "NT144 last in every case").

Job names: `fp-018gpu-ntl-{l18,l36,l72,anch,glo2,glo1,ghi,n144}{m,e}`; logs
`<wt>/logs/slurm/slurm-<jobname>-<jobid>.{out,err}`.

## Acceptance

Per rung: complete all steps, finite CT, bounded monitor04 Γ/σ², gate_rc=0;
banner must show linegauss, correct N/NT/rlxf/pps, guard floor/ceil per table,
omission 1/41 @0.12R, split+stretch, mesh ct4. Ladder read: CT slope over
NT36→72→144 at fixed model-def labels; ANCH−A5 = linegauss effect at champion
config; L72−ANCH = N-effect (2→4) at NT72; guard arms vs ANCH = floor
dose-response (CT + circulation distribution).

## Job IDs (submitted 2026-09-23 ~22:10)

| Arm | m13h | eng |
|---|---|---|
| L18 | 13878867 | 13878880 |
| L36 | 13878868 | 13878879 |
| L72 | 13878869 | 13878878 |
| ANCH | 13878870 | 13878877 |
| G-lo2 | 13878871 | 13878876 |
| G-lo1 | 13878872 | 13878875 |
| G-hi | 13878873 | 13878874 |
| NT144 | 13878881 (dep) | 13878882 (dep) |

NT144 twins `--dependency=afterany:` all 14 other job IDs. Pair watcher
`/home/rander39/campaigns/p018-ntladder-20260923/ntl_pair_watcher.sh`
detached on orc login (PID 2415612 at launch — verify via
`pgrep -af ntl_pair_watcher`, not the PID), log
`<wt>/logs/slurm/ntl_pair_watcher.log`. Loser cancellation pre-authorized;
nothing else is.

## Round-1 FAILURE + r2 resubmission (2026-09-23 ~22:30–23:00)

**All 7 non-NT144 arms died at init**: `DAS_ARC_TABLE=data/
p018_cs_l3p4_rs1_te_downwash_te.csv` is an UNTRACKED asset (a symlink in the
p032 worktree → `~/projects/FLOWPanel.jl/data/`); the fresh worktree lacked
it, so every started eng twin failed the driver's isfile guard (gate_rc=1,
"GPU source-influence path never ran") ~1 min after start; the watcher had
already cancelled the m13h twins as designed. **Lesson: fresh campaign
worktrees need untracked data assets staged — check every relative-path env
value in the case rows before submitting.** Fix: CSV copied into
`<wt>/data/` (md5 08375291ed2b542ea946d09730e0b629, verified).

**NT144-eng 13878882 survived and is RUNNING (started 22:41:27)**: its
afterany deps cleared as the failures terminated, it started before a hold
landed, and the CSV copy beat its isfile check (banner clean: linegauss,
N=4, NT144, pps 3, no ERRORs). It therefore runs FIRST, not last —
Ryan-gated decision whether to keep or kill it. 13878881 (n144m)
watcher-cancelled.

**r2 resubmission** (same env/case per arm, run names `..._lg_r2_{m13h,eng}`,
job names `fp-018gpu-ntl2-*`, submitted m13h forward / eng reverse ~22:58):

| Arm | m13h | eng |
|---|---|---|
| L18 | 13879081 | 13879094 |
| L36 | 13879082 | 13879093 |
| L72 | 13879083 | 13879092 |
| ANCH | 13879084 | 13879091 |
| G-lo2 | 13879085 | 13879090 |
| G-lo1 | 13879086 | 13879089 |
| G-hi | 13879087 | 13879088 |

Watcher `ntl2_pair_watcher.sh` (PID 2570978 at launch), log
`<wt>/logs/slurm/ntl2_pair_watcher.log`. The failed round-1 runs left
partial `data/<runname>_eng/` dirs (merge_events.csv only, no monitors) —
ignore/clean at harvest.

## Results (harvest as arms finish, 2026-09-24)

Twin winners so far: all **m13h** for L18/L36/G-hi (eng twins
watcher-cancelled); NT144 winner = **eng** 13878882 (ran first, out of
order — **Ryan ruled KEEP 2026-09-24**, out-of-order execution accepted).

| Arm | Job | Verdict | Steps | gate_rc | CT cycle-mean | mon04 max / end | particles peak / final | min_sigma | floor_clamp_cum |
|---|---|---|---|---|---|---|---|---|---|
| L18 (expl.) | 13879081 m13h | **PASS** | 540/540 | 0 | 0.068488 ± 5.65e-5 | 109.4 / 82.1 | 295,153 / 295,153 | 0.001197 | 13,471 |
| L36 | 13879082 m13h | **PASS** | 1080/1080 | 0 | 0.070457 ± 6.86e-5 | 112.0 / 71.1 | 330,562 / 330,562 | 0.001194 | 247,922 |
| L72 | 13879083 m13h | **PASS** | 2160/2160 | 0 | 0.070656 ± 8.32e-5 | 86.2 / 86.2 | 729,022 / 729,022 | 0.00119 | 4,353,691 |
| G-hi | 13879087 m13h | **PASS** | 2160/2160 | 0 | 0.0714012 ± 7.21e-5 | 66.8 / 66.8 | 675,801 / 675,801 | 0.00238 | 74,255,077 |
| G-lo1 | 13879089 eng | **PASS** | 2160/2160 | 0 | 0.0713961 ± 1.24e-4 | 92.2 / 92.2 | 871,041 / 871,041 | 0.000595 | 11,760 |
| NT144 | 13878882 eng | **PASS** | 4320/4320 | 0 | 0.0733267 ± 1.51e-4 | 93.5 / 93.5 | 766,435 / 766,435 | 0.00119 | 2,764,836 |
| G-lo2 | 13879090 eng | **PASS** | 2160/2160 | 0 | 0.0714079 ± 1.08e-4 | 164.1 / 78.0 | 986,443 / 986,443 | 0.000381 | 0 |
| ANCH | 13879091 eng | **PASS** | 2160/2160 | 0 | 0.0714104 ± 8.58e-5 | 111.9 / 76.4 | 878,034 / 878,034 | 0.00119 | 7,114,406 |

Banner values verified against the arm matrix for both (linegauss,
NT/N/rlxf/pps, guard 0.25/0.030, omission 1/41 @0.12R, split+stretch,
45_185_ct4, H200, 64 threads, pinned paths). mon04 bounded (A5 reference
max 154.8). Particle counts still growing at end-of-run (peak=final), well
under 1.5M cap. floor_clamp_cum L36/L18 ≈ 18.4× (more rows + lower rlxf).
Early ladder read: L18→L36 CT 0.0685→0.0705, toward A5 anchor 0.07136.

Harvest 2 (2026-09-24 ~PM): L72/G-hi/G-lo1/NT144 rows added above. Banner
verification PASS for all four against the arm matrix (verbatim banners
pulled: linegauss pinned, correct N/NT/rlxf/pps, SIGMA_FLOOR_FRAC and
SIGMA_CEIL per table — L72 0.25/0.030, G-hi 0.5/Inf, G-lo1 0.125/Inf,
NT144 0.25/Inf — omission 1/41 @0.12R, split+stretch, ct4, 64 threads,
SFS=DynamicSFS(rlxf=0.0025031, maxC=1.0, alpha=0.999,
clippings=backscatter, controls=none, nostatic=false)). An earlier
monitor pass mis-parsed sigma_floor=0; verbatim banners refute it.
Twin outcomes: L72/G-hi winners m13h; G-lo1/G-lo2 winners eng
(13879086/85/92 watcher-cancelled as designed). Caveat: harvest-2 mon04
max==end for all four rows (monitor still at its running max at end of
run, unlike L18/L36 which decayed from peak) — all values well under the
A5 reference max 154.8, but worth a glance at the mon04 traces.
Driver's `CONVERGED (Phase 2e criterion)`: true for L72 and G-lo1, false
for NT144 (p2p 0.0255 vs 0.02 tol) and G-hi (0.0226) — informational,
not an acceptance gate.

Key reads so far: ladder CT 0.0685 (NT18) → 0.070457 (NT36) → 0.070656
(NT72) → 0.073327 (NT144); NT36→72 nearly flat (+0.28%), NT72→144 +3.8%
(NOTE ceil confound: L18/L36/L72 run ceil 0.030, NT144 ceil off).
Floor dose at NT72/N=2: G-hi (floor 0.5) 0.0714012 vs G-lo1 (0.125)
0.0713961 — Δ 0.007% despite floor_clamp_cum 74.3M vs 11.8k: CT is
floor-insensitive over this range. ANCH pending for N-effect (L72−ANCH)
and linegauss effect (ANCH−A5 0.07136).

Still live at harvest 2: G-lo2 13879090 eng RUNNING (banner PASS, floor
0.0625/ceil Inf, ~4.2 s/step early); ANCH twins 13879084/91 PENDING.

Harvest 3 (2026-09-25 AM) — **CAMPAIGN TERMINAL, 8/8 arms PASS**. G-lo2
13879090 + ANCH 13879091 (both eng winners; ANCH m13h 13879084 manually
killed by Ryan 2026-09-24, acknowledged) finished 2160/2160, gate_rc=0,
clean .err. Verbatim banners verified against the arm matrix: ANCH
N=2 NT=72 rlxf=0.16334 pps=6 floor 0.25 (0.0011899 m) / ceil 0.030,
linegauss, omission 1/41 @0.12R, split+stretch (clamp=[0.0011899,0.03]),
ct4, SFS=DynamicSFS(rlxf=0.0025031, maxC=1.0, alpha=0.999,
clippings=backscatter, controls=none, nostatic=false); G-lo2 identical
except floor 0.0625 (0.0002975 m) / ceil Inf (clamp=[0.0002975,off]).
mon04 column identified this harvest as `max_gamma_over_sigma2` (col 7 of
monitor04_wake_health_system1.csv); all prior rows consistent. **G-lo2
mon04 peaked at 164.1 — first arm to EXCEED the A5 reference max 154.8 —
but decayed to 78.0 by end (bounded, no ignition; CT unaffected).**
G-lo2 floor_clamp_cum=0 (floor 0.0625 never binds; min_sigma 0.000381 >
floor 0.0002975). Driver CONVERGED (Phase 2e): G-lo2 true (p2p 0.0133),
ANCH false (p2p 0.0243 vs 0.02) — informational, same as NT144/G-hi.
Pair watcher exited cleanly ("all pairs resolved", 2026-09-25 05:03).
Non-critical warns only (FastMultipole constant redefinition, FMM sigma
adequacy ratio→1.0 late in run — expected at these particle counts).

**Key reads (ladder complete):**

1. **Linegauss effect ≈ NULL**: ANCH 0.0714104 ± 8.6e-5 vs A5 (legacy
   filament reg) 0.07136 ± 5.6e-5 → Δ +5.0e-5 (+0.071%, ~0.5σ combined).
2. **N-effect (pure N: 2→4 at equal pps=6) at NT72**: L72 0.070656 vs ANCH
   0.0714104 → −7.5e-4 (−1.06%), floors/ceils matched (0.25/0.030). The
   largest surviving lever besides NT itself.
3. **Floor dose-response FLAT at NT72/N=2**: floors 0.0625/0.125/0.25/0.5
   (G-lo2/G-lo1/ANCH/G-hi) → CT 0.0714079/0.0713961/0.0714104/0.0714012;
   full spread 0.02% while floor_clamp_cum spans 0 → 74.3M. ANCH (ceil
   0.030) sits inside the ceil-Inf band → **ceil also null at NT72/N=2**,
   which weakens (but does not eliminate — N=4 untested) the ceil-confound
   caveat on the NT72→144 +3.8% jump.
4. **CT vs NT (labeled model-def rungs: N=1/2/4/4, pps=24/12/6/3, ceil
   0.030 except NT144 off)**: 0.068488 (18)
   → 0.070457 (36) → 0.070656 (72) → 0.073327 (144); +2.9%, +0.28%,
   +3.8%. Non-monotone convergence pattern: near-flat 36→72 then a jump
   72→144 — not a clean Richardson ladder; NT144 also failed the p2p
   criterion (0.0255), so its cycle-mean carries more within-rev
   variation.

## Round 3 STAGED (Ryan order 2026-09-25): SFS-rlxf correction + N=8 rung

**Motivation (Ryan finding 2026-09-25):** DynamicSFS `rlxf` = Δt/T (FLOWVPM
docstring, `FLOWVPM_subfilterscale.jl:714`) — the Lagrangian-average
relaxation of the dynamic coefficient C. The r1/r2 ladder froze
SFS_RLXF=0.0025031 (itself the exact-rate NT72 conversion of the driver
default 0.005: 1−(1−0.005)^{1/2}) across all rungs, so the physical
averaging time T varied 4× down the ladder (NT18 ~22 revs … NT144 ~2.8
revs vs the intended ~5.6). This is an unflagged confound on the NT axis
(wake rlxf got the exact-rate treatment; SFS rlxf did not). alpha=0.999 is
the filter-level ratio α_τ (spatial), correctly Δt-independent.

**r3 arm matrix (4 arms × m13h/eng twins = 8 jobs).** SFS rlxf per the
exact-rate law from base 0.005@NT36 (same derivation as the r2 value):

| Arm | Case | N | NT | pps | SFS_RLXF | floor/ceil | max_p | wall | run name (…_{m13h,eng}) |
|---|---|---|---|---|---|---|---|---|---|
| S18 | p018_csarc_n1_nt18_l3p0 | 1 | 18 | 24 | 0.009975 | 0.25/0.030 | 1.5M | 6 h | …n1_nt18…_lg_sfsr3 |
| S36 | p018_csarc_n2_nt36_l3p0 | 2 | 36 | 12 | 0.005 | 0.25/0.030 | 1.5M | 8 h | …n2_nt36…_lg_sfsr3 |
| S144 | p018_csarc_n4_nt144_l3p0 | 4 | 144 | 3 | 0.0012524 | 0.25/off | 3M | 48 h | …n4_nt144…g25nc…_lg_sfsr3 |
| N8-144 | p018_csarc_n8_nt144_l3p0 (NEW case) | 8 | 144 | 3 | 0.0012524 | 0.25/off | 3M | 48 h | …n8_nt144…g25nc…_lg_sfsr3 |

All other env identical to r2 (common-env block above; SIGMA_FLOOR_FRAC=0.25,
linegauss, omission, split+stretch, worktree/env overrides). N8-144 extends
the N-prop-NT progression (1@18, 2@36, 4@72 → 8@144); NT×pps=432 held.
Key r3 reads: S144−NT144(r2) = SFS-timescale effect at NT144;
S18/S36 vs L18/L36 = same at coarse rungs; N8-144−S144 = N-effect at NT144
(matched corrected SFS); corrected-SFS ladder slope S18→S36→L72→S144
(L72's frozen 0.0025031 IS the corrected NT72 value, so L72 is reusable).

**Staged assets:** submit script
`/home/rander39/campaigns/p018-ntladder-20260923/p018_submit_ntl3.sh`
(syntax-checked; twins m13h qos gpu / eng qos eng, intel+H200, 64 cpu,
192G, no-requeue; NT144-last via afterany on the four short-arm jobids;
generates+launches `ntl3_pair_watcher.sh` from the ntl2 watcher). The
script self-guards: refuses to run unless the DAS arc table AND the n8
case row are present in the worktree.

**SUBMITTED 2026-09-26 (Ryan approved; N8-144 wall raised to 72 h,
S144 kept at 48 h):** n8 case row committed in the campaign worktree as
**6e9a815** on top of 3d7fba4, tag **campaign/p018-ntladder-r3-20260925**
(annotated; not pushed — pushes Ryan-gated). Jobs:

| Arm | m13h | eng |
|---|---|---|
| S18 | 13898944 | 13898945 |
| S36 | 13898946 | 13898947 |
| S144 | 13898948 | 13898949 |
| N8-144 | 13898950 | 13898951 |

NT144-last enforced: 13898948–51 carry `afterany` on all four short-arm
jobids. Pair watcher `ntl3_pair_watcher.sh` detached on orc login (launch
PID 2629745 — verify via `pgrep -af ntl3_pair_watcher`, not the PID), log
`<wt>/logs/slurm/ntl3_pair_watcher.log`; loser cancellation pre-authorized,
nothing else is.

Archiver (032 round-2, detached 09-24): COMPLETED, VERIFY_FAIL_COUNT=0,
TOTAL_FREED_MB=283,929 (~284 GB), DF_AFTER=467G — still above the 400 G
cap (expected ~377 G; ~90 G attributed to live NT-ladder run growth).
Another storage pass will be needed once the ladder finishes (apply is
Ryan-gated).

## Round 3 results (harvest as arms finish, 2026-09-26)

Babysit events 09-26: watcher cancelled eng losers 13898945 (08:38) and
13898947 (10:01) after their m13h twins started, then the watcher process
died (log intact, process gone). Relaunched detached on orc login
~mid-day 09-26 as PID 1693433 (verified via pgrep) to cover the two NT144
pairs; script is idempotent (done pairs just log "both gone -> done").
S18/S36 m13h banners verified against the r3 matrix (linegauss, exact-rate
wake rlxf, per-arm SFS_RLXF, floor 0.25/ceil 0.030 on, omission 1/41
@0.12R, split+stretch, ct4, correct N/NT/pps).

Same columns/conventions as r1/r2 §Results (CT = monitor02 −CFx cycle
mean; mon04 = col 7 `max_gamma_over_sigma2`):

| Arm | Job | Verdict | Steps | gate_rc | CT cycle-mean | mon04 max / end | particles peak / final | min_sigma | floor_clamp_cum |
|---|---|---|---|---|---|---|---|---|---|
| S18 | 13898944 m13h | **PASS** | 540/540 | 0 | 0.069629 ± 3.94e-5 | 91.5 / 76.5 | 253,349 / 253,349 | 0.001197 | 3,274 |
| S36 | 13898946 m13h | **PASS** | 1080/1080 | 0 | 0.070816 ± 8.04e-5 | 107.2 / 84.9 | 477,501 / 477,501 | 0.00119 | 552,530 |
| S144 | 13898948 m13h | **PASS** | 4320/4320 | 0 | 0.072402 ± 8.67e-5 | 168.7 / 82.2 | 1,069,455 / 1,069,455 | 0.00119 | 11,614,098 |
| N8-144 | 13898951 eng | **PASS** | 4320/4320 | 0 | 0.071396 ± 9.40e-5 | 126.6 / 98.3 | 768,685 / 768,685 | 0.00119 | 6,780,726 |

S18 auxiliary reads: wall 1 h 23 m (08:38–10:01); p2p CT ripple last rev
0.43% (< 2% Phase 2e criterion, CONVERGED=true); mon04 peaked mid-run and
decayed (no ignition); floor clamping minimal (3,274 vs L18's 13,471).

**Key read — S18−L18 = +1.67%** (0.069629 vs 0.068488): correcting the
frozen SFS averaging time at NT18 (rlxf 0.009975 vs 0.0025031, i.e. T
~22 revs → ~5.6 revs) raises CT by +1.67%, moving the coarsest rung
*toward* the NT72/NT144 values — the SFS-timescale confound flattens the
apparent ladder slope at the coarse end. Slope conclusion awaits S36/S144.

S36 (harvested 09-26 ~PM): wall 2 h 59 m (10:01–13:00); p2p CT ripple
last rev 1.16% (< 2%, CONVERGED=true); mon04 decayed 107.2→84.9 (no
ignition). **S36−L36 = +0.51%** (0.070816 vs 0.070457) — same sign as
S18−L18, smaller as expected (NT36 frozen-T error was 2× vs NT18's 4×).
NOTE: harvester subagents misreported both r3 deltas (S18 "+0.076%",
S36 "+1.04%") and the S36 ripple ("0.0116%"); values here recomputed
inline — always recompute subagent-derived percentages.

Corrected-SFS ladder so far (S18, S36, L72 reusable): 0.069629 →
0.070816 → 0.070656, i.e. NT18→36 +1.70%, NT36→72 −0.23% — vs the frozen
ladder's +2.9% / +0.28%. The corrected coarse-end slope is roughly half
the frozen one and NT36→72 is now flat-to-slightly-negative; the
NT72→144 rung (S144, running) decides whether the ladder has converged.

**Dose–response + S144 prediction (09-26).** rlxf = Δt/T so the frozen
0.0025031 implies T = 1/(NT·rlxf) revs: 22.2 / 11.1 / 5.6 / 2.8 at
NT18/36/72/144 vs the intended 5.6 everywhere. The measured corrections
scale with the T-distortion being removed (4× → +1.67%, 2× → +0.51%,
1× → 0 by construction) — a monotone dose–response supporting a real
SFS-timescale effect over noise. At NT144 the frozen T errs in the
OPPOSITE direction (2.8 revs, 2× too short), so the mechanism predicts
**S144 < NT144(r2) = 0.0733**, pulling the ladder top toward L72's
0.0707. S144 landing below 0.0733 confirms the mechanism; the size of
the drop sets the corrected NT72→144 slope (ladder-convergence verdict).

**S144 + N8-144 harvested (2026-09-28; rows above, CORRECTED same day —
see convention note).** Both arms completed 4320/4320, gate
gpu_gemv=4320 / cpu_gemv=0 / 0 NaN / rc=0, no .err fatals; judged by
outputs. NOTE: the harvester subagent's first CT pass was again wrong
(bad awk window + broken std → mean 0.07176, phantom "min=0" rows); an
inline recompute then used the WRONG WINDOW (final rev only, 0.072267 /
0.071486). The r1/r2 and S18/S36 rows all equal the driver's
`case_metadata.toml` `CT_window_mean` (revs 20–30 window; ± =
`CT_cycle_std`), so the rows above now use that same convention:
**S144 = 0.072402 ± 8.67e-5, N8-144 = 0.071396 ± 9.40e-5** (window p2p
`CT_ptp_rel`: 3.11% / 2.57%). Always harvest CT from case_metadata, not
hand-windowed monitor02. Auxiliary:

- **S144**: wall 17 h 07 m (09-26 13:00:32 → 09-27 06:07:39), node
  m13h-1-2. Window p2p ripple 3.11% — FAILS the 2% Phase 2e criterion,
  slightly worse than r2 NT144's 2.55% (the fine rung still resolves
  ≥2% within-rev unsteadiness; informational, not a gate). mon04 peaked
  168.7 mid-run and decayed to 82.2 (no ignition). Particles 1,069,455
  peak=final (under 3M cap; +40% vs r2 NT144's 766k — corrected SFS
  retains more wake). floor_clamp_cum 11.6M (r2 NT144: 2.76M).
- **N8-144**: wall 14 h 44 m (09-26 13:31:28 → 09-27 04:15:17), node
  eng-1-1. Window p2p 2.57% (fails 2%, like every NT144 arm). mon04
  126.6 → 98.3 (no ignition). Particles 768,685 peak=final (−28% vs
  S144). floor_clamp_cum 6.78M.

**Campaign-deciding reads (recomputed inline 2026-09-28, window
convention):**

1. **Prediction CONFIRMED: S144 = 0.072402 < 0.0733.** S144−NT144(r2) =
   (0.072402−0.0733267)/0.0733267 = **−1.26%** — the SFS-timescale
   correction lowers CT at NT144, opposite sign vs the coarse rungs
   (+1.67% @NT18, +0.51% @NT36), as the T-distortion mechanism
   predicts. Sign confirmed at all three off-reference rungs; magnitude
   NOT symmetric with NT36 (both 2× distortions: +0.51% vs −1.26%) —
   plausibly amplified by the finer wake at NT144 (+40% particles).
   Caveat: this read is cross-partition (S144 m13h vs r2 NT144 eng);
   twin equivalence was never measured (losers cancelled).
2. **Corrected-SFS ladders — TWO consistent readings (framing fixed
   2026-09-29 after Ryan pushback; the earlier "N is a confound" wording
   applied only to the mixed S144 ladder):**
   - **Constant-handoff ladder (N ∝ NT, the campaign's design intent;
     N/NT = 0.056 revs at every rung): S18 → S36 → L72 → N8-144** =
     0.069629 → 0.070816 → 0.070656 → 0.071396, i.e. **+1.70% / −0.23% /
     +1.05%**. Here N is part of the scaling law, not a confound.
   - **Fixed-N ladder (pure NT at matched N): S36 → ANCH (N=2) =
     +0.84%; L72 → S144 (N=4) = +2.47%.**
   - The as-run S18/S36/L72/**S144** sequence (+1.70/−0.23/+2.47) mixes
     the two policies (N∝NT for two rungs, then N frozen at 4) — do not
     quote it as "the" ladder.
   - Verdict: NEITHER ladder closes (top rungs +1.05% and +2.47%), and
     the two policies DIVERGE with refinement — endpoints 0.071396 vs
     0.072402 differ by 1.39% = the N-effect, which GREW from NT72 to
     NT144 (read 3). Constant-handoff is the flatter, better-behaved
     family. (Frozen-SFS as-run ladder for reference: +2.9/+0.28/+3.78.)
     All corrected rungs m13h except N8-144 (eng).
3. **N-effect at NT144 (matched corrected SFS): N8-144−S144 =
   (0.071396−0.072402)/0.072402 = −1.39%** vs −1.06% at NT72 (L72−ANCH,
   N 2→4). Same sign, ~30% LARGER at the finer rung — the N axis is not
   converged either, and the growth means fixed-N and N∝NT ladders do
   not meet at these resolutions (cross-partition caveat: N8-144 ran on
   eng, S144 on m13h).

**Settings audit (Ryan request 2026-09-28): banner-diff of L36/L72/
NT144(r2)/S36/S144/N8-144 logs + case_metadata + mon04.** Identical
across all rungs: mesh 45_185_ct4, RPM, formulation, depth 4R, settle
22, Das λ3.0 arc table/src, overlap 2.75, merge_r 0.0055, conversion
legacy (conv_overlap 1.3), attribution upstream, omission 1/41 @0.12R,
split params (f_visc 0.587 / f_comp 0.73 / f_elong 0.3 / mo4 / every 1),
linegauss, FMM knobs (body 17/0.7/109, wake 16/0.6/38), σ floor
0.00119 m, SFS alpha/maxC/clippings, integrator, relax_filter offR,
sigma_chord 0.313. Differences beyond the intended NT/Δt, exact-rate
wake rlxf, SFS rlxf ∝1/NT, pps ∝1/NT:

1. **N: the as-run S-ladder is 1/2/4/4** — N∝NT for two rungs, then
   frozen at the top. Not a stray setting: with N8-144 the data supports
   both a clean N∝NT (constant-handoff) ladder and clean fixed-N pairs;
   only the mixed S144 sequence should not be read as a ladder (read 2).
2. **σ ceiling 0.030 at NT18/36/72 vs OFF at NT144** (guard AND split
   clamp) — **empirically moot**: run-max max_sigma ≈ 0.0174 in both
   L72 and S144 (and 0.015 in N8-144), never within 1.7× of the cap;
   the ceiling never binds anywhere on the ladder.
3. MAX_PARTICLES 1.5M vs 3M — never binds (peaks 477k/729k/1.07M).
4. Partition: corrected ladder all-m13h (clean); only the S144−NT144(r2)
   and N8-144−S144 reads are cross-partition (flagged above).

No other inconsistency found; the ladder's non-convergence is not a
settings artifact.
