# 018 mechanism-isolation wave 1 — campaign provenance (2026-09-15)

Follow-on to `gamma_distribution_status_20260915.md` (far-wake, substantially
temporal Σ|Γ| divergence) and `rlxfscaled_provenance_20260915.md` (g25 guarded
reference pair, +1.27%/doubling). Matrix approved by Ryan (chat, 2026-09-15;
see `mechanism_tests_reset_prompt_20260915.md`): T1 (extend both g25 rungs to
rev 60 by RESTART, gated on a restart-fidelity check), T2 (NT72 merge cadence
rate-matched to NT36), T4 (ConstantSFS both rungs), T5 (Inviscid both rungs).
T3 (merging off) explicitly deferred.

## Pins (annotated tags `campaign/p018-mech-tests-20260915`, cluster repos)

| repo | worktree (cluster) | commit | note |
|---|---|---|---|
| FLOWPanel.jl | `~/wt018/FLOWPanel-mechtests` | `5cdf058` | base `d5dd772` (= expguard pin) + examples commit `5d8711c` + warmstart skew fix (below) |
| FLOWVPM.jl | `~/wt018/FLOWVPM-expguard` | `7468712` | unchanged; alias tag on the expguard pin |
| FastMultipole | `~/wt018/FastMultipole-expguard` | `3da58a1a` | unchanged; alias tag on the expguard pin |

FLOWPanel branch `p018-mech-tests-20260915`, two commits (worktree verified
clean after each; the tag was MOVED from `5d8711c` to `5cdf058` after the
incident below — nothing ran successfully under the first tag placement):

Commit `5cdf058` — **warmstart version-skew fix (src change), forced by the
first gate attempt**. Job 13712091 (first T1 chain) FAILED in
`simulate_warmstart!`: `src/FLOWPanel_warmstart.jl` was written against a
newer FLOWVPM whose `SplittingState` carries `dsigma2_visc`/`dsigma2_rvpm`;
the pinned FLOWVPM `7468712` lacks both fields and its VTPs carry only the
four `split_{sigma_0,H_chi,hold,cooldown}` arrays (verified in the g25 VTP
headers), so both `_clear_splitting_state!` (crash site, :374) and the
all-or-nothing `split_*` restore block (:337–349, would have thrown the
partial-set `ArgumentError` next) were guarded with
`hasproperty(st, :dsigma2_visc)`. Cold-start behavior is untouched (the g25
production pair never entered this path — the skew was latent). No other
skew-prone accesses found (FilamentEdgeGraph fields all match the pin).
Slurm auto-cancelled the first T1b (13712092, dependency-never-satisfied);
both were resubmitted (see Submitted). The T1 restart-fidelity gate remains
the empirical arbiter of restart validity per the reset prompt.

Commit `5d8711c` (examples-only):

1. `examples/rotor_hover_pressure_comparison.jl:778` — merge cadence exposed:
   `every = merge_particles ? parse(Int, get(ENV,"MERGE_EVERY","1")) : 0`
   (default 1 = behavior identical to the pin).
2. `examples/run_dji9443_hover_ct_hpc.slurm.sh` — banner now prints
   `merge_every:${MERGE_EVERY:-1}`; two new case arms `p018_csarc_l3p0_nv` and
   `p018_csarc_n2_nt72_l3p0_nv`, identical to their parents except
   `CORE_SPREADING_ACTIVE=false` (new arms per the standing ops rule: never
   env-override a knob a case arm exports unconditionally).

The worktree's `data/` is a symlink to `~/projects/FLOWPanel.jl/data` (same as
the expguard worktree), so run outputs land in the standard data root and the
untracked DAS arc table `data/p018_cs_l3p4_rs1_te_downwash_te.csv` resolves.

## Julia environment

`~/p018wtenv-expguard-gh200` reused (depot `~/fm052depot-gh200`, julia 1.11.7
aarch64; stack proven by smoke 13610775 and production 13704962/13704963);
after the warmstart fix its Manifest dev-path for FLOWPanel was re-pointed
from `~/wt018/FLOWPanel-expguard` to `~/wt018/FLOWPanel-mechtests` (manual
path edit, same UUID/version — all three dev-paths now at campaign
worktrees). Timing nuance: T2 (13712093) started while the dev-path still
read the expguard worktree; immaterial for T2 — the two worktrees' `src/`
differ only in the warmstart module, which cold starts never load-path
through, and T2's driver/case script always ran from the mechtests worktree
via `P018_REPO_OVERRIDE`.

## Restart mechanism (T1)

Driver env knobs `RESTART_STEP/RESTART_NAME/RESTART_PATH`
(`examples/rotor_hover_pressure_comparison.jl:1368-1370` →
`pnl.simulate_warmstart!`, `src/FLOWPanel_warmstart.jl`): reconstructs from
`<name>_body1.<S>.vtu` + `<name>_wake1.{1,2}.<S>.vts` +
`<name>_wake1_particles.<S>.vtp`, then continues at `start_step=S+1`. Complete
restart quadruplets verified present at steps 1044 and 1079 (NT36) and 2159
(NT72). The metadata flag `restart_reconstruct_required=true` was audited: it
is a generic TOML-serialization marker (`src/FLOWPanel_metadata.jl:15`), never
read back — no bearing on restart validity. Fidelity is instead gated
empirically (below).

**Restart-fidelity gate** (`p018_csarc_l3p0_3r_g25_rgate`): restart NT36 from
step 1044, rerun to 1080; scripted check
(`orc:~/p018_mech_tests_20260915/gate_check.py`) compares the per-step force
monitor (thrust component) against the original over the overlapping steps,
skipping the first 3 after the seam; PASS = max relative ΔCF ≤ 2e-4
(reset-prompt tolerance ΔCT ≲ 1e-4 rel, doubled for the seam transient), plus
an informational particle-count / Σ|Γ| comparison at step 1079. The T1 long
restarts run ONLY downstream of a PASS: T1a is staged after the check inside
the same chained job; T1b is submitted `--dependency=afterok` on that job.

## Dynamic-Cs extraction for T4

No monitor logs the DynamicSFS coefficient; it lives in the particle VTPs as
per-particle array `C` (component 1 of 3; components 2–3 are
Lagrangian-average storage). Phase-averaged over 12 snapshots spanning rev
29→30 (script `orc:~/p018_gamma_dist_20260915/p018_sfs_C_avg.py`, log
`sfs_C_avg.log`):

| rung | mean C | \|Γ\|-weighted mean C | frac C=0 (clipped) | frac C≥0.999 |
|---|---|---|---|---|
| NT36 g25 | 0.14455 | 0.23976 | 0.542 | 0.0050 |
| NT72 g25 | 0.13078 | 0.23532 | 0.558 | 0.0087 |

**Chosen `SFS_CONST_CS = 0.14`** (NT36 run-average, rounded; one value for
BOTH rungs — the readout is NT-invariance, not level). Caveat disclosed:
ConstantSFS applies Cs uniformly where the dynamic model had clipped ~54% of
particles to zero, so T4's dissipation level is not matched to the dynamic
runs; `clipping_backscatter` remains active (driver passes `clippings` to
`ConstantSFS`, `rotor_hover_pressure_comparison.jl:507`). NoSFS deliberately
NOT used (052 ignition lesson; ruling 9).

## Arms

Common env (all runs): `WAKE_EXPINT=false`, `SIGMA_FLOOR_FRAC=0.25` (floor
0.00119 m), `SIGMA_CEIL=0.030`, `TRUNCATION_RADIUS_R=3.0`,
`MAX_PARTICLES=1500000` (headroom: NT72 peaked ~340k at rev 30; ~490k linear
projection at rev 60), `P018_REPO_OVERRIDE=~/wt018/FLOWPanel-mechtests`,
`P018_PROJECT_OVERRIDE=~/p018wtenv-expguard-gh200`. Launcher
`examples/run_dji9443_hover_ct_gpu.slurm.sh gh200 <case>`, submitted from the
mechtests worktree with `--constraint=arm` (mgh, 72 cpu, 192G, 1×GH200).
Baselines for every comparison: `p018_csarc_l3p0_3r_g25` (job 13704962) and
`p018_csarc_n2_nt72_l3p0_3r_srlx_g25` (job 13704963).

| test | run name | case tag | NT | extra env | wall |
|---|---|---|---|---|---|
| gate | `p018_csarc_l3p0_3r_g25_rgate` | `p018_csarc_l3p0` | 36 | `SFS_RLXF=0.005`, `P018_SETTLE_REVS=22`, `RESTART_STEP=1044`, `RESTART_NAME/PATH=…_3r_g25` | stage of 9 h chain |
| T1a | `p018_csarc_l3p0_3r_g25_s2` | `p018_csarc_l3p0` | 36 | `SFS_RLXF=0.005`, `P018_SETTLE_REVS=52`, `RESTART_STEP=1079`, same source | stage of 9 h chain |
| T1b | `p018_csarc_n2_nt72_l3p0_3r_srlx_g25_s2` | `p018_csarc_n2_nt72_l3p0` | 72 | `SFS_RLXF=0.0025031`, `P018_SETTLE_REVS=52`, `RESTART_STEP=2159`, `RESTART_NAME/PATH=…_srlx_g25` | 16:00 (afterok on chain) |
| T2 | `p018_csarc_n2_nt72_l3p0_3r_srlx_mrg2_g25` | `p018_csarc_n2_nt72_l3p0` | 72 | `SFS_RLXF=0.0025031`, `MERGE_EVERY=2`, `P018_SETTLE_REVS=22` | 14:00 |
| T4a | `p018_csarc_l3p0_3r_csfs_g25` | `p018_csarc_l3p0` | 36 | `SFS_CONST_CS=0.14`, `P018_SETTLE_REVS=22` (no `SFS_RLXF`: inert under ConstantSFS) | 05:00 |
| T4b | `p018_csarc_n2_nt72_l3p0_3r_csfs_g25` | `p018_csarc_n2_nt72_l3p0` | 72 | `SFS_CONST_CS=0.14`, `P018_SETTLE_REVS=22` | 14:00 |
| T5a | `p018_csarc_l3p0_3r_nv_g25` | `p018_csarc_l3p0_nv` | 36 | `SFS_RLXF=0.005`, `P018_SETTLE_REVS=22` | 05:00 |
| T5b | `p018_csarc_n2_nt72_l3p0_3r_srlx_nv_g25` | `p018_csarc_n2_nt72_l3p0_nv` | 72 | `SFS_RLXF=0.0025031`, `P018_SETTLE_REVS=22` | 14:00 |

Run-name collisions checked: all eight `data/<run>` paths free before
submission. mgh capacity note: 2 nodes × 1 GH200; jobs queue behind two
in-flight 022g jobs (walls 4 h / 10 h at submission).

## Scoring (per reset prompt)

Per completed run: (1) per-rev near/far Σ|Γ| time series
(`~/p018_gamma_dist_20260915/p018_gamma_ts.py`, 4-phase avg); (2)
phase-averaged per-bin ratios at rev 29→30 (12 phases,
`p018_gamma_dist_avg.py`); (3) CT̄ split windows via
`python3 scripts/p018_analyze.py m1 --revs {21 25, 26 30, 21 30}` (+ 31–60
windows for T1). Primary discriminator: the NT72 late far-wake re-acceleration
and NT36/NT72 far-wake divergence; a mechanism is implicated if the signature
shrinks >~50%.

## Submitted

2026-09-15 (~17:15 MDT), from `~/wt018/FLOWPanel-mechtests`, all
`--constraint=arm`. The T1 chain job runs gate → fidelity check → T1a as
stages of ONE job (`orc:~/p018_mech_tests_20260915/t1chain.slurm.sh`); T1b is
`--dependency=afterok` on it, so the long restarts cannot start unless the
scripted fidelity check passes. All queued behind two running 022g jobs
(mgh = 2 nodes × 1 GH200); banners to be verified on start (mandatory).

| job | test | run name / wrapper | status |
|---|---|---|---|
| 13712091 | gate + T1a chain (take 1) | `t1chain.slurm.sh` | FAILED in warmstart (dsigma2 skew, see Pins); banner otherwise fully correct (NT36, guard on, DynamicSFS rlxf=0.005, resume from step 1044, merge_every:1) |
| 13712092 | T1b take 1 | — | auto-CANCELLED by Slurm (dependency never satisfied) |
| 13712154 | gate + T1a chain (take 2, post-fix) | `t1chain.slurm.sh` → `…_rgate`, then `p018_csarc_l3p0_3r_g25_s2` | queued |
| 13712155 | T1b (afterok:13712154) | `p018_csarc_n2_nt72_l3p0_3r_srlx_g25_s2` | queued |
| 13712093 | T2 | `p018_csarc_n2_nt72_l3p0_3r_srlx_mrg2_g25` | RUNNING mgh-1-1; banner VERIFIED: NT:72, rlxf:0.16334, **merge_every:2**, guard=on, DynamicSFS(rlxf=0.0025031) |
| 13712094 | T4a | `p018_csarc_l3p0_3r_csfs_g25` | queued |
| 13712095 | T4b | `p018_csarc_n2_nt72_l3p0_3r_csfs_g25` | queued |
| 13712096 | T5a | `p018_csarc_l3p0_3r_nv_g25` | queued |
| 13712097 | T5b | `p018_csarc_n2_nt72_l3p0_3r_srlx_nv_g25` | queued |

Banner verification for the queued jobs happens as each starts (T4: expect
`SFS=ConstantSFS(Cs=0.14…)`; T5: `visc:false`/Inviscid; T1 chain: resume
lines from steps 1044/1079/2159).
