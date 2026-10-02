# 018 rlxf-scaled guarded reference pair — campaign provenance (2026-09-15)

Guarded, no-expint reference NT36/NT72 pair with **exact-rate SFS relaxation**,
to convert the `_3r_srlx` single-rung inference (exact-rate `SFS_RLXF` removes
39% of the NT climb: +1.51% → +0.92%) into a measured corrected slope on one
internally consistent stack and architecture. Authorised by Ryan 2026-09-15
("launch the small GPU campaign").

## Rationale / design

- `sfs_rlxf` enters `DynamicSFS` as a per-step relaxation with no other `dt`
  dependence; the campaign default 0.005 was never scaled with NT while the
  vortex `RELAX_RLXF` is exact-rate scaled. Exact-rate values:
  NT36 → 0.005 (nominal definition), NT72 → $1-\sqrt{0.995}$ = 0.0025031
  (bit-identical to `_3r_srlx`).
- Both rungs run guarded (`SIGMA_FLOOR_FRAC=0.25`, `SIGMA_CEIL=0.030`): the
  guard shifts CT (+0.39% on the expguard control), so the pair must carry it
  identically; the readout is the **slope**, not the levels vs the unguarded
  reference.
- Both rungs on the SAME arch (gh200/mgh) to exclude the untested
  h200-vs-gh200 FMM-kernel confound (the old reference pair and `_3r_srlx`
  mixed arms across arches).
- Pre-launch audit (2026-09-15, this session): on `_3r_sv_s1p5` NT36/72/144 at
  matched rev ≈30, birth σ (particles within 0.1R of the disk) is constant to
  <1% (6.92/6.96/6.95e-3 m) and nominal σ params are identical, but **gross
  sheds/rev is NOT constant**: 30,993 / 34,222 / 38,357 (+10.4%/doubling at
  NT36→72; the `max(1,ceil)` station quantization floor). This is common-mode
  between these arms and the reference family, so the rlxf slope comparison
  remains valid — but the NT axis carries a ~10%/doubling wake-refinement
  ride-along (recorded for the error-budget discussion).

## Pins (annotated tags `campaign/p018-rlxfscaled-20260915`, cluster repos)

Same commits as `campaign/p018-expguard-20260908`; the new tag is an alias for
discoverability. Worktrees reused from the expguard campaign (verified clean at
their pins, `git status --porcelain` empty, no queued/running job uses them).

| repo | worktree (cluster) | commit | note |
|---|---|---|---|
| FLOWVPM.jl | `~/wt018/FLOWVPM-expguard` | `7468712` | guard on euler (052c native) + euler_exp |
| FLOWPanel.jl | `~/wt018/FLOWPanel-expguard` | `d5dd772` | base `67b7383` + guard forwarding |
| FastMultipole | `~/wt018/FastMultipole-expguard` | `3da58a1a` | unmodified `unified-052` |

Local development counterparts carry `campaign/p018-expguard-20260908`
(FLOWVPM `21eeaaa`, FLOWPanel `7dee1ab`).

Verified before launch: `sigma_guard` is native on the plain
(`WAKE_EXPINT=false`) ReformulatedVPM euler path in FLOWVPM `7468712`
(`src/FLOWVPM_timeintegration.jl`, `_sigma_guard_params` / `_euler_cpu_reformulated!`).

## Julia environment

`~/p018wtenv-expguard-gh200` (dev-paths → the three worktrees above), depot
`~/fm052depot-gh200`, julia `~/julia/julia-1.11.7/bin/julia` (aarch64).
Unchanged from the expguard campaign; stack proven by smoke job 13610775 and
completed jobs 13610777/13610778.

## Arms

Common env: `WAKE_EXPINT=false`, `SIGMA_FLOOR_FRAC=0.25` (floor 0.00119 m),
`SIGMA_CEIL=0.030`, `TRUNCATION_RADIUS_R=3.0`, `MAX_PARTICLES=1500000`,
`P018_SETTLE_REVS=22`, `P018_REPO_OVERRIDE=~/wt018/FLOWPanel-expguard`,
`P018_PROJECT_OVERRIDE=~/p018wtenv-expguard-gh200`. Launcher
`examples/run_dji9443_hover_ct_gpu.slurm.sh gh200 <case_tag>`, submitted from
the FLOWPanel worktree; arch gh200 (partition mgh, 72 cpu, 192G).

| run name | case tag | NT | `SFS_RLXF` | wall |
|---|---|---|---|---|
| `p018_csarc_l3p0_3r_g25` | `p018_csarc_l3p0` | 36 | 0.005 | 05:00:00 |
| `p018_csarc_n2_nt72_l3p0_3r_srlx_g25` | `p018_csarc_n2_nt72_l3p0` | 72 | 0.0025031 | 14:00:00 |

New `_g25` run names — no clobbering of existing data dirs (checked).
Wall sizing: expguard NT36 gh200 ran 2h01; `_3r_srlx` NT72 h200 ran 5h05.

## Results (2026-09-15, scored same session)

Both COMPLETED (13704962: 2h16; 13704963: 4h41). `p018_analyze.py m1 --revs 21 30`:

| rung | CT̄ (revs 21–30) | 95% CI | per-rev std |
|---|---|---|---|
| NT36 `_3r_g25` | 0.070519 | [0.070478, 0.070577] | 0.000106 |
| NT72 `_3r_srlx_g25` | 0.071417 | [0.071306, 0.071548] | 0.000196 |

**Corrected climb = +1.27%/doubling** (vs unguarded reference +1.51%, and vs
the srlx-inferred +0.92%). Block drift 0.000%, non-monotone, both rungs.

**Verdict: exact-rate `SFS_RLXF` is NOT the carrier.** On a clean guarded
same-arch pair it removes only ~16% of the climb (0.24 points, marginal vs the
~±0.15-point CI on the slope), far less than the 39% inferred from the
cross-stack srlx comparison (which mixed guard state, stack, and arch). The
srlx-inferred 39% is superseded by this measurement. Residual +1.27%/doubling
remains unexplained; the shed-count quantization excess (+10.4%/doubling, see
audit above) and the untested FMM near-set adequacy lead are now the front
runners. Guard level effect also flipped sign vs the exp arm: guarded NT36
reference sits −0.36% below the unguarded 0.070775 (exp arm: +0.39%) — one
more reason levels across stacks must not be compared, only slopes.

Full per-step VTP sets (1080 / 2160 files) were still live in both run dirs at
scoring time — harvest or protect (`vtk_protect_list.txt`) before the sweeper
culls to newest-36.

## Scoring

`scripts/p018_analyze.py`, settled window revs 21–30 per `decision_rules.md`.
Readout: guarded corrected climb NT36→NT72, compared against the unguarded
reference climb +1.51% and the srlx-inferred +0.92%.

## Submitted

2026-09-15, from `~/wt018/FLOWPanel-expguard`, both `--constraint=arm`
(required by the mgh submit plugin; first attempt without it was rejected,
nothing submitted).

| job | run name | state at submit |
|---|---|---|
| 13704962 | `p018_csarc_l3p0_3r_g25` | RUNNING mgh-1-1; banner `WAKE_EXPINT=false`, `guard=on` (floor 0.00119 m, ceil 0.03 m), `DynamicSFS(rlxf=0.005)`; ~4.9 s/step |
| 13704963 | `p018_csarc_n2_nt72_l3p0_3r_srlx_g25` | RUNNING mgh-1-2; banner `WAKE_EXPINT=false`, `guard=on` (floor 0.00119 m, ceil 0.03 m), `DynamicSFS(rlxf=0.0025031)`; ~4.8 s/step |
