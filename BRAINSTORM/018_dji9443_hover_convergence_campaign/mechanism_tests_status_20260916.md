# 018 mechanism tests wave 1 — status (2026-09-16)

Monitoring/scoring session per `mechanism_tests_reset_prompt_20260915b.md`.
Provenance: `mechanism_tests_provenance_20260915.md` (pins unchanged, tag
`campaign/p018-mech-tests-20260915`, FLOWPanel `5cdf058`). Baselines
throughout: g25 pair `p018_csarc_l3p0_3r_g25` (NT36, 13704962) /
`p018_csarc_n2_nt72_l3p0_3r_srlx_g25` (NT72, 13704963).

## Job outcomes

| job | test | run | outcome |
|---|---|---|---|
| 13712093 | T2 mrg2 | `…_srlx_mrg2_g25` | COMPLETED 5:30 (vs 4:41 baseline; +17%, consistent with less merging). Banner ✓ (merge_every:2). GATE clean. |
| 13712094 | T4a | `…l3p0_3r_csfs_g25` | COMPLETED 1:53. Banner ✓ `SFS=ConstantSFS(Cs=0.14, clippings=backscatter)`. GATE clean. |
| 13712095 | T4b | `…nt72_l3p0_3r_csfs_g25` | COMPLETED 3:45. Banner ✓ ConstantSFS(Cs=0.14). GATE clean. |
| 13712096 | T5a | `…l3p0_3r_nv_g25` | COMPLETED 3:25. Banner ✓ `visc:false`, no CoreSpreading line, DynamicSFS rlxf=0.005. GATE clean. |
| 13712097 | T5b | `…nt72_l3p0_3r_srlx_nv_g25` | **NODE_FAIL** at 0:38 (step ~379/2159). Banner was ✓ (visc:false, NT72, rlxf 0.0025031). Not requeued (`--no-requeue`). Partial dir moved to `…_nv_g25_nodefail13712097`; **resubmitted as 13733314** (identical submit line). |
| 13712154 | T1 chain take 2 | rgate → T1a | Stage 1 (rgate) COMPLETED cleanly (GATE clean, resume from 1044 ✓). Stage 2 **FAILED for infra reasons only**: job-env `python3` has no numpy (`ModuleNotFoundError`, see `.err`). Chain exited 42 before T1a. |
| 13712155 | T1b take 2 | — | auto-CANCELLED (dependency). |

## T1 restart-fidelity gate: **PASS** (manual rerun of the scripted check)

`gate_check.py` run manually on the login node (numpy available there)
against the completed rgate data — output:
thrust comp CFx, overlap 35 steps (1045–1079), **max rel ΔCF (skip first 3)
= 2.107e-5, tol 2e-4** — over an order of magnitude inside tolerance, seam
transient decaying (6e-10 at step 1045 growing smoothly to 2.1e-5 at 1079,
i.e. slow drift, not a jump); wake at step 1079 identical between original
and rgate: n=311310, Σ|Γ|=5.9700 both. Restart machinery is faithful ⇒ per
charter, T1 proceeds.

**Resubmissions (2026-09-16 ~09:30 MDT, from the mechtests worktree):**
T1a as standalone job **13733312** (`~/p018_mech_tests_20260915/t1a.slurm.sh`
= stage 3 of the chain verbatim, wall 08:00, gate stage dropped since PASS is
established); T1b **13733313** (original submit line minus the dependency);
T5b **13733314**. All pending at write time; banners to be verified on start.

## Scoring — completed runs vs g25 baselines

Scripts: `~/p018_mech_tests_20260915/p018_gamma_ts2.py` / `p018_gamma_dist_avg2.py`
(parameterized copies of the 0915 scripts, RUNS via argv; logs
`gamma_ts_wave1.log`, `gamma_avg_wave1.log`, `ct_windows_wave1.log`,
`ct_windows_baselines.log` in `~/p018_mech_tests_20260915/`). Full VTP sets
were intact (sweeper had not culled).

### 1. Far-wake (x ≥ 1.5R) Σ|Γ| per-rev time series (4-phase avg)

| rev | NT36 base | NT72 base | T2 (NT72 mrg2) | T4a (NT36 cSFS) | T4b (NT72 cSFS) | T5a (NT36 inv) |
|---|---|---|---|---|---|---|
| 16 | 3.33 | 3.18 | 3.14 | 3.61 | 3.52 | 4.24 |
| 20 | 3.53 | 3.53 | 3.24 | 3.27 | 3.20 | 4.41 |
| 24 | 3.75 | 3.72 | 3.64 | 3.58 | 3.70 | 4.48 |
| 26 | 3.86 | 3.95 | 3.75 | 3.37 | 4.39 | 4.64 |
| 28 | 3.94 | 4.17 | 3.86 | 3.62 | 4.61 | 4.96 |
| 29 | 3.90 | 4.33 | 3.91 | 3.61 | 4.94 | 5.05 |

- **T2**: NT72 far-wake excess over NT36 at rev 29 drops from **+0.43 to
  +0.01 (−98%)**; late slope (revs 24→29) drops from +0.12/rev to +0.05/rev
  (NT36 baseline: +0.03/rev). The rev-28 crossing is gone. Particle count
  monotone (~316k, vs baseline's 336k→323k→339k non-monotone).
- **T4**: divergence **amplified** — NT72−NT36 far excess at rev 29 is
  **+1.34** (3× baseline), T4b late slope ~+0.25/rev. ConstantSFS also
  shrinks wake counts ~40% (NT36 rev 29: 184k vs 311k) — uniform Cs=0.14
  removes far more than the dynamic model (which clips ~54% of particles
  to C=0).
- **T5a** (NT36 only so far): whole far wake inflated (~+29% at rev 29,
  Σ|Γ| total 7.30 vs 5.96) and still climbing — viscous decay clearly a
  major sink; pair verdict awaits T5b.

### 2. Phase-avg per-bin NT72/NT36 ratios, rev 29→30 (12 phases)

TOTAL-row (count / Σ|Γ|) and worst far bins:

| pair | count | Σ\|Γ\| | 1.5–2R Σ\|Γ\| ratio | 2.5–3R Σ\|Γ\| ratio |
|---|---|---|---|---|
| baseline g25 | 1.09 | 1.05 | 1.33 | 1.26 |
| T2 vs NT36 base | 1.02 | 0.98 | 1.02 | 1.14 |
| T4b vs T4a | 1.17 | 1.29 | 1.53 | 1.96 |

Same story: T2 flattens the far-bin excess to near-parity; ConstantSFS
roughly doubles it. (Earlier windows 15→16 / 23→24 in the log show T2 ≈
baseline there — the fix is specifically to the late-onset excess.)

### 3. CT̄ windows (`p018_analyze.py m1`)

| run | 21–25 | 26–30 | 21–30 | climb vs NT36 partner (21–30) |
|---|---|---|---|---|
| NT36 base | 0.070571 | 0.070454 | 0.070519 | — |
| NT72 base | 0.071394 | 0.071446 | 0.071417 | **+1.27%** |
| T2 (NT72 mrg2) | 0.071403 | 0.071358 | 0.071383 | **+1.23%** (vs NT36 base) |
| T4a (NT36 cSFS) | 0.070661 | 0.070767 | 0.070708 | — |
| T4b (NT72 cSFS) | 0.071303 | 0.071316 | 0.071308 | **+0.85%** (vs T4a) |
| T5a (NT36 inv) | 0.070846 | 0.070345 | 0.070623 | (pair pending T5b) |

## Verdicts (wave 1, partial — T1/T5b pending)

- **T2 / merging cadence: IMPLICATED for the far-wake Γ signature, NOT for
  the CT climb.** Rate-matching NT72's merge cadence to NT36 (36 merges/rev)
  removes ~98% of the late far-wake Σ|Γ| excess and the count non-monotonicity
  — yet CT̄ is statistically unchanged (+1.23% vs +1.27% climb, −0.05% on
  CT̄ itself). **The far-wake Σ|Γ| divergence and the CT climb are decoupled:
  the Γ-distribution signature is a merge-cadence artifact, but it is not the
  carrier of the +1.27%/doubling.**
- **T4 / SFS dynamics: NOT the far-wake divergence mechanism (it was damping
  it), but implicated in ~1/3 of the CT climb** (1.27% → 0.85%/doubling).
  Consistent with the dynamic-Cs NT-dependence measured in provenance (mean C
  0.145 @ NT36 vs 0.131 @ NT72 — less dissipation at NT72). Caveat: the
  uniform Cs=0.14 is not dissipation-matched to the dynamic runs (54%
  clipping), and it changes the wake globally (counts −40%), so this is
  directional evidence, not a clean attribution. Below the >50% threshold ⇒
  "partially implicated".
- **T5 / viscous:** NT36 arm shows core spreading is a first-order wake-Γ
  sink (+29% far wake when off) but small on CT (+0.15% on NT36). Divergence
  verdict needs T5b (13733314).
- **T1 / transient length:** restart fidelity PASS; long runs in flight
  (13733312/13733313).

## Wave-2 recommendation

T2 positive ⇒ **T3 (merging off) is triggered per the pre-approval** — it is
the clean limit of the cadence axis and would confirm the merge-removal
asymmetry as the Γ-signature driver. Feasibility flag before submitting:
baseline NT72 peaked ~340k particles at rev 30 *with* per-step merging; with
merging fully off expect several× that — check against `MAX_PARTICLES=1.5e6`
and GH200 memory, and consider NT36-only or shorter horizon if projections
are marginal. However, given T2 shows the Γ signature does **not** carry the
CT climb, the higher-value wave-2 axis is arguably the SFS one (e.g. a
dissipation-matched ConstantSFS pair, or dynamic-Cs with NT-matched
averaging) plus the T1 long-window readout — worth deciding with Ryan once
T1/T5b land.

## Housekeeping / pending

- T2/T4a/T4b/T5a full VTP sets still on disk at scoring time; Γ scripts run.
  If tails matter later, harvest to
  `/nobackup/archive/usr/rander39/FLOWPanel_runs/` before the sweeper culls.
- rgate run dir `p018_csarc_l3p0_3r_g25_rgate` retained (35 VTPs).
- T5b node-failed partial: `…_nv_g25_nodefail13712097` (delete after 13733314
  lands cleanly).
- `gate_check.py` numpy failure in job env: if a scripted in-job gate is ever
  reused, load a numpy-capable python (or venv) in the slurm script first.
- To verify on start (mandatory): 13733312 (T1a — expect resume from step
  1079, NT36, rlxf 0.3/DynamicSFS 0.005, guard=on), 13733313 (T1b — resume
  from 2159, NT72, rlxf 0.16334/DynamicSFS 0.0025031), 13733314 (T5b —
  visc:false, no CoreSpreading).
- Notebook/ledger backlog (stale since 2026-09-08) now additionally includes:
  wave-1 launch + warmstart-skew incident (0915), gate PASS, T5b node
  failure + resubmits, and these wave-1 partial verdicts. Awaiting Ryan's
  approval to draft.

---

# 2026-09-16 evening session — free analyses 1–3 (Task B) + wake-transplant audit (Task C.4 prep)

Per `mechanism_tests_reset_prompt_20260916.md`. In-flight at session time:
13733312 (T1a) RUNNING mgh-1-1 (step ~1356/2159 at 00:39 elapsed, log clean —
the only "nan" matches are the banner's unused `das_chord/das_uniform/das_beta`
knobs), 13733313 (T1b) RUNNING mgh-1-2 (step ~2421/4319, log clean), 13733314
(T5b resubmit) PENDING behind T1a. Banner verification for 13733314 still owed
at its start.

## Staging note (sweeper had culled the baselines)

Both baseline run dirs are marked `ARCHIVED.txt` and hold only their last 5
particle VTPs (NT36: 1075–1079; NT72: 2155–2159) — not enough for the
12-phase rev-29→30 average. The 24 needed VTPs (NT36 steps 1044+3j, NT72
2088+6j, j=0..11) were re-extracted from the archive tarballs
(`/nobackup/archive/usr/rander39/FLOWPanel_runs/projects_FLOWPanel.jl/`, 26G/54G)
to `orc:~/p018_mech_tests_20260915/extracted/` (script
`extract_rev29_vtps.sh`, log `extract_rev29.log`). Delete after wave-2 scoring
settles.

## B3 — spanwise loading split (monitor03 bound circulation)

No sectional *force* monitor exists (monitors are force/bound_circulation/
wake_health), but `monitor03_bound_circulation` carries per-blade, per-section
`circulation_te` at every step (`circulation_slice` is all-NaN — unused).
Script `orc:~/p018_mech_tests_20260915/p018_spanwise_loading.py`, log
`spanwise_loading.log`. Window-averaged NT72/NT36 ratios of section-mean Γ_te
(both blades symmetric; representative sections):

| r/R | revs 21–25 | revs 26–30 |
|---|---|---|
| 0.11 | 1.007 | 1.032–1.035 |
| 0.22 | 1.012 | 1.023–1.026 |
| 0.40 | 1.009 | 1.018–1.020 |
| 0.60 | 1.001–1.002 | 1.009–1.010 |
| 0.70 (local dip) | 0.994–0.997 | 1.002 |
| 0.79–0.81 (bump) | 1.023–1.031 | 1.024–1.029 |
| 0.85–0.90 | 1.003–1.011 | 1.008–1.011 |
| 0.94–0.99 | 0.994–0.996 | 0.996–0.999 |
| **Σ\|Γ_te\| (integrated)** | **1.0068** | **1.0146** |

Findings: (1) the integrated bound-circulation excess **doubles from the
21–25 to the 26–30 window (+0.68% → +1.46%)**, tracking the late-onset CT
climb (CT̄ ratio +1.17% → +1.41% over the same windows) — bound circulation is
a faithful spanwise proxy for the carrier. (2) The **late-developing part is
inboard/mid-board weighted** (largest ratio growth at r/R ≲ 0.4, ~+2.5 pts;
mid-board ~+1 pt), while the tip (r/R > 0.9) is flat and slightly *below*
parity, and the r/R ≈ 0.80 bump (~+2.5%) is **window-independent** (already
present in 21–25, not part of the late onset). A uniform-inflow change would
scale sections more evenly; an inboard-weighted growth is what a slowly
developing change in induced inflow over the inner disk looks like —
consistent with the wake-geometry hypothesis, not with tip-loading noise.

## B1 + B2 — wake geometry profiling (12-phase avg, rev 29→30, matched phase)

Script `orc:~/p018_mech_tests_20260915/p018_wake_geometry.py`, log
`wake_geometry.log`. Cylindrical r about axial axis x, weights |Γ|.

(a) x-slice structure (Δx = 0.25R; Γ-weighted mean radius r̄/R and weighted
percentiles of r/R):

| x/R | NT36 r̄ | NT72 r̄ | NT36 p50 | NT72 p50 | NT36 p75 | NT72 p75 |
|---|---|---|---|---|---|---|
| 0.12 | 0.843 | 0.731 | 0.777 | 0.767 | 0.818 | 0.802 |
| 0.38 | 0.917 | 0.672 | 0.733 | 0.715 | 0.821 | 0.757 |
| 0.62 | 0.867 | 0.691 | 0.725 | 0.715 | 0.808 | 0.773 |
| 0.88 | 0.895 | 0.737 | 0.751 | 0.748 | 0.989 | 0.839 |
| 1.12 | 0.957 | 0.871 | 0.805 | 0.814 | 1.057 | 0.997 |
| 1.38 | 1.049 | 1.041 | 0.866 | 0.976 | 1.270 | 1.323 |
| 1.88 | 0.987 | 1.101 | 0.903 | 0.997 | 1.196 | 1.414 |
| 2.37 | 1.143 | 1.264 | 1.037 | 1.232 | 1.462 | 1.594 |

Readout: **within the first tip passage (x ≲ 1R) the two rungs' wakes are
geometrically different.** NT72's slipstream is modestly more contracted at
the median (p50 lower by 0.010–0.018 R) and markedly tighter in the outer
quartile (p75 lower by 0.05–0.15 R); the Γ-weighted mean differs far more
(r̄ lower by 0.11–0.25 R) because NT36 carries a substantial high-radius Γ
tail near the disk plane that NT72 lacks. The pattern **reverses downstream**:
for 1.4 ≲ x/R ≲ 2.5 NT72 sits wider (p50 higher by 0.09–0.20 R). Caveat: r̄
is tail-sensitive and Σ|Γ|-based weights are discretization-dependent (pps 12
vs 6, different merge cadence), so the robust statement is the percentile one.

(b) tip-band axial |Γ| profile (r ∈ [0.85, 1.15]R, Δx = 0.02R, matched phase
j=0): no clean discrete spiral peaks are resolvable in either rung (profiles
rise monotonically with noise — passages are already diffused/merged at rev
29), so **first-passage axial spacing is not extractable** from these
snapshots. The band-integrated |Γ| is ~2× larger for NT72 over
0.75 ≲ x/R ≲ 1.5, consistent with (c) below.

(c) radial rebin of Σ|Γ| restricted to x < 1.5R (phase-avg):

| r/R bin | NT36 | NT72 | ratio |
|---|---|---|---|
| 0.0–0.2 | 0.0769 | 0.0663 | 0.86 |
| 0.2–0.4 | 0.1016 | 0.1074 | 1.06 |
| 0.4–0.6 | 0.0384 | 0.0500 | 1.30 |
| 0.6–0.8 | 0.8916 | 0.8477 | 0.95 |
| 0.8–1.0 | 0.4813 | 0.5081 | 1.06 |
| 1.0–1.2 | 0.0911 | 0.1371 | **1.51** |
| 1.2–1.4 | 0.0712 | 0.0844 | 1.19 |
| total | 1.7522 | 1.8010 | 1.028 |

Near-total Σ|Γ| is nearly NT-invariant (+2.8%; the earlier 4-phase rev-29
series had it closer to parity — sampling difference, both small), but there
is clear **radial redistribution**: NT72 has less Σ|Γ| in the root (0–0.2R,
−14%) and in the main sheet band (0.6–0.8R, −5%), and substantially more in
the outboard band beyond the tip radius (1.0–1.4R, +19…+51%).

**Interpretation w.r.t. the geometry hypothesis.** The surviving suspect
class (wake-borne, slowly developing, not total far-wake strength) now has a
concrete signature: the two rungs place their near/mid wake differently —
NT72's slipstream is tighter/more contracted within the first passage while
carrying more Γ in the 1.0–1.4R outboard band of the near field, and its
mid/far wake rides wider. The disk response matches: the late-developing
loading excess is inboard-weighted (B3), which is where induced-inflow changes
from near-wake geometry act. This strengthens the geometry/induced-inflow
candidacy for the CT carrier and sharpens what T1 must answer: if the
geometry difference (and CT gap) saturates by rev ~45–60 it is a
transient-length effect; if it persists it is a resolution-linked steady-state
difference. Sign attribution (which geometry feature raises inboard loading)
is left open — that is what test 4 (wake transplant) would discriminate
directly.

## Task C.4 — wake-transplant cross-restart feasibility audit (prep only, NOT submitted)

Audited `src/FLOWPanel_warmstart.jl` + `FLOWPanel_replay.jl` +
`FLOWPanel_metadata.jl` at the campaign pin (`5cdf058`). Corrections to the
charter's priors: **both** arms write both wake surfaces (`wake1.1` and
`wake1.2` VTS exist for NT36 and NT72 — the "NT36 has only wake1.1" caveat is
wrong), and the conversion fingerprints are **identical** across arms
(`LegacyEdgeJumpConversion`, parameter-free) so the fingerprint gate passes.
The real arm difference is wake *rows*: NT36 VTS extent `0 1 0 40` (1 wake
row), NT72 `0 2 0 40` (2 rows); spanwise node count (41) matches.

Mechanics that matter (file:line at pin): `restart_step` is interpreted in the
**target arm's** step/time coordinates (kinematic replay `warmstart.jl:527-531`,
`start_step=restart_step+1` :694), so a donor checkpoint must be staged under
the **matched-azimuth target step index**: NT36 step 1079 = NT72 step 2158
(both 29.972 revs — NT72's final step 2159 has NO NT36 counterpart; use 2158).
Metadata restores key on `i_step == restart_step` and hard-check
`step_identity == restart_step` (`replay.jl:534-535`, `kutta.jl:1954`), so the
staged `.metadata.toml` needs those keys rewritten. The PVD manifest may be
stale — an explicit `RESTART_STEP` whose body VTU exists proceeds with a warn
(`warmstart.jl:469-477`). Frame state: target replays its own kinematics and
the manifest cross-check either matches (same azimuth) or falls back to replay.
Particle field restore is arm-agnostic (capacity-checked only; all 9 required
fields + the 4-array `split_*` set are present in this pin's VTPs).

- **NT72 ← NT36 (donor NT36 rev-30): feasible with staging only, no code
  change.** Panel wake loads 1 donor row into the 2-row buffer
  (`nwakes[]=1`, matches donor metadata `active_row_count=1`); the second row
  regenerates by shedding within a step.
- **NT36 ← NT72 (mirror): one real blocker.** Donor VTS has 2 rows; the
  loader copies `nodes[:,1:dim1,:]` with dim1=3 into NT36's 2-node-row buffer
  → BoundsError (`warmstart.jl:196`). Options: (i) stage a row-truncated VTS
  (keep the newest row) + rewrite `active_row_count/live_rows` — fiddly,
  needs row-order confirmation; (ii) a small load-time row clamp (code
  change, wave-2 approval); (iii) run the transplant one-directional
  (NT72←NT36 already discriminates evaluation-carrier vs wake-carrier).

Proposed staging + submit (NT72←NT36 direction; **awaiting Ryan's approval**):

1. Stage `data/p018_xr72f36_stage/` with symlinks renaming the NT36 step-1079
   quadruplet (`_body1.1079.vtu`, `_wake1.{1,2}.1079.vts`,
   `_wake1_particles.1079.vtp`) to the staged name at step **2158**, copy
   `_body1.pvd`, and write a key-rewritten `.metadata.toml`
   (i_step/step_identity 1079→2158 in the step records).
2. Submit (modeled on the T1b submit line; 3 revs past the seam → settle 25 ⇒
   33-rev t_range, 2376 steps, ~1.5 h of stepping — wall 04:00):
   `sbatch --parsable --constraint=arm --job-name=fp-018mech-xr72f36
   --time=04:00:00 --export=ALL,WAKE_EXPINT=false,SIGMA_FLOOR_FRAC=0.25,
   SIGMA_CEIL=0.030,TRUNCATION_RADIUS_R=3.0,MAX_PARTICLES=1500000,
   P018_REPO_OVERRIDE=$HOME/wt018/FLOWPanel-mechtests,
   P018_PROJECT_OVERRIDE=$HOME/p018wtenv-expguard-gh200,
   SFS_RLXF=0.0025031,P018_SETTLE_REVS=25,RESTART_STEP=2158,
   RESTART_NAME=p018_xr72f36_stage,RESTART_PATH=data/p018_xr72f36_stage,
   P018_RUN_NAME=p018_csarc_n2_nt72_l3p0_3r_srlx_xr72f36
   examples/run_dji9443_hover_ct_gpu.slurm.sh gh200 p018_csarc_n2_nt72_l3p0`
   Readout: per-step CFx from the seam onward — snaps to NT72's level within
   ~a rev ⇒ per-step evaluation carrier; starts at NT36's level and drifts
   over revs ⇒ wake-evolution carrier.
3. Mirror direction: decision needed on options (i)/(ii)/(iii) above.

Open validation before real submission: dry-check that the driver accepts
`RESTART_NAME≠P018_RUN_NAME` staging (T1 used same-name restart; the env knobs
`RESTART_STEP/NAME/PATH` at `rotor_hover_pressure_comparison.jl:1368-1370`
pass through generically) and that `P018_SETTLE_REVS=25` yields a 33-rev
t_range on this case (T1 pattern: settle = total − 8).

## Task C.5 / C.6 — prep (proposed lines only, awaiting approval)

**C.5 single-knob restart perturbations** (same-arm restarts from own rev-30
state — none of the cross-arm staging above is needed; RESTART_STEP=1079 (NT36)
/ 2159 (NT72), settle 27 ⇒ 35-rev t_range, 5 revs past the seam, wall 03:00 /
05:00). `RELAX_RLXF` and `SFS_RLXF` are plain env knobs — override directly;
**pps/nwakerows are case-arm constants** ⇒ per the standing ops rule those two
flips need a new-case-arms examples commit first (wave-2 code change, not
prepped here). Highest-value pair given T4: rlxf-scaling checks —
NT72 restart with `RELAX_RLXF=0.3` (NT36's unscaled value; tests whether the
½-power srlx compensation leaves residual NT-dependence) and NT72 restart with
`SFS_RLXF=0.005` (same, SFS channel). Template = the T1b submit line with
`P018_SETTLE_REVS=27`, `--time=05:00:00`, the flipped knob, and run names
`…_srlx_g25_rx30` / `…_srlx_g25_sfsrx5`. Watch the late-window CT response
(revs 31–35 vs baseline's 26–30).

**C.6 per-rung matched ConstantSFS**: reuse the 13712094/13712095 submit lines
verbatim with `SFS_CONST_CS=0.145` (NT36 arm) and `SFS_CONST_CS=0.131` (NT72
arm), run names `…_csfs145_g25` / `…nt72…_csfs131_g25`, walls 05:00/14:00.
Readout: climb returns to ~1.27% ⇒ T4's reduction was SFS *level*; stays
~0.85% ⇒ *dynamics/fluctuations*.

## 2026-09-16 ~22:15 MDT — T1a landed clean; T5b started, banner PASS; T1a scored

**T1a (13733312, NT36 revs 30→60)** COMPLETED ~22:08 MDT. GATE clean
(`gpu_gemv=1080 cpu_gemv=0 nan_lines=0 dispatcher_rc=0`), launcher rc=0.
**T5b (13733314)** started immediately on the freed node (mgh-1-1); the
automated banner watcher verified and recorded **PASS** at 22:08:16
(`visc:false`, no CoreSpreading line, NT:72, rlxf 0.16334, DynamicSFS
rlxf=0.0025031, guard=on, mechtests repo —
`~/p018_mech_tests_20260915/t5b_banner_check.txt`). T1b still RUNNING
(expected to land by ~11:20 MDT 09-17).

**T1a CT̄ windows** (`ct_windows_t1a.log`; baseline windows repeated for
context):

| window | NT36 CT̄ | 95% CI |
|---|---|---|
| 21–25 (base) | 0.070571 | — |
| 26–30 (base) | 0.070454 | — |
| 31–40 | 0.070329 | [0.070241, 0.070406] |
| 41–50 | 0.070429 | [0.070343, 0.070485] |
| 51–59 | 0.070495 | [0.070417, 0.070565] |

**The baseline's "NT36 drifts down" REVERSES in the extension**: CT̄ bottoms
in revs 31–40 and recovers monotonically through rev 59, ending near the
26–30 level (0.070495 vs 0.070454; the 31–40 → 51–59 recovery of +0.00017
clears the CIs). NT36 is still slowly equilibrating at rev ~50 with an
oscillatory (undershoot-recover) approach, settling near ≈0.0705.

**T1a far-wake Σ|Γ| (x ≥ 1.5R, 4-phase avg, `gamma_ts_t1a.log`)**: 3.91 (rev
30) → 4.63 (40) → 4.75 (50) → 4.84 (59). Slope collapses from +0.072/rev
(revs 30–40) to ~+0.010/rev (50–59): the NT36 far wake **largely saturates by
rev ~45** with mild residual creep. Near wake stable (1.87–1.95).

Implication so far (pair verdict awaits T1b): the NT36 half of the two-sided
late-onset divergence looks like a **transient** — its rev-21–30 downdrift
was an undershoot during far-wake fill-in, not a persistent trend. Whether
the +1.27% gap closes now rests on NT72's late windows (does 0.07145 hold, or
does it relax back toward NT36's ≈0.0705?).

Housekeeping: T1a's rev-44 and rev-59 12-phase VTP windows (24 files) were
harvested to `~/p018_mech_tests_20260915/extracted/p018_csarc_l3p0_3r_g25_s2/`
before the sweeper culls (pair per-bin scoring vs T1b needs rev 44; the
newest-36 cull retains only steps 2124–2159). T1b/T5b completion watcher
running locally; T1b Γ time series must be scored immediately at land (its
rev-30–58 VTPs die in the newest-36 cull).

---

# 2026-09-17 — T1b and T5b landed; T1/T5 verdicts. **T1b IGNITED at rev ~47.5.**

Both scored automatically at landing by the ORC-side watchers
(`t1b_land_watch.sh`/`t5b_land_watch.sh`; logs `gamma_ts_t1b.log`,
`gamma_avg_t1pair.log`, `ct_windows_t1b.log`, `gamma_ts_t5b.log`,
`gamma_avg_t5pair.log`, `ct_windows_t5b.log` in
`~/p018_mech_tests_20260915/`).

## T1b (13733313, NT72 revs 30→60): COMPLETED 03:54 MDT (8:27) — but the run
**destabilized at rev ~47.5**. GATE printed clean (`gpu_gemv=2160
nan_lines=0`) — a NaN-free ignition: exit code and GATE are NOT health, the
monitors are. Signature is the 052-style stretching runaway
(monitor04: max |Γ|/σ² 8.0e2 → 1.4e3 (rev 45) → 1.6e4 (rev 47) → 2.2e9
(rev 48); max_u 25 → 55 → 2.4e6 m/s; particle count collapses 411k → 24k as
the guard/truncation culls). Per-rev CT is clean through rev 47 and garbage
from rev 48 (±1e3). SIGMA_CEIL guard was ON.

**Valid-window CT̄ (revs ≤ 47) vs T1a:**

| window | NT36 (T1a) | NT72 (T1b) | gap |
|---|---|---|---|
| 26–30 (base) | 0.070454 | 0.071446 | +1.41% |
| 31–40 | 0.070329 | 0.071244 | **+1.30%** |
| 41–50 / 41–47 | 0.070429 | 0.071268 (41–47) | **+1.19%** |
| 51–59 | 0.070495 | (post-ignition, invalid) | — |

**Far wake (x ≥ 1.5R, 4-phase avg, pre-ignition):** NT72 climbs 4.41 (rev
30) → 5.27 (40) → 6.13 (46), slope **+0.11/rev at rev 46, no saturation** —
vs NT36's saturation at +0.010/rev by rev ~45. Per-bin pair ratio at rev
44→45 (12 phases): NT72/NT36 total Σ|Γ| 1.21; far bins 2.5–3R 1.26, 3–3.5R
1.15 — the divergence keeps compounding right up to ignition.

**T1 verdict:** restart fidelity PASS (gate, previously); the extension
answers the central arbiter question: **the CT climb is a PERSISTENT carrier,
not transient length.** NT36's baseline downdrift was a fill-in transient
(undershoot then recovery to ≈0.0705); NT72 holds ≈0.0712–0.0713 with no
relaxation toward NT36 through rev 47. The gap narrows only marginally
(+1.41% → +1.19%). **New finding: the NT72 g25 rung is metastable at long
horizon** — its far-wake Σ|Γ| accumulation (the T2-diagnosed merge-cadence
artifact) never equilibrates and the run ignites at rev ~47.5. Corollary
prediction (cheap wave-2 falsifiable): a long mrg2 run (rate-matched merging)
should NOT ignite, since T2 removed the far-wake accumulation.

## T5b (13733314, inviscid NT72): COMPLETED 03:09 MDT (5:03), banner PASS
(auto-verified at start), GATE clean, monitors healthy through rev 30.

**T5 pair CT̄** (windows as reported by m1; T5b's last window is 26–29):

| window | T5a (NT36 inv) | T5b (NT72 inv) | climb |
|---|---|---|---|
| 21–25 | 0.070846 | 0.071460 | +0.87% |
| 26–30 / 26–29 | 0.070345 | 0.071208 | +1.23% |
| 21–30 / 21–29 | 0.070623 | 0.071348 | **+1.03%** |

**T5 far wake:** inviscid inflates both rungs (T5b rev-29 far Σ|Γ| 5.94 vs
viscous 4.33) and the NT72/NT36 far excess persists (rev-29 far total +17.6%
vs baseline's +11%; per-bin rev 29→30 ratios 1.5–2R 1.24, 2.5–3R 1.39 vs
baseline 1.33/1.26).

**T5 verdict: viscous machinery (core spreading) is NOT the CT carrier and
NOT the divergence mechanism.** The climb persists inviscid at +1.03%/doubling
(vs +1.27% viscous; reduction well below the 50% implication threshold), the
late-onset window pattern persists (+0.87% → +1.23%), and the far-wake
divergence persists. Core spreading is a first-order Σ|Γ| sink (+29–42%
far-wake inflation when off) but nearly CT-neutral.

## Wave-1 synthesis (all tests now scored)

| axis | verdict on CT climb |
|---|---|
| T2 merge cadence | NOT the carrier (explains far-wake Γ signature only) |
| T4 SFS dynamics | partially implicated (~1/3; level-vs-dynamics open → C.6) |
| T5 viscous | not implicated (+1.03% persists inviscid) |
| T1 transient length | REJECTED as explanation — carrier is persistent to rev 47 |

Remainder: a persistent, wake-borne, geometry-correlated carrier (the 09-16
free analyses: inboard-weighted late loading growth + near-wake contraction /
outboard-band redistribution) — the transplant (C.4), n0 convert-at-shed
ladder, and per-rung Cs (C.6) are the live discriminators, plus the new
stability axis (mrg2 long run).

## Housekeeping

- Harvested to `~/p018_mech_tests_20260915/extracted/…_srlx_g25_s2/`: T1b
  12-phase windows revs 44, 46, 47, 48, 59 (60 files; 46–48 straddle
  ignition for post-mortem; 59 is post-blowup). If a full ignition
  post-mortem is wanted, tar T1b's VTP set to the archive before the sweeper
  culls (~2160 files on disk at scoring time).
- T5b node-fail partial `…_nv_g25_nodefail13712097`: 13733314 landed clean so
  it is now deletable per the 09-16 note, but the delete was blocked by the
  local permission classifier — **left for Ryan**
  (`rm -rf ~/projects/FLOWPanel.jl/data/p018_csarc_n2_nt72_l3p0_3r_srlx_nv_g25_nodefail13712097`).
- Notebook/ledger backlog additionally: T1a/T1b/T5b outcomes, the T1b
  ignition incident, and the wave-1 final verdicts. Awaiting approval.
