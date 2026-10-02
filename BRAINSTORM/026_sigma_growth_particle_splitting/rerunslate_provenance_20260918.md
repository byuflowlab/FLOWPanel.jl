# 026 §22.3 rerun slate — provenance (2026-09-18)

Ryan-approved 2026-09-18 (AskUserQuestion): 5-arm slate; NT144 cap-ladder
rung PARKED until this slate reads out.

## Purpose

A/B the §22 merge-σ redesign (second-moment merged σ + merge-path
lineage, incl. axis ruling) against the wave-2 s9 deaths, and validate
caps as a non-binding tripwire rather than a pump compensator. Also the
first live-sim exercise of the axis-lineage change (unit-tested only;
smoke was pre-axis).

## Pins (annotated tag `campaign/p026-rerunslate-20260918`)

| repo | commit | note |
|---|---|---|
| FLOWPanel.jl (`fastmultipole`) | `6e52628` | 026 doc bundle `8eff54e` + RUN_NAME_OVERRIDE launcher hook |
| FLOWVPM.jl (`flowpanel`) | `8d4a3b4` | §22 merge redesign (σ law + merge lineage + axis) |
| FastMultipole | `ac7230a6` | unchanged; carries `campaign/p026-derisk-20260914` |

Worktrees: `orc:~/campaigns/p026-derisk-20260914/{FLOWPanel.jl,FLOWVPM.jl,FastMultipole}`
checked out detached at the tags above, `dirty=0` verified 2026-09-18.
Env: `orc:~/campaigns/p026-derisk-20260914/env` — Manifest dev-paths
verified pointing at the three campaign worktrees (unchanged from wave-2).
Tags pushed to `orc:projects/{FLOWPanel.jl,FLOWVPM.jl}` by name (branches
NOT pushed — divergence ruling still parked with Ryan).

## Slate

All s9 wave-2 case defs unchanged; new code generation distinguished by
`RUN_NAME_OVERRIDE` (`scr_p026s9r2_*` run dirs — no harvest collision
with wave-2). FLOOR = `SIGMA_FLOOR_FRAC=0.1`; SPLIT =
`WAKE_SPLIT_VISCOUS=true WAKE_SPLIT_FRAC_VISCOUS=0.587
WAKE_SPLIT_FRAC_COMPRESS=0.73 WAKE_SPLIT_FRAC_ELONGATE=0.3`. Wall 12 h
(wave-2 ctrl fs/split hit the 8 h wall at 355/370). H200 via
`run_p018_screen_gpu052.slurm.sh h200 <case>`, cwd = campaign FLOWPanel
worktree, `P018_REPO_OVERRIDE`/`P018_PROJECT_OVERRIDE` at the campaign
worktree/env.

| arm | case | run name | mech env | extra | wave-2 death it A/Bs |
|---|---|---|---|---|---|
| A1 | scr_p026s9_explg_fs | scr_p026s9r2_explg_fs | FLOOR+SPLIT | — | overflow @295 (np 497k) |
| A2 | scr_p026s9_explg_floor | scr_p026s9r2_explg_floor | FLOOR | — | DomainError dt·\|L\|=2133 @274 |
| A3 | scr_p026s9_ctrllg_fs | scr_p026s9r2_ctrllg_fs | FLOOR+SPLIT | — | \|CFx\|>1 @259 (ctrl channel, pump-independence probe) |
| A4 | scr_p026s9_explg_fs | scr_p026s9r2_explg_fs_cap3 | FLOOR+SPLIT | `WAKE_SPLIT_SIGMA_MAX=0.0071 SIGMA_CEIL=0.0071` (k=3 × shed σ 0.00237; adequacy ≈0.008) | — (tripwire vs A1) |
| A5 | scr_p026s9_explg_split | scr_p026s9r2_explg_split | SPLIT | — | DomainError @280 (replicate of A2 class) |

## Acceptance / read-out

- A1/A2/A5 complete 467/467 with σ_max staying near smoke behavior
  (~1.1× shed σ, i.e. ≪ 2.5–3.4× wave-2 death band) → pump-attribution
  confirmed; caps not needed as compensators on s9.
- A4: zero σ_max/SIGMA_CEIL clamp engagements is the pass criterion
  (cap validated as non-binding tripwire); any binding → k=3 too tight
  under new law, report.
- A3: death in the ~260–290 window again → ctrl-channel confirmed
  pump-independent (separate investigation); completion → wave-2 ctrl
  blow-up was pump-coupled after all.
- Banner verification mandatory per wave-2 checklist (reg linegauss,
  mechanism knobs per arm, merge-event log + sigma telemetry lines).
- Judge by outputs, not sacct.

## Submissions (2026-09-18, cwd = campaign FLOWPanel worktree, wall 12 h, h200/m13h)

| job | arm | case | RUN_NAME_OVERRIDE |
|---|---|---|---|
| 13758586 | A1 | scr_p026s9_explg_fs | scr_p026s9r2_explg_fs |
| ~~13758587~~ → 13763819 | A2 | scr_p026s9_explg_floor | scr_p026s9r2_explg_floor |
| ~~13758588~~ → 13763820 | A3 | scr_p026s9_ctrllg_fs | scr_p026s9r2_ctrllg_fs |
| ~~13758589~~ → 13763821 | A4 | scr_p026s9_explg_fs | scr_p026s9r2_explg_fs_cap3 |
| ~~13758590~~ → 13763822 | A5 | scr_p026s9_explg_split | scr_p026s9r2_explg_split |

Mechanism env exactly as the Slate table (FLOOR/SPLIT strings verbatim in
the submission `--export=ALL,...`). Slurm `.out` files land in the
campaign FLOWPanel worktree under `logs/slurm/slurm-fp-052-scr-gpu-<jobid>.out`.
Banner verification owed on start (wave-2 checklist + `RUN_NAME` = the r2
override in the "Artifacts:" line).

### eng requeue (2026-09-18 ~18:15, Ryan: "we get priority on eng")

m13h start estimates for A3-A5 were 05:30-09:00 next morning; `sbatch
--test-only` on eng said 18:20 same day (eng-1-1: 7 idle H200). Cancelled
the four PENDING arms 13758587-90 (no compute lost) and resubmitted
identically but `-p eng --qos=eng --gres=gpu:h200:1` (wall/env/case/run
names unchanged): A2=13763819, A3=13763820, A4=13763821, A5=13763822.
A2 was RUNNING on eng-1-1 within a minute. Note: each job asks 64 CPUs,
so eng-1-1 serializes the arms; still hours ahead of m13h. A1 stays on
m13h (13758586, already running) — same H200 model, backend-matched.

## Outcomes

| job | arm | outcome |
|---|---|---|
| 13758586 | A1 explg_fs uncapped | **DIED step 328/467** (2026-09-18 18:24): DomainError dt·\|L\|=3461.6 exceeds euler_exp substep budget (`FLOWVPM_timeintegration.jl:713` via `propagate!`). Same death CLASS as wave-2 explg_floor (2133 @274) / explg_split (2229 @280); wave-2 explg_fs itself died of 500k overflow @295. Redesign DELAYED death (+33 steps vs explg_fs, +48–54 vs the DomainError arms) but did not prevent exp-family ignition. Precursors: radix FMM sigma adequacy warnings steps 310–328 (ratio 0.999→0.984→0.919; **adequacy limit had shrunk to 0.004463**, count-driven), merge activity ramping (87,287 merge_events rows, 12+/step late). No NaN, no overflow, GPU mem stable. No CT_vs_rev.csv (checkpoint never fired); force monitor stops at 327. Data: campaign worktree `data/scr_p026s9r2_explg_fs/`. Early read: σ-pump removal is NOT sufficient for exp-family s9 survival — consistent with the Γ-side ignition channel; A2/A5 (floor-only / split-only) will discriminate. |

| 13763819 | A2 explg_floor uncapped | **DIED step 274/467** (2026-09-18, elapsed 37 min): DomainError dt·\|L\|=5251.8 exceeds euler_exp substep budget — **exactly the same step as the wave-2 twin (2133 @274), larger magnitude**. Zero delay from the merge redesign on the floor-only arm ⇒ the Γ-side ignition channel in the floor arm is fully independent of the merge σ-pump; A1's +33-step delay presumably came from the split interaction. No adequacy warnings, no NaN, sigma growth nominal 1.195 at death, 2754 merge events. GOTCHA: the traceback lands in `logs/slurm/slurm-fp-052-scr-gpu-<jobid>.err`, NOT the .out (the .out ends with the GATE line, `dispatcher_rc=1`) — first-pass autopsy of the .out alone mislabels this a launcher failure. |

| 13763820 | A3 ctrllg_fs | **COMPLETED 467/467 but PHYSICALLY BLEW UP (2026-09-19 05:30)** — sacct COMPLETED is misleading (judge by outputs): CT_per_rev shows onset in rev block 8 = steps 253–288 (CT_ptp 1535, escalating to 4272 in block 9), i.e. **the same ~260–290 window as the wave-2 ctrl \|CFx\|>1 death @259**; final readout CT −49.4 ± 23.9%, CONVERGED=false, GATE technically clean (no NaN, gpu_gemv=468). Unlike wave-2 the run never tripped a fatal guard and marched to 467. **Slate read-out: ctrl-channel blow-up recurs in-window under the merge redesign → ctrl channel confirmed pump-INDEPENDENT** (separate investigation per acceptance). |

| 13763821 | A4 explg_fs_cap3 | **DIED step ~316/467 (2026-09-19): PARTICLE OVERFLOW (500k cap)** — different terminal symptom than A1 (which RODE the cap from ~318 with merge maintenance holding, then DomainError @328). np 497.6k and climbing at step 315 (wake-health), max-col telemetry surging 460→590 — same-class late-run particle surge. Merge events 62,421 through ~315 vs A1's 87,287 through 327. Read: with σ_max=SIGMA_CEIL=0.0071 the maintenance could NOT hold np at the cap where A1's uncapped merging could — i.e. **the k=3 cap appears to BIND through the merge law** (per acceptance: "any binding → k=3 too tight under new law, report"). Caveat: no per-event clamp-engagement telemetry exists (known gap), so binding is inferred from the overflow-vs-ride contrast, not counted directly. |

| 13763822 | A5 explg_split | **DIED step ~306/467 (2026-09-19): DomainError dt·\|L\|=5.59e7** — same euler_exp guard class as A1/A2, +26 steps vs its wave-2 twin (@280). Wake-health at 305: np 447k, extreme tail ignited in ONE step (max Γ/σ² 889→3.08e7, min σ collapsing 1.6e-4→2.9e-5 with floor OFF — this arm has no σ floor). Classic Γ-ignition, fastest-blowing variant (no floor to arrest σ collapse). |

### SLATE VERDICT (all five arms terminal, 2026-09-19)

**The §22 merge-σ redesign does NOT rescue any s9 arm.** All five died or
blew up, four of them in/near the ~274–328 band: A1 DomainError @328, A2
DomainError @274 (exact wave-2 twin step), A3 physical blow-up onset
253–288 (sacct-COMPLETED with CT −49.4), A4 overflow @~316, A5 DomainError
@~306. Read-outs vs acceptance:

- **Pump attribution REFUTED as the death driver**: σ_max stayed nominal
  where measured; deaths are Γ-side ignition (localized, fountain-region)
  in every exp arm and the ctrl channel — all **pump-independent**.
- **A4 (cap3): the cap appears to BIND through the merge law** (overflow
  where A1's uncapped merging held np at the cap) → k=3 too tight under
  the new σ law, per the "any binding" clause.
- **The productive lever is BRAINSTORM 032**: same A1 config with
  `PARTICLE_OMIT_ROOT_R_OVER_R=0.15/0.20` completed 467/467 converged
  (CT 0.0716/0.0713, clip-width-insensitive) — root-shed circulation is
  the fountain Γ-ignition driver (`032_rerun_provenance_20260918.md`).
  B-floor (A2-config + omission, 13772015) running as of this entry.
- Candidate follow-up for Ryan: ctrl-channel + omission discriminator
  (out of current directive scope).

### Harvest (2026-09-19; run dirs MOVEd to `orc:~/projects/FLOWPanel.jl/data/` + symlinks back at the worktree paths; VTK preserved for 032 comparison; archiving deferred)

| arm | last step | final np | max Γ/σ² (last 10 rows) | min σ (last 10) | CT | size |
|-----|-----------|----------|-------------------------|-----------------|----|------|
| A1 | 327 | 499,377 | 9.08e4 | 2.37e-4 (floor) | no readout | 11 G |
| A2 | 273 | 266,288 | 1.44e6 | 2.41e-4 | no readout | 5.2 G |
| A3 | 467 | **17,169** | **4.69e9** | 2.37e-4 | −49.4 ± 24% (per-rev CSV; harvester missed it) | 6.3 G |
| A4 | 315 | 497,589 | 590 (overflow, not ignition-terminal) | 2.37e-4 | no readout | 9.3 G |
| A5 | 305 | 447,320 | 3.08e7 | 2.94e-5 (no floor) | no readout | 8.5 G |

A3 detail: np COLLAPSED to 17k by the end (from ~250k-class mid-run) with
max Γ/σ² 4.69e9 — the same merge-frenzy/self-destruction signature as
wave-2's ctrllg_floor (np 246k→54k). Total slate footprint 40.3 G.

### Death-character correction (Ryan ParaView review, 2026-09-18 evening)

Ryan inspected A1 steps 278–327 in ParaView: field looks fine, a surge
of particles near the root in fountain flow, "not fully unstable (yet)"
— and the telemetry agrees. "Blow-up" was an overstatement:

- The DomainError is a **numerical tractability guard**, max-over-
  particles: euler_exp subdivides to θ=dt·\|L\|≤0.5/substep with
  EXP_MAX_SUBSTEPS=4096, so it throws at worst-particle dt·\|L\|>2048
  (`FLOWVPM_timeintegration.jl:647-713`; the elementwise bound also
  over-estimates \|L\| ≤3×). One bad close pair trips it while the rest
  of the field is healthy.
- A1 monitor04 tail: global metrics stayed sane (mean σ flat 0.0019,
  forces nominal, np smooth) while the EXTREME TAIL ran away over the
  last ~15 steps: min_σ pinned at the viscous floor (2.374e-4) from
  ~315, max Γ/σ² 826→9.1e4, max_u 49→919, max_dtZ 1.3→107 (steps
  315→327). Classic **localized Γ-side ignition at the σ floor** (root/
  fountain region per ParaView), not field-wide instability.
- A1 also **rode the 500k particle cap from ~step 318** (np pinned
  ~499k, maintenance holding) — wave-2 explg_fs died OF that cap @295;
  the new merge law kept the run alive at the ceiling.
- Same refinement applies to A2 (σ-growth nominal 1.195 at death).
  Whether either run would have gone globally unstable is untested —
  the guard terminates first, by design.

### A1 offending-particle localization (step-327 VTP analysis, 2026-09-18)

Frame: rotor axis = x (`omega_axis=[-1,0,0]`, thrust −x, wake convects
+x to ~+2.8R; radial = √(y²+z²)). Per-particle guard statistic
θ = dt·Σ|J| recomputed from the saved `velocity_gradient` (dt=3.086e-4):
max θ at step 327 = **3461.6 — exactly the DomainError value**, so the
step-328 throw used this state.

- **Offender: depth x = −0.37R (ABOVE the rotor plane, thrust side),
  radial r = 0.344R** (x=−0.0437, y=0.0136, z=+0.0386 m; σ=2.99e-4).
- The entire top-10 θ cluster is a ~1 mm knot at depth −0.36 to −0.46R,
  radial 0.34–0.38R — the fountain-flow recirculation region Ryan
  identified in ParaView.
- Driver vs victims: the guard-tripping particles carry tiny \|Γ\|
  (1.7e-5–2.3e-5) — they are close-by VICTIMS of the max-Γ/σ² particle
  at the same spot (\|Γ\|=5.38e-3 at σ=2.43e-4 ≈ floor, Γ/σ²=9.1e4),
  0.44 mm = 1.47σ from the offender. Close-pair singularity: floor-σ,
  Γ-concentrated particle imposing a huge local gradient on neighbors.
- Fountain context: 29,611 particles (5.9%) sit above the rotor plane;
  12,146 of them in the 0.25–0.45R radial band.
- Tooling note: step VTPs carry `velocity_gradient` (9 comp) — the
  guard statistic is exactly reconstructible offline; parse the raw
  appended VTP binary directly (meshio cannot read .vtp).

## Banner verification

| job | arm | verdict |
|---|---|---|
| 13758586 | A1 | **PASS 9/9** (2026-09-18): 12 revs/468 steps, linegauss, floor 0.1 guard=on, f_visc/f_comp/f_elong = 0.587/0.73/0.3, expint on, SIGMA_CEIL=Inf, merge log + sigma telemetry at `data/scr_p026s9r2_explg_fs/`, no errors; step 277/467 at check — already past the wave-2 exp fs death step (295 was the overflow death) with no overflow signs |
| 13763819 | A2 | **PASS 9/9** (2026-09-18, eng-1-1): 12 revs/468 steps, linegauss, floor 0.1, split knobs OFF as intended (floor-only arm), expint on, SIGMA_CEIL=Inf, merge log + telemetry at `data/scr_p026s9r2_explg_floor/`, no errors through step 17 (spinup) |
| 13763820 | A3 | **PASS 9/9** (2026-09-18, eng-1-1, started ~18:58): 12 revs/468 steps, linegauss, floor 0.1 guard=on, split trio 0.587/0.73/0.3, WAKE_EXPINT=false (ctrl arm — correct), SIGMA_CEIL=Inf, merge log + telemetry at `data/scr_p026s9r2_ctrllg_fs/`, no errors at step 0 (das_eta:nan is benign/expected) |
| 13763821 | A4 | **PASS** (2026-09-19 ~06:15, eng): `SIGMA_CEIL=0.0071 m (guard=on)` and split `clamp=[0.0002374, 0.0071]` — cap3 knobs confirmed; split trio 0.587/0.73/0.3, expint on, floor 0.1, CONVERSION=legacy, linegauss. Prediction on record: cap does NOT bind the Γ-driver (its σ sits at the FLOOR) — watch for a same-class death ~274–328. |
| 13763822 | A5 | **PASS** (2026-09-19, eng): split-only confirmed — split trio 0.587/0.73/0.3 ACTIVE with `clamp=[off, off]`, `SIGMA_FLOOR_FRAC=0.0 (floor=-Inf)`, SIGMA_CEIL=Inf, expint on, CONVERSION=legacy. Watch the ~274–280 window (wave-2 twin died @280). |
