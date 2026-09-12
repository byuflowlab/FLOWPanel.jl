# RESET BRIEF — NT-divergence discriminators, worktree era (2026-09-05)

Supersedes `smooth_ladders_handoff_20260902.md` as entry point. Guardrails from
that file's §5 chain still apply (sacct authoritative, `ssh orc` needs
`bash -lc` + live ControlMaster, NODE_FAIL restarts are fresh, exit 124 →
`_s2` chain, don't touch `fp-il-*`/`fp052*` jobs). NEW binding rules:
**silos are dead** (Ryan 2026-09-05, memory `feedback_no_more_silos_worktrees`)
— every submission runs from a pinned git worktree; pins are **annotated tags**
not bare SHAs (Campaign Reproducibility section of `~/.claude/CLAUDE.md`).

## 0. Where the investigation stands

CT-vs-NT **anti-converges** and the per-doubling delta roughly **doubles**
each doubling (ratios Δ72→144/Δ36→72 = 1.68–1.92 across all five completed
arms) — the fingerprint of a fixed bias per STEP application integrated over
∝NT steps/rev, NOT truncation (would halve) and not simply a per-step-scaled
time constant (would stay constant). Formal CT_per_rev-window scores
(2026-09-05, all 30 revs, no zero-fill; match log-tails to ±0.07%):

| Arm (λ3.0, 3R) | NT36 | NT72 | NT144 | Δ36→72 | Δ72→144 |
|---|---|---|---|---|---|
| Reference | 0.070775 | 0.071844 | 0.073829 | +1.51% | +2.76% |
| Merge OFF | 0.070811 | 0.071934 | 0.074132 | +1.59% | +3.06% |
| Fixed rlxf=0.3 | 0.070775 | 0.071709 | 0.073300 | +1.32% | +2.22% |
| Traditional Pedrizzetti | 0.071706 | 0.072730 | 0.074594 | +1.43% | +2.56% |
| tp + rlxf=0.3 | — | 0.072424 | 0.074154 | — | +2.39% |

Merge, relaxation family, and relaxation dose are all cleared. λ converges.

## 1. LIVE jobs (submitted 2026-09-05; ladders SHORTENED to NT36+NT72 by Ryan)

All from the pin-`f46c3fe` worktrees (see §2). NT36→mgh, NT72→eng
(`--qos=eng --constraint=intel`; eng preemption privilege held). NT144 rungs
cancelled for turnaround.

| Arm | NT36 (mgh) | NT72 (eng) | Tests |
|---|---|---|---|
| 1 `_3r_sv` smooth conv, overlap 2.75 | 13592738 | ~~13592723~~ → 13593011 | unsteady interface Γ carrier |
| 2 `_3r_sv_h2p0` smooth, overlap 2.0 | ~~13592759~~ | ~~13592726~~ | **DROPPED** (Ryan 2026-09-05) |
| 3 `_3r_sfs3nb` legacy shed, `SFS_THREELEVEL=true` | 13592760 | ~~13592729~~ DEAD | SFS procedure/backscatter |
| 4 `_3r_sv_s1p5` smooth, 1.5σ | 13592731 **DONE 0.072222** | 13592732 | smoothing width |
| 5 `_3r_exp` `WAKE_EXPINT=true`, legacy shed | 13592894 | 13592895 | stretching-integrator stiffness |
| 6 `_3r_srlx` `SFS_RLXF=0.0025031` | = ref 13507289 (0.070775) | 13592896 | per-step SFS coefficient memory |

### CURRENT FLEET (2026-09-05 evening — authoritative over the table above)

| Job | Arm/rung | State at reset |
|---|---|---|
| 13592732 | arm 4 NT72 (s1p5) | RUNNING eng-1-1 (~5.5 h in, 14 h wall) |
| 13593011 | arm 1 NT72 (sv, resubmit) | pending eng (Priority) |
| 13592896 | arm 6 NT72 (srlx) | pending eng (Priority) |
| 13593711 | arm 3 NT36 (sfs3nb, resubmit after NODE_FAIL) | pending mgh |
| 13593717 | arm 5 NT36 (exp, GUARD-OFF resubmit) | pending mgh |
| 13593718 | arm 5 NT72 (exp, GUARD-OFF resubmit) | pending eng |

Scored: arm 4 NT36 = 0.072222 (windowed, 10 rev blocks). Reference NT36 =
0.070775, NT72 = 0.071844 (§0 table; the flatten-vs-climb test compares each
arm's NT36→72 delta to the reference +1.51%).

Rulings this session (Ryan): arm 2 (h2p0) DROPPED entirely; arm 1 = option
(c) hold — decide the NT36 retry after 13593011's outcome; arm 5 resubmitted
WITHOUT SIGMA_CEIL (euler_exp rejects sigma_guard; note guard-on/off confound
vs other arms when scoring). Monitor via hpc-monitor; judge by outputs.

Storage: hpc-storage's detached --apply archive pass was launched ~17:00 on
the cluster (logs on /home); one checkout skipped due to a stale lock from a
dead 2026-09-03 session (lock left for a human). Verify it completed.

Session log (2026-09-05, details behind the table above):
- 13592723 died at launch on the `isfile(DAS_ARC_TABLE)` check
  (`rotor_hover_pressure_comparison.jl:814`): started 13:22:57, the worktree
  `data` symlink landed 13:25 — launch-time race, not a bug. Resubmitted
  identical → **13593011**.
- Arm 2 NT72 (13592726) died of WakeGeometryError ~57 min in (~step 660);
  Ryan ruled: **drop the overlap-2.0 ladder**. CORRECTION: NT36 13592759
  was NOT cancelled-while-pending — it had already run 2h16 and died of
  WakeGeometryError at step 1076/1079 (99.6%, σ-growth 1.19 vs 1e9
  critical) before the scancel landed. Its `_CT_per_rev.csv` is likely
  scoreable without a rerun if the dropped arm is ever wanted as data.
- Arm 3 NT72 (13592729) died at step 1551/2159 (rev ~21.5, just before
  settle=22): FMM regularized-nearfield adequacy gate, `sigma_max=0.2036`
  vs `SIGMA_CEIL=0.030` — σ blow-up under SFS_THREELEVEL that the σ-guard
  ceil did NOT clamp (units mismatch or guard-bypassing growth path —
  UNRESOLVED). CF excursion −0.069→−0.117 on the final step. Not resubmitted
  (deterministic). Possibly a discriminator datum: reference smooth arm
  survived this rung.
- Arm 4 NT36 (13592731) COMPLETED: windowed mean CT = 0.072222 (10 rev
  blocks) = +2.0% level shift vs reference NT36; verdict waits on 13592732
  (RUNNING since ~19:00).
- Arm 1 NT36 (13592738) died at step 536/1079 (rev ~14.9, 1h23 elapsed):
  WakeGeometryError, "wake panel folded or inverted, normal reverses at
  (xi,eta)=(0,1)" — same open family as arm 2 NT72 + the five wave-3
  deaths. NOT resubmitted (presumed deterministic); NEEDS RYAN RULING.
  Note the smooth-conversion family is now 2-for-4 on WGE deaths this wave
  (arm1 NT36, arm2 NT72 dead; arm4 NT36 done, arm4 NT72 running) while
  legacy-shed arms have zero WGE so far (arm3 died of σ/FMM gate instead).
- SIGMA_CEIL confirmed in METERS (driver prints "SIGMA_CEIL=... m"), so
  arm 3's sigma_max=0.2036 m vs ceil 0.030 m is a real guard bypass, not
  a units mismatch.
- Arm 1 ruling (Ryan): option (c) — HOLD, wait for 13593011's outcome
  before deciding on an NT36 retry.
- Arm 3 NT36 (13592760) NODE_FAIL at 1:23 (step 816/1079, log healthy) —
  genuine node death; fresh resubmit per policy → **13593711**.
- Arm 5 NT36 (13592894) failed at startup in 4 min: "sigma_guard is not
  supported by the euler_exp integrator" (FLOWPanel_wake.jl:2227) — the
  wave-common SIGMA_CEIL=0.030 export conflicts with WAKE_EXPINT=true.
  NT72 twin 13592895 (pending) will die identically. NEEDS RYAN: rerun
  arm 5 with guard off (confound: other arms run guard-on), or drop arm.
- σ-guard log grep (suspect #3) is a DEAD END: the guard in
  `FLOWVPM_timeintegration.jl` clamps silently (no print on CPU or GPU
  path) — measuring activations needs instrumentation or post-hoc σ from
  snapshots.

As of last check: 13592731 RUNNING past first shed (~2 s/step — **GPU-path
risk RETIRED**: SurfaceVorticityConversion works on the CuArray pfield);
13592726 RUNNING on eng; others pending (mgh shares with an fp052c job).
Common exports every job: `SIGMA_CEIL=0.030,TRUNCATION_RADIUS_R=3.0,
MAX_PARTICLES=1500000,P018_SETTLE_REVS=22` + arm extras +
`P018_REPO_OVERRIDE=$HOME/wt018/FLOWPanel-pin-<arch>,
P018_PROJECT_OVERRIDE=$HOME/p018wtenv-<arch>,P018_RUN_NAME=<case>_<suffix>`;
payload `examples/run_dji9443_hover_ct_gpu.slurm.sh <arch> <case_tag>`; logs
`~/projects/FLOWPanel.jl/logs/slurm/` (absolute).

**Readout per arm**: `<run>_CT_per_rev.csv`, mean CT_mean over
`in_convergence_window==true` (fallback monitor02_force for restart carriers).
Decision: an arm whose NT36→72 climb flattens vs reference (+1.51%) fingers
its mechanism. Arm 6 is single-rung: if its NT72 lands near 0.070775 instead
of 0.071844, per-step SFS rlxf is the driver. If ALL still climb → per-step-
bias discriminators in §4.

## 2. Worktree/tag structure (replaces silos)

Cluster `~/wt018/`, envs `~/p018wtenv-*`. Running campaign pin = tag
`campaign/p018-smooth-ladders-20260905` (FLOWPanel `f46c3fe`, exact snapshot
of the 018-gpu silo code all completed arms ran + data-symlink commit; dep
silo snapshots tagged `...-{h200,gh200}` in FLOWVPM/FastMultipole repos —
deps genuinely differ per arch: FMM radix + CUDA translate kernels).
Full SHAs + saga (rk3-kwarg mirror disaster, data-symlink fix) in
`smooth_ladders_provenance_20260905.md`.

NEXT-launch base (prepped, UNUSED, not GPU-smoked): tag
`campaign/p018-nt-20260905` in all three repos = local FLOWPanel `0e08ab4`
(launcher constraint fix + `scripts/prep_campaign_worktree.sh`) + unified
FLOWVPM `3315b22` + FastMultipole `3da58a1a`. Worktrees
`~/wt018/{FLOWPanel-nt-h200,FLOWPanel-nt-gh200,FLOWVPM-nt,FastMultipole-nt}`,
envs `~/p018wtenv-nt-{h200,gh200}`. Caveats: smoke NT36 before trusting
walltimes; unified deps drop per-arch kernel patches (kills the cross-arch
confound going forward, but new gh200 runs ≠ bit-comparable to old mgh arms).
New-campaign recipe: `git show <tag>:scripts/prep_campaign_worktree.sh | bash
-s -- <tag> <dir>` after pushing the tag.

## 3. Suspect ranking (audit rounds 1 AND 2 done)

1. Per-step ops with fixed bias per application (matches doubling signature):
   FMM eval bias (tuned knobs' CT null only verified at one NT — 023/025),
   iterative-solver per-step tolerance, f32 accumulation in GPU path.
2. DynamicSFS per-step rlxf (`FLOWVPM_subfilterscale.jl:954-955`, memory time
   dt/0.005 shrinks 4× NT36→144) + backscatter clip duty cycle — being tested
   NOW by arms 3+6.
3. `sigma_guard` ceil (`FLOWVPM_timeintegration.jl:24-28`) documents a prior
   NT144-specific σ-compression blowup — grep completed arms' logs for guard
   activation counts by NT (cheap, not yet done).
4. Cross-arch confound: NT36 always ran ARM/gh200 with different FMM kernels
   than NT72/144 (x86) — one NT36 rung on h200 (~2.5 h) would size it.
5. Viscous CoreSpreading cleared (per-physical-time trigger). Shed
   quantization floor plausible-secondary (floor engages inboard at NT144).

Round-2 audit (2026-09-05) RULED OUT for the production config (velocity
formulation + Backslash + pinned FMM knobs, verified in driver defaults):
per-step iterative-solver tolerance (solver is direct Backslash; Krylov/FGS
`atol=1e-6`-per-solve WOULD match the doubling shape if ever used —
`src/FLOWPanel_solver.jl:1001-1002` — check 021-benchmark variants);
FMM autotune cadence (all autotune off, knobs static; fixed relative FMM
error is multiplicative → integrates NT-independent over a fixed rev, only
NT-coupled if the near/far split geometry shifts with particle spacing —
unprofiled); GPU Float32 (`FLOAT_TYPE=Float64` governs state; the f32 knob
is VTP output only — but f32 VTP **warm-starts** would inject it, flagged in
`FLOWPanel_warmstart.jl:263-264`); Das refresh (off), GreenReconstruction
staleness (velocity formulation doesn't call it), TE filament (exact, no
truncation), handoff blending (nwakerows=1 → inactive). Residual round-2
leads: verify all completed rungs truly ran identical solver/formulation
config (only driver defaults + ops_reference.md were checked), and profile
`FLOWVPM_timeintegration.jl`'s per-particle update for fixed non-dt-scaled
corrections (not yet done line-by-line).

## 4. Next actions

1. Monitor the 11 live jobs (hpc-monitor; judge by outputs not sacct state).
   NT36 ~2.5–8 h, NT72 ~5.4–14 h walltimes. On exit-124 → `_s2` chain; on
   NODE_FAIL → fresh resubmit (never chain).
2. Kick `hpc-storage` (VTK-writing jobs live; not yet kicked this wave).
3. Harvest per §1 readout when rungs land; extend the §0 table.
4. If all arms climb: backend-matched NT36 (h200) + conservative-FMM NT72
   discriminators; also the σ-guard log grep and metadata check for particle
   precision (f32 vs f64) in the GPU runs.
5. Owed to Ryan: notebook entry for the whole 09-05 arc (formal scores table,
   worktree migration, discriminator design — ask verbosity first); ledger.md
   still stale (pre-08-29); INDEX.md + item-file RESET BRIEF stale (~08-21);
   mgh-1-2 triple-NODE_FAIL unticketed.

## 5. Local repo state

`fastmultipole` branch at `0e08ab4` (clean tree at reset; tag
`campaign/p018-nt-20260905`). Local scratch worktrees `wt593`, `wt_*` under
the session scratchpad may linger — `git worktree prune` after deleting.
Cluster repo has local branch pushed as `p018-nt-local`. Formal-score files
in the 09-05 session scratchpad are disposable (scores are in §0 and the
provenance file).

## 6. Fleet update 2026-09-06 (context-reset addendum)

hpc-monitor check after the 13593711 completion event:

| Job | Arm | State | Notes |
|---|---|---|---|
| 13593711 | sfs3nb NT36 (resubmit) | COMPLETED clean | 1080/1080; cycle-mean CT 0.06938 ± 0.0001 (±0.14% / 10 revs); Phase-2e gate PASS rc=0 |
| 13592731 | s1p5 NT36 | COMPLETED clean | 1080/1080, converged (scored §0) |
| 13592760 | sfs3nb | RUNNING | 816/1079 |
| 13592732 | s1p5 NT72 | RUNNING | ~1838/2159 |
| 13592894 | exp | FAILED step 1 | `sigma_guard is not supported by euler_exp` |
| 13593717 | exp (GUARD-OFF resubmit) | FAILED step 3 | **same sigma_guard error** — the SIGMA_CEIL drop never took |
| 13592759 | sv-h2p0 | FAILED 1076/1080 | GPU dispatch error after convergence window; per-rev CSV likely scoreable |

**First action on resume: re-launch the exp arms without sigma_guard, for real
this time.** 13593717 was launched as the guard-off resubmit yet died on the
identical `ArgumentError`, so the `SIGMA_CEIL=0.030` export is still reaching
the euler_exp dispatcher — check the GPU launcher / wave env script and confirm
(e.g. grep the job's `.err`/env dump) that `SIGMA_CEIL` is unset/empty for
`WAKE_EXPINT=true` arms *before* submitting. Remember the guard-on/off confound
note (§ triage): the expint arm runs unguarded relative to the rest of the
ladder — Ryan approved this.
