# RESET PROMPT — 018 NT ladders: harvest, then investigate SFS Cd vs NT (2026-09-08)

You are picking up BRAINSTORM item 018 (DJI 9443 hover CT convergence campaign).
**Entry point: this file.** Then read, in order:
`nt_discriminators_handoff_20260906.md` §§14–17 (recent sessions), §§11–12
(recipes + gotchas learned in anger), and for the science background the sibling
`nt_discriminators_handoff_20260905.md` §0 (scores) / §2 (worktrees) / §3
(suspect ranking). Do NOT read whole BRAINSTORM files inline — use the
`brainstorm-scout` subagent.

## Ryan's task for this session, in order

1. **Score the two completed rungs** — 13605974 (exp-nt NT72 on-model) and
   13605983 (Ladder S Cs=0.18 NT36). Windowed CT, rev coverage stated, the
   exp-nt NT36→NT72 climb vs the reference +1.51%, logged in the handoff. Offer
   Ryan a notebook entry (ask verbosity first — never write to the notebook
   unasked).
2. **Investigate the SFS machinery**: determine whether NT has a clear,
   accountable impact on the SFS model coefficient Cd — and if so, what the
   campaign should do about it. Detail in §"The investigation" below.
3. **Diagnose the Ladder C blow-up mechanism** — is σ running away, where, and
   is it the reformulated-VPM stretching term acting in an anti-stretching
   region? See the ADDENDUM at the end of this file, which is Ryan's own framing
   of the question and the most specific brief here.

## Standing guardrails (violating these has cost this campaign whole days)

- **sacct state is NOT evidence.** Judge every run by its outputs. `FAILED` is
  routinely a controlled launcher exit with clean physics; `COMPLETED` can hide
  an off-model run.
- **Read the `.err` before the `.out`.** The GPU launcher writes them separately;
  every Julia stacktrace and the `ERROR: dispatcher exited rc=` line go to `.err`.
- **`ssh orc` needs `bash -lc` + a live ControlMaster** (2FA otherwise).
- **NODE_FAIL → fresh resubmit, never a chain.**
- **Worktrees + annotated tags only; no silos.** Never edit a worktree that has a
  queued or running job. Campaign worktrees must carry no uncommitted tracked state.
- **Another agent may share the local checkout and the queue.** Re-check
  `git status` and `squeue` before acting; never `git checkout .` or `git stash`.
- **No further 018 submissions** beyond what Ryan authorises. Ryan's standing
  remedy for the σ blow-up is 026 particle splitting (needs 026 FLOWPanel commits
  4–6 wired + a fresh pin); the σ-cap and all-direct-fallback remedies are PARKED.
- Local jobs: never more than 4 threads (HPC exempt).

## Scoring recipe

`scripts/p018_harvest_ct.py` — windowed CT = mean over revs with
`in_convergence_window == true` in `*_CT_per_rev.csv`; force-monitor fallback
(`monitors/*_monitor02_force_system1.csv`) for runs killed before the end-of-run
CSV. **References: NT36 0.070775, NT72 0.071844** (reference-family climb +1.51%).

⚠ **Window mismatch is a live trap.** A run that dies before rev 30 makes the
script fall back to a *different, less-settled* rev span. Comparing that against a
rung scored on revs 21–30 measures the window shift, not the physics — this
already produced a bogus +2.03% climb on 2026-09-08 (§14). Always report the rev
coverage with the number, and refuse to compute a climb across mismatched windows.

⚠ **Verify CSV timestamps against the job** before scoring: several data dirs mix
old and new runs, and the base-case dirs were clobbered on 09-07.

## Verified-good numbers so far (do not re-derive)

| Arm | NT36 | NT72 | Climb | Verdict |
|---|---|---|---|---|
| reference family | 0.070775 | 0.071844 | +1.51% | baseline |
| `_3r_sv_s1p5` (smoothing width) | 0.072222 (+2.00%) | 0.073285 (+2.01%) | +1.47% | **NULL** |
| `_3r_srlx` (SFS relax 0.0025031) | — | 0.071429 (−0.58%) | — | partial (~40% of a doubling) |
| `_3r_sfs3nb` (three-level SFS) | 0.069382 (−1.97%) | died @1550 | — | NT72 owed |
| `_3r_sfs3nb_m2` (pairs-only merge) | 0.069260 (−2.14%) | died @1550 | — | merge chaining ≠ sole σ source |
| `_3r_nosfs` (SFS off) | 0.069431 (−1.90%) | blew up @rev 20.24 | — | NT72 will not be rerun (Ryan) |
| `_3r_exp_nt` (stretching integrator) | 0.070536 (−0.34%) | see fleet table | — | provisional null |

## Fleet at reset — ALL TERMINAL, nothing running, nothing pending

No 018 jobs are in the queue. Full triage is in handoff §17; summary:

| Job | Arm | Status | What to do |
|---|---|---|---|
| **13605974** | exp-nt NT72 on-model | COMPLETE 2159/2159, gate_rc=0 | **SCORE IT.** Pairs with NT36 0.070536 (13603853) → the formal stretching-integrator slope. Both on-model, same pin lineage, both full-window — a clean comparison |
| **13605983** | Ladder S `_3r_cs0p18` NT36 | COMPLETE 1079/1079, gate_rc=0 | **SCORE IT.** First static-Cd rung; compare vs ref NT36 0.070775 and vs dynamic-Cd s1p5-family |
| 13605984 | Ladder S `_3r_cs0p18` NT72 | DIED rev 26.06 — FMM adequacy gate | Partial harvest possible via force-monitor fallback (reaches rev 26 > window open ~rev 20), but **window-mismatched** — do not use it for a slope |
| 13605985 | Ladder C `_3r_cs0p002_exp_nt` NT36 | DIED rev 24.9 — force blow-up (`DomainError` 8.6e9) | Not scoreable; and see the confound below |
| 13605986 | Ladder C `_3r_cs0p34_exp_nt` NT36 | DIED rev 10.0 — force blow-up | Not scoreable |

Run dirs are under `~/projects/FLOWPanel.jl/data/<run name>/`; logs under each
worktree's `logs/slurm/` (worktree per job listed in handoff §16).

**⚠ Ladder C is confounded — say so before anyone reads it as a result.** It needs
`WAKE_EXPINT=true`, and euler_exp rejects `sigma_guard`, so both rungs ran
UNGUARDED. Ladder S NT36 used the same static-Cs mechanism *with* the guard and
completed fine. So the Ladder C blow-ups cannot be attributed to freezing Cd
rather than to the missing σ ceiling. If Ryan wants the Cs sensitivity, it must be
rerun on the WAVE stack (guarded, no expint) against a same-stack dynamic control.

## The investigation: does NT have a clear, accountable impact on Cd?

This is Ryan's second task and the substantive one. What is already established
(handoff §15, do NOT re-derive):

- Cd is the **dynamic** model: per-particle, recomputed every step by
  `dynamicprocedure_pseudo3level_*` in `FLOWVPM_subfilterscale_models.jl`
  (~lines 744–1024). Stored at `C_INDEX = 37:39` (`C[1]`=coefficient,
  `C[2]/C[3]`=Lagrangian nume/deno), `FLOWVPM_particlefield.jl:485`.
- It is **already written to every wake `.vtp`** by
  `FLOWPanel_wake.jl:2427` — no new instrumentation needed. No monitor
  aggregates it, so it exists only as raw per-particle data.
- Extractors written and working: `~/vtp_C.py` and `~/vtp_C2.py` on the cluster
  (also in this session's scratchpad). They parse raw-appended VTK XML directly
  — **meshio cannot read these files**. `vtp_C2.py` compares at **equal
  revolution** (rev = step/NT), which is essential: NT36 step 1079 and NT72 step
  1550 are rev 30.0 and rev 21.5.
- **The distribution is bimodal**: ~50–55% of particles carry C exactly 0
  (`clipping_backscatter` zeroes them), so median-over-all is 0.0 everywhere.
  Use `mC` (mean over all active) and `mC+` (mean over the C>0 subset).
  Saturation at the `maxC=1.0` clamp is rare (≤0.6%), so maxC is not binding.
- **Cd already looks NT-dependent**: merge2 pair at matched revs 4/8/12/16 gives
  NT72/NT36 ratios 1.33, 1.40, 1.17, 1.29; `_3r_sv_s1p5` at rev 30 gives +14.6%
  (mC+ +19.9%). **The 2.1–2.4× at revs 20–21.4 is CONFOUNDED** — the NT72 merge2
  rung is entering its σ blow-up there. ~1.2–1.4× is the trustworthy figure.
- The three-level preset (`SFS_THREELEVEL`, α=0.667) runs Cd ~**200× smaller**
  than the default two-level `DynamicSFS` (α=0.999): 0.0008 vs 0.158 at matched
  rev. That is a far larger model difference than anything else on the ladder.

### What the investigation still needs to answer

1. **Is the NT dependence of Cd a real property of the dynamic procedure, or an
   artifact of the wake state differing between rungs?** Read the actual formula
   (`nume = (Γ·(Γ·∇)dU/dσ)(3α−2)`, `deno = (Γ·dSFS/dσ)/(ζ₀/σ³)`, Lagrangian
   relaxation by `rlxf`) and work out its explicit and implicit NT dependence:
   what in it scales with Δt, with particles-shed-per-step, with helix spacing,
   with σ? A derivation is worth more here than more runs.
2. **Is the `rlxf` Lagrangian averaging the carrier?** `rlxf` is a per-STEP
   relaxation, so at fixed physical time a doubled NT applies it twice as often —
   the averaging window in *physical time* halves with NT. That is a concrete,
   checkable NT-coupling candidate, and it is consistent with `_3r_srlx`
   (rlxf 0.005→0.0025031, i.e. ~halved) landing at −0.58%, ~40% of one doubling
   step (handoff §8). **If rlxf is the carrier, the fix is to make rlxf scale with
   Δt** so the physical-time averaging window is NT-invariant. Check whether that
   exactly reproduces the srlx result.
3. **Quantify cleanly, without the blow-up confound.** Use only rev ranges where
   both rungs are healthy (σ-growth flat, CF flat). Prefer arms with full VTK
   retained: `_3r_sfs3nb_m2` (1080 / 1550 files) and `_3r_exp_nt` NT36 (1080).
   Note most other runs kept only the newest ~288 snapshots, so early revs are
   gone.
4. **Then: what should the campaign account for?** Options to evaluate, not
   presume — a Δt-scaled `rlxf`; quoting Cd with its rev and arm rather than as a
   constant; or treating Cd's NT drift as a term in the CT error budget.

### Useful extra data now available
13605974 (exp-nt NT72, complete, VTK on) and 13605983 (Cs=0.18 NT36, complete)
both have fresh particle output. 13605983 is a **constant-Cd control**: its `C`
field should be exactly 0.18 on unclipped particles and 0 elsewhere — a good
end-to-end validation of the extractor and of the clipping-fraction assumption
(expect fr0 ≈ 0.5 if the ~50% clipping rate is genuinely Cs-independent).

## Ryan's open decisions (do not act on these unprompted)
- WGE-arm reruns and the `_3r_sv` resubmit (unblocked by the s1p5 null since 09-07).
- Where a merge-result σ cap lands in FLOWVPM (pair merges still uncapped at ×1.26).
- 026 particle-splitting wiring (FLOWPanel commits 4–6) — his designated remedy
  for the σ blow-up; needs a fresh pin before any `_m2s`-style rung.
- Whether to rerun the Cs sweep guarded on the wave stack (see the Ladder C confound).
- **Disk: `/home/rander39` is at 430 G against the 400 G cap.** A strip belongs to
  Ryan's OTHER session — do not start a second archiver without checking with him.
- **Owed**: the base-case reference CSVs (0.070775 / 0.071844) were clobbered on
  09-07 and survive only in handoff tables; recover them from the archive tarballs.
- **Notebook**: entries owed for the whole 09-07/09-08 arc (Cd finding, both new
  ladders, exp-nt slope). ASK RYAN for verbosity before writing; append-only,
  under a `# YYYYMMDD` header he has approved.

## Local repo state at reset
Branch `fastmultipole`. Commit **720947a** (mine, this session): `SFS_CONST_CS`
branch + `SFS=` diagnostics label in
`examples/rotor_hover_pressure_comparison.jl`.
Uncommitted/untracked in the shared checkout include other agents' work on
BRAINSTORM 026 and 018 handoff files — **do not `git checkout .` or `git stash`**.

---

# ADDENDUM (Ryan, 2026-09-08): investigate the Ladder C blow-up mechanism

Ryan wants the Ladder C deaths diagnosed, not just labelled as confounded. Both
tasks stand: **score the two completed rungs (13605974, 13605983) first**, then
this. Treat this as a forensic study of existing output — it needs no new runs.

## The hypothesis chain to test, in order

1. **Is σ getting large before the blow-up?** If yes —
2. **Where in the field?** Which particles: tip vortex, root, near-disk, far
   wake, the young sheet or the old wake? Spatial and age localisation both matter.
3. **Is the σ growth driven by the reformulated-VPM stretching term?** Read the
   actual σ ODE rather than assuming its form — the reformulated-VPM σ dynamics
   and the guard clamp live in `FLOWVPM_timeintegration.jl` (the clamp is around
   lines 327–330 and it SKIPS static particles; the euler_exp path is separate,
   lines ~33–34/108, and rejects `sigma_guard` entirely, which is why Ladder C
   ran unguarded).
4. **Does it make physical sense?** Ryan's specific prediction: for the stretching
   term to *grow* σ, those particles must sit in a region of **anti-stretching** —
   vorticity aligned with a negative velocity gradient along its own direction,
   i.e. axial compression of the vortex tube. Test it directly: compute the
   normalised stretching rate along the vortex axis,
   $$ S = \hat{\Gamma} \cdot (\nabla u) \cdot \hat{\Gamma} $$
   and check whether the large-σ population is concentrated in $S < 0$.
5. **If that holds, 026 particle splitting is a candidate fix** — splitting
   large-σ particles back toward the target resolution attacks the σ source
   directly, which is exactly Ryan's standing remedy for the step-1550/1876
   adequacy deaths. A confirmed anti-stretching mechanism would tie the two
   failure modes together and strengthen the case for 026. **Do not implement or
   submit anything** on the strength of this — report the finding and let Ryan rule.

## The data is already there — no instrumentation needed

Both Ladder C runs wrote **every step** right up to the death, under
`~/projects/FLOWPanel.jl/data/`:

| Run | Snapshots | Last step | Death |
|---|---|---|---|
| `p018_csarc_l3p0_3r_cs0p002_exp_nt` | 899 | 898 | rev 24.9, `DomainError` 8.6e9 |
| `p018_csarc_l3p0_3r_cs0p34_exp_nt` | 362 | 361 | rev 10.0, `DomainError` 8.2e4 |

Neither wrote a `CT_per_rev.csv`, but both have
`monitors/*_monitor02_force_system1.csv` — use that for the CF time series to
locate the exact onset step, then work backwards through the VTPs from there.

**Every field you need is already in the `.vtp` point data**: `Points`, `gamma`,
`sigma`, `vol`, `circulation`, `velocity`, `vorticity`, `C`, `SFS`, and crucially
**`velocity_gradient`** — so $S = \hat{\Gamma}\cdot\nabla u\cdot\hat{\Gamma}$ is
computable directly from saved output. Extend `~/vtp_C.py` (its `read_vtp(path,
want)` takes a set of array names and already handles the 3- and 9-component
arrays; confirm `velocity_gradient`'s component count and row/column convention
before trusting any sign — **a sign error here inverts the entire conclusion**).

## Suggested shape of the analysis

- σ time series: percentiles (median, p95, p99, max) vs step for both runs, and
  the count above the wave ceiling 0.030. Establish whether σ runs away before,
  with, or after the CF blow-up — the ordering is the whole question.
- Spatial: for the top-σ percentile at a few pre-death steps, report radial and
  axial position (rotor radius R ≈ 0.119 m) and particle age if recoverable.
- Stretching: joint distribution of σ (or dσ) against $S$; specifically, what
  fraction of the high-σ tail sits in $S<0$ versus the population baseline.
- **Control**: do the same on `p018_csarc_l3p0_3r_cs0p18` (13605983, 1080
  snapshots, guarded, COMPLETED, never blew up). This is the decisive comparison
  — same static-Cs mechanism, guard on, healthy. If its σ tail also sits in
  anti-stretching regions but stays bounded, the guard is what's holding it and
  the anti-stretching is ordinary; if only the unguarded runs develop the tail,
  the guard is masking a real σ source. Also available:
  `p018_csarc_n2_nt72_l3p0_3r_cs0p18` (1876 snapshots, guarded, died at the FMM
  adequacy gate) — a third point linking σ growth to the adequacy deaths.

## Cautions

- Ladder C is **unguarded** (no `SIGMA_CEIL`) while Ladder S is guarded. That
  confound is the reason for the control above; do not report a Ladder-C-only
  finding as a property of static Cs.
- Cs=0.34 died at rev 10 and Cs=0.002 at rev 24.9. Tempting to read as
  "larger Cs destabilises faster", but with n=2, unguarded, and two different Cs
  values spanning 170×, it is a lead for the σ forensics to explain — not a result.
- Reading VTPs is I/O-heavy (each is tens of MB). Subsample steps first to find
  the onset, then go dense only near it; delegate bulk sweeps to a subagent so
  raw output stays out of the main context.
