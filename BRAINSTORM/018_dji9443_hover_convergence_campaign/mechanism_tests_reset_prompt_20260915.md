# RESET PROMPT — 018 Γ-distribution mechanism tests (2026-09-15, post-Γ-analysis)

You are a clean-context agent working BRAINSTORM item 018 (DJI-9443 hover
convergence campaign). Skim the header + standing rulings of
`../018_dji9443_hover_convergence_campaign.md`, then `ops_reference.md` and
`decision_rules.md`. Never tune to CT_exp; ask Ryan before writing anything to
the lab notebook or `ledger.md` (both stale since 2026-09-08; notebook owes
entries for Cd transient, rlxf derivation, Ladder C forensics, expguard arc,
the 2026-09-15 audit + rlxf-scaled pair, and the 2026-09-15 Γ-distribution
result).

## Where the campaign stands (do NOT re-derive)

- NT anti-converges: unguarded reference climb +1.51%/doubling; guarded
  same-arch g25 pair with exact-rate-scaled relaxation still climbs
  **+1.27%/doubling** (`rlxfscaled_provenance_20260915.md`):
  - job 13704962 `p018_csarc_l3p0_3r_g25` (NT36): CT̄(21–30) = 0.070519
  - job 13704963 `p018_csarc_n2_nt72_l3p0_3r_srlx_g25` (NT72): CT̄ = 0.071417
- **2026-09-15 Γ-distribution result** (`gamma_distribution_status_20260915.md`
  — read it, it has all tables): the Γ distribution difference vs NT is a
  **far-field, wake-evolution, substantially TEMPORAL effect**:
  - Near-disk (<1.5R) Σ|Γ| is NT-invariant (≲3–4%) in both stacks at all revs
    after ~12. Shedding/conversion class is ruled out for this signal.
  - Far wake (≥1.5R): NT36's Σ|Γ| saturates (~3.9–4.0 by rev 29) while NT72's
    plateaus revs 18–22 then RE-ACCELERATES (3.72 → 4.33 over revs 24–29,
    still growing at run end, crossing NT36 near rev 28). NT72's wake is not
    statistically steady anywhere in the rev-21–30 scoring window.
  - CT window split agrees: climb widens +1.17% (revs 21–25) → +1.41% (26–30);
    NT36 CT drifts down late while NT72 holds/rises.
  - Count-normalized: NT72 has ~6–9% more particles; far-field median |Γ|
    drops ~30–50% per doubling (more, weaker particles) — but far-field Σ|Γ|
    itself differs at 10%+ so it is NOT purely benign discretization.
  - Cross-check `_3r_sv_s1p5` NT36/72/144 at rev 30: near-disk invariant
    within 2.7% over 4× NT; far field differs 12.5% with the OPPOSITE sign
    (NT36 > NT72 ≈ NT144) — consistent with rungs sitting at different points
    on NT-dependent equilibration curves.
  - Single matched-rev snapshots are USELESS at the ±10% level (structure
    aliasing ≈ signal); always phase-average ≥12 snapshots over a full rev.
    Analysis scripts live on ORC in `~/p018_gamma_dist_20260915/`
    (`p018_gamma_dist_avg.py`, `p018_gamma_ts.py`, logs `gdist_avg.log`,
    `gdist_ts.log`).

## Mechanism landscape (from run metadata + FLOWVPM source survey, 2026-09-15)

Run metadata: `<rundir>/<run>.metadata.toml` (pfield_optargs section) and
`<run>_case_metadata.toml`. Facts established:

1. **Merging is ON, per-step, absolute-radius**: `[[wake.particle_maintenance
   .functional_policies]] type="MergeParticles", every=1, r=6.545e-4,
   r_hash=2.38e-3, sigma_relative=false, max_sigma_ratio=2.0`. Per-step
   cadence ⇒ NT72 gets 2× merge passes per unit time. Merging sums Γ vectors
   ⇒ cancellation destroys Σ|Γ|, preferentially in the crowded old wake.
   FLOWVPM exposes it as `run_vpm!` kwargs `merge_every` (0 = off) and
   `merge_kwargs` (`src/FLOWVPM_merging.jl`, `merge_particles!` sig at :440);
   the p018 driver exposes `merge_particles = true` in case metadata — find
   how the driver maps config → `merge_every`/`merge_kwargs` before assuming.
2. **Merging (+ the geometric trim) is the ONLY wake eater**: FLOWVPM has NO
   strength-threshold pruning anywhere (verified by source survey; removal =
   `remove_particle` primitive called by merging, static-particle
   bookkeeping, and the GlobalCylinder trim: radius 0.357 m, extrude
   [0.476,0,0], origin [-0.0595,0,0] — NT-independent, gives the observed
   x_max = 3.47R). So the documented removal-side asymmetry ("old wake eaten
   ~2× faster at NT36") is almost certainly carried by merging. **Merging is
   the front-runner.**
3. **Pedrizzetti relaxation is already exact-rate-matched in the g25 pair**:
   NT36 rlxf=0.3, NT72 rlxf=0.16334, (1−0.16334)² = 0.700 = 1−0.3. Excluded
   beyond the ~16% already measured.
4. **SFS = DynamicSFS pseudo-3-level with clipping_backscatter**
   (`src/FLOWVPM_subfilterscale.jl`; DynamicSFS at :277, clipping at :435).
   Dynamic-coefficient dynamics (test filter, per-step clipping events,
   Lagrangian-avg rlxf) are not guaranteed NT-invariant. Set at
   ParticleField construction (driver-side), not a run_vpm! kwarg.
5. **Viscous = CoreSpreading, nu=1.4334e-5, beta=1e9** ⇒ core reset never
   triggers; σ growth is a pure per-unit-time ODE. Weak candidate; cheap
   toggle (`Inviscid`).
6. Shed quantization (`max(1,ceil)` stations, +10% counts/doubling) lives on
   the FLOWPanel/driver side, not FLOWVPM. Not in this wave (near-disk Σ|Γ|
   is matched, so it is not the Σ|Γ| carrier), but it drives the benign
   "more, weaker" shape change.

## Your task (Ryan, 2026-09-15): launch mechanism-isolation wave 1 on ORC

All runs: g25 stack (same arch/gh200 pattern as jobs 13704962/13704963 —
copy their sbatch/driver configs as the template), guarded, rev length as
specified. Ryan approved this matrix (chat, 2026-09-15):

- **T1 — extend both g25 rungs to rev 60, RESTART from rev 30** (Ryan chose
  restart over fresh). Both full per-step VTP sets are safe: harvested to
  `/nobackup/archive/usr/rander39/FLOWPanel_runs/
  p018_csarc_l3p0_3r_g25_wake1_particles_full_20260915.tar` (33 GB, 1080
  members, verified) and `..._n2_nt72_l3p0_3r_srlx_g25_..._full_20260915.tar`
  (70 GB, 2160 members); originals still live in the run dirs.
  - **Restart-fidelity gate first**: metadata flags integration/kernel/
    formulation `restart_reconstruct_required = true`. Locate the p018
    driver's restart path (how it reconstructs the pfield from saved state —
    likely reads a particle VTP + case config). Gate stage: restart NT36
    from step ~1044 and rerun to 1079; compare CT per step and a wake
    snapshot vs the original (tolerance: CT trace visually indistinguishable,
    per-step ΔCT ≲ 1e-4 relative after the first couple steps). Only if the
    gate passes, launch the two rev-30→60 restarts. If restart is unreliable
    or unsupported, STOP and report to Ryan (he accepted restart only as the
    cheaper option; fresh 60-rev runs are the fallback but need his sign-off).
  - Watch maxparticles headroom: NT72 counts were still growing at rev 30
    (~340k avg). Check the driver's maxparticles setting before submitting.
- **T2 — NT72 merge-cadence rate-match**: clone of
  `p018_csarc_n2_nt72_l3p0_3r_srlx_g25` with merging every 2 steps
  (`merge_every = 2` equivalent in the driver config), 30 revs. Suggested
  name: `p018_csarc_n2_nt72_l3p0_3r_srlx_mrg2_g25`. This matches merge passes
  per unit time to NT36. No NT36 counterpart needed (NT36 side unchanged).
- **T4 — SFS frozen (ConstantSFS), BOTH rungs**, 30 revs. Replace DynamicSFS
  with `ConstantSFS` at a fixed Cs. Choice of Cs: extract the run-averaged
  dynamic coefficient from an existing g25 run if the driver/monitors expose
  it; otherwise propose a standard value to Ryan before submitting. Do NOT
  use NoSFS (052 lesson: stretching runaway/ignition risk in unguarded
  no-SFS wakes).
- **T5 — Inviscid, BOTH rungs**, 30 revs. Swap CoreSpreading → Inviscid.
  Cheap stage; keep SFS as-is.
- T3 (merging fully off) was explicitly deferred — do not run it this wave.

Staging: per the cluster-jobs rule, combine ready GPU stages into one sbatch
where envs/resources match (multi-stage `cuda_048_run.sh` pattern); the
restart-fidelity gate must complete before the T1 long restarts. T2/T4/T5 are
independent of T1 and can go in parallel stages/jobs.

**Campaign rules apply** (these are cited runs): commit any driver changes,
pin with annotated tags (`campaign/p018-mech-tests-20260915` convention),
run from worktrees created from the tags, Manifest dev-paths at the
worktrees, write a provenance file
(`mechanism_tests_provenance_20260915.md`) in this directory with tags(+SHAs),
job IDs, and per-run configs BEFORE submitting. Run outputs to the standard
data root, not the worktree.

## Scoring / analysis (same for every test)

For each completed run, compute — reusing/adapting the ORC scripts in
`~/p018_gamma_dist_20260915/`:

1. Per-rev far-wake time series (near <1.5R vs far ≥1.5R Σ|Γ|, 4-phase avg):
   `p018_gamma_ts.py`. The diagnostic signature to kill: NT72's late
   far-wake re-acceleration and the NT36/NT72 far-wake divergence.
2. Phase-averaged per-bin ratios at rev 29→30 (12 phases):
   `p018_gamma_dist_avg.py`.
3. CT̄ with split windows: `python3 scripts/p018_analyze.py m1 --revs A B
   <run>` from `~/projects/FLOWPanel.jl` (one run per invocation; use 21 25,
   26 30, 21 30; for T1 also 31 60 windows, e.g. 41 50 / 51 60).

Verdicts sought:
- T1: does NT72's far wake saturate by rev ~45–60, and does the CT climb
  shrink as the window moves late? If yes → residual carrier is (partly) a
  transient-length artifact; report the asymptotic climb.
- T2: does rate-matched merging collapse the far-wake Σ|Γ| divergence and/or
  the climb? Quantify: corrected climb and far-wake curve vs the g25 pair.
- T4/T5: same comparison; a mechanism is implicated if the divergence
  signature shrinks by >~50% with the toggle.

## Data and tooling gotchas (ORC, `ssh orc`, user rander39)

- Run dirs `~/projects/FLOWPanel.jl/data/<run>/`; particles in
  `<run>_wake1_particles/<run>_wake1_particles.<step>.vtp` (unpadded steps,
  step 0 exists; NT36 30 revs = steps 0–1079, NT72 = 0–2159).
- Reader `~/vtp_C.py` (`read_vtp(path, {"Points","gamma","sigma"})`,
  `run_files(rundir)` — takes the RUN dir, not the particles subdir;
  `step_of`). Raw-appended VTK XML; meshio cannot read them.
- Remote python/slurm need a login shell: `ssh orc 'bash -lc "…"'`. The
  banner's ANSI codes glue to the first stdout line — start with a bare
  `echo`, never anchor greps to line starts on line 1.
- Print summaries only; never cat VTP/CSV bytes into context. Axial axis is
  x (disk at x≈0, wake → +x, R = 0.11995 m); do NOT pick the axis by max
  extent.
- VTK sweeper culls non-live runs to newest-36 VTPs (`vtk_protect_list.txt`
  is Ryan's — agents read, never write). New test runs you launch are live
  while running; if any full VTP set matters afterward, harvest to
  `/nobackup/archive/usr/rander39/FLOWPanel_runs/` like the g25 tars above.
- The g25 sbatch scripts / driver configs for jobs 13704962/13704963 are the
  template for every clone — locate them (likely under
  `~/projects/FLOWPanel.jl` scripts or the job's submit dir via
  `sacct -j 13704962 --format=WorkDir%200`) and diff-minimally.

## Deliverable

Provenance file before submission; after runs complete, a dated status file
in this directory with: the three metrics per test run, side-by-side with the
g25 baselines; a verdict per mechanism (merging cadence / SFS dynamics /
viscous / transient-length) with the far-wake divergence signature as the
primary discriminator; and a recommendation for wave 2 (T3 merging-off is the
pre-approved follow-up if T2 is positive). Offer (don't write) notebook
entries — the backlog list is at the top of this prompt.
