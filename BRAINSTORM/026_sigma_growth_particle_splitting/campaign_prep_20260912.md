# 026 campaign prep package — commit-7 arms (2026-09-12, DRAFT for Ryan gating)

Prepared per Ryan's 2026-09-12 three-task directive
(`sigma_check_gpusplit_campaignprep_reset_prompt_20260912.md`). Nothing here
is launched; tags/worktrees are staged only after Ryan approves the commit
sequencing in §4.

## 1. Split fractions (fraction space)

§21 already rules the operating points post-§19 fractional gating:
`WAKE_SPLIT_FRAC_COMPRESS = 0.73`, `WAKE_SPLIT_FRAC_ELONGATE = 0.3`,
adaptive elongation at default Φ_t (= shedding `OVERLAP`). The remaining
§19 obligation was the cap-arm adequacy check ("splits fire before
particles sit long at the clamp"):

- Compress trigger fires at attempted σ = (1 + f_comp)·σ₀ = **1.73 σ₀**.
- Cap arms shed σ₀ = 0.313·c(r) (`SIGMA_CHORD_FRACTION=0.313`); DJI 9443
  max chord ≈ 0.025 m ⇒ σ₀ ≤ ~0.0078 m ⇒ trigger σ ≤ **~0.0135 m**, below
  BOTH clamps (cap018 = 0.018 with ~25% margin at the widest chord; cap030 =
  0.030 with ~2.2×). Compress splits fire before either clamp at every
  station, and §19 ruling 2 (attempted-Δσ² accrual) bounds clamp dwell even
  for outliers. *(Chord bound is a literature value — verify at submission
  from the logged `sigma_max / current_adequacy_limit` census, which the §5
  ruling already requires per radix rebuild.)*
- Shrink arms: the old D4 threshold re-derivation is superseded by the §21
  ruling f_elong = 0.3 (elongate fires when attempted σ < 0.7 σ₀; m = 3
  children/fire confirmed in the 2026-09-10 driver smoke and smoke A/C).

## 2. Arm matrix (14 gated arms, defined in `run_p018_screen_hpc.slurm.sh`)

| group | arms | submission knobs (beyond dispatcher defs) |
|---|---|---|
| §8.4 grow-side caps | `scr_p026sp_nt144_cap030`, `_cap018` | `WAKE_SPLIT_FRAC_COMPRESS=0.73 WAKE_SPLIT_FRAC_ELONGATE=0.3`; warm-start from archived `p018_csarc_n5_nt144_l2p4_s2gpu` steps 2200–2248 |
| §9 s020v shrink matrix | `scr_p026s9_{ctrl,exp,ctrllg,explg}_{floor,split,fs}` (12) | floor arms `SIGMA_FLOOR_FRAC=0.25`; split/fs arms `WAKE_SPLIT_FRAC_ELONGATE=0.3`; fs = both; warm-start at bracket steps (ctrl 224/226, exp 209/211, ctrl-lg 229/231, exp-lg 284/286) |

Per-arm overlap/merge proposals:

- `WAKE_SPLIT_ELONGATE_OVERLAP`: leave at default Φ_t = shedding `OVERLAP`
  (2.4 for s9 arms, 2.75 for cap arms) per §21. Both satisfy the merge
  stability constraint Φ_merge = 3.5 > Φ_t (children placed at overlap Φ_t
  do not meet the merge gate at birth — confirmed empirically in Task 1:
  sibling re-merges were 0–3% of merge events in smoke B).
- **First merge discriminator (per §21): `MERGE_OVERLAP=3.5` vs
  production-absolute, exp bracket first** — i.e. add
  `scr_p026s9_exp_split` twin runs differing only in `MERGE_OVERLAP=3.5`
  vs the dispatcher's `MERGE_R_FACTOR=0.00275` absolute radius.
  Regime note from Task 1 (see §5 below): at the s9 arms' shed σ ≈ 0.0381R
  the crossover σ* = 3.5·r_abs sits ≈ 2× above shed σ, so the overlap gate
  merges young particles LESS and aged particles MORE than production —
  the theory-doc §4 σ-pump A/B is live in this regime. (The 09-12 smoke B
  ran a coarse-σ config, σ_shed = 0.26R, where the overlap gate is wider
  than absolute for every particle — its 2.7× retention drop does NOT
  transfer to the campaign arms.)
- Backend: GPU for the particle-heavy arms (Ryan 2026-09-05 preference),
  contingent on §4 commits + CUDA smoke. Keep discriminator pairs
  backend-matched (both GPU or both CPU).

## 3. Reproducibility staging (per global campaign policy)

- Annotated tags `campaign/p026-splitmerge-20260912` in FLOWPanel, FLOWVPM,
  FastMultipole; worktrees created from the tags; campaign Julia env
  Manifest dev-paths pointed at the worktrees; pins (tag + SHA) recorded in
  a provenance file in the campaign data root.
- Pin state (updated 2026-09-14, all committed AND pushed):
  - FLOWVPM `flowpanel`: `bf88806` (device accumulators `edc9d95` +
    device-safe lifecycle hooks).
  - FLOWPanel `fastmultipole`: `035f50b` (seam sync `a804a95` +
    verification/mo35 case defs + gate-safe banner).
  - FastMultipole `flowpanel-20260817`: `ac7230a6`.
- Tagging proceeds once the v3 CUDA smokes (13689273/74/75) pass.

## 4. Commit sequencing proposal (needs Ryan approval)

1. **FLOWVPM commit** (after reverting the Task 1 merge-logger diagnostic):
   parameterized `ResolutionSplitState` (backing arrays match the pfield's
   storage type), broadcast accumulation twins in the three integrators +
   three CoreSpreading variants, t12 CPU-vs-broadcast parity testset.
   Verified: `runtests_resolution_split.jl` 1055/1055.
2. **FLOWPanel commit**: `_gpu_sync_rsplit!` seam sync (widened H2D against
   stale device tails), device-side `enable_resolution_split!` in
   `_apply_particle_maintenance_device!`, updated design comments.
   Verified: `runtests_unit_wake.jl` clean, `runtests_unit_replay.jl`
   148/148 (pre-existing warmstart first-testset failure untouched).
3. **HPC CUDA smoke** (verification, not campaign): clone
   `scr_p026ph1b_expgpu_smoke` + `WAKE_SPLIT_STRETCH=true
   WAKE_SPLIT_FRAC_COMPRESS=0.73 WAKE_SPLIT_FRAC_ELONGATE=0.3`, plus a
   `MERGE_OVERLAP=3.5` companion; accept on nonzero split counters, sane σ
   census, split-event rates comparable to a backend-matched CPU run, no
   per-step cost blowup (rs sync adds 7/45 more rows to the ~2×67 MB
   maintenance traffic).
4. Then tag + worktrees + provenance (§3).

## 5. Task 1 verdict (σ-distribution check, smoke B) — FINAL

Evidence: smoke B surviving VTK (steps 322–467; 0–321 were overwritten by
the local smoke C rerun — B's 322–467 backed up to session scratchpad) +
an instrumented probe (merge-event logger patched into merge_particles!,
uncommitted and since reverted; B config rerun locally to step 236, np at
233 = 8934 vs B's 8876, 0.7% ≈ RNG-level — faithful reproduction; 21.9k
logged merge events).

**Verdict: B's 2.7× lower retention is split→merge churn of a specific,
bounded kind — the overlap gate consuming freshly split children across
overlapping filament FAMILIES — not healthy aged-wake σ-thinning, and not
either pathological channel:**

- Sibling re-merge (children of one split event): 0–3% of merges. The
  Φ_merge = 3.5 > Φ_t placement margin works as designed.
- At-release shed merging: 0.3–0.8% of merge members. Dead channel.
- Dominant flux (stationary across passes 30–233, ~96–103 merges/step):
  **67% split-child × split-child pairs, 72–78% of members split children**,
  pair σ_min ≈ 0.47–0.88 σ_shed (median 0.65–0.68), **98% of members fresh
  (σ ≈ σ₀)** — children merge soon after birth, with other families'
  children in braided/dense regions.
- Aged-wake σ-pump (merge products re-merging): ~5–9% of flux; visible as
  a slow σ-tail extension in the snapshots (σ q90 0.0434 → 0.0481 m over
  steps 322→457) but a minor removal channel.

**Regime caveat (critical for §2):** this smoke sheds at σ = 0.26R
(overlap_pps), where the overlap gate radius (σ/3.5 ≈ 0.0089 m at shed) is
~3.7× WIDER than the production absolute radius everywhere — the retention
drop is baked into the gate width. The campaign s9 arms shed at
σ ≈ 0.0381R where the crossover σ* = 3.5·r_abs sits ~2× ABOVE shed σ: the
overlap gate is then TIGHTER than production for young particles and wider
only for aged ones. Additionally the campaign child overlap Φ_t = 2.4/2.75
(vs the smoke's 3.0) leaves 1.46×/1.27× relative margin to the 3.5 gate
(vs 1.17×), so cross-family child churn should be substantially weaker in
the campaign arms. Smoke B therefore neither validates nor indicts
MERGE_OVERLAP=3.5 for the campaign — it confirms the mechanism inventory
and that the A/B (per §21, exp bracket first) is the right discriminator.

## Launch rulings (Ryan, 2026-09-14)

- **3-arm de-risk first**: `scr_p026sp_nt144_cap030` + `scr_p026s9_exp_split`
  + its merge-A/B twin `scr_p026s9_exp_split_mo35` (MERGE_OVERLAP=3.5,
  new case def). Remaining 11 arms follow if the trio passes.
- **Cold-start** (GPU is much faster): no RESTART_* env; exp arms run past
  the ignition bracket (~step 210) from step 0; cap030 runs through the
  original cliff region (~step 2250, NT144).
- **`f_visc` enabled** ("I don't think it will trigger"): value not ruled
  in §21 — adopted `f_visc = 4^(1/3) − 1 ≈ 0.587`, the count-matched tetra4
  analog of the ruled `f_comp = √3 − 1` (Ryan picked 0.587 from the
  AskUserQuestion options, 2026-09-14). Pairs with
  `WAKE_SPLIT_VISCOUS=true` (driver fail-fast requires the pair).
- **Guard floor at 0.1 for the de-risk trio** (`SIGMA_FLOOR_FRAC=0.1`
  instead of 0.25): if the trio passes, the rest of the campaign keeps
  0.1; if it fails, bump back to 0.25.
