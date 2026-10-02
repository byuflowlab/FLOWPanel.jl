# P2 rerun matrix — σ-guard rescue of the NT72 GPU class (2026-09-21)

Ryan-directed (2026-09-21) rerun of the failed 032-reopen arm P2
(`p018_csarc_n2_nt72_l3p0_3r_sfs3nb_om15`, job 13774449) after the revised
ignition diagnosis. DRAFT — job ids and the omission threshold are filled in
at submission; see the submission block below.

## Revised P2 diagnosis motivating this matrix (established 2026-09-21)

Tiebreaker harvest on the P2 run dir (correct x-axis cylinder convention,
VTP totals cross-validated against monitor04 exactly):

- Single Γ-ignition: last stable step 1729 (rev 24.0, CT_bernoulli 0.133,
  already elevated); step 1750 CT → −1084.7. Max |Γ| grows 0 → 13.8 (1729,
  rank-1 ~2000× the 99.9th pct) → 1563 (1740) → 6487 (1750).
- Particle count 158k (1700) → 115k (1740) → 7.3k (1750): the wake is
  blasted upstream + radially outward and the 3R truncation cylinder culls
  it. Culling is a SYMPTOM; the `3r` lever is exculpated as cause.
- The earlier "frozen offender at 3.95R outside truncation" was a
  coordinate artifact (sqrt(x²+y²) instead of sqrt(y²+z²)); recomputed the
  offender is inside the domain. No culling bug (trim removes correctly,
  `FLOWPanel_wake.jl` `apply_particle_maintenance!`).
- Ignition seeds: deep wake x≈2.8–3.5R, including a root-radius column
  particle (r_perp/R 0.169) — omission suppresses shedding, not downstream
  accumulation. Seed σ 7.3e-5–1.8e-3 m (σ/R down to 0.0006, ~40× below
  shed σ) and P2 ran SIGMA_FLOOR_FRAC=0 / SIGMA_CEIL=Inf → the small-σ
  Γ-ignition channel had no arrest. **Lever under test here: the g25 σ
  guard.** (SIGMA_CEIL rerun stays deprioritized: max σ was stable −2.9%
  pre-ignition.)
- Secondary observation: ignition at rev 24.3, two revs after
  P018_SETTLE_REVS=22 ended — freestream-withdrawal transient is a
  plausible trigger; noted, not tested in this matrix.

## Rulings applied (Ryan 2026-09-21)

- Innermost-only shed omission for the new-config arms:
  PARTICLE_OMIT_ROOT_R_OVER_R=0.12 masks exactly 1/41 stations per blade
  (verified by local driver probe 2026-09-21; re-verify in the submission
  banner).
- All new-config arms use standard DynamicSFS with exact-rate rlxf
  (SFS_RLXF=0.0025031 at NT72 per r(NT)=1−(1−r36)^(36/NT), r36=0.005) —
  no SFS_THREELEVEL.
- Matrix extended with the minimal-delta control (A1) whose sole purpose
  is rescue attribution: it keeps P2's sfs3nb + om15 and adds ONLY the
  guard.
- Merge-overlap gate arm at MERGE_OVERLAP=4 (Φ_merge=4 > child/shedding
  overlap 2.75, satisfying the split-stability constraint).
- Splitting arms use the 026 §22.3 wave-2 quartet (rerunslate provenance):
  WAKE_SPLIT_VISCOUS=true, FRAC_VISCOUS=0.587, FRAC_COMPRESS=0.73,
  FRAC_ELONGATE=0.3. SIGMA_CEIL stays at the guard's 0.030 m (k=3 cap
  retune still owed — NOT reused here).

## Pins (unchanged; no code changes needed — all knobs present at pin)

- FLOWPanel: tag `campaign/p032-rootomit-20260918` (`af92740`), worktree
  `orc:/home/rander39/campaigns/p032-rootomit-20260918/FLOWPanel.jl`, env
  `.../p032-rootomit-20260918/env` (GPU, h200) — same repo/env as failed
  P2, so A1 is exact-P2-plus-guard.
- FLOWVPM: `8d4a3b4` (new merge law, production default per Ryan
  2026-09-19) — same as P2.

## Arm matrix

Base (all arms): case `p018_csarc_n2_nt72_l3p0` (champion shedding law:
sigma_overlap + SIGMA_CHORD_FRACTION=0.313, OVERLAP=2.75, NWAKEROWS=2,
DAS λ=3.0 arc-placed steady table, MERGE_R_FACTOR=0.0055, RELAX_RLXF
0.16334, CoreSpreading β=1e9), h200 GPU, TRUNCATION_RADIUS_R=3.0,
MAX_PARTICLES=1500000, P018_SETTLE_REVS=22, **g25 guard:
SIGMA_FLOOR_FRAC=0.25, SIGMA_CEIL=0.030** (guard CT offset +0.39% known —
compare slopes vs unguarded history, not levels).

| Arm | Run name | Delta vs base | Job |
|---|---|---|---|
| A1 ctrl | `p018_csarc_n2_nt72_l3p0_3r_sfs3nb_om15_g25` | exact failed P2 + guard: SFS_THREELEVEL=true, OMIT_ROOT=0.15 (3 stations) | 13842791 |
| A2 base | `p018_csarc_n2_nt72_l3p0_3r_srlx_g25_omi1` | SFS_RLXF=0.0025031, innermost-only omission (0.12 → 1/41) | 13842792 |
| A3 | `..._omi1_mo4` | A2 + MERGE_OVERLAP=4 | 13842793 |
| A4 | `..._omi1_split` | A2 + 026 split quartet (absolute merge_r) | 13842794 |
| A5 | `..._omi1_split_mo4` | A2 + split quartet + MERGE_OVERLAP=4 | 13842795 |

Submitted 2026-09-21 (Ryan "Go"); walls 16 h (A1–A3) / 24 h (A4–A5);
job names fp-018gpu-p2rr-a1..a5. Banner verification owed once running:
omission count (3/41 for A1, 1/41 for A2–A5), guard=on with floor
0.25/ceil 0.030, SFS label (threelevel for A1, DynamicSFS rlxf=0.0025031
for A2–A5), MERGE_OVERLAP and WAKE_SPLIT lines on A3–A5, and the GPU-path
gate.

Questions: A1 = does the σ floor alone rescue NT72? A2 = new baseline
(guard + champion SFS + minimal omission). A3 = overlap-gate merging at
NT72. A4 = splitting viability once ignition is suppressed. A5 =
split×merge-gate interaction.

## Acceptance

- Rescue: arm survives past step 1750 (P2 ignition) and completes 2160
  steps (30 revs) with finite CT and no Γ-ignition signature (monitor04
  max|Γ|/σ² stays bounded).
- If A1 survives: σ guard is established as the rescue lever for the NT72
  class → per the 09-19 reopen conditional, the NT144 offer question can
  be re-raised with Ryan (still gated on him).
- If A1 dies but A2 survives: rescue is confounded between SFS and
  omission width — escalate before further arms.
- Splitting arms judged additionally on particle-count trajectory vs
  MAX_PARTICLES=1.5M and merged-σ behavior.

## Submission block (filled at submit)

Template (mirrors the recovered P2 submit line; sacct SubmitLine of
13774449):

```
sbatch --job-name=<name> --partition=eng --qos=eng --constraint=intel \
  --gres=gpu:h200:1 --cpus-per-task=64 --no-requeue --mem=192G \
  --time=<wall> \
  --output=/home/rander39/campaigns/p032-rootomit-20260918/FLOWPanel.jl/logs/slurm/slurm-%x-%j.out \
  --error=/home/rander39/campaigns/p032-rootomit-20260918/FLOWPanel.jl/logs/slurm/slurm-%x-%j.err \
  --export=ALL,<arm env>,P018_REPO_OVERRIDE=/home/rander39/campaigns/p032-rootomit-20260918/FLOWPanel.jl,P018_PROJECT_OVERRIDE=/home/rander39/campaigns/p032-rootomit-20260918/env,P018_RUN_NAME=<run name> \
  examples/run_dji9443_hover_ct_gpu.slurm.sh h200 p018_csarc_n2_nt72_l3p0
```

Wall: 16 h for A1–A3 (P2 completed 2160 steps within 14 h), 24 h for the
split arms (particle growth). Guard env common to all arms:
`SIGMA_FLOOR_FRAC=0.25,SIGMA_CEIL=0.030,TRUNCATION_RADIUS_R=3.0,MAX_PARTICLES=1500000,P018_SETTLE_REVS=22`.
