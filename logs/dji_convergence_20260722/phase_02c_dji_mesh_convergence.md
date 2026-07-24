# Phase 2c Log — DJI 9443 Chordwise Mesh Convergence

Plan: [Phase 2c — DJI 9443 Chordwise Mesh Convergence](../../plans/dji_convergence_20260722/phase_02c_dji_mesh_convergence.md)

Dashboard: [DJI convergence progress](../../plans/dji_convergence_20260722/dji_convergence_20260722.md)

## Current snapshot

Status: **Phase 2c — COMPLETE; decision recorded, awaiting Ryan's Phase 3 approval.**

Goal: verify on the actual DJI 9443 rotor mesh that the Dirichlet–Neumann
bound-circulation gap converges under chordwise refinement, as the Phase 2b oracle
predicted.

**Decision: FINER MESH DESIRABLE (attribution direction supported, not converged; one
outlier rung flagged).** All 6 solves clean (monitor/direct error 0.0, blade symmetry
≤1.5e-13). The three trustworthy signals support the Phase 2b attribution:

- **Neumann reference flat/converged** across all rungs (∫Γ Neu per-rung Δ ≤ 0.30%) —
  the oracle's "open Neumann chordwise-converged by ~60" behavior.
- **Outboard (|r/R|≥0.7) gap shrinks monotonically** 3.325% → 2.954% → 2.686%, Dirichlet
  climbing toward Neumann (the oracle mechanism) — but **still 2.69% at n=121, not ≤1%.**
- **Full ∫Γ gap is NON-monotonic** (3.764% → 1.011% → 2.975%) solely because the
  **n_airfoil=97 capped/Dirichlet mesh is an inboard outlier**: its Dirichlet *inboard*
  circulation spikes to the Neumann level (3.774→**3.915**→3.779) then recedes, while
  Neumann inboard stays flat (3.926→3.911→3.898) and Dirichlet outboard stays monotone.
  Same single-rung pattern as the documented rotor-hover 56_57 outlier.

Sanity check passed: the n=81 gap (3.76% full / 3.33% outboard) is below Phase 1's
4.85% at n=57, consistent with the attribution.

**Spanwise addendum (2026-07-24, Ryan-supplied 60_97 meshes):** the n=97 outlier is a
30-spanwise **mesh defect, not a resolution or attribution problem**. Solved 60_97 capped
(24160 panels) + uncapped (22656 panels), compared to the 30_97 pair on 35 matched bins.
**Spanwise is converged**: the outlier-free indicators — the Neumann solve (∫Γ change
30→60 = +0.065%, outboard −0.285%) and the outboard integral (Dir −0.20%) — are
spanwise-flat (≤0.3%). The large Dirichlet full-∫Γ move (−2.15%) is the 30_97 capped
**inboard** defect being corrected: the 60_97 Dirichlet ∫Γ = 0.017825 rejoins the
chordwise family (81c 0.017786 / 121c 0.017824), and the inboard Dir–Neu gap opens from
the anomalous −0.20% back to the family +3.57%. With 30_97c replaced by 60_97c the
Dir–Neu ∫Γ gap is **monotonic**: 3.76% (81) → 3.30% (97) → 2.98% (121). The Phase 2c
non-monotonicity was purely the 30_97 mesh.

**Recommendation to Ryan (updated):** the 30-spanwise grid is adequate — spanwise is
converged. **Regenerate/replace the 30_97 capped mesh (or use 60_97).** The gap is still
~3% at n=121, so the **finer-chordwise** recommendation stands (e.g. 30_145, 30_161).
Attribution NOT challenged. **Do NOT proceed to Phase 3/5 until Ryan approves.** Artifacts:
`spanwise_report.md`, `spanwise_metrics.csv`; driver mode `PHASE2C_MODE=spanwise`.

**Extended chordwise addendum (2026-07-24, Ryan-supplied 45-span 145/185/201/249):** the
DJI Dir–Neu gap **converges monotonically toward the Neumann reference under chordwise
refinement — the Phase 2b oracle mechanism, now confirmed on the real DJI mesh — but
slowly.** Best-point-per-n_airfoil ladder (n=97 from the corrected 60-span solve):

| n_airfoil | ∫Γ gap % | outboard gap % |
|---:|---:|---:|
| 81  | 3.76 | 3.33 |
| 97  | 3.30 | 2.87 |
| 121 | 2.98 | 2.69 |
| 145 | 2.73 | 2.55 |
| 185/201/249 | pending HPC | pending HPC |

Neumann is flat (Δ ≤ 0.32%/rung), Dirichlet now nearly flat (Δ 0.11% at 145), gap closing
~0.24 pts/rung. Not yet ≤1% at n=145; extrapolation suggests many more chordwise panels
to reach 1% (decay is slow). Only n=145 solved locally (27608/25344 panels, ~6 GB dense
Backslash, 413/112 s, error 0.0). **185/201/249 deferred to HPC** (dense Backslash needs
10/12/18 GB > the 17 GB laptop; 249 exceeds RAM entirely — Ryan chose "145 local, rest
HPC"). Artifacts: `extended_report.md`, `extended_fixed_bin.csv`; driver mode
`PHASE2C_MODE=extend`.

**Mesh-format note:** Ryan's first 45_145/45_185 exports were MSH 2.2 (GeoIO needs ≥4.1);
he re-exported them as 4.1 and added 45_201/45_249 (all 4.1). Driver uses the original
filenames. A transient meshio 4.1 conversion was used and then removed.

### HPC submission for 185/201/249 (Ryan approved 2026-07-24)

**Submitted: Slurm job `12890990`** (`fp-dji9443-mesh-ext`) on BYU `orc`
(ssh alias, `ssh.rc.byu.edu`), repo `/home/rander39/projects/FLOWPanel.jl` (branch
`fastmultipole`). 1 node, 64 threads, 32 GB, 12:00:00.

Deployed via `scp` (no git ops on the cluster — its worktree had unrelated dirty files,
left untouched): the driver `examples/dji9443_mesh_convergence.jl`, the TE detector
`examples/dji9443_trailing_edge.jl`, the launcher
`examples/run_dji9443_mesh_convergence_hpc.slurm.sh`, and the six meshes
`dji9443_20260723_45_{185,201,249}_{capped,uncapped}.msh` (all confirmed MSH 4.1 on the
cluster). The launcher runs each case in a **fresh Julia process**
(`PHASE2C_MODE=solve PHASE2C_CASE=<tag>`, tags `dji45_{185,201,249}{c,u}`) so each large
dense solve gets clean memory (249 capped peak ~18 GB). New driver mode `hpc_ext` also
solves exactly this deferred set in one process if preferred.

**Results retrieval (when the job finishes):**

1. `scp orc:/home/rander39/projects/FLOWPanel.jl/data/dji_convergence_20260722/phase_02c_dji_mesh_convergence/raw/dji45_{185,201,249}{c,u}_circulation.csv`
   into the local `.../raw/`.
2. `PHASE2C_MODE=extend julia --project examples/dji9443_mesh_convergence.jl` locally —
   the pending 185/201/249 rows fill in automatically (driver checkpoint-skips present
   raws; the local-RAM guard only defers solves, not analysis).
3. Re-read the fine-tail trend and decide converged vs finer-still; update this snapshot
   and the dashboard.

Handoff facts carried from Phase 2b (approved 2026-07-24):

- **Attribution**: the Dir–Neu gap is a **chordwise (n_airfoil) section-resolution
  discretization difference** between the Morino source+doublet (Dirichlet) and
  vortex-ring (Neumann) formulations. Open Neumann is chordwise-converged by ~60
  panels/section; capped Dirichlet climbs *up* to meet it.
- **Oracle convergence**: clean solves agree to **0.11% at 120 chordwise
  panels/section** (Dir Δ/rung 0.32→0.19→0.13%, Neu ~flat at ≈0.01%). Below the 1% and
  the stricter 0.25% criteria.
- **Excluded**: solve-formulation/Kutta-trace (0.74%), extraction/mix (loop-circ = Γ_TE
  to 3e-5), thinness (thinner → smaller gap). The earlier "~1% caps term" was a
  drop-one-panel artifact; genuine cap effect ≲0.1–0.2%.
- **Phase 1 reference gaps** (uncapped/Neumann above capped/Dirichlet): +6.519% at
  n_airfoil=41, +4.853% at n_airfoil=57 (40/56 spanwise family — non-matched).

Driver: `examples/dji9443_mesh_convergence.jl` (modes `setup`/`solve`/`all`/`analyze`),
adapted from Phase 1 `dji9443_circulation_audit.jl`, seeds from the TE detector, bins
trimmed to common station support.

## Working records

| Date | Command or file | Purpose/result |
|---|---|---|
| 2026-07-24 | `examples/dji9443_mesh_convergence.jl` | New Phase 2c driver; all 6 meshes pass `setup` (panels 10528/9280/12640/11136/15800/13920; 29 shedding sections/blade; \|r/R\| support [0.1159, 0.9817] covers the full 35-bin 0.125:0.025:0.975 grid, no trimming). |
| 2026-07-24 | `PHASE2C_MODE=all` | 6 FMM solves coarse→fine (67/12/87/15/134/25 s); monitor/direct error 0.0 all cases. |
| 2026-07-24 | `PHASE2C_MODE=analyze` | Convergence tables + `phase_02c_report.md`. Decision: FINER MESH DESIRABLE; n_airfoil=97 capped is an inboard outlier. |
| 2026-07-24 | `PHASE2C_MODE=spanwise` | Ryan-supplied 60_97 meshes. Solved 60_97c (24160 panels, 306 s) + 60_97u (22656 panels, 80 s), error 0.0. `spanwise_report.md`: spanwise CONVERGED (Neu/outboard ≤0.3%); 30_97 capped confirmed a mesh defect; corrected ladder monotonic 3.76→3.30→2.98%. |
| 2026-07-24 | `PHASE2C_MODE=extend` | Ryan-supplied 45-span 145/185/201/249 (first 145/185 were MSH 2.2 → he re-exported 4.1). Added `CHORDWISE_EXT`/`EXTENDED_LADDER` + `extend` mode with local-RAM guard. Solved 145 locally (27608/25344 panels, 413/112 s, error 0.0); 185/201/249 deferred to HPC (10/12/18 GB dense > 17 GB). `extended_report.md`: ∫Γ gap 3.76→3.30→2.98→2.73%, monotonic; CONVERGING (provisional). |

## Dated entries

### 2026-07-24 — Phase 2c staged and setup validated

- Ryan approved Phase 2b (attribution: chordwise section-resolution difference; oracle
  converged to 0.11% at 120 panels/section). Phase 2c inserted before Phase 3.
- Built `examples/dji9443_mesh_convergence.jl` from the Phase 1 driver pieces. Seeds come
  from `find_dji9443_trailing_edge_indices` per mesh (not hardcoded); two-stage body
  build with `ensure_winding=true`; hardcoded Phase-1 count checks dropped in favor of
  recording observed counts.
- `setup` passed on all six meshes; observed panel counts match the plan's inventory
  exactly. Common station support covers the full fixed-bin grid, so no bin trimming was
  needed (35 matched bins).
### 2026-07-24 — Solves + analysis complete; decision recorded

- `PHASE2C_MODE=all` solved all six cases coarse→fine (10528/9280/12640/11136/15800/13920
  panels; wall 67.4/12.1/86.5/15.3/133.7/25.1 s). Monitor-vs-direct max relative error
  0.0 for every case; blade-symmetry relative error ≤1.5e-13.
- `PHASE2C_MODE=analyze` wrote `case_metrics.csv`, `fixed_bin_circulation.csv`, and
  `phase_02c_report.md` on 35 matched bins (no trimming). Convergence table (full ∫Γ,
  gap % = 100·(Neu−Dir)/|Dir|):

  | n_airfoil | Dir ∫Γ | Neu ∫Γ | gap % |
  |---:|---:|---:|---:|
  | 81  | 0.017786 | 0.018455 | +3.764 |
  | 97  | 0.018216 | 0.018400 | +1.011 |
  | 121 | 0.017824 | 0.018354 | +2.975 |

  Outboard gap (|r/R|≥0.7) is monotone: +3.325% → +2.954% → +2.686%.
- Root cause of the full-∫Γ non-monotonicity localized by splitting inboard/outboard
  (`fixed_bin_circulation.csv`): the n_airfoil=97 capped Dirichlet *inboard* sum spikes
  (3.774→3.915→3.779) while Neumann inboard is flat (3.926→3.911→3.898) and Dirichlet
  outboard is monotone (2.285→2.289→2.293). The 97-capped mesh is a single-rung inboard
  outlier — not a failure of the attribution. Mirrors the rotor-hover 56_57 outlier.
- Decision: **FINER MESH DESIRABLE** (see snapshot). Reported to Ryan; Phase 3 stays
  blocked pending approval.

### 2026-07-24 — Extended chordwise ladder (45-span 145/185/201/249)

- Ryan supplied 45-span meshes at n_airfoil 145/185/201/249 (capped+uncapped). First
  145/185 exports were MSH format 2.2 (GeoIO requires ≥4.1) → I converted to 4.1 sidecars
  via meshio to unblock, Ryan then re-exported the originals as 4.1 and added 201/249; I
  reverted the driver to the original filenames and removed the sidecars.
- Added `CHORDWISE_EXT` (the 8 new cases), `EXTENDED_LADDER` (best point per n_airfoil,
  n=97 = corrected 60-span), and `PHASE2C_MODE=extend`. Because `Backslash` is a dense
  solve and this Mac has 17 GB RAM, the mode solves only n≤145 locally (`LOCAL_EXT_TAGS`,
  overridable via `PHASE2C_LOCAL_EXT`) and analyzes whatever raws exist, marking the rest
  pending. Panel→memory: 145=27608 (6.1 GB), 185=35288 (10 GB), 201=38356 (11.8 GB),
  249=47560 (18.1 GB > RAM). Ryan chose "145 local, rest HPC".
- Solved dji45_145c/u locally (413/112 s, monitor/direct error 0.0). Extended ∫Γ gap
  ladder: 3.76% (81) → 3.30% (97) → 2.98% (121) → 2.73% (145), monotonic; outboard
  3.33→2.87→2.69→2.55%. Neumann flat (Δ≤0.32%/rung), Dirichlet near-flat (Δ0.11% at 145).
  **CONVERGING toward the Neumann reference (oracle mechanism confirmed on DJI), but
  slowly — not yet ≤1% at n=145.** Provisional pending the HPC tail.
- 185/201/249 prepared for HPC (six per-case Backslash solves), NOT submitted — study
  rule requires Ryan's approval before deployment. Handoff steps in the snapshot above.

### 2026-07-24 — Ryan approved HPC; submitted job 12890990

- Ryan approved submission. Deployed to BYU `orc` (`/home/rander39/projects/FLOWPanel.jl`,
  branch `fastmultipole`) via `scp` — driver, TE detector, launcher, and the six 45-span
  185/201/249 meshes (verified MSH 4.1 on the cluster). No git operations on the cluster;
  its pre-existing dirty worktree files were left untouched.
- Added launcher `examples/run_dji9443_mesh_convergence_hpc.slurm.sh` (64 threads, 32 GB,
  12 h; one fresh Julia process per case via `PHASE2C_MODE=solve`) and driver mode
  `hpc_ext`. Submitted with `sbatch` → **job 12890990**, PENDING (Priority).
- Retrieval + local `extend` re-analysis procedure recorded in the snapshot's HPC section.
