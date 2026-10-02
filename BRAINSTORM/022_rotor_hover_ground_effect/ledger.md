# 022 ledger — single running results table

One row per run. CT = M1 cycle-mean ± cycle-std (thrust convention, negated
axial Bernoulli channel) over the quoted window. Health = banner verified +
GS converged (gs_nonconverged=0) + tangency bounded (n/a for OGE).

| tag | mesh | ground | job | submitted | landed | window (revs) | CT̄ ± std | health | notes |
|---|---|---|---|---|---|---|---|---|---|
| smoke_ige | 40_40 | h/R 1, 4R/0.15R disc, trunc 3R, none | local | 2026-08-18 | 2026-08-18 | n/a (24 steps) | n/a | PASS (GS 6–7/50 iters, 0 non-conv) | Phase 0 |
| smoke_oge | 40_40 | off | local | 2026-08-18 | 2026-08-18 | n/a | n/a | PASS | Phase 0 |
| p022_oge_fine | 45_185_ct4 | off | 13207679 | 2026-08-18 | 2026-08-20 (all 1007 steps) | final 10 revs, in hover | **0.07480 ± 0.00152** | all_finite, window_in_hover; strict per-rev criterion false (cycle-mean is headline per decision rules) | Phase 1 — CT_OGE fine anchor |
| p022_oge_coarse | 56_57 | off | 13207680 | 2026-08-18 | — | — | — | banner VERIFIED 2026-08-18 | Phase 1, 24 h |
| p022_ige_fine | 45_185_ct4 | h/R 1, 4R/0.15R disc, trunc 3R, none | 13207681 | 2026-08-18 | 2026-08-21 (all 1007 steps) | final 10 revs, in hover | **0.07934 ± 0.00234** | gs_nonconverged=0, tangency RMS max 0.577, below-ground max 4679 (plateau); survived rev 17.6 ⇒ blow-up NOT mesh-independent | Phase 2 — CT_IGE fine anchor |
| p022_ige_coarse_damp | 56_57 | h/R 1, 4R/0.15R disc, trunc 3R, none + damp band 0.1R | 13246557 | 2026-08-20 | — | — | — | banner VERIFIED 2026-08-20 (policy none, damp 0.1R, depth 4R) | Phase 3 A/B "with" arm, 48 h |
| p022_ige_coarse_trunc | 56_57 | h/R 1, 4R/0.15R disc, floor@ground (1.457R) + damp 0.1R | 13246558 | 2026-08-20 | — | — | — | banner VERIFIED 2026-08-20 (policy truncate, damp 0.1R, depth 1.457R) | Phase 3 production candidate, 48 h |
| p022_ige_coarse | 56_57 | h/R 1, 4R/0.15R disc, trunc 3R, none | 13207682 | 2026-08-18 | 2026-08-19 (cancelled) | n/a — blew up at rev 17.6 | n/a | **FAIL (blow-up)** | Phase 2, 24 h; see below |

Deploy note (2026-08-18): driver+launcher md5-verified on cluster. src NOT
pushed: local carries uncommitted 021 solver WIP; the cluster src (live 018
state) already has everything 022 needs, and the tuple block-GS `solve!` was
diff-verified byte-identical local↔remote. Remote Manifest pins 1.11.7 ✓.
Queue at submission: 18 active (12 R) — 022 adds 4.

## p022_ige_coarse blow-up (2026-08-19)

Cancelled at 16 h 28 m after diverging at step 636/1007 (rev ≈ 17.64). Signature:

| quantity | steps 240–630 (healthy) | step 630 | step 636 |
|---|---|---|---|
| CF_x | −0.074 ± 0.004 | — | **−371.9** (then +2.32) |
| ground RMS(U·n) | 0.12 → 0.22, flat after step 340 | **0.842** | — |
| ground max\|U·n\| | 3.4 → 5.4 | **20.17** | — |
| below-ground sum\|Γ\| | 0.31 → 0.38 | **9.859** | — |
| below-ground count | 456 → 6362 (monotone) | 6813 | — |
| walltime/step | 100–160 s | 189 s | 413 s → 633 s |

Only steps 636–637 are corrupted; everything up to 635 is clean and was used for
the preliminary read below. The below-ground census grew monotonically from
step 240 under the leave-be policy (ruling 2) while sum|Γ| stayed ~0.35 — so
accumulation alone was tolerated for ~390 steps and then something tipped. The
leading suspect is the below-ground particle population; this is the trigger the
contingency chain names for moving Phase 4's cull ablation earlier. Not yet known
to be mesh-independent: p022_ige_fine (13207681) was at step 445 with comparable
tangency (RMS ≈ 0.24, 2714 below-ground) and no excursion when checked.

## Preliminary matched read (2026-08-19, NOT an M1 harvest)

Revs 7–17, entirely post-freestream-withdrawal (hover starts rev 6.5) and
entirely pre-blow-up, but inside the settling transient (`SETTLE_REVS=20`
intends revs 18–28). Computed by `scripts/p022_harvest.py ratio`:

| run | CT̄ | cycle-std | SE of mean |
|---|---|---|---|
| p022_oge_coarse | 0.075516 | 0.003093 (4.10%) | 0.000978 (1.30%) |
| p022_ige_coarse | 0.078888 | 0.005376 (6.82%) | 0.001700 (2.16%) |

Ratio IGE/OGE = **1.045** ± 0.083 (cycle-std propagated) / ± 0.026 (SE
propagated). Right direction, consistent with both the momentum-theory anchor
(1.067) and experiment (1.078 ± 0.008), but not discriminating. See the
resolvability analysis in `phase_05_hr_sweep.md`.

## Headline (Phase 2, fill at harvest)

| rung | CT_OGE | CT_IGE | ratio IGE/OGE | anchor (momentum theory ≈1.07) |
|---|---|---|---|---|
| 45_185_ct4 | 0.07480 ± 0.00152 | 0.07934 ± 0.00234 | **1.061 ± 0.038** (cycle-std propagated) | consistent (1.067; experiment 1.078 ± 0.008) |
| 56_57 | walled @ rev 25 (no M1 window) | blew up rev 17.6 | prelim 1.045 ± 0.083 (revs 7–17, NOT M1) | — |

Fine-rung numbers are cycle-means over the final 10 in-hover revs from
`*_case_metadata.toml` (harvested 2026-08-21). Both fine runs fail the strict
per-rev criterion (`converged=false`) but the cycle-mean is the headline per
`decision_rules.md`. Phase-2 exit still needs the clean-context verify pass +
Ryan's read (and note the coarse rung contributes nothing harvestable).

## Storage cycle 2026-09-21 (hpc-storage watchdog)

Two passes on `/home/rander39` (400 G cap), both APPLY mode, byte-for-byte
verified before any delete, newest 5 restartable steps left on `/home`.
Protect list, CSVs, TOMLs, `monitors/`, and the job queue untouched.

Pass 1 (routine sweep) — archived 2 finished 032 runs from the
`projects_FLOWPanel.jl` checkout to
`/nobackup/archive/usr/rander39/FLOWPanel_runs/projects_FLOWPanel.jl/`:

| run | src MB | tar MB | steps kept | freed MB |
|---|---|---|---|---|
| `scr_p032om15_ctrllg_fs` | 15934 | 12019 | 463–467 | 15643 |
| `scr_p032om15_explg_fs_oldlaw` | 15987 | 12341 | 463–467 | 15691 |

`/home` 287.3 G → 256.9 G (freed 30.4 G). `VERIFY_FAIL=0`, `STALE=0`,
`CHECKOUT_LOCKED=0`. No VTK sweep needed. 20 candidate directories deduped by
`realpath(data/)` to 3 real checkouts (`projects/FLOWPanel.jl`,
`flowpanel-021/FLOWPanel.jl`, `wt052/FLOWPanel.jl`); the `wt018/`, `wt021/`,
`wt026/`, and GPU-clone paths share one `data/` tree and were ALIAS-SKIPped.

Pass 2 (Ryan-approved RECENT run, `--include-recent --only`) — archived
`p022lg_hr10_om15`: src 47544 MB → tar 37499 MB
(`p022lg_hr10_om15.tar.zst`, 39,320,612,716 B), newest 5 steps [1003–1007]
kept on `/home`, freed 47239 MB, `VERIFY_FAIL=0`, `ARCHIVED.txt` breadcrumb
written. `scr_p020r_geom_s020v_om15` (live job 13829223) left untouched.
`/home` 241469 → 216309 MB (235.8 G → 211.2 G). CAVEAT: the "before" `du` for
this pass overlapped the apply, so 235.8 G is not a clean pre-action baseline —
other live runs grew during the ~25 min window; the 47239 MB action figure is
the reliable one.

Archive quota after both passes: 4.549 T, 93,402 files. Two items flagged for
Ryan, no action taken: (1) `lfs quota` on `/nobackup/archive` now reports no
enforced block/file limit (`0k`), contradicting the 20 TiB / 1 M baseline noted
in the 018 ledger; (2) restarting `p022lg_hr10_om15` from any step before 1003
now requires unpacking the tarball — steps 1003–1007 restart directly off
`/home` as usual.
