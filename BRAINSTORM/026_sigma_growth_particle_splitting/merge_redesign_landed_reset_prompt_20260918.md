# 026 reset prompt — merge redesign LANDED (uncommitted); commit gate, then rerun slate (2026-09-18)

You are picking up BRAINSTORM 026 in
`/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (branch
**fastmultipole**) + sibling `/Users/ryan/Dropbox/research/projects/FLOWVPM.jl`
(branch **flowpanel**). Read `CLAUDE.md` and the policies it names.
Predecessor: `merge_sigma_redesign_reset_prompt_20260917.md` — its tasks
are DONE. This file supersedes it.

## State: what is implemented and verified (all UNCOMMITTED)

The §22 merge redesign (`particle_splitting_design.md` §22, incl. the
2026-09-18 implementation-review rulings block at the end of §22.2) is
fully implemented in FLOWVPM and verified. Durable record:
`p026_derisk_20260914_provenance.md` §"Merge-σ redesign implementation
+ smoke (2026-09-17/18)" — read that section for the full numbers.
Short form:

- **§22.1**: merged σ = second-moment rule (σ_new² = ⟨σ²⟩_w +
  (1/3)⟨|xᵢ−x̄|²⟩_w, w=|α|), in `_finalize_merged_particle!`
  (`FLOWVPM.jl/src/FLOWVPM_merging.jl`; the new `ext` loop computes all
  weighted means about the actual placement). Coincident equal pair →
  σ (pump gone); equal pair at d → σ² + d²/12.
- **§22.2**: merge-path `_rsplit_reset_slot!` replaced by
  `_rsplit_merge_lineage!` (`FLOWVPM.jl/src/FLOWVPM_resolution_split.jl`):
  σ₀²/dvisc/drvpm := same |α|-weighted means, separation term
  (1/3)⟨|Δx|²⟩_w added to drvpm. PLUS the 2026-09-18 axis ruling:
  axis/weight := sign-aligned |α|-weighted means (no more zeroing).
  Host + device dispatch mirrors `_rsplit_reset_slot!`
  (`_rsplit_slot_value` readers, view-broadcast writes; no scalar
  indexing on device storage). Split-path resets untouched.
- Rulings 2026-09-18: axis lineage = DO IT (done); merged vol = LEAVE
  as Σvol (traced: vol inert in production dynamics — integrator/SFS
  use σ³ directly; only PSE/smooth-conversion/filament ops read it);
  zero-total-|Γ| fallback = reviewed, fine.
- **Tests green**: FLOWVPM merging 78/78, resolution-split 1058/1058,
  filament-edge-graph 477/477 (run `julia --project=test --threads 4
  test/runtests_<x>.jl` from FLOWVPM.jl); FLOWPanel
  unit_wake/unit_replay/unit_simulate exit 0. Three old cbrt-law
  assertions were updated to the new law (expected, documented in the
  provenance append). Known UNRELATED pre-existing failure:
  FLOWPanel `runtests_unit_warmstart.jl` first testset — ignore.
- **Smoke PASSED** (`data/smoke_mergesigma2m_20260917/` in FLOWPanel,
  467 steps, wave-2 fs env, exit 0): 4,763 merge events; per-event σ
  growth median 0.24% (old law on same pairs: 25.8%); max σ ended at
  1.109× shed σ vs wave-2 deaths at 2.5–3.4×. Note: smoke ran
  PRE-axis-lineage (σ law + ledger only); axis change is
  direction-state only and fully unit-tested — do not rerun the smoke
  unless Ryan asks.

## Your task 1 — commit gate (Ryan-gated: prepare, ask, then commit)

When Ryan approves, commit as a bundle:
- FLOWVPM (branch `flowpanel`): `src/FLOWVPM_merging.jl`,
  `src/FLOWVPM_resolution_split.jl`, `test/runtests_merging.jl`,
  `test/runtests_resolution_split.jl`, `test/runtests_filament_edge_graph.jl`.
  (Untracked `examples/p026_ring_split_test_out/` is NOT ours — leave.)
- FLOWPanel (branch `fastmultipole`), the 026 doc bundle:
  `particle_splitting_design.md` (§22 + §22.2 rulings block),
  `p026_derisk_20260914_provenance.md` (Wave-2 appends + merge-redesign
  append), `ledger.md`, `wave2_inflight_reset_prompt_20260915.md`,
  `merge_sigma_redesign_reset_prompt_20260917.md`, this file. Do NOT
  sweep in other arcs' in-flight files (018 gamma_distribution_*/
  mechanism_tests_*/rlxfscaled_*, 021 files, expguard_provenance,
  examples/run_rotor_multi_ground_effect_gpu.slurm.sh, the
  rotor_hover_pressure_comparison metadata.toml).
- The smoke run dir `data/smoke_mergesigma2m_20260917/` (~VTK-heavy):
  ask Ryan whether to keep/commit CSVs only/delete; default = leave on
  disk uncommitted.

## Your task 2 — after the commit: re-derive the §22.3 rerun slate

Parked until the redesign landed; now unblocked. s9 caps were specified
at k≈2–3 shed multiples to compensate for the σ-pump — the smoke shows
the pump is gone (max σ 1.109× shed unforced), so re-derive from the
new-law behavior rather than reusing k blindly. Slate design needs
Ryan's approval before any HPC submission; campaigns need worktrees +
annotated tags per global CLAUDE.md.

## Owed / parked (carried)

- Notebook entry for the whole 026 arc still owed (Ryan "not yet" ×2) —
  offer, don't write.
- orc `fastmultipole`/`flowpanel` branch-divergence ruling (tags carry
  the campaign pins; do NOT force-push).
- scr_p026gpuv_split archive retry via hpc-storage once quiet ≥24 h.
- 021 silo cleanup owed (021 package §7); 018/022 queue jobs are other
  sessions'.

## Ground rules (carry-over)

Local ≤4 threads. `ssh orc` needs a live ControlMaster socket (ask Ryan
to run `! ssh orc echo ok` if 2FA blocks). MOTD contaminates ssh output
— filter. Judge runs by outputs, never sacct state. For any ORC
execution read `BYU_ORC_AGENTS.md` first. Commits are Ryan-gated.
