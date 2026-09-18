# 026 reset prompt — implement 2nd-moment merged σ + coincident-limit growth lineage (2026-09-17)

You are picking up BRAINSTORM 026 in
`/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (branch
**fastmultipole**) + sibling `/Users/ryan/Dropbox/research/projects/FLOWVPM.jl`
(branch **flowpanel**). Read `CLAUDE.md` and the policies it names.
Predecessor: `wave2_inflight_reset_prompt_20260915.md` — its tasks are DONE
(banner verification, harvest, data moved to shared root; record in
`p026_derisk_20260914_provenance.md` §Wave-2-banner / §Wave-2-harvest).
This file supersedes it.

**Your task: implement the two Ryan-ruled merge-system changes in design
doc §22 (`particle_splitting_design.md`, end of file — READ §22 FIRST;
also §19 for trigger semantics, §21 for the merge-A/B settings).** Both
changes are in FLOWVPM. Do NOT launch any HPC runs; the rerun slate is
parked until this lands (§22.3). Commits are Ryan-gated.

## Why (one paragraph)

Wave-2 autopsy: the volume-conserving merged σ (`cbrt(σᵢ³+σⱼ³)`) grows
+26% per equal pair regardless of overlap — a pure-bookkeeping σ-pump
(coincident identical particles should merge to σ, not 1.26σ). And the
merge path resets the resolution-split ledger (σ₀ := merged σ,
accumulators := 0), so merge-driven growth is invisible to the split
triggers. Together these produced the ctrl-family blow-ups (merge frenzy
~1500 pairs/step, np 246k→54k collapse, CT→270) and the exp-family σ
growth to the FMM adequacy limit. Full evidence: provenance §Wave-2.

## Change 1 — second-moment merged σ

File: `FLOWVPM.jl/src/FLOWVPM_merging.jl`. Current rule at :166
(`sigma = cbrt(sigma3_sum)`, accumulated at :241 and in the
acc[18]/acc[20] path ~:296–:337). Merging is pairs-only per pass
(verified: `paired[]` marking at :525/:537/:576-577; every union-find
root has exactly 2 members) — but implement over the existing
cluster-accumulation structure so it stays correct for n members.

New rule (weights w_i = |α_i|, the vector-strength magnitude; x̄ = the
representative's final position — compute the variance about where the
merged particle is actually placed):

    sigma_new^2 = ⟨σ²⟩_w + (1/3)·⟨|x_i − x̄|²⟩_w

Equal pair at distance d: σ_new² = σ² + d²/12 (exact unit-test value).
Coincident limit: σ_new = σ (pump gone). Keep unchanged: pair gates
(`max_sigma_ratio=2.0`, `gamma_align_cos`), strength/position
combination, merge_events.csv schema (`step,np,sigma_i,sigma_j,dist`),
the σ-weighted scalar circulation at :167 (unless it references
sigma3_sum — then adapt minimally and note it).

Note: check what weighting the existing code uses to place the
representative (centroid). Do NOT change particle placement; if it is
not |α|-weighted, still compute the new σ variance about the actual
placement and FLAG the inconsistency in your report for Ryan.

## Change 2 — coincident-limit ledger lineage on merge

Files: `FLOWVPM.jl/src/FLOWVPM_merging.jl` :184–:190 (the
`_rsplit_reset_slot!` call — this reset is what you are replacing on
the MERGE path only) and `FLOWVPM.jl/src/FLOWVPM_resolution_split.jl`
(state: `sigma_0` + per-mechanism attempted-Δσ² accumulators, e.g.
`dvisc` ~:78, drvpm; read the state struct first).

On merge, instead of resetting, set every ledger line to the SAME
|α|-weighted mean over members, with the separation term omitted:

    sigma_0,new² = ⟨σ₀²⟩_w        Δσ²_k,new = ⟨Δσ²_k⟩_w   for every mechanism k

then ADD the separation term (1/3)·⟨|x_i−x̄|²⟩_w to `drvpm` (grow side —
Ryan's default ruling §22.2; it routes merge-coarsening to the
compress/tri3 trigger; do not add a new accumulator). Everything is
σ²-additive with one weight definition, so the ledger identity
σ² ≈ σ₀² + ΣΔσ² survives the merge exactly, and the invariants hold:
merging equals is maturity-neutral; merging unequals inherits
weighted-mean maturity; only genuine coarsening advances the split
clock. Split children keep their existing fresh-σ₀ reset (correct —
do not touch the split path). Anti-livelock fence is pre-existing
(child overlap 2.4 < Φ_merge 3.5), no action.

GPU caveat: the resolution-split state has a device-backed branch (see
`enable_resolution_split!`); the GPU splitting work (2026-09-11
directive) added device accumulators. Make the merge-path lineage work
for BOTH array types the way `_rsplit_reset_slot!` already does —
follow its dispatch pattern; no scalar indexing on device storage.

## Tests (local, ≤4 threads)

Existing (all must stay green; run FLOWVPM first):
- FLOWVPM: merging testset (was 50/50), resolution-split (1055/1055),
  filament-edge-graph (477/477)
- FLOWPanel: `test/runtests_unit_wake.jl` (730), `runtests_unit_replay.jl`
  (142), `runtests_unit_simulate.jl`
- Known UNRELATED pre-existing failure: `runtests_unit_warmstart.jl`
  first testset — ignore. Some merging tests may assert the OLD cbrt
  law — update those asserted values to the new law (that is expected,
  not a regression; say so explicitly in your report).

New unit tests to add (FLOWVPM merging testset):
1. Coincident equal pair → σ_new == σ (tolerance eps-level).
2. Equal pair at distance d → σ_new² == σ² + d²/12.
3. Maturity neutrality: two particles with σ₀ set so each is at growth
   ratio g merge → merged σ/σ₀ == g when coincident.
4. Separation credit: merged drvpm gained exactly (1/3)⟨|Δx|²⟩_w.
5. Unequal-w pair: all merged ledger lines equal the |α|-weighted means.

Then one end-to-end smoke: the local CPU driver smoke used for wave-2
telemetry (467-step `rotor_hover_pressure_comparison.jl` config — see
provenance §Wave-2 telemetry; MERGE_EVENT_LOG on) — verify it runs,
merge events log, and report the max σ census vs the wave-2 baseline
(expect max σ growth visibly slower).

## Report / gates

- Report: diffs summary, test tallies, smoke σ-census comparison, any
  flagged inconsistencies (centroid weighting, circulation rule).
- Commits are Ryan-gated: prepare but do NOT commit. When Ryan
  approves, the bundle should also include the still-uncommitted 026
  docs: provenance §Wave-2 appends, ledger lines, design-doc §22,
  wave2_inflight + this reset prompt.
- Notebook remains owed for the whole 026 arc (Ryan "not yet" ×2) —
  offer, don't write.

## Owed / parked (carried)

- orc `fastmultipole`/`flowpanel` branch-divergence ruling (see
  provenance NOTE; tags carry the campaign pins).
- Rerun slate (s9 with caps in shed-σ multiples k≈2–3 + cap ladder) —
  re-derive AFTER this redesign lands (§22.3).
- scr_p026gpuv_split archive retry via hpc-storage once quiet ≥24 h.
- 021 silo cleanup owed (021 package §7); 018/022 queue jobs are other
  sessions'.

## Ground rules (carry-over)

Local ≤4 threads. `ssh orc` needs a live ControlMaster socket (ask Ryan
to run `! ssh orc echo ok` if 2FA blocks) — but this task is local-only.
MOTD contaminates ssh output — filter. Judge runs by outputs, never
sacct state. For any ORC execution read `BYU_ORC_AGENTS.md` first.
