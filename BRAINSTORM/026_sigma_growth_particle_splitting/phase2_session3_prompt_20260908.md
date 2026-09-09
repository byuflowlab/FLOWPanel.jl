# 026 Phase 2 Session 3 prompt — ship resolution splitting AS the FLOWVPM splitting system; remove the legacy path

Written 2026-09-08 at the end of Session 2 (commits 4–6). Copy-paste the
prompt below to the next agent. Authorization on record: Ryan 2026-09-08 —
"ship this to FLOWVPM (not FLOWPanel). We can remove existing particle
splitting code in favor of this." That supersedes the Phase-2 byte-identity
rule for `src/FLOWVPM_splitting.jl` (which existed to protect the old path
only while the new one was being built).

---

## Prompt

Ship BRAINSTORM 026 resolution splitting as THE particle-splitting system of
FLOWVPM (`~/Dropbox/research/projects/FLOWVPM.jl`, branch `flowpanel`), and
REMOVE the legacy experimental splitting path in its favor (Ryan authorization
2026-09-08). Read FIRST, in order: (1) this file's "State inventory" section
below — it is the map of every file you will touch and was verified against
the working trees on 2026-09-08; (2)
`BRAINSTORM/026_sigma_growth_particle_splitting/phase2_impl_handoff_20260907.md`
(FLOWPanel repo) — Session 1 + Session 2 reports document the new system's
API and seams; (3) the plan
`~/.claude/plans/shimmering-wobbling-sunrise.md` for design rationale only —
do NOT re-open its closed decisions.

The deliverable is FLOWVPM commits. FLOWPanel edits are allowed ONLY as
minimal companion commits that delete legacy consumers so the two-repo dev
stack still builds (FLOWPanel dev-paths FLOWVPM, so a FLOWVPM-only removal
would break it) — no new FLOWPanel features, and list every companion commit
separately in your report.

Scope, in commit order (gate each before the next):

1. **FLOWVPM: first-class integration.** Make the new system usable
   standalone: wire an optional resolution-split step into `run_vpm!`
   (pattern-match the existing `merge_every`/`merge_kwargs` kwargs: e.g.
   `split_every::Int=0`, `split_opts=nothing`; apply after merging, order per
   W3) with an on-merge state reset for standalone users (see decision D-A
   below). Keep exports (`split_particles!`, `ResolutionSplitOpts`,
   `ResolutionSplitState`, `enable_resolution_split!`) and add docs to the
   `split_particles!` docstring noting it is the ONLY splitting system.
   Gate: `runtests_resolution_split.jl` green + a new small testset for the
   `run_vpm!` wiring.

2. **FLOWVPM: legacy removal.** Delete the legacy path:
   - `src/FLOWVPM_splitting.jl` (2659 lines: `SplittingState`,
     `SplitOptions`, `SplitDirection`, all `SplitTrigger`s, `_do_split!`,
     `accumulate_H_chi!`, `should_split`/`severity`) and its `"splitting"`
     entry in the include loop (`src/FLOWVPM.jl:91`).
   - The `splitting_state::SplittingState{R}` field on `ParticleField`
     (`FLOWVPM_particlefield.jl:248`) + its ctor init + lockstep hooks
     (:375, :665) + the `accumulate_H_chi!` call in `nextstep` (:789).
   - The legacy W2 `dsigma2_visc/rvpm` writes at the six sites
     (`FLOWVPM_viscous.jl:167,196,218` and
     `FLOWVPM_timeintegration.jl:575,1118` + the euler site near :289) —
     **KEEP the colocated `_rsplit_accumulate_dsigma2!` mirror calls**, they
     are the new system's attribution and are separate statements.
   - The legacy wholesale reset block in `_finalize_merged_particle!`
     (`FLOWVPM_merging.jl:187-194`) — see decision D-A.
   - Legacy exports (`SplitOptions`, `SplitDirection`, `STRENGTH`,
     `STREAMLINE`, `STRAIN1`, `SplitTrigger`, `AllTrigger`, `AnyTrigger`,
     `HoldTrigger`, `GammaMagTrigger`, `ZTrigger`, `StretchTrigger`,
     `SeparationTrigger`, `SigmaShrinkTrigger` — `src/FLOWVPM.jl:48-53`;
     keep `split_particles!` itself).
   - `test/runtests_dsigma2_accumulators.jl` (tests the legacy accumulators;
     port any coverage of the NEW mirrors into
     `runtests_resolution_split.jl` first — do not lose dvisc/drvpm
     attribution coverage) and its include in `test/runtests.jl:21`.
   - Grep-audit afterwards: `SplittingState|SplitOptions|H_chi|hold_counter|
     cooldown_counter|dsigma2_` must have zero live references in
     `FLOWVPM/src` and `FLOWVPM/test`.
   Gate: FULL FLOWVPM suite green (`julia --project=test test/runtests.jl`
   from the FLOWVPM root — `--project=.` trips TestEnv; max 4 threads).

3. **FLOWPanel companion (minimal deletions only, branch `fastmultipole`).**
   Delete the legacy consumers so FLOWPanel builds against the new FLOWVPM:
   - `SplitParticles` policy struct + its `apply_particle_policy!` +
     the mutual-exclusion guard mentions (`src/FLOWPanel_wake.jl` :1644ff,
     :1730ff area — the ResolutionSplit policy and its
     `_resolution_split_merge_hook`/`_heal_unseeded_rsplit_slots!` STAY).
   - The `FLOWVPM.accumulate_H_chi!(w.pfield, dt)` call
     (`FLOWPanel_wake.jl:2351`).
   - The `split_*` six-field VTP writer block (`_write_particles_vtp`
     `split_state` kwarg + block, `FLOWPanel_wake.jl` ~:2437-2444) and the
     `split_state=w.pfield.splitting_state` argument at the call site — the
     new `rsplit_*` block stays.
   - The `split_*` loader block + `splitting_state` lines in
     `_clear_splitting_state!` (`FLOWPanel_warmstart.jl` — keep the
     filament-edge-graph clearing and all `rsplit_*` code). NOTE the
     warm-start COMPATIBILITY decision D-B below before touching this.
   - `splitting_state` in the GPU side-buffer copy list
     (`FLOWPanel_gpu_wake.jl:~66-72`; `filament_edge_graph` stays).
   - `SplitParticles` replay serialization branch (`FLOWPanel_replay.jl`).
   - Tests: the 38-test "warm start SplittingState persistence (026 W1)"
     testset (`runtests_unit_warmstart.jl:440-585`), the
     `SplitParticles(nothing)` mutual-exclusion line in the new
     ResolutionSplit testset (`runtests_unit_wake.jl` — reduce that test to
     construction-only or drop just that assertion), and any `SplitParticles`
     mention in `runtests_unit_replay.jl`.
   Gate: wake + replay suites green; warm-start 026 testsets green run in
   isolation (the file aborts earlier at the KNOWN pre-existing Julia 1.12
   WeakKeyDict/`WarmstartNoopSolver` failure — not yours); IGE suite
   113/113 in isolation; splitting-off short-march smoke vs pre-removal
   HEAD bit-identical (Session 2 report describes the exact recipe:
   scratchpad cwd, `NREVS=0.25 FREESTREAM_RAMP_REVS=0.1
   FREESTREAM_HOLD_REVS=0.05 FREESTREAM_WITHDRAW_REVS=0.05
   SETTLE_REVS=0.05`, compare md5 of all outputs; only the wake-health
   `wall_s` column may differ).

Decisions you must resolve (recommendations on record; ask Ryan only if you
disagree):

- **D-A (merge reset ownership).** With the legacy wholesale reset deleted
  from `_finalize_merged_particle!`, standalone FLOWVPM users calling
  `merge_particles!` would get stale `ResolutionSplitState` on merged
  representatives. RECOMMENDED: add a one-branch guarded reset in
  `_finalize_merged_particle!` (`rs = pfield.resolution_split; rs ===
  nothing || _rsplit_reset_slot!(rs, representative, sigma)`) and KEEP the
  `on_representative` hook (other consumers may use it). This supersedes the
  2026-09-07 "merging stays ignorant" ruling, which protected old-path
  byte-identity that no longer applies; the FLOWPanel closure becomes
  redundant-but-harmless (leave it, or drop it in the companion commit —
  either way say which in the report).
- **D-B (checkpoint compatibility).** Existing on-disk checkpoints carry
  `split_*` fields (e.g. the `scr_p026ef_*` warm-start brackets). After
  removal, the loader must IGNORE unknown `split_*` fields silently (they
  become inert extra point data) rather than erroring — verify a
  Session-2-era checkpoint still loads (the 026 W1 testset's writer can
  fabricate one before you delete it). Do not write `split_*` anymore.

Ground rules:
- FLOWVPM branch `flowpanel` currently at `65247ee` (on top of Ryan's
  `119fe23` merge-runaway guard and `21eeaaa` euler_exp sigma_guard — his
  runtests_merging additions must stay green); FLOWPanel `fastmultipole` at
  `e663e59`. Both live checkouts are shared with other agents — `git status`
  before committing, stage ONLY your files, leave unrelated dirty files
  alone.
- The new system's files: `src/FLOWVPM_resolution_split.jl` (included BEFORE
  `particlefield` — struct-field dependency; keep it that way),
  `test/runtests_resolution_split.jl` (812 tests),
  `examples/p026_ring_split_test.jl`. Do not rename any of the new API.
- Max 4 local threads. Commits 7–8 of the Phase-2 plan (campaign launches,
  SIGMA_CEIL removal) stay Ryan-gated — untouched.

Deliverable: FLOWVPM commits (+ minimal FLOWPanel companion commits, listed
separately) with gates green; per-commit summary, the D-A/D-B outcomes, a
grep-audit statement, and exact SHAs appended to
`BRAINSTORM/026_sigma_growth_particle_splitting/phase2_impl_handoff_20260907.md`
under "## Session 3 report". Offer — do not write — a notebook entry.

---

## State inventory (verified 2026-09-08)

New system (KEEP, FLOWVPM): `src/FLOWVPM_resolution_split.jl` —
`ResolutionSplitState{R}` (sigma_0/axis/weight/exposure/dvisc/drvpm),
`ResolutionSplitOpts{R}` (hand-written kwarg ctor, no `@kwdef`),
`enable_resolution_split!`, lockstep `_rsplit_{init,reset,swap,zero}_slot!`,
integrator-inline `_rsplit_accumulate!` + `_rsplit_accumulate_dsigma2!`
mirrors, `_rsplit_direction`, kernels `_split_viscous_tetra4!` /
`_split_compress_tri3!` / `_split_elongate_pair2!`,
`split_particles!(pfield, ::ResolutionSplitOpts; verbose=false, dt=nothing)`
(dt accepted+ignored), `_radix_log_sigma_adequacy` in `FLOWVPM_fmm_radix.jl`.

New system (KEEP, FLOWPanel): `ResolutionSplit{TO}` policy +
`_resolution_split_merge_hook` + `_heal_unseeded_rsplit_slots!`
(`FLOWPanel_wake.jl`), `rsplit_*` VTP writer block + all-or-nothing loader +
`_zero_resolution_split_state!` (`FLOWPanel_wake.jl`/`FLOWPanel_warmstart.jl`),
GPU-seam doc in `_gpu_copy_side_buffers!`, ResolutionSplit replay
serialize-then-drop, `WAKE_SPLIT_*` driver knobs + `scr_p026sp_*`/
`scr_p026s9_*` dispatcher arms, tests in `runtests_unit_wake.jl` (17),
`runtests_unit_warmstart.jl` ("ResolutionSplitState persistence", 18),
`runtests_unit_replay.jl` (6), opt-in CUDA seam testset.

Legacy system (REMOVE): see scope items 2–3 above; grep counts as of today —
FLOWVPM: `FLOWVPM_splitting.jl` (whole file), `FLOWVPM_particlefield.jl`
:248/:375/:665/:789, `FLOWVPM_merging.jl` :187-194, `FLOWVPM_viscous.jl`
:167/:196/:218/:282, `FLOWVPM_timeintegration.jl` :289/:575/:1118,
`FLOWVPM.jl` exports :48-53 + include :91, `test/runtests_dsigma2_accumulators.jl`.
FLOWPanel: `FLOWPanel_wake.jl` (11 refs), `FLOWPanel_warmstart.jl` (6),
`FLOWPanel_replay.jl` (3), `FLOWPanel_gpu_wake.jl` (1),
`runtests_unit_warmstart.jl` (12), `runtests_unit_wake.jl` (1),
`runtests_unit_replay.jl` (1).

Session 2 SHAs (FLOWPanel): housekeeping `9d2c0f1`, policy `c78f6d1`,
persistence `ebbbf3b`, driver/dispatcher `e7681d2`, report `e663e59`.
Session 1 SHAs (FLOWVPM): `9d63578`, `99f4d54`, `5583443`; ring test
`65247ee`.
