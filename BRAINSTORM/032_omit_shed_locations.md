# 032 — Omit shed locations from a panel object

Opened 2026-09-18 (Ryan directive, via 026 §22.3 rerun-slate autopsy). Status:
**feature implemented, unit-tested; A/B reruns pending Ryan-gated submission.**

## Motivation

Ryan (2026-09-18, verbatim intent): *make a new item to add a feature to omit
shed locations from a panel object, and then try doing that at the root-most
shed locations. I suspect the strong root particle strength is causing the
instability here, and it isn't even strictly physical.*

Evidence (026 rerun slate, `026_sigma_growth_particle_splitting/rerunslate_provenance_20260918.md`):

- A1 (`scr_p026s9r2_explg_fs`, job 13758586) died step 328/467 at the
  euler_exp substep guard (dt·|L| > 2048). Offender localized from the
  step-327 VTP: depth x = −0.37R (thrust side), r = 0.344R — the fountain-flow
  recirculation region. Driver/victim: a σ-at-floor particle with |Γ|=5.4e-3
  (Γ/σ² = 9.1e4, field max) imposing the gradient on near-zero-Γ neighbors
  1.5σ away. Global field healthy (Ryan, ParaView).
- A2 (floor-only, 13763819) died step 274 — the exact step of its wave-2 twin,
  so this ignition channel is independent of the (fixed) merge σ-pump.
- Hypothesis: root-shed circulation feeds the fountain region; blade-root
  shedding is itself of dubious physicality (root cutout / hub interference
  in reality).

## What already existed

`examples/rotor_hover_pressure_comparison.jl` has since the 018 TE-trace
rework carried an explicit **modeling root clip**: trace the full TE
(anchored outer→inner by `end_node`), then `clip_shedding_root` drops
shedding-matrix columns whose edge midpoint has |r|/R < `SHEDDING_R_OVER_R`
(env knob, default **0.1** for 018 comparability). On the stock DJI blade the
TE runs from r/R 0.111 outward, so the 0.1 clip is inert for the *bladed* TE —
but it does clip the true-root TE segment down at r/R 0.0095, i.e. edge
omission is already exercised in every 018/026 production run. What was
missing: the filter lived driver-side (radial-only, untested) and the knob was
recorded nowhere in run metadata.

## Feature (implemented 2026-09-18)

Package-level, general omission in `src/FLOWPanel_liftingbody.jl` (exported):

- `shedding_edge_midpoint(nodes, cells, shedding, j)` — midpoint of column
  `j`, resolved through the shedding panel's cell-local node slots.
- `filter_shedding(shedding, keep::Vector{Bool})` — column mask form.
- `filter_shedding(keep_edge, nodes, cells, shedding)` — predicate form,
  `keep_edge(midpoint, j) -> Bool`.

Both return a fresh `Matrix{Int}` with retained columns in order. Driver's
`clip_shedding_root` now delegates to the predicate form with the identical
criterion (edge-midpoint |r|/R ≥ cutoff) — stock behavior bit-identical.

Winding invariant preserved by construction: the filter is applied *after*
`calc_shedding_from_seed` on the constructed base body's `.nodes`/`.cells`
and before the rebuild, exactly where the clip already sat (CLAUDE.md
critical invariant). Omitted edges cannot silently reappear — the rebuilt
body's `shedding`/`shedding_full`/`Das`/`velocity_te` are all sized from the
filtered matrix (unit-tested).

Plumbing:

- `SHEDDING_R_OVER_R` is the user knob (kept; no duplicate env var added).
  Now exported with default in `examples/run_p018_screen_hpc.slurm.sh` so
  sbatch `--export` A/Bs work, and recorded in the case-metadata TOML along
  with `blade_root_r_over_R` and traced/retained edge counts per blade. The
  driver banner already printed the requested clip and per-blade root
  midpoints.

## Feature B (implemented 2026-09-18, Ryan directive — PRIMARY A/B lever)

Ryan's refinement: *don't just move the root filament outboard — eliminate it.
Keep shedding Das (rigid panels) and the free wake-panel row, but don't shed
the inner-most particles.* Implemented as `OmitStations` in
`src/FLOWPanel_wake.jl`:

- `OmitStations(method, omit)` wraps any station-resolvable
  `WakeSheddingMethod`; `omit[i_surf][j]` masks wake-node column `j` (same
  vertex-based ordering as `StationSigmaOverlap`/`Das`: station j = edge j's
  nib node, station n+1 = last edge's nia node). Masked trailing filaments —
  including the chain-closing terminal filament — route to an accounting sink
  instead of shedding particles; an unsteady filament is masked only when
  BOTH its stations are. The solve, Das/Kutta closure, and wake-panel rows
  are completely untouched: the difference between arms is purely which
  vorticity enters the free particle field.
- Deleted circulation is measured, not hidden: cumulative
  `omitted_circulation[]` (Σ|Γ|·Δl) and `omitted_filaments[]` counters,
  printed post-march and written to the metadata TOML.
- Legacy conversion only (`SurfaceVorticityConversion` rejects explicit line
  policies, so the guard is structural). Wraps compose with
  `StationSigmaOverlap` (chord–σ co-scaling) transparently.
- Driver knob `PARTICLE_OMIT_ROOT_R_OVER_R` (default 0.0 = off,
  bit-identical): masks stations with |r|/R below the value, via the
  driver's existing vertex-based `station_radii`. Refuses armed-but-inert
  (0 stations) and fully-masked configurations (018 silent-clip history).
  Exported in `run_p018_screen_hpc.slurm.sh`; banner prints masked counts.

Why B over A for the mechanism test: edge omission (A) relocates the closing
filament to the new terminal station where bound Γ is generally *larger*;
station omission (B) deletes it. B's physics cost is explicit circulation
non-conservation at the handoff (shed vortex lines end at the omission
boundary — a ∇·ω source VPM tolerates but never repairs), justified as
hub/root-cutout diffusion of an unphysically concentrated root vortex;
precedent: FLOWUnsteady's `no_shedding_Rthreshold`. Feature A is retained as
general-purpose infrastructure.

## Conservation / bookkeeping notes (flagged, not hidden)

- **Omission relocates, not deletes, the chain-closing root filament.** In
  `_convert_to_particles!` (LegacyEdgeJumpConversion,
  `src/FLOWPanel_wake.jl`), a non-wrapping chain sheds a terminal trailing
  filament carrying the full wake-column circulation −Γ at the last station.
  Raising the clip moves that station outboard and its strength becomes Γ at
  the *new* terminal column. Since bound circulation generally rises from
  root cutout toward mid-span, the relocated root filament can be *stronger*
  than before; the mechanism bet is that releasing it further from the
  hub/fountain recirculation (and away from the smallest-chord, σ-at-floor
  stations) breaks the ignition feedback, not that the net root vortex
  vanishes. Kelvin is respected either way — the shed sheet still closes.
- Interior gaps (non-contiguous masks) would split one chain into pieces that
  each close with their own terminal filaments; if ever intended, pass the
  pieces as separate shedding matrices. The 032 use is a contiguous root clip.
- Kutta/jump closure: an omitted TE edge simply has no attached-wake influence
  and no wake column — the same regime every non-TE panel and the already
  clipped r/R < 0.1 stock-root segment lives in. No solver special-casing.
- `BoundCirculationMonitor` sizes its TE stations from the constructed body's
  shedding, so clipped stations drop out of `circulation_te` rather than
  reading as zeros; Γ(r/R) comparisons across arms must align stations by
  radius, not by column index.

## Validation

- Unit tests: `test/runtests_unit_liftingbody.jl`, testset
  "shedding omission (BRAINSTORM 032)" — midpoint resolution, identity mask
  (knob-off regression), mask and predicate forms, constructed-body
  consistency (`nsheddings`, `shedding_full == -1` on omitted panels, `Das`
  sizing), and error cases. Wing/plate regressions are covered by the
  existing suites (the filter is inert unless invoked).
- Driver regression: default `SHEDDING_R_OVER_R=0.1` retains the exact
  edge set of the previous inline clip (same criterion, same ordering).

## A/B plan (Ryan-gated submission)

Rerun the A1-class case (`scr_p026s9_explg_fs`, campaign
`campaign/p026-rerunslate-20260918` lineage) with root-most shedding omitted:

1. Local smoke first (few steps, banner shows the clip and edge counts;
   metadata TOML records it).
2. HPC arms (feature B is the primary lever): `PARTICLE_OMIT_ROOT_R_OVER_R`
   ∈ {0.15, 0.20} vs the 0.0 baseline (stock TE stations start at r/R 0.111,
   so these mask the innermost ~few stations *including the terminal
   root-closing filament* — record exact masked counts from the banner in
   provenance, plus the run-end deleted-circulation total). Optionally one
   `SHEDDING_R_OVER_R=0.2` arm (feature A) as a discriminator between
   "delete the root vortex" and "move the root vortex". Extend to A2
   (floor-only) if the first A/B moves the failure.
3. Success signal: fountain-region Γ concentration and the ~274–328-step
   guard trips disappear or move out; CT change small and explainable
   (record the ΔCT as a modeling-clip cost the way 018 budgets carry the
   relaxation term).
4. Campaign rules apply if it graduates beyond a smoke: worktrees, annotated
   tags, provenance file, Manifest dev-paths at the worktrees.

## Cross-references

- 026 `rerunslate_provenance_20260918.md` (A1/A2 autopsies, slate table),
  `root_shed_omission_reset_prompt_20260918.md` (this item's charter).
- 018 TE-trace rework comment block in
  `examples/rotor_hover_pressure_comparison.jl` (~line 340) — why trace and
  clip are separate jobs; cap-wrap failure history.
- CLAUDE.md critical invariant: shedding from the *constructed* body's cells.

## Log

- 2026-09-18: item opened; `filter_shedding`/`shedding_edge_midpoint`
  implemented + exported; driver clip refactored onto them (bit-identical
  default); knob added to launcher env + metadata TOML; unit testset added.
  A/B pending Ryan gate.
- 2026-09-18: unit suite PASS (liftingbody 56/56, 14 new). Local smoke PASS
  (40_40 mesh, NT=36, 5 steps each): clip 0.1 → 37→35 edges/blade, root
  midpoint r/R 0.124; clip 0.2 → 37→31 edges/blade, root midpoint r/R 0.225;
  both runs stepped cleanly, metadata TOML records all new keys. Scratch
  outputs in `data/scratch_p032_smoke{A,B}` (disposable).
- 2026-09-18 (later): Feature B (`OmitStations` particle-side omission)
  implemented per Ryan's refinement; wake + liftingbody unit suites PASS
  (755 + 56, incl. new OmitStations testset: all-false-mask bit-identity vs
  same-method reference, terminal-filament deletion = exactly 2 particles at
  station y=3 with Σ|Γ|·Δl = 0.5, boundary-spanning unsteady still sheds).
  First failure was a test bug (compared against the OverlapPPS golden
  reference from a SigmaOverlap fixture), not a feature bug. Driver smokes
  C (knob 0.15) / D (knob off) launched — verify results in
  `data/scratch_p032_smoke{C,D}` if this session ended before they reported.
  Rerun charter → `032_reset_prompt_20260918.md`.
- 2026-09-18 (rerun session): predecessor's smokes C/D both died with the
  session pre-post-march; rerun to completion. GOTCHA: `NREVS` alone does
  not shorten a smoke (`required_revs = max(nrevs, schedule_revs)`, schedule
  = 13 revs) — zero `FREESTREAM_{RAMP,HOLD,WITHDRAW}_REVS`/`SETTLE_REVS`
  too. smokeC (knob 0.15) PASS: banner 2/36 stations per blade, totals
  Σ|Γ|·Δl = 1.593e-3 m³/s / 12 filaments, `particle_omit_*` keys in
  `*_case_metadata.toml` (NOT `*.metadata.toml`). smokeD (off) PASS: no
  banner, knob 0.0 recorded, no station keys. Committed `af92740`; tag
  `campaign/p032-rootomit-20260918` pushed to orc; worktree + env at
  `orc:~/campaigns/p032-rootomit-20260918/` (FLOWVPM/FMM pins unchanged,
  p026 worktrees reused). B15=13771353, B20=13771354 submitted m13h,
  both RUNNING 22:32. Provenance: `032_rerun_provenance_20260918.md`.
- 2026-09-18 (late): **BOTH ARMS RESCUE A1.** B15 and B20 completed
  467/467 with zero guard trips (~71/75 min wall, np peaked ~370k vs A1's
  500k-cap ride), CT converged 0.07162/0.07132 (Δ −0.42% for 48% more
  deleted circulation: 0.486 vs 0.722 m³/s) — the fountain Γ-ignition is
  gone and the result is insensitive to clip width. A2 extension B-floor
  (`scr_p032om15_explg_floor`, knob 0.15, eng H200) = job 13772015,
  PENDING behind the 026 slate arms. Banner + outcomes tables in the
  provenance file.
- 2026-09-19: **FINAL VERDICT — omission rescues every tested class.**
  B-floor (A2-config) also COMPLETED 467/467 (CT 0.071654 ± 0.014%,
  deleted 0.4865 m³/s, min σ never below 4.6× floor). Fountain
  quantitative before/after: B15 max Γ/σ² flat ~46 through A1's death
  window (A1: 826→9.1e4). ΔCT vs cap018 −3.15% (no clean twin — caveats
  in provenance). ParaView side-by-side staged locally:
  `~/scr_p026s9r2_explg_fs_last50steps/` vs
  `~/scr_p032om15_explg_fs_steps278_327/`. All three run dirs moved to
  consolidated root + symlinks. 026 slate simultaneously ALL TERMINAL,
  redesign rescues nothing (see 026 provenance slate-verdict section) —
  032 omission is the standing lever. Ryan-gated offers: ctrl+omission
  discriminator, feature-A arm, baseline-for-ΔCT arm, notebook entry,
  archiving.
- 2026-09-19 (wrap): merge-burden analysis added to provenance (omission
  collapses fs-family merges ~30×, 84–87k→~3k → merge-law choice likely
  near-inert under omission; floor family never had a frenzy; np −30%
  relieves the count-driven FMM adequacy gate → NT ladders candidate).
  A/B CLOSED. Next session charter (026 next-runs recommendation +
  Ryan-gated launches + notebook entry w/ side-by-side GIF) =
  `032_reset_prompt_20260919b.md`.

- 2026-09-19 (later session): Follow-up arms HARVESTED — C-ctrl 13773490
  COMPLETED 467/467 CT 0.0714872 ± 0.067% (omission is a GENERAL rescue,
  ctrl channel included); C-oldlaw 13773491 COMPLETED CT 0.0716405 ±
  0.0014%, inside B15's noise band (§22 merge redesign unnecessary under
  omission — prediction held). Details:
  `032_followup_provenance_20260919.md` Outcomes. Notebook follow-up
  written (Ryan-approved). Ryan directive: NEW MERGE LAW = PRODUCTION
  DEFAULT henceforth (flowpanel→unified-052 merge in flight, commit
  d896145 pending test verdict). Re-open preps approved (022 hr10 /
  018 NT72 _3r / 020 Phase-2R discriminator; NT144 conditional on P2):
  `032_omission_reopen_prep_20260919.md`.
- 2026-09-19T16:52Z ARCHIVE (hpc-storage): 3 terminal root-omission arms
  (om15_explg_fs 16005→12346MB, om20_explg_fs 15773→12154MB,
  om15_explg_floor 12716→9704MB) tarred+verified to
  /nobackup/archive/usr/rander39/FLOWPanel_runs/projects_FLOWPanel.jl/,
  newest 5 restartable steps kept on /home; freed 43682MB. DEFERRED
  (RECENT-HOT <2h quiet): scr_p032om15_ctrllg_fs + explg_fs_oldlaw
  (~31G) — archive on a later pass. /home 476.9→396.4 GiB of 400 cap
  (3.6 GiB headroom); 212GB of p018_csarc_*_3r_* identified reclaimable
  but NOT authorized this cycle.
