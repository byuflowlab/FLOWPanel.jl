# 026 Phase 2 — implementation handoff (2026-09-07, split into two sessions)

Session 1 = FLOWVPM (plan commits 1–3). Session 2 = FLOWPanel (plan commits 4–6),
launched only after Session 1 reports green. Copy-paste the matching prompt below.

---

## Session 1 prompt — FLOWVPM (commits 1–3)

Implement **commits 1–3** of the BRAINSTORM 026 particle-splitting plan at
`~/.claude/plans/shimmering-wobbling-sunrise.md`. Read that file FIRST and follow it
exactly — it is context-complete (struct definitions, file:line targets for every edit,
commit sequencing with test gates, and Ryan's rulings dated 2026-09-05/07). Do not
re-plan or re-open closed decisions listed in its Context section.

Scope — FLOWVPM only (`~/Dropbox/research/projects/FLOWVPM.jl`, branch `flowpanel`):
1. Commit 1: `ResolutionSplitState` + ParticleField field + lockstep + integrator
   accumulation + direction resolution (plan Stage 1). Gate: 8.1 tests 6–7 + existing
   FLOWVPM suite green.
2. Commit 2: the three kernels `_split_viscous_tetra4!` / `_split_compress_tri3!` /
   `_split_elongate_pair2!` + random orientation draws (Stage 2), plus a short design-doc
   §3b amendment noting shrink/elongation events now route to the 2-child in-line kernel
   (third 2026-09-07 ruling). Gate: 8.1 tests 1–5, 9.
3. Commit 3: `split_particles!(pfield, ::ResolutionSplitOpts)` — trigger check, routing,
   single-loop main pass, counters — + FMM radix adequacy logging (Stage 3). Gate: 8.1
   test 8; `runtests_dsigma2_accumulators.jl` green.
Work each commit to its gate before starting the next. Do NOT touch FLOWPanel source
(commits 4–6 are a separate session), and do NOT start commits 7–8 (Ryan-gated).

Ground rules:
- The old splitting path must stay byte-identical: `src/FLOWVPM_splitting.jl` and
  `src/FLOWVPM_merging.jl` untouched; all new code goes in a NEW file
  `src/FLOWVPM_resolution_split.jl` (+ the minimal ParticleField/timeintegration edits
  the plan specifies). The existing `SplittingState` W2 accumulator writes stay as-is.
- Before your first commit, `git status` the repo; if pre-existing uncommitted 026 edits
  exist, commit them as their own housekeeping commit; leave unrelated files alone.
- New tests go in `test/runtests_resolution_split.jl`, included from runtests.jl (plan
  8.1). Max 4 local threads.
- The plan's open knobs D2/D2b/D3/D4 are env-settable defaults — implement the proposed
  defaults, don't block on them. D5 (SIGMA_CEIL) is commit 8: skip.

Deliverable: commits 1–3 landed with gates green, a compact per-commit summary
(what/where/test result), a list of any deviations from the plan (with justification),
and the exact FLOWVPM commit SHAs — Session 2 needs them. Append that summary to
`BRAINSTORM/026_sigma_growth_particle_splitting/phase2_impl_handoff_20260907.md`
(FLOWPanel repo) under "## Session 1 report". Offer — do not write — a notebook entry.

---

## Session 2 prompt — FLOWPanel (commits 4–6; launch AFTER Session 1 reports green)

Implement **commits 4–6** of the BRAINSTORM 026 particle-splitting plan at
`~/.claude/plans/shimmering-wobbling-sunrise.md`. Read that file FIRST, then the
"Session 1 report" section at the bottom of
`BRAINSTORM/026_sigma_growth_particle_splitting/phase2_impl_handoff_20260907.md` —
Session 1 already landed the FLOWVPM side (state struct, kernels,
`split_particles!(::ResolutionSplitOpts)`); honor any deviations it recorded. Do not
re-plan or re-open closed decisions.

Scope — FLOWPanel only (`~/Dropbox/research/projects/FLOWPanel.jl`, branch `fastmultipole`):
4. Commit 4: `ResolutionSplit` policy + merge `on_representative` reset closure +
   replay drop + GPU-seam comment (plan Stage 4). Gate: 8.2 policy/seam tests.
5. Commit 5: warm-start persistence — the all-or-nothing `rsplit_*` VTP block
   (Stage 5). Gate: 8.2 warm-start tests; IGE suite green.
6. Commit 6: driver env knobs + dispatcher cases (Stage 6) + FLOWVPM ring collective
   script `examples/p026_ring_split_test.jl` (Stage 8.3). Gate: ring clean;
   splitting-off smoke bit-identical.
Then run the plan's Stage 7 reconciliation checklist and report each line's status.
Do NOT start commits 7–8 (campaign launches, SIGMA_CEIL removal — Ryan-gated).

Ground rules:
- Read `FLOWPanel.jl/CLAUDE.md` and `agent_policies/WORKFLOW.md` + `TESTING.md` before
  editing; run tests via the test-runner subagent where the change maps to the
  TESTING.md matrix. Max 4 local threads.
- FLOWPanel has uncommitted state (026 design-doc/dispatcher edits, 018 handoff files) —
  before your first commit, commit the pre-existing 026 dispatcher/design-doc edits as
  their own housekeeping commit(s); leave 018/others' files alone.
- `src/FLOWVPM_merging.jl` stays untouched (the reset closure uses the existing
  `on_representative` hook from the FLOWPanel side).
- Known pre-existing local test breakage (NOT yours to fix): Julia 1.12 WeakKeyDict
  finalizer errors on testsets using immutable `WarmstartNoopSolver`.

Deliverable: commits 4–6 landed with gates green, per-commit summary, Stage-7 checklist
status, deviations list. Offer — do not write — a notebook entry covering the full
Phase-2 implementation (both sessions).

---

Provenance: plan approved by Ryan in-session 2026-09-07 after three revision rounds
(independent structs; no cooldown/skip_static/max_fraction/spacing-guard; physics
naming; two-regime stretch mechanism with new elongate-pair2 kernel). Supersedes
`phase2_handoff_prompt_1.md` / `phase2_handoff_prompt_2.md`.

## Session 1 report (commits 1–3, FLOWVPM, 2026-09-07)

All three FLOWVPM commits landed on `flowpanel` with gates green. **Session 2
needs these SHAs:**

| # | SHA | what | gate result |
|---|---|---|---|
| 1 | `9d63578` | `ResolutionSplitState{R}` + `ResolutionSplitOpts{R}` + `ParticleField.resolution_split` field (lazy, `nothing` ⇒ no-op) + add/remove lockstep + integrator-inline accumulation (euler, euler_exp `S=L*G0`, RK3 final-stage `b==8/15`) + dvisc/drvpm mirrors at the W2 sites + `_rsplit_direction` (coherence gate, Γ̂ fallback) | 8.1 t6–t7 (48 tests) + dsigma2 suite 44/44 + FULL FLOWVPM suite green |
| 2 | `99f4d54` | `_split_viscous_tetra4!` (σ_c=σ_p·4^(−1/3), a=1.3503σ_p, Shoemake SO(3) draw), `_split_compress_tri3!` (σ_c=σ_p/√3, ring 0.6σ_p, random in-plane angle), `_split_elongate_pair2!` (σ_c=σ_p, ±0.5σ_p in-line, no draw), shared `_rsplit_emit_child!` (Γ/m exact, vol=(4/3)πσ_c³ per W5, zeroed U/ω/J/PSE/M/C/SFS/U_prev, fresh state) | 8.1 t1–t5+t9 (788 tests cumulative). t4 regression pins on default offset ratios: L2 support residual tetra4 ∈ (0.50,0.68), tri3 ∈ (0.92,1.18), pair2 ∈ (0.095,0.112) — measured 0.57–0.61 / 1.03–1.07 / ≈0.103 |
| 3 | `5583443` | `split_particles!(pfield, ::ResolutionSplitOpts; verbose, dt)` — `_rsplit_check` (NaN-guarded; grow wins; floor pin with `sigma_0>floor` anti-refire guard), routing (shrink→pair2, grow→`dvisc≥drvpm ? tetra4 : tri3`, disabled ⇒ skip+count, never reroute), single serial pass, counters NamedTuple + verbose skip telemetry; `_radix_log_sigma_adequacy` at every radix coupling (re)build (`@info` ratio, `@warn` >0.9 with σ-cap guidance — observed firing live) | 8.1 t8 (812 tests total, all green); dsigma2 suite 44/44; radix FMM host suite 63/63 |

Old-path byte-identity audit: `git diff 9d63578^..5583443 -- src/FLOWVPM_splitting.jl
src/FLOWVPM_merging.jl` is EMPTY. New code is entirely in
`src/FLOWVPM_resolution_split.jl` + minimal particlefield/timeintegration/viscous/
fmm_radix edits. Tests in `test/runtests_resolution_split.jl`, included from
`runtests.jl` after the dsigma2 suite.

Deviations from the plan (all judged forced/minor):

1. **Include order**: `resolution_split` is included BEFORE `particlefield`, not
   after `splitting` — the `ParticleField` struct declares
   `resolution_split::Union{Nothing,ResolutionSplitState{R}}`, so the type must
   exist first (same reason `SplittingState` lives ahead of the struct).
   Consequence: functions in the new file leave `pfield` unannotated (dispatch is
   on the opts/state types); documented in the file header.
2. **No `Base.@kwdef`** for `ResolutionSplitOpts`: `@kwdef` on the parametric
   struct auto-generates a zero-arg constructor that collides with the
   `ResolutionSplitOpts() = ResolutionSplitOpts{FLOAT_TYPE}()` convenience method
   (precompile-time method-overwrite error). Hand-written kwarg constructor with
   the same defaults instead.
3. **`split_particles!` accepts `dt=nothing` (ignored)**: Stage 4's policy seam
   calls `FLOWVPM.split_particles!(pfield, opts; dt)`; accumulation is
   integrator-inline so dt is unused here — kept as a documented seam kwarg so
   Session 2 can follow the plan text verbatim.
4. **`_rsplit_emit_child!` also freshens the (inactive) legacy `SplittingState`
   slot** for the in-place child (appended children get this from `add_particle`),
   keeping lockstep bookkeeping coherent even though the two split policies must
   never be co-active. Zero behavior change for the old path (its files are
   untouched and its tests unmodified).
5. **dvisc/drvpm mirror sites are six, not "two"**: the landed W2 writes live at
   three viscous sites (`FLOWVPM_viscous.jl` euler/euler_exp/RK3 branches) and
   three rVPM sites (`FLOWVPM_timeintegration.jl` euler/euler_exp/RK3) — mirrored
   at all six, colocated as the plan intends.
6. **Design doc §3b amendment** added (blockquote under the §3b heading) noting the
   third 2026-09-07 ruling: shrink/elongation events route to the 2-child in-line
   kernel. Left UNCOMMITTED in FLOWPanel alongside the pre-existing §17–18
   working-tree edits (this session commits nothing in FLOWPanel per scope).

Env-knob defaults (D2/D2b/D3/D4) implemented as proposed: compress 0.6,
elongate 0.5, viscous 1.3503, `log_stretch_max` NaN-disabled until §9 arms.
D5 (SIGMA_CEIL removal) untouched (commit 8, Ryan-gated). RBF-reset site
untouched per ruling. Commits 4–6 (FLOWPanel policy/persistence/driver) and 7–8
remain for Session 2 / Ryan.

## Session 2 report (commits 4–6, FLOWPanel + FLOWVPM ring script, 2026-09-08)

All three commits landed with gates green, on `fastmultipole` (FLOWPanel) and
`flowpanel` (FLOWVPM). Housekeeping commit `9d2c0f1` first captured the
pre-existing 026 working-tree state (design-doc §17–18 + §3b amendment,
linegauss dispatcher default + rk3 arm, phase2 handoff/prompt files).

| # | SHA (repo) | what | gate result |
|---|---|---|---|
| 4 | `c78f6d1` (FLOWPanel) | `ResolutionSplit{TO}` policy (lazy `enable_resolution_split!` on every application so accumulation runs from step 1 at any cadence); mutual exclusion with `SplitParticles` enforced in the `ParticleMaintenance` tuple ctor; W3 merge interplay via an `on_representative` closure built in `apply_particle_policies!` at the MergeParticles application site (`FLOWVPM_merging.jl` untouched); `_heal_unseeded_rsplit_slots!` for device-shed particles; GPU-seam documentation in `_gpu_copy_side_buffers!` (NO new mirror entries — on device-backed wakes the state is canonical to the host mirror, device field stays `nothing`); replay serialize-then-drop | wake unit suite green incl. 17 new policy tests; replay suite 142+6 green; opt-in `FLOWPANEL_TEST_RESOLUTION_SPLIT_CUDA` seam testset (self-skips, no local GPU) |
| 5 | `ebbbf3b` (FLOWPanel) | rsplit_* VTP persistence: writer block (six fields, only when enabled, host-mirror-aware, series precision); all-or-nothing loader (all six → enable + exact restore; none → `nothing`, no warn; partial → `ArgumentError`); `_clear_splitting_state!` extended | new 18-test warm-start testset green; IGE suite 113/113 green (run in isolation — see known breakage below) |
| 6 | `e7681d2` (FLOWPanel) + `65247ee` (FLOWVPM) | WAKE_SPLIT_* env knobs per the Stage 6 table (enabled iff mechanism AND trigger; errors on mech-without-trigger, trigger-without-mech, SIGMA_CEIL+splitting, ON_FLOOR without SIGMA_FLOOR_FRAC>0); splice after MergeParticles; dispatcher arms `scr_p026sp_nt144_cap030`/`_cap018` + twelve `scr_p026s9_*_{floor,split,fs}` §9 matrix arms; self-contained ring collective test `FLOWVPM.jl/examples/p026_ring_split_test.jl` | all four error combos verified firing; enabled setup prints ACTIVE with (Merge, ResolutionSplit) order; splitting-off short-march smoke vs pre-commit-6 driver **bit-identical** (59 files; sole diff = wall_s wall-clock column of the wake-health CSV); ring script runs clean — control conserves to machine precision, split arms ≤0.8% circ / ≤0.4% impulse / ≤3% enstrophy drift over ~1 convective time |

Stage 7 reconciliation checklist:

- Merge hook resets ALL six fields, sigma_0 := merged σ — DONE (closure →
  `_rsplit_reset_slot!`; asserted by the W3 policy test).
- Merge→Split order documented — DONE (policy docstring, driver comment,
  splice order; tuple order preserved by `_split_particle_policies`).
- Every trigger self-limiting on fresh children — AUDITED (ratio restarts at
  1 via fresh sigma_0; grow caps: tetra4/tri3 birth σ_c < σ_p, pair2 σ_c=σ_p
  unreachable on the grow side because grow wins ties; exposure resets to 0;
  floor guard `sigma_0 > floor` disarms floor-born lineages — Session 1 t8
  anti-refire tests cover all of these).
- Floor never GATES splitting — AUDITED (grep: `sigma_floor` appears in
  `_rsplit_check` only as a shrink TRIGGER; the 052c guard clamps σ in the
  integrator and never consults split state).
- Skip telemetry returned + printed — DONE (counters NamedTuple + verbose
  skip println; WAKE_SPLIT_VERBOSE default true).
- Split-vs-merge churn observable — DONE (per-step `n_split_*` from the
  policy verbose line vs merge counts; the deferred-cooldown diagnostic).
- Old splitting path byte-identical / `FLOWVPM_merging.jl` untouched — HOLDS
  for this session: zero FLOWVPM `src/` edits (only `examples/` added).

Deviations / notes:

1. **FLOWVPM moved under us**: Ryan landed `119fe23` (merge runaway guard —
   touches `FLOWVPM_merging.jl`, hook signature unchanged) and `21eeaaa`
   (euler_exp sigma_guard) mid-session. All Session 2 work was built and
   tested against that newer HEAD; nothing in them conflicts with the split
   seams.
2. **GPU seam heal (small addition beyond plan text)**: particles shed on the
   device field between maintenance passes miss the add_particle lockstep
   hook on the host mirror (sigma_0 = 0 would read as infinite growth ratio).
   `_heal_unseeded_rsplit_slots!` seeds those slots from the current σ before
   each split application; host-backed wakes are a no-op scan. Documented
   device-path limitations (accepted): device integrator twins skip
   accumulation → axis/weight/exposure stay zero (Γ̂ fallback, exposure
   trigger inert) and dvisc/drvpm stay zero (grow ties route to the viscous
   kernel).
3. **§9 matrix thresholds are submission knobs**: SIGMA_FLOOR_FRAC (floor/fs
   arms) and WAKE_SPLIT_LOG_STRETCH_MAX (split arms) are deliberately NOT
   baked into the dispatcher — D4 is Ryan-open. The driver fails fast if an
   arm is submitted without its threshold (mechanism-without-trigger error),
   so no arm can silently run unarmed. Γ̂-comparison arms via
   WAKE_SPLIT_STRETCH_AXIS=false at submission.
4. **Replay drop**: ResolutionSplit is serialized (opaque opts, same
   convention as SplitParticles) and dropped without warning on
   deserialization — replay reads particle states from disk and never
   re-runs splitting.
5. **Known pre-existing breakage (NOT this session's)**: the first testset of
   `runtests_unit_warmstart.jl` errors on Julia 1.12 (WeakKeyDict finalizer
   on immutable `WarmstartNoopSolver`) and aborts the file; the 026 and IGE
   testsets were run in isolation (green). Also pre-existing dirty files left
   alone: 018 provenance/handoff docs, `scripts/p018_harvest_ct.py`,
   formulation src/test edits.
6. Ring-test output (`examples/p026_ring_split_test_out/`) is generated data,
   left untracked in FLOWVPM.

Commits 7–8 (campaign launches, SIGMA_CEIL removal) remain Ryan-gated.

## Session 3 report (ship resolution splitting as THE system; legacy removal, 2026-09-08)

Ryan authorization 2026-09-08 executed: resolution splitting is now THE
particle-splitting system of FLOWVPM; the legacy experimental path is
removed. Merging (`merge_particles!`) and the filament edge graph are
untouched features (per Ryan's mid-session reminder, only splitting was
removed).

FLOWVPM commits (branch `flowpanel`, on top of `65247ee`):

| # | SHA | what | gate result |
|---|---|---|---|
| 1 | `2a1b970` | first-class `run_vpm!` wiring: `split_every::Int=0` / `split_opts=nothing` kwargs (pattern-matched on `merge_every`/`merge_kwargs`), applied AFTER merging (W3), `enable_resolution_split!` called before step 1 so triggers integrate at any cadence, error on `split_every>0` without opts; D-A native merge reset in `_finalize_merged_particle!` (`rs === nothing || _rsplit_reset_slot!(rs, representative, sigma)`); docstrings updated (only splitting system) | `runtests_resolution_split.jl` 829/829 incl. new s3 wiring/merge-reset testset |
| 2 | `8b0b70d` | legacy removal: legacy machinery deleted (`SplittingState`/`SplitOptions`/`SplitDirection`/trigger tree/`_do_split!`/`accumulate_H_chi!` + exports + include); `ParticleField` drops `splitting_state`/`splitting_workspace`/`track_H_chi`/`H_chi_axis`/`H_chi_clip_positive` + lockstep hooks + `nextstep` hook; six legacy `dsigma2_*` writes dropped (mirrors kept) + legacy RBF-reset accumulator clear dropped (new mirrors untouched there per the standing Session-1 ruling); legacy merge-reset block replaced by D-A; `runtests_dsigma2_accumulators.jl` deleted after porting its invariants onto the `dvisc`/`drvpm` mirrors | `runtests_resolution_split.jl` 861/861; FULL suite `julia --project=test test/runtests.jl` exit 0, 0 failures (merging + filament suites green) |

**Inventory correction found during commit 2**: `src/FLOWVPM_splitting.jl`
was two subsystems in one file — legacy splitting (top ~650 lines) AND the
whole filament-edge-graph machinery (`add_edge!` … `filament_calibration_sweep`,
exported, 477+ tests). Deleting the file wholesale broke the filament suite;
resolved by moving the filament machinery verbatim (plus its four geometric
helpers `_strain_tensor`/`_eSe`/`_leading_eig_sym3`/`_unit_strength`/
`_unit_streamline`/`_filament_axis_unit`) to new `src/FLOWVPM_filament_edges.jl`
(included after `merging`), then deleting the file. Exports unchanged.

FLOWPanel companion commit (branch `fastmultipole`, on top of `ef239de`),
listed separately per scope:

| SHA | what | gate result |
|---|---|---|
| `7dd1dc7` | delete `SplitParticles` policy + apply + mutual-exclusion guard; `accumulate_H_chi!` call in `propagate!`; `split_*` VTP writer block + kwarg; `split_*` loader block + legacy `_clear_splitting_state!` lines; `splitting_state` GPU side-buffer entry; `SplitParticles` replay branch; W1 testset removed, new D-B testset added | wake unit suite green (16 policy tests); replay 142+6; warm-start 026 testset 23/23 isolated; IGE 113/113 isolated; smoke below |

**D-A outcome**: adopted as recommended — native guarded reset in
`_finalize_merged_particle!`; `on_representative` hook KEPT. The FLOWPanel
merge-hook closure (`_resolution_split_merge_hook`) is now
redundant-but-harmless and was LEFT IN PLACE (minimal-companion principle;
it double-resets the same slot idempotently).

**D-B outcome**: verified by a new permanent testset
("D: stale legacy split_* fields ignored silently") in
`runtests_unit_warmstart.jl`: a fabricated Session-2-era checkpoint carrying
all six `split_*` arrays plus a full `rsplit_*` block loads with no
warning/error, the stale arrays are inert extra point data, and the
`rsplit_*` state restores exactly. New checkpoints no longer write `split_*`.

**Grep audit**: `SplittingState|SplitOptions|H_chi|hold_counter|cooldown_counter|dsigma2_`
→ zero live references in FLOWVPM `src/`+`test/`+`examples/` and FLOWPanel
`src/`+`test/`. Sole intentional hits: three string literals inside the D-B
test that fabricate the stale legacy field names.

**Bit-identity smoke** (splitting off, Session-2 recipe: scratchpad cwd,
`NREVS=0.25 FREESTREAM_RAMP_REVS=0.1 FREESTREAM_HOLD_REVS=0.05
FREESTREAM_WITHDRAW_REVS=0.05 SETTLE_REVS=0.05`, `julia -t1`,
`examples/rotor_hover_pressure_comparison.jl`): baseline = git worktrees
FLOWPanel `ef239de` + FLOWVPM `65247ee`; candidate = post-removal trees.
108 output files, identical file lists, 99/108 md5-identical. The 9 diffs
are fully accounted for: wake-health CSV differs ONLY in the `wall_s`
column (allowed); case-metadata TOML only in absolute paths and wall-clock
timings; the 7 particle VTPs are bit-identical in every shared point-data
field (verified field-by-field via ReadVTK) with only the six removed
`split_*` arrays absent from the candidate — the intended format change.

Notes: the CoreSpreading RBF σ-reset no longer clears any Δσ² attribution
(the legacy clear was deleted with the legacy state; the new mirrors were
deliberately never cleared there per the Session-1 ruling — flag to Ryan if
CoreSpreading+splitting is ever co-armed). Commits 7–8 (campaign launches,
SIGMA_CEIL removal) remain Ryan-gated and untouched.
