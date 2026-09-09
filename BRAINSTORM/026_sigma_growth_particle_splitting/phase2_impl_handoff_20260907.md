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
