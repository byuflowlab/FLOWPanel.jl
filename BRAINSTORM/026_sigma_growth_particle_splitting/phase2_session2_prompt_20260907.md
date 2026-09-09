# Session 2 prompt — BRAINSTORM 026 particle splitting, commits 4–6 (FLOWPanel + ring script)

(Prepared by Session 1, 2026-09-07. Session 1 landed FLOWVPM commits 1–3; its full
report is in `phase2_impl_handoff_20260907.md` § "Session 1 report".)

---

Implement commits 4–6 of the BRAINSTORM 026 particle-splitting plan at
`~/.claude/plans/shimmering-wobbling-sunrise.md`. Read that file FIRST and follow it
exactly — it is context-complete (struct definitions, file:line targets, commit
sequencing with test gates, Ryan's rulings dated 2026-09-05/07). Do not re-plan or
re-open closed decisions listed in its Context section. Then read
`BRAINSTORM/026_sigma_growth_particle_splitting/phase2_impl_handoff_20260907.md`
§ "Session 1 report" for what already landed and the API deviations you must match.

State you inherit (verified 2026-09-07):

- FLOWVPM (`~/Dropbox/research/projects/FLOWVPM.jl`, branch `flowpanel`, clean tree)
  has commits `9d63578` → `99f4d54` → `5583443`: ResolutionSplitState/Opts, the three
  kernels, and `split_particles!(pfield, opts::ResolutionSplitOpts; verbose=false,
  dt=nothing)` (dt accepted and ignored — accumulation is integrator-inline).
  `enable_resolution_split!(pfield)` is lazy; state lives HOST-side only. Merge
  interplay is a closure: pass `on_representative = i ->
  FLOWVPM._rsplit_reset_slot!(rs, i, sigma_merged)` at the FLOWPanel merge-policy
  site — FLOWVPM_merging.jl was NOT edited. Tests: `test/runtests_resolution_split.jl`
  (812 pass), included from runtests.jl.
- FLOWPanel (`~/Dropbox/research/projects/FLOWPanel.jl`, branch `fastmultipole`) has
  PRE-EXISTING uncommitted edits (design-doc §17–18 + §3b amendment, 018/021 docs,
  dispatcher, formulation src/test) that predate this session. Commit the
  026-related doc/handoff files as their own housekeeping commit before your first
  code commit; leave unrelated dirty files (018/021 docs, formulation src/test,
  dispatcher) alone unless commit 6 touches the dispatcher, in which case commit its
  pre-existing hunk separately from yours or ask Ryan.

Scope — commits 4–6 (plan Stages 4–6, gates from its commit table):

1. Commit 4 (FLOWPanel): `ResolutionSplit{TO}` policy in `src/FLOWPanel_wake.jl`
   next to `SplitParticles` (~:1644-1651); `apply_particle_policy!` calls
   `FLOWVPM.split_particles!(pfield, opts; ...)` on cadence with lazy
   `enable_resolution_split!`; Merge-before-ResolutionSplit ordering; the
   `on_representative` reset closure at the merge application site (~:1717-1728)
   when a ResolutionSplit policy is in the tuple; GPU seam comment in
   `FLOWPanel_gpu_wake.jl` (~:109-123 — host-side state needs NO new mirror
   entries; verify + comment); replay drop in `FLOWPanel_replay.jl` (~:735ff).
   Gate: 8.2 policy tests in `runtests_unit_wake.jl` (fires through maintenance,
   order, cadence, lazy enable) + CUDA-guarded seam testset (self-skips without
   GPU) + replay-drop test.
2. Commit 5 (FLOWPanel): persistence — VTP writer block in `FLOWPanel_wake.jl`
   (~:2436-2443): `rsplit_sigma_0`, `rsplit_axis` (3×np), `rsplit_weight`,
   `rsplit_exposure`, `rsplit_dvisc`, `rsplit_drvpm`, written only when
   `resolution_split !== nothing`; loader in `FLOWPanel_warmstart.jl` (~:328-360):
   ONE all-or-nothing field-presence block (all six → enable + exact restore; none
   → leave `nothing`, no warn; partial → typed error); extend
   `_clear_splitting_state!` (~:367-384). No TOML changes. Gate: 8.2 warm-start
   tests (`runtests_unit_warmstart.jl` ~:436-590 pattern) + IGE suite green.
3. Commit 6 (FLOWPanel driver/dispatcher + FLOWVPM ring script): env knobs in
   `examples/rotor_hover_pressure_comparison.jl` near :649 per the plan's Stage 6
   table (`WAKE_SPLIT_VISCOUS/STRETCH/SIGMA_MAX/SIGMA_GROWTH_RATIO_MAX/
   LOG_STRETCH_MAX/ON_FLOOR/STRETCH_AXIS/*_OFFSET_RATIO/EVERY/VERBOSE`; enabled iff
   mechanism AND trigger armed; error on inconsistent combos; error if SIGMA_CEIL
   set with splitting; splice `(GlobalCylinder, MergeParticles, maybe_split...)` at
   ~:728-735); dispatcher arms `scr_p026sp_nt144_cap030` (+ optional `_cap018`, §9
   s020v matrix) in `examples/run_p018_screen_hpc.slurm.sh`; ring collective test
   `FLOWVPM.jl/examples/p026_ring_split_test.jl` (§8.3: control vs force-split-all
   per kernel, CSV out). Gate: ring script runs clean; splitting-off smoke
   bit-identical.

Work each commit to its gate before starting the next. Do NOT start commits 7–8
(campaign launches and SIGMA_CEIL removal are Ryan-gated; D5 = commit 8, skip).

Ground rules:

- FLOWVPM's `src/FLOWVPM_splitting.jl` and `src/FLOWVPM_merging.jl` stay
  byte-identical. Prefer zero FLOWVPM src edits this session; if one proves
  unavoidable, it must be minimal, plan-consistent, and called out in the report.
- Plan file:line anchors may have drifted a few lines — re-grep before editing;
  do not trust them blindly.
- The two split policies (`SplitParticles`, `ResolutionSplit`) must never be
  active together — enforce/document at the policy layer per plan Stage 3.3 note.
- Max 4 local threads (`JULIA_NUM_THREADS=4`). FLOWPanel test env quirk: run suites
  the way `agent_policies/TESTING.md` says; for FLOWVPM use
  `julia --project=test test/runtests.jl` (`--project=.` trips TestEnv).
- Open knobs D2/D2b/D3/D4 are env-settable defaults — already implemented in
  FLOWVPM; don't block on them.

Deliverable: commits 4–6 landed with gates green; a compact per-commit summary
(what/where/test result), deviations with justification, and the exact FLOWPanel
commit SHAs, appended to `phase2_impl_handoff_20260907.md` under
"## Session 2 report". Offer — do not write — a notebook entry.
