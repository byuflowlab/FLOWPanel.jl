# Handoff prompt — rigid-motion tree/cache reuse (copy-paste for a fresh agent)

Prepared 2026-08-19 at context reset. Ryan launches a fresh agent with the
prompt below; the launch itself is the execution go the item's status line
asks for.

---

Implement the rigid-motion tree/cache reuse item for BRAINSTORM 021. The
complete staged plan is
BRAINSTORM/021_rotor_hover_solver_benchmarks/rigid_motion_tree_reuse_item.md
— read it first and follow its "Sketch (execution order)" exactly; it
contains the idea, the invariance arguments, the payoffs, and a CORRECTNESS
FLAG that must be verified FIRST. My launching you with this prompt is the
execution go — flip that file's status line to EXECUTING with today's date.

Required policy reads before touching src: agent_policies/WORKFLOW.md and
agent_policies/TESTING.md. HPC or anything >20 min local: read
agent_policies/HPC.md (you should not need the cluster for this task — see
constraints).

Repo state you inherit (do not "fix" it):
- FastMultipole dev checkout at ../FastMultipole, branch flowpanel-20260817,
  HEAD 1ec0af9. The five 2026-08-19 near-field-cache commits (72d3f3d cache
  core, 02f071c fmm!/FmmPlan integration, 0ef4e83 estimator + max_build_time
  guard, 204188a cached-path tuning, 1ec0af9 provided-cache leaf correction)
  are prerequisites this item builds on. scripts/20250404_prediction_accuracy.jl
  is dirty from unrelated work — never stage or commit it.
- FLOWPanel's working tree is heavily dirty and UNCOMMITTED by design:
  src edits (FLOWPanel_fmm.jl cache_nearfield plumbing, FLOWPanel_solver.jl
  KrylovSolver cache_tree/cache_nearfield, FLOWPanel_instrumentation.jl,
  FLOWPanel_metadata.jl), test additions, and benchmark scripts
  (rotor_hover_solver_phase2.jl nfcache configs,
  rotor_hover_solver_phase2b_nearfield_cache.jl). Do not revert, stash, or
  commit them; build on top. Commit only in FastMultipole, staging only
  files you touched.

Hard constraints:
- HPC p2b chains (nearfield-cache benchmarks, separate cluster checkout
  ~/projects/p2b/) and possibly Phase-1 chains are in flight: rsync NOTHING
  to the cluster, do not edit benchmark/rotor_hover_solver_phase1* or
  anything under benchmark/results/.
- Local runs ≤4 threads TOTAL; nothing here should need >20 min locally.

Execution notes (the item's sketch, condensed — follow the file's order):
1. FGS unsteady staleness discriminator FIRST (short R1 unsteady run, FGS
   vs krylov_gmres per-step strength divergence vs rotation angle; the
   unsteady driver pattern is benchmark/rotor_hover_solver_unsteady.jl with
   RHPC_SETUP_ONLY include — but do NOT edit that script; write your own
   small driver if needed). This establishes the stakes and is the item's
   pre-registered correctness check.
2. FastMultipole `transform_tree!(tree, R, t)` (rotate/translate branch
   centers; audit consumers of stored axis-aligned box dims before trusting
   adaptive-P paths) + unit test: fmm! on a rotated system with a
   transformed tree == fresh-tree fmm! rtol 1e-12, INCLUDING a panel system
   with vertices (vertices/positions refresh through buffers; verify).
3. `transform_plan!` + NearfieldInfluenceCache validity under rigid motion:
   Dirichlet scalar cache blocks are expected EXACTLY invariant — test
   cache-reused-after-rotation == fresh build (rtol 1e-12, shedding panels
   included); velocity/gradient outputs either get the per-block G→R·G
   handling or a loud v1 refusal (document which).
4. FGS: transform its trees via the SAME machinery (Ryan directive: shared
   parts); dense self/nonself matrices stay untouched (rotation-invariant);
   unsteady FGS-vs-Krylov A/B must collapse the divergence measured in
   step 1.
5. FLOWPanel surface: per-step transform hook from rigid-frame kinematics +
   the `persistent_plan` opt-in (the deferred nearfield-cache commit 4;
   its contract sketch is in nearfield_matrix_cache_plan.md §3 point 5 —
   kerneloffset state restored before every solve makes persistence sound).
6. Note (docs only): run-amortized build changes the *_nfcache knob
   economics — flag for the HPC re-examination, do not rerun benchmarks.

Testing conventions: every test needs non-vacuous premise guards (campaign
convention). Verification: FastMultipole full suite
(julia --project=test -t 4 test/runtests.jl from its repo root) and
FLOWPanel test/runtests_unit_solver.jl + runtests_unit_fmm.jl +
runtests_unit_simulate.jl at 4T all green. Gotchas: standalone FM test
runners must define the Gravitational compat overloads that live at the top
of test/solve_test.jl; Julia stack traces don't flush reliably to redirected
logs (silent log end = error, rerun foreground); don't pipe verdict-bearing
commands through grep/tail without checking pipestatus.

Log as you go: dated entry in BRAINSTORM/021_rotor_hover_solver_benchmarks/
phase_02_single_step_benchmarks.md Log + pointer in log.md (newest-first) +
update memory project_021_solver_benchmarks.md at the end. Judgment calls:
make them, log them, flag them for Ryan — don't block. Do not write to the
lab notebook; offer an entry instead.

Deliverable: staleness verdict with measured divergence numbers;
transform_tree!/transform_plan!/FGS-transform landed as FastMultipole
commits (suggested split: (1) transform_tree! + tests, (2) plan+cache
persistence + tests, (3) FGS + unsteady A/B fix); FLOWPanel hook left
uncommitted in the working tree like the rest; all suites green; docs +
memory updated. Stop before any HPC benchmark reruns — report instead.
