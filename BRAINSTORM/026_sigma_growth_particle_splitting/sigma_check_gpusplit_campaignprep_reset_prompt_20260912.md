# 026 reset prompt — σ-distribution check → GPU splitting → campaign prep (2026-09-12 evening)

You are picking up BRAINSTORM 026 (resolution-preserving particle splitting)
in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (+ sibling
FLOWVPM.jl, FastMultipole). Read `CLAUDE.md` and the policies it names.
Ryan directed this exact three-task sequence on 2026-09-12.

**Where the item stands** (full detail:
`hpc_smokes_harvest_reset_prompt_20260912.md`, same directory, FINAL
section; design context: `adaptive_elongation_reset_prompt_20260911.md`):
all three verification smokes PASSED (A adaptive default; B
MERGE_OVERLAP=3.5, low retention ruled a PASS as σ-scaled merging; C legacy
pair2, exactly 2 children/event). Commits DONE: FLOWVPM `f51f4ee` on
`flowpanel` (splitting impl + tests), FLOWPanel side in WIP snapshot
`004ce84` on `fastmultipole`. The HPC silo used for B/C is deleted. Working
trees are clean of 026 obligations except the perpetually-rewritten
`data/rotor_hover_pressure_comparison/rotor_hover_pressure_comparison.metadata.toml`
(leave uncommitted).

## Task 1 — σ-distribution check on smoke B (small, do first)

Question: is B's 2.7× lower particle retention (8879 vs A's 24275 @ step
233) healthy σ-scaled thinning of the aged wake, or split-merge churn
(children re-merging)? This ruling shapes the campaign's per-arm
MERGE_OVERLAP choice, so it precedes Task 3's knob proposals.

Mechanism to test: with `sigma_relative=true`, `merge_particles!` uses
`r_pair = σ_min/Φ_merge` (FLOWVPM_merging.jl:562, Φ_merge=3.5); legacy uses
absolute 0.02R (R=0.119 m). Crossover σ* = 3.5·0.02R ≈ 0.0083 m ≈ 2.2× the
shed σ (~0.0381R ≈ 0.00453 m). σ-thinning hypothesis ⇒ merging concentrated
in aged/large-σ particles (σ > σ*); churn hypothesis ⇒ merging of
just-split children near shed σ.

Data on disk (local, all smoke-B written — B ran full 467 steps and
overwrote A's VTK in the shared dir):
- B per-step VTK: `data/rotor_hover_pressure_comparison/rotor_hover_pressure_comparison_wake1_particles/*.{step}.vtp`
  (f32, appended binary; σ is a point-data array — read with ReadVTK.jl or
  a small Julia script, not grep).
- A per-step particle counts: session-A scratchpad
  `/private/tmp/claude-502/-Users-ryan-Dropbox-research-projects-FLOWPanel-jl/74eae7d8-c925-45f8-83ab-41a854dafb9e/scratchpad/smoke_a_particle_counts.csv`
  (+ `smoke_b_mergeoverlap.log`, `hpc_smoke_b_particle_counts.csv`,
  `hpc_smoke_b_13662583.log` there; copy anything you need out of the
  scratchpad early — it is session-scoped).
  NOTE: A's VTK is gone (overwritten); only its counts CSV survives.

Deliverable: σ histograms/quantiles at ~step 233 and 467, fraction of
particles above σ*, and (if feasible from B's VTK alone) an age/σ profile of
where the population deficit vs A sits. If churn can't be excluded from
snapshots, say so and propose the cheapest instrumented probe (e.g. a
~50-step rerun logging merges of particles younger than N steps) — but run
it only within the ≤4-thread budget. Report the verdict to Ryan with the
evidence before Task 3's MERGE_OVERLAP proposals.

## Task 2 — GPU splitting implementation

Entry point: `gpu_split_merge_reset_prompt_20260911.md` (this directory) —
read it FIRST and follow it; per its status the remaining gap is
**device-side accumulators only** (host-mirror path already works; see the
device-path comment in `src/FLOWPanel_gpu_wake.jl`). Notes since it was
written:
- The FLOWVPM splitting base it builds on is now COMMITTED (`f51f4ee`), so
  work from a clean tree; the CPU tests are
  `test/runtests_resolution_split.jl` (t10 adaptive elongation, t11 circ).
- Any new commits need Ryan's approval (the 09-12 commit ruling covered only
  what is already committed).
- GPU smoke/validation runs: HPC per `agent_policies/HPC.md` (GPU pools,
  arch-keyed launcher pattern, worktree/no-silo rules — the 09-12 silo was a
  one-off Ryan exception, already deleted; do not create silos without his
  say-so).

## Task 3 — campaign prep (commit-7 arms; NO launches without Ryan)

Prepare, do not launch:
1. Re-derive the split fractions in fraction space (old values predate §19
   fractional gating; design doc `particle_splitting_design.md` §19/§20).
2. Draft the arm matrix: §8.4 cap030/cap018 arms + §9 s020v matrix, with
   per-arm `WAKE_SPLIT_ELONGATE_OVERLAP` / `MERGE_OVERLAP` proposals. First
   discriminator = MERGE_OVERLAP=3.5 vs absolute-radius merging A/B
   (σ-pump, theory doc §4) — informed by Task 1's verdict.
3. Stage reproducibility mechanics per global policy: annotated tags
   (`campaign/p026-<slug>-YYYYMMDD`) + worktrees for the
   FLOWPanel/FLOWVPM/FastMultipole triple, Manifest dev-paths at the
   worktrees, pins recorded in a provenance file. FLOWVPM is committed
   (`f51f4ee`); FLOWPanel `fastmultipole` tip is `004ce84`; check
   FastMultipole for uncommitted state before tagging (if dirty, ask Ryan
   what to include).
4. Present the whole package (fractions, arms, knobs, pins, cost estimate)
   to Ryan for gating. GPU backends preferred for particle-heavy arms where
   sensible (Ryan 2026-09-05), which is why Task 2 precedes this.

## Ground rules (carry-over)
Local runs ≤4 threads TOTAL (check the machine first). Harness kill bug:
long local jobs via `nohup ... & disown` + log polling, never harness
run_in_background. `ssh orc` needs the live ControlMaster socket (2FA
otherwise); slurm needs `PATH=$PATH:/apps/slurm/latest/bin` in non-login
shells. Commits and campaign launches Ryan-gated. Design doc append-only;
theory doc rewritten in place. Notebook writes need Ryan's approval. Known
unrelated failure: `runtests_unit_warmstart.jl` first testset — not ours.
