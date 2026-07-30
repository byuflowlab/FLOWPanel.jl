# Progress Log — carrying out plans/20260721_panel_wake.md

Owner note: this log lets a fresh agent resume. Read `plans/20260721_panel_wake.md`
(the authoritative plan — its STATUS block mirrors this log) first, then this.
All local code + smoke deliverables are DONE; resume at the TODO checklist below
(HPC deploy onward).

**Update this log as you go:** whenever a TODO item completes (or partially
completes), move it to a DONE entry right away with concrete results — commands
run verbatim, Slurm job IDs, per-case wall times, output paths, and any
deviations from the plan with the reason. This log is the single source of truth
for execution state.

## Status snapshot (2026-07-21, updated after plan review)

### DONE — DirectWakePotential formulation (production)
Implemented in `src/FLOWPanel_formulation.jl`, exported in `src/FLOWPanel.jl`.
Productionizes the Task-3 diagnostic (`debug/dirichlet_solve/task3.md`):
`DirectWakePotential(; recompute_interval=1)` + `DirectWakePotentialState{TF}`;
`initialize_formulation` calls `_validate_formulation_common` then rejects
non-`PanelWake` wakes, `include_final_filament=true`, and any source with
`has_vector_potential`; **no DirectBackend guard** (runs on FMM, deliberate);
`solve_formulation!` does `clear_wake_correction!` → `_split_sigma!` → σ0 →
`_wake_potential!` → q_f → `_source_potential!` → Sσ0 → `rhs = −Sσ0 − q_f` →
`ldiv!(view(body.strength,:,2), solver.Glu, rhs)` → store σ0 in strength col 1.
**Key design fact:** no extra factorization needed — the `Backslash` solver's
`Glu` IS the finite-body G_Δ LU that Task-3 uses.

### DONE — tests in `test/formulation_test.jl`
Stages 8/9/10 appended; **ALL STAGES PASS (verified twice)**.
- Stage 8 Task-3 manual-solve equivalence 8/8; single-shot rel residual 8.6e-16.
- Stage 9 invalid-config rejections 8/8 (+ FMM backend accepted).
- Stage 10 ConstantDoublet vs VortexRing 5/5: q_f agrees 1.7e-13, γ 3.4e-14;
  induced-velocity differs 1.8e-4 (expected discrete-kernel difference; the
  formulation consumes only q_f).

### DONE — `examples/rotor_panel_wake_study.jl` (driver) + local smoke
Standalone driver, complete: winding-safe source+ConstantDoublet body (TE seeds
1614,1574,45 / 3324,3284,1755 ParaView 0-based → +1); finite
`PanelWake(nwakerows=WAKE_ROWS, core_size=1e-3, include_final_filament=false)`;
Backslash solver; FMM body order 8/acc 0.4/leaf 20, wake order 4/acc 0.4/leaf 50;
monitors PressureBernoulli(unsteady=true, allow_partial=false), pressure
ForceMonitor+RotorNormalization, SpanwiseLoadingMonitor (after force monitor),
KuttaJoukowskiForce, BoundCirculationMonitor, NO PressureLaplace; derived
`NREVS = ceil(WAKE_ROWS/NT)+FREESTREAM_DECREASE_REVS+SETTLE_REVS`; CT-vs-rev
CSV, metadata TOML, COMPLETED marker (only if all-finite), final-rev stats
(ptp≤5%, drift≤2.5%).
**Schedule note (deliberate, keep):** freestream ramp starts at t=0 then holds —
NOT fill-first — so the final buffer is entirely terminal-freestream-shed. See
the plan's "Schedule ordering" bullet for the rationale; do not change this.
**Smoke PASSED:** `data/dji9443_panelwake_smoke/` — COMPLETED marker,
`all_finite=true`, 12 finite steps, VTK + metadata + all 4 monitor CSVs;
CT_bernoulli_final_mean=0.00756, wall 204 s local (NT=12, 3 rows, 4 threads).
`stable=false` is expected at 1 rev. Do not re-run unless code changes.

### DONE — `examples/run_dji9443_panel_wake_convergence.jl` (adaptive driver)
Complete: subprocess-per-case; completion-marker + metadata checkpointing
(`is_complete` requires COMPLETED + `all_finite=true` + finite CT); axial ladder
72..432 (constant 4.0 m/s, settle 3); hover 0.25 ladder from
0.5/1.0/1.5× the axial cap (decrease_revs default 4, env `HOVER_DECREASE_REVS`),
extended to 432 **in steps equal to the ladder's initial spacing** (equal-sized
refinements — reviewed fix, do not revert to 36-row steps); zero-flow single
confirmation at the accepted hover cap (gentler ramp, `ZEROFLOW_DECREASE_REVS`
default hover+2); two-consecutive-refinement gate CT≤2% AND loading≤2% with
stability; time-budget aware (`BUDGET_S`, 15-min safety reserve, clean stop,
resume on resubmission); Stage-4 aggregate report with CCBlade references.
Parse-checked after the 2026-07-21 review edits (ladder step + dead-line cleanup).

### DONE — `scripts/dji9443_panel_wake_convergence.slurm.sh`
Matches the plan template: 1 node, `--ntasks=64`, `--mem=32G`, `--time=12:00:00`,
repo-dir guard, thread env exports, no cd/srun. `--time` is a first estimate;
re-tune from the first axial rung's wall time (checkpoint-resume makes
underestimation safe).

### IN PROGRESS — HPC deploy (2026-07-21)
- (a/b) All six study files were **already on the cluster and byte-identical**
  (sha256 match local↔remote for `src/FLOWPanel.jl`, `src/FLOWPanel_formulation.jl`,
  `test/formulation_test.jl`, both `examples/*.jl` drivers, and the launcher);
  remote git status shows only these expected modifications/untracked files.
  Remote checkout: branch `fastmultipole` @ `b094d02` (matches local base).
- (c) Assets verified on cluster: `examples/data/dji9443_new_40_40.msh` present
  (311269 bytes); `data/rotor_axial_j0187_ccblade/` complete (sectional ncrit4/9
  CSVs, polar CSV, validation CSV all present). No prior `data/dji9443_panelwake/`
  or job logs; queue empty.
- **DEVIATION (justified): julia is NOT on PATH in non-interactive ssh shells**,
  and the shared `/apps/juliaup` launcher errors (permission denied on its config
  db). The Manifest pins julia 1.11.7; the site spack module
  `julia/1.11.7-6bmogfl` provides a working binary at
  `/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin/julia`.
  Added a guarded PATH fallback (`command -v julia || export PATH=<that bin>:$PATH`)
  to `scripts/dji9443_panel_wake_convergence.slurm.sh` (local+remote in sync,
  sha256 `83a1797e...`), since the batch job inherits the non-interactive submit env.
- (d) On cluster: `bash -n` launcher OK; `Meta.parseall` parse checks OK for the
  two drivers + formulation file; focused `test/formulation_test.jl` **ALL 10
  STAGES PASS on the cluster** (julia 1.11.7, 2 threads; Stage 8 single-shot rel
  residual 5.03e-16; Stage 9 rejections 8/8; Stage 10 q_f 1.85e-13, γ 4.17e-14;
  full run ~5.5 min, dominated by Stage 7).

### DONE — HPC deploy (2026-07-21). Submission command (recorded verbatim, plan
"Execution Boundary" session-scoped authorization):
```
ssh orc 'cd /home/rander39/projects/FLOWPanel.jl && sbatch scripts/dji9443_panel_wake_convergence.slurm.sh'
```
Reuse this identical command for any resubmission after allocation timeout.
**Note:** non-interactive ssh lacks `sbatch` on PATH too — the working invocation
(used for the actual submission) wraps it in a login shell:
```
ssh orc "bash -lc 'cd /home/rander39/projects/FLOWPanel.jl && sbatch scripts/dji9443_panel_wake_convergence.slurm.sh'"
```

### SUBMITTED — Slurm job **12860375** (2026-07-21). Logs:
`/home/rander39/projects/FLOWPanel.jl/slurm-fp-dji9443-panelwake-12860375.{out,err}`.

### FAILED — job 12860375 (2.6 min): nested-RUN_NAME path bug; FIXED locally
First case `axial_rows072` crashed at the end-of-run metadata write:
`SystemError: opening file "data/dji9443_panelwake/axial_rows072/dji9443_panelwake/axial_rows072.metadata.toml"`.
Root cause: the adaptive driver passes `RUN_NAME=dji9443_panelwake/axial_rows072`
(slash nests the save dir, per plan), but `examples/rotor_panel_wake_study.jl`
used the full `run_name` as the `simulate!` `name` and as its CSV/metadata file
prefix → paths nest a second time inside `save_path`. simulate!'s metadata writer
doesn't mkpath, so it crashed; monitor CSVs also landed in a nested
`monitors/dji9443_panelwake/` subdir where the adaptive driver's non-recursive
`readdir` spanwise reader would never find them. (The original smoke passed
because its RUN_NAME had no slash. The adaptive driver already expected
`basename(run_name)` filenames — the study driver was the one at fault.)
**Fix (study driver only):** added `base_name = basename(run_name)` and use it
for `simulate!(name=...)` and the `_CT_vs_rev.csv` / `_study_metadata.toml`
prefixes. Verified by a local nested-RUN_NAME smoke
(`RUN_NAME=dji9443_panelwake_smokenest/axial_rows003`, NT=12, 3 rows).
Physics/schedule untouched — the time march itself ran fine for all 180 steps.
Before resubmission: rm the incomplete remote `data/dji9443_panelwake/` attempt.
Nested smoke PASSED locally (COMPLETED marker, flat filenames, flat
`monitors/axial_rows003_*.csv`, CT_bern mean 0.007556 — matches the original
smoke 0.00756; wall 207 s). Fixed driver redeployed (sha256 `79a99eee...`,
local=remote), incomplete remote `data/dji9443_panelwake/` removed.

### RESUBMITTED — Slurm job **12860413** (2026-07-21), same recorded command.
Logs: `slurm-fp-dji9443-panelwake-12860413.{out,err}` at the repo top level.

### DIAGNOSIS (2026-07-21) — axial rows072 CT way low; localized to the SOLVE
First rung (axial_rows072): CT_bern 0.00806 (ptp 6.6e-5, stable), CT_kj 0.00317 —
vs CCBlade 0.0325–0.0420. Subagent compared final-step body strengths vs the
trusted particle-wake run at the same condition (orc:data/rotor_axial_j0187_ccblade,
FINER mesh 28,432 cells vs 7,288 — caveat, but prior steady study showed mesh
moves CT only ~0.25%):
- sigma: new run equals -(Uinf-Omega×r)·n to machine precision — source BC fine.
- doublet jump (bound circulation): new/ref = 0.44 integrated, SAME sign, but
  deficit grows outboard — 0.66 @ r/R=0.29 → 0.26 @ tip; loading peak pushed
  inboard (0.51 vs 0.64). Not a 4π/sign artifact (ratio varies with radius).
- Cross-check: CT_kj(new)/CT_bern(ref) = 0.43 ≈ circulation ratio 0.44 (KJ linear
  in Γ) — internally consistent: the solved μ is ~2.3× low overall, ~4× at tip.
- SEPARATE flag: new-run CT_bern (0.00806) is 8% ABOVE ref despite 2.3×-low μ →
  Bernoulli pressure recovery inconsistent with the solved μ in this formulation.
Artifacts in scratchpad: strength_comparison_summary.md, spanwise_doublet_jump.csv,
spanwise_binned_strengths.csv (+ scripts, copied .vtu files).
**MAJOR CORRECTION (steady-Bernoulli check, same day):** the ref CT 0.00743 used
in the cross-check above is the *pre-aa38710 BROKEN* Bernoulli — the ref run's
trusted axial CT is **≈0.050** (CT_laplace_matderiv 0.0512, CT_laplace_lamb
0.0506, corrected bernoulli_replay 0.0503; see
`data/rotor_axial_j0187_ccblade/bernoulli_replay/bernoulli_forensic.csv`, ratio
fixed/run ≈6.8). Steady-Bernoulli recomputed from the stored VTU velocity gives:
NEW 0.0028 (≈ its KJ 0.0032), REF 0.0075 (= its broken in-run Bernoulli, NOT
0.050) — i.e. the naive velocity-based Bernoulli path underestimates ~7× on
rotating axial cases in BOTH runs; both VTUs' `velocity` fields equally violate
rotating-frame impermeability (|n·u_rel| rms ratio 0.75/0.78; 16 root cells in
NEW have |u| up to 11,600 m/s). Consistency within the velocity-based family:
new/ref steady-Bern 0.37 ≈ tip-weighted μ ratio — everything scales with the
2.3×-low μ. Bottom line: (1) μ deficit (0.44) is real and is the solve-side bug
signal; (2) NONE of the new run's in-run CT numbers (Bern 0.00806, KJ 0.00317)
are trustworthy absolute thrust — the trusted path on this branch is
PressureLaplace / corrected-Bernoulli replay, which the study driver deliberately
omitted (plan said no PressureLaplace). Script: scratchpad/steady_bernoulli_ct.py.
### AXIAL LADDER RESULT (job 12860413) — CONVERGED FLAT AT THE LOW ANSWER
rows=72 CT 0.0080554 (45.4 min) → 144 CT 0.0079971 (82.4 min, dCT 0.73%,
dLoad 0.85%) → 216 CT 0.0079785 (128.2 min, dCT 0.23%, dLoad 0.26%);
ladder formally CONVERGED at cap 216. (CT = the in-run Bernoulli path; flatness,
not the absolute value, is the signal.) **Truncation hypothesis is DEAD** — the
2.3×-low μ is wake-length-converged ⇒ defect in how the marched wake's q_f
enters the DirectWakePotential solve (prime suspects: production shed-strength
transfer to the wake rows — the Stage-8/10 tests used prescribed frozen
strengths — or the missing TE/final-filament potential-jump representation).
All three axial checkpoints (COMPLETED markers) preserved on orc.

**Job 12860413 then FAILED on its own at 04:16:07** immediately after "LADDER
CONVERGED": `cannot assign a value to imported variable Base.step` at
`run_dji9443_panel_wake_convergence.jl:249` (top-level variable `step` collides
with `Base.step` after use as a function). My scancel raced its death; hover
stages never started (would have been wasted anyway pending the formulation
fix). **FIXED:** renamed to `rung_step`, parse-checked, redeployed to orc.
DO NOT resubmit until the q_f/shedding defect is diagnosed and fixed.

### SHED-STRENGTH AUDIT (subagent, 2026-07-21) — formulation machinery EXONERATED
Verified numerically on a 12-step marched small wing: (1) TE-jump→wake-row
transfer exact (ratio 1.0000, right sign, all 7 edges; `shed_wake!`
src/FLOWPanel_simulate.jl:957-992 → `_get_wakestrength_mu`
src/FLOWPanel_liftingbody.jl:1603-1608; one benign convection-step lag);
(2) q_f covers all rows, Direct-vs-FMM wake backend 2.9e-16, FMM march
reproduces Direct γ 1.0000; (3) production RHS matches Task-3 term-for-term, no
σ₀/q_f double-removal; (4) **bug does NOT reproduce on a marched wing** —
DirectWakePotential vs VelocityThroughSources μ differ only ~5-7% (direct
slightly HIGHER). Benign caveat: evaluating wake influence between shed_wake!
and update_TE! hits degenerate row-1 nodes → NaN (production order never does).
Artifacts: scratchpad/shed_audit/{shed_audit.jl,qf_parity.jl,results.txt}.
Ranked rotor-specific hypotheses:
1. **Confounded comparison** — reference is a PARTICLE wake; deficit may belong
   to the finite PanelWake wake model, not the formulation.
2. **Wrong-side q_f branch** — doublet-sheet potential jumps by γ across the
   sheet; preceding blade's sheet passing close under the next blade could give
   O(γ) wrong-branch q_f at nearby control points (deterministic, converged,
   tip-weighted — matches 4× tip deficit). Untested regime (all validations had
   downstream-trailing wakes).

### DIAGNOSTIC SUBMITTED — Slurm job **12864823** (fp-dji9443-veldiag)
`scripts/dji9443_velocity_diag.slurm.sh` (new, deployed): identical axial
rows=72 case but `FORMULATION=velocity` + same finite PanelWake →
`data/dji9443_panelwake_diag/axial_velocity_rows072/`. Discriminates hypothesis
1 (velocity run ALSO ~2.3× low γ ⇒ wake model) vs 2 (recovers reference γ ⇒
potential-pathway defect). ~2.5 h allocation (rows072 direct took 45 min; the
velocity formulation adds wake-induced-velocity→σ work each step).
**RESULT (job 12864823, COMPLETED 42 min): HYPOTHESIS 1 CONFIRMED — the
formulation is VINDICATED; the deficit belongs to the finite PanelWake model.**
Velocity-formulation run: CT_bern 0.00800, CT_kj 0.00536, stable=true.
Spanwise doublet-jump comparison (scratchpad, direct vs velocity vs particle
ref): **direct and velocity μ agree to 3-4 significant digits at every radial
station** (jump ratios to ref identical: 0.77@r/R 0.25 → 0.45@0.55 → 0.25@0.95),
i.e. with the SAME finite panel wake, VelocityThroughSources solves essentially
the SAME 2.3×-low circulation as DirectWakePotential. The wrong-side-q_f
hypothesis is dead too (velocity path has no q_f). Note: KJ CT differs between
formulations (0.00317 vs 0.00536) despite identical μ — KJ samples local
velocity, which differs by pathway; another reason not to trust in-run KJ/Bern
absolutes here. **New lever: why does the finite PanelWake induce so much more
downwash (or mis-convect) vs the particle wake at the same condition?**
Candidate angles for next session: wake geometry comparison (panel-wake node
positions vs particle positions in VTK), core_size=1e-3 smoothing difference,
paired-shedding filament structure the particle wake has that rows lack, and
the semi-infinite-baseline steady solve (CT 0.0505) as the third reference.

**Open hypothesis discriminator (job left running deliberately — OVERTAKEN,
see above; kept for context):** in the
potential formulation a truncated doublet sheet lacks the far-wake potential jump,
which *lowers* μ (opposite of velocity-formulation truncation). 3-row smoke CT
0.00756 → 72-row 0.00806 is a slow rise. If rows 144/216 climb materially, the
deficit is (slow) wake-length convergence, not a code bug; if flat, it's a bug in
DirectWakePotential's marched-wake q_f (e.g. shed-strength transfer / missing
final-filament compensation). Axial ladder is self-limiting: if CT keeps moving
>2%, no cap is accepted and hover stages are skipped.

## TODO (remaining work — resume here)
1. ~~**HPC deploy**~~ DONE (see above) (plan "HPC Launcher and Deployment"; read `agent_policies/HPC.md`
   first). Cluster checkout `/home/rander39/projects/FLOWPanel.jl` (host `orc`).
   (a) Compare local vs remote for every file to be copied; abort on unrelated
   remote modifications. (b) Copy ONLY the study files: `src/FLOWPanel.jl`,
   `src/FLOWPanel_formulation.jl`, `test/formulation_test.jl`,
   `examples/rotor_panel_wake_study.jl`,
   `examples/run_dji9443_panel_wake_convergence.jl`,
   `scripts/dji9443_panel_wake_convergence.slurm.sh` — never a deleting/whole-tree
   sync; preserve unrelated remote state. (c) Verify on cluster:
   `examples/data/dji9443_new_40_40.msh` and `data/rotor_axial_j0187_ccblade/`
   CSVs exist. (d) From the cluster checkout: Julia parse/include checks,
   `julia --project=. test/formulation_test.jl`,
   `bash -n scripts/dji9443_panel_wake_convergence.slurm.sh`.
2. **Submit** (session-scoped sbatch authorization is recorded in the plan's
   "Execution Boundary" section). Record, then run, exactly:
   `ssh orc 'cd /home/rander39/projects/FLOWPanel.jl && sbatch scripts/dji9443_panel_wake_convergence.slurm.sh'`
3. **Monitor / resubmit.** Check `slurm-fp-dji9443-panelwake-<jobid>.{out,err}`
   at the repo top level; on allocation timeout, resubmit the identical command —
   the adaptive driver resumes at the first incomplete case via completion
   markers under `data/dji9443_panelwake/`. After the first axial rung reports
   wall time, re-estimate and (if warranted) adjust `--time` for resubmissions.
4. **Reports.** Axial: CT histories + final-cycle stats, direct-potential vs
   Kutta–Joukowski CT, spanwise loading + integrated CT vs CCBlade ncrit=4
   (CT 0.0419673) and ncrit=9 (CT 0.0325479), wake-row convergence table/plots.
   Hover: same + comparison (not equation) with experiment ≈0.072, viscous BEM
   ≈0.068, particle wake ≈0.062, steady rigid wake ≈0.0505; report 0-flow if it
   holds stability, else 0.25 m/s labeled J≈0.0117 with the destabilization
   freestream documented; record the decrease-rate used for each accepted result.
5. **Final reviewed conclusion** separating: solver/formulation correctness;
   wake-row convergence; mesh independence (cited from steady study); stability
   under terminal freestream; physical agreement. Apply the CLAUDE.md
   conclude→review→report loop.

## Useful references already located
- Task-3 solve mechanics: `debug/dirichlet_solve/dirichlet_solve.jl:1100-1151`
  (`_task3_direct_solve`).
- Formulation wiring in simulate: `_scalar_potential_sources` /
  `_collect_wake_scalar_sources` at `src/FLOWPanel_simulate.jl:269-279`;
  `formulation_prewake!` call at :363; `solve_formulation!` call at :454.
- Winding-safe RigidWakeBody: `examples/rotor_hover_convergence.jl`; CLAUDE.md
  invariant.
- Freestream schedule implementation: `examples/rotor_panel_wake_study.jl`
  (`magVinf`, ramp from t=0; NREVS derivation near the top of the env-config
  block).
