# Phase 1 — Consistency: same solution, matched tolerances, no wake

**Status:** IN PROGRESS 2026-08-13 — tentative ladder selected (see Log); meshes to
generate, then cost-check and freeze. **Approval:** [ ]

## Prior-phase handoff

Phase 0 delivers: all six configs runnable, shared true-residual evaluator, history capture,
CSV schema, smoke harness.

## Scope

No-wake solves (single `steady!`; skeleton = `examples/rotor_hover_convergence.jl`, minding the
shedding-from-constructed-cells invariant in CLAUDE.md) on the mesh ladder. Two deliverables:

1. **Ladder selection (then frozen)**: pick ≥4 rungs from `examples/data/dji9443_*` spanning
   the available panel counts (candidates: 30_81 → 30_97 → 30_121 → 45_145 → 45_185); record
   exact files + panel counts in `ledger.md`. All later phases use exactly these rungs.
2. **Agreement + calibration**:
   - Reference = `Backslash` (direct; residual at round-off scale).
   - **FMM residual floor first (ruling 11c)**: per rung, evaluate the Backslash solution's
     true residual through the FMM operator (`apply_G!` with `FastMultipoleBackend`) — that
     mismatch is the floor any matrix-free config can reach in *direct-evaluated* residual.
     Set the matched-residual targets above the measured floor; if the floor sits above the
     provisional 1e-6-relative target, tune `expansion_order` up (recording the cost) before
     freezing targets.
   - **Right-preconditioning only (ruling 11a)**: all preconditioned Krylov configs run with
     `N=` routing so the stopping criterion is the true residual. (The Jacobi config's
     Phase 0 left-routing stopped 20× worse in true residual than plain GMRES at the same
     nominal tolerance — that comparison was not honest.)
   - For each iterative config, find settings (tol, `itmax`, `preconditioner_cell_size`, FGS
     `max_iterations`/`inner_iterations`/`rlx`/`expansion_order`/`leaf_size`, FMM knobs) that
     achieve a **matched true-residual level** (targets in `decision_rules.md`) on every rung.
     Smoke-scale settings (`rtol=1e-4`, `itmax=400`) are runtime conveniences, never carried
     forward (ruling 11b).
   - Agreement metrics vs reference: relative L2 and max norm of the strength-vector
     difference, plus CT from the three pressure monitors as physical cross-checks.
   - **A mismatch at matched residual is a bug**: diagnose and fix before proceeding (root
     cause first, per policy). Known suspect class: sign/convention drift between backends.

## Exit criteria

1. Agreement table (solver × rung × threading mode): strength-difference norms + CT deltas,
   all within `decision_rules.md` thresholds.
2. Frozen per-solver settings table — copied into `decision_rules.md` and used verbatim by
   Phases 2–4.
3. Any fixes landed with tests.

## Log

(dated entries appended newest-last)

### 2026-08-13 — Tentative ladder selected (constant aspect ratio, ×2/rung, cost-ceiling drop-out)

Phase 1 opened per Ryan's ruling (Phase 0 clear-context approved same day). Ryan's design
directives for the ladder, replacing the doc's original 30_81…45_185 candidate list:

1. **Refine spanwise and chordwise together** so panel aspect ratio stays constant (the
   original candidates drift chord/span panel ratio 2.76→4.18 across rungs);
2. **Chordwise parity chosen so control points sit adjacent across the x-y plane** —
   `scripts/generate_dji9443_mesh.sh` already enforces odd `Tess_W`, which gives an even
   chordwise panel count with equal upper/lower tessellation, i.e. control points mirror
   pairwise across the section plane;
3. **Span N far enough that matrix-constructing methods become genuinely expensive**, then
   drop them (per config, not per rung) once cost exceeds a ceiling — demonstrating that
   matrix-free configs reach problem sizes matrix-ful ones cannot. Tentative ladder first,
   adjusted from measured cost.

Mesh facts (verified): cells N = 4(U−1)(W−1) for uncapped two-blade meshes (45_145 →
25,344 ✓, 30_81 → 9,280 ✓); generation is local and cheap via
`scripts/generate_dji9443_mesh.sh U W uncapped`.

**Tentative ladder** (in `ledger.md`): anchor the aspect ratio on the workhorse 45_145
((U−1,W−1) = (44,144)); scale both by √2 per rung → N doubles per rung. R1 23_73 (6,336)
· R2 31_103 (12,648) · R3 45_145 (25,344, exists) · R4 63_205 (50,592) · R5 89_289
(101,376) · R6 125_409 (202,368) · R7 177_577 (405,504; added per Ryan 2026-08-13).
Dense-G memory 8N²: 0.3 / 1.3 / 5.1 / 20 / 82 / 328 GB / 1.3 TB — Backslash drops out
mid-ladder by construction; the exact ceiling is set from measured cost once R1–R3
numbers exist. All W odd ✓.

Ryan rulings recorded same day: Phase 1 development runs local (≤6 threads, small rungs),
the production sweep runs on HPC; **Phase 1's frozen settings + agreement tables are
treated as published numbers (ruling 5: HPC-only)**. Next: generate R1/R2/R4–R6 meshes,
verify N and watertightness, cost-check, then freeze the ladder in `ledger.md`.

### 2026-08-13 — Meshes generated; OpenVSP chordwise rounding found; R2/R4 revised

All missing rungs generated locally (`scripts/generate_dji9443_mesh.sh U W uncapped`,
seconds each) and verified with a throwaway script (cell count, TE detection incl.
`iswatertight` topology check, shedding traced on the *constructed* body per the
CLAUDE.md invariant; every rung sheds the full span, U−1 edges per blade).

**Finding: OpenVSP rounds the chordwise panel count up to a multiple of 8.** R1/R3/
R5–R7 already had (W−1) ≡ 0 (mod 8) and came out exact; the tentative R2 31_103 and
R4 63_205 did not — they produced 12,480 and 51,584 cells (chordwise 102→104,
204→208), i.e. mislabeled files. Every previously existing DJI mesh also satisfies
the mod-8 rule, corroborating. **Revision**: require (W−1) ≡ 0 (mod 8) and re-balance
U to hold the anchor AR ≈ 3.27 → R2 = 33_105 (13,312 cells, AR 3.25), R4 = 65_209
(53,248, AR 3.25); mislabeled meshes deleted, replacements generated + verified exact.
AR spread across the ladder is now 3.25–3.29 (±0.7% about the anchor); N ratios per
rung 2.10/1.90/2.10/1.90/2.00/2.00 (geometric mean 2.0). Ladder table updated in
`ledger.md`, still TENTATIVE pending the R1–R3 cost check.

### 2026-08-13 — Stale DJI conventions found in rotor_hover.jl / Phase 0 smoke (scaling + thrust axis)

While adapting the Phase 1 driver: `examples/rotor_hover.jl` scales nodes by
`maximum(nodes[1,:])` and reports thrust as `force[2,:]` — both Phantom-mesh-era
lines (y axial, x radial) left stale when commit 5615ada switched the file to the
DJI mesh (x axial, y radial; `git log -L` confirms). Effect: geometry ~20× the
intended radius relative to the fixed absolute kerneloffset/Das lengths, and a
radial force component labeled CT — this is why the Phase 0 smoke CTs were ~1e-5.
Phase 0's conclusions (availability, relative solver timings at 6,240 panels) are
unaffected — it disclaimed physics — but the numbers must not be reused as physical.
The Phase 1 driver (`benchmark/rotor_hover_solver_phase1.jl`) uses the corrected
conventions from `examples/rotor_hover_convergence.jl` (scale by radial dim 2,
rotation about x, thrust = −force_x). `examples/rotor_hover.jl` itself NOT edited
(fix pending Ryan's ruling, since downstream scripts copy it).

### 2026-08-13 — FMM residual floor measured (R1–R3): p-saturated, kerneloffset-dominated

Driver built: `benchmark/rotor_hover_solver_phase1.jl` (per-rung, one process per
rung; RUNG env; appends `benchmark/results/phase1/<mode>/runs.csv` + `floor.csv`,
never truncates). Frozen apparent-velocity RHS, Backslash reference with separately
timed assembly/factorization/solve, then ruling-11c floor sweep over expansion
order p ∈ {4,6,8,10,12,14} at production knobs (MAC 0.4, leaf 20, koff = R·1e-3).
Local dev runs, 4 threads (multi mode), commit b251071-dirty.

**Floor (rel = RMS(r)/RMS(b), direct-evaluated Backslash residual ~1e-13 rms):**

| Rung | N | rel floor (p ≥ 8, saturated) | max residual |
| --- | --- | --- | --- |
| R1 | 6,336 | 8.5e-7 | 2.3e-4 |
| R2 | 13,312 | 1.8e-6 | 4.1e-4 |
| R3 (July file) | 25,344 | 5.0e-5 | 8.3e-3 |
| R3 (fresh 20260813 file) | 25,344 | 1.3e-5 | 2.8e-3 |

Two findings, both decision-relevant:

1. **The floor saturates in p** (identical to 3 digits from p = 8 up) — ruling 11c's
   "tune `expansion_order` up" contingency does nothing. Diagnosis (R1
   localization + sensitivity): residual energy is mid-blade (r/R 0.7–0.9), NOT on
   shedding panels (0.05% of ‖r‖²); floor scales strongly with kerneloffset
   (koff 1e-2R → 9.0e-3 rel; 1e-3R → 8.5e-7; 1e-4R → 7.4e-9), with MAC
   (0.25/0.4/0.6 → 6.6e-8/8.5e-7/8.0e-6) and leaf size (20/50/150 →
   8.5e-7/2.9e-7/7.3e-9). **Root cause: the direct backend evaluates the
   offset-regularized kernel for every pair, while multipole expansions represent
   the unregularized kernel — every expansion-mediated pair carries an O(koff)
   discrepancy that expansion order cannot remove.** The floor also grows with N at
   fixed absolute koff (more relatively-close pairs enter expansion range), so the
   provisional 1e-6-relative target is unreachable at R3+ with production knobs
   (leaf 20, MAC 0.4). Consequences: (a) matched-residual targets must be set
   against the floor measured with each config's OWN FMM knobs (e.g. FGS at
   leaf 150 has a ~100× lower floor than the leaf-20 backend); (b) the tolerance
   freeze needs a ruling: raise leaf/lower MAC (changes the benchmark subject's
   cost) vs. set per-rung floor-relative targets. koff is case physics
   (production R·1e-3) — not a legitimate lever.
2. **The July 45_145 file is dirtier than a fresh regeneration** (4× higher floor,
   same N and topology). Ladder frozen on the fresh 20260813 file for
   generation-pipeline uniformity (all seven rungs same-day, same pipeline) —
   flagged for Ryan since the July file is the historical workhorse.

**Matrix-ful cost check (4-thread local, min-of-3):** assembly 1.03/3.73/13.3 s and
LU 0.92/7.59/55.1 s for R1/R2/R3 → fitted exponents ≈ N^2.0 (assembly) and N^2.9
(LU). Projected: R4 ≈ 70 s + 8.5 min, 23 GiB; R5 ≈ 4.7 min + 1.1 h, 82 GiB; R6
≈ 19 min + 8 h, 328 GiB; R7 1.3 TiB. **Ceiling set: matrix-ful configs run
R1–R5 (dense G fits a large HPC node and LU ≤ ~1 h), drop out at R6–R7** —
recorded in `ledger.md`; exact HPC timings land in Phase 2.

### 2026-08-13 evening — Ryan rulings + implementation: radius fix, Dirichlet case switch, tune_fmm

Ryan rulings (plan approved, plan file
`~/.claude/plans/we-ll-need-to-tune-playful-kettle.md`): (1) fold the
regularization reach into the radius FastMultipole sees, ON by default; (2) tune
the FMM per case with `tune_fmm` + a one-at-a-time perturbation double-check;
(3) **the campaign case matches `rotor_hover_pressure_comparison.jl` (RHPC), i.e.
Dirichlet** — capped captess4 watertight meshes, `Union{ConstantSource,VortexRing}`,
DBC=true, kerneloffset_panel=R·1e-10 / targets=1e-3, rho 1.179, RPM 6000;
(4) `examples/rotor_hover.jl`'s stale lines left as-is for now.

**Radius fix (src, on by default).** New kernel-dispatched
`radius_inflation(kernel, koff, tol)` (`src/FLOWPanel_elements_fmm.jl`, beside
`regularize`): source/doublet → exactly `koff` (compact support); VortexRing →
`koff·(2/tol)^(1/4)` (Vatistas quartic error). Applied in
`FastMultipole.source_system_to_buffer!` for panel bodies
(`src/FLOWPanel_abstractbody.jl`) using the ACTIVE `body.kerneloffset` (buffers
fill per pass, and `_set_kerneloffsets!` selects the pass's offset immediately
before each influence call — solve passes therefore pay only the panel offset,
target passes the targets offset), plus a one-time `@warn` when the inflation
exceeds 10× the geometric radius. Tolerance = module-global
`pnl.FMM_RADIUS_TOL[]` (default 1e-6; `Inf` disables → pre-fix radii). Same fix
for `PanelWake` panels (`+ radius_inflation(TK, core_size, tol)`) and the
filament writer (was `+ core_size`, under-inflated). Unit tests added
(`test/runtests_unit_fmm.jl`, testset "kerneloffset-aware FMM radius"): formula
values, buffer-radius assertions incl. active-offset policy and `Inf` escape,
PanelWake radius, and a large-koff (0.4) FMM-vs-direct regression with a
premise guard (uninflated MUST fail the bound — it does; inflated passes). All
17 pass; pre-existing FMM testsets green.
**A/B on the old Neumann case (koff=R·1e-3, the case that exposed the floor) —
MEASURED**: `RADIUS_TOL=Inf` reproduces the morning's saturated floor exactly
(rel 8.53e-7, p-flat); the default (1e-6) drops it **32× to 2.67e-8** with
p-dependence restored (7.2e-8 at p=6 → 2.7e-8 at p≥10), at ~+0.1 s per apply
(0.27→0.38 s at R1). On the new Dirichlet case (panel offset R·1e-10) the A/B
measured IDENTICAL floors with inflation on/off, as expected — the inflation is
inert when the active offset is tiny.

**Ladder re-frozen capped** (ledger.md): same U_W pairs, `_capped_captess4`
files, measured N = blade + caps: 8,016 / 15,760 / 28,752 / 58,192 / 108,240 /
212,108 / 419,276; watertight verified; RHPC shedding recipe (end_node-anchored
full TE trace, bbox for blade separation only, root clip r/R 0.1, cap-wrap
guard) — 22 edges → 20 after clip at R1.

**Driver rewritten for Dirichlet** (`benchmark/rotor_hover_solver_phase1.jl`):
RHPC parameters; frozen b captured per the `assemble_rhs!` contract (sources
from BC → direct source-potential assembly → b = −potential; solution =
strength col 2); reference x_ref via `_solve!` against exactly that b (direct
true residual ~6e-17 rms ✓); timed solves via the production tuple path with the
note that Dirichlet `solve!` reassembles the source potential internally (~6.5 s
of the 11.5 s at R1). `RADIUS_TOL` env knob; `floor.csv` gains a `radius_tol`
column (old Neumann rows rotated to `floor_neumann_uncapped_20260813.csv`).
New `TUNE=1` stage: `tune_fmm` at `PowerAbsolutePotential(0.1·target_rel·rms_b)`
+ 7-variant perturbation table → `tune.csv`. New availability script
(`benchmark/rotor_hover_solver_phase1_availability.jl`) re-proves the 7-config
roster on the Dirichlet system (B5; Phase 0 proved Neumann only).

**Two FastMultipole bugs found & fixed by the tune stage** (dev checkout,
commit 5adde3b + these): (a) `fmm.jl` tune-path interaction counting builds a
zero-step range when `direct_list` is empty (`n = max(..., 1)` fix); (b)
`tune_fmm`'s cache-preallocation call did not forward `kwargs...`, so the cache
buffer layout mismatched any tuning call with non-default derivative flags
(`autotune.jl` fix). Flagged for Ryan — these live in the FastMultipole repo.

**Dirichlet floors (capped ladder, MAC 0.4 / leaf 20 / koff_panel R·1e-10)**:

| Rung | N | rel floor (p ≥ 8) | behavior |
| --- | --- | --- | --- |
| R1 | 8,016 | 7.4e-7 | p-dependent to p≈8, then plateau |
| R2 | 15,760 | 2.24e-6 | plateau |
| R3 | 28,752 | 2.41e-5 | plateau; max-norm 3.07e-5 flat from p=4 (localized) |

The koff mechanism is gone (A/B-inert), but a **second, smaller p-saturated
mechanism remains** and already exceeds the provisional 1e-6 relative target at
R2. It is NOT kerneloffset (offset ≈ 0 here) — origin unknown (suspects: panel
multipole moments vs true panel integrals, self/near geometry at cap junctions).
**Tuner interaction**: at R2 the requested absolute tolerance (1.18e-9) sits
below this floor, so `tune_fmm` saturates at max P=20 with warnings and returns
knobs (leaf 12, MAC 0.5) whose measured floor (1.2e-5) is WORSE than the
production knobs' 2.24e-6 — tolerance-constrained tuning is degenerate below the
operator floor. At R1 (target above floor) the tuner returned p=17/MAC 0.5/leaf
21, and the perturbation check showed MAC 0.6 17% faster and leaf 10 15% faster
at identical measured floor — the double-check Ryan mandated is doing exactly its
job. **Tolerance freeze remains blocked on a Ryan ruling**: per-rung
floor-relative targets vs root-causing the residual plateau first.

Dirichlet cost check (4-thread local): assembly 5.7/31.6 s, LU 1.7/18.2 s at
R1/R2 — Dirichlet assembly is ~5× the Neumann cost at matched N (2-kernel
influence), LU ~2× (bigger N from caps); R1–R5/R6–R7 ceiling split unchanged.
**R3 completed 2026-08-14 00:3x** (after three background launches were killed by
the environment at startup; a foreground-rolled-to-background run succeeded):
floor rel 2.41e-5 (table above); assembly 71.6 s, LU 71.9 s, timed solve 143.8 s
(incl. ~72 s Dirichlet source-potential assembly), dense G 6.16 GiB, reference
residual 1.2e-16. The R3 TUNE stage was skipped (tuner is degenerate below the
floor anyway — see R2); rerun `RUNG=R3 TUNE=1` after the floor ruling if wanted.
**Secondary-floor growth across the ladder: 7.4e-7 → 2.24e-6 → 2.41e-5 (~N^1.8
in rel RMS)** — at R3 it is 24× the provisional target; the max-norm error is
p-flat from p=4 (localized), consistent with a few worst panels dominating.

**Roster availability on Dirichlet R1 (B5, measured; `availability.csv`)**:
`backslash_iterative` PASS (rel 5e-15); `krylov_gmres` / `krylov_jacobi` /
`krylov_ilu` PASS at smoke tolerance (rel 5e-5–9e-5). **Findings for Ryan**:
(1) `fgs` and `fgmres_fgs` FAIL — `BoundsError` in a 2×147 Int32 matrix inside
the FastGaussSeidel machinery: it assumes the Neumann single-strength layout and
cannot operate on the 2-strength DBC system as-is; (2) `backslash_coupled` runs
but its solution measures rel_rms = 1.0 against the frozen b (x ≈ 0 in strength
col 2) — likely a solution-column/coupling-layout mismatch on DBC bodies;
treated as NOT proven, needs root-cause. Roster reshaping (or FGS Dirichlet
support) is a Ryan decision, per the plan.

**Suites**: `runtests_unit_fmm.jl` (incl. 17 new radius tests) and
`runtests_unit_wake.jl` both green.

### 2026-08-14 — "Secondary floor" ROOT-CAUSED AND FIXED: it was a direct-kernel defect, not the FMM

Per Ryan's directive (single-panel isolation), the koff-independent p-saturated
floor was chased down on R1 Dirichlet:

1. Localization: full-residual error concentrated on a handful of targets (worst:
   tiny tip panel 7934); single-source columns initially looked exact, but a
   binary search over source subsets isolated the entire worst-target error to
   ONE source panel (6652, a TE-strip sliver, AR ≈ 24, NOT a shedding panel —
   attached-wake hypothesis refuted by `suppress_attached_wake` A/B).
2. The (6652 → 7934) pair: error = −1.0003 × the pair's direct influence
   (pair effectively missing), deterministic, bit-identical across leaf sizes,
   gone at MAC 0.3, p-converged. Tree/list inspection: the pair IS covered by a
   legitimately admitted M2L; 99 of 101 leaf-mates deliver correctly through the
   same expansion chain.
3. **Ground truth by quadrature inverted the picture**: the target sits 2.4 μm
   from the sliver's PLANE, 0.034 m away along its long-edge extension; the
   refined-quadrature dipole potential is +2.9255e-9 — the FMM matches it to
   7 digits, while the ANALYTIC DIRECT KERNEL returned −1.0466e-5, wrong by
   3,600×. The dense G (and hence the Backslash reference and the RHS) inherits
   the same wrong entries — the measured "floor" was the direct/dense operator's
   error, with the multipole expansion CORRECT.
4. **Root cause** (`src/FLOWPanel_elements_fmm.jl`, `compute_source_dipole`, 3
   sites): the tan_term guard `abs(abs(R_dot_s) − ri·ds) < 1e-12` uses an
   ABSOLUTE tolerance in units of length², i.e. a trigger tube of half-angle
   ~√(1e-12/(ri·ds)) rad around every panel-side extension line — ~1e-4 rad for
   TE slivers at 0.03 m range. Inside it, one edge's solid-angle term is zeroed,
   breaking the three-edge cancellation. Explains everything measured:
   p-independence, MAC dependence (both operators share the defect when the pair
   is direct → discrepancy invisible), the on/off behavior between adjacent
   slivers, and the N-growth (θ_crit ∝ 1/√(edge·dist) grows as panels shrink).
5. **Fix**: relative tolerance `≤ 1e-12·ri·ds` (`<=` so the degenerate
   target-on-vertex case, ri·ds = 0, still triggers; the zero limit is exact ON
   the extension line). Regression test added
   (`runtests_unit_fmm.jl` testset "sliver panel near-plane far-field
   potential"): DJI-like sliver + near-plane extension target vs in-test
   quadrature — analytic now matches to 6 digits (was 10³× off).
6. **Result**: the floor is GONE on all three rungs (driver reruns recorded in
   floor.csv, 2026-08-14; MAC 0.4, leaf 20, radius_tol 1e-6):

   | rel RMS | p=6 | p=8 | p=10 | p=12 | p=14 |
   | --- | --- | --- | --- | --- | --- |
   | R1 (8,016) | 6.4e-6 | 2.9e-7 | 1.4e-8 | 2.2e-9 | 3.4e-10 |
   | R2 (15,760) | 3.0e-6 | 2.1e-7 | 1.8e-8 | 1.1e-9 | 8.0e-11 |
   | R3 (28,752) | 4.5e-6 | 2.8e-7 | 2.3e-8 | 2.2e-9 | 2.2e-10 |

   Clean p-decay, rung-INDEPENDENT (the old N-growth is gone with the defect).
   The provisional 1e-6-relative target clears at p≈8 with headroom; p=10 gives
   ~50× margin. Dirichlet cost check (4-thread local): assembly
   8.5/34/101 s, LU 2.8/18/105 s at R1/R2/R3. Suites green: kernel_gradient,
   fmm (18/18 incl. new test), solver, liftingbody (with helpers),
   instrumentation, wake.
   NOTE: all previous floor tables above measured THIS defect, not FMM error;
   the tuner-degeneracy observations remain valid in principle (tolerance below
   any true floor) but the specific numbers are superseded. The defect also sat
   in every historical dense/direct evaluation (same guard feeds the velocity
   z-component); the post-fix shift in global quantities looks tiny (rms_b moved
   at the 12th digit at R1) but any past result hinging on near-plane
   edge-extension pairs deserves a spot check.

### 2026-08-14 — FGS-on-Dirichlet failure diagnosed (fix suggested, not applied)

Ryan pointer confirmed: FGS runs Dirichlet fine in the duct examples
(`duct_vpm.jl` etc., `Union{ConstantSource,VortexRing}`, default DBC=true) — the
break is NOT the BC. Reproduced with full stack: BoundsError in
`FastMultipole.sort_by_target` called from `add_self_interactions` inside the
`FastGaussSeidel` constructor (`FastMultipole/src/solve.jl:529,585`;
`interaction_list.jl:612`). Mechanism, verified by building both trees for the
R1 Dirichlet body (leaf 150): **target tree 229 branches / source tree 147** for
the same bodies — FastMultipole inflates SOURCE tree bounding boxes by per-body
radii (buffer row 4) while TARGET trees use bare positions, so the octrees
diverge; `FastGaussSeidel` implicitly assumes structurally identical trees (it
pushes source-leaf indices into the TARGET slot of the direct list for its self
blocks, and sizes the sort counter by `source_tree.branches` while the list's
target indices refer to the target tree, up to 229 > 147 → BoundsError at 166).
The duct/Neumann cases worked because the trees happened to coincide (Phase 0's
matched-residual check confirms those solves were consistent); the capped DJI
mesh's tiny tip-cap panels broke the coincidence. **Suggested fix (Ryan's
call, FastMultipole internals)**: make `FastGaussSeidel` single-tree by
construction — reuse the source tree('s structure) for the target side, or
build the interaction lists with the source tree on both sides since
target ≡ source bodies; a mere counter-resize would silence the crash but leave
silent block mis-association whenever the trees differ. Until then the roster's
`fgs`/`fgmres_fgs` configs stay unavailable on the campaign case.

**STAGED (Ryan 2026-08-14): fix = "one topology, two shrink results" (option 2
of the discussed set), NEXT ACTION for this item.** Design agreed with Ryan:

- In `FastGaussSeidel` (`../FastMultipole/src/solve.jl:585` ctor), the target
  tree must share the SOURCE tree's octree topology (same root box, same
  subdivision decisions, same `bodies_index` partitions), then each copy runs
  its role-appropriate shrink pass so the error machinery keeps tight per-role
  centers/radii (naively reusing the source tree would make target radii loose —
  panel radius + regularization inflation — degrading MAC/error-bound tightness;
  that is why plain "single tree" was rejected).
- Implementation route: build the source tree first, then construct the target
  tree by REPLAYING its splits instead of re-deriving them from positions
  (structure copy + target-buffer fill + target shrink). Check first whether
  `Branch` already carries separate source/target radius slots (its ctor
  appends two `zero(radius)` fields, containers.jl:181) — if so, a shared
  structure holding both radii is the same idea with less duplication.
- Either way, assert the invariant in the ctor (equal branch counts +
  `bodies_index` equality) so divergence errors loudly.
- Rationale trail: crash at `sort_by_target` (interaction_list.jl:612) from
  `add_self_interactions` (solve.jl:529) mixing tree index spaces; R1-Dirichlet
  capped body at leaf 150 measures 229 target vs 147 source branches (asymmetry:
  source boxes are radius-inflated, tree.jl:1428; Barba `exceeds` is
  role-independent, tree.jl:896).
- Regression test: build `FGSSolver` on the capped R1 Dirichlet ladder body
  (`dji9443_20260813_23_73_capped_captess4.msh`, construction recipe in
  `benchmark/rotor_hover_solver_phase1.jl` / `_availability.jl`) — the mesh
  that provably splits the trees — then one solve + true-residual check vs the
  frozen b; rerun the availability script to flip `fgs`/`fgmres_fgs` to PASS.
- This is FastMultipole (Ryan's package, dev checkout `../FastMultipole`,
  branch main, dirty): land there with its own unit test, mirroring how the
  two approved tune_fmm fixes were handled.
- **PENDING CHECK (Ryan 2026-08-14): re-examine error / `tune_fmm` after the
  tree fix verifies** — the FGS target tree now replays the (shallower,
  radius-limited) SOURCE topology instead of an independently subdivided
  target tree, so `tune_fmm`'s error model / knob search may no longer behave
  as intended on FGS configs (tuning assumptions about tree depth, MAC/error
  bounds, and near-field block sizes can shift). Check once the availability
  rerun lands, before freezing any FGS tolerances.

### 2026-08-14 — FGS tree fix LANDED (FastMultipole) + availability VERIFIED: fgs/fgmres_fgs PASS

Implemented the staged "one topology, two shrink results" fix in
`../FastMultipole` (uncommitted working-tree changes alongside the two tune_fmm
fixes, per convention; fm_commit 5adde3b-dirty):

- **True divergence mechanism (sharper than the earlier diagnosis)**: not the
  root-box inflation — `center_box` (def. tree.jl:~1368 pre-fix, called from the
  ctor) is positions-only for both roles and Barba `exceeds` (post-fix :949) is
  role-independent. The cause is the **source-only subdivision stop guard** in
  `branch!` (post-fix tree.jl:759; also the multithread variant :796 and the
  `(target || child_radius >= max_body_radius)`
  forms in the `child_branches!` variants :512/:583/:678/:711): subdivision
  halts when the child half-radius would drop below the branch's largest
  per-body radius (census buffer row 4). Radius-inflated panels therefore make
  the source tree strictly SHALLOWER (R1: 147 branches) than the position-only
  target tree (229). Corollary: within-octant sort order is position-only, so
  the source `sort_index` is a valid target ordering — the replay is clean.
- **Fix**: new replay constructor `Tree(source_tree, target_systems, switches;
  shrink, recenter)` in `tree.jl` (+53) — copies branches/levels_index/
  leaf_index/sort indices verbatim, fills target buffers directly in
  source-sorted order via `target_to_buffer!`, then runs the role-appropriate
  `update_min_influence!` + `shrink_recenter_target!` so target centers/radii
  stay tight for the error machinery ("two shrink results"). `FastGaussSeidel`
  ctor (solve.jl, +32/−3) now builds the target tree this way and calls
  `assert_shared_topology` — a loud `error()` (not compilable-out `@assert`)
  on branch-count or per-branch `bodies_index`/`branch_index` mismatch.
- **Tests**: new testset "Fast Gauss Seidel: shared source/target tree
  topology" in `test/solve_test.jl` with a NON-VACUOUS premise guard (synthetic
  Gravitational system, radius_factor tuned to 4.0 so independent source/target
  builds provably diverge; the first attempt 0.3 did NOT diverge and the guard
  caught it). Full FastMultipole suite at 4 threads: exit 0, zero failures.
- **Design note for Ryan**: the shared topology is the (shallower) SOURCE tree,
  so FGS near-field blocks are larger than an independent target tree would
  give at the same leaf_size — inherent to the agreed design. The invariant
  check is O(branches) per ctor (negligible, permanent).

**Verification on the campaign case (B5 availability rerun, multi @4 threads,
2026-08-14 06:37–06:41, `availability.csv` rows of record)**: `fgs` PASS
rel_rms 3.54e-6 (rms residual 4.07e-8), `fgmres_fgs` PASS rel_rms 6.29e-5 —
both formerly BoundsError — at smoke settings (rtol 1e-4/itmax 400); all other
configs unchanged (backslash_iterative 5.4e-15, Krylov trio 5e-5–9e-5).
**`backslash_coupled` still measures rel_rms = 1.0 (x ≈ 0)** — the script's
PASS only means it ran; the coupling-layout suspicion stands, still item (3)
in the next-actions order. Full seven-config roster is now available on the
Dirichlet campaign case.

**Clear-context review 2026-08-14: APPROVED, zero blocking findings** (fresh
subagent, code + diff + live single-threaded testset run 5/5 + CSV
cross-check). Verified: stop-guard mechanism claims; replay-ctor correctness
(Branch immutability makes tree aliasing safe; shrink=false role-appropriate
since unshrunk cubic geometry is role-independent; buffers zero-initialized;
expansions sized per tree); premise guard genuinely non-vacuous (pre-fix the
branch-count test would fail); evidence numbers match `availability.csv`
exactly; FLOWPanel src untouched by the fix; both prior tune_fmm fixes
undisturbed. Non-blocking notes: (a) doc line citations corrected above
(post-diff shift +53); (b) the testset's final `solve!` converges in 0
iterations (strengths already exact after direct!-based setup) — exercises
the old BoundsError path but is not a nontrivial solve; strengthening is
optional; (c) `assert_shared_topology` is tautological on the replay path by
construction — its value is guarding future regressions; the pre-existing
`target_systems === source_systems` assert closes the system-mismatch hole;
(d) housekeeping for Ryan pre-commit: `scripts/20250404_prediction_accuracy.jl`
is dirty in the FastMultipole working tree from unrelated error-method
experimentation — segregate before committing the three pending fixes;
`test/Project.toml` gains FastMultipole as a test dep (benign, needed).

### 2026-08-14 — TUNE=1 rerun R1–R3 against clean floors; tuner healthy; tolerance-freeze PROPOSAL (judgment call)

Reruns (multi @4 threads, `tune.csv`/`floor.csv` rows 2026-08-14 06:5x–07:3x).
Floors reconfirmed identical to the post-sliver-fix table (e.g. R3 rel 2.34e-8
at p=10, production knobs MAC 0.4/leaf 20). **The tuner is no longer
degenerate**: pre-fix it saturated P=20 at R2 even at MAC 0.4/0.5; now it
converges at MAC 0.3–0.5 on every rung and only fails its bound at MAC
0.6/0.7 (P=20 warnings), which it correctly skips.

| Rung | tuned (p, MAC, leaf) | t_apply | rel_rms | notable perturbations (all target-meeting) |
| --- | --- | --- | --- | --- |
| R1 | 17, 0.5, 21 | 1.60 s | 1.20e-10 | MAC 0.6: −16.6% @ 4.7e-8; leaf 10: −14.3% |
| R2 | 17, 0.5, 12 | 3.50 s | 5.01e-10 | MAC 0.6: −19.4% @ 2.7e-8; leaf 6: −8.9% |
| R3 | 18, 0.5, 18 | 9.39 s | 2.00e-10 | MAC 0.6: −21.1% @ 1.4e-8; leaf 9: −13.4% |

Tuner requests abs tol = 0.1·1e-6·rms(b) (rel 1e-7); achieved rel_rms lands
200–500× below the 1e-7 request — the error bound is conservative ~2 orders.
Consistently on all rungs, MAC 0.6 measures 17–21% faster at rel_rms 1.4–4.7e-8
(still ≤1e-7) but is BOUND-UNCERTIFIED (P=20 warning). Note the tuned points
also beat the production knobs outright: at R3, tuned apply 9.4 s vs 12.3 s
(p=10/MAC 0.4/leaf 20) at 100× lower error — tuned knobs dominate.

**Judgment call (autonomous, Ryan 2026-08-14 directive "make your own best
judgment, log, report"; override freely):**

1. **Matched true-residual target frozen at rel 1e-6 per rung** (the
   provisional value) — post-fix floors clear it by ≥40× at production knobs
   and ~5000× at tuned knobs; every iterative config can honestly stop on it
   (right-preconditioned, ruling 11a).
2. **FMM-backend apply knobs for R1–R3 = the per-rung TUNED values above**
   (they dominate production knobs in both cost and error). R4–R7: run the
   same TUNE procedure on HPC per rung (tuner needs no dense reference) and
   adopt its output under the same rules.
3. **Bound-uncertified faster points (MAC 0.6) REJECTED** despite measuring
   in-tolerance at R1–R3: at R6–R7 no dense reference exists, so the error
   bound is the only guard; the campaign operates only at bound-certified
   points. Cost: ~20% apply time. Ryan may override (e.g. accept MAC 0.6 with
   per-rung empirical floors where measurable, R1–R5).
4. **FGS (and fgmres_fgs preconditioner) knobs/tolerances NOT frozen** —
   gated on the PENDING CHECK above (tune_fmm error model vs the
   replayed-source-topology FGS trees). backslash_coupled remains
   quarantined until its x≈0 root-cause (item 3).

Copied to `decision_rules.md` as a PROPOSED block (marked pending Ryan).

### 2026-08-15 — Ryan rulings: BC-error metric + per-rung/per-solver tuning; plan STAGED

Ryan (in session, 2026-08-15): (1) primary metric = BC error (induced
potential at CPs vs the BC), rel L2 target 1e-6, FMM-evaluated with the
error-tolerance feature (reference-free — supersedes the strength-agreement
thresholds and moots the relL2 3e-5 proposal); (2) NO production-knob
freeze — tune FMM knobs per rung per solver for optimal performance
(supersedes the 2026-08-14 PROPOSED freeze); (3) share one tuned apply-knob
set across all Krylov configs, and one FGS knob set between solver and
preconditioner roles, the latter selected by a τ-ladder (tune FGS at a range
of coarser tolerances, measure end-to-end FGMRES runtime per τ, pick best);
(4) FGS 1e-6 landing: revisit after tuning {p, MAC, leaf, inner}; if it
still overshoots, report both 1e-6-target cost and the last iterate above
1e-6 (error snapped to nearest half decade). Consequence: the MAC≥0.6
bound-certification rejection is withdrawn (per-run certified BC verification
is the guard at every rung). Full staged plan (approved to stage only, NOT
executed): `phase_01_tuning_plan.md`. Baseline production-knob per-solve
timings completed same day (`solvetime.csv`: gmres 94/226/651 s,
jacobi 31/49/151 s, ilu 11/35/104 s at R1/R2/R3; fgs 2.6/5.7/20 s and
fgmres_fgs 17/33/101 s in fgs_check.csv) — kept as the untuned baseline.
Pending small item: backslash_ldiv rows (config added to the solvetime
script post-launch, not yet run).

### 2026-08-14 evening — backslash_coupled "x≈0" ROOT-CAUSED: script state leak + a real latent rhs-accumulation defect (both fixed, test-gated)

Ryan away; autonomous session per his 2026-08-14 directive (judgments logged
here, morning report drafted, no notebook writes).

The availability rows' `rel_rms = 1.0` was NOT x≈0 and NOT a
solution-column/coupling-layout mismatch. Probe on the R1 Dirichlet case
(scratchpad `probe_backslash_coupled.jl`, 1 thread, direct-evaluated
residuals vs the frozen b; b-capture + Backslash x_ref exactly as the
driver):

- **A (availability reproduction)**: entry `rotor.potential` = the frozen
  source-potential assembly → `rel_rms = 1.0` with **‖x‖/‖x_ref‖ = 2.0 and
  ‖x − 2·x_ref‖/‖x_ref‖ = 0.0** — the solution was exactly DOUBLED, and
  G·(2x_ref) − b = b makes rel_rms exactly 1.0 (the availability
  signature). Mechanism: the coupled solvers (`BackslashCoupled` solve!
  ~line 1709, `KrylovCoupled` alike) treat entry `body.potential` as an
  EXTERNAL potential and add it to their internally re-assembled source
  potential; the availability script's `reset!` never zeroed
  `rotor.potential`, which still held the b-capture's assembly → Dirichlet
  RHS doubled. Single-body solve paths zero the workspace themselves,
  which is why only `backslash_coupled` was affected.
- **B (contract respected)**: entry potential zeroed, fresh solver →
  rel_rms 5.4e-15, x ≡ x_ref (‖xB − x_ref‖ = 0). **The coupled
  solution-column/coupling layout on DBC bodies is CORRECT.**
- **C (latent src defect, real)**: SECOND solve on the same solver (entry
  potential 0) → x = 2·x_ref again. Root cause: `calc_bc_dirichlet`
  ACCUMULATED (`RHS .-= self.potential`) into the coupled solvers'
  persistent `solver.rhs`, which is never zeroed between solves — every
  solve after the first doubled the Dirichlet rows (the Neumann path
  assigns, `RHS[i] = ...`, and was safe). Affects `BackslashCoupled` AND
  `KrylovCoupled` on Dirichlet bodies in any repeated-solve use
  (time-marching, warmstart studies).

**Fixes landed** (judgment: root cause unambiguous — probe predictions
matched to machine precision — so fix now rather than quarantine, per the
work-queue rule):

1. src: `calc_bc_dirichlet` now assigns (`RHS .= .-self.potential`),
   matching Neumann semantics (`src/FLOWPanel_solver.jl`). Only call site
   is the coupled boundary_condition! path.
2. Regression test: new testset "Coupled Dirichlet rhs idempotence
   (regression)" in `test/runtests_unit_solver.jl` — BackslashCoupled and
   KrylovCoupled on a DBC octahedron, solve twice, both solves must equal
   the single-body Backslash reference (pre-fix, solve 2 = 2× solve 1
   fails it).
3. Availability script: `reset!` now zeroes `rotor.potential` (with a
   comment recording the contract). The 2026-08-13/14 `backslash_coupled`
   availability rows measured this script bug — annotate, keep rows.

Suite runs + availability rerun queued behind the FGS-check rungs (4-thread
budget). backslash_coupled QUARANTINE LIFTED pending green tests + a clean
availability row.

**Verification landed (same evening):** full `runtests_unit_solver.jl` suite
GREEN incl. the new regression testset (5/5); premise guard run with the
pre-fix `.-=` temporarily restored measured ‖x₂‖/‖x₁‖ = 2.0 exactly (test is
non-vacuous), fix re-applied; availability rerun (2026-08-14 ~19:5x rows):
**backslash_coupled PASS rel_rms 5.43e-15** — identical to
backslash_iterative, i.e. exact — with all six other configs unchanged.
Quarantine LIFTED; the seven-config roster is fully clean on the Dirichlet
campaign case.

### 2026-08-14 evening — PENDING CHECK RESOLVED: tune_fmm-vs-FGS-topology, empirical R1–R3; FGS knobs/tolerances selected (PROPOSED)

The check Ryan requested (does tune_fmm's output transfer to FGS configs now
that FastGaussSeidel replays the shallower SOURCE topology on the target
side?). New driver `benchmark/rotor_hover_solver_phase1_fgs_check.jl`
(construction mirrors the availability script; frozen b; direct-evaluated
true residuals; multi @4 threads; rows in `fgs_check.csv` 2026-08-14
19:07–19:4x). Design judgments: three knob candidates per rung — `prod`
(p=10/MAC 0.4/leaf 150, the Phase 0 production FGS knobs), `tuned` (the
per-rung tune_fmm output incl. its small leaf), `hybrid` (tuned p+MAC,
leaf 150 — because FGS's leaf doubles as the dense near-field block size,
which tune_fmm knows nothing about); FGS internal tolerance set ABSOLUTE
max-abs = 1e-6·rms(b) (its metric; max ≥ rms makes it conservative for the
RMS target); fgmres rtol=1e-6/atol=1e-14, memory=50 (avoids restart
confounds; 9–11 its used), FGS solve times exclude the N² Dirichlet
source-potential assembly (frozen assembly restored, not rebuilt).

**Result: ALL 18 runs (3 rungs × {fgs, fgmres_fgs} × 3 candidates) meet the
rel-1e-6 direct-evaluated target — the source-topology replay introduces NO
new residual floor**; internal convergence honestly implies the direct
target with ≥5% margin everywhere (≥17× margin for all `fgs` rows).

| Rung | config | p/MAC/leaf | iters | t solve (s) | direct rel_rms |
| --- | --- | --- | --- | --- | --- |
| R1 | fgs prod | 10/0.4/150 | 5 | 2.56 | 4.7e-8 |
| R1 | fgs tuned | 17/0.5/21 | 6 | 3.12 | 3.9e-8 |
| R1 | fgs hybrid | 17/0.5/150 | 6 | 2.68 | 7.4e-9 |
| R2 | fgs prod | 10/0.4/150 | 6 | 5.68 | 3.2e-8 |
| R2 | fgs tuned | 17/0.5/12 | 12 | 9.40 | 1.8e-7 |
| R2 | fgs hybrid | 17/0.5/150 | 7 | 6.39 | 1.7e-8 |
| R3 | fgs prod | 10/0.4/150 | 8 | 19.98 | 3.5e-8 |
| R3 | fgs tuned | 18/0.5/18 | 10 | 27.22 | 5.9e-8 |
| R3 | fgs hybrid | 18/0.5/150 | 9 | 19.33 | 4.6e-8 |
| R1 | fgmres_fgs prod | 10/0.4/150 | 9 | 17.17 | 6.7e-7 |
| R1 | fgmres_fgs tuned | 17/0.5/21 | 9 | 17.93 | 5.7e-7 |
| R1 | fgmres_fgs hybrid | 17/0.5/150 | 9 | 18.39 | 6.9e-7 |
| R2 | fgmres_fgs prod | 10/0.4/150 | 9 | 33.49 | 7.9e-7 |
| R2 | fgmres_fgs tuned | 17/0.5/12 | 10 | 41.19 | 4.1e-7 |
| R2 | fgmres_fgs hybrid | 17/0.5/150 | 9 | 33.96 | 8.2e-7 |
| R3 | fgmres_fgs prod | 10/0.4/150 | 10 | 101.11 | 9.5e-7 |
| R3 | fgmres_fgs tuned | 18/0.5/18 | 11 | 115.48 | 4.4e-7 |
| R3 | fgmres_fgs hybrid | 18/0.5/150 | 11 | 113.40 | 4.0e-7 |

(fgmres_fgs times include per-iteration tree rebuilds in the Krylov apply —
the known Phase 2b lever; precondition column = FGSPreconditioner knobs,
apply backend = per-rung tuned knobs throughout.)

Findings: (1) correctness transfers — tune_fmm knobs are SAFE inside FGS on
every rung; (2) but they are not optimal for FGS: the tuner's small leaf
degrades FGS convergence (R2: 12 outer its at leaf 12 vs 6 at leaf 150 —
leaf is also the Gauss–Seidel block size) and costs 20–65% more wall time;
(3) production FGS knobs are fastest or tied at every rung, with `hybrid`
statistically indistinguishable (±3%); (4) FGS-family error is guarded
per-run anyway: every campaign solve records a direct-evaluated true
residual (ruling 3), which is an O(N²) apply (no dense G), feasible at all
rungs incl. R6–R7 — so FGS needs no tuner bound certification.

**Judgment (autonomous, PROPOSED pending Ryan): FGS-family freeze** —

- `fgs`: p=10, MAC=0.4, leaf=150, inner_iterations=20, rlx=1.0,
  shrink=true, recenter=false, tolerance = 1e-6·rms(b) per rung (absolute,
  max-abs metric), max_iterations=300 cap (≤12 used). Same knobs at R4–R7
  (they were never rung-sensitive R1–R3); guard = the per-run
  direct-evaluated residual, target rel 1e-6.
- `fgmres_fgs`: preconditioner = the same FGS knobs (sweeps=1,
  inner_iterations=2); apply backend = the per-rung TUNED knobs (already
  frozen); rtol=1e-6, atol=1e-14, memory=50, itmax=500.
- tune_fmm's role for FGS configs: none beyond the apply backend — do NOT
  re-tune FGS-internal knobs per rung.

Copied into the `decision_rules.md` PROPOSED block.

### 2026-08-14 late evening — Agreement stage vs Backslash (exit criterion 1): R1–R3, all 7 configs, production steady! path

New stage `benchmark/rotor_hover_solver_phase1_agreement.jl` → `agreement.csv`
(multi @4T). Design judgments (autonomous): reference = `backslash_ref`
(Backslash, backend_solve=Direct) through the SAME `steady!` path (skeleton =
rotor_hover_convergence.jl); all configs share one post-processing/monitor
backend (per-rung tuned knobs) so CT deltas isolate the solved strengths;
backend_solve=Direct for every config's Dirichlet rhs assembly; iterative
settings = the frozen/PROPOSED values; CT from
PressureLaplace(:material_derivative, :edge_difference)+Force, steady
PressureBernoulli+Force, KuttaJoukowski — fresh monitors per config; steady-
Bernoulli's moving-body caveat accepted because the metric is the
cross-config DELTA at identical recovery; per-config direct-evaluated
residual vs the frozen driver b recorded (matched-residual verification).

Operational notes: (1) the previous config's steady! target pass leaves
kerneloffset_targets active — the script now restores the panel offset before
each solver ctor (else ILU/FGS ctor-time trees see radius-inflated panels;
the ILU direct-list pattern guard tripped at R1). (2) R3 required chunked
one-process-per-config runs: long background jobs on this machine get killed
(three kills, twice at the two-dense-G peak; 16 GB RAM) — reference solution
persisted to `agreement_xref_R3.bin` so heavy configs run standalone.
(3) **ILU pattern guard legitimately trips at R3** (default 512N; the
leaf-10/MAC-1.0 Barba direct list exceeds 512 entries/panel as the
sliver-heavy mesh refines — Phase 2 knob consideration for R4+); agreement
runs use max_pattern_entries=2048N, judgment logged. (4) PressureLaplace CG
hits its 1000-iteration cap at R3 (rel residual ~3e-6) — identical for every
config, so deltas stay comparable; raise itmax for the HPC sweep.

**Results (agreement.csv):** CT agreement is uniformly excellent — worst
|dCT| ≈ 8e-4 % (backslash_coupled), 120× inside the 0.1% threshold, across
all three monitors and all rungs. Strength agreement vs provisional
thresholds (relL2 ≤ 1e-5 / relmax ≤ 1e-4): jacobi/ilu/fgs/fgmres_fgs PASS
with large margins (fgs best at relL2 ~2e-8); krylov_gmres relL2 =
1.2/1.5/1.6e-5 (R1/R2/R3) sits just above 1e-5 — exactly the
residual→solution amplification κ_eff·rtol with κ_eff ≈ 12–16 (its true
residual honestly meets 1e-6). **Judgment (PROPOSED): calibrate the final
threshold to relL2 ≤ 3e-5 at the rel-1e-6 target (per the decision_rules
calibration contract), relmax ≤ 1e-4 and dCT ≤ 0.1% unchanged.** Alternative:
tighten the matched target to 3e-7 if the 1e-5 solution bar is wanted.

**FLAG → ROOT-CAUSED same night: backslash_coupled inside steady! is a
FORMULATION-GAUGE inconsistency, not a solver bug.** Observed: direct
residual vs the frozen b 4.6–4.8e-5 with relL2 ≈ 1.3e-5 on every rung, while
the SAME solver called directly is exact (5e-15). Probe chain (scratchpad
`probe_coupled_steady.jl` → `probe_coupled_das.jl` → `probe_coupled_rhs.jl`):

1. Trisection: fresh-body steady! reproduces 4.778681231e-5 BIT-identically;
   prior-config state irrelevant; deterministic in the steady!+coupled
   combination.
2. Internals diff: Das unchanged, **G/LU factors bit-identical**, entry
   velocity and BC sources bit-identical — the ENTIRE difference is the rhs,
   and `solver.phi_ext` (entry potential snapshot) is nonzero with rel norm
   exactly the rhs delta.
3. Writer identified: `_steady_aerodynamics!` calls `apply_freestream!`
   (FLOWPanel_simulate.jl:489 → FLOWPanel_abstractbody.jl:1062), which adds
   the freestream scalar potential φ∞ = U∞·x into `body.potential`
   (magnitude check: magVinf 1e-4 × CP coords → norm 4.92e-5 ✓ exact).

Mechanism: the single-body Dirichlet solve! wrapper ZEROES `body.potential`
(perturbation-potential gauge — discards φ∞; this is the production
convention every historical result uses), while the coupled solvers honor
entry potential as an external field (total-potential-ish gauge — φ∞ enters
the Dirichlet rhs). Both are self-consistent formulations; they differ by
the φ∞-driven μ-gauge response (relL2 1.3e-5 at magVinf=1e-4, scaling with
|U∞|; CT deltas 8e-4 % = gauge×discretization). At hover (magVinf→0) the
difference vanishes. **Disposition = Ryan ruling**: unify the two paths'
Dirichlet external-potential convention (either zero φ∞ in the coupled path
inside steady!/simulate!, or honor external potential in the single-body
wrapper — the latter shifts production history). Not a Phase 2 blocker:
agreement is quantified and benign at campaign thresholds.

### 2026-08-15 evening — Tuning plan EXECUTION begins (autonomous session): Stage 1 evaluator built + validated

Executing `phase_01_tuning_plan.md` per Ryan's handoff (staged plan granted
execution). Ryan away; judgment calls logged here per standing directive.

**Judgment calls (this entry):**

1. **Tuning objective confirmed as carried judgment** (plan ruling 5): per-solve
   wall time subject to certified BC ≤ target, setup timed separately. Ryan not
   in session to confirm — proceeding under the staged judgment; override
   propagates mechanically (objective is one line in the tuner).
2. **Shared case construction factored into `benchmark/phase1_case.jl`** —
   the RHPC/Dirichlet case block was duplicated verbatim across four scripts
   ("keep in sync" comments); the three new Stage 1–3 scripts and any Stage 5
   regeneration include the shared file instead. The four pre-existing scripts
   are untouched (their CSVs are evidence of record).
3. **Correction to the 2026-08-15 reset-point language**: solvetime.csv's
   Krylov rows are NOT "production-knob" timings — their apply backend already
   used the per-rung TUNED knobs (notes column: `apply backend tuned(...)`).
   They are the tuned-apply baseline; only the FGS-family knobs and the
   backslash_ldiv rows remain to be (re)generated for the Stage 5 table.
4. **bc_error! gate logic fix after first R1 run**: the first validation run
   (bcerror.csv rows 2026-08-15 ~17:4x) mislabeled the krylov_ilu point
   (direct rel 8.29e-7, FMM 8.30e-7, disc 0.11%) as "below evaluator
   resolution" (cutoff was direct ≥ 1e-6) and then failed it against the
   floor criterion. Gate rewritten: PASS = error_success AND (disc ≤ 10% OR
   fmm ≤ direct + safety·target [additive requested evaluation tolerance]).
   R1 rerun with the fixed gate; both runs' rows kept in bcerror.csv.

**Stage 1 (certified BC-error evaluator) — implementation of record:**
`bc_error!` in `benchmark/common.jl`: one influence pass with BOTH strength
columns (σ re-derived via `set_strengths!` from the frozen apparent velocity;
col 2 = candidate x); control-point potential = φ_σ + G_μ x = G x − b = the
BC residual; rel L2 normalized by rms_b (= RMS φ_σ). FMM path:
`FastMultipole.fmm!` under `PowerAbsolutePotential(0.1·target·rms_b)` with
dynamic-P cap 20, per-rung tuned MAC/leaf; **certification = the fmm! return's
`error_success` flag** (every M2L met its bound at P ≤ cap). Direct path
(:direct) = O(N²) validation reference. Side-effect-free (full state
save/restore incl. kerneloffset forced to panel value).

**R1 validation result (bcerror.csv, first run — numbers identical in rerun):**
6 points spanning 1e-15…2e-2 (backslash, krylov_ilu, fgs, fgs stopped at
1/2/3 outer iterations): discrepancy FMM vs direct ≤ 0.88% at every
resolvable point (0.11% at the 1e-6 scale), error_success true everywhere;
the evaluator's own floor measured 4.0e-9 rel on the backslash point (25×
below the 1e-7 request — conservative-bound behavior consistent with
tune.csv). Direct-vs-true_residual! cross-check column agrees (same operator
algebraically). R2/R3 validation runs queued.

### 2026-08-15 late evening — Stage 1 gate PASS R1–R3; Stage 4 (decision_rules) + Stage 6 (HPC spec) delivered

**Stage 1 COMPLETE** (bcerror.csv, logs bcerror_R{1,2,3}.log): all three rungs
GATE PASS with the corrected gate. Worst FMM-vs-direct discrepancy across all
18 points: 3.72% (R2 fgs point at 3.2e-8 — near the evaluator floor); at
every point ≥1e-5 the discrepancy is ≤0.04%; at the 1e-6 scale
(krylov_ilu points 3.3e-7–8.3e-7) ≤0.11%. `error_success=true` on every
evaluation (P cap 20 never exhausted). Evaluator self-floor (backslash probe,
true BC ~1e-14): 4.0/6.0/4.9e-9 rel at R1/R2/R3 — rung-independent, 20–25×
below the 1e-7 request, 200× below target. The certified evaluator is hereby
the standard post-solve column (decision_rules.md rewritten accordingly).

**Stage 4 DONE**: decision_rules.md — 2026-08-14 PROPOSED freeze block struck,
replaced by: certified-BC-metric definition (target 1e-6/rung) + evaluator
contract; per-rung tuning procedures (tune_fmm for the shared Krylov apply;
purpose-built FGS τ-ladder); ruling-3 shared-knob policy; ruling-4 dual-row
FGS reporting; certified per-run guard replacing bound-certification
restrictions (MAC ≥ 0.6 rejection withdrawn per Ryan). Agreement-stage
strength thresholds retired as historical (dCT ≤ 0.1% cross-check retained).

**Stage 6 SPEC drafted**: `phase_01_hpc_procedure.md` — per-rung R4–R7 recipe
(order: case/frozen-b → tune_fmm → R4-only direct spot-validation → fgstune →
fgsprecond → table), knobs that grow with N (ILU max_pattern_entries,
PressureLaplace itmax, FGS MAXIT/LEAF_SET), dense-config drop at R6–R7, and
the uncommitted-FastMultipole rsync + src md5 caution.

**Interpretation logged (Stage 2/3 tension)**: plan Stage 2 says "solver role
= the τ=1e-6 rung" while ruling 3 says ONE shared FGS knob set selected by
Stage 3's end-to-end winner. Ruling 3 (Ryan's words) controls: shared set =
Stage 3 winner, run to 1e-6 in the solver role; `..._table.jl` also emits an
`fgs_tuned1e6` comparison row whenever the τ=1e-6-tuned config differs, so
Ryan can override from data.

### 2026-08-15 night — Stages 2–3 R1–R2: τ-ladder tuned, preconditioner selected; big FGS speedups

**Stage 2 (fgstune) R1+R2 GATE PASS** (fgstune_{staircase,selected,verify}.csv;
34/35 candidates per rung, coordinate descent + coarse-τ seeds, MAXIT=40):

- R1 τ=1e-6 winner: p=6/MAC 0.3/leaf 150/**inner 5** — BC 8.7e-7 in 10 outer
  iterations, 1.39 s staircase (fresh verify 1.66 s / 12 it) vs 2.56 s at the
  Phase 0 production knobs (p=10/inner 20): **~1.8× faster**. Low p + few
  inner sweeps dominate once the certified BC guard replaces bound margins;
  p=4 stalls above 1e-6 (apply-accuracy-limited fixed point) but WINS the
  coarse rungs τ=1e-1/1e-4.
- R2 τ=1e-6 winner: p=8/MAC 0.4/leaf 100/inner 10 — 2.66 s staircase vs
  5.68 s prod (**~2.1×**). Verify PASS at 2.4e-7 but needed 15 it/6.0 s vs
  the staircase's 7 it — the known FGS multithread order-nondeterminism makes
  near-threshold stopping noisy at R2; Stage 5 timings carry this scatter
  (HPC k_reps ≥ 5 will bound it).
- Ladders are clean: crossing iteration grows ~1 per decade (R1 sweeps
  1→9, R2 1→6-9 across τ=1e-1→1e-6); every per-τ winner re-verified
  BC ≤ τ on a fresh cold solve.

**Stage 3 (fgsprecond) R1+R2** (fgsprecond.csv; FIXED-sweep linear applies,
shrink=true judgment, tuned outer apply):

- R1 winner τ=1e-6: sweeps=9 → FGMRES converges in 1 iteration, 2.87 s
  end-to-end (all six τ beat Phase 0's 17.2 s; plain FGS at 1.4 s remains
  the fastest config overall).
- R2 winner τ=1e-4: p=8/MAC 0.5/leaf 150/inner 5, sweeps=6, niter=2,
  9.81 s end-to-end vs Phase 0's 33.5 s (**3.4×**). Shared-set ≠
  τ=1e-6-tuned config at R2 → the Stage 5 table will emit the
  `fgs_tuned1e6` comparison row (ruling-3 interpretation, logged above).

R3 tuning in flight.

### 2026-08-15 night — STAGED (Ryan directive): FGS multithreading determinism + performance plan

Ryan (in session): stage work to examine FGS multithreading and improve
determinism + performance. Staged (NOT executed):
`fgs_determinism_performance_plan.md`. Code survey done during staging (so
the executing agent starts warm): solve.jl's GS loop is serial at the Julia
level (only BLAS inside leaf `mat \ rhs` threads); the in-loop fmm! passes an
empty direct list; assign_m2l!/make_direct_assignments! respect target
boundaries; the threaded interaction-list build merges in assignment-index
order — all deterministic-by-construction, so the observed run-to-run
variation (Phase 0: 267 vs 500 iters @4T, 1T reproducible; Stage 2 R2:
staircase 7 iters vs verify 15 in-process) is NOT explained by the obvious
suspects. Plan: Part A localization experiments (same-solver repeats,
ctor-vs-ctor bitwise diffs, thread matrix, M2M parent-boundary audit, BLAS
leaf-solve repeats) → Part B fixes gated on attribution (list-order
canonicalization / M2M partition fix / BLAS pinning) + an ungated
harness-side mitigation (margin-based stopping tol in the τ-ladder tuner) →
Part C performance (headline: cache leaf LUs — currently refactorized every
sweep; then deterministically-colored parallel sweeps + batched
update_nonself_influence!, both gated on Ryan as solver-variant changes).

### 2026-08-15 late night — Stage 2 R3 done (selection healthy); zero-margin verify FAIL confirms the margin fix; replay staged

R3 fgstune (old-script run, 35 candidates): τ=1e-6 winner p=6/MAC 0.3/
leaf 100/inner 5, 7.65 s staircase vs ~20 s at prod knobs (**2.6×**); ladder
sweeps 1→13. Zero-margin verification: **GATE FAIL at τ=1e-5** — fresh solve
stopped EARLY (9 it, BC 7.9e-5 > τ; staircase needed 12) because the
nondeterministic internal max-abs residual dipped under the zero-margin
tolerance prematurely; τ=1e-6 scattered the other way (27 vs 14 it, PASS at
2.0e-7). This is the exact failure mode the tuner edit (2026-08-15, during
the R3 run) addresses: verification stopping tolerance = geometric mean of
the crossing residual and its first LOWER successor (correctly directed —
stops at-or-after the crossing, never before), output renamed
fgstune_margin_verify.csv. The R1/R2/R3 runs all executed the pre-edit
script (their zero-margin rows live in fgstune_verify.csv, kept as
evidence).

Action (token-conscious replay): new `..._margin_verify.jl` re-runs ONLY the
per-τ verification solves (winners from fgstune_selected.csv, residual
staircases from fgstune_staircase.csv; MAXIT 60), writing
fgstune_margin_verify.csv — the 35-candidate pools are NOT re-run.
`..._table.jl` updated to derive the FGS production stopping tolerance by
the same margin rule (was: zero-margin crossing residual). The selection
CSVs are margin-independent (winners chosen on staircase cost-to-τ) and
stand as-is. Queue: fgsprecond R3 (in flight) → margin_verify R1/R2/R3 →
Stage 5 tables.

### 2026-08-15 ~23:30 — Margin verification GATE PASS R1–R3; Stage 2 closed; R3 staircase-extension mechanism added

fgstune_margin_verify.csv complete for all 18 (rung × τ) winners, all PASS.
Highlights: the R2 τ=1e-6 case that scattered 7→15 iterations under
zero-margin stopping now stops at 7 it / 2.74 s, matching its staircase —
the margin rule removes the threshold-stopping cost scatter. The R3 τ=1e-6
verify initially errored ("no finite residual below the crossing"): the
winner's instrumented solve had stopped AT the crossing iterate, so no
successor residual existed. (Ops correction: the two "killed" background
attempts were in fact this loud error — the stacktrace hadn't flushed to the
nohup log when I tailed it; foreground rerun revealed it. The
long-background-job kill hypothesis was wrong here.) Fix per the tuner
contract: `..._margin_verify.jl` now auto-EXTENDS a too-short staircase (one
instrumented cold solve of the winner knobs at 0.01·target·rms_b internal
tolerance, rows appended to fgstune_staircase.csv, same schema), then
re-derives the margin. R3 τ=1e-6 after extension: tol 8.77e-8, 15 it,
8.95 s, BC 3.27e-7 PASS. Production stopping tolerances for Stage 5 come
from these margin rows (`..._table.jl` mirrors the same margin_tol rule).

**Stage 3 summary R1–R3 (fgsprecond.csv)**: winners τ=1e-6 / 1e-4 / 1e-6;
end-to-end FGMRES 2.87 / 9.81 / 19.41 s vs Phase 0 fgmres_fgs
17.2 / 33.5 / 101 s (**6.0× / 3.4× / 5.2×**). R1+R3 winners coincide with
the τ=1e-6-tuned config (no comparison row); R2's differs (p8/.5/150/5
sweeps 6 vs 1e-6-tuned p8/.4/100/10) → table emits fgs_tuned1e6 row at R2.

### 2026-08-16 ~00:00 — Stage 5 R1–R2 tables done; shared-set ADMISSIBILITY constraint added (judgment, flag for Ryan)

solvetable.csv R1+R2 complete (all rows certified). R1: ldiv 0.019 s
(setup 7.5), gmres 94.4 s/59 it, jacobi 30.6/19, ilu 11.3/7, fgs 1.88/11 it
(single row — staircase crossing bc 8.7e-7 is within a half decade of
target, no ruling-4 dual row), fgmres_fgs 2.84 s/1 it.

**Finding + judgment (R2)**: the Stage 3 end-to-end winner (τ=1e-4 config
p8/MAC.5/leaf150/inner5) is INADMISSIBLE as the shared set: its staircase
STALLS at BC ≈ 3.57e-6 (apply-accuracy-limited fixed point — internal
residual falls to 5.6e-10 while true BC error floors above target), so it
can never serve the 1e-6 solver role. Ruling 3 demands ONE set for both
roles ⇒ `stage3_winner()` in `..._table.jl` now restricts to configs whose
staircase crosses 1e-6, picking the fastest admissible end-to-end row and
logging skips loudly. R2 shared set ⇒ τ=1e-5 winner p8/MAC.4/leaf100/
inner5, sweeps 9 (11.05 s end-to-end, +1.24 s vs the inadmissible pick).
R2 fgs rows: shared set 3.30 s/13 it (bc 1.6e-7); fgs_tuned1e6 comparison
2.67 s/7 it (bc 6.0e-7) — the two-set alternative is 0.6 s faster in the
solver role; RYAN CALL whether one-set simplicity is worth it. fgmres_fgs
11.15 s/2 it. (Ops note: nohup logs do NOT get julia stacktraces flushed on
error — silent log ends were errors, not kills; foreground reruns reveal
them. Applied to R2/R3 table chunking.)

### 2026-08-16 — Stage 5 COMPLETE: tuned solver × rung table R1–R3 (solvetable.csv); plan stages 1–6 all delivered

All rows certified (`bc_certified=true`, BC ≤ 1e-6 everywhere; per-solve
seconds, local 4-thread, k_reps=1, setup separate; R3 chunked per config
group):

| config (row_kind) | R1 t_solve | R2 | R3 | R1 setup | R2 | R3 |
| --- | --- | --- | --- | --- | --- | --- |
| backslash_ldiv | 0.019 | 0.077 | 0.226 | 7.5 | 34.1 | 144.1 |
| krylov_gmres | 94.4 | 226.9 | 651.3 | 0.2 | 0.2 | 0.3 |
| krylov_jacobi | 30.6 | 49.3 | 154.9 | 6.9 | 24.8 | 77.2 |
| krylov_ilu | 11.3 | 35.1 | 106.0 | 3.6 | 6.6 | 14.9 |
| fgs (shared set) | 1.88 | 3.30 | 8.64 | 10.4 | 21.5 | 65.0 |
| fgs_tuned1e6 (comparison, R2 only) | — | 2.67 | — | — | 18.2 | — |
| fgmres_fgs (shared set) | 2.84 | 11.1 | 17.1 | 9.4 | 18.1 | 62.0 |

Shared FGS sets: R1 τ=1e-6 p6/MAC.3/leaf150/inner5 sweeps9; R2 τ=1e-5
p8/.4/100/5 sweeps9 (admissibility-filtered, see previous entry); R3 τ=1e-6
p6/.3/100/5 sweeps13. No ruling-4 dual rows triggered (no staircase
crossing overshot by more than a half decade). FGS-family speedups vs
Phase 0 configs: solver ~1.4–2.3×, preconditioned ~3–6×. krylov rows
reproduce solvetime.csv within noise (gmres 651.3 vs 651.3; jacobi
154.9 vs 150.8; ilu 106.0 vs 103.5) — the earlier baseline stands.

**Phase 1 tuning plan (2026-08-15): ALL SIX STAGES DELIVERED.** Remaining
Phase 1 work = HPC R4–R7 per `phase_01_hpc_procedure.md` (fresh agent) +
Ryan reviews: (a) one-set-vs-two FGS tradeoff at R2 (0.6 s cost of ruling-3
simplicity); (b) the admissibility judgment; (c) the staged FGS
determinism/performance plan (`fgs_determinism_performance_plan.md`).

### 2026-08-17 — FGS multithread determinism FIXED (Ryan, FastMultipole): M2L owner overlap via non-owner-major list

Ryan root-caused and fixed the FGS multithread nondeterminism: the
interaction list was deterministic but not owner-major, so threaded M2L
assignment could overlap target owners; FastMultipole now canonicalizes
(target, source) order (stable counting sorts) before threaded assignment
and dependent matrix/index construction. Verification (Ryan,
`results/fgs_determinism/summary.md`): full four-way {Julia 1T/4T} ×
{BLAS 1T/4T} matrix bit-exact within every configuration; timing +0.87%
(passes the ≤5% gate); all 18 R1–R3 post-crossing margins re-pass certified
BC. Plan doc marked COMPLETE (Parts B1+B4 implemented). Consequences:
τ-ladder stopping tolerances now transfer between runs; FGS timing scatter
eliminated locally; the margin rule stays (still correct, now un-stressed).
FastMultipole dev checkout now carries FOUR uncommitted campaign fixes.

### 2026-08-17 — Ryan rulings implemented: FGS set = per-rung 1e-6-tuned config (both roles); post-determinism-fix FGS rows refreshed

Rulings (Ryan in session): (a) the shared FGS knob set per rung = the
τ=1e-6-tuned config for BOTH roles (supersedes end-to-end selection of the
knob set; only the preconditioner SWEEP COUNT stays end-to-end-selected,
now over the 1e-6 config's own staircase crossings); (b) the admissibility
filter is confirmed vestigial under (a) — kept as a tripwire, rationale
(apply-accuracy BC floor of coarse-tuned configs) recorded in
decision_rules.md. Implementation: fgsprecond.jl `SWEEP_LADDER_1E6=1` mode
(sweep ladder of the 1e-6 config); table.jl `stage3_winner` pins knobs to
the 1e-6 config.

Sweep ladders R1–R3 (fgsprecond.csv, post-determinism-fix): winners
R1 sweeps 9 (2.15 s end-to-end), R2 sweeps 6 (4.57 s), R3 sweeps 13
(14.64 s). Refreshed solvetable.csv FGS-family rows (all certified;
BC values reproduce BIT-FOR-BIT pre/post fix — trajectories unchanged,
only threading order canonicalized):

| row | R1 | R2 | R3 | pre-fix R1/R2/R3 |
| --- | --- | --- | --- | --- |
| fgs t_solve (s) | 0.98 | 1.52 | 4.30 | 1.88 / 3.30 / 8.64 |
| fgmres_fgs t_solve (s) | 2.06 | 4.75 | 13.05 | 2.84 / 11.1 / 17.1 |

The ~2× FGS wall-time drop is the owner-major M2L canonicalization's
locality/balance side effect plus (R2) the ruling's knob switch — far
beyond the +0.87% fixed-ten delta Ryan measured on the isolated apply;
worth confirming on HPC. FGS shared sets of record: R1 p6/MAC.3/leaf150/
inner5 sw9; R2 p8/.4/100/10 sw6; R3 p6/.3/100/5 sw13. Superseded rows
(pre-fix fgs/fgmres/fgs_tuned1e6) remain in solvetable.csv — filter by
date/commit; latest rows are of record.

### 2026-08-17 — Attribution correction + τ-ladder re-tune launched (Ryan: leaf-LU caching landed)

Ryan implemented FastMultipole leaf-LU caching (Part C1) — it was in the
dirty working tree during today's refreshed runs, which CORRECTS the
previous entry's speculation: the ~2× FGS wall-time drop is the LU cache
(leaf solves O(L³)→O(L²) per sweep; bit-neutral, hence the bit-identical
BC values; outside the isolated-apply +0.87% measurement), not an
owner-major locality effect. Consequence (pre-registered in the plan):
the 08-15 knob tuning is STALE — leaf/inner tradeoffs changed — so the
τ-ladders R1–R3 are being re-run (fgstune, then margin verify with replay
fallback, then SWEEP_LADDER_1E6 sweep ladders, then table refresh). All
mechanical scripts, no re-design. FastMultipole now carries FIVE
uncommitted fixes — commit hygiene before HPC rsync remains flagged.

### 2026-08-18 — Campaign moved to HPC (Ryan: local machine too loaded for honest benchmarks); R1–R4 chains submitted

Local post-LU-cache re-tune KILLED per Ryan (R1 had finished GATE PASS —
τ=1e-6 verify 0.55 s / 7 it with the cache — but rows are load-contaminated;
HPC supersedes ALL local timing rows; local CSVs remain as pipeline-shakeout
evidence, distinguishable by hardware_tag). Sync verified: FastMultipole
branch flowpanel-20260817 @ a9b734a BOTH ends (Ryan committed the LU cache
after the md5 check caught it uncommitted; src md5 MATCH), FLOWPanel src
rsynced + md5 MATCH (LC_ALL=C — locale sort produced a false mismatch
first), ILUZero resolved into the cluster Manifest, R4–R7 capped meshes
rsynced. Infrastructure landed: phase1_case.jl (LADDER R4–R7; TUNED
csv-first w/ hardcoded R1–R3 fallback; CACHE_B frozen-b disk cache;
KNOBS_MODE for single-mode jobs), dense-free tune stage
(rotor_hover_solver_phase1_tune_hpc.jl: FGS state + certified check +
tune_fmm), Slurm launchers benchmark/slurm/p1_{tune,table}.sh (CONFIGS
colon-separated — commas collide with sbatch --export).

Submitted (2026-08-18, per-rung afterok chains tune → table-multi(64T,
k_reps 5, all six configs) → table-single(1T)): R1 13193383/84/85,
R2 13193386/87/88, R3 13193389/90/91, R4 13193392/93/94. R5–R7 HELD until
R4's tune validates the pipeline at scale (then same chains; R6/R7 drop
backslash_ldiv per ledger). Monitoring at ≥5-min cadence (BYU rule).

### 2026-08-18 — Hardware pinned (Ryan consistency question exposed two gaps); R1–R3 mesh miss fixed; all chains resubmitted

Ryan asked which machine the HPC benchmarks request → audit found the first
submission specified NO constraint and NO --exclusive (jobs landed on m12
zen3 nodes only by luck, and SHARED). Also the R1–R3 tune jobs FAILED at
~5 min: the R1–R3 capped ladder meshes were never rsynced (only R4–R7
were; the cluster's 9 pre-existing captess4 files are other DJI meshes).
Fixes: meshes rsynced; ALL p1 jobs now submitted with
`--constraint=zen3 --exclusive` (m12-class: 2×64-core AMD Zen3, 512 GB —
the campaign's pinned hardware; ruling 5's dedicated-node requirement now
enforced, incl. tune jobs since knob selection is timing-based); the
originally-running shared-node R4 tune killed rather than trusted. Jobs of
record: R1 13193401/02/03, R2 13193404/05/06, R3 13193407/08/09,
R4 13193412/13/14 (tune/table-multi/table-single per rung). hardware_tag
column records the node; banner shows julia 1.12.6, fm a9b734a. R5–R7
still held on R4 tune validation.

### 2026-08-18 — Ryan cost rulings implemented: direct column retired, adaptive k, gmres 1T drop-out at R6–R7, coarse-τ seeds deleted

Re-estimate showed R7 table-single ≈ 33–35 h even without the direct
column (1T gmres ~12 h/solve dominating) → per Ryan's conditional, the
full package landed: (1) per-row direct O(N²) residual RETIRED
(DIRECT_RESID=1 re-enables; certified BC column is the metric of record);
(2) adaptive min-of-k in table.jl (5/<60 s, 3/<10 min, 2 above; k recorded
per row) + cost-ceiling drop-out: krylov_gmres out of SINGLE mode at
R6–R7; (3) coarse-τ seed candidates deleted from fgstune (vestigial under
the 08-17 sweep-ladder ruling); (4) hierarchical knob extrapolation +
warm-started pattern search HELD. decision_rules.md timing protocol
amended. Revised projections: R5 chain ~1 day, R6 ~1 day, R7 ~1–1.5 days
(tuner becomes R7's long pole at ~10–16 h). Protocol note: pending R1–R4
TABLE jobs pick up the new table.jl at start (uniform protocol across all
rungs); the RUNNING R1–R4 tune jobs still carry the coarse-τ seeds
(harmless — extra candidates only).

### 2026-08-18 — NFS concurrent-append clobbering found + fixed (per-rung dirs); all chains resubmitted; first Zen3 calibration points

R1's table-multi FAILED with the "impossible" tripwire → root cause: the
four rung pipelines run concurrently on different nodes APPENDING TO THE
SAME CSVs over NFS; non-atomic appends clobber (R1's 32-candidate
staircase reduced to 17 foreign rows). Local runs never hit it
(sequential). Fix: `PER_RUNG_DIR=1` in the Slurm launchers →
results/phase1/<mode>/<rung>/ isolates every writer (phase1_case.jl);
polluted CSVs quarantined on the cluster
(multi/quarantine_nfs_clobber_20260818/). All chains resubmitted:
R1 13193470/71/72, R2 13193473/74/75, R3 13193476/77/78,
R4 13193479/80/81 (zen3+exclusive, adaptive-k table protocol, seedless
tuner).

Zen3 calibration from the clobbered-but-informative first round: R1 tune
16 min wall; table-multi krylov rows (k=5) gmres 34.0 s / jacobi 11.1 s /
ilu 4.08 s (vs 94/31/11 s at 4T local); **R1 1e-6 FGS winner on Zen3 =
p14/MAC 0.6/leaf 150/inner 10** (vs p6/MAC 0.3 on the Mac) — hardware-
dependent tuning confirmed in the strongest way; per-rung/per-hardware
re-tuning (Ryan ruling) was the right call.

### 2026-08-18 — Round-2 verdicts: R1/R3/R4 clean; R2 tune died on a transient depot error; R5–R7 submitted (single-mode tables split serially)

Chain outcomes (13193470–81): **R1, R3 chains COMPLETE; R4 tune + table-multi
COMPLETE, table-single running** (13193481, backslash row landed). Every row
of every finished solvetable.csv is `bc_certified=true` (R1 7 rows incl. the
FGS dual row, R3/R4 6 rows). Headline multi-mode (64T) solve times, R4
(58 192 panels): backslash ldiv 1.3 s / gmres 207 s / jacobi 54 s / ilu 44 s /
**fgs 13.6 s** / fgmres_fgs 23.7 s — FGS fastest iterative at every rung both
modes; ILU best true-Krylov. Single-mode R3: gmres 5 033 s vs fgs 7.9 s.

**R2 tune 13193473 FAILED** — NOT the tuner: the τ-ladder finished with
verification GATE PASS (fgstune CSVs complete on disk); the chain died when
the *fgsprecond* stage's fresh Julia process failed to load PythonPlot:
`InitError: could not find source path for package MbedTLS_jll` (transient
depot/precompile state on m12-3-9 at 06:50; R1 before and R3/R4 after loaded
the identical stack cleanly). afterok cancelled R2's table jobs. Recovery:
new `benchmark/slurm/p1_tune_s3.sh` (stage-3-only launcher — the stage
scripts APPEND rather than checkpoint-skip, so rerunning the full tune would
duplicate tune.csv/fgstune rows) → chain 13206067 (s3) → 13206068
(table-multi) → 13206069 (table-single), all six configs.

**R5–R7 submitted** per the drop-out recipe (R5 all six; R6/R7 drop
backslash_ldiv both modes + krylov_gmres single mode): R5 13206072/73/74/75,
R6 13206076/77/78/79, R7 13206080/81/82/83, all zen3+exclusive, PER_RUNG_DIR.
**Ops judgment (logged): each rung's table-single is SPLIT into two SERIAL
afterok jobs** (bulk configs, then the slowest config alone — R5 gmres;
R6/R7 jacobi): extrapolating measured 1T costs (~3×/rung; gmres R3-single
84 min/solve) puts a monolithic R5/R6/R7 single job near or past the 3-day
walltime; the split keeps each job inside 3-00:00:00 and a timeout loses only
one config's rows. Serial chaining preserves the NFS single-writer rule.
Mem/time per job sized from round-2 measurements (R5 dense G ≈ 94 GB → 256G
nodes for R5 tables). Residual risk flagged: R7-single-jacobi extrapolates to
~69 h ≈ its 3-day limit — may time out; rerunnable alone if so.
