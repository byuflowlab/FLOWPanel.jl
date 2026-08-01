# Phase 4 — Testing, Verification, and Validation

**Status:** IN PROGRESS — V0 complete; V1 still blocked at the wing
formulation-sanity gate, which was rebuilt as a scripted experiment on
2026-07-30. The gate now fails *informatively*: the leading suspect (open-tip
bodies) is eliminated and the gap is split into a reconstruction part and an
unsteady-circulation part.

**Prerequisite:** approved Phase 3 implementation (approved 2026-07-29)

**Phase approval:** [ ]

## Prior-phase handoff

Phase 1 approved a four-combination attribution matrix: Route A retains rigid
`Das`, Route B is TE-anchored, and each route independently supports jump or
pressure closure. Route B testing begins with B1 Euler timestep refinement,
then B2 persistent subpanels after an equal-strength cancellation null gate;
Katz–Plotkin \(k=0.2,0.25,0.3\) remains a comparator. The pressure residual is
paired panel-centroid pressure continuity, so sharp-edge convergence must be
demonstrated rather than assumed.

Phase 2 must define the approved interfaces and diagnostics. Phase 3 must
provide the approved implementation and focused mechanical-test record. Read
the hub and the approved [Phase 3 implementation](phase_03_implementation.md)
before beginning validation.

## Scope

Run unit/default regressions, pressure-provider and failure-path tests,
cross-solver and backend comparisons, mesh/tolerance/initial-guess studies,
the attached-wake length/`Das` sweep, and the planned symmetric-wing,
finite/swept-wing, and rotor-hover validations. Use at most six local threads.

The validation must distinguish whether pressure closure and/or Route B
removes, reduces, or merely relocates wake-geometry sensitivity. Route A must
never be described as `Das`-free, and Route B conclusions require
timestep/wake-resolution convergence.

## Governing decisions

- Evaluate all four wake-attachment/Kutta-closure combinations independently.
- Treat default preservation and deterministic failure behavior as gates.
- Require both B1 temporal convergence and the B2 cancellation/refinement
  sequence before drawing Route B accuracy conclusions.
- Report negative results and relocated sensitivities as valid outcomes.

## Deliverables

- Reproducible verification commands, configurations, data, and results.
- Default-regression and solver/backend/provider compatibility evidence.
- Nonlinear robustness studies over tolerance and initial guess.
- Four-combination attribution results.
- B1 timestep-refinement/Richardson results; B2 equal-strength null and
  \(N_s=2,4,8\) interpolation results; B3 comparator results if implemented.
- Symmetric-wing, finite/swept-wing, and rotor-hover validation conclusions.

## Acceptance gate

Phase 4 and Item 015 are complete only after all required evidence satisfies
the approved criteria and the user explicitly approves the conclusions. Only
then may the item-level approval columns in `BRAINSTORM/INDEX.md` be checked.

## Campaign plan (staged 2026-07-29, from the session-approved plan)

- **V0 — regression & mechanical gates** (local): unit regression + kutta
  suites; default-preservation and failure-path evidence (cite testsets).
- **V1 — four-combination attribution matrix** (core deliverable): diamond
  fixture record → symmetric wing (item-014 corrected SSW mesh,
  `assert_ssw_mesh_symmetry`; zero-lift symmetry then finite-α) → finite/swept
  wing steady. Metrics: bound-Γ distribution, CL, paired TE Δp, accepted `c`,
  iteration counts. A/jump vs A/pressure isolates the closure; A/jump vs
  B/jump isolates the geometry; B/pressure tests composition.
- **V2 — sensitivity relocation** (the item's reason to exist): Route A `Das`
  sweep under BOTH closures (η-ladder machinery from
  `examples/ssw_eta_convergence.jl`); sharp-edge convergence for pressure
  closure (chordwise/TE refinement + kernel-offset study — Phase 1 §4
  mandates this).
- **V3 — Route B temporal convergence**: B1 unsteady SSW, ≥3 halved
  timesteps, observed order + Richardson vs the finite-AR Wagner oracle
  (Jones NACA-681); B2 equal-strength null gate then N_s=2,4,8 (NOTE: B2
  blocks not constructible in Phase 3 — becomes a Phase 4 implementation
  ticket or documented deferral); B3 Katz–Plotkin k=0.2/0.25/0.3 comparator
  if cheap.
- **V4 — nonlinear robustness**: tolerance sweep, initial-guess perturbation,
  Broyden-vs-fallback on a physical case, density invariance of accepted c.
- **V5 — rotor hover** (HPC, mission tie-in): A/pressure (and B if validated)
  on DJI9443 at production settings vs the η ladder (floor-free η=0.2
  0.0522 … η=1.0 0.0713) and the 0.068–0.072 band.

Constraints: ≤6 local threads; >20 min runs go to HPC per
`agent_policies/HPC.md` (user submits); negative/relocated-sensitivity
results are valid outcomes; Route A is never described as `Das`-free.

## Phase 4 progress and approval log

- **2026-07-29 — Phase seeded.** Blocked pending explicit Phase 3 approval; no
  validation has started.
- **2026-07-29 — Phase 3 approved; V0 complete.** Full legacy unit regression
  (fmm, kernel_gradient, body, solver, liftingbody, wake, postprocess,
  simulate, added_mass, warmstart, replay) green single-threaded, plus the
  hardened kutta suites (658 + 62 assertions). Default-preservation evidence:
  testset "default pair is the legacy path" (bitwise `==`);
  rejection-before-mutation evidence: "support-boundary validation" (settled
  body/wake pair unchanged through a rejected call); failure-path
  determinism: "failure rollback (:error)", "failure rollback at a settled
  step (:error)", "commit failure is deterministic (:commit)", "explicit jump
  fallback (:jump)" (bitwise `==` vs legacy), "snapshot restore round trip".
  All four combinations converge on the diamond fixture with genuinely
  nonzero corrections ("pressure closure changes the solution": max|c| ≈
  8.6e-3·C_s, A/pressure moves strengths ≈2.2% vs A/jump, independent Δp
  recompute within 2× the gate while the jump baseline violates continuity by
  ≫10× the gate). Next: V1 wing-level attribution.

## 2026-07-29 — V1 attempted; diamond complete, wing matrix blocked

> **SUPERSEDED 2026-07-30.** The 2026-07-30 clear-context review found this
> section's diagnosis unreproducible: of the numbers quoted below, only the
> \(t^*=2\) gate row was produced by a script and written to CSV. The
> \(t^*=8\) row, the semi-infinite row, both Kutta–Joukowski CLs, and the
> paired-TE Γ ranges were collected by hand and are not regenerable, and they
> are mutually inconsistent (Γ agreeing to ~3% while KJ CL differs 2.5×, and
> Neumann's own KJ and pressure CLs differing 3.3×). **The localization claim
> is withdrawn**, and the diamond table's metrics have been renormalized. Read
> the 2026-07-30 section for the current record; the material below is retained
> only as the history of a blocked attempt.

### Reproducible driver and configuration changes

Added `examples/kutta_v1_attribution.jl`, selected with
`KUTTAV1_STAGE=diamond|ssw0|ssw_alpha|sweptwing` and writing restartable CSV
records under `data/kutta_v1_attribution/`. The driver refuses more than six
Julia threads. It records summary metrics, full SSW CL histories, paired-TE
bound circulation, and swept-wing sectional loading when the prerequisites
permit those stages.

`SSWConfig` now has `bodytype=:neumann|:dirichlet`; the Dirichlet option is
`RigidWakeBody{Union{ConstantSource,ConstantDoublet}}`. The no-shedding probe
and final body use the same selected type, and shedding remains derived from
the constructed probe body's rewound cells. `bodytype` is included in the
resume tag. The classic swept-wing example now selects the same Dirichlet
source+doublet type. These defaults preserve the legacy SSW Neumann case.

Commands used:

```bash
JULIA_NUM_THREADS=2 KUTTAV1_STAGE=diamond KUTTAV1_FORCE=true \
    julia --project=. examples/kutta_v1_attribution.jl
JULIA_NUM_THREADS=2 KUTTAV1_STAGE=ssw0 KUTTAV1_FORCE=true \
    julia --project=. examples/kutta_v1_attribution.jl
```

### Diamond four-combination record

Source: `data/kutta_v1_attribution/diamond_summary.csv` and
`diamond_spanwise_gamma.csv`.

> **SUPERSEDED 2026-07-30 (metrics).** The \(\widehat{\Delta p}\) columns below
> were each divided by *that cell's own* frozen kutta `pressure_scale`, which is
> ~116 for A/pressure and ~122 for B/pressure but is undefined (NaN, so the code
> substituted the freestream \(q\)) for the two jump cells. Four cells were
> therefore normalized three different ways and the "pre" and "post" columns of
> the same row by different scales, so no two numbers in the table were
> comparable. See the renormalized table in the 2026-07-30 section. The mean
> iteration counts below also averaged the deterministic `startup_jump` step in.

| cell | fixture CL | max paired \(|\widehat{\Delta p}|\), jump/pre | max paired \(|\widehat{\Delta p}|\), post | max \(|c|/C_s\) | mean outer iterations | body solves | disposition |
|---|---:|---:|---:|---:|---:|---:|---|
| A/jump | 33.1714 | 4.3362e-2 | 4.3362e-2 | 0 | 0 | 0 | jump |
| A/pressure | 33.8502 | 4.3362e-2 | 1.7474e-7 | 8.8513e-3 | 2.0 | 16 | converged |
| B/jump | 30.9702 | 1.8133e-1 | 1.8133e-1 | 0 | 0 | 0 | jump |
| B/pressure | 38.5168 | 1.8133e-1 | 1.1534e-7 | 2.4455e-1 | 4.5 | 26 | startup_jump, then converged |

"Fixture CL" is **not** a lift coefficient: the diamond `ForceMonitor` runs with
`NoNormalization` and the driver divides only by \(q\), never by a reference
area, so the column is area-dimensional lift \(L/q\) on an ad-hoc fixture. It is
comparable across the four cells and meaningless in absolute terms.

This formalizes the V0 smoke. Pressure closure is not a no-op and suppresses
the independent paired-pressure residual on both geometry routes. Route B's
first `startup_jump` is the designed cold-start disposition.

### Hard-prerequisite failure and diagnosis

The required AR=6, alpha=5-degree, corrected-mesh A/jump comparison failed
before the zero-lift or finite-alpha matrix was run:

| comparison | Neumann CL | Dirichlet CL | relative gap |
|---|---:|---:|---:|
| settled near \(t^*=2\) | 0.356613 | 0.273324 | 23.36% |
| settled near \(t^*=8\) | 0.411381 | 0.302816 | 26.39% |
| semi-infinite steady | 0.421219 | 0.359815 | 14.58% |

The mismatch is not premature transient sampling. In the steady diagnostic,
the paired-TE circulation ranges are much closer:
Neumann \([-0.23091,-0.12384]\), Dirichlet
\([-0.22391,-0.12004]\), about a 3.0% change at both ends (3.03% and 3.07%),
while pressure-integrated CL differs by 14.6%. Kutta--Joukowski force is also
formulation-dependent on this coarse surface (Neumann CL 0.12740, Dirichlet
0.31323).

> **WITHDRAWN 2026-07-30.** These numbers are not reproducible and are mutually
> inconsistent: the two KJ CLs differ by 2.5× while the circulations behind them
> differ by 3%, and Neumann's own KJ CL (0.12740) and pressure CL (0.42122)
> differ by 3.3×, which no consistent extractor permits. The scripted
> Kutta–Joukowski extractor added on 2026-07-30 reproduces the Neumann steady
> pressure CL to 5.2% and the lifting-line anchor to 2.9%, and it contradicts
> both hand-collected KJ values. The localization sentence that followed —
> attributing the gap to global force reconstruction rather than circulation —
> is therefore **withdrawn, not merely softened**: the diagnostics it rested on
> were mutually inconsistent and localization was never established. The
> 2026-07-30 section replaces it with a scripted split.

Per the pre-registered hard gate, `ssw0` stopped with an error and V1 did
**not** run or interpret the SSW zero-lift matrix, finite-alpha matrix, or
swept-wing comparison. V1 remains blocked; this is a negative gate result, not
evidence against either Route B or pressure continuity.

### Regression record

- `runtests_unit_kutta.jl`: 658/658 pass.
- `runtests_unit_kutta_routeb.jl`: 62/62 pass.
- `runtests_unit_simulate.jl`: not green in the existing dirty item-014 wake
  worktree. It errors in `FLOWVPM.Estr_direct_multithreaded` for a zero-particle
  field (`ArgumentError: step cannot be zero`) in the backend-split test. The V1
  change touches no `src/` file and does not enter that particle-wake path; the
  failure is retained as an unresolved verification limitation.

Next action is to resolve or explicitly replace the Dirichlet--Neumann CL
sanity oracle. Until then, V2--V5 and all approval boxes remain gated.

> **Note added 2026-07-30 (V0/V1 apparent contradiction, resolved).** The V0
> entry above reports the full legacy unit regression green, including
> `runtests_unit_simulate.jl`; this V1 entry reports it failing. Both are
> correct: the difference is **thread count**, not any V1 change. See the
> reconciliation in the 2026-07-30 section.

## 2026-07-30 — V1 review remediation: gate rebuilt as a scripted experiment

Plan: `plans/20260730_kutta_v1_review_remediation.md`. This section supersedes
the 2026-07-29 gate diagnosis in full. Every number quoted below is written by
`examples/kutta_v1_attribution.jl` into `data/kutta_v1_attribution/`, and the
whole battery reruns in about three minutes on six local threads; the run was
repeated and reproduced bitwise.

### What the review found, and what changed

1. **Unreproducible, self-contradicting diagnosis.** Only the \(t^*=2\) gate row
   had a script behind it. The rest was hand-collected and internally
   inconsistent. The localization claim is withdrawn (see the annotations in the
   2026-07-29 section) and every gate number is now scripted.
2. **The gate sampled a transient and measured only CL.** It compared at
   \(t^*\le 2\) and recorded no strength-level comparison, so it could not
   separate the solve from the velocity reconstruction. The gate now settles
   over \(t^*\in[7,8]\) (9 samples; peak-to-peak spread 0.09–0.17% of CL, versus
   a curve still climbing 15% between \(t^*=2\) and \(t^*=8\)) and records
   spanwise TE \(\Delta\)strength, a semi-infinite steady solve, and a
   Kutta–Joukowski circulation diagnostic per cell.
3. **Open-tip Dirichlet bodies (the decisive new fact).**
   `suddenly_started_wing.jl` builds a "Triangular open-tip extrusion" with
   `watertight=false`, and `simplewing` lofts sections with no caps either. The
   interior-Dirichlet Green identity assumes a **closed** surface; the solver
   warns for the Neumann analogue of this mistake
   (`src/FLOWPanel_solver.jl:261`) but is silent for non-watertight Dirichlet.
   Flat centroid-fan tip caps are now the Dirichlet configuration, via a shared
   `add_flat_tip_caps` in `examples/helper_functions.jl`, and uncapped Dirichlet
   is retained as a gate row so the caps' effect is *measured*.
4. **`examples/sweptwing.jl` had been silently flipped to Dirichlet.** Its
   default is back to the legacy `RigidWakeBody{VortexRing,1,Float64,false}`,
   the Dirichlet variant is opt-in via `SWEPTWING_BODYTYPE=dirichlet` (which
   also enables flat caps), and a formulation token now appears in `grid_tag`,
   `run_name`, `save_path`, and the mirrored-Cp cache load paths, so the tracked
   Neumann-era artifacts can never mix with Dirichlet output.
5. **Diamond metric blemishes and doc/log gaps** — addressed below and in the
   annotations above.

### Cap utility and mesh checks

`add_flat_tip_caps(nodes, cells)` chains the singly-referenced edges into
boundary loops, adds one centroid node per loop, and fans it to every loop edge,
traversing each boundary edge in the direction opposite to the surface cell so
the winding stays consistent. Nodes and cells are **appended**, so existing
indices — and any shedding matrix expressed in them — survive unchanged. The
centroid fan inherits the contour's symmetry: the centroid maps to the centroid
and each fan face to the fan face of the image edge.

`SSWConfig` gains `caps::Symbol = :auto`, resolved by `ssw_caps` to `:flat` for
`bodytype=:dirichlet` and `:none` for `:neumann` — the only pairing consistent
with both formulations, since a watertight Neumann `RigidWakeBody` is
rank-deficient. `:auto` rather than a literal `:none` default is what lets an
explicit `caps=:none` on a Dirichlet body be distinguished from the default,
which the gate's third row requires. The resolved value is a `_ssw_case_tag`
token, so capped and uncapped histories can never resume across each other.

Checks at AR=6, `n_span=12`, `n_airfoil=21` (all passing):

| check | result |
|---|---|
| uncapped mesh watertight | `false` (unchanged) |
| capped mesh watertight | `true` |
| nodes / cells, uncapped → capped | 260 → 262, 480 → 520 (= +2·20 fan faces) |
| appended-only (prefix bitwise equal) | yes, nodes and cells |
| `assert_ssw_mesh_symmetry`, capped | passes (both reflections, incl. cap faces) |
| shedding pairing, uncapped vs capped | bitwise identical, 12 pairs |
| Neumann legacy tag and mesh | unchanged (`caps=:none` resolves for `:neumann`) |
| spanwise \(\Gamma\) symmetry, all six gate cells | 0 exactly |

### Gate battery: six cells, pre-registered criterion

Three formulation/cap configurations × two `grad_mu` reconstructions, all at
A/jump (no attachment or closure runtime involved), AR=6, \(\alpha=5°\),
`n_span=12`, `n_airfoil=21`, `dt*=0.125`, direct backend, `t_end_star=8`.
`:tri` is the legacy `SSW_GRAD_MU_OPTIONS` setting and the repo's known-divergent
configuration on an all-triangle mesh; `:tri_robust` is `basis=:tri` with
`tri_robust=true`. Quad bases are N/A — the mesh has no quads to recover.

Pre-registered before the runs: pass requires the settled-CL gap
\(|\mathrm{Dir_{capped}} - \mathrm{Neu}|/\mathrm{Neu} \le 5\%\) **and** spanwise
TE \(\Delta\)strength agreement (no sign flips, shape correlation \(\ge 0.98\)).

Source: `gate_summary.csv`.

| cell | caps | basis | settled CL, \(t^*\in[7,8]\) | steady semi-inf. CL | \(C_L^{KJ}\) unsteady | \(C_L^{KJ}\) steady | cells |
|---|---|---|---:|---:|---:|---:|---:|
| Neumann | none | `:tri` | 0.411270 | 0.421219 | −0.389578 | −0.399434 | 480 |
| Neumann | none | `:tri_robust` | 0.367931 | 0.377169 | −0.389578 | −0.399434 | 480 |
| Dirichlet | flat | `:tri` | 0.306989 | 0.364783 | −0.339535 | −0.387432 | 520 |
| Dirichlet | flat | `:tri_robust` | 0.301138 | 0.358840 | −0.339535 | −0.387432 | 520 |
| Dirichlet | none | `:tri` | 0.302757 | 0.359815 | −0.339542 | −0.387462 | 480 |
| Dirichlet | none | `:tri_robust` | 0.297236 | 0.354244 | −0.339542 | −0.387462 | 480 |

The Kutta–Joukowski column is \(C_L = 2\int\Gamma\,dy/(U S)\) over the paired TE
strength jump (\(\Delta\gamma\) for the vortex-ring body, \(\Delta\mu\) for the
source+doublet body), with Voronoi spanwise weights — exact for this uniform
mesh. It is identical across the two bases by construction, which is itself a
check. **Extractor validation** (`gate_kj_validation.csv`), run before KJ is
admitted as evidence: on the Neumann steady case it gives −0.399434 against that
case's own pressure-integrated CL of 0.421219 — a 5.2% magnitude error, inside
the pre-registered 10% tolerance — and against the lifting-line anchor
\(2\pi\alpha\,AR/(AR+2) = 0.411234\), a 2.9% error. KJ is therefore admitted.
This also resolves the 2026-07-29 internal contradiction: the hand-collected KJ
values there (Neumann 0.12740, Dirichlet 0.31323) are both wrong, and the two
formulations' KJ CLs do **not** differ by 2.5×.

### Gate verdict: FAIL, and the caps are not the cause

`ssw_sanity_gate.csv`:

| quantity | value |
|---|---:|
| Neumann settled CL (`:tri`) | 0.411270 |
| Dirichlet **flat-capped** settled CL | 0.306989 |
| Dirichlet **uncapped** settled CL | 0.302757 |
| relative gap, capped | **25.36%** (criterion ≤ 5%) |
| relative gap, uncapped | 26.38% |
| fraction of the gap closed by caps | **3.90%** |
| spanwise sign flips | 0 |
| spanwise shape correlation | 0.999815 |
| spanwise least-squares level ratio Dir/Neu | 0.870195 |
| spanwise max relative difference | 14.02% |
| max basis sensitivity across cells | 10.54% |
| **pass** | **false** |

The plan's leading suspect is **eliminated**: closing the tips moves the
Dirichlet settled CL from 0.302757 to 0.306989, closing 3.9% of a 26% gap. The
open-tip Green-identity violation was real and worth fixing on its own terms —
the Dirichlet body is now watertight and inside its formulation's assumptions —
but it is not what the gate was detecting. The user's reported observation that
flat caps match Neumann-no-caps best is reproduced in *direction* (capped is
closer than uncapped) and is negligible in *magnitude* at this resolution.

The spanwise comparison is the other half of the story: the two formulations'
TE circulation distributions have the same shape to a correlation of 0.9998 with
no sign flips, and differ by a nearly uniform level factor of 0.870. The
disagreement is a circulation **level**, not a spanwise redistribution.

### Where the gap actually lives

`gate_localization.csv` splits each gap two ways. The KJ CL depends only on the
TE strength jump, so it is blind to `grad_mu` and to pressure integration; the
`ForceMonitor` CL depends on both.

| comparison | basis | unsteady settled gap | unsteady KJ gap | steady pressure gap | steady KJ gap |
|---|---|---:|---:|---:|---:|
| Dir capped vs Neu | `:tri` | 25.36% | 12.85% | 13.40% | 3.00% |
| Dir capped vs Neu | `:tri_robust` | 18.15% | 12.85% | **4.86%** | 3.00% |
| Dir uncapped vs Neu | `:tri` | 26.38% | 12.84% | 14.58% | 3.00% |
| Dir uncapped vs Neu | `:tri_robust` | 19.21% | 12.84% | 6.08% | 3.00% |
| Neumann, `:tri_robust` vs `:tri` | — | 10.54% | 0 | 10.46% | 0 |
| Dirichlet capped, `:tri_robust` vs `:tri` | — | 1.91% | 0 | 1.63% | 0 |
| Dirichlet uncapped, `:tri_robust` vs `:tri` | — | 1.82% | 0 | 1.55% | 0 |

Three findings, all scripted:

1. **In the steady semi-infinite solve the two formulations agree on
   circulation to 3.0%** — which is what the withdrawn 2026-07-29 "~3%" number
   was actually measuring, and it was right about the steady case. With the
   robust reconstruction the steady pressure-CL gap is 4.86%, i.e. **the steady
   comparison would pass the 5% criterion**. The residual steady gap is
   dominated by `:tri`, which alone moves the Neumann steady CL by 10.5% while
   moving the Dirichlet one by 1.6%.
2. **The reconstruction sensitivity is almost entirely on the Neumann body**
   (10.5% vs 1.6–1.9%). This is consistent with the repo's prior finding that
   `basis=:tri` diverges on all-triangle meshes and with the sweptwing-mirror
   investigation, which traced a CL gap to triangulation-sensitive velocity
   evaluation rather than to the strength solution.
3. **The circulation gap is created by the unsteady finite wake, not by the
   formulation as such**: 3.0% with a semi-infinite wake, 12.8% with the
   marching `PanelWake`. Consistently, each body's unsteady pressure CL sits at
   1.056 (Neumann `:tri`) / 0.944 (Neumann `:tri_robust`) / 0.904 / 0.887
   (Dirichlet) times its own KJ CL. The 25.4% settled gap is the product of the
   two effects: \(0.870 \times (0.904/1.056) = 0.745\).

That relocates the open question squarely inside item 015's own subject —
near-wake attachment for a finite marching wake — rather than in force
reconstruction, which is where the withdrawn 2026-07-29 note put it.

### Consequence: V1 stays blocked; surviving suspects

Per the pre-registered hard gate, `ssw0` still stops and the SSW zero-lift
matrix, finite-\(\alpha\) matrix, and swept-wing comparison were **not** run or
interpreted. This is a negative gate result, not evidence against Route B or
pressure continuity. Surviving suspects, in the order the evidence ranks them:

1. **First-wake-row offset / near-wake attachment.** The circulation gap is
   absent in the semi-infinite steady solve and appears only under the marching
   `PanelWake` with `Das = \eta\,\Delta t\,U` (\(\eta = 0.25\)). Both bodies get
   the same rigid offset, but a constant-doublet sheet and a vortex-ring sheet
   need not respond to it identically. Next test: sweep \(\eta\) for both
   formulations and check whether the 12.8% KJ gap is an \(\eta\) artifact —
   which would also connect directly to item 014's \(\eta\) ladder.
2. **`grad_mu` reconstruction on an all-triangle mesh.** Quantified, not
   eliminated: 10.5% on the Neumann body. The gate should probably adopt
   `:tri_robust` as its reference basis, but that is a change to the *oracle*
   and should be decided deliberately rather than because it narrows the gap.
3. **Dirichlet influence-matrix regularization.** `Backslash` assembles \(G\)
   at `body.kerneloffset_panel` (`src/FLOWPanel_solver.jl:229`), and `SSWConfig`
   sets `kerneloffset_over_c = 1e-6`, whereas the DJI campaign's convention is
   that the Dirichlet \(G\) is unregularized by design. Untested here.
4. **Mesh resolution.** The battery ran at a single resolution
   (`n_span=12`, `n_airfoil=21`, 480 cells). Nothing in this record shows
   whether the gap refines away; a two-rung ladder is the cheapest way to find
   out and should precede any structural conclusion.

> **All four were tested the same day** — see the gate diagnosis battery
> section below. Outcome: (1) confirmed dominant, (2) largely a resolution
> artifact, (3) eliminated, (4) confirmed as a real second contributor.

### Driver metric fixes (diamond re-run)

Rerun command: `KUTTAV1_STAGE=diamond KUTTAV1_FORCE=true julia --project -t 2
examples/kutta_v1_attribution.jl`. Source: `diamond_summary.csv`.

| cell | \(L/q\) (see below) | \(\widehat{\Delta p}\) pre (from jump run) | \(\widehat{\Delta p}\) post | max \(|c|/C_s\) | mean outer its. | body solves | startup steps | disposition |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| A/jump | 33.1714 | 4.3362e−2 | 4.3362e−2 | 0 | 0 | 0 | 0 | jump |
| A/pressure | 33.8502 | 4.3362e−2 | 3.0967e−5 | 8.8513e−3 | 2.0 | 16 | 0 | converged |
| B/jump | 30.9702 | 1.8133e−1 | 1.8133e−1 | 0 | 0 | 0 | 0 | jump |
| B/pressure | 38.5168 | 1.8133e−1 | 2.1525e−5 | 2.4455e−1 | 6.0 | 26 | 1 | startup_jump, then converged |

Changes, all of which alter numbers in the 2026-07-29 table:

- **One \(\widehat{\Delta p}\) normalization for every cell and both columns**:
  freestream \(q = \tfrac12\rho\|U_\infty\|^2 = 0.654\), now published as its own
  `dp_hat_scale` column. The closure's frozen `pressure_scale` (115.90 for
  A/pressure, 122.05 for B/pressure, undefined for the jump cells) is recorded
  separately as `kutta_pressure_scale` instead of silently normalizing the
  table. The post-closure residuals therefore read 3.10e−5 / 2.15e−5 rather than
  1.75e−7 / 1.15e−7; the physics is unchanged — those are the same residuals
  divided by 0.654 instead of ~120 — but the four cells are now comparable.
- **`max_dp_hat_pre` renamed `max_dp_hat_pre_from_jump_run`.** It was never
  measured inside the pressure run; it is the *jump* run's post-closure residual
  on the same route, planted as the closure's starting point. The name now says
  so.
- **Iteration means exclude `:startup_jump` steps**, which run no nonlinear
  iteration by design, and the excluded steps are reported as `startup_steps`
  and `startup_body_solves`. B/pressure's mean outer iterations is 6.0, not the
  4.5 that averaged the cold start in.
- **"Fixture CL" is not a CL.** The diamond `ForceMonitor` uses
  `NoNormalization` and the driver divides only by \(q\), never by a reference
  area, so the column is area-dimensional lift \(L/q\): comparable across cells,
  meaningless in absolute terms.
- **Dead `KUTTAV1_RHO` removed.**

### Relocated-sensitivity signal to carry into V2/V3

Route B with pressure closure needs a correction of 24.5% of the correction
scale on the diamond fixture (`c_inf_scaled` = 2.4455e−1, versus 8.9e−3 for
Route A), and its jump baseline violates paired-pressure continuity 4.2× worse
than Route A's (1.8133e−1 vs 4.3362e−2). Both are consistent with the closure
doing much more work to repair a TE-anchored geometry than a rigid-`Das` one.
This is exactly the "relocated rather than removed" outcome Phase 1 asked to be
watched for, and it is tracked here as an open V2/V3 question, not a conclusion:
a fixture with a 3.4% reference thickness is not where it should be decided.

### Regression record, and the V0/V1 `runtests_unit_simulate.jl` contradiction

- `runtests_unit_kutta.jl`: 658/658 pass. `runtests_unit_kutta_routeb.jl`:
  62/62 pass. No `src/` file was changed in this remediation — only
  `examples/kutta_v1_attribution.jl`, `examples/suddenly_started_wing.jl`,
  `examples/sweptwing.jl`, and `examples/helper_functions.jl`.
- **The V0 "green" and V1 "not green" reports of `runtests_unit_simulate.jl`
  are both correct.** Reproduced on the same tree by varying only the thread
  count:

  ```bash
  cd test
  julia --project=.. -t 1 -e 'using Test; import FLOWPanel as pnl; using WriteVTK; using LinearAlgebra; include("test_helpers.jl"); include("runtests_unit_simulate.jl")'
  julia --project=.. -t 6 -e 'using Test; import FLOWPanel as pnl; using WriteVTK; using LinearAlgebra; include("test_helpers.jl"); include("runtests_unit_simulate.jl")'
  ```

  At `-t 1` the "simulate! backend split and validation" testset passes 17/17.
  At `-t 6` it errors with `ArgumentError: step cannot be zero` inside
  `FLOWVPM.Estr_direct_multithreaded` on a zero-particle field
  (`runtests_unit_simulate.jl:138`) and reports 11 pass / 1 error. V0 ran the
  regression single-threaded and V1 ran it multithreaded; the failure is a
  pre-existing multithreaded zero-particle bug in the dependency, is unrelated
  to item 015, and is a genuine verification limitation only for
  multithreaded runs.
- **Unverified in this environment: the swept-wing cap path.**
  `simplewing_mirrored` (and therefore `make_sweptwing_body` and
  `examples/sweptwing.jl`) fails inside GeometricTools' Python spline path with
  `Python: IndexError: list index out of range` during `gt.generate_loft`,
  **identically with `caps=:none` and `caps=:flat`**, so the failure is
  pre-existing and independent of this change. The swept-wing caps wiring and
  the restored Neumann default are therefore code-reviewed and parse-checked but
  not executed. The V1 swept-wing stage is gated off behind the failed sanity
  gate regardless.
- Note that `runtests_unit_simulate.jl` cannot be included standalone: it
  relies on names (`cross`, `test_helpers.jl`) that `test/runtests.jl` brings
  into scope, and a bare `include` fails with `UndefVarError: cross` in an
  unrelated testset. Use the preamble above.

### Reproducibility and regression record

Commands (all local, ≤6 threads):

```bash
KUTTAV1_STAGE=diamond KUTTAV1_FORCE=true julia --project -t 2 \
    examples/kutta_v1_attribution.jl
KUTTAV1_STAGE=gate KUTTAV1_FORCE=true julia --project -t 6 \
    examples/kutta_v1_attribution.jl
```

The `gate` stage is new and can be run on its own; `ssw0` runs it if
`ssw_sanity_gate.csv` is absent and hard-stops if it did not pass. Gate CSVs
are now parsed **by header name** rather than by last-column position, and a
**resume no longer masks a failed gate**: `main_kuttav1` re-verifies the stage's
prerequisite gate records before honoring a checkpoint, so a stage that ran and
failed can never present itself as complete.

Files written by the gate: `gate_summary.csv`, `gate_cl_history.csv`,
`gate_spanwise_gamma.csv`, `gate_spanwise_gamma_steady.csv`,
`gate_kj_validation.csv`, `gate_localization.csv`, `ssw_sanity_gate.csv`. Every
number in this section comes from one of them or from arithmetic on their
published columns.

### Out-of-scope riders (noted, not fixed)

- **Non-watertight Dirichlet bodies are accepted silently.** `solve!` warns when
  a watertight `RigidWakeBody` is solved with Neumann
  (`src/FLOWPanel_solver.jl:261-265`) but says nothing when a *non*-watertight
  body is solved with Dirichlet, which is the symmetric error and is what let
  this campaign run open-tip Dirichlet wings unremarked. A symmetric warning is
  proposed as a follow-up; no `src/` change was made in V1.
- SSW example `n_span=1` breakage against the even-`n_span` guard → item 014.
- `TraceCorrected` deprecation follow-ups.

## 2026-07-30 — gate diagnosis battery (η, mesh, kerneloffset)

Plan: `plans/20260730_kutta_v1_review_remediation.md` follow-up. **Reporting
only.** This battery does not touch `ssw_sanity_gate.csv`, the pass criterion,
or `KUTTAV1_GATE_REFERENCE_BASIS`; V1 remains blocked. Per the user, the
`:tri` vs `:tri_robust` oracle choice and any ∇μ-method investigation are
**deferred** and are not decided here.

Command (single-threaded, as requested — also makes the runs deterministic):

```bash
KUTTAV1_STAGE=gate_diagnosis KUTTAV1_FORCE=true julia --project -t 1 \
    examples/kutta_v1_attribution.jl
```

Sources: `gate_diagnosis_summary.csv` (44 cells, same row schema as
`gate_summary.csv`), `gate_diagnosis_gaps.csv`, `gate_diagnosis_cl_history.csv`,
`gate_diagnosis_spanwise_gamma.csv`. The gap arithmetic is `formulation_gap`,
factored out of `run_gate`, so the gate and the battery cannot drift apart.

### Harness self-checks (run before reading any result)

| check | result |
|---|---|
| η=0.25 / 480-cell / 1e-6 cells reproduce `gate_summary.csv` **bitwise** (`settled_CL`, `CL_kj_unsteady`, `steady_semiinfinite_CL`) | PASS |
| `steady_semiinfinite_CL` and `CL_kj_steady` constant across η within each (formulation, basis) group | PASS |
| `CL_kj_*` identical across `:tri` / `:tri_robust` at fixed knobs | PASS |
| re-running `KUTTAV1_STAGE=gate` leaves `ssw_sanity_gate.csv` and `gate_localization.csv` byte-identical | PASS |

The second check is not decoration: `_ssw_steady_cl` sets `Das` to the unit
freestream direction, so the semi-infinite columns are η-invariant *by
construction*. Had they moved with η, the harness would have been wrong.

### A — η ladder: the circulation gap is a near-wake offset artifact

Neumann-uncapped vs Dirichlet-flat-capped, A/jump, 480 cells, `Das = η·Δt·U`.
`unsteady_kj_gap` is basis-independent by construction (verified above), so one
column serves both bases. The semi-infinite steady KJ gap, 3.00%, is the anchor:
it is what the two formulations agree to when the wake is *not* a marching
finite sheet.

| η | KJ gap (circulation) | settled-CL gap, `:tri` | settled-CL gap, `:tri_robust` | spanwise level ratio Dir/Neu |
|---:|---:|---:|---:|---:|
| 0.0625 | 20.34% | 33.54% | 26.97% | 0.7906 |
| 0.125 | 17.03% | 29.89% | 23.11% | 0.8271 |
| 0.25 *(gate)* | 12.85% | 25.36% | 18.15% | 0.8702 |
| 0.5 | 9.19% | 21.22% | 13.57% | 0.9073 |
| 1.0 | 6.56% | 18.05% | 10.02% | 0.9341 |
| 2.0 | 4.88% | 15.89% | 7.60% | 0.9511 |
| 4.0 | 3.90% | 14.59% | **6.15%** | 0.9609 |
| — | **3.00%** (semi-infinite anchor) | 13.40% | 4.86% | — |

**This is the answer to the question the battery was built for.** The
circulation gap is monotone in η and converges to the semi-infinite anchor:
20.34% → 3.90% against an anchor of 3.00%, with the spanwise level ratio walking
from 0.79 to 0.96 toward unity. The gap the gate measured at η=0.25 is
overwhelmingly a **first-wake-row offset artifact**, not a property of the two
formulations. Route A is not `Das`-free, and this is a direct, quantified
instance of that — the same lever item 014 found worth 16% of rotor-hover \(C_T\).

The settled-CL columns fall with η too but flatten well above the KJ column,
because the second factor (the pressure/KJ offset difference) does not depend on
η. At η=4 with `:tri_robust` the settled gap is 6.15% against a 5% criterion.

### B — mesh ladder: the gap refines, and the basis question largely dissolves

η=0.25 throughout; `dt_star` fixed at 0.125, so this is a spatial ladder only.

| cells | KJ gap | settled gap, `:tri` | settled gap, `:tri_robust` | steady pressure gap, `:tri` | steady pressure gap, `:tri_robust` | steady KJ gap |
|---:|---:|---:|---:|---:|---:|---:|
| 480 (`n_span`=12, `n_airfoil`=21) | 12.85% | 25.36% | 18.15% | 13.40% | 4.86% | 3.00% |
| 1920 (`n_span`=24, `n_airfoil`=41) | 9.77% | 11.74% | 12.09% | **1.37%** | **1.45%** | 1.75% |

Three things happen at once, and the third is the interesting one:

1. **The unsteady gap roughly halves** at `:tri` (25.36% → 11.74%). It is not a
   fixed formulation offset.
2. **The steady semi-infinite comparison converges properly**: the pressure gap
   falls to ~1.4% and the KJ gap to 1.75% at both bases. On the steady problem
   the two formulations agree in the limit, as they must — a genuine
   verification result, and the strongest evidence yet that nothing is
   structurally wrong with either body.
3. **The `:tri` vs `:tri_robust` spread collapses from 7.2 pp to 0.35 pp**
   (11.74% vs 12.09%). The 10.5% basis sensitivity that dominated the gate at
   480 cells is a *resolution* artifact, not a standing property of the
   reconstruction. This does not decide the oracle question — deferred per the
   user — but it does suggest the question may be close to moot at adequate
   resolution, which is worth knowing before spending effort on the ∇μ method.

### C — kerneloffset ladder: the Dirichlet regularization is inert, and 0 fails silently

Dirichlet-capped and Neumann-uncapped, `:tri`, η=0.25, 480 cells.

| `kerneloffset_over_c` | settled gap | KJ gap | status |
|---:|---:|---:|---|
| 1e−4 | 25.36% | 12.85% | ok |
| 1e−6 *(current default)* | 25.36% | 12.85% | ok |
| 1e−8 | 25.35% | 12.84% | ok |
| 0 | — | — | Neumann non-finite |

Over four orders of magnitude the settled CL is identical to 8–9 significant
figures (0.306989163 at 1e−4 and 1e−6; 0.306993025 at 1e−8). **The Dirichlet
`G` regularization is not the cause of the gap** — suspect eliminated.

Two things are worth recording about the `0` rung:

- **The Neumann body returns NaN without throwing.** The vortex-ring near-wake
  self-influence is singular at zero offset and the NaN propagates silently
  through the solve to the reported CL. The battery originally recorded this
  cell as `status="ok"`; `_run_gate_cell` now checks finiteness explicitly and
  reports `status="nonfinite"`. A silent NaN is worse than an exception, and
  this is exactly the kind of thing that would have been quoted as a result.
- **The Dirichlet body runs at zero offset but is a different problem**: settled
  CL 0.237341 (vs 0.306989) and \(C_L^{KJ}\) −0.371968 (vs −0.339535), while its
  *steady* CL is bit-identical at 0.364783. So `kerneloffset` never enters the
  Dirichlet `G` self-term at all; it acts only through the near-wake panel
  influence. That the effect is confined to the unsteady near wake is one more
  independent pointer at the same place as batteries A and B.

### What this leaves

Settled by the battery: the near-wake offset dominates the circulation gap;
mesh resolution is a real second contributor; caps (3.9%, earlier) and Dirichlet
`G` regularization (inert) are both eliminated.

Not answered here, and deliberately not run: **the combined test**, η large
*and* the refined mesh. η=4 alone gives 6.15% (`:tri_robust`) and the refinement
alone gives ~12%, and both act on partly different factors, so the combination
is the one number that would say whether the gate is attainable at converged
settings. That is a decision to take deliberately, not a knob to keep turning
until the gate passes.

Deferred per the user: the `:tri`/`:tri_robust` oracle choice and the ∇μ-method
investigation — though B suggests refinement may matter more than the choice.

**V1 remains blocked and the gate record is untouched.** Reopening it is a
separate decision, and nothing in this battery was allowed to alter it.

## 2026-07-30 — gate diagnosis II: small-Das attribution

Plan: `plans/20260730_kutta_small_das_attribution.md`. **Reporting only.** This
stage does not touch `ssw_sanity_gate.csv`, the 5% criterion, or
`KUTTAV1_GATE_REFERENCE_BASIS`; V1 remains blocked, and reopening it is a
separate decision. T4 was **not** run — it is gated on explicit user approval.

Command (single-threaded, per the standing directive; also deterministic):

```bash
KUTTAV1_STAGE=gate_diagnosis2 KUTTAV1_FORCE=true julia --project -t 1 \
    examples/kutta_v1_attribution.jl
```

Sources: `gate_diagnosis2_summary.csv` (21 cells), `_recovery.csv`, `_gaps.csv`,
`_strip_jump.csv`, `_checks.csv`, `_cl_history.csv`, `_spanwise_gamma.csv`, all
in `data/kutta_v1_attribution/`. Runtime ~40 min.

### Why this battery exists

The previous battery reported "η is the dominant lever" because the gap falls
monotonically toward the semi-infinite anchor as η grows. **That framing is
nearly tautological** — `Das` *is* the finite stand-in for the semi-infinite
row, so agreement at large `Das` is what the construction guarantees. The
non-trivial finding is the opposite end: measured against **its own**
semi-infinite anchor, *both* bodies lose bound circulation as Das→0, and the
Dirichlet body loses ~3.5× more. That **differential** sensitivity is what
needed a cause. The metric here is therefore circulation recovery
$|C_L^{KJ}(\eta)| / |C_L^{KJ,\text{steady}}|$, per body against its own anchor,
not the Dirichlet-vs-Neumann gap.

### Harness self-checks (`gate_diagnosis2_checks.csv`, run before any reading)

| check | value | result |
|---|---|---|
| η=0.25 `neumann` and `dirichlet_vts` reproduce `gate_diagnosis_summary.csv` **bitwise** (`settled_CL`, `CL_kj_unsteady`, `steady_semiinfinite_CL`) | 0 mismatched columns | PASS |
| `CL_kj_steady` invariant per body across all 21 cells (13 Dirichlet, 8 Neumann) | spread exactly 0 | PASS |
| Neumann/Dirichlet meshes identical (nodes, cells) on the uncapped mesh | 0 | PASS |
| Neumann/Dirichlet **control points** identical | max abs diff 0 (scale 2.833) | PASS |
| `Das` parallel to U∞, max over 21 cells | 0 | PASS |
| `Das` unchanged from step 0 to final step, max over 21 cells | 0 | PASS |

The anchors are −0.399434 (Neumann, uncapped) and −0.387432 (Dirichlet,
flat-capped), unchanged from battery I.

### The two questions that needed code reading, not experiment

**Is `Das` a suspect through its direction or through drift?** No.
`_set_ssw_Das!(wing, config.eta*dt*full_uinf)` makes it freestream-parallel by
construction, and `set_Das_eta_freestream=NaN` makes `initialize_Das!` return
immediately and disables the per-step refresh, leaving only `rotate_Das!` — and
this body never rotates. Both are now *measured* rather than asserted: parallel
error and drift are exactly 0 on every cell.

**Could the discrepancy live in wake→control-point evaluation?** No.
`calc_controlpoints!` places the control point at the triangle centroid with
**zero offset** for both formulations (the "exterior for Neumann, interior for
Dirichlet" line in the `Backslash` docstring is stale; self-influence uses
side-aware closed-form limits selected by `DBC`). With identical geometry and
identical wake strengths the wake-induced velocity at the control points is
identical *by construction* — verified above, max abs CP difference exactly 0.
The discrepancy must therefore live in what each formulation **does** with that
velocity. That is what T1 tests.

### T1 — formulation axis: CONFIRMED

Three arms, `:tri`, 480 cells, A/jump. `dirichlet_green` is the same Dirichlet
flat-capped body under `GreenReconstruction` (defaults: `gauge=:area_mean`,
`recompute_interval=1`, `green_solver=nothing`), which transmits the shed wake
as a reconstructed potential trace $q$ — solve $(I-B)q = S\sigma$, then
$G\mu_E = -S\sigma_0 - q$ — instead of the default `VelocityThroughSources`,
which represents it by an equivalent source sheet $\sigma = -u_{wake}\cdot n$ and
drops its direct contribution to the interior-potential condition (the exact RHS
term is $-\phi_{wake}$).

Circulation recovery $|C_L^{KJ}(\eta)| / |C_L^{KJ,\text{steady}}|$:

| η | Neumann | Dirichlet VTS | Dirichlet **Green** |
|---:|---:|---:|---:|
| 0.0625 | 0.9342 | 0.7672 | **0.9232** |
| 0.125 | 0.9619 | 0.8228 | **0.9363** |
| 0.25 *(gate)* | 0.9753 | 0.8764 | **0.9494** |
| 0.5 | 0.9824 | 0.9198 | **0.9621** |
| 1.0 | 0.9865 | 0.9504 | **0.9729** |

Read as circulation **loss** (1 − recovery), and as the Dirichlet-minus-Neumann
**differential**, which is the quantity that needed explaining:

| η | Neumann loss | VTS loss | Green loss | differential, VTS | differential, Green |
|---:|---:|---:|---:|---:|---:|
| 0.0625 | 6.58% | 23.28% | 7.68% | 16.70 pp | **1.10 pp** |
| 0.125 | 3.81% | 17.72% | 6.37% | 13.91 pp | **2.56 pp** |
| 0.25 | 2.47% | 12.36% | 5.06% | 9.89 pp | **2.59 pp** |
| 0.5 | 1.76% | 8.02% | 3.79% | 6.26 pp | **2.03 pp** |
| 1.0 | 1.35% | 4.96% | 2.71% | 3.61 pp | **1.36 pp** |

**The hypothesis is confirmed.** Under `VelocityThroughSources` the differential
grows 3.61 → 16.70 pp as η falls by 16×, i.e. it varies by 13.1 pp across the
ladder. Under `GreenReconstruction` it is **flat at 1.1–2.6 pp**, varying by
1.5 pp, with no systematic trend in η. At the smallest offset — where the
substitution should be worst, and is — `GreenReconstruction` removes **93%** of
the differential (16.70 → 1.10 pp).

The formulation gap tells the same story from the other side
(`gate_diagnosis2_gaps.csv`):

| η | KJ gap, VTS | KJ gap, **Green** | settled-CL gap, VTS | settled-CL gap, **Green** |
|---:|---:|---:|---:|---:|
| 0.0625 | 20.34% | **4.15%** | 33.54% | 14.70% |
| 0.125 | 17.03% | **5.59%** | 29.89% | 15.86% |
| 0.25 | 12.85% | **5.58%** | 25.36% | 15.79% |
| 0.5 | 9.19% | **5.02%** | 21.22% | 15.25% |
| 1.0 | 6.56% | **4.34%** | 18.05% | 14.62% |
| — | 3.00% (semi-infinite anchor) | — | 13.40% | — |

Under `GreenReconstruction` the circulation gap loses its η-dependence entirely
and sits at 4–6% against the 3.00% semi-infinite anchor at every offset. **So the
Dirichlet arm's η-sensitivity is substantially an artifact of the default
wake→body transfer scheme, not physics.**

Two limits worth stating plainly:

- **Green does not close the settled-CL gap.** It flattens it (14.6–15.9% across
  the whole ladder, vs 18.1–33.5% for VTS) but leaves a ~15% floor. That floor is
  the *force-reconstruction* factor, which battery I already localized: the
  steady semi-infinite pressure gap is 13.40% at 480 cells and falls to **1.37%
  at 1920 cells**. The residue is resolution, not formulation.
  *(2026-07-30: this was an inference across two batteries when written. It is
  now a measurement — see "gate diagnosis II(b)" below, where Green's unsteady
  reconstruction ratio equals the steady one to ≤0.11% at three mesh rungs and
  the floor's contribution falls to 0.28 pp at 4320 cells. The circulation half,
  however, is NOT converged.)*
- **The Neumann body's own small-Das loss is real but small** (6.58% at
  η=0.0625), and T2 shows most of it is wake rollup.

### T2 — geometry control: the rollup confound is exonerated

Battery I ran with `shed_with_induced_velocity=true, freestream_convection=false`,
so each body convected its wake with its **own** induced velocity and the two
wakes rolled up differently from step 1. `freestream_convection=true` translates
every row rigidly by `dt·U∞`, making the sheet straight, uniform, and identical
between the two runs.

| η | Neumann recovery, rolled up → straight | Dirichlet VTS recovery, rolled up → straight | differential, rolled up → straight |
|---:|---|---|---|
| 0.0625 | 0.9342 → **0.9786** | 0.7672 → **0.7701** | 16.70 → **20.85 pp** |
| 0.25 | 0.9753 → **0.9847** | 0.8764 → **0.8784** | 9.89 → **12.16 pp** |
| 1.0 | 0.9865 → **0.9884** | 0.9504 → **0.9514** | 3.61 → **4.86 pp** |

**Removing the wake-shape difference makes the differential larger, not
smaller.** Two-thirds of the Neumann body's small-Das loss was rollup (6.58% →
2.14% at η=0.0625); essentially none of the Dirichlet body's was (23.28% →
22.99%). The residual η-sensitivity survives on a straight, uniform, identical
sheet — the regime in which the cancelling-interior-edge argument
(`examples/suddenly_started_wing.jl:757-766`) says `Das` should be nearly inert.
**Wake shape cannot explain the asymmetry**, and it was partly masking it.

### T3 — strip / wake-row-1 strength jump: eliminated

The attached `Das` strip carries the **current** bound circulation (it is in the
LHS); wake row 1 carries the **previous** step's, planted by `shed_wake!`. A
read-only monitor recorded $\Gamma_{strip} - \Gamma_{row1}$ per shedding edge per
step, before `propagate!`/`shed_wake!`, on all 21 runs. Settled window
t*∈[7,8], normalized by the settled-window maximum |Γ| of that cell (one scale
per cell, so near-tip edges where Γ→0 cannot manufacture a ratio):

| η | Neumann | Dirichlet VTS | Dirichlet Green |
|---:|---:|---:|---:|
| 0.0625 | 8.1e−4 | 5.3e−4 | 8.0e−4 |
| 0.25 | 6.0e−4 | 4.7e−4 | 5.5e−4 |
| 1.0 | 4.1e−4 | 3.7e−4 | 3.9e−4 |

(mean |jump|/Γ_max over edges; per-cell maxima are 5.4e−4…1.4e−3.)

The jump is **three orders of magnitude below** the 10–20% circulation gap, and
its ordering across arms is *anti*-correlated with the gap: the Dirichlet VTS
body, which has by far the largest gap, has the **smallest** jump at every η,
while `GreenReconstruction` — which fixes the gap — moves the jump back up
toward Neumann's. It does fall with η, in the same direction as the gap, but it
cannot carry an effect of that size and does not order correctly.
**Eliminated as the cause.**

### Verdict

The differential small-Das sensitivity is attributed to the **default wake→body
transfer scheme for the Dirichlet formulation**: `set_strengths!` represents the
shed wake by an equivalent source sheet derived from its normal velocity
($\sigma = -u_{wake}\cdot n$) and drops its direct contribution to the
interior-potential condition, whose exact RHS term is $-\phi_{wake}$. As Das→0
the shed sheet approaches the TE, its potential contribution grows, and the
normal-velocity proxy degrades — on the Dirichlet side only, because the Neumann
body consumes $u\cdot n$ directly: that *is* its boundary condition, unconverted.
Replacing the proxy with `GreenReconstruction`'s reconstructed potential trace
removes 93% of the differential at the smallest offset and all of its
η-dependence.

Eliminated by this battery: wake rollup / wake shape (T2, makes it worse), the
strip/row-1 strength jump (T3, three orders too small and wrongly ordered),
`Das` direction and drift (measured 0), and wake→control-point evaluation
(control points identical to the bit). Previously eliminated: flat tip caps
(3.9% of the gap) and Dirichlet `G` regularization (inert over 1e-8…1e-4).

**Not run:** T4 (T4a local operator error, T4b global closure identity), gated on
explicit user approval. T4a costs ~1 min, needs no solve, and would measure the
operator error itself rather than its effect on CL — the most direct test of the
mechanism above. T5 (strip aspect ratio) stays dropped; the plan records why.

### Amendment to the 2026-07-30 battery-I section

The battery-I headline — "the circulation gap is a near-wake offset artifact",
"the gap the gate measured at η=0.25 is overwhelmingly a first-wake-row offset
artifact" — **stands as an arithmetic description and is superseded as an
attribution.** The η-dependence is real, but for the Dirichlet arm it is
substantially a property of `VelocityThroughSources`, not of the near-wake offset
per se: at fixed mesh and fixed `Das`, changing only the formulation flattens the
circulation gap to 4–6% at every η. Battery I's numbers are unchanged and were
reproduced bitwise; only their interpretation is corrected.

**V1 remains blocked and the gate record is untouched** — `ssw_sanity_gate.csv`
and `gate_localization.csv` re-verified byte-identical after re-running
`KUTTAV1_STAGE=gate`. `runtests_unit_kutta.jl` (658) and
`runtests_unit_kutta_routeb.jl` (62) pass single-threaded.

### Driver changes (`examples/kutta_v1_attribution.jl`, no `src/` edits)

- `_ssw_combo_run` gains `formulation` and `extra_monitors`; the formulation is
  forwarded to `pnl.simulate!`, extra monitors are **appended** so the existing
  `audit_monitors` ordering contract is untouched.
- `_run_gate_cell` reads optional `formulation` and `freestream_convection` off
  the cell, records both in the row schema, and **errors** (never coerces) if a
  non-default formulation is paired with a non-Dirichlet body.
- `StripJumpRecorder` + `_settled_strip_jump`: the T3 monitor and its reduction.
- `_das_control` / `_controlpoint_control`: the two control assertions.
  `_controlpoint_control` calls `calc_normals!`/`calc_controlpoints!` first —
  control points are populated by the solve, not the constructor, so comparing
  freshly built bodies would otherwise compare two zero arrays and pass
  vacuously.
- New stage `KUTTAV1_STAGE=gate_diagnosis2` → `run_gate_diagnosis2`. **No entry
  in `KUTTAV1_STAGE_GATES`**: it gates nothing and is gated by nothing.

### 2026-07-30 — gate diagnosis II(b): mesh ladder under GreenReconstruction

The section above attributed the ~15% settled-CL floor that survives
`GreenReconstruction` to "the force-reconstruction factor, which battery I
already localized". **That was an inference stitched from two separate
batteries, not a measurement.** This sub-battery measures it.

```bash
KUTTAV1_STAGE=gate_diagnosis2 KUTTAV1_DIAG2_BATTERIES=mesh KUTTAV1_FORCE=true \
    SSW_NO_PLOT=true julia --project -t 1 examples/kutta_v1_attribution.jl
```

9 runs (3 rungs × 3 arms) at η=0.25, `:tri`, ~55 min single-threaded. Sources:
`gate_diagnosis2_*_mesh.csv`. **Three** rungs, not two — a two-rung agreement is
not convergence (the DJI Phase 2c lesson).

#### The decomposition being tested

The settled-CL ratio is the *exact* product of two factors (`gap_factorization`,
`examples/kutta_v1_attribution.jl`; the identity residual is a published check,
max 1.11e−16 over 6 pairs):

$$\frac{C_{L,\text{dir}}}{C_{L,\text{neu}}}
= \underbrace{\frac{|C^{KJ}_{\text{dir}}|}{|C^{KJ}_{\text{neu}}|}}_{\text{circulation}}
\times
\underbrace{\frac{C_{L,\text{dir}}/|C^{KJ}_{\text{dir}}|}{C_{L,\text{neu}}/|C^{KJ}_{\text{neu}}|}}_{\text{reconstruction}}$$

The comparator is the same reconstruction ratio built from the **steady
semi-infinite** solve, which no formulation kwarg can reach (`_ssw_steady_cl`
runs `steady!` on a `semiinfinite_wake=true` body, where `GreenReconstruction`
is rejected by `_validate_formulation_common`). An unsteady reconstruction ratio
equal to it is the steady discretization error and nothing more.

**Pre-registered prediction** (written before the runs): under
`GreenReconstruction` the unsteady reconstruction ratio tracks the steady one at
every mesh, and the settled-CL gap collapses onto the circulation gap alone.

#### Result: CONFIRMED

| cells | arm | settled gap | KJ gap | circulation | **reconstruction** | **steady reconstruction** |
|---:|---|---:|---:|---:|---:|---:|
| 480 | VTS | 25.36% | 12.85% | 0.8715 | 0.8565 | 0.8928 |
| 480 | **Green** | 15.79% | 5.58% | 0.9442 | **0.8919** | **0.8928** |
| 1920 | VTS | 11.74% | 9.77% | 0.9023 | 0.9782 | 1.0039 |
| 1920 | **Green** | 4.89% | 5.22% | 0.9478 | **1.0035** | **1.0039** |
| 4320 | VTS | 9.09% | 7.71% | 0.9229 | 0.9851 | 1.0031 |
| 4320 | **Green** | 4.15% | 4.42% | 0.9558 | **1.0028** | **1.0031** |

Green's unsteady reconstruction ratio tracks the steady comparator to
**0.110% / 0.040% / 0.034%** across the three rungs. `VelocityThroughSources`
does not — it departs by **4.08% / 2.56% / 1.80%**, i.e. the default transfer
scheme injects a *spurious* reconstruction error on top of its circulation
error, and that contamination shrinks with refinement but has not vanished at
4440 cells.

Both halves of the prediction hold:

1. **The floor is force reconstruction, measured.** Under Green the
   reconstruction factor goes 0.8919 → 1.0035 → 1.0028, i.e. its contribution to
   the gap falls from 10.8 pp at 480 cells to **0.28 pp** at 4320 — and it is
   the steady discretization error at every rung, to ≤0.11%.
2. **The settled gap collapses onto the circulation gap.** At the two fine rungs
   settled ≈ KJ (4.89% vs 5.22%; 4.15% vs 4.42%), where at 480 cells they were
   15.79% vs 5.58%.

Green also removes the *mesh* sensitivity of circulation recovery, not just the
η sensitivity: `kj_recovery` is 0.9494 / 0.9448 / 0.9479 (flat to ±0.5%) while
VTS climbs 0.8764 → 0.8994 → 0.9152 and is still climbing.

#### What has NOT converged, and should not be glossed

The reconstruction half is converged — the steady comparator moves only 0.08%
between the two fine rungs (1.0039 → 1.0031). **The circulation half is not.**
Under Green the KJ gap is 5.58% → 5.22% → 4.42%, still falling, and its own
steady anchor (the steady KJ gap, from `steady_circulation_ratio`) is
3.00% → 1.75% → 1.19%, also still falling. At the finest mesh the unsteady
circulation gap is still ~3.7× its steady counterpart. So the residue is now
**entirely circulation, and open** — this battery moved the question, it did not
close it.

#### Consequence for the gate — reported, not acted on

At 1920 and 4320 cells the Dirichlet-Green settled gap is **4.89% and 4.15%**,
below the gate's 5% criterion. Stated because it is what the numbers say.
**The gate record is untouched**: `ssw_sanity_gate.csv` and
`gate_localization.csv` re-verified byte-identical (md5) after re-running
`KUTTAV1_STAGE=gate`, the criterion and `KUTTAV1_GATE_REFERENCE_BASIS` are
unchanged, and V1 remains blocked. Per the user's decision this session,
reopening the gate at converged settings is a **separate decision** and is
deliberately not taken here — not least because the circulation half is not
converged, so a pass at these two rungs would be a pass on a still-moving number.

#### Amendment to the section above

The sentence "That floor is the *force-reconstruction* factor, which battery I
already localized" was an inference; it is now a **measurement**, made at three
rungs, with the steady semi-infinite reconstruction ratio as the independent
comparator. The supporting figure quoted there (steady pressure gap
13.4% → 1.4% over 480→1920) is battery I's and stands.

#### Harness self-checks (all PASS, `gate_diagnosis2_checks_mesh.csv`)

| check | value |
|---|---|
| rung 1 (480) reproduces the T1 η=0.25 cells **bitwise**, all three arms | 0 mismatched columns each |
| rungs 1–2 `neumann` / `dirichlet_vts` reproduce battery I's meshladder rows **bitwise** | 0 mismatched columns each |
| `gap_factorization` identity, max over 6 pairs | 1.11e−16 |
| `CL_kj_steady` invariant **per (body, mesh)** — 6 groups | spread exactly 0 |
| `Das` parallel error / drift, max over 9 cells | 0 / 0 |
| Neumann vs Dirichlet control points identical | 0 (scale 2.833) |

The steady-anchor check was regrouped for this battery: it previously keyed on
`(bodytype, caps)` alone, which would have failed spuriously the moment the mesh
varied, because `CL_kj_steady` legitimately refines (Neumann −0.399434 →
−0.396583 → −0.395941; Dirichlet −0.387432 → −0.389629 → −0.391219).

#### Driver changes

- `gap_factorization(n, d)` next to `formulation_gap`, so the gate and both
  batteries keep sharing one piece of arithmetic; its columns
  (`circulation_ratio`, `reconstruction_ratio`, `steady_reconstruction_ratio`,
  `steady_circulation_ratio`, `settled_ratio`, `factorization_residual`) are
  published in `gate_diagnosis2_gaps*.csv` rather than left as prose arithmetic.
  `pressure_per_kj` / `steady_pressure_per_kj` added per cell to
  `gate_diagnosis2_recovery*.csv` so the split can be re-derived from the
  recovery table alone.
- `run_gate_diagnosis2` is battery-selectable
  (`KUTTAV1_DIAG2_BATTERIES=formulation,convection,mesh`) with **suffixed output
  files** when a subset is selected, mirroring `run_gate_diagnosis`, so the mesh
  ladder cannot overwrite the completed T1/T2 record.
- `_repro_check` / `_prior_rows` factor the bitwise-reproduction comparison,
  which now covers three families of cells.
