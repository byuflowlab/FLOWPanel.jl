# Phase 2b Log — Formulation/Topology Attribution

Plan: [Phase 2b — Formulation/Topology Attribution](../../plans/dji_convergence_20260722/phase_02b_formulation_attribution.md)

Dashboard: [DJI convergence progress](../../plans/dji_convergence_20260722/dji_convergence_20260722.md)

## Current snapshot

Status: **COMPLETE (pending Ryan's Phase 3 approval).** Attribution reviewed and
confirmed, and **convergence obtained**: the Dirichlet-Neumann circulation gap is a
**chordwise (n_airfoil) section-resolution discretization difference** between the
Morino source+doublet (Dirichlet) and vortex-ring (Neumann) formulations. Under
chordwise refinement the clean solves converge to **0.11% at 120 panels/section**
(open Neumann converges by ~60; Dirichlet climbs up to meet it). NOT an extraction
artifact, NOT a solve-formulation/Kutta-trace effect, NOT explained by thinness. The
earlier "~1% physical caps term" was CORRECTED — it was the drop-one-panel artifact;
the genuine cap effect is ≲0.1–0.2%.

### Extended chordwise ladder (2026-07-23) — CONVERGENCE OBTAINED; caps claim corrected

Added rungs xfine/xxfine/ultra = n_airfoil 81/101/121 (= 80/100/120 chordwise
panels/section) at fixed n_span=11, n_endcap=7 (chordwise is the gap lever).

Oracle wing geometry: rectangular NACA 0015, chord 1, span 4, **aspect ratio = 4**.
Panel counts per rung ("chordwise/section" = panels around the section perimeter =
2× per-surface chordwise; "spanwise" = TE shedding stations = n_span):

| rung | spanwise | chordwise/surface | chordwise/section | lifting-surf panels | total (capped) |
|---|---:|---:|---:|---:|---:|
| coarse | 5  | 15 | 30  | 300  | 524  |
| medium | 7  | 20 | 40  | 560  | 864  |
| fine   | 11 | 30 | 60  | 1320 | 2016 |
| xfine  | 11 | 40 | 80  | 1760 | 2696 |
| xxfine | 11 | 50 | 100 | 2200 | 3376 |
| ultra  | 11 | 60 | 120 | 2640 | 4056 |

(open Neumann = lifting-surf panels only, no caps; capped total adds the tip caps.
"120 chordwise panels/section" at convergence = 60 panels along the chord on each of
the upper and lower surfaces.)

| nsec | cap+Dir ∫Γ | open+Neu ∫Γ | production gap | Dir Δ∫Γ/rung | Neu Δ∫Γ/rung |
|---:|---:|---:|---:|---:|---:|
| 60  | −66.380 | −66.895 | 0.78% | (−4.2% span) | (−3.5% span) |
| 80  | −66.594 | −66.891 | 0.45% | −0.32% | +0.007% |
| 100 | −66.721 | −66.884 | 0.24% | −0.19% | +0.010% |
| 120 | −66.805 | −66.878 | **0.11%** | −0.13% | +0.010% |

(production gap = open+Neumann vs capped+Dirichlet — the clean, trustworthy solves,
and the direct analogue of Phase-1's uncapped/Neumann vs capped/Dirichlet.)

**Convergence obtained.** Open Neumann is chordwise-converged by ~60 panels/section
(flat at −66.88, Δ≈0.01%); capped Dirichlet (Morino) converges more slowly and climbs
*up* to meet it (Δ 0.32→0.19→0.13%). The two clean solves agree to **0.11% at 120
chordwise panels/section** — below the 1% and the stricter 0.25% criteria. The Dir-Neu
gap is fully explained by chordwise under-resolution of the Dirichlet formulation, which
converges to the Neumann value.

**CORRECTION to the earlier "~1% physical caps/topology term":** that was mostly an
artifact. It came only from the fragile drop-one-panel capped-Neumann (near-singular).
Since capped Dirichlet ≈ open Neumann to 0.11% in the converged limit, the genuine
physical cap effect on ∫Γ is small (≲0.1–0.2%); the persistent ~0.8% `capped-Neumann
(drop1)` offset (BConly 0.92% at ultra) is the drop-regularization artifact, not clean
physics, and is disregarded in favor of the production comparison.

**Caveat (absolute ∫Γ):** the large early per-rung ∫Γ changes (−5.5%, −4.2%) were mostly
*spanwise* (n_span 5→7→11), not chordwise; at fixed n_span=11 the pure-chordwise ∫Γ
change plateaus (0.32→0.19→0.13%). The *gap* is converged; certifying the absolute ∫Γ
would need a separate spanwise sweep (not required for the attribution).

Artifacts: `oracle_stations_{xfine,xxfine,ultra}.csv`, updated `oracle_summary.csv`.

### Confirmations (2026-07-23) — thickness screen + refinement-dimension separation

**Thickness screen** (`thickness` mode, medium refinement (7,41), semi-infinite wake):

| t/c | cap+Dir | open+Neu | cap+Neu(d1) | BConly% | topo% | extract_res |
|---:|---:|---:|---:|---:|---:|---:|
| 0.06 | −61.043 | −61.564 | −62.247 | −1.972 | +1.098 | 3.1e-5 |
| 0.12 | −62.818 | −63.603 | −64.260 | −2.295 | +1.022 | 3.1e-5 |
| 0.18 | −64.585 | −65.632 | −66.335 | −2.711 | +1.061 | 3.0e-5 |

|BConly%| **grows with thickness** (thinner → smaller gap). So global thinness does
NOT explain DJI's ~5% > oracle's ~2.5%; thinner sections would predict a *smaller*
gap. `topo%` ~1% and extraction residual 3e-5 are thickness-independent.

**Refinement-dimension separation** (medium baseline (7,41), n_endcap=5 fixed):

| case | BConly% | topo% |
|---|---:|---:|
| baseline (n_span=7, n_airfoil=41) | −2.494 | +1.036 |
| chordwise refine (7, 61) | **−1.602** | +0.780 |
| spanwise refine (11, 41) | −2.690 | +1.219 |
| both (11, 61) | −1.700 | +0.879 |

**Chordwise resolution (`n_airfoil`) is the entire convergence lever**: it closes
BConly ~0.9%; spanwise refinement alone does not help (slightly widens); spanwise
on top of chordwise adds nothing. The Dir-Neu difference lives in the chordwise
resolution of the airfoil *section* surface (LE suction peak + TE closure), where
the source+doublet and vortex-ring representations differ until resolved. This
reconciles the thickness result: at fixed `n_airfoil`, a thicker section is
relatively under-resolved chordwise → bigger gap.

Robustness: one refinement step per axis, but the direction is large and
unambiguous (chordwise −0.9%, spanwise +0.2%). Extraction residual stayed 3e-5.
`formulation_test.jl` passes all 10 stages after the `pitching_wing.jl` caps/thickness
refactor (default behaviour unchanged).

### Oracle Gate A ladder + tiebreaker result (2026-07-23) — the gap is numerically unresolved

Refactored `examples/pitching_wing.jl`: `pitching_wing_mesh` now takes `caps::Bool`
(open lifting surface vs capped, on identical lifting-surface coordinates — verified
open == capped main node/cell block exactly) and `thickness::Real`; threaded
`thickness` into `build_pitching_wing_body`. Existing default behaviour (caps=true,
t/c=0.15) unchanged. Driver `oracle` + `analyze` modes implemented. All three oracle
bodies solved through one matched **semi-infinite attached wake** via
`steady!`/VelocityThroughSources, DirectBackend:

- capped+Dirichlet: `RigidWakeBody{Union{ConstantSource,ConstantDoublet},2,Float64,true}`
- open+Neumann: `RigidWakeBody{VortexRing,1,Float64,false}`, caps=false, watertight=false
- capped+Neumann(drop-one-panel): same VortexRing type, mid panel (ncells÷2) dropped

| refine | ncells(cap) | cap+Dir ∫Γ | open+Neu ∫Γ | cap+Neu(d1) ∫Γ | BConly% | topo% |
|---|---:|---:|---:|---:|---:|---:|
| coarse | 524  | −60.395 | −61.613 | −62.314 | −3.177 | +1.125 |
| medium | 864  | −63.707 | −64.619 | −65.295 | −2.494 | +1.036 |
| fine   | 2016 | −66.380 | −66.895 | −67.562 | −1.780 | +0.987 |

`BConly%` = capped Neumann vs capped Dirichlet (topology fixed → pure boundary
condition). `topo%` = open vs capped Neumann (BC fixed → pure caps/closure).

**Findings:**

1. **Not an extraction/mix artifact.** Gate A velocity-based loop-integral
   circulation equals `Gamma_TE = μ_upper−μ_lower` to `max|loopcirc−Γ|/|Γ| ≈ 3e-5`
   for **all** bodies (Dirichlet and both Neumann), at every rung. The Dirichlet
   circulation faithfully represents its exterior field; the Dir-Neu difference is a
   genuine field difference, not a weakly-observable Dirichlet mix. Confirms the
   Step 0 prior conditioning evidence. → representation/extraction branch CLOSED.

2. **Pure-BC gap shrinks monotonically under refinement:** BConly −3.177% → −2.494%
   → −1.780% (524→864→2016 cells), no plateau. The production-like open-Neu-vs-Dir
   gap likewise shrinks −2.02% → −1.43% → −0.78%. This is the Phase-1 pattern
   (Neumann drifting toward Dirichlet under refinement). → the Dir-Neu gap is
   **dominantly numerically unresolved** (source+doublet Morino vs vortex-ring
   Neumann discretizations converging toward each other).

3. **Caps/topology contribute a small, nearly-converged ~1%** (topo% +1.125 →
   +1.036 → +0.987), separable and persistent — physical (caps close the tips).

4. **Drop-one-panel tiebreaker is drop-location-sensitive** (fragility finding):
   mid-panel drop clean & consistent, but an alternate drop (coarse, panel 63 near a
   TE/strip-2 location) swung ∫Γ by −47.6% — the closed VortexRing system is
   near-singular and only conditionally regularized. **Robustness safeguard:** the
   mid-drop BConly% matches the drop-free decomposition (open-vs-Dir%) − (topo%) to
   ~0.02% at all three rungs, so the attribution does NOT hinge on the fragile drop
   case; **open+Neumann is the trustworthy (genuinely full-rank) Neumann reference.**

**Reviewed attribution (exit-gate answer): MIXED, quantified —**
- **dominant: numerically-unresolved discretization difference** between the
  Dirichlet source+doublet (Morino) and Neumann vortex-ring boundary conditions,
  converging toward zero under refinement (residual at 2016 cells: BC-only −1.78%,
  production −0.78%; no plateau observed);
- **small separable physical caps/topology term ~+1%** (nearly converged);
- **excluded:** solve-formulation/Kutta-trace channel (Gate A0, 0.74%) and Dirichlet
  representation/extraction (loop-circ = Γ_TE to 3e-5).

Refinement dimension: `n_span` and `n_airfoil` were refined together; which
dominates is not yet separated. Oracle t/c=0.15 (NACA0015); DJI sections are
thinner and the DJI gap (~5%) exceeds the oracle's (~2.5%), so global thickness is
the leading candidate to explain the magnitude difference — motivates the gated
thickness screen and an `n_airfoil`-vs-`n_span` separation as the recommended
confirmations before the Phase 3 handoff.

Artifacts: `oracle_summary.csv`, `oracle_stations_{coarse,medium,fine}.csv`.

### Gate A0 formsweep result (2026-07-23) — formulation channel is NOT the ~5% lever

Driver `examples/dji9443_formulation_attribution.jl` (`smoke` + `formsweep` modes)
built, reusing `test/formulation_test.jl` harness (`small_body`/`flat_wake`/
`run_aero!`/`kutta_map`). `smoke` passed (all machinery). Gate A0 on the unchanged
capped/Dirichlet pitching-wing oracle (864 cells, 7 shedding edges), identical
frozen flat wake across formulations, DirectBackend:

| formulation | ∫Γ dy | Δ vs VTS |
|---|---:|---:|
| VelocityThroughSources | −32.0087 | — |
| TraceCorrected(:green) | −32.2464 | −0.743% |
| TraceCorrected(:line_integral, simpson) | −32.2529 | −0.763% |
| GreenReconstruction | −32.2464 | −0.743% |

- **Disambiguation (important):** "formulation" here = the `AbstractSolveFormulation`
  (wake-potential/Kutta-trace *solve scheme*), which is the same for all three tested
  cases — they are **all Dirichlet solves**. Gate A0 did NOT run Neumann and says
  nothing about the Dirichlet-vs-Neumann (boundary-condition) question, which remains
  the open primary question for the oracle tiebreaker. What Gate A0 excludes is only
  the solve-scheme sub-hypothesis (that the Dirichlet trace correction `c` would lift
  circulation ~5% toward Neumann).
- The affine Kutta correction `c` shifts ∫Γ by only **~0.74%** (toward Neumann,
  i.e. larger magnitude), ~1/7 of the Phase-1 ~5% Dirichlet-Neumann gap. **The
  formulation / wake-potential / Kutta-trace channel is excluded as the dominant
  cause.** Does NOT trigger the plan's ≥5% short-circuit.
- `GreenReconstruction == TraceCorrected(:green)` to displayed precision — the
  expected consistency identity (Stage 6 of formulation_test), a sanity check.
- Corroborates prior `debug/dirichlet_solve/` Task 4/5: trace estimator is not
  the lift/circulation lever (quadrature changed lift ~1e-9 there; here the whole
  `c` channel is ~0.74%).
- Robustness/caveats: (a) oracle NACA0015 finite wing, not the DJI rotor mesh —
  but the affine-`c` math is geometry-general and 0.74% is unlikely to become 5%
  on DJI; plan permits a DJI Gate A0 only "if the oracle result is material,"
  which it is not. (b) fabricated flat wake seeded from a VTS pass (mean edge
  Γ≈−5.48); `c` depends on wake velocity, but prior evidence shows `c` is a small
  TE-trace refinement, not a 5% lever. Conclusion held after method review.
- **Attribution consequence:** the ~5% gap is NOT formulation-dominant. Leading
  remaining candidates are topology/closure (caps) vs a genuine
  Neumann-vs-Dirichlet operator difference — to be separated by the oracle
  open+Neumann / capped+Dirichlet / capped+Neumann(drop-one-panel) tiebreaker
  and Gate A exterior-velocity + loop-integral circulation.

Artifacts: `formsweep_summary.csv`, `formsweep_stations.csv`, `smoke_stations.csv`
under `data/dji_convergence_20260722/phase_02b_formulation_attribution/`.

Open item for oracle build: the loop-integral self-check read ~0 circulation on
the tiny 3-station smoke body (relative spread meaningless near zero) — validate/
debug the loop-integral geometry (wake-sheet crossing + sign) on the real oracle
before trusting Gate A loop circulation.


### Step 0 prior-evidence intake result (2026-07-23, read-only, no solves)

Read `debug/dirichlet_solve/{task5.md, dirichlet_solve.md, green_identity_analysis.md
takeaways}` and `debug/debug_bc_comparison.jl`. Note the prior suite targeted the
*wake-route* discrepancy (semi-infinite vs finite rolled free wake; CL up to −3.44%
at 3.94°, −9.3% at 45°), not directly the Phase-1 bound-circulation Dirichlet-vs-Neumann
gap — but its conditioning result is strong prior evidence here.

Mechanisms **excluded / implicated** by the prior work:

- **Nonuniqueness / weak-mix branch — largely excluded a priori.** `(I−B)` has a
  *single* isolated near-null mode = the constant (`const_projection=1.0`,
  σ_min≈4e-11 at 6688 panels); area gauge removes it, σ_2≈9.5e-3, gauge-fixed
  cond≈100 stable across resolution, **no tiny-σ cluster**. The constant cancels in
  μ_upper−μ_lower, so it cannot explain a Gamma_TE gap. ⇒ Dirichlet extraction is
  well-determined; the mix/nonuniqueness branch is low-prior unless Gate A contradicts.
- **Gauge/solve failure — excluded.** The trace degradation is a *forward
  discretization residual* `r = Sσ − (I−B)q_direct`, centroid-sampled constant-panel
  representation error, converging algebraically (far-field ~h^1.1–1.4; near-field
  wake→body ~h^0.5–1.1, the slow part).
- **Velocity-kernel regularization inconsistency — excluded.** `core_size=0`
  reproduces every residual bit-for-bit (core inert; no body CP within
  ~1e-3 m ≈0.003c of a wake panel).
- **Trace-estimator error — excluded as a lift-gap cause.** Straight Simpson vs
  deformed Simpson vs trapezoid change the trace metric strongly but lift only
  ~1e-9; the remaining Task 3 lift gap is shared fixed-wake geometry error.
- **Defect flag:** the (201,16,11)/10360-panel NACA0015 mesh is defective
  (constant-mode defect 6.6e-3, 8 orders high) — avoid that resolution.

Tooling Phase 2b will **reuse**:

- drop-one-panel closed Neumann build (`debug_bc_comparison.jl:22–66`): drop
  `size(cells,2)÷2` mid-mesh panel → full-rank pure-VortexRing system; reuse for the
  capped+Neumann tiebreaker;
- Green-identity residual `r = Sσ − (I−B)q` (body-only, wake columns excluded) as a
  closed-surface consistency check;
- source/doublet exterior-velocity split;
- `TraceCorrected` (estimator ∈ {:green, :line_integral}) and `GreenReconstruction`
  formulation hooks, both `DirectBackend`-only, for Gate A0.

Implication for phase ordering: Gate A0 (`TraceCorrected` affine `c` shift) remains the
cheapest decisive test; the mix/nonuniqueness follow-ups are now *lower* prior given the
prior conditioning result and should stay gated behind an explicit Gate A contradiction.

---

Status (pre-execution snapshot retained below): **Staged — awaiting Ryan approval**

- Goal: determine why Dirichlet solves produce lower integrated circulation
  than Neumann solves.
- This phase is inserted between Phase 2 and Phase 3.
- Phase 2 accepted sharp capped/Dirichlet TE adequacy only for the tested
  local-failure mode; it did not explain the Phase 1 `~5%` Dirichlet-Neumann
  circulation gap.
- 2026-07-23 review revised the phase: it now leads with (1) a read-only
  prior-evidence intake of `debug/dirichlet_solve/` and
  `debug/debug_bc_comparison.jl` (commit `be46062`, "Morino formulation
  velocity-only discrepancy"), and (2) a Gate A0 formulation sweep
  (`VelocityThroughSources` vs `TraceCorrected` vs `GreenReconstruction` on an
  unchanged capped/Dirichlet body) — `TraceCorrected`'s affine Kutta
  correction `γ = μ_upper − μ_lower − c` is a broad-span circulation shift
  matching the observed gap shape, and is the current leading hypothesis.
- Then the simple generated mesh oracle runs as a default three-point
  refinement ladder (not gated; Phase 1 showed Neumann drifting −1.314%
  toward converged Dirichlet under refinement) with the cheap Gate A
  field-circulation test: exterior velocity equivalence, loop-integral
  circulation (with a two-radii/two-quadrature self-check), and Kutta map
  consistency.
- Oracle case list revised: open+Dirichlet dropped (ill-posed Morino BC on an
  open surface); capped+Neumann with drop-one-panel regularization
  (`debug/debug_bc_comparison.jl:27–34`) promoted to the designated
  topology-versus-formulation tiebreaker.
- Naming fact recorded: one `Backslash` solver dispatches on the body `DBC`
  type parameter; `steady!` does not yet accept `formulation=` (only
  `simulate!` does), so Gate A0 needs minor wiring or a one-step frozen-wake
  `simulate!`.
- Gate A decides whether the gap is physical/formulation, topology/closure,
  refinement-unresolved, or likely Dirichlet representation/extraction before
  running expensive follow-ups.
- Thickness/TE closure sweeps, exterior-invisible mix modes, Green identity
  closure, wake/Kutta decomposition, forcing-versus-operator attribution, and
  DJI bridge cases are conditional rather than automatic.
- A DJI 2x2 topology/formulation bridge is conditional, not automatic.
- The geometry-only aerodynamic audit is intentionally excluded from this
  phase.
- No scripts, meshes, or solves have been run for Phase 2b yet.
- Phase 3 should remain blocked until Phase 2b is complete or Ryan explicitly
  cancels this inserted phase.

Phase 2 handoff: capped/Dirichlet new-mesh refinement changed integrated
circulation by only `+0.255%`, but uncapped/Neumann circulation was
approximately `5–7%` higher than capped/Dirichlet. Phase 2 found no sharp-TE
operator-fragility trigger, no material sensitivity growth, no TE-local
off-collocation outlier, and no panel-kernel-offset sensitivity.

## Working records

### Implementation and validation

| Date | Command or file | Purpose/result |
|---|---|---|
| 2026-07-23 | `plans/dji_convergence_20260722/phase_02b_formulation_attribution.md` | Added Dirichlet mix/nonuniqueness diagnostics to improve attribution of the Dirichlet-Neumann integrated-circulation gap. |
| 2026-07-23 | `plans/dji_convergence_20260722/phase_02b_formulation_attribution.md` | Revised Phase 2b into a gated diagnostic tree: Gate A field-circulation outputs first; thickness, mix, Green, wake/Kutta, forcing/operator, and DJI bridge follow-ups now conditional. |
| 2026-07-23 | `plans/dji_convergence_20260722/phase_02b_formulation_attribution.md` | Review pass: added prior-evidence intake (`debug/dirichlet_solve/`, `debug/debug_bc_comparison.jl`) and Gate A0 formulation sweep (`TraceCorrected`/`GreenReconstruction`); dropped ill-posed open+Dirichlet; promoted capped+Neumann drop-one-panel tiebreaker; made 3-point oracle refinement ladder default; specified loop-integral procedure + self-check; recorded constant-gauge cancellation and `Backslash`/`DBC` dispatch facts. |

### Oracle cases

| Case | Mesh/topology | Formulation | Refinement | Integrated circulation | Decision note |
|---|---|---|---:|---:|---|
| oracle | capped+Dirichlet (SD) | VTS, semi-inf wake | coarse/med/fine | −60.40/−63.71/−66.38 | baseline |
| oracle | open+Neumann (VR) | plain solve, semi-inf wake | coarse/med/fine | −61.61/−64.62/−66.90 | topology+BC ref (trustworthy) |
| oracle | capped+Neumann drop1 (VR) | plain solve, semi-inf wake | coarse/med/fine | −62.31/−65.30/−67.56 | pure-BC (mid-drop; drop-sensitive) |
| Gate A0 | capped+Dirichlet (SD) | VTS/TraceCorrected/Green, flat wake | ~medium | −32.01…−32.25 | formulation channel only 0.74% |

### Cross-residual audit

| Case | Primary residual checked | Spatial concentration | Relation to `Delta Gamma` | Decision note |
|---|---|---|---|---|
| all oracle bodies, all rungs | Gate A extraction max\|loopcirc−Γ\|/\|Γ\| | uniform ~3e-5 | ≪ the ~2.5% gap | faithful representation → real field difference, not extraction/mix |
| capped Neumann drop1 | drop-panel sensitivity | mid clean; alt panel-63 swung −47.6% | — | closed VR near-singular; open+Neumann is the reference |

### Wake/Kutta decomposition

| Case | Body-only contribution | Attached-wake/Kutta contribution | Closure/root/tip contribution | Decision note |
|---|---:|---:|---:|---|

### Conditional DJI bridge

| Case | Mesh source | Topology | Formulation | Integrated circulation | Attribution note |
|---|---|---|---|---:|---|

### Decisions and next-phase handoff

**Phase 2b attribution (reviewed, exit-gate answer): MIXED, quantified.**

1. **Dominant mechanism — chordwise section-resolution discretization difference.**
   The lower Dirichlet integrated circulation vs Neumann is a numerically-unresolved
   difference between the Morino source+doublet (Dirichlet) and vortex-ring (Neumann)
   discretizations of the airfoil section. It converges toward zero under **chordwise
   (`n_airfoil`) refinement** (BConly −3.18%→−2.49%→−1.78% on the joint ladder;
   chordwise-only closes ~0.9% while spanwise-only does not help). No plateau
   observed; residual at the finest tested rung ~1.8% BC-only / ~0.8% production-like.
2. **Convergence obtained:** the clean solves (capped Dirichlet, open Neumann)
   converge to a **0.11% gap at 120 chordwise panels/section** (open Neumann
   chordwise-converged by ~60; Dirichlet climbs up to it). Below the 1% and 0.25%
   criteria. The genuine physical cap effect on ∫Γ is small (≲0.1–0.2%); the earlier
   "~1% caps term" was the drop-one-panel regularization artifact (corrected).
3. **Excluded:** solve-formulation / Kutta-trace channel (Gate A0: 0.74%); Dirichlet
   representation/extraction/mix artifact (loop-integral velocity circulation =
   `Γ_TE` to 3e-5 for all bodies, all rungs); global thinness (thinner → *smaller*
   gap, opposite of what DJI's larger gap would require).

**Handoff to Phase 3 / Phase 5.** The Phase-1 DJI ~5% Dirichlet-Neumann gap is
predicted to be **chordwise-under-resolution**, not a physical or formulation defect.
Actionable lever for Phase 5 CT convergence: **refine DJI blade chordwise (section)
panel resolution**; expect the Dir-Neu gap and CT to converge together. Required
comparisons to carry forward: run the DJI capped/Dirichlet CT convergence with a
chordwise-resolution sweep, and (if a Neumann DJI case is built) confirm it drifts
toward Dirichlet under the same chordwise refinement. The DJI 2x2 bridge was NOT
needed — the oracle isolated the mechanism cleanly. Spanwise refinement and
timestep are not the Dir-Neu lever.

**Driver/asset state:** `examples/dji9443_formulation_attribution.jl` (modes
`smoke`, `formsweep`, `oracle`, `thickness`, `analyze`; `dji_bridge` stub, not
needed). `examples/pitching_wing.jl` gained `caps::Bool`/`thickness::Real`
(default behaviour unchanged; regression-tested). Artifacts under
`data/dji_convergence_20260722/phase_02b_formulation_attribution/`:
`formsweep_{summary,stations}.csv`, `oracle_summary.csv`,
`oracle_stations_{coarse,medium,fine}.csv`, `thickness_summary_medium.csv`,
`smoke_stations.csv`.

## Dated entries

Append exact commands, changed files and why, generated mesh descriptions,
artifact paths, failures/root causes, quantitative results, and reviewed
conclusions here. Keep the snapshot above current after every material action.

### 2026-07-24 — Ryan approved Phase 2b; handoff to Phase 2c

Ryan approved the Phase 2b attribution on 2026-07-24: the Dirichlet–Neumann
bound-circulation gap is a chordwise (n_airfoil) section-resolution discretization
difference (Morino source+doublet vs vortex-ring), with the oracle converging to a
0.11% gap at 120 chordwise panels/section.

Phase 2c was inserted before Phase 3 to verify the attribution on the actual DJI 9443
rotor mesh (chordwise ladder, fixed 30 spanwise, n_airfoil ∈ {81, 97, 121}). See
`plans/dji_convergence_20260722/phase_02c_dji_mesh_convergence.md` and its log. Phase 3
remains blocked until Ryan approves the Phase 2c decision.

### 2026-07-23 — Review pass on the staged Phase 2b plan (no solves run)

Exploration verified the plan's assumptions against the code and found:

- No velocity-based circulation check exists in Phases 1–2; Gate A is
  genuinely the missing test. Phase 1 extraction was `μ_upper − μ_lower` only
  (TE monitor and direct jump agreed to 0.0 relative; blade symmetry <6e-9).
- Commit `be46062` already added a "Morino formulation velocity-only
  discrepancy" diagnostic suite (`debug/dirichlet_solve/`: task1–5,
  `green_identity_analysis.md`, `unregularized_velocity_plan.md`; CL
  differences up to −3.4% between wake routes) plus
  `debug/debug_bc_comparison.jl` (working Dirichlet-vs-Neumann harness,
  drop-one-panel full-rank closed Neumann). The plan now intakes this first.
- `TraceCorrected` (`src/FLOWPanel_formulation.jl:112`) activates
  `γ = μ_upper − μ_lower − c` via `set_wake_correction!`
  (`src/FLOWPanel_liftingbody.jl:1628`) — a broad-span circulation shift
  matching the observed gap shape (integrated ≈ outboard percentage in
  Phase 1). Added as Gate A0, the cheapest decisive test.
  `TraceCorrected` requires `DirectBackend`; `steady!` lacks the
  `formulation=` kwarg (only `simulate!`, `src/FLOWPanel_simulate.jl:332`) —
  minor wiring noted as a driver task.
- Confound fix: open+Dirichlet is ill-posed (no interior region) and was
  dropped; capped+Neumann with drop-one-panel regularization is now the
  topology-vs-formulation tiebreaker.
- Refinement gating inverted for the cheap oracle: 3-point ladder now
  default; Phase 1 refinement asymmetry (Dirichlet +0.255% converged,
  Neumann −1.314% toward Dirichlet) recorded as evidence.
- Gauge scoping: a constant gauge mode cancels in `μ_upper − μ_lower`, so the
  `(I−B)` constant nullspace alone cannot explain the gap; mix-mode search
  restricted to non-constant directions. Off-body probes bypass `grad_mu`
  (induced-kernel path), so Gate A is clean of the triangulation-sensitive
  surface-velocity issue.
- Loop-integral procedure specified (station-plane circuit crossing the wake
  sheet once, offset probes, two-radii × two-quadrature self-check at 0.1%).
- Naming corrected: single `Backslash` solver dispatching on body `DBC`;
  Neumann = `VortexRing, DBC=false, watertight=false` (column 1);
  Dirichlet = `Union{ConstantSource,VortexRing}, DBC=true, watertight=true`
  (column 2, RHS `-self.potential`).

Phase 2b remains staged and awaiting Ryan's approval; no scripts, meshes, or
solves have been run.
