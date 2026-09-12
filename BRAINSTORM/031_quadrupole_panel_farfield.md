# 031 — Quadrupole (higher-order) far-field representation of panel kernels

**Opened:** 2026-09-11 (Ryan directive, on retiring 030 Phase 3b's dipole prototype).
**Status:** staged; exploration authorized, any production adoption Ryan-gated.

**Item-level approvals:** Technical [ ]; clear-context [ ]; user [ ]

## RESET BRIEF

**Question:** would a per-panel far-field multipole expansion of the Hess-Smith
panel kernels — quadrupole and possibly higher — ever be useful for assembling
influence blocks (near-field cache, and any future dense-block consumer), given
that 030 Phase 3b showed the monopole/dipole truncation is both too inaccurate
and, on the fixtures where the panel kernel actually dominates, rarely
applicable?

**What 030 Phase 3b established (measured 2026-09-10, retired 2026-09-11 —
prototype recoverable from FLOWPanel `030-block-assembly` commit `bc139c3`,
reverted by `0c444a4`):**

- A centroid point source σA + point dipole μA·n̂ reproduces `_induced` exactly
  in the far limit: φ = −A/4π(σ/r + μ n̂·r/r³), U = +∇φ, GT right-hand-rule
  normal; VortexRing ≡ constant-doublet panel (ratio 1.0000 measured).
  Attached-wake term must ALWAYS stay exact (semi-infinite filaments do not
  decay with panel distance).
- Error law (sphere Neumann cache, max block-normalized entry error):
  ≈ 0.14/η² for pairs beyond η·L (L = longest edge) — clean (L/r)². Note the
  centroid-centered expansion already has a vanishing first-moment correction
  (constant σ ⇒ ∫(x′−c)dS = 0), so the measured (L/r)² IS the quadrupole term:
  adding quadrupole buys ONE power, → O((L/r)³).
- Speed: point kernel 76 ns/pair vs 370–620 ns full; sphere cache build
  12.2×/5.8×/1.9× at η=2/4/8 (far fractions 97/88/50%).
- Blockers that retired the dipole version: (1) FMM-p=8-commensurate accuracy
  (~1e-8) needs η~4×10³ at (L/r)² — unreachable; quadrupole only moves this to
  η~(1e8)^(1/3)≈460, still unreachable; roughly (1/η)^p = 1e-8 needs p≈12 terms
  at η=4. (2) On the all-direct shedding fixture (Dirichlet diamond), core_size
  radius inflation puts essentially every direct pair NEAR: 1.2% far at η=2, 0%
  at η≥3 — the panels whose µs kernel motivates the work don't qualify.

**So the exploration must answer, in order:**

1. **Accuracy ladder:** derive closed-form centroid-centered moments for
   constant-strength flat triangles (source: monopole exact, quadrupole =
   second area moments ∫(x′ᵢx′ⱼ)dS — polynomial, exact, cheap per-column;
   doublet/VortexRing: dipole μA n̂ plus its second-moment correction). Measure
   the error-vs-η law with quadrupole (expect ~C/η³; measure C) and, if cheap,
   octupole. The 030 sweep harness (`eta_sweep.jl` pattern, error metrics,
   solve-level acceptance) is the template — resurrect from `bc139c3`.
2. **Utility criterion:** for each accuracy tier (1e-3 / 1e-5 / 1e-8 operator
   error), what η achieves it and what fraction of direct-list pairs qualify on
   (a) sphere-like bodies (cheap ~80 ns kernels) and (b) shedding/rotor bodies
   (µs kernels)? A tier is "useful" only if ≥~30% of pairs qualify on a fixture
   class whose kernel cost matters. 030's data says (b) fails at ANY η ≥ 3 —
   check whether that is intrinsic (MAC/leaf geometry) or an artifact of
   core_size radius inflation in the direct-list construction.
3. **The structural alternative (weigh before building anything big):**
   FastMultipole already forms panel multipole expansions at BRANCH level
   (`bodytomultipole.jl`) to arbitrary order — a per-pair per-panel expansion
   is a p≈2 branch of the same math. If the far pairs exist only because
   core_size inflation forces them into the direct list, the better lever may
   be radius-aware interaction-list construction (letting existing branch
   expansions absorb them at native p) rather than a new per-pair
   approximation. If instead a relaxed-accuracy niche is the target (e.g. FGS
   preconditioner blocks, where a 1e-4 operator perturbation is plausibly
   harmless — solve-level strength perturbation measured 1.9e-4 at η=4 with
   plain dipole), the dipole version may already suffice and quadrupole is
   unnecessary. Deliverable: a recommendation memo, not code adoption.

**Constraints carried over from 030:** exact assembly stays the default behind
a knob; rtol-1e-12 exactness tests untouched at the exact setting;
attached-wake and hessian rows exact; sign conventions calibrated numerically
against `_induced`, never re-derived on paper alone; prototype at the
assembly-hook level only (`FastMultipole.assemble_influence_block!` overload in
`src/FLOWPanel_abstractbody.jl`), never in `direct!` evaluation paths.

## Phases (proposed)

- **Phase 1 (theory + harness):** closed-form quadrupole moments for
  ConstantSource / ConstantDoublet / VortexRing triangles; numeric calibration
  script (extend the 030 pattern); measured error-vs-η law at dipole vs
  +quadrupole on sphere + diamond fixtures.
- **Phase 2 (utility study):** far-fraction and build-speedup maps per fixture
  class and accuracy tier; direct-list anatomy on a production-like rotor body
  (why are pairs direct: MAC, leaf, or core_size inflation?).
- **Phase 3 (recommendation):** adopt / park / pivot-to-interaction-list memo
  for Ryan. Production adoption is a separate Ryan gate regardless.

## Log

- 2026-09-11: opened on Ryan's directive while retiring 030 Phase 3 (route A)
  and Phase 3b (dipole prototype). All quantitative context above is from the
  030 session-3 sweep (see 030 Log, 2026-09-10 session-3 entry).
