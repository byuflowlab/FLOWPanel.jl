# 028 — Non-uniform filter-width error in the SFS model

**Opened:** 2026-08-31 (Ryan). **Status:** staged; theory and measurement plan only.

**Item-level approvals:** Technical [ ]; clear-context [ ]; user [ ]

## RESET BRIEF

**Question:** how much error is introduced when the SFS derivation assumes a
spatially uniform filter width, while VPM particles carry non-uniform core sizes
$\sigma_p$, and which mitigation is accurate enough without defeating VPM locality?
**Status:** staged on 2026-08-31; no code, tests, or runs have been performed.
**Standing distinction:** source-attached particle cores $\sigma_q$ are not automatically
equivalent to one normalized observation-point filter $\sigma(\boldsymbol{x},t)$; Phase 0
must identify the actual FLOWVPM operation before assigning a commutation correction.
**Primary error:** variable-width filtering does not commute with spatial or temporal
differentiation. It can add a missing term to $\partial_jT_{ij}$ and can make the
filtered vorticity non-solenoidal even when $\partial_j\omega_j=0$.
**Relationship:** item 020 owns the proposed sigma-aware subgrid closure; item 028
isolates the variable-filter-width consistency error and may amend 020's equations or
recommend decoupling the LES filter from the particle core.
**Next action:** Phase 0 — audit the exact filter/core semantics used by the SFS and VPM
operators, then freeze the commutator definition and manufactured tests before examining
production snapshots.
**Gate:** do not implement a mitigation until Phases 0–2 quantify the error relative to
the modeled SFS term, resolved stretching, and observed tail-runaway events.

## Objective

Quantify the error caused by non-uniform particle/filter widths in the current SFS
formulation, determine whether it is negligible, localized, or dynamically important,
and select the least intrusive mitigation that restores the required consistency.

The investigation must distinguish three effects that are easy to conflate:

1. a true variable-width LES-filter commutation error;
2. quadrature and partition-of-unity error from irregular, unequally sized particles; and
3. the use of source-attached regularization cores in the VPM reconstruction, which need
   not define a single observation-point filter at all.

## Why this matters

The usual filtered vorticity derivation assumes

$$
\begin{align}
\partial_j\overline{f}
&= \overline{\partial_j f} \,.
\end{align}
$$

That identity holds for a homogeneous convolution kernel, but particle core growth,
stretching, splitting, merging, and remeshing can create a field of unequal $\sigma_p$.
If the particle core is interpreted as the LES filter, the uniform-filter assumption may
then omit terms at precisely the sharp core-size transitions and small-$\sigma$ tails that
matter to SFS activation and stability.

Item 020 asks whether an additional sigma-aware closure is needed. Before attributing a
measured discrepancy to missing subfilter physics, item 028 must determine whether part of
the discrepancy is instead a variable-filter commutation or reconstruction error.

## Theory seed: the missing commutator

For an observation-point filter whose width varies in space and time,

$$
\begin{align}
\overline{f}(\boldsymbol{x},t)
&= \int_V
   G_{\sigma(\boldsymbol{x},t)}(\boldsymbol{x}-\boldsymbol{y})\,
   f(\boldsymbol{y},t)\,
   \mathrm{d}\boldsymbol{y} \,.
\end{align}
$$

Define the spatial and temporal commutators by

$$
\begin{align}
\mathcal{C}_j(f)
&\equiv \partial_j\overline{f}
 - \overline{\partial_j f} \, , \\
\mathcal{C}_t(f)
&\equiv \partial_t\overline{f}
 - \overline{\partial_t f} \,.
\end{align}
$$

For a smooth, normalized kernel in an unbounded or periodic domain,

$$
\begin{align}
\mathcal{C}_j(f)
&= \left(\partial_j\sigma\right)
   \int_V
   \frac{\partial G_\sigma}{\partial\sigma}
   (\boldsymbol{x}-\boldsymbol{y})\,
   f(\boldsymbol{y})\,
   \mathrm{d}\boldsymbol{y} \\
&= \left(\partial_j\sigma\right)
   \frac{\partial\overline{f}}{\partial\sigma} \, , \\
\mathcal{C}_t(f)
&= \left(\partial_t\sigma\right)
   \frac{\partial\overline{f}}{\partial\sigma} \,.
\end{align}
$$

Consequently, the earlier uniform-filter identity for the SFS tensor becomes

$$
\begin{align}
\partial_j T_{ij}
&= \overline{\omega_j\,\partial_j u_i}
 - \overline{\omega_j}\,\partial_j\overline{u_i}
 + \mathcal{C}_j(u_i\omega_j)
 - \overline{u_i}\,\mathcal{C}_j(\omega_j) \,.
\end{align}
$$

The filtered vorticity also need not remain solenoidal:

$$
\begin{align}
\partial_j\overline{\omega_j}
&= \mathcal{C}_j(\omega_j) \,.
\end{align}
$$

These equations are a theory seed, not yet the FLOWVPM correction. They apply only if the
actual particle operation can be represented by a normalized observation-point filter.
For a Gaussian filter, the identity

$$
\begin{align}
\frac{\partial\overline{f}}{\partial\sigma}
&= \sigma\,\nabla^2\overline{f}
\end{align}
$$

provides a useful analytic reference, giving

$$
\begin{align}
\mathcal{C}_j(f)
&= \sigma\left(\partial_j\sigma\right)
   \nabla^2\overline{f} \,.
\end{align}
$$

The production VPM reconstruction instead has source-attached widths,

$$
\begin{align}
\omega_{h,j}(\boldsymbol{x})
&= \sum_{q=1}^{N}
   \Gamma_{q,j}\,
   \zeta_{\sigma_q}(\boldsymbol{x}-\boldsymbol{X}_q) \,.
\end{align}
$$

Here $\sigma_q$ is constant with respect to the evaluation coordinate
$\boldsymbol{x}$ for each summand. Differentiating this sum does not directly produce an
observation-point $\nabla\sigma$ term. Any correction must therefore follow from the
specific filter interpretation, particle-to-field reconstruction, and evolution law—not
from inserting $\sigma_p$ into the continuum commutator formula by analogy.

## Questions to answer

1. What operation is called “filtering” in the current SFS implementation: the particle
   basis itself, a test filter, a neighbor average, or a combination?
2. Is its width attached to the source particle, target particle, pair, or reconstructed
   spatial field?
3. Does the discrete operator reproduce constants and low-order moments when neighboring
   particles have unequal $\sigma$ and volume?
4. What are the spatial, material-time, and divergence commutators of the actual discrete
   operator?
5. How large are those terms relative to resolved vortex stretching, the modeled SFS
   contribution, viscous core spreading, and the residual of the filtered equation?
6. Are the largest errors concentrated in the small-$\sigma$ tail, shedding/handoff zone,
   splitting/merging transitions, or wake regions with poor kernel overlap?
7. Do error bursts precede or merely accompany the instability signatures studied in
   items 019 and 020?

## Phase 0 — operator and semantics audit

No model changes. Trace the current implementation and write down one discrete equation for
each relevant operation:

- particle vorticity reconstruction;
- velocity and velocity-gradient evaluation;
- SFS test filtering and coefficient evaluation;
- sigma evolution under stretching and diffusion;
- particle splitting, merging, and remeshing; and
- any normalization, volume weighting, or pair symmetrization.

The Phase-0 deliverable is a table stating, for every operator, which width is used and
whether it is source-, target-, pair-, or observation-point attached. It must decide whether
the continuum variable-filter commutator applies directly, requires a discrete analogue, or
is inapplicable because regularization and LES filtering are separate operations.

## Phase 1 — manufactured-field quantification

Construct smooth, divergence-free manufactured vorticity fields with known velocity and
prescribe width fields ranging from uniform to sharply graded. Include at least:

- a linear width gradient;
- a smooth localized width transition;
- a particle cloud with the same mean $\sigma$ but randomized neighbor-to-neighbor widths;
  and
- a source-attached-width reconstruction matched against an observation-point filter.

Evaluate derivatives in two ways: differentiate the filtered field directly, and filter
the analytic derivative. Verify convergence against high-accuracy quadrature and measure
the commutator as particle spacing is refined independently of the width variation.

For Gaussian filtering, first recover the analytic relation

$$
\begin{align}
\mathcal{C}_j(f)
&= \sigma\left(\partial_j\sigma\right)
   \nabla^2\overline{f}
\end{align}
$$

before using less ideal particle distributions. This separates implementation error from
the physical variable-filter term.

## Phase 2 — production-snapshot error budget

Use frozen snapshots from a representative stable interval and from a pre-instability
small-$\sigma$ tail interval. Do not begin with full reruns. For each particle, record:

$$
\begin{align}
r_{\mathrm{SFS}}
&= \frac{\left\|\boldsymbol{E}_{\mathrm{comm}}\right\|}
        {\left\|\boldsymbol{S}_{\mathrm{SFS}}\right\|+s_0} \, , \\
r_{\mathrm{stretch}}
&= \frac{\left\|\boldsymbol{E}_{\mathrm{comm}}\right\|}
        {\left\|\left(\overline{\boldsymbol{\omega}}\!\cdot\!\boldsymbol{\nabla}\right)
        \overline{\boldsymbol{u}}\right\|+s_0} \, , \\
e_{\nabla\cdot\omega}
&= \frac{\sigma\left|\boldsymbol{\nabla}\!\cdot
        \overline{\boldsymbol{\omega}}\right|}
        {\left\|\overline{\boldsymbol{\omega}}\right\|+\omega_0} \,.
\end{align}
$$

Here $\boldsymbol{E}_{\mathrm{comm}}$ is the complete missing spatial/material commutator
for the audited discrete operator, and $s_0$ and $\omega_0$ are preregistered floors used
only to prevent meaningless ratios in quiescent regions.

Report distributions, not only means: median, 90th, 99th, and maximum values, volume- or
circulation-weighted counterparts, and conditional statistics versus

$$
\begin{align}
\chi_{\sigma,p}
&= \max_{q\in\mathcal{N}_p}
   \frac{|\sigma_q-\sigma_p|}{\sigma_p} \, , \\
\chi_{h,p}
&= \frac{\sigma_p}{h_p} \,.
\end{align}
$$

Also compare the timing of commutator-error bursts with SFS activation, margin loss,
extreme strain, and sigma-tail collapse. Correlation alone is not causation; any apparent
lead-lag relation must be tested by a controlled mitigation A/B in Phase 3.

## Phase 3 — mitigation ladder

Test mitigations from least to most intrusive. Preserve a no-mitigation control and change
one mechanism at a time.

### A. Diagnose and bound width variation

- Add monitors for neighbor-width variation, overlap, partition of unity, and the discrete
  commutator residual.
- Impose a conservative limiter on neighbor-to-neighbor width ratios only if Phase 2 shows
  a resolved error threshold.
- Use splitting, merging, or remeshing to remove pathological width tails while preserving
  circulation and selected vorticity moments.

### B. Restore discrete consistency

- Renormalize kernel sums to reproduce constants and apply first-moment corrections where
  irregular particle volumes or supports break consistency.
- Use conservative pair-symmetric widths, such as a preregistered symmetric function
  $\sigma_{pq}=\sigma_{qp}$, if the audited operator currently treats a pair asymmetrically.
- Derive a discrete commutator correction by differentiating the actual normalized particle
  operator, including derivatives of both its kernel weights and normalization denominator.

### C. Add an explicit commutation model

- For an observation-point Gaussian filter, evaluate the exact
  $\sigma\,\partial_j\sigma\,\nabla^2\overline{f}$ correction where derivatives are
  sufficiently resolved.
- If direct derivatives are too noisy, test a scale-similarity or approximate-deconvolution
  estimate at a wider test-filter scale.
- Include the material-time contribution proportional to
  $D\sigma/Dt$; correcting only the spatial gradient would leave an inconsistent evolving
  filter.

### D. Decouple numerical core size from LES filter width

- Retain variable $\sigma_p$ for VPM regularization, but evaluate the SFS model with a
  separate, explicitly defined, locally uniform test filter.
- Choose that test-filter width independently of the smallest particle cores and verify
  that the resulting SFS model is insensitive to admissible changes in the particle-core
  distribution.
- This is the cleanest conceptual separation, but it adds filtering cost and may reduce
  locality or resolution if the explicit width must exceed the largest nearby core.

### E. Reformulate the closure for variable resolution

If the preceding approaches fail, derive the SFS equations from the chosen variable-width
filter from the outset, retaining spatial and temporal commutators and enforcing the
filtered-vorticity divergence constraint. This is the most principled and most invasive
option; it requires new conservation, consistency, and stability tests.

## Decision rules

Before Phase 2, preregister numerical thresholds for “negligible,” “model-relevant,” and
“dominant” based on ratios to the SFS and resolved-stretching terms. Do not choose the
thresholds after viewing the production distributions.

- If the commutator is small in both bulk and tails and does not change a controlled A/B,
  document it as a bounded modeling error and make no production change.
- If the error is primarily quadrature/normalization error, fix discrete consistency rather
  than adding a continuum commutation closure.
- If source-attached core variation is not mathematically the LES filter, decouple the
  concepts and amend item 020's derivation and terminology.
- If the error is localized to width discontinuities, prefer conservative width management
  or pair consistency over a global closure.
- If a resolved commutator is comparable to or larger than the existing SFS term and its
  mitigation materially changes the registered instability or fidelity metrics, promote
  the selected correction into item 020's validation matrix.

## Acceptance criteria

Technical completion requires all of the following:

- an explicit definition of the actual discrete filtering operator and its width semantics;
- a derived continuum or discrete commutator with sign and limiting-case checks;
- manufactured-field verification showing the expected convergence and uniform-width
  zero limit;
- a stable-versus-tail production error budget relative to SFS and resolved dynamics;
- at least one controlled mitigation A/B if the error is not negligible;
- conservation, divergence, and kernel-consistency checks for the selected mitigation; and
- a concise ruling that either bounds the error, amends item 020, or recommends a specific
  implementation phase with preregistered validation gates.

## Relationship to other items and notes

- `020_sigma_aware_subgrid_closure.md` owns the sigma-aware closure and its validation gate.
- `012_robust_resolution_overlap_management.md` owns broader overlap and resolution
  management; mitigation involving splitting, merging, or remeshing should be coordinated
  there rather than duplicated.
- `026_sigma_growth_particle_splitting/` owns particle splitting under sigma growth and is
  a likely source of both mitigation machinery and controlled width transitions.
- `les_vorticity_filter_stress_derivation.md` contains the seed uniform-filter derivation
  and the relationship between convolution filtering and VPM reconstruction.

## Log

- 2026-08-31: item opened at Ryan's request. Staged theory, manufactured-field,
  snapshot-error-budget, and mitigation phases. No implementation or run authorized.
