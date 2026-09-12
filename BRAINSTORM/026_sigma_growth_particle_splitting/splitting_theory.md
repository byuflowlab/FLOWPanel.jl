# Resolution-preserving particle splitting — theory

This is the living theory document for the FLOWVPM splitting system
(`FLOWVPM.jl/src/FLOWVPM_resolution_split.jl`). It is maintained as a single
current draft — always the latest design, no change history. For the decision
record and campaign history see `particle_splitting_design.md`.

## 1. Conventions

- **Overlap** $\Phi \equiv \sigma/h$, where $h$ is inter-particle spacing.
  Shedding realizes this as $\sigma = h\,\Phi$
  (`FLOWPanel_wake.jl`: `sigma = dist·overlap/p_per_step`,
  `h = sigma/overlap`). Driver `OVERLAP` default 3.0; the 018 hover campaign
  convention is $\Phi = 2.75$. A particle born on a filament tiles length

  $$\ell_0 = \sigma_0/\Phi .$$

- **Per-particle split state** (`ResolutionSplitState`): reference radius
  $\sigma_0$ (stamped at creation / last split / last merge), sign-invariant
  averaged stretch axis + coherence weight, and two *attempted* (pre-clamp)
  $\Delta\sigma^2$ accumulators: `dvisc` $\ge 0$ (viscous spreading) and
  `drvpm` (signed net rVPM area evolution; compression $+$, elongation $-$).

- **Triggers** (per-mechanism growth fractions vs $\sigma_0$): mechanism $k$
  fires when $\sqrt{\sigma_0^2 + \Delta\sigma^2_k}/\sigma_0$ leaves
  $[\,1-f_{\rm elong},\ 1+f_k\,]$ on its side. The trigger IS the mechanism.
  $\sigma_{\min}/\sigma_{\max}$ are emission clamps, never triggers.

## 2. The tube picture: what a particle represents, and how it stretches

A vortex-tube particle represents a segment of tube of length $\ell$, core
radius $\sigma$, circulation $\Gamma_{\rm circ}$, and vector strength
$|\boldsymbol\Gamma| \approx \Gamma_{\rm circ}\,\ell$.

**The rVPM channel conserves $\sigma^2\ell$ exactly.** The rVPM closure
(Alvarez 2022; re-derived in `020_sigma_aware_subgrid_closure/phase_01_theory.md`)
evolves, with $Z = h_\sigma s$ the effective axial rate ($s$ = resolved
stretching rate along $\hat\Gamma$, $h_\Gamma = 2h_\sigma$ at the production
closure $(f,g)=(0,1/5)$):

$$\dot\sigma = -\sigma Z, \qquad
\frac{\mathrm d|\boldsymbol\Gamma|}{\mathrm dt}\bigg/|\boldsymbol\Gamma| = h_\Gamma s = 2Z .$$

Circulation is materially conserved (Kelvin), so $|\boldsymbol\Gamma| \propto
\ell$: the **represented segment length grows at rate $2Z$ while $\sigma^2$
decays at rate $2Z$**, hence $\sigma^2\ell = \mathrm{const}$ for the rVPM
channel, exactly. The length-stretch factor since the last split is therefore
read directly off the attempted accumulator:

$$\boxed{\ \lambda = \frac{\ell}{\ell_0}
   = \frac{\sigma_0^2}{\sigma_0^2 + \texttt{drvpm}}\ } \qquad
(\texttt{drvpm} < 0 \Rightarrow \lambda > 1),$$

and at the elongation trigger $\lambda^* = (1-f_{\rm elong})^{-2}$
(2.04 at $f_{\rm elong}=0.3$).

Two caveats, both acceptable for split sizing:

1. *Represented vs material stretch.* The material line element stretches as
   $e^{\int s\,dt}$; the represented tube as $e^{2\int Z\,dt} =
   e^{2h_\sigma\int s\,dt}$. At $h_\sigma = 1/5$ only $2/5$ of kinematic
   stretching becomes represented length (the rest lengthens
   $|\boldsymbol\Gamma|$ against fixed circulation bookkeeping). $\lambda$
   measures the represented tube — the correct quantity for
   re-discretization, since it is the represented support that must stay
   resolved.
2. *Attempted accumulation under clamps.* Each step's $\Delta\sigma^2$ is
   computed off the realized (clamped) $\sigma$, so `drvpm` chains
   linearized increments rather than following the free trajectory:
   floor-pinned at $\sigma_f$,
   $\lambda_{\rm att} = [1 - 2(\sigma_f/\sigma_0)^2\!\int\!Z\,dt]^{-1}$
   vs the free $e^{2\int Z\,dt}$ — equal to first order, diverging only
   over long pinned windows. Splits fire at small fractions (short
   windows), so the linearization is adequate.

## 3. Elongation kernel (in-line, adaptive child count)

Elongation ($\texttt{drvpm}<0$): the tube lengthens and thins; the
cross-section stays resolved, so the split re-discretizes the LENGTH. Children
lie on a line through the parent along the averaged stretch axis (Γ̂
fallback), each with $\sigma_c = \sigma_p$ (realized parent σ, cross-section
untouched), $\boldsymbol\Gamma_c = \boldsymbol\Gamma_p/m$, and `circulation`
unchanged — cutting a tube crosswise preserves each piece's circulation.

**Child count and spacing.** Choose a target child overlap $\Phi_t$ (knob
`elongate_overlap`; default = the shedding convention). The parent represents
$\lambda\,\ell_0$ of tube with $\ell_0 = \sigma_0/\Phi_t$; children of core
$\sigma_c$ should be spaced $s_t = \sigma_c/\Phi_t$ to retain overlap. Tiling
gives

$$m^\* = \frac{\lambda\,\ell_0}{s_t}
       = \lambda\,\frac{\sigma_0}{\sigma_c}
\;\;\xrightarrow{\ \sigma_c = (1-f)\sigma_0\ \text{(unclamped)}\ }\;
(1-f_{\rm elong})^{-3} = \lambda^{3/2},$$

implemented as $m = \mathrm{clamp}\!\big(\mathrm{round}(\lambda\,
\sigma_0/\sigma_c),\, 2,\, m_{\max}\big)$ with spacing
$s = \lambda\,\ell_0/m$ (children exactly tile the stretched length; after
integer rounding $s \approx s_t$). Offsets are centroid-symmetric,
$\big(k - \tfrac{m+1}{2}\big)s$ for $k = 1{:}m$. $m_{\max}$ (knob
`elongate_m_max`, default 4) caps capacity consumption per event and covers
$f_{\rm elong} \le 0.37$ without clamping; $f_{\rm elong}=0.5$ would want
$m^\*=8$.

**Properties** (exact up to integer rounding):

1. **Volume rule / merge inverse**: $m^\*\sigma_c^3 = \sigma_0^3$ — the
   overlap-matched tiling is precisely the inverse of merging's
   $\sigma = \sqrt[3]{\Sigma\sigma^3}$.
2. **No per-particle $\ell_0$ state needed**: a child's implied birth tile
   $\sigma_c/\Phi_t$ equals its true tile $\lambda\ell_0/m^\*$ identically —
   including when the parent is floor-pinned ($\sigma_c = \sigma_f$ enters
   both sides). The single knob $\Phi_t$ carries the discretization.
3. **Composition**: $m_1 m_2 = \lambda_1\tfrac{\sigma_0}{\sigma_1}\cdot
   \lambda_2\tfrac{\sigma_1}{\sigma_2} = \lambda_{\rm tot}\,
   \tfrac{\sigma_0}{\sigma_2}$ — repeated small splits reproduce the child
   count, σ, and tiled span of one large split.
4. **Conservation**: total $\boldsymbol\Gamma$, centroid, and linear impulse
   $\tfrac12\rho\,\Sigma\,\mathbf x_i\times\boldsymbol\Gamma_i$ are exact
   for any $m$ (symmetric offsets with equal $\boldsymbol\Gamma_c$ cancel
   pairwise; odd $m$ leaves the middle child at the parent position).
   Angular impulse error is quadratic in the offsets — the ± symmetry does
   not cancel it — bounded by $a^2|\boldsymbol\Gamma|/3$ with
   $a = (m{-}1)s/2$ the outermost offset, and vanishing when the split axis
   is parallel to $\boldsymbol\Gamma$ (the coherent-tube case the averaged
   axis targets).
5. Limits: $\lambda \to 1^+$ gives $m \to 2$ at spacing $\to \ell_0/2$
   (minimal re-discretization); the pathological
   $\sigma_0^2 + \texttt{drvpm} \le 0$ (attempted total collapse) saturates
   at $m = m_{\max}$.

With `elongate_overlap = NaN` the kernel falls back to the legacy fixed
2-child form: $\pm\,b$ with $b = \texttt{elongate\_offset\_ratio}\cdot
\sigma_p$ (default spacing $1.0\,\sigma_p$). The fixed form is dimensionally
arbitrary: its over/under-coverage ratio is $2\Phi(1-f_{\rm elong})^3$
(≈ 1.9× over-separated at $\Phi{=}2.75$, $f{=}0.3$; crossover at
$f \approx 0.43$), which is why the adaptive form replaces it.

## 4. Interaction with merging: the overlap gate

Splitting and merging are inverse operations, and their thresholds must not
overlap or they ping-pong. Worse than churn: split leaves σ unchanged while
merge returns $\sigma = m^{1/3}\sigma_c$, so a split→merge cycle is a
**σ-pump** (×$m^{1/3}$ per cycle, both split and merge state re-armed fresh
each time).

`merge_particles!` supports an overlap criterion natively: with
`sigma_relative=true` it merges a pair when
$\mathrm{dist} < r_{\rm merge}\cdot\sigma_{\min}$, i.e. when the pair overlap
$\Phi_{\rm pair} = \sigma_{\min}/\mathrm{dist}$ exceeds $1/r_{\rm merge}$.
**Stability condition**: the merge overlap threshold must exceed the child
target overlap,

$$\Phi_{\rm merge} > \Phi_t,$$

so split children (at $\Phi_t$) and the shed lattice (at $\Phi$) are never
merge candidates, at any σ (both sides scale with σ, so the guarantee
survives floor-pinning). Adopted values: $\Phi_t = 2.75$,
$\Phi_{\rm merge} = 3.5$ ($r_{\rm merge} = 1/3.5 \approx 0.286$, relative) —
tight structural pairs (shed pairs sit near $\sigma/h \approx 4.2$) still
merge; the filament lattice does not.

Reference: production hover merging currently uses an ABSOLUTE radius
$0.02R = 0.525\,\sigma^*$ (at the campaign $\sigma^* = 0.0381R$), which
exceeds the shed spacing $\sigma^*/\Phi = 0.0139R$ — i.e. it coarsens even
the fresh lattice, and would re-merge overlap-matched children immediately.
The overlap gate is therefore wired as an env knob (`MERGE_OVERLAP`,
default off) pending a production A/B: it merges strictly less than the
absolute radius does, so particle counts and cost rise.

## 5. Compression kernel (tri3)

Compression ($\texttt{drvpm}>0$): the tube shortens and fattens
($\lambda < 1$; cross-section area grows by $1/\lambda = (1+f_{\rm comp})^2$
at trigger). The fattened element is a bundle of thinner parallel filaments,
so the split re-discretizes the CROSS-SECTION: $m = 3$ children on an
equilateral triangle in the plane normal to the averaged stretch axis,
$\sigma_c = \sigma_p/\sqrt3$ (mass-per-length rule $3\sigma_c^2 =
\sigma_p^2$), $\boldsymbol\Gamma_c = \boldsymbol\Gamma_p/3$ parallel to the
parent, and `circulation`$/3$ — dividing a tube lengthwise splits the
vorticity flux among the filaments (contrast the crosswise elongation cut,
which preserves it).

**Count matching.** Re-discretizing the fattened cross-section into
birth-sized cores wants $m \approx \sigma_{\rm att}^2/\sigma_0^2 =
(1+f_{\rm comp})^2$ filaments: the fixed $m=3$ triangle is matched to
$f_{\rm comp} \approx \sqrt3 - 1 \approx 0.73$. The kernel keeps $m = 3$
regardless (an $m=2$ pair would impose an artificial transverse anisotropy
on what is physically an axisymmetric fattening; the triangle is the
smallest isotropic-in-the-plane arrangement), which effectively pairs the
mechanism with $f_{\rm comp} \approx 0.73$. Larger accumulated compression
could be matched by higher $m$ with a more complicated child arrangement
(e.g. ring + center, or two-ring patterns) — deliberately not pursued yet.

**Ring radius.** $a = \texttt{compress\_offset\_ratio}\cdot\sigma_p$
(default 0.6, spacing $1.8\,\sigma_c$; the transverse second-moment match is
$a = 1.155\,\sigma_p$ — knob-selectable, decision open pending the kernel-fit
study). Because $a$ scales with the realized $\sigma_p$, it already tracks
accumulated compression except when ceiling-pinned
($\sigma_p < \sigma_{\rm att}$); substituting
$\sigma_{\rm att} = \sqrt{\sigma_0^2 + \texttt{drvpm}}$ would restore the
attempted picture if that regime ever matters. Grow-side children at the
default spacing are comfortably outside both merge criteria under production
numbers.

## 6. Viscous kernel (tetra4)

Isotropic spreading ($\mathrm d\sigma^2/\mathrm dt = 2\nu$): 4 children on a
randomly oriented regular tetrahedron, $\sigma_c = \sigma_p\,4^{-1/3}$
(volume rule = merge inverse), vertex radius $a = 1.3503\,\sigma_p$ (per-axis
second-moment match: $\sigma_p^2 = \sigma_c^2 + a^2/3$),
$\boldsymbol\Gamma_c = \boldsymbol\Gamma_p/4$, `circulation`$/4$ (bundle
division, as in §5). Nearest-child spacing $\approx 3.5\,\sigma_c$ is
under-overlapped; the child geometry remains open pending a pointwise/$L_2$
kernel-fit and induced-field study (larger child sets on the table). The
mechanism has zero observed production events and stays evidence-gated.

## 7. Conservation and verification obligations

Non-negotiable per split event, any kernel, any $m$: total
$\boldsymbol\Gamma$, Γ-weighted centroid, and linear impulse exact; angular
impulse error bounded by $a^2|\boldsymbol\Gamma|/3$ ($a$ = outermost child
offset) and pinned by test for every kernel. Far-field velocity/strain equivalence and combined-support
$L_2$ residuals are regression-pinned per kernel geometry; any geometry
change re-derives the pins rather than loosening them. The adaptive-$m$
elongation kernel additionally pins the composition property (§3.3) and the
child-overlap target. Test suite:
`FLOWVPM.jl/test/runtests_resolution_split.jl`.
