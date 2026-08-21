# 020 Phase 1 — Theory: a σ-aware subgrid closure for the rVPM

**Status: SUPERSEDED IN PART (2026-08-12 Phase-2R audit).** Sections 1–8
remain the historical Phase-1 registration, but the corrections immediately
below govern all subsequent work. Author: agent.

### 2026-08-12 verification corrections (binding)

The Phase-2R independent audit found that several conclusions were stated more
strongly than the derivation supports:

1. The retained-$E_{\mathrm{adv}}$ calculation establishes the standard
   core-spreading route only under a locally uniform diffusivity, an isolated
   Gaussian blob, and the localized/no-interparticle-transfer approximation.
   Global conservation of $\sum_p\boldsymbol\Gamma_p$ does not uniquely force
   every $\dot{\boldsymbol\Gamma}_p=0$; conservative interparticle strength
   transfer is also admissible. Thus $\dot\sigma=\nu_{\mathrm{sgs}}/\sigma$ at
   fixed particle strength is a well-motivated candidate, not a unique closure.
2. Hereafter $\nu$ is molecular viscosity,
   $\nu_{\mathrm{sgs}}\ge0$ is the modeled addition, and
   $\nu_{\mathrm{eff}}=\nu+\nu_{\mathrm{sgs}}$. Expressions below that write
   `$\nu_t=\nu+\cdots$` use $\nu_t$ as the *total* effective viscosity; they
   must not then be added to $\nu$ a second time.
3. The fixed-filter form and $\kappa=1/\sqrt5$ are resolution-cutoff modeling
   choices, not consequences of the filtered equation or physical constants.
   The value $1/\sqrt5$ applies only when the closure uses positive aligned
   stretching $s_+=\max(\hat\Gamma^T D\hat\Gamma,0)$, for which $s_+=5Z$ in
   the live aligned case. A norm such as $\sqrt{2D:D}$ has a different
   normalization and therefore a different equilibrium constant.
4. The unconditional result in §4.5 is withdrawn. Positivity is provable for
   the corrected frozen-gradient geometric update, but bounded closed-loop
   stability is conditional on the scalar Kelvin-filament assumption
   $Z\sim\Gamma_v/(2\pi\sigma^2)$, bounded ambient strain, and simultaneous
   strainers not invalidating that bound. The failed registered
   $\Gamma_{\mathrm{implied}}$ ladder criterion prevents elevating this model
   to a general proof.
5. A pure core-size update at fixed position and strength preserves total
   vorticity and linear impulse. It does **not** preserve angular impulse:
   for an isotropic blob, angular impulse contains a contribution proportional
   to $\sigma^2\boldsymbol\Gamma$. The original §5.2/T4 angular-impulse claim
   is withdrawn.
6. The 1.9–3.3-decade result is a scenario estimate conditional on the assumed
   $r_c/c=0.01$–$0.05$ bracket and vortex-profile coefficient, not a validated
   DJI9443 physics target. Squire supports apparent viscosity proportional to
   circulation; the empirical coefficient is model-dependent over a broad
   literature range and does not determine $\kappa$.
7. The original `euler_exp` froze the initial stretching vector as a forcing.
   It was therefore not the frozen-gradient rVPM solution: under aligned strain
   its strength multiplier tends to $5/3$, rather than $e^{2\Delta tZ}$.
   Original Stage A and Stage B mechanistic conclusions are invalidated and
   must be replaced by the corrected Phase-2R integrator and rerun.

Every claim carries a confidence label per the 018 house style:
**CONFIRMED** = directly measured/verified in code or data, evidence cited ·
**DERIVED** = follows from stated equations by algebra reproduced here ·
**SUPPORTED** = consistent multi-source evidence · **HYPOTHESIS** = stated, untested.

Scope guard: this document is Phase 1 of
`BRAINSTORM/020_sigma_aware_subgrid_closure.md` — theory only. It derives and
surveys; it runs nothing and edits no source. Phases 2–3 protocols are
pre-registered in §7–§8 *before* any of their data is examined.

---

## 1. The live discrete system, stated exactly (CONFIRMED, code-cited)

### 1.1 Continuous governing equations (rVPM)

FLOWVPM (dev'd at `~/Dropbox/research/projects/FLOWVPM.jl`, v4.1.0) integrates
the reformulated VPM of Alvarez & Ning: the LES-filtered vorticity equation

$$
\frac{\mathrm{d}}{\mathrm{d}t}\overline{\boldsymbol\omega}
= (\overline{\boldsymbol\omega}\cdot\nabla)\overline{\mathbf u}
+ \nu\nabla^2\overline{\boldsymbol\omega}
- \mathbf E_{\mathrm{adv}} - \mathbf E_{\mathrm{str}},
$$

discretized with particles $\overline{\boldsymbol\omega}(\mathbf x) \approx
\sum_p \boldsymbol\Gamma_p\, \zeta_{\sigma_p}(\mathbf x-\mathbf x_p)$, where the
filter kernel doubles as the particle basis and the filter width **is** the
core size $\sigma_p$, a dynamic per-particle quantity. The governing equations
with the free parameters $(f,g)$ (`docs/src/rVPM.md`, Eqs. 6–8 of the excerpt;
`FLOWVPM_formulation.jl:31-37`):

$$
\begin{align}
\frac{\mathrm{d}\boldsymbol\Gamma_p}{\mathrm{d}t}
&= (\boldsymbol\Gamma_p\cdot\nabla)\overline{\mathbf u}
 - \frac{g+f}{\tfrac13+f}\Big[\big((\boldsymbol\Gamma_p\cdot\nabla)\overline{\mathbf u}\big)\cdot\hat{\boldsymbol\Gamma}_p\Big]\hat{\boldsymbol\Gamma}_p
 - \frac{C_d}{\zeta_{\sigma_p}(0)}\Big[\mathbf E_{\mathrm{str}} - \frac{f}{\tfrac13+f}\big(\mathbf E_{\mathrm{str}}\cdot\hat{\boldsymbol\Gamma}_p\big)\hat{\boldsymbol\Gamma}_p\Big],
\\[4pt]
\frac{\mathrm{d}\sigma_p}{\mathrm{d}t}
&= -\frac{g+f}{1+3f}\,\frac{\sigma_p}{\Vert\boldsymbol\Gamma_p\Vert}\Big[\big((\boldsymbol\Gamma_p\cdot\nabla)\overline{\mathbf u}\big)\cdot\hat{\boldsymbol\Gamma}_p\Big]
 + \frac{f}{1+3f}\,\frac{\sigma_p}{\Vert\boldsymbol\Gamma_p\Vert}\,\frac{C_d}{\zeta_{\sigma_p}(0)}\,\mathbf E_{\mathrm{str}}\cdot\hat{\boldsymbol\Gamma}_p .
\end{align}
$$

The live default is $f=0,\ g=1/5$ (`FLOWVPM.jl:99,133-134`, `formulation_rVPM`),
chosen so a spherical fluid element conserves mass and momentum. Available
alternatives (`FLOWVPM.jl:97-116`): classic VPM (no σ equation), cVPM $(0,0)$,
tube-continuity $(1/2,0)$, tube-momentum $(1/4,1/4)$, sphere-momentum
$(0,1/5{+}10^{-8})$. Note the SFS term enters the σ equation **only** through
the $f/(1+3f)$ prefactor — identically zero at $f=0$. This is the
**channel gap**: in the live formulation no SFS model, at any coefficient, can
write σ (CONFIRMED: `FLOWVPM_timeintegration.jl:147`; the SFS test-filter
σ-scalings are exactly undone, `FLOWVPM_subfilterscale.jl:613/674, 890/904`;
018 corpse: σ marched to the laminar floor with SFS on).

### 1.2 Notation for the map analysis

Define the **axial resolved stretching rate** (units 1/s) and the live Z:

$$
s \equiv \frac{\big[(\boldsymbol\Gamma\cdot\nabla)\overline{\mathbf u}\big]\cdot\hat{\boldsymbol\Gamma}}{\Vert\boldsymbol\Gamma\Vert},
\qquad
Z \equiv \frac{g+f}{1+3f}\,s \;\underset{f=0,\,g=1/5}{=}\; \frac{s}{5}.
$$

$Z$ is the code's `MM4` (`FLOWVPM_timeintegration.jl:143-151`) and the exact
per-step readout is $\Delta\sigma/\sigma = -\Delta t\,Z$ (the 019 `max_dtZ`
monitor column). The transposed-scheme switch changes which gradient enters
$s$, not this structure.

### 1.3 The Euler-step map, all three channels (CONFIRMED at `FLOWVPM_timeintegration.jl:103-173`)

Per particle, per step (SFS terms omitted here; they add to Γ only at $f=0$):

1. $\mathbf x \leftarrow \mathbf x + \Delta t\,(\overline{\mathbf u}+\mathbf u_\infty)$;
2. $\boldsymbol\Gamma \leftarrow \boldsymbol\Gamma + \Delta t\,\big(\mathbf S - 3Z\,\boldsymbol\Gamma\big)$ with $\mathbf S \equiv (\boldsymbol\Gamma\cdot\nabla)\overline{\mathbf u}$ (`:153-158`);
3. $\sigma \leftarrow \sigma\,(1-\Delta t Z)$ (`:160-161`);
4. viscous (CoreSpreading only): $\sigma \leftarrow \sqrt{\sigma^2 + 2\nu\Delta t}$ (`FLOWVPM_viscous.jl:161`), the exact integral of $\dot\sigma = \nu/\sigma$, applied post-stretch. `Inviscid()` is a no-op (`:58`).

Splitting Γ into axial and transverse parts,
$\boldsymbol\Gamma = \Gamma_\parallel\hat{\boldsymbol\Gamma} + \boldsymbol\Gamma_\perp$
(with $\Gamma_\perp = 0$ instantaneously but sourced by $\mathbf S_\perp$), and
using $\mathbf S\cdot\hat{\boldsymbol\Gamma} = s\,\Vert\boldsymbol\Gamma\Vert = 5Z\Vert\boldsymbol\Gamma\Vert$:

$$
\begin{align}
\Gamma_\parallel &\leftarrow \Gamma_\parallel\,\big(1 + \Delta t\,(s - 3Z)\big) = \Gamma_\parallel\,(1 + 2\Delta t Z), \\
\boldsymbol\Gamma_\perp &\leftarrow (1 - 3\Delta t Z)\,\boldsymbol\Gamma_\perp + \Delta t\,\mathbf S_\perp, \\
\sigma^2 = y &\leftarrow (1-\Delta t Z)^2\,y + 2\nu\Delta t .
\end{align}
$$

**Instability thresholds of the frozen-Z map (DERIVED):**

| channel | multiplier | $|{\cdot}|>1$ when | character |
| --- | --- | --- | --- |
| $\Gamma_\parallel$ | $1+2\Delta t Z$ | always ($Z>0$) | **physical** growth; Euler under-shoots the exact $e^{2\Delta t Z}$ |
| $\boldsymbol\Gamma_\perp$ | $1-3\Delta t Z$ | $\Delta t Z > 2/3$ | **spurious** sign-flipping amplification |
| $\sigma$ (linear) | $1-\Delta t Z$ | $\Delta t Z > 2$ | spurious flip at $\Delta t Z>1$ (negative σ, inviscid), geometric beyond 2 |

**Correction to the discussion record (DERIVED, load-bearing).** Item 3 of the
020 discussion record states the local dynamics as
$\sigma \leftarrow \sigma e^{-\Delta t Z}$, $\Gamma \leftarrow \Gamma e^{-3\Delta t Z}$
and concludes the velocity scale $\Gamma/\sigma^2 \sim e^{-Z t}$ *decays*.
That caricature drops the stretching source $\mathbf S$: the actual code
retains it (`:156-158`), the axial multiplier is $(1+2\Delta t Z)$, and the
continuous local dynamics under sustained aligned stretching gives

$$
\frac{\mathrm{d}}{\mathrm{d}t}\ln\frac{\Vert\boldsymbol\Gamma\Vert}{\sigma^2}
= 2Z + 2Z = 4Z \;>\;0 :
$$

the velocity scale **grows** (exponentially in time) — vortex-stretching
intensification, which is physical (Kelvin: tube circulation fixed, core
thins, peak swirl rises), and is saturated in reality only by (turbulent)
viscosity at the Burgers scale. So the claim "the blow-up is discrete, not
physical" must be **split**: the $\Delta t Z > 2/3$ transverse flip and the
$\Delta t Z>1$ negative-σ flip are purely discrete artifacts an exact local
integrator removes (§3f); the underlying exponential intensification at rate
$\sim 4Z$ with $Z \propto 1/\sigma^2$ through the field coupling is
*continuous physics that molecular ν cannot stop at these $Re_\Gamma$* — that
part only a closure (or DNS-resolution) addresses. This sharpened statement
replaces discussion-record item 3 and is testable in Phase 2 Stage A/B (§7).

The two regimes measured in 018 (`018_.../sigma_blowup_mechanism.md` §2–§3)
are the y-map rows of this table plus the CoreSpreading offset: ambient
$0<\Delta t Z<2$ has the stable fixed point
$\sigma^* = \sqrt{2\nu/[Z(2-\Delta t Z)]} \to \sigma_{eq}=\sqrt{\nu/Z}$; the
measured ignition route is Γ-side amplification at floor-pinned σ
($|\Gamma|/\sigma^2$ ×3.5/step at the $\sqrt{2\nu\Delta t}$ floor), consistent
with the $\Gamma_\perp$/$\Gamma_\parallel$ rows at raw $\mathrm{d}t\!\cdot\!Z$
p95 > 2.4.

### 1.4 What physics the σ channel is missing (SUPPORTED)

At the operating resolution the entire vortex core is one particle
($\sigma \approx$ core radius): there is no resolved field inside the core,
so core-interior turbulence is wholly subfilter. Its dominant effect is
**radial transport** — turbulent mixing diffuses core vorticity outward,
setting core growth (Squire 1965; Bhagwat & Leishman's generalized viscous
vortex model: effective viscosity $\nu_{\mathrm{eff}} = \delta\,\nu$ with
$\delta = 1 + a_1 Re_v$, $Re_v \equiv \Gamma_v/\nu$, $a_1 \approx 6.5\times10^{-5}$;
Lamb–Oseen growth $r_c = \sqrt{r_0^2 + 4\alpha_L\delta\nu\,t}$,
$\alpha_L = 1.25643$). In VPM variables that is a σ-channel effect:
$\dot\sigma = \nu_{\mathrm{eff}}/\sigma$. The live model has **no** such term
beyond molecular CoreSpreading, and the live SFS is structurally Γ-side
(§1.1). Honest scale note (DERIVED from the model constants): at this rotor's
$Re_v \approx 0.278/1.43\times10^{-5} \approx 1.9\times10^4$,
$\delta \approx 1 + 6.5\times10^{-5}\cdot 1.9\times10^4 \approx 2.3$ — the
Squire/BL correction here is a factor ~2, **not** "tens of ×" (that figure
belongs to full-scale rotors at $Re_v\sim10^6$). The discussion-record's
"measured tip-vortex ν_eff runs tens of × molecular" is thereby **scoped, not
refuted**: the physical subgrid viscosity available from the BL anchor at this
scale is modest, which §4 shows has structural consequences (no
physically-constanted eddy viscosity arrests filament self-thinning; the
cutoff must come from the resolution scale, not from δ).

---

## 2. The σ-equation subfilter term, derived (targeted extension)

Citations: Alvarez (2022) dissertation ("A", printed pages) and Alvarez–Ning
AIAA J ("B"). This section starts from A's filtered equations and collocation
machinery — cited, not reproduced — and re-derives only the σ equation with
the dropped terms retained.

### 2.1 How the particle equations arise: collocation, not moment matching (CONFIRMED, A/B)

A does **not** match moments. The route [A Eq. (1.8)–(1.14), pp. 17–18;
B (18)–(32)]: substitute
$\overline{\boldsymbol\omega} \approx \sum_q \boldsymbol\Gamma_q\zeta_{\sigma_q}$
into the *inviscid* filtered vorticity equation, expand the Lagrangian
derivative of the basis [A App. A.1, p. 158], and **collocate at
$\mathbf x = \mathbf x_p$**, using
$\tfrac{\partial\zeta_\sigma}{\partial t}(0) = -\tfrac{3}{\sigma}\dot\sigma\,\zeta_\sigma(0)$
[A (1.13)] to obtain the general governing equation [A (1.14)]

$$
\frac{\mathrm d\boldsymbol\Gamma_p}{\mathrm dt}
= (\boldsymbol\Gamma_p\cdot\nabla)\overline{\mathbf u}(\mathbf x_p)
+ 3\,\boldsymbol\Gamma_p\,\frac{\dot\sigma_p}{\sigma_p}
+ \frac{-M^0_p + M^1_p + M^2_p}{\zeta_{\sigma_p}(0)}
- \frac{\mathbf E_{\mathrm{adv}}(\mathbf x_p) + \mathbf E_{\mathrm{str}}(\mathbf x_p)}{\zeta_{\sigma_p}(0)} .
$$

One equation, two unknowns $(\boldsymbol\Gamma_p, \sigma_p)$ — underdetermined.
The neighbor terms $M^{0,1,2}_p$ are dropped by the **localized-vorticity
assumption** [A pp. 19–21], with the explicit remark that the neglected terms
"become part of the subfilter-scale contributions." The closure is the
one-parameter-family ansatz [A (2.17), p. 28]

$$
\frac{\mathrm d\sigma_p}{\mathrm dt}
= -f\,\frac{\sigma_p}{\Vert\boldsymbol\Gamma_p\Vert}\,\frac{\mathrm d\boldsymbol\Gamma_p}{\mathrm dt}\cdot\hat{\boldsymbol\Gamma}_p
\;-\; g\,\frac{\sigma_p}{\Vert\boldsymbol\Gamma_p\Vert}\,\big[(\boldsymbol\Gamma_p\cdot\nabla)\overline{\mathbf u}\big]\cdot\hat{\boldsymbol\Gamma}_p ,
$$

($f$: feedback from the *realized* strength change; $g$: direct coupling to
resolved stretching), whose back-substitution [A (2.18)–(2.20), p. 29] yields
the §1.1 equations with $h_\Gamma = \tfrac{1-3g}{1+3f}$,
$h_\sigma = \tfrac{g+f}{1+3f}$, and the consistency constraint
$h_\Gamma + 3h_\sigma = 1$: every $(f,g)$ member differs only in how
stretching is split between Γ lengthening and σ shrinking. $f{=}0,\,g{=}1/5$
follows from angular-momentum conservation of a spherical element
[A §2.2.1, p. 25]; Table 2.1 [A p. 30] catalogs the alternatives
(over-stretching for $h_\Gamma>1$, anti-stretching for $h_\sigma>1/3$).
Note in these variables the live axial rates are $h_\Gamma s = 2Z$ and
$-h_\sigma s = -Z$ — the §1.3 map multipliers, cross-checked.

**Also CONFIRMED and relevant: no stability analysis exists in A or B.**
B p. 4, verbatim: "Formal analysis will be needed in future work to better
understand the numerical properties of the rVPM." Stability is established
empirically (round jet, hover rotor; rVPM and SFS *both* required). The 018
§2 map analysis and §4 of this document are, to our knowledge, that missing
formal analysis for the σ/Γ subsystem.

### 2.2 The three doors through which σ-channel physics was lost (CONFIRMED, cited verbatim)

1. **Viscosity — removed by operator splitting *before* collocation**
   [A §1.3, p. 14; explicit at A p. 18 and B p. 5: core spreading "must not
   be included" in the inviscid $\dot\sigma$]. Viscous diffusion re-enters as
   the separate scheme $\dot\sigma = \nu/\sigma$ (exact Gaussian solution
   $\mathrm d\sigma^2/\mathrm dt = 2\nu$ [A (4.10), p. 55]; additive
   composition [A p. 61] — exactly the live code). **Molecular ν throughout;
   no turbulent viscosity is ever added.**
2. **$\mathbf E_{\mathrm{adv}}$ — carried symbolically, then zeroed in one
   sentence** [A §3.6, p. 40: "we also assume $E_{adv}=0$ since our interest
   is in testing and validating our SFS vortex-stretching model"; B §V.A:
   "it is recommended that $E_{adv}$ is added in other cases, where advection
   is the physical mechanism for the development of turbulence"].
3. **$f=0$ closes the only SFS route into σ.** In the ansatz, SFS reaches σ
   only through the $f$-weighted feedback of
   $\dot{\boldsymbol\Gamma}_p\cdot\hat{\boldsymbol\Gamma}_p$, giving the
   $\tfrac{f}{1+3f}$ prefactor [A (2.19)] — so the prefactor **is derived**
   (from the ansatz structure), and its vanishing at $f=0$ is a **silent
   structural consequence, not an argued approximation**: no sentence in A or
   B analyzes it. (Plus, zeroth: the M-terms — the un-modeled discrete
   analogs of $E_{\mathrm{adv}}/E_{\mathrm{str}}$ [A §3.2, p. 35] — were
   already absorbed "into the subfilter-scale contributions.")

**Verdict on discussion-record item 1 (derivation-level): SUPPORTED and
sharpened.** The missing channel is not merely "SFS can't write σ at $f=0$";
it is that the *physics* of turbulent core transport — subfilter velocity
fluctuations advecting core vorticity radially — lives in
$\mathbf E_{\mathrm{adv}}$ (advective subfilter flux), which was zeroed for
scope, in a flow class (concentrated strained cores) that B's own caveat says
calls for it. The channels are complementary by construction:
$\mathbf E_{\mathrm{str}}$ (modeled, Γ-side) is enstrophy transfer by
stretching; $\mathbf E_{\mathrm{adv}}$ (zeroed) is transport. A σ-channel
closure therefore does not compete with the existing SFS model — it fills the
term the existing model explicitly set aside.

### 2.3 Re-derivation: retaining $\mathbf E_{\mathrm{adv}}$ under a gradient-diffusion closure (DERIVED — the core Phase-1 result)

Model the subfilter advective flux by gradient diffusion (the standard
eddy-diffusivity closure of $\overline{u_j'\omega_i'}$):

$$
\mathbf E_{\mathrm{adv}} \;\approx\; -\,\nabla\cdot\big(\nu_t\,\nabla\overline{\boldsymbol\omega}\big)
\;=\; -\,\nu_t\,\nabla^2\overline{\boldsymbol\omega} \quad (\nu_t \text{ locally uniform}),
$$

with $\nu_t \ge 0$ a subfilter diffusivity whose *scale* is deliberately left
open here (that choice is §3's subject; the derivation constrains the form,
not the constant). Retain this term through the collocation. For a single
blob (localized-vorticity, as in A's single-blob route [A (1.17), p. 20]),
$\overline{\boldsymbol\omega}(\mathbf x) = \boldsymbol\Gamma_p\zeta_{\sigma_p}(\mathbf x-\mathbf x_p)$
and, for the Gaussian basis,

$$
\nabla^2\zeta_{\sigma}(\mathbf 0) = -\frac{3}{\sigma^2}\,\zeta_\sigma(\mathbf 0),
$$

so the retained term contributes, at the collocation point,

$$
-\frac{\mathbf E_{\mathrm{adv}}(\mathbf x_p)}{\zeta_{\sigma_p}(0)}
= \frac{\nu_t\,\nabla^2\overline{\boldsymbol\omega}(\mathbf x_p)}{\zeta_{\sigma_p}(0)}
= -\,\frac{3\,\nu_t}{\sigma_p^2}\,\boldsymbol\Gamma_p .
$$

The collocated balance [A (1.14) with only this term and the
$3\boldsymbol\Gamma_p\dot\sigma_p/\sigma_p$ term active] is again
underdetermined — the diffusive decay of the centerline vorticity can be
booked as strength decay at fixed σ, as core growth at fixed Γ, or any mix:

$$
\frac{\mathrm d\boldsymbol\Gamma_p}{\mathrm dt}\bigg|_{\mathrm{adv}}
- 3\,\boldsymbol\Gamma_p\,\frac{\dot\sigma_p}{\sigma_p}\bigg|_{\mathrm{adv}}
= -\,\frac{3\,\nu_t}{\sigma_p^2}\,\boldsymbol\Gamma_p .
$$

*(Sign bookkeeping: in A (1.14) the $+3\Gamma\dot\sigma/\sigma$ term sits on
the RHS of the Γ equation; moving the split to the left as here makes the
routing explicit.)* For the isolated-blob core-spreading realization, impose
fixed particle strength. This is compatible with the conservation law
**total represented vorticity** $\int\overline{\boldsymbol\omega}\,\mathrm d\mathbf x
= \sum_p\boldsymbol\Gamma_p$ must be invariant under pure diffusion (diffusion
transports, it does not create or destroy circulation). It excludes local
strength decay for the isolated blob, although global conservation also admits
conservative interparticle transfer. The local realization gives

$$
\boxed{\ \frac{\mathrm d\boldsymbol\Gamma_p}{\mathrm dt}\bigg|_{\mathrm{adv}} = 0,
\qquad
\frac{\mathrm d\sigma_p}{\mathrm dt}\bigg|_{\mathrm{adv}} = \frac{\nu_t}{\sigma_p}\ }
$$

— which also reproduces, exactly, the analytic diffusion of a Gaussian blob
($\sigma^2 \mathrel{+}= 2\nu_t t$ at fixed Γ; the same argument that makes
CoreSpreading the exact molecular scheme [A (4.10)]). Cross-check at the
collocation point: $\overline{\omega}(\mathbf x_p) = \Gamma\zeta_\sigma(0)$
with $\dot\sigma = \nu_t/\sigma$ gives
$\dot{\overline\omega}(\mathbf x_p) = -\tfrac{3\nu_t}{\sigma^2}\overline\omega(\mathbf x_p)
= \nu_t\nabla^2\overline\omega(\mathbf x_p)$ ✓.

**Corrected result (SUPPORTED CANDIDATE).** Retaining the σ-equation subfilter
term under the standard gradient-diffusion closure, collocation, and the
isolated-blob/fixed-strength realization produces

$$
\frac{\mathrm d\sigma_p}{\mathrm dt} = -\,h_\sigma\, s\,\sigma_p \;+\; \frac{\nu + \nu_{\mathrm{sgs}}(\mathbf x_p)}{\sigma_p},
\qquad
\frac{\mathrm d\boldsymbol\Gamma_p}{\mathrm dt} = \text{(unchanged inviscid + }E_{\mathrm{str}}\text{ terms)} .
$$

That is: **the resulting candidate is a per-particle turbulent viscosity
entering the existing CoreSpreading channel** — the
$\dot\sigma = \nu_{\mathrm{eff}}/\sigma$ form assumed by candidates 3a–3c,
where $\nu_{\mathrm{eff}}=\nu+\nu_{\mathrm{sgs}}$ — implemented at
the same operator-splitting stage as molecular ν (answering placeholder Q3:
$\nu \to \nu + \nu_{\mathrm{sgs}}$ preserves the derivation's structure term-by-term,
because the split precedes the collocation). This isolated-blob realization
leaves the Γ equation untouched; global conservation alone does not forbid
conservative interparticle Γ transfer. What the derivation does **not** fix is
the scale inside $\nu_t$ — the eddy-diffusivity coefficient is a model input.
§3's R1/R2 then act as the selection principle among scales, and the
fixed-filter choice (3b) is the one consistent with reading
$\zeta_{\sigma_0}$ as the LES filter whose subfilter content is being
modeled. (Alvarez's stated reason for rejecting Smagorinsky-type models —
"overly diffusive in simulations with coherent vortical structures"
[A §3.1, p. 34] — concerns the **Γ-side vorticity-stress** model for
$\mathbf E_{\mathrm{str}}$; it does not bear on a σ-channel diffusivity,
whose action is confined to core-scale transport and which vanishes with
$|S|$ at resolved coherent structures. Note also Cottet's remark, quoted at
A p. 16/B p. 4, attributing $\nabla\cdot\overline{\boldsymbol\omega}$ growth
to "unphysical small scales… which could be properly addressed with
subfilter-scale diffusion" — currently palliated by Pedrizzetti relaxation:
independent literature support for a diffusive subfilter σ-channel.)

---

## 3. Candidate closure survey

Each candidate: provenance, form, the constants' meaning, and its action on
the discrete map (full analysis in §4). Two structural results organize the
survey; both are proved in §4:

- **(R1) Scale-invariance no-go (DERIVED):** any closure of the form
  $\nu_t = (C\,\sigma_p)^2\,F(\text{resolved strain})$ — i.e. using the
  particle's **own dynamic σ** as the mixing length — is scale-free: because
  both the strain-thinning term and $\nu_t/\sigma$ then scale identically
  with σ, the closure only rescales the laminar equilibrium
  $\sigma_{eq} \to \sigma_{eq}/\sqrt{1 - C^2 F/Z}$ (for $C^2F<Z$) or tips the
  core into unconditional growth (for $C^2F>Z$). It can never introduce a
  *cutoff length*: the $\Delta t Z$ ceiling moves at most by a constant
  factor, with a sensitive pole at $C^2F \to Z$. Standard LES escapes this
  because its filter width is *fixed*; in rVPM the filter width is the
  dynamical variable under attack.
- **(R2) A closure must carry its own scale (DERIVED, corollary of R1):**
  admissible σ-channel closures therefore import a scale from outside the
  particle's own state: a fixed resolution length (the shed core size
  $\sigma_0$), a circulation-based velocity scale ($a_1|\Gamma_p|$, which is
  σ-independent), or the molecular ν (the laminar limit).

### 3a. Smagorinsky-σ with the dynamic core as mixing length — REJECTED by R1

$$\nu_t = \nu + (C_s\,\sigma_p)^2\,|S|,\qquad \dot\sigma_{\mathrm{visc}} = \nu_t/\sigma_p,$$

with $|S|$ the resolved strain-rate magnitude at the particle (available from
the existing velocity-gradient evaluation). Provenance: direct transplant of
the Smagorinsky eddy viscosity with $\ell = C_s\sigma$. Verdict (DERIVED, §4.2):
falls to R1 — it rescales $\sigma_{eq}$ without bounding $\Delta t Z$, and its
behavior is discontinuously sensitive at $C_s^2|S| \approx Z$. Retained in the
survey as the instructive failure that motivates R2.

### 3b. Fixed-filter Smagorinsky: the shed scale as mixing length — PRIMARY CANDIDATE

$$\boxed{\;\nu_t(\mathbf x_p) = \nu + (\kappa\,\sigma_0)^2\,|S(\mathbf x_p)|\;}
\qquad \dot\sigma_{\mathrm{visc}} = \nu_t/\sigma_p,$$

where $\sigma_0$ is the **shed** (target-resolution) core size — a constant of
the discretization, not the particle's dynamic σ — and $\kappa = O(1/2)$ a
dimensionless constant. Provenance: this is the LES stance made local — it is
exactly 019's remedy 2, $\nu_{\mathrm{eff}} = (\sigma_{\mathrm{target}}/2)^2\bar Z$,
with the global mean strain replaced by the local resolved strain, recovering
uniform-ν_eff as its spatial average. Constants' meaning: $\kappa\sigma_0$ is
the subgrid cutoff — the smallest core the model maintains against strain.
Key property (DERIVED, §4.3): the strained-core equilibrium becomes
$\sigma^* = \kappa\sigma_0\sqrt{|S|/Z} \ge \sqrt5\,\kappa\sigma_0$ (using
$|S| \ge s = 5Z$; the ratio $|S|/Z$ is an alignment factor, not a magnitude)
— **independent of the strain magnitude**: every core equilibrates near the
resolution scale, strong strain or weak, which is precisely "sub-σ
strain-thinning balances sub-σ turbulent mixing at the filter scale."
Laminar-limit recovery: $|S|\to 0 \Rightarrow \nu_t \to \nu$ ✓; quiescent
no-op up to molecular CoreSpreading ✓.

### 3c. Re_Γ-scaled effective viscosity (Squire / Bhagwat–Leishman) — physical anchor, insufficient alone

$$\nu_t = \nu\,\big(1 + a_1\,Re_v\big) = \nu + a_1|\Gamma_v|,\qquad a_1 \approx 6.5\times10^{-5}.$$

Provenance: Squire (1965) hypothesis; Bhagwat & Leishman's generalized
viscous vortex model, empirically validated on measured tip-vortex core
growth; per-particle version $\nu_t = \nu + a_1\Vert\boldsymbol\Gamma_p\Vert/(\text{—})$
requires a circulation proxy (the particle strength is not a tube
circulation; see §5 consistency note). σ-independent per particle, so it
passes R2 and creates an equilibrium $\sigma^*_t = \sqrt{\nu_t/Z}$ that is
$\sqrt{\delta}\approx1.5\times$ the laminar one at this scale. Verdict
(DERIVED, §4.4): **too weak to be the stabilizing closure at this flow's
constants** — arresting filament self-thinning needs
$\nu_t \gtrsim \Gamma_v/2\pi \approx 0.044\ \mathrm{m^2/s} \approx 3000\nu$,
while BL provides $\approx 2.3\nu$. Its role: the **validation anchor** —
whatever closure is adopted must reproduce BL core growth rates on an isolated
aged tip vortex (a Phase-4 check), and the honest physical statement that the
*true* subgrid viscosity at this scale is small: the stabilizing cutoff is a
*resolution* statement (3b), not a claim about physical turbulence intensity.

### 3d. Winckelmans CSM / CSM+PSE hybrids and the in-code precedent

Provenance: core-spreading method with RBF spatial adaptation (Barba);
Winckelmans & Leonard 1993; `zuhaldung_2014` CSM. FLOWVPM already ships the
machinery: CoreSpreading carries a growth tracker and, past
$\beta_{\mathrm{cur}} \ge \beta$, **resets** σ to `sgm0` and re-solves particle
strengths by RBF conjugate-gradient to preserve the vorticity field
(`FLOWVPM_viscous.jl:191-212, 314`). Relevance: (i) precedent that σ can be
re-set at fixed represented vorticity — the conservative reset operator a
closure-driven σ *decrease* would also use; (ii) the natural home for the
closure implementation is a `ViscousScheme` subtype (per-particle $\nu_t$ in
`viscousdiffusion`), since ν currently enters only as one scalar per scheme
and **no per-particle or strain-dependent viscosity hook exists anywhere**
(CONFIRMED, grep). Not itself a closure — a transport mechanism the closure
parameterizes.

### 3e. Γ-side consistency: must the closure also modify the Γ equation?

The $(f,g)$ machinery ties $\dot\sigma$ and $\dot{\boldsymbol\Gamma}$ through
the conservation constraints of the moment matching (§2). The question is
whether adding $\dot\sigma = \nu_t/\sigma$ demands a compensating Γ term.
Answer (DERIVED, §2.3): **no compensating Γ term is required**. Under the
explicitly adopted isolated-blob, localized/no-transfer approximation, the
retained subfilter transport is routed into σ at fixed Γ. Global
total-vorticity conservation alone does not forbid redistribution among
particles. Consistently, viscous
core spreading changes σ at *fixed*
$\boldsymbol\Gamma_p$, exactly as Lamb–Oseen diffusion conserves circulation
while spreading the core, and exactly as the live CoreSpreading already
composes with the inviscid $(f,g)$ update by viscous splitting
(`timeintegration.jl:161` → `viscous.jl:161`). Replacing
$\nu \to \nu_t(\mathbf x_p)$ inherits that structure term-by-term. What the
closure must **not** do is touch the $-3Z\boldsymbol\Gamma$ term: that term is
the inviscid moment-matching consequence of $\dot\sigma_{\mathrm{stretch}}$,
and §4.5 shows the ignition threshold it sets ($\Delta t Z = 2/3$, transverse)
is removed by integration (3f), not by modeling. The total-vorticity and
impulse budgets under this position are checked in §5.

### 3f. The stiff local integrator — mandatory non-model companion

Replace the forward-Euler local update at frozen velocity gradient
$\nabla\overline{\mathbf u}$ with the exact solution of the *linear* local
system over the step. Since
$\mathbf S = (\nabla\overline{\mathbf u})^{\mathsf T}\boldsymbol\Gamma$ (or the
transposed variant), the local ODE is

$$
\dot{\boldsymbol\Gamma} = A(t)\,\boldsymbol\Gamma
\ \text{ with }\
A = (\nabla\overline{\mathbf u})^{\mathsf T} - 3Z(\boldsymbol\Gamma)\,I,
\qquad
\dot\sigma = -Z(\boldsymbol\Gamma)\,\sigma,
$$

nonlinear only through the scalar $Z(\boldsymbol\Gamma)$. Exact-in-Z
integration (matrix exponential of the frozen gradient plus the scalar norm
correction, or per-particle sub-stepping as the cheap equivalent) gives
multipliers $e^{\lambda\Delta t}$: **no sign flips at any $\Delta t Z$**, σ>0
preserved exactly, axial growth at the physical rate $e^{2\Delta t Z}$ instead
of the Euler undershoot, and the spurious $\Delta t Z\in\{2/3, 1, 2\}$
thresholds gone. What it fixes: ignition as a *discrete overshoot* — the
transverse flip route and negative σ. What it cannot fix (DERIVED, §4.5): the
field-coupled contraction runaway — with lagged Z and $Z \propto 1/\sigma^2$,
$\Delta t Z_{n+1} \approx \Delta t Z_n\,e^{+2\Delta t Z_n}$ grows
super-exponentially once contraction outruns diffusion, and the exactly
integrated core still collapses to the (laminar) viscous floor, carrying
$Z_{\max}$ far beyond $2/\Delta t$ *in the field* even though each local step
is stable. Molecular CoreSpreading is composed with contraction through the
constant-effective-$Z$ exact update for $y=\sigma^2$; this replaces the old
full-step viscous Lie split and remains positivity-preserving for $\nu\ge0$.
Consistency in the smooth limit: $e^{x} = 1+x+O(x^2)$ reproduces the Euler
update to first order — bit-drift only at $O(\Delta t^2 Z^2)$, i.e. nowhere in
a healthy wake (ambient $\Delta t Z \sim 10^{-3}$).

**Design conclusion of the survey (conditional scalar-model result in §4.5): closure and
integrator are complements, not alternatives.** The fixed-filter closure (3b)
bounds the *equilibrium* strain the field can present
($Z \lesssim \Gamma_v/2\pi(\kappa\sigma_0)^2$ in the scalar filament model)
but a transient can still
overshoot a threshold before σ re-equilibrates; the exact integrator removes
every finite-$\Delta t Z$ threshold but leaves the laminar runaway. Together:
no scalar-model threshold to cross and no scalar-model runaway to feed. This
motivates, but does not prove generally, the item's mandate that 3f accompany
any closure.

---

## 4. The closed-loop discrete map with closure

*(Analysis section — the algebra below is checked numerically in the
scratchpad companion; see the Log entry.)*

### 4.1 Setup

Augment step 4 of the map (§1.3) with a per-particle $\nu_t$:

$$
y_{n+1} = (1-\Delta t Z_n)^2\,y_n + 2\,\nu_t(y_n, |S|_n)\,\Delta t,
\qquad y \equiv \sigma^2 .
$$

Field coupling (the "closed loop"): Z is lagged, and the worst-case strainer
is the compact filament measured in 018/019 (margin-curve exponent $p\approx2$
expected; $\Gamma_{\mathrm{implied}} \approx \Gamma_v$):

$$
Z_n = \frac{\Gamma_v}{2\pi\,\max(y_n^{\mathrm{strainer}},\, r_c^2)} + \bar Z,
$$

where the strainer's own core obeys the same map — self-consistent
contraction — and $r_c$ is the physical core radius (unreachable without DNS;
kept only to mark where physics saturates).

**Kelvin bookkeeping (DERIVED; a subtlety that changes the loop's verdict).**
The circulation in the numerator is the **tube** circulation
$\Gamma_v$ — a Kelvin invariant of the strain source — *not* the particle
strength $\Vert\boldsymbol\Gamma_p\Vert$. Under stretching,
$\Vert\boldsymbol\Gamma_p\Vert$ grows as *segment lengthening*
($\Gamma_p \approx \Gamma_{\mathrm{tube}}\,\Delta\ell_p$) while the strain a
neighbor feels scales with $\Gamma_{\mathrm{tube}}$, which is constant. A
closed-loop model that feeds the growing $\Vert\boldsymbol\Gamma_p\Vert$
back into Z manufactures a spurious mutual-intensification runaway that no
closure can bound (checked numerically: that mis-specified loop diverges at
any sub-stepping). This is also exactly why 019's
$\Gamma_{\mathrm{implied}} = 2\pi\sigma^2 M/\Delta t \approx \Gamma_v \approx$
const across the ladder is the *validation* of the filament picture — the
strain source presents a conserved circulation, not a growing strength.

### 4.2 R1, proved (dynamic-σ closures are scale-free)

With $\nu_t = \nu + C^2 y\,F$ ($F$ any function of resolved strain,
$Z = \alpha F$ with alignment factor $\alpha \le g$ fixed by geometry):

$$
y_{n+1} = y_n\big[(1-\Delta t Z)^2 + 2C^2 F\,\Delta t\big] + 2\nu\Delta t .
$$

The closure moved from the offset into the **multiplier**: the map is still
linear in y with an inhomogeneous term $2\nu\Delta t$ set by molecular ν
alone. Continuous limit $\dot y = 2y(C^2F - Z) + 2\nu$:

- $C^2 F < Z$: equilibrium $y^* = \nu/(Z - C^2F)$ — the **laminar**
  $\sigma_{eq}$ inflated by $(1 - C^2F/Z)^{-1/2}$. Since $F/Z = 1/\alpha$ is a
  pure alignment ratio, independent of σ and of strain magnitude, the
  inflation is a *constant factor*: no new length scale, $Z_{\max}$ at the
  equilibrium still $\propto 1/y^*\propto Z$ — the ceiling structure is
  untouched.
- $C^2 F > Z$: y grows exponentially wherever strain exists — unconditional
  fattening, saturated only by field-level smoothing; and the boundary
  between the two behaviors is an alignment condition, not a tunable margin.

∎ (R1). A dynamic-σ mixing length cannot set a cutoff. Note this is special
to rVPM: in grid LES, $\Delta$ is not a dynamical variable.

### 4.3 The fixed-filter closure (3b): equilibrium at the resolution scale

With $\nu_t = \nu + (\kappa\sigma_0)^2|S|$ and $|S| \ge s = 5Z$ (equality when
strain is fully axial):

$$
\dot y = -2Zy + 2\nu + 2(\kappa\sigma_0)^2|S|
\quad\Longrightarrow\quad
y^* = \frac{\nu}{Z} + (\kappa\sigma_0)^2\,\frac{|S|}{Z}
\;\ge\; \frac{\nu}{Z} + 5\,(\kappa\sigma_0)^2 .
$$

**The strained-core equilibrium is bounded below by
$\sigma^* \ge \sqrt5\,\kappa\sigma_0$ regardless of strain magnitude**
(DERIVED): stronger strain thins faster but also mixes faster, in exactly
compensating proportion, so cores park at the resolution scale — the wanted
subgrid-cutoff behavior. Stability of the fixed point: multiplier
$(1-\Delta t Z)^2 < 1$ for $\Delta t Z < 2$; convergence e-fold time
$1/(2Z)$ — fast precisely where strain is strong.

Closed loop: the strainer's core is also bounded, $y_s \ge 5(\kappa\sigma_0)^2$,
so the resolvable strain is capped:

$$
Z_{\max} = \frac{\Gamma_v}{2\pi\,y_s^{\min}} = \frac{\Gamma_v}{10\pi(\kappa\sigma_0)^2},
\qquad
(\Delta t Z)_{\max} = \frac{1}{5}\Big(\frac{\sigma_{\mathrm{stab}}}{\kappa\sigma_0}\Big)^2,
\quad \sigma_{\mathrm{stab}} \equiv \sqrt{\frac{\Gamma_v\Delta t}{2\pi}} .
$$

**Consequence (DERIVED, honest):** with the Euler update retained, freedom
from the transverse threshold $(\Delta t Z)_{\max} < 2/3$ requires
$\kappa\sigma_0 > \sqrt{3/10}\,\sigma_{\mathrm{stab}} \approx 0.55\,\sigma_{\mathrm{stab}}$
— for $\kappa = 1/2$, $\sigma_0 \gtrsim 1.1\,\sigma_{\mathrm{stab}}$: the
closure alone buys **little** operating range below today's boundary (and the
bound above is an *equilibrium* statement — transients can overshoot it).
This is the quantitative form of "a σ floor alone is ruled out" (018 §6) and
the reason 3f is mandatory: the closure's value is fidelity (cores equilibrate
at a chosen, meaningful scale instead of the laminar floor), while the
*stability* freedom comes from the integrator (§4.5).

### 4.4 The BL/Squire closure (3c): the physical-constants no-arrest result

Filament self-thinning with any σ-independent $\nu_t$: put
$Z = \Gamma_v/(2\pi y)$ (partner filaments at matched σ) into
$\dot y = -2Zy + 2\nu_t$:

$$
\dot y = -\frac{\Gamma_v}{\pi}\Big(1 - \frac{2\pi \nu_t}{\Gamma_v}\Big)
\quad\text{(}y\text{-independent!)}
\qquad\Longrightarrow\qquad
y(t) = y_0 - \frac{\Gamma_v}{\pi}\Big(1 - \frac{2\pi\nu_t}{\Gamma_v}\Big)\,t .
$$

The contraction is *linear in* $y$-space and reaches any floor in finite time
unless $\nu_t \ge \Gamma_v/2\pi$. At this rotor,
$\Gamma_v/2\pi = 0.044\ \mathrm{m^2/s} \approx 3000\,\nu$, while BL supplies
$\nu_t = \delta\nu \approx 2.3\,\nu$: the physically-anchored eddy viscosity
**slows the collapse by ~0.1% and cannot arrest it** (DERIVED). This is the
precise sense in which "reaching $r_c$ faithfully means DNS of the core"
(discussion item 4): the physical core is maintained by *resolved* core-scale
turbulence, not by an effective viscosity a one-particle core can host.
Verdict: 3c is an anchor and a fidelity target, not the stabilizing closure —
and any claim that a "physical" closure stabilizes this flow at its true
constants should be treated as refuted by this bound.

### 4.5 The combined result: closure + exact integrator

**Historical conditional claim (WITHDRAWN as a general proof by the 2026-08-12
audit).** Under (i) the fixed-filter
closure $\nu_t = \nu + (\kappa\sigma_0)^2|S|$ applied through the
positivity-preserving viscous update
$\sigma \leftarrow \sqrt{\sigma^2 + 2\nu_t\Delta t}$, and (ii) the exact-in-Z
local inviscid update (3f), the per-particle map satisfies, for **all**
$\Delta t > 0$, $\Gamma$, σ > 0:

1. $\sigma_{n+1} > 0$ (no sign flips: the inviscid factor is $e^{-\Delta t Z}>0$,
   the viscous update adds a positive quantity under a square root);
2. $y_n \ge \min\!\big(y_0,\; 2\nu_t^{\min}\Delta t\big)$ and, in sustained
   strain, $y_n \to y^* \ge 5(\kappa\sigma_0)^2 + \nu/Z$: cores equilibrate at
   the resolution scale, not the laminar floor;
3. **on the attractor** the field-coupled strain is bounded,
   $Z \le \Gamma_v/(10\pi(\kappa\sigma_0)^2) + \bar Z$ (Kelvin-constant
   source, §4.1), hence $\Delta t\,Z$ is bounded for fixed
   $(\Delta t, \sigma_0)$; transients present whatever strain the initial
   state implies (a deep-collapsed core reads a large $\Delta t Z$ until it
   re-equilibrates), but because no discrete threshold exists under (ii),
   **no finite $\Delta t Z$ is catastrophic**: regime 2 is unreachable in
   the map-level sense — no update's multiplier magnitude departs from its
   exact value $e^{\lambda\Delta t}$, every state converges to the bounded
   attractor, and nothing irreversible happens on the way;
4. neither ingredient suffices alone: dropping (ii) leaves the transverse
   flips wherever the §4.3 cap exceeds $2/3$ (numerically confirmed at
   $\sigma_0 = 0.02R$, cap $\approx 2.4$: the Euler+closure run survives
   *with* sign flips); dropping (i) leaves the super-exponential lagged-Z
   runaway $\Delta t Z_{n+1} \approx \Delta t Z_n e^{2\Delta t Z_n}$ ending
   at the laminar floor with unbounded strain (numerically: $\Delta t Z$
   reaches $10^3$ and beyond).

*Proof sketch and numerical exercise* — 1 and 2 are immediate from the
update forms; 3 follows by induction: if every strainer core satisfies 2 at
step n, the induced Z at step n+1 obeys the cap, which sustains 2. The
scratchpad companion iterates the Kelvin-correct closed loop
(§4.1) through all four (integrator × closure) quadrants at the campaign
constants ($\sigma_0 = 0.02R$, below today's ignition boundary): only
exact+closure is flip-free and bounded, with peak $\Delta t Z$ = the §4.3
cap (2.43 vs predicted 2.41) and settled $\sigma/\sigma_0 = 1.05$
(sub-stepped composed update) — the §4.6 $\kappa$ design landing where
intended. Two implementation-relevant numerical findings, recorded for
Phase 3: (a) operator-splitting the exact update at full $\Delta t$
(inviscid exponential, then viscous add) parks the equilibrium ~45% above
the composed-ODE value at $\Delta t Z \sim 2$ — benign (extra margin) but
the Phase-3 T2 tolerance must be stated against the *scheme's own* fixed
point, not the continuous one; (b) the formal weak point of the induction —
simultaneity of many strainers raising $\bar Z$ — is not exercised by the
scalar loop and is flagged as the one place Phase-4 empirics outrank
Phase-1 algebra.

**Interpretation.** The composed design replaces "clip the update" with: the
integrator makes the *discretization* faithful at any local stiffness, and
the closure makes the *model* possess the equilibrium that physics (turbulent
core transport) provides and molecular ν does not — which is exactly the
item's objective statement.

### 4.6 Constants

Two constants, both with stated meaning and defaults to be agreed at the gate:

| constant | meaning | default proposal | falsifiable via |
| --- | --- | --- | --- |
| $\kappa$ | subgrid cutoff in units of shed σ: cores equilibrate at $\sqrt5\kappa\sigma_0$ | $\kappa = 1/\sqrt5 \approx 0.45$ (equilibrium exactly at $\sigma_0$: shed scale = maintained scale) | σ/σ_shed settled histograms (fid144-class runs) |
| $a_1$ (anchor, not a knob) | Squire parameter in the BL validation target | $6.5\times10^{-5}$ (literature) | isolated aged-vortex core-growth check (Phase 4) |

No dynamic procedure is proposed for $\kappa$ in Phase 1: the Γ-side SFS keeps
its own $C_d$ machinery untouched, and $\kappa$ is a resolution statement
(like $\sigma_0$ itself), not a flow constant. Revisit only if Phase-4
cross-family validation fails at fixed $\kappa$ (the item's no-retuning rule).

---

## 5. Conservation and consistency constraints

The closure must satisfy, and Phase 3 must test (tolerances in §8):

1. **Total vorticity / strength invariance (viscous channel):**
   $\dot\sigma = \nu_t/\sigma$ at fixed $\boldsymbol\Gamma_p$ conserves
   $\sum_p\boldsymbol\Gamma_p$ trivially and, because
   $\int\zeta_\sigma\,\mathrm d\mathbf x = 1$ for all σ, conserves total
   represented vorticity $\int\overline{\boldsymbol\omega}\,\mathrm d\mathbf x$
   exactly (DERIVED). Same argument as molecular CoreSpreading.
2. **Impulse correction:** linear impulse
   $\tfrac12\sum_p\mathbf x_p\times\boldsymbol\Gamma_p$ is untouched by any
   pure-σ update (positions and strengths fixed) and is conserved exactly.
   Angular impulse is **not** invariant because an isotropic blob contributes
   a term proportional to $\sigma^2\boldsymbol\Gamma$; its predicted change
   must be tested instead.
3. **Energy/enstrophy budget sign:** the closure must be dissipative on
   enstrophy at sub-σ scales. Spreading a fixed-strength core strictly lowers
   $\int|\overline{\boldsymbol\omega}|^2$ (Gaussian: $\propto \sigma^{-3}$ at
   fixed Γ) — so $\nu_t \ge 0$ suffices; **no backscatter channel is
   introduced** (backscatter remains the Γ-side SFS's clipped business).
4. **Burgers/laminar recovery:** $|S| \to 0 \Rightarrow \nu_t \to \nu$ and
   the map reduces to the live CoreSpreading map with its
   $\sigma_{eq} = \sqrt{\nu/Z}$ fixed point (the 018 §2 physics) — exact by
   construction of 3b.
5. **Zero-strain no-op:** quiescent wake ($|S| = 0$): closure adds nothing
   beyond molecular CoreSpreading; inviscid-mode ($\nu = 0$) with closure off:
   bit-identical to the live code (the Phase-3 off-state contract).
6. **Γ-equation non-interference (from 3e):** the closure writes only the
   viscous σ channel; the inviscid $(f,g)$ σ/Γ coupling is untouched, so the
   collocation-and-ansatz structure that fixes $f=0, g=1/5$ is preserved
   term-by-term. Under the isolated-blob, localized/no-transfer approximation,
   the viscous/turbulent split precedes the collocation and needs no
   compensating Γ term; this is not a uniqueness claim from global
   conservation alone.
7. **BL anchor (fidelity, Phase 4):** an isolated aged vortex under the
   closure must grow its core within the BL envelope
   $r_c(t) = \sqrt{r_0^2 + 4\alpha_L\delta\nu t}$ scaled to this $Re_v$ —
   i.e. the closure must not *over*-diffuse a resolved, weakly-strained
   vortex (the known Smagorinsky failure mode in coherent-structure flows
   that Alvarez's Γ-side SFS was designed to avoid; the fixed-filter form
   inherits protection because $|S|$ at a resolved isolated vortex core is
   small by construction — its strain lives at radius $\sim r_c$, resolved).

---

## 6. Falsifiable predictions (written before any Phase-2 data is examined)

P-1 (**strain ceiling**): the 019 margin curve will show $M(\sigma)\propto
\sigma^{-p}$ with $p \approx 2$ and $\Gamma_{\mathrm{implied}} =
2\pi\sigma^2 M/\Delta t \approx \Gamma_v$ across the ladder — today's max
resolvable strain is stability-set, not physics-set. (Already the 019
pre-registered interpretation; 020 inherits it as its gap premise.)

P-2 (**integrator, Stage A**): re-integrating the recorded collapse
trajectories under the exact local map removes every sign flip and finite-step
divergence *pointwise* (no negative σ, no transverse overshoot), while still
tracking σ into deep contraction (the laminar floor) — i.e. Stage A shows
stability without fidelity.

P-3 (**integrator, Stage B**): a below-boundary run (σ ≈ 0.02R viscous class)
with the exact integrator alone (no closure) does **not** ignite in the
discrete-overshoot sense but still develops unbounded-growing
$|\Gamma|/\sigma^2$ tails / M readings far above ε as cores grind to the
floor — the boundary *moves* but a fidelity collapse remains. (This is the
sharpened item-3 claim of §1.3; if instead the run is healthy and bounded,
the field-runaway half of the claim is refuted and the closure's stability
motivation weakens to fidelity-only — recorded either way.)

P-4 (**SFS null**): SFS-on vs SFS-off at a collapsing σ shows σ-trajectory
and ignition unchanged (the $f/(1+3f)=0$ structural null, plus measured
Γ-side unimportance in collapse).

P-5 (**closure, Phase 4**): with 3b at $\kappa = 1/\sqrt5$ and the integrator,
previously-igniting 019 regime-map cells (σ/R ∈ {0.015…0.030}) run stable
with settled median σ/σ_shed ≈ 1 (not the laminar-floor migration), M bounded
per §4.3's cap; and the σ-iterate procedure of 019 converges at smaller σ*.
Loads: CT and Γ(r/R) move *toward* the 018 B0/L1 anchors or stay within their
error budget — a large load shift at matched resolution would falsify the
"stability without aerodynamic contamination" premise (cf. the measured NULL
ν-delta at L1).

P-6 (**constant-revs reconciliation**): the `scr_ufdt_nt{36,72,144}`
constant-physical-time ignitions (019 open question) are explained by the
field-runaway mechanism (contraction rate $\dot\sigma/\sigma = -Z$ is
Δt-independent, so the grind to the floor takes fixed physical time; the
Euler threshold merely dates the death, the contraction drives it). The
closure removes the grind itself ⇒ prediction: with the **closure+integrator
pair** ON, nt36/72/144 at 0.0291R all survive past 8 revs. (Both ingredients
are required by the prediction: at $\sigma_0 = 0.0291R$, $\kappa = 1/\sqrt5$,
the §4.3 cap is $\approx 1.1 > 2/3$ — closure-only would leave Euler
transverse flips in play.) A Δt-dependent survival pattern under the pair
would falsify the mechanism attribution.

---

## 7. Phase-2 protocol — PRE-REGISTERED (gap demonstration)

Per Ryan's scoping decision (2026-08-06): registered now, against the
in-flight 019 Campaign E runs; their outcomes test this registration.
Metrics and thresholds fixed here; any deviation at execution time is
recorded in the item Log, not patched.

### 7.1 Leg 1 — the strain ceiling, measured (mostly reuse; no new runs expected)

- **Data:** 019 pipeline outputs (margin curve + regime map) from Campaign E
  jobs 13060963-66 (`s015v/s020v/s025v/s030v`), 13061166-69
  (`s015/s025/s038/s038v`), 13061089 (`fid144`), 13064696 (`sstab`), plus the
  018 reuse points already in the 019 REGISTRY. Source of record: monitors
  CSVs via `scripts/p019_sigma_procedure.py` — never stdout (018 rule).
- **Metrics:** fitted $p$ from $M(\sigma)\propto\sigma^{-p}$ (screen-class,
  per viscosity row, as 019-registered); $\Gamma_{\mathrm{implied}}(\sigma) =
  2\pi\sigma^2M/\Delta t$; the measured ignition boundary
  $\sigma_{\mathrm{ign}}$.
- **Thresholds (data-blind):** the ceiling is "stability-set" if
  $p \ge 1.5$ over the ladder and
  $\Gamma_{\mathrm{implied}} \in [0.5, 2]\,\Gamma_v$ for ≥ 3 consecutive
  rungs. Physics target: $Z_{\mathrm{phys}} = \Gamma_v/(2\pi r_c^2)$ with
  $r_c$ from literature tip-vortex cores at matched $Re_v$
  ($r_c \approx 0.01$–$0.05\,c$; both endpoints reported). The gap statement
  is the decade count $\log_{10}(Z_{\mathrm{phys}}/Z_{\max}^{\mathrm{resolved}})$
  at the boundary, reported with the $r_c$ bracket.

### 7.2 Leg 2 — the stiff-integration test

- **Stage A (offline, no src change):** re-integrate
  `data/p018_L2_visc_forensics/death_trajectory.csv` and the
  `p018_ufront_n1{,_visc}` wake-health rows under (a) the live Euler map and
  (b) the exact local map of §3f (axial/transverse split, frozen per-step Z
  from the recorded data), composed with the CoreSpreading update where the
  source run was viscous. **Pass criterion for P-2:** under (b), zero
  negative-σ events and zero per-step multiplier magnitudes > its exact
  $e^{|\lambda|\Delta t}$; σ still reaches ≤ 1.1× the recorded floor
  (fidelity failure persists). Deliverable: overlay trajectories, one figure.
- **Stage B (in-code prototype; default-off/bit-identical/test-gated, the
  standard src contract — implementation happens in Phase 2 only as this
  single gated prototype):** one screen-class run at σ/R = 0.02, viscous,
  NT=36, N=1, D=3.4, uniform-d_front recipe, `max_dtZ` on, integrator ON,
  closure OFF. **Registered readouts:** survival to 8 revs; M trajectory;
  settled min/median σ/σ_shed; max $|\Gamma|/\sigma^2$ tail. P-3 verdict
  rules: "threshold moves" = survives 8 revs with no negative σ and no
  merge-OOM; "fidelity collapse remains" = settled median σ/σ_shed < 0.85 or
  M > ε=0.2 sustained. All four outcome quadrants have stated meaning (item
  file, Phase-2 objective); none is a protocol failure.

### 7.3 Leg 3 — the SFS null

- **Derivation leg:** §1.1's structural null restated standalone (done in
  this document; cite 018 corpse evidence).
- **Measured leg:** reuse first — the 018 record already carries SFS-on
  collapses (campaign runs had SFS on via `sfs_choice`,
  `rotor_hover_pressure_comparison.jl:420-509`) and the code-audit null. If
  no clean SFS-off twin exists at a collapsing σ in the 018/019 ledgers
  (check at execution time), run ONE screen pair at σ/R = 0.025 viscous
  (the mid-ladder igniter), SFS on/off, identical otherwise. **Threshold:**
  time-to-ignition within ±15% and min-σ trajectory overlay within line
  width ⇒ SFS viscous-null CONFIRMED at collapse scale.

### 7.4 Evidence pack

One self-contained document in the item subdir composing the three legs, with
publication-standard figures per the repo figure policy (standalone TikZ +
same-named CSV dirs), closing with the item's registered claim: *a subfilter
viscous (σ-channel) model is the missing piece that would let the simulation
approach physics-predicted strain rates.* Gate: Ryan review; explicit approval
before any Phase-3 work.

---

## 8. Phase-3 consistency tests — PRE-REGISTERED (tolerances fixed now)

All tests live in FLOWVPM's test suite (upstream repo), mirror the 019
`max_dtZ` test pattern (exact map replication), and gate the Phase-3 PR.
Implementation contract: default-off, bit-identical when off, env-knob
exposed in the FLOWPanel driver, `_case_metadata.toml` + banner carry it.

| # | test | tolerance |
| --- | --- | --- |
| T1 | single-particle exact map vs closed-form: prescribed constant $\nabla\overline{\mathbf u}$, N steps of integrator+closure vs the analytic solution of §3f/§4.3 (axial, transverse, σ) | rtol 1e-10 |
| T2 | Burgers equilibrium: sustained axial strain Z, closure ON ⇒ $y \to y^*$ = the **implemented scheme's own** discrete fixed point (splitting parks it above the continuous $\nu/Z + (\kappa\sigma_0)^2|S|/Z$ by the $(1-\Delta t Z/2)^{-1}$-type factor — state the scheme's formula in the test) | rtol 1e-6 after 20 e-folds |
| T3 | laminar recovery: $|S|\to0$ run reproduces live CoreSpreading trajectory | rtol 1e-12 |
| T4 | moments under a pure closure update: $\sum\boldsymbol\Gamma_p$ and linear impulse invariant; angular impulse changes by the analytic $\sigma^2\boldsymbol\Gamma$ contribution | rtol 1e-12 |
| T5 | isolated Gaussian blob: closure contribution lowers analytic self-enstrophy; no whole-field monotonicity claim for arbitrary overlap/spatially varying diffusivity without a proof | sign check, 0 violations |
| T6 | off-state bit-identity: closure+integrator flags off ⇒ byte-identical trajectories vs current master on the vortex-ring example | exact (== on state vectors) |
| T7 | ignition fixture: the 019 σ-ladder ignition fixture (recorded collapse ICs) under integrator+closure: no negative σ, bounded M per §4.3's cap | 0 sign flips; $M \le 1.2\times$ predicted cap |
| T8 | stiffness stress: single particle, $\Delta t Z \in \{0.5, 1, 2, 5, 50\}$ frozen: integrator multipliers match $e^{\lambda\Delta t}$ | rtol 1e-10 |

Deviations discovered while implementing are recorded in the item Log with
the failing tolerance — never retuned silently (the 018/019 house rule).

---

## 9. References

- E. J. Alvarez (2022), *Reformulated Vortex Particle Method and Meshless
  Large Eddy Simulation of Multirotor Aircraft*, PhD dissertation, BYU.
  [local: `~/Dropbox/research/literature/alvarez_2022_reformulated…aa.pdf`]
- E. J. Alvarez & A. Ning (2024), *Stable Vortex Particle Method Formulation
  for Meshless Large-Eddy Simulation*, AIAA J., DOI 10.2514/1.J063045.
  [local: `alvarezning_2023_stable…aa.pdf`]
- H. B. Squire (1965), *The growth of a vortex in turbulent flow*, Aero. Q. 16.
- M. J. Bhagwat & J. G. Leishman (2002), *Generalized viscous vortex model
  for application to free-vortex wake and aeroacoustic calculations*, AHS 58.
- S. Ananthan & J. G. Leishman (2004) — $\alpha_L = 1.25643$ convention.
- G. S. Winckelmans & A. Leonard (1993), JCP 109. [local]
- P. Ploumhans & G. S. Winckelmans (1999/2000). [local]
- L. A. Barba et al. (2005), IJNMF 47 (RBF spatial adaptation).
- Repo: `018_.../sigma_blowup_mechanism.md`;
  `019_sigma_selection_procedure.md` (P1 + addendum);
  `FLOWVPM_timeintegration.jl`, `FLOWVPM_viscous.jl`,
  `FLOWVPM_subfilterscale{,_models}.jl`, `FLOWVPM_formulation.jl`,
  `docs/src/rVPM.md`.

Web sources for the Squire/BL constants:
[Bhagwat & Leishman, Generalized viscous vortex model (ResearchGate)](https://www.researchgate.net/publication/255470975_Generalized_viscous_vortex_model_for_application_to_free-vortex_wake_and_aeroacoustic_calculations) ·
[A Reynolds Number-Based Blade Tip Vortex Model](https://www.researchgate.net/publication/233593771_A_Reynolds_Number-Based_Blade_Tip_Vortex_Model) ·
[WES 2019 actuator-line tip correction (α_L, δ usage)](https://wes.copernicus.org/articles/4/369/2019/)

## Log

- 2026-08-06 — draft written; §2 merged from the dissertation/paper
  extraction (collocation structure, the three doors, the $E_{\mathrm{adv}}$
  re-derivation). Numerical check of every §4 claim: ALL PASS (session
  scratchpad `check_020_map.py`, 19 checks — fixed points, R1 inflation
  factor, the §4.3 cap and κ landing, the 4.4 linear-contraction rate to
  1e-9, the four-quadrant Kelvin-correct loop, the lagged-Z recursion, the
  §2.3 Laplacian identity). The check run surfaced and corrected two draft
  defects, recorded here per house rule: (1) a spurious ½ factor in the
  §4.4 contraction rate; (2) the first closed-loop model fed particle
  strength into Z — diverges unphysically; replaced by the Kelvin-constant
  loop (§4.1 remark). Known open points for the gate: κ default; the §4.5
  simultaneous-strainer caveat; the choice of $|S|$ estimator (full
  strain-rate magnitude vs the already-computed axial $s$) for the 3b
  closure.
