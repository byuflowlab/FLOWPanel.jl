# 025 Phase 0 — Peak Comparison and Default Decision

**Date:** 2026-08-20. **Script:** `phase_00_peaks.jl` (this directory; runs in
seconds, single thread). **Decision rule (Ryan, item ruling 1):** the family
with the lower maximum velocity and/or velocity gradient wins the default; if
Gaussian is lower than the existing compact-support family, Gaussian becomes
the new codebase default.

## Setup

All three candidates regularize the transverse profile of the straight
vortex-filament kernel. For an infinite filament of circulation $\Gamma$ at
perpendicular distance $h$, in units of $\Gamma/(2\pi)$:

$$
u_{\rm sing}(h) = \frac{1}{h}, \qquad
u_{\rm Vat}(h) = \frac{h}{\sqrt{h^4 + r_c^4}}, \qquad
u_{\rm cpt}(h) = \frac{h}{h^2 + \delta(h)}, \qquad
u_{\rm gau}(h) = \frac{1 - e^{-h^2/(2\sigma^2)}}{h},
$$

with $\delta(h) = (h - r_c)^2$ for $h < r_c$ and $0$ otherwise — the same
compact-support `regularize` already used by the source/doublet kernels
(`src/FLOWPanel_elements_fmm.jl:954`), transplanted to the filament; the
Gaussian is the Lamb–Oseen swirl profile.

**Fair-comparison convention (item ruling 3):** each family's core parameter
is set so its far-field matching distance — the $h$ beyond which the relative
error against $1/h$ is $\le$ tol — equals a common $\Delta r$ (normalized to
1):

$$
d_{\rm Vat} = r_c\left(\tfrac{2}{\rm tol}\right)^{1/4}, \qquad
d_{\rm cpt} = r_c \ \text{(exact beyond } r_c\text{)}, \qquad
d_{\rm gau} = \sigma\sqrt{2\ln(1/{\rm tol})}.
$$

Peaks are then meaningfully comparable: every candidate buys the same FMM
radius inflation. (Correction 2026-08-20, review finding 5: the Vatistas rule
$r_c(2/{\rm tol})^{1/4}$ is CONSERVATIVE — it yields velocity error
$\approx$ tol/4 at $\Delta r$; the asymptotically tight distance is
$r_c(1/(2\,{\rm tol}))^{1/4}$. Back-solving $r_c$ from the conservative
radius therefore forces a SMALLER core than the convention requires, i.e. it
PENALIZES Vatistas in this matched-$\Delta r$ framing — the original text
stated the bias direction backwards. A tight-radius back-solve would raise
the Vatistas $r_c$ by $(4)^{1/4} \approx 1.41\times$ and lower its peaks by
the same factor; the ranking below is unchanged.)

## Results (script output, 2×10⁶-point grid)

(Metric labels corrected 2026-08-20, review finding 4: the derivative metric
is the RADIAL derivative $|du/dh|$ of the transverse profile, not a full
gradient norm. The infinite-filament transverse gradient has two nonzero
entry types — $du/dh$ (radial) and $u/h$ (curvature) — so the operator norm
is the pointwise max of the two; the script now computes both, and for all
three families and all tolerances the operator-norm maximum COINCIDES with
$\max|du/dh|$ (the radial term dominates everywhere), so the rankings are
unchanged under the correct label.)

| tol | family | core param/Δr | u_max·Δr | at h/Δr | max\|du/dh\|·Δr² | max‖∇u‖_op·Δr² |
|---|---|---|---:|---:|---:|---:|
| 1e-4 | Vatistas n=2 | 0.0841 | 8.41 | 0.084 | 141.4 | 141.4 |
| 1e-4 | compact | 1 | **1.207** | 0.707 | **2.55** | **2.55** |
| 1e-4 | Gaussian | 0.233 | 1.94 | 0.369 | 9.21 | 9.21 |
| 1e-5 | Vatistas n=2 | 0.0473 | 14.95 | 0.047 | 447.2 | 447.2 |
| 1e-5 | compact | 1 | **1.207** | 0.707 | **2.55** | **2.55** |
| 1e-5 | Gaussian | 0.208 | 2.17 | 0.330 | 11.51 | 11.51 |
| 1e-6 | Vatistas n=2 | 0.0266 | 26.59 | 0.027 | 1414.2 | 1414.2 |
| 1e-6 | compact | 1 | **1.207** | 0.707 | **2.55** | **2.55** |
| 1e-6 | Gaussian | 0.190 | 2.37 | 0.302 | 13.82 | 13.82 |

Far-field check at $h=\Delta r$: Vatistas rel. err tol/4, compact exactly 0,
Gaussian exactly tol — all within convention.

## Findings

1. **Compact-support has the lowest peaks at every tolerance, on both
   metrics** — u_max·Δr = 1.207 (vs Gaussian 1.9–2.4, Vatistas 8.4–26.6) and
   max|du/dh|·Δr² = 2.55 (vs Gaussian 9.2–13.8, Vatistas 141–1414) — and its
   peaks are **tolerance-independent** (exactly singular beyond $r_c$, so the
   matched-Δr scaling never squeezes its core).
2. The Vatistas columns quantify why the current default conditions poorly
   under matched inflation: to buy the same far-field agreement its core must
   shrink like $\mathrm{tol}^{1/4}$, driving peaks up by the same factor.
3. The compact profile is $C^1$ at $r_c$ (interior slope $-1/r_c^2$ at
   $h\to r_c^-$ matches the singular $-1/h^2$; see Phase 1), so no gradient
   jump is introduced at the support boundary.

## DECISION (per ruling 1 + contingency chain)

Gaussian is **not** lower than compact-support ⇒ the default does NOT go to
Gaussian. Per the item contingency ("Phase 0 finds compact-support has the
lower peaks → compact-support becomes the default for ALL offset-regularized
panel kernels incl. VortexRing"):

- **New default: compact-support** for the VortexRing filament kernel
  (source/doublet already use it — the codebase becomes uniform).
- **Gaussian implemented as a selectable option** (marginal cost is small:
  the same single-switch mechanism carries all three families).
- **Vatistas n=2 retained as the selectable legacy option** for A/B and
  rollback.

FMM radius-inflation consequence (`radius_inflation`, tol = `FMM_RADIUS_TOL[]`
= 1e-6, kerneloffset 1e-3): Δr = 0.0376 m (Vatistas) → **0.001 m (compact,
tol-independent)** → the 023 body-pass forced-direct mechanism is removed at
its source.

## Addendum (2026-08-20, Ryan convention ruling): matched-CORE-SIZE table and revised decision

The tables above hold the far-field matching distance Δr fixed and back-solve
each family's core size rc — under that convention compact wins. Ryan's
review identified that the operationally pinned quantity is the core size
itself (`kerneloffset`), not Δr. At matched rc the peaks are
tolerance-independent for every family (exact, verified by converting the
tables above with rc/Δr = (tol/2)^{1/4}, 1, and 1/√(2 ln 1/tol)):

| family | u_max·rc [Γ/2π] | max|du/dh|·rc² [Γ/2π] | Δr(tol=1e-6) |
|---|---:|---:|---:|
| Vatistas n=2 | 0.71 | 1.00 | 37.6·rc |
| compact | 1.21 | 2.55 | 1.0·rc |
| **Gaussian** | **0.45** | **0.50** | 5.3·rc |

(Metric note 2026-08-20, review finding 4: the derivative column above is
the radial $|du/dh|$; the operator-norm maxima under matched rc coincide
with it — Gaussian 0.50, Vatistas 1.00, compact 2.55 in $\Gamma/2\pi\,r_c^2$
units — so the ranking is identical under either derivative metric.)

**REVISED DECISION (Ryan 2026-08-20): the matched-rc convention governs;
Gaussian has the lowest peak velocity AND gradient, so per ruling 1 the
default is GAUSSIAN.** Its Δr ≈ 5.3·rc at tol 1e-6 retains nearly all of the
FMM radius-inflation win (0.0376 m → 0.0053 m at production kerneloffset).
Compact remains selectable (it is the joint-tradeoff winner if rc were free
to grow to Δr); Vatistas remains the legacy option.
