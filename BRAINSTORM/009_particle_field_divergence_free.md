# Is the Particle-Induced Velocity Field Divergence-Free?

## Context

The rotor-hover wake (items 002–008) is carried by a FLOWVPM vortex-particle
field. Several downstream constructions — the pressure-Poisson monitor
(`PressureLaplace`), the unsteady Bernoulli term, and the overlap-coupled
strength evolution proposed in item 008 — implicitly assume the particle-induced
velocity field is incompressible, i.e. `\nabla\cdot u = 0`. If the regularized
(mollified) particle kernel quietly injected spurious dilatation, that would
show up as an error source in the pressure reconstruction and in any vorticity
balance that uses `(\omega\cdot\nabla)u = J\omega` with a non-trace-free `J`.

This item settles the question both analytically and numerically. The headline
claim to be verified is: **a vortex-particle field built from a radially
symmetric regularization is divergence-free everywhere and for arbitrary
overlap, exactly, not just in the far field.** The numerical check uses
ForwardDiff to differentiate *through the induced-velocity kernel itself* (not a
finite difference), so the measured divergence is limited only by floating-point
round-off, giving an unambiguous pass/fail.

This is a foundational sanity item: it does not by itself move CT, but it
removes (or exposes) an assumption that items 004–008 all lean on.

## Theory

### Setup and the kernel

FLOWVPM's default kernel is the Gaussian-error-function kernel
(`zeta_gauserf` / `g_gauserf`). The regularized Biot–Savart velocity induced by
a set of particles with positions `X_p`, vector strengths `\Gamma_p`, and core
sizes `\sigma_p` is

```math
u(x) = -\frac{1}{4\pi}\sum_p
        \frac{g\!\left(\rho_p\right)}{r_p^{3}}\,
        (x - X_p)\times\Gamma_p,
\qquad
r_p = |x - X_p|,\quad \rho_p = r_p/\sigma_p,
```

where `g` is the dimensionless regularization function. For the Gaussian-erf
kernel,

```math
g(\rho) = \operatorname{erf}\!\left(\tfrac{\rho}{\sqrt2}\right)
          - \sqrt{\tfrac{2}{\pi}}\,\rho\,e^{-\rho^2/2}.
```

The only structural facts the derivation uses are that `g` depends on position
through the scalar radius `r_p` alone, and that `g(\rho)/\rho^3` is finite and
smooth as `\rho\to0` (the regularization removes the singular `1/r^3` so the
field is `C^1` at every particle center).

### Each particle's contribution is divergence-free

Write the single-particle velocity as `u_p(x) = K_p(x)\times\Gamma_p` with the
purely radial vector field

```math
K_p(x) = -\frac{1}{4\pi}\,f(r_p)\,(x-X_p),
\qquad
f(r_p) = \frac{g(r_p/\sigma_p)}{r_p^{3}},
```

so `K_p = f(r_p)\,r` (up to the constant) with `r = x - X_p` and `f` a scalar
function of `r_p = |r|` only. Using the vector identity, for constant `\Gamma_p`,

```math
\nabla\cdot(K_p\times\Gamma_p)
 = (\nabla\times K_p)\cdot\Gamma_p - K_p\cdot(\nabla\times\Gamma_p)
 = (\nabla\times K_p)\cdot\Gamma_p,
```

since `\Gamma_p` is constant. The curl of any radial field of the form
`f(|r|)\,r` vanishes:

```math
\nabla\times\bigl(f(r_p)\,r\bigr)
 = \nabla f(r_p)\times r + f(r_p)\,(\nabla\times r)
 = f'(r_p)\,\frac{r}{r_p}\times r + 0
 = 0,
```

because `r\times r = 0` and `\nabla\times r = 0`. Therefore

```math
\boxed{\;\nabla\cdot u_p(x) = 0\;}
```

for **every** particle, at **every** field point, for **any** radial
regularization `g` (Gaussian-erf, Winckelmans algebraic, Gaussian, singular
`g\equiv1`, …) and any core size `\sigma_p`. By linearity the total field
`u=\sum_p u_p` is divergence-free as well. Overlap is irrelevant: the result
holds pointwise per particle, so it cannot be broken by summing overlapping
particles.

### Why this is true at a deeper level

The regularized velocity is the curl of a regularized vector potential. With the
mollified Green's function `G_\sigma = \zeta_\sigma * G` (where `G` is the
Laplace Green's function and `\zeta_\sigma` the radial mollifier),

```math
u(x) = \nabla\times\!\sum_p \Gamma_p\, G_\sigma(x - X_p),
\qquad
\nabla\cdot u = \nabla\cdot(\nabla\times A) \equiv 0.
```

So divergence-freeness is structural — it survives regularization precisely
because the mollifier acts on the streamfunction and the velocity is still an
exact curl. This also tells us what *would* break it: a kerneloffset / control
point treatment that evaluates the velocity at a *shifted* argument relative to
the curl operator, a per-component (anisotropic) core, a non-radial mollifier,
or a `\sigma` that varies in space at the evaluation point. Those are the
failure modes the numerical test should be able to detect if introduced.

### Relation to the regularized vorticity

Note `\nabla\cdot u = 0` does **not** require the regularized vorticity to equal
the curl of the regularized velocity particle-by-particle. The blob vorticity is
`\omega(x)=\sum_p\Gamma_p\zeta_\sigma(x-X_p)`, and `\nabla\times u` recovers a
*smoothed* version of it; the two differ by overlap-dependent terms (this is
exactly the subject of item 008). Divergence, by contrast, is clean: it is
identically zero regardless of overlap. The test must therefore check the
*divergence* (`\operatorname{tr} J`), not the curl, against zero.

## Numerical Verification (ForwardDiff)

`examples/particle_divergence_free_check.jl` implements the test. It defines the
Gaussian-erf induced-velocity kernel in pure Julia (so ForwardDiff differentiates
the analytic kernel cleanly, independent of FLOWVPM's mutating internals) and
forms the velocity Jacobian `J_{ab}=\partial u_a/\partial x_b` via
`ForwardDiff.jacobian`. The reported divergence is `\operatorname{tr} J`.

Test matrix — **2 particles**, three overlap regimes set by center separation
`d` relative to core size `\sigma` (equal cores `\sigma=1`):

- **(a) overlap a lot:** `d = 0.25\,\sigma`,
- **(b) overlap a little:** `d = 2\,\sigma`,
- **(c) far apart:** `d = 20\,\sigma`.

Strengths are non-parallel 3-vectors (e.g. `\Gamma_1=(1,0,0)`,
`\Gamma_2=(0,1,0.5)`) so the field is genuinely 3-D and `J` is generically full.

Evaluation points, for each regime:

- the two **particle centers** `X_1`, `X_2` (the regularization makes `J` finite
  here; the singular kernel would be the only one needing care), and
- a few **hundred random points** drawn in a box scaled to the particle cloud,
  including points inside the cores, in the overlap zone, and in the far field.

Pass criterion: `|\operatorname{tr} J| / \|J\|_F` (and the bare
`|\operatorname{tr} J|`) at the level of floating-point round-off,
`\lesssim 10^{-10}` in Float64 (cancellation across overlapping particles lifts
the floor a little above `\epsilon`), at *every* evaluation point in *every*
regime. A simultaneous **sanity assertion** confirms the test has teeth: the
curl `\|\nabla\times u\| = \|(\,J_{32}-J_{23},\,J_{13}-J_{31},\,J_{21}-J_{12}\,)\|`
is many orders of magnitude larger than the divergence, i.e. the field is not
trivially zero or irrotational — divergence is small because of structure, not
because nothing is happening.

To prove the diagnostic can fail, the script also runs a **deliberately broken
control**: a non-radial mollifier (anisotropic core, `g` evaluated on a
coordinate-scaled radius). That control should produce `\operatorname{tr} J`
well above round-off, demonstrating the metric is not vacuously zero.

### Run

ForwardDiff is a dependency of the **test** environment, so run the script under
that project:

```bash
julia --project=test examples/particle_divergence_free_check.jl
```

(The script imports only `ForwardDiff`/`SpecialFunctions`, not FLOWPanel itself,
and defines the kernel inline; it falls back to a self-contained `erf` series if
`SpecialFunctions` is unavailable.)

The script prints, per regime, max/mean `|\operatorname{tr} J|`, the relative
divergence, and the curl magnitude over the particle centers plus random-point
cloud, then a final `PASS`/`FAIL`. Output table is also written to
`data/particle_divergence_free/divergence_report.csv`.

## Acceptance

Technically complete when the script shows, for all three overlap regimes and
both at centers and at the random-point cloud:

- `|\operatorname{tr} J|/\|J\|_F \lesssim 10^{-10}` (Float64 round-off), and
- the irrotationality sanity check fails on purpose (curl well above divergence), and
- the broken-control case reports a divergence many orders of magnitude larger,

confirming the metric is meaningful. The expected outcome is a clean positive
result consistent with the theory; a *negative* result (non-zero divergence at
round-off-beating magnitude for the standard radial kernel) would instead be a
significant finding pointing at a kernel/kerneloffset bug and should be escalated
to items 004/008.

## Caveats

- The proof assumes the velocity is evaluated at the *same* argument the curl is
  taken with respect to. FLOWPanel's panel↔particle coupling uses a
  `kerneloffset` / control-point offset; if a downstream consumer differentiates
  the offset-evaluated velocity w.r.t. the *un*-offset point, the field it sees
  is still divergence-free in its own argument, but care is needed not to
  misattribute an offset chain-rule term as a divergence. The test deliberately
  uses the bare kernel to isolate the structural claim.
- This item certifies the *velocity* divergence only. It says nothing about
  whether the discrete strength-update conserves anything (item 008) or about
  vorticity-field accuracy. Keep the two questions separate.
- Float64 round-off floor scales with the dynamic range of `J`; very deep overlap
  (`d\ll\sigma`) with large strengths can lift the floor slightly. Report the
  relative metric, not just the absolute, so the round-off interpretation holds.
- AD at the particle center is subtle: the velocity is analytic there, but the
  naive `g(\rho)/r^3` form with `\rho=|r|` is `0/0` and `|r|` is not
  differentiable at `r=0`, so ForwardDiff returns `NaN`. The kernel must be
  written through `\rho^2=r\cdot r` with a small-radius Taylor branch
  `\phi(\rho)=g(\rho)/\rho^3=\sqrt{2/\pi}(1/3-\rho^2/10+\rho^4/56-\cdots)` to be
  AD-clean at centers (this is implemented in the script). This is a property of
  the *expression*, not the field.

## 2026-06-23 Verification Run

Ran `julia --project=test examples/particle_divergence_free_check.jl`. All three
overlap regimes pass; the broken anisotropic-core control fails on purpose,
confirming the metric is non-vacuous.

| regime | npts | max\|tr J\| | max rel div | mean\|curl\| | result |
| --- | --- | --- | --- | --- | --- |
| overlap a lot (d=0.25σ) | 402 | 1.9e-14 | 9.5e-12 | 8.0e-3 | PASS |
| overlap a little (d=2σ) | 402 | 5.8e-14 | 2.2e-11 | 8.2e-3 | PASS |
| far apart (d=20σ) | 402 | 9.4e-15 | 7.2e-12 | 2.5e-4 | PASS |
| **broken** anisotropic (d=2σ) | 402 | 3.3e-2 | **1.06e+0** | 6.9e-3 | div≠0 ✓ |

(npts = 2 particle centers + 400 random points spanning cores → far field.)

**Conclusion:** the Gaussian-erf vortex-particle velocity field is
divergence-free to floating-point round-off (relative divergence `~10^{-11}`) for
all overlap regimes and at every sampled point, including the particle centers —
exactly as the theory predicts, since each radial particle kernel is the curl of
a regularized vector potential. Meanwhile the curl is `O(10^{-2})`, ~9–12 orders
of magnitude above the divergence, so the field is genuinely rotational and the
near-zero divergence is structural, not vacuous. The deliberately non-radial
control gives relative divergence `~1`, proving the AD diagnostic has teeth. The
incompressibility assumption used by `PressureLaplace`, the unsteady Bernoulli
term, and item 008's `J\omega` algebra is therefore validated for the velocity
field; any observed dilatation downstream would point at an offset/chain-rule or
discretization artifact, not the kernel itself. Report `data/particle_divergence_free/divergence_report.csv`.

## 2026-06-24 Clear-Context Review

Reviewed the theory, script, and recorded output, then reran
`julia --project=test examples/particle_divergence_free_check.jl`. The rerun
matches the 2026-06-23 table exactly enough for approval:

- radial-kernel proof is pointwise and independent of particle overlap;
- ForwardDiff samples include both particle centers and 400 random points for
  each overlap regime;
- max relative divergence remains `2.3e-11`, below the stated `1e-10` threshold;
- the broken anisotropic control reports relative divergence `1.06`, so the
  diagnostic is not vacuous.

Approved. This item validates the velocity-field incompressibility assumption
for the radial Gaussian-erf particle kernel. It does not validate discrete
strength evolution or overlap-corrected vorticity recovery; those remain item
008 concerns.
