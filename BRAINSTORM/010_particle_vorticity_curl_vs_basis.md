# Particle Vorticity: Curl-of-J vs Zeta-Basis Evaluation

## Context

FLOWVPM particle states expose two different-looking estimates of the
particle-induced vorticity field at particle locations:

1. **Curl-of-gradient estimate** from the stored velocity gradient `J`:

   ```math
   \omega_J(x) = \nabla\times u(x),
   ```

   available cheaply from the FLOWVPM `velocity_gradient` layout
   (column-major `J_{ab}=\partial u_a/\partial x_b`):

   ```math
   \omega_J =
   \left(
       J_6-J_8,\,
       J_7-J_3,\,
       J_2-J_4
   \right).
   ```

   This is exactly the vorticity FLOWVPM's Pedrizzetti relaxation already reads
   (`relax_pedrizzetti`, `FLOWVPM_relaxation.jl:41`).

2. **Direct zeta-basis estimate** from the particle vorticity ansatz:

   ```math
   \omega_\zeta(x)
   =
   \sum_i \Gamma_i\,\zeta_{\sigma_i}(x-X_i),
   ```

   computed by the FLOWVPM basis evaluators `zeta_direct` / `zeta_fmm`
   (`FLOWVPM_viscous.jl:488`).

Item 008's saved-state check reported a nonzero difference between these two
(`basis ω vs curl-of-J ≈ 0.186` on the settled rotor-hover wake) and used it as a
**normalization calibration**. This item shows that interpretation is wrong: the
two estimates are **not** mathematically equal, the difference is a structural
consequence of the regularized-particle representation, and `0.186` is a physical
discrepancy — not slop confirming the kernel constant.

The questions this item settles are: *why* are they unequal (gauge / Helmholtz),
*how big* is it (a closed-form self-term, plus an open overlap question), and
*what can be done about it efficiently* while keeping the Gaussian kernel (whose
viscous core-spreading law `\sigma^2 \mathrel{+}= 2\nu\,dt` is especially clean,
`FLOWVPM_viscous.jl:160`).

## Theory

### Exact identity

Use the same regularized representation as items 008/009. Each particle has
position `X_i`, vector strength `\Gamma_i`, core size `\sigma_i`. The two kernels
are a matched Green's pair: with the mollified Laplace Green's function `G_\sigma`
and the radial mollifier `\zeta_\sigma`,

```math
-\nabla^2 G_{\sigma}(x-X_i) = \zeta_{\sigma}(x-X_i).
```

The regularized velocity is an **exact curl** of a mollified vector potential
(this is the FLOWVPM Biot–Savart kernel; see item 009 for the Gaussian-erf
form):

```math
u(x) = \nabla\times A(x),
\qquad
A(x)=\sum_i \Gamma_i\,G_{\sigma_i}(x-X_i).
```

Then, with no approximation,

```math
\omega_J(x)
= \nabla\times u
= \nabla\times\nabla\times A
= \nabla(\nabla\cdot A) - \nabla^2 A.
```

Since `\nabla^2 A=\sum_i\Gamma_i\nabla^2 G_{\sigma_i}=-\sum_i\Gamma_i\zeta_{\sigma_i}=-\omega_\zeta`,

```math
\boxed{\;
\omega_J(x) = \omega_\zeta(x) + \nabla(\nabla\cdot A)(x)
\;}
```

The two estimates coincide **only** if `\nabla(\nabla\cdot A)=0`. They do not in
general: `\nabla\cdot A=\sum_i\Gamma_i\cdot\nabla G_{\sigma_i}(x-X_i)\ne 0` for
constant vector strengths and a radial `G_\sigma`. VPM never imposes the Coulomb
gauge `\nabla\cdot A=0`, so `\omega_J\ne\omega_\zeta` is the generic situation.

This is exact for any matched radial pair, hence for FLOWVPM's Gaussian-erf
kernel `\zeta_{\mathrm{gauserf}}(r)=\mathrm{const1}\,e^{-r^2/2}`,
`\mathrm{const1}=1/(2\pi)^{3/2}` (`FLOWVPM_kernel.jl:50`). No sign or
normalization remains "to be verified": the relation follows from
`u=\nabla\times A` and `-\nabla^2 G=\zeta` alone.

### Solenoidality is the root cause

The cleanest way to see *what* the two fields are:

- `\omega_J=\nabla\times u` is **automatically divergence-free**, because it is a
  curl: `\nabla\cdot\omega_J\equiv 0` for any velocity field whatsoever.
- `\omega_\zeta=\sum_i\Gamma_i\zeta_i` is **not** divergence-free:

  ```math
  \nabla\cdot\omega_\zeta=\sum_i\Gamma_i\cdot\nabla\zeta_i\ne 0.
  ```

**The blunt version: a scalar Gaussian blob is not the curl of anything.** A field
can be written as `\nabla\times u` for some `u` *only if* it is divergence-free,
because `\nabla\cdot(\nabla\times u)\equiv 0` for every `u`. Since
`\nabla\cdot(\Gamma\,\zeta_\sigma)=\Gamma\cdot\nabla\zeta_\sigma\ne 0`, a single
scalar-Gaussian blob `\Gamma\,\zeta_\sigma` is **not the curl of any velocity
field at all** — and the same holds for the generic multi-particle sum
`\omega_\zeta`. The obstruction is not a gauge choice on `u` (you can always make
`u` incompressible); it is that the target vorticity itself is not solenoidal, so
no velocity has it as its curl. **This is not specific to the Gaussian:** any
scalar radial kernel times a constant vector has `\nabla\cdot(\Gamma\zeta)=
\Gamma\cdot\nabla\zeta\ne 0`. The culprit is the **scalar-kernel × constant-vector
ansatz**, not the Gaussian shape — which is why the fix does not require
abandoning the Gaussian: M4 keeps the Gaussian heat-kernel/core-spreading
structure (giving up only locality, not the Gaussian).

What FLOWVPM actually does is declare `\omega_\zeta` to be the vorticity and feed
it to regularized Biot–Savart. The resulting `u` *is* divergence-free (item 009),
but its true curl `\omega_J=\nabla\times u\ne\omega_\zeta`. So `u` is the velocity
of the *nearest solenoidal* vorticity field, the model carries `\omega_\zeta` as
strengths, and the physics (everything fed by `\nabla\times u`) silently uses
`\omega_J`.

Writing `\delta=\nabla(\nabla\cdot A)` and `\varphi=\nabla\cdot A`, the boxed
identity rearranges to

```math
\omega_\zeta = \omega_J - \nabla\varphi,
\qquad
\nabla^2\varphi = -\nabla\cdot\omega_\zeta.
```

This is exactly the **Helmholtz decomposition** of the blob field:
`\omega_J` is the solenoidal (divergence-free) part and `-\nabla\varphi` is the
irrotational part. Therefore:

> **`\omega_J` is the exact Helmholtz divergence-free projection of the blob
> vorticity `\omega_\zeta`. The discrepancy `\delta=\nabla(\nabla\cdot A)` is
> precisely the spurious irrotational vorticity that the scalar-blob ansatz
> carries but physical vorticity cannot have.**

This is consistent with — and complementary to — item 009. Item 009 proved
`\nabla\cdot u\equiv 0` (the divergence of a curl), which is *silent* about
`\nabla\times u`. The curl question is settled **here**: `\nabla\times u` does not
recover the blob vorticity; it recovers its solenoidal projection. (009's
"Relation to the regularized vorticity" pointer, which currently sends the reader
to 008, really points to this item.)

### The single-particle self-term (headline result)

The discrepancy is **not primarily an overlap effect**. It is already large for a
*single isolated particle*, evaluated at its own center.

For one particle, `\delta=\nabla(\nabla\cdot A)=(\Gamma\cdot\nabla)\nabla G_\sigma
=\mathrm{Hess}(G_\sigma)\,\Gamma`. The Hessian of a radial function at its center
is isotropic,

```math
\mathrm{Hess}(G_\sigma)(0)=\tfrac13\nabla^2 G_\sigma(0)\,I
=-\tfrac13\zeta_\sigma(0)\,I,
```

so at the particle center `X`:

```math
\boxed{\;
\omega_J(X)=\zeta_\sigma(0)\Gamma-\tfrac13\zeta_\sigma(0)\Gamma
=\tfrac{2}{3}\,\zeta_\sigma(0)\,\Gamma
=\tfrac23\,\omega_\zeta(X)
\;}
```

A **33% discrepancy at zero overlap**. The factor `2/3` is universal for any
matched radial pair (it uses only isotropy at the origin and `-\nabla^2 G=\zeta`),
hence exact for the Gaussian-erf kernel. The curl **redistributes** the blob's
vorticity (lowering the center, adding compensating lobes and a long projected
tail) rather than applying a local rescaling. The whole-space moment of that
projected field needs the M4 caveat below; do not infer from the local self-term
that `\int\mathrm{Hess}(G_\sigma)\Gamma\,dV=0`.

### A locality subtlety

It is tempting to say the discrepancy "is" `\nabla\cdot\omega_\zeta`. That is too
glib pointwise. The discrepancy is `\delta=\nabla(\nabla\cdot A)`, governed by the
gradient of the *potential's* divergence, with only the integral relation
`\nabla\cdot\delta=-\nabla\cdot\omega_\zeta`. At an isolated particle center
`\nabla\cdot\omega_\zeta=\Gamma\cdot\nabla\zeta(0)=0` (radial kernel peaks there),
yet `\delta\ne 0`. The numerical check below must therefore compare `\delta`
against the analytic `\mathrm{Hess}(G_\sigma)\Gamma`, not against
`\nabla\cdot\omega_\zeta`.

### Why VPM still works despite this

In a converged, well-overlapped discretization, the blob field `\omega_\zeta`
approximates a smooth, genuinely divergence-free vorticity field, so
`\nabla\cdot\omega_\zeta\to 0` and `\omega_J\approx\omega_\zeta`. The isolated
particle is the **maximally under-resolved** case: it has no neighbors to "fill
in" the redistributed center, so it shows the worst self-term. Hence the
discrepancy `\delta` is large wherever `\omega_\zeta` is far from a smooth
solenoidal field — under-resolved regions, heavy *non-aligned* overlap, or noisy
strengths — which is exactly the settled rotor-hover regime examined in item 008
(≈22.5k neighbors per particle, non-diagonally-dominant `M`, non-parallel
strengths).

### Answers to the original questions

- **Equivalent under any common assumption?** No. They differ by
  `\nabla(\nabla\cdot A)` exactly; equality needs the Coulomb gauge, which VPM
  does not enforce.
- **Is the difference the projection term or a convention issue?** The projection
  term `\nabla(\nabla\cdot A)`, derived above with fixed signs — not an
  implementation convention.
- **Does it vanish in the no-overlap limit at a center?** **No.** A single
  regularized particle already shows the closed-form `2/3` self-term at its
  center.
- **Effect of core size / off-center / non-parallel strengths?** The `2/3`
  center factor is independent of `\sigma`; off-center, `\delta` follows
  `\mathrm{Hess}(G_\sigma)\Gamma` (sign-changing lobes); overlap and
  non-parallel strengths add neighbor Hessian contributions on top of the self
  term.
- **Which estimate should item 008 use?** See the directive/discussion below —
  this is deliberately left as a question to explore, because the self-vs-overlap
  split in a real wake is not yet quantified.

## Relationship to item 008 — to explore, not yet decided

Item 008's overlap equation is posed in the `\zeta` basis: its left-hand side is
`\sum_i\zeta_i(x)\dot\Gamma_i`, so the *consistent* vorticity to use on its
right-hand side is the basis field `\omega(X_k)=\sum_j\Gamma_j\zeta_j(X_k)=
\omega_\zeta`. Using `\omega_J` there would be inconsistent with the basis on the
left. That argues 008 should keep `\omega_\zeta` as its operational vorticity and
treat `\omega_J` only as a divergence diagnostic / projection target.

What must **not** be carried forward from 008 is the claim that
`basis ω vs curl-of-J = 0.186` "confirms normalization." By the identity above,
that `0.186` is the physical self-plus-overlap discrepancy
`\|\nabla(\nabla\cdot A)\|/\|\omega_\zeta\|`, not kernel-constant slop.

However — and this is the open point — we have **not** decomposed how much of a
real-wake `\omega_J-\omega_\zeta` gap is the isolated self-term versus genuine
inter-particle overlap. At an isolated center it is *pure* self-term (the
closed-form `1/3`); in the settled hover wake the two are mixed and the `0.186`
sits *below* the `1/3` single-center value, suggesting partial cancellation from
neighbors and/or field-norm averaging over off-center points. **Quantifying this
split is future work** (a decomposition experiment: per-particle self-term
`\mathrm{Hess}(G_{\sigma_i})\Gamma_i` versus the neighbor remainder over a saved
state). Until it is done, 010 records the *interpretation* correction to 008 but
does not issue a hard directive on which field 008's solver must use.

**M2 update (2026-06-25): M2 is not a replacement for 008 — it is the same overlap
solve with `\Gamma_{\mathrm{eff}}` coefficients.** Now that M2 is worked out, the
"which field does 008 use" tension above has a clean resolution, and the three tracks
form one machine with different coefficients:

- **008** is the per-step overlap solve `M\dot\Gamma=b` in the scalar `\zeta` basis,
  solving for the strength *rate* with overlap accounted for. Per the consistency
  argument above, as posed it must keep `\omega_\zeta` (not `\omega_J`) on its RHS.
- **M2's full law** `M\dot\Gamma_{\mathrm{eff}}=b_{\mathrm{eff}}` is *that same overlap
  solve* — same `M`, same neighbor sums, same CG — but with the coefficients
  reinterpreted as `\Gamma_{\mathrm{eff}}` so the represented field is the solenoidal
  `\omega_J=M\Gamma_{\mathrm{eff}}`. This is the **legitimate** way to put `\omega_J`
  into an 008-style scalar solve: not by swapping the RHS (which the LHS basis forbids),
  but by retargeting the carried strength. M2 therefore **refines** 008, reusing all of
  its machinery (collocation matvec, neighbor lists, residual diagnostic), rather than
  replacing it.
- **M2 option (1)** (drift-corrected self-only stretching) is *not* a cheaper bypass
  around the 008 solve: its `\dot T\,\Gamma` term contains `M^{-1}`, so in the
  heavy-overlap wake it inherits the same per-step CG and collapses to ~full-solve cost
  (see the cost analysis in the M2 section). Once you are paying for `M^{-1}` anyway you
  would solve the exact `M\dot\Gamma_{\mathrm{eff}}=b_{\mathrm{eff}}` rather than the
  approximation — option (1) is dominated (same cost, less exact). It is a conceptual
  bridge, not a shippable shortcut. The only genuinely *cheaper-than-008* tier is the
  self-only `\tfrac23\Gamma`+relaxation scheme, which is approximate (drifts).
- **M6** is what genuinely *supersedes* 008 in spirit: the same overlap solve redone
  exactly in the divergence-free matrix basis `K_\sigma` — "the 008-style overlap-aware
  strength update in that basis."

So: would option (1) replace 008? No. 008 is the overlap-solve framework; M2 is that
framework with `\Gamma_{\mathrm{eff}}` coefficients (the correct route for `\omega_J`),
and M6 is its exact matrix-basis successor. The cheap-and-approximate end is self-only;
the accurate-and-expensive end is the 008/M6 solve. (Light-overlap caveat: where `M^{-1}`
truncates in a few terms, option (1) *would* be a cheap consistent alternative to a full
008 solve — but the settled hover wake is the expensive end.)

## Numerical Task (N=1)

A single-particle synthetic check is sufficient to confirm the headline result
unambiguously and is independent of rotor-hover data. Create
`examples/particle_vorticity_curl_vs_basis_check.jl`, reusing the AD-clean
Gaussian-erf kernel infrastructure from item 009's
`examples/particle_divergence_free_check.jl` (in particular the `\rho^2` Taylor
branch that makes `g(\rho)/\rho^3` differentiable at the center).

### Setup

One particle: `\sigma=1`, `\Gamma=(1,\,0.2,\,0)`. Evaluate at:

- the particle center `X` (the decisive point), and
- a short radial line sample `r/\sigma\in[0,4]` to expose the spatial structure
  of `\delta`.

### Computed quantities

For each evaluation point compute:

```math
\omega_\zeta=\Gamma\,\zeta_\sigma(x-X),
\qquad
\omega_J=\nabla\times u
\ \text{(ForwardDiff through the Gaussian-erf velocity kernel)},
```

```math
\delta=\omega_J-\omega_\zeta,
\qquad
\delta_{\text{analytic}}=\mathrm{Hess}(G_\sigma)(x-X)\,\Gamma,
\qquad
\nabla\cdot\omega_J=\operatorname{tr}\nabla\omega_J.
```

### Acceptance gates (these give the test teeth)

- **Self-term:** `\|\omega_J-\omega_\zeta\|/\|\omega_\zeta\| \to 1/3` at the
  center (closed-form prediction `\approx 0.3333`).
- **Structure:** `\delta` matches `\delta_{\text{analytic}}=
  \mathrm{Hess}(G_\sigma)\Gamma` to floating-point round-off along the line —
  confirming the discrepancy is the projection term, not numerical noise.
- **Solenoidality:** `\nabla\cdot\omega_J\approx` round-off (ties to item 009;
  the discrepancy lives entirely in the curl, not the divergence).

Run under the test project (ForwardDiff/SpecialFunctions live there):

```bash
julia --project=test examples/particle_vorticity_curl_vs_basis_check.jl
```

## Model-Improvement Brainstorm

The blob field carries a spurious irrotational vorticity `-\nabla\varphi`. The
goal is a **computationally efficient** way to make the particle model
self-consistent with respect to it, **preferably keeping the Gaussian** so the
viscous core-spreading model stays clean.

A shared efficiency fact unlocks all options: because `A` is itself a particle
sum of Gaussians, *every* Helmholtz piece — `\nabla\cdot A`,
`\nabla(\nabla\cdot A)`, `\varphi`, `\nabla\varphi`, and the per-particle
`\mathrm{Hess}(G_\sigma)` — is an **analytic Gaussian kernel sum on the same
particles**. So the correction is `O(N)` / FMM-able and needs **no auxiliary grid
and no global Poisson solve**.

This item explores **M1 (the divergence-correction term)**, **M2 (effective
strength `Γ_eff`)**, **M4 (the Leray-projected Gaussian blob)**, and **M6
(projected-Gaussian overlap evolution)** in depth. M1 is the diagnostic /
physics-term route that leaves the carried strengths untouched, M2 is the
coefficient-only route, M4 is the structural representation fix, and M6 is the
operational synthesis: use M4's projected basis as the represented vorticity and
re-derive the 008-style overlap-aware strength update in that basis. The correct
divergence-free basis (the projected blob
`K_\sigma=\zeta_\sigma I+\mathrm{Hess}(G_\sigma)`, which equals `\omega_J`)
already hides inside the model; M1 evaluates it on demand as a correction to
`\omega_\zeta`, M2 tries to fake it with the existing scalar kernel and adjusted
coefficients, M4 adopts it outright, and M6 asks how to evolve it.

### M1 (explored): the divergence-correction term

**Definition.** Keep the carried strengths `\Gamma_i` exactly as the model already
evolves them, and — wherever a *true* (solenoidal) vorticity is required — evaluate

```math
\omega_J(x)
=\omega_\zeta(x)+\sum_i \mathrm{Hess}(G_{\sigma_i})(x-X_i)\,\Gamma_i
=\sum_i K_{\sigma_i}(x-X_i)\,\Gamma_i,
\qquad
K_\sigma=\zeta_\sigma I+\mathrm{Hess}(G_\sigma),
```

as one extra analytic Gaussian-erf sum on the same particles. M1 is the
**diagnostic face of M4's `K_\sigma`**: it adds the anisotropic Hessian part to the
isotropic `\zeta_\sigma I` already carried, recovering the exact Leray-projected
field `\omega_J=\nabla\times u` without introducing any new state or changing the
strength evolution. It is the simplest member of the `K_\sigma` family: M2 best-fits
`\omega_J` with the scalar kernel and adjusted coefficients, M4 adopts `K_\sigma` as
the carried basis, M6 evolves it; M1 just *evaluates* it where the physics asks for
a divergence-free vorticity.

**Closed-form analytic Hessian (no ForwardDiff, no grid).** For the radial
Gaussian-erf Green's function `G_\sigma(r)=\mathrm{erf}(r/(\sqrt2\,\sigma))/(4\pi r)`,
write `y=x-X`, `s=|y|^2/\sigma^2`, and let `\phi(s)=g(\rho)/\rho^3`
(`\rho=\sqrt s`, `g(\rho)=\mathrm{erf}(\rho/\sqrt2)-\sqrt{2/\pi}\,\rho\,e^{-\rho^2/2}`)
be the same radial profile FLOWVPM's velocity kernel already uses. Then

```math
\boxed{\;
\mathrm{Hess}(G_\sigma)(y)\,\Gamma
= A\,\Gamma + B\,(y\!\cdot\!\Gamma)\,y,
\qquad
A=-\frac{\phi(s)}{4\pi\sigma^3},
\qquad
B=-\frac{\phi'(s)}{2\pi\sigma^5},
\;}
```

with the radial derivative

```math
\phi'(s)=\frac{\sqrt{2/\pi}\,\rho^3 e^{-\rho^2/2}-3\,g(\rho)}{2\rho^5}
\quad\Bigl(\text{Taylor at }s\!\to\!0:\;
\phi'(s)=\sqrt{\tfrac{2}{\pi}}\bigl(-\tfrac1{10}+\tfrac{s}{28}-\tfrac{s^2}{144}+\cdots\bigr)\Bigr).
```

Two closed-form consistency checks fix all signs/constants:

- **Center / self-term.** At `y=0`, `A(0)=-\tfrac13\zeta_\sigma(0)` and the `B` term
  vanishes, so `\mathrm{Hess}(G_\sigma)(0)\Gamma=-\tfrac13\zeta_\sigma(0)\Gamma` and
  `K_\sigma(0)=\tfrac23\zeta_\sigma(0)I` — exactly the `2/3` self-term confirmed by
  the N=1 check.
- **Trace.** `\mathrm{tr}\,\mathrm{Hess}(G_\sigma)=\nabla^2 G_\sigma=-\zeta_\sigma`,
  so the isotropic part of `K_\sigma` reproduces `-\nabla^2 G=\zeta` and the
  represented `\omega_J=\sum_i K_{\sigma_i}\Gamma_i` is the solenoidal projection.

The Taylor branches near `s=0` (mirroring `phi_gauserf`'s) keep the kernel clean at
and near the particle center, where the naive `g/\rho^3`, `g/\rho^5` forms are
`0/0`. **Numerically verified:** the closed form matches the ForwardDiff
`\mathrm{Hess}(G_\sigma)\Gamma` (`delta_analytic`) to `\le 8\times10^{-13}` relative
error along the radial line (`M1_SELFTEST=1
examples/particle_vorticity_curl_vs_basis_check.jl`), so M1 needs **no automatic
differentiation in production** — only `\mathrm{erf}` and `\exp`.

**Cost / reuse.** M1 is one additional Gaussian-erf particle sum over the same
particles (and the same neighbor lists / FMM tree the velocity evaluation already
uses). The kernel `\mathrm{Hess}(G_\sigma)` is the *same analytic Gaussian Hessian
FLOWVPM's FMM already forms for the velocity gradient* `J`
(`FLOWVPM_fmm.jl:114`), so M1 is `O(N)`/FMM-able with **no auxiliary grid and no
global Poisson solve**. The minimal implementation is the brute-force form in
`examples/particle_vorticity_curl_vs_basis_check.jl`
(`hess_green_times_gamma`, `omega_J_corrected`); the production form swaps the
inner loop for the hard-cutoff neighbor list / FMM Hessian sum already built in
`examples/particle_overlap_residual.jl`.

**Where it applies — and where it does not.** M1 is for any consumer that needs a
*true solenoidal* vorticity from a saved/live particle state: divergence
diagnostics, pressure / unsteady-Bernoulli source terms, force cross-checks, and
the item-008 overlap residual when a divergence-free target is wanted. It is **not**
needed by relaxation/stretching, which already read `\omega_J=\nabla\times u`
directly from the stored `J` (`relax_pedrizzetti`, `FLOWVPM_relaxation.jl:41`; the
`J`-stretching in `FLOWVPM_timeintegration.jl`) — those paths use the projected
vorticity already; the gap is only between that and the *stored* basis field
`\omega_\zeta` (`VORTICITY_INDEX`) that monitors/post-processing read.

**Self-vs-overlap decomposition (the item's open question — M1 is its vehicle).**
Because M1 forms the full sum `\sum_i\mathrm{Hess}(G_{\sigma_i})\Gamma_i`
explicitly, it splits cleanly at each evaluation particle `k` into the closed-form
**self-term** `\mathrm{Hess}(G_{\sigma_k})(0)\Gamma_k=-\tfrac13\zeta_{\sigma_k}(0)\Gamma_k`
and the **neighbor remainder** `\sum_{i\ne k}\mathrm{Hess}(G_{\sigma_i})(X_k-X_i)\Gamma_i`.
Quantifying the two contributions over a settled rotor-hover state — and thereby
explaining why the item-008 `0.186` sits *below* the isolated `1/3` (neighbor
partial cancellation vs field-norm averaging) — is the deferred numerical test for
M1; the `.vtp` loader and neighbor lists in
`examples/particle_overlap_residual.jl` already supply the machinery.

**Viscous compatibility.** `\mathrm{Hess}(G_\sigma)` is a Gaussian-derivative
kernel, so the clean core-spreading law `\sigma^2\mathrel{+}=2\nu\,dt` applies to it
unchanged (the same Fourier `P(k)`-commutes-with-`e^{-\nu|k|^2t}` argument as M4).
M1 evaluates a correction; it adds no diffusion bookkeeping of its own.

**Limitation (state, do not hide).** M1 makes the vorticity *that is evaluated*
solenoidal at the query points, but the carried strengths and their evolution are
untouched — the stored `\omega_\zeta` and the dynamics that stretch it remain
non-solenoidal. M1 is therefore a **diagnostic / physics-term correction**, not a
self-consistent model fix; closing the loop in the evolution is exactly M2's
best-fit retargeting and M4/M6's exact projected-basis dynamics. M1's value is that
it is the cheapest, zero-state way to *consume* the correct `\omega_J` wherever the
downstream physics needs it, and the closed-form kernel it provides is the shared
building block all of M2/M4/M6 reuse.

### M2 (explored): evolve an effective strength `Γ_eff` directly

**Definition.** Choose effective coefficients `\Gamma_{\mathrm{eff},i}` whose
basis sum best-represents the *solenoidal* vorticity:

```math
\sum_i \Gamma_{\mathrm{eff},i}\,\zeta_i(x)\ \approx\ \omega_J(x)
=\sum_j\bigl(\zeta_j I+\mathrm{Hess}(G_\sigma)_j\bigr)\Gamma_j.
```

In Galerkin/collocation form with overlap matrix
`M_{ki}=\zeta_{\sigma_i}(X_k-X_i)` and the analytic Hessian-overlap operator
`H_{kj}=\mathrm{Hess}(G_{\sigma_j})(X_k-X_j)`,

```math
M\,\Gamma_{\mathrm{eff}}=(M+H)\,\Gamma
\quad\Longrightarrow\quad
\Gamma_{\mathrm{eff}}=\Gamma+M^{-1}H\,\Gamma.
```

**Honest representation caveat (state, do not hide).** The scalar-`\zeta` ×
vector-`\Gamma` basis sum is *itself* non-solenoidal, with the same
`\nabla(\nabla\cdot A)` structure — a scalar Gaussian blob is not the curl of
anything. So `\Gamma_{\mathrm{eff}}` makes the field a **best fit** to `\omega_J`,
not an exactly divergence-free field; exactness requires the divergence-free
matrix basis developed as M4 below. M2's payoff is therefore **evolution
self-consistency**, not a perfectly clean field:
today FLOWVPM stretches the coefficients of `\omega_\zeta` using the `J` whose
curl is `\omega_J\ne\omega_\zeta`. Evolving `\Gamma_{\mathrm{eff}}` — a strength
that represents `\omega_J` — under that same `J` removes the mismatch between the
stretched quantity and the stretching operator.

**Operational contract (important): M2 is a two-field, auxiliary-strength model by
default.** The physical particle strength `\Gamma` remains the circulation vector and
the source used by regularized Biot–Savart to generate `u` and `J`. The effective
strength `\Gamma_{\mathrm{eff}}` is an auxiliary coefficient set used to represent the
projected vorticity in the scalar `\zeta` basis for RHS/stretching diagnostics and any
future projection/core-reset experiment. Replacing `\Gamma` with
`\Gamma_{\mathrm{eff}}` in the velocity solve is **not** part of M2 as written: it
would change the velocity field, change the fixed point whose curl is being fitted, and
can shift discrete total circulation because `\sum_i(M^{-1}H\Gamma)_i` need not vanish
even though the continuous `\int\mathrm{Hess}(G)\Gamma\,dx=0`. Any
`Γ→Γ_eff` replacement is a separate model requiring its own velocity/circulation proof.
With the auxiliary contract, "drift" means `\Gamma_{\mathrm{eff}}` drifting away from
the instantaneous projection of the live physical `\Gamma`; deciding whether stretching
advances `\Gamma`, `\Gamma_{\mathrm{eff}}`, or both is the unresolved live-design
question.

**Re-derive `\dot\Gamma_{\mathrm{eff}}` directly (the key task).** Rather than
carry `\Gamma`, project to `\Gamma_{\mathrm{eff}}` each step, and differentiate,
work the entire governing equation in terms of the effective strength: derive its
transport directly from the vorticity equation applied to the solenoidal field
`\omega_J` and projected onto the `\zeta` basis. The aim is to bake the
consistency into the time-derivative `\dot\Gamma_{\mathrm{eff}}` so that **no
separate relaxation or projection pass is needed** — the corrected dynamics are
the dynamics.

**Carried-out derivation.** Represent the *physical solenoidal* vorticity in the
scalar `\zeta` basis with the effective coefficients,
`\omega_J(x,t)\approx\sum_i\Gamma_{\mathrm{eff},i}(t)\,\zeta_{\sigma_i}(x-X_i(t))`,
and apply the inviscid transport `D\omega/Dt=(\omega\cdot\nabla)u` to that field.
Taking the material derivative (particles convecting with their own velocity `u_i`)
and collocating at the moving centers `X_k` gives the same basis-motion structure as
item 008's scalar overlap equation, but with `\omega_J` as the transported field:

```math
\boxed{\;
M\,\dot\Gamma_{\mathrm{eff}}=b_{\mathrm{eff}},
\qquad
b_{\mathrm{eff},k}
= \bigl[(\omega_J\!\cdot\!\nabla)u\bigr](X_k)
-\sum_i\Gamma_{\mathrm{eff},i}\,(u(X_k)-u_i)\!\cdot\!\nabla\zeta_i(X_k)
\;}
```

with `M_{ki}=\zeta_{\sigma_i}(X_k-X_i)`. The stretching term uses the *represented*
vorticity `\omega_J=M\Gamma_{\mathrm{eff}}`: with the FLOWVPM layout
`J_{ab}=\partial u_a/\partial x_b`,
`[(\omega_J\!\cdot\!\nabla)u](X_k)=J(X_k)\,\omega_J(X_k)=J(X_k)\,(M\Gamma_{\mathrm{eff}})_k`
(classic scheme; the transposed FLOWVPM default uses `J^\top`, the same toggle as the
`_euler` stretching `G[1]+=dt(J_1 G_1+J_4 G_2+J_7 G_3)` vs the transposed branch in
`FLOWVPM_timeintegration.jl`). The single subtraction term is the convection of the
basis through the nonuniform velocity field.

**Why there is no core-spreading term in the strength equation.** Viscous diffusion
is *not* a separate forcing of `\dot\Gamma_{\mathrm{eff}}`: it is the core-spreading
`\sigma`-update applied with the strengths held fixed (FLOWVPM's split — stretch, then
`viscousdiffusion`). The Gaussian basis is a heat kernel, so `\sigma^2\mathrel{+}=2\nu\,dt`
gives exactly `\dot\sigma_i\,\partial_{\sigma_i}\zeta_i=\nu\nabla^2\zeta_i`, i.e. the
basis-drift from spreading *is* the modeled `\nu\nabla^2\omega`. In the full viscous
transport `D\omega/Dt=(\omega\cdot\nabla)u+\nu\nabla^2\omega` that drift cancels the
diffusion source term, so it never enters `b_{\mathrm{eff}}` (inviscidly, `\dot\sigma=0`
and the question is moot). Sanity check: under pure diffusion (`J=0`, uniform `u`) the
strength must be conserved, `\dot\Gamma_{\mathrm{eff}}=0`; the cancellation gives exactly
that, whereas keeping a standalone `-\sum_i\Gamma_{\mathrm{eff},i}\dot\sigma_i\partial_\sigma\zeta_i`
term would spuriously grow circulation by `+3(\dot\sigma/\sigma)\Gamma`.

**No-overlap reduction (the teeth).** In the self-dominant limit
`M\approx\mathrm{diag}(\zeta_{\sigma_k}(0))` the basis-motion terms vanish (a particle
does not convect relative to itself) and the overlap matrix cancels on both sides,
leaving the *exact functional form of the present isolated stretching* now applied to
the effective strength:

```math
\dot\Gamma_{\mathrm{eff},k}=J_k\,\Gamma_{\mathrm{eff},k}
\qquad(\text{self-only; classic}).
```

Together with the self-only projection `\Gamma_{\mathrm{eff},k}\approx\tfrac23\Gamma_k`
this is the cheap closed-form scheme: a scalar `2/3` rescale at initialization plus
**unchanged** `J`-stretching — cheaper than Pedrizzetti's per-particle
rotate-and-renormalize, and self-consistent *in the isolated limit* (the stretched
quantity `\Gamma_{\mathrm{eff}}` represents the same `\omega_J=\nabla\times u` whose
gradient `J` does the stretching). **Caveat, now backed by data:** the 2026-06-25
saved-state diagnostic finds that in the settled rotor-hover wake the neighbor Hessian is
`\approx3.1\times` the self-term and only weakly aligned with it (`\cos\approx0.34`), so
the self-only `\tfrac23` rescale is **not** a useful standalone model there — it is
actually *worse* than carrying raw `\Gamma` (residual `0.295` vs `0.186`). The cheap limit
is a clean isolated-particle result, not a wake-ready scheme; in heavy overlap a real
local/global projection (or M4/M6) is required. See that dated section.

**Projection-preservation subtlety (state honestly).** The cheap self-only update
`\dot\Gamma_{\mathrm{eff},k}=J_k\Gamma_{\mathrm{eff},k}` is also drift-prone, and it
is worth seeing exactly why, because it is the crux of whether M2's "the corrected
dynamics *are* the dynamics" aspiration holds.

*Two things move at once.* Write the projection as
`\Gamma_{\mathrm{eff}}=(I+T)\Gamma` with `T=M^{-1}H`. The operator `T` is **not
constant**: both `M_{ki}=\zeta_{\sigma_i}(X_k-X_i)` and
`H_{kj}=\mathrm{Hess}(G_{\sigma_j})(X_k-X_j)` depend on the *instantaneous geometry*
— the positions `X_i` and the core sizes `\sigma_i`. As particles convect through
strain and their cores spread (`\sigma^2\mathrel{+}=2\nu\,dt`), the geometry changes,
so `T=T(t)` and `\dot T\ne 0`.

*Ordinary stretching is blind to `\dot T`.* The standard VPM update `\dot\Gamma=S\Gamma`
(with `S` the per-particle `J`-contraction) advances the *coefficients* but treats
each blob as a fixed-shape stamp; it carries no information about `T`. Ask what
evolution `\Gamma_{\mathrm{eff}}` must obey for the relation `\Gamma_{\mathrm{eff}}=(I+T)\Gamma`
to stay true. Differentiating,

```math
\dot\Gamma_{\mathrm{eff}}
=(I+T)\dot\Gamma+\dot T\,\Gamma
=(I+T)\,S\,\Gamma+\dot T\,\Gamma,
```

and re-expressing in the carried state `\Gamma=(I+T)^{-1}\Gamma_{\mathrm{eff}}`,

```math
\dot\Gamma_{\mathrm{eff}}
=\underbrace{(I+T)\,S\,(I+T)^{-1}\Gamma_{\mathrm{eff}}}_{\ne\,S\,\Gamma_{\mathrm{eff}}\text{ unless }[S,T]=0}
+\;\dot T\,(I+T)^{-1}\Gamma_{\mathrm{eff}}.
```

The consistent law therefore differs from "just stretch `\Gamma_{\mathrm{eff}}`" by
two pieces: a **conjugation mismatch** `[(I+T)S(I+T)^{-1}-S]` (zero only when `S` and
`T` commute, which they do not in general) and an explicit **geometry-drift term**
`\dot T\,\Gamma`. Drop them — as the cheap self-only update does — and a
`\Gamma_{\mathrm{eff}}` that starts as the exact projection of `\Gamma` is, one step
later, **no longer** the best fit to `\omega_J` for the cloud's new arrangement: the
projection operator moved underneath the update and the cheap law ignored it. (This is
the same kind of representation degradation that forces FLOWVPM to relax / core-reset
at all; M2's drift is that phenomenon in the `\Gamma_{\mathrm{eff}}` map.)

*Two consistent realizations — pay the cost continuously or in bulk.*
1. **Carry the drift correction `\dot T\,\Gamma`** (the dominant, intuitive piece;
   the conjugation term is the same order). `\dot T=\tfrac{d}{dt}(M^{-1}H)` contains
   `M^{-1}`, which couples every particle; if it is finite-differenced against the
   previous step's stored projection it is nearly free *given* the projection, but you
   still owe the `M^{-1}H\Gamma` solve each step to have that projection. Whether `\dot T`
   can be cheaply localized depends on the overlap regime (see cost below): it can in
   lightly-overlapped regions, but **not** in the heavy-overlap settled wake.
2. **Re-project at the core-reset cadence.** Let `\Gamma_{\mathrm{eff}}` drift under
   ordinary stretching and periodically re-solve `\Gamma_{\mathrm{eff}}=\Gamma+M^{-1}H\Gamma`
   from scratch, wiping the accumulated drift. FLOWVPM **already runs a periodic global
   RBF solve at exactly this cadence** (`rbf_conjugategradient`, the core-spreading
   core-reset), so this piggybacks on existing machinery and amortizes the expensive
   global solve over many cheap steps.

*Cost of the two realizations (overlap-dependent — the punchline corrects a naive
"option (1) is cheap" reading).* Measure work in **near-field kernel matvecs**: one pass
over the neighbor lists, `\sim N\bar n` kernel evaluations (`\bar n` = average neighbors
in the `\sim4\sigma` cutoff), roughly the near/P2P portion of the one FMM
velocity+gradient evaluation already paid each step. The cost is governed entirely by the
overlap `\bar n` and the conditioning of `M`; item 008's settled rotor-hover state pins
both: `\bar n\approx 2.25\times10^4` (within the `4\cdot\sigma_{\max}` cutoff) and a
diagonal-dominance ratio `\approx 0.006` (so `M` is *far* from diagonally dominant).

> **Cost-count caveat (item 011 overlap re-characterization).** The `\bar n\approx2.25\times10^4`
> here is the neighbor count within `4\cdot\sigma_{\max}`, and `\sigma_{\max}` is the *global*
> max core size (`\sigma_{\max}/\sigma_{\min}\approx15`), so that cutoff `\approx1.7\,R`
> engulfs half the wake. The **true** overlap is far smaller: `\sigma/h\approx4.2`,
> `\approx394` within `2\sigma_{\mathrm{local}}`, kernel-effective `\approx181` (item 011,
> `eff_overlap_sample`). The **relative** cost ordering below is unaffected — it is set by the
> CG *iteration* count (conditioning), not `\bar n` — but the **absolute** per-pass cost is
> pessimistic: a `\sigma_{\mathrm{local}}`-based cutoff would cut the P2P work `\sim50\text{–}100\times`
> while leaving the operators (and the conditioning verdict) essentially unchanged, because
> the far ζ entries it drops are `\sim e^{-8}` and negligible. So the "heavy cost" is partly a
> cutoff-choice artifact; the intrinsic obstruction is the conditioning, not the count.

- **Free tier — self-only `\tfrac23\Gamma` + standard stretching:** `\sim 0` extra units,
  cheaper than Pedrizzetti's rotate-and-renormalize, but drifts.
- **Option (1), carry `\dot T\,\Gamma`:** the killer is `M^{-1}`. Forming the projection
  once costs one `H\Gamma` Hessian neighbor sum (`\sim1` unit) **plus a CG solve**
  (`rbf_conjugategradient` runs `\mathrm{itmax}\approx15`, tol `10^{-3}` ⇒ `\sim15`
  units); the conjugation `[S,T]` piece needs `T` applied a few more times, roughly
  doubling it. So option (1) is `\sim25\text{–}30` near-field passes **per step**, order
  `10\times` a baseline VPM step — the *same order as just doing the full M6/008 per-step
  solve*. The "neighbor-truncate `M^{-1}`" shortcut is precisely what fails here: a
  local/Jacobi/short-Neumann inverse only converges when `M` is near-diagonally-dominant,
  and at ratio `0.006` it is not. (In lightly-overlapped regions `M^{-1}` truncates in a
  few terms and option (1) genuinely *is* a couple of cheap passes — but the wake we care
  about is the expensive end.)
- **Option (2), re-project every `K` steps:** the same `\sim16`-unit solve amortized,
  `\sim16/K` units/step, reusing FLOWVPM's existing periodic RBF core-reset; drift
  accumulates between resets.

So the honest ordering in the heavy-overlap regime is: option (1) `\approx` full solve
(not cheap); option (2) is the cost-sensible realization (amortized); self-only is the
free tier but, per the saved-state data below, *inaccurate* there (worse than raw — not
just drift-prone). The verdict is regime-dependent; the 2026-06-25 saved-state diagnostic
(see that section) measured `\bar n\approx2.25\times10^4`, diagonal-dominance `\approx0.006`,
and a neighbor Hessian `\approx3.1\times` the self-term — confirming these
order-of-magnitude figures and that a truncated-Jacobi `M^{-1}` is unstable (`\sim27\times`
row/probe residual). A true local projection fits its own collocation rows to round-off,
but held-out probes remain poor (`\sim0.54` residual), so the local solve is a solve check
rather than live-model evidence.

**Numerical confirmation of the projection.** The static projection
`\Gamma_{\mathrm{eff}}=\Gamma+M^{-1}H\Gamma` is implemented and checked in
`examples/particle_vorticity_curl_vs_basis_check.jl` (`M2_SELFTEST=1`), reusing the M1
Hessian kernel for `H`: (1) a single particle collapses to
`\Gamma_{\mathrm{eff}}=\tfrac23\Gamma` to round-off; (2) two overlapping
equal-`\sigma` particles agree between the dense solve and a CG mirroring
`rbf_conjugategradient`, with collocation residual `\sim10^{-17}`; (3) **off-node**,
the `\Gamma_{\mathrm{eff}}` field is `\approx 8.6\times` closer to `\omega_J` than the
raw `\omega_\zeta` field (relative residual `0.527\to0.061`), with a nonzero floor; and
(4) the field divergence drops `0.642\to0.360` — the direct test of the "best-fit a
*solenoidal* field" claim. Note (4)'s reduction (`\approx1.8\times`) is far milder than
(3)'s closeness gain (`8.6\times`): the scalar best-fit moves the field toward `\omega_J`
but leaves a substantial irreducible divergence floor — exactly the honest caveat, and a
concrete reason the *exact* fix needs M4's matrix basis (closing the floor to zero), not
better scalar coefficients.

**Hunt for a cheap closed-form limit competitive with Pedrizzetti.** The
self-term alone gives a per-particle rescale

```math
\Gamma_{\mathrm{eff},i}\approx \tfrac23\,\Gamma_i
\quad\text{(self-only, collocation at the center)},
```

which is `O(1)` per particle — *cheaper* than Pedrizzetti's per-particle
rotation. The open quantitative question is how much the neighbor-overlap Hessian
sum `M^{-1}H` adds on top, and whether a neighbor-truncated `M^{-1}H` (the same
neighbor lists 008 already builds) keeps `\dot\Gamma_{\mathrm{eff}}` at or below
Pedrizzetti cost. Target deliverable: a per-particle or short-neighbor-list update
rule for `\dot\Gamma_{\mathrm{eff}}`.

**Reuse existing machinery.**
- The projection solve `M\Gamma_{\mathrm{eff}}=(M+H)\Gamma` maps directly onto
  FLOWVPM's RBF/CG core-reset solver (`rbf_conjugategradient`,
  `FLOWVPM_viscous.jl:314`), which already re-fits `\Gamma` to a target
  vorticity — today targeting `\omega_\zeta`; M2 retargets it to `\omega_J`.
- **Pedrizzetti relaxation already rotates `\Gamma` toward `\nabla\times u=
  \omega_J`** (`FLOWVPM_relaxation.jl:41`). M2 is the magnitude- and
  overlap-correct generalization of that crude direction-only rotation. This ties
  back to item 005's finding that relaxation does real, load-bearing work: M2
  explains *why* (it is fighting the irrotational `-\nabla\varphi`) and proposes
  doing it consistently rather than as an after-the-fact damping.

**Viscous compatibility.** `\Gamma_{\mathrm{eff}}` is still a Gaussian-basis
coefficient, so the clean core-spreading law `\sigma^2\mathrel{+}=2\nu\,dt`
continues to apply to the basis. The projection `M^{-1}H` is linear; note and
check any concern about whether projection and core-spreading commute over a step,
but there is no obvious obstruction to retaining the Gaussian viscous model.

### M4 (explored): the Leray-projected Gaussian blob

**M4 in four claims.**
1. M4 means carrying the **Leray-projected Gaussian blob**
   `K_\sigma=P(\zeta_\sigma I)=\zeta_\sigma I+\nabla\nabla G_\sigma` (sign per the
   `-\Delta G_\sigma=\zeta_\sigma` convention) as the represented vorticity — which
   equals `\omega_J` — **not** a local compact divergence-free RBF.
2. Gaussian viscous core spreading **survives**, because the projection commutes
   with heat diffusion.
3. `K_\sigma` is exactly solenoidal and velocity-consistent, but its algebraic
   Hessian tail makes whole-space moments more subtle than the scalar blob's.
4. This fixes representation/velocity consistency, but does **not** by itself solve
   the time-evolution problem — deriving the projected-basis strength update is the
   real remaining work.

**Key recognition — the correct divergence-free basis is the Leray projection of
the existing blob, and it is already implicit in the model.** Apply the Leray
(Helmholtz) projector `P` to a scalar-Gaussian blob:

```math
\omega_i(x)=P\bigl(\Gamma_i\,\zeta_{\sigma_i}\bigr)
=\bigl(\zeta_{\sigma_i}I+\nabla\nabla G_{\sigma_i}\bigr)\Gamma_i
=K_{\sigma_i}(x-X_i)\,\Gamma_i,
\qquad
K_\sigma=\zeta_\sigma I+\mathrm{Hess}(G_\sigma),
```

using `-\Delta G_\sigma=\zeta_\sigma`. Summed over particles this is **exactly
`\omega_J`, with the same coefficients `\Gamma_i`** — precisely the vorticity
whose velocity FLOWVPM already computes. Two settled properties make this the
*physically compatible* divergence-free Gaussian basis:

- **Exact solenoidality.** Because `K_\sigma` is the Leray projection of the scalar
  blob, `\nabla\cdot(K_\sigma\Gamma)\equiv0`; equivalently, it is the curl of the
  same regularized velocity field whose gradient FLOWVPM already stores.
- **Velocity is self-consistent.** Under the usual whole-space, decaying,
  no-harmonic-field assumptions — the relevant setting for free vortex particles —
  Biot–Savart of `\omega_J` returns the same solenoidal `u` (uniqueness of the
  decaying solenoidal field with a given curl), so adopting `K_\sigma` does **not**
  change the velocity kernel; it relabels the vorticity to be consistent with the
  velocity FLOWVPM already produces.

Relation to the rest of the menu is unchanged: `\omega_\zeta` uses only the
**isotropic part** `\zeta_\sigma I`; adding the **anisotropic**
`\mathrm{Hess}(G_\sigma)` recovers the full `K_\sigma` and makes the represented
vorticity **exactly solenoidal**. **M1** is its diagnostic face; **M2** best-fits
it with the isotropic part only (and — because a scalar Gaussian blob is never
divergence-free — can only ever be a best fit).

**What you give up is locality, not the Gaussian.** `K_\sigma` is built from the
*Green's function* `G_\sigma` (a `1/r` mollified potential), so
`\mathrm{Hess}(G_\sigma)` carries an algebraic `\sim 1/r^3` tail. The represented
vorticity keeps a Gaussian **core** but is no longer compactly Gaussian-decaying;
it has a nonlocal, sign-changing tail whose volume moments require care. The
Biot–Savart machinery is unchanged (the velocity already has its own far tail).

**Whole-space moment caveat — do not overclaim scalar-blob circulation
preservation.** The tempting argument `\int\mathrm{Hess}(G_\sigma)\,dx=0` is false
for this Green-function Hessian because `G_\sigma\sim1/(4\pi r)` and
`\nabla\nabla G_\sigma\sim r^{-3}` has a non-vanishing surface contribution at
infinity. A symmetric finite-volume integral of the pointwise kernel tends toward
the angular average of the Leray projector, `\int K_\sigma\Gamma\,dx\to
\tfrac23\Gamma`, not `\Gamma`. Thus M4 preserves the velocity-relevant projected
vorticity and the carried coefficient still labels the particle's source
strength in the vector potential, but the scalar blob's naive total-vorticity
moment does **not** carry over unchanged. This does not affect the local identity
`K_\sigma\Gamma=\omega_J` or viscous diffusion, but it is a real representation
caveat to keep explicit before any live model change.

**The viscous model survives — cleanly, via Fourier (the decisive property,
cleaner than M2).** In Fourier space,

```math
\hat\omega_i(k)=P(k)\,\Gamma_i\,e^{-\sigma_i^2|k|^2/2},
\qquad
P(k)=I-\frac{k k^\top}{|k|^2}.
```

Viscous diffusion multiplies by `e^{-\nu|k|^2 t}`, and `P(k)` is `t`-independent,
so it **commutes** with diffusion:

```math
\hat\omega_i(k,t)=P(k)\,\Gamma_i\,e^{-(\sigma_i^2+2\nu t)|k|^2/2}
\;\Longrightarrow\;
\sigma_i^2\leftarrow\sigma_i^2+2\nu\,dt
```

(the usual Gaussian core-spreading law, up to the convention-dependent factor).
The projection does not obstruct diffusion — what is preserved is the **Gaussian
heat kernel / Gaussian core spreading**; the result is a *projected Gaussian
basis*, not a local Gaussian one.

**Alternative flavor — a purely local divergence-free Gaussian RBF, with a
serious caveat.** One could instead use the Narcowich–Ward / Lowitzsch matrix
kernel `(\nabla\nabla^\top-\Delta I)\phi` built on a *local* Gaussian `\phi` (not
the Green's function) — compactly decaying and exactly divergence-free. But each
such basis function has **zero volume integral**
(`\int(\nabla\nabla^\top-\Delta I)\phi\,dx=0`), hence zero net circulation, so
`\Gamma_i` would lose its meaning as a vortex-blob circulation vector. Clean
mathematically, but not obviously the right particle model — flagged, not the
default.

**Remaining open task for M4 (projected flavor):**
- **Strength evolution.** Re-pose `D\omega/Dt=(\omega\cdot\nabla)u` for the
  projected basis: how does the stretching term project onto `K_\sigma`
  coefficients, and does `\dot\Gamma_i` stay as local as the present scalar
  update? (Velocity self-consistency and viscous diffusion are settled above; this
  is the main remaining unknown.) **Cost:** `K_\sigma` adds the
  `\mathrm{Hess}(G_\sigma)` sum — the same analytic Gaussian Hessian FLOWVPM's FMM
  already evaluates for velocity gradients (`FLOWVPM_fmm.jl:114`) — plus the
  longer-range `1/r^3` tail; quantify the overhead.

**M2 vs M4.** M4-projected is the **exact** representation fix
(`\nabla\cdot\omega\equiv0`, `\omega=\nabla\times u`, velocity-consistent) and
shares M2's Hessian/FMM machinery; its costs are the nonlocal tail, the
whole-space moment caveat above, and a re-derived strength update. M2 is the
minimal-change, coefficient-only route that leaves a best-fit residual. Develop
them together.

### M6 (explored): projected-Gaussian overlap evolution

**Definition.** Promote the M4 projected blob from a representation observation to
the actual particle model:

```math
\omega(x,t)=\sum_i K_{\sigma_i}(x-X_i)\Gamma_i,
\qquad
K_\sigma=P(\zeta_\sigma I)=\zeta_\sigma I+\nabla\nabla G_\sigma.
```

Then derive the strength equation in this **same matrix-valued basis**, rather
than deriving an overlap solve in the scalar `\zeta` basis and later trying to
repair its non-solenoidal residual. In collocation form, the M6 analogue of item
008 is

```math
\sum_i K_{\sigma_i}(X_k-X_i)\,\dot\Gamma_i = b_k^{(K)},
```

where `b_k^{(K)}` is the projected-basis right-hand side from the vorticity
equation:

```math
\frac{D\omega}{Dt}=(\omega\cdot\nabla)u
```

plus the same model terms item 008 keeps explicit: convection of particle
centers through a nonuniform velocity field, core spreading through
`\dot\sigma_i\,\partial_{\sigma_i}K_{\sigma_i}`, SFS/baroclinic terms if retained,
and any chosen local truncation or regularization.

**Why this is better than the M1-M4 menu alone.**
- It eliminates the velocity/vorticity discrepancy exactly, because the carried
  vorticity is `\omega_J=\nabla\times u` by construction.
- It keeps `\Gamma_i` as the coefficient of the velocity-consistent projected
  basis, but must inherit M4's whole-space moment convention instead of assuming
  the scalar-blob identity `\int K_\sigma\,dx=I`.
- It preserves clean viscous core spreading, because `P(k)` commutes with the heat
  multiplier.
- It avoids M2's unavoidable best-fit residual: the left-hand side and the
  represented vorticity live in the same divergence-free basis.
- It turns M4's "remaining open task" into a concrete deliverable rather than a
  caveat.

**No-overlap sanity check.** At a particle center,
`K_\sigma(0)=\tfrac23\zeta_\sigma(0)I`. If only particle `i` contributes at
`X_i`, the collocation equation becomes

```math
K_{\sigma_i}(0)\dot\Gamma_i
=
J_i\,K_{\sigma_i}(0)\Gamma_i,
```

so the scalar factor cancels and the model recovers the standard isolated update
`\dot\Gamma_i=J_i\Gamma_i`. This is the key consistency check M6 must pass before
touching real rotor-hover states.

**Projected-basis RHS to derive carefully.** The scalar-basis item 008 equation
contains basis-motion terms

```math
\sum_i \Gamma_i
  \left(u(x)-u_i\right)\cdot\nabla\zeta_i
+\sum_i\Gamma_i\dot\sigma_i\partial_{\sigma_i}\zeta_i.
```

M6 needs the matrix-kernel analogue:

```math
\sum_i \left[\left(u(x)-u_i\right)\cdot\nabla K_i\right]\Gamma_i
+\sum_i \dot\sigma_i\,(\partial_{\sigma_i}K_i)\Gamma_i,
```

with `K_i=K_{\sigma_i}(x-X_i)`. The stretching term should use the represented
vorticity `\omega=\sum_iK_i\Gamma_i`, not `\omega_\zeta`. Sign, transpose, and
Jacobian conventions must be checked against FLOWVPM's `J_{ab}=\partial u_a /
\partial x_b` layout.

**Implementation shape.** M6 likely reuses the same ingredients already identified
for M2/M4:
- matrix-vector products with `K_\sigma=\zeta I+\mathrm{Hess}(G_\sigma)`;
- neighbor-list or FMM-accelerated Hessian sums;
- local sparse/CG solves analogous to the RBF core-reset machinery, but with a
  `3\times3` block kernel instead of scalar `M_{ki}I`;
- Pedrizzetti/corrected-Pedrizzetti and rVPM `Z`-coupling as baselines, not as
  parts of the definition.

**Acceptance target for M6.** Produce either a defensible projected-basis
`\dot\Gamma` equation with no-overlap reduction, preserved core-spreading term,
and an offline residual diagnostic on saved particle states, or a clear negative
result showing that the block-kernel solve is too ill-conditioned, too nonlocal,
or redundant with existing FLOWVPM machinery. Report cost separately for
self-only, neighbor-truncated, and FMM/global variants.

**Status (2026-06-25): carried out — defensible equation + a clear negative
diagnostic result.** The projected-basis `\dot\Gamma` equation is derived below
and its no-overlap reduction (`\dot\Gamma_k=J_k\Gamma_k`), block-solve structure,
viscous core-spreading cancellation, and `\partial_\sigma K` term are
self-tested at N=1 (`M6_SELFTEST=1
examples/particle_vorticity_curl_vs_basis_check.jl`). The offline saved-state
diagnostic (`examples/particle_m6_projected_basis_diag.jl`) returns the **negative
branch**: in the settled rotor-hover wake the block operator is *more* ill-conditioned
than the scalar overlap `M` and a truncated local block solve produces meaningless
rates. See the dated section "2026-06-25 M6 — projected-basis strength evolution
derived + diagnostic".

### Enumerated alternatives (not explored in this item)

- **M3 — periodic Leray cleaning.** Project `\omega_\zeta\leftarrow
  \omega_\zeta-\nabla\varphi` with `\nabla^2\varphi=\nabla\cdot\omega_\zeta`; for
  Gaussians `\varphi` is a closed-form Gaussian sum (kernel sum, not a grid
  solve). Offline/maintenance cadence (e.g., piggyback on core reset).
- **M5 — existing formulation levers as baselines.** rVPM `Z`-coupling
  (`FLOWVPM_formulation.jl:31`, `FLOWVPM_timeintegration.jl:103`) and
  corrected-Pedrizzetti as already-shipped partial mitigations to benchmark any
  new `\Gamma_{\mathrm{eff}}` / `K_\sigma` scheme against.

### Additional staged ideas

M7-M11 are not acceptance scope for this item. They are parked here to avoid
losing viable directions while keeping 010 focused on M2/M4/M6. **M12 is
different:** it is a required theory-clarification stage in item 010's completion
path; it resolves the representation no-go/interpretation question while leaving
the projected-basis dynamics and cost to M6.

- **M7 — projected evolution with scalar storage.** Keep scalar-blob coefficients
  and existing maintenance machinery, but project every physics residual before
  fitting coefficient updates:

  ```math
  \omega=P(\omega_\zeta),
  \qquad
  \partial_t\omega
  =
  P\left[-(u\cdot\nabla)\omega_\zeta
  +(\omega_\zeta\cdot\nabla)u
  +\nu\Delta\omega_\zeta+\cdots\right].
  ```

  This is less exact than M6 because the stored basis remains scalar and
  non-solenoidal, but it may be a cheaper bridge that keeps all dynamic residuals
  on the solenoidal side before refitting.
- **M8 — carry the vector potential as the primary state.** Treat the object
  FLOWVPM already differentiates as the model state:

  ```math
  A(x)=\sum_i\Gamma_iG_{\sigma_i}(x-X_i),
  \qquad
  u=\nabla\times A,
  \qquad
  \omega=\nabla\times\nabla\times A.
  ```

  This is conceptually clean and preserves Gaussian core spreading, but it moves
  the hard derivation to an evolution law for `A` plus a gauge choice.
- **M9 — split solenoidal and gauge-error coefficient content.** Decompose the
  scalar representation into parts used for curl/velocity physics and parts that
  only fit scalar-basis residue, e.g. `\Gamma=\Gamma^S+\Gamma^G`, then measure,
  damp, discard, or constrain `\Gamma^G`. Mostly diagnostic, but it could explain
  when relaxation is suppressing nonphysical residue.
- **M10 — moment-corrected projected particles.** Extend M4/M6 with local moments
  built from derivatives of the projected Gaussian heat kernel:

  ```math
  \omega_i=K_\sigma\Gamma_i+\nabla K_\sigma:M_i+\cdots.
  ```

  Projection and diffusion remain clean if all moments are Gaussian-derivative
  based. This may reduce the burden on large overlap solves, at the cost of more
  state and stability questions.
- **M11 — local divergence-free cloud reconstruction.** During maintenance/core
  reset, reconstruct a divergence-free local polynomial/Taylor model of `\omega`
  on each particle neighborhood, then refit projected Gaussian particles to that
  model. This is closer to meshfree remeshing than pure VPM, so reserve it for a
  future maintenance-focused item.
#### M12 — no local circulation-monopole kernel; use velocity-primary pairing

M12's literal first hope was a new particle representation that is simultaneously
local/no-algebraic-tail, exactly divergence-free, Gaussian/heat-kernel compatible,
FMM-efficient, and still carries circulation in the same `\int\omega\,dV` sense as
the scalar blob. That target is impossible in the precise zeroth-moment sense.

**No-go theorem (generator-independent).** Let `K(x)` be a matrix vorticity kernel
whose columns are exactly divergence-free and local enough that `K\in L^1`, so
`\widehat K(k)` is continuous at `k=0`. Columnwise solenoidality gives

```math
k^\top\widehat K(k)=0\qquad(k\ne0).
```

Taking `k=t e_a\to0` along each coordinate axis gives
`\widehat K_{ab}(t e_a)=0`, hence by continuity `\widehat K_{ab}(0)=0` for every
row `a` and column `b`. Therefore

```math
\boxed{\;\widehat K(0)=\int_{\mathbb R^3}K(x)\,dx=0\;}
```

for every local exactly divergence-free matrix kernel. The zero integral of the
Narcowich-Ward / Lowitzsch-style local divergence-free Gaussian RBF is therefore
not an accident of that generator; it is the required consequence of locality plus
exact solenoidality. The impossible triad is:

```text
local / L1 tail  +  exact divergence-free  +  nonzero zeroth vorticity moment
```

**Trilemma.** A particle model must relax one corner:

- **Relax locality:** M4's Leray-projected blob `K_\sigma=\zeta_\sigma I+
  \mathrm{Hess}(G_\sigma)` is exact-solenoidal and velocity-consistent, but has the
  nonlocal `1/r^3` Hessian tail. This tail is FMM-benign because the same analytic
  Hessian structure is already used for velocity-gradient evaluation.
- **Relax exact solenoidality:** the scalar Gaussian blob keeps the simple
  `\int\omega_\zeta\,dV=\Gamma` moment, but this is exactly the inconsistency item
  010 exposes.
- **Relax `\int\omega\,dV` as the circulation definition:** keep the velocity
  kernel primary and define represented vorticity as its curl,
  `\omega_J=\nabla\times u=K_\sigma\Gamma`. Then `\Gamma` keeps its physical
  circulation meaning through velocity circulation / vorticity flux, not through
  the full-volume vector moment.

**Recommended M12 resolution: velocity-primary paired representation.** State the
particle as `(X,\Gamma,\sigma)` with the existing Gaussian-erf Biot-Savart velocity
kernel as primary. The represented vorticity is defined by the actual curl of that
velocity:

```math
u_\Gamma=\text{existing regularized Biot-Savart kernel},\qquad
\omega_\Gamma=\nabla\times u_\Gamma=K_\sigma\Gamma.
```

This preserves the existing velocity field, FMM path, and Gaussian core-spreading
model; it also makes the represented vorticity exactly divergence-free. The key
semantic correction is that physical circulation is

```math
\oint_{\partial S}u\cdot dl=\int_S\omega\cdot n\,dS,
```

not the full-volume vector moment `\int\omega\,dV`. The symmetric M4 volume moment
`\int\omega_J\,dV\to\tfrac23\Gamma` is a `k=0`/far-field convention diagnostic, not
loss of Stokes circulation. Thus `\Gamma` remains a circulation/velocity-strength
coefficient even though it is not the scalar blob's naive vorticity monopole.

This also clarifies what a local exactly divergence-free zero-mean particle can
represent. Such a kernel cannot converge as a vorticity distribution to an isolated
open filament segment with nonzero `\int\omega\,dV`. It can naturally represent
closed vortex structures, where the leading information is in impulse/dipole-like
moments (for a small ring, roughly `\rho\Gamma A`) rather than the zeroth moment.
Angular momentum or impulse comparisons should therefore be made through the
converged velocity field or the appropriate vorticity moments, not through
`\int\omega\,dV` alone.

M12's theory result is therefore a clear negative result for the literal local
circulation-monopole target and a viable representation path: velocity-primary M4
pairing. The only open piece is the strength dynamics, which remains M6's
projected-basis `\dot\Gamma` derivation and cost analysis.

## Acceptance

This item is **not complete** until it delivers **both**:

(a) **the discrepancy thoroughly exposed in math and interpreting prose** — the
exact identity `\omega_J=\omega_\zeta+\nabla(\nabla\cdot A)` with fixed signs for
the Gaussian-erf kernel; the gauge + Helmholtz framing (`\omega_J` is the
solenoidal projection of `\omega_\zeta`); the closed-form single-particle
self-term `\omega_J(X)=\tfrac23\zeta_\sigma(0)\Gamma`; the locality sharpening
(`\delta=\nabla(\nabla\cdot A)`, not pointwise `\nabla\cdot\omega_\zeta`); the
reconciliation with why VPM still works; and the N=1 numerical confirmation
(`\to 1/3` at center, `\delta` matching `\mathrm{Hess}(G_\sigma)\Gamma` to
round-off, `\nabla\cdot\omega_J` at round-off); **and**

(b) **a thorough, in-depth exploration of M1 (the divergence-correction term),
M2 (effective `\Gamma_{\mathrm{eff}}`),
M4 (the Leray-projected Gaussian blob
`K_\sigma=\zeta_\sigma I+\mathrm{Hess}(G_\sigma)=\omega_J`), and M6
(projected-Gaussian overlap evolution)** — for M1: its definition (carry `\Gamma`
unchanged, evaluate `\omega_J=\omega_\zeta+\sum_i\mathrm{Hess}(G_{\sigma_i})\Gamma_i`
on demand), the closed-form analytic Hessian
`\mathrm{Hess}(G_\sigma)\Gamma=A\Gamma+B(y\!\cdot\!\Gamma)y` (numerically matched to
the ForwardDiff term to round-off), the cost/FMM reuse, the diagnostic-vs-physics
role, the self-vs-overlap decomposition framing, and viscous compatibility; for M2:
its
definition, the honest best-fit/representation caveat, the direct
`\dot\Gamma_{\mathrm{eff}}` re-derivation, the cheap-limit (`\approx\tfrac23\Gamma`
self-term) cost comparison against Pedrizzetti, reuse of the RBF/relaxation
machinery, and Gaussian-viscous compatibility; for M4: the recognition that
`K_\sigma` is already implicit in the model, the exact-solenoidality, velocity
self-consistency, preserved-diffusion (Fourier: `P(k)` commutes with
`e^{-\nu|k|^2t}`), and whole-space moment caveat (the Green-function Hessian tail
prevents the naive `\int K_\sigma=I` claim) results, the locality tradeoff
(nonlocal `1/r^3` tail, not a local basis), the zero-moment caveat of the
local-RBF flavor; for M6: the projected-basis `\dot\Gamma` equation,
the block-kernel overlap solve, the no-overlap reduction to
`\dot\Gamma_i=J_i\Gamma_i`, the `\dot\sigma\,\partial_\sigma K` viscous
core-spreading term, and offline residual/cost diagnostics — with M3, M5
enumerated as alternatives only, plus the open self-vs-overlap decomposition
question framed as future work; and **M12** as the theory
clarification that proves the local/exact-solenoidal/nonzero-`\int\omega`
triad impossible, distinguishes Stokes circulation from the full-volume
vorticity moment, and recommends the velocity-primary paired representation
(`\omega=\nabla\times u`) while leaving strength dynamics to M6.

## Caveats

- Do not use rotor-hover saved states as the first evidence; the N=1 synthetic
  test isolates the self-term cleanly and the source of the discrepancy is
  controlled.
- Do not treat `\omega_J` (curl-of-J) as a normalization calibration of `\zeta`.
  By the identity above the two genuinely differ by `\nabla(\nabla\cdot A)`; the
  008 `0.186` is that physical discrepancy, not kernel-constant slop.
- The split between the isolated self-term and genuine inter-particle overlap in
  a real wake is **unquantified**; treat any single-number wake mismatch as a mix
  until the decomposition experiment is run.
- Keep the M2 best-fit caveat explicit: an effective *scalar*-Gaussian strength
  reduces but cannot eliminate the basis divergence error — only the matrix-valued
  projected-Gaussian basis (M4) makes `\nabla\cdot\omega` identically zero.
- For M4, distinguish the two flavors: the **Leray-projected blob**
  (`K_\sigma=\zeta_\sigma I+\mathrm{Hess}(G_\sigma)`) is exact-solenoidal and
  velocity-consistent but has a nonlocal `1/r^3` tail and a nontrivial
  whole-space moment convention; a **local** divergence-free Gaussian RBF is
  compact but has zero net circulation, so it does not obviously preserve the
  physical meaning of `\Gamma_i`. Default to the projected blob.
- For M6, do not silently reuse item 008's scalar `M\dot\Gamma=b` equation. The
  point is to move the whole strength-evolution derivation into the projected
  matrix basis `K_\sigma`; otherwise M6 collapses back to M2/M3 and loses exact
  solenoidality.
- M7-M11 are only parking-lot ideas here. Do not require or evaluate them for
  item 010 technical completion unless a future item explicitly promotes one.
- M12 is required and now resolved at the theory level: the literal
  local/no-algebraic-tail circulation-monopole target is impossible, and the
  viable path is the velocity-primary paired representation. Do not treat this as
  a live dynamics solution; M6 still owns the `\dot\Gamma` derivation.

## 2026-06-24 N=1 Numerical Verification

Implemented `examples/particle_vorticity_curl_vs_basis_check.jl`, a
self-contained ForwardDiff check using the same AD-clean Gaussian-erf velocity
kernel structure as item 009. The script uses one particle with `σ=1`,
`Γ=(1,0.2,0)`, evaluates the center plus a short radial line, and compares:

- `ω_ζ=Γ ζ_σ(x-X)`;
- `ω_J=∇×u`, from ForwardDiff through the regularized Biot-Savart kernel;
- `δ=ω_J-ω_ζ`;
- `δ_analytic=Hess(G_σ)Γ`, with an AD-clean Taylor branch for `G_σ` at the
  particle center; and
- `∇·ω_J=tr(∇ω_J)`, from a nested ForwardDiff Jacobian.

Run:

```bash
julia --project=test examples/particle_vorticity_curl_vs_basis_check.jl
```

Result:

| r/σ | \|\|δ\|\|/\|\|ωζ\|\| | δ vs Hess(G)Γ relerr | rel div ωJ | 2/3 center err |
| ---: | ---: | ---: | ---: | ---: |
| 0.000 | 3.333333e-01 | 1.797172e-16 | 0.000000e+00 | 0.000000e+00 |
| 0.125 | 3.332564e-01 | 6.114427e-14 | 0.000000e+00 | 1.504558e-03 |
| 0.250 | 3.330649e-01 | 1.215124e-14 | 1.912471e-17 | 6.058795e-03 |
| 0.500 | 3.329341e-01 | 6.401928e-15 | 0.000000e+00 | 2.490143e-02 |
| 1.000 | 3.457911e-01 | 1.792448e-15 | 0.000000e+00 | 1.114773e-01 |
| 2.000 | 8.058638e-01 | 1.890840e-15 | 0.000000e+00 | 7.517007e-01 |
| 4.000 | 8.353614e+01 | 6.305829e-13 | 1.259344e-16 | 8.355327e+01 |

The decisive center gate passes:

```text
center discrepancy: 3.3333333333333337e-01 (target 1/3 = 3.3333333333333331e-01)
max delta analytic relative error: 6.306e-13
max relative div(omega_J): 1.259e-16
RESULT: PASS -- curl-of-J is the projected vorticity, with the 2/3 self-term at the center.
```

Conclusion: the N=1 synthetic test confirms the theory. The discrepancy between
`ω_J` and `ω_ζ` is not a normalization/sign artifact; it is the analytic
projection term `Hess(G_σ)Γ`. At the particle center, `ω_J=(2/3)ω_ζ` exactly to
roundoff, while `∇·ω_J` remains roundoff-zero. Report written to
`data/particle_vorticity_curl_vs_basis/curl_vs_basis_report.csv`.

## 2026-06-24 M1 (divergence-correction term) — explored

Promoted **M1** from an enumerated alternative to a full "(explored)" subsection
(now placed before M2 as the diagnostic face / simplest member of the `K_σ`
family). Headline deliverable: a **closed-form analytic `Hess(G_σ)` kernel** that
needs no ForwardDiff and no grid,

```
Hess(G_σ)(y)Γ = A·Γ + B·(y·Γ)·y,
A = −φ(s)/(4π σ³),   B = −φ′(s)/(2π σ⁵),   s = |y|²/σ², y = x−X,
φ′(s) = [√(2/π) ρ³ e^{−ρ²/2} − 3 g(ρ)] / (2 ρ⁵)   (Taylor branch near s=0),
```

with `φ`, `g` the existing Gaussian-erf profiles. Center reduces to
`A(0)Γ = −⅓ζ_σ(0)Γ` (the `2/3` self-term) and `tr Hess(G_σ) = −ζ_σ`, fixing all
signs/constants.

Implemented in `examples/particle_vorticity_curl_vs_basis_check.jl` as
`phi_prime_gauserf`, `hess_green_times_gamma`, and the M1 diagnostic
`omega_J_corrected` (ω_ζ + Σ Hess(G_σᵢ)Γᵢ over a particle set; brute-force here,
production path swaps in the neighbor-list/FMM Hessian sum from
`examples/particle_overlap_residual.jl`). An opt-in `M1_SELFTEST=1` path asserts the
closed form matches the ForwardDiff `delta_analytic`:

```text
M1 self-test max relerr = 7.687e-13 -- PASS
```

The default `main()` gates are unchanged and still PASS. **Deferred** (per scope):
the numerical self-vs-overlap decomposition on a saved rotor-hover state — split
`Σᵢ Hess(G_σᵢ)Γᵢ` at each particle into the closed-form self-term
`−⅓ζ_σ(0)Γ_k` and the neighbor remainder, quantifying why the item-008 `0.186`
sits below the isolated `1/3`. The `.vtp` loader and neighbor lists already exist in
`examples/particle_overlap_residual.jl`.

## 2026-06-25 M2 (effective strength Γ_eff) — derivation carried out + projection implemented

Filled the two gaps in the M2 "(explored)" prose: the `Γ̇_eff` evolution law was only
*asserted*, and there was no runnable projection.

**Derivation carried out.** Representing the solenoidal `ω_J` in the scalar ζ basis
with coefficients `Γ_eff(t)` and projecting the inviscid transport
`Dω/Dt=(ω·∇)u` onto the basis at the convecting centers gives

```
M Γ̇_eff = b_eff,
b_eff,k = J(X_k)(M Γ_eff)_k − Σ_i Γ_eff,i (u(X_k)−u_i)·∇ζ_i(X_k),
```

(no core-spreading term: for the Gaussian heat kernel `σ̇ ∂_σζ = ν∇²ζ`, so viscous
diffusion is the separate core-spreading σ-update with `Γ` held fixed and cancels the
`ν∇²ω` source — it does not force `Γ̇_eff`.)

— item-008 basis-motion structure with `ω_J=MΓ_eff` as the stretched field (classic;
transposed uses `Jᵀ`). Self-dominant limit: the overlap matrix cancels and
`Γ̇_eff,k = J_k Γ_eff,k`, i.e. **unchanged stretching on `Γ_eff`**; with the self-only
`Γ_eff≈⅔Γ` projection this is an `O(1)`/particle scheme cheaper than Pedrizzetti.
Honest subtlety (expanded in the M2 section): the projection `Γ_eff=(I+T)Γ`, `T=M⁻¹H`,
has a geometry-dependent `T(t)` (positions + core sizes), so `Ṫ≠0`. The consistent
`Γ̇_eff` carries a drift term `Ṫ Γ` (plus a `[S,T]` conjugation piece) that ordinary
stretching drops, so cheap stretching makes `Γ_eff` drift off the best fit. Realize M2
either by carrying that drift correction every step or by re-projecting at the RBF
core-reset cadence. Cost (overlap-dependent, expanded in the M2 section): in item 008's
heavy-overlap settled wake (`n̄≈2.25e4`, diagonal-dominance `≈0.006`) the drift term
inherits the per-step `M⁻¹` CG solve (`~15` matvecs, conjugation roughly doubling it), so
option (1) ≈ a full M6/008 solve — *not* cheap, and the local-`M⁻¹` shortcut fails since
`M` is far from diagonally dominant. Option (2) amortizes that solve over `K` steps and is
the cost-sensible realization; self-only `⅔Γ` is the free tier but *inaccurate* in heavy
overlap (the 2026-06-25 saved-state section below shows it is worse than raw, neighbor
Hessian ≈3.1× the self-term), so not a wake-ready standalone scheme.

**Projection implemented** in `examples/particle_vorticity_curl_vs_basis_check.jl`
(`overlap_matrix`, `omega_J_nodes`, `gamma_eff` dense, `gamma_eff_cg`, `zeta_field`),
reusing the M1 `hess_green_times_gamma` for `H`. Opt-in `M2_SELFTEST=1`:

```text
  (1) single particle: ||Geff - (2/3)Gamma||/||Gamma|| = 1.122e-16
  (2) two particles: CG vs dense = 2.737e-16, linear residual = 7.938e-17
  (3) off-node field vs omega_J: raw Gamma = 5.273e-01, Gamma_eff = 6.131e-02 (floor)
  (4) off-node rel divergence: raw Gamma = 6.424e-01, Gamma_eff = 3.596e-01 (floor)
M2 self-test -- PASS
```

So `Γ_eff` is the `2/3` self-limit exactly, the CG (mirroring `rbf_conjugategradient`)
matches the dense solve, the `Γ_eff` field is ~8.6× closer to `ω_J` off-node, and its
divergence drops ~1.8× — both with a nonzero floor (the scalar-basis best-fit limit;
closing it to zero is M4). The milder divergence gain vs the field-closeness gain is the
honest signature that scalar coefficients cannot make the field solenoidal. Default N=1
gates unchanged and still PASS. **Deferred** (per scope): `Γ_eff`/`Γ̇_eff` on a saved
rotor-hover state and any retargeting of the live FLOWVPM core-reset.

## 2026-06-25 M2 clear-context check

Checked the M2 derivation and implementation against the INDEX approval workflow,
without changing item-level gates. The derivation is internally consistent: the
material derivative at moving centers gives the item-008 basis-motion terms with
`Γ_eff` replacing `Γ`, the stretching term uses the represented vorticity
`ω_J=MΓ_eff`, and the self-dominant limit correctly cancels the scalar overlap
factor to recover `Γ̇_eff,k=J_kΓ_eff,k`. The projection-preservation caveat is also
necessary and correctly stated; ordinary stretching alone does not maintain
`Γ_eff=(I+M⁻¹H)Γ` as `M`, `H`, and particle positions evolve.

**Correction (later review, 2026-06-25):** this check missed one error in the first
write-up — `b_eff` originally carried a standalone core-spreading term
`−Σ_i Γ_eff,i σ̇_i ∂_σζ_i`. That term is wrong in the strength equation: the derivation is
inviscid (`σ̇=0`), and viscously the Gaussian heat-kernel identity `σ̇ ∂_σζ = ν∇²ζ` makes
it the modeled `ν∇²ω`, realized by the separate core-spreading σ-update with `Γ` fixed, so
it cancels rather than forcing `Γ̇_eff` (a pure-diffusion check `J=0` ⇒ `Γ̇=0` confirms;
the standalone term would give `+3(σ̇/σ)Γ`). It has been removed from the boxed law and
both restatements; the no-overlap reduction `Γ̇_eff,k=J_kΓ_eff,k` is now exact (it had
silently dropped that non-vanishing term). The projection code is unaffected (it only does
the static projection).

Re-ran:

```bash
julia --project=test examples/particle_vorticity_curl_vs_basis_check.jl
M1_SELFTEST=1 M2_SELFTEST=1 julia --project=test examples/particle_vorticity_curl_vs_basis_check.jl
```

Both passed. The M2 self-test reproduced the documented values:

```text
(1) single particle: ||Geff - (2/3)Gamma||/||Gamma|| = 1.122e-16
(2) two particles: CG vs dense = 2.737e-16, linear residual = 7.938e-17
(3) off-node field vs omega_J: raw Gamma = 5.273e-01, Gamma_eff = 6.131e-02 (floor)
M2 self-test -- PASS
```

No M2 blocker found. The remaining M2 limitations are already explicit and should
stay deferred rather than folded into this check: saved-state quantification of
the neighbor part of `M⁻¹H`, and any live FLOWVPM core-reset or strength-evolution
wiring. The CG mirror is correctly restricted to equal-`σ`/SPD cases; unequal core
sizes still require the dense/nonsymmetric interpretation or a separate solver
choice.

## 2026-06-25 M2 saved-state diagnostic — self-only is not enough

Implemented `examples/particle_m2_effective_strength_diag.jl`, an offline/read-only
diagnostic over saved particle VTP states. It reuses the item-008 VTP loader and
CSR neighbor grid, evaluates the analytic M1 Hessian kernel on sampled target
particles, and decomposes

```math
H\Gamma(X_k)=\sum_i \mathrm{Hess}(G_{\sigma_i})(X_k-X_i)\Gamma_i
=H_{\mathrm{self},k}+H_{\mathrm{neighbor},k},
\qquad
H_{\mathrm{self},k}=-\tfrac13\zeta_{\sigma_k}(0)\Gamma_k.
```

The shared Gaussian-erf Hessian/Taylor helpers now live in
`examples/particle_gauserf_hessian_helpers.jl`, included by both the synthetic
curl-vs-basis check and the saved-state M2 diagnostic, so future kernel fixes do not
diverge between scripts.

The default saved-state mode now samples 500 target particles per state using a
fixed `RNG_SEED` (default `1`) so sample selection is reproducible without aliasing
against particle ordering. Exhaustive all-target Hessian evaluation is much slower
than item 008's Gaussian-only pass because each neighbor pair evaluates the
erf-based Hessian. The script keeps exact synthetic gates and adds the real-state
closure metric

```math
\frac{\|(M\Gamma+H\Gamma)-\mathrm{curl}(J)\|}{\|\mathrm{curl}(J)\|},
```

which measures how well the truncated analytic target reconstructs the saved
`velocity_gradient` curl. This closure error includes Hessian truncation, saved-`J`
/ FMM error, and any convention mismatch; it is not a pure truncation error.

```text
[gate 1] isolated particle gives H_neighbor=0 and self residual ~0
[gate 2] two-particle local solve improves over raw/self/Jacobi
```

Verification / settled rotor-hover run:

```bash
LOCAL_SOLVE_NSAMPLE=20 LOCAL_SOLVE_MAX_NEIGHBORS=160 RNG_SEED=1 \
  julia --project examples/particle_m2_effective_strength_diag.jl \
  rotor_hover_pressure_comparison 340:359
```

Summary over 20 settled states (`data/rotor_hover_pressure_comparison/particle_m2_effective_strength_diag.csv`):

```text
basis omega vs curl-of-J rel-err : mean 2.408e-01
target (omega+H) vs curl-of-J    : mean 1.139e-01
||H Gamma|| / ||M Gamma||        : mean 1.686e-01
||H_neighbor|| / ||H_self||      : mean 3.097e+00
cos(H_total, H_self)             : mean 0.340
residual raw / self              : 1.854e-01 / 2.955e-01
local row residual raw/self/Jac/check: 4.215e-01 / 1.916e-01 / 2.694e+01 / 3.062e-08
local probe residual raw/self/Jac/solve: 4.182e-01 / 1.891e-01 / 2.670e+01 / 5.397e-01
diag-dominance median            : mean 5.620e-03
neighbors mean / max             : 22442.8 / 34927
```

Cutoff sensitivity on step 359 with `RNG_SEED=1`, `TARGET_NSAMPLE=150`, and
`LOCAL_SOLVE_NSAMPLE=0`:

| cutoff factor | target vs curl(J) | H/ω | H_neighbor/H_self | neighbors mean |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 2.661e-01 | 2.511e-01 | 2.457 | 4,931 |
| 4 | 2.123e-01 | 1.826e-01 | 1.682 | 21,999 |
| 8 | 2.143e-01 | 1.765e-01 | 1.633 | 39,169 |
| 12 | 2.144e-01 | 1.764e-01 | 1.633 | 39,394 |

Interpretation: adding the analytic `H\Gamma` target roughly halves the saved-state
`M\Gamma` vs `curl(J)` gap in the settled-state run, but it does **not** collapse it.
The cutoff sweep shows the closure error plateaus once the neighbor set is effectively
the whole particle cloud, so the remaining `~0.21` single-state closure floor is not a
simple `4σ` truncation artifact. It is instead a real-state closure mismatch involving
saved `J`/FMM accuracy, possible panel/body or self-term convention content in the saved
gradient, and/or analytic-vs-FLOWVPM Hessian convention differences. Therefore the
saved-state M2 numbers should be treated as **sampled diagnostic evidence**, not as a
fully closed proof of `ω_J=MΓ+HΓ` on the saved VTP states.

The random-sampled settled run still supports the practical M2 conclusion: the wake's
`ω_J-ω_ζ` gap is not mostly the isolated `2/3` self-term. The neighbor Hessian is about
`3.1×` the self contribution in norm and only weakly aligned with it (`cos≈0.34`). A raw
scalar-basis field is closer to the sampled target than the uniform self-only `2/3Γ`
field (`0.185` vs `0.296` residual), so the cheap self-only rescale is **not** a useful
standalone M2 model for this rotor-hover state.

The local 160-neighbor solve is now explicitly a row-solve check, not model-quality
evidence: it fits its own collocation rows to roundoff, but held-out local probes have
`~0.54` residual. The meaningful negative result remains that truncated Jacobi is
unstable (`~27×` row/probe residual) and should not be promoted as a live shortcut. If
M2 is pursued, use a proper projection at core-reset cadence or move to M4/M6; do not
wire a self-only `2/3Γ` or naive Jacobi correction into live FLOWVPM.

## 2026-06-25 M4 projected-Gaussian representation check

Worked M4 as a docs/theory item, not live FLOWVPM wiring. The default M4 object is the
Leray-projected Gaussian blob

```math
K_\sigma=P(\zeta_\sigma I)=\zeta_\sigma I+\mathrm{Hess}(G_\sigma),
```

so the represented vorticity is `\omega(x)=\sum_iK_{\sigma_i}(x-X_i)\Gamma_i`.
This is the same field as `\omega_J=\nabla\times u`; M4 therefore fixes the
representation/velocity mismatch by changing the represented basis, while leaving the
existing velocity kernel conceptually unchanged. It does **not** derive a live
strength-evolution update; that remains M6.

Updated the M4 prose to make four points explicit:

- `K_\sigma\Gamma` is exactly solenoidal and equals the curl of the existing
  Gaussian-erf velocity field.
- Gaussian core spreading remains compatible because the Fourier Leray projector
  `P(k)=I-kk^\top/|k|^2` commutes with the heat multiplier.
- M4 gives up locality, not the Gaussian: the core remains Gaussian, but
  `\mathrm{Hess}(G_\sigma)` contributes a nonlocal `1/r^3` tail.
- The earlier `\int K_\sigma\,dx=I` / circulation-preservation claim was too strong.
  For the Green-function Hessian, the algebraic tail has a surface contribution at
  infinity; a symmetric finite-volume integral of the pointwise kernel tends toward
  `\tfrac23\Gamma`, not `\Gamma`. Treat the carried `\Gamma` as the vector-potential
  source coefficient unless a future live model supplies a more careful moment
  convention.

Added opt-in `M4_SELFTEST=1` to
`examples/particle_vorticity_curl_vs_basis_check.jl`. It verifies:

```text
  (1) K_sigma center factor relerr = 1.657e-16
  (2) K_sigma Gamma vs curl(u) max relerr = 1.828e-14
  (3) relative div(K_sigma Gamma) max = 7.578e-14
  (4) finite-cube integral vs (2/3)Gamma relerr = 2.810e-09
M4 self-test -- PASS
```

Command:

```bash
M4_SELFTEST=1 julia --project=test examples/particle_vorticity_curl_vs_basis_check.jl
```

Default N=1 gates still pass unchanged. The M4 outcome is therefore stronger than
"prose only" for representation consistency, but it is deliberately still not a
production implementation.

## 2026-06-25 M12 — no-go theorem + velocity-primary pairing

Resolved M12 at the theory level. The original target — a local/no-algebraic-tail,
exactly divergence-free heat-kernel vorticity particle that still carries a nonzero
scalar-blob-style `\int\omega\,dV` moment — is impossible. The proof is independent
of the Gaussian generator: if a matrix kernel `K\in L^1` is exactly divergence-free
columnwise, then `k^\top\widehat K(k)=0` for `k\ne0`; taking `k\to0` along each
coordinate axis and using continuity forces `\widehat K(0)=\int K\,dx=0`.

This sharpens the representation tradeoff:

- M4 keeps exact solenoidality and velocity consistency by relaxing locality
  (the `1/r^3` projected tail).
- The scalar Gaussian keeps a nonzero volume moment by relaxing exact
  solenoidality.
- A local exactly divergence-free kernel must relax the idea that `\int\omega\,dV`
  is the circulation carried by one particle.

The recommended path is the velocity-primary paired representation: keep the
existing Gaussian-erf Biot-Savart velocity kernel as primary and define the
represented vorticity as `\omega=\nabla\times u=K_\sigma\Gamma`. Then `\Gamma`
keeps its physical meaning through velocity circulation / vorticity flux,

```math
\oint_{\partial S}u\cdot dl=\int_S\omega\cdot n\,dS,
```

not through the full-volume vector moment `\int\omega\,dV`. The M4
`\tfrac23\Gamma` result is therefore a symmetric-volume moment caveat, not a loss
of Stokes circulation. Closed vortex structures can still be represented by local
zero-mean divergence-free kernels through impulse/dipole-like moments, but an
isolated open filament segment with nonzero `\int\omega\,dV` cannot be represented
as a vorticity distribution by finite sums of such zero-mean kernels.

Added opt-in `M12_SELFTEST=1` to
`examples/particle_vorticity_curl_vs_basis_check.jl`. It verifies the theorem's
numerical signatures and the circulation/moment distinction:

```text
  (1a) scalar zeta volume moment relerr to Gamma = 4.214e-09
  (1b) local div-free Gaussian volume moment / ||Gamma|| = 1.074e-07
  (1c) M4 symmetric volume moment relerr to (2/3)Gamma = 2.810e-09
  (2) finite-loop circulation vs Stokes flux relerr = 1.374e-05
M12 self-test -- PASS
```

Command:

```bash
M12_SELFTEST=1 julia --project=test examples/particle_vorticity_curl_vs_basis_check.jl
```

No live FLOWVPM wiring is implied. M12 closes the representation theorem and
interpretation; M6 still owns the projected-basis strength evolution.

### Pedrizzetti relaxation is an overlap correction (2026-06-25)

A corollary of the `K_\sigma` self-term, and a sharpening of why M4/M6 as a
*representation* change do not retire relaxation. Pedrizzetti relaxation realigns
each particle's carried strength `\Gamma_i` with the local curl of velocity.
**Source verified** (`FLOWVPM.jl` dev v4.1.0,
`src/FLOWVPM_relaxation.jl:41`, the active backend per `Manifest.toml`):
`relax_pedrizzetti` reads `\omega` as the **curl-of-J**
`(J_6-J_8,\,J_7-J_3,\,J_2-J_4)=\nabla\times u=\omega_J` — *not* the scalar basis
field `\omega_\zeta` — and applies the magnitude-preserving rotation

```math
\Gamma_i \leftarrow (1-f)\,\Gamma_i + f\,\|\Gamma_i\|\,\hat\omega_{J}(X_i),
\qquad \hat\omega_J=\omega_J/\|\omega_J\|
```

(corrected-Pedrizzetti adds a `/\sqrt{b^2}` renormalization that fixes Pedrizzetti's
strength-decay bug; it reads the same `\omega_J`). So the alignment target is the
curl field `\omega_J(X_i)`.

**The correction angle is created entirely by overlap (airtight geometry).**
Decompose the target into its self- and neighbor parts:

```math
\omega_J(X_i)=K_{\sigma_i}(0)\,\Gamma_i+\sum_{j\ne i}K_{\sigma_j}(X_i-X_j)\,\Gamma_j
=\underbrace{\tfrac23\zeta_{\sigma_i}(0)\,\Gamma_i}_{\parallel\,\Gamma_i}
+\ \text{neighbor sum}.
```

`K_\sigma(0)=\tfrac23\zeta_\sigma(0)I` is **isotropic** (scalar × identity), so the
self-term is exactly parallel to `\Gamma_i` and contributes **zero angle**. Each
`K_{\sigma_j}(X_i-X_j)`, `j\ne i`, is a full anisotropic `3\times3` block (the
Hessian tilts vectors), so only the **neighbor (overlap) sum** can rotate
`\omega_J(X_i)` away from `\Gamma_i`. Therefore the misalignment
`\angle(\Gamma_i,\omega_J(X_i))` that Pedrizzetti removes is **entirely
overlap-induced**: in the no-overlap limit `\omega_J(X_i)\parallel\Gamma_i` and
Pedrizzetti is idle. The neighbor sum's component **perpendicular** to `\Gamma_i`
is the direction error Pedrizzetti fixes; its **parallel** component is the
magnitude error M2/M6 additionally fix — the same direction-only-vs-overlap-correct
ladder noted in the M2 section.

**Source vs. medium (state honestly).** Overlap does not *cause* the drift; the
discrete stretching dynamics do (per-particle `\Gamma_i` and the collective field
`\omega_J(X_i)=(\nabla\times u)(X_i)` follow related-but-different discrete updates
and fail to preserve the continuous `\Gamma\parallel\omega` Helmholtz identity).
Overlap is the *medium* that expresses and sizes the resulting angle.

**Consequence for M4/M6.** Adopting the velocity-primary projected basis
`\omega=K_\sigma\Gamma` makes the *represented field* solenoidal **trivially** — a
curl is divergence-free regardless of how misaligned the `\Gamma_i` are — but it does
**not** rotate `\Gamma_i` toward `\omega_J(X_i)`. The dynamics still act on the
carried `\Gamma_i`, so keeping that coefficient tracking the local solenoidal
vorticity (Pedrizzetti's actual job) is **not** retired by the representation
change; it is precisely the overlap reconciliation (M2's `\Gamma_{\mathrm{eff}}`
or M6's block solve). Hence M6 *with the solve set aside* does not eliminate
relaxation: the solenoidal readout was never the operative defect — coefficient-vs-
local-field tracking is, and that is irreducibly an overlap operation. This is the
same conclusion the merging discussion reaches from the other side: reduce overlap
and Pedrizzetti's workload shrinks (less neighbor-induced tilt) for the identical
reason the overlap solve becomes better-conditioned. **Deferred numerical check:**
the per-particle angle `\angle(\Gamma_i,\omega_J(X_i))` is directly computable from
the M6 diagnostic data (`\Gamma_k` and `\omega_J(X_k)`); measuring its distribution
on settled states would quantify Pedrizzetti's actual workload and how much is
overlap.

## 2026-06-25 M6 — projected-basis strength evolution derived + diagnostic

Carried M6 out: a defensible projected-basis `\dot\Gamma` equation with its
no-overlap reduction, viscous core-spreading cancellation, and `\partial_\sigma K`
term (all N=1 self-tested), plus an offline saved-state diagnostic. The diagnostic
lands on M6's **negative branch**: in the settled rotor-hover wake the block-kernel
overlap operator is *more* ill-conditioned than the scalar overlap `M` and a
truncated local block solve produces meaningless rates. So M6 is exact in
representation but, like M2, blocked operationally by the same heavy-overlap regime.

### Derivation

Represent the (solenoidal) vorticity in the matrix basis
`\omega(x,t)=\sum_i K_{\sigma_i}(x-X_i)\Gamma_i`,
`K_\sigma=\zeta_\sigma I+\mathrm{Hess}(G_\sigma)` (=`\omega_J=\nabla\times u`).
Each blob convects with its own velocity `u_i` and spreads with `\dot\sigma_i`, so
`\partial_t K_i=-(u_i\cdot\nabla)K_i+\dot\sigma_i\,\partial_{\sigma_i}K_i`.
Substituting `\partial_t\omega+(u\cdot\nabla)\omega=(\omega\cdot\nabla)u+\nu\nabla^2\omega`
and collocating at the moving centers `X_k` gives the block-kernel analogue of
item 008:

```math
\sum_i K_{\sigma_i}(X_k-X_i)\,\dot\Gamma_i
= J(X_k)\,\omega(X_k)
- \sum_i\bigl[(u(X_k)-u_i)\cdot\nabla K_i(X_k)\bigr]\Gamma_i
- \sum_i\dot\sigma_i\,[\partial_{\sigma_i}K_i(X_k)]\Gamma_i
+ \nu\nabla^2\omega(X_k),
```

with `\omega=\sum_iK_i\Gamma_i` the represented vorticity (classic stretching
`J\omega`; the FLOWVPM transposed default uses `J^\top`). This is a `3n_p\times3n_p`
linear system with `3\times3` blocks `\mathcal K_{ki}=K_{\sigma_i}(X_k-X_i)` (each
block symmetric; the assembled operator non-symmetric through the source-`\sigma_i`
dependence), versus item 008's scalar `M_{ki}I`.

- **No-overlap reduction.** With only particle `k` at `X_k`: `K_{\sigma_k}(0)=\tfrac23\zeta_{\sigma_k}(0)I`,
  the convective term vanishes (`u(X_k)-u_k=0`), `\omega(X_k)=K_{\sigma_k}(0)\Gamma_k`,
  and the scalar `\tfrac23\zeta(0)` factor cancels both sides ⟹
  `\dot\Gamma_k=J_k\Gamma_k`, the standard isolated VPM stretch. This is M6's
  decisive consistency check, recovered to round-off (below).
- **Viscous core-spreading term.** Same argument as M2: for the Gaussian heat
  kernel the explicit `-\sum_i\dot\sigma_i\,\partial_\sigma K_i\,\Gamma_i` is exactly
  the `\nu\nabla^2\omega` source realized by the separate `\sigma^2\mathrel{+}=2\nu\,dt`
  update with `\Gamma` fixed, so it **cancels** rather than forcing `\dot\Gamma`
  (inviscidly `\dot\sigma=0`). The term is present and verified
  (`\partial_\sigma(K_\sigma\Gamma)`, AD vs FD relerr `2.2\times10^{-11}`), not used
  in the inviscid hover diagnostic.
- **Convective basis-motion term.** The matrix analogue
  `(u-u_i)\cdot\nabla K_i` needs the third-derivative tensor `\nabla\mathrm{Hess}(G_\sigma)`;
  it is **deferred** in the saved-state diagnostic (the "reduced" RHS `J\omega`,
  exactly as item 008 staged reduced-vs-sampled). It does not change the
  conditioning conclusion below.

### N=1 self-test (`M6_SELFTEST=1 examples/particle_vorticity_curl_vs_basis_check.jl`)

```text
  (1) K(0) vs (2/3)zeta(0) I relerr = 1.639e-16
  (1) no-overlap Gdot vs J Gamma relerr = 0.000e+00
  (2) block solve residual = 1.066e-16, cond(K_block) = 8.892e+00, cond(M_scalar) = 5.259e+00
  (2) well-separated Gdot vs J Gamma relerr: d=20s 9.167e-04, d=40s 1.146e-04 (ratio 8.00 ~ 8 for 1/r^3 tail)
  (3) d/dsigma(K Gamma) AD vs FD relerr = 2.234e-11 (term present; cancels viscously)
M6 self-test -- PASS
```

The well-separated residual decays as `(\sigma/d)^3` (ratio 8.00 across `d=20\sigma\to40\sigma`),
i.e. the block coupling vanishes only **algebraically** through M4's `1/r^3` Hessian
tail, never to machine zero like a compact kernel — the honest M6/M4 nonlocality,
which is FMM-benign (same Hessian the velocity-gradient FMM already forms).

### Saved-state diagnostic (`examples/particle_m6_projected_basis_diag.jl`)

Read-only over saved VTP states. Reuses item 008's loader / CSR neighbor grid and
the M1 analytic Hessian helpers. For sampled targets it forms the represented
`\omega=M\Gamma+H\Gamma`, the reduced RHS `b_k=J_k\omega(X_k)`, and a local
neighbor-truncated block solve, comparing three rate estimates at the center —
**self** `J_k\Gamma_k`, **block-Jacobi** `K(0)^{-1}b_k=\tfrac{3}{2\zeta(0)}J_k\omega(X_k)`,
and the local **block solve** — plus block/scalar conditioning and block
diagonal-dominance. Synthetic gates pass (isolated ⟹ block solve `=J\Gamma` to
round-off, `cond=1`; two-particle solve fits its rows).

```bash
TARGET_NSAMPLE=80 BLOCK_SOLVE_NSAMPLE=12 BLOCK_MAX_NEIGHBORS=160 RNG_SEED=1 \
  julia --project examples/particle_m6_projected_basis_diag.jl \
  rotor_hover_pressure_comparison 340:359
```

Summary over 20 settled states
(`data/rotor_hover_pressure_comparison/particle_m6_projected_basis_diag.csv`):

```text
basis omega vs curl-of-J          : mean 2.209e-01
represented (omega+H) vs curl-J   : mean 6.959e-02
||H||/||omega||                   : mean 1.837e-01
Gdot solve-vs-self / jac-vs-self  : 1.198e+14 / 8.965e+01
Gdot solve-vs-jac                 : 1.000e+00
block solve residual (solve check): 4.830e-04
block cond / scalar cond (median) : 9.159e+14 / 4.603e+13
block diag-dominance (median)     : 9.778e-03
represented omega vs isolated     : 8.881e-01
scalar diag-dominance median      : 6.833e-03
neighbors mean / max              : 21258.9 / 34086
```

Interpretation (the negative result, with cost tiers):

- **Adding the matrix basis does not help conditioning — it hurts.** The block
  operator's median condition number `\approx 9.2\times10^{14}` is `\sim20\times`
  *worse* than the scalar overlap `M` (`4.6\times10^{13}`); both are effectively
  singular in this heavy-overlap regime (block diagonal-dominance `\approx0.010` vs
  scalar `0.0068`). The `3\times3`-block coupling adds anisotropic off-diagonal
  content without restoring diagonal dominance. (The `\bar n\approx2.1\times10^4`
  neighbor count quoted here is *within* the `4\cdot\sigma_{\max}` cutoff, which
  spans `\approx1.7\,R`; per item 011's direct measurement the **true** overlap is
  `\sigma/h\approx4.2`, `\approx394` within `2\sigma_{\mathrm{local}}`, kernel-effective
  `\approx181` — heavy, but not the cutoff count. The conditioning verdict rests on
  the kernel-weighted diagonal-dominance, not the raw count, so it stands.)
- **The local block solve is numerically meaningless.** `\|\dot\Gamma_{\rm solve}-\dot\Gamma_{\rm self}\|/\|\dot\Gamma_{\rm self}\|\approx10^{14}`:
  the near-null directions of `\mathcal K` amplify the RHS into absurd rates, even
  though the truncated solve is internally consistent (linear residual `\sim5\times10^{-4}`,
  itself degraded from round-off by the conditioning). This is the block analogue of
  M2's "fits its rows but the model is not trustworthy", but far more severe.
- **Even the well-posed block-Jacobi rate is far from self.** `\dot\Gamma_{\rm jac}`
  differs from `\dot\Gamma_{\rm self}` by `\approx90\times`, because the represented
  `\omega(X_k)` is `\approx0.89` (relative) away from the isolated `K(0)\Gamma_k` —
  overlap dominates the vorticity itself. So overlap **would** materially change the
  rate (consistent with item 008's GO), but the operator that should resolve it is
  unusable.
- **Cost tiers.** *Self-only* (`\dot\Gamma_k=J_k\Gamma_k`): `O(1)`/particle, free,
  but ignores the `\approx0.89` overlap content. *Neighbor-truncated block solve*:
  a dense `3n_{\rm loc}\times3n_{\rm loc}` factorization — `\sim27\times` the scalar
  local solve at the same neighbor count and `9\times` the storage — and **unstable**
  here (cond `10^{15}`), so the truncation is not a usable shortcut (the same verdict
  M2 reached for truncated Jacobi). *FMM/global block solve*: would need a
  block-preconditioned Krylov solve over the `3\times3`-block `\mathcal K`; with median
  cond `\sim10^{15}` it is the same near-singular system globally and offers no
  conditioning advantage over the scalar M2 solve while costing `\sim9\times` per matvec.

**M6 conclusion.** The projected matrix basis `K_\sigma` is the correct *exact*
divergence-free representation (M4) and yields a clean, self-consistent strength
equation (no-overlap reduction + viscous cancellation verified). But operationally
it is **redundant-to-worse** than the scalar overlap solve in the settled hover wake:
the block-kernel collocation operator is even more ill-conditioned than the already
near-singular scalar `M`, and neither a truncated nor a global block solve is a
usable live update there. As with M1/M2, the binding constraint is the heavy-overlap
regime, not the basis choice — pointing back to resolution/overlap control (items
004/005/006) rather than a richer particle kernel. This satisfies M6's acceptance
("a clear negative result … too ill-conditioned … report cost for self-only,
neighbor-truncated, and FMM/global variants").
