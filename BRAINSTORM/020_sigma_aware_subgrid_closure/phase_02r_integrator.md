# 020 Phase 2R — Corrected frozen-gradient geometric integrator

**Status:** implementation and analytic verification complete (2026-08-12);
decisive rotor screen 13154223 terminal and rescored in the Phase-2 pack.

## Governing map

For the live $f=0$ rVPM, hold the velocity-gradient stretching operator $L$
fixed over one step and solve

$$
\dot q=Lq,
\qquad
q(0)=\Gamma_n.
$$

Let

$$
q_{n+1}=e^{L\Delta t}\Gamma_n,
\qquad
r=\frac{\lVert q_{n+1}\rVert}{\lVert\Gamma_n\rVert}.
$$

The homogeneous geometric rVPM update is

$$
\Gamma_{n+1}=q_{n+1}r^{-3g},
\qquad
\sigma_{n+1}=\sigma_n r^{-g}.
$$

For $g=1/5$ and aligned strain $L\Gamma=s\Gamma$, $r=e^{s\Delta t}$,
$Z=gs=s/5$, and therefore

$$
\Gamma_{n+1}=e^{2Z\Delta t}\Gamma_n,
\qquad
\sigma_{n+1}=e^{-Z\Delta t}\sigma_n.
$$

This retains the physical aligned stretching missing from the superseded
prototype. The existing SFS contribution is an additive modeled forcing and is
applied afterward as an explicit first-order Lie split.

## Coupled molecular diffusion

Define the step's constant-effective contraction rate

$$
Z_{\mathrm{eff}}=\frac{g\log r}{\Delta t}.
$$

After the homogeneous map, CoreSpreading adds the exact diffusion contribution
for $y=\sigma^2$ under constant $Z_{\mathrm{eff}}$:

$$
y_{n+1}=e^{-2Z_{\mathrm{eff}}\Delta t}y_n+
\frac{\nu}{Z_{\mathrm{eff}}}
\left(1-e^{-2Z_{\mathrm{eff}}\Delta t}\right),
$$

with the continuous limit $y_{n+1}=y_n+2\nu\Delta t$ at
$Z_{\mathrm{eff}}=0$. No SGS closure is present in Phase 2R.

## Verified properties and limits

- Exact aligned-strain amplification and contraction through
  $\Delta tZ\in\{0.5,1,2,5,50\}$.
- General constant-gradient agreement with an independent matrix-exponential
  reference for both FLOWVPM stretching-operator conventions.
- Positive $\sigma$ for every finite gradient/timestep represented without
  floating-point overflow.
- First-order agreement with the stock Euler step in the smooth limit.
- First-order convergence for a prescribed time-varying aligned gradient
  against its closed-form solution.
- Exact coupled CoreSpreading formula and its $Z\to0$ limit.
- Explicit rejection of reformulations with $f\ne0$; their geometry requires a
  separate derivation.
- SFS remains first-order split, not exactly integrated. Thus “exact” refers
  only to the homogeneous frozen-gradient geometry and constant-effective-$Z$
  molecular diffusion.

Verification gates:

- Full local FLOWVPM suite passes, including rings/leapfrog, merging 36/36,
  filament storage 477/477, calibration 40/40, SFS 505/505, vorticity storage
  66/66, relaxation 11/11, and the corrected integrator 28/28.
- FLOWPanel simulation/constructor coverage, including the opt-in integrator
  path, passes.
- Cluster Julia 1.11 corrected-integrator gate passes 28/28 after checksum
  identity was established for all staged FLOWPanel/FLOWVPM sources.
