# Diagnosis: SFS blow-up in `rotor_hover_pressure_comparison`

## TL;DR

The instability is **not** in the FLOWPanel coupling (body-on-particle, solved
Kutta wake row, or free-wake row). It is intrinsic to the FLOWVPM SFS model
`SFS_Cd_twolevel_nobackscatter`. In the trailing-edge wake region, the dynamic
coefficient `C[1]` is silently nullified by the `clipping_backscatter` strategy
on every particle where the stretching SFS points opposite to Γ — which is
precisely where the model is supposed to regularize. With `C[1]=0` the SFS term
in the RK3 update collapses, leaving inertial vortex stretching `J·Γ` unopposed.
Γ then roughly doubles per step on those particles and the runaway spreads
through the neighborhood.

The decisive line of code is
[`FLOWVPM/src/FLOWVPM_subfilterscale.jl:233`](../../FLOWVPM.jl/src/FLOWVPM_subfilterscale.jl):

```julia
if clipping(pfield, i)
    # Clip SFS model by nullifying the model coefficient
    pfield.particles[C_INDEX[1],i] *= 0
end
```

with `clipping = clipping_backscatter` defined at lines 282–296 of the same file.

## Symptoms (from canonical run + paraview)

In `data/rotor_hover_pressure_comparison/` the simulation diverges over the
final ~5 saved snapshots:

1. `get_SFS` grows first
2. `get_J` (velocity gradient) grows next
3. `get_Gamma` is last, starting in a single particle near the blade trailing
   edge and infecting neighbors over the last few steps

Sanity controls reported by the user:

| viscous | SFS                              | outcome              |
| ------- | -------------------------------- | -------------------- |
| on      | off                              | stable               |
| off     | off                              | marginally unstable, finite outcome |
| on      | `SFS_Cd_twolevel_nobackscatter`  | diverges (the case)  |

## Diagnostic procedure (what was run)

Two artifacts were created to do a forward-march restart of the last ~10
saved snapshots with full particle-state instrumentation:

- [`debug/sfs_blowup_tracer.jl`](sfs_blowup_tracer.jl) — a `SFSBlowupTracer`
  monitor that scans `wake.pfield` each step, tracks running medians/max of
  `|SFS|`, `|J|_F`, `|Γ|`, and on the first particle that exceeds
  `SFS_DEBUG_RATIO × running median` (default 10×):

    - prints the seed plus its k-NN neighborhood in full;
    - writes the seed + 16 nearest neighbors with complete particle state
      (`X`, `σ`, `Γ`, full 3×3 `J`, `SFS`, `C`, plus derived `J·Γ` and
      `Γ·SFS`) to `debug/logs/sfs_seed_step<N>_idx<I>.csv`;
    - throws `SFSBlowupDetected` to unwind `simulate!` cleanly before the
      OS reaps the process for OOM.

  Behavior is gated by `SFS_DEBUG=1`. Other env knobs: `SFS_DEBUG_VERBOSE`,
  `SFS_DEBUG_K`, `SFS_DEBUG_THRESHOLD_{SFS,J,GAMMA}`, `SFS_DEBUG_RATIO`,
  `SFS_DEBUG_HALT` (default 1), `SFS_DEBUG_FOCUS_INDEX`, `SFS_DEBUG_DUMP_DIR`.

- [`examples/rotor_hover_pressure_comparison_warmstart_debug.jl`](../examples/rotor_hover_pressure_comparison_warmstart_debug.jl)
  — driver that includes the canonical
  [`examples/rotor_hover_pressure_comparison.jl`](../examples/rotor_hover_pressure_comparison.jl)
  under `RHPC_SETUP_ONLY=1` (to source the geometry/frames/wake/monitor setup
  without firing the canonical `simulate!`), then calls
  `pnl.simulate_warmstart!` from the canonical PVD with the tracer appended
  to the monitor tuple. Output goes to a fresh
  `data/rotor_hover_pressure_comparison_warmstart_debug/` directory so the
  canonical run is preserved.

The single command that produced the diagnosis:

```bash
NREVS=10 SFS_DEBUG=1 SFS_DEBUG_RATIO=10 \
    julia --project examples/rotor_hover_pressure_comparison_warmstart_debug.jl
```

The run halted at step 42 with seed particle `idx=68` and wrote
`debug/logs/sfs_seed_step42_idx68.csv` (17 particles: seed + 16 NN).

## What the seed snapshot showed

Reading `debug/logs/sfs_seed_step42_idx68.csv`:

| quantity                | seed (idx=68) | top neighbors      |
| ----------------------- | -------------:| ------------------:|
| `\|SFS\|`               | 2.10e5        | 1.0e5–1.7e5        |
| `\|J\|_F`               | 1950          | 1700–2900          |
| `\|Γ\|`                 | 5.57e-4       | mostly 1e-4, several already 5e-4–1.3e-3 |
| `C[1]`                  | **0.0**       | **0.0** on every anomalous-SFS particle |
| `Γ·SFS`                 | **−66.2**     | strongly negative (−7 to −158) on every anomalous-SFS particle |
| `J·Γ` (per component)   | ~0.03–0.45    | ~0.01–0.5          |

Two facts stand out:

1. **`C[1] = 0` exactly** on every particle with anomalously high `|SFS|`.
   That is not a coincidence: it is the post-clip value.
2. **`Γ·SFS < 0`** on every one of those particles. The SFS stretching term
   `Estr_fmm` is pointing in the *backscatter* direction relative to Γ.

There is no small denominator. There is no particle pair at near-zero
separation. The failure is logical, not numerical.

## The mechanism, step by step

### Step 1 — `Estr_fmm` produces a backscatter-pointing SFS

`Estr_fmm` (FLOWVPM_subfilterscale_models.jl:16–40) populates each particle's
`get_SFS(P)` from the classic vortex-stretching tensor times neighbors'
Γ-weighted kernel ζ\_σ. In the high-strain region just downstream of the
trailing edge, the resulting vector field has `Γ·SFS < 0` for a large fraction
of particles — physically, the model wants to inject energy back to small
scales in this turbulent region.

### Step 2 — `pseudo3level_positive_afterUJ` computes `C[1]` and forces it ≥ 0

In `dynamicprocedure_pseudo3level_afterUJ`
(FLOWVPM_subfilterscale.jl:540–656):

```julia
# Store model coefficient
C_p[1] = C_p[2] / C_p[3]                       # nume / deno
...
# Force the coefficient to be positive
C_p[1] *= sign(C_p[1])^force_positive          # → |C_p[1]| when force_positive=true
```

`SFS_Cd_twolevel_nobackscatter` uses `pseudo3level_positive_afterUJ`, so
`force_positive=true` and `C_p[1]` is left as `|nume/deno| ≥ 0`.

### Step 3 — `clipping_backscatter` zeros `C[1]` on backscatter particles

The `DynamicSFS(::AfterUJ)` dispatch
(FLOWVPM_subfilterscale.jl:216–268) iterates the `clippings` tuple:

```julia
for clipping in SFS.clippings
    ...
    if clipping(pfield, i)
        # Clip SFS model by nullifying the model coefficient
        pfield.particles[C_INDEX[1],i] *= 0      # ← line 233 / 243
    end
    ...
end
```

`SFS_Cd_twolevel_nobackscatter` is constructed with
`clippings=(clipping_backscatter,)`. The predicate
(FLOWVPM_subfilterscale.jl:287–296) is:

```julia
function clipping_backscatter(pfield, i::Int)
    C  = pfield.particles[C_INDEX[1], i]
    G1, G2, G3 = pfield.particles[GAMMA_INDEX[1..3], i]
    S1, S2, S3 = pfield.particles[SFS_INDEX[1..3],  i]
    return C*(G1*S1 + G2*S2 + G3*S3) < 0
end
```

Combined with `force_positive=true`, the product
`C · (Γ·SFS)` is negative iff `Γ·SFS < 0` (with `C > 0`). On every such
particle the next statement nullifies `C[1]`. After this point,
`get_C(p)[1] == 0` on the entire backscatter region of the wake.

### Step 4 — RK3 RHS collapses to pure inertial stretching

`FLOWVPM_timeintegration.jl:276`–`300`:

```julia
C::R = get_C(p)[1]            # ← 0 for the offending particles
...
# Classic scheme (Γ⋅∇)U - Cϵ
M[4] = a*M[4] + dt*(J[1]*G[1]+J[4]*G[2]+J[7]*G[3] - C*get_SFS1(p)*get_sigma(p)[]^3/zeta0)
M[5] = a*M[5] + dt*(J[2]*G[1]+J[5]*G[2]+J[8]*G[3] - C*get_SFS2(p)*get_sigma(p)[]^3/zeta0)
M[6] = a*M[6] + dt*(J[3]*G[1]+J[6]*G[2]+J[9]*G[3] - C*get_SFS3(p)*get_sigma(p)[]^3/zeta0)
```

With `C=0` the SFS term vanishes and Γ updates purely from inertial vortex
stretching `J·Γ`.

### Step 5 — Inertial stretching dominates dt and Γ doubles each step

At the seed neighborhood: `|J| ≈ 2000`, `|Γ| ≈ 10⁻⁴` →
`|J·Γ|·dt ≈ 0.2 × dt`. With `dt = 1/3600 s ≈ 2.8×10⁻⁴`, that is
`ΔΓ ≈ 6×10⁻⁵` per step, comparable to the current `|Γ|` itself. Without an
SFS or viscous brake of comparable magnitude, Γ doubles roughly every step on
the affected particles.

### Step 6 — The instability spreads

Growing `|Γ|` on one particle linearly amplifies the FMM-induced velocity
gradient `|J|` it produces on its neighbors. Neighbors' RK3 RHS therefore
sees larger `|J|`, their `|J·Γ|` rises in turn, and the runaway propagates.
This matches the user's paraview observation that the blow-up starts in one
particle near the trailing edge and infects its surroundings.

### Step 7 — Why viscous alone is enough but viscous + this SFS is not

The CoreSpreading viscous model
(`viscous=CoreSpreading(wake_nu, wake_core_size, zeta_fmm; beta=1.5)`)
grows `σ` by `σ ← √(σ² + 2ν·dt)`. At trailing-edge particles σ is still
small (5–8×10⁻³) so the diffusion timescale is slower than the
stretching-doubling timescale; viscous alone *would* lose this race, but
without SFS interference the stretching is too small to start the chain.
With `SFS_Cd_twolevel_nobackscatter`, the model's clip removes the only
fast regularization on exactly the particles that needed it, and viscous
cannot catch up in time.

This is why the user's three sanity runs behave the way they do:

- **viscous on + SFS off**: clip path never executes; the model has no
  bug-prone surface; viscous handles regularization. Stable.
- **viscous off + SFS off**: stretching is weakly unregularized everywhere
  but still small in absolute terms. Marginally unstable, finite outcome.
- **viscous on + `SFS_Cd_twolevel_nobackscatter`**: the clip silently
  disables SFS in the backscatter region, while viscous is too slow.
  Diverges.

## Recommended fixes (not implemented)

In order of cleanliness:

1. **FLOWVPM**: replace the "nullify C" clip in `clipping_backscatter` with
   "project SFS onto the forward-scatter component" — subtract the
   `min(Γ·SFS, 0)·Γ̂` piece from SFS so the model still regularizes the
   forward-scatter part instead of abdicating. This preserves the no-backscatter
   intent of the model without the all-or-nothing failure mode.
2. **User workaround**: configure a DynamicSFS variant *without*
   `clipping_backscatter` in its `clippings` tuple (e.g. allow backscatter)
   and let the dynamic coefficient procedure self-regulate. The implementation
   already enforces magnitude bounds `minC`/`maxC` independently.
3. **User workaround**: keep SFS off and rely on CoreSpreading viscous alone
   — already known to be stable for this case.

The linear-doublet vs constant-doublet panel question raised separately is
orthogonal here: the seed particle is in the wake, not on a panel edge, and
the failure is purely intra-VPM. Linear panels would not change this.

## Artifacts

- [`debug/sfs_blowup_tracer.jl`](sfs_blowup_tracer.jl) — the tracer monitor
- [`examples/rotor_hover_pressure_comparison_warmstart_debug.jl`](../examples/rotor_hover_pressure_comparison_warmstart_debug.jl) — warmstart driver
- [`examples/rotor_hover_pressure_comparison.jl`](../examples/rotor_hover_pressure_comparison.jl) — patched with `RHPC_SETUP_ONLY` gate
- `debug/logs/sfs_seed_step42_idx68.csv` — seed snapshot from the
  diagnosis run (regenerated by re-running the warmstart script)

## Reproducing

From the FLOWPanel.jl repo root, with the `fastmultipole` branch checked out:

```bash
NREVS=10 SFS_DEBUG=1 SFS_DEBUG_RATIO=10 \
    julia --project examples/rotor_hover_pressure_comparison_warmstart_debug.jl
```

The run will warmstart from ~10 saved snapshots before the end of the
canonical PVD, march forward, and halt at the first seed detection with a
CSV dump in `debug/logs/`.

## Stage 1 / 2 ablation outcomes

Two follow-up runs were executed against the same warmstart driver and
tracer to isolate (a) whether the panel-wake-row → particle Hessian is a
meaningful driver of `|J|` (relevant to the constant-doublet vs
linear-doublet question), and (b) whether removing
`clipping_backscatter` is a viable fix in isolation. Env knobs:

- `WAKEROW_NO_HESSIAN_TO_PARTICLES=1` — splits the
  `_steady_aerodynamics!` wake-source `influence!` call into panel-part
  (rings + closing filaments) with `velocity_gradient=false` for
  `ParticleField` targets, and pfield-part with the standard VPM
  Hessian. Implemented as an env-gated branch around
  `FLOWPanel_simulate.jl:251`.
- `SFS_NO_BACKSCATTER_CLIP=1` — at the
  `pnl.PanelParticleWake(...; SFS=...)` call in the canonical script,
  replaces `SFS_Cd_twolevel_nobackscatter` with an equivalent
  `DynamicSFS(Estr_fmm, pseudo3level_beforeUJ, pseudo3level_positive_afterUJ; alpha=0.999, clippings=())`.
  No FLOWVPM source edit needed.

### Step-42 trigger comparison

All three configurations trip the tracer at step 42 with the same seed
particle `idx=68` and the same `|SFS|≈2.1×10⁵`. They differ in the
global `|J|` max and in the long-term trajectory:

| run                                   | step-42 `\|J\|` max | catastrophic step | notes                                  |
| ------------------------------------- | ------------------:| -----------------:| -------------------------------------- |
| baseline (clip ON, wake-H ON)         | 11451              | ~50               | OOM-killed near step 50                |
| Stage 1 (clip ON, wake-H OFF)         | 2056 (5.5× ↓)      | ~89               | OOM-killed at step 90                  |
| Stage 2 (clip OFF, wake-H ON), HALT=0 | 11467              | ~45               | 3-step ignition: \|J\|→58k, \|SFS\|→1.7×10⁶ |

### What Stage 1 says about the linear-panel question

Disabling the panel-wake-row → particle Hessian cut the **global** `|J|`
max at step 42 by 5.5× (11451 → 2056), confirming that some particles
near panel-row edges do experience the `~1/r²` filament singularity that
linear/sheet panels would smooth. But the **seed particle's own**
`|J|≈1949` is unchanged (it sits in pfield-on-pfield interaction, not
near a panel edge), so the seed signature at step 42 is identical.

Net effect on the blow-up: not eliminated, but **delayed by ~45 steps**.
Catastrophic ignition moves from step ~50 to step ~89. Linear panels
would therefore buy time but not fix the underlying instability — that
requires addressing the SFS-clip / `force_positive` interaction.

### What Stage 2 says about the clipping_backscatter ablation

Disabling the clip restores `C[1]` to non-zero values
(seed `C[1]=1.0` saturating at `maxC=1`; neighbors `0.0005–0.48`). But
because `Γ·SFS < 0` is the same pathology (the SFS field still points
in the backscatter direction) and `force_positive=true` forces `C ≥ 0`,
the RK3 SFS term `−C·SFS·σ³/ζ0` now has a component aligned **with** Γ
(since SFS ≈ −α·Γ gives `−C·SFS ≈ +C·α·Γ`). That is energy injection on
top of the inertial stretching `J·Γ`.

Result: the blow-up gets **worse**, not better. Within 3 steps of the
step-42 detection, `|J|` max jumps to 58k and `|SFS|` max to 1.7×10⁶.

### Revised recommendation

The two "user workaround" entries in the previous Recommended Fixes
section are now disproven by Stage 2:

- ~~Configure DynamicSFS without `clipping_backscatter`~~ — makes it
  worse for this case, not better.
- Keeping SFS off + viscous on is the only known-stable configuration
  (unchanged).

The real fix has to live in the SFS model itself. The decisive failure
is **not** the clip in isolation; it is the **combination** of:

1. `Estr_fmm` producing a backscatter-direction SFS in the trailing-edge
   wake (`Γ·SFS<0` on a wide neighborhood);
2. `pseudo3level_positive_afterUJ` forcing `C[1]=|C[1]|` (line 655);
3. the all-or-nothing `clipping_backscatter` strategy (line 233).

With (2) in force, neither setting of (3) is safe: clip on → no
regularization on backscatter particles; clip off → backscatter
injection. The clean fix needs to remove `force_positive=true` and
allow `C[1]` to be negative when the model wants backscatter, *or*
project SFS onto the forward-scatter component so the magnitude is
preserved but the direction always dissipates Γ.

### Linear panels — revised verdict

Stage 1 shows linear panels are not as orthogonal to this instability as
the earlier note in this file claimed. They wouldn't eliminate the
blow-up (the seed sits in pfield-on-pfield interaction), but they could
significantly slow it by removing panel-edge `1/r²` `|J|` spikes on
particles that drift close to ring boundaries. Whether that buys enough
time for viscous diffusion to catch up depends on the resolved
configuration — not answerable from this diagnostic alone. A separate
plan would scope the kernel refactor.

### Artifacts added by these stages

- `FLOWPanel_simulate.jl:256` — env-gated wake-source split branch
- `examples/rotor_hover_pressure_comparison.jl` (around line 142) — env-gated SFS variant constructor
- `/tmp/sfs_stage1.log`, `/tmp/sfs_stage2.log`,
  `/tmp/sfs_stage2_halt0.log` — tracer outputs (transient, not checked in)
- `debug/logs/sfs_seed_step42_idx68.csv` — last-write was Stage 2 (`C[1]`
  non-zero, `Γ·SFS<0`); rerun under desired configuration to refresh
