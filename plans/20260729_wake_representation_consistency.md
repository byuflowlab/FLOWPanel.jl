# Wake-representation consistency: freestream-convected `PanelWake` vs semi-infinite wake

**Status:** ready to execute. Written 2026-07-29 for a fresh agent — everything
needed is inline; no exploration required.

**Owning item:** `BRAINSTORM/014_first_wake_row_offset_selection.md`
**Related:** `BRAINSTORM/006_rotor_hover_converge_stable_near_reference.md`,
`logs/dji_convergence_20260722/phase_02e_unsteady_ct.md`

---

## 1. Why this exists

`Das` is the offset of the first wake row from the trailing edge, set as
`Das = max(η·dt·|V_te|, min_displacement)`.

It is the **largest sensitivity found anywhere in the rotor-hover study**. On the
DJI9443 production mesh (36,752 panels, 5400 RPM, NT=36, truncation 4R,
unfiltered, velocity formulation), settled 10-rev cycle-mean CT over revs 10–19:

| `Das` η | CT | note |
|---:|---:|---|
| 0.2 (floor removed) | 0.05215 ± 0.00078 | |
| 0.2 (with 0.01R floor) | 0.06148 ± 0.00152 | the long-standing default |
| 0.5 | 0.06942 ± 0.00185 | |
| 1.0 | 0.07133 ± 0.00159 | in 006's 0.068–0.072 band |
| 4.0 | 0.07190 ± 0.00101 | flat vs η=1.0 (+0.8%, < scatter) |

**+37% floor-free** (0.05215 → 0.07133), versus < cycle-scatter for mesh
refinement, `GreenReconstruction`, truncation depth and relaxation strength. That
single parameter moved CT from ~15% below experiment into the acceptance band.
**If `Das` is a free parameter, CT = 0.0713 is tuned, not predicted.**

Three mechanisms were proposed and **all refuted** (do not re-tread):
1. *"Larger `Das` pushes the near wake out of influence"* — refuted: the body
   always carries an attached doublet panel on `[TE, TE+Das]`.
2. *"η<1 seeds particles inside the attached panel's span"* — refuted:
   `update_TE!` places wake row 1 exactly at TE+`Das`; particles abut it.
3. *"`Das` < σ makes the shed particle core straddle the TE"* — refuted by
   measurement: CT is flat from η=1 where the core still exceeds `Das` 3.7×.

A Γ_TE diagnostic showed the effect acts through the **bound circulation**
(Γ tracks CT in lockstep; concentrated inboard where `Das`/chord is smallest),
which rules out "wake acts downstream on unchanged loading" but does not separate
a numerical (Kutta/paneling) route from a physical (induced-velocity) one.

### Why we do NOT lead with a `Das` sweep

A steady `Das` sweep with a uniform-Γ wake is **a theorem, not an experiment**.
Constant-strength doublet panels have **exactly cancelling interior edges** — a
uniform-Γ sheet is equivalent to a vortex loop around its outer boundary only.
So "attached strip on `[0,Das]` ∪ wake rows on `[Das,L]`" is analytically
identical to one uniform sheet on `[0,L]`, independent of where the internal
boundary sits, and tends to the semi-infinite wake as `L→∞`. Γ(`Das`) would be
constant by construction.

That identity exposes the sharper question: **do the code's two wake
representations actually agree?** Everything downstream — the whole η ladder —
presumes they do, and it has never been checked.

---

## 2. The experiment

Give `PanelWake` a mode where **all** rows convect with the freestream only, so
its geometry is *exactly* the semi-infinite wake's (a straight sheet along `U∞`).
Then check that it converges to the semi-infinite answer as wake length grows.

- **Converges** ⇒ representations are consistent; the `Das` question can be
  pursued inside that framework (gated work in §7).
- **Does not converge** ⇒ the discrepancy is a **direct candidate for the +37%**,
  and localizing it (§6) takes priority over every queued CT run.

This is also `BRAINSTORM/013`'s recommended first verification rung: an
impulsively-started wing with a freestream-convected wake *is* the Wagner
problem, whose CL must rise toward the steady value as the starting vortex
convects away.

---

## 3. `src/` change — new convection mode

**File:** `src/FLOWPanel_wake.jl`

### 3a. Current state (verbatim, so you need not go looking)

Struct, lines ~92-104:
```julia
struct PanelWake{TK,NK,TF} <: AbstractFreeWake
    nwakes::Array{Int, 0}
    nodes::Vector{Array{TF, 3}}
    strength::Vector{Array{TF, 3}}
    velocity::Vector{Array{TF, 3}}
    freestream::Vector{TF}
    core_size::Float64
    overflowed::Array{Bool, 0}
    shed_with_induced_velocity::Bool
    unsteady_filament::Bool
    include_final_filament::Bool
end
```

Constructor, line 124:
```julia
function PanelWake(shedding::Vector{Matrix{Int}}, kernel, TF=Float64;
        core_size=1e-3, nwakerows=100, shed_with_induced_velocity=true,
        unsteady_filament=true, include_final_filament=true)
```
(plus a convenience `PanelWake(body::AbstractLiftingBody, kernel=...; nwakerows, kwargs...)`)

`propagate!`, lines 487-501 — **two modes today, neither adequate**:
```julia
function propagate!(wake::PanelWake, dt; step=0, frames=nothing)
    for i_surf in eachindex(wake.nodes)
        if wake.shed_with_induced_velocity
            view(wake.velocity[i_surf], :, 1:wake.nwakes[]+1, :) .*= dt
            view(wake.nodes[i_surf], :, 1:wake.nwakes[]+1, :) .+= view(wake.velocity[i_surf], :, 1:wake.nwakes[]+1, :)
            view(wake.velocity[i_surf], :, 1:wake.nwakes[]+1, :) ./= dt
        else
            view(wake.nodes[i_surf], :, 1, :) .+= dt .* wake.freestream
            if wake.nwakes[] >= 1
                view(wake.velocity[i_surf], :, 2:wake.nwakes[]+1, :) .*= dt
                view(wake.nodes[i_surf], :, 2:wake.nwakes[]+1, :) .+= view(wake.velocity[i_surf], :, 2:wake.nwakes[]+1, :)
                view(wake.velocity[i_surf], :, 2:wake.nwakes[]+1, :) ./= dt
            end
        end
    end
end
```
`true` → all rows roll up. `false` → **only row 1** uses freestream; rows 2+ still
roll up. Neither keeps the whole sheet straight.

### 3b. The change

1. Add field `freestream_convection::Bool` to the struct (append last, update the
   inner constructor call at the end of `PanelWake(...)`).
2. Add kwarg `freestream_convection=false` — **default preserves behaviour
   exactly**; verify no existing test/example changes result.
3. New leading branch in `propagate!`:
```julia
if wake.freestream_convection
    view(wake.nodes[i_surf], :, 1:wake.nwakes[]+1, :) .+=
        dt .* reshape(wake.freestream, 3, 1, 1)
elseif wake.shed_with_induced_velocity
    ...unchanged...
else
    ...unchanged...
end
```
**Note the `reshape`**: the existing row-1 line broadcasts a 3-vector against a
`3×C` view; a `3×R×C` view needs explicit singleton dims.
4. Check `src/FLOWPanel_warmstart.jl:166-218` (restores `nwakes[]`/`overflowed[]`)
   and any other site that reconstructs a `PanelWake`, so the new field is carried.

### 3c. Unit test — `test/runtests_unit_wake.jl`

That file already hand-builds wakes (lines 30-135); follow its style. With
`freestream_convection=true`, after `k` calls to `propagate!(wake, dt)` every node
must equal its initial position `+ k·dt·freestream` **to round-off**, and rows
must remain exactly coplanar. Also assert the two legacy modes are byte-identical
to before the change.

---

## 4. Driver — extend `examples/suddenly_started_wing.jl`

**Do not write a new example.** This file already implements ~90% of the
experiment: impulsively-started wing, `PanelWake`, a `Das` knob, a semi-infinite
steady reference, and Wagner comparison. Existing results live in
`data/suddenly_started_wing/` (`convergence.csv`, `na21_ns*_dt0p125*`).

Relevant existing pieces:

| Location | What it does |
|---|---|
| `build_suddenly_started_wing(config; semiinfinite_wake)` | NACA 0012 open-tip wing, `bodytype = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}` (Neumann — matches the rotor study) |
| `_set_ssw_Das!(body, displacement)` (line 188) | `for Das in body.Das, j in axes(Das,2); Das[:,j] .= displacement; end` |
| `prepare_suddenly_started_wing(config)` (line ~222) | Builds the unsteady case. Currently: `_set_ssw_Das!(wing, config.eta*dt*full_uinf)`; `PanelWake(wing; nwakerows=length(t_range)-1, core_size=..., include_final_filament=false, shed_with_induced_velocity=true)` |
| `initialize_wake_convection!` monitor (line ~251) | At `i_step==0` seeds `velocity[:,1,:] .= full_uinf` to avoid a degenerate first wake panel |
| `_ssw_steady_cl(config, backend; ...)` (line ~264) | **The semi-infinite reference**: `semiinfinite_wake=true`, `_set_ssw_Das!(wing, full_uinf/norm(full_uinf))` (**unit** Das), `pnl.steady!`, returns `(; cl, cd, ...)` |
| `_ssw_solver(body, backend, max_backslash)` | `Backslash` if `ncells <= max_backslash` else `KrylovSolver` |
| `pnl.simulate!(...)` call (line ~349) | with `set_Das_eta_freestream=NaN`, `grad_mu_options=SSW_GRAD_MU_OPTIONS` (`(; basis=:tri)`) |

Existing ENV knobs: `SSW_MODE` (single|coarse|convergence), `SSW_BACKEND`
(direct|fmm), `SSW_NAIRFOIL`, `SSW_NSPAN`, `SSW_AR` (default 100),
`SSW_ETA` (**the `Das` knob**, default 0.25), `SSW_DT_STAR` (default 0.125),
`SSW_T_END_STAR` (default 7), `SSW_OUTPUT`, `SAVE_VTK`.

### Changes to make

1. Add `SSW_FREESTREAM_CONVECTION` (default `false`) to `SSWConfig`, threaded into
   the `PanelWake(...)` call in `prepare_suddenly_started_wing`.
2. Add `SSW_MODE=wake_consistency`: sweep wake length and compare against the
   semi-infinite reference from `_ssw_steady_cl`.
3. Write `data/suddenly_started_wing/wake_consistency.csv` with one row per
   wake length: `n_rows, L_over_c, cl_panel, cl_semiinf, rel_err, cd_panel,
   cd_semiinf, eta, freestream_convection, backend, elapsed_s`.

### Run matrix

Wake length grows one row per step, so **number of rows = number of steps**.
With `SSW_DT_STAR = 0.125`, `L/c = 0.125 · n_rows`. Sweep
`SSW_T_END_STAR ∈ {1, 2, 4, 8, 16, 32, 64}` ⇒ `L/c` from 1 to 64
(`n_rows` 8 → 512).

Run **both** `SSW_FREESTREAM_CONVECTION ∈ {false, true}`:
- `false` (rolled-up wake, current behaviour) — expected NOT to match the
  semi-infinite reference exactly, since a rolled-up wake is a different (more
  physical) geometry. This is the control.
- `true` (straight sheet) — **must** converge to the semi-infinite reference.
  This is the test.

Use `SSW_AR=100` (near-2D, cleanest Wagner comparison) for the headline, and
repeat the converged endpoint at `SSW_AR=6` to confirm it is not AR-specific.

### Backend: use `SSW_BACKEND=direct`

**Critical.** `src/FLOWPanel_liftingbody.jl:996` inflates the FMM panel radius by
`|Das|`:
```julia
if !system.semiinfinite_wake
    buffer[4, i_buffer] += update_radius   # update_radius = max(|Da|, |Db|)
end
```
so the multipole acceptance criterion is an **explicit function of `Das`** — FMM
would sweep its own error alongside the variable of interest.

**This also means every unsteady rotor η run used `FastMultipoleBackend` and
therefore carries this as an untested contributor to the +37%.** Record that in
`BRAINSTORM/014` as a fourth candidate mechanism; cheapest check is two points of
the η ladder rerun with `DirectBackend`.

---

## 5. Success criterion

With `freestream_convection=true`, `cl_panel(L) → cl_semiinf` monotonically, with
`rel_err < 1%` by `L/c ≈ 32` and still falling. Report the convergence order.

---

## 6. If they do NOT converge — diagnostics, ranked

1. **Streamwise Γ uniformity.** At the largest `L`, print
   `wake.strength[i][1, :, j]` along the streamwise index. At steady state every
   row must carry the same Γ (Kelvin). Non-uniformity ⇒ wake not converged, or
   shedding inconsistent.
2. **Attached strip vs first wake row strength — the prime suspect.** These come
   from *different code paths*: the strip radiates
   `get_strength_doublet(source_system, i_source)` + `wake_strength_shift[i_source]`
   inside the influence routine (`src/FLOWPanel_elements_fmm.jl:1030-1038`),
   while wake rows are set from `μ_upper − μ_lower` via
   `pnl._get_wakestrength_Gamma(body, i_shed, i_surf)`
   (`src/FLOWPanel_liftingbody.jl:1608`, used by `shed_wake!` at
   `src/FLOWPanel_simulate.jl:1076-1084`). Compare per shedding edge.
   **A mismatch is the finding** — it makes `Das` a lever on how much area
   carries the wrong Γ, explaining the +37% directly.
3. **Final-filament / starting-vortex treatment.** The current SSW config uses
   `include_final_filament=false` (strictly finite, panel-only). The docstring at
   `src/FLOWPanel_wake.jl:88-91` warns wake-length convergence must then be
   checked explicitly — which is exactly what this study does. Also try
   `include_final_filament=true` with `unsteady_filament=false`, which makes
   `_final_filament_strength` (`:1291-1294`) cancel the last row and emulate a
   semi-infinite closure. Assert
   `FastMultipole.get_n_bodies(pnl.FilamentWrapper(wake)) == n_shedding_edges` —
   it returns **0** unless `wake.overflowed[] == true` (`:1343`), silently
   reinstating a starting vortex.
4. **Wagner cross-check.** Plot `cl(L)/cl_semiinf` vs convected distance; it
   should rise from ~0.5 toward 1. Plateauing below 1 ⇒ missing far wake;
   overshoot ⇒ excess near-field.
5. **Geometry identity.** Assert `PanelWake` nodes lie on the ray `TE + s·d̂` to
   round-off — confirms §3 works and both representations describe one surface.

---

## 7. GATED follow-ons (do not start until §5 reports)

- **If it converges:** run the `Das` sweep inside this framework as an identity
  check with **predicted answer zero**, reporting `max_j|ΔΓ_j|/max_j|Γ_j|` as a
  multiple of the numerical floor. Calibrate the floor with a null test (perturb
  something that provably cannot matter; Γ must be bitwise identical). A clean
  null **proves** the +37% is not a steady paneling effect.
- **Either way — the experiment that can explain the +37%** is unsteady: sweep
  `Das` at fixed *large* `nwakerows`, so particles sit far downstream and the near
  wake is all panels at near-uniform Γ where the cancellation theorem applies.
  Collapse of the sensitivity ⇒ it is the **panel→particle representation change**
  (measured σ ≈ 1.5 local chords at shedding), routing to `BRAINSTORM/012`.
  `p2e_nrows4_das1p0` (η=1.0, nrows=4, CT ≈ 0.0744 vs 0.0713 at nrows=1) is
  already one point; it needs **η=0.2 at nrows=4** to form the pair.

---

## 8. Traps (each has already bitten this project)

1. **`Backslash` caches a `Das`-dependent operator.** `Backslash(body)`
   (`src/FLOWPanel_solver.jl:222-233`) assembles **and LU-factorizes** `G` at
   construction, and `_G!` → `induced(...)` walks the `TE→TE+Das` strip, so `G`
   depends on `Das`. **Rebuild the solver after any `Das` change.** This exact
   defect invalidated the unsteady `das1p0_refresh` run (job 12927924): the driver
   skipped `initialize_Das!`, so `G` was factorized with `Das = 0`, collapsing
   bound circulation 21.7× and CT to 0.0076.
2. **Semi-infinite requires a UNIT `Das`.** `_phi_semiinfinite`
   (`src/FLOWPanel_elements_fmm.jl:1219-1229`) *errors* if
   `|‖d‖²−1| > 2eps()`. Normalize and assert before assigning; `cosd/sind`
   round-off can trip it.
3. **`_accumulate_Das!` is `+=`-based** (`src/FLOWPanel_simulate.jl:3-19`). Zero
   with `_das_zero!` before re-deriving, or offsets accumulate.
4. **Shedding must come from the *constructed* body**, not the raw mesh
   (CLAUDE.md invariant). The SSW builder already handles this.
5. **Don't judge a running job.** Only terminal states count — a job that has
   passed a previous failure point is not evidence of survival.
6. **FMM vs `Das`** — see §4.

---

## 9. Verification before interpreting anything

1. Unit test of §3c passes; the two legacy convection modes are unchanged.
2. `_ssw_steady_cl` reproduces the existing `data/suddenly_started_wing/`
   results for the same config — validates geometry, shedding, monitors.
3. `SSW_FREESTREAM_CONVECTION=false` reproduces the existing Wagner curve
   byte-for-byte (default-preserving change).
4. Spanwise symmetry of Γ at solver tolerance.
5. Γ_TE cross-check: `pnl._get_wakestrength_Gamma` vs the direct `body.strength`
   form (see `examples/dji9443_circulation_audit.jl:230-245`) vs
   `BoundCirculationMonitor` — all agree to round-off; pins the sign convention.

## 10. Cost

`SSW_AR=100`, `SSW_NAIRFOIL=21`, `SSW_NSPAN=12` is small (well under the
`Backslash` cutoff). The `L/c=64` case is 512 steps — the long pole. Whole sweep
(7 lengths × 2 convection modes + references): **~1–3 h locally** with
`julia --project -t 6`. Start with `L/c ∈ {1,2,4,8}` to confirm the trend before
paying for the long cases.

Run:
```
SSW_MODE=wake_consistency SSW_BACKEND=direct SSW_FREESTREAM_CONVECTION=true \
  julia --project -t 6 examples/suddenly_started_wing.jl
```

## 11. Report into

- `BRAINSTORM/014_first_wake_row_offset_selection.md` — headline verdict, and add
  the FMM-radius-vs-`Das` coupling as a fourth candidate mechanism.
- `logs/dji_convergence_20260722/phase_02e_unsteady_ct.md` — run ledger.
- Update `BRAINSTORM/INDEX.md` outcome cell for 014.
