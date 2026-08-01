#=##############################################################################
# code_audit Task 1 — Steady wake-model consistency (spanwise lift)
#
# Same wing as examples/pitching_wing.jl at fixed alpha (validated range),
# steady conditions, STEADY PressureBernoulli everywhere, so only the wake
# model varies:
#   a. semiinfinite_wake=true, steady solve            (trusted baseline)
#   b. unsteady march, wake_model=:panel, semiinfinite_wake=false
#   c. unsteady march, wake_model=:particle
#
# Outputs: code_audit/results/task1/{task1_spanwise.csv, task1_cl_history.csv,
#          task1_overlay.png, task1_cl_history.png}
# Run:  julia --project code_audit/scripts/task1_steady_wake_consistency.jl
=###############################################################################

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
import LinearAlgebra: cross, dot, norm
import Printf: @printf
import Statistics: mean

# Reuse the pitching-wing generator and helpers (safe to include; main is guarded)
include(joinpath(@__DIR__, "..", "..", "examples", "pitching_wing.jl"))

const RESULTS = joinpath(@__DIR__, "..", "results", "task1")
mkpath(RESULTS)

# ------------------------------------------------------------------ parameters
alpha_deg    = 4.5                    # fixed alpha in the validated range
c_ft         = 1.0
aspect_ratio = 4.0
dims  = _resolve_pitching_wing_dimensions(c_ft, aspect_ratio, nothing, nothing)
c, b  = dims.c, dims.b
U     = 330.2 * FT_TO_M
rho   = 1.225
Sref  = b * c
qinf  = 0.5 * rho * U^2

n_span, n_airfoil, n_endcap = 15, 31, 5
n_bins   = 30
c_per_dt = 0.5
nsteps   = 80                          # wake convects ~nsteps*c_per_dt = 40 chords
dt       = c / U * c_per_dt
t_range  = range(0.0, step=dt, length=nsteps + 1)
das      = 0.05 * c                    # das_chord_fraction=0.05 as in the example
pivot    = SVector{3}(0.25 * c, 0.0, 0.0)

backend() = pnl.FastMultipoleBackend(expansion_order=8,
    multipole_acceptance=0.4, leaf_size=40)

# --------------------------------------------------------- case a: semiinfinite
"""Steady solve, semiinfinite wake — mirrors run_pitching_wing_static_polar."""
function run_case_a()
    body = build_pitching_wing_body(c, b;
        n_span, n_airfoil, n_endcap, semiinfinite_wake=true)
    Uinf = _uinf_from_alpha(U, alpha_deg)
    Lhat, Dhat, Shat = _lift_drag_span_directions(Uinf)
    set_wake_Das!(body, Uinf)

    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false)
    force = pnl.ForceMonitor(1, 1;
        i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c),
        correct_kuttacondition=false,
        verbose=false,
    )
    spanwise = pnl.SpanwiseLoadingMonitor(n_bins, 1;
        components=(lift=Lhat, drag=Dhat),
        span_axis=Shat,
        per_length=true,
        normalization=pnl.NoSectionalNormalization(),
        file=false,
        binning=:span_overlap,
    )

    pnl.steady!(body, pnl.ReferenceFrame(body), Uinf;
        body_solvers=pnl.Backslash(body),
        backend=backend(),
        monitors=(pressure, force, spanwise),
        path=nothing,
        name="task1_static",
        verbose=false,
    )

    eta = spanwise.bin_center ./ (b / 2)
    cl  = spanwise.load_components[1, :] ./ (qinf * c)
    CL  = dot(SVector{3}(force.force[:, 1]), Lhat)
    return (; eta, cl, CL)
end

# --------------------------------------------- cases b/c: unsteady march, fixed alpha
"""Unsteady march at constant Uinf/alpha with steady PressureBernoulli.
Mirrors prepare_pitching_wing's unsteady body path (rotated wing, horizontal
freestream) but with zero pitching rate."""
function run_case_unsteady(wake_model::Symbol; panel_wake_rows::Int,
        max_particles::Int=100_000)
    wing = build_pitching_wing_body(c, b;
        n_span, n_airfoil, n_endcap, semiinfinite_wake=false)
    frames = pitching_wing_frame(wing, pivot, deg2rad(alpha_deg))
    Uinf(t) = SVector{3}(U, 0.0, 0.0)
    set_wake_Das!(wing, Uinf(first(t_range)); magnitude=das)
    Lhat = SVector{3}(0.0, 0.0, 1.0)
    Dhat = SVector{3}(1.0, 0.0, 0.0)
    Shat = SVector{3}(0.0, 1.0, 0.0)

    steady_maneuver!(frames, systems, wakes, t) = nothing

    wake = if wake_model == :panel
        pnl.PanelWake(wing; nwakerows=panel_wake_rows, include_final_filament=false)
    elseif wake_model == :particle
        pnl.PanelParticleWake(wing;
            nwakerows=panel_wake_rows,
            max_particles,
            method_trailing=pnl.OverlapPPS(1.3, 2),
            method_unsteady=pnl.OverlapPPS(1.3, 2),
            particle_maintenance=pnl.ParticleMaintenance(),
        )
    else
        throw(ArgumentError("wake_model must be :panel or :particle"))
    end

    # STEADY PressureBernoulli (unsteady=false default) — kinematics are zero,
    # so the body-relative trace is exact and only the wake model varies.
    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false)
    force = pnl.ForceMonitor(length(t_range), 1;
        i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c),
        correct_kuttacondition=false,
        verbose=false,
    )
    spanwise = pnl.SpanwiseLoadingMonitor(n_bins, 1;
        components=(lift=Lhat, drag=Dhat),
        span_axis=Shat,
        per_length=true,
        normalization=pnl.NoSectionalNormalization(),
        file=false,
        binning=:span_overlap,
    )

    pnl.simulate!((wing,), (wake,), frames, steady_maneuver!, Uinf, t_range;
        body_solvers=(pnl.Backslash(wing),),
        backend=backend(),
        monitors=(pressure, force, spanwise),
        path=nothing,
        name="task1_$(wake_model)",
        set_Das_eta_freestream=NaN,
        verbose=false,
    )

    eta = spanwise.bin_center ./ (b / 2)
    cl  = spanwise.load_components[1, :] ./ (qinf * c)   # last-step distribution
    CL_hist = collect(force.force[3, :])

    # wake extent diagnostics (chords downstream of the TE)
    te_x = (1 - 0.25) * c + pivot[1]    # approx TE x-location after rotation
    panel_extent = try
        pw = wake_model == :panel ? wake : wake.panel_wake
        nw = pw.nwakes[]
        nw > 0 ? maximum(maximum(nodes[1, 1:nw+1, :]) for nodes in pw.nodes) : NaN
    catch
        NaN
    end
    particle_extent = if wake_model == :particle
        try
            pf = wake.pfield
            np = pnl.FLOWVPM.get_np(pf)
            np > 0 ? maximum(pnl.FLOWVPM.get_X(pf, i)[1] for i in 1:np) : NaN
        catch
            NaN
        end
    else
        NaN
    end
    n_particles = wake_model == :particle ?
        (try pnl.FLOWVPM.get_np(wake.pfield) catch; -1 end) : 0

    return (; eta, cl, CL=CL_hist[end], CL_hist,
        panel_extent_chords=(panel_extent - te_x) / c,
        particle_extent_chords=(particle_extent - te_x) / c,
        n_particles)
end

# ------------------------------------------------------------------------- run
println("Task 1 — steady wake-model consistency")
@printf("  wing: c=%.4g m, b=%.4g m (AR=%.3g), n_span=%d n_airfoil=%d n_endcap=%d\n",
    c, b, b / c, n_span, n_airfoil, n_endcap)
@printf("  alpha=%.2f deg, U=%.4g m/s, dt=%.4g s (%.2g chords/step), nsteps=%d\n",
    alpha_deg, U, dt, c_per_dt, nsteps)

println("\n[a] semiinfinite wake, steady solve...")
ta = @elapsed a = run_case_a()
@printf("    CL_a = %.8g  (%.1f s)\n", a.CL, ta)

println("\n[b] unsteady march, panel wake (no truncation: nwakerows=$(nsteps + 2))...")
tb = @elapsed bres = run_case_unsteady(:panel; panel_wake_rows=nsteps + 2)
@printf("    CL_b(final) = %.8g  (%.1f s); wake extent = %.1f chords\n",
    bres.CL, tb, bres.panel_extent_chords)

println("\n[c] unsteady march, particle wake (panel buffer 10 rows)...")
tc = @elapsed cres = run_case_unsteady(:particle; panel_wake_rows=10)
@printf("    CL_c(final) = %.8g  (%.1f s); panel extent = %.1f chords, particle extent = %.1f chords, np = %d\n",
    cres.CL, tc, cres.panel_extent_chords, cres.particle_extent_chords, cres.n_particles)

# ------------------------------------------------------------ settledness check
function settled_stats(CL_hist, nwindow=10)
    tail = CL_hist[end-nwindow+1:end]
    drift = (tail[end] - tail[1]) / abs(tail[end])
    ripple = (maximum(tail) - minimum(tail)) / abs(tail[end])
    return drift, ripple
end
db, rb = settled_stats(bres.CL_hist)
dc, rc = settled_stats(cres.CL_hist)
@printf("\nSettledness (last 10 steps): panel drift=%.3e ripple=%.3e | particle drift=%.3e ripple=%.3e\n",
    db, rb, dc, rc)

# --------------------------------------------------------------------- compare
# bins are identical across runs (same monitor construction + span axis)
@assert isapprox(a.eta, bres.eta; rtol=1e-10) && isapprox(a.eta, cres.eta; rtol=1e-10)
inboard = findall(eta -> abs(eta) <= 0.8, a.eta)
clmax = maximum(abs, a.cl)
dev_ba = maximum(abs.(bres.cl[inboard] .- a.cl[inboard])) / clmax
dev_cb = maximum(abs.(cres.cl[inboard] .- bres.cl[inboard])) / clmax
dev_ca = maximum(abs.(cres.cl[inboard] .- a.cl[inboard])) / clmax
@printf("\nTotal CL:  a=%.6g  b=%.6g  c=%.6g\n", a.CL, bres.CL, cres.CL)
@printf("CL rel diff: b vs a = %+.3f%%,  c vs a = %+.3f%%,  c vs b = %+.3f%%\n",
    100 * (bres.CL - a.CL) / a.CL, 100 * (cres.CL - a.CL) / a.CL,
    100 * (cres.CL - bres.CL) / bres.CL)
@printf("Max sectional deviation (|eta|<=0.8, normalized by max cl_a):\n")
@printf("  b vs a = %.3f%%   c vs b = %.3f%%   c vs a = %.3f%%\n",
    100dev_ba, 100dev_cb, 100dev_ca)

# ------------------------------------------------------------------------- CSV
open(joinpath(RESULTS, "task1_spanwise.csv"), "w") do io
    println(io, "eta,cl_a_semiinfinite,cl_b_panel,cl_c_particle")
    for i in eachindex(a.eta)
        println(io, "$(a.eta[i]),$(a.cl[i]),$(bres.cl[i]),$(cres.cl[i])")
    end
end
open(joinpath(RESULTS, "task1_cl_history.csv"), "w") do io
    println(io, "step,time,CL_b_panel,CL_c_particle,CL_a_semiinfinite")
    for i in eachindex(t_range)
        println(io, "$(i-1),$(t_range[i]),$(bres.CL_hist[i]),$(cres.CL_hist[i]),$(a.CL)")
    end
end
println("\nWrote CSVs under $(RESULTS)")

# ------------------------------------------------------------------------ plot
import PythonPlot as plt

fig, ax = plt.subplots(figsize=(7.2, 4.8))
ax.plot(a.eta, a.cl, "-k"; linewidth=2.0, label="a: semiinfinite (steady)")
ax.plot(bres.eta, bres.cl, "--"; color="C0", linewidth=1.6, label="b: panel wake (settled)")
ax.plot(cres.eta, cres.cl, ":"; color="C3", linewidth=1.8, label="c: particle wake (settled)")
ax.set_xlabel("η = 2y/b"); ax.set_ylabel("cl")
ax.grid(true, alpha=0.35); ax.legend(fontsize=9)
ax.set_title("Task 1: steady wake-model consistency, α=$(alpha_deg)°")
fig.tight_layout()
fig.savefig(joinpath(RESULTS, "task1_overlay.png"), dpi=170)
plt.pyplot.close(fig)

fig, ax = plt.subplots(figsize=(7.2, 4.8))
steps = 0:nsteps
ax.plot(steps, bres.CL_hist; color="C0", label="b: panel wake")
ax.plot(steps, cres.CL_hist; color="C3", label="c: particle wake")
ax.axhline(a.CL; color="k", linestyle="--", label="a: semiinfinite")
ax.set_xlabel("step"); ax.set_ylabel("CL")
ax.grid(true, alpha=0.35); ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(joinpath(RESULTS, "task1_cl_history.png"), dpi=170)
plt.pyplot.close(fig)
println("Wrote plots under $(RESULTS)")
