#=##############################################################################
# code_audit Task 2 — Unsteady panel vs particle wake (per-station amplitude/phase)
#
# Same wing as examples/pitching_wing.jl, gentle sinusoidal pitching (the
# validation case's own forcing: alpha = 3.94 + 1.99 sin(2*pi*4.01 t) deg,
# k = omega*c/(2U) ~ 0.038), identical dt/discretization for both wake models,
# unsteady PressureBernoulli. Uses prepare_pitching_wing() verbatim so the
# monitor stack matches the failing validation case; only wake_model differs.
#
# Per-station lift amplitude/phase extracted by least-squares sinusoid fit at
# the forcing frequency over the last full cycle (cycle 3); cycle-2 fit kept as
# washout/settledness evidence.
#
# Outputs: code_audit/results/task2/{task2_station_fits.csv, task2_histories.csv,
#          task2_cl_history.png, task2_amp_phase.png}
# Run:  julia --project code_audit/scripts/task2_unsteady_panel_vs_particle.jl
=###############################################################################

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
import Printf: @printf
import Statistics: mean

include(joinpath(@__DIR__, "..", "..", "examples", "pitching_wing.jl"))

const RESULTS = joinpath(@__DIR__, "..", "results", "task2")
mkpath(RESULTS)

# ------------------------------------------------------------------ parameters
# Discretization matches Task 1 (n_span=15, n_airfoil=31, n_endcap=5); forcing,
# dt, das, wake length match the pitching-wing validation defaults.
n_span, n_airfoil, n_endcap = 15, 31, 5
n_section_bins = 30
section_eta = [0.15, 0.30, 0.45, 0.60, 0.75, 0.90]
alpha_mean_deg = 3.94
alpha_amp_deg  = 1.99
frequency_hz   = 4.01
n_cycles       = 3
c_per_dt       = 0.5
wake_length_spans = 2.0

function run_model(wake_model::Symbol)
    sim = prepare_pitching_wing(;
        wake_model,
        n_span, n_airfoil, n_endcap, n_section_bins,
        section_eta,
        alpha_mean_deg, alpha_amp_deg, frequency_hz,
        n_cycles, c_per_dt, wake_length_spans,
        include_static_polar=false,
        save_vtk=false,
    )
    pnl.simulate!((sim.wing,), (sim.wake,), sim.frames, sim.maneuver!,
        sim.Uinf, sim.t_range;
        body_solvers=(sim.solver,),
        backend=sim.backend,
        monitors=sim.monitors,
        path=nothing,
        name="task2_$(wake_model)",
        set_Das_eta_freestream=NaN,
        verbose=false,
    )
    np = wake_model == :particle ?
        (try pnl.FLOWVPM.get_np(sim.wake.pfield) catch; -1 end) : 0
    return (; sim, setup=sim.setup,
        t=collect(sim.t_range),
        CL=collect(sim.force_monitor.force[3, :]),
        section=copy(sim.section_monitor.cl),
        n_particles=np)
end

# --------------------------------------------------------- sinusoid fit helper
"""Fit y ~ a0 + A*sin(omega*t + phi) over indices idx (least squares).
Returns (mean=a0, amp=A, phase_deg, resid_rel)."""
function fit_sinusoid(t, y, omega, idx)
    tt = t[idx]; yy = y[idx]
    M = hcat(ones(length(tt)), sin.(omega .* tt), cos.(omega .* tt))
    coef = M \ yy
    a0, as, ac = coef
    A = hypot(as, ac)
    phase = atan(ac, as)                    # y = a0 + A sin(wt + phase)
    resid = M * coef .- yy
    resid_rel = sqrt(sum(abs2, resid) / length(yy)) / A
    return (; mean=a0, amp=A, phase_deg=rad2deg(phase), resid_rel)
end

wrapdeg(x) = mod(x + 180, 360) - 180

# ------------------------------------------------------------------------- run
println("Task 2 — unsteady panel vs particle wake")
println("  forcing: alpha = $(alpha_mean_deg) + $(alpha_amp_deg) sin(2*pi*$(frequency_hz) t) deg")

tp = @elapsed panel = run_model(:panel)
@printf("  [panel]    done in %.1f s (%d steps)\n", tp, length(panel.t))
tq = @elapsed particle = run_model(:particle)
@printf("  [particle] done in %.1f s (%d steps, np=%d)\n", tq, length(particle.t),
    particle.n_particles)

setup = panel.setup
omega, period = setup.omega, setup.period
@assert panel.t == particle.t
t = panel.t
@printf("  steps/cycle = %.1f, k = omega*c/(2U) = %.4f\n",
    period / setup.dt, omega * setup.c / (2 * setup.U))

cyc(n) = findall(ti -> (n - 1) * period <= ti <= n * period, t)   # cycle-n window
idx2, idx3 = cyc(2), cyc(3)

# per-station + total fits, both models, cycles 2 and 3
labels = vcat(["CL"], ["cl(eta=$(eta))" for eta in section_eta])
series(m) = vcat([m.CL], [m.section[i, :] for i in eachindex(section_eta)])

rows = []
for (lab, yp, yq) in zip(labels, series(panel), series(particle))
    fp2, fp3 = fit_sinusoid(t, yp, omega, idx2), fit_sinusoid(t, yp, omega, idx3)
    fq2, fq3 = fit_sinusoid(t, yq, omega, idx2), fit_sinusoid(t, yq, omega, idx3)
    push!(rows, (; lab, fp2, fp3, fq2, fq3))
end

println("\n=== Cycle-3 fits (y = mean + A sin(wt + phi); phase rel. to alpha forcing) ===")
@printf("%-14s | %10s %10s %8s | %10s %10s %8s | %8s %8s\n",
    "station", "A_panel", "A_part", "dA %", "ph_panel", "ph_part", "dph deg",
    "mean_pnl", "mean_prt")
for r in rows
    dA = 100 * (r.fq3.amp - r.fp3.amp) / r.fp3.amp
    dph = wrapdeg(r.fq3.phase_deg - r.fp3.phase_deg)
    @printf("%-14s | %10.6f %10.6f %+8.3f | %10.4f %10.4f %+8.4f | %8.5f %8.5f\n",
        r.lab, r.fp3.amp, r.fq3.amp, dA, r.fp3.phase_deg, r.fq3.phase_deg, dph,
        r.fp3.mean, r.fq3.mean)
end

println("\n=== Washout evidence: cycle-2 -> cycle-3 change (same model) ===")
@printf("%-14s | panel: dA %% %9s dph deg | particle: dA %% %9s dph deg | resid3 pnl/prt\n", "station", "", "")
for r in rows
    @printf("%-14s | %+12.4f %12.5f | %+15.4f %12.5f | %.2e %.2e\n",
        r.lab,
        100 * (r.fp3.amp - r.fp2.amp) / r.fp3.amp, wrapdeg(r.fp3.phase_deg - r.fp2.phase_deg),
        100 * (r.fq3.amp - r.fq2.amp) / r.fq3.amp, wrapdeg(r.fq3.phase_deg - r.fq2.phase_deg),
        r.fp3.resid_rel, r.fq3.resid_rel)
end

# ------------------------------------------------------------------------- CSV
open(joinpath(RESULTS, "task2_station_fits.csv"), "w") do io
    println(io, "station,amp_panel,amp_particle,damp_pct,phase_panel_deg," *
        "phase_particle_deg,dphase_deg,mean_panel,mean_particle," *
        "amp_panel_c2,amp_particle_c2,phase_panel_c2_deg,phase_particle_c2_deg," *
        "resid_rel_panel,resid_rel_particle")
    for r in rows
        dA = 100 * (r.fq3.amp - r.fp3.amp) / r.fp3.amp
        dph = wrapdeg(r.fq3.phase_deg - r.fp3.phase_deg)
        println(io, join([r.lab, r.fp3.amp, r.fq3.amp, dA, r.fp3.phase_deg,
            r.fq3.phase_deg, dph, r.fp3.mean, r.fq3.mean,
            r.fp2.amp, r.fq2.amp, r.fp2.phase_deg, r.fq2.phase_deg,
            r.fp3.resid_rel, r.fq3.resid_rel], ","))
    end
end
open(joinpath(RESULTS, "task2_histories.csv"), "w") do io
    hdr = vcat(["step", "time", "alpha_deg"],
        ["$(m)_$(l)" for m in ("panel", "particle") for l in
            vcat(["CL"], ["cl_eta$(eta)" for eta in section_eta])])
    println(io, join(hdr, ","))
    for i in eachindex(t)
        alpha = alpha_mean_deg + alpha_amp_deg * sin(omega * t[i])
        vals = Any[i - 1, t[i], alpha]
        append!(vals, vcat([panel.CL[i]], panel.section[:, i]))
        append!(vals, vcat([particle.CL[i]], particle.section[:, i]))
        println(io, join(vals, ","))
    end
end
println("\nWrote CSVs under $(RESULTS)")

# ------------------------------------------------------------------------ plot
import PythonPlot as plt

fig, axs = plt.subplots(2, 1, figsize=(7.6, 6.4), sharex=true)
axs[0].plot(t ./ period, panel.CL; color="C0", label="panel wake")
axs[0].plot(t ./ period, particle.CL; color="C3", linestyle="--", label="particle wake")
axs[0].set_ylabel("CL"); axs[0].grid(true, alpha=0.35); axs[0].legend(fontsize=9)
axs[0].set_title("Task 2: pitching wing, panel vs particle wake")
axs[1].plot(t ./ period, particle.CL .- panel.CL; color="k")
axs[1].set_ylabel("CL particle - panel"); axs[1].set_xlabel("t / T")
axs[1].grid(true, alpha=0.35)
fig.tight_layout()
fig.savefig(joinpath(RESULTS, "task2_cl_history.png"), dpi=170)
plt.pyplot.close(fig)

amp_p = [r.fp3.amp for r in rows[2:end]]
amp_q = [r.fq3.amp for r in rows[2:end]]
ph_p = [r.fp3.phase_deg for r in rows[2:end]]
ph_q = [r.fq3.phase_deg for r in rows[2:end]]
fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.2))
axs[0].plot(section_eta, amp_p, "o-"; color="C0", label="panel")
axs[0].plot(section_eta, amp_q, "s--"; color="C3", label="particle")
axs[0].set_xlabel("η = 2y/b"); axs[0].set_ylabel("cl amplitude (cycle 3)")
axs[0].grid(true, alpha=0.35); axs[0].legend(fontsize=9)
axs[1].plot(section_eta, ph_p, "o-"; color="C0", label="panel")
axs[1].plot(section_eta, ph_q, "s--"; color="C3", label="particle")
axs[1].set_xlabel("η = 2y/b"); axs[1].set_ylabel("phase re. α forcing (deg)")
axs[1].grid(true, alpha=0.35); axs[1].legend(fontsize=9)
fig.tight_layout()
fig.savefig(joinpath(RESULTS, "task2_amp_phase.png"), dpi=170)
plt.pyplot.close(fig)
println("Wrote plots under $(RESULTS)")
