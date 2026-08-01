#=##############################################################################
# code_audit Task 2 addendum — particle wake with a SHORT panel buffer
#
# The base task2 script used the example's default panel_wake_rows
# (wake_length/das = 160 rows ~ 80 chords), so the particle-model run kept the
# entire trimmed wake as panels and culled every converted particle at the
# downstream boundary (np=62 at the end): the particle path was not exercised.
# This addendum reruns ONLY the particle model with panel_wake_rows=10 (a
# 5-chord buffer, matching Task 1 case c), so the 5->8.75-chord band is genuine
# particles influencing the body, and compares per-station amplitude/phase
# against the panel-wake histories already stored in task2_histories.csv.
#
# Outputs: code_audit/results/task2/{task2b_station_fits.csv,
#          task2b_histories.csv}
# Run:  julia --project code_audit/scripts/task2b_particle_short_buffer.jl
=###############################################################################

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
import Printf: @printf
import DelimitedFiles: readdlm

include(joinpath(@__DIR__, "..", "..", "examples", "pitching_wing.jl"))

const RESULTS = joinpath(@__DIR__, "..", "results", "task2")

# parameters identical to task2_unsteady_panel_vs_particle.jl except the buffer
n_span, n_airfoil, n_endcap = 15, 31, 5
n_section_bins = 30
section_eta = [0.15, 0.30, 0.45, 0.60, 0.75, 0.90]
alpha_mean_deg = 3.94
alpha_amp_deg  = 1.99
frequency_hz   = 4.01
n_cycles       = 3
c_per_dt       = 0.5
wake_length_spans = 2.0
panel_wake_rows = 10          # 5-chord buffer (rows are 0.5c long)

sim = prepare_pitching_wing(;
    wake_model=:particle,
    n_span, n_airfoil, n_endcap, n_section_bins,
    section_eta,
    alpha_mean_deg, alpha_amp_deg, frequency_hz,
    n_cycles, c_per_dt, wake_length_spans,
    panel_wake_rows,
    include_static_polar=false,
    save_vtk=false,
)

telapsed = @elapsed pnl.simulate!((sim.wing,), (sim.wake,), sim.frames, sim.maneuver!,
    sim.Uinf, sim.t_range;
    body_solvers=(sim.solver,),
    backend=sim.backend,
    monitors=sim.monitors,
    path=nothing,
    name="task2b_particle_short",
    set_Das_eta_freestream=NaN,
    verbose=false,
)
np = try pnl.FLOWVPM.get_np(sim.wake.pfield) catch; -1 end
t = collect(sim.t_range)
CL = collect(sim.force_monitor.force[3, :])
section = copy(sim.section_monitor.cl)
@printf("particle short-buffer run done in %.1f s (%d steps, np=%d)\n",
    telapsed, length(t), np)

# ------------------------------------------------- load stored panel histories
hist, hdr = readdlm(joinpath(RESULTS, "task2_histories.csv"), ','; header=true)
hdr = vec(hdr)
col(name) = findfirst(==(name), hdr)
t_ref = Float64.(hist[:, col("time")])
@assert isapprox(t_ref, t; atol=1e-10) "time grids differ from stored panel run"
panel_series = Dict{String,Vector{Float64}}()
panel_series["CL"] = Float64.(hist[:, col("panel_CL")])
for eta in section_eta
    panel_series["cl(eta=$(eta))"] = Float64.(hist[:, col("panel_cl_eta$(eta)")])
end

# ------------------------------------------------------------------------ fits
function fit_sinusoid(t, y, omega, idx)
    tt = t[idx]; yy = y[idx]
    M = hcat(ones(length(tt)), sin.(omega .* tt), cos.(omega .* tt))
    coef = M \ yy
    a0, as, ac = coef
    A = hypot(as, ac)
    phase = atan(ac, as)
    resid_rel = sqrt(sum(abs2, M * coef .- yy) / length(yy)) / A
    return (; mean=a0, amp=A, phase_deg=rad2deg(phase), resid_rel)
end
wrapdeg(x) = mod(x + 180, 360) - 180

omega, period = sim.setup.omega, sim.setup.period
cyc(n) = findall(ti -> (n - 1) * period <= ti <= n * period, t)
idx3 = cyc(3)

labels = vcat(["CL"], ["cl(eta=$(eta))" for eta in section_eta])
this_series = vcat([CL], [section[i, :] for i in eachindex(section_eta)])

println("\n=== Cycle-3 fits: panel (160-row buffer, stored) vs particle (10-row buffer) ===")
@printf("%-14s | %10s %10s %8s | %10s %10s %8s | %8s %8s\n",
    "station", "A_panel", "A_part10", "dA %", "ph_panel", "ph_part10", "dph deg",
    "mean_pnl", "mean_p10")
rows = []
for (lab, yq) in zip(labels, this_series)
    fp = fit_sinusoid(t, panel_series[lab], omega, idx3)
    fq = fit_sinusoid(t, yq, omega, idx3)
    dA = 100 * (fq.amp - fp.amp) / fp.amp
    dph = wrapdeg(fq.phase_deg - fp.phase_deg)
    @printf("%-14s | %10.6f %10.6f %+8.3f | %10.4f %10.4f %+8.4f | %8.5f %8.5f\n",
        lab, fp.amp, fq.amp, dA, fp.phase_deg, fq.phase_deg, dph, fp.mean, fq.mean)
    push!(rows, (; lab, fp, fq, dA, dph))
end

open(joinpath(RESULTS, "task2b_station_fits.csv"), "w") do io
    println(io, "station,amp_panel,amp_particle10,damp_pct,phase_panel_deg," *
        "phase_particle10_deg,dphase_deg,mean_panel,mean_particle10," *
        "resid_rel_panel,resid_rel_particle10")
    for r in rows
        println(io, join([r.lab, r.fp.amp, r.fq.amp, r.dA, r.fp.phase_deg,
            r.fq.phase_deg, r.dph, r.fp.mean, r.fq.mean,
            r.fp.resid_rel, r.fq.resid_rel], ","))
    end
end
open(joinpath(RESULTS, "task2b_histories.csv"), "w") do io
    hdr2 = vcat(["step", "time"],
        vcat(["particle10_CL"], ["particle10_cl_eta$(eta)" for eta in section_eta]))
    println(io, join(hdr2, ","))
    for i in eachindex(t)
        vals = Any[i - 1, t[i], CL[i]]
        append!(vals, section[:, i])
        println(io, join(vals, ","))
    end
end
@printf("\nnp at end = %d; wrote task2b CSVs under %s\n", np, RESULTS)
