#=##############################################################################
# code_audit Task 1 follow-up — attribute the b-vs-a ~11% CL deficit
#
# Checks:
#  1. Orientation confound: steady semiinfinite solve on the ROTATED body with
#     horizontal freestream (exact configuration of cases b/c) — should
#     reproduce CL_a=0.29386 if orientation is not the issue.
#  2. Das-magnitude sensitivity: panel-wake march (40 steps) with das x0.5 and
#     x2 — if CL moves materially, implicates Das/first-row wake geometry.
#  3. Circulation level: compare mid-span doublet-strength (mu) jump at the TE
#     (shed circulation) between the steady semiinfinite solve and the settled
#     panel-wake march — locates the deficit in the solve vs post-processing.
#
# Run:  julia --project code_audit/scripts/task1_followup_attribution.jl
=###############################################################################

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
import LinearAlgebra: cross, dot, norm
import Printf: @printf

include(joinpath(@__DIR__, "..", "..", "examples", "pitching_wing.jl"))

alpha_deg    = 4.5
dims  = _resolve_pitching_wing_dimensions(1.0, 4.0, nothing, nothing)
c, b  = dims.c, dims.b
U     = 330.2 * FT_TO_M
rho   = 1.225
Sref  = b * c
qinf  = 0.5 * rho * U^2
n_span, n_airfoil, n_endcap = 15, 31, 5
n_bins   = 30
c_per_dt = 0.5
dt       = c / U * c_per_dt
das0     = 0.05 * c
pivot    = SVector{3}(0.25 * c, 0.0, 0.0)
backend() = pnl.FastMultipoleBackend(expansion_order=8,
    multipole_acceptance=0.4, leaf_size=40)

# shed-circulation summary: total |Delta mu| along the TE from body wake strength
function shed_strength_summary(body)
    i_gamma = pnl.get_Gammai(body)
    vals = Float64[]
    for shedding in body.shedding
        for j in axes(shedding, 2)
            pi_upper = shedding[1, j]      # shedding panel index
            pi_lower = shedding[4, j]      # partner panel index (may be -1)
            mu_u = body.strength[pi_upper, i_gamma]
            mu_l = pi_lower > 0 ? body.strength[pi_lower, i_gamma] : 0.0
            push!(vals, mu_u - mu_l)
        end
    end
    return vals
end

function steady_rotated()
    body = build_pitching_wing_body(c, b;
        n_span, n_airfoil, n_endcap, semiinfinite_wake=true)
    frames = pitching_wing_frame(body, pivot, deg2rad(alpha_deg))
    Uinf = SVector{3}(U, 0.0, 0.0)
    set_wake_Das!(body, Uinf)
    Lhat = SVector{3}(0.0, 0.0, 1.0); Dhat = SVector{3}(1.0, 0.0, 0.0)
    Shat = SVector{3}(0.0, 1.0, 0.0)
    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false)
    force = pnl.ForceMonitor(1, 1; i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c),
        correct_kuttacondition=false, verbose=false)
    pnl.steady!(body, frames, Uinf;
        body_solvers=pnl.Backslash(body), backend=backend(),
        monitors=(pressure, force), path=nothing, verbose=false)
    CL = dot(SVector{3}(force.force[:, 1]), Lhat)
    return CL, shed_strength_summary(body)
end

function march_panel(nsteps; das=das0)
    wing = build_pitching_wing_body(c, b;
        n_span, n_airfoil, n_endcap, semiinfinite_wake=false)
    frames = pitching_wing_frame(wing, pivot, deg2rad(alpha_deg))
    Uinf(t) = SVector{3}(U, 0.0, 0.0)
    set_wake_Das!(wing, Uinf(0.0); magnitude=das)
    t_range = range(0.0, step=dt, length=nsteps + 1)
    steady_maneuver!(frames, systems, wakes, t) = nothing
    wake = pnl.PanelWake(wing; nwakerows=nsteps + 2, include_final_filament=false)
    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false)
    force = pnl.ForceMonitor(length(t_range), 1; i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c),
        correct_kuttacondition=false, verbose=false)
    pnl.simulate!((wing,), (wake,), frames, steady_maneuver!, Uinf, t_range;
        body_solvers=(pnl.Backslash(wing),), backend=backend(),
        monitors=(pressure, force), path=nothing,
        set_Das_eta_freestream=NaN, verbose=false)
    return collect(force.force[3, :]), shed_strength_summary(wing)
end

println("[1] Orientation cross-check: rotated body + horizontal Uinf, semiinfinite steady")
CL_a2, dmu_a = steady_rotated()
@printf("    CL_a2 = %.8g   (case-a baseline was 0.29385887)\n", CL_a2)

println("\n[2] Das-magnitude sensitivity, panel-wake march, 40 steps")
CLh_ref, dmu_b = march_panel(40; das=das0)
CLh_half, _ = march_panel(40; das=0.5 * das0)
CLh_x2, _ = march_panel(40; das=2.0 * das0)
@printf("    das=0.05c : CL(40) = %.8g\n", CLh_ref[end])
@printf("    das=0.025c: CL(40) = %.8g  (%+.3f%%)\n", CLh_half[end],
    100 * (CLh_half[end] - CLh_ref[end]) / CLh_ref[end])
@printf("    das=0.10c : CL(40) = %.8g  (%+.3f%%)\n", CLh_x2[end],
    100 * (CLh_x2[end] - CLh_ref[end]) / CLh_ref[end])

println("\n[3] Shed circulation (Delta-mu along TE), steady vs marched (40 steps)")
n = length(dmu_a)
mid = cld(n, 2)
@printf("    stations: %d | mid-span Delta-mu: steady=%.6g marched=%.6g (ratio %.4f)\n",
    n, dmu_a[mid], dmu_b[mid], dmu_b[mid] / dmu_a[mid])
@printf("    sum|Delta-mu|: steady=%.6g marched=%.6g (ratio %.4f)\n",
    sum(abs, dmu_a), sum(abs, dmu_b), sum(abs, dmu_b) / sum(abs, dmu_a))
