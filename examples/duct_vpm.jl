#=##############################################################################
# DESCRIPTION
    Validation script: compare panel-only wake (nwakerows=10) with hybrid
    panel+VPM wake (nwakerows=3 + vortex particles). After 10 timesteps the
    body strengths should match within 5%.

# AUTHORSHIP
  * Author    : Ryan Anderson
  * Created   : Mar 2026
  * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
import GeometricTools as gt
include(joinpath(pnl.examples_path, "helper_functions.jl"))
import CSV
import DataFrames: DataFrame
using FLOWPanel.FastMultipole.StaticArrays
import LinearAlgebra: norm

run_name = "vpm"
save_path = joinpath("data", run_name)

# ----------------- GEOMETRY (shared) -----------------------------------------
filename   = joinpath(pnl.examples_path, "data", "naca662015.csv")
contour    = CSV.read(filename, DataFrame)

aspectratio = 0.6
d           = 2 * 0.835
magVinf     = 30.0
AOA         = 5.0       # degrees
n_rfl       = 6
NDIVS_theta = 21

NDIVS_rfl = [(0.25, n_rfl, 10.0, false),
             (0.50, n_rfl,  1.0, true),
             (0.25, n_rfl,  0.1, false)]

kernel   = Union{pnl.ConstantSource, pnl.VortexRing}
bodytype = pnl.RigidWakeBody{kernel}

function make_body()
    xs, ys = gt.rediscretize_airfoil(contour[:, 1], contour[:, 2],
                                          NDIVS_rfl, NDIVS_rfl; verify_spline=false)
    ys[end] = ys[1]
    xs *= d * aspectratio
    ys *= d * aspectratio
    ys .+= d / 2
    points = hcat(xs, ys)
    return generate_revolution_liftbody(bodytype, points, NDIVS_theta;
                bodyoptargs=(cp_outer=true, core_size=1e-2,
                             kernelcutoff=1e-14,
                             characteristiclength=(args...)->d*aspectratio,
                             semiinfinite_wake=false))
end

# Shared parameters
Vinf    = magVinf * [cosd(AOA), 0.0, sind(AOA)]
Uinf(t) = Vinf
eta     = 0.3
l       = d * aspectratio
dt_val  = magVinf / l / (n_rfl * 500)
n_steps = 401
t_range = range(0.0, step=dt_val, length=n_steps)

expansion_order      = 10
multipole_acceptance = 0.4
leaf_size            = 20
backend = pnl.FastMultipoleBackend(; expansion_order, multipole_acceptance, leaf_size)
# backend = pnl.DirectBackend()
maneuver = (args...; optargs...) -> nothing

# No rotation for this test
function make_frames(body)
    return pnl.ReferenceFrame(body;
        origin = SVector{3}(0.0, 0.0, 0.0),
        v = SVector{3}(0.0, 0.0, 0.0),
        # ω_axis = SVector{3}(0.0, 1.0, 0.0),
        # ω = 0.0,
        ω_axis = SVector{3}(0.0, 1.0, 0.0),
        ω = 0.1 * 2 * pi,
        R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name = "vehicle",
        child_index = Int[],
        dependent_index = [1]
    )
end

# ======================== (A) Panel-only reference ===========================
# println("="^60)
# println("(A) Panel-only wake with nwakerows=10")
# println("="^60)

# body_A = make_body()
# body_A.Das[1] .= repeat(Vinf * dt_val * eta, 1, size(body_A.Das[1], 2))
# wake_A = pnl.PanelWake(body_A; nwakerows=100)
# frames_A = make_frames(body_A)
# solver_A = pnl.Backslash(body_A)

# @time pnl.simulate!((body_A,), (wake_A,), frames_A, maneuver, Uinf, t_range;
#     body_solvers=(solver_A,), backend, verbose=true, path="panel_wake", name="test_stability")

# strength_A = copy(body_A.strength)

# ======================== (B) Hybrid panel + VPM =============================
println("\n" * "="^60)
println("(B) ParticleWake — immediate particle shedding")
println("="^60)

body_B = make_body()
body_B.Das[1] .= repeat(Vinf * dt_val * eta, 1, size(body_B.Das[1], 2))
wake_B = pnl.PanelParticleWake(body_B; max_particles=50000, nwakerows=1,
                            method_trailing=pnl.OverlapPPS(2.0, 1),
                            method_unsteady=pnl.OverlapPPS(2.0, 3))
frames_B = make_frames(body_B)
# solver_B = pnl.Backslash(body_B)
solver_B = pnl.FGSSolver(body_B;
            max_iterations=500,
            tolerance=1.0e-6,
            rlx=1.0,
            expansion_order=14,
            multipole_acceptance,
            leaf_size=150,
            shrink=true,
            recenter=false,
            inner_iterations=20,
            reverse_pass=false,
            verbose=false
        )

# t_range2 = range(0.0, step=dt_val/2, length=5)#n_steps*20)
l = d * aspectratio
dt = magVinf / l / (n_rfl * 500)
# t_range2 = range(0.0, step=dt, length=11) 
@time pnl.simulate!((body_B,), (wake_B,), frames_B, maneuver, Uinf, t_range;
    body_solvers=(solver_B,), backend, verbose=true, path=save_path, name="fgs")

strength_B = copy(body_B.strength)

# # ======================== Compare ============================================
# println("\n" * "="^60)
# println("Comparison")
# println("="^60)

# println("Particles shed: $(wake_B.pfield.np)")

# # Relative error (exclude near-zero strengths where relative error is meaningless)
# strength_scale = maximum(abs.(strength_A))
# significant = abs.(strength_A) .> 0.01 * strength_scale
# if any(significant)
#     rel_err = abs.(strength_A[significant] .- strength_B[significant]) ./ abs.(strength_A[significant])
#     sorted_err = sort(rel_err)
#     max_rel = maximum(rel_err)
#     mean_rel = sum(rel_err) / length(rel_err)
#     p95 = sorted_err[max(1, round(Int, 0.95 * length(sorted_err)))]
#     p99 = sorted_err[max(1, round(Int, 0.99 * length(sorted_err)))]

#     println("Significant points: $(sum(significant))/$(length(significant))")
#     println("Max  relative error: $(round(max_rel * 100; digits=2))%")
#     println("Mean relative error: $(round(mean_rel * 100; digits=2))%")
#     println("95th percentile:     $(round(p95 * 100; digits=2))%")
#     println("99th percentile:     $(round(p99 * 100; digits=2))%")

#     # Also compare with panel-only nwakerows=3 baseline
#     println("\n--- Panel-only nwakerows=3 baseline ---")
#     body_C = make_body()
#     body_C.Das[1] .= repeat(Vinf * dt_val * eta, 1, size(body_C.Das[1], 2))
#     wake_C = pnl.PanelWake(body_C; nwakerows=3)
#     frames_C = make_frames(body_C)
#     pnl.simulate!((body_C,), (wake_C,), frames_C, maneuver, Uinf, t_range;
#         body_solvers=(pnl.Backslash(body_C),), backend, verbose=false, path=nothing)
#     strength_C = body_C.strength
#     rel_err_C = abs.(strength_A[significant] .- strength_C[significant]) ./ abs.(strength_A[significant])
#     sorted_C = sort(rel_err_C)
#     println("Mean relative error: $(round(sum(rel_err_C)/length(rel_err_C)*100; digits=2))%")
#     println("95th percentile:     $(round(sorted_C[max(1, round(Int, 0.95*length(sorted_C)))]*100; digits=2))%")

#     if p95 < 0.05
#         println("\nPASS: 95th percentile relative error < 5%")
#     else
#         println("\nVPM improved mean error from $(round(sum(rel_err_C)/length(rel_err_C)*100; digits=1))% to $(round(mean_rel*100; digits=1))%")
#     end
# else
#     println("No significant strengths to compare")
# end
