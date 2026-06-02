#=##############################################################################
# DESCRIPTION
    Validation script: compare solvers BackslashDirichlet and Fast Gauss-Seidel.
    After 10 timesteps the wake strengths should match within 1e-12, if the leaf size is large enough.

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
                bodyoptargs=(cp_outer=true, kerneloffset=1e-2,
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
n_steps = 2
t_range = range(0.0, step=dt_val, length=n_steps)

expansion_order      = 10
multipole_acceptance = 0.4
leaf_size            = 20
backend = pnl.FastMultipoleBackend(; expansion_order, multipole_acceptance, leaf_size)
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

# ======================== (A) backslash dirichlet =============================
# println("\n" * "="^60)
# println("(A) ParticleWake — immediate particle shedding")
# println("="^60)

body_A = make_body()
body_A.Das[1] .= repeat(Vinf * dt_val * eta, 1, size(body_A.Das[1], 2))
wake_A = pnl.PanelParticleWake(body_A; max_particles=50000, nwakerows=1,
                            method_trailing=pnl.OverlapPPS(2.0, 1),
                            method_unsteady=pnl.OverlapPPS(2.0, 3))
frames_A = make_frames(body_A)
solver_A = pnl.Backslash(body_A)

pnl.simulate!((body_A,), (wake_A,), frames_A, maneuver, Uinf, t_range;
    body_solvers=(solver_A,), backend, verbose=false, path=save_path, name="bsd")

strength_A = copy(wake_A.panel_wake.strength[1])

# ======================== (B) Hybrid panel + VPM =============================
# println("\n" * "="^60)
# println("(B) ParticleWake — immediate particle shedding")
# println("="^60)

body_B = make_body()
body_B.Das[1] .= repeat(Vinf * dt_val * eta, 1, size(body_B.Das[1], 2))
wake_B = pnl.PanelParticleWake(body_B; max_particles=50000, nwakerows=1,
                            method_trailing=pnl.OverlapPPS(2.0, 1),
                            method_unsteady=pnl.OverlapPPS(2.0, 3))
frames_B = make_frames(body_B)
solver_B = pnl.FGSSolver(body_B;
            max_iterations=500,
            tolerance=1.0e-6,
            rlx=1.0,
            expansion_order=14,
            multipole_acceptance,
            leaf_size=15000,
            shrink=true,
            recenter=false,
            inner_iterations=20,
            reverse_pass=false,
            verbose=false
        )

pnl.simulate!((body_B,), (wake_B,), frames_B, maneuver, Uinf, t_range;
    body_solvers=(solver_B,), backend, verbose=false, path=save_path, name="fgs")

strength_B = copy(wake_B.panel_wake.strength[1])

# error 
err = norm(strength_A - strength_B) / norm(strength_A)
println("Relative error: ", err)
