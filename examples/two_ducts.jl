#=##############################################################################
# DESCRIPTION
    Two fan ducts in close proximity, solved with the multi-body iterative
    solve!. Uses the same geometry as duct.jl (Hill 1978 / Lewis 1991).

    Tests:
    1) BackslashDirichlet solver per body with outer fixed-point iteration
    2) FGSSolver (leaf_size=10000) for comparison
    Both should yield the same RMS flow tangency residual.

# AUTHORSHIP
  * Created   : Apr 2026
  * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
import CSV
import DataFrames: DataFrame
import LinearAlgebra: norm

# ----------------- SIMULATION PARAMETERS --------------------------------------
AOA             = 5                             # (deg) angle of attack
magVinf         = 30.0                          # (m/s) freestream velocity
rho             = 1.225                         # (kg/m^3) air density
Vinf            = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)]

# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
filename        = joinpath(pnl.examples_path, "data", "naca662015.csv")
contour         = CSV.read(filename, DataFrame)

aspectratio     = 0.6                           # Duct trailing edge aspect ratio l/d
d               = 2*0.835                       # (m) duct diameter

# ----------------- SOLVER PARAMETERS ------------------------------------------
NDIVS_theta     = 20                            # Number of azimuthal panels
n_rfl           = 6                             # Chordwise panel control

NDIVS_rfl_up    = [(0.25, n_rfl, 10.0, false),
                   (0.50, n_rfl,  1.0, true),
                   (0.25, n_rfl,  0.1, false)]
NDIVS_rfl_lo    = NDIVS_rfl_up

# ----------------- HELPER: GENERATE A DUCT ------------------------------------
function make_duct(;
        kernel = Union{pnl.ConstantSource, pnl.VortexRing},
        bodytype = pnl.RigidWakeBody
    )
    bodytype = bodytype{kernel}

    xs, ys = pnl.gt.rediscretize_airfoil(contour[:, 1], contour[:, 2],
                                          NDIVS_rfl_up, NDIVS_rfl_lo;
                                          verify_spline=false)
    ys[end] = ys[1]
    xs *= d*aspectratio
    ys *= d*aspectratio
    ys .+= d/2
    points = hcat(xs, ys)

    body = pnl.generate_revolution_liftbody(bodytype, points, NDIVS_theta;
                bodyoptargs = (
                    CPoffset=1e-12,
                    kerneloffset=1e-2,
                    kernelcutoff=1e-14,
                    characteristiclength=(args...)->d*aspectratio,
                ))
    return body
end

# ----------------- GENERATE TWO DUCTS ----------------------------------------
body1 = make_duct()
# body2 = make_duct()
body2 = make_duct(;
            kernel = pnl.ConstantSource,
            bodytype = pnl.NonLiftingBody
        )

# Offset second duct in z by 1.5 diameters
pnl.rotate!(body2, 0, 0, 0; translation=[0.0, 0.0, 1.5*d])

println("Body 1 panels: $(body1.ncells)")
println("Body 2 panels: $(body2.ncells)")

# ----------------- HELPER: SETUP & POST-PROCESS ------------------------------
function setup_bodies!(body, Vinf, magVinf)
    # Reset fields
    pnl.reset!(body)

    # Apply freestream
    pnl.apply_freestream!(body, Vinf)

    # Set wake directions
    if typeof(body) <: pnl.RigidWakeBody
        body.Das[1] .= repeat(Vinf/magVinf, 1, size(body.Das[1], 2))
        body.semiinfinite_wake = false
    end
end

function postprocess!(bodies, Vinf, magVinf, rho, backend)
    pnl.calcfield_U!(bodies; backend)
    pnl.apply_freestream!(bodies, Vinf)
    pnl.calcfield_Cp!(bodies, magVinf; correct_kuttacondition=fill(true, length(bodies)))
    pnl.calcfield_F!(bodies, magVinf, rho)
end

function report_tangency(bodies, label)
    println("\n--- $label ---")
    for (i, body) in enumerate(bodies)
        Udotn = sum(body.velocity .* body.normals, dims=1)
        rms   = sqrt(sum(Udotn .^ 2) / body.ncells)
        maxr  = maximum(abs.(Udotn))
        println("  Body $i: RMS flow tangency = $rms, max = $maxr")
    end
end

backend_direct = pnl.DirectBackend()

# ============================================================================
# TEST 1: BackslashDirichlet
# ============================================================================
println("\n" * "="^60)
println("TEST 1: BackslashDirichlet multi-body solve")
println("="^60)

setup_bodies!(body1, Vinf, magVinf)
setup_bodies!(body2, Vinf, magVinf)
body1.strength .= 0.0
body2.strength .= 0.0

solver1 = pnl.Backslash(body1)
solver2 = pnl.Backslash(body2)

println("\nSolving...")
@time pnl.solve!((body1, body2), (solver1, solver2);
    backend=fill(backend_direct, 2),
    max_outer_iterations=50,
    outer_tolerance=1e-8,
    verbose=true)

println("\nPost-processing...")
postprocess!((body1, body2), Vinf, magVinf, rho, backend_direct)
report_tangency((body1, body2), "BackslashDirichlet")

# Save for comparison
bd_strengths1 = copy(body1.strength)
bd_strengths2 = copy(body2.strength)

# ============================================================================
# TEST 2: FGSSolver (leaf_size=10000)
# ============================================================================
println("\n" * "="^60)
println("TEST 2: FGSSolver multi-body solve (leaf_size=10000)")
println("="^60)

setup_bodies!(body1, Vinf, magVinf)
setup_bodies!(body2, Vinf, magVinf)
body1.strength .= 0.0
body2.strength .= 0.0

println("\nInitializing FGSSolvers (one per body)...")
@time fgs1 = pnl.FGSSolver(body1;
    leaf_size=10000,
    expansion_order=10,
    multipole_acceptance=0.4,
    max_iterations=500,
    inner_iterations=20,
    tolerance=1e-6,
    rlx=1.0,
    shrink=true,
    reverse_pass=false,
    verbose=false
)
@time fgs2 = pnl.FGSSolver(body2;
    leaf_size=10000,
    expansion_order=10,
    multipole_acceptance=0.4,
    max_iterations=500,
    inner_iterations=20,
    tolerance=1e-6,
    rlx=1.0,
    shrink=true,
    reverse_pass=false,
    verbose=false
)

backend_fmm = pnl.FastMultipoleBackend(
    expansion_order=10,
    multipole_acceptance=0.4,
    leaf_size=10000
)

println("\nSolving...")
@time pnl.solve!(body2, fgs2;
    backend=backend_direct)

pnl.reset!(body2)
pnl.apply_freestream!(body2, Vinf)
pnl.calcfield_U!(body2; backend=backend_direct)
udotn = sum(body2.velocity .* body2.normals, dims=1)
rms = sqrt(sum(udotn .^ 2) / body2.ncells)
maxr = maximum(abs.(udotn))
println("  Body 2 (FGS only): RMS flow tangency = $rms, max = $maxr")

   
@time pnl.solve!((body1, body2), (fgs1, fgs2);
    backend=fill(backend_direct, 2),
    max_outer_iterations=50,
    outer_tolerance=1e-8,
    verbose=true)

println("\nPost-processing...")
postprocess!((body1, body2), Vinf, magVinf, rho, backend_direct)
report_tangency((body1, body2), "FGSSolver (per-body)")

# Save for comparison
fgs_strengths1 = copy(body1.strength)
fgs_strengths2 = copy(body2.strength)

backend_fmm = pnl.FastMultipoleBackend(
    expansion_order=10,
    multipole_acceptance=0.4,
    leaf_size=10000
)

# # ============================================================================
# # TEST 3: Single coupled FGSSolver (both bodies in one solver)
# # ============================================================================
# println("\n" * "="^60)
# println("TEST 3: Single coupled FGSSolver (leaf_size=10000)")
# println("="^60)

# setup_bodies!(body1, Vinf, magVinf)
# setup_bodies!(body2, Vinf, magVinf)
# body1.strength .= 0.0
# body2.strength .= 0.0

# println("\nInitializing coupled FGSSolver...")
# @time fgs_coupled = pnl.FGSSolver((body1, body2);
#     leaf_size=10000,
#     expansion_order=10,
#     multipole_acceptance=0.4,
#     max_iterations=500,
#     inner_iterations=20,
#     tolerance=1e-8,
#     rlx=1.0,
#     shrink=true,
#     reverse_pass=false,
#     verbose=false
# )

# println("\nSolving...")
# @time pnl.solve!((body1, body2), fgs_coupled; backend=backend_fmm)

# println("\nPost-processing...")
# postprocess!((body1, body2), Vinf, magVinf, rho, backend_fmm)
# report_tangency((body1, body2), "FGSSolver (coupled)")

# ============================================================================
# COMPARISON
# ============================================================================
println("\n" * "="^60)
println("COMPARISON")
println("="^60)
println("BackslashDirichlet vs FGSSolver (per-body):")
println("  Max strength diff (body1): ", maximum(abs.(fgs_strengths1 .- bd_strengths1)))
println("  Max strength diff (body2): ", maximum(abs.(fgs_strengths2 .- bd_strengths2)))
# println("BackslashDirichlet vs FGSSolver (coupled):")
# println("  Max strength diff (body1): ", maximum(abs.(body1.strength .- bd_strengths1)))
# println("  Max strength diff (body2): ", maximum(abs.(body2.strength .- bd_strengths2)))
