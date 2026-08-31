#=##############################################################################
# DESCRIPTION
    Fan duct replicating the experiment reported by V. P. Hill (1978), "A
    Surface Vorticity Theory for Propeller Ducts and Turbofan Engine Cowls in
    Non-Axisymmetric Incompressible Flow." The same experiment is also
    discussed in Sections 4.4 and 6.3.1 of Lewis, R. (1991), "Vortex Element
    Methods for Fluid Dynamic Analysis of Engineering Systems."

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jan 2023
  * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
import GeometricTools as gt
include(joinpath(pnl.examples_path, "helper_functions.jl"))
import CSV
import DataFrames: DataFrame
using LinearAlgebra

include(joinpath(pnl.examples_path, "duct_postprocessing.jl"))

run_name        = "duct-hill00"             # Name of this run

save_path       = joinpath("data", run_name) # Where to save outputs
fluiddomain     = false                     # Whether to generate fluid domain
paraview        = true                      # Whether to visualize with Paraview
call_paraview   = false                     # Whether to call Paraview at the end

save_plots      = false                     # Whether to save plots or not
# Where to save plots (default to re-generating the figures that are used
# in the docs)
fig_path = joinpath(pnl.examples_path, "..", "docs", "resources", "images")


# ----------------- SIMULATION PARAMETERS --------------------------------------
AOAs            = [0, 5, 15]                # (deg) angles of attack to evaluate
# AOAs            = [5]                # (deg) angles of attack to evaluate
magVinf         = 30.0                      # (m/s) freestream velocity
rho             = 1.225                     # (kg/m^3) air density


# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
# Read duct contour (table in figure 7.4 of Lewis 1991)
filename        = joinpath(pnl.examples_path, "data", "naca662015.csv")
contour         = CSV.read(filename, DataFrame)

aspectratio     = 0.6                       # Duct trailing edge aspect ratio l/d
d               = 2*0.835                   # (m) duct diameter


# ----------------- SOLVER PARAMETERS ------------------------------------------
# Discretization
NDIVS_theta     = 21                        # Number of azimuthal panels

# NOTE: NDIVS is the number of divisions (panels) in each dimension. This can be
#       either an integer, or an array of tuples as shown below

n_rfl           = 6                        # This controls the number of chordwise panels

NDIVS_rfl_up = [                            # Discretization of airfoil upper surface
            # 0 to 0.25 of the airfoil has `n_rfl` panels at a geometric expansion of 10 that is not central
                # (0.25, n_rfl >> 1,   1.0, false),
                (0.25, n_rfl,   10.0, false),
            # 0.25 to 0.75 of the airfoil has `n_rfl` panels evenly spaced
                (0.50, n_rfl,    1.0, true),
            # 0.75 to 1.00 of the airfoil has `n_rfl` panels at a geometric expansion of 0.1 that is not central
                # (0.25, n_rfl >> 1,    1.0, false)]
                (0.25, n_rfl,    0.1, false)]

NDIVS_rfl_lo = NDIVS_rfl_up                 # Discretization of airfoil lower surface

# NOTE: A geometric expansion of 10 that is not central means that the last
#       panel is 10 times larger than the first panel. If central, the
#       middle panel is 10 times larger than the peripheral panels.

# Solver: Vortex-ring least-squares
# bodytype        = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64} # Elements and wake model
# bodytype        = pnl.RigidWakeBody{pnl.ConstantDoublet, 1, Float64} # Elements and wake model
# kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}
kernel = Union{pnl.ConstantSource, pnl.VortexRing}
# kernel = pnl.VortexRing
bodytype = pnl.RigidWakeBody{kernel} # Elements and wake model


# ----------------- GENERATE BODY ----------------------------------------------
# Re-discretize the contour of the body of revolution according to NDIVS
xs, ys = gt.rediscretize_airfoil(contour[:, 1], contour[:, 2],
                                        NDIVS_rfl_up, NDIVS_rfl_lo;
                                        verify_spline=false)

# Make sure that the contour is closed
ys[end] = ys[1]
ys .*= 1.0

# Scale contour by duct length
xs *= d*aspectratio
ys *= d*aspectratio

# Move contour to the radial position
ys .+= d/2

# Collect points that make the contour of the body of revolution
points = hcat(xs, ys)

# Generate body of revolution
generate_body() = generate_revolution_liftbody(bodytype, points, NDIVS_theta;
                                        bodyoptargs = (
                                                        cp_outer=true,
                                                        core_size=1e-10,
                                                        kernelcutoff=1e-14,
                                                        characteristiclength=(args...)->d*aspectratio,
                                                        semiinfinite_wake=false
                                            )
                                        )
body_fgs = generate_body()
body_bsd = generate_body()

# body.shedding[1] = body.shedding[1][:,1:1]
# body.nsheddings = 1
# body.shedding_full .= -1
# idx_shed = body.shedding[1, 1]
# body.shedding_full[1:2, idx_shed] .= view(body.shedding[1], 2:3, 1)
# body.shedding_full[3:4, idx_shed] .= 1
# idx_shed = body.shedding[1][4, 1]
# body.shedding_full[1:2, idx_shed] .= view(body.shedding[1], 5:6, 1)
# body.shedding_full[3:4, idx_shed] .= 1

println("Number of panels:\t$(body_fgs.ncells)")

vtks = save_path*"/"                        # String with VTK output files


# ----------------- CALL SOLVER ------------------------------------------------
global solver = nothing
# for (i, AOA) in enumerate(AOAs)             # Sweep over angle of attack
i = 2
AOA = AOAs[i]
    
# println("Solving body...")

# Freestream vector
Vinf = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)]

# Freestream at every control point
Uinfs = repeat(Vinf, 1, body_fgs.ncells)
leaf_size = 100
expansion_order = 10
multipole_acceptance = 0.4
backend = pnl.FastMultipoleBackend(;
                                expansion_order,
                                multipole_acceptance,
                                leaf_size=20
                            )

dt_val = 0.00998003992015968
eta = 0.3
body_bsd.Das[1] .= repeat(Vinf * dt_val * eta, 1, size(body_bsd.Das[1], 2))
body_fgs.Das[1] .= repeat(Vinf * dt_val * eta, 1, size(body_fgs.Das[1], 2))

body_bsd.velocity .= repeat(Vinf, 1, body_bsd.ncells)
body_fgs.velocity .= repeat(Vinf, 1, body_fgs.ncells)

println("Initializaing...")
solver_fgs = pnl.FGSSolver(body_fgs;
    max_iterations=500,         # Maximum number of iterations
    tolerance=1.0e-6,            # Convergence tolerance
    rlx=1.0,                  # Relaxation factor
    expansion_order,
    multipole_acceptance,
    leaf_size=100,
    shrink=true,
    recenter=false,
    inner_iterations=20,
    reverse_pass=false,
    verbose=false
)
solver_bsd = pnl.Backslash(body_bsd)

println("solving...")
pnl.solve!(body_bsd, solver_bsd; backend)
pnl.solve!(body_fgs, solver_fgs; backend)

err = norm(body_bsd.strength - body_fgs.strength) / norm(body_bsd.strength)
println("Relative error: ", err)
