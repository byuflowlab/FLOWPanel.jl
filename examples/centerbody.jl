#=##############################################################################
# DESCRIPTION
    Centerbody (body of revolution) replicating the experiment reported in
    Section Section 4.3.2 of Lewis, R. (1991), "Vortex Element Methods for
    Fluid Dynamic Analysis of Engineering Systems."

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Dec 2022
  * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
import CSV
import DataFrames: DataFrame

run_name        = "centerbody-lewis00"      # Name of this run

save_path       = ""                        # Where to save outputs
paraview        = true                      # Whether to visualize with Paraview


# ----------------- SIMULATION PARAMETERS --------------------------------------
AOA             = 0.0                       # (deg) angle of attack
magVinf         = 30.0                      # (m/s) freestream velocity
Vinf            = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)] # Freestream

rho             = 1.225                     # (kg/m^3) air density


# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
# Read body contour (table 4.2 in Lewis 1991)
filename        = joinpath(pnl.examples_path, "data",
                                "centerbody-lewis-table4p2.csv")
contour_lewis   = CSV.read(filename, DataFrame)

R               = maximum(contour_lewis[:, 2]) # (m) max radius

# ----------------- SOLVER PARAMETERS ------------------------------------------
# Discretization
NDIVS_theta     = 60                        # Number of azimuthal panels

# Solver
# kerneltype = pnl.ConstantDoublet          # Kernel type to use
# kerneltype = pnl.VortexRing               # Kernel type to use
kerneltype = pnl.ConstantSource               # Kernel type to use
# kerneltype = Union{pnl.ConstantSource, pnl.ConstantDoublet}               # Kernel type to use
bodytype        = pnl.NonLiftingBody{kerneltype}    # Elements and wake model
# bodytype2        = pnl.NonLiftingBody{pnl.ConstantDoublet}    # Elements and wake model


# ----------------- GENERATE BODY ----------------------------------------------
# Generate grid of body of revolution
# holeradius      = 0.10*R                    # Hole in centerbody
                                            # Points of contour to revolve
points = Matrix(contour_lewis[2:end-1, :])
# points[1, 2] += holeradius

grid = pnl.gt.surface_revolution(points, NDIVS_theta;
                                    # Loop the azimuthal dimension to close the surface
                                    loop_dim=2,
                                    # Rotate the axis of rotation to align with x-axis
                                    axis_angle=90
                                )

# Rotate the body of revolution to align centerline with x-axis
Oaxis = pnl.gt.rotation_matrix2(0, 0, 90)          # Rotation matrix
O = zeros(3)                                       # Translation of coordinate system
pnl.gt.lintransform!(grid, Oaxis, O)

# Triangular grid (splits quadrangular panels into triangular panels)
split_dim = 1                               # Dimension to split into triangles
trigrid = pnl.gt.GridTriangleSurface(grid, split_dim)

# Generate body to be solved
body = bodytype(trigrid)
# body2 = bodytype2(trigrid)

println("Number of panels:\t$(body.ncells)")


# ----------------- CALL SOLVER ------------------------------------------------
println("Solving body...")

# Freestream at every control point
Uinfs = repeat(Vinf, 1, body.ncells)

#--- debugging ---#

G = zeros(body.ncells, body.ncells)
normals = pnl._calc_normals(body)
CPs = pnl._calc_controlpoints(body, normals; off=1e-14)
backend = pnl.FastMultipoleBackend(leaf_size=1000000)
pnl._G_phi!(body, kerneltype, G, CPs, backend; kerneloffset=1.0e-8)

body.strength .= 1.0
test_rhs = zeros(body.ncells)
pnl._phi!(body, CPs, test_rhs, backend; kerneloffset=1.0e-8)
debug_rhs = G * body.strength[:,1]
@show maximum(abs.(test_rhs .- debug_rhs))



# Solve body (panel strengths) giving `Uinfs` as boundary conditions
# @time pnl.solve(body, Uinfs)
# strengths1 = deepcopy(body.strength)
@time begin
    # solver = pnl.Backslash(body; least_squares=false)
    # solver = pnl.KrylovSolver(body;
    #     method=:gmres,
    #     itmax=20,
    #     atol=1e-4,
    #     rtol=1e-4,
        backend=pnl.FastMultipoleBackend(
                                    expansion_order=7,
                                    multipole_acceptance=0.4,
                                    leaf_size=10
                                )
    # )

    # G, rhs = pnl.solve(body, Uinfs)
    # G2, rhs2 = pnl.solve(body2, Uinfs)

    solver = pnl.Backslash(body; least_squares=false)
    pnl.solve2!(body, Uinfs, solver; backend)
end
strengths2 = deepcopy(body.strength)

# ----------------- POST PROCESSING --------------------------------------------
println("Post processing...")

# Calculate surface velocity on the body with the FMM backend
# backend = pnl.FastMultipoleBackend(
#                                     expansion_order=7,
#                                     multipole_acceptance=0.4,
#                                     leaf_size=10
#                                 )
println("Running FMM:")
@time Us = pnl.calcfield_U(body, body; backend)
Us_fmm = deepcopy(Us) # save for comparison/verification
# str = pnl.save(body, "debug_fmm"; path="./", debug=true, backend)

# uncomment to profile
# using Profile
# using PProf
# @profile Us = pnl.calcfield_U(body, body; backend)
# Profile.clear()
# @profile Us = pnl.calcfield_U(body, body; backend)
# pprof()

Us .= zero(eltype(Us))
backend = pnl.DirectBackend()
println("Running Direct:")
@time Us = pnl.calcfield_U(body, body; backend)
Us_direct = deepcopy(Us)
# str = pnl.save(body, "debug_direct"; path="./", debug=true, backend)

# calculate error
max_err = maximum(abs.(Us_direct - Us_fmm))
println("Max velocity error between Direct and FMM normalized by Vinf: $(max_err/magVinf)")

if kerneltype == pnl.ConstantDoublet || kerneltype == pnl.VortexRing || kerneltype == Union{pnl.ConstantSource, pnl.ConstantDoublet}
    # Calculate velocity induced by doublet strength gradient
    if kerneltype == Union{pnl.ConstantSource, pnl.ConstantDoublet}
        Gammai = 2
    else
        Gammai = 1
    end
    @time UDeltaGamma = pnl.calcfield_Ugradmu(body; sharpTE=true, force_cellTE=false, Gammai)

    # Add both velocities together
    pnl.addfields(body, "Ugradmu", "U")
end

Us_tot = pnl.get_field(body, "U")["field_data"] # save for comparison/verification
Us_tot = [Us_tot[j][i] for i=1:3, j=1:size(Us_tot,1)]

# check normal flow condition
normals = pnl._calc_normals(body)
Udotn = sum(Us_tot .* normals, dims=1)
resid = maximum(abs.(Udotn))
println("Max flow tangency residual: $resid")

# Calculate pressure coefficient
@time Cps = pnl.calcfield_Cp(body, magVinf)

# Calculate the force of each panel
@time Fs = pnl.calcfield_F(body, magVinf, rho)


# ----------------- VISUALIZATION ----------------------------------------------
if paraview
    str = save_path*"/"

    # Save body as a VTK
    str *= pnl.save(body, run_name; path=save_path, debug=true)

    # Call Paraview
    run(`paraview --data=$(str)`)
end


# ----------------- COMPARISON TO EXPERIMENTAL DATA ----------------------------
include(joinpath(pnl.examples_path, "centerbody_postprocessing.jl"))

save_outputs = !true                        # Whether to save outputs or not

# Where to save figures (default to re-generating the figures that are used
# in the docs)
fig_path = joinpath(pnl.examples_path, "..", "docs", "resources", "images")
outdata_path = joinpath(pnl.examples_path, "..", "docs", "resources", "data")

if save_outputs
    fig.savefig(joinpath(fig_path, "$(run_name)-velocity-source.png"),
                                                dpi=300, transparent=true)
end


# ----------------- VORTEX RING SOLVER -----------------------------------------
# include(joinpath(pnl.examples_path, "centerbody_vortexring.jl"))
