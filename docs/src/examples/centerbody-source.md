```@raw html
<center>
  <img src="../../assets/images/centerbody-viz00.png" alt="Pic here" style="width: 100%;"/>
</center>
```

In this example we solve the flow around a body of revolution resembling
the centerbody (hub) of a ducted fan.

# Source Elements
First we run this example with source elements, which are especially
accurate for non-lifting bodies.


```julia
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
bodytype        = pnl.NonLiftingBody{pnl.ConstantSource}    # Elements and wake model


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

println("Number of panels:\t$(body.ncells)")


# ----------------- CALL SOLVER ------------------------------------------------
println("Solving body...")

# Freestream at every control point
Uinfs = repeat(Vinf, 1, body.ncells)

# Solve body (panel strengths) giving `Uinfs` as boundary conditions
@time pnl.solve(body, Uinfs)


# ----------------- POST PROCESSING --------------------------------------------
println("Post processing...")

# Calculate surface velocity on the body
@time Us = pnl.calcfield_U(body, body)

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


```
(see the complete example under
[examples/centerbody.jl](https://github.com/byuflowlab/FLOWPanel.jl/blob/master/examples/centerbody.jl)
)

