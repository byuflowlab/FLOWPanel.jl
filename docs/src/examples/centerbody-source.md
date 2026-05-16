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
include(joinpath(pnl.examples_path, "helper_functions.jl"))
import GeometricTools as gt
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

body = generate_revolution(bodytype, points, NDIVS_theta;
                            # Loop the azimuthal dimension to close the surface
                            loop_dim=2,
                            # Rotate the axis of rotation to align with x-axis
                            axis_angle=90,
                            # Split quadrangular panels into triangles along dimension 1
                            dimsplit=1,
                            # Rotate centerline to align with x-axis
                            gridprocessing = grid -> gt.lintransform!(grid,
                                                gt.rotation_matrix2(0, 0, 90), zeros(3))
                          )

println("Number of panels:\t$(body.ncells)")


# ----------------- CALL SOLVER ------------------------------------------------
println("Solving body...")

@time begin
    solver = pnl.FGSSolver(body;
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
    backend = pnl.FastMultipoleBackend(
        expansion_order=7,
        multipole_acceptance=0.4,
        leaf_size=10,
    )

    pnl.solve!(body, solver; backend)
end


# ----------------- POST PROCESSING --------------------------------------------
println("Post processing...")

# Calculate surface velocity on the body
@time pnl.calcfield_U!(body; backend)

# Calculate gauge pressure
eachcol(body.velocity) .+= Ref(Vinf)
@time pnl.calcfield_P!(body, magVinf, rho)

# Calculate the force of each panel
@time pnl.calcfield_F!(body)


# ----------------- VISUALIZATION ----------------------------------------------
if paraview
    pnl.write_vtk(joinpath(save_path, "centerbody"), body, 0, 0)
end


```
(see the complete example under
[examples/centerbody.jl](https://github.com/byuflowlab/FLOWPanel.jl/blob/master/examples/centerbody.jl)
)

