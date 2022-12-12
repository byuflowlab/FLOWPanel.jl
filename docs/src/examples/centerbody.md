# Centerbody

```@raw html
<center>
  <img src="../../assets/images/centerbody-viz00.png" alt="Pic here" style="width: 100%;"/>
</center>
```

In this example we solve the flow around a body of revolution resembling
the centerbody (hub) of a ducted fan.

## Source Elements
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
# holeradius      = 0.10*R                    # Whole in centerbody
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

## Slice

FLOWPanel provides the following function to obtain the solution field
along a slice along a body:

```@docs
FLOWPanel.slicefield
```

Now we process the solution to plot the surface velocity along a slice
of the body of revolution.


```julia
import PyPlot as plt
import LaTeXStrings: @L_str

# ----------------- COMPARISON TO EXPERIMENTAL DATA ----------------------------
#=
    NOTE: Here we take a slice of the body and plot the velocity distribution
    along the slice.
=#

# Get a slice of the body
position        = 0.0        # Position of slice (slice along origin)
direction       = [0, 1, 0]  # Direction of slice (slice along the xz-plane)
row             = false      # If true, it slices along azimuth; centerline if false

slicepoints, sliceCps = pnl.slicefield(body, "Cp", position, direction, row)
slicepoints, sliceUs = pnl.slicefield(body, "U", position, direction, row)

# Plot experimental surface velocity distribution (figure 4.6 in Lewis 1991)
fig = plt.figure(figsize=[7, 5*0.8]*2/3)
ax = fig.gca()

filename = joinpath(pnl.examples_path, "data",
                                "centerbody-lewis-fig4p6.csv")
VoVinf_lewis = CSV.read(filename, DataFrame)

ax.plot(VoVinf_lewis[:, 1], VoVinf_lewis[:, 2], "ok",
                            markersize=5, label="Experimental")

# Plot surface velocity distribution of FLOWPanel
ax.plot(slicepoints[1, :], pnl.norm.(sliceUs)/magVinf, "-", color="cyan",
                            linewidth=2.0, alpha=0.9, label="FLOWPanel")

# Plot contour of centerbody
ax2 = ax.twinx()
xs = vcat(slicepoints[1, :], reverse(slicepoints[1, :]), slicepoints[1, 1])
ys = vcat(slicepoints[3, :], -reverse(slicepoints[3, :]), slicepoints[3, 1])
ax2.plot(xs, ys, "-k", alpha=0.25)

```
```@raw html
<center>
    <br><b>Surface velocity</b><br>
    <img src="../../assets/images/centerbody-lewis00-velocity-source.png" alt="Pic here" style="width: 60%;"/>
</center>
```


## Vortex Ring Elements

While source elements are physically adequate to model a non-lifting body,
in some circumstances it may be benefitial to use all vortex ring elements.
A thick body with only vortex ring elements leads to a surface velocity
that is inaccurate at the exact surface of the body, but that
approximates the physical solution away from the surface. For this
reason, we probe the velocity used to calculate Cp slightly away from
the body.

Here we repeat the example using only vortex ring elements.


```julia
bodytype = pnl.RigidWakeBody{pnl.VortexRing}    # Elements and wake model

body = bodytype(trigrid;            # Body to be solved
                CPoffset=1e-14,     # Offset control points slightly in the positive normal direction
                characteristiclength=(args...)->R,   # Characheristic length for control point offset
                kerneloffset=1e-8,  # Offset of kernel to avoid singularities
                kernelcutoff=1e-14  # Cutoff of kernel to avoid singularities
                )


# ----------------- CALL SOLVER ------------------------------------------------
println("Solving body...")

# Freestream at every control point
Uinfs = repeat(Vinf, 1, body.ncells)

# Unitary direction of semi-infinite vortex at points `a` and `b` of each
# trailing edge panel
# NOTE: In this case they are empty arrays since there is no wake
Das = repeat(Vinf/magVinf, 1, body.nsheddings)
Dbs = repeat(Vinf/magVinf, 1, body.nsheddings)

# Solve body (panel strengths) giving `Uinfs` as boundary conditions and
# `Das` and `Dbs` as trailing edge rigid wake direction
@time pnl.solve(body, Uinfs, Das, Dbs)


# ----------------- POST PROCESSING --------------------------------------------
println("Post processing...")
# NOTE: A thick body with only vortex ring elements leads to a surface velocity
#       that is inaccurate at the exact surface of the body, but that
#       approximates the physical solution away from the surface. For this
#       reason, we probe the velocity used to calculate Cp slightly away from
#       the body

# Calculate surface velocity on the body
Us = pnl.calcfield_U(body, body,
                            offset=0.08, characteristiclength=(args...)->R)

# Calculate pressure coefficient
Cps = pnl.calcfield_Cp(body, magVinf)

# Calculate the force of each panel
Fs = pnl.calcfield_F(body, magVinf, rho)


```
```@raw html
<center>
    <br><b>Surface velocity</b><br>
    <img src="../../assets/images/centerbody-lewis00-velocity-vortexring.png" alt="Pic here" style="width: 60%;"/>
</center>
```

