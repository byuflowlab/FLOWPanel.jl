envpath = joinpath(pwd(), "environment-testfmm")

#= -------------------- CREATE THE ENVIRONMENT -----------------------------------
---------------------------------------------------------------------------------- =#
if isdir(envpath)
    rm(envpath, recursive=true)
end
mkpath(envpath)

using Pkg
Pkg.activate(envpath)

Pkg.add(url="https://github.com/byuflowlab/FastMultipole")
Pkg.add(url="https://github.com/byuflowlab/FLOWPanel.jl", rev="fastmultipole-mvp")

Pkg.add("CSV")
Pkg.add("DataFrames")

Pkg.status()


# -------------------- IMPORT MODULES --------------------------------------------
import FLOWPanel as pnl
import FLOWPanel: dot, norm

# Check if PythonPlot is available
try
    import FLOWPanel: plt, @L_str
    global pyplot = true
catch
    global pyplot = false
end



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

run_name        = "centerbody"      # Name of this run

save_path       = "centerbody-lewis-baseline00/"   # Where to save outputs
paraview        = false                      # Whether to visualize with Paraview


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
bodytype        = pnl.RigidWakeBody{pnl.VortexRing}    # Elements and wake model


# ----------------- GENERATE BODY ----------------------------------------------
# Generate grid of body of revolution
# holeradius      = 0.10*R                    # Hole in centerbody
                                            # Points of contour to revolve
points = Matrix(contour_lewis[2:end-1, :])
# points[1, 2] += holeradius

# Force body to have no wake
shedding = zeros(Int, 6, 0)

body = pnl.generate_revolution_liftbody(bodytype, points, NDIVS_theta;
                                        # Loop the azimuthal dimension to close the surface
                                        loop_dim=2,
                                        # Rotate the axis of rotation to align with x-axis
                                        axis_angle=90,
                                        # Indicate that this body is open at the trailing edge
                                        closed_contour=false,
                                        overwrite_shedding=shedding
                                        )

println("Number of panels:\t$(body.ncells)")


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

# Calculate surface velocity on the body
@time Us = pnl.calcfield_U(body, body)

Us_old = deepcopy(Us)
strength_old = deepcopy(body.strength)

# Calculate surface velocity U_∇μ due to the gradient of the doublet strength
UDeltaGamma = pnl.calcfield_Ugradmu(body)

# Add both velocities together
pnl.addfields(body, "Ugradmu", "U")

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

# display(Us_old)
# display(strength_old)

if pyplot

    # ----------------- COMPARISON TO EXPERIMENTAL DATA ------------------------
    #=
        NOTE: Here we take a slice of the body and plot the velocity distribution
        along the slice.
    =#

    # Get a slice of the body
    local position  = 0.0        # Position of slice (slice along origin)
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

    Usurfs_old = pnl.norm.(sliceUs)/magVinf

    # Plot contour of centerbody
    ax2 = ax.twinx()
    xs = vcat(slicepoints[1, :], reverse(slicepoints[1, :]), slicepoints[1, 1])
    ys = vcat(slicepoints[3, :], -reverse(slicepoints[3, :]), slicepoints[3, 1])
    ax2.plot(xs, ys, "-k", alpha=0.25)

    ax2.set_aspect("equal")

    # Beautify the plot
    xlims = [-0.1, 1.5]
    ylims = [-0.5, 1.5]
    ax.set_xlim(xlims)
    ax.set_xticks(0:0.5:1.5)
    ax.set_ylim(ylims)
    ax.set_yticks(0:0.5:1.5)

    ax2.set_aspect(1.0)
    ax2.set_ylim([-0.62, 1])
    ax2.set_yticks([])
    ax2.plot(xlims, zeros(2), ":k", alpha=0.25, linewidth=1)

    ax.set_xlabel(L"$x$-position (m)")
    ax.set_ylabel(L"Velocity $u/u_\infty$")
    ax.legend(loc="best", frameon=false, fontsize=8)

    for a in [ax, ax2]
        a.spines["right"].set_visible(false)
        a.spines["top"].set_visible(false)
    end

    fig.tight_layout()

end

# Convert body to FMM body
body_fmm = pnl.body2fmmbody(body)

# Set freestream at each FMM panel
for (i, Uinf) in zip(eachindex(body_fmm.velocity), eachcol(Uinfs))
    body_fmm.velocity[i] -= body_fmm.velocity[i]
    body_fmm.velocity[i] += Uinf
end

# Evaluate mutually-induced velocity of panels for benchmark
@time tree_unstructured = pnl.FMM.FastMultipole.fmm!(body_fmm; expansion_order=14, leaf_size=50, multipole_threshold=0.4,
                                                    unsort_bodies=true)

# Save FMM body for debugging purposes
if !isnothing(save_path)

    pnl.gt.create_path(save_path, false)
    files = pnl.FMM.vtk(run_name, body_fmm; save_path=save_path)

    if paraview
        str = files[1]
        for fname in files[2:end]
            str *= "; "*splitdir(fname)[end]
        end

        run(`paraview --data=$(str)`)
    end
end

Us_fmm = [V[i] for i in 1:3, V in body_fmm.velocity]

println("Maximum error:\t$(maximum(abs.(Us_fmm .- Us_old)))")
println("Maximum error:\t$(maximum(abs.(Us_fmm .- Us_old) ./ abs.(Us_old)) * 100) %")

# Convert body to FMM body
body_fmm = pnl.body2fmmbody(body)

# Set freestream at each FMM panel
for (i, Uinf) in zip(eachindex(body_fmm.velocity), eachcol(Uinfs))
    body_fmm.velocity[i] -= body_fmm.velocity[i]
    body_fmm.velocity[i] += Uinf
end

# Evaluate mutually-induced velocity of panels for benchmark
@time tree_unstructured = pnl.FMM.FastMultipole.direct!(body_fmm)

# Save FMM body for debugging purposes
if !isnothing(save_path)

    pnl.gt.create_path(save_path, false)
    files = pnl.FMM.vtk(run_name, body_fmm; save_path=save_path)

    if paraview
        str = files[1]
        for fname in files[2:end]
            str *= "; "*splitdir(fname)[end]
        end

        run(`paraview --data=$(str)`)
    end
end

Us_fmm = [V[i] for i in 1:3, V in body_fmm.velocity]
# display(Us_old)

println("Maximum error:\t$(maximum(abs.(Us_fmm .- Us_old)))")
println("Maximum error:\t$(maximum(abs.(Us_fmm .- Us_old) ./ abs.(Us_old)) * 100) %")
