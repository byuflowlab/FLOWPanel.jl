#=
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
=#

# -------------------- IMPORT MODULES --------------------------------------------
import FLOWPanel as pnl
import FLOWPanel: dot, norm
using Test
using StaticArrays

# Check if PyPlot is available
# try
#     import FLOWPanel: plt, @L_str
#     global pyplot = true
# catch
    global pyplot = false
# end



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
# bodytype        = pnl.RigidWakeBody{pnl.VortexRing}    # Elements and wake model
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
# if paraview
#     str = save_path*"/"
#
#     # Save body as a VTK
#     str *= pnl.save(body, run_name; path=save_path, debug=true)
#
#     # Call Paraview
#     run(`paraview --data=$(str)`)
# end

# display(Us_old)
# display(strength_old)

# if pyplot
#
#     # ----------------- COMPARISON TO EXPERIMENTAL DATA ------------------------
#     #=
#         NOTE: Here we take a slice of the body and plot the velocity distribution
#         along the slice.
#     =#
#
#     # Get a slice of the body
#     local position  = 0.0        # Position of slice (slice along origin)
#     direction       = [0, 1, 0]  # Direction of slice (slice along the xz-plane)
#     row             = false      # If true, it slices along azimuth; centerline if false
#
#     slicepoints, sliceCps = pnl.slicefield(body, "Cp", position, direction, row)
#     slicepoints, sliceUs = pnl.slicefield(body, "U", position, direction, row)
#
#     # Plot experimental surface velocity distribution (figure 4.6 in Lewis 1991)
#     fig = plt.figure(figsize=[7, 5*0.8]*2/3)
#     ax = fig.gca()
#
#     filename = joinpath(pnl.examples_path, "data",
#                                     "centerbody-lewis-fig4p6.csv")
#     VoVinf_lewis = CSV.read(filename, DataFrame)
#
#     ax.plot(VoVinf_lewis[:, 1], VoVinf_lewis[:, 2], "ok",
#                                 markersize=5, label="Experimental")
#
#     # Plot surface velocity distribution of FLOWPanel
#     ax.plot(slicepoints[1, :], pnl.norm.(sliceUs)/magVinf, "-", color="cyan",
#                                 linewidth=2.0, alpha=0.9, label="FLOWPanel")
#
#     Usurfs_old = pnl.norm.(sliceUs)/magVinf
#
#     # Plot contour of centerbody
#     ax2 = ax.twinx()
#     xs = vcat(slicepoints[1, :], reverse(slicepoints[1, :]), slicepoints[1, 1])
#     ys = vcat(slicepoints[3, :], -reverse(slicepoints[3, :]), slicepoints[3, 1])
#     ax2.plot(xs, ys, "-k", alpha=0.25)
#
#     ax2.set_aspect("equal")
#
#     # Beautify the plot
#     xlims = [-0.1, 1.5]
#     ylims = [-0.5, 1.5]
#     ax.set_xlim(xlims)
#     ax.set_xticks(0:0.5:1.5)
#     ax.set_ylim(ylims)
#     ax.set_yticks(0:0.5:1.5)
#
#     ax2.set_aspect(1.0)
#     ax2.set_ylim([-0.62, 1])
#     ax2.set_yticks([])
#     ax2.plot(xlims, zeros(2), ":k", alpha=0.25, linewidth=1)
#
#     ax.set_xlabel(L"$x$-position (m)")
#     ax.set_ylabel(L"Velocity $u/u_\infty$")
#     ax.legend(loc="best", frameon=false, fontsize=8)
#
#     for a in [ax, ax2]
#         a.spines["right"].set_visible(false)
#         a.spines["top"].set_visible(false)
#     end
#
#     fig.tight_layout()
#
# end

#--- use FastMultipole to compute induced velocity ---#


println("\n--- INDUCED VELOCITY: FMM ---\n")

# Convert body to FMM body
body_fmm = pnl.body2fmmbody(body)

# Set freestream at each FMM panel
for (i, Uinf) in zip(eachindex(body_fmm.velocity), eachcol(Uinfs))
    body_fmm.velocity[i] -= body_fmm.velocity[i]
    body_fmm.velocity[i] += Uinf
end

# Evaluate mutually-induced velocity of panels for benchmark
expansion_order, leaf_size_source, multipole_threshold = 6, 20, 0.4
# @time pnl.FMM.FastMultipole.fmm!(body_fmm; expansion_order, leaf_size_source, multipole_threshold)

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

#--- direct calcs just in case ---#

println("\n--- INDUCED VELOCITY: DIRECT ---\n")

# Convert body to FMM body
body_fmm = pnl.body2fmmbody(body)

# Set freestream at each FMM panel
for (i, Uinf) in zip(eachindex(body_fmm.velocity), eachcol(Uinfs))
    body_fmm.velocity[i] -= body_fmm.velocity[i]
    body_fmm.velocity[i] += Uinf
end

# Evaluate mutually-induced velocity of panels for benchmark
# @time tree_unstructured = pnl.FMM.FastMultipole.direct!(body_fmm)

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


function save_strengths(body)
    strengths = zeros(length(body.strengths))
    for i in eachindex(strengths)
        strengths[i] = body.strengths[i][1]
    end
    return strengths
end

#--- test direct solver ---#
println("\n--- SOLVER: ---\n")

function convert_strength(strength, TK::Type{pnl.FMM.ConstantNormalDoublet}, new_kernel::Type{pnl.FMM.UniformSourceNormalDoublet})
    return SVector{2,eltype(strength)}(0.0, strength[1])
end

function convert_strength(strength, TK::Type{pnl.FMM.UniformSource}, new_kernel::Type{pnl.FMM.UniformSourceNormalDoublet})
    return SVector{2,eltype(strength)}(strength[1], 0.0)
end

function change_kernel(panel_array::pnl.FMM.UnstructuredGrid{TK,TF,NK,NS}, new_kernel) where {TK,TF,NK,NS}
    (; points, meshcells, control_points, normals, strengths, potential, velocity, sigma, panels, wake_points) = panel_array
    new_NK = pnl.FMM.kernel_multiplicity(new_kernel)
    new_panel_type = pnl.FMM.Panel{TF,new_NK,NS}
    new_panels = Vector{new_panel_type}(undef, length(panels))
    for i in eachindex(panels)
        old_panel = panels[i]
        (; vertices, control_point, normal, strength, radius) = old_panel
        # new_strength = convert_strength(strength, TK, new_kernel)
        new_strength = zero(SVector{new_NK,TF})
        new_panels[i] = new_panel_type(vertices, control_point, normal, new_strength, radius)
    end
    new_strengths = fill(zero(SVector{new_NK, TF}), length(panels))

    return pnl.FMM.UnstructuredGrid{new_kernel, TF, new_NK, NS}(points, meshcells, control_points, normals, new_strengths, potential, velocity, sigma, new_panels, wake_points)
end



#------- try FastMultipole FGS solver -------#

# reset strengths
for i in eachindex(body_fmm.strengths)
    body_fmm.strengths[i] = ones(eltype(body_fmm.strengths))
end
pnl.FMM.grid_2_panels_strength!(body_fmm)
# pnl.FMM.solve_fgs!(body_fmm, Vinf; tolerance=1e-9)
# pnl.FMM.vtk("test_fgs", body_fmm; save_path=".")

#--- try source+dipole panels ---#

panels = change_kernel(body_fmm, pnl.FMM.UniformSourceNormalDoublet)
fgs = pnl.FMM.solve_fgs!(panels, Vinf; tolerance=1e-9, leaf_size=400)
pnl.FMM.vtk("test_fgs_sourcedipole", panels; save_path=".")





# solver = pnl.FMM.BlockGaussSeidel(1000, 1000, 1.0)
# solver = pnl.FMM.LUSolver()
# resid, A, rhs, σ, μ = pnl.FMM.solve_panels!(panels, solver; watertight=true)

# function zero_strength!(panels, i)
#     mask = SVector{2}(j!=i for j in 1:2)
#     for i in eachindex(panels.panels)
#         (; vertices, control_point, normal, strength, radius) = panels.panels[i]
#         panels.panels[i] = eltype(panels.panels)(vertices, control_point, normal, strength .* mask, radius)
#     end
# end

# pnl.FMM.reset_potential_velocity!(panels)

# zero_strength!(panels, 1)
# pnl.FMM.fmm!(panels)



# pnl.FMM.grid_2_panels_strength!(panels)
# zero_strength!(panels, 2)
# pnl.FMM.fmm!(panels)

# for i in eachindex(panels.panels)
#     panels.velocity[i] = panels.velocity[i] + Vinf
# end

# resid = [dot(panels.panels[i].normal, panels.velocity[i]) for i in eachindex(panels.panels)]






# println("\n--- SOLVER: LU Decomposition ---\n")

# strengths_fgs = save_strengths(body_fmm)
# pnl.FMM.reset_potential_velocity!(body_fmm)
# pnl.FMM.apply_freestream!(body_fmm, pnl.FMM.StaticArrays.SVector{3}(Uinfs[:,1]))
# solver = pnl.FMM.LUDecomposition(body_fmm, pnl.FMM.Scheme{pnl.FMM.DirectNeumann, pnl.FMM.FlowTangency})
# pnl.FMM.solve!(body_fmm, solver)
# pnl.FMM.grid_2_panels_strength!(body_fmm)
# pnl.FMM.FastMultipole.direct!(body_fmm)

# pnl.FMM.vtk("test_fgs", body_fmm; save_path=".")

# # test normal velocity
# for (velocity,panel) in zip(body_fmm.velocity, body_fmm.panels)
#     @test isapprox(dot(panel.normal, velocity), 0.0; atol=1e-6)
# end

# strengths_lu = save_strengths(body_fmm)







#--- test FGS solver ---#

# println("\n--- SOLVER: FGS ---\n")

# strengths_lu = save_strengths(body_fmm)

# pnl.FMM.reset_potential_velocity!(body_fmm)
# pnl.FMM.apply_freestream!(body_fmm, pnl.FMM.StaticArrays.SVector{3}(Uinfs[:,1]))
# solver = pnl.FMM.FastGaussSeidel(body_fmm, pnl.FMM.Scheme{pnl.FMM.DirectNeumann,pnl.FMM.FlowTangency}; leaf_size, multipole_threshold, expansion_order)
# pnl.FMM.solve!(body_fmm, solver; max_iterations=10000, tolerance=1e-3, verbose=true, history=true, relaxation=1.4)
# pnl.FMM.grid_2_panels_strength!(body_fmm)
# pnl.FMM.FastMultipole.fmm!(body_fmm; expansion_order, leaf_size, multipole_threshold)

# # test normal velocity
# Us_fmm = [V[i] for i in 1:3, V in body_fmm.velocity]
# println("Maximum error:\t$(maximum(abs.(Us_fmm .- Us_old)))")
# println("Maximum error:\t$(maximum(abs.(Us_fmm .- Us_old) ./ abs.(Us_old)) * 100) %")
# for (velocity,panel) in zip(body_fmm.velocity, body_fmm.panels)
#     @test isapprox(dot(panel.normal, velocity), 0.0; atol=1e-6)
# end
# strengths_fgs = save_strengths(body_fmm)
