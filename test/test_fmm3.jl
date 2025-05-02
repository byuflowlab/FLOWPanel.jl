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

# Check if PyPlot is available
# try
#     import FLOWPanel: plt, @L_str
#     global pyplot = true
# catch
#     global pyplot = false
# end
import PyPlot as plt
using LaTeXStrings
const pyplot = false



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
if paraview
    str = save_path*"/"

    # Save body as a VTK
    str *= pnl.save(body, run_name; path=save_path, debug=true)

    # Call Paraview
    run(`paraview --data=$(str)`)
end

#--- use FastMultipole to compute induced velocity ---#

sigma = 1e-2
println("\n--- INDUCED VELOCITY: FMM ---\n")

# Convert body to FMM body
body_fmm = pnl.body2fmmbody(body; sigma)

# Set freestream at each FMM panel
body_fmm.potential .= 0.0
for (i, Uinf) in zip(eachindex(body_fmm.velocity), eachcol(Uinfs))
    body_fmm.velocity[i] -= body_fmm.velocity[i]
    body_fmm.velocity[i] += Uinf
end

# Evaluate mutually-induced velocity of panels for benchmark
expansion_order, leaf_size_source, multipole_threshold = 30, 30, 0.5
ε_tol=pnl.FMM.FastMultipole.PowerAbsoluteVelocity(1e-8, false)
# ε_tol=nothing
@time _, _, target_tree, source_tree, m2l_list, direct_list, _ = pnl.FMM.FastMultipole.fmm!(body_fmm; expansion_order, leaf_size_source, multipole_threshold, ε_tol)
@show length(m2l_list) length(direct_list)
# pnl.FMM.FastMultipole.direct!(body_fmm)

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

# println("Maximum error:\t$(maximum(abs.(Us_fmm .- Us_old)))")
# println("Maximum error:\t$(maximum(abs.(Us_fmm .- Us_old) ./ abs.(Us_old)) * 100) %")

#--- direct calcs just in case ---#

println("\n--- INDUCED VELOCITY: DIRECT ---\n")

# Convert body to FMM body
body_fmm = pnl.body2fmmbody(body; sigma)

# Set freestream at each FMM panel
for (i, Uinf) in zip(eachindex(body_fmm.velocity), eachcol(Uinfs))
    body_fmm.velocity[i] -= body_fmm.velocity[i]
    body_fmm.velocity[i] += Uinf
end

# Evaluate mutually-induced velocity of panels for benchmark
@time pnl.FMM.FastMultipole.direct!(body_fmm)

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

Us_direct = [V[i] for i in 1:3, V in body_fmm.velocity]
# display(Us_old)

println("Maximum error:\t$(maximum(abs.(Us_fmm .- Us_direct)))")
println("Maximum error:\t$(maximum(abs.(Us_fmm .- Us_direct) ./ abs.(Us_direct)) * 100) %")


function save_strengths(body)
    strengths = zeros(length(body.strengths))
    for i in eachindex(strengths)
        strengths[i] = body.strengths[i][1]
    end
    return strengths
end

# #--- test direct solver ---#
# println("\n--- SOLVER: LU DECOMPOSITION ---\n")

# strengths_old = save_strengths(body_fmm)
# pnl.FMM.reset_potential_velocity!(body_fmm)
# pnl.FMM.apply_freestream!(body_fmm, pnl.FMM.StaticArrays.SVector{3}(Uinfs[:,1]))
# solver = pnl.FMM.LUDecomposition(body_fmm, pnl.FMM.Scheme{pnl.FMM.DirectNeumann, pnl.FMM.FlowTangency})
# pnl.FMM.solve!(body_fmm, solver)
# pnl.FMM.grid_2_panels_strength!(body_fmm)
# pnl.FMM.FastMultipole.direct!(body_fmm)

# strengths_lu = save_strengths(body_fmm)

# # test normal velocity
# for (velocity,panel) in zip(body_fmm.velocity, body_fmm.panels)
#     @test isapprox(dot(panel.normal, velocity), 0.0; atol=1e-6)
# end


#=

#--- test FGS solver ---#

println("\n--- SOLVER: FGS ---\n")

strengths_lu = save_strengths(body_fmm)

pnl.FMM.reset_potential_velocity!(body_fmm)
pnl.FMM.apply_freestream!(body_fmm, pnl.FMM.StaticArrays.SVector{3}(Uinfs[:,1]))
solver = pnl.FMM.FastGaussSeidel(body_fmm, pnl.FMM.Scheme{pnl.FMM.DirectNeumann,pnl.FMM.FlowTangency}; leaf_size, multipole_threshold, expansion_order)
pnl.FMM.solve!(body_fmm, solver; max_iterations=10000, tolerance=1e-3, verbose=true, history=true, relaxation=1.4)
pnl.FMM.grid_2_panels_strength!(body_fmm)
pnl.FMM.FastMultipole.fmm!(body_fmm; expansion_order, leaf_size, multipole_threshold)

# test normal velocity
Us_fmm = [V[i] for i in 1:3, V in body_fmm.velocity]
println("Maximum error:\t$(maximum(abs.(Us_fmm .- Us_old)))")
println("Maximum error:\t$(maximum(abs.(Us_fmm .- Us_old) ./ abs.(Us_old)) * 100) %")
for (velocity,panel) in zip(body_fmm.velocity, body_fmm.panels)
    @test isapprox(dot(panel.normal, velocity), 0.0; atol=1e-6)
end
strengths_fgs = save_strengths(body_fmm)
=#
