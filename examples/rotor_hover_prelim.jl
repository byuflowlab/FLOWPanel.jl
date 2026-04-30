## Plate, rotor, with no ground effect
# Author: Timothy Harlow
# Created: December 5, 2025

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
import LinearAlgebra: norm, I
import FLOWPanel.GeometricTools.Meshes
import GeoIO
using VSPGeom
using Plots

## =========================================================
# SIMULATION PARAMETERS
# ==========================================================
magVinf = 0.0001    # Freestream velocity magnitude (m/s)
AOA     = 0.0       # Angle of attack (degrees)
rho     = 1.225     # Air density (kg/m^3)

Vinf    = magVinf * [cosd(AOA), 0.0, sind(AOA)]
eta     = 0.3

nrevs   = 10        # Number of revolutions
nt      = 36        # Number of time steps per revolution
dt      = 60 / RPM / nt
n_steps = nt * nrevs
t_range = range(0.0, step=dt, length=n_steps)

## =========================================================
# ROTOR GEOMETRY
# ==========================================================
read_path   = joinpath(pnl.examples_path, "data")
mesh_file   = joinpath(read_path, "phantom_3_sharp.msh")
# stl_file    = joinpath(read_path, "phantom_3_sharp.STL")

R       = 0.12      # Rotor radius
RPM     = 6000      # Rotation speed (rpm)

msh = GeoIO.load(mesh_file).geometry
msh = msh |> Meshes.Scale(0.001)
grid = pnl.gt.GridTriangleSurface(msh)

te_file = "/Users/Timothy/Library/CloudStorage/OneDrive-BrighamYoungUniversity/Capstone/flow-panel/meshes/phantom_3_sharp_te1.msh"
te_msh = GeoIO.load(te_file).geometry
te_msh = te_msh |> Meshes.Scale(1.0)

nte       = length(te_msh.vertices) # Number of trailing edge points
te_matrix = zeros(3, nte)

for i in 1:nte
    te_matrix[:, i] .= te_msh.vertices[i].coords[1:3]
end

# n       = length(msh.vertices) # Number of trailing edge points
# matrix = zeros(3, n)

# for i in 1:n
#     matrix[:, i] .= msh.vertices[i].coords[1:3]
# end

#  Test finding indices
# tolerance = R
# points = eachcol(te_matrix)

#     # Calculate the index of the closed point in line to each node
#     indices = Tuple{Int, Int}[] # (index of node, index of line point)

#     for (nodei, node) in enumerate(eachcol(grid._nodes))

#         # Find closest point
#         distance, pointi = findmin(X -> norm(node - X), points)
#         @show(distance, pointi)

#         # Check if it is close enough
#         if distance <= 3
#             push!(indices, (nodei, pointi))
#         end
#     end

#     # Sort the indices by `pointi`
#     sort!(indices, by = x -> x[2] )


# scatter(te_matrix[1,:], te_matrix[2,:],
#     xlabel = "x",
#     ylabel = "y",
#     xlims  = [-120,-20],
#     ylims  = [-2.6, 5.5],
#     markersize = 3
# )
# scatter!(matrix[1,:], matrix[2,:], markersize=3)

# shedding = pnl.calc_shedding(grid._nodes, pnl.grid2cells(grid), te_matrix; tolerance=0.3 * R)
shedding = pnl.noshedding

# --- Construct RigidWakeBody ---
kernel = Union{pnl.ConstantSource, pnl.VortexRing}
rotor = pnl.RigidWakeBody{kernel}(grid, shedding;
            CPoffset=1e-14,
            kerneloffset=1e-2,
            kernelcutoff=1e-14,
            semiinfinite_wake=false,
            watertight=false)

pnl.write_vtk("rotor_hover_prelim", rotor)

# update shedding
bbox = (pnl.SVector{3}(-R*1.2, -1.0, -1.0), pnl.SVector{3}(-R*0.1, 1.0, 1.0))
shedding1 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, 102, 758; bbox, end_node=nothing, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
bbox = (pnl.SVector{3}(R*0.1, -1.0, -1.0), pnl.SVector{3}(R*1.2, 1.0, 1.0))
shedding2 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, 47, 383; bbox, end_node=nothing, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

rotor = pnl.RigidWakeBody{kernel}(grid, [shedding1, shedding2],
            CPoffset=1e-14,
            kerneloffset=1e-2,
            kernelcutoff=1e-14,
            semiinfinite_wake=false,
            watertight=true)
pnl.write_vtk("rotor_hover_prelim", rotor)

# --- Construct RigidWakeBody ---
# kernel = Union{pnl.ConstantSource, pnl.VortexRing}
# rotor = pnl.RigidWakeBody{kernel}(grid, shedding;
#             CPoffset=1e-14,
#             kerneloffset=1e-2,
#             kernelcutoff=1e-14,
#             semiinfinite_wake=false,
#             watertight=true)

println("Rotor: $(rotor.nnodes) nodes, $(rotor.ncells) panels, $(rotor.nsheddings) shedding edges")

## =========================================================
# WAKE SETUP
# ==========================================================
das_offset = 0.05                               # Das scale (fraction of unit Vinf direction)

rotor.Das[1] .= repeat(Vinf / magVinf * das_offset, 1, size(rotor.Das[1], 2))

p_per_step = 2
overlap    = 1.3

wake_rotor = pnl.PanelParticleWake(rotor;
                nwakerows=3,
                max_particles=20000,
                method_trailing=pnl.OverlapPPS(overlap, p_per_step),
                method_unsteady=pnl.OverlapPPS(overlap, p_per_step))

## =========================================================
# SIMULATION SETUP
# ==========================================================
Uinf(t) = Vinf

expansion_order      = 10
multipole_acceptance = 0.4
leaf_size            = 20
backend = pnl.FastMultipoleBackend(; expansion_order, multipole_acceptance, leaf_size)

# solver_rotor = pnl.FGSSolver(rotor;
#             max_iterations=500,
#             tolerance=1.0e-6,
#             rlx=1.0,
#             expansion_order=14,
#             multipole_acceptance,
#             leaf_size=150,
#             shrink=true,
#             recenter=false,
#             inner_iterations=20,
#             reverse_pass=false,
#             verbose=false
#         )

solver_rotor = pnl.BackslashDirichlet(rotor)

backend = pnl.FastMultipoleBackend()

# Reference frame
frames = pnl.ReferenceFrame(rotor;
    origin = SVector{3}(0.0, 0.0, 0.0),
    v = SVector{3}(0.0, 0.0, 0.0),
    ω_axis = SVector{3}(0.0, 0.0, 1.0),
    ω = 2*pi * RPM*60,
    R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    name = "vehicle",
    child_index = Int[],
    dependent_index = [1]
)

# Maneuver
maneuver!(frames, systems, wakes, t) = nothing

## =========================================================
# RUN SIMULATION
# ==========================================================
println("\nBegin rotor hover simulation ($(n_steps) steps)...")
@time pnl.simulate!(rotor, wake_rotor, frames, maneuver!, Uinf, t_range;
    body_solvers, backend, rho, verbose=true,
    path="rotor_hover", name="rotor_hover"
)