## Plate, rotor, with no ground effect
# Author: Timothy Harlow
# Created: December 5, 2025

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
using FLOWPanel.FastMultipole.StaticArrays
using VSPGeom
import GeoIO

## =========================================================
# SIMULATION PARAMETERS
# ==========================================================
magVinf = 0.0001    # Freestream velocity magnitude (m/s)
AOA     = 0.0       # Angle of attack (degrees)
rho     = 1.225     # Air density (kg/m^3)
RPM     = 6000      # Rotation speed (rpm)
Vinf    = magVinf * [0.0, -cosd(AOA), sind(AOA)]
R       = 0.119      # Rotor radius (m)
nrevs   = 1        # Number of revolutions
nt      = 36        # Number of time steps per revolution
dt      = 60 / RPM / nt
n_steps = nt * nrevs
t_range = range(0.0, step=dt, length=n_steps)[1:3]

# ==========================================================
# Sensitivity parameters
# ==========================================================
CPoffset     = R * 1e-6
kerneloffset = R * 1e-3
kernelcutoff = R * 1e-13
p_per_step   = 2
overlap      = 2.0
merge_r_factor = 0.02
merge_r_hash_factor = 0.04
init_Das_eta_kinematic = 0.2
p_correct_kuttacondition_flag = false
wake_core_size = parse(Float64, get(ENV, "WAKE_CORE_SIZE", "1e-3"))

## =========================================================
# ROTOR GEOMETRY
# ==========================================================
read_path   = joinpath(pnl.examples_path, "data")
# stl_file   = joinpath(read_path, "phantom_3_mod3_rev5.stl")

# phantom_3_rebuild_r2.msh
msh_file  = joinpath(read_path, "phantom_3_rebuild_r2.msh")
te_indices_1 = [9, 175, 127]
te_indices_2 = [13, 286, 238]

# # phantom_3_rebuild_r3.msh
# msh_file  = joinpath(read_path, "phantom_3_rebuild_r3.msh")
# te_indices_1 = [8, 523, 223] .+ 1
# te_indices_2 = [12, 997, 697] .+ 1

# # phantom_3_rebuild_r4.msh
# msh_file  = joinpath(read_path, "phantom_3_rebuild_r4.msh")
# te_indices_1 = [7, 952, 4] .+ 1
# te_indices_2 = [3, 478, 0] .+ 1

# STL file
# mesh = VSPGeom.readSTL(stl_file)[1]
# scale = 1/1000 # convert to meters
# radius = 119.38 * scale
# for point in mesh.points
#     point .*= scale
# end

# MSH file
msh = GeoIO.load(msh_file).geometry
nodes, cells = pnl.meshes2nodes_cells(msh)

# scale to proper radius
nodes .*= R / maximum(nodes[1, :])

# place-holder shedding
shedding = pnl.noshedding

# --- Construct RigidWakeBody ---
kernel = Union{pnl.ConstantSource, pnl.VortexRing}
# kernel = pnl.VortexRing
DBC = kernel == pnl.VortexRing ? false : true
rotor = pnl.RigidWakeBody{kernel}(nodes, cells, shedding;
            CPoffset,
            kerneloffset,
            kernelcutoff,
            semiinfinite_wake=false,
            watertight=true,
            DBC)

pnl.write_vtk("rotor_hover_check", rotor)

# update shedding
bbox = (pnl.SVector{3}(-R*1.2, -1.0, -1.0), pnl.SVector{3}(-R*0.1, 1.0, 1.0))
bbox = nothing
shedding1 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_1[1], te_indices_1[2]; bbox, end_node=te_indices_1[3], normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
bbox = (pnl.SVector{3}(R*0.1, -1.0, -1.0), pnl.SVector{3}(R*1.2, 1.0, 1.0))
bbox = nothing
shedding2 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_2[1], te_indices_2[2]; bbox, end_node=te_indices_2[3], normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

rotor = pnl.RigidWakeBody{kernel}(rotor.nodes, rotor.cells, [shedding1, shedding2];
                        CPoffset,
                        kerneloffset,
                        kernelcutoff,
                        semiinfinite_wake=false,
                        watertight=true,
                        ensure_winding=true,
                        DBC)

pnl.write_vtk("rotor_hover", rotor)

println("Rotor: $(rotor.nnodes) nodes, $(rotor.ncells) panels, $(rotor.nsheddings) shedding edges")

## =========================================================
# WAKE SETUP
# ==========================================================

wake_rotor = pnl.PanelParticleWake(rotor;
                nwakerows=1,
                core_size=wake_core_size,
                max_particles=100000,
                method_trailing=pnl.OverlapPPS(overlap, p_per_step),
                method_unsteady=pnl.OverlapPPS(overlap, p_per_step),
                particle_maintenance=pnl.ParticleMaintenance((
                    pnl.MergeParticles(every=1,
                        r=R*merge_r_factor,
                        r_hash=R*merge_r_hash_factor,
                        sigma_relative=false,
                        max_sigma_ratio=2.0,
                        skip_static=true,
                        check_neighboring_cells=false),
                )))

## =========================================================
# SIMULATION SETUP
# ==========================================================
Uinf(t) = Vinf

# Reference frame
frames = pnl.ReferenceFrame(rotor;
    origin = SVector{3}(0.0, 0.0, 0.0),
    v = SVector{3}(0.0, 0.0, 0.0),
    ω_axis = SVector{3}(0.0, 1.0, 0.0),
    ω = 2*pi * RPM/60, # rad/s
    R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    name = "vehicle",
    child_index = Int[],
    dependent_index = [1]
)

pnl.initialize_Das!((rotor,), frames, Uinf, t_range[1], t_range[2] - t_range[1];
    set_Das_eta_kinematic=init_Das_eta_kinematic)

# solver_rotor = pnl.FGSSolver(rotor;
#             max_iterations=500,
#             tolerance=1.0e-6,
#             rlx=1.0,
#             expansion_order=10,
#             multipole_acceptance=0.4,
#             leaf_size=150,
#             shrink=true,
#             recenter=false,
#             inner_iterations=20,
#             reverse_pass=false,
#             verbose=false
#         )
solver_rotor = pnl.Backslash(rotor)

backend = pnl.FastMultipoleBackend()

# Maneuver
maneuver!(frames, systems, wakes, t) = nothing

## =========================================================
# RUN SIMULATION
# ==========================================================
systems      = (rotor,)
wakes        = (wake_rotor,)
body_solvers = (solver_rotor,)
monitors = (pnl.PressureBernoulli(rho; unsteady=true,
                    correct_kuttacondition=p_correct_kuttacondition_flag),
            pnl.ForceMonitor(length(t_range), 1; # un-normalized, global frame
                    i_frame=-1,
                    normalization=pnl.RotorNormalization(rho, 2*R, 1),
                    correct_kuttacondition=p_correct_kuttacondition_flag,
                    # normalization=pnl.NoNormalization(),
                    verbose=false
                ),
            # pnl.KuttaJoukowskiForce(rotor, length(t_range), 1;
            #         rho, backend,
            #         normalization=pnl.RotorNormalization(rho, 2*R, 1)
            #     )
            )

println("\nBegin rotor hover simulation ($(n_steps) steps)...")
name = "rotor_hover"
@time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
    # Das was initialized before constructing the matrixful solver.
    set_Das_eta_kinematic=NaN,
    # set_Das_eta_freestream=0.1,
    monitors,
    body_solvers, backend, verbose=true,
    path="rotor_hover", name,
)

println("Thrust Coefficient: ", monitors[2].force[2,:])
