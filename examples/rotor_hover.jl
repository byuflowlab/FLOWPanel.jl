## Plate, rotor, with no ground effect
# Author: Timothy Harlow
# Created: December 5, 2025

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
using VSPGeom

## =========================================================
# SIMULATION PARAMETERS
# ==========================================================
magVinf = 0.0001    # Freestream velocity magnitude (m/s)
AOA     = 0.0       # Angle of attack (degrees)
rho     = 1.225     # Air density (kg/m^3)
RPM     = 5400      # Rotation speed (rpm)

Vinf    = magVinf * [0.0, -cosd(AOA), sind(AOA)]
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
stl_file   = joinpath(read_path, "phantom_3_mod3_rev5.stl")

R       = 0.12      # Rotor radius
RPM     = 6000      # Rotation speed (rpm)

mesh = VSPGeom.readSTL(stl_file)[1]
scale = 1/1000 # convert to meters
radius = 119.38 * scale
for point in mesh.points
    point .*= scale
end
shedding = pnl.noshedding

# --- Construct RigidWakeBody ---
kernel = Union{pnl.ConstantSource, pnl.VortexRing}
rotor = pnl.RigidWakeBody{kernel}(mesh, shedding;
            CPoffset=1e-14,
            kerneloffset=1e-2,
            kernelcutoff=1e-14,
            semiinfinite_wake=false,
            watertight=true)

pnl.write_vtk("rotor_hover", rotor)

# update shedding
bbox = (pnl.SVector{3}(-radius*1.2, -1.0, -1.0), pnl.SVector{3}(-radius*0.1, 1.0, 1.0))
shedding1 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, 1991, 1989; bbox, end_node=nothing, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
bbox = (pnl.SVector{3}(radius*0.1, -1.0, -1.0), pnl.SVector{3}(radius*1.2, 1.0, 1.0))
shedding2 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, 778, 776; bbox, end_node=nothing, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

rotor = pnl.RigidWakeBody{kernel}(mesh, [shedding1, shedding2],
            CPoffset=1e-14,
            kerneloffset=1e-2,
            kernelcutoff=1e-14,
            semiinfinite_wake=false,
            watertight=true)

println("Rotor: $(rotor.nnodes) nodes, $(rotor.ncells) panels, $(rotor.nsheddings) shedding edges")

## =========================================================
# WAKE SETUP
# ==========================================================

p_per_step = 2
overlap    = 2.0

wake_rotor = pnl.PanelParticleWake(rotor;
                nwakerows=1,
                max_particles=20000,
                method_trailing=pnl.OverlapPPS(overlap, p_per_step),
                method_unsteady=pnl.OverlapPPS(overlap, p_per_step),
                merge_every=1, # merge every step
                merge_r=radius * 0.02, # r_merge for merging particles
                merge_sigma_relative=false, # use relative sigma for merging
                merge_max_sigma_ratio=2.0, # prevents particles of very different strengths from merging
                merge_skip_static=true, # skip merging static particles
                merge_verbose=true)

## =========================================================
# SIMULATION SETUP
# ==========================================================
Uinf(t) = Vinf

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
solver_rotor = pnl.BackslashDirichlet(rotor)

backend = pnl.FastMultipoleBackend()

# Reference frame
frames = pnl.ReferenceFrame(rotor;
    origin = SVector{3}(0.0, 0.0, 0.0),
    v = SVector{3}(0.0, 0.0, 0.0),
    ω_axis = SVector{3}(0.0, 1.0, 0.0),
    ω = 2*pi * RPM/60,
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
systems      = (rotor,)
wakes        = (wake_rotor,)
body_solvers = (solver_rotor,)

println("\nBegin rotor hover simulation ($(n_steps) steps)...")
name = Threads.nthreads() > 1 ? "rotor_hover_mt" : "rotor_hover"
@time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
    set_Das_eta_kinematic=0.1,
    # set_Das_eta_freestream=0.1,
    body_solvers, backend, rho, verbose=true,
    path="rotor_hover", name
)
