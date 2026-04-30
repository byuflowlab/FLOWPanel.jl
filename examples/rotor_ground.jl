## Plate, rotor, with no ground effect
# Author: Timothy Harlow
# Created: December 5, 2025

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
using VSPGeom
using DelimitedFiles
using GeoIO
using Meshes

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
mesh_file   = joinpath(read_path, "phantom_3_sharp.msh")

# mesh = VSPGeom.readSTL(stl_file)[1]
scale = 1/1000 # convert to meters
radius = 119.38 * scale
# for point in mesh.points
#     point .*= scale
# end

msh = GeoIO.load(mesh_file).geometry
msh = msh |> Meshes.Scale(scale)
mesh = pnl.gt.GridTriangleSurface(msh)

shedding = pnl.noshedding

# --- Construct RigidWakeBody ---
kernel = Union{pnl.ConstantSource, pnl.VortexRing}
rotor = pnl.RigidWakeBody{kernel}(mesh, shedding;
            CPoffset=1e-14,
            kerneloffset=1e-2,
            kernelcutoff=1e-14,
            semiinfinite_wake=false,
            watertight=true)

pnl.write_vtk("rotor_ground", rotor)

# update shedding
bbox = nothing#(pnl.SVector{3}(-radius*1.2, -1.0, -1.0), pnl.SVector{3}(-radius*0.1, 1.0, 1.0))
shedding1 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, 64, 585; bbox, end_node=nothing, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
bbox = nothing#(pnl.SVector{3}(radius*0.1, -1.0, -1.0), pnl.SVector{3}(radius*1.2, 1.0, 1.0))
shedding2 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, 119, 1520; bbox, end_node=nothing, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

rotor = pnl.RigidWakeBody{kernel}(mesh, [shedding1, shedding2],
            CPoffset=1e-14,
            kerneloffset=1e-2,
            kernelcutoff=1e-14,
            semiinfinite_wake=false,
            watertight=true)

println("Rotor: $(rotor.nnodes) nodes, $(rotor.ncells) panels, $(rotor.nsheddings) shedding edges")

## =============================================================================
# GROUND GEOMETRY — flat ground plane below the rotor
# =============================================================================
ground_y        = -0.5 * radius                 # (m) ground plane y-position
ground_center   = [0.0, ground_y, 0.0]          # center
ground_normal   = [0.0, 1.0, 0.0]               # pointing up toward rotor
ground_radius   = 0.5                           # (m) radius of ground disc
ground_panel_length = ground_radius/20          # (m) side length of equilateral triangles

ground = pnl.FlatGround(ground_center, ground_normal, ground_radius;
                         panel_length=ground_panel_length)

println("Ground: $(ground.nnodes) nodes, $(ground.ncells) panels")


## =========================================================
# WAKE SETUP
# ==========================================================
p_per_step = 2
overlap    = 2.125

wake_rotor = pnl.PanelParticleWake(rotor;
                nwakerows=1,
                max_particles=100000,
                method_trailing=pnl.OverlapPPS(overlap, p_per_step),
                method_unsteady=pnl.OverlapPPS(overlap, p_per_step),
                merge_every=1, # merge every step
                merge_r=radius * 0.00422, # r_merge for merging particles
                merge_sigma_relative=false, # use relative sigma for merging
                merge_max_sigma_ratio=2.0, # prevents particles of very different strengths from merging
                merge_skip_static=true, # skip merging static particles
                merge_verbose=true)

wake_ground = nothing

## =========================================================
# SIMULATION SETUP
# ==========================================================
Uinf(t) = Vinf

solver_rotor = pnl.FGSSolver(rotor;
            max_iterations=500,
            tolerance=1.0e-6,
            rlx=1.0,
            expansion_order=10,
            multipole_acceptance=0.4,
            leaf_size=150,
            shrink=true,
            recenter=false,
            inner_iterations=20,
            reverse_pass=false,
            verbose=false
        )
solver_rotor = pnl.BackslashDirichlet(rotor)
solver_ground = pnl.FlatGroundSolver(ground)

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

## =============================================================================
# FLOW TANGENCY MONITOR — check U·n on ground at every step
# =============================================================================
function tangency_monitor(systems, wakes, i_step)
    ground = systems[2]
    normals = pnl.calc_normals(ground)
    Udotn = sum(ground.velocity .* normals, dims=1)
    rms = sqrt(sum(Udotn .^ 2) / ground.ncells)
    maxerr = maximum(abs.(Udotn))
    println("    Ground tangency step $i_step: RMS(U·n) = $(round(rms; sigdigits=4)),  " *
            "max|U·n| = $(round(maxerr; sigdigits=4)),  " *
            "RMS/Vinf = $(round(rms/magVinf*100; sigdigits=4))%")
end

## =========================================================
# RUN SIMULATION
# ==========================================================
systems      = (rotor, ground)
wakes        = (wake_rotor, wake_ground)
body_solvers = (solver_rotor, solver_ground)
monitors = (pnl.ForceMonitor(length(t_range), 1; # un-normalized, global frame
                i_frame=-1, rho=1.0, Sref=1.0, Lref=1.0, TF=Float64),
            )

println("\nBegin rotor + ground simulation ($(n_steps) steps)...")
name = Threads.nthreads() > 1 ? "rotor_ground_mt" : "rotor_ground"
@time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
    set_Das_eta_kinematic=0.1,
    # set_Das_eta_freestream=0.1,
    monitors,
    body_solvers, backend, rho, verbose=true,
    path="rotor_ground", name
)

## =============================================================================
# RESULTS AND FINAL FLOW TANGENCY VERIFICATION
# =============================================================================
println("\n=== RESULTS ===")

# Forces
F_rotor = sum(rotor.F, dims=2)
F_ground = sum(ground.F, dims=2)
println("Rotor total force:   ", F_rotor)
println("Ground total force: ", F_ground)

# Final flow tangency check for both bodies
println("\n--- Flow Tangency (from last simulation step) ---")

for (i, (body, label)) in enumerate(zip(systems, ("Rotor", "Ground")))
    normals = pnl.calc_normals(body)
    Udotn = sum(body.velocity .* normals, dims=1)
    rms = sqrt(sum(Udotn .^ 2) / body.ncells)
    maxr = maximum(abs.(Udotn))
    println("  $label: RMS(U·n) = $(round(rms; sigdigits=4)),  " *
            "max|U·n| = $(round(maxr; sigdigits=4)),  " *
            "RMS/Vinf = $(round(rms/magVinf*100; sigdigits=4))%")
end

out_name = "CT_v_t_hover_ground.csv"
writedlm(out_name, monitors[1].CF, ',')