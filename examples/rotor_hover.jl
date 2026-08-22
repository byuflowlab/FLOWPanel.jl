## Plate, rotor, with no ground effect
# Author: Timothy Harlow
# Created: December 5, 2025

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
include(joinpath(pnl.examples_path, "dji9443_trailing_edge.jl"))
using FLOWPanel.FastMultipole.StaticArrays
using VSPGeom
using LinearAlgebra: norm

run_name = get(ENV, "RUN_NAME", "rotor_hover")
save_path = joinpath("data", run_name)

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
n_steps = parse(Int, get(ENV, "N_STEPS", string(nt * nrevs)))
n_steps > 0 || error("N_STEPS must be a positive integer, got $(n_steps)")
t_range = range(0.0, step=dt, length=n_steps)# [1:3]

# ==========================================================
# Sensitivity parameters
# ==========================================================
cp_outer=true
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
msh_file, watertight = (
    joinpath(read_path, "dji9443_20260723_30_121_uncapped.msh"),
    false,
)
te_indices_1, te_indices_2 =
    find_dji9443_trailing_edge_indices(msh_file; watertight)

# MSH file
msh = pnl.read_gmsh(msh_file)
nodes, cells = pnl.meshes2nodes_cells(msh)

# scale to proper radius
println("Scaling x-component of geometry to radius R = $R m")
nodes .*= R / maximum(nodes[1, :])

# place-holder shedding
shedding = pnl.noshedding

# --- Construct RigidWakeBody ---
kernel = watertight ? Union{pnl.ConstantSource, pnl.VortexRing} : pnl.VortexRing
DBC = watertight
rotor = pnl.RigidWakeBody{kernel}(nodes, cells, shedding;
            kerneloffset,
            kernelcutoff,
            semiinfinite_wake=false,
            watertight,
            DBC)

# update shedding
shedding1 = pnl.calc_shedding_from_seed(
    rotor.nodes, rotor.cells, te_indices_1[1], te_indices_1[2];
    end_node=te_indices_1[3], normal_jump_tol=0.2,
    max_turn_angle=pi/3, debug=false)
shedding2 = pnl.calc_shedding_from_seed(
    rotor.nodes, rotor.cells, te_indices_2[1], te_indices_2[2];
    end_node=te_indices_2[3], normal_jump_tol=0.2,
    max_turn_angle=pi/3, debug=false)

rotor = pnl.RigidWakeBody{kernel}(rotor.nodes, rotor.cells, [shedding1, shedding2];
                        kerneloffset,
                        kernelcutoff,
                        semiinfinite_wake=false,
                        watertight,
                        ensure_winding=true,
                        DBC)

pnl.write_vtk(joinpath(save_path, run_name), rotor)

println("Rotor: $(rotor.nnodes) nodes, $(rotor.ncells) panels, $(rotor.nsheddings) shedding edges")

## =========================================================
# WAKE SETUP
# ==========================================================

wake_rotor = pnl.PanelWake(rotor;
                nwakerows=n_steps + 1,
                core_size=wake_core_size)

## =========================================================
# SIMULATION SETUP
# ==========================================================
Uinf(t) = Vinf

# Reference frame
frames = pnl.ReferenceFrame(rotor;
    origin = SVector{3}(0.0, 0.0, 0.0),
    v = SVector{3}(0.0, 0.0, 0.0),
    ω_axis = SVector{3}(1.0, 0.0, 0.0),
    ω = -2*pi * RPM/60, # rad/s
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
pressure_bernoulli = pnl.PressureBernoulli(rho;
                        unsteady=true,
                        allow_partial=true,
                        correct_kuttacondition=p_correct_kuttacondition_flag,
                        backend=backend)
force_monitor_bernoulli = pnl.ForceMonitor(length(t_range), 1;
                            i_frame=1,
                            normalization=pnl.RotorNormalization(rho, 2*R, 1),
                            correct_kuttacondition=p_correct_kuttacondition_flag,
                            verbose=false
                        )
pressure_laplace = pnl.PressureLaplace(rotor, rho; verbose=false,
                        unsteady=true)
force_monitor_laplace = pnl.ForceMonitor(length(t_range), 1;
                            i_frame=1,
                            normalization=pnl.RotorNormalization(rho, 2*R, 1),
                            correct_kuttacondition=p_correct_kuttacondition_flag,
                            verbose=false
                        )
# Circulation-based cross-check against the pressure-integrated force history.
kj_monitor = pnl.KuttaJoukowskiForce(rotor, length(t_range), 1;
                rho, backend,
                i_frame=1,
                normalization=pnl.RotorNormalization(rho, 2*R, 1),
                verbose=false
            )
monitors = (
    pressure_laplace,
    force_monitor_laplace,
    pressure_bernoulli,
    force_monitor_bernoulli,
    kj_monitor,
)

println("\nBegin rotor hover simulation ($(n_steps) steps)...")
name = run_name
@time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
    # Das was initialized before constructing the matrixful solver.
    set_Das_eta_kinematic=NaN,
    # set_Das_eta_freestream=0.1,
    monitors,
    body_solvers, backend, verbose=true,
    path=save_path, name,
)

println("Thrust Coefficient (PressureBernoulli + ForceMonitor): ", force_monitor_bernoulli.force[2, :])
println("Thrust Coefficient (PressureLaplace + ForceMonitor): ", force_monitor_laplace.force[2, :])
println("Thrust Coefficient (KuttaJoukowskiForce): ", kj_monitor.force[2, :])

CT_bernoulli = force_monitor_bernoulli.force[2, :]
CT_laplace = force_monitor_laplace.force[2, :]
CT_kj = kj_monitor.force[2, :]

relative_difference(a, b) = abs(a - b) / max(abs(b), eps())

println("\nstep | CT Bernoulli | CT Laplace | CT KJ | rel(B-L) | rel(B-KJ)")
for k in eachindex(CT_bernoulli)
    cb = CT_bernoulli[k]
    cl = CT_laplace[k]
    ck = CT_kj[k]
    println("  $k  |  $(round(cb, sigdigits=6))  |  $(round(cl, sigdigits=6))  |  $(round(ck, sigdigits=6))  |  $(round(relative_difference(cb, cl), sigdigits=4))  |  $(round(relative_difference(cb, ck), sigdigits=4))")
end

bern_lap_rel = norm(CT_bernoulli - CT_laplace) / max(norm(CT_bernoulli), eps())
bern_kj_rel = norm(CT_bernoulli - CT_kj) / max(norm(CT_bernoulli), eps())
lap_kj_rel = norm(CT_laplace - CT_kj) / max(norm(CT_laplace), eps())

println("\nRelative CT history differences:")
println("  Bernoulli vs Laplace: $(bern_lap_rel)")
println("  Bernoulli vs KuttaJoukowski: $(bern_kj_rel)")
println("  Laplace vs KuttaJoukowski: $(lap_kj_rel)")
