## Plate, rotor, with no ground effect
# C_T vs. number of cells
# Author: Timothy Harlow
# Created: December 5, 2025

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
using VSPGeom
import GeoIO
using DelimitedFiles

## =========================================================
# SIMULATION PARAMETERS
# ==========================================================
magVinf = 0.0001    # Freestream velocity magnitude (m/s)
AOA     = 0.0       # Angle of attack (degrees)
rho     = 1.225     # Air density (kg/m^3)
RPM     = 5400      # Rotation speed (rpm)
Vinf    = magVinf * [0.0, -cosd(AOA), sind(AOA)]
nrevs   = 10        # Number of revolutions
nt      = 36        # Number of time steps per revolution
dt      = 60 / RPM / nt
n_steps = nt * nrevs
t_range = range(0.0, step=dt, length=n_steps)

## =========================================================
# ROTOR GEOMETRY
# ==========================================================
read_path    = joinpath(pnl.examples_path, "data")
msh_file_1k  = joinpath(read_path, "phantom_3_1k.msh")
msh_file_4k  = joinpath(read_path, "phantom_3_4k.msh")
msh_file_8k  = joinpath(read_path, "phantom_3_8k.msh")

R       = 0.119      # Rotor radius

te_1k_1 = [10, 73]
te_1k_2 = [13, 144]

te_4k_1 = [9, 175]
te_4k_2 = [13, 286]

te_8k_1 = [10, 196]
te_8k_2 = [13, 406]

function build_rotor(msh_file, te_1, te_2, R=0.119; find_te=false)

    # MSH files
    msh  = GeoIO.load(msh_file).geometry
    mesh = pnl.gt.GridTriangleSurface(msh)

    # scale to proper radius
    mesh._nodes .*= R / maximum(mesh._nodes[1, :])

    # place-holder shedding
    shedding = pnl.noshedding

    # --- Construct RigidWakeBody ---
    kernel = Union{pnl.ConstantSource, pnl.VortexRing}
    rotor = pnl.RigidWakeBody{kernel}(mesh, shedding;
                CPoffset=1e-14,
                kerneloffset=1e-2,
                kernelcutoff=1e-14,
                semiinfinite_wake=false,
                watertight=true)

    if find_te
        pnl.write_vtk("rotor_hover", rotor)
        println("Find te indices in Paraview")
        return nothing
    end

    # update shedding
    # bbox = nothing
    bbox1 = (pnl.SVector{3}(R*0.1, -1.0, -1.0), pnl.SVector{3}(R*1.2, 1.0, 1.0))
    shedding1 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_1[1], te_1[2]; bbox=bbox1, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
    # bbox = nothing
    bbox2 = (pnl.SVector{3}(-R*1.2, -1.0, -1.0), pnl.SVector{3}(-R*0.1, 1.0, 1.0))
    shedding2 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_2[1], te_2[2]; bbox=bbox2, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

    rotor = pnl.RigidWakeBody{kernel}(rotor.nodes, rotor.cells, [shedding1, shedding2],
                            CPoffset=1e-6,
                            kerneloffset=R*0.001,
                            kernelcutoff=1e-14,
                            semiinfinite_wake=false,
                            watertight=true,
                            ensure_winding=true)

    pnl.write_vtk("rotor_hover", rotor)

    println("Rotor: $(rotor.nnodes) nodes, $(rotor.ncells) panels, $(rotor.nsheddings) shedding edges")

    return rotor
end

rotor_1k = build_rotor(msh_file_1k, te_1k_1, te_1k_2)
rotor_4k = build_rotor(msh_file_4k, te_4k_1, te_4k_2)
rotor_8k = build_rotor(msh_file_8k, te_8k_1, te_8k_2)

## =========================================================
# WAKE SETUP
# ==========================================================
p_per_step = 2
overlap    = 2.0

wake_rotor = pnl.PanelParticleWake(rotor_8k;
                nwakerows=2,
                max_particles=100000,
                method_trailing=pnl.OverlapPPS(overlap, p_per_step),
                method_unsteady=pnl.OverlapPPS(overlap, p_per_step),
                particle_maintenance=pnl.ParticleMaintenance((
                    pnl.MergeParticles(every=1,
                        r=R*0.02,
                        r_hash=R*0.04,
                        sigma_relative=false,
                        max_sigma_ratio=2.0,
                        skip_static=true,
                        check_neighboring_cells=false),
                )))

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
solver_rotor = pnl.Backslash(rotor_8k)

backend = pnl.FastMultipoleBackend()

# Reference frame
frames = pnl.ReferenceFrame(rotor_8k;
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
systems      = (rotor_8k,)
wakes        = (wake_rotor,)
body_solvers = (solver_rotor,)
monitors = (pnl.ForceMonitor(length(t_range), 1; # un-normalized, global frame
                    i_frame=-1, 
                    normalization=pnl.RotorNormalization(rho, 2*R, 1)
                ),
            )

nc = rotor_8k.ncells # number of cells on this run

println("\nBegin rotor hover simulation ($(n_steps) steps, $(nc) cells)...")
name = "rotor_hover_8k"
@time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
    set_Das_eta_kinematic=0.1,
    # set_Das_eta_freestream=0.1,
    monitors,
    body_solvers, backend, rho, verbose=true,
    path="rotor_hover_8k", name
)


out_name = "CT_v_t_hover_RPM"*"$RPM"*"_nc$nc"*"_pps$p_per_step"*"_nt$nt"*"_overlap$overlap"*"_kerneloff$kerneloffset"*".csv"
writedlm(out_name, monitors[1].force, ',')