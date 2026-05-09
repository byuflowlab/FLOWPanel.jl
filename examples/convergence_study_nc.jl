## Plate, rotor, with no ground effect
# C_T vs. number of cells
# Author: Timothy Harlow
# Created: December 5, 2025

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
using VSPGeom
import GeoIO
using DelimitedFiles

function build_rotor(msh_file, te_1, te_2, R=0.119; find_te=false, kerneloffset=R*0.001)
    # ==========================================================
    # Sensitivity parameters
    # ==========================================================
    CPoffset     = R * 1e-6
    kernelcutoff = R * 1e-13
    p_per_step   = 2
    overlap      = 2.0
    merge_r_factor = 0.02
    merge_r_hash_factor = 0.04
    init_Das_eta_kinematic = 0.2
    p_correct_kuttacondition_flag = false
    wake_core_size = parse(Float64, get(ENV, "WAKE_CORE_SIZE", "1e-3"))

    # MSH files
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

    if find_te
        pnl.write_vtk("rotor_hover_check", rotor)
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

    return rotor
end


function run_simulation(rotor, run_name; R=0.119, Vinf=zeros(3), RPM=5400, nt=36, nrevs=10)
    dt      = 60 / RPM / nt
    n_steps = nt * nrevs
    t_range = range(0.0, step=dt, length=n_steps)
    
    # =========================================================
    # WAKE SETUP
    # ==========================================================
    p_per_step = 2
    overlap    = 2.0
    wake_rotor = pnl.PanelParticleWake(rotor;
                    nwakerows=1,
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
    solver_rotor = pnl.Backslash(rotor)

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
    monitors = (pnl.ForceMonitor(length(t_range), 1; # un-normalized, global frame
                        i_frame=-1, 
                        normalization=pnl.RotorNormalization(rho, 2*R, 1)
                    ),
                )

    nc = rotor.ncells # number of cells on this run

    println("\nBegin rotor hover simulation ($(n_steps) steps, $(nc) cells)...")
    name = run_name
    @time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
        set_Das_eta_kinematic=0.1,
        # set_Das_eta_freestream=0.1,
        monitors,
        body_solvers, backend, rho, verbose=true,
        path=run_name, name
    )

    out_name = "CT_v_t_hover_RPM"*"$RPM"*"_nc$nc"*"_pps$p_per_step"*"_nt$nt"*"_overlap$overlap"*"_kerneloff0.00012"*".csv"
    writedlm(out_name, monitors[1].force, ',')

end

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

# =========================================================
# ROTOR GEOMETRY
# ==========================================================
read_path    = joinpath(pnl.examples_path, "data")
R       = 0.119      # Rotor radius

# ~uniform cell size
msh_file_1k  = joinpath(read_path, "phantom_3_1k.msh")
msh_file_4k  = joinpath(read_path, "phantom_3_4k.msh")
msh_file_8k  = joinpath(read_path, "phantom_3_8k.msh")
msh_file_16k  = joinpath(read_path, "phantom_3_16k.msh")

te_1k_1 = [10, 73]
te_1k_2 = [13, 144]

te_4k_1 = [9, 175]
te_4k_2 = [13, 286]

te_8k_1 = [10, 196]
te_8k_2 = [13, 406]

te_16k_1 = [10, 268]
te_16k_2 = [13, 562]

# rotor_1k  = build_rotor(msh_file_1k, te_1k_1, te_1k_2)
# rotor_4k  = build_rotor(msh_file_4k, te_4k_1, te_4k_2)
# rotor_8k  = build_rotor(msh_file_8k, te_8k_1, te_8k_2)
# rotor_16k = build_rotor(msh_file_16k, te_16k_1, te_16k_2)

# adaptive cell size (smaller toward LE and TE)
msh_adapt_12k  = joinpath(read_path, "phantom_3_adaptive_12k.msh")
msh_adapt_24k  = joinpath(read_path, "phantom_3_adaptive_24k.msh")

te_a_12k_1 = [10, 398]
te_a_12k_2 = [13, 994]

te_a_24k_1 = [10, 553]
te_a_24k_2 = [13, 1396]

rotor_adapt_12k = build_rotor(msh_adapt_12k, te_a_12k_1, te_a_12k_2)
rotor_adapt_24k = build_rotor(msh_adapt_24k, te_a_24k_1, te_a_24k_2)

# =========================================================
# RUN SIMULATION
# ==========================================================
run_simulation(rotor_adapt_12k, "rotor_hover_12k"; Vinf=Vinf)
run_simulation(rotor_adapt_24k, "rotor_hover_24k"; Vinf=Vinf)