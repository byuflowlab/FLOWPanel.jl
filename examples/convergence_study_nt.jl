## Plate, rotor, with no ground effect
# Author: Timothy Harlow
# Created: December 5, 2025

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
using VSPGeom
import GeoIO
using DelimitedFiles

function run_convergence(p_per_step, overlap, nt, kerneloffset)
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
    dt      = 60 / RPM / nt
    n_steps = nt * nrevs
    t_range = range(0.0, step=dt, length=n_steps)

    ## =========================================================
    # ROTOR GEOMETRY
    # ==========================================================
    read_path   = joinpath(pnl.examples_path, "data")
    stl_file   = joinpath(read_path, "phantom_3_mod3_rev5.stl")
    msh_file  = joinpath(read_path, "phantom_3_rebuild_r2.msh")

    R       = 0.12      # Rotor radius
    # RPM     = 6000      # Rotation speed (rpm)

    # STL file
    # mesh = VSPGeom.readSTL(stl_file)[1]
    # scale = 1/1000 # convert to meters
    # radius = 119.38 * scale
    # for point in mesh.points
    #     point .*= scale
    # end

    # MSH file
    msh = GeoIO.load(msh_file).geometry
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

    # pnl.write_vtk("rotor_hover", rotor)

    te_indices_1 = [9, 175, 127]
    te_indices_2 = [13, 286, 238]

    # update shedding
    bbox = (pnl.SVector{3}(-R*1.2, -1.0, -1.0), pnl.SVector{3}(-R*0.1, 1.0, 1.0))
    bbox = nothing
    shedding1 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_1[1], te_indices_1[2]; bbox, end_node=te_indices_1[3], normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
    bbox = (pnl.SVector{3}(R*0.1, -1.0, -1.0), pnl.SVector{3}(R*1.2, 1.0, 1.0))
    bbox = nothing
    shedding2 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_2[1], te_indices_2[2]; bbox, end_node=te_indices_2[3], normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

    rotor = pnl.RigidWakeBody{kernel}(rotor.nodes, rotor.cells, [shedding1, shedding2],
                            CPoffset=1e-6,
                            kerneloffset=kerneloffset,
                            kernelcutoff=1e-14,
                            semiinfinite_wake=false,
                            watertight=true,
                            ensure_winding=true)

    pnl.write_vtk("rotor_hover", rotor)

    println("Rotor: $(rotor.nnodes) nodes, $(rotor.ncells) panels, $(rotor.nsheddings) shedding edges")

    ## =========================================================
    # WAKE SETUP
    # ==========================================================
    wake_rotor = pnl.PanelParticleWake(rotor;
                    nwakerows=1,
                    max_particles=100000,
                    method_trailing=pnl.OverlapPPS(overlap, p_per_step),
                    method_unsteady=pnl.OverlapPPS(overlap, p_per_step),
                    merge_every=1, # merge every step
                    merge_r=R*0.02, # r_merge for merging particles
                    merge_r_hash=R*0.04, # r_merge for hashing particles
                    merge_sigma_relative=false, # use relative sigma for merging
                    merge_max_sigma_ratio=2.0, # prevents particles of very different strengths from merging
                    merge_skip_static=true, # skip merging static particles
                    check_neighboring_cells=false, # check neighboring cells for merging (prevents merging across large gaps)
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

    println("\nBegin rotor hover simulation ($(n_steps) steps)...")
    name = "rotor_convergence"
    @time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
        set_Das_eta_kinematic=0.1,
        # set_Das_eta_freestream=0.1,
        monitors,
        body_solvers, backend, rho, verbose=true,
        path="rotor_convergence", name
    )

    out_name = "CT_v_t_hover_RPM"*"$RPM"*"_pps$p_per_step"*"_nt$nt"*"_overlap$overlap"*"_kerneloff$kerneloffset"*".csv"
    writedlm(out_name, monitors[1].force, ',')
end

## =========================================================
# RUN CONVERGENCE STUDY
# ==========================================================
R = 0.12
np = 120 # number of particles per revolution

# Constant values
# nt = 36
# p_per_step = 2
kerneloffset = R*0.001
overlap0 = 2.0


# Variable values
nts         = [10, 20, 40, 60, 120]
p_per_steps = Int.(np ./ nts)

for i in eachindex(p_per_steps)
    nt = nts[i]
    p_per_step = p_per_steps[i]
    sigma0 = overlap0 * (2*pi*R) / (2*nt*p_per_step) # maintain sigma constant through each run
    overlap    = sigma0 * (nt*p_per_step) * 2 ./ (2*pi*R)
    println("Simulation settings:\n p_per_step: $p_per_step\n nt: $nt\n sigma: $sigma0\n λ: $overlap")
    @show Threads.nthreads()

    run_convergence(p_per_step, overlap, nt, kerneloffset)
end