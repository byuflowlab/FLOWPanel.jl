## Plate, rotor, with no ground effect
# C_T vs. number of cells
# Author: Timothy Harlow
# Created: December 5, 2025

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
using VSPGeom
import GeoIO
using LinearAlgebra: norm
using DelimitedFiles

function build_rotor(msh_file, te_indices_1, te_indices_2, R=0.119; find_te=false, kerneloffset=R*0.001)
    # ==========================================================
    # Sensitivity parameters
    # ==========================================================
    run_name  = "rotor_hover_convergence_nc"
    save_path = joinpath("data", run_name)

    kerneloffset_panel = parse(Float64, get(ENV, "KERNELOFFSET_PANEL", string(R * 1e-10)))
    kerneloffset_targets = parse(Float64, get(ENV, "KERNELOFFSET_TARGETS", get(ENV, "KERNELOFFSET", "1e-3")))
    kernelcutoff = R * 1e-13
    shedding_r_over_R = 0.1

    radial_dimension = occursin("dji9443", msh_file) ? 2 : 1 # this might be wrong for non-dji9443

    msh = GeoIO.load(msh_file).geometry
    nodes, cells = pnl.meshes2nodes_cells(msh)
    nodes .*= R / maximum(nodes[radial_dimension, :])

    shedding = pnl.noshedding
    kernel = Union{pnl.ConstantSource, pnl.VortexRing}
    DBC = kernel == pnl.VortexRing ? false : true
    rotor = pnl.RigidWakeBody{kernel}(nodes, cells, shedding;
        kerneloffset=kerneloffset_panel, kerneloffset_panel, kerneloffset_targets, kernelcutoff,
        semiinfinite_wake=false, watertight=true, DBC)

    0.0 <= shedding_r_over_R <= 1.0 || error("shedding_r_over_R must be between 0 and 1")

    function make_shedding_bbox(nodes, seed_nodes, radial_dimension, R, shedding_r_over_R)
        radial_midpoint = sum(nodes[radial_dimension, seed_nodes]) / length(seed_nodes)
        radial_sign = sign(radial_midpoint)
        radial_sign == 0 && error("Seed edge lies on the rotor axis; cannot determine shedding side")

        lower = [minimum(nodes[i, :]) for i in 1:size(nodes, 1)]
        upper = [maximum(nodes[i, :]) for i in 1:size(nodes, 1)]
        padding = max(sqrt(eps(eltype(nodes))) * R, R * 1e-6)
        lower .-= padding
        upper .+= padding

        radial_cutoff = shedding_r_over_R * R
        if radial_sign > 0
            lower[radial_dimension] = radial_cutoff - padding
        else
            upper[radial_dimension] = -radial_cutoff + padding
        end

        return (pnl.SVector{3}(lower...), pnl.SVector{3}(upper...))
    end

    # save vtk file to inspect for shedding nodes
    if find_te
        vtk_path = joinpath(save_path, "rotor_initial")
        pnl.write_vtk(vtk_path, rotor)
        println("Find TE nodes in Paraview")
        return nothing
    end

    bbox1 = make_shedding_bbox(rotor.nodes, te_indices_1[1:2], radial_dimension, R, shedding_r_over_R)
    shedding1 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_1[1], te_indices_1[2];
        bbox=bbox1, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
    bbox2 = make_shedding_bbox(rotor.nodes, te_indices_2[1:2], radial_dimension, R, shedding_r_over_R)
    shedding2 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_2[1], te_indices_2[2];
        bbox=bbox2, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

    rotor = pnl.RigidWakeBody{kernel}(rotor.nodes, rotor.cells, [shedding1, shedding2];
        kerneloffset=kerneloffset_panel, kerneloffset_panel, kerneloffset_targets, kernelcutoff,
        semiinfinite_wake=false, watertight=true,
        ensure_winding=true, DBC)

    function shedding_root_r_over_R(nodes, shedding, cells, radial_dimension, R)
        isempty(shedding) && return NaN
        root_edge = shedding[:, end]
        pi, nia, nib = root_edge[1], root_edge[2], root_edge[3]
        edge_nodes = cells[[nia, nib], pi]
        midpoint = (nodes[:, edge_nodes[1]] + nodes[:, edge_nodes[2]]) / 2
        return midpoint[radial_dimension] / R
    end

    println("Requested shedding root at |r/R| >= $(shedding_r_over_R)")
    println("  shedding1 root midpoint r/R = $(shedding_root_r_over_R(rotor.nodes, shedding1, rotor.cells, radial_dimension, R))")
    println("  shedding2 root midpoint r/R = $(shedding_root_r_over_R(rotor.nodes, shedding2, rotor.cells, radial_dimension, R))")

    println("Rotor: $(rotor.nnodes) nodes, $(rotor.ncells) panels, $(rotor.nsheddings) shedding edges")

    return rotor
end


function run_simulation(rotor, run_name, msh_file; R=0.119, Vinf=zeros(3), RPM=5400, nt=36, nrevs=10)
    run_name = "rotor_hover_convergence_nc"
    save_path = parse(Bool, get(ENV, "SAVE_VTK", "true")) ? joinpath("data", run_name) : nothing

    rho = 1.179 # from NASA paper
    dt      = 60 / RPM / nt
    spinup_revs = max(0.0, parse(Float64, get(ENV, "SPINUP_REVS", "0.0")))
    spinup_start_fraction = clamp(parse(Float64, get(ENV, "SPINUP_START_FRACTION", "0.05")), eps(), 1.0)
    spinup_duration = spinup_revs * 60 / RPM
    spinup_steps = ceil(Int, spinup_duration / dt)
    n_steps = round(Int, nt * nrevs) + spinup_steps
    t_range = range(0.0, step=dt, length=n_steps)

    ramp_nrev = 5.0
    ramp_t = ramp_nrev * 60 / RPM
    magVinf_start = 0.0
    cylinder_depth = 4R

    cp_outer=true
    kerneloffset_panel = parse(Float64, get(ENV, "KERNELOFFSET_PANEL", string(R * 1e-10)))
    kerneloffset_targets = parse(Float64, get(ENV, "KERNELOFFSET_TARGETS", get(ENV, "KERNELOFFSET", "1e-3")))
    kernelcutoff = R * 1e-13
    p_per_step = 2
    overlap = 3.0
    merge_r_factor = 0.02
    merge_r_hash_factor = 0.02
    merge_sigma_relative = false
    merge_particles = parse(Bool, get(ENV, "MERGE_PARTICLES", "true"))
    init_Das_eta_kinematic = 0.2
    set_Das_min_kinematic_displacement = 0.01 * R
    p_correct_kuttacondition_flag = false
    wake_core_size = parse(Float64, get(ENV, "WAKE_CORE_SIZE", string(kerneloffset_targets)))
    # wake_nu_default = 1.85508e-5 / rho # from FLOWUnsteady docs
    wake_nu_default = 1.69e-5 / rho # from NASA paper
    wake_nu = parse(Float64, get(ENV, "WAKE_NU", string(wake_nu_default)))
    # wake_nu = parse(Float64, get(ENV, "WAKE_NU", "1.5e-5"))
    wake_core_beta = parse(Float64, get(ENV, "WAKE_CORE_BETA", "1.5"))
    run_kj = parse(Bool, get(ENV, "RUN_KJ", "false"))
    lamb_only = parse(Bool, get(ENV, "LAMB_ONLY", "false"))

    axial_dimension = occursin("dji9443", msh_file) ? 1 : 2 # DJI9443 geometry is rotated compared to typical rotor convention
    radial_dimension = occursin("dji9443", msh_file) ? 2 : 1 # this might be wrong for non-dji9443

    # wake_rotor = pnl.PanelWake(rotor; nwakerows=12, core_size=wake_core_size)
    wake_rotor = pnl.PanelParticleWake(rotor;
        nwakerows=1, max_particles=500_000, core_size=wake_core_size,
        particle_kerneloffset=kerneloffset_targets,
        viscous=pnl.FLOWVPM.CoreSpreading(wake_nu, wake_core_size, pnl.FLOWVPM.zeta_fmm;
            beta=wake_core_beta),
        SFS=pnl.FLOWVPM.SFS_Cd_twolevel_nobackscatter,
        # method_trailing=pnl.SigmaOverlap(R*0.05, 4.0),
        method_trailing=pnl.OverlapPPS(overlap, 2),
        method_unsteady=pnl.NoShed(),
        # method_unsteady=pnl.OverlapPPS(1.3, 2),
        unsteady_filament=false, # should be false if method_unsteady is NoShed
        shed_with_induced_velocity=false,
        particle_maintenance=pnl.ParticleMaintenance((
            pnl.GlobalCylinder([-0.5R, 0.0, 0.0], [cylinder_depth, 0.0, 0.0], 1.5R),
            pnl.MergeParticles(;
                every=merge_particles ? 1 : 0,
                r=merge_sigma_relative ? merge_r_factor : merge_r_factor * R,
                r_hash=merge_sigma_relative ? merge_r_hash_factor : merge_r_hash_factor * R,
                sigma_relative=merge_sigma_relative),
        ))
        )

    ramp_magVinf(t) = t <= ramp_t ? magVinf_start + (magVinf - magVinf_start) * (t / ramp_t) : magVinf
    Uinf(t) = ramp_magVinf(t) * Vinf_direction
    # Uinf(t) = magVinf * Vinf_direction

    smoothstep(x) = x <= 0 ? zero(x) : x >= 1 ? one(x) : x * x * (3 - 2 * x)
    function spinup_fraction(t)
        spinup_revs <= 0 && return 1.0
        return spinup_start_fraction + (1 - spinup_start_fraction) * smoothstep(t / spinup_duration)
    end

    ω_full = 2 * pi * RPM / 60
    function set_frame_omega!(frames, i_frame, ω)
        frame = frames[i_frame]
        frames[i_frame] = typeof(frame)(
            frame.x, frame.v, frame.ω_axis, ω, frame.R, frame.Rp2g,
            frame.name, frame.parent_index, frame.child_index, frame.dependent_index,
        )
    end

    frames = pnl.ReferenceFrame(rotor;
        origin=SVector{3}(0.0, 0.0, 0.0),
        v=SVector{3}(0.0, 0.0, 0.0),
        ω_axis=occursin("dji9443", msh_file) ? SVector{3}(-1.0, 0.0, 0.0) : SVector{3}(0.0, 1.0, 0.0),
        ω=ω_full * spinup_fraction(t_range[1]),
        R=SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name="vehicle",
        child_index=Int[],
        dependent_index=[1]
    )

    pnl.initialize_Das!((rotor,), frames, Uinf, t_range[1], t_range[2] - t_range[1];
        set_Das_eta_kinematic=init_Das_eta_kinematic,
        set_Das_min_kinematic_displacement)

    solver_rotor = pnl.Backslash(rotor)
    backend = pnl.FastMultipoleBackend(;
        expansion_order=parse(Int, get(ENV, "FMM_EXPANSION_ORDER", "8")),
        multipole_acceptance=parse(Float64, get(ENV, "FMM_ACCEPTANCE", "0.4")),
        leaf_size=parse(Int, get(ENV, "FMM_LEAF_SIZE", "20")),
    )
    backend_wake = pnl.FastMultipoleBackend(;
        expansion_order=parse(Int, get(ENV, "FMM_EXPANSION_ORDER", "4")),
        multipole_acceptance=parse(Float64, get(ENV, "FMM_ACCEPTANCE", "0.4")),
        leaf_size=parse(Int, get(ENV, "FMM_LEAF_SIZE", "50")),
    )
    kj_backend = pnl.FastMultipoleBackend(;
        expansion_order=parse(Int, get(ENV, "KJ_FMM_EXPANSION_ORDER", "3")),
        multipole_acceptance=parse(Float64, get(ENV, "KJ_FMM_ACCEPTANCE", "0.4")),
        leaf_size=parse(Int, get(ENV, "KJ_FMM_LEAF_SIZE", "1000")),
    )
    function maneuver!(frames, systems, wakes, t)
        set_frame_omega!(frames, 1, ω_full * spinup_fraction(t))
        return nothing
    end

    systems = (rotor,)
    wakes = (wake_rotor,)
    body_solvers = (solver_rotor,)

    pressure_bernoulli = pnl.PressureBernoulli(rho;
        unsteady=false,
        correct_kuttacondition=p_correct_kuttacondition_flag,
        backend=backend)
    force_monitor_bernoulli = pnl.ForceMonitor(length(t_range), 1;
        i_frame=1,
        normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        correct_kuttacondition=p_correct_kuttacondition_flag,
        verbose=true)

    pressure_laplace_matderiv = pnl.PressureLaplace(rotor, rho;
        acceleration_form=:material_derivative, verbose=true)
    force_monitor_laplace_matderiv = pnl.ForceMonitor(length(t_range), 1;
        i_frame=1,
        normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        correct_kuttacondition=p_correct_kuttacondition_flag,
        verbose=true)

    pressure_laplace_lamb = pnl.PressureLaplace(rotor, rho;
        acceleration_form=:lamb_vector, verbose=true)
    force_monitor_laplace_lamb = pnl.ForceMonitor(length(t_range), 1;
        i_frame=1,
        normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        correct_kuttacondition=p_correct_kuttacondition_flag,
        verbose=true)

    kj_monitor = run_kj ? pnl.KuttaJoukowskiForce(rotor, length(t_range), 1;
            rho, backend=kj_backend,
            i_frame=1,
            normalization=pnl.RotorNormalization(rho, 2 * R, 1),
            verbose=true) : nothing
    bound_circulation = pnl.BoundCirculationMonitor(rotor, length(t_range), 1;
        i_frame=1,
        radial_dimension,
        R)

    monitors = run_kj ? (
            pressure_laplace_matderiv,
            force_monitor_laplace_matderiv,
            pressure_laplace_lamb,
            force_monitor_laplace_lamb,
            pressure_bernoulli,
            force_monitor_bernoulli,
            kj_monitor,
            bound_circulation,
        ) : lamb_only ? (
            pressure_laplace_lamb,
            force_monitor_laplace_lamb,
            bound_circulation,
        ) : (
            pressure_laplace_matderiv,
            force_monitor_laplace_matderiv,
            pressure_laplace_lamb,
            force_monitor_laplace_lamb,
            pressure_bernoulli,
            force_monitor_bernoulli,
            bound_circulation,
        )

    if spinup_revs > 0
        println("\nRotor spin-up enabled: $(spinup_revs) nominal revs from $(spinup_start_fraction)×RPM to full RPM")
    end
    println("\nBegin rotor hover pressure comparison ($(length(t_range)) steps)...")
    name = run_name
    @time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
        set_Das_eta_kinematic=NaN,
        set_Das_min_kinematic_displacement,
        monitors,
        body_solvers, backend, backend_wake,
        verbose=true,
        path=save_path, name,
    )

    CT_bernoulli   = lamb_only ? fill(NaN, length(t_range)) : force_monitor_bernoulli.force[axial_dimension, :]
    CT_laplace_md  = lamb_only ? fill(NaN, length(t_range)) : force_monitor_laplace_matderiv.force[axial_dimension, :]
    CT_laplace_lv  = force_monitor_laplace_lamb.force[axial_dimension, :]
    CT_kj          = run_kj ? kj_monitor.force[axial_dimension, :] : fill(NaN, length(t_range))

    function relative_difference(a, b)
        denom = max(abs(b), eps())
        return abs(a - b) / denom
    end

    println("\nstep | CT Bernoulli | CT Laplace(∇u) | CT Laplace(λ) | CT KJ | rel(B-∇u) | rel(B-λ) | rel(B-KJ)")
    for k in 1:length(t_range)
        cb  = CT_bernoulli[k]
        cmd = CT_laplace_md[k]
        clv = CT_laplace_lv[k]
        ck  = CT_kj[k]
        println("  $k  |  $(round(cb, sigdigits=6))  |  $(round(cmd, sigdigits=6))  |  $(round(clv, sigdigits=6))  |  $(round(ck, sigdigits=6))  |  $(round(relative_difference(cb, cmd), sigdigits=4))  |  $(round(relative_difference(cb, clv), sigdigits=4))  |  $(round(relative_difference(cb, ck), sigdigits=4))")
    end

    bern_md_rel  = norm(CT_bernoulli - CT_laplace_md) / max(norm(CT_bernoulli), eps())
    bern_lv_rel  = norm(CT_bernoulli - CT_laplace_lv) / max(norm(CT_bernoulli), eps())
    bern_kj_rel  = norm(CT_bernoulli - CT_kj)         / max(norm(CT_bernoulli), eps())
    md_lv_rel    = norm(CT_laplace_md - CT_laplace_lv) / max(norm(CT_laplace_md), eps())

    println("\nRelative history differences (L2 norm):")
    println("  Bernoulli vs Laplace(∇u):  $(bern_md_rel)")
    println("  Bernoulli vs Laplace(λ):   $(bern_lv_rel)")
    println("  Bernoulli vs KJ:           $(bern_kj_rel)")
    println("  Laplace(∇u) vs Laplace(λ): $(md_lv_rel)")

    nc = rotor.ncells
    out_name = "CT_v_t_hover_RPM"*"$RPM"*"_nc$nc"*"_pps$p_per_step"*"_nt$nt"*"_overlap$overlap"*".csv"
    CT_all   = ["CT_bernoulli" "CT_laplace_md" "CT_laplace_lv" "CT_kj";
                 CT_bernoulli   CT_laplace_md   CT_laplace_lv   CT_kj] # combine different force calcs into one matrix
    writedlm(out_name, CT_all, ',')

end

## =========================================================
# SIMULATION PARAMETERS
# ==========================================================
magVinf = 0.0001    # Freestream velocity magnitude (m/s)
AOA     = 0.0       # Angle of attack (degrees)
rho     = 1.225     # Air density (kg/m^3)
RPM     = 5400      # Rotation speed (rpm)
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

# rotor_adapt_12k = build_rotor(msh_adapt_12k, te_a_12k_1, te_a_12k_2)
# rotor_adapt_24k = build_rotor(msh_adapt_24k, te_a_24k_1, te_a_24k_2)

# OpenVSP model
te_vsp_40_1 = [1614, 1574, 45] .+ 1 # (or 45 instead of 0) convert from 0-based to 1-based indexing
te_vsp_40_2 = [3324, 3284, 1755] .+ 1 # (or 1755 instead of 1711) convert from 0-based to 1-based indexing

te_vsp_56_1 = [6370, 6314, 3569] .+ 1 
te_vsp_56_2 = [3117, 3061, 316] .+ 1

te_vsp_80_1 = [6351, 6271, 590] .+ 1 
te_vsp_80_2 = [12898, 12818, 7137] .+ 1 

msh_40 = joinpath(read_path, "dji9443_new_40_40.msh")
msh_56 = joinpath(read_path, "dji9443_56_57.msh")
msh_80 = joinpath(read_path, "dji9443_80_81.msh")

rotor_vsp_40 = build_rotor(msh_40, te_vsp_40_1, te_vsp_40_2)
rotor_vsp_56 = build_rotor(msh_56, te_vsp_56_1, te_vsp_56_2)
rotor_vsp_80 = build_rotor(msh_80, te_vsp_80_1, te_vsp_80_2)


# =========================================================
# RUN SIMULATION
# ==========================================================
Vinf_direction = occursin("dji9443", msh_40) ? [cosd(AOA), sind(AOA), 0.0] : [0.0, -cosd(AOA), sind(AOA)]
Vinf = magVinf * Vinf_direction
run_simulation(rotor_vsp_40, "rotor_hover_40", msh_40; Vinf=Vinf)

Vinf_direction = occursin("dji9443", msh_56) ? [cosd(AOA), sind(AOA), 0.0] : [0.0, -cosd(AOA), sind(AOA)]
Vinf = magVinf * Vinf_direction
run_simulation(rotor_vsp_56, "rotor_hover_56", msh_56; Vinf=Vinf)

Vinf_direction = occursin("dji9443", msh_80) ? [cosd(AOA), sind(AOA), 0.0] : [0.0, -cosd(AOA), sind(AOA)]
Vinf = magVinf * Vinf_direction
run_simulation(rotor_vsp_80, "rotor_hover_80", msh_80; Vinf=Vinf)