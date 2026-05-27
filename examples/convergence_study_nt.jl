## Plate, rotor, with no ground effect
# Author: Timothy Harlow
# Created: December 5, 2025

import FLOWPanel as pnl
# include(joinpath(pnl.examples_path, "helper_functions.jl"))
using FLOWPanel.FastMultipole.StaticArrays
using VSPGeom
import GeoIO
using LinearAlgebra: norm
using DelimitedFiles

function run_convergence(p_per_step, overlap, nt)
    run_name = "rotor_hover_pressure_comparison_nt$nt"
    save_path = parse(Bool, get(ENV, "SAVE_VTK", "true")) ? joinpath("data", run_name) : nothing

    magVinf = 0.0001 # + 10
    AOA = 0.0
    # rho = 1.225
    rho = 1.071778
    RPM = 5400
    R = 0.119
    shedding_r_over_R = 0.1
    nrevs = parse(Float64, get(ENV, "NREVS", "10"))
    # nt = parse(Int, get(ENV, "NT", "72"))
    dt = 60 / RPM / nt
    spinup_revs = max(0.0, parse(Float64, get(ENV, "SPINUP_REVS", "0.0")))
    spinup_start_fraction = clamp(parse(Float64, get(ENV, "SPINUP_START_FRACTION", "0.05")), eps(), 1.0)
    spinup_duration = spinup_revs * 60 / RPM
    spinup_steps = ceil(Int, spinup_duration / dt)
    n_steps = round(Int, nt * nrevs) + spinup_steps
    t_range = range(0.0, step=dt, length=n_steps)

    CPoffset = R * 1e-8
    kerneloffset_panel = parse(Float64, get(ENV, "KERNELOFFSET_PANEL", string(R * 1e-10)))
    kerneloffset_targets = parse(Float64, get(ENV, "KERNELOFFSET_TARGETS", get(ENV, "KERNELOFFSET", "3e-3")))
    kernelcutoff = R * 1e-13
    # p_per_step = 2
    # overlap = 2.0
    merge_r_factor = 0.02
    merge_r_hash_factor = 0.02
    merge_sigma_relative = false
    merge_particles = parse(Bool, get(ENV, "MERGE_PARTICLES", "true"))
    init_Das_eta_kinematic = 0.2
    set_Das_min_kinematic_displacement = 0.01 * R
    p_correct_kuttacondition_flag = false
    wake_core_size = parse(Float64, get(ENV, "WAKE_CORE_SIZE", string(kerneloffset_targets)))
    wake_nu_default = 1.85508e-5 / rho
    wake_nu = parse(Float64, get(ENV, "WAKE_NU", string(wake_nu_default)))
    # wake_nu = parse(Float64, get(ENV, "WAKE_NU", "1.5e-5"))
    wake_core_beta = parse(Float64, get(ENV, "WAKE_CORE_BETA", "1.5"))
    run_kj = parse(Bool, get(ENV, "RUN_KJ", "true"))
    lamb_only = parse(Bool, get(ENV, "LAMB_ONLY", "true"))

    read_path = joinpath(pnl.examples_path, "data")
    # msh_file = joinpath(read_path, "phantom_3_rebuild_r2.msh")
    # te_indices_1 = [9, 175, 127]
    # te_indices_2 = [13, 286, 238]
    msh_file = joinpath(read_path, "dji9443_new_40_40.msh")
    te_indices_1 = [1614, 1574, 45] .+ 1 # (or 45 instead of 0) convert from 0-based to 1-based indexing
    te_indices_2 = [3324, 3284, 1755] .+ 1 # (or 1755 instead of 1711) convert from 0-based to 1-based indexing

    axial_dimension = occursin("dji9443", msh_file) ? 1 : 2 # DJI9443 geometry is rotated compared to typical rotor convention
    radial_dimension = occursin("dji9443", msh_file) ? 2 : 1 # this might be wrong for non-dji9443

    Vinf_direction = occursin("dji9443", msh_file) ? [cosd(AOA), sind(AOA), 0.0] : [0.0, -cosd(AOA), sind(AOA)]
    Vinf = magVinf * Vinf_direction

    msh = GeoIO.load(msh_file).geometry
    nodes, cells = pnl.meshes2nodes_cells(msh)
    nodes .*= R / maximum(nodes[radial_dimension, :])

    shedding = pnl.noshedding
    kernel = Union{pnl.ConstantSource, pnl.VortexRing}
    DBC = kernel == pnl.VortexRing ? false : true
    rotor = pnl.RigidWakeBody{kernel}(nodes, cells, shedding;
        CPoffset, kerneloffset=kerneloffset_panel, kerneloffset_panel, kerneloffset_targets, kernelcutoff,
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
    # vtk_path = joinpath(save_path, "rotor_initial")
    # pnl.write_vtk(vtk_path, rotor)

    bbox1 = make_shedding_bbox(rotor.nodes, te_indices_1[1:2], radial_dimension, R, shedding_r_over_R)
    shedding1 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_1[1], te_indices_1[2];
        bbox=bbox1, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
    bbox2 = make_shedding_bbox(rotor.nodes, te_indices_2[1:2], radial_dimension, R, shedding_r_over_R)
    shedding2 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_2[1], te_indices_2[2];
        bbox=bbox2, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

    rotor = pnl.RigidWakeBody{kernel}(rotor.nodes, rotor.cells, [shedding1, shedding2];
        CPoffset, kerneloffset=kerneloffset_panel, kerneloffset_panel, kerneloffset_targets, kernelcutoff,
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

    # wake_rotor = pnl.PanelWake(rotor; nwakerows=12, core_size=wake_core_size)
    wake_rotor = pnl.PanelParticleWake(rotor;
        nwakerows=1, max_particles=500_000, core_size=wake_core_size,
        particle_kerneloffset=kerneloffset_targets,
        viscous=pnl.FLOWVPM.CoreSpreading(wake_nu, wake_core_size, pnl.FLOWVPM.zeta_fmm;
            beta=wake_core_beta),
        SFS=pnl.FLOWVPM.SFS_Cd_twolevel_nobackscatter,
        method_trailing=pnl.SigmaOverlap(R*0.05, 4.0),
        # method_trailing=pnl.OverlapPPS(1.3, 2),
        method_unsteady=pnl.NoShed(),
        # method_unsteady=pnl.OverlapPPS(1.3, 2),
        unsteady_filament=false, # should be false if method_unsteady is NoShed
        shed_with_induced_velocity=false,
        particle_maintenance=pnl.MergeParticles(;
            every=merge_particles ? 1 : 0,
            r=merge_sigma_relative ? merge_r_factor : merge_r_factor * R,
            r_hash=merge_sigma_relative ? merge_r_hash_factor : merge_r_hash_factor * R,
            sigma_relative=merge_sigma_relative)
        )

    ramp_t = 0.5 * 60 / RPM
    magVinf_start = 0.0
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
        unsteady=true,
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

    monitors = run_kj ? (
            pressure_laplace_matderiv,
            force_monitor_laplace_matderiv,
            pressure_laplace_lamb,
            force_monitor_laplace_lamb,
            pressure_bernoulli,
            force_monitor_bernoulli,
            kj_monitor,
        ) : lamb_only ? (
            pressure_laplace_lamb,
            force_monitor_laplace_lamb,
        ) : (
            pressure_laplace_matderiv,
            force_monitor_laplace_matderiv,
            pressure_laplace_lamb,
            force_monitor_laplace_lamb,
            pressure_bernoulli,
            force_monitor_bernoulli,
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
        body_solvers, backend, verbose=true,
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

    out_name = "CT_v_t_hover_RPM"*"$RPM"*"_pps$p_per_step"*"_nt$nt"*"_overlap$overlap"*".csv"
    CT_all   = ["CT_bernoulli" "CT_laplace_md" "CT_laplace_lv" "CT_kj";
                 CT_bernoulli   CT_laplace_md   CT_laplace_lv   CT_kj] # combine different force calcs into one matrix
    writedlm(out_name, CT_all, ',')
end

## =========================================================
# RUN CONVERGENCE STUDY
# ==========================================================
R = 0.119
np = 144 # number of particles per revolution

# Constant values
# nt = 36
# p_per_step = 2
overlap = 2.0

# Variable values
nts         = [36, 72, 144]
p_per_steps = Int.(np ./ nts)

for i in eachindex(p_per_steps)
    nt = nts[i]
    p_per_step = p_per_steps[i]
    println("Simulation settings:\n p_per_step: $p_per_step\n nt: $nt\n λ: $overlap")
    @show Threads.nthreads()

    run_convergence(p_per_step, overlap, nt)
end