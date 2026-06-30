## Rotor hover: compare PressureBernoulli and PressureLaplace.
## KuttaJoukowskiForce is available as an opt-in diagnostic with RUN_KJ=true.

import FLOWPanel as pnl
# include(joinpath(pnl.examples_path, "helper_functions.jl"))
using FLOWPanel.FastMultipole.StaticArrays
using VSPGeom
import GeoIO
using LinearAlgebra: norm

run_name = get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison")
save_path = parse(Bool, get(ENV, "SAVE_VTK", "true")) ? joinpath("data", run_name) : nothing

AOA = 0.0
# rho = 1.225
# rho = 1.071778 # from FLOWUnsteady docs
rho = 1.179 # from NASA paper
RPM = parse(Float64, get(ENV, "RPM", "6000"))
R = 0.119
shedding_r_over_R = 0.1
nrevs = parse(Float64, get(ENV, "NREVS", "10"))
nt = parse(Int, get(ENV, "NT", "36"))
dt = 60 / RPM / nt
spinup_revs = max(0.0, parse(Float64, get(ENV, "SPINUP_REVS", "0.0")))
spinup_start_fraction = clamp(parse(Float64, get(ENV, "SPINUP_START_FRACTION", "0.05")), eps(), 1.0)
spinup_duration = spinup_revs * 60 / RPM
spinup_steps = ceil(Int, spinup_duration / dt)
n_steps = round(Int, nt * nrevs) + spinup_steps
t_range = range(0.0, step=dt, length=n_steps)

# --- Item 005: staged-startup convecting-freestream pulse ---------------------
# A four-phase smoothstep profile for magVinf(t), measured in rotor revolutions:
#   ramp-up   : magVinf_start -> magVinf_peak  over FREESTREAM_RAMP_REVS
#   hold      : magVinf_peak                    for FREESTREAM_HOLD_REVS
#   withdraw  : magVinf_peak -> magVinf_end     over FREESTREAM_WITHDRAW_REVS
#   hover     : magVinf_end                     for the remainder
# magVinf_peak defaults to ~2x the hover induced velocity v_h = sqrt(T/(2 rho A))
# (with CT~0.06, R=0.119, RPM=6000: v_h ~ 4.6 m/s) so the freestream initially
# dominates self-induction and sweeps the shed wake clear, then is slowly
# withdrawn so wake self-induction sustains the column into hover.
sec_per_rev = 60 / RPM
magVinf_start  = parse(Float64, get(ENV, "MAGVINF_START",  "0.0"))
magVinf_peak   = parse(Float64, get(ENV, "MAGVINF_PEAK",   "10.0"))
magVinf_end    = parse(Float64, get(ENV, "MAGVINF_END",    "0.0"))
freestream_ramp_revs     = parse(Float64, get(ENV, "FREESTREAM_RAMP_REVS",     "2.0"))
freestream_hold_revs     = parse(Float64, get(ENV, "FREESTREAM_HOLD_REVS",     "3.0"))
freestream_withdraw_revs = parse(Float64, get(ENV, "FREESTREAM_WITHDRAW_REVS", "4.0"))
settle_revs              = parse(Float64, get(ENV, "SETTLE_REVS",              "4.0"))
t_ramp_up   = freestream_ramp_revs     * sec_per_rev
t_hold      = freestream_hold_revs     * sec_per_rev
t_withdraw  = freestream_withdraw_revs * sec_per_rev
cylinder_depth = parse(Float64, get(ENV, "TRUNCATION_DEPTH_R", "4")) * R

# Ensure the run is long enough to cover ramp-up + hold + withdraw + a settle
# tail (so it never ends mid-withdrawal and the plateau ripple is observable).
schedule_revs = freestream_ramp_revs + freestream_hold_revs +
    freestream_withdraw_revs + settle_revs
required_revs = max(nrevs, schedule_revs)
n_steps = round(Int, nt * required_revs) + spinup_steps
t_range = range(0.0, step=dt, length=n_steps)

cp_outer=true
kerneloffset_panel = parse(Float64, get(ENV, "KERNELOFFSET_PANEL", string(R * 1e-10)))
kerneloffset_targets = parse(Float64, get(ENV, "KERNELOFFSET_TARGETS", get(ENV, "KERNELOFFSET", "1e-3")))
kernelcutoff = R * 1e-13
p_per_step = parse(Int, get(ENV, "P_PER_STEP", "2"))
overlap = parse(Float64, get(ENV, "OVERLAP", "3.0"))
particle_shedding = lowercase(get(ENV, "PARTICLE_SHEDDING", "overlap_pps"))
merge_r_factor = parse(Float64, get(ENV, "MERGE_R_FACTOR", "0.02"))
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
# Item 005 E4: multiplicative bump to the viscous core-spreading rate for the
# damping sweep (WAKE_NU_FACTOR=3, 10, ...); default 1.0 leaves wake_nu unchanged.
wake_nu *= parse(Float64, get(ENV, "WAKE_NU_FACTOR", "1.0"))
# wake_nu = parse(Float64, get(ENV, "WAKE_NU", "1.5e-5"))
wake_core_beta = parse(Float64, get(ENV, "WAKE_CORE_BETA", "1.5"))
run_kj = parse(Bool, get(ENV, "RUN_KJ", "false"))
lamb_only = parse(Bool, get(ENV, "LAMB_ONLY", "false"))
# Item 005 tuning: steady-Bernoulli-only monitor set. Skips both CG pressure
# solves and the per-step FMM Hessian, so iterations are far cheaper while still
# giving a CT history to read plateau ripple from.
bernoulli_only = parse(Bool, get(ENV, "BERNOULLI_ONLY", "false"))
run_monitors = parse(Bool, get(ENV, "RUN_MONITORS", "true"))

read_path = joinpath(pnl.examples_path, "data")

# msh_file = joinpath(read_path, "phantom_3_rebuild_r2.msh")
# te_indices_1 = [9, 175, 127]
# te_indices_2 = [13, 286, 238]

rhpc_mesh = lowercase(get(ENV, "RHPC_MESH", "80_81"))
if rhpc_mesh == "40_40"
    msh_file = joinpath(read_path, "dji9443_new_40_40.msh")
    te_indices_1 = [1614, 1574, 45] .+ 1 # (or 45 instead of 0) convert from 0-based to 1-based indexing
    te_indices_2 = [3324, 3284, 1755] .+ 1 # (or 1755 instead of 1711) convert from 0-based to 1-based indexing
elseif rhpc_mesh == "56_57"
    msh_file = joinpath(read_path, "dji9443_56_57.msh")
    te_indices_1 = [6370, 6314, 3255] .+ 1 # (or 3255 instead of 3253) convert from 0-based to 1-based indexing
    te_indices_2 = [3117, 3061, 0] .+ 1 # (or 0 instead of 1) convert from 0-based to 1-based indexing
elseif rhpc_mesh == "80_81"
    msh_file = joinpath(read_path, "dji9443_80_81.msh")
    te_indices_1 = [12898, 12818, 6549] .+ 1 # (or 6549 instead of 6547) convert from 0-based to 1-based indexing
    te_indices_2 = [6351, 6271, 3] .+ 1 # (or 0 instead of 3) convert from 0-based to 1-based indexing
else
    error("Unknown RHPC_MESH=$(repr(rhpc_mesh)); use 40_40, 56_57, or 80_81")
end

axial_dimension = occursin("dji9443", msh_file) ? 1 : 2 # DJI9443 geometry is rotated compared to typical rotor convention
radial_dimension = occursin("dji9443", msh_file) ? 2 : 1 # this might be wrong for non-dji9443

Vinf_direction = occursin("dji9443", msh_file) ? [cosd(AOA), sind(AOA), 0.0] : [0.0, -cosd(AOA), sind(AOA)]

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
# println("Saving initial at $(save_path)...")
# vtk_path = joinpath(save_path, "rotor_initial")
# pnl.write_vtk(vtk_path, rotor)
# sherlock

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

# Diagnostic gate: ablate the clipping_backscatter strategy by replacing the
# canonical SFS with an equivalent DynamicSFS that omits the clip. Default
# off, leaves the run bit-identical to the canonical configuration.
sfs_no_backscatter_clip    = parse(Bool,    get(ENV, "SFS_NO_BACKSCATTER_CLIP",    "false"))
sfs_no_backscatter_project = parse(Bool,    get(ENV, "SFS_NO_BACKSCATTER_PROJECT", "false"))
sfs_backscatter_signed     = parse(Bool,    get(ENV, "SFS_BACKSCATTER_SIGNED",     "false"))
sfs_magnitude_control      = parse(Bool,    get(ENV, "SFS_MAGNITUDE_CONTROL",      "false"))
sfs_directional_control    = parse(Bool,    get(ENV, "SFS_DIRECTIONAL_CONTROL",    "false"))
sfs_threelevel             = parse(Bool,    get(ENV, "SFS_THREELEVEL",             "false"))
sfs_nostatic               = parse(Bool,    get(ENV, "SFS_NOSTATIC",               "false"))
sfs_maxC                   = parse(Float64, get(ENV, "SFS_MAXC",                   "1.0"))
sfs_rlxf                   = parse(Float64, get(ENV, "SFS_RLXF",                   "0.005"))
panel_wake_hessian_to_particles = parse(Bool, get(ENV, "PANEL_WAKE_HESSIAN_TO_PARTICLES",
    string(!parse(Bool, get(ENV, "WAKEROW_NO_HESSIAN_TO_PARTICLES", "true")))))
wakerow_no_hessian_to_particles = !panel_wake_hessian_to_particles
body_hessian_to_particles       = parse(Bool, get(ENV, "BODY_HESSIAN_TO_PARTICLES",        "false"))
body_on_wake                    = parse(Bool, get(ENV, "BODY_ON_WAKE",                     "true"))
panel_wake_on_particles         = parse(Bool, get(ENV, "PANEL_WAKE_VELOCITY_TO_PARTICLES",
    get(ENV, "PANEL_WAKE_ON_PARTICLES", "true")))
particle_hessian_self           = parse(Bool, get(ENV, "PARTICLE_HESSIAN_SELF",            "true"))
particle_relax                  = parse(Bool, get(ENV, "PARTICLE_RELAX",                   "true"))
bound_strength_rlx              = parse(Float64, get(ENV, "BOUND_STRENGTH_RLX",            "1.0"))  # E4.8 body-strength low-pass (1.0=off)
diagnose_particle_gamma         = parse(Bool, get(ENV, "DIAGNOSE_PARTICLE_GAMMA",          "false"))
diagnose_particle_influence     = parse(Bool, get(ENV, "DIAGNOSE_PARTICLE_INFLUENCE",      "false"))
particle_diagnostic_vertical    = ntuple(i -> i == axial_dimension ? 1.0 : 0.0, 3)
sfs_off                         = parse(Bool, get(ENV, "SFS_OFF",                          "false"))

# Item 006: optional spatially-filtered relaxation. When RELAX_FILTER_DOWNSTREAM_R is
# set, apply Pedrizzetti relaxation only to particles that have propagated at least
# RELAX_FILTER_DOWNSTREAM_R*R downstream (+axial) of the rotor plane, leaving the
# near-rotor band unrelaxed. Unset/NaN => unfiltered full-wake relaxation (default).
relax_filter_downstream_R = parse(Float64, get(ENV, "RELAX_FILTER_DOWNSTREAM_R", "NaN"))
# Relaxation factor (rlxf) of the corrected-Pedrizzetti scheme. Defaults to the FLOWVPM
# stock value so unset behavior is unchanged; override with RELAX_RLXF.
stock_relaxation = pnl.FLOWVPM.relaxation_correctedpedrizzetti
relax_rlxf = parse(Float64, get(ENV, "RELAX_RLXF", string(stock_relaxation.rlxf)))
base_relaxation = pnl.FLOWVPM.Relaxation(stock_relaxation.relax,
    stock_relaxation.nsteps_relax, relax_rlxf)
relaxation_scheme = if isnan(relax_filter_downstream_R)
    base_relaxation
else
    d = relax_filter_downstream_R * R
    plane_point  = SVector{3,Float64}(ntuple(i -> i == axial_dimension ? d   : 0.0, 3))
    plane_normal = SVector{3,Float64}(ntuple(i -> i == axial_dimension ? 1.0 : 0.0, 3))
    pnl.plane_filtered_relaxation(base_relaxation, plane_point, plane_normal; i_frame=1)
end

sfs_choice = if sfs_off
    pnl.FLOWVPM.noSFS
elseif sfs_backscatter_signed
    pnl.FLOWVPM.SFS_Cd_twolevel_backscatter_signed
elseif sfs_no_backscatter_project
    pnl.FLOWVPM.SFS_Cd_twolevel_nobackscatter_projection
elseif sfs_threelevel
    pnl.FLOWVPM.SFS_Cd_threelevel_nobackscatter
else
    sfs_controls = ()
    sfs_magnitude_control   && (sfs_controls = (sfs_controls..., pnl.FLOWVPM.control_magnitude))
    sfs_directional_control && (sfs_controls = (sfs_controls..., pnl.FLOWVPM.control_directional))
    sfs_clippings = sfs_no_backscatter_clip ? () : (pnl.FLOWVPM.clipping_backscatter,)
    sfs_model = sfs_nostatic ?
        ((pfield; optargs...) -> pnl.FLOWVPM.E_nostaticparticles(pfield; E=pnl.FLOWVPM.Estr_fmm, optargs...)) :
        pnl.FLOWVPM.Estr_fmm
    pnl.FLOWVPM.DynamicSFS(sfs_model,
        pnl.FLOWVPM.pseudo3level_beforeUJ,
        pnl.FLOWVPM.pseudo3level_positive_afterUJ;
        alpha=0.999, clippings=sfs_clippings, controls=sfs_controls,
        maxC=sfs_maxC, rlxf=sfs_rlxf)
end

# wake_rotor = pnl.PanelWake(rotor; nwakerows=12, core_size=wake_core_size)
FV = pnl.FLOWVPM

method_trailing = if particle_shedding == "overlap_pps"
    pnl.OverlapPPS(overlap, p_per_step)
elseif particle_shedding == "sigma_overlap"
    tip_sigma = 2 * pi * R / nt * overlap / p_per_step
    pnl.SigmaOverlap(tip_sigma, overlap)
else
    error("Unknown PARTICLE_SHEDDING=$(repr(particle_shedding)); use overlap_pps or sigma_overlap")
end

wake_rotor = pnl.PanelParticleWake(rotor;
    nwakerows=1, max_particles=500_000, core_size=wake_core_size,
    particle_kerneloffset=kerneloffset_targets,
    viscous=pnl.FLOWVPM.CoreSpreading(wake_nu, wake_core_size, pnl.FLOWVPM.zeta_fmm;
        beta=wake_core_beta),
    SFS=sfs_choice,
    relaxation=relaxation_scheme,
    method_trailing,
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

smoothstep(x) = x <= 0 ? zero(x) : x >= 1 ? one(x) : x * x * (3 - 2 * x)

# Item 005: four-phase smoothstep freestream pulse (see parameter block above).
function magVinf_pulse(t)
    if t <= t_ramp_up
        return magVinf_start + (magVinf_peak - magVinf_start) * smoothstep(t_ramp_up > 0 ? t / t_ramp_up : 1.0)
    elseif t <= t_ramp_up + t_hold
        return magVinf_peak
    elseif t <= t_ramp_up + t_hold + t_withdraw
        s = t_withdraw > 0 ? (t - t_ramp_up - t_hold) / t_withdraw : 1.0
        return magVinf_peak + (magVinf_end - magVinf_peak) * smoothstep(s)
    else
        return magVinf_end
    end
end
Uinf(t) = magVinf_pulse(t) * Vinf_direction

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

monitors = !run_monitors ? () : run_kj ? (
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
    ) : bernoulli_only ? (
        pressure_bernoulli,
        force_monitor_bernoulli,
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
let
    r0 = freestream_ramp_revs
    r1 = r0 + freestream_hold_revs
    r2 = r1 + freestream_withdraw_revs
    println("\nFreestream pulse (revs): ramp-up [0, $(round(r0,digits=2))] -> peak=$(magVinf_peak); " *
        "hold [$(round(r0,digits=2)), $(round(r1,digits=2))]; " *
        "withdraw [$(round(r1,digits=2)), $(round(r2,digits=2))] -> end=$(magVinf_end); " *
        "hover/settle [$(round(r2,digits=2)), $(round(required_revs,digits=2))]")
    println("Total run length: $(round(required_revs,digits=2)) revs ($(length(t_range)) steps), truncation depth=$(round(cylinder_depth/R,digits=2))R")
end
println("\nBegin rotor hover pressure comparison ($(length(t_range)) steps)...")
println("Particle diagnostics: PARTICLE_SHEDDING=$(particle_shedding), RUN_MONITORS=$(run_monitors), BODY_HESSIAN_TO_PARTICLES=$(body_hessian_to_particles), PANEL_WAKE_HESSIAN_TO_PARTICLES=$(panel_wake_hessian_to_particles), PANEL_WAKE_VELOCITY_TO_PARTICLES=$(panel_wake_on_particles), PARTICLE_HESSIAN_SELF=$(particle_hessian_self), PARTICLE_RELAX=$(particle_relax), DIAGNOSE_PARTICLE_GAMMA=$(diagnose_particle_gamma), DIAGNOSE_PARTICLE_INFLUENCE=$(diagnose_particle_influence), diagnostic_vertical=$(particle_diagnostic_vertical)")
name = run_name

# Allow other scripts to `include` this file purely for setup (geometry,
# frames, wake, monitors) without executing the time-marching call.
rhpc_setup_only = parse(Bool, get(ENV, "RHPC_SETUP_ONLY", "false"))

if !rhpc_setup_only

@time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
    set_Das_eta_kinematic=NaN,
    set_Das_min_kinematic_displacement,
    monitors,
    body_solvers, backend, backend_wake,
    wakerow_no_hessian_to_particles,
    body_hessian_to_particles,
    body_on_wake,
    panel_wake_on_particles,
    particle_hessian_self,
    particle_relax,
    bound_strength_rlx,
    diagnose_particle_gamma,
    diagnose_particle_influence,
    diagnostic_vertical=particle_diagnostic_vertical,
    verbose=true,
    path=save_path, name,
)

if !run_monitors
    println("\nRUN_MONITORS=false; skipping pressure/force history summary.")
else

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

# Item 005 primary metric: peak-to-peak CT ripple over the final settle window
# (after the freestream is fully withdrawn), with the residual freestream so a
# plateau under nonzero convection is not mistaken for hover.
let
    settle_window_revs = min(settle_revs, required_revs)
    k_start = max(1, length(t_range) - round(Int, nt * settle_window_revs) + 1)
    tail_b  = filter(isfinite, CT_bernoulli[k_start:end])
    tail_md = filter(isfinite, CT_laplace_md[k_start:end])
    tail_lv = filter(isfinite, CT_laplace_lv[k_start:end])
    ptp(v) = isempty(v) ? NaN : maximum(v) - minimum(v)
    mean(v) = isempty(v) ? NaN : sum(v) / length(v)
    residual_magVinf = magVinf_pulse(t_range[end])
    println("\nItem 005 plateau diagnostics (final $(round(settle_window_revs,digits=2)) revs, steps $(k_start):$(length(t_range))):")
    println("  residual magVinf at readout = $(residual_magVinf)  (hover requires ≈ 0)")
    println("  CT Bernoulli   plateau mean=$(round(mean(tail_b), sigdigits=5))  peak-to-peak=$(round(ptp(tail_b), sigdigits=4))")
    println("  CT Laplace(∇u) plateau mean=$(round(mean(tail_md),sigdigits=5))  peak-to-peak=$(round(ptp(tail_md),sigdigits=4))")
    println("  CT Laplace(λ)  plateau mean=$(round(mean(tail_lv),sigdigits=5))  peak-to-peak=$(round(ptp(tail_lv),sigdigits=4))")
end

if save_path !== nothing
    isdir(save_path) || mkpath(save_path)
    csv_path = joinpath(save_path, "$(run_name)_CT_vs_rev.csv")
    open(csv_path, "w") do io
        println(io, "step,revolution,CT_bernoulli,CT_laplace_matderiv,CT_laplace_lamb,CT_kj")
        for k in 1:length(t_range)
            rev = (k - 1) * dt * RPM / 60
            println(io, "$k,$rev,$(-CT_bernoulli[k]),$(-CT_laplace_md[k]),$(-CT_laplace_lv[k]),$(-CT_kj[k])")
        end
    end
    println("\nWrote CT vs revolution CSV: $csv_path")
end

end # run_monitors

end # !rhpc_setup_only
