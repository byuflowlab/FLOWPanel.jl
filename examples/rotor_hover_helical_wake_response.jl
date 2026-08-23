## Frozen-geometry strength-response sweep for the prescribed helical wake probe.
##
## Conclusive diagnostic for BRAINSTORM/002: the coupled strength relaxation in
## examples/rotor_hover_prescribed_helical_wake.jl diverges monotonically. This
## script freezes the iteration-1 helical wake GEOMETRY and measures the
## strength-response map J(Gw): the solved trailing-edge doublet jump J as a
## function of the externally imposed wake strength Gw = alpha * J0.
##
## If J(Gw) is a line of slope ~+1 lying above the diagonal J = Gw (never
## crossing it), the fixed-point iteration Gw <- relax*J + (1-relax)*Gw has NO
## finite fixed point and must diverge -- conclusively pinning the cause on the
## decoupled, non-contractive Kutta/strength relaxation rather than on missing
## hover wake physics or wake length.
##
## Run:
##   RUN_NAME=/private/tmp/helix_response SAVE_VTK=false julia --project examples/rotor_hover_helical_wake_response.jl

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
using LinearAlgebra: norm
using Printf: @sprintf

run_name  = get(ENV, "RUN_NAME", "rotor_hover_helical_wake_response")
save_path = isabspath(run_name) ? run_name : joinpath("data", run_name)
isdir(save_path) || mkpath(save_path)

save_vtk          = get(ENV, "SAVE_VTK", "true") == "true"
wake_revs         = parse(Float64, get(ENV, "WAKE_REVS", "1.25"))
rows_per_rev      = parse(Int, get(ENV, "WAKE_ROWS_PER_REV", "36"))
wake_core_size    = parse(Float64, get(ENV, "WAKE_CORE_SIZE", "1e-3"))
initial_inflow    = parse(Float64, get(ENV, "INITIAL_INFLOW", "0.08"))
near_wake_mode    = Symbol(get(ENV, "NEAR_WAKE_MODE", "te"))

# strength scale factors alpha applied to the baseline (zero-wake) TE jump J0
alpha_list = [parse(Float64, s) for s in split(get(ENV, "ALPHA_LIST",
    "0.0,0.25,0.5,0.75,1.0,1.5,2.0"), ",")]

# ----- fixed operating point / geometry (DJI9443 40_40, matches audit) -------
msh_file = joinpath(pnl.examples_path, "data", "dji9443_new_40_40.msh")
isfile(msh_file) || error("Mesh file not found: $(msh_file)")
mesh_tag = splitext(basename(msh_file))[1]

te_indices_1 = [1614, 1574, 45]   .+ 1
te_indices_2 = [3324, 3284, 1755] .+ 1

magVinf = 0.0001
AOA     = 0.0
rho     = 1.179
RPM     = parse(Float64, get(ENV, "RPM", "5400"))
R       = 0.119
shedding_r_over_R = 0.1

core_size_panel   = R * 1e-10
core_size_targets = parse(Float64, get(ENV, "CORE_SIZE", get(ENV, "KERNELOFFSET", "1e-3")))
kernelcutoff         = R * 1e-13
init_Das_eta_kinematic = 0.2
set_Das_min_kinematic_displacement = 0.01 * R

axial_dimension  = 1
radial_dimension = 2
omega_axis       = SVector{3}(-1.0, 0.0, 0.0)
default_helix_winding_sign = -sign(omega_axis[axial_dimension])
default_axial_wake_sign = -sign(omega_axis[axial_dimension])
helix_winding_sign = parse(Float64, get(ENV, "HELIX_WINDING_SIGN",
    string(default_helix_winding_sign)))
axial_wake_sign = parse(Float64, get(ENV, "WAKE_AXIAL_SIGN",
    string(default_axial_wake_sign)))
Vinf             = magVinf * [cosd(AOA), sind(AOA), 0.0]
Uinf(t)          = Vinf
omega            = 2 * pi * RPM / 60
dt               = 60 / RPM / 36
tip_speed        = abs(omega) * R

msh = pnl.read_gmsh(msh_file)
base_nodes, base_cells = pnl.meshes2nodes_cells(msh)
base_nodes .*= R / maximum(base_nodes[radial_dimension, :])

kernel = Union{pnl.ConstantSource, pnl.VortexRing}
DBC    = kernel == pnl.VortexRing ? false : true

function make_shedding_bbox(nodes, seed_nodes)
    radial_midpoint = sum(nodes[radial_dimension, seed_nodes]) / length(seed_nodes)
    radial_sign = sign(radial_midpoint)
    lower = [minimum(nodes[i, :]) for i in 1:size(nodes, 1)]
    upper = [maximum(nodes[i, :]) for i in 1:size(nodes, 1)]
    padding = max(sqrt(eps(eltype(nodes))) * R, R * 1e-6)
    lower .-= padding
    upper .+= padding
    radial_cutoff = shedding_r_over_R * R
    radial_sign > 0 ? (lower[radial_dimension] = radial_cutoff - padding) :
                      (upper[radial_dimension] = -radial_cutoff + padding)
    return (SVector{3}(lower...), SVector{3}(upper...))
end

const BASE_ROTOR = pnl.RigidWakeBody{kernel}(base_nodes, base_cells, pnl.noshedding;
    core_size=core_size_panel, core_size_panel, core_size_targets,
    kernelcutoff, semiinfinite_wake=false, watertight=true, DBC)
const SHEDDING1 = pnl.calc_shedding_from_seed(BASE_ROTOR.nodes, BASE_ROTOR.cells,
    te_indices_1[1], te_indices_1[2];
    bbox=make_shedding_bbox(BASE_ROTOR.nodes, te_indices_1[1:2]),
    normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
const SHEDDING2 = pnl.calc_shedding_from_seed(BASE_ROTOR.nodes, BASE_ROTOR.cells,
    te_indices_2[1], te_indices_2[2];
    bbox=make_shedding_bbox(BASE_ROTOR.nodes, te_indices_2[1:2]),
    normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

function build_rotor_and_frames()
    rotor = pnl.RigidWakeBody{kernel}(copy(BASE_ROTOR.nodes), copy(BASE_ROTOR.cells),
        [copy(SHEDDING1), copy(SHEDDING2)];
        core_size=core_size_panel, core_size_panel, core_size_targets,
        kernelcutoff, semiinfinite_wake=false, watertight=true,
        ensure_winding=true, DBC)

    frames = pnl.ReferenceFrame(rotor;
        origin=SVector{3}(0.0, 0.0, 0.0), v=SVector{3}(0.0, 0.0, 0.0),
        ω_axis=omega_axis, ω=omega,
        R=SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name="vehicle", child_index=Int[], dependent_index=[1])

    pnl.initialize_Das!((rotor,), frames, Uinf, 0.0, dt;
        set_Das_eta_kinematic=init_Das_eta_kinematic,
        set_Das_min_kinematic_displacement)
    return rotor, frames
end

rotor, frames = build_rotor_and_frames()
n_wake_rows = max(1, ceil(Int, wake_revs * rows_per_rev))
wake = pnl.PanelWake(rotor; nwakerows=n_wake_rows, core_size=wake_core_size,
    unsteady_filament=false)
backend = pnl.FastMultipoleBackend(;
    expansion_order=parse(Int, get(ENV, "FMM_EXPANSION_ORDER", "8")),
    multipole_acceptance=parse(Float64, get(ENV, "FMM_ACCEPTANCE", "0.4")),
    leaf_size=parse(Int, get(ENV, "FMM_LEAF_SIZE", "20")),
)

rotation_x(theta) = SMatrix{3,3}(
    1.0, 0.0,         0.0,
    0.0, cos(theta), -sin(theta),
    0.0, sin(theta),  cos(theta),
)

function te_vertices(system, shedding)
    out = zeros(eltype(system.nodes), 3, size(shedding, 2) + 1)
    for i_shed in axes(shedding, 2)
        i_panel = shedding[1, i_shed]
        idx = shedding[3, i_shed]
        out[:, i_shed] .= view(system.nodes, :, system.cells[idx, i_panel])
    end
    i_panel = shedding[1, end]
    idx = shedding[2, end]
    out[:, end] .= view(system.nodes, :, system.cells[idx, i_panel])
    return out
end

function validate_near_wake_mode(mode::Symbol)
    mode in (:te, :das_offset) ||
        error("Invalid NEAR_WAKE_MODE=$(mode). Expected te or das_offset.")
    return mode
end

function near_wake_vertices(system, i_surf; mode)
    validate_near_wake_mode(mode)
    out = te_vertices(system, system.shedding[i_surf])
    if mode == :das_offset
        out .+= system.Das[i_surf]
    end
    return out
end

function pin_near_wake_row!(wake, system; mode)
    for i_surf in eachindex(wake.nodes)
        wake.nodes[i_surf][:, 1, :] .= near_wake_vertices(system, i_surf; mode)
    end
    return nothing
end

function populate_helical_nodes!(wake, system; induced_inflow_ratio, mode)
    vi = max(induced_inflow_ratio, 0.0) * tip_speed
    for i_surf in eachindex(wake.nodes)
        seed = near_wake_vertices(system, i_surf; mode)
        nodes = wake.nodes[i_surf]
        fill!(nodes, zero(eltype(nodes)))
        for jrow in 1:(wake.nwakes[] + 1)
            age_revs = (jrow - 1) / rows_per_rev
            theta = 2 * pi * age_revs
            q = rotation_x(theta * helix_winding_sign)
            dx = axial_wake_sign * vi * theta / max(abs(omega), eps())
            for jcol in axes(nodes, 3)
                r = q * SVector{3}(seed[1, jcol], seed[2, jcol], seed[3, jcol])
                nodes[:, jrow, jcol] .= r .+ SVector{3}(dx, 0.0, 0.0)
            end
        end
    end
    return nothing
end

function candidate_strengths(system)
    out = Vector{Vector{Float64}}(undef, length(system.shedding))
    for i_surf in eachindex(system.shedding)
        out[i_surf] = [pnl._get_wakestrength_Gamma(system, ished, i_surf)
                       for ished in axes(system.shedding[i_surf], 2)]
    end
    return out
end

function strength_stats(values)
    total = 0.0; maxval = 0.0; n = 0
    for surf_values in values, gamma in surf_values
        agamma = abs(gamma)
        total += agamma
        maxval = max(maxval, agamma)
        n += 1
    end
    return (mean = n == 0 ? 0.0 : total / n, max = maxval)
end

# Impose wake strength = scale * target[i_surf][ished] on EVERY age-row of each
# shedding column (the uniform-sheet convention used by the relaxation probe).
function impose_strengths!(wake, target; scale)
    for i_surf in eachindex(wake.strength)
        str = wake.strength[i_surf]
        for ished in axes(str, 3)
            value = scale * target[i_surf][ished]
            str[1, :, ished] .= value
            if size(str, 1) > 1
                str[2:end, :, ished] .= zero(eltype(str))
            end
        end
    end
    return nothing
end

# Mean SIGNED wake-induced axial velocity at the blade control points, isolated
# by a dedicated wake->body influence pass (sign reveals downwash vs upwash).
function wake_induced_axial_on_body(system, wake)
    pnl.reset!(system)
    pnl.reset!(wake)
    pnl.influence!((system,), (wake,), backend; precalc=true,
        scalar_potential=false, velocity=true, velocity_gradient=false)
    total = 0.0; n = 0
    for i in axes(system.velocity, 2)
        total += system.velocity[axial_dimension, i]
        n += 1
    end
    return n == 0 ? NaN : total / n
end

function run_steady_with_explicit_wake!(system, wake, frames)
    pressure = pnl.PressureBernoulli(rho;
        unsteady=false, correct_kuttacondition=false, backend=backend, file=false)
    force = pnl.ForceMonitor(1, 1; i_frame=1,
        normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        correct_kuttacondition=false, verbose=false, file=false)
    vorticity_force = pnl.SurfaceVorticityForce(system, 1, 1; rho=rho, i_frame=1,
        normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        correct_kuttacondition=false, verbose=false, file=false)
    kj = pnl.KuttaJoukowskiForce(system, 1, 1; rho=rho, backend=backend,
        i_frame=1, normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        verbose=false, file=false)
    monitors = (pressure, force, vorticity_force, kj)

    systems_tuple = (system,)
    wakes_tuple = (wake,)
    pnl.audit_monitors(monitors)
    for sys in systems_tuple
        sys.needs_velocity_gradient[] = any(pnl.monitor_requires_body_hessian, monitors)
        pnl.calc_normals!(sys)
        pnl.calc_controlpoints!(sys)
    end

    pnl._steady_aerodynamics!(systems_tuple, systems_tuple, wakes_tuple, frames, Vinf,
        (pnl.Backslash(system),);
        backend_wake=backend, backend_solve=backend, backend_system=backend,
        needs_induced_vorticity=any(pnl.monitor_requires_induced_vorticity, monitors),
        update_trailing_edges=false)

    ctx = pnl.MonitorContext()
    pnl.monitor_set_time!(ctx, 0.0)
    for monitor in monitors
        pnl._run_monitor!(monitor, ctx, systems_tuple, wakes_tuple, frames, Vinf, 0, 1.0)
    end

    return (
        CT_bernoulli = -force.force[axial_dimension, 1],
        CT_vorticity = -vorticity_force.force[axial_dimension, 1],
        CT_kj = -kj.force[axial_dimension, 1],
    )
end

function main()
    wake.nwakes[] = n_wake_rows
    wake.overflowed[] = true
    # freeze the iteration-1 helical geometry for the whole sweep
    populate_helical_nodes!(wake, rotor; induced_inflow_ratio=initial_inflow,
        mode=near_wake_mode)
    pin_near_wake_row!(wake, rotor; mode=near_wake_mode)

    println("\nFrozen-geometry strength-response sweep - $(mesh_tag), $(RPM) RPM")
    println(@sprintf("  wake %.2f revs, %d rows/rev, %d panel rows (geometry frozen)",
        wake_revs, rows_per_rev, wake.nwakes[]))
    println(@sprintf("  near wake: %s, helix winding sign %.1f, axial wake sign %.1f",
        string(near_wake_mode), helix_winding_sign, axial_wake_sign))
    println("  References: steady rigid wake CT ~= 0.0505; VPM ~= 0.062; experiment ~= 0.072")

    # baseline: zero wake strength -> J0 (the TE jump with no shed sheet)
    fill!.(wake.strength, 0.0)
    base_result = run_steady_with_explicit_wake!(rotor, wake, frames)
    J0 = candidate_strengths(rotor)
    J0_stats = strength_stats(J0)
    println(@sprintf("\n  Baseline (Gw=0): CT_bern=%.5f  meanJ0=%.5e\n", base_result.CT_bernoulli, J0_stats.mean))

    println("  alpha   imposed_Gw   solved_J     J - Gw       CT(Bern)  CT(vort)  CT(KJ)   wake_uax")
    println("  " * "-"^96)

    rows = NamedTuple[]
    for alpha in alpha_list
        impose_strengths!(wake, J0; scale=alpha)
        Gw_stats = strength_stats([[alpha * g for g in s] for s in J0])
        uax = wake_induced_axial_on_body(rotor, wake)  # measured before solve, on frozen wake
        # re-impose (reset! inside the velocity probe cleared system fields only)
        impose_strengths!(wake, J0; scale=alpha)
        result = run_steady_with_explicit_wake!(rotor, wake, frames)
        J = candidate_strengths(rotor)
        J_stats = strength_stats(J)
        diff = J_stats.mean - Gw_stats.mean
        push!(rows, (alpha=alpha, Gw=Gw_stats.mean, J=J_stats.mean, diff=diff,
            CT_bernoulli=result.CT_bernoulli, CT_vorticity=result.CT_vorticity,
            CT_kj=result.CT_kj, wake_uax=uax))
        println(@sprintf("  %5.2f  %10.5e  %10.5e  %10.5e  %8.5f  %8.5f  %8.5f  %+.4e",
            alpha, Gw_stats.mean, J_stats.mean, diff,
            result.CT_bernoulli, result.CT_vorticity, result.CT_kj, uax))
    end

    # least-squares slope/offset of J vs Gw
    n = length(rows)
    sx = sum(r.Gw for r in rows); sy = sum(r.J for r in rows)
    sxx = sum(r.Gw^2 for r in rows); sxy = sum(r.Gw * r.J for r in rows)
    denom = n * sxx - sx^2
    slope = denom == 0 ? NaN : (n * sxy - sx * sy) / denom
    offset = denom == 0 ? NaN : (sy - slope * sx) / n
    crosses = any(r.J <= r.Gw + 1e-12 for r in rows if r.Gw > 0)

    csv_path = joinpath(save_path, "response_sweep.csv")
    open(csv_path, "w") do io
        println(io, "near_wake_mode,helix_winding_sign,axial_wake_sign,alpha,imposed_Gw,solved_J,J_minus_Gw,CT_bernoulli,CT_vorticity,CT_kj,wake_induced_axial")
        for r in rows
            println(io, @sprintf("%s,%.1f,%.1f,%.4f,%.8e,%.8e,%.8e,%.8f,%.8f,%.8f,%.8e",
                string(near_wake_mode), helix_winding_sign, axial_wake_sign,
                r.alpha, r.Gw, r.J, r.diff, r.CT_bernoulli, r.CT_vorticity, r.CT_kj, r.wake_uax))
        end
    end
    println("\nWrote $(csv_path)")

    println("\nDiagnosis:")
    println(@sprintf("  J(Gw) least-squares fit: slope = %.4f, offset = %.5e", slope, offset))
    if !crosses && slope > 0.5
        println("  J(Gw) is a positive-slope line that NEVER crosses the diagonal J = Gw:")
        println("  the strength relaxation has NO finite fixed point and MUST diverge.")
        println("  => Cause = non-contractive decoupled Kutta/strength relaxation (positive feedback),")
        println("     NOT missing hover wake physics or wake length.")
    elseif crosses
        println("  J(Gw) crosses the diagonal => a finite fixed point exists; divergence")
        println("  would then point elsewhere (geometry relaxation, wake length, etc.).")
    else
        println("  Inconclusive slope; inspect response_sweep.csv.")
    end
end

main()
