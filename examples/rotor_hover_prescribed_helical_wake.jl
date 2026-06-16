## Bagai/Leishman-style finite free-wake relaxation probe for DJI9443 hover.
##
## This is a diagnostic harness, not a supported wake model. It reuses the
## DJI9443 setup from rotor_hover_force_method_audit.jl, fills a PanelWake with
## a finite helical sheet as the initial condition, then pseudo-implicitly
## relaxes the wake geometry and TE-jump strengths together.
##
## Quick smoke:
##   RUN_NAME=/private/tmp/bl_freewake_smoke SAVE_VTK=false MAX_ITER=2 WAKE_REVS=0.25 WAKE_ROWS_PER_REV=12 julia --project examples/rotor_hover_prescribed_helical_wake.jl
##
## Intended diagnostic:
##   RUN_NAME=/private/tmp/bl_freewake_default SAVE_VTK=false MAX_ITER=20 julia --project examples/rotor_hover_prescribed_helical_wake.jl

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
import GeoIO
using LinearAlgebra: norm
using Printf: @sprintf

run_name  = get(ENV, "RUN_NAME", "rotor_hover_prescribed_helical_wake")
save_path = isabspath(run_name) ? run_name : joinpath("data", run_name)
isdir(save_path) || mkpath(save_path)

save_vtk          = get(ENV, "SAVE_VTK", "true") == "true"
wake_revs         = parse(Float64, get(ENV, "WAKE_REVS", "1.25"))
rows_per_rev      = parse(Int, get(ENV, "WAKE_ROWS_PER_REV", "36"))
max_iter          = parse(Int, get(ENV, "MAX_ITER", "20"))
wake_relax        = parse(Float64, get(ENV, "WAKE_RELAX", "0.15"))
strength_relax    = parse(Float64, get(ENV, "STRENGTH_RELAX", "0.25"))
pseudo_dt_factor  = parse(Float64, get(ENV, "WAKE_PSEUDO_DT_FACTOR", "1.0"))
max_node_step_R   = parse(Float64, get(ENV, "MAX_NODE_STEP_R", "0.03"))
gamma_tol         = parse(Float64, get(ENV, "GAMMA_TOL", "5e-3"))
node_tol_R        = parse(Float64, get(ENV, "NODE_TOL_R", "2e-3"))
wake_core_size    = parse(Float64, get(ENV, "WAKE_CORE_SIZE", "1e-3"))
initial_inflow    = parse(Float64, get(ENV, "INITIAL_INFLOW", "0.08"))
ct_tol            = parse(Float64, get(ENV, "CT_TOL", "2e-4"))
sample_r_min      = parse(Float64, get(ENV, "SAMPLE_R_MIN", "0.25"))
sample_r_max      = parse(Float64, get(ENV, "SAMPLE_R_MAX", "0.95"))
near_wake_mode    = Symbol(get(ENV, "NEAR_WAKE_MODE", "te"))

# ----- fixed operating point / geometry (DJI9443 40_40, matches audit) -------
msh_file = joinpath(pnl.examples_path, "data", "dji9443_new_40_40.msh")
isfile(msh_file) || error("Mesh file not found: $(msh_file)")
mesh_tag = splitext(basename(msh_file))[1]

# 0-based ParaView TE seed point IDs (see rotor_hover_pressure_comparison.jl)
te_indices_1 = [1614, 1574, 45]   .+ 1
te_indices_2 = [3324, 3284, 1755] .+ 1

magVinf = 0.0001
AOA     = 0.0
rho     = 1.179
RPM     = parse(Float64, get(ENV, "RPM", "5400"))
R       = 0.119
shedding_r_over_R = 0.1

kerneloffset_panel   = R * 1e-10
kerneloffset_targets = parse(Float64, get(ENV, "KERNELOFFSET", "1e-3"))
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
wake_row_dt      = 60 / RPM / rows_per_rev
pseudo_dt        = pseudo_dt_factor * wake_row_dt
tip_speed        = abs(omega) * R

msh = GeoIO.load(msh_file).geometry
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
    kerneloffset=kerneloffset_panel, kerneloffset_panel, kerneloffset_targets,
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
        kerneloffset=kerneloffset_panel, kerneloffset_panel, kerneloffset_targets,
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
        idx = shedding[3, i_shed] # nib
        out[:, i_shed] .= view(system.nodes, :, system.cells[idx, i_panel])
    end
    i_panel = shedding[1, end]
    idx = shedding[2, end] # final nia
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
    total = 0.0
    maxval = 0.0
    n = 0
    for surf_values in values, gamma in surf_values
        agamma = abs(gamma)
        total += agamma
        maxval = max(maxval, agamma)
        n += 1
    end
    return (mean = n == 0 ? 0.0 : total / n, max = maxval)
end

function relax_strengths!(wake, candidates; relax)
    gamma_num = 0.0
    gamma_den = 0.0
    relaxed = Vector{Vector{Float64}}(undef, length(wake.strength))
    for i_surf in eachindex(wake.strength)
        str = wake.strength[i_surf]
        relaxed[i_surf] = similar(candidates[i_surf])
        for ished in axes(str, 3)
            current = str[1, 1, ished]
            candidate = candidates[i_surf][ished]
            next_gamma = (1 - relax) * current + relax * candidate
            relaxed[i_surf][ished] = next_gamma
            gamma_num += abs(next_gamma - current)
            gamma_den += abs(current)
            str[1, :, ished] .= next_gamma
            if size(str, 1) > 1
                str[2:end, :, ished] .= zero(eltype(str))
            end
        end
    end
    return relaxed, gamma_num / max(gamma_den, eps())
end

function wake_row_inflow_ratio(wake)
    total = 0.0
    n = 0
    for (isurf, nodes) in enumerate(wake.nodes)
        vel = wake.velocity[isurf]
        for jrow in 2:(wake.nwakes[] + 1), jcol in axes(nodes, 3)
            y = nodes[2, jrow, jcol]
            z = nodes[3, jrow, jcol]
            rr = sqrt(y^2 + z^2) / R
            if sample_r_min <= rr <= sample_r_max
                induced_axial = abs(vel[axial_dimension, jrow, jcol] - Vinf[axial_dimension])
                total += induced_axial
                n += 1
            end
        end
    end
    return n == 0 ? NaN : max(total / n / tip_speed, 0.0)
end

function relax_wake_nodes!(wake, system; relax, pseudo_dt, max_step)
    pin_near_wake_row!(wake, system; mode=near_wake_mode)
    total = 0.0
    maxdisp = 0.0
    n = 0
    for i_surf in eachindex(wake.nodes)
        nodes = wake.nodes[i_surf]
        vel = wake.velocity[i_surf]
        for jrow in 2:(wake.nwakes[] + 1), jcol in axes(nodes, 3)
            # A steady wake has nonzero convection velocity. Relax toward the
            # age-marched position from the upstream row so node motion has a
            # fixed point, instead of adding U*dt to every node indefinitely.
            upstream = SVector{3}(nodes[1, jrow - 1, jcol],
                                  nodes[2, jrow - 1, jcol],
                                  nodes[3, jrow - 1, jcol])
            current = SVector{3}(nodes[1, jrow, jcol],
                                 nodes[2, jrow, jcol],
                                 nodes[3, jrow, jcol])
            u_up = SVector{3}(vel[1, jrow - 1, jcol],
                              vel[2, jrow - 1, jcol],
                              vel[3, jrow - 1, jcol])
            u_here = SVector{3}(vel[1, jrow, jcol],
                                vel[2, jrow, jcol],
                                vel[3, jrow, jcol])
            target = upstream + 0.5 * (u_up + u_here) * pseudo_dt
            residual = target - current
            residual_norm = norm(residual)
            capped = residual_norm > max_step ? residual * (max_step / residual_norm) : residual
            disp = relax * capped
            nodes[:, jrow, jcol] .+= disp
            disp_norm = norm(disp)
            total += disp_norm
            maxdisp = max(maxdisp, disp_norm)
            n += 1
        end
    end
    pin_near_wake_row!(wake, system; mode=near_wake_mode)
    wake.nwakes[] = n_wake_rows
    wake.overflowed[] = true
    return (mean_R = n == 0 ? 0.0 : total / n / R, max_R = maxdisp / R)
end

function run_steady_with_explicit_wake!(system, wake, frames, iteration; write_state::Bool)
    pressure = pnl.PressureBernoulli(rho;
        unsteady=false, correct_kuttacondition=false, backend=backend, file=false)
    force = pnl.ForceMonitor(1, 1; i_frame=1,
        normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        correct_kuttacondition=false, verbose=false, file=false)
    pressure_laplace_lamb = pnl.PressureLaplace((system,), rho;
        unsteady=false, acceleration_form=:lamb_vector, file=false)
    force_laplace_lamb = pnl.ForceMonitor(1, 1; i_frame=1,
        normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        correct_kuttacondition=false, verbose=false, file=false)
    vorticity_force = pnl.SurfaceVorticityForce(system, 1, 1; rho=rho, i_frame=1,
        normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        correct_kuttacondition=false, verbose=false, file=false)
    kj = pnl.KuttaJoukowskiForce(system, 1, 1; rho=rho, backend=backend,
        i_frame=1, normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        verbose=false, file=false)
    monitors = (pressure, force, pressure_laplace_lamb, force_laplace_lamb,
        vorticity_force, kj)

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
    pnl.monitor_set_time!(ctx, Float64(iteration - 1))
    for monitor in monitors
        pnl._run_monitor!(monitor, ctx, systems_tuple, wakes_tuple, frames, Vinf, 0, 1.0)
    end

    if write_state
        tag = "$(basename(run_name))_iter$(lpad(iteration, 2, '0'))"
        pnl.write_vtk(joinpath(save_path, tag * "_body1"), system, iteration - 1, iteration - 1;
            monitors=monitors, i_system=1, overwrite=iteration == 1, compress=true)
        pnl.write_vtk(joinpath(save_path, tag * "_wake"), wake, iteration - 1, iteration - 1;
            overwrite=iteration == 1, compress=true)
    end

    return (
        CT_bernoulli = -force.force[axial_dimension, 1],
        CT_laplace_lamb = -force_laplace_lamb.force[axial_dimension, 1],
        CT_vorticity = -vorticity_force.force[axial_dimension, 1],
        CT_kj = -kj.force[axial_dimension, 1],
        wake_row_inflow = wake_row_inflow_ratio(wake),
    )
end

function write_iteration_table(rows)
    csv_path = joinpath(save_path, "iteration_table.csv")
    open(csv_path, "w") do io
        println(io, "iteration,wake_revs,rows_per_rev,nwake_panels,CT_bernoulli,CT_laplace_lamb,CT_vorticity,CT_kj,delta_CT,mean_abs_gamma,max_abs_gamma,mean_abs_gamma_candidate,max_abs_gamma_candidate,delta_gamma_rel,mean_node_step_R,max_node_step_R,wake_row_inflow_ratio,converged,stopping_reason")
        for r in rows
            println(io, @sprintf("%d,%.8f,%d,%d,%.8f,%.8f,%.8f,%.8f,%.8f,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8f,%s,%s",
                r.iteration, r.wake_revs, r.rows_per_rev, r.nwake_panels,
                r.CT_bernoulli, r.CT_laplace_lamb, r.CT_vorticity, r.CT_kj,
                r.delta_CT,
                r.mean_abs_gamma, r.max_abs_gamma,
                r.mean_abs_gamma_candidate, r.max_abs_gamma_candidate,
                r.delta_gamma_rel, r.mean_node_step_R, r.max_node_step_R,
                r.wake_row_inflow_ratio, string(r.converged), r.stopping_reason))
        end
    end
    return csv_path
end

function main()
    wake.nwakes[] = n_wake_rows
    wake.overflowed[] = true
    populate_helical_nodes!(wake, rotor; induced_inflow_ratio=initial_inflow,
        mode=near_wake_mode)
    fill!.(wake.strength, 0.0)
    pin_near_wake_row!(wake, rotor; mode=near_wake_mode)

    println("\nBagai/Leishman-style PanelWake relaxation probe - $(mesh_tag), $(RPM) RPM")
    println(@sprintf("  wake %.2f revs, %d rows/rev, %d panel rows, core %.3e m",
        wake_revs, rows_per_rev, wake.nwakes[], wake_core_size))
    println(@sprintf("  near wake: %s, helix winding sign %.1f, axial wake sign %.1f",
        string(near_wake_mode), helix_winding_sign, axial_wake_sign))
    println(@sprintf("  relax: wake %.3f, strength %.3f, pseudo dt %.3e s, max node step %.4f R",
        wake_relax, strength_relax, pseudo_dt, max_node_step_R))
    println("  References: steady rigid wake CT ~= 0.0505; VPM ~= 0.062; BEM ~= 0.068; experiment ~= 0.072")
    println("\n  iter  CT(Bern)  CT(Lamb)  CT(vort)  CT(KJ)    dCT      dGamma   nodeMax/R  gammaMean  candMean  rowInflow  status")
    println("  " * "-"^123)

    rows = NamedTuple[]
    previous_ct = NaN
    previous_growth = (ct = NaN, gamma = NaN, node = NaN)
    unstable_streak = 0
    converged = false
    stopping_reason = "max_iter"

    for iter in 1:max_iter
        wake.nwakes[] = n_wake_rows
        result = run_steady_with_explicit_wake!(rotor, wake, frames, iter; write_state=save_vtk)
        candidates = candidate_strengths(rotor)
        candidate_stats = strength_stats(candidates)
        relaxed_strengths, delta_gamma_rel = relax_strengths!(wake, candidates; relax=strength_relax)
        relaxed_stats = strength_stats(relaxed_strengths)
        node_stats = relax_wake_nodes!(wake, rotor;
            relax=wake_relax, pseudo_dt=pseudo_dt, max_step=max_node_step_R * R)

        delta_ct = isfinite(previous_ct) ? abs(result.CT_bernoulli - previous_ct) : Inf
        row_inflow = result.wake_row_inflow

        converged = iter > 1 &&
            delta_ct < ct_tol &&
            delta_gamma_rel < gamma_tol &&
            node_stats.max_R < node_tol_R
        if converged
            stopping_reason = "converged"
        end

        ct_abs = abs(result.CT_bernoulli)
        grows_ct = isfinite(previous_growth.ct) && ct_abs > 1.35 * max(previous_growth.ct, eps())
        grows_gamma = isfinite(previous_growth.gamma) && relaxed_stats.mean > 1.35 * max(previous_growth.gamma, eps())
        grows_node = isfinite(previous_growth.node) && node_stats.max_R > 1.35 * max(previous_growth.node, eps())
        if grows_ct || grows_gamma || grows_node || !isfinite(result.CT_bernoulli) || !isfinite(delta_gamma_rel)
            unstable_streak += 1
        else
            unstable_streak = 0
        end
        unstable = unstable_streak >= 3
        if unstable
            stopping_reason = "unstable_or_divergent"
        end

        push!(rows, (
            iteration=iter, wake_revs=wake_revs, rows_per_rev=rows_per_rev,
            nwake_panels=sum(size(s, 3) for s in wake.strength) * wake.nwakes[],
            CT_bernoulli=result.CT_bernoulli,
            CT_laplace_lamb=result.CT_laplace_lamb,
            CT_vorticity=result.CT_vorticity, CT_kj=result.CT_kj,
            delta_CT=delta_ct,
            mean_abs_gamma=relaxed_stats.mean, max_abs_gamma=relaxed_stats.max,
            mean_abs_gamma_candidate=candidate_stats.mean,
            max_abs_gamma_candidate=candidate_stats.max,
            delta_gamma_rel=delta_gamma_rel,
            mean_node_step_R=node_stats.mean_R, max_node_step_R=node_stats.max_R,
            wake_row_inflow_ratio=row_inflow,
            converged=converged, stopping_reason=stopping_reason,
        ))

        status = converged ? "converged" : (unstable ? "unstable" : "relaxing")
        println(@sprintf("  %4d  %8.5f  %8.5f  %8.5f  %8.5f  %7.4f  %8.3e  %9.4f  %9.3e  %8.3e  %9.4f  %s",
            iter, result.CT_bernoulli, result.CT_laplace_lamb,
            result.CT_vorticity, result.CT_kj,
            isfinite(delta_ct) ? delta_ct : NaN, delta_gamma_rel,
            node_stats.max_R, relaxed_stats.mean, candidate_stats.mean,
            isfinite(row_inflow) ? row_inflow : NaN, status))

        previous_ct = result.CT_bernoulli
        previous_growth = (ct = ct_abs, gamma = relaxed_stats.mean, node = node_stats.max_R)
        (converged || unstable) && break
    end

    if !isempty(rows) && rows[end].stopping_reason == "max_iter" && length(rows) >= max_iter
        rows[end] = merge(rows[end], (stopping_reason="max_iter",))
    end

    csv_path = write_iteration_table(rows)
    last = rows[end]
    println("\nWrote $(csv_path)")
    if save_vtk
        println("Wrote body/wake VTK snapshots under $(save_path)/")
    end

    println("\nVerdict:")
    if last.stopping_reason == "unstable_or_divergent"
        println(@sprintf("  Unstable/divergent: CT %.5f, delta_gamma %.3e, max node step %.4f R.",
            last.CT_bernoulli, last.delta_gamma_rel, last.max_node_step_R))
    elseif last.converged && 0.060 <= last.CT_bernoulli <= 0.080
        println(@sprintf("  Converged at CT %.5f (Laplace Lamb %.5f), in the target VPM/BEM/experiment neighborhood.",
            last.CT_bernoulli, last.CT_laplace_lamb))
    elseif last.converged && last.CT_bernoulli < 0.060
        println(@sprintf("  Converged at CT %.5f (Laplace Lamb %.5f), below the target VPM/BEM/experiment neighborhood.",
            last.CT_bernoulli, last.CT_laplace_lamb))
    elseif last.converged
        println(@sprintf("  Converged at CT %.5f (Laplace Lamb %.5f), above the target VPM/BEM/experiment neighborhood.",
            last.CT_bernoulli, last.CT_laplace_lamb))
    else
        println(@sprintf("  Not converged within MAX_ITER=%d: CT %.5f (Laplace Lamb %.5f), delta_gamma %.3e, max node step %.4f R.",
            max_iter, last.CT_bernoulli, last.CT_laplace_lamb,
            last.delta_gamma_rel, last.max_node_step_R))
    end
end

main()
