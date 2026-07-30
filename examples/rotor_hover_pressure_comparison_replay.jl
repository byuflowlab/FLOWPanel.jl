## Rotor hover pressure comparison replay.
##
## Replays VTK output from rotor_hover_pressure_comparison.jl and reruns
## pressure/force monitors without advancing the simulation.
## Any Lamb-vector variants below are deprecated diagnostic channels only.

import FLOWPanel as pnl
using Dates
using LinearAlgebra: norm
using Printf

function env_bool(name, default)
    return parse(Bool, get(ENV, name, string(default)))
end

function env_float(name, default)
    return parse(Float64, get(ENV, name, string(default)))
end

function env_int(name, default)
    return parse(Int, get(ENV, name, string(default)))
end

function parse_recompute(value::AbstractString)
    v = lowercase(strip(value))
    v == "auto" && return (:auto,)
    v == "none" && return ()
    v == "all" && return (:all,)
    allowed = Set(["velocity", "potential", "velocity_gradient", "induced_vorticity", "p", "f"])
    parts = filter(!isempty, strip.(split(v, ",")))
    all(in(allowed), parts) || error("RECOMPUTE must be auto, none, all, or a comma-separated list from $(sort(collect(allowed)))")
    return Tuple(Symbol(p == "p" ? "P" : p == "f" ? "F" : p) for p in parts)
end

function parse_steps(value::AbstractString)
    v = lowercase(strip(value))
    v == "all" && return :all

    selected = Int[]
    for part in split(v, ",")
        token = strip(part)
        isempty(token) && continue
        if occursin(":", token)
            bounds = split(token, ":")
            length(bounds) in (2, 3) || error("Invalid STEPS range: $(token)")
            first_step = parse(Int, bounds[1])
            if length(bounds) == 2
                last_step = parse(Int, bounds[2])
                append!(selected, collect(first_step:last_step))
            else
                step = parse(Int, bounds[2])
                last_step = parse(Int, bounds[3])
                step != 0 || error("Invalid zero stride in STEPS range: $(token)")
                append!(selected, collect(first_step:step:last_step))
            end
        else
            push!(selected, parse(Int, token))
        end
    end
    isempty(selected) && error("STEPS did not select any saved step")
    return selected
end

function read_saved_steps(save_path, run_name)
    pvd_path = joinpath(save_path, run_name * "_body1.pvd")
    return pnl._read_pvd_steps(pvd_path)
end

function default_replay_steps(save_path, run_name, nsteps_replay)
    _, idxs = read_saved_steps(save_path, run_name)
    n = min(nsteps_replay, length(idxs))
    n > 0 || error("No saved steps found in $(joinpath(save_path, run_name * "_body1.pvd"))")
    return idxs[(end - n + 1):end]
end

function selected_replay_steps(save_path, run_name, nsteps_replay)
    if haskey(ENV, "STEPS")
        return parse_steps(ENV["STEPS"])
    end
    return default_replay_steps(save_path, run_name, nsteps_replay)
end

function step_description(steps)
    steps === :all && return "all"
    isempty(steps) && return "none"
    contiguous = length(steps) == 1 || all(diff(steps) .== 1)
    if contiguous
        return string(first(steps), ":", last(steps), " (", length(steps), " samples)")
    end
    return string(steps, " (", length(steps), " samples)")
end

function print_ct_table(result, histories, history_order, axial_dimension, rpm)
    names = [name for name in history_order if haskey(histories, name)]
    print("\n")
    @printf "%8s  %12s  %12s" "step" "time" "azimuth"
    for name in names
        @printf "  %16s" String(name)
    end
    print("\n")

    for (i, step) in enumerate(result.steps)
        t = result.t_range[i]
        azimuth_deg = 360 * t * rpm / 60
        @printf "%8d  %12.6g  %12.3f" step t azimuth_deg
        for name in names
            @printf "  %16.8g" histories[name][axial_dimension, i]
        end
        print("\n")
    end
end

function print_summary(histories, history_order, axial_dimension, ct_experiment)
    println("\nSummary over replayed samples:")
    for name in history_order
        haskey(histories, name) || continue
        ct = histories[name][axial_dimension, :]
        mean_ct = sum(ct) / length(ct)
        final_ct = ct[end]
        @printf "  %-16s mean CT = %.8g, final CT = %.8g" String(name) mean_ct final_ct
        if !isnan(ct_experiment)
            @printf ", |mean-exp| = %.8g, |final-exp| = %.8g" abs(mean_ct - ct_experiment) abs(final_ct - ct_experiment)
        end
        print("\n")
    end
end

mutable struct Stage1ReplayState
    rows::Vector{NamedTuple}
    seen_steps::Set{Int}
end

Stage1ReplayState() = Stage1ReplayState(NamedTuple[], Set{Int}())

struct Stage1ScaleMonitor
    state::Stage1ReplayState
    t_range::Vector{Float64}
    rho::Float64
    R::Float64
    rpm::Float64
    axial_dimension::Int
end

struct Stage1PressureSnapshot
    state::Stage1ReplayState
    label::Symbol
    t_range::Vector{Float64}
    rho::Float64
    R::Float64
    axial_dimension::Int
    correct_kuttacondition::Bool
end

struct RotorRelativeBernoulli
    rho::Float64
    reference_mode::Symbol
    correct_kuttacondition::Bool
end

pnl.monitor_provides(::RotorRelativeBernoulli) = (:P,)

struct LaplaceRHSComparison
    state::Stage1ReplayState
    raw::pnl.PressureLaplace
    surface::pnl.PressureLaplace
    t_range::Vector{Float64}
end

pnl.monitor_provides(::LaplaceRHSComparison) = ()
pnl.monitor_requires_body_hessian(::LaplaceRHSComparison) = true
pnl.monitor_requires_induced_vorticity(::LaplaceRHSComparison) = true

struct LaplacePressureComparison
    state::Stage1ReplayState
    md::pnl.PressureLaplace
    lamb::pnl.PressureLaplace
    t_range::Vector{Float64}
    rho::Float64
    R::Float64
    axial_dimension::Int
    correct_kuttacondition::Bool
end

pnl.monitor_provides(::LaplacePressureComparison) = ()
pnl.monitor_requires_body_hessian(::LaplacePressureComparison) = true
pnl.monitor_requires_induced_vorticity(::LaplacePressureComparison) = true

function vector_norm_summary(A)
    values = [norm(view(A, :, j)) for j in axes(A, 2)]
    return scalar_summary(values)
end

function scalar_summary(values)
    vals = collect(skipmissing(vec(values)))
    isempty(vals) && return (min=NaN, q05=NaN, mean=NaN, median=NaN, q95=NaN, max=NaN, rms=NaN)
    sort!(vals)
    n = length(vals)
    q(p) = vals[clamp(ceil(Int, p * n), 1, n)]
    rms = sqrt(sum(abs2, vals) / n)
    return (min=first(vals), q05=q(0.05), mean=sum(vals) / n, median=q(0.50), q95=q(0.95), max=last(vals), rms=rms)
end

function pearson_correlation(x, y)
    length(x) == length(y) || throw(ArgumentError("correlation inputs must have matching lengths"))
    isempty(x) && return NaN
    mx = sum(x) / length(x)
    my = sum(y) / length(y)
    dx = x .- mx
    dy = y .- my
    denom = sqrt(sum(abs2, dx) * sum(abs2, dy))
    denom == 0 && return NaN
    return sum(dx .* dy) / denom
end

function correlation_from_summands(x, y)
    return pearson_correlation(collect(x), collect(y))
end

function affine_fit(x, y)
    length(x) == length(y) || throw(ArgumentError("affine fit inputs must have matching lengths"))
    isempty(x) && return (scale=NaN, offset=NaN, rms_residual=NaN)
    mx = sum(x) / length(x)
    my = sum(y) / length(y)
    dx = x .- mx
    denom = sum(abs2, dx)
    scale = denom == 0 ? NaN : sum(dx .* (y .- my)) / denom
    offset = my - scale * mx
    residual = y .- (scale .* x .+ offset)
    return (scale=scale, offset=offset, rms_residual=sqrt(sum(abs2, residual) / length(residual)))
end

function push_summary!(state::Stage1ReplayState, step, time, category, quantity, units, summary; path=:state, scale=:none)
    push!(state.rows, (;
        step,
        time,
        category=String(category),
        path=String(path),
        quantity=String(quantity),
        units=String(units),
        scale=String(scale),
        min=summary.min,
        q05=summary.q05,
        mean=summary.mean,
        median=summary.median,
        q95=summary.q95,
        max=summary.max,
        rms=summary.rms,
    ))
end

function push_scalar!(state::Stage1ReplayState, step, time, category, quantity, units, value; path=:state, scale=:none)
    summary = (min=value, q05=value, mean=value, median=value, q95=value, max=value, rms=abs(value))
    push_summary!(state, step, time, category, quantity, units, summary; path, scale)
end

function push_relation!(state::Stage1ReplayState, step, time, category, quantity, values; path=:laplace_rhs_diag)
    push_summary!(state, step, time, category, quantity, "1", scalar_summary(values); path)
end

function rotor_reference_values(frames, rho, R, i_frame)
    omega = frames[i_frame].ω
    omega_r = abs(omega) * R
    q_rotor = 0.5 * rho * omega_r^2
    n = omega / (2 * pi)
    ct_ref = rho * n^2 * (2 * R)^4
    return omega, omega_r, q_rotor, ct_ref
end

function (m::RotorRelativeBernoulli)(systems, wakes, frames, uinf, i_step::Int, dt::Real)
    body = systems[1]
    fill!(body.P, zero(eltype(body.P)))

    if m.reference_mode == :zero
        pnl.calcfield_P!(body.P, body, body.velocity, 0.0, m.rho, nothing;
            correct_kuttacondition=m.correct_kuttacondition)
    elseif m.reference_mode == :kinematic
        half_rho = 0.5 * m.rho
        for i in axes(body.velocity, 2)
            u_rel2 = sum(abs2, view(body.velocity, :, i))
            u_ref2 = sum(abs2, view(body.velocity_kinematic, :, i))
            body.P[i] = half_rho * (u_ref2 - u_rel2)
        end
        if m.correct_kuttacondition && body isa pnl.AbstractLiftingBody
            for shedding in body.shedding
                for (pi, nia, nib, pj, nja, njb) in eachcol(shedding)
                    if pj != -1
                        ave = (body.P[pi] + body.P[pj]) / 2
                        body.P[pi] = ave
                        body.P[pj] = ave
                    end
                end
            end
        end
    else
        throw(ArgumentError("Unknown RotorRelativeBernoulli reference_mode $(m.reference_mode)."))
    end
    return nothing
end

function rotor_relative_label(reference_mode::Symbol)
    reference_mode == :zero && return :bernoulli_rotor_relative
    reference_mode == :kinematic && return :bernoulli_local_omega
    return Symbol(:bernoulli_, reference_mode)
end

function laplace_variant_label(gradient_mode::Symbol, acceleration_form::Symbol)
    gradient_part = gradient_mode == :raw_hessian ? :raw : :surface
    rhs_part = acceleration_form == :material_derivative ? :md : :lamb
    return Symbol(:laplace_, gradient_part, :_, rhs_part)
end

function (m::Stage1ScaleMonitor)(systems, wakes, frames, uinf, i_step::Int, dt::Real)
    step = i_step + 1
    step in m.state.seen_steps && return nothing
    push!(m.state.seen_steps, step)

    body = systems[1]
    time = m.t_range[step]
    omega, omega_r, q_rotor, ct_ref = rotor_reference_values(frames, m.rho, m.R, 1)
    q_uinf = 0.5 * m.rho * sum(abs2, uinf)

    push_scalar!(m.state, step, time, :frame, :omega, "rad/s", omega)
    push_scalar!(m.state, step, time, :frame, :omega_r, "m/s", omega_r)
    push_scalar!(m.state, step, time, :frame, :q_rotor, "Pa", q_rotor)
    push_scalar!(m.state, step, time, :frame, :q_uinf, "Pa", q_uinf)
    push_scalar!(m.state, step, time, :frame, :ct_force_reference, "N", ct_ref)
    push_scalar!(m.state, step, time, :frame, :uinf_norm, "m/s", norm(uinf))
    push_summary!(m.state, step, time, :velocity, :rotor_relative, "m/s", vector_norm_summary(body.velocity))
    push_summary!(m.state, step, time, :velocity, :kinematic, "m/s", vector_norm_summary(body.velocity_kinematic))
    push_summary!(m.state, step, time, :velocity, :inertial_estimate, "m/s", vector_norm_summary(body.velocity .+ body.velocity_kinematic))
    push_summary!(m.state, step, time, :kinematics, :angular_velocity_body, "rad/s", scalar_summary(body.angular_velocity))

    for (i, name) in enumerate(pnl.strength_names(body))
        push_summary!(m.state, step, time, :strength, Symbol(name), "solver units",
            scalar_summary(view(body.strength, :, i)))
    end
    return nothing
end

function (m::Stage1PressureSnapshot)(systems, wakes, frames, uinf, i_step::Int, dt::Real)
    body = systems[1]
    step = i_step + 1
    time = m.t_range[step]
    _, _, q_rotor, ct_ref = rotor_reference_values(frames, m.rho, m.R, 1)
    q_uinf = 0.5 * m.rho * sum(abs2, uinf)
    pressure = scalar_summary(body.P)
    force_global = pnl.calcfield_Ftot(body)
    _, R_f2g = pnl.frame_global_transform(frames, 1)
    force_local = transpose(R_f2g) * force_global
    force_axial = force_local[m.axial_dimension]
    ct_axial = force_axial / ct_ref

    push_summary!(m.state, step, time, :pressure, :gauge_pressure, "Pa", pressure;
        path=m.label, scale=:absolute)
    q_rotor > 0 && push_summary!(m.state, step, time, :pressure, :gauge_pressure_over_q_rotor, "1",
        (min=pressure.min / q_rotor, q05=pressure.q05 / q_rotor, mean=pressure.mean / q_rotor,
         median=pressure.median / q_rotor, q95=pressure.q95 / q_rotor, max=pressure.max / q_rotor,
         rms=pressure.rms / q_rotor); path=m.label, scale=:q_rotor)
    q_uinf > 0 && push_summary!(m.state, step, time, :pressure, :gauge_pressure_over_q_uinf, "1",
        (min=pressure.min / q_uinf, q05=pressure.q05 / q_uinf, mean=pressure.mean / q_uinf,
         median=pressure.median / q_uinf, q95=pressure.q95 / q_uinf, max=pressure.max / q_uinf,
         rms=pressure.rms / q_uinf); path=m.label, scale=:q_uinf)
    push_scalar!(m.state, step, time, :force, :axial_dimensional, "N", force_axial; path=m.label)
    push_scalar!(m.state, step, time, :force, :ct_axial_global, "1", ct_axial; path=m.label)

    areas = pnl.calc_areas(body)
    panel_force_global = zeros(eltype(body.P), 3, body.ncells)
    pnl.calcfield_F!(panel_force_global, body, areas, body.normals, body.P;
        correct_kuttacondition=m.correct_kuttacondition)
    panel_force_local = transpose(R_f2g) * panel_force_global
    ct_panel = vec(panel_force_local[m.axial_dimension, :] ./ ct_ref)
    ct_positive = sum(max(c, 0.0) for c in ct_panel)
    ct_negative = sum(min(c, 0.0) for c in ct_panel)
    ct_abs_sum = sum(abs, ct_panel)
    ct_net = sum(ct_panel)
    cancellation_ratio = abs(ct_net) > eps() ? ct_abs_sum / abs(ct_net) : Inf
    axial_pressure_weight = ct_panel ./ vec(body.P)
    finite_weight = isfinite.(axial_pressure_weight)

    push_scalar!(m.state, step, time, :force_cancellation, :ct_positive_sum, "1", ct_positive; path=m.label)
    push_scalar!(m.state, step, time, :force_cancellation, :ct_negative_sum, "1", ct_negative; path=m.label)
    push_scalar!(m.state, step, time, :force_cancellation, :ct_abs_sum, "1", ct_abs_sum; path=m.label)
    push_scalar!(m.state, step, time, :force_cancellation, :ct_net_panel_sum, "1", ct_net; path=m.label)
    push_scalar!(m.state, step, time, :force_cancellation, :cancellation_ratio, "1", cancellation_ratio; path=m.label)
    push_scalar!(m.state, step, time, :force_correlation, :corr_pressure_ct_panel, "1",
        pearson_correlation(vec(body.P), ct_panel); path=m.label)
    push_scalar!(m.state, step, time, :force_correlation, :corr_abs_pressure_abs_ct_panel, "1",
        pearson_correlation(abs.(vec(body.P)), abs.(ct_panel)); path=m.label)
    push_scalar!(m.state, step, time, :force_correlation, :corr_pressure_axial_weight, "1",
        any(finite_weight) ? pearson_correlation(vec(body.P)[finite_weight], axial_pressure_weight[finite_weight]) : NaN;
        path=m.label)
    return nothing
end

function raw_curl_at(body, i)
    G = body.velocity_gradient
    return (
        G[3, 2, i] - G[2, 3, i],
        G[1, 3, i] - G[3, 1, i],
        G[2, 1, i] - G[1, 2, i],
    )
end

function edge_rhs_components(pl::pnl.PressureLaplace, body, i_body)
    edges = pl.edges[i_body]
    u_inertial = pl.u_inertial[i_body]
    acceleration = pl.acceleration[i_body]
    velocity_dot = pl.velocity_dot[i_body]
    nedge = size(edges, 2)
    udot = zeros(Float64, nedge)
    md_edge = zeros(Float64, nedge)
    md_raw_accel = zeros(Float64, nedge)
    lamb_induced = zeros(Float64, nedge)
    lamb_rawcurl = zeros(Float64, nedge)
    kinetic = zeros(Float64, nedge)
    kinetic_relative = zeros(Float64, nedge)
    lamb_induced_cross = zeros(Float64, nedge)
    lamb_rawcurl_cross = zeros(Float64, nedge)

    @inbounds for k in axes(edges, 2)
        _, _, i, j = edges[:, k]
        r1 = body.controlpoints[1, j] - body.controlpoints[1, i]
        r2 = body.controlpoints[2, j] - body.controlpoints[2, i]
        r3 = body.controlpoints[3, j] - body.controlpoints[3, i]

        udot[k] = 0.5 * (
            (velocity_dot[1, i] + velocity_dot[1, j]) * r1 +
            (velocity_dot[2, i] + velocity_dot[2, j]) * r2 +
            (velocity_dot[3, i] + velocity_dot[3, j]) * r3)

        nx_i, ny_i, nz_i = body.normals[1, i], body.normals[2, i], body.normals[3, i]
        nx_j, ny_j, nz_j = body.normals[1, j], body.normals[2, j], body.normals[3, j]
        ui_n = body.velocity[1, i] * nx_i + body.velocity[2, i] * ny_i + body.velocity[3, i] * nz_i
        uj_n = body.velocity[1, j] * nx_j + body.velocity[2, j] * ny_j + body.velocity[3, j] * nz_j
        urel1 = 0.5 * (body.velocity[1, i] - ui_n * nx_i + body.velocity[1, j] - uj_n * nx_j)
        urel2 = 0.5 * (body.velocity[2, i] - ui_n * ny_i + body.velocity[2, j] - uj_n * ny_j)
        urel3 = 0.5 * (body.velocity[3, i] - ui_n * nz_i + body.velocity[3, j] - uj_n * nz_j)
        ui_t1 = body.velocity[1, i] - ui_n * nx_i
        ui_t2 = body.velocity[2, i] - ui_n * ny_i
        ui_t3 = body.velocity[3, i] - ui_n * nz_i
        uj_t1 = body.velocity[1, j] - uj_n * nx_j
        uj_t2 = body.velocity[2, j] - uj_n * ny_j
        uj_t3 = body.velocity[3, j] - uj_n * nz_j

        du1 = body.velocity[1, j] - body.velocity[1, i]
        du2 = body.velocity[2, j] - body.velocity[2, i]
        du3 = body.velocity[3, j] - body.velocity[3, i]
        md_edge[k] = urel1 * du1 + urel2 * du2 + urel3 * du3

        md_raw_accel[k] = 0.5 * (
            (acceleration[1, i] + acceleration[1, j]) * r1 +
            (acceleration[2, i] + acceleration[2, j]) * r2 +
            (acceleration[3, i] + acceleration[3, j]) * r3)

        qi = 0.5 * (u_inertial[1, i]^2 + u_inertial[2, i]^2 + u_inertial[3, i]^2)
        qj = 0.5 * (u_inertial[1, j]^2 + u_inertial[2, j]^2 + u_inertial[3, j]^2)
        kinetic[k] = qj - qi
        qi_rel = 0.5 * (ui_t1^2 + ui_t2^2 + ui_t3^2)
        qj_rel = 0.5 * (uj_t1^2 + uj_t2^2 + uj_t3^2)
        kinetic_relative[k] = qj_rel - qi_rel

        wx_i, wy_i, wz_i = body.induced_vorticity[1, i], body.induced_vorticity[2, i], body.induced_vorticity[3, i]
        wx_j, wy_j, wz_j = body.induced_vorticity[1, j], body.induced_vorticity[2, j], body.induced_vorticity[3, j]
        wx = 0.5 * (wx_i + wx_j)
        wy = 0.5 * (wy_i + wy_j)
        wz = 0.5 * (wz_i + wz_j)
        lamb1 = wy * urel3 - wz * urel2
        lamb2 = wz * urel1 - wx * urel3
        lamb3 = wx * urel2 - wy * urel1
        lamb_induced_cross[k] = lamb1 * r1 + lamb2 * r2 + lamb3 * r3

        rwx_i, rwy_i, rwz_i = raw_curl_at(body, i)
        rwx_j, rwy_j, rwz_j = raw_curl_at(body, j)
        rwx = 0.5 * (rwx_i + rwx_j)
        rwy = 0.5 * (rwy_i + rwy_j)
        rwz = 0.5 * (rwz_i + rwz_j)
        rlamb1 = rwy * urel3 - rwz * urel2
        rlamb2 = rwz * urel1 - rwx * urel3
        rlamb3 = rwx * urel2 - rwy * urel1
        lamb_rawcurl_cross[k] = rlamb1 * r1 + rlamb2 * r2 + rlamb3 * r3

        lamb_induced[k] = kinetic[k] + lamb_induced_cross[k]
        lamb_rawcurl[k] = kinetic[k] + lamb_rawcurl_cross[k]
    end

    return (; udot, md_edge, md_raw_accel, lamb_induced, lamb_rawcurl,
        kinetic, kinetic_relative, lamb_induced_cross, lamb_rawcurl_cross)
end

function log_rhs_components!(state, step, time, label, c)
    md_edge_total = c.udot .+ c.md_edge
    md_accel_convective = c.md_raw_accel .- c.udot
    lamb_induced_total = c.udot .+ c.lamb_induced
    lamb_rawcurl_total = c.udot .+ c.lamb_rawcurl

    push_summary!(state, step, time, :laplace_rhs_diag, :velocity_dot_edge, "m2/s2",
        scalar_summary(c.udot); path=label)
    push_summary!(state, step, time, :laplace_rhs_diag, :md_edge, "m2/s2",
        scalar_summary(c.md_edge); path=label)
    push_summary!(state, step, time, :laplace_rhs_diag, :md_edge_total, "m2/s2",
        scalar_summary(md_edge_total); path=label)
    push_summary!(state, step, time, :laplace_rhs_diag, :md_acceleration_total, "m2/s2",
        scalar_summary(c.md_raw_accel); path=label)
    push_summary!(state, step, time, :laplace_rhs_diag, :md_acceleration_convective, "m2/s2",
        scalar_summary(md_accel_convective); path=label)
    push_summary!(state, step, time, :laplace_rhs_diag, :lamb_induced_vorticity, "m2/s2",
        scalar_summary(c.lamb_induced); path=label)
    push_summary!(state, step, time, :laplace_rhs_diag, :lamb_induced_vorticity_total, "m2/s2",
        scalar_summary(lamb_induced_total); path=label)
    push_summary!(state, step, time, :laplace_rhs_diag, :lamb_raw_hessian_curl, "m2/s2",
        scalar_summary(c.lamb_rawcurl); path=label)
    push_summary!(state, step, time, :laplace_rhs_diag, :lamb_raw_hessian_curl_total, "m2/s2",
        scalar_summary(lamb_rawcurl_total); path=label)
    push_summary!(state, step, time, :laplace_rhs_diag, :kinetic_jump, "m2/s2",
        scalar_summary(c.kinetic); path=label)
    push_summary!(state, step, time, :laplace_rhs_diag, :kinetic_relative_jump, "m2/s2",
        scalar_summary(c.kinetic_relative); path=label)
    push_summary!(state, step, time, :laplace_rhs_diag, :lamb_induced_cross, "m2/s2",
        scalar_summary(c.lamb_induced_cross); path=label)
    push_summary!(state, step, time, :laplace_rhs_diag, :lamb_rawcurl_cross, "m2/s2",
        scalar_summary(c.lamb_rawcurl_cross); path=label)

    eps_scale = sqrt(eps(Float64))
    push_relation!(state, step, time, :laplace_rhs_relation, :rms_lamb_induced_over_md_edge,
        [scalar_summary(c.lamb_induced).rms / max(scalar_summary(c.md_edge).rms, eps_scale)]; path=label)
    push_relation!(state, step, time, :laplace_rhs_relation, :rms_lamb_rawcurl_over_md_edge,
        [scalar_summary(c.lamb_rawcurl).rms / max(scalar_summary(c.md_edge).rms, eps_scale)]; path=label)
    push_relation!(state, step, time, :laplace_rhs_relation, :rms_md_accel_conv_over_md_edge,
        [scalar_summary(md_accel_convective).rms / max(scalar_summary(c.md_edge).rms, eps_scale)]; path=label)
    push_relation!(state, step, time, :laplace_rhs_relation, :rms_md_accel_total_over_md_edge_total,
        [scalar_summary(c.md_raw_accel).rms / max(scalar_summary(md_edge_total).rms, eps_scale)]; path=label)
    push_relation!(state, step, time, :laplace_rhs_relation, :corr_lamb_induced_md_edge,
        [correlation_from_summands(c.lamb_induced, c.md_edge)]; path=label)
    push_relation!(state, step, time, :laplace_rhs_relation, :corr_lamb_rawcurl_md_edge,
        [correlation_from_summands(c.lamb_rawcurl, c.md_edge)]; path=label)
    push_relation!(state, step, time, :laplace_rhs_relation, :corr_lamb_induced_total_md_edge_total,
        [correlation_from_summands(lamb_induced_total, md_edge_total)]; path=label)
    push_relation!(state, step, time, :laplace_rhs_relation, :corr_lamb_rawcurl_total_md_edge_total,
        [correlation_from_summands(lamb_rawcurl_total, md_edge_total)]; path=label)
    return nothing
end

function (m::LaplaceRHSComparison)(systems, wakes, frames, uinf, i_step::Int, dt::Real)
    body = systems[1]
    step = i_step + 1
    time = m.t_range[step]

    pnl._pressure_fill_inertial_surface_velocity!(m.raw.u_inertial[1], body)
    pnl._pressure_velocity_dot!(m.raw, body, 1, i_step, dt)
    pnl._pressure_material_acceleration!(m.raw, body, 1)
    log_rhs_components!(m.state, step, time, :raw_hessian, edge_rhs_components(m.raw, body, 1))

    pnl._pressure_fill_inertial_surface_velocity!(m.surface.u_inertial[1], body)
    pnl._pressure_velocity_dot!(m.surface, body, 1, i_step, dt)
    pnl._pressure_material_acceleration!(m.surface, body, 1)
    log_rhs_components!(m.state, step, time, :surface_velocity, edge_rhs_components(m.surface, body, 1))

    return nothing
end

function pressure_force_ct(body, pressure, frames, axial_dimension, ct_ref, correct_kuttacondition)
    areas = pnl.calc_areas(body)
    panel_force_global = zeros(eltype(pressure), 3, body.ncells)
    pnl.calcfield_F!(panel_force_global, body, areas, body.normals, pressure;
        correct_kuttacondition)
    _, R_f2g = pnl.frame_global_transform(frames, 1)
    panel_force_local = transpose(R_f2g) * panel_force_global
    ct_panel = vec(panel_force_local[axial_dimension, :] ./ ct_ref)
    return ct_panel
end

function log_pressure_force_decomposition!(state, step, time, label, pressure, ct_panel, q_rotor)
    push_summary!(state, step, time, :pressure_compare, Symbol(label, :_pressure), "Pa",
        scalar_summary(pressure); path=:laplace_md_lamb)
    q_rotor > 0 && push_summary!(state, step, time, :pressure_compare, Symbol(label, :_pressure_over_q_rotor), "1",
        scalar_summary(pressure ./ q_rotor); path=:laplace_md_lamb)
    push_summary!(state, step, time, :pressure_force_compare, Symbol(label, :_ct_panel), "1",
        scalar_summary(ct_panel); path=:laplace_md_lamb)
    push_scalar!(state, step, time, :pressure_force_compare, Symbol(label, :_ct_net), "1",
        sum(ct_panel); path=:laplace_md_lamb)
    push_scalar!(state, step, time, :pressure_force_compare, Symbol(label, :_ct_positive_sum), "1",
        sum(max(c, 0.0) for c in ct_panel); path=:laplace_md_lamb)
    push_scalar!(state, step, time, :pressure_force_compare, Symbol(label, :_ct_negative_sum), "1",
        sum(min(c, 0.0) for c in ct_panel); path=:laplace_md_lamb)
    push_scalar!(state, step, time, :pressure_force_compare, Symbol(label, :_ct_abs_sum), "1",
        sum(abs, ct_panel); path=:laplace_md_lamb)
    return nothing
end

function log_vector_compare!(state, step, time, category, label, x, y; units="1", path=:laplace_md_lamb)
    dy = y .- x
    fit = affine_fit(x, y)
    y_affine = fit.scale .* x .+ fit.offset
    residual = y .- y_affine
    sx = scalar_summary(x)
    sy = scalar_summary(y)
    sd = scalar_summary(dy)
    sr = scalar_summary(residual)

    push_summary!(state, step, time, category, Symbol(label, :_x), units, sx; path)
    push_summary!(state, step, time, category, Symbol(label, :_y), units, sy; path)
    push_summary!(state, step, time, category, Symbol(label, :_delta), units, sd; path)
    push_summary!(state, step, time, category, Symbol(label, :_affine_residual), units, sr; path)
    push_scalar!(state, step, time, Symbol(category, :_relation), Symbol(label, :_corr), "1",
        pearson_correlation(x, y); path)
    push_scalar!(state, step, time, Symbol(category, :_relation), Symbol(label, :_affine_scale), "1",
        fit.scale; path)
    push_scalar!(state, step, time, Symbol(category, :_relation), Symbol(label, :_affine_offset), units,
        fit.offset; path)
    push_scalar!(state, step, time, Symbol(category, :_relation), Symbol(label, :_affine_residual_rms), units,
        fit.rms_residual; path)
    push_scalar!(state, step, time, Symbol(category, :_relation), Symbol(label, :_rms_delta_over_x), "1",
        sd.rms / max(sx.rms, sqrt(eps(Float64))); path)
    push_scalar!(state, step, time, Symbol(category, :_relation), Symbol(label, :_rms_residual_over_x), "1",
        sr.rms / max(sx.rms, sqrt(eps(Float64))); path)
    return (; fit, delta=dy, affine=y_affine, residual)
end

function edge_values_to_rhs(pl::pnl.PressureLaplace, body, i_body, values)
    edges = pl.edges[i_body]
    b = zeros(Float64, body.ncells)
    @inbounds for k in axes(edges, 2)
        edge_a, edge_b, i, j = edges[1, k], edges[2, k], edges[3, k], edges[4, k]
        w = pnl._pressure_edge_conormal_weight(body, edge_a, edge_b, i, j)[1]
        flux = pl.rho * w * values[k]
        b[i] += flux
        b[j] -= flux
    end
    b[pl.reference_panel] = 0.0
    return b
end

function pressure_from_rhs(pl::pnl.PressureLaplace, b)
    return pl.L[1] \ b
end

function log_rhs_component_solution!(state, step, time, label, pl, body, values,
                                     frames, axial_dimension, ct_ref,
                                     correct_kuttacondition)
    b = edge_values_to_rhs(pl, body, 1, values)
    p = pressure_from_rhs(pl, b)
    ct = pressure_force_ct(body, p, frames, axial_dimension, ct_ref,
        correct_kuttacondition)
    push_summary!(state, step, time, :rhs_component, Symbol(label, :_b), "Pa/m",
        scalar_summary(b); path=:laplace_md_lamb)
    push_summary!(state, step, time, :rhs_component, Symbol(label, :_p), "Pa",
        scalar_summary(p); path=:laplace_md_lamb)
    push_scalar!(state, step, time, :rhs_component_force, Symbol(label, :_ct_net), "1",
        sum(ct); path=:laplace_md_lamb)
    push_scalar!(state, step, time, :rhs_component_force, Symbol(label, :_ct_positive_sum), "1",
        sum(max(c, 0.0) for c in ct); path=:laplace_md_lamb)
    push_scalar!(state, step, time, :rhs_component_force, Symbol(label, :_ct_negative_sum), "1",
        sum(min(c, 0.0) for c in ct); path=:laplace_md_lamb)
    push_scalar!(state, step, time, :rhs_component_force, Symbol(label, :_ct_abs_sum), "1",
        sum(abs, ct); path=:laplace_md_lamb)
    return (; b, p, ct)
end

function (m::LaplacePressureComparison)(systems, wakes, frames, uinf, i_step::Int, dt::Real)
    body = systems[1]
    step = i_step + 1
    time = m.t_range[step]
    _, _, q_rotor, ct_ref = rotor_reference_values(frames, m.rho, m.R, 1)

    m.md(systems, wakes, frames, uinf, i_step, dt)
    p_md = copy(body.P)
    b_md = copy(m.md.b[1])

    m.lamb(systems, wakes, frames, uinf, i_step, dt)
    p_lamb = copy(body.P)
    b_lamb = copy(m.lamb.b[1])

    dp = p_lamb .- p_md
    fit = affine_fit(p_md, p_lamb)
    p_md_centered = p_md .- sum(p_md) / length(p_md)
    p_lamb_centered = p_lamb .- sum(p_lamb) / length(p_lamb)

    ct_md = pressure_force_ct(body, p_md, frames, m.axial_dimension, ct_ref,
        m.correct_kuttacondition)
    ct_lamb = pressure_force_ct(body, p_lamb, frames, m.axial_dimension, ct_ref,
        m.correct_kuttacondition)
    ct_delta = ct_lamb .- ct_md
    p_affine = fit.scale .* p_md .+ fit.offset
    p_residual = p_lamb .- p_affine
    ct_affine = pressure_force_ct(body, p_affine, frames, m.axial_dimension, ct_ref,
        m.correct_kuttacondition)
    ct_residual = pressure_force_ct(body, p_residual, frames, m.axial_dimension, ct_ref,
        m.correct_kuttacondition)

    log_pressure_force_decomposition!(m.state, step, time, :md, p_md, ct_md, q_rotor)
    log_pressure_force_decomposition!(m.state, step, time, :lamb, p_lamb, ct_lamb, q_rotor)
    log_pressure_force_decomposition!(m.state, step, time, :delta, dp, ct_delta, q_rotor)
    log_pressure_force_decomposition!(m.state, step, time, :affine, p_affine, ct_affine, q_rotor)
    log_pressure_force_decomposition!(m.state, step, time, :residual, p_residual, ct_residual, q_rotor)

    rhs_compare = log_vector_compare!(m.state, step, time, :rhs_compare, :b_lamb_from_b_md,
        b_md, b_lamb; units="Pa/m", path=:laplace_md_lamb)
    L = m.md.L[1]
    rhs_from_pressure_residual = L * p_residual
    rhs_from_pressure_delta = L * dp
    log_vector_compare!(m.state, step, time, :rhs_pressure_residual_compare,
        :rhs_residual_from_pressure, rhs_compare.residual, rhs_from_pressure_residual;
        units="Pa/m", path=:laplace_md_lamb)
    log_vector_compare!(m.state, step, time, :rhs_pressure_delta_compare,
        :rhs_delta_from_pressure, b_lamb .- b_md, rhs_from_pressure_delta;
        units="Pa/m", path=:laplace_md_lamb)

    c = edge_rhs_components(m.lamb, body, 1)
    md_component = log_rhs_component_solution!(m.state, step, time, :md_edge_component,
        m.lamb, body, c.md_edge, frames, m.axial_dimension, ct_ref,
        m.correct_kuttacondition)
    kinetic_component = log_rhs_component_solution!(m.state, step, time, :lamb_kinetic_component,
        m.lamb, body, c.kinetic, frames, m.axial_dimension, ct_ref,
        m.correct_kuttacondition)
    kinetic_relative_component = log_rhs_component_solution!(m.state, step, time,
        :relative_kinetic_component, m.lamb, body, c.kinetic_relative, frames,
        m.axial_dimension, ct_ref, m.correct_kuttacondition)
    cross_component = log_rhs_component_solution!(m.state, step, time, :lamb_cross_component,
        m.lamb, body, c.lamb_induced_cross, frames, m.axial_dimension, ct_ref,
        m.correct_kuttacondition)
    lamb_component_b = kinetic_component.b .+ cross_component.b
    log_vector_compare!(m.state, step, time, :rhs_component_consistency,
        :lamb_components_from_b_lamb, lamb_component_b, b_lamb;
        units="Pa/m", path=:laplace_md_lamb)
    log_vector_compare!(m.state, step, time, :rhs_component_consistency,
        :md_component_from_b_md, md_component.b, b_md;
        units="Pa/m", path=:laplace_md_lamb)
    log_vector_compare!(m.state, step, time, :rhs_component_compare,
        :relative_kinetic_from_md_edge, md_component.b, kinetic_relative_component.b;
        units="Pa/m", path=:laplace_md_lamb)
    log_vector_compare!(m.state, step, time, :rhs_component_compare,
        :inertial_kinetic_from_md_edge, md_component.b, kinetic_component.b;
        units="Pa/m", path=:laplace_md_lamb)

    push_scalar!(m.state, step, time, :pressure_compare_relation, :corr_p_lamb_p_md, "1",
        pearson_correlation(p_md, p_lamb); path=:laplace_md_lamb)
    push_scalar!(m.state, step, time, :pressure_compare_relation, :corr_centered_p_lamb_p_md, "1",
        pearson_correlation(p_md_centered, p_lamb_centered); path=:laplace_md_lamb)
    push_scalar!(m.state, step, time, :pressure_compare_relation, :affine_scale_lamb_from_md, "1",
        fit.scale; path=:laplace_md_lamb)
    push_scalar!(m.state, step, time, :pressure_compare_relation, :affine_offset_lamb_from_md, "Pa",
        fit.offset; path=:laplace_md_lamb)
    push_scalar!(m.state, step, time, :pressure_compare_relation, :affine_residual_rms, "Pa",
        fit.rms_residual; path=:laplace_md_lamb)
    q_rotor > 0 && push_scalar!(m.state, step, time, :pressure_compare_relation,
        :affine_residual_rms_over_q_rotor, "1", fit.rms_residual / q_rotor;
        path=:laplace_md_lamb)
    push_scalar!(m.state, step, time, :pressure_compare_relation, :rms_delta_over_md, "1",
        scalar_summary(dp).rms / max(scalar_summary(p_md).rms, sqrt(eps(Float64)));
        path=:laplace_md_lamb)
    push_scalar!(m.state, step, time, :pressure_compare_relation, :ct_delta_net, "1",
        sum(ct_delta); path=:laplace_md_lamb)
    push_scalar!(m.state, step, time, :pressure_compare_relation, :ct_delta_abs_sum, "1",
        sum(abs, ct_delta); path=:laplace_md_lamb)
    return nothing
end

function write_stage1_log(path, state::Stage1ReplayState)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "step,time,category,path,quantity,units,scale,min,q05,mean,median,q95,max,rms")
        for r in state.rows
            @printf(io, "%d,%.17g,%s,%s,%s,%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                r.step, r.time, r.category, r.path, r.quantity, r.units, r.scale,
                r.min, r.q05, r.mean, r.median, r.q95, r.max, r.rms)
        end
    end
    return path
end

function append_stage1_ct_histories!(state::Stage1ReplayState, result, histories, history_order, axial_dimension)
    for name in history_order
        haskey(histories, name) || continue
        for (i, step) in enumerate(result.steps)
            ct = histories[name][axial_dimension, i]
            push_scalar!(state, i, result.t_range[i], :force, :ct_axial_monitor, "1", ct; path=name)
        end
    end
    return state
end

# run_name = get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison")
# run_name = get(ENV, "RUN_NAME", "rotor_hover_pressure_cacceleration_form=:lamb_vectoromparison_nt36_2")
run_name = get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison_warmstart_debug")
# run_name = get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison_nt36")
# run_name = get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison_nt72_sfs")
save_path = get(ENV, "SAVE_PATH", joinpath("data", run_name))

R = env_float("R", 0.119)
# rho = env_float("RHO", 1.225)
rho = env_float("RHO", 1.071778)
rpm = env_float("RPM", 6000)
axial_dimension = env_int("AXIAL_DIMENSION", 1)
nsteps_replay = env_int("NSTEPS_REPLAY", 72)
ct_experiment = env_float("CT_EXPERIMENT", -0.072)
p_correct_kuttacondition = env_bool("P_CORRECT_KUTTA", false)

run_bernoulli = env_bool("RUN_BERNOULLI", true)
run_laplace_md = env_bool("RUN_LAPLACE_MD", false)
run_laplace_lamb = env_bool("RUN_LAPLACE_LAMB", true)
run_vortex_force = env_bool("RUN_VORTEX_FORCE", false)
run_kj = env_bool("RUN_KJ", false)
run_stage2_laplace = env_bool("RUN_STAGE2_LAPLACE", false)
run_stage3_bernoulli = env_bool("RUN_STAGE3_BERNOULLI", false)
run_stage1_log = env_bool("RUN_STAGE1_LOG", false)
run_laplace_rhs_diag = env_bool("RUN_LAPLACE_RHS_DIAG", false)
run_laplace_pressure_compare = env_bool("RUN_LAPLACE_PRESSURE_COMPARE", false)
stage1_log_path = get(ENV, "STAGE1_LOG_PATH",
    joinpath("debug", "logs", "rotor_hover_stage1_$(Dates.format(now(), "yyyymmdd_HHMMSS")).csv"))

recompute = parse_recompute(get(ENV, "RECOMPUTE", "auto"))
steps = selected_replay_steps(save_path, run_name, nsteps_replay)

backend = pnl.FastMultipoleBackend(;
    expansion_order=env_int("FMM_EXPANSION_ORDER", 8),
    multipole_acceptance=env_float("FMM_ACCEPTANCE", 0.4),
    leaf_size=env_int("FMM_LEAF_SIZE", 20),
)
kj_backend = pnl.FastMultipoleBackend(;
    expansion_order=env_int("KJ_FMM_EXPANSION_ORDER", 4),
    multipole_acceptance=env_float("KJ_FMM_ACCEPTANCE", 0.4),
    leaf_size=env_int("KJ_FMM_LEAF_SIZE", 50),
)

monitor_store = Dict{Symbol,Any}()
stage1_state = Stage1ReplayState()
history_order = [:bernoulli, :bernoulli_rotor_relative, :bernoulli_local_omega,
    :laplace_md, :laplace_lamb,
    :laplace_raw_md, :laplace_surface_md, :laplace_raw_lamb, :laplace_surface_lamb,
    :vortex_force, :kutta_joukowski]
normalization = pnl.RotorNormalization(rho, 2 * R, 1)

monitor_factory = (systems, wakes, frames, t_range) -> begin
    nt = length(t_range)
    monitors = Any[]
    empty!(monitor_store)
    run_stage1_log && push!(monitors, Stage1ScaleMonitor(stage1_state, collect(t_range), rho, R, rpm, axial_dimension))

    if run_bernoulli
        pressure = pnl.PressureBernoulli(rho;
            unsteady=false,
            correct_kuttacondition=p_correct_kuttacondition,
            backend)
        force = pnl.ForceMonitor(nt, 1;
            i_frame=1,
            normalization,
            correct_kuttacondition=p_correct_kuttacondition,
            verbose=true)
        push!(monitors, pressure, force)
        run_stage1_log && push!(monitors, Stage1PressureSnapshot(stage1_state, :bernoulli, collect(t_range), rho, R, axial_dimension, p_correct_kuttacondition))
        monitor_store[:bernoulli] = force
    end

    if run_stage3_bernoulli
        for reference_mode in (:zero, :kinematic)
            label = rotor_relative_label(reference_mode)
            pressure = RotorRelativeBernoulli(rho, reference_mode, p_correct_kuttacondition)
            force = pnl.ForceMonitor(nt, 1;
                i_frame=1,
                normalization,
                correct_kuttacondition=p_correct_kuttacondition,
                verbose=true)
            push!(monitors, pressure, force)
            run_stage1_log && push!(monitors, Stage1PressureSnapshot(stage1_state, label, collect(t_range), rho, R, axial_dimension, p_correct_kuttacondition))
            monitor_store[label] = force
        end
    end

    if run_laplace_md
        pressure = pnl.PressureLaplace(systems[1], rho;
            acceleration_form=:material_derivative,
            verbose=false)
        force = pnl.ForceMonitor(nt, 1;
            i_frame=1,
            normalization,
            correct_kuttacondition=p_correct_kuttacondition,
            verbose=true)
        push!(monitors, pressure, force)
        run_stage1_log && push!(monitors, Stage1PressureSnapshot(stage1_state, :laplace_md, collect(t_range), rho, R, axial_dimension, p_correct_kuttacondition))
        monitor_store[:laplace_md] = force
    end

    if run_laplace_lamb
        pressure = pnl.PressureLaplace(systems[1], rho;
            acceleration_form=:lamb_vector,
            verbose=false)
        force = pnl.ForceMonitor(nt, 1;
            i_frame=1,
            normalization,
            correct_kuttacondition=p_correct_kuttacondition,
            verbose=true)
        push!(monitors, pressure, force)
        run_stage1_log && push!(monitors, Stage1PressureSnapshot(stage1_state, :laplace_lamb, collect(t_range), rho, R, axial_dimension, p_correct_kuttacondition))
        monitor_store[:laplace_lamb] = force
    end

    if run_stage2_laplace
        for acceleration_form in (:material_derivative, :lamb_vector)
            for gradient_mode in (:raw_hessian, :surface_velocity)
                label = laplace_variant_label(gradient_mode, acceleration_form)
                pressure = pnl.PressureLaplace(systems[1], rho;
                    acceleration_form,
                    gradient_mode,
                    verbose=false)
                force = pnl.ForceMonitor(nt, 1;
                    i_frame=1,
                    normalization,
                    correct_kuttacondition=p_correct_kuttacondition,
                    verbose=true)
                push!(monitors, pressure, force)
                run_stage1_log && push!(monitors, Stage1PressureSnapshot(stage1_state, label, collect(t_range), rho, R, axial_dimension, p_correct_kuttacondition))
                monitor_store[label] = force
            end
        end
    end

    if run_laplace_rhs_diag
        raw = pnl.PressureLaplace(systems[1], rho;
            acceleration_form=:material_derivative,
            gradient_mode=:raw_hessian,
            verbose=false)
        surface = pnl.PressureLaplace(systems[1], rho;
            acceleration_form=:material_derivative,
            gradient_mode=:surface_velocity,
            verbose=false)
        push!(monitors, LaplaceRHSComparison(stage1_state, raw, surface, collect(t_range)))
    end

    if run_laplace_pressure_compare
        md = pnl.PressureLaplace(systems[1], rho;
            acceleration_form=:material_derivative,
            verbose=false)
        lamb = pnl.PressureLaplace(systems[1], rho;
            acceleration_form=:lamb_vector,
            verbose=false)
        push!(monitors, LaplacePressureComparison(stage1_state, md, lamb,
            collect(t_range), rho, R, axial_dimension, p_correct_kuttacondition))
    end

    if run_vortex_force
        force = pnl.SurfaceVorticityForce(systems[1], nt, 1;
            rho,
            i_frame=1,
            normalization,
            correct_kuttacondition=p_correct_kuttacondition,
            verbose=true)
        push!(monitors, force)
        monitor_store[:vortex_force] = force
    end

    if run_kj
        force = pnl.KuttaJoukowskiForce(systems[1], nt, 1;
            rho,
            backend=kj_backend,
            i_frame=1,
            normalization,
            verbose=true)
        push!(monitors, force)
        monitor_store[:kutta_joukowski] = force
    end

    isempty(monitors) && error("No replay monitors enabled; set at least one RUN_* flag to true")
    return Tuple(monitors)
end

println("Rotor hover pressure comparison replay")
println("  run_name:  $(run_name)")
println("  save_path: $(save_path)")
println("  steps:     $(step_description(steps))")
println("  recompute: $(recompute)")
println("  stage1:    $(run_stage1_log ? stage1_log_path : "disabled")")
println("  stage2:    $(run_stage2_laplace ? "PressureLaplace variants enabled" : "disabled")")
println("  stage3:    $(run_stage3_bernoulli ? "rotor-relative Bernoulli enabled" : "disabled")")
println("  rhs diag:  $(run_laplace_rhs_diag ? "enabled" : "disabled")")
println("  p compare: $(run_laplace_pressure_compare ? "enabled" : "disabled")")

result = pnl.replay(save_path, run_name;
    monitor_factory,
    recompute,
    steps,
    backend,
    backend_wake=backend,
    backend_system=backend,
    verbose=true)

histories = Dict(name => monitor.force for (name, monitor) in monitor_store)

print_ct_table(result, histories, history_order, axial_dimension, rpm)
print_summary(histories, history_order, axial_dimension, ct_experiment)

if run_stage1_log || run_laplace_rhs_diag || run_laplace_pressure_compare
    append_stage1_ct_histories!(stage1_state, result, histories, history_order, axial_dimension)
    write_stage1_log(stage1_log_path, stage1_state)
    println("\nStage 1/diagnostic log written to $(stage1_log_path)")
end
