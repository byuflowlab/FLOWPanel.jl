## Rotor hover pressure comparison replay.
##
## Replays VTK output from rotor_hover_pressure_comparison.jl and reruns
## pressure/force monitors without advancing the simulation.

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
# run_name = get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison_nt36")
run_name = get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison_nt72_sfs")
save_path = get(ENV, "SAVE_PATH", joinpath("data", run_name))

R = env_float("R", 0.119)
# rho = env_float("RHO", 1.225)
rho = env_float("RHO", 1.071778)
rpm = env_float("RPM", 6000)
axial_dimension = env_int("AXIAL_DIMENSION", 1)
nsteps_replay = env_int("NSTEPS_REPLAY", 3)
ct_experiment = env_float("CT_EXPERIMENT", -0.072)
p_correct_kuttacondition = env_bool("P_CORRECT_KUTTA", false)

run_bernoulli = env_bool("RUN_BERNOULLI", true)
run_laplace_md = env_bool("RUN_LAPLACE_MD", true)
run_laplace_lamb = env_bool("RUN_LAPLACE_LAMB", true)
run_kj = env_bool("RUN_KJ", true)
run_stage2_laplace = env_bool("RUN_STAGE2_LAPLACE", false)
run_stage3_bernoulli = env_bool("RUN_STAGE3_BERNOULLI", false)
run_stage1_log = env_bool("RUN_STAGE1_LOG", false)
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
    :kutta_joukowski]
normalization = pnl.RotorNormalization(rho, 2 * R, 1)

monitor_factory = (systems, wakes, frames, t_range) -> begin
    nt = length(t_range)
    monitors = Any[]
    empty!(monitor_store)
    run_stage1_log && push!(monitors, Stage1ScaleMonitor(stage1_state, collect(t_range), rho, R, rpm, axial_dimension))

    if run_bernoulli
        pressure = pnl.PressureBernoulli(rho;
            unsteady=true,
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

if run_stage1_log
    append_stage1_ct_histories!(stage1_state, result, histories, history_order, axial_dimension)
    write_stage1_log(stage1_log_path, stage1_state)
    println("\nStage 1 scale log written to $(stage1_log_path)")
end
