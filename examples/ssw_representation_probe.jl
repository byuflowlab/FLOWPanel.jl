# Offline sheet-to-particle representation diagnostic for BRAINSTORM 017.
# This file deliberately uses only example-local helpers and DirectBackend.

if !isdefined(@__MODULE__, :SSWConfig)
    include(joinpath(@__DIR__, "suddenly_started_wing.jl"))
end

import Serialization

struct SSWFrozenWake <: pnl.AbstractFreeWake
    sources::Tuple
end
pnl.get_sources(wake::SSWFrozenWake) = wake.sources
pnl.get_probes(::SSWFrozenWake) = ()
pnl.reset!(::SSWFrozenWake) = nothing
pnl.apply_freestream!(::SSWFrozenWake, uinf) = nothing

function ssw_subset_wake(body, source::pnl.PanelWake, first_row::Int,
        last_row::Int)
    1 <= first_row <= last_row <= source.nwakes[] ||
        throw(ArgumentError("row range $first_row:$last_row is outside 1:$(source.nwakes[])"))
    nrows = last_row - first_row + 1
    subset = pnl.PanelWake(body; nwakerows=nrows,
        core_size=source.core_size, include_final_filament=false,
        shed_with_induced_velocity=source.shed_with_induced_velocity,
        freestream_convection=source.freestream_convection)
    for surface in eachindex(source.nodes)
        subset.nodes[surface][:, 1:nrows+1, :] .=
            source.nodes[surface][:, first_row:last_row+1, :]
        subset.strength[surface][:, 1:nrows+1, :] .=
            source.strength[surface][:, first_row:last_row+1, :]
    end
    subset.nwakes[] = nrows
    return subset
end

function _sswrp_add_segment!(pfield, r1, r2, circulation, method,
        expected_gamma, expected_impulse)
    circulation == 0 && return nothing
    direction = r2 - r1
    gamma = circulation * direction
    midpoint = 0.5 * (r1 + r2)
    expected_gamma .+= gamma
    expected_impulse .+= LA.cross(midpoint, gamma)
    pnl._shed_particles!(pfield, r1, r2, circulation, method)
    return nothing
end

"""
    ssw_convert_rows(wake, first_row, last_row, method; max_particles)

Reproduce the legacy trailing, unsteady, and open-tip deposition for an
arbitrary contiguous row range. Returns the fresh particle field and exact
vector-circulation/impulse ledgers.
"""
function ssw_convert_rows(wake::pnl.PanelWake, first_row::Int, last_row::Int,
        method::pnl.WakeSheddingMethod; max_particles::Int=2_000_000)
    1 <= first_row <= last_row <= wake.nwakes[] ||
        throw(ArgumentError("invalid conversion range $first_row:$last_row"))
    pfield = pnl.FLOWVPM.ParticleField(max_particles, Float64;
        fmm=pnl.FLOWVPM.FMM(autotune_reg_error=false))
    expected_gamma = zeros(3)
    expected_impulse = zeros(3)
    family_counts = Dict(:trailing => 0, :unsteady => 0, :tip => 0)
    # Deposit complete rings for the selected shell. Adjacent edges cancel
    # algebraically to the production trailing/unsteady/tip deposition; the two
    # range interfaces remain explicitly closed because the retained panel
    # buffer/final filament that closes them in production is absent offline.
    for row in first_row:last_row, surface in eachindex(wake.nodes)
        nodes = wake.nodes[surface]
        strength = wake.strength[surface]
        ncols = size(nodes, 3) - 1
        for column in 1:ncols
            gamma = strength[1, row, column]
            upstream_left = SVector{3}(nodes[:, row, column])
            downstream_left = SVector{3}(nodes[:, row + 1, column])
            downstream_right = SVector{3}(nodes[:, row + 1, column + 1])
            upstream_right = SVector{3}(nodes[:, row, column + 1])
            before = pfield.np
            _sswrp_add_segment!(pfield, upstream_left, downstream_left, gamma,
                method, expected_gamma, expected_impulse)
            family_counts[:trailing] += pfield.np - before
            before = pfield.np
            _sswrp_add_segment!(pfield, downstream_left, downstream_right, gamma,
                method, expected_gamma, expected_impulse)
            family_counts[:unsteady] += pfield.np - before
            before = pfield.np
            _sswrp_add_segment!(pfield, downstream_right, upstream_right, gamma,
                method, expected_gamma, expected_impulse)
            family_counts[:tip] += pfield.np - before
            _sswrp_add_segment!(pfield, upstream_right, upstream_left, gamma,
                method, expected_gamma, expected_impulse)
        end
    end
    actual_gamma = zeros(3)
    actual_impulse = zeros(3)
    for particle in 1:pfield.np
        gamma = pnl.FLOWVPM.get_Gamma(pfield, particle)
        position = pnl.FLOWVPM.get_X(pfield, particle)
        actual_gamma .+= gamma
        actual_impulse .+= LA.cross(position, gamma)
    end
    return (; pfield, expected_gamma, actual_gamma, expected_impulse,
        actual_impulse, family_counts)
end

function _sswrp_probe(points, sources)
    probes = pnl.FastMultipole.ProbeSystem(length(points), Float64)
    probes.position .= SVector{3,Float64}.(points)
    pnl.influence!((probes,), sources, pnl.DirectBackend();
        precalc=false, scalar_potential=false, velocity=true,
        velocity_gradient=true)
    return (; velocity=copy(probes.gradient), gradient=copy(probes.hessian))
end

function ssw_compare_fields(points, sheet::pnl.PanelWake, pfield;
        velocity_scale=1.0, gradient_scale=1.0)
    sheet_field = _sswrp_probe(points, pnl.get_sources(sheet))
    particle_field = _sswrp_probe(points, (pfield,))
    velocity_error = [LA.norm(a - b) / velocity_scale for
        (a, b) in zip(sheet_field.velocity, particle_field.velocity)]
    gradient_error = [LA.norm(a - b) / gradient_scale for
        (a, b) in zip(sheet_field.gradient, particle_field.gradient)]
    return (; velocity_rms=LA.norm(velocity_error) / sqrt(length(velocity_error)),
        velocity_max=maximum(velocity_error),
        gradient_rms=LA.norm(gradient_error) / sqrt(length(gradient_error)),
        gradient_max=maximum(gradient_error), sheet_field, particle_field)
end

function ssw_probe_families(sim)
    body = sim.wing
    control = [SVector{3}(body.controlpoints[:, i]) for i in 1:body.ncells]
    # A small normal offset avoids evaluating the straight sheet's singular,
    # side-dependent on-sheet limit.
    zvals = vcat(collect(range(-0.5, -0.1; length=5)),
        collect(range(0.1, 0.5; length=5)))
    chord_plane = [SVector(x, 0.0, z) for
        x in range(-0.25, 1.75; length=17),
        z in zvals]
    ymin, ymax = extrema(body.nodes[2, :])
    # Keep this line at least 15 core widths from the σ/c=0.01 validation
    # shell; closer locations measure the deliberately finite regularization
    # error instead of validating the converter's trivial limit.
    span_line = [SVector(1.25, y, 0.15) for
        y in range(ymin, ymax; length=49)]
    return (; body=vec(control), chord_plane=vec(chord_plane),
        aft_span=vec(span_line))
end

function _sswrp_write_metrics(rows, path)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "arm,selection,sigma_over_c,pps,buffer_over_c,probe,velocity_rms,velocity_max,gradient_rms,gradient_max,np,gamma_ledger_error,impulse_ledger_error")
        for row in rows
            @printf(io, "%s,%s,%.8g,%d,%.8g,%s,%.16e,%.16e,%.16e,%.16e,%d,%.16e,%.16e\n",
                row.arm, row.selection, row.sigma_over_c, row.pps, row.buffer_over_c,
                row.probe, row.velocity_rms, row.velocity_max,
                row.gradient_rms, row.gradient_max, row.np,
                row.gamma_ledger_error, row.impulse_ledger_error)
        end
    end
    return path
end

function _sswrp_complement_sources(body, wake, first_row, last_row)
    sources = ()
    if first_row > 1
        upstream = ssw_subset_wake(body, wake, 1, first_row - 1)
        sources = (sources..., pnl.get_sources(upstream)...)
    end
    if last_row < wake.nwakes[]
        downstream = ssw_subset_wake(body, wake, last_row + 1, wake.nwakes[])
        sources = (sources..., pnl.get_sources(downstream)...)
    end
    return sources
end

function ssw_frozen_resolve(snapshots, config; first_row, last_row,
        method=nothing)
    length(snapshots) == 3 ||
        throw(ArgumentError("frozen BDF2 replay requires exactly three snapshots"))
    backend = pnl.DirectBackend()
    solver = pnl.Backslash(first(snapshots).wing)
    pressure = pnl.PressureBernoulli(config.rho; unsteady=true,
        correct_kuttacondition=false, allow_partial=!isnothing(method), backend)
    force = pnl.ForceMonitor(3, 1; normalization=pnl.NoNormalization(),
        i_frame=-1, correct_kuttacondition=false, verbose=false)
    circulation = pnl.BoundCirculationMonitor(first(snapshots).wing, 3, 1;
        i_frame=1, radial_dimension=2, R=config.AR * config.c / 2,
        file=false, verbose=false)
    context = pnl.MonitorContext()
    dt = config.dt_star * config.c / config.U
    for (index, snapshot) in enumerate(snapshots)
        body = deepcopy(snapshot.wing)
        wake = snapshot.wake
        snapshot_last_row = min(last_row, wake.nwakes[])
        sources = if isnothing(method)
            pnl.get_sources(wake)
        else
            converted = ssw_convert_rows(wake, first_row, snapshot_last_row, method)
            (_sswrp_complement_sources(body, wake, first_row, snapshot_last_row)...,
                converted.pfield)
        end
        frozen = SSWFrozenWake(sources)
        frames = pnl.ReferenceFrame(body)
        pnl._steady_aerodynamics!(body, (body,), (frozen,), frames,
            config.U .* _ssw_directions(config)[1], solver;
            backend_wake=backend, backend_solve=backend,
            backend_system=backend, update_trailing_edges=false,
            grad_mu_options=SSW_GRAD_MU_OPTIONS, i_step=index - 1)
        pnl._run_monitor!(pressure, context, (body,), (frozen,), frames,
            config.U .* _ssw_directions(config)[1], index - 1, dt,
            snapshot.time_star * config.c / config.U)
        pnl._run_monitor!(force, context, (body,), (frozen,), frames,
            config.U .* _ssw_directions(config)[1], index - 1, dt,
            snapshot.time_star * config.c / config.U)
        pnl._run_monitor!(circulation, context, (body,), (frozen,), frames,
            config.U .* _ssw_directions(config)[1], index - 1, dt,
            snapshot.time_star * config.c / config.U)
    end
    qS = 0.5 * config.rho * config.U^2 * config.AR * config.c^2
    lift_hat = _ssw_directions(config)[2]
    CL = LA.dot(force.force[:, end], lift_hat) / qS
    gamma = vec(circulation.circulation_te[:, 1, end])
    return (; CL, gamma)
end

function _sswrp_write_predictions(rows, path)
    isfile(path) && throw(ArgumentError(
        "refusing to overwrite immutable Phase B prediction at $path"))
    open(path, "w") do io
        println(io, "sigma_over_c,pps,buffer_over_c,predicted_delta_CL_percent,predicted_gamma_error_percent")
        for row in rows
            @printf(io, "%.8g,%d,%.8g,%.16e,%.16e\n", row.sigma_over_c,
                row.pps, row.buffer_over_c, row.delta_CL_percent,
                row.gamma_error_percent)
        end
    end
    return path
end

function _sswrp_write_overlap_predictions(rows, path)
    isfile(path) && throw(ArgumentError(
        "refusing to overwrite immutable Phase B overlap prediction at $path"))
    open(path, "w") do io
        println(io, "sigma_over_c,overlap,buffer_over_c,predicted_delta_CL_percent,predicted_gamma_error_percent")
        for row in rows
            @printf(io, "%.8g,%.8g,%.8g,%.16e,%.16e\n", row.sigma_over_c,
                row.overlap, row.buffer_over_c, row.delta_CL_percent,
                row.gamma_error_percent)
        end
    end
    return path
end

function run_ssw_representation_probe(; trivial_only=false)
    output = get(ENV, "SSWRP_OUTPUT", joinpath("data", "ssw_representation_probe"))
    base = SSWConfig(AR=6.0, n_span=24, n_airfoil=21, dt_star=1 / 8,
        t_end_star=20.0, eta=1.0, backend_kind=:direct, save_vtk=false,
        verbose=_envbool("SSW_VERBOSE", false), output_root=output)
    sigmas = trivial_only ? (0.01, 0.005) :
        (0.05, 0.08, 0.15, 0.3, 0.6, 1.2, 1.5)
    pps_values = trivial_only ? (32, 64) : (1, 2, 4, 8)
    buffers = trivial_only ? (0.25,) : (0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
    rows = NamedTuple[]
    prediction_only = _envbool("SSWRP_PREDICTION_ONLY", false)
    for freestream_convection in (true, false)
        arm = freestream_convection ? :straight : :rolledup
        config = ssw_with(base; freestream_convection)
        state_dir = joinpath(output, "states")
        final_state = joinpath(state_dir, "$(arm)_t20p0.jls")
        if _envbool("SSWRP_REUSE_STATES", true) && isfile(final_state)
            state_paths = sort(filter(path -> startswith(basename(path),
                "$(arm)_t"), readdir(state_dir; join=true)))
            snapshots = Serialization.deserialize.(state_paths[end-2:end])
            sim = (; wing=last(snapshots).wing, wake=last(snapshots).wake)
        else
            sim = prepare_suddenly_started_wing(config)
            snapshots = NamedTuple[]
            function record_snapshot!(systems, wakes, frames, uinf, i_step, dt)
                if i_step >= length(sim.t_range) - 3
                    push!(snapshots, (; step=i_step,
                        time_star=sim.t_range[i_step + 1] * config.U / config.c,
                        wing=deepcopy(systems[1]), wake=deepcopy(wakes[1])))
                end
                return nothing
            end
            pnl.simulate!((sim.wing,), (sim.wake,), sim.frames, sim.maneuver!,
                sim.Uinf, sim.t_range; body_solvers=(sim.solver,),
                backend=sim.backend, monitors=(sim.monitors..., record_snapshot!),
                path=nothing, set_Das_eta_freestream=NaN,
                grad_mu_options=SSW_GRAD_MU_OPTIONS, verbose=false)
            mkpath(state_dir)
            for snapshot in snapshots
                token = replace(string(snapshot.time_star), "." => "p")
                Serialization.serialize(joinpath(state_dir,
                    "$(arm)_t$(token).jls"), merge(snapshot, (; config)))
            end
        end
        prediction_only && continue
        probes = ssw_probe_families(sim)
        selections = trivial_only ?
            [(:shell, first(buffers))] :
            vcat([(:shell, b) for b in buffers],
                [(:tail, b) for b in (0.25, 1.0, 4.0)])
        for (selection, buffer) in selections
            last_row = sim.wake.nwakes[]
            shell_row = clamp(round(Int, buffer / base.dt_star), 1, last_row)
            first_row = shell_row
            selected_last_row = selection == :shell ? shell_row : last_row
            sheet = ssw_subset_wake(sim.wing, sim.wake, first_row,
                selected_last_row)
            sigma_pps = trivial_only ? zip(sigmas, pps_values) :
                Iterators.product(sigmas, pps_values)
            for (sigma, pps) in sigma_pps
                converted = ssw_convert_rows(sim.wake, first_row,
                    selected_last_row,
                    pnl.SigmaPPS(sigma * base.c, pps))
                gamma_error = LA.norm(converted.actual_gamma -
                    converted.expected_gamma)
                impulse_error = LA.norm(converted.actual_impulse -
                    converted.expected_impulse)
                for (probe_name, points) in pairs(probes)
                    metrics = ssw_compare_fields(points, sheet, converted.pfield;
                        velocity_scale=base.U, gradient_scale=base.U / base.c)
                    push!(rows, (; arm, selection, sigma_over_c=sigma, pps,
                        buffer_over_c=buffer, probe=probe_name,
                        velocity_rms=metrics.velocity_rms,
                        velocity_max=metrics.velocity_max,
                        gradient_rms=metrics.gradient_rms,
                        gradient_max=metrics.gradient_max,
                        np=converted.pfield.np,
                        gamma_ledger_error=gamma_error,
                        impulse_ledger_error=impulse_error))
                end
            end
        end
    end
    path = prediction_only ? joinpath(output, "field_metrics.csv") :
        _sswrp_write_metrics(rows, joinpath(output, "field_metrics.csv"))
    if trivial_only
        maximum(max(r.velocity_max, r.gradient_max) for r in rows) < 1e-3 ||
            error("Phase A trivial-limit gate failed; see $path")
    end
    prediction_path = nothing
    overlap_prediction_path = nothing
    if !trivial_only
        prediction_path = joinpath(output, "phase_b_prediction.csv")
        overlap_prediction_path =
            joinpath(output, "phase_b_prediction_sigma_overlap.csv")
        if isfile(prediction_path) && isfile(overlap_prediction_path)
            return (; rows, path, prediction_path, overlap_prediction_path)
        end
        rolled_paths = sort(filter(path -> startswith(basename(path),
            "rolledup_t"), readdir(joinpath(output, "states"); join=true)))
        rolled = Serialization.deserialize.(rolled_paths[end-2:end])
        control = ssw_frozen_resolve(rolled, last(rolled).config;
            first_row=1, last_row=1)
        if !isfile(prediction_path)
            predictions = NamedTuple[]
            for buffer in buffers, sigma in sigmas, pps in pps_values
                row = clamp(round(Int, buffer / base.dt_star), 1,
                    last(rolled).wake.nwakes[])
                prediction = ssw_frozen_resolve(rolled, last(rolled).config;
                    first_row=row, last_row=last(rolled).wake.nwakes[],
                    method=pnl.SigmaPPS(sigma * base.c, pps))
                dcl = 100 * (prediction.CL - control.CL) / abs(control.CL)
                dgamma = 100 * maximum(abs.(prediction.gamma - control.gamma)) /
                    max(maximum(abs, control.gamma), eps(Float64))
                push!(predictions, (; sigma_over_c=sigma, pps,
                    buffer_over_c=buffer, delta_CL_percent=dcl,
                    gamma_error_percent=dgamma))
            end
            _sswrp_write_predictions(predictions, prediction_path)
        end
        if !isfile(overlap_prediction_path)
            overlap_predictions = NamedTuple[]
            overlap = 1.3
            for buffer in (0.25, 1.0, 4.0),
                    sigma in (0.08, 0.15, 0.3, 0.6, 1.2)
                row = clamp(round(Int, buffer / base.dt_star), 1,
                    last(rolled).wake.nwakes[])
                prediction = ssw_frozen_resolve(rolled, last(rolled).config;
                    first_row=row, last_row=last(rolled).wake.nwakes[],
                    method=pnl.SigmaOverlap(sigma * base.c, overlap))
                dcl = 100 * (prediction.CL - control.CL) / abs(control.CL)
                dgamma = 100 * maximum(abs.(prediction.gamma - control.gamma)) /
                    max(maximum(abs, control.gamma), eps(Float64))
                push!(overlap_predictions, (; sigma_over_c=sigma, overlap,
                    buffer_over_c=buffer, delta_CL_percent=dcl,
                    gamma_error_percent=dgamma))
            end
            _sswrp_write_overlap_predictions(overlap_predictions,
                overlap_prediction_path)
        end
    end
    return (; rows, path, prediction_path, overlap_prediction_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_ssw_representation_probe(;
        trivial_only=_envbool("SSWRP_TRIVIAL_ONLY", true))
end
