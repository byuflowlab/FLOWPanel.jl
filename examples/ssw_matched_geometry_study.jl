# Matched-geometry dynamic discriminator for BRAINSTORM 017.

if !isdefined(@__MODULE__, :_sps_case_row)
    include(joinpath(@__DIR__, "ssw_sheet_particle_split.jl"))
end

function _sswmc_csv(path)
    lines = readlines(path)
    header = split(first(lines), ',')
    return [Dict(header .=> split(line, ',')) for line in lines[2:end]]
end

_sswmc_float(row, name) = parse(Float64, row[name])
_sswmc_int(row, name) = parse(Int, row[name])
_sswmc_bool(row, name) = lowercase(row[name]) == "true"

function _sswmc_prediction_maps(root)
    pps = Dict{Tuple{Float64,Int,Float64},Tuple{Float64,Float64}}()
    for row in _sswmc_csv(joinpath(root, "phase_b_prediction.csv"))
        pps[(_sswmc_float(row, "sigma_over_c"),
            _sswmc_int(row, "pps"),
            _sswmc_float(row, "buffer_over_c"))] =
            (_sswmc_float(row, "predicted_delta_CL_percent"),
             _sswmc_float(row, "predicted_gamma_error_percent"))
    end
    overlap = Dict{Tuple{Float64,Float64,Float64},Tuple{Float64,Float64}}()
    for row in _sswmc_csv(joinpath(root,
            "phase_b_prediction_sigma_overlap.csv"))
        overlap[(_sswmc_float(row, "sigma_over_c"),
            _sswmc_float(row, "overlap"),
            _sswmc_float(row, "buffer_over_c"))] =
            (_sswmc_float(row, "predicted_delta_CL_percent"),
             _sswmc_float(row, "predicted_gamma_error_percent"))
    end
    growing = Dict{Tuple{Float64,Float64,Float64},Tuple{Float64,Float64}}()
    growing_path = joinpath(root,
        "phase_b_prediction_sigma_overlap_growing.csv")
    if isfile(growing_path)
        for row in _sswmc_csv(growing_path)
            growing[(_sswmc_float(row, "sigma_over_c"),
                _sswmc_float(row, "overlap"),
                _sswmc_float(row, "buffer_over_c"))] =
                (_sswmc_float(row, "predicted_delta_CL_percent"),
                 _sswmc_float(row, "predicted_gamma_error_percent"))
        end
    end
    return (; pps, overlap, growing)
end

function _sswmc_frozen(row, predictions)
    method = Symbol(row["method"])
    if method == :sigma_pps
        return predictions.pps[(_sswmc_float(row, "sigma_over_c"),
            _sswmc_int(row, "pps"), _sswmc_float(row, "buffer_over_c"))]
    end
    table = method == :sigma_overlap_growing ?
        predictions.growing : predictions.overlap
    return table[(_sswmc_float(row, "sigma_over_c"),
        _sswmc_float(row, "overlap"), _sswmc_float(row, "buffer_over_c"))]
end

function run_ssw_matched_geometry_study()
    split_root = get(ENV, "SSWMC_SPLIT_ROOT",
        joinpath("data", "ssw_sheet_particle_split"))
    prediction_root = get(ENV, "SSWMC_PREDICTION_ROOT",
        joinpath("data", "ssw_representation_probe"))
    output = get(ENV, "SSWMC_OUTPUT",
        joinpath(split_root, "matched_geometry"))
    rows = _sswmc_csv(joinpath(split_root, "split_summary.csv"))
    predictions = _sswmc_prediction_maps(prediction_root)

    sigma_pps = filter(r -> r["method"] == "sigma_pps", rows)
    isempty(sigma_pps) && error("matched controller requires SigmaPPS observations")
    worst_cl = sigma_pps[argmax(abs.(
        _sswmc_float.(sigma_pps, Ref("delta_CL_percent"))))]
    worst_gamma = sigma_pps[argmax(
        _sswmc_float.(sigma_pps, Ref("gamma_error_percent")))]
    admissible = filter(r -> _sswmc_bool(r, "admissible"), rows)
    finest = isempty(admissible) ? nothing :
        sort(admissible; by=r -> (_sswmc_float(r, "sigma_over_c"),
            -_sswmc_float(r, "buffer_over_c")))[1]

    roles = [(role="worst_cl", row=worst_cl),
        (role="worst_gamma", row=worst_gamma)]
    !isnothing(finest) && push!(roles, (role="finest_admissible", row=finest))
    selected = Dict{Tuple,NamedTuple}()
    for item in roles
        row = item.row
        key = (row["method"], row["sigma_over_c"], row["pps"],
            row["overlap"], row["buffer_over_c"])
        if haskey(selected, key)
            selected[key] = merge(selected[key],
                (; role=selected[key].role * ";" * item.role))
        else
            selected[key] = (; role=item.role, row)
        end
    end

    base = SSWConfig(AR=6.0, n_span=24, n_airfoil=21, dt_star=1 / 8,
        t_end_star=20.0, eta=1.0, backend_kind=:direct, save_vtk=false,
        verbose=false, output_root=split_root)
    control = load_ssw_result(ssw_with(base; wake_model=:panel))
    output_rows = NamedTuple[]
    for item in values(selected)
        row = item.row
        method = Symbol(row["method"])
        sigma = _sswmc_float(row, "sigma_over_c")
        pps = _sswmc_int(row, "pps")
        overlap = _sswmc_float(row, "overlap")
        buffer = _sswmc_float(row, "buffer_over_c")
        config = ssw_with(base; wake_model=:particle,
            panel_rows=max(1, round(Int, buffer / base.dt_star)),
            shed_method=method == :sigma_overlap_growing ?
                :sigma_overlap : method, sigma_over_c=sigma, pps_n=pps,
            pps_overlap=overlap, output_root=output)
        controller = SSWMatchedGeometry(
            config.U .* _ssw_directions(config)[1])
        maintenance = pnl.ParticleMaintenance((controller,))
        result = run_suddenly_started_wing(config;
            particle_maintenance=maintenance)
        while !(result.tail_CL.settled && result.tail_gamma.settled) &&
                config.t_end_star < 40
            config = ssw_with(config;
                t_end_star=min(40.0, config.t_end_star + 5.0))
            controller = SSWMatchedGeometry(
                config.U .* _ssw_directions(config)[1])
            maintenance = pnl.ParticleMaintenance((controller,))
            result = run_suddenly_started_wing(config;
                particle_maintenance=maintenance)
        end
        matched = _sps_case_row(result, control, method, sigma, pps,
            overlap, buffer)
        frozen_cl, frozen_gamma = _sswmc_frozen(row, predictions)
        default_cl = _sswmc_float(row, "delta_CL_percent")
        default_gamma = _sswmc_float(row, "gamma_error_percent")
        stats = ssw_matched_geometry_stats(controller)
        push!(output_rows, (; role=item.role, method, sigma_over_c=sigma,
            pps, overlap, buffer_over_c=buffer,
            default_delta_CL_percent=default_cl,
            matched_delta_CL_percent=matched.delta_CL_percent,
            frozen_delta_CL_percent=frozen_cl,
            M3b_CL_percent=default_cl - matched.delta_CL_percent,
            M3a_bound_CL_percent=matched.delta_CL_percent - frozen_cl,
            default_gamma_error_percent=default_gamma,
            matched_gamma_error_percent=matched.gamma_error_percent,
            frozen_gamma_error_percent=frozen_gamma,
            M3b_gamma_percent=default_gamma - matched.gamma_error_percent,
            M3a_bound_gamma_percent=matched.gamma_error_percent - frozen_gamma,
            settled=matched.settled,
            mean_position_correction=stats.mean_position_correction,
            max_position_correction=stats.max_position_correction,
            n_particles=matched.n_particles, tag=matched.tag))
    end
    mkpath(output)
    path = _write_sps_summary(output_rows,
        joinpath(output, "matched_geometry_summary.csv"))
    return (; output_rows, path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_ssw_matched_geometry_study()
end
