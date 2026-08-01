# Conditional frozen prediction for the h ∝ sigma^1.5 refinement path.

if !isdefined(@__MODULE__, :ssw_frozen_resolve)
    include(joinpath(@__DIR__, "ssw_representation_probe.jl"))
end

function _sswgo_read_csv(path)
    lines = readlines(path)
    header = split(first(lines), ',')
    return [Dict(header .=> split(line, ',')) for line in lines[2:end]]
end

function _sswgo_triggered(summary)
    for buffer in (0.25, 1.0, 4.0)
        row_at(sigma) = only([r for r in summary
            if r["method"] == "sigma_overlap" &&
               parse(Float64, r["buffer_over_c"]) == buffer &&
               parse(Float64, r["sigma_over_c"]) == sigma])
        r008, r015, r030 = row_at(0.08), row_at(0.15), row_at(0.3)
        cl(r) = abs(parse(Float64, r["delta_CL_percent"]))
        gamma(r) = parse(Float64, r["gamma_error_percent"])
        cl(r015) < cl(r030) && gamma(r015) < gamma(r030) &&
            cl(r008) < cl(r015) && gamma(r008) < gamma(r015) ||
            return true
    end
    return false
end

function run_ssw_growing_overlap_prediction()
    split_path = get(ENV, "SSWGO_SUMMARY",
        joinpath("data", "ssw_sheet_particle_split", "split_summary.csv"))
    output = get(ENV, "SSWRP_OUTPUT",
        joinpath("data", "ssw_representation_probe"))
    summary = _sswgo_read_csv(split_path)
    _sswgo_triggered(summary) || return (; triggered=false, path=nothing)
    path = joinpath(output,
        "phase_b_prediction_sigma_overlap_growing.csv")
    isfile(path) && return (; triggered=true, path)
    rolled_paths = sort(filter(p -> startswith(basename(p), "rolledup_t"),
        readdir(joinpath(output, "states"); join=true)))
    rolled = Serialization.deserialize.(rolled_paths[end-2:end])
    config = last(rolled).config
    control = ssw_frozen_resolve(rolled, config; first_row=1, last_row=1)
    predictions = NamedTuple[]
    for buffer in (0.25, 1.0, 4.0),
            sigma in (0.08, 0.15, 0.3, 0.6, 1.2)
        overlap = 1.3 * sqrt(0.3 / sigma)
        row = clamp(round(Int, buffer / config.dt_star), 1,
            last(rolled).wake.nwakes[])
        prediction = ssw_frozen_resolve(rolled, config; first_row=row,
            last_row=last(rolled).wake.nwakes[],
            method=pnl.SigmaOverlap(sigma * config.c, overlap))
        dcl = 100 * (prediction.CL - control.CL) / abs(control.CL)
        dgamma = 100 * maximum(abs.(prediction.gamma - control.gamma)) /
            max(maximum(abs, control.gamma), eps(Float64))
        push!(predictions, (; sigma_over_c=sigma, overlap,
            buffer_over_c=buffer, delta_CL_percent=dcl,
            gamma_error_percent=dgamma))
    end
    _sswrp_write_overlap_predictions(predictions, path)
    return (; triggered=true, path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = run_ssw_growing_overlap_prediction()
    println(result.triggered ? "growing-overlap prediction: $(result.path)" :
        "fixed-overlap refinement passed; growing path not triggered")
end
