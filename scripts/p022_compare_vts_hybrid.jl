#!/usr/bin/env julia

using TOML
using Statistics: mean

# These fields define the marched physical/numerical case. Formulation-specific
# controls and execution metadata are intentionally excluded. A mismatch is a
# hard error: this tool compares two independent runs, never two formulations
# marched through one mutable simulation state.
const MATCH_KEYS = (
    "mesh_key", "mesh_file", "RPM", "R", "rho", "NT", "dt", "total_steps",
    "requested_revs", "required_revs", "truncation_depth_R", "nwakerows",
    "p_per_step", "overlap", "particle_shedding", "conversion",
    "wake_core_size", "wake_nu", "wake_core_beta", "core_size_panel",
    "core_size_targets", "merge_particles", "particle_relax", "spinup_revs",
    "spinup_start_fraction", "magVinf_peak", "magVinf_end",
    "freestream_ramp_revs", "freestream_hold_revs",
    "freestream_withdraw_revs", "settle_revs", "ground_enable", "ground_h_r",
    "ground_radius_r", "ground_panel_length_r", "ground_particle_policy",
    "ground_damp_band_r", "nrotors", "rotor_spacing_r", "rotor_directions",
    "backend_body_order", "backend_wake_order", "vpm_arraytype",
    "gpu_influence_mode", "gpu_allow_fallback", "shared_pfield",
    "shared_rotor_operator", "body_hessian_to_particles",
    "panel_wake_hessian_to_particles", "panel_wake_velocity_to_particles",
    "particle_hessian_self", "particle_body_overlap_action",
    "particle_body_overlap_core_ratio", "particle_body_overlap_every",
    "das_eta_kinematic", "das_chord_fraction", "das_uniform_dsigma",
    "sigma_chord_fraction", "sigma_floor_r", "das_sigma_lambda",
    "das_curvature_beta", "das_arc_placed", "das_min_displacement",
    "das_kinematic_arc", "das_refresh", "shed_with_induced_velocity",
)

function final_circulation_means(metadata_path)
    monitor_dir = joinpath(dirname(metadata_path), "monitors")
    isdir(monitor_dir) || return Dict{String,Float64}()
    files = sort(filter(f -> occursin("bound_circulation_system", f) &&
                            endswith(f, ".csv"), readdir(monitor_dir)))
    result = Dict{String,Float64}()
    for file in files
        lines = readlines(joinpath(monitor_dir, file))
        length(lines) > 1 || continue
        rows = split.(lines[2:end], ',')
        final_step = maximum(parse(Int, row[1]) for row in rows)
        vals = [parse(Float64, row[6]) for row in rows
                if parse(Int, row[1]) == final_step &&
                   isfinite(parse(Float64, row[6]))]
        isempty(vals) || (result[file] = mean(vals))
    end
    return result
end

function compare_vts_hybrid(vts_path, hybrid_path, output_path)
    vts_path = abspath(vts_path)
    hybrid_path = abspath(hybrid_path)
    output_path = abspath(output_path)
    vts = TOML.parsefile(vts_path)
    hybrid = TOML.parsefile(hybrid_path)

    get(vts, "formulation", "") == "velocity" || error(
        "VTS metadata must have formulation=velocity: $vts_path")
    get(hybrid, "formulation", "") == "hybrid" || error(
        "hybrid metadata must have formulation=hybrid: $hybrid_path")
    for key in MATCH_KEYS
        haskey(vts, key) || error("VTS metadata lacks matched-setting key $key")
        haskey(hybrid, key) || error("hybrid metadata lacks matched-setting key $key")
        isequal(vts[key], hybrid[key]) || error(
            "matched-setting mismatch for $key: $(vts[key]) != $(hybrid[key])")
    end

    result = Dict{String,Any}(
        "artifact" => "p022_vts_vs_hybrid",
        "vts_metadata" => vts_path,
        "hybrid_metadata" => hybrid_path,
        "matched_settings_verified" => true,
        "nrotors" => vts["nrotors"],
        "ground_enable" => vts["ground_enable"],
    )
    for i in 1:Int(vts["nrotors"])
        for metric in ("CT_window_mean_r", "CQ_window_mean_r")
            key = metric * string(i)
            haskey(vts, key) && haskey(hybrid, key) || continue
            result["vts_" * key] = vts[key]
            result["hybrid_" * key] = hybrid[key]
            result["delta_" * key] = hybrid[key] - vts[key]
        end
    end
    result["vts_circulation_final_mean"] = final_circulation_means(vts_path)
    result["hybrid_circulation_final_mean"] = final_circulation_means(hybrid_path)

    mkpath(dirname(output_path))
    open(output_path, "w") do io
        TOML.print(io, result; sorted=true)
    end
    println("Wrote matched independent-run VTS-vs-hybrid artifact: $output_path")
    return result
end

function main(args=ARGS)
    length(args) in (2, 3) || error(
        "usage: julia p022_compare_vts_hybrid.jl VTS_METADATA HYBRID_METADATA [OUTPUT]")
    output_path = length(args) == 3 ? args[3] :
        joinpath(dirname(abspath(args[2])), "p022_vts_vs_hybrid.toml")
    compare_vts_hybrid(args[1], args[2], output_path)
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
