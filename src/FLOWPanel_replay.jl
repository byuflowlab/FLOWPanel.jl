#=##############################################################################
# DESCRIPTION
#   Replay saved VTK simulation output through post-processing monitors without
#   advancing the solver, wake, or kinematics.
=###############################################################################

struct ReplayResult{TS,TW,TF,TR,TM}
    systems::TS
    wakes::TW
    frames::TF
    t_range::TR
    steps::Vector{Int}
    monitors::TM
end

struct ReplayLoadedFields
    fields::Set{Symbol}
    pressure::Union{Nothing, Vector{Float64}}
    force::Union{Nothing, Matrix{Float64}}
end

function _body_vtu_path(path::AbstractString, body_name::AbstractString, idx::Int)
    return joinpath(path, body_name, "$(body_name).$(idx).vtu")
end

function _read_vtu_cells(vtk)
    vcells = ReadVTK.get_cells(vtk)
    all(Int(t) == 5 for t in vcells.types) ||
        throw(ArgumentError("Replay supports triangular VTU cells only; got VTK cell types $(unique(vcells.types))."))
    ncells = length(vcells.offsets)
    cells = Matrix{Int}(undef, 3, ncells)
    start = 1
    for i in 1:ncells
        stop = vcells.offsets[i]
        stop - start + 1 == 3 ||
            throw(ArgumentError("Replay supports triangular VTU cells only; cell $(i) has $(stop - start + 1) nodes."))
        cells[:, i] .= vcells.connectivity[start:stop]
        start = stop + 1
    end
    return cells
end

function _strength_names_from_cell_data(cell_data)
    names = String[]
    "sigma" in keys(cell_data) && push!(names, "sigma")
    "mu" in keys(cell_data) && push!(names, "mu")
    "gamma" in keys(cell_data) && push!(names, "gamma")
    isempty(names) && throw(ArgumentError("Cannot reconstruct replay body: missing strength field; expected one of sigma, mu, gamma."))
    return Tuple(names)
end

_kernel_from_strength_names(names) =
    names == ("sigma",) ? ConstantSource :
    names == ("mu",) ? ConstantDoublet :
    names == ("gamma",) ? VortexRing :
    names == ("sigma", "mu") ? Union{ConstantSource, ConstantDoublet} :
    names == ("sigma", "gamma") ? Union{ConstantSource, VortexRing} :
    throw(ArgumentError("Cannot infer replay body kernel from strength fields $(names)."))

function _body_kind_string(body::AbstractBody)
    body isa RigidWakeBody && return "RigidWakeBody"
    body isa NonLiftingBody && return "NonLiftingBody"
    return string(nameof(typeof(body)))
end

function _wake_shedding_manifest(method)
    if method isa NoShed
        return Dict{String, Any}("type" => "NoShed")
    elseif method isa SigmaPPS
        return Dict{String, Any}(
            "type" => "SigmaPPS",
            "sigma" => method.sigma,
            "p_per_step" => method.p_per_step,
        )
    elseif method isa SigmaOverlap
        return Dict{String, Any}(
            "type" => "SigmaOverlap",
            "sigma" => method.sigma,
            "overlap" => method.overlap,
        )
    elseif method isa OverlapPPS
        return Dict{String, Any}(
            "type" => "OverlapPPS",
            "overlap" => method.overlap,
            "p_per_step" => method.p_per_step,
        )
    elseif method isa StationSigmaOverlap
        return Dict{String, Any}(
            "type" => "StationSigmaOverlap",
            "sigmas" => [collect(sig) for sig in method.sigmas],
            "overlap" => method.overlap,
        )
    else
        return _metadata_unsupported_dict(typeof(method))
    end
end

function _particle_policy_manifest(policy)
    if policy isa MinGamma
        return Dict{String, Any}("type" => "MinGamma", "threshold" => policy.threshold)
    elseif policy isa RelativeMinGamma
        return Dict{String, Any}("type" => "RelativeMinGamma", "fraction" => policy.fraction)
    elseif policy isa MaxGamma
        return Dict{String, Any}("type" => "MaxGamma", "threshold" => policy.threshold)
    elseif policy isa GlobalBox
        return Dict{String, Any}(
            "type" => "GlobalBox",
            "xmin" => collect(policy.xmin),
            "xmax" => collect(policy.xmax),
        )
    elseif policy isa GlobalCylinder
        return Dict{String, Any}(
            "type" => "GlobalCylinder",
            "origin" => collect(policy.origin),
            "extrude" => collect(policy.extrude),
            "radius" => policy.radius,
        )
    elseif policy isa FrameBox
        return Dict{String, Any}(
            "type" => "FrameBox",
            "i_frame" => policy.i_frame,
            "xmin" => collect(policy.xmin),
            "xmax" => collect(policy.xmax),
        )
    elseif policy isa MergeParticles
        return Dict{String, Any}(
            "type" => "MergeParticles",
            "every" => policy.every,
            "r" => policy.r,
            "r_hash" => policy.r_hash,
            "sigma_relative" => policy.sigma_relative,
            "max_sigma_ratio" => policy.max_sigma_ratio,
            "skip_static" => policy.skip_static,
        )
    elseif policy isa SplitParticles
        return Dict{String, Any}(
            "type" => "SplitParticles",
            "every" => policy.every,
            "verbose" => policy.verbose,
            # Split options are opaque FLOWVPM objects; record that manual
            # reconstruction is still required instead of guessing defaults.
            "opts" => _metadata_unsupported_dict(typeof(policy.opts)),
        )
    else
        # Prepared* variants are runtime products of prepare_particle_policy;
        # manifests record the original unprepared maintenance tuple.
        return _metadata_unsupported_dict(typeof(policy))
    end
end

function _particle_maintenance_manifest(maintenance::ParticleMaintenance)
    return Dict{String, Any}(
        "type" => "ParticleMaintenance",
        "trim_policies" => [_particle_policy_manifest(policy) for policy in maintenance.trim_policies],
        "functional_policies" => [_particle_policy_manifest(policy) for policy in maintenance.functional_policies],
    )
end

function _singleton_tag(mod::Module, obj)::Union{String, Nothing}
    for name in names(mod; all=true)
        occursin('#', String(name)) && continue
        isdefined(mod, name) || continue
        getfield(mod, name) === obj && return "$(nameof(mod)).$(name)"
    end
    return nothing
end

function _resolve_singleton(tag::AbstractString, expected::Type)
    prefix = "FLOWVPM."
    startswith(tag, prefix) ||
        throw(ArgumentError("replay: singleton tag \"$(tag)\" must start with \"$(prefix)\"."))
    sym = Symbol(tag[lastindex(prefix) + 1:end])
    isdefined(FLOWVPM, sym) ||
        throw(ArgumentError("replay: FLOWVPM binding for tag \"$(tag)\" does not exist; this may be a FLOWVPM version mismatch. Try migrate_metadata_toml on the metadata file."))
    x = getfield(FLOWVPM, sym)
    x isa expected ||
        throw(ArgumentError("replay: tag \"$(tag)\" resolved to $(typeof(x)), expected $(expected)."))
    return x
end

function _viscous_manifest(viscous)
    if viscous isa FLOWVPM.Inviscid
        return Dict{String, Any}("type" => "FLOWVPM.Inviscid", "nu" => viscous.nu)
    elseif viscous isa FLOWVPM.CoreSpreading
        kernel_tag = _singleton_tag(FLOWVPM, viscous.zeta)
        return Dict{String, Any}(
            "type" => "FLOWVPM.CoreSpreading",
            "nu" => viscous.nu,
            "core_size" => viscous.sgm0,
            "kernel" => isnothing(kernel_tag) ? "FLOWVPM.zeta_fmm" : kernel_tag,
            "beta" => viscous.beta,
        )
    elseif viscous isa FLOWVPM.ParticleStrengthExchange
        return Dict{String, Any}(
            "type" => "FLOWVPM.ParticleStrengthExchange",
            "nu" => viscous.nu,
            "recalculate_vols" => viscous.recalculate_vols,
        )
    elseif viscous === nothing
        return Dict{String, Any}("type" => "nothing")
    else
        return _metadata_unsupported_dict(typeof(viscous))
    end
end

function _sfs_manifest(sfs)
    if sfs === nothing
        return Dict{String, Any}("type" => "nothing")
    end
    tag = _singleton_tag(FLOWVPM, sfs)
    return isnothing(tag) ? _metadata_unsupported_dict(typeof(sfs)) : Dict{String, Any}("type" => tag)
end

function _relaxation_filter_manifest(filter)
    if filter === FLOWVPM.relax_filter_all
        return Dict{String, Any}("type" => "all")
    elseif filter isa RelaxationPlaneFilter
        return Dict{String, Any}(
            "type" => "RelaxationPlaneFilter",
            "point" => collect(filter.point_local),
            "normal" => collect(filter.normal_local),
            "i_frame" => filter.i_frame,
        )
    else
        return _metadata_unsupported_dict(typeof(filter))
    end
end

function _relaxation_manifest(relaxation::FLOWVPM.Relaxation)
    # Identify the scheme by its relaxation method function. Pedrizzetti and
    # CorrectedPedrizzetti share the same Relaxation struct, differing only in
    # this function, so the method identity is the only reliable discriminator.
    if relaxation.relax === FLOWVPM.relax_pedrizzetti
        type = "FLOWVPM.relax_pedrizzetti"
    elseif relaxation.relax === FLOWVPM.relax_correctedpedrizzetti
        type = "FLOWVPM.relax_correctedpedrizzetti"
    elseif relaxation.nsteps_relax < 0
        type = "none"
    else
        type = string(nameof(typeof(relaxation.relax)))
    end
    return Dict{String, Any}(
        "type" => type,
        "nsteps_relax" => relaxation.nsteps_relax,
        "rlxf" => relaxation.rlxf,
        "filter" => _relaxation_filter_manifest(relaxation.filter),
    )
end

# Generic per-value serialization for the FLOWVPM particle-field optargs NamedTuple.
# Scalars pass through; known non-scalar FLOWVPM types get tagged dicts; anything
# else falls back to an unsupported marker (record-only is acceptable here).
# _singleton_tag/_resolve_singleton are the reusable mechanism for future
# FLOWVPM singleton-tag round trips.
_vpm_optarg_manifest(val::Union{Bool, Integer, AbstractFloat, AbstractString, Symbol}) = val
_vpm_optarg_manifest(val::FLOWVPM.Relaxation) = _relaxation_manifest(val)
_vpm_optarg_manifest(val::FLOWVPM.ViscousScheme) = _viscous_manifest(val)
_vpm_optarg_manifest(val::FLOWVPM.SubFilterScale) = _sfs_manifest(val)
function _vpm_optarg_manifest(val)
    return _metadata_unsupported_dict(typeof(val))
end

function _vpm_optargs_manifest(nt::NamedTuple)
    d = Dict{String, Any}()
    for (k, v) in pairs(nt)
        d[string(k)] = _vpm_optarg_manifest(v)
    end
    return d
end

function _body_manifest_dict(body::AbstractBody, i::Int)
    d = Dict{String, Any}(
        "i" => i,
        "kind" => _body_kind_string(body),
        "strength_names" => collect(strength_names(body)),
        "dbc" => has_dirichlet_bc(body),
        "core_size_panel" => body.core_size_panel,
        "core_size_targets" => body.core_size_targets,
        "kernelcutoff" => body.kernelcutoff,
        "watertight" => body.watertight,
    )
    if body isa RigidWakeBody
        d["semiinfinite_wake"] = body.semiinfinite_wake
        d["shedding"] = [_matrix_to_rows(s) for s in body.shedding]
        d["Das"] = [_matrix_to_rows(Das) for Das in body.Das]
    end
    return d
end

_conversion_manifest(::LegacyEdgeJumpConversion) = Dict{String,Any}(
    "version" => 1, "type" => "LegacyEdgeJumpConversion")

function _conversion_manifest(c::SurfaceVorticityConversion)
    return Dict{String,Any}(
        "version" => 1,
        "type" => "SurfaceVorticityConversion",
        "sigma" => c.sigma,
        "overlap" => c.overlap,
        "rank_rtol" => c.rank_rtol,
        "geometry_rtol" => c.geometry_rtol,
        "attribution" => String(c.attribution),
        "diagnose_nearfield" => c.diagnose_nearfield,
    )
end

_conversion_fingerprint(c::AbstractPanelParticleConversion) =
    sprint(io -> TOML.print(io, Dict("conversion" => _conversion_manifest(c)); sorted=true))

function _require_conversion_key(meta, key)
    haskey(meta, key) || throw(ArgumentError(
        "replay: SurfaceVorticityConversion metadata is missing required key \"$(key)\""))
    return meta[key]
end

"""Strict decoder for the version-1 panel/particle conversion schema."""
function _deserialize_conversion(meta)
    meta isa AbstractDict || throw(ArgumentError(
        "replay: conversion metadata must be a table"))
    haskey(meta, "version") || throw(ArgumentError(
        "replay: conversion metadata is missing required key \"version\""))
    version = try Int(meta["version"]) catch; throw(ArgumentError(
        "replay: conversion version must be integer 1")) end
    version == 1 || throw(ArgumentError(
        "replay: unsupported conversion metadata version $(version)"))
    haskey(meta, "type") || throw(ArgumentError(
        "replay: conversion metadata is missing required key \"type\""))
    ctype = String(meta["type"])
    if ctype == "LegacyEdgeJumpConversion"
        return LegacyEdgeJumpConversion()
    elseif ctype == "SurfaceVorticityConversion"
        sigma = Float64(_require_conversion_key(meta, "sigma"))
        overlap = Float64(_require_conversion_key(meta, "overlap"))
        rank_rtol = Float64(_require_conversion_key(meta, "rank_rtol"))
        geometry_rtol = Float64(_require_conversion_key(meta, "geometry_rtol"))
        attribution = Symbol(String(_require_conversion_key(meta, "attribution")))
        diagnose_nearfield = Bool(_require_conversion_key(meta, "diagnose_nearfield"))
        return SurfaceVorticityConversion(sigma; overlap, rank_rtol,
            geometry_rtol, attribution, diagnose_nearfield)
    end
    throw(ArgumentError("replay: unsupported conversion metadata type \"$(ctype)\""))
end

function _wake_manifest_dict(wake, i::Int)
    d = Dict{String, Any}("i" => i)
    if isnothing(wake)
        d["type"] = "nothing"
    elseif wake isa PanelParticleWake
        d["type"] = "PanelParticleWake"
        # logical value: 0 for a convert-at-shed wake (BRAINSTORM 024), whose
        # N=1 storage is an implementation detail the constructor recreates
        d["nwakerows"] = _logical_nwakerows(wake.panel_wake)
        d["max_particles"] = size(wake.pfield.particles, 2)
        d["core_size"] = wake.panel_wake.core_size
        d["shed_with_induced_velocity"] = wake.panel_wake.shed_with_induced_velocity
        d["unsteady_filament"] = wake.panel_wake.unsteady_filament
        d["freestream_convection"] = wake.panel_wake.freestream_convection
        d["particle_core_size"] = wake.particle_core_size
        d["conversion"] = _conversion_manifest(wake.conversion)
        if wake.conversion isa LegacyEdgeJumpConversion
            # Preserve the historical line-policy representation exactly. The
            # smooth strategy has no such policies and must not fabricate them.
            d["method_trailing"] = _wake_shedding_manifest(wake.method_trailing)
            d["method_unsteady"] = _wake_shedding_manifest(wake.method_unsteady)
        end
        d["particle_maintenance"] = _particle_maintenance_manifest(wake.particle_maintenance)
        # Resolved FLOWVPM particle-field construction options (viscous, SFS,
        # relaxation scheme, formulation, kernel). Serialized generically so new
        # scalar optargs are captured without bespoke code.
        d["pfield_optargs"] = _vpm_optargs_manifest(wake.pfield_optargs)
    elseif wake isa PanelWake
        d["type"] = "PanelWake"
        d["nwakerows"] = _logical_nwakerows(wake)
        d["core_size"] = wake.core_size
        d["shed_with_induced_velocity"] = wake.shed_with_induced_velocity
        d["unsteady_filament"] = wake.unsteady_filament
        d["include_final_filament"] = wake.include_final_filament
        d["freestream_convection"] = wake.freestream_convection
    else
        d["type"] = string(nameof(typeof(wake)))
    end
    return d
end

function _construct_body_from_metadata(nodes, cells, body_meta, cell_data)
    names = haskey(body_meta, "strength_names") ? Tuple(String.(body_meta["strength_names"])) :
            _strength_names_from_cell_data(cell_data)
    E = _kernel_from_strength_names(names)
    dbc = Bool(get(body_meta, "dbc", false))
    kwargs = (;
        DBC=dbc,
        # NOTE: replay metadata written before 2026-08-22 spelled these keys
        # "kerneloffset_panel"/"kerneloffset_targets"/"kerneloffset". The old
        # keys stay accepted (new keys win) so existing restart datasets load.
        core_size_panel=Float64(get(body_meta, "core_size_panel", get(body_meta, "kerneloffset_panel", get(body_meta, "core_size", get(body_meta, "kerneloffset", 1e-8))))),
        core_size_targets=Float64(get(body_meta, "core_size_targets", get(body_meta, "kerneloffset_targets", get(body_meta, "core_size", get(body_meta, "kerneloffset", 1e-8))))),
        kernelcutoff=Float64(get(body_meta, "kernelcutoff", 1e-14)),
        watertight=Bool(get(body_meta, "watertight", false)),
        ensure_winding=false,
    )
    kind = String(get(body_meta, "kind", "NonLiftingBody"))
    if kind == "RigidWakeBody"
        haskey(body_meta, "shedding") ||
            throw(ArgumentError("Cannot reconstruct RigidWakeBody for replay: missing shedding metadata. Provide a replay manifest or reconstruct callback."))
        shedding = [_rows_to_matrix(rows, Int, 6) for rows in body_meta["shedding"]]
        body = RigidWakeBody{E}(copy(nodes), copy(cells), shedding;
            check_mesh=false,
            semiinfinite_wake=Bool(get(body_meta, "semiinfinite_wake", true)),
            kwargs...)
        if haskey(body_meta, "Das")
            Das = [_rows_to_matrix(rows, eltype(body.nodes), 3) for rows in body_meta["Das"]]
            length(Das) == length(body.Das) && foreach((a, b) -> (a .= b), body.Das, Das)
        end
        return body
    elseif kind == "NonLiftingBody"
        return NonLiftingBody{E}(copy(nodes), copy(cells); kwargs...)
    else
        throw(ArgumentError("Cannot reconstruct replay body kind $(kind). Provide reconstruct callback."))
    end
end

function _load_optional_body_fields!(body::AbstractBody, cell_data)
    loaded = Set{Symbol}()
    pressure = nothing
    force = nothing
    for (i, sname) in enumerate(strength_names(body))
        sname in keys(cell_data) || throw(ArgumentError("Cannot load replay body: missing strength field $(sname)."))
        body.strength[:, i] .= ReadVTK.get_data(cell_data[sname])
        push!(loaded, Symbol(sname))
    end
    if "potential" in keys(cell_data)
        body.potential .= ReadVTK.get_data(cell_data["potential"])
        push!(loaded, :potential)
    end
    if "velocity" in keys(cell_data)
        body.velocity .= ReadVTK.get_data(cell_data["velocity"])
        push!(loaded, :velocity)
    end
    if "gauge pressure" in keys(cell_data)
        pressure = Vector{Float64}(ReadVTK.get_data(cell_data["gauge pressure"]))
        push!(loaded, :P)
    end
    if "F" in keys(cell_data)
        force = Matrix{Float64}(ReadVTK.get_data(cell_data["F"]))
        push!(loaded, :F)
    end
    if "normals" in keys(cell_data)
        body.normals .= ReadVTK.get_data(cell_data["normals"])
        push!(loaded, :normals)
    end
    return ReplayLoadedFields(loaded, pressure, force)
end

function _load_body_step!(body::AbstractBody, path, body_name, idx)
    vtk = ReadVTK.VTKFile(_body_vtu_path(path, body_name, idx))
    points = ReadVTK.get_points(vtk)
    size(points) == size(body.nodes) ||
        throw(ArgumentError("Replay body $(body_name) step $(idx) has point size $(size(points)); expected $(size(body.nodes))."))
    body.nodes .= points
    cell_data = ReadVTK.get_cell_data(vtk)
    loaded = _load_optional_body_fields!(body, cell_data)
    calc_normals!(body)
    calc_controlpoints!(body)
    return loaded
end

function _load_replay_wake_step!(w::PanelWake, path, wake_name, idx)
    vts_path = joinpath(path, wake_name, "$(wake_name).1.$(idx).vts")
    if !isfile(vts_path)
        w.nwakes[] = 0
        w.overflowed[] = false
        return w
    end
    return _load_panel_wake_vtk!(w, path, wake_name, idx)
end

function _load_replay_wake_step!(w::PanelParticleWake, path, wake_name, idx)
    vts_path = joinpath(path, wake_name, "$(wake_name).1.$(idx).vts")
    if !isfile(vts_path)
        w.panel_wake.nwakes[] = 0
        w.panel_wake.overflowed[] = false
        return w
    end
    vtp_path = joinpath(path, wake_name * "_particles", "$(wake_name)_particles.$(idx).vtp")
    isfile(vtp_path) || return _load_replay_wake_step!(w.panel_wake, path, wake_name, idx)
    return _load_panel_particle_wake_vtk!(w, path, wake_name, idx)
end

function _metadata_step_record(metadata, idx::Int)
    metadata === nothing && return nothing
    steps = get(metadata, "step", Any[])
    k = findfirst(s -> Int(get(s, "i_step", typemin(Int))) == idx, steps)
    return isnothing(k) ? nothing : steps[k]
end

function _continuation_record(step, i_wake::Int)
    step === nothing && return nothing
    records = get(step, "wake_continuation", Any[])
    matches = filter(r -> Int(get(r, "i", -1)) == i_wake, records)
    step_label = get(step, "i_step", "?")
    length(matches) <= 1 || throw(WakeContinuationStateError(
        "step $(step_label) has multiple continuation records for wake $(i_wake)"))
    return isempty(matches) ? nothing : only(matches)
end

function _required_continuation(record, key, i_wake, idx)
    haskey(record, key) || throw(WakeContinuationStateError(
        "wake $(i_wake) step $(idx) continuation state is missing \"$(key)\""))
    return record[key]
end

"""
Restore non-VTK panel/particle continuation state after a wake VTK load. Old
metadata that selects the legacy conversion retains the historical inferred
defaults; smooth conversion never guesses an absent or ambiguous handoff.
"""
function _restore_wake_continuation!(wakes::Tuple, metadata, idx::Int)
    step = _metadata_step_record(metadata, idx)
    for (i, wake) in enumerate(wakes)
        wake isa PanelParticleWake || continue
        record = _continuation_record(step, i)
        if record === nothing
            wake.conversion isa SurfaceVorticityConversion &&
                throw(WakeContinuationStateError(
                    "smooth wake $(i) step $(idx) has no per-step continuation state"))
            continue
        end

        phase = String(_required_continuation(record, "snapshot_phase", i, idx))
        phase == "pre_end_of_step_shedding" || throw(WakeContinuationStateError(
            "wake $(i) step $(idx) has unsupported snapshot phase \"$(phase)\""))
        step_identity = Int(_required_continuation(record, "step_identity", i, idx))
        step_identity == idx || throw(WakeContinuationStateError(
            "wake $(i) continuation step identity $(step_identity) does not match loaded step $(idx)"))
        fingerprint = String(_required_continuation(record,
            "conversion_fingerprint", i, idx))
        expected_fingerprint = _conversion_fingerprint(wake.conversion)
        fingerprint == expected_fingerprint || throw(WakeContinuationStateError(
            "wake $(i) step $(idx) conversion fingerprint does not match the reconstructed conversion"))

        active = Int(_required_continuation(record, "active_row_count", i, idx))
        active == wake.panel_wake.nwakes[] || throw(WakeContinuationStateError(
            "wake $(i) step $(idx) VTK has $(wake.panel_wake.nwakes[]) active rows but metadata records $(active)"))
        capacity = size(wake.panel_wake.nodes[1], 2) - 1
        0 <= active <= capacity || throw(WakeContinuationStateError(
            "wake $(i) step $(idx) active row count $(active) is outside 0:$(capacity)"))
        live_rows = Int(_required_continuation(record, "live_rows", i, idx))
        0 <= live_rows <= active || throw(WakeContinuationStateError(
            "wake $(i) step $(idx) live row count $(live_rows) is incompatible with $(active) active rows"))
        count = Int(_required_continuation(record, "conversion_count", i, idx))
        count >= 0 || throw(WakeContinuationStateError(
            "wake $(i) step $(idx) has negative conversion count $(count)"))
        handoff = Bool(_required_continuation(record, "handoff_active", i, idx))
        handoff == (count > 0) || throw(WakeContinuationStateError(
            "wake $(i) step $(idx) handoff flag is inconsistent with conversion count $(count)"))
        weight = Float64(_required_continuation(record, "handoff_weight", i, idx))
        isfinite(weight) && 0 <= weight <= 1 || throw(WakeContinuationStateError(
            "wake $(i) step $(idx) has invalid handoff weight $(weight)"))

        terminal = get(record, "terminal_strength", nothing)
        if terminal === nothing
            wake.conversion isa SurfaceVorticityConversion &&
                throw(WakeContinuationStateError(
                    "smooth wake $(i) step $(idx) continuation state is missing \"terminal_strength\""))
        else
            length(terminal) == length(wake.panel_wake.strength) ||
                throw(WakeContinuationStateError(
                    "wake $(i) step $(idx) terminal-strength surface count does not match the wake"))
            for (isurf, rows) in enumerate(terminal)
                s = wake.panel_wake.strength[isurf]
                restored = try
                    _rows_to_matrix(rows, eltype(s), size(s, 1))
                catch err
                    throw(WakeContinuationStateError(
                        "wake $(i) step $(idx) terminal strength $(isurf) is malformed: $(sprint(showerror, err))"))
                end
                size(restored) == (size(s,1), size(s,3)) ||
                    throw(WakeContinuationStateError(
                        "wake $(i) step $(idx) terminal strength $(isurf) has size $(size(restored)); expected $((size(s,1), size(s,3)))"))
                s[:, active + 1, :] .= restored
            end
        end

        pw = wake.panel_wake
        pw.nwakes[] = active
        pw.overflowed[] = Bool(_required_continuation(record, "overflowed", i, idx))
        pw.live_rows[] = live_rows
        pw.live_step_id[] = Int(_required_continuation(record, "live_step_id", i, idx))
        pw.particle_handoff_active[] = handoff
        pw.particle_handoff_weight[] = weight
        wake.conversion_count[] = count
    end
    return wakes
end

function _read_body_metadata(path, run_name, idx, manifest)
    body_pvds = sort(filter(f -> occursin(Regex("^" * run_name * "_body\\d+\\.pvd\$"), f), readdir(path)))
    isempty(body_pvds) && throw(ArgumentError("Replay found no body PVD files matching $(run_name)_body*.pvd in $(path)."))
    bodies = Any[]
    metadata = Dict{String, Any}(
        "path" => path,
        "run_name" => run_name,
        "manifest" => manifest,
        "missing" => String[],
    )
    body_metas = manifest === nothing ? Any[] : get(manifest, "body", Any[])
    for (i, _) in enumerate(body_pvds)
        body_name = run_name * "_body$(i)"
        vtu = _body_vtu_path(path, body_name, idx)
        isfile(vtu) || throw(ArgumentError("Replay body VTU not found: $(vtu)."))
        vtk = ReadVTK.VTKFile(vtu)
        nodes = Matrix{Float64}(ReadVTK.get_points(vtk))
        cells = _read_vtu_cells(vtk)
        cell_data = ReadVTK.get_cell_data(vtk)
        body_meta = i <= length(body_metas) ? body_metas[i] : Dict{String, Any}()
        if String(get(body_meta, "kind", "NonLiftingBody")) == "RigidWakeBody" && !haskey(body_meta, "shedding")
            push!(metadata["missing"], "body$(i): shedding metadata")
            body_meta = copy(body_meta)
            body_meta["kind"] = "NonLiftingBody"
        end
        push!(bodies, _construct_body_from_metadata(nodes, cells, body_meta, cell_data))
    end
    metadata["bodies"] = bodies
    return Tuple(bodies), metadata
end

function _construct_wakes_from_manifest(systems::Tuple, manifest)
    manifest === nothing && return Tuple(nothing for _ in systems)
    wake_metas = get(manifest, "wake", Any[])
    wakes = Any[]
    for i in eachindex(systems)
        wmeta = i <= length(wake_metas) ? wake_metas[i] : Dict{String, Any}("type" => "nothing")
        wtype = String(get(wmeta, "type", "nothing"))
        if wtype == "nothing"
            push!(wakes, nothing)
        elseif wtype == "PanelWake"
            systems[i] isa AbstractLiftingBody ||
                throw(ArgumentError("Cannot reconstruct PanelWake for non-lifting body $(i)."))
            push!(wakes, PanelWake(systems[i];
                nwakerows=Int(get(wmeta, "nwakerows", 100)),
                core_size=Float64(get(wmeta, "core_size", 1e-3)),
                shed_with_induced_velocity=Bool(get(wmeta, "shed_with_induced_velocity", true)),
                unsteady_filament=Bool(get(wmeta, "unsteady_filament", true)),
                # Preserve the historical default for manifests written before
                # finite panel-only wakes recorded this setting explicitly.
                include_final_filament=Bool(get(wmeta, "include_final_filament", true)),
                freestream_convection=Bool(get(wmeta, "freestream_convection", false))))
        elseif wtype == "PanelParticleWake"
            systems[i] isa AbstractLiftingBody ||
                throw(ArgumentError("Cannot reconstruct PanelParticleWake for non-lifting body $(i)."))
            # A wholly absent conversion block is the historical legacy
            # default. Once present, the version/type/parameters are strict and
            # never use the warning/fallback shedding decoder.
            conversion = haskey(wmeta, "conversion") ?
                _deserialize_conversion(wmeta["conversion"]) :
                LegacyEdgeJumpConversion()
            particle_maintenance = _deserialize_particle_maintenance(get(wmeta, "particle_maintenance", Dict{String, Any}("type" => "ParticleMaintenance")))
            # FLOWVPM particle-field optargs now live under "pfield_optargs"; fall
            # back to legacy top-level "viscous"/"SFS" keys for older metadata files.
            pf_optargs = get(wmeta, "pfield_optargs", wmeta)
            viscous = _deserialize_viscous(get(pf_optargs, "viscous", Dict{String, Any}("type" => "FLOWVPM.Inviscid")))
            sfs = _deserialize_sfs(get(pf_optargs, "SFS", Dict{String, Any}("type" => "FLOWVPM.SFS_default")))
            relaxation = _deserialize_relaxation(get(pf_optargs, "relaxation",
                Dict{String, Any}("type" => "FLOWVPM.relax_correctedpedrizzetti", "nsteps_relax" => 1, "rlxf" => 0.3)))
            common = (;
                nwakerows=Int(get(wmeta, "nwakerows", 3)),
                max_particles=Int(get(wmeta, "max_particles", 10000)),
                core_size=Float64(get(wmeta, "core_size", 1e-3)),
                shed_with_induced_velocity=Bool(get(wmeta, "shed_with_induced_velocity", true)),
                unsteady_filament=Bool(get(wmeta, "unsteady_filament", true)),
                freestream_convection=Bool(get(wmeta, "freestream_convection", false)),
                particle_maintenance=particle_maintenance,
                # old key "particle_kerneloffset" accepted for legacy restarts
                particle_core_size=Float64(get(wmeta, "particle_core_size", get(wmeta, "particle_kerneloffset", NaN))),
                viscous=viscous,
                SFS=sfs,
                relaxation=relaxation,
                conversion=conversion)
            if conversion isa LegacyEdgeJumpConversion
                method_trailing = _deserialize_wake_shedding(get(wmeta,
                    "method_trailing", Dict{String, Any}("type" => "OverlapPPS")))
                method_unsteady = _deserialize_wake_shedding(get(wmeta,
                    "method_unsteady", Dict{String, Any}("type" => "OverlapPPS")))
                push!(wakes, PanelParticleWake(systems[i]; common...,
                    method_trailing, method_unsteady))
            else
                push!(wakes, PanelParticleWake(systems[i]; common...))
            end
        else
            throw(ArgumentError("Cannot reconstruct replay wake type $(wtype). Provide reconstruct callback."))
        end
    end
    return Tuple(wakes)
end

function _deserialize_wake_shedding(meta)
    wtype = String(get(meta, "type", "unsupported"))
    if wtype == "NoShed"
        return NoShed()
    elseif wtype == "SigmaPPS"
        return SigmaPPS(Float64(get(meta, "sigma", 1.0)), Int(get(meta, "p_per_step", 1)))
    elseif wtype == "SigmaOverlap"
        return SigmaOverlap(Float64(get(meta, "sigma", 1.0)), Float64(get(meta, "overlap", 1.0)))
    elseif wtype == "OverlapPPS"
        return OverlapPPS(Float64(get(meta, "overlap", 1.0)), Int(get(meta, "p_per_step", 1)))
    elseif wtype == "StationSigmaOverlap"
        sigmas = [Float64.(sig) for sig in get(meta, "sigmas", Any[])]
        return StationSigmaOverlap(sigmas, Float64(get(meta, "overlap", 1.0)))
    elseif wtype == "nothing"
        return NoShed()
    else
        @warn "replay: unrecognized wake shedding tag \"$(wtype)\"; falling back to NoShed (result may differ)"
        return NoShed()
    end
end

function _deserialize_particle_policy(meta)
    ptype = String(get(meta, "type", "unsupported"))
    if ptype == "MinGamma"
        return MinGamma(Float64(get(meta, "threshold", 0.0)))
    elseif ptype == "RelativeMinGamma"
        return RelativeMinGamma(Float64(get(meta, "fraction", 0.0)))
    elseif ptype == "MaxGamma"
        return MaxGamma(Float64(get(meta, "threshold", 0.0)))
    elseif ptype == "GlobalBox"
        return GlobalBox(get(meta, "xmin", [0.0, 0.0, 0.0]), get(meta, "xmax", [0.0, 0.0, 0.0]))
    elseif ptype == "GlobalCylinder"
        return GlobalCylinder(
            get(meta, "origin", [0.0, 0.0, 0.0]),
            get(meta, "extrude", [1.0, 0.0, 0.0]),
            Float64(get(meta, "radius", 0.0)),
        )
    elseif ptype == "FrameBox"
        return FrameBox(Int(get(meta, "i_frame", 1)), get(meta, "xmin", [0.0, 0.0, 0.0]), get(meta, "xmax", [0.0, 0.0, 0.0]))
    elseif ptype == "MergeParticles"
        return MergeParticles(
            every=Int(get(meta, "every", 0)),
            r=Float64(get(meta, "r", 0.5)),
            r_hash=Float64(get(meta, "r_hash", -1.0)),
            sigma_relative=Bool(get(meta, "sigma_relative", true)),
            max_sigma_ratio=Float64(get(meta, "max_sigma_ratio", 2.0)),
            skip_static=Bool(get(meta, "skip_static", true)),
        )
    else
        return nothing
    end
end

function _deserialize_particle_maintenance(meta)
    if String(get(meta, "type", "unsupported")) != "ParticleMaintenance"
        return ParticleMaintenance()
    end
    trim = [_deserialize_particle_policy(p) for p in get(meta, "trim_policies", Any[]) if !isnothing(_deserialize_particle_policy(p))]
    functional = [_deserialize_particle_policy(p) for p in get(meta, "functional_policies", Any[]) if !isnothing(_deserialize_particle_policy(p))]
    return ParticleMaintenance((trim..., functional...))
end

function _deserialize_viscous(meta)
    vtype = String(get(meta, "type", "unsupported"))
    if vtype == "FLOWVPM.Inviscid"
        return FLOWVPM.Inviscid()
    elseif vtype == "FLOWVPM.CoreSpreading"
        kernel_fn = _resolve_singleton(String(get(meta, "kernel", "FLOWVPM.zeta_fmm")), Function)
        return FLOWVPM.CoreSpreading(
            Float64(get(meta, "nu", 0.0)),
            Float64(get(meta, "core_size", 0.0)),
            kernel_fn;
            beta=Float64(get(meta, "beta", 1.5)),
        )
    elseif vtype == "FLOWVPM.ParticleStrengthExchange"
        return FLOWVPM.ParticleStrengthExchange{Float64}(
            Float64(get(meta, "nu", 0.0));
            recalculate_vols=Bool(get(meta, "recalculate_vols", true)),
        )
    elseif vtype == "nothing"
        return FLOWVPM.Inviscid()
    else
        @warn "replay: unrecognized viscous tag \"$(vtype)\"; falling back to Inviscid (result may differ)"
        return FLOWVPM.Inviscid()
    end
end

function _deserialize_sfs(meta)
    stype = String(get(meta, "type", "unsupported"))
    if stype == "nothing" || stype == "unsupported"
        @warn "replay: unrecognized SFS tag \"$(stype)\"; falling back to SFS_none (result may differ)"
        return FLOWVPM.SFS_none
    end
    return _resolve_singleton(stype, FLOWVPM.SubFilterScale)
end

const _LEGACY_TAG_REMAP = Dict{String, String}()

function _canonicalize_flowvpm_tag(tag::AbstractString, unresolved::Vector{String})
    startswith(tag, "FLOWVPM.") || return tag
    resolved_tag = tag
    try
        obj = _resolve_singleton(resolved_tag, Any)
        canonical = _singleton_tag(FLOWVPM, obj)
        return isnothing(canonical) ? resolved_tag : canonical
    catch
        if haskey(_LEGACY_TAG_REMAP, tag)
            resolved_tag = _LEGACY_TAG_REMAP[tag]
            try
                obj = _resolve_singleton(resolved_tag, Any)
                canonical = _singleton_tag(FLOWVPM, obj)
                return isnothing(canonical) ? resolved_tag : canonical
            catch
            end
        end
        push!(unresolved, tag)
        return tag
    end
end

function _canonicalize_metadata_tags!(x, unresolved::Vector{String})
    if x isa AbstractDict
        for (k, v) in collect(x)
            String(k) == "julia_type" && continue
            if v isa AbstractString
                x[k] = _canonicalize_flowvpm_tag(v, unresolved)
            else
                _canonicalize_metadata_tags!(v, unresolved)
            end
        end
    elseif x isa AbstractVector
        for i in eachindex(x)
            v = x[i]
            if v isa AbstractString
                x[i] = _canonicalize_flowvpm_tag(v, unresolved)
            else
                _canonicalize_metadata_tags!(v, unresolved)
            end
        end
    end
    return x
end

"""
    migrate_metadata_toml(path::AbstractString; backup=true) -> path

Rewrite a replay metadata TOML file to the current layout and canonical FLOWVPM
binding tags. A backup is written to `path * ".bak"` by default.
"""
function migrate_metadata_toml(path::AbstractString; backup=true)
    data = TOML.parsefile(path)
    unresolved = String[]
    for wmeta in get(data, "wake", Any[])
        wmeta isa AbstractDict || continue
        pf_optargs = get!(wmeta, "pfield_optargs", Dict{String, Any}())
        for key in ("viscous", "SFS", "relaxation")
            if haskey(wmeta, key)
                haskey(pf_optargs, key) || (pf_optargs[key] = wmeta[key])
                delete!(wmeta, key)
            end
        end
        _canonicalize_metadata_tags!(wmeta, unresolved)
    end
    if !isempty(unresolved)
        unique_tags = unique(unresolved)
        throw(ArgumentError("migrate_metadata_toml could not resolve FLOWVPM tags: $(join(unique_tags, ", ")). Add renamed bindings to _LEGACY_TAG_REMAP if needed."))
    end
    backup && cp(path, path * ".bak"; force=true)
    open(path, "w") do io
        TOML.print(io, data)
    end
    return path
end

function _deserialize_relaxation_filter(meta)
    ftype = String(get(meta, "type", "all"))
    if ftype == "all" || ftype == "FLOWVPM.relax_filter_all"
        return FLOWVPM.relax_filter_all
    elseif ftype == "RelaxationPlaneFilter"
        return RelaxationPlaneFilter(
            get(meta, "point", [0.0, 0.0, 0.0]),
            get(meta, "normal", [1.0, 0.0, 0.0]);
            i_frame=Int(get(meta, "i_frame", 0)))
    else
        return FLOWVPM.relax_filter_all
    end
end

function _deserialize_relaxation(meta)
    rtype = String(get(meta, "type", "unsupported"))
    nsteps = Int(get(meta, "nsteps_relax", 1))
    rlxf = Float64(get(meta, "rlxf", 0.3))
    filter = _deserialize_relaxation_filter(get(meta, "filter", Dict{String, Any}("type" => "all")))
    if rtype == "FLOWVPM.relax_pedrizzetti"
        return FLOWVPM.Relaxation(FLOWVPM.relax_pedrizzetti, nsteps, rlxf, filter)
    elseif rtype == "FLOWVPM.relax_correctedpedrizzetti"
        return FLOWVPM.Relaxation(FLOWVPM.relax_correctedpedrizzetti, nsteps, rlxf, filter)
    elseif rtype == "none"
        return FLOWVPM.relaxation_none
    else
        # Unknown/legacy metadata: fall back to the PanelParticleWake default.
        return FLOWVPM.relaxation_correctedpedrizzetti
    end
end

function _frame_from_state_and_static(state, static, default_name, dependent_index)
    TF = Float64
    name = static === nothing ? String(get(state, "name", default_name)) : String(get(static, "name", default_name))
    parent_index = static === nothing ? Int(get(state, "parent_index", -1)) : Int(get(static, "parent_index", -1))
    child_index = static === nothing ? Int.(get(state, "child_index", Int[])) : Int.(get(static, "child_index", Int[]))
    dep_index = static === nothing ? Int.(get(state, "dependent_index", dependent_index)) : Int.(get(static, "dependent_index", dependent_index))
    return ReferenceFrame{TF}(
        FastMultipole.SVector{3,TF}(state["x"]...),
        FastMultipole.SVector{3,TF}(state["v"]...),
        FastMultipole.SVector{3,TF}(state["omega_axis"]...),
        TF(state["omega"]),
        FastMultipole.SMatrix{3,3,TF,9}(state["R"]...),
        FastMultipole.SMatrix{3,3,TF,9}(state["Rp2g"]...),
        name, parent_index, child_index, dep_index)
end

function _frames_from_metadata(path, name, idx, systems, metadata)
    step_frames = _metadata_step_frames(metadata, idx)
    if step_frames === nothing
        frame_file = _frames_toml_path(path, name)
        if isfile(frame_file)
            data = TOML.parsefile(frame_file)
            steps = get(data, "step", Any[])
            i_toml = findfirst(s -> Int(s["i_step"]) == idx, steps)
            isnothing(i_toml) && throw(ArgumentError("Replay frame state for step $(idx) not found in $(frame_file)."))
            step_frames = steps[i_toml]["frame"]
            static_frames = get(metadata === nothing ? data : metadata, "frame", Any[])
        else
            if metadata !== nothing && haskey(metadata, "frame")
                steps_available = [Int(s["i_step"]) for s in get(metadata, "step", Any[])]
                throw(ArgumentError("Replay frame state for step $(idx) not found in $(_metadata_toml_path(path, name)). " *
                    "Frame states available for steps: $(steps_available). " *
                    "This may be caused by a bug in simulate! that overwrote earlier step frame states; re-running the simulation will fix this for future runs. " *
                    "To replay existing data, restrict STEPS to the available steps listed above."))
            end
            return ReferenceFrame(systems[1]; dependent_index=collect(eachindex(systems)))
        end
    else
        static_frames = get(metadata, "frame", Any[])
    end

    frames = ReferenceFrame(systems[1]; dependent_index=collect(eachindex(systems)))
    empty!(frames)
    for (i, fd) in enumerate(step_frames)
        static = i <= length(static_frames) ? static_frames[i] : nothing
        push!(frames, _frame_from_state_and_static(fd, static, "frame$(i)", collect(eachindex(systems))))
    end
    return frames
end

function _apply_reconstruct_callback(reconstruct, path, run_name, metadata, systems, wakes, frames)
    isnothing(reconstruct) && return systems, wakes, frames
    out = reconstruct(path, run_name, metadata)
    out === nothing && return systems, wakes, frames
    if out isa NamedTuple
        systems = get(out, :systems, systems)
        wakes = get(out, :wakes, wakes)
        frames = get(out, :frames, frames)
    elseif out isa Tuple && length(out) == 3
        systems, wakes, frames = out
    else
        systems = out
    end
    systems = _systems_tuple(systems)
    wakes = _wakes_tuple(systems, wakes)
    return systems, wakes, frames
end

function _selected_replay_steps(path, run_name, steps)
    timesteps, idxs = _read_pvd_steps(joinpath(path, run_name * "_body1.pvd"))
    selected = steps === :all ? collect(idxs) :
               steps isa Integer ? [Int(steps)] :
               collect(Int, steps)
    missing = setdiff(selected, idxs)
    isempty(missing) || throw(ArgumentError("Replay steps $(missing) not found in $(run_name)_body1.pvd; available steps are $(idxs)."))
    return timesteps, idxs, selected
end

function _recompute_set(recompute)
    supported = Set([:velocity, :potential, :velocity_gradient, :induced_vorticity])
    recompute === :all && return copy(supported)
    recompute isa Tuple || throw(ArgumentError("recompute must be a tuple of symbols, (:auto,), (:all,), (), or :all."))
    :all in recompute && return copy(supported)
    :auto in recompute && length(recompute) == 1 && return Set{Symbol}([:auto])
    fields = Set(Symbol.(recompute))
    unsupported = setdiff(fields, supported)
    isempty(unsupported) || throw(ArgumentError(
        "Unsupported replay recompute field(s): $(join(sort!(String.(collect(unsupported))), ", ")). " *
        "Supported fields are :velocity, :potential, :velocity_gradient, and :induced_vorticity."))
    return fields
end

function _monitor_replay_requirements(monitors)
    req = Set{Symbol}()
    for m in monitors
        union!(req, Symbol.(monitor_requires(m)))
        if m isa PressureBernoulli
            push!(req, :velocity)
        elseif m isa PressureLaplace
            push!(req, :velocity)
            monitor_requires_body_hessian(m) && push!(req, :velocity_gradient)
            monitor_requires_induced_vorticity(m) && push!(req, :induced_vorticity)
        elseif m isa SurfaceVorticityForce
            push!(req, :velocity)
        end
    end
    return req
end

function _recompute_replay_fields!(systems::Tuple, wakes::Tuple, frames, uinf, fields::AbstractSet, backend_wake, backend_system;
        grad_mu_options=(;))
    if (:velocity in fields) || (:potential in fields) || (:velocity_gradient in fields) || (:induced_vorticity in fields)
        recompute_velocity = :velocity in fields
        recompute_potential = :potential in fields
        recompute_velocity_gradient = :velocity_gradient in fields
        recompute_induced_vorticity = :induced_vorticity in fields
        for body in systems
            recompute_velocity && (body.velocity .= 0)
            recompute_potential && (body.potential .= 0)
            recompute_velocity_gradient && (body.velocity_gradient .= 0)
            recompute_induced_vorticity && (body.induced_vorticity .= 0)
            if recompute_velocity
                body.velocity_kinematic .= 0
                body.angular_velocity .= 0
            end
        end
        wake_sources = _collect_wake_sources(wakes)
        if length(wake_sources) > 0
            influence!(systems, wake_sources, backend_wake; precalc=true,
                scalar_potential=false,
                velocity=recompute_velocity,
                velocity_gradient=Tuple(recompute_velocity_gradient for _ in systems),
                extra_outputs=_induced_vorticity_extra_outputs(systems, recompute_induced_vorticity))
        end
        _set_core_sizes!(systems, :core_size_targets)
        scalar_sources = _filter_scalar_potential_sources(systems)
        if recompute_potential && length(scalar_sources) > 0
            influence!(systems, scalar_sources, backend_system; scalar_potential=true, velocity=false)
        end
        # Mirror simulate!: add κ = n × ∇sμ on top of the wake-induced curl
        # before the body-on-body pass, so replay produces the same
        # induced_vorticity = wake_induced + κ that simulate! does. simulate!
        # normalizes its grad_mu_options with default_basis=:quad
        # (_steady_aerodynamics!), so the same default is applied here — the
        # bare _add_bound_surface_vorticity! default is :tri, which was
        # observed to shift the replayed lamb-vector CT by ~2.5e-3 (4%) on the
        # rotor-hover case relative to the pressures the original simulation
        # wrote.
        recompute_induced_vorticity && _add_bound_surface_vorticity!(systems;
            grad_mu_options=_normalize_grad_mu_options(grad_mu_options;
                default_basis=:quad))
        if recompute_velocity || recompute_velocity_gradient || recompute_induced_vorticity
            influence!(systems, systems, backend_system; precalc=false,
                scalar_potential=false,
                velocity=recompute_velocity,
                velocity_gradient=Tuple(recompute_velocity_gradient for _ in systems),
                extra_outputs=_induced_vorticity_extra_outputs(systems, recompute_induced_vorticity),
                direct_conditioning=_self_panel_core_size_conditioning())
        end
        if recompute_velocity
            apply_freestream!(systems, uinf)
            kinematic_velocity!(systems, frames)
        elseif recompute_potential
            for body in systems
                for i in axes(body.controlpoints, 2)
                    body.potential[i] += uinf[1] * body.controlpoints[1, i] +
                                         uinf[2] * body.controlpoints[2, i] +
                                         uinf[3] * body.controlpoints[3, i]
                end
            end
        end
        if !recompute_velocity && (recompute_velocity_gradient || recompute_induced_vorticity)
            _refresh_replay_kinematic_sidecars!(systems, frames)
        end
    end

    return nothing
end

function _missing_replay_fields(loaded::Vector{ReplayLoadedFields}, required::Set{Symbol})
    missing = String[]
    for (i, loaded_fields) in enumerate(loaded)
        fields = loaded_fields.fields
        for r in required
            r in fields || push!(missing, "body$(i):$(r)")
        end
    end
    return missing
end

function _replay_step_uinf(metadata, idx::Int, t, Uinf)
    uinf = _metadata_step_uinf(metadata, idx)
    uinf !== nothing && return uinf
    u = Uinf(t)
    length(u) == 3 || throw(ArgumentError("Replay Uinf fallback returned length $(length(u)); expected 3."))
    return FastMultipole.SVector{3, Float64}(u[1], u[2], u[3])
end

function _seed_replay_monitor_context!(ctx::MonitorContext, provided::Set{Symbol},
                                      loaded::Vector{ReplayLoadedFields}, fields::AbstractSet)
    for (i, lf) in enumerate(loaded)
        if lf.pressure !== nothing && !(:P in fields)
            monitor_register!(ctx, :P, i, lf.pressure)
            push!(provided, :P)
        end
        if lf.force !== nothing && !(:F in fields)
            monitor_register!(ctx, :F, i, lf.force)
            push!(provided, :F)
        end
    end
    return nothing
end

function _refresh_replay_kinematic_sidecars!(systems::Tuple, frames)
    saved_velocity = [copy(body.velocity) for body in systems]
    for body in systems
        body.velocity_kinematic .= 0
        body.angular_velocity .= 0
    end
    kinematic_velocity!(systems, frames)
    for (body, velocity) in zip(systems, saved_velocity)
        body.velocity .= velocity
    end
    return nothing
end

"""
    replay(path, run_name; monitors=(), monitor_factory=nothing, reconstruct=nothing,
           recompute=(:auto,), Uinf=t -> FastMultipole.SVector{3, Float64}(0.0, 0.0, 0.0),
           backend=FastMultipoleBackend(), backend_wake=backend,
           backend_system=backend, steps=:all, grad_mu_options=(;), verbose=false)

Load saved VTK body/wake/frame state and run monitors for selected saved steps
without solving, propagating wakes, shedding wakes, or advancing kinematics.
"""
function replay(path::AbstractString, run_name::AbstractString;
        monitors=(),
        monitor_factory=nothing,
        reconstruct=nothing,
        recompute=(:auto,),
        Uinf=t -> FastMultipole.SVector{3, Float64}(0.0, 0.0, 0.0),
        backend=FastMultipoleBackend(),
        backend_wake=backend,
        backend_system=backend,
        steps=:all,
        grad_mu_options=(;),
        verbose=false)
    timesteps, idxs, selected = _selected_replay_steps(path, run_name, steps)
    metadata = _read_metadata_toml(path, run_name)
    # BRAINSTORM 015: replay has no wake_attachment/kutta_closure plumbing, so
    # a run recorded with a non-default configuration would silently replay as
    # a legacy run and reconstruct wrong fields. Refuse until supported.
    if !isnothing(metadata) && haskey(metadata, "kutta")
        throw(ArgumentError("this run was recorded with a non-default "*
            "wake_attachment/kutta_closure configuration (manifest has a "*
            "[kutta] table); replay does not yet support it and would "*
            "silently reconstruct legacy fields."))
    end
    first_step = first(selected)
    systems, body_metadata = _read_body_metadata(path, run_name, first_step, metadata)
    wakes = _construct_wakes_from_manifest(systems, metadata)
    frames = _frames_from_metadata(path, run_name, first_step, systems, metadata)
    systems, wakes, frames = _apply_reconstruct_callback(reconstruct, path, run_name, body_metadata, systems, wakes, frames)
    if !isempty(body_metadata["missing"]) && isnothing(reconstruct)
        throw(ArgumentError("Replay cannot reconstruct all state from VTK: $(join(body_metadata["missing"], ", ")). Provide a replay manifest or reconstruct callback."))
    end

    if any(!isnothing, wakes)
        for (i, w) in enumerate(wakes)
            isnothing(w) && continue
            wake_name = run_name * "_wake$(i)"
            if w isa PanelParticleWake
                _load_replay_wake_step!(w, path, wake_name, first_step)
            elseif w isa PanelWake
                _load_replay_wake_step!(w, path, wake_name, first_step)
            else
                throw(ArgumentError("Replay does not support wake type $(typeof(w))."))
            end
        end
        _restore_wake_continuation!(wakes, metadata, first_step)
    end

    t_range = [timesteps[findfirst(==(idx), idxs)] for idx in selected]
    monitors = isnothing(monitor_factory) ? monitors : monitor_factory(systems, wakes, frames, t_range)
    req = _monitor_replay_requirements(monitors)
    # An unsteady PressureBernoulli converts the panel-following Dphi/Dt into an
    # Eulerian phi_dot using velocity_kinematic, which the loaded-fields fast
    # path never populates unless a recompute already refreshed it. Refresh the
    # kinematic sidecars explicitly in that case (see also the equivalent
    # trigger inside _recompute_replay_fields!).
    needs_kinematic_sidecars = any(
        m -> m isa PressureBernoulli && m.unsteady, monitors)
    rset = _recompute_set(recompute)
    auto = :auto in rset
    forced = auto ? Set{Symbol}() : rset

    for (out_i, idx) in enumerate(selected)
        loaded = ReplayLoadedFields[]
        for (i, body) in enumerate(systems)
            push!(loaded, _load_body_step!(body, path, run_name * "_body$(i)", idx))
        end
        for (i, w) in enumerate(wakes)
            isnothing(w) && continue
            wake_name = run_name * "_wake$(i)"
            if w isa PanelParticleWake
                _load_replay_wake_step!(w, path, wake_name, idx)
            elseif w isa PanelWake
                _load_replay_wake_step!(w, path, wake_name, idx)
            end
        end
        _restore_wake_continuation!(wakes, metadata, idx)
        frames = _frames_from_metadata(path, run_name, idx, systems, metadata)
        uinf = _replay_step_uinf(metadata, idx, t_range[out_i], Uinf)
        dt = out_i < length(selected) ? t_range[out_i + 1] - t_range[out_i] :
             out_i > 1 ? t_range[out_i] - t_range[out_i - 1] : 1.0

        needed = union(req, forced)
        if auto
            to_recompute = Set{Symbol}()
            for r in needed
                any(!(r in f.fields) for f in loaded) && push!(to_recompute, r)
            end
        else
            to_recompute = forced
            missing = _missing_replay_fields(loaded, setdiff(req, forced))
            isempty(missing) || throw(ArgumentError("Replay missing required monitor fields and recompute is disabled for them: $(join(missing, ", "))."))
        end

        _recompute_replay_fields!(systems, wakes, frames, uinf, to_recompute, backend_wake, backend_system;
            grad_mu_options)
        if needs_kinematic_sidecars && !(:velocity in to_recompute) &&
                !(:velocity_gradient in to_recompute) && !(:induced_vorticity in to_recompute)
            _refresh_replay_kinematic_sidecars!(systems, frames)
        end

        monitor_context = MonitorContext()
        provided = Set{Symbol}()
        _seed_replay_monitor_context!(monitor_context, provided, loaded, to_recompute)
        for (i, m) in enumerate(monitors)
            missing = String[]
            for r in monitor_requires(m)
                if !(Symbol(r) in provided)
                    push!(missing, String(r))
                end
            end
            isempty(missing) || throw(ArgumentError("Replay monitor $(nameof(typeof(m))) at position $(i) is missing fields $(join(missing, ", "))."))
            _run_monitor!(m, monitor_context, systems, wakes, frames, uinf, out_i - 1, dt)
            union!(provided, Symbol.(monitor_provides(m)))
        end
        verbose && println("replay: processed step $(idx)")
    end

    return ReplayResult(systems, wakes, frames, t_range, selected, monitors)
end
