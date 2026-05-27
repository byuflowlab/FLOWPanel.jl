#=##############################################################################
# DESCRIPTION
#   Replay saved VTK simulation output through post-processing monitors without
#   advancing the solver, wake, or kinematics.
=###############################################################################

import ReadVTK

struct ReplayResult{TS,TW,TF,TR,TM}
    systems::TS
    wakes::TW
    frames::TF
    t_range::TR
    steps::Vector{Int}
    monitors::TM
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
    else
        return _metadata_unsupported_dict(typeof(method))
    end
end

function _particle_policy_manifest(policy)
    if policy isa MinGamma
        return Dict{String, Any}("type" => "MinGamma", "threshold" => policy.threshold)
    elseif policy isa MaxGamma
        return Dict{String, Any}("type" => "MaxGamma", "threshold" => policy.threshold)
    elseif policy isa GlobalBox
        return Dict{String, Any}(
            "type" => "GlobalBox",
            "xmin" => collect(policy.xmin),
            "xmax" => collect(policy.xmax),
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
    else
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

function _viscous_manifest(viscous)
    if viscous isa FLOWVPM.Inviscid
        return Dict{String, Any}("type" => "FLOWVPM.Inviscid", "nu" => viscous.nu)
    elseif viscous isa FLOWVPM.CoreSpreading
        return Dict{String, Any}(
            "type" => "FLOWVPM.CoreSpreading",
            "nu" => viscous.nu,
            "core_size" => viscous.sgm0,
            "kernel" => "FLOWVPM.zeta_fmm",
            "beta" => viscous.beta,
        )
    elseif viscous === nothing
        return Dict{String, Any}("type" => "nothing")
    else
        return _metadata_unsupported_dict(typeof(viscous))
    end
end

function _sfs_manifest(sfs)
    if sfs === FLOWVPM.SFS_default
        return Dict{String, Any}("type" => "FLOWVPM.SFS_default")
    elseif sfs === FLOWVPM.SFS_Cd_twolevel_nobackscatter
        return Dict{String, Any}("type" => "FLOWVPM.SFS_Cd_twolevel_nobackscatter")
    elseif sfs === nothing
        return Dict{String, Any}("type" => "nothing")
    else
        return _metadata_unsupported_dict(typeof(sfs))
    end
end

function _body_manifest_dict(body::AbstractBody, i::Int)
    d = Dict{String, Any}(
        "i" => i,
        "kind" => _body_kind_string(body),
        "strength_names" => collect(strength_names(body)),
        "dbc" => has_dirichlet_bc(body),
        "cp_offset" => body.CPoffset,
        "kerneloffset_panel" => body.kerneloffset_panel,
        "kerneloffset_targets" => body.kerneloffset_targets,
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

function _wake_manifest_dict(wake, i::Int)
    d = Dict{String, Any}("i" => i)
    if isnothing(wake)
        d["type"] = "nothing"
    elseif wake isa PanelParticleWake
        d["type"] = "PanelParticleWake"
        d["nwakerows"] = size(wake.panel_wake.nodes[1], 2) - 1
        d["max_particles"] = size(wake.pfield.particles, 2)
        d["core_size"] = wake.panel_wake.core_size
        d["shed_with_induced_velocity"] = wake.panel_wake.shed_with_induced_velocity
        d["unsteady_filament"] = wake.panel_wake.unsteady_filament
        d["particle_kerneloffset"] = wake.particle_kerneloffset
        d["method_trailing"] = _wake_shedding_manifest(wake.method_trailing)
        d["method_unsteady"] = _wake_shedding_manifest(wake.method_unsteady)
        d["particle_maintenance"] = _particle_maintenance_manifest(wake.particle_maintenance)
        d["viscous"] = _viscous_manifest(wake.pfield.viscous)
        d["SFS"] = _sfs_manifest(wake.pfield.SFS)
    elseif wake isa PanelWake
        d["type"] = "PanelWake"
        d["nwakerows"] = size(wake.nodes[1], 2) - 1
        d["core_size"] = wake.core_size
        d["shed_with_induced_velocity"] = wake.shed_with_induced_velocity
        d["unsteady_filament"] = wake.unsteady_filament
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
        CPoffset=Float64(get(body_meta, "cp_offset", E === VortexRing ? 1e-6 : 1e-14)),
        kerneloffset_panel=Float64(get(body_meta, "kerneloffset_panel", get(body_meta, "kerneloffset", 1e-8))),
        kerneloffset_targets=Float64(get(body_meta, "kerneloffset_targets", get(body_meta, "kerneloffset", 1e-8))),
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
        body.P .= ReadVTK.get_data(cell_data["gauge pressure"])
        push!(loaded, :P)
    end
    if "F" in keys(cell_data)
        body.F .= ReadVTK.get_data(cell_data["F"])
        push!(loaded, :F)
    end
    if "normals" in keys(cell_data)
        body.normals .= ReadVTK.get_data(cell_data["normals"])
        push!(loaded, :normals)
    end
    return loaded
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
    calc_controlpoints!(body; off=abs(body.CPoffset))
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
                unsteady_filament=Bool(get(wmeta, "unsteady_filament", true))))
        elseif wtype == "PanelParticleWake"
            systems[i] isa AbstractLiftingBody ||
                throw(ArgumentError("Cannot reconstruct PanelParticleWake for non-lifting body $(i)."))
            method_trailing = _deserialize_wake_shedding(get(wmeta, "method_trailing", Dict{String, Any}("type" => "OverlapPPS")))
            method_unsteady = _deserialize_wake_shedding(get(wmeta, "method_unsteady", Dict{String, Any}("type" => "OverlapPPS")))
            particle_maintenance = _deserialize_particle_maintenance(get(wmeta, "particle_maintenance", Dict{String, Any}("type" => "ParticleMaintenance")))
            viscous = _deserialize_viscous(get(wmeta, "viscous", Dict{String, Any}("type" => "FLOWVPM.Inviscid")))
            sfs = _deserialize_sfs(get(wmeta, "SFS", Dict{String, Any}("type" => "FLOWVPM.SFS_default")))
            push!(wakes, PanelParticleWake(systems[i];
                nwakerows=Int(get(wmeta, "nwakerows", 3)),
                max_particles=Int(get(wmeta, "max_particles", 10000)),
                core_size=Float64(get(wmeta, "core_size", 1e-3)),
                shed_with_induced_velocity=Bool(get(wmeta, "shed_with_induced_velocity", true)),
                unsteady_filament=Bool(get(wmeta, "unsteady_filament", true)),
                method_trailing=method_trailing,
                method_unsteady=method_unsteady,
                particle_maintenance=particle_maintenance,
                particle_kerneloffset=Float64(get(wmeta, "particle_kerneloffset", NaN)),
                viscous=viscous,
                SFS=sfs))
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
    elseif wtype == "nothing"
        return NoShed()
    else
        return NoShed()
    end
end

function _deserialize_particle_policy(meta)
    ptype = String(get(meta, "type", "unsupported"))
    if ptype == "MinGamma"
        return MinGamma(Float64(get(meta, "threshold", 0.0)))
    elseif ptype == "MaxGamma"
        return MaxGamma(Float64(get(meta, "threshold", 0.0)))
    elseif ptype == "GlobalBox"
        return GlobalBox(get(meta, "xmin", [0.0, 0.0, 0.0]), get(meta, "xmax", [0.0, 0.0, 0.0]))
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
        kernel = get(meta, "kernel", "FLOWVPM.zeta_fmm")
        kernel_fn = kernel == "FLOWVPM.zeta_fmm" ? FLOWVPM.zeta_fmm : FLOWVPM.zeta_fmm
        return FLOWVPM.CoreSpreading(
            Float64(get(meta, "nu", 0.0)),
            Float64(get(meta, "core_size", 0.0)),
            kernel_fn;
            beta=Float64(get(meta, "beta", 1.5)),
        )
    elseif vtype == "nothing"
        return FLOWVPM.Inviscid()
    else
        return FLOWVPM.Inviscid()
    end
end

function _deserialize_sfs(meta)
    stype = String(get(meta, "type", "unsupported"))
    if stype == "FLOWVPM.SFS_default"
        return FLOWVPM.SFS_default
    elseif stype == "FLOWVPM.SFS_Cd_twolevel_nobackscatter"
        return FLOWVPM.SFS_Cd_twolevel_nobackscatter
    elseif stype == "nothing"
        return FLOWVPM.SFS_default
    else
        return FLOWVPM.SFS_default
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
    recompute === :all && return Set([:velocity, :potential, :velocity_gradient, :induced_vorticity, :P, :F])
    recompute isa Tuple || throw(ArgumentError("recompute must be a tuple of symbols, (:auto,), (:all,), (), or :all."))
    :all in recompute && return Set([:velocity, :potential, :velocity_gradient, :induced_vorticity, :P, :F])
    :auto in recompute && length(recompute) == 1 && return Set{Symbol}([:auto])
    return Set(Symbol.(recompute))
end

function _monitor_replay_requirements(monitors)
    req = Set{Symbol}()
    for m in monitors
        for r in monitor_requires(m)
            push!(req, Symbol(r))
        end
        if m isa PressureBernoulli
            push!(req, :velocity)
            m.unsteady && push!(req, :potential)
        elseif m isa PressureLaplace
            push!(req, :velocity)
            monitor_requires_body_hessian(m) && push!(req, :velocity_gradient)
            monitor_requires_induced_vorticity(m) && push!(req, :induced_vorticity)
        end
    end
    return req
end

function _recompute_replay_fields!(systems::Tuple, wakes::Tuple, frames, uinf, fields::AbstractSet, backend_wake, backend_system)
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
        _set_kerneloffsets!(systems, :kerneloffset_targets)
        scalar_sources = _filter_scalar_potential_sources(systems)
        if recompute_potential && length(scalar_sources) > 0
            influence!(systems, scalar_sources, backend_system; scalar_potential=true, velocity=false)
        end
        if recompute_velocity || recompute_velocity_gradient || recompute_induced_vorticity
            influence!(systems, systems, backend_system; precalc=false,
                scalar_potential=false,
                velocity=recompute_velocity,
                velocity_gradient=Tuple(recompute_velocity_gradient for _ in systems),
                extra_outputs=_induced_vorticity_extra_outputs(systems, recompute_induced_vorticity),
                direct_conditioning=_self_panel_kerneloffset_conditioning())
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
    if :F in fields
        for body in systems
            fill!(body.F, zero(eltype(body.F)))
        end
        calcfield_F!(systems)
    end
    return nothing
end

function _missing_replay_fields(loaded::Vector{Set{Symbol}}, required::Set{Symbol})
    missing = String[]
    for (i, fields) in enumerate(loaded)
        for r in required
            r in fields || push!(missing, "body$(i):$(r)")
        end
    end
    return missing
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
           recompute=(:auto,), backend=FastMultipoleBackend(), backend_wake=backend,
           backend_system=backend, steps=:all, verbose=false)

Load saved VTK body/wake/frame state and run monitors for selected saved steps
without solving, propagating wakes, shedding wakes, or advancing kinematics.
"""
function replay(path::AbstractString, run_name::AbstractString;
        monitors=(),
        monitor_factory=nothing,
        reconstruct=nothing,
        recompute=(:auto,),
        backend=FastMultipoleBackend(),
        backend_wake=backend,
        backend_system=backend,
        steps=:all,
        verbose=false)
    timesteps, idxs, selected = _selected_replay_steps(path, run_name, steps)
    metadata = _read_metadata_toml(path, run_name)
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
    end

    t_range = [timesteps[findfirst(==(idx), idxs)] for idx in selected]
    monitors = isnothing(monitor_factory) ? monitors : monitor_factory(systems, wakes, frames, t_range)
    req = _monitor_replay_requirements(monitors)
    rset = _recompute_set(recompute)
    auto = :auto in rset
    forced = auto ? Set{Symbol}() : rset

    for (out_i, idx) in enumerate(selected)
        loaded = Set{Symbol}[]
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
        frames = _frames_from_metadata(path, run_name, idx, systems, metadata)
        uinf = FastMultipole.SVector{3, Float64}(0.0, 0.0, 0.0)
        dt = out_i < length(selected) ? t_range[out_i + 1] - t_range[out_i] :
             out_i > 1 ? t_range[out_i] - t_range[out_i - 1] : 1.0

        needed = union(req, forced)
        if auto
            to_recompute = Set{Symbol}()
            for r in needed
                any(!(r in f) for f in loaded) && push!(to_recompute, r)
            end
        else
            to_recompute = forced
            missing = _missing_replay_fields(loaded, setdiff(req, forced))
            isempty(missing) || throw(ArgumentError("Replay missing required monitor fields and recompute is disabled for them: $(join(missing, ", "))."))
        end

        (:P in to_recompute) && throw(ArgumentError("Replay cannot generically recompute :P without a pressure monitor; include PressureBernoulli or PressureLaplace before consumers, or load gauge pressure from VTK."))
        _recompute_replay_fields!(systems, wakes, frames, uinf, to_recompute, backend_wake, backend_system)

        provided = Set{Symbol}()
        for (i, m) in enumerate(monitors)
            missing = String[]
            for r in monitor_requires(m)
                if !(Symbol(r) in provided) && any(!(Symbol(r) in f) for f in loaded) && !(Symbol(r) in to_recompute)
                    push!(missing, String(r))
                end
            end
            isempty(missing) || throw(ArgumentError("Replay monitor $(nameof(typeof(m))) at position $(i) is missing fields $(join(missing, ", "))."))
            m(systems, wakes, frames, uinf, out_i - 1, dt)
            union!(provided, Symbol.(monitor_provides(m)))
        end
        verbose && println("replay: processed step $(idx)")
    end

    return ReplayResult(systems, wakes, frames, t_range, selected, monitors)
end
