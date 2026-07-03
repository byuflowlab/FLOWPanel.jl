#=##############################################################################
# DESCRIPTION
#   Unified restart/replay metadata helpers.
=###############################################################################

_metadata_toml_path(path, name) = joinpath(path, name * ".metadata.toml")
_replay_manifest_path(path, name) = joinpath(path, name * ".replay.toml")
_frames_toml_path(path, name) = joinpath(path, name * ".frames.toml")

_matrix_to_rows(A::AbstractMatrix) = [collect(A[:, j]) for j in axes(A, 2)]

_metadata_unsupported_dict(julia_type) = Dict{String, Any}(
    "type" => "unsupported",
    "julia_type" => string(julia_type),
    "restart_reconstruct_required" => true,
)

function _rows_to_matrix(rows, ::Type{T}=Int) where T
    if rows isa AbstractMatrix
        return T.(rows)
    elseif isempty(rows)
        return zeros(T, 0, 0)
    elseif first(rows) isa AbstractVector
        cols = [T.(r) for r in rows]
        M = hcat(cols...)
        return ndims(M) == 1 ? reshape(M, length(M), 1) : M
    else
        v = T.(rows)
        return reshape(v, length(v), 1)
    end
end

function _rows_to_matrix(rows, ::Type{T}, nrows::Int) where T
    if rows isa AbstractMatrix
        return T.(rows)
    elseif isempty(rows)
        return zeros(T, nrows, 0)
    elseif first(rows) isa AbstractVector
        cols = [T.(r) for r in rows]
        M = hcat(cols...)
        return ndims(M) == 1 ? reshape(M, length(M), 1) : M
    else
        v = T.(rows)
        length(v) % nrows == 0 || throw(ArgumentError("Cannot reshape $(length(v)) values into $(nrows) rows."))
        return reshape(v, nrows, :)
    end
end

_frame_static_dict(frame::ReferenceFrame, i::Int) = Dict{String, Any}(
    "i" => i,
    "name" => frame.name,
    "parent_index" => frame.parent_index,
    "child_index" => collect(frame.child_index),
    "dependent_index" => collect(frame.dependent_index),
)

_frame_state_dict(frame::ReferenceFrame, i::Int) = Dict{String, Any}(
    "i" => i,
    "x" => collect(frame.x),
    "v" => collect(frame.v),
    "omega" => frame.ω,
    "omega_axis" => collect(frame.ω_axis),
    "R" => collect(reshape(frame.R, 9)),
    "Rp2g" => collect(reshape(frame.Rp2g, 9)),
)

function _step_dict(frames, i_step::Int, t::Real; uinf=nothing)
    d = Dict{String, Any}(
        "i_step" => i_step,
        "t" => float(t),
        "frame" => [_frame_state_dict(frames[i], i) for i in eachindex(frames)],
    )
    if uinf !== nothing
        d["uinf"] = [float(uinf[1]), float(uinf[2]), float(uinf[3])]
    end
    return d
end

function _simulation_metadata_dict(t_range, start_step::Int, set_Das_eta_kinematic,
        set_Das_eta_freestream, set_Das_min_kinematic_displacement, clean_files::Bool;
        solver_options::NamedTuple=(;))
    n_steps = max(length(t_range) - 1, 0)
    dt = length(t_range) > 1 ? float(t_range[2] - t_range[1]) : NaN
    d = Dict{String, Any}(
        "t_start" => float(first(t_range)),
        "t_stop" => float(last(t_range)),
        "dt" => dt,
        "n_steps" => n_steps,
        "start_step" => start_step,
        "clean_files" => clean_files,
        "set_Das_eta_kinematic" => set_Das_eta_kinematic,
        "set_Das_eta_freestream" => set_Das_eta_freestream,
        "set_Das_min_kinematic_displacement" => set_Das_min_kinematic_displacement,
    )
    # Run-affecting simulate! toggles (coupling flags, relaxation on/off, body
    # strength relaxation). Scalars, recorded so a run can be reproduced.
    for (k, v) in pairs(solver_options)
        d[string(k)] = v
    end
    return d
end

function _backend_metadata_dict(backend)
    if backend isa Tuple
        return [_backend_metadata_dict(b) for b in backend]
    elseif backend isa DirectBackend
        return Dict{String, Any}("type" => "DirectBackend")
    elseif backend isa FastMultipoleBackend
        return Dict{String, Any}(
            "type" => "FastMultipoleBackend",
            "expansion_order" => backend.expansion_order,
            "multipole_acceptance" => backend.multipole_acceptance,
            "leaf_size" => backend.leaf_size,
        )
    elseif backend === nothing
        return Dict{String, Any}("type" => "nothing")
    else
        return _metadata_unsupported_dict(typeof(backend))
    end
end

function _solver_metadata_dict(solver)
    if solver isa Tuple
        return [_solver_metadata_dict(s) for s in solver]
    elseif solver isa Backslash
        return Dict{String, Any}("type" => "Backslash")
    elseif solver isa KrylovSolver
        return Dict{String, Any}(
            "type" => "KrylovSolver",
            "method" => String(solver.method),
            "itmax" => solver.itmax,
            "atol" => solver.atol,
            "rtol" => solver.rtol,
            "backend" => _backend_metadata_dict(solver.backend),
        )
    elseif solver isa FGSSolver
        return Dict{String, Any}(
            "type" => "FGSSolver",
            "expansion_order" => solver.expansion_order,
            "leaf_size" => solver.leaf_size,
            "multipole_acceptance" => solver.multipole_acceptance,
            "max_iterations" => solver.max_iterations,
            "inner_iterations" => solver.inner_iterations,
            "tolerance" => solver.tolerance,
            "rlx" => solver.rlx,
            "reverse_pass" => solver.reverse_pass,
            "verbose" => solver.verbose,
            "solution_history_length" => solver.solution_history_length,
            "project_solution" => solver.project_solution,
            "project_solution_order" => solver.project_solution_order,
        )
    elseif solver isa BackslashCoupled
        return Dict{String, Any}("type" => "BackslashCoupled")
    elseif solver === nothing
        return Dict{String, Any}("type" => "nothing")
    else
        return _metadata_unsupported_dict(typeof(solver))
    end
end

function _monitor_normalization_metadata(norm)
    if norm isa WingNormalization
        return Dict{String, Any}(
            "type" => "WingNormalization",
            "rho" => norm.rho,
            "Sref" => norm.Sref,
            "Lref" => norm.Lref,
        )
    elseif norm isa NoNormalization
        return Dict{String, Any}("type" => "NoNormalization")
    elseif norm isa RotorNormalization
        return Dict{String, Any}(
            "type" => "RotorNormalization",
            "rho" => norm.rho,
            "D" => norm.D,
            "i_frame" => norm.i_frame,
        )
    elseif norm isa RotorNormalization2
        return Dict{String, Any}(
            "type" => "RotorNormalization2",
            "rho" => norm.rho,
            "D" => norm.D,
            "i_frame" => norm.i_frame,
        )
    elseif norm isa NoSectionalNormalization
        return Dict{String, Any}("type" => "NoSectionalNormalization")
    elseif norm isa FreestreamSectionalNormalization
        return Dict{String, Any}(
            "type" => "FreestreamSectionalNormalization",
            "rho" => norm.rho,
            "Lref" => norm.Lref,
        )
    elseif norm isa RotorSectionalNormalization
        return Dict{String, Any}(
            "type" => "RotorSectionalNormalization",
            "rho" => norm.rho,
            "R" => norm.R,
            "i_frame" => norm.i_frame,
            "omega_scale" => string(norm.omega_scale),
        )
    else
        return _metadata_unsupported_dict(typeof(norm))
    end
end

function _monitor_metadata(m)
    if m isa PressureBernoulli
        return Dict{String, Any}(
            "type" => "PressureBernoulli",
            "rho" => m.rho,
            "unsteady" => m.unsteady,
            "correct_kuttacondition" => m.correct_kuttacondition,
            "clip" => m.clip === nothing ? "nothing" : string(m.clip),
            "backend" => _backend_metadata_dict(m.backend),
            "vtk_fields" => collect(string.(m.vtk_fields)),
            "file" => m.file,
        )
    elseif m isa PressureLaplace
        return Dict{String, Any}(
            "type" => "PressureLaplace",
            "rho" => m.rho,
            "unsteady" => m.unsteady,
            "atol" => m.atol,
            "rtol" => m.rtol,
            "itmax" => m.itmax,
            "preconditioner" => _metadata_unsupported_dict(typeof(m.preconditioner)),
            "reference_panel" => m.reference_panel,
            "reference_pressure" => m.reference_pressure,
            "rebuild_every_step" => m.rebuild_every_step,
            "verbose" => m.verbose,
            "gradient_mode" => string(m.gradient_mode),
            "acceleration_form" => string(m.acceleration_form),
            "vtk_fields" => collect(string.(m.vtk_fields)),
            "file" => m.file,
        )
    elseif m isa ForceMonitor
        return Dict{String, Any}(
            "type" => "ForceMonitor",
            "i_system" => m.i_system,
            "i_frame" => m.i_frame,
            "normalization" => _monitor_normalization_metadata(m.normalization),
            "correct_kuttacondition" => m.correct_kuttacondition,
            "verbose" => m.verbose,
            "vtk_fields" => collect(string.(m.vtk_fields)),
            "file" => m.file,
        )
    elseif m isa SurfaceVorticityForce
        return Dict{String, Any}(
            "type" => "SurfaceVorticityForce",
            "i_system" => m.i_system,
            "i_frame" => m.i_frame,
            "rho" => m.rho,
            "normalization" => _monitor_normalization_metadata(m.normalization),
            "correct_kuttacondition" => m.correct_kuttacondition,
            "grad_mu_options" => Dict(string(k) => string(v) for (k, v) in pairs(m.grad_mu_options)),
            "verbose" => m.verbose,
            "vtk_fields" => collect(string.(m.vtk_fields)),
            "file" => m.file,
        )
    elseif m isa SpanwiseLoadingMonitor
        return Dict{String, Any}(
            "type" => "SpanwiseLoadingMonitor",
            "nbins" => m.nbins,
            "i_system" => m.i_system,
            "i_frame" => m.i_frame,
            "component_names" => collect(string.(m.component_names)),
            "span_axis" => collect(m.span_axis),
            "normalization" => _monitor_normalization_metadata(m.normalization),
            "per_length" => m.per_length,
            "binning" => string(m.binning),
            "file" => m.file,
            "verbose" => m.verbose,
            "vtk_fields" => collect(string.(m.vtk_fields)),
        )
    elseif m isa KuttaJoukowskiForce
        return Dict{String, Any}(
            "type" => "KuttaJoukowskiForce",
            "i_system" => m.i_system,
            "i_frame" => m.i_frame,
            "rho" => m.rho,
            "backend" => _backend_metadata_dict(m.backend),
            "normalization" => _monitor_normalization_metadata(m.normalization),
            "verbose" => m.verbose,
            "file" => m.file,
        )
    elseif m isa BoundCirculationMonitor
        return Dict{String, Any}(
            "type" => "BoundCirculationMonitor",
            "i_system" => m.i_system,
            "i_frame" => m.i_frame,
            "radial_dimension" => m.radial_dimension,
            "R" => m.radius,
            "section_tol" => m.section_tol === nothing ? "nothing" : m.section_tol,
            "verbose" => m.verbose,
            "file" => m.file,
        )
    elseif m isa CylindricalFieldProbeMonitor
        return Dict{String, Any}(
            "type" => "CylindricalFieldProbeMonitor",
            "i_frame" => m.i_frame,
            "axial_axis" => string(m.axial_axis),
            "n_axial" => m.n_axial,
            "n_radial" => m.n_radial,
            "n_azimuth" => m.n_azimuth,
            "save_path" => m.save_path === nothing ? "nothing" : m.save_path,
            "name" => m.name,
            "verbose" => m.verbose,
        )
    else
        return _metadata_unsupported_dict(typeof(m))
    end
end

function _metadata_manifest_dict(name, systems::Tuple, wakes::Tuple, frames,
        t_range, body_solvers, backend_wake, backend_solve, backend_system,
        monitors; start_step::Int=0,
        set_Das_eta_kinematic=NaN,
        set_Das_eta_freestream=NaN,
        set_Das_min_kinematic_displacement=0.0,
        clean_files::Bool=true,
        solver_options::NamedTuple=(;))
    return Dict{String, Any}(
        "meta" => Dict{String, Any}(
            "schema_version" => 1,
            "purpose" => "FLOWPanel restart/replay metadata",
            "run_name" => name,
            "n_bodies" => length(systems),
            "n_wakes" => length(wakes),
            "n_frames" => length(frames),
        ),
        "simulation" => _simulation_metadata_dict(t_range, start_step,
            set_Das_eta_kinematic, set_Das_eta_freestream,
            set_Das_min_kinematic_displacement, clean_files; solver_options),
        "backends" => Dict{String, Any}(
            "wake" => _backend_metadata_dict(backend_wake),
            "system" => _backend_metadata_dict(backend_system),
            "solve" => _backend_metadata_dict(backend_solve),
            "body_solvers" => _solver_metadata_dict(body_solvers),
        ),
        "body" => [_body_manifest_dict(body, i) for (i, body) in enumerate(systems)],
        "wake" => [_wake_manifest_dict(wake, i) for (i, wake) in enumerate(wakes)],
        "frame" => [_frame_static_dict(frame, i) for (i, frame) in enumerate(frames)],
        "monitor" => [_monitor_metadata(m) for m in monitors],
    )
end

function _write_metadata_toml(path, name, systems::Tuple, wakes::Tuple, frames,
        t_range, body_solvers, backend_wake, backend_solve, backend_system,
        monitors; start_step::Int=0,
        set_Das_eta_kinematic=NaN,
        set_Das_eta_freestream=NaN,
        set_Das_min_kinematic_displacement=0.0,
        clean_files::Bool=true,
        solver_options::NamedTuple=(;))
    mkpath(path)
    file = _metadata_toml_path(path, name)
    open(file, "w") do io
        TOML.print(io, _metadata_manifest_dict(name, systems, wakes, frames,
            t_range, body_solvers, backend_wake, backend_solve, backend_system,
            monitors; start_step, set_Das_eta_kinematic, set_Das_eta_freestream,
            set_Das_min_kinematic_displacement, clean_files, solver_options))
    end
    return file
end

function _append_metadata_step_toml(path, name, frames, i_step::Int, t::Real; uinf=nothing)
    file = _metadata_toml_path(path, name)
    step = _step_dict(frames, i_step, t; uinf)
    data = isfile(file) ? TOML.parsefile(file) : Dict{String, Any}()
    steps = get(data, "step", Any[])
    existing = findfirst(s -> Int(s["i_step"]) == i_step, steps)
    if isnothing(existing)
        push!(steps, step)
    else
        steps[existing] = step
    end
    data["step"] = steps
    open(file, "w") do io
        TOML.print(io, data)
    end
    return file
end

function _read_metadata_toml(path, name)
    file = _metadata_toml_path(path, name)
    if isfile(file)
        return TOML.parsefile(file)
    end

    legacy = _replay_manifest_path(path, name)
    if isfile(legacy)
        data = TOML.parsefile(legacy)
        frame_file = _frames_toml_path(path, name)
        if isfile(frame_file)
            data["step"] = get(TOML.parsefile(frame_file), "step", Any[])
        end
        return data
    end

    return nothing
end

function _metadata_step_frames(data, restart_step::Int)
    data === nothing && return nothing
    steps = get(data, "step", Any[])
    idx_in_toml = findfirst(s -> Int(s["i_step"]) == restart_step, steps)
    isnothing(idx_in_toml) && return nothing
    return steps[idx_in_toml]["frame"]
end

function _metadata_step_uinf(data, restart_step::Int)
    data === nothing && return nothing
    steps = get(data, "step", Any[])
    idx_in_toml = findfirst(s -> Int(s["i_step"]) == restart_step, steps)
    isnothing(idx_in_toml) && return nothing
    step = steps[idx_in_toml]
    haskey(step, "uinf") || return nothing
    u = step["uinf"]
    length(u) == 3 || throw(ArgumentError("Metadata step $(restart_step) has uinf with length $(length(u)); expected 3."))
    return FastMultipole.SVector{3, Float64}(u[1], u[2], u[3])
end
