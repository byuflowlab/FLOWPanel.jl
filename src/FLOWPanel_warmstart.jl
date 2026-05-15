#=##############################################################################
# DESCRIPTION
    Warm-start support for `simulate!`. Resume a time-marching panel-method
    simulation from VTK files produced by a previous run, plus a small companion
    `{name}.frames.toml` that captures `ReferenceFrame` state (which the VTK
    output cannot represent).
=###############################################################################

import ReadVTK
import TOML

#--- frame-state TOML helpers ---#

function _frame_to_dict(frame::ReferenceFrame, i::Int)
    return Dict{String, Any}(
        "i"          => i,
        "name"       => frame.name,
        "x"          => collect(frame.x),
        "v"          => collect(frame.v),
        "omega"      => frame.ω,
        "omega_axis" => collect(frame.ω_axis),
        "R"          => collect(reshape(frame.R, 9)),
        "Rp2g"       => collect(reshape(frame.Rp2g, 9)),
    )
end

function _step_dict(frames, i_step::Int, t::Real)
    return Dict{String, Any}(
        "i_step" => i_step,
        "t"      => float(t),
        "frame"  => [_frame_to_dict(frames[i], i) for i in eachindex(frames)],
    )
end

function _frames_toml_path(path, name)
    return joinpath(path, name * ".frames.toml")
end

"""
Write or append the per-step ReferenceFrame state to `{path}/{name}.frames.toml`.
At `truncate=true`, the file is overwritten with `[meta]` and the first
`[[step]]` block; otherwise a single `[[step]]` block is appended.
"""
function _write_frame_state_toml(path, name, frames, i_step::Int, t::Real; truncate::Bool=false)
    isnothing(frames) && return nothing

    file = _frames_toml_path(path, name)
    if truncate
        open(file, "w") do io
            meta = Dict{String, Any}(
                "sim_name"    => name,
                "n_frames"    => length(frames),
                "frame_names" => [f.name for f in frames],
            )
            TOML.print(io, Dict("meta" => meta))
            println(io)
            TOML.print(io, Dict("step" => [_step_dict(frames, i_step, t)]))
        end
    else
        open(file, "a") do io
            TOML.print(io, Dict("step" => [_step_dict(frames, i_step, t)]))
        end
    end
    return file
end

"""
Load `frames` state from `{path}/{name}.frames.toml` for the entry whose
`i_step == restart_step`. Replaces each `frames[i]` in place.
"""
function _load_frame_state_toml!(frames, path, name, restart_step::Int)
    isnothing(frames) && error("simulate_warmstart! requires frames; ReferenceFrame state is not available from VTK output alone.")

    file = _frames_toml_path(path, name)
    isfile(file) || error("Warm-start frames file not found: $(file). The original simulate! run must produce this file.")

    data = TOML.parsefile(file)
    steps = data["step"]
    idx_in_toml = findfirst(s -> Int(s["i_step"]) == restart_step, steps)
    isnothing(idx_in_toml) && error("restart_step=$(restart_step) not found in $(file)")
    step = steps[idx_in_toml]

    framedicts = step["frame"]
    @assert length(framedicts) == length(frames) "frames file has $(length(framedicts)) frames, but `frames` argument has $(length(frames))"

    for (i, fd) in enumerate(framedicts)
        old = frames[i]
        TF = eltype(old.x)
        x          = FastMultipole.SVector{3,TF}(fd["x"]...)
        v          = FastMultipole.SVector{3,TF}(fd["v"]...)
        ω_axis     = FastMultipole.SVector{3,TF}(fd["omega_axis"]...)
        ω          = TF(fd["omega"])
        Rvec       = fd["R"]
        Rp2gvec    = fd["Rp2g"]
        R          = FastMultipole.SMatrix{3,3,TF,9}(Rvec...)
        Rp2g       = FastMultipole.SMatrix{3,3,TF,9}(Rp2gvec...)
        frames[i] = ReferenceFrame{TF}(x, v, ω_axis, ω, R, Rp2g,
            old.name, old.parent_index, old.child_index, old.dependent_index)
    end
    return frames
end

#--- PVD parsing ---#

"""
Read a `.pvd` collection and return `(timesteps::Vector{Float64}, idxs::Vector{Int})`,
where `idxs[k]` is the integer step index parsed from the k-th entry's filename.
"""
function _read_pvd_steps(pvd_path::String)
    isfile(pvd_path) || error("PVD file not found: $(pvd_path)")
    pvd = ReadVTK.PVDFile(pvd_path)
    idxs = Int[]
    for fname in pvd.vtk_filenames
        # filename like ".../{base}.{idx}.vtm" → match the integer between the last two dots
        m = match(r"\.(\d+)\.[a-zA-Z]+$", fname)
        isnothing(m) && error("Could not parse step index from VTK filename: $(fname)")
        push!(idxs, parse(Int, m.captures[1]))
    end
    return pvd.timesteps, idxs
end

#--- Body VTK loading ---#

"""
Load body state at step `idx` from `{path}/{name}/{name}.{idx}.vtu` into the
already-constructed `body`. Updates `body.nodes`, `body.strength`, and
`body.potential`. Velocity / pressure / forces are not restored — they are
recomputed at the start of the next step.
"""
function _load_body_vtk!(body::AbstractBody, path::String, name::String, idx::Int)
    vtu_path = joinpath(path, name, "$(name).$(idx).vtu")
    isfile(vtu_path) || error("Body VTU not found: $(vtu_path)")

    vtk = ReadVTK.VTKFile(vtu_path)

    # nodes (3 × n_points)
    points = ReadVTK.get_points(vtk)
    @assert size(points) == size(body.nodes) "Loaded VTU points size $(size(points)) != body.nodes size $(size(body.nodes))"
    body.nodes .= points

    # cell data
    cell_data = ReadVTK.get_cell_data(vtk)

    # strengths (one cell-data array per element type)
    for (i, sname) in enumerate(strength_names(body))
        sname in keys(cell_data) || error("Strength field '$(sname)' missing from $(vtu_path)")
        arr = ReadVTK.get_data(cell_data[sname])
        body.strength[:, i] .= arr
    end

    # potential (cell)
    if "potential" in keys(cell_data)
        body.potential .= ReadVTK.get_data(cell_data["potential"])
    end

    return body
end

#--- Panel wake VTK loading ---#

"""
Load `PanelWake` state at step `idx` from `{path}/{wake_name}/{wake_name}.{i_surf}.{idx}.vts`
files (one per surface). Restores `wake.nodes`, `wake.strength`, `wake.velocity`,
and `wake.nwakes[]`. `wake.overflowed[]` is set when the loaded panel buffer is
at the maximum row count.
"""
function _load_panel_wake_vtk!(wake::PanelWake, path::String, wake_name::String, idx::Int)
    nwakes_loaded = -1
    for i_surf in eachindex(wake.nodes)
        vts_path = joinpath(path, wake_name, "$(wake_name).$(i_surf).$(idx).vts")
        isfile(vts_path) || error("Wake VTS not found: $(vts_path)")
        vtk = ReadVTK.VTKFile(vts_path)

        # points: structured grid of shape (3, nrows+1, ncols+1, 1) -> reshape
        points = ReadVTK.get_points(vtk)
        # ReadVTK returns points as a 3 x n_points matrix; we need to map to
        # (3, nrows+1, ncols+1) using whole_extent.
        dims, _ = ReadVTK.get_wholeextent(vtk.xml_file)
        dim1, dim2, dim3 = dims

        @assert dim3 == 1 "Wake VTS expected to have third structured dim of 1, got $(dim3)"

        nwakes_this = dim1 - 1
        if nwakes_loaded < 0
            nwakes_loaded = nwakes_this
        else
            @assert nwakes_loaded == nwakes_this "Inconsistent nwakes across wake surfaces ($(nwakes_loaded) vs $(nwakes_this))"
        end

        # reshape points (3 × n_points) into (3, dim1, dim2)
        nodes_flat = reshape(points, 3, dim1, dim2)
        @assert size(wake.nodes[i_surf], 3) == dim2 "Wake VTS surface $(i_surf) has $(dim2) cols, but wake.nodes has $(size(wake.nodes[i_surf], 3))"
        wake.nodes[i_surf][:, 1:dim1, :] .= nodes_flat

        # velocity (point data)
        point_data = ReadVTK.get_point_data(vtk)
        if "velocity" in keys(point_data)
            vel_arr = ReadVTK.get_data(point_data["velocity"])
            vel_flat = reshape(vel_arr, 3, dim1, dim2)
            wake.velocity[i_surf][:, 1:dim1, :] .= vel_flat
        end

        # strength (cell data) — shape (dim_strength, dim1-1, dim2-1)
        cell_data = ReadVTK.get_cell_data(vtk)
        if "strength" in keys(cell_data)
            str_arr = ReadVTK.get_data(cell_data["strength"])
            dim_strength = size(wake.strength[i_surf], 1)
            str_flat = reshape(str_arr, dim_strength, dim1-1, dim2-1)
            wake.strength[i_surf][:, 1:dim1-1, :] .= str_flat
        end
    end

    wake.nwakes[] = max(nwakes_loaded, 0)
    n_rows_max = size(wake.nodes[1], 2)
    wake.overflowed[] = (wake.nwakes[] >= n_rows_max - 1)
    return wake
end

#--- Particle wake VTK loading ---#

"""
Load `PanelParticleWake` state at step `idx`. First loads the panel-wake
component, then loads the particle field from
`{path}/{wake_name}_particles/{wake_name}_particles.{idx}.vtp`.
"""
function _load_panel_particle_wake_vtk!(wake::PanelParticleWake, path::String, wake_name::String, idx::Int)
    _load_panel_wake_vtk!(wake.panel_wake, path, wake_name, idx)

    vtp_dir  = joinpath(path, wake_name * "_particles")
    vtp_path = joinpath(vtp_dir, "$(wake_name)_particles.$(idx).vtp")
    isfile(vtp_path) || error("Particles VTP not found: $(vtp_path)")
    vtk = ReadVTK.VTKFile(vtp_path)

    # number of particles
    np = vtk.n_points
    pf = wake.pfield
    @assert np <= size(pf.particles, 2) "Loaded $(np) particles but pfield has capacity $(size(pf.particles, 2))"

    # first reset any existing particles
    while pf.np > 0
        FLOWVPM.remove_particle(pf, pf.np)
    end

    pf.np = 0
    if np == 0
        return wake
    end

    # positions
    points = ReadVTK.get_points(vtk)  # 3 × np
    pf.particles[FLOWVPM.X_INDEX, 1:np] .= points

    # per-particle fields
    point_data = ReadVTK.get_point_data(vtk)
    if "gamma" in keys(point_data)
        pf.particles[FLOWVPM.GAMMA_INDEX, 1:np] .= ReadVTK.get_data(point_data["gamma"])
    end
    if "sigma" in keys(point_data)
        pf.particles[FLOWVPM.SIGMA_INDEX, 1:np] .= ReadVTK.get_data(point_data["sigma"])
    end
    if "vol" in keys(point_data)
        pf.particles[FLOWVPM.VOL_INDEX, 1:np] .= ReadVTK.get_data(point_data["vol"])
    end
    if "circulation" in keys(point_data)
        pf.particles[FLOWVPM.CIRCULATION_INDEX, 1:np] .= ReadVTK.get_data(point_data["circulation"])
    end
    if "velocity" in keys(point_data)
        pf.particles[FLOWVPM.U_INDEX, 1:np] .= ReadVTK.get_data(point_data["velocity"])
    end
    if "vorticity" in keys(point_data)
        pf.particles[FLOWVPM.VORTICITY_INDEX, 1:np] .= ReadVTK.get_data(point_data["vorticity"])
    end
    if "velocity_gradient" in keys(point_data)
        J_arr = ReadVTK.get_data(point_data["velocity_gradient"])
        # written as reshape(view(..., J_INDEX, 1:np), 3, 3, np); ReadVTK gives back as 9 × np or similar.
        pf.particles[FLOWVPM.J_INDEX, 1:np] .= reshape(J_arr, 9, np)
    end

    pf.np = np
    return wake
end

#--- public API ---#

"""
    simulate_warmstart!(systems, wakes, frames, maneuver!, Uinf, t_range; restart_path, restart_name, restart_step=-1, body_solvers, backend_wake=backend, backend_solve=backend, backend_system=backend, optargs...)

Resume a simulation from VTK output and a `{name}.frames.toml` companion file
written by a previous `simulate!` call. The body and wake objects must be
constructed identically to the original run; their state is overwritten from
disk.

Reads the last entry of `{restart_path}/{restart_name}_body1.pvd` (or the entry
matching `restart_step` if given), restores body + wake + frame state, replays
the end-of-step propagate / kinematics / shed sequence (which `simulate!` skips
on the final step of a run), then forwards to `simulate!` with
`start_step=restart_step+1` so that subsequent steps append to the same VTK and
TOML output.
"""
function simulate_warmstart!(systems, wakes, frames, maneuver!::Function, Uinf::Function, t_range;
        name="default_sim", path="./default_simulation",
        restart_path=nothing, restart_name=nothing,
        restart_step::Int=-1,
        body_solvers,
        backend=FastMultipoleBackend(;
                expansion_order=10,
                multipole_acceptance=0.4,
                leaf_size=100,
            ),
        backend_wake=backend,
        backend_solve=backend,
        backend_system=backend,
        monitors=(),
        set_Das_eta_kinematic=NaN,
        set_Das_eta_freestream=NaN,
        verbose=false,
    )
    systems_tuple = _systems_tuple(systems)
    wakes_tuple = _wakes_tuple(systems, wakes)
    _validate_body_solvers(systems, body_solvers)
    _validate_influence_backend(:backend_wake, backend_wake)
    _validate_influence_backend(:backend_system, backend_system)
    _validate_solve_backend(systems, body_solvers, backend_solve)
    audit_monitors(monitors)

    rpath = isnothing(restart_path) ? path : restart_path
    rname = isnothing(restart_name) ? name : restart_name

    # 1. Resolve restart_step from the body PVD
    body_pvd = joinpath(rpath, rname * "_body1.pvd")
    timesteps, idxs = _read_pvd_steps(body_pvd)
    if restart_step < 0
        restart_step = idxs[end]
    end
    @assert restart_step in idxs "restart_step=$(restart_step) not found in $(body_pvd) (available: $(idxs))"
    @assert restart_step + 1 < length(t_range) "restart_step=$(restart_step) leaves no steps to simulate (t_range has $(length(t_range)) entries)"

    if verbose
        println("simulate_warmstart!: resuming from step $(restart_step) (file count $(length(idxs)))")
    end

    # 2. Mirror simulate!'s pre-loop initialization on the freshly-constructed bodies.
    for sys in systems_tuple
        calc_normals!(sys)
        calc_controlpoints!(sys; off=abs(sys.CPoffset))
    end

    if !isnan(set_Das_eta_freestream) || !isnan(set_Das_eta_kinematic)
        dt0 = t_range[2] - t_range[1]
        if !isnan(set_Das_eta_freestream)
            uinf0 = Uinf(t_range[1])
            for sys in systems_tuple
                extra_reset!(sys)
                extra_apply_freestream!(sys, uinf0)
                _accumulate_Das!(sys, dt0 * set_Das_eta_freestream)
            end
        end
        if !isnan(set_Das_eta_kinematic)
            for sys in systems_tuple
                extra_reset!(sys)
            end
            kinematic_velocity!(systems_tuple, frames)
            for sys in systems_tuple
                _accumulate_Das!(sys, dt0 * set_Das_eta_kinematic)
            end
        end
        for sys in systems_tuple
            reset!(sys)
        end
    end

    # 3. Load on-disk state at restart_step.
    for (i, sys) in enumerate(systems_tuple)
        body_name = rname * "_body$(i)"
        _load_body_vtk!(sys, rpath, body_name, restart_step)
    end

    for (i, w) in enumerate(wakes_tuple)
        isnothing(w) && continue
        wake_name = rname * "_wake$(i)"
        if w isa PanelParticleWake
            _load_panel_particle_wake_vtk!(w, rpath, wake_name, restart_step)
        elseif w isa PanelWake
            _load_panel_wake_vtk!(w, rpath, wake_name, restart_step)
        else
            error("simulate_warmstart! does not yet support wake type $(typeof(w))")
        end
    end

    _load_frame_state_toml!(frames, rpath, rname, restart_step)

    # 4. Refresh derived geometry that wasn't persisted.
    for sys in systems_tuple
        calc_normals!(sys)
        calc_controlpoints!(sys; off=abs(sys.CPoffset))
    end

    # 5. Replay the end-of-step-`restart_step` actions that simulate! skipped
    #    because that step was the final step of the previous run. Use the dt
    #    that simulate! itself would use at that step.
    dt_end = t_range[restart_step + 2] - t_range[restart_step + 1]

    for w in wakes_tuple
        !isnothing(w) && propagate!(w, dt_end; step=restart_step, frames)
    end
    propagate_kinematics!(systems_tuple, frames, dt_end)
    for sys in systems_tuple
        calc_normals!(sys)
        calc_controlpoints!(sys; off=abs(sys.CPoffset))
    end
    for (sys, w) in zip(systems_tuple, wakes_tuple)
        !isnothing(w) && shed_wake!(w, sys)
    end

    # 6. Forward to simulate! with start_step pointing at the next step.
    simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
        name=name, path=path,
        body_solvers=body_solvers,
        backend=backend,
        backend_wake=backend_wake,
        backend_solve=backend_solve,
        backend_system=backend_system,
        monitors=monitors,
        set_Das_eta_kinematic=set_Das_eta_kinematic,
        set_Das_eta_freestream=set_Das_eta_freestream,
        start_step=restart_step + 1,
        verbose=verbose,
    )

    return systems, wakes, frames
end
