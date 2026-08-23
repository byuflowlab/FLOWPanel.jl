"""
    CylindricalFieldProbeMonitor

Sample induced velocity and velocity gradient on a cylindrical structured grid
laid out in the local axes of `frames[i_frame]`, evaluating one source entity
at a time so each entity's contribution can be inspected independently in
ParaView. Outputs one `.vts` file per entity per step at
`<save_path>/<name>_<entity_label>.<i_step>.vts`.

Fields written (in the global frame): `velocity` (3-vector) and
`velocity_gradient` (3×3 tensor, flattened as 9 components per point).

Entities are derived from the simulation state at call time:

- one entity per body in `systems` (label `"body<i>"`)
- for each `PanelParticleWake`: `"wake<i>_panels"` (panel/filament sources)
  and `"wake<i>_particles"` (the vortex particle field)
- for each `PanelWake`: `"wake<i>_panels"`

No "total" field is written; the per-entity files sum to the total in
ParaView via the Calculator filter.
"""
struct CylindricalFieldProbeMonitor{TF, TB} <: AbstractMonitor
    i_frame::Int
    axial_axis::Symbol
    n_axial::Int
    n_radial::Int
    n_azimuth::Int
    positions_local::Array{FastMultipole.SVector{3, TF}, 3}  # (n_axial, n_radial, n_azimuth)
    probes::FastMultipole.ProbeSystem{TF}
    backend::TB
    save_path::Union{Nothing, String}
    name::String
    verbose::Bool
end

monitor_provides(::CylindricalFieldProbeMonitor) = ()
monitor_requires(::CylindricalFieldProbeMonitor) = ()
monitor_requires_body_hessian(::CylindricalFieldProbeMonitor) = false

function CylindricalFieldProbeMonitor(;
        i_frame::Int = 1,
        axial_axis::Symbol = :x,
        axial_range = (-1.0, 1.0),
        radial_range = (0.0, 1.0),
        azimuth_range = (0.0, 2pi),
        n_axial::Int = 20,
        n_radial::Int = 8,
        n_azimuth::Int = 24,
        backend = FastMultipoleBackend(),
        save_path::Union{Nothing, AbstractString} = nothing,
        name::AbstractString = "fieldprobe",
        TF = Float64,
        verbose::Bool = false,
    )
    axial_axis in (:x, :y, :z) || throw(ArgumentError("axial_axis must be :x, :y, or :z; got $(axial_axis)"))
    n_axial   >= 2 || throw(ArgumentError("n_axial must be >= 2"))
    n_radial  >= 1 || throw(ArgumentError("n_radial must be >= 1"))
    n_azimuth >= 3 || throw(ArgumentError("n_azimuth must be >= 3"))

    z_vals = collect(range(TF(axial_range[1]),  TF(axial_range[2]);  length=n_axial))
    r_vals = collect(range(TF(radial_range[1]), TF(radial_range[2]); length=n_radial))
    θ_vals = collect(range(TF(azimuth_range[1]), TF(azimuth_range[2]); length=n_azimuth))

    positions_local = Array{FastMultipole.SVector{3, TF}, 3}(undef, n_axial, n_radial, n_azimuth)
    for k in 1:n_azimuth, j in 1:n_radial, i in 1:n_axial
        z = z_vals[i]; r = r_vals[j]; θ = θ_vals[k]
        x = r * cos(θ); y = r * sin(θ)
        p = axial_axis === :x ? (z, x, y) :
            axial_axis === :y ? (x, z, y) :
                                (x, y, z)
        positions_local[i, j, k] = FastMultipole.SVector{3, TF}(p[1], p[2], p[3])
    end

    n_total = n_axial * n_radial * n_azimuth
    probes = FastMultipole.ProbeSystem(n_total, TF)

    if save_path !== nothing
        isdir(save_path) || mkpath(save_path)
    end

    return CylindricalFieldProbeMonitor{TF, typeof(backend)}(
        i_frame, axial_axis, n_axial, n_radial, n_azimuth,
        positions_local, probes, backend,
        save_path === nothing ? nothing : String(save_path),
        String(name), verbose)
end

# ---- internals ----

_fieldprobe_entities(::Nothing, ::Int) = ()
function _fieldprobe_entities(w::PanelParticleWake, iw::Int)
    return (("wake$(iw)_panels",    get_sources(w.panel_wake)),
            ("wake$(iw)_particles", (w.pfield,)))
end
function _fieldprobe_entities(w::PanelWake, iw::Int)
    return (("wake$(iw)_panels", get_sources(w)),)
end
_fieldprobe_entities(w, iw::Int) = throw(ArgumentError(
    "CylindricalFieldProbeMonitor does not know how to split wake of type $(typeof(w))"))

function _fieldprobe_collect_entities(systems, wakes)
    entities = Tuple{String, Tuple}[]
    for (ib, body) in pairs(systems)
        push!(entities, ("body$(ib)", (body,)))
    end
    for (iw, w) in pairs(wakes)
        for e in _fieldprobe_entities(w, iw)
            push!(entities, e)
        end
    end
    return entities
end

function _fieldprobe_set_global_positions!(probes::FastMultipole.ProbeSystem{TF},
                                            positions_local::Array{<:FastMultipole.SVector{3, TF}, 3},
                                            origin_global, R_f2g) where {TF}
    n_total = length(positions_local)
    @inbounds for k in 1:n_total
        probes.position[k] = origin_global + R_f2g * positions_local[k]
    end
    return nothing
end

function _fieldprobe_zero!(probes::FastMultipole.ProbeSystem{TF}) where {TF}
    zero_v = zero(FastMultipole.SVector{3, TF})
    zero_h = zero(FastMultipole.SMatrix{3, 3, TF, 9})
    @inbounds for k in eachindex(probes.gradient)
        probes.gradient[k] = zero_v
        probes.hessian[k]  = zero_h
        probes.scalar_potential[k] = zero(TF)
    end
    return nothing
end

function _fieldprobe_write_vts(m::CylindricalFieldProbeMonitor{TF}, label::String,
                                i_step::Int) where {TF}
    m.save_path === nothing && return nothing
    na, nr, nθ = m.n_axial, m.n_radial, m.n_azimuth
    pos = m.probes.position
    grad = m.probes.gradient
    hess = m.probes.hessian

    xs = Array{TF, 3}(undef, na, nr, nθ)
    ys = Array{TF, 3}(undef, na, nr, nθ)
    zs = Array{TF, 3}(undef, na, nr, nθ)
    vel = Array{TF, 4}(undef, 3, na, nr, nθ)
    vgrad = Array{TF, 4}(undef, 9, na, nr, nθ)

    flat = 0
    @inbounds for k in 1:nθ, j in 1:nr, i in 1:na
        flat += 1
        p = pos[flat]
        xs[i, j, k] = p[1]; ys[i, j, k] = p[2]; zs[i, j, k] = p[3]
        v = grad[flat]
        vel[1, i, j, k] = v[1]; vel[2, i, j, k] = v[2]; vel[3, i, j, k] = v[3]
        H = hess[flat]
        vgrad[1, i, j, k] = H[1, 1]; vgrad[2, i, j, k] = H[2, 1]; vgrad[3, i, j, k] = H[3, 1]
        vgrad[4, i, j, k] = H[1, 2]; vgrad[5, i, j, k] = H[2, 2]; vgrad[6, i, j, k] = H[3, 2]
        vgrad[7, i, j, k] = H[1, 3]; vgrad[8, i, j, k] = H[2, 3]; vgrad[9, i, j, k] = H[3, 3]
    end

    fname = joinpath(m.save_path, "$(m.name)_$(label)_step$(lpad(i_step, 4, '0'))")
    vtk_grid(fname, xs, ys, zs) do vtk
        vtk["velocity"] = vel
        vtk["velocity_gradient"] = vgrad
    end
    return nothing
end

function (m::CylindricalFieldProbeMonitor{TF})(systems, wakes,
                                                frames::AbstractVector{<:ReferenceFrame},
                                                uinf, i_step::Int, dt::Real) where {TF}
    origin_global, R_f2g = frame_global_transform(frames, m.i_frame)
    _fieldprobe_set_global_positions!(m.probes, m.positions_local, origin_global, R_f2g)

    # Set every body's core_size to its off-body (target) value while we
    # treat panels as sources for off-body probes.
    saved_offsets = [sys.core_size for sys in systems]
    try
        for sys in systems
            sys.core_size = sys.core_size_targets
        end
        for (label, source_tuple) in _fieldprobe_collect_entities(systems, wakes)
            _fieldprobe_zero!(m.probes)
            isempty(source_tuple) && (_fieldprobe_write_vts(m, label, i_step); continue)
            influence!((m.probes,), source_tuple, m.backend;
                precalc=false,
                scalar_potential=false,
                gradient=true,
                hessian=true)
            _fieldprobe_write_vts(m, label, i_step)
            if m.verbose
                vmax = zero(TF)
                @inbounds for v in m.probes.gradient
                    s = v[1]*v[1] + v[2]*v[2] + v[3]*v[3]
                    s > vmax && (vmax = s)
                end
                println("\t\tCylindricalFieldProbeMonitor[step=$(i_step+1), $(label)]: max |velocity| = $(sqrt(vmax))")
            end
        end
    finally
        for (sys, off) in zip(systems, saved_offsets)
            sys.core_size = off
        end
    end
    return nothing
end
