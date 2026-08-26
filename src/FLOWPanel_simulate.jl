#--- Das initialization helpers ---#

function _accumulate_Das!(sys::AbstractLiftingBody, eta; min_displacement=0.0)
    for ishedding in eachindex(sys.Das)
        Das = sys.Das[ishedding]
        Vte = sys.velocity_te[ishedding]
        for j in axes(Das, 2)
            speed = sqrt(Vte[1, j]^2 + Vte[2, j]^2 + Vte[3, j]^2)
            if speed > zero(speed)
                displacement_length = max(abs(eta) * speed, min_displacement)
                signed_length = eta < zero(eta) ? -displacement_length : displacement_length
                scale = signed_length / speed
                Das[1, j] += Vte[1, j] * scale
                Das[2, j] += Vte[2, j] * scale
                Das[3, j] += Vte[3, j] * scale
            end
        end
    end
end

_accumulate_Das!(::AbstractBody, eta; min_displacement=0.0) = nothing

# Prescribed per-station |Das| (BRAINSTORM/018 chord-proportional offset):
# station_lengths[k][j] is the offset magnitude for the j-th station of the k-th
# shedding edge set; direction is the local kinematic TE tangent (velocity_te
# must already be populated by kinematic_velocity!). Accumulates like
# `_accumulate_Das!` so a freestream contribution laid down earlier survives.
function _set_Das_station_lengths!(sys::AbstractLiftingBody, station_lengths;
        min_displacement=0.0)
    length(station_lengths) == length(sys.Das) || error(
        "station_lengths has $(length(station_lengths)) shedding entries, " *
        "body has $(length(sys.Das))")
    for k in eachindex(sys.Das)
        Das = sys.Das[k]
        Vte = sys.velocity_te[k]
        lens = station_lengths[k]
        length(lens) == size(Das, 2) || error(
            "station_lengths[$(k)] has $(length(lens)) stations, " *
            "shedding $(k) has $(size(Das, 2))")
        for j in axes(Das, 2)
            speed = sqrt(Vte[1, j]^2 + Vte[2, j]^2 + Vte[3, j]^2)
            # zero TE velocity gives no direction: leave the station untouched,
            # matching `_accumulate_Das!`'s handling of stationary nodes
            speed > zero(speed) || continue
            scale = max(abs(lens[j]), min_displacement) / speed
            Das[1, j] += Vte[1, j] * scale
            Das[2, j] += Vte[2, j] * scale
            Das[3, j] += Vte[3, j] * scale
        end
    end
    return nothing
end
_set_Das_station_lengths!(::AbstractBody, station_lengths; min_displacement=0.0) = nothing

# F1b (BRAINSTORM/018 Phase 16, Route B): endpoint-on-arc variant of the
# station-lengths path. The offset endpoint is placed ON the local frozen-wake
# path instead of along the straight TE tangent: the material that left the TE
# a lag τ ago sits at the TE's backward swept-arc position
# (`_rigid_back_displacement`, summed over the frame tree exactly like
# `accumulate_Das_arc!`) PLUS the induced-drift convection u_j·τ, where u_j is
# a prescribed per-station drift velocity (global frame; typically the
# steady-state downwash measured from a settled wake). τ_j is found by
# ARC-LENGTH INTEGRATION: sub-step backward along the path, advancing the lag
# by δτ = δs/|dp/dτ| with the closed-form instantaneous backward velocity,
# until the accumulated path length reaches the prescribed |Das|_j. This is
# self-limiting near stagnation points (the endpoint stops accumulating length
# and lands short — reported, never extrapolated), needs no speed floor, and
# handles in-plane drift components exactly. The stored Das column is the
# CHORD to the endpoint (shorter than the arc length — intended; the
# admissibility bands are in arc length). Reduces to the pure
# `_rigid_back_displacement` arc as u → 0. Callers not passing drifts get the
# legacy tangent path — zero behavior change.
const _DAS_ARC_NSUB = 16                 # sub-steps per station (init-time only)
const _DAS_ARC_STAGNATION_RTOL = 1e-3    # stagnation detector (relative)

# Instantaneous backward-path velocity of one frame's contribution at lag τ:
# d/dτ of `_rigid_back_displacement(te, origin, v, ω, τ)`.
function _rigid_back_velocity(te, origin, v, ω, τ)
    d = te - origin
    ωmag = sqrt(ω[1]^2 + ω[2]^2 + ω[3]^2)
    θ = ωmag * τ
    if abs(θ) > eps(typeof(θ))^(1/3)
        axis = ω / ωmag
        drot = Rodrigues(axis, -θ) * d
        return -cross(ω, drot) - v
    else
        return -cross(ω, d) - v
    end
end

function _set_Das_station_arc_te!(body::AbstractLiftingBody, motions,
        station_lengths, station_drifts, min_displacement)
    for ishedding in eachindex(body.Das)
        Das = body.Das[ishedding]
        shedding = body.shedding[ishedding]
        lens = station_lengths[ishedding]
        drift = station_drifts[ishedding]
        length(lens) == size(Das, 2) || error(
            "station_lengths[$(ishedding)] has $(length(lens)) stations, " *
            "shedding $(ishedding) has $(size(Das, 2))")
        size(drift) == size(Das) || error(
            "station_drifts[$(ishedding)] is $(size(drift)), Das is $(size(Das))")
        for j in axes(Das, 2)
            # TE node for column j: nib for 1..nshed, nia for the last column
            # (mirrors _kinematic_velocity_te!)
            if j <= size(shedding, 2)
                node_idx = body.cells[shedding[3, j], shedding[1, j]]
            else
                node_idx = body.cells[shedding[2, end], shedding[1, end]]
            end
            te = FastMultipole.SVector{3}(body.nodes[1, node_idx],
                body.nodes[2, node_idx], body.nodes[3, node_idx])
            u = FastMultipole.SVector{3}(drift[1, j], drift[2, j], drift[3, j])
            L = max(abs(lens[j]), min_displacement)
            L > zero(L) || continue
            δs = L / _DAS_ARC_NSUB
            τ = zero(L)
            # reference speed for the stagnation detector: the total backward
            # speed at zero lag (kinematic + drift)
            V0 = u
            for m in motions
                V0 += _rigid_back_velocity(te, m.origin, m.v, m.ω, τ)
            end
            speed0 = sqrt(V0[1]^2 + V0[2]^2 + V0[3]^2)
            # stationary node with no drift: no direction, leave untouched
            # (matches the legacy paths' handling)
            speed0 > zero(speed0) || continue
            stalled = false
            for _ in 1:_DAS_ARC_NSUB
                V = u
                for m in motions
                    V += _rigid_back_velocity(te, m.origin, m.v, m.ω, τ)
                end
                speed = sqrt(V[1]^2 + V[2]^2 + V[3]^2)
                if speed < _DAS_ARC_STAGNATION_RTOL * speed0
                    stalled = true
                    break
                end
                τ += δs / speed
            end
            stalled && @warn "endpoint-on-arc Das: backward path stagnated " *
                "before reaching the prescribed length at shedding " *
                "$(ishedding), station $(j); endpoint lands short" maxlog = 8
            Δ = u * τ
            for m in motions
                Δ += _rigid_back_displacement(te, m.origin, m.v, m.ω, τ)
            end
            Das[1, j] += Δ[1]
            Das[2, j] += Δ[2]
            Das[3, j] += Δ[3]
        end
    end
    return nothing
end
_set_Das_station_arc_te!(::AbstractBody, motions, lens, drifts, mindisp) = nothing

# Collect each system's affecting frame motions (global-frame v, ω, origin)
# by walking the frame tree exactly like `accumulate_Das_arc!`.
function _collect_frame_motions!(motions, frames::AbstractVector{ReferenceFrame{TF}},
        i_frame::Int, dx_parent_to_global, R_parent_to_global) where TF
    frame = frames[i_frame]
    origin_global = R_parent_to_global * frame.x + dx_parent_to_global
    v_global = R_parent_to_global * frame.v
    ω_global = R_parent_to_global * frame.ω_axis * frame.ω
    for isurf in frame.dependent_index
        push!(motions[isurf], (v=v_global, ω=ω_global, origin=origin_global))
    end
    dx_parent_to_global = origin_global
    R_parent_to_global = R_parent_to_global * frame.R
    for i in frame.child_index
        _collect_frame_motions!(motions, frames, i, dx_parent_to_global,
            R_parent_to_global)
    end
    return nothing
end

function _set_Das_station_arc!(systems_tuple::Tuple,
        frames::AbstractVector{ReferenceFrame{TF}}, station_lengths,
        station_drifts; min_displacement=0.0) where TF
    length(station_lengths) == length(systems_tuple) || error(
        "station_lengths has $(length(station_lengths)) entries, " *
        "systems tuple has $(length(systems_tuple))")
    length(station_drifts) == length(systems_tuple) || error(
        "station_drifts has $(length(station_drifts)) entries, " *
        "systems tuple has $(length(systems_tuple))")
    motions = [NamedTuple{(:v, :ω, :origin),
        NTuple{3, FastMultipole.SVector{3,TF}}}[] for _ in systems_tuple]
    _collect_frame_motions!(motions, frames, 1,
        zero(FastMultipole.SVector{3,TF}),
        FastMultipole.SMatrix{3,3,TF,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    for (isys, sys) in enumerate(systems_tuple)
        _set_Das_station_arc_te!(sys, motions[isys], station_lengths[isys],
            station_drifts[isys], min_displacement)
    end
    return nothing
end

function _set_Das_station_lengths!(systems_tuple::Tuple, station_lengths;
        min_displacement=0.0)
    length(station_lengths) == length(systems_tuple) || error(
        "station_lengths has $(length(station_lengths)) entries, " *
        "systems tuple has $(length(systems_tuple))")
    for (sys, lens) in zip(systems_tuple, station_lengths)
        _set_Das_station_lengths!(sys, lens; min_displacement)
    end
    return nothing
end

# --- arc-following kinematic Das ---
#
# `_accumulate_Das!` lays the offset along the trailing-edge velocity, i.e. along
# the *tangent* to the path the trailing edge sweeps. That is only first-order
# accurate in θ = |ω|τ: for a rotor it lands at radius r√(1+θ²) rather than r.
# `accumulate_Das_arc!` (FLOWPanel_frames.jl) follows the actual arc instead.

_das_stash(sys::AbstractLiftingBody) = [copy(d) for d in sys.Das]
_das_stash(::AbstractBody) = nothing

function _das_zero!(sys::AbstractLiftingBody)
    for d in sys.Das
        fill!(d, zero(eltype(d)))
    end
    return nothing
end
_das_zero!(::AbstractBody) = nothing

# Apply the minimum-displacement floor to the freshly computed kinematic
# contribution, then add back whatever was already in `Das` (e.g. a freestream
# contribution accumulated earlier). This preserves `_accumulate_Das!`'s
# semantics, where the floor applies to the kinematic increment alone.
function _das_floor_restore!(sys::AbstractLiftingBody, saved, min_displacement)
    for (k, Das) in enumerate(sys.Das)
        Vte = sys.velocity_te[k]
        for j in axes(Das, 2)
            len = sqrt(Das[1,j]^2 + Das[2,j]^2 + Das[3,j]^2)
            if len < min_displacement
                if len > zero(len)
                    scale = min_displacement / len
                    for i in 1:3
                        Das[i,j] *= scale
                    end
                else
                    # degenerate (stationary node): fall back to the TE velocity
                    speed = sqrt(Vte[1,j]^2 + Vte[2,j]^2 + Vte[3,j]^2)
                    if speed > zero(speed)
                        scale = min_displacement / speed
                        for i in 1:3
                            Das[i,j] = Vte[i,j] * scale
                        end
                    end
                end
            end
            for i in 1:3
                Das[i,j] += saved[k][i,j]
            end
        end
    end
    return nothing
end
_das_floor_restore!(::AbstractBody, saved, min_displacement) = nothing

"""
    _accumulate_Das_kinematic!(systems_tuple, frames, τ; min_displacement, arc)

Accumulate the kinematic first-wake-row offset over time `τ`. With `arc=true`
(default) the offset follows the trailing edge's swept arc; with `arc=false` it
uses the legacy tangent-vector construction. The two agree to first order in
`|ω|τ` and are identical for purely translating bodies.
"""
function _accumulate_Das_kinematic!(systems_tuple::Tuple, frames, τ;
        min_displacement=0.0, arc::Bool=true)
    if !arc
        for sys in systems_tuple
            _accumulate_Das!(sys, τ; min_displacement)
        end
        return nothing
    end

    saved = map(_das_stash, systems_tuple)
    for sys in systems_tuple
        _das_zero!(sys)
    end
    accumulate_Das_arc!(systems_tuple, frames, τ)
    for (i, sys) in enumerate(systems_tuple)
        _das_floor_restore!(sys, saved[i], min_displacement)
    end
    return nothing
end

#--- wake tuple helpers ---#

function _systems_tuple(systems::Tuple)
    for system in systems
        system isa AbstractBody || throw(ArgumentError("Each system must be an AbstractBody, got $(typeof(system))"))
    end
    return systems
end
_systems_tuple(system::AbstractBody) = (system,)
_systems_tuple(systems) = throw(ArgumentError("systems must be an AbstractBody or a Tuple of AbstractBody objects, got $(typeof(systems))"))

function _wakes_tuple(systems::Tuple, wakes)
    wakes isa Tuple || throw(ArgumentError("wakes must be a Tuple when systems is a Tuple, got $(typeof(wakes))"))
    length(wakes) == length(systems) || throw(ArgumentError("Number of wakes ($(length(wakes))) must match number of systems ($(length(systems)))"))
    for wake in wakes
        (isnothing(wake) || wake isa AbstractFreeWake) || throw(ArgumentError("Each wake must be an AbstractFreeWake or nothing, got $(typeof(wake))"))
    end
    return wakes
end

function _wakes_tuple(system::AbstractBody, wake)
    (isnothing(wake) || wake isa AbstractFreeWake) || throw(ArgumentError("wake must be an AbstractFreeWake or nothing when systems is an AbstractBody, got $(typeof(wake))"))
    return (wake,)
end

_induced_vorticity_extra_outputs(targets::Tuple, enabled::Bool) =
    Tuple(enabled && target isa AbstractBody ? 3 : 0 for target in targets)

function _bound_surface_vorticity_te_info(body::AbstractBody)
    return hasproperty(body, :shedding_full) ?
        view(body.shedding_full, 1:2, :) :
        zeros(Int, 2, body.ncells)
end

function _add_bound_surface_vorticity!(systems::Tuple; grad_mu_options=(;))
    for body in systems
        _add_bound_surface_vorticity!(body; grad_mu_options)
    end
    return nothing
end

"""
    _bound_surface_vorticity!(out, body; grad_mu_options=(;))

Write the bound surface vorticity κ = n × ∇sμ of `body` into `out`
(3 × ncells, overwritten). `out` is zeroed and returned unchanged for bodies
without a doublet/Gamma strength (`has_grad_mu(body) == false`).
"""
function _bound_surface_vorticity!(out::AbstractMatrix, body::AbstractBody;
        grad_mu_options=(;))
    fill!(out, zero(eltype(out)))
    has_grad_mu(body) || return out
    normalized_grad_mu_options = _normalize_grad_mu_options(grad_mu_options;
        default_basis=:tri)

    # Build ∇sμ in a scratch buffer so compute_mu_gradient! (which writes,
    # not accumulates) does not clobber caller state, then form n × ∇sμ.
    TF = eltype(out)
    grad_mu = zeros(TF, 3, body.ncells)
    compute_mu_gradient!(grad_mu, body.controlpoints, body.normals,
        body.cells, body.neighbor,
        view(body.strength, :, get_Gammai(body)),
        _bound_surface_vorticity_te_info(body);
        scale=-1.0,
        nodes=body.nodes,
        grad_mu_options=normalized_grad_mu_options)

    @inbounds for i in axes(grad_mu, 2)
        nx, ny, nz = body.normals[1, i], body.normals[2, i], body.normals[3, i]
        gx = grad_mu[1, i]
        gy = grad_mu[2, i]
        gz = grad_mu[3, i]
        out[1, i] = ny * gz - nz * gy
        out[2, i] = nz * gx - nx * gz
        out[3, i] = nx * gy - ny * gx
    end

    return out
end

"""
    _axpy_bound_surface_vorticity!(out, body, a; grad_mu_options=(;))

Accumulate `a * κ` (κ = n × ∇sμ) into `out` (3 × ncells).
"""
function _axpy_bound_surface_vorticity!(out::AbstractMatrix, body::AbstractBody,
        a::Real; grad_mu_options=(;))
    has_grad_mu(body) || return out
    kappa = zeros(eltype(out), 3, body.ncells)
    _bound_surface_vorticity!(kappa, body; grad_mu_options)
    @inbounds out .+= a .* kappa
    return out
end

_add_bound_surface_vorticity_into!(out::AbstractMatrix, body::AbstractBody; optargs...) =
    _axpy_bound_surface_vorticity!(out, body, 1; optargs...)
_subtract_bound_surface_vorticity!(out::AbstractMatrix, body::AbstractBody; optargs...) =
    _axpy_bound_surface_vorticity!(out, body, -1; optargs...)

function _add_bound_surface_vorticity!(body::AbstractBody; grad_mu_options=(;))
    # Accumulate the bound surface vorticity κ = n × ∇sμ into
    # body.induced_vorticity on top of any wake-induced contribution already
    # there.
    _add_bound_surface_vorticity_into!(body.induced_vorticity, body; grad_mu_options)
    return nothing
end

function _validate_body_solvers(systems::Tuple, body_solvers)
    if body_solvers isa Tuple
        length(body_solvers) == length(systems) || throw(ArgumentError("Number of body_solvers ($(length(body_solvers))) must match number of systems ($(length(systems)))"))
    elseif !(body_solvers isa AbstractSolver)
        throw(ArgumentError("body_solvers must be a Tuple or an AbstractSolver, got $(typeof(body_solvers))"))
    end
    return body_solvers
end

function _validate_body_solvers(system::AbstractBody, body_solvers)
    if body_solvers isa Tuple
        throw(ArgumentError("body_solvers must be an AbstractSolver when systems is an AbstractBody, got $(typeof(body_solvers))"))
    elseif !(body_solvers isa AbstractSolver)
        throw(ArgumentError("body_solvers must be a Tuple or an AbstractSolver, got $(typeof(body_solvers))"))
    end
    return body_solvers
end

function _validate_influence_backend(name::Symbol, backend)
    backend isa AbstractBackend || throw(ArgumentError("$(name) must be an AbstractBackend for influence! calls, got $(typeof(backend))"))
    return backend
end

function _validate_solve_backend(systems, body_solvers, backend)
    if systems isa AbstractBody && (backend isa Tuple || backend isa AbstractVector)
        throw(ArgumentError("backend_solve must be a single backend when systems is an AbstractBody, got $(typeof(backend))"))
    elseif body_solvers isa Tuple && (backend isa Tuple || backend isa AbstractVector)
        length(backend) == length(_systems_tuple(systems)) || throw(ArgumentError("backend_solve length ($(length(backend))) must match number of systems ($(length(_systems_tuple(systems))))"))
    elseif body_solvers isa AbstractSolver && (backend isa Tuple || backend isa AbstractVector)
        throw(ArgumentError("backend_solve must be a single backend when body_solvers is an AbstractSolver, got $(typeof(backend))"))
    end
    return backend
end

function _collect_wake_probes(wakes::Tuple)
    result = ()
    for w in wakes
        if !isnothing(w)
            result = (result..., get_probes(w)...)
        end
    end
    return result
end

_collect_wake_probes(wake::AbstractFreeWake) = _collect_wake_probes((wake,))
_collect_wake_probes(::Nothing) = ()

function _set_core_sizes!(systems::Tuple, field::Symbol)
    for system in systems
        system.core_size = getfield(system, field)
    end
    return nothing
end

function _self_panel_core_size_conditioning()
    before! = function (source_buffer, source_system, i_source_system, target_buffer, i_target_system)
        source_system.core_size = source_system.core_size_panel
        return nothing
    end
    after! = function (source_buffer, source_system, i_source_system, target_buffer, i_target_system)
        source_system.core_size = source_system.core_size_targets
        return nothing
    end
    return FastMultipole.DirectConditioningRule(FastMultipole.SelfPairs(), before!, after!)
end

function _collect_wake_sources(wakes::Tuple)
    result = ()
    for w in wakes
        if !isnothing(w)
            result = (result..., get_sources(w)...)
        end
    end
    return result
end

_collect_wake_sources(wake::AbstractFreeWake) = _collect_wake_sources((wake,))
_collect_wake_sources(::Nothing) = ()

function _particle_j_rms(pfield::FLOWVPM.ParticleField)
    np = pfield.np
    np == 0 && return 0.0
    J = view(pfield.particles, FLOWVPM.J_INDEX, 1:np)
    return sqrt(sum(abs2, J) / np)
end

function _set_particle_j_zero!(pfield::FLOWVPM.ParticleField)
    pfield.np == 0 && return nothing
    view(pfield.particles, FLOWVPM.J_INDEX, 1:pfield.np) .= 0
    return nothing
end

function _diagnostic_particle_j_from_sources!(pfield::FLOWVPM.ParticleField, sources::Tuple,
        backend; velocity_gradient::Bool=true, optargs...)
    length(sources) == 0 && return 0.0
    pfield.np == 0 && return 0.0

    particles = copy(view(pfield.particles, :, 1:pfield.np))
    try
        _set_particle_j_zero!(pfield)
        influence!((pfield,), sources, backend;
            scalar_potential=false,
            velocity=false,
            velocity_gradient=(velocity_gradient,),
            precalc=false,
            postcalc=false,
            optargs...)
        return _particle_j_rms(pfield)
    finally
        view(pfield.particles, :, 1:size(particles, 2)) .= particles
    end
end

function _diagnose_particle_influence!(wakes_tuple::Tuple, systems_tuple::Tuple,
        backend_wake, backend_system; needs_induced_vorticity::Bool=false,
        particle_hessian_self::Bool=true, diagnostic_vertical=(0.0, 0.0, 1.0))
    for (iw, w) in enumerate(wakes_tuple)
        w isa PanelParticleWake || continue
        pfield = w.pfield
        pfield.np == 0 && continue

        wake_sources = _collect_wake_sources((w,))
        panel_sources = Tuple(s for s in wake_sources if !(s isa FLOWVPM.ParticleField))
        particle_sources = Tuple(s for s in wake_sources if s isa FLOWVPM.ParticleField)

        j_total = _particle_j_rms(pfield)
        j_body = _diagnostic_particle_j_from_sources!(pfield, systems_tuple, backend_system;
            direct_conditioning=_self_panel_core_size_conditioning())
        j_panelwake = _diagnostic_particle_j_from_sources!(pfield, panel_sources, backend_wake)
        j_particles = _diagnostic_particle_j_from_sources!(pfield, particle_sources, backend_wake;
            velocity_gradient=particle_hessian_self)

        println("particle influence wake=$(iw) np=$(pfield.np) |J_total|=$(j_total) |J_body|=$(j_body) |J_panelwake|=$(j_panelwake) |J_particles|=$(j_particles)")
    end
    return nothing
end

# Wake sources that have a well-defined scalar potential. Excludes vortex
# particle fields (which carry only a vector potential).
_scalar_potential_sources(w::PanelWake) = get_sources(w)
_scalar_potential_sources(w::PanelParticleWake) = get_sources(w.panel_wake)
_scalar_potential_sources(::Nothing) = ()

function _collect_wake_scalar_sources(wakes::Tuple)
    result = ()
    for w in wakes
        result = (result..., _scalar_potential_sources(w)...)
    end
    return result
end

function initialize_Das!(systems, frames, Uinf::Function, t0, dt0;
        set_Das_eta_kinematic=NaN,
        set_Das_eta_freestream=NaN,
        set_Das_min_kinematic_displacement=0.0,
        set_Das_kinematic_arc::Bool=true,
        set_Das_station_lengths=nothing,
        set_Das_station_drifts=nothing)
    if isnan(set_Das_eta_freestream) && isnan(set_Das_eta_kinematic) &&
            isnothing(set_Das_station_lengths)
        isnothing(set_Das_station_drifts) || error(
            "set_Das_station_drifts requires set_Das_station_lengths " *
            "(the endpoint-on-arc path needs the per-station lengths)")
        return systems
    end

    systems_tuple = _systems_tuple(systems)

    if !isnan(set_Das_eta_freestream)
        uinf0 = Uinf(t0)
        for sys in systems_tuple
            extra_reset!(sys)
            extra_apply_freestream!(sys, uinf0)
            _accumulate_Das!(sys, dt0 * set_Das_eta_freestream)
        end
    end

    if !isnan(set_Das_eta_kinematic)
        isnothing(set_Das_station_lengths) || error(
            "set_Das_eta_kinematic and set_Das_station_lengths are mutually " *
            "exclusive: both set the kinematic first-row offset")
        for sys in systems_tuple
            extra_reset!(sys)
        end
        kinematic_velocity!(systems_tuple, frames)
        _accumulate_Das_kinematic!(systems_tuple, frames, dt0 * set_Das_eta_kinematic;
            min_displacement=set_Das_min_kinematic_displacement,
            arc=set_Das_kinematic_arc)
    end

    if !isnothing(set_Das_station_lengths)
        # Prescribed per-station |Das| (e.g. a fraction of the local chord,
        # BRAINSTORM/018): direction is the local kinematic TE tangent, exactly
        # as the eta path, but the magnitude at station j is the given length
        # instead of eta*dt*|V_te,j|. dt-independent by construction, so it is
        # inherently safe inside dt-refinement studies.
        # With `set_Das_station_drifts` (F1b, Phase 16 Route B) the endpoint is
        # instead placed ON the frozen-wake path (backward swept arc + induced
        # drift) at the prescribed arc length — see _set_Das_station_arc!.
        for sys in systems_tuple
            extra_reset!(sys)
        end
        kinematic_velocity!(systems_tuple, frames)
        if isnothing(set_Das_station_drifts)
            _set_Das_station_lengths!(systems_tuple, set_Das_station_lengths;
                min_displacement=set_Das_min_kinematic_displacement)
        else
            _set_Das_station_arc!(systems_tuple, frames,
                set_Das_station_lengths, set_Das_station_drifts;
                min_displacement=set_Das_min_kinematic_displacement)
        end
    end

    # reset velocity fields modified during Das computation
    for sys in systems_tuple
        reset!(sys)
    end

    return systems
end

# ------------------------------------------------------------------------
# Stage helpers of _steady_aerodynamics! (pure code motion; BRAINSTORM 015
# Phase 3). The legacy function calls them in the identical statement order;
# the non-default Kutta runtime (_kutta_step! in FLOWPanel_kutta.jl) reuses
# them so both paths share one implementation of each stage.
# ------------------------------------------------------------------------

"Collect wake probes/targets/sources (pure; stage helper of
`_steady_aerodynamics!`)."
function _sa_collect(systems_tuple::Tuple, wakes_tuple::Tuple)
    wake_probes = _collect_wake_probes(wakes_tuple)
    targets = (systems_tuple..., wake_probes...)
    wake_sources = _collect_wake_sources(wakes_tuple)
    return wake_probes, targets, wake_sources
end

"Reset wakes and bodies, apply the freestream, and add kinematic velocities
(stage helper of `_steady_aerodynamics!`)."
function _sa_reset_freestream_kinematic!(systems_tuple::Tuple,
        wakes_tuple::Tuple, frames, uinf)
    for w in wakes_tuple
        !isnothing(w) && reset!(w)
    end
    for sys in systems_tuple
        reset!(sys)
    end

    apply_freestream!(systems_tuple, uinf)
    for w in wakes_tuple
        !isnothing(w) && apply_freestream!(w, uinf)
    end

    kinematic_velocity!(systems_tuple, frames)
    return nothing
end

function _steady_aerodynamics!(systems, systems_tuple::Tuple, wakes_tuple::Tuple,
        frames, uinf, body_solvers; backend_wake=nothing, backend_solve,
        backend_system, needs_induced_vorticity::Bool=false,
        update_trailing_edges::Bool=false,
        wakerow_no_hessian_to_particles::Bool=false,
        body_hessian_to_particles::Bool=false,
        body_gradient_core_size::Float64=NaN,
        body_on_wake::Bool=true,
        panel_wake_on_particles::Bool=true,
        particle_hessian_self::Bool=true,
        diagnose_particle_influence::Bool=false,
        diagnostic_vertical=(0.0, 0.0, 1.0),
        grad_mu_options=(;),
        formulation::AbstractSolveFormulation=VelocityThroughSources(),
        formulation_state=nothing,
        i_step::Int=0)
    normalized_grad_mu_options = _normalize_grad_mu_options(grad_mu_options;
        default_basis=:quad)

    wake_probes, targets, wake_sources = _step_timer_measure(:remaining_aerodynamics) do
        collected = _sa_collect(systems_tuple, wakes_tuple)
        _sa_reset_freestream_kinematic!(systems_tuple, wakes_tuple, frames, uinf)

        if update_trailing_edges
            for (sys, w) in zip(systems_tuple, wakes_tuple)
                !isnothing(w) && update_TE!(w, sys)
            end
        end

        # snapshot pre-wake control-point velocity for formulations that isolate
        # the wake-only contribution afterwards (no-op for the default)
        formulation_prewake!(formulation, formulation_state, systems_tuple)
        collected
    end

    _step_timer_measure(:wake_influence) do
        _sa_wake_influence!(targets, wake_sources, backend_wake;
            needs_induced_vorticity, wakerow_no_hessian_to_particles,
            panel_wake_on_particles, particle_hessian_self)
    end

    _step_timer_measure(:solve) do
        _set_core_sizes!(systems_tuple, :core_size_panel)
        solve_formulation!(formulation, formulation_state, systems, systems_tuple,
            wakes_tuple, body_solvers; backend_solve, backend_wake, i_step)
    end

    _step_timer_measure(:remaining_aerodynamics) do
        needs_induced_vorticity && _add_bound_surface_vorticity!(systems_tuple;
            grad_mu_options=normalized_grad_mu_options)
        _set_core_sizes!(systems_tuple, :core_size_targets)
    end

    _step_timer_measure(:body_influence) do
        _sa_body_influence!(targets, systems_tuple, backend_system;
            needs_induced_vorticity, body_on_wake, body_hessian_to_particles,
            body_gradient_core_size)
    end

    _step_timer_measure(:remaining_aerodynamics) do
        if diagnose_particle_influence
            _diagnose_particle_influence!(wakes_tuple, systems_tuple, backend_wake, backend_system;
                needs_induced_vorticity, particle_hessian_self, diagnostic_vertical)
        end

        _sa_half_jump!(systems_tuple, normalized_grad_mu_options)
    end

    return nothing
end

"Frozen wake-influence stage of `_steady_aerodynamics!` (pure code motion):
wake sources → bodies and wake probes, with the legacy diagnostic gates."
function _sa_wake_influence!(targets::Tuple, wake_sources::Tuple, backend_wake;
        needs_induced_vorticity::Bool=false,
        wakerow_no_hessian_to_particles::Bool=false,
        panel_wake_on_particles::Bool=true,
        particle_hessian_self::Bool=true)
    if length(wake_sources) > 0
        # Diagnostic gates:
        #   wakerow_no_hessian_to_particles: ablate the panel-wake-row ->
        #     particle Hessian (vortex rings + closing filaments produce
        #     |∇U| ~ 1/r^2 at edges that linear panels would smooth into a sheet).
        #   panel_wake_on_particles=false: drop panel-wake-row -> particle
        #     velocity entirely; particles still receive particle-on-particle
        #     velocity, and panel wake still acts on bodies and panel-wake nodes.
        if wakerow_no_hessian_to_particles || !panel_wake_on_particles
            panel_sources  = Tuple(s for s in wake_sources if !(s isa FLOWVPM.ParticleField))
            pfield_sources = Tuple(s for s in wake_sources if   s isa FLOWVPM.ParticleField)
            panel_targets = !panel_wake_on_particles ?
                Tuple(t for t in targets if !(t isa FLOWVPM.ParticleField)) : targets
            if length(panel_sources) > 0 && length(panel_targets) > 0
                influence!(panel_targets, panel_sources, backend_wake; precalc=true,
                    scalar_potential=false,
                    velocity=true,
                    velocity_gradient=Tuple(sys isa FLOWVPM.ParticleField ? false : requires_hessian(sys) for sys in panel_targets),
                    extra_outputs=_induced_vorticity_extra_outputs(panel_targets, needs_induced_vorticity))
            end
            if length(pfield_sources) > 0
                particle_targets = Tuple(t for t in targets if t isa FLOWVPM.ParticleField)
                nonparticle_targets = Tuple(t for t in targets if !(t isa FLOWVPM.ParticleField))
                if particle_hessian_self || length(particle_targets) == 0
                    influence!(targets, pfield_sources, backend_wake; precalc=true,
                        postcalc=true,
                        scalar_potential=false,
                        velocity=true,
                        velocity_gradient=Tuple(requires_hessian(sys) for sys in targets),
                        extra_outputs=_induced_vorticity_extra_outputs(targets, needs_induced_vorticity))
                else
                    if length(nonparticle_targets) > 0
                        influence!(nonparticle_targets, pfield_sources, backend_wake; precalc=true,
                            postcalc=false,
                            scalar_potential=false,
                            velocity=true,
                            velocity_gradient=Tuple(requires_hessian(sys) for sys in nonparticle_targets),
                            extra_outputs=_induced_vorticity_extra_outputs(nonparticle_targets, needs_induced_vorticity))
                    end
                    influence!(particle_targets, pfield_sources, backend_wake; precalc=true,
                        postcalc=true,
                        scalar_potential=false,
                        velocity=true,
                        velocity_gradient=Tuple(false for _ in particle_targets),
                        extra_outputs=_induced_vorticity_extra_outputs(particle_targets, needs_induced_vorticity))
                end
            end
        elseif !particle_hessian_self && any(source isa FLOWVPM.ParticleField for source in wake_sources)
            panel_sources = Tuple(s for s in wake_sources if !(s isa FLOWVPM.ParticleField))
            pfield_sources = Tuple(s for s in wake_sources if s isa FLOWVPM.ParticleField)
            if length(panel_sources) > 0
                influence!(targets, panel_sources, backend_wake; precalc=true,
                    scalar_potential=false,
                    velocity=true,
                    velocity_gradient=Tuple(requires_hessian(sys) for sys in targets),
                    extra_outputs=_induced_vorticity_extra_outputs(targets, needs_induced_vorticity))
            end
            particle_targets = Tuple(t for t in targets if t isa FLOWVPM.ParticleField)
            nonparticle_targets = Tuple(t for t in targets if !(t isa FLOWVPM.ParticleField))
            if length(pfield_sources) > 0
                if length(nonparticle_targets) > 0
                    influence!(nonparticle_targets, pfield_sources, backend_wake; precalc=true,
                        postcalc=false,
                        scalar_potential=false,
                        velocity=true,
                        velocity_gradient=Tuple(requires_hessian(sys) for sys in nonparticle_targets),
                        extra_outputs=_induced_vorticity_extra_outputs(nonparticle_targets, needs_induced_vorticity))
                end
                if length(particle_targets) > 0
                    influence!(particle_targets, pfield_sources, backend_wake; precalc=true,
                        postcalc=true,
                        scalar_potential=false,
                        velocity=true,
                        velocity_gradient=Tuple(false for _ in particle_targets),
                        extra_outputs=_induced_vorticity_extra_outputs(particle_targets, needs_induced_vorticity))
                end
            end
        else
            wake_postcalc = any(source isa FLOWVPM.ParticleField for source in wake_sources)
            influence!(targets, wake_sources, backend_wake; precalc=true,
                postcalc=wake_postcalc,
                scalar_potential=false,
                velocity=true,
                velocity_gradient=Tuple(requires_hessian(sys) for sys in targets),
                extra_outputs=_induced_vorticity_extra_outputs(targets, needs_induced_vorticity))
        end
    end
    return nothing
end

"Post-solve body-influence stage of `_steady_aerodynamics!` (pure code
motion): body → (bodies, wake probes) at `core_size_targets`, with the
legacy `body_on_wake` and split-gradient gates. The caller sets the target
kernel offsets first."
function _sa_body_influence!(targets::Tuple, systems_tuple::Tuple,
        backend_system;
        needs_induced_vorticity::Bool=false,
        body_on_wake::Bool=true,
        body_hessian_to_particles::Bool=false,
        body_gradient_core_size::Float64=NaN)
    if !body_on_wake
        # body-on-body only; skip the body-on-wake-probes pass so the wake
        # never receives body-induced velocity this step.
        influence!(systems_tuple, systems_tuple, backend_system; precalc=false,
            scalar_potential=false,
            velocity=true,
            velocity_gradient=Tuple(requires_hessian(sys) for sys in systems_tuple),
            extra_outputs=_induced_vorticity_extra_outputs(systems_tuple, needs_induced_vorticity),
            direct_conditioning=_self_panel_core_size_conditioning())
    else
        # Optional split regularization: evaluate the body->particle velocity
        # GRADIENT with a larger kernel offset than the velocity itself
        # (body_gradient_core_size; NaN = disabled). Not strictly physical —
        # it smooths the |∇U| "bumpiness" of piecewise-constant doublet panels
        # felt by nearby particles, while leaving the advecting velocity at the
        # sharper core_size_targets.
        split_grad = body_hessian_to_particles && !isnan(body_gradient_core_size) &&
            any(t isa FLOWVPM.ParticleField for t in targets)
        if !split_grad
            influence!(targets, systems_tuple, backend_system; precalc=false,
                scalar_potential=false,
                velocity=true,
                velocity_gradient=Tuple((sys isa FLOWVPM.ParticleField && !body_hessian_to_particles) ? false : requires_hessian(sys) for sys in targets),
                extra_outputs=_induced_vorticity_extra_outputs(targets, needs_induced_vorticity),
                direct_conditioning=_self_panel_core_size_conditioning())
        else
            # pass 1: velocity for all targets (+ gradient for non-particle
            # targets) at core_size_targets
            influence!(targets, systems_tuple, backend_system; precalc=false,
                scalar_potential=false,
                velocity=true,
                velocity_gradient=Tuple(sys isa FLOWVPM.ParticleField ? false : requires_hessian(sys) for sys in targets),
                extra_outputs=_induced_vorticity_extra_outputs(targets, needs_induced_vorticity),
                direct_conditioning=_self_panel_core_size_conditioning())
            # pass 2: gradient for particle targets at the larger offset. The
            # velocity computed alongside (backends may not support
            # hessian-only) is discarded via snapshot/restore so particles keep
            # the pass-1 velocity. No body self-pairs here, so the self-panel
            # conditioning rule is unnecessary (its after! hook would also
            # clobber the gradient offset).
            particle_targets = Tuple(t for t in targets if t isa FLOWVPM.ParticleField)
            for sys in systems_tuple
                sys.core_size = body_gradient_core_size
            end
            saved_U = [copy(pf.particles[FLOWVPM.U_INDEX, 1:pf.np]) for pf in particle_targets]
            influence!(particle_targets, systems_tuple, backend_system; precalc=false,
                scalar_potential=false,
                velocity=true,
                velocity_gradient=Tuple(true for _ in particle_targets))
            for (pf, U0) in zip(particle_targets, saved_U)
                pf.particles[FLOWVPM.U_INDEX, 1:pf.np] .= U0
            end
            _set_core_sizes!(systems_tuple, :core_size_targets)
        end
    end
    return nothing
end

"Exterior half-jump stage of `_steady_aerodynamics!` (pure code motion): add
the +½∇μ tangential half-jump on each surface so body.velocity is the
EXTERIOR surface limit (matching OLD calcfield_U!). The kernel-induced
velocity at the on-surface centroid is the PV (continuous through the doublet
sheet); the exterior limit requires this extra tangential half-jump that
depends on neighbor strengths, which _self_limit cannot supply locally."
function _sa_half_jump!(systems_tuple::Tuple, normalized_grad_mu_options)
    for body in systems_tuple
        if has_grad_mu(body)
            compute_mu_gradient!(body.velocity, body.controlpoints, body.normals,
                body.cells, body.neighbor,
                view(body.strength, :, get_Gammai(body)),
                _bound_surface_vorticity_te_info(body);
                scale=0.5,
                nodes=body.nodes,
                grad_mu_options=normalized_grad_mu_options)
        end
    end
    return nothing
end

function _clean_monitor_csv_files!(path, name)
    isnothing(path) && return nothing
    dir = joinpath(path, "monitors")
    isdir(dir) || return nothing
    prefix = "$(name)_monitor"
    for file in readdir(dir; join=false)
        if startswith(file, prefix) && endswith(file, ".csv")
            rm(joinpath(dir, file); force=true)
        end
    end
    return nothing
end

"""
    steady!(systems, frames, uinf; body_solvers, backend=FastMultipoleBackend(...),
            backend_solve=backend, backend_system=backend, monitors=(), i_run=1,
            dt=1.0, path=nothing, name="default_steady", verbose=false)

Run a one-shot steady panel solve for one body or a tuple of bodies. The solve
uses wake geometry already stored on the bodies, if any, and does not create,
propagate, or shed external free wakes.
"""
function steady!(systems, frames, uinf;
        name="default_steady",
        path=nothing,
        body_solvers,
        backend=FastMultipoleBackend(;
                expansion_order=10,
                multipole_acceptance=0.4,
                leaf_size=100,
            ),
        backend_solve=backend,
        backend_system=backend,
        monitors=(),
        i_run::Int=1,
        dt::Real=1.0,
        clean_files::Bool=true,
        compress_vtk::Bool=true,
        wakerow_no_hessian_to_particles::Bool=false,
        body_hessian_to_particles::Bool=false,
        body_gradient_core_size::Float64=NaN,
        body_on_wake::Bool=true,
        panel_wake_on_particles::Bool=true,
        particle_hessian_self::Bool=true,
        diagnose_particle_influence::Bool=false,
        diagnostic_vertical=(0.0, 0.0, 1.0),
        grad_mu_options=(;),
        wake_attachment::AbstractWakeAttachment=RigidTransitionAttachment(),
        kutta_closure::AbstractKuttaClosure=JumpKutta(),
        verbose=false
    )
    i_run >= 1 || throw(ArgumentError("i_run must be >= 1, got $(i_run)."))
    dt > 0 || throw(ArgumentError("dt must be positive, got $(dt)."))

    systems_tuple = _systems_tuple(systems)
    wakes_tuple = Tuple(nothing for _ in systems_tuple)
    _validate_body_solvers(systems, body_solvers)
    _validate_influence_backend(:backend_system, backend_system)
    _validate_solve_backend(systems, body_solvers, backend_solve)
    audit_monitors(monitors)

    # Non-default wake-attachment/Kutta-closure configuration (BRAINSTORM 015):
    # steady! supports Route A only, validated before any state is mutated.
    kutta_is_legacy = _is_legacy_kutta(wake_attachment, kutta_closure)
    if !kutta_is_legacy
        _validate_kutta_configuration(:steady, systems_tuple, wakes_tuple,
            body_solvers, VelocityThroughSources(), backend_system,
            wake_attachment, kutta_closure)
    end

    i_step = i_run - 1
    verbose && println("\tsteady run $(i_run)")

    needs_grad = any(monitor_requires_body_hessian, monitors)
    needs_induced_vorticity = any(monitor_requires_induced_vorticity, monitors)
    for sys in systems_tuple
        sys.needs_velocity_gradient[] = needs_grad
    end

    if !isnothing(path) && !isdir(path)
        mkpath(path)
    end
    clean_files && i_run == 1 && _clean_monitor_csv_files!(path, name)

    for sys in systems_tuple
        calc_normals!(sys)
        calc_controlpoints!(sys)
    end

    if kutta_is_legacy
        _steady_aerodynamics!(systems, systems_tuple, wakes_tuple, frames, uinf,
            body_solvers; backend_solve, backend_system, needs_induced_vorticity,
            wakerow_no_hessian_to_particles,
            body_hessian_to_particles,
            body_gradient_core_size,
            body_on_wake,
            panel_wake_on_particles,
            particle_hessian_self,
            diagnose_particle_influence,
            diagnostic_vertical,
            grad_mu_options)
    else
        kutta_runtime = _initialize_kutta(:steady, systems_tuple[1],
            _single_body_solver(body_solvers), nothing,
            wake_attachment, kutta_closure)
        _kutta_step!(kutta_runtime, systems_tuple, wakes_tuple, frames, uinf;
            backend_solve, backend_system, needs_induced_vorticity,
            grad_mu_options, i_step,
            wakerow_no_hessian_to_particles,
            body_hessian_to_particles,
            body_gradient_core_size,
            body_on_wake,
            panel_wake_on_particles,
            particle_hessian_self)
    end

    monitor_context = MonitorContext()
    monitor_set_time!(monitor_context, i_step * dt)
    monitor_csv_dir = isnothing(path) ? nothing : joinpath(path, "monitors")
    for (i_monitor, monitor) in enumerate(monitors)
        _run_monitor!(monitor, monitor_context, systems_tuple, wakes_tuple, frames, uinf, i_step, dt)
        if !isnothing(monitor_csv_dir)
            write_monitor_csv!(monitor, monitor_csv_dir, name, i_monitor,
                               monitor_context, systems_tuple, i_step, dt;
                               overwrite=i_run == 1)
        end
    end

    if !isnothing(path)
        t = i_step * dt
        for (i, sys) in enumerate(systems_tuple)
            body_name = name * "_body$(i)"
            write_vtk(joinpath(path, body_name), sys, i_step, t;
                      monitors=monitors, i_system=i, overwrite=i_run==1,
                      compress=compress_vtk)
        end
    end

    return systems
end

"""
    simulate!(systems, wakes, frames, maneuver!, Uinf, t_range; body_solvers, backend=FastMultipoleBackend(...), backend_wake=backend, backend_solve=backend, backend_system=backend, monitors=(), ...)

Advance one or more coupled body-wake systems through `t_range`, solving the
aerodynamics, updating wakes, optionally writing VTK output, and calling any
registered monitors.
"""
function simulate!(systems, wakes, frames, maneuver!::Function, Uinf::Function, t_range;
        name="default_sim", path=joinpath("data", "default_simulation"),
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
        set_Das_min_kinematic_displacement=0.0,
        set_Das_kinematic_arc::Bool=true,
        set_Das_refresh::Bool=false,
        start_step::Int=0,
        clean_files::Bool=true,
        compress_vtk::Bool=true,
        wakerow_no_hessian_to_particles::Bool=false,
        body_hessian_to_particles::Bool=false,
        body_gradient_core_size::Float64=NaN,
        body_on_wake::Bool=true,
        panel_wake_on_particles::Bool=true,
        particle_hessian_self::Bool=true,
        particle_relax::Bool=true,
        bound_strength_rlx::Real=1.0,
        diagnose_particle_gamma::Bool=false,
        diagnose_particle_influence::Bool=false,
        diagnostic_vertical=(0.0, 0.0, 1.0),
        grad_mu_options=(;),
        formulation::AbstractSolveFormulation=VelocityThroughSources(),
        wake_attachment::AbstractWakeAttachment=RigidTransitionAttachment(),
        kutta_closure::AbstractKuttaClosure=JumpKutta(),
        verbose=false
    )
    @assert 0 <= start_step < length(t_range) "start_step ($(start_step)) must be in [0, $(length(t_range))-1)"
    systems_tuple = _systems_tuple(systems)
    wakes_tuple = _wakes_tuple(systems, wakes)
    _validate_body_solvers(systems, body_solvers)
    _validate_influence_backend(:backend_wake, backend_wake)
    _validate_influence_backend(:backend_system, backend_system)
    _validate_solve_backend(systems, body_solvers, backend_solve)
    audit_monitors(monitors)

    # Non-default wake-attachment/Kutta-closure configurations (BRAINSTORM 015)
    # are validated in full before any body, wake, or solver state is mutated.
    # The exact legacy pair skips this entirely and branches into the
    # pre-existing call sequence with no new allocation.
    kutta_is_legacy = _is_legacy_kutta(wake_attachment, kutta_closure)
    if !kutta_is_legacy
        _validate_kutta_configuration(:simulate, systems_tuple, wakes_tuple,
            body_solvers, formulation, backend_system, wake_attachment,
            kutta_closure;
            bound_strength_rlx, set_Das_eta_kinematic, set_Das_eta_freestream,
            set_Das_min_kinematic_displacement, set_Das_kinematic_arc,
            set_Das_refresh)
    end

    # Flip body.needs_velocity_gradient based on monitor contracts so that
    # requires_hessian(body) propagates the right HS flag into the per-step
    # influence! calls. Done once before the time loop; bodies are mutable
    # so the Ref persists across steps.
    needs_grad = any(monitor_requires_body_hessian, monitors)
    needs_induced_vorticity = any(monitor_requires_induced_vorticity, monitors)
    for sys in systems_tuple
        sys.needs_velocity_gradient[] = needs_grad
    end

    # create save path if it does not exist
    if !isnothing(path) && !isdir(path)
        mkpath(path)
    end
    clean_files && start_step == 0 && _clean_monitor_csv_files!(path, name)

    # update control points and normals according to Neumann/Dirichlet BCs
    for sys in systems_tuple
        calc_normals!(sys)
        calc_controlpoints!(sys)
    end

    if !isnan(set_Das_eta_freestream) || !isnan(set_Das_eta_kinematic)
        initialize_Das!(systems_tuple, frames, Uinf, t_range[1], t_range[2] - t_range[1];
            set_Das_eta_kinematic, set_Das_eta_freestream,
            set_Das_min_kinematic_displacement, set_Das_kinematic_arc)
    end

    # Validate the wake→body solve formulation and build its runtime state
    # (one-time operator assembly; `nothing` for the default formulation).
    # Placed after geometry/Das initialization so cached operators see final
    # trailing-edge wake directions.
    formulation_state = initialize_formulation(formulation, systems_tuple,
        wakes_tuple, body_solvers, backend_solve, backend_system)

    # Build the non-default Kutta runtime after geometry/Das initialization so
    # Route A's one-time attachment operator sees the final trailing-edge wake
    # directions (nothing on the legacy default path).
    kutta_runtime = kutta_is_legacy ? nothing :
        _initialize_kutta(:simulate, systems_tuple[1],
            _single_body_solver(body_solvers), wakes_tuple[1],
            wake_attachment, kutta_closure)

    # Body bound-circulation low-pass (item 005 E4.8): when bound_strength_rlx < 1
    # we blend each step's freshly solved strength with the previous (relaxed)
    # strength, Γ_n = (1-α)Γ_{n-1} + α Γ_solve, to artificially damp the body↔wake
    # feedback loop. α = 1 (default) is a no-op. State is one buffer per body.
    apply_bound_rlx = bound_strength_rlx != 1
    prev_strengths = apply_bound_rlx ? [copy(sys.strength) for sys in systems_tuple] : nothing
    have_prev_strength = false
    # The TraceCorrected Kutta-trace correction c is blended with the same α so
    # the affine wake strength γ = C·μ̃ − c stays consistent with the relaxed μ̃.
    prev_c = apply_bound_rlx && formulation isa TraceCorrected ?
        copy(formulation_state.c) : nothing

    # begin simulation
    i_step = start_step
    _t_wall_prev = NaN
    for t in @view t_range[start_step+1:end]
        _step_timer_token = _step_timer_begin_step!()
        # When full step timing is explicitly armed, make WakeHealthMonitor's
        # first wall sample cover the first continuation step instead of its
        # historical NaN sentinel. Off-state behavior remains unchanged.
        if _step_timer_token !== nothing
            for monitor in monitors
                if monitor isa WakeHealthMonitor && isnan(monitor.t_last)
                    monitor.t_last = time()
                end
            end
        end
        if verbose
            # Wall-clock stamp and flush: stdout is block-buffered when
            # redirected to a Slurm log, so an unflushed per-step line leaves a
            # long run looking hung and gives no timing at all. The only timing
            # otherwise available is a single `@time` around the whole call,
            # which is why early-vs-late per-step cost has never been measured
            # (BRAINSTORM/018 S0a). WakeHealthMonitor records `wall_s` in CSV;
            # this is the same information for runs without that monitor.
            _t_wall_now = time()
            _dt_wall = isnan(_t_wall_prev) ? NaN : _t_wall_now - _t_wall_prev
            _t_wall_prev = _t_wall_now
            println("\tstep $(i_step)/$(length(t_range)-1) at time $(t)" *
                    (isnan(_dt_wall) ? "" : "  [$(round(_dt_wall, sigdigits=3)) s]"))
            flush(stdout)
        end

        dynamics_toggle, uinf, dt = _step_timer_measure(:controls_setup) do
            #------- controls -------#

            # update frames based on maneuver
            # (RPMs, tilting systems, prescribed trajectory, etc.)
            dynamics_toggle_local = maneuver!(frames, systems_tuple, wakes_tuple, t)

            #------- aerodynamics -------#

            uinf_local = Uinf(t)

            # dt for this step
            dt_local = i_step < length(t_range) - 1 ? t_range[i_step+2] - t_range[i_step+1] : t_range[i_step+1] - t_range[i_step]

            # Re-derive the first-wake-row offset from the *current* kinematic state
            # (BRAINSTORM 014). By default Das is frozen at its t=0 magnitude and
            # only rotated by propagate_kinematics!, so it does not track the
            # operating condition (e.g. it retains a spin-up-fraction RPM). Zeroing
            # first is mandatory: the accumulate helpers are `+=`-based.
            if set_Das_refresh &&
                    (!isnan(set_Das_eta_freestream) || !isnan(set_Das_eta_kinematic))
                for sys in systems_tuple
                    _das_zero!(sys)
                end
                initialize_Das!(systems_tuple, frames, Uinf, t, dt_local;
                    set_Das_eta_kinematic, set_Das_eta_freestream,
                    set_Das_min_kinematic_displacement, set_Das_kinematic_arc)
            end
            (dynamics_toggle_local, uinf_local, dt_local)
        end

        if isnothing(kutta_runtime)
            _steady_aerodynamics!(systems, systems_tuple, wakes_tuple, frames, uinf,
                body_solvers; backend_wake, backend_solve, backend_system,
                needs_induced_vorticity, update_trailing_edges=true,
                wakerow_no_hessian_to_particles,
                body_hessian_to_particles,
                body_gradient_core_size,
                body_on_wake,
                panel_wake_on_particles,
                particle_hessian_self,
                diagnose_particle_influence,
                diagnostic_vertical,
                grad_mu_options,
                formulation,
                formulation_state,
                i_step)
        else
            _kutta_step!(kutta_runtime, systems_tuple, wakes_tuple, frames, uinf;
                backend_wake, backend_solve, backend_system,
                needs_induced_vorticity, grad_mu_options, i_step,
                wakerow_no_hessian_to_particles,
                body_hessian_to_particles,
                body_gradient_core_size,
                body_on_wake,
                panel_wake_on_particles,
                particle_hessian_self)
        end

        # body bound-circulation low-pass (item 005 E4.8): damp body↔wake feedback
        # by under-relaxing the solved strength before it is shed into the wake.
        _step_timer_measure(:remaining_aerodynamics) do
            if apply_bound_rlx
                α = bound_strength_rlx
                for (sys, prev) in zip(systems_tuple, prev_strengths)
                    if have_prev_strength
                        @. sys.strength = (1 - α) * prev + α * sys.strength
                    end
                    prev .= sys.strength
                end
                if !isnothing(prev_c)
                    if have_prev_strength
                        @. formulation_state.c = (1 - α) * prev_c + α * formulation_state.c
                        set_wake_correction!(systems_tuple[1], formulation_state.c)
                    end
                    prev_c .= formulation_state.c
                end
                have_prev_strength = true
            end
        end

        #------- other solvers -------#

        # e.g. structures, acoustics, dynamics, etc.

        #------- update state -------#


        #------- monitors -------#

        # Run monitors before VTK write so monitor-owned fields can be passed to
        # downstream monitors and output files.
        _step_timer_measure(:monitors) do
            monitor_context = MonitorContext()
            monitor_set_time!(monitor_context, t)
            monitor_csv_dir = isnothing(path) ? nothing : joinpath(path, "monitors")
            for (i_monitor, monitor) in enumerate(monitors)
                _run_monitor!(monitor, monitor_context, systems_tuple, wakes_tuple, frames, uinf, i_step, dt, t)
                if !isnothing(monitor_csv_dir)
                    write_monitor_csv!(monitor, monitor_csv_dir, name, i_monitor,
                                       monitor_context, systems_tuple, i_step, dt;
                                       overwrite=i_step == start_step)
                end
            end
        end

        #------- save state -------#

        _step_timer_measure(:io) do
        if !isnothing(path)
            metadata_path = _metadata_toml_path(path, name)
            if i_step == start_step || !isfile(metadata_path)
                _write_metadata_toml(path, name, systems_tuple, wakes_tuple, frames, t_range,
                    body_solvers, backend_wake, backend_solve, backend_system, monitors;
                    start_step=start_step,
                    set_Das_eta_kinematic=set_Das_eta_kinematic,
                    set_Das_eta_freestream=set_Das_eta_freestream,
                    set_Das_min_kinematic_displacement=set_Das_min_kinematic_displacement,
                    clean_files=clean_files,
                    solver_options=(;
                        set_Das_kinematic_arc,
                        set_Das_refresh,
                        particle_relax,
                        body_on_wake,
                        panel_wake_on_particles,
                        particle_hessian_self,
                        body_hessian_to_particles,
                        body_gradient_core_size,
                        wakerow_no_hessian_to_particles,
                        bound_strength_rlx,
                    ),
                    kutta=isnothing(kutta_runtime) ? nothing :
                        _kutta_manifest_dict(kutta_runtime))
            end

            for (i, sys) in enumerate(systems_tuple)
                body_name = name * "_body$(i)"
                write_vtk(joinpath(path, body_name), sys, i_step, t;
                          monitors=monitors, i_system=i, overwrite=i_step==0,
                          compress=compress_vtk)
            end

            for (i, w) in enumerate(wakes_tuple)
                if !isnothing(w)
                    wake_name = name * "_wake$(i)"
                    write_vtk(joinpath(path, wake_name), w, i_step, t;
                              overwrite=i_step==0, compress=compress_vtk)
                end
            end

            _append_metadata_step_toml(path, name, frames, i_step, t; uinf,
                wakes=wakes_tuple,
                kutta=isnothing(kutta_runtime) ? nothing :
                    _kutta_step_dict(kutta_runtime))
        end
        end

        #------- propagate system -------#

        if i_step < length(t_range) - 1

            #--- state evolution ---#

            # propagate wake
            _step_timer_measure(:wake_propagation_maintenance) do
                for w in wakes_tuple
                    if w isa PanelParticleWake
                        propagate!(w, dt; relax=particle_relax,
                            step=i_step, frames, diagnose_particle_gamma, diagnostic_vertical)
                    elseif !isnothing(w)
                        propagate!(w, dt; step=i_step, frames)
                    end
                end
            end

            # Propagate rigid-body kinematics first. Persistent solver target
            # buffers consume control points, so mirror the rigid delta only
            # AFTER normals/control points have been refreshed below; doing it
            # here leaves those buffers one timestep behind the moved nodes.
            _step_timer_measure(:rigid_kinematics) do
                step_transforms = propagate_kinematics!(systems_tuple, frames, dt)

                # update control points and normals according to Neumann/Dirichlet BCs
                for sys in systems_tuple
                    calc_normals!(sys)
                    calc_controlpoints!(sys)
                end

                # Mirror the same rigid delta into persistent FMM state after all
                # kernel-consumed target geometry is current.
                transform_body_solvers!(body_solvers, systems_tuple, step_transforms)
            end

            #--- shed new wake ---#

            _step_timer_measure(:shedding) do
                for (sys, w) in zip(systems_tuple, wakes_tuple)
                    !isnothing(w) && shed_wake!(w, sys)
                end

                # Route B topology advancement bookkeeping (BRAINSTORM 015): the
                # accepted live block was just shifted into old-wake storage and
                # the fresh row-1 deposit is the reserved next live slot.
                isnothing(kutta_runtime) ||
                    _kutta_advance_topology!(kutta_runtime, i_step)
            end

        end

        _step_timer_finish_step!(_step_timer_token, i_step)

        # increment step
        i_step += 1
    end

    verbose && println()
end

get_Gammai(::AbstractBody{TK,NK,TF}) where {TK, NK, TF} = NK==2 ? 2 : 1
has_grad_mu(::AbstractBody{TK,NK,TF}) where {TK, NK, TF} = TK == ConstantDoublet || TK == VortexRing || TK == Union{ConstantSource, ConstantDoublet} || TK == Union{ConstantSource, VortexRing}

#------- wake shedding -------#

function update_TE!(wake::PanelWake, system::AbstractBody)

    # update first row based on system
    for i_surf in eachindex(wake.nodes)
        nodes = wake.nodes[i_surf]
        shedding = system.shedding[i_surf]
        Das = system.Das[i_surf]

        # loop over shedding panels
        for i_shed in axes(shedding, 2)
            i_panel = shedding[1, i_shed]
            idx_1 = shedding[3, i_shed] # nib (second node of the shedding edge)
            v1 = FastMultipole.SVector{3}(
                system.nodes[1, system.cells[idx_1, i_panel]],
                system.nodes[2, system.cells[idx_1, i_panel]],
                system.nodes[3, system.cells[idx_1, i_panel]],
            )
            nodes[:, 1, i_shed] .= v1 .+ view(Das, :, i_shed) # nib = vertex i_shed in Das
        end

        # final node of this edge (nia of last edge = last vertex in Das)
        i_panel = shedding[1, end]
        idx_2 = shedding[2, end] # nia (first node of the shedding edge)
        v2 = FastMultipole.SVector{3}(
            system.nodes[1, system.cells[idx_2, i_panel]],
            system.nodes[2, system.cells[idx_2, i_panel]],
            system.nodes[3, system.cells[idx_2, i_panel]],
        )
        nodes[:, 1, end] .= v2 .+ view(Das, :, size(Das, 2)) # nia of last edge = last vertex
    end
end

update_TE!(wake::PanelParticleWake, system::AbstractBody) = update_TE!(wake.panel_wake, system)

# update_TE!(w::ParticleWake, sys) = update_TE!(w.panel_wake, sys)

function shed_wake!(wake::PanelWake, system::AbstractBody)

    # check storage dimensions
    @assert length(system.Das) == length(system.shedding) == length(wake.nodes) "Length of system.Das ($(length(system.Das))) must match number of surfaces in wake ($(length(wake.nodes)))"

    # shift panels back a row
    n_rows = size(wake.nodes[1], 2)
    for i_surf in eachindex(wake.nodes)
        nodes = wake.nodes[i_surf]
        strength = wake.strength[i_surf]

        for j_row in min(wake.nwakes[]+1, n_rows-1):-1:1
            nodes[:, j_row+1, :] .= nodes[:, j_row, :]
        end
        for j_row in min(wake.nwakes[], n_rows-1):-1:1
            strength[:, j_row+1, :] .= strength[:, j_row, :]
        end
    end

    # update first row of strengths
    for i_surf in eachindex(wake.strength)
        for i_shed in axes(wake.strength[i_surf], 3)
            strengthi, strengthj = _get_wakestrength_mu(system, i_shed, i_surf)
            # Negate the strength because we are swapping the order of the nodes
            # to form a contiguous sequence (v1=nib, v4=nia instead of v1=nia, v4=nib)
            wake.strength[i_surf][1, 1, i_shed] = strengthi - strengthj
        end
    end

    # update nwakes
    if wake.nwakes[] == n_rows - 1 # about to overflow
        wake.overflowed[] = true
    end
    wake.nwakes[] = min(wake.nwakes[] + 1, n_rows - 1) # ensure we don't exceed storage

end

#------- PanelParticleWake shedding -------#

function shed_wake!(wake::PanelParticleWake, system::AbstractBody)
    pw = wake.panel_wake

    if pw.convert_at_shed
        # Convert-at-shed, nwakerows = 0 (BRAINSTORM 024): shed the row and
        # convert it to particles in the same call, so no free sheet survives
        # into the next solve and particles appear at the TE+Das line.
        pw.live_rows[] == 0 || error(
            "convert-at-shed (nwakerows = 0) does not support a reserved " *
            "live row block (Kutta Route B / TEAnchoredAttachment)")
        # 1. Strength history: the fresh row's downstream (unsteady) face jump
        #    and the retained-filament bookkeeping read row 2 as the previous
        #    shed's strength, but shed_wake!'s own strength shift is empty at
        #    nwakes[] == 0, so carry it explicitly before shedding.
        for i_surf in eachindex(pw.strength)
            s = pw.strength[i_surf]
            s[:, 2, :] .= s[:, 1, :]
        end
        # 2. Shift nodes (row 2 <- this step's convected row-1 line) and write
        #    the fresh row-1 strengths from the just-solved mu jump.
        shed_wake!(pw, system)
        # 3. Pin the fresh row's upstream edge to the current TE+Das line: the
        #    transient row spans the Das line -> last step's convected line.
        #    (For nwakerows >= 1 this pinning happens at the next solve's
        #    update_TE!; here the row is consumed before that.)
        update_TE!(pw, system)
        # 4. The single-row buffer "overflows" at every shed by construction;
        #    downstream filament/VTK guards key on this flag.
        pw.overflowed[] = true
        # 5. Convert the just-shed row (nwakes[] == capacity == 1). The smooth
        #    strategy's upstream (rigid-row) face jump is identically zero
        #    here — the rigid Das row carries the same just-solved mu — so
        #    :downstream attribution deposits the whole unsteady face and the
        #    retained filament on the Das line cancels the full row strength.
        _convert_to_particles!(wake, system)
        # 6. No free sheet survives into the next solve.
        pw.nwakes[] = 0
        return
    end

    n_rows = size(pw.nodes[1], 2)
    buffer_full = pw.nwakes[] >= n_rows - 1

    if buffer_full
        # new particles (the smooth strategy needs the body to complete its
        # streamwise stencil on a single-row wake; the legacy one ignores it)
        _convert_to_particles!(wake, system)
    end

    # Shift panel rows (existing PanelWake method)
    shed_wake!(pw, system)
end
