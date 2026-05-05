"""
    WingNormalization(rho, Sref, Lref)

Default normalization for `ForceMonitor`: divides force by `0.5 ρ |U∞|² Sref`
and moment by `0.5 ρ |U∞|² Sref Lref`.
"""
struct WingNormalization{TF}
    rho::TF
    Sref::TF
    Lref::TF
end

function (n::WingNormalization)(CF, CM, systems, frames, uinf)
    qinf = n.rho / 2 * (uinf[1]^2 + uinf[2]^2 + uinf[3]^2)
    return CF / (qinf * n.Sref), CM / (qinf * n.Sref * n.Lref)
end

"""
    NoNormalization()

Pass-through normalization for `ForceMonitor`: returns dimensional force and
moment unchanged.
"""
struct NoNormalization end

function (::NoNormalization)(CF, CM, systems, frames, uinf)
    return CF, CM
end

"""
    RotorNormalization(rho, D, i_frame)

Rotor thrust/torque coefficient normalization.  At each call, shaft speed `n`
(rev/s) is read from `frames[i_frame].ω / 2π`.  Force is divided by `ρ n² D⁴`
(giving CT) and moment by `ρ n² D⁵` (giving CQ), where `D` is rotor diameter.
"""
struct RotorNormalization{TF}
    rho::TF
    D::TF       # rotor diameter
    i_frame::Int
end

function (c::RotorNormalization)(CF, CM, systems, frames, uinf)
    n = frames[c.i_frame].ω / (2 * pi)   # rev/s
    ref = c.rho * n^2 * c.D^4
    return CF / ref, CM / (ref * c.D)
end

struct RotorNormalization2{TF}
    rho::TF
    D::TF       # rotor diameter
    i_frame::Int
end

function (c::RotorNormalization2)(CF, CM, systems, frames, uinf)
    n = frames[c.i_frame].ω # rad/s
    ref = c.rho * pi * c.D^2 * 0.25 * (n * c.D / 2)^2
    return CF / ref, CM / (ref * c.D*0.5)
end

"""
    ForceMonitor(nt, i_system; i_frame=-1, rho=1.0, Sref=1.0, Lref=1.0, TF=Float64, normalization=WingNormalization(...), verbose=false)

Storage monitor for integrated force and moment coefficient histories over `nt`
simulation steps.  Forces and moments on `systems[i_system]` are computed in
the coordinate system of `frames[i_frame]`.

`normalization` is called as `normalization(CF, CM, systems, frames, uinf)` and
must return a `(CF_norm, CM_norm)` tuple.  The default is `WingNormalization`
which divides by `0.5 ρ |U∞|² Sref` (and `… Lref` for moments).

If `verbose=true`, the normalized CF and CM for each step are printed to stdout
with a single `\\t` indent.
"""
struct ForceMonitor{TF, TN}
    force::Matrix{TF}
    moment::Matrix{TF}
    i_system::Int
    i_frame::Int
    normalization::TN
    verbose::Bool
end

"""
    KuttaJoukowskiForce(body, nt, i_system; rho=1.225, backend=DirectBackend(),
                        TF=Float64, normalization=NoNormalization(), verbose=false)

Diagnostic force monitor that integrates the Kutta–Joukowski force directly on
each panel edge:  ``F = ρ Σᵢ Σⱼ γᵢ \\, (Δsᵢⱼ × Vᵢⱼ)`` summed over panels `i`
and their three edges `j`, where `γᵢ` is the panel circulation (column of
`body.strength` named `"gamma"` or `"mu"`), `Δsᵢⱼ` is the directed edge vector,
and `Vᵢⱼ` is the inertial-frame fluid velocity (induced + freestream) at the
edge midpoint, evaluated via a `FastMultipole.ProbeSystem`.

Provides an independent cross-check against the pressure-integral force
returned by `ForceMonitor`/`calcfield_F!`. The body must use a `ConstantDoublet`
or `VortexRing` kernel.

If `verbose=true`, the normalized CF for each step is printed to stdout with a
single `\\t` indent.
"""
struct KuttaJoukowskiForce{TF, TB, TN}
    force::Matrix{TF}
    i_system::Int
    i_strength::Int
    rho::TF
    backend::TB
    probes::FastMultipole.ProbeSystem{TF}
    edge_node_a::Vector{Int}
    edge_node_b::Vector{Int}
    panel_of_probe::Vector{Int}
    normalization::TN
    verbose::Bool
end

function KuttaJoukowskiForce(body::AbstractBody, nt::Int, i_system::Int;
                              rho::Real=1.225,
                              backend::AbstractBackend=DirectBackend(),
                              TF=Float64,
                              normalization=NoNormalization(),
                              verbose::Bool=false)

    names = strength_names(body)
    i_strength = something(findfirst(==("gamma"), names),
                            findfirst(==("mu"), names),
                            0)
    i_strength == 0 && throw(ArgumentError(
        "KuttaJoukowskiForce requires a body with a ConstantDoublet or "*
        "VortexRing kernel; got strengths $(names)."))

    cells = body.cells
    @assert size(cells, 1) == 3 "KuttaJoukowskiForce assumes triangular panels (3 nodes per cell), got $(size(cells, 1))."
    ncells = body.ncells
    n_probes = 3 * ncells

    edge_node_a = Vector{Int}(undef, n_probes)
    edge_node_b = Vector{Int}(undef, n_probes)
    panel_of_probe = Vector{Int}(undef, n_probes)
    @inbounds for i in 1:ncells
        n1 = cells[1, i]; n2 = cells[2, i]; n3 = cells[3, i]
        k = 3*(i-1)
        edge_node_a[k+1] = n1; edge_node_b[k+1] = n2
        edge_node_a[k+2] = n2; edge_node_b[k+2] = n3
        edge_node_a[k+3] = n3; edge_node_b[k+3] = n1
        panel_of_probe[k+1] = i
        panel_of_probe[k+2] = i
        panel_of_probe[k+3] = i
    end

    probes = FastMultipole.ProbeSystem(n_probes, TF)
    force = zeros(TF, 3, nt)

    return KuttaJoukowskiForce{TF, typeof(backend), typeof(normalization)}(
        force, i_system, i_strength, TF(rho), backend, probes,
        edge_node_a, edge_node_b, panel_of_probe, normalization, verbose)
end

function (m::KuttaJoukowskiForce{TF})(systems, wakes,
                                       frames::AbstractVector{<:ReferenceFrame},
                                       uinf, i_step::Int) where {TF}
    body = systems[m.i_system]

    # 1. Update probe positions to edge midpoints; reset accumulators.
    zero_v = zero(SVector{3, TF})
    zero_h = zero(SMatrix{3, 3, TF, 9})
    @inbounds for k in eachindex(m.edge_node_a)
        a = m.edge_node_a[k]; b = m.edge_node_b[k]
        m.probes.position[k] = SVector{3, TF}(
            0.5 * (body.nodes[1, a] + body.nodes[1, b]),
            0.5 * (body.nodes[2, a] + body.nodes[2, b]),
            0.5 * (body.nodes[3, a] + body.nodes[3, b]),
        )
        m.probes.gradient[k] = zero_v
        m.probes.scalar_potential[k] = zero(TF)
        m.probes.hessian[k] = zero_h
    end

    # 2. Collect all sources: every body, plus every wake's source list.
    wake_sources = _collect_wake_sources(wakes)
    all_sources  = (systems..., wake_sources...)

    # 3. Compute induced velocity at probes.
    influence!((m.probes,), all_sources, m.backend;
                precalc=false,
                scalar_potential=false,
                gradient=true,
                hessian=(false,))

    # 4. Add freestream manually (no apply_freestream! method for ProbeSystem).
    uinf_sv = SVector{3, TF}(uinf[1], uinf[2], uinf[3])
    @inbounds for k in eachindex(m.probes.gradient)
        m.probes.gradient[k] += uinf_sv
    end

    # 5. Sum Kutta–Joukowski contributions: F = ρ Σ γ (Δs × V).
    Fsum = zero(SVector{3, TF})
    @inbounds for k in eachindex(m.edge_node_a)
        a = m.edge_node_a[k]; b = m.edge_node_b[k]
        Δs = SVector{3, TF}(
            body.nodes[1, b] - body.nodes[1, a],
            body.nodes[2, b] - body.nodes[2, a],
            body.nodes[3, b] - body.nodes[3, a],
        )
        γ = body.strength[m.panel_of_probe[k], m.i_strength]
        V = m.probes.gradient[k]
        Fsum += m.rho * γ * cross(Δs, V)
    end

    # 6. Optional normalization (moment is not computed; pass zero).
    Mzero = zero(SVector{3, TF})
    CF, _ = m.normalization(Fsum, Mzero, systems, frames, uinf)
    m.force[:, i_step + 1] .= CF

    if m.verbose
        println("\t\tKuttaJoukowskiForce[i_system=$(m.i_system), step=$(i_step+1)]:")
        println("\t\t\tCF = ($(round(CF[1], sigdigits=4)), $(round(CF[2], sigdigits=4)), $(round(CF[3], sigdigits=4)))")
    end
end


function ForceMonitor(nt::Int, i_system::Int;
                       i_frame=-1, TF=Float64,
                       normalization=WingNormalization(TF(rho), TF(Sref), TF(Lref)),
                       verbose::Bool=false)
    force = zeros(TF, 3, nt)
    moment = zeros(TF, 3, nt)
    return ForceMonitor{TF, typeof(normalization)}(force, moment, i_system, i_frame, normalization, verbose)
end

function (monitor::ForceMonitor)(systems, wakes,
                                  frames::AbstractVector{<:ReferenceFrame},
                                  uinf, i_step::Int)
    systems_tuple = _systems_tuple(systems)
    body = systems_tuple[monitor.i_system]

    # total force in global frame (body.F already populated by calcfield_F!)
    Ftot = calcfield_Ftot(body)
    Fvec = FastMultipole.SVector{3}(Ftot[1], Ftot[2], Ftot[3])

    if monitor.i_frame < 0
        # global frame: moment about the origin
        Mtot = calcfield_Mtot(body, zeros(3))
        Mvec = FastMultipole.SVector{3}(Mtot[1], Mtot[2], Mtot[3])
    else
        # frame-local: moment about frame origin, rotated into frame axes
        origin_global, R_f2g = frame_global_transform(frames, monitor.i_frame)
        Mtot = calcfield_Mtot(body, collect(origin_global))
        R_g2f = transpose(R_f2g)
        Fvec = R_g2f * Fvec
        Mvec = R_g2f * FastMultipole.SVector{3}(Mtot[1], Mtot[2], Mtot[3])
    end

    # normalise to coefficients
    CF_norm, CM_norm = monitor.normalization(Fvec, Mvec, systems_tuple, frames, uinf)
    monitor.force[:, i_step + 1] .= CF_norm
    monitor.moment[:, i_step + 1] .= CM_norm

    if monitor.verbose
        println("\t\tForceMonitor[i_system=$(monitor.i_system), step=$(i_step+1)]:")
        println("\t\t\tCF = ($(round(CF_norm[1], sigdigits=4)), $(round(CF_norm[2], sigdigits=4)), $(round(CF_norm[3], sigdigits=4)))")
        println("\t\t\tCM = ($(round(CM_norm[1], sigdigits=4)), $(round(CM_norm[2], sigdigits=4)), $(round(CM_norm[3], sigdigits=4)))")
    end
end
