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

"""
    ForceMonitor(nt, i_system; i_frame=-1, rho=1.0, Sref=1.0, Lref=1.0, TF=Float64, normalization=WingNormalization(...))

Storage monitor for integrated force and moment coefficient histories over `nt`
simulation steps.  Forces and moments on `systems[i_system]` are computed in
the coordinate system of `frames[i_frame]`.

`normalization` is called as `normalization(CF, CM, systems, frames, uinf)` and
must return a `(CF_norm, CM_norm)` tuple.  The default is `WingNormalization`
which divides by `0.5 ρ |U∞|² Sref` (and `… Lref` for moments).
"""
struct ForceMonitor{TF, TN}
    force::Matrix{TF}
    moment::Matrix{TF}
    i_system::Int
    i_frame::Int
    normalization::TN
end

function ForceMonitor(nt::Int, i_system::Int;
                       i_frame=-1, TF=Float64,
                       normalization=WingNormalization(TF(rho), TF(Sref), TF(Lref)))
    force = zeros(TF, 3, nt)
    moment = zeros(TF, 3, nt)
    return ForceMonitor{TF, typeof(normalization)}(force, moment, i_system, i_frame, normalization)
end

function (monitor::ForceMonitor)(systems::Tuple, wakes::Tuple,
                                  frames::AbstractVector{<:ReferenceFrame},
                                  uinf, i_step::Int)
    body = systems[monitor.i_system]

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
    CF_norm, CM_norm = monitor.normalization(Fvec, Mvec, systems, frames, uinf)
    monitor.force[:, i_step + 1] .= CF_norm
    monitor.moment[:, i_step + 1] .= CM_norm
end
