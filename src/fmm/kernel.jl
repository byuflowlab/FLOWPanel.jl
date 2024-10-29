function induced(target, panel, kernel::AbstractRotatedKernel, derivatives_switch=DerivativesSwitch(true,true,true))

    Rprime, Rxprime, Ryprime, Rzprime = rotate_to_panel(panel)

    potential, velocity, velocity_gradient = _induced(target, panel, kernel, Rprime, Rxprime, Ryprime, Rzprime, derivatives_switch)

    return potential, velocity, velocity_gradient
end

function induced(target, panel, kernel::AbstractUnrotatedKernel, derivatives_switch=DerivativesSwitch(true,true,true))

    potential, velocity, velocity_gradient = _induced(target, panel, kernel, derivatives_switch)

    return potential, velocity, velocity_gradient
end

@inline function get_vertices_xyR(Rxyprime::SVector{3,TF}, vertices::SVector{NS,<:SVector}, centroid) where {NS,TF}
    vec1 = SVector{NS,TF}(
        Rxyprime' * (vertices[i] - centroid) for i in 1:NS
    )
    return SVector{NS+1,TF}(vec1..., vec1[1])
end

@inline function get_dxys(vertices_dxyR::SVector{NSP1,TF}) where {NSP1,TF}
    return SVector{NSP1-1,TF}(
        vertices_dxyR[i+1] - vertices_dxyR[i] for i in 1:NSP1-1
    )
end

@inline function get_ds(dxs::SVector{NS,TF}, dys::SVector{NS,TF}) where {TF,NS}
    return SVector{NS,TF}(
        sqrt(dxs[i]^2 + dys[i]^2) for i in 1:NS
    )
end

@inline function get_ms(dxs::SVector{NS,TF}, dys::SVector{NS,TF}) where {NS,TF}
    return SVector{NS,TF}(
        dys[i] / dxs[i] for i in 1:NS
    )
end

@inline function get_rxys(target_Rxy::TFR, vertices_xyR::SVector{NSP1,TF}) where {TFR,NSP1,TF}
    return SVector{NSP1,promote_type(TFR,TF)}(
        target_Rxy - vertices_xyR[i] for i in 1:NSP1
    )
end

@inline function get_es_hs(rxs::SVector{NSP1,TF}, rys::SVector{NSP1,TF}, dz2) where {NSP1,TF}
    es = SVector{NSP1,TF}(
        rxs[i]^2 + dz2 for i in 1:NSP1
    )

    hs = SVector{NSP1,TF}(
        rxs[i] * rys[i] for i in 1:NSP1
    )

    return es, hs
end

@inline function get_rs(es::SVector{NSP1,TFT}, target_Ry::TFT, vertices_yR::SVector{NSP1,TF}) where {NSP1,TF,TFT}
    return SVector{NSP1,promote_type(TF,TFT)}(
        sqrt(es[i] + (vertices_yR[i] - target_Ry)^2) for i in 1:NSP1
    )
end

@inline function get_rxy_over_rs(rxys::SVector{NSP1,TF}, rs::SVector{NSP1,TF}) where {NSP1,TF}
    return SVector{NSP1,TF}(
        rxys[i] / rs[i] for i in 1:NSP1
    )
end

@inline function rotate_to_panel(panel::Panel{TFP,<:Any,<:Any}) where TFP
    # unpack panel
    # (; normal, vertices) = panel
    normal = panel.normal
    vertices = panel.vertices

    # rotate into panel frame
    new_z = normal
    new_x = vertices[3] - vertices[1]
    new_x /= norm(new_x)
    new_y = cross(new_z, new_x)
    R = hcat(new_x, new_y, new_z)
    Rprime = R'

    Rxprime = SVector{3,TFP}(Rprime[1,1],Rprime[1,2],Rprime[1,3])
    Ryprime = SVector{3,TFP}(Rprime[2,1],Rprime[2,2],Rprime[2,3])
    Rzprime = SVector{3,TFP}(Rprime[2,1],Rprime[2,2],Rprime[2,3])

    return R, Rxprime, Ryprime, Rzprime
end

@inline function source_dipole_preliminaries(target, panel, TFT, TFP, R, Rxprime, Ryprime)
    # unpack panel
    # (;control_point, normal, vertices, strength) = panel
    control_point = panel.control_point
    vertices = panel.vertices
    strength = panel.strength
    centroid = control_point

    # promote types
    TF = promote_type(TFT,TFP)

    # rotate target
    target_R = R' * (target-centroid)

    # induced potential, velocity, gradient
    potential = zero(TF)
    velocity = @SVector zeros(TF,3)
    velocity_gradient = @SMatrix zeros(TF,3,3)

    # intermediate quantities
    vertices_xR = get_vertices_xyR(Rxprime, vertices, centroid)
    vertices_yR = get_vertices_xyR(Ryprime, vertices, centroid)

    dxs = get_dxys(vertices_xR)
    dys = get_dxys(vertices_yR)

    ds = get_ds(dxs, dys)
    ms = get_ms(dxs, dys)

    dz = target_R[3]# - Rzprime' * centroid
    iszero(dz) && (dz = eps())

    dz2 = dz^2

    rxs = get_rxys(target_R[1], vertices_xR)
    rys = get_rxys(target_R[2], vertices_yR)

    es, hs = get_es_hs(rxs, rys, dz2)

    rs = get_rs(es, target_R[2], vertices_yR)
    rx_over_rs = get_rxy_over_rs(rxs, rs)
    ry_over_rs = get_rxy_over_rs(rys, rs)

    return strength, TF, potential, velocity, velocity_gradient, dxs, dys, ds, ms, dz, dz2, rxs, rys, es, hs, rs, rx_over_rs, ry_over_rs
end

#####
##### constant source
#####
struct Source{M} <: AbstractRotatedKernel{M} end
const ConstantSource = Source{1}

@inline function _induced(target::AbstractVector{TFT}, panel::Panel{TFP,<:Any,NS}, kernel::ConstantSource, R, Rxprime, Ryprime, Rzprime, derivatives_switch::DerivativesSwitch{PS,VS,GS}) where {TFT,TFP,NS,PS,VS,GS}
    # prelimilary computations
    strength, TF, potential, velocity, velocity_gradient, dxs, dys, ds, ms, dz, dz2, rxs, rys, es, hs, rs, rx_over_rs, ry_over_rs = source_dipole_preliminaries(target, panel, TFT, TFP, R, Rxprime, Ryprime)

    # loop over side contributions
    for i in 1:NS
        # intermediate quantities
        # singularity if probing on a side [ SOLVED ]; (easy way out is to perturb the evaluation point slightly)
        num = rs[i] + rs[i+1] - ds[i]
        iszero(num) && (num += eps())
        log_term = log(num / (rs[i] + rs[i+1] + ds[i]))
        # singularity if d_z=0 [ SOLVED ] or probing at a vertex [ SOLVED ]; (easy way out is to perturb the evaluation point slightly)
        ri = rs[i]
        iszero(ri) && (ri += eps())
        rip1 = rs[i+1]
        iszero(rip1) && (rip1 += eps())
        tan_term = atan((ms[i] * es[i] - hs[i]) / dz / ri) - atan((ms[i] * es[i+1] - hs[i+1]) / dz / rip1)
        dx = dxs[i]
        dy = dys[i]

        if PS# && !isinf(ms[i])
            potential += (rxs[i] * dy - rys[i] * dx) / ds[i] * log_term
            potential += dz * tan_term
        end

        if VS
            velocity += SVector{3}(
                dy / ds[i] * log_term,
                -dx / ds[i] * log_term,
                tan_term
            )
        end

        if GS
            # intermediate values
            d2 = ds[i]^2
            r_plus_rp1 = rs[i] + rs[i+1]
            r_plus_rp1_2 = r_plus_rp1^2
            r_times_rp1 = rs[i] * rs[i+1]

            rho = r_times_rp1 + rxs[i] * rxs[i+1] + rys[i] * rys[i+1] + dz2
            lambda = rxs[i] * rys[i+1] - rxs[i+1] * rys[i]

            val1 = r_plus_rp1_2 - d2
            val2 = rx_over_rs[i] + rx_over_rs[i+1]
            val3 = ry_over_rs[i] + ry_over_rs[i+1]
            val4 = r_plus_rp1 / (r_times_rp1 * rho)

            # construct velocity_gradient
            phi_xx = 2 * dy / val1 * val2
            phi_xy = -2 * dx / val1 * val2
            phi_xz = dz * dy * val4
            phi_yy = -2 * dx / val1 * val3
            phi_yz = -dz * dx * val4
            phi_zz = lambda * val4
            velocity_gradient += SMatrix{3,3,eltype(velocity_gradient),9}(
                phi_xx, phi_xy, phi_xz,
                phi_xy, phi_yy, phi_yz,
                phi_xz, phi_yz, phi_zz
            )
        end

    end

    if PS
        potential *= strength[1] * ONE_OVER_4PI
    end
    if VS
        velocity = -strength[1] * ONE_OVER_4PI * R * velocity
    end
    if GS
        velocity_gradient = -strength[1] * ONE_OVER_4PI * R * velocity_gradient * R'
    end

    return potential, velocity, velocity_gradient
end


#####
##### constant normal doublet
#####

struct NormalDoublet{M} <: AbstractRotatedKernel{M} end
const ConstantNormalDoublet = NormalDoublet{1}

function _induced(target::AbstractVector{TFT}, panel::Panel{TFP,<:Any,NS}, kernel::ConstantNormalDoublet, R, Rxprime, Ryprime, Rzprime, derivatives_switch::DerivativesSwitch{PS,VS,GS}) where {TFT,TFP,NS,PS,VS,GS}
    # prelimilary computations
    strength, TF, potential, velocity, velocity_gradient, dxs, dys, ds, ms, dz, dz2, rxs, rys, es, hs, rs, rx_over_rs, ry_over_rs = source_dipole_preliminaries(target, panel, TFT, TFP, R, Rxprime, Ryprime)

    # loop over side contributions
    for i in 1:NS
        # intermediate quantities
        # singularity if probing on a side [ SOLVED ]; (easy way out is to perturb the evaluation point slightly)
        num = rs[i] + rs[i+1] - ds[i]
        iszero(num) && (num += eps())
        log_term = log(num / (rs[i] + rs[i+1] + ds[i]))
        # singularity if d_z=0 [ SOLVED ] or probing at a vertex [ SOLVED ]; (easy way out is to perturb the evaluation point slightly)
        ri = rs[i]
        iszero(ri) && (ri += eps())
        rip1 = rs[i+1]
        iszero(rip1) && (rip1 += eps())
        tan_term = atan((ms[i] * es[i] - hs[i]) / dz / ri) - atan((ms[i] * es[i+1] - hs[i+1]) / dz / rip1)

        if PS# && !isinf(ms[i])
            potential -= tan_term
        end

        if VS
            r_plus_rp1 = rs[i] + rs[i+1]
            r_times_rp1 = rs[i] * rs[i+1]
            rho = r_times_rp1 + rxs[i] * rxs[i+1] + rys[i] * rys[i+1] + dz2
            lambda = rxs[i] * rys[i+1] - rxs[i+1] * rys[i]
            val4 = r_plus_rp1 / (r_times_rp1 * rho)
            velocity += SVector{3}(
                dz * dys[i] * val4,
                -dz * dxs[i] * val4,
                lambda * val4
            )
        end

        if GS
            # intermediate values
            d2 = ds[i]^2
            r_plus_rp1 = rs[i] + rs[i+1]
            r_plus_rp1_2 = r_plus_rp1^2
            r_times_rp1 = rs[i] * rs[i+1]

            rho = r_times_rp1 + rxs[i] * rxs[i+1] + rys[i] * rys[i+1] + dz2
            lambda = rxs[i] * rys[i+1] - rxs[i+1] * rys[i]

            val1 = r_times_rp1 * r_plus_rp1_2 + rho * rs[i+1]^2
            val1 /= rho * rs[i] * r_plus_rp1
            val2 = r_times_rp1 * r_plus_rp1_2 + rho * rs[i]^2
            val2 /= rho * rs[i+1] * r_plus_rp1
            val3 = r_plus_rp1 / (rho * r_times_rp1^2)

            # construct velocity_gradient
            psi_xx = dz * dys[i] * val3 * (rxs[i] * val1 + rxs[i+1] * val2)
            psi_xy = dz * dys[i] * val3 * (rys[i] * val1 + rys[i+1] * val2)
            psi_yy = -dz * dxs[i] * val3 * (rys[i] * val1 + rys[i+1] * val2)
            val4 = r_plus_rp1_2 / rho
            val5 = (rs[i]^2 - r_times_rp1 + rs[i+1]^2) / r_times_rp1
            val6 = dz * (val4 + val5)
            psi_zz = lambda * val3 * val6
            val7 = r_times_rp1 - dz * val6
            val8 = val3 * val7
            psi_xz = -dys[i] * val8
            psi_yz = dxs[i] * val8
            velocity_gradient += SMatrix{3,3,eltype(velocity_gradient),9}(
                psi_xx, psi_xy, psi_xz,
                psi_xy, psi_yy, psi_yz,
                psi_xz, psi_yz, psi_zz
            )
        end

    end

    if PS
        potential *= strength[1] * ONE_OVER_4PI
    end
    if VS
        velocity = strength[1] * ONE_OVER_4PI * R * velocity
    end
    if GS
        velocity_gradient = -strength[1] * ONE_OVER_4PI * R * velocity_gradient * R'
    end

    # dx = target - panel.control_point
    # if isapprox(dx' * dx, 0.0; atol=2*eps())
    #     velocity /= 2
    #     velocity_gradient /= 2
    # end

    return potential, velocity, velocity_gradient
end


#####
##### constant source plus constant normal doublet
#####

struct SourceNormalDoublet{M} <: AbstractRotatedKernel{M} end
const ConstantSourceNormalDoublet = SourceNormalDoublet{2}

function _induced(target::AbstractVector{TFT}, panel::Panel{TFP,<:Any,NS}, kernel::ConstantSourceNormalDoublet, R, Rxprime, Ryprime, Rzprime, derivatives_switch::DerivativesSwitch{PS,VS,GS}) where {TFT,TFP,NS,PS,VS,GS}
    # prelimilary computations;
    # note that strength[1] is the source strength and strength[2] is the dipole strength
    strength, TF, potential, velocity, velocity_gradient, dxs, dys, ds, ms, dz, dz2, rxs, rys, es, hs, rs, rx_over_rs, ry_over_rs = source_dipole_preliminaries(target, panel, TFT, TFP, R, Rxprime, Ryprime)

    # loop over side contributions
    for i in 1:NS
        # intermediate quantities
        # singularity if probing on a side [ SOLVED ]; (easy way out is to perturb the evaluation point slightly)
        num = rs[i] + rs[i+1] - ds[i]
        iszero(num) && (num += eps())
        log_term = log(num / (rs[i] + rs[i+1] + ds[i]))
        # singularity if d_z=0 [ SOLVED ] or probing at a vertex [ SOLVED ]; (easy way out is to perturb the evaluation point slightly)
        ri = rs[i]
        iszero(ri) && (ri += eps())
        rip1 = rs[i+1]
        iszero(rip1) && (rip1 += eps())
        tan_term = atan((ms[i] * es[i] - hs[i]) / dz / ri) - atan((ms[i] * es[i+1] - hs[i+1]) / dz / rip1)

        if PS# && !isinf(ms[i])
            potential += strength[1] * ((rxs[i] * dys[i] - rys[i] * dxs[i]) / ds[i] * log_term + dz * tan_term)
            potential -= strength[2] * tan_term
        end

        if VS
            r_plus_rp1 = rs[i] + rs[i+1]
            r_times_rp1 = rs[i] * rs[i+1]
            rho = r_times_rp1 + rxs[i] * rxs[i+1] + rys[i] * rys[i+1] + dz2
            lambda = rxs[i] * rys[i+1] - rxs[i+1] * rys[i]
            val4 = r_plus_rp1 / (r_times_rp1 * rho)
            velocity -= strength[1] * SVector{3}(
                dys[i] / ds[i] * log_term,
                -dxs[i] / ds[i] * log_term,
                tan_term
            )
            velocity += strength[2] * SVector{3}(
                dz * dys[i] * val4,
                -dz * dxs[i] * val4,
                lambda * val4
            )
        end

        if GS
            # intermediate values
            d2 = ds[i]^2
            r_plus_rp1 = rs[i] + rs[i+1]
            r_plus_rp1_2 = r_plus_rp1^2
            r_times_rp1 = rs[i] * rs[i+1]

            rho = r_times_rp1 + rxs[i] * rxs[i+1] + rys[i] * rys[i+1] + dz2
            lambda = rxs[i] * rys[i+1] - rxs[i+1] * rys[i]

            val1 = r_plus_rp1_2 - d2
            val2 = rx_over_rs[i] + rx_over_rs[i+1]
            val2 *= strength[1]
            val3 = ry_over_rs[i] + ry_over_rs[i+1]
            val3 *= strength[1]
            val4 = r_plus_rp1 / (r_times_rp1 * rho)
            val4 *= strength[1]

            # construct velocity_gradient
            phi_xx = 2 * dys[i] / val1 * val2
            phi_xy = -2 * dxs[i] / val1 * val2
            phi_xz = dz * dys[i] * val4
            phi_yy = -2 * dxs[i] / val1 * val3
            phi_yz = -dz * dxs[i] * val4
            phi_zz = lambda * val4
            velocity_gradient += SMatrix{3,3,eltype(velocity_gradient),9}(
                phi_xx, phi_xy, phi_xz,
                phi_xy, phi_yy, phi_yz,
                phi_xz, phi_yz, phi_zz
            )

            val1 = r_times_rp1 * r_plus_rp1_2 + rho * rs[i+1]^2
            val1 /= rho * rs[i] * r_plus_rp1
            val2 = r_times_rp1 * r_plus_rp1_2 + rho * rs[i]^2
            val2 /= rho * rs[i+1] * r_plus_rp1
            val3 = r_plus_rp1 / (rho * r_times_rp1^2)

            # construct velocity_gradient
            psi_xx = dz * dys[i] * val3 * (rxs[i] * val1 + rxs[i+1] * val2)
            psi_xy = dz * dys[i] * val3 * (rys[i] * val1 + rys[i+1] * val2)
            psi_yy = -dz * dxs[i] * val3 * (rys[i] * val1 + rys[i+1] * val2)
            val4 = r_plus_rp1_2 / rho
            val5 = (rs[i]^2 - r_times_rp1 + rs[i+1]^2) / r_times_rp1
            val6 = dz * (val4 + val5)
            psi_zz = lambda * val3 * val6
            val7 = r_times_rp1 - dz * val6
            val8 = val3 * val7
            psi_xz = -dys[i] * val8
            psi_yz = dxs[i] * val8
            velocity_gradient += strength[2] * SMatrix{3,3,eltype(velocity_gradient),9}(
                psi_xx, psi_xy, psi_xz,
                psi_xy, psi_yy, psi_yz,
                psi_xz, psi_yz, psi_zz
            )
        end

    end

    PS && (potential *= ONE_OVER_4PI)
    VS && (velocity = ONE_OVER_4PI * R * velocity)
    GS && (velocity_gradient = -ONE_OVER_4PI * R * velocity_gradient * R')

    return potential, velocity, velocity_gradient
end


#####
##### vortex ring panel
#####

struct Vortex{M} <: AbstractUnrotatedKernel{M} end
const VortexRing = Vortex{1}


function _induced(target::AbstractVector{TFT}, panel::Panel{TFP,<:Any,NS}, kernel::VortexRing, derivatives_switch::DerivativesSwitch{PS,VS,GS}) where {TFT,TFP,NS,PS,VS,GS}
    TF = promote_type(TFT,TFP)
    corner_vectors = SVector{NS,SVector{3,TF}}(corner - target for corner in panel.vertices)
    velocity = zero(SVector{3,TF})
    gradient = zero(SMatrix{3,3,TF,9})

    # finite core settings
    finite_core = false
    core_size = 1e-3

    # evaluate velocity/gradient
    for i in 1:NS-1
        r1 = panel.vertices[i] - target
        r2 = panel.vertices[i+1] - target

        # parameters
        r1norm = sqrt(r1'*r1)
        r2norm = sqrt(r2'*r2)
        r1normr2norm = r1norm*r2norm
        rcross = cross(r1, r2)
        rdot = dot(r1, r2)
        ONE_OVER_4PI = 1/4/pi

        if VS
            # velocity
            v = _bound_vortex_velocity(r1norm, r2norm, r1normr2norm, rcross, rdot, finite_core, core_size; epsilon=10*eps())
            velocity += v
        end
        if GS
            # velocity gradient
            g = _bound_vortex_gradient(r1, r2, r1norm, r2norm, r1normr2norm, rcross, rdot; epsilon=10*eps())
            gradient += g
        end
    end

    # wrap vertex
    r1 = panel.vertices[NS] - target
    r2 = panel.vertices[1] - target

    # parameters
    r1norm = sqrt(r1'*r1)
    r2norm = sqrt(r2'*r2)
    r1normr2norm = r1norm*r2norm
    rcross = cross(r1, r2)
    rdot = dot(r1, r2)
    ONE_OVER_4PI = 1/4/pi

    if VS
        # velocity
        v = _bound_vortex_velocity(r1norm, r2norm, r1normr2norm, rcross, rdot, finite_core, core_size; epsilon=10*eps())
        velocity += v
    end
    if GS
        # velocity gradient
        g = _bound_vortex_gradient(r1, r2, r1norm, r2norm, r1normr2norm, rcross, rdot; epsilon=10*eps())
        gradient += g
    end

    return zero(TF), velocity, gradient
end

function _bound_vortex_velocity(r1norm::TF, r2norm, r1normr2norm, rcross, rdot, finite_core::Bool, core_size::Number; epsilon=10*eps()) where TF
    # check if evaluation point is colinear with the bound vortex
    if norm(rcross) < epsilon # colinear if true
        if isapprox(rdot, -r1normr2norm; atol=epsilon) # at the midpoint, so return zero
            return zero(SVector{3,TF})
        elseif rdot <= 0.0 && finite_core # coincident with the filament so use the finite core model
            r1s, r2s, εs = r1norm^2, r2norm^2, core_size^2
            f1 = rcross/(r1s*r2s - rdot^2 + εs*(r1s + r2s - 2*r1normr2norm))
            f2 = (r1s - rdot)/sqrt(r1s + εs) + (r2s - rdot)/sqrt(r2s + εs)
            velocity = (f1*f2)/(4*pi)
            return velocity
        end
    end

    # otherwise, use singular kernel
    f1 = rcross/(r1normr2norm + rdot)
    f2 = (1/r1norm + 1/r2norm)

    velocity = (f1*f2) * ONE_OVER_4PI

    return velocity
end

function _bound_vortex_gradient(r1, r2, r1norm::TF, r2norm, r1normr2norm, rcross, rdot; epsilon=10*eps()) where TF
    # zeta
    t1 = 1/(r1norm*r2norm + rdot)
    t2 = 1/r1norm + 1/r2norm
    z = t1*t2*ONE_OVER_4PI

    # zeta gradient
    r1norm3 = r1norm^3
    r2norm3 = r2norm^3
    t4 = SVector{3,TF}(r1[i]/r1norm^3 + r2[i]/r2norm^3 for i in 1:3)
    t5 = SVector{3,TF}(r1norm/r2norm*r2[i] + r2norm/r1norm*r1[i] + r1[i] + r2[i] for i in 1:3)
    zgrad = ONE_OVER_4PI*(-t1*t4 - t2*t5*t1^2)

    # Omega
    o = cross(r1,r2)

    # Omega gradient
    ograd = SMatrix{3,3,TF,9}(
        0.0,# 1,1
        r1[3]-r2[3], # 2,1
        r2[2]-r1[2], # 3,1
        r2[3]-r1[3], # 1,2
        0.0, # 2,2
        r1[1]-r2[1], # 3,2
        r1[2]-r2[2], # 1,3
        r2[1]-r1[1], # 2,3
        0.0 # 3,3
    )
    gradient = transpose(zgrad * transpose(o)) + z * ograd

    return gradient
end
