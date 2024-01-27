function induced(target, panel, kernel; toggle_potential=true, toggle_velocity=true, toggle_hessian=true)
    
    Rprime, Rxprime, Ryprime, Rzprime = rotate_to_panel(panel)

    potential, velocity, velocity_gradient = _induced(target, panel, kernel, Rprime, Rxprime, Ryprime, Rzprime; toggle_potential, toggle_velocity, toggle_hessian)
    
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
    (; normal, vertices) = panel
    
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

#####
##### constant source
#####

struct ConstantSource <: AbstractKernel end

@inline function _induced(target::AbstractVector{TFT}, panel::Panel{TFP,<:Any,NS}, kernel::ConstantSource, R, Rxprime, Ryprime, Rzprime; toggle_potential=true, toggle_velocity=true, toggle_hessian=true) where {TFT,TFP,NS}
    # unpack panel
    (;control_point, normal, vertices, strength) = panel
    centroid = control_point
    strength = strength[1]

    # promote types
    TF = promote_type(TFT,TFP)

    # rotate target
    target_R = R' * (target-centroid)

    # induced potential, velocity, gradient
    potential = zero(TF)
    velocity = @SVector zeros(TF,3)
    hessian = @SMatrix zeros(TF,3,3)

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
        
        if toggle_potential# && !isinf(ms[i])
            potential += (rxs[i] * dy - rys[i] * dx) / ds[i] * log_term
            potential += dz * tan_term
        end

        if toggle_velocity
            velocity += SVector{3}(
                dy / ds[i] * log_term,
                -dx / ds[i] * log_term,
                tan_term
            )
        end

        if toggle_hessian
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

            # construct hessian
            phi_xx = 2 * dy / val1 * val2
            phi_xy = -2 * dx / val1 * val2
            phi_xz = dz * dy * val4
            phi_yy = -2 * dx / val1 * val3
            phi_yz = -dz * dx * val4
            phi_zz = lambda * val4
            hessian += SMatrix{3,3,eltype(hessian),9}(
                phi_xx, phi_xy, phi_xz,
                phi_xy, phi_yy, phi_yz,
                phi_xz, phi_yz, phi_zz
            )
        end

    end

    potential *= strength * ONE_OVER_4PI
    velocity = -strength * ONE_OVER_4PI * R * velocity
    velocity_gradient = -strength * ONE_OVER_4PI * R * hessian * R'

    return potential, velocity, velocity_gradient
end

@inline kernel_multiplicity(::ConstantSource) = 1

#####
##### constant normal doublet
#####

struct ConstantNormalDoublet <: AbstractKernel end

function _induced(target::AbstractVector{TFT}, panel::Panel{TFP,<:Any,NS}, kernel::ConstantNormalDoublet, R, Rxprime, Ryprime, Rzprime; toggle_potential=true, toggle_velocity=true, toggle_hessian=true) where {TFT,TFP,NS}
    # unpack panel
    (;control_point, normal, vertices, strength) = panel
    centroid = control_point
    strength = strength[1]

    # promote types
    TF = promote_type(TFT,TFP)
    
    # rotate target
    target_R = R' * target

    # induced potential, velocity, gradient
    potential = zero(TF)
    velocity = @SVector zeros(TF,3)
    hessian = @SMatrix zeros(TF,3,3)
    
    # intermediate quantities
    dz = target_R[3] - Rzprime' * centroid
    iszero(dz) && (dz = eps())
    
    vertices_xR = get_vertices_xyR(Rxprime, vertices)
    vertices_yR = get_vertices_xyR(Ryprime, vertices)
    
    dxs = get_dxys(vertices_xR)
    dys = get_dxys(vertices_yR)

    ds = get_ds(dxs, dys)
    ms = get_ms(dxs, dys)
    
    dz2 = dz^2
    
    rxs = get_rxys(target_R[1], vertices_xR)
    rys = get_rxys(target_R[2], vertices_yR)

    es, hs = get_es_hs(rxs, rys, dz2)

    rs = get_rs(es, target_R[2], vertices_yR)
    rx_over_rs = get_rxy_over_rs(rxs, rs)
    ry_over_rs = get_rxy_over_rs(rys, rs)

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
        
        if toggle_potential# && !isinf(ms[i])
            potential -= tan_term
        end
        
        if toggle_velocity
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
            
        if toggle_hessian
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
            
            # construct hessian
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
            hessian += SMatrix{3,3,eltype(hessian),9}(
                psi_xx, psi_xy, psi_xz,
                psi_xy, psi_yy, psi_yz,
                psi_xz, psi_yz, psi_zz
            )
        end

    end

    potential *= strength * ONE_OVER_4PI
    velocity = strength * ONE_OVER_4PI * R * velocity
    velocity_gradient = -strength * ONE_OVER_4PI * R * hessian * R'

    return potential, velocity, velocity_gradient
end

@inline kernel_multiplicity(::ConstantNormalDoublet) = 1

#####
##### constant source plus constant normal doublet
#####

struct ConstantSourceNormalDoublet <: AbstractKernel end

function _induced(target::AbstractVector{TFT}, panel::Panel{TFP,<:Any,NS}, kernel::ConstantSourceNormalDoublet, R, Rxprime, Ryprime, Rzprime; toggle_potential=true, toggle_velocity=true, toggle_hessian=true) where {TFT,TFP,NS}
    # unpack panel
    (;control_point, normal, vertices, strength) = panel
    centroid = control_point
    source_strength, doublet_strength = strength

    # promote types
    TF = promote_type(TFT,TFP)
    
    # rotate target
    target_R = R' * target

    # induced potential, velocity, gradient
    potential = zero(TF)
    velocity = @SVector zeros(TF,3)
    hessian = @SMatrix zeros(TF,3,3)
    
    # intermediate quantities
    dz = target_R[3] - Rzprime' * centroid
    iszero(dz) && (dz = eps())
    
    vertices_xR = get_vertices_xyR(Rxprime, vertices)
    vertices_yR = get_vertices_xyR(Ryprime, vertices)
    
    dxs = get_dxys(vertices_xR)
    dys = get_dxys(vertices_yR)

    ds = get_ds(dxs, dys)
    ms = get_ms(dxs, dys)
    
    dz2 = dz^2
    
    rxs = get_rxys(target_R[1], vertices_xR)
    rys = get_rxys(target_R[2], vertices_yR)

    es, hs = get_es_hs(rxs, rys, dz2)

    rs = get_rs(es, target_R[2], vertices_yR)
    rx_over_rs = get_rxy_over_rs(rxs, rs)
    ry_over_rs = get_rxy_over_rs(rys, rs)

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
        
        if toggle_potential# && !isinf(ms[i])
            potential += source_strength * ((rxs[i] * dys[i] - rys[i] * dxs[i]) / ds[i] * log_term + dz * tan_term)
            potential -= doublet_strength * tan_term
        end
        
        if toggle_velocity
            r_plus_rp1 = rs[i] + rs[i+1]
            r_times_rp1 = rs[i] * rs[i+1]
            rho = r_times_rp1 + rxs[i] * rxs[i+1] + rys[i] * rys[i+1] + dz2
            lambda = rxs[i] * rys[i+1] - rxs[i+1] * rys[i]
            val4 = r_plus_rp1 / (r_times_rp1 * rho)
            velocity -= source_strength * SVector{3}(
                dys[i] / ds[i] * log_term,
                -dxs[i] / ds[i] * log_term,
                tan_term
            )
            velocity += doublet_strength * SVector{3}(
                dz * dys[i] * val4,
                -dz * dxs[i] * val4,
                lambda * val4
            )
        end
            
        if toggle_hessian
            # intermediate values
            d2 = ds[i]^2
            r_plus_rp1 = rs[i] + rs[i+1]
            r_plus_rp1_2 = r_plus_rp1^2
            r_times_rp1 = rs[i] * rs[i+1]
            
            rho = r_times_rp1 + rxs[i] * rxs[i+1] + rys[i] * rys[i+1] + dz2
            lambda = rxs[i] * rys[i+1] - rxs[i+1] * rys[i]

            val1 = r_plus_rp1_2 - d2
            val2 = rx_over_rs[i] + rx_over_rs[i+1]
            val2 *= source_strength
            val3 = ry_over_rs[i] + ry_over_rs[i+1]
            val3 *= source_strength
            val4 = r_plus_rp1 / (r_times_rp1 * rho)
            val4 *= source_strength

            # construct hessian
            phi_xx = 2 * dys[i] / val1 * val2
            phi_xy = -2 * dxs[i] / val1 * val2
            phi_xz = dz * dys[i] * val4
            phi_yy = -2 * dxs[i] / val1 * val3
            phi_yz = -dz * dxs[i] * val4
            phi_zz = lambda * val4
            hessian += SMatrix{3,3,eltype(hessian),9}(
                phi_xx, phi_xy, phi_xz,
                phi_xy, phi_yy, phi_yz,
                phi_xz, phi_yz, phi_zz
            )
            
            val1 = r_times_rp1 * r_plus_rp1_2 + rho * rs[i+1]^2
            val1 /= rho * rs[i] * r_plus_rp1
            val2 = r_times_rp1 * r_plus_rp1_2 + rho * rs[i]^2
            val2 /= rho * rs[i+1] * r_plus_rp1
            val3 = r_plus_rp1 / (rho * r_times_rp1^2)
            
            # construct hessian
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
            hessian += doublet_strength * SMatrix{3,3,eltype(hessian),9}(
                psi_xx, psi_xy, psi_xz,
                psi_xy, psi_yy, psi_yz,
                psi_xz, psi_yz, psi_zz
            )
        end

    end

    potential *= ONE_OVER_4PI
    velocity = ONE_OVER_4PI * R * velocity
    velocity_gradient = -ONE_OVER_4PI * R * hessian * R'

    return potential, velocity, velocity_gradient
end

@inline kernel_multiplicity(::ConstantSourceNormalDoublet) = 2