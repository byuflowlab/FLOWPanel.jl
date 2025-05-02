function induced(target, panel, kernel::AbstractRotatedKernel, derivatives_switch=DerivativesSwitch(true,true,true); sigma=1.0e-2)

    R = rotate_to_panel(panel)

    potential, velocity, velocity_gradient = _induced(target, panel.vertices, panel.control_point, panel.strength, kernel, sigma, R, derivatives_switch)

    # sigma = panel.radius * 0.1
    # potential, velocity, velocity_gradient = regularize(potential, velocity, velocity_gradient, target, panel.vertices, kernel, sigma)

    return potential, velocity, velocity_gradient
end

function induced(target, panel, kernel::AbstractUnrotatedKernel, derivatives_switch=DerivativesSwitch(true,true,true); sigma=1.0e-2)

    potential, velocity, velocity_gradient = _induced(target, panel, kernel, sigma, derivatives_switch)

    # sigma = panel.radius * 0.1
    # potential, velocity, velocity_gradient = regularize(potential, velocity, velocity_gradient, target, panel, kernel, sigma)

    return potential, velocity, velocity_gradient
end

# Function to calculate the distance from point P to the line segment AB
function minimum_distance(A, B, target)
    # compute vectos AB and Atarget
    ABx = B[1] - A[1]
    ABy = B[2] - A[2]
    ABz = B[3] - A[3]

    Atargetx = target[1] - A[1]
    Atargety = target[2] - A[2]
    Atargetz = target[3] - A[3]

    # Compute dot products directly (AB · AB and Atarget · AB)
    AB_dot_AB = ABx*ABx + ABy*ABy + ABz*ABz
    Atarget_dot_AB = Atargetx*ABx + Atargety*ABy + Atargetz*ABz

    # targetrojection length
    projection_length = Atarget_dot_AB / AB_dot_AB

    # Clamp projection_length to [0, 1]
    projection_length = clamp(projection_length, 0.0, 1.0)

    # Closest point calculation
    closest_x = A[1] + projection_length * ABx
    closest_y = A[2] + projection_length * ABy
    closest_z = A[3] + projection_length * ABz

    # Distance from target to the closest point
    dx = target[1] - closest_x
    dy = target[2] - closest_y
    dz = target[3] - closest_z

    return sqrt(dx*dx + dy*dy + dz*dz)
end

function minimum_distance(target, vertices)
    dist = eltype(target)(Inf)
    for i_vertex in 1:length(vertices)-1
        A = vertices[i_vertex]
        B = vertices[i_vertex+1]
        dist = min(dist, minimum_distance(A, B, target))
    end
    A = vertices[end]
    B = vertices[1]
    dist = min(dist, minimum_distance(A, B, target))

    return dist
end

@inline function rotate_to_panel(panel::Panel{TFP,<:Any,<:Any}) where TFP
    # unpack panel
    # (; normal, vertices) = panel
    normal = panel.normal
    vertices = panel.vertices

    return rotate_to_panel(normal, vertices)
end

@inline function rotate_to_panel(normal, vertices)
    # rotate into panel frame
    new_z = normal
    new_x = vertices[3] - vertices[1]
    new_x /= norm(new_x)
    new_y = cross(new_z, new_x)
    R = hcat(new_x, new_y, new_z)

    return R
end

@inline mysign(val) = val >= zero(val) ? 1 : -1

@inline function source_dipole_preliminaries(TFT, TFP, target, centroid, R)

    # promote types
    TF = promote_type(TFT,TFP)

    # target in the rotated frame
    target_Rx, target_Ry, target_Rz = transpose(R) * (target-centroid)
    if abs(target_Rz) < 10*eps()
        target_Rz = typeof(target_Rx)(10*eps()) * mysign(target_Rz)
    end

    # preallocate results
    potential = zero(TF)
    velocity = zero(SVector{3,TF})
    velocity_gradient = zero(SMatrix{3,3,TF,9})

    return potential, velocity, velocity_gradient, target_Rx, target_Ry, target_Rz
end

@inline function recurse_source_dipole(target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1)
    # eip1, hip1, rip1
    dx = target_Rx - vx_ip1
    dy = target_Ry - vy_ip1
    eip1 = dx * dx + target_Rz * target_Rz
    hip1 = dx * dy
    rip1 = sqrt(eip1 + dy * dy)

    # ei, hi, ri
    dx = target_Rx - vx_i
    dy = target_Ry - vy_i
    ei = dx * dx + target_Rz * target_Rz
    hi = dx * dy
    ri = sqrt(ei + dy * dy)

    # ds
    dx = vx_ip1 - vx_i
    dy = vy_ip1 - vy_i
    ds = sqrt(dx * dx + dy * dy)

    # mi
    mi = dy / dx

    return eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy
end

#------- constant source, normal doublet, source + normal doublet -------#

struct Source{M} <: AbstractRotatedKernel{M} end
const ConstantSource = Source{1}

struct NormalDoublet{M} <: AbstractRotatedKernel{M} end
const ConstantNormalDoublet = NormalDoublet{1}

struct SourceNormalDoublet{M} <: AbstractRotatedKernel{M} end
const ConstantSourceNormalDoublet = SourceNormalDoublet{2}

@inline function compute_source_dipole(::DerivativesSwitch{PS,VS,GS}, target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1, eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, strength::AbstractVector{TF}, ::ConstantSource, reg_term) where {PS,VS,GS,TF}

    #--- compute values ---#

    # intermediate quantities
    # singularity if probing on a side [ SOLVED ]; (easy way out is to perturb the evaluation point slightly)
    num = ri + rip1 - ds + reg_term # try adding the regularization term to the numerator
    # abs(num) < 5 * eps() && (num = typeof(num)(5*eps()))
    log_term = log(num / (ri + rip1 + ds))
    # singularity if d_z=0 [ SOLVED ] or probing at a vertex [ SOLVED ]; (easy way out is to perturb the evaluation point slightly)
    # iszero(ri) && (ri += eps())
    # iszero(rip1) && (rip1 += eps())
    # TODO: is reg_term needed here? maybe just for the gradient?
    tan_term = atan((mi * ei - hi) / (target_Rz * ri + reg_term * 0.0)) - atan((mi * eip1 - hip1) / (target_Rz * rip1 + reg_term * 0.0))

    potential = zero(TF)
    if PS# && !isinf(mi)
        potential += strength[1] * ((target_Rx - vx_i) * dy - (target_Ry - vy_i) * dx) / ds * log_term
        potential += strength[1] * target_Rz * tan_term
    end

    velocity = zero(SVector{3,TF})
    if VS
        velocity += strength[1] * SVector{3}(
            dy / ds * log_term,
            -dx / ds * log_term,
            tan_term
        )
    end

    velocity_gradient = zero(SMatrix{3,3,TF,9})
    if GS
        # intermediate values
        d2 = ds * dx
        r_plus_rp1 = ri + rip1
        r_plus_rp1_2 = r_plus_rp1 * r_plus_rp1
        r_times_rp1 = ri * rip1

        rho = r_times_rp1 + (target_Rx - vx_i) * (target_Rx - vx_ip1) + (target_Ry - vy_i) * (target_Ry - vy_ip1) + target_Rz * target_Rz
        lambda = (target_Rx - vx_i) * (target_Ry - vy_ip1) - (target_Rx - vx_ip1) * (target_Ry - vy_i)

        ri_inv = 1/ri
        rip1_inv = 1/rip1

        val1 = r_plus_rp1_2 - d2
        val2 = (target_Rx - vx_i) * ri_inv + (target_Rx - vx_ip1) * rip1_inv
        val3 = (target_Ry - vy_i) * ri_inv + (target_Ry - vy_ip1) * rip1_inv
        val4 = r_plus_rp1 / (r_times_rp1 * rho)

        # construct velocity_gradient
        phi_xx = 2 * dy / val1 * val2
        phi_xy = -2 * dx / val1 * val2
        phi_xz = target_Rz * dy * val4
        phi_yy = -2 * dx / val1 * val3
        phi_yz = -target_Rz * dx * val4
        phi_zz = lambda * val4
        velocity_gradient += strength[1] * SMatrix{3,3,eltype(velocity_gradient),9}(
            phi_xx, phi_xy, phi_xz,
            phi_xy, phi_yy, phi_yz,
            phi_xz, phi_yz, phi_zz
        )
    end

    return potential, velocity, velocity_gradient
end

@inline function compute_source_dipole(::DerivativesSwitch{PS,VS,GS}, target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1, eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, strength::AbstractVector{TF}, ::ConstantNormalDoublet, reg_term) where {PS,VS,GS,TF}

    # intermediate quantities
    # singularity if probing on a side [ SOLVED ]; (easy way out is to perturb the evaluation point slightly)
    # num = ri + rip1 - ds
    # abs(num) < 5*eps() && (num = typeof(num)(eps()))
    # log_term = log(num / (ri + rip1 + ds))
    # singularity if d_z=0 [ SOLVED ] or probing at a vertex [ SOLVED ]; (easy way out is to perturb the evaluation point slightly)
    # iszero(ri) && (ri += eps())
    # rip1 = rip1
    # iszero(rip1) && (rip1 += eps())
    tan_term = atan((mi * ei - hi) / (target_Rz * ri)) - atan((mi * eip1 - hi) / (target_Rz * rip1))
    # arg1 = (mi * es - hs) / (target_Rz * ri)
    # arg2 = (mi * esp1 - hs) / (target_Rz * rip1)
    # tan_term = atan((arg1 + arg2) / (1 - arg1 * arg2))
    # if arg1 * arg2 > 1.0
    #     if arg1 > 0
    #         tan_term += π
    #     else
    #         tan_term -= π
    #     end
    # end

    potential = -strength[1] * tan_term

    velocity = zero(SVector{3,TF})
    if VS
        r_plus_rp1 = ri + rip1
        r_times_rp1 = ri * rip1
        rho = r_times_rp1 + (target_Rx - vx_i) * (target_Rx - vx_ip1) + (target_Ry - vy_i) * (target_Ry - vy_ip1) + target_Rz * target_Rz
        lambda = (target_Rx - vx_i) * (target_Ry - vy_ip1) - (target_Rx - vx_ip1) * (target_Ry - vy_i)
        val4 = r_plus_rp1 / (r_times_rp1 * rho + reg_term)
        velocity += strength[1] * SVector{3}(
            target_Rz * dy * val4,
            -target_Rz * dx * val4,
            lambda * val4
        )
    end

    velocity_gradient = zero(SMatrix{3,3,TF,9})
    if GS
        # intermediate values
        d2 = ds * ds
        r_plus_rp1 = ri + rip1
        r_plus_rp1_2 = r_plus_rp1 * r_plus_rp1
        r_times_rp1 = ri * rip1

        rho = r_times_rp1 + (target_Rx - vx_i) * (target_Rx - vx_ip1) + (target_Ry - vy_i) * (target_Ry - vy_ip1) + target_Rz * target_Rz
        lambda = (target_Rx - vx_i) * (target_Ry - vy_ip1) - (target_Rx - vx_ip1) * (target_Ry - vy_i)

        val1 = r_times_rp1 * r_plus_rp1_2 + rho * rip1 * rip1
        val1 /= rho * ri * r_plus_rp1 + reg_term
        val2 = r_times_rp1 * r_plus_rp1_2 + rho * ri * ri
        val2 /= rho * rip1 * r_plus_rp1 + reg_term
        val3 = r_plus_rp1 / (rho * r_times_rp1 * r_times_rp1 + reg_term)

        # construct velocity_gradient
        psi_xx = target_Rz * dy * val3 * ((target_Rx - vx_i) * val1 + (target_Rx - vx_ip1) * val2)
        psi_xy = target_Rz * dy * val3 * ((target_Ry - vy_i) * val1 + (target_Ry - vy_ip1) * val2)
        psi_yy = -target_Rz * dx * val3 * ((target_Ry - vy_i) * val1 + (target_Ry - vy_ip1) * val2)
        val4 = r_plus_rp1_2 / rho
        val5 = (ri * ri - r_times_rp1 + rip1 * rip1) / (r_times_rp1 + reg_term)
        val6 = target_Rz * (val4 + val5)
        psi_zz = lambda * val3 * val6
        val7 = r_times_rp1 - target_Rz * val6
        val8 = val3 * val7
        psi_xz = -dy * val8
        psi_yz = dx * val8
        velocity_gradient += strength[1] * SMatrix{3,3,eltype(velocity_gradient),9}(
            psi_xx, psi_xy, psi_xz,
            psi_xy, psi_yy, psi_yz,
            psi_xz, psi_yz, psi_zz
        )
    end

    return potential, velocity, velocity_gradient
end

@inline function compute_source_dipole(::DerivativesSwitch{PS,VS,GS}, target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1, eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, strength::AbstractVector{TF}, ::ConstantSourceNormalDoublet, reg_term) where {PS,VS,GS,TF}

    # intermediate quantities
    # singularity if probing on a side [ SOLVED ]; (easy way out is to perturb the evaluation point slightly)
    num = ri + rip1 - ds + reg_term
    # abs(num) < 5*eps() && (num = typeof(num)(5*eps()))
    log_term = log(num / (ri + rip1 + ds))
    # singularity if d_z=0 [ SOLVED ] or probing at a vertex [ SOLVED ]; (easy way out is to perturb the evaluation point slightly)
    # ri = ri
    # iszero(ri) && (ri += eps())
    # rip1 = rip1
    # iszero(rip1) && (rip1 += eps())
    tan_term = atan((mi * ei - hi) / (target_Rz * ri)) - atan((mi * eip1 - hip1) / (target_Rz * rip1))

    potential = zero(TF)
    if PS# && !isinf(mi)
        potential += strength[1] * (((target_Rx - vx_i) * dy - (target_Ry - vy_i) * dx) / ds * log_term + target_Rz * tan_term)
        potential -= strength[2] * tan_term
    end

    velocity = zero(SVector{3,TF})
    if VS
        r_plus_rp1 = ri + rip1
        r_times_rp1 = ri * rip1
        rho = r_times_rp1 + (target_Rx - vx_i) * (target_Rx - vx_ip1) + (target_Ry - vy_i) * (target_Ry - vy_ip1) + target_Rz * target_Rz
        lambda = (target_Rx - vx_i) * (target_Ry - vy_ip1) - (target_Rx - vx_ip1) * (target_Ry - vy_i)
        val4 = r_plus_rp1 / (r_times_rp1 * rho + reg_term)
        velocity -= strength[1] * SVector{3}(
            dy / ds * log_term,
            -dx / ds * log_term,
            tan_term
        )
        velocity += strength[2] * SVector{3}(
            target_Rz * dy * val4,
            -target_Rz * dx * val4,
            lambda * val4
        )
    end

    velocity_gradient = zero(SMatrix{3,3,TF,9})
    if GS
        # intermediate values
        d2 = ds * ds
        r_plus_rp1 = ri + rip1
        r_plus_rp1_2 = r_plus_rp1 * r_plus_rp1
        r_times_rp1 = ri * rip1

        rho = r_times_rp1 + (target_Rx - vx_i) * (target_Rx - vx_ip1) + (target_Ry - vy_i) * (target_Ry - vy_ip1) + target_Rz * target_Rz
        lambda = (target_Rx - vx_i) * (target_Ry - vy_ip1) - (target_Rx - vx_ip1) * (target_Ry - vy_i)

        ri_inv = 1/(ri + reg_term)
        rip1_inv = 1/(rip1 + reg_term)

        val1 = r_plus_rp1_2 - d2
        val2 = (target_Rx - vx_i) * ri_inv + (target_Rx - vx_ip1) * rip1_inv
        val2 *= strength[1]
        val3 = (target_Ry - vy_i) * ri_inv + (target_Ry - vy_ip1) * rip1_inv
        val3 *= strength[1]
        val4 = r_plus_rp1 / (r_times_rp1 * rho + reg_term)
        val4 *= strength[1]

        # construct velocity_gradient
        phi_xx = 2 * dy / val1 * val2
        phi_xy = -2 * dx / val1 * val2
        phi_xz = target_Rz * dy * val4
        phi_yy = -2 * dx / val1 * val3
        phi_yz = -target_Rz * dx * val4
        phi_zz = lambda * val4
        velocity_gradient += SMatrix{3,3,eltype(velocity_gradient),9}(
            phi_xx, phi_xy, phi_xz,
            phi_xy, phi_yy, phi_yz,
            phi_xz, phi_yz, phi_zz
        )

        val1 = r_times_rp1 * r_plus_rp1_2 + rho * rip1 * rip1
        val1 /= rho * ri * r_plus_rp1
        val2 = r_times_rp1 * r_plus_rp1_2 + rho * ri * ri
        val2 /= rho * rip1 * r_plus_rp1
        val3 = r_plus_rp1 / (rho * r_times_rp1 * r_times_rp1 + reg_term)

        # construct velocity_gradient
        psi_xx = target_Rz * dy * val3 * ((target_Rx - vx_i) * val1 + (target_Rx - vx_ip1) * val2)
        psi_xy = target_Rz * dy * val3 * ((target_Ry - vy_i) * val1 + (target_Ry - vy_ip1) * val2)
        psi_yy = -target_Rz * dx * val3 * ((target_Ry - vy_i) * val1 + (target_Ry - vy_ip1) * val2)
        val4 = r_plus_rp1_2 / rho
        val5 = (ri * ri - r_times_rp1 + rip1 * rip1) / (r_times_rp1 + reg_term)
        val6 = target_Rz * (val4 + val5)
        psi_zz = lambda * val3 * val6
        val7 = r_times_rp1 - target_Rz * val6
        val8 = val3 * val7
        psi_xz = -dy * val8
        psi_yz = dx * val8
        velocity_gradient += strength[2] * SMatrix{3,3,eltype(velocity_gradient),9}(
            psi_xx, psi_xy, psi_xz,
            psi_xy, psi_yy, psi_yz,
            psi_xz, psi_yz, psi_zz
        )
    end

    return potential, velocity, velocity_gradient
end

function _induced(target::AbstractVector{TFT}, vertices::SVector{NS,<:Any}, centroid::AbstractVector{TFP}, strength, kernel::Union{ConstantSource, ConstantNormalDoublet, ConstantSourceNormalDoublet}, core_radius, R, derivatives_switch::DerivativesSwitch{PS,VS,GS}) where {TFT,TFP,NS,PS,VS,GS}
    #--- prelimilary computations ---#

    # note that target_Rz is ensured to be nonzero in the source_dipole_preliminaries function
    potential, velocity, velocity_gradient, target_Rx, target_Ry, target_Rz = source_dipole_preliminaries(TFT, TFP, target, centroid, R)

    #--- first recursive quantities ---#

    # current vertex locations
    vertex_ip1 = vertices[1] - centroid
    vx_ip1 = R[1,1] * vertex_ip1[1] + R[2,1] * vertex_ip1[2] + R[3,1] * vertex_ip1[3]
    vy_ip1 = R[1,2] * vertex_ip1[1] + R[2,2] * vertex_ip1[2] + R[3,2] * vertex_ip1[3]

    # loop over side contributions
    for i in 1:NS-1
        
        #--- recurse values ---#
        
        # current vertex locations
        vertex_i = vertex_ip1
        vx_i = vx_ip1
        vy_i = vy_ip1
        
        # the next vertex locations
        vertex_ip1 = vertices[i+1] - centroid
        vx_ip1 = R[1,1] * vertex_ip1[1] + R[2,1] * vertex_ip1[2] + R[3,1] * vertex_ip1[3]
        vy_ip1 = R[1,2] * vertex_ip1[1] + R[2,2] * vertex_ip1[2] + R[3,2] * vertex_ip1[3]

        #--- regularization term based on the minimum distance to this side ---#

        m_dist = minimum_distance(vertex_i, vertex_ip1, target)
        reg_term = regularize(m_dist, core_radius)
        
        #--- the rest ---#

        eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy = recurse_source_dipole(target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1)

        p, v, vg = compute_source_dipole(derivatives_switch, target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1, eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, strength, kernel, reg_term)
        if PS
            potential += p
        end
        if VS
            velocity += v
        end
        if GS
            velocity_gradient += vg
        end

    end

    #--- recurse values ---#

    # current vertex locations
    vertex_i = vertex_ip1
    vx_i = vx_ip1
    vy_i = vy_ip1

    # the next vertex locations
    vertex_ip1 = vertices[1] - centroid
    vx_ip1 = R[1,1] * vertex_ip1[1] + R[2,1] * vertex_ip1[2] + R[3,1] * vertex_ip1[3]
    vy_ip1 = R[1,2] * vertex_ip1[1] + R[2,2] * vertex_ip1[2] + R[3,2] * vertex_ip1[3]

    #--- regularization term based on the minimum distance to this side ---#

    m_dist = minimum_distance(vertex_i, vertex_ip1, target)
    reg_term = regularize(m_dist, core_radius)
    
    #--- the rest ---#

    eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy = recurse_source_dipole(target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1)

    #--- compute values ---#

    p, v, vg = compute_source_dipole(derivatives_switch, target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1, eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, strength, kernel, reg_term)

    #--- return result ---#

    if PS
        potential += p
        potential *= ONE_OVER_4PI
    end
    if VS
        velocity += v
        velocity = ONE_OVER_4PI * R * velocity
    end
    if GS
        velocity_gradient += vg
        velocity_gradient = -ONE_OVER_4PI * R * velocity_gradient * transpose(R)
    end

    return potential, velocity, velocity_gradient
end

"version for use with FastMultipole"
function _induced(target::AbstractVector{TFT}, vertices, centroid::AbstractVector{TFP}, strength, kernel::Union{ConstantSource, ConstantNormalDoublet, ConstantSourceNormalDoublet}, core_radius, R, derivatives_switch::DerivativesSwitch{PS,VS,GS}) where {TFT,TFP,NS,PS,VS,GS}
    #--- prelimilary computations ---#

    potential, velocity, velocity_gradient, target_Rx, target_Ry, target_Rz = source_dipole_preliminaries(TFT, TFP, target, centroid, R)

    #--- first recursive quantities ---#

    # current vertex locations
    vertex_ip1 = vertices[1] - centroid
    vx_ip1 = R[1,1] * vertex_ip1[1] + R[2,1] * vertex_ip1[2] + R[3,1] * vertex_ip1[3]
    vy_ip1 = R[1,2] * vertex_ip1[1] + R[2,2] * vertex_ip1[2] + R[3,2] * vertex_ip1[3]

    # loop over side contributions
    for i in 1:NS-1

        #--- recurse values ---#

        # current vertex locations
        vertex_i = vertex_ip1
        vx_i = vx_ip1
        vy_i = vy_ip1

        # the next vertex locations
        vertex_ip1 = vertices[i+1] - centroid
        vx_ip1 = R[1,1] * vertex_ip1[1] + R[2,1] * vertex_ip1[2] + R[3,1] * vertex_ip1[3]
        vy_ip1 = R[1,2] * vertex_ip1[1] + R[2,2] * vertex_ip1[2] + R[3,2] * vertex_ip1[3]

        # the rest
        eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy = recurse_source_dipole(target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1)

        p, v, vg = compute_source_dipole(derivatives_switch, target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1, eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, strength, kernel)
        if PS
            potential += p
        end
        if VS
            velocity += v
        end
        if GS
            velocity_gradient += vg
        end

    end

    #--- recurse values ---#

    # current vertex locations
    vertex_i = vertex_ip1
    vx_i = vx_ip1
    vy_i = vy_ip1

    # the next vertex locations
    vertex_ip1 = vertices[1] - centroid
    vx_ip1 = R[1,1] * vertex_ip1[1] + R[2,1] * vertex_ip1[2] + R[3,1] * vertex_ip1[3]
    vy_ip1 = R[1,2] * vertex_ip1[1] + R[2,2] * vertex_ip1[2] + R[3,2] * vertex_ip1[3]

    # the rest
    eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy = recurse_source_dipole(target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1)

    #--- compute values ---#

    p, v, vg = compute_source_dipole(derivatives_switch, target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1, eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, strength, kernel)

    #--- return result ---#

    if PS
        potential += p
        potential *= ONE_OVER_4PI
    end
    if VS
        velocity += v
        velocity = ONE_OVER_4PI * R * velocity
    end
    if GS
        velocity_gradient += vg
        velocity_gradient = -ONE_OVER_4PI * R * velocity_gradient * transpose(R)
    end

    return potential, velocity, velocity_gradient
end


#------- vortex ring panel -------#

struct Vortex{M} <: AbstractUnrotatedKernel{M} end
const VortexRing = Vortex{1}


function _induced(target::AbstractVector{TFT}, panel::Panel{TFP,<:Any,NS}, kernel::VortexRing, core_radius, derivatives_switch::DerivativesSwitch{PS,VS,GS}) where {TFT,TFP,NS,PS,VS,GS}
    TF = promote_type(TFT,TFP)
    corner_vectors = SVector{NS,SVector{3,TF}}(corner - target for corner in panel.vertices)
    velocity = zero(SVector{3,TF})
    gradient = zero(SMatrix{3,3,TF,9})

    # finite core settings
    finite_core = true
    core_size = core_radius

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

#------- regularization functions -------#

# function regularize(potential::TF, velocity, velocity_gradient, target, vertices, kernel::Union{NormalDoublet, SourceNormalDoublet}, sigma) where TF
#     distance = minimum_distance(target, vertices)
#     r_over_σ = distance / sigma
#     factor = one(TF) - exp(-r_over_σ * r_over_σ * r_over_σ)

#     # TODO: use the product rule to regularize the gradient also

#     return potential, velocity * factor, velocity_gradient * factor
# end

# function regularize(potential, velocity, velocity_gradient, target, panel, kernel, sigma)
#     return potential, velocity, velocity_gradient
# end

@inline function regularize(distance, core_size)
    δ = distance < core_size ? (distance-core_size) * (distance-core_size) : zero(distance)

    return δ
end