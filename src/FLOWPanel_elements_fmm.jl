#=##############################################################################
# DESCRIPTION
    Definition of panel elements

# AUTHORSHIP
  * Created by  : Ryan Anderson
  * Email       : Ry.M.Anderson@gmail.com
  * Date        : Jan 2026
  * License     : GNU Public License
=###############################################################################

#-------- panel kernels -------#

function rotate_to_panel(source_system::AbstractBody{<:Any,NK,<:Any}, source_buffer::Matrix{TF}, i_source::Int) where {TF,NK}

    #--- rotate into panel frame ---#

    # get vertices
    v1x = source_buffer[5+NK, i_source]
    v1y = source_buffer[6+NK, i_source]
    v1z = source_buffer[7+NK, i_source]
    v2x = source_buffer[8+NK, i_source]
    v2y = source_buffer[9+NK, i_source]
    v2z = source_buffer[10+NK, i_source]
    v3x = source_buffer[11+NK, i_source]
    v3y = source_buffer[12+NK, i_source]
    v3z = source_buffer[13+NK, i_source]

    return rotate_to_panel(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z)
end

function rotate_to_panel(v1x::TF, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z) where TF
    # normal (new z axis)
    # explicit cross(v2-v1, v3-v1)
    dx1 = v2x - v1x; dy1 = v2y - v1y; dz1 = v2z - v1z
    dx2 = v3x - v1x; dy2 = v3y - v1y; dz2 = v3z - v1z
    new_z_x = dy1 * dz2 - dz1 * dy2
    new_z_y = dz1 * dx2 - dx1 * dz2
    new_z_z = dx1 * dy2 - dy1 * dx2
    new_z_norm = sqrt(new_z_x*new_z_x + new_z_y*new_z_y + new_z_z*new_z_z)
    new_z_x /= new_z_norm
    new_z_y /= new_z_norm
    new_z_z /= new_z_norm

    # new x axis: (v3 - v1) / norm(v3-v1)
    new_x_x = v3x - v1x
    new_x_y = v3y - v1y
    new_x_z = v3z - v1z
    new_x_norm = sqrt(new_x_x*new_x_x + new_x_y*new_x_y + new_x_z*new_x_z)
    new_x_x /= new_x_norm
    new_x_y /= new_x_norm
    new_x_z /= new_x_norm

    # new y axis: (orthogonal to x and z)
    new_y_x = new_z_y*new_x_z - new_z_z*new_x_y
    new_y_y = new_z_z*new_x_x - new_z_x*new_x_z
    new_y_z = new_z_x*new_x_y - new_z_y*new_x_x
    R = FastMultipole.StaticArrays.SMatrix{3,3,TF,9}(new_x_x, new_x_y, new_x_z, 
                                                     new_y_x, new_y_y, new_y_z, 
                                                     new_z_x, new_z_y, new_z_z)

    # assemble vertices
    v1 = FastMultipole.StaticArrays.SVector{3,TF}(v1x, v1y, v1z)
    v2 = FastMultipole.StaticArrays.SVector{3,TF}(v2x, v2y, v2z)
    v3 = FastMultipole.StaticArrays.SVector{3,TF}(v3x, v3y, v3z)
    return R, v1, v2, v3
end

function rotate_to_panel(source_system::AbstractBody, i_source::Int)
    # get vertices
    v1i, v2i, v3i = source_system.cells[1, i_source], source_system.cells[2, i_source], source_system.cells[3, i_source]
    v1x, v1y, v1z = source_system.nodes[1, v1i], source_system.nodes[2, v1i], source_system.nodes[3, v1i]
    v2x, v2y, v2z = source_system.nodes[1, v2i], source_system.nodes[2, v2i], source_system.nodes[3, v2i]
    v3x, v3y, v3z = source_system.nodes[1, v3i], source_system.nodes[2, v3i], source_system.nodes[3, v3i]

    return rotate_to_panel(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z)
end

function get_vertices(source_system::AbstractBody{<:Any,NK,<:Any}, source_buffer::Matrix{TF}, i_source::Int) where {TF,NK}
    # get vertices
    v1x = source_buffer[5+NK, i_source]
    v1y = source_buffer[6+NK, i_source]
    v1z = source_buffer[7+NK, i_source]
    v2x = source_buffer[8+NK, i_source]
    v2y = source_buffer[9+NK, i_source]
    v2z = source_buffer[10+NK, i_source]
    v3x = source_buffer[11+NK, i_source]
    v3y = source_buffer[12+NK, i_source]
    v3z = source_buffer[13+NK, i_source]

    # assemble vertices
    v1 = FastMultipole.StaticArrays.SVector{3,TF}(v1x, v1y, v1z)
    v2 = FastMultipole.StaticArrays.SVector{3,TF}(v2x, v2y, v2z)
    v3 = FastMultipole.StaticArrays.SVector{3,TF}(v3x, v3y, v3z)
    return v1, v2, v3
end

function get_vertices(source_system::AbstractBody, i_source::Int)
    # get vertices
    v1i, v2i, v3i = source_system.cells[1, i_source], source_system.cells[2, i_source], source_system.cells[3, i_source]
    v1x, v1y, v1z = source_system.nodes[1, v1i], source_system.nodes[2, v1i], source_system.nodes[3, v1i]
    v2x, v2y, v2z = source_system.nodes[1, v2i], source_system.nodes[2, v2i], source_system.nodes[3, v2i]
    v3x, v3y, v3z = source_system.nodes[1, v3i], source_system.nodes[2, v3i], source_system.nodes[3, v3i]

    # assemble vertices
    TF = eltype(source_system.nodes)
    v1 = FastMultipole.StaticArrays.SVector{3,TF}(v1x, v1y, v1z)
    v2 = FastMultipole.StaticArrays.SVector{3,TF}(v2x, v2y, v2z)
    v3 = FastMultipole.StaticArrays.SVector{3,TF}(v3x, v3y, v3z)
    return v1, v2, v3
end

"Rotated kernels"
function induced(target::AbstractVector{TF}, source_system::AbstractBody{TK,NK,<:Any}, source_buffer::Matrix, i_source, derivatives_switch=FastMultipole.DerivativesSwitch(false,true,false); kerneloffset=1.0e-3) where {TF,TK,NK}

    # get vertices, rotation matrix
    R, v1, v2, v3 = rotate_to_panel(source_system, source_buffer, i_source)

    # get control point and strength
    control_point = FastMultipole.get_position(source_buffer, i_source)
    # strength = FastMultipole.get_strength(source_buffer, source_system, i_source)
    strength = FastMultipole.StaticArrays.SVector{NK,TF}(view(source_buffer, 5:4+NK, i_source))

    # evaluate influence
    potential, velocity, velocity_gradient = _induced(target, (v1, v2, v3), control_point, strength, TK, kerneloffset, R, derivatives_switch)

    # check for wake (if any)
    p, v, vg = _induced_wake(target, (v1, v2, v3), source_system, source_buffer, i_source, derivatives_switch)

    # return potential, velocity, velocity_gradient
    return potential+p, velocity+v, velocity_gradient+vg
end

function induced(target::AbstractVector{TF}, source_system::AbstractBody{TK,NK,<:Any}, i_source::Int, derivatives_switch=FastMultipole.DerivativesSwitch(false,true,false); kerneloffset=1.0e-3) where {TF,TK,NK}

    # get vertices, rotation matrix
    R, v1, v2, v3 = rotate_to_panel(source_system, i_source)

    # get control point and strength
    control_point = (v1 + v2 + v3) * 0.3333333333333333
    strength = FastMultipole.StaticArrays.SVector{NK,TF}(view(source_system.strength, i_source, :))

    # evaluate influence
    potential, velocity, velocity_gradient = _induced(target, (v1, v2, v3), control_point, strength, TK, kerneloffset, R, derivatives_switch)

    # check for wake (if any)
    p, v, vg = _induced_wake(target, (v1, v2, v3), source_system, i_source, derivatives_switch)

    # return potential, velocity, velocity_gradient
    return potential+p, velocity+v, velocity_gradient+vg
end

"Overload for non-rotated kernels"
function induced(target::AbstractVector{TF}, source_system::AbstractBody{VortexRing,NK,<:Any}, source_buffer::Matrix, i_source, derivatives_switch=FastMultipole.DerivativesSwitch(false,true,false); kerneloffset=1.0e-3) where {TF,NK}
    
    # get vertices
    v1, v2, v3 = get_vertices(source_system, source_buffer, i_source)
    
    # strength = FastMultipole.get_strength(source_buffer, source_system, i_source)
    strength = FastMultipole.StaticArrays.SVector{NK,TF}(view(source_buffer, 5:4+NK, i_source))

    potential, velocity, velocity_gradient = _induced(target, (v1, v2, v3), strength, VortexRing, kerneloffset, derivatives_switch)

    return potential, velocity, velocity_gradient
end

function induced(target::AbstractVector{TF}, source_system::AbstractBody{VortexRing,NK,<:Any}, i_source::Int, derivatives_switch=FastMultipole.DerivativesSwitch(false,true,false); kerneloffset=1.0e-3) where {TF,NK}
    
    # get vertices
    v1, v2, v3 = get_vertices(source_system, i_source)
    
    strength = FastMultipole.StaticArrays.SVector{NK,TF}(view(source_system.strength, i_source, :))

    potential, velocity, velocity_gradient = _induced(target, (v1, v2, v3), strength, VortexRing, kerneloffset, derivatives_switch)

    return potential, velocity, velocity_gradient
end

# Function to calculate the distance from point P to the line segment AB
function minimum_distance(A, B, target)
    # compute vectors AB and Atarget
    ABx = B[1] - A[1]
    ABy = B[2] - A[2]
    ABz = B[3] - A[3]

    Atargetx = target[1] - A[1]
    Atargety = target[2] - A[2]
    Atargetz = target[3] - A[3]

    # Compute dot products directly (AB · AB and Atarget · AB)
    AB_dot_AB = ABx*ABx + ABy*ABy + ABz*ABz
    Atarget_dot_AB = Atargetx*ABx + Atargety*ABy + Atargetz*ABz

    # target projection length
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

mysign(val) = val >= zero(val) ? 1 : -1

function source_dipole_preliminaries(TFT, TFP, target, centroid, R)

    # promote types
    TF = promote_type(TFT,TFP)

    # target in the rotated frame
    target_Rx, target_Ry, target_Rz = transpose(R) * (target-centroid)
    # if abs(target_Rz) < 10*eps()
    #     target_Rz = typeof(target_Rx)(10*eps()) * mysign(target_Rz)
    # end

    # preallocate results
    potential = zero(TF)
    velocity = zero(FastMultipole.StaticArrays.SVector{3,TF})
    velocity_gradient = zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})

    return potential, velocity, velocity_gradient, target_Rx, target_Ry, target_Rz
end

function recurse_source_dipole(target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1)
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

    # R_dot_s
    R_dot_s = dx * (target_Rx - vx_i) + dy * (target_Ry - vy_i)

    # mi
    mi = dy / dx

    return eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, R_dot_s
end

#------- constant source, normal doublet, source + normal doublet -------#

function compute_source_dipole(::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1, eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, strength::AbstractVector{TF}, ::Type{ConstantSource}, R_dot_s, reg_term) where {PS,VS,GS,TF}

    #--- compute values ---#

    # singularity if probing on a side
    num = max(eps(typeof(ri)), ri + rip1 - ds)
    log_term = log(num / (ri + rip1 + ds))

    # singularity at extension of the panel side
    if abs(abs(R_dot_s) - ri * ds) < 1e-12
        tan_term = zero(target_Rz)
    else
        # remove the singularity as much as possible
        arg1 = (dy * ei - hi * dx) / ri
        arg2 = (dy * eip1 - hip1 * dx) / rip1
        num = dx * target_Rz * (arg1 - arg2)
        den = target_Rz * target_Rz * dx * dx + arg1 * arg2
        tan_term = atan(num, den)
    end

    potential = zero(TF)
    if PS# && !isinf(mi)
        potential += strength[1] * ((target_Rx - vx_i) * dy - (target_Ry - vy_i) * dx) / ds * log_term
        potential += strength[1] * target_Rz * tan_term
    end

    velocity = zero(FastMultipole.StaticArrays.SVector{3,TF})
    if VS
        velocity += strength[1] * FastMultipole.StaticArrays.SVector{3}(
            dy / ds * log_term,
            -dx / ds * log_term,
            tan_term
        )
    end

    velocity_gradient = zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
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
        velocity_gradient += strength[1] * FastMultipole.StaticArrays.SMatrix{3,3,eltype(velocity_gradient),9}(
            phi_xx, phi_xy, phi_xz,
            phi_xy, phi_yy, phi_yz,
            phi_xz, phi_yz, phi_zz
        )
    end

    return potential, velocity, velocity_gradient
end

function compute_source_dipole(::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1, eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, strength::AbstractVector{TF}, ::Type{ConstantDoublet}, R_dot_s, reg_term) where {PS,VS,GS,TF}
    
    # singularity at extension of the panel side
    if abs(abs(R_dot_s) - ri * ds) < 1e-12
        tan_term = zero(target_Rz)
    else
        # remove the singularity as much as possible
        arg1 = (dy * ei - hi * dx) / ri
        arg2 = (dy * eip1 - hip1 * dx) / rip1
        num = dx * target_Rz * (arg1 - arg2)
        den = target_Rz * target_Rz * dx * dx + arg1 * arg2
        tan_term = atan(num, den)
    end

    potential = -strength[1] * tan_term

    velocity = zero(FastMultipole.StaticArrays.SVector{3,TF})
    if VS
        r_plus_rp1 = ri + rip1
        r_times_rp1 = ri * rip1
        rho = r_times_rp1 + (target_Rx - vx_i) * (target_Rx - vx_ip1) + (target_Ry - vy_i) * (target_Ry - vy_ip1) + target_Rz * target_Rz
        lambda = (target_Rx - vx_i) * (target_Ry - vy_ip1) - (target_Rx - vx_ip1) * (target_Ry - vy_i)
        val4 = r_plus_rp1 / (r_times_rp1 * rho + reg_term)
        velocity -= strength[1] * FastMultipole.StaticArrays.SVector{3}(
            target_Rz * dy * val4,
            -target_Rz * dx * val4,
            lambda * val4
        )
    end

    velocity_gradient = zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
    if GS
        # intermediate values
        d2 = ds * ds
        r_plus_rp1 = ri + rip1
        r_plus_rp1_2 = r_plus_rp1 * r_plus_rp1
        r_times_rp1 = ri * rip1

        rho = r_times_rp1 + (target_Rx - vx_i) * (target_Rx - vx_ip1) + (target_Ry - vy_i) * (target_Ry - vy_ip1) + target_Rz * target_Rz
        lambda = (target_Rx - vx_i) * (target_Ry - vy_ip1) - (target_Rx - vx_ip1) * (target_Ry - vy_i)

        val1 = r_times_rp1 * r_plus_rp1_2 + rho * rip1 * rip1
        val1 /= rho * ri * r_plus_rp1# + reg_term
        val2 = r_times_rp1 * r_plus_rp1_2 + rho * ri * ri
        val2 /= rho * rip1 * r_plus_rp1# + reg_term
        val3 = r_plus_rp1 / (rho * r_times_rp1 * r_times_rp1)# + reg_term)

        # construct velocity_gradient
        psi_xx = target_Rz * dy * val3 * ((target_Rx - vx_i) * val1 + (target_Rx - vx_ip1) * val2)
        psi_xy = target_Rz * dy * val3 * ((target_Ry - vy_i) * val1 + (target_Ry - vy_ip1) * val2)
        psi_yy = -target_Rz * dx * val3 * ((target_Ry - vy_i) * val1 + (target_Ry - vy_ip1) * val2)
        val4 = r_plus_rp1_2 / rho
        val5 = (ri * ri - r_times_rp1 + rip1 * rip1) / (r_times_rp1)# + reg_term)
        val6 = target_Rz * (val4 + val5)
        psi_zz = lambda * val3 * val6
        val7 = r_times_rp1 - target_Rz * val6
        val8 = val3 * val7
        psi_xz = -dy * val8
        psi_yz = dx * val8
        velocity_gradient += strength[1] * FastMultipole.StaticArrays.SMatrix{3,3,eltype(velocity_gradient),9}(
            psi_xx, psi_xy, psi_xz,
            psi_xy, psi_yy, psi_yz,
            psi_xz, psi_yz, psi_zz
        )
    end

    return potential, velocity, velocity_gradient
end

function compute_source_dipole(::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1, eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, strength::AbstractVector{TF}, ::Type{Union{ConstantSource, ConstantDoublet}}, R_dot_s, reg_term) where {PS,VS,GS,TF}

    # singularity if probing on a side
    # println("\nTESTING...")
    num = max(eps(typeof(ri)), ri + rip1 - ds)
    log_term = log(num / (ri + rip1 + ds))

    # singularity at extension of the panel side
    if abs(abs(R_dot_s) - ri * ds) < 1e-12
        tan_term = zero(target_Rz)
    else
        # println("NOT HERE")
        # remove the singularity as much as possible
        arg1 = (dy * ei - hi * dx) / ri
        arg2 = (dy * eip1 - hip1 * dx) / rip1
        num = dx * target_Rz * (arg1 - arg2)
        den = target_Rz * target_Rz * dx * dx + arg1 * arg2
        tan_term = atan(num, den)
    end

    potential = zero(TF)
    if PS# && !isinf(mi)
        potential += strength[1] * (((target_Rx - vx_i) * dy - (target_Ry - vy_i) * dx) / ds * log_term + target_Rz * tan_term)
        potential -= strength[2] * tan_term
    end

    velocity = zero(FastMultipole.StaticArrays.SVector{3,TF})
    if VS
        r_plus_rp1 = ri + rip1
        r_times_rp1 = ri * rip1
        rho = r_times_rp1 + (target_Rx - vx_i) * (target_Rx - vx_ip1) + (target_Ry - vy_i) * (target_Ry - vy_ip1) + target_Rz * target_Rz
        lambda = (target_Rx - vx_i) * (target_Ry - vy_ip1) - (target_Rx - vx_ip1) * (target_Ry - vy_i)

        val4 = r_plus_rp1 / (r_times_rp1 * rho + reg_term)
        velocity += strength[1] * FastMultipole.StaticArrays.SVector{3}(
            dy / ds * log_term,
            -dx / ds * log_term,
            tan_term
        )
        velocity -= strength[2] * FastMultipole.StaticArrays.SVector{3}(
            target_Rz * dy * val4,
            -target_Rz * dx * val4,
            lambda * val4
        )
    end

    velocity_gradient = zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
    if GS
        # intermediate values
        d2 = ds * ds
        r_plus_rp1 = ri + rip1
        r_plus_rp1_2 = r_plus_rp1 * r_plus_rp1
        r_times_rp1 = ri * rip1

        rho = r_times_rp1 + (target_Rx - vx_i) * (target_Rx - vx_ip1) + (target_Ry - vy_i) * (target_Ry - vy_ip1) + target_Rz * target_Rz
        lambda = (target_Rx - vx_i) * (target_Ry - vy_ip1) - (target_Rx - vx_ip1) * (target_Ry - vy_i)

        ri_inv = 1/(ri)# + reg_term)
        rip1_inv = 1/(rip1)# + reg_term)

        val1 = r_plus_rp1_2 - d2
        val2 = (target_Rx - vx_i) * ri_inv + (target_Rx - vx_ip1) * rip1_inv
        val2 *= strength[1]
        val3 = (target_Ry - vy_i) * ri_inv + (target_Ry - vy_ip1) * rip1_inv
        val3 *= strength[1]
        val4 = r_plus_rp1 / (r_times_rp1 * rho)# + reg_term)
        val4 *= strength[1]

        # construct velocity_gradient
        phi_xx = 2 * dy / val1 * val2
        phi_xy = -2 * dx / val1 * val2
        phi_xz = target_Rz * dy * val4
        phi_yy = -2 * dx / val1 * val3
        phi_yz = -target_Rz * dx * val4
        phi_zz = lambda * val4
        velocity_gradient += FastMultipole.StaticArrays.SMatrix{3,3,eltype(velocity_gradient),9}(
            phi_xx, phi_xy, phi_xz,
            phi_xy, phi_yy, phi_yz,
            phi_xz, phi_yz, phi_zz
        )

        val1 = r_times_rp1 * r_plus_rp1_2 + rho * rip1 * rip1
        val1 /= rho * ri * r_plus_rp1
        val2 = r_times_rp1 * r_plus_rp1_2 + rho * ri * ri
        val2 /= rho * rip1 * r_plus_rp1
        val3 = r_plus_rp1 / (rho * r_times_rp1 * r_times_rp1)# + reg_term)

        # construct velocity_gradient
        psi_xx = target_Rz * dy * val3 * ((target_Rx - vx_i) * val1 + (target_Rx - vx_ip1) * val2)
        psi_xy = target_Rz * dy * val3 * ((target_Ry - vy_i) * val1 + (target_Ry - vy_ip1) * val2)
        psi_yy = -target_Rz * dx * val3 * ((target_Ry - vy_i) * val1 + (target_Ry - vy_ip1) * val2)
        val4 = r_plus_rp1_2 / rho
        val5 = (ri * ri - r_times_rp1 + rip1 * rip1) / (r_times_rp1)# + reg_term)
        val6 = target_Rz * (val4 + val5)
        psi_zz = lambda * val3 * val6
        val7 = r_times_rp1 - target_Rz * val6
        val8 = val3 * val7
        psi_xz = -dy * val8
        psi_yz = dx * val8
        velocity_gradient += strength[2] * FastMultipole.StaticArrays.SMatrix{3,3,eltype(velocity_gradient),9}(
            psi_xx, psi_xy, psi_xz,
            psi_xy, psi_yy, psi_yz,
            psi_xz, psi_yz, psi_zz
        )
    end

    return potential, velocity, velocity_gradient
end

function _induced(target, vertices::NTuple{NS}, centroid::AbstractVector{TFP}, strength, kernel::Union{Type{ConstantSource}, Type{ConstantDoublet}, Type{Union{ConstantSource, ConstantDoublet}}}, core_radius, R, derivatives_switch::FastMultipole.DerivativesSwitch{PS,VS,GS}) where {TFP,NS,PS,VS,GS}
    #--- prelimilary computations ---#

    # note that target_Rz is ensured to be nonzero in the source_dipole_preliminaries function
    TFT = eltype(target)
    potential, velocity, velocity_gradient, target_Rx, target_Ry, target_Rz = source_dipole_preliminaries(TFT, TFP, target, centroid, R)

    # check if we're on the centroid
    to_centroid = target_Rx * target_Rx + target_Ry * target_Ry + target_Rz * target_Rz
    if to_centroid < eps(eltype(target))
        target_Rz += eps(eltype(target))
    end

    #--- first recursive quantities ---#

    # current vertex locations
    vertex_ip1 = vertices[1] - centroid
    vx_ip1 = R[1,1] * vertex_ip1[1] + R[2,1] * vertex_ip1[2] + R[3,1] * vertex_ip1[3]
    vy_ip1 = R[1,2] * vertex_ip1[1] + R[2,2] * vertex_ip1[2] + R[3,2] * vertex_ip1[3]

    # loop over side contributions
    for i in 1:NS
        
        #--- recurse values ---#
        
        # current vertex locations
        vertex_i = vertex_ip1
        vx_i = vx_ip1
        vy_i = vy_ip1
        
        # the next vertex locations
        ip1 = i < NS ? i+1 : 1
        vertex_ip1 = vertices[ip1] - centroid
        vx_ip1 = R[1,1] * vertex_ip1[1] + R[2,1] * vertex_ip1[2] + R[3,1] * vertex_ip1[3]
        vy_ip1 = R[1,2] * vertex_ip1[1] + R[2,2] * vertex_ip1[2] + R[3,2] * vertex_ip1[3]

        #--- regularization term based on the minimum distance to this side ---#

        m_dist = minimum_distance(vertex_i, vertex_ip1, target - centroid)
        reg_term = regularize(m_dist, core_radius)
        
        #--- the rest ---#

        eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, R_dot_s = recurse_source_dipole(target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1)

        p, v, vg = compute_source_dipole(derivatives_switch, target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1, eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, strength, kernel, R_dot_s, reg_term)
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

    # #--- recurse values ---#

    # # current vertex locations
    # vertex_i = vertex_ip1
    # vx_i = vx_ip1
    # vy_i = vy_ip1

    # # the next vertex locations
    # vertex_ip1 = vertices[1] - centroid
    # vx_ip1 = R[1,1] * vertex_ip1[1] + R[2,1] * vertex_ip1[2] + R[3,1] * vertex_ip1[3]
    # vy_ip1 = R[1,2] * vertex_ip1[1] + R[2,2] * vertex_ip1[2] + R[3,2] * vertex_ip1[3]

    # #--- regularization term based on the minimum distance to this side ---#

    # m_dist = minimum_distance(vertex_i, vertex_ip1, target - centroid)
    # reg_term = regularize(m_dist, core_radius)
    
    # #--- the rest ---#

    # eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, R_dot_s = recurse_source_dipole(target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1)

    # #--- compute values ---#

    # p, v, vg = compute_source_dipole(derivatives_switch, target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1, eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, strength, kernel, R_dot_s, reg_term)

    #--- return result ---#

    if PS
        # potential += p
        potential *= -ONE_OVER_4PI
    end
    if VS
        # velocity += v
        velocity = -ONE_OVER_4PI * R * velocity
    end
    if GS
        # velocity_gradient += vg
        velocity_gradient = -ONE_OVER_4PI * R * velocity_gradient * transpose(R)
    end

    return potential, velocity, velocity_gradient
end


#------- vortex ring panel -------#

function _induced(target, vertices::NTuple{NS}, strength, ::Type{VortexRing}, core_size, derivatives_switch::FastMultipole.DerivativesSwitch{PS,VS,GS}) where {NS,PS,VS,GS}
    TFT = eltype(target)
    TFP = eltype(strength)
    TF = promote_type(TFT,TFP)
    potential = zero(TF)
    velocity = zero(FastMultipole.StaticArrays.SVector{3,TF})
    gradient = zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})

    # finite core toggle
    finite_core = true
    
    # evaluate potential
    if PS
        @assert NS == 3
        v1, v2, v3 = vertices
        control_point = (v1 + v2 + v3) * 0.3333333333333333
        R, _ = rotate_to_panel(v1[1], v1[2], v1[3], v2[1], v2[2], v2[3], v3[1], v3[2], v3[3])
        p, _ = _induced(target, vertices, control_point, strength, ConstantDoublet, core_size, R, FastMultipole.DerivativesSwitch{true, false, false}())
        potential += p
    end

    # evaluate velocity/gradient
    if VS || GS
        for i in 1:NS
            ip1 = i < NS ? i+1 : 1

            # filament
            r1 = vertices[i] - target
            r2 = vertices[ip1] - target

            if VS
                # velocity
                # println("\nEvaluating Vortex Velocity")
                # @show vertices[i], vertices[ip1], target, strength[1]
                if strength[1] > 0.0 && false
                    DEBUG[] = true
                    @show vertices[i], vertices[ip1], target
                end
                v = _bound_vortex_velocity(r1, r2, finite_core, core_size)
                velocity += v
            end
            if GS
                # velocity gradient
                g = _bound_vortex_gradient(r1, r2, finite_core, core_size)
                gradient += g
            end
        end
    end

    return potential, velocity * strength[1], gradient * strength[1]
end

_induced(target, vertices, centroid, strength, kernel::Type{VortexRing}, core_size, R, derivatives_switch) =
    _induced(target, vertices, strength, kernel, core_size, derivatives_switch)


function _bound_vortex_velocity(r1::SVector{3,TF}, r2::SVector{3,TF}, finite_core, core_size) where TF
    # Vatistas n=2 core model: 1/h^2 → 1/sqrt(h^4 + rc^4)
    nr1 = norm(r1)
    nr2 = norm(r2)

    # target coincides with a filament endpoint: induced velocity is zero in the limit
    if nr1 < 5*eps(TF) || nr2 < 5*eps(TF)
        return zero(r1)
    end

    num = cross(r1, r2)
    r0 = r1 - r2
    dotrixrj = dot(num, num)        # |r1×r2|^2
    r0sqr = dot(r0, r0)             # |r0|^2
    rijdothat = dot(r0, r1/nr1 - r2/nr2)

    if finite_core
        V = num * rijdothat / sqrt(dotrixrj*dotrixrj + core_size*core_size*core_size*core_size * r0sqr*r0sqr) / (4*pi)
    else
        # singular kernel (no regularization)
        V = num * rijdothat / dotrixrj / (4*pi)
    end

    # if norm(V) > 500.0 && DEBUG[]
    #     @warn "V large!"
    #     @show V, nr1, nr2, finite_core, core_size
    #     stop
    #     println("============================================================================")
    # end

    return V
end

function _bound_vortex_gradient(r1::AbstractVector{TF}, r2, finite_core, core_size; epsilon=10*eps()) where TF
    # preliminaries
    r1norm = norm(r1)
    r2norm = norm(r2)

    # target coincides with a filament endpoint: induced gradient is zero in the limit
    if r1norm < 5*eps(TF) || r2norm < 5*eps(TF)
        return zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
    end

    rdot = dot(r1, r2)

    # zeta (with Vatistas-style core regularization: add core_size² to denominator)
    denom = r1norm*r2norm + rdot
    if finite_core
        denom += core_size*core_size
    end
    t1 = 1/denom
    t2 = 1/r1norm + 1/r2norm
    z = t1*t2*ONE_OVER_4PI

    # zeta gradient
    r1norm3 = r1norm*r1norm*r1norm
    r2norm3 = r2norm*r2norm*r2norm
    t4 = FastMultipole.StaticArrays.SVector{3,TF}(r1[i]/r1norm3 + r2[i]/r2norm3 for i in 1:3)
    t5 = FastMultipole.StaticArrays.SVector{3,TF}(r1norm/r2norm*r2[i] + r2norm/r1norm*r1[i] + r1[i] + r2[i] for i in 1:3)
    zgrad = ONE_OVER_4PI*(-t1*t4 - t2*t5*t1^2)

    # Omega
    o = cross(r1,r2)

    # Omega gradient
    ograd = FastMultipole.StaticArrays.SMatrix{3,3,TF,9}(
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

    gradient = zero(SMatrix{3,3,TF,9})

    return gradient
end

function _induced(target, vertices::NTuple, centroid::AbstractVector, strength, kernel::Type{Union{ConstantSource, VortexRing}}, core_radius, R, derivatives_switch::FastMultipole.DerivativesSwitch)

    # source influence
    p, v, vg = _induced(target, vertices, centroid, SVector{1}(strength[1]), ConstantSource, core_radius, R, derivatives_switch)

    # vortex ring
    p_vr, v_vr, vg_vr = _induced(target, vertices, SVector{1}(strength[2]), VortexRing, core_radius, derivatives_switch)

    return p + p_vr, v + v_vr, vg + vg_vr
end

#------- regularization functions -------#

@inline function regularize(distance, core_size)
    δ = distance < core_size ? (distance-core_size) * (distance-core_size) : zero(distance)

    return δ
end

#------- semi-infinite panels -------#

function get_strength_doublet(source_system::AbstractBody{Union{ConstantSource, ConstantDoublet}, 2, <:Any}, source_buffer, i_source)
    # get the strength of the doublet
    return source_buffer[6, i_source]
end

function get_strength_doublet(source_system::AbstractBody{<:Union{ConstantDoublet, VortexRing}, 1, <:Any}, source_buffer, i_source)
    # get the strength of the doublet
    return source_buffer[5, i_source]
end

function get_strength_doublet(source_system::AbstractBody{<:Union{ConstantSource, UniformVortexSheet}, 1, <:Any}, source_buffer, i_source)
    # get the strength of the doublet
    @warn "get_strength_doublet requested or a system that does not have doublets: returning zero"
    return zero(eltype(source_buffer))
end

function _induced_wake(target, (v1, v2, v3), source_system::AbstractBody, source_buffer, i_source, derivatives_switch)
    TF = eltype(target)
    return zero(TF), zero(FastMultipole.StaticArrays.SVector{3,TF}), zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
end

function get_strength_doublet(source_system::AbstractBody{Union{ConstantSource, ConstantDoublet}, 2, <:Any}, i_source::Int)
    # get the strength of the doublet
    return source_system.strength[i_source, 2]
end

function get_strength_doublet(source_system::AbstractBody{Union{ConstantSource, VortexRing}, 2, <:Any}, i_source::Int)
    # get the strength of the doublet
    return source_system.strength[i_source, 2]
end

function get_strength_doublet(source_system::AbstractBody{<:Union{ConstantDoublet, VortexRing}, 1, <:Any}, i_source::Int)
    # get the strength of the doublet
    return source_system.strength[i_source, 1]
end

function get_strength_doublet(source_system::AbstractBody{<:Union{ConstantSource, UniformVortexSheet}, 1, <:Any}, i_source::Int)
    # get the strength of the doublet
    @warn "get_strength_doublet requested or a system that does not have doublets: returning zero"
    return zero(eltype(source_system.strength))
end

function _induced_wake(target, (v1, v2, v3), source_system::AbstractBody, i_source::Int, derivatives_switch)
    TF = eltype(target)
    return zero(TF), zero(FastMultipole.StaticArrays.SVector{3,TF}), zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
end

get_wake_kernel(::AbstractBody{Union{ConstantSource, ConstantDoublet}}) = ConstantDoublet
get_wake_kernel(::AbstractBody{Union{ConstantSource, VortexRing}}) = VortexRing
get_wake_kernel(::AbstractBody{VortexRing}) = VortexRing
get_wake_kernel(::AbstractBody{ConstantDoublet}) = ConstantDoublet

function _induced_wake(target::AbstractVector{TF}, vertices::Tuple, source_system::RigidWakeBody, i_source::Int, derivatives_switch::FastMultipole.DerivativesSwitch{PS,VS,GS}) where {TF,PS,VS,GS}
    # check if this panel has a wake
    idx_1 = source_system.shedding_full[1, i_source]
    if idx_1 > 0
        # check which vertices are used
        v1i, v2i, v3i = source_system.cells[1, i_source], source_system.cells[2, i_source], source_system.cells[3, i_source]
        nodes_idx = (v1i, v2i, v3i)

        TE1 = nodes_idx[idx_1]
        v1 = FastMultipole.StaticArrays.SVector{3,TF}(source_system.nodes[1, TE1], source_system.nodes[2, TE1], source_system.nodes[3, TE1])
        v1x, v1y, v1z = v1[1], v1[2], v1[3]

        idx_2 = source_system.shedding_full[2, i_source]
        TE2 = nodes_idx[idx_2]
        v2 = FastMultipole.StaticArrays.SVector{3,TF}(source_system.nodes[1, TE2], source_system.nodes[2, TE2], source_system.nodes[3, TE2])
        v2x, v2y, v2z = v2[1], v2[2], v2[3]

        # get the wake shedding direction using correct Das columns
        i_surf = source_system.shedding_full[3, i_source]
        das_col_1 = source_system.shedding_full[5, i_source]  # Das column for TE1
        das_col_2 = source_system.shedding_full[6, i_source]  # Das column for TE2
        Dax, Day, Daz = source_system.Das[i_surf][1, das_col_1], source_system.Das[i_surf][2, das_col_1], source_system.Das[i_surf][3, das_col_1]

        # get strength
        strength = get_strength_doublet(source_system, i_source)
        TK = get_wake_kernel(source_system)

        # evaluate potential
        if source_system.semiinfinite_wake
            return induced_semiinfinite(target, TK, v1x, v1y, v1z, v2x, v2y, v2z, Dax, Day, Daz, strength, derivatives_switch; kerneloffset=source_system.kerneloffset)
        else
            # wake node connected to the first vertex (TE1)
            v1w_x = v1x + Dax
            v1w_y = v1y + Day
            v1w_z = v1z + Daz
            vw1 = FastMultipole.StaticArrays.SVector{3,TF}(v1w_x, v1w_y, v1w_z)

            # influence of first triangle
            control_point = (v1 + v2 + vw1) * 0.333333333333333
            strength_vec = FastMultipole.StaticArrays.SVector{1,TF}(strength)
            R, _ = rotate_to_panel(v1x, v1y, v1z, v2x, v2y, v2z, v1w_x, v1w_y, v1w_z)
            p, v, g = _induced(target, (v1, v2, vw1), control_point, strength_vec, TK, source_system.kerneloffset, R, derivatives_switch)

            # wake node connected to the second vertex (TE2)
            Dbx, Dby, Dbz = source_system.Das[i_surf][1, das_col_2], source_system.Das[i_surf][2, das_col_2], source_system.Das[i_surf][3, das_col_2]
            v2w_x = v2x + Dbx
            v2w_y = v2y + Dby
            v2w_z = v2z + Dbz
            vw2 = FastMultipole.StaticArrays.SVector{3,TF}(v2w_x, v2w_y, v2w_z)
            
            # influence of the second triangle
            control_point = (vw1 + v2 + vw2) * 0.333333333333333
            R, _ = rotate_to_panel(v1w_x, v1w_y, v1w_z, v2x, v2y, v2z, v2w_x, v2w_y, v2w_z)
            dp, dv, dg = _induced(target, (vw1, v2, vw2), control_point, strength_vec, TK, source_system.kerneloffset, R, derivatives_switch)
            if PS
                p += dp
            end
            if VS
                v += dv
            end
            if GS
                g += dg
            end

            return p, v, g
        end
    else
        # no wake for this panel
        return zero(TF), zero(FastMultipole.StaticArrays.SVector{3,TF}), zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
    end
end

function _induced_wake(target::AbstractVector{TF}, vertices::Tuple, source_system::RigidWakeBody{<:Any,NK,<:Any}, source_buffer::Matrix, i_source::Int, derivatives_switch::FastMultipole.DerivativesSwitch{PS,VS,GS}) where {TF,NK,PS,VS,GS}
    # Buffer layout (rows are 1-indexed, NK = number of element types):
    #   1-3:       center (cx, cy, cz)
    #   4:         radius
    #   5..4+NK:   strengths
    #   5+NK..13+NK: vertices v1, v2, v3 (3 rows each)
    #   14+NK:     idx1  -- local vertex index (1/2/3) of first TE node, or -1
    #   15+NK..17+NK: Da (wake direction for first TE node)
    #   18+NK:     idx2  -- local vertex index of second TE node, or -1
    #   19+NK..21+NK: Db (wake direction for second TE node)

    idx_1 = round(Int, source_buffer[14+NK, i_source])
    if idx_1 > 0
        # vertex row start: vertex i starts at row 2+NK+3*i
        r1 = 2 + NK + 3*idx_1
        v1x = source_buffer[r1,   i_source]
        v1y = source_buffer[r1+1, i_source]
        v1z = source_buffer[r1+2, i_source]

        idx_2 = round(Int, source_buffer[18+NK, i_source])
        r2 = 2 + NK + 3*idx_2
        v2x = source_buffer[r2,   i_source]
        v2y = source_buffer[r2+1, i_source]
        v2z = source_buffer[r2+2, i_source]

        v1 = FastMultipole.StaticArrays.SVector{3,TF}(v1x, v1y, v1z)
        v2 = FastMultipole.StaticArrays.SVector{3,TF}(v2x, v2y, v2z)

        # wake direction from buffer
        Dax = source_buffer[15+NK, i_source]
        Day = source_buffer[16+NK, i_source]
        Daz = source_buffer[17+NK, i_source]

        # doublet strength from buffer (NK-th element = last strength row)
        strength = source_buffer[4+NK, i_source]
        TK = get_wake_kernel(source_system)

        if source_system.semiinfinite_wake
            return induced_semiinfinite(target, TK, v1x, v1y, v1z, v2x, v2y, v2z, Dax, Day, Daz, strength, derivatives_switch; kerneloffset=source_system.kerneloffset)
        else
            # wake node connected to the first vertex
            v1w_x = v1x + Dax
            v1w_y = v1y + Day
            v1w_z = v1z + Daz
            vw1 = FastMultipole.StaticArrays.SVector{3,TF}(v1w_x, v1w_y, v1w_z)

            # wake direction for second node (stored at rows 19+NK..21+NK)
            Dbx = source_buffer[19+NK, i_source]
            Dby = source_buffer[20+NK, i_source]
            Dbz = source_buffer[21+NK, i_source]
            v2w_x = v2x + Dbx
            v2w_y = v2y + Dby
            v2w_z = v2z + Dbz
            vw2 = FastMultipole.StaticArrays.SVector{3,TF}(v2w_x, v2w_y, v2w_z)
            
            # strength
            strength_vec = FastMultipole.StaticArrays.SVector{1,TF}(strength)

            # induced influence due to the quad
            return _induced_quad(target, (v1, v2, vw2, vw1), strength_vec, TK, source_system.kerneloffset, derivatives_switch)

        end
    else
        return zero(TF), zero(FastMultipole.StaticArrays.SVector{3,TF}), zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
    end
end

function _induced_quad(target, vertices, strength, kernel::Type{ConstantDoublet}, kerneloffset, derivatives_switch)
    # influence of first triangle
    v1 = vertices[1]
    v2 = vertices[2]
    vw1 = vertices[4]
    control_point = (v1 + v2 + vw1) * 0.333333333333333
    R, _ = rotate_to_panel(v1x, v1y, v1z, v2x, v2y, v2z, v1w_x, v1w_y, v1w_z)
    p, vel, g = _induced(target, (v1, v2, vw1), control_point, strength_vec, TK, kerneloffset, R, derivatives_switch)

    # influence of the second triangle
    vw2 = vertices[3]
    control_point = (vw1 + v2 + vw2) * 0.333333333333333
    R, _ = rotate_to_panel(v1w_x, v1w_y, v1w_z, v2x, v2y, v2z, v2w_x, v2w_y, v2w_z)
    dp, dvel, dg = _induced(target, (vw1, v2, vw2), control_point, strength_vec, TK, kerneloffset, R, derivatives_switch)

    return p+dp, vel+dvel, g+dg
end

function _induced_quad(target, vertices, strength, kernel::Type{VortexRing}, kerneloffset, derivatives_switch::FastMultipole.DerivativesSwitch{PS,VS,GS}) where {PS,VS,GS}
    if PS
        # influence of first triangle
        v1 = vertices[1]
        v2 = vertices[2]
        vw1 = vertices[4]
        control_point = (v1 + v2 + vw1) * 0.333333333333333
        R, _ = rotate_to_panel(v1[1], v1[2], v1[3], v2[1], v2[2], v2[3], vw1[1], vw1[2], vw1[3])
        p, vel, g = _induced(target, (v1, v2, vw1), control_point, strength, kernel, kerneloffset, R, derivatives_switch)

        # influence of the second triangle
        vw2 = vertices[3]
        control_point = (vw1 + v2 + vw2) * 0.333333333333333
        R, _ = rotate_to_panel(vw1[1], vw1[2], vw1[3], v2[1], v2[2], v2[3], vw2[1], vw2[2], vw2[3])
        dp, dvel, dg = _induced(target, (vw1, v2, vw2), control_point, strength, kernel, kerneloffset, R, derivatives_switch)

        return p+dp, vel+dvel, g+dg
    else
        return _induced(target, vertices, strength, kernel, kerneloffset, derivatives_switch)
    end
end

function induced_semiinfinite(target::AbstractVector, TK::Type{VortexRing}, args...; kerneloffset)
    return induced_semiinfinite(target, ConstantDoublet, args...; kerneloffset)
end

function induced_semiinfinite(target::AbstractVector{TF}, TK::Type{ConstantDoublet}, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, strength, ::FastMultipole.DerivativesSwitch{PS,VS,GS}; kerneloffset) where {TF,PS,VS,GS}
    potential = zero(TF)
    velocity = zero(FastMultipole.StaticArrays.SVector{3,TF})
    gradient = zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
    if PS
        potential += _phi_semiinfinite(target, TK, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, strength; kerneloffset)
    end
    if VS
        velocity += _U_semiinfinite(target, TK, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, strength; kerneloffset)
    end

    return potential, velocity, gradient
end

function _phi_semiinfinite(target::AbstractVector{TF}, TK::Type{ConstantDoublet}, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, strength::Number; kerneloffset=1e-8) where TF

    # initialize result
    phi = zero(TF)
    scaling_factor = 1/(4*π)

    if abs(d1*d1 + d2*d2 + d3*d3 - 1) > 2*eps()
        error("Found non-unitary semi-infinite direction"*
                " norm([d1, d2, d3]) = norm($([d1,d2,d3])) = $(norm((d1,d2,d3)))")
    end

    # Split panel into a bound panel and a semi-infinite panel
    @inbounds begin
        p_i1, p_i2, p_i3 = v1x, v1y, v1z
        p_j1, p_j2, p_j3 = v2x, v2y, v2z

        # p_a = p_i + [(p_j-p_i)⋅d]d
        pijdotd = (p_j1 - p_i1)*d1 + (p_j2 - p_i2)*d2 + (p_j3 - p_i3)*d3
        pa1 = p_i1 + pijdotd*d1
        pa2 = p_i2 + pijdotd*d2
        pa3 = p_i3 + pijdotd*d3

        # Panel local coordinate system
        x1, x2, x3 = d1, d2, d3
        nrmpbpa = sqrt( (p_j1-pa1)^2 + (p_j2-pa2)^2 + (p_j3-pa3)^2  )
        y1, y2, y3 = (p_j1-pa1)/nrmpbpa, (p_j2-pa2)/nrmpbpa, (p_j3-pa3)/nrmpbpa
        z1, z2, z3 = x2*y3-x3*y2, x3*y1-x1*y3, x1*y2-x2*y1

        O1, O2, O3 = p_i1, p_i2, p_i3                   # Origin
        R = FastMultipole.StaticArrays.SMatrix{3,3,TF,9}(
            x1, x2, x3,
            y1, y2, y3,
            z1, z2, z3
        )                                               # Transformation matrix
    end

    # Target position in panel coordinate system
    # X = Oaxis*(targets[:, ti]-O)
    @inbounds begin
        x = x1*(target[1]-O1) + x2*(target[2]-O2) + x3*(target[3]-O3)
        y = y1*(target[1]-O1) + y2*(target[2]-O2) + y3*(target[3]-O3)
        z = z1*(target[1]-O1) + z2*(target[2]-O2) + z3*(target[3]-O3)
    end

    # ------------ Potential of bound panel
    if abs(pijdotd) > 2*eps()               # <--- Avoids the case that there is no bound panel

        # assemble vertices
        v1 = FastMultipole.StaticArrays.SVector{3,TF}(p_i1, p_i2, p_i3)
        v2 = FastMultipole.StaticArrays.SVector{3,TF}(p_j1, p_j2, p_j3)
        v3 = FastMultipole.StaticArrays.SVector{3,TF}(pa1, pa2, pa3)
        control_point = (v1 + v2 + v3) * 0.3333333333333333

        # assemble strength
        this_strength = FastMultipole.StaticArrays.SVector{1,TF}(strength)

        # other arguments
        derivatives_switch = FastMultipole.DerivativesSwitch{true, false, false}()

        # compute potential
        potential, _ = _induced(target, (v1, v2, v3), control_point, this_strength, TK, kerneloffset, R, derivatives_switch)
        phi += potential

    end

    # ------------ Potential of semi-infinite panel
    val = zero(TF)

    # Convert nodes to panel coordinate system
    xa = x1*(pa1-O1) + x2*(pa2-O2) + x3*(pa3-O3)
    ya = y1*(pa1-O1) + y2*(pa2-O2) + y3*(pa3-O3)
    # za = z1*(pa1-O1) + z2*(pa2-O2) + z3*(pa3-O3)
    xb = x1*(p_j1-O1) + x2*(p_j2-O2) + x3*(p_j3-O3)
    yb = y1*(p_j1-O1) + y2*(p_j2-O2) + y3*(p_j3-O3)
    # zb = z1*(p_j1-O1) + z2*(p_j2-O2) + z3*(p_j3-O3)

    # TODO: What is the domain of evaluation of this atan function in the theory?
    val += atan((yb-y) / z) + atan( (yb-y)*(x-xb) / (z*sqrt((x-xb)^2 + (yb-y)^2 + z^2)) )
    val -= atan((ya-y) / z) + atan( (ya-y)*(x-xa) / (z*sqrt((x-xa)^2 + (ya-y)^2 + z^2)) )
    # val += atan(yb-y, z) + atan( (yb-y)*(x-xa), z*sqrt((x-xa)^2 + (yb-y)^2 + z^2) )
    # val -= atan(ya-y, z) + atan( (ya-y)*(x-xa), z*sqrt((x-xa)^2 + (ya-y)^2 + z^2) )

    phi += strength * scaling_factor * val

    return phi
end

function _U_semiinfinite(target::AbstractVector{TF}, ::Type{ConstantDoublet}, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, strength::Number; kerneloffset=1e-8) where TF
    Ux1, Uy1, Uz1 = _U_semiinfinite_vortex(v1x, v1y, v1z,
                                    d1, d2, d3,
                                    -strength,
                                    target; offset=kerneloffset)
    
    Ux2, Uy2, Uz2 = _U_semiinfinite_vortex(v2x, v2y, v2z,
                                    d1, d2, d3,
                                    strength,
                                    target; offset=kerneloffset)
    
    Uxb, Uyb, Uzb = _U_boundvortex(v1x, v1y, v1z,
                                    v2x, v2y, v2z,
                                    strength,
                                    target; offset=kerneloffset)

    # Uxb, Uyb, Uzb = _bound_vortex_velocity(
    #                     FastMultipole.StaticArrays.SVector{3,TF}(v1x, v1y, v1z),
    #                     FastMultipole.StaticArrays.SVector{3,TF}(v2x, v2y, v2z),
    #                     true,
    #                     kerneloffset
    #                 ) .* strength
    
    # combine contributions
    Ux = Ux1 + Ux2 + Uxb
    Uy = Uy1 + Uy2 + Uyb
    Uz = Uz1 + Uz2 + Uzb

    return FastMultipole.StaticArrays.SVector{3,TF}(-Ux, -Uy, -Uz)
end

function _U_boundvortex( pa1::TF1, pa2::Number, pa3::Number,
                        pb1::Number, pb2::Number, pb3::Number,
                        strength::TF2,
                        target::AbstractVector{TF3};
                        cutoff=1e-14, offset=1e-8,
                       ) where{TF1, TF2, TF3}

    TF = promote_type(TF1, TF2, TF3)
    Ux = zero(TF)
    Uy = zero(TF)
    Uz = zero(TF)
    scaling_factor = 1/(4*π) # boundary

    # rij = pj - pi
    rij1 = pb1 - pa1
    rij2 = pb2 - pa2
    rij3 = pb3 - pa3

    @inbounds begin
        # ri = x - pi
        ri1 = target[1] - pa1
        ri2 = target[2] - pa2
        ri3 = target[3] - pa3

        # rj = x - pj
        rj1 = target[1] - pb1
        rj2 = target[2] - pb2
        rj3 = target[3] - pb3
    end

    # ri × rj
    rixrj1 = ri2*rj3 - ri3*rj2
    rixrj2 = ri3*rj1 - ri1*rj3
    rixrj3 = ri1*rj2 - ri2*rj1

    # ‖ ri × rj ‖^2
    dotrixrj = rixrj1^2 + rixrj2^2 + rixrj3^2

    # rij ⋅ (hat{ri} - hat{rj})
    normri = sqrt(ri1*ri1 + ri2*ri2 + ri3*ri3)
    normrj = sqrt(rj1*rj1 + rj2*rj2 + rj3*rj3)
    rijdothat = rij1*(ri1/normri - rj1/normrj) + rij2*(ri2/normri - rj2/normrj) + rij3*(ri3/normri - rj3/normrj)

    if dotrixrj > cutoff^2 # This makes the self induced velocity zero

        # Vatistas n=2 core model: 1/|r1×r2|^2 → 1/sqrt(|r1×r2|^4 + rc^4*|r0|^4)
        r0sqr = rij1^2 + rij2^2 + rij3^2
        aux = strength * rijdothat / sqrt(dotrixrj^2 + offset^4 * r0sqr^2)

        # NOTE: Negative sign is not needed since we defined rij = rj - ri
        Ux += aux * rixrj1 * scaling_factor
        Uy += aux * rixrj2 * scaling_factor
        Uz += aux * rixrj3 * scaling_factor
    end

    return -Ux, -Uy, -Uz
end

function _U_semiinfinite_vortex(p1::TF1, p2::Number, p3::Number,
                                d1::Number, d2::Number, d3::Number,
                                strength::TF2,
                                target::AbstractVector{TF3};
                                cutoff=1e-14, offset=1e-8) where {TF1,TF2,TF3}


    TF = promote_type(TF1,TF2,TF3)
    scaling_factor = 1/(4*π) # boundary, with negative to match Morino sign convention

    if abs(d1*d1 + d2*d2 + d3*d3 - 1) > 2*eps()
        error("Found non-unitary semi-infinite direction"*
                " norm([d1, d2, d3]) = norm($([d1,d2,d3])) = $(norm((d1,d2,d3)))")
    end

    vx, vy, vz = zero(TF), zero(TF), zero(TF)
    tx, ty, tz = target[1], target[2], target[3]

    # Split vortex into bound and semi-infinite sections
    # p0 = p + [(x-p)⋅d]d
    @inbounds xmpdotd = (tx - p1)*d1 + (ty - p2)*d2 + (tz - p3)*d3
    p01 = p1 + xmpdotd*d1
    p02 = p2 + xmpdotd*d2
    p03 = p3 + xmpdotd*d3

    # ----------------- Bound Vortex ---------------------------------------
    if (p01-p1) * (p01-p1) + (p02-p2) * (p02-p2) + (p03-p3) * (p03-p3) > offset^2 # Check that there is a bound section

        # rij = pj - pi
        rij1 = p01 - p1
        rij2 = p02 - p2
        rij3 = p03 - p3

        @inbounds begin
            # ri = x - pi
            ri1 = tx - p1
            ri2 = ty - p2
            ri3 = tz - p3

            # rj = x - pj
            rj1 = tx - p01
            rj2 = ty - p02
            rj3 = tz - p03
        end

        # ri × rj
        rixrj1 = ri2*rj3 - ri3*rj2
        rixrj2 = ri3*rj1 - ri1*rj3
        rixrj3 = ri1*rj2 - ri2*rj1

        # ‖ ri × rj ‖^2
        dotrixrj = rixrj1^2 + rixrj2^2 + rixrj3^2

        # rij ⋅ (hat{ri} - hat{rj})
        normri = sqrt(ri1^2 + ri2^2 + ri3^2)
        normrj = sqrt(rj1^2 + rj2^2 + rj3^2)
        rijdothat = rij1*(ri1/normri - rj1/normrj) + rij2*(ri2/normri - rj2/normrj) + rij3*(ri3/normri - rj3/normrj)

        if dotrixrj > cutoff^2 # This makes the self induced velocity zero

            # Vatistas n=2 core model: 1/|r1×r2|^2 → 1/sqrt(|r1×r2|^4 + rc^4*|r0|^4)
            r0sqr = rij1^2 + rij2^2 + rij3^2
            aux = strength * rijdothat / sqrt(dotrixrj^2 + offset^4 * r0sqr^2)

            # NOTE: Negative sign is not needed since we defined rij = rj - ri
            vx += aux * rixrj1 * scaling_factor
            vy += aux * rixrj2 * scaling_factor
            vz += aux * rixrj3 * scaling_factor
        end

    end

    # ----------------- Semi-Infinite Vortex -------------------------------
    # h = ‖x - p0‖
    # @inbounds h = sqrt( (tx - p01)^2 + (targets[2, ti] - p02)^2 + (targets[3, ti] - p03)^2 )
    hsqr = (tx - p01)^2 + (ty - p02)^2 + (tz - p03)^2

    # hhat = (x - p0) / ‖x - p0‖
    # h1 = (tx - p01)/h
    # h2 = (targets[2, ti] - p02)/h
    # h3 = (targets[3, ti] - p03)/h
    h1 = (tx - p01)
    h2 = (ty - p02)
    h3 = (tz - p03)

    # nhat = dhat × hhat
    n1 = d2*h3 - d3*h2
    n2 = d3*h1 - d1*h3
    n3 = d1*h2 - d2*h1

    # if h > cutoff # This makes the self induced velocity zero
    if hsqr > cutoff^2

        # Vatistas n=2 core model: 1/h^2 → 1/sqrt(h^4 + rc^4)
        aux = strength / sqrt(hsqr^2 + offset^4)

        vx += aux * n1 * scaling_factor
        vy += aux * n2 * scaling_factor
        vz += aux * n3 * scaling_factor

    end

    return -vx, -vy, -vz

end
