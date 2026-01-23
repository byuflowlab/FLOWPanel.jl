#=##############################################################################
# DESCRIPTION
    Definition of panel elements

# AUTHORSHIP
  * Created by  : Ryan Anderson
  * Email       : Ry.M.Anderson@gmail.com
  * Date        : Jan 2026
  * License     : GNU Public License
=###############################################################################

#-------- high-level interface -------#

function _Uind!(self::AbstractBody, targets, out, backend::FastMultipoleBackend; optargs...)
    # wrap targets in a probe system
    TF = eltype(targets)
    potential = Vector{TF}(undef, 0) # unused
    hessian = Array{TF, 3}(undef, 0, 0, 0)  # unused
    probe_system = FastMultipole.ProbeSystemArray(targets, potential, out, hessian)

    # perform N-body calculation
    FastMultipole.fmm!(probe_system, _unpack_fmm(self); expansion_order=backend.expansion_order,
                                        multipole_acceptance=backend.multipole_acceptance,
                                        leaf_size_source=backend.leaf_size,
                                        hessian=false,
                                        gradient=true, 
                                        scalar_potential=false)

    _rigid_wake_U!(self, targets, out; optargs...)

    return nothing
end

function _phi!(self::AbstractBody, targets, out, backend::FastMultipoleBackend; optargs...)
    # wrap targets in a probe system
    TF = eltype(targets)
    velocity = Array{TF, 2}(undef, 0, 0)  # unused
    hessian = Array{TF, 3}(undef, 0, 0, 0)  # unused
    probe_system = FastMultipole.ProbeSystemArray(targets, out, velocity, hessian)
    
    # perform N-body calculation
    FastMultipole.fmm!(probe_system, _unpack_fmm(self); expansion_order=backend.expansion_order,
                                        multipole_acceptance=backend.multipole_acceptance,
                                        leaf_size_source=backend.leaf_size,
                                        hessian=false,
                                        gradient=false, 
                                        scalar_potential=true)

    _rigid_wake_phi!(self, targets, out; optargs...)

    return nothing
end

"Defaults to do nothing."
_rigid_wake_phi!(self::AbstractBody, targets, out; optargs...) = nothing

"Defaults to do nothing."
_rigid_wake_U!(self::AbstractBody, targets, out; optargs...) = nothing

"Recursively applies rigid wake contributions from each member body."
function _rigid_wake_U!(self::MultiBody, targets, out; optargs...)
    for body in self.bodies
        _rigid_wake_U!(body, targets, out; optargs...)
    end
end

"Defaults to do nothing."
_unpack_fmm(self::AbstractBody) = self

"Recursively unpacks multibody objects."
function _unpack_fmm(self::MultiBody)
    return Tuple(_unpack_fmm.(self.bodies))
end

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

"Rotated kernels"
function induced(target::AbstractVector{TF}, source_system::AbstractBody{TK,NK,<:Any}, source_buffer::Matrix, i_source, derivatives_switch=FastMultipole.DerivativesSwitch(false,true,false); kerneloffset=1.0e-3) where {TF,TK,NK}

    # get vertices, rotation matrix
    R, v1, v2, v3 = rotate_to_panel(source_system, source_buffer, i_source)

    # get control point and strength
    control_point = FastMultipole.get_position(source_buffer, i_source)
    # strength = FastMultipole.get_strength(source_buffer, source_system, i_source)
    strength = FastMultipole.StaticArrays.SVector{NK,TF}(view(source_buffer, 5:4+NK, i_source))

    potential, velocity, velocity_gradient = _induced(target, (v1, v2, v3), control_point, strength, TK, kerneloffset, R, derivatives_switch)

    return potential, velocity, velocity_gradient
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

    return -potential, -velocity, -velocity_gradient
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
        velocity += strength[1] * FastMultipole.StaticArrays.SVector{3}(
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
        velocity -= strength[1] * FastMultipole.StaticArrays.SVector{3}(
            dy / ds * log_term,
            -dx / ds * log_term,
            tan_term
        )
        velocity += strength[2] * FastMultipole.StaticArrays.SVector{3}(
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
        potential *= ONE_OVER_4PI
    end
    if VS
        # velocity += v
        velocity = ONE_OVER_4PI * R * velocity
    end
    if GS
        # velocity_gradient += vg
        velocity_gradient = -ONE_OVER_4PI * R * velocity_gradient * transpose(R)
    end

    return potential, velocity, velocity_gradient
end


#------- vortex ring panel -------#

function _induced(target, vertices::NTuple{NS}, strength, kernel::Type{VortexRing}, core_radius, derivatives_switch::FastMultipole.DerivativesSwitch{PS,VS,GS}) where {NS,PS,VS,GS}
    TFT = eltype(target)
    TFP = eltype(strength)
    TF = promote_type(TFT,TFP)
    velocity = zero(FastMultipole.StaticArrays.SVector{3,TF})
    gradient = zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})

    # finite core settings
    finite_core = true
    core_size = core_radius

    # evaluate velocity/gradient
    for i in 1:NS
        ip1 = i < NS ? i+1 : 1

        # filament
        r1 = vertices[i] - target
        r2 = vertices[ip1] - target

        if VS
            # velocity
            v = _bound_vortex_velocity(r1, r2, finite_core, core_size)
            velocity += v
        end
        if GS
            # velocity gradient
            g = _bound_vortex_gradient(r1, r2, finite_core, core_size)
            gradient += g
        end
    end

    return zero(TF), velocity * strength[1], gradient * strength[1]
end

@inline function get_δ(distance, core_size)
    δ = distance < core_size ? (distance-core_size) * (distance-core_size) : zero(distance)
    return δ
end

function _bound_vortex_velocity(r1, r2, finite_core, core_size)
    # intermediate values
    nr1 = norm(r1)
    nr2 = norm(r2)
    num = cross(r1, r2)
    denom = nr1 * nr2 + dot(r1, r2)

    if finite_core
        # core size comes into play here
        # nr3 = norm(r1 - r2) # length of the segment
        # distance_line = norm(num) / nr3
        δ0 = get_δ(denom, core_size)
        δ1 = get_δ(nr1, core_size) # distance to first endpoint
        δ2 = get_δ(nr2, core_size) # distance to second endpoint
        # nr3 = norm(r1 - r2) # length of the filament
        # distance_3 = sqrt(max(2*nr1*nr1 + 2*nr2*nr2 - nr3*nr3, zero(nr3)))*0.5 # distance to midpoint
        # δ3 = get_δ(distance_3, core_size)
    else
        δ0 = zero(denom)
        δ1 = zero(nr1)
        δ2 = zero(nr2)
    end

    # desingularized terms
    f1 = num/(denom + δ0)
    f2 = 1/(nr1+δ1)
    f3 = 1/(nr2+δ2)

    # evaluate vector field
    V = (f1*(f2+f3))/(4*pi)

    return V
end

function _bound_vortex_gradient(r1::AbstractVector{TF}, r2, finite_core, core_size; epsilon=10*eps()) where TF
    # preliminaries
    r1norm = norm(r1)
    r2norm = norm(r2)
    # r1normr2norm = r1norm * r2norm
    rdot = dot(r1, r2)
    # rcross = cross(r1, r2)

    # zeta
    t1 = 1/(r1norm*r2norm + rdot)
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

    return gradient
end

#------- regularization functions -------#

@inline function regularize(distance, core_size)
    δ = distance < core_size ? (distance-core_size) * (distance-core_size) : zero(distance)

    return δ
end