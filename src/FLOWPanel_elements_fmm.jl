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

"""
    rotate_to_panel(source_system, source_buffer, i_source)
    rotate_to_panel(source_system, i_source)
    rotate_to_panel(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z)

Return the panel-local rotation matrix and vertices for source panel `i_source`.
"""
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
    new_z_x /= (new_z_norm + eps(new_z_norm))
    new_z_y /= (new_z_norm + eps(new_z_norm))
    new_z_z /= (new_z_norm + eps(new_z_norm))

    # new x axis: (v3 - v1) / norm(v3-v1)
    new_x_x = v3x - v1x
    new_x_y = v3y - v1y
    new_x_z = v3z - v1z
    new_x_norm = sqrt(new_x_x*new_x_x + new_x_y*new_x_y + new_x_z*new_x_z)
    new_x_x /= (new_x_norm + eps(new_x_norm))
    new_x_y /= (new_x_norm + eps(new_x_norm))
    new_x_z /= (new_x_norm + eps(new_x_norm))

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

"""
    get_vertices(source_system, source_buffer, i_source)
    get_vertices(source_system, i_source)

Return the three panel vertices for source panel `i_source`.
"""
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

"""
    _self_limit(phi_raw, u_raw, ug_raw, ::Type{TK}, strength, n_gt)

Override the on-centroid `_induced` result with FLOWPanel's fixed self limit.

The evaluation point is NOT moved (`_induced` does not nudge target_Rz). At
`target_Rz == 0`, `compute_source_dipole` returns the principal value of the
solid-angle term (forces `tan_term = 0`), so the kernel outputs at a self pair
are clean PVs:
* doublet/ring potential PV = 0           (jump = ∓μ across the panel)
* source `u·n_GT` PV = 0                   (jump = ±σ across the panel)
* source potential PV: finite, from log_term integral (continuous)
* doublet/ring tangential & normal velocity PV: finite (continuous)

`_self_limit` therefore applies the fixed physical convention used by
FLOWPanel's exact control-point path:
* self-pair velocity is always the exterior surface limit, so source
  `u_s·n_GT = +σ/2`;
* self-pair potential is always the interior surface limit, so
  `φ_d = +μ/2` or `φ_Γ = +Γ/2`;
* continuous components keep their principal-value result.

For a self panel, the exterior doublet/ring potential can be recovered from the
returned interior value by subtracting the panel's doublet or vortex strength.

Consumers (`_G!`, Krylov matvec, post-processing) receive the actual
surface-limit value and do not add their own ±0.5 jump. Wake contributions from
`_induced_wake` are still added on top by the caller. The velocity gradient
self-limit is zeroed (out of scope).
"""
@inline function _self_limit(phi_raw, u_raw, ug_raw, ::Type{ConstantSource}, strength, n_gt)
    # source phi is continuous (PV in phi_raw is correct); source u·n_GT has the jump.
    sigma = strength[1]
    half_n = sigma * 0.5
    u_n_raw = u_raw[1]*n_gt[1] + u_raw[2]*n_gt[2] + u_raw[3]*n_gt[3]
    u_self = u_raw + (half_n - u_n_raw) * n_gt
    return phi_raw, u_self, zero(ug_raw)
end

@inline function _self_limit(phi_raw, u_raw, ug_raw, ::Type{ConstantDoublet}, strength, n_gt)
    # doublet velocity at the centroid is continuous (u_raw is the clean PV);
    # potential uses the interior half-jump (PV = 0).
    mu = strength[1]
    phi_self = mu * 0.5
    return phi_self, u_raw, zero(ug_raw)
end

@inline function _self_limit(phi_raw, u_raw, ug_raw, ::Type{VortexRing}, strength, n_gt)
    gamma = strength[1]
    phi_self = gamma * 0.5
    return phi_self, u_raw, zero(ug_raw)
end

@inline function _self_limit(phi_raw, u_raw, ug_raw, ::Type{Union{ConstantSource, ConstantDoublet}}, strength, n_gt)
    # phi_raw = source PV (clean) + 0 (doublet PV).  Add interior doublet half-jump.
    # u_raw   = doublet PV (clean) + source tangential PV (clean) + 0 (source normal PV).
    # Add only the exterior source half-jump so the continuous doublet velocity is preserved.
    sigma = strength[1]
    mu    = strength[2]
    phi_self = phi_raw + mu * 0.5
    half_n   = sigma * 0.5
    u_self   = u_raw + half_n * n_gt
    return phi_self, u_self, zero(ug_raw)
end

@inline function _self_limit(phi_raw, u_raw, ug_raw, ::Type{Union{ConstantSource, VortexRing}}, strength, n_gt)
    sigma = strength[1]
    gamma = strength[2]
    phi_self = phi_raw + gamma * 0.5
    half_n   = sigma * 0.5
    u_self   = u_raw + half_n * n_gt
    return phi_self, u_self, zero(ug_raw)
end

# Relative tolerance for self-pair detection: ||target - centroid|| < ε_rel * √A.
const SELF_PAIR_EPS_REL = 1.0e-12

# On-plane snap for the solid-angle PV branch (see _induced): plain floats
# snap so the ±2π branch side is deterministic across host/device/exact
# arithmetic; AD duals pass through unchanged so their partials survive.
@inline _onplane_snap(tRz::T, L2) where {T<:Union{Float32,Float64}} =
    ifelse(tRz*tRz <= 1e-24 * L2, zero(T), tRz)
@inline _onplane_snap(tRz, L2) = tRz

@inline function _is_self_pair(target, control_point, vertices)
    v1, v2, v3 = vertices
    e1 = v2 - v1
    e2 = v3 - v1
    nx = e1[2]*e2[3] - e1[3]*e2[2]
    ny = e1[3]*e2[1] - e1[1]*e2[3]
    nz = e1[1]*e2[2] - e1[2]*e2[1]
    area = 0.5 * sqrt(nx*nx + ny*ny + nz*nz)
    char_len_sq = area  # √A squared = A
    dx = target[1] - control_point[1]
    dy = target[2] - control_point[2]
    dz = target[3] - control_point[3]
    return (dx*dx + dy*dy + dz*dz) < (SELF_PAIR_EPS_REL * SELF_PAIR_EPS_REL) * char_len_sq
end

@inline function _panel_normal_gt(vertices)
    v1, v2, v3 = vertices
    e1 = v2 - v1
    e2 = v3 - v1
    nx = e1[2]*e2[3] - e1[3]*e2[2]
    ny = e1[3]*e2[1] - e1[1]*e2[3]
    nz = e1[1]*e2[2] - e1[2]*e2[1]
    inv_n = 1.0 / sqrt(nx*nx + ny*ny + nz*nz)
    return FastMultipole.StaticArrays.SVector{3,Float64}(nx*inv_n, ny*inv_n, nz*inv_n)
end

"""
    induced(target, source_system, source_buffer, i_source, derivatives_switch=FastMultipole.DerivativesSwitch(false, true, false); core_size=1.0e-3)
    induced(target, source_system, i_source, derivatives_switch=FastMultipole.DerivativesSwitch(false, true, false); core_size=1.0e-3)

Evaluate the panel-induced potential, velocity, and optional gradient at
`target` for source panel `i_source`.
"""
function induced(target::AbstractVector{TF}, source_system::AbstractBody{TK,NK,<:Any}, source_buffer::Matrix, i_source, derivatives_switch=FastMultipole.DerivativesSwitch(false,true,false),
        fam::Val=Val(FILAMENT_REGULARIZATION[]); core_size=1.0e-3) where {TF,TK,NK}

    # get vertices, rotation matrix
    R, v1, v2, v3 = rotate_to_panel(source_system, source_buffer, i_source)

    # get control point and strength
    control_point = FastMultipole.get_position(source_buffer, i_source)
    # strength = FastMultipole.get_strength(source_buffer, source_system, i_source)
    strength = FastMultipole.StaticArrays.SVector{NK,TF}(view(source_buffer, 5:4+NK, i_source))

    # evaluate influence
    potential, velocity, velocity_gradient = _induced(target, (v1, v2, v3), control_point, strength, TK, core_size, R, derivatives_switch, fam)

    # self-pair short-circuit: exterior velocity limit and interior potential limit.
    if _is_self_pair(target, control_point, (v1, v2, v3))
        n_gt = _panel_normal_gt((v1, v2, v3))
        potential, velocity, velocity_gradient = _self_limit(potential, velocity, velocity_gradient, TK, strength, n_gt)
    end

    # check for wake (if any)
    p, v, vg = _induced_wake(target, (v1, v2, v3), source_system, source_buffer, i_source, derivatives_switch)

    return potential+p, velocity+v, velocity_gradient+vg
end

function induced(target::AbstractVector{TF}, source_system::AbstractBody{TK,NK,<:Any}, i_source::Int, derivatives_switch=FastMultipole.DerivativesSwitch(false,true,false),
        fam::Val=Val(FILAMENT_REGULARIZATION[]); core_size=1.0e-3) where {TF,TK,NK}

    # get vertices, rotation matrix
    R, v1, v2, v3 = rotate_to_panel(source_system, i_source)

    # get control point and strength
    control_point = (v1 + v2 + v3) * 0.3333333333333333
    strength = FastMultipole.StaticArrays.SVector{NK,TF}(view(source_system.strength, i_source, :))

    # evaluate influence
    potential, velocity, velocity_gradient = _induced(target, (v1, v2, v3), control_point, strength, TK, core_size, R, derivatives_switch, fam)

    # self-pair short-circuit: exterior velocity limit and interior potential limit.
    if _is_self_pair(target, control_point, (v1, v2, v3))
        n_gt = _panel_normal_gt((v1, v2, v3))
        potential, velocity, velocity_gradient = _self_limit(potential, velocity, velocity_gradient, TK, strength, n_gt)
    end

    # check for wake (if any)
    p, v, vg = _induced_wake(target, (v1, v2, v3), source_system, i_source, derivatives_switch)

    # isnan(p) && println("Warning: NaN wake-induced potential at target $(target) from source panel $i_source with vertices $v1, $v2, $v3, core_size $core_size and strength $strength")

    return potential+p, velocity+v, velocity_gradient+vg
end

"Overload for non-rotated kernels"
function induced(target::AbstractVector{TF}, source_system::AbstractBody{VortexRing,NK,<:Any}, source_buffer::Matrix, i_source, derivatives_switch=FastMultipole.DerivativesSwitch(false,true,false),
        fam::Val=Val(FILAMENT_REGULARIZATION[]); core_size=1.0e-3) where {TF,NK}

    # get vertices
    v1, v2, v3 = get_vertices(source_system, source_buffer, i_source)

    # strength = FastMultipole.get_strength(source_buffer, source_system, i_source)
    strength = FastMultipole.StaticArrays.SVector{NK,TF}(view(source_buffer, 5:4+NK, i_source))

    # influence
    potential, velocity, velocity_gradient = _induced(target, (v1, v2, v3), strength, VortexRing, core_size, derivatives_switch, fam)

    # self-pair short-circuit
    control_point = (v1 + v2 + v3) * 0.3333333333333333
    if _is_self_pair(target, control_point, (v1, v2, v3))
        n_gt = _panel_normal_gt((v1, v2, v3))
        potential, velocity, velocity_gradient = _self_limit(potential, velocity, velocity_gradient, VortexRing, strength, n_gt)
    end

    # check for wake (if any)
    p, v, vg = _induced_wake(target, (v1, v2, v3), source_system, source_buffer, i_source, derivatives_switch)

    return potential+p, velocity+v, velocity_gradient+vg
end

function induced(target::AbstractVector{TF}, source_system::AbstractBody{VortexRing,NK,<:Any}, i_source::Int, derivatives_switch=FastMultipole.DerivativesSwitch(false,true,false),
        fam::Val=Val(FILAMENT_REGULARIZATION[]); core_size=1.0e-3) where {TF,NK}

    # get vertices
    v1, v2, v3 = get_vertices(source_system, i_source)

    strength = FastMultipole.StaticArrays.SVector{NK,TF}(view(source_system.strength, i_source, :))

    potential, velocity, velocity_gradient = _induced(target, (v1, v2, v3), strength, VortexRing, core_size, derivatives_switch, fam)

    # self-pair short-circuit
    control_point = (v1 + v2 + v3) * 0.3333333333333333
    if _is_self_pair(target, control_point, (v1, v2, v3))
        n_gt = _panel_normal_gt((v1, v2, v3))
        potential, velocity, velocity_gradient = _self_limit(potential, velocity, velocity_gradient, VortexRing, strength, n_gt)
    end

    # check for wake (if any)
    p, v, vg = _induced_wake(target, (v1, v2, v3), source_system, i_source, derivatives_switch)

    return potential+p, velocity+v, velocity_gradient+vg
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

# fam-discarding forwarder: the generic `induced` passes the filament
# regularization Val unconditionally; source/doublet kernels ignore it
_induced(target, vertices::NTuple, centroid::AbstractVector, strength,
    kernel::Union{Type{ConstantSource}, Type{ConstantDoublet}, Type{Union{ConstantSource, ConstantDoublet}}},
    core_radius, R, derivatives_switch::FastMultipole.DerivativesSwitch, ::Val) =
    _induced(target, vertices, centroid, strength, kernel, core_radius, R, derivatives_switch)


function compute_source_dipole(::FastMultipole.DerivativesSwitch{PS,VS,GS,NO,NM}, target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1, eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, strength::AbstractVector{TF}, ::Type{ConstantSource}, R_dot_s, reg_term) where {PS,VS,GS,NO,NM,TF}

    #--- compute values ---#

    # singularity if probing on a side
    num = max(eps(typeof(ri)), ri + rip1 - ds)
    log_term = log(num / (ri + rip1 + ds))

    # singularity at extension of the panel side; also force principal value
    # (tan_term = 0) when the target is on the panel plane (target_Rz == 0).
    # Self-pair (target == centroid in panel-local coords) ⇒ PV of solid-angle
    # integral is 0. Also force tan_term = 0 on the panel-side extension singularity.
    if target_Rz == zero(target_Rz) ||   # tRz snapped to an exact zero at roundoff scale in _induced (on-plane PV; keeps host/device/exact on the same ±2π side)
       abs(abs(R_dot_s) - ri * ds) <= 1e-12 * ri * ds  # (<=: ri*ds==0, target on a vertex, must trigger)
       # RELATIVE tol: R_dot_s and ri*ds carry units of length², so an
       # absolute window is geometry-dependent — it zeroed the tan_term of far targets within ~√(1e-12/(ri·ds)) rad of a
       # side extension, breaking the per-edge solid-angle cancellation for near-plane sliver pairs (spurious potential
       # ~3600× the truth on the DJI TE strip; found 2026-08-14 as the "p-saturated FMM floor", which was really this
       # direct-kernel defect — the multipole expansion was correct). On the extension line itself the limit IS zero.
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
        d2 = ds * ds
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

function compute_source_dipole(::FastMultipole.DerivativesSwitch{PS,VS,GS,NO,NM}, target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1, eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, strength::AbstractVector{TF}, ::Type{ConstantDoublet}, R_dot_s, reg_term) where {PS,VS,GS,NO,NM,TF}

    # singularity at extension of the panel side; also force principal value
    # (tan_term = 0) when the target is on the panel plane (target_Rz == 0).
    # Self-pair (target == centroid in panel-local coords) ⇒ PV of solid-angle
    # integral is 0. Also force tan_term = 0 on the panel-side extension singularity.
    if target_Rz == zero(target_Rz) ||   # tRz snapped to an exact zero at roundoff scale in _induced (on-plane PV; keeps host/device/exact on the same ±2π side)
       abs(abs(R_dot_s) - ri * ds) <= 1e-12 * ri * ds  # (<=: ri*ds==0, target on a vertex, must trigger)
       # RELATIVE tol: R_dot_s and ri*ds carry units of length², so an
       # absolute window is geometry-dependent — it zeroed the tan_term of far targets within ~√(1e-12/(ri·ds)) rad of a
       # side extension, breaking the per-edge solid-angle cancellation for near-plane sliver pairs (spurious potential
       # ~3600× the truth on the DJI TE strip; found 2026-08-14 as the "p-saturated FMM floor", which was really this
       # direct-kernel defect — the multipole expansion was correct). On the extension line itself the limit IS zero.
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

function compute_source_dipole(::FastMultipole.DerivativesSwitch{PS,VS,GS,NO,NM}, target_Rx, target_Ry, target_Rz, vx_i, vy_i, vx_ip1, vy_ip1, eip1, hip1, rip1, ei, hi, ri, ds, mi, dx, dy, strength::AbstractVector{TF}, ::Type{Union{ConstantSource, ConstantDoublet}}, R_dot_s, reg_term) where {PS,VS,GS,NO,NM,TF}

    # singularity if probing on a side
    # println("\nTESTING...")
    num = max(eps(typeof(ri)), ri + rip1 - ds)
    log_term = log(num / (ri + rip1 + ds))

    # singularity at extension of the panel side; also force principal value
    # (tan_term = 0) when the target is on the panel plane (target_Rz == 0).
    # Self-pair (target == centroid in panel-local coords) ⇒ PV of solid-angle
    # integral is 0. Also force tan_term = 0 on the panel-side extension singularity.
    if target_Rz == zero(target_Rz) ||   # tRz snapped to an exact zero at roundoff scale in _induced (on-plane PV; keeps host/device/exact on the same ±2π side)
       abs(abs(R_dot_s) - ri * ds) <= 1e-12 * ri * ds  # (<=: ri*ds==0, target on a vertex, must trigger)
       # RELATIVE tol: R_dot_s and ri*ds carry units of length², so an
       # absolute window is geometry-dependent — it zeroed the tan_term of far targets within ~√(1e-12/(ri·ds)) rad of a
       # side extension, breaking the per-edge solid-angle cancellation for near-plane sliver pairs (spurious potential
       # ~3600× the truth on the DJI TE strip; found 2026-08-14 as the "p-saturated FMM floor", which was really this
       # direct-kernel defect — the multipole expansion was correct). On the extension line itself the limit IS zero.
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

function _induced(target, vertices::NTuple{NS}, centroid::AbstractVector{TFP}, strength, kernel::Union{Type{ConstantSource}, Type{ConstantDoublet}, Type{Union{ConstantSource, ConstantDoublet}}}, core_radius, R, derivatives_switch::FastMultipole.DerivativesSwitch{PS,VS,GS,NO,NM}) where {TFP,NS,PS,VS,GS,NO,NM}
    #--- prelimilary computations ---#

    TFT = eltype(target)
    potential, velocity, velocity_gradient, target_Rx, target_Ry, target_Rz = source_dipole_preliminaries(TFT, TFP, target, centroid, R)
    # No on-centroid nudge: when target_Rz == 0 the kernel returns the principal
    # value (compute_source_dipole forces tan_term = 0 at target_Rz == 0), and
    # `_self_limit` applies the fixed exterior velocity and interior potential limits at self pairs.
    # On-plane snap: roundoff-scale target_Rz means the target IS on the panel
    # plane; snapping to an exact zero makes every edge take the PV branch of
    # the solid-angle tan_term, so the ±2π side cannot follow the sign of
    # arithmetic junk (host/device FMA divergence, job 13309929, 2026-08-22).
    # Mirrored in FastMultipole _rect_tri_source_doublet with the same 1e-24
    # relative-tolerance-squared vs the panel scale L² = Σᵢ|vᵢ - centroid|².
    # Plain-float types only: an AD dual must keep flowing through the smooth
    # atan branch or zero(target_Rz) erases its partials (∇φ ≡ 0 at on-plane
    # evaluation points, e.g. the control-point convergence tests).
    L2 = zero(promote_type(TFT, TFP))
    for v in vertices
        L2 += (v[1]-centroid[1])^2 + (v[2]-centroid[2])^2 + (v[3]-centroid[3])^2
    end
    target_Rz = _onplane_snap(target_Rz, L2)

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

function _induced(target, vertices::NTuple{NS}, strength, ::Type{VortexRing}, core_size, derivatives_switch::FastMultipole.DerivativesSwitch{PS,VS,GS,NO,NM},
        fam::Val=Val(FILAMENT_REGULARIZATION[])) where {NS,PS,VS,GS,NO,NM}
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
        p, _ = _induced(target, vertices, control_point, strength, ConstantDoublet, core_size, R, FastMultipole.DerivativesSwitch(true, false, false))
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
                v = _bound_vortex_velocity(r1, r2, finite_core, core_size, fam)
                velocity += v
            end
            if GS
                # velocity gradient
                g = _bound_vortex_gradient(r1, r2, finite_core, core_size, fam)
                gradient += g
            end
        end
    end

    return potential, velocity * strength[1], gradient * strength[1]
end

_induced(target, vertices, centroid, strength, kernel::Type{VortexRing}, core_size, R, derivatives_switch,
    fam::Val=Val(FILAMENT_REGULARIZATION[])) =
    _induced(target, vertices, strength, kernel, core_size, derivatives_switch, fam)


"""
    FilamentRegularization

Selectable regularization family for the bound-vortex filament kernel
(BRAINSTORM 025). The three point-blob-profile families share the numerator
`c*q` and differ only in the scalar denominator `D` (velocity) and
`∇D = κ ∇A` (gradient); derivations in
`BRAINSTORM/025_kernel_regularization_update/phase_01_theory.md`.
`LineGaussRegularization` is instead the exact line-convolved kernel — its
gradient does not fit the `(D, κ∇A)` shape and uses a dedicated cylindrical
assembly (052d DERIVATION.md §5).

- `GaussianRegularization` (default; Ryan ruling 2026-08-20, matched-CORE-
  SIZE convention): Lamb–Oseen transverse profile (`σ ≡ core_size`) — lowest
  peak velocity AND gradient at fixed core size. Radius inflation is
  GRADIENT-AWARE: `Δr = core_size·√(2z*)` with `z*` solving
  `e^(-z)(1+2z) = tol` (≈ `5.90·core_size` at tol 1e-6); the velocity-only
  radius `√(2 ln(1/tol))` would leave gradient error `tol·(1+2 ln(1/tol))`.
- `CompactRegularization`: the same compact-support family the
  doublet-velocity kernel uses (`regularize`) transplanted to the filament
  transverse profile — exactly singular (velocity AND gradient) beyond
  `core_size`, so `Δr = core_size`, tolerance-independent.
- `VatistasRegularization`: the legacy Vatistas n=2 core
  (`1/h² → 1/√(h⁴ + rc⁴)`); `Δr = core_size·(2/tol)^(1/4)` (velocity-derived,
  legacy-pinned; leaves gradient error ≤ 1.25·tol at that radius — see
  `radius_inflation`).
- `LineGaussRegularization`: exact closed form of the singular segment kernel
  convolved with the FLOWVPM Gaussian blob (`σ ≡ core_size`; FastMultipole
  `MATRIX_OPERATOR_REFACTOR/prototypes/052d_compact_kernel/DERIVATION.md`,
  2026-08-28). `GaussianRegularization` is exactly its infinite-line limit;
  the deviation from the singular kernel decays as `poly·e^(−d²/2σ²)` in the
  distance `d` to the SEGMENT — the Gaussian family's open along-line error
  channel closes by construction, so the radius inflation truly bounds the
  direct/expansion mismatch in the geometry the MAC measures. Inflation is
  the Gaussian gradient-aware fixed point plus a measured `0.35σ` pad
  (`≈ 6.25σ` at tol 1e-6). Costs 4 erf + 1 exp per edge (vs 1 expm1).

Select via [`set_filament_regularization!`](@ref) or the
`FLOWPANEL_FILAMENT_REG` environment variable (`compact`/`gaussian`/
`vatistas`/`linegauss`, read at package load) — the env hook exists so frozen drivers can
pin a family without code changes. The FMM stays aligned with the direct
kernel by construction: [`radius_inflation`](@ref) for `VortexRing` reads the
same global.
"""
@enum FilamentRegularization begin
    VatistasRegularization
    CompactRegularization
    GaussianRegularization
    LineGaussRegularization
end

"Active filament regularization family (see [`FilamentRegularization`](@ref))."
# Default: Gaussian (Ryan ruling 2026-08-20, matched-CORE-SIZE convention:
# at fixed core_size the Gaussian has the lowest peak velocity (0.45 vs
# compact 1.21, Vatistas 0.71, units Gamma/2pi/rc) AND lowest peak gradient
# (0.50 vs 2.55 / 1.00); its radius inflation rc*sqrt(2 ln 1/tol) ~ 5rc still
# removes the Vatistas 37.6rc pathology. See BRAINSTORM/025 phase_00 doc.
const FILAMENT_REGULARIZATION = Ref(GaussianRegularization)

"Set the active filament regularization family (type or Symbol
`:compact`/`:gaussian`/`:vatistas`/`:linegauss`)."
set_filament_regularization!(family::FilamentRegularization) =
    (FILAMENT_REGULARIZATION[] = family)
function set_filament_regularization!(family::Symbol)
    family === :compact && return set_filament_regularization!(CompactRegularization)
    family === :gaussian && return set_filament_regularization!(GaussianRegularization)
    family === :vatistas && return set_filament_regularization!(VatistasRegularization)
    family === :linegauss && return set_filament_regularization!(LineGaussRegularization)
    throw(ArgumentError("unknown filament regularization $(repr(family)); " *
        "use :compact, :gaussian, :vatistas, or :linegauss"))
end

#------- LineGauss closed-form kernel (052d, 2026-08-28) -------#
# Exact blob-line convolution of the segment Biot–Savart kernel with the
# FLOWVPM Gaussian core g(t) = erf(t/√2) − √(2/π)·t·e^(−t²/2). Derivation and
# validation live in FastMultipole `MATRIX_OPERATOR_REFACTOR/prototypes/
# 052d_compact_kernel/` (DERIVATION.md; k01/k03 harnesses). All lengths in
# the helpers are σ-scaled (σ ≡ core_size). Guard thresholds and series
# truncations are Float64-derived; a Float32 device port needs re-derivation.

const _LG_SQ2OPI = sqrt(2 / pi)
const _lg_erf = FLOWVPM.erf   # SpecialFunctions.erf via the FLOWVPM dep

# blob velocity function g(t) (odd in t), series-guarded as t → 0
@inline function _lg_gfun(t)
    at = abs(t)
    at >= 9.3 && return copysign(one(t), t)   # deviation < 2e-18
    if at < 0.125
        t2 = t * t
        term = t * t2 / 3
        s = term
        for m in 1:12
            term *= -t2 * (2m + 1) / (2m * (2m + 3))
            s += term
        end
        return _LG_SQ2OPI * s
    end
    return _lg_erf(t / sqrt(2)) - _LG_SQ2OPI * t * exp(-t * t / 2)
end

# on-axis antiderivative ψ (odd, ψ(0) = 0): M_axis = ψ(ẑ1) − ψ(ẑ2)
@inline function _lg_psi(z)
    z == 0 && return zero(z)
    if abs(z) < 0.125
        z2 = z * z
        return _LG_SQ2OPI * z * (1 / 3 - z2 / 30 + z2^2 / 280 -
                                 z2^3 / 3024 + z2^4 / 38016 - z2^5 / 549120)
    end
    return _LG_SQ2OPI * (z / 2) * exp(-z * z / 2) -
           _lg_gfun(z) / (2 * z * z) + _lg_gfun(z) / 2
end

# g(R)/R³ including its finite R = 0 limit (gradient axial factor)
@inline function _lg_kfun(R)
    if R < 0.125
        r2 = R * R
        return _LG_SQ2OPI * (1 / 3 - r2 / 10 + r2^2 / 56 -
                             r2^3 / 432 + r2^4 / 4224 - r2^5 / 49920)
    end
    return _lg_gfun(R) / R^3
end

# Fixed threshold at the error crossover: the axis limit truncates O(ĥ²)
# terms (rel err ≈ 0.25·ĥ² for a long mid-span segment) while the general
# branch loses ~eps/ĥ² to cancellation in N = ĥ²·M — both ≤ ~2.5e-8 at
# ĥ² = 1e-7. Scaling the threshold by min ẑ² (pre-2026-08-28 form) let the
# guard fire at physically large ĥ for long segments (ẑ ~ 4·10³ ⇒ ĥ² < 0.16),
# silently dropping the 1% O(ĥ²) correction at ĥ = 0.2.
@inline _lg_axis_guard(ĥ2, ẑ1, ẑ2) =
    ĥ2 < 1e-7 && ẑ1 != 0 && ẑ2 != 0

# wholly-small configuration: term-by-term integral of the convolution;
# returns M and the radial-gradient factor D = M + 2ĥ²·∂M/∂ĥ²
function _lg_small_radius_MD(ẑ1, ẑ2, ĥ2)
    M = zero(ĥ2)
    dM = zero(ĥ2)
    coeff = 1 / 3
    for m in 0:12
        Im = zero(ĥ2)
        dIm = zero(ĥ2)
        for k in 0:m
            p = m - k
            dzpow = (ẑ1^(2k + 1) - ẑ2^(2k + 1)) / (2k + 1)
            bc = binomial(m, k)
            Im += bc * ĥ2^p * dzpow
            p > 0 && (dIm += bc * p * ĥ2^(p - 1) * dzpow)
        end
        M += coeff * Im
        dM += coeff * dIm
        coeff *= -(2m + 3) / (2 * (m + 1) * (2m + 5))
    end
    M *= _LG_SQ2OPI
    return M, M + 2ĥ2 * _LG_SQ2OPI * dM
end

# one endpoint in the small-radius region: split the integral at |ẑ| = SMALL_R
const _LG_SMALL_R = 0.125
@inline _lg_endpoint_split_guard(ĥ2, ẑ1, ẑ2, R̂1, R̂2) =
    ĥ2 < 1e-8 * _LG_SMALL_R^2 && min(abs(ẑ1), abs(ẑ2)) < _LG_SMALL_R &&
    max(R̂1, R̂2) >= _LG_SMALL_R

function _lg_endpoint_split_MD(ẑ1, ẑ2, ĥ2)
    if abs(ẑ1) < _LG_SMALL_R
        split = -_LG_SMALL_R
        Mc, Dc = _lg_small_radius_MD(ẑ1, split, ĥ2)
        Mf = _lg_psi(split) - _lg_psi(ẑ2)
    else
        split = _LG_SMALL_R
        Mc, Dc = _lg_small_radius_MD(split, ẑ2, ĥ2)
        Mf = _lg_psi(ẑ1) - _lg_psi(split)
    end
    return Mf + Mc, Mf + Dc
end

# M = N/ĥ² with u = c·M/(4π σ² L); guarded near the axis and endpoints
function _lg_M(ẑ1, ẑ2, ĥ2, R̂1, R̂2)
    if ĥ2 == 0
        return _lg_psi(ẑ1) - _lg_psi(ẑ2)
    elseif max(R̂1, R̂2) < _LG_SMALL_R
        M, _ = _lg_small_radius_MD(ẑ1, ẑ2, ĥ2)
        return M
    elseif _lg_endpoint_split_guard(ĥ2, ẑ1, ẑ2, R̂1, R̂2)
        M, _ = _lg_endpoint_split_MD(ẑ1, ẑ2, ĥ2)
        return M
    elseif _lg_axis_guard(ĥ2, ẑ1, ẑ2)
        return _lg_psi(ẑ1) - _lg_psi(ẑ2)
    end
    # cancellation-reduced closed form: the endpoint Gaussians inside the
    # four g functions cancel exactly — 4 erf + 1 exp per edge
    G = exp(-ĥ2 / 2)
    N = ẑ1 * _lg_erf(R̂1 / sqrt(2)) / R̂1 -
        ẑ2 * _lg_erf(R̂2 / sqrt(2)) / R̂2 -
        G * (_lg_erf(ẑ1 / sqrt(2)) - _lg_erf(ẑ2 / sqrt(2)))
    return N / ĥ2
end

@inline _lg_skewmat(t::SVector{3,TF}) where TF = SMatrix{3,3,TF,9}(
    zero(TF), t[3], -t[2],
    -t[3], zero(TF), t[1],
    t[2], -t[1], zero(TF))

# LineGauss ∂u_i/∂x_j per unit Γ (cylindrical assembly, DERIVATION.md §5);
# same index convention as _bound_vortex_gradient (pinned by k01 T3c)
function _linegauss_gradient(r1::AbstractVector{TF}, r2, σ) where TF
    Z = zero(SMatrix{3,3,TF,9})
    s = r1 - r2
    B = dot(s, s)
    L = sqrt(B)
    L < 5 * eps(TF) && return Z
    that = -s / L
    ẑ1 = -dot(that, r1) / σ
    ẑ2 = ẑ1 - L / σ
    c = cross(r1, r2)
    ĥ2 = dot(c, c) / (B * σ * σ)
    R̂1 = norm(r1) / σ
    R̂2 = norm(r2) / σ
    M = _lg_M(ẑ1, ẑ2, ĥ2, R̂1, R̂2)
    C = 1 / (4 * pi * σ * σ)
    if ĥ2 == 0
        # on the segment axis: the regularized transverse derivative is finite
        return (C * M) * _lg_skewmat(that)
    end
    ĥ = sqrt(ĥ2)
    hvec = -r1 - (σ * ẑ1) * that      # h n̂ = (x − P1) − z1 t̂
    nh = norm(hvec)
    if nh <= 1e-10 * σ * max(R̂1, R̂2)
        # transverse direction lost to projection roundoff (can be exactly
        # zero → NaN): collapse to the deterministic axis-limit skew form
        return (C * M) * _lg_skewmat(that)
    end
    n̂ = hvec / nh
    b̂ = cross(that, n̂)
    k1 = _lg_kfun(R̂1)
    k2 = _lg_kfun(R̂2)
    duθdz = C * ĥ * (k1 - k2)
    if max(R̂1, R̂2) < _LG_SMALL_R
        _, radial = _lg_small_radius_MD(ẑ1, ẑ2, ĥ2)
        duθdh = C * radial
    elseif _lg_endpoint_split_guard(ĥ2, ẑ1, ẑ2, R̂1, R̂2)
        _, radial = _lg_endpoint_split_MD(ẑ1, ẑ2, ĥ2)
        duθdh = C * radial
    elseif _lg_axis_guard(ĥ2, ẑ1, ẑ2)
        duθdh = C * M          # bracket → 2M on the axis, so brk − M → M
    else
        G = exp(-ĥ2 / 2)
        brk = -ẑ1 * k1 + ẑ2 * k2 +
              G * (_lg_erf(ẑ1 / sqrt(2)) - _lg_erf(ẑ2 / sqrt(2)))
        duθdh = C * (brk - M)
    end
    uθ_h = C * M
    return duθdh * (b̂ * n̂') + duθdz * (b̂ * that') - uθ_h * (n̂ * b̂')
end

# Performance contract (BRAINSTORM 025 regression fix, 2026-08-20): the hot
# direct! loops must NEVER read FILAMENT_REGULARIZATION[] per edge — the
# non-const Ref load + 3-way branch inside the innermost kernels measured
# +34-49% on the production body influence pass (65.0 s vs 43.5-48.5 s,
# cluster A/B). The family is read ONCE per direct!-level call and crossed
# through a function barrier as `Val(family)`, so these `::Val{F}` methods
# compile with zero runtime family branches. The Val-less methods below are
# thin Ref-reading fallbacks for cold call sites (tests, probes) only.
@inline function _bound_vortex_velocity(r1::SVector{3,TF}, r2::SVector{3,TF}, finite_core, core_size,
        ::Val{F}) where {TF, F}
    # regularized filament kernel u = c*q/(4π D); D per the active
    # FilamentRegularization family (phase_01_theory.md)
    nr1 = norm(r1)
    nr2 = norm(r2)

    # target coincides with a filament endpoint: induced velocity is zero in the limit
    if nr1 < 5*eps(TF) || nr2 < 5*eps(TF)
        return zero(r1)
    end

    num = cross(r1, r2)
    r0 = r1 - r2
    dotrixrj = dot(num, num)        # A = |r1×r2|^2
    r0sqr = dot(r0, r0)             # B = |r0|^2
    rijdothat = dot(r0, r1/nr1 - r2/nr2)

    if !finite_core
        # singular kernel (no regularization)
        return num * rijdothat / dotrixrj / (4*pi)
    end

    if F === VatistasRegularization
        # 1/h^2 → 1/sqrt(h^4 + rc^4)
        V = num * rijdothat / sqrt(dotrixrj*dotrixrj + core_size*core_size*core_size*core_size * r0sqr*r0sqr) / (4*pi)
    elseif F === CompactRegularization
        # 1/h^2 → 1/(h^2 + δ(h)), δ = (h-rc)^2 inside the support, 0 beyond;
        # D = A + δB is exactly A for h ≥ rc and Brc² > 0 at h = 0 (no guard needed)
        h = sqrt(dotrixrj / r0sqr)
        D = h < core_size ? dotrixrj + (h - core_size)*(h - core_size) * r0sqr : dotrixrj
        V = num * rijdothat / D / (4*pi)
    elseif F === LineGaussRegularization
        # u = c·M/(4π σ² L): exact blob-line convolution (LineGauss section
        # above); matches the singular kernel to tol beyond ~6σ of the SEGMENT
        L = sqrt(r0sqr)
        σ = core_size
        ẑ1 = dot(r0, r1) / (L * σ)
        M = _lg_M(ẑ1, ẑ1 - L / σ, dotrixrj / (r0sqr * σ * σ), nr1 / σ, nr2 / σ)
        V = num * (M / (4*pi * σ * σ * L))
    else # GaussianRegularization
        # u = c*q*g(h)/(4π A), g = 1 - exp(-h²/2rc²); evaluated as
        # g/A = (g/x²)/(B rc²) with x² = (h/rc)² so the h → 0 limit is exact
        x2 = dotrixrj / (r0sqr * core_size * core_size)
        gscaled = x2 < 1e-12 ? TF(0.5) : TF(-expm1(-x2/2) / x2)
        V = num * rijdothat * gscaled / (r0sqr * core_size * core_size) / (4*pi)
    end

    return V
end

# cold-path fallback: one Ref read + one dynamic dispatch per CALL (not per
# edge inside a hot loop) — hot paths pass Val explicitly via the barriers
_bound_vortex_velocity(r1::SVector{3,<:Any}, r2::SVector{3,<:Any}, finite_core, core_size) =
    _bound_vortex_velocity(r1, r2, finite_core, core_size, Val(FILAMENT_REGULARIZATION[]))

@inline function _bound_vortex_gradient(r1::AbstractVector{TF}, r2, finite_core, core_size,
        ::Val{F}) where {TF, F}
    nr1 = norm(r1)
    nr2 = norm(r2)

    # target coincides with a filament endpoint: induced gradient is zero in the limit
    if nr1 < 5*eps(TF) || nr2 < 5*eps(TF)
        return zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
    end

    # LineGauss does not fit the (D, κ∇A) shape — dedicated cylindrical assembly
    if finite_core && F === LineGaussRegularization
        return _linegauss_gradient(r1, r2, core_size)
    end

    c = cross(r1, r2)
    s = r1 - r2
    A = dot(c, c)
    B = dot(s, s)
    q = dot(s, r1/nr1 - r2/nr2)

    # per-family D and ∇D = κ ∇A with ∇A = 2 s×c (B is target-independent);
    # see phase_01_theory.md for the κ derivations
    if finite_core
        if F === VatistasRegularization
            rc4 = core_size * core_size * core_size * core_size
            D = sqrt(A*A + rc4 * B*B)
            D == zero(D) && return zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
            κ = A / D
        elseif F === CompactRegularization
            B == zero(B) && return zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
            h = sqrt(A / B)
            if h < core_size
                D = A + (h - core_size)*(h - core_size) * B
                # κ = 2 - rc/h; the h → 0 clamp keeps κ finite where ∇A → 0
                # anyway (κ∇A has a finite limit; measure-zero perturbation)
                κ = 2 - core_size / max(h, eps(TF)*core_size)
            else
                D = A
                κ = one(TF)
            end
            D == zero(D) && return zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
        else # GaussianRegularization
            B == zero(B) && return zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
            x2 = A / (B * core_size * core_size)     # (h/rc)^2
            if x2 < 1e-12
                # series limits: D → 2 B rc², κ → 1/2
                D = 2 * B * core_size * core_size
                κ = TF(0.5)
            else
                g = -expm1(-x2/2)
                D = A / g
                κ = (1 - x2 * exp(-x2/2) / (2*g)) / g
            end
        end
        dD_coeff = κ * (2 * cross(s, c))
    else
        D = A
        D == zero(D) && return zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
        dD_coeff = 2 * cross(s, c)
    end

    r1hat = r1 / nr1
    r2hat = r2 / nr2
    dq_coeff = -((I - r1hat * transpose(r1hat)) * s) / nr1 +
                ((I - r2hat * transpose(r2hat)) * s) / nr2

    dc_dx = FastMultipole.StaticArrays.SMatrix{3,3,TF,9}(
        zero(TF), -s[3],    s[2],
           s[3], zero(TF), -s[1],
          -s[2],    s[1], zero(TF),
    )
    df_coeff = dq_coeff / D - q * dD_coeff / (D*D)

    return ONE_OVER_4PI * (dc_dx * (q / D) + c * transpose(df_coeff))
end

# cold-path fallback (see _bound_vortex_velocity note)
_bound_vortex_gradient(r1::AbstractVector, r2, finite_core, core_size) =
    _bound_vortex_gradient(r1, r2, finite_core, core_size, Val(FILAMENT_REGULARIZATION[]))

function _induced(target, vertices::NTuple, centroid::AbstractVector, strength, kernel::Type{Union{ConstantSource, VortexRing}}, core_radius, R, derivatives_switch::FastMultipole.DerivativesSwitch,
        fam::Val=Val(FILAMENT_REGULARIZATION[]))

    # source influence
    p, v, vg = _induced(target, vertices, centroid, SVector{1}(strength[1]), ConstantSource, core_radius, R, derivatives_switch)

    # vortex ring
    p_vr, v_vr, vg_vr = _induced(target, vertices, SVector{1}(strength[2]), VortexRing, core_radius, derivatives_switch, fam)

    return p + p_vr, v + v_vr, vg + vg_vr
end

#------- regularization functions -------#

@inline function regularize(distance, core_size)
    δ = distance < core_size ? (distance-core_size) * (distance-core_size) : zero(distance)

    return δ
end

"""
    radius_inflation(kernel, core_size, tol)

Distance beyond which the offset-regularized kernel matches the singular kernel
within relative tolerance `tol`. Added to the geometric panel radius written
into the FastMultipole source buffer (`source_system_to_buffer!`), so the
multipole-acceptance criterion only admits expansions — which represent the
*unregularized* kernel — where they agree with the regularized direct kernel to
`tol`. Without this term the direct/FMM operator mismatch saturates with
expansion order (021 Phase 1 finding, 2026-08-13).

- Source/doublet kernels: [`regularize`](@ref) is compactly supported — the
  regularized kernel is exactly singular beyond `core_size` — so the
  inflation is `core_size`, independent of `tol`.
- `VortexRing`: per the active [`FilamentRegularization`](@ref) family
  (BRAINSTORM 025; derivations in
  `BRAINSTORM/025_kernel_regularization_update/phase_01_theory.md`):
  - Gaussian (default): GRADIENT-AWARE radius. The velocity relative error is
    `e^(-z)` with `z = h²/2rc²`, but the gradient relative error is
    `e^(-z)(1+2z)`, so the velocity-derived radius `rc·√(2 ln(1/tol))` leaves
    gradient error `tol·(1+2 ln(1/tol))` (28.6× at 1e-6). The inflation
    solves `e^(-z*)(1+2z*) = tol` by the fixed point `z ← ln((1+2z)/tol)`
    (contraction rate `2/(1+2z)` ≈ 0.06; 5 iterations from `z₀ = ln(1/tol)`)
    ⇒ `Δr = rc·√(2z*)` ≈ `4.99rc / 5.47rc / 5.90rc` at tol `1e-4/1e-5/1e-6`.
  - compact-support: exactly singular — velocity AND gradient — beyond `rc`
    ⇒ inflation `rc`, independent of `tol` (matches the source/doublet rule).
  - LineGauss: the Gaussian fixed-point radius plus a `0.35rc` pad, calibrated
    against the measured segment-distance matching radii `5.25/5.75/6.25rc`
    at tol `1e-4/1e-5/1e-6` (052d k01 T7/T7b dense L/direction scans; the
    bare fixed point `4.99/5.47/5.90rc` is slightly non-conservative for the
    finite segment's polynomial prefactor). Semantic upgrade: for LineGauss
    the deviation from the singular kernel is bounded by distance to the
    SEGMENT — exactly the geometry the multipole-acceptance sphere measures —
    so the inflated radius is a true bound (the Gaussian family's line-
    distance `h` caveat does not apply).
  - Vatistas n=2 (`1/h² → 1/√(h⁴+rc⁴)`, legacy): velocity relative error
    ≈ ½(rc/h)⁴ ⇒ `rc·(2/tol)^(1/4)`. This shipped rule is velocity-derived
    and pinned by legacy-reproduction tests — the gradient relative error
    coefficient is 2.5(rc/h)⁴, so at the shipped radius the gradient error is
    ≤ 2.5·(tol/2) = 1.25·tol, absorbed by the multipole-acceptance margin
    (clearance ≥ Δ(1/MAC − 1) beyond the summed radii; thin for MAC > 0.5).

`tol = Inf` disables the inflation (pre-2026-08-13 behavior, for A/B runs).
"""
@inline radius_inflation(::Type{ConstantSource}, core_size, tol) =
    isinf(tol) ? zero(core_size) : core_size
@inline radius_inflation(::Type{ConstantDoublet}, core_size, tol) =
    isinf(tol) ? zero(core_size) : core_size
@inline radius_inflation(::Type{Union{ConstantSource, ConstantDoublet}}, core_size, tol) =
    isinf(tol) ? zero(core_size) : core_size
@inline function radius_inflation(::Type{VortexRing}, core_size, tol)
    isinf(tol) && return zero(core_size)
    family = FILAMENT_REGULARIZATION[]
    family == CompactRegularization && return core_size * one(tol)
    if family == GaussianRegularization || family == LineGaussRegularization
        # gradient-aware: solve e^(-z)(1+2z) = tol (see docstring)
        z = log(1 / tol)
        for _ in 1:5
            z = log((1 + 2z) / tol)
        end
        # LineGauss: +0.35rc pad calibrated against the measured SEGMENT-
        # distance matching radii 5.25/5.75/6.25rc at tol 1e-4/1e-5/1e-6
        # (052d k01 T7/T7b dense scans) vs the bare fixed point
        # 4.99/5.47/5.90rc; unlike the Gaussian family the bound is by
        # segment distance — exactly what the MAC sphere geometry measures
        family == LineGaussRegularization &&
            return core_size * (sqrt(2z) + oftype(sqrt(2z), 0.35))
        return core_size * sqrt(2z)
    end
    return core_size * (2 / tol)^0.25    # Vatistas (legacy, velocity-derived)
end
@inline radius_inflation(::Type{Union{ConstantSource, VortexRing}}, core_size, tol) =
    radius_inflation(VortexRing, core_size, tol)

#------- semi-infinite panels -------#

function get_strength_doublet(source_system::AbstractBody{Union{ConstantSource, ConstantDoublet}, 2, <:Any}, source_buffer, i_source)
    # get the strength of the doublet
    return source_buffer[6, i_source]
end

function get_strength_doublet(source_system::AbstractBody{Union{ConstantSource, VortexRing}, 2, <:Any}, source_buffer, i_source)
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

function _induced_wake(target::AbstractVector{TF}, vertices::Tuple, source_system::RigidWakeBody, i_source::Int, derivatives_switch::FastMultipole.DerivativesSwitch{PS,VS,GS,NO,NM}) where {TF,PS,VS,GS,NO,NM}
    # body-only operator assembly (e.g. the Green-system B matrix) skips the
    # attached wake entirely
    if source_system.suppress_attached_wake[]
        return zero(TF), zero(FastMultipole.StaticArrays.SVector{3,TF}), zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
    end

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

        # get strength; in physical mode the affine attached-wake correction
        # shifts the transition-panel strength (∓c/2 across a paired edge, so
        # the net attached circulation is γ = μ_u − μ_l − c)
        strength = get_strength_doublet(source_system, i_source)
        if source_system.wake_correction_active[]
            strength += source_system.wake_strength_shift[i_source]
        end
        TK = get_wake_kernel(source_system)

        # evaluate potential
        if source_system.semiinfinite_wake
            return induced_semiinfinite(target, TK, v1x, v1y, v1z, v2x, v2y, v2z, Dax, Day, Daz, strength, derivatives_switch; core_size=source_system.core_size)
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
            p, v, g = _induced(target, (v1, v2, vw1), control_point, strength_vec, TK, source_system.core_size, R, derivatives_switch)

            # wake node connected to the second vertex (TE2)
            Dbx, Dby, Dbz = source_system.Das[i_surf][1, das_col_2], source_system.Das[i_surf][2, das_col_2], source_system.Das[i_surf][3, das_col_2]
            v2w_x = v2x + Dbx
            v2w_y = v2y + Dby
            v2w_z = v2z + Dbz
            vw2 = FastMultipole.StaticArrays.SVector{3,TF}(v2w_x, v2w_y, v2w_z)
            
            # influence of the second triangle
            control_point = (vw1 + v2 + vw2) * 0.333333333333333
            R, _ = rotate_to_panel(v1w_x, v1w_y, v1w_z, v2x, v2y, v2z, v2w_x, v2w_y, v2w_z)
            dp, dv, dg = _induced(target, (vw1, v2, vw2), control_point, strength_vec, TK, source_system.core_size, R, derivatives_switch)
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

function _induced_wake(target::AbstractVector{TF}, vertices::Tuple, source_system::RigidWakeBody{<:Any,NK,<:Any}, source_buffer::Matrix, i_source::Int, derivatives_switch::FastMultipole.DerivativesSwitch{PS,VS,GS,NO,NM}) where {TF,NK,PS,VS,GS,NO,NM}
    # body-only operator mode (e.g. Green-system B products): no attached wake
    if source_system.suppress_attached_wake[]
        return zero(TF), zero(FastMultipole.StaticArrays.SVector{3,TF}), zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
    end

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
            return induced_semiinfinite(target, TK, v1x, v1y, v1z, v2x, v2y, v2z, Dax, Day, Daz, strength, derivatives_switch; core_size=source_system.core_size)
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
            return _induced_quad(target, (v1, v2, vw2, vw1), strength_vec, TK, source_system.core_size, derivatives_switch)

        end
    else
        return zero(TF), zero(FastMultipole.StaticArrays.SVector{3,TF}), zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
    end
end

function _induced_quad(target, vertices, strength, kernel::Type{ConstantDoublet}, core_size, derivatives_switch)
    # influence of first triangle
    v1 = vertices[1]
    v2 = vertices[2]
    vw1 = vertices[4]
    control_point = (v1 + v2 + vw1) * 0.333333333333333
    R, _ = rotate_to_panel(v1[1], v1[2], v1[3], v2[1], v2[2], v2[3], vw1[1], vw1[2], vw1[3])
    p, vel, g = _induced(target, (v1, v2, vw1), control_point, strength, kernel, core_size, R, derivatives_switch)

    # influence of the second triangle
    vw2 = vertices[3]
    control_point = (vw1 + v2 + vw2) * 0.333333333333333
    R, _ = rotate_to_panel(vw1[1], vw1[2], vw1[3], v2[1], v2[2], v2[3], vw2[1], vw2[2], vw2[3])
    dp, dvel, dg = _induced(target, (vw1, v2, vw2), control_point, strength, kernel, core_size, R, derivatives_switch)
    
    return p+dp, vel+dvel, g+dg
end

function _induced_quad(target, vertices, strength, kernel::Type{VortexRing}, core_size, derivatives_switch::FastMultipole.DerivativesSwitch{PS,VS,GS,NO,NM}) where {PS,VS,GS,NO,NM}
    if PS
        # influence of first triangle
        v1 = vertices[1]
        v2 = vertices[2]
        vw1 = vertices[4]
        control_point = (v1 + v2 + vw1) * 0.333333333333333
        R, _ = rotate_to_panel(v1[1], v1[2], v1[3], v2[1], v2[2], v2[3], vw1[1], vw1[2], vw1[3])
        p, vel, g = _induced(target, (v1, v2, vw1), control_point, strength, kernel, core_size, R, derivatives_switch)

        # influence of the second triangle
        vw2 = vertices[3]
        control_point = (vw1 + v2 + vw2) * 0.333333333333333
        R, _ = rotate_to_panel(vw1[1], vw1[2], vw1[3], v2[1], v2[2], v2[3], vw2[1], vw2[2], vw2[3])
        dp, dvel, dg = _induced(target, (vw1, v2, vw2), control_point, strength, kernel, core_size, R, derivatives_switch)

        return p+dp, vel+dvel, g+dg
    else
        return _induced(target, vertices, strength, kernel, core_size, derivatives_switch)
    end
end

function induced_semiinfinite(target::AbstractVector, TK::Type{VortexRing}, args...; core_size)
    return induced_semiinfinite(target, ConstantDoublet, args...; core_size)
end

function induced_semiinfinite(target::AbstractVector{TF}, TK::Type{ConstantDoublet}, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, strength, ::FastMultipole.DerivativesSwitch{PS,VS,GS,NO,NM}; core_size) where {TF,PS,VS,GS,NO,NM}
    potential = zero(TF)
    velocity = zero(FastMultipole.StaticArrays.SVector{3,TF})
    gradient = zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
    if PS
        potential += _phi_semiinfinite(target, TK, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, strength; core_size)
    end
    if VS
        velocity += _U_semiinfinite(target, TK, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, strength; core_size)
    end
    if GS
        gradient += _U_semiinfinite_gradient(target, TK, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, strength; core_size)
    end

    return potential, velocity, gradient
end

function _phi_semiinfinite(target::AbstractVector{TF}, TK::Type{ConstantDoublet}, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, strength::Number; core_size=1e-8) where TF

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
        derivatives_switch = FastMultipole.DerivativesSwitch(true, false, false)

        # compute potential
        potential, _ = _induced(target, (v1, v2, v3), control_point, this_strength, TK, core_size, R, derivatives_switch)
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

function _U_semiinfinite(target::AbstractVector{TF}, ::Type{ConstantDoublet}, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, strength::Number; core_size=1e-8) where TF
    Ux1, Uy1, Uz1 = _U_semiinfinite_vortex(v1x, v1y, v1z,
                                    d1, d2, d3,
                                    -strength,
                                    target; offset=core_size)
    
    Ux2, Uy2, Uz2 = _U_semiinfinite_vortex(v2x, v2y, v2z,
                                    d1, d2, d3,
                                    strength,
                                    target; offset=core_size)
    
    Uxb, Uyb, Uzb = _U_boundvortex(v1x, v1y, v1z,
                                    v2x, v2y, v2z,
                                    strength,
                                    target; offset=core_size)

    # Uxb, Uyb, Uzb = _bound_vortex_velocity(
    #                     FastMultipole.StaticArrays.SVector{3,TF}(v1x, v1y, v1z),
    #                     FastMultipole.StaticArrays.SVector{3,TF}(v2x, v2y, v2z),
    #                     true,
    #                     core_size
    #                 ) .* strength
    
    # combine contributions
    Ux = Ux1 + Ux2 + Uxb
    Uy = Uy1 + Uy2 + Uyb
    Uz = Uz1 + Uz2 + Uzb

    return FastMultipole.StaticArrays.SVector{3,TF}(-Ux, -Uy, -Uz)
end

function _U_semiinfinite_gradient(target::AbstractVector{TF}, ::Type{ConstantDoublet}, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, strength::Number; core_size=1e-8) where TF
    g1 = _U_semiinfinite_vortex_gradient(v1x, v1y, v1z,
                                    d1, d2, d3,
                                    -strength,
                                    target; offset=core_size)
    
    g2 = _U_semiinfinite_vortex_gradient(v2x, v2y, v2z,
                                    d1, d2, d3,
                                    strength,
                                    target; offset=core_size)

    gb = _U_boundvortex_gradient(v1x, v1y, v1z,
                                    v2x, v2y, v2z,
                                    strength,
                                    target; offset=core_size)

    return -(g1 + g2 + gb)
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

function _U_boundvortex_gradient(pa1::TF1, pa2::Number, pa3::Number,
                        pb1::Number, pb2::Number, pb3::Number,
                        strength::TF2,
                        target::AbstractVector{TF3};
                        offset=1e-8,
                       ) where{TF1, TF2, TF3}
    TF = promote_type(TF1, TF2, TF3)
    pa = FastMultipole.StaticArrays.SVector{3,TF}(pa1, pa2, pa3)
    pb = FastMultipole.StaticArrays.SVector{3,TF}(pb1, pb2, pb3)
    x = FastMultipole.StaticArrays.SVector{3,TF}(target[1], target[2], target[3])
    return -strength * _bound_vortex_gradient(pa - x, pb - x, true, offset)
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

function _U_semiinfinite_vortex_gradient(p1::TF1, p2::Number, p3::Number,
                                d1::Number, d2::Number, d3::Number,
                                strength::TF2,
                                target::AbstractVector{TF3};
                                cutoff=1e-14, offset=1e-8) where {TF1,TF2,TF3}
    TF = promote_type(TF1, TF2, TF3)
    scaling_factor = 1/(4*π)

    if abs(d1*d1 + d2*d2 + d3*d3 - 1) > 2*eps()
        error("Found non-unitary semi-infinite direction"*
                " norm([d1, d2, d3]) = norm($([d1,d2,d3])) = $(norm((d1,d2,d3)))")
    end

    p = FastMultipole.StaticArrays.SVector{3,TF}(p1, p2, p3)
    d = FastMultipole.StaticArrays.SVector{3,TF}(d1, d2, d3)
    x = FastMultipole.StaticArrays.SVector{3,TF}(target[1], target[2], target[3])
    y = x - p
    a = dot(y, d)
    h = y - a*d
    H = dot(h, h)

    if H <= cutoff^2
        return zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})
    end

    y_norm = norm(y)
    y_norm == zero(y_norm) && return zero(FastMultipole.StaticArrays.SMatrix{3,3,TF,9})

    a2 = a*a
    active_bound = a2 > offset^2 && a2*H > cutoff^2
    E = sqrt(H*H + offset*offset*offset*offset)
    n = cross(d, h)
    factor = one(TF) + (active_bound ? a / y_norm : zero(TF))
    f = factor / E

    grad_inv_E = -2 * H / (E*E*E) * h
    grad_factor = active_bound ? d / y_norm - a * y / (y_norm*y_norm*y_norm) :
                                 zero(FastMultipole.StaticArrays.SVector{3,TF})
    grad_f = grad_factor / E + factor * grad_inv_E

    dn_dx = FastMultipole.StaticArrays.SMatrix{3,3,TF,9}(
        zero(TF),   d[3], -d[2],
          -d[3], zero(TF),  d[1],
           d[2],  -d[1], zero(TF),
    )

    return -strength * scaling_factor * (dn_dx * f + n * transpose(grad_f))
end
