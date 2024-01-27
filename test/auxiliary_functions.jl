struct Quaternion{TF}
    real::TF
    pure::SVector{3,TF}
end

function rotate(v, q::Quaternion)
    t = 2*cross(q.pure, v)
    return v + q.real*t + cross(q.pure,t)
end

"perform the reverse rotation"
function antirotate(v, q::Quaternion)
    q = Quaternion(q.real, -q.pure)
    return rotate(v, q)
end

function get_quaternion(axis, angle)
    @assert isapprox(norm(axis), 1.0; atol=1e-12) "axis must be a unit vector"
    stheta_over_2, ctheta_over_2 = sincos(angle/2)
    real_part = ctheta_over_2
    pure_part = axis * stheta_over_2
    quat = Quaternion(real_part, pure_part)
    return quat
end

function rotate(v, axis, angle)
    return rotate(v, get_quaternion(axis, angle))
end

function antirotate(v, axis, angle)
    return rotate(v, get_quaternion(axis, -angle))
end

function create_panel(;
        vertices = SVector{4}(
            SVector{3}(-0.5,-0.5,0.0),
            SVector{3}(0.5,-0.5,0.0),
            SVector{3}(0.5,0.5,0.0),
            SVector{3}(-0.5,0.5,0.0),
        ),
        control_point = SVector{3}(0.0,0.0,0.0),
        normal = SVector{3}(0,0,1.0),
        strength = SVector{1}(1.0),
        invert_normals = false,
    )

    radius = FLOWPanel.get_radius(control_point, vertices...)
    
    if invert_normals
        vertices = reverse(vertices)
        normal = -1 * normal
    end

    return Panel(vertices, control_point, normal, strength, radius)
end

function test_velocity(target, panel, kernel=ConstantSource())
    this_potential(target) = induced(target, panel, kernel)[1]
    return -ForwardDiff.gradient(this_potential, target)
end

function test_gradient(target, panel, kernel=ConstantSource())
    this_velocity(target) = induced(target, panel, kernel)[2]
    return ForwardDiff.jacobian(this_velocity, target)
end
