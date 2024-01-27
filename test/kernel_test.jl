@testset "single constant source panel" begin

# create single panel
panel = create_panel()
(; vertices, control_point, normal, strength) = panel
strength = strength[1]

#####
##### potential function
#####

probe(target) = induced(target, panel, ConstantSource(); toggle_potential=true, toggle_velocity=true, toggle_hessian=true)

pots = []
vs = []
grads = []
scls = range(0, stop=50, length=100)
scls = [100.0]
probe_direction = SVector{3,Float64}(0.0,0.0,1.0)
for scl in scls
    pot, vel, grad = probe(control_point + probe_direction * scl)
    push!(pots, pot)
    push!(vs, vel)
    push!(grads, grad)
end

# what it should be
x, y, z = 0.0, 0.0, scls[end]
xks = [v[1] for v in vertices]
xkp1 = vcat(xks[2:end],xks[1])
yks = [v[2] for v in vertices]
ykp1 = vcat(yks[2:end],yks[1])
ds = [sqrt((x-xp1)^2+(y-yp1)^2) for (x,y,xp1,yp1) in zip(xks,yks,xkp1,ykp1)]
ms = [(ykp1[i]-yks[i])/(xkp1[i]-xks[i]) for i in 1:4]
rs = [sqrt(x^2 + y^2 + z^2) for (x,y) in zip(xks,yks)]
rp1s = vcat(rs[2:end],rs[1])
es = [x^2 + z^2 for x in xks]
ep1s = vcat(es[2:end],es[1])
hs = [-xks[i] * -yks[i] for i in 1:4]
hp1s = vcat(hs[2:end],hs[1])

phi_test = 0.0
phi_xx_test = 0.0
for i in 1:4
    f1 = (-xks[i] * (ykp1[i]-yks[i]) - (-yks[i]) * (xkp1[i] - xks[i])) / ds[i]
    log_term = log((rs[i] + rp1s[i] - ds[i]) / (rs[i] + rp1s[i] + ds[i]))
    atan_term = z * ( atan((ms[i]*es[i] - hs[i]) / (z*rs[i])) - atan((ms[i]*ep1s[i] - hp1s[i]) / (z*rp1s[i])) )
    phi_test += f1*log_term + atan_term
    # global phi_test
    # global phi_xx_test
    phi_xx_test += 2*(ykp1[i]-yks[i])/((rs[i]+rp1s[i])^2 - ds[i]^2) * ((x-xks[i])/rs[i] + (x-xkp1[i])/rp1s[i])
end
phi_test *= 1/4/pi
phi_xx_test *= -1/4/pi

@test isapprox(phi_test, 1/scls[end]/4/pi; atol=1e-6)
@test isapprox(pots[end], 1/scls[end]/4/pi; atol=1e-6)

# velocity is the negative gradient of the potential, or the negative partial derivative in the z direction
@test isapprox(vs[end][3], 1/scls[end]^2/4/pi; atol=1e-6)
x, y, z = probe_direction * scls[end]
x2y2z2_to_52 = (x^2+y^2+z^2)^(5/2)
test_hessian = -SMatrix{3,3,Float64,9}(
    (2*x^2 - (y^2+z^2)) / x2y2z2_to_52, 3*x*y/x2y2z2_to_52, 3*x*z/x2y2z2_to_52,
    3*x*y/x2y2z2_to_52, (2*y^2 - x^2 - z^2)/x2y2z2_to_52, 3*y*z/x2y2z2_to_52,
    3*x*z/x2y2z2_to_52, 3*y*z/x2y2z2_to_52, (2*z^2-x^2-y^2)/x2y2z2_to_52
)
for i in 1:9
    @test isapprox(grads[end][i], test_hessian[i]; atol=2e-5)
end

#####
##### potential function- inverted normal
#####
panel_inverted = create_panel(; invert_normals=true)

probe2(target) = induced(target, panel_inverted, ConstantSource())

pots = []
vs = []
grads = []
scls = range(1, stop=50, length=100)
scls = [50.0]
probe_direction = panel_inverted.normal
for scl in scls
    target = control_point + probe_direction * scl
    pot, vel, grad = probe2(target)
    push!(pots, pot)
    push!(vs, vel)
    push!(grads, grad)
end

@test isapprox(pots[end], 1/scls[end]/4/pi; atol=1e-6)
# velocity is the negative gradient of the potential, or the negative partial derivative in the z direction
@test isapprox(vs[end][3], -1/scls[end]^2/4/pi; atol=1e-6)
x, y, z = probe_direction * scls[end]
x2y2z2_to_52 = (x^2+y^2+z^2)^(5/2)
test_hessian = -SMatrix{3,3,Float64,9}(
    (2*x^2 - (y^2+z^2)) / x2y2z2_to_52, 3*x*y/x2y2z2_to_52, 3*x*z/x2y2z2_to_52,
    3*x*y/x2y2z2_to_52, (2*y^2 - x^2 - z^2)/x2y2z2_to_52, 3*y*z/x2y2z2_to_52,
    3*x*z/x2y2z2_to_52, 3*y*z/x2y2z2_to_52, (2*z^2-x^2-y^2)/x2y2z2_to_52
    )
for i in 1:9
    @test isapprox(grads[end][i], test_hessian[i]; atol=2e-5)
end
    
end


@testset "constant source velocity/gradient" begin

# simulation options
n_targets = 10
nondimensional_targets = range(0.0, stop=5, length=n_targets)

# set random seed
Random.seed!(123)

panel = create_panel(; invert_normals=false)

centroid = panel.control_point
normal = panel.normal
vertex_1, vertex_2, vertex_3, vertex_4 = panel.vertices
vertex_1 += SVector{3}(0.1 * (rand()-0.5), 0.1 * (rand()-0.5), 0.0)
vertex_2 += SVector{3}(0.1 * (rand()-0.5), 0.1 * (rand()-0.5), 0.0)
vertex_3 += SVector{3}(0.1 * (rand()-0.5), 0.1 * (rand()-0.5), 0.0)
vertex_4 += SVector{3}(0.1 * (rand()-0.5), 0.1 * (rand()-0.5), 0.0)

axis, angle = rand(SVector{3,Float64}), pi/6
axis /= norm(axis)
vertex_1 = rotate(vertex_1, axis, angle)
vertex_2 = rotate(vertex_2, axis, angle)
vertex_3 = rotate(vertex_3, axis, angle)
vertex_4 = rotate(vertex_4, axis, angle)
normal = rotate(normal, axis, angle)

new_panel = Panel(SVector{4}(vertex_1, vertex_2, vertex_3, vertex_4), centroid, normal, panel.strength, panel.radius)

check_normal = cross(vertex_2-vertex_1, vertex_4-vertex_1)
check_normal /= norm(check_normal)
@assert isapprox(check_normal, normal; atol=1e-12) "normal no longer matches"

# create targets
Random.seed!(123)
direction = rand(3)
direction = panel.normal
direction /= norm(direction)
direction = SVector{3,Float64}(direction)
targets = [centroid + direction * nondimensional_targets[i] for i in eachindex(nondimensional_targets)]

for target in targets[2:end]
    v = induced(target, new_panel, ConstantSource())[2]
    v_check = test_velocity(target, new_panel)
    for j in eachindex(v)
        @test isapprox(v, v_check; atol=1e-12)
    end
end

vgrads = []
vgradchecks = []
for target in targets[2:end]
    vgrad = induced(target, new_panel, ConstantSource())[3]
    vgrad_check = test_gradient(target, new_panel)
    push!(vgrads, vgrad)
    push!(vgradchecks, vgrad_check)
    for j in eachindex(vgrad)
        @test isapprox(vgrad[j], vgrad_check[j]; atol=1e-12)
    end
end

end

@testset "single constant doublet panel" begin

panel = create_panel()

target = rand(SVector{3}) * 2
vz_source = induced(target, panel, ConstantSource())[2][3]
psi_doublet, v_doublet, grad_doublet = induced(target, panel, ConstantNormalDoublet())

@test isapprox(vz_source, psi_doublet; atol=1e-10)

v_doublet_test = test_velocity(target, panel, ConstantNormalDoublet())

for i in eachindex(v_doublet)
    @test isapprox(v_doublet[i], v_doublet_test[i]; atol=1e-10)
end

grad_doublet_test = test_gradient(target, panel, ConstantNormalDoublet())

for i in eachindex(grad_doublet)
    @test isapprox(grad_doublet[i], grad_doublet_test[i]; atol=1e-10)
end

end

@testset "single constant source plus doublet panel" begin

panel = create_panel(; strength=SVector{2}(1.0,1.0))

target = rand(SVector{3}) * 2

psi_source, v_source, grad_source = induced(target, panel, ConstantSource())
psi_doublet, v_doublet, grad_doublet = induced(target, panel, ConstantNormalDoublet())
psi_combined = psi_source + psi_doublet
v_combined = v_source + v_doublet
grad_combined = grad_source + grad_doublet

psi_source_doublet, v_source_doublet, grad_source_doublet = induced(target, panel, ConstantSourceNormalDoublet())

@test isapprox(psi_combined, psi_source_doublet; atol=1e-10)

for i in eachindex(v_combined)
    @test isapprox(v_combined[i], v_source_doublet[i]; atol=1e-10)
end

for i in eachindex(grad_doublet)
    @test isapprox(grad_combined[i], grad_source_doublet[i]; atol=1e-10)
end

end