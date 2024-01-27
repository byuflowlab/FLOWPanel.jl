# @testset "solver" begin

# create ground panels
nx, dx, z = 10, 0.5, 0.0
expansion_order, n_per_branch, theta, shrink_recenter = 7, 10, 0.4, true
panels = initialize_panels(nx, dx, z; expansion_order_panels=expansion_order, n_per_branch_panels=n_per_branch, theta_panels=theta, shrink_recenter_panels=shrink_recenter, invert_normals=false)

# apply external velocity
Random.seed!(123)
for ipanel in eachindex(panels.velocity)
    panels.velocity[ipanel] = -rand(SVector{3,Float64})
    panels.velocity[ipanel] = SVector{3,Float64}(0,0,1.0)
end

# solve for panel strengths
solve_panels!(panels)

# apply panel-induced velocity
self_induced_direct!(panels)

# test for flow tangency
global v_normal = 0.0
for v in panels.velocity
    global  v_normal += v[3]
end
v_normal /= length(panels.velocity)
@test isapprox(v_normal, 0.0; atol=1e-10)

#####
##### inverted normals
#####
panels = initialize_panels(nx, dx, z; expansion_order_panels=expansion_order, n_per_branch_panels=n_per_branch, theta_panels=theta, shrink_recenter_panels=shrink_recenter, invert_normals=true)

# apply external velocity
Random.seed!(123)
for ipanel in eachindex(panels.velocity)
    panels.velocity[ipanel] = -rand(SVector{3,Float64})
    panels.velocity[ipanel] = SVector{3,Float64}(0,0,1.0)
end

# save external velocity
# vtk(panels; run_name="external_velocity")

# solve for panel strengths
solve_panels!(panels)

# apply panel-induced velocity
self_induced_direct!(panels)

# save net velocity
# vtk(panels; run_name="solved_velocity")

# test for flow tangency
global v_normal = 0.0
for v in panels.velocity
    global  v_normal += v[3]
end
v_normal /= length(panels.velocity)
@test isapprox(v_normal, 0.0; atol=1e-10)

# end

#=
@testset "fmm- single panel" begin

# create single panel
nx, dx, z = 1, 0.5, 0.0
expansion_order = 5
n_per_branch, theta, shrink_recenter = 1, 0.4, false
single_panel = initialize_panels(nx, dx, z; expansion_order_panels=expansion_order, n_per_branch_panels=n_per_branch, theta_panels=theta, shrink_recenter_panels=shrink_recenter, invert_normals=true)
strength = 1.0
single_panel.strengths[1] = strength

# multipole expansion
single_panel_wrapped = fmm.SortWrapper(single_panel)
tree = fmm.Tree(single_panel_wrapped)
fmm.upward_pass_single_thread!(tree.branches, single_panel_wrapped, expansion_order)

# probe potential/velocity/gradient
centroid = single_panel.centroid_grid[1]
normal = single_panel.normals[1]
vertices = get_vertices(single_panel, 1)
probe(target) = potential_velocity_gradient_source_panel(target, centroid, normal, strength, vertices...)

# direct evaluation
pots = []
vs = []
grads = []
scls = range(1.0, stop=50, length=100)
for scl in scls
    pot, vel, grad = probe(centroid + normal * scl)
    push!(pots, pot)
    push!(vs, vel)
    push!(grads, grad)
end

# expansion
function evaluate_expansion(target, i_branch, tree)
    target_potential = zeros(4)
    fmm.M2B!(target_potential, target, i_branch, tree)

    return target_potential
end

pots_expansion = []
for scl in scls
    target = centroid + normal * scl
    i_branch = 1
    pot = evaluate_expansion(target, i_branch, tree)[1]
    push!(pots_expansion, pot)
end

for i in 1:length(pots)
    @test isapprox(pots[i], pots_expansion[i]; atol=1e-6)
end

#####
##### test velocity/velocity gradient
#####


end

@testset "fmm- multiple panels" begin

# create ground panels
nx, dx, z = 10, 0.5, 0.0
expansion_order, n_per_branch, theta, shrink_recenter = 6, 1, 0.25, true
panels = initialize_panels(nx, dx, z; expansion_order_panels=expansion_order, n_per_branch_panels=n_per_branch, theta_panels=theta, shrink_recenter_panels=shrink_recenter, invert_normals=true)

# random strengths
for i in eachindex(panels.strengths)
    panels.strengths[i] = rand()
end

# direct interaction
self_induced_direct!(panels)
potential_direct = deepcopy(panels.potential)
velocity_direct = deepcopy(panels.velocity)
velocity_gradient_direct = deepcopy(panels.velocity_gradient)

# fmm interaction
reset_panels!(panels)
panels_wrapped = fmm.SortWrapper(panels)
self_induced_fmm!(panels_wrapped)
potential_fmm = deepcopy(panels.potential)
velocity_fmm = deepcopy(panels.velocity)
velocity_gradient_fmm = deepcopy(panels.velocity_gradient)

# test
for i in eachindex(potential_direct)
    @test isapprox(potential_direct[i], potential_fmm[i]; atol=1e-12)
    for j in 1:3
        @test isapprox(velocity_direct[i][j], velocity_fmm[i][j]; atol=1e-12)
    end
    for j in 1:9
        @test isapprox(velocity_gradient_direct[i][j], velocity_gradient_fmm[i][j]; atol=1e-12)
    end
end

end
=#