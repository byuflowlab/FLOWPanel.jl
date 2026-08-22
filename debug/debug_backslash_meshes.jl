using Test
import FLOWPanel as pnl
using LinearAlgebra
import Meshes
using StaticArrays: SVector

include(joinpath(@__DIR__, "..", "test", "test_helpers.jl"))

# ---------------- Same setup as test/runtests_unit_solver.jl:507 -------------
run_names = ["nasa_wing.msh", "nasa_surface_spaced_repaired.msh"]
files = [joinpath(pnl.examples_path, "wing_aileron", name) for name in run_names]
nodes1 = [42, 19]
nodes2 = [34, 3]

m              = 0.0254
magVinf        = 117.3 * m * 12

c_body1 = 10 * m
b       = 60 * m
c_body2 = 2.0
chords  = [c_body1, c_body2]
scaling = 1.0

kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}
bodytype = pnl.RigidWakeBody{kernel}

Vinf = magVinf * [1.0, 0.0, 0.0]

bodies = tuple([generate_body(file, chord, b, bodytype, scaling, 1, Vinf, firstnode, secondnode)
                for (file, chord, firstnode, secondnode) in zip(files, chords, nodes1, nodes2)]...)

backend = pnl.DirectBackend()

println("\n== sizes ==")
for (i, bi) in enumerate(bodies)
    println("body $i: ncells=$(bi.ncells)  watertight=$(bi.watertight)  DBC=$(typeof(bi).parameters[4])")
end

npanels = [bi.ncells for bi in bodies]
offsets = cumsum(vcat(0, npanels))
N       = sum(npanels)

# --- Step A: Assemble G ourselves via _G!, no solve! call ---
println("\n== Assembling G directly via _G! ==")
TF = Float64
G  = zeros(TF, N, N)

# Set cp_outer/normals/controlpoints the way BackslashCoupled would
for bi in bodies
    bi.cp_outer = !pnl.has_dirichlet_bc(bi)   # Dirichlet ⇒ false (interior)
    pnl.calc_normals!(bi)
    pnl.calc_controlpoints!(bi)
end

for (sj, src) in enumerate(bodies)
    cols = offsets[sj]+1 : offsets[sj+1]
    for (ti, tgt) in enumerate(bodies)
        rows = offsets[ti]+1 : offsets[ti+1]
        pnl._G!(view(G, rows, cols), tgt, src; update_geometry=false)
    end
end

println("size(G)       = ", size(G))
println("norm(G)       = ", norm(G))
println("norm(G - G')  = ", norm(G - G'))

println("\n== SVD of true G ==")
sv = svdvals(G)
println("svd top 5   : ", sv[1:5])
println("svd bottom 5: ", sv[end-4:end])
println("cond(G)     = ", sv[1] / sv[end])

tol = max(N, N) * eps(sv[1])
rk  = count(>(tol), sv)
println("numerical rank = $rk / $N  (tol=$tol)")

println("\n== Smallest right-singular vectors ==")
F = svd(G)
for k in 0:min(2, N-1)
    v = F.V[:, end-k]
    println("k=$k  sigma=$(F.S[end-k])  max|v|=$(maximum(abs,v))  min|v|=$(minimum(abs,v))  mean(v)=$(sum(v)/length(v))  std(v)=$(sqrt(sum((v .- sum(v)/length(v)).^2)/length(v)))")
end

# --- Step B: Build rhs the way BackslashCoupled does ---
println("\n== Build rhs ==")
# save freestream potential at CPs
phi_ext = zeros(TF, N)
for (bi, body) in enumerate(bodies)
    rng = offsets[bi]+1 : offsets[bi+1]
    phi_ext[rng] .= body.potential
end

for body in bodies
    pnl.set_strengths!(body)   # σ = -V·n, μ = 0
    body.potential .= 0.0
end
pnl.influence!(bodies, bodies, backend;
    scalar_potential=[pnl.has_dirichlet_bc(b) for b in bodies],
    velocity=[!pnl.has_dirichlet_bc(b) for b in bodies])

for (bi, body) in enumerate(bodies)
    rng = offsets[bi]+1 : offsets[bi+1]
    body.potential .+= phi_ext[rng]
end

rhs = zeros(TF, N)
for (bi, body) in enumerate(bodies)
    rng = offsets[bi]+1 : offsets[bi+1]
    pnl.boundary_condition!(body, view(rhs, rng), backend)
end
println("norm(rhs) = ", norm(rhs))

# --- Step C: solve and compare to backslash_meshes residual ---
println("\n== Solve G \\ rhs ==")
x = G \ rhs
println("max|x|         = ", maximum(abs, x))
println("median|x|      = ", sort(abs.(x))[div(N+1, 2)])
println("norm(G*x - rhs)/norm(rhs) = ", norm(G*x - rhs) / norm(rhs))

# --- Step D: matvec ≡ influence! check ---
println("\n== matvec G*x vs influence! ==")
x_test = sin.(LinRange(0, 2π, N))

y_mat = G * x_test

for (bi, body) in enumerate(bodies)
    rng = offsets[bi]+1 : offsets[bi+1]
    body.strength[:, 1] .= 0.0
    body.strength[:, 2] .= x_test[rng]
    body.potential .= 0.0
    body.cp_outer = false
end
pnl.influence!(bodies, bodies, backend;
    scalar_potential=[true for _ in bodies],
    velocity=[false for _ in bodies])
y_phi = vcat([copy(b.potential) for b in bodies]...)

diff = y_mat - y_phi
println("max|G*x|       = ", maximum(abs, y_mat))
println("max|phi|       = ", maximum(abs, y_phi))
println("max|G*x - phi| = ", maximum(abs, diff))
println("max rel diff   = ", maximum(abs.(diff) ./ (abs.(y_mat) .+ 1e-30)))
ix = argmax(abs.(diff))
println("argmax diff at i=$ix:  G*x=$(y_mat[ix])  phi=$(y_phi[ix])  diff=$(diff[ix])")

# --- Step E: Single-panel probes ---
println("\n== Single-panel probes ==")
for k_probe in (1, div(npanels[1], 2), npanels[1]+1)
    e = zeros(N); e[k_probe] = 1.0
    y_mat_e = G * e

    for (bi, body) in enumerate(bodies)
        rng = offsets[bi]+1 : offsets[bi+1]
        body.strength[:, 1] .= 0.0
        body.strength[:, 2] .= e[rng]
        body.potential .= 0.0
        body.cp_outer = false
    end
    pnl.influence!(bodies, bodies, backend;
        scalar_potential=[true for _ in bodies],
        velocity=[false for _ in bodies])
    y_phi_e = vcat([copy(b.potential) for b in bodies]...)

    diff_e = y_mat_e - y_phi_e
    println("probe panel $k_probe:")
    println("  G[$k_probe,$k_probe]   = $(G[k_probe, k_probe])    phi[$k_probe] = $(y_phi_e[k_probe])")
    println("  max|G*e - phi| = $(maximum(abs, diff_e))  argmax at $(argmax(abs.(diff_e)))")
end

println("\ndone.")
