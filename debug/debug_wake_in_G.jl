using Test
import FLOWPanel as pnl
using LinearAlgebra
import Meshes
using StaticArrays: SVector
import GeoIO

include(joinpath(@__DIR__, "..", "test", "test_helpers.jl"))

run_names = ["nasa_wing.msh", "nasa_surface_spaced_repaired.msh"]
files = [joinpath(pnl.examples_path, "wing_aileron", name) for name in run_names]
nodes1 = [42, 19]
nodes2 = [34, 3]
m = 0.0254
magVinf = 117.3 * m * 12
c_body1 = 10 * m
b = 60 * m
c_body2 = 2.0
chords = [c_body1, c_body2]
kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}
bodytype = pnl.RigidWakeBody{kernel}
Vinf = magVinf * [1.0, 0.0, 0.0]

bodies = tuple([generate_body(file, chord, b, bodytype, scaling, 1, Vinf, firstnode, secondnode)
                for (file, chord, scaling, firstnode, secondnode) in zip(files, chords, [1.0, 1.0], nodes1, nodes2)]...)

backend = pnl.DirectBackend()
for bi in bodies
    bi.cp_outer = false
    pnl.calc_normals!(bi)
    pnl.calc_controlpoints!(bi)
end

# count TE panels per body
for (i, bi) in enumerate(bodies)
    te = count(>(0), bi.shedding_full[1, :])
    println("body $i: ncells=$(bi.ncells), TE panels (shedding_full[1,j]>0): $te")
end

# pick one TE panel from body 1, one non-TE from body 1
b1 = bodies[1]
te_panels = findall(>(0), b1.shedding_full[1, :])
non_te    = setdiff(1:b1.ncells, te_panels)
println("\nbody1 TE panels (first 8): ", te_panels[1:min(8,end)])
println("body1 non-TE example: ", non_te[1])

# --- compare column of G for TE vs non-TE built via _G! (which calls induced→_induced_wake) ---
npanels = [bi.ncells for bi in bodies]
offsets = cumsum(vcat(0, npanels))
N = sum(npanels)
G = zeros(N, N)
for (sj, src) in enumerate(bodies)
    cols = offsets[sj]+1 : offsets[sj+1]
    for (ti, tgt) in enumerate(bodies)
        rows = offsets[ti]+1 : offsets[ti+1]
        pnl._G!(view(G, rows, cols), tgt, src; update_geometry=false)
    end
end

# Now build G_nowake: same thing but with wake suppressed (set shedding_full[1, .] = 0 temporarily)
sf1 = copy(bodies[1].shedding_full)
sf2 = copy(bodies[2].shedding_full)
bodies[1].shedding_full[1, :] .= 0
bodies[2].shedding_full[1, :] .= 0
G_nowake = zeros(N, N)
for (sj, src) in enumerate(bodies)
    cols = offsets[sj]+1 : offsets[sj+1]
    for (ti, tgt) in enumerate(bodies)
        rows = offsets[ti]+1 : offsets[ti+1]
        pnl._G!(view(G_nowake, rows, cols), tgt, src; update_geometry=false)
    end
end
bodies[1].shedding_full .= sf1
bodies[2].shedding_full .= sf2

D = G - G_nowake   # pure wake contribution per source column
println("\nnorm(G - G_nowake) = ", norm(D))
println("max|D|             = ", maximum(abs, D))

# For each TE panel, check that its column in D is nonzero
for j in te_panels[1:min(4,end)]
    println("D[:, $j]:  norm=$(norm(D[:,j]))  max|.|=$(maximum(abs, D[:,j]))")
end
# For non-TE panel, column of D should be zero
println("D[:, non-TE $(non_te[1])]:  norm=$(norm(D[:,non_te[1]]))")

# Are TE-panel wake influences in body2 columns affecting body1 rows (and vice versa)?
te2 = findall(>(0), bodies[2].shedding_full[1, :])
j_te2 = offsets[2] + te2[1]
println("D[body1 rows, body2 TE col $(j_te2)] norm = ", norm(D[1:offsets[2], j_te2]))

# --- Now test if G with wake included is full rank ---
println("\n== SVD comparison ==")
sv     = svdvals(G)
sv_nw  = svdvals(G_nowake)
println("with wake    : smallest 4 = ", sv[end-3:end])
println("without wake : smallest 4 = ", sv_nw[end-3:end])
println("with wake    : cond = ", sv[1]/sv[end])
println("without wake : cond = ", sv_nw[1]/sv_nw[end])

# Probe: how does G act on the uniform-doublet-on-body-1 vector?
v1 = zeros(N)
v1[1:offsets[2]] .= 1.0
println("\n== G * (uniform μ on body 1) ==")
y_wake   = G * v1
y_nowake = G_nowake * v1
y_wake_only = D * v1
println("max|G * v1|         = ", maximum(abs, y_wake))
println("max|G_nowake * v1|  = ", maximum(abs, y_nowake))
println("max|D * v1|         = ", maximum(abs, y_wake_only))
println("(D*v1) at body1 rows: max|.|=$(maximum(abs, y_wake_only[1:offsets[2]])), at body2 rows: max|.|=$(maximum(abs, y_wake_only[offsets[2]+1:end]))")

# Same for body 2 uniform
v2 = zeros(N); v2[offsets[2]+1:end] .= 1.0
println("max|G  * v2|        = ", maximum(abs, G * v2))
println("max|D  * v2|        = ", maximum(abs, D * v2))
