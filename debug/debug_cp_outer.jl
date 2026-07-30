using Test
import FLOWPanel as pnl
using LinearAlgebra
import Meshes
using StaticArrays: SVector
import GeoIO

include(joinpath(@__DIR__, "..", "test", "test_helpers.jl"))

run_names = ["nasa_wing.msh"]
files = [joinpath(pnl.examples_path, "wing_aileron", name) for name in run_names]
m = 0.0254
magVinf = 117.3 * m * 12
c_body1 = 10 * m
b = 60 * m
kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}
bodytype = pnl.RigidWakeBody{kernel}
Vinf = magVinf * [1.0, 0.0, 0.0]

body = generate_body(files[1], c_body1, b, bodytype, 1.0, 1, Vinf, 42, 19)

println("After construction: body.cp_outer = ", body.cp_outer)
println("has_dirichlet_bc(body) = ", pnl.has_dirichlet_bc(body))

# Probe: what does induced return for a self pair, at cp_outer=false vs true?
pnl.calc_normals!(body)

# pick a panel, place a unit doublet on just that panel
i = 56
body.strength .= 0.0
body.strength[i, 2] = 1.0

# disable wake for this isolated test
sf_save = copy(body.shedding_full)
body.shedding_full[1, :] .= 0

# evaluate induced at the centroid of panel i with cp_outer=false
v1 = body.nodes[:, body.cells[1, i]]
v2 = body.nodes[:, body.cells[2, i]]
v3 = body.nodes[:, body.cells[3, i]]
centroid = (v1 + v2 + v3) / 3
target = SVector(centroid[1], centroid[2], centroid[3])

body.cp_outer = false
ds_phi = FastMultipole.DerivativesSwitch(true, false, false)
phi_int, _, _ = pnl.induced(target, body, i, ds_phi; kerneloffset=body.kerneloffset)
println("\nself-pair, cp_outer=false → phi = ", phi_int, "   (expected interior limit +μ/2 = +0.5)")

body.cp_outer = true
phi_ext, _, _ = pnl.induced(target, body, i, ds_phi; kerneloffset=body.kerneloffset)
println("self-pair, cp_outer=true  → phi = ", phi_ext, "   (expected exterior limit -μ/2 = -0.5)")

# Also check what the raw _induced returns (bypassing _self_limit) at the centroid
# by setting cp_outer = false AND moving the target slightly perpendicular
body.cp_outer = false
n = body.normals[:, i]
δ = 1e-6 * sqrt(0.5)  # small offset
println("\nslightly OUTSIDE (centroid + δ*n): cp_outer=false")
t_out = SVector(centroid[1]+δ*n[1], centroid[2]+δ*n[2], centroid[3]+δ*n[3])
phi_o, _, _ = pnl.induced(t_out, body, i, ds_phi; kerneloffset=body.kerneloffset)
println("  phi = $phi_o   (true exterior, should be -0.5)")

println("slightly INSIDE  (centroid - δ*n): cp_outer=false")
t_in = SVector(centroid[1]-δ*n[1], centroid[2]-δ*n[2], centroid[3]-δ*n[3])
phi_i, _, _ = pnl.induced(t_in, body, i, ds_phi; kerneloffset=body.kerneloffset)
println("  phi = $phi_i   (true interior, should be +0.5)")

# restore
body.shedding_full .= sf_save

# Also: what is the GeometricTools normal direction vs outward?
println("\nbody.normals[:, $i] = ", body.normals[:, i])
println("centroid = ", centroid)
