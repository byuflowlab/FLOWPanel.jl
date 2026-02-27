# Example: construct a NonLiftingBody with a single triangular panel
# Vertices: [0,0,0], [1,0,0], [0,1,0]

using FLOWPanel
const pnl = FLOWPanel
using FLOWPanel.gt.Meshes
import FLOWPanel.PyPlot as plt

function make_single_triangle_body()
    pnl = FLOWPanel

    # Create vertices for a single triangle
    vertices = [
        Point(0.0, 0.0, 0.0),
        Point(1.0, 0.0, 0.0),
        Point(0.0, 1.0, 0.0)
        # Point(0.4972949032316967, -8.097796542024505e-17, -0.9098827194222483),
        # Point(0.36953825743300084, 0.18883245189104392, -0.8883868388132481),
        # Point(0.36953825743300084, -8.859888522845575e-17, -0.9082339292625988),
    ]

    # Create a single triangle using Meshes.jl
    triangle = [connect((1,2,3))]
    mesh = SimpleMesh(vertices, triangle)

    # Wrap it as a GridTriangleSurface
    grid = FLOWPanel.gt.GridTriangleSurface(mesh)

    # Construct NonLiftingBody passing `cells` explicitly to avoid calling gt-based grid2cells
    # body = pnl.NonLiftingBody{pnl.ConstantSource}(grid)
    body = pnl.NonLiftingBody{pnl.ConstantDoublet}(grid)
    # body = pnl.NonLiftingBody{Union{pnl.ConstantSource, pnl.ConstantDoublet}}(grid)

    return body
end

# generate body
body = make_single_triangle_body()

# set unit strength
str = 1.0
body.strength[1,1] = str
# body.strength[1,2] = str

normals = pnl._calc_normals(body)
@show normals # outward facing
cps = pnl._calc_controlpoints(body, normals; off=1e-14)

# compute velocity at a line of points above the triangle
npoints = 1000
dz = range(-1.0, stop=1.0, length=npoints)
points = hcat([ cps[:,1] .+ [0.0,0.0,z] for z in dz ]...)  # points above centroid of triangle
# points = hcat([ cps[:,1] .+ z .* normals[:,1] for z in dz ]...)  # points above centroid of triangle
# points = hcat([[0.33; 0.33; 0.0] .+ [z; z; z] for z in dz ]...)  # points in y from the centroid

# set "solved" flag
FLOWPanel._solvedflag(body, true)

# direct backend
backend = pnl.FastMultipoleBackend(
                                    expansion_order=20,
                                    multipole_acceptance=1.0,
                                    leaf_size=10000,
                                )
out_direct = zeros(3, npoints)
pnl.Uind!(body, points, out_direct, backend)
# out_direct .*= normals[:,1]
phi_direct = zeros(npoints)
pnl.phi!(body, points, phi_direct, backend)

# get semi-infinite panel potential and velocities
phi_semiinfinite = zeros(npoints)
U_semiinfinite = zeros(3, npoints)
v1x, v1y, v1z = 1.0, 0.0, 0.0
v2x, v2y, v2z = 0.0, 1.0, 0.0
d1, d2, d3 = 0.0, -1.0, 0.0
for i in 1:npoints
    phi_semiinfinite[i] = pnl._phi_semiinfinite(view(points,:,i), pnl.ConstantDoublet, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, str; kerneloffset=1e-8)
    U_semiinfinite[:,i] .= pnl._U_semiinfinite(view(points,:,i), pnl.ConstantDoublet, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, str; kerneloffset=1e-8)
end

# fast multipole backend
backend_fmm = pnl.FastMultipoleBackend(
                                    expansion_order=20,
                                    multipole_acceptance=1.0,
                                    leaf_size=1,
                                )
out_fmm = zeros(3, npoints)
pnl.Uind!(body, points, out_fmm, backend_fmm)
# out_fmm .*= normals[:,1]
phi_fmm = zeros(npoints)
pnl.phi!(body, points, phi_fmm, backend_fmm)

# plot results
out_direct_z = out_direct[3, :]
out_fmm_z = out_fmm[3, :]
fig = plt.figure("verify panel")
fig.clf()
fig.add_subplot(121, xlabel="z", ylabel="phi")
fig.add_subplot(122, xlabel="z", ylabel="Uind_z")
axs = fig.get_axes()
axs[2].plot(dz, out_direct_z, label="direct")
axs[2].plot(dz, out_fmm_z, "--", label="FMM")
axs[2].plot(dz, U_semiinfinite[3,:], ":", label="semi-infinite panel")
axs[2].legend()
axs[1].plot(dz, phi_direct, label="direct")
axs[1].plot(dz, phi_fmm, "--", label="FMM")
axs[1].plot(dz, phi_semiinfinite, ":", label="semi-infinite panel")
axs[1].legend()
plt.tight_layout()
plt.show()
# plt.plot(dz, abs.(out_direct_z - out_fmm_z), label="error")
# plt.yscale("log")
# plt.xlabel("z (above triangle)")
# plt.ylabel("Uind_z (direct - FMM)")
# # plt.plot(dz, out_direct_z, label="direct")
# # plt.plot(dz, out_fmm_z, label="FMM", "--")
# # plt.legend()
# plt.tight_layout()
# plt.show()
