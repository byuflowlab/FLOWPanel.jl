# Example: construct a NonLiftingBody with a single triangular panel
# Vertices: [0,0,0], [1,0,0], [0,1,0]

using FLOWPanel
using FLOWPanel.gt.Meshes
import FLOWPanel.PyPlot as plt
import LinearAlgebra: norm

const pnl = FLOWPanel

function make_single_triangle_body()
    pnl = FLOWPanel

    # Create vertices for a single triangle
    vertices = [
        Point(0.0, 0.0, 0.0),
        Point(0.0, 1.0, 0.0),
        Point(-0.5, 0.5, 0.5),
        Point(-0.5, 0.5, -0.5)
    ]

    # Create a single triangle using Meshes.jl
    triangles = [connect((1,2,3)), connect((2,1,4))]
    mesh = SimpleMesh(vertices, triangles)

    # Wrap it as a GridTriangleSurface
    grid = FLOWPanel.gt.GridTriangleSurface(mesh)

    # Construct NonLiftingBody passing `cells` explicitly to avoid calling gt-based grid2cells
    body = pnl.NonLiftingBody{pnl.ConstantSource}(grid)

    return body
end

function make_single_triangle_lifting_body()
    pnl = FLOWPanel

    # Create vertices for a single triangle
    vertices = [
        Point(0.0, 0.0, 0.0),
        Point(0.0, 1.0, 0.0),
        Point(-0.5, 0.5, 0.5),
        Point(-0.5, 0.5, -0.5)
    ]

    # Create a single triangle using Meshes.jl
    triangles = [connect((1,2,3)), connect((2,1,4))]
    mesh = SimpleMesh(vertices, triangles)

    # Wrap it as a GridTriangleSurface
    grid = FLOWPanel.gt.GridTriangleSurface(mesh)

    # shedding info
    shedding = zeros(Int, 6, 1)
    shedding[1, 1] = 1  # panel 1
    shedding[2, 1] = 1  # vertex 1
    shedding[3, 1] = 2  # vertex 2
    shedding[4, 1] = 2  # panel 2
    shedding[5, 1] = 1  # vertex 1
    shedding[6, 1] = 2  # vertex 2

    # Construct NonLiftingBody passing `cells` explicitly to avoid calling gt-based grid2cells
    body = pnl.RigidWakeBody{Union{pnl.ConstantDoublet, pnl.ConstantSource}, 2, Float64}(grid, shedding)

    return body
end

# generate body
body = make_single_triangle_body()

# set unit strength
body.strength[1,1] = 1.9
body.strength[2,1] = 1.2

# compute velocity at a line of points above the triangle
npoints = 100
dz = range(0.0, stop=10.0, length=npoints)
points = hcat([ [0.33; 0.33; z] for z in dz ]...)  # points above centroid of triangle
# points = hcat([[0.33; 0.33; 0.0] .+ [z; z; z] for z in dz ]...)  # points in y from the centroid

# set "solved" flag
FLOWPanel._solvedflag(body, true)

# direct backend
backend = pnl.DirectBackend()
out_direct = zeros(3, npoints)
pnl.Uind!(body, points, out_direct, backend)

# fast multipole backend
backend_fmm = pnl.FastMultipoleBackend(
                                    expansion_order=20,
                                    multipole_acceptance=1.0,
                                    leaf_size=1,
                                )
out_fmm = zeros(3, npoints)
pnl.Uind!(body, points, out_fmm, backend_fmm)

# add solutions to body
# pnl.addfield!(body, "Uind_direct", collect(eachcol(out_direct)), "cell")
# pnl.addfield!(body, "Uind_fmm", collect(eachcol(out_fmm)), "cell")

# plot results
out_direct_z = out_direct[1, :]
out_fmm_z = out_fmm[1, :]
plt.figure("verify panel")
plt.cla()
plt.plot(dz, abs.(out_direct_z - out_fmm_z), label="error")
plt.yscale("log")
plt.xlabel("z (above triangle)")
plt.ylabel("Uind_z (direct - FMM)")
# plt.plot(dz, out_direct_z, label="direct")
# plt.plot(dz, out_fmm_z, label="FMM", "--")
# plt.legend()
plt.tight_layout()
plt.show()

#--- debug wake ---#

pnl.EXTRA_FARFIELD[1] = false
body = make_single_triangle_lifting_body()

# get freestream
Uinfs = zeros(3, body.ncells)
alpha = 25.0 * pi/180
Uinfs[1, :] .= cos(alpha)
Uinfs[3, :] .= sin(alpha)
body.Das .= repeat(Uinfs[:, 1] ./ norm(Uinfs[:, 1]), 1, body.nsheddings)
body.Dbs .= repeat(Uinfs[:, 1] ./ norm(Uinfs[:, 1]), 1, body.nsheddings)

# get solver
solver = pnl.BackslashDirichlet(body)

# solve body
backend = pnl.FastMultipoleBackend(
                                    expansion_order=7,
                                    multipole_acceptance=0.4,
                                    leaf_size=200000
                                )
pnl.solve2!(body, Uinfs, solver; backend)

# induced velocity at control points
Us = pnl.calcfield_U(body, body; backend)

# visualize results
save_path="two_triangle_body"
name="two_triangle_body"
# normals = pnl.calc_normals(body)
# add_field(body, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
# add_field(body, "Da", "vector", collect(eachcol(body.Das)), "system")
# add_field(body, "Db", "vector", collect(eachcol(body.Dbs)), "system")
# add_field(body, "sigma", "scalar", view(body.strength, :, 1), "cell")
# add_field(body, "mu", "scalar", view(body.strength, :, 2), "cell")
# add_field(body, "normals", "vector", collect(eachcol(normals)), "cell")
pnl.save(body, name; path=save_path, wake_panel=true, debug=true)
