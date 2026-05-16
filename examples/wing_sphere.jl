#=##############################################################################
# DESCRIPTION
    Multi-body unsteady simulation: a low-AR capped wing (RigidWakeBody) with
    PanelParticleWake interacting with a downstream sphere (NonLiftingBody).
    Verifies flow tangency for both bodies after simulation.

# AUTHORSHIP
  * Author    : Ryan Anderson and Claude (AI assistant)
  * Created   : Apr 2026
  * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
import GeometricTools as gt
using FLOWPanel.FastMultipole.StaticArrays
import LinearAlgebra: norm, I
import Meshes
import GeoIO

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================
AOA             = 5.0                           # (deg) angle of attack
magVinf         = 30.0                          # (m/s) freestream velocity
rho             = 1.225                         # (kg/m^3) air density
Vinf            = magVinf * [cosd(AOA), 0.0, sind(AOA)]

# =============================================================================
# WING GEOMETRY — import the capped wing mesh used by simple_wing_capped.jl
# =============================================================================

chord           = 2.0                           # (m) root chord length
AR              = 4.0                           # aspect ratio b/c
b               = AR * chord                    # (m) span
read_path       = joinpath(pnl.examples_path, "data")
meshfile        = joinpath(read_path, "wing_ar4_naca0016_5.msh")

msh = GeoIO.load(meshfile).geometry
msh = msh |> Meshes.Scale(1.0)
grid = gt.GridTriangleSurface(msh)

nte = 10000
trailingedge = zeros(3, nte)
trailingedge[1, :] .= chord
trailingedge[2, :] .= range(-b/2, stop=b/2, length=nte)
trailingedge[3, :] .= 0.0

# --- Construct RigidWakeBody ---
kernel = Union{pnl.ConstantSource, pnl.VortexRing}
wing = pnl.RigidWakeBody{kernel}(grid;
            CPoffset=1e-14,
            kerneloffset=1e-2,
            kernelcutoff=1e-14,
            semiinfinite_wake=false,
            watertight=true)
shedding = pnl.calc_shedding(wing.nodes, wing.cells, trailingedge; tolerance=0.001 * b)
wing = pnl.RigidWakeBody{kernel}(wing.nodes, wing.cells, shedding;
            CPoffset=1e-14,
            kerneloffset=1e-2,
            kernelcutoff=1e-14,
            semiinfinite_wake=false,
            watertight=true,
            ensure_winding=false)

println("Wing: $(wing.nnodes) nodes, $(wing.ncells) panels, $(wing.nsheddings) shedding edges")

# =============================================================================
# SPHERE GEOMETRY — NonLiftingBody downstream of the wing
# =============================================================================
R_sphere        = 0.5                           # (m) sphere radius
x_sphere        = chord + 3 * R_sphere          # (m) sphere center x-position

# Parametric grid in (theta, phi)
theta_pad = 0.15                                # avoid polar singularity
P_min = [theta_pad, 0.0, 0.0]
P_max = [π - theta_pad, 2π, 0.0]
NDIVS_sphere = [20, 40, 0]

sphere_grid = gt.Grid(P_min, P_max, NDIVS_sphere; loop_dim=2)

# Transform to Cartesian spherical coordinates
gt.transform!(sphere_grid, X -> gt.spherical3D(vcat(R_sphere, X[1:2])))

# Translate downstream
Oaxis_sphere = Matrix{Float64}(I, 3, 3)
O_sphere = [x_sphere, 0.0, 0.0]
gt.lintransform!(sphere_grid, Oaxis_sphere, O_sphere)

# Triangulate and create body
triang_sphere = gt.GridTriangleSurface(sphere_grid, 1)
sphere = pnl.NonLiftingBody{pnl.ConstantSource}(triang_sphere)

println("Sphere: $(sphere.nnodes) nodes, $(sphere.ncells) panels, center at $O_sphere")

# =============================================================================
# WAKE SETUP — PanelParticleWake for wing, nothing for sphere
# =============================================================================
dt_val = chord / magVinf / 5                    # timestep (~0.2 chord per step)
das_offset = 0.05                               # Das scale (fraction of unit Vinf direction)

wing.Das[1] .= repeat(Vinf / magVinf * das_offset, 1, size(wing.Das[1], 2))

wake_wing = pnl.PanelParticleWake(wing;
                nwakerows=3,
                max_particles=20000,
                method_trailing=pnl.OverlapPPS(1.3, 2),
                method_unsteady=pnl.OverlapPPS(1.3, 2))

wake_sphere = nothing

# =============================================================================
# SIMULATION SETUP
# =============================================================================
Uinf(t) = Vinf

n_steps = 21
t_range = range(0.0, step=dt_val, length=n_steps)

# Solvers
solver_wing   = pnl.Backslash(wing)
solver_sphere = pnl.Backslash(sphere)

# Backend
backend = pnl.FastMultipoleBackend()

# Reference frame (static — no rotation or translation)
frames = pnl.ReferenceFrame(wing;
    origin = SVector{3}(0.0, 0.0, 0.0),
    v = SVector{3}(0.0, 0.0, 0.0),
    ω_axis = SVector{3}(0.0, 1.0, 0.0),
    ω = 0.0,
    R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    name = "vehicle",
    child_index = Int[],
    dependent_index = [1, 2]
)

# Maneuver (no-op for static case)
maneuver!(frames, systems, wakes, t) = nothing

# =============================================================================
# RUN SIMULATION
# =============================================================================
systems = (wing, sphere)
wakes   = (wake_wing, wake_sphere)
body_solvers = (solver_wing, solver_sphere)

println("\nBegin wing+sphere simulation ($(n_steps) steps)...")
@time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
    body_solvers, backend, verbose=true,
    path="wing_sphere", name="wing_sphere"
)

# =============================================================================
# RESULTS AND FLOW TANGENCY VERIFICATION
# =============================================================================
println("\n=== RESULTS ===")

# Forces
F_wing = sum(wing.F, dims=2)
F_sphere = sum(sphere.F, dims=2)
println("Wing total force:   ", F_wing)
println("Sphere total force: ", F_sphere)

# Flow tangency — check velocity field left by the last simulation step.
# simulate! leaves body.velocity = freestream + kinematic + wake + body influence,
# which should satisfy U·n ≈ 0 if the solve converged.
println("\n--- Flow Tangency (from last simulation step) ---")

for (i, (body, label)) in enumerate(zip(systems, ("Wing", "Sphere")))
    Udotn = sum(body.velocity .* body.normals, dims=1)
    rms = sqrt(sum(Udotn .^ 2) / body.ncells)
    maxr = maximum(abs.(Udotn))
    println("  $label: RMS(U·n) = $(round(rms; digits=6)),  " *
            "max|U·n| = $(round(maxr; digits=6)),  " *
            "RMS/Vinf = $(round(rms/magVinf*100; digits=4))%")
end
