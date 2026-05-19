#=##############################################################################
# DESCRIPTION
    Multi-body unsteady simulation: a capped wing (RigidWakeBody) flying in
    ground effect over a flat ground plane (NonLiftingBody with FlatGroundSolver).
    Verifies no-flow-through on the ground at every time step via a monitor.

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

run_name = "flat_ground"
save_path = joinpath("data", run_name)

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================
AOA             = 5.0                           # (deg) angle of attack
magVinf         = 30.0                          # (m/s) freestream velocity
rho             = 1.225                         # (kg/m^3) air density
Vinf            = magVinf * [cosd(AOA), 0.0, sind(AOA)]

# =============================================================================
# WING GEOMETRY — import capped wing mesh (same as simple_wing_capped.jl)
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
# GROUND GEOMETRY — flat ground plane below the wing
# =============================================================================
ground_z        = -3.0                          # (m) ground plane z-position
ground_center   = [chord/2, 0.0, ground_z]     # center below wing mid-chord
ground_normal   = [0.0, 0.0, 1.0]              # pointing up toward wing
ground_radius   = 20.0                          # (m) radius of ground disc
ground_panel_length = 1.0                       # (m) side length of equilateral triangles

ground = pnl.FlatGround(ground_center, ground_normal, ground_radius;
                         panel_length=ground_panel_length)

println("Ground: $(ground.nnodes) nodes, $(ground.ncells) panels")

# =============================================================================
# WAKE SETUP — PanelParticleWake for wing, nothing for ground
# =============================================================================
dt_val = chord / magVinf / 5                    # timestep (~0.2 chord per step)
das_offset = 0.05                               # Das scale (fraction of unit Vinf direction)

wing.Das[1] .= repeat(Vinf / magVinf * das_offset, 1, size(wing.Das[1], 2))

wake_wing = pnl.PanelParticleWake(wing;
                nwakerows=3,
                max_particles=20000,
                method_trailing=pnl.OverlapPPS(1.3, 2),
                method_unsteady=pnl.OverlapPPS(1.3, 2))

wake_ground = nothing

# =============================================================================
# SIMULATION SETUP
# =============================================================================
Uinf(t) = Vinf

n_steps = 21
t_range = range(0.0, step=dt_val, length=n_steps)

# Solvers
solver_wing   = pnl.Backslash(wing)
solver_ground = pnl.FlatGroundSolver(ground)

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
# FLOW TANGENCY MONITOR — check U·n on ground at every step
# =============================================================================
function tangency_monitor(systems, wakes, i_step)
    ground = systems[2]
    normals = pnl.calc_normals(ground)
    Udotn = sum(ground.velocity .* normals, dims=1)
    rms = sqrt(sum(Udotn .^ 2) / ground.ncells)
    maxerr = maximum(abs.(Udotn))
    println("    Ground tangency step $i_step: RMS(U·n) = $(round(rms; sigdigits=4)),  " *
            "max|U·n| = $(round(maxerr; sigdigits=4)),  " *
            "RMS/Vinf = $(round(rms/magVinf*100; sigdigits=4))%")
end

# =============================================================================
# RUN SIMULATION
# =============================================================================
systems = (wing, ground)
wakes   = (wake_wing, wake_ground)
body_solvers = (solver_wing, solver_ground)

println("\nBegin wing + ground simulation ($(n_steps) steps)...")
@time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
    body_solvers, backend, verbose=true,
    monitors=(tangency_monitor,),
    path=save_path, name=run_name
)

# =============================================================================
# RESULTS AND FINAL FLOW TANGENCY VERIFICATION
# =============================================================================
println("\n=== RESULTS ===")

# Forces
F_wing = sum(wing.F, dims=2)
F_ground = sum(ground.F, dims=2)
println("Wing total force:   ", F_wing)
println("Ground total force: ", F_ground)

# Final flow tangency check for both bodies
println("\n--- Flow Tangency (from last simulation step) ---")

for (i, (body, label)) in enumerate(zip(systems, ("Wing", "Ground")))
    normals = pnl.calc_normals(body)
    Udotn = sum(body.velocity .* normals, dims=1)
    rms = sqrt(sum(Udotn .^ 2) / body.ncells)
    maxr = maximum(abs.(Udotn))
    println("  $label: RMS(U·n) = $(round(rms; sigdigits=4)),  " *
            "max|U·n| = $(round(maxr; sigdigits=4)),  " *
            "RMS/Vinf = $(round(rms/magVinf*100; sigdigits=4))%")
end
