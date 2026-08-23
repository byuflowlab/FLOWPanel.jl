## Rotor hover: compare CT from PressureLaplace vs. KuttaJoukowskiForce
# Derived from rotor_hover.jl. The pressure-monitor swap replaces PressureBernoulli
# with PressureLaplace (surface pressure Poisson solve). KuttaJoukowskiForce is
# enabled as an independent cross-check that bypasses the pressure field entirely.
# Both forces are normalized as rotor CT via RotorNormalization(rho, 2R, 1).

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
using FLOWPanel.FastMultipole.StaticArrays
using VSPGeom

run_name = "rotor_hover_pressurelaplace"
save_path = joinpath("data", run_name)

## =========================================================
# SIMULATION PARAMETERS
# ==========================================================
magVinf = 0.0001    # Freestream velocity magnitude (m/s)
AOA     = 0.0       # Angle of attack (degrees)
rho     = 1.225     # Air density (kg/m^3)
RPM     = 6000      # Rotation speed (rpm)
Vinf    = magVinf * [0.0, -cosd(AOA), sind(AOA)]
R       = 0.119      # Rotor radius (m)
nrevs   = 1        # Number of revolutions
nt      = 36        # Number of time steps per revolution
dt      = 60 / RPM / nt
n_steps = nt * nrevs
t_range = range(0.0, step=dt, length=n_steps)[1:3]

# ==========================================================
# Sensitivity parameters
# ==========================================================
cp_outer=true
core_size = R * 1e-3
kernelcutoff = R * 1e-13
p_per_step   = 2
overlap      = 2.0
merge_r_factor = 0.02
merge_r_hash_factor = 0.04
init_Das_eta_kinematic = 0.2
p_correct_kuttacondition_flag = false
wake_core_size = parse(Float64, get(ENV, "WAKE_CORE_SIZE", "1e-3"))

## =========================================================
# ROTOR GEOMETRY
# ==========================================================
read_path   = joinpath(pnl.examples_path, "data")
msh_file  = joinpath(read_path, "phantom_3_rebuild_r2.msh")
te_indices_1 = [9, 175, 127]
te_indices_2 = [13, 286, 238]

msh = pnl.read_gmsh(msh_file)
nodes, cells = pnl.meshes2nodes_cells(msh)
nodes .*= R / maximum(nodes[1, :])

shedding = pnl.noshedding

kernel = Union{pnl.ConstantSource, pnl.VortexRing}
DBC = kernel == pnl.VortexRing ? false : true
rotor = pnl.RigidWakeBody{kernel}(nodes, cells, shedding;
            core_size, kernelcutoff,
            semiinfinite_wake=false, watertight=true, DBC)

shedding1 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_1[1], te_indices_1[2];
                bbox=nothing, end_node=te_indices_1[3], normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
shedding2 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_2[1], te_indices_2[2];
                bbox=nothing, end_node=te_indices_2[3], normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

rotor = pnl.RigidWakeBody{kernel}(rotor.nodes, rotor.cells, [shedding1, shedding2];
                        core_size, kernelcutoff,
                        semiinfinite_wake=false, watertight=true,
                        ensure_winding=true, DBC)

pnl.write_vtk(joinpath(save_path, run_name), rotor)
println("Rotor: $(rotor.nnodes) nodes, $(rotor.ncells) panels, $(rotor.nsheddings) shedding edges")

## =========================================================
# WAKE SETUP
# ==========================================================
wake_rotor = pnl.PanelParticleWake(rotor;
                nwakerows=1,
                core_size=wake_core_size,
                max_particles=100000,
                method_trailing=pnl.OverlapPPS(overlap, p_per_step),
                method_unsteady=pnl.OverlapPPS(overlap, p_per_step),
                particle_maintenance=pnl.ParticleMaintenance((
                    pnl.MergeParticles(every=1,
                        r=R*merge_r_factor,
                        r_hash=R*merge_r_hash_factor,
                        sigma_relative=false,
                        max_sigma_ratio=2.0,
                        skip_static=true),
                )))

## =========================================================
# SIMULATION SETUP
# ==========================================================
Uinf(t) = Vinf

frames = pnl.ReferenceFrame(rotor;
    origin = SVector{3}(0.0, 0.0, 0.0),
    v = SVector{3}(0.0, 0.0, 0.0),
    ω_axis = SVector{3}(0.0, 1.0, 0.0),
    ω = 2*pi * RPM/60,
    R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    name = "vehicle",
    child_index = Int[],
    dependent_index = [1]
)

pnl.initialize_Das!((rotor,), frames, Uinf, t_range[1], t_range[2] - t_range[1];
    set_Das_eta_kinematic=init_Das_eta_kinematic)

solver_rotor = pnl.Backslash(rotor)
backend = pnl.FastMultipoleBackend()
maneuver!(frames, systems, wakes, t) = nothing

## =========================================================
# MONITORS — PressureLaplace vs. Kutta–Joukowski cross-check
# ==========================================================
systems      = (rotor,)
wakes        = (wake_rotor,)
body_solvers = (solver_rotor,)

# PressureLaplace populates body.P each step; ForceMonitor then integrates body.F
# from body.P. KuttaJoukowskiForce is independent of body.P entirely — it sums
# ρ Σ γ (Δs × V) at panel edges using a FastMultipole.ProbeSystem.
pressure_laplace = pnl.PressureLaplace(rotor, rho; verbose=false)
force_monitor    = pnl.ForceMonitor(length(t_range), 1;
                        i_frame=-1,
                        normalization=pnl.RotorNormalization(rho, 2*R, 1),
                        correct_kuttacondition=p_correct_kuttacondition_flag,
                        verbose=false)
kj_monitor       = pnl.KuttaJoukowskiForce(rotor, length(t_range), 1;
                        rho, backend,
                        normalization=pnl.RotorNormalization(rho, 2*R, 1),
                        verbose=false)
monitors = (pressure_laplace, force_monitor, kj_monitor)

## =========================================================
# RUN SIMULATION
# ==========================================================
println("\nBegin rotor hover simulation ($(length(t_range)) steps)...")
name = "rotor_hover_pressurelaplace"
@time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
    set_Das_eta_kinematic=NaN,
    monitors,
    body_solvers, backend, verbose=true,
    path=save_path, name,
)

## =========================================================
# COMPARE CT
# ==========================================================
# Rotor axis is +y (ω_axis = [0,1,0]); thrust is the y-component of force.
CT_pressure = force_monitor.force[2, :]
CT_kj       = kj_monitor.force[2, :]

println("\n========== CT comparison ==========")
println("step | CT (PressureLaplace+ForceMonitor) | CT (Kutta–Joukowski) | abs diff | rel diff")
for k in 1:length(t_range)
    cp = CT_pressure[k]
    ck = CT_kj[k]
    ad = abs(cp - ck)
    rd = abs(ck) > 0 ? ad / abs(ck) : NaN
    println("  $k  |  $(round(cp, sigdigits=6))  |  $(round(ck, sigdigits=6))  |  $(round(ad, sigdigits=4))  |  $(round(rd, sigdigits=4))")
end

println("\nFull force arrays (3 × nt) printed for inspection:")
println("  PressureLaplace+ForceMonitor force:")
display(force_monitor.force)
println("  Kutta–Joukowski force:")
display(kj_monitor.force)
