## Replay the last 10 saved steps of rotor_hover_pressure_comparison with a
## CylindricalFieldProbeMonitor, writing per-entity velocity / velocity-gradient
## fields on a cylindrical grid for ParaView inspection of the blowup.

import FLOWPanel as pnl

ENV["RHPC_SETUP_ONLY"] = "true"
include(joinpath(pnl.examples_path, "rotor_hover_pressure_comparison.jl"))

replay_save_path = joinpath("data", "rotor_hover_pressure_comparison_fieldprobe")
isdir(replay_save_path) || mkpath(replay_save_path)

# Light default grid: ~4k probes
probe_monitor = pnl.CylindricalFieldProbeMonitor(;
    i_frame      = 1,
    axial_axis   = :x,                       # dji9443 spins about global x
    axial_range  = (-0.5R, 4R),              # matches GlobalCylinder maintenance bbox
    radial_range = (0.0, 1.5R),
    azimuth_range = (0.0, 2pi),
    n_axial      = 20,
    n_radial     = 8,
    n_azimuth    = 24,
    backend      = backend_wake,
    save_path    = replay_save_path,
    name         = "fieldprobe",
    verbose      = true,
)

source_path = joinpath("data", "rotor_hover_pressure_comparison")
# Bodies stopped saving at step 79 (blowup); particle field went to 99.
# Replay the last 10 body steps (70..79) so all entities have valid state.
replay_steps = 70:79
println("\nReplaying steps $(first(replay_steps))..$(last(replay_steps)) of $(source_path) with CylindricalFieldProbeMonitor...")

result = pnl.replay(source_path, run_name;
    monitors     = (probe_monitor,),
    Uinf         = Uinf,
    backend      = backend,
    backend_wake = backend_wake,
    steps        = replay_steps,
    verbose      = true,
)

println("\nWrote per-entity probe VTKs to $(replay_save_path)")
