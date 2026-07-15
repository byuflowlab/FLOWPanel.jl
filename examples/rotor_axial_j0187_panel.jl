## DJI 9443 axial-flow panel case for the CCBlade comparison (J = 0.1867).
##
## This is deliberately a small entry point over rotor_hover_pressure_comparison:
## it preserves the established 80_81 geometry, wake, monitors, and VTK/replay
## metadata while replacing the startup/withdrawal pulse by V∞ = 4 m/s.

function setdefaultenv!(key, value)
    haskey(ENV, key) || (ENV[key] = string(value))
end

setdefaultenv!("RUN_NAME", "rotor_axial_j0187_ccblade")
setdefaultenv!("RPM", 5400)
setdefaultenv!("RHPC_MESH", "80_81")
setdefaultenv!("NT", 36)
setdefaultenv!("NREVS", 12)
setdefaultenv!("SAVE_VTK", true)
setdefaultenv!("MAGVINF_START", 4.0)
setdefaultenv!("MAGVINF_PEAK", 4.0)
setdefaultenv!("MAGVINF_END", 4.0)
setdefaultenv!("FREESTREAM_RAMP_REVS", 0.0)
setdefaultenv!("FREESTREAM_HOLD_REVS", 0.0)
setdefaultenv!("FREESTREAM_WITHDRAW_REVS", 0.0)
setdefaultenv!("SETTLE_REVS", 0.0)

println("Axial CCBlade comparison: V∞=4 m/s, RPM=5400, J=$(4 / ((5400 / 60) * (2 * 0.119)))")
include(joinpath(@__DIR__, "rotor_hover_pressure_comparison.jl"))
