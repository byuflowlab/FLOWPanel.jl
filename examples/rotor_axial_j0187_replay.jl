## Replay the axial comparison and accept only a flat final full revolution.

function setdefaultenv!(key, value)
    haskey(ENV, key) || (ENV[key] = string(value))
end

setdefaultenv!("RUN_NAME", "rotor_axial_j0187_ccblade")
setdefaultenv!("SAVE_PATH", joinpath("data", "rotor_axial_j0187_ccblade"))
setdefaultenv!("OUT_DIR", joinpath("data", "rotor_axial_j0187_ccblade", "spanwise_loading_replay"))
setdefaultenv!("RPM", 5400)
setdefaultenv!("NREVS_AVG", 1)
setdefaultenv!("REQUIRE_FINAL_FULL_REV", true)
setdefaultenv!("CCBLADE_CSV", joinpath("data", "rotor_axial_j0187_ccblade", "rotor_hover_ccblade_sectional_ncrit4_Vc4_J0p1867.csv"))
include(joinpath(@__DIR__, "rotor_hover_spanwise_loading_replay.jl"))
