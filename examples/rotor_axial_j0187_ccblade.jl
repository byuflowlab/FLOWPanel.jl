## Explicit DJI 9443 CCBlade operating point used by the axial panel comparison.

function setdefaultenv!(key, value)
    haskey(ENV, key) || (ENV[key] = string(value))
end

setdefaultenv!("SAVE_PATH", joinpath(@__DIR__, "..", "data", "rotor_axial_j0187_ccblade"))
setdefaultenv!("RPM", 5400)
setdefaultenv!("VC_LIST", "4.0")
setdefaultenv!("TARGET_VC", 4.0)
setdefaultenv!("OP_TAG", "Vc4_J0p1867")
setdefaultenv!("NCRIT_LIST", "4,9")
setdefaultenv!("STRICT_ALPHA_GATE", true)
include(joinpath(@__DIR__, "rotor_hover_ccblade.jl"))
