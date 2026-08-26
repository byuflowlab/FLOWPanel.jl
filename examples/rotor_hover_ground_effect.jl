## Rotor hover IN GROUND EFFECT (BRAINSTORM/022).
##
## Fork of examples/rotor_hover_pressure_comparison.jl at commit b251071
## (2026-08-18); RHPC itself is frozen for the live 018 campaign and must not
## be modified. This driver adds a flat ground plane (existing machinery:
## FlatGround + FlatGroundSolver + block Gauss-Seidel multi-body solve!) at
## GROUND_H_R rotor radii downstream of the rotor plane, plus ground-specific
## diagnostics. GROUND_ENABLE=false collapses to the single-body RHPC path so
## the out-of-ground-effect (OGE) baseline comes from this identical file.
##
## Operating point defaults (022, overriding 018): RPM=6000, RHO=1.16,
## ROTOR_R=0.1195. Downstream (wake) direction is +axial (+x for dji9443
## meshes): the truncation GlobalCylinder extrudes +x and the ground sits at
## x = x_rotor + GROUND_H_R*R with its normal pointing back at the rotor.
##
## KuttaJoukowskiForce is available as an opt-in diagnostic with RUN_KJ=true.

import FLOWPanel as pnl
# include(joinpath(pnl.examples_path, "helper_functions.jl"))
using FLOWPanel.FastMultipole.StaticArrays
using VSPGeom
using LinearAlgebra: norm, dot, cross
import LinearAlgebra

# BLAS threads: honor BLAS_NUM_THREADS, else OMP_NUM_THREADS, else Julia's
# thread count. Env vars alone are not reliable across BLAS vendors/load
# order, so set it explicitly here.
let n = tryparse(Int, get(ENV, "BLAS_NUM_THREADS", get(ENV, "OMP_NUM_THREADS", string(Threads.nthreads()))))
    if n !== nothing && n > 0
        LinearAlgebra.BLAS.set_num_threads(n)
    end
    println("BLAS threads: $(LinearAlgebra.BLAS.get_num_threads()) (Julia threads: $(Threads.nthreads()))")
end

run_name = get(ENV, "RUN_NAME", "rotor_hover_ground_effect")
save_path = parse(Bool, get(ENV, "SAVE_VTK", "true")) ? joinpath("data", run_name) : nothing

AOA = 0.0
# 022 operating point (explicit override of 018's 1.179 / 0.119):
rho = parse(Float64, get(ENV, "RHO", "1.16"))
RPM = parse(Float64, get(ENV, "RPM", "6000"))
R = parse(Float64, get(ENV, "ROTOR_R", "0.1195"))
# Modeling root clip: the trailing edge inboard of this radius does not shed.
# Kept at 0.1 for comparability across the whole 018 campaign; inert on hub
# meshes, whose blade root is outboard of it. NOT cap protection -- see the
# end_node anchoring below.
shedding_r_over_R = parse(Float64, get(ENV, "SHEDDING_R_OVER_R", "0.1"))
nrevs = parse(Float64, get(ENV, "NREVS", "10"))
nt = parse(Int, get(ENV, "NT", "36"))
dt = 60 / RPM / nt
spinup_revs = max(0.0, parse(Float64, get(ENV, "SPINUP_REVS", "0.0")))
spinup_start_fraction = clamp(parse(Float64, get(ENV, "SPINUP_START_FRACTION", "0.05")), eps(), 1.0)
spinup_duration = spinup_revs * 60 / RPM
spinup_steps = ceil(Int, spinup_duration / dt)
n_steps = round(Int, nt * nrevs) + spinup_steps
t_range = range(0.0, step=dt, length=n_steps)

# --- Item 005: staged-startup convecting-freestream pulse ---------------------
# A four-phase smoothstep profile for magVinf(t), measured in rotor revolutions:
#   ramp-up   : magVinf_start -> magVinf_peak  over FREESTREAM_RAMP_REVS
#   hold      : magVinf_peak                    for FREESTREAM_HOLD_REVS
#   withdraw  : magVinf_peak -> magVinf_end     over FREESTREAM_WITHDRAW_REVS
#   hover     : magVinf_end                     for the remainder
# magVinf_peak defaults to ~2x the hover induced velocity v_h = sqrt(T/(2 rho A))
# (with CT~0.06, R=0.119, RPM=6000: v_h ~ 4.6 m/s) so the freestream initially
# dominates self-induction and sweeps the shed wake clear, then is slowly
# withdrawn so wake self-induction sustains the column into hover.
sec_per_rev = 60 / RPM
magVinf_start  = parse(Float64, get(ENV, "MAGVINF_START",  "0.0"))
magVinf_peak   = parse(Float64, get(ENV, "MAGVINF_PEAK",   "10.0"))
magVinf_end    = parse(Float64, get(ENV, "MAGVINF_END",    "0.0"))
freestream_ramp_revs     = parse(Float64, get(ENV, "FREESTREAM_RAMP_REVS",     "2.0"))
freestream_hold_revs     = parse(Float64, get(ENV, "FREESTREAM_HOLD_REVS",     "3.0"))
freestream_withdraw_revs = parse(Float64, get(ENV, "FREESTREAM_WITHDRAW_REVS", "4.0"))
settle_revs              = parse(Float64, get(ENV, "SETTLE_REVS",              "4.0"))
t_ramp_up   = freestream_ramp_revs     * sec_per_rev
t_hold      = freestream_hold_revs     * sec_per_rev
t_withdraw  = freestream_withdraw_revs * sec_per_rev
cylinder_depth = parse(Float64, get(ENV, "TRUNCATION_DEPTH_R", "4")) * R

# --- BRAINSTORM/022: ground-plane knobs ---------------------------------------
# GROUND_ENABLE=false gives the OGE baseline from this same file (all other
# knobs identical), which is the matched-settings comparison contract.
ground_enable         = parse(Bool,    get(ENV, "GROUND_ENABLE",         "true"))
ground_h_r            = parse(Float64, get(ENV, "GROUND_H_R",            "1.0"))  # rotor plane -> ground, in R
ground_radius_r       = parse(Float64, get(ENV, "GROUND_RADIUS_R",       "4.0"))  # ground disc radius, in R
ground_panel_length_r = parse(Float64, get(ENV, "GROUND_PANEL_LENGTH_R", "0.15")) # triangle side, in R
# Below-ground particle policy (Ryan 2026-08-18 ruling 2; refined 2026-08-20
# after the p022_ige_coarse blow-up, phase_03_particle_policy.md):
#   none     : leave particles that sneak past the ground alone (v1)
#   cull     : GlobalBox trim at the ground plane
#   truncate : raise the truncation-cylinder FLOOR to coincide with the ground
#              (same deletion set as cull, expressed through the cylinder)
ground_particle_policy = lowercase(get(ENV, "GROUND_PARTICLE_POLICY", "none"))
ground_particle_policy in ("none", "cull", "truncate") ||
    error("Unknown GROUND_PARTICLE_POLICY=$(repr(ground_particle_policy)); use none, cull, or truncate")
# Near-ground vertical velocity cutoff (Ryan 2026-08-20): for particles within
# GROUND_DAMP_BAND_R*R above the ground, the GROUND-WARD axial velocity
# component is scaled by the linear factor f = d/band (d = height above the
# ground), reaching 0 at the plane; receding motion is never damped, so
# strays below the plane can recover. Applied to the stored particle U just
# before the Euler position update (driver-local propagate! overwrite).
# 0 disables (bit-identical to the Phase 1/2 runs). Start value R/10.
ground_damp_band_r = parse(Float64, get(ENV, "GROUND_DAMP_BAND_R", "0.0"))
ground_damp_band_r >= 0 || error("GROUND_DAMP_BAND_R must be >= 0")
# Radial extent of the truncation cylinder. In ground effect the wake becomes
# a radial wall jet, so the OGE default 1.5R would amputate it; default 3R
# with the ground on (explicit convergence axis of 022 Phase 3).
trunc_radius_r = parse(Float64, get(ENV, "TRUNC_RADIUS_R", ground_enable ? "3.0" : "1.5"))
# Block Gauss-Seidel outer-loop instrumentation (rotor <-> ground coupling).
# Uses solve!'s OWN verbose/history kwargs (src/FLOWPanel_solver.jl) — no new
# iteration machinery — surfaced via a driver-level solve_formulation! method.
gs_log       = parse(Bool,    get(ENV, "GS_LOG",       "true"))
gs_verbose   = parse(Bool,    get(ENV, "GS_VERBOSE",   "false"))
gs_max_outer = parse(Int,     get(ENV, "GS_MAX_OUTER", "50"))
gs_tol       = parse(Float64, get(ENV, "GS_TOL",       "1e-8"))
# --- Multi-rotor knobs (022 extension, Ryan's specs 2026-08-25) ---------------
# NROTORS=1 is the exact legacy single-rotor path (bit-identical construction,
# frames, monitors, and outputs). 2 = side-by-side pair, 4 = quadcopter square,
# both at ROTOR_SPACING_R center-to-center and counter-rotating in adjacent
# pairs (mirrored blade geometry + negated frame omega; see build_rotor below).
nrotors = parse(Int, get(ENV, "NROTORS", "1"))
nrotors in (1, 2, 4) || error("NROTORS must be 1, 2, or 4; got $(nrotors)")
rotor_spacing_r = parse(Float64, get(ENV, "ROTOR_SPACING_R", "2.7"))
rotor_directions_env = get(ENV, "ROTOR_DIRECTIONS", "")
# -----------------------------------------------------------------------------

# Ensure the run is long enough to cover ramp-up + hold + withdraw + a settle
# tail (so it never ends mid-withdrawal and the plateau ripple is observable).
schedule_revs = freestream_ramp_revs + freestream_hold_revs +
    freestream_withdraw_revs + settle_revs
required_revs = max(nrevs, schedule_revs)
n_steps = round(Int, nt * required_revs) + spinup_steps
t_range = range(0.0, step=dt, length=n_steps)

cp_outer=true
core_size_panel = parse(Float64, get(ENV, "CORE_SIZE_PANEL", get(ENV, "KERNELOFFSET_PANEL", string(R * 1e-10))))
core_size_targets = parse(Float64, get(ENV, "CORE_SIZE_TARGETS", get(ENV, "KERNELOFFSET_TARGETS", get(ENV, "CORE_SIZE", get(ENV, "KERNELOFFSET", "1e-3")))))
kernelcutoff = R * 1e-13
p_per_step = parse(Int, get(ENV, "P_PER_STEP", "2"))
overlap = parse(Float64, get(ENV, "OVERLAP", "3.0"))
particle_shedding = lowercase(get(ENV, "PARTICLE_SHEDDING", "overlap_pps"))
# BRAINSTORM/016 opt-in smooth surface-vorticity conversion. "legacy" (default)
# is the historical edge-jump conversion and is bit-identical to before this
# knob existed. "smooth" selects SurfaceVorticityConversion, which replaces
# method_trailing/method_unsteady entirely and REQUIRES unsteady_filament=true
# (supplying either method alongside it is a configuration error the wake
# constructor rejects). CONVERSION_SIGMA defaults to the same tip sigma the
# sigma_overlap trailing policy uses, so "smooth" sheds at the ladder's sigma.
conversion_mode = lowercase(get(ENV, "CONVERSION", "legacy"))
conversion_overlap = parse(Float64, get(ENV, "CONVERSION_OVERLAP", "1.3"))
conversion_sigma_env = get(ENV, "CONVERSION_SIGMA", "")
conversion_attribution = Symbol(lowercase(get(ENV, "ATTRIBUTION", "upstream")))
conversion_diagnose = lowercase(get(ENV, "CONVERSION_DIAGNOSE", "false")) == "true"
conversion_mode in ("legacy", "smooth") ||
    error("Unknown CONVERSION=$(repr(conversion_mode)); use legacy or smooth")
merge_r_factor = parse(Float64, get(ENV, "MERGE_R_FACTOR", "0.02"))
merge_r_hash_factor = 0.02
merge_sigma_relative = false
merge_particles = parse(Bool, get(ENV, "MERGE_PARTICLES", "true"))
# Das (first shed-row offset) = max(eta*dt*|V_te|, min_displacement), so eta is
# the fraction of one timestep's trailing-edge travel. Note the offset displaces
# the whole attached wake row downstream of the TE, but does NOT set the shed
# particle sigma: `_shed_particles!` spans node row 1 -> row 2, whose separation
# is one step's TE travel (the Das cancels between the rows).
init_Das_eta_kinematic = parse(Float64, get(ENV, "DAS_ETA_KINEMATIC", "0.2"))
# Chord-proportional first-row offset (BRAINSTORM/018, 2026-08-04): when
# finite, |Das| at each shedding station = DAS_CHORD_FRACTION * local chord,
# direction along the local kinematic TE tangent. This replaces the eta
# parameterization (Das = eta*dt*|V_te| ~ r), whose wedge-shaped Das/c_local
# distribution cannot place the whole span on the 014 plateau (0.25-1.5c) at
# any single eta. dt-independent by construction; the 0.4x spin-up freeze
# factor and the min_displacement floor become irrelevant. NaN = off (legacy).
das_chord_fraction = parse(Float64, get(ENV, "DAS_CHORD_FRACTION", "NaN"))
# Span-uniform handoff clearance (BRAINSTORM/018 phase_13 section 4b,
# 2026-08-05): when finite, per-station |Das|_j = D*sigma - (N-1)*|w x r_j|*dt
# with D = DAS_UNIFORM_DSIGMA, so the TE -> front-of-first-free-row distance
# d_front,j = |Das_j| + (N-1)*travel_j equals D*sigma IDENTICALLY at every
# station and timestep -- Das absorbs the dt dependence by construction (at
# N=1 there is no travel term: |Das| = D*sigma, span-uniform absolute).
# sigma is the SigmaOverlap shed sigma (tip_sigma_default). NaN = off.
# Mutually exclusive with DAS_CHORD_FRACTION and DAS_REFRESH.
das_uniform_dsigma = parse(Float64, get(ENV, "DAS_UNIFORM_DSIGMA", "NaN"))
# Chord–sigma co-scaling (BRAINSTORM/018 Phase 16, 2026-08-14): when
# SIGMA_CHORD_FRACTION = s* is finite, the shed sigma becomes PER-STATION,
# sigma_j = max(s* * c_local_j, SIGMA_FLOOR_R * R), replacing the span-uniform
# PARTICLE_SHEDDING sigma law (StationSigmaOverlap; OVERLAP still sets the
# particle count per filament). DAS_SIGMA_LAMBDA = lambda then sets
# |Das|_j = lambda * sigma_j, so Das/c AND Das/sigma are span-uniform
# simultaneously — the unique one-parameter family satisfying the 014 chord
# band, the Phase-12 clearance band, and dt-independence at once (see
# phase_16_chord_sigma_coscaling.md). NaN = off (both).
sigma_chord_fraction = parse(Float64, get(ENV, "SIGMA_CHORD_FRACTION", "NaN"))
sigma_floor_r = parse(Float64, get(ENV, "SIGMA_FLOOR_R", "0.0"))
das_sigma_lambda = parse(Float64, get(ENV, "DAS_SIGMA_LAMBDA", "NaN"))
# Phase 16 F1: curvature cap |Das|_j = min(lambda*sigma_j, beta*r_j/N) so the
# rigid extent never subtends more than beta radians of the local helix
# circle (theta_j = N*|Das|_j/r_j <= beta). NaN = off (pure co-scaling).
das_curvature_beta = parse(Float64, get(ENV, "DAS_CURVATURE_BETA", "NaN"))
# Phase 16 F1b (Route B): place the Das ENDPOINT on the frozen-wake path
# (backward swept arc + induced drift) at arc length lambda*sigma_j, instead
# of along the straight TE tangent. DAS_ARC_HELIX_SOURCE=steady reads the
# per-station induced-drift table DAS_ARC_TABLE (CSV written by the TE
# downwash probe, components in the local TE basis); =kinematic uses zero
# drift (pure backward arc).
das_arc_placed = parse(Bool, get(ENV, "DAS_ARC_PLACED", "false"))
das_arc_source = lowercase(get(ENV, "DAS_ARC_HELIX_SOURCE", "steady"))
das_arc_table = get(ENV, "DAS_ARC_TABLE", "")
set_Das_min_kinematic_displacement = parse(Float64, get(ENV, "DAS_MIN_DISPLACEMENT_R", "0.01")) * R
# Lay the offset along the arc the TE actually sweeps (true) rather than along
# the tangent to it (false, legacy). The two differ by ~r*θ²/2 with θ = eta*dψ:
# negligible at the default eta=0.2 (6e-4 R) but 0.22R at eta=4, where the
# tangent construction flings the first wake row outside the rotor disk.
set_Das_kinematic_arc = parse(Bool, get(ENV, "DAS_KINEMATIC_ARC", "true"))
# Re-derive Das from the *current* kinematic state each step (BRAINSTORM 014).
# Default false = legacy behavior: Das frozen at its t=0 magnitude, which is
# computed at the spin-up fraction (0.4x RPM), so eta_eff = 0.4*eta.
set_Das_refresh = parse(Bool, get(ENV, "DAS_REFRESH", "false"))
# Panel-wake rows between the Das row and the particle handoff (BRAINSTORM 014
# Proposal 1): row 1 is rigid (re-placed at TE+Das each step), rows 2..N convect
# freely, particles shed from rows N -> N+1. Handoff distance ~ Das+(N-1)*travel.
nwakerows = parse(Int, get(ENV, "NWAKEROWS", "1"))
# Shed the wake at the local induced velocity rather than the kinematic one.
shed_with_induced_velocity = parse(Bool, get(ENV, "SHED_WITH_INDUCED_VELOCITY", "false"))
p_correct_kuttacondition_flag = false
wake_core_size = parse(Float64, get(ENV, "WAKE_CORE_SIZE", string(core_size_targets)))
# wake_nu_default = 1.85508e-5 / rho # from FLOWUnsteady docs
wake_nu_default = 1.69e-5 / rho # from NASA paper
wake_nu = parse(Float64, get(ENV, "WAKE_NU", string(wake_nu_default)))
# Item 005 E4: multiplicative bump to the viscous core-spreading rate for the
# damping sweep (WAKE_NU_FACTOR=3, 10, ...); default 1.0 leaves wake_nu unchanged.
wake_nu *= parse(Float64, get(ENV, "WAKE_NU_FACTOR", "1.0"))
# wake_nu = parse(Float64, get(ENV, "WAKE_NU", "1.5e-5"))
wake_core_beta = parse(Float64, get(ENV, "WAKE_CORE_BETA", "1.5"))
run_kj = parse(Bool, get(ENV, "RUN_KJ", "false"))
lamb_only = parse(Bool, get(ENV, "LAMB_ONLY", "false"))
# Item 005 tuning: steady-Bernoulli-only monitor set. Skips both CG pressure
# solves and the per-step FMM Hessian, so iterations are far cheaper while still
# giving a CT history to read plateau ripple from.
bernoulli_only = parse(Bool, get(ENV, "BERNOULLI_ONLY", "false"))
run_monitors = parse(Bool, get(ENV, "RUN_MONITORS", "true"))
# Wake->body solve formulation (src/FLOWPanel_formulation.jl).
#   velocity : VelocityThroughSources() -- production default
#   green    : GreenReconstruction(...) -- reconstructs the wake body-trace from
#              sampled wake velocities; valid for a particle wake under FMM.
formulation_name = lowercase(get(ENV, "RHPC_FORMULATION", "velocity"))
green_recompute_interval = parse(Int, get(ENV, "GREEN_RECOMPUTE_INTERVAL", "1"))
green_gauge = Symbol(lowercase(get(ENV, "GREEN_GAUGE", "area_mean")))

read_path = joinpath(pnl.examples_path, "data")

# msh_file = joinpath(read_path, "phantom_3_rebuild_r2.msh")
# te_indices_1 = [9, 175, 127]
# te_indices_2 = [13, 286, 238]

# Mesh selection. RHPC_MESH names a known mesh key; RHPC_MESH_FILE overrides it
# with an explicit file (absolute, or relative to examples/data). The three legacy
# meshes carry hard-coded trailing-edge seed triples; every other mesh gets its
# seeds detected automatically by examples/dji9443_trailing_edge.jl.
rhpc_mesh = lowercase(get(ENV, "RHPC_MESH", "40_40"))
rhpc_mesh_file = get(ENV, "RHPC_MESH_FILE", "")

known_meshes = Dict(
    "40_40"      => "dji9443_new_40_40.msh",
    "56_57"      => "dji9443_56_57.msh",
    "80_81"      => "dji9443_80_81.msh",
    # Phase 2d production recipe: flat root cap + round CapUMinTess=4 tip cap.
    # Do NOT substitute a round-ct3 "capped" mesh — that recipe carries a tip
    # circulation artifact under Dirichlet.
    "45_145_ct4" => "dji9443_20260725_45_145_capped_captess4.msh",
    "45_185_ct4" => "dji9443_20260725_45_185_capped_captess4.msh",
)

# 0-based seed triples from the original hand inspection, converted to 1-based.
legacy_te_indices = Dict(
    "40_40" => ([1614, 1574, 45] .+ 1,    [3324, 3284, 1755] .+ 1),
    "56_57" => ([6370, 6314, 3255] .+ 1,  [3117, 3061, 0] .+ 1),
    "80_81" => ([12898, 12818, 6549] .+ 1, [6351, 6271, 3] .+ 1),
)

msh_file = if !isempty(rhpc_mesh_file)
    isabspath(rhpc_mesh_file) ? rhpc_mesh_file : joinpath(read_path, rhpc_mesh_file)
elseif haskey(known_meshes, rhpc_mesh)
    joinpath(read_path, known_meshes[rhpc_mesh])
else
    error("Unknown RHPC_MESH=$(repr(rhpc_mesh)); known keys: " *
          "$(join(sort(collect(keys(known_meshes))), ", ")). " *
          "Set RHPC_MESH_FILE to use another mesh.")
end
isfile(msh_file) || error("Mesh file does not exist: $(msh_file)")

if isempty(rhpc_mesh_file) && haskey(legacy_te_indices, rhpc_mesh)
    te_indices_1, te_indices_2 = legacy_te_indices[rhpc_mesh]
    te_seed_source = "hardcoded"
else
    include(joinpath(pnl.examples_path, "dji9443_trailing_edge.jl"))
    te_indices_1, te_indices_2 = find_dji9443_trailing_edge_indices(msh_file; watertight=true)
    te_seed_source = "auto"
end
println("Mesh: $(msh_file)")
println("  TE seeds ($(te_seed_source)): blade1=$(te_indices_1)  blade2=$(te_indices_2)")

axial_dimension = occursin("dji9443", msh_file) ? 1 : 2 # DJI9443 geometry is rotated compared to typical rotor convention
radial_dimension = occursin("dji9443", msh_file) ? 2 : 1 # this might be wrong for non-dji9443

Vinf_direction = occursin("dji9443", msh_file) ? [cosd(AOA), sind(AOA), 0.0] : [0.0, -cosd(AOA), sind(AOA)]

# --- Multi-rotor layout (022 extension) ---------------------------------------
# The mirror dimension is the lateral axis that is neither axial nor radial
# (dims are a permutation of 1:3): negating it flips blade handedness while
# leaving every axial and radial coordinate -- and hence the TE seed indices,
# root-clip radii, and shedding bboxes -- untouched. A mirrored blade spun with
# negated omega thrusts in the same axial direction.
mirror_dimension = 6 - axial_dimension - radial_dimension
rotor_centers = if nrotors == 1
    [zeros(3)]
elseif nrotors == 2
    s = rotor_spacing_r * R
    [begin c = zeros(3); c[radial_dimension] = sgn * s / 2; c end
     for sgn in (1.0, -1.0)]
else  # 4: square in the lateral plane, ordered around the perimeter
    s = rotor_spacing_r * R
    [begin c = zeros(3); c[radial_dimension] = sa * s / 2; c[mirror_dimension] = sb * s / 2; c end
     for (sa, sb) in ((1.0, 1.0), (-1.0, 1.0), (-1.0, -1.0), (1.0, -1.0))]
end
# Counter-rotating pairs: 2r = (+1, -1); 4r = quadcopter alternation, adjacent
# rotors opposite, diagonal pairs co-rotating (dir = sign(a)*sign(b)).
rotor_directions = if !isempty(rotor_directions_env)
    dirs = parse.(Int, split(rotor_directions_env, ','))
    length(dirs) == nrotors || error("ROTOR_DIRECTIONS has $(length(dirs)) " *
        "entries; NROTORS=$(nrotors)")
    all(d -> d in (-1, 1), dirs) || error("ROTOR_DIRECTIONS entries must be +-1")
    dirs
elseif nrotors == 1
    [1]
elseif nrotors == 2
    [1, -1]
else
    [Int(sign(c[radial_dimension]) * sign(c[mirror_dimension])) for c in rotor_centers]
end
max_rotor_offset = nrotors == 1 ? 0.0 : maximum(norm(c) for c in rotor_centers)
# Multi-rotor supports the production carrier only; the exotic Das/sigma laws,
# Laplace-family monitors, and warm start are single-rotor machinery.
if nrotors > 1
    bernoulli_only || error("NROTORS>1 requires BERNOULLI_ONLY=true " *
        "(PressureLaplace monitors are single-body)")
    (run_kj || lamb_only) && error("NROTORS>1 is incompatible with RUN_KJ / LAMB_ONLY")
    isnan(sigma_chord_fraction) || error("NROTORS>1 does not support SIGMA_CHORD_FRACTION")
    isnan(das_chord_fraction) || error("NROTORS>1 does not support DAS_CHORD_FRACTION")
    das_arc_placed && error("NROTORS>1 does not support DAS_ARC_PLACED")
    conversion_mode == "legacy" || error("NROTORS>1 requires CONVERSION=legacy")
    formulation_name == "velocity" || error("NROTORS>1 requires RHPC_FORMULATION=velocity")
end
# -----------------------------------------------------------------------------

msh = pnl.read_gmsh(msh_file)
nodes0, cells0 = pnl.meshes2nodes_cells(msh)
nodes0 .*= R / maximum(nodes0[radial_dimension, :])

kernel = Union{pnl.ConstantSource, pnl.VortexRing}
DBC = kernel == pnl.VortexRing ? false : true

0.0 <= shedding_r_over_R <= 1.0 || error("shedding_r_over_R must be between 0 and 1")

function make_shedding_bbox(nodes, seed_nodes, radial_dimension, R, shedding_r_over_R)
    radial_midpoint = sum(nodes[radial_dimension, seed_nodes]) / length(seed_nodes)
    radial_sign = sign(radial_midpoint)
    radial_sign == 0 && error("Seed edge lies on the rotor axis; cannot determine shedding side")

    lower = [minimum(nodes[i, :]) for i in 1:size(nodes, 1)]
    upper = [maximum(nodes[i, :]) for i in 1:size(nodes, 1)]
    padding = max(sqrt(eps(eltype(nodes))) * R, R * 1e-6)
    lower .-= padding
    upper .+= padding

    radial_cutoff = shedding_r_over_R * R
    if radial_sign > 0
        lower[radial_dimension] = radial_cutoff - padding
    else
        upper[radial_dimension] = -radial_cutoff + padding
    end

    return (pnl.SVector{3}(lower...), pnl.SVector{3}(upper...))
end

# save vtk file to inspect for shedding nodes
# println("Saving initial at $(save_path)...")
# vtk_path = joinpath(save_path, "rotor_initial")
# pnl.write_vtk(vtk_path, rotor)
# sherlock

# Anchor the trailing-edge walk at BOTH ends. `find_dji9443_trailing_edge_indices`
# already identifies the complete TE chain intrinsically -- edges that are
# predominantly radial AND sharp -- validates it as a simple open path, and
# returns `[outer, second_outer, inner]`. Passing `inner` as `end_node` makes
# `calc_shedding_from_seed` trace exactly that chain and ERROR if it cannot
# reach it.
#
# This replaces an earlier, wrong approach. `shedding_r_over_R` had been doing
# three unrelated jobs at once: separating the two blades, clipping the root as
# a modeling choice, and -- only by coincidence -- stopping the walk before it
# wrapped onto the ROOT CAP. Any mesh that moves the blade root breaks that
# coincidence, and it failed SILENTLY in both directions: too small and the
# chain wrapped the whole cap perimeter (measured on the 0.15R hub mesh: 136
# shedding edges of which 92 were cap edges; that run diverged in 17 steps,
# BRAINSTORM/018 phase_05 job 13031568), too large and it quietly discarded real
# trailing edge (raising the cutoff to root+0.05 dropped 3 genuine TE edges at
# the hub). With `end_node` the walk is anchored by geometry instead, and a
# wrong turn is loud rather than silent.
#
# `shedding_r_over_R` is now ONLY what its name says: an explicit modeling root
# clip, plus blade separation via the bbox half-space. It is deliberately left
# at 0.1 for comparability with the whole 018 campaign on the stock blade
# (which sheds from r/R 0.111, not from the true root 0.0095). On any hub
# variant the blade root is OUTBOARD of 0.1, so the clip is inert and the blade
# sheds its full trailing edge down to the root -- which is the desired
# behavior, and needs no special-casing.
# The two jobs must be separated, because they conflict: the modeling clip
# deliberately truncates the TE on the stock blade (whose true root TE node sits
# at r/R 0.0095, inboard of the 0.1 clip), so a walk confined by the clip can
# never reach `end_node`. So: TRACE the full trailing edge, anchored end-to-end
# and therefore cap-free, then CLIP it afterwards as an explicit modeling choice.
# The bbox now only separates the two blades (cutoff 0).
te_end_1 = length(te_indices_1) >= 3 ? te_indices_1[3] : nothing  # legacy hardcoded
te_end_2 = length(te_indices_2) >= 3 ? te_indices_2[3] : nothing  # seeds give only 2

# Modeling root clip, applied to the traced chain. Identical criterion to the old
# bbox test (edge MIDPOINT radius), so the retained set -- and hence stock
# behavior -- is unchanged; it is only applied after tracing instead of during.
function clip_shedding_root(nodes, shedding, cells, radial_dimension, R, clip_r_over_R)
    keep = Int[]
    for j in axes(shedding, 2)
        p, nia, nib = shedding[1, j], shedding[2, j], shedding[3, j]
        na, nb = cells[nia, p], cells[nib, p]
        mid = (nodes[radial_dimension, na] + nodes[radial_dimension, nb]) / 2
        abs(mid) / R >= clip_r_over_R && push!(keep, j)
    end
    return shedding[:, keep]
end

# Regression net: no shed edge may run CIRCUMFERENTIALLY. A trailing-edge edge
# runs radially (spanwise); a root-cap perimeter edge runs around the cap, i.e.
# perpendicular to radial. Testing edge DIRECTION rather than edge RADIUS is
# what makes this work now that the chain legitimately reaches the blade root --
# a radius test could not tell a root TE edge from a cap edge. Threshold is
# loose (0.3 vs the detector's 0.7) because this only needs to catch cap wrap
# (ratio ~0), not re-do the detector's selection.
function count_circumferential_shedding(nodes, shedding, cells, radial_dimension; tol=0.3)
    n = 0
    for j in axes(shedding, 2)
        p, nia, nib = shedding[1, j], shedding[2, j], shedding[3, j]
        na, nb = cells[nia, p], cells[nib, p]
        delta = nodes[:, na] - nodes[:, nb]
        len = sqrt(sum(abs2, delta))
        len > 0 && abs(delta[radial_dimension]) / len < tol && (n += 1)
    end
    return n
end

function shedding_root_r_over_R(nodes, shedding, cells, radial_dimension, R)
    isempty(shedding) && return NaN
    root_edge = shedding[:, end]
    pi, nia, nib = root_edge[1], root_edge[2], root_edge[3]
    edge_nodes = cells[[nia, nib], pi]
    midpoint = (nodes[:, edge_nodes[1]] + nodes[:, edge_nodes[2]]) / 2
    return midpoint[radial_dimension] / R
end

# --- Per-rotor build (022 multi-rotor extension) ------------------------------
# The CLAUDE.md shedding invariant is preserved PER ROTOR: construct a
# noshedding body first (the constructor re-winds cells), trace the trailing
# edge from THAT body's nodes/cells, then rebuild with the shedding.
# direction == -1 mirrors the lateral (non-axial, non-radial) coordinate BEFORE
# the noshedding build, flipping blade handedness while leaving the TE seed
# indices and all radial bookkeeping valid; the constructor's outward-winding
# pass re-orients the reflected cells. Translation to `center` happens only
# AFTER the shedding is traced (a rigid translation is winding- and
# shedding-invariant), so every radius/chord computation stays rotor-local --
# downstream station math must use `local_nodes`, never the translated body's.
function build_rotor(i, center, direction)
    lab = nrotors == 1 ? "" : "rotor$(i) "
    local_nodes0 = copy(nodes0)
    direction == -1 && (local_nodes0[mirror_dimension, :] .*= -1)
    # copy(cells0) is LOAD-BEARING: the constructor re-winds `cells` IN PLACE
    # and may store the array without copying, so sharing one cells matrix
    # across builds lets a later (mirrored) build flip an earlier rotor's
    # already-constructed cells behind its back (measured: stock rotor gamma
    # 2.4x off, CT collapsed 250x, with zero rotor-rotor interaction).
    base = pnl.RigidWakeBody{kernel}(local_nodes0, copy(cells0), pnl.noshedding;
        core_size=core_size_panel, core_size_panel, core_size_targets, kernelcutoff,
        semiinfinite_wake=false, watertight=true, DBC)

    blade_root_r_over_R = minimum(abs.(base.nodes[radial_dimension, :])) / R

    bbox1 = make_shedding_bbox(base.nodes, te_indices_1[1:2], radial_dimension, R, 0.0)
    shedding1_full = pnl.calc_shedding_from_seed(base.nodes, base.cells, te_indices_1[1], te_indices_1[2];
        bbox=bbox1, end_node=te_end_1, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
    bbox2 = make_shedding_bbox(base.nodes, te_indices_2[1:2], radial_dimension, R, 0.0)
    shedding2_full = pnl.calc_shedding_from_seed(base.nodes, base.cells, te_indices_2[1], te_indices_2[2];
        bbox=bbox2, end_node=te_end_2, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

    shedding1 = clip_shedding_root(base.nodes, shedding1_full, base.cells, radial_dimension, R, shedding_r_over_R)
    shedding2 = clip_shedding_root(base.nodes, shedding2_full, base.cells, radial_dimension, R, shedding_r_over_R)
    println("$(lab)Trailing edge traced: blade1 $(size(shedding1_full, 2)) edges -> $(size(shedding1, 2)) after root clip at r/R $(shedding_r_over_R) (blade root r/R $(round(blade_root_r_over_R, digits=4)))")

    for (k, shed) in enumerate((shedding1, shedding2))
        nbad = count_circumferential_shedding(base.nodes, shed, base.cells, radial_dimension)
        nbad == 0 || error("$(lab)shedding$(k) includes $(nbad) of $(size(shed, 2)) circumferential " *
            "edges -- the trailing-edge chain has wrapped onto the ROOT CAP. The walk should " *
            "be anchored by end_node (te_indices[3]); check that the TE detector returned it.")
    end

    final_nodes = all(iszero, center) ? base.nodes : base.nodes .+ center
    body = pnl.RigidWakeBody{kernel}(final_nodes, base.cells, [shedding1, shedding2];
        core_size=core_size_panel, core_size_panel, core_size_targets, kernelcutoff,
        semiinfinite_wake=false, watertight=true,
        ensure_winding=true, DBC)

    println("$(lab)Requested shedding root at |r/R| >= $(shedding_r_over_R)")
    println("  shedding1 root midpoint r/R = $(shedding_root_r_over_R(base.nodes, shedding1, base.cells, radial_dimension, R))")
    println("  shedding2 root midpoint r/R = $(shedding_root_r_over_R(base.nodes, shedding2, base.cells, radial_dimension, R))")
    nrotors > 1 && println("  $(lab)center/R = $(round.(center ./ R, digits=4)), " *
        "direction = $(direction > 0 ? "+1 (stock)" : "-1 (mirrored, reversed)")")

    return (; body, local_nodes=base.nodes, sheddings=(shedding1, shedding2), center, direction)
end

rotor_builds = [build_rotor(i, rotor_centers[i], rotor_directions[i]) for i in 1:nrotors]
rotors = [rb.body for rb in rotor_builds]
rotor = rotors[1]   # legacy alias: single-rotor code paths below use it directly

# Diagnostic gate: ablate the clipping_backscatter strategy by replacing the
# canonical SFS with an equivalent DynamicSFS that omits the clip. Default
# off, leaves the run bit-identical to the canonical configuration.
sfs_no_backscatter_clip    = parse(Bool,    get(ENV, "SFS_NO_BACKSCATTER_CLIP",    "false"))
sfs_no_backscatter_project = parse(Bool,    get(ENV, "SFS_NO_BACKSCATTER_PROJECT", "false"))
sfs_backscatter_signed     = parse(Bool,    get(ENV, "SFS_BACKSCATTER_SIGNED",     "false"))
sfs_magnitude_control      = parse(Bool,    get(ENV, "SFS_MAGNITUDE_CONTROL",      "false"))
sfs_directional_control    = parse(Bool,    get(ENV, "SFS_DIRECTIONAL_CONTROL",    "false"))
sfs_threelevel             = parse(Bool,    get(ENV, "SFS_THREELEVEL",             "false"))
sfs_nostatic               = parse(Bool,    get(ENV, "SFS_NOSTATIC",               "false"))
sfs_maxC                   = parse(Float64, get(ENV, "SFS_MAXC",                   "1.0"))
sfs_rlxf                   = parse(Float64, get(ENV, "SFS_RLXF",                   "0.005"))
panel_wake_hessian_to_particles = parse(Bool, get(ENV, "PANEL_WAKE_HESSIAN_TO_PARTICLES",
    string(!parse(Bool, get(ENV, "WAKEROW_NO_HESSIAN_TO_PARTICLES", "true")))))
wakerow_no_hessian_to_particles = !panel_wake_hessian_to_particles
body_hessian_to_particles       = parse(Bool, get(ENV, "BODY_HESSIAN_TO_PARTICLES",        "false"))
# Split regularization for the body->particle influence: evaluate the velocity
# GRADIENT with this (larger) kernel offset while the velocity keeps
# CORE_SIZE_TARGETS. NaN disables (single-pass, shared offset). Only active
# when BODY_HESSIAN_TO_PARTICLES=true. Not strictly physical — smooths the
# |∇U| bumpiness of piecewise-constant doublet panels felt by nearby particles.
body_gradient_core_size      = parse(Float64, get(ENV, "BODY_GRADIENT_CORE_SIZE",    get(ENV, "BODY_GRADIENT_KERNELOFFSET", "NaN")))
body_on_wake                    = parse(Bool, get(ENV, "BODY_ON_WAKE",                     "true"))
panel_wake_on_particles         = parse(Bool, get(ENV, "PANEL_WAKE_VELOCITY_TO_PARTICLES",
    get(ENV, "PANEL_WAKE_ON_PARTICLES", "true")))
particle_hessian_self           = parse(Bool, get(ENV, "PARTICLE_HESSIAN_SELF",            "true"))
particle_relax                  = parse(Bool, get(ENV, "PARTICLE_RELAX",                   "true"))
bound_strength_rlx              = parse(Float64, get(ENV, "BOUND_STRENGTH_RLX",            "1.0"))  # E4.8 body-strength low-pass (1.0=off)
diagnose_particle_gamma         = parse(Bool, get(ENV, "DIAGNOSE_PARTICLE_GAMMA",          "false"))
diagnose_particle_influence     = parse(Bool, get(ENV, "DIAGNOSE_PARTICLE_INFLUENCE",      "false"))
particle_diagnostic_vertical    = ntuple(i -> i == axial_dimension ? 1.0 : 0.0, 3)
sfs_off                         = parse(Bool, get(ENV, "SFS_OFF",                          "false"))

# Item 006: optional spatially-filtered relaxation. When RELAX_FILTER_DOWNSTREAM_R is
# set, apply Pedrizzetti relaxation only to particles that have propagated at least
# RELAX_FILTER_DOWNSTREAM_R*R downstream (+axial) of the rotor plane, leaving the
# near-rotor band unrelaxed. Unset/NaN => unfiltered full-wake relaxation (default).
relax_filter_downstream_R = parse(Float64, get(ENV, "RELAX_FILTER_DOWNSTREAM_R", "NaN"))
# Relaxation factor (rlxf) of the corrected-Pedrizzetti scheme. Defaults to the FLOWVPM
# stock value so unset behavior is unchanged; override with RELAX_RLXF.
stock_relaxation = pnl.FLOWVPM.relaxation_correctedpedrizzetti
relax_rlxf = parse(Float64, get(ENV, "RELAX_RLXF", string(stock_relaxation.rlxf)))
base_relaxation = pnl.FLOWVPM.Relaxation(stock_relaxation.relax,
    stock_relaxation.nsteps_relax, relax_rlxf)
relaxation_scheme = if isnan(relax_filter_downstream_R)
    base_relaxation
else
    d = relax_filter_downstream_R * R
    plane_point  = SVector{3,Float64}(ntuple(i -> i == axial_dimension ? d   : 0.0, 3))
    plane_normal = SVector{3,Float64}(ntuple(i -> i == axial_dimension ? 1.0 : 0.0, 3))
    pnl.plane_filtered_relaxation(base_relaxation, plane_point, plane_normal; i_frame=1)
end

sfs_choice = if sfs_off
    pnl.FLOWVPM.noSFS
elseif sfs_backscatter_signed
    pnl.FLOWVPM.SFS_Cd_twolevel_backscatter_signed
elseif sfs_no_backscatter_project
    pnl.FLOWVPM.SFS_Cd_twolevel_nobackscatter_projection
elseif sfs_threelevel
    pnl.FLOWVPM.SFS_Cd_threelevel_nobackscatter
else
    sfs_controls = ()
    sfs_magnitude_control   && (sfs_controls = (sfs_controls..., pnl.FLOWVPM.control_magnitude))
    sfs_directional_control && (sfs_controls = (sfs_controls..., pnl.FLOWVPM.control_directional))
    sfs_clippings = sfs_no_backscatter_clip ? () : (pnl.FLOWVPM.clipping_backscatter,)
    sfs_model = sfs_nostatic ?
        ((pfield; optargs...) -> pnl.FLOWVPM.E_nostaticparticles(pfield; E=pnl.FLOWVPM.Estr_fmm, optargs...)) :
        pnl.FLOWVPM.Estr_fmm
    pnl.FLOWVPM.DynamicSFS(sfs_model,
        pnl.FLOWVPM.pseudo3level_beforeUJ,
        pnl.FLOWVPM.pseudo3level_positive_afterUJ;
        alpha=0.999, clippings=sfs_clippings, controls=sfs_controls,
        maxC=sfs_maxC, rlxf=sfs_rlxf)
end

# wake_rotor = pnl.PanelWake(rotor; nwakerows=12, core_size=wake_core_size)
FV = pnl.FLOWVPM

# Local chord at each Das STATION. Das/velocity_te are vertex-based: station
# j (1..nshed) is edge j's nib node, station nshed+1 is the last edge's nia
# node (`_kinematic_velocity_te!`, src/FLOWPanel_frames.jl). Chord = max
# distance from the station's TE node to any same-blade node within a radius
# band (one ring spacing) — the same 3D LE-TE section measure
# `scripts/p018_mesh_profile.py` reports for interior rings; chord varies
# slowly spanwise so band smoothing is benign. (Defined here, before the
# shedding-method construction, because Phase 16's per-station sigma needs it.)
function station_chords(nodes, shedding, cells, radial_dimension; band)
    nshed = size(shedding, 2)
    station_nodes = Vector{Int}(undef, nshed + 1)
    for j in 1:nshed
        station_nodes[j] = cells[shedding[3, j], shedding[1, j]]   # nib
    end
    station_nodes[nshed + 1] = cells[shedding[2, nshed], shedding[1, nshed]]  # final nia
    chords = zeros(nshed + 1)
    for (j, node_idx) in enumerate(station_nodes)
        te = nodes[:, node_idx]
        rj = te[radial_dimension]
        best = 0.0
        for k in axes(nodes, 2)
            rk = nodes[radial_dimension, k]
            sign(rk) == sign(rj) || continue      # same blade only
            abs(rk - rj) <= band || continue
            d2 = (nodes[1, k] - te[1])^2 + (nodes[2, k] - te[2])^2 +
                 (nodes[3, k] - te[3])^2
            d2 > best && (best = d2)
        end
        best > 0 || error("station_chords: no nodes within band at station $(j)")
        chords[j] = sqrt(best)
    end
    return chords
end

# Per-station TE radii along a shedding edge set, vertex-based to match
# Das/velocity_te (station j = edge j's nib node, station nshed+1 = last nia).
function station_radii(nodes, shedding, cells, radial_dimension)
    nshed = size(shedding, 2)
    r = zeros(nshed + 1)
    for j in 1:nshed
        r[j] = abs(nodes[radial_dimension, cells[shedding[3, j], shedding[1, j]]])
    end
    r[nshed + 1] = abs(nodes[radial_dimension, cells[shedding[2, nshed], shedding[1, nshed]]])
    return r
end

method_trailing = if particle_shedding == "overlap_pps"
    pnl.OverlapPPS(overlap, p_per_step)
elseif particle_shedding == "sigma_overlap"
    tip_sigma = 2 * pi * R / nt * overlap / p_per_step
    pnl.SigmaOverlap(tip_sigma, overlap)
else
    error("Unknown PARTICLE_SHEDDING=$(repr(particle_shedding)); use overlap_pps or sigma_overlap")
end

# Phase 16 chord–sigma co-scaling: per-station sigma_j = s* * c_local_j
# (floored at SIGMA_FLOOR_R * R) replaces the span-uniform sigma law above.
# The floor, where it binds, re-introduces sigma/c non-uniformity outboard of
# the crossing station — reported below, not glossed.
station_sigmas = nothing
if !isnan(sigma_chord_fraction)
    particle_shedding == "sigma_overlap" || error(
        "SIGMA_CHORD_FRACTION replaces the shed-sigma law and expects " *
        "PARTICLE_SHEDDING=sigma_overlap (got $(repr(particle_shedding)))")
    conversion_mode == "smooth" && error(
        "SIGMA_CHORD_FRACTION is incompatible with CONVERSION=smooth " *
        "(SurfaceVorticityConversion owns shedding with a single scalar sigma)")
    band_cs = 0.02 * R  # ~ one spanwise ring spacing on the 45-ring blade
    station_sigmas = [
        max.(sigma_chord_fraction .*
                station_chords(rotor.nodes, shed, rotor.cells, radial_dimension;
                               band=band_cs),
             sigma_floor_r * R)
        for shed in rotor.shedding]
    method_trailing = pnl.StationSigmaOverlap(station_sigmas, overlap)
    for (k, sig) in enumerate(station_sigmas)
        n_floored = count(s -> s == sigma_floor_r * R, sig)
        println("Sigma chord mode: shedding$(k) sigma/R range " *
            "$(round(minimum(sig)/R, digits=4))-$(round(maximum(sig)/R, digits=4)) " *
            "(s* = $(sigma_chord_fraction), floor $(sigma_floor_r)R binds at " *
            "$(n_floored)/$(length(sig)) stations)")
    end
end

tip_sigma_default = 2 * pi * R / nt * overlap / p_per_step
conversion_sigma = isempty(conversion_sigma_env) ? tip_sigma_default :
    parse(Float64, conversion_sigma_env)

# Viscous scheme (BRAINSTORM 018, 2026-08-03 finding). The CoreSpreading
# passed here was silently inert in every prior run: the pfield declared
# integration=rungekutta3 while FLOWPanel steps with _euler, so
# viscousdiffusion no-opped (no spreading, no beta resets). That src mismatch
# is now fixed (PanelParticleWake declares integration=euler), which would
# have turned CoreSpreading ON and changed physics under every existing case
# tag. To keep default runs bit-identical to the campaign baseline, the de
# facto behavior (inviscid wake) is now the EXPLICIT default; set
# CORE_SPREADING_ACTIVE=true to get working core spreading. sgm0 (the reset
# core size and beta-trigger reference) defaults to the shed sigma, not
# wake_core_size=1e-3: resets stamp EVERY particle to sgm0, so it must be the
# ladder's sigma or a reset would collapse the resolution to 1e-3 m.
core_spreading_active = parse(Bool, get(ENV, "CORE_SPREADING_ACTIVE", "false"))
core_spreading_sgm0 = parse(Float64,
    get(ENV, "CORE_SPREADING_SGM0", string(tip_sigma_default)))
# iterror=false: an RBF (beta-reset) solve that misses tolerance must WARN,
# not kill a 40 h run — the projection over heavily-overlapped particles is
# known to be ill-conditioned (BRAINSTORM 008/M6), so best-effort resets are
# accepted and disclosed. itmax raised from the FLOWVPM default 15.
core_spreading_itmax = parse(Int, get(ENV, "CORE_SPREADING_ITMAX", "50"))
core_spreading_tol = parse(Float64, get(ENV, "CORE_SPREADING_TOL", "1e-3"))
core_spreading_iterror = parse(Bool, get(ENV, "CORE_SPREADING_ITERROR", "false"))
viscous_scheme = core_spreading_active ?
    pnl.FLOWVPM.CoreSpreading(wake_nu, core_spreading_sgm0,
        pnl.FLOWVPM.zeta_fmm; beta=wake_core_beta,
        itmax=core_spreading_itmax, tol=core_spreading_tol,
        iterror=core_spreading_iterror, verbose=true) :
    pnl.FLOWVPM.Inviscid()
# Phase 16 guard: a CoreSpreading beta-reset stamps EVERY particle to the
# scalar sgm0, which would destroy a per-station sigma distribution in one
# event. Production runs spreading-only (beta=1e9, resets unreachable); any
# finite-beta viscous config is incompatible with chord–sigma co-scaling.
!isnan(sigma_chord_fraction) && core_spreading_active && wake_core_beta < 1e6 &&
    error("SIGMA_CHORD_FRACTION requires spreading-only viscosity " *
          "(WAKE_CORE_BETA >= 1e6): a beta-reset would stamp the uniform " *
          "sgm0 = $(core_spreading_sgm0) over the per-station sigmas")

# Frozen-gradient geometric local integrator, BRAINSTORM 020 Phase 2R.
# Default false = stock forward Euler (bit-identical off-state).
wake_expint = parse(Bool, get(ENV, "WAKE_EXPINT", "false"))

# The two conversions need mutually exclusive wake options, so build the
# differing kwargs once rather than duplicating the constructor call.
# legacy: method_trailing/method_unsteady drive shedding, no unsteady filament.
# smooth: SurfaceVorticityConversion owns shedding; passing either method is a
#         configuration error, and unsteady_filament MUST be true.
conversion_kwargs = if conversion_mode == "smooth"
    (conversion = pnl.SurfaceVorticityConversion(conversion_sigma;
         overlap = conversion_overlap,
         attribution = conversion_attribution,
         diagnose_nearfield = conversion_diagnose),
     unsteady_filament = true)
else
    (method_trailing = method_trailing,
     method_unsteady = pnl.NoShed(),
     unsteady_filament = false)  # should be false if method_unsteady is NoShed
end

# --- BRAINSTORM/022: ground plane geometry ------------------------------------
# Built from the CONSTRUCTED rotor's nodes (post-rebuild), before the wake so
# the optional cull policy below knows the plane location. Downstream is the
# +axial direction (the truncation cylinder extrudes +axial); the ground disc
# sits GROUND_H_R*R downstream of the rotor plane with its normal pointing
# back at the rotor. Rotor-plane axial station = mean axial node coordinate
# (blades are thin; the value is printed for audit).
axial_unit = [i == axial_dimension ? 1.0 : 0.0 for i in 1:3]
x_rotor_plane = sum(rotor.nodes[axial_dimension, :]) / size(rotor.nodes, 2)
ground_x = x_rotor_plane + ground_h_r * R
# Multi-rotor: the disc (and the truncation cylinder below) keep their
# per-rotor margins by growing with the layout footprint -- effective radius =
# nominal radius + the largest rotor-center offset (0 at NROTORS=1).
ground_radius_eff = ground_radius_r * R + max_rotor_offset
ground, solver_ground = if ground_enable
    ground_center = axial_unit .* ground_x
    ground_normal = -axial_unit                      # toward the rotor
    g = pnl.FlatGround(ground_center, ground_normal, ground_radius_eff;
        panel_length=ground_panel_length_r * R)
    println("Ground plane: x_rotor_plane/R = $(round(x_rotor_plane / R, digits=4)), " *
        "ground at axial/R = $(round(ground_x / R, digits=4)) (h/R = $(ground_h_r)), " *
        "disc radius $(ground_radius_r)R" *
        (nrotors > 1 ? " (+$(round(max_rotor_offset / R, digits=3))R layout offset -> " *
            "$(round(ground_radius_eff / R, digits=3))R effective)" : "") *
        ", panel_length $(ground_panel_length_r)R " *
        "-> $(g.ncells) panels, $(g.nnodes) nodes")
    g, pnl.FlatGroundSolver(g)
else
    println("Ground plane DISABLED (OGE baseline); rotor plane axial/R = " *
        "$(round(x_rotor_plane / R, digits=4))")
    nothing, nothing
end

# Trim policies: truncation cylinder (radial extent widened in ground effect,
# see TRUNC_RADIUS_R above) + optional below-ground cull. GROUND_PARTICLE_POLICY
# =none (default) deliberately leaves below-ground particles alone.
# policy=truncate raises the cylinder FLOOR to the ground plane: the cylinder
# starts 0.5R upstream of the origin, so its length becomes ground_x + 0.5R.
effective_cylinder_depth = (ground_enable && ground_particle_policy == "truncate") ?
    (ground_x + 0.5R) : cylinder_depth
# The cylinder is ALSO the horizontal domain truncation: particles past
# TRUNC_RADIUS_R are deleted every step, so keep it INSIDE the paneled ground
# disc or the wall jet convects over unpaneled ground before deletion.
ground_enable && trunc_radius_r > ground_radius_r && @warn(
    "TRUNC_RADIUS_R ($(trunc_radius_r)R) exceeds GROUND_RADIUS_R " *
    "($(ground_radius_r)R): particles will travel beyond the paneled ground " *
    "support before the truncation cylinder deletes them. Widen the disc or " *
    "shrink the truncation radius.")
trim_policies = (pnl.GlobalCylinder(-0.5R .* axial_unit,
    effective_cylinder_depth .* axial_unit, trunc_radius_r * R + max_rotor_offset),)
if ground_enable && ground_particle_policy == "cull"
    box_lo = [-1e3 * R, -1e3 * R, -1e3 * R]
    box_hi = [1e3 * R, 1e3 * R, 1e3 * R]
    box_hi[axial_dimension] = ground_x
    trim_policies = (trim_policies..., pnl.GlobalBox(box_lo, box_hi))
end
# -----------------------------------------------------------------------------

wakes_rotor = [pnl.PanelParticleWake(r;
    nwakerows, max_particles=500_000, core_size=wake_core_size,
    particle_core_size=core_size_targets,
    viscous=viscous_scheme,
    SFS=sfs_choice,
    relaxation=relaxation_scheme,
    expint=wake_expint,
    conversion_kwargs...,
    shed_with_induced_velocity,
    particle_maintenance=pnl.ParticleMaintenance((
            trim_policies...,
            pnl.MergeParticles(;
                every=merge_particles ? 1 : 0,
                r=merge_sigma_relative ? merge_r_factor : merge_r_factor * R,
                r_hash=merge_sigma_relative ? merge_r_hash_factor : merge_r_hash_factor * R,
                sigma_relative=merge_sigma_relative),
        ))
    ) for r in rotors]
wake_rotor = wakes_rotor[1]   # legacy alias

smoothstep(x) = x <= 0 ? zero(x) : x >= 1 ? one(x) : x * x * (3 - 2 * x)

# Item 005: four-phase smoothstep freestream pulse (see parameter block above).
function magVinf_pulse(t)
    if t <= t_ramp_up
        return magVinf_start + (magVinf_peak - magVinf_start) * smoothstep(t_ramp_up > 0 ? t / t_ramp_up : 1.0)
    elseif t <= t_ramp_up + t_hold
        return magVinf_peak
    elseif t <= t_ramp_up + t_hold + t_withdraw
        s = t_withdraw > 0 ? (t - t_ramp_up - t_hold) / t_withdraw : 1.0
        return magVinf_peak + (magVinf_end - magVinf_peak) * smoothstep(s)
    else
        return magVinf_end
    end
end
Uinf(t) = magVinf_pulse(t) * Vinf_direction

function spinup_fraction(t)
    spinup_revs <= 0 && return 1.0
    return spinup_start_fraction + (1 - spinup_start_fraction) * smoothstep(t / spinup_duration)
end

ω_full = 2 * pi * RPM / 60
function set_frame_omega!(frames, i_frame, ω)
    frame = frames[i_frame]
    frames[i_frame] = typeof(frame)(
        frame.x, frame.v, frame.ω_axis, ω, frame.R, frame.Rp2g,
        frame.name, frame.parent_index, frame.child_index, frame.dependent_index,
    )
end

rotor_omega_axis = occursin("dji9443", msh_file) ? SVector{3}(-1.0, 0.0, 0.0) : SVector{3}(0.0, 1.0, 0.0)
# NROTORS=1: the legacy single rotating root frame, bit-identical to before.
# NROTORS>1: the root is INERT (omega = 0, v = 0) and each rotor gets its own
# child frame at its center with omega = direction * omega_full. add_frame!
# registers each rotor as a dependent of the root too, but an inert frame
# contributes zero kinematic velocity, so there is no co-rotation (the
# phase_00 trap requires a SPINNING ancestor). The ground stays in no frame.
frames = if nrotors == 1
    pnl.ReferenceFrame(rotor;
        origin=SVector{3}(0.0, 0.0, 0.0),
        v=SVector{3}(0.0, 0.0, 0.0),
        ω_axis=rotor_omega_axis,
        ω=ω_full * spinup_fraction(t_range[1]),
        R=SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name="vehicle",
        child_index=Int[],
        dependent_index=[1]
    )
else
    f = pnl.ReferenceFrame(rotor;
        origin=SVector{3}(0.0, 0.0, 0.0),
        v=SVector{3}(0.0, 0.0, 0.0),
        ω_axis=rotor_omega_axis,
        ω=0.0,
        R=SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name="vehicle",
        child_index=Int[],
        dependent_index=Int[]
    )
    for i in 1:nrotors
        pnl.add_frame!(f, "rotor$(i)", 1,
            SVector{3}(rotor_centers[i][1], rotor_centers[i][2], rotor_centers[i][3]),
            [i];
            ω_axis=rotor_omega_axis,
            ω=rotor_directions[i] * ω_full * spinup_fraction(t_range[1]))
    end
    f
end
rotor_frame_index(i) = nrotors == 1 ? 1 : 1 + i

# With DAS_REFRESH the pre-init is skipped: simulate! re-derives Das at every
# step (including the first) from the same eta, so pre-initializing here would
# double-accumulate through simulate!'s own initialize_Das! call.
!isnan(das_chord_fraction) && set_Das_refresh && error(
    "DAS_CHORD_FRACTION is incompatible with DAS_REFRESH (chord-proportional " *
    "Das is frozen by construction)")
!isnan(das_uniform_dsigma) && set_Das_refresh && error(
    "DAS_UNIFORM_DSIGMA is incompatible with DAS_REFRESH (uniform-clearance " *
    "Das is frozen by construction)")
!isnan(das_uniform_dsigma) && !isnan(das_chord_fraction) && error(
    "DAS_UNIFORM_DSIGMA and DAS_CHORD_FRACTION are mutually exclusive Das laws")
!isnan(das_sigma_lambda) && isnan(sigma_chord_fraction) && error(
    "DAS_SIGMA_LAMBDA requires SIGMA_CHORD_FRACTION (|Das|_j = lambda * sigma_j " *
    "needs the per-station sigma law)")
!isnan(das_sigma_lambda) && set_Das_refresh && error(
    "DAS_SIGMA_LAMBDA is incompatible with DAS_REFRESH (co-scaled Das is " *
    "frozen by construction)")
!isnan(das_sigma_lambda) && (!isnan(das_chord_fraction) || !isnan(das_uniform_dsigma)) && error(
    "DAS_SIGMA_LAMBDA, DAS_CHORD_FRACTION, and DAS_UNIFORM_DSIGMA are " *
    "mutually exclusive Das laws")
!isnan(das_curvature_beta) && isnan(das_sigma_lambda) && error(
    "DAS_CURVATURE_BETA requires DAS_SIGMA_LAMBDA (the F1 cap applies to " *
    "the co-scaled law)")
!isnan(das_curvature_beta) && das_curvature_beta <= 0 && error(
    "DAS_CURVATURE_BETA must be > 0 (radians of local helix arc)")
das_arc_placed && isnan(das_sigma_lambda) && error(
    "DAS_ARC_PLACED requires DAS_SIGMA_LAMBDA (the endpoint-on-arc placement " *
    "applies to the co-scaled law)")
das_arc_placed && !isnan(das_curvature_beta) && error(
    "DAS_ARC_PLACED and DAS_CURVATURE_BETA are mutually exclusive Phase-16 " *
    "fallbacks (arc placement vs F1 cap)")
das_arc_placed && !(das_arc_source in ("steady", "kinematic")) && error(
    "DAS_ARC_HELIX_SOURCE must be steady or kinematic (live is deferred: " *
    "Kutta attachment operator is assembled once per run)")
das_arc_placed && das_arc_source == "steady" && !isfile(das_arc_table) && error(
    "DAS_ARC_HELIX_SOURCE=steady requires DAS_ARC_TABLE=<csv> (per-station " *
    "induced-drift table from the TE downwash probe); got " *
    "$(repr(das_arc_table))")
theta_max_cs = NaN   # curvature diagnostic; finite only under the lambda law
if !set_Das_refresh
    if !isnan(das_sigma_lambda)
        # Phase 16: |Das|_j = lambda * sigma_j = lambda * s* * c_local_j, so
        # Das/c AND Das/sigma are span-uniform at once. Curvature diagnostic:
        # theta_j = N * |Das|_j / r_j is the arc the rigid extent subtends on
        # the local helix circle — under chord-scaling the binding station is
        # INBOARD (large chord, small r), the reverse of the eta law.
        das_station_lengths = (Tuple(
            begin
                lens_k = das_sigma_lambda .* station_sigmas[k]
                if !isnan(das_curvature_beta)
                    # F1 (2026-08-17 ladder verdict): cap the rigid extent at
                    # beta radians of local helix arc, theta_j <= beta.
                    r_k = station_radii(rotor.nodes, rotor.shedding[k],
                                        rotor.cells, radial_dimension)
                    lens_k = min.(lens_k, das_curvature_beta .* r_k ./ nwakerows)
                end
                lens_k
            end
            for k in eachindex(rotor.shedding)),)
        println("Das co-scaling mode: |Das|_j = $(das_sigma_lambda) * sigma_j " *
            "(= $(round(das_sigma_lambda * sigma_chord_fraction, digits=3)) * " *
            "c_local; DAS_ETA_KINEMATIC ignored)")
        if !isnan(das_curvature_beta)
            ncap = count(das_station_lengths[1][1] .<
                         das_sigma_lambda .* station_sigmas[1] .- eps())
            println("  F1 curvature cap ACTIVE: theta_j <= beta = " *
                "$(das_curvature_beta) rad (blade1: cap binds at $(ncap) of " *
                "$(length(station_sigmas[1])) stations)")
        end
        let lens = das_station_lengths[1][1],
            r = station_radii(rotor.nodes, rotor.shedding[1], rotor.cells, radial_dimension),
            sig = station_sigmas[1],
            c = station_chords(rotor.nodes, rotor.shedding[1], rotor.cells,
                               radial_dimension; band=0.02 * R)
            theta = nwakerows .* lens ./ r
            println("  blade1 station table (r/R, c/R, sigma/R, |Das|/R, theta rad):")
            for j in eachindex(lens)
                println("    $(round(r[j]/R, digits=4)) $(round(c[j]/R, digits=4)) " *
                    "$(round(sig[j]/R, digits=4)) $(round(lens[j]/R, digits=4)) " *
                    "$(round(theta[j], digits=4))")
            end
        end
        theta_max_cs = maximum(
            maximum(nwakerows .* das_station_lengths[1][k] ./
                    station_radii(rotor.nodes, rotor.shedding[k], rotor.cells,
                                  radial_dimension))
            for k in eachindex(rotor.shedding))
        println("  theta_max (all blades) = $(round(theta_max_cs, digits=4)) rad")
        das_station_drifts = nothing
        if das_arc_placed
            # F1b Route B: per-station drift vectors in the GLOBAL frame,
            # reconstructed from the table's local-TE-basis components using
            # the current geometry (basis co-rotates with the body by
            # construction). Basis: e_t = kinematic TE tangent -ω̂×d̂,
            # e_r = radial (spanwise) unit, e_n = e_t × e_r.
            axis_hat = frames[1].ω_axis / norm(frames[1].ω_axis)
            frame_origin = frames[1].x
            table_r = Float64[]; table_ut = Float64[]
            table_ur = Float64[]; table_un = Float64[]
            if das_arc_source == "steady"
                for line in eachline(das_arc_table)
                    sline = strip(line)
                    (isempty(sline) || startswith(sline, '#') ||
                        occursin("r_over_R", sline)) && continue
                    parts = split(sline, ',')
                    length(parts) >= 4 || continue
                    push!(table_r, parse(Float64, parts[1]))
                    push!(table_ut, parse(Float64, parts[2]))
                    push!(table_ur, parse(Float64, parts[3]))
                    push!(table_un, parse(Float64, parts[4]))
                end
                length(table_r) >= 2 || error("DAS_ARC_TABLE $(das_arc_table) " *
                    "has $(length(table_r)) usable rows; need >= 2")
                if !issorted(table_r)
                    perm = sortperm(table_r)
                    table_r = table_r[perm]; table_ut = table_ut[perm]
                    table_ur = table_ur[perm]; table_un = table_un[perm]
                end
                issorted(table_r) || error("DAS_ARC_TABLE r_over_R not sortable")
            end
            _interp(xs, ys, x) = begin   # linear, flat extrapolation
                x <= xs[1] && return ys[1]
                x >= xs[end] && return ys[end]
                i = searchsortedlast(xs, x)
                t = (x - xs[i]) / (xs[i+1] - xs[i])
                (1 - t) * ys[i] + t * ys[i+1]
            end
            das_station_drifts = (Tuple(
                begin
                    shed = rotor.shedding[k]
                    nshed = size(shed, 2)
                    drift = zeros(3, nshed + 1)
                    for j in 1:(nshed + 1)
                        node_idx = j <= nshed ?
                            rotor.cells[shed[3, j], shed[1, j]] :
                            rotor.cells[shed[2, nshed], shed[1, nshed]]
                        p = SVector{3}(rotor.nodes[1, node_idx],
                            rotor.nodes[2, node_idx], rotor.nodes[3, node_idx])
                        d = p - frame_origin
                        d_perp = d - axis_hat * dot(axis_hat, d)
                        rj = norm(d_perp)
                        rj > 0 || continue        # on-axis station: no basis
                        e_r = d_perp / rj
                        vte = -cross(axis_hat, d_perp)   # ∝ kinematic TE tangent
                        nv = norm(vte)
                        nv > 0 || continue
                        e_t = vte / nv
                        e_n = cross(e_t, e_r)
                        if das_arc_source == "steady"
                            roR = rj / R
                            u = _interp(table_r, table_ut, roR) .* e_t .+
                                _interp(table_r, table_ur, roR) .* e_r .+
                                _interp(table_r, table_un, roR) .* e_n
                            drift[:, j] .= u
                        end # kinematic: leave zero
                    end
                    drift
                end
                for k in eachindex(rotor.shedding)),)
            umax = maximum(maximum(sqrt.(sum(abs2, d; dims=1))) for d in das_station_drifts[1])
            println("  F1b endpoint-on-arc ACTIVE: source=$(das_arc_source)" *
                (das_arc_source == "steady" ? " table=$(basename(das_arc_table))" : "") *
                ", max |u| = $(round(umax, sigdigits=4)) m/s " *
                "($(round(umax / (ω_full * R), sigdigits=3)) of tip speed)")
        end
        pnl.initialize_Das!(Tuple(rotors), frames, Uinf, t_range[1], t_range[2] - t_range[1];
            set_Das_station_lengths=das_station_lengths,
            set_Das_station_drifts=das_station_drifts,
            set_Das_min_kinematic_displacement)
    elseif !isnan(das_uniform_dsigma)
        # phase_13 section 4b: |Das|_j = D*sigma - (N-1)*travel_j so that
        # d_front = D*sigma uniformly. travel at FULL rpm (the frozen-Das
        # convention would otherwise re-introduce the 0.4x spin-up factor
        # through the velocity magnitude; lengths here are absolute).
        dt_das = t_range[2] - t_range[1]
        # Station radii come from each rotor's LOCAL (untranslated) nodes: on a
        # translated body abs(node[radial_dimension]) is not the blade radius.
        lens_uniform_by_body = [map(rb.body.shedding) do shed
            r = station_radii(rb.local_nodes, shed, rb.body.cells, radial_dimension)
            das_uniform_dsigma * tip_sigma_default .- (nwakerows - 1) .* ω_full .* r .* dt_das
        end for rb in rotor_builds]
        lens_uniform = lens_uniform_by_body[1]
        for (ib, lens_body) in enumerate(lens_uniform_by_body), (k, lens) in enumerate(lens_body)
            minimum(lens) > 0 || error("DAS_UNIFORM_DSIGMA=$(das_uniform_dsigma) " *
                "infeasible: rotor$(ib) shedding$(k) min |Das| = $(minimum(lens)) <= 0 at " *
                "N=$(nwakerows) (tip travel exceeds D*sigma; raise D or lower N)")
        end
        band = 0.02 * R
        chords1 = station_chords(rotor_builds[1].local_nodes, rotor.shedding[1], rotor.cells, radial_dimension; band)
        dasc = lens_uniform[1] ./ chords1
        println("Das uniform-clearance mode: |Das|_j = $(das_uniform_dsigma)*sigma" *
            " - (N-1)*travel_j (DAS_ETA_KINEMATIC ignored), sigma = " *
            "$(round(tip_sigma_default, sigdigits=5)) m, N = $(nwakerows)")
        println("  blade1 |Das|/R range $(round(minimum(lens_uniform[1])/R, digits=4))" *
            "-$(round(maximum(lens_uniform[1])/R, digits=4)); Das/c_local range " *
            "$(round(minimum(dasc), digits=3))-$(round(maximum(dasc), digits=3)) " *
            "(014 band 0.25-1.5)")
        pnl.initialize_Das!(Tuple(rotors), frames, Uinf, t_range[1], t_range[2] - t_range[1];
            set_Das_station_lengths=Tuple(Tuple(lens_body) for lens_body in lens_uniform_by_body),
            set_Das_min_kinematic_displacement)
    elseif !isnan(das_chord_fraction)
        band = 0.02 * R  # ~ one spanwise ring spacing on the 45-ring blade
        chords1 = station_chords(rotor.nodes, rotor.shedding[1], rotor.cells, radial_dimension; band)
        chords2 = station_chords(rotor.nodes, rotor.shedding[2], rotor.cells, radial_dimension; band)
        das_station_lengths = ((das_chord_fraction .* chords1,
                                das_chord_fraction .* chords2),)
        println("Das chord mode: |Das| = $(das_chord_fraction) * local chord " *
            "(DAS_ETA_KINEMATIC ignored)")
        println("  blade1 chords/R: min $(round(minimum(chords1)/R, digits=4)) " *
            "max $(round(maximum(chords1)/R, digits=4)); " *
            "|Das|/R range $(round(das_chord_fraction*minimum(chords1)/R, digits=4))" *
            "-$(round(das_chord_fraction*maximum(chords1)/R, digits=4))")
        pnl.initialize_Das!(Tuple(rotors), frames, Uinf, t_range[1], t_range[2] - t_range[1];
            set_Das_station_lengths=das_station_lengths,
            set_Das_min_kinematic_displacement)
    else
        pnl.initialize_Das!(Tuple(rotors), frames, Uinf, t_range[1], t_range[2] - t_range[1];
            set_Das_eta_kinematic=init_Das_eta_kinematic,
            set_Das_min_kinematic_displacement,
            set_Das_kinematic_arc)
    end
end

# --- Panel->particle handoff clearance profile (BRAINSTORM/018 S0c) ---------
# The distance from the TE at which the panel wake hands off to particles is
# NOT a knob: it is emergent from (Das, nwakerows, dt). Row 1 is rigidly pinned
# at TE+Das; rows 2..N convect freely and are separated by one step of TE
# travel (the Das cancels between consecutive rows); particles are deposited on
# the segment spanning node rows N -> N+1 (`_convert_to_particles!`,
# src/FLOWPanel_wake.jl). Hence, per station j,
#
#     d_front,j = |Das_j| + (N-1) * |w x r_j| * dt
#
# and the quantity that governs whether the blade sees a particle singularity
# is d_front/sigma -- BRAINSTORM/018 Phase 12 C1 measured the admissible value
# as d/sigma* ~ 2.6 (median) / 3.4 (p95). This block reports the profile so the
# clearance is auditable from the run's own banner instead of being inferred.
#
# NOTE the span dependence: |w x r| ~ r while sigma is span-uniform under
# SigmaOverlap, so d/sigma rises outboard and the BINDING station is the
# innermost one. Judge the pin on the minimum, not the tip.
#
# Placed before `Backslash` deliberately: that constructor assembles and
# LU-factors a dense ncells^2 matrix (~10 GB at 45_185_ct4), so gating an early
# exit here makes the profile a seconds-long query.
function das_handoff_profile(nodes, shedding, cells, Das, radial_dimension,
                             omega_full, dt, nrows, sigma)
    nshed = size(shedding, 2)
    station_nodes = Vector{Int}(undef, nshed + 1)
    for j in 1:nshed
        station_nodes[j] = cells[shedding[3, j], shedding[1, j]]   # nib
    end
    station_nodes[nshed + 1] = cells[shedding[2, nshed], shedding[1, nshed]]  # final nia
    r = [abs(nodes[radial_dimension, n]) for n in station_nodes]
    das = [sqrt(Das[1, j]^2 + Das[2, j]^2 + Das[3, j]^2) for j in axes(Das, 2)]
    travel = omega_full .* r .* dt          # one step of TE travel at FULL rpm
    d_front = das .+ (nrows - 1) .* travel
    return (; station_nodes, r, das, travel, d_front, dsigma=d_front ./ sigma)
end

let sigma_handoff = tip_sigma_default, dt_handoff = t_range[2] - t_range[1]
    println("Handoff clearance (d_front = |Das| + (N-1)*|w x r|*dt, N=$(nwakerows)):")
    println("  sigma = $(round(sigma_handoff, sigdigits=5)) m  " *
            "(sigma/R = $(round(sigma_handoff / R, digits=5)), " *
            "one-step tip travel = $(round(ω_full * R * dt_handoff, sigdigits=5)) m)")
    csv_path = get(ENV, "RHPC_HANDOFF_CSV", "")
    rows = String[]
    for (ib, rb) in enumerate(rotor_builds), (k, shed) in enumerate(rb.body.shedding)
        # local_nodes, not the translated body's: station radii are blade radii
        p = das_handoff_profile(rb.local_nodes, shed, rb.body.cells, rb.body.Das[k],
                                radial_dimension, ω_full, dt_handoff, nwakerows,
                                sigma_handoff)
        jmin = argmin(p.dsigma)
        println("  $(nrotors == 1 ? "" : "rotor$(ib) ")shedding$(k): d/sigma min $(round(p.dsigma[jmin], digits=3)) " *
                "at r/R $(round(p.r[jmin] / R, digits=3)); " *
                "max $(round(maximum(p.dsigma), digits=3)) " *
                "at r/R $(round(p.r[argmax(p.dsigma)] / R, digits=3)); " *
                "median $(round(sort(p.dsigma)[cld(length(p.dsigma), 2)], digits=3))")
        for j in eachindex(p.dsigma)
            push!(rows, "$((ib - 1) * length(rb.body.shedding) + k),$(j)," *
                        "$(p.r[j] / R),$(p.das[j]),$(p.travel[j])," *
                        "$(p.d_front[j]),$(p.dsigma[j])")
        end
    end
    if !isempty(csv_path)
        mkpath(dirname(csv_path))
        open(csv_path, "w") do io
            println(io, "i_shed,station,r_over_R,das,step_travel,d_front,d_over_sigma")
            for row in rows; println(io, row); end
        end
        println("  wrote $(csv_path)")
    end
end
parse(Bool, get(ENV, "RHPC_HANDOFF_PROFILE_ONLY", "false")) && exit(0)

solvers_rotor = [pnl.Backslash(r) for r in rotors]
solver_rotor = solvers_rotor[1]   # legacy alias
# RHPC_BACKEND=direct is the BRAINSTORM/018 Phase-2 discriminator for the FMM
# |Das|-panel-radius coupling (src/FLOWPanel_liftingbody.jl); fmm is production.
rhpc_backend = get(ENV, "RHPC_BACKEND", "fmm")
# BRAINSTORM/023 + 025 (ported from rotor_hover_pressure_comparison.jl): the
# body and wake passes tune independently, so each takes its own knob triple.
# The shared FMM_* names remain the fallback, so any existing submission that
# sets only FMM_EXPANSION_ORDER / FMM_ACCEPTANCE / FMM_LEAF_SIZE behaves
# exactly as before.
rhpc_fmm_knob(specific, shared, default) = get(ENV, specific, get(ENV, shared, default))
backend, backend_wake = if rhpc_backend == "direct"
    pnl.DirectBackend(), pnl.DirectBackend()
elseif rhpc_backend == "fmm"
    pnl.FastMultipoleBackend(;
        expansion_order=parse(Int, rhpc_fmm_knob("FMM_BODY_EXPANSION_ORDER", "FMM_EXPANSION_ORDER", "8")),
        multipole_acceptance=parse(Float64, rhpc_fmm_knob("FMM_BODY_ACCEPTANCE", "FMM_ACCEPTANCE", "0.4")),
        leaf_size=parse(Int, rhpc_fmm_knob("FMM_BODY_LEAF_SIZE", "FMM_LEAF_SIZE", "20")),
    ),
    pnl.FastMultipoleBackend(;
        expansion_order=parse(Int, rhpc_fmm_knob("FMM_WAKE_EXPANSION_ORDER", "FMM_EXPANSION_ORDER", "4")),
        multipole_acceptance=parse(Float64, rhpc_fmm_knob("FMM_WAKE_ACCEPTANCE", "FMM_ACCEPTANCE", "0.4")),
        leaf_size=parse(Int, rhpc_fmm_knob("FMM_WAKE_LEAF_SIZE", "FMM_LEAF_SIZE", "50")),
    )
else
    error("Unknown RHPC_BACKEND=$(repr(rhpc_backend)); use fmm or direct")
end
kj_backend = pnl.FastMultipoleBackend(;
    expansion_order=parse(Int, get(ENV, "KJ_FMM_EXPANSION_ORDER", "3")),
    multipole_acceptance=parse(Float64, get(ENV, "KJ_FMM_ACCEPTANCE", "0.4")),
    leaf_size=parse(Int, get(ENV, "KJ_FMM_LEAF_SIZE", "1000")),
)
formulation = if formulation_name == "velocity"
    pnl.VelocityThroughSources()
elseif formulation_name == "green"
    pnl.GreenReconstruction(; gauge=green_gauge,
        recompute_interval=green_recompute_interval, green_solver=nothing)
else
    error("Unknown RHPC_FORMULATION=$(repr(formulation_name)); use velocity or green")
end

function maneuver!(frames, systems, wakes, t)
    for i in 1:nrotors
        set_frame_omega!(frames, rotor_frame_index(i),
            rotor_directions[i] * ω_full * spinup_fraction(t))
    end
    return nothing
end

# BRAINSTORM/022: with the ground on, simulate! receives per-body solvers and
# the block Gauss-Seidel outer loop in solve!(bodies::Tuple, solvers::Tuple)
# iterates the rotor <-> ground coupling to convergence. The ground is a
# NonLiftingBody with wake `nothing`, and is deliberately in NO frame's
# dependent_index: propagate_kinematics!/kinematic_velocity! only touch frame
# dependents, so the ground stays static with freestream-only onset flow.
if ground_enable
    systems = (rotors..., ground)
    wakes = (wakes_rotor..., nothing)
    body_solvers = (solvers_rotor..., solver_ground)
else
    systems = Tuple(rotors)
    wakes = Tuple(wakes_rotor)
    body_solvers = Tuple(solvers_rotor)
end

# --- GS outer-loop instrumentation (022 Phase 0 requirement) ------------------
# solve_formulation!(::VelocityThroughSources) is a one-liner that forwards to
# solve! WITHOUT the verbose/history kwargs, and solve!'s non-convergence is
# silent unless verbose (src/FLOWPanel_solver.jl). Overwrite that one method
# here (driver-local; src untouched) to pass a driver-owned ConvergenceHistory
# into the EXISTING outer loop and log per-step iteration counts. No new
# iteration logic. Julia prints a method-overwrite warning; expected.
gs_iters = Int[]
gs_final_delta = Float64[]
if (ground_enable || nrotors > 1) && gs_log
    formulation_name == "velocity" || error(
        "GS logging (and the 022 ground path generally) assumes " *
        "RHPC_FORMULATION=velocity; got $(repr(formulation_name)). " *
        "GreenReconstruction with a ground body is untested — do not combine.")
    gs_history = pnl.ConvergenceHistory(:blockgs_maxdelta)
    function pnl.solve_formulation!(::pnl.VelocityThroughSources, state, systems,
            systems_tuple, wakes_tuple, body_solvers;
            backend_solve, backend_wake, i_step::Int=0)
        pnl.solve!(systems, body_solvers; backend=backend_solve,
            max_outer_iterations=gs_max_outer, outer_tolerance=gs_tol,
            verbose=gs_verbose, history=gs_history)
        n = length(gs_history.iter)
        final = n == 0 ? NaN : gs_history.residual_internal[end]
        push!(gs_iters, n)
        push!(gs_final_delta, final)
        if n >= gs_max_outer && !(final < gs_tol)
            @warn "Block Gauss-Seidel outer loop hit max_outer_iterations " *
                "without converging" i_step n final gs_tol
        end
        return nothing
    end
end
# -----------------------------------------------------------------------------

# --- Near-ground vertical velocity cutoff (022 Phase 3, Ryan 2026-08-20) -----
# Driver-local overwrite of propagate!(::PanelParticleWake, ...): body copied
# verbatim from src/FLOWPanel_wake.jl @ b251071 with ONE insertion — just
# before the Euler position update, scale the ground-ward axial velocity of
# every particle within the damping band by the linear factor f = d/band
# (d = height above the ground; f = 0 at/below the plane). Receding motion is
# untouched. _euler convects with the STORED particle U (pfield.Uinf is
# nofreestream here; apply_freestream! already folded the freestream into U),
# so this caps exactly the velocity that moves particles. With dt = 2.8e-4 s
# and band = R/10, a particle cannot cross the plane below ~43 m/s axial.
# Counters feed the ground diagnostics monitor for the pass-through A/B.
ground_damp_last_n = Ref(0)        # particles damped in the latest step
ground_damp_last_inband = Ref(0)   # particles inside the band (above ground)
ground_damp_total = Ref(0)         # cumulative damped-particle-steps
if ground_enable && ground_damp_band_r > 0
    damp_band = ground_damp_band_r * R
    function pnl.propagate!(w::pnl.PanelParticleWake, dt; relax=true, step=0,
            frames=nothing, diagnose_particle_gamma::Bool=false,
            diagnostic_vertical=(0.0, 0.0, 1.0))

        # panel wake
        pnl.propagate!(w.panel_wake, dt)

        gamma_before = diagnose_particle_gamma && w.pfield.np > 0 ?
            copy(view(w.pfield.particles, pnl.FLOWVPM.GAMMA_INDEX, 1:w.pfield.np)) : nothing
        if diagnose_particle_gamma
            println("particle gamma step=$(step) phase=before_euler " *
                pnl._particle_gamma_direction_stats(w.pfield; vertical=diagnostic_vertical))
        end

        # refresh any frame-tracking relaxation filter to the current frame pose
        pnl.refresh_relaxation_filter!(w.pfield.relaxation.filter, frames)

        # >>> 022 insertion: linear vertical velocity cutoff near the ground.
        n_damped = 0
        n_inband = 0
        u_ax_index = pnl.FLOWVPM.U_INDEX[axial_dimension]
        x_ax_index = pnl.FLOWVPM.X_INDEX[axial_dimension]
        for i in 1:w.pfield.np
            d = ground_x - w.pfield.particles[x_ax_index, i]  # height above ground
            d < damp_band || continue
            0 <= d && (n_inband += 1)
            u_ax = w.pfield.particles[u_ax_index, i]
            if u_ax > 0                       # moving toward the ground (+axial)
                f = clamp(d / damp_band, 0.0, 1.0)
                w.pfield.particles[u_ax_index, i] = f * u_ax
                n_damped += 1
            end
        end
        ground_damp_last_n[] = n_damped
        ground_damp_last_inband[] = n_inband
        ground_damp_total[] += n_damped
        # <<< end 022 insertion

        # convect particles (`pfield.integration` is the single source of truth:
        # stock forward Euler unless the wake was built with `expint=true`)
        if w.pfield.integration === pnl.FLOWVPM.euler_exp
            pnl.FLOWVPM._euler_exp(w.pfield, dt; relax)
        else
            pnl.FLOWVPM._euler(w.pfield, dt; relax)
        end

        if diagnose_particle_gamma
            println("particle gamma step=$(step) phase=after_euler relax=$(relax) " *
                pnl._particle_gamma_direction_stats(w.pfield;
                    vertical=diagnostic_vertical, before_gamma=gamma_before))
        end

        # particle maintenance
        pnl.apply_particle_maintenance!(w.pfield, w.particle_maintenance,
            pnl.ParticleMaintenanceContext(frames, step, dt))
    end
end
# -----------------------------------------------------------------------------

pressure_bernoulli = pnl.PressureBernoulli(rho;
    unsteady=false,
    correct_kuttacondition=p_correct_kuttacondition_flag,
    backend=backend)
# One ForceMonitor per rotor (i_system = i), each normalized by its OWN frame's
# omega (RotorNormalization uses n^2, so counter-rotation is sign-safe). At
# NROTORS=1 this is exactly the legacy single monitor.
force_monitors_bernoulli = [pnl.ForceMonitor(length(t_range), i;
    i_frame=rotor_frame_index(i),
    normalization=pnl.RotorNormalization(rho, 2 * R, rotor_frame_index(i)),
    correct_kuttacondition=p_correct_kuttacondition_flag,
    verbose=true) for i in 1:nrotors]
force_monitor_bernoulli = force_monitors_bernoulli[1]   # legacy alias

if nrotors == 1
pressure_laplace_matderiv = pnl.PressureLaplace(rotor, rho;
    acceleration_form=:material_derivative, verbose=true)
force_monitor_laplace_matderiv = pnl.ForceMonitor(length(t_range), 1;
    i_frame=1,
    normalization=pnl.RotorNormalization(rho, 2 * R, 1),
    correct_kuttacondition=p_correct_kuttacondition_flag,
    verbose=true)

# DIAGNOSTIC ONLY: deprecated Lamb mode is retained for this historical comparison.
pressure_laplace_lamb = pnl.PressureLaplace(rotor, rho;
    acceleration_form=:lamb_vector, verbose=true)
force_monitor_laplace_lamb = pnl.ForceMonitor(length(t_range), 1;
    i_frame=1,
    normalization=pnl.RotorNormalization(rho, 2 * R, 1),
    correct_kuttacondition=p_correct_kuttacondition_flag,
    verbose=true)

kj_monitor = run_kj ? pnl.KuttaJoukowskiForce(rotor, length(t_range), 1;
        rho, backend=kj_backend,
        i_frame=1,
        normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        verbose=true) : nothing
bound_circulation = pnl.BoundCirculationMonitor(rotor, length(t_range), 1;
    i_frame=1,
    radial_dimension,
    R)

monitors = !run_monitors ? () : run_kj ? (
        pressure_laplace_matderiv,
        force_monitor_laplace_matderiv,
        pressure_laplace_lamb,
        force_monitor_laplace_lamb,
        pressure_bernoulli,
        force_monitor_bernoulli,
        kj_monitor,
        bound_circulation,
    ) : lamb_only ? (
        pressure_laplace_lamb,
        force_monitor_laplace_lamb,
        bound_circulation,
    ) : bernoulli_only ? (
        pressure_bernoulli,
        force_monitor_bernoulli,
        bound_circulation,
    ) : (
        pressure_laplace_matderiv,
        force_monitor_laplace_matderiv,
        pressure_laplace_lamb,
        force_monitor_laplace_lamb,
        pressure_bernoulli,
        force_monitor_bernoulli,
        bound_circulation,
    )
else
# NROTORS>1: bernoulli_only is enforced above; the Laplace/KJ family and the
# BoundCirculationMonitor (whose radius bookkeeping assumes an on-axis rotor)
# are not constructed. Order matters: pressure before the force monitors.
monitors = !run_monitors ? () : (pressure_bernoulli, force_monitors_bernoulli...)
end

# Wake tripwire (BRAINSTORM/018 S0a). Appended LAST so existing monitor indices
# -- and therefore existing monitor CSV filenames -- are unchanged. Diagnostic
# only: provides/requires nothing and mutates no simulation state, so runs stay
# bit-identical; it adds one CSV. Default ON deliberately: a divergence tripwire
# you have to remember to enable is not a tripwire.
wake_health_active = parse(Bool, get(ENV, "WAKE_HEALTH", "true"))
# Default-off extra column: max per-step core contraction max_p dt*Z_p (see
# WakeHealthMonitor docstring). Off keeps the CSV bit-identical to legacy.
wake_health_dtz = parse(Bool, get(ENV, "WAKE_HEALTH_DTZ", "false"))
# min_sr attribution columns (BRAINSTORM/018 phase 15): p1 percentile of
# sigma/sigma_shed + argmin position, appended AFTER the legacy columns so
# positional readers of the existing columns are unaffected. Default ON so
# campaign runs carry the attribution for free.
wake_health_attribution = parse(Bool, get(ENV, "WAKE_HEALTH_ATTRIBUTION", "true"))
if run_monitors && wake_health_active
    monitors = (monitors..., pnl.WakeHealthMonitor(dtz=wake_health_dtz,
                                 attribution=wake_health_attribution))
end

# Banded wake inventory (BRAINSTORM/018 phase 15 drift-source study): per-step
# per-band count, sum|Gamma|, vector-sum Gamma, and sigma quantiles. Appended
# LAST so existing monitor CSV names are unchanged; diagnostic only. Bands end
# at the MEASURED deletion boundary 3.5R downstream (the GlobalCylinder above
# starts 0.5R upstream, so its 4R length = 3.5R downstream extent). Default ON
# so the dt-ladder rungs double as the H1/H3 measurement runs.
wake_inventory_active = parse(Bool, get(ENV, "WAKE_INVENTORY", "true"))
if run_monitors && wake_inventory_active
    monitors = (monitors..., pnl.WakeInventoryMonitor(R))
end

# --- BRAINSTORM/022 ground diagnostics monitor --------------------------------
# (a) Flow tangency on the ground: RMS and max of U·n over ground control
#     points (flat_ground.jl pattern) — corroborates that the coupled solve
#     actually enforced no-flow-through each step.
# (b) Below-ground particle census: count and sum|Gamma| of particles past the
#     plane — measures the cost of the leave-be policy before deciding on cull.
# Appended LAST (after WakeHealth/WakeInventory) so existing monitor CSV
# indices are unchanged. Diagnostic only; provides/requires nothing.
ground_steps = Int[]
ground_tangency_rms = Float64[]
ground_tangency_max = Float64[]
ground_below_count = Int[]
ground_below_gamma = Float64[]
function ground_diagnostics_monitor(systems, wakes, frames, uinf, i_step, dt)
    g = systems[nrotors + 1]
    normals = pnl.calc_normals(g)
    Udotn = sum(g.velocity .* normals, dims=1)
    rms = sqrt(sum(abs2, Udotn) / g.ncells)
    maxerr = maximum(abs.(Udotn))
    nbelow = 0
    gbelow = 0.0
    for w in wakes[1:nrotors]
        pfield = w.pfield
        for i in 1:pfield.np
            x = pnl.FLOWVPM.get_X(pfield, i)
            if x[axial_dimension] > ground_x
                nbelow += 1
                G = pnl.FLOWVPM.get_Gamma(pfield, i)
                gbelow += sqrt(G[1]^2 + G[2]^2 + G[3]^2)
            end
        end
    end
    push!(ground_steps, i_step)
    push!(ground_tangency_rms, rms)
    push!(ground_tangency_max, maxerr)
    push!(ground_below_count, nbelow)
    push!(ground_below_gamma, gbelow)
    if i_step % 10 == 0
        println("    Ground step $i_step: RMS(U·n)=$(round(rms; sigdigits=4)) " *
            "max|U·n|=$(round(maxerr; sigdigits=4)) " *
            "below-ground particles=$nbelow sum|Γ|=$(round(gbelow; sigdigits=4)) " *
            "damped=$(ground_damp_last_n[]) in-band=$(ground_damp_last_inband[])")
    end
    # Incremental CSVs (2026-08-20 lesson: the post-run dump loses everything
    # on a walled/cancelled/diverged run — 13207682's census survived only in
    # the slurm log). Appended every step; header written at first call.
    if save_path !== nothing
        isdir(save_path) || mkpath(save_path)
        gd_csv = joinpath(save_path, "$(run_name)_ground_diagnostics.csv")
        write_header = !isfile(gd_csv) || length(ground_steps) == 1
        open(gd_csv, write_header ? "w" : "a") do io
            write_header && println(io, "i_step,tangency_rms,tangency_max," *
                "below_ground_count,below_ground_sum_gamma,n_damped,n_inband")
            println(io, "$i_step,$rms,$maxerr,$nbelow,$gbelow," *
                "$(ground_damp_last_n[]),$(ground_damp_last_inband[])")
        end
        if gs_log
            gs_csv = joinpath(save_path, "$(run_name)_gs_convergence.csv")
            gs_header = !isfile(gs_csv) || gs_csv_written[] == 0
            open(gs_csv, gs_header ? "w" : "a") do io
                gs_header && println(io, "solve_index,outer_iterations,final_max_delta")
                for k in (gs_csv_written[] + 1):length(gs_iters)
                    println(io, "$k,$(gs_iters[k]),$(gs_final_delta[k])")
                end
            end
            gs_csv_written[] = length(gs_iters)
        end
    end
    return nothing
end
gs_csv_written = Ref(0)
if run_monitors && ground_enable
    monitors = (monitors..., ground_diagnostics_monitor)
end

if spinup_revs > 0
    println("\nRotor spin-up enabled: $(spinup_revs) nominal revs from $(spinup_start_fraction)×RPM to full RPM")
end
let
    r0 = freestream_ramp_revs
    r1 = r0 + freestream_hold_revs
    r2 = r1 + freestream_withdraw_revs
    println("\nFreestream pulse (revs): ramp-up [0, $(round(r0,digits=2))] -> peak=$(magVinf_peak); " *
        "hold [$(round(r0,digits=2)), $(round(r1,digits=2))]; " *
        "withdraw [$(round(r1,digits=2)), $(round(r2,digits=2))] -> end=$(magVinf_end); " *
        "hover/settle [$(round(r2,digits=2)), $(round(required_revs,digits=2))]")
    println("Total run length: $(round(required_revs,digits=2)) revs ($(length(t_range)) steps), truncation depth=$(round(cylinder_depth/R,digits=2))R")
end
println("\nBegin rotor hover ground effect run ($(length(t_range)) steps)...")
println("Mesh=$(rhpc_mesh) file=$(basename(msh_file)) formulation=$(formulation_name) " *
        "RPM=$(RPM) rho=$(rho) R=$(R) NT=$(nt) truncation_depth=$(round(cylinder_depth/R,digits=3))R " *
        "trunc_radius=$(trunc_radius_r)R nwakerows=$(nwakerows) das_refresh=$(set_Das_refresh)")
nrotors > 1 && println("Rotors: NROTORS=$(nrotors), spacing=$(rotor_spacing_r)R " *
    "(center-to-center), directions=$(rotor_directions), centers/R=" *
    "$([round.(c ./ R, digits=3) for c in rotor_centers]), " *
    "max_offset=$(round(max_rotor_offset / R, digits=3))R")
println("Ground: GROUND_ENABLE=$(ground_enable)" * (ground_enable ?
        ", h/R=$(ground_h_r), disc_radius=$(ground_radius_r)R, " *
        "panel_length=$(ground_panel_length_r)R, ncells=$(ground.ncells), " *
        "particle_policy=$(ground_particle_policy), " *
        "damp_band=$(ground_damp_band_r)R, " *
        "effective_depth=$(round(effective_cylinder_depth/R, digits=3))R, " *
        "GS_LOG=$(gs_log), GS_MAX_OUTER=$(gs_max_outer), GS_TOL=$(gs_tol)" :
        " (OGE baseline)"))
println("Particle diagnostics: PARTICLE_SHEDDING=$(particle_shedding), CONVERSION=$(conversion_mode)$(conversion_mode == "smooth" ? ", CONVERSION_SIGMA=$(conversion_sigma), CONVERSION_OVERLAP=$(conversion_overlap), ATTRIBUTION=$(conversion_attribution)" : ""), RUN_MONITORS=$(run_monitors), BODY_HESSIAN_TO_PARTICLES=$(body_hessian_to_particles), PANEL_WAKE_HESSIAN_TO_PARTICLES=$(panel_wake_hessian_to_particles), PANEL_WAKE_VELOCITY_TO_PARTICLES=$(panel_wake_on_particles), PARTICLE_HESSIAN_SELF=$(particle_hessian_self), PARTICLE_RELAX=$(particle_relax), DIAGNOSE_PARTICLE_GAMMA=$(diagnose_particle_gamma), DIAGNOSE_PARTICLE_INFLUENCE=$(diagnose_particle_influence), diagnostic_vertical=$(particle_diagnostic_vertical), WAKE_HEALTH=$(wake_health_active), WAKE_HEALTH_DTZ=$(wake_health_dtz), WAKE_HEALTH_ATTRIBUTION=$(wake_health_attribution), WAKE_INVENTORY=$(wake_inventory_active), WAKE_EXPINT=$(wake_expint)")
name = run_name

# Allow other scripts to `include` this file purely for setup (geometry,
# frames, wake, monitors) without executing the time-marching call.
rhpc_setup_only = parse(Bool, get(ENV, "RHPC_SETUP_ONLY", "false"))

if !rhpc_setup_only

# Warm-start mode: RESTART_STEP >= 0 resumes from a previous run's saved state
# (default: this example's own baseline output) and continues with the CURRENT
# ENV configuration — the basis for knob-perturbation experiments. The body,
# wake, and frames objects must be construction-compatible with the restart
# run (same RHPC_MESH, NT, RPM); construction-time knobs (RELAX_RLXF,
# CORE_SIZE_*, ...) and simulate!-kwarg knobs (BODY_HESSIAN_TO_PARTICLES,
# ...) both take effect in the continuation.
restart_step = parse(Int, get(ENV, "RESTART_STEP", "-1"))
restart_step >= 0 && nrotors > 1 && error(
    "Warm-start with NROTORS>1 is UNTESTED (multi-body restart state); run cold")
restart_step >= 0 && ground_enable && error(
    "Warm-start with the ground body is UNTESTED (multi-body restart state); " *
    "run cold or set GROUND_ENABLE=false")
restart_name = get(ENV, "RESTART_NAME", "rotor_hover_ground_effect")
restart_path = get(ENV, "RESTART_PATH", joinpath("data", restart_name))

sim_wall_start = time()

if restart_step >= 0

@time pnl.simulate_warmstart!(systems, wakes, frames, maneuver!, Uinf, t_range;
    restart_path, restart_name, restart_step,
    set_Das_eta_kinematic=set_Das_refresh ? init_Das_eta_kinematic : NaN,
    set_Das_min_kinematic_displacement,
    set_Das_kinematic_arc,
    set_Das_refresh,
    monitors,
    formulation,
    body_solvers, backend, backend_wake,
    wakerow_no_hessian_to_particles,
    body_hessian_to_particles,
    body_gradient_core_size,
    body_on_wake,
    panel_wake_on_particles,
    particle_hessian_self,
    particle_relax,
    bound_strength_rlx,
    diagnose_particle_gamma,
    diagnose_particle_influence,
    diagnostic_vertical=particle_diagnostic_vertical,
    verbose=true,
    path=save_path, name,
)

else

@time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
    set_Das_eta_kinematic=set_Das_refresh ? init_Das_eta_kinematic : NaN,
    set_Das_min_kinematic_displacement,
    set_Das_kinematic_arc,
    set_Das_refresh,
    monitors,
    formulation,
    body_solvers, backend, backend_wake,
    wakerow_no_hessian_to_particles,
    body_hessian_to_particles,
    body_gradient_core_size,
    body_on_wake,
    panel_wake_on_particles,
    particle_hessian_self,
    particle_relax,
    bound_strength_rlx,
    diagnose_particle_gamma,
    diagnose_particle_influence,
    diagnostic_vertical=particle_diagnostic_vertical,
    verbose=true,
    path=save_path, name,
)

end

sim_wall_seconds = time() - sim_wall_start
println("\nTime marching wall time: $(round(sim_wall_seconds, digits=1)) s")

# --- BRAINSTORM/022 post-run diagnostics --------------------------------------
gs_iters_max = isempty(gs_iters) ? 0 : maximum(gs_iters)
gs_iters_mean = isempty(gs_iters) ? NaN : sum(gs_iters) / length(gs_iters)
gs_nonconverged = count(i -> gs_iters[i] >= gs_max_outer &&
    !(gs_final_delta[i] < gs_tol), eachindex(gs_iters))
if (ground_enable || nrotors > 1) && gs_log
    println("\nBlock Gauss-Seidel outer-loop summary ($(length(gs_iters)) solves):")
    println("  iterations: max=$(gs_iters_max) mean=$(round(gs_iters_mean, digits=2)) " *
        "cap=$(gs_max_outer)")
    println("  non-converged solves: $(gs_nonconverged)")
    # CSVs are written incrementally by the ground diagnostics monitor
    # (walled-run lesson, 2026-08-20); top up any GS entries logged after the
    # final monitor call.
    if save_path !== nothing && gs_csv_written[] < length(gs_iters)
        gs_csv = joinpath(save_path, "$(run_name)_gs_convergence.csv")
        open(gs_csv, gs_csv_written[] == 0 ? "w" : "a") do io
            gs_csv_written[] == 0 &&
                println(io, "solve_index,outer_iterations,final_max_delta")
            for k in (gs_csv_written[] + 1):length(gs_iters)
                println(io, "$k,$(gs_iters[k]),$(gs_final_delta[k])")
            end
        end
        gs_csv_written[] = length(gs_iters)
    end
end
if ground_enable && ground_damp_band_r > 0
    println("Ground damping summary: band=$(ground_damp_band_r)R, cumulative " *
        "damped particle-steps=$(ground_damp_total[])")
end
# -----------------------------------------------------------------------------

# BRAINSTORM/016: external conservation evidence for the smooth conversion.
# Sibling file (same shape as the SSW driver's conversion_diagnostics.csv);
# legacy runs carry `conversion_diagnostics === nothing` and emit nothing, so
# their outputs stay bit-identical. On a restart-chained run this ledger covers
# only the final segment (in-memory accumulator).
let wake = wakes[1]
    if save_path !== nothing && wake isa pnl.PanelParticleWake &&
            wake.conversion_diagnostics !== nothing
        d = wake.conversion_diagnostics
        ledger_path = joinpath(save_path, "$(run_name)_conversion_diagnostics.csv")
        open(ledger_path, "w") do io
            println(io, "conversions,particles,expected_x,expected_y," *
                "expected_z,deposited_x,deposited_y,deposited_z,residual_norm")
            println(io, "$(wake.conversion_count[]),$(d.total_particles)," *
                "$(d.total_expected[1]),$(d.total_expected[2])," *
                "$(d.total_expected[3]),$(d.total_deposited[1])," *
                "$(d.total_deposited[2]),$(d.total_deposited[3])," *
                "$(norm(d.total_residual))")
        end
        println("Wrote conversion ledger: $ledger_path " *
            "(residual_norm = $(norm(d.total_residual)))")
    end
end

if !run_monitors
    println("\nRUN_MONITORS=false; skipping pressure/force history summary.")
else

if nrotors == 1
CT_bernoulli   = lamb_only ? fill(NaN, length(t_range)) : force_monitor_bernoulli.force[axial_dimension, :]
CT_laplace_md  = lamb_only ? fill(NaN, length(t_range)) : force_monitor_laplace_matderiv.force[axial_dimension, :]
CT_laplace_lv  = force_monitor_laplace_lamb.force[axial_dimension, :]
CT_kj          = run_kj ? kj_monitor.force[axial_dimension, :] : fill(NaN, length(t_range))
CT_bernoulli_by_rotor = [CT_bernoulli]
else
# Per-rotor CT plus the layout total. Every rotor shares rho, D, and |omega|,
# so the per-rotor coefficients share one reference and the total is their sum.
CT_bernoulli_by_rotor = [m.force[axial_dimension, :] for m in force_monitors_bernoulli]
CT_bernoulli   = sum(CT_bernoulli_by_rotor)
CT_laplace_md  = fill(NaN, length(t_range))
CT_laplace_lv  = fill(NaN, length(t_range))
CT_kj          = fill(NaN, length(t_range))
end

function relative_difference(a, b)
    denom = max(abs(b), eps())
    return abs(a - b) / denom
end

println("\nstep | CT Bernoulli | CT Laplace(∇u) | CT Laplace(λ) | CT KJ | rel(B-∇u) | rel(B-λ) | rel(B-KJ)")
for k in 1:length(t_range)
    cb  = CT_bernoulli[k]
    cmd = CT_laplace_md[k]
    clv = CT_laplace_lv[k]
    ck  = CT_kj[k]
    println("  $k  |  $(round(cb, sigdigits=6))  |  $(round(cmd, sigdigits=6))  |  $(round(clv, sigdigits=6))  |  $(round(ck, sigdigits=6))  |  $(round(relative_difference(cb, cmd), sigdigits=4))  |  $(round(relative_difference(cb, clv), sigdigits=4))  |  $(round(relative_difference(cb, ck), sigdigits=4))")
end

bern_md_rel  = norm(CT_bernoulli - CT_laplace_md) / max(norm(CT_bernoulli), eps())
bern_lv_rel  = norm(CT_bernoulli - CT_laplace_lv) / max(norm(CT_bernoulli), eps())
bern_kj_rel  = norm(CT_bernoulli - CT_kj)         / max(norm(CT_bernoulli), eps())
md_lv_rel    = norm(CT_laplace_md - CT_laplace_lv) / max(norm(CT_laplace_md), eps())

println("\nRelative history differences (L2 norm):")
println("  Bernoulli vs Laplace(∇u):  $(bern_md_rel)")
println("  Bernoulli vs Laplace(λ):   $(bern_lv_rel)")
println("  Bernoulli vs KJ:           $(bern_kj_rel)")
println("  Laplace(∇u) vs Laplace(λ): $(md_lv_rel)")

# Item 005 primary metric: peak-to-peak CT ripple over the final settle window
# (after the freestream is fully withdrawn), with the residual freestream so a
# plateau under nonzero convection is not mistaken for hover.
let
    settle_window_revs = min(settle_revs, required_revs)
    k_start = max(1, length(t_range) - round(Int, nt * settle_window_revs) + 1)
    tail_b  = filter(isfinite, CT_bernoulli[k_start:end])
    tail_md = filter(isfinite, CT_laplace_md[k_start:end])
    tail_lv = filter(isfinite, CT_laplace_lv[k_start:end])
    ptp(v) = isempty(v) ? NaN : maximum(v) - minimum(v)
    mean(v) = isempty(v) ? NaN : sum(v) / length(v)
    residual_magVinf = magVinf_pulse(t_range[end])
    println("\nItem 005 plateau diagnostics (final $(round(settle_window_revs,digits=2)) revs, steps $(k_start):$(length(t_range))):")
    println("  residual magVinf at readout = $(residual_magVinf)  (hover requires ≈ 0)")
    println("  CT Bernoulli   plateau mean=$(round(mean(tail_b), sigdigits=5))  peak-to-peak=$(round(ptp(tail_b), sigdigits=4))")
    println("  CT Laplace(∇u) plateau mean=$(round(mean(tail_md),sigdigits=5))  peak-to-peak=$(round(ptp(tail_md),sigdigits=4))")
    println("  CT Laplace(λ)  plateau mean=$(round(mean(tail_lv),sigdigits=5))  peak-to-peak=$(round(ptp(tail_lv),sigdigits=4))")
end

if save_path !== nothing
    isdir(save_path) || mkpath(save_path)
    csv_path = joinpath(save_path, "$(run_name)_CT_vs_rev.csv")
    open(csv_path, "w") do io
        # Legacy header/columns at NROTORS=1; per-rotor columns appended after
        # the legacy ones for N>1 (CT_bernoulli is then the layout TOTAL).
        extra_cols = nrotors == 1 ? "" :
            join(",CT_bernoulli_r$(i)" for i in 1:nrotors)
        println(io, "step,revolution,CT_bernoulli,CT_laplace_matderiv,CT_laplace_lamb,CT_kj" * extra_cols)
        for k in 1:length(t_range)
            rev = (k - 1) * dt * RPM / 60
            extra = nrotors == 1 ? "" :
                join(",$(-CT_bernoulli_by_rotor[i][k])" for i in 1:nrotors)
            println(io, "$k,$rev,$(-CT_bernoulli[k]),$(-CT_laplace_md[k]),$(-CT_laplace_lv[k]),$(-CT_kj[k])" * extra)
        end
    end
    println("\nWrote CT vs revolution CSV: $csv_path")
end

# --- Phase 2e convergence metric: per-revolution blocks over the final revs ----
# Ryan's criterion: CT settles to a small amplitude with a mean that changes
# little over 5 revolutions. Operationalized as, over the final
# CONVERGENCE_REVS complete revolutions of self-sustained hover:
#   (a) every per-rev mean within CONVERGENCE_MEAN_TOL of their average,
#   (b) every within-rev peak-to-peak <= CONVERGENCE_PTP_TOL of that average,
#   (c) all CT samples in the window finite.
# Thrust convention matches the CSV: CT_thrust = -(axial force channel).
convergence_revs     = parse(Int,     get(ENV, "CONVERGENCE_REVS",     "5"))
convergence_mean_tol = parse(Float64, get(ENV, "CONVERGENCE_MEAN_TOL", "0.005"))
convergence_ptp_tol  = parse(Float64, get(ENV, "CONVERGENCE_PTP_TOL",  "0.02"))

CT_thrust = -CT_bernoulli
nsteps_total = length(t_range)

# ALL complete revolutions, oldest first. The filtered/low-relaxation
# configurations do not plateau -- they settle into a ~9-rev limit cycle -- so
# the whole series is needed to see the cycle and to feed
# examples/analyze_stable_wake_oscillation.jl.
all_block_ranges = Vector{UnitRange{Int}}()
for stop in nt:nt:nsteps_total
    push!(all_block_ranges, (stop - nt + 1):stop)
end
all_block_means = [sum(CT_thrust[r]) / length(r) for r in all_block_ranges]
all_block_ptps  = [maximum(CT_thrust[r]) - minimum(CT_thrust[r]) for r in all_block_ranges]

# The convergence window is the final CONVERGENCE_REVS of those.
block_ranges = length(all_block_ranges) <= convergence_revs ? all_block_ranges :
    all_block_ranges[(end - convergence_revs + 1):end]

block_means = [sum(CT_thrust[r]) / length(r) for r in block_ranges]
block_ptps  = [maximum(CT_thrust[r]) - minimum(CT_thrust[r]) for r in block_ranges]
block_finite = all(r -> all(isfinite, CT_thrust[r]), block_ranges)
window_mean = isempty(block_means) ? NaN : sum(block_means) / length(block_means)
mean_spread_rel = isempty(block_means) ? NaN :
    maximum(abs.(block_means .- window_mean)) / max(abs(window_mean), eps())
ptp_rel = isempty(block_ptps) ? NaN : maximum(block_ptps) / max(abs(window_mean), eps())

# The window must lie entirely after the freestream has been withdrawn, or the
# "plateau" is a forced-convection artifact rather than hover.
hover_start_rev = freestream_ramp_revs + freestream_hold_revs + freestream_withdraw_revs
window_start_rev = isempty(block_ranges) ? NaN : (first(block_ranges[1]) - 1) * dt * RPM / 60
residual_magVinf_end = magVinf_pulse(t_range[end])
window_in_hover = !isnan(window_start_rev) && window_start_rev >= hover_start_rev

converged = block_finite && window_in_hover &&
    mean_spread_rel <= convergence_mean_tol && ptp_rel <= convergence_ptp_tol

# Cycle-mean: the headline number when the wake settles into a limit cycle
# rather than a plateau. Quote as CT = cycle_mean +/- cycle_std over >= 1 period
# (~9-10 revs); the strict per-rev spread test above cannot fire on a limit
# cycle and is retained only as a diagnostic.
cycle_mean = window_mean
cycle_std = length(block_means) < 2 ? NaN :
    sqrt(sum((block_means .- cycle_mean) .^ 2) / (length(block_means) - 1))
cycle_std_rel = isfinite(cycle_std) ? cycle_std / max(abs(cycle_mean), eps()) : NaN

println("\nPhase 2e per-revolution CT blocks (final $(length(block_ranges)) revs, thrust convention):")
println("  rev_block | steps | mean CT | peak-to-peak | rel. dev. from window mean")
for (i, r) in enumerate(block_ranges)
    dev = abs(block_means[i] - window_mean) / max(abs(window_mean), eps())
    println("  $(i) | $(first(r)):$(last(r)) | $(round(block_means[i], sigdigits=6)) | " *
            "$(round(block_ptps[i], sigdigits=4)) | $(round(dev, sigdigits=3))")
end
println("  CYCLE-MEAN CT (headline)    = $(round(cycle_mean, sigdigits=6)) ± $(round(cycle_std, sigdigits=3))" *
        "  (±$(round(100*cycle_std_rel, sigdigits=3))%, over $(length(block_means)) revs)")
println("  max per-rev mean spread     = $(round(mean_spread_rel, sigdigits=3)) (tol $(convergence_mean_tol))")
println("  max within-rev p-p / mean   = $(round(ptp_rel, sigdigits=3)) (tol $(convergence_ptp_tol))")
println("  window starts at rev        = $(round(window_start_rev, digits=2)) (hover begins at rev $(round(hover_start_rev, digits=2)))")
println("  residual magVinf at readout = $(residual_magVinf_end)")
println("  all finite                  = $(block_finite)")
println("  CONVERGED (Phase 2e criterion) = $(converged)")

if save_path !== nothing
    perrev_path = joinpath(save_path, "$(run_name)_CT_per_rev.csv")
    open(perrev_path, "w") do io
        println(io, "rev_block,step_start,step_stop,rev_start,rev_stop,CT_mean,CT_ptp,rel_dev_from_window_mean,in_convergence_window")
        for (i, r) in enumerate(all_block_ranges)
            rev_start = (first(r) - 1) * dt * RPM / 60
            rev_stop  = (last(r)  - 1) * dt * RPM / 60
            dev = abs(all_block_means[i] - window_mean) / max(abs(window_mean), eps())
            in_window = i > length(all_block_ranges) - length(block_ranges)
            println(io, "$i,$(first(r)),$(last(r)),$rev_start,$rev_stop,$(all_block_means[i]),$(all_block_ptps[i]),$dev,$in_window")
        end
    end
    println("Wrote per-revolution block stats: $perrev_path")

    meta_path = joinpath(save_path, "$(run_name)_case_metadata.toml")
    open(meta_path, "w") do io
        println(io, "run_name = \"$(run_name)\"")
        println(io, "formulation = \"$(formulation_name)\"")
        println(io, "green_gauge = \"$(green_gauge)\"")
        println(io, "green_recompute_interval = $(green_recompute_interval)")
        println(io, "mesh_key = \"$(rhpc_mesh)\"")
        println(io, "mesh_file = \"$(msh_file)\"")
        println(io, "te_seed_source = \"$(te_seed_source)\"")
        println(io, "te_indices_1 = $(collect(te_indices_1))")
        println(io, "te_indices_2 = $(collect(te_indices_2))")
        println(io, "ncells = $(size(rotor.cells, 2))")
        println(io, "RPM = $(RPM)")
        println(io, "R = $(R)")
        println(io, "rho = $(rho)")
        println(io, "NT = $(nt)")
        println(io, "dt = $(dt)")
        println(io, "n_steps = $(nsteps_total)")
        println(io, "required_revs = $(required_revs)")
        println(io, "truncation_depth_R = $(cylinder_depth / R)")
        println(io, "nwakerows = $(nwakerows)")
        println(io, "p_per_step = $(p_per_step)")
        println(io, "overlap = $(overlap)")
        println(io, "particle_shedding = \"$(particle_shedding)\"")
        println(io, "conversion = \"$(conversion_mode)\"")
        if conversion_mode == "smooth"
            println(io, "conversion_sigma = $(conversion_sigma)")
            println(io, "conversion_overlap = $(conversion_overlap)")
            println(io, "conversion_attribution = \"$(conversion_attribution)\"")
        end
        println(io, "wake_core_size = $(wake_core_size)")
        println(io, "wake_nu = $(wake_nu)")
        println(io, "wake_core_beta = $(wake_core_beta)")
        println(io, "core_spreading_active = $(core_spreading_active)")
        println(io, "core_spreading_sgm0 = $(core_spreading_sgm0)")
        println(io, "wake_expint = $(wake_expint)")
        println(io, "core_size_panel = $(core_size_panel)")
        println(io, "core_size_targets = $(core_size_targets)")
        println(io, "merge_particles = $(merge_particles)")
        println(io, "relax_rlxf = $(relax_rlxf)")
        println(io, "relax_filter_downstream_R = $(relax_filter_downstream_R)")
        println(io, "particle_relax = $(particle_relax)")
        println(io, "spinup_revs = $(spinup_revs)")
        println(io, "spinup_start_fraction = $(spinup_start_fraction)")
        println(io, "magVinf_peak = $(magVinf_peak)")
        println(io, "magVinf_end = $(magVinf_end)")
        println(io, "freestream_ramp_revs = $(freestream_ramp_revs)")
        println(io, "freestream_hold_revs = $(freestream_hold_revs)")
        println(io, "freestream_withdraw_revs = $(freestream_withdraw_revs)")
        println(io, "settle_revs = $(settle_revs)")
        println(io, "bernoulli_only = $(bernoulli_only)")
        println(io, "restart_step = $(restart_step)")
        println(io, "backend_body_order = $(backend.expansion_order)")
        println(io, "backend_wake_order = $(backend_wake.expansion_order)")
        println(io, "julia_threads = $(Threads.nthreads())")
        println(io, "wall_time_s = $(sim_wall_seconds)")
        println(io, "convergence_revs = $(length(block_ranges))")
        println(io, "convergence_mean_tol = $(convergence_mean_tol)")
        println(io, "convergence_ptp_tol = $(convergence_ptp_tol)")
        println(io, "das_eta_kinematic = $(init_Das_eta_kinematic)")
        println(io, "das_chord_fraction = $(das_chord_fraction)")
        println(io, "das_uniform_dsigma = $(das_uniform_dsigma)")
        println(io, "sigma_chord_fraction = $(sigma_chord_fraction)")
        println(io, "sigma_floor_r = $(sigma_floor_r)")
        println(io, "das_sigma_lambda = $(das_sigma_lambda)")
        println(io, "das_curvature_beta = $(das_curvature_beta)")
        println(io, "das_arc_placed = $(das_arc_placed)")
        println(io, "das_arc_source = $(repr(das_arc_source))")
        println(io, "das_arc_table = $(repr(basename(das_arc_table)))")
        println(io, "theta_max = $(theta_max_cs)")
        println(io, "das_min_displacement = $(set_Das_min_kinematic_displacement)")
        println(io, "das_kinematic_arc = $(set_Das_kinematic_arc)")
        println(io, "das_refresh = $(set_Das_refresh)")
        println(io, "shed_with_induced_velocity = $(shed_with_induced_velocity)")
        println(io, "n_revs_recorded = $(length(all_block_ranges))")
        println(io, "CT_cycle_mean = $(cycle_mean)")
        println(io, "CT_cycle_std = $(cycle_std)")
        println(io, "CT_cycle_std_rel = $(cycle_std_rel)")
        println(io, "CT_cycle_revs = $(length(block_means))")
        println(io, "CT_window_mean = $(window_mean)")
        println(io, "CT_mean_spread_rel = $(mean_spread_rel)")
        println(io, "CT_ptp_rel = $(ptp_rel)")
        println(io, "window_start_rev = $(window_start_rev)")
        println(io, "hover_start_rev = $(hover_start_rev)")
        println(io, "window_in_hover = $(window_in_hover)")
        println(io, "residual_magVinf = $(residual_magVinf_end)")
        println(io, "all_finite = $(block_finite)")
        println(io, "converged = $(converged)")
        # BRAINSTORM/022 ground keys
        println(io, "ground_enable = $(ground_enable)")
        println(io, "ground_h_r = $(ground_h_r)")
        println(io, "ground_radius_r = $(ground_radius_r)")
        println(io, "ground_panel_length_r = $(ground_panel_length_r)")
        println(io, "ground_ncells = $(ground_enable ? ground.ncells : 0)")
        println(io, "ground_x = $(ground_x)")
        println(io, "x_rotor_plane = $(x_rotor_plane)")
        println(io, "ground_particle_policy = \"$(ground_particle_policy)\"")
        println(io, "ground_damp_band_r = $(ground_damp_band_r)")
        println(io, "ground_damp_total = $(ground_damp_total[])")
        println(io, "effective_cylinder_depth_R = $(effective_cylinder_depth / R)")
        println(io, "trunc_radius_r = $(trunc_radius_r)")
        println(io, "gs_log = $(gs_log)")
        println(io, "gs_max_outer = $(gs_max_outer)")
        println(io, "gs_tol = $(gs_tol)")
        println(io, "gs_iters_max = $(gs_iters_max)")
        println(io, "gs_iters_mean = $(gs_iters_mean)")
        println(io, "gs_nonconverged = $(gs_nonconverged)")
        # 022 multi-rotor keys
        println(io, "nrotors = $(nrotors)")
        println(io, "rotor_spacing_r = $(rotor_spacing_r)")
        println(io, "rotor_directions = $(rotor_directions)")
        println(io, "rotor_centers_r = [" *
            join(("[" * join(round.(c ./ R, digits=6), ", ") * "]" for c in rotor_centers), ", ") * "]")
        println(io, "max_rotor_offset_r = $(max_rotor_offset / R)")
        println(io, "ground_radius_eff_r = $(ground_radius_eff / R)")
        if nrotors > 1
            for (i, ct) in enumerate(CT_bernoulli_by_rotor)
                ct_thrust_i = -ct
                wmean = isempty(block_ranges) ? NaN :
                    sum(sum(ct_thrust_i[r]) / length(r) for r in block_ranges) / length(block_ranges)
                println(io, "CT_window_mean_r$(i) = $(wmean)")
            end
        end
        if ground_enable && run_monitors && !isempty(ground_tangency_rms)
            println(io, "ground_tangency_rms_final = $(ground_tangency_rms[end])")
            println(io, "ground_tangency_rms_max = $(maximum(ground_tangency_rms))")
            println(io, "below_ground_count_final = $(ground_below_count[end])")
            println(io, "below_ground_count_max = $(maximum(ground_below_count))")
        end
    end
    # Julia prints non-finite floats as NaN/Inf/-Inf, which are invalid bare
    # TOML values; rewrite them as the TOML-valid lowercase floats
    write(meta_path, replace(read(meta_path, String),
        r"= NaN$"m => "= nan", r"= Inf$"m => "= inf", r"= -Inf$"m => "= -inf"))
    println("Wrote case metadata: $meta_path")
end

end # run_monitors

end # !rhpc_setup_only
