## Minimal reproducer for the FastMultipole OOM seen in rotor_hover_pressure_comparison.jl.
##
## Sets up the same rotor + PanelWake state as step 2 of that example, then issues
## a single FastMultipole.fmm! call with the same probes/sources as
## KuttaJoukowskiForce. Designed for fast iteration when bisecting which phase
## of fmm! (tree build, nearfield, farfield) allocates pathologically.
##
## Env vars:
##   FMM_NEARFIELD=true|false  (default true)
##   FMM_FARFIELD=true|false   (default true)
##   FMM_EXPANSION_ORDER       (default 10)
##   FMM_ACCEPTANCE            (default 0.4)
##   FMM_LEAF_SIZE             (default 20)
##   FMM_SHRINK=true|false     (default true)
##   REPRO_SOURCES             "all" (default), "rotor", "wake", "filament",
##                             "rotor+wake", "rotor+filament", "wake+filament"

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
using FLOWPanel.FastMultipole
using FLOWPanel.FastMultipole.StaticArrays
using VSPGeom
import GeoIO
using LinearAlgebra: norm

# ----- match the example geometry/kinematics exactly -----

magVinf = 0.0001
AOA     = 0.0
rho     = 1.225
RPM     = 6000
Vinf    = magVinf * [0.0, -cosd(AOA), sind(AOA)]
R       = 0.119
nrevs   = 1
nt      = 36
dt      = 60 / RPM / nt
n_steps = nt * nrevs
t_range = range(0.0, step=dt, length=n_steps)

cp_outer=true
kerneloffset   = R * 1e-3
kernelcutoff   = R * 1e-13
wake_core_size = 1e-3

read_path = joinpath(pnl.examples_path, "data")
_use_phantom = get(ENV, "REPRO_MESH", "phantom") == "phantom"
if _use_phantom
    msh_file     = joinpath(read_path, "phantom_3_rebuild_r2.msh")
    te_indices_1 = [9, 175, 127]
    te_indices_2 = [13, 286, 238]
else
    msh_file     = joinpath(read_path, "dji9443_new_40_40.msh")
    te_indices_1 = [1614, 1574, 0]    .+ 1
    te_indices_2 = [3323, 3284, 1711] .+ 1
end

msh = GeoIO.load(msh_file).geometry
nodes, cells = pnl.meshes2nodes_cells(msh)
nodes .*= R / maximum(nodes[1, :])

shedding = pnl.noshedding
kernel   = Union{pnl.ConstantSource, pnl.VortexRing}
DBC      = kernel == pnl.VortexRing ? false : true

rotor = pnl.RigidWakeBody{kernel}(nodes, cells, shedding;
    kerneloffset, kernelcutoff,
    semiinfinite_wake=false, watertight=true, DBC)

shedding1 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_1[1], te_indices_1[2];
    bbox=nothing, end_node=te_indices_1[3], normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
shedding2 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_2[1], te_indices_2[2];
    bbox=nothing, end_node=te_indices_2[3], normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

rotor = pnl.RigidWakeBody{kernel}(rotor.nodes, rotor.cells, [shedding1, shedding2];
    kerneloffset, kernelcutoff,
    semiinfinite_wake=false, watertight=true,
    ensure_winding=true, DBC)

wake_rotor = pnl.PanelWake(rotor; nwakerows=12, core_size=wake_core_size)

Uinf(t) = Vinf

frames = pnl.ReferenceFrame(rotor;
    origin       = SVector{3}(0.0, 0.0, 0.0),
    v            = SVector{3}(0.0, 0.0, 0.0),
    ω_axis       = SVector{3}(0.0, 1.0, 0.0),
    ω            = 2 * pi * RPM / 60,
    R            = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    name         = "vehicle",
    child_index  = Int[],
    dependent_index = [1])

pnl.initialize_Das!((rotor,), frames, Uinf, t_range[1], t_range[2] - t_range[1];
    set_Das_eta_kinematic = 0.2)

solver_rotor = pnl.Backslash(rotor)
backend      = pnl.FastMultipoleBackend()  # matches the example
maneuver!(frames, systems, wakes, t) = nothing

# ----- advance through one full step so the wake has one shed row -----
println("Running one step of simulate! (no monitors) to populate wake to nwakes[]=1...")
@time pnl.simulate!((rotor,), (wake_rotor,), frames, maneuver!, Uinf, t_range[1:2];
    set_Das_eta_kinematic = NaN,
    monitors   = (),
    body_solvers = (solver_rotor,),
    backend, verbose=true,
    path = nothing)  # skip VTK
println("wake_rotor.nwakes[] = ", wake_rotor.nwakes[])

# At this point the next thing simulate! would do is step 1's main fmm! calls.
# We mimic the call KuttaJoukowskiForce would make at the top of step 1.

# Build the probe set at edge midpoints, same as KuttaJoukowskiForce constructor
# (src/FLOWPanel_simulate_monitors.jl:1115-1155).
edge_node_a = Int[]
edge_node_b = Int[]
for shed in rotor.shedding
    for k in axes(shed, 2)
        i_panel = shed[1, k]
        idx_a   = shed[2, k]
        idx_b   = shed[3, k]
        push!(edge_node_a, rotor.cells[idx_a, i_panel])
        push!(edge_node_b, rotor.cells[idx_b, i_panel])
    end
end
n_probes = length(edge_node_a)
probes   = FastMultipole.ProbeSystem(n_probes, Float64)
@inbounds for k in 1:n_probes
    a, b = edge_node_a[k], edge_node_b[k]
    probes.position[k] = SVector{3, Float64}(
        0.5 * (rotor.nodes[1, a] + rotor.nodes[1, b]),
        0.5 * (rotor.nodes[2, a] + rotor.nodes[2, b]),
        0.5 * (rotor.nodes[3, a] + rotor.nodes[3, b]),
    )
end
println("n_probes = $n_probes")

# Assemble sources per env var.
wake_sources_all = pnl.get_sources(wake_rotor)  # = (wake_rotor, FilamentWrapper(wake_rotor))
println("wake_sources length = $(length(wake_sources_all)); types = $(typeof.(wake_sources_all))")

sources_key = get(ENV, "REPRO_SOURCES", "all")
all_sources = if sources_key == "rotor"
    (rotor,)
elseif sources_key == "wake"
    (wake_sources_all[1],)
elseif sources_key == "filament"
    (wake_sources_all[2],)
elseif sources_key == "rotor+wake"
    (rotor, wake_sources_all[1])
elseif sources_key == "rotor+filament"
    (rotor, wake_sources_all[2])
elseif sources_key == "wake+filament"
    wake_sources_all
else
    (rotor, wake_sources_all...)
end
println("Using sources ($sources_key): $(typeof.(all_sources))")

# FMM parameters.
expansion_order      = parse(Int,     get(ENV, "FMM_EXPANSION_ORDER", "10"))
multipole_acceptance = parse(Float64, get(ENV, "FMM_ACCEPTANCE",      "0.4"))
leaf_size_source     = parse(Int,     get(ENV, "FMM_LEAF_SIZE",       "20"))
shrink               = parse(Bool,    get(ENV, "FMM_SHRINK",          "true"))
nearfield            = parse(Bool,    get(ENV, "FMM_NEARFIELD",       "true"))
farfield             = parse(Bool,    get(ENV, "FMM_FARFIELD",        "true"))

println("FMM kwargs: expansion_order=$expansion_order multipole_acceptance=$multipole_acceptance leaf_size_source=$leaf_size_source shrink=$shrink nearfield=$nearfield farfield=$farfield")

# Make sure rotor's source-side state is fresh (pre-evaluate, as
# FLOWPanel.influence! does when precalc=true).
pnl.pre_evaluate_influence!(rotor)

println("\n--- Issuing FastMultipole.fmm! ---")
@time FastMultipole.fmm!((probes,), all_sources;
    expansion_order,
    multipole_acceptance,
    leaf_size_source,
    scalar_potential = false,
    gradient         = true,
    hessian          = (false,),
    extra_farfield   = false,
    shrink,
    nearfield,
    farfield,
)
println("--- fmm! returned successfully ---")
