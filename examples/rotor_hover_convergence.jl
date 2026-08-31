## Rotor hover mesh-convergence helper.
##
## Iterative workflow (mesh generation and refinement done by the user):
##   1. Generate a .msh file for the rotor and point MSH_FILE at it.
##   2. Run without TE indices -> writes data/<run_name>/rotor_initial VTK to
##      inspect in ParaView and pick trailing-edge seed nodes for calc_shedding.
##   3. Rerun with TE_INDICES_1/TE_INDICES_2 (0-based ParaView point IDs)
##      -> single steady! evaluation with a rotor thrust-coefficient monitor.
##   4. Refine the mesh and repeat; CT history accumulates in
##      data/<run_name>/convergence_history.csv.
##
## Example:
##   MSH_FILE=dji9443_new_40_40.msh \
##   TE_INDICES_1="1614,1574,45" TE_INDICES_2="3324,3284,1755" \
##   julia --project examples/rotor_hover_convergence.jl

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
using LinearAlgebra: norm

run_name = get(ENV, "RUN_NAME", "rotor_hover_convergence")
save_path = joinpath("data", run_name)

msh_input = get(ENV, "MSH_FILE", "")
isempty(msh_input) && error("""
    MSH_FILE is required. Provide the rotor mesh, e.g.:
        MSH_FILE=dji9443_new_40_40.msh julia --project examples/rotor_hover_convergence.jl
    Paths are resolved relative to $(joinpath(pnl.examples_path, "data")) unless they exist as given.""")
msh_file = isfile(msh_input) ? msh_input : joinpath(pnl.examples_path, "data", msh_input)
isfile(msh_file) || error("Mesh file not found: $(msh_input) (also tried $(msh_file))")
mesh_tag = splitext(basename(msh_file))[1]

parse_te_indices(str) = [parse(Int, s) for s in split(str, ",")] .+ 1 # 0-based ParaView IDs to 1-based
te_str_1 = get(ENV, "TE_INDICES_1", "")
te_str_2 = get(ENV, "TE_INDICES_2", "")
inspect_only = isempty(te_str_1) || isempty(te_str_2)

magVinf = 0.0001
AOA = 0.0
rho = 1.179 # from NASA paper
RPM = parse(Float64, get(ENV, "RPM", "5400"))
R = 0.119
shedding_r_over_R = parse(Float64, get(ENV, "SHEDDING_R_OVER_R", "0.1"))
0.0 <= shedding_r_over_R <= 1.0 || error("shedding_r_over_R must be between 0 and 1")

core_size_panel = parse(Float64, get(ENV, "CORE_SIZE_PANEL", get(ENV, "KERNELOFFSET_PANEL", string(R * 1e-10))))
core_size_targets = parse(Float64, get(ENV, "CORE_SIZE_TARGETS", get(ENV, "KERNELOFFSET_TARGETS", get(ENV, "CORE_SIZE", get(ENV, "KERNELOFFSET", "1e-3")))))
kernelcutoff = R * 1e-13
init_Das_eta_kinematic = 0.2
set_Das_min_kinematic_displacement = 0.01 * R

# Surface-velocity (∇μ) reconstruction options, built from env vars. A key is
# included only when its env var is set, so anything unset falls back to the
# production default inside steady!. Valid keys match _normalize_grad_mu_options
# in src/FLOWPanel_postprocess.jl; bad values throw a clear ArgumentError there.
function grad_mu_options_from_env()
    opts = NamedTuple()
    haskey(ENV, "GRAD_MU_BASIS")               && (opts = merge(opts, (; basis = Symbol(ENV["GRAD_MU_BASIS"]))))
    haskey(ENV, "GRAD_MU_TRI_ROBUST")          && (opts = merge(opts, (; tri_robust = parse(Bool, ENV["GRAD_MU_TRI_ROBUST"]))))
    haskey(ENV, "GRAD_MU_QUAD_GROW")           && (opts = merge(opts, (; quad_grow = parse(Bool, ENV["GRAD_MU_QUAD_GROW"]))))
    haskey(ENV, "GRAD_MU_QUAD_GROW_STOP")      && (opts = merge(opts, (; quad_grow_stop = Symbol(ENV["GRAD_MU_QUAD_GROW_STOP"]))))
    haskey(ENV, "GRAD_MU_QUAD_GROW_MAX_DEPTH") && (opts = merge(opts, (; quad_grow_max_depth = parse(Int, ENV["GRAD_MU_QUAD_GROW_MAX_DEPTH"]))))
    haskey(ENV, "GRAD_MU_QUAD_GROW_COND_MAX")  && (opts = merge(opts, (; quad_grow_cond_max = parse(Float64, ENV["GRAD_MU_QUAD_GROW_COND_MAX"]))))
    return opts
end
grad_mu_options = grad_mu_options_from_env()
grad_mu_tag = haskey(grad_mu_options, :basis) ? string(grad_mu_options.basis) : "default"

axial_dimension = occursin("dji9443", msh_file) ? 1 : 2 # DJI9443 geometry is rotated compared to typical rotor convention
radial_dimension = occursin("dji9443", msh_file) ? 2 : 1 # this might be wrong for non-dji9443
omega_axis = occursin("dji9443", msh_file) ? SVector{3}(-1.0, 0.0, 0.0) : SVector{3}(0.0, 1.0, 0.0)
Vinf_direction = occursin("dji9443", msh_file) ? [cosd(AOA), sind(AOA), 0.0] : [0.0, -cosd(AOA), sind(AOA)]
Vinf = magVinf * Vinf_direction
Uinf(t) = Vinf

msh = pnl.read_gmsh(msh_file)
nodes, cells = pnl.meshes2nodes_cells(msh)
nodes .*= R / maximum(nodes[radial_dimension, :])

kernel = Union{pnl.ConstantSource, pnl.VortexRing}
DBC = kernel == pnl.VortexRing ? false : true
rotor = pnl.RigidWakeBody{kernel}(nodes, cells, pnl.noshedding;
    core_size=core_size_panel, core_size_panel, core_size_targets, kernelcutoff,
    semiinfinite_wake=false, watertight=true, DBC)

isdir(save_path) || mkpath(save_path)

if inspect_only
    vtk_path = joinpath(save_path, "rotor_initial")
    pnl.write_vtk(vtk_path, rotor)
    println("""
        Mesh $(mesh_tag): $(size(rotor.nodes, 2)) nodes, $(rotor.ncells) cells.
        Wrote inspection VTK: $(vtk_path)

        Open it in ParaView and, for each blade, pick the trailing-edge seed:
        two adjacent TE point IDs (tip-most first) plus the root-end point ID.
        Then rerun with the 0-based ParaView point IDs, e.g.:
            MSH_FILE=$(msh_input) TE_INDICES_1="i,j,k" TE_INDICES_2="i,j,k" \\
            julia --project examples/rotor_hover_convergence.jl""")
    exit()
end

te_indices_1 = parse_te_indices(te_str_1)
te_indices_2 = parse_te_indices(te_str_2)

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

bbox1 = make_shedding_bbox(rotor.nodes, te_indices_1[1:2], radial_dimension, R, shedding_r_over_R)
shedding1 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_1[1], te_indices_1[2];
    bbox=bbox1, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
bbox2 = make_shedding_bbox(rotor.nodes, te_indices_2[1:2], radial_dimension, R, shedding_r_over_R)
shedding2 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_2[1], te_indices_2[2];
    bbox=bbox2, normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

rotor = pnl.RigidWakeBody{kernel}(rotor.nodes, rotor.cells, [shedding1, shedding2];
    core_size=core_size_panel, core_size_panel, core_size_targets, kernelcutoff,
    semiinfinite_wake=false, watertight=true,
    ensure_winding=true, DBC)

function shedding_root_r_over_R(nodes, shedding, cells, radial_dimension, R)
    isempty(shedding) && return NaN
    root_edge = shedding[:, end]
    pi, nia, nib = root_edge[1], root_edge[2], root_edge[3]
    edge_nodes = cells[[nia, nib], pi]
    midpoint = (nodes[:, edge_nodes[1]] + nodes[:, edge_nodes[2]]) / 2
    return midpoint[radial_dimension] / R
end

println("Requested shedding root at |r/R| >= $(shedding_r_over_R)")
println("  shedding1 root midpoint r/R = $(shedding_root_r_over_R(rotor.nodes, shedding1, rotor.cells, radial_dimension, R))")
println("  shedding2 root midpoint r/R = $(shedding_root_r_over_R(rotor.nodes, shedding2, rotor.cells, radial_dimension, R))")

omega = 2 * pi * RPM / 60
frames = pnl.ReferenceFrame(rotor;
    origin=SVector{3}(0.0, 0.0, 0.0),
    v=SVector{3}(0.0, 0.0, 0.0),
    ω_axis=omega_axis,
    ω=omega,
    R=SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    name="vehicle",
    child_index=Int[],
    dependent_index=[1]
)

dt = 60 / RPM / 36 # nominal step for Das initialization
pnl.initialize_Das!((rotor,), frames, Uinf, 0.0, dt;
    set_Das_eta_kinematic=init_Das_eta_kinematic,
    set_Das_min_kinematic_displacement)

backend = pnl.FastMultipoleBackend(;
    expansion_order=parse(Int, get(ENV, "FMM_EXPANSION_ORDER", "8")),
    multipole_acceptance=parse(Float64, get(ENV, "FMM_ACCEPTANCE", "0.4")),
    leaf_size=parse(Int, get(ENV, "FMM_LEAF_SIZE", "20")),
)

pressure_bernoulli = pnl.PressureBernoulli(rho;
    unsteady=false,
    correct_kuttacondition=false,
    backend=backend)
force_monitor = pnl.ForceMonitor(1, 1;
    i_frame=1,
    normalization=pnl.RotorNormalization(rho, 2 * R, 1),
    correct_kuttacondition=false,
    verbose=true)

println("\nSteady solve of $(mesh_tag) ($(rotor.ncells) cells) at $(RPM) RPM...")
println("  grad_mu_options = $(isempty(grad_mu_options) ? "(default)" : grad_mu_options); core_size_targets = $(core_size_targets)")
@time pnl.steady!((rotor,), frames, Vinf;
    body_solvers=(pnl.Backslash(rotor),),
    backend,
    monitors=(pressure_bernoulli, force_monitor),
    grad_mu_options,
    path=save_path,
    name=run_name * "_" * mesh_tag,
    verbose=true)

CT = -force_monitor.force[axial_dimension, 1]
println("\nCT (Bernoulli, steady) = $(CT)")

# Accumulate the convergence history across runs, one row per mesh.
csv_path = joinpath(save_path, "convergence_history.csv")
header = "mesh,nnodes,ncells,CT,grad_mu_basis,core_size_targets,RPM"
rows = isfile(csv_path) ? readlines(csv_path)[2:end] : String[]
filter!(row -> split(row, ",")[1] != mesh_tag, rows)
push!(rows, "$(mesh_tag),$(size(rotor.nodes, 2)),$(rotor.ncells),$(CT),$(grad_mu_tag),$(core_size_targets),$(RPM)")
sort!(rows, by=row -> parse(Int, split(row, ",")[3]))
open(csv_path, "w") do io
    println(io, header)
    foreach(row -> println(io, row), rows)
end

println("\nConvergence history ($(csv_path)):")
println("  " * header)
foreach(row -> println("  " * row), rows)
