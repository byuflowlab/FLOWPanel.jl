## Rotor hover: compare PressureBernoulli, PressureLaplace, and KuttaJoukowskiForce
## Short unsteady run intended as a pressure/force cross-check.

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
using FLOWPanel.FastMultipole.StaticArrays
using VSPGeom
import GeoIO
using LinearAlgebra: norm

magVinf = 0.0001
AOA = 0.0
rho = 1.225
RPM = 6000
Vinf = magVinf * [0.0, -cosd(AOA), sind(AOA)]
R = 0.119
nrevs = 1
nt = 36
dt = 60 / RPM / nt
n_steps = nt * nrevs
t_range = range(0.0, step=dt, length=n_steps)[1:12]

CPoffset = R * 1e-6
kerneloffset = R * 1e-3
kernelcutoff = R * 1e-13
p_per_step = 2
overlap = 2.0
merge_r_factor = 0.02
merge_r_hash_factor = 0.04
init_Das_eta_kinematic = 0.2
p_correct_kuttacondition_flag = false
wake_core_size = parse(Float64, get(ENV, "WAKE_CORE_SIZE", "1e-3"))

read_path = joinpath(pnl.examples_path, "data")
msh_file = joinpath(read_path, "phantom_3_rebuild_r2.msh")
te_indices_1 = [9, 175, 127]
te_indices_2 = [13, 286, 238]

msh = GeoIO.load(msh_file).geometry
nodes, cells = pnl.meshes2nodes_cells(msh)
nodes .*= R / maximum(nodes[1, :])

shedding = pnl.noshedding
kernel = Union{pnl.ConstantSource, pnl.VortexRing}
DBC = kernel == pnl.VortexRing ? false : true
rotor = pnl.RigidWakeBody{kernel}(nodes, cells, shedding;
    CPoffset, kerneloffset, kernelcutoff,
    semiinfinite_wake=false, watertight=true, DBC)

shedding1 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_1[1], te_indices_1[2];
    bbox=nothing, end_node=te_indices_1[3], normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
shedding2 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_2[1], te_indices_2[2];
    bbox=nothing, end_node=te_indices_2[3], normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

rotor = pnl.RigidWakeBody{kernel}(rotor.nodes, rotor.cells, [shedding1, shedding2];
    CPoffset, kerneloffset, kernelcutoff,
    semiinfinite_wake=false, watertight=true,
    ensure_winding=true, DBC)

wake_rotor = pnl.PanelWake(rotor; nwakerows=12, core_size=wake_core_size)

Uinf(t) = Vinf

frames = pnl.ReferenceFrame(rotor;
    origin=SVector{3}(0.0, 0.0, 0.0),
    v=SVector{3}(0.0, 0.0, 0.0),
    ω_axis=SVector{3}(0.0, 1.0, 0.0),
    ω=2 * pi * RPM / 60,
    R=SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    name="vehicle",
    child_index=Int[],
    dependent_index=[1]
)

pnl.initialize_Das!((rotor,), frames, Uinf, t_range[1], t_range[2] - t_range[1];
    set_Das_eta_kinematic=init_Das_eta_kinematic)

solver_rotor = pnl.Backslash(rotor)
backend = pnl.FastMultipoleBackend()
maneuver!(frames, systems, wakes, t) = nothing

systems = (rotor,)
wakes = (wake_rotor,)
body_solvers = (solver_rotor,)

pressure_bernoulli = pnl.PressureBernoulli(rho;
    unsteady=true,
    correct_kuttacondition=p_correct_kuttacondition_flag,
    backend=backend)
force_monitor_bernoulli = pnl.ForceMonitor(length(t_range), 1;
    i_frame=1,
    normalization=pnl.RotorNormalization(rho, 2 * R, 1),
    correct_kuttacondition=p_correct_kuttacondition_flag,
    verbose=false)

pressure_laplace = pnl.PressureLaplace(rotor, rho; verbose=false)
force_monitor_laplace = pnl.ForceMonitor(length(t_range), 1;
    i_frame=1,
    normalization=pnl.RotorNormalization(rho, 2 * R, 1),
    correct_kuttacondition=p_correct_kuttacondition_flag,
    verbose=false)

kj_monitor = pnl.KuttaJoukowskiForce(rotor, length(t_range), 1;
    rho, backend,
    i_frame=1,
    normalization=pnl.RotorNormalization(rho, 2 * R, 1),
    verbose=false)

monitors = (
    pressure_laplace,
    force_monitor_laplace,
    pressure_bernoulli,
    force_monitor_bernoulli,
    kj_monitor,
)

println("\nBegin rotor hover pressure comparison ($(length(t_range)) steps)...")
name = "rotor_hover_pressure_comparison"
@time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
    set_Das_eta_kinematic=NaN,
    monitors,
    body_solvers, backend, verbose=true,
    path=name, name,
)

CT_bernoulli = force_monitor_bernoulli.force[1, :]
CT_laplace = force_monitor_laplace.force[1, :]
CT_kj = kj_monitor.force[1, :]

function relative_difference(a, b)
    denom = max(abs(b), eps())
    return abs(a - b) / denom
end

println("\nstep | CT Bernoulli | CT Laplace | CT KJ | rel(B-L) | rel(B-KJ)")
for k in 1:length(t_range)
    cb = CT_bernoulli[k]
    cl = CT_laplace[k]
    ck = CT_kj[k]
    println("  $k  |  $(round(cb, sigdigits=6))  |  $(round(cl, sigdigits=6))  |  $(round(ck, sigdigits=6))  |  $(round(relative_difference(cb, cl), sigdigits=4))  |  $(round(relative_difference(cb, ck), sigdigits=4))")
end

bern_lap_rel = norm(CT_bernoulli - CT_laplace) / max(norm(CT_bernoulli), eps())
bern_kj_rel = norm(CT_bernoulli - CT_kj) / max(norm(CT_bernoulli), eps())

println("\nRelative history differences:")
println("  Bernoulli vs Laplace: $(bern_lap_rel)")
println("  Bernoulli vs KuttaJoukowski: $(bern_kj_rel)")
