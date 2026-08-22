## DJI 9443 direct-potential finite panel-wake study driver.
##
## Standalone driver for the `DirectWakePotential` formulation on the DJI 9443
## rotor (see plans/20260721_panel_wake.md). Unlike the particle-wake pressure
## comparison, this uses a strictly finite scalar-potential `PanelWake`
## (include_final_filament=false) and a source+ConstantDoublet body so every
## wake source exposes a scalar potential. The only schedule is a single
## monotone freestream *decrease* at constant RPM — no spin-up, no up-ramp.
##
## Configuration is entirely through env vars (see the parameter block). NREVS
## is DERIVED so the rolling wake-row buffer fills under the benign starting
## freestream, the freestream eases down, and the final revolution is settled:
##   NREVS = ceil(WAKE_ROWS/NT) + FREESTREAM_DECREASE_REVS + SETTLE_REVS
##
## Usage (smoke):
##   RUN_NAME=dji9443_panelwake_smoke RHPC_MESH=40_40 RPM=5400 NT=12 WAKE_ROWS=3 \
##   SETTLE_REVS=0 START_VINF=4.0 TERMINAL_VINF=4.0 FREESTREAM_DECREASE_REVS=0 \
##   FORMULATION=direct SAVE_VTK=true \
##   julia --project=. -t 4 examples/rotor_panel_wake_study.jl

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
using LinearAlgebra: norm
import LinearAlgebra

# BLAS threads: honor BLAS_NUM_THREADS, else OMP_NUM_THREADS, else Julia's count.
let n = tryparse(Int, get(ENV, "BLAS_NUM_THREADS", get(ENV, "OMP_NUM_THREADS", string(Threads.nthreads()))))
    if n !== nothing && n > 0
        LinearAlgebra.BLAS.set_num_threads(n)
    end
    println("BLAS threads: $(LinearAlgebra.BLAS.get_num_threads()) (Julia threads: $(Threads.nthreads()))")
end

# ---------------------------------------------------------------- parameters --
run_name  = get(ENV, "RUN_NAME", "dji9443_panelwake")
save_vtk   = parse(Bool, get(ENV, "SAVE_VTK", "true"))
save_path  = joinpath("data", run_name)   # always used for CSV/metadata/marker
# RUN_NAME may contain "/" to nest the save dir; filenames (simulate! name,
# CSV/metadata prefixes) must stay slash-free or paths nest again inside save_path.
base_name  = basename(run_name)

rho = 1.179                      # NASA paper
R   = 0.119
AOA = 0.0
shedding_r_over_R = 0.1

RPM = parse(Float64, get(ENV, "RPM", "5400"))
nt  = parse(Int,     get(ENV, "NT",  "36"))
dt  = 60 / RPM / nt

wake_rows = parse(Int,     get(ENV, "WAKE_ROWS", "72"))
settle_revs             = parse(Float64, get(ENV, "SETTLE_REVS",             "3"))
start_vinf              = parse(Float64, get(ENV, "START_VINF",              "4.0"))
terminal_vinf           = parse(Float64, get(ENV, "TERMINAL_VINF",           "4.0"))
freestream_decrease_revs = parse(Float64, get(ENV, "FREESTREAM_DECREASE_REVS", "3"))

# NREVS is derived; allow an explicit override for debugging only.
fill_revs = ceil(wake_rows / nt)
derived_nrevs = fill_revs + freestream_decrease_revs + settle_revs
nrevs = parse(Float64, get(ENV, "NREVS", string(derived_nrevs)))

n_steps = round(Int, nt * nrevs)
t_range = range(0.0, step=dt, length=n_steps)
sec_per_rev = 60 / RPM

formulation_name = lowercase(get(ENV, "FORMULATION", "direct"))
recompute_interval = parse(Int, get(ENV, "RECOMPUTE_INTERVAL", "1"))

kerneloffset_panel   = parse(Float64, get(ENV, "KERNELOFFSET_PANEL", string(R * 1e-10)))
kerneloffset_targets = parse(Float64, get(ENV, "KERNELOFFSET_TARGETS", get(ENV, "KERNELOFFSET", "1e-3")))
kernelcutoff = R * 1e-13
wake_core_size = parse(Float64, get(ENV, "WAKE_CORE_SIZE", "1e-3"))
nbins = parse(Int, get(ENV, "NBINS", "30"))

read_path = joinpath(pnl.examples_path, "data")
rhpc_mesh = lowercase(get(ENV, "RHPC_MESH", "40_40"))
if rhpc_mesh == "40_40"
    msh_file = joinpath(read_path, "dji9443_new_40_40.msh")
    te_indices_1 = [1614, 1574, 45]   .+ 1   # 0-based ParaView IDs -> 1-based Julia
    te_indices_2 = [3324, 3284, 1755] .+ 1
else
    error("Unknown RHPC_MESH=$(repr(rhpc_mesh)); this study supports only 40_40")
end

axial_dimension  = 1   # DJI9443 geometry rotated vs typical rotor convention
radial_dimension = 2
Vinf_direction   = [cosd(AOA), sind(AOA), 0.0]

println("""
DJI 9443 direct-potential panel-wake study
  run_name=$(run_name)  formulation=$(formulation_name)  mesh=$(rhpc_mesh)
  RPM=$(RPM)  NT=$(nt)  WAKE_ROWS=$(wake_rows)  core_size=$(wake_core_size)
  START_VINF=$(start_vinf)  TERMINAL_VINF=$(terminal_vinf)  DECREASE_REVS=$(freestream_decrease_revs)
  SETTLE_REVS=$(settle_revs)  fill_revs=$(fill_revs)  NREVS=$(nrevs)  steps=$(n_steps)
""")

# ---------------------------------------------------------- winding-safe body --
# Source + ConstantDoublet (Dirichlet) so the wake kernel is ConstantDoublet and
# every wake source exposes a scalar potential (required by DirectWakePotential).
kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}
DBC = true

msh = pnl.read_gmsh(msh_file)
nodes, cells = pnl.meshes2nodes_cells(msh)
nodes .*= R / maximum(nodes[radial_dimension, :])

# (1) no-shedding body so the constructor re-winds cells; compute shedding from
# ITS rewound .nodes/.cells (CLAUDE.md critical invariant).
rotor = pnl.RigidWakeBody{kernel}(nodes, cells, pnl.noshedding;
    kerneloffset=kerneloffset_panel, kerneloffset_panel, kerneloffset_targets, kernelcutoff,
    semiinfinite_wake=false, watertight=true, DBC)

function make_shedding_bbox(nodes, seed_nodes, radial_dimension, R, shedding_r_over_R)
    radial_midpoint = sum(nodes[radial_dimension, seed_nodes]) / length(seed_nodes)
    radial_sign = sign(radial_midpoint)
    radial_sign == 0 && error("Seed edge lies on the rotor axis; cannot determine shedding side")
    lower = [minimum(nodes[i, :]) for i in 1:size(nodes, 1)]
    upper = [maximum(nodes[i, :]) for i in 1:size(nodes, 1)]
    padding = max(sqrt(eps(eltype(nodes))) * R, R * 1e-6)
    lower .-= padding; upper .+= padding
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

# (3) rebuild with the derived shedding.
rotor = pnl.RigidWakeBody{kernel}(rotor.nodes, rotor.cells, [shedding1, shedding2];
    kerneloffset=kerneloffset_panel, kerneloffset_panel, kerneloffset_targets, kernelcutoff,
    semiinfinite_wake=false, watertight=true, ensure_winding=true, DBC)

# --------------------------------------------------------- wake + solver setup --
# Strictly finite, scalar-potential-only panel wake (no trailing filament).
wake_rotor = pnl.PanelWake(rotor;
    nwakerows=wake_rows, core_size=wake_core_size, include_final_filament=false,
    shed_with_induced_velocity=true)

solver_rotor = pnl.Backslash(rotor)
backend = pnl.FastMultipoleBackend(;              # body-only solve/system
    expansion_order=parse(Int,     get(ENV, "FMM_EXPANSION_ORDER", "8")),
    multipole_acceptance=parse(Float64, get(ENV, "FMM_ACCEPTANCE", "0.4")),
    leaf_size=parse(Int,           get(ENV, "FMM_LEAF_SIZE", "20")))
backend_wake = pnl.FastMultipoleBackend(;         # wake evaluation
    expansion_order=parse(Int,     get(ENV, "FMM_WAKE_EXPANSION_ORDER", "4")),
    multipole_acceptance=parse(Float64, get(ENV, "FMM_WAKE_ACCEPTANCE", "0.4")),
    leaf_size=parse(Int,           get(ENV, "FMM_WAKE_LEAF_SIZE", "50")))
kj_backend = pnl.FastMultipoleBackend(;
    expansion_order=3, multipole_acceptance=0.4, leaf_size=1000)

formulation = if formulation_name == "direct"
    pnl.DirectWakePotential(; recompute_interval)
elseif formulation_name == "velocity"
    pnl.VelocityThroughSources()
elseif formulation_name == "green"
    pnl.GreenReconstruction()
else
    error("Unknown FORMULATION=$(repr(formulation_name)); use direct, velocity, or green")
end

# ---------------------------------------------------------- freestream schedule --
# Constant RPM; single monotone freestream decrease START_VINF -> TERMINAL_VINF
# over FREESTREAM_DECREASE_REVS, then held.
t_decrease = freestream_decrease_revs * sec_per_rev
function magVinf(t)
    if t_decrease <= 0
        return terminal_vinf
    elseif t >= t_decrease
        return terminal_vinf
    else
        s = t / t_decrease                                   # linear 0->1
        return start_vinf + (terminal_vinf - start_vinf) * s
    end
end
Uinf(t) = magVinf(t) * Vinf_direction

# ------------------------------------------------------------------- frames -----
ω_full = 2 * pi * RPM / 60
frames = pnl.ReferenceFrame(rotor;
    origin=SVector{3}(0.0, 0.0, 0.0),
    v=SVector{3}(0.0, 0.0, 0.0),
    ω_axis=SVector{3}(-1.0, 0.0, 0.0),
    ω=ω_full,
    R=SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    name="vehicle",
    child_index=Int[],
    dependent_index=[1])

pnl.initialize_Das!((rotor,), frames, Uinf, t_range[1], t_range[2] - t_range[1];
    set_Das_eta_kinematic=0.2,
    set_Das_min_kinematic_displacement=0.01 * R)

maneuver!(frames, systems, wakes, t) = nothing   # constant RPM, no spin-up

systems      = (rotor,)
wakes        = (wake_rotor,)
body_solvers = (solver_rotor,)

# ------------------------------------------------------------------ monitors ----
norm_rotor = pnl.RotorNormalization(rho, 2 * R, 1)

pressure_bernoulli = pnl.PressureBernoulli(rho;
    unsteady=true, allow_partial=false, backend=backend)
force_monitor_bernoulli = pnl.ForceMonitor(length(t_range), 1;
    i_frame=1, normalization=norm_rotor, verbose=true)

# spanwise pressure-force loading on the +y blade (thrust dir -x)
one_blade = cp -> cp[2] > 0
span_bernoulli = pnl.SpanwiseLoadingMonitor(nbins, 1;
    i_frame=1, span_axis=[0.0, 1.0, 0.0],
    components=(thrust=[-1.0, 0.0, 0.0], tangential=[0.0, 0.0, 1.0]),
    per_length=true, select=one_blade, verbose=false)

# independent total-load cross-check
kj_monitor = pnl.KuttaJoukowskiForce(rotor, length(t_range), 1;
    rho, backend=kj_backend, i_frame=1, normalization=norm_rotor, verbose=true)

bound_circulation = pnl.BoundCirculationMonitor(rotor, length(t_range), 1;
    i_frame=1, radial_dimension, R)

# spanwise monitor must immediately follow the force monitor whose :F it reads.
monitors = (
    pressure_bernoulli, force_monitor_bernoulli, span_bernoulli,
    kj_monitor,
    bound_circulation,
)

# --------------------------------------------------------------- time march -----
if !isdir(save_path); mkpath(save_path); end
vtk_path = save_vtk ? save_path : nothing

println("\nBegin time march ($(length(t_range)) steps, formulation=$(formulation_name))...")
elapsed = @elapsed pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
    set_Das_eta_kinematic=NaN,
    set_Das_min_kinematic_displacement=0.01 * R,
    monitors,
    body_solvers, backend, backend_wake,
    formulation,
    verbose=true,
    path=vtk_path, name=base_name)
println("Time march wall time: $(round(elapsed, digits=1)) s")

# ------------------------------------------------------------- validation -------
finite_checks = Dict{String,Bool}(
    "body_strength"  => all(isfinite, rotor.strength),
    "body_velocity"  => all(isfinite, rotor.velocity),
    "wake_strength"  => all(all(isfinite, s) for s in wake_rotor.strength),
    "wake_nodes"     => all(all(isfinite, n) for n in wake_rotor.nodes),
    "CT_bernoulli"   => all(isfinite, force_monitor_bernoulli.force[axial_dimension, :]),
    "CT_kj"          => all(isfinite, kj_monitor.force[axial_dimension, :]),
)
all_finite = all(values(finite_checks))
for (k, v) in finite_checks
    v || println("  NON-FINITE: $k")
end

# CT sign convention matches the pressure-comparison driver (thrust = −axial).
CT_bernoulli = -force_monitor_bernoulli.force[axial_dimension, :]
CT_kj        = -kj_monitor.force[axial_dimension, :]

# final complete revolution statistics
k_final_start = max(1, length(t_range) - nt + 1)
tail_b = filter(isfinite, CT_bernoulli[k_final_start:end])
mean(v) = isempty(v) ? NaN : sum(v) / length(v)
ptp(v)  = isempty(v) ? NaN : maximum(v) - minimum(v)
# relative linear drift over the final rev: slope*window / mean
function linear_drift(v)
    n = length(v); n < 2 && return NaN
    xs = collect(0:n-1); xbar = sum(xs)/n; ybar = sum(v)/n
    num = sum((xs .- xbar) .* (v .- ybar)); den = sum((xs .- xbar).^2)
    slope = den == 0 ? 0.0 : num/den
    return abs(slope * (n-1)) / max(abs(ybar), eps())
end
CT_final_mean = mean(tail_b)
CT_final_ptp_rel   = ptp(tail_b) / max(abs(CT_final_mean), eps())
CT_final_drift_rel = linear_drift(tail_b)
CT_kj_final_mean = mean(filter(isfinite, CT_kj[k_final_start:end]))

residual_vinf = magVinf(t_range[end])
println("""
Final-revolution statistics (steps $(k_final_start):$(length(t_range))):
  residual freestream magVinf = $(round(residual_vinf, sigdigits=5)) m/s  (J = $(round(residual_vinf/(RPM/60*2R), sigdigits=4)))
  CT (Bernoulli) mean = $(round(CT_final_mean, sigdigits=6))
  CT (Bernoulli) rel peak-to-peak = $(round(CT_final_ptp_rel, sigdigits=4))  (gate <=0.05)
  CT (Bernoulli) rel linear drift = $(round(CT_final_drift_rel, sigdigits=4))  (gate <=0.025)
  CT (Kutta-Joukowski) mean = $(round(CT_kj_final_mean, sigdigits=6))
""")

stable = all_finite &&
    isfinite(CT_final_ptp_rel) && CT_final_ptp_rel <= 0.05 &&
    isfinite(CT_final_drift_rel) && CT_final_drift_rel <= 0.025

# --------------------------------------------------------------- outputs --------
# CT vs revolution CSV
csv_path = joinpath(save_path, "$(base_name)_CT_vs_rev.csv")
open(csv_path, "w") do io
    println(io, "step,revolution,CT_bernoulli,CT_kj")
    for k in 1:length(t_range)
        rev = (k - 1) * dt * RPM / 60
        println(io, "$k,$rev,$(CT_bernoulli[k]),$(CT_kj[k])")
    end
end
println("Wrote CT vs revolution CSV: $csv_path")

# metadata TOML
meta_path = joinpath(save_path, "$(base_name)_study_metadata.toml")
open(meta_path, "w") do io
    println(io, "run_name = \"$(run_name)\"")
    println(io, "formulation = \"$(formulation_name)\"")
    println(io, "recompute_interval = $(recompute_interval)")
    println(io, "mesh = \"$(rhpc_mesh)\"")
    println(io, "RPM = $(RPM)")
    println(io, "NT = $(nt)")
    println(io, "wake_rows_cap = $(wake_rows)")
    println(io, "active_wake_rows = $(wake_rotor.nwakes[])")
    println(io, "NREVS = $(nrevs)")
    println(io, "fill_revs = $(fill_revs)")
    println(io, "settle_revs = $(settle_revs)")
    println(io, "start_vinf = $(start_vinf)")
    println(io, "terminal_vinf = $(terminal_vinf)")
    println(io, "freestream_decrease_revs = $(freestream_decrease_revs)")
    println(io, "residual_vinf = $(residual_vinf)")
    println(io, "n_steps = $(n_steps)")
    println(io, "wall_time_s = $(elapsed)")
    println(io, "backend_body_order = $(backend.expansion_order)")
    println(io, "backend_wake_order = $(backend_wake.expansion_order)")
    println(io, "CT_bernoulli_final_mean = $(CT_final_mean)")
    println(io, "CT_bernoulli_final_ptp_rel = $(CT_final_ptp_rel)")
    println(io, "CT_bernoulli_final_drift_rel = $(CT_final_drift_rel)")
    println(io, "CT_kj_final_mean = $(CT_kj_final_mean)")
    println(io, "all_finite = $(all_finite)")
    println(io, "stable = $(stable)")
end
println("Wrote study metadata: $meta_path")

# completion marker — written ONLY after all validation and reporting succeed
if all_finite
    marker_path = joinpath(save_path, "COMPLETED")
    open(marker_path, "w") do io
        println(io, "run_name=$(run_name)")
        println(io, "stable=$(stable)")
        println(io, "CT_bernoulli_final_mean=$(CT_final_mean)")
        println(io, "wall_time_s=$(elapsed)")
    end
    println("Wrote completion marker: $marker_path  (stable=$(stable))")
else
    println("NON-FINITE fields detected; completion marker NOT written.")
end
