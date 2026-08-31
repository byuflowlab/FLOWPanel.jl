## Plate, rotor, with no ground effect
# Author: Timothy Harlow
# Created: December 5, 2025

import FLOWPanel as pnl
using Printf
include(joinpath(pnl.examples_path, "helper_functions.jl"))
using FLOWPanel.FastMultipole.StaticArrays
using VSPGeom

run_name = "rotor_hover_mudiag"
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
t_range = range(0.0, step=dt, length=n_steps)[1:10]

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
# stl_file   = joinpath(read_path, "phantom_3_mod3_rev5.stl")

# phantom_3_rebuild_r2.msh
msh_file  = joinpath(read_path, "phantom_3_rebuild_r2.msh")
te_indices_1 = [9, 175, 127]
te_indices_2 = [13, 286, 238]

# # phantom_3_rebuild_r3.msh
# msh_file  = joinpath(read_path, "phantom_3_rebuild_r3.msh")
# te_indices_1 = [8, 523, 223] .+ 1
# te_indices_2 = [12, 997, 697] .+ 1

# # phantom_3_rebuild_r4.msh
# msh_file  = joinpath(read_path, "phantom_3_rebuild_r4.msh")
# te_indices_1 = [7, 952, 4] .+ 1
# te_indices_2 = [3, 478, 0] .+ 1

# STL file
# mesh = VSPGeom.readSTL(stl_file)[1]
# scale = 1/1000 # convert to meters
# radius = 119.38 * scale
# for point in mesh.points
#     point .*= scale
# end

# MSH file
msh = pnl.read_gmsh(msh_file)
nodes, cells = pnl.meshes2nodes_cells(msh)

# scale to proper radius
nodes .*= R / maximum(nodes[1, :])

# place-holder shedding
shedding = pnl.noshedding

# --- Construct RigidWakeBody ---
kernel = Union{pnl.ConstantSource, pnl.VortexRing}
# kernel = pnl.VortexRing
DBC = kernel == pnl.VortexRing ? false : true
rotor = pnl.RigidWakeBody{kernel}(nodes, cells, shedding;
            core_size,
            kernelcutoff,
            semiinfinite_wake=false,
            watertight=true,
            DBC)

pnl.write_vtk(joinpath(save_path, "rotor_hover_check"), rotor)

# update shedding
bbox = (pnl.SVector{3}(-R*1.2, -1.0, -1.0), pnl.SVector{3}(-R*0.1, 1.0, 1.0))
bbox = nothing
shedding1 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_1[1], te_indices_1[2]; bbox, end_node=te_indices_1[3], normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)
bbox = (pnl.SVector{3}(R*0.1, -1.0, -1.0), pnl.SVector{3}(R*1.2, 1.0, 1.0))
bbox = nothing
shedding2 = pnl.calc_shedding_from_seed(rotor.nodes, rotor.cells, te_indices_2[1], te_indices_2[2]; bbox, end_node=te_indices_2[3], normal_jump_tol=0.2, max_turn_angle=pi/3, debug=false)

rotor = pnl.RigidWakeBody{kernel}(rotor.nodes, rotor.cells, [shedding1, shedding2];
                        core_size,
                        kernelcutoff,
                        semiinfinite_wake=false,
                        watertight=true,
                        ensure_winding=true,
                        DBC)

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

# Reference frame
frames = pnl.ReferenceFrame(rotor;
    origin = SVector{3}(0.0, 0.0, 0.0),
    v = SVector{3}(0.0, 0.0, 0.0),
    ω_axis = SVector{3}(0.0, 1.0, 0.0),
    ω = 2*pi * RPM/60, # rad/s
    R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    name = "vehicle",
    child_index = Int[],
    dependent_index = [1]
)

pnl.initialize_Das!((rotor,), frames, Uinf, t_range[1], t_range[2] - t_range[1];
    set_Das_eta_kinematic=init_Das_eta_kinematic)

# solver_rotor = pnl.FGSSolver(rotor;
#             max_iterations=500,
#             tolerance=1.0e-6,
#             rlx=1.0,
#             expansion_order=10,
#             multipole_acceptance=0.4,
#             leaf_size=150,
#             shrink=true,
#             recenter=false,
#             inner_iterations=20,
#             reverse_pass=false,
#             verbose=false
#         )
solver_rotor = pnl.Backslash(rotor)

backend = pnl.FastMultipoleBackend()

# Maneuver
maneuver!(frames, systems, wakes, t) = nothing

## =========================================================
# RUN SIMULATION
# ==========================================================
systems      = (rotor,)
wakes        = (wake_rotor,)
body_solvers = (solver_rotor,)
mu_log = Tuple{Int,Float64,Float64}[]
function mu_monitor(systems, wakes, frames, uinf, i_step)
    body = systems[1]
    mu = view(body.strength, :, 2)
    push!(mu_log, (i_step, maximum(abs, mu), sum(abs, mu) / length(mu)))
end

# wake-induced velocity probe at body control points, split panel vs particle
wake_v_log = Tuple{Int,Float64,Float64,Float64,Float64,Float64,Float64,Float64}[]
const wake_v_probes = pnl.FastMultipole.ProbeSystem(rotor.ncells, Float64)
function _reset_probes_at_cps!(probes, body)
    zero_v = zero(SVector{3, Float64})
    zero_h = zero(SMatrix{3, 3, Float64, 9})
    @inbounds for k in 1:body.ncells
        probes.position[k] = SVector{3, Float64}(
            body.controlpoints[1, k],
            body.controlpoints[2, k],
            body.controlpoints[3, k],
        )
        probes.gradient[k] = zero_v
        probes.scalar_potential[k] = 0.0
        probes.hessian[k] = zero_h
    end
end
function _vmag_stats(probes, body)
    n_dot_v_sum = 0.0
    v_mag_sum = 0.0
    v_mag_max = 0.0
    @inbounds for k in 1:body.ncells
        V = probes.gradient[k]
        nrm = SVector{3, Float64}(body.normals[1, k], body.normals[2, k], body.normals[3, k])
        n_dot_v_sum += abs(V[1]*nrm[1] + V[2]*nrm[2] + V[3]*nrm[3])
        m = sqrt(V[1]^2 + V[2]^2 + V[3]^2)
        v_mag_sum += m
        v_mag_max = max(v_mag_max, m)
    end
    return n_dot_v_sum / body.ncells, v_mag_sum / body.ncells, v_mag_max
end
function wake_v_monitor(systems, wakes, frames, uinf, i_step)
    body = systems[1]
    wake_sources = pnl._collect_wake_sources(wakes)

    # split sources: PanelParticleWake's get_sources is (panel_wake_sources..., pfield)
    # so all but the last are panels, and the last is the particle field
    panel_sources = isempty(wake_sources) ? () : wake_sources[1:end-1]
    particle_source = isempty(wake_sources) ? () : (wake_sources[end],)

    # 1) panel-only contribution
    _reset_probes_at_cps!(wake_v_probes, body)
    if !isempty(panel_sources)
        pnl.influence!((wake_v_probes,), panel_sources, backend;
            precalc=false, scalar_potential=false, gradient=true, hessian=(false,))
    end
    p_ndv_mean, p_v_mean, p_v_max = _vmag_stats(wake_v_probes, body)

    # find the body CP with worst panel-induced |V|, and distance to nearest wake vertex
    if i_step > 0 && !isempty(panel_sources)
        worst_k = 0
        worst_v = 0.0
        @inbounds for k in 1:body.ncells
            V = wake_v_probes.gradient[k]
            m = sqrt(V[1]^2 + V[2]^2 + V[3]^2)
            if m > worst_v
                worst_v = m; worst_k = k
            end
        end
        cp = SVector{3,Float64}(body.controlpoints[1, worst_k], body.controlpoints[2, worst_k], body.controlpoints[3, worst_k])
        # nearest wake-panel vertex; PanelParticleWake.panel_wake.nodes is a Vector{Array{TF,3}} indexed by surface
        min_d = Inf
        pw = wakes[1].panel_wake
        for i_surf in eachindex(pw.nodes)
            nodes_surf = pw.nodes[i_surf]   # (3, nrows+1, nshed+1)
            for j in axes(nodes_surf, 2), k2 in axes(nodes_surf, 3)
                vx = nodes_surf[1, j, k2]; vy = nodes_surf[2, j, k2]; vz = nodes_surf[3, j, k2]
                d = sqrt((vx-cp[1])^2 + (vy-cp[2])^2 + (vz-cp[3])^2)
                if d < min_d; min_d = d; end
            end
        end
        @printf("    [step %d] worst body CP idx=%d at (%.4f,%.4f,%.4f), |V|=%.2f, nearest wake vertex at d=%.4f mm\n",
                i_step, worst_k, cp[1], cp[2], cp[3], worst_v, min_d*1000)
    end

    # 2) particle-only contribution
    _reset_probes_at_cps!(wake_v_probes, body)
    if !isempty(particle_source)
        pnl.influence!((wake_v_probes,), particle_source, backend;
            precalc=false, scalar_potential=false, gradient=true, hessian=(false,))
    end
    pa_ndv_mean, pa_v_mean, pa_v_max = _vmag_stats(wake_v_probes, body)

    # reference: ω·R_tip
    ωmag = abs(frames[1].ω)
    rmax = 0.0
    @inbounds for k in 1:body.ncells
        rmax = max(rmax, sqrt(body.controlpoints[1, k]^2 + body.controlpoints[3, k]^2))
    end
    omega_R_tip = ωmag * rmax

    push!(wake_v_log, (i_step, p_ndv_mean, p_v_mean, p_v_max, pa_ndv_mean, pa_v_mean, pa_v_max, omega_R_tip))
end

monitors = (pnl.PressureBernoulli(rho; unsteady=true,
                    allow_partial=true,
                    correct_kuttacondition=p_correct_kuttacondition_flag),
            pnl.ForceMonitor(length(t_range), 1; # un-normalized, global frame
                    i_frame=-1,
                    normalization=pnl.RotorNormalization(rho, 2*R, 1),
                    correct_kuttacondition=p_correct_kuttacondition_flag,
                    # normalization=pnl.NoNormalization(),
                    verbose=false
                ),
            mu_monitor,
            wake_v_monitor,
            # pnl.KuttaJoukowskiForce(rotor, length(t_range), 1;
            #         rho, backend,
            #         normalization=pnl.RotorNormalization(rho, 2*R, 1)
            #     )
            )

println("\nBegin rotor hover simulation ($(n_steps) steps)...")
name = "rotor_hover"
@time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
    # Das was initialized before constructing the matrixful solver.
    set_Das_eta_kinematic=NaN,
    # set_Das_eta_freestream=0.1,
    monitors,
    body_solvers, backend, verbose=true,
    path=save_path, name,
)

println("Thrust Coefficient: ", monitors[2].force[2,:])
println("\n|μ| evolution (step, max, mean):")
for (i, mx, mn) in mu_log
    println("  step $i: max=$(round(mx, sigdigits=5))  mean=$(round(mn, sigdigits=5))")
end
println("\nWake-induced velocity at body CPs (panels vs particles; ωR_tip ≈ $(round(wake_v_log[end][8], sigdigits=4))):")
println("  step | panel: n·V mean   |V| mean   |V| max  || particle: n·V mean   |V| mean   |V| max")
for (i, p_ndv, p_vm, p_vx, pa_ndv, pa_vm, pa_vx, _) in wake_v_log
    @printf("  %4d |        %7.3f   %7.3f   %7.3f  ||           %7.3f   %7.3f   %7.3f\n",
        i, p_ndv, p_vm, p_vx, pa_ndv, pa_vm, pa_vx)
end
