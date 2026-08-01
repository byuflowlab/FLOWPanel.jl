#=##############################################################################
# Task 3a — Sphere added-mass gate (code_audit/plan.md)
#
# NonLiftingBody{ConstantSource} sphere, prescribed smooth-ramp translation
# U(t) x̂ via frames/maneuver, unsteady PressureBernoulli + ForceMonitor,
# no wake/shedding. Isolates the φ probe, exterior-trace handling, BE/BDF2
# ladder, ALE w·∇φ term, and inertial-trace kinetic term.
#
# Exact solution (sphere of radius R translating at U(t) in still fluid,
# θ measured from the direction of motion):
#   φ_ext = −½ U(t) R³ cosθ / r²
#   p − p∞ |surface = (ρR/2) U̇ cosθ + (ρU²/8)(9cos²θ − 5)
#   F = −½ ρ Vol U̇ x̂,  Vol = (4/3)πR³  (added mass of a sphere = ½ρVol)
#
# Gates:
#   G1: per-panel P within a few % of the analytic surface pressure away from
#       the mesh poles at the peak-acceleration step.
#   G2: force error < ~2–5% once BDF2 is active (step ≥ 2) during the ramp.
#   G3: after the ramp (U̇=0), |Fx| decays to ≈0 relative to the peak force.
#
# Run:  julia --project code_audit/scripts/task3a_sphere_added_mass.jl
# Outputs: code_audit/results/task3a/
=###############################################################################

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
import Statistics: mean, median
import Printf: @printf, @sprintf

const RESULTS_DIR = joinpath(@__DIR__, "..", "results", "task3a")
mkpath(RESULTS_DIR)

#------- parameters -------#

R        = 1.0                  # (m) sphere radius
rho      = 1.225                # (kg/m^3)
Umax     = 1.0                  # (m/s) terminal translation speed
T_ramp   = 1.0                  # (s) ramp duration
dt       = parse(Float64, get(ENV, "TASK3A_DT", "0.02"))    # (s) time step
t_end    = 1.4 * T_ramp
n_theta  = parse(Int, get(ENV, "TASK3A_NTHETA", "16"))      # polar panel count
n_phi    = parse(Int, get(ENV, "TASK3A_NPHI", "32"))        # azimuthal panel count

U_of(t)    = t <= T_ramp ? Umax * 0.5 * (1 - cos(pi * t / T_ramp)) : Umax
dUdt_of(t) = t <= T_ramp ? Umax * 0.5 * (pi / T_ramp) * sin(pi * t / T_ramp) : 0.0

t_range = range(0.0, stop=t_end, step=dt)
nt = length(t_range)

#------- geometry: UV sphere, outward winding -------#

function sphere_mesh(R, n_theta, n_phi)
    np = n_phi
    nnodes = 2 + (n_theta - 1) * np
    nodes = zeros(3, nnodes)
    nodes[:, 1] .= (0.0, 0.0, R)     # north pole (+z)
    nodes[:, 2] .= (0.0, 0.0, -R)    # south pole (−z)
    node_id(j, k) = 2 + (j - 1) * np + mod(k - 1, np) + 1
    for j in 1:n_theta-1
        θ = j * pi / n_theta
        for k in 1:np
            ϕ = 2pi * (k - 1) / np
            nodes[:, node_id(j, k)] .= (R * sin(θ) * cos(ϕ), R * sin(θ) * sin(ϕ), R * cos(θ))
        end
    end

    ncells = 2 * np + 2 * np * (n_theta - 2)
    cells = zeros(Int, 3, ncells)
    c = 0
    # north cap
    for k in 1:np
        c += 1
        cells[:, c] .= (1, node_id(1, k), node_id(1, k + 1))
    end
    # bands
    for j in 1:n_theta-2
        for k in 1:np
            uk, ukp = node_id(j, k), node_id(j, k + 1)
            lk, lkp = node_id(j + 1, k), node_id(j + 1, k + 1)
            c += 1
            cells[:, c] .= (uk, lk, lkp)
            c += 1
            cells[:, c] .= (uk, lkp, ukp)
        end
    end
    # south cap
    for k in 1:np
        c += 1
        cells[:, c] .= (2, node_id(n_theta - 1, k + 1), node_id(n_theta - 1, k))
    end
    @assert c == ncells
    return nodes, cells
end

nodes, cells = sphere_mesh(R, n_theta, n_phi)

ElementTypes = Union{pnl.ConstantSource}
body = pnl.NonLiftingBody{ElementTypes}(nodes, cells)

pnl.calc_normals!(body)
pnl.calc_controlpoints!(body)

# verify outward normals (auto-flip if the constructor re-wound them inward)
let
    center = vec(mean(body.nodes, dims=2))
    n_out = count(1:body.ncells) do i
        n = SVector{3}(body.normals[1,i], body.normals[2,i], body.normals[3,i])
        cp = SVector{3}(body.controlpoints[1,i], body.controlpoints[2,i], body.controlpoints[3,i])
        pnl.dot(n, cp - center) > 0
    end
    if n_out < body.ncells / 2
        global body = pnl.NonLiftingBody{ElementTypes}(nodes, cells; flip_normals=true)
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        n_out2 = count(1:body.ncells) do i
            n = SVector{3}(body.normals[1,i], body.normals[2,i], body.normals[3,i])
            cp = SVector{3}(body.controlpoints[1,i], body.controlpoints[2,i], body.controlpoints[3,i])
            pnl.dot(n, cp - center) > 0
        end
        println("flipped normals: outward $(n_out2)/$(body.ncells)")
        @assert n_out2 == body.ncells "normals not consistently outward after flip"
    else
        @assert n_out == body.ncells "normals not consistently outward ($(n_out)/$(body.ncells))"
        println("normals outward: $(n_out)/$(body.ncells)")
    end
end

# discrete volume via divergence theorem (for added-mass discretization ratio)
areas = pnl.calc_areas(body)
V_mesh = sum(1:body.ncells) do i
    n = SVector{3}(body.normals[1,i], body.normals[2,i], body.normals[3,i])
    cp = SVector{3}(body.controlpoints[1,i], body.controlpoints[2,i], body.controlpoints[3,i])
    pnl.dot(cp, n) * areas[i] / 3
end
V_exact = 4/3 * pi * R^3
println("panels = $(body.ncells);  V_mesh/V_exact = $(V_mesh / V_exact)")

#------- frames / maneuver -------#

frames = pnl.ReferenceFrame(body;
    origin=SVector{3}(0.0, 0.0, 0.0),
    v=SVector{3}(0.0, 0.0, 0.0),
    name="sphere",
)

function sphere_maneuver!(frames, systems, wakes, t)
    frame = frames[1]
    frames[1] = typeof(frame)(
        frame.x,
        SVector{3}(U_of(t), 0.0, 0.0),
        frame.ω_axis,
        0.0,
        frame.R,
        frame.Rp2g,
        frame.name,
        frame.parent_index,
        frame.child_index,
        frame.dependent_index,
    )
    return nothing
end

Uinf(t) = SVector{3}(0.0, 0.0, 0.0)

#------- monitors -------#

backend = pnl.FastMultipoleBackend(expansion_order=8, multipole_acceptance=0.4, leaf_size=40)

pressure_monitor = pnl.PressureBernoulli(rho; unsteady=true, backend=backend, file=false)
force_monitor = pnl.ForceMonitor(nt, 1;
    normalization=pnl.NoNormalization(),
    i_frame=-1,
    correct_kuttacondition=false,
    verbose=false,
    file=false,
)

# recorder: per-step copies of pressure, phi_dot, geometry, and frame velocity
record = (;
    t=Float64[],
    U=Float64[],
    P=Vector{Float64}[],
    phi_dot=Vector{Float64}[],
    cps=Matrix{Float64}[],
    center=SVector{3,Float64}[],
)
function recorder(systems, wakes, frames, uinf, i_step, dt)
    b = systems isa Tuple ? systems[1] : systems
    push!(record.t, i_step * dt)
    push!(record.U, frames[1].v[1])
    push!(record.P, copy(pressure_monitor.pressure[1]))
    push!(record.phi_dot, copy(pressure_monitor.phi_dot[1]))
    push!(record.cps, copy(b.controlpoints))
    push!(record.center, SVector{3}(vec(mean(b.nodes, dims=2))...))
    return nothing
end

monitors = (pressure_monitor, force_monitor, recorder)

#------- run -------#

solver = pnl.Backslash(body)

elapsed = @elapsed pnl.simulate!((body,), (nothing,), frames, sphere_maneuver!, Uinf, t_range;
    body_solvers=(solver,),
    backend,
    monitors,
    path=nothing,
    verbose=false,
)
println("simulate! finished in $(round(elapsed, digits=1)) s ($(nt) steps)")

#------- analysis -------#

# identify mesh-pole cap panels (skinny triangles touching the ±z poles)
cap_panel = [any(cells[k, i] in (1, 2) for k in 1:3) for i in 1:body.ncells]

P_exact_of(cosθ, U, dUdt) = rho * R / 2 * dUdt * cosθ + rho * U^2 / 8 * (9 * cosθ^2 - 5)

# --- Gate 1: per-panel pressure at the peak-acceleration step ---
i_peak = argmin(abs.(collect(t_range) .- T_ramp / 2))   # 1-based index into records
t_peak = record.t[i_peak]
U_pk, A_pk = U_of(t_peak), dUdt_of(t_peak)

cps = record.cps[i_peak]
ctr = record.center[i_peak]
cosθs = [ (cps[1, p] - ctr[1]) / pnl.norm(SVector{3}(cps[1,p]-ctr[1], cps[2,p]-ctr[2], cps[3,p]-ctr[3])) for p in 1:body.ncells ]
P_num = record.P[i_peak]
P_ex  = [P_exact_of(cosθs[p], U_pk, A_pk) for p in 1:body.ncells]
P_scale = maximum(abs, P_ex)

rel_err = abs.(P_num .- P_ex) ./ P_scale
sel = .!cap_panel
g1_median = median(rel_err[sel])
g1_p90 = sort(rel_err[sel])[ceil(Int, 0.90 * count(sel))]
g1_max = maximum(rel_err[sel])

println("\n--- Gate 1: per-panel pressure @ t=$(round(t_peak, digits=3)) (peak U̇) ---")
@printf("  U = %.4f m/s, dU/dt = %.4f m/s^2, P_scale = %.4f Pa\n", U_pk, A_pk, P_scale)
@printf("  |P_num - P_exact|/max|P_exact| (non-cap panels): median = %.4f, p90 = %.4f, max = %.4f\n",
        g1_median, g1_p90, g1_max)
@printf("  cap panels only: median = %.4f, max = %.4f\n",
        median(rel_err[cap_panel]), maximum(rel_err[cap_panel]))

# --- Gate 2: force during the ramp, steps >= 2 (BDF2 active) ---
Fx = force_monitor.force[1, :]
Fy = force_monitor.force[2, :]
Fz = force_monitor.force[3, :]
F_exact = [-0.5 * rho * V_exact * dUdt_of(t) for t in t_range]
F_peak = maximum(abs, F_exact)

gate2_idx = [i for (i, t) in enumerate(t_range)
             if i >= 3 && t <= T_ramp && abs(dUdt_of(t)) >= 0.2 * maximum(abs, dUdt_of.(t_range))]
g2_rel = [abs(Fx[i] - F_exact[i]) / abs(F_exact[i]) for i in gate2_idx]
g2_max = maximum(g2_rel)
g2_median = median(g2_rel)

println("\n--- Gate 2: force vs -0.5*rho*Vol*dU/dt (steps >= 2, ramp, |U̇| >= 20% peak) ---")
@printf("  F_peak(exact) = %.4f N;  V_mesh/V_exact = %.4f\n", F_peak, V_mesh / V_exact)
@printf("  rel force error: median = %.4f, max = %.4f over %d steps\n", g2_median, g2_max, length(gate2_idx))
@printf("  lateral force: max|Fy| = %.2e, max|Fz| = %.2e (vs F_peak %.4f)\n",
        maximum(abs, Fy), maximum(abs, Fz), F_peak)

# --- Gate 3: post-ramp force decay ---
post_idx = [i for (i, t) in enumerate(t_range) if t > T_ramp + 2 * dt]
g3 = maximum(abs, Fx[post_idx]) / F_peak
println("\n--- Gate 3: post-ramp |Fx|/F_peak = $(round(g3, digits=5)) ---")

#------- outputs -------#

open(joinpath(RESULTS_DIR, "task3a_force_history.csv"), "w") do io
    println(io, "t,U,dUdt,Fx_num,Fy_num,Fz_num,Fx_exact,rel_err")
    for (i, t) in enumerate(t_range)
        re = abs(F_exact[i]) > 1e-12 ? abs(Fx[i] - F_exact[i]) / abs(F_exact[i]) : NaN
        @printf(io, "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                t, U_of(t), dUdt_of(t), Fx[i], Fy[i], Fz[i], F_exact[i], re)
    end
end

open(joinpath(RESULTS_DIR, "task3a_panel_pressure_peak.csv"), "w") do io
    println(io, "cos_theta,P_num,P_exact,rel_err,is_cap")
    for p in 1:body.ncells
        @printf(io, "%.6e,%.6e,%.6e,%.6e,%d\n",
                cosθs[p], P_num[p], P_ex[p], rel_err[p], cap_panel[p] ? 1 : 0)
    end
end

#------- verdict -------#

pass1 = g1_p90 < 0.05
pass2 = g2_max < 0.05
pass3 = g3 < 0.02
println("\n================ TASK 3a VERDICT ================")
println("  Gate 1 (per-panel P, p90 < 5%):        $(pass1 ? "PASS" : "FAIL") (p90 = $(round(g1_p90, digits=4)))")
println("  Gate 2 (force err < 5%, BDF2 active):  $(pass2 ? "PASS" : "FAIL") (max = $(round(g2_max, digits=4)))")
println("  Gate 3 (post-ramp force ~ 0):          $(pass3 ? "PASS" : "FAIL") (ratio = $(round(g3, digits=5)))")
println("  OVERALL: $(pass1 && pass2 && pass3 ? "PASS" : "FAIL")")
println("=================================================")
