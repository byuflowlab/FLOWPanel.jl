# Sphere added-mass gate (code_audit Task 3a).
#
# A closed NonLiftingBody{ConstantSource} sphere undergoes a prescribed
# smooth-ramp translation U(t) x̂ through frames/maneuver with no wake and no
# shedding, monitored by unsteady PressureBernoulli + ForceMonitor. This
# isolates the unsteady-Bernoulli phi-dot machinery: the scalar-potential
# probe, the exterior-trace handling, the zero/BE/variable-step-BDF2 startup
# ladder, the ALE w·∇φ term, and the inertial-trace kinetic term.
#
# Exact solution (sphere radius R translating at U(t) in still fluid, θ from
# the direction of motion):
#   surface pressure  p − p∞ = (ρR/2) U̇ cosθ + (ρU²/8)(9cos²θ − 5)
#   total force       F = −½ ρ Vol U̇ x̂  (added mass of a sphere = ½ρVol)
#
# Tolerances were calibrated in
# code_audit/scripts/task3a_sphere_added_mass.jl (2026-07-17). The discrete
# added mass overshoots the continuum by a mesh-dependent constant that
# converges under refinement (force error 4.9% at 528 panels, 3.7% at 960,
# 2.7% at 2208), so the raw force gate is a coarse-mesh bound (7%) and the
# tight gate checks that the ratio Fx/F_exact is constant in time (spread
# ~0.4% measured): any BE/BDF2-ladder, ALE, or trace-sign defect shows up as
# a time-varying ratio, while a flat ratio isolates the residual to spatial
# discretization. Non-cap per-panel pressure p90 measured ~5.0% at 528 panels.

using Test
import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
import Statistics: mean, median

function _added_mass_sphere_mesh(R::Real, n_theta::Int, n_phi::Int)
    np = n_phi
    nodes = zeros(3, 2 + (n_theta - 1) * np)
    nodes[:, 1] .= (0.0, 0.0, R)     # north pole
    nodes[:, 2] .= (0.0, 0.0, -R)    # south pole
    node_id(j, k) = 2 + (j - 1) * np + mod(k - 1, np) + 1
    for j in 1:n_theta-1
        θ = j * pi / n_theta
        for k in 1:np
            ϕ = 2pi * (k - 1) / np
            nodes[:, node_id(j, k)] .= (R*sin(θ)*cos(ϕ), R*sin(θ)*sin(ϕ), R*cos(θ))
        end
    end

    ncells = 2 * np + 2 * np * (n_theta - 2)
    cells = zeros(Int, 3, ncells)
    c = 0
    for k in 1:np
        cells[:, c += 1] .= (1, node_id(1, k), node_id(1, k + 1))
    end
    for j in 1:n_theta-2, k in 1:np
        uk, ukp = node_id(j, k), node_id(j, k + 1)
        lk, lkp = node_id(j + 1, k), node_id(j + 1, k + 1)
        cells[:, c += 1] .= (uk, lk, lkp)
        cells[:, c += 1] .= (uk, lkp, ukp)
    end
    for k in 1:np
        cells[:, c += 1] .= (2, node_id(n_theta - 1, k + 1), node_id(n_theta - 1, k))
    end
    return nodes, cells
end

@testset verbose=true "Sphere added-mass gate (unsteady Bernoulli)" begin
    R, rho, Umax, T_ramp = 1.0, 1.225, 1.0, 1.0
    dt, t_end = 0.04, 1.2
    n_theta, n_phi = 12, 24

    U_of(t)    = t <= T_ramp ? Umax * 0.5 * (1 - cos(pi * t / T_ramp)) : Umax
    dUdt_of(t) = t <= T_ramp ? Umax * 0.5 * (pi / T_ramp) * sin(pi * t / T_ramp) : 0.0
    t_range = range(0.0, stop=t_end, step=dt)
    nt = length(t_range)

    nodes, cells = _added_mass_sphere_mesh(R, n_theta, n_phi)
    body = pnl.NonLiftingBody{Union{pnl.ConstantSource}}(nodes, cells)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)

    # normals must be consistently outward for the exterior problem
    center0 = vec(mean(body.nodes, dims=2))
    n_outward = count(1:body.ncells) do i
        n = SVector{3}(body.normals[1,i], body.normals[2,i], body.normals[3,i])
        cp = SVector{3}(body.controlpoints[1,i], body.controlpoints[2,i], body.controlpoints[3,i])
        pnl.dot(n, cp - center0) > 0
    end
    @test n_outward == body.ncells

    areas = pnl.calc_areas(body)
    V_mesh = sum(1:body.ncells) do i
        n = SVector{3}(body.normals[1,i], body.normals[2,i], body.normals[3,i])
        cp = SVector{3}(body.controlpoints[1,i], body.controlpoints[2,i], body.controlpoints[3,i])
        pnl.dot(cp, n) * areas[i] / 3
    end
    V_exact = 4/3 * pi * R^3
    @test isapprox(V_mesh, V_exact; rtol=0.05)

    frames = pnl.ReferenceFrame(body;
        origin=SVector{3}(0.0, 0.0, 0.0),
        v=SVector{3}(0.0, 0.0, 0.0),
        name="sphere")

    function sphere_maneuver!(frames, systems, wakes, t)
        frame = frames[1]
        frames[1] = typeof(frame)(frame.x, SVector{3}(U_of(t), 0.0, 0.0),
            frame.ω_axis, 0.0, frame.R, frame.Rp2g, frame.name,
            frame.parent_index, frame.child_index, frame.dependent_index)
        return nothing
    end
    Uinf(t) = SVector{3}(0.0, 0.0, 0.0)

    backend = pnl.FastMultipoleBackend(expansion_order=8, multipole_acceptance=0.4, leaf_size=40)
    pressure_monitor = pnl.PressureBernoulli(rho; unsteady=true, backend=backend, file=false)
    force_monitor = pnl.ForceMonitor(nt, 1;
        normalization=pnl.NoNormalization(), i_frame=-1,
        correct_kuttacondition=false, file=false)

    # record the per-panel pressure and geometry at the peak-acceleration step
    i_peak = argmin(abs.(collect(t_range) .- T_ramp / 2)) - 1   # i_step of peak U̇
    peak = (; P=Float64[], cps=zeros(3, 0), center=zeros(3))
    function peak_recorder(systems, wakes, frames, uinf, i_step, dt)
        i_step == i_peak || return nothing
        b = systems isa Tuple ? systems[1] : systems
        peak = (; P=copy(pressure_monitor.pressure[1]),
                  cps=copy(b.controlpoints),
                  center=vec(mean(b.nodes, dims=2)))
        return nothing
    end

    solver = pnl.Backslash(body)
    pnl.simulate!((body,), (nothing,), frames, sphere_maneuver!, Uinf, t_range;
        body_solvers=(solver,),
        backend,
        monitors=(pressure_monitor, force_monitor, peak_recorder),
        path=nothing,
        verbose=false)

    #--- Gate 1: per-panel surface pressure at peak acceleration ---#
    @test length(peak.P) == body.ncells
    t_peak = t_range[i_peak + 1]
    U_pk, A_pk = U_of(t_peak), dUdt_of(t_peak)
    cap_panel = [any(cells[k, i] in (1, 2) for k in 1:3) for i in 1:body.ncells]

    P_ex = zeros(body.ncells)
    rel_err = zeros(body.ncells)
    for p in 1:body.ncells
        d = SVector{3}(peak.cps[1,p] - peak.center[1],
                       peak.cps[2,p] - peak.center[2],
                       peak.cps[3,p] - peak.center[3])
        cosθ = d[1] / pnl.norm(d)
        P_ex[p] = rho * R / 2 * A_pk * cosθ + rho * U_pk^2 / 8 * (9 * cosθ^2 - 5)
    end
    P_scale = maximum(abs, P_ex)
    rel_err .= abs.(peak.P .- P_ex) ./ P_scale
    noncap = rel_err[.!cap_panel]
    @test median(noncap) < 0.04
    @test sort(noncap)[ceil(Int, 0.90 * length(noncap))] < 0.07

    #--- Gate 2: total force vs added-mass reaction (BDF2 active, ramp) ---#
    Fx = force_monitor.force[1, :]
    A_max = maximum(abs, dUdt_of.(t_range))
    gate_idx = [i for (i, t) in enumerate(t_range)
                if i >= 3 && t <= T_ramp && abs(dUdt_of(t)) >= 0.2 * A_max]
    @test !isempty(gate_idx)
    ratios = [Fx[i] / (-0.5 * rho * V_exact * dUdt_of(t_range[i])) for i in gate_idx]
    err_cont = maximum(r -> abs(r - 1), ratios)
    ratio_spread = (maximum(ratios) - minimum(ratios)) / mean(ratios)
    @test err_cont < 0.07        # continuum added mass, coarse-mesh bound
    @test ratio_spread < 0.02    # flat-in-time ratio isolates the phi-dot ladder

    # lateral force components stay at machine-precision levels of the peak
    F_peak = 0.5 * rho * V_exact * A_max
    @test maximum(abs, force_monitor.force[2, :]) < 1e-8 * F_peak
    @test maximum(abs, force_monitor.force[3, :]) < 1e-8 * F_peak

    #--- Gate 3: post-ramp (U̇ = 0) force decays to zero ---#
    post_idx = [i for (i, t) in enumerate(t_range) if t > T_ramp + 2 * dt]
    @test !isempty(post_idx)
    @test maximum(abs, Fx[post_idx]) < 0.02 * F_peak
end
