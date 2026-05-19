## Capped Dirichlet wing: compare PressureBernoulli, PressureLaplace, and KuttaJoukowskiForce
## This is an exploratory cross-check on a solved capped body.

import FLOWPanel as pnl
using LinearAlgebra: norm, dot

include(joinpath(@__DIR__, "simple_wing_capped_dirichlet.jl"))

function expected_negative_tangent_velocity(body)
    expected = similar(body.velocity)
    for p in 1:body.ncells
        n = body.normals[:, p]
        u = body.velocity[:, p]
        expected[:, p] .= -(u .- dot(u, n) .* n .+ body.velocity_kinematic[:, p])
    end
    return expected
end

function prepare_pressure_state!(body, Vinf; backend=pnl.DirectBackend(), velocity_gradient=true)
    body.needs_velocity_gradient[] = velocity_gradient
    body.velocity .= 0.0
    body.velocity_gradient .= 0.0
    pnl.calcfield_U!(body, Vinf; backend)
    if velocity_gradient
        pnl.influence!((body,), (body,), backend;
            scalar_potential=false,
            velocity=false,
            velocity_gradient=true)
    end
    return body
end

function compare_pressure_models(;
    meshfile=joinpath(pnl.examples_path, "data", "naca0012_nc101_nw26.msh"),
    AOA=5.88,
    magVinf=56.0,
    rho=1.225,
    backend=pnl.DirectBackend(),
    dt=1.0,
    laplace_gradient_mode=:surface_velocity,
)
    body_ref, Vinf = build_simple_wing_capped_dirichlet(; meshfile, AOA, magVinf)
    solve_simple_wing_capped_dirichlet!(body_ref, :backslash; backend)
    prepare_pressure_state!(body_ref, Vinf; backend,
        velocity_gradient=laplace_gradient_mode == :raw_hessian)

    body_b = deepcopy(body_ref)
    body_l = deepcopy(body_ref)
    body_kj = deepcopy(body_ref)

    frames_b = pnl.ReferenceFrame(body_b)
    frames_l = pnl.ReferenceFrame(body_l)
    frames_kj = pnl.ReferenceFrame(body_kj)

    chord = maximum(body_ref.nodes[1, :]) - minimum(body_ref.nodes[1, :])
    span = maximum(body_ref.nodes[2, :]) - minimum(body_ref.nodes[2, :])
    Sref = chord * span
    normalization = pnl.WingNormalization(rho, Sref, chord)

    pressure_bernoulli = pnl.PressureBernoulli(rho;
        unsteady=false,
        correct_kuttacondition=false,
        backend=backend)
    pressure_laplace = pnl.PressureLaplace((body_l,), rho;
        reference_panel=1,
        reference_pressure=0.0,
        gradient_mode=laplace_gradient_mode,
        verbose=false)
    force_bernoulli = pnl.ForceMonitor(1, 1;
        i_frame=-1,
        normalization=normalization,
        correct_kuttacondition=false,
        verbose=false)
    force_laplace = pnl.ForceMonitor(1, 1;
        i_frame=-1,
        normalization=normalization,
        correct_kuttacondition=false,
        verbose=false)
    kj_force = pnl.KuttaJoukowskiForce(body_kj, 1, 1;
        rho,
        backend,
        i_frame=-1,
        normalization=normalization,
        verbose=false)

    pressure_bernoulli((body_b,), (nothing,), frames_b, Vinf, 0, dt)
    force_bernoulli((body_b,), (nothing,), frames_b, Vinf, 0, dt)

    pressure_laplace.velocity_dot[1] .= expected_negative_tangent_velocity(body_l)
    pressure_laplace((body_l,), (nothing,), frames_l, Vinf, 0, dt)
    force_laplace((body_l,), (nothing,), frames_l, Vinf, 0, dt)

    kj_force((body_kj,), (nothing,), frames_kj, Vinf, 0, dt)

    p_b = body_b.P .- body_b.P[1]
    p_l = body_l.P .- body_l.P[1]
    pressure_abs_diff = maximum(abs.(p_b .- p_l))
    pressure_rel_diff = norm(p_b .- p_l) / max(norm(p_b), eps())

    cf_b = vec(force_bernoulli.force[:, 1])
    cf_l = vec(force_laplace.force[:, 1])
    cf_kj = vec(kj_force.force[:, 1])

    bernoulli_laplace_rel = norm(cf_b .- cf_l) / max(norm(cf_b), eps())
    bernoulli_kj_rel = norm(cf_b .- cf_kj) / max(norm(cf_b), eps())

    return (
        pressure_abs_diff=pressure_abs_diff,
        pressure_rel_diff=pressure_rel_diff,
        bernoulli_force=cf_b,
        laplace_force=cf_l,
        kj_force=cf_kj,
        bernoulli_laplace_rel=bernoulli_laplace_rel,
        bernoulli_kj_rel=bernoulli_kj_rel,
        npanels=body_ref.ncells,
        laplace_gradient_mode=laplace_gradient_mode,
    )
end

function main()
    println()
    println("#===== SIMPLE WING CAPPED PRESSURE COMPARISON =====#")

    result = compare_pressure_models()

    println("Panels: $(result.npanels)")
    println("PressureLaplace gradient mode: $(result.laplace_gradient_mode)")
    println("PressureBernoulli vs PressureLaplace:")
    println("\tmax |Δp| = $(result.pressure_abs_diff)")
    println("\trel ||Δp|| = $(result.pressure_rel_diff)")
    println()
    println("Integrated force coefficients [Cx, Cy, Cz]:")
    println("\tPressureBernoulli: $(result.bernoulli_force)")
    println("\tPressureLaplace:   $(result.laplace_force)")
    println("\tKuttaJoukowski:   $(result.kj_force)")
    println()
    println("Relative force differences:")
    println("\tBernoulli vs Laplace: $(result.bernoulli_laplace_rel)")
    println("\tBernoulli vs KuttaJoukowski: $(result.bernoulli_kj_rel)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
