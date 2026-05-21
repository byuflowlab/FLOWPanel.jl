## Sweep `gradient_ar_threshold` to see whether a tighter mask closes the
## PressureLaplace vs PressureBernoulli gap on the unrefined capped wing.

import FLOWPanel as pnl
using LinearAlgebra: norm, dot
using Printf

include(joinpath(pnl.examples_path, "simple_wing_capped_dirichlet.jl"))

function expected_negative_tangent_velocity(body)
    expected = similar(body.velocity)
    for p in 1:body.ncells
        n = body.normals[:, p]
        u = body.velocity[:, p]
        expected[:, p] .= -(u .- dot(u, n) .* n .+ body.velocity_kinematic[:, p])
    end
    return expected
end

function run_one(; thresh, AOA=5.88, magVinf=56.0, rho=1.225,
                   backend=pnl.DirectBackend(), dt=1.0,
                   meshfile=joinpath(pnl.examples_path, "data", "wing_ar4_naca0016_5.msh"))

    body_ref, Vinf = build_simple_wing_capped_dirichlet(; meshfile, AOA, magVinf)
    solve_simple_wing_capped_dirichlet!(body_ref, :backslash; backend)

    body_ref.needs_velocity_gradient[] = true
    body_ref.velocity .= 0.0
    body_ref.velocity_gradient .= 0.0
    pnl.calcfield_U!(body_ref, Vinf; backend, gradient_ar_threshold=thresh)
    pnl.influence!((body_ref,), (body_ref,), backend;
        scalar_potential=false, velocity=false, velocity_gradient=true)

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

    pressure_bernoulli = pnl.PressureBernoulli(rho; unsteady=false,
        correct_kuttacondition=false, backend=backend)
    pressure_laplace = pnl.PressureLaplace((body_l,), rho;
        reference_panel=1, reference_pressure=0.0, verbose=false)
    force_bernoulli = pnl.ForceMonitor(1, 1; i_frame=-1, normalization, correct_kuttacondition=false, verbose=false)
    force_laplace   = pnl.ForceMonitor(1, 1; i_frame=-1, normalization, correct_kuttacondition=false, verbose=false)
    kj_force = pnl.KuttaJoukowskiForce(body_kj, 1, 1; rho, backend, i_frame=-1, normalization, verbose=false)

    pressure_bernoulli((body_b,), (nothing,), frames_b, Vinf, 0, dt)
    force_bernoulli((body_b,), (nothing,), frames_b, Vinf, 0, dt)
    pressure_laplace.velocity_dot[1] .= 0.0
    pressure_laplace((body_l,), (nothing,), frames_l, Vinf, 0, dt)
    force_laplace((body_l,), (nothing,), frames_l, Vinf, 0, dt)
    kj_force((body_kj,), (nothing,), frames_kj, Vinf, 0, dt)

    # residual diagnostic
    L = pressure_laplace.L[1]
    b_L = pressure_laplace.b[1]
    P_B = body_b.P .- body_b.P[1]
    rel_resid = norm(L * P_B .- b_L) / max(norm(b_L), eps())

    cf_b = vec(force_bernoulli.force[:, 1])
    cf_l = vec(force_laplace.force[:, 1])
    cf_kj = vec(kj_force.force[:, 1])

    mask = pnl.panel_aspect_ratio_mask(body_ref.nodes, body_ref.cells; threshold=thresh)
    return (thresh=thresh, nbad=count(mask), rel_resid=rel_resid,
            cf_b=cf_b, cf_l=cf_l, cf_kj=cf_kj,
            bl_rel=norm(cf_b .- cf_l) / max(norm(cf_b), eps()),
            bkj_rel=norm(cf_b .- cf_kj) / max(norm(cf_b), eps()))
end

function main()
    println()
    println("#===== AR threshold sweep =====#")
    @printf("%-8s %-6s %-12s %-12s %-12s %-12s %-12s\n",
        "thresh", "nbad", "rel_resid", "CzB", "CzL", "BL_rel", "BKJ_rel")
    for thr in (3.0, 5.0, 7.0, 10.0, 15.0, Inf)
        r = run_one(; thresh=thr)
        @printf("%-8.2f %-6d %-12.3e %-12.4f %-12.4f %-12.3e %-12.3e\n",
            thr, r.nbad, r.rel_resid, r.cf_b[3], r.cf_l[3], r.bl_rel, r.bkj_rel)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
