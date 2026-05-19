## Compare Bernoulli vs Laplace vs Kutta-Joukowski force on the capped wing
## across BC formulations:
##   1. Dirichlet (Source+Doublet), gradient_robust=false
##   2. Dirichlet (Source+Doublet), gradient_robust=true (AR threshold=7)
##   3. Neumann   (VortexRing)
##
## Tests whether the Bernoulli<->Laplace gap is driven by `compute_mu_gradient!`
## sliver contamination (the current strongest hypothesis logged in
## DEBUG_PRESSURE_LAPLACE.md). Under Neumann the velocity is built directly
## from the vortex-ring kernel without `compute_mu_gradient!`, so if the
## hypothesis is right the gap should collapse.

import FLOWPanel as pnl
import GeoIO
using LinearAlgebra: norm, dot
using Printf

include(joinpath(pnl.examples_path, "simple_wing_capped_dirichlet.jl"))

const MESHFILE = joinpath(pnl.examples_path, "data", "naca0012_nc101_nw26.msh")

function build_neumann_body(; meshfile=MESHFILE, AOA=5.88, magVinf=56.0,
                              cp_offset=1e-10, kernel_offset=1e-3)
    msh = GeoIO.load(meshfile).geometry
    nodes, cells = pnl.meshes2nodes_cells(msh)

    # Remove one panel to make the closed surface non-watertight, so the
    # pure-VortexRing Neumann influence matrix is full rank. Drop a mid-mesh
    # panel (avoids cap/TE panels where dropping would distort the wake).
    drop_idx = size(cells, 2) ÷ 2
    cells = hcat(cells[:, 1:drop_idx-1], cells[:, drop_idx+1:end])
    println("\tNeumann: dropped panel index $drop_idx; new ncells = $(size(cells, 2))")

    xte = maximum(nodes[1, :])
    ystart = minimum(nodes[2, :])
    ystop = maximum(nodes[2, :])
    span = ystop - ystart

    trailingedge = zeros(eltype(nodes), 3, 10_000)
    trailingedge[1, :] .= xte
    trailingedge[2, :] .= range(ystart, stop=ystop, length=size(trailingedge, 2))

    Vinf = magVinf .* [cosd(AOA), 0.0, sind(AOA)]

    kernel = pnl.VortexRing
    DBC = pnl.kernel_dim(kernel) == 2
    bodytype = pnl.RigidWakeBody{kernel, pnl.kernel_dim(kernel), Float64, DBC}

    shedding = pnl.calc_shedding(nodes, cells, trailingedge; tolerance=0.001 * span)
    body = bodytype(nodes, cells, shedding;
        CPoffset=cp_offset, kerneloffset=kernel_offset,
        ensure_winding=true, semiinfinite_wake=true)
    shedding = pnl.calc_shedding(body.nodes, body.cells, trailingedge; tolerance=0.001 * span)
    body = bodytype(body.nodes, body.cells, shedding;
        CPoffset=cp_offset, kerneloffset=kernel_offset,
        ensure_winding=true, semiinfinite_wake=true)

    wake_direction = reshape(Vinf ./ norm(Vinf), :, 1)
    for i in eachindex(body.Das)
        body.Das[i] .= repeat(wake_direction, 1, size(body.Das[i], 2))
    end

    pnl.apply_freestream!(body, Vinf)

    return body, Vinf
end

function normal_velocity_stats(body)
    max_rel = 0.0
    sum_sq_rel = 0.0
    n = 0
    for p in 1:body.ncells
        u = view(body.velocity, :, p)
        nrm = view(body.normals, :, p)
        umag = norm(u)
        umag < eps() && continue
        un = abs(dot(u, nrm))
        rel = un / umag
        max_rel = max(max_rel, rel)
        sum_sq_rel += rel^2
        n += 1
    end
    mean_sq = n > 0 ? sum_sq_rel / n : 0.0
    return max_rel, mean_sq
end

function run_case(; case_name, body, Vinf, gradient_robust, gradient_ar_threshold=7.0,
                    rho=1.225, backend=pnl.DirectBackend(), dt=1.0)

    body.needs_velocity_gradient[] = true
    body.velocity .= 0.0
    body.velocity_gradient .= 0.0
    pnl.calcfield_U!(body, Vinf; backend,
        gradient_robust=gradient_robust,
        gradient_ar_threshold=gradient_ar_threshold)
    pnl.influence!((body,), (body,), backend;
        scalar_potential=false, velocity=false, velocity_gradient=true)

    max_un_rel, mean_un_sq = normal_velocity_stats(body)

    body_b = deepcopy(body)
    body_l = deepcopy(body)
    body_kj = deepcopy(body)
    frames_b = pnl.ReferenceFrame(body_b)
    frames_l = pnl.ReferenceFrame(body_l)
    frames_kj = pnl.ReferenceFrame(body_kj)

    chord = maximum(body.nodes[1, :]) - minimum(body.nodes[1, :])
    span = maximum(body.nodes[2, :]) - minimum(body.nodes[2, :])
    Sref = chord * span
    normalization = pnl.WingNormalization(rho, Sref, chord)

    pressure_bernoulli = pnl.PressureBernoulli(rho; unsteady=false,
        correct_kuttacondition=false, backend=backend)
    pressure_laplace = pnl.PressureLaplace((body_l,), rho;
        reference_panel=1, reference_pressure=0.0, verbose=false)
    force_bernoulli = pnl.ForceMonitor(1, 1; i_frame=-1, normalization,
        correct_kuttacondition=false, verbose=false)
    force_laplace   = pnl.ForceMonitor(1, 1; i_frame=-1, normalization,
        correct_kuttacondition=false, verbose=false)
    kj_force = pnl.KuttaJoukowskiForce(body_kj, 1, 1; rho, backend,
        i_frame=-1, normalization, verbose=false)

    pressure_bernoulli((body_b,), (nothing,), frames_b, Vinf, 0, dt)
    force_bernoulli((body_b,), (nothing,), frames_b, Vinf, 0, dt)
    pressure_laplace.velocity_dot[1] .= 0.0
    pressure_laplace((body_l,), (nothing,), frames_l, Vinf, 0, dt)
    force_laplace((body_l,), (nothing,), frames_l, Vinf, 0, dt)
    kj_force((body_kj,), (nothing,), frames_kj, Vinf, 0, dt)

    L = pressure_laplace.L[1]
    b_L = pressure_laplace.b[1]
    P_B = body_b.P .- body_b.P[1]
    rel_resid = norm(L * P_B .- b_L) / max(norm(b_L), eps())

    cf_b  = vec(force_bernoulli.force[:, 1])
    cf_l  = vec(force_laplace.force[:, 1])
    cf_kj = vec(kj_force.force[:, 1])

    return (case=case_name,
            CzB=cf_b[3], CzL=cf_l[3], CzKJ=cf_kj[3],
            BL_rel=norm(cf_b .- cf_l) / max(norm(cf_b), eps()),
            BKJ_rel=norm(cf_b .- cf_kj) / max(norm(cf_b), eps()),
            max_un_rel=max_un_rel, mean_un_sq=mean_un_sq,
            rel_resid=rel_resid)
end

function main()
    backend = pnl.DirectBackend()

    println()
    println("#===== BC comparison on $(basename(MESHFILE)) =====#")

    # Dirichlet body (built once, then deep-copied per case)
    body_dir1, Vinf_dir = build_simple_wing_capped_dirichlet(; meshfile=MESHFILE)
    solve_simple_wing_capped_dirichlet!(body_dir1, :backslash; backend)

    body_dir2 = deepcopy(body_dir1)

    # Neumann body
    body_neu, Vinf_neu = build_neumann_body(; meshfile=MESHFILE)
    solver_neu = pnl.Backslash(body_neu)
    pnl.solve!(body_neu, solver_neu; backend)

    results = [
        run_case(; case_name="Dirichlet",       body=body_dir1, Vinf=Vinf_dir,
                   gradient_robust=false, backend),
        run_case(; case_name="Dirichlet+AR=3",  body=body_dir2, Vinf=Vinf_dir,
                   gradient_robust=true, gradient_ar_threshold=3.0, backend),
        run_case(; case_name="Neumann",         body=body_neu,  Vinf=Vinf_neu,
                   gradient_robust=false, backend),
    ]

    @printf("\n%-18s %-10s %-10s %-10s %-10s %-10s %-14s %-10s\n",
        "case", "CzB", "CzL", "CzKJ", "BL_rel", "BKJ_rel", "max|u.n|/|u|", "rel_resid")
    for r in results
        @printf("%-18s %-10.4f %-10.4f %-10.4f %-10.3e %-10.3e %-14.3e %-10.3e\n",
            r.case, r.CzB, r.CzL, r.CzKJ, r.BL_rel, r.BKJ_rel, r.max_un_rel, r.rel_resid)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
