#=##############################################################################
# DESCRIPTION
#   Zero-new-code test of the existing sliver fallback in compute_mu_gradient!:
#   does flagging TE slivers with panel_aspect_ratio_mask and using the
#   BFS-gathered healthy-only LS (with a depth large enough to escape the
#   flagged TE band) arrest the 1/h_min divergence of the tri half-jump?
#   See "Chordwise divergence root cause" in data/sweptwing_sweeploft/findings.md.
#
#   Run single-threaded:  julia --project examples/sweptwing_bfsgrad.jl
#
#   Environment knobs:
#     FLOWPANEL_BFSGRAD_CASES="40,80,156"
#     FLOWPANEL_BFSGRAD_NSPAN=48
#     FLOWPANEL_BFSGRAD_AR=10.0          aspect-ratio flag threshold
#     FLOWPANEL_BFSGRAD_DEPTH=60         bfs_max_depth (default 4 is too
#                                        shallow: ~26 flagged rows at n_ch=156)
#     FLOWPANEL_BFSGRAD_TARGET=8         bfs_target_healthy
=###############################################################################

ENV["FLOWPANEL_CHORDDIV_RUN"] = "false"
include(joinpath(@__DIR__, "sweptwing_chorddivergence.jl"))

function run_bfs_case(n_ch::Int, n_span::Int; swap_diagonals::Bool,
                      ar_threshold::Real, bfs_max_depth::Int,
                      bfs_target_healthy::Int, variant_name::String="")
    t0 = time()
    body = sweeploft_wing(n_ch, n_span; mirror_diagonals=true, swap_diagonals)
    N = body.ncells
    npts = 2*n_ch

    solver = pnl.Backslash(body)
    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false)
    force = pnl.ForceMonitor(1, 1; i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c_ref),
        correct_kuttacondition=false, verbose=false)
    Dhat = Vinf/LA.norm(Vinf)
    Lhat = LA.cross(Dhat, [0.0, 1.0, 0.0])

    pnl.steady!(body, pnl.ReferenceFrame(body), Vinf;
        body_solvers=solver, backend=pnl.FastMultipoleBackend(),
        monitors=(pressure, force), path=nothing, verbose=false)
    CL_solver = LA.dot(force.force[:, 1], Lhat)
    solver = nothing; GC.gc()

    mu = view(body.strength, :, pnl.get_Gammai(body))
    te_info = view(body.shedding_full, 1:2, :)

    # Production reconstruction (diverges) and the BFS-fallback variant under test
    Jtri = pnl.compute_mu_gradient!(zeros(3, N), body.controlpoints,
        body.normals, body.cells, body.neighbor, mu, te_info;
        scale=0.5, nodes=body.nodes, grad_mu_options=(; basis=:tri))
    Jbfs = pnl.compute_mu_gradient!(zeros(3, N), body.controlpoints,
        body.normals, body.cells, body.neighbor, mu, te_info;
        scale=0.5, nodes=body.nodes,
        grad_mu_options=(; basis=:tri, tri_robust=true,
            tri_robust_ar_threshold=ar_threshold,
            tri_robust_max_depth=bfs_max_depth,
            tri_robust_target_healthy=bfs_target_healthy))

    U0 = body.velocity .- Jtri                       # = U∞ + PV
    geom = quad_geometry(body)
    res_tri = quad_velocity_CL(geom, U0, 2 .* quad_average(geom, Jtri), +1.0,
                               Lhat, Dhat)
    res_bfs = quad_velocity_CL(geom, U0, 2 .* quad_average(geom, Jbfs), +1.0,
                               Lhat, Dhat)

    strip(ci) = mod1((ci + 1) ÷ 2, npts)
    te_cells = [ci for ci in 1:N if strip(ci) == 1 || strip(ci) == npts]
    Jtri_te_mean, Jtri_te_max = colnorm_stats(Jtri, te_cells)
    Jbfs_te_mean, Jbfs_te_max = colnorm_stats(Jbfs, te_cells)
    Jbfs_all_mean, Jbfs_all_max = colnorm_stats(Jbfs, 1:N)
    h_min = c*(1 - cos(pi/n_ch))/2

    row = (; variant=variant_name, n_ch, n_span, panels=N, h_min,
        ar_threshold, bfs_max_depth, bfs_target_healthy,
        n_flagged=count(mask), flagged_frac=count(mask)/N,
        CL_solver, CL_tri=res_tri.CL, CL_bfs=res_bfs.CL,
        Jtri_te_mean, Jtri_te_max, Jbfs_te_mean, Jbfs_te_max,
        Jbfs_all_mean, Jbfs_all_max, elapsed=time() - t0)
    return row, res_bfs.U
end

function main_bfsgrad()
    cases = [parse(Int, s) for s in
        split(get(ENV, "FLOWPANEL_BFSGRAD_CASES", "40,80,156"), ",")]
    n_span = parse(Int, get(ENV, "FLOWPANEL_BFSGRAD_NSPAN", "48"))
    ar_threshold = parse(Float64, get(ENV, "FLOWPANEL_BFSGRAD_AR", "10.0"))
    bfs_max_depth = parse(Int, get(ENV, "FLOWPANEL_BFSGRAD_DEPTH", "60"))
    bfs_target_healthy = parse(Int, get(ENV, "FLOWPANEL_BFSGRAD_TARGET", "8"))
    out_csv = joinpath("data", "sweptwing_sweeploft", "bfsgrad_test.csv")

    println("BFS-fallback ∇μ test: n_ch = ", cases, " x n_span = ", n_span,
            ", AR threshold = ", ar_threshold, ", depth = ", bfs_max_depth,
            ", target_healthy = ", bfs_target_healthy)

    rows = []
    for n_ch in cases
        U = Dict{String, Any}()
        for (vname, swap) in (("diagAC", false), ("diagBD", true))
            r, Ubfs = run_bfs_case(n_ch, n_span; swap_diagonals=swap,
                ar_threshold, bfs_max_depth, bfs_target_healthy,
                variant_name=vname)
            @printf("%7s %4d flagged %6d (%4.1f%%) CLslv %+.6f CLtri %+.6f CLbfs %+.6f | Jte_max tri %9.3e bfs %9.3e (bfs*h %9.3e) t %.0fs\n",
                r.variant, r.n_ch, r.n_flagged, 100*r.flagged_frac,
                r.CL_solver, r.CL_tri, r.CL_bfs,
                r.Jtri_te_max, r.Jbfs_te_max, r.Jbfs_te_max*r.h_min, r.elapsed)
            flush(stdout)
            push!(rows, r)
            CSV.write(out_csv, DataFrame(rows))
            U[vname] = Ubfs
            GC.gc()
        end
        @printf("        %4d cross-variant max|dU_bfs|/Vinf: %.3f\n",
            n_ch, maxdiff(U["diagAC"], U["diagBD"])/magVinf)
        flush(stdout)
    end
    println("CSV written to ", out_csv)
    return rows
end

if !isinteractive() && get(ENV, "FLOWPANEL_BFSGRAD_RUN", "true") == "true"
    bfs_rows = main_bfsgrad()
end
