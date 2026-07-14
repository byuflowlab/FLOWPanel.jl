## Per-panel and per-edge diagnostic VTU writer for the capped Dirichlet wing.
##
## Writes <out_dir>/panel_diag.vtu with per-cell:
##   - aspect_ratio                  (max_edge / min_edge)
##   - min_edge, max_edge, area
##   - u_dot_n, abs_u_dot_n, rel_u_dot_n  (tangency residual)
##   - velocity_magnitude
##   - P_bernoulli, P_laplace        (both shifted to P[ref]=0)
##   - pressure_abs_diff             (|P_B - P_L|)
##   - laplace_panel_residual        (|L*P_B - b_L|_panel)
##
## And <out_dir>/edge_diag.vtu with per-edge:
##   - edge_length
##   - dihedral_angle_deg            (acos(n_i · n_j))
##   - conormal_magnitude            |ν_ij|
##   - bernoulli_edge_flux, laplace_edge_flux, edge_flux_diff
##   - panel_i_AR, panel_j_AR, max_panel_AR
##
## Open both in Paraview, color by aspect_ratio and abs_u_dot_n, and confirm
## that high-AR, high-tangency, and high-pressure-residual regions coincide.

import FLOWPanel as pnl
using LinearAlgebra: norm, dot
using WriteVTK

include(joinpath(pnl.examples_path, "simple_wing_capped_pressure_comparison.jl"))

cross3(a, b) = [a[2]*b[3] - a[3]*b[2],
                a[3]*b[1] - a[1]*b[3],
                a[1]*b[2] - a[2]*b[1]]

## Tangent-plane stencil conditioning for compute_mu_gradient!.
## Replicates the 1-ring LS stencil (without TE upstream expansion), projects
## the neighbor offsets into the panel tangent plane, and reports the 2x2 Gram
## matrix condition number. High κ means one tangent direction of ∇μ is
## unreliable — exactly the failure mode expected on sliver triangles.
function tangent_stencil_condition(body)
    ncells = body.ncells
    kappa = zeros(ncells)
    nstencil = zeros(Int, ncells)
    aniso = zeros(ncells)  # ratio max/min singular value of the 2D stencil
    for i in 1:ncells
        ni = body.normals[:, i]
        # build orthonormal tangent basis (t1, t2)
        ref = abs(ni[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
        t1 = ref .- dot(ref, ni) .* ni
        t1 ./= max(norm(t1), eps())
        t2 = cross3(ni, t1)

        G = zeros(2, 2)
        ns = 0
        for k in 1:size(body.neighbor, 1)
            j = body.neighbor[k, i]
            j <= 0 && continue
            d = body.controlpoints[:, j] .- body.controlpoints[:, i]
            u = dot(d, t1); v = dot(d, t2)
            G[1, 1] += u*u; G[1, 2] += u*v
            G[2, 1] += v*u; G[2, 2] += v*v
            ns += 1
        end
        nstencil[i] = ns
        if ns >= 2
            # eigenvalues of 2x2 SPD G
            tr = G[1,1] + G[2,2]
            det_ = G[1,1]*G[2,2] - G[1,2]*G[2,1]
            disc = max(tr*tr/4 - det_, 0.0)
            lam_max = tr/2 + sqrt(disc)
            lam_min = tr/2 - sqrt(disc)
            kappa[i] = lam_max / max(lam_min, eps())
            aniso[i] = sqrt(lam_max / max(lam_min, eps()))
        else
            kappa[i] = Inf; aniso[i] = Inf
        end
    end
    return kappa, aniso, nstencil
end

function panel_edge_metrics(body)
    ncells = body.ncells
    AR = zeros(ncells)
    emin = zeros(ncells); emax = zeros(ncells); area = zeros(ncells)
    for p in 1:ncells
        ns = body.cells[:, p]
        v1 = body.nodes[:, ns[1]]; v2 = body.nodes[:, ns[2]]; v3 = body.nodes[:, ns[3]]
        e1 = norm(v2 .- v1); e2 = norm(v3 .- v2); e3 = norm(v1 .- v3)
        emin[p] = min(e1, e2, e3)
        emax[p] = max(e1, e2, e3)
        AR[p] = emin[p] > 0 ? emax[p] / emin[p] : Inf
        area[p] = 0.5 * norm(cross3(v2 .- v1, v3 .- v1))
    end
    return AR, emin, emax, area
end

function write_diagnostics(; meshfile=joinpath(pnl.examples_path, "data", "wing_ar4_naca0016_5.msh"),
                             AOA=5.88, magVinf=56.0, rho=1.225, dt=1.0,
                             kernel_offset=1e-3,
                             backend=pnl.DirectBackend(),
                             out_dir=joinpath(@__DIR__, "results", "panel_diagnostics"))

    mkpath(out_dir)

    body = build_pressure_comparison_wing(; kernel_offset)
    Vinf = magVinf .* [cosd(AOA), 0.0, sind(AOA)]
    set_pressure_comparison_wake!(body, Vinf)
    pnl.steady!(body, pnl.ReferenceFrame(body), Vinf;
        body_solvers=pnl.Backslash(body), backend, verbose=false)

    body.needs_velocity_gradient[] = true
    body.velocity .= 0.0
    body.velocity_gradient .= 0.0
    pnl.calcfield_U!(body, Vinf; backend)
    pnl.influence!((body,), (body,), backend;
        scalar_potential=false, velocity=false, velocity_gradient=true)

    # --- Bernoulli copy
    body_b = deepcopy(body)
    frames_b = pnl.ReferenceFrame(body_b)
    pb_mon = pnl.PressureBernoulli(rho; unsteady=false,
        correct_kuttacondition=false, backend=backend)
    pb_mon((body_b,), (nothing,), frames_b, Vinf, 0, dt)

    # --- Laplace copy
    body_l = deepcopy(body)
    frames_l = pnl.ReferenceFrame(body_l)
    pl_mon = pnl.PressureLaplace((body_l,), rho;
        reference_panel=1, reference_pressure=0.0, verbose=false)
    for p in 1:body_l.ncells
        n = @view body_l.normals[:, p]
        u = @view body_l.velocity[:, p]
        tang = u .- dot(u, n) .* n
        pl_mon.velocity_dot[1][:, p] .= -(tang .+ body_l.velocity_kinematic[:, p])
    end
    pl_mon((body_l,), (nothing,), frames_l, Vinf, 0, dt)

    # --- Per-panel scalars
    AR, emin, emax, area = panel_edge_metrics(body)
    stencil_kappa, stencil_aniso, stencil_nbrs = tangent_stencil_condition(body)
    udotn = zeros(body.ncells); umag = zeros(body.ncells)
    for p in 1:body.ncells
        n = @view body.normals[:, p]
        u = @view body.velocity[:, p]
        udotn[p] = dot(u, n)
        umag[p] = norm(u)
    end
    abs_udotn = abs.(udotn)
    rel_udotn = abs_udotn ./ max.(umag, eps())

    P_B = body_b.P .- body_b.P[1]
    P_L = body_l.P .- body_l.P[1]
    L = pl_mon.L[1]; b_L = pl_mon.b[1]
    LP = L * P_B
    panel_resid = abs.(LP .- b_L)
    pressure_abs_diff = abs.(P_B .- P_L)

    # --- Write panel VTU
    vtk_path = joinpath(out_dir, "panel_diag")
    vtk_grid(vtk_path, body.nodes, body.vtk_cells) do vtk
        vtk["aspect_ratio", VTKCellData()] = AR
        vtk["min_edge", VTKCellData()] = emin
        vtk["max_edge", VTKCellData()] = emax
        vtk["area", VTKCellData()] = area
        vtk["u_dot_n", VTKCellData()] = udotn
        vtk["abs_u_dot_n", VTKCellData()] = abs_udotn
        vtk["rel_u_dot_n", VTKCellData()] = rel_udotn
        vtk["velocity_magnitude", VTKCellData()] = umag
        vtk["normals", VTKCellData()] = body.normals
        vtk["velocity", VTKCellData()] = body.velocity
        vtk["P_bernoulli", VTKCellData()] = P_B
        vtk["P_laplace", VTKCellData()] = P_L
        vtk["pressure_abs_diff", VTKCellData()] = pressure_abs_diff
        vtk["laplace_panel_residual", VTKCellData()] = panel_resid
        vtk["stencil_kappa", VTKCellData()] = stencil_kappa
        vtk["log10_stencil_kappa", VTKCellData()] = log10.(max.(stencil_kappa, eps()))
        vtk["stencil_anisotropy", VTKCellData()] = stencil_aniso
        vtk["stencil_nbr_count", VTKCellData()] = stencil_nbrs
    end
    println("Wrote $(vtk_path).vtu")

    # --- Per-edge data: build a polyline VTU
    edges = pl_mon.edges[1]
    nedges = size(edges, 2)
    edge_length = zeros(nedges)
    dihedral_deg = zeros(nedges)
    conormal_mag = zeros(nedges)
    bern_flux = zeros(nedges)
    lapl_flux = zeros(nedges)
    panel_i_AR = zeros(nedges); panel_j_AR = zeros(nedges); max_panel_AR = zeros(nedges)

    # Recompute the same acceleration the monitor used: stored in pl_mon.acceleration[1]
    a_panel = pl_mon.acceleration[1]

    edge_points = zeros(3, 2 * nedges)
    edge_cells = WriteVTK.MeshCell[]
    for k in 1:nedges
        na, nb, pi, pj = edges[1, k], edges[2, k], edges[3, k], edges[4, k]
        xa = body_l.nodes[:, na]; xb = body_l.nodes[:, nb]
        edge_points[:, 2k-1] = xa
        edge_points[:, 2k]   = xb
        push!(edge_cells, MeshCell(VTKCellTypes.VTK_LINE, [2k-1, 2k]))

        edge_vec = xb .- xa
        edge_length[k] = norm(edge_vec)

        ni = body_l.normals[:, pi]; nj = body_l.normals[:, pj]
        nij = 0.5 .* (ni .+ nj)
        # in-plane co-normal: edge tangent crossed with averaged normal
        t = edge_vec ./ max(edge_length[k], eps())
        nu = cross3(t, nij)
        # orient from i to j
        rij = body_l.controlpoints[:, pj] .- body_l.controlpoints[:, pi]
        if dot(nu, rij) < 0
            nu = -nu
        end
        conormal_mag[k] = norm(nu)

        cosang = clamp(dot(ni, nj), -1.0, 1.0)
        dihedral_deg[k] = acosd(cosang)

        # Bernoulli flux across edge (gradient-flux): w * (P_i - P_j) where w is the
        # operator weight (read directly from L for consistency).
        wij = -L[pi, pj]
        bern_flux[k] = wij * (P_B[pi] - P_B[pj])

        # Laplace RHS flux: ρ ℓ * 0.5*(a_i + a_j) · ν_ij
        a_avg = 0.5 .* (a_panel[:, pi] .+ a_panel[:, pj])
        lapl_flux[k] = rho * edge_length[k] * dot(a_avg, nu)

        panel_i_AR[k] = AR[pi]; panel_j_AR[k] = AR[pj]
        max_panel_AR[k] = max(AR[pi], AR[pj])
    end

    edge_path = joinpath(out_dir, "edge_diag")
    vtk_grid(edge_path, edge_points, edge_cells) do vtk
        vtk["edge_length", VTKCellData()] = edge_length
        vtk["dihedral_angle_deg", VTKCellData()] = dihedral_deg
        vtk["conormal_magnitude", VTKCellData()] = conormal_mag
        vtk["bernoulli_edge_flux", VTKCellData()] = bern_flux
        vtk["laplace_edge_flux", VTKCellData()] = lapl_flux
        vtk["edge_flux_diff", VTKCellData()] = bern_flux .- lapl_flux
        vtk["abs_edge_flux_diff", VTKCellData()] = abs.(bern_flux .- lapl_flux)
        vtk["panel_i_AR", VTKCellData()] = panel_i_AR
        vtk["panel_j_AR", VTKCellData()] = panel_j_AR
        vtk["max_panel_AR", VTKCellData()] = max_panel_AR
    end
    println("Wrote $(edge_path).vtu")

    # --- Quick text summary
    println()
    println("Summary:")
    println("  panels = $(body.ncells), edges = $nedges")
    println("  kerneloffset = $kernel_offset")
    println("  median AR = $(median(AR)), p99 AR = $(quantile(AR, 0.99)), max AR = $(maximum(AR))")
    println("  median |u·n| = $(median(abs_udotn)), p99 = $(quantile(abs_udotn, 0.99))")
    println("  median rel |u·n| = $(median(rel_udotn)), p99 = $(quantile(rel_udotn, 0.99))")
    println("  ||L*P_B - b_L|| / ||b_L|| = $(norm(LP .- b_L) / max(norm(b_L), eps()))")

    # --- AR vs rel_u_dot_n vs stencil_kappa correlation
    finite_k = isfinite.(stencil_kappa)
    if any(finite_k)
        k_finite = stencil_kappa[finite_k]
        r_finite = rel_udotn[finite_k]
        ar_finite = AR[finite_k]
        # Spearman-ish: rank correlation via sortperm
        function rankvec(x)
            o = sortperm(x); r = zeros(length(x)); r[o] .= 1:length(x); return r
        end
        rk = rankvec(k_finite); rr = rankvec(r_finite); ra = rankvec(ar_finite)
        function pearson(a, b)
            ma = sum(a)/length(a); mb = sum(b)/length(b)
            num = sum((a .- ma) .* (b .- mb))
            den = sqrt(sum((a .- ma).^2) * sum((b .- mb).^2))
            return num / max(den, eps())
        end
        println("  rank corr(rel_u·n, stencil_kappa) = $(pearson(rr, rk))")
        println("  rank corr(rel_u·n, AR)            = $(pearson(rr, ra))")
        println("  rank corr(AR, stencil_kappa)      = $(pearson(ra, rk))")
        println("  stencil_kappa quantiles: median=$(quantile(k_finite, 0.5)), p90=$(quantile(k_finite, 0.9)), p99=$(quantile(k_finite, 0.99)), max=$(maximum(k_finite))")
    end
    return nothing
end

using Statistics: median, quantile

if abspath(PROGRAM_FILE) == @__FILE__
    write_diagnostics()
end
