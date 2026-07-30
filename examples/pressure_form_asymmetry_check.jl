#=##############################################################################
# DESCRIPTION (T1 of plans/20260709_pressure_monitor_reliability.md)
#   Wake-free discrimination of the two PressureLaplace RHS forms
#   (:material_derivative vs :lamb_vector) on steady analytic bodies where the
#   true fluid vorticity at the control points is 0, so any pressure difference
#   between the forms is discretization/projection asymmetry (mechanisms 2+3 of
#   BRAINSTORM/007) plus, on the lifting body only, the bound-sheet kappa
#   injection (mechanism 1c):
#     - sphere (NonLiftingBody, ConstantSource): exact analytic Cp available,
#       so each form also gets an absolute error;
#     - 45-degree swept NACA0012 wing (RigidWakeBody, VortexRing, AR=5,
#       AOA=4.2 deg): lifting case, forms compared against each other and
#       via CL.
#   NOTE: the plan called for the Weber RAE101 wing via
#   examples/helper_functions.jl, but that helper requires GeometricTools,
#   which is not installed in this workspace; the wing here is an equivalent
#   analytic loft (same AR, sweep, AOA) that serves the same purpose.
#
#   Mesh sweep (4,16) / (8,32) / (16,64); errors reported vs ncells. All
#   pressures are gauge fields (Laplace pins panel 1 to 0), so comparisons are
#   made after removing the mean of each field.
#
# Usage: julia -t auto --project examples/pressure_form_asymmetry_check.jl
=###############################################################################

import FLOWPanel as pnl
import LinearAlgebra as LA
using Printf

const OUTDIR = joinpath("data", "rotor_hover_pressure_reliability")
mkpath(OUTDIR)

demean(v) = v .- (sum(v) / length(v))
l2(v) = sqrt(sum(abs2, v) / length(v))
linf(v) = maximum(abs, v)

# ---------------- shared runner ------------------------------------------------
# One steady! per body with both Laplace forms + force monitors in a single
# monitor tuple (matches how they coexist in the rotor replay).

function run_both_forms(body, frames, Vinf, rho, normalization)
    pl_mat = pnl.PressureLaplace((body,), rho; reference_panel=1,
        reference_pressure=0.0, acceleration_form=:material_derivative,
        verbose=false, file=false)
    fm_mat = pnl.ForceMonitor(1, 1; i_frame=-1, normalization,
        correct_kuttacondition=false, verbose=false, file=false, vtk_fields=())
    # DIAGNOSTIC ONLY: deprecated Lamb mode is the subject of this asymmetry audit.
    pl_lamb = pnl.PressureLaplace((body,), rho; reference_panel=1,
        reference_pressure=0.0, acceleration_form=:lamb_vector,
        verbose=false, file=false)
    fm_lamb = pnl.ForceMonitor(1, 1; i_frame=-1, normalization,
        correct_kuttacondition=false, verbose=false, file=false, vtk_fields=())

    # tri basis for the bound-vorticity mu gradient: the steady! default (quad)
    # produces empty stencils at the sharp trailing edge of the analytic loft.
    pnl.steady!(body, frames, Vinf;
        body_solvers=pnl.Backslash(body),
        backend=pnl.DirectBackend(),
        monitors=(pl_mat, fm_mat, pl_lamb, fm_lamb),
        grad_mu_options=(; basis=:tri),
        path=nothing, verbose=false)

    return (; p_mat=copy(pl_mat.p[1]), p_lamb=copy(pl_lamb.p[1]),
        F_mat=copy(fm_mat.force[:, 1]), F_lamb=copy(fm_lamb.force[:, 1]),
        omega_rms=l2(vec(sqrt.(sum(abs2, body.induced_vorticity; dims=1)))))
end

# ---------------- geometry builders ---------------------------------------------

"Lat-long triangulated sphere of radius R centered at the origin (open polar
caps like examples/sphere.jl, theta in [0.15, pi-0.15])."
function sphere_nodes_cells(R, n_theta, n_phi)
    th = range(0.15, pi - 0.15; length=n_theta + 1)
    ph = range(0, 2pi; length=n_phi + 1)[1:end - 1]
    nid(i, j) = (i - 1) * n_phi + j          # i=1..n_theta+1, j=1..n_phi
    nodes = zeros(3, (n_theta + 1) * n_phi)
    for i in 1:n_theta + 1, j in 1:n_phi
        s, c = sincos(th[i])
        nodes[:, nid(i, j)] .= (R * c, R * s * cos(ph[j]), R * s * sin(ph[j]))
    end
    cells = zeros(Int, 3, 2 * n_theta * n_phi)
    k = 0
    for i in 1:n_theta, j in 1:n_phi
        jp = j == n_phi ? 1 : j + 1
        a, b, c, d = nid(i, j), nid(i, jp), nid(i + 1, jp), nid(i + 1, j)
        # (a, c, b)/(a, d, c) ordering gives outward normals (verified against
        # the exterior solution: |u|max = 1.5 Vinf at the equator)
        cells[:, k += 1] .= (a, c, b)
        cells[:, k += 1] .= (a, d, c)
    end
    return nodes, cells
end

naca0012_halfthickness(xc) = 5 * 0.12 * (0.2969 * sqrt(xc) - 0.1260 * xc -
    0.3516 * xc^2 + 0.2843 * xc^3 - 0.1036 * xc^4)

"Closed chordwise loop of 2nc nodes (sharp TE): TE -> upper -> LE -> lower."
function naca_loop(nc)
    xs = [0.5 * (1 + cos(pi * i / nc)) for i in 0:nc]   # 1 -> 0, TE to LE
    loop = zeros(2, 2 * nc)
    for i in 0:nc                                        # upper surface
        i == nc && continue                              # LE handled below
        loop[:, i + 1] .= (xs[i + 1], naca0012_halfthickness(xs[i + 1]))
    end
    loop[:, nc + 1] .= (0.0, 0.0)                        # LE
    for i in nc-1:-1:1                                   # lower surface LE->TE
        loop[:, nc + 1 + (nc - i)] .= (xs[i + 1], -naca0012_halfthickness(xs[i + 1]))
    end
    return loop
end

"45-deg swept, untapered, untwisted NACA0012 wing, AR=5, open tips."
function sweptwing_nodes_cells(b, c, sweep_deg, nc, nspan_full)
    loop = naca_loop(nc)
    nloop = size(loop, 2)
    ys = [b / 2 * cos(pi * (1 - k / nspan_full)) for k in 0:nspan_full]
    nsec = length(ys)
    nid(s, k) = (s - 1) * nloop + k
    nodes = zeros(3, nsec * nloop)
    tansw = tand(sweep_deg)
    for s in 1:nsec, k in 1:nloop
        x_le = tansw * abs(ys[s])
        nodes[:, nid(s, k)] .= (x_le + c * loop[1, k], ys[s], c * loop[2, k])
    end
    cells = zeros(Int, 3, 2 * (nsec - 1) * nloop)
    m = 0
    for s in 1:nsec - 1, k in 1:nloop
        kp = k == nloop ? 1 : k + 1
        a, bb, cc, d = nid(s, k), nid(s, kp), nid(s + 1, kp), nid(s + 1, k)
        cells[:, m += 1] .= (a, cc, bb)
        cells[:, m += 1] .= (a, d, cc)
    end
    return nodes, cells, ys
end

# ---------------- sphere case ---------------------------------------------------

function run_sphere(n_theta::Int, n_phi::Int)
    R = 1.0
    magVinf = 10.0
    Vinf = magVinf * [1.0, 0.0, 0.0]
    rho = 1.225
    nodes, cells = sphere_nodes_cells(R, n_theta, n_phi)
    body = pnl.NonLiftingBody{Union{pnl.ConstantSource}}(nodes, cells)
    frames = pnl.ReferenceFrame(body)
    normalization = pnl.WingNormalization(rho, pi * R^2, R)

    res = run_both_forms(body, frames, Vinf, rho, normalization)

    # analytic gauge pressure at each control point: p = q (1 - 9/4 sin^2 th)
    q = 0.5 * rho * magVinf^2
    p_exact = [begin
            x = body.controlpoints[1, i]
            r = sqrt(sum(abs2, view(body.controlpoints, :, i)))
            costh = x / r
            q * (1 - 2.25 * (1 - costh^2))
        end for i in 1:body.ncells]

    d_forms = demean(res.p_mat) .- demean(res.p_lamb)
    e_mat = demean(res.p_mat) .- demean(p_exact)
    e_lamb = demean(res.p_lamb) .- demean(p_exact)
    return (; case="sphere", n1=n_theta, n2=n_phi, ncells=body.ncells, q,
        res.omega_rms,
        d_forms_l2=l2(d_forms), d_forms_linf=linf(d_forms),
        e_mat_l2=l2(e_mat), e_mat_linf=linf(e_mat),
        e_lamb_l2=l2(e_lamb), e_lamb_linf=linf(e_lamb),
        dF=res.F_mat .- res.F_lamb, F_mat=res.F_mat, F_lamb=res.F_lamb)
end

# ---------------- swept-wing case ------------------------------------------------

function run_wing(n_rfl::Int, n_span::Int)
    AOA = 4.2
    magVinf = 30.0
    Vinf = magVinf * [cosd(AOA), 0.0, sind(AOA)]
    rho = 1.225
    b = 98 * 0.0254
    ar = 5.0
    c = b / ar
    Sref = b * c
    q = 0.5 * rho * magVinf^2
    nc = 3 * n_rfl                    # chordwise panels per side
    nspan_full = 2 * n_span           # spanwise panels across the full span

    nodes, cells, _ = sweptwing_nodes_cells(b, c, 45.0, nc, nspan_full)
    bodytype = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}

    # Shedding must be computed from the *constructed* (re-wound) cells; build a
    # noshedding body first, then rebuild with shedding from its nodes/cells.
    body0 = bodytype(nodes, cells, pnl.noshedding;
        ensure_winding=true, watertight=false)
    # Analytic distance-to-TE criterion: the loft's TE nodes lie exactly on
    # x = c + tan(sweep)|y|, z = 0, so a tight tolerance selects only true TE
    # edges (a sampled polyline with a loose tolerance also catches the first
    # cosine-clustered off-TE node row on fine meshes).
    tansw = tand(45.0)
    te_criterion(X) = (sqrt((X[1] - (c + tansw * abs(X[2])))^2 + X[3]^2), X[2])
    shedding = pnl.calc_shedding(body0.nodes, body0.cells, te_criterion;
        tolerance=1e-6 * c)
    body = bodytype(body0.nodes, body0.cells, shedding;
        ensure_winding=true, watertight=false)
    for i in eachindex(body.Das)
        body.Das[i] .= repeat(Vinf ./ magVinf, 1, size(body.Das[i], 2))
    end

    frames = pnl.ReferenceFrame(body)
    normalization = pnl.WingNormalization(rho, Sref, c)

    res = run_both_forms(body, frames, Vinf, rho, normalization)

    Dhat = Vinf / magVinf
    Lhat = LA.cross(Dhat, [0.0, 1.0, 0.0])
    d_forms = demean(res.p_mat) .- demean(res.p_lamb)
    return (; case="wing", n1=n_rfl, n2=n_span, ncells=body.ncells, q,
        res.omega_rms,
        d_forms_l2=l2(d_forms), d_forms_linf=linf(d_forms),
        e_mat_l2=NaN, e_mat_linf=NaN, e_lamb_l2=NaN, e_lamb_linf=NaN,
        dF=res.F_mat .- res.F_lamb, F_mat=res.F_mat, F_lamb=res.F_lamb,
        CL_mat=LA.dot(res.F_mat, Lhat), CL_lamb=LA.dot(res.F_lamb, Lhat))
end

# ---------------- sweep ----------------------------------------------------------

resolutions = [(4, 16), (8, 32), (16, 64)]

rows = []
for (n1, n2) in resolutions
    println("sphere ($n1,$n2) ...")
    push!(rows, run_sphere(n1, n2))
end
for (n1, n2) in resolutions
    println("wing ($n1,$n2) ...")
    push!(rows, run_wing(n1, n2))
end

println("\n=== form asymmetry sweep (pressures de-meaned; normalized by q) ===")
@printf("%-7s %10s %7s %10s %13s %13s %12s %12s %12s %12s\n",
    "case", "(n1,n2)", "ncells", "rms|omega|", "d_forms_L2/q", "d_forms_Loo/q",
    "e_mat_L2/q", "e_lamb_L2/q", "|dF|", "|F_mat|")
for r in rows
    @printf("%-7s %10s %7d %10.3e %13.4e %13.4e %12.4e %12.4e %12.4e %12.4e\n",
        r.case, "($(r.n1),$(r.n2))", r.ncells, r.omega_rms,
        r.d_forms_l2 / r.q, r.d_forms_linf / r.q,
        r.e_mat_l2 / r.q, r.e_lamb_l2 / r.q,
        LA.norm(r.dF), LA.norm(r.F_mat))
end

println("\n=== wing CL by form ===")
for r in rows
    r.case == "wing" || continue
    @printf("wing (%2d,%2d) ncells=%6d  CL_mat=%+.5f  CL_lamb=%+.5f  dCL=%+.5f\n",
        r.n1, r.n2, r.ncells, r.CL_mat, r.CL_lamb, r.CL_mat - r.CL_lamb)
end

csv_path = joinpath(OUTDIR, "T1_form_asymmetry_sweep.csv")
open(csv_path, "w") do io
    println(io, "case,n1,n2,ncells,q,omega_rms,d_forms_l2,d_forms_linf,e_mat_l2,e_mat_linf,e_lamb_l2,e_lamb_linf,dFx,dFy,dFz,CL_mat,CL_lamb")
    for r in rows
        cl_mat = r.case == "wing" ? r.CL_mat : NaN
        cl_lamb = r.case == "wing" ? r.CL_lamb : NaN
        println(io, join([r.case, r.n1, r.n2, r.ncells, r.q, r.omega_rms,
            r.d_forms_l2, r.d_forms_linf, r.e_mat_l2, r.e_mat_linf,
            r.e_lamb_l2, r.e_lamb_linf, r.dF[1], r.dF[2], r.dF[3],
            cl_mat, cl_lamb], ","))
    end
end
println("\nWrote $(csv_path)")
