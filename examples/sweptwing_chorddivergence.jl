#=##############################################################################
# DESCRIPTION
#   Root-cause analysis of the chordwise CL divergence of the triangulated
#   sweep-lofted Weber wing (see data/sweptwing_sweeploft/PLAN_chord_divergence.md
#   and findings.md).
#
#   Chord ladder at fixed n_span with component decomposition of the surface
#   velocity (PV vs tri half-jump ½∇μ), four ∇μ reconstructions (tri-LS jump,
#   structured quad LS, area-weighted P1, angle-weighted P1), Kutta-Joukowski
#   CL from γ alone, condition number of the influence matrix, flow-tangency
#   residuals, and TE-row sliver diagnostics.
#
#   Run single-threaded:  julia --project examples/sweptwing_chorddivergence.jl
#   (do NOT use julia -t: steady! calls FastMultipole.fmm!, whose threaded
#   path crashes with a MethodError on extra_farfield!)
#
#   Environment knobs:
#     FLOWPANEL_CHORDDIV_CASES="40,56,80,112,156"   chord levels n_ch
#     FLOWPANEL_CHORDDIV_NSPAN=48
#     FLOWPANEL_CHORDDIV_VARIANTS="diagAC,diagBD"
#     FLOWPANEL_CHORDDIV_KERNELOFFSETS=""           e.g. "1e-8,1e-12" enables
#                                                   the section-C probe
#     FLOWPANEL_CHORDDIV_KO_CASES="80,112"          chord levels for the probe
=###############################################################################

ENV["FLOWPANEL_SWEEPLOFT_RUN"] = "false"
include(joinpath(@__DIR__, "sweptwing_sweeploft.jl"))

# Keeping a second copy of G (for the true algebraic residual ||G*gamma-rhs||)
# only fits in RAM up to this many panels (2 dense matrices: 2*8*N^2 bytes;
# 16 GB machine -> cap at N=23000, ~8.5 GB for the pair).
const RESID_COPY_MAX_PANELS = 23000

# ----------------- ANGLE-WEIGHTED P1 GRADIENT (plan section B) ----------------
"""
Interior angle of the triangle at vertex `p` (edges p->q, p->r).
"""
function vertex_angle(p, q, r)
    u = q .- p
    v = r .- p
    nu, nv = LA.norm(u), LA.norm(v)
    (nu <= 0 || nv <= 0) && return 0.0
    return acos(clamp(LA.dot(u, v)/(nu*nv), -1.0, 1.0))
end

"""
    p1_mu_gradient_angle(body, n_ch, n_span) -> 3 x Nq

Same as `p1_mu_gradient` but the nodal averaging weight is the triangle's
vertex ANGLE at that node instead of the triangle AREA, so a sliver triangle
contributes ~0 weight at its two far (acute) vertices instead of polluting
them with its full area. TE-seam upper/lower register split and the
per-triangle P1 gradient + area-weighted quad averaging are unchanged.
"""
function p1_mu_gradient_angle(body, n_ch::Int, n_span::Int)
    npts = 2*n_ch
    nst = n_span + 1
    Nq = body.ncells ÷ 2
    areas = pnl.calc_areas(body)
    gamma = view(body.strength, :, pnl.get_Gammai(body))

    strip(ci) = mod1((ci + 1) ÷ 2, npts)
    is_upper(ci) = strip(ci) <= n_ch
    te_node = falses(size(body.nodes, 2))
    for k in 1:nst
        te_node[(k - 1)*npts + 1] = true
    end

    nnodes = size(body.nodes, 2)
    val = zeros(nnodes); wt = zeros(nnodes)            # regular register
    val_lo = zeros(nnodes); wt_lo = zeros(nnodes)      # lower-side register (TE nodes)
    for ci in 1:body.ncells
        g = gamma[ci]
        n1, n2, n3 = body.cells[1, ci], body.cells[2, ci], body.cells[3, ci]
        p1 = view(body.nodes, :, n1)
        p2 = view(body.nodes, :, n2)
        p3 = view(body.nodes, :, n3)
        angles = (vertex_angle(p1, p2, p3),
                  vertex_angle(p2, p3, p1),
                  vertex_angle(p3, p1, p2))
        for (ni, a) in zip((n1, n2, n3), angles)
            if te_node[ni] && !is_upper(ci)
                val_lo[ni] += a*g; wt_lo[ni] += a
            else
                val[ni] += a*g; wt[ni] += a
            end
        end
    end
    mu_node = val ./ max.(wt, eps())
    mu_node_lo = ifelse.(wt_lo .> 0, val_lo ./ max.(wt_lo, eps()), mu_node)

    grad = zeros(3, Nq)
    for q in 1:Nq
        Aq = 0.0
        for t in (2q - 1, 2q)
            n1, n2, n3 = body.cells[1, t], body.cells[2, t], body.cells[3, t]
            f(ni) = (te_node[ni] && !is_upper(t)) ? mu_node_lo[ni] : mu_node[ni]
            gt = tri_p1_gradient(view(body.nodes, :, n1), view(body.nodes, :, n2),
                                 view(body.nodes, :, n3), f(n1), f(n2), f(n3))
            grad[:, q] .+= areas[t] .* gt
            Aq += areas[t]
        end
        grad[:, q] ./= Aq
    end
    return grad
end

# ----------------- HELPERS ----------------------------------------------------
"""Area-weighted average of a 3 x ncells tri field onto quads (3 x Nq)."""
function quad_average(geom, X)
    Nq = length(geom.A)
    out = zeros(3, Nq)
    for q in 1:Nq
        a1, a2 = geom.areas[2q-1], geom.areas[2q]
        out[:, q] .= (a1 .* view(X, :, 2q-1) .+ a2 .* view(X, :, 2q)) ./ (a1 + a2)
    end
    return out
end

"""Mean and max of column norms of `X` restricted to cells in `idx`."""
function colnorm_stats(X, idx)
    isempty(idx) && return (0.0, 0.0)
    s, mx = 0.0, 0.0
    for i in idx
        v = sqrt(X[1, i]^2 + X[2, i]^2 + X[3, i]^2)
        s += v
        mx = max(mx, v)
    end
    return s/length(idx), mx
end

"""Mean and max of |x| over `idx`."""
function abs_stats(x, idx)
    isempty(idx) && return (0.0, 0.0)
    s, mx = 0.0, 0.0
    for i in idx
        v = abs(x[i])
        s += v
        mx = max(mx, v)
    end
    return s/length(idx), mx
end

"""Max column-norm difference between two 3 x Nq fields."""
maxdiff(A, B) = maximum(sqrt.(vec(sum(abs2, A .- B; dims=1))))

# ----------------- ONE CASE ---------------------------------------------------
"""
Solve one (n_ch, n_span, variant) case with a manually assembled dense
Backslash solver (so cond1 and the algebraic residual are recoverable), then
decompose body.velocity into U0 = U∞+PV and the tri half-jump Jtri, evaluate
all four ∇μ reconstructions, and collect the plan's A/B/B2 diagnostics.
Returns `(row::NamedTuple, Ufields::NamedTuple)`.
"""
function run_case(n_ch::Int, n_span::Int; swap_diagonals::Bool,
                  core_size::Real=1e-10, variant_name::String="")
    t0 = time()
    body = sweeploft_wing(n_ch, n_span; mirror_diagonals=true, swap_diagonals,
                          core_size)
    N = body.ncells
    npts = 2*n_ch

    # --- B2.2: assemble G ourselves so the 1-norm survives the in-place lu! ---
    # (pnl.Backslash's constructor destroys G: Glu = lu!(G).)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    G = zeros(N, N)
    pnl._G!(G, body, body; core_size=body.core_size_panel,
            update_geometry=false)
    anorm = LA.opnorm(G, 1)
    Gres = N <= RESID_COPY_MAX_PANELS ? copy(G) : nothing
    Glu = LA.lu!(G)
    rcond = LA.LAPACK.gecon!('1', Glu.factors, anorm)
    cond1 = 1/rcond
    solver = pnl.Backslash{Float64, typeof(Glu)}(G, Glu, zeros(N),
                                                 zeros(3, N), zeros(N))

    Dhat = Vinf/LA.norm(Vinf)
    Shat = [0.0, 1.0, 0.0]
    Lhat = LA.cross(Dhat, Shat)

    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false)
    force = pnl.ForceMonitor(1, 1; i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c_ref),
        correct_kuttacondition=false, verbose=false)

    pnl.steady!(body, pnl.ReferenceFrame(body), Vinf;
        body_solvers=solver,
        backend=pnl.FastMultipoleBackend(),
        monitors=(pressure, force),
        path=nothing, verbose=false)

    CL_solver = LA.dot(force.force[:, 1], Lhat)
    CD_solver = LA.dot(force.force[:, 1], Dhat)
    gamma = collect(view(body.strength, :, pnl.get_Gammai(body)))

    # --- B2.3a: algebraic residual ||G*gamma - rhs||inf (true G, when it fits)
    resid_inf, resid_rel = NaN, NaN
    if Gres !== nothing
        r = Gres*gamma .- solver.rhs
        resid_inf = LA.norm(r, Inf)
        resid_rel = resid_inf/max(LA.norm(solver.rhs, Inf), eps())
        Gres = nothing
    end
    solver = nothing; Glu = nothing; G = nothing
    GC.gc()

    # --- A: recover PV without any FMM call (matches the call at
    #     src/FLOWPanel_simulate.jl:441, validated bitwise earlier) ----------
    Jtri = pnl.compute_mu_gradient!(zeros(3, N), body.controlpoints,
        body.normals, body.cells, body.neighbor,
        view(body.strength, :, pnl.get_Gammai(body)),
        view(body.shedding_full, 1:2, :); scale=0.5, nodes=body.nodes)
    U0 = body.velocity .- Jtri                       # = U∞ + PV
    decomp_err = maximum(abs.((U0 .+ Jtri) .- body.velocity))

    geom = quad_geometry(body)
    grad_quad = quad_mu_gradient(geom, n_ch, n_span)
    grad_p1_area = p1_mu_gradient(body, n_ch, n_span)
    grad_p1_angle = p1_mu_gradient_angle(body, n_ch, n_span)
    Jtri_quad = quad_average(geom, Jtri)

    # Jump-sign convention (findings.md): body.velocity = U∞+PV+Jtri with
    # Jtri ≈ -½∇μ in the quad-gradient convention, so quad/P1 gradients use
    # jump_sign=-1 and the tri jump itself enters as 2*Jtri_quad with +1.
    res_tri  = quad_velocity_CL(geom, U0, 2 .* Jtri_quad, +1.0, Lhat, Dhat)
    res_quad = quad_velocity_CL(geom, U0, grad_quad, -1.0, Lhat, Dhat)
    res_p1a  = quad_velocity_CL(geom, U0, grad_p1_area, -1.0, Lhat, Dhat)
    res_p1g  = quad_velocity_CL(geom, U0, grad_p1_angle, -1.0, Lhat, Dhat)

    # Production auto-detecting quad-pitch reconstruction (basis=:quad):
    # should reproduce CL_quadgrad without knowing (n_ch, n_span). Jmd is the
    # per-triangle half-jump (-½∇μ_agg, same convention as Jtri); quad-average it
    # and feed as the jump exactly like res_tri.
    Jmd = pnl.compute_mu_gradient!(zeros(3, N), body.controlpoints, body.normals,
        body.cells, body.neighbor, view(body.strength, :, pnl.get_Gammai(body)),
        view(body.shedding_full, 1:2, :);
        scale=0.5, nodes=body.nodes, grad_mu_options=(; basis=:quad))
    res_mudiff = quad_velocity_CL(geom, U0, 2 .* quad_average(geom, Jmd), +1.0, Lhat, Dhat)

    # --- B2.1: Kutta-Joukowski CL from gamma alone (quad-averaged gamma of the
    #     TE rows; upper TE = chordwise strip j=1, lower = j=npts) -----------
    nid(k, j) = (k - 1)*npts + j
    CL_KJ = 0.0
    for k in 1:n_span
        Gamma_k = geom.mu[nid(k, 1)] - geom.mu[nid(k, npts)]
        dy_k = body.nodes[2, nid(k + 1, 1)] - body.nodes[2, nid(k, 1)]
        CL_KJ += Gamma_k*dy_k
    end
    # Sign flip: with this body's vortex-ring γ convention, positive lift gives
    # γ_upperTE − γ_lowerTE < 0 (verified at 40x48), so flip to report CL_KJ > 0.
    CL_KJ *= -rho*magVinf/(0.5*rho*magVinf^2*Sref)

    # --- TE-row sliver diagnostics (cells of chordwise strips 1/npts and
    #     2/npts-1) -----------------------------------------------------------
    strip(ci) = mod1((ci + 1) ÷ 2, npts)
    te_cells   = [ci for ci in 1:N if strip(ci) == 1 || strip(ci) == npts]
    row2_cells = [ci for ci in 1:N if strip(ci) == 2 || strip(ci) == npts - 1]
    dU0 = U0 .- Vinf                                  # pure PV
    Jtri_te_mean,  Jtri_te_max  = colnorm_stats(Jtri, te_cells)
    Jtri_r2_mean,  Jtri_r2_max  = colnorm_stats(Jtri, row2_cells)
    dU0_te_mean,   dU0_te_max   = colnorm_stats(dU0, te_cells)
    dU0_r2_mean,   dU0_r2_max   = colnorm_stats(dU0, row2_cells)

    # --- B2.3b: physical flow-tangency residual U.n/V∞ -----------------------
    tang = [LA.dot(view(body.velocity, :, i), view(body.normals, :, i))
            for i in 1:N] ./ magVinf
    tang_mean, tang_max = abs_stats(tang, 1:N)
    tang_te_mean, tang_te_max = abs_stats(tang, te_cells)

    h_min = c*(1 - cos(pi/n_ch))/2                    # TE sliver width

    row = (; variant=variant_name, n_ch, n_span, panels=N, h_min, core_size,
        CL_solver, CD_solver,
        CL_tri=res_tri.CL, CL_quadgrad=res_quad.CL, CL_mudiff=res_mudiff.CL,
        CL_p1_area=res_p1a.CL, CL_p1_angle=res_p1g.CL, CL_KJ,
        cond1, cond1_eps=cond1*eps(Float64), resid_inf, resid_rel,
        tang_mean, tang_max, tang_te_mean, tang_te_max,
        Jtri_te_mean, Jtri_te_max, Jtri_row2_mean=Jtri_r2_mean,
        Jtri_row2_max=Jtri_r2_max,
        dU0_te_mean, dU0_te_max, dU0_row2_mean=dU0_r2_mean,
        dU0_row2_max=dU0_r2_max,
        decomp_err, elapsed=time() - t0)
    Ufields = (; U_tri=res_tri.U, U_quad=res_quad.U,
               U_p1_area=res_p1a.U, U_p1_angle=res_p1g.U)
    return row, Ufields
end

# ----------------- DRIVER -------------------------------------------------------
function print_row(r)
    @printf("%7s %4d ko=%-7.0e CLslv %+.6f CLtri %+.6f CLquad %+.6f CLp1a %+.6f CLp1g %+.6f CLkj %+.6f cond1 %.2e res %.1e t %.0fs\n",
        r.variant, r.n_ch, r.core_size, r.CL_solver, r.CL_tri, r.CL_quadgrad,
        r.CL_p1_area, r.CL_p1_angle, r.CL_KJ, r.cond1, r.resid_inf, r.elapsed)
    flush(stdout)
end

function main_chorddiv()
    chord_levels = [parse(Int, s) for s in
        split(get(ENV, "FLOWPANEL_CHORDDIV_CASES", "40,56,80,112,156"), ",")]
    n_span = parse(Int, get(ENV, "FLOWPANEL_CHORDDIV_NSPAN", "48"))
    variants = [(name=String(s), swap=(s == "diagBD")) for s in
        split(get(ENV, "FLOWPANEL_CHORDDIV_VARIANTS", "diagAC,diagBD"), ",")]

    out_csv = joinpath("data", "sweptwing_sweeploft", "chord_divergence.csv")
    xvar_csv = joinpath("data", "sweptwing_sweeploft",
                        "chord_divergence_crossvariant.csv")
    mkpath(dirname(out_csv))

    println("Chord-divergence ladder: n_ch = ", chord_levels, " x n_span = ",
            n_span, ", variants = ", [v.name for v in variants])

    rows = []
    xrows = []
    for n_ch in chord_levels
        Uf = Dict{String, Any}()
        for v in variants
            r, Ufields = run_case(n_ch, n_span; swap_diagonals=v.swap,
                                  variant_name=v.name)
            print_row(r)
            push!(rows, r)
            CSV.write(out_csv, DataFrame(rows))
            Uf[v.name] = Ufields
            GC.gc()
        end
        if length(variants) == 2
            a, b = Uf[variants[1].name], Uf[variants[2].name]
            xr = (; n_ch, n_span,
                dU_tri=maxdiff(a.U_tri, b.U_tri)/magVinf,
                dU_quadgrad=maxdiff(a.U_quad, b.U_quad)/magVinf,
                dU_p1_area=maxdiff(a.U_p1_area, b.U_p1_area)/magVinf,
                dU_p1_angle=maxdiff(a.U_p1_angle, b.U_p1_angle)/magVinf)
            push!(xrows, xr)
            CSV.write(xvar_csv, DataFrame(xrows))
            @printf("        %4d cross-variant max|dU|/Vinf: tri %.3f quad %.3f p1area %.3f p1angle %.3f\n",
                n_ch, xr.dU_tri, xr.dU_quadgrad, xr.dU_p1_area, xr.dU_p1_angle)
            flush(stdout)
        end
    end

    # ----- summary tables ----------------------------------------------------
    println("\n#===== CL COLUMNS (divergence discriminator) =====#")
    @printf("%7s %5s %8s %10s %10s %10s %10s %10s %10s %10s\n", "variant", "n_ch",
        "panels", "CL_solver", "CL_tri", "CL_quad", "CL_mudiff", "CL_p1area", "CL_p1angl", "CL_KJ")
    for r in rows
        @printf("%7s %5d %8d %+10.6f %+10.6f %+10.6f %+10.6f %+10.6f %+10.6f %+10.6f\n",
            r.variant, r.n_ch, r.panels, r.CL_solver, r.CL_tri, r.CL_quadgrad,
            r.CL_mudiff, r.CL_p1_area, r.CL_p1_angle, r.CL_KJ)
    end

    println("\n#===== SOLVE QUALITY (B2) =====#")
    @printf("%7s %5s %10s %10s %10s %10s %10s %11s %11s\n", "variant", "n_ch",
        "cond1", "cond1*eps", "resid_inf", "resid_rel", "tang_mean", "tang_max",
        "tang_te_max")
    for r in rows
        @printf("%7s %5d %10.3e %10.3e %10.3e %10.3e %10.3e %11.3e %11.3e\n",
            r.variant, r.n_ch, r.cond1, r.cond1_eps, r.resid_inf, r.resid_rel,
            r.tang_mean, r.tang_max, r.tang_te_max)
    end

    println("\n#===== TE-ROW SLIVER SCALING (x h_min, x h_min^2) =====#")
    @printf("%7s %5s %9s | %9s %9s %9s | %9s %9s %9s\n", "variant", "n_ch",
        "h_min", "Jte_max", "*h", "*h^2", "dU0te_max", "*h", "*h^2")
    for r in rows
        @printf("%7s %5d %9.3e | %9.3e %9.3e %9.3e | %9.3e %9.3e %9.3e\n",
            r.variant, r.n_ch, r.h_min,
            r.Jtri_te_max, r.Jtri_te_max*r.h_min, r.Jtri_te_max*r.h_min^2,
            r.dU0_te_max, r.dU0_te_max*r.h_min, r.dU0_te_max*r.h_min^2)
    end

    println("\nCSV written to ", out_csv, " and ", xvar_csv)
    return rows, xrows
end

# ----------------- SECTION C: CORE_SIZE PROBE (run only if A implicates PV)
function core_size_probe()
    offsets = [parse(Float64, s) for s in
        split(ENV["FLOWPANEL_CHORDDIV_KERNELOFFSETS"], ",") if !isempty(strip(s))]
    ko_cases = [parse(Int, s) for s in
        split(get(ENV, "FLOWPANEL_CHORDDIV_KO_CASES", "80,112"), ",")]
    n_span = parse(Int, get(ENV, "FLOWPANEL_CHORDDIV_NSPAN", "48"))
    out_csv = joinpath("data", "sweptwing_sweeploft",
                       "chord_divergence_core_size.csv")

    println("\nCore_size probe (diagAC): n_ch = ", ko_cases,
            ", core_size = ", offsets)
    rows = []
    for n_ch in ko_cases, ko in offsets
        r, _ = run_case(n_ch, n_span; swap_diagonals=false, core_size=ko,
                        variant_name="diagAC")
        h_min = c*(1 - cos(pi/n_ch))/2
        @printf("  n_ch %4d ko %.0e (ko/h_min %.1e): CLslv %+.6f CLtri %+.6f dU0te_max %.3e\n",
            n_ch, ko, ko/h_min, r.CL_solver, r.CL_tri, r.dU0_te_max)
        flush(stdout)
        push!(rows, r)
        CSV.write(out_csv, DataFrame(rows))
        GC.gc()
    end
    return rows
end

# ----------------- GRAD_MU_OPTIONS CONVERGENCE SWEEP ---------------------------
# Six grad_mu surface-velocity reconstructions, swept across a chord ladder at a
# cheap n_span to see which keep CL bounded (converge) and which blow up as the
# TE slivers sharpen (diverge). Each config feeds steady! directly via
# grad_mu_options, so CL is the production ForceMonitor value for body.velocity.
const GRADMU_CONFIGS = [
    ("tri",          (; basis=:tri)),
    ("tri_robust",   (; basis=:tri, tri_robust=true)),
    ("quad_nogrow",  (; basis=:quad, quad_grow=false)),
    ("quad_default", (; basis=:quad)),
    ("quad_depth2",  (; basis=:quad, quad_grow=true, quad_grow_stop=:depth,
                         quad_grow_max_depth=2)),
    ("quad_cond12",  (; basis=:quad, quad_grow_cond_max=12.0)),
]

"""
    run_gradmu_case(n_ch, n_span; swap_diagonals, variant_name) -> Vector{NamedTuple}

Solve one (n_ch, n_span, variant) sweep-lofted wing ONCE (Backslash, LU factored
a single time), then for each grad_mu_options config in `GRADMU_CONFIGS` re-run the
production steady! velocity reconstruction + force integration and record CL/CD.
gamma is identical across configs (same prefactored system); only the half-jump
reconstruction and force integration repeat.
"""
function run_gradmu_case(n_ch::Int, n_span::Int; swap_diagonals::Bool,
                         variant_name::String="")
    body = sweeploft_wing(n_ch, n_span; mirror_diagonals=true, swap_diagonals,
                          core_size=1e-10)
    N = body.ncells
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    solver = pnl.Backslash(body)               # LU factored once, reused below

    Dhat = Vinf/LA.norm(Vinf)
    Shat = [0.0, 1.0, 0.0]
    Lhat = LA.cross(Dhat, Shat)

    rows = NamedTuple[]
    for (name, opts) in GRADMU_CONFIGS
        t0 = time()
        pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false)
        force = pnl.ForceMonitor(1, 1; i_frame=-1,
            normalization=pnl.WingNormalization(rho, Sref, c_ref),
            correct_kuttacondition=false, verbose=false)
        pnl.steady!(body, pnl.ReferenceFrame(body), Vinf;
            body_solvers=solver,
            backend=pnl.FastMultipoleBackend(),
            monitors=(pressure, force),
            grad_mu_options=opts,
            path=nothing, verbose=false)
        CL = LA.dot(force.force[:, 1], Lhat)
        CD = LA.dot(force.force[:, 1], Dhat)
        push!(rows, (; variant=variant_name, n_ch, n_span, panels=N,
                       config=name, CL, CD, elapsed=time() - t0))
        @printf("  %7s n_ch %4d %-12s CL %+.6f CD %+.6f  t %.0fs\n",
            variant_name, n_ch, name, CL, CD, time() - t0)
        flush(stdout)
    end
    GC.gc()
    return rows
end

function main_gradmu_sweep()
    chord_levels = [parse(Int, s) for s in
        split(get(ENV, "FLOWPANEL_CHORDDIV_GRADMU_CASES",
                  "16,24,40,56,80,112,156"), ",")]
    n_span = parse(Int, get(ENV, "FLOWPANEL_CHORDDIV_NSPAN", "12"))
    variants = [(name="diagAC", swap=false), (name="diagBD", swap=true)]
    blowup_factor = parse(Float64,
        get(ENV, "FLOWPANEL_CHORDDIV_GRADMU_BLOWUP", "10.0"))
    abs_floor = 5.0

    out_csv = joinpath("data", "sweptwing_sweeploft",
                       "chord_divergence_gradmu.csv")
    mkpath(dirname(out_csv))

    println("grad_mu_options sweep: n_ch = ", chord_levels, " x n_span = ", n_span,
            ", variants = ", [v.name for v in variants],
            ", configs = ", [c[1] for c in GRADMU_CONFIGS])

    rows = NamedTuple[]
    reference_CL = NaN            # set from first n_ch level (median |CL|)
    for n_ch in chord_levels
        level_CL = Float64[]
        for v in variants
            crows = run_gradmu_case(n_ch, n_span; swap_diagonals=v.swap,
                                    variant_name=v.name)
            append!(rows, crows)
            append!(level_CL, [abs(r.CL) for r in crows])
            CSV.write(out_csv, DataFrame(rows))
        end
        if isnan(reference_CL)
            s = sort(level_CL)
            n = length(s)
            reference_CL = isodd(n) ? s[(n+1)÷2] : (s[n÷2] + s[n÷2+1])/2
        end
        worst = maximum(level_CL)
        thresh = max(blowup_factor*reference_CL, abs_floor)
        if worst > thresh
            @printf("\nCL blowing up at n_ch=%d (max|CL| %.3f > %.3f); ending ladder early.\n",
                n_ch, worst, thresh)
            flush(stdout)
            break
        end
    end

    # ----- pivot tables: rows = config, cols = n_ch, cells = CL -----
    ran_levels = sort(unique(r.n_ch for r in rows))
    for v in variants
        println("\n#===== CL(n_ch) | variant = ", v.name,
                " (rows: config) =====#")
        @printf("%-13s", "config")
        for nc in ran_levels
            @printf(" %10d", nc)
        end
        println()
        for (name, _) in GRADMU_CONFIGS
            @printf("%-13s", name)
            for nc in ran_levels
                idx = findfirst(r -> r.variant == v.name && r.n_ch == nc &&
                                     r.config == name, rows)
                if idx === nothing
                    @printf(" %10s", "-")
                else
                    @printf(" %+10.6f", rows[idx].CL)
                end
            end
            println()
        end
    end

    println("\nCSV written to ", out_csv)
    return rows
end

if !isinteractive() && get(ENV, "FLOWPANEL_CHORDDIV_GRADMU", "false") == "true"
    gradmu_rows = main_gradmu_sweep()
elseif !isinteractive() && get(ENV, "FLOWPANEL_CHORDDIV_RUN", "true") == "true"
    rows, xrows = main_chorddiv()
    if !isempty(get(ENV, "FLOWPANEL_CHORDDIV_KERNELOFFSETS", ""))
        ko_rows = core_size_probe()
    end
end
