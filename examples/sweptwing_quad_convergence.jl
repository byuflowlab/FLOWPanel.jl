# Grid-convergence study of CL for the sweep-lofted Weber swept wing using the
# quad surface-velocity reconstruction grad_mu_options = (; basis=:quad,
# quad_grow=true, quad_grow_stop=:depth, quad_grow_max_depth=2).
#
# Two stages:
#   Stage 1 (span):  hold n_ch, sweep n_span -> pick converged spanwise resolution
#   Stage 2 (chord): hold n_span, sweep n_ch  (run later)
#
# Run with 3 threads:
#   julia --project -t3 examples/sweptwing_quad_convergence.jl
#
# Env knobs:
#   FLOWPANEL_QUADCONV_STAGE   = span | chord            (default span)
#   FLOWPANEL_QUADCONV_NCH     = 48                       (Stage 1 fixed chord)
#   FLOWPANEL_QUADCONV_NSPANS  = 16,24,32,48,64,96        (Stage 1 ladder)
#   FLOWPANEL_QUADCONV_NSPAN_FIXED = <int>                (Stage 2 fixed span; required)
#   FLOWPANEL_QUADCONV_NCHS    = 24,32,48,64,96,128       (Stage 2 ladder)
#   FLOWPANEL_QUADCONV_TOL     = 5e-4                      (|dCL| convergence flag)

ENV["FLOWPANEL_SWEEPLOFT_RUN"] = "false"
include(joinpath(@__DIR__, "sweptwing_sweeploft.jl"))

using GeometricTools: PyPlot as plt

# Fixed quad reconstruction (quad_grow_max_depth may be raised later).
const GRAD_MU = (; basis=:quad, quad_grow=true, quad_grow_stop=:depth,
                   quad_grow_max_depth=2)

const QUADCONV_DIR = joinpath("data", "sweptwing_quad_convergence")

# ----------------- ONE CASE ---------------------------------------------------
"""
    solve_cl(n_ch, n_span; grad_mu_options, backslash_max_panels, path, name)

Solve one diagAC sweep-lofted wing with the production steady! path and the given
grad_mu_options, returning a NamedTuple with CL/CD (Bernoulli + WingNormalization).
Mirrors solve_case in sweptwing_sweeploft.jl but threads grad_mu_options through and
defaults to a lower Backslash cutoff so the larger cases use Krylov+FMM.
"""
function solve_cl(n_ch::Int, n_span::Int; grad_mu_options=GRAD_MU,
                  backslash_max_panels::Int=25000, path=nothing, name="quadconv")
    body = sweeploft_wing(n_ch, n_span; mirror_diagonals=true, swap_diagonals=false)

    solver = body.ncells <= backslash_max_panels ?
        pnl.Backslash(body) :
        pnl.KrylovSolver(body; backend=pnl.FastMultipoleBackend(),
                         method=:gmres, atol=1e-9, rtol=1e-9, itmax=1000)
    solver_name = body.ncells <= backslash_max_panels ? "backslash" : "krylov"

    Dhat = Vinf/LA.norm(Vinf)
    Shat = [0.0, 1.0, 0.0]
    Lhat = LA.cross(Dhat, Shat)

    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false)
    force = pnl.ForceMonitor(1, 1; i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c_ref),
        correct_kuttacondition=false, verbose=false)

    elapsed = @elapsed pnl.steady!(body, pnl.ReferenceFrame(body), Vinf;
        body_solvers=solver,
        backend=pnl.FastMultipoleBackend(),
        monitors=(pressure, force),
        grad_mu_options=grad_mu_options,
        path=path, name=name, verbose=false)

    CL = LA.dot(force.force[:, 1], Lhat)
    CD = LA.dot(force.force[:, 1], Dhat)
    return (; n_ch, n_span, panels=body.ncells, solver=solver_name,
              CL, CD, CL_error_pct=100*(CL - CLexp)/CLexp, elapsed)
end

# ----------------- PLOT -------------------------------------------------------
function plot_convergence(rows, xsym::Symbol, xlabel::AbstractString, stage::String;
                          xtransform=identity)
    x  = [xtransform(getproperty(r, xsym)) for r in rows]
    CL = [r.CL for r in rows]
    CD = [r.CD for r in rows]

    fig = plt.figure("quadconv_$(stage)", figsize=(11, 4.5))
    plt.clf()
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.plot(x, CL; marker="o", color="C0")
    ax1.axhline(CLexp; color="k", ls="--", lw=1, label="CLexp = $(CLexp)")
    ax1.set_xlabel(xlabel); ax1.set_ylabel("C_L")
    ax1.set_title("C_L convergence ($stage)")
    ax1.grid(true, alpha=0.3); ax1.legend(loc="best", fontsize=9)

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.plot(x, CD; marker="o", color="C3")
    ax2.set_xlabel(xlabel); ax2.set_ylabel("C_D")
    ax2.set_title("C_D convergence ($stage)")
    ax2.grid(true, alpha=0.3)

    plt.tight_layout()
    out = joinpath(QUADCONV_DIR, "$(stage)_convergence.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    println("Saved ", out)
end

# ----------------- DRIVERS ----------------------------------------------------
function _run_ladder(cases::Vector{Tuple{Int,Int}}, stage::String, xsym::Symbol,
                     xlabel::AbstractString; xtransform=identity)
    cl_tol = parse(Float64, get(ENV, "FLOWPANEL_QUADCONV_TOL", "5e-4"))
    out_csv = joinpath(QUADCONV_DIR, "$(stage)_convergence.csv")
    mkpath(QUADCONV_DIR)

    println("quad convergence ($stage): cases (n_ch,n_span) = ", cases,
            ", grad_mu = ", GRAD_MU)
    @printf("%5s %5s %8s %10s %+11s %10s %10s %8s\n",
        "n_ch", "nspan", "panels", "solver", "CL", "dCL", "err%", "t[s]")

    rows = NamedTuple[]
    CL_prev = NaN
    for (n_ch, n_span) in cases
        casetag = "$(stage)_nch$(n_ch)_ns$(n_span)"
        path = joinpath(QUADCONV_DIR, casetag)
        r = solve_cl(n_ch, n_span; path=path, name=casetag)
        push!(rows, r)
        CSV.write(out_csv, DataFrame(rows))

        dCL = isnan(CL_prev) ? NaN : r.CL - CL_prev
        flag = (!isnan(dCL) && abs(dCL) < cl_tol) ? "  <-- converged" : ""
        @printf("%5d %5d %8d %10s %+11.6f %+10.6f %+9.3f %8.2f%s\n",
            r.n_ch, r.n_span, r.panels, r.solver, r.CL, dCL, r.CL_error_pct,
            r.elapsed, flag)
        flush(stdout)
        CL_prev = r.CL
    end

    plot_convergence(rows, xsym, xlabel, stage; xtransform)
    println("\nCSV written to ", out_csv)
    return rows
end

function converge_span()
    n_ch = parse(Int, get(ENV, "FLOWPANEL_QUADCONV_NCH", "48"))
    n_spans = [parse(Int, s) for s in
        split(get(ENV, "FLOWPANEL_QUADCONV_NSPANS", "16,24,32,48,64,96"), ",")]
    cases = [(n_ch, ns) for ns in n_spans]
    return _run_ladder(cases, "span", :n_span, "n_span (spanwise panels)")
end

function converge_chord()
    haskey(ENV, "FLOWPANEL_QUADCONV_NSPAN_FIXED") ||
        error("Stage chord requires FLOWPANEL_QUADCONV_NSPAN_FIXED=<int>")
    n_span = parse(Int, ENV["FLOWPANEL_QUADCONV_NSPAN_FIXED"])
    n_chs = [parse(Int, s) for s in
        split(get(ENV, "FLOWPANEL_QUADCONV_NCHS", "24,32,48,64,96"), ",")]
    cases = [(nc, n_span) for nc in n_chs]
    return _run_ladder(cases, "chord", :n_ch, "n_ch (chordwise panels)")
end

# ----------------- LOADING PLOTS (one discretization) -------------------------
"""
    loading_plots(n_ch, n_span)

Solve one diagAC case with the quad reconstruction, then draw (a) the spanwise
sectional lift/drag loading from a SpanwiseLoadingMonitor and (b) chordwise Cp at a
few spanwise stations, extracted from quad_geometry + body.P. Saves two PNGs under
data/sweptwing_quad_convergence/.
"""
function loading_plots(n_ch::Int, n_span::Int; backslash_max_panels::Int=25000)
    body = sweeploft_wing(n_ch, n_span; mirror_diagonals=true, swap_diagonals=false)
    solver = body.ncells <= backslash_max_panels ?
        pnl.Backslash(body) :
        pnl.KrylovSolver(body; backend=pnl.FastMultipoleBackend(),
                         method=:gmres, atol=1e-9, rtol=1e-9, itmax=1000)

    Dhat = Vinf/LA.norm(Vinf)
    Shat = [0.0, 1.0, 0.0]
    Lhat = LA.cross(Dhat, Shat)

    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false)
    force = pnl.ForceMonitor(1, 1; i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c_ref),
        correct_kuttacondition=false, verbose=false)
    spanwise = pnl.SpanwiseLoadingMonitor(n_span, 1;
        components=(lift=Lhat, drag=Dhat), span_axis=Shat, per_length=true,
        normalization=pnl.NoSectionalNormalization(), file=false)

    pnl.steady!(body, pnl.ReferenceFrame(body), Vinf;
        body_solvers=solver, backend=pnl.FastMultipoleBackend(),
        monitors=(pressure, force, spanwise),
        grad_mu_options=GRAD_MU, path=nothing, verbose=false)

    qinf = 0.5*rho*magVinf^2
    mkpath(QUADCONV_DIR)

    # --- (a) spanwise sectional loading -------------------------------------
    y = spanwise.bin_center
    eta = y ./ (b/2)
    cl = spanwise.load_components[1, :] ./ (qinf*c)   # sectional cl (untapered c)
    cd = spanwise.load_components[2, :] ./ (qinf*c)

    fig = plt.figure("quadconv_spanloading", figsize=(7.5, 4.5)); plt.clf()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(eta, cl; marker="o", ms=3, color="C0", label="sectional c_l")
    ax.plot(eta, cd; marker="s", ms=3, color="C3", label="sectional c_d")
    ax.axhline(0.0; color="k", lw=0.6, alpha=0.5)
    ax.set_xlabel("2y/b"); ax.set_ylabel("sectional coefficient")
    ax.set_title("Spanwise loading  ($(n_ch)x$(n_span), quad_depth2)")
    ax.grid(true, alpha=0.3); ax.legend(loc="best")
    plt.tight_layout()
    out1 = joinpath(QUADCONV_DIR, "loading_spanwise.png")
    plt.savefig(out1, dpi=150, bbox_inches="tight"); println("Saved ", out1)

    # --- (b) chordwise Cp at a few stations ---------------------------------
    geom = quad_geometry(body)
    npts = 2*n_ch
    P = pressure.pressure[1]                            # per-cell pressure (Pa)
    Cp_quad = [(geom.areas[2q-1]*P[2q-1] + geom.areas[2q]*P[2q]) /
               geom.A[q] / qinf for q in 1:(body.ncells ÷ 2)]
    strip_y = [sum(geom.centroid[2, (k-1)*npts+1 : k*npts]) / npts for k in 1:n_span]

    fig2 = plt.figure("quadconv_chordloading", figsize=(7.5, 4.5)); plt.clf()
    ax2 = fig2.add_subplot(1, 1, 1)
    for (ci, target) in enumerate((0.2, 0.5, 0.8))     # 2y/b stations (+y side)
        k = argmin(abs.(strip_y ./ (b/2) .- target))
        cols = (k-1)*npts+1 : k*npts
        xs = geom.centroid[1, cols]
        xLE, xTE = minimum(xs), maximum(xs)
        xc = (xs .- xLE) ./ (xTE - xLE)
        upper = 1:n_ch                                  # cols 1..n_ch one surface
        lower = n_ch+1:npts
        col = "C$(ci-1)"
        ax2.plot(xc[upper], [Cp_quad[c] for c in cols[upper]];
                 color=col, lw=1.2, label="2y/b≈$(round(strip_y[k]/(b/2), digits=2))")
        ax2.plot(xc[lower], [Cp_quad[c] for c in cols[lower]];
                 color=col, lw=1.2, ls="--")
    end
    ax2.invert_yaxis()                                  # suction (−Cp) up
    ax2.set_xlabel("x/c"); ax2.set_ylabel("C_p")
    ax2.set_title("Chordwise Cp  ($(n_ch)x$(n_span); solid=upper, dashed=lower)")
    ax2.grid(true, alpha=0.3); ax2.legend(loc="best", fontsize=8)
    plt.tight_layout()
    out2 = joinpath(QUADCONV_DIR, "loading_chordwise.png")
    plt.savefig(out2, dpi=150, bbox_inches="tight"); println("Saved ", out2)

    return body, (; spanwise, pressure, force)
end

# Joint refinement: scale n_ch and n_span together so chordwise and spanwise
# truncation errors (which have OPPOSITE signs here) shrink simultaneously,
# revealing the genuine continuum CL. x-axis is sqrt(panels) ~ 1/h.
function converge_joint()
    cases = parse_case_list(get(ENV, "FLOWPANEL_QUADCONV_JOINT",
        "24:16,32:24,48:32,64:48,96:64"))
    return _run_ladder(cases, "joint", :panels, "sqrt(panels)  (~ 1/h)";
                       xtransform=sqrt)
end

if !isinteractive()
    stage = get(ENV, "FLOWPANEL_QUADCONV_STAGE", "span")
    rows = stage == "span"  ? converge_span()  :
           stage == "chord" ? converge_chord() :
           stage == "joint" ? converge_joint() :
           stage == "loading" ? loading_plots(
               parse(Int, get(ENV, "FLOWPANEL_QUADCONV_NCH", "96")),
               parse(Int, get(ENV, "FLOWPANEL_QUADCONV_NSPAN_FIXED", "64"))) :
           error("FLOWPANEL_QUADCONV_STAGE must be span, chord, joint, or loading, got $(stage)")
end
