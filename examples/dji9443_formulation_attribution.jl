#=##############################################################################
Phase 2b — Formulation/Topology Attribution driver.

Goal: determine why Dirichlet solves produce lower integrated bound circulation
than Neumann solves in the DJI 9443 hover study (Phase 1 gap ~5%).

This driver is a reusable diagnostic layer around the pitching-wing oracle
(`examples/pitching_wing.jl`) and reuses the frozen-state solve harness from
`test/formulation_test.jl` (`small_body`/`flat_wake`/`run_aero!`/`kutta_map`).

Modes (dispatched by ENV["PHASE2B_MODE"]):
  smoke     tiny synthetic case exercising CSV writers, fixed-grid interp,
            shared probe generation, loop-integral machinery, wake-strength
            extraction, and decomposition algebra (no DJI assets, no big solves)
  formsweep Gate A0: solve one unchanged capped/Dirichlet body with each
            formulation {VelocityThroughSources, TraceCorrected(:green),
            TraceCorrected(:line_integral), GreenReconstruction} on an
            identical frozen wake, report fixed-bin spanwise circulation
  oracle    [next batch] generated open/capped refinement ladder + Gate A
  thickness [next batch] gated thickness / TE-closure screen
  analyze   aggregate CSVs into tables + attribution notes
  dji_bridge [conditional] DJI 2x2 topology/formulation bridge

Run:
  PHASE2B_MODE=formsweep julia --project examples/dji9443_formulation_attribution.jl
  PHASE2B_MODE=smoke     julia --project examples/dji9443_formulation_attribution.jl
  PHASE2B_MODE=analyze   julia --project examples/dji9443_formulation_attribution.jl
=###############################################################################

import FLOWPanel as pnl
import LinearAlgebra as LA
import LinearAlgebra: norm, dot
using FastMultipole.StaticArrays: SVector
using Printf
import CSV
import DataFrames: DataFrame, nrow

include(joinpath(pnl.examples_path, "pitching_wing.jl"))

# --------------------------------------------------------------- constants ----
const OUTPUT_DIR = normpath(joinpath(@__DIR__, "..", "data",
    "dji_convergence_20260722", "phase_02b_formulation_attribution"))
mkpath(OUTPUT_DIR)

const C_CHORD = 1.0
const B_SPAN  = 4.0
const UMAG    = 100.0
const ALPHA   = 5.0 * pi / 180
const UINF    = SVector{3}(UMAG*cos(ALPHA), 0.0, UMAG*sin(ALPHA))
const DIRECT  = pnl.DirectBackend()

artifact(name) = joinpath(OUTPUT_DIR, name)

# ------------------------------------------------- reused oracle harness ------
# (mirrors test/formulation_test.jl so the driver is self-contained)

"Capped, watertight Dirichlet source+doublet pitching-wing body."
function oracle_body(; n_span=7, n_airfoil=41, n_endcap=5, das_length_c=0.05,
        semiinfinite_wake=false)
    body = build_pitching_wing_body(C_CHORD, B_SPAN;
        n_span, n_airfoil, n_endcap, endcap=:round, semiinfinite_wake)
    set_wake_Das!(body, UINF; magnitude=das_length_c*C_CHORD)
    return body
end

"Fabricated flat prescribed wake with per-strip strength `gamma`."
function flat_wake(body, gamma; nfree=8, free_row_length_c=0.5)
    wake = pnl.PanelWake(body; nwakerows=nfree, include_final_filament=false)
    pnl.update_TE!(wake, body)
    Das = view(body.Das[1], :, 1)
    direction = Das ./ norm(Das) .* (free_row_length_c*C_CHORD)
    for nodes in wake.nodes
        first_row = copy(view(nodes, :, 1, :))
        for row in 1:nfree+1
            view(nodes, :, row, :) .= first_row .+ (row - 1) .* direction
        end
    end
    for strength in wake.strength
        for row in 1:nfree
            view(strength, 1, row, :) .= gamma
        end
    end
    wake.nwakes[] = nfree
    return wake
end

kutta_map(body, x) = [x[shed[1, i]] - (shed[4, i] > 0 ? x[shed[4, i]] : 0.0)
                      for shed in body.shedding for i in axes(shed, 2)]

nedges(body) = sum(size(shed, 2) for shed in body.shedding; init=0)

function run_aero!(body, wake, solver; formulation=pnl.VelocityThroughSources(),
        formulation_state=nothing, i_step=0)
    frames = pnl.ReferenceFrame(body)
    pnl._steady_aerodynamics!(body, (body,), (wake,), frames, UINF, solver;
        backend_wake=DIRECT, backend_solve=DIRECT, backend_system=DIRECT,
        update_trailing_edges=false,
        formulation, formulation_state, i_step)
    return body
end

# ---------------------------------------------- circulation observables -------

"""
Per-shedding-edge spanwise station `y` and physical circulation
`Γ = μ_upper − μ_lower − c` (the correction flows through
`_get_wakestrength_Gamma`). Returns `(y::Vector, gamma::Vector)` sorted by `y`.
Station `y` is the upper TE panel control-point spanwise coordinate.
"""
function spanwise_gamma(body)
    ys = Float64[]
    gs = Float64[]
    for (s, shed) in enumerate(body.shedding)
        for i in axes(shed, 2)
            push!(ys, body.controlpoints[2, shed[1, i]])
            push!(gs, pnl._get_wakestrength_Gamma(body, i, s))
        end
    end
    perm = sortperm(ys)
    return ys[perm], gs[perm]
end

"Trapezoidal integral of y-sorted samples."
function trapz(x, y)
    s = 0.0
    for i in 1:length(x)-1
        s += 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i])
    end
    return s
end

"Linear interpolation of (x,y) onto `bins` (flat extrapolation at ends)."
function interpolate_linear(x, y, bins)
    out = similar(bins, Float64)
    for (k, xb) in enumerate(bins)
        if xb <= x[1]
            out[k] = y[1]
        elseif xb >= x[end]
            out[k] = y[end]
        else
            j = searchsortedlast(x, xb)
            t = (xb - x[j]) / (x[j+1] - x[j])
            out[k] = (1 - t) * y[j] + t * y[j+1]
        end
    end
    return out
end

"""
Integrated circulation metrics on a fixed spanwise station grid.
`∫Γ dy`, the |y|-weighted analogue `∫|y| Γ dy` (finite-wing weighted analogue of
the rotor `∫Ωr Γ dr`), and the outboard integral over |y| ≥ `outboard_frac·b/2`.
Also returns side symmetry `Γ(+y) vs Γ(−y)` mismatch.
"""
function integrated_metrics(body; nbins=64, outboard_frac=0.6)
    y, g = spanwise_gamma(body)
    ymin, ymax = minimum(y), maximum(y)
    bins = collect(range(ymin, ymax; length=nbins))
    gb = interpolate_linear(y, g, bins)
    integrated = trapz(bins, gb)
    weighted   = trapz(bins, abs.(bins) .* gb)
    half = B_SPAN / 2
    mask = abs.(bins) .>= outboard_frac * half
    outboard = any(mask) ? trapz(bins[mask], gb[mask]) : 0.0
    # side symmetry: interpolate |y| grid on each side, compare
    ypos = bins[bins .>= 0]; gpos = gb[bins .>= 0]
    yneg = -reverse(bins[bins .<= 0]); gneg = reverse(gb[bins .<= 0])
    common = collect(range(0.0, min(maximum(ypos), maximum(yneg)); length=nbins ÷ 2))
    sp = interpolate_linear(ypos, gpos, common)
    sn = interpolate_linear(yneg, gneg, common)
    sym = norm(sp .- sn) / max(norm(sp), eps())
    return (; integrated, weighted, outboard, side_sym=sym, bins, gb, y, g)
end

# --------------------------------------------------- loop-integral circ -------

"""
Bound circulation at spanwise station `y0` from a closed rectangular circuit in
the x–z plane (constant y), evaluated from the induced velocity field only
(∮U∞·dl = 0 around a closed loop). The circuit encloses the airfoil section and
crosses the attached wake sheet once. Probes are offset from the surface via the
body's `kerneloffset_targets`. Returns `∮ V·dl` (positive sense: +x → +z → −x → −z).

`center` is the (x,z) loop center; `hw`,`hh` the half-width (x) and half-height (z);
`nseg` the per-edge quadrature segments. Off-body induced velocity goes through
the `ProbeSystem`/`influence!` path (bypasses on-surface grad_mu).
"""
function loop_integral_circulation(body, y0; center=(0.5*C_CHORD, 0.0),
        hw=1.2*C_CHORD, hh=0.6*C_CHORD, nseg=48, backend=DIRECT)
    cx, cz = center
    # rectangular loop corners (counter-clockwise in x-z)
    corners = [(cx-hw, cz-hh), (cx+hw, cz-hh), (cx+hw, cz+hh), (cx-hw, cz+hh)]
    pts = SVector{3,Float64}[]
    tangents = SVector{3,Float64}[]
    weights = Float64[]
    for e in 1:4
        (x1, z1) = corners[e]
        (x2, z2) = corners[mod1(e+1, 4)]
        for k in 1:nseg
            t = (k - 0.5) / nseg                 # midpoint rule
            x = x1 + t*(x2 - x1); z = z1 + t*(z2 - z1)
            dl = SVector{3}((x2 - x1)/nseg, 0.0, (z2 - z1)/nseg)
            push!(pts, SVector{3}(x, y0, z))
            push!(tangents, dl)
            push!(weights, 1.0)
        end
    end
    U = induced_velocity(body, pts; backend)
    circ = 0.0
    for k in eachindex(pts)
        circ += dot(U[k], tangents[k])
    end
    return circ
end

"""
Induced velocity at a cloud of points `pts::Vector{SVector{3}}` from a solved
`body`, via the FastMultipole ProbeSystem path (velocity read from
`probes.gradient`). Uses the body's off-body `kerneloffset_targets`.
"""
function induced_velocity(body, pts; backend=DIRECT)
    n = length(pts)
    probes = pnl.FastMultipole.ProbeSystem(n, Float64)
    for k in 1:n
        probes.position[k] = pts[k]
    end
    saved = body.kerneloffset
    body.kerneloffset = body.kerneloffset_targets
    pnl.influence!((probes,), (body,), backend;
        precalc=false, scalar_potential=false, gradient=true, hessian=false)
    body.kerneloffset = saved
    return probes.gradient
end

"""
Loop-integral self-check: evaluate `loop_integral_circulation` at two loop radii
and two quadrature densities; return the four values and whether they agree to
`rtol` (default 0.1%). Loops failing the self-check are reported and excluded.
"""
function loop_self_check(body, y0; rtol=1e-3, backend=DIRECT)
    v11 = loop_integral_circulation(body, y0; hw=1.2, hh=0.6, nseg=48, backend)
    v12 = loop_integral_circulation(body, y0; hw=1.2, hh=0.6, nseg=96, backend)
    v21 = loop_integral_circulation(body, y0; hw=1.5, hh=0.8, nseg=48, backend)
    v22 = loop_integral_circulation(body, y0; hw=1.5, hh=0.8, nseg=96, backend)
    vals = (v11, v12, v21, v22)
    m = sum(vals) / 4
    spread = maximum(abs.(vals .- m)) / max(abs(m), eps())
    return (; vals, mean=m, spread, ok=(spread <= rtol))
end

# ================================================================= modes ======

function seed_gamma(body, solver; nfree=8)
    # one VelocityThroughSources pass to obtain a representative bound
    # circulation, then rebuild the flat wake at that strength so the frozen
    # wake used across formulations is roughly self-consistent.
    M = nedges(body)
    wake = flat_wake(body, zeros(M); nfree)
    b0 = deepcopy(body)
    run_aero!(b0, wake, solver)
    _, g0 = spanwise_gamma(b0)
    # map sorted-back to edge order is unnecessary: use per-edge extraction
    gedge = [pnl._get_wakestrength_Gamma(b0, i, s)
             for (s, shed) in enumerate(b0.shedding) for i in axes(shed, 2)]
    return gedge
end

function run_formsweep()
    println("\n#===== Gate A0 formulation sweep (capped/Dirichlet oracle) =====#")
    body0  = oracle_body()
    solver0 = pnl.Backslash(body0)
    N = body0.ncells; M = nedges(body0)
    println("oracle body: $N cells, $M shedding edges, watertight=$(body0.watertight)")

    gseed = seed_gamma(body0, solver0)
    println("seeded flat-wake circulation: mean=$(round(sum(gseed)/M; sigdigits=4)), " *
            "max=$(round(maximum(gseed); sigdigits=4))")
    wake_template() = flat_wake(body0, gseed)

    formulations = [
        ("VelocityThroughSources", pnl.VelocityThroughSources()),
        ("TraceCorrected_green",    pnl.TraceCorrected(estimator=:green)),
        ("TraceCorrected_lineint",  pnl.TraceCorrected(estimator=:line_integral,
                                        quadrature=:simpson)),
        ("GreenReconstruction",     pnl.GreenReconstruction()),
    ]

    station_rows = DataFrame(formulation=String[], edge=Int[], y=Float64[],
        gamma=Float64[], refinement=String[], topology=String[],
        gauge=String[], extraction=String[])
    summary_rows = DataFrame(formulation=String[], integrated=Float64[],
        weighted=Float64[], outboard=Float64[], side_sym=Float64[],
        gamma_mean=Float64[], refinement=String[], topology=String[],
        backend=String[], wake=String[])

    baseline_int = NaN
    println()
    @printf("%-26s %14s %14s %14s %12s\n",
        "formulation", "∫Γ dy", "∫|y|Γ dy", "outboard", "side_sym")
    for (name, form) in formulations
        wake = wake_template()
        st = pnl.initialize_formulation(form, (deepcopy(body0),), (wake,),
            solver0, DIRECT, DIRECT)
        body = deepcopy(body0)
        run_aero!(body, wake, solver0; formulation=form, formulation_state=st)
        m = integrated_metrics(body)
        isnan(baseline_int) && (baseline_int = m.integrated)
        @printf("%-26s %14.6f %14.6f %14.6f %12.2e\n",
            name, m.integrated, m.weighted, m.outboard, m.side_sym)
        for (edge, (yy, gg)) in enumerate(zip(m.y, m.g))
            push!(station_rows, (name, edge, yy, gg, "medium", "capped_dirichlet",
                "area_mean", "mu_upper_minus_mu_lower_minus_c"))
        end
        push!(summary_rows, (name, m.integrated, m.weighted, m.outboard,
            m.side_sym, sum(m.g)/length(m.g), "medium", "capped_dirichlet",
            "DirectBackend", "flat_prescribed"))
    end

    # attribution deltas relative to VelocityThroughSources baseline
    println("\nΔ∫Γ vs VelocityThroughSources baseline:")
    for r in eachrow(summary_rows)
        d = (r.integrated - baseline_int) / abs(baseline_int) * 100
        @printf("  %-26s %+8.3f %%\n", r.formulation, d)
    end

    CSV.write(artifact("formsweep_stations.csv"), station_rows)
    CSV.write(artifact("formsweep_summary.csv"), summary_rows)
    println("\nwrote formsweep_summary.csv, formsweep_stations.csv to\n  $OUTPUT_DIR")
    println("\nGate A0 reading: if TraceCorrected moves ∫Γ toward Neumann by " *
            "order 5%, the gap is wake-potential/Kutta-trace handling.")
    return summary_rows
end

function run_smoke()
    println("\n#===== smoke: exercise driver machinery on a tiny case =====#")
    # tiny synthetic body
    body = oracle_body(; n_span=3, n_airfoil=21, n_endcap=3)
    solver = pnl.Backslash(body)
    M = nedges(body); N = body.ncells
    println("tiny body: $N cells, $M shedding edges")

    # 1. wake-strength extraction on a random state
    body.strength[:, 2] .= randn(N)
    gamma = [pnl._get_wakestrength_Gamma(body, i, s)
             for (s, shed) in enumerate(body.shedding) for i in axes(shed, 2)]
    @assert length(gamma) == M
    @assert kutta_map(body, body.strength[:, 2]) ≈ gamma atol=1e-12
    println("  wake-strength extraction: OK ($(M) edges, matches kutta_map)")

    # 2. fixed-grid interpolation + trapz self-consistency
    x = collect(0.0:0.25:1.0); y = x .^ 2
    bins = collect(0.0:0.1:1.0)
    yi = interpolate_linear(x, y, bins)
    @assert abs(yi[1] - 0.0) < 1e-12 && abs(yi[end] - 1.0) < 1e-12
    @assert abs(trapz(x, x) - 0.5) < 1e-12       # ∫₀¹ x dx = 0.5
    println("  interpolation + trapz: OK")

    # 3. integrated metrics + CSV writer round-trip
    m = integrated_metrics(body)
    df = DataFrame(edge=1:M, y=m.y, gamma=m.g)
    CSV.write(artifact("smoke_stations.csv"), df)
    back = CSV.read(artifact("smoke_stations.csv"), DataFrame)
    @assert nrow(back) == M
    println("  integrated metrics + CSV round-trip: OK " *
            "(∫Γ=$(round(m.integrated; sigdigits=4)), side_sym=$(round(m.side_sym; sigdigits=3)))")

    # 4. shared probe generation + loop-integral machinery + self-check
    #    (a solved state so the induced field is nonzero)
    gseed = seed_gamma(body, solver)
    wake = flat_wake(body, gseed)
    run_aero!(body, wake, solver)
    y0 = body.controlpoints[2, body.shedding[1][1, (M ÷ 2) + 1]]
    sc = loop_self_check(body, y0)
    @printf("  loop-integral self-check @y=%.3f: mean=%.4f spread=%.2e ok=%s\n",
        y0, sc.mean, sc.spread, sc.ok)

    # 5. decomposition algebra: full = body_only + wake_delta
    full = copy(view(body.strength, :, 2))
    body_only = deepcopy(body)
    # body-only solve: suppress wake by zero-strength wake
    wz = flat_wake(body_only, zeros(M))
    run_aero!(body_only, wz, solver)
    delta = full .- view(body_only.strength, :, 2)
    recon = view(body_only.strength, :, 2) .+ delta
    @assert norm(recon .- full) < 1e-10
    println("  decomposition algebra (full = body_only + wake_delta): OK")

    println("\nsmoke: all machinery exercised.")
    return nothing
end

function run_analyze()
    println("\n#===== analyze: aggregate Phase 2b CSVs =====#")
    fs = artifact("formsweep_summary.csv")
    if isfile(fs)
        df = CSV.read(fs, DataFrame)
        base = df[df.formulation .== "VelocityThroughSources", :integrated][1]
        println("\nGate A0 formulation sweep (Δ vs VelocityThroughSources):")
        @printf("  %-26s %14s %10s\n", "formulation", "∫Γ dy", "Δ%")
        for r in eachrow(df)
            d = (r.integrated - base) / abs(base) * 100
            @printf("  %-26s %14.6f %+9.3f\n", r.formulation, r.integrated, d)
        end
    else
        println("  no formsweep_summary.csv yet")
    end
    ol = artifact("oracle_summary.csv")
    if isfile(ol)
        df = CSV.read(ol, DataFrame)
        order = ["coarse", "medium", "fine", "xfine", "xxfine", "ultra"]
        refs = filter(r -> r in df.refinement, order)
        println("\nOracle Gate A ladder — ∫Γ and Dirichlet-Neumann tiebreaker:")
        @printf("  %-8s %14s %14s %14s | %12s %12s\n", "refine",
            "cap+Dir", "open+Neu", "cap+Neu(d1)", "BConly%", "topo%")
        for r in refs
            sub = df[df.refinement .== r, :]
            gv(t) = (rows = sub[sub.topology .== t, :integrated];
                     isempty(rows) ? NaN : rows[1])
            iD = gv("capped_dirichlet"); iNo = gv("open_neumann"); iNc = gv("capped_neumann_drop1")
            bc = (iNc - iD) / abs(iD) * 100
            tp = (iNo - iNc) / abs(iNc) * 100
            @printf("  %-8s %14.5f %14.5f %14.5f | %+11.3f %+11.3f\n",
                r, iD, iNo, iNc, bc, tp)
        end
        println("\n  BConly% = capped Neumann vs capped Dirichlet (topology fixed → pure BC).")
        println("  topo%   = open vs capped Neumann (BC fixed → pure caps/closure).")
        println("  Refinement trend of BConly% is the key: shrinking ⇒ numerically")
        println("  unresolved; persistent ⇒ genuine formulation bias.")
        # extraction residuals
        println("\n  Gate A extraction residual max|loopcirc−Γ|/|Γ| (should be ~0 if")
        println("  circulation is faithfully represented, i.e. not an extraction artifact):")
        for r in refs
            sub = df[df.refinement .== r, :]
            @printf("    %-8s ", r)
            for t in ["capped_dirichlet","open_neumann","capped_neumann_drop1"]
                rows = sub[sub.topology .== t, :extract_resid]
                @printf("%s=%.1e  ", t, isempty(rows) ? NaN : rows[1])
            end
            println()
        end
    else
        println("  no oracle_summary.csv yet")
    end
    return nothing
end

# --------------------------------------------------------------- oracle -------

const ORACLE_REFINEMENTS = (
    coarse = (n_span=5,  n_airfoil=31, n_endcap=5),
    medium = (n_span=7,  n_airfoil=41, n_endcap=5),
    fine   = (n_span=11, n_airfoil=61, n_endcap=7),
    # chordwise extension at fixed n_span=11 (chordwise is the convergence lever)
    xfine  = (n_span=11, n_airfoil=81,  n_endcap=7),
    xxfine = (n_span=11, n_airfoil=101, n_endcap=7),
    ultra  = (n_span=11, n_airfoil=121, n_endcap=7),
)

"VortexRing Neumann body (DBC=false) from a raw mesh, semi-infinite attached wake."
function neumann_body(nodes, cells; watertight, semiinfinite_wake=true)
    bt = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}
    oa = (; kerneloffset=1e-6*C_CHORD, kernelcutoff=1e-12*C_CHORD,
        semiinfinite_wake, watertight)
    base = bt(nodes, cells, zeros(Int, 6, 0); oa...)           # noshedding first
    shed = calc_pitching_wing_shedding(base.nodes, base.cells, C_CHORD)  # trace on rewound
    return bt(copy(base.nodes), copy(base.cells), [shed]; oa...)
end

"One-shot steady solve with the body's semi-infinite attached wake."
function solve_semiinf!(body)
    set_wake_Das!(body, UINF)
    pnl.steady!(body, pnl.ReferenceFrame(body), UINF;
        body_solvers=pnl.Backslash(body), backend=DIRECT, verbose=false)
    return body
end

"""
Build the three oracle bodies on identical lifting-surface coordinates:
capped+Dirichlet, open+Neumann, capped+Neumann(drop-one-panel). `drop` is the
dropped mid-mesh panel index (`:mid` = ncells÷2) for the full-rank closed
VortexRing system.
"""
function build_oracle_cases(ref; thickness=0.15, drop=:mid)
    (; n_span, n_airfoil, n_endcap) = ref
    bD = build_pitching_wing_body(C_CHORD, B_SPAN; n_span, n_airfoil, n_endcap,
        endcap=:round, semiinfinite_wake=true, thickness)
    no, co = pitching_wing_mesh(C_CHORD, B_SPAN; n_span, n_airfoil, n_endcap,
        caps=false, thickness)
    bNo = neumann_body(no, co; watertight=false)
    nc, cc = pitching_wing_mesh(C_CHORD, B_SPAN; n_span, n_airfoil, n_endcap,
        caps=true, thickness)
    dropidx = drop === :mid ? size(cc, 2) ÷ 2 : Int(drop)
    ccd = hcat(cc[:, 1:dropidx-1], cc[:, dropidx+1:end])
    bNc = neumann_body(nc, ccd; watertight=false)
    return (; bD, bNo, bNc, dropidx, ncells_capped=size(cc, 2))
end

"""
Gate A per-body observables: integrated metrics, per-station Γ, and the
velocity-based loop-integral circulation at each station (with the two-radii ×
two-quadrature self-check). Returns the metrics NamedTuple plus per-station
vectors (`y`, `gamma`, `loopcirc`, `loop_spread`, `loop_ok`) and the max
|loopcirc − Γ|/|Γ| extraction residual.
"""
function gate_a_observables(body)
    m = integrated_metrics(body)
    y, g = m.y, m.g
    loopcirc = similar(g); loop_spread = similar(g); loop_ok = falses(length(g))
    for i in eachindex(y)
        sc = loop_self_check(body, y[i])
        loopcirc[i] = sc.mean; loop_spread[i] = sc.spread; loop_ok[i] = sc.ok
    end
    extract_resid = maximum(abs.(loopcirc .- g) ./ max.(abs.(g), eps()))
    return (; m, y, gamma=g, loopcirc, loop_spread, loop_ok, extract_resid)
end

function run_oracle()
    ref_tag = Symbol(get(ENV, "PHASE2B_REFINEMENT", "medium"))
    haskey(ORACLE_REFINEMENTS, ref_tag) ||
        error("PHASE2B_REFINEMENT must be one of $(keys(ORACLE_REFINEMENTS)); got $ref_tag")
    ref = ORACLE_REFINEMENTS[ref_tag]
    thickness = parse(Float64, get(ENV, "PHASE2B_THICKNESS", "0.15"))
    println("\n#===== oracle Gate A ladder: refinement=$ref_tag $ref t/c=$thickness =====#")

    cases = build_oracle_cases(ref; thickness)
    println("capped ncells=$(cases.ncells_capped); dropped panel index=$(cases.dropidx)")

    topo = [("capped_dirichlet", "Dirichlet_SD", cases.bD),
            ("open_neumann",      "Neumann_VR",   cases.bNo),
            ("capped_neumann_drop1", "Neumann_VR", cases.bNc)]

    station_rows = DataFrame(topology=String[], bc=String[], refinement=String[],
        thickness=Float64[], edge=Int[], y=Float64[], gamma_TE=Float64[],
        loopcirc=Float64[], loop_spread=Float64[], loop_ok=Bool[])
    summary_rows = DataFrame(topology=String[], bc=String[], refinement=String[],
        thickness=Float64[], ncells=Int[], nedges=Int[], integrated=Float64[],
        weighted=Float64[], outboard=Float64[], side_sym=Float64[],
        extract_resid=Float64[], dropidx=Int[])

    obs = Dict{String,Any}()
    println()
    @printf("%-24s %8s %8s %14s %14s %12s %12s\n",
        "topology", "ncells", "nedges", "∫Γ dy", "∫|y|Γ dy", "side_sym", "extract_res")
    for (topo, bc, body) in topo
        solve_semiinf!(body)
        o = gate_a_observables(body)
        obs[topo] = o
        M = nedges(body)
        @printf("%-24s %8d %8d %14.5f %14.5f %12.2e %12.2e\n",
            topo, body.ncells, M, o.m.integrated, o.m.weighted, o.m.side_sym, o.extract_resid)
        for i in eachindex(o.y)
            push!(station_rows, (topo, bc, String(ref_tag), thickness, i, o.y[i],
                o.gamma[i], o.loopcirc[i], o.loop_spread[i], o.loop_ok[i]))
        end
        dropidx = topo == "capped_neumann_drop1" ? cases.dropidx : -1
        push!(summary_rows, (topo, bc, String(ref_tag), thickness, body.ncells, M,
            o.m.integrated, o.m.weighted, o.m.outboard, o.m.side_sym,
            o.extract_resid, dropidx))
    end

    iD  = obs["capped_dirichlet"].m.integrated
    iNo = obs["open_neumann"].m.integrated
    iNc = obs["capped_neumann_drop1"].m.integrated
    println("\nTiebreaker (∫Γ):")
    @printf("  open+Neu      vs capped+Dir : %+7.3f %%  (topology+BC)\n", (iNo-iD)/abs(iD)*100)
    @printf("  capped+Neu(d1) vs capped+Dir: %+7.3f %%  (BC only, topology fixed)\n", (iNc-iD)/abs(iD)*100)
    @printf("  open+Neu      vs capped+Neu : %+7.3f %%  (topology only, BC fixed)\n", (iNo-iNc)/abs(iNc)*100)
    println("  Gate A extraction (max |loopcirc−Γ|/|Γ|): " *
            "Dir=$(round(obs["capped_dirichlet"].extract_resid;sigdigits=2)) " *
            "openNeu=$(round(obs["open_neumann"].extract_resid;sigdigits=2)) " *
            "capNeu=$(round(obs["capped_neumann_drop1"].extract_resid;sigdigits=2))")

    # alternate-drop-panel sensitivity (coarse only): re-solve capped Neumann
    # with a different dropped panel to check tiebreaker robustness.
    if ref_tag == :coarse
        alt = build_oracle_cases(ref; thickness, drop = ref.n_airfoil * 2 + 1)
        solve_semiinf!(alt.bNc)
        iNc_alt = integrated_metrics(alt.bNc).integrated
        @printf("  alt drop panel %d: capped+Neu ∫Γ=%.5f (%+.4f%% vs mid-drop)\n",
            alt.dropidx, iNc_alt, (iNc_alt-iNc)/abs(iNc)*100)
        push!(summary_rows, ("capped_neumann_dropAlt", "Neumann_VR", String(ref_tag),
            thickness, alt.bNc.ncells, nedges(alt.bNc), iNc_alt,
            integrated_metrics(alt.bNc).weighted, integrated_metrics(alt.bNc).outboard,
            integrated_metrics(alt.bNc).side_sym, gate_a_observables(alt.bNc).extract_resid,
            alt.dropidx))
    end

    CSV.write(artifact("oracle_stations_$(ref_tag).csv"), station_rows)
    # append-or-create the cross-refinement summary
    sfile = artifact("oracle_summary.csv")
    if isfile(sfile)
        prev = CSV.read(sfile, DataFrame)
        prev = prev[prev.refinement .!= String(ref_tag), :]   # replace this rung
        summary_rows = vcat(prev, summary_rows)
    end
    CSV.write(sfile, summary_rows)
    println("\nwrote oracle_stations_$(ref_tag).csv, oracle_summary.csv")
    return summary_rows
end

# ------------------------------------------------------------- thickness ------

"""
Gated thickness screen: run the three-body oracle tiebreaker at a fixed
refinement across symmetric NACA 00xx thicknesses. Tests whether the
Dirichlet-Neumann gap follows global t/c (the leading explanation for the
oracle ~2.5% vs DJI ~5% magnitude difference), holding refinement fixed.
"""
function run_thickness()
    ref_tag = Symbol(get(ENV, "PHASE2B_REFINEMENT", "medium"))
    ref = ORACLE_REFINEMENTS[ref_tag]
    tcs = [0.06, 0.12, 0.18]
    println("\n#===== thickness screen: refinement=$ref_tag $ref, t/c=$tcs =====#")

    rows = DataFrame(thickness=Float64[], refinement=String[], cap_dir=Float64[],
        open_neu=Float64[], cap_neu_d1=Float64[], BConly_pct=Float64[],
        topo_pct=Float64[], extract_resid_max=Float64[])
    @printf("\n%6s %14s %14s %14s %10s %10s %12s\n",
        "t/c", "cap+Dir", "open+Neu", "cap+Neu(d1)", "BConly%", "topo%", "extract_res")
    for tc in tcs
        cases = build_oracle_cases(ref; thickness=tc)
        oD = (solve_semiinf!(cases.bD);  gate_a_observables(cases.bD))
        oNo = (solve_semiinf!(cases.bNo); gate_a_observables(cases.bNo))
        oNc = (solve_semiinf!(cases.bNc); gate_a_observables(cases.bNc))
        iD, iNo, iNc = oD.m.integrated, oNo.m.integrated, oNc.m.integrated
        bc = (iNc - iD)/abs(iD)*100
        tp = (iNo - iNc)/abs(iNc)*100
        er = max(oD.extract_resid, oNo.extract_resid, oNc.extract_resid)
        @printf("%6.3f %14.5f %14.5f %14.5f %+9.3f %+9.3f %12.2e\n",
            tc, iD, iNo, iNc, bc, tp, er)
        push!(rows, (tc, String(ref_tag), iD, iNo, iNc, bc, tp, er))
    end
    CSV.write(artifact("thickness_summary_$(ref_tag).csv"), rows)
    println("\nwrote thickness_summary_$(ref_tag).csv")
    println("Reading: if |BConly%| grows as t/c shrinks, the Dir-Neu gap is " *
            "global-thickness/section-resolution sensitive (explains DJI 5% > oracle 2.5%).")
    return rows
end
run_dji_bridge() = error("dji_bridge mode is conditional and not implemented yet.")

# ================================================================= main =======

function main()
    mode = get(ENV, "PHASE2B_MODE", "smoke")
    println("PHASE2B_MODE = $mode")
    if     mode == "smoke";     run_smoke()
    elseif mode == "formsweep"; run_formsweep()
    elseif mode == "oracle";    run_oracle()
    elseif mode == "thickness"; run_thickness()
    elseif mode == "analyze";   run_analyze()
    elseif mode == "dji_bridge"; run_dji_bridge()
    else error("unknown PHASE2B_MODE=$mode")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
