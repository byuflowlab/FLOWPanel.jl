#=##############################################################################
BRAINSTORM 021 Phase 3 — warmstart analysis (exit criteria 1 and 2).

Reads the per-step rows written by benchmark/rotor_hover_solver_unsteady.jl
(results/phase3/<mode>/<rung>/unsteady.csv) and emits:

  1. iteration-count and wall-time reduction, warmstart type x solver x rung,
     as wake-developed mean +/- spread against that config's cold baseline.
     The iteration metric is `niter_first` -- the FIRST inner solve of each
     step, i.e. the one that consumes the step-to-step initial guess. Legacy
     `niter` (the LAST inner solve of the step) is reported as a diagnostic
     column only: sampling it after the step is what made the pre-fix smoke
     look like a warmstart reduction (021 Phase 3, 2026-08-23);
  2. break-even step count -- steps needed for cumulative warm savings to
     amortize the warmstart machinery's setup/bookkeeping overhead;
  3. the order sweep for the extrapolated arms;
  4. a standalone TikZ figure + same-named data dir, per the figure policy.

The wake-developed window drops each row set's leading `skip_steps` steps (the
post-restart history-fill transient); the count is carried in the CSV rather
than assumed here.

Validity guards folded into the `guard` column: the cold-vs-warm trajectory
identity check, the count of steps whose inner solves did not converge, and --
for these one-body rotor runs -- the requirement that `nsolves == 1` on every
scored row. `nsolves == -1` is "unavailable" (direct/generic solver), reported
as such rather than as a pass. Inputs predating the Phase 3 schema (no
`niter_first` / `nsolves`) are refused.

The validity check was REBUILT 2026-08-25 (Ryan). It no longer compares the
arms' SOLUTIONS to each other. Instead it asks the direct question -- does each
arm still satisfy the BOUNDARY CONDITIONS it promised to? -- from the driver's
per-step `bc_error!` pass: for this Dirichlet body the control-point potential
IS the BC residual, measured under a requested FMM error 10x below the arm's own
stopping tolerance, with certification reported. PASS = every measured step has
|phi|_inf < that tolerance and certified.

Two superseded designs, both recorded so neither is re-attempted: the original
`sum(abs, strength)` vs a hard 1e-8 (read a residual tolerance as a solution
tolerance, and the L1 aggregate cancels), and a relative-L2-between-arms bar
calibrated against a tolerance-null arm (correct but needed a conditioning
argument and an extra arm to say what a residual check says directly).

Identity is REPORTED, not asserted: CT per step is the physical identity signal
and particle counts are shown for gross divergence. The wake is chaotic, so tiny
differences compound over many steps and exact trajectory identity is not
expected -- no threshold is attached to either.

Run:
  julia --project benchmark/phase3_analysis.jl
Env:
  PHASE_DIR  results subdir to read (default benchmark/results/phase3)
  MODE       threading mode subdir  (default multi)
  FIG_DIR    figure output dir (default BRAINSTORM/021_.../figures)
=###############################################################################

using LinearAlgebra: norm

const RESULTS = get(ENV, "PHASE_DIR",
                    joinpath(@__DIR__, "results", get(ENV, "PHASE", "phase3")))
const MODE = get(ENV, "MODE", "multi")
const FIG_DIR = get(ENV, "FIG_DIR", normpath(joinpath(@__DIR__, "..", "BRAINSTORM",
    "021_rotor_hover_solver_benchmarks", "figures")))


# --- minimal CSV reader (the unsteady schema quotes only trailing `notes`) ---
function read_csv(path)
    isfile(path) || return nothing
    lines = filter(!isempty, readlines(path))
    length(lines) >= 2 || return nothing
    cols = Dict(c => i for (i, c) in enumerate(split(lines[1], ",")))
    rows = [split(l, ",") for l in lines[2:end]]
    return cols, rows
end

cell(cols, r, name) = haskey(cols, name) && length(r) >= cols[name] ? r[cols[name]] : ""
fnum(cols, r, name) = (v = cell(cols, r, name); isempty(v) ? NaN : something(tryparse(Float64, v), NaN))

mean(v) = isempty(v) ? NaN : sum(v) / length(v)
function stdev(v)
    length(v) < 2 && return 0.0
    m = mean(v)
    return sqrt(sum((x - m)^2 for x in v) / (length(v) - 1))
end

# --- gather every rung's rows ------------------------------------------------
struct Arm
    rung::String
    config::String
    warmstart::String
    order::Int
    niter_first::Vector{Float64}  # FIRST inner solve of each step (the metric)
    niter::Vector{Float64}        # LAST inner solve (Phase 1/2 legacy diagnostic)
    nsolves::Vector{Float64}      # per-body solves per step (1 for a one-body run)
    t_solve::Vector{Float64}
    t_setup::Float64
    n_particles::Float64
    checksum::Float64        # sum(abs, strength); legacy, pre-rebuild rows only
    unconverged::Int
    nsteps::Int
    # Per-step BC check (2026-08-25). `ratio` is |phi|_inf / this arm's own
    # absolute acceptance level, so <1 means the arm kept the promise it made.
    bc_ratio::Vector{Float64}
    bc_spread::Vector{Float64}   # |phi|_inf / median|phi|, the distribution's tail
    bc_cert::Vector{Bool}
    bc_eps_frac::Vector{Float64} # requested FMM error / tol: the measurement's
                                 # own certified uncertainty, as a fraction
    ct::Vector{Float64}          # per-step CT (identity metric)
    npart::Vector{Float64}       # per-step particle count (REPORTED, not asserted)
end

function collect_arms()
    arms = Arm[]
    root = joinpath(RESULTS, MODE)
    isdir(root) || (@warn "no results under $root"; return arms)
    for rung in sort(readdir(root))
        parsed = read_csv(joinpath(root, rung, "unsteady.csv"))
        parsed === nothing && continue
        cols, rows = parsed
        groups = Dict{Tuple{String,String,String,Int},Vector{Any}}()
        for r in rows
            key = (cell(cols, r, "rung"), cell(cols, r, "config"),
                   cell(cols, r, "warmstart"),
                   round(Int, something(tryparse(Float64, cell(cols, r, "warmstart_order")), 0.0)))
            push!(get!(groups, key, Any[]), r)
        end

        # Phase 3 requires the widened schema. Without these columns the only
        # available iteration count is `niter` (the LAST inner solve), which is
        # what made the voided smoke table look like a warmstart reduction.
        for c in ("niter_first", "nsolves")
            haskey(cols, c) || error(
                "$(joinpath(root, rung, "unsteady.csv")) predates the Phase 3 " *
                "schema (missing column `$c`). Re-run the arm with the current " *
                "benchmark/rotor_hover_solver_unsteady.jl; pre-Phase-3 rows " *
                "cannot be scored as warmstart comparisons.")
        end
        for (key, grp) in groups
            skip = round(Int, maximum(fnum(cols, r, "skip_steps") for r in grp))
            steps = [round(Int, fnum(cols, r, "step")) for r in grp]
            keep = [i for (i, st) in enumerate(steps) if st > skip]
            isempty(keep) && continue
            niter = [fnum(cols, grp[i], "niter") for i in keep]
            nfirst = [fnum(cols, grp[i], "niter_first") for i in keep]
            nsolv = [fnum(cols, grp[i], "nsolves") for i in keep]
            tsol = [fnum(cols, grp[i], "t_solve") for i in keep]
            unconv = count(i -> lowercase(cell(cols, grp[i], "solved")) == "false", keep)
            setup = maximum(x -> isnan(x) ? 0.0 : x, (fnum(cols, r, "t_setup") for r in grp))
            npart_end = maximum(x -> isnan(x) ? 0.0 : x, (fnum(cols, r, "n_particles_end") for r in grp))
            csum = maximum(x -> isnan(x) ? 0.0 : x, (fnum(cols, r, "strength_checksum") for r in grp))
            # BC rows only where the step was actually measured (BCERR_EVERY
            # can subsample); a blank cell is "not measured", never a pass.
            meas = [i for i in keep if !isnan(fnum(cols, grp[i], "bcerr_max"))]
            ratio = [fnum(cols, grp[i], "bcerr_max") / fnum(cols, grp[i], "bcerr_tol")
                     for i in meas]
            spread = [(m = fnum(cols, grp[i], "bcerr_med");
                       m == 0 ? NaN : fnum(cols, grp[i], "bcerr_max") / m) for i in meas]
            cert = [lowercase(cell(cols, grp[i], "bcerr_certified")) == "true" for i in meas]
            epsf = [fnum(cols, grp[i], "bcerr_eps") / fnum(cols, grp[i], "bcerr_tol")
                    for i in meas]
            ct = filter(!isnan, [fnum(cols, grp[i], "CT") for i in keep])
            npart = filter(!isnan, [fnum(cols, grp[i], "n_particles") for i in keep])
            push!(arms, Arm(key[1], key[2], key[3], key[4], nfirst, niter,
                            nsolv, tsol, setup, npart_end, csum, unconv, length(grp),
                            ratio, spread, cert, epsf, ct, npart))        end
    end
    return arms
end

arms = collect_arms()
isempty(arms) && (println("No Phase 3 rows found under $(joinpath(RESULTS, MODE)) — nothing to analyze."); exit(0))

baseline(rung, config) = findfirst(a -> a.rung == rung && a.config == config &&
                                        a.warmstart == "cold", arms)


# ---------------------------------------------------------------------------
# Validity guard: does each arm still satisfy the BOUNDARY CONDITIONS?
# (Ryan 2026-08-25 — supersedes both the strength-checksum guard and the
# relative-L2-between-arms design that briefly replaced it.)
# ---------------------------------------------------------------------------
# WHAT A WARM START CAN AND CANNOT DO. It changes only WHERE the iteration
# begins. The solver's acceptance test is a RESIDUAL bound, absolute and
# initial-guess-independent (FastMultipole/src/solve.jl:1222), so a converged
# warm arm is in exactly the same acceptance set as the cold arm.
#
# The old guard compared the ARMS' SOLUTIONS and asked whether they agreed to
# 1e-8. That is the wrong question twice over: it read a residual tolerance as a
# solution tolerance (implicitly ||A^-1|| = 1, where Phase 1 measures ~12x at
# R1), and `sum(abs, strength)` cancels equal-and-opposite changes. It reported
# the campaign's own extrapolated arm as DIVERGED at 2.4e-7 — roughly three
# decades INSIDE the solver's own solution ambiguity.
#
# The direct question is whether each arm satisfies the BCs it promised to. The
# driver answers it per step with one `bc_error!` pass: for this Dirichlet body
# the control-point potential IS the BC residual phi_sigma + G_mu x, and the
# pass requests an absolute FMM error 10x below the arm's own stopping tolerance
# and reports whether it certified that bound. No reference solution, no dense
# G, no conditioning argument.
#
#   bc_ratio  = |phi|_inf / (that arm's own absolute acceptance level)
#   PASS      = every measured step has bc_ratio < 1 AND certified
#
# IDENTITY IS REPORTED, NOT ASSERTED. The wake is chaotic: tiny differences
# compound over many steps, so demanding identical particle counts or
# trajectories would flag correct runs. CT is the physical quantity the campaign
# cares about and is the identity signal; `n_particles` is shown per step to
# make gross divergence visible. Neither carries a threshold.

"Order statistics of a vector, NaN-safe, for compact reporting."
_mx(v) = isempty(v) ? NaN : maximum(v)
_md(v) = isempty(v) ? NaN : (w = sort(v); w[clamp(cld(length(w), 2), 1, length(w))])

guard_of = Dict{Tuple{String,String,String,Int},String}()
for a in arms
    key = (a.rung, a.config, a.warmstart, a.order)
    g = if isempty(a.bc_ratio)
        # blank bcerr_* means the step was never measured. Report that, never
        # a pass — the same convention as `nsolves == -1` below.
        "no BC rows (pre-2026-08-25 driver, or BCERR_EVERY skipped every step)"
    else
        worst = _mx(a.bc_ratio)
        nunc = count(!, a.bc_cert)
        # The arm stopped on ITS OWN residual estimate (leaf-local, at its tuned
        # expansion order); `bc_error!` re-evaluates the same quantity
        # independently and more sharply. The two legitimately disagree by about
        # the error we REQUESTED of the measurement, so the pass band is
        # 1 + eps/tol -- the measurement's own certified uncertainty, not a
        # fudge factor. Ratios inside that band are reported as marginal, so a
        # drift toward the bar stays visible instead of being rounded to "OK".
        band = 1 + (isempty(a.bc_eps_frac) ? 0.0 : _mx(a.bc_eps_frac))
        r3 = round(worst; sigdigits=3)
        base = worst <= 1 ?
            "BC OK (max |φ|/tol = $r3 over $(length(a.bc_ratio)) steps)" :
            worst <= band ?
            "BC OK, marginal (max |φ|/tol = $r3, within the ±$(round(100*(band-1); sigdigits=2))% measurement band)" :
            "**BC VIOLATED** (max |φ|/tol = $r3 > band $(round(band; sigdigits=3)))"
        nunc == 0 ? base : base * " / **$nunc UNCERTIFIED**"
    end
    a.unconverged > 0 && (g *= " / **$(a.unconverged) UNCONVERGED STEPS**")
    guard_of[key] = g
end

# --- exit criterion 1: reduction table --------------------------------------
println("\n## Warmstart reduction, wake-developed window (mode=$MODE)\n")
println("`niter_first` (FIRST inner solve of each step) is the warmstart metric: " *
        "it is the solve that consumes the step-to-step guess. Legacy `niter` " *
        "(LAST inner solve) is shown only as a diagnostic.\n")
println("| rung | config | warmstart | order | steps | niter_first mean ± sd | Δiter % | t_solve mean ± sd [s] | Δt % | break-even [steps] | legacy niter | guard |")
println("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")

fig_rows = Tuple{String,String,String,Int,Float64,Float64}[]

for a in sort(arms, by = x -> (x.rung, x.config, x.warmstart, x.order))
    ib = baseline(a.rung, a.config)
    bi = ib === nothing ? nothing : arms[ib]
    mi, si = mean(a.niter_first), stdev(a.niter_first)
    mt, st = mean(a.t_solve), stdev(a.t_solve)
    d_iter = bi === nothing || mean(bi.niter_first) == 0 ? NaN :
        100 * (mi - mean(bi.niter_first)) / mean(bi.niter_first)
    d_t = bi === nothing || mean(bi.t_solve) == 0 ? NaN : 100 * (mt - mean(bi.t_solve)) / mean(bi.t_solve)

    # Break-even: the warmstart machinery's overhead is the extra setup this arm
    # paid over the cold arm (history allocation + projection bookkeeping);
    # per-step saving is the cold minus warm solve time. Report the step count
    # at which cumulative saving overtakes that overhead.
    be = if bi === nothing || a.warmstart == "cold"
        "—"
    else
        saving = mean(bi.t_solve) - mt
        overhead = a.t_setup - bi.t_setup
        saving <= 0 ? "never (no saving)" :
            overhead <= 0 ? "0 (no overhead)" : string(ceil(Int, overhead / saving))
    end

    # Trajectory-identity verdict, computed above against a MEASURED bar.
    guard = guard_of[(a.rung, a.config, a.warmstart, a.order)]

    # One-body solve-count check (021 Phase 3). These are 1-tuple rotor runs, so
    # a valid warmstart comparison needs exactly one per-body solve per step;
    # more means the step re-solved and niter_first is not a step-to-step
    # measurement. nsolves == -1 is "unavailable" (direct/generic solver), which
    # is explicitly NOT a pass. Not generalized to multibody: real block-GS outer
    # counts legitimately vary per timestep there.
    if isempty(a.nsolves) || any(isnan, a.nsolves)
        guard *= " / nsolves missing"
    elseif all(==(-1.0), a.nsolves)
        guard *= " / nsolves n/a"
    else
        nbad = count(!=(1.0), a.nsolves)
        nbad == 0 || (guard *= " / **$nbad ROWS nsolves != 1 (INVALID)**")
    end

    println("| $(a.rung) | $(a.config) | $(a.warmstart) | $(a.order) | $(length(a.niter)) | " *
            "$(round(mi; digits=2)) ± $(round(si; digits=2)) | $(round(d_iter; digits=1)) | " *
            "$(round(mt; digits=4)) ± $(round(st; digits=4)) | $(round(d_t; digits=1)) | $be | " *
            "$(round(mean(a.niter); digits=2)) | $guard |")
    push!(fig_rows, (a.rung, a.config, a.warmstart, a.order, mi, mt))
end

# --- BC satisfaction and identity report -------------------------------------
println("\n## BC satisfaction and identity (mode=$MODE)\n")
println("Validity: each arm's own per-step BC residual |φ| = |φ_σ + G_μ x| " *
        "against the absolute tolerance THAT ARM promised, measured by one " *
        "`bc_error!` pass per step at a requested FMM error 10× below that " *
        "tolerance. `certified` means every M2L pair proved its own error " *
        "bound. The arm stopped on its own coarser residual estimate, so the " *
        "pass band is 1 + (requested error)/tol — the measurement's own " *
        "certified uncertainty; the `band` column states it per arm.\n")
println("Identity: CT and particle count are REPORTED, not asserted — the wake " *
        "is chaotic, so small differences compound and exact trajectory " *
        "identity is not expected. Δ columns are against the cold arm.\n")
println("| rung | config | arm | BC steps | max \\|φ\\|/tol | band | med \\|φ\\|/tol | tail max/med | certified | CT final | ΔCT % | n_part final | Δn_part |")
println("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
for a in sort(arms, by = x -> (x.rung, x.config, x.warmstart, x.order))
    ib = baseline(a.rung, a.config)
    bi = ib === nothing ? nothing : arms[ib]
    arm = a.warmstart * (a.order > 0 ? string(a.order) : "")
    ct_f = isempty(a.ct) ? NaN : a.ct[end]
    ct_b = bi === nothing || isempty(bi.ct) ? NaN : bi.ct[end]
    dct = isnan(ct_f) || isnan(ct_b) || ct_b == 0 ? NaN : 100 * (ct_f - ct_b) / ct_b
    np_f = isempty(a.npart) ? NaN : a.npart[end]
    np_b = bi === nothing || isempty(bi.npart) ? NaN : bi.npart[end]
    fmt(v) = isnan(v) ? "—" : string(round(v; sigdigits=4))
    println("| $(a.rung) | $(a.config) | $arm | $(length(a.bc_ratio)) | " *
            "$(fmt(_mx(a.bc_ratio))) | " *
            "$(isempty(a.bc_eps_frac) ? "—" : fmt(1 + _mx(a.bc_eps_frac))) | " *
            "$(fmt(_md(a.bc_ratio))) | " *
            "$(fmt(_md(a.bc_spread))) | " *
            "$(isempty(a.bc_cert) ? "—" : "$(count(a.bc_cert))/$(length(a.bc_cert))") | " *
            "$(fmt(ct_f)) | $(fmt(dct)) | $(fmt(np_f)) | " *
            "$(isnan(np_f) || isnan(np_b) ? "—" : string(round(Int, np_f - np_b))) |")
end
println("\n`tail max/med` is |φ|_∞ / median|φ| — how far the worst panel sits " *
        "above the typical one. A stable ratio across arms says the warm start " *
        "did not concentrate error anywhere new.\n")

# --- exit criterion 2: order sweep ------------------------------------------
orders = sort(unique(a.order for a in arms if a.warmstart == "extrap"))
if length(orders) > 1
    println("\n## Extrapolation order sweep\n")
    println("| rung | config | order | niter_first mean | t_solve mean [s] |")
    println("| --- | --- | --- | --- | --- |")
    for a in sort(filter(x -> x.warmstart == "extrap", arms), by = x -> (x.rung, x.config, x.order))
        println("| $(a.rung) | $(a.config) | $(a.order) | " *
                "$(round(mean(a.niter_first); digits=2)) | " *
                "$(round(mean(a.t_solve); digits=4)) |")
    end
else
    println("\n(Only one extrapolation order present — order sweep needs runs at " *
            "additional WARMSTART_ORDER values before exit criterion 2 is met.)")
end

# --- figure: standalone TikZ + same-named data dir ---------------------------
mkpath(FIG_DIR)
figname = "warmstart_iterations"
datadir = joinpath(FIG_DIR, figname)
mkpath(datadir)

configs = sort(unique(r[2] for r in fig_rows))
for ws in sort(unique(r[3] for r in fig_rows))
    open(joinpath(datadir, "$(ws).csv"), "w") do io
        println(io, "x,config,rung,order,niter_first_mean,t_solve_mean")
        for (i, cfg) in enumerate(configs)
            for r in fig_rows
                r[3] == ws && r[2] == cfg || continue
                println(io, "$i,$(r[2]),$(r[1]),$(r[4]),$(r[5]),$(r[6])")
            end
        end
    end
end

open(joinpath(FIG_DIR, "$(figname).tex"), "w") do io
    println(io, raw"""
\documentclass[tikz,border=2mm]{standalone}
\usepackage{pgfplots}
\pgfplotsset{compat=1.18}
\begin{document}
\begin{tikzpicture}
\begin{axis}[
    width=12cm, height=7cm,
    ylabel={first-solve iterations to target (wake-developed mean)},
    xtick=data, xticklabel style={rotate=30, anchor=east},
    xticklabels={""" * join(configs, ",") * raw"""},
    legend pos=north east, ymin=0, grid=major,
]""")
    for ws in sort(unique(r[3] for r in fig_rows))
        println(io, "\\addplot+[only marks] table[x=x, y=niter_first_mean, col sep=comma] " *
                    "{$(figname)/$(ws).csv};")
        println(io, "\\addlegendentry{$ws}")
    end
    println(io, raw"""
\end{axis}
\end{tikzpicture}
\end{document}""")
end

println("\nFigure written: $(joinpath(FIG_DIR, figname * ".tex")) (+ $(datadir)/)")
