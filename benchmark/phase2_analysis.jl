#=##############################################################################
BRAINSTORM 021 Phase 2 — analysis + figure pipeline (rulings 4, 8, 9).

Reads solvetable.csv / phase2.csv / history sidecars across rungs (flat local
layout AND per-rung HPC dirs merged; latest row per (rung, config, row_kind)
wins by date) and emits, per the figure policy (standalone TikZ .tex + a
same-named CSV data dir, pdflatex-compilable):

  figures/solver_scaling.tex + solver_scaling/<config>.csv + fits.csv
      t_solve and t_setup vs N per config (log-log) + fitted exponents
      t ∝ N^p (ruling 9)
  figures/convergence_<rung>.tex + convergence_<rung>/<config>.csv
      residual vs wall-clock per iterative config on shared axes (ruling 4),
      from the phase2 history sidecars
  results/phase2_analysis/itertime_validation.csv
      equal-time-per-iteration check (ruling 4): mean per-iteration wall
      time + max drift %, flagged when drift > 10%
  results/phase2_analysis/memory.csv
      solver_state_bytes vs N per config where phase2.csv rows exist
      (ruling 8; the memory figure is added once HPC rows span the ladder)

The pipeline is the deliverable — rerun after the HPC sweeps land to
regenerate everything from the CSVs of record. Local rows are dev-only.

Run:  ANALYSIS_MODE=multi julia --project benchmark/phase2_analysis.jl
=###############################################################################

mode = get(ENV, "ANALYSIS_MODE", "multi")
repo = normpath(joinpath(@__DIR__, ".."))
item_dir = joinpath(repo, "BRAINSTORM", "021_rotor_hover_solver_benchmarks")
fig_dir = joinpath(item_dir, "figures")
ana_dir = joinpath(item_dir, "results", "phase2_analysis")
mkpath(fig_dir); mkpath(ana_dir)

RUNG_ORDER = ["R1", "R2", "R3", "R4", "R5", "R6", "R7"]
CONFIGS = ["backslash_ldiv", "krylov_gmres", "krylov_jacobi", "krylov_ilu",
           "fgs", "fgmres_fgs"]

# ---- CSV plumbing (schema-tolerant: rows keyed by header names) ------------
function read_dict_rows(path)
    isfile(path) || return NamedTuple[]
    lines = readlines(path)
    # header = first non-comment line (history sidecars carry a "# metric =" line)
    i_header = findfirst(l -> !isempty(strip(l)) && !startswith(l, "#"), lines)
    i_header === nothing && return NamedTuple[]
    cols = [String(strip(c)) for c in split(lines[i_header], ",")]
    out = []
    for l in lines[i_header+1:end]
        (isempty(strip(l)) || startswith(l, "#")) && continue
        vals = split(l, ",")
        length(vals) == length(cols) || continue
        push!(out, Dict(zip(cols, String.(vals))))
    end
    return out
end

"All rows of `csvname` for this mode: flat dir + per-rung subdirs merged."
function gather_rows(base, csvname)
    root = joinpath(@__DIR__, "results", base, mode)
    isdir(root) || return []
    paths = [joinpath(root, csvname);
             [joinpath(root, r, csvname) for r in RUNG_ORDER]]
    return reduce(vcat, (read_dict_rows(p) for p in paths); init=[])
end

"Latest row per (rung, config, row_kind) by the date column."
function latest_rows(rows)
    best = Dict{Tuple{String,String,String},Any}()
    for r in rows
        key = (r["rung"], r["config"], get(r, "row_kind", "standard"))
        if !haskey(best, key) || get(r, "date", "") >= get(best[key], "date", "")
            best[key] = r
        end
    end
    return best
end

fnum(r, k) = something(tryparse(Float64, get(r, k, "")), NaN)

solvetable = latest_rows(gather_rows("phase1", "solvetable.csv"))
phase2 = latest_rows(gather_rows("phase2", "phase2.csv"))

# prefer phase2 rows (componentized), fall back to solvetable
function cost_row(rung, config)
    for (tbl, kinds) in ((phase2, ("standard",)),
                         (solvetable, ("standard", "target_1e-6")))
        for kind in kinds
            haskey(tbl, (rung, config, kind)) && return tbl[(rung, config, kind)]
        end
    end
    return nothing
end

# ---- 1. scaling fits + figure data (ruling 9) ------------------------------
"Least-squares fit log t = log a + p log N; returns (p, a)."
function fit_powerlaw(Ns, ts)
    n = length(Ns)
    x = log10.(Float64.(Ns)); y = log10.(ts)
    mx = sum(x) / n; my = sum(y) / n
    p = sum((x .- mx) .* (y .- my)) / sum((x .- mx) .^ 2)
    return p, 10.0^(my - p * mx)
end

scaling_dir = joinpath(fig_dir, "solver_scaling")
mkpath(scaling_dir)
fits = String[]
present_configs = String[]
for config in CONFIGS
    rows = [(r, cost_row(r, config)) for r in RUNG_ORDER]
    filter!(x -> x[2] !== nothing, rows)
    length(rows) < 2 && continue
    push!(present_configs, config)
    Ns = [parse(Int, x[2]["n_panels"]) for x in rows]
    tsolve = [fnum(x[2], "t_solve_min") for x in rows]
    tsetup = [max(fnum(x[2], haskey(x[2], "t_setup_total") ? "t_setup_total" :
                             "t_setup"), 1e-6) for x in rows]
    open(joinpath(scaling_dir, "$(config).csv"), "w") do io
        println(io, "N,t_solve,t_setup")
        for (N, ts, tu) in zip(Ns, tsolve, tsetup)
            println(io, "$N,$ts,$tu")
        end
    end
    p_solve, _ = fit_powerlaw(Ns, tsolve)
    p_setup, _ = fit_powerlaw(Ns, tsetup)
    push!(fits, "$config,$(length(Ns)),$(round(p_solve; digits=3))," *
                "$(round(p_setup; digits=3))")
end
open(joinpath(scaling_dir, "fits.csv"), "w") do io
    println(io, "config,n_rungs,p_solve,p_setup")
    foreach(l -> println(io, l), fits)
end

plotstyle = Dict(
    "backslash_ldiv" => "color=black, mark=square*",
    "krylov_gmres"   => "color=red!70!black, mark=o",
    "krylov_jacobi"  => "color=orange!80!black, mark=triangle*",
    "krylov_ilu"     => "color=blue!70!black, mark=diamond*",
    "fgs"            => "color=green!50!black, mark=*",
    "fgmres_fgs"     => "color=violet, mark=pentagon*",
)
legname = Dict("backslash_ldiv" => "Backslash (ldiv)",
    "krylov_gmres" => "GMRES", "krylov_jacobi" => "GMRES+Jacobi",
    "krylov_ilu" => "GMRES+ILU", "fgs" => "FGS",
    "fgmres_fgs" => "FGMRES+FGS")

open(joinpath(fig_dir, "solver_scaling.tex"), "w") do io
    print(io, """
% 021 Phase 2 — per-solve cost vs panel count (ruling 9). Data:
% solver_scaling/<config>.csv (regenerated by benchmark/phase2_analysis.jl);
% fitted exponents in solver_scaling/fits.csv.
\\documentclass[tikz]{standalone}
\\usepackage{pgfplots}
\\pgfplotsset{compat=1.17}
\\begin{document}
\\begin{tikzpicture}
\\begin{loglogaxis}[width=11cm, height=8cm,
    xlabel={panel count \$N\$}, ylabel={per-solve time [s]},
    legend pos=north west, legend style={font=\\footnotesize},
    grid=both, minor grid style={gray!15}, major grid style={gray!35}]
""")
    for config in present_configs
        println(io, "\\addplot[$(plotstyle[config]), thick] table " *
            "[col sep=comma, x=N, y=t_solve] {solver_scaling/$(config).csv};")
        println(io, "\\addlegendentry{$(legname[config])}")
    end
    print(io, """
\\end{loglogaxis}
\\end{tikzpicture}
\\end{document}
""")
end

# ---- 2. convergence vs wall-clock per rung (ruling 4) ----------------------
# history sidecars: iter,t_wall,residual_internal,residual_true (phase2 dirs)
function history_paths(rung, config)
    root = joinpath(@__DIR__, "results", "phase2", mode)
    cands = [joinpath(root, "history_$(config)_$(rung).csv"),
             joinpath(root, rung, "history_$(config)_$(rung).csv")]
    return filter(isfile, cands)
end

conv_rungs = String[]
for rung in RUNG_ORDER
    dir = joinpath(fig_dir, "convergence_$(rung)")
    written = String[]
    for config in CONFIGS
        paths = history_paths(rung, config)
        isempty(paths) && continue
        rows = read_dict_rows(paths[end])
        isempty(rows) && continue
        mkpath(dir)
        open(joinpath(dir, "$(config).csv"), "w") do io
            println(io, "iter,t_wall,residual")
            for r in rows
                res = fnum(r, "residual_internal")
                isnan(res) && (res = fnum(r, "residual_true"))
                (isnan(res) || res <= 0) && continue
                println(io, "$(r["iter"]),$(fnum(r, "t_wall")),$res")
            end
        end
        push!(written, config)
    end
    isempty(written) && continue
    push!(conv_rungs, rung)
    open(joinpath(fig_dir, "convergence_$(rung).tex"), "w") do io
        print(io, """
% 021 Phase 2 — solver-internal residual vs wall-clock, $(rung) (ruling 4).
% Metrics differ per solver (Krylov preconditioned norm; FGS max-abs) — the
% shared axis is TIME; matched stopping is certified separately via the BC
% metric. Data: convergence_$(rung)/<config>.csv (benchmark/phase2_analysis.jl).
\\documentclass[tikz]{standalone}
\\usepackage{pgfplots}
\\pgfplotsset{compat=1.17}
\\begin{document}
\\begin{tikzpicture}
\\begin{semilogyaxis}[width=11cm, height=8cm,
    xlabel={wall-clock time [s]}, ylabel={solver-internal residual},
    legend pos=north east, legend style={font=\\footnotesize},
    grid=both, minor grid style={gray!15}, major grid style={gray!35}]
""")
        for config in written
            println(io, "\\addplot[$(plotstyle[config]), thick] table " *
                "[col sep=comma, x=t_wall, y=residual] " *
                "{convergence_$(rung)/$(config).csv};")
            println(io, "\\addlegendentry{$(legname[config])}")
        end
        print(io, """
\\end{semilogyaxis}
\\end{tikzpicture}
\\end{document}
""")
    end
end

# ---- 3. equal-time-per-iteration validation (ruling 4) ---------------------
open(joinpath(ana_dir, "itertime_validation.csv"), "w") do io
    println(io, "rung,config,n_iters,mean_iter_s,max_drift_pct,flat_within_10pct")
    for rung in RUNG_ORDER, config in CONFIGS
        paths = history_paths(rung, config)
        isempty(paths) && continue
        rows = read_dict_rows(paths[end])
        ts = [fnum(r, "t_wall") for r in rows]
        length(ts) < 3 && continue
        dts = diff(ts)
        m = sum(dts) / length(dts)
        drift = maximum(abs.(dts .- m)) / m * 100
        println(io, "$rung,$config,$(length(ts)),$(round(m; sigdigits=4))," *
            "$(round(drift; digits=1)),$(drift <= 10)")
    end
end

# ---- 4. memory vs N (ruling 8) ---------------------------------------------
open(joinpath(ana_dir, "memory.csv"), "w") do io
    println(io, "rung,config,n_panels,solver_state_bytes,solve_alloc_bytes")
    for rung in RUNG_ORDER, config in CONFIGS
        haskey(phase2, (rung, config, "standard")) || continue
        r = phase2[(rung, config, "standard")]
        println(io, "$rung,$config,$(r["n_panels"])," *
            "$(get(r, "solver_state_bytes", ""))," *
            "$(get(r, "solve_alloc_bytes", ""))")
    end
end

println("phase2 analysis written:")
println("  $(joinpath(scaling_dir, "fits.csv")) ($(length(fits)) configs)")
println("  figures: solver_scaling.tex" *
        (isempty(conv_rungs) ? "" : ", convergence_{$(join(conv_rungs, ","))}.tex"))
println("  $(joinpath(ana_dir, "itertime_validation.csv"))")
println("  $(joinpath(ana_dir, "memory.csv"))")
