#=##############################################################################
BRAINSTORM 021 — shared knob-selection + timing helpers, extracted verbatim
from rotor_hover_solver_phase1_table.jl (2026-08-18, pure code motion) so the
Phase 2 scripts reuse the frozen Phase-1 knob logic instead of duplicating it.

Include AFTER common.jl and phase1_case.jl. Call-time global requirements:
  rung, knobs_dir   (from phase1_case.jl)
  target_rel        (defined by the including script before calling
                     stage3_winner; 1e-6 for the campaign)
=###############################################################################

# ---- CSV helpers (read latest matching rows from the Stage 2/3 outputs) ----
function read_rows(path)
    isfile(path) || return nothing
    lines = readlines(path)
    cols = Dict(String(c) => i for (i, c) in enumerate(split(lines[1], ",")))
    rows = [split(l, ",") for l in lines[2:end] if !isempty(strip(l))]
    return (cols, rows)
end

"""
Shared FGS set (Ryan ruling 2026-08-17): the knob set is PINNED to the
per-rung τ=1e-6-tuned config; only the preconditioner sweep count is
selected end-to-end — fastest meets_target fgsprecond row whose knobs match
that config. The staircase-crossing check stays as a tripwire (coarse-tuned
configs have apply-accuracy BC floors that can sit above 1e-6; measured:
R2 τ=1e-4 config floors at BC ≈ 3.6e-6 — impossible for the 1e-6-tuned
config by construction).
"""
function stage3_winner()
    cfg6 = stage2_selected(1e-6)
    cfg6 === nothing && error("no τ=1e-6 row in fgstune_selected.csv for $rung")
    parsed = read_rows(joinpath(knobs_dir, "fgsprecond.csv"))
    parsed === nothing && return nothing
    cands = NamedTuple[]
    cols, rows = parsed
    for c in rows
        c[cols["rung"]] == rung || continue
        lowercase(c[cols["meets_target"]]) == "true" || continue
        knobs = (parse(Int, c[cols["p"]]), parse(Float64, c[cols["mac"]]),
                 parse(Int, c[cols["leaf"]]), parse(Int, c[cols["inner"]]))
        knobs == (cfg6.p, cfg6.mac, cfg6.leaf, cfg6.inner) || continue
        push!(cands, (; tau=parse(Float64, c[cols["tau"]]),
                 p=knobs[1], mac=knobs[2], leaf=knobs[3], inner=knobs[4],
                 sweeps=parse(Int, c[cols["sweeps"]]),
                 t=parse(Float64, c[cols["t_solve_min"]])))
    end
    isempty(cands) && error("no fgsprecond rows match the 1e-6-tuned config " *
        "($(cfg6.p)/$(cfg6.mac)/$(cfg6.leaf)/$(cfg6.inner)) — run fgsprecond " *
        "with SWEEP_LADDER_1E6=1")
    sort!(cands; by=c -> c.t)
    for c in cands
        sc = staircase_for(c.p, c.mac, c.leaf, c.inner)
        if sc !== nothing && any(t -> t[4] <= target_rel, sc)
            return c
        end
        println("stage3 row τ=$(c.tau) SKIPPED: staircase never crosses " *
                "BC $(target_rel) (tripwire — should be impossible for the " *
                "1e-6-tuned config)")
    end
    return nothing
end

"Stage 2 selected row for a given τ (latest row wins)."
function stage2_selected(tau)
    parsed = read_rows(joinpath(knobs_dir, "fgstune_selected.csv"))
    parsed === nothing && return nothing
    cols, rows = parsed
    out = nothing
    for c in rows
        c[cols["rung"]] == rung || continue
        parse(Float64, c[cols["tau"]]) == tau || continue
        out = (; p=parse(Int, c[cols["p"]]), mac=parse(Float64, c[cols["mac"]]),
                leaf=parse(Int, c[cols["leaf"]]), inner=parse(Int, c[cols["inner"]]))
    end
    return out
end

"Staircase for given knobs: sorted (iter, t_wall, mse, bc) tuples (latest run)."
function staircase_for(p, mac, leaf, inner)
    parsed = read_rows(joinpath(knobs_dir, "fgstune_staircase.csv"))
    parsed === nothing && return nothing
    cols, rows = parsed
    pts = Tuple{Int,Float64,Float64,Float64}[]
    for c in rows
        c[cols["rung"]] == rung || continue
        parse(Int, c[cols["p"]]) == p || continue
        parse(Float64, c[cols["mac"]]) == mac || continue
        parse(Int, c[cols["leaf"]]) == leaf || continue
        parse(Int, c[cols["inner"]]) == inner || continue
        it = parse(Int, c[cols["iter"]])
        # duplicate iters (rerun) — keep latest by overwriting below
        filter!(t -> t[1] != it, pts)
        push!(pts, (it, parse(Float64, c[cols["t_wall"]]),
                    something(tryparse(Float64, c[cols["mse_internal"]]), NaN),
                    parse(Float64, c[cols["bc_rel_l2"]])))
    end
    sort!(pts; by=first)
    return isempty(pts) ? nothing : pts
end

snap_half_decade(x) = 10.0^(round(2 * log10(x)) / 2)

"Post-crossing margin stopping tolerance (2026-08-15 fgstune edit): geometric
mean of the crossing residual and its first LOWER successor — threshold
stopping can neither stop early above target nor sit on the crossing value.
`sc` = staircase tuples (iter, t_wall, mse, bc); `i_cross` = BC crossing."
function margin_tol(sc, i_cross)
    r_cross = sc[i_cross][3]
    (isfinite(r_cross) && r_cross > 0) ||
        error("no positive finite residual at the BC crossing — extend fgstune MAXIT")
    for j in i_cross+1:length(sc)
        r = sc[j][3]
        isfinite(r) && 0 < r < r_cross && return sqrt(r_cross * r)
    end
    error("no finite residual below the crossing residual — extend fgstune MAXIT")
end

"""
Adaptive min-of-k (Ryan ruling 2026-08-18, amends ruling 5's fixed k>=5):
run one untimed-choice warmup solve (timed to CHOOSE k, excluded from the
min), then k timed reps: k=5 below 60 s, k=3 below 10 min, k=2 above.
Returns (t_min, k_used).
"""
function adaptive_min_of_k(f; setup! = nothing)
    setup! === nothing || setup!()
    t0 = time_ns(); f(); t1 = (time_ns() - t0) / 1e9
    k = t1 < 60 ? 5 : (t1 < 600 ? 3 : 2)
    t_min = Inf
    for _ in 1:k
        setup! === nothing || setup!()
        t0 = time_ns(); f()
        t_min = min(t_min, (time_ns() - t0) / 1e9)
    end
    return t_min, k
end
