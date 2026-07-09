#=##############################################################################
# DESCRIPTION
#   Post-process the settled 2.0R rotor-hover cycle (revs 29-39) without
#   re-simulating: re-integrate the three saved gauge-pressure fields
#   (Laplace material-derivative, Laplace lamb-vector, Bernoulli) and evaluate
#   a Kutta-Joukowski cross-check, each with and without hub panels, over a
#   chain of warm-start run directories. Built on `pnl.replay`, which
#   reconstructs body/wake/frame state from the runs' VTK + metadata manifest.
#
#   Only the KJ monitors do N-body work (one induced-velocity evaluation per
#   step at panel-edge midpoints); the pressure-based thrust is a pure
#   re-integration of saved panel pressures.
#
# ENV knobs:
#   REPLAY_MODE          "full" (default: replay chain + CSV + plot),
#                        "validate" (5-step gate vs the original runs' CSVs),
#                        "plot" (re-plot from existing CSV), "none" (define only)
#   HUB_CUTOFF_R_OVER_R  cylindrical-radius cutoff for hub exclusion (default 0.1,
#                        the driver's blade-shedding root r/R)
#   REPLAY_OUT           output dir (default data/rotor_hover_replay2p0_forces)
#
# Usage: julia --project examples/rotor_hover_replay_forces.jl
=###############################################################################

import FLOWPanel as pnl
const ReadVTK = pnl.ReadVTK
using Printf

# PythonPlot is only needed by `plot_replay_ct`; import it at top level (a
# runtime `require` hits world-age errors) but skip it for compute-only modes.
if get(ENV, "REPLAY_MODE", "full") in ("full", "plot")
    import PythonPlot
end

# ---------------- constants matching the original runs -----------------------
# (dji9443 rotor, RHPC_MESH=40_40, NT=36, RPM=6000; see
#  examples/rotor_hover_pressure_comparison.jl and data/rotor_hover_sweeps.md)
const ROTOR_R = 0.119               # rotor radius [m]
const RHO = 1.179                   # air density [kg/m^3]
const NT_PER_REV = 36               # steps per revolution
const AXIAL_DIMENSION = 1           # thrust axis in frame 1 (dji9443 mesh)
const EXPERIMENT_CT = 0.072         # dji9443 hover reference

const DEFAULT_CHAIN = [
    (joinpath("data", "rotor_hover_relaxfilter2p0_ws"),
        "rotor_hover_relaxfilter2p0_ws", 1044:1223),
    (joinpath("data", "rotor_hover_relaxfilter2p0_ws_ext"),
        "rotor_hover_relaxfilter2p0_ws_ext", 1224:1403),
]

# VTU cell-data names of the saved gauge-pressure fields (fixed by the monitor
# tuple order of the original runs; verified on the saved files).
const PRESSURE_FIELDS = (
    laplace_matderiv = "gauge pressure",
    laplace_lamb     = "gauge pressure (PressureLaplace #3)",
    bernoulli        = "gauge pressure (PressureBernoulli #5)",
)

# ---------------- saved-pressure register monitor ----------------------------

"""
Replay-only monitor that reads a named gauge-pressure cell array from the saved
body VTU of the current step and registers it as the `:P` context field, so a
downstream `ForceMonitor` integrates that field instead of the default one.
`selected` maps the replay-local step counter to the on-disk step index.
"""
struct SavedPressureRegister
    path::String
    body_name::String
    field::String
    selected::Vector{Int}
    i_system::Int
end

pnl.monitor_provides(::SavedPressureRegister) = (:P,)

function pnl._run_monitor!(m::SavedPressureRegister, ctx::pnl.MonitorContext,
                           systems, wakes,
                           frames::AbstractVector{<:pnl.ReferenceFrame},
                           uinf, i_step::Int, dt::Real, t=nothing)
    idx = m.selected[i_step + 1]
    vtu = joinpath(m.path, m.body_name, "$(m.body_name).$(idx).vtu")
    cell_data = ReadVTK.get_cell_data(ReadVTK.VTKFile(vtu))
    m.field in keys(cell_data) || error(
        "Pressure field '$(m.field)' not found in $(vtu); " *
        "available: $(collect(keys(cell_data)))")
    P = Vector{Float64}(ReadVTK.get_data(cell_data[m.field]))
    pnl.monitor_register!(ctx, :P, m.i_system, P)
    return nothing
end

# ---------------- monitor tuple for one chain link ---------------------------

nohub_select(hub_cutoff_r_over_R) =
    cp -> sqrt(cp[2]^2 + cp[3]^2) >= hub_cutoff_r_over_R * ROTOR_R

"""
Build the replay monitor tuple for one chain link: three
(register, ForceMonitor-all, ForceMonitor-nohub) triplets — one per saved
pressure field — followed by two KuttaJoukowskiForce instances (all / no-hub).
"""
function build_replay_monitors(path, name, selected, systems;
                               hub_cutoff_r_over_R)
    nt = length(selected)
    body_name = name * "_body1"
    nohub = nohub_select(hub_cutoff_r_over_R)
    normalization() = pnl.RotorNormalization(RHO, 2 * ROTOR_R, 1)
    make_fm(select) = pnl.ForceMonitor(nt, 1;
        i_frame=1, normalization=normalization(),
        correct_kuttacondition=false,   # matches the original runs
        select, file=false)
    make_reg(field) = SavedPressureRegister(path, body_name, field, selected, 1)
    kj_backend = pnl.FastMultipoleBackend(;
        expansion_order=3, multipole_acceptance=0.4, leaf_size=1000)
    make_kj(select) = pnl.KuttaJoukowskiForce(systems[1], nt, 1;
        rho=RHO, backend=kj_backend, i_frame=1,
        normalization=normalization(), select, file=false)

    return (
        make_reg(PRESSURE_FIELDS.laplace_matderiv),
        make_fm(nothing), make_fm(nohub),
        make_reg(PRESSURE_FIELDS.laplace_lamb),
        make_fm(nothing), make_fm(nohub),
        make_reg(PRESSURE_FIELDS.bernoulli),
        make_fm(nothing), make_fm(nohub),
        make_kj(nothing), make_kj(nohub),
    )
end

# Column order of the consolidated CSV; maps CSV name -> monitor tuple index in
# `build_replay_monitors` output.
const CT_COLUMNS = (
    (:CT_laplace_matderiv, 2), (:CT_laplace_matderiv_nohub, 3),
    (:CT_laplace_lamb, 5), (:CT_laplace_lamb_nohub, 6),
    (:CT_bernoulli, 8), (:CT_bernoulli_nohub, 9),
    (:CT_kj, 10), (:CT_kj_nohub, 11),
)

# ---------------- chain replay ------------------------------------------------

"""
    replay_saved_forces(chain=DEFAULT_CHAIN; hub_cutoff_r_over_R=0.1, verbose=true)

Replay every step of `chain` (`[(path, name, steprange), ...]`) through the
force monitors and return `(steps, cts)` where `steps::Vector{Int}` are the
on-disk step indices concatenated across the chain and `cts` is a NamedTuple
of CT vectors (sign convention: thrust positive, i.e. `-force[axial]`).
"""
function replay_saved_forces(chain=DEFAULT_CHAIN;
        hub_cutoff_r_over_R=parse(Float64, get(ENV, "HUB_CUTOFF_R_OVER_R", "0.1")),
        verbose=true)
    steps = Int[]
    cts = NamedTuple{Tuple(first.(CT_COLUMNS))}(ntuple(_ -> Float64[], length(CT_COLUMNS)))
    hub_reported = false
    for (path, name, steprange) in chain
        selected = collect(steprange)
        verbose && println("replay_saved_forces: $(name), steps $(first(selected))-$(last(selected)) ($(length(selected)) steps)")
        factory = (systems, wakes, frames, t_range) ->
            build_replay_monitors(path, name, selected, systems;
                                  hub_cutoff_r_over_R)
        elapsed = @elapsed result = pnl.replay(path, name;
            monitor_factory=factory, steps=selected, verbose)
        verbose && @printf("replay_saved_forces: %s done in %.1f s\n", name, elapsed)
        if !hub_reported
            mask = result.monitors[3].select_mask
            println("replay_saved_forces: hub cutoff r/R=$(hub_cutoff_r_over_R) " *
                "excludes $(count(==(false), mask)) of $(length(mask)) panels")
            hub_reported = true
        end
        append!(steps, selected)
        for (colname, i_monitor) in CT_COLUMNS
            append!(getproperty(cts, colname),
                    -result.monitors[i_monitor].force[AXIAL_DIMENSION, :])
        end
    end
    return steps, cts
end

function write_replay_csv(steps, cts, out_dir)
    isdir(out_dir) || mkpath(out_dir)
    csv_path = joinpath(out_dir, "replay_CT_vs_rev.csv")
    colnames = Tuple(first.(CT_COLUMNS))
    open(csv_path, "w") do io
        println(io, "step,revolution," * join(String.(colnames), ","))
        for (k, step) in enumerate(steps)
            rev = step / NT_PER_REV
            vals = join((getproperty(cts, c)[k] for c in colnames), ",")
            println(io, "$step,$rev,$vals")
        end
    end
    println("Wrote consolidated CT history: $csv_path")
    return csv_path
end

# ---------------- validation gate ---------------------------------------------

"""
    validate_replay(; chain=DEFAULT_CHAIN, steprange=1050:1054, rtol=1e-6)

Cheap correctness gate: replayed full-panel CT from the saved pressures must
reproduce the original run's per-monitor force CSVs at the same steps (same P,
same geometry — validates loaders + integration before trusting new numbers).
Also checks `select=cp->false` yields zero force and no-hub != all-panels.
"""
function validate_replay(; chain=DEFAULT_CHAIN, steprange=1050:1054, rtol=1e-6)
    path, name, _ = chain[1]
    selected = collect(steprange)
    hub_cutoff = parse(Float64, get(ENV, "HUB_CUTOFF_R_OVER_R", "0.1"))
    monitors = nothing
    zero_fm = pnl.ForceMonitor(length(selected), 1;
        i_frame=1, normalization=pnl.RotorNormalization(RHO, 2 * ROTOR_R, 1),
        correct_kuttacondition=false, select=cp -> false, file=false)
    factory = (systems, wakes, frames, t_range) -> begin
        monitors = build_replay_monitors(path, name, selected, systems;
                                         hub_cutoff_r_over_R=hub_cutoff)
        return (monitors..., zero_fm)
    end
    pnl.replay(path, name; monitor_factory=factory, steps=selected, verbose=true)

    # Reference: the original run's per-monitor force CSVs, keyed by absolute
    # step index. monitor02 = Laplace matderiv, 04 = Laplace lamb, 06 = Bernoulli.
    read_force_csv(i_monitor) = begin
        f = joinpath(path, "monitors",
            @sprintf("%s_monitor%02d_force_system1.csv", name, i_monitor))
        d = Dict{Int, Float64}()
        for line in Iterators.drop(eachline(f), 1)
            parts = split(line, ",")
            d[parse(Int, parts[1])] = parse(Float64, parts[3])  # CFx
        end
        d
    end
    refs = ((2, "CT_laplace_matderiv"), (4, "CT_laplace_lamb"), (6, "CT_bernoulli"))
    ok = true
    for ((i_csv, label), i_mon) in zip(refs, (2, 5, 8))
        ref = read_force_csv(i_csv)
        for (k, step) in enumerate(selected)
            got = monitors[i_mon].force[AXIAL_DIMENSION, k]
            want = ref[step]
            pass = isapprox(got, want; rtol)
            pass || (ok = false)
            @printf("  %-22s step %d: replay=% .12f original=% .12f  %s\n",
                    label, step, got, want, pass ? "OK" : "MISMATCH")
        end
    end
    maxzero = maximum(abs, zero_fm.force)
    println("  select=cp->false max |CF| = $(maxzero) $(maxzero == 0 ? "OK" : "FAIL")")
    maxzero == 0 || (ok = false)
    for (label, i_all, i_nohub) in (("laplace_matderiv", 2, 3), ("kj", 10, 11))
        dct = maximum(abs, monitors[i_all].force .- monitors[i_nohub].force)
        println("  $(label): max |CT_all - CT_nohub| = $(dct) " *
                (dct > 0 ? "OK (mask active)" : "FAIL (mask inert)"))
        dct > 0 || (ok = false)
    end
    println(ok ? "validate_replay: PASSED" : "validate_replay: FAILED")
    return ok
end

# ---------------- statistics + plot -------------------------------------------

function read_replay_csv(csv_path)
    lines = readlines(csv_path)
    header = Symbol.(split(lines[1], ","))
    cols = Dict(h => Float64[] for h in header)
    for line in lines[2:end]
        for (h, v) in zip(header, split(line, ","))
            push!(cols[h], parse(Float64, v))
        end
    end
    return header, cols
end

cycle_stats(rev, ct, rev_lo, rev_hi) = begin
    sel = [v for (r, v) in zip(rev, ct) if rev_lo <= r <= rev_hi && isfinite(v)]
    isempty(sel) && return (NaN, NaN)
    m = sum(sel) / length(sel)
    s = sqrt(sum((v - m)^2 for v in sel) / max(length(sel) - 1, 1))
    return (m, s)
end

"""
    plot_replay_ct(csv_path; experiment=0.072, warmup_revs=1)

Plot CT vs revolution for all monitors (solid = all panels, dashed = no-hub),
with the experimental reference line. Laplace curves are masked for the first
`warmup_revs` after the chain start (PressureLaplace du/dt warm-up after a warm
start). Legend reports cycle means over revs 30-39. Saves a PNG next to the CSV.
"""
function plot_replay_ct(csv_path; experiment=EXPERIMENT_CT, warmup_revs=1.0,
                        stats_rev_lo=30.0, stats_rev_hi=Inf)
    plt = PythonPlot
    _, cols = read_replay_csv(csv_path)
    rev = cols[:revolution]
    rev0 = minimum(rev)
    monitors = (
        (:CT_laplace_matderiv, "Laplace (Du/Dt)", "tab:blue", true),
        (:CT_laplace_lamb, "Laplace (lamb)", "tab:green", true),
        (:CT_bernoulli, "Bernoulli", "tab:orange", false),
        (:CT_kj, "Kutta-Joukowski", "tab:red", false),
    )
    fig, ax = plt.subplots(figsize=(9, 5.5))
    stats = Dict{Symbol, NTuple{2, Float64}}()
    for (col, label, color, is_laplace) in monitors
        for (suffix, style) in (("", "-"), ("_nohub", "--"))
            name = Symbol(string(col) * suffix)
            ct = copy(cols[name])
            if is_laplace   # mask warm-up
                for (k, r) in enumerate(rev)
                    r < rev0 + warmup_revs && (ct[k] = NaN)
                end
            end
            m, s = cycle_stats(rev, ct, stats_rev_lo, stats_rev_hi)
            stats[name] = (m, s)
            ax.plot(rev, ct; color, linestyle=style, linewidth=1.2,
                label=@sprintf("%s%s: %.4f +/- %.4f", label,
                               suffix == "" ? "" : " (no hub)", m, s))
        end
    end
    # KJ carries a sharp fixed-azimuth 1/rev spike (probes passing near the
    # wake); overlay per-rev means so its cycle behavior stays readable.
    for (name, color) in ((:CT_kj, "tab:red"), (:CT_kj_nohub, "darkred"))
        ct = cols[name]
        revbins = sort(unique(floor.(Int, rev)))
        means = [let sel = [v for (r, v) in zip(rev, ct) if floor(Int, r) == b]
                     sum(sel) / length(sel)
                 end for b in revbins]
        ax.plot(revbins .+ 0.5, means; color, linestyle="none", marker="o",
                markersize=4, label="$(name) per-rev mean")
    end
    ax.axhline(experiment; color="k", linestyle=":", linewidth=1.5,
               label="experiment $(experiment)")
    ax.set_ylim(0.05, 0.08)   # KJ spikes go off-scale; std in legend reflects them
    ax.set_xlabel("revolution")
    ax.set_ylabel("CT")
    ax.set_title("Rotor hover CT, replayed 2.0R-filter cycle (means over revs $(stats_rev_lo)+)")
    ax.grid(true, alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    png_path = splitext(csv_path)[1] * ".png"
    fig.savefig(png_path, dpi=150)
    println("Wrote plot: $png_path")

    println("\nCycle statistics (revs $(stats_rev_lo)-$(stats_rev_hi == Inf ? "end" : stats_rev_hi)):")
    for (col, label, _, _) in monitors
        for suffix in ("", "_nohub")
            name = Symbol(string(col) * suffix)
            m, s = stats[name]
            @printf("  %-32s CT = %.4f +/- %.4f  (gap to %.3f: %+.4f)\n",
                    string(name), m, s, experiment, m - experiment)
        end
    end
    return stats
end

# ---------------- main ---------------------------------------------------------

function main()
    mode = get(ENV, "REPLAY_MODE", "full")
    out_dir = get(ENV, "REPLAY_OUT", joinpath("data", "rotor_hover_replay2p0_forces"))
    if mode == "validate"
        validate_replay()
    elseif mode == "full"
        steps, cts = replay_saved_forces()
        csv_path = write_replay_csv(steps, cts, out_dir)
        plot_replay_ct(csv_path)
    elseif mode == "plot"
        plot_replay_ct(joinpath(out_dir, "replay_CT_vs_rev.csv"))
    elseif mode == "none"
        # define-only (for interactive use)
    else
        error("Unknown REPLAY_MODE=$(mode); use full, validate, plot, or none")
    end
end

if get(ENV, "REPLAY_MODE", "full") != "none"
    main()
end
