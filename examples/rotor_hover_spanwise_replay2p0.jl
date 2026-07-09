#=##############################################################################
# DESCRIPTION
#   Spanwise loading percentile statistics for the settled 2.0R-filter
#   rotor-hover cycle (revs 30-39), replayed from saved VTK output without
#   re-simulating and without any N-body work: for each of the three saved
#   gauge-pressure fields (Laplace Du/Dt, Laplace lamb, Bernoulli), re-register
#   the saved pressures, re-integrate per-panel forces, and bin blade thrust
#   loading dT/dr per blade per step. Writes per-bin quantile stats (min/q25/
#   median/q75/max across the 324-step window) to CSV and one percentile plot
#   per pressure source.
#
#   Reuses SavedPressureRegister/PRESSURE_FIELDS/DEFAULT_CHAIN from
#   rotor_hover_replay_forces.jl (REPLAY_MODE=none) and the snapshot/quantile/
#   plotting machinery from rotor_hover_spanwise_loading_replay.jl
#   (SPANWISE_REPLAY_RUN=false).
#
# ENV knobs:
#   SPAN_OUT          output dir (default data/rotor_hover_replay2p0_forces/spanwise)
#   NBINS             spanwise bins per blade (default 35 = shedding stations)
#   SPANWISE_BINNING  control_point | span_overlap (default)
#   STEPS1 / STEPS2   step ranges for the two chain links, "a:b" format
#                     (defaults 1080:1223 and 1224:1403 = revs 30-39);
#                     STEPS2="" skips link 2 (for cheap gating)
#   CCBLADE_CSV       optional BEM overlay CSV (skipped if absent)
#
# Usage: julia -t auto --project examples/rotor_hover_spanwise_replay2p0.jl
=###############################################################################

ENV["REPLAY_MODE"] = "none"
include("rotor_hover_replay_forces.jl")          # SavedPressureRegister, PRESSURE_FIELDS, constants
ENV["SPANWISE_REPLAY_RUN"] = "false"
include("rotor_hover_spanwise_loading_replay.jl") # snapshot/quantiles/CSV/plotting
import PythonPlot as plt

# ---------------- configuration ----------------------------------------------

span_out = get(ENV, "SPAN_OUT",
    joinpath("data", "rotor_hover_replay2p0_forces", "spanwise"))
nbins = env_int("NBINS", 35)
binning_sym = parse_spanwise_binning(get(ENV, "SPANWISE_BINNING", "span_overlap"))
ccblade_csv = get(ENV, "CCBLADE_CSV",
    joinpath("data", "rotor_hover_ccblade", "rotor_hover_ccblade_sectional_ncrit4.csv"))
replay_ct_csv = joinpath("data", "rotor_hover_replay2p0_forces", "replay_CT_vs_rev.csv")

parse_link_steps(value) = begin
    v = strip(value)
    isempty(v) && return Int[]
    parse_steps(v)
end
link_steps = [
    parse_link_steps(get(ENV, "STEPS1", "1080:1223")),
    parse_link_steps(get(ENV, "STEPS2", "1224:1403")),
]

sources = keys(PRESSURE_FIELDS)   # (:laplace_matderiv, :laplace_lamb, :bernoulli)

# Shared binning: blades along +/-y in frame 1 (dji9443 2-blade rotor), thrust
# along -x; r_min/r_max auto-computed on the first snapshot call and shared
# across links so bin edges are identical everywhere.
e_axial = unit_axis(AXIAL_DIMENSION)
thrust_dir = -e_axial
blade_dirs = fallback_blade_dirs(2, AXIAL_DIMENSION, 2)
binning = SpanwiseReplayBinning(blade_dirs, e_axial, thrust_dir, ROTOR_R,
    fill(NaN, 2), fill(NaN, 2), nbins, 1, binning_sym)

println("Rotor hover spanwise loading replay (2.0R-filter chain)")
println("  chain:   ", join(("$(name) $(step_description(sel))"
    for ((_, name, _), sel) in zip(DEFAULT_CHAIN, link_steps) if !isempty(sel)), " + "))
println("  nbins:   $(nbins) ($(binning_sym))")
println("  sources: ", join(String.(sources), ", "))
println("  output:  $(span_out)")

# ---------------- replay chain ------------------------------------------------

link_snapshots = Vector{Dict{Symbol, SpanwiseReplaySnapshot}}()

for ((path, name, _), selected) in zip(DEFAULT_CHAIN, link_steps)
    isempty(selected) && continue
    snaps = Dict{Symbol, SpanwiseReplaySnapshot}()
    factory = (systems, wakes, frames, t_range) -> begin
        nt = length(t_range)
        normalization = pnl.RotorNormalization(RHO, 2 * ROTOR_R, 1)
        monitors = Any[]
        for source in sources
            reg = SavedPressureRegister(path, name * "_body1",
                getproperty(PRESSURE_FIELDS, source), selected, 1)
            fm = pnl.ForceMonitor(nt, 1; i_frame=1, normalization,
                correct_kuttacondition=false, verbose=false,
                file=false, vtk_fields=())
            snap = SpanwiseReplaySnapshot(source, binning, nt)
            push!(monitors, reg, fm, snap)
            snaps[source] = snap
        end
        return Tuple(monitors)
    end
    elapsed = @elapsed pnl.replay(path, name;
        monitor_factory=factory, steps=selected, verbose=true)
    @printf("link %s done in %.1f s (%.2f s/step)\n", name, elapsed,
            elapsed / length(selected))
    push!(link_snapshots, snaps)
end

isempty(link_snapshots) && error("No chain links selected (STEPS1/STEPS2 both empty)")
all_steps = vcat((sel for sel in link_steps if !isempty(sel))...)

# ---------------- concatenate links across the window --------------------------

combined = Dict{Symbol, SpanwiseReplaySnapshot}()
for source in sources
    vals = vcat((snaps[source].values for snaps in link_snapshots)...)
    cnts = vcat((snaps[source].counts for snaps in link_snapshots)...)
    combined[source] = SpanwiseReplaySnapshot(source, binning, vals, cnts)
end

nsteps_total = size(first(values(combined)).values, 1)
window_idxs = 1:nsteps_total
println("\nCombined window: $(nsteps_total) steps (disk steps $(first(all_steps)):$(last(all_steps)))")

# ---------------- stats CSV + plots --------------------------------------------

mkpath(span_out)
stats_path = joinpath(span_out, "spanwise_loading_stats.csv")
stats = write_stats_csv(stats_path, combined, collect(sources), window_idxs,
                        ROTOR_R, length(blade_dirs))
println("Wrote $(stats_path)")

ccblade = load_ccblade(ccblade_csv)
isnothing(ccblade) &&
    println("(no CCBlade CSV at $(ccblade_csv); plots omit BEM overlay)")

for source in sources
    png = joinpath(span_out, "spanwise_loading_$(source).png")
    plot_source(png, stats, source, ccblade)
    println("Wrote $(png)")
end

# ---------------- integral consistency check -----------------------------------
# Integrate the binned per-step loading back to total thrust and compare its
# window mean against the CT cycle-mean from the forces-replay CSV at the same
# steps. Expect ratio ~ 1 (small deficit: panels with degenerate radial
# projection or outside [r_min, r_max] are skipped by the binner).

n_rev_per_s = 6000 / 60                                    # RPM=6000
thrust_norm = RHO * n_rev_per_s^2 * (2 * ROTOR_R)^4        # rho n^2 D^4 [N]

ct_ref = if isfile(replay_ct_csv)
    _, cols = read_replay_csv(replay_ct_csv)
    stepset = Dict(Int(s) => k for (k, s) in enumerate(cols[:step]))
    idxs = [stepset[s] for s in all_steps if haskey(stepset, s)]
    length(idxs) == length(all_steps) ||
        @warn "Only $(length(idxs)) of $(length(all_steps)) window steps found in $(replay_ct_csv)"
    Dict(source => mean(cols[Symbol("CT_", source)][idxs]) for source in sources)
else
    println("(no reference CT CSV at $(replay_ct_csv); skipping consistency check)")
    nothing
end

println("\nIntegral consistency (binned dT/dr vs forces-replay CT, same steps):")
for source in sources
    snap = combined[source]
    T_steps = map(1:nsteps_total) do k
        T = 0.0
        for blade in axes(snap.values, 2)
            width = (binning.r_max[blade] - binning.r_min[blade]) / nbins
            for bin in 1:nbins
                v = snap.values[k, blade, bin]
                isfinite(v) && (T += v * width)
            end
        end
        T
    end
    T_mean = mean(T_steps)
    if isnothing(ct_ref)
        @printf("  %-17s integrated T = %.4f N (no reference)\n", String(source), T_mean)
    else
        T_ref = ct_ref[source] * thrust_norm
        @printf("  %-17s integrated T = %.4f N, reference T = %.4f N (CT %.5f), ratio = %.4f\n",
                String(source), T_mean, T_ref, ct_ref[source], T_mean / T_ref)
    end
end
