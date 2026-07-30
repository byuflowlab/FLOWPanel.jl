#!/usr/bin/env julia
# Aggregate the rotor-hover wake-decay audit (BRAINSTORM/004).
#
# Scans every run directory data/<AUDIT_PREFIX>_<variant>/ produced by
# examples/run_rotor_hover_wake_decay_audit.sh, then:
#   * reads each run's <run_name>_CT_vs_rev.csv  (CT already stored as positive thrust)
#   * extracts a per-step particle-count history from the .vtp wake files
#     (NumberOfPoints header), writing <run_dir>/<run_name>_particle_count.csv
#   * writes data/<AUDIT_PREFIX>_aggregate/<AUDIT_PREFIX>_summary.csv with, per variant:
#         name, peak_CT, rev_at_peak, final_CT, drift_peak_to_final, crossed_0072,
#         final_particle_count
#   * plots a CT-vs-rev overlay (with the 0.072 reference line) and a
#     particle-count-vs-rev overlay across all variants.
#
# ENV overrides:
#   AUDIT_PREFIX  run-name prefix              (default wake_decay_audit)
#   RPM, NT       timing, to tag .vtp steps    (defaults 6000, 36; match the harness)
#   CT_SERIES     which CT column to overlay   (default CT_bernoulli)
#   CT_TARGET     reference CT line            (default 0.072)

using DelimitedFiles
using PythonPlot

prefix    = get(ENV, "AUDIT_PREFIX", "wake_decay_audit")
RPM       = parse(Float64, get(ENV, "RPM", "6000"))
nt        = parse(Int, get(ENV, "NT", "36"))
dt        = 60 / RPM / nt
ct_series = get(ENV, "CT_SERIES", "CT_bernoulli")
ct_target = parse(Float64, get(ENV, "CT_TARGET", "0.072"))

agg_dir = joinpath("data", "$(prefix)_aggregate")
isdir(agg_dir) || mkpath(agg_dir)

# ---- discover variant run directories ---------------------------------------
# data/<prefix>_<variant>/, excluding the aggregate dir itself.
all_dirs = filter(isdir, readdir("data"; join=true))
run_dirs = filter(all_dirs) do d
    b = basename(d)
    startswith(b, "$(prefix)_") && b != "$(prefix)_aggregate"
end
sort!(run_dirs)
isempty(run_dirs) && error("No run directories matching data/$(prefix)_* found.")

variant_name(run_dir) = replace(basename(run_dir), "$(prefix)_" => "")

# ---- CT CSV reader ----------------------------------------------------------
function read_ct(run_dir)
    run_name = basename(run_dir)
    csv = joinpath(run_dir, "$(run_name)_CT_vs_rev.csv")
    isfile(csv) || return nothing
    raw, header = readdlm(csv, ',', header=true)
    header = vec(string.(header))
    idx = findfirst(==(ct_series), header)
    idx === nothing && error("Column $(ct_series) not in $(csv)")
    rev = Float64.(raw[:, findfirst(==("revolution"), header)])
    ct  = Float64.(raw[:, idx])
    return rev, ct
end

# ---- particle-count history from .vtp headers (reused pattern) --------------
step_index(fname) = parse(Int, split(basename(fname), '.')[end-1])
function read_num_points(path)
    io = open(path, "r")
    try
        m = match(r"NumberOfPoints=\"(\d+)\"", String(read(io, 2048)))
        m === nothing && error("Could not find NumberOfPoints in $(path)")
        return parse(Int, m.captures[1])
    finally
        close(io)
    end
end

function read_particle_counts(run_dir)
    run_name = basename(run_dir)
    pdir = joinpath(run_dir, "$(run_name)_wake1_particles")
    isdir(pdir) || return nothing
    vtps = filter(f -> endswith(f, ".vtp"), readdir(pdir; join=true))
    isempty(vtps) && return nothing
    sort!(vtps; by=step_index)
    steps  = step_index.(vtps)
    counts = read_num_points.(vtps)
    revs   = steps .* dt .* RPM ./ 60
    # persist per-run count history
    open(joinpath(run_dir, "$(run_name)_particle_count.csv"), "w") do io
        println(io, "step,revolution,particle_count")
        for i in eachindex(steps)
            println(io, "$(steps[i]),$(revs[i]),$(counts[i])")
        end
    end
    return revs, counts
end

# ---- gather, summarize, plot ------------------------------------------------
ct_fig, ct_ax = PythonPlot.subplots(figsize=(9, 5.5))
pc_fig, pc_ax = PythonPlot.subplots(figsize=(9, 5.5))

summary_rows = Vector{String}()
push!(summary_rows,
    "name,peak_CT,rev_at_peak,final_CT,drift_peak_to_final,crossed_0072,final_particle_count")

for run_dir in run_dirs
    name = variant_name(run_dir)
    ct_data = read_ct(run_dir)
    if ct_data === nothing
        @warn "skipping $(name): no CT CSV"
        continue
    end
    rev, ct = ct_data
    if all(isnan, ct)
        @warn "skipping $(name): CT series $(ct_series) is all NaN"
        continue
    end

    peak_ct   = maximum(skipmissing(filter(!isnan, ct)))
    peak_k    = findfirst(==(peak_ct), ct)
    rev_peak  = rev[peak_k]
    final_ct  = ct[findlast(!isnan, ct)]
    drift     = peak_ct - final_ct
    crossed   = any(c -> !isnan(c) && c >= ct_target, ct)

    pc = read_particle_counts(run_dir)
    final_pc = pc === nothing ? -1 : pc[2][end]

    push!(summary_rows,
        "$(name),$(peak_ct),$(rev_peak),$(final_ct),$(drift),$(crossed),$(final_pc)")

    ct_ax.plot(rev, ct, label=name, linewidth=1.3)
    if pc !== nothing
        pc_ax.plot(pc[1], pc[2], label=name, linewidth=1.3)
    end
end

# CT overlay
ct_ax.axhline(ct_target, color="k", linestyle="--", linewidth=1.0,
    label="target $(ct_target)")
ct_ax.set_xlabel("Revolution")
ct_ax.set_ylabel("Thrust coefficient \$C_T\$")
ct_ax.set_title("Wake-decay audit: \$C_T\$ vs revolution ($(ct_series))")
ct_ax.grid(true, alpha=0.3)
ct_ax.legend(fontsize=8, ncol=2)
ct_fig.tight_layout()
ct_png = joinpath(agg_dir, "$(prefix)_CT_vs_rev.png")
ct_fig.savefig(ct_png, dpi=150)

# particle-count overlay
pc_ax.set_xlabel("Revolution")
pc_ax.set_ylabel("Particle count")
pc_ax.set_title("Wake-decay audit: particle count vs revolution")
pc_ax.grid(true, alpha=0.3)
pc_ax.legend(fontsize=8, ncol=2)
pc_fig.tight_layout()
pc_png = joinpath(agg_dir, "$(prefix)_particles_vs_rev.png")
pc_fig.savefig(pc_png, dpi=150)

# summary table
summary_csv = joinpath(agg_dir, "$(prefix)_summary.csv")
open(summary_csv, "w") do io
    for row in summary_rows
        println(io, row)
    end
end

println("Wrote:")
println("  $(ct_png)")
println("  $(pc_png)")
println("  $(summary_csv)")
println()
foreach(println, summary_rows)
