## Predict the computational savings from particle merging in the rotor-hover run.
##
## Reads the per-step particle VTK (.vtp) files written by
## examples/rotor_hover_pressure_comparison.jl and extracts the actual particle
## count at each timestep (merging ON) straight from the PolyData header
## (NumberOfPoints="..."). It then computes the no-merge baseline as pure linear
## accumulation of shed particles:
##
##     particles_merge_off[k] = n_shed_locations * p_per_step * (shed events so far)
##
## where n_shed_locations counts shedding edges across BOTH blades and is derived
## empirically from the first nonzero merged count (= n_shed_locations * pps).
## No culling is applied: in this run the count plateau seen with merging on is
## driven by wake roll-up making merging more effective, not by the GlobalCylinder
## trim policy, so the merge-off field simply grows linearly.
##
## Output: a CSV with one row per timestep reporting merge-on and merge-off counts.
##
## ENV overrides:
##   RUN_DIR   -- directory of the completed run (default data/rotor_hover_pressure_comparison)
##   RUN_NAME  -- run name prefix          (default rotor_hover_pressure_comparison)
##   PPS       -- particles per step        (default 2)
##   OUT_CSV   -- output CSV path            (default <RUN_DIR>/particles_merging_savings.csv)

run_name = get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison")
run_dir  = get(ENV, "RUN_DIR", joinpath("data", run_name))
pps      = parse(Int, get(ENV, "PPS", "2"))

# Simulation timing (mirrored from rotor_hover_pressure_comparison.jl) so we can
# tag each step with its rotor revolution without loading the solver.
RPM = 6000
nt  = 36
dt  = 60 / RPM / nt

particles_dir = joinpath(run_dir, "$(run_name)_wake1_particles")
isdir(particles_dir) || error("Particle directory not found: $(particles_dir)")

# Collect the per-step .vtp files and sort numerically by their step index, which
# appears as the second-to-last dot-delimited token, e.g. "...particles.137.vtp".
step_index(fname) = parse(Int, split(basename(fname), '.')[end-1])
vtp_files = filter(f -> endswith(f, ".vtp"), readdir(particles_dir; join=true))
isempty(vtp_files) && error("No .vtp files found in $(particles_dir)")
sort!(vtp_files; by=step_index)

# Pull NumberOfPoints out of the XML header. The file tail is raw binary, so read
# only the leading chunk and regex-match (equivalent to grep -a -m1).
function read_num_points(path)
    io = open(path, "r")
    try
        header = String(read(io, 2048))
        m = match(r"NumberOfPoints=\"(\d+)\"", header)
        m === nothing && error("Could not find NumberOfPoints in $(path)")
        return parse(Int, m.captures[1])
    finally
        close(io)
    end
end

steps     = step_index.(vtp_files)
merge_on  = read_num_points.(vtp_files)

# Per-step shed count = n_shed_locations * pps, taken from the first nonzero
# merged count (the first shedding event, before merging removes anything).
first_nonzero_idx = findfirst(>(0), merge_on)
first_nonzero_idx === nothing && error("All merged particle counts are zero; nothing to compare.")
shed_per_step    = merge_on[first_nonzero_idx]
n_shed_locations = shed_per_step / pps
s0               = steps[first_nonzero_idx]   # step at which shedding becomes visible

println("Run directory:      $(run_dir)")
println("Particle files:     $(length(vtp_files)) (steps $(first(steps))..$(last(steps)))")
println("p_per_step:         $(pps)")
println("shed_per_step:      $(shed_per_step)  (first nonzero merged count, step $(s0))")
println("n_shed_locations:   $(n_shed_locations)  (both blades)")

# No-merge baseline: pure accumulation, onset aligned to the data at step s0.
function merge_off_count(step)
    n_events = step - s0 + 1
    return n_events <= 0 ? 0 : shed_per_step * n_events
end
merge_off = merge_off_count.(steps)

out_csv = get(ENV, "OUT_CSV", joinpath(run_dir, "particles_merging_savings.csv"))
open(out_csv, "w") do io
    println(io, "step,revolution,particles_merge_on,particles_merge_off,savings_ratio")
    for i in eachindex(steps)
        rev   = steps[i] * dt * RPM / 60
        ratio = merge_off[i] / max(merge_on[i], 1)
        println(io, "$(steps[i]),$(rev),$(merge_on[i]),$(merge_off[i]),$(ratio)")
    end
end

final_ratio = merge_off[end] / max(merge_on[end], 1)
println("\nFinal step $(last(steps)): merge_on=$(merge_on[end]), merge_off=$(merge_off[end]), " *
        "ratio=$(round(final_ratio, digits=3))")
println("Wrote CSV: $(out_csv)")
