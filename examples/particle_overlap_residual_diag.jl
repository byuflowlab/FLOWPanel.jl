# =============================================================================
# Driver: particle-overlap vorticity-evolution residual diagnostic (item 008)
#
# Offline, read-only over SAVED particle VTP states. Runs validation gates, then
# the per-state collocation mat-vec residual diagnostic, writes a per-state CSV,
# and prints a go/no-go summary.
#
# Usage:
#   julia --project examples/particle_overlap_residual_diag.jl [RUN_NAME] [STEPS]
#   ENV overrides: RUN_NAME, STEPS (e.g. "320:359" or "350"), CUTOFF_FACTOR,
#                  NSTATES (last-N fallback when STEPS unset), TRANSPOSED.
#
# Defaults to RUN_NAME=rotor_hover_pressure_comparison, last NSTATES=40 steps.
# =============================================================================

include(joinpath(@__DIR__, "particle_overlap_residual.jl"))
const POR = ParticleOverlapResidual

import Printf: @printf, @sprintf
import Statistics: mean

# ---------------------------------------------------------------------------
# Resolve run / step selection
# ---------------------------------------------------------------------------
run_name      = length(ARGS) >= 1 ? ARGS[1] : get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison")
steps_arg     = length(ARGS) >= 2 ? ARGS[2] : get(ENV, "STEPS", "")
cutoff_factor = parse(Float64, get(ENV, "CUTOFF_FACTOR", "4.0"))
nstates       = parse(Int,     get(ENV, "NSTATES", "40"))
transposed    = parse(Bool,    get(ENV, "TRANSPOSED", "true"))

data_root  = joinpath(@__DIR__, "..", "data")
part_dir   = joinpath(data_root, run_name, "$(run_name)_wake1_particles")
isdir(part_dir) || error("Particle dir not found: $(part_dir)")

vtp_path(idx) = joinpath(part_dir, "$(run_name)_wake1_particles.$(idx).vtp")

# discover available step indices
avail = Int[]
for f in readdir(part_dir)
    m = match(Regex("$(run_name)_wake1_particles\\.(\\d+)\\.vtp\$"), f)
    m === nothing || push!(avail, parse(Int, m.captures[1]))
end
sort!(avail)
isempty(avail) && error("No particle VTP files in $(part_dir)")

steps =
    if !isempty(steps_arg)
        if occursin(':', steps_arg)
            lo, hi = parse.(Int, split(steps_arg, ':'))
            filter(s -> lo <= s <= hi, avail)
        else
            [parse(Int, steps_arg)]
        end
    else
        avail[max(1, end-nstates+1):end]
    end
isempty(steps) && error("No steps selected (requested '$(steps_arg)' from $(length(avail)) available).")

println("="^78)
println("Particle-overlap residual diagnostic (item 008)")
println("  run         = $(run_name)")
println("  steps       = $(first(steps))…$(last(steps))  ($(length(steps)) states; $(length(avail)) available)")
println("  cutoff_fac  = $(cutoff_factor)   transposed = $(transposed)")
println("="^78)

# ---------------------------------------------------------------------------
# Validation gates
# ---------------------------------------------------------------------------
println("\n[gate 1] no-overlap limit (isolated particles ⇒ collocation = isolated update)")
let
    # two well-separated particles; residual must be ~machine-eps
    X = [0.0 10.0; 0.0 0.0; 0.0 0.0]
    gamma = [0.3 -0.2; 0.5 0.7; -0.4 0.1]
    sigma = [0.1, 0.13]
    J = repeat([0.2, -0.1, 0.05, 0.3, -0.2, 0.1, -0.15, 0.25, 0.05], 1, 2)
    vel = zeros(3,2); vort = POR.overlap_matvec_bruteforce((; X, sigma, np=2), gamma)
    st = (; X, gamma, sigma, velocity=vel, vorticity=vort, J, np=2)
    res = POR.analyze(st; cutoff_factor, transposed)
    @printf("    r_red_phys=%.2e  r_red_fvpm=%.2e  (expect ~0)\n", res.r_red_phys, res.r_red_fvpm)
    @assert res.r_red_phys < 1e-10 && res.r_red_fvpm < 1e-10 "no-overlap limit failed"
    println("    PASS")
end

println("\n[gate 2] mat-vec (hash-grid) vs brute-force O(N²) reference, on a real state")
let
    st = POR.load_state(vtp_path(steps[end]))
    # subsample for cheap O(N²) reference if large
    nsub = min(st.np, 4000)
    idx = round.(Int, range(1, st.np; length=nsub))
    sub = (; X=st.X[:,idx], gamma=st.gamma[:,idx], sigma=st.sigma[idx],
             velocity=st.velocity[:,idx], vorticity=st.vorticity[:,idx],
             J=st.J[:,idx], np=nsub)
    cutoff = cutoff_factor * maximum(sub.sigma)
    grid = POR.build_grid(sub.X, sub.np, cutoff)
    fast = POR.overlap_matvec(grid, sub, sub.gamma)
    ref  = POR.overlap_matvec_bruteforce(sub, sub.gamma)
    rel  = POR.relnorm(fast .- ref, ref)
    @printf("    np_sub=%d  cutoff=%.3g  rel(matvec - bruteforce)=%.2e  (= Gaussian-tail truncation error)\n", nsub, cutoff, rel)
    # The only difference is the hard cutoff dropping kernel tails; for the
    # Gaussian-erf at 4·σmax this is ~1e-5, far below any meaningful residual.
    @assert rel < 1e-3 "hash-grid mat-vec truncation error unexpectedly large"
    println("    PASS (truncation at $(cutoff_factor)·σmax negligible vs residual signals)")
end

# ---------------------------------------------------------------------------
# Per-state diagnostic + CSV
# ---------------------------------------------------------------------------
csv_path = joinpath(data_root, run_name, "particle_overlap_residual.csv")
header = ["step","np","cutoff","sig_min","sig_max",
          "basis_curl_relerr","basis_savedvort_relerr","saved_vort_norm","conv_rel",
          "r_red_phys","r_red_fvpm","r_smp_phys","r_smp_fvpm",
          "domin_min","domin_med","nbr_mean","nbr_max","elapsed_s"]

rows = Vector{NamedTuple}()
println("\n[diagnostic] per-state residuals")
@printf("  %5s %7s %10s %10s %10s %10s %10s %9s\n",
        "step","np","basisCurl","r_smp_phy","r_smp_fvp","r_red_phy","r_red_fvp","nbr_mean")
for s in steps
    st = POR.load_state(vtp_path(s))
    res = POR.analyze(st; cutoff_factor, transposed)
    res === nothing && (println("  step $s: 0 particles, skipped"); continue)
    push!(rows, (; step=s, res...))
    @printf("  %5d %7d %10.3e %10.3e %10.3e %10.3e %10.3e %9.0f\n",
            s, res.np, res.basis_curl_relerr, res.r_smp_phys, res.r_smp_fvpm,
            res.r_red_phys, res.r_red_fvpm, res.nbr_mean)
end

open(csv_path, "w") do io
    println(io, join(header, ","))
    for r in rows
        println(io, join((r.step, r.np, r.cutoff, r.sig_min, r.sig_max,
                          r.basis_curl_relerr, r.basis_savedvort_relerr, r.saved_vort_norm, r.conv_rel,
                          r.r_red_phys, r.r_red_fvpm, r.r_smp_phys, r.r_smp_fvpm,
                          r.domin_min, r.domin_med, r.nbr_mean, r.nbr_max, r.elapsed), ","))
    end
end

# ---------------------------------------------------------------------------
# Summary verdict
# ---------------------------------------------------------------------------
if !isempty(rows)
    mean_field(f) = mean(getfield.(rows, f))
    println("\n" * "="^78)
    println("SUMMARY over $(length(rows)) states")
    @printf("  basis ω vs curl-of-J rel-err     : mean %.3e  (kernel-normalization calibration)\n", mean_field(:basis_curl_relerr))
    if all(r -> r.saved_vort_norm == 0, rows)
        println("  saved `vorticity` point-field    : EMPTY in all states (calibration uses curl-of-J)")
    end
    @printf("  convective-correction rel size   : mean %.3e\n", mean_field(:conv_rel))
    @printf("  r_sampled,phys (GO/NO-GO)        : mean %.3e\n", mean_field(:r_smp_phys))
    @printf("  r_sampled,fvpm (GO/NO-GO)        : mean %.3e\n", mean_field(:r_smp_fvpm))
    @printf("  r_reduced,phys                   : mean %.3e\n", mean_field(:r_red_phys))
    @printf("  r_reduced,fvpm                   : mean %.3e\n", mean_field(:r_red_fvpm))
    @printf("  diag-dominance (min / median)    : %.3f / %.3f\n", mean_field(:domin_min), mean_field(:domin_med))
    @printf("  neighbors (mean / max)           : %.1f / %.0f\n", mean_field(:nbr_mean), mean_field(:nbr_max))
    println("  CSV: $(csv_path)")
    println("="^78)
end
