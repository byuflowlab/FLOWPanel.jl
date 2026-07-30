## Rotor hover drag-polar replay + inflow-direction comparison.
##
## Replays saved VTK output from a rotor_hover run and adds sectional viscous
## (profile-drag) forces from an airfoil drag polar via DragPolarMonitor. The
## drag monitor determines the effective sectional inflow direction two ways —
## an area-weighted average of the surface velocity, and an average over a ring
## of off-body probe points around the quarter-chord — and this script reports
## how well the two agree, together with the sectional drag distribution and the
## viscous correction to CT/CQ.
##
## Switch the source run with DRAGPOLAR_DATA_PATH (and DRAGPOLAR_RUN_NAME if the
## run name differs from the directory basename). Output is written to
## data/rotor_hover_dragpolar_replay/ (overwritten in place each run).
##
## Key env overrides:
##   DRAGPOLAR_DATA_PATH   source run dir (default data/rotor_hover_pressure_comparison)
##   DRAGPOLAR_RUN_NAME    run name       (default basename of DATA_PATH)
##   DRAGPOLAR_NSTEPS      #steps to replay from the end (default 36 = 1 rev)
##   DRAGPOLAR_NBINS       spanwise bins  (default 30)
##   DRAGPOLAR_NCRIT       XFOIL ncrit for the polar (default 4)
##   DRAGPOLAR_PROBE_RADIUS probe-circle radius in chords (default 1.0)
##   DRAGPOLAR_INFLOW      :both | :surface | :probe (default both)

import FLOWPanel as pnl
import CSV
import DataFrames
using DataFrames: DataFrame
using LinearAlgebra: cross, dot, norm
using Printf: @printf, @sprintf
using Statistics: mean, median, std
using FLOWPanel.FastMultipole.StaticArrays: SVector

# optional plotting
const HAVE_PLT = try
    @eval import PythonPlot as plt
    true
catch err
    @warn "PythonPlot unavailable; skipping plots" exception=err
    false
end

# ------------------------------------------------------------------ config ----
env_str(name, default)   = get(ENV, name, default)
env_int(name, default)   = parse(Int, get(ENV, name, string(default)))
env_float(name, default) = parse(Float64, get(ENV, name, string(default)))

DATA_PATH  = env_str("DRAGPOLAR_DATA_PATH", joinpath("data", "rotor_hover_pressure_comparison"))
RUN_NAME   = env_str("DRAGPOLAR_RUN_NAME", basename(DATA_PATH))
N_STEPS    = env_int("DRAGPOLAR_NSTEPS", 36)
NBINS      = env_int("DRAGPOLAR_NBINS", 30)
NCRIT      = env_int("DRAGPOLAR_NCRIT", 4)
PROBE_RADIUS_CHORDS = env_float("DRAGPOLAR_PROBE_RADIUS", 1.0)
INFLOW_METHOD = Symbol(env_str("DRAGPOLAR_INFLOW", "both"))
SAVE_PATH  = env_str("DRAGPOLAR_SAVE_PATH", joinpath("data", "rotor_hover_dragpolar_replay"))
CCBLADE_DIR = env_str("DRAGPOLAR_CCBLADE_DIR", joinpath("data", "rotor_hover_ccblade"))

rho = env_float("DRAGPOLAR_RHO", 1.179)
R   = env_float("DRAGPOLAR_R", 0.119)
rpm = env_float("DRAGPOLAR_RPM", 6000.0)
axial_dimension  = env_int("DRAGPOLAR_AXIAL_DIM", 1)
radial_dimension = env_int("DRAGPOLAR_RADIAL_DIM", 2)
hub_cut_rR = env_float("DRAGPOLAR_HUBCUT", 0.15)

n = rpm / 60          # rev/s
D = 2 * R

# reset output dir
isdir(SAVE_PATH) && rm(SAVE_PATH; recursive=true, force=true)
mkpath(SAVE_PATH)

unit_axis(dim::Int) = SVector{3,Float64}(ntuple(i -> i == dim ? 1.0 : 0.0, 3))

# --------------------------------------------------------------- step select --
function last_n_steps(data_path, run_name, nsteps)
    _, idxs = pnl._read_pvd_steps(joinpath(data_path, run_name * "_body1.pvd"))
    isempty(idxs) && error("No saved steps found for $(run_name) in $(data_path)")
    n = min(nsteps, length(idxs))
    return idxs[(end - n + 1):end]
end

# ---------------------------------------------------------- polar + chord ------
# Build cd(cl) per section from the CCBlade polars CSV (pre-stall increasing-cl
# branch, linearly interpolated and clamped) plus chord(r/R) from the sectional
# CSV. Returns callables polar(cl, span_coord) and chord(span_coord) with the
# span coordinate in meters (frame units).
function build_polar_and_chord(ccblade_dir, ncrit, R)
    polar_csv = joinpath(ccblade_dir, "rotor_hover_ccblade_polars.csv")
    sect_csv  = joinpath(ccblade_dir, "rotor_hover_ccblade_sectional_ncrit$(ncrit).csv")
    isfile(polar_csv) || error("Polars CSV not found: $(polar_csv)")
    isfile(sect_csv)  || error("Sectional CSV not found: $(sect_csv)")

    pol = CSV.read(polar_csv, DataFrame)
    sub = pol[pol.ncrit .== ncrit, :]
    isempty(sub) && error("No rows with ncrit=$(ncrit) in $(polar_csv)")

    # per-section (increasing-cl branch) tables, sorted by r_over_R
    secs = sort(unique(sub.r_over_R))
    sec_cl = Vector{Vector{Float64}}(undef, length(secs))
    sec_cd = Vector{Vector{Float64}}(undef, length(secs))
    for (i, rr) in enumerate(secs)
        s = sub[sub.r_over_R .== rr, :]
        order = sortperm(s.alpha_deg)
        alpha = s.alpha_deg[order]; cl = s.cl[order]; cd = s.cd[order]
        # restrict to the branch from cl_min (lowest alpha) up to cl_max (stall)
        imax = argmax(cl)
        cl_b = cl[1:imax]; cd_b = cd[1:imax]
        # enforce strictly increasing cl for interpolation
        keep = [true; diff(cl_b) .> 0]
        sec_cl[i] = cl_b[keep]; sec_cd[i] = cd_b[keep]
    end

    lininterp(xs, ys, x) = begin
        x <= xs[1]  && return ys[1]
        x >= xs[end] && return ys[end]
        j = searchsortedlast(xs, x)
        t = (x - xs[j]) / (xs[j+1] - xs[j])
        (1 - t) * ys[j] + t * ys[j+1]
    end

    cd_of_section(i, cl) = lininterp(sec_cl[i], sec_cd[i], cl)

    function polar(cl, span_coord)
        rr = abs(span_coord) / R
        # bracketing sections
        if rr <= secs[1];   return cd_of_section(1, cl); end
        if rr >= secs[end]; return cd_of_section(length(secs), cl); end
        j = searchsortedlast(secs, rr)
        w = (rr - secs[j]) / (secs[j+1] - secs[j])
        return (1 - w) * cd_of_section(j, cl) + w * cd_of_section(j+1, cl)
    end

    # chord distribution
    sec = CSV.read(sect_csv, DataFrame)
    order = sortperm(sec.r_over_R)
    rR_c = sec.r_over_R[order]; chord_m = sec.chord_m[order]
    chord(span_coord) = lininterp(rR_c, chord_m, abs(span_coord) / R)

    return polar, chord, secs
end

polar, chord, polar_sections = build_polar_and_chord(CCBLADE_DIR, NCRIT, R)

# ------------------------------------------------------------- monitors --------
backend = pnl.FastMultipoleBackend(;
    expansion_order = env_int("FMM_EXPANSION_ORDER", 8),
    multipole_acceptance = env_float("FMM_ACCEPTANCE", 0.4),
    leaf_size = env_int("FMM_LEAF_SIZE", 20),
)

# handles captured by the factory so we can read results after replay
force_inviscid = Ref{pnl.ForceMonitor}()
force_total    = Ref{pnl.ForceMonitor}()
dragpolar_ref  = Ref{pnl.DragPolarMonitor}()

monitor_factory = (systems, wakes, frames, t_range) -> begin
    nt = length(t_range)
    e_axial = unit_axis(axial_dimension)
    thrust_dir = -e_axial                                   # thrust = -F_axial
    tan_dir = cross(e_axial, unit_axis(radial_dimension))
    tan_dir = tan_dir / norm(tan_dir)                       # blade-1 tangential
    select_blade1 = cp_frame -> cp_frame[radial_dimension] > hub_cut_rR * R

    pressure = pnl.PressureBernoulli(rho; unsteady=false)
    force = pnl.ForceMonitor(nt, 1; i_frame=1, rho,
        normalization = pnl.NoNormalization(), file=false, vtk_fields=())
    spanwise = pnl.SpanwiseLoadingMonitor(NBINS, 1; i_frame=1,
        span_axis = unit_axis(radial_dimension),
        components = (thrust=thrust_dir, tangential=tan_dir),
        per_length = false, binning = :span_overlap,
        normalization = pnl.NoSectionalNormalization(),
        select = select_blade1, file=false, vtk_fields=(), verbose=false)
    dragpolar = pnl.DragPolarMonitor(nt, 1, polar, chord; i_frame=1, rho,
        inflow_method = INFLOW_METHOD, backend = pnl.DirectBackend(),
        n_probe = env_int("DRAGPOLAR_NPROBE", 16),
        probe_radius_chords = PROBE_RADIUS_CHORDS, file=false, verbose=true)
    total = pnl.ForceMonitor(nt, 1; i_frame=1, rho, source=:context_force,
        normalization = pnl.RotorNormalization(rho, D, 1), file=false, vtk_fields=())

    force_inviscid[] = force
    force_total[]    = total
    dragpolar_ref[]  = dragpolar
    return (pressure, force, spanwise, dragpolar, total)
end

# ------------------------------------------------------------- run replay ------
steps = last_n_steps(DATA_PATH, RUN_NAME, N_STEPS)
println("Rotor hover drag-polar replay")
println("  data:      $(DATA_PATH)")
println("  run_name:  $(RUN_NAME)")
println("  steps:     $(first(steps)):$(last(steps)) ($(length(steps)) samples)")
println("  rho/R/RPM: $(rho) / $(R) / $(rpm)")
println("  inflow:    $(INFLOW_METHOD), probe radius = $(PROBE_RADIUS_CHORDS) chords")
println("  output:    $(SAVE_PATH)")

result = pnl.replay(DATA_PATH, RUN_NAME;
    monitor_factory,
    recompute = (:auto,),
    steps,
    backend,
    backend_wake = backend,
    backend_system = backend,
    verbose = true)

dp = dragpolar_ref[]
ncols = length(steps)                     # replay writes columns 1:ncols
cols = 1:ncols

# ------------------------------------------------- inflow-direction comparison -
# angle_diff is nbins × nt; only the first `ncols` columns were populated.
angle = dp.inflow_angle_diff[:, cols]
bin_r = dp.bin_center_history ./ R        # r/R per bin (last step; bins fixed)

mean_deg = [mean(filter(!isnan, angle[b, :])) for b in 1:dp.nbins]
std_deg  = [let v = filter(!isnan, angle[b, :]); length(v) > 1 ? std(v) : 0.0 end for b in 1:dp.nbins]
max_deg  = [let v = filter(!isnan, angle[b, :]); isempty(v) ? NaN : maximum(v) end for b in 1:dp.nbins]

inflow_df = DataFrame(r_over_R=bin_r, mean_deg=mean_deg, std_deg=std_deg, max_deg=max_deg)
CSV.write(joinpath(SAVE_PATH, "inflow_angle_diff.csv"), inflow_df)

allang = filter(!isnan, vec(angle))
if isempty(allang)
    println("\nInflow comparison: no probe/surface pairs (inflow_method=$(INFLOW_METHOD))")
else
    worst_b = argmax(replace(mean_deg, NaN => -Inf))
    @printf("\nInflow direction comparison (surface vs probe circle):\n")
    @printf("  overall  mean = %.3f deg   median = %.3f deg   max = %.3f deg\n",
            mean(allang), median(allang), maximum(allang))
    @printf("  largest per-bin mean disagreement: %.3f deg at r/R = %.3f\n",
            mean_deg[worst_b], bin_r[worst_b])
end

# ------------------------------------------------------- sectional drag / cl,cd
last_col = ncols
sectional_df = DataFrame(
    r_over_R = bin_r,
    cl = dp.cl_history[:, last_col],
    cd = dp.cd_history[:, last_col],
    Dx = dp.drag_force[1, :],
    Dy = dp.drag_force[2, :],
    Dz = dp.drag_force[3, :],
    n_probe_skipped = dp.n_probe_skipped,
)
CSV.write(joinpath(SAVE_PATH, "sectional_drag.csv"), sectional_df)

# ------------------------------------------------------------- force totals ----
# inviscid ForceMonitor uses NoNormalization → dimensional; normalize by rho n^2 D^4.
ct_ref = rho * n^2 * D^4
cq_ref = rho * n^2 * D^5
CT_inv = -force_inviscid[].force[axial_dimension, cols] ./ ct_ref
# tangential moment about the axial axis → torque; RotorNormalization already
# divided total ForceMonitor by rho n^2 D^4 (CF) and D^5 (CM).
CT_tot = -force_total[].force[axial_dimension, cols]
CQ_tot =  abs.(force_total[].moment[axial_dimension, cols])

CT_inviscid = mean(CT_inv)
CT_total    = mean(CT_tot)
# viscous-only totals (dimensional, frame coords) as a cross-check
Fvisc = dp.force[:, cols]
Mvisc = dp.moment[:, cols]
dCT_visc = mean(-Fvisc[axial_dimension, :]) / ct_ref
dCQ_visc = mean(abs.(Mvisc[axial_dimension, :])) / cq_ref

@printf("\nForce totals (averaged over %d replayed steps):\n", ncols)
@printf("  inviscid  CT = %.5f\n", CT_inviscid)
@printf("  corrected CT = %.5f   (ΔCT = %+.5f)\n", CT_total, CT_total - CT_inviscid)
@printf("  corrected CQ = %.6f\n", mean(CQ_tot))
@printf("  viscous-only ΔCT = %+.5f   ΔCQ = %+.6f  (cross-check)\n", dCT_visc, dCQ_visc)

# --------------------------------------------------------------- plots ---------
if HAVE_PLT
    # inflow angle vs r/R
    if !isempty(allang)
        fig, ax = plt.subplots(figsize=(7.0, 4.5))
        ax.plot(bin_r, mean_deg, "-o"; color="tab:blue", markersize=3, label="mean")
        ax.fill_between(bin_r, mean_deg .- std_deg, mean_deg .+ std_deg;
                        color="tab:blue", alpha=0.25, label="±1σ")
        ax.set_xlabel("r/R"); ax.set_ylabel("inflow direction difference (deg)")
        ax.set_title("Surface-average vs probe-circle inflow"); ax.legend(); ax.grid(true, alpha=0.3)
        fig.tight_layout(); fig.savefig(joinpath(SAVE_PATH, "inflow_angle_vs_rR.png"), dpi=150)
        plt.close()
    end

    # sectional cl, cd, dD/dr
    dDdr = [norm(dp.drag_force[:, b]) for b in 1:dp.nbins]  # bin_width baked into drag_force
    fig, axs = plt.subplots(3, 1; figsize=(7.0, 9.0), sharex=true)
    axs[0].plot(bin_r, dp.cl_history[:, last_col], "-o"; markersize=3, color="tab:green")
    axs[0].set_ylabel("cl"); axs[0].grid(true, alpha=0.3)
    axs[1].plot(bin_r, dp.cd_history[:, last_col], "-o"; markersize=3, color="tab:red")
    axs[1].set_ylabel("cd"); axs[1].grid(true, alpha=0.3)
    axs[2].plot(bin_r, dDdr, "-o"; markersize=3, color="tab:purple")
    axs[2].set_ylabel("|drag per bin| (N)"); axs[2].set_xlabel("r/R"); axs[2].grid(true, alpha=0.3)
    axs[0].set_title("Sectional loading (last replayed step)")
    fig.tight_layout(); fig.savefig(joinpath(SAVE_PATH, "sectional_drag.png"), dpi=150)
    plt.close()
    println("\nPlots written to $(SAVE_PATH)")
end

println("\nDone. CSVs and plots in $(SAVE_PATH)")
