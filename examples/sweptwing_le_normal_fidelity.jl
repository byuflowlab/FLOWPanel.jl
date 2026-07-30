# Compares per-panel normals against analytic RAE101 surface normals for both
# mirrored builders. Tests hypothesis D: body_neg's diagonal choice yields
# systematically poorer LE surface approximation than body_pos's.

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
include(joinpath(pnl.examples_path, "sweptwing_te_topology.jl"))  # builds both bodies, solves
import CSV, DataFrames
using GeometricTools: PyPlot as plt

# ---- RAE101 airfoil normals -----------------------------------------------
df = CSV.read(joinpath(pnl.examples_path, "data", "airfoil-rae101.csv"), DataFrames.DataFrame)
xs_raw = collect(df[:, 1]); zs_raw = collect(df[:, 2])
# File goes TE→upper→LE→lower→TE. Split at LE (x=0).
le_idx = findfirst(x -> x == 0.0, xs_raw)
upper_x = reverse(xs_raw[1:le_idx]);  upper_z = reverse(zs_raw[1:le_idx])  # 0..1 ascending
lower_x = xs_raw[le_idx:end];          lower_z = zs_raw[le_idx:end]         # 0..1 ascending

function linterp(xs, ys, x)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x)
    t = (x - xs[i]) / (xs[i+1] - xs[i])
    return (1-t)*ys[i] + t*ys[i+1]
end

# Build dense (x/c, z/c, dz/dx) samples
ns = 2001
xs_dense = collect(range(0.0, 1.0; length=ns))
zu_dense = [linterp(upper_x, upper_z, x) for x in xs_dense]
zl_dense = [linterp(lower_x, lower_z, x) for x in xs_dense]
function central_diff(xs, ys)
    n = length(xs); d = similar(ys)
    d[1] = (ys[2] - ys[1]) / (xs[2] - xs[1])
    d[end] = (ys[end] - ys[end-1]) / (xs[end] - xs[end-1])
    for i in 2:n-1
        d[i] = (ys[i+1] - ys[i-1]) / (xs[i+1] - xs[i-1])
    end
    return d
end
dzu_dense = central_diff(xs_dense, zu_dense)
dzl_dense = central_diff(xs_dense, zl_dense)

# Analytic 3D normal at (x/c, y, side) for this swept wing.
# Surface r(s, y) = (x_LE(y) + s*c, y, c*z_surface(s)), x_LE(y) = |y|*tan(λ)
# ∂r/∂s × ∂r/∂y = (-c*z'(s), c*z'(s)*sign(y)*tan(λ), c)
function analytic_normal(xc::Float64, y::Float64, side::Symbol)
    zprime = side == :upper ? linterp(xs_dense, dzu_dense, xc) : linterp(xs_dense, dzl_dense, xc)
    sgn = sign(y)
    nx = -zprime
    ny = zprime * sgn * tan(lambda * pi / 180)
    nz = 1.0
    nrm = sqrt(nx^2 + ny^2 + nz^2)
    return (nx / nrm, ny / nrm, nz / nrm)
end

c_root = b / ar
function panel_normal_error(body, p)
    xc_world = body.controlpoints[1, p]
    yc_world = body.controlpoints[2, p]
    zc_world = body.controlpoints[3, p]
    nz = body.normals[3, p]
    side = nz > 0 ? :upper : :lower
    x_over_c = (xc_world - abs(yc_world) * tan(lambda * pi / 180)) / c_root
    # clamp to [0,1] (tiny excursions at TE/LE numerical)
    x_over_c = clamp(x_over_c, 0.0, 1.0)
    nan, nay, naz = analytic_normal(x_over_c, yc_world, side)
    # ensure analytic normal points same side as panel (upper vs lower)
    if side == :lower
        nan, nay, naz = -nan, -nay, -naz
    end
    nx, ny, nz_b = body.normals[:, p]
    dotp = nx*nan + ny*nay + nz_b*naz
    dotp = clamp(dotp, -1.0, 1.0)
    return (x_over_c, side, acos(dotp))  # radians
end

function histogram_errors(body, label)
    le_errs = Float64[]; mid_errs = Float64[]; te_errs = Float64[]
    for p in 1:body.ncells
        x, side, err = panel_normal_error(body, p)
        err_deg = rad2deg(err)
        if x < 0.1
            push!(le_errs, err_deg)
        elseif x < 0.7
            push!(mid_errs, err_deg)
        else
            push!(te_errs, err_deg)
        end
    end
    println("  $label   LE (x/c<0.1): n=$(length(le_errs))  RMS=$(round(sqrt(sum(abs2,le_errs)/length(le_errs)),sigdigits=4))°  max=$(round(maximum(le_errs),sigdigits=4))°")
    println("  $label   MID (0.1≤x/c<0.7): n=$(length(mid_errs))  RMS=$(round(sqrt(sum(abs2,mid_errs)/length(mid_errs)),sigdigits=4))°  max=$(round(maximum(mid_errs),sigdigits=4))°")
    println("  $label   TE  (x/c≥0.7): n=$(length(te_errs))  RMS=$(round(sqrt(sum(abs2,te_errs)/length(te_errs)),sigdigits=4))°  max=$(round(maximum(te_errs),sigdigits=4))°")
    return (le=le_errs, mid=mid_errs, te=te_errs)
end

println("\n#===== LE normal fidelity comparison =====#")
errs_pos = histogram_errors(body_pos, "body_pos")
errs_neg = histogram_errors(body_neg, "body_neg")

# ---- vs-x/c scatter plot ---------------------------------------------------
xs_pos = Float64[]; es_pos = Float64[]
xs_neg = Float64[]; es_neg = Float64[]
for p in 1:body_pos.ncells
    x, _, e = panel_normal_error(body_pos, p)
    push!(xs_pos, x); push!(es_pos, rad2deg(e))
end
for p in 1:body_neg.ncells
    x, _, e = panel_normal_error(body_neg, p)
    push!(xs_neg, x); push!(es_neg, rad2deg(e))
end

fig = plt.figure(figsize=(8, 4.5))
ax = fig.subplots(1, 1)
ax.scatter(xs_pos, es_pos, s=4, alpha=0.4, label="body_pos (+y lofter)")
ax.scatter(xs_neg, es_neg, s=4, alpha=0.4, label="body_neg (−y lofter)")
ax.set_xlabel("x/c")
ax.set_ylabel("|panel normal − analytic normal|  (deg)")
ax.set_title("Surface-normal fidelity, n_rfl=$n_rfl, n_span=$n_span")
ax.set_yscale("log")
ax.legend()
ax.grid(true, alpha=0.3)
fig.tight_layout()
fig.savefig(joinpath(@__DIR__, "..", "sweptwing_normal_fidelity.png"), dpi=150)
println("Saved normal-fidelity plot to sweptwing_normal_fidelity.png")
