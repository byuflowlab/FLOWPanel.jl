# BRAINSTORM 016 Phase 3 Stage 4c -- discriminating attribution criterion.
#
# Reproduce the pre-registered attribution evidence. Run from the repository
# root with at most six threads:
#
#     julia --project data/surface_vorticity_conversion_static/generate_attribution.jl
#
# The decision rules, thresholds, and dispositions were committed to
# BRAINSTORM/016_surface_vorticity_particle_shedding/phase_03_implementation.md
# ("Stage 4c -- pre-registration") BEFORE this script was ever run. Nothing here
# time-marches, sheds, solves, or mutates the checked-in evidence of Stage 4b.
#
# Why this exists: Stage 4b did not select `:upstream`; it failed to displace it.
# Its near-field metric scored the hybrid against `u_sheet` -- the edge-jump
# representation this item exists to replace -- so it structurally rewards the
# artefact. Every metric below is measured against a *refined reference* or is
# reference-free.
using LinearAlgebra
using StaticArrays
using Printf
import FLOWPanel as pnl
import FLOWVPM
import FastMultipole
include(joinpath(@__DIR__, "..", "..", "test", "test_helpers.jl"))

const OUT      = @__DIR__
const SIGMA0   = 0.06
const OVERLAP  = 1.3
const NSPAN0   = 12
# Standoffs. Stage 4b stopped at 8; the alpha difference is a dipole layer of
# thickness `h_row`, so d/sigma <= 8 sits INSIDE it (d < h_row) at every
# resolution swept here and cannot show far-field decay. The extra decades are
# what separate "alpha does not vanish in the near field" from "alpha does not
# vanish at all".
const RATIOS   = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0)
const NSTATION = 8                       # T1 pre-registration: >= 8 stations
const ALPHA_OF = (downstream = 0.0, split = 0.5, upstream = 1.0)
const MODES    = (:downstream, :split, :upstream)

#--- The prescribed sheet: identical maps and strength field to Stage 4b's
#    `generate.jl`, so the two bodies of evidence describe the same object.

mapx(s) = 3expm1(1.1s) / expm1(1.1)
mapy(q) = 1.5sinh(0.8q) / sinh(0.8)
mapz(s, q) = 0.10sin(pi * s) * cos(pi * q / 2)
sheet_point(s, q) = SVector(mapx(s), mapy(q), mapz(s, q))

muhat(x, y, z) = 0.7 + 0.35x + 0.08x^2 + 0.18sin(1.2y) + 0.04x * y + 0.1z
muhat(p::SVector{3}) = muhat(p[1], p[2], p[3])
grad_muhat(x, y, z) = SVector(0.35 + 0.16x + 0.04y,
                              0.216cos(1.2y) + 0.04x,
                              0.1)
grad_muhat(p::SVector{3}) = grad_muhat(p[1], p[2], p[3])

"Outward normal of the parametric sheet at `(s,q)`, by central differences."
function sheet_normal(s, q; h = 1e-6)
    ds = (sheet_point(min(s + h, 1.0), q) - sheet_point(max(s - h, 0.0), q))
    dq = (sheet_point(s, min(q + h, 1.0)) - sheet_point(s, max(q - h, -1.0)))
    n = cross(ds, dq)
    return n / norm(n)
end

"""
Deterministic prescribed-sheet `PanelParticleWake` at `nrows` streamwise rows.

The physical sheet, its node map, and the physical `muhat` field are INVARIANT
under `nrows`: node row `r` sits at `s = (r-1)/nrows`, so refining `nrows` shrinks
the row extent without moving or reshaping the sheet, and strength is assigned
from each cell's centroid *position*. That invariance is what makes the T0
collapse sweep a single-variable experiment.
"""
function prescribed_sheet(mode; nrows = 6, nspan = NSPAN0, sigma = SIGMA0)
    body = make_dirichlet_diamond_body(; nspan)
    conversion = pnl.SurfaceVorticityConversion(sigma; overlap = OVERLAP,
                                                attribution = mode)
    wake = pnl.PanelParticleWake(body; nwakerows = nrows, max_particles = 2_000_000,
                                 core_size = sigma, conversion = conversion)
    nodes = wake.panel_wake.nodes[1]
    strength = wake.panel_wake.strength[1]
    nc = size(nodes, 3) - 1
    for r in axes(nodes, 2), c in axes(nodes, 3)
        s = (r - 1) / nrows
        q = 2(c - 1) / nc - 1
        nodes[:, r, c] .= sheet_point(s, q)
    end
    for r in axes(strength, 2), c in axes(strength, 3)
        rr = min(r, size(nodes, 2) - 1)
        p = 0.25 * (SVector{3}(nodes[:, rr, c]) + SVector{3}(nodes[:, rr + 1, c]) +
                    SVector{3}(nodes[:, rr + 1, c + 1]) + SVector{3}(nodes[:, rr, c + 1]))
        # The stored row beyond the active ones is the smooth continuation of the
        # field, not a starting vortex: this fixture is a steady-state sheet.
        r == size(strength, 2) &&
            (p += 0.5 * (SVector{3}(nodes[:, end, c]) - SVector{3}(nodes[:, end - 1, c])))
        strength[1, r, c] = muhat(p)
    end
    wake.panel_wake.nwakes[] = nrows
    # `FastMultipole.get_n_bodies(::FilamentWrapper{<:PanelWake})` returns ZERO
    # unless `overflowed[]` is set (`src/FLOWPanel_wake.jl:2641`), and only
    # `shed_wake!` ever sets it. A statically built fixture therefore emits a
    # FilamentWrapper that contains no bodies, so the retained final filament --
    # the one alpha-dependent source in the whole comparison -- radiates nothing.
    # A live run is `overflowed[] == true` from the first row shift onward, which
    # is the steady state emulated here, so setting it is what makes the fixture
    # represent the real wake. (Stage 4b's `generate.jl` does NOT set it; see the
    # Stage 4c log, finding R4.)
    wake.panel_wake.overflowed[] = true
    return body, wake
end

#--- Field evaluation ------------------------------------------------------

function _probe_system(positions)
    probes = FastMultipole.ProbeSystem(length(positions), Float64)
    for i in eachindex(positions)
        probes.position[i] = positions[i]
        probes.scalar_potential[i] = 0.0
        probes.gradient[i] = zero(SVector{3,Float64})
        probes.hessian[i] = zero(SMatrix{3,3,Float64,9})
    end
    return probes
end

function direct_velocity(positions, sources)
    probes = _probe_system(positions)
    pnl.influence!((probes,), sources, pnl.DirectBackend();
                   scalar_potential = false, gradient = true, hessian = (false,))
    return copy(probes.gradient)
end

"Velocity-gradient tensors. In this package `gradient` is velocity and `hessian`
is the velocity gradient."
function direct_velocity_gradient(positions, sources)
    probes = _probe_system(positions)
    pnl.influence!((probes,), sources, pnl.DirectBackend();
                   scalar_potential = false, gradient = false, hessian = (true,))
    return copy(probes.hessian)
end

#--- Probe set: fixed in PHYSICAL space, independent of `nrows` -------------
#
# Anchored on the sheet's own aft edge (s = 1), which does not move under the
# `nrows` sweep -- unlike the panel/particle partition at s = 1 - 1/nrows, which
# does, and therefore cannot carry an `nrows`-invariant probe. The `:control`
# family sits at mid-chord (s = 0.5), far from the converted row; it is an
# instrument self-check, not evidence: attribution must not move it.

function probe_set(; sigma = SIGMA0)
    positions = SVector{3,Float64}[]
    family = Symbol[]
    ratio_of = Float64[]
    qs = range(-0.96, 0.96; length = NSTATION)
    for (k, q) in enumerate(qs), (s, base) in ((1.0, :aft), (0.5, :control))
        fam = base === :control ? :control :
              k == 1 ? :root : k == NSTATION ? :tip : :aft
        x0 = sheet_point(s, q)
        n = sheet_normal(s, q)
        for ratio in RATIOS, side in (-1, 1)
            push!(positions, x0 + side * ratio * sigma * n)
            push!(family, fam)
            push!(ratio_of, ratio)
        end
    end
    return positions, family, ratio_of
end

"""
Hybrid field: retained panel rows + the alpha-weighted retained filament + the
particles deposited from the outgoing row.

`particle_handoff_active` is pre-set so the conversion takes its STEADY-STATE
branch. That is the branch attribution actually governs, and it is the one
consistent with this fixture, whose stored terminal row is a smooth continuation
of the field rather than a starting vortex. Post-shift is emulated exactly as
Stage 4b's `generate.jl` does (`nwakes[] = N-1`), which makes
`_final_filament_strength` return `-(alpha*muhat_A + (1-alpha)*muhat_G)`.
"""
function hybrid(mode; nrows, nspan = NSPAN0, sigma = SIGMA0, positions)
    body, wake = prescribed_sheet(mode; nrows, nspan, sigma)
    pw = wake.panel_wake
    N = pw.nwakes[]
    # Steady state: the previous conversion already committed this alpha, so the
    # incoming retained filament must carry it too. Without this the incoming aft
    # face would not equal `beta * downstream` and conservation would fail for
    # every mode but `:upstream` -- an artefact of the emulation, not the code.
    pw.particle_handoff_active[] = true
    pw.particle_handoff_weight[] = ALPHA_OF[mode]
    ledger_before = filament_ledger(pw)
    pnl._convert_to_particles!(wake, body)
    pw.nwakes[] = N - 1
    ledger_after = filament_ledger(pw)
    u = direct_velocity(positions, (pnl.get_sources(pw)..., wake.pfield))
    return (; u, wake, pw, ledger_before, ledger_after)
end

"Full prescribed sheet, no conversion: the reference object."
function full_sheet(; nrows, nspan = NSPAN0, sigma = SIGMA0, positions)
    _, wake = prescribed_sheet(:upstream; nrows, nspan, sigma)
    return direct_velocity(positions, pnl.get_sources(wake.panel_wake))
end

"""
`S = sum_j filament_j * edge_j` over the retained final filament. Every panel
ring's perimeter vector sums to zero, so this is the whole field-relevant
content of the panel wake; the particles' gain must equal its loss. Derived from
`get_sources`, so a wake that emits no `FilamentWrapper` contributes nothing --
the blind spot that hid review finding R2.
"""
function filament_ledger(pw)
    # Asking only whether a FilamentWrapper is *present* is the R2 blind spot in
    # a subtler guise: the wrapper can be present and still emit zero bodies.
    # Require that it actually radiates.
    emits = any(s -> s isa pnl.FilamentWrapper && FastMultipole.get_n_bodies(s) > 0,
                pnl.get_sources(pw))
    emits || return zero(SVector{3,Float64})
    total = zero(SVector{3,Float64})
    N = pw.nwakes[]
    N >= 1 || return total
    for is in eachindex(pw.nodes)
        nd = pw.nodes[is]
        for j in 1:(size(nd, 3) - 1)
            f = pnl._final_filament_strength(pw, is, N, j)
            e = SVector{3}(nd[1, N+1, j+1], nd[2, N+1, j+1], nd[3, N+1, j+1]) -
                SVector{3}(nd[1, N+1, j],   nd[2, N+1, j],   nd[3, N+1, j])
            total += f * e
        end
    end
    return total
end

net_gamma(pf) = SVector{3}(sum(pf.particles[FLOWVPM.GAMMA_INDEX[k], 1:pf.np])
                           for k in 1:3)

#--- Small statistics helpers ---------------------------------------------

function pctile(v, p)
    isempty(v) && return NaN
    o = sort(v)
    return o[clamp(ceil(Int, p * length(o)), 1, length(o))]
end
rms(v) = isempty(v) ? NaN : sqrt(sum(abs2, v) / length(v))
function stdev(v)
    n = length(v)
    n < 2 && return NaN
    m = sum(v) / n
    return sqrt(sum(x -> (x - m)^2, v) / (n - 1))
end
"Least-squares slope of log(y) against log(x): the observed order."
function fit_order(x, y)
    keep = findall(i -> y[i] > 0 && x[i] > 0, eachindex(y))
    length(keep) < 2 && return NaN
    lx = log.(x[keep]); ly = log.(y[keep])
    mx = sum(lx) / length(lx); my = sum(ly) / length(ly)
    return sum((lx .- mx) .* (ly .- my)) / sum(abs2, lx .- mx)
end

writecsv(name, header, rows) = open(joinpath(OUT, name), "w") do io
    println(io, header)
    for r in rows
        println(io, join(r, ','))
    end
end

#=============================================================================
T0 -- is alpha a truncation knob or a modelling choice?
=============================================================================#

function run_T0(positions, family, ratio_of, uref_mag)
    println("\n== T0: alpha spread vs row extent ==")
    nrows_list = (4, 6, 8, 12, 16)
    rows = Any[]
    per_family = Dict{Tuple{Symbol,Float64},Vector{Tuple{Float64,Float64}}}()

    for nrows in nrows_list
        # Physical extent of the converted (final) row, along the sheet midline.
        h_row = mapx(1.0) - mapx(1 - 1 / nrows)
        fields = Dict(m => hybrid(m; nrows, positions).u for m in MODES)
        for fam in (:aft, :root, :tip, :control), ratio in RATIOS
            ids = findall(i -> family[i] === fam && ratio_of[i] == ratio,
                          eachindex(positions))
            isempty(ids) && continue
            spread = [maximum(norm(fields[m1][i] - fields[m2][i])
                              for m1 in MODES, m2 in MODES) for i in ids]
            s_raw = rms(spread)
            s_norm = s_raw / rms([uref_mag[i] for i in ids])
            push!(rows, (nrows, h_row, SIGMA0 / h_row, fam, ratio,
                         ratio * SIGMA0 / h_row, s_raw, s_norm, maximum(spread)))
            push!(get!(per_family, (fam, ratio), Tuple{Float64,Float64}[]),
                  (h_row, s_raw))
        end
    end

    order_rows = Any[]
    for ((fam, ratio), pts) in sort(collect(per_family), by = first)
        p = fit_order(first.(pts), last.(pts))
        push!(order_rows, (fam, ratio, p))
        fam === :aft && @printf("  aft d/sigma=%-5.2f  observed order p = %.3f\n", ratio, p)
    end

    writecsv("attribution_collapse.csv",
             "nrows,h_row,sigma_over_h_row,family,d_over_sigma,d_over_h_row," *
             "S_raw_rms,S_norm_rms,S_raw_max", rows)
    writecsv("attribution_collapse_order.csv", "family,d_over_sigma,observed_order",
             order_rows)
    return order_rows
end

#=============================================================================
T1 -- accuracy against a REFINED reference (the selector)
=============================================================================#

function run_T1(positions, family, ratio_of)
    println("\n== T1: accuracy vs refined reference ==")
    # (a) reference convergence -- the reference must be converged to well below
    #     the alpha differences it is asked to resolve, or T1 is inconclusive.
    # Richardson-extrapolate the two finest rungs: the raw `nrows = 128` field is
    # not converged to below the alpha differences it must resolve, which the
    # pre-registration makes a disqualifying condition for T1.
    fine = (64, 128, 256, 512)
    uref = Dict(n => full_sheet(; nrows = n, positions) for n in fine)
    d21 = [norm(uref[256][i] - uref[128][i]) for i in eachindex(uref[128])]
    d32 = [norm(uref[512][i] - uref[256][i]) for i in eachindex(uref[256])]
    order = log2(rms(d21) / rms(d32))
    gold = uref[512] .+ (uref[512] .- uref[256]) ./ (2^order - 1)
    gold_mag = norm.(gold)
    resid = rms([norm(gold[i] - uref[512][i]) for i in eachindex(gold)])
    conv_rows = Any[]
    for (a, b) in ((64, 128), (128, 256), (256, 512))
        d = [norm(uref[b][i] - uref[a][i]) for i in eachindex(gold)]
        push!(conv_rows, (a, b, rms(d), maximum(d), rms(d) / rms(gold_mag)))
        @printf("  reference nrows %3d -> %3d : rms rel change = %.3e\n",
                a, b, rms(d) / rms(gold_mag))
    end
    @printf("  observed reference order %.2f; Richardson residual = %.3e (rel %.3e)\n",
            order, resid, resid / rms(gold_mag))
    writecsv("attribution_reference_convergence.csv",
             "nrows_coarse,nrows_fine,rms_abs,max_abs,rms_rel", conv_rows)

    # (b) paired accuracy of each mode at the production-like coarse resolution
    nrows_coarse = 6
    E = Dict(m => [norm(hybrid(m; nrows = nrows_coarse, positions).u[i] - gold[i])
                   for i in eachindex(gold)] for m in MODES)
    rows = Any[]
    for fam in (:aft, :root, :tip, :control), ratio in RATIOS
        ids = findall(i -> family[i] === fam && ratio_of[i] == ratio, eachindex(positions))
        isempty(ids) && continue
        paired = [E[:split][i] - E[:upstream][i] for i in ids]
        common = rms([E[:upstream][i] for i in ids])
        push!(rows, (fam, ratio, common,
                     rms([E[:upstream][i] for i in ids]),
                     rms([E[:split][i] for i in ids]),
                     rms([E[:downstream][i] for i in ids]),
                     sum(paired) / length(paired), pctile(paired, 0.5),
                     pctile(paired, 0.95), maximum(paired),
                     count(<(0), paired), length(paired)))
    end
    writecsv("attribution_accuracy.csv",
             "family,d_over_sigma,common_mode_rms,E_upstream_rms,E_split_rms," *
             "E_downstream_rms,paired_mean,paired_median,paired_p95,paired_max," *
             "n_split_better,n_probes", rows)

    aft = findall(i -> family[i] in (:aft, :root, :tip), eachindex(positions))
    paired_all = [E[:split][i] - E[:upstream][i] for i in aft]
    frac = count(<(0), paired_all) / length(paired_all)
    mean_paired = sum(paired_all) / length(paired_all)
    @printf("  over %d sheet-adjacent probes: :split better at %.1f%%, mean paired diff %+.3e\n",
            length(paired_all), 100frac, mean_paired)
    conclusive = resid < 0.1 * rms(abs.(paired_all))
    @printf("  reference residual %.3e vs |paired| rms %.3e -> T1 is %s\n",
            resid, rms(abs.(paired_all)),
            conclusive ? "CONCLUSIVE" : "INCONCLUSIVE (T3 fires)")
    @printf("  selector: %s\n", mean_paired > 0 ? ":upstream (lower error)" :
                                ":split (lower error)")
    return (; gold_mag, frac, paired_all, E, resid, conclusive, mean_paired)
end

#=============================================================================
T2 -- the hazard, judged on the hazard (the veto)
=============================================================================#

function run_T2()
    println("\n== T2: panel-induced stretching at deposited particles ==")
    rows = Any[]
    paired = Dict(:interior => Float64[], :roottip => Float64[], :perimeter => Float64[])

    for sigfac in (0.5, 1.0, 2.0), nrows in (6, 12), nspan in (8, 12, 16)
        sigma = sigfac * SIGMA0
        p95 = Dict{Symbol,Dict{Symbol,Float64}}()
        for mode in (:upstream, :split)
            body, wake = prescribed_sheet(mode; nrows, nspan, sigma)
            pw = wake.panel_wake
            N = pw.nwakes[]
            pw.particle_handoff_active[] = true
            pw.particle_handoff_weight[] = ALPHA_OF[mode]
            pnl._convert_to_particles!(wake, body)
            pw.nwakes[] = N - 1
            pf = wake.pfield
            X = [SVector{3}(pf.particles[FLOWVPM.X_INDEX, i]) for i in 1:pf.np]
            G = [SVector{3}(pf.particles[FLOWVPM.GAMMA_INDEX, i]) for i in 1:pf.np]
            cls = wake.conversion_workspace.classes
            J = direct_velocity_gradient(X, pnl.get_sources(pw))
            buckets = Dict(:interior => Float64[], :roottip => Float64[],
                           :perimeter => Float64[])
            for i in 1:pf.np
                nG = norm(G[i])
                nG > 0 || continue
                key = cls[i] == pnl.InteriorSurfaceParticle ? :interior :
                      cls[i] == pnl.RootTipSurfaceParticle ? :roottip : :perimeter
                push!(buckets[key], norm(J[i] * G[i]) / nG)
            end
            p95[mode] = Dict(k => pctile(v, 0.95) for (k, v) in buckets)
            for (k, v) in buckets
                push!(rows, (sigfac, nrows, nspan, mode, k, length(v),
                             isempty(v) ? NaN : sum(v) / length(v),
                             pctile(v, 0.5), pctile(v, 0.95),
                             isempty(v) ? NaN : maximum(v)))
            end
        end
        for k in keys(paired)
            a = get(p95[:split], k, NaN); b = get(p95[:upstream], k, NaN)
            isfinite(a) && isfinite(b) && push!(paired[k], a - b)
        end
    end

    writecsv("attribution_hazard.csv",
             "sigma_factor,nrows,nspan,mode,class,count,mean,median,p95,max", rows)

    verdict = Any[]
    for k in (:interior, :roottip, :perimeter)
        d = paired[k]
        isempty(d) && continue
        m = sum(d) / length(d); s = stdev(d)
        rejected = m > 0 && abs(m) > 2s
        @printf("  %-10s paired p95 diff (split - upstream): mean %+.4e, sd %.4e -> %s\n",
                k, m, s, rejected ? "REJECT :split" : "no resolved increase")
        push!(verdict, (k, length(d), m, s, m / s, rejected))
    end
    writecsv("attribution_hazard_verdict.csv",
             "class,n_replicates,paired_mean,paired_sd,ratio,split_rejected", verdict)
    return verdict
end

#=============================================================================
alpha* -- diagnostic only, no new source mode
=============================================================================#

function run_alpha_star(; nrows = 6, nspan = NSPAN0, sigma = SIGMA0)
    println("\n== alpha*: error-minimizing streamwise weight (diagnostic) ==")
    kap = Dict{Symbol,Vector{SVector{3,Float64}}}()
    local nodes
    for mode in MODES
        body, wake = prescribed_sheet(mode; nrows, nspan, sigma)
        pw = wake.panel_wake
        pw.particle_handoff_active[] = true
        pnl._convert_to_particles!(wake, body)
        rec = wake.conversion_diagnostics.records[end]
        kap[mode] = [p.kappa_conservative for p in rec.panels]
        nodes = pw.nodes[1]
    end
    N = nrows
    rows = Any[]
    affinity = 0.0
    for j in eachindex(kap[:upstream])
        k0 = kap[:downstream][j]; k1 = kap[:upstream][j]
        affinity = max(affinity, norm(kap[:split][j] - 0.5 * (k0 + k1)) /
                                 max(norm(k0), norm(k1), eps()))
        v1 = SVector{3}(nodes[:, N, j]);     v2 = SVector{3}(nodes[:, N+1, j])
        v3 = SVector{3}(nodes[:, N+1, j+1]); v4 = SVector{3}(nodes[:, N, j+1])
        c = 0.25 * (v1 + v2 + v3 + v4)
        n = cross(v2 - v1, v4 - v1); n /= norm(n)
        kexact = -cross(n, grad_muhat(c))
        # kappa(alpha) is affine, so ||kappa(alpha) - kexact||^2 is a quadratic
        # in alpha with an exact minimizer.
        d = k1 - k0
        r0 = k0 - kexact
        alpha_opt = norm(d)^2 > 0 ? -dot(r0, d) / dot(d, d) : NaN
        err(a) = norm(k0 + a * d - kexact) / norm(kexact)
        # Theoretical second-order weight b^2/(a^2+b^2) from the streamwise
        # centroid separations. Upstream neighbour: a = (h_{N-1}+h_N)/2. The
        # downstream sample is this fixture's terminal strength row, which
        # `prescribed_sheet` places exactly half a row extent beyond the aft
        # edge, so b = h_N/2 -- fixture-specific, not a general wake property.
        hm = mapx((N - 1) / N) - mapx((N - 2) / N)
        hN = mapx(1.0) - mapx((N - 1) / N)
        a_sep = (hm + hN) / 2; b_sep = hN / 2
        alpha_theory = b_sep^2 / (a_sep^2 + b_sep^2)
        # Consistency weight: alpha*a + (1-alpha)*b = h_N, i.e. the alpha that
        # makes the leading term reproduce -mu'. On a graded mesh this exceeds 1,
        # so NO admissible alpha in [0,1] is consistent -- the pre-registered
        # theory result, now with a number attached.
        alpha_consistency = (hN - b_sep) / (a_sep - b_sep)
        push!(rows, (j, hN / hm, alpha_opt, alpha_theory, alpha_consistency,
                     err(0.0), err(0.5), err(1.0),
                     isnan(alpha_opt) ? NaN : err(alpha_opt)))
    end
    @printf("  affinity check |kappa(1/2) - mean(kappa(0),kappa(1))| / |kappa| = %.2e\n",
            affinity)
    opts = [r[3] for r in rows if !isnan(r[3])]
    interior = [r[3] for r in rows[2:end-1] if !isnan(r[3])]
    @printf("  error-minimizing alpha over %d panels: median %.3f, range [%.3f, %.3f]\n",
            length(opts), pctile(opts, 0.5), minimum(opts), maximum(opts))
    @printf("  interior-only median %.3f; consistency alpha %.3f; 2nd-order alpha %.3f\n",
            pctile(interior, 0.5), rows[1][5], rows[1][4])
    @printf("  kappa rel err at alpha = 0 / 0.5 / 1 (interior panel 2): %.4f / %.4f / %.4f\n",
            rows[2][6], rows[2][7], rows[2][8])
    writecsv("attribution_kappa_alpha.csv",
             "panel,h_ratio,alpha_optimal,alpha_theory_2nd_order,alpha_consistency," *
             "rel_err_alpha0,rel_err_alpha_half,rel_err_alpha1,rel_err_alpha_opt",
             rows)
    return affinity, opts
end

#=============================================================================
Instrument self-checks -- run BEFORE any conclusion is drawn
=============================================================================#

function self_checks(positions, family)
    println("== instrument self-checks ==")

    # 1. Conservation, derived from `get_sources`, at every alpha.
    worst = 0.0
    for mode in MODES, nrows in (4, 8)
        h = hybrid(mode; nrows, positions)
        lost = h.ledger_before - h.ledger_after
        gained = net_gamma(h.wake.pfield)
        worst = max(worst, norm(lost - gained))
    end
    @printf("  conservation |filament loss - particle gain|      = %.3e\n", worst)
    worst < 1e-10 || error("Stage 4c: conservation self-check failed")

    # 2. The probe set must not move with `nrows` (T0 is void otherwise).
    p1, _, _ = probe_set()
    p2, _, _ = probe_set()
    @printf("  probe positions independent of nrows              = %s\n",
            all(p1 .== p2) ? "yes (analytic)" : "NO")

    # 3. The T2 instrument must be alpha-SENSITIVE. Stage 5's preflight probe is
    #    not (it runs before commit and before the shift, so the retained
    #    filament is still in its alpha-independent legacy form); repeating that
    #    blind spot would make T2 meaningless.
    sens = Float64[]
    for mode in (:upstream, :split)
        body, wake = prescribed_sheet(mode; nrows = 6)
        pw = wake.panel_wake; N = pw.nwakes[]
        pw.particle_handoff_active[] = true
        pnl._convert_to_particles!(wake, body)
        pw.nwakes[] = N - 1
        pf = wake.pfield
        X = [SVector{3}(pf.particles[FLOWVPM.X_INDEX, i]) for i in 1:pf.np]
        J = direct_velocity_gradient(X, pnl.get_sources(pw))
        push!(sens, sum(norm, J) / length(J))
    end
    rel = abs(sens[1] - sens[2]) / max(sens...)
    @printf("  T2 probe alpha-sensitivity (relative)             = %.3e %s\n",
            rel, rel > 1e-6 ? "(sensitive)" : "(BLIND -- unusable)")
    rel > 1e-6 || error("Stage 4c: T2 instrument is alpha-blind")
    return nothing
end

#--- driver ----------------------------------------------------------------

function main()
    positions, family, ratio_of = probe_set()
    @printf("probes: %d (families: aft/root/tip + control at mid-chord)\n",
            length(positions))
    self_checks(positions, family)

    gold = full_sheet(; nrows = 128, positions)
    uref_mag = norm.(gold)

    run_T0(positions, family, ratio_of, uref_mag)
    run_T1(positions, family, ratio_of)
    run_T2()
    run_alpha_star()
    println("\nwrote CSVs to $(OUT)")
end

main()
