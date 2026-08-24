#=##############################################################################
# fm051_pass3_attribution.jl
#
# Three-way attribution of the pass-3 body self-influence deviation observed in
# H200 job 13306588 (body_velocity diff_to_scale 1.235e-2 vs the 1e-3 phase
# gate): the parity harness compares the GPU rect seam (exact brute force)
# against CPU fmm! (approximate), so a gate FAIL cannot distinguish
#   (a) a seam defect                                    from
#   (b) the CPU arm's own fmm truncation error on the self-body farfield.
# Both arms evaluate self-body pairs at core_size_panel (p018: R*1e-10, i.e.
# effectively the singular kernel), so the regularized-vs-singular term is
# nulled and hypothesis (b) predicts: seam == exact direct to roundoff, CPU
# fmm off by the observed ~1.2e-2.
#
# This file provides the comparison/reporting layer. The EXACT arm itself is
# produced by the caller running the same `_sa_body_influence!` pass with
# `pnl.DirectBackend()` and body-only targets (production conditioning
# semantics are preserved: standalone FastMultipole.direct! honors
# direct_conditioning; with body-only targets the SelfPairs rule covers every
# block). Cost: n_panels^2 pair evaluations -- ~1.35e9 for p018's 36752
# panels, tens of seconds multithreaded; see fm051_attribution_debug.jl for a
# small-problem check and a measured per-pair extrapolation.
#
# Used by:
#   - benchmark/fm051_attribution_debug.jl  (standalone small-problem debug)
#   - benchmark/fm051_pass_parity.jl        (full-mode stage after pass 3)
=###############################################################################

"""
    attr_class_stats(ref, tst; gate=1e-3) -> NamedTuple

Identical metric definitions to the parity harness's `class_parity_stats`
(duplicated here so this file is self-contained for the debug driver; keep the
two in sync): per-target relative error with a `1e-6*scale` floor, `scale` =
largest per-target reference norm, `diff_to_scale` = max per-target |diff| /
scale, plus over-gate/divergent counts and a top-5 offender list.
"""
function attr_class_stats(ref::Matrix{Float64}, tst::Matrix{Float64}; gate::Float64=1e-3)
    (size(ref, 2) == 0 || size(ref) != size(tst)) &&
        return (worst_rel=0.0, iw=0, dn_w=0.0, rn_w=0.0, scale=0.0,
                diff_to_scale=0.0, n_over=0, n_div=0, n=0,
                top=NTuple{4,Float64}[])
    n = size(ref, 2)
    scale = 0.0
    @inbounds for i in 1:n
        s = 0.0
        for c in axes(ref, 1)
            s += ref[c, i]^2
        end
        scale = max(scale, sqrt(s))
    end
    floorv = 1e-6 * scale
    rels = Vector{Float64}(undef, n)
    dns = Vector{Float64}(undef, n)
    rns = Vector{Float64}(undef, n)
    dmax = 0.0
    n_over = 0
    n_div = 0
    @inbounds for i in 1:n
        dn = 0.0
        rn = 0.0
        for c in axes(ref, 1)
            dn += (tst[c, i] - ref[c, i])^2
            rn += ref[c, i]^2
        end
        dn = sqrt(dn); rn = sqrt(rn)
        den = max(rn, floorv)
        e = den > 0 ? dn / den : (dn > 0 ? Inf : 0.0)
        rels[i] = e; dns[i] = dn; rns[i] = rn
        dmax = max(dmax, dn)
        e > gate && (n_over += 1)
        e > 0.5 && (n_div += 1)
    end
    order = partialsortperm(rels, 1:min(5, n); rev=true)
    top = [(Float64(i), rels[i], dns[i], rns[i]) for i in order]
    iw = isempty(order) ? 0 : order[1]
    return (worst_rel=(iw == 0 ? 0.0 : rels[iw]), iw=iw,
            dn_w=(iw == 0 ? 0.0 : dns[iw]), rn_w=(iw == 0 ? 0.0 : rns[iw]),
            scale=scale, diff_to_scale=(scale > 0 ? dmax / scale : (dmax > 0 ? Inf : 0.0)),
            n_over=n_over, n_div=n_div, n=n, top=top)
end

function _attr_print_arm(io, name, st; nshow=5)
    Printf.@printf(io, "  %-22s diff/scale = %.3e   worst rel = %.3e (target %d: |diff| %.3e, |ref| %.3e)\n",
        name, st.diff_to_scale, st.worst_rel, st.iw, st.dn_w, st.rn_w)
    Printf.@printf(io, "  %-22s over-1e-3 %d / divergent %d of %d   (class scale %.3e)\n",
        "", st.n_over, st.n_div, st.n, st.scale)
    for (i, r, dn, rn) in st.top[1:min(nshow, length(st.top))]
        r > 1e-3 || break
        Printf.@printf(io, "      offender target %-8d rel %.3e  |diff| %.3e  |ref| %.3e\n",
            Int(i), r, dn, rn)
    end
end

"""
    attr_report(dexact, dcpu, dseam; gate_seam=1e-10, io=stdout, label="pass3_attribution")
        -> (; seam, cpu, cross, seam_ok)

Three-way attribution of one output class (matrices are per-target CONTRIBUTION
columns, all from the same pre-pass state):

  - `seam`  = seam arm vs exact direct  -- the seam-correctness verdict; gated
              at `gate_seam` on diff_to_scale (F64 summation-order roundoff
              headroom; the arms sum the same pair values in different orders).
  - `cpu`   = CPU fmm! arm vs exact direct -- NOT gated: this is the measured
              legacy-path truncation error, reported as a finding.
  - `cross` = CPU vs seam (should reproduce the parity-stage number; printed
              as a consistency check).

Returns the stats plus `seam_ok`. The caller owns recording/gating.
"""
function attr_report(dexact::Matrix{Float64}, dcpu::Matrix{Float64},
        dseam::Matrix{Float64}; gate_seam::Float64=1e-10, io=stdout,
        label::AbstractString="pass3_attribution")
    st_seam = attr_class_stats(dexact, dseam)
    st_cpu = attr_class_stats(dexact, dcpu)
    st_cross = attr_class_stats(dseam, dcpu)
    println(io, "\n--- $(label): body_velocity three-way vs exact direct ---")
    _attr_print_arm(io, "seam vs exact", st_seam)
    _attr_print_arm(io, "cpu-fmm vs exact", st_cpu)
    _attr_print_arm(io, "cpu-fmm vs seam", st_cross; nshow=0)
    seam_ok = st_seam.diff_to_scale <= gate_seam
    verdict = seam_ok ?
        (st_cpu.diff_to_scale > 10 * st_seam.diff_to_scale ?
            "seam matches exact direct; the parity-stage deviation is the CPU arm's fmm truncation error" :
            "seam matches exact direct; CPU fmm error is at the same (small) level -- no attribution needed") :
        "SEAM DOES NOT MATCH EXACT DIRECT -- genuine seam defect, investigate before touching any gate"
    println(io, "  verdict: $(verdict)")
    return (; seam=st_seam, cpu=st_cpu, cross=st_cross, seam_ok)
end
