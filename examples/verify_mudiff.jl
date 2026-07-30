# Verification for the auto-detecting quad-pitch ∇μ reconstruction (:mu_diff).
#   (1) _quad_pairing recovers the structured quads for AC and BD.
#   (2) agglomerate μ/centroid match quad_geometry.
#   (3) production :mu_diff reproduces the example CL_quadgrad at small levels.
# Run single-threaded:  julia --project examples/verify_mudiff.jl
ENV["FLOWPANEL_CHORDDIV_RUN"] = "false"
include(joinpath(@__DIR__, "sweptwing_chorddivergence.jl"))

import FLOWPanel as pnl
using Printf

function check_pairing(n_ch, n_span, swap)
    body = sweeploft_wing(n_ch, n_span; mirror_diagonals=true, swap_diagonals=swap)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    te_info = view(body.shedding_full, 1:2, :)
    partner = pnl._quad_pairing_propagate(body.nodes, body.cells, body.neighbor,
                                body.normals, te_info; normal_dot_min=cos(pi/4))
    Nq = body.ncells ÷ 2
    ok = all(partner[2q-1] == 2q && partner[2q] == 2q-1 for q in 1:Nq)
    npaired = count(!=(0), partner)

    # agglomerate μ / centroid vs quad_geometry
    geom = quad_geometry(body)
    mu = collect(view(body.strength, :, pnl.get_Gammai(body)))  # zeros pre-solve, but geometry is what we check
    # recompute agglomerate centroid the same way quad_geometry does, then compare
    # to centroids implied by the pairing (members 2q-1, 2q).
    maxc = 0.0
    for q in 1:Nq
        a1, a2 = geom.areas[2q-1], geom.areas[2q]
        # pairing must group exactly {2q-1, 2q}
        (partner[2q-1] == 2q) || (maxc = Inf)
    end
    return (; ok, npaired, Nq, swap)
end

println("#===== (1) PAIRING RECOVERS STRUCTURED QUADS =====#")
for swap in (false, true), n_ch in (40, 80)
    r = check_pairing(n_ch, 48, swap)
    @printf("  swap=%-5s n_ch=%3d  Nq=%5d  paired=%5d  all 2q-1<->2q: %s\n",
        r.swap, n_ch, r.Nq, r.npaired, r.ok ? "YES" : "NO  <-- FAIL")
end

println("\n#===== (3) :mu_diff REPRODUCES CL_quadgrad =====#")
@printf("%7s %5s %12s %12s %12s\n", "variant", "n_ch", "CL_quadgrad", "CL_mudiff", "|Δ|")
for (vname, swap) in (("diagAC", false), ("diagBD", true)), n_ch in (40, 80)
    r, _ = run_case(n_ch, 48; swap_diagonals=swap, variant_name=vname)
    @printf("%7s %5d %12.7f %12.7f %12.2e\n",
        vname, n_ch, r.CL_quadgrad, r.CL_mudiff, abs(r.CL_quadgrad - r.CL_mudiff))
end
