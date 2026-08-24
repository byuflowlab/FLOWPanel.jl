#=##############################################################################
# fm051_attribution_debug.jl
#
# Small-problem debug/validation for the pass-3 attribution stage (see
# fm051_pass3_attribution.jl): builds the mini diamond RigidWakeBody at a
# larger span count with a FORCED-farfield fmm backend (small leaf), sets
# deterministic panel strengths, and evaluates the body self-influence pass
# (`_sa_body_influence!`, body-only targets, production conditioning) three
# ways:
#
#   exact   DirectBackend()               -- brute-force reference
#   cpu     FastMultipoleBackend(P,...)   -- for each P in FM051_DEBUG_PS
#   seam    :host rect functor + fmm backend (seam intercepts)
#
# Expected on a healthy tree: seam vs exact at F64 summation roundoff
# (diff_to_scale <~ 1e-13); cpu-fmm vs exact error nonzero and DECREASING with
# P (truncation). Also times the exact arm and extrapolates the per-pair cost
# to the p018 size (36752 panels) so the full-run stage cost is known before
# submission.
#
# Run from the FLOWPanel.jl root (4 threads max per local-run rules):
#   FM051_DEBUG_NSPAN=40 julia --project=. --threads=4 benchmark/fm051_attribution_debug.jl
=###############################################################################

import Printf
using Printf
import FLOWPanel as pnl
import FastMultipole

include(joinpath(@__DIR__, "..", "test", "test_helpers.jl"))
include(joinpath(@__DIR__, "fm051_pass3_attribution.jl"))

const NSPAN = parse(Int, get(ENV, "FM051_DEBUG_NSPAN", "40"))
const NCHORD = parse(Int, get(ENV, "FM051_DEBUG_NCHORD", "20"))
const PS = [parse(Int, s) for s in split(get(ENV, "FM051_DEBUG_PS", "4,10"), ",")]
const MAC = parse(Float64, get(ENV, "FM051_DEBUG_MAC", "0.5"))
const LEAF = parse(Int, get(ENV, "FM051_DEBUG_LEAF", "16"))
const GATE_SEAM = parse(Float64, get(ENV, "FM051_DEBUG_GATE_SEAM", "1e-10"))
const P018_NPANELS = 36752

println("FM051 pass-3 attribution debug -- nspan=$(NSPAN) nchord=$(NCHORD) " *
        "(=> $(4*NSPAN*NCHORD) panels), fmm P=$(PS) mac=$(MAC) leaf=$(LEAF), " *
        "threads=$(Threads.nthreads())")

# ---- body with deterministic strengths (no solve needed: influence is linear
# in strength, so any smooth deterministic field exercises the same code) ----
#
# Chordwise-REFINED diamond (test_helpers' make_dirichlet_diamond_body has only
# 4 chordwise cells, so its panels span half the chord; the source octree then
# cannot subdivide below panel-radius scale and collapses to ONE leaf -- no m2l
# pairs exist and the "fmm" arm is trivially exact (probe_lists.jl evidence:
# src tree n_branches=1 at every MAC/method). Small panels restore a deep tree
# and a real farfield, which is the whole point of this debug case.)
function make_refined_diamond_body(; nspan::Int, nchord::Int, thick=0.06, das=0.15)
    ys = range(0, 1; length=nspan + 1)
    nperim = 2 * nchord                    # closed perimeter loop per station
    zup(x) = x <= 0.5 ? thick * 2x : thick * (2 - 2x)
    nodes = Float64[]
    for y in ys, k in 0:nperim-1
        if k <= nchord
            x = k / nchord                 # upper surface LE -> TE
            append!(nodes, [x, y, zup(x)])
        else
            x = (nperim - k) / nchord      # lower surface TE -> LE
            append!(nodes, [x, y, -zup(x)])
        end
    end
    nodes = reshape(nodes, 3, :)
    nid(j, k) = (j - 1) * nperim + mod(k, nperim) + 1
    cells = Int[]
    for j in 1:nspan, k in 0:nperim-1
        a1, b1 = nid(j, k), nid(j, k + 1)
        a2, b2 = nid(j + 1, k), nid(j + 1, k + 1)
        append!(cells, [a1, b1, b2])
        append!(cells, [a1, b2, a2])
    end
    cells = reshape(cells, 3, :)
    shedding = [pnl.calc_shedding_from_seed(nodes, cells, nid(1, nchord), nid(2, nchord))]
    body = pnl.RigidWakeBody{Union{pnl.ConstantSource,pnl.VortexRing}}(
        nodes, cells, shedding; check_mesh=false, watertight=false,
        semiinfinite_wake=false)
    for i in eachindex(body.Das)
        body.Das[i] .= repeat([das, 0.0, 0.0], 1, size(body.Das[i], 2))
    end
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    return body
end

body = make_refined_diamond_body(; nspan=NSPAN, nchord=NCHORD, thick=0.06, das=0.15)
body.needs_velocity_gradient[] = false
for j in axes(body.strength, 2), i in axes(body.strength, 1)
    body.strength[i, j] = sin(0.7 * i + 1.3 * j) + 0.1
end
# p018-like core-size split: effectively-singular self block, physical targets core
body.core_size_panel = 1e-10
body.core_size_targets = 1e-3
n = body.ncells
@assert n == size(body.velocity, 2)

function run_arm!(body, backend; seam::Symbol=:off)
    result = (0.0, 0, zeros(0, 0))
    for rep in 1:2                          # rep 1 warms up (JIT); rep 2 is timed
        body.velocity .= 0.0
        pnl._set_core_sizes!((body,), :core_size_targets)
        pnl.set_gpu_influence!(seam)
        h0 = pnl.GPU_INFLUENCE_HITS[]
        t = @elapsed pnl._sa_body_influence!((body,), (body,), backend;
            needs_induced_vorticity=false, body_on_wake=true,
            body_hessian_to_particles=false, body_gradient_core_size=NaN)
        hits = pnl.GPU_INFLUENCE_HITS[] - h0
        pnl.set_gpu_influence!(:off)
        result = (t, hits, copy(body.velocity))
    end
    t, hits, vel = result
    return vel, t, hits
end

# ---- arms ----
dexact, t_exact, _ = run_arm!(body, pnl.DirectBackend())
@printf("exact direct: %.3f s for %d^2 pairs (%.3e s/pair)\n",
    t_exact, n, t_exact / n^2)

fmm_backend(P) = pnl.FastMultipoleBackend(P, MAC, LEAF)

dseam, t_seam, hits_seam = run_arm!(body, fmm_backend(PS[end]); seam=:host)
hits_seam > 0 || error("seam (:host) did not accept the pass-3-shaped call -- " *
    "attribution debug is vacuous; check seam eligibility for this body")

ok = true
for P in PS
    dcpu, t_cpu, hits_cpu = run_arm!(body, fmm_backend(P))
    hits_cpu == 0 || error("seam hit during a CPU arm (P=$(P)) -- seam not disarmed")
    r = attr_report(dexact, dcpu, dseam; gate_seam=GATE_SEAM,
        label="debug P=$(P) mac=$(MAC) leaf=$(LEAF)")
    @printf("  walls: exact %.3f s | cpu-fmm %.3f s | seam %.3f s (hits %d)\n",
        t_exact, t_cpu, t_seam, hits_seam)
    global ok &= r.seam_ok
end

# ---- p018 extrapolation (exact arm is O(n^2) with the same per-pair kernel) --
t_p018 = t_exact * (P018_NPANELS / n)^2
@printf("\np018 exact-arm extrapolation: %.1f s at %d threads (%d panels; scales ~1/threads)\n",
    t_p018, Threads.nthreads(), P018_NPANELS)

println(ok ? "\nDEBUG PASS: seam matches exact direct at every P" :
             "\nDEBUG FAIL: seam-vs-exact exceeded the gate (see verdict lines)")
exit(ok ? 0 : 1)
