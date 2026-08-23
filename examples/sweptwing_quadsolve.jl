#=##############################################################################
# DESCRIPTION
#   Quad-collocation vortex-ring solve for the sweep-lofted Weber wing
#   ("channel 2" of the triangulation-sensitivity investigation).
#
#   The triangle-based solves converge to two different CL depending on the
#   quad-diagonal handedness (diagAC 0.2678 / diagBD 0.2459 at 80x112) even
#   though the node set is identical. Here the triangulation is removed from
#   the SOLVE itself: one vortex-ring element per quad (4 edges), one Neumann
#   BC collocation per quad centroid, semi-infinite rigid-wake legs on the TE
#   rows — exactly mirroring the conventions of the triangle Backslash solver
#   (G[i,j] = u_j(cp_i)·n̂_i at unit γ, rhs = -U∞·n̂, core_size 1e-10).
#
#   Everything here is direct N-body (no FastMultipole.fmm!), so the script is
#   safe and intended to run threaded:  julia -t 4 --project examples/sweptwing_quadsolve.jl
#
#   ENV knobs:
#     FLOWPANEL_QUADSOLVE_CASE="80:112"   n_ch:n_span of the main case
#     FLOWPANEL_QUADSOLVE_BUDGET_S=3600   abort main case if estimate exceeds this
=###############################################################################

ENV["FLOWPANEL_SWEEPLOFT_RUN"] = "false"
include(joinpath(@__DIR__, "sweptwing_sweeploft.jl"))

import StaticArrays: SVector
import Statistics: mean

const SWITCH = pnl.FastMultipole.DerivativesSwitch(false, true, false)
const CORE = core_size                  # 1e-10, matches the triangle study

# ----------------- QUAD VORTEX-RING KERNEL ------------------------------------
"""
Velocity induced at `target` by a unit-strength vortex ring over the 4 nodes
(in winding order) — identical math to the triangle VortexRing kernel
(`pnl._bound_vortex_velocity` per edge, Vatistas n=2 core), one more edge.
"""
@inline function vortex_ring_quad(target::SVector{3,Float64}, v1, v2, v3, v4,
                                  core::Float64)
    u = pnl._bound_vortex_velocity(v1 - target, v2 - target, true, core)
    u += pnl._bound_vortex_velocity(v2 - target, v3 - target, true, core)
    u += pnl._bound_vortex_velocity(v3 - target, v4 - target, true, core)
    u += pnl._bound_vortex_velocity(v4 - target, v1 - target, true, core)
    return u
end

"Semi-infinite wake legs of a TE quad at unit strength."
@inline function wake_velocity(target::SVector{3,Float64}, te1, te2, Da,
                               core::Float64)
    _, v, _ = pnl.induced_semiinfinite(target, pnl.VortexRing,
        te1[1], te1[2], te1[3], te2[1], te2[2], te2[3],
        Da[1], Da[2], Da[3], 1.0, SWITCH; core_size=core)
    return v
end

# ----------------- QUAD MESH DATA ----------------------------------------------
"""
Build all per-quad data for the solve from the sweep-loft node array:
vertex tuples (winding order A,D,C,B), control points (4-node mean), normals
(diagonal cross — triangulation-free), areas (mean of both diagonal splits),
and TE wake edges (TE1,TE2) per TE quad (reverse of ring traversal, matching
the triangle shedding convention).
"""
function quad_solver_data(nodes, n_ch, n_span; wake_reverse::Bool=true)
    quads = quad_mesh(n_ch, n_span)
    Nq = size(quads, 2)
    npts = 2*n_ch

    node(i) = SVector{3,Float64}(nodes[1, i], nodes[2, i], nodes[3, i])
    verts = Vector{NTuple{4,SVector{3,Float64}}}(undef, Nq)
    cps = Vector{SVector{3,Float64}}(undef, Nq)
    normals = Vector{SVector{3,Float64}}(undef, Nq)
    areas = zeros(Nq)
    for q in 1:Nq
        A, D, C, B = (node(quads[r, q]) for r in 1:4)
        verts[q] = (A, D, C, B)
        cps[q] = (A + D + C + B)/4
        n = LA.cross(C - A, B - D)
        normals[q] = n/LA.norm(n)
        areas[q] = (LA.norm(LA.cross(D - A, C - A))/2 +    # split on diagonal A-C
                    LA.norm(LA.cross(C - A, B - A))/2 +
                    LA.norm(LA.cross(D - A, B - A))/2 +    # split on diagonal D-B
                    LA.norm(LA.cross(C - D, B - D))/2)/2
    end

    # TE quads: strips j=1 (upper; TE edge traversed A->D) and j=npts (lower;
    # TE edge traversed C->B, i.e. naturally opposite). wake_reverse picks the
    # reverse of the traversal as (TE1,TE2).
    te_quads = Int[]
    te_edges = Tuple{SVector{3,Float64},SVector{3,Float64}}[]
    for k in 1:n_span
        for (j, (ia, ib)) in ((1, (1, 2)), (npts, (3, 4)))  # (A,D) or (C,B) slots
            q = (k - 1)*npts + j
            e1, e2 = verts[q][ia], verts[q][ib]
            push!(te_quads, q)
            push!(te_edges, wake_reverse ? (e2, e1) : (e1, e2))
        end
    end

    return (; quads, verts, cps, normals, areas, te_quads, te_edges, Nq, npts)
end

# ----------------- DIRECT INFLUENCE & SOLVER -----------------------------------
"""
    influence_quad!(U, targets, data, gammas, Da; core=CORE)

Accumulate the direct N-body velocity of all quad rings (+ TE wake legs) with
strengths `gammas` at `targets` (vector of SVector). Threaded over targets.
"""
function influence_quad!(U::Matrix{Float64}, targets, data, gammas, Da;
                         core::Float64=CORE)
    Threads.@threads for i in eachindex(targets)
        t = targets[i]
        acc = zero(SVector{3,Float64})
        for q in 1:data.Nq
            v = data.verts[q]
            acc += gammas[q]*vortex_ring_quad(t, v[1], v[2], v[3], v[4], core)
        end
        for (q, (te1, te2)) in zip(data.te_quads, data.te_edges)
            acc += gammas[q]*wake_velocity(t, te1, te2, Da, core)
        end
        U[1, i], U[2, i], U[3, i] = acc[1], acc[2], acc[3]
    end
    return U
end

"""
    solve_quads!(data, Uinf_vec, Da; core=CORE) -> gammas

Assemble the dense influence matrix G[i,j] = u_j(cp_i)·n̂_i (unit-γ quad ring
+ wake legs for TE quads) and solve G γ = -U∞·n̂ with the backslash operator
(LU). Threaded over columns. Returns γ (and G for reuse diagnostics).
"""
function solve_quads!(data, Uinf_vec, Da; core::Float64=CORE, verbose=true)
    Nq = data.Nq
    G = Matrix{Float64}(undef, Nq, Nq)
    is_te = zeros(Int, Nq)                       # 0 or index into te_edges
    for (idx, q) in enumerate(data.te_quads)
        is_te[q] = idx
    end

    t_asm = @elapsed Threads.@threads for j in 1:Nq
        v = data.verts[j]
        tej = is_te[j]
        for i in 1:Nq
            u = vortex_ring_quad(data.cps[i], v[1], v[2], v[3], v[4], core)
            if tej > 0
                te1, te2 = data.te_edges[tej]
                u += wake_velocity(data.cps[i], te1, te2, Da, core)
            end
            G[i, j] = LA.dot(u, data.normals[i])
        end
    end

    rhs = [-LA.dot(Uinf_vec, data.normals[i]) for i in 1:Nq]
    t_lu = @elapsed begin
        Glu = LA.lu!(G)
        gammas = Glu \ rhs
    end
    verbose && @printf("  assembly %.1f s, LU+solve %.1f s (Nq=%d, %d threads)\n",
        t_asm, t_lu, Nq, Threads.nthreads())
    return gammas
end

# ----------------- BENCHMARK GATE ----------------------------------------------
"""
Time the ring kernel and extrapolate the full job (assembly + LU + 1 post
influence pass). Returns estimated seconds.
"""
function estimate_cost(data; nsample=200_000)
    v = data.verts[1]
    targets = data.cps[mod1.(1:nsample, data.Nq)]
    acc = zero(SVector{3,Float64})
    # warmup
    for i in 1:1000
        acc += vortex_ring_quad(targets[i], v[1], v[2], v[3], v[4], CORE)
    end
    t = @elapsed for i in 1:nsample
        acc += vortex_ring_quad(targets[i], v[1], v[2], v[3], v[4], CORE)
    end
    t_eval = t/nsample
    nthreads = Threads.nthreads()
    pair_evals = (Float64(data.Nq)^2 + 2*data.Nq*length(data.te_quads))
    t_asm = pair_evals*t_eval/nthreads
    t_lu = (2/3)*Float64(data.Nq)^3/50e9       # assume ~50 GFLOP/s LU
    t_post = 2*Float64(data.Nq)^2*t_eval/nthreads
    total = t_asm + t_lu + t_post
    @printf("Cost estimate (Nq=%d, %d threads): kernel %.0f ns, assembly %.0f s, LU %.0f s, post %.0f s => total ~%.0f s\n",
        data.Nq, nthreads, t_eval*1e9, t_asm, t_lu, t_post, total)
    return total + (acc[1] == Inf ? 1.0 : 0.0)  # keep acc live
end

# ----------------- CASE DRIVER --------------------------------------------------
function load_vtu_gamma(body_name)
    vtu_path = joinpath(save_path, body_name, "$(body_name).0.vtu")
    isfile(vtu_path) || return nothing
    vtk = pnl.ReadVTK.VTKFile(vtu_path)
    cell_data = pnl.ReadVTK.get_cell_data(vtk)
    return Vector{Float64}(pnl.ReadVTK.get_data(cell_data["gamma"]))
end

function run_quadsolve(n_ch, n_span; budget_s=3600.0, write_vtk=true,
                       enforce_budget=true)
    println("\n#===== QUAD-COLLOCATION SOLVE $(n_ch)x$(n_span) =====#")
    body = sweeploft_wing(n_ch, n_span)        # nodes (+ areas for cross-checks)
    data = quad_solver_data(body.nodes, n_ch, n_span)
    Da = SVector{3,Float64}((Vinf./magVinf)...)
    Uinf_vec = SVector{3,Float64}(Vinf...)

    est = estimate_cost(data)
    if enforce_budget && est > budget_s
        println("Estimated cost $(round(est, digits=0)) s exceeds budget $(budget_s) s; skipping this case.")
        return nothing
    end

    gammas = solve_quads!(data, Uinf_vec, Da)

    # symmetry check: strips k and n_span+1-k mirror each other at equal j
    sym_err = 0.0
    for k in 1:(n_span ÷ 2), j in 1:data.npts
        q1 = (k - 1)*data.npts + j
        q2 = (n_span - k)*data.npts + j
        sym_err = max(sym_err, abs(gammas[q1] - gammas[q2]))
    end
    println("  max |γ(k,j) - γ(mirror k,j)| = ", sym_err)

    # post-solve velocity: U∞ + PV + (-1)*0.5*∇μ
    U0 = Matrix{Float64}(undef, 3, data.Nq)
    influence_quad!(U0, data.cps, data, gammas, Da)
    U0 .+= Vinf

    centroid = Matrix{Float64}(undef, 3, data.Nq)
    normal = Matrix{Float64}(undef, 3, data.Nq)
    for q in 1:data.Nq
        centroid[:, q] .= data.cps[q]
        normal[:, q] .= data.normals[q]
    end
    geom = (; A=data.areas, centroid, normal, mu=gammas, areas=Float64[])
    gradmu = quad_mu_gradient(geom, n_ch, n_span)

    Dhat = Vinf/LA.norm(Vinf)
    Lhat = LA.cross(Dhat, [0.0, 1.0, 0.0])
    res = quad_velocity_CL(geom, U0, gradmu, -1.0, Lhat, Dhat)

    @printf("  CL = %.6f   CD = %.6f\n", res.CL, res.CD)
    @printf("  vs triangle solves at this mesh: diagAC %.4f | diagBD %.4f | experiment %.3f\n",
        0.2678, 0.2459, CLexp)

    # cross-check γ against the saved triangle solves (quad-averaged), if available
    for (label, body_name) in (("diagAC", run_name*"_diagAC_body1"),
                               ("diagBD", run_name*"_diagBD_body1"))
        g_tri = load_vtu_gamma(body_name)
        (g_tri === nothing || length(g_tri) != 2*data.Nq) && continue
        g_tri_quad = (g_tri[1:2:end] .+ g_tri[2:2:end])./2
        rel = LA.norm(gammas .- g_tri_quad)/LA.norm(g_tri_quad)
        println("  rms relative γ difference vs $(label) (quad-averaged): ",
            round(rel, sigdigits=4))
    end

    if write_vtk
        cells = [pnl.WriteVTK.MeshCell(pnl.WriteVTK.VTKCellTypes.VTK_QUAD,
                                       data.quads[:, q]) for q in 1:data.Nq]
        vtk_path = joinpath(save_path, "quadsolve_$(n_ch)x$(n_span)")
        pnl.WriteVTK.vtk_grid(vtk_path, body.nodes, cells) do vtk
            vtk["gamma"] = gammas
            vtk["gauge pressure"] = res.p
            vtk["velocity"] = res.U
            vtk["Cp"] = res.p ./ (0.5*rho*magVinf^2)
        end
        println("  wrote $(vtk_path).vtu")
    end

    return (; gammas, res, data, sym_err)
end

function main_quadsolve()
    case = get(ENV, "FLOWPANEL_QUADSOLVE_CASE", "80:112")
    n_ch, n_span = parse.(Int, split(case, ":"))
    budget_s = parse(Float64, get(ENV, "FLOWPANEL_QUADSOLVE_BUDGET_S", "3600"))

    # smoke case first (small, also serves as JIT warmup)
    smoke = run_quadsolve(24, 24; budget_s, write_vtk=false, enforce_budget=false)
    @assert smoke.sym_err < 1e-8 "Smoke case γ not mirror-symmetric: $(smoke.sym_err)"
    @assert 0.15 < smoke.res.CL < 0.35 "Smoke case CL implausible: $(smoke.res.CL)"

    return run_quadsolve(n_ch, n_span; budget_s)
end

if !isinteractive() && get(ENV, "FLOWPANEL_QUADSOLVE_RUN", "true") == "true"
    quadsolve_result = main_quadsolve()
end
