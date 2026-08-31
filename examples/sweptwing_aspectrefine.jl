#=##############################################################################
# DESCRIPTION
#   Mesh-surgery fix for the triangulated chordwise CL divergence (see
#   "Chordwise divergence root cause" in data/sweptwing_sweeploft/findings.md):
#   recursively bisect high-aspect-ratio triangles (Rivara longest-edge
#   propagation, conforming, O(N)) until every triangle has edge-ratio <= tol,
#   then rerun the chord ladder and check that CL converges and the
#   diagAC/diagBD split collapses.
#
#   The refinement is MESH-AGNOSTIC: new nodes are plain edge midpoints, so
#   the refined mesh represents the IDENTICAL piecewise-linear surface (no
#   surface/CAD knowledge needed, no geometry change at all).
#
#   Run single-threaded:  julia --project examples/sweptwing_aspectrefine.jl
#
#   Environment knobs:
#     FLOWPANEL_ASPECTREF_CASES="40,56,80,112,156"
#     FLOWPANEL_ASPECTREF_NSPAN=48
#     FLOWPANEL_ASPECTREF_TOL=10.0        target max edge-ratio
#     FLOWPANEL_ASPECTREF_MAX_PANELS=80000
#     FLOWPANEL_ASPECTREF_BACKSLASH_MAX=30000   (Krylov+FMM above this)
#     FLOWPANEL_ASPECTREF_KRYLOV_TOL=1e-9
=###############################################################################

ENV["FLOWPANEL_SWEEPLOFT_RUN"] = "false"
ENV["FLOWPANEL_CHORDDIV_RUN"] = "false"
include(joinpath(@__DIR__, "sweptwing_chorddivergence.jl"))

# ----------------- O(N) CONFORMING LONGEST-EDGE BISECTION ---------------------
"""
Growable triangle-mesh state. `neighbor[k, ci]` is the cell across edge k of
cell ci, where edge k connects local nodes (k, mod1(k+1, 3)) — the same
convention as `pnl.calc_neighbors`.
"""
mutable struct RefineState
    nodes::Matrix{Float64}
    nnodes::Int
    cells::Matrix{Int}
    neighbor::Matrix{Int}
    ncells::Int
end

function RefineState(nodes::Matrix{Float64}, cells::Matrix{Int}, cell_cap::Int)
    nn, nc = size(nodes, 2), size(cells, 2)
    node_cap = nn + (cell_cap - nc) ÷ 2 + 8   # one node per 2 added cells
    N = zeros(3, node_cap); N[:, 1:nn] .= nodes
    C = zeros(Int, 3, cell_cap); C[:, 1:nc] .= cells
    NB = zeros(Int, 3, cell_cap)
    NB[:, 1:nc] .= pnl.calc_neighbors(cells)  # guaranteed edge convention
    return RefineState(N, nn, C, NB, nc)
end

function edge_len2(st::RefineState, a::Int, b::Int)
    dx = st.nodes[1, a] - st.nodes[1, b]
    dy = st.nodes[2, a] - st.nodes[2, b]
    dz = st.nodes[3, a] - st.nodes[3, b]
    return dx*dx + dy*dy + dz*dz
end

"""Local index of the longest edge of cell ci, and its squared length."""
function longest_edge(st::RefineState, ci::Int)
    v1, v2, v3 = st.cells[1, ci], st.cells[2, ci], st.cells[3, ci]
    l1, l2, l3 = edge_len2(st, v1, v2), edge_len2(st, v2, v3), edge_len2(st, v3, v1)
    return l1 >= l2 ? (l1 >= l3 ? (1, l1) : (3, l3)) : (l2 >= l3 ? (2, l2) : (3, l3))
end

"""Squared edge-ratio (max/min edge length squared) of cell ci."""
function edge_ratio2(st::RefineState, ci::Int)
    v1, v2, v3 = st.cells[1, ci], st.cells[2, ci], st.cells[3, ci]
    l1, l2, l3 = edge_len2(st, v1, v2), edge_len2(st, v2, v3), edge_len2(st, v3, v1)
    return max(l1, l2, l3)/min(l1, l2, l3)
end

"""Replace the pointer to `old` with `new` in the neighbor slots of cell ci."""
function repoint!(st::RefineState, ci::Int, old::Int, new::Int)
    ci == 0 && return nothing
    for k in 1:3
        if st.neighbor[k, ci] == old
            st.neighbor[k, ci] = new
            return nothing
        end
    end
    error("Cell $ci has no neighbor pointer to $old")
end

"""
Bisect the conforming pair (i, j=neighbor[k,i]) through the midpoint of edge k
of cell i, splitting BOTH cells (always conforming, no hanging nodes).
Children overwrite the parents' slots and two new cells are appended; thanks
to the slot layout only two external neighbor pointers need repointing.
Returns the four child cell indices.
"""
function bisect_pair!(st::RefineState, i::Int, k::Int)
    j = st.neighbor[k, i]
    j > 0 || error("Edge $k of cell $i is a boundary edge; mesh not watertight")
    st.ncells + 2 <= size(st.cells, 2) ||
        error("Cell capacity $(size(st.cells, 2)) exceeded; raise FLOWPANEL_ASPECTREF_MAX_PANELS")

    # Rotate i so the split edge is (v1, v2); v3 opposite
    v1 = st.cells[k, i]
    v2 = st.cells[mod1(k + 1, 3), i]
    v3 = st.cells[mod1(k + 2, 3), i]
    B_i = st.neighbor[mod1(k + 1, 3), i]   # across (v2, v3)
    C_i = st.neighbor[mod1(k + 2, 3), i]   # across (v3, v1)

    # Locate the shared edge in j: consistent winding => directed edge (v2, v1)
    kj = 0
    for kk in 1:3
        if st.cells[kk, j] == v2 && st.cells[mod1(kk + 1, 3), j] == v1
            kj = kk
            break
        end
    end
    kj > 0 || error("Cells $i and $j do not share directed edge ($v2,$v1); winding inconsistent")
    u3 = st.cells[mod1(kj + 2, 3), j]      # j's opposite vertex
    B_j = st.neighbor[mod1(kj + 1, 3), j]  # across (v1, u3)
    C_j = st.neighbor[mod1(kj + 2, 3), j]  # across (u3, v2)

    # Midpoint node (ON the existing edge: geometry unchanged)
    st.nnodes += 1
    m = st.nnodes
    m <= size(st.nodes, 2) || error("Node capacity exceeded")
    for d in 1:3
        st.nodes[d, m] = 0.5*(st.nodes[d, v1] + st.nodes[d, v2])
    end

    i2 = st.ncells + 1
    j2 = st.ncells + 2
    st.ncells += 2

    # child1 in slot i: (v1, m, v3) — keeps i's external edge (v3, v1) -> C_i
    st.cells[1, i] = v1; st.cells[2, i] = m; st.cells[3, i] = v3
    st.neighbor[1, i] = j2; st.neighbor[2, i] = i2; st.neighbor[3, i] = C_i
    # child2 appended: (m, v2, v3) — takes i's external edge (v2, v3) -> B_i
    st.cells[1, i2] = m; st.cells[2, i2] = v2; st.cells[3, i2] = v3
    st.neighbor[1, i2] = j; st.neighbor[2, i2] = B_i; st.neighbor[3, i2] = i
    # child3 in slot j: (v2, m, u3) — keeps j's external edge (u3, v2) -> C_j
    st.cells[1, j] = v2; st.cells[2, j] = m; st.cells[3, j] = u3
    st.neighbor[1, j] = i2; st.neighbor[2, j] = j2; st.neighbor[3, j] = C_j
    # child4 appended: (m, v1, u3) — takes j's external edge (v1, u3) -> B_j
    st.cells[1, j2] = m; st.cells[2, j2] = v1; st.cells[3, j2] = u3
    st.neighbor[1, j2] = i; st.neighbor[2, j2] = B_j; st.neighbor[3, j2] = j

    repoint!(st, B_i, i, i2)
    repoint!(st, B_j, j, j2)

    return (i, i2, j, j2)
end

"""
Single-sided bisection of cell i through its boundary edge k (no cell on the
other side, e.g. the open wingtip rings). Always conforming.
"""
function bisect_boundary!(st::RefineState, i::Int, k::Int)
    st.neighbor[k, i] == 0 || error("Edge $k of cell $i is not a boundary edge")
    st.ncells + 1 <= size(st.cells, 2) ||
        error("Cell capacity $(size(st.cells, 2)) exceeded; raise FLOWPANEL_ASPECTREF_MAX_PANELS")

    v1 = st.cells[k, i]
    v2 = st.cells[mod1(k + 1, 3), i]
    v3 = st.cells[mod1(k + 2, 3), i]
    B_i = st.neighbor[mod1(k + 1, 3), i]
    C_i = st.neighbor[mod1(k + 2, 3), i]

    st.nnodes += 1
    m = st.nnodes
    m <= size(st.nodes, 2) || error("Node capacity exceeded")
    for d in 1:3
        st.nodes[d, m] = 0.5*(st.nodes[d, v1] + st.nodes[d, v2])
    end

    i2 = st.ncells + 1
    st.ncells += 1
    st.cells[1, i] = v1; st.cells[2, i] = m; st.cells[3, i] = v3
    st.neighbor[1, i] = 0; st.neighbor[2, i] = i2; st.neighbor[3, i] = C_i
    st.cells[1, i2] = m; st.cells[2, i2] = v2; st.cells[3, i2] = v3
    st.neighbor[1, i2] = 0; st.neighbor[2, i2] = B_i; st.neighbor[3, i2] = i
    repoint!(st, B_i, i, i2)
    return (i, i2)
end

const BISECTION_COUNT = Ref(0)

"""
Refine cell `s` until every affected cell satisfies edge_ratio <= tol, by
FORCED conforming bisection: split the over-tolerance cell's own longest
edge directly, splitting the (single, found by direct `neighbor[k,i]`
lookup) cell across that edge through the same midpoint — always conforming,
no hanging nodes — and re-enqueue any of the children still over tolerance.

Note: classic Rivara LEPP (bisect only mutual-longest "terminal pairs") was
tried first and cascades globally on this uniformly anisotropic grid — every
bulk quad's diagonal is its longest edge, so LEPP chains bisect compliant
bulk cells wholesale (observed 7.7k -> 57k+ cells without terminating).
Forced bisection stays local: a collateral split at most ~doubles the
neighbor's aspect ratio, so propagation decays geometrically away from the
sliver bands. Termination is guarded by the capacity cap (+ state dump).
"""
function force_refine!(st::RefineState, s::Int, tol2::Float64, stack::Vector{Int};
                       frozen=nothing)
    empty!(stack)
    push!(stack, s)
    while !isempty(stack)
        i = pop!(stack)
        edge_ratio2(st, i) <= tol2 && continue
        k, _ = longest_edge(st, i)
        if frozen !== nothing &&
           frozen(st.nodes, st.cells[k, i], st.cells[mod1(k + 1, 3), i])
            continue   # leave this cell over-tol rather than split a frozen edge
        end
        j = st.neighbor[k, i]
        kids = j == 0 ? bisect_boundary!(st, i, k) : bisect_pair!(st, i, k)
        for ci in kids
            edge_ratio2(st, ci) > tol2 && push!(stack, ci)
        end
        BISECTION_COUNT[] += 1
        if BISECTION_COUNT[] % 25_000 == 0
            println("    [refine] $(BISECTION_COUNT[]) bisections, ",
                    "ncells = $(st.ncells)")
            flush(stdout)
        end
    end
    return nothing
end

"""
    refine_aspect(nodes, cells; tol, max_panels) -> (nodes, cells, nsplits)

Conforming longest-edge bisection until every triangle has edge-ratio <= tol.
Mesh-agnostic (midpoint nodes only). Requires a watertight, consistently
wound triangle mesh.
"""
function refine_aspect(nodes::Matrix{Float64}, cells::Matrix{Int};
                       tol::Real, max_panels::Int, frozen=nothing)
    st = RefineState(nodes, cells, max_panels)
    tol2 = Float64(tol)^2
    nc0 = st.ncells
    BISECTION_COUNT[] = 0
    stack = Int[]
    try
        # Children are handled recursively via the stack; slots reused by
        # children are re-tested on pop, so stale entries are harmless.
        for ci in 1:nc0
            force_refine!(st, ci, tol2, stack; frozen)
        end
    catch err
        _dump_refine_state(st, tol)
        rethrow(err)
    end
    return (st.nodes[:, 1:st.nnodes], st.cells[:, 1:st.ncells],
            st.neighbor[:, 1:st.ncells])
end

"""Diagnostic dump on refinement failure: AR histogram and x-position of cells."""
function _dump_refine_state(st::RefineState, tol)
    nbad = 0
    ar_buckets = zeros(Int, 6)   # <tol, tol-2tol, 2-4, 4-8, 8-16, >16tol
    xbad = Float64[]
    for ci in 1:st.ncells
        ar = sqrt(edge_ratio2(st, ci))
        b = ar <= tol ? 1 : min(6, 2 + floor(Int, log2(ar/tol)))
        ar_buckets[b] += 1
        if ar > tol
            nbad += 1
            length(xbad) < 200_000 &&
                push!(xbad, sum(st.nodes[1, st.cells[k, ci]] for k in 1:3)/3)
        end
    end
    println("    [refine DUMP] ncells=$(st.ncells), over-tol=$(nbad), ",
            "AR buckets (<tol,1-2,2-4,4-8,8-16,>16 xtol): ", ar_buckets)
    if !isempty(xbad)
        # x - |y|tanλ is the local chord coordinate proxy; just use raw x here
        println("    [refine DUMP] over-tol cell x-centroid range: ",
                extrema(xbad), ", mean: ", sum(xbad)/length(xbad))
    end
    flush(stdout)
end

# ----------------- BODY CONSTRUCTION -------------------------------------------
"""
Build the refined RigidWakeBody. Returns `(body, info)` where info has the
panel counts and post-refinement max edge-ratio.
"""
function build_refined_body(n_ch::Int, n_span::Int; swap_diagonals::Bool,
                            tol::Real, max_panels::Int, freeze_te::Bool=false,
                            freeze_band::Real=0.0, allow_annulus=nothing)
    base = sweeploft_wing(n_ch, n_span; mirror_diagonals=true, swap_diagonals)
    nodes0, cells0 = copy(base.nodes), copy(base.cells)
    panels_base = size(cells0, 2)

    # freeze_te: never split a TE-chain edge (diagnostic mode — leaves the
    # TE-touching slivers unrefined so the wake topology is unchanged)
    # freeze_band: never split an edge with BOTH endpoints within chord
    # fraction `freeze_band` of the TE or LE line (diagnostic mode — refines
    # only away from the sharp edges; band cells stay over-tol)
    frozen = nothing
    if allow_annulus !== nothing
        # Inverse mode: ONLY edges fully inside the chordwise annulus
        # xt ∈ [lo, hi] (mirrored at the TE side) may be split. Confines the
        # conforming cascade to the annulus — minimal, edge-avoiding surgery.
        lo, hi = allow_annulus
        in_ann = function (nodes, ni)
            xt = clamp((nodes[1, ni] - abs(nodes[2, ni])*tan(lambda_rad))/c, 0.0, 1.0)
            return (lo <= xt <= hi) || (1 - hi <= xt <= 1 - lo)
        end
        frozen = (nodes, a, b) -> !(in_ann(nodes, a) && in_ann(nodes, b))
    elseif freeze_band > 0
        delta = Float64(freeze_band)
        in_band = function (nodes, ni)
            xt = clamp((nodes[1, ni] - abs(nodes[2, ni])*tan(lambda_rad))/c, 0.0, 1.0)
            return xt < delta || xt > 1 - delta
        end
        frozen = (nodes, a, b) -> in_band(nodes, a) && in_band(nodes, b)
    elseif freeze_te
        tolg0 = 1e-9*c
        iste = [abs(nodes0[1, ni] - (abs(nodes0[2, ni])*tan(lambda_rad) + c)) < tolg0 &&
                abs(nodes0[3, ni]) < tolg0 for ni in axes(nodes0, 2)]
        # nodes added during refinement (index > base count) are never TE nodes
        frozen = (nodes, a, b) -> a <= length(iste) && b <= length(iste) &&
                                  iste[a] && iste[b]
    end

    nodes, cells, neighbor = refine_aspect(nodes0, cells0; tol, max_panels, frozen)
    panels = size(cells, 2)

    # Audit: adjacency consistent with a fresh rebuild, watertight, outward
    neighbor == pnl.calc_neighbors(cells) ||
        error("Incrementally maintained adjacency differs from rebuild")
    # The base wing has open tips (no caps), so watertightness must simply be
    # UNCHANGED by refinement (boundary edges bisect into boundary edges).
    watertight, _ = pnl.iswatertight(nodes, cells)
    watertight_base, _ = pnl.iswatertight(nodes0, cells0)
    watertight == watertight_base ||
        error("Refinement changed watertightness: $watertight_base -> $watertight")
    vol = 0.0
    for ci in axes(cells, 2)
        a = view(nodes, :, cells[1, ci])
        bb = view(nodes, :, cells[2, ci])
        cc = view(nodes, :, cells[3, ci])
        vol += LA.dot(a, LA.cross(bb, cc))/6
    end
    vol > 0 || error("Refined mesh signed volume $(vol) <= 0")
    ar_max = 0.0
    for ci in axes(cells, 2)
        v1, v2, v3 = cells[1, ci], cells[2, ci], cells[3, ci]
        p1, p2, p3 = view(nodes, :, v1), view(nodes, :, v2), view(nodes, :, v3)
        l1, l2, l3 = LA.norm(p2 .- p1), LA.norm(p3 .- p2), LA.norm(p1 .- p3)
        ar_max = max(ar_max, max(l1, l2, l3)/min(l1, l2, l3))
    end
    frozen === nothing && ar_max > tol + 1e-9 &&
        error("Max edge-ratio $(ar_max) > tol $(tol) after refinement")
    n_over = count(ci -> begin
            v1, v2, v3 = cells[1, ci], cells[2, ci], cells[3, ci]
            p1, p2, p3 = view(nodes, :, v1), view(nodes, :, v2), view(nodes, :, v3)
            l1, l2, l3 = LA.norm(p2 .- p1), LA.norm(p3 .- p2), LA.norm(p1 .- p3)
            max(l1, l2, l3)/min(l1, l2, l3) > tol + 1e-9
        end, axes(cells, 2))

    # TE line of this loft: x = |y| tanλ + c, z = 0 (mesh-specific; a general
    # tool would carry the TE chain through the splits instead)
    tolg = 1e-9*c
    te = [ni for ni in axes(nodes, 2)
          if abs(nodes[1, ni] - (abs(nodes[2, ni])*tan(lambda_rad) + c)) < tolg &&
             abs(nodes[3, ni]) < tolg]
    sort!(te; by=ni -> nodes[2, ni])
    shedding = pnl.calc_shedding(nodes, cells, te, zeros(Float64, 3, 0))
    size(shedding, 2) == length(te) - 1 ||
        error("Shedding has $(size(shedding, 2)) edges, expected $(length(te) - 1)")

    body = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}(
        nodes, cells, [shedding];
        watertight, ensure_winding=false, core_size=1e-10)
    wake_direction = reshape(Vinf ./ magVinf, :, 1)
    for i in eachindex(body.Das)
        body.Das[i] .= repeat(wake_direction, 1, size(body.Das[i], 2))
    end
    info = (; panels_base, panels, growth=panels/panels_base - 1, ar_max,
            n_over, n_te_nodes=length(te))
    return body, info
end

# ----------------- ONE CASE -----------------------------------------------------
function run_refined_case(n_ch::Int, n_span::Int; swap_diagonals::Bool,
                          tol::Real, max_panels::Int, backslash_max::Int,
                          krylov_tol::Real, variant_name::String="",
                          freeze_band::Real=0.0, allow_annulus=nothing)
    t0 = time()
    body, info = build_refined_body(n_ch, n_span; swap_diagonals, tol,
                                    max_panels, freeze_band, allow_annulus)
    N = body.ncells

    use_dense = N <= backslash_max
    solver = use_dense ? pnl.Backslash(body) :
        pnl.KrylovSolver(body; backend=pnl.FastMultipoleBackend(),
                         method=:gmres, atol=krylov_tol, rtol=krylov_tol,
                         itmax=200)

    Dhat = Vinf/LA.norm(Vinf)
    Lhat = LA.cross(Dhat, [0.0, 1.0, 0.0])
    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false)
    force = pnl.ForceMonitor(1, 1; i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c_ref),
        correct_kuttacondition=false, verbose=false)

    pnl.steady!(body, pnl.ReferenceFrame(body), Vinf;
        body_solvers=solver, backend=pnl.FastMultipoleBackend(),
        monitors=(pressure, force), path=nothing, verbose=false)
    CL = LA.dot(force.force[:, 1], Lhat)
    CD = LA.dot(force.force[:, 1], Dhat)
    solver = nothing; GC.gc()

    # Flow-tangency audit (mandatory sanity check on the Krylov levels)
    tang = [LA.dot(view(body.velocity, :, i), view(body.normals, :, i))
            for i in 1:N] ./ magVinf
    te_node = falses(size(body.nodes, 2))
    tolg = 1e-9*c
    for ni in axes(body.nodes, 2)
        te_node[ni] = abs(body.nodes[1, ni] -
                          (abs(body.nodes[2, ni])*tan(lambda_rad) + c)) < tolg &&
                      abs(body.nodes[3, ni]) < tolg
    end
    te_cells = [ci for ci in 1:N if any(te_node[body.cells[k, ci]] for k in 1:3)]
    tang_mean, tang_max = abs_stats(tang, 1:N)
    _, tang_te_max = abs_stats(tang, te_cells)

    # Kutta-Joukowski CL from gamma alone (solve-vs-reconstruction
    # discriminator; sign flip per the body's gamma convention, as in
    # sweptwing_chorddivergence.jl)
    gam = view(body.strength, :, pnl.get_Gammai(body))
    sh = body.shedding[1]
    CL_KJ = 0.0
    for e in axes(sh, 2)
        na = body.cells[sh[2, e], sh[1, e]]
        nb = body.cells[sh[3, e], sh[1, e]]
        dy = abs(body.nodes[2, nb] - body.nodes[2, na])
        gpj = sh[4, e] > 0 ? gam[sh[4, e]] : 0.0
        CL_KJ += (gam[sh[1, e]] - gpj)*dy
    end
    CL_KJ *= -rho*magVinf/(0.5*rho*magVinf^2*Sref)

    # TE-row jump stats: bounded-vs-1/h_min discriminator
    Jtri = pnl.compute_mu_gradient!(zeros(3, N), body.controlpoints,
        body.normals, body.cells, body.neighbor,
        view(body.strength, :, pnl.get_Gammai(body)),
        view(body.shedding_full, 1:2, :); scale=0.5, nodes=body.nodes)
    Jte_mean, Jte_max = colnorm_stats(Jtri, te_cells)

    return (; variant=variant_name, n_ch, n_span, tol, freeze_band,
        annulus=allow_annulus === nothing ? "" :
                "$(allow_annulus[1]):$(allow_annulus[2])",
        n_over=info.n_over,
        panels_base=info.panels_base, panels=info.panels,
        growth_pct=100*info.growth, ar_max=info.ar_max,
        solver=use_dense ? "backslash" : "krylov", CL, CD, CL_KJ,
        gam_max=maximum(abs.(gam)),
        tang_mean, tang_max, tang_te_max,
        Jte_mean, Jte_max, h_min=c*(1 - cos(pi/n_ch))/2,
        elapsed=time() - t0)
end

# ----------------- DRIVER -------------------------------------------------------
function main_aspectrefine()
    cases = [parse(Int, s) for s in
        split(get(ENV, "FLOWPANEL_ASPECTREF_CASES", "40,56,80,112,156"), ",")]
    n_span = parse(Int, get(ENV, "FLOWPANEL_ASPECTREF_NSPAN", "48"))
    tol = parse(Float64, get(ENV, "FLOWPANEL_ASPECTREF_TOL", "10.0"))
    max_panels = parse(Int, get(ENV, "FLOWPANEL_ASPECTREF_MAX_PANELS", "80000"))
    backslash_max = parse(Int, get(ENV, "FLOWPANEL_ASPECTREF_BACKSLASH_MAX", "30000"))
    krylov_tol = parse(Float64, get(ENV, "FLOWPANEL_ASPECTREF_KRYLOV_TOL", "1e-9"))
    freeze_band = parse(Float64, get(ENV, "FLOWPANEL_ASPECTREF_FREEZE_BAND", "0"))
    ann_s = get(ENV, "FLOWPANEL_ASPECTREF_ANNULUS", "")
    allow_annulus = isempty(ann_s) ? nothing :
        Tuple(parse.(Float64, split(ann_s, ":")))
    out_csv = joinpath("data", "sweptwing_sweeploft",
        (freeze_band > 0 || allow_annulus !== nothing) ?
            "aspectrefine_band.csv" : "aspectrefine_convergence.csv")

    println("Aspect-refined chord ladder: n_ch = ", cases, " x n_span = ", n_span,
            ", tol = ", tol, ", freeze_band = ", freeze_band,
            ", annulus = ", allow_annulus)

    rows = []
    for n_ch in cases, (vname, swap) in (("diagAC", false), ("diagBD", true))
        r = try
            run_refined_case(n_ch, n_span; swap_diagonals=swap, tol, max_panels,
                             backslash_max, krylov_tol, variant_name=vname,
                             freeze_band, allow_annulus)
        catch err
            println("SKIPPED $vname $n_ch: ", sprint(showerror, err))
            flush(stdout)
            continue
        end
        @printf("%7s %4d %6d->%6d (+%4.0f%%) AR %5.1f %9s CL %+.6f CLkj %+.6f maxg %7.3f CD %+.6f tangmax %.1e (TE %.1e) Jte_max %9.3e t %.0fs\n",
            r.variant, r.n_ch, r.panels_base, r.panels, r.growth_pct, r.ar_max,
            r.solver, r.CL, r.CL_KJ, r.gam_max, r.CD, r.tang_max,
            r.tang_te_max, r.Jte_max, r.elapsed)
        flush(stdout)
        push!(rows, r)
        CSV.write(out_csv, DataFrame(rows))
        GC.gc()
    end

    println("\n#===== REFINED vs UNREFINED CL (unrefined CL_solver from chord_divergence.csv) =====#")
    unref = Dict{Tuple{String, Int}, Float64}()
    cd_csv = joinpath("data", "sweptwing_sweeploft", "chord_divergence.csv")
    if isfile(cd_csv)
        for r in CSV.File(cd_csv)
            unref[(r.variant, r.n_ch)] = r.CL_solver
        end
    end
    @printf("%7s %5s %8s %12s %12s %10s\n", "variant", "n_ch", "panels",
        "CL_refined", "CL_unref", "Jte_max")
    for r in rows
        cu = get(unref, (r.variant, r.n_ch), NaN)
        @printf("%7s %5d %8d %+12.6f %+12.6f %10.3e\n",
            r.variant, r.n_ch, r.panels, r.CL, cu, r.Jte_max)
    end
    for v in ("diagAC", "diagBD")
        cls = [r.CL for r in rows if r.variant == v]
        length(cls) > 1 &&
            println("  $v ΔCL per step: ", join(map(x -> @sprintf("%+.6f", x),
                                                    diff(cls)), " "))
    end
    println("CSV written to ", out_csv)
    return rows
end

"""Mesh-only test: refinement + audits, no solve (fast iteration on the mesher)."""
function meshtest()
    cases = [parse(Int, s) for s in
        split(get(ENV, "FLOWPANEL_ASPECTREF_CASES", "40,56,80,112,156"), ",")]
    n_span = parse(Int, get(ENV, "FLOWPANEL_ASPECTREF_NSPAN", "48"))
    tol = parse(Float64, get(ENV, "FLOWPANEL_ASPECTREF_TOL", "10.0"))
    max_panels = parse(Int, get(ENV, "FLOWPANEL_ASPECTREF_MAX_PANELS", "80000"))
    freeze_band = parse(Float64, get(ENV, "FLOWPANEL_ASPECTREF_FREEZE_BAND", "0"))
    ann_s = get(ENV, "FLOWPANEL_ASPECTREF_ANNULUS", "")
    allow_annulus = isempty(ann_s) ? nothing :
        Tuple(parse.(Float64, split(ann_s, ":")))
    for n_ch in cases, (vname, swap) in (("diagAC", false), ("diagBD", true))
        t = @elapsed body, info = build_refined_body(n_ch, n_span;
            swap_diagonals=swap, tol, max_panels, freeze_band, allow_annulus)
        @printf("%7s %4d: %6d -> %6d (+%5.1f%%), max AR %6.2f, over-tol %5d, TE nodes %4d, %.1fs\n",
            vname, n_ch, info.panels_base, info.panels, 100*info.growth,
            info.ar_max, info.n_over, info.n_te_nodes, t)
        flush(stdout)
    end
end

if !isinteractive() && get(ENV, "FLOWPANEL_ASPECTREF_RUN", "true") == "true"
    if get(ENV, "FLOWPANEL_ASPECTREF_MESHTEST", "false") == "true"
        meshtest()
    else
        aspectref_rows = main_aspectrefine()
    end
end
