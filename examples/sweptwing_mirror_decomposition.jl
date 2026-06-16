#=##############################################################################
# DESCRIPTION
    Diagnostic for the sweptwing mirrored-discretization discrepancy: loads the
    saved 40x48 VTK state for the +y-half-mirrored and -y-half-mirrored meshes
    (no re-solve), decomposes the surface velocity into the panel-induced (PV)
    part and the ∇μ half-jump part, and compares them across the two meshes at
    the level of the shared structured quads. Also compares the quad pairs
    chosen by compute_mu_gradient!'s quad-consistent pass and runs a controlled
    same-μ gradient experiment to isolate triangulation sensitivity.

    Run: julia --project examples/sweptwing_mirror_decomposition.jl
=###############################################################################

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
import LinearAlgebra
const _norm = LinearAlgebra.norm

# ----------------- PARAMETERS (must match examples/sweptwing.jl) --------------
AOA             = 4.2
aoa_tag         = replace(string(AOA), "." => "p")
magVinf         = 30.0
Vinf            = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)]

b               = 98*0.0254
ar              = 5.0
tr              = 1.0
twist_root      = 0
twist_tip       = 0
lambda          = 45
gamma           = 0
airfoil         = "airfoil-rae101.csv"
airfoil_path    = joinpath(pnl.examples_path, "data")

n_rfl           = 40
NDIVS_rfl = [ (0.25, n_rfl,   10.0, false),
              (0.50, n_rfl,    1.0, true),
              (0.25, n_rfl, 1/10.0, false)]
n_span_full     = 48
NDIVS_span      = [(1.0, n_span_full, 1.0, true)]

grid_tag        = "nrf$(n_rfl)_nspan$(n_span_full)"
run_name        = "sweptwing000_" * grid_tag
load_path       = joinpath("data", run_name)
positive_body_name = run_name * "_bernoulli_AOA$(aoa_tag)_body1"
negative_body_name = run_name * "_negative_half_mirror_bernoulli_AOA$(aoa_tag)"
out_path        = joinpath("data", "sweptwing_mirror_decomposition")

bodytype = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}
bodyoptargs = (;)

# Copied from examples/sweptwing.jl
function simplewing_mirrored_from_negative(b, ar, tr, twist_root, twist_tip, lambda, gamma;
                                           bodytype, bodyoptargs=(;),
                                           airfoil_root, airfoil_tip, airfoil_path,
                                           rfl_NDIVS, span_NDIVS, delim=",",
                                           mirror_tol=100eps(Float64),
                                           reference_nodes=nothing)
    half = simplewing(b, ar, tr, twist_root, twist_tip, lambda, gamma;
                      bodytype=bodytype, bodyoptargs=bodyoptargs,
                      airfoil_root=airfoil_root, airfoil_tip=airfoil_tip,
                      airfoil_path=airfoil_path,
                      rfl_NDIVS=rfl_NDIVS,
                      delim=delim,
                      span_NDIVS=span_NDIVS,
                      b_low=-1.0, b_up=0.0)

    half_nodes = half.nodes
    half_cells = half.cells
    mirror_index = Vector{Int}(undef, size(half_nodes, 2))
    nodes = copy(half_nodes)
    for ni in axes(half_nodes, 2)
        if abs(half_nodes[2, ni]) <= mirror_tol
            mirror_index[ni] = ni
        else
            nodes = hcat(nodes, [half_nodes[1, ni], -half_nodes[2, ni], half_nodes[3, ni]])
            mirror_index[ni] = size(nodes, 2)
        end
    end

    half_centers_y = [sum(half_nodes[2, half_cells[:, ci]]) / 3 for ci in axes(half_cells, 2)]
    neg_order = sort(collect(axes(half_cells, 2)); by=ci -> half_centers_y[ci])
    pos_order = sort(collect(axes(half_cells, 2)); by=ci -> -half_centers_y[ci])

    cells = Matrix{Int}(undef, 3, 2 * size(half_cells, 2))
    out_ci = 0
    for ci in neg_order
        out_ci += 1
        cells[:, out_ci] .= half_cells[:, ci]
    end
    for ci in pos_order
        out_ci += 1
        cells[:, out_ci] .= reverse(mirror_index[half_cells[:, ci]])
    end

    te_nodes = Int[]
    for col in eachcol(half.shedding[1])
        pi, nia, nib = col[1], col[2], col[3]
        push!(te_nodes, half_cells[nia, pi])
        push!(te_nodes, half_cells[nib, pi])
    end
    full_te_nodes = unique(vcat(te_nodes, mirror_index[te_nodes]))

    if !isnothing(reference_nodes)
        size(reference_nodes) == size(nodes) ||
            error("Cannot reindex negative-half mirror: reference node size $(size(reference_nodes)) differs from generated node size $(size(nodes)).")

        old_to_reference = zeros(Int, size(nodes, 2))
        reference_used = falses(size(reference_nodes, 2))
        for old_i in axes(nodes, 2)
            match_i = 0
            for ref_i in axes(reference_nodes, 2)
                reference_used[ref_i] && continue
                if maximum(abs.(view(nodes, :, old_i) .- view(reference_nodes, :, ref_i))) <= mirror_tol
                    match_i = ref_i
                    break
                end
            end
            match_i != 0 ||
                error("Cannot reindex negative-half mirror: generated node $old_i has no matching reference node within $(mirror_tol).")
            old_to_reference[old_i] = match_i
            reference_used[match_i] = true
        end

        all(reference_used) ||
            error("Cannot reindex negative-half mirror: not all reference nodes were matched.")

        nodes = copy(reference_nodes)
        cells = old_to_reference[cells]
        full_te_nodes = unique(old_to_reference[full_te_nodes])
    end

    sort!(full_te_nodes; by=ni -> nodes[2, ni])
    shedding = pnl.calc_shedding(nodes, cells, full_te_nodes, zeros(eltype(nodes), 3, 0))

    watertight, _ = pnl.iswatertight(nodes, cells)
    final_bodyoptargs = merge((ensure_winding=false,), bodyoptargs)
    return bodytype(nodes, cells, [shedding]; watertight, final_bodyoptargs...)
end

# Copied from examples/sweptwing.jl
function cell_scalar_to_nodes(body, cell_values)
    length(cell_values) == body.ncells ||
        error("Expected one cell value per panel; got $(length(cell_values)) values for $(body.ncells) panels.")
    areas = pnl.calc_areas(body.nodes, body.cells)
    values = zeros(promote_type(eltype(cell_values), eltype(areas)), size(body.nodes, 2))
    weights = zeros(eltype(areas), size(body.nodes, 2))
    for ci in 1:body.ncells
        area = areas[ci]
        value = cell_values[ci]
        for ni in view(body.cells, :, ci)
            values[ni] += area * value
            weights[ni] += area
        end
    end
    unused = findall(iszero, weights)
    isempty(unused) || error("Cannot interpolate cell values to nodes; $(length(unused)) nodes have no incident panels.")
    values ./= weights
    return values
end

# ----------------- BUILD GEOMETRY AND LOAD SAVED STATE ------------------------
println("Generating +y-mirrored body...")
@time body_pos = simplewing_mirrored(b, ar, tr, twist_root, twist_tip, lambda, gamma;
                                     bodytype=bodytype, bodyoptargs=bodyoptargs,
                                     airfoil_root=airfoil, airfoil_tip=airfoil,
                                     airfoil_path=airfoil_path,
                                     rfl_NDIVS=NDIVS_rfl,
                                     delim=",",
                                     span_NDIVS=NDIVS_span)

println("Generating -y-mirrored body (reindexed to shared nodes)...")
@time body_neg = simplewing_mirrored_from_negative(
    b, ar, tr, twist_root, twist_tip, lambda, gamma;
    bodytype=bodytype, bodyoptargs=bodyoptargs,
    airfoil_root=airfoil, airfoil_tip=airfoil,
    airfoil_path=airfoil_path,
    rfl_NDIVS=NDIVS_rfl,
    delim=",",
    span_NDIVS=NDIVS_span,
    reference_nodes=body_pos.nodes)

# Semi-infinite wake direction along the freestream, as in examples/sweptwing.jl
wake_direction = reshape(Vinf ./ magVinf, :, 1)
for bd in (body_pos, body_neg), i in eachindex(bd.Das)
    bd.Das[i] .= repeat(wake_direction, 1, size(bd.Das[i], 2))
end

node_mismatch = maximum(abs.(body_pos.nodes .- body_neg.nodes))
println("Max node coordinate mismatch: $(node_mismatch)")
node_mismatch <= 1e-10 || error("Node arrays differ between meshes; quad-level comparison invalid.")

println("Loading saved VTK state (no re-solve)...")
loaded_pos = pnl._load_body_step!(body_pos, load_path, positive_body_name, 0)
loaded_neg = pnl._load_body_step!(body_neg, load_path, negative_body_name, 0)

for (label, body) in (("+y", body_pos), ("-y", body_neg))
    speeds = [_norm(c) for c in eachcol(body.velocity)]
    gammas = view(body.strength, :, pnl.get_Gammai(body))
    println("  $label: |velocity| in [$(round(minimum(speeds), sigdigits=4)), " *
            "$(round(maximum(speeds), sigdigits=4))] m/s, max|gamma|=$(round(maximum(abs, gammas), sigdigits=4))")
    all(isfinite, body.velocity) || error("Loaded velocity for $label mesh has non-finite entries.")
    maximum(abs, gammas) > 0 || error("Loaded gamma for $label mesh is all zero.")
end

# ----------------- QUESTION 2: QUAD-PAIRING COMPARISON ------------------------
te_info(body) = view(body.shedding_full, 1:2, :)

function quad_pairs(body; normal_dot_min=cos(pi/4))
    pnl.calc_normals!(body)
    cells, neighbors = body.cells, body.neighbor
    nodes, normals = body.nodes, body.normals
    te = te_info(body)
    M = size(cells, 2)
    paired = falses(M)
    pairs = Tuple{Int,Int}[]
    for i in 1:M
        paired[i] && continue
        j = pnl._quad_candidate_neighbor(nodes, cells, neighbors, normals, te, i;
                                         normal_dot_min=normal_dot_min)
        j > 0 || continue
        paired[j] && continue
        pnl._quad_candidate_neighbor(nodes, cells, neighbors, normals, te, j;
                                     normal_dot_min=normal_dot_min) == i || continue
        push!(pairs, (i, j))
        paired[i] = true
        paired[j] = true
    end
    return pairs, paired
end

quad_key(cells, i, j) = Tuple(sort(unique(vcat(cells[:, i], cells[:, j]))))

pairs_pos, paired_pos = quad_pairs(body_pos)
pairs_neg, paired_neg = quad_pairs(body_neg)

keys_pos = Set(quad_key(body_pos.cells, i, j) for (i, j) in pairs_pos)
keys_neg = Set(quad_key(body_neg.cells, i, j) for (i, j) in pairs_neg)

function describe_quads(label, keys_only, nodes, b)
    isempty(keys_only) && return
    println("  Quads only in $label pairing: $(length(keys_only))")
    ys = Float64[]; xs = Float64[]
    for key in keys_only
        push!(ys, sum(nodes[2, k] for k in key) / 4 * 2 / b)
        push!(xs, sum(nodes[1, k] for k in key) / 4)
    end
    println("    2y/b range: [$(round(minimum(ys), sigdigits=3)), $(round(maximum(ys), sigdigits=3))], " *
            "x range: [$(round(minimum(xs), sigdigits=3)), $(round(maximum(xs), sigdigits=3))]")
end

println("\n#===== QUAD-PAIRING COMPARISON (compute_mu_gradient! quad-consistent pass) =====#")
println("  +y mesh: $(length(pairs_pos)) pairs, $(count(.!paired_pos)) unpaired panels of $(body_pos.ncells)")
println("  -y mesh: $(length(pairs_neg)) pairs, $(count(.!paired_neg)) unpaired panels of $(body_neg.ncells)")
common_keys = intersect(keys_pos, keys_neg)
println("  Identical quads (same 4 nodes) in both pairings: $(length(common_keys))")
describe_quads("+y", setdiff(keys_pos, keys_neg), body_pos.nodes, b)
describe_quads("-y", setdiff(keys_neg, keys_pos), body_neg.nodes, b)

# ----------------- CROSS-MESH STRUCTURED-QUAD CORRESPONDENCE ------------------
# Independent of the gradient pairing: A-mesh triangles i,j sharing an edge form
# a structured quad iff the B mesh contains the two complement triangles of the
# same 4-node set (B splits the quad by the other diagonal).
function structured_quads(bodyA, bodyB)
    triB = Dict{NTuple{3,Int},Int}()
    for ci in axes(bodyB.cells, 2)
        triB[Tuple(sort(bodyB.cells[:, ci]))] = ci
    end
    quads = Tuple{NTuple{4,Int},Tuple{Int,Int},Tuple{Int,Int}}[]
    usedA = falses(bodyA.ncells)
    for i in axes(bodyA.cells, 2)
        usedA[i] && continue
        for k in 1:3
            j = bodyA.neighbor[k, i]
            (j > 0 && j != i && !usedA[j]) || continue
            key = quad_key(bodyA.cells, i, j)
            length(key) == 4 || continue
            triplesA = (Tuple(sort(bodyA.cells[:, i])), Tuple(sort(bodyA.cells[:, j])))
            hits = Int[]
            for combo in ((1,2,3), (1,2,4), (1,3,4), (2,3,4))
                triple = (key[combo[1]], key[combo[2]], key[combo[3]])
                triple in triplesA && continue
                haskey(triB, triple) && push!(hits, triB[triple])
            end
            length(hits) == 2 || continue
            push!(quads, (key, (i, j), (hits[1], hits[2])))
            usedA[i] = true
            usedA[j] = true
            break
        end
    end
    return quads, usedA
end

quads, usedA = structured_quads(body_pos, body_neg)
coverage = 2 * length(quads) / body_pos.ncells
println("\nStructured-quad correspondence: $(length(quads)) quads, covering " *
        "$(round(100 * coverage, digits=2))% of panels")
coverage > 0.99 || @warn "Incomplete quad coverage; uncovered panels excluded from comparison."

# ----------------- QUESTION 1: VELOCITY DECOMPOSITION AT QUAD LEVEL -----------
function half_jump(body; basis=:quad)
    jump = zeros(3, body.ncells)
    pnl.compute_mu_gradient!(jump, body.controlpoints, body.normals,
        body.cells, body.neighbor,
        view(body.strength, :, pnl.get_Gammai(body)),
        te_info(body);
        scale=0.5,
        nodes=basis === :tri ? nothing : body.nodes,
        grad_mu_options=(; basis))
    return jump
end

jump_pos = half_jump(body_pos)
jump_neg = half_jump(body_neg)

# Verify the loaded-velocity reconstruction with a single N-body evaluation on
# the +y mesh, replicating the post-solve influence pass of steady!
# (src/FLOWPanel_simulate.jl): zero velocity, panel convolution, +half-jump,
# +freestream. This must reproduce the loaded velocity.
println("\nVerifying loaded velocity with one N-body evaluation on the +y mesh...")
velocity_loaded_pos = copy(body_pos.velocity)
backend = pnl.FastMultipoleBackend()
body_pos.velocity .= 0
pnl._set_kerneloffsets!((body_pos,), :kerneloffset_targets)
@time pnl.influence!((body_pos,), (body_pos,), backend; precalc=false,
    scalar_potential=false, velocity=true,
    direct_conditioning=pnl._self_panel_kerneloffset_conditioning())
body_pos.velocity .+= jump_pos
pnl.apply_freestream!(body_pos, Vinf)
verr = [_norm(body_pos.velocity[:, i] .- velocity_loaded_pos[:, i]) for i in 1:body_pos.ncells]
println("  recomputed vs loaded velocity: rms=$(round(sqrt(sum(abs2, verr)/length(verr)), sigdigits=4)) " *
        "max=$(round(maximum(verr), sigdigits=4)) m/s (loaded |V| rms=" *
        "$(round(sqrt(sum(abs2, velocity_loaded_pos)/body_pos.ncells), sigdigits=4)))")
body_pos.velocity .= velocity_loaded_pos

pv_pos = body_pos.velocity .- jump_pos .- Vinf
pv_neg = body_neg.velocity .- jump_neg .- Vinf
totalind_pos = body_pos.velocity .- Vinf
totalind_neg = body_neg.velocity .- Vinf

areas_pos = pnl.calc_areas(body_pos.nodes, body_pos.cells)
areas_neg = pnl.calc_areas(body_neg.nodes, body_neg.cells)

function quad_average(field, areas, (i, j))
    f = field isa AbstractVector ? reshape(field, 1, :) : field
    w = areas[i] + areas[j]
    return (areas[i] .* f[:, i] .+ areas[j] .* f[:, j]) ./ w
end

function quad_diffs(fieldA, fieldB)
    diffs = Vector{Float64}(undef, length(quads))
    for (qi, (_, pairA, pairB)) in enumerate(quads)
        diffs[qi] = _norm(quad_average(fieldA, areas_pos, pairA) .-
                          quad_average(fieldB, areas_neg, pairB))
    end
    return diffs
end

stats(d) = (mean=sum(d)/length(d), rms=sqrt(sum(abs2, d)/length(d)), max=maximum(d))

gamma_pos_v = collect(view(body_pos.strength, :, pnl.get_Gammai(body_pos)))
gamma_neg_v = collect(view(body_neg.strength, :, pnl.get_Gammai(body_neg)))

diff_fields = (
    ("gamma",               quad_diffs(gamma_pos_v, gamma_neg_v)),
    ("PV induced velocity", quad_diffs(pv_pos, pv_neg)),
    ("half-jump (grad mu)", quad_diffs(jump_pos, jump_neg)),
    ("total induced vel",   quad_diffs(totalind_pos, totalind_neg)),
)

println("\n#===== QUAD-LEVEL +y vs -y DIFFERENCES (area-weighted quad averages) =====#")
for (label, d) in diff_fields
    s = stats(d)
    println(rpad(label, 22), " mean=$(round(s.mean, sigdigits=5))",
            " rms=$(round(s.rms, sigdigits=5))",
            " max=$(round(s.max, sigdigits=5))")
end

# Spanwise profile of rms differences
quad_y2b = [sum(body_pos.nodes[2, k] for k in key) / 4 * 2 / b for (key, _, _) in quads]
nbins = 24
println("\nSpanwise rms profiles (bin-center 2y/b : gamma, PV, jump, total):")
for bi in 1:nbins
    lo, hi = -1 + 2*(bi-1)/nbins, -1 + 2*bi/nbins
    sel = findall(y -> lo <= y < hi, quad_y2b)
    isempty(sel) && continue
    vals = [sqrt(sum(abs2, d[sel]) / length(sel)) for (_, d) in diff_fields]
    println("  $(rpad(round((lo+hi)/2, digits=3), 7)) : " *
            join((round(v, sigdigits=4) for v in vals), "  "))
end

# ----------------- CONTROLLED SAME-mu GRADIENT EXPERIMENT ---------------------
# Same nodal mu field averaged onto both triangulations: any quad-level gradient
# difference is purely triangulation sensitivity of compute_mu_gradient!.
mu_nodes = cell_scalar_to_nodes(body_pos, gamma_pos_v)
mu_tri(body) = [sum(mu_nodes[body.cells[k, ci]] for k in 1:3) / 3 for ci in axes(body.cells, 2)]
mu_pos = mu_tri(body_pos)
mu_neg = mu_tri(body_neg)

function gradient_of(body, mu; basis)
    g = zeros(3, body.ncells)
    pnl.compute_mu_gradient!(g, body.controlpoints, body.normals,
        body.cells, body.neighbor, mu, te_info(body);
        scale=0.5,
        nodes=basis === :tri ? nothing : body.nodes,
        grad_mu_options=(; basis))
    return g
end

println("\n#===== SAME-mu CONTROLLED GRADIENT EXPERIMENT (scale=0.5) =====#")
for basis in (:quad, :tri)
    gA = gradient_of(body_pos, mu_pos; basis)
    gB = gradient_of(body_neg, mu_neg; basis)
    s = stats(quad_diffs(gA, gB))
    ref = stats(quad_diffs(gA, zero(gA)))
    println("  basis=$(basis): diff mean=$(round(s.mean, sigdigits=5))" *
            " rms=$(round(s.rms, sigdigits=5)) max=$(round(s.max, sigdigits=5))" *
            "   (|grad| rms=$(round(ref.rms, sigdigits=5)))")
end

# ----------------- CL WITH QUAD-CONSISTENT PASS ON vs OFF ---------------------
# The loaded velocity was computed with the quad-consistent pass ON. Toggling it
# off only changes the local half-jump term, so no N-body evaluation is needed:
# U_off = U_loaded - jump_on + jump_off.
rho  = 1.225
Sref = b^2 / ar
Dhat = Vinf ./ magVinf
Lhat = LinearAlgebra.cross(Dhat, [0.0, 1.0, 0.0])

# Replicates the example's monitor stack: PressureBernoulli(rho) (steady, TE
# pair-averaged) -> ForceMonitor(correct_kuttacondition=false) -> WingNormalization.
function compute_CL(body, U, areas)
    P = zeros(body.ncells)
    pnl.calcfield_P!(P, body, U, magVinf, rho, nothing; correct_kuttacondition=true)
    F = zeros(3, body.ncells)
    pnl.calcfield_F!(F, body, areas, body.normals, P; correct_kuttacondition=false)
    CF = vec(sum(F, dims=2)) ./ (0.5 * rho * magVinf^2 * Sref)
    return LinearAlgebra.dot(CF, Lhat)
end

function CL_from_pressure(body, P, areas)
    F = zeros(3, body.ncells)
    pnl.calcfield_F!(F, body, areas, body.normals, P; correct_kuttacondition=false)
    return LinearAlgebra.dot(vec(sum(F, dims=2)) ./ (0.5 * rho * magVinf^2 * Sref), Lhat)
end

jump_off_pos = half_jump(body_pos; basis=:tri)
jump_off_neg = half_jump(body_neg; basis=:tri)
U_off_pos = body_pos.velocity .- jump_pos .+ jump_off_pos
U_off_neg = body_neg.velocity .- jump_neg .+ jump_off_neg

CL_on_pos  = compute_CL(body_pos, body_pos.velocity, areas_pos)
CL_on_neg  = compute_CL(body_neg, body_neg.velocity, areas_neg)
CL_off_pos = compute_CL(body_pos, U_off_pos, areas_pos)
CL_off_neg = compute_CL(body_neg, U_off_neg, areas_neg)

println("\n#===== CL: QUAD-CONSISTENT PASS ON vs OFF =====#")
println("  Sanity (CL from loaded gauge pressure): +y=$(round(CL_from_pressure(body_pos, loaded_pos.pressure, areas_pos), digits=5)), " *
        "-y=$(round(CL_from_pressure(body_neg, loaded_neg.pressure, areas_neg), digits=5))")
println(rpad("  mesh", 12), rpad("quad ON", 12), "quad OFF")
println(rpad("  +y", 12), rpad(round(CL_on_pos, digits=5), 12), round(CL_off_pos, digits=5))
println(rpad("  -y", 12), rpad(round(CL_on_neg, digits=5), 12), round(CL_off_neg, digits=5))
println("  ΔCL(+y vs -y): ON=$(round(abs(CL_on_pos - CL_on_neg), sigdigits=4)) " *
        "OFF=$(round(abs(CL_off_pos - CL_off_neg), sigdigits=4))")

# ----------------- GAMMA-TRANSPLANT CL EXPERIMENT -----------------------------
# Transplant the +y mesh's solved strengths onto the -y mesh and recompute CL
# (one N-body evaluation per case, no solve). If CL(-y, transplanted gamma)
# lands near the +y value, the mesh-to-mesh CL gap is attributable to the
# panel-strength solution rather than post-processing.
#
# Two interpolation variants, each with a +y self-interpolation control so the
# smoothing artifact of the interpolation itself is separated from the actual
# solution difference:
#   - node roundtrip: gamma -> shared nodes (area-weighted) -> target triangles
#     (mean of 3 node values). Caveat: TE nodes are shared by upper and lower
#     surfaces, so this smears the TE circulation jump.
#   - quad average: each target triangle gets the area-weighted gamma average of
#     the two +y triangles in its parent structured quad (never mixes surfaces).

# Node roundtrip reuses mu_nodes = cell_scalar_to_nodes(body_pos, gamma_pos_v)
# from the same-mu experiment.
gamma_nr_pos = mu_tri(body_pos)   # +y self-roundtrip control
gamma_nr_neg = mu_tri(body_neg)   # -y <- +y transplant

gamma_qa_pos = copy(gamma_pos_v)
gamma_qa_neg = similar(gamma_neg_v)
for (_, (i, j), (k, l)) in quads
    qgamma = (areas_pos[i] * gamma_pos_v[i] + areas_pos[j] * gamma_pos_v[j]) /
             (areas_pos[i] + areas_pos[j])
    gamma_qa_pos[i] = qgamma
    gamma_qa_pos[j] = qgamma
    gamma_qa_neg[k] = qgamma
    gamma_qa_neg[l] = qgamma
end

# Nearest-centroid transplant within each structured quad: each -y triangle
# takes the gamma of the closer of the quad's two +y triangles. No averaging at
# all, so the TE jump and the chordwise gradient survive up to a half-panel
# shift; the +y self counterpart is exactly the identity (control = baseline).
gamma_nc_neg = similar(gamma_neg_v)
for (_, (i, j), (k, l)) in quads
    for tk in (k, l)
        di = sum(abs2, body_neg.controlpoints[:, tk] .- body_pos.controlpoints[:, i])
        dj = sum(abs2, body_neg.controlpoints[:, tk] .- body_pos.controlpoints[:, j])
        gamma_nc_neg[tk] = gamma_pos_v[di <= dj ? i : j]
    end
end

# Quantify the TE-jump smearing of the node roundtrip on the +y mesh
te_panels = unique(vcat([vec(shedding[[1, 4], :]) for shedding in body_pos.shedding]...))
filter!(p -> p > 0, te_panels)
rms(v) = sqrt(sum(abs2, v) / length(v))
println("\n#===== GAMMA-TRANSPLANT CL EXPERIMENT =====#")
println("  Node-roundtrip gamma change on +y mesh: rms(all)=" *
        "$(round(rms(gamma_nr_pos .- gamma_pos_v), sigdigits=4)), " *
        "rms(TE panels)=$(round(rms(gamma_nr_pos[te_panels] .- gamma_pos_v[te_panels]), sigdigits=4))")
println("  Transplanted vs own gamma on -y mesh: node-roundtrip rms=" *
        "$(round(rms(gamma_nr_neg .- gamma_neg_v), sigdigits=4)), quad-average rms=" *
        "$(round(rms(gamma_qa_neg .- gamma_neg_v), sigdigits=4)), nearest-centroid rms=" *
        "$(round(rms(gamma_nc_neg .- gamma_neg_v), sigdigits=4))")

function CL_with_gamma(bd, gamma_new, areas)
    Gi = pnl.get_Gammai(bd)
    gamma_stash = copy(view(bd.strength, :, Gi))
    velocity_stash = copy(bd.velocity)
    view(bd.strength, :, Gi) .= gamma_new

    bd.velocity .= 0
    pnl._set_kerneloffsets!((bd,), :kerneloffset_targets)
    pnl.influence!((bd,), (bd,), backend; precalc=false,
        scalar_potential=false, velocity=true,
        direct_conditioning=pnl._self_panel_kerneloffset_conditioning())
    bd.velocity .+= half_jump(bd)
    pnl.apply_freestream!(bd, Vinf)
    CL = compute_CL(bd, bd.velocity, areas)

    view(bd.strength, :, Gi) .= gamma_stash
    bd.velocity .= velocity_stash
    return CL
end

println("  Evaluating 4 transplant cases (one N-body evaluation each)...")
@time CL_nr_pos = CL_with_gamma(body_pos, gamma_nr_pos, areas_pos)
@time CL_nr_neg = CL_with_gamma(body_neg, gamma_nr_neg, areas_neg)
@time CL_qa_pos = CL_with_gamma(body_pos, gamma_qa_pos, areas_pos)
@time CL_qa_neg = CL_with_gamma(body_neg, gamma_qa_neg, areas_neg)
@time CL_nc_neg = CL_with_gamma(body_neg, gamma_nc_neg, areas_neg)

println()
println(rpad("  gamma source", 34), rpad("CL(+y)", 12), rpad("CL(-y)", 12), "ΔCL")
println(rpad("  own solve (baseline)", 34), rpad(round(CL_on_pos, digits=5), 12),
        rpad(round(CL_on_neg, digits=5), 12), round(abs(CL_on_pos - CL_on_neg), sigdigits=4))
println(rpad("  +y gamma, node roundtrip", 34), rpad(round(CL_nr_pos, digits=5), 12),
        rpad(round(CL_nr_neg, digits=5), 12), round(abs(CL_nr_pos - CL_nr_neg), sigdigits=4))
println(rpad("  +y gamma, quad average", 34), rpad(round(CL_qa_pos, digits=5), 12),
        rpad(round(CL_qa_neg, digits=5), 12), round(abs(CL_qa_pos - CL_qa_neg), sigdigits=4))
println(rpad("  +y gamma, nearest centroid", 34), rpad(round(CL_on_pos, digits=5), 12),
        rpad(round(CL_nc_neg, digits=5), 12), round(abs(CL_on_pos - CL_nc_neg), sigdigits=4))

# ----------------- VTK OUTPUT FOR PARAVIEW ------------------------------------
rm(out_path; force=true, recursive=true)
mkpath(out_path)

# Broadcast each quad's diff onto its two +y-mesh triangles so the diffs render
# as a surface in ParaView. Uncovered panels get NaN.
function quad_to_cells(diffs)
    vals = fill(NaN, body_pos.ncells)
    for (qi, (_, (i, j), _)) in enumerate(quads)
        vals[i] = diffs[qi]
        vals[j] = diffs[qi]
    end
    return vals
end

triangle_cells = [collect(c) .- 1 for c in eachcol(body_pos.cells)]
pnl._write_vtk_points_or_lines(joinpath(out_path, "mirror_decomposition_diffs"),
    body_pos.nodes;
    cells=triangle_cells,
    cell_data=Tuple(Dict("field_name" => "diff " * label, "field_data" => quad_to_cells(d))
                    for (label, d) in diff_fields),
    override_cell_type=pnl.WriteVTK.VTKCellTypes.VTK_TRIANGLE)
println("\nSaved quad-level diff surface to $(out_path)/mirror_decomposition_diffs.vtu")
