# fm051_pass_parity.jl — task 051 stage 2 parity harness.
#
# Builds a reduced 018-like configuration (diamond RigidWakeBody{ConstantSource+
# VortexRing} with a finite attached TE wake, a populated PanelWake sheet, and a
# synthetic gaussianerf particle cloud) and evaluates the two cross-influence
# passes of `_steady_aerodynamics!` three ways through their REAL call sites
# (`_sa_wake_influence!` / `_sa_body_influence!`):
#
#   1. unarmed FastMultipoleBackend (the production fmm! path)
#   2. armed host rectangular seam  (FLOWPANEL_GPU_INFLUENCE / set_gpu_influence!)
#   3. unarmed DirectBackend        (brute-force reference)
#
# and reports relative RMS differences per pass, per target output, per filament
# regularization family. Expected gates:
#   armed vs direct : ~1e-12 (same math, different loop order)
#   armed vs fmm    : <= ~1e-3 (fmm! is approximate)
#
# Run (CPU, from the FLOWPanel.jl root):
#   julia --project=. --threads=4 examples/fm051_pass_parity.jl
#
# Task 052: FM051_MODE=cuda arms the seam in :cuda mode instead of :host,
# validating the device direct_rectangular! kernels through the same gates
# (requires functional CUDA; the run FAILS if the seam falls back to host).

armed_mode = lowercase(get(ENV, "FM051_MODE", "host")) == "cuda" ? :cuda : :host

using LinearAlgebra: norm
import Random
import FLOWPanel as pnl
import FLOWPanel.FLOWVPM as vpm

Random.seed!(51)

# ------------------------------------------------------------------ geometry
# diamond wing body (copied from test/test_helpers.jl `make_dirichlet_diamond_body`
# to keep this harness independent of the mesh-file fixtures)
function make_diamond_body(; nspan::Int=3, thick=0.06, das=0.3)
    ys = range(0, 1; length=nspan+1)
    nodes = Float64[]
    for y in ys
        append!(nodes, [0.0, y, 0.0])
        append!(nodes, [0.5, y, thick])
        append!(nodes, [1.0, y, 0.0])
        append!(nodes, [0.5, y, -thick])
    end
    nodes = reshape(nodes, 3, :)
    idx(j, k) = (j-1)*4 + k
    cells = Int[]
    for j in 1:nspan
        le1, up1, te1, lo1 = idx(j, 1), idx(j, 2), idx(j, 3), idx(j, 4)
        le2, up2, te2, lo2 = idx(j+1, 1), idx(j+1, 2), idx(j+1, 3), idx(j+1, 4)
        append!(cells, [le1, up1, up2]); append!(cells, [le1, up2, le2])
        append!(cells, [up1, te1, te2]); append!(cells, [up1, te2, up2])
        append!(cells, [le1, le2, lo2]); append!(cells, [le1, lo2, lo1])
        append!(cells, [lo1, lo2, te2]); append!(cells, [lo1, te2, te1])
    end
    cells = reshape(cells, 3, :)
    shedding = [pnl.calc_shedding_from_seed(nodes, cells, idx(1, 3), idx(2, 3))]
    body = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.VortexRing}}(
        nodes, cells, shedding; check_mesh=false, watertight=false,
        semiinfinite_wake=false)
    for i in eachindex(body.Das)
        body.Das[i] .= repeat([das, 0.0, 0.0], 1, size(body.Das[i], 2))
    end
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    return body
end

nspan = 4
nwakerows = 6
n_particles = 400

body = make_diamond_body(; nspan=nspan)                       # 8*nspan panels
body.strength .= 0.1 .* randn(body.ncells, 2)
body.needs_velocity_gradient[] = true                          # exercise H path

wake = pnl.PanelParticleWake(body; nwakerows=nwakerows, max_particles=8000)
pw = wake.panel_wake

# populate the wake sheet: rows trail downstream of the TE (x = 1) with mild
# span/depth waviness, affine strengths
nod = pw.nodes[1]
str = pw.strength[1]
for irow in 1:size(nod, 2), icol in 1:size(nod, 3)
    nod[1, irow, icol] = 1.0 + 0.35 * (irow - 1)
    nod[2, irow, icol] = (icol - 1) / (size(nod, 3) - 1)
    nod[3, irow, icol] = 0.05 * sin(0.7 * irow + 0.9 * icol)
end
for irow in 1:size(str, 2), icol in 1:size(str, 3)
    str[1, irow, icol] = 0.1 * irow + 0.03 * icol
end
pw.nwakes[] = nwakerows

# synthetic particle cloud around/behind the wake
pfield = wake.pfield
for _ in 1:n_particles
    X = [1.0 + 2.5*rand(), rand(), 0.6*(rand() - 0.5)]
    Gamma = 0.05 .* randn(3)
    sigma = 0.05 + 0.1*rand()
    vpm.add_particle(pfield, X, Gamma, sigma)
end

# assemble the pass tuples exactly as _steady_aerodynamics! does
wake_probes, targets, wake_sources = pnl._sa_collect((body,), (wake,))
println("targets      = ", map(typeof, targets))
println("wake_sources = ", map(typeof, wake_sources))

# ------------------------------------------------------------------ machinery
fmm_backend = pnl.FastMultipoleBackend(; expansion_order=10,
    multipole_acceptance=0.4, leaf_size=20)
direct_backend = pnl.DirectBackend()

function reset_outputs!()
    body.velocity .= 0
    body.velocity_gradient .= 0
    body.potential .= 0
    for v in pw.velocity
        v .= 0
    end
    pfield.particles[vpm.U_INDEX, 1:pfield.np] .= 0
    pfield.particles[vpm.J_INDEX, 1:pfield.np] .= 0
    return nothing
end

snapshot() = (
    bodyU = copy(body.velocity),
    bodyH = copy(body.velocity_gradient),
    wakeU = reduce(vcat, vec.(pw.velocity)),
    partU = copy(pfield.particles[vpm.U_INDEX, 1:pfield.np]),
    partJ = copy(pfield.particles[vpm.J_INDEX, 1:pfield.np]),
)

relrms(a, b) = norm(vec(a) .- vec(b)) / max(norm(vec(b)), eps())

function run_pass1!(backend)
    pnl._sa_wake_influence!(targets, wake_sources, backend)
    return snapshot()
end

function run_pass3!(backend)
    pnl._set_core_sizes!((body,), :core_size_targets)
    pnl._sa_body_influence!(targets, (body,), backend)
    return snapshot()
end

function evaluate(pass!, label)
    pnl.set_gpu_influence!(:off)
    reset_outputs!()
    ref_fmm = pass!(fmm_backend)
    reset_outputs!()
    ref_dir = pass!(direct_backend)
    hits0 = pnl.GPU_INFLUENCE_HITS[]
    pnl.set_gpu_influence!(armed_mode)
    if armed_mode === :cuda
        pnl._gpu_device() === :cuda || error(
            "FM051_MODE=cuda but the CUDA seam is not functional: " *
            pnl.FastMultipole.cuda_radix_status())
    end
    reset_outputs!()
    got = pass!(fmm_backend)
    pnl.set_gpu_influence!(:off)
    armed = pnl.GPU_INFLUENCE_HITS[] - hits0
    armed >= 1 || error("$label: the rectangular seam did not arm (hits=$armed)")
    println("\n== $label (seam handled $armed influence! call(s)) ==")
    for k in keys(got)
        r_dir = relrms(got[k], ref_dir[k])
        r_fmm = relrms(got[k], ref_fmm[k])
        println(rpad(String(k), 6), "  vs direct: ", r_dir, "   vs fmm: ", r_fmm)
    end
    return nothing
end

# ------------------------------------------------------------------ run
for fam in (:vatistas, :compact, :gaussian)
    pnl.set_filament_regularization!(fam)
    println("\n#---------------- filament regularization: $fam ----------------#")
    evaluate(run_pass1!, "pass 1 (wake -> bodies + probes + particles) [$fam]")
    evaluate(run_pass3!, "pass 3 (bodies -> bodies + probes + particles) [$fam]")
end
pnl.set_filament_regularization!(:gaussian)   # restore the working-tree default
println("\ndone.")
