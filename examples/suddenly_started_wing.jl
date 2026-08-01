#=##############################################################################
# DESCRIPTION
#   Wagner-function validation using a suddenly started, open-tip NACA 0012
#   wing and FLOWPanel's unsteady Neumann vortex-ring formulation.
#
#   The reference case is VortexLattice.jl/test/suddenly_started_wing_uvlm.jl.
#   Its executable timestep is U*dt/c = 1/8 (the 1/16 header comment is stale).
#
# RUN
#   SSW_MODE=single julia --project examples/suddenly_started_wing.jl
#   SSW_MODE=convergence julia --project -t auto examples/suddenly_started_wing.jl
#   SSW_MODE=wake_consistency SSW_BACKEND=direct julia --project -t 6 \
#       examples/suddenly_started_wing.jl
#
# ENVIRONMENT
#   SSW_MODE       single | coarse | convergence | wake_consistency
#                  (default: single)
#   SSW_OUTPUT     output root (default: data/suddenly_started_wing)
#   SSW_BACKEND    direct | fmm (default: fmm)
#   SSW_FMM_ORDER  FMM expansion order (default: 10)
#   SSW_FMM_THETA  FMM multipole acceptance threshold (default: 0.4)
#   SSW_FMM_LEAF   FMM source leaf size (default: 100)
#   SSW_NAIRFOIL   contour request for single mode (default: 21)
#   SSW_NSPAN      spanwise strips for single mode (default: 12)
#   SSW_AR         wing aspect ratio (default: 100; use 6 or 3 to compare
#                  against R.T. Jones' finite-AR indicial functions).
#                  WARNING: AR=100 is numerically broken with this mesh — the
#                  steady semi-infinite solve returns cl~1e5 (not ~0.5) at
#                  n_span=12, and still ~1e2 at n_span=48, because the spanwise
#                  panels are several chords wide. Absolute forces there are
#                  meaningless; only ratios of equally-corrupted quantities
#                  survive. wake_consistency mode therefore defaults to AR=6.
#   SSW_ETA        first-wake-row offset Das = eta*dt*U (default: 0.25)
#   SSW_BODYTYPE   neumann | dirichlet (default: neumann)
#   SSW_CAPS       auto | none | flat (default: auto -> flat for dirichlet,
#                  none for neumann; interior Dirichlet assumes a closed
#                  surface, watertight Neumann is rank-deficient)
#   SSW_DT_STAR    U*dt/c for single mode (default: 0.125)
#   SSW_T_END_STAR final tU/c (default: 7)
#   SSW_BACKSLASH_MAX_PANELS dense-solver cutoff (default: 10000)
#   SSW_RESUME     reuse completed convergence histories (default: true)
#   SSW_VERBOSE    true | false (default: true)
#   SAVE_VTK       true | false (default: true; default false in
#                  wake_consistency mode, where 512-step wakes are large)
#   SSW_FREESTREAM_CONVECTION
#                  convect every PanelWake row with the freestream only, so the
#                  sheet stays straight along Uinf and its geometry matches the
#                  semi-infinite wake (default: false)
#   SSW_WAKE_LENGTHS
#                  comma-separated t_end_star values swept by wake_consistency
#                  (default: 1,2,4,8,16,32,64 -> L/c at dt*=0.125)
#   SSW_WAKE_MODES comma-separated freestream_convection flags swept by
#                  wake_consistency (default: false,true)
#   SSW_WAKE_MODEL panel | particle (default: panel)
#   SSW_PANEL_ROWS retained panel rows in particle mode (default: 8)
#   SSW_SHED_METHOD overlap_pps | sigma_pps | sigma_overlap
#   SSW_SIGMA_OVER_C required positive sigma/c for sigma methods
#   SSW_PPS / SSW_OVERLAP particle spacing controls (defaults: 2 / 1.3)
#   SSW_WAKE_CORE wake core/c (default: 0.001)
#   SSW_MAX_PARTICLES particle capacity (default: 200000)
=###############################################################################

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
import LinearAlgebra as LA
import Printf: @printf
get!(ENV, "MPLBACKEND", "Agg")
const SSW_PLOTTING = !(lowercase(get(ENV, "SSW_NO_PLOT", "false")) in
    ("1", "true", "yes", "on"))
SSW_PLOTTING && (@eval import PythonPlot as plt)

# `add_flat_tip_caps` / `open_boundary_loops` (shared with the swept-wing path).
isdefined(@__MODULE__, :add_flat_tip_caps) ||
    include(joinpath(@__DIR__, "helper_functions.jl"))

const DEFAULT_SSW_OUTPUT = joinpath("data", "suddenly_started_wing")
const SSW_GRAD_MU_OPTIONS = (; basis=:tri)

# Mesh revision, carried in the case tag so a mesh change can never silently
# resume histories computed on the previous geometry. Bump this whenever
# `suddenly_started_wing_mesh` changes.
#   1 — pre-2026-07-29: every quad split along the same diagonal. The
#       triangulation was invariant under NEITHER reflection: TE panels were not
#       chordwise mirrors, and all cells violated y -> -y, putting ~10%
#       asymmetry into shed circulation that did not refine away.
#   2 — 2026-07-29: splitting diagonal keyed on the XOR of the chordwise and
#       spanwise half-tests, so it flips under each reflection.
const SSW_MESH_REVISION = 2

Base.@kwdef struct SSWConfig
    c::Float64 = 1.0
    AR::Float64 = 100.0
    U::Float64 = 1.0
    alpha_deg::Float64 = 5.0
    rho::Float64 = 1.0
    t_end_star::Float64 = 7.0
    dt_star::Float64 = 1 / 8
    n_airfoil::Int = 21
    n_span::Int = 12
    eta::Float64 = 0.25
    # Boundary-condition/body formulation. The legacy suddenly-started-wing
    # validation remains Neumann; :dirichlet enables the source+doublet body
    # required by the pressure-continuity Kutta attribution campaign.
    bodytype::Symbol = :neumann
    # Tip closure. `:auto` resolves to :flat for :dirichlet and :none for
    # :neumann (see `ssw_caps`): the interior-Dirichlet Green identity assumes a
    # closed surface, whereas the legacy Neumann validation is an open-tip
    # extrusion and a watertight Neumann body is rank-deficient. Set explicitly
    # to :none or :flat to override either default.
    caps::Symbol = :auto
    # Convect every wake row with the freestream only (no rollup), giving the
    # PanelWake exactly the semi-infinite wake's straight-sheet geometry. This
    # is the control knob for the wake-representation consistency study; the
    # default reproduces the rolled-up (physical) wake.
    freestream_convection::Bool = false
    # Sheet-vs-particle split (BRAINSTORM 014): :panel keeps the pure PanelWake
    # spanning every step; :particle uses a PanelParticleWake whose sheet buffer
    # is `panel_rows` rows, with particles carrying the wake beyond it. In
    # particle mode the final-edge filament stays ACTIVE (constructor default)
    # and unsteady content is shed via OverlapPPS — NoShed would silently delete
    # the starting vortex at the first conversion.
    wake_model::Symbol = :panel
    panel_rows::Int = 8
    max_particles::Int = 200_000
    shed_method::Symbol = :overlap_pps
    sigma_over_c::Float64 = NaN
    pps_overlap::Float64 = 1.3
    pps_n::Int = 2
    wake_core_over_c::Float64 = 1e-3
    kerneloffset_over_c::Float64 = 1e-6
    kernelcutoff_over_c::Float64 = 1e-12
    output_root::String = DEFAULT_SSW_OUTPUT
    save_vtk::Bool = true
    backend_kind::Symbol = :fmm
    fmm_expansion_order::Int = 10
    fmm_multipole_acceptance::Float64 = 0.4
    fmm_leaf_size::Int = 100
    verbose::Bool = true
    # Dense factorization is substantially faster and more robust than the
    # current unpreconditioned matrix-free solve for this extreme-AR geometry.
    backslash_max_panels::Int = 10_000
end

function _ssw_shedding_method(config::SSWConfig)
    config.shed_method in (:overlap_pps, :sigma_pps, :sigma_overlap) ||
        throw(ArgumentError("shed_method must be :overlap_pps, :sigma_pps, or " *
            ":sigma_overlap; got $(config.shed_method)"))
    config.pps_n > 0 || throw(ArgumentError("pps_n must be positive"))
    config.pps_overlap > 0 || throw(ArgumentError("pps_overlap must be positive"))
    if config.shed_method == :overlap_pps
        return pnl.OverlapPPS(config.pps_overlap, config.pps_n)
    end
    isfinite(config.sigma_over_c) && config.sigma_over_c > 0 ||
        throw(ArgumentError("sigma_over_c must be finite and positive for " *
            "$(config.shed_method)"))
    sigma = config.sigma_over_c * config.c
    return config.shed_method == :sigma_pps ?
        pnl.SigmaPPS(sigma, config.pps_n) :
        pnl.SigmaOverlap(sigma, config.pps_overlap)
end

mutable struct SSWRelaxationRecorder
    lock::ReentrantLock
    sum_relative_change::Float64
    max_relative_change::Float64
    count::Int
end
SSWRelaxationRecorder() = SSWRelaxationRecorder(ReentrantLock(), 0.0, 0.0, 0)

function (recorder::SSWRelaxationRecorder)(rlxf::Real, particle)
    gamma = pnl.FLOWVPM.get_Gamma(particle)
    before = copy(gamma)
    pnl.FLOWVPM.relax_correctedpedrizzetti(rlxf, particle)
    change = LA.norm(gamma - before) / max(LA.norm(before), eps(Float64))
    lock(recorder.lock) do
        recorder.sum_relative_change += change
        recorder.max_relative_change = max(recorder.max_relative_change, change)
        recorder.count += 1
    end
    return nothing
end

ssw_relaxation_stats(recorder::SSWRelaxationRecorder) =
    (; mean_relative_change=recorder.count == 0 ? 0.0 :
            recorder.sum_relative_change / recorder.count,
        max_relative_change=recorder.max_relative_change,
        samples=recorder.count)

ssw_recording_relaxation(recorder::SSWRelaxationRecorder=SSWRelaxationRecorder()) =
    (pnl.FLOWVPM.Relaxation(recorder,
        pnl.FLOWVPM.relaxation_correctedpedrizzetti.nsteps_relax,
        pnl.FLOWVPM.relaxation_correctedpedrizzetti.rlxf), recorder)

"""
Jones' two-exponential approximation to Wagner's function. Argument is chord
time t* = t·U/c (Jones' original coefficients 0.0455/0.3 are in semichord time
s = 2t*, hence the doubled exponents here).
"""
ssw_wagner(tstar::Real) = 1 - 0.165exp(-0.09tstar) - 0.335exp(-0.6tstar)

"""
    ssw_caps(config) -> :none | :flat

Resolve `config.caps`. `:auto` (the default) gives flat tip caps to the
Dirichlet body and leaves the Neumann body open, which is the only combination
consistent with both formulations: interior Dirichlet needs a closed surface,
and a watertight Neumann `RigidWakeBody` has a rank-deficient influence matrix
(`FLOWPanel_solver.jl` warns about exactly this).
"""
function ssw_caps(config::SSWConfig)
    config.caps === :auto || config.caps === :none || config.caps === :flat ||
        throw(ArgumentError("caps must be :auto, :none, or :flat; got $(config.caps)"))
    config.caps === :auto || return config.caps
    return config.bodytype === :dirichlet ? :flat : :none
end

# R.T. Jones' one-exponential indicial lift functions for *elliptic* wings of
# finite aspect ratio (NACA Rept. 681, 1940; tabulated in Bisplinghoff et al.,
# "Aeroelasticity", and Leishman). Converted to chord time t* = s/2 like
# ssw_wagner. The SSW mesh is rectangular, not elliptic, so these are an
# external sanity band, NOT an exactness reference — the primary convergence
# criterion for finite-AR studies is self-convergence of the computed curves.
ssw_jones_ar6(tstar::Real) = 1 - 0.361exp(-0.762tstar)
ssw_jones_ar3(tstar::Real) = 1 - 0.283exp(-1.080tstar)

"""
    ssw_reference(config) -> (fun, name)

Select the indicial-lift reference for `config.AR`: Wagner for near-2D aspect
ratios, Jones' elliptic-wing functions near AR=6 and AR=3, otherwise Wagner
with a warning (no published reference at that AR).
"""
function ssw_reference(config::SSWConfig)
    AR = config.AR
    AR >= 20 && return (ssw_wagner, "Wagner (AR=Inf)")
    isapprox(AR, 6.0; atol=0.5) && return (ssw_jones_ar6, "Jones AR=6 (elliptic)")
    isapprox(AR, 3.0; atol=0.5) && return (ssw_jones_ar3, "Jones AR=3 (elliptic)")
    @warn "no indicial reference at AR=$(AR); using Wagner (AR=Inf)" maxlog=1
    return (ssw_wagner, "Wagner (AR=Inf)")
end

function ssw_time_range(config::SSWConfig)
    config.dt_star > 0 || throw(ArgumentError("dt_star must be positive"))
    nsteps = round(Int, config.t_end_star / config.dt_star)
    isapprox(nsteps * config.dt_star, config.t_end_star; atol=100eps(Float64), rtol=0) ||
        throw(ArgumentError("t_end_star must be an integer multiple of dt_star"))
    dt = config.dt_star * config.c / config.U
    return range(0.0; step=dt, length=nsteps + 1)
end

function _ssw_directions(config::SSWConfig)
    alpha = deg2rad(config.alpha_deg)
    drag = SVector{3}(cos(alpha), 0.0, sin(alpha))
    span = SVector{3}(0.0, 1.0, 0.0)
    lift = LA.cross(drag, span)
    return drag, lift
end

"Closed-trailing-edge NACA 0012 contour, ordered TE-lower-LE-upper-TE."
function naca0012_contour(n::Integer=21; thickness::Real=0.12)
    n >= 21 || throw(ArgumentError("n_airfoil must be at least 21"))
    isodd(n) || throw(ArgumentError("n_airfoil must be odd"))
    n_half = cld(n, 2)
    beta = range(0.0, pi; length=n_half)
    x = 0.5 .* (1 .- cos.(beta))
    yt = 5thickness .* (0.2969 .* sqrt.(x) .- 0.1260 .* x .-
        0.3516 .* x.^2 .+ 0.2843 .* x.^3 .- 0.1036 .* x.^4)
    lower = hcat(reverse(x), -reverse(yt))
    upper = hcat(x[2:end-1], yt[2:end-1])
    return vcat(lower, upper)
end

"""
Triangular open-tip extrusion of a NACA 0012 contour.

The node cloud is symmetric under both reflections of this problem: z -> -z
(node `i` <-> node `n_section+2-i`, so quad `i` <-> quad `n_section+1-i`) and
y -> -y (node column `j` <-> column `n_span+2-j`, so quad column `j` <->
column `n_span+1-j`). The *triangulation* must respect both, and it does not
come for free: splitting every quad along the diagonal `n11-n22` is preserved
by neither reflection, because each maps that diagonal onto the opposite
diagonal of the image quad. The consequences are real, not cosmetic:

  * chordwise: the upper- and lower-surface triangles meeting at the trailing
    edge are not reflections, so their control points sit at different
    chordwise stations and the Kutta condition is enforced asymmetrically;
  * spanwise: shed circulation is asymmetric on a symmetric wing by ~10%, and
    because the bias is topological it does NOT refine away (measured 1.01e-1 /
    1.13e-1 / 1.10e-1 at n_span = 12 / 24 / 48).

Each reflection is satisfied by making the split *flip* under it, so both are
satisfied by keying the choice on the XOR of the two half-tests. Both orderings
below traverse the quad loop n11->n21->n22->n12 in the same rotational sense,
so face normals stay consistent.

Exact spanwise mirroring needs an even `n_span`: an odd one leaves the centre
column as its own mirror image, which would have to be its own opposite.

`caps=:flat` closes both tips with a centroid fan (`add_flat_tip_caps`),
producing the watertight surface the interior-Dirichlet formulation assumes.
The fan inherits both reflections from the contour, so the symmetry assertion
still holds. `caps=:none` (the default) reproduces the legacy open-tip mesh
byte for byte, which is why `SSW_MESH_REVISION` does not move.
"""
function suddenly_started_wing_mesh(c::Real, b::Real;
        n_span::Integer=12, n_airfoil::Integer=21, caps::Symbol=:none)
    caps === :none || caps === :flat ||
        throw(ArgumentError("caps must be :none or :flat; got $caps"))
    n_span >= 1 || throw(ArgumentError("n_span must be positive"))
    iseven(n_span) || throw(ArgumentError(
        "n_span must be even so the triangulation can mirror about y=0; got $n_span"))
    contour = naca0012_contour(n_airfoil)
    n_section = size(contour, 1)
    ys = range(-b / 2, b / 2; length=n_span + 1)
    node_index(i, j) = i + (j - 1) * n_section

    nodes = zeros(Float64, 3, n_section * length(ys))
    for (j, y) in enumerate(ys), i in 1:n_section
        nodes[:, node_index(i, j)] .= (c * contour[i, 1], y, c * contour[i, 2])
    end

    cells = zeros(Int, 3, 2 * n_section * n_span)
    k = 0
    for j in 1:n_span, i in 1:n_section
        ip = i == n_section ? 1 : i + 1
        n11, n21 = node_index(i, j), node_index(ip, j)
        n12, n22 = node_index(i, j + 1), node_index(ip, j + 1)
        # Flip the splitting diagonal under each reflection: XOR of the
        # chordwise and spanwise half-tests.
        if xor(2i <= n_section, 2j <= n_span)
            cells[:, k += 1] .= (n11, n21, n22)
            cells[:, k += 1] .= (n11, n22, n12)
        else
            cells[:, k += 1] .= (n11, n21, n12)
            cells[:, k += 1] .= (n21, n22, n12)
        end
    end
    caps === :flat || return nodes, cells
    capped_nodes, capped_cells, _, _ = add_flat_tip_caps(nodes, cells)
    return capped_nodes, capped_cells
end

"""
    assert_ssw_mesh_symmetry(c, b; n_span, n_airfoil, caps)

Assert that the triangulation — not merely the node cloud — is invariant under
both reflections of this problem. A symmetric node cloud with an asymmetric
connectivity silently biases the solution: it produced a ~10% spanwise
asymmetry in shed circulation that did not refine away. O(ncells).

With `caps=:flat` the two appended centroid nodes are checked geometrically
(each lies on the z=0 plane and the two are y-reflections of each other, within
`cap_tol`) and then folded into the same face-set reflection test, so cap faces
must mirror exactly like surface faces.
"""
function assert_ssw_mesh_symmetry(c::Real, b::Real; n_span::Integer=12,
        n_airfoil::Integer=21, caps::Symbol=:none,
        cap_tol::Real=1e-10 * max(abs(c), abs(b)))
    nodes, cells = suddenly_started_wing_mesh(c, b; n_span, n_airfoil, caps)
    n_section = size(naca0012_contour(n_airfoil), 1)
    n_base = n_section * (n_span + 1)
    # node index <-> (i, j)
    ij(n) = (mod1(n, n_section), cld(n, n_section))
    idx(i, j) = i + (j - 1) * n_section

    # Cap centroids are appended after the structured block; pair them by y.
    cap_ids = collect((n_base + 1):size(nodes, 2))
    cap_partner = Dict{Int, Int}()
    if !isempty(cap_ids)
        length(cap_ids) == 2 || error("expected two tip caps, got $(length(cap_ids))")
        lo, hi = cap_ids[argmin(nodes[2, cap_ids])], cap_ids[argmax(nodes[2, cap_ids])]
        cap_partner[lo], cap_partner[hi] = hi, lo
        for n in cap_ids
            abs(nodes[3, n]) <= cap_tol || error("tip-cap centroid $n is off the " *
                "z=0 plane by $(abs(nodes[3, n])) > $cap_tol; the cap breaks z -> -z")
        end
        for d in (1, 3)
            abs(nodes[d, lo] - nodes[d, hi]) <= cap_tol ||
                error("tip-cap centroids differ in coordinate $d by " *
                      "$(abs(nodes[d, lo] - nodes[d, hi])) > $cap_tol; " *
                      "the caps break y -> -y")
        end
        abs(nodes[2, lo] + nodes[2, hi]) <= cap_tol ||
            error("tip-cap centroids are not y-reflections of each other")
    end

    # i=1 is the closed TE node, its own mirror; likewise the LE node.
    mirror_z(n) = n > n_base ? n :
        ((i, j) = ij(n); idx(i == 1 ? 1 : n_section + 2 - i, j))
    mirror_y(n) = n > n_base ? cap_partner[n] :
        ((i, j) = ij(n); idx(i, n_span + 2 - j))

    faces = Set(Tuple(sort(cells[:, k])) for k in axes(cells, 2))
    for (name, m) in (("z -> -z", mirror_z), ("y -> -y", mirror_y))
        for k in axes(cells, 2)
            image = Tuple(sort([m(n) for n in cells[:, k]]))
            image in faces || error("suddenly_started_wing_mesh is not symmetric " *
                "under $name: cell $k maps to a triangle that is not in the mesh")
        end
    end
    return true
end

function _ssw_shedding(nodes, cells, c)
    tol = max(1e-8 * abs(c), 100eps(Float64) * max(abs(c), 1.0))
    te_nodes = findall(i -> isapprox(nodes[1, i], c; atol=tol, rtol=0), axes(nodes, 2))
    sort!(te_nodes; by=i -> nodes[2, i])
    length(te_nodes) >= 2 || error("could not identify the trailing-edge chain")
    lower = [c - tol, minimum(nodes[2, te_nodes]) - tol, minimum(nodes[3, te_nodes]) - tol]
    upper = [c + tol, maximum(nodes[2, te_nodes]) + tol, maximum(nodes[3, te_nodes]) + tol]
    return pnl.calc_shedding_from_seed(nodes, cells, te_nodes[1], te_nodes[2];
        end_node=te_nodes[end], bbox=(lower, upper), normal_jump_tol=1.0,
        max_turn_angle=pi / 2)
end

function build_suddenly_started_wing(config::SSWConfig; semiinfinite_wake::Bool=false)
    b = config.AR * config.c
    nodes, cells = suddenly_started_wing_mesh(config.c, b;
        n_span=config.n_span, n_airfoil=config.n_airfoil, caps=ssw_caps(config))
    bodytype = if config.bodytype == :neumann
        pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}
    elseif config.bodytype == :dirichlet
        pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.ConstantDoublet}}
    else
        throw(ArgumentError(
            "bodytype must be :neumann or :dirichlet; got $(config.bodytype)"))
    end
    # Uncapped meshes are open by construction; capped ones must actually close,
    # and `watertight=true` is what makes the constructor orient the component
    # outward rather than merely consistently.
    watertight, _ = pnl.iswatertight(nodes, cells)
    ssw_caps(config) === :flat && !watertight &&
        error("flat tip caps did not close the surface; check add_flat_tip_caps")
    options = (;
        kerneloffset=config.kerneloffset_over_c * config.c,
        kernelcutoff=config.kernelcutoff_over_c * config.c,
        semiinfinite_wake,
        watertight,
    )

    # The constructor may re-wind cells. Derive shedding from the constructed
    # base body, never from the raw mesh.
    base = bodytype(nodes, cells, zeros(Int, 6, 0); options...)
    shedding = _ssw_shedding(base.nodes, base.cells, config.c)
    return bodytype(copy(base.nodes), copy(base.cells), [shedding]; options...)
end

function _set_ssw_Das!(body, displacement)
    for Das in body.Das, j in axes(Das, 2)
        Das[:, j] .= displacement
    end
    return body
end

function _ssw_backend(config::SSWConfig)
    config.backend_kind == :direct && return pnl.DirectBackend()
    if config.backend_kind == :fmm
        config.fmm_expansion_order > 0 ||
            throw(ArgumentError("fmm_expansion_order must be positive"))
        config.fmm_multipole_acceptance >= 0 ||
            throw(ArgumentError("fmm_multipole_acceptance must be nonnegative"))
        config.fmm_leaf_size > 0 || throw(ArgumentError("fmm_leaf_size must be positive"))
        return pnl.FastMultipoleBackend(
            expansion_order=config.fmm_expansion_order,
            multipole_acceptance=config.fmm_multipole_acceptance,
            leaf_size=config.fmm_leaf_size,
        )
    end
    throw(ArgumentError(
        "backend_kind must be :direct or :fmm; got $(config.backend_kind)"))
end

function _ssw_solver(body, backend, max_backslash::Int)
    if body.ncells <= max_backslash
        return pnl.Backslash(body), :backslash
    end
    return pnl.KrylovSolver(body; backend, method=:gmres, atol=1e-9,
        rtol=1e-9, itmax=1000), :krylov
end

function prepare_suddenly_started_wing(config::SSWConfig; relaxation=nothing,
        particle_maintenance=pnl.ParticleMaintenance())
    t_range = ssw_time_range(config)
    drag_hat, lift_hat = _ssw_directions(config)
    full_uinf = config.U * drag_hat
    Uinf(t) = t <= first(t_range) ? zero(full_uinf) : full_uinf
    maneuver!(frames, systems, wakes, t) = nothing

    wing = build_suddenly_started_wing(config; semiinfinite_wake=false)
    dt = step(t_range)
    _set_ssw_Das!(wing, config.eta * dt * full_uinf)
    wake = if config.wake_model == :particle
        config.panel_rows >= 1 || throw(ArgumentError(
            "panel_rows must be >= 1 in particle mode"))
        shedding_method = _ssw_shedding_method(config)
        # include_final_filament stays at its default (true): the final-edge
        # filament carries the interface circulation where the sheet hands off
        # to particles (unsteady_filament=true pairing, matching the
        # pitching_wing.jl particle branch).
        particle_kwargs = isnothing(relaxation) ? (;) : (; relaxation)
        pnl.PanelParticleWake(wing;
            nwakerows=config.panel_rows,
            max_particles=config.max_particles,
            core_size=config.wake_core_over_c * config.c,
            shed_with_induced_velocity=true,
            freestream_convection=config.freestream_convection,
            method_trailing=shedding_method,
            method_unsteady=shedding_method, particle_kwargs...,
            particle_maintenance)
    elseif config.wake_model == :panel
        pnl.PanelWake(wing;
            nwakerows=length(t_range) - 1,
            core_size=config.wake_core_over_c * config.c,
            include_final_filament=false,
            shed_with_induced_velocity=true,
            freestream_convection=config.freestream_convection)
    else
        throw(ArgumentError("wake_model must be :panel or :particle"))
    end
    frames = pnl.ReferenceFrame(wing)
    backend = _ssw_backend(config)
    solver, solver_name = _ssw_solver(wing, backend, config.backslash_max_panels)
    # Particles (and the active final filament) carry no scalar potential, so
    # the unsteady Bernoulli must be allowed to run on partial potential there.
    pressure = pnl.PressureBernoulli(config.rho; unsteady=true,
        backend, correct_kuttacondition=false,
        allow_partial=(config.wake_model == :particle))
    force = pnl.ForceMonitor(length(t_range), 1;
        normalization=pnl.NoNormalization(), i_frame=-1,
        correct_kuttacondition=false, verbose=false)
    spanwise = pnl.SpanwiseLoadingMonitor(config.n_span, 1;
        components=(lift=lift_hat, drag=drag_hat),
        span_axis=(0.0, 1.0, 0.0), per_length=true,
        normalization=pnl.NoSectionalNormalization(), file=false,
        binning=:span_overlap, verbose=false)
    circulation = pnl.BoundCirculationMonitor(wing, length(t_range), 1;
        i_frame=1, radial_dimension=2, R=config.AR * config.c / 2,
        file=false, verbose=false)
    spanwise_history = Vector{Float64}[]
    function record_spanwise!(systems, wakes, frames, uinf, i_step, dt)
        push!(spanwise_history, copy(spanwise.load_components[1, :]))
        return nothing
    end

    # FLOWPanel solves at every sample, including t=0. Prime Bernoulli history
    # with a genuinely zero-flow body solve, then give the unattached wake row
    # the first interval's full convection velocity. After propagation/shedding,
    # update_TE! resets the new upstream edge while the downstream edge retains
    # this displacement, avoiding a degenerate first wake panel.
    function initialize_wake_convection!(systems, wakes, frames, uinf, i_step, dt)
        i_step == 0 || return nothing
        panel_wake = wakes[1] isa pnl.PanelParticleWake ?
            wakes[1].panel_wake : wakes[1]
        if config.freestream_convection
            # This mode ignores `velocity` entirely and translates every row by
            # dt*wake.freestream, which is still zero at t=0. Prime the
            # freestream instead, or the first shed panel collapses to zero area
            # and that degenerate row is convected downstream forever.
            panel_wake.freestream .= full_uinf
        else
            for velocity in panel_wake.velocity
                velocity[:, 1, :] .= full_uinf
            end
        end
        return nothing
    end
    monitors = (pressure, force, spanwise, circulation, record_spanwise!,
        initialize_wake_convection!)
    return (; wing, wake, frames, maneuver!, Uinf, full_uinf, t_range,
        drag_hat, lift_hat, backend, solver, solver_name, pressure, force,
        spanwise, circulation, spanwise_history, record_spanwise!,
        initialize_wake_convection!, monitors)
end

function _ssw_steady_cl(config::SSWConfig, backend; path=nothing, name="steady",
        grad_mu_options=SSW_GRAD_MU_OPTIONS)
    wing = build_suddenly_started_wing(config; semiinfinite_wake=true)
    drag_hat, lift_hat = _ssw_directions(config)
    full_uinf = config.U * drag_hat
    _set_ssw_Das!(wing, full_uinf / LA.norm(full_uinf))
    solver, solver_name = _ssw_solver(wing, backend, config.backslash_max_panels)
    pressure = pnl.PressureBernoulli(config.rho; correct_kuttacondition=false,
        backend)
    force = pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization(),
        i_frame=-1, correct_kuttacondition=false, verbose=false)
    elapsed = @elapsed pnl.steady!(wing, pnl.ReferenceFrame(wing), full_uinf;
        body_solvers=solver, backend, monitors=(pressure, force), path, name,
        grad_mu_options, verbose=false)
    qS = 0.5 * config.rho * config.U^2 * config.AR * config.c^2
    cl = LA.dot(force.force[:, 1], lift_hat) / qS
    cd = LA.dot(force.force[:, 1], drag_hat) / qS
    return (; cl, cd, elapsed, solver_name, wing)
end

_ssw_num_tag(x) = replace(string(isinteger(x) ? Int(x) : x), "." => "p", "-" => "m")

function _ssw_case_tag(config::SSWConfig)
    # AR and eta are part of the tag so sweeps over them do not collide in the
    # resume cache. Histories written before these tokens existed will simply
    # not be found and will re-run.
    tag = "ar$(_ssw_num_tag(config.AR))_eta$(_ssw_num_tag(config.eta))_" *
        "na$(config.n_airfoil)_ns$(config.n_span)_" *
        "dt$(_ssw_num_tag(config.dt_star))_$(config.backend_kind)_" *
        "$(config.bodytype)"
    # Only appended when non-default, so tags of previously-run cases are
    # unchanged and their histories still resume.
    config.t_end_star == 7.0 || (tag *= "_t$(_ssw_num_tag(config.t_end_star))")
    # The mesh itself changes with caps, so this token is mandatory; :none
    # reproduces every pre-caps tag.
    ssw_caps(config) === :none || (tag *= "_caps$(ssw_caps(config))")
    config.freestream_convection && (tag *= "_fc")
    if config.wake_model == :particle
        _ssw_shedding_method(config)
        tag *= "_pw$(config.panel_rows)"
        config.max_particles == 200_000 ||
            (tag *= "_mp$(config.max_particles)")
        if config.shed_method == :overlap_pps
            tag *= "_opp"
        elseif config.shed_method == :sigma_pps
            tag *= "_spp_sig$(_ssw_num_tag(config.sigma_over_c))"
        elseif config.shed_method == :sigma_overlap
            tag *= "_sov_sig$(_ssw_num_tag(config.sigma_over_c))"
        end
        tag *= "_ov$(_ssw_num_tag(config.pps_overlap))_pn$(config.pps_n)"
    end
    config.wake_core_over_c == 1e-3 ||
        (tag *= "_core$(_ssw_num_tag(config.wake_core_over_c))")
    # Always appended: a history is only comparable to a run on the same mesh.
    tag *= "_m$(SSW_MESH_REVISION)"
    return tag
end

function ssw_settled_stats(values, time_star; window=1.0)
    length(values) == length(time_star) ||
        throw(ArgumentError("values and time_star must have equal length"))
    isempty(values) && throw(ArgumentError("settled statistics require samples"))
    selected = findall(>=(maximum(time_star) - window), time_star)
    tail = Float64.(values[selected])
    scale = max(abs(sum(tail) / length(tail)), eps(Float64))
    x = Float64.(time_star[selected])
    x0 = sum(x) / length(x)
    y0 = sum(tail) / length(tail)
    denom = sum(abs2, x .- x0)
    slope = iszero(denom) ? 0.0 : sum((x .- x0) .* (tail .- y0)) / denom
    drift = abs(slope) * window / scale
    ripple = (maximum(tail) - minimum(tail)) / scale
    return (; mean=y0, min=minimum(tail), max=maximum(tail), drift, ripple,
        settled=drift <= 1e-3 && ripple <= 2.5e-3, indices=selected)
end

function _ssw_station_tail_stats(values::AbstractMatrix, indices)
    return map(axes(values, 1)) do station
        tail = view(values, station, indices)
        meanval = sum(tail) / length(tail)
        scale = max(abs(meanval), eps(Float64))
        (; mean=meanval, min=minimum(tail), max=maximum(tail),
            ripple=(maximum(tail) - minimum(tail)) / scale)
    end
end

function ssw_matrix_settled_stats(values::AbstractMatrix, time_star; window=1.0)
    size(values, 2) == length(time_star) ||
        throw(ArgumentError("matrix history columns must match time_star"))
    selected = findall(>=(maximum(time_star) - window), time_star)
    tail = values[:, selected]
    means = vec(sum(tail; dims=2) ./ length(selected))
    scale = max(maximum(abs, means), eps(Float64))
    x = Float64.(time_star[selected])
    x0 = sum(x) / length(x)
    denom = sum(abs2, x .- x0)
    slopes = iszero(denom) ? zeros(size(values, 1)) :
        vec(sum((tail .- reshape(means, :, 1)) .*
            reshape(x .- x0, 1, :); dims=2)) ./ denom
    drift = maximum(abs, slopes) * window / scale
    ripple = maximum(maximum(view(tail, i, :)) - minimum(view(tail, i, :))
        for i in axes(tail, 1)) / scale
    return (; drift, ripple, settled=drift <= 1e-3 && ripple <= 2.5e-3,
        indices=selected)
end

function _write_ssw_tail_artifacts(result)
    mkpath(result.path)
    gamma_path = joinpath(result.path, "gamma_te.csv")
    gamma_stats = _ssw_station_tail_stats(result.gamma_te, result.tail_indices)
    open(gamma_path, "w") do io
        println(io, "y_over_semispan,mean,min,max,ripple")
        for i in eachindex(result.gamma_y_over_semispan)
            s = gamma_stats[i]
            @printf(io, "%.16e,%.16e,%.16e,%.16e,%.16e\n",
                result.gamma_y_over_semispan[i], s.mean, s.min, s.max, s.ripple)
        end
    end
    loading_path = joinpath(result.path, "spanwise_loading.csv")
    loading_stats = _ssw_station_tail_stats(result.spanwise_loading,
        result.tail_indices)
    open(loading_path, "w") do io
        println(io, "y_over_semispan,mean,min,max,ripple")
        for i in eachindex(result.loading_y_over_semispan)
            s = loading_stats[i]
            @printf(io, "%.16e,%.16e,%.16e,%.16e,%.16e\n",
                result.loading_y_over_semispan[i], s.mean, s.min, s.max, s.ripple)
        end
    end
    settledness_path = joinpath(result.path, "settledness.csv")
    open(settledness_path, "w") do io
        println(io, "quantity,drift,ripple,settled")
        println(io, "CL,$(result.tail_CL.drift),$(result.tail_CL.ripple),$(result.tail_CL.settled)")
        println(io, "gamma_te,$(result.tail_gamma.drift),$(result.tail_gamma.ripple),$(result.tail_gamma.settled)")
    end
    diagnostics_path = joinpath(result.path, "case_diagnostics.csv")
    relaxation = :relaxation_stats in keys(result) ? result.relaxation_stats :
        (; mean_relative_change=NaN, max_relative_change=NaN, samples=0)
    open(diagnostics_path, "w") do io
        println(io, "wake_rows,n_particles,wake_extent_over_c,relaxation_mean_relative_change,relaxation_max_relative_change,relaxation_samples")
        println(io, "$(result.wake_rows),$(result.n_particles),$(result.wake_extent_over_c),$(relaxation.mean_relative_change),$(relaxation.max_relative_change),$(relaxation.samples)")
    end
    return (; gamma_path, loading_path, settledness_path, diagnostics_path)
end

function ssw_with(config::SSWConfig; kwargs...)
    names = fieldnames(SSWConfig)
    values = ntuple(i -> getfield(config, names[i]), length(names))
    base = NamedTuple{names}(values)
    return SSWConfig(; merge(base, (; kwargs...))...)
end

function _write_ssw_case_csv(result)
    mkpath(result.path)
    filename = joinpath(result.path, "history.csv")
    open(filename, "w") do io
        println(io, "time_star,CL,CD,CL_over_CLsteady,wagner,error")
        for i in eachindex(result.time_star)
            @printf(io, "%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                result.time_star[i], result.CL[i], result.CD[i], result.ratio[i],
                result.wagner[i], result.ratio[i] - result.wagner[i])
        end
    end
    return filename
end

function load_ssw_result(config::SSWConfig)
    tag = _ssw_case_tag(config)
    path = joinpath(config.output_root, tag)
    csv_path = joinpath(path, "history.csv")
    isfile(csv_path) || throw(ArgumentError("no completed history at $csv_path"))
    rows = [parse.(Float64, split(line, ',')) for line in readlines(csv_path)[2:end]]
    isempty(rows) && error("empty history at $csv_path")
    time_star = getindex.(rows, 1)
    CL = getindex.(rows, 2)
    CD = getindex.(rows, 3)
    ratio = getindex.(rows, 4)
    wagner = getindex.(rows, 5)
    error = ratio .- wagner
    steady_candidates = [cl / r for (cl, r) in zip(CL, ratio)
        if isfinite(cl) && isfinite(r) && abs(r) > eps(Float64)]
    steady_CL = isempty(steady_candidates) ? NaN : first(steady_candidates)
    tail_CL = ssw_settled_stats(CL, time_star)
    function read_tail_artifact(name)
        artifact = joinpath(path, name)
        isfile(artifact) || return nothing
        data = [parse.(Float64, split(line, ',')) for line in readlines(artifact)[2:end]]
        return (; coordinate=getindex.(data, 1), mean=getindex.(data, 2),
            min=getindex.(data, 3), max=getindex.(data, 4),
            ripple=getindex.(data, 5), path=artifact)
    end
    gamma_artifact = read_tail_artifact("gamma_te.csv")
    loading_artifact = read_tail_artifact("spanwise_loading.csv")
    settledness_path = joinpath(path, "settledness.csv")
    tail_gamma = if isfile(settledness_path)
        fields = split(readlines(settledness_path)[3], ',')
        (; drift=parse(Float64, fields[2]), ripple=parse(Float64, fields[3]),
            settled=lowercase(fields[4]) == "true")
    else
        (; drift=Inf, ripple=Inf, settled=false)
    end
    diagnostics_path = joinpath(path, "case_diagnostics.csv")
    diagnostic_fields = isfile(diagnostics_path) ?
        split(readlines(diagnostics_path)[2], ',') : String[]
    n_particles = length(diagnostic_fields) >= 2 ?
        parse(Int, diagnostic_fields[2]) : -1
    wake_extent_over_c = length(diagnostic_fields) >= 3 ?
        parse(Float64, diagnostic_fields[3]) : config.t_end_star
    wake_extent_estimated = length(diagnostic_fields) < 3
    relaxation_stats = length(diagnostic_fields) >= 6 ?
        (; mean_relative_change=parse(Float64, diagnostic_fields[4]),
            max_relative_change=parse(Float64, diagnostic_fields[5]),
            samples=parse(Int, diagnostic_fields[6])) :
        (; mean_relative_change=NaN, max_relative_change=NaN, samples=0)
    panels = 2 * (config.n_airfoil - 1) * config.n_span
    solver_name = panels <= config.backslash_max_panels ? :backslash : :krylov
    return (; config, tag, path, panels, solver_name, time_star, CL, CD, ratio,
        wagner, steady_CL, steady_CD=NaN, tail_CL, tail_gamma,
        gamma_artifact, loading_artifact,
        rms_error=LA.norm(error) / sqrt(length(error)),
        max_error=maximum(abs, error), elapsed=NaN, steady_elapsed=NaN,
        wake_rows=length(time_star), n_particles, wake_extent_over_c,
        wake_extent_estimated, relaxation_stats, csv_path)
end

function run_suddenly_started_wing(config::SSWConfig=SSWConfig(); relaxation=nothing,
        relaxation_recorder=nothing,
        particle_maintenance=pnl.ParticleMaintenance())
    sim = prepare_suddenly_started_wing(config; relaxation, particle_maintenance)
    tag = _ssw_case_tag(config)
    case_path = joinpath(config.output_root, tag)
    vtk_path = config.save_vtk ? case_path : nothing

    config.verbose && println("Suddenly-started wing: $tag, panels=$(sim.wing.ncells), solver=$(sim.solver_name)")
    steady = _ssw_steady_cl(config, sim.backend;
        path=vtk_path, name=tag * "_steady")
    elapsed = @elapsed pnl.simulate!((sim.wing,), (sim.wake,), sim.frames,
        sim.maneuver!, sim.Uinf, sim.t_range;
        body_solvers=(sim.solver,), backend=sim.backend,
        monitors=sim.monitors, path=vtk_path, name=tag,
        set_Das_eta_freestream=NaN, grad_mu_options=SSW_GRAD_MU_OPTIONS,
        verbose=config.verbose)

    qS = 0.5 * config.rho * config.U^2 * config.AR * config.c^2
    CL_all = vec(sim.lift_hat' * sim.force.force) ./ qS
    CD_all = vec(sim.drag_hat' * sim.force.force) ./ qS
    # VortexLattice reports one result per completed interval, at dt:dt:t_end.
    indices = 2:length(sim.t_range)
    time_star = collect(sim.t_range[indices]) .* config.U ./ config.c
    CL = CL_all[indices]
    CD = CD_all[indices]
    ratio = CL ./ steady.cl
    ref_fun, ref_name = ssw_reference(config)
    wagner = ref_fun.(time_star)   # reference curve (Wagner or Jones finite-AR)
    error = ratio .- wagner
    rms_error = LA.norm(error) / sqrt(length(error))
    max_error = maximum(abs, error)
    tail_CL = ssw_settled_stats(CL, time_star)
    monitor_indices = indices
    gamma_te = Matrix(sim.circulation.circulation_te[:, 1, monitor_indices])
    spanwise_loading = hcat(sim.spanwise_history[monitor_indices]...)
    gamma_y_over_semispan = vec(sim.circulation.r_over_R[:, 1])
    loading_y_over_semispan =
        sim.spanwise.bin_center ./ (config.AR * config.c / 2)
    tail_gamma = ssw_matrix_settled_stats(gamma_te, time_star)
    body_te = maximum(LA.dot(SVector{3}(sim.wing.nodes[:, i]), sim.drag_hat)
        for i in axes(sim.wing.nodes, 2))
    panel_wake = sim.wake isa pnl.PanelParticleWake ?
        sim.wake.panel_wake : sim.wake
    wake_extent = 0.0
    for nodes in panel_wake.nodes, row in 1:panel_wake.nwakes[] + 1,
            column in axes(nodes, 3)
        wake_extent = max(wake_extent,
            LA.dot(SVector{3}(nodes[:, row, column]), sim.drag_hat) - body_te)
    end
    if sim.wake isa pnl.PanelParticleWake
        for i in 1:sim.wake.pfield.np
            wake_extent = max(wake_extent,
                LA.dot(SVector{3}(pnl.FLOWVPM.get_X(sim.wake.pfield, i)),
                    sim.drag_hat) - body_te)
        end
    end
    result = (; config, tag, path=case_path, panels=sim.wing.ncells,
        solver_name=sim.solver_name, time_star, CL, CD, ratio, wagner,
        reference_name=ref_name,
        steady_CL=steady.cl, steady_CD=steady.cd, rms_error, max_error,
        elapsed, steady_elapsed=steady.elapsed,
        tail_CL, tail_gamma, tail_indices=tail_CL.indices, gamma_te,
        gamma_y_over_semispan, spanwise_loading, loading_y_over_semispan,
        wake_rows=(sim.wake isa pnl.PanelParticleWake ?
            sim.wake.panel_wake : sim.wake).nwakes[],
        n_particles=sim.wake isa pnl.PanelParticleWake ? sim.wake.pfield.np : 0,
        wake_extent_over_c=wake_extent / config.c)
    if !isnothing(relaxation_recorder)
        result = merge(result, (; relaxation_stats=
            ssw_relaxation_stats(relaxation_recorder)))
    end
    csv_path = _write_ssw_case_csv(result)
    artifacts = _write_ssw_tail_artifacts(result)

    @printf("  CLsteady=%+.8f, %s RMS=%8.4e, max=%8.4e, elapsed=%.2fs\n",
        result.steady_CL, ref_name, result.rms_error, result.max_error,
        result.elapsed)
    return merge(result, (; csv_path), artifacts)
end

function ssw_curve_change(coarse, fine)
    coarse.time_star == fine.time_star || throw(ArgumentError("curve times differ"))
    delta = fine.ratio .- coarse.ratio
    scale = max.(abs.(fine.ratio), sqrt(eps(Float64)))
    return (; max_absolute=maximum(abs, delta),
        max_relative=maximum(abs.(delta) ./ scale),
        relative_l2=LA.norm(delta) / max(LA.norm(fine.ratio), eps(Float64)))
end

function _write_ssw_summary(results, output_root)
    mkpath(output_root)
    filename = joinpath(output_root, "convergence.csv")
    open(filename, "w") do io
        println(io, "tag,backend,n_airfoil,n_span,dt_star,panels,solver,CLsteady,wagner_rms,wagner_max,successive_max_absolute,successive_max_relative,successive_relative_l2,elapsed")
        previous = Dict{Tuple{Symbol, Float64}, Any}()
        for r in sort(collect(results); by=x -> (string(x.config.backend_kind), x.config.dt_star, x.panels))
            key = (r.config.backend_kind, r.config.dt_star)
            change = haskey(previous, key) && previous[key].time_star == r.time_star ?
                ssw_curve_change(previous[key], r) :
                (; max_absolute=NaN, max_relative=NaN, relative_l2=NaN)
            @printf(io, "%s,%s,%d,%d,%.8g,%d,%s,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.8g\n",
                r.tag, string(r.config.backend_kind), r.config.n_airfoil,
                r.config.n_span, r.config.dt_star, r.panels, string(r.solver_name),
                r.steady_CL, r.rms_error, r.max_error, change.max_absolute,
                change.max_relative, change.relative_l2, r.elapsed)
            previous[key] = r
        end
    end
    return filename
end

function plot_ssw_results(results; output_root=first(results).config.output_root,
        emphasize_tag=nothing)
    isempty(results) && return nothing
    mkpath(output_root)
    fig = plt.figure("ssw_wagner", figsize=(11, 4.5))
    fig.clear()
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    tfine = range(0.0, maximum(first(results).time_star); length=400)
    ref_fun, ref_name = ssw_reference(first(results).config)
    ax1.plot(tfine, ref_fun.(tfine), "k:"; linewidth=2, label=ref_name)
    for (i, r) in enumerate(results)
        label = "$(r.config.n_airfoil)/$(r.config.n_span), $(r.config.backend_kind), dt*=$(r.config.dt_star), eta=$(r.config.eta)"
        emphasized = isnothing(emphasize_tag) ? i == length(results) : r.tag == emphasize_tag
        ax1.plot(r.time_star, r.ratio; linewidth=emphasized ? 2.2 : 1.0,
            label)
        ax2.plot(r.time_star, r.ratio .- r.wagner;
            linewidth=emphasized ? 2.2 : 1.0, label)
    end
    ax1.set(xlabel="t*=t Uinf/c", ylabel="CL/CLsteady",
        title="Suddenly-started NACA 0012 wing")
    ax2.set(xlabel="t*=t Uinf/c", ylabel="FLOWPanel - reference",
        title="Validation error vs $ref_name")
    for ax in (ax1, ax2)
        ax.grid(true; alpha=0.3)
        ax.legend(fontsize=7, loc="best")
    end
    fig.tight_layout()
    filename = joinpath(output_root, "wagner_comparison.png")
    fig.savefig(filename; dpi=180, bbox_inches="tight")
    return filename
end

function plot_ssw_convergence(results; output_root=first(results).config.output_root,
        backend_kind=first(results).config.backend_kind,
        dt_star=first(results).config.dt_star)
    spatial = [r for r in results if r.config.backend_kind == backend_kind &&
        r.config.dt_star == dt_star]
    sort!(spatial; by=r -> r.panels)
    # Keep only the finest-dimensional case at duplicate panel counts.
    spatial = [spatial[i] for i in eachindex(spatial)
        if i == length(spatial) || spatial[i + 1].panels != spatial[i].panels]
    isempty(spatial) && return nothing
    panels = getproperty.(spatial, :panels)
    rms = getproperty.(spatial, :rms_error)
    maxerr = getproperty.(spatial, :max_error)
    successive_absolute = fill(NaN, length(spatial))
    successive_l2 = fill(NaN, length(spatial))
    for i in 2:length(spatial)
        spatial[i - 1].time_star == spatial[i].time_star || continue
        change = ssw_curve_change(spatial[i - 1], spatial[i])
        successive_absolute[i] = change.max_absolute
        successive_l2[i] = change.relative_l2
    end

    fig = plt.figure("ssw_convergence", figsize=(10, 4.5))
    fig.clear()
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    ax1.loglog(panels, rms, "o-"; label="Wagner RMS")
    ax1.loglog(panels, maxerr, "s-"; label="Wagner max")
    ax2.semilogx(panels[2:end], 100 .* successive_absolute[2:end], "o-";
        label="successive max absolute")
    ax2.semilogx(panels[2:end], 100 .* successive_l2[2:end], "s-";
        label="successive relative L2")
    ax2.axhline(2.0; color="k", linestyle=":", label="2% target")
    ax1.set(xlabel="surface panels", ylabel="error", title="Wagner error")
    ax2.set(xlabel="surface panels", ylabel="percent", title="Mesh convergence")
    for ax in (ax1, ax2)
        ax.grid(true; alpha=0.3)
        ax.legend(fontsize=8, loc="best")
    end
    fig.tight_layout()
    filename = joinpath(output_root, "convergence.png")
    fig.savefig(filename; dpi=180, bbox_inches="tight")
    return filename
end

function run_ssw_convergence(config::SSWConfig=SSWConfig(); tolerance=0.02,
        span_values=(12, 24, 48, 96, 192), airfoil_values=(21, 41, 81),
        resume::Bool=true)
    cache = Dict{Tuple{Int, Int, Float64, Symbol}, Any}()
    results = Any[]
    run_case(na, ns, dtstar=config.dt_star, kind=config.backend_kind) =
            get!(cache, (na, ns, Float64(dtstar), kind)) do
        c = ssw_with(config; n_airfoil=na, n_span=ns, dt_star=dtstar,
            backend_kind=kind)
        history_path = joinpath(c.output_root, _ssw_case_tag(c), "history.csv")
        r = resume && isfile(history_path) ? load_ssw_result(c) :
            run_suddenly_started_wing(c)
        push!(results, r)
        _write_ssw_summary(results, config.output_root)
        plot_ssw_results(results; output_root=config.output_root)
        r
    end

    direct_crosscheck = nothing
    backend_change = nothing
    if config.backend_kind == :fmm
        direct_crosscheck = run_case(first(airfoil_values), first(span_values),
            config.dt_star, :direct)
        fmm_coarse = run_case(first(airfoil_values), first(span_values))
        backend_change = ssw_curve_change(direct_crosscheck, fmm_coarse)
    end

    previous = nothing
    accepted_span = last(span_values)
    for ns in span_values
        current = run_case(first(airfoil_values), ns)
        if previous !== nothing
            change = ssw_curve_change(previous, current)
            @printf("  span %d -> %d: max relative=%7.3f%%, L2=%7.3f%%\n",
                previous.config.n_span, ns, 100 * change.max_relative,
                100 * change.relative_l2)
            if change.max_absolute <= tolerance && change.relative_l2 <= tolerance
                accepted_span = ns
                break
            end
        end
        previous = current
    end

    previous = nothing
    accepted_airfoil = last(airfoil_values)
    for na in airfoil_values
        current = run_case(na, accepted_span)
        if previous !== nothing
            change = ssw_curve_change(previous, current)
            @printf("  airfoil %d -> %d: max relative=%7.3f%%, L2=%7.3f%%\n",
                previous.config.n_airfoil, na, 100 * change.max_relative,
                100 * change.relative_l2)
            if change.max_absolute <= tolerance && change.relative_l2 <= tolerance
                accepted_airfoil = na
                break
            end
        end
        previous = current
    end

    accepted = run_case(accepted_airfoil, accepted_span)
    confirm = run_case(accepted_airfoil, 2 * accepted_span)
    joint_change = ssw_curve_change(accepted, confirm)
    joint_converged = joint_change.max_absolute <= tolerance &&
        joint_change.relative_l2 <= tolerance
    accepted = confirm
    if !joint_converged
        @warn "joint span confirmation exceeded tolerance" joint_change
    end

    temporal = run_case(accepted.config.n_airfoil, accepted.config.n_span,
        accepted.config.dt_star / 2)
    temporal_ratio = temporal.ratio[2:2:end]
    length(temporal_ratio) == length(accepted.ratio) ||
        error("half-timestep history does not align with baseline")
    temporal_delta = temporal_ratio .- accepted.ratio
    temporal_change = (; max_absolute=maximum(abs, temporal_delta),
        max_relative=maximum(abs.(temporal_delta) ./
            max.(abs.(temporal_ratio), sqrt(eps(Float64)))),
        relative_l2=LA.norm(temporal_delta) / LA.norm(temporal_ratio))

    summary_path = _write_ssw_summary(results, config.output_root)
    plot_path = plot_ssw_results(results; output_root=config.output_root,
        emphasize_tag=accepted.tag)
    convergence_plot_path = plot_ssw_convergence(results;
        output_root=config.output_root, backend_kind=config.backend_kind,
        dt_star=config.dt_star)
    return (; results, accepted, confirm, temporal, joint_change,
        joint_converged, temporal_change, direct_crosscheck, backend_change,
        summary_path, plot_path, convergence_plot_path)
end

#--- wake-representation consistency study -------------------------------------
#
# A constant-strength doublet sheet has exactly cancelling interior edges, so
# "attached strip on [0,Das] + wake rows on [Das,L]" is analytically one uniform
# sheet on [0,L]. With `freestream_convection=true` the PanelWake sheet is
# straight along Uinf, i.e. *geometrically identical* to the semi-infinite wake
# truncated at L. The panel-wake CL must therefore converge to the semi-infinite
# CL as L grows. If it does not, the two wake representations disagree, which is
# a direct candidate for the +37% Das sensitivity seen in the rotor study.
# See plans/20260729_wake_representation_consistency.md.

function _write_ssw_wake_consistency_csv(rows, output_root)
    mkpath(output_root)
    filename = joinpath(output_root, "wake_consistency.csv")
    open(filename, "w") do io
        println(io, "n_rows,L_over_c,cl_panel,cl_semiinf,rel_err,cd_panel," *
            "cd_semiinf,eta,freestream_convection,AR,backend,elapsed_s,tag")
        for r in rows
            @printf(io, "%d,%.8g,%.16e,%.16e,%.16e,%.16e,%.16e,%.8g,%s,%.8g,%s,%.8g,%s\n",
                r.n_rows, r.L_over_c, r.cl_panel, r.cl_semiinf, r.rel_err,
                r.cd_panel, r.cd_semiinf, r.eta, string(r.freestream_convection),
                r.AR, string(r.backend), r.elapsed_s, r.tag)
        end
    end
    return filename
end

"""
Observed convergence order of |rel_err| against wake length, from a least-squares
fit of log|rel_err| vs log(L/c). Returns NaN with fewer than two usable points.
"""
function _ssw_convergence_order(L, err)
    keep = [i for i in eachindex(L) if isfinite(err[i]) && err[i] > 0 && L[i] > 0]
    length(keep) >= 2 || return NaN
    x = log.(L[keep])
    y = log.(abs.(err[keep]))
    xbar, ybar = sum(x) / length(x), sum(y) / length(y)
    denom = sum((x .- xbar).^2)
    denom > 0 || return NaN
    return -sum((x .- xbar) .* (y .- ybar)) / denom
end

function plot_ssw_wake_consistency(rows; output_root)
    isempty(rows) && return nothing
    mkpath(output_root)
    fig = plt.figure("ssw_wake_consistency", figsize=(10, 4.5))
    fig.clear()
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    for fc in unique(getproperty.(rows, :freestream_convection))
        sub = sort([r for r in rows if r.freestream_convection == fc];
            by=r -> r.L_over_c)
        isempty(sub) && continue
        label = fc ? "freestream-convected (straight sheet)" : "rolled-up (control)"
        marker = fc ? "o-" : "s--"
        ax1.semilogx(getproperty.(sub, :L_over_c), getproperty.(sub, :cl_panel),
            marker; label)
        ax2.loglog(getproperty.(sub, :L_over_c),
            abs.(getproperty.(sub, :rel_err)), marker; label)
    end
    semiinf = first(rows).cl_semiinf
    ax1.axhline(semiinf; color="k", linestyle=":", label="semi-infinite wake")
    ax2.axhline(0.01; color="k", linestyle=":", label="1%")
    ax1.set(xlabel="wake length L/c", ylabel="CL at final step",
        title="PanelWake vs semi-infinite wake")
    ax2.set(xlabel="wake length L/c", ylabel="|CL/CLsemiinf - 1|",
        title="Wake-representation discrepancy")
    for ax in (ax1, ax2)
        ax.grid(true; alpha=0.3, which="both")
        ax.legend(fontsize=8, loc="best")
    end
    fig.tight_layout()
    filename = joinpath(output_root, "wake_consistency.png")
    fig.savefig(filename; dpi=180, bbox_inches="tight")
    return filename
end

function run_ssw_wake_consistency(config::SSWConfig=SSWConfig();
        t_end_values=(1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0),
        convection_modes=(false, true), resume::Bool=true)

    config.backend_kind == :direct || @warn """
        FMM inflates the wake-panel radius by |Das| (FLOWPanel_liftingbody.jl),
        so its acceptance criterion is an explicit function of the variable
        under study. Use SSW_BACKEND=direct for this experiment.""" maxlog=1

    # A connectivity asymmetry here would bias Gamma spanwise without refining
    # away, so check before spending hours on the sweep.
    assert_ssw_mesh_symmetry(config.c, config.AR * config.c;
        n_span=config.n_span, n_airfoil=config.n_airfoil, caps=ssw_caps(config))

    rows = Any[]
    for fc in convection_modes, t_end in t_end_values
        c = ssw_with(config; freestream_convection=fc, t_end_star=t_end)
        history_path = joinpath(c.output_root, _ssw_case_tag(c), "history.csv")
        r = resume && isfile(history_path) ? load_ssw_result(c) :
            run_suddenly_started_wing(c)
        n_rows = length(r.time_star)
        push!(rows, (; n_rows, L_over_c=n_rows * c.dt_star,
            cl_panel=r.CL[end], cl_semiinf=r.steady_CL,
            rel_err=r.CL[end] / r.steady_CL - 1,
            cd_panel=r.CD[end], cd_semiinf=r.steady_CD,
            eta=c.eta, freestream_convection=fc, AR=c.AR,
            backend=c.backend_kind, elapsed_s=r.elapsed, tag=r.tag))
        csv_path = _write_ssw_wake_consistency_csv(rows, config.output_root)
        @printf("  L/c=%6.2f fc=%-5s CL=%+.8f  semiinf=%+.8f  rel_err=%+9.3f%%\n",
            rows[end].L_over_c, string(fc), rows[end].cl_panel,
            rows[end].cl_semiinf, 100 * rows[end].rel_err)
        config.verbose && println("  -> $csv_path")
    end

    csv_path = _write_ssw_wake_consistency_csv(rows, config.output_root)
    plot_path = plot_ssw_wake_consistency(rows; output_root=config.output_root)

    orders = Dict{Bool, Float64}()
    for fc in convection_modes
        sub = sort([r for r in rows if r.freestream_convection == fc];
            by=r -> r.L_over_c)
        orders[fc] = _ssw_convergence_order(getproperty.(sub, :L_over_c),
            abs.(getproperty.(sub, :rel_err)))
        @printf("  freestream_convection=%-5s observed order in L/c: %+.3f\n",
            string(fc), orders[fc])
    end
    return (; rows, orders, csv_path, plot_path)
end

_envbool(name, default) = lowercase(get(ENV, name, string(default))) in
    ("1", "true", "yes", "on")

function _ssw_config_from_env()
    backend = Symbol(lowercase(get(ENV, "SSW_BACKEND", "fmm")))
    return SSWConfig(
        n_airfoil=parse(Int, get(ENV, "SSW_NAIRFOIL", "21")),
        n_span=parse(Int, get(ENV, "SSW_NSPAN", "12")),
        AR=parse(Float64, get(ENV, "SSW_AR", "100.0")),
        eta=parse(Float64, get(ENV, "SSW_ETA", "0.25")),
        bodytype=Symbol(lowercase(get(ENV, "SSW_BODYTYPE", "neumann"))),
        caps=Symbol(lowercase(get(ENV, "SSW_CAPS", "auto"))),
        freestream_convection=_envbool("SSW_FREESTREAM_CONVECTION", false),
        wake_model=Symbol(lowercase(get(ENV, "SSW_WAKE_MODEL", "panel"))),
        panel_rows=parse(Int, get(ENV, "SSW_PANEL_ROWS", "8")),
        shed_method=Symbol(lowercase(get(ENV, "SSW_SHED_METHOD", "overlap_pps"))),
        sigma_over_c=parse(Float64, get(ENV, "SSW_SIGMA_OVER_C", "NaN")),
        pps_n=parse(Int, get(ENV, "SSW_PPS", "2")),
        pps_overlap=parse(Float64, get(ENV, "SSW_OVERLAP", "1.3")),
        wake_core_over_c=parse(Float64, get(ENV, "SSW_WAKE_CORE", "0.001")),
        max_particles=parse(Int, get(ENV, "SSW_MAX_PARTICLES", "200000")),
        dt_star=parse(Float64, get(ENV, "SSW_DT_STAR", "0.125")),
        t_end_star=parse(Float64, get(ENV, "SSW_T_END_STAR", "7.0")),
        output_root=get(ENV, "SSW_OUTPUT", DEFAULT_SSW_OUTPUT),
        save_vtk=_envbool("SAVE_VTK", true),
        backend_kind=backend,
        fmm_expansion_order=parse(Int, get(ENV, "SSW_FMM_ORDER", "10")),
        fmm_multipole_acceptance=parse(Float64,
            get(ENV, "SSW_FMM_THETA", "0.4")),
        fmm_leaf_size=parse(Int, get(ENV, "SSW_FMM_LEAF", "100")),
        verbose=_envbool("SSW_VERBOSE", true),
        backslash_max_panels=parse(Int,
            get(ENV, "SSW_BACKSLASH_MAX_PANELS", "10000")),
    )
end

function main()
    mode = Symbol(lowercase(get(ENV, "SSW_MODE", "single")))
    config = _ssw_config_from_env()
    if mode == :single
        result = run_suddenly_started_wing(config)
        plot_ssw_results([result]; output_root=config.output_root)
        return result
    elseif mode == :coarse
        coarse = ssw_with(config; n_airfoil=21, n_span=12)
        result = run_suddenly_started_wing(coarse)
        plot_ssw_results([result]; output_root=config.output_root)
        return result
    elseif mode == :convergence
        return run_ssw_convergence(config; resume=_envbool("SSW_RESUME", true))
    elseif mode == :wake_consistency
        t_end_values = Tuple(parse.(Float64,
            split(get(ENV, "SSW_WAKE_LENGTHS", "1,2,4,8,16,32,64"), ',')))
        convection_modes = Tuple(lowercase(strip(s)) in ("1", "true", "yes", "on")
            for s in split(get(ENV, "SSW_WAKE_MODES", "false,true"), ','))
        # AR=100 is unusable here: its steady semi-infinite CL is O(1e5) (see
        # the header note), so cl_panel/cl_semiinf would be a ratio of garbage.
        # Default this mode to AR=6, where the steady solve gives cl≈0.39.
        AR = haskey(ENV, "SSW_AR") ? config.AR : 6.0
        # 512-step wakes make VTK output dominate runtime and disk.
        c = ssw_with(config; AR, save_vtk=_envbool("SAVE_VTK", false))
        return run_ssw_wake_consistency(c; t_end_values, convection_modes,
            resume=_envbool("SSW_RESUME", true))
    end
    error("SSW_MODE must be single, coarse, convergence, or wake_consistency; " *
        "got $mode")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
