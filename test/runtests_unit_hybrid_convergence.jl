using Test
import FLOWPanel as pnl
import FastMultipole
import LinearAlgebra as LA
using StaticArrays: SVector

include("test_helpers.jl")

const DIRECT_HYBRID_ORACLE = pnl.DirectBackend()
const UINF_HYBRID_ORACLE = SVector(1.0, 0.07, 0.12)

# Stretched, warped, non-uniform sheet and non-affine strength field reused
# from runtests_unit_wake.jl's static sheet-to-hybrid acceptance fixture.
_hybrid_oracle_nodes(irow, icol) =
    [0.4 * (1.35^(irow - 1) - 1) / 0.35,
     0.8 * (icol - 1) + 0.15 * (icol - 1)^2,
     0.05 * (irow - 1) * (icol - 1)]
_hybrid_oracle_mu(irow, icol) =
    0.1irow + 0.3icol + 0.07 * irow * icol - 0.04 * icol^2

function hybrid_oracle_pair(; nwakerows=3, sigma=0.12, h=0.06,
        attribution=:upstream)
    overlap = sigma / h
    conversion = pnl.SurfaceVorticityConversion(sigma; overlap, attribution)
    hybrid = make_conversion_fixture(; nwakerows, wrap=false,
        max_particles=200_000, conversion,
        node_fun=_hybrid_oracle_nodes, strength_fun=_hybrid_oracle_mu)
    hpw = hybrid.panel_wake

    # The velocity oracle must retain the final filament: it is part of both
    # the full panel wake and the converted representation.  The formulation
    # comparison below makes a scalar-only copy without that filament because
    # DirectWakePotential deliberately rejects vector-only wake sources.
    body = make_dirichlet_diamond_body(; nspan=3)
    pure = pnl.PanelWake(body; nwakerows, include_final_filament=true)
    for (dst, src) in zip(pure.nodes, hpw.nodes)
        dst .= src
    end
    for (dst, src) in zip(pure.strength, hpw.strength)
        dst .= src
    end
    pure.nwakes[] = nwakerows

    pnl._convert_to_particles!(hybrid, body)
    hpw.nwakes[] = nwakerows - 1
    return body, pure, hybrid
end

function hybrid_oracle_probe_positions(pure, distance)
    positions = SVector{3,Float64}[]
    nodes = pure.nodes[1]
    row = pure.nwakes[]
    for col in axes(pure.strength[1], 3)
        v1 = SVector{3}(nodes[:, row, col])
        v2 = SVector{3}(nodes[:, row + 1, col])
        v3 = SVector{3}(nodes[:, row + 1, col + 1])
        v4 = SVector{3}(nodes[:, row, col + 1])
        center = (v1 + v2 + v3 + v4) / 4
        normal = LA.cross(v2 - v1, v4 - v1)
        normal /= LA.norm(normal)
        for side in (-1, 1)
            push!(positions, center + side * distance * normal)
        end
    end
    return positions
end

function hybrid_oracle_velocity(sources, positions)
    probes = FastMultipole.ProbeSystem(length(positions), Float64)
    for i in eachindex(positions)
        probes.position[i] = positions[i]
        probes.scalar_potential[i] = 0.0
        probes.gradient[i] = zero(SVector{3,Float64})
        probes.hessian[i] = zero(FastMultipole.SMatrix{3,3,Float64,9})
    end
    pnl.influence!((probes,), sources, DIRECT_HYBRID_ORACLE;
        scalar_potential=false, gradient=true, hessian=(false,))
    return copy(probes.gradient)
end

function hybrid_oracle_field_metrics(; nwakerows=3, sigma=0.12, h=0.06,
        distance=0.24, attribution=:upstream)
    _, pure, hybrid = hybrid_oracle_pair(;
        nwakerows, sigma, h, attribution)
    positions = hybrid_oracle_probe_positions(pure, distance)
    reference = hybrid_oracle_velocity(pnl.get_sources(pure), positions)
    converted = hybrid_oracle_velocity(
        (pnl.get_sources(hybrid.panel_wake)..., hybrid.pfield), positions)
    delta = [LA.norm(converted[i] - reference[i]) for i in eachindex(reference)]
    reference_norm = [LA.norm(u) for u in reference]
    converted_norm = [LA.norm(u) for u in converted]
    rms = sqrt(sum(abs2, delta) / length(delta)) /
        max(sqrt(sum(abs2, reference_norm) / length(reference_norm)), eps())
    maximum_relative = maximum(delta) / max(maximum(reference_norm), eps())
    reference_rms = sqrt(sum(abs2, reference_norm) / length(reference_norm))
    converted_rms = sqrt(sum(abs2, converted_norm) / length(converted_norm))
    correlation = sum(LA.dot(reference[i], converted[i])
        for i in eachindex(reference)) /
        max(sqrt(sum(abs2, reference_norm) * sum(abs2, converted_norm)), eps())
    return (; rms, maximum_relative, reference_rms, converted_rms,
        correlation, np=hybrid.pfield.np)
end

function run_hybrid_oracle_aero!(body, wake, solver, formulation, state)
    frames = pnl.ReferenceFrame(body)
    pnl._steady_aerodynamics!(body, (body,), (wake,), frames,
        UINF_HYBRID_ORACLE, solver;
        backend_wake=DIRECT_HYBRID_ORACLE,
        backend_solve=DIRECT_HYBRID_ORACLE,
        backend_system=DIRECT_HYBRID_ORACLE,
        update_trailing_edges=false,
        formulation, formulation_state=state, i_step=1)
    return body
end

function hybrid_oracle_formulation_metrics(; sigma=0.12, h=0.06,
        attribution=:upstream)
    body_template, pure, hybrid = hybrid_oracle_pair(;
        nwakerows=3, sigma, h, attribution)

    pure_scalar = pnl.PanelWake(body_template; nwakerows=pure.nwakes[],
        include_final_filament=false)
    for (dst, src) in zip(pure_scalar.nodes, pure.nodes)
        dst .= src
    end
    for (dst, src) in zip(pure_scalar.strength, pure.strength)
        dst .= src
    end
    pure_scalar.nwakes[] = pure.nwakes[]

    body_direct = deepcopy(body_template)
    solver_direct = pnl.Backslash(body_direct)
    direct = pnl.DirectWakePotential()
    direct_state = pnl.initialize_formulation(direct, (body_direct,), (pure_scalar,),
        solver_direct, DIRECT_HYBRID_ORACLE, DIRECT_HYBRID_ORACLE)
    run_hybrid_oracle_aero!(body_direct, pure_scalar, solver_direct, direct,
        direct_state)

    body_hybrid = deepcopy(body_template)
    solver_hybrid = pnl.Backslash(body_hybrid)
    formulation = pnl.HybridWakePotential(outer_tolerance=1e-11,
        require_outer_convergence=true)
    hybrid_state = pnl.initialize_formulation(formulation, (body_hybrid,),
        (hybrid,), solver_hybrid, DIRECT_HYBRID_ORACLE, DIRECT_HYBRID_ORACLE)
    run_hybrid_oracle_aero!(body_hybrid, hybrid, solver_hybrid, formulation,
        hybrid_state)

    bs = hybrid_state.bodies[1]
    areas = pnl.calc_areas(body_hybrid)
    dq = bs.q_total .- direct_state.q_wake
    gauge_offset = LA.dot(areas, dq) / sum(areas)
    dq .-= gauge_offset
    q_reference = direct_state.q_wake .-
        LA.dot(areas, direct_state.q_wake) / sum(areas)
    q_error = LA.norm(dq) / max(LA.norm(q_reference), eps())
    source_error = LA.norm(view(body_hybrid.strength, :, 1) .-
        view(body_direct.strength, :, 1)) /
        max(LA.norm(view(body_direct.strength, :, 1)), eps())
    doublet_error = LA.norm(view(body_hybrid.strength, :, 2) .-
        view(body_direct.strength, :, 2)) /
        max(LA.norm(view(body_direct.strength, :, 2)), eps())
    status = pnl.block_gs_status((solver_hybrid,))
    return (; q_error, source_error, doublet_error,
        green_residual=bs.green_residual[], gauge_defect=bs.gauge_defect[],
        normalized_residual=status.normalized_residual,
        converged=status.converged)
end

@testset verbose=true "Panel-wake-to-particle convergence oracle (052b B)" begin
    # Collect and emit the independent sweeps before applying the calibrated
    # gates.  Calibration baseline (2026-08-27, Float64/DirectBackend): the
    # finest h case gave RMS/max 0.6422/0.6058 and correlation 0.869; all nine
    # handoff cases stayed below 0.783 RMS; formulation q and doublet errors
    # stayed below 0.104 and 0.091.  The limits retain roughly 2--20% headroom
    # while still detecting a lost filament, sign reversal, stalled spatial
    # refinement, or broken Green solve.  Near-core d=0.25 is diagnostic only:
    # a regularized particle sheet is not expected to reproduce the singular
    # panel kernel there.
    h_sweep = [(h, hybrid_oracle_field_metrics(; sigma=0.12, h,
        distance=1.0)) for h in (0.24, 0.12, 0.06)]
    sigma_sweep = [(sigma, hybrid_oracle_field_metrics(; sigma, h=0.03,
        distance=1.0)) for sigma in (0.24, 0.12, 0.06)]
    distance_sweep = [(distance, hybrid_oracle_field_metrics(;
        sigma=0.12, h=0.06, distance)) for distance in (0.25, 0.5, 1.0)]
    handoff_sweep = [((nwakerows, attribution),
        hybrid_oracle_field_metrics(; nwakerows, sigma=0.12, h=0.06,
            distance=1.0, attribution))
        for nwakerows in (2, 3, 4),
            attribution in (:upstream, :split, :downstream)]

    println("\n052b hybrid convergence diagnostics")
    println("  h sweep (fixed sigma=0.12, d=1.0): ", h_sweep)
    println("  sigma sweep (fixed h=0.03, d=1.0): ", sigma_sweep)
    println("  distance sweep (fixed sigma=0.12, h=0.06): ", distance_sweep)
    println("  handoff row/attribution sweep: ", vec(handoff_sweep))

    all_field_metrics = [last(x) for x in vcat(h_sweep, sigma_sweep,
        distance_sweep, vec(handoff_sweep))]
    @test all(m -> all(isfinite, (m.rms, m.maximum_relative,
        m.reference_rms, m.converted_rms, m.correlation)) && m.np > 0,
        all_field_metrics)

    # Fixed-core spatial quadrature converges monotonically as h decreases.
    h_metrics = last.(h_sweep)
    @test h_metrics[3].rms < h_metrics[2].rms < h_metrics[1].rms
    @test h_metrics[3].maximum_relative < h_metrics[2].maximum_relative <
        h_metrics[1].maximum_relative
    @test h_metrics[3].rms <= 0.66
    @test h_metrics[3].maximum_relative <= 0.62
    @test h_metrics[3].correlation >= 0.86

    # At a standoff more than eight cores away, changing sigma at fixed h is
    # correctly negligible: this isolates core bias from particle spacing.
    sigma_metrics = last.(sigma_sweep)
    @test maximum(m.rms for m in sigma_metrics) -
        minimum(m.rms for m in sigma_metrics) <= 5e-5
    @test all(m -> m.rms <= 0.66 && m.correlation >= 0.86, sigma_metrics)

    # The converted representation approaches the direct sheet as physical
    # standoff grows; the near-core endpoint remains visible in diagnostics.
    distance_metrics = last.(distance_sweep)
    @test distance_metrics[3].rms < distance_metrics[2].rms <
        distance_metrics[1].rms
    @test distance_metrics[3].maximum_relative <
        distance_metrics[2].maximum_relative <
        distance_metrics[1].maximum_relative
    @test distance_metrics[3].rms <= 0.66

    # Moving the handoff and exercising all attribution modes must remain
    # bounded.  This is not an attribution selector: the separate static
    # near-field decision rule in runtests_unit_wake.jl retains :upstream.
    handoff_metrics = last.(vec(handoff_sweep))
    @test all(m -> m.rms <= 0.80 && m.maximum_relative <= 0.75,
        handoff_metrics)

    formulation_metrics = [(attribution,
        hybrid_oracle_formulation_metrics(; attribution))
        for attribution in (:upstream, :split, :downstream)]
    println("  formulation Direct-vs-Hybrid: ", formulation_metrics)
    formulation_values = last.(formulation_metrics)
    @test all(m -> all(isfinite, (m.q_error, m.source_error,
        m.doublet_error, m.green_residual, m.gauge_defect,
        m.normalized_residual)) && m.converged, formulation_values)
    @test all(m -> m.q_error <= 0.12, formulation_values)
    @test all(m -> m.source_error <= 1e-13, formulation_values)
    @test all(m -> m.doublet_error <= 0.10, formulation_values)
    @test all(m -> m.green_residual <= 1e-12, formulation_values)
    @test all(m -> m.gauge_defect <= 1e-12, formulation_values)
    @test all(m -> m.normalized_residual <= 1e-12, formulation_values)
end
