using Test
import FLOWPanel as pnl
import FLOWVPM
const TOML_FW = pnl.TOML
import FastMultipole
using LinearAlgebra: norm
using StaticArrays: SVector
using Random

if !isdefined(@__MODULE__, :make_dirichlet_diamond_body)
    include("test_helpers.jl")
end

#=
A6 unit tests for TrailingFilamentSheet / FilamentParticleWake
(BRAINSTORM 018 r4). Numbering follows the execution plan:
 (1) index-map round-trips, (2) ring-decomposition identity (LOAD-BEARING,
 pins the kernel sign convention), (3) steady-DeltaGamma invariance,
 (4) FMM vs Direct with per-station cores, (5) handoff exactness,
 (6) guards, (7) replay manifest + warmstart round-trips.
=#

# Deterministic TrailingFilamentSheet with hand-written rows (sibling of
# make_conversion_fixture, which builds the PanelParticleWake analogue).
function make_filament_fixture(; nwakerows::Int, wrap::Bool, max_particles=2000,
        method_trailing=pnl.SigmaOverlap(0.2, 2.75), strength_fun=nothing,
        node_fun=nothing, core_size=0.05, optargs...)
    body = make_dirichlet_diamond_body(; nspan=3)
    wake = pnl.FilamentParticleWake(body; nwakerows=nwakerows,
        max_particles=max_particles, method_trailing=method_trailing,
        core_size=core_size, optargs...)
    sheet = wake.sheet
    nodes = sheet.nodes[1]
    strength = sheet.strength[1]
    n_node_rows = size(nodes, 2)
    n_node_cols = size(nodes, 3)

    for irow in 1:n_node_rows, icol in 1:n_node_cols
        if node_fun !== nothing
            nodes[:, irow, icol] .= node_fun(irow, icol)
        elseif wrap
            theta = 2pi * (icol - 1) / (n_node_cols - 1)
            nodes[1, irow, icol] = cos(theta)
            nodes[2, irow, icol] = sin(theta)
            nodes[3, irow, icol] = 0.25 * (irow - 1)
        else
            nodes[1, irow, icol] = 0.5 * (irow - 1)
            nodes[2, irow, icol] = icol - 1
            nodes[3, irow, icol] = 0.0
        end
    end
    wrap && (nodes[:, :, n_node_cols] .= nodes[:, :, 1])

    f = strength_fun === nothing ? (irow, icol) -> 0.1 * irow + 0.3 * icol : strength_fun
    for irow in 1:size(strength, 2), icol in 1:size(strength, 3)
        strength[1, irow, icol] = f(irow, icol)
    end

    sheet.nwakes[] = nwakerows
    # ctor wraps are degenerate-true on zero nodes; re-detect from the real rows
    pnl._refresh_wraps!(sheet, nothing)
    return wake
end

function make_probes(targets)
    probes = FastMultipole.ProbeSystem(length(targets), Float64)
    for i in eachindex(targets)
        probes.position[i] = SVector{3,Float64}(targets[i])
        probes.scalar_potential[i] = 0.0
        probes.gradient[i] = zero(SVector{3,Float64})
        probes.hessian[i] = zero(FastMultipole.SMatrix{3,3,Float64,9})
    end
    return probes
end

function probe_field(sources, targets; backend=pnl.DirectBackend(), hessian=true)
    probes = make_probes(targets)
    pnl.influence!((probes,), sources, backend;
        scalar_potential=false, gradient=true, hessian=(hessian,))
    return copy(probes.gradient), copy(probes.hessian)
end

relrms_fw(a, b) = norm(vec(a) .- vec(b)) / max(norm(vec(b)), eps())

# Spanwise completion of the trailing decomposition, built as a TRANSPOSED
# TrailingFilamentSheet: its "trailing" segments run along the original
# spanwise direction. With stored strengths Gamma'(c, rho) = -Gamma(rho, c),
# its station jumps telescope to exactly the spanwise segment strengths of the
# ring decomposition: node row rho carries Gamma(rho-1,c) - Gamma(rho,c) in
# the +col direction, with Gamma(0,.) = 0 and the bottom boundary
# rho = nwakes+1 carrying +Gamma(nwakes,c).
function make_spanwise_closure(nodes, strengthP, nwakes, ncols, core_size)
    tsheet = pnl.TrailingFilamentSheet([zeros(Int, 3, nwakes)], Float64;
        nwakerows=ncols, core_size=core_size)
    tnodes = tsheet.nodes[1]
    tstrength = tsheet.strength[1]
    for a in 1:(ncols+1), b in 1:(nwakes+1)
        tnodes[:, a, b] .= view(nodes, :, b, a)
    end
    for a in 1:ncols, rho in 1:nwakes
        tstrength[1, a, rho] = -strengthP[1, rho, a]
    end
    tsheet.nwakes[] = ncols
    pnl._refresh_wraps!(tsheet, nothing)
    @assert !tsheet.wraps[1] "spanwise-closure sheet must be an open chain"
    return tsheet
end

@testset verbose=true "FilamentParticleWake (018 r4)" begin

    #--- (1) index-map round-trips ---#
    @testset "index map round-trips (open, wrapped, two-surface, partial fill)" begin
        nwakerows = 3
        shedding = [zeros(Int, 3, 4), zeros(Int, 3, 2)]
        sigmas = [collect(0.01 .* (1:5)), collect(0.1 .* (1:3))]
        sheet = pnl.TrailingFilamentSheet(shedding, Float64;
            nwakerows=nwakerows, core_size=0.03, filament_core_size=sigmas)

        # surface 1 open, surface 2 wrapped (last node column == first)
        for (isurf, nodes) in enumerate(sheet.nodes)
            for irow in axes(nodes, 2), icol in axes(nodes, 3)
                nodes[1, irow, icol] = 100.0 * isurf + irow
                nodes[2, irow, icol] = Float64(icol)
                nodes[3, irow, icol] = 0.1 * irow + 0.01 * icol
            end
        end
        sheet.nodes[2][:, :, end] .= sheet.nodes[2][:, :, 1]
        for (isurf, strength) in enumerate(sheet.strength)
            for irow in axes(strength, 2), icol in axes(strength, 3)
                strength[1, irow, icol] = 10.0 * isurf + irow + 0.1 * icol
            end
        end

        for nwakes in (2, 3) # partial fill, then full
            sheet.nwakes[] = nwakes
            pnl._refresh_wraps!(sheet, nothing)
            @test sheet.wraps == [false, true]
            n_st = (pnl._n_stations(sheet, 1), pnl._n_stations(sheet, 2))
            @test n_st == (5, 2)
            nb = FastMultipole.get_n_bodies(sheet)
            @test nb == nwakes * sum(n_st)

            seen = Set{Tuple{Int,Int,Int}}()
            for i in 1:nb
                isurf, irow, j = pnl.global_to_matrix_index(sheet, i)
                @test 1 <= isurf <= 2
                @test 1 <= irow <= nwakes
                @test 1 <= j <= n_st[isurf]
                @test pnl.matrix_to_global_index(sheet, isurf, irow, j) == i
                push!(seen, (isurf, irow, j))

                nodes = sheet.nodes[isurf]
                mid = 0.5 .* (nodes[:, irow, j] .+ nodes[:, irow+1, j])
                @test FastMultipole.get_position(sheet, i) ≈ SVector{3}(mid...)

                buffer = zeros(12, 1)
                FastMultipole.source_system_to_buffer!(buffer, 1, sheet, i)
                @test buffer[1:3, 1] ≈ mid
                @test buffer[6:8, 1] == nodes[:, irow, j]
                @test buffer[9:11, 1] == nodes[:, irow+1, j]
                @test buffer[12, 1] == sigmas[isurf][j]
                @test buffer[5, 1] ==
                    pnl._trailing_filament_strength(sheet, isurf, irow, j)
            end
            @test length(seen) == nb # bijective
        end

        # wrapped station-1 jump closes the chain
        s2 = sheet.strength[2]
        @test pnl._trailing_filament_strength(sheet, 2, 1, 1) ==
            s2[1, 1, 1] - s2[1, 1, 2]
        # open boundary stations
        s1 = sheet.strength[1]
        @test pnl._trailing_filament_strength(sheet, 1, 2, 1) == s1[1, 2, 1]
        @test pnl._trailing_filament_strength(sheet, 1, 2, 5) == -s1[1, 2, 4]

        # probe map (wake nodes as convection probes)
        sheet.nwakes[] = 2
        pw = pnl.ProbeWrapper(sheet)
        npr = FastMultipole.get_n_bodies(pw)
        @test npr == (2 + 1) * (5 + 3)
        for i in 1:npr
            isurf, irow, icol = pnl.global_to_matrix_index(pw, i)
            @test pnl.matrix_to_global_index(pw, isurf, irow, icol) == i
            @test FastMultipole.get_position(pw, i) ≈
                SVector{3}(sheet.nodes[isurf][:, irow, icol]...)
        end
    end

    #--- (2) ring-decomposition identity (pins the kernel sign convention) ---#
    @testset "ring decomposition identity: PanelWake rings == trailing + spanwise" begin
        # DO NOT "fix" the _direct_filaments! argument asymmetry (target - v
        # for velocity, v - target for gradient/vorticity): this test pins it.
        Random.seed!(1234)
        sigma0 = 0.08
        nwakes, ncols = 3, 5

        for wrap in (false, true)
            ring = pnl.PanelWake([zeros(Int, 3, ncols)], pnl.VortexRing, Float64;
                core_size=sigma0, nwakerows=nwakes,
                include_final_filament=false, unsteady_filament=false)
            nodes = ring.nodes[1]
            strengthP = ring.strength[1]
            for irow in axes(nodes, 2), icol in axes(nodes, 3)
                if wrap
                    theta = 2pi * (icol - 1) / (size(nodes, 3) - 1)
                    nodes[1, irow, icol] = cos(theta) + 0.05 * randn()
                    nodes[2, irow, icol] = sin(theta) + 0.05 * randn()
                    nodes[3, irow, icol] = 0.4 * (irow - 1) + 0.05 * randn()
                else
                    nodes[1, irow, icol] = 0.6 * (irow - 1) + 0.05 * randn()
                    nodes[2, irow, icol] = 1.0 * (icol - 1) + 0.05 * randn()
                    nodes[3, irow, icol] = 0.05 * randn()
                end
            end
            wrap && (nodes[:, :, end] .= nodes[:, :, 1])
            for irow in axes(strengthP, 2), icol in axes(strengthP, 3)
                strengthP[1, irow, icol] = randn()
            end
            ring.nwakes[] = nwakes
            @test pnl.get_sources(ring) == (ring,)

            tsheet = pnl.TrailingFilamentSheet([zeros(Int, 3, ncols)], Float64;
                nwakerows=nwakes, core_size=sigma0)
            tsheet.nodes[1] .= nodes
            tsheet.strength[1] .= strengthP
            tsheet.nwakes[] = nwakes
            pnl._refresh_wraps!(tsheet, nothing)
            @test tsheet.wraps[1] == wrap

            ssheet = make_spanwise_closure(nodes, strengthP, nwakes, ncols, sigma0)

            # segment bookkeeping: rings decompose into exactly these segments
            @test FastMultipole.get_n_bodies(tsheet) ==
                nwakes * (wrap ? ncols : ncols + 1)
            @test FastMultipole.get_n_bodies(ssheet) == ncols * (nwakes + 1)

            # probes well off the sheet (identity is analytic; stay away from
            # the segments so cancellation noise cannot inflate the residual)
            targets = [SVector(4.0 * rand() - 1.0, 5.0 * rand() - 1.5,
                               (0.6 + 2.0 * rand()) * (rand() < 0.5 ? -1 : 1) .+
                               (wrap ? 0.5 : 0.0)) for _ in 1:30]

            Uring, Hring = probe_field((ring,), targets)
            Ufil, Hfil = probe_field((tsheet, ssheet), targets)

            @test isapprox(Ufil, Uring; rtol=1e-12)
            @test isapprox(Hfil, Hring; rtol=1e-12)

            # the identity must BREAK if the trailing half is removed
            Uspan, _ = probe_field((ssheet,), targets)
            @test relrms_fw(Uspan, Uring) > 1e-3
        end
    end

    #--- (3) steady loading: conversion is step-invariant ---#
    @testset "steady DeltaGamma = 0: identical particle sets, invariant sheet field" begin
        body = make_dirichlet_diamond_body(; nspan=3)
        body.strength .= reshape(0.01 .* (1:length(body.strength)),
            size(body.strength))
        wake = pnl.FilamentParticleWake(body; nwakerows=2, max_particles=500,
            core_size=0.05, method_trailing=pnl.SigmaOverlap(0.2, 2.75),
            freestream_convection=true)
        sheet = wake.sheet
        pnl.apply_freestream!(sheet, SVector(1.0, 0.0, 0.0))
        dt = 0.1
        targets = (SVector(1.8, 0.5, 0.7), SVector(0.9, -0.6, -0.5),
                   SVector(2.5, 1.4, 0.2))

        np_at = Int[]
        fields = Vector{Vector{SVector{3,Float64}}}()
        for step in 1:8
            pnl.update_TE!(wake, body)
            # probe at the solve-time state (post-pin, pre-shed): the
            # post-shift state transiently holds coincident rows 1-2 that a
            # solve never sees
            U, _ = probe_field((sheet,), collect(targets); hessian=false)
            push!(fields, U)
            pnl.shed_wake!(wake, body)
            pnl.propagate!(sheet, dt)
            push!(np_at, wake.pfield.np)
        end

        # conversions start once the 2-row buffer fills (step 3) and shed a
        # fixed, nonzero number of particles every step after
        born = diff(np_at)
        @test np_at[1] == np_at[2] == 0
        @test born[2] > 0
        @test all(born[2:end] .== born[2])

        # consecutive conversions deposit byte-identical particle sets
        nper = born[2]
        pf = wake.pfield
        for k in 2:length(born)-1
            a = (np_at[k]+1):np_at[k+1]
            b = (np_at[k+1]+1):np_at[k+2]
            @test pf.particles[1:3, a] == pf.particles[1:3, b]
            @test pf.particles[FLOWVPM.GAMMA_INDEX, a] ==
                  pf.particles[FLOWVPM.GAMMA_INDEX, b]
            @test pf.particles[FLOWVPM.SIGMA_INDEX, a] ==
                  pf.particles[FLOWVPM.SIGMA_INDEX, b]
        end

        # sheet-induced velocity at fixed probes is step-invariant once full
        for step in 4:8
            @test isapprox(fields[step], fields[3]; rtol=1e-13, atol=1e-15)
        end
        @test any(norm.(fields[3]) .> 0)
    end

    #--- (4) FMM vs Direct with per-station cores ---#
    @testset "FMM agrees with Direct; per-station sigma enters the kernel" begin
        Random.seed!(99)
        nwakes, ncols = 8, 16
        sigmas = [0.02 * 10.0^((j - 1) / ncols) for j in 1:(ncols+1)] # 10x span
        sheet = pnl.TrailingFilamentSheet([zeros(Int, 3, ncols)], Float64;
            nwakerows=nwakes, core_size=0.05, filament_core_size=[sigmas])
        nodes = sheet.nodes[1]
        for irow in axes(nodes, 2), icol in axes(nodes, 3)
            nodes[1, irow, icol] = 0.35 * (irow - 1) + 0.03 * randn()
            nodes[2, irow, icol] = 0.25 * (icol - 1) + 0.03 * randn()
            nodes[3, irow, icol] = 0.05 * randn()
        end
        strength = sheet.strength[1]
        for irow in axes(strength, 2), icol in axes(strength, 3)
            strength[1, irow, icol] = randn()
        end
        sheet.nwakes[] = nwakes
        pnl._refresh_wraps!(sheet, nothing)
        @test !sheet.wraps[1]

        targets = [SVector(6.0 * rand() - 1.5, 7.0 * rand() - 1.5,
                           (0.5 + 3.0 * rand()) * (rand() < 0.5 ? -1 : 1))
                   for _ in 1:40]

        Udir, Hdir = probe_field((sheet,), targets)
        Ufmm, Hfmm = probe_field((sheet,), targets;
            backend=pnl.FastMultipoleBackend(; expansion_order=12,
                multipole_acceptance=0.4, leaf_size=5))
        @test relrms_fw(Ufmm, Udir) < 1e-5
        @test relrms_fw(Hfmm, Hdir) < 1e-4

        # canary: collapsing the sigma ladder to its mean must change the
        # direct field — proves slot 12 actually reaches the kernel
        uniform = pnl.TrailingFilamentSheet([zeros(Int, 3, ncols)], Float64;
            nwakerows=nwakes, core_size=sum(sigmas) / length(sigmas))
        uniform.nodes[1] .= nodes
        uniform.strength[1] .= strength
        uniform.nwakes[] = nwakes
        pnl._refresh_wraps!(uniform, nothing)
        Uuni, _ = probe_field((uniform,), targets)
        @test relrms_fw(Uuni, Udir) > 1e-6
    end

    #--- (5) handoff exactness ---#
    @testset "handoff: buffer strengths == shed particles (incl. zero elision)" begin
        for wrap in (false, true)
            nwakerows = 2
            overlap = 2.75
            sigmas = [0.11, 0.17, 0.23, 0.31] # per node column (ncols + 1 = 4)
            method = pnl.StationSigmaOverlap([sigmas], overlap)
            # station 2 jump exactly zero on every row (s(1) == s(2))
            sfun = (irow, icol) -> 0.1 * irow + (icol <= 2 ? 0.3 : 0.3 * icol)
            wake = make_filament_fixture(; nwakerows, wrap,
                method_trailing=method, strength_fun=sfun)
            sheet = wake.sheet
            nodes = sheet.nodes[1]
            n_st = pnl._n_stations(sheet, 1)
            @test n_st == (wrap ? 3 : 4)
            # sigma-coupling rule: hold-sigma == shed-sigma by construction
            @test sheet.filament_core_size[1] == sigmas

            # FMM buffer strengths on the outgoing row
            nb = FastMultipole.get_n_bodies(sheet)
            buffer = zeros(12, nb)
            for i in 1:nb
                FastMultipole.source_system_to_buffer!(buffer, i, sheet, i)
            end
            outrow = sheet.nwakes[]
            expected = Dict{Int,Float64}()
            for j in 1:n_st
                i_body = pnl.matrix_to_global_index(sheet, 1, outrow, j)
                G = pnl._trailing_filament_strength(sheet, 1, outrow, j)
                @test buffer[5, i_body] === G
                expected[j] = G
            end
            @test expected[2] == 0.0 # the crafted zero jump
            @test count(!iszero, values(expected)) == n_st - 1

            pnl._convert_to_particles!(wake)
            pf = wake.pfield

            # reconstruct the exact particle sequence _shed_particles! deposits
            exp_X = Vector{Float64}[]
            exp_G = Vector{Float64}[]
            exp_sig = Float64[]
            for j in 1:n_st
                G = expected[j]
                G == 0 && continue # elided: a zero-strength filament induces nothing
                r1 = SVector{3}(nodes[:, outrow, j]...)
                r2 = SVector{3}(nodes[:, outrow+1, j]...)
                dist = norm(r2 - r1)
                p = max(1, ceil(Int, overlap * dist / sigmas[j]))
                dl = (r2 - r1) / p
                X = r1 + 0.5 * dl
                for _ in 1:p
                    push!(exp_X, collect(X))
                    push!(exp_G, collect(G * dl))
                    push!(exp_sig, sigmas[j])
                    X += dl
                end
            end
            @test pf.np == length(exp_X)
            for i in 1:pf.np
                @test pf.particles[1:3, i] == exp_X[i]
                @test pf.particles[FLOWVPM.GAMMA_INDEX, i] == exp_G[i]
                @test pf.particles[FLOWVPM.SIGMA_INDEX, i][1] == exp_sig[i]
            end

            # accounting: sum of particle Gamma per station == DeltaGamma * segment
            i0 = 0
            for j in 1:n_st
                G = expected[j]
                G == 0 && continue
                r1 = SVector{3}(nodes[:, outrow, j]...)
                r2 = SVector{3}(nodes[:, outrow+1, j]...)
                p = max(1, ceil(Int, overlap * norm(r2 - r1) / sigmas[j]))
                tot = sum(pf.particles[FLOWVPM.GAMMA_INDEX, i0+1:i0+p]; dims=2)
                @test vec(tot) ≈ G .* (r2 - r1) rtol = 1e-14
                i0 += p
            end
        end
    end

    #--- (6) guards ---#
    @testset "guards" begin
        body = make_dirichlet_diamond_body(; nspan=3)

        @test_throws ArgumentError pnl.FilamentParticleWake(body)
        @test_throws ArgumentError pnl.FilamentParticleWake(body; nwakerows=0)
        @test_throws ArgumentError pnl.FilamentParticleWake(body; nwakerows=2,
            conversion=pnl.LegacyEdgeJumpConversion())
        @test_throws ArgumentError pnl.FilamentParticleWake(body; nwakerows=2,
            method_unsteady=pnl.NoShed())

        # filament_core_size shape/positivity
        @test_throws ArgumentError pnl.TrailingFilamentSheet(body; nwakerows=2,
            filament_core_size=[[0.1, 0.2]]) # needs ncols + 1 = 4
        @test_throws ArgumentError pnl.TrailingFilamentSheet(body; nwakerows=2,
            filament_core_size=[[0.1, -0.2, 0.3, 0.4]])
        @test_throws ArgumentError pnl.TrailingFilamentSheet(body.shedding, Float64;
            nwakerows=2, filament_core_size=Vector{Float64}[])
        @test_throws ArgumentError pnl.TrailingFilamentSheet(body.shedding, Float64)

        # sigma-coupling: mismatch rejected, exact match accepted
        sig = [0.1, 0.2, 0.3, 0.4]
        method = pnl.StationSigmaOverlap([sig], 2.0)
        @test_throws ArgumentError pnl.FilamentParticleWake(body; nwakerows=2,
            method_trailing=method, filament_core_size=[sig .* 2])
        wok = pnl.FilamentParticleWake(body; nwakerows=2,
            method_trailing=method, filament_core_size=[copy(sig)])
        @test wok.sheet.filament_core_size[1] == sig
        # OmitStations-wrapped StationSigmaOverlap couples too
        wom = pnl.FilamentParticleWake(body; nwakerows=2,
            method_trailing=pnl.OmitStations(method, [[true, false, false, false]]))
        @test wom.sheet.filament_core_size[1] == sig
        # non-station-based method: scalar fallback fill
        wsc = pnl.FilamentParticleWake(body; nwakerows=2, core_size=0.07)
        @test wsc.sheet.filament_core_size[1] == fill(0.07, 4)

        # a FilamentParticleWake never contributes scalar-potential sources
        @test pnl._scalar_potential_sources(wsc) == ()
        @test FastMultipole.has_vector_potential(wsc.sheet)

        # live-row block (Kutta Route B) is rejected at conversion
        wlr = make_filament_fixture(; nwakerows=2, wrap=false)
        wlr.sheet.live_rows[] = 1
        @test_throws ErrorException pnl._convert_to_particles!(wlr)

        # mid-run topology change: cached wraps vs geometric test fails loudly
        wtp = make_filament_fixture(; nwakerows=2, wrap=false)
        wtp.sheet.wraps[1] = true
        @test_throws ErrorException pnl._convert_to_particles!(wtp)

        # Kutta Route B configuration names the new types in its rejection
        wk = pnl.FilamentParticleWake(body; nwakerows=2, max_particles=100)
        @test_throws "FilamentParticleWake" pnl._validate_kutta_configuration(
            :simulate, (body,), (wk,), (pnl.Backslash(body),),
            pnl.VelocityThroughSources(), pnl.DirectBackend(),
            pnl.TEAnchoredAttachment(), pnl.JumpKutta())
    end

    #--- (7) replay manifest + warmstart round-trips ---#
    @testset "replay manifest write -> reconstruct field-by-field" begin
        body = make_dirichlet_diamond_body(; nspan=3)
        sig = [0.11, 0.17, 0.23, 0.31]
        wake = pnl.FilamentParticleWake(body; nwakerows=3, max_particles=128,
            core_size=0.07,
            method_trailing=pnl.StationSigmaOverlap([sig], 2.75),
            shed_with_induced_velocity=false, freestream_convection=true,
            particle_core_size=0.09)

        meta = pnl._wake_manifest_dict(wake, 1)
        @test meta["type"] == "FilamentParticleWake"
        @test meta["nwakerows"] == 3
        @test meta["core_size"] == 0.07
        @test meta["filament_core_size"] == [sig]
        @test meta["shed_with_induced_velocity"] == false
        @test meta["freestream_convection"] == true
        @test meta["particle_core_size"] == 0.09
        @test meta["method_trailing"]["type"] == "StationSigmaOverlap"

        # shedding-method (de)serialization round-trips
        m2 = pnl._deserialize_wake_shedding(meta["method_trailing"])
        @test m2 isa pnl.StationSigmaOverlap
        @test m2.sigmas == [sig]
        @test m2.overlap == 2.75
        om = pnl.OmitStations(pnl.StationSigmaOverlap([sig], 2.75),
            [[true, false, false, false]])
        om2 = pnl._deserialize_wake_shedding(pnl._wake_shedding_manifest(om))
        @test om2 isa pnl.OmitStations
        @test om2.method isa pnl.StationSigmaOverlap
        @test om2.method.sigmas == [sig]
        @test om2.omit == om.omit

        # full metadata replay reconstructs the wake
        frames = pnl.ReferenceFrame(body)
        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0; overwrite=true)
        pnl.write_vtk(joinpath(path, "run_wake1"), wake, 0, 0.0; overwrite=true)
        pnl._write_metadata_toml(path, "run", (body,), (wake,), frames,
            [0.0, 0.1], (pnl.Backslash(body),), pnl.DirectBackend(),
            pnl.DirectBackend(), pnl.DirectBackend(), ())
        pnl._append_metadata_step_toml(path, "run", frames, 0, 0.0)

        result = pnl.replay(path, "run"; steps=0, recompute=())
        w2 = result.wakes[1]
        @test w2 isa pnl.FilamentParticleWake
        @test pnl._logical_nwakerows(w2.sheet) == 3
        @test w2.sheet.core_size == 0.07
        @test w2.sheet.filament_core_size == [sig]
        @test w2.sheet.shed_with_induced_velocity == false
        @test w2.sheet.freestream_convection == true
        @test w2.particle_core_size == 0.09
        @test w2.method_trailing isa pnl.StationSigmaOverlap
        @test w2.method_trailing.sigmas == [sig]
        @test w2.method_trailing.overlap == 2.75

        # a manifest without nwakerows must fail loudly (no default)
        broken = pnl._wake_manifest_dict(wake, 1)
        delete!(broken, "nwakerows")
        toml_path = joinpath(path, "run.metadata.toml")
        md = TOML_FW.parsefile(toml_path)
        md["wake"] = [broken]
        open(toml_path, "w") do io
            TOML_FW.print(io, md)
        end
        @test_throws ArgumentError pnl.replay(path, "run"; steps=0, recompute=())
    end

    @testset "warmstart .vts save/load round-trip" begin
        wake = make_filament_fixture(; nwakerows=3, wrap=false,
            max_particles=256)
        sheet = wake.sheet
        sheet.overflowed[] = true
        for i in 1:4
            FLOWVPM.add_particle(wake.pfield, (0.3i, 0.1, 0.2),
                (0.01, 0.0, 0.02i), 0.1 + 0.01i)
        end
        np = wake.pfield.np

        path = mktempdir()
        pnl.write_vtk(joinpath(path, "w1"), wake, 5, 0.25)
        # the filament polyline companion is written alongside the .vts sheet
        @test isdir(joinpath(path, "w1_filaments"))

        wake2 = make_filament_fixture(; nwakerows=3, wrap=false,
            max_particles=256)
        wake2.sheet.nodes[1] .= 0
        wake2.sheet.strength[1] .= 0
        wake2.sheet.nwakes[] = 0
        pnl._load_panel_particle_wake_vtk!(wake2, path, "w1", 5)

        @test wake2.sheet.nwakes[] == sheet.nwakes[]
        @test wake2.sheet.nodes[1] ≈ sheet.nodes[1]
        # the .vts stores active rows only; the stale strength row nwakes+1 is
        # never read by the filament sheet (stations read rows 1:nwakes) and
        # restores as zeros — same lossy-stale-row semantics as PanelWake
        nw = sheet.nwakes[]
        @test wake2.sheet.strength[1][:, 1:nw, :] ≈ sheet.strength[1][:, 1:nw, :]
        @test wake2.pfield.np == np
        @test wake2.pfield.particles[1:3, 1:np] ≈
              wake.pfield.particles[1:3, 1:np]
        @test wake2.pfield.particles[FLOWVPM.GAMMA_INDEX, 1:np] ≈
              wake.pfield.particles[FLOWVPM.GAMMA_INDEX, 1:np]
        @test wake2.pfield.particles[FLOWVPM.SIGMA_INDEX, 1:np] ≈
              wake.pfield.particles[FLOWVPM.SIGMA_INDEX, 1:np]
    end
end
