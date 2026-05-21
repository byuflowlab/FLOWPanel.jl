using Test
import FLOWPanel as pnl
import FLOWVPM
using LinearAlgebra: norm
using StaticArrays: SVector

if !isdefined(@__MODULE__, :make_plate_vortex_body)
    include("test_helpers.jl")
end

@testset verbose=true "Free Wakes" begin
    @testset "SigmaOverlap chooses particle count from fixed sigma" begin
        pfield = FLOWVPM.ParticleField(10, Float64;
            fmm=FLOWVPM.FMM(autotune_reg_error=false))
        r1 = SVector(0.0, 0.0, 0.0)
        r2 = SVector(1.0, 0.0, 0.0)
        sigma = 0.2
        target_overlap = 0.7

        pnl._shed_particles!(pfield, r1, r2, 1.0, pnl.SigmaOverlap(sigma, target_overlap))

        @test pfield.np == 4
        @test all(pfield.particles[FLOWVPM.SIGMA_INDEX, 1:pfield.np] .== sigma)
        @test sigma * pfield.np / norm(r2 - r1) >= target_overlap
        @test pfield.particles[1, 1:pfield.np] ≈ [0.125, 0.375, 0.625, 0.875]
    end

    @testset "PanelWake first row shedding velocity option" begin
        body = make_plate_vortex_body()
        dt = 0.25
        uinf = [1.0, 2.0, 3.0]
        induced_first = [10.0, 20.0, 30.0]
        induced_downstream = [100.0, 200.0, 300.0]

        wake_total = pnl.PanelWake(body; nwakerows=2)
        wake_total.nwakes[] = 1
        pnl.apply_freestream!(wake_total, uinf)
        wake_total.velocity[1][:, 1, :] .+= induced_first
        wake_total.velocity[1][:, 2, :] .+= induced_downstream
        pnl.propagate!(wake_total, dt)

        @test wake_total.nodes[1][:, 1, 1] ≈ dt .* (uinf .+ induced_first)
        @test wake_total.nodes[1][:, 2, 1] ≈ dt .* (uinf .+ induced_downstream)

        wake_freestream = pnl.PanelWake(body; nwakerows=2,
            shed_with_induced_velocity=false)
        wake_freestream.nwakes[] = 1
        pnl.apply_freestream!(wake_freestream, uinf)
        wake_freestream.velocity[1][:, 1, :] .+= induced_first
        wake_freestream.velocity[1][:, 2, :] .+= induced_downstream
        pnl.propagate!(wake_freestream, dt)

        @test wake_freestream.nodes[1][:, 1, 1] ≈ dt .* uinf
        @test wake_freestream.nodes[1][:, 2, 1] ≈ dt .* (uinf .+ induced_downstream)
    end

    @testset "PanelParticleWake forwards shedding velocity option" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body; shed_with_induced_velocity=false)

        @test !wake.panel_wake.shed_with_induced_velocity
    end

    @testset "ParticleMaintenance separates mixed policies" begin
        maintenance = pnl.ParticleMaintenance((
            pnl.MinGamma(1e-3),
            pnl.MergeParticles(every=2, r_hash=0.25),
            pnl.GlobalBox([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]),
        ))

        @test maintenance.trim_policies isa Tuple{<:pnl.MinGamma,<:pnl.GlobalBox}
        @test maintenance.functional_policies isa Tuple{<:pnl.MergeParticles}
        @test_throws ArgumentError pnl.ParticleMaintenance((pnl.MinGamma(1e-3), :not_a_policy))
    end

    @testset "PanelParticleWake particle maintenance trims by strength" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.MinGamma(1e-3),)))

        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [1e-4, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [1.0, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [2.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [3.0, 0.0, 0.0], [5e-4, 0.0, 0.0], 1.0)

        pnl.propagate!(wake, 0.0; relax=false, step=1)

        @test wake.pfield.np == 2
        @test collect(FLOWVPM.get_Gamma(wake.pfield, 1)) == [2e-3, 0.0, 0.0]
        @test collect(FLOWVPM.get_Gamma(wake.pfield, 2)) == [1e-3, 0.0, 0.0]
    end

    @testset "PanelParticleWake particle maintenance trims by strength range" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.MinGamma(1e-3), pnl.MaxGamma(3e-3))))

        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [5e-4, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [1.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [2.0, 0.0, 0.0], [4e-3, 0.0, 0.0], 1.0)

        pnl.propagate!(wake, 0.0; relax=false, step=1)

        @test wake.pfield.np == 1
        @test collect(FLOWVPM.get_Gamma(wake.pfield, 1)) == [2e-3, 0.0, 0.0]
    end

    @testset "PanelParticleWake particle maintenance trims by global box" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.GlobalBox([0.0, -1.0, -1.0], [2.0, 1.0, 1.0]),)))

        FLOWVPM.add_particle(wake.pfield, [-1.0, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [2.0, 0.0, 0.0], [3e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [3.0, 0.0, 0.0], [4e-3, 0.0, 0.0], 1.0)

        pnl.propagate!(wake, 0.0; relax=false, step=1)

        @test wake.pfield.np == 2
        @test sort([FLOWVPM.get_X(wake.pfield, i)[1] for i in 1:wake.pfield.np]) == [0.0, 2.0]
    end

    @testset "PanelParticleWake particle maintenance trims by frame box" begin
        body = make_plate_vortex_body()
        frames = pnl.ReferenceFrame(body; origin=SVector(10.0, 0.0, 0.0))
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.FrameBox(1, [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]),)))

        FLOWVPM.add_particle(wake.pfield, [10.5, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [12.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)

        @test_throws ArgumentError pnl.propagate!(wake, 0.0; relax=false, step=1)

        pnl.propagate!(wake, 0.0; relax=false, step=1, frames)

        @test wake.pfield.np == 1
        @test collect(FLOWVPM.get_X(wake.pfield, 1)) == [10.5, 0.0, 0.0]
    end

    @testset "PanelParticleWake particle maintenance merges particles" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.MergeParticles(every=1, r=0.5, r_hash=0.5),)))

        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [0.1, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)

        pnl.propagate!(wake, 0.0; relax=false, step=1)

        @test wake.pfield.np == 1
    end
end
