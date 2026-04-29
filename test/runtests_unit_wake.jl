using Test
import FLOWPanel as pnl
import FLOWVPM

@testset verbose=true "Free Wakes" begin
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
