using Test
import FLOWPanel as pnl
using StaticArrays, LinearAlgebra

const FD = pnl.FD

kernel_switch(ps, vs, gs) = pnl.FastMultipole.DerivativesSwitch(ps, vs, gs)

function assert_velocity_is_potential_gradient(eval_kernel, targets; atol=1e-10, rtol=1e-9)
    for x in targets
        _, velocity, _ = eval_kernel(x, false, true, false)
        grad_phi = FD.gradient(t -> eval_kernel(t, true, false, false)[1], x)

        @test isapprox(velocity, grad_phi; atol, rtol)
    end
end

function assert_velocity_gradient_is_velocity_jacobian(eval_kernel, targets; atol=1e-10, rtol=1e-9)
    for x in targets
        _, _, gradient = eval_kernel(x, false, false, true)
        jac_velocity = FD.jacobian(t -> eval_kernel(t, false, true, false)[2], x)

        @test isapprox(Matrix(gradient), jac_velocity; atol, rtol)
    end
end

@testset verbose=true "Kernel potential/velocity/gradient consistency" begin
    # Only scalar-potential kernels assert velocity == grad(potential). Vortex-ring
    # kernels use a vector potential in the FMM formulation; their scalar
    # potential path is a doublet-panel surrogate and is not the velocity
    # potential. Vortex rings are still checked for velocity-gradient consistency.

    vertices = (
        SVector{3,Float64}(0.0, 0.0, 0.0),
        SVector{3,Float64}(1.2, 0.1, 0.0),
        SVector{3,Float64}(0.2, 1.1, 0.3),
    )
    control_point = (vertices[1] + vertices[2] + vertices[3]) / 3
    R, _, _, _ = pnl.rotate_to_panel(vertices[1][1], vertices[1][2], vertices[1][3],
                                     vertices[2][1], vertices[2][2], vertices[2][3],
                                     vertices[3][1], vertices[3][2], vertices[3][3])
    panel_targets = (
        SVector{3,Float64}( 0.4,  0.2, 0.8),
        SVector{3,Float64}(-0.5,  0.7, 1.1),
        SVector{3,Float64}( 1.4, -0.3, 0.6),
    )

    @testset "ConstantSource triangle" begin
        strength = SVector{1,Float64}(1.3)
        eval_kernel(x, ps, vs, gs) =
            pnl._induced(x, vertices, control_point, strength, pnl.ConstantSource,
                         1e-8, R, kernel_switch(ps, vs, gs))

        assert_velocity_is_potential_gradient(eval_kernel, panel_targets)
        assert_velocity_gradient_is_velocity_jacobian(eval_kernel, panel_targets)
    end

    @testset "ConstantDoublet triangle" begin
        strength = SVector{1,Float64}(-0.7)
        eval_kernel(x, ps, vs, gs) =
            pnl._induced(x, vertices, control_point, strength, pnl.ConstantDoublet,
                         1e-8, R, kernel_switch(ps, vs, gs))

        assert_velocity_is_potential_gradient(eval_kernel, panel_targets)
        assert_velocity_gradient_is_velocity_jacobian(eval_kernel, panel_targets)
    end

    @testset "ConstantSource + ConstantDoublet triangle" begin
        strength = SVector{2,Float64}(1.3, -0.7)
        eval_kernel(x, ps, vs, gs) =
            pnl._induced(x, vertices, control_point, strength,
                         Union{pnl.ConstantSource,pnl.ConstantDoublet},
                         1e-8, R, kernel_switch(ps, vs, gs))

        assert_velocity_is_potential_gradient(eval_kernel, panel_targets)
        assert_velocity_gradient_is_velocity_jacobian(eval_kernel, panel_targets)
    end

    @testset "VortexRing triangle" begin
        strength = SVector{1,Float64}(0.9)
        eval_kernel(x, ps, vs, gs) =
            pnl._induced(x, vertices, control_point, strength, pnl.VortexRing,
                         1e-8, R, kernel_switch(ps, vs, gs))

        assert_velocity_gradient_is_velocity_jacobian(eval_kernel, panel_targets)
    end

    @testset "VortexRing control-point velocity converges to scalar-potential gradient as core shrinks" begin
        strength = SVector{1,Float64}(1.0)
        target = control_point
        core_sizes = (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5)
        errors = map(core_sizes) do core_size
            _, velocity, _ = pnl._induced(target, vertices, control_point, strength,
                                          pnl.VortexRing, core_size, R,
                                          kernel_switch(false, true, false))
            grad_phi = FD.gradient(
                t -> pnl._induced(t, vertices, control_point, strength,
                                  pnl.VortexRing, core_size, R,
                                  kernel_switch(true, false, false))[1],
                target,
            )
            norm(velocity - grad_phi) / max(norm(grad_phi), eps())
        end

        @test errors[1] > 1e-3
        @test errors[5] < 1e-9
        @test errors[end] < 1e-12
        @test all(errors[i+1] <= 0.02 * errors[i] for i in 1:6)
    end

    @testset "ConstantSource + VortexRing triangle" begin
        strength = SVector{2,Float64}(1.3, 0.9)
        eval_kernel(x, ps, vs, gs) =
            pnl._induced(x, vertices, control_point, strength,
                         Union{pnl.ConstantSource,pnl.VortexRing},
                         1e-8, R, kernel_switch(ps, vs, gs))

        assert_velocity_gradient_is_velocity_jacobian(eval_kernel, panel_targets)
    end

    @testset "ConstantDoublet wake quad" begin
        quad_vertices = (
            SVector{3,Float64}(0.0, 0.0, 0.0),
            SVector{3,Float64}(1.0, 0.2, 0.1),
            SVector{3,Float64}(1.2, 0.6, 0.8),
            SVector{3,Float64}(0.1, 0.4, 0.9),
        )
        strength = SVector{1,Float64}(0.8)
        eval_kernel(x, ps, vs, gs) =
            pnl._induced_quad(x, quad_vertices, strength, pnl.ConstantDoublet,
                              1e-8, kernel_switch(ps, vs, gs))

        assert_velocity_is_potential_gradient(eval_kernel, panel_targets)
        assert_velocity_gradient_is_velocity_jacobian(eval_kernel, panel_targets)
    end

    @testset "VortexRing wake quad" begin
        quad_vertices = (
            SVector{3,Float64}(0.0, 0.0, 0.0),
            SVector{3,Float64}(1.0, 0.2, 0.1),
            SVector{3,Float64}(1.2, 0.6, 0.8),
            SVector{3,Float64}(0.1, 0.4, 0.9),
        )
        strength = SVector{1,Float64}(0.8)
        eval_kernel(x, ps, vs, gs) =
            pnl._induced_quad(x, quad_vertices, strength, pnl.VortexRing,
                              1e-8, kernel_switch(ps, vs, gs))

        assert_velocity_gradient_is_velocity_jacobian(eval_kernel, panel_targets)
    end

    @testset "VortexRing filament (_bound_vortex_gradient)" begin
        va = SVector{3,Float64}( 1.0,  0.0,  0.0)
        vb = SVector{3,Float64}( 0.0,  1.0,  0.0)
        core = 0.05

        # Off-singularity targets, on-edge proximity, and a point inside the core radius.
        targets = (
            SVector{3,Float64}( 0.3,  0.3,  0.5),       # well off filament
            SVector{3,Float64}( 0.3,  0.3,  0.005),     # very close to the filament plane
            SVector{3,Float64}( 0.5,  0.5,  0.02),      # within ~core of the segment
            SVector{3,Float64}(-2.0,  3.0, -1.5),       # far away
        )

        for finite_core in (true, false), x in targets
            J_an = pnl._bound_vortex_gradient(va - x, vb - x, finite_core, core)
            J_fd = FD.jacobian(t -> pnl._bound_vortex_velocity(va - t, vb - t, finite_core, core), x)
            @test isapprox(Matrix(J_an), J_fd; atol=1e-10, rtol=1e-9)
        end
    end

    @testset "Semi-infinite ConstantDoublet wake" begin
        geom = (0.0, 0.0, 0.0,
                1.0, 0.2, 0.1,
                1.0, 0.0, 0.0,
                1.2)
        targets = (
            SVector{3,Float64}( 0.35, -0.4, 0.8),
            SVector{3,Float64}( 0.20,  1.1, 0.7),
            SVector{3,Float64}(-0.50,  0.6, 1.2),
        )
        eval_kernel(x, ps, vs, gs) =
            pnl.induced_semiinfinite(x, pnl.ConstantDoublet, geom...,
                                     kernel_switch(ps, vs, gs); kerneloffset=1e-8)

        assert_velocity_is_potential_gradient(eval_kernel, targets)
        assert_velocity_gradient_is_velocity_jacobian(eval_kernel, targets)
    end

    @testset "Semi-infinite VortexRing wake alias" begin
        geom = (0.0, 0.0, 0.0,
                1.0, 0.2, 0.1,
                1.0, 0.0, 0.0,
                1.2)
        targets = (
            SVector{3,Float64}( 0.35, -0.4, 0.8),
            SVector{3,Float64}(-0.50,  0.6, 1.2),
        )
        eval_kernel(x, ps, vs, gs) =
            pnl.induced_semiinfinite(x, pnl.VortexRing, geom...,
                                     kernel_switch(ps, vs, gs); kerneloffset=1e-8)

        assert_velocity_gradient_is_velocity_jacobian(eval_kernel, targets)
    end

end
