using LinearAlgebra
using StaticArrays
using Test

import FLOWPanel as pnl

const FD = pnl.FD

function erf_approx(x)
    # Abramowitz-Stegun 7.1.26; sufficient for this numerical divergence check
    # and compatible with ForwardDiff dual numbers.
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911

    s = ifelse(x < zero(x), -one(x), one(x))
    z = abs(x)
    t = inv(1 + p * z)
    y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-z^2)

    return s * y
end

function gaussian_vortex_particle_velocity(x, xp, Gamma, sigma)
    r = x - xp
    rho = norm(r)
    eta = rho / sigma

    F = (
        erf_approx(eta / sqrt(2)) -
        sqrt(2 / pi) * eta * exp(-eta^2 / 2)
    ) / (4 * pi * rho^3)

    return cross(Gamma, r) * F
end

function finite_core_vortex_ring_velocity(x, vertices, strength, core_size)
    velocity = zero(x)
    nvertices = length(vertices)

    for i in 1:nvertices
        va = vertices[i]
        vb = vertices[mod1(i + 1, nvertices)]
        velocity += pnl._bound_vortex_velocity(va - x, vb - x, true, core_size) * strength
    end

    return velocity
end

function divergence(velocity, x)
    J = FD.jacobian(velocity, x)
    return tr(Matrix(J))
end

@testset verbose=true "Regularized vortex kernel divergence" begin
    @testset "Gaussian vortex particle" begin
        xp = SVector(0.2, -0.1, 0.4)
        Gamma = SVector(0.7, -1.1, 0.5)
        sigma = 0.35
        targets = (
            SVector(0.8, 0.3, -0.2),
            SVector(-0.4, 0.9, 0.6),
            SVector(0.25, -0.05, 0.55),
            SVector(2.0, -1.0, 1.4),
        )

        for x in targets
            divu = divergence(t -> gaussian_vortex_particle_velocity(t, xp, Gamma, sigma), x)
            @info "Gaussian vortex particle divergence" x divu
            @test abs(divu) < 1e-12
        end
    end

    @testset "Vatistas finite-core vortex ring panel" begin
        vertices = (
            SVector(0.0, 0.0, 0.0),
            SVector(1.2, 0.1, 0.0),
            SVector(1.1, 1.0, 0.2),
            SVector(-0.1, 0.8, -0.1),
        )
        strength = 0.8
        core_size = 0.05
        targets = (
            SVector(0.3, 0.3, 0.7),
            SVector(0.5, 0.4, 0.04),
            SVector(-0.8, 0.2, 0.5),
            SVector(2.0, -1.5, 1.1),
        )

        for x in targets
            divu = divergence(t -> finite_core_vortex_ring_velocity(t, vertices, strength, core_size), x)
            @info "Vatistas finite-core vortex ring divergence" x divu
            @test abs(divu) < 1e-10
        end
    end
end
