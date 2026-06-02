using Test
import FLOWPanel as pnl
using Statistics: mean

@testset verbose=true "Analytical Validation" begin
    @testset "Sphere in uniform flow - Neumann" begin
        body = make_sphere_source_body(ntheta=18, nphi=36)
        body, pvals, Fs = solve_source_body!(body)

        cps = body.controlpoints
        xcoords = view(cps, 1, :)
        ycoords = view(cps, 2, :)
        zcoords = view(cps, 3, :)

        # With rho=1, Uinf=1: q = 0.5, so P_stag = q*1.0 = 0.5, P_equator = q*(-1.25) = -0.625
        i_stag = argmax(xcoords)
        @test isapprox(pvals[i_stag], 0.5; atol=0.025)

        eq_mask = abs.(xcoords) .< 0.05
        @test any(eq_mask)
        p_equator = mean(pvals[eq_mask])
        @test isapprox(p_equator, -0.625; atol=0.1)

        q = 0.5 * pi
        ftot = pnl.calcfield_Ftot!(zeros(3), body, Fs)
        @test abs(ftot[1]) / q < 0.05
        @test abs(ftot[3]) / q < 0.05

        pos_y = pvals[ycoords .> 0]
        neg_y = pvals[ycoords .< 0]
        @test !isempty(pos_y)
        @test !isempty(neg_y)
        @test abs(mean(pos_y) - mean(neg_y)) < 0.01

        pos_z = pvals[zcoords .> 0]
        neg_z = pvals[zcoords .< 0]
        @test !isempty(pos_z)
        @test !isempty(neg_z)
        @test abs(mean(pos_z) - mean(neg_z)) < 0.01

        @testset "Cp pointwise comparison" begin
            q = 0.5 * 1.0 * 1.0^2  # 0.5*rho*Uinf^2
            Cp = pvals ./ q
            # Analytic: Cp = 1 - (9/4)*sin²θ; on the unit sphere sin²θ = y² + z²
            Cp_analytic = [1 - (9/4) * (cps[2, i]^2 + cps[3, i]^2) for i in 1:size(cps, 2)]
            abs_err = abs.(Cp .- Cp_analytic)
            @test maximum(abs_err) < 0.15
            @test mean(abs_err) < 0.05
            @show maximum(abs_err)
            @show mean(abs_err)
            @show minimum(abs_err)
        end
    end
end
