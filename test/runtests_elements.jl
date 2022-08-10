using Test
import ForwardDiff as FD
import Printf: @printf
import FLOWPanel as pnl

verbose = true
v_lvl = 0


@testset verbose=verbose "Element Tests" begin

    # --------------- ∇ϕ = u TESTS ---------------------------------------------
    if verbose
        println("\t"^(v_lvl)*"∇ϕ = u element tests")
        @printf "%s%10.10s %-26s %-26s\t%-10s\n" "\t"^(v_lvl+1) "Element" "∇ϕ (automatic diff)" "u (analytic)" "‖∇ϕ - u‖/‖u‖"
    end

    for (lbl, (phi, U)) in (
                            ("source", (pnl.phi_constant_source, pnl.U_constant_source)),
                            ("doublet", (pnl.phi_constant_doublet, pnl.U_constant_doublet)),
                           )
        @test begin

            # Nodes of planar panel
            nodes = Array{Float64, 1}[
                        [0, 0, 0],
                        [0, 1, 0],
                        [0.4, 1.2, 0],
                        [1, 0.8, 0],
                        [0.8, 0, 0]
                    ]

            # Translate and rotate panel
            O = ones(3)
            Oaxis = pnl.gt.rotation_matrix(30, 10, 40)
            nodes = [Oaxis'*node + O for node in nodes]

            # Unitary panel's strength
            strength = 1.0

            function fwrap(X::Array{R, 1}) where {R}
                targets = reshape(X, 3, 1)
                phis = zeros(R, 1)
                phi(nodes, strength, targets, phis)
                return phis
            end

            Nin = 3                        # Number of inputs
            Nout = 1                       # Number of outputs

            X = rand(3)                    # Point to evaluate
            J = zeros(Nout, Nin)           # Stores the Jacobian here, J[i,j]=dfi/dxj

            cfg = FD.JacobianConfig(nothing, X, FD.Chunk{Nin}())
            FD.jacobian!(J, fwrap, X, cfg)

            Udiff = reshape(J, 3, 1)          # Velocity from automatic differentiation

            # Analytic velocity
            Uana = zeros(3, 1)
            U(nodes, strength, reshape(X, 3, 1), Uana)

            if verbose
                @printf "%s%10.10s [% 6.3g, % 6.3g, % 6.3g] [% 6.3g, % 6.3g, % 6.3g]\t%4.3g ﹪\n" "\t"^(v_lvl+1) lbl Udiff... Uana... pnl.norm(Uana-Udiff)/pnl.norm(Uana)*100
            end

            res = abs.(Udiff .- Uana ) ./ abs.(Uana) *100 .<= 2e0

            # Test result
            prod(res)
        end
    end

end
