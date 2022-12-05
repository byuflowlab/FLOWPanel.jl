#=##############################################################################
# DESCRIPTION
    Unit tests of elements
=###############################################################################

using Test
import ForwardDiff as FD
import Printf: @printf
import FLOWPanel as pnl

try
    verbose
catch
    global verbose = true
end
v_lvl = 0

run_gradphiu = true
run_memoryalloc = true
run_divcurl = true
run_properties = true

elements_to_test = (
                        ( "source", (pnl.phi_constant_source, pnl.U_constant_source), ((offset=0,), (offset=0, )) ),
                        ( "doublet", (pnl.phi_constant_doublet, pnl.U_constant_doublet), ((), (cutoff=0, offset=0)) ),
                  )


@testset verbose=verbose "Element Tests" begin

    # --------------- ∇ϕ = u TESTS ---------------------------------------------
    if verbose
        println("\t"^(v_lvl)*"∇ϕ = u element tests")
        @printf "%s%10.10s %-24s %-26s\t%-10s\n" "\t"^(v_lvl+1) "Element" "∇ϕ (automatic diff)" "u (analytic)" "max((∇ϕ .- u)./u)"
    end

    for (lbl, (phi, U), (phioptargs, Uoptargs)) in elements_to_test[1:(run_gradphiu ? end : -1)]
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

            # Store the nodes as a 3xn matrix
            nodes = hcat(nodes...)

            # Define panel
            nnodes = size(nodes, 2)
            panel = 1:nnodes

            # Unitary panel's strength
            strength = 1.0

            function fwrap(X::Array{R, 1}) where {R}
                targets = reshape(X, 3, 1)
                phis = zeros(R, 1)
                phi(nodes, panel, strength, targets, phis; phioptargs...)
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
            U(nodes, panel, strength, reshape(X, 3, 1), Uana; Uoptargs...)

            if verbose
                @printf "%s%10.10s [% 6.3f, % 6.3f, % 6.3f] [% 6.3f, % 6.3f, % 6.3f]\t%4.3g ﹪\n" "\t"^(v_lvl+1) lbl Udiff... Uana... maximum(abs.(Uana.-Udiff)./abs.(Uana))*100
            end

            res = abs.(Udiff .- Uana ) ./ abs.(Uana) *100 .<= 1e-8

            # Test result
            prod(res)
        end
    end



    # --------------- MEMORY-ALLOCATION TESTS ----------------------------------
    if verbose
        println("\t"^(v_lvl)*"Memory-allocation element tests")
        @printf "%s%10.10s %10s %16s\n" "\t"^(v_lvl+1) "Element" "ϕ" "u"
    end

    for (lbl, (phi, U)) in elements_to_test[1:(run_memoryalloc ? end : -1)]
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

            # Store the nodes as a 3xn matrix
            nodes = hcat(nodes...)

            # Define panel
            nnodes = size(nodes, 2)
            panel = 1:nnodes

            # Unitary panel's strength
            strength = 1.0

            X = rand(3)                   # Point to evaluate
            phis = zeros(1)
            Us = zeros(3, 1)
            targets = reshape(X, 3, 1)

            # Warm-up runs
            phi(nodes, panel, strength, targets, phis)
            U(nodes, panel, strength, targets, Us)

            phi_alloc = @allocated phi(nodes, panel, strength, targets, phis)
            U_alloc = @allocated U(nodes, panel, strength, targets, Us)


            if verbose
                @printf "%s%10.10s %10i bytes %10i bytes\n" "\t"^(v_lvl+1) lbl phi_alloc U_alloc
            end

            # Test result
            phi_alloc==0 && U_alloc==0
        end
    end



    # --------------- SOLENOIDAL AND IRROTATIONAL TESTS ------------------------
    if verbose
        println("\t"^(v_lvl)*"Solenoidal and irrotational tests")
        @printf "%s%10.10s %-10s %16s\n" "\t"^(v_lvl+1) "Element" "    ∇ ⋅u" "∇ × u"
    end

    for (lbl, (phi, U), (phioptargs, Uoptargs)) in elements_to_test[1:(run_divcurl ? end : -1)]
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

            # Store the nodes as a 3xn matrix
            nodes = hcat(nodes...)

            # Define panel
            nnodes = size(nodes, 2)
            panel = 1:nnodes

            # Unitary panel's strength
            strength = 1.0

            function fwrap(X::Arr1) where {R, Arr1<:AbstractArray{R}}
                targets = reshape(X, 3, 1)
                Us = zeros(R, 3, 1)
                U(nodes, panel, strength, targets, Us; Uoptargs...)
                return Us
            end

            Nin = 3                        # Number of inputs
            Nout = 3                       # Number of outputs

            X = rand(3)                    # Point to evaluate
            J = zeros(Nout, Nin)           # Stores the Jacobian here, J[i,j]=dfi/dxj

            cfg = FD.JacobianConfig(nothing, X, FD.Chunk{Nin}())
            FD.jacobian!(J, fwrap, X, cfg)

            div = sum([J[i,i] for i in 1:3])
            curl = [J[3,2]-J[2,3], J[1,3]-J[3,1], J[2,1]-J[1,2]]

            if verbose
                @printf "%s%10.10s % 5.3e [% 6.1e, % 6.1e, % 6.1e]\n" "\t"^(v_lvl+1) lbl div curl[1] curl[2] curl[3]
            end

            # Test result
            abs(div) < 1e-14
        end
    end
end



if run_properties
    @testset verbose=verbose "Element Properties" begin

        # Nodes of planar panel
        nodes = [
                    -1 -1 0;
                     1 -1 0;
                     1  1 0;
                    -1  1 0;
                ]'

        # Define panel
        nnodes = size(nodes, 2)
        panel = 1:nnodes

        # Translate and rotate panel
        O = zeros(3)
        Oaxis = pnl.gt.rotation_matrix(0, 0, 0)
        # O = (rand(3) .- 0.5)*1.99
        # Oaxis = pnl.gt.rotation_matrix(rand(0:360), rand(0:360), rand(0:360))
        nodes = Oaxis'*nodes + repeat(O, 1, nnodes)

        # Unitary panel's strength
        strength = 1.0


        # Evaluate potential field at z=0 and z=inf
        hinf = rand([-1, 1])*1e16
        h = rand([-1, 1])*1e-8
        normal = Oaxis'*[0, 0, 1]
        offset = vcat((round.(rand(2), digits=4) .- 0.5)*1.99, 0)
        targets = hcat(0*normal+O+offset, hinf*normal+O+offset, h*normal+O+offset)


        # --------------- SOURCE PANEL PROPERTIES ------------------------------
        # Test dphi/dz(±0) == ±sigma/2
        if verbose
            println("\n"*"\t"^(v_lvl)*"Test source properties")
            println("\t"^(v_lvl+1)*"xinf:\t\t\t$(targets[:, 2])\n")
        end

        phis = zeros(3)
        pnl.phi_constant_source(nodes, panel, strength, targets, phis; offset=0)

        phiinf = phis[2]
        dphidz = (phis[3] - phis[1])/h

        phiinfres = abs(phiinf) < 1e-8
        dphidzres = abs(dphidz - sign(h)*strength/2) < 1e-6

        if verbose
            hinfsign = sign(h)==1 ? "+" : "-"
            hsign = sign(h)==1 ? "+" : "-"
            println("\t"^(v_lvl+1)*"ϕ($(hinfsign)∞):\t\t\t$(phis[2])")
            println("\t"^(v_lvl+1)*"ϕ($(hinfsign)∞) == 0 ? \t\t$(phiinfres)\n")

            println("\t"^(v_lvl+1)*"dϕ/dz(0, 0, $(hsign)0):\t$(dphidz)")
            println("\t"^(v_lvl+1)*"σ:\t\t\t$(strength)")
            println("\t"^(v_lvl+1)*"dϕ/dz($(hsign)0) == $(hsign)σ/2 ? \t$(dphidzres)")
        end

        @test phiinfres==true
        @test dphidzres==true

        # --------------- DOUBLET PANEL PROPERTIES -----------------------------
        # Test phi(z=±0) == ∓mu/2
        if verbose
            println("\n"*"\t"^(v_lvl)*"Test doublet properties")
            println("\t"^(v_lvl+1)*"xinf:\t\t\t$(targets[:, 2])\n")
        end

        phis = zeros(3)
        pnl.phi_constant_doublet(nodes, panel, strength, targets, phis)

        phiinf = phis[2]
        phidz = phis[3]

        phiinfres = abs(phiinf) < 1e-8
        phidzres = abs(phidz - sign(-h)*strength/2) < 1e-6

        if verbose
            hinfsign = sign(h)==1 ? "+" : "-"
            hsign, musign = sign(h)==1 ? ("+", "-") : ("-", "+")
            println("\t"^(v_lvl+1)*"ϕ($(hinfsign)∞):\t\t\t$(phis[2])")
            println("\t"^(v_lvl+1)*"ϕ($(hinfsign)∞) == 0 ? \t\t$(phiinfres)\n")

            println("\t"^(v_lvl+1)*"ϕ(0, 0, $(hsign)0):\t\t$(phidz)")
            println("\t"^(v_lvl+1)*"μ:\t\t\t$(strength)")
            println("\t"^(v_lvl+1)*"ϕ(0, 0, $(hsign)0) == $(musign)μ/2 ? \t$(phidzres)")
        end

        @test phiinfres
        @test phidzres

    end
end
