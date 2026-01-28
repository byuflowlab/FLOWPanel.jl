#=##############################################################################
# DESCRIPTION
    Auxiliary functions for optimization

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Jan 2026
  * License     : MIT License
=###############################################################################


"""
Arctan2 of complex numbers to accomodate complex-step derivative approximation.
"""
function Base.atan(y::Complex, x::Complex)
    
    re = atan(real(y), real(x))
    ima = (real(x)*imag(y) - real(y)*imag(x)) / (real(y)*real(y) + real(x)*real(x))
    
    return Complex(re, ima)
end

"""
Comparison of complex numbers to accomodate complex-step derivative approximation
"""
Base.isless(x::Complex, y::Real) = isless(real(x), y)
Base.isless(x::Real, y::Complex) = isless(x, real(y))
Base.isless(x::Complex, y::Complex) = isless(real(x), real(y))

"""
Override `abs(::Complex)` to be consistent with CSDA instead of the complex 
norm. This is important in order to get NonLinearSolve to converge the CSDA 
primal at the same accuracy than Dual.
"""
Base.abs(x::Complex) = sqrt(x^2)

"""
Wrapper function that generates a lifting line (or morphs an existing one found
under the `cache` optional argument), and evaluates and returns forces, moments, 
and stability derivatives.
"""
function run_liftingline(;
        stability_derivatives = false,                  # Whether to calculate stability derivatives
        
        save_path       = nothing,                      # Where to save outputs
        
        run_name        = "liftingline",                # Name of this run
        airfoil_path    = joinpath(examples_path, "data"), # Where to find 2D polars
        
        verbose         = true,
        paraview        = false,                        # Whether to visualize with Paraview
        
        
        # ----------------- SIMULATION PARAMETERS --------------------------------------
        alpha::R1       = 4.2,                          # (deg) angle of attack
        beta::R2        = 0.0,                          # (deg) sideslip angle
        
        magUinf::R3     = 49.7,                         # (m/s) freestream velocity magnitude
        ground_distance = Inf,                          # (m) distance to ground
        
        rho::R4         = 1.225,                        # (kg/m^3) air density
        
        # ------------------ GEOMETRY PARAMETERS ---------------------------------------
        
        # High-level parameters
        b::R5           = 98*0.0254,                    # (m) wing span
        
        # Discretization parameters
        nelements       = 40,                           # Number of stripwise elements per semi-span
        discretization  = [                             # Multi-discretization of wing (seg length, ndivs, expansion, central)
                            (1.00,  nelements, 1/10, false),
                          ],
        symmetric       = true,                         # Whether the wing is symmetric
        
        # Chord distribution (nondim y-position, nondim chord)
        chord_distribution::AbstractMatrix{R6} = [
        #   2*y/b   c/b
            0.0     0.2
            1.0     0.2
        ],
        
        # Twist distribution (nondim y-position, twist)
        twist_distribution::AbstractMatrix{R7} = [
        #   2*y/b   twist (deg)
            0.0     0.0
            1.0     0.0
        ],
        
        # Sweep distribution (nondim y-position, sweep)
        sweep_distribution::AbstractMatrix{R8} = [
        #   2*y/b   sweep (deg)
            0.0     45.0
            1.0     45.0
        ],
        
        # Dihedral distribution (nondim y-position, dihedral)
        dihedral_distribution::AbstractMatrix{R9} = [
        #   2*y/b   dihedral (deg)
            0.0     0.0
            1.0     0.0
        ],
        
        # Span-axis distribution: chordwise point about which the wing is twisted, 
        # swept, and dihedralized (nondim y-position, nondim chord-position)
        spanaxis_distribution::AbstractMatrix{R10} = [
        #   2*y/b   x/c
            0.0     0.25
            1.0     0.25
        ],
        
        # Airfoil contour distribution (nondim y-position, polar, airfoil type)
        airfoil_distribution = [
        #    2*y/b  polar file                        airfoil type
            (0.00, "rae101-Re1p7e6-smooth180-2.csv",  SimpleAirfoil),
            (1.00, "rae101-Re1p7e6-smooth180-2.csv",  SimpleAirfoil)
        ],
        
        element_optargs = (;    path = airfoil_path,
                                extrapolatepolar = false,   # Whether to extrapolate the 2D polars to ±180 deg
                                plot_polars = !true,
                                verbose
                            ),
        
        X0              = [0.25 * chord_distribution[1, 2]*b, 0, 0], # Center about which to calculate moments
        
        
        # ------------------ SOLVER PARAMETERS -----------------------------------------
        deltasb         = 1.0,                          # Blending distance, deltasb = 2*dy/b
        deltajoint      = 1.0,                          # Joint distance, deltajoint = dx/c
        
        sigmafactor     = 0.0,                          # Dragging line amplification factor (set to -1.0 for rebust post-stall method)
        sigmaexponent   = 4.0,                          # Dragging line amplification exponent (no effects if `sigmafactor==0.0`)
        
                                                        # Nonlinear solver
        # solver        = SimpleNonlinearSolve.SimpleDFSane(),             # Indifferent to initial guess, but somewhat not robust post stall   <---- NOT COMPATIBLE WITH FORWARDDIFF (but compatible with CSDA and optimization converges well)
        # solver        = SimpleNonlinearSolve.SimpleTrustRegion(),        # Trust region needs a good initial guess, but solver converges very reliably post stall, not compatible with CSDA nor ForwardDiff
        # solver        = NonlinearSolve.SimpleBroyden(),                  # Optimization converges well while being compatible with ForwardDiff (not compatible with CSDA). EXTREMELY ROBUST ACROSS LINEAR, MILD STALL, AND DEEP STALL (it returns an answer, but it might be noise and not accurate)
        # solver        = NonlinearSolve.NonlinearSolveQuasiNewton.Broyden(; autodiff = ADTypes.AutoForwardDiff(),  init_jacobian=Val(:true_jacobian)) # Also extremely robust across regions (it returns an answer, but it might be noise and not accurate)
        # solver        = NonlinearSolve.NLsolveJL(method = :trust_region),# Optimization converges very well with ForwardDiff, not compatible with CSDA. Solver converges slowly but realibly in linear and mild stall regions, does not converge post stall
        # solver        = NonlinearSolve.SIAMFANLEquationsJL(method = :newton, autodiff=ADTypes.AutoForwardDiff()), # Also robust in linear and mild stall regions, but much faster
        # solver        = NonlinearSolve.NonlinearSolve.FastShortcutNLLSPolyalg(; autodiff = ADTypes.AutoForwardDiff(), vjp_autodiff = ADTypes.AutoForwardDiff(), jvp_autodiff = ADTypes.AutoForwardDiff()) # Very robust, but it can take a long time. Not ForwardDiff nor CSDA compatible
        # solver        = NonlinearSolve.NonlinearSolve.FastShortcutNonlinearPolyalg(; autodiff = ADTypes.AutoForwardDiff(), vjp_autodiff = ADTypes.AutoForwardDiff(), jvp_autodiff = ADTypes.AutoForwardDiff(), prefer_simplenonlinearsolve = Val(true)) # 100% convergence, and super fast, but it might return the wrong solution. ForwardDiff compatible (though solver might be a bit noise, so set optimizer tol ~5e-5)
        solver          = optimization_solver,
        
        solver_optargs  = (; 
                            abstol = 1e-13,  
                            maxiters = 800,
                            ),
        
        align_joints_with_Uinfs = false,                # Whether to align joint bound vortices with the freestream
        
        use_Uind_for_force = true,                      # Whether to use Uind as opposed to selfUind for force postprocessing
                                                        # (`true` for more accurate spanwise cd distribution, but worse integrated CD)

        distributions = false,                          # Whether to output spanwise distributions

        cache = Dict(),                                 # Model cache
        
        plot_convergence = false,                       # Whether to plot solver convergence
        
) where {R1, R2, R3, R4, R5, R6, R7, R8, R9, R10}

    # Populate cache fields
    for field in ("ll", "Uinfs", "solver_cache")
        if !(field in keys(cache))
            cache[field] = Dict()
        end
    end

    if !("fcalls" in keys(cache))
        cache["fcalls"] = 0
    end

    # Increase function call counter
    cache["fcalls"] += 1
                                                    
    # Number type for LiftingLine
    NumType = promote_type(R1, R2, R3, R4, R5, R6, R7, R8, R9, R10)  
    
    # NOTE: In these dual numbers we will implicitely defined the first partial to be 
    #       the derivative w.r.t. angle of attack and the second partial to be
    #       the derivative w.r.t. sideslip angle
    if stability_derivatives

        # Promote the primals of alpha and beta
        alpha = NumType(alpha)
        beta = NumType(beta)

        # Convert alpha and beta into Dual numbers for stability derivatives
        tag = FD.Tag{:stabilityderivative, NumType}
        alpha = FD.Dual{tag}(alpha, FD.Partials((1.0, 0.0)))    # Convert angle of attack into dual number for automatic differentiation
        beta  = FD.Dual{tag}(beta,  FD.Partials((0.0, 1.0)))    # Convert sideslip angle into dual number for automatic differentiation
        
        # Number type for LiftingLine
        NumType = promote_type(typeof(alpha), typeof(beta))
        
    end
    
    Uinf            = magUinf*direction(; alpha, beta) # Freestream vector

    Dhat            = direction(; alpha)        # Drag direction
    Shat            = [0, 1, 0]                     # Span direction
    Lhat            = cross(Dhat, Shat)             # Lift direction
    
    lhat            = Dhat                          # Rolling direction
    mhat            = Shat                          # Pitching direction
    nhat            = Lhat                          # Yawing direction
    
    cref            = chord_distribution[1, 2]*b    # (m) reference chord

    # Initialize solver cache
    if !(NumType in keys(cache["solver_cache"]))
        cache["solver_cache"][NumType] = Dict()
    end 


    # ------------------ GENERATE GEOMETRY -------------------------------------

    # Arguments for geometry generation
    geom_optargs = (; b, chord_distribution, twist_distribution,
                                        sweep_distribution, dihedral_distribution,
                                        spanaxis_distribution,
                                        discretization,
                                        symmetric,
                                        plot_discretization = !true,)

    # Generate or fetch Lifting Line
    if NumType in keys(cache["ll"])
        
        ll = cache["ll"][NumType]
        
    else
        
        ll = LiftingLine{NumType}(  airfoil_distribution; 
                                        geom_optargs...,
                                        deltasb, deltajoint, sigmafactor, sigmaexponent,
                                        element_optargs,
                                        )
        verbose && display(ll)

        cache["ll"][NumType] = ll
    end

    # Morph Lifting Line into its shape
    remorph!(   ll; 
                    geom_optargs...,
                    deltasb, deltajoint,
                    )


    # ------------------ CALL NONLINEAR SOLVER -------------------------------------

    # Fetch cache
    if NumType in keys(cache["Uinfs"])
        Uinfs = cache["Uinfs"][NumType]
    else
        Uinfs = repeat(Uinf, 1, ll.nelements)
        cache["Uinfs"][NumType] = Uinfs
    end
    
    # Freestream velocity at each stripwise element
    for U in eachcol(Uinfs)
        U .= Uinf
    end

    # Set ground distance
    if isfinite(ground_distance)
        set_ground(ll, ground_distance)
    end
    
    # Run solver
    result, solver_cache = solve(ll, Uinfs; 
                                        debug=plot_convergence,      # `true` returns the residual rms
                                        aoas_initial_guess=alpha, 
                                        align_joints_with_Uinfs, 
                                        solver, solver_optargs,
                                        solver_cache=cache["solver_cache"][NumType]
                                        )

    # Check solver success
    success = SimpleNonlinearSolve.SciMLBase.successful_retcode(result)

    if !success
        @warn "Lifting line solver did not converge!"
        @show success
        @show solver_cache[:fcalls]
        @show result.retcode
        # display(result)
    end

    # (Optional) plot residual convergence
    if plot_convergence
        fig = plt.figure(figsize=[7, 0.5*5]*7/9)
        ax = fig.gca()
        
        ax.plot(1:solver_cache[:fcalls], solver_cache[:residual_rms], "-k", linewidth=1)
        
        ax.set_xlabel("Iteration")
        ax.set_xlim([0, solver_cache[:fcalls]])
        
        ax.set_ylabel("Residual RMS")
        ax.set_yscale("log")
        
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)
        
        fig.tight_layout()
    end


    # ------------------ POSTPROCESSING --------------------------------------------

    distributions = distributions ? [] : nothing

    # Calculate dimensional forces and moments (vectors)
    forcesmoments = calc_forcesmoments(ll, Uinfs, Uinf, rho; 
                                                X0,
                                                use_Uind_for_force,
                                                Dhat, Shat, Lhat,
                                                lhat, mhat, nhat,
                                                distributions,
                                                )
    
    # Calculate force and moment coefficients (scalars)
    coefficients = calc_forcemoment_coefficients(ll, Uinfs, Uinf,
                                                rho, cref, b;
                                                forcesmoments,
                                                distributions,
                                                )
    
    # Unpack calculations
    (; lift, drag, side) = forcesmoments        # Forces (vectors)
    (; roll, pitch, yaw) = forcesmoments        # Moments (vectors)
    (; X0) = forcesmoments                      # Point about with the moments where calculated
    
    (; Ftot, Mtot) = forcesmoments              # Total force and moment vectors
    
    (; CD, CY, CL) = coefficients               # Drag, side, and lift forces
    (; Cl, Cm, Cn) = coefficients               # Roll, pitch, and yawing moment
    
    # (; Dhat, Shat, Lhat) = coefficients       # Direction of each force
    # (; lhat, mhat, nhat) = coefficients       # Direction of each moment
    
    (; q, Aref, bref, cref) = coefficients      # Reference dynamic pressure, area, span, and chord

    # Calculate scalar forces and moments from coefficients
    D = CD * (q*Aref)
    Y = CY * (q*Aref)
    L = CL * (q*Aref)
    l = Cl * (q*Aref*cref)
    m = Cm * (q*Aref*cref)
    n = Cn * (q*Aref*cref)

    if stability_derivatives
    
        # Fetch stability derivatives
        dCDdα = FD.partials(CD)[1]
        dCYdα = FD.partials(CY)[1]
        dCLdα = FD.partials(CL)[1]
        dCldα = FD.partials(Cl)[1]
        dCmdα = FD.partials(Cm)[1]
        dCndα = FD.partials(Cn)[1]
        
        dCDdβ = FD.partials(CD)[2]
        dCYdβ = FD.partials(CY)[2]
        dCLdβ = FD.partials(CL)[2]
        dCldβ = FD.partials(Cl)[2]
        dCmdβ = FD.partials(Cm)[2]
        dCndβ = FD.partials(Cn)[2]
    
        # Get primals of forces and moments
        CD = FD.value(CD)
        CY = FD.value(CY)
        CL = FD.value(CL)
        Cl = FD.value(Cl)
        Cm = FD.value(Cm)
        Cn = FD.value(Cn)
        
        D = FD.value(D)
        Y = FD.value(Y)
        L = FD.value(L)
        l = FD.value(l)
        m = FD.value(m)
        n = FD.value(n)
        
        lift = FD.value.(lift)
        drag = FD.value.(drag)
        side = FD.value.(side)
        roll = FD.value.(roll)
        pitch = FD.value.(pitch)
        yaw = FD.value.(yaw)
        
        Ftot = FD.value.(Ftot)
        Mtot = FD.value.(Mtot)
    
        q = FD.value(q)
        Aref = FD.value(Aref)
        bref = FD.value(bref)
        cref = FD.value(cref)
    
        Dhat = FD.value.(Dhat)
        Shat = FD.value.(Shat)
        Lhat = FD.value.(Lhat)
    
        lhat = FD.value.(lhat)
        mhat = FD.value.(mhat)
        nhat = FD.value.(nhat)

        X0 = FD.value.(X0)
    else
        dCDdα = NaN
        dCYdα = NaN
        dCLdα = NaN
        dCldα = NaN
        dCmdα = NaN
        dCndα = NaN
        
        dCDdβ = NaN
        dCYdβ = NaN
        dCLdβ = NaN
        dCldβ = NaN
        dCmdβ = NaN
        dCndβ = NaN
        
    end
    
    # ------------------ OUTPUT SOLUTION -------------------------------------------
    if !isnothing(save_path)
        
        str = save(ll, run_name; path=save_path, debug=!true) # Use `debug=true` to output the effective horseshoes
    
        if paraview
            run(`paraview --data=$(joinpath(save_path, str))`)
        end
        
    end
    
    return (;   CD, CY, CL, Cl, Cm, Cn,                       # Force and moment coefficients
        
                dCDdα, dCYdα, dCLdα, dCldα, dCmdα, dCndα,     # Stability derivatives
                dCDdβ, dCYdβ, dCLdβ, dCldβ, dCmdβ, dCndβ,
        
                q, Aref, bref, cref,                          # Non-dimensionalization parameters
        
                D, Y, L, l, m, n,                             # Dimensional forces and moments (scalars)
        
                Dhat, Shat, Lhat,                             # Direction of forces and moments (vectors)
                lhat, mhat, nhat,
        
                lift, drag, side, roll, pitch, yaw,           # Dimensional forces and moments (vectors)
                Ftot, Mtot,                                   # Total forces and moment vectors

                X0,                                           # Point about which the moments where calculated
        
                ll,                                           # Lifting line object that was used to obtain these results
                cache,                                        # Model cache

                distributions = isnothing(distributions) ? NaN : distributions[end], # Spanwise distributions

                success                                       # Whether the nonlinear solver was successful
                )
end