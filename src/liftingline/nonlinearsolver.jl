#=##############################################################################
# DESCRIPTION
    Nonlinear solver for lifting line method with sweep and dihedral.

    This solver is based on the following references:

    * Martinez-Tossas, L. A., Allaerts, D., Branlard, E., and Churchfield, M. J. 
        (2025), ASME. J. Fluids Eng., "A Solution Method for the Filtered 
        Lifting Line Theory."

    * Cory D. Goates and Douglas F. Hunsaker (2023), Journal of Aircraft, 
        "Modern Implementation and Evaluation of Lifting-Line Theory for Complex 
        Geometries" 

    * Jackson T. Reid (2020), PhD Dissertation, "A General Approach to 
        Lifting-Line Theory, Applied to Wings With Sweep"

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Nov 2025
  * License     : MIT License
=###############################################################################

"""
Robust, fast, and ForwardDiff compatible (though solver might be a bit noise, so 
set optimizer tol ~5e-5). Often it returns the secondary solution that is 
unphysical post stall. Not CSDA compatible.
"""
const optimization_solver = NonlinearSolve.FastShortcutNonlinearPolyalg(; 
                                autodiff = ADTypes.AutoForwardDiff(), 
                                vjp_autodiff = ADTypes.AutoForwardDiff(), jvp_autodiff = ADTypes.AutoForwardDiff(), 
                                prefer_simplenonlinearsolve = Val(true)
                                )

"""
Very robust and physically accurate, but it can take a long time. 
Not ForwardDiff nor CSDA compatible.
"""
const analysis_solver = NonlinearSolve.FastShortcutNLLSPolyalg(; 
                                autodiff = ADTypes.AutoForwardDiff(), 
                                # vjp_autodiff = ADTypes.AutoForwardDiff(), 
                                # jvp_autodiff = ADTypes.AutoForwardDiff(), 
                                )


function solve(self::LiftingLine, Uinfs; 
                        solver_cache=Dict(),
                        debug=false,
                        update_states=tuple,
                        optargs...)


    # Generate residual function of form `residual!(r, y, x, p)` for ImplicitAD
    residual! = generate_f_residual(self, Uinfs, update_states; 
                                        cache=solver_cache, debug)

    result = _solve(self, Uinfs, residual!; update_states, optargs...)
   
    return result, solver_cache
end

                                
function _solve(self::LiftingLine, Uinf::AbstractVector, 
                                            args...; optargs...) 
    return _solve(self, repeat(Uinf, 1, self.nelements), args...; optargs...)
end

function _solve(self::LiftingLine{R}, 
                        Uinfs::AbstractMatrix,
                        residual!::Function;
                        aoas_initial_guess=0.0,
                        align_joints_with_Uinfs=false,
                        addfields=true, raise_warn=false,
                        solver=SimpleNonlinearSolve.SimpleDFSane(),
                        solver_optargs=(; abstol = 1e-9),
                        solver_cache=Dict(),
                        Dinfs=Uinfs,
                        update_states=tuple,
                        optargs...
                        ) where {R<:Number}

    # Align joint nodes with freestream
    if align_joints_with_Uinfs

        # jointerize!(self)
        align_joints_with_Uinfs!(self, Uinfs)
        

    # NOTE: we comment jointerize! out to avoid recomputing Geffnowake every 
    #       time, but this also means that if align_joints_with_Uinfs! is ever
    #       called on the lifting line, the joints will stay aligned with Uinf
    #       afterward

    # else
    #     jointerize!(self)
    
    end

    # Update semi-infinite wake to align with freestream
    calc_Dinfs!(self, Dinfs)

    # # Generate residual function of form `residual!(r, y, x, p)` for ImplicitAD
    # residual! = generate_f_residual(self, Uinfs, update_states; 
    #                                     cache=solver_cache, debug)

    # Wrap residual for NonlinearSolve of form `f!(du, u, p)`. Note that the
    # `p` parameter used by NonlinearSolve is the `x` parameter expected by
    # ImplictAD, while there is no `p` parameter for ImplicitAD (in reality,
    # `p` in ImplicitAD is the entire geometry inside the lifting line that is 
    # already being used when calling `generate_f_residual`).
    f!(du, u, p) = residual!(du, u, p, nothing)

    # Define solver initial guess
    u0 = zeros(R, self.nelements)
    u0 .= aoas_initial_guess

    # Define nonlinear problem, using p = Uinfs
    isinplace = true
    problem = SimpleNonlinearSolve.NonlinearProblem{isinplace}(f!, u0, Uinfs)

    # Call nonlinear solver
    result = SimpleNonlinearSolve.solve(problem, solver; solver_optargs...)

    # Set solved AOA
    self.aoas .= result.u

    # Update element states from the AOA solution
    update_states(self.elements_settings, self, self.aoas, Uinfs)

    # Calculate Gamma from AOA
    # NOTE: We should use U instead of Uinf, but there isn't a clear way of 
    #       calculating U without iterating on Gamma until convergence.
    #       Just using Uinf might be good enough of an approximation?
    # UPDATE: Not needed since we already backtrack Uind from Uinf and the 
    #       effective AOA?
    # calc_Gammas!(self.Gammas, self, self.aoas, self.Us, self.elements_settings)
    calc_Gammas!(self.Gammas, self, self.aoas, Uinfs, self.elements_settings)

    # Calculate velocity at lifting-line midpoints
    self.Us .= Uinfs
    # Uind!(self, self.midpoints, self.Us)
    selfUind!(self, self.Us)

    if addfields
        gt.add_field(self.grid, "Uinf", "vector", collect(eachcol(Uinfs)), "cell"; raise_warn)
        gt.add_field(self.grid, "Gamma", "scalar", self.Gammas, "cell"; raise_warn)
        gt.add_field(self.grid, "angleofattack", "scalar", self.aoas, "cell"; raise_warn)
    end
    

    return result
end

"""
Residual = aoa_input - aoa_effective
"""
function calc_residuals!(residuals::AbstractVector, 
                            ll::LiftingLine, 
                            Uinfs::AbstractMatrix, 
                            aoas::AbstractVector, 
                            Gammas::AbstractVector, sigmas::AbstractVector, 
                            Us::AbstractMatrix,
                            elements_settings::AbstractMatrix,
                            update_states
                            )

    # NOTE: remember to update initial guess aoas and initial state Us before 
    # calling this function

    # --------- Steps 1 and 2: input aoa -> lookup cl, convert cl to Gamma -----
    # calc_Gammas!(Gammas, ll, aoas, Us, elements_settings)
    calc_Gammas!(Gammas, ll, aoas, Uinfs, elements_settings) # <--- We can use Uinf instead of U to reduce nonlinearity

    # --------- Step 3: compute inflow velocity U at lifting line --------------
    
    # Dragging line component
    for ei in 1:ll.nelements
    # Threads.@threads for ei in 1:ll.nelements

        sweep = calc_sweep(ll, ei)

        # Calculate drag coefficient
        cd = calc_cd(ll.elements[ei], aoas[ei], view(elements_settings, ei, :)...)

        # Project the velocity onto the filament direction
        UsΛ = Uinfs[1, ei]*ll.lines[1, ei] + Uinfs[2, ei]*ll.lines[2, ei] + Uinfs[3, ei]*ll.lines[3, ei]

        # Substract filament-component from the velocity
        # NOTE: We can use Uinf instead of U to reduce nonlinearity
        UinfΛ1 = Uinfs[1, ei] - UsΛ*ll.lines[1, ei]
        UinfΛ2 = Uinfs[2, ei] - UsΛ*ll.lines[2, ei]
        UinfΛ3 = Uinfs[3, ei] - UsΛ*ll.lines[3, ei]
        magUinfΛ = sqrt(UinfΛ1^2 + UinfΛ2^2 + UinfΛ3^2)

        # Calculate local geometric twist relative to freestream
        beta = calc_sweptaoa(ll, UinfΛ1, UinfΛ2, UinfΛ3, ei)

        # Calculate inflow angle
        phi = aoas[ei] - beta

        # Calculate velocity magnitude from the inflow angle as in Algorithm 1 
        # of Martinez 2025 ASME paper
        # NOTE: this assumes that the lifting line induces no velocity in the 
        #       direction of the freestream, which doesn't really hold for
        #       dihedral or winglets
        magUΛ = magUinfΛ / cosd(phi)

        # Lifting filament
        dl1 = ll.horseshoes[1, 3, ei] - ll.horseshoes[1, 2, ei]
        dl2 = ll.horseshoes[2, 3, ei] - ll.horseshoes[2, 2, ei]
        dl3 = ll.horseshoes[3, 3, ei] - ll.horseshoes[3, 2, ei]
        magdl = sqrt(dl1^2 + dl2^2 + dl3^2)

        # Velocity counter-projected on the filament direction
        Uxdl1 = UinfΛ2*dl3 - UinfΛ3*dl2
        Uxdl2 = UinfΛ3*dl1 - UinfΛ1*dl3
        Uxdl3 = UinfΛ1*dl2 - UinfΛ2*dl1
        magUxdl = sqrt(Uxdl1^2 + Uxdl2^2 + Uxdl3^2)

        # Apply the same assumption to calculate the total velocity counter-projected
        magUxdl /= math.abs_smooth(cosd(phi), 0.01)

        # Area of this section
        area = ll.chords[ei] * math.abs_smooth( dl1*ll.spans[1, ei] + dl2*ll.spans[2, ei] + dl3*ll.spans[3, ei], 0.01)

        # Calculate Lambda just like Gamma (using Eq. 41 in Goates 2022 JoA paper)
        Lambda = cd * 0.5*magUΛ^2*area / magUxdl    # Source filament strength
        sigma = Lambda/math.abs_smooth(ll.chords[ei]*cosd(sweep), 0.01)   # Equivalent constant source panel strength

        for i in 1:3
            # Here we approximate the velocity induced by the dragging line
            # by using an approximation of only the self-induced velocity
            Us[i, ei] = -0.5 * sigma^ll.sigmaexponent/2 * ll.swepttangents[i, ei]
            Us[i, ei] *= ll.sigmafactor
        end

        sigmas[ei] = sigma

    end

    # Freestream component
    Us .+= Uinfs
    # Us .= Uinfs

    # Lifting line component
    # Uind!(ll, Gammas, ll.midpoints, Us)
    selfUind!(ll, Gammas, Us)

    # Convert U to U_Λ (effective velocity seen by the swept sections)
    calc_UΛs!(ll, Us)


    for ei in 1:ll.nelements

        # --------- Step 4: compute effective swept aoa from swept inflow velocity ---
        aoa_effective = calc_sweptaoa(ll, Us, ei)

        # --------- Step 5: residual = input aoa - effective aoa ---------------
        aoa_input = aoas[ei]

        residuals[ei] = aoa_input - aoa_effective


        # # Martinez' residual

        # # Calculate local geometric twist relative to freestream
        # # beta = calc_aoa(ll, Uinfs, ei)

        # # magU = sqrt(Us[1, ei]^2 + Us[2, ei]^2 + Us[3, ei]^2)
        # # magUinf = sqrt(Uinfs[1, ei]^2 + Uinfs[2, ei]^2 + Uinfs[3, ei]^2)

        # # residuals[ei] = magU*cosd(phi) - magUinf*sind(phi)

        # # Project the velocity onto the filament direction
        # UinfsΛ = Uinfs[1, ei]*ll.lines[1, ei] + Uinfs[2, ei]*ll.lines[2, ei] + Uinfs[3, ei]*ll.lines[3, ei]
        # UinfΛ1 = Uinfs[1, ei] - UinfsΛ*ll.lines[1, ei]
        # UinfΛ2 = Uinfs[2, ei] - UinfsΛ*ll.lines[2, ei]
        # UinfΛ3 = Uinfs[3, ei] - UinfsΛ*ll.lines[3, ei]
        # magUinfΛ = sqrt(UinfΛ1^2 + UinfΛ2^2 + UinfΛ3^2)

        # # magUΛind = sqrt((Us[1, ei]-UinfΛ1)^2 + (Us[2, ei]-UinfΛ2)^2 + (Us[3, ei]-UinfΛ3)^2)
        # magUΛind = sqrt(Us[1, ei]^2 + Us[2, ei]^2 + Us[3, ei]^2)

        # beta = calc_sweptaoa(ll, UinfΛ1, UinfΛ2, UinfΛ3, ei)

        # # Calculate inflow angle
        # phi = aoas[ei] - beta

        # residuals[ei] = magUΛind*cosd(phi) - magUinfΛ*sind(phi)

    end

    # --------- Step 6: update tightly-coupled states --------------
    update_states(elements_settings, ll, aoas, Uinfs)

end

"""
Calculate effective angle of attack seen by the ei-th stripwise element
"""
function calc_aoa(U1::Number, U2::Number, U3::Number, 
                    normals::AbstractMatrix, tangents::AbstractMatrix,
                    ei::Int)

    Uy = U1*normals[1, ei] + U2*normals[2, ei]  + U3*normals[3, ei]
    Ux = U1*tangents[1, ei] + U2*tangents[2, ei]  + U3*tangents[3, ei]

    return atand(Uy, Ux)

end
function calc_aoa(Us::AbstractMatrix, normals::AbstractMatrix, tangents::AbstractMatrix, ei::Int)
    return calc_aoa(Us[1, ei], Us[2, ei], Us[3, ei], normals, tangents, ei)
end
function calc_aoa(ll::LiftingLine, U1::Number, U2::Number, U3::Number, ei::Int)
    return calc_aoa(U1, U2, U3, ll.normals, ll.tangents, ei)
end
function calc_aoa(ll::LiftingLine, Us::AbstractMatrix, args...)
    return calc_aoa(Us, ll.normals, ll.tangents, args...)
end

function calc_aoas!(aoas::AbstractVector, ll::LiftingLine, Us::AbstractMatrix)

    for ei in 1:ll.nelements
        aoas[ei] = calc_aoa(ll, Us, ei)
    end

end

"""
Calculate effective swept angle of attack at the ei-th stripwise element as in
Goates 2022, Eq. (24).
This function must be called after calc_UΛs! (Us must have been converted to 
swept velocities).
"""
function calc_sweptaoa(ll::LiftingLine, UΛs::AbstractMatrix, ei::Int)

    Uy = UΛs[1, ei]*ll.sweptnormals[1, ei] + UΛs[2, ei]*ll.sweptnormals[2, ei]  + UΛs[3, ei]*ll.sweptnormals[3, ei]
    Ux = UΛs[1, ei]*ll.swepttangents[1, ei] + UΛs[2, ei]*ll.swepttangents[2, ei]  + UΛs[3, ei]*ll.swepttangents[3, ei]

    return atand(Uy, Ux)
end
function calc_sweptaoa(ll::LiftingLine, UΛ1::Number, UΛ2::Number, UΛ3::Number, ei::Int)

    Uy = UΛ1*ll.sweptnormals[1, ei] + UΛ2*ll.sweptnormals[2, ei]  + UΛ3*ll.sweptnormals[3, ei]
    Ux = UΛ1*ll.swepttangents[1, ei] + UΛ2*ll.swepttangents[2, ei]  + UΛ3*ll.swepttangents[3, ei]

    return atand(Uy, Ux)
end

"""
Calculate nonlinear lifting line strengths from the given AOAs and inflow 
velocities.

* `aoas` are the effective AOAs seen by the swept sections (𝛼_Λ in Goates' 2022 
    notation).

* `Uinfs` is the undisturbed freestream seen by each unswept section.

"""
function calc_Gammas!(Gammas::AbstractVector, ll::LiftingLine, 
                        aoas::AbstractVector, Uinfs::AbstractMatrix,
                        elements_settings::AbstractMatrix)

    for ei in 1:ll.nelements

        sweep = calc_sweep(ll, ei)

        # Calculate swept sectional cl (C_𝐿Λ in Goates 2022, Eq. (28))
        clΛ = calc_sweptcl(ll.elements[ei], sweep, aoas[ei], view(elements_settings, ei, :)...; claero=true)

        # Project the velocity onto the filament direction
        UsΛ = Uinfs[1, ei]*ll.lines[1, ei] + Uinfs[2, ei]*ll.lines[2, ei] + Uinfs[3, ei]*ll.lines[3, ei]

        # Substract filament-component from the velocity
        UinfΛ1 = Uinfs[1, ei] - UsΛ*ll.lines[1, ei]
        UinfΛ2 = Uinfs[2, ei] - UsΛ*ll.lines[2, ei]
        UinfΛ3 = Uinfs[3, ei] - UsΛ*ll.lines[3, ei]
        magUinfΛ = sqrt(UinfΛ1^2 + UinfΛ2^2 + UinfΛ3^2)

        # Calculate local geometric twist relative to freestream
        # NOTE: Do I need to convert Uinfs to swept Uinfs and calculate beta as
        #       a swept AOA?
        # beta = calc_aoa(ll, Uinfs, ei)

        beta = calc_sweptaoa(ll, UinfΛ1, UinfΛ2, UinfΛ3, ei)

        # Calculate inflow angle
        phi = aoas[ei] - beta

        # Calculate velocity magnitude from the inflow angle as in Algorithm 1 
        # of Martinez 2025 ASME paper
        # NOTE: this assumes that the lifting line induces no velocity in the 
        #       direction of the freestream, which doesn't really hold for
        #       dihedral or winglets
        magUΛ = magUinfΛ / cosd(phi)

        # Lifting filament
        dl1 = ll.horseshoes[1, 3, ei] - ll.horseshoes[1, 2, ei]
        dl2 = ll.horseshoes[2, 3, ei] - ll.horseshoes[2, 2, ei]
        dl3 = ll.horseshoes[3, 3, ei] - ll.horseshoes[3, 2, ei]
        magdl = sqrt(dl1^2 + dl2^2 + dl3^2)

        # Velocity counter-projected on the filament direction
        # Uxdl1 = Uinfs[2, ei]*dl3 - Uinfs[3, ei]*dl2
        # Uxdl2 = Uinfs[3, ei]*dl1 - Uinfs[1, ei]*dl3
        # Uxdl3 = Uinfs[1, ei]*dl2 - Uinfs[2, ei]*dl1
        Uxdl1 = UinfΛ2*dl3 - UinfΛ3*dl2
        Uxdl2 = UinfΛ3*dl1 - UinfΛ1*dl3
        Uxdl3 = UinfΛ1*dl2 - UinfΛ2*dl1
        magUxdl = sqrt(Uxdl1^2 + Uxdl2^2 + Uxdl3^2)

        # Apply the same assumption to calculate the total velocity counter-projected
        magUxdl /= cosd(phi)

        # Area of this section
        area = ll.chords[ei] * math.abs_smooth( dl1*ll.spans[1, ei] + dl2*ll.spans[2, ei] + dl3*ll.spans[3, ei], 0.01 )

        # Calculate Gamma using Eq. 41 in Goates 2022 JoA paper
        Gammas[ei] = clΛ * 0.5*magUΛ^2*area / magUxdl

    end
end

"""
Calculate the effective velocity seeing by the swept section as in Goates 2022,
Eq. 22. This is calculated in-place, converting the velocities given in Us.
"""
calc_UΛs!(ll::LiftingLine) = calc_UΛs!(ll, ll.Us)

function calc_UΛs!(ll::LiftingLine, Us::AbstractMatrix)
    for ei in 1:ll.nelements
        calc_UΛs!(Us, ll.lines, ei)
    end
end

function calc_UΛs!(Us::AbstractMatrix, lines::AbstractMatrix, ei::Int)

    # Project the velocity onto the filament direction
    UsΛ = Us[1, ei]*lines[1, ei] + Us[2, ei]*lines[2, ei] + Us[3, ei]*lines[3, ei]

    # Substract filament-component from the velocity
    Us[1, ei] -= UsΛ*lines[1, ei]
    Us[2, ei] -= UsΛ*lines[2, ei]
    Us[3, ei] -= UsΛ*lines[3, ei]

end

"""
Generate residual wrapper for NonlinerSolver methods
"""
function generate_f_residual(ll::LiftingLine{T1},
                                Uinfs::AbstractMatrix, update_states; 
                                cache=Dict(), debug=false
                                ) where T1<:Number

    cache[:fcalls] = 0

    if debug
        cache[:residual_rms] = []
    end

    reset_cache(cache, ll.elements_settings)

    # f_residual follows the format `residual!(r, y, x, p)` of ImplicitAD
    function f_residual!(du, u::AbstractVector{T2}, Uinfs, p; 
                                cache=cache, update_states=update_states
                                ) where T2<:Number

        T = promote_type(T1, T2)

        # Fetch AOAs from input variables
        aoas = u

        # Initiate cache
        if !(T in keys(cache))
            cache[T] = (;   residuals = zeros(T, ll.nelements), 
                            Gammas = zeros(T, ll.nelements),
                            sigmas = zeros(T, ll.nelements),
                            Us = zeros(T, 3, ll.nelements),
                            fcalls = [0],
                            elements_settings = zeros(T, size(ll.elements_settings)),
                        )

            # Set Uinf as the initial velocity
            cache[T].Us .= Uinfs

            # Set initial element settings
            cache[T].elements_settings .= ll.elements_settings
        end

        # Increase function call counter
        cache[:fcalls] += 1

        cache[T].Us .= Uinfs  # <-- Force it to use only Uinfs to dimensionalize cl and cd reducing nonlinearity

        # Calculate residual
        calc_residuals!(cache[T].residuals, ll, Uinfs, 
                        aoas, cache[T].Gammas, cache[T].sigmas, cache[T].Us,
                        cache[T].elements_settings,
                        update_states)

        # Set residual as state
        du .= cache[T].residuals

        if debug
            push!( cache[:residual_rms], FD.value(sqrt(mean(du.^2))))
        end

    end

    return f_residual!

end

function reset_cache(cache, elements_settings)

    for (T, data) in cache
        if T != :fcalls && T != :residual_rms

            data.residuals .= 0
            data.Gammas .= 0
            data.sigmas .= 0
            data.Us .= 0
            data.elements_settings .= elements_settings

        end
    end

end
