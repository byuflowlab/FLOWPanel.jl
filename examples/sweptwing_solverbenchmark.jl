#=##############################################################################
# DESCRIPTION
    Benchmark of linear solvers

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Dec 2022
  * License   : MIT License
=###############################################################################

import Printf: @sprintf

function calc_lift_drag(body, b, ar, Vinf, magVinf, rho; verbose=true, lbl="")

    CLexp = 0.238
    CDexp = 0.005

    str = ""

    if verbose
        str *= @sprintf "| %15.15s | %-7s | %-7s |\n" "Solver" "CL" "CD"
        str *= "| --------------: | :-----: | :-----: |\n"
    end

    # Calculate surface velocity induced by the body on itself
    Us = pnl.calcfield_U(body, body)

    # NOTE: Since the boundary integral equation of the potential flow has a
    #       discontinuity at the boundary, we need to add the gradient of the
    #       doublet strength to get an accurate surface velocity

    # Calculate surface velocity U_∇μ due to the gradient of the doublet strength
    UDeltaGamma = pnl.calcfield_Ugradmu(body)

    # Add both velocities together
    pnl.addfields(body, "Ugradmu", "U")

    # Calculate pressure coefficient
    Cps = pnl.calcfield_Cp(body, magVinf)

    # Calculate the force of each panel
    Fs = pnl.calcfield_F(body, magVinf, rho)
    
    # Calculate total force of the vehicle decomposed as lfit, drag, and sideslip
    Dhat = Vinf/pnl.norm(Vinf)        # Drag direction
    Shat = [0, 1, 0]              # Span direction
    Lhat = pnl.cross(Dhat, Shat)      # Lift direction

    LDS = pnl.calcfield_LDS(body, Lhat, Dhat)

    L = LDS[:, 1]
    D = LDS[:, 2]

    # Force coefficients
    nondim = 0.5*rho*magVinf^2*b^2/ar   # Normalization factor
    CL = sign(pnl.dot(L, Lhat)) * pnl.norm(L) / nondim
    CD = sign(pnl.dot(D, Dhat)) * pnl.norm(D) / nondim
    err = abs(CL-CLexp)/CLexp

    if verbose
        str *= @sprintf "| %15.15s | %-7.4f | %-7.4f |\n" lbl CL CD
        str *= @sprintf "| %15.15s | %-7s | %-7s |\n" "Experimental" CLexp CDexp
        str *= @sprintf "\n\tCL Error:\t%4.3g﹪\n" err*100
    end

    return CL, CD, str

end


# ----- Linear solver test: backslash \ operator
t = @elapsed pnl.solve(body, Uinfs, Das, Dbs; solver=pnl.solve_backslash!)

CL, CD, str = calc_lift_drag(body, b, ar, Vinf, magVinf, rho; lbl="Backslash")

str *= @sprintf "\tRun time:\t%4.2f seconds\n" t

open(joinpath(outdata_path, run_name*"-backslash.md"), "w") do f
    println(f, str)
end

println(str)


# ----- Linear solver test: LU decomposition + div
t = @elapsed pnl.solve(body, Uinfs, Das, Dbs; solver=pnl.solve_ludiv!)

CL, CD, str = calc_lift_drag(body, b, ar, Vinf, magVinf, rho; lbl="LUdiv")

str *= @sprintf "\tRun time:\t%4.2f seconds\n" t

open(joinpath(outdata_path, run_name*"-ludiv.md"), "w") do f
    println(f, str)
end

println(str)


# ----- Linear solver test: GMRES tol=1e-8
stats = []                      # Stats of GMRES get stored here

t = @elapsed pnl.solve(body, Uinfs, Das, Dbs;
                        solver = pnl.solve_gmres!,
                        solver_optargs = (atol=1e-8, rtol=1e-8, out=stats))

CL, CD, str = calc_lift_drag(body, b, ar, Vinf, magVinf, rho; lbl="GMRES tol=1e-8")

str *= replace("\n$(stats[1])", "\n"=>"\n\t")
str *= @sprintf "\n\tRun time:\t%4.2f seconds\n" t

open(joinpath(outdata_path, run_name*"-gmres8.md"), "w") do f
    println(f, str)
end

println(str)


# ----- Linear solver test: GMRES tol=1e-2
stats = []                      # Stats of GMRES get stored here

t = @elapsed pnl.solve(body, Uinfs, Das, Dbs;
                        solver = pnl.solve_gmres!,
                        solver_optargs = (atol=1e-2, rtol=1e-2, out=stats))

CL, CD, str = calc_lift_drag(body, b, ar, Vinf, magVinf, rho; lbl="GMRES tol=1e-2")

str *= replace("\n$(stats[1])", "\n"=>"\n\t")
str *= @sprintf "\n\tRun time:\t%4.2f seconds\n" t

open(joinpath(outdata_path, run_name*"-gmres2.md"), "w") do f
    println(f, str)
end

println(str)



# ----- Linear solver test: GMRES tol=1e-1
stats = []                      # Stats of GMRES get stored here

t = @elapsed pnl.solve(body, Uinfs, Das, Dbs;
                        solver = pnl.solve_gmres!,
                        solver_optargs = (atol=1e-1, rtol=1e-1, out=stats))

CL, CD, str = calc_lift_drag(body, b, ar, Vinf, magVinf, rho; lbl="GMRES tol=1e-1")

str *= replace("\n$(stats[1])", "\n"=>"\n\t")
str *= @sprintf "\n\tRun time:\t%4.2f seconds\n" t

open(joinpath(outdata_path, run_name*"-gmres1.md"), "w") do f
    println(f, str)
end

println(str)
