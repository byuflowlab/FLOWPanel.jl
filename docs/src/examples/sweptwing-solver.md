# Solver Benchmark

The problem of solving the panel strengths that satisfy the
no-flow-through condition poses a linear system of equation of the form

```math
\begin{align*}
        A y = b
,\end{align*}
```

where $A$ is the matrix containing the geometry of the panels and wake,
$y$ is the vector of panel strengths, and $b$ is the vector of boundary
conditions.
This is trivially solved as

```math
\begin{align*}
        y = A^{-1}b
,\end{align*}
```

however, depending on the size of $A$ (which depends on the number of
panels) it can become inefficient or even unfeasible to explicitely
calculate the inverse of $A$.
Multiple linear solvers are available in FLOWPanel that avoid
explicitely inverting $A$, which are described and benchmarked as
follows.

> The complete code is available at [examples/sweptwing_solverbenchmark.jl](https://github.com/byuflowlab/FLOWPanel.jl/blob/master/examples/sweptwing_solverbenchmark.jl) but you should also be able to copy and paste these lines after running the first section of this example.


## Backslash operator `\`

Most programming languages implement an operator `\` that directly
calculates the matrix-vector product $A^{-1}b$.
This is more efficient than directly inverting $A$ and then multiplying
by $b$, without loosing any accuracy.
This linear solver is available under this function:
```@docs
FLOWPanel.solve_backslash!
```
and is used as follows:

```julia
import Printf: @sprintf

function calc_lift_drag(body, b, ar, Vinf, magVinf, rho; verbose=true, lbl="")

    CLexp = 0.238
    CDexp = 0.005

    str = ""

    if verbose
        str *= @sprintf "| %15.15s | %-7s | %-7s |\n" "Solver" "CL" "CD"
        str *= "| --------------: | :-----: | :-----: |\n"
    end

    # Calculate velocity away from the body
    Us = pnl.calcfield_U(body, body; fieldname="Uoff",
                            offset=0.02, characteristiclength=(args...)->b/ar)

    # Calculate surface velocity U_∇μ due to the gradient of the doublet strength
    UDeltaGamma = pnl.calcfield_Ugradmu(body)

    # Add both velocities together
    pnl.addfields(body, "Ugradmu", "Uoff")

    # Calculate pressure coeffiecient
    Cps = pnl.calcfield_Cp(body, magVinf; U_fieldname="Uoff")

    # Calculate the force of each panel
    Fs = pnl.calcfield_F(body, magVinf, rho; U_fieldname="Uoff")
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
```


|          Solver | CL      | CD      |
| --------------: | :-----: | :-----: |
|       Backslash | 0.2335  | 0.0132  |
|    Experimental | 0.238   | 0.005   |

	CL Error:	1.89﹪
	Run time:	0.44 seconds

## LU decomposition

Pre-calculating and re-using the
[LU decomposition](https://en.wikipedia.org/wiki/LU_decomposition) of $A$
is advantageous when the linear system needs to be solved for multiple
boundary conditions.

```@docs
FLOWPanel.solve_ludiv!
```

Running the solver:

```julia
t = @elapsed pnl.solve(body, Uinfs, Das, Dbs; solver=pnl.solve_ludiv!)

CL, CD, str = calc_lift_drag(body, b, ar, Vinf, magVinf, rho; lbl="LUdiv")
```
|          Solver | CL      | CD      |
| --------------: | :-----: | :-----: |
|           LUdiv | 0.2335  | 0.0132  |
|    Experimental | 0.238   | 0.005   |

	CL Error:	1.89﹪
	Run time:	0.44 seconds

## Iterative Krylov Solver

Iterative Krylov subspace solvers converge to the right solution
rather than directly solving the system of equations.
This allows the user to trade off accuracy for computational speed by
tuning the tolerance of the solver.

The [generalized minimal residual](https://en.wikipedia.org/wiki/Generalized_minimal_residual_method)
(GMRES) method provided by
[Krylov.jl](https://juliasmoothoptimizers.github.io/Krylov.jl/stable/solvers/unsymmetric/#GMRES)
is available in FLOWPanel.

```@docs
FLOWPanel.solve_gmres!
```

Running the solver with tolerance $10^{-8}$:
```julia
stats = []                      # Stats of GMRES get stored here

t = @elapsed pnl.solve(body, Uinfs, Das, Dbs;
                        solver = pnl.solve_gmres!,
                        solver_optargs = (atol=1e-8, rtol=1e-8, out=stats))

CL, CD, str = calc_lift_drag(body, b, ar, Vinf, magVinf, rho; lbl="GMRES tol=1e-8")
```
|          Solver | CL      | CD      |
| --------------: | :-----: | :-----: |
|  GMRES tol=1e-8 | 0.2335  | 0.0132  |
|    Experimental | 0.238   | 0.005   |

	CL Error:	1.89﹪

	Simple stats
	 niter: 286
	 solved: true
	 inconsistent: false
	 residuals: []
	 Aresiduals: []
	 κ₂(A): []
	 status: solution good enough given atol and rtol
	
	Run time:	0.71 seconds

Running the solver with tolerance $10^{-2}$:
```julia
stats = []                      # Stats of GMRES get stored here

t = @elapsed pnl.solve(body, Uinfs, Das, Dbs;
                        solver = pnl.solve_gmres!,
                        solver_optargs = (atol=1e-2, rtol=1e-2, out=stats))

CL, CD, str = calc_lift_drag(body, b, ar, Vinf, magVinf, rho; lbl="GMRES tol=1e-2")
```
|          Solver | CL      | CD      |
| --------------: | :-----: | :-----: |
|  GMRES tol=1e-2 | 0.2331  | 0.0129  |
|    Experimental | 0.238   | 0.005   |

	CL Error:	2.06﹪

	Simple stats
	 niter: 88
	 solved: true
	 inconsistent: false
	 residuals: []
	 Aresiduals: []
	 κ₂(A): []
	 status: solution good enough given atol and rtol
	
	Run time:	0.39 seconds

Running the solver with tolerance $10^{-1}$:
```julia
stats = []                      # Stats of GMRES get stored here

t = @elapsed pnl.solve(body, Uinfs, Das, Dbs;
                        solver = pnl.solve_gmres!,
                        solver_optargs = (atol=1e-1, rtol=1e-1, out=stats))

CL, CD, str = calc_lift_drag(body, b, ar, Vinf, magVinf, rho; lbl="GMRES tol=1e-1")
```
|          Solver | CL      | CD      |
| --------------: | :-----: | :-----: |
|  GMRES tol=1e-1 | 0.2625  | 0.0133  |
|    Experimental | 0.238   | 0.005   |

	CL Error:	10.3﹪

	Simple stats
	 niter: 25
	 solved: true
	 inconsistent: false
	 residuals: []
	 Aresiduals: []
	 κ₂(A): []
	 status: solution good enough given atol and rtol
	
	Run time:	0.29 seconds

