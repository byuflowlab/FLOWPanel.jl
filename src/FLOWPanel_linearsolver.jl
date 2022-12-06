#=##############################################################################
# DESCRIPTION
    Methods to solve a linear system of equations.

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Dec 2022
  * License     : MIT License
=###############################################################################


################################################################################
# SOLVERS
################################################################################
"""
    solve_backslash!(y::AbstractVector, A::AbstractMatrix, b::AbstractVector)

Solves a linear system of equations of the form Ay = b using the \\ operator.
The solution is stored under `y`.
"""
function solve_backslash!(y::AbstractVector,
                                A::AbstractMatrix, b::AbstractVector)

    y .= A\b

    return y
end

"""
    solve_ludiv!(y::AbstractVector,
                    A::AbstractMatrix, b::AbstractVector; Alu=nothing)

Solves a linear system of equations of the form Ay = b using the LU
decomposition of `A` provided under `Alu` and `LinearAlgebra.ldiv!`. If `Alu`
is not provided, it will be automatically calculated using `LinearAlgebra.lu`.

This method is useful when the system needs to be solved multiple times for
different `b` vectors since `Alu` can be precomputed and re-used. We
recommend you use `calc_Alu` to compute `Alu` since it has overloaded to also
handle Dual and TrackedReal numbers. `solve_ludiv!` has also been overloaded
with ImplicitAD to efficiently differentiate the linear solver as needed.

The solution is stored under `y`.


```@example
import FLOWPanel as pnl

Alu = pnl.calc_lu(A)
solve_ludiv!(y, A, b; Alu=Alu)
```
"""
function solve_ludiv!(y::AbstractVector,
                        A::AbstractMatrix{T}, b::AbstractVector; Alu=nothing
                        ) where {T}


    if isnothing(Alu)           # Case: No LU decomposition provided

        # Calculate LU decomposition
        _Alu = calc_Alu(A)

        # Recursive step
        solve_ludiv!(y, A, b; Alu=_Alu)

    else                        # Case: LU decomposition provided

        # Solve the linear system using ldiv!
        y .= IAD.implicit_linear(A, b, LA.ldiv!, Alu)

    end

    return y
end




"""
    solve_gmres!(y::AbstractVector,
                    A::AbstractMatrix, b::AbstractVector; Avalue=nothing,
                    atol=1e-8, rtol=1e-8, optargs...)

Solves a linear system of equations of the form Ay = b through the generalized
minimal residual (GMRES) method, which is an iterative method in the Krylov
subspace.

This iterative method is more efficient than a direct method (`solve_backslack!`
or `solve_ludiv!`) when `A` is larger than 3000x3000 or so. Also, iterative
methods can trade off accuracy for speed by lowering the tolerance (`atol` and
`rtol`). Optional arguments `optargs...` will be passed to `Krylov.gmres`.

Differentiating through the solver will require extracting the primal values of
`A`, which can be provided through the argument `Avalue` (this is calculated
automatically if not already provided).

The solution is stored under `y`.

```@example
import FLOWPanel as pnl

Avalue = pnl.calc_Avalue(A)
solve_gmres!(y, A, b; Avalue=Avalue)
```
"""
function solve_gmres!(y::AbstractVector,
                        A::AbstractMatrix, b::AbstractVector; Avalue=nothing,
                        optargs...)

    if isnothing(Avalue)        # Case: No primal values provided

        # Extract primal values
        _Avalue = calc_Avalue(A)

        # Recursive step
        solve_gmres!(y, A, b; Avalue=_Avalue, optargs...)

    else                        # Case: Primal values provided

        # Define GMRES function with these optional arguments
        this_gmres(A, b) = gmres(A, b; optargs...)

        # Solve the linear system through GMRES
        y .= IAD.implicit_linear(A, b, this_gmres, Avalue)

    end

    return y
end


function gmres(A, b::AbstractVector; memory=50, history=false,
                                        atol=sqrt(eps()), rtol=sqrt(eps()),
                                        out=[],
                                        optargs...)

    x, stats = Krylov.gmres(A, b, memory=memory, history=history;
                                        atol=atol, rtol=rtol, optargs...)

    push!(out, stats)

    return x
end

function gmres(A, b::AbstractMatrix; optargs...)

    return hcat(
            collect(
                gmres(A, Vector(bcol); optargs...) for bcol in eachcol(b)
                )...
            )

end
#### END OF SOLVERS ############################################################





################################################################################
# PREPROCESSING OF LINEAR SYSTEM
################################################################################
"""
    calc_Alu!(Apivot::AbstractMatrix, A::AbstractMatrix) -> Alu

Returns the LU decomposition of `A` using `Apivot` as storage memory to pivot
leaving `A` unchanged.
"""
function calc_Alu!(Apivot, A::AbstractMatrix{T}) where {T}

    # Prepare pivot array
    calc_Avalue!(Apivot, A)

    # LU decomposition
    Alu = LA.lu!(Apivot)

    return Alu
end


"""
    calc_Alu!(A::AbstractMatrix) -> Alu

Returns the LU decomposition of `A`. If `A` does not carry Dual nor TrackedReal
numbers, computation will be done in-place using `A`; hence `A` should not be
reused for multiple solves or for implicit differentiation (use `calc_Alu(A)`
and `calc_Alu!(Apivot, A)` instead).
"""
function calc_Alu!(A::AbstractMatrix{T}) where {T}

    # Allocate memory for pivot
    if T<:FD.Dual || T<:RD.TrackedReal  # Automatic differentiation case

        Tprimal = T.parameters[T<:FD.Dual ? 2 : 1]
        Apivot = zeros(Tprimal, size(A))

        # LU decomposition
        Alu = calc_Alu!(Apivot, A)

    else
        # LU decomposition
        Alu =  LA.lu!(A)
    end

    # LU decomposition
    return calc_Alu!(Apivot, A)
end

"""
    calc_Alu(A::AbstractMatrix) -> Alu

Returns the LU decomposition of `A`.
"""
function calc_Alu(A::AbstractMatrix{T}) where {T}

    # Allocate memory for pivot
    if T<:FD.Dual || T<:RD.TrackedReal  # Automatic differentiation case

        Tprimal = T.parameters[T<:FD.Dual ? 2 : 1]
        Apivot = zeros(Tprimal, size(A))

    else
        Apivot = zeros(T, size(A))
    end

    # LU decomposition
    return calc_Alu!(Apivot, A)
end

"""
    calc_Avalue!(Avalue::AbstractMatrix, A::AbstractMatrix)

Copy the primal values of `A` into `Avalue`.
"""
function calc_Avalue!(Avalue, A::AbstractMatrix{T}) where {T}

    if T<:FD.Dual || T<:RD.TrackedReal  # Automatic differentiation case

        # Extract primal values of A
        value = T<:FD.Dual ? FD.value : RD.value
        map!(value, Avalue, A)

    else                                # Normal case
        # Deep copy A
        Avalue .= A
    end

    return Avalue
end

"""
    calc_Avalue(A::AbstractMatrix)

Return the primal values of matrix `A`, which is simply `A` if the elements
of `A` are not Dual nor TrackedReal numbers.
"""
function calc_Avalue(A::AbstractMatrix{T}) where {T}

    if T<:FD.Dual || T<:RD.TrackedReal  # Automatic differentiation case

        Tprimal = T.parameters[T<:FD.Dual ? 2 : 1]
        Avalue = zeros(Tprimal, size(A))
        calc_Avalue!(Avalue, A)

        return Avalue
    else                                # Normal case

        return A
    end
end
#### END OF PREPROCESSING ######################################################






################################################################################
# BOUNDARY CONDITIONS
################################################################################

"""
    calc_bc_noflowthrough!(RHS::AbstractVector,
                            Us::AbstractMatrix, normals::AbstractMatrix)

Given the velocity and normals at/of each control point, it calculates the
right-hand side of the linear system of equations that imposes the
no-flow-through condition, which is stored under `RHS`.
"""
function calc_bc_noflowthrough!(RHS::AbstractVector,
                                Us::AbstractMatrix, normals::AbstractMatrix)

    @assert size(Us)==size(normals) ""*
        "Invalid matrices `Us` and `normals`;"*
        " expected to have the same size, got $(size(Us)) and $(size(normals))."
    @assert length(RHS)==size(Us, 2)
        "Invalid vector `RHS`;"*
        " expected length $(size(Us, 2)), got $(length(RHS))."

    for (i, (U, normal)) in enumerate(zip(eachcol(Us), eachcol(normals)))
        RHS[i] = -dot(U, normal)
    end

    return RHS
end

"""
    calc_bc_noflowthrough(Us::AbstractMatrix, normals::AbstractMatrix) -> RHS

Given the velocity and normals at/of each control point, it calculates and
returns the right-hand side of the linear system of equations that imposes the
no-flow-through condition.
"""
function calc_bc_noflowthrough(Us::AbstractMatrix{T1},
                                normals::AbstractMatrix{T2}) where {T1, T2}
    T = promote_type(T1, T2)
    RHS = zeros(T, size(Us, 2))
    calc_bc_noflowthrough!(RHS, Us, normals)

    return RHS
end
#### END OF BOUNDARY CONDITIONS ################################################
