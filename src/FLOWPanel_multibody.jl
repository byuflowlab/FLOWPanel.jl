#=##############################################################################
# DESCRIPTION
    Definition of a multi-body type made out of multiple bodies with different
    element types.

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Sep 2022
  * License     : MIT License
=###############################################################################




################################################################################
# MULTIBODY TYPE
################################################################################
struct MultiBody{E, N, B<:AbstractBody, TF<:Number} <: AbstractBody{E, N, TF}

    # User inputs
    bodies::Array{B, 1}                       # Array of bodies
    names::Array{String, 1}                   # Name of each body

    # Properties
    nbodies::Int                              # Number of bodies
    nnodes::Int                               # Number of nodes
    nsheddings::Int                           # Number of shedding edges
    ncells::Int                               # Number of cells
    solved::Bool                              # Whether it has been solved
    Oaxis::Array{TF,2}                  # Coordinate system orientation
    O::Array{TF,1}                      # Coordinate system origin
    Cps::Vector{TF}                           # Pressure coefficient at each cell

    function MultiBody{E, N, B, TF}(    bodies, names;
                                    nbodies=length(bodies),
                                    nnodes=sum(body.nnodes for body in bodies; init=0),
                                    ncells=sum(body.ncells for body in bodies; init=0),
                                    nsheddings=_calc_nsheddings(bodies),
                                    solved=false,
                                    Oaxis=Array{TF,2}(1.0I, 3, 3), O=zeros(TF,3)
                                    ) where {E, N, B, TF}

            @assert length(bodies)==length(names) ""*
                "Found different number of bodies than names"

            @assert length(unique(names)) == length(names) ""*
                "Got repeated names: $(names)"

            return new(
                        bodies, names,
                        nbodies,
                        nnodes,
                        nsheddings,
                        ncells,
                        solved,
                        Oaxis, O,
                        zeros(TF, ncells)
                      )
    end
end

function MultiBody(bodies::Array{B, 1}, args...; optargs...) where {B<:Union{AbstractBody, AbstractLiftingBody}}
    # Create Union of all body types
    ellist = [typeof(body).parameters[1] for body in bodies]
    E = Union{ellist...}

    countEbodies = length(bodies) > 0 ? maximum(typeof(body).parameters[2] for body in bodies) : 0

    N = max(countEbodies, _count(E))

    TF = promote_type([eltype(body.grid._nodes) for body in bodies]...)

    return MultiBody{E, N, B, TF}(bodies, args...; optargs...)
end

function MultiBody(bodies::Array{<:Union{AbstractBody, AbstractLiftingBody}, 1}; optargs...)
    names = ["Body$(bi)" for (bi, body) in enumerate(bodies)]
    return MultiBody(bodies, names; optargs...)
end

# Empty initializer
MultiBody() = MultiBody(AbstractBody[], String[])

"Returns the requested body"
function get_body(self::MultiBody, bodyname::String)

    bodyfound = findall(x->x==bodyname, self.names)

    if length(bodyfound)==0
        error("$bodyname not found. Bodies available are $(self.names).")
    elseif length(bodyfound)!=1
        error("Logic error: Found more than one grid $bodyname.")
    end

    return get_body(self, bodyfound[1])
end

"Returns the requested body"
function get_body(self::MultiBody, bodyindex::Int)
    if bodyindex>self.nbodies
        error("Requested invalid body index $bodyindex; max is $(self.nbodies).")
    end

    return self.bodies[bodyindex]
end

function get_strength(self::MultiBody, i)

    counter = 0

    for body in self.bodies
        offset = get_nstrengths(body)

        if i>counter && i<=counter+offset
            return get_strength(body, i-counter)
        end

        counter += offset
    end
end

get_nstrengths(self::MultiBody) = sum(get_nstrengths(body) for body in self.bodies, init=0)

function set_strength(self::MultiBody, i::Int, val)

    counter = 0

    for body in self.bodies
        offset = get_nstrengths(body)

        if i>counter && i<=counter+offset
            return set_strength(body, i-counter, val)
        end

        counter += offset
    end
end

function check_solved(self::MultiBody)
    if check_field(self, "solved")
        return get_fieldval(self, "solved", 1, "system")
    else
        return false
    end
end

function save(multibody::MultiBody, filename::String, args...; optargs...)

    str = ""

    for (name, body) in zip(multibody.names, multibody.bodies)
        str *= save(body, filename*"_"*name, args...; optargs...)
    end

    return str
end

function rotatetranslate!(multibody::MultiBody,
                            M::AbstractMatrix, T::AbstractVector; optargs...)

    # Recursively update bodies
    for body in multibody.bodies
        rotatetranslate!(body, M, T; optargs...)
    end

    # Update position and coordinate system of the overall system
    multibody.O .+= T
    multibody.Oaxis .= M*multibody.Oaxis

    nothing
end

function set_coordinatesystem(multibody::MultiBody,
                                O::AbstractVector, Oaxis::AbstractMatrix;
                                optargs...)

    for body in multibody.bodies
        set_coordinatesystem(body, O, Oaxis; optargs...)
    end

    # Undo its current coordinate system
    rotatetranslate!(body, inv(body.Oaxis), -body.O; optargs...)

    # Set new coordinate system
    rotatetranslate!(body, Oaxis, O; optargs...)

    return nothing
end



################################################################################
# PURE VORTEX RING SOLVER
################################################################################
function solve(self::MultiBody{VortexRing, 1},
                Uinfs::AbstractMatrix{T1},
                Das::AbstractMatrix{T2}; optargs...) where {T1, T2}

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # Compute geometric matrix (left-hand-side influence matrix)
    Gdims = _get_Gdims(self)
    G = zeros(promote_type(T1, T2), Gdims)

    solve!(self, G, normals, CPs, Uinfs, Das; optargs...)
end

function solve!(self::MultiBody{VortexRing, 1},
                G::AbstractMatrix{T1},
                normals::AbstractMatrix{T2}, CPs::AbstractMatrix,
                Uinfs::AbstractMatrix{T3}, Das::AbstractMatrix;
                solver=solve_ludiv!, solver_optargs=(), optargs...
                ) where {T1, T2, T3}

    if size(Uinfs) != (3, self.ncells)
        error("Invalid Uinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Uinfs))")
    elseif size(Das) != (3, self.nsheddings+1)
        error("Invalid Das;"*
              " expected size (3, $(self.nsheddings+1)), got $(size(Das))")
    end

    # Compute geometric matrix (left-hand-side influence matrix)
    _G_U!(self, G, CPs, normals, Das; optargs...)

    # Calculate boundary conditions (right-hand side of system of equations)
    RHS = calc_bc_noflowthrough(Uinfs, normals)

    # Solve system of equations
    Gamma = zeros(promote_type(T1, T2, T3), size(G, 1))
    solver(Gamma, G, RHS; solver_optargs...)

    # Save solution
    _set_strength(self, Gamma)

    _solvedflag(self, true)
    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(self, "Da", "vector", collect(eachcol(Das)), "system")
    add_field(self, "Gamma", "scalar", Gamma, "cell")
end

function solve2!(self::MultiBody,
                Uinfs::Matrix{T1},
                solver::AbstractMatrixfulSolver{false};
                    update_G::Bool=true,
                    solver_optargs=(),
                    elprescribe="automatic",
                    optargs...
                    # elprescribe_index::Int=1, elprescribe_value=0,
                    # weight_gammat=0, weight_gammao=1
                ) where T1<:Real

    if size(Uinfs) != (3, self.ncells)
        error("Invalid Uinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Uinfs))")
    end

    # get CPs and normals
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # Compute geometric matrix (left-hand-side influence matrix)
    solver.G .= zero(eltype(solver.G))
    _G_U!(self, solver.G, CPs, normals; optargs...)

    # Calculate boundary conditions (right-hand side of system of equations)
    RHS = calc_bc_noflowthrough(Uinfs, normals)

    # Solve system of equations
    Gamma = zeros(T1, size(solver.G, 1))
    solve_matrix!(Gamma, solver.G, RHS, solver; solver_optargs...)

    # Save solution
    set_solution(self, Gamma, Gamma, Tuple{Int,Float64}[], Uinfs)
end




################################################################################
# PURE VORTEX-RING LEAST-SQUARE SOLVER
################################################################################
function solve(self::MultiBody{VortexRing, 2},
                Uinfs::AbstractMatrix{T1},
                Das::AbstractMatrix{T2};
                elprescribe="automatic",
                optargs...
                ) where {T1, T2}

    T = promote_type(T1, T2)

    # Determine prescribed elements
    _elprescribe = elprescribe=="automatic" ? calc_elprescribe(self) : elprescribe

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # Allocate solver memory
    (; Gamma, Gammals, G, Gred, tGred, gpuGred, Gls, RHS, RHSls) = allocate_solver(self, _elprescribe, T)

    # Solve the least-squares problem to calculate strengths
    solve!(self, Gamma, Gammals,
            G, Gred, tGred, gpuGred, Gls, RHS, RHSls,
            normals, CPs,
            Uinfs, Das; elprescribe=_elprescribe, optargs...)
end

function allocate_solver(self::AbstractBody, _elprescribe, T::Type)

    n = self.ncells
    npres = length(_elprescribe)

    Gamma = zeros(T, n)
    Gammals = zeros(T, n-npres)
    G = zeros(T, n, n)
    Gred = zeros(T, n, n-npres)
    tGred = zeros(T, n-npres, n)
    gpuGred = zeros(T, n, n-npres)
    Gls = zeros(T, n-npres, n-npres)
    RHS = zeros(T, n)
    RHSls = zeros(T, n-npres)

    return (; Gamma, Gammals, G, Gred, tGred, gpuGred, Gls, RHS, RHSls)
end

function solve!(self::MultiBody{VortexRing, 2},
                Gamma, Gammals,
                G::AbstractMatrix, Gred, tGred, gpuGred, Gls, RHS, RHSls,
                normals::AbstractMatrix, CPs::AbstractMatrix,
                Uinfs::AbstractMatrix, Das::AbstractMatrix;
                solver=solve_ludiv!, solver_optargs=(),
                elprescribe=Tuple{Int, Float64}[],
                optargs...
                )

    n = self.ncells
    npres = length(elprescribe)

    if size(Uinfs) != (3, self.ncells)
        error("Invalid Uinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Uinfs))")
    elseif size(Das) != (3, self.nsheddings+1)
        error("Invalid Das;"*
              " expected size (3, $(self.nsheddings+1)), got $(size(Das))")
    end

    @assert length(Gammals)==n-npres ""*
        "Invalid Gammals; expected length $(n-npres), got $(length(Gammals))"
    @assert length(Gamma)==n ""*
        "Invalid Gamma; expected length $(n), got $(length(Gamma))"

    # Compute geometric matrix (left-hand-side influence matrix) and boundary
    # conditions (right-hand-side) converted into a least-squares problem
    _G_U_RHS!(self, G, Gred, tGred, gpuGred, Gls, RHS, RHSls,
                        Uinfs, CPs, normals, Das, elprescribe; optargs...)

    # Solve system of equations
    solver(Gammals, Gls, RHSls; solver_optargs...)

    # Save solution
    set_solution(self, Gamma, Gammals, elprescribe, Uinfs)
end

function solve2!(self::MultiBody,
                Uinfs::AbstractMatrix{T1},
                solver::AbstractMatrixfulSolver{true},
                Das::AbstractMatrix{T2};
                    update_G::Bool=true,
                    solver_optargs=(),
                    elprescribe="automatic",
                    optargs...
                    # elprescribe_index::Int=1, elprescribe_value=0,
                    # weight_gammat=0, weight_gammao=1
                ) where {T1, T2}

    if size(Uinfs) != (3, self.ncells)
        error("Invalid Uinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Uinfs))")
    elseif size(Das) != (3, self.nsheddings+1)
        error("Invalid Das;"*
              " expected size (3, $(self.nsheddings+1)), got $(size(Das))")
    end

    # Determine prescribed elements
    _elprescribe = elprescribe=="automatic" ? calc_elprescribe(self) : elprescribe

    # Allocate solver memory
    T = promote_type(T1, T2)
    (; Gamma, Gammals, G, Gred, tGred, gpuGred, Gls, RHS, RHSls) = allocate_solver(self, _elprescribe, T)

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # Compute geometric matrix (left-hand-side influence matrix) and boundary
    # conditions (right-hand-side) converted into a least-squares problem
    _G_U_RHS!(self, G, Gred, tGred, gpuGred, Gls, RHS, RHSls,
                        Uinfs, CPs, normals, Das, _elprescribe; optargs...)

    # Solve system of equations
    solve_matrix!(Gammals, Gls, RHSls, solver; solver_optargs...)

    # Save solution
    set_solution(self, Gamma, Gammals, _elprescribe, Uinfs)
end

function numtype(self::MultiBody)
    Tlist = [numtype(body) for body in self.bodies]
    return promote_type(Tlist...)
end

"Assumes the solved strengths are all the first columns of the strength arrays."
function set_solution(self::AbstractBody,
                        Gamma, Gammals, elprescribe,
                        Uinfs)

    if length(elprescribe)==0
        Gamma .= Gammals
    else
        prev_eli = 0
        for (i, (eli, elval)) in enumerate(elprescribe)

            @show size(Gamma) size(Gammals) (prev_eli+1):(eli-1) (prev_eli+2-i):(eli-i)
            Gamma[(prev_eli+1):(eli-1)] .= view(Gammals, (prev_eli+2-i):(eli-i))

            Gamma[eli] = elval

            if i==length(elprescribe) && eli!=length(Gamma)
                Gamma[eli+1:end] .= view(Gammals, (eli-i+1):length(Gammals))
            end

            prev_eli = eli
        end
    end
    _set_strength(self, Gamma)

    _solvedflag(self, true)
    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field_D(self)
    add_field(self, "scalar", Gamma, "cell")
end

function add_field_D(self::MultiBody)
    for body in self.bodies
        add_field_D(body)
    end
end

function add_field_D(self::AbstractBody)
    return nothing
end

function add_field_D(self::RigidWakeBody)
    add_field(self, "Da", "vector", collect(eachcol(self.Das)), "system")
end

function set_solution(self::MultiBody,
                        Gamma, Gammals, elprescribe,
                        Uinfs)

    if length(elprescribe)==0
        Gamma .= Gammals
    else
        prev_eli = 0
        for (i, (eli, elval)) in enumerate(elprescribe)

            Gamma[(prev_eli+1):(eli-1)] .= view(Gammals, (prev_eli+2-i):(eli-i))

            Gamma[eli] = elval

            if i==length(elprescribe) && eli!=length(Gamma)
                Gamma[eli+1:end] .= view(Gammals, (eli-i+1):length(Gammals))
            end

            prev_eli = eli
        end
    end
    _set_strength(self, Gamma)
    
    _solvedflag(self, true)
    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field_D(self)
    add_field(self, "scalar", Gamma, "cell")
end

function _G_U_RHS(self::MultiBody{VortexRing, 1}, args...; optargs...)
    return _G_U_RHS_leastsquares(self, args...; optargs...)
end

function _G_U_RHS!(self::MultiBody{VortexRing, 1}, args...; optargs...)
    return _G_U_RHS_leastsquares!(self, args...; optargs...)
end


################################################################################
# COMMON SOLVER FUNCTIONS
################################################################################
function _G_U!(self::MultiBody,
                    G::Arr1, CPs::Arr2, normals::Arr3;
                    optargs...
               ) where{ T1, Arr1<:AbstractArray{T1, 2},
                        T2, Arr2<:AbstractArray{T2, 2},
                        T3, Arr3<:AbstractArray{T3, 2}}

    ncells = 0
    # nsheddings = 0

    for body in self.bodies

        # if body isa AbstractLiftingBody || body isa MultiBody

            _G_U!(body, view(G, :, (1:body.ncells) .+ ncells), CPs, normals; optargs...)

        #     nsheddings += body.nsheddings

        # else
        #     _G_U!(body, view(G, :, (1:body.ncells) .+ ncells), CPs, normals; optargs...)
        # end

        ncells += body.ncells

    end

end

function _G_Uvortexring!(self::MultiBody, G, CPs, normals; optargs...)

    ncells = 0
    # nsheddings = 0

    for body in self.bodies

        # if body isa AbstractLiftingBody || body isa MultiBody

            _G_Uvortexring!(body, view(G, :, (1:body.ncells) .+ ncells), CPs, normals; optargs...)

        #     nsheddings += body.nsheddings

        # else
        #     _G_Uvortexring!(body, view(G, :, (1:body.ncells) .+ ncells), CPs, normals; optargs...)
        # end

        ncells += body.ncells

    end

end

"""
    `calc_elprescribe(multibody::MultiBody; indices=[1], values=[0.0])`

Automatically calculate the elements to prescribe in a MultiBody, recognizing
all bodies that are watertight and prescribing the strength of elements of
indices `indices` in each body the values `values`. Returns array `elprescribe`
that is used by `solve!`.
"""
function calc_elprescribe(self::MultiBody;
                            indices::AbstractArray{Int}=[1],
                            values::AbstractArray{<:Number}=[0.0],
                            elprescribe=Tuple{Int, Float64}[], ncells0=0)

    @assert length(indices)==length(values) ""*
        "Got unequal number of indices ($(length(indices)) and "*
        "values ($(length(values)))"

    ncells = ncells0

    for body in self.bodies

        # Recursive case
        if body isa MultiBody

            calc_elprescribe(body; elprescribe=elpresecribe, ncells0=ncells)

        # Base case: watertight body
        elseif body.watertight

            for (ind, val) in zip(indices, values)
                push!(elprescribe, (ncells+ind, val))
            end

        # Base case: non-watertight body
        else
            nothing
        end

        ncells += body.ncells
    end

    return elprescribe
end

function _set_strength(self::MultiBody{VortexRing}, Gamma)

    ncells = 0

    for body in self.bodies

        strength = view(Gamma, (1:body.ncells) .+ ncells)

        if body isa MultiBody
            _set_strength(body, strength)
        else
            body.strength[:, 1] .= strength
        end

        ncells += body.ncells
    end
end


#### INTERNAL FUNCTIONS ########################################################
function _Uind!(self::MultiBody, targets, out, backend::DirectBackend, args...; optargs...)
    for body in self.bodies
        _Uind!(body, targets, out, backend, args...; optargs...)
    end
end

function _phi!(self::MultiBody, targets, out, backend::DirectBackend, args...; optargs...)
    for body in self.bodies
        _phi!(body, targets, out, backend, args...; optargs...)
    end
end

function _get_Gdims(self::MultiBody)
    n, m = 0, 0
    for body in self.bodies
        Gdims = _get_Gdims(body)
        n += Gdims[1]
        m += Gdims[2]
    end

    return (n, m)
end

function _calc_nsheddings(bodies)
    nsheddings = 0

    for body in bodies
        if body isa MultiBody
            nsheddings += _calc_nsheddings(body)
        elseif body isa AbstractLiftingBody
            nsheddings += body.nsheddings
        end
    end

    return nsheddings
end

function _calc_controlpoints(self::MultiBody, args...; optargs...)
    controlpoints = zeros(3, self.ncells)
    _calc_controlpoints!(self, controlpoints, args...; optargs...)
    return controlpoints
end
function _calc_controlpoints!(self::MultiBody, controlpoints, normals; optargs...)
    ncells = 0
    for body in self.bodies

        # Determine direction of control point offset based on Body's default
        # offset (if offset was requested)
        if body isa MultiBody
            this_optargs = optargs
        else
            this_optargs = ((key, key==:off ? sign(body.CPoffset)*val : val) for (key, val) in optargs)
        end

        rng = (1:body.ncells) .+ ncells
        _calc_controlpoints!(body, view(controlpoints, :, rng), view(normals, :, rng); this_optargs...)

        ncells += body.ncells
    end
end

function _calc_tangents(self::MultiBody)
    tangents = zeros(3, self.ncells)
    _calc_tangents!(self, tangents)
    return tangents
end
function _calc_tangents!(self::MultiBody, tangents)
    ncells = 0
    for body in self.bodies
        _calc_tangents!(body, view(tangents, :, (1:body.ncells) .+ ncells))
        ncells += body.ncells
    end
end

function _calc_obliques(self::MultiBody)
    obliques = zeros(3, self.ncells)
    _calc_obliques!(self, obliques)
    return obliques
end
function _calc_obliques!(self::MultiBody, obliques)
    ncells = 0
    for body in self.bodies
        _calc_obliques!(body, view(obliques, :, (1:body.ncells) .+ ncells))
        ncells += body.ncells
    end
end

function _calc_areas(self::MultiBody)
    areas = zeros(self.ncells)
    _calc_areas!(self, areas)
    return areas
end
function _calc_areas!(self::MultiBody, areas)
    ncells = 0
    for body in self.bodies
        _calc_areas!(body, view(areas, (1:body.ncells) .+ ncells))
        ncells += body.ncells
    end
end

function _solvedflag(self::MultiBody, val::Bool)

    # Remove all existing fields
    for field in Iterators.reverse(self.fields)
        remove_field(self, field)
    end

    # Add solved flag
    add_field(self, "solved", "scalar", [val], "system")

    for body in self.bodies
        _solvedflag(body, val)
    end
end
#### END OF MULTIBODY ##########################################################
