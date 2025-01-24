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
struct MultiBody{E, N, B<:Union{AbstractBody, AbstractLiftingBody}} <: AbstractBody{E, N}

    # User inputs
    bodies::Array{B, 1}                       # Array of bodies
    names::Array{String, 1}                   # Name of each body

    # Properties
    nbodies::Int                              # Number of bodies
    nnodes::Int                               # Number of nodes
    ncells::Int                               # Number of cells
    nsheddings::Int                           # Number of shedding edges
    fields::Array{String, 1}                  # Available fields (solutions)
    Oaxis::Array{<:Number,2}                  # Coordinate system orientation
    O::Array{<:Number,1}                      # Coordinate system origin

    function MultiBody{E, N, B}(    bodies, names;
                                    nbodies=length(bodies),
                                    nnodes=sum(body.nnodes for body in bodies; init=0),
                                    ncells=sum(body.ncells for body in bodies; init=0),
                                    nsheddings=_calc_nsheddings(bodies),
                                    fields=Array{String,1}(),
                                    Oaxis=Array(1.0I, 3, 3), O=zeros(3)
                                    ) where {E, N, B}

            @assert length(bodies)==length(names) ""*
                "Found different number of bodies than names"

            @assert length(unique(names)) == length(names) ""*
                "Got repeated names: $(names)"

            return new(
                        bodies, names,
                        nbodies,
                        nnodes,
                        ncells,
                        nsheddings,
                        fields,
                        Oaxis, O,
                      )
    end
end

function MultiBody(bodies::Array{B, 1}, args...; optargs...) where {B<:Union{AbstractBody, AbstractLiftingBody}}
    # Create Union of all body types
    ellist = [typeof(body).parameters[1] for body in bodies]
    E = Union{ellist...}

    countEbodies = length(bodies) > 0 ? maximum(typeof(body).parameters[2] for body in bodies) : 0

    N = max(countEbodies, _count(E))

    return MultiBody{E, N, B}(bodies, args...; optargs...)
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

# """
#   `get_ndivscells(body::MultiBody)`
#
# Returns a tuple with the number of cells in each parametric dimension
# """
# get_ndivscells(self::MultiBody) = Tuple(sum(collect(get_ndivscells(body)) for body in self.bodies))

# """
#   `get_ndivnodes(body::MultiBody)`
#
# Returns a tuple with the number of nodes in each parametric dimension
# """
# get_ndivsnodes(self::MultiBody) = Tuple(sum(collect(get_ndivsnodes(body)) for body in self.bodies))

get_ndivscells(self::MultiBody) = error("`get_ndivscells(...)` not implemented yet.")
get_ndivnodes(self::MultiBody) = error("`get_ndivnodes(...)` not implemented yet.")
get_cart2lin_cells(self::MultiBody) = error("`get_cart2lin_cells(...)` not implemented yet.")
get_cart2lin_nodes(self::MultiBody) = error("`get_cart2lin_nodes(...)` not implemented yet.")

function add_field(self::MultiBody{E, N, B}, field_name::String,
                    field_type::String, field_data, entry_type::String;
                    raise_warn=false) where {E, N, B<:Union{AbstractBody, AbstractLiftingBody}}

    if !(field_name in self.fields)
        push!(self.fields, field_name)
    end

    counter = 0
    prop = entry_type=="node" ? :nnodes : :ncells

    for body in self.bodies

        offset = getproperty(body, prop)

        if entry_type=="system"
            data_slice = field_data
        else
            data_slice = view(collect(field_data), (1:offset) .+ counter)
        end

        add_field(body, field_name, field_type, data_slice, entry_type;
                                                          raise_warn=raise_warn)
        counter += offset
    end

end

function remove_field(self::MultiBody, field_name)
    if check_field(self, field_name)

        i = findfirst(name->name==field_name, self.fields)

        splice!(self.fields, i)

        for body in self.bodies
            remove_field(body, field_name)
        end
    end
end


function get_field(self::MultiBody, field_name::String)

    # Collect all fields
    fields = [get_field(body, field_name) for body in self.bodies]

    # Concatenate field data
    if fields[1]["entry_type"]=="system"
        field_data = mean([field["field_data"] for field in fields])
    else
        field_data = vcat([field["field_data"] for field in fields]...)
    end

    # Create a new field
    field = Dict( "field_name" => field_name,
                  "field_type" => fields[1]["field_type"],
                  "entry_type" => fields[1]["entry_type"],
                  "field_data" => field_data)

    return field
end

function get_fieldval(self::MultiBody, field_name::String, i::Int,
                        entry_type::String; _check::Bool=true)

    prop = entry_type=="node"   ? :nnodes :
           entry_type=="cell"   ? :ncells :
           entry_type=="system" ? :ncells :
           error("Invalid entry type $(entry_type)")

    if entry_type!="system" && (i<=0 || i>getproperty(self, prop))
        error("Invalid index $(i) requested. Valid range is [1, $(getproperty(self, prop))]")
    end

    counter = 0

    for body in self.bodies
        offset = getproperty(body, prop)

        if i>counter && i<=counter+offset
            if body isa MultiBody
                return get_fieldval(body, field_name, i-counter, entry_type; _check=_check)
            else
                return get_fieldval(body, field_name, i-counter; _check=_check)
            end
        end

        counter += offset
    end

    error("Logic error!")
end

function get_fieldval(self::MultiBody, field_name::String, i::Int; optargs...)

    error("Invalid method for MultiBody."*
          " Try `get_fieldval(multibody, field_name, i, entry_type)` instead")

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
                Das::AbstractMatrix{T2},
                Dbs::AbstractMatrix{T3}; optargs...) where {T1, T2, T3}

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # Compute geometric matrix (left-hand-side influence matrix)
    Gdims = _get_Gdims(self)
    G = zeros(promote_type(T1, T2, T3), Gdims)

    solve!(self, G, normals, CPs, Uinfs, Das, Dbs; optargs...)
end

function solve!(self::MultiBody{VortexRing, 1},
                G::AbstractMatrix{T1},
                normals::AbstractMatrix{T2}, CPs::AbstractMatrix,
                Uinfs::AbstractMatrix{T3}, Das::AbstractMatrix, Dbs::AbstractMatrix;
                solver=solve_ludiv!, solver_optargs=(), optargs...
                ) where {T1, T2, T3}

    if size(Uinfs) != (3, self.ncells)
        error("Invalid Uinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Uinfs))")
    elseif size(Das) != (3, self.nsheddings)
        error("Invalid Das;"*
              " expected size (3, $(self.nsheddings)), got $(size(Das))")
    elseif size(Dbs) != (3, self.nsheddings)
        error("Invalid Dbs;"*
              " expected size (3, $(self.nsheddings)), got $(size(Dbs))")
    end

    # Compute geometric matrix (left-hand-side influence matrix)
    _G_U!(self, G, CPs, normals, Das, Dbs; optargs...)

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
    add_field(self, "Db", "vector", collect(eachcol(Dbs)), "system")
    add_field(self, "Gamma", "scalar", Gamma, "cell")
end





################################################################################
# PURE VORTEX-RING LEAST-SQUARE SOLVER
################################################################################
function solve(self::MultiBody{VortexRing, 2},
                Uinfs::AbstractMatrix{T1},
                Das::AbstractMatrix{T2},
                Dbs::AbstractMatrix{T3};
                elprescribe="automatic",
                optargs...
                ) where {T1, T2, T3}

    T = promote_type(T1, T2, T3)

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
            Uinfs, Das, Dbs; elprescribe=_elprescribe, optargs...)
end

function allocate_solver(self::Union{MultiBody{VortexRing, 2}, RigidWakeBody{VortexRing, 2}},
                                                                        _elprescribe, T::Type)

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
                Uinfs::AbstractMatrix, Das::AbstractMatrix, Dbs::AbstractMatrix;
                solver=solve_ludiv!, solver_optargs=(),
                elprescribe=Tuple{Int, Float64}[],
                optargs...
                )

    n = self.ncells
    npres = length(elprescribe)

    if size(Uinfs) != (3, self.ncells)
        error("Invalid Uinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Uinfs))")
    elseif size(Das) != (3, self.nsheddings)
        error("Invalid Das;"*
              " expected size (3, $(self.nsheddings)), got $(size(Das))")
    elseif size(Dbs) != (3, self.nsheddings)
        error("Invalid Dbs;"*
              " expected size (3, $(self.nsheddings)), got $(size(Dbs))")
    end

    @assert length(Gammals)==n-npres ""*
        "Invalid Gammals; expected length $(n-npres), got $(length(Gammals))"
    @assert length(Gamma)==n ""*
        "Invalid Gamma; expected length $(n), got $(length(Gamma))"

    # Compute geometric matrix (left-hand-side influence matrix) and boundary
    # conditions (right-hand-side) converted into a least-squares problem
    _G_U_RHS!(self, G, Gred, tGred, gpuGred, Gls, RHS, RHSls,
                        Uinfs, CPs, normals, Das, Dbs, elprescribe; optargs...)

    # Solve system of equations
    solver(Gammals, Gls, RHSls; solver_optargs...)

    # Save solution
    set_solution(self, Gamma, Gammals, elprescribe, Uinfs, Das, Dbs)
end

function set_solution(self::MultiBody{VortexRing, 2},
                        Gamma, Gammals, elprescribe,
                        Uinfs, Das, Dbs)

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
    add_field(self, "Da", "vector", collect(eachcol(Das)), "system")
    add_field(self, "Db", "vector", collect(eachcol(Dbs)), "system")
    add_field(self, "Gamma", "scalar", Gamma, "cell")
end

function _G_U_RHS(self::MultiBody{VortexRing, 2}, args...; optargs...)
    return _G_U_RHS_leastsquares(self, args...; optargs...)
end

function _G_U_RHS!(self::MultiBody{VortexRing, 2}, args...; optargs...)
    return _G_U_RHS_leastsquares!(self, args...; optargs...)
end


################################################################################
# COMMON SOLVER FUNCTIONS
################################################################################
function _G_U!(self::MultiBody,
                    G::Arr1, CPs::Arr2, normals::Arr3, Das, Dbs;
                    optargs...
               ) where{ T1, Arr1<:AbstractArray{T1, 2},
                        T2, Arr2<:AbstractArray{T2, 2},
                        T3, Arr3<:AbstractArray{T3, 2}}

    ncells = 0
    nsheddings = 0

    for body in self.bodies

        if body isa AbstractLiftingBody || body isa MultiBody

            _G_U!(body, view(G, :, (1:body.ncells) .+ ncells), CPs, normals,
                        view(Das, :, (1:body.nsheddings) .+ nsheddings),
                        view(Dbs, :, (1:body.nsheddings) .+ nsheddings); optargs...)

            nsheddings += body.nsheddings

        else
            _G_U!(body, view(G, :, (1:body.ncells) .+ ncells), CPs, normals; optargs...)
        end

        ncells += body.ncells

    end

end

function _G_Uvortexring!(self::MultiBody, G, CPs, normals, Das, Dbs; optargs...)

    ncells = 0
    nsheddings = 0

    for body in self.bodies

        if body isa AbstractLiftingBody || body isa MultiBody

            _G_Uvortexring!(body, view(G, :, (1:body.ncells) .+ ncells), CPs, normals,
                        view(Das, :, (1:body.nsheddings) .+ nsheddings),
                        view(Dbs, :, (1:body.nsheddings) .+ nsheddings); optargs...)

            nsheddings += body.nsheddings

        else
            _G_Uvortexring!(body, view(G, :, (1:body.ncells) .+ ncells), CPs, normals; optargs...)
        end

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
        elseif body isa RigidWakeBody{VortexRing, 2}

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
function _Uind!(self::MultiBody, args...; optargs...)
    for body in self.bodies
        _Uind!(body, args...; optargs...)
    end
end

function _phi!(self::MultiBody, args...; optargs...)
    for body in self.bodies
        _phi!(body, args...; optargs...)
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

function _calc_normals(self::MultiBody; optargs...)
    normals = zeros(3, self.ncells)
    _calc_normals!(self, normals; optargs...)
    return normals
end
function _calc_normals!(self::MultiBody, normals; optargs...)
    ncells = 0
    for body in self.bodies
        _calc_normals!(body, view(normals, :, (1:body.ncells) .+ ncells);
                                                                    optargs...)
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
