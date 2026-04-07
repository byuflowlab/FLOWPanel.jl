# ### ADDED
# function calc_normals!(bodies::Tuple; flipbyCPoffset=fill(false, length(bodies)))
#     for (i, body) in enumerate(bodies)
#         calc_normals!(body; flipbyCPoffset=flipbyCPoffset[i])
#     end

#     return nothing
# end

# ### ADDED
# function calc_controlpoints!(bodies::Tuple)
#     for body in bodies
#         calc_controlpoints!(body)
#     end

#     return nothing
# end

# ### ADDED
# function influence!(bodies::Tuple, targets::Tuple, backend::DirectBackend; optargs...)
#     for body in bodies
#         for target in targets
#             influence!(body, target, backend; optargs...)
#         end
#     end
# end

### ADDED
# function solve!(bodies::Tuple, solver::FGSSolver; backend = FastMultipoleBackend(
#         expansion_order=solver.expansion_order,
#         multipole_acceptance=solver.multipole_acceptance,
#         leaf_size=solver.leaf_size
#     ), optargs...) where TF

#     npanels = map(body -> body.ncells, bodies)
#     offsets = cumsum(vcat(0, npanels))

#     # save external velocity and potential
#     for (bi, body) in enumerate(bodies)
#         valid = body isa RigidWakeBody{<:Union{ConstantSource, ConstantDoublet, VortexRing}, 2}
#         @assert valid "All bodies must be RigidWakeBody{<:Union{ConstantSource, ConstantDoublet, VortexRing}, 2}"

#         r = offsets[bi]+1 : offsets[bi+1]
#         @views solver.Uext[:, r] .= body.velocity
#         @views solver.phi_ext[r] .= body.potential
#     end

#     # Store CP offsets
#     CPoffset_old = map(body -> body.CPoffset, bodies)

#     # Flip CP offset
#     for body in bodies
#         body.CPoffset = -abs(body.CPoffset)
#     end
    
#     # get normals and control points (inside)
#     calc_normals!(bodies)
#     calc_controlpoints!(bodies)

#     # Geometry + source strengths
#     for body in bodies
#         body.strength[:, 1] .= 0.0
#         for d in 1:3
#             @views body.strength[:, 1] .-= body.velocity[d, :] .* body.normals[d, :]
#         end
#         body.strength[:, 2] .= 0.0
#         body.potential .= 0
#     end

#     sigma = map(body -> copy(view(body.strength, :, 1)), bodies)
    
#     influence!(bodies, bodies, backend; scalar_potential=true, optargs...)
    
#     # run fgs solver
#     FastMultipole.solve!(bodies, solver.fgs;
#         max_iterations=solver.max_iterations,
#         inner_iterations=solver.inner_iterations,
#         tolerance=solver.tolerance,
#         rlx=solver.rlx,
#         derivatives_switches=FastMultipole.DerivativesSwitch(true, false, false, to_tuple(_unpack_fmm(bodies))),
#         reverse_pass=solver.reverse_pass,
#         verbose=solver.verbose,
#         final_update=false
#     )

#     # Save solution
#     # restore CPoffset, external velocity and potential
#     for (bi, b) in enumerate(bodies)
#         _solvedflag(b, true)
#         r = offsets[bi]+1 : offsets[bi+1]

#         b.strength[:, 1] .= sigma[bi]        
#         @views b.CPoffset = CPoffset_old[bi]
#         @views b.velocity .= solver.Uext[:, r]
#         @views b.potential .= solver.phi_ext[r]
#     end
# end

### ADDED
struct BackslashCoupled{TF}
    G::Matrix{TF}
    Glu::LU{TF, Matrix{TF}}
    rhs::Vector{TF}
    Uext::Array{TF, 3}
    phi_ext::Vector{TF}
    boundary::Symbol
end

function BackslashCoupled(bodies::Tuple{<:AbstractBody{<:Any,<:Any,TF}}) where TF
    ncs = sum(b -> b.ncells, bodies)

    G       = zeros(TF, ncs, ncs)
    rhs     = zeros(TF, ncs)
    Uext    = zeros(TF, 3, ncs)
    phi_ext = zeros(TF, ncs)

    # infer boundary
    if bodies[1] isa AbstractLiftingBody
        boundary = :dirichlet
    else
        boundary = :neumann
    end

    Glu = lu!(G)  # dummy init; will be overwritten on first update_G=true

    BackslashCoupled{TF}(G, Glu, rhs, Uext, phi_ext, boundary)
end

Backslash(bodies::Tuple) = BackslashCoupled(bodies)

### ADDED
function solve!(bodies::Tuple, solver::BackslashCoupled; backend=DirectBackend(), update_G::Bool=false, optargs...)

    # Sizes
    npanels = map(b -> b.ncells, bodies)
    offsets = cumsum(vcat(0, npanels))

    # save external fields
    for (bi, body) in enumerate(bodies)
        r = offsets[bi]+1 : offsets[bi+1]
        @views solver.Uext[:, r] .= body.velocity
        @views solver.phi_ext[r]  .= body.potential
    end

    # flip CP offset
    CPoffset_old = map(b -> b.CPoffset, bodies)
    for b in bodies
        b.CPoffset = -abs(b.CPoffset)
    end

    calc_normals!(bodies)
    calc_controlpoints!(bodies)

    # build RHS depending on boundary type
    if solver.boundary === :dirichlet
        for b in bodies
            b.strength[:,1] .= 0
            for d in 1:3
                @views b.strength[:,1] .-= b.velocity[d,:] .* b.normals[d,:]
            end
            b.strength[:,2] .= 0
            b.potential .= 0
        end

        influence!(bodies, bodies, backend; scalar_potential=true)
        for (bi, b) in enumerate(bodies)
            r = offsets[bi]+1 : offsets[bi+1]
            @views solver.rhs[r] .= -b.potential
        end

    elseif solver.boundary === :neumann
        for (bi, b) in enumerate(bodies)
            r = offsets[bi]+1 : offsets[bi+1]
            rhs_local = zeros(eltype(solver.rhs), length(r))
            for d in 1:3
                @views rhs_local .+= solver.Uext[d,r] .* b.normals[d,:]
            end
            @views solver.rhs[r] .= -rhs_local
        end

        # zero strengths
        for b in bodies
            b.strength .= 0
        end
    end

    if update_G
        CPs_all     = hcat(map(b -> b.controlpoints, bodies)...)
        normals_all = hcat(map(b -> b.normals, bodies)...)

        fill!(solver.G, 0)
        _G!(bodies, solver, solver.G, CPs_all, normals_all, backend;
            kerneloffset=bodies[1].kerneloffset)

        solver.Glu = lu!(solver.G)
    end

    # solve with cached LU
    sol = similar(solver.rhs)
    ldiv!(sol, solver.Glu, solver.rhs)
    
    # write solution back
    for (bi, b) in enumerate(bodies)
        r = offsets[bi]+1 : offsets[bi+1]
        if solver.boundary === :dirichlet
            @views b.strength[:,2] .= sol[r]
        else
            @views b.strength[:,1] .= sol[r]
        end

        b.CPoffset = CPoffset_old[bi]
        @views b.velocity  .= solver.Uext[:, r]
        @views b.potential .= solver.phi_ext[r]
        _solvedflag(b, true)
    end
end


### ADDED
# function FastMultipole.buffer_to_target_system!(target_system::RigidWakeBody, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
#     if PS
#         phi = target_buffer[4, i_buffer]
#         target_system.potential[i_target] += phi
#     end

#     if VS
#         vx, vy, vz = target_buffer[5, i_buffer], target_buffer[6, i_buffer], target_buffer[7, i_buffer]
#         target_system.velocity[1, i_target] += vx
#         target_system.velocity[2, i_target] += vy
#         target_system.velocity[3, i_target] += vz
#     end
# end

function _G!(bodies::Tuple, solver, G, CPs, normals, backend; kerneloffset=bodies[1].kerneloffset)
    if solver.boundary === :dirichlet
        # Dirichlet → potential influence of doublets
        _G_phi!(bodies, ConstantDoublet, G, CPs, backend;
                kerneloffset=kerneloffset)

    elseif solver.boundary === :neumann
        # Neumann → normal-velocity influence of sources
        _G_neumann!(bodies, ConstantSource, G, CPs, normals, backend;
                    kerneloffset=kerneloffset)

    else
        error("Unknown boundary type $(solver.boundary). Expected :dirichlet or :neumann.")
    end
end


function _G_sigma!(bodies::Tuple{<:AbstractBody{<:Any,NK,TF}},
                    kernel, G, CPs, normals,
                    backend::AbstractBackend=DirectBackend();
                    strength_index=1,  # sources usually in column 1
                    kerneloffset=bodies[1].kerneloffset,
                    optargs...) where {NK,TF}

    N = sum(body -> body.ncells, bodies)
    M = size(CPs, 2)

    if size(G, 1)!=M || size(G, 2)!=N
        error("Matrix G with invalid dimensions; got $(size(G)), expected ($M, $N).")
    end

    # only velocity
    derivatives_switch = FastMultipole.DerivativesSwitch(false,true,false)

    # store old strengths and set to unit sources
    old_strengths = map(body -> body.strength, bodies)

    for body in bodies
        body.strength .= zero(eltype(body.strength))
        body.strength[:, strength_index] .= 1.0
    end

    Threads.@threads for i_target in 1:M
        # target position
        tx, ty, tz = CPs[1, i_target], CPs[2, i_target], CPs[3, i_target]
        target = FastMultipole.StaticArrays.SVector{3,TF}(tx, ty, tz)

        # target normal
        nx, ny, nz = normals[1, i_target], normals[2, i_target], normals[3, i_target]

        col_offset = 0
        for body in bodies
            for i_source in 1:body.ncells
                _, u, _ = induced(target, body, i_source, derivatives_switch; kerneloffset=kerneloffset)

                # normal velocity = u ⋅ n
                un = u[1]*nx + u[2]*ny + u[3]*nz

                G[i_target, col_offset + i_source] = un
            end
            col_offset += body.ncells
        end
    end

    # restore strengths
    for (body, old_strength) in zip(bodies, old_strengths)
        body.strength .= old_strength
    end
end

### ADDED
function _G_phi!(bodies::Tuple{<:AbstractBody{<:Any,NK,TF}}, kernel, G, CPs, backend::AbstractBackend=DirectBackend(); strength_index=kernel==ConstantDoublet || kernel==VortexRing && NK>1 ? 2 : 1, kerneloffset=bodies[1].kerneloffset, optargs...) where {NK,TF}
    N = sum(body -> body.ncells, bodies)
    M = size(CPs, 2)

    if size(G, 1)!=M || size(G, 2)!=N
        error("Matrix G with invalid dimensions;"*
              " got $(size(G)), expected ($M, $N).")
    end

    # Build geometric matrix
    derivatives_switch = FastMultipole.DerivativesSwitch(true,false,false) # only potential
    
    # store old strength and set to unit
    old_strengths = map(body -> body.strength, bodies)
    
    for body in bodies
        body.strength .= zero(eltype(body.strength))
        if strength_index > 0
            body.strength[:, strength_index] .= 1.0
        else
            for i in 1:NK
                body.strength[:, i] .= 1.0
            end
        end
    end

    Threads.@threads for i_target in 1:M   ## Control points index
        # get target
        tx, ty, tz = CPs[1, i_target], CPs[2, i_target], CPs[3, i_target]
        target = FastMultipole.StaticArrays.SVector{3,TF}(tx, ty, tz)

        col_offset = 0
        for body in bodies
            for i_source in 1:body.ncells
                # compute influence
                phi, _, _ = induced(target, body, i_source, derivatives_switch; kerneloffset=kerneloffset)

                # update G
                G[i_target, col_offset + i_source] = phi
            end
            col_offset += body.ncells
        end
    end

    # restore strengths
    for (body, old_strength) in zip(bodies, old_strengths)
        body.strength .= old_strength
    end
end

### ADDED
# function calcfield_U!(targetbody::Tuple, sourcebody::Tuple;
#     backend::AbstractBackend=DirectBackend(),
#     reset=true,
#     convolve_panels=true,
#     doublet_gradient=true
# )

#     # ERROR CASES
#     for sb in sourcebody
#         @assert check_solved(sb) "Source body hasn't been solved yet. Please call `solve(...)` first."
#     end

#     # reset + geometry prep
#     for tb in targetbody
#         reset && (tb.velocity .= zero(eltype(tb.velocity)))
#     end

#     calc_normals!(targetbody)
#     calc_controlpoints!(targetbody)

#     convolve_panels && influence!(sourcebody, targetbody, backend; scalar_potential=false)

#     # doublet gradient
#     if doublet_gradient
#         for tb in targetbody
#             if has_grad_mu(tb)
#                 compute_mu_gradient!(
#                     tb.velocity,
#                     tb.controlpoints,
#                     tb.normals,
#                     tb.cells,
#                     tb.neighbor,
#                     view(tb.strength, :, get_Gammai(tb)),
#                     view(tb.shedding_full, 1:2, :),
#                     scale=0.5
#                 )
#             end
#         end
#     end

#     # return (match original behavior: single vs tuple)
#     return map(tb -> tb.velocity, targetbody)
# end

# ### ADDED
# function calcfield_Cp!(out::Arr1,
#                         bodies::Tuple{<:RigidWakeBody{<:Union{ConstantSource, ConstantDoublet, VortexRing}, 2}},
#                         Us::Arr2, Uref::Number;
#                         dphidt::Union{Nothing, AbstractVector}=nothing,
#                         correct_kuttacondition=true,
#                         clip::Union{Nothing, Function}=nothing,
#                         ) where {Arr1<:AbstractArray{<:Number,1},
#                                  Arr2<:AbstractArray{<:Number,2}}

#     @assert size(Us, 1) == 3 "Velocity matrix `Us` must have 3 rows (for 3 velocity components); got $(size(Us))."
#     @assert size(Us, 2) == length(out) "Number of columns in `Us` must match length of `out`; got $(size(Us)), $(length(out))."

#     # Calculate pressure coefficient
#     for (i, U) in enumerate(eachcol(Us))
#         out[i] += 1 - (norm(U)/Uref)^2
#     end

#     # Unsteady Bernoulli term: -2(∂φ/∂t) / V∞²
#     if !isnothing(dphidt)
#         inv_Uref2 = 2.0 / Uref^2
#         for i in eachindex(out)
#             out[i] -= inv_Uref2 * dphidt[i]
#         end
#     end

#     if correct_kuttacondition 
#         for body in bodies
#             if body isa AbstractLiftingBody
#                 for shedding in body.shedding
#                     for (pi, nia, nib, pj, nja, njb) in eachcol(shedding)

#                         if pj != -1
#                             ave = (out[pi] + out[pi+1] + out[pj] + out[pj-1]) / 4
#                             out[pi] = ave
#                             out[pi+1] = ave
#                             out[pj] = ave
#                             out[pj-1] = ave
#                         else
#                             ave = (out[pi] + out[pi+1] ) / 2
#                             out[pi] = ave
#                             out[pi+1] = ave
#                         end

#                     end
#                 end
#             end
#         end 
#     end

#     # Clip values if requested
#     if !isnothing(clip)
#         for (i, Cp) in enumerate(out)
#             out[i] = clip(Cp)
#         end
#     end

#     return out
# end

# ### ADDED
# function calcfield_F!(bodies::Tuple, Uinf::Number, rho::Number; correct_kuttacondition=true)
#     for body in bodies
#         calcfield_F!(body.F, body, calc_areas(body), body.normals, body.Cp, Uinf, rho; correct_kuttacondition=correct_kuttacondition)
#     end
#     return map(b -> b.F, bodies) 
# end

# ### ADDED
# function calcfield_LDS!(out::AbstractMatrix, bodies::Tuple,
#                         Fs::Vector{<:AbstractMatrix},
#                         Lhat::AbstractVector, Dhat::AbstractVector,
#                         Shat::AbstractVector;
#                         addfield=true)

#     fill!(out[:,3], 0.0)

#     for Fbody in Fs
#         out[:,3] .+= sum(Fbody, dims=2)[:]
#     end

#     out[:,1] = Lhat .* dot(out[:,3], Lhat)
#     out[:,2] = Dhat .* dot(out[:,3], Dhat)
#     out[:,3] = Shat .* dot(out[:,3], Shat)

#     return out
# end