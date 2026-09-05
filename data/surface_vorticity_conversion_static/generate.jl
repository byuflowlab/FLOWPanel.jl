# Reproduce the simulation-free Stage 4b probe and ledger evidence. Run from
# the repository root with at most six threads:
#
#     julia --project data/surface_vorticity_conversion_static/generate.jl
#
# The script prints CSV rows; it does not time-march, shed, solve, or modify the
# checked-in evidence files.
using LinearAlgebra
using StaticArrays
import FLOWPanel as pnl
import FLOWVPM
import FastMultipole
include(joinpath(@__DIR__, "..", "..", "test", "test_helpers.jl"))

const SIGMA = 0.06
const OVERLAP = 1.3
const RATIOS = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0)

function prescribed_sheet(strategy; nrows=6, nspan=12)
    body = make_dirichlet_diamond_body(; nspan)
    conversion = strategy === :legacy ? pnl.LegacyEdgeJumpConversion() :
        pnl.SurfaceVorticityConversion(SIGMA; overlap=OVERLAP,
            attribution=strategy)
    wake = strategy === :legacy ?
        pnl.PanelParticleWake(body; nwakerows=nrows, max_particles=200_000,
            core_size=SIGMA, conversion,
            method_trailing=pnl.SigmaOverlap(SIGMA, OVERLAP),
            method_unsteady=pnl.SigmaOverlap(SIGMA, OVERLAP)) :
        pnl.PanelParticleWake(body; nwakerows=nrows, max_particles=200_000,
            core_size=SIGMA, conversion)
    nodes = wake.panel_wake.nodes[1]
    strength = wake.panel_wake.strength[1]
    nc = size(nodes,3)-1
    mapx(s) = 3expm1(1.1s)/expm1(1.1)
    mapy(q) = 1.5sinh(0.8q)/sinh(0.8)
    for r in axes(nodes,2), c in axes(nodes,3)
        s=(r-1)/nrows; q=2(c-1)/nc-1
        nodes[:,r,c] .= (mapx(s), mapy(q), 0.10sin(pi*s)*cos(pi*q/2))
    end
    for r in axes(strength,2), c in axes(strength,3)
        rr=min(r,size(nodes,2)-1)
        p=0.25*(SVector{3}(nodes[:,rr,c])+SVector{3}(nodes[:,rr+1,c])+
                SVector{3}(nodes[:,rr+1,c+1])+SVector{3}(nodes[:,rr,c+1]))
        r == size(strength,2) &&
            (p += 0.5*(SVector{3}(nodes[:,end,c])-SVector{3}(nodes[:,end-1,c])))
        x,y,z=p
        strength[1,r,c]=0.7+0.35x+0.08x^2+0.18sin(1.2y)+0.04x*y+0.1z
    end
    wake.panel_wake.nwakes[]=nrows
    return body,wake
end

function direct_velocity(positions,sources)
    probes=FastMultipole.ProbeSystem(length(positions),Float64)
    for i in eachindex(positions)
        probes.position[i]=positions[i]
        probes.scalar_potential[i]=0.0
        probes.gradient[i]=zero(SVector{3,Float64})
        probes.hessian[i]=zero(SMatrix{3,3,Float64,9})
    end
    pnl.influence!((probes,),sources,pnl.DirectBackend();
        scalar_potential=false,gradient=true,hessian=(false,))
    return copy(probes.gradient)
end

function evidence(strategy)
    body,wake=prescribed_sheet(strategy)
    pw=wake.panel_wake; N=pw.nwakes[]; nodes=pw.nodes[1]
    positions=SVector{3,Float64}[]; family=Symbol[]; ratio_of=Float64[]
    for ratio in RATIOS
        for (label,q) in ((:root,-0.96),(:interior,-0.35),(:interior,0.0),
                          (:interior,0.35),(:tip,0.96))
            j=clamp(round(Int,1+(size(nodes,3)-1)*(q+1)/2),1,size(nodes,3)-1)
            v1=SVector{3}(nodes[:,N,j]); v2=SVector{3}(nodes[:,N+1,j])
            v3=SVector{3}(nodes[:,N+1,j+1]); v4=SVector{3}(nodes[:,N,j+1])
            x=0.25*(v1+v2+v3+v4); n=cross(v2-v1,v4-v1); n/=norm(n)
            for side in (-1,1)
                push!(positions,x+side*ratio*SIGMA*n)
                push!(family,label); push!(ratio_of,ratio)
            end
        end
        for j in (max(1,div(size(nodes,3)-1,3)),
                  max(1,div(2*(size(nodes,3)-1),3)))
            x=0.5*(SVector{3}(nodes[:,N,j])+SVector{3}(nodes[:,N,j+1]))
            a=SVector{3}(nodes[:,N+1,j])-SVector{3}(nodes[:,N,j])
            b=SVector{3}(nodes[:,N,j+1])-SVector{3}(nodes[:,N,j])
            n=cross(a,b); n/=norm(n)
            for side in (-1,1)
                push!(positions,x+side*ratio*SIGMA*n)
                push!(family,:handoff); push!(ratio_of,ratio)
            end
        end
    end
    sheet=direct_velocity(positions,pnl.get_sources(pw))
    strategy === :legacy ? pnl._convert_to_particles!(wake) :
        pnl._convert_to_particles!(wake,body)
    pw.nwakes[]=N-1
    hybrid=direct_velocity(positions,(pnl.get_sources(pw)...,wake.pfield))
    for ratio in RATIOS
        ids=findall(i->family[i]===:handoff && ratio_of[i]==ratio,
                    eachindex(positions))
        delta=[norm(hybrid[i]-sheet[i]) for i in ids]
        rms=sqrt(sum(abs2,delta)/length(ids))/
            sqrt(sum(norm(sheet[i])^2 for i in ids)/length(ids))
        mx=maximum(delta)/maximum(norm(sheet[i]) for i in ids)
        println(join((strategy,ratio,rms,mx,wake.pfield.np),','))
    end
    G=zero(SVector{3,Float64}); I=zero(SVector{3,Float64})
    for i in 1:wake.pfield.np
        x=SVector{3}(wake.pfield.particles[FLOWVPM.X_INDEX,i])
        g=SVector{3}(wake.pfield.particles[FLOWVPM.GAMMA_INDEX,i])
        G+=g; I+=cross(x,g)
    end
    println("ledger,",strategy,",",wake.pfield.np,",",join(G,','),",",join(I,','))
end

println("strategy,d_over_sigma,handoff_normalized_rms,handoff_normalized_max,particle_count")
foreach(evidence,(:legacy,:upstream,:split))
