using ReadVTK, Printf, Statistics
const DIR = "/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/data/rotor_hover_pressure_comparison/rotor_hover_pressure_comparison_wake1_particles"
const SIG_SHED = 0.031154; const PHI = 3.5; const MAXRATIO = 2.0
fname(s) = joinpath(DIR, "rotor_hover_pressure_comparison_wake1_particles.$s.vtp")
function load(s)
    vtk = VTKFile(fname(s)); pd = get_point_data(vtk)
    get_points(vtk), Float64.(vec(get_data(pd["sigma"]))), Float64.(vec(get_data(pd["rsplit_sigma_0"])))
end
function would_merge_pairs(X, sig)
    N = length(sig); cell = (1/PHI)*mean(sig)
    o = (minimum(X[1,:]), minimum(X[2,:]), minimum(X[3,:]))
    cells = Dict{NTuple{3,Int},Vector{Int}}()
    for i in 1:N
        k = (floor(Int,(X[1,i]-o[1])/cell), floor(Int,(X[2,i]-o[2])/cell), floor(Int,(X[3,i]-o[3])/cell))
        push!(get!(cells,k,Int[]), i)
    end
    paired = falses(N); pairs = Tuple{Int,Int,Float64}[]
    for (_,idxs) in cells
        length(idxs)<2 && continue
        for (ai,ia) in enumerate(idxs)
            paired[ia] && continue
            nearest=0; nd2=Inf
            for b in ai+1:length(idxs)
                ib=idxs[b]; paired[ib] && continue
                smin=min(sig[ia],sig[ib]); smax=max(sig[ia],sig[ib])
                smax/smin>MAXRATIO && continue
                d2=(X[1,ib]-X[1,ia])^2+(X[2,ib]-X[2,ia])^2+(X[3,ib]-X[3,ia])^2
                (d2<(smin/PHI)^2 && d2<nd2) && (nearest=ib; nd2=d2)
            end
            nearest!=0 && (paired[ia]=true; paired[nearest]=true; push!(pairs,(ia,nearest,sqrt(nd2))))
        end
    end
    pairs
end
for s in (322,400,467)
    X,sig,s0 = load(s)
    pr = would_merge_pairs(X,sig)
    n = length(pr)
    sib = [(i,j,d) for (i,j,d) in pr if abs(s0[i]-s0[j]) < 1e-5*max(s0[i],s0[j])]
    sibsplit = [(i,j,d) for (i,j,d) in sib if s0[i] < 0.9*SIG_SHED]
    @printf("\nSTEP %d: pairs=%d  equal-sigma0 pairs=%d (%.0f%%)  of which split-sigma0=%d\n",
        s, n, length(sib), 100*length(sib)/max(n,1), length(sibsplit))
    if !isempty(sibsplit)
        # distance relative to placement spacing sigma0/3.0 and gate sigma/3.5
        rel = [d/(s0[i]/3.0) for (i,j,d) in sibsplit]
        @printf("  sibling-split pair dist / (sigma_c/3.0 placement): q(5,50,95)= %.2f %.2f %.2f\n",
            quantile(rel,0.05), quantile(rel,0.5), quantile(rel,0.95))
        gr = vcat([sig[i]/s0[i] for (i,_,_) in sibsplit],[sig[j]/s0[j] for (_,j,_) in sibsplit])
        @printf("  sibling members sigma/sigma_0: q(5,50,95)= %.4f %.4f %.4f\n",
            quantile(gr,0.05), quantile(gr,0.5), quantile(gr,0.95))
    end
    # non-sibling pairs composition
    nonsib = [(i,j) for (i,j,d) in pr if abs(s0[i]-s0[j]) >= 1e-5*max(s0[i],s0[j])]
    if !isempty(nonsib)
        s0p = vcat([s0[i] for (i,_) in nonsib],[s0[j] for (_,j) in nonsib])
        @printf("  non-sibling members: split-child=%.2f  ~shed=%.2f  merge-product=%.2f\n",
            count(x->x<0.9*SIG_SHED,s0p)/length(s0p),
            count(x->abs(x-SIG_SHED)<0.02*SIG_SHED,s0p)/length(s0p),
            count(x->x>1.1*SIG_SHED,s0p)/length(s0p))
    end
end
