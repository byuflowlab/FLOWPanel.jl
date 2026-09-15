using ReadVTK, Printf, Statistics
const DIR = "/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/data/rotor_hover_pressure_comparison/rotor_hover_pressure_comparison_wake1_particles"
const SIG_SHED = 0.031154
const PHI = 3.5
const MAXRATIO = 2.0
fname(s) = joinpath(DIR, "rotor_hover_pressure_comparison_wake1_particles.$s.vtp")
function load(s)
    vtk = VTKFile(fname(s)); pd = get_point_data(vtk)
    get_points(vtk), Float64.(vec(get_data(pd["sigma"]))), Float64.(vec(get_data(pd["rsplit_sigma_0"])))
end
qs(v) = quantile(v, [0.05,0.25,0.5,0.75,0.95])

# emulate merge_particles! same-cell disjoint nearest pairing
function would_merge_pairs(X, sig)
    N = length(sig)
    cell = (1/PHI) * mean(sig)
    o = (minimum(X[1,:]), minimum(X[2,:]), minimum(X[3,:]))
    key(i) = (floor(Int,(X[1,i]-o[1])/cell), floor(Int,(X[2,i]-o[2])/cell), floor(Int,(X[3,i]-o[3])/cell))
    cells = Dict{NTuple{3,Int},Vector{Int}}()
    for i in 1:N
        push!(get!(cells, key(i), Int[]), i)
    end
    paired = falses(N); pairs = Tuple{Int,Int}[]
    for (_, idxs) in cells
        length(idxs) < 2 && continue
        for (ai,ia) in enumerate(idxs)
            paired[ia] && continue
            nearest = 0; nd2 = Inf
            for b in ai+1:length(idxs)
                ib = idxs[b]; paired[ib] && continue
                smin = min(sig[ia],sig[ib]); smax = max(sig[ia],sig[ib])
                smax/smin > MAXRATIO && continue
                r = smin/PHI
                d2 = (X[1,ib]-X[1,ia])^2 + (X[2,ib]-X[2,ia])^2 + (X[3,ib]-X[3,ia])^2
                (d2 < r*r && d2 < nd2) && (nearest = ib; nd2 = d2)
            end
            nearest != 0 && (paired[ia]=true; paired[nearest]=true; push!(pairs,(ia,nearest)))
        end
    end
    return pairs
end

for s in (322, 400, 467)
    X, sig, s0 = load(s)
    pr = would_merge_pairs(X, sig)
    println("\n=== STEP $s (np=$(length(sig))): would-merge pairs next pass = $(length(pr))")
    isempty(pr) && continue
    smins = [min(sig[i],sig[j]) for (i,j) in pr]
    s0p = vcat([s0[i] for (i,_) in pr], [s0[j] for (_,j) in pr])
    gp = vcat([sig[i]/s0[i] for (i,_) in pr], [sig[j]/s0[j] for (_,j) in pr])
    @printf("  pair sigma_min/shed quantiles (5/25/50/75/95): %s\n", join([@sprintf("%.2f",q) for q in qs(smins./SIG_SHED)], " "))
    @printf("  member sigma_0 composition: ~shed(2%%)=%.3f split-child(<0.9shed)=%.3f merge-product(>1.1shed)=%.3f\n",
        count(x->abs(x-SIG_SHED)<0.02*SIG_SHED, s0p)/length(s0p),
        count(x->x<0.9*SIG_SHED, s0p)/length(s0p),
        count(x->x>1.1*SIG_SHED, s0p)/length(s0p))
    @printf("  member growth sigma/sigma_0 q: %s\n", join([@sprintf("%.3f",q) for q in qs(gp)], " "))
    # sigma_min histogram
    edges = collect(0.0:0.25:3.5)
    h = [count(x-> e1 <= x/SIG_SHED < e2, smins) for (e1,e2) in zip(edges[1:end-1], edges[2:end])]
    println("  hist sigma_min/shed @0.25 bins from 0: ", h)
end
