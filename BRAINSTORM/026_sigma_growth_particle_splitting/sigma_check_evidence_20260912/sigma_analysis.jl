using ReadVTK, Printf, Statistics

const DIR = "/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/data/rotor_hover_pressure_comparison/rotor_hover_pressure_comparison_wake1_particles"
const R = 0.119
const SIG_SHED = 0.031154
const PHI = 3.5
const SIG_STAR = PHI * 0.02 * R   # crossover vs legacy absolute radius
const MAXRATIO = 2.0

fname(s) = joinpath(DIR, "rotor_hover_pressure_comparison_wake1_particles.$s.vtp")

function load(s)
    vtk = VTKFile(fname(s))
    pd = get_point_data(vtk)
    X = get_points(vtk)                       # 3×N
    sig = vec(get_data(pd["sigma"]))
    s0  = vec(get_data(pd["rsplit_sigma_0"]))
    return X, Float64.(sig), Float64.(s0)
end

qs(v) = quantile(v, [0.05,0.25,0.5,0.75,0.95])

# imminent-merge pairs: dist < sigma_min/PHI, sigma ratio <= MAXRATIO
function imminent_pairs(X, sig)
    N = length(sig)
    pairs = Tuple{Int,Int}[]
    # cell binning at max gate width
    rmax = maximum(sig)/PHI
    for i in 1:N-1
        xi, yi, zi = X[1,i], X[2,i], X[3,i]
        si = sig[i]
        for j in i+1:N
            sj = sig[j]
            smin = min(si,sj); smax = max(si,sj)
            smax/smin > MAXRATIO && continue
            r = smin/PHI
            dx = X[1,j]-xi; dy = X[2,j]-yi; dz = X[3,j]-zi
            d2 = dx*dx+dy*dy+dz*dz
            d2 < r*r && push!(pairs, (i,j))
        end
    end
    return pairs
end

function detail(s)
    X, sig, s0 = load(s)
    N = length(sig)
    println("\n===== STEP $s  (np=$N) =====")
    @printf("sigma quantiles [m] (5/25/50/75/95%%): %s\n", join([@sprintf("%.5f",q) for q in qs(sig)], " "))
    @printf("sigma/sig_shed quantiles: %s\n", join([@sprintf("%.2f",q) for q in qs(sig./SIG_SHED)], " "))
    @printf("frac sigma > sig* (%.5f m): %.4f\n", SIG_STAR, count(>(SIG_STAR), sig)/N)
    @printf("frac sigma within 10%% of shed sigma: %.4f\n", count(x->abs(x-SIG_SHED)<0.1*SIG_SHED, sig)/N)
    # birth-sigma composition
    @printf("sigma_0 quantiles [m]: %s\n", join([@sprintf("%.5f",q) for q in qs(s0)], " "))
    fshed = count(x->abs(x-SIG_SHED)<0.02*SIG_SHED, s0)/N
    fsmall = count(x-> x < 0.9*SIG_SHED, s0)/N        # split children born smaller
    fbig = count(x-> x > 1.1*SIG_SHED, s0)/N          # merge products born bigger
    @printf("sigma_0 composition: ~shed(±2%%)=%.3f  <0.9*shed(split-child)=%.3f  >1.1*shed(merge-product)=%.3f\n", fshed, fsmall, fbig)
    # age proxy: growth since last event
    g = sig ./ s0
    @printf("sigma/sigma_0 (growth since last event) quantiles: %s\n", join([@sprintf("%.3f",q) for q in qs(g)], " "))
    # imminent merge pairs
    pr = imminent_pairs(X, sig)
    println("imminent-merge pairs: ", length(pr))
    if !isempty(pr)
        smins = [min(sig[i],sig[j]) for (i,j) in pr]
        gpair = vcat([g[i] for (i,_) in pr], [g[j] for (_,j) in pr])
        s0pair = vcat([s0[i] for (i,_) in pr], [s0[j] for (_,j) in pr])
        @printf("  pair sigma_min quantiles [m]: %s\n", join([@sprintf("%.5f",q) for q in qs(smins)], " "))
        @printf("  pair sigma_min/sig_shed quantiles: %s\n", join([@sprintf("%.2f",q) for q in qs(smins./SIG_SHED)], " "))
        @printf("  pair-member growth sigma/sigma_0 quantiles: %s\n", join([@sprintf("%.3f",q) for q in qs(gpair)], " "))
        @printf("  pair-member sigma_0 composition: ~shed=%.3f  split-child=%.3f  merge-product=%.3f\n",
            count(x->abs(x-SIG_SHED)<0.02*SIG_SHED, s0pair)/length(s0pair),
            count(x-> x<0.9*SIG_SHED, s0pair)/length(s0pair),
            count(x-> x>1.1*SIG_SHED, s0pair)/length(s0pair))
        # fresh-event members (sigma within 5% of sigma_0) = just split/merged/shed
        ffresh = count(x-> x < 1.05, gpair)/length(gpair)
        @printf("  frac pair-members fresh (sigma/sigma_0<1.05): %.3f\n", ffresh)
    end
    # histogram (10 bins in sigma/shed)
    edges = range(0.5, 3.5, length=13)
    h = [count(x-> e1 <= x/SIG_SHED < e2, sig) for (e1,e2) in zip(edges[1:end-1], edges[2:end])]
    println("hist sigma/sig_shed bins ", round.(collect(edges[1:end-1]),digits=2), ":")
    println("  ", h)
end

for s in (322, 400, 467)
    detail(s)
end

# time series
println("\n===== TIME SERIES (B surviving steps) =====")
println("step,np,q10,q50,q90,frac_fresh,imminent_pairs")
for s in 322:15:467
    X, sig, s0 = load(s)
    g = sig ./ s0
    pr = imminent_pairs(X, sig)
    q = quantile(sig, [0.1,0.5,0.9])
    @printf("%d,%d,%.5f,%.5f,%.5f,%.3f,%d\n", s, length(sig), q[1], q[2], q[3], count(x->x<1.05,g)/length(g), length(pr))
end
