# Offline iteration on the dumped R4 dagteam graph (fgs_dag_L_graph/edges CSVs).
# (a) task counts for the theta-cutoff variants; (b) greedy critical-path-driven
# selective split: repeatedly split every leaf on the current critical path.
dir = ARGS[1]
lines = readlines(joinpath(dir, "fgs_dag_L_graph_R4.csv"))[2:end]
n = length(lines)
nof = zeros(Int, n); ptot = zeros(Int, n)
for l in lines
    p = split(l, ','); i = parse(Int, p[1])
    nof[i] = parse(Int, p[2]); ptot[i] = parse(Int, p[3])
end
preds = [Int[] for _ in 1:n]
for l in readlines(joinpath(dir, "fgs_dag_L_edges_R4.csv"))[2:end]
    p = split(l, ','); push!(preds[parse(Int, p[1])], parse(Int, p[2]))
end
sz = 4
g = [sz*nof[i]*ptot[i] for i in 1:n]
lu = [sz*nof[i]^2 for i in 1:n]
red = [sz*nof[i]*length(preds[i]) for i in 1:n]
edge(i,j) = sz*nof[i]*nof[j]

function Lpath(S::Vector{Bool})
    fin = zeros(n); pre = zeros(Int, n)
    for i in 1:n
        if S[i]
            t, ti = 0.0, 0
            for j in preds[i]
                v = fin[j] + edge(i,j); v > t && (t = v; ti = j)
            end
            fin[i] = t + red[i] + lu[i]; pre[i] = ti
        else
            a, ai = 0.0, 0
            for j in preds[i]
                fin[j] > a && (a = fin[j]; ai = j)
            end
            fin[i] = a + g[i] + lu[i]; pre[i] = ai
        end
    end
    L, iend = findmax(fin)
    path = Int[]; i = iend
    while i != 0; push!(path, i); i = pre[i]; end
    return L, path
end

mb(x) = round(x/1e6, digits=1)
S = falses(n) |> collect
L0, _ = Lpath(S)
println("| round | leaves split | tasks/sweep | L (MB) | L_node/L |")
println("|---|---|---|---|---|")
for r in 0:40
    L, path = Lpath(S)
    nt = n + sum(length(preds[i]) for i in 1:n if S[i]; init=0)
    r % 2 == 0 && println("| $r | $(count(S)) | $nt | $(mb(L)) | $(round(L0/L, digits=2)) |")
    added = false
    for i in path
        S[i] || (S[i] = true; added = true)
    end
    added || (println("(converged at round $r)"); break)
end

# theta-cutoff task counts (edges >= theta become tasks; rest aggregate per leaf)
println("\n| theta | split edges | tasks/sweep |")
println("|---|---|---|")
for th in (0, 4096, 16384, 65536, 262144)
    ne = sum(count(j -> edge(i,j) >= th, preds[i]) for i in 1:n)
    println("| $(div(th,1024))KB | $ne | $(n + ne) |")
end
