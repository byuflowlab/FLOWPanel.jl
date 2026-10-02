#=##############################################################################
BRAINSTORM 033 B-R2: paired A/B of FGS influence-matrix setup — legacy serial
probe (`setup_threads=0`, the untouched default) vs the threaded per-leaf
populate path (`setup_threads=nthreads()`, toggled via the new
`threaded_setup` kwarg on FastGaussSeidel / pnl.FGSSolver).

Both arms run in ONE process on the SAME trees/direct list (built once via
the constructor-replay pattern of fgs_setup_profile.jl), so the comparison is
paired. The last pass of each arm's matrices are compared ELEMENT-WISE
(expected bitwise identical: the threaded path executes the identical
reset!/direct!/influence! per-column sequence on worker-private buffers, and
matrix writes are disjoint per source leaf).

Optionally (CERT_SOLVE=1, requires SKIP_B unset/0) builds full pnl.FGSSolver
solvers both ways and runs one cold FGS solve each, certifying identical
iteration counts and bitwise-identical solutions.

Usage (per (rung, threads) process):
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 THREADING_MODE=multi \
  EXPECT_JULIA_THREADS=4 BENCH_BLAS_THREADS=1 RUNG=R4 SKIP_B=1 \
  AB_K=2 ARMS=old,new CERT_SOLVE=0 \
    julia --project=. -t 4 benchmark/fgs_setup_ab.jl <outdir>

Env knobs:
  RUNG        R1..R7 (champion R4 knobs transplanted; see fgs_setup_profile.jl)
  ARMS        comma list of {old,new} (default "old,new"; cert needs both)
  AB_K        timed passes per arm (default 2; pass 0 = compile warmup, flagged)
  CERT_SOLVE  1 = full-ctor + cold-solve certification (needs SKIP_B=0)
  SKIP_B      1 = geometry-only (no RHS assembly; phases + matrix cert only)

Every pass is a CSV row; judge takes min per (arm, phase) over warmup=0 rows.
=###############################################################################

import TOML

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

const FM = pnl.FastMultipole

outdir = isempty(ARGS) ? pwd() : ARGS[1]
mkpath(outdir)
open(joinpath(outdir, "banner_ab_$(rung)_j$(banner.julia_threads).txt"), "w") do io
    println(io, banner.text)
end

champ = TOML.parsefile(joinpath(@__DIR__, "retained_r4_champion.toml"))
champ["rung"] == rung || @warn "champion TOML is $(champ["rung"]); transplanting its knobs onto $rung"

P         = champ["P"]
MAC       = Float64(champ["MAC"])
leaf      = champ["leaf"]
precision = Symbol(champ["dagteam_precision"])
shrink_   = true
recenter_ = false

ab_k       = parse(Int, get(ENV, "AB_K", "2"))
arms       = split(get(ENV, "ARMS", "old,new"), ",")
cert_solve = get(ENV, "CERT_SOLVE", "0") == "1"

################################################################################
# Shared trees + direct list (built once; populate is a pure function of them)
################################################################################

target_systems = source_systems = (rotor,)
leaf_size_v = FM.to_vector(leaf, length(source_systems))
interaction_list_method = FM.Barba()
derivatives_switches = FM.DerivativesSwitch(true, true, false, target_systems)

switches = FM.DerivativesSwitch(true, true, true, source_systems)
source_tree = FM.Tree(source_systems, FM.SourceTree(), switches;
    expansion_order=P, leaf_size=leaf_size_v, shrink=shrink_, recenter=recenter_,
    interaction_list_method)
switches = FM.DerivativesSwitch(true, true, true, target_systems)
target_tree = FM.Tree(source_tree, target_systems, switches;
    shrink=shrink_, recenter=recenter_)
FM.assert_shared_topology(target_tree, source_tree)

farfield, nearfield, self_induced = true, true, false
_m2l, direct_list = FM.build_interaction_lists(target_tree.branches,
    source_tree.branches, leaf_size_v, MAC, farfield, nearfield,
    self_induced, interaction_list_method)
direct_list = FM.sort_by_target(FM.sort_by_source(direct_list, source_tree.branches),
                                target_tree.branches)
println("$rung: trees + lists built ($(length(direct_list)) direct-list entries)")

################################################################################
# Arm passes
################################################################################

setup_threads_for(arm) = arm == "new" ? Threads.nthreads() : 0

kept = Dict{String,Any}()   # arm => (nonself_matrices, sorted_list, self_matrices)

csv_path = joinpath(outdir, "fgs_setup_ab_$(rung)_j$(banner.julia_threads).csv")
io = open(csv_path, "w")
println(io, "rung,n_panels,julia_threads,blas_threads,arm,setup_threads,pass,warmup,phase,t_s,gc_bytes")

for arm in arms
    st = setup_threads_for(arm)
    for pass in 0:ab_k
        GC.gc(); GC.gc()

        t0 = time_ns(); a0 = Base.gc_bytes()
        nonself_matrices, sorted_list = FM.nonself_influence_matrices(
            target_tree.buffers, source_tree.buffers, source_systems,
            target_tree, source_tree, direct_list, derivatives_switches;
            setup_threads=st)
        t_nonself = (time_ns()-t0)/1e9; b_nonself = Int(Base.gc_bytes()-a0)

        t0 = time_ns(); a0 = Base.gc_bytes()
        self_matrices = FM.self_influence_matrices(target_tree.buffers,
            source_tree.buffers, source_systems, target_tree, source_tree,
            derivatives_switches; setup_threads=st)
        t_self = (time_ns()-t0)/1e9; b_self = Int(Base.gc_bytes()-a0)

        for (phase, t, b) in (("nonself", t_nonself, b_nonself),
                              ("self", t_self, b_self))
            println(io, "$rung,$(rotor.ncells),$(banner.julia_threads)," *
                "$(banner.blas_threads),$arm,$st,$pass,$(pass == 0 ? 1 : 0)," *
                "$phase,$t,$b")
        end
        println("arm=$arm pass=$pass: nonself $(round(t_nonself, digits=2)) s, " *
                "self $(round(t_self, digits=2)) s")
        flush(io)

        if pass == ab_k
            kept[arm] = (nonself_matrices, sorted_list, self_matrices)
        end
    end
end

################################################################################
# Matrix equivalence certification (element-wise; expect bitwise)
################################################################################

function maxabsdiff(a::AbstractVector, b::AbstractVector)
    length(a) == length(b) || return Inf
    m = 0.0
    @inbounds for i in eachindex(a)
        m = max(m, abs(Float64(a[i]) - Float64(b[i])))
    end
    return m
end

if haskey(kept, "old") && haskey(kept, "new")
    no, slo, so = kept["old"]
    nn, sln, sn = kept["new"]
    cert = [
        ("sorted_list_equal",      slo == sln),
        ("nonself_sizes_equal",    no.sizes == nn.sizes),
        ("nonself_data_bitwise",   no.data == nn.data),
        ("nonself_rhs_bitwise",    no.rhs == nn.rhs),
        ("self_sizes_equal",       so.sizes == sn.sizes),
        ("self_data_bitwise",      so.data == sn.data),
        ("self_rhs_bitwise",       so.rhs == sn.rhs),
    ]
    d_nonself = maxabsdiff(no.data, nn.data)
    d_self    = maxabsdiff(so.data, sn.data)
    for (name, ok) in cert
        println(io, "$rung,$(rotor.ncells),$(banner.julia_threads)," *
            "$(banner.blas_threads),cert,0,0,0,$name,$(ok ? 1 : 0),0")
        println("CERT $name = $ok")
    end
    println(io, "$rung,$(rotor.ncells),$(banner.julia_threads)," *
        "$(banner.blas_threads),cert,0,0,0,nonself_max_abs_diff,$d_nonself,0")
    println(io, "$rung,$(rotor.ncells),$(banner.julia_threads)," *
        "$(banner.blas_threads),cert,0,0,0,self_max_abs_diff,$d_self,0")
    println("CERT nonself max|Δ| = $d_nonself, self max|Δ| = $d_self")
    all(last.(cert)) || @error "MATRIX CERTIFICATION FAILED — see CSV"
end

# release phase-arm matrices before the ctor certification
empty!(kept); GC.gc(); GC.gc()

################################################################################
# Optional full-ctor + cold-solve certification (SKIP_B=0 rungs only)
################################################################################

if cert_solve
    b_skipped && error("CERT_SOLVE=1 requires SKIP_B=0 (a solvable case)")

    function build_solver(threaded::Bool)
        GC.gc(); GC.gc()
        t0 = time_ns()
        solver = pnl.FGSSolver(rotor; expansion_order=P,
            multipole_acceptance=MAC, leaf_size=leaf,
            inner_iterations=champ["inner"], max_iterations=champ["max_iterations"],
            tolerance=champ["tolerance"], rlx=champ["rlx"], shrink=shrink_,
            recenter=recenter_, reverse_pass=false, cache_leaf_lu=true,
            sweep_order=:dagteam, dagteam_precision=precision,
            verbose=false, project_solution=false, solution_history_length=0,
            threaded_setup=threaded)
        t = (time_ns()-t0)/1e9
        println(io, "$rung,$(rotor.ncells),$(banner.julia_threads)," *
            "$(banner.blas_threads),$(threaded ? "new" : "old")," *
            "$(threaded ? Threads.nthreads() : 0),1,0,full_ctor,$t,0")
        println("full ctor (threaded_setup=$threaded): $(round(t, digits=2)) s")
        flush(io)
        return solver
    end

    function cold_solve!(solver)
        reset_cold!()
        rotor.core_size = rotor.core_size_panel
        solver.niter = 0; solver.solved = false
        t = @elapsed pnl._solve!(rotor, solver)
        return (t, copy(rotor.strength[:, solution_column]), solver.niter, solver.solved)
    end

    solver_old = build_solver(false)
    solver_new = build_solver(true)

    # solver-internal matrices must also match bitwise
    m_ok = solver_old.fgs.nonself_matrices.data == solver_new.fgs.nonself_matrices.data &&
           solver_old.fgs.self_matrices.data == solver_new.fgs.self_matrices.data
    # B-I2: the dagteam plan's repacked split storage too (threaded repack
    # under threaded_setup must be bitwise-identical to the serial pass)
    d_o, d_n = solver_old.fgs.dagteam, solver_new.fgs.dagteam
    plan_ok = d_o === nothing ? d_n === nothing :
        (d_o.Lmat == d_n.Lmat && d_o.Umat == d_n.Umat)
    m_ok &= plan_ok
    println("CERT solver_matrices_bitwise = $m_ok (dagteam Lmat/Umat = $plan_ok)")
    println(io, "$rung,$(rotor.ncells),$(banner.julia_threads)," *
        "$(banner.blas_threads),cert,0,0,0,solver_matrices_bitwise,$(m_ok ? 1 : 0),0")

    t_old, x_old, n_old, ok_old = cold_solve!(solver_old)
    t_new, x_new, n_new, ok_new = cold_solve!(solver_new)
    x_bit = x_old == x_new
    rel = norm(x_old - x_new) / max(norm(x_old), eps())
    println("CERT solve: niter old=$n_old new=$n_new solved old=$ok_old new=$ok_new " *
            "x bitwise=$x_bit rel=$rel (t $(round(t_old,digits=2)) vs $(round(t_new,digits=2)) s)")
    println(io, "$rung,$(rotor.ncells),$(banner.julia_threads)," *
        "$(banner.blas_threads),cert,0,0,0,solve_niter_old,$n_old,0")
    println(io, "$rung,$(rotor.ncells),$(banner.julia_threads)," *
        "$(banner.blas_threads),cert,0,0,0,solve_niter_new,$n_new,0")
    println(io, "$rung,$(rotor.ncells),$(banner.julia_threads)," *
        "$(banner.blas_threads),cert,0,0,0,solve_x_bitwise,$(x_bit ? 1 : 0),0")
    println(io, "$rung,$(rotor.ncells),$(banner.julia_threads)," *
        "$(banner.blas_threads),cert,0,0,0,solve_x_reldiff,$rel,0")
    (m_ok && x_bit && n_old == n_new && ok_old == ok_new) ||
        @error "SOLVE CERTIFICATION FAILED — see CSV"
end

close(io)
println("CSV written to $csv_path")
