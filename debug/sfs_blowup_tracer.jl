## SFSBlowupTracer
##
## A diagnostic monitor that watches the wake's FLOWVPM ParticleField each
## timestep for the first signs of an SFS / velocity-gradient / Gamma blow-up.
## Designed to be inserted as the LAST entry in the simulate!/simulate_warmstart!
## monitor tuple so it runs after pressure/force monitors have read state.
##
## Behavior is gated by env flags so the file can live alongside the canonical
## examples without affecting normal runs:
##
##   SFS_DEBUG=1                 -- enable any tracer output (off by default)
##   SFS_DEBUG_VERBOSE=1         -- per-step ranked top-K dump
##   SFS_DEBUG_K=8               -- how many particles to dump per call
##   SFS_DEBUG_THRESHOLD_SFS=... -- |SFS| threshold (absolute, default 1.0)
##   SFS_DEBUG_THRESHOLD_J=...   -- |J|_F threshold (default 1.0e6)
##   SFS_DEBUG_THRESHOLD_GAMMA=. -- |Γ| threshold (default 1.0)
##   SFS_DEBUG_RATIO=10          -- "anomaly" = > RATIO × running median
##   SFS_DEBUG_FOCUS_INDEX=...   -- once seed is known, lock focal particle
##                                  index and dump full RK3 RHS each step
##   SFS_DEBUG_HALT=1            -- throw SFSBlowupDetected on first seed
##                                  detection so simulate! unwinds before OOM
##                                  (default: 1)
##   SFS_DEBUG_DUMP_DIR=...      -- directory for the seed-snapshot CSV
##                                  (default: debug/logs)
##
## The tracer never mutates pfield state.

struct SFSBlowupDetected <: Exception
    msg::String
end
Base.showerror(io::IO, e::SFSBlowupDetected) = print(io, "SFSBlowupDetected: ", e.msg)

using LinearAlgebra: norm, dot
import FLOWPanel as pnl
using FLOWPanel.FLOWVPM: ParticleField, get_X, get_Gamma, get_sigma, get_J,
                        get_C, get_SFS, get_SFS1, get_SFS2, get_SFS3

env_int(name, default)   = parse(Int,    get(ENV, name, string(default)))
env_float(name, default) = parse(Float64,get(ENV, name, string(default)))
env_bool(name, default)  = parse(Bool,   get(ENV, name, string(default)))

mutable struct SFSBlowupTracer
    enabled::Bool
    verbose::Bool
    k::Int
    thresh_sfs::Float64
    thresh_J::Float64
    thresh_gamma::Float64
    ratio::Float64
    focus_index::Int    # -1 = no focus locked; otherwise dump RHS each step
    halt::Bool
    dump_dir::String
    # running medians, persisted across calls
    last_sfs_median::Float64
    last_J_median::Float64
    last_gamma_median::Float64
    seed_step::Int
    seed_index::Int
end

SFSBlowupTracer() = SFSBlowupTracer(
    env_bool("SFS_DEBUG", false),
    env_bool("SFS_DEBUG_VERBOSE", false),
    env_int("SFS_DEBUG_K", 8),
    env_float("SFS_DEBUG_THRESHOLD_SFS", 1.0),
    env_float("SFS_DEBUG_THRESHOLD_J", 1.0e6),
    env_float("SFS_DEBUG_THRESHOLD_GAMMA", 1.0),
    env_float("SFS_DEBUG_RATIO", 10.0),
    env_int("SFS_DEBUG_FOCUS_INDEX", -1),
    env_bool("SFS_DEBUG_HALT", true),
    get(ENV, "SFS_DEBUG_DUMP_DIR", joinpath("debug", "logs")),
    NaN, NaN, NaN, -1, -1,
)

pnl.monitor_provides(::SFSBlowupTracer) = ()
pnl.monitor_requires(::SFSBlowupTracer) = ()

# scan a ParticleField, return (idx, mag) tuples sorted descending by mag
function _topk_by(pf::ParticleField, mag::Function, k::Int)
    n = pf.np
    n == 0 && return Tuple{Int,Float64}[]
    pairs = Vector{Tuple{Int,Float64}}(undef, n)
    @inbounds for i in 1:n
        pairs[i] = (i, mag(pf, i))
    end
    sort!(pairs; by=p -> -p[2])
    return pairs[1:min(k, n)]
end

_mag_sfs(pf, i)   = norm(get_SFS(pf, i))
_mag_J(pf, i)     = norm(get_J(pf, i))       # 9-vector view; Euclidean = Frobenius
_mag_gamma(pf, i) = norm(get_Gamma(pf, i))

# median of a small N is not in stdlib without sort allocation — quick version
function _median_mag(pf::ParticleField, mag::Function)
    n = pf.np
    n == 0 && return NaN
    vals = [mag(pf, i) for i in 1:n]
    sort!(vals)
    return vals[(n + 1) ÷ 2]
end

function _knn(pf::ParticleField, i0::Int, k::Int)
    n = pf.np
    n <= 1 && return Int[]
    x0 = get_X(pf, i0)
    d2 = [sum((get_X(pf, j) .- x0).^2) for j in 1:n]
    d2[i0] = Inf
    order = sortperm(d2)
    return order[1:min(k, n - 1)]
end

# Write a CSV with one row per particle in `indices`, recording everything
# relevant (X, sigma, Gamma, J, SFS, C, |.|, J·Γ, Γ·SFS).
function _dump_seed_csv(path::String, pf::ParticleField, step::Int,
                       seed::Int, indices::AbstractVector{Int})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "kind,step,idx,X1,X2,X3,sigma,",
            "G1,G2,G3,absG,",
            "J11,J12,J13,J21,J22,J23,J31,J32,J33,absJ,",
            "SFS1,SFS2,SFS3,absSFS,",
            "C1,C2,C3,",
            "JdotG1,JdotG2,JdotG3,GdotSFS")
        for i in indices
            kind = (i == seed) ? "seed" : "nbr"
            X = collect(get_X(pf, i))
            σ = get_sigma(pf, i)[]
            Γ = collect(get_Gamma(pf, i))
            J = collect(get_J(pf, i))
            SFSv = collect(get_SFS(pf, i))
            C = collect(get_C(pf, i))
            absG = norm(Γ); absJ = norm(J); absSFS = norm(SFSv)
            JG1 = J[1]*Γ[1] + J[2]*Γ[2] + J[3]*Γ[3]
            JG2 = J[4]*Γ[1] + J[5]*Γ[2] + J[6]*Γ[3]
            JG3 = J[7]*Γ[1] + J[8]*Γ[2] + J[9]*Γ[3]
            GS  = Γ[1]*SFSv[1] + Γ[2]*SFSv[2] + Γ[3]*SFSv[3]
            print(io, kind, ",", step, ",", i)
            for v in X;   print(io, ",", v); end
            print(io, ",", σ)
            for v in Γ;   print(io, ",", v); end
            print(io, ",", absG)
            for v in J;   print(io, ",", v); end
            print(io, ",", absJ)
            for v in SFSv; print(io, ",", v); end
            print(io, ",", absSFS)
            for v in C;   print(io, ",", v); end
            print(io, ",", JG1, ",", JG2, ",", JG3, ",", GS)
            println(io)
        end
    end
    return path
end

function _dump_particle(io::IO, pf::ParticleField, i::Int; tag::String="")
    X     = get_X(pf, i)
    Γ     = get_Gamma(pf, i)
    σ     = get_sigma(pf, i)[]
    J     = get_J(pf, i)
    SFSv  = get_SFS(pf, i)
    C     = get_C(pf, i)
    Jmat  = reshape(collect(J), 3, 3)
    Γ_dot_SFS = Γ[1]*SFSv[1] + Γ[2]*SFSv[2] + Γ[3]*SFSv[3]
    # RK3 stretching contribution J·Γ
    JΓ1 = J[1]*Γ[1] + J[2]*Γ[2] + J[3]*Γ[3]
    JΓ2 = J[4]*Γ[1] + J[5]*Γ[2] + J[6]*Γ[3]
    JΓ3 = J[7]*Γ[1] + J[8]*Γ[2] + J[9]*Γ[3]
    println(io, "  [$tag] idx=$i  X=$(round.(X; digits=6))  σ=$(σ)  |Γ|=$(norm(Γ))  |J|=$(norm(J))  |SFS|=$(norm(SFSv))  C=$(C)  Γ·SFS=$(Γ_dot_SFS)  J·Γ=($(JΓ1),$(JΓ2),$(JΓ3))")
    println(io, "    Γ=$(Γ)")
    println(io, "    SFS=$(SFSv)")
    println(io, "    J=\n$(Jmat)")
end

function (m::SFSBlowupTracer)(systems, wakes, frames, uinf, i_step::Int, dt::Real)
    m.enabled || return nothing

    # Locate the particle wake — first PanelParticleWake in wakes
    pwake = nothing
    for w in wakes
        if w isa pnl.PanelParticleWake
            pwake = w
            break
        end
    end
    pwake === nothing && return nothing
    pf = pwake.pfield
    pf.np == 0 && return nothing

    io = stdout
    step = i_step + 1

    sfs_top   = _topk_by(pf, _mag_sfs,   m.k)
    J_top     = _topk_by(pf, _mag_J,     m.k)
    gamma_top = _topk_by(pf, _mag_gamma, m.k)

    sfs_med   = _median_mag(pf, _mag_sfs)
    J_med     = _median_mag(pf, _mag_J)
    gamma_med = _median_mag(pf, _mag_gamma)

    sfs_max   = first(sfs_top)[2]
    J_max     = first(J_top)[2]
    gamma_max = first(gamma_top)[2]

    println(io, "[SFSBlowupTracer step=$step  np=$(pf.np)]  |SFS| med=$(sfs_med) max=$(sfs_max)  |J| med=$(J_med) max=$(J_max)  |Γ| med=$(gamma_med) max=$(gamma_max)")

    # detect first anomaly
    sfs_anomaly   = sfs_max   > m.thresh_sfs   || (isfinite(m.last_sfs_median)   && sfs_max   > m.ratio * m.last_sfs_median)
    J_anomaly     = J_max     > m.thresh_J     || (isfinite(m.last_J_median)     && J_max     > m.ratio * m.last_J_median)
    gamma_anomaly = gamma_max > m.thresh_gamma || (isfinite(m.last_gamma_median) && gamma_max > m.ratio * m.last_gamma_median)

    if m.verbose || sfs_anomaly || J_anomaly || gamma_anomaly
        println(io, "  top |SFS|:")
        for (idx, mag) in sfs_top
            _dump_particle(io, pf, idx; tag="SFS")
        end
        if m.verbose
            println(io, "  top |J|:")
            for (idx, mag) in J_top
                _dump_particle(io, pf, idx; tag="J")
            end
            println(io, "  top |Γ|:")
            for (idx, mag) in gamma_top
                _dump_particle(io, pf, idx; tag="Γ")
            end
        end
    end

    # First-detection: lock seed
    if m.seed_step < 0 && (sfs_anomaly || J_anomaly || gamma_anomaly)
        seed = first(sfs_top)[1]
        if J_anomaly && first(J_top)[2] >= sfs_max  # tie-break to whichever blew up earliest in magnitude
            seed = first(J_top)[1]
        end
        if gamma_anomaly && first(gamma_top)[2] >= max(sfs_max, J_max)
            seed = first(gamma_top)[1]
        end
        m.seed_step = step
        m.seed_index = seed
        # widen neighborhood for the dump beyond the per-step verbose k
        kdump = max(m.k, 16)
        nbrs = _knn(pf, seed, kdump)
        all_indices = vcat([seed], nbrs)
        println(io, "  *** SEED DETECTED: step=$step  particle_index=$seed  |SFS|=$(sfs_max) |J|=$(J_max) |Γ|=$(gamma_max) ***")
        println(io, "  $(length(nbrs)) nearest neighbors:")
        for j in nbrs
            _dump_particle(io, pf, j; tag="nbr")
        end
        # Dump full state to a CSV for later analysis
        csv_path = joinpath(m.dump_dir, "sfs_seed_step$(step)_idx$(seed).csv")
        try
            _dump_seed_csv(csv_path, pf, step, seed, all_indices)
            println(io, "  *** Seed snapshot written to $(csv_path) ***")
        catch err
            println(io, "  WARNING: failed to write seed CSV: $(err)")
        end
        if m.halt
            throw(SFSBlowupDetected(
                "Seed detected at step=$step particle_index=$seed; " *
                "|SFS|=$(sfs_max) |J|=$(J_max) |Γ|=$(gamma_max). " *
                "Snapshot: $(csv_path). Set SFS_DEBUG_HALT=0 to continue past."))
        end
    end

    # If user has locked a focal particle, dump full state always
    if m.focus_index > 0 && m.focus_index <= pf.np
        println(io, "  --- focal particle (FOCUS_INDEX=$(m.focus_index)) ---")
        _dump_particle(io, pf, m.focus_index; tag="FOCUS")
        nbrs = _knn(pf, m.focus_index, m.k)
        for j in nbrs
            _dump_particle(io, pf, j; tag="FOCUS-nbr")
        end
    end

    m.last_sfs_median = sfs_med
    m.last_J_median = J_med
    m.last_gamma_median = gamma_med
    return nothing
end
