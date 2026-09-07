#=
BRAINSTORM 021 — synthesize the per-job CSVs into flat tables.

Ryan's ruling, 2026-08-24: "every job writes its own CSV; synthesize
afterwards." Both campaign steps run one job per rung, each owning its own
directory under `results/`, because the alternative — many nodes appending to
one file over NFS — is a non-atomic write that CLOBBERS rows (observed
2026-08-18: R1's staircase rows destroyed by the R2–R4 writers).

The point of merging in a script rather than by `cat` is the ASSERTIONS. Every
expected row is either present exactly once or reported as MISSING by name, so
a job that was preempted, crashed or never ran is a *detectable absence* rather
than a silent hole in a published table.

    julia --project benchmark/p021_merge.jl                 # multi, all rungs
    RUNGS=R1:R2 THREADING_MODE=multi julia --project benchmark/p021_merge.jl
    STRICT=1 julia --project benchmark/p021_merge.jl        # missing row -> error

Inputs   results/phase2/<mode>/<rung>/tune_phase2.csv
         results/phase1/<mode>/<rung>/agreement.csv
         results/phase1/<mode>/<rung>/agreement_spread.csv
Outputs  results/phase2/<mode>/tune_phase2_merged.csv
         results/phase1/<mode>/agreement_merged.csv
         results/phase1/<mode>/agreement_spread_merged.csv

Merged files are DERIVED: they are rewritten from scratch on every run and must
never be appended to by a driver.
=#

const MODE = get(ENV, "THREADING_MODE", "multi")
const ALL_RUNGS = ["R1", "R2", "R3", "R4", "R5", "R6", "R7"]
const RUNGS = haskey(ENV, "RUNGS") ?
    String.(split(ENV["RUNGS"], r"[:,]")) : ALL_RUNGS
const STRICT = get(ENV, "STRICT", "0") == "1"

results(phase) = joinpath(@__DIR__, "results", phase, MODE)

problems = String[]
note(msg) = (push!(problems, msg); println("  ! $msg"))

"Split a CSV line, honouring the doubled-quote quoting `_csv_cell` emits."
function split_csv(line)
    cells = String[]
    buf = IOBuffer()
    inq = false
    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if inq
            if c == '"'
                if i < lastindex(line) && line[nextind(line, i)] == '"'
                    write(buf, '"'); i = nextind(line, i)
                else
                    inq = false
                end
            else
                write(buf, c)
            end
        elseif c == '"'
            inq = true
        elseif c == ','
            push!(cells, String(take!(buf)))
        else
            write(buf, c)
        end
        i = nextind(line, i)
    end
    push!(cells, String(take!(buf)))
    return cells
end

"""
Concatenate `<dir>/<rung>/<name>` over `RUNGS` into `<dir>/<out>`.

Every part must carry the SAME header — a schema drift between rungs means the
drivers were not all at the same commit, which would make the merged table a
quiet mixture of two definitions, so it is an error rather than a warning.
`expect(rung, rows, hdr) -> Vector{Tuple{key, label}}` names the rows that rung is
supposed to have contributed; each is checked for presence and uniqueness.
`keyof(cells, header) -> key` extracts a row's identity.
`dedup(rows, hdr, rung) -> (keep_indices, infos)`, when given, may collapse rows
that are duplicates BY DESIGN (see the resume-replay rule below) before the
uniqueness check runs; collapsed rows are excluded from the merged output and
reported as info lines, not problems.
"""
function merge_csv(dir, name, out; keyof, expect,
                   gate = (c, h, rung) -> String[], dedup = nothing)
    println("\n$name  ($MODE)")
    header = nothing
    all_rows = String[]
    found = Dict{String, Vector{Vector{String}}}()
    for rung in RUNGS
        path = joinpath(dir, rung, name)
        if !isfile(path) || filesize(path) == 0
            note("$rung: $(relpath(path, @__DIR__)) is missing or empty")
            continue
        end
        lines = [l for l in eachline(path) if !isempty(strip(l))]
        if header === nothing
            header = lines[1]
        elseif lines[1] != header
            error("$path has a different header from the earlier parts — the " *
                  "rung jobs were not all at the same schema, so merging them " *
                  "would produce a mixture of two definitions:\n  earlier: " *
                  "$header\n  here:    $(lines[1])")
        end
        rows = [split_csv(l) for l in Iterators.drop(lines, 1)]
        raw = lines[2:end]
        if dedup !== nothing
            keepidx, infos = dedup(rows, split_csv(lines[1]), rung)
            for m in infos
                println("  · $m")
            end
            rows = rows[keepidx]
            raw = raw[keepidx]
        end
        found[rung] = rows
        append!(all_rows, raw)
        println("  $rung: $(length(rows)) row(s)")
    end
    if header === nothing
        note("$name: nothing to merge — no rung produced a file")
        return 0
    end
    hdr = split_csv(header)
    for rung in RUNGS
        rows = get(found, rung, nothing)
        rows === nothing && continue
        keys_here = [keyof(r, hdr) for r in rows]
        for (k, label) in expect(rung, rows, hdr)
            n = count(==(k), keys_here)
            n == 1 && continue
            note(n == 0 ? "$rung: MISSING row $label" :
                          "$rung: $n duplicate rows for $label")
        end
    end
    # ---- quality gate --------------------------------------------------
    # Row PRESENCE is not row VALIDITY. The drivers already record why a row
    # may be unpublishable, but only in columns nobody is obliged to read, and
    # the accompanying @warn goes to a .err file that in practice goes unread
    # (2026-08-25: a healthy 14 h job and a wedged one had identical logs).
    # Promoting those columns to reported problems is what stops a best-so-far
    # row from reaching a published table as though it were a tuned optimum.
    for rung in RUNGS
        rows = get(found, rung, nothing)
        rows === nothing && continue
        for r in rows, msg in gate(r, hdr, rung)
            note(msg)
        end
    end

    outpath = joinpath(dir, out)
    open(outpath, "w") do io
        println(io, header)
        for r in all_rows
            println(io, r)
        end
    end
    println("  -> $(relpath(outpath, @__DIR__))  ($(length(all_rows)) rows)")
    return length(all_rows)
end

col(cells, hdr, name) = (i = findfirst(==(name), hdr);
    i === nothing ? error("no column '$name' in the merged header") : cells[i])

# ---- Step A: Phase 2 tuning -------------------------------------------------
# identity is (rung, mem_budget_gib). The expected budget set is not knowable
# from outside the job (MEM_BUDGETS is a submit-line choice), so it is taken
# from BUDGETS if given, else from the rows themselves — which still catches
# duplicates, the clobber signature.
const BUDGETS = haskey(ENV, "BUDGETS") ?
    String.(split(ENV["BUDGETS"], r"[:,]")) : nothing
_istrue(s) = lowercase(strip(s)) == "true"
# Resume-replay rule (observed 2026-09-07, jobs 13593015-21). A tuner that
# resumes from its trace APPENDS a fresh row for every budget it had already
# completed — marked `resumed_from_trace n_replayed=N` in `notes`, with a
# near-zero t_tune (pure replay, no descent) and the same winner knobs. That
# second row is not a clobber and not new information; the ORIGINAL row is the
# one with the real t_tune. So: within a (rung, mem_budget_gib) group, replay
# rows are collapsed onto the original iff exactly one original exists AND the
# winner knobs agree. Two originals, or a replay that disagrees on the knobs,
# is NOT collapsed — that is either the NFS-clobber signature or a descent that
# continued past the trace, and both must reach the duplicate detector. A group
# that is all replays (the original was lost) keeps only the newest row.
_is_replay(c, h) = occursin("resumed_from_trace", col(c, h, "notes"))
_knobs(c, h) = Tuple(col(c, h, k) for k in
    ("cached", "expansion_order", "multipole_acceptance", "leaf_size"))
function replay_dedup(rows, hdr, rung)
    keep = trues(length(rows))
    infos = String[]
    groups = Dict{Any, Vector{Int}}()
    for (i, r) in enumerate(rows)
        key = (col(r, hdr, "rung"), col(r, hdr, "mem_budget_gib"))
        push!(get!(groups, key, Int[]), i)
    end
    for (key, idxs) in sort!(collect(groups); by = first)
        length(idxs) == 1 && continue
        at = "$rung budget $(key[2]) GiB"
        orig = [i for i in idxs if !_is_replay(rows[i], hdr)]
        replays = [i for i in idxs if _is_replay(rows[i], hdr)]
        if length(orig) == 1
            agree = [i for i in replays if _knobs(rows[i], hdr) ==
                                           _knobs(rows[orig[1]], hdr)]
            for i in agree
                keep[i] = false
            end
            length(agree) == length(replays) && push!(infos,
                "$at: collapsed $(length(agree)) resume-replay row(s) onto " *
                "the original")
        elseif isempty(orig) && length(replays) > 1
            for i in replays[1:end-1]
                keep[i] = false
            end
            push!(infos, "$at: all $(length(replays)) rows are resume " *
                "replays (original lost); kept the newest")
        end
    end
    return findall(keep), infos
end

# A below-floor row legitimately carries cache_capped=true and
# bc_certified=false: it is the driver's explanatory answer for a budget
# cheaper than the rung's minimum configuration (e.g. 16 GiB at R6), not a
# tuned winner. Gating it would fire on every such row by construction.
_below_floor(c, h) = occursin("BELOW this rung's floor", col(c, h, "notes"))

merge_csv(results("phase2"), "tune_phase2.csv", "tune_phase2_merged.csv";
    dedup = replay_dedup,
    gate = function (c, h, rung)
        _below_floor(c, h) && return String[]
        at = "$rung budget $(col(c, h, "mem_budget_gib")) GiB"
        msgs = String[]
        _istrue(col(c, h, "tune_timed_out")) && push!(msgs,
            "$at: descent TIMED OUT — knobs are BEST-SO-FAR, not a tuned " *
            "optimum; do not publish as one (n_candidates=" *
            "$(col(c, h, "n_candidates")))")
        _istrue(col(c, h, "cache_capped")) && push!(msgs,
            "$at: winner is CACHE-CAPPED — bounded by the memory budget, " *
            "not by cost; not an unconstrained optimum")
        _istrue(col(c, h, "bc_certified")) || push!(msgs,
            "$at: row is NOT bc_certified — it is not a result")
        msgs
    end,
    keyof = (c, h) -> (col(c, h, "rung"), col(c, h, "mem_budget_gib")),
    expect = function (rung, rows, hdr)
        want = BUDGETS === nothing ?
            unique(col(r, hdr, "mem_budget_gib") for r in rows) : BUDGETS
        [((rung, b), "(rung=$rung, mem_budget_gib=$b)") for b in want]
    end)

# ---- Step B: Phase 1 agreement ---------------------------------------------
# identity is (rung, config, residual_backend). CONFIGS is per-rung at submit
# time (backslash_* are dropped above R4), so absent an explicit list the
# expectation is again "each config that did appear, appeared once".
#
# residual_backend IS PART OF THE KEY (added 2026-08-25). The O(N^2) BC
# cross-check runs the SAME configs at the SAME rung and writes to the SAME
# per-rung agreement.csv, because RESIDUAL_BACKEND is a CSV COLUMN, not a path
# component. Job 13469159 (RESIDUAL_BACKEND=direct) wrote R1's krylov_gmres and
# fgs rows; when p1-agree-R1 is re-run it appends fmm rows for those same two
# configs to that same file. Keying on (rung, config) alone would report those
# as DUPLICATE ROWS — the NFS-clobber signature — on data that is correct and
# independent. Keying on the backend too keeps a real clobber detectable while
# letting the two backends coexist.
const CONFIGS = haskey(ENV, "CONFIGS") ?
    String.(split(ENV["CONFIGS"], r"[:,]")) : nothing
# A pre-2026-08-25 CSV may predate the residual_backend column; treat a missing
# column as the default backend rather than erroring on an older file.
_backend(c, h) = "residual_backend" in h ? col(c, h, "residual_backend") : "fmm"
merge_csv(results("phase1"), "agreement.csv", "agreement_merged.csv";
    gate = function (c, h, rung)
        _istrue(col(c, h, "bc_certified")) ? String[] :
            ["$rung config $(col(c, h, "config")) " *
             "(residual_backend=$(_backend(c, h))): NOT bc_certified — " *
             "it is not a result (bc_rel_l2=$(col(c, h, "bc_rel_l2")))"]
    end,
    keyof = (c, h) -> (col(c, h, "rung"), col(c, h, "config"), _backend(c, h)),
    expect = function (rung, rows, hdr)
        # Duplicate check: every (config, backend) pair that DID appear must
        # have appeared exactly once. Built from observed pairs, because the
        # direct cross-check deliberately runs only a SUBSET of the configs
        # (13469159 ran krylov_gmres,fgs only) — demanding every config in
        # every backend would report absences that are by design.
        seen = unique([(col(r, hdr, "config"), _backend(r, hdr)) for r in rows])
        out = [((rung, cf, bk),
                "(rung=$rung, config=$cf, residual_backend=$bk)")
               for (cf, bk) in seen]
        # Coverage check: an explicitly requested config must appear under at
        # least one backend. The sentinel key cannot match any row, so a config
        # that never ran is reported as MISSING exactly once.
        if CONFIGS !== nothing
            ran = Set(cf for (cf, _) in seen)
            for cf in CONFIGS
                cf in ran || push!(out, ((rung, cf, "\0none"),
                    "(rung=$rung, config=$cf) under any residual_backend"))
            end
        end
        out
    end)

# the spread file carries one row per config plus a summary row, so config is
# not a unique key there; concatenate and check only for an empty part.
merge_csv(results("phase1"), "agreement_spread.csv",
    "agreement_spread_merged.csv";
    keyof = (c, h) -> nothing, expect = (rung, rows, hdr) -> ())

println()
if isempty(problems)
    println("MERGE CLEAN: every expected row present exactly once")
else
    println("MERGE INCOMPLETE — $(length(problems)) problem(s):")
    for p in problems
        println("  - $p")
    end
    STRICT && error("p021_merge: $(length(problems)) problem(s); see above " *
                    "(STRICT=1)")
end
