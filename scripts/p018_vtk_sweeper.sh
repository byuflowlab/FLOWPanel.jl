#!/usr/bin/env bash
#
# p018_vtk_sweeper.sh -- restart-set-aware VTK sweeper for the 400 G home budget.
#
# Encodes ruling 10 as amended (BRAINSTORM/018_.../ops_reference.md) mechanically so
# a watchdog never has to do per-file bookkeeping by hand:
#
#   * CSVs, *_case_metadata.toml, *.metadata.toml and monitors/ are NEVER touched.
#   * Every run keeps its final $KEEP_STEPS RESTARTABLE steps, with ALL files at
#     those steps (.vtm indices and their pieces together) so the retained state is
#     both warm-startable and openable in ParaView.  A restartable step S is one
#     where all four paths the warmstart loader reads exist:
#         <run>_body1/<run>_body1.<S>.vtu
#         <run>_wake1/<run>_wake1.1.<S>.vts
#         <run>_wake1/<run>_wake1.2.<S>.vts
#         <run>_wake1_particles/<run>_wake1_particles.<S>.vtp
#     (src/FLOWPanel_warmstart.jl; filaments are not read but are retained at the
#     kept steps for debugging.)
#   * A run with NO complete restartable step is skipped, never guessed at.
#   * Runs named in the protect list are skipped entirely.
#   * LIVE runs that are not protected ARE swept (Ryan, 2026-08-05: "for any run
#     not in the protect list, you can delete files while it runs, so long as you
#     leave the newest timesteps (keep count: $KEEP_STEPS, 288 since 2026-08-24) -- no need to wait for the run to
#     complete").  Files younger than $LIVE_QUIET_SEC are still spared so the
#     step currently being written is never touched, and the $KEEP_STEPS most recent
#     restartable steps are kept exactly as for a closed run.  Together these
#     mean the restart file that killed job 13036477 on 2026-08-04 cannot be the
#     one removed.  Pass --skip-live for the old conservative behaviour.
#
# Dry-run by default.  Pass --apply to actually delete.
#
# Usage:
#   scripts/p018_vtk_sweeper.sh [--apply] [--skip-live] [--only RUN] [--keep N]
#                               [--last-resort]
#
# Retention ladder (Ryan, 2026-08-24; all G tiers dropped 100 G by Ryan 2026-08-25): 288 steps normally; 72 at the >300 G
# escalation tier; 36 is the absolute floor and is a genuine last resort,
# reachable only with --last-resort.  Anything below 36, or a non-integer,
# is refused outright -- see the keep-value validation block below.
#
# Exit codes: 2 bad usage, 3 no protect list, 4 bad --keep, 5 unmatched --only,
# 6 another sweeper already holds the lock.
#
# Run from the top level of the FLOWPanel.jl checkout.

set -euo pipefail

# ------------------------------------------------------------------ locking
# Two concurrent --apply runs would race on the same tree; a watchdog cron can
# fire while the previous cycle is still sweeping.  Re-exec under flock so only
# one sweeper touches data/ at a time.  Guarded on flock(1) being present so the
# script stays runnable on macOS for syntax and fixture testing.
LOCK="${TMPDIR:-/tmp}/.p018_vtk_sweeper.$(id -u).lock"
if command -v flock >/dev/null 2>&1 && [[ -z "${SWEEPER_LOCK_HELD:-}" ]]; then
    export SWEEPER_LOCK_HELD=1
    exec flock -n -E 6 "$LOCK" "$0" "$@"
fi

KEEP_STEPS="${KEEP_STEPS:-288}"     # restartable steps retained per run (10->36 Ryan 2026-08-20; 36->288 Ryan 2026-08-24)
LIVE_MTIME_MIN="${LIVE_MTIME_MIN:-20}"   # newer than this => run is LIVE
LIVE_QUIET_SEC="${LIVE_QUIET_SEC:-300}"  # never delete files younger than this
DATA_DIR="${DATA_DIR:-data}"
PROTECT_FILE="${PROTECT_FILE:-BRAINSTORM/018_dji9443_hover_convergence_campaign/vtk_protect_list.txt}"

APPLY=false
SKIP_LIVE=false
ONLY=""
LAST_RESORT=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --apply)     APPLY=true; shift ;;
        --skip-live) SKIP_LIVE=true; shift ;;
        --emergency) shift ;;   # accepted for compatibility; live sweeping is now default
        --only)      ONLY="$2"; shift 2 ;;
        --keep)      KEEP_STEPS="$2"; shift 2 ;;
        --last-resort) LAST_RESORT=true; shift ;;
        -h|--help)   sed -n '2,44p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

[[ -d "$DATA_DIR" ]] || { echo "no $DATA_DIR/ here -- run from the checkout root" >&2; exit 2; }

# ------------------------------------------------------- keep-value validation
# --keep is only ever passed under time pressure (the >300 G escalation tier),
# and bad values are silently catastrophic: a KEEP_STEPS of 0, "" or negative
# makes the awk top-list below retain NOTHING, deleting 100% of a run's VTK with
# no warning.  Verified against a synthetic fixture 2026-08-24.  Validated here,
# after arg parsing, so it covers the KEEP_STEPS environment variable too.
HARD_FLOOR=36    # absolute floor; last-resort tier (Ryan, 2026-08-24)
SOFT_FLOOR=72    # >300 G escalation tier; below this needs --last-resort

if [[ ! "$KEEP_STEPS" =~ ^[0-9]+$ ]]; then
    echo "FATAL: --keep must be a positive integer, got '$KEEP_STEPS'" >&2
    exit 4
fi
if (( KEEP_STEPS < HARD_FLOOR )); then
    echo "FATAL: --keep $KEEP_STEPS is below the hard floor of $HARD_FLOOR steps" >&2
    exit 4
fi
if (( KEEP_STEPS < SOFT_FLOOR )) && ! $LAST_RESORT; then
    echo "FATAL: --keep $KEEP_STEPS is below the ${SOFT_FLOOR}-step escalation floor." >&2
    echo "       Pass --last-resort only if the 400 G cap is still breached at 72." >&2
    exit 4
fi

# An --only that matches nothing must not look like a successful empty sweep:
# during incremental escalation that reads as "swept, freed nothing" when the
# truth is "never looked at it".
if [[ -n "$ONLY" && ! -d "$DATA_DIR/$ONLY" ]]; then
    echo "FATAL: --only '$ONLY' matches no directory under $DATA_DIR/" >&2
    exit 5
fi

# ---------------------------------------------------------------- protect list
# Read fresh on every invocation so it can be edited mid-flight.  This script
# only ever READS this file.
protected_runs=""
if [[ -f "$PROTECT_FILE" ]]; then
    protected_runs="$(sed -e 's/#.*//' -e 's/[[:space:]]//g' "$PROTECT_FILE" | grep -v '^$' || true)"
else
    echo "FATAL: protect list not found at $PROTECT_FILE -- refusing to sweep" >&2
    exit 3
fi

is_protected() {
    local run="$1" pat
    while IFS= read -r pat; do
        [[ -z "$pat" ]] && continue
        # shellcheck disable=SC2053  -- glob match is intentional
        [[ "$run" == $pat ]] && return 0
    done <<< "$protected_runs"
    return 1
}

# ------------------------------------------------------------------- liveness
# Belt-and-braces: mtime OR a matching squeue job name.  Either signal => LIVE.
# Slurm truncates job names, so the name test is a normalised substring test,
# which errs toward calling things live.
job_names=""
if command -v squeue >/dev/null 2>&1; then
    job_names="$(squeue -u "$USER" -h -o '%j' 2>/dev/null | tr '-' '_' | sed 's/^fp_//' || true)"
fi

is_live() {
    local run="$1" dir="$2" jn
    [[ -n "$(find "$dir" -type f -mmin "-$LIVE_MTIME_MIN" -print -quit 2>/dev/null)" ]] && return 0
    while IFS= read -r jn; do
        [[ -z "$jn" ]] && continue
        [[ "$run" == *"$jn"* ]] && return 0
    done <<< "$job_names"
    return 1
}

# ---------------------------------------------------------------------- sweep
NOW="$(date +%s)"
total_freed=0
manifest="$(mktemp)"
trap 'rm -f "$manifest"' EXIT

echo "# p018_vtk_sweeper  keep=$KEEP_STEPS steps  apply=$APPLY  skip_live=$SKIP_LIVE"
echo "# protect list: $PROTECT_FILE"

for dir in "$DATA_DIR"/*/; do
    run="$(basename "$dir")"
    [[ -n "$ONLY" && "$run" != "$ONLY" ]] && continue

    if is_protected "$run"; then
        echo "PROTECTED $run"
        continue
    fi

    live=false
    is_live "$run" "$dir" && live=true

    if $live && $SKIP_LIVE; then
        echo "LIVE-SKIP $run"
        continue
    fi

    # One find per run; awk does the step algebra and emits deletion candidates.
    : > "$manifest"
    find "$dir" -type f -name '*.vt[upsm]' -printf '%s\t%T@\t%p\n' 2>/dev/null |
    awk -F'\t' -v run="$run" -v keep="$KEEP_STEPS" -v live="$live" \
               -v now="$NOW" -v quiet="$LIVE_QUIET_SEC" '
    {
        size = $1; mt = int($2); path = $3
        base = path; sub(/.*\//, "", base)
        stem = base; sub(/\.[^.]*$/, "", stem)      # drop extension
        step = stem; sub(/.*\./, "", step)          # step = last dot-field
        if (step !~ /^[0-9]+$/) next

        n++; P[n] = path; S[n] = step; Z[n] = size; M[n] = mt
        seen[step] = 1

        # Mark the four roles that together make a step restartable.  Separate
        # arrays rather than a bitmask: or() is a gawk extension.
        if (base == run "_body1." step ".vtu")                 hb[step]  = 1
        else if (base == run "_wake1.1." step ".vts")          hw1[step] = 1
        else if (base == run "_wake1.2." step ".vts")          hw2[step] = 1
        else if (base == run "_wake1_particles." step ".vtp")  hp[step]  = 1
    }
    END {
        # Select the `keep` highest restartable steps by insertion into a small
        # descending top-list -- avoids asort(), also a gawk extension.
        nk = 0
        nc = 0
        for (s in seen) {
            if (!(hb[s] && hw1[s] && hw2[s] && hp[s])) continue
            nc++
            v = s + 0
            pos = nk + 1
            for (i = 1; i <= nk; i++) if (v > top[i]) { pos = i; break }
            if (pos > keep) continue
            for (i = (nk < keep ? nk : keep - 1); i >= pos; i--) top[i+1] = top[i]
            top[pos] = v
            if (nk < keep) nk++
        }
        if (nc == 0) { print "NO-RESTART-SET"; exit }

        kept_list = ""
        for (i = nk; i >= 1; i--) {
            kept[top[i] ""] = 1
            kept_list = kept_list (kept_list == "" ? "" : ",") top[i]
        }
        printf "KEEP\t%d\t%s\n", nc, kept_list

        for (i = 1; i <= n; i++) {
            if (kept[S[i] ""] == 1) continue
            if (live == "true" && now - M[i] < quiet) continue
            printf "DEL\t%d\t%s\n", Z[i], P[i]
        }
    }' > "$manifest" || true

    if [[ ! -s "$manifest" ]]; then
        echo "EMPTY $run"
        continue
    fi

    if grep -q '^NO-RESTART-SET' "$manifest"; then
        echo "NO-RESTART-SET $run  (skipped, not guessing)"
        continue
    fi

    ncomplete="$(awk -F'\t' '/^KEEP/{print $2}' "$manifest")"
    keptsteps="$(awk -F'\t' '/^KEEP/{print $3}' "$manifest")"
    nfiles="$(grep -c '^DEL' "$manifest" || true)"
    bytes="$(awk -F'\t' '/^DEL/{s+=$2} END{print s+0}' "$manifest")"
    mb=$(( bytes / 1048576 ))

    state="CLOSED"; $live && state="LIVE"
    echo "SWEEP $run  state=$state  restartable_steps=$ncomplete  keeping=[$keptsteps]  delete=${nfiles}files/${mb}MB"

    if $APPLY && [[ "${nfiles:-0}" -gt 0 ]]; then
        awk -F'\t' '/^DEL/{print $3}' "$manifest" | tr '\n' '\0' | xargs -0 -r rm -f
        echo "APPLIED $run  freed=${mb}MB"
    fi
    total_freed=$(( total_freed + bytes ))
done

echo "TOTAL_FREED_MB=$(( total_freed / 1048576 ))"
echo "DF_AFTER=$(df -BG "$HOME" | tail -1 | awk '{print $3}')"
echo "SWEEP_MODE=$( $APPLY && echo APPLY || echo DRY-RUN )"
