#!/usr/bin/env bash
#
# p018_vtk_sweeper.sh -- restart-set-aware VTK sweeper for the 200 G home budget.
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
#     leave the newest timesteps (keep count: $KEEP_STEPS, 36 since 2026-08-20) -- no need to wait for the run to
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
#
# Run from the top level of the FLOWPanel.jl checkout.

set -euo pipefail

KEEP_STEPS="${KEEP_STEPS:-36}"      # restartable steps retained per run (10->36, Ryan 2026-08-20)
LIVE_MTIME_MIN="${LIVE_MTIME_MIN:-20}"   # newer than this => run is LIVE
LIVE_QUIET_SEC="${LIVE_QUIET_SEC:-300}"  # never delete files younger than this
DATA_DIR="${DATA_DIR:-data}"
PROTECT_FILE="${PROTECT_FILE:-BRAINSTORM/018_dji9443_hover_convergence_campaign/vtk_protect_list.txt}"

APPLY=false
SKIP_LIVE=false
ONLY=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --apply)     APPLY=true; shift ;;
        --skip-live) SKIP_LIVE=true; shift ;;
        --emergency) shift ;;   # accepted for compatibility; live sweeping is now default
        --only)      ONLY="$2"; shift 2 ;;
        --keep)      KEEP_STEPS="$2"; shift 2 ;;
        -h|--help)   sed -n '2,40p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

[[ -d "$DATA_DIR" ]] || { echo "no $DATA_DIR/ here -- run from the checkout root" >&2; exit 2; }

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
