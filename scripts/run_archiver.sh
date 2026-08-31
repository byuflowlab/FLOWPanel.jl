#!/usr/bin/env bash
#
# run_archiver.sh -- archive-first reclaim for the 400 G /home/rander39 budget.
#
# Finished runs are TARRED TO ARCHIVE, not deleted.  This is the primary reclaim
# mechanism; scripts/p018_vtk_sweeper.sh is now the fallback and applies only to
# runs that are still writing (a live run cannot be consistently tarred).
#
# Policy encoded here (Ryan, 2026-08-27):
#
#   * FINISHED = no matching Slurm job AND newest file mtime older than
#     $QUIET_HOURS (24 h).  Both signals are required.  A live run is never
#     touched by this script.
#   * The TARBALL IS SELF-CONTAINED: the whole run directory goes in, including
#     CSVs, *_case_metadata.toml, *.metadata.toml, *.pvd and monitors/.
#   * On /home ONLY PARAVIEW FILES ARE REMOVED, and only those OUTSIDE the
#     newest $KEEP_STEPS restartable steps (default 5, Ryan 2026-08-27).  An
#     archived run therefore stays warm-startable straight off /home -- no
#     unpacking a tarball to continue it.  The CSV/TOML/pvd/monitors residue
#     stays exactly where it was too: it is small, it is the scientific output,
#     and /home is backed up daily where the archive is NOT.
#   * NOTHING IS DELETED UNTIL THE TARBALL HAS BEEN WRITTEN AND VERIFIED by
#     reading it back.  Verification failure aborts that run and deletes nothing.
#
# Archive tier facts that this script depends on (https://rc.byu.edu/wiki/?id=Storage):
#   /nobackup/archive/usr/$USER -- 20 TiB but only 1 MILLION FILES, no backup,
#   no snapshots, no auto-delete, and mounted on LOGIN NODES ONLY.  A Slurm job
#   cannot read or write it.  One tarball per run is what keeps per-timestep VTK
#   from eating the inode budget.  It is also slow by design.
#
# Dry-run by default.  Pass --apply to actually write and delete.
#
# There are SEVERAL FLOWPanel checkouts under $HOME, each with its own data/
# tree, and all of them are in scope.  Archives are filed per source checkout so
# the origin of a tarball is never in doubt:
#
#   $ARCHIVE_DIR/<checkout-slug>/<run>.tar.zst
#
# where <checkout-slug> is the checkout path relative to $HOME with '/' -> '_'
# (e.g. projects_FLOWPanel.jl, FLOWPanel-018-gpu-h200,
# flowpanel-021_FLOWPanel.jl).  Two clones may hold a run of the same name; the
# slug directory keeps them apart.  Provenance is recorded THREE ways so it
# survives any one of them being lost:
#
#   1. the slug directory the tarball sits in;
#   2. ARCHIVE_SOURCE.txt INSIDE the tarball -- absolute source path, host, git
#      commit, date, sizes -- so a tarball moved elsewhere still says where it
#      came from;
#   3. $ARCHIVE_DIR/INDEX.tsv, one append-only row per archived run, so
#      "which clone did this come from?" is one grep across everything.
#
# Usage:
#   scripts/run_archiver.sh [--apply] [--root DIR]... [--all-checkouts]
#                           [--only RUN] [--keep N] [--quiet-hours N]
#                           [--compress zstd|gzip|none] [--fast-verify]
#   scripts/run_archiver.sh --restore RUN [--root DIR] [--apply]
#   scripts/run_archiver.sh --resume-delete RUN [--root DIR] [--apply]
#
# With no --root, the current checkout (cwd) is the only one processed -- the
# same behaviour as before this option existed.
#
# Concurrency: workers on DIFFERENT --root checkouts may run in parallel
# (they are disjoint in every path they touch); locks live on the SHARED
# filesystem under $LOCKROOT (default $HOME/.cache/flowpanel/locks) so they
# hold across login nodes.  No lock auto-expires: a stale lock after a crash
# is removed by a HUMAN, and the error message names the exact rm -rf.
#
# Exit codes: 2 bad usage, 3 no protect list, 5 unmatched --only/--restore,
#             6 lock conflict: a sweeper holds the shared lock, or a checkout
#               was skipped because another archiver holds its lock,
#             7 archive dir
#             unreachable or not writable (you are probably on a compute node),
#             8 a tarball failed verification (those runs deleted NOTHING) or a
#               --restore left members missing on disk (tarball kept intact),
#             9 ARCHIVED-STALE runs found -- ask Ryan, do not proceed alone.
#
# Run from the top level of the FLOWPanel.jl checkout.

set -euo pipefail

ARCHIVE_DIR="${ARCHIVE_DIR:-/nobackup/archive/usr/$USER/FLOWPanel_runs}"
DATA_DIR="${DATA_DIR:-data}"
CHECKOUT_GLOBS="${CHECKOUT_GLOBS:-$HOME/* $HOME/*/*}"   # searched by --all-checkouts
PROTECT_FILE="${PROTECT_FILE:-BRAINSTORM/018_dji9443_hover_convergence_campaign/vtk_protect_list.txt}"
QUIET_HOURS="${QUIET_HOURS:-24}"        # finished => newest mtime older than this
RECENT_HOT_HOURS="${RECENT_HOT_HOURS:-2}"  # below this, not even Ryan may approve
KEEP_STEPS="${KEEP_STEPS:-5}"           # restartable steps left on /home after archiving
COMPRESS="${COMPRESS:-zstd}"            # zstd | gzip | none
# Archive is mounted on LOGIN NODES ONLY, so compression necessarily runs on a
# shared interactive node.  Cap the thread count -- `zstd -T0` would take every
# core on a machine other people are using.
ZSTD_THREADS="${ZSTD_THREADS:-4}"

APPLY=false
ONLY=""
RESTORE=""
RESUME_DELETE=""
FAST_VERIFY=false
ALL_CHECKOUTS=false
INCLUDE_RECENT=false
ROOTS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --apply)          APPLY=true; shift ;;
        --only)           ONLY="${2:-}"; shift 2 ;;
        --root)           ROOTS+=("${2:-}"); shift 2 ;;
        --all-checkouts)  ALL_CHECKOUTS=true; shift ;;
        --include-recent) INCLUDE_RECENT=true; shift ;;
        --keep)           KEEP_STEPS="${2:-}"; shift 2 ;;
        --restore)        RESTORE="${2:-}"; shift 2 ;;
        --resume-delete)  RESUME_DELETE="${2:-}"; shift 2 ;;
        --quiet-hours)    QUIET_HOURS="${2:-}"; shift 2 ;;
        --compress)       COMPRESS="${2:-}"; shift 2 ;;
        --fast-verify)    FAST_VERIFY=true; shift ;;
        -h|--help)        sed -n '2,40p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

# ------------------------------------------------------------------ locking
# Cross-NODE mutual exclusion on the SHARED filesystem.  The old flock lived
# in ${TMPDIR:-/tmp}, which is node-local: a sweeper on login1 and an archiver
# on login2 shared no lock at all.  The primitive is mkdir(2), atomic on both
# VAST NFS and Lustre; flock(2) is NOT used because on Lustre mounted with
# -o localflock it SUCCEEDS while protecting nothing.
#
# Scheme (mirrored in scripts/p018_vtk_sweeper.sh):
#   * each archiver publishes a READER marker, then re-checks that no sweeper
#     holds sweeper.excl.  Publish-then-check on BOTH sides is what makes the
#     race safe: whichever order two parties act in, at least one sees the
#     other and backs off;
#   * the sweeper takes sweeper.excl exclusively and refuses to start while
#     any reader exists.  A sweep and an archive must never overlap: a sweeper
#     deleting VTK from a run that has no verified tarball yet is a permanent,
#     unarchived loss, and CHANGED-DURING-ARCHIVE cannot catch it (deleting
#     old files does not raise the tree's newest mtime);
#   * archivers exclude EACH OTHER per checkout via checkout.<slug> locks in
#     the main loop: two archivers on the SAME checkout can overwrite each
#     other's tarball after one has already stripped /home, while archivers on
#     DIFFERENT checkouts are disjoint in every path they touch and may run
#     concurrently (one worker per checkout is the supported parallelism).
#
# NO lock auto-expires.  A worker that goes quiet is usually hung on NFS or
# Lustre (D-state, routinely minutes), not dead; stealing its lock puts two
# workers on one run, which is the catastrophic case.  After a SIGKILL or a
# node crash a HUMAN removes the stale lock with the rm -rf the message names.
LOCKROOT="${LOCKROOT:-$HOME/.cache/flowpanel/locks}"
READERS="$LOCKROOT/readers"
mkdir -p "$READERS"
LOCK_TOKEN="$(hostname -s 2>/dev/null || echo unknown).$$.$(date +%s)"
READER_MARK="$READERS/archiver.$LOCK_TOKEN"
CHECKOUT_LOCKS=()
TMPD=""
release_locks() {
    local d
    for d in ${CHECKOUT_LOCKS[@]+"${CHECKOUT_LOCKS[@]}"}; do rm -rf "$d"; done
    rm -rf "$READER_MARK"
    [[ -n "$TMPD" ]] && rm -rf "$TMPD"
    return 0
}
trap release_locks EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
mkdir "$READER_MARK"
printf 'host %s\npid %s\nstarted_utc %s\n' \
    "$(hostname -s 2>/dev/null || echo unknown)" "$$" \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$READER_MARK/owner" 2>/dev/null || true
if [[ -d "$LOCKROOT/sweeper.excl" ]]; then
    echo "FATAL: a sweeper holds $LOCKROOT/sweeper.excl -- an archive and a sweep must never overlap." >&2
    sed 's/^/       sweeper: /' "$LOCKROOT/sweeper.excl/owner" >&2 2>/dev/null || true
    echo "       If that sweeper is dead (host/pid gone), a human may: rm -rf $LOCKROOT/sweeper.excl" >&2
    exit 6
fi

# Per-checkout exclusive lock among archivers.  Returns 1 (silently) when
# another archiver holds the checkout.  Released at the end of that checkout's
# loop iteration, and again -- harmlessly -- by the EXIT trap.
acquire_checkout_lock() {
    CKLOCK="$LOCKROOT/checkout.$1"
    mkdir "$CKLOCK" 2>/dev/null || return 1
    CHECKOUT_LOCKS+=("$CKLOCK")
    printf 'host %s\npid %s\nstarted_utc %s\n' \
        "$(hostname -s 2>/dev/null || echo unknown)" "$$" \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$CKLOCK/owner" 2>/dev/null || true
    return 0
}
checkout_locked_msg() {
    echo "CHECKOUT-LOCKED $1  -- another archiver holds $LOCKROOT/checkout.$1; skipping"
    sed 's/^/    holder: /' "$LOCKROOT/checkout.$1/owner" 2>/dev/null || true
    echo "    if that archiver is dead (host/pid gone), a human may: rm -rf $LOCKROOT/checkout.$1"
}

TMPD="$(mktemp -d)"

# ------------------------------------------------------------- checkout roots
# A directory counts as a FLOWPanel checkout only if its Project.toml declares
# the FLOWPanel package UUID and it has a data/ directory.
#
# Weaker tests were tried and are wrong:
#   * Project.toml + data/ alone matched 19 directories on the cluster, of which
#     8 were per-checkout `test/` environments (a [deps]-only Project.toml plus a
#     couple of scratch outputs) and 3 were OTHER packages entirely --
#     FastMultipole-026/046 and FLOWVPM-046.  None of those hold FLOWPanel runs.
#   * matching name = "FLOWPanel" would work today but breaks on a rename; the
#     UUID is the package's canonical identity.
FLOWPANEL_UUID="${FLOWPANEL_UUID:-6be8c882-484d-4309-b349-d23112750151}"
is_checkout() {
    [[ -f "$1/Project.toml" && -d "$1/data" ]] || return 1
    grep -q "$FLOWPANEL_UUID" "$1/Project.toml" 2>/dev/null
}

if $ALL_CHECKOUTS; then
    [[ ${#ROOTS[@]} -eq 0 ]] || { echo "FATAL: --all-checkouts and --root are mutually exclusive" >&2; exit 2; }
    for c in $CHECKOUT_GLOBS; do
        [[ -d "$c" ]] || continue
        is_checkout "$c" && ROOTS+=("$c")
    done
    [[ ${#ROOTS[@]} -gt 0 ]] || { echo "FATAL: --all-checkouts found no FLOWPanel checkout under $CHECKOUT_GLOBS" >&2; exit 2; }
fi

if [[ ${#ROOTS[@]} -eq 0 ]]; then
    # Default: this checkout only, exactly as before --root existed.
    [[ -d "$DATA_DIR" ]] || { echo "no $DATA_DIR/ here -- run from the checkout root, or pass --root" >&2; exit 2; }
    ROOTS=(".")
fi

for r in "${ROOTS[@]}"; do
    [[ -d "$r/data" ]] || { echo "FATAL: --root '$r' has no data/ directory" >&2; exit 2; }
done

# Slug identifying which checkout a tarball came from: path relative to $HOME
# with '/' -> '_'.  Basename alone is NOT enough -- projects/FLOWPanel.jl and
# flowpanel-021/FLOWPanel.jl would collide, and so would their run names.
slug_for() {
    local abs
    abs="$(cd "$1" && pwd -P)"
    local rel="${abs#"$HOME"/}"
    [[ "$rel" == "$abs" ]] && rel="${abs#/}"    # outside $HOME: use the full path
    echo "${rel//\//_}"
}

if [[ ! "$QUIET_HOURS" =~ ^[0-9]+$ ]] || (( QUIET_HOURS < 1 )); then
    echo "FATAL: --quiet-hours must be a positive integer, got '$QUIET_HOURS'" >&2
    exit 2
fi

# Unlike the sweeper's --keep this has no hard floor: the tarball is verified
# before anything is removed, so a small (even zero) window is recoverable
# rather than catastrophic.  0 means "strip all VTK; restarting needs a
# --restore first" -- legal, but say so out loud rather than typing it by
# accident.
if [[ ! "$KEEP_STEPS" =~ ^[0-9]+$ ]]; then
    echo "FATAL: --keep must be a non-negative integer, got '$KEEP_STEPS'" >&2
    exit 2
fi

case "$COMPRESS" in
    zstd) command -v zstd >/dev/null 2>&1 || COMPRESS=gzip ;;
    gzip|none) ;;
    *) echo "FATAL: --compress must be zstd, gzip or none (got '$COMPRESS')" >&2; exit 2 ;;
esac

case "$COMPRESS" in
    zstd) EXT="tar.zst"; COMP_CMD=(zstd -3 -q -T"$ZSTD_THREADS" -c); DECOMP_CMD=(zstd -dcq) ;;
    gzip) EXT="tar.gz";  COMP_CMD=(gzip -1 -c);        DECOMP_CMD=(gzip -dc) ;;
    none) EXT="tar";     COMP_CMD=(cat);               DECOMP_CMD=(cat) ;;
esac

# stat(1) is not portable; resolve the "size path" form once.
if stat -c '%s %n' . >/dev/null 2>&1; then
    STAT_FMT=(stat -c '%s %n')          # GNU
    STAT_MTIME=(stat -c '%Y')
else
    STAT_FMT=(stat -f '%z %N')          # BSD/macOS
    STAT_MTIME=(stat -f '%m')
fi

# Emits "<size> <path>" for NUL-separated paths on stdin, and NEVER fails.
# Both guards are load-bearing:
#   * GNU xargs runs the command once even on EMPTY input, so a run with no
#     ParaView files made `stat` run with no arguments -> exit 123 -> under
#     `set -e` the whole --all-checkouts pass died partway through, silently
#     skipping every later checkout (observed 2026-08-27 on FLOWPanel-052-h200).
#   * a live job deleting a file between find and stat produces the same 123.
stat_sizes() { xargs -0 "${STAT_FMT[@]}" 2>/dev/null || true; }

# ParaView output -- the ONLY class this script ever removes from /home.
# Deliberately excludes *.pvd, *.csv and *.toml, and never descends monitors/.
pv_find() {
    find "$1" -type f -not -path '*/monitors/*' \
        \( -name '*.vt[upsmri]' -o -name '*.vtk' -o -name '*.pvt[up]' \) "${@:2}"
}

# ------------------------------------------------------- restart-window select
# Emits, for one run directory, the VTK that may be removed from /home once the
# tarball is verified -- i.e. everything OUTSIDE the newest $KEEP_STEPS
# RESTARTABLE steps.  Step S is restartable only if all four paths the warmstart
# loader reads exist (src/FLOWPanel_warmstart.jl), the same definition
# scripts/p018_vtk_sweeper.sh uses:
#     <run>_body1/<run>_body1.<S>.vtu
#     <run>_wake1/<run>_wake1.1.<S>.vts   and   .2.<S>.vts
#     <run>_wake1_particles/<run>_wake1_particles.<S>.vtp
# ALL files at a kept step are kept -- .vtm indices together with their pieces --
# so the retained state is both warm-startable and ParaView-openable.  Leaving an
# index stub without its pieces is what killed job 13036477 on 2026-08-04.
#
# Output: one "KEEP\t<n_restartable>\t<comma list>" line, then "DEL\t<size>\t<path>".
# A run with NO restartable step yields KEEP with an empty list and every file as
# DEL: there is no restart to preserve, and unlike the sweeper's refusal to guess,
# here the verified tarball still holds all of it.
select_deletable() {
    local dir="$1" run="$2" keep="$3"
    pv_find "$1" -print0 2>/dev/null | stat_sizes |
    awk -v run="$run" -v keep="$keep" '
    {
        size = $1; path = $0; sub(/^[0-9]+ /, "", path)
        base = path; sub(/.*\//, "", base)
        stem = base; sub(/\.[^.]*$/, "", stem)      # drop extension
        step = stem; sub(/.*\./, "", step)          # step = last dot-field
        if (step !~ /^[0-9]+$/) next

        n++; P[n] = path; S[n] = step; Z[n] = size
        seen[step] = 1

        if (base == run "_body1." step ".vtu")                 hb[step]  = 1
        else if (base == run "_wake1.1." step ".vts")          hw1[step] = 1
        else if (base == run "_wake1.2." step ".vts")          hw2[step] = 1
        else if (base == run "_wake1_particles." step ".vtp")  hp[step]  = 1
    }
    END {
        # keep highest `keep` restartable steps, via a small descending top-list
        nk = 0; nc = 0
        for (s in seen) {
            if (!(hb[s] && hw1[s] && hw2[s] && hp[s])) continue
            nc++
            if (keep == 0) continue
            v = s + 0
            pos = nk + 1
            for (i = 1; i <= nk; i++) if (v > top[i]) { pos = i; break }
            if (pos > keep) continue
            for (i = (nk < keep ? nk : keep - 1); i >= pos; i--) top[i+1] = top[i]
            top[pos] = v
            if (nk < keep) nk++
        }
        kept_list = ""
        for (i = nk; i >= 1; i--) {
            kept[top[i] ""] = 1
            kept_list = kept_list (kept_list == "" ? "" : ",") top[i]
        }
        printf "KEEP\t%d\t%s\n", nc, kept_list
        for (i = 1; i <= n; i++) {
            if (kept[S[i] ""] == 1) continue
            printf "DEL\t%d\t%s\n", Z[i], P[i]
        }
    }'
}

# ---------------------------------------------------------------- protect list
# One protect list governs ALL checkouts: run names are meaningful campaign-wide,
# not per-clone.  Resolve to an absolute path now so it still points at the same
# file when we are working inside a different checkout.
if [[ -f "$PROTECT_FILE" ]]; then
    PROTECT_FILE="$(cd "$(dirname "$PROTECT_FILE")" && pwd -P)/$(basename "$PROTECT_FILE")"
fi
protected_runs=""
if [[ -f "$PROTECT_FILE" ]]; then
    protected_runs="$(sed -e 's/#.*//' -e 's/[[:space:]]//g' "$PROTECT_FILE" | grep -v '^$' || true)"
else
    echo "FATAL: protect list not found at $PROTECT_FILE -- refusing to archive" >&2
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
# Same belt-and-braces test as the sweeper: mtime OR a matching squeue job name.
# Slurm truncates and mangles names (fp-018-L1visc vs p018_L1_visc), so the name
# test is a normalised substring test and errs toward calling things live.
job_names=""
squeue_ok=false
if command -v squeue >/dev/null 2>&1; then
    if job_names="$(squeue -u "$USER" -h -o '%j' 2>/dev/null | tr '-' '_' | sed 's/^fp_//')"; then
        squeue_ok=true
    fi
fi

# The two liveness signals are reported SEPARATELY rather than OR-ed, because
# they disagree often and in both directions (measured on the cluster
# 2026-08-27):
#
#   * FALSE POSITIVE: job `fp-018-csarc_n5_nt144_l2p4` normalises to a string
#     that is a substring of SEVEN different run directories (_s2gpu, _s3gpu,
#     _sortdiag, _sortdiag2, _bisect_ell8, and the bare name in two checkouts).
#     Six of them are not being written at all.
#   * FALSE NEGATIVE: job `fp-018gpu-n2_nt72_l3p0` normalises to
#     `018gpu_n2_nt72_l3p0`, which is NOT a substring of the directory it is
#     actively writing, `p018_csarc_n2_nt72_l3p0`.
#
# So the queue is treated as advisory and mtime as the real safety signal.
# A run with a queue match is LIVE and untouchable.  A run with no queue match
# but recent activity is RECENT: it looks finished, but the quiet window has not
# elapsed, so it needs Ryan's say-so (Ryan, 2026-08-27).

# Echoes the matching job name, or nothing.
queue_match() {
    local run="$1" jn
    while IFS= read -r jn; do
        [[ -z "$jn" ]] && continue
        [[ "$run" == *"$jn"* ]] && { echo "$jn"; return 0; }
    done <<< "$job_names"
    return 1
}

# Newest mtime under a run directory, excluding this script's own breadcrumbs --
# they are written by us, not the simulation, and counting them would make a run
# we just touched look active for the next $QUIET_HOURS.  That once left a run
# whose archive failed un-retryable for a day.
# `|| true` is load-bearing, not defensive clutter.  These pipelines stat tens of
# thousands of files in a directory a live job may be writing; a file that
# disappears between find and stat makes xargs exit 123, and under `set -e` that
# aborted an entire --all-checkouts pass partway through (observed 2026-08-27 on
# a run with 23,609 VTK files).  `head -1` closing the pipe early can also raise
# SIGPIPE upstream.  Every stat-over-find pipeline in this script is guarded.
newest_mtime() {
    { find "$1" -type f \
        ! -name 'ARCHIVED.txt' ! -name 'RESTORED.txt' ! -name 'ARCHIVE_SOURCE.txt' \
        -print0 2>/dev/null | xargs -0 "${STAT_MTIME[@]}" 2>/dev/null | sort -rn | head -1; } || true
}

quiet_hours_of() {
    local m; m="$(newest_mtime "$1")"
    [[ -z "$m" ]] && { echo 99999; return; }
    echo $(( ( $(date +%s) - m ) / 3600 ))
}

# ---------------------------------------------------------------- verification
# Reads the finished tarball back and checks BOTH the member path set and the
# total content byte count against the source manifest.  One decompression pass:
# `tar -xvO` writes contents to stdout and member names to stderr.
#
# manifest format: "<size>\t<relpath>" sorted by relpath (paths relative to
# $DATA_DIR, i.e. they start with "<run>/").
verify_tarball() {
    local tarball="$1" manifest="$2" run="$3"
    local names bytes_seen bytes_want

    names="$TMPD/names"
    if $FAST_VERIFY; then
        "${DECOMP_CMD[@]}" < "$tarball" | tar -tf - > "$names" 2>/dev/null || {
            echo "VERIFY-FAIL $run  (tar listing failed)"; rm -f "$names"; return 1; }
        bytes_seen=""
    else
        bytes_seen="$("${DECOMP_CMD[@]}" < "$tarball" | tar -xvOf - 2>"$names" | wc -c | tr -cd '0-9')"
    fi

    # Compare member path sets.  Normalise the two tar dialects: bsdtar prefixes
    # extracted names with "x " where GNU tar does not, and either may emit a
    # leading "./".  Drop directory entries (trailing /).
    local want got
    want="$TMPD/want"; got="$TMPD/got"
    cut -f2- "$manifest" | sort > "$want"
    sed -e 's|^x ||' -e 's|^\./||' -e '/\/$/d' "$names" | grep -v '^$' | sort -u > "$got"

    if ! cmp -s "$want" "$got"; then
        echo "VERIFY-FAIL $run  (member list differs from source: $(comm -23 "$want" "$got" | wc -l | tr -cd '0-9') missing, $(comm -13 "$want" "$got" | wc -l | tr -cd '0-9') extra)"
        
        return 1
    fi

    if ! $FAST_VERIFY; then
        bytes_want="$(awk -F'\t' '{s+=$1} END{print s+0}' "$manifest")"
        if [[ "$bytes_seen" != "$bytes_want" ]]; then
            echo "VERIFY-FAIL $run  (content bytes $bytes_seen != source $bytes_want)"
        
            return 1
        fi
    fi

        
    return 0
}

# Build "<size>\t<relpath>" for every regular file under $DATA_DIR/$run.
build_manifest() {
    local run="$1" out="$2"
    ( cd "$DATA_DIR" && find "$run" -type f -print0 2>/dev/null | stat_sizes ) |
        sed -e 's/^\([0-9][0-9]*\) /\1\t/' | sort -t"$(printf '\t')" -k2 > "$out"
}

# ------------------------------------------------- single-checkout resolution
# --restore and --resume-delete act on exactly one checkout.  With one --root
# that is unambiguous; with several (or --all-checkouts) it is not, and guessing
# would restore into or delete from the wrong clone.
single_root() {
    if [[ ${#ROOTS[@]} -ne 1 ]]; then
        echo "FATAL: $1 needs exactly one checkout; got ${#ROOTS[@]}. Pass --root DIR." >&2
        exit 2
    fi
    echo "${ROOTS[0]}"
}

# Locate a run's tarball: the slug subdir first, then the flat legacy layout
# used before archives were filed per checkout.
find_tarball() {
    local sub="$1" run="$2" e
    for e in tar.zst tar.gz tar; do
        [[ -f "$sub/$run.$e" ]]          && { echo "$sub/$run.$e"; return 0; }
    done
    for e in tar.zst tar.gz tar; do
        [[ -f "$ARCHIVE_DIR/$run.$e" ]]  && { echo "$ARCHIVE_DIR/$run.$e"; return 0; }
    done
    return 1
}

# ------------------------------------------------------------------- restore
if [[ -n "$RESTORE" ]]; then
    root="$(single_root --restore)"
    DATA_DIR="$root/data"
    SLUG="$(slug_for "$root")"
    # A restore writes into data/ that a concurrent archiver may be tarring.
    acquire_checkout_lock "$SLUG" || { checkout_locked_msg "$SLUG"; exit 6; }
    tarball="$(find_tarball "$ARCHIVE_DIR/$SLUG" "$RESTORE")" || {
        echo "FATAL: no tarball for '$RESTORE' under $ARCHIVE_DIR/$SLUG (checkout $(cd "$root" && pwd -P))" >&2
        echo "       If it came from a different clone, pass that clone's --root." >&2
        exit 5; }

    case "$tarball" in
        *.tar.zst) DECOMP_CMD=(zstd -dcq) ;;
        *.tar.gz)  DECOMP_CMD=(gzip -dc) ;;
        *)         DECOMP_CMD=(cat) ;;
    esac

    echo "# run_archiver RESTORE  $tarball -> $DATA_DIR/$RESTORE  (checkout $SLUG)"
    if ! $APPLY; then
        echo "RESTORE-PLAN $RESTORE  members=$("${DECOMP_CMD[@]}" < "$tarball" | tar -tf - | wc -l | tr -cd '0-9')"
        echo "ARCHIVE_MODE=DRY-RUN"
        exit 0
    fi
    # -k refuses to overwrite existing files, so a restore can never clobber
    # newer /home state (notably the CSV/monitor residue left behind on archive).
    # The `|| true` is unavoidable: -k makes tar exit nonzero on the EXPECTED
    # collisions with that residue, indistinguishable from a real failure.  So
    # tar's exit status proves nothing here -- completeness is verified below
    # by comparing the member list against what actually landed on disk.
    "${DECOMP_CMD[@]}" < "$tarball" | tar -xkf - -C "$DATA_DIR" 2>/dev/null || true

    # A truncated tarball, a zstd checksum failure, or ENOSPC would otherwise
    # be swallowed and reported RESTORED.  Verify the same way verify_tarball
    # does -- one decompression pass, member names on stderr, content bytes on
    # stdout -- then require every member to exist on disk AND the members'
    # total on-disk size to match the archived content.  The byte check matters:
    # a truncated file left by an interrupted earlier restore EXISTS, and -k
    # would skip it forever; existence alone would call that restore complete.
    # (pipefail makes the substitution fail if zstd or tar dies mid-stream.)
    rnames="$TMPD/rnames"
    bytes_want="$("${DECOMP_CMD[@]}" < "$tarball" | tar -xvOf - 2>"$rnames" | wc -c | tr -cd '0-9')" || {
        echo "RESTORE-INCOMPLETE $RESTORE  (tarball unreadable: $tarball)"
        echo "         nothing certified restored -- fix the tarball before trusting this run"
        exit 8; }
    sed -e 's|^x ||' -e 's|^\./||' -e '/\/$/d' "$rnames" | grep -v '^$' | sort -u > "$TMPD/members"
    ( cd "$DATA_DIR" && find "$RESTORE" -type f 2>/dev/null ) | sort -u > "$TMPD/ondisk"
    n_missing="$(comm -23 "$TMPD/members" "$TMPD/ondisk" | wc -l | tr -cd '0-9')"
    bytes_disk="$( ( cd "$DATA_DIR" && tr '\n' '\0' < "$TMPD/members" | stat_sizes ) | awk '{s+=$1} END{print s+0}' )"
    if [[ "${n_missing:-0}" -gt 0 || "$bytes_disk" != "$bytes_want" ]]; then
        # A partial restore DID land VTK, so leave the marker: the run will
        # (correctly) read ARCHIVED-STALE on the next pass, annotated with why.
        # ARCHIVED.txt is NOT cleared -- the tarball is still authoritative.
        printf 'restored %s INCOMPLETE\nfrom %s\ndate %s\nmissing %s\n' \
            "$RESTORE" "$tarball" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$n_missing" > "$DATA_DIR/$RESTORE/RESTORED.txt"
        if [[ "${n_missing:-0}" -gt 0 ]]; then
            echo "RESTORE-INCOMPLETE $RESTORE  ($n_missing members absent on disk; first: $(comm -23 "$TMPD/members" "$TMPD/ondisk" | head -1))"
        else
            echo "RESTORE-INCOMPLETE $RESTORE  (members present but ${bytes_disk} bytes on disk != ${bytes_want} archived -- a file is truncated, likely from an interrupted earlier restore; -k will never overwrite it)"
        fi
        echo "         tarball intact at $tarball -- fix the cause, remove the bad files, and re-run --restore"
        exit 8
    fi

    rm -f "$DATA_DIR/$RESTORE/ARCHIVED.txt"
    # A restore leaves tarball + VTK-on-/home, which is byte-for-byte the
    # interrupted-delete state.  Leave a marker so the next archiver pass can
    # say WHY it is stale.  The run is still reported ARCHIVED-STALE and still
    # goes to Ryan -- the marker annotates the alert, it never suppresses it.
    printf 'restored %s\nfrom %s\ndate %s\n' \
        "$RESTORE" "$tarball" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$DATA_DIR/$RESTORE/RESTORED.txt"
    echo "RESTORED $RESTORE  files=$(find "$DATA_DIR/$RESTORE" -type f | wc -l | tr -cd '0-9')"
    echo "ARCHIVE_MODE=APPLY"
    exit 0
fi

# -------------------------------------------------------------- resume-delete
# For an ARCHIVED-STALE run: the tarball landed but the /home ParaView delete
# did not complete.  Ryan authorises this per run; it is never batched.  It
# honours the same $KEEP_STEPS restart window as a fresh archive, so a resumed
# delete leaves the run warm-startable in place too.  Safety: every file it is
# about to remove must be present in the tarball at the same size first.
if [[ -n "$RESUME_DELETE" ]]; then
    run="$RESUME_DELETE"
    root="$(single_root --resume-delete)"
    DATA_DIR="$root/data"
    SLUG="$(slug_for "$root")"
    acquire_checkout_lock "$SLUG" || { checkout_locked_msg "$SLUG"; exit 6; }
    dir="$DATA_DIR/$run"
    [[ -d "$dir" ]] || { echo "FATAL: --resume-delete '$run' matches no directory under $DATA_DIR/" >&2; exit 5; }
    tarball="$(find_tarball "$ARCHIVE_DIR/$SLUG" "$run")" || {
        echo "FATAL: no tarball for '$run' under $ARCHIVE_DIR/$SLUG (checkout $(cd "$root" && pwd -P))" >&2; exit 5; }
    case "$tarball" in
        *.tar.zst) DECOMP_CMD=(zstd -dcq) ;;
        *.tar.gz)  DECOMP_CMD=(gzip -dc) ;;
        *)         DECOMP_CMD=(cat) ;;
    esac

    # Paths inside the tarball are "<run>/...", which is exactly what pv_find
    # yields from inside $DATA_DIR, so the two sets are directly comparable.
    remaining="$TMPD/remaining"; rem_paths="$TMPD/rem_paths"; intar="$TMPD/intar"; missing_f="$TMPD/missing"
    sel="$TMPD/sel"
    ( cd "$DATA_DIR" && select_deletable "$run" "$run" "$KEEP_STEPS" ) > "$sel" 2>/dev/null || true
    keptsteps="$(awk -F'\t' '/^KEEP/{print $3}' "$sel")"
    awk -F'\t' '/^DEL/{printf "%s\t%s\n", $2, $3}' "$sel" | sort -t"$(printf '\t')" -k2 > "$remaining"
    cut -f2- "$remaining" | sort > "$rem_paths"
    if [[ ! -s "$remaining" ]]; then
        echo "RESUME-DELETE $run  nothing outside the ${KEEP_STEPS}-step window; already complete"
        echo "ARCHIVE_MODE=$( $APPLY && echo APPLY || echo DRY-RUN )"
        exit 0
    fi
    if ! "${DECOMP_CMD[@]}" < "$tarball" | tar -tf - | sed 's|^\./||' | sort -u > "$intar"; then
        echo "VERIFY-FAIL $run  (tarball unreadable)"; exit 8
    fi

    nrem="$(wc -l < "$rem_paths" | tr -cd '0-9')"
    comm -23 "$rem_paths" "$intar" > "$missing_f"
    missing="$(wc -l < "$missing_f" | tr -cd '0-9')"
    if (( missing > 0 )); then
        echo "VERIFY-FAIL $run  ($missing of $nrem files due for deletion are absent from $tarball)"
        echo "         first missing: $(head -1 "$missing_f")"
        exit 8
    fi

    # Membership is not enough -- the archived copy must also carry the same
    # bytes.  One pass, restricted to exactly the members about to be deleted.
    want_bytes="$(awk -F'\t' '{s+=$1} END{print s+0}' "$remaining")"
    got_bytes="$("${DECOMP_CMD[@]}" < "$tarball" | tar -xOf - -T "$rem_paths" 2>/dev/null | wc -c | tr -cd '0-9')"
    if [[ "$got_bytes" != "$want_bytes" ]]; then
        echo "VERIFY-FAIL $run  (archived bytes $got_bytes != on-disk $want_bytes for the $nrem files due for deletion)"
        exit 8
    fi

    mb=$(( want_bytes / 1048576 ))
    echo "RESUME-DELETE $run  verified=${nrem}files/${mb}MB byte-for-byte in $tarball  keeping=[$keptsteps]"
    if $APPLY; then
        ( cd "$DATA_DIR" && tr '\n' '\0' < "$rem_paths" | xargs -0 rm -f )
        printf 'archived %s\ncheckout_slug %s\nsource_checkout %s\ntarball %s\nfiles %s\ndate %s\nkept_steps %s\nnote resumed delete; kept steps are warm-startable in place\nrestore: scripts/run_archiver.sh --root %s --restore %s --apply\n' \
            "$run" "$SLUG" "$(cd "$root" && pwd -P)" "$tarball" "$nrem" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$keptsteps" "$(cd "$root" && pwd -P)" "$run" > "$dir/ARCHIVED.txt"
        rm -f "$dir/RESTORED.txt"
        echo "APPLIED $run  freed=${mb}MB"
    fi
    echo "ARCHIVE_MODE=$( $APPLY && echo APPLY || echo DRY-RUN )"
    exit 0
fi

# ------------------------------------------------------------------- archive
# --only takes one run or a comma-separated list, so a batch of RECENT runs can
# be approved in one go without ever becoming a blanket "archive everything".
in_only() {
    [[ -z "$ONLY" ]] && return 0
    local x
    IFS=',' read -r -a _only <<< "$ONLY"
    for x in "${_only[@]}"; do [[ "$1" == "$x" ]] && return 0; done
    return 1
}

if [[ -n "$ONLY" ]]; then
    IFS=',' read -r -a _only_check <<< "$ONLY"
    for o in "${_only_check[@]}"; do
        # A name must resolve to EXACTLY ONE checkout.  Run names are not
        # unique across clones: approving `fm052d_gpu_1080` once also selected
        # an identically-named, unapproved run in a different checkout
        # (observed 2026-08-27; only the dry run caught it).  Refusing here
        # turns that class of mistake into a hard error instead of a
        # vigilance test -- scope with a single --root to disambiguate.
        matches=()
        for r in "${ROOTS[@]}"; do [[ -d "$r/data/$o" ]] && matches+=("$r"); done
        if [[ ${#matches[@]} -eq 0 ]]; then
            echo "FATAL: --only '$o' matches no directory under any --root's data/" >&2
            exit 5
        fi
        if [[ ${#matches[@]} -gt 1 ]]; then
            echo "FATAL: --only '$o' matches runs in ${#matches[@]} checkouts: ${matches[*]}" >&2
            echo "       Run names are not unique across clones.  Scope with a single --root to name which one is approved." >&2
            exit 2
        fi
    done
fi

# RECENT runs are archived ONLY when Ryan has named them.  Without this guard
# --include-recent would silently become "archive everything that is not in the
# queue", which is precisely the judgement call the RECENT label exists to
# escalate.
if $INCLUDE_RECENT && [[ -z "$ONLY" ]]; then
    echo "FATAL: --include-recent must be scoped with --only RUN[,RUN...]." >&2
    echo "       It exists to act on runs Ryan has approved by name, not to archive every recent run." >&2
    exit 2
fi

if $APPLY; then
    mkdir -p "$ARCHIVE_DIR" 2>/dev/null || true
    if [[ ! -d "$ARCHIVE_DIR" || ! -w "$ARCHIVE_DIR" ]]; then
        echo "FATAL: $ARCHIVE_DIR is not a writable directory." >&2
        echo "       Archive is mounted on LOGIN NODES ONLY -- a compute node cannot see it." >&2
        exit 7
    fi
fi

echo "# run_archiver  keep=$KEEP_STEPS steps on /home  quiet_hours=$QUIET_HOURS  compress=$COMPRESS  apply=$APPLY  fast_verify=$FAST_VERIFY"
(( KEEP_STEPS == 0 )) && echo "# WARNING: --keep 0 -- archived runs will NOT be warm-startable without a --restore first"
echo "# archive dir:  $ARCHIVE_DIR   (one <checkout-slug>/ subdir per source checkout)"
echo "# checkouts:    ${#ROOTS[@]}"
for r in "${ROOTS[@]}"; do echo "#   $(slug_for "$r")  <-  $(cd "$r" && pwd -P)"; done
echo "# protect list: $PROTECT_FILE"
$squeue_ok || echo "# NOTE: squeue unavailable -- liveness rests on mtime alone"

total_src=0
total_tar=0
total_freed=0
stale_count=0
recent_count=0
verify_fail=0
locked_count=0

for root in "${ROOTS[@]}"; do
DATA_DIR="$root/data"
SLUG="$(slug_for "$root")"
if ! acquire_checkout_lock "$SLUG"; then
    checkout_locked_msg "$SLUG"
    locked_count=$(( locked_count + 1 ))
    continue
fi
ARCH_SUB="$ARCHIVE_DIR/$SLUG"
SRC_ABS="$(cd "$root" && pwd -P)"
GIT_REV="$(git -C "$root" rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "## checkout $SLUG  ($SRC_ABS, git $GIT_REV)"
$APPLY && mkdir -p "$ARCH_SUB" 2>/dev/null || true

for dir in "$DATA_DIR"/*/; do
    [[ -d "$dir" ]] || continue
    run="$(basename "$dir")"
    in_only "$run" || continue

    if is_protected "$run"; then
        echo "PROTECTED $run"
        continue
    fi

    existing=""
    for e in tar.zst tar.gz tar; do
        [[ -f "$ARCH_SUB/$run.$e" ]] && { existing="$ARCH_SUB/$run.$e"; break; }
    done

    if [[ -n "$existing" ]]; then
        # A tarball already here means the archive step ran.  If ParaView files
        # are ALSO still on /home the delete step did not complete -- report it
        # and stop; do not re-tar and do not silently delete.  Ryan decides,
        # then authorises --resume-delete for that run.
        # An archived run is EXPECTED to still hold its newest $KEEP_STEPS steps
        # on /home -- that is the point of the window.  Stale means VTK survives
        # OUTSIDE that window, i.e. the delete did not finish.
        sel="$TMPD/sel"
        select_deletable "$dir" "$run" "$KEEP_STEPS" > "$sel" 2>/dev/null || true
        del_n="$(grep -c '^DEL' "$sel" || true)"
        if [[ "${del_n:-0}" -gt 0 ]]; then
            del_mb=$(( $(awk -F'\t' '/^DEL/{s+=$2} END{print s+0}' "$sel") / 1048576 ))
            why="delete may have been interrupted"
            [[ -f "$dir/RESTORED.txt" ]] && why="marker says deliberately restored $(awk '/^date/{print $2}' "$dir/RESTORED.txt")"
            echo "ARCHIVED-STALE $run  tarball=$existing  outside_keep_window=${del_n}files/${del_mb}MB  -- ASK RYAN, $why"
            stale_count=$(( stale_count + 1 ))
        else
            echo "ARCHIVED-ALREADY $run  tarball=$existing  keeping=[$(awk -F'\t' '/^KEEP/{print $3}' "$sel")]"
        fi
        continue
    fi

    jm="$(queue_match "$run" || true)"
    if [[ -n "$jm" ]]; then
        echo "LIVE $run  (queue job '$jm' -- archiver never touches a live run)"
        continue
    fi

    # Checked before the RECENT branch: a run with no ParaView files has nothing
    # to reclaim, so putting it in Ryan's approval queue would be pure noise.
    pv_n="$(pv_find "$dir" | wc -l | tr -cd '0-9')"
    if [[ "${pv_n:-0}" -eq 0 ]]; then
        echo "NO-VTK $run  (nothing to reclaim)"
        continue
    fi

    qh="$(quiet_hours_of "$dir")"
    if (( qh < QUIET_HOURS )); then
        # A run written in the last couple of hours is almost certainly still
        # being written -- and the queue name match cannot be relied on to say
        # so (it failed to match `p018_csarc_n2_nt72_l3p0` while its job was
        # actively writing 55 GB into it, 2026-08-27).  Refuse approval outright
        # rather than let a hand-approved mistake reach the tar step.
        if (( qh < RECENT_HOT_HOURS )); then
            pv_mb=$(( $(pv_find "$dir" -print0 2>/dev/null | stat_sizes | awk '{s+=$1} END{print s+0}') / 1048576 ))
            echo "RECENT-HOT $run  quiet=${qh}h (<${RECENT_HOT_HOURS}h), vtk=${pv_mb}MB  -- too recently written to archive; NOT approvable, wait for it to go quiet"
            recent_count=$(( recent_count + 1 ))
            continue
        fi
        if $INCLUDE_RECENT && in_only "$run"; then
            echo "RECENT-APPROVED $run  quiet=${qh}h (<${QUIET_HOURS}h) -- archiving on Ryan's explicit --only"
        else
            pv_mb=$(( $(pv_find "$dir" -print0 2>/dev/null | stat_sizes | awk '{s+=$1} END{print s+0}') / 1048576 ))
            echo "RECENT $run  quiet=${qh}h (<${QUIET_HOURS}h), no queue job, vtk=${pv_mb}MB  -- ASK RYAN: looks finished but was touched recently"
            recent_count=$(( recent_count + 1 ))
            continue
        fi
    fi

    # Provenance goes INSIDE the tarball so a tarball that is moved, renamed, or
    # found years later still names the clone it came from.  It is staged in
    # $TMPD and tarred with a second -C, NOT written into the run directory:
    # touching /home would bump the run's mtime and make it look LIVE on the
    # next pass.
    manifest="$TMPD/manifest"
    build_manifest "$run" "$manifest"
    prov_dir="$TMPD/prov/$run"
    prov_rel="$run/ARCHIVE_SOURCE.txt"
    if $APPLY; then
        mkdir -p "$prov_dir"
        printf 'run %s\nsource_checkout %s\nsource_data_dir %s\ncheckout_slug %s\ngit_rev %s\nhost %s\narchived_by run_archiver.sh\narchived_utc %s\n' \
            "$run" "$SRC_ABS" "$SRC_ABS/data/$run" "$SLUG" "$GIT_REV" \
            "$(hostname -s 2>/dev/null || echo unknown)" \
            "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$prov_dir/ARCHIVE_SOURCE.txt"
        # it is a tar member like any other, so verification must expect it
        printf '%s\t%s\n' \
            "$("${STAT_FMT[@]}" "$prov_dir/ARCHIVE_SOURCE.txt" | awk '{print $1}')" \
            "$prov_rel" >> "$manifest"
    fi
    nfiles="$(wc -l < "$manifest" | tr -cd '0-9')"
    src_bytes="$(awk -F'\t' '{s+=$1} END{print s+0}' "$manifest")"
    src_mb=$(( src_bytes / 1048576 ))

    sel="$TMPD/sel"
    select_deletable "$dir" "$run" "$KEEP_STEPS" > "$sel" 2>/dev/null || true
    nrestart="$(awk -F'\t' '/^KEEP/{print $2}' "$sel")"
    keptsteps="$(awk -F'\t' '/^KEEP/{print $3}' "$sel")"
    del_n="$(grep -c '^DEL' "$sel" || true)"
    del_bytes="$(awk -F'\t' '/^DEL/{s+=$2} END{print s+0}' "$sel")"
    del_mb=$(( del_bytes / 1048576 ))
    [[ "${nrestart:-0}" -eq 0 ]] && keptsteps="none: no restartable step"

    if ! $APPLY; then
        echo "ARCHIVE $run  files=$nfiles  src=${src_mb}MB  restartable=${nrestart:-0}  keeping=[$keptsteps]  delete=${del_n}files/${del_mb}MB (would)"
        total_src=$(( total_src + src_bytes ))
        total_freed=$(( total_freed + del_bytes ))
        continue
    fi

    mtime_before="$(newest_mtime "$dir")"
    if [[ -z "$mtime_before" ]]; then
        # newest_mtime is ||-true-guarded (load-bearing, see its comment), so
        # a transient find/stat failure surfaces as an EMPTY string.  Comparing
        # "" against a later "" would silently PASS the changed-during-archive
        # guard; an unreadable snapshot is a hard no-go for this run instead.
        echo "MTIME-UNREADABLE $run  (cannot snapshot pre-archive mtime; nothing archived, nothing deleted)"
        verify_fail=$(( verify_fail + 1 ))
        continue
    fi
    part="$ARCH_SUB/$run.$EXT.partial"
    final="$ARCH_SUB/$run.$EXT"
    rm -f "$part"

    # Stream straight to archive -- never stage a temp copy on /home, which is
    # the filesystem under pressure.  The .partial name means an interrupted
    # transfer can never be mistaken for a complete archive.
    if ! tar -cf - -C "$DATA_DIR" "$run" -C "$TMPD/prov" "$prov_rel" | "${COMP_CMD[@]}" > "$part"; then
        echo "VERIFY-FAIL $run  (tar/compress failed; leaving $part for inspection)"
        verify_fail=$(( verify_fail + 1 ))
        continue
    fi
    mv "$part" "$final"

    # A run that is being written while we tar it would produce an inconsistent
    # archive.  The mtime snapshot taken before the tar is the check: if anything
    # moved, we do NOT delete, regardless of what the classification said.  This
    # is the backstop for a queue-match false negative and for any RECENT run
    # approved by hand that turned out to still be active.
    # Quarantine names carry host+pid+start-time, not bare $$: PIDs collide
    # across login nodes, and one node overwriting another's preserved
    # evidence is how a bad tarball's forensics vanish.
    mtime_after="$(newest_mtime "$dir")"
    if [[ -z "$mtime_after" || "$mtime_after" != "$mtime_before" ]]; then
        mv "$final" "$final.changed.$LOCK_TOKEN" 2>/dev/null || true
        verify_fail=$(( verify_fail + 1 ))
        echo "CHANGED-DURING-ARCHIVE $run  the run was written while being tarred (or its mtime became unreadable); nothing deleted, tarball moved to $final.changed.$LOCK_TOKEN"
        continue
    fi

    if ! verify_tarball "$final" "$manifest" "$run"; then
        # Rename it out of the way: a file sitting at the canonical archive name
        # is indistinguishable from a good archive on the next pass, and would be
        # reported ARCHIVED-ALREADY.  Renamed, it is preserved for inspection and
        # the next run re-archives from scratch.
        mv "$final" "$final.corrupt.$LOCK_TOKEN" 2>/dev/null || true
        verify_fail=$(( verify_fail + 1 ))
        echo "         nothing deleted for $run; bad tarball moved to $final.corrupt.$LOCK_TOKEN"
        continue
    fi

    tar_bytes="$("${STAT_FMT[@]}" "$final" | awk '{print $1}')"
    tar_mb=$(( tar_bytes / 1048576 ))

    # The verify read-back can take many minutes on a large run, so the world
    # may have changed since the pre-tar snapshot.  Every no-go signal is
    # re-checked IMMEDIATELY before the only rm in this script:
    #   * mtime: a write landing during verify would otherwise be deleted with
    #     its new content unarchived (the tarball holds the OLD bytes);
    #   * protect list: read FRESH, so "quick, protect that run!" still works
    #     mid-pass on a long-lived worker;
    #   * squeue: a job that started mid-pass makes the run live again.
    # In every no-go case the verified tarball stays (or is quarantined) and
    # NOTHING is deleted; the run surfaces as ARCHIVED-STALE on the next pass,
    # which is the existing ask-Ryan path.
    mtime_final="$(newest_mtime "$dir")"
    if [[ -z "$mtime_final" || "$mtime_final" != "$mtime_before" ]]; then
        mv "$final" "$final.changed.$LOCK_TOKEN" 2>/dev/null || true
        verify_fail=$(( verify_fail + 1 ))
        echo "CHANGED-BEFORE-DELETE $run  written between tar and delete; nothing deleted, tarball moved to $final.changed.$LOCK_TOKEN"
        continue
    fi
    if [[ -f "$PROTECT_FILE" ]]; then
        protected_runs="$(sed -e 's/#.*//' -e 's/[[:space:]]//g' "$PROTECT_FILE" | grep -v '^$' || true)"
    fi
    if is_protected "$run"; then
        echo "PROTECTED-LATE $run  protected mid-pass; tarball kept at $final, NOTHING deleted -- will read ARCHIVED-STALE next pass, ask Ryan"
        stale_count=$(( stale_count + 1 ))
        continue
    fi
    if $squeue_ok; then
        fresh_names="$(squeue -u "$USER" -h -o '%j' 2>/dev/null | tr '-' '_' | sed 's/^fp_//')" && job_names="$fresh_names"
    fi
    jm="$(queue_match "$run" || true)"
    if [[ -n "$jm" ]]; then
        echo "LIVE-LATE $run  (queue job '$jm' appeared mid-pass; tarball kept at $final, NOTHING deleted -- will read ARCHIVED-STALE next pass)"
        stale_count=$(( stale_count + 1 ))
        continue
    fi

    # Verified.  Now, and only now, remove the VTK outside the restart window.
    if [[ "${del_n:-0}" -gt 0 ]]; then
        awk -F'\t' '/^DEL/{print $3}' "$sel" | tr '\n' '\0' | xargs -0 rm -f
    fi
    printf 'archived %s\ncheckout_slug %s\nsource_checkout %s\ntarball %s\nfiles %s\nsrc_bytes %s\ntar_bytes %s\ndate %s\nkept_steps %s\nnote the kept steps are warm-startable in place; no unpacking needed\nrestore: scripts/run_archiver.sh --root %s --restore %s --apply\n' \
        "$run" "$SLUG" "$SRC_ABS" "$final" "$nfiles" "$src_bytes" "$tar_bytes" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$keptsteps" "$SRC_ABS" "$run" > "$dir/ARCHIVED.txt"

    # Append-only index at the archive root: one row per archived run, so
    # "which clone did this come from?" is a single grep over everything.
    # HEADER VIA >> ONLY: with concurrent workers, `>` after a stale [[ -f ]]
    # check truncates rows another worker appended in between.  Appending can
    # at worst duplicate the header line; it can never lose a row.
    idx="$ARCHIVE_DIR/INDEX.tsv"
    [[ -f "$idx" ]] || printf 'archived_utc\tcheckout_slug\trun\ttarball\tsrc_bytes\ttar_bytes\tfiles\tkept_steps\tgit_rev\tsource_checkout\n' >> "$idx" 2>/dev/null || true
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$SLUG" "$run" "${final#"$ARCHIVE_DIR"/}" \
        "$src_bytes" "$tar_bytes" "$nfiles" "${keptsteps:-none}" "$GIT_REV" "$SRC_ABS" >> "$idx" 2>/dev/null || true

    echo "ARCHIVE $run  files=$nfiles  src=${src_mb}MB  tar=${tar_mb}MB  keeping=[$keptsteps]  freed=${del_mb}MB  -> $SLUG/$(basename "$final")"
    total_src=$(( total_src + src_bytes ))
    total_tar=$(( total_tar + tar_bytes ))
    total_freed=$(( total_freed + del_bytes ))
done
rm -rf "$CKLOCK"    # free this checkout for other workers as soon as we leave it
done

echo "TOTAL_ARCHIVED_MB=$(( total_tar / 1048576 ))"
echo "TOTAL_SOURCE_MB=$(( total_src / 1048576 ))"
echo "TOTAL_FREED_MB=$(( total_freed / 1048576 ))"
echo "STALE_COUNT=$stale_count"
echo "RECENT_COUNT=$recent_count"
echo "VERIFY_FAIL_COUNT=$verify_fail"
echo "CHECKOUT_LOCKED_COUNT=$locked_count"
echo "DF_AFTER=$(df -BG "$HOME" 2>/dev/null | tail -1 | awk '{print $3}')"
echo "ARCHIVE_MODE=$( $APPLY && echo APPLY || echo DRY-RUN )"

# Exit status carries the conditions a caller must not miss.  Note this
# breaks the sweeper's "non-zero means nothing happened" convention: runs that
# archived cleanly in this pass DID archive.  8 means the runs that failed
# verification deleted nothing; 9 means a human decision is owed; 6 means a
# checkout was skipped because another archiver holds its lock.  A launcher
# running several workers MUST collect and OR every worker's exit status --
# losing an 8 or a 9 in one of N scrollbacks means "do not trust this pass"
# went unread.
if (( verify_fail > 0 )); then exit 8; fi
if (( stale_count > 0 )); then exit 9; fi
if (( locked_count > 0 )); then exit 6; fi
