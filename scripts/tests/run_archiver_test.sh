#!/usr/bin/env bash
# Behaviour suite for scripts/run_archiver.sh against a synthetic data/ tree.
#
# Every check here exists because the corresponding failure is unrecoverable in
# production: this script is the only thing standing between a verified archive
# and an irreversible delete of ParaView state.  Run it after ANY edit to
# run_archiver.sh, and before mirroring the script to the cluster.
#
#   bash scripts/tests/run_archiver_test.sh
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
S="$REPO/scripts/run_archiver.sh"
SW="$REPO/scripts/p018_vtk_sweeper.sh"
WORK="${WORK:-${TMPDIR:-/tmp}/run_archiver_test.$$}"
mkdir -p "$WORK" && cp "$REPO/scripts/tests/run_archiver_fixture.sh" "$WORK/mkfx.sh"
cd "$WORK"
trap 'rm -rf "$WORK"' EXIT
ZSTD="$(command -v zstd)"
FAILED=0
pass(){ echo "  PASS $1"; }; fail(){ echo "  ** FAIL $1"; FAILED=1; }
nvtk(){ find "data/$1" -name '*.vt[upsm]' -not -path '*/monitors/*' | wc -l | tr -d ' '; }
# steps still present on /home, newest first
steps(){ ls "data/$1/$1_wake1_particles" 2>/dev/null | sed 's/.*\.\([0-9]*\)\.vtp/\1/' | sort -rn | tr '\n' ',' ; }
# LOCKROOT lives inside the fixture so tests never touch ~/.cache and a
# crashed test run cannot leave a stale lock that blocks the next one.
fresh(){ bash mkfx.sh; cd fx || exit 1; mkdir -p locks; export ARCHIVE_DIR="$PWD/arch" PROTECT_FILE="$PWD/BR/protect.txt" LOCKROOT="$PWD/locks"; }

fresh
echo "T1 dry-run classification"
out=$(bash "$S" 2>&1)
for e in "ARCHIVE finished_a" "ARCHIVE finished_b" "RECENT-HOT live_d" "RECENT recent_f" "NO-VTK novtk_e" "PROTECTED protected_c" "ARCHIVE_MODE=DRY-RUN"; do
  echo "$out" | grep -q "$e" && pass "$e" || fail "$e"; done
[ -z "$(ls -A arch)" ] && pass "dry-run wrote nothing" || fail "dry-run wrote nothing"
diff -rq ref data >/dev/null && pass "dry-run deleted nothing" || fail "dry-run deleted nothing"

echo "T2 apply: vtk migrates, RESTART WINDOW stays, residue stays"
bash "$S" --apply >/dev/null 2>&1; rc=$?
[ "$rc" = 0 ] && pass "exit 0" || fail "exit 0 (got $rc)"
# finished_a has 8 steps; the newest 5 (4..8) must survive in full, 1..3 must go
[ "$(nvtk finished_a)" = 25 ] && pass "5 steps x 5 files kept" || fail "5 steps x 5 files kept (got $(nvtk finished_a))"
[ "$(steps finished_a)" = "8,7,6,5,4," ] && pass "kept the NEWEST 5 steps" || fail "kept newest 5 (got $(steps finished_a))"
for st in 4 5 6 7 8; do
  n=$(find "data/finished_a" -name "*.$st.vt*" -not -path '*/monitors/*' | wc -l | tr -d ' ')
  [ "$n" = 5 ] && pass "step $st complete (5 files incl. .vtm)" || fail "step $st complete (got $n)"
done
[ "$(find data/finished_a -name '*.3.vt*' -not -path '*/monitors/*' | wc -l | tr -d ' ')" = 0 ] && pass "step 3 removed" || fail "step 3 removed"
# finished_b has only 3 steps -- fewer than the window, so nothing is deleted
[ "$(nvtk finished_b)" = 15 ] && pass "short run keeps everything" || fail "short run keeps everything"
for f in finished_a_CT_vs_rev.csv finished_a.metadata.toml finished_a_case_metadata.toml \
         finished_a_body1.pvd monitors/finished_a_force.csv monitors/finished_a_never_touch.vtu; do
  [ -f "data/finished_a/$f" ] && pass "kept $f" || fail "kept $f"; done
diff -rq ref/live_d      data/live_d      >/dev/null && pass "recent untouched"    || fail "recent untouched"
diff -rq ref/protected_c data/protected_c >/dev/null && pass "protected untouched" || fail "protected untouched"

echo "T3 restore is byte-identical"
rm -rf data/finished_a; bash "$S" --restore finished_a --apply >/dev/null 2>&1
rm -f data/finished_a/RESTORED.txt data/finished_a/ARCHIVE_SOURCE.txt
diff -r ref/finished_a data/finished_a >/dev/null && pass "round-trip identical" || fail "round-trip identical"

echo "T4 ARCHIVED-STALE: vtk present OUTSIDE the keep window"
out=$(bash "$S" 2>&1); rc=$?
echo "$out" | grep -q 'ARCHIVED-STALE finished_a'   && pass "stale flagged"    || fail "stale flagged"
echo "$out" | grep -q 'ARCHIVED-ALREADY finished_b' && pass "clean => ALREADY"  || fail "clean => ALREADY"
echo "$out" | grep -q 'STALE_COUNT=1'               && pass "STALE_COUNT=1"     || fail "STALE_COUNT=1"
echo "$out" | grep -q 'ASK RYAN'                    && pass "escalation text"   || fail "escalation text"
[ "$rc" = 9 ] && pass "exit 9" || fail "exit 9 (got $rc)"
[ "$(nvtk finished_a)" = 40 ] && pass "stale run untouched" || fail "stale run untouched"

echo "T5 resume-delete refuses when bytes disagree"
echo "TAMPERED-and-longer-than-the-original" > data/finished_a/finished_a_body1/finished_a_body1.2.vtu
bash "$S" --resume-delete finished_a --apply >/dev/null 2>&1; rc=$?
[ "$rc" = 8 ] && pass "exit 8" || fail "exit 8 (got $rc)"
[ "$(nvtk finished_a)" = 40 ] && pass "deleted nothing" || fail "deleted nothing"

echo "T6 resume-delete completes an interrupted delete"
cd ..; fresh
bash "$S" --apply >/dev/null 2>&1
rm -rf data/finished_a; bash "$S" --restore finished_a --apply >/dev/null 2>&1
out=$(bash "$S" --resume-delete finished_a --apply 2>&1)
echo "$out" | grep -q 'verified=15files' && pass "verified 15 files (outside window)" || fail "verified 15 files"
[ "$(nvtk finished_a)" = 25 ]            && pass "window restored to 5 steps" || fail "window restored (got $(nvtk finished_a))"
echo "$out" | grep -q 'keeping=\[4,5,6,7,8\]' && pass "resume keeps newest 5 (label ascending, as the sweeper does)" || fail "resume keeps newest 5"
[ -f data/finished_a/finished_a_CT_vs_rev.csv ] && pass "csv kept" || fail "csv kept"
[ ! -f data/finished_a/RESTORED.txt ]    && pass "restore marker cleared" || fail "restore marker cleared"

# T7 also guards a real regression: writing a breadcrumb into the run directory
# once bumped its mtime, so a run whose archive failed looked LIVE and could not
# be retried for 24 h.  If "recovery exit 0" starts reporting LIVE, that is back.
echo "T7 corrupt tarball: quarantined, nothing deleted, next pass recovers"
cd ..; fresh
mkdir -p shim
printf '#!/usr/bin/env bash\nfor a in "$@"; do [ "$a" = "-dcq" ] && exec %s "$@"; done\n%s "$@" | head -c 200\n' "$ZSTD" "$ZSTD" > shim/zstd
chmod +x shim/zstd
PATH="$PWD/shim:$PATH" bash "$S" --only finished_a --apply >/dev/null 2>&1; rc=$?
[ "$rc" = 8 ] && pass "exit 8" || fail "exit 8 (got $rc)"
find arch -name '*corrupt*' | grep -q . && pass "quarantined"        || fail "quarantined"
find arch -name 'finished_a.tar.zst' | grep -q . && fail "no valid name left" || pass "no valid name left"
[ "$(nvtk finished_a)" = 40 ]               && pass "deleted nothing"        || fail "deleted nothing"
bash "$S" --only finished_a --apply >/dev/null 2>&1; rc=$?
[ "$rc" = 0 ] && pass "recovery exit 0" || fail "recovery exit 0 (got $rc)"
[ "$(nvtk finished_a)" = 25 ] && pass "recovered to window" || fail "recovered to window"

echo "T8 --keep window is honoured"
cd ..; fresh
bash "$S" --only finished_a --keep 2 --apply >/dev/null 2>&1
[ "$(steps finished_a)" = "8,7," ] && pass "--keep 2 keeps steps 8,7" || fail "--keep 2 (got $(steps finished_a))"
cd ..; fresh
out=$(bash "$S" --only finished_a --keep 0 --apply 2>&1)
[ "$(nvtk finished_a)" = 0 ] && pass "--keep 0 strips all vtk" || fail "--keep 0 strips all vtk"
echo "$out" | grep -q 'WARNING: --keep 0' && pass "--keep 0 warns" || fail "--keep 0 warns"
[ -f data/finished_a/finished_a_CT_vs_rev.csv ] && pass "--keep 0 still keeps csv" || fail "--keep 0 still keeps csv"

echo "T9 RECENT: looks finished but touched recently -> ask Ryan"
cd ..; fresh
out=$(bash "$S" 2>&1)
echo "$out" | grep -q 'RECENT recent_f' && pass "recent flagged" || fail "recent flagged"
echo "$out" | grep -q 'ASK RYAN'        && pass "escalation text" || fail "escalation text"
echo "$out" | grep -q 'RECENT_COUNT=3'  && pass "RECENT_COUNT=3"  || fail "RECENT_COUNT (got $(echo "$out" | grep RECENT_COUNT))"
echo "$out" | grep -q 'RECENT-HOT live_d' && pass "sub-2h run is RECENT-HOT" || fail "sub-2h run is RECENT-HOT"
echo "$out" | grep -q 'NOT approvable'     && pass "hot run marked not approvable" || fail "hot run marked not approvable"
bash "$S" --apply >/dev/null 2>&1
[ "$(nvtk recent_f)" = 35 ] && pass "recent NOT archived by default" || fail "recent NOT archived by default"
[ -z "$(find arch -name 'recent_f.tar.zst' 2>/dev/null)" ] && pass "no tarball for recent" || fail "no tarball for recent"
# blanket approval is refused; per-run approval works
bash "$S" --include-recent --apply >/dev/null 2>&1
[ $? = 2 ] && pass "--include-recent without --only refused" || fail "--include-recent without --only refused"
out=$(bash "$S" --include-recent --only recent_f --apply 2>&1)
echo "$out" | grep -q 'RECENT-APPROVED recent_f' && pass "named recent run approved" || fail "named recent run approved"
[ "$(nvtk recent_f)" = 25 ] && pass "approved run trimmed to the 5-step window" || fail "approved run trimmed (got $(nvtk recent_f))"
[ -n "$(find arch -name 'recent_f.tar.zst')" ] && pass "tarball written for approved run" || fail "tarball written"
[ "$(nvtk live_d)" = 15 ] && pass "unapproved run untouched" || fail "unapproved run untouched"
# --only takes a list
cd ..; fresh
out=$(bash "$S" --include-recent --only recent_f,recent_g --apply 2>&1)
[ "$(echo "$out" | grep -c 'RECENT-APPROVED')" = 2 ] && pass "--only accepts a list" || fail "--only accepts a list"
# a sub-2h run cannot be approved even when named explicitly
cd ..; fresh
out=$(bash "$S" --include-recent --only live_d --apply 2>&1)
echo "$out" | grep -q 'RECENT-HOT live_d' && pass "hot run refused despite --only" || fail "hot run refused despite --only"
[ -z "$(find arch -name 'live_d.tar.zst' 2>/dev/null)" ] && pass "no tarball for hot run" || fail "no tarball for hot run"

echo "T10 a run written DURING the archive is never deleted"
cd ..; fresh
# Deterministic: a slow compressor shim holds the tar open long enough for the
# background write to land inside the window the mtime guard watches.
mkdir -p shim
printf '#!/usr/bin/env bash\nfor a in "$@"; do [ "$a" = "-dcq" ] && exec %s "$@"; done\nsleep 3\nexec %s "$@"\n' "$ZSTD" "$ZSTD" > shim/zstd
chmod +x shim/zstd
( sleep 1; echo "a late write from a still-running job" >> data/finished_a/finished_a_body1/finished_a_body1.1.vtu ) &
out=$(PATH="$PWD/shim:$PATH" bash "$S" --only finished_a --apply 2>&1); rc=$?
wait
echo "$out" | grep -q 'CHANGED-DURING-ARCHIVE finished_a' && pass "change detected" || fail "change detected"
[ "$rc" = 8 ] && pass "exit 8" || fail "exit 8 (got $rc)"
[ "$(nvtk finished_a)" = 40 ] && pass "nothing deleted" || fail "nothing deleted (got $(nvtk finished_a))"
find arch -name '*.changed.*' | grep -q . && pass "tarball quarantined" || fail "tarball quarantined"

echo "T11 multi-checkout: per-source organization and provenance"
cd ..; fresh
bash "$S" --root . --root ../fx2 --apply >/dev/null 2>&1; rc=$?
[ "$rc" = 0 ] && pass "exit 0 across 2 checkouts" || fail "exit 0 across 2 checkouts (got $rc)"
SL1=$(basename "$(ls -d arch/*_fx)"); SL2=$(basename "$(ls -d arch/*_fx2)")
[ -f "arch/$SL1/finished_a.tar.zst" ] && pass "checkout 1 filed under its own slug" || fail "checkout 1 slug dir"
[ -f "arch/$SL2/finished_a.tar.zst" ] && pass "checkout 2 filed under its own slug" || fail "checkout 2 slug dir"
[ -f "arch/$SL2/other_only.tar.zst" ] && pass "second checkout's own run archived" || fail "second checkout's own run"
# the same-named runs must NOT have clobbered each other
a=$(zstd -dcq "arch/$SL1/finished_a.tar.zst" | tar -tf - | wc -l | tr -d ' ')
b=$(zstd -dcq "arch/$SL2/finished_a.tar.zst" | tar -tf - | wc -l | tr -d ' ')
[ "$a" != "$b" ] && pass "same-named runs kept distinct ($a vs $b members)" || fail "same-named runs kept distinct"
# provenance inside the tarball names the right clone
src=$(zstd -dcq "arch/$SL2/finished_a.tar.zst" | tar -xOf - finished_a/ARCHIVE_SOURCE.txt | awk '/^source_checkout/{print $2}')
case "$src" in *?/fx2) pass "in-tarball provenance names the source clone";; *) fail "in-tarball provenance (got $src)";; esac
# index answers "which clone?" in one grep
[ -f arch/INDEX.tsv ] && pass "INDEX.tsv written" || fail "INDEX.tsv written"
[ "$(grep -c 'finished_a' arch/INDEX.tsv)" = 2 ] && pass "index has both finished_a rows" || fail "index rows"
head -1 arch/INDEX.tsv | grep -q 'checkout_slug' && pass "index has a header" || fail "index header"
# protect list spans checkouts
[ ! -f "arch/$SL1/protected_c.tar.zst" ] && pass "protect list honoured across roots" || fail "protect list across roots"
# --only must not guess which clone either: the same name in two checkouts once
# selected an unapproved 13.8 G run alongside the approved one (2026-08-27).
bash "$S" --root . --root ../fx2 --only finished_a >/dev/null 2>&1
[ $? = 2 ] && pass "--only refuses a name matching two checkouts" || fail "--only refuses ambiguous name"
out=$(bash "$S" --root . --root ../fx2 --only finished_a 2>&1); rc=$?
echo "$out" | grep -q 'Scope with a single --root' && pass "ambiguous --only names the remedy" || fail "ambiguous --only remedy text"
bash "$S" --root ../fx2 --only finished_a >/dev/null 2>&1
[ $? = 0 ] && pass "--only + single --root still accepted" || fail "--only + single --root accepted"
bash "$S" --root . --root ../fx2 --only other_only >/dev/null 2>&1
[ $? = 0 ] && pass "unique name across roots still accepted" || fail "unique name across roots accepted"
# restore must not guess which clone
bash "$S" --root . --root ../fx2 --restore finished_a --apply >/dev/null 2>&1
[ $? = 2 ] && pass "restore refuses ambiguous checkout" || fail "restore refuses ambiguous checkout"
rm -rf ../fx2/data/finished_a
bash "$S" --root ../fx2 --restore finished_a --apply >/dev/null 2>&1
[ "$(ls ../fx2/data/finished_a/finished_a_wake1_particles | wc -l | tr -d ' ')" = 6 ] && pass "restore pulled the CORRECT clone's copy" || fail "restore pulled correct clone"

echo "T12 a VTK-less run must not abort the pass"
cd ..; fresh
# novtk_e has no ParaView files at all.  GNU xargs runs `stat` even on empty
# input, and its exit 123 once killed an --all-checkouts pass partway through,
# silently skipping every later checkout.  The summary lines are the proof the
# pass reached the end.
out=$(bash "$S" --root . --root ../fx2 2>&1); rc=$?
[ "$rc" = 0 ] && pass "exit 0 with a VTK-less run present" || fail "exit 0 (got $rc)"
echo "$out" | grep -q 'NO-VTK novtk_e'   && pass "VTK-less run dismissed early" || fail "VTK-less run dismissed"
echo "$out" | grep -q 'RECENT novtk_e'   && fail "VTK-less run kept out of approval queue" || pass "VTK-less run kept out of approval queue"
echo "$out" | grep -q 'ARCHIVE_MODE='    && pass "pass reached its summary" || fail "pass reached its summary"
[ "$(echo "$out" | grep -c '^## checkout ')" = 2 ] && pass "both checkouts still reached" || fail "both checkouts reached"

echo "T13 --all-checkouts discovery"
cd ..; fresh
out=$(HOME="$(cd .. && pwd -P)" bash "$S" --all-checkouts 2>&1)
[ "$(echo "$out" | grep -c '^## checkout ')" = 2 ] && pass "found both checkouts" || fail "found both checkouts"
echo "$out" | grep -q 'ARCHIVE other_only' && pass "reached the second checkout's runs" || fail "reached second checkout"
mkdir -p ../notacheckout/data && touch ../notacheckout/data/x
out=$(HOME="$(cd .. && pwd -P)" bash "$S" --all-checkouts 2>&1)
[ "$(echo "$out" | grep -c '^## checkout ')" = 2 ] && pass "data/ without Project.toml is not a checkout" || fail "rejects non-checkout"
# a different Julia package (FastMultipole/FLOWVPM) must not be swept in
mkdir -p ../otherpkg/data
printf 'name = "FastMultipole"\nuuid = "ce07d0d3-2b9f-49ba-89eb-12c800257c85"\n' > ../otherpkg/Project.toml
# nor a checkout's own test/ environment ([deps]-only Project.toml)
mkdir -p ../fx2/test/data && printf '[deps]\nTest = "8dfed614-e22c-5e08-85e1-65c5234f0b40"\n' > ../fx2/test/Project.toml
out=$(HOME="$(cd .. && pwd -P)" bash "$S" --all-checkouts 2>&1)
[ "$(echo "$out" | grep -c '^## checkout ')" = 2 ] && pass "other packages and test/ envs excluded" || fail "other packages and test/ envs excluded (got $(echo "$out" | grep -c '^## checkout '))"
bash "$S" --all-checkouts --root . >/dev/null 2>&1; [ $? = 2 ] && pass "--all-checkouts + --root rejected" || fail "--all-checkouts + --root"

echo "T14 argument and environment guards"
cd ..; fresh
bash "$S" --only nope                      >/dev/null 2>&1; [ $? = 5 ] && pass "5 bad --only"      || fail "5 bad --only"
bash "$S" --quiet-hours 0                  >/dev/null 2>&1; [ $? = 2 ] && pass "2 bad hours"       || fail "2 bad hours"
bash "$S" --compress bz2                   >/dev/null 2>&1; [ $? = 2 ] && pass "2 bad compress"    || fail "2 bad compress"
bash "$S" --keep -1                        >/dev/null 2>&1; [ $? = 2 ] && pass "2 bad --keep"       || fail "2 bad --keep"
bash "$S" --keep two                       >/dev/null 2>&1; [ $? = 2 ] && pass "2 non-integer keep" || fail "2 non-integer keep"
bash "$S" --root /nope/nope                >/dev/null 2>&1; [ $? = 2 ] && pass "2 --root without data/" || fail "2 --root without data/"
bash "$S" --restore nope                   >/dev/null 2>&1; [ $? = 5 ] && pass "5 bad --restore"   || fail "5 bad --restore"
bash "$S" --resume-delete nope             >/dev/null 2>&1; [ $? = 5 ] && pass "5 bad --resume"    || fail "5 bad --resume"
PROTECT_FILE=/nope bash "$S"               >/dev/null 2>&1; [ $? = 3 ] && pass "3 no protect list" || fail "3 no protect list"
ARCHIVE_DIR=/proc/nope/x bash "$S" --apply >/dev/null 2>&1; [ $? = 7 ] && pass "7 archive unwritable (compute node)" || fail "7 archive unwritable"

echo "T15 restore verifies completeness instead of trusting tar"
# tar -xk exits nonzero on EXPECTED residue collisions, so its status is
# discarded -- which once meant a truncated tarball or ENOSPC still printed
# RESTORED.  The restore must now prove every member landed at full size.
cd ..; fresh
bash "$S" --only finished_a --apply >/dev/null 2>&1
rm -rf data/finished_a
tb=$(find arch -name 'finished_a.tar.zst')
cp "$tb" "$tb.good"
head -c "$(( $(wc -c < "$tb") / 2 ))" "$tb.good" > "$tb"
out=$(bash "$S" --restore finished_a --apply 2>&1); rc=$?
[ "$rc" = 8 ] && pass "truncated tarball -> exit 8" || fail "truncated tarball exit 8 (got $rc)"
echo "$out" | grep -q 'RESTORE-INCOMPLETE finished_a' && pass "RESTORE-INCOMPLETE reported" || fail "RESTORE-INCOMPLETE reported"
echo "$out" | grep -q '^RESTORED finished_a' && fail "no false RESTORED" || pass "no false RESTORED"
mv "$tb.good" "$tb"
rm -rf data/finished_a
bash "$S" --restore finished_a --apply >/dev/null 2>&1; rc=$?
[ "$rc" = 0 ] && pass "intact tarball restores, exit 0" || fail "intact restore exit 0 (got $rc)"
# a file half-written by an interrupted restore EXISTS, so -k skips it forever;
# existence alone would certify the re-restore -- the byte check must refuse
v=$(find data/finished_a -name '*.vtu' | head -1)
printf 'x' > "$v"
out=$(bash "$S" --restore finished_a --apply 2>&1); rc=$?
[ "$rc" = 8 ] && pass "truncated on-disk file -> exit 8" || fail "truncated on-disk file (got $rc)"
echo "$out" | grep -q 'RESTORE-INCOMPLETE finished_a' && pass "size mismatch reported" || fail "size mismatch reported"
# recovery: remove the bad file, restore fills the hole
rm -f "$v"
bash "$S" --restore finished_a --apply >/dev/null 2>&1; rc=$?
[ "$rc" = 0 ] && pass "restore after removing bad file, exit 0" || fail "restore after removing bad file (got $rc)"
rm -f data/finished_a/RESTORED.txt data/finished_a/ARCHIVE_SOURCE.txt
diff -r ref/finished_a data/finished_a >/dev/null && pass "recovery round-trip identical" || fail "recovery round-trip identical"

echo "T16 --fast-verify still catches a truncated tarball"
cd ..; fresh
mkdir -p shim
printf '#!/usr/bin/env bash\nfor a in "$@"; do [ "$a" = "-dcq" ] && exec %s "$@"; done\n%s "$@" | head -c 200\n' "$ZSTD" "$ZSTD" > shim/zstd
chmod +x shim/zstd
PATH="$PWD/shim:$PATH" bash "$S" --fast-verify --only finished_a --apply >/dev/null 2>&1; rc=$?
[ "$rc" = 8 ] && pass "exit 8" || fail "exit 8 (got $rc)"
find arch -name '*corrupt*' | grep -q . && pass "quarantined" || fail "quarantined"
[ "$(nvtk finished_a)" = 40 ] && pass "deleted nothing" || fail "deleted nothing"
bash "$S" --fast-verify --only finished_a --apply >/dev/null 2>&1; rc=$?
[ "$rc" = 0 ] && pass "fast-verify recovery exit 0" || fail "fast-verify recovery (got $rc)"
[ "$(nvtk finished_a)" = 25 ] && pass "recovered to window" || fail "recovered to window (got $(nvtk finished_a))"

# T17-T19 are reserved for the Phase 2 queue/claim machinery (see the plan);
# they are added together with that code, not before it.

echo "T20 late-change and unreadable-mtime guards before the strip (L5/L6)"
cd ..; fresh
mkdir -p shim
# Decompression (the verify read-back) is slowed so a write can land AFTER the
# post-tar mtime guard but BEFORE the strip; only the pre-delete re-check can
# catch it.
printf '#!/usr/bin/env bash\nfor a in "$@"; do [ "$a" = "-dcq" ] && { sleep 4; exec %s "$@"; }; done\nexec %s "$@"\n' "$ZSTD" "$ZSTD" > shim/zstd
chmod +x shim/zstd
( sleep 2; echo "a late write during verify" >> data/finished_a/finished_a_body1/finished_a_body1.1.vtu ) &
out=$(PATH="$PWD/shim:$PATH" bash "$S" --only finished_a --apply 2>&1); rc=$?
wait
echo "$out" | grep -q 'CHANGED-BEFORE-DELETE finished_a' && pass "pre-delete recheck caught the write" || fail "pre-delete recheck (got: $(echo "$out" | grep finished_a | head -2 | tr '\n' ';'))"
[ "$rc" = 8 ] && pass "exit 8" || fail "exit 8 (got $rc)"
[ "$(nvtk finished_a)" = 40 ] && pass "nothing deleted" || fail "nothing deleted (got $(nvtk finished_a))"
find arch -name '*.changed.*' | grep -q . && pass "tarball quarantined" || fail "tarball quarantined"
# L6: a stat that cannot read mtimes must hard-abort the run, not silently
# compare "" == "" and pass the changed-during-archive guard.
cd ..; fresh
mkdir -p shim2
REALSTAT="$(command -v stat)"
printf '#!/usr/bin/env bash\nfor a in "$@"; do case "$a" in *%%m*|*%%Y*) exit 1;; esac; done\nexec %s "$@"\n' "$REALSTAT" > shim2/stat
chmod +x shim2/stat
out=$(PATH="$PWD/shim2:$PATH" bash "$S" --only finished_a --apply 2>&1); rc=$?
echo "$out" | grep -q 'MTIME-UNREADABLE finished_a' && pass "empty mtime -> hard abort" || fail "empty mtime hard abort (got: $(echo "$out" | grep finished_a | head -2 | tr '\n' ';'))"
[ "$rc" = 8 ] && pass "exit 8" || fail "exit 8 (got $rc)"
[ "$(nvtk finished_a)" = 40 ] && pass "nothing deleted" || fail "nothing deleted (got $(nvtk finished_a))"
[ -z "$(find arch -name 'finished_a.tar.zst' 2>/dev/null)" ] && pass "no canonical tarball written" || fail "no canonical tarball"

echo "T21 shared-FS locking: archiver readers vs sweeper exclusivity"
cd ..; fresh
SLUGME="$(pwd -P | sed -e "s|^$HOME/||" -e 's|^/||' -e 's|/|_|g')"
# a sweeper holding sweeper.excl keeps every archiver out
mkdir -p "$LOCKROOT/sweeper.excl"
bash "$S" >/dev/null 2>&1; rc=$?
[ "$rc" = 6 ] && pass "archiver backs off while sweeper.excl held" || fail "archiver vs sweeper.excl (got $rc)"
rm -rf "$LOCKROOT/sweeper.excl"
# another archiver's checkout lock: that checkout is skipped, exit 6
mkdir -p "$LOCKROOT/checkout.$SLUGME"
out=$(bash "$S" --apply 2>&1); rc=$?
echo "$out" | grep -q 'CHECKOUT-LOCKED' && pass "locked checkout reported" || fail "locked checkout reported"
[ "$rc" = 6 ] && pass "exit 6 on locked checkout" || fail "exit 6 on locked checkout (got $rc)"
[ -z "$(ls -A arch 2>/dev/null)" ] && pass "locked checkout archived nothing" || fail "locked checkout archived nothing"
rm -rf "$LOCKROOT/checkout.$SLUGME"
# an archiver reader keeps the sweeper out
mkdir -p "$LOCKROOT/readers/archiver.test.1"
bash "$SW" >/dev/null 2>&1; rc=$?
[ "$rc" = 6 ] && pass "sweeper backs off while a reader exists" || fail "sweeper vs reader (got $rc)"
rm -rf "$LOCKROOT/readers/archiver.test.1"
# locks are released on normal exit
bash "$S" --apply >/dev/null 2>&1
[ -z "$(ls -A "$LOCKROOT/readers" 2>/dev/null)" ] && pass "reader released" || fail "reader released"
[ ! -d "$LOCKROOT/checkout.$SLUGME" ] && pass "checkout lock released" || fail "checkout lock released"
bash "$SW" >/dev/null 2>&1; rc=$?
[ "$rc" = 0 ] && pass "sweeper runs and exits 0 when unlocked" || fail "sweeper unlocked (got $rc)"
[ ! -d "$LOCKROOT/sweeper.excl" ] && pass "sweeper.excl released" || fail "sweeper.excl released"
# locks are released on usage-error exits too (trap fires on every path)
bash "$S" --quiet-hours 0 >/dev/null 2>&1
[ -z "$(ls -A "$LOCKROOT/readers" 2>/dev/null)" ] && pass "reader released on error exit" || fail "reader released on error exit"

echo; [ "$FAILED" = 0 ] && echo "ALL PASS" || { echo "SOME FAILED"; exit 1; }
