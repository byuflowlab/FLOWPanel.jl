#!/bin/bash
# Poll 018 ladder jobs every 10 min. Exit (waking the agent) on:
#  - any job reaching a terminal state
#  - first RUNNING of any job still needing verification (all 5 initially)
#  - a free mgh-1-1 GH200 while the merge-off NT36 (13518484) is PENDING
# 13518931 = NT144 nm resubmit (3rd attempt, mgh-1-2, after NODE_FAILs of
# 13518482 and 13518894 both on mgh-1-2; mgh-1-1 still leaked 61.6 GB)
# 13518932 = NT144 nm eng/h200 backup (run name *_nm_eng).
# m13h doubles (qos=gpu, -C intel, run names *_m13h): 13518979 NT72 rlxf,
# 13518980 NT144 rlxf (trimmed to 48 cpus for the 192-cpu user cap),
# 13518981 NT144 nm. Ryan: keep whichever copy finishes, cancel losers.
# 13518483 (NT72 nm) COMPLETED 07:55 exit 0:0 — CT̄(10-rev) 0.0719342 — dropped
JOBS="13518931,13518932,13518484,13518485,13518486,13518979,13518980,13518981"
NT36=13518484
# jobs whose start still needs metadata/banner verification (space-separated)
# excluded (verification handled by dedicated polls): 13518931, 13518484,
# 13518979
UNVERIFIED="13518932 13518485 13518486 13518981"
while true; do
  OUT=$(ssh orc 'bash -lc "sacct -j '"$JOBS"' -X -P -n -o JobID,State; echo ---SINFO---; sinfo -p mgh -N -h -O NodeList,GresUsed"' 2>/dev/null | sed $'s/\x1b\\[[0-9;]*m//g')
  echo "== $(date '+%F %T') =="
  echo "$OUT"
  SACCT=$(echo "$OUT" | sed -n '/---SINFO---/q;p')
  SINFO=$(echo "$OUT" | sed -n '/---SINFO---/,$p' | tail -n +2)
  # terminal state?
  if echo "$SACCT" | grep -Eq 'FAILED|COMPLETED|CANCELLED|TIMEOUT|OUT_OF_ME|NODE_FAIL|PREEMPTED'; then
    echo "WAKE: terminal state detected"; exit 0
  fi
  # unverified job now RUNNING?
  for j in $UNVERIFIED; do
    if echo "$SACCT" | grep -q "^$j|RUNNING"; then
      echo "WAKE: job $j RUNNING (needs verification)"; exit 0
    fi
  done
  # mgh-1-1 hop DISABLED 2026-08-31: GPU has ~61.6 GB leaked memory with no
  # processes (caused 13518892 OOM-fail); do not hop anything there.
  sleep 1800
done
