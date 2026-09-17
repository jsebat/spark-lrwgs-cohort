#!/usr/bin/env bash
# wait_job.sh <jobid> [max_seconds] [poll_seconds]
# Poll a Slurm job until it leaves the queue, then print its terminal state.
# JS rule 2026-09-05: for a job expected to take <= 5 min, poll every minute
# rather than reporting a state that is already stale by the time it is read.
# Default: poll every 60 s, give up after 15 min (then say so - never claim
# an unknown state is "pending").
set -uo pipefail
J="${1:?usage: wait_job.sh <jobid> [max_seconds] [poll_seconds]}"
MAX="${2:-900}"
POLL="${3:-60}"
t=0
while [ "$t" -lt "$MAX" ]; do
  st=$(sacct -j "$J" -X -n -P --format=State 2>/dev/null | head -1 | tr -d ' ')
  case "$st" in
    COMPLETED|FAILED|CANCELLED*|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|BOOT_FAIL|DEADLINE)
      el=$(sacct -j "$J" -X -n -P --format=Elapsed 2>/dev/null | head -1 | tr -d ' ')
      echo "JOB $J TERMINAL: $st after ${el:-?} (waited ${t}s)"
      exit 0
      ;;
  esac
  sleep "$POLL"
  t=$((t + POLL))
done
st=$(sacct -j "$J" -X -n -P --format=State 2>/dev/null | head -1 | tr -d ' ')
case "$st" in
  COMPLETED|FAILED|CANCELLED*|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL) echo "JOB $J TERMINAL: $st (reached at the ${MAX}s deadline)"; exit 0;;
esac
echo "JOB $J STILL ${st:-UNKNOWN} after ${MAX}s - not terminal, reporting as still running"
exit 1
