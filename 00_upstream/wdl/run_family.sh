#!/bin/bash
# Per-family launcher.  Submits the miniwdl driver AND its watchdog as Slurm
# jobs (EXPANSE.md rule 1: no workflow drivers on login nodes).
# Double-launch guard uses squeue, not pgrep, so it works from either login node.
#   run_family.sh [--driver-only] <FAMILY_ID> [RUNROOT]
set -uo pipefail
# --- site configuration: every path and account comes from config/upstream.env (copy upstream.env.example) ---
UP_HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
UP_ENV="${UPSTREAM_ENV:-$UP_HERE/config/upstream.env}"
[ -f "$UP_ENV" ] && source "$UP_ENV"
: "${WDL_DIR:?set WDL_DIR (HiFi-human-WGS-WDL checkout, v3.3.1) in config/upstream.env}"
: "${INPUTS_DIR:?}" "${RUNROOT:?}" "${SLURM_ACCOUNT:?}" "${SLURM_PARTITION:?}" "${SLURM_QOS:?}"
NOTIFY_EMAIL="${NOTIFY_EMAIL:-}"
DRIVER_ONLY=0
[ "${1:-}" = "--driver-only" ] && { DRIVER_ONLY=1; shift; }
FAMILY="${1:?usage: run_family.sh [--driver-only] <FAMILY_ID> [RUNROOT]}"
RUNROOT="${2:-${RUNROOT:?}}"
C="$UP_HERE/wdl"
WDL="$WDL_DIR/workflows/family.wdl"
INPUTS="$INPUTS_DIR/$FAMILY.inputs.json"
RUNDIR="$RUNROOT/run_$FAMILY"
ACCT="$SLURM_ACCOUNT"; PART="$SLURM_PARTITION"; QOS="$SLURM_QOS"
DRV_TIME=${DRV_TIME:-48:00:00}; DRV_MEM=${DRV_MEM:-4G}
WD_TIME=${WD_TIME:-50:00:00};   WD_MEM=${WD_MEM:-2G}   # the watchdog must OUTLIVE a driver killed by TIMEOUT to restart it; equal walltimes cannot

[ -f "$INPUTS" ] || { echo "FATAL: no inputs at $INPUTS (run make_inputs.py)" >&2; exit 2; }
[ -f "$WDL" ]    || { echo "FATAL: no WDL at $WDL" >&2; exit 2; }

# squeue-based double-launch guard (replaces pgrep; login-node independent)
existing=$(squeue -u "$USER" -h -n "drv_$FAMILY" -o "%i %T" 2>/dev/null)
if [ -n "$existing" ]; then
  echo "FATAL: $FAMILY already has a driver job: $existing" >&2; exit 3
fi

# Global cap: refuse to add another driver past MAX_DRIVERS, whoever calls us.
# Protects against a stale pgrep-era submit_wave.sh that cannot see drv_* jobs
# and would otherwise launch every remaining family in a loop.
# The cap applies to NEW families only. A restart or resume of a family that
# already has work on disk is EXEMPT -- refusing one strands completed work and
# burns the watchdog's restart budget, which is what happened to <family> on
# 2026-08-31 (D34). The per-family double-launch guard above still applies, so
# an exemption can never produce two drivers for the same family.
IS_RESTART=0
[ "$DRIVER_ONLY" -eq 1 ] && IS_RESTART=1
ls -d "$RUNDIR"/*_humanwgs_family >/dev/null 2>&1 && IS_RESTART=1

if [ "$IS_RESTART" -eq 1 ]; then
  echo "$FAMILY: restart/resume of existing work -- exempt from MAX_DRIVERS"
else
  MAX_DRIVERS=${MAX_DRIVERS:-3}
  ndrv=$(squeue -u "$USER" -h -o "%j" 2>/dev/null | grep -c '^drv_')
  if [ "$ndrv" -ge "$MAX_DRIVERS" ]; then
    echo "REFUSING $FAMILY: $ndrv driver jobs already active (MAX_DRIVERS=$MAX_DRIVERS)" >&2
    exit 5
  fi
fi

mkdir -p "$RUNDIR"
rm -f "$RUNDIR"/.watchdog.* "$RUNDIR/.stall_notified"

DRV=$(sbatch --parsable \
  --job-name="drv_$FAMILY" --account=$ACCT --partition=$PART --qos=$QOS \
  --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=$DRV_MEM --time=$DRV_TIME \
  --output="$RUNDIR/driver.%j.out" ${NOTIFY_EMAIL:+--mail-type=END,FAIL --mail-user=$NOTIFY_EMAIL} \
  --export=ALL,FAMILY="$FAMILY",RUNROOT="$RUNROOT" \
  "$C/driver_family.sb") || { echo "FATAL: driver sbatch failed" >&2; exit 4; }
echo "$(date) submitted driver job $DRV for $FAMILY"

if [ "$DRIVER_ONLY" -eq 0 ]; then
  WD=$(sbatch --parsable \
    --job-name="wd_$FAMILY" --account=$ACCT --partition=$PART --qos=$QOS \
    --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=$WD_MEM --time=$WD_TIME \
    --output="$RUNDIR/watchdog.%j.out" \
    --export=ALL,FAMILY="$FAMILY",RUNROOT="$RUNROOT" \
    "$C/watchdog_family.sb") && echo "$(date) submitted watchdog job $WD for $FAMILY" \
    || echo "WARNING: watchdog sbatch failed for $FAMILY" >&2
fi
