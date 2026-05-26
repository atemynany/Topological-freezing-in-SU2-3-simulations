#!/bin/bash
# ==============================================================================
# sync_results.sh — Sync analysis data from Fuchs cluster to local machine
# ==============================================================================
# Copies only missing .dat files (topcharge, plaquette, etc.) and input.txt
# files — skips binary configs and does not overwrite local files by default.
# Run this from your LOCAL machine.
#
# Usage:
#   ./cluster/sync_results.sh             # missing files only
#   ./cluster/sync_results.sh --refresh   # update/overwrite changed files
# ==============================================================================

set -euo pipefail

REMOTE="barros@fuchs.hhlr-gu.de"
#REMOTE_DIR="/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/data/results/"
REMOTE_DIR="/home/mesonqcd/barros/SU2_Calc/Topological-freezing-in-SU2-3-simulations/data/results/"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCAL_DIR="$(dirname "$SCRIPT_DIR")/data/results/"

RSYNC_EXISTING_FLAG="--ignore-existing"
if [ "${1:-}" = "--refresh" ]; then
  RSYNC_EXISTING_FLAG=""
elif [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
  sed -n '1,12p' "$0"
  exit 0
elif [ $# -gt 0 ]; then
  echo "Unknown option: $1" >&2
  echo "Usage: $0 [--refresh]" >&2
  exit 1
fi

rsync -avm \
  $RSYNC_EXISTING_FLAG \
  --include="*/" \
  --include="*.dat" \
  --include="input.txt" \
  --exclude="*" \
  "$REMOTE:$REMOTE_DIR" \
  "$LOCAL_DIR"
