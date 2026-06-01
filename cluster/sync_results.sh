#!/bin/bash
# ==============================================================================
# sync_results.sh — Sync analysis data from Fuchs cluster to local machine
# ==============================================================================
# Copies only missing .dat files from both remote result locations and does not
# overwrite local files by default.
# Run this from your LOCAL machine.
#
# Usage:
#   ./cluster/sync_results.sh             # missing files only
#   ./cluster/sync_results.sh --refresh   # update/overwrite changed files
# ==============================================================================

set -euo pipefail

REMOTE="barros@fuchs.hhlr-gu.de"
REMOTE_HOME_DIR="/home/mesonqcd/barros/SU2_Calc/Topological-freezing-in-SU2-3-simulations/data/results"
REMOTE_WORK_DIR="/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/data/results"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
LOCAL_HOME_DIR="$PROJECT_DIR/data/results_home"
LOCAL_WORK_DIR="$PROJECT_DIR/data/results_work"

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

sync_dat_results_dir() {
  local remote_dir="$1"
  local local_dir="$2"
  local label="$3"

  mkdir -p "$local_dir"

  echo "Syncing $label .dat results:"
  echo "  $REMOTE:$remote_dir/ -> $local_dir/"

  rsync -avm \
    $RSYNC_EXISTING_FLAG \
    --include="*/" \
    --include="*.dat" \
    --exclude="*" \
    "$REMOTE:$remote_dir/" \
    "$local_dir/"
}

sync_dat_results_dir "$REMOTE_HOME_DIR" "$LOCAL_HOME_DIR" "home"
sync_dat_results_dir "$REMOTE_WORK_DIR" "$LOCAL_WORK_DIR" "work"
