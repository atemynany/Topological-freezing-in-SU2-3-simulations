#!/bin/bash
# ==============================================================================
# sync_results.sh — Sync analysis data from Fuchs cluster to local machine
# ==============================================================================
# Copies only .dat files (topcharge, plaquette, etc.) — skips binary configs.
# Run this from your LOCAL machine.
#
# Usage:
#   ./cluster/sync_results.sh
# ==============================================================================

REMOTE="barros@fuchs.hhlr-gu.de"
#REMOTE_DIR="/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/data/results/"
REMOTE_DIR="/home/mesonqcd/barros/SU2_Calc/Topological-freezing-in-SU2-3-simulations/data/results"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCAL_DIR="$(dirname "$SCRIPT_DIR")/data/results/"

rsync -avm \
  --include="*/" \
  --include="*.dat" \
  --include="input.txt" \
  --exclude="*" \
  "$REMOTE:$REMOTE_DIR" \
  "$LOCAL_DIR"
