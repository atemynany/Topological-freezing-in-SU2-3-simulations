#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
PYTHON="/home/alex/miniconda3/envs/master_thesis/bin/python"

PLAQ_FILE=""
Q_FILE=""
SMEAR=""
THERM="0"

while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--plaq) PLAQ_FILE="$2"; shift 2;;
        -q|--topcharge) Q_FILE="$2"; shift 2;;
        --smear) SMEAR="$2"; shift 2;;
        --therm) THERM="$2"; shift 2;;
        *) echo "Unknown: $1"; exit 1;;
    esac
done

ARGS=""
[ -n "$PLAQ_FILE" ] && ARGS="$ARGS --plaq $PLAQ_FILE --therm $THERM"
[ -n "$Q_FILE" ] && ARGS="$ARGS --topcharge $Q_FILE"
[ -n "$SMEAR" ] && ARGS="$ARGS --smear $SMEAR"

$PYTHON "${PROJECT_ROOT}/analysis/autocorrelation.py" $ARGS
