#!/bin/bash
# ==============================================================================
# Lattice QCD Simulation - Main Entry Point
# ==============================================================================
# This is the main script to run lattice gauge theory simulations.
# It wraps the beta_scan workflow which handles:
#   1. Building the project
#   2. Running MC heatbath simulations for each beta value
#   3. Auto-detecting thermalization
#   4. Measuring topological charge
#   5. Running analysis and generating plots
#
# Usage:
#   ./run_simulation.sh [--su2|--su3] [options]
#
# Options:
#   --su2           Use SU(2) gauge group (default)
#   --su3           Use SU(3) gauge group
#   --dry-run       Show what would be run without executing
#   --skip-build    Skip the build step (use existing binaries)
#   --open          Use open boundary conditions (uses *_open.txt params)
#   --beta-file     Custom file with beta values (one per line)
#   --params-file   Custom base parameters file
#   -h, --help      Show this help message
#
# Examples:
#   ./run_simulation.sh                    # Run SU(2) with default betas
#   ./run_simulation.sh --su3              # Run SU(3) with default betas
#   ./run_simulation.sh --su2 --open       # Run SU(2) with open boundaries
#   ./run_simulation.sh --dry-run          # Preview without executing
#
# Output:
#   Results are saved to data/results/T{T}_L{L}_b{beta}_{boundary}_seed{seed}/
#   Analysis plots go to output/figures_beta_scan_{su2|su3}/
#
# Configuration:
#   Edit input/base_params_{su2|su3}.txt to change simulation parameters
#   Edit input/beta_scan_{su2|su3}.txt to change beta values
#
# ==============================================================================

set -e

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default settings
GAUGE_GROUP="su2"
OPEN_BC=false
EXTRA_ARGS=()

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --su2)
            GAUGE_GROUP="su2"
            shift
            ;;
        --su3)
            GAUGE_GROUP="su3"
            shift
            ;;
        --open)
            OPEN_BC=true
            shift
            ;;
        -h|--help)
            head -45 "$0" | tail -40
            exit 0
            ;;
        *)
            EXTRA_ARGS+=("$1")
            shift
            ;;
    esac
done

# Select appropriate parameter files
if [ "$OPEN_BC" = true ]; then
    PARAMS_FILE="input/base_params_${GAUGE_GROUP}_open.txt"
    BETA_FILE="input/beta_scan_${GAUGE_GROUP}_open.txt"
else
    PARAMS_FILE="input/base_params_${GAUGE_GROUP}.txt"
    BETA_FILE="input/beta_scan_${GAUGE_GROUP}.txt"
fi

# Check if files exist, fall back to non-open versions if needed
if [ ! -f "$BETA_FILE" ]; then
    BETA_FILE="input/beta_scan_${GAUGE_GROUP}.txt"
fi

# Run the beta scan workflow
exec "$SCRIPT_DIR/scripts/run_beta_scan.sh" \
    "--${GAUGE_GROUP}" \
    --params-file "$PARAMS_FILE" \
    --beta-file "$BETA_FILE" \
    "${EXTRA_ARGS[@]}"
