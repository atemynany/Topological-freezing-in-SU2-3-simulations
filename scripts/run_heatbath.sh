#!/bin/bash
# ==============================================================================
# run_heatbath.sh
# ==============================================================================
# Job script for running Monte Carlo heatbath configuration generation.
# Generates gauge field configurations for SU(2) lattice gauge theory.
#
# Usage:
#   ./scripts/run_heatbath.sh <config_file>
#   or with job scheduler:
#   sbatch scripts/run_heatbath.sh <config_file>
#
# Author: Alexander de Barros Noll
# Date: January 2026
# ==============================================================================

#SBATCH --job-name=su2_heatbath
#SBATCH --output=logs/heatbath_%j.out
#SBATCH --error=logs/heatbath_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -e

# ==============================================================================
# Configuration
# ==============================================================================

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="${PROJECT_ROOT}/build"
EXECUTABLE="${BUILD_DIR}/bin/mc_heatbath"

# Load configuration file or use defaults
CONFIG_FILE="${1:-}"

if [ -n "$CONFIG_FILE" ] && [ -f "$CONFIG_FILE" ]; then
    echo "Loading configuration from: $CONFIG_FILE"
    source "$CONFIG_FILE"
else
    # Default parameters
    OUTPUT_DIR="${PROJECT_ROOT}/data/configs"
    BETA="2.4"
    T="16"
    L="16"
    SEED="12345"
    START_TYPE="cold"      # cold or hot
    NUM_SWEEPS="10000"     # Total number of MC sweeps
    OUTPUT_INTERVAL="500"  # Save configuration every N sweeps
    BOUNDARY="periodic"    # periodic or open
fi

# ==============================================================================
# Setup
# ==============================================================================

# Create output and log directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${PROJECT_ROOT}/logs"

# Check if executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    echo "Please run ./scripts/build.sh first"
    exit 1
fi

# ==============================================================================
# Run Simulation
# ==============================================================================

echo "=============================================="
echo "SU(2) Heatbath Monte Carlo Simulation"
echo "=============================================="
echo "Output directory: ${OUTPUT_DIR}"
echo "Beta: ${BETA}"
echo "Lattice: ${T} x ${L}^3"
echo "Seed: ${SEED}"
echo "Start type: ${START_TYPE}"
echo "MC sweeps: ${NUM_SWEEPS}"
echo "Output interval: ${OUTPUT_INTERVAL}"
echo "Boundary conditions: ${BOUNDARY}"
echo "=============================================="
echo ""

# Create output subdirectory with parameters
RUN_DIR="${OUTPUT_DIR}/T${T}_L${L}_b${BETA}"
mkdir -p "${RUN_DIR}"

# Record parameters
cat > "${RUN_DIR}/parameters.txt" << EOF
# Simulation Parameters
# Generated: $(date)
beta=${BETA}
T=${T}
L=${L}
seed=${SEED}
start_type=${START_TYPE}
num_sweeps=${NUM_SWEEPS}
output_interval=${OUTPUT_INTERVAL}
boundary=${BOUNDARY}
EOF

# Run the simulation
echo "Starting simulation at $(date)"
START_TIME=$(date +%s)

"$EXECUTABLE" "${RUN_DIR}/" "$BETA" "$T" "$L" "$SEED" "$START_TYPE" "$NUM_SWEEPS" "$OUTPUT_INTERVAL" "$BOUNDARY"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "=============================================="
echo "Simulation complete!"
echo "=============================================="
echo "Elapsed time: ${ELAPSED} seconds"
echo "Configurations saved in: ${RUN_DIR}"
echo "=============================================="
