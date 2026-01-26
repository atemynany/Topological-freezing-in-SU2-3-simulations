#!/bin/bash
# ==============================================================================
# run_full_analysis.sh
# ==============================================================================
# Complete workflow script: Generate configurations, measure topological charge,
# and create plots.
#
# Usage:
#   ./scripts/run_full_analysis.sh <parameter_file>
#
# Author: Alexander de Barros Noll
# Date: January 2026
# ==============================================================================

#SBATCH --job-name=su2_full_analysis
#SBATCH --output=logs/full_analysis_%j.out
#SBATCH --error=logs/full_analysis_%j.err
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -e

# ==============================================================================
# Configuration
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="${PROJECT_ROOT}/build"

# Load parameter file or use defaults
PARAM_FILE="${1:-}"

if [ -n "$PARAM_FILE" ] && [ -f "$PARAM_FILE" ]; then
    echo "Loading parameters from: $PARAM_FILE"
    source "$PARAM_FILE"
else
    # Default parameters for a test run
    BETA="2.4"
    T="16"
    L="16"
    SEED="12345"
    START_TYPE="cold"
    NUM_SWEEPS="5000"
    OUTPUT_INTERVAL="500"
    BOUNDARY="periodic"
    SMEAR_STEPS="200"
    SMEAR_INTERVAL="20"
    SMEAR_ALPHA="0.5"
fi

# Derived paths
DATA_DIR="${PROJECT_ROOT}/data"
CONFIG_DIR="${DATA_DIR}/configs/T${T}_L${L}_b${BETA}"
RESULTS_DIR="${DATA_DIR}/results/T${T}_L${L}_b${BETA}"
PLOTS_DIR="${PROJECT_ROOT}/plots/T${T}_L${L}_b${BETA}"

# ==============================================================================
# Setup
# ==============================================================================

mkdir -p "${CONFIG_DIR}"
mkdir -p "${RESULTS_DIR}"
mkdir -p "${PLOTS_DIR}"
mkdir -p "${PROJECT_ROOT}/logs"

echo "=============================================="
echo "SU(2) Full Analysis Pipeline"
echo "=============================================="
echo "Parameters:"
echo "  Beta: ${BETA}"
echo "  Lattice: ${T} x ${L}^3"
echo "  MC sweeps: ${NUM_SWEEPS}"
echo "  Smearing steps: ${SMEAR_STEPS}"
echo "=============================================="
echo ""

# ==============================================================================
# Step 1: Build
# ==============================================================================

echo "[1/4] Building project..."
"${SCRIPT_DIR}/build.sh" release
echo ""

# ==============================================================================
# Step 2: Generate Configurations
# ==============================================================================

echo "[2/4] Generating gauge configurations..."
echo "      This may take a while..."

"${BUILD_DIR}/bin/mc_heatbath" \
    "${CONFIG_DIR}/" \
    "$BETA" \
    "$T" \
    "$L" \
    "$SEED" \
    "$START_TYPE" \
    "$NUM_SWEEPS" \
    "$OUTPUT_INTERVAL" \
    "$BOUNDARY"

echo "Configurations saved to: ${CONFIG_DIR}"
echo ""

# ==============================================================================
# Step 3: Measure Topological Charge
# ==============================================================================

echo "[3/4] Measuring topological charge..."

# Create input file for measurement
INPUT_FILE="${RESULTS_DIR}/topcharge_input.txt"
cat > "$INPUT_FILE" << EOF
# Topological charge measurement input
# Generated: $(date)
config_dir          ${CONFIG_DIR}
output_dir          ${RESULTS_DIR}
beta                ${BETA}
T                   ${T}
L                   ${L}
start_conf          ${OUTPUT_INTERVAL}
end_conf            ${NUM_SWEEPS}
conf_step           ${OUTPUT_INTERVAL}
smear_steps         ${SMEAR_STEPS}
smear_output_interval ${SMEAR_INTERVAL}
smear_alpha         ${SMEAR_ALPHA}
seed                ${SEED}
EOF

"${BUILD_DIR}/bin/meas_topcharge" -i "$INPUT_FILE"

echo "Results saved to: ${RESULTS_DIR}"
echo ""

# ==============================================================================
# Step 4: Generate Plots
# ==============================================================================

echo "[4/4] Generating plots..."

# Find the output data file
DATA_FILE=$(ls "${RESULTS_DIR}"/topcharge_*.dat 2>/dev/null | head -1)

if [ -n "$DATA_FILE" ] && [ -f "$DATA_FILE" ]; then
    python3 "${PROJECT_ROOT}/analysis/plot_topcharge.py" \
        --input "$DATA_FILE" \
        --output "${PLOTS_DIR}" \
        --title "SU(2) ${T}x${L}^3 beta=${BETA}"
    
    echo "Plots saved to: ${PLOTS_DIR}"
else
    echo "Warning: No data file found for plotting"
fi

echo ""
echo "=============================================="
echo "Full Analysis Complete!"
echo "=============================================="
echo "Outputs:"
echo "  Configurations: ${CONFIG_DIR}"
echo "  Results: ${RESULTS_DIR}"
echo "  Plots: ${PLOTS_DIR}"
echo "=============================================="
