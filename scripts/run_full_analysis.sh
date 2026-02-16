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

# Parse command-line arguments
INPUT_FILE=""
PARAM_FILE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        *)
            PARAM_FILE="$1"
            shift
            ;;
    esac
done

# Helper function to parse input file
parse_input_file() {
    local file="$1"
    local key="$2"
    grep -E "^${key}[[:space:]]" "$file" | awk '{print $2}' | head -1
}

# Load from input file (-i flag) or parameter file (sourced) or use defaults
if [ -n "$INPUT_FILE" ] && [ -f "$INPUT_FILE" ]; then
    echo "Loading parameters from input file: $INPUT_FILE"
    # Parse the unified input file format
    OUTPUT_DIR=$(parse_input_file "$INPUT_FILE" "output_dir")
    CONFIG_DIR_FROM_FILE=$(parse_input_file "$INPUT_FILE" "config_dir")
    T=$(parse_input_file "$INPUT_FILE" "T")
    L=$(parse_input_file "$INPUT_FILE" "L")
    BETA=$(parse_input_file "$INPUT_FILE" "beta")
    SEED=$(parse_input_file "$INPUT_FILE" "seed")
    START_TYPE=$(parse_input_file "$INPUT_FILE" "start_type")
    BOUNDARY=$(parse_input_file "$INPUT_FILE" "boundary")
    NUM_SWEEPS=$(parse_input_file "$INPUT_FILE" "num_sweeps")
    OUTPUT_INTERVAL=$(parse_input_file "$INPUT_FILE" "save_interval")
    SMEAR_STEPS=$(parse_input_file "$INPUT_FILE" "smear_steps")
    SMEAR_INTERVAL=$(parse_input_file "$INPUT_FILE" "smear_interval")
    SMEAR_ALPHA=$(parse_input_file "$INPUT_FILE" "smear_alpha")
    START_CONF=$(parse_input_file "$INPUT_FILE" "start_conf")
    END_CONF=$(parse_input_file "$INPUT_FILE" "end_conf")
    CONF_STEP=$(parse_input_file "$INPUT_FILE" "conf_step")
    USE_INPUT_FILE=true
elif [ -n "$PARAM_FILE" ] && [ -f "$PARAM_FILE" ]; then
    echo "Loading parameters from: $PARAM_FILE"
    source "$PARAM_FILE"
    USE_INPUT_FILE=false
else
    echo "No parameter file provided, using defaults"
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
    USE_INPUT_FILE=false
fi

# Derived paths (use config_dir from file if provided, otherwise construct it)
DATA_DIR="${PROJECT_ROOT}/data"
if [ -n "$CONFIG_DIR_FROM_FILE" ]; then
    CONFIG_DIR="${PROJECT_ROOT}/${CONFIG_DIR_FROM_FILE}"
else
    CONFIG_DIR="${DATA_DIR}/configs/T${T}_L${L}_b${BETA}"
fi
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

# Create a temporary input file for mc_heatbath (it now uses -i flag)
HEATBATH_INPUT="${PROJECT_ROOT}/logs/heatbath_input_temp.txt"
cat > "$HEATBATH_INPUT" << EOF
# Heatbath input file (auto-generated)
# Generated: $(date)
output_dir          ${PROJECT_ROOT}/output/
config_dir          ${CONFIG_DIR}/
T                   ${T}
L                   ${L}
beta                ${BETA}
seed                ${SEED}
start_type          ${START_TYPE}
boundary            ${BOUNDARY}
num_sweeps          ${NUM_SWEEPS}
save_interval       ${OUTPUT_INTERVAL}
EOF

"${BUILD_DIR}/bin/mc_heatbath" -i "$HEATBATH_INPUT"

echo "Configurations saved to: ${CONFIG_DIR}"
echo ""

# ==============================================================================
# Step 3: Measure Topological Charge
# ==============================================================================

echo "[3/4] Measuring topological charge..."

# Create input file for measurement
MEAS_INPUT_FILE="${RESULTS_DIR}/topcharge_input.txt"

# Use values from input file if provided, otherwise compute from generation params
if [ "$USE_INPUT_FILE" = true ] && [ -n "$START_CONF" ]; then
    START_CONF_VAL="$START_CONF"
    END_CONF_VAL="${END_CONF:-$NUM_SWEEPS}"
    CONF_STEP_VAL="${CONF_STEP:-$OUTPUT_INTERVAL}"
else
    START_CONF_VAL="$OUTPUT_INTERVAL"
    END_CONF_VAL="$NUM_SWEEPS"
    CONF_STEP_VAL="$OUTPUT_INTERVAL"
fi

cat > "$MEAS_INPUT_FILE" << EOF
# Topological charge measurement input
# Generated: $(date)
config_dir          ${CONFIG_DIR}
output_dir          ${RESULTS_DIR}
beta                ${BETA}
T                   ${T}
L                   ${L}
start_conf          ${START_CONF_VAL}
end_conf            ${END_CONF_VAL}
conf_step           ${CONF_STEP_VAL}
smear_steps         ${SMEAR_STEPS}
smear_interval      ${SMEAR_INTERVAL}
smear_alpha         ${SMEAR_ALPHA}
seed                ${SEED}
EOF

"${BUILD_DIR}/bin/meas_topcharge" -i "$MEAS_INPUT_FILE"

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
