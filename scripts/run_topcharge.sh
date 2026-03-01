#!/bin/bash
# ==============================================================================
# run_topcharge.sh
# ==============================================================================
# Job script for measuring topological charge on gauge configurations.
# Applies APE smearing and measures Q as a function of smearing steps.
#
# Usage:
#   ./scripts/run_topcharge.sh <input_file>              # SU(2)
#   ./scripts/run_topcharge.sh --su3 <input_file>        # SU(3)
#   or with job scheduler:
#   sbatch scripts/run_topcharge.sh <input_file>
#
# Author: Alexander de Barros Noll
# Date: January 2026
# ==============================================================================

#SBATCH --job-name=topcharge
#SBATCH --output=logs/topcharge_%j.out
#SBATCH --error=logs/topcharge_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -e

# ==============================================================================
# Configuration
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="${PROJECT_ROOT}/build"

# Parse arguments
GAUGE_GROUP="su2"
INPUT_FILE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --su3)
            GAUGE_GROUP="su3"
            shift
            ;;
        --su2)
            GAUGE_GROUP="su2"
            shift
            ;;
        *)
            INPUT_FILE="$1"
            shift
            ;;
    esac
done

# Select executable based on gauge group
if [ "$GAUGE_GROUP" = "su3" ]; then
    EXECUTABLE="${BUILD_DIR}/bin/meas_topcharge_su3"
    GROUP_NAME="SU(3)"
else
    EXECUTABLE="${BUILD_DIR}/bin/meas_topcharge"
    GROUP_NAME="SU(2)"
fi

# ==============================================================================
# Validation
# ==============================================================================

if [ -z "$INPUT_FILE" ]; then
    echo "Usage: $0 [--su2|--su3] <input_file>"
    echo ""
    echo "Options:"
    echo "  --su2    Use SU(2) measurement (default)"
    echo "  --su3    Use SU(3) measurement"
    echo ""
    echo "Input file format:"
    echo "  config_dir        <path>    # Directory with gauge configurations"
    echo "  output_dir        <path>    # Output directory"
    echo "  beta              <value>   # Coupling constant"
    echo "  T                 <value>   # Temporal lattice extent"
    echo "  L                 <value>   # Spatial lattice extent"
    echo "  start_conf        <value>   # First configuration number"
    echo "  end_conf          <value>   # Last configuration number"
    echo "  conf_step         <value>   # Step between configurations"
    echo "  smear_steps       <value>   # Number of smearing steps"
    echo "  smear_alpha       <value>   # APE smearing parameter"
    exit 1
fi

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found: $INPUT_FILE"
    exit 1
fi

if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    echo "Please run ./scripts/build.sh first"
    exit 1
fi

# ==============================================================================
# Setup
# ==============================================================================

mkdir -p "${PROJECT_ROOT}/logs"

# ==============================================================================
# Run Measurement
# ==============================================================================

echo "=============================================="
echo "${GROUP_NAME} Topological Charge Measurement"
echo "=============================================="
echo "Gauge group: ${GROUP_NAME}"
echo "Executable:  ${EXECUTABLE}"
echo "Input file:  ${INPUT_FILE}"
echo "=============================================="
echo ""

echo "Starting measurement at $(date)"
START_TIME=$(date +%s)

"$EXECUTABLE" -i "$INPUT_FILE"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "=============================================="
echo "Measurement complete!"
echo "=============================================="
echo "Elapsed time: ${ELAPSED} seconds"
echo "=============================================="
