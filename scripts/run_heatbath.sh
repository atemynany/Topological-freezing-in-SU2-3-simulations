#!/bin/bash
# ==============================================================================
# run_heatbath.sh
# ==============================================================================
# Standalone heatbath runner: generates gauge configs and shows thermalization
# plots (plaquette + topological charge) for a single beta value.
#
# Usage:
#   ./scripts/run_heatbath.sh [--su2|--su3] [--beta <val>] [--params-file <file>]
#
# Options:
#   --su2            SU(2) gauge group (default)
#   --su3            SU(3) gauge group
#   --beta <val>     Override beta value
#   --params-file    Override base params file
#   --skip-build     Skip build step
#
#SBATCH --job-name=su2_heatbath
#SBATCH --output=logs/heatbath_%j.out
#SBATCH --error=logs/heatbath_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="${PROJECT_DIR}/build"
RESULTS_DIR="${PROJECT_DIR}/data/results"
CONDA_ENV="master_thesis"

GAUGE_GROUP="su2"
BETA_OVERRIDE=""
PARAMS_FILE=""
SKIP_BUILD=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --su2)        GAUGE_GROUP="su2";  shift ;;
        --su3)        GAUGE_GROUP="su3";  shift ;;
        --beta)       BETA_OVERRIDE="$2"; shift 2 ;;
        --params-file) PARAMS_FILE="$2";  shift 2 ;;
        --skip-build) SKIP_BUILD=true;    shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [ "$GAUGE_GROUP" = "su2" ]; then
    PARAMS_FILE="${PARAMS_FILE:-$PROJECT_DIR/input/base_params_su2.txt}"
    HEATBATH_BIN="mc_heatbath"
    PLAQ_FILE="plaquette.dat"
    SU3_FLAG=""
else
    PARAMS_FILE="${PARAMS_FILE:-$PROJECT_DIR/input/base_params_su3.txt}"
    HEATBATH_BIN="mc_heatbath_su3"
    PLAQ_FILE="plaquette_su3.dat"
    SU3_FLAG="--su3"
fi

EXECUTABLE="${BUILD_DIR}/bin/${HEATBATH_BIN}"

# Helper: read a parameter from the params file
read_param() {
    grep -E "^${1}\s+" "$2" | awk '{print $2}' | head -1
}

# Generate unique seed
generate_seed() {
    python3 -c "import time; print(int(time.time()*1e9))" 2>/dev/null | cut -c1-10 \
        || echo $(($(date +%s) * $$ % 2147483647))
}

# ==============================================================================
# Build
# ==============================================================================
if [ "$SKIP_BUILD" = false ]; then
    echo "[INFO] Building $HEATBATH_BIN..."
    mkdir -p "${BUILD_DIR}"
    cd "${BUILD_DIR}"
    cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
    make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4) "$HEATBATH_BIN" 2>&1 | tail -3
    cd "${PROJECT_DIR}"
    echo "[OK] Build complete"
fi

if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found: $EXECUTABLE"
    echo "Run ./scripts/build.sh first or use --skip-build if already built."
    exit 1
fi

# ==============================================================================
# Read parameters
# ==============================================================================
T=$(read_param "T" "$PARAMS_FILE")
L=$(read_param "L" "$PARAMS_FILE")
BETA="${BETA_OVERRIDE:-$(read_param "beta" "$PARAMS_FILE")}"
# Fall back to reading first beta from beta scan file if not in params
if [ -z "$BETA" ]; then
    BETA_SCAN_FILE="${PROJECT_DIR}/input/beta_scan_${GAUGE_GROUP}.txt"
    BETA=$(grep -v '^#' "$BETA_SCAN_FILE" | grep -v '^$' | head -1)
fi
START_TYPE=$(read_param "start_type" "$PARAMS_FILE")
BOUNDARY=$(read_param "boundary" "$PARAMS_FILE")
NUM_SWEEPS=$(read_param "num_sweeps" "$PARAMS_FILE")
SAVE_INTERVAL=$(read_param "save_interval" "$PARAMS_FILE")
OVERRELAX=$(read_param "overrelax_steps" "$PARAMS_FILE")
END_CONF=$(read_param "end_conf" "$PARAMS_FILE")
CONF_STEP=$(read_param "conf_step" "$PARAMS_FILE")
SMEAR_STEPS=$(read_param "smear_steps" "$PARAMS_FILE")
SMEAR_ALPHA=$(read_param "smear_alpha" "$PARAMS_FILE")
EXCLUDE_BC=$(read_param "exclude_boundary_slices" "$PARAMS_FILE")

SEED=$(generate_seed)
BOUNDARY="${BOUNDARY:-periodic}"
RUN_NAME="T${T}_L${L}_b${BETA}_${BOUNDARY}_seed${SEED}"
RUN_DIR="${RESULTS_DIR}/${RUN_NAME}"
RUN_OUTPUT_DIR="${RUN_DIR}/output"
RUN_CONFIG_DIR="${RUN_DIR}/configs"

echo "=============================================="
echo "${GAUGE_GROUP^^} Heatbath - Config Generation"
echo "=============================================="
echo "Params file:   $PARAMS_FILE"
echo "Lattice:       ${T} x ${L}^3"
echo "Beta:          ${BETA}"
echo "Seed:          ${SEED}"
echo "Boundary:      ${BOUNDARY}"
echo "Sweeps:        ${NUM_SWEEPS}"
echo "Save interval: ${SAVE_INTERVAL}"
echo "Output:        ${RUN_DIR}"
echo "=============================================="

mkdir -p "${RUN_OUTPUT_DIR}"
mkdir -p "${RUN_CONFIG_DIR}"

# ==============================================================================
# Create input file
# ==============================================================================
INPUT_FILE="${RUN_DIR}/input.txt"
cat > "$INPUT_FILE" << EOF
# Auto-generated input - run_heatbath.sh
# Generated: $(date)

T                   ${T:-16}
L                   ${L:-16}
beta                ${BETA}
seed                ${SEED}
start_type          ${START_TYPE:-cold}
boundary            ${BOUNDARY}
num_sweeps          ${NUM_SWEEPS:-1000}
save_interval       ${SAVE_INTERVAL:-10}
EOF

if [ "$GAUGE_GROUP" = "su3" ] && [ -n "$OVERRELAX" ]; then
    echo "overrelax_steps     $OVERRELAX" >> "$INPUT_FILE"
fi

cat >> "$INPUT_FILE" << EOF

start_conf          10
end_conf            ${END_CONF:-1000}
conf_step           ${CONF_STEP:-10}
smear_steps         ${SMEAR_STEPS:-20}
smear_alpha         ${SMEAR_ALPHA:-0.3}
exclude_boundary_slices ${EXCLUDE_BC:-0}

output_dir          ${RUN_OUTPUT_DIR}/
config_dir          ${RUN_CONFIG_DIR}/
EOF

if [ "$GAUGE_GROUP" = "su3" ]; then
    echo "output_file         ${RUN_OUTPUT_DIR}/topcharge_su3.dat" >> "$INPUT_FILE"
fi

# ==============================================================================
# Run heatbath
# ==============================================================================
echo ""
echo "[STEP] Running heatbath..."
START_TIME=$(date +%s)

"$EXECUTABLE" -i "$INPUT_FILE" 2>&1 | tee "${RUN_DIR}/heatbath.log"

END_TIME=$(date +%s)
echo ""
echo "[OK] Heatbath complete in $((END_TIME - START_TIME))s"
echo "     Configs: ${RUN_CONFIG_DIR}"
echo "     Plaquette: ${RUN_OUTPUT_DIR}/${PLAQ_FILE}"

# ==============================================================================
# Plot thermalization
# ==============================================================================
echo ""
echo "[STEP] Generating thermalization plot..."

_run_python() {
    if command -v conda &> /dev/null && conda run -n "$CONDA_ENV" python3 --version &> /dev/null 2>&1; then
        conda run -n "$CONDA_ENV" python3 "$@"
    elif command -v python3 &> /dev/null; then
        python3 "$@"
    else
        echo "[WARN] Python not found, skipping plot"
        return 1
    fi
}

# ==============================================================================
# Save run info
# ==============================================================================
cat > "${RUN_DIR}/run_info.txt" << EOF
# Run Information - run_heatbath.sh
# Generated: $(date)

gauge_group         ${GAUGE_GROUP^^}
beta                ${BETA}
seed                ${SEED}
T                   ${T}
L                   ${L}
boundary            ${BOUNDARY}
save_interval       ${SAVE_INTERVAL}
run_dir             ${RUN_DIR}

# Files
plaquette_file      ${RUN_OUTPUT_DIR}/${PLAQ_FILE}
config_dir          ${RUN_CONFIG_DIR}
EOF

echo ""
echo "=============================================="
echo "Done! To measure topological charge, run:"
echo "  ./build/bin/meas_topcharge${GAUGE_GROUP/su2/}$([ "$GAUGE_GROUP" = "su3" ] && echo "_su3") -i ${INPUT_FILE}"
echo "Then re-run analysis:"
echo "  python3 analysis/analysis.py --${GAUGE_GROUP} --results-dir ${RESULTS_DIR}"
echo "=============================================="
