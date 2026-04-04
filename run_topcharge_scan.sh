#!/bin/bash
# ==============================================================================
# run_topcharge_scan.sh  —  Phase 2: measure topological charge on existing configs
# ==============================================================================
# Finds run dirs created by run_heatbath_scan.sh and runs the topcharge
# measurement on each. Uses the input.txt (with start_conf) already saved
# by the heatbath phase.
#
# Usage:
#   ./run_topcharge_scan.sh [--su2|--su3] [options]
#
# Options:
#   --su2           SU(2) gauge group (default)
#   --su3           SU(3) gauge group
#   --open          Open boundary conditions (used for dir matching)
#   --dry-run       Preview without executing
#   --skip-build    Skip build step
#   --beta-file     Custom beta scan file
#   --params-file   Custom base parameters file
#
# Requires run dirs from:
#   ./run_heatbath_scan.sh [--su2|--su3] [same options]
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"
BUILD_DIR="$PROJECT_DIR/build"
RESULTS_DIR="$PROJECT_DIR/data/results"

GAUGE_GROUP="su2"
OPEN_BC=false
DRY_RUN=false
SKIP_BUILD=false
BETA_FILE=""
LATTICE_FILE=""
TOPCHARGE_FILE_PARAMS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --su2)              GAUGE_GROUP="su2";          shift ;;
        --su3)              GAUGE_GROUP="su3";          shift ;;
        --open)             OPEN_BC=true;               shift ;;
        --dry-run)          DRY_RUN=true;               shift ;;
        --skip-build)       SKIP_BUILD=true;            shift ;;
        --beta-file)        BETA_FILE="$2";             shift 2 ;;
        --lattice-file)     LATTICE_FILE="$2";          shift 2 ;;
        --topcharge-file)   TOPCHARGE_FILE_PARAMS="$2"; shift 2 ;;
        -h|--help)          head -25 "$0" | tail -22;   exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ==============================================================================
# File selection
# ==============================================================================
if [ -z "$LATTICE_FILE" ]; then
    LATTICE_FILE="input/lattice_${GAUGE_GROUP}.txt"
    [ ! -f "$LATTICE_FILE" ] && LATTICE_FILE="input/example_base_params.txt"
fi

if [ -z "$TOPCHARGE_FILE_PARAMS" ]; then
    TOPCHARGE_FILE_PARAMS="input/topcharge_params_${GAUGE_GROUP}.txt"
fi

if [ -z "$BETA_FILE" ]; then
    if [ "$OPEN_BC" = true ]; then
        BETA_FILE="input/beta_scan_${GAUGE_GROUP}_open.txt"
        [ ! -f "$BETA_FILE" ] && BETA_FILE="input/beta_scan_${GAUGE_GROUP}.txt"
    else
        BETA_FILE="input/beta_scan_${GAUGE_GROUP}.txt"
    fi
    [ ! -f "$BETA_FILE" ] && BETA_FILE="input/example_beta_scan.txt"
fi

if [ "$GAUGE_GROUP" = "su3" ]; then
    TOPCHARGE_BIN="meas_topcharge_su3"
    TOPCHARGE_FILE="topcharge_su3.dat"
else
    TOPCHARGE_BIN="meas_topcharge"
    TOPCHARGE_FILE="topcharge.dat"
fi

GAUGE_GROUP_UPPER=$(echo "$GAUGE_GROUP" | tr '[:lower:]' '[:upper:]')

# ==============================================================================
# Helpers
# ==============================================================================
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; MAGENTA='\033[0;35m'; NC='\033[0m'

log_info()    { echo -e "${CYAN}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_warn()    { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error()   { echo -e "${RED}[ERROR]${NC} $1"; }
log_step()    { echo -e "${MAGENTA}[STEP]${NC} $1"; }

read_param() { grep -E "^${1}\s+" "$2" | awk '{print $2}' | head -1; }

# ==============================================================================
# Validation
# ==============================================================================
echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}     ${GAUGE_GROUP_UPPER} Topcharge Scan — Phase 2/2${NC}"
echo -e "${GREEN}======================================================${NC}"
echo ""

[ ! -f "$BETA_FILE" ]            && { log_error "Beta file not found: $BETA_FILE";             exit 1; }
[ ! -f "$LATTICE_FILE" ]         && { log_error "Lattice file not found: $LATTICE_FILE";         exit 1; }
[ ! -f "$TOPCHARGE_FILE_PARAMS" ] && { log_error "Topcharge params not found: $TOPCHARGE_FILE_PARAMS"; exit 1; }

T_SIZE=$(read_param "T" "$LATTICE_FILE")
L_SIZE=$(read_param "L" "$LATTICE_FILE")
BOUNDARY=$(read_param "boundary" "$LATTICE_FILE"); BOUNDARY=${BOUNDARY:-periodic}

log_info "Gauge group:       ${GAUGE_GROUP_UPPER}"
log_info "Beta file:         $BETA_FILE"
log_info "Lattice file:      $LATTICE_FILE"
log_info "Topcharge params:  $TOPCHARGE_FILE_PARAMS"
log_info "Lattice:       ${T_SIZE}x${L_SIZE}^3"
log_info "Results dir:   $RESULTS_DIR"
echo ""

BETAS=($(grep -v '^#' "$BETA_FILE" | grep -v '^$'))
log_info "Beta values: ${BETAS[*]}"
echo ""

# ==============================================================================
# Build
# ==============================================================================
if [ "$SKIP_BUILD" = false ] && [ "$DRY_RUN" = false ]; then
    log_step "Building $TOPCHARGE_BIN..."
    mkdir -p "$BUILD_DIR" && cd "$BUILD_DIR"
    cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
    make -j$(nproc 2>/dev/null || echo 4) "$TOPCHARGE_BIN" 2>&1 | tail -3
    cd "$PROJECT_DIR"
    log_success "Build complete"
    echo ""
fi

# ==============================================================================
# Main loop
# ==============================================================================
RUN_COUNT=0
for BETA in "${BETAS[@]}"; do
    RUN_COUNT=$((RUN_COUNT + 1))
    echo ""
    echo -e "${GREEN}------------------------------------------------------${NC}"
    echo -e "${GREEN}  Run $RUN_COUNT/${#BETAS[@]}: beta = $BETA${NC}"
    echo -e "${GREEN}------------------------------------------------------${NC}"

    # Find the run dir created by the heatbath phase
    RUN_DIR=$(ls -d "$RESULTS_DIR"/T${T_SIZE}_L${L_SIZE}_b${BETA}_${BOUNDARY}_seed* 2>/dev/null | head -1)

    if [ -z "$RUN_DIR" ]; then
        log_error "No existing run dir found for beta=$BETA (T=${T_SIZE}, L=${L_SIZE}, boundary=${BOUNDARY})"
        log_error "Run ./run_heatbath_scan.sh --${GAUGE_GROUP} first"
        continue
    fi

    RUN_OUTPUT_DIR="$RUN_DIR/output"
    TEMP_INPUT="$RUN_DIR/input.txt"
    RUN_NAME=$(basename "$RUN_DIR")

    log_info "Found run: $RUN_NAME"

    if [ ! -f "$TEMP_INPUT" ]; then
        log_error "input.txt missing in $RUN_DIR — was heatbath phase completed?"
        continue
    fi

    if [ "$DRY_RUN" = true ]; then
        echo "  Would run topcharge on: $RUN_DIR"
        echo "  Input: $TEMP_INPUT"
        continue
    fi

    # -------------------------------------------------------------------------
    # Patch input.txt with topcharge params from topcharge_params file
    # -------------------------------------------------------------------------
    END_CONF=$(read_param "end_conf"                "$TOPCHARGE_FILE_PARAMS")
    CONF_STEP=$(read_param "conf_step"              "$TOPCHARGE_FILE_PARAMS")
    SMEAR_STEPS=$(read_param "smear_steps"          "$TOPCHARGE_FILE_PARAMS")
    SMEAR_ALPHA=$(read_param "smear_alpha"          "$TOPCHARGE_FILE_PARAMS")
    SMEAR_INTERVAL=$(read_param "smear_interval"    "$TOPCHARGE_FILE_PARAMS")
    EXCLUDE_BC=$(read_param "exclude_boundary_slices" "$TOPCHARGE_FILE_PARAMS")

    for key in end_conf conf_step smear_steps smear_alpha smear_interval exclude_boundary_slices; do
        sed -i "s/^${key}[[:space:]].*//" "$TEMP_INPUT"
    done

    cat >> "$TEMP_INPUT" << EOF

# Topcharge params (from $TOPCHARGE_FILE_PARAMS)
end_conf            ${END_CONF:-1000}
conf_step           ${CONF_STEP:-10}
smear_steps         ${SMEAR_STEPS:-20}
smear_alpha         ${SMEAR_ALPHA:-0.3}
exclude_boundary_slices ${EXCLUDE_BC:-0}
EOF
    [ -n "$SMEAR_INTERVAL" ] && echo "smear_interval      $SMEAR_INTERVAL" >> "$TEMP_INPUT"

    # -------------------------------------------------------------------------
    # Run topological charge measurement
    # -------------------------------------------------------------------------
    log_step "Measuring topological charge..."
    cd "$PROJECT_DIR"
    "$BUILD_DIR/bin/$TOPCHARGE_BIN" -i "$TEMP_INPUT" 2>&1 | tee "$RUN_DIR/topcharge.log"

    TOPCHARGE_PATH="$RUN_OUTPUT_DIR/$TOPCHARGE_FILE"
    if [ ! -f "$TOPCHARGE_PATH" ] && [ -f "$PROJECT_DIR/output/$TOPCHARGE_FILE" ]; then
        mv "$PROJECT_DIR/output/$TOPCHARGE_FILE" "$TOPCHARGE_PATH"
    fi

    # Update run_info.txt with topcharge output
    cat >> "$RUN_DIR/run_info.txt" << EOF

# Topcharge phase
topcharge_file      $RUN_OUTPUT_DIR/$TOPCHARGE_FILE
EOF

    log_success "Beta=$BETA done — results in $RUN_OUTPUT_DIR"
done

echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}  Topcharge scan complete. To run analysis:${NC}"
echo -e "${GREEN}  python3 analysis/analysis.py --${GAUGE_GROUP}${NC}"
echo -e "${GREEN}======================================================${NC}"
