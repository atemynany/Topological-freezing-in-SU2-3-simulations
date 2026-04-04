#!/bin/bash
# ==============================================================================
# run_heatbath_scan.sh  —  Phase 1: generate gauge configs for all beta values
# ==============================================================================
# For each beta: auto-seeds, runs heatbath, detects thermalization, saves
# input.txt with start_conf so run_topcharge_scan.sh can pick up from here.
#
# Usage:
#   ./run_heatbath_scan.sh [--su2|--su3] [options]
#
# Options:
#   --su2           SU(2) gauge group (default)
#   --su3           SU(3) gauge group
#   --open          Open boundary conditions
#   --dry-run       Preview without executing
#   --skip-build    Skip build step
#   --beta-file     Custom beta scan file
#   --params-file   Custom base parameters file
#
# After this completes, run:
#   ./run_topcharge_scan.sh [--su2|--su3] [same options]
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"
BUILD_DIR="$PROJECT_DIR/build"
RESULTS_DIR="$PROJECT_DIR/data/results"
CONDA_ENV="master_thesis"

GAUGE_GROUP="su2"
OPEN_BC=false
DRY_RUN=false
SKIP_BUILD=false
BETA_FILE=""
LATTICE_FILE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --su2)           GAUGE_GROUP="su2"; shift ;;
        --su3)           GAUGE_GROUP="su3"; shift ;;
        --open)          OPEN_BC=true;      shift ;;
        --dry-run)       DRY_RUN=true;      shift ;;
        --skip-build)    SKIP_BUILD=true;   shift ;;
        --beta-file)     BETA_FILE="$2";    shift 2 ;;
        --lattice-file)  LATTICE_FILE="$2"; shift 2 ;;
        -h|--help)       head -25 "$0" | tail -22; exit 0 ;;
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

if [ -z "$BETA_FILE" ]; then
    if [ "$OPEN_BC" = true ]; then
        BETA_FILE="input/beta_scan_${GAUGE_GROUP}_open.txt"
        [ ! -f "$BETA_FILE" ] && BETA_FILE="input/beta_scan_${GAUGE_GROUP}.txt"
    else
        BETA_FILE="input/beta_scan_${GAUGE_GROUP}.txt"
    fi
    [ ! -f "$BETA_FILE" ] && BETA_FILE="input/example_beta_scan.txt"
fi

PARAMS_FILE="$LATTICE_FILE"

if [ "$GAUGE_GROUP" = "su3" ]; then
    HEATBATH_BIN="mc_heatbath_su3"
    PLAQ_FILE="plaquette_su3.dat"
else
    HEATBATH_BIN="mc_heatbath"
    PLAQ_FILE="plaquette.dat"
fi

GAUGE_GROUP_UPPER=$(echo "$GAUGE_GROUP" | tr '[:lower:]' '[:upper:]')

# ==============================================================================
# Helpers (copied from run_beta_scan.sh)
# ==============================================================================
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; MAGENTA='\033[0;35m'; NC='\033[0m'

log_info()    { echo -e "${CYAN}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_warn()    { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error()   { echo -e "${RED}[ERROR]${NC} $1"; }
log_step()    { echo -e "${MAGENTA}[STEP]${NC} $1"; }

read_param() { grep -E "^${1}\s+" "$2" | awk '{print $2}' | head -1; }

generate_seed() {
    python3 -c "import time; print(int(time.time()*1e9) % 2147483646 + 1)" 2>/dev/null \
        || echo $(($(date +%s) * $$ % 2147483647))
}

create_input_file() {
    local output_file="$1" beta="$2" seed="$3" start_conf="$4"
    local run_output_dir="$5" run_config_dir="$6"

    local T=$(read_param "T" "$PARAMS_FILE")
    local L=$(read_param "L" "$PARAMS_FILE")
    local start_type=$(read_param "start_type" "$PARAMS_FILE")
    local boundary=$(read_param "boundary" "$PARAMS_FILE")
    local num_sweeps=$(read_param "num_sweeps" "$PARAMS_FILE")
    local save_interval=$(read_param "save_interval" "$PARAMS_FILE")
    local end_conf=$(read_param "end_conf" "$PARAMS_FILE")
    local conf_step=$(read_param "conf_step" "$PARAMS_FILE")
    local smear_steps=$(read_param "smear_steps" "$PARAMS_FILE")
    local smear_alpha=$(read_param "smear_alpha" "$PARAMS_FILE")
    local exclude_boundary_slices=$(read_param "exclude_boundary_slices" "$PARAMS_FILE")
    local overrelax_steps=$(read_param "overrelax_steps" "$PARAMS_FILE")
    local smear_interval=$(read_param "smear_interval" "$PARAMS_FILE")

    cat > "$output_file" << EOF
# Auto-generated input file
# Beta scan: beta=$beta, seed=$seed
# Generated: $(date)

T                   ${T:-16}
L                   ${L:-16}
beta                $beta
seed                $seed
start_type          ${start_type:-cold}
boundary            ${boundary:-periodic}
num_sweeps          ${num_sweeps:-1000}
save_interval       ${save_interval:-10}
EOF

    if [ "$GAUGE_GROUP" = "su3" ] && [ -n "$overrelax_steps" ]; then
        echo "overrelax_steps     $overrelax_steps" >> "$output_file"
    fi

    cat >> "$output_file" << EOF

start_conf          ${start_conf:-50}
end_conf            ${end_conf:-1000}
conf_step           ${conf_step:-10}
smear_steps         ${smear_steps:-20}
smear_alpha         ${smear_alpha:-0.3}
exclude_boundary_slices ${exclude_boundary_slices:-0}
EOF

    if [ "$GAUGE_GROUP" = "su2" ] && [ -n "$smear_interval" ]; then
        echo "smear_interval      $smear_interval" >> "$output_file"
    fi

    cat >> "$output_file" << EOF

output_dir          ${run_output_dir}/
config_dir          ${run_config_dir}/
EOF

    if [ "$GAUGE_GROUP" = "su3" ]; then
        echo "output_file         ${run_output_dir}/topcharge_su3.dat" >> "$output_file"
    fi
}

# ==============================================================================
# Validation
# ==============================================================================
echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}     ${GAUGE_GROUP_UPPER} Heatbath Scan — Phase 1/2${NC}"
echo -e "${GREEN}======================================================${NC}"
echo ""

[ ! -f "$BETA_FILE" ]    && { log_error "Beta file not found: $BETA_FILE";    exit 1; }
[ ! -f "$LATTICE_FILE" ] && { log_error "Lattice file not found: $LATTICE_FILE"; exit 1; }

SAVE_INTERVAL=$(read_param "save_interval" "$LATTICE_FILE")
T_SIZE=$(read_param "T" "$LATTICE_FILE")
L_SIZE=$(read_param "L" "$LATTICE_FILE")

log_info "Gauge group:   ${GAUGE_GROUP_UPPER}"
log_info "Beta file:     $BETA_FILE"
log_info "Lattice file:  $LATTICE_FILE"
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
    log_step "Building $HEATBATH_BIN..."
    mkdir -p "$BUILD_DIR" && cd "$BUILD_DIR"
    cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
    make -j$(nproc 2>/dev/null || echo 4) "$HEATBATH_BIN" 2>&1 | tail -3
    cd "$PROJECT_DIR"
    log_success "Build complete"
    echo ""
fi

mkdir -p "$RESULTS_DIR"

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

    SEED=$(generate_seed)
    sleep 0.1
    BOUNDARY=$(read_param "boundary" "$PARAMS_FILE"); BOUNDARY=${BOUNDARY:-periodic}

    RUN_NAME="T${T_SIZE}_L${L_SIZE}_b${BETA}_${BOUNDARY}_seed${SEED}"
    RUN_DIR="$RESULTS_DIR/$RUN_NAME"
    RUN_OUTPUT_DIR="$RUN_DIR/output"
    RUN_CONFIG_DIR="$RUN_DIR/configs"

    log_info "Run name:  $RUN_NAME"
    log_info "Seed:      $SEED"

    if [ "$DRY_RUN" = true ]; then
        echo "  Would create: $RUN_DIR"
        echo "  Would run heatbath with beta=$BETA, seed=$SEED"
        echo "  Would detect thermalization"
        continue
    fi

    mkdir -p "$RUN_OUTPUT_DIR" "$RUN_CONFIG_DIR"
    TEMP_INPUT="$RUN_DIR/input.txt"
    create_input_file "$TEMP_INPUT" "$BETA" "$SEED" "10" "$RUN_OUTPUT_DIR" "$RUN_CONFIG_DIR"

    # -------------------------------------------------------------------------
    # Step 1: Heatbath
    # -------------------------------------------------------------------------
    log_step "Running MC heatbath..."
    cd "$PROJECT_DIR"
    "$BUILD_DIR/bin/$HEATBATH_BIN" -i "$TEMP_INPUT" 2>&1 | tee "$RUN_DIR/heatbath.log"
    log_success "Heatbath complete"

    # -------------------------------------------------------------------------
    # Step 2: Detect thermalization
    # -------------------------------------------------------------------------
    log_step "Detecting thermalization point..."
    PLAQ_PATH="$RUN_OUTPUT_DIR/$PLAQ_FILE"
    if [ ! -f "$PLAQ_PATH" ] && [ -f "$PROJECT_DIR/output/$PLAQ_FILE" ]; then
        mv "$PROJECT_DIR/output/$PLAQ_FILE" "$PLAQ_PATH"
    fi
    if [ ! -f "$PLAQ_PATH" ]; then
        log_warn "Plaquette file not found, using default start_conf=50"
        START_CONF=50
    else
        START_CONF=$(conda run -n "$CONDA_ENV" python3 \
            "$SCRIPT_DIR/scripts/detect_thermalization.py" "$PLAQ_PATH" "$SAVE_INTERVAL" 2>/dev/null \
            || echo "50")
        [[ "$START_CONF" =~ ^[0-9]+$ ]] || START_CONF=50
    fi
    log_info "start_conf = $START_CONF"

    # -------------------------------------------------------------------------
    # Step 3: Update input.txt with start_conf (ready for topcharge scan)
    # -------------------------------------------------------------------------
    create_input_file "$TEMP_INPUT" "$BETA" "$SEED" "$START_CONF" "$RUN_OUTPUT_DIR" "$RUN_CONFIG_DIR"

    cat > "$RUN_DIR/run_info.txt" << EOF
# Run Information — heatbath phase
# Generated: $(date)

gauge_group         ${GAUGE_GROUP_UPPER}
beta                $BETA
seed                $SEED
T                   $T_SIZE
L                   $L_SIZE
boundary            $BOUNDARY
save_interval       $SAVE_INTERVAL
start_conf          $START_CONF
run_dir             $RUN_DIR

# Files
plaquette_file      $RUN_OUTPUT_DIR/$PLAQ_FILE
config_dir          $RUN_CONFIG_DIR
EOF

    log_success "Beta=$BETA done — configs in $RUN_CONFIG_DIR"
done

echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}  Heatbath scan complete. To measure topcharge run:${NC}"
echo -e "${GREEN}  ./run_topcharge_scan.sh --${GAUGE_GROUP}${NC}"
echo -e "${GREEN}======================================================${NC}"
