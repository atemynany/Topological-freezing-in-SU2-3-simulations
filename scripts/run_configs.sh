#!/bin/bash
# ==============================================================================
# run_configs.sh — Config Generation (Step 1)
# ==============================================================================
# For each beta value: runs MC heatbath, detects thermalization, generates
# thermalization plots. Does NOT measure topological charge.
#
# Usage:
#   ./scripts/run_configs.sh [--su2|--su3] [options]
#
# Options:
#   --su2           SU(2) gauge group (default)
#   --su3           SU(3) gauge group
#   --dry-run       Preview without executing
#   --skip-build    Skip build step
#   --beta-file     Custom beta scan file
#   --params-file   Custom base parameters file
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

GAUGE_GROUP="su2"
DRY_RUN=false
SKIP_BUILD=false
BETA_FILE=""
PARAMS_FILE=""
CONDA_ENV="master_thesis"

while [[ $# -gt 0 ]]; do
    case $1 in
        --su2)        GAUGE_GROUP="su2";  shift ;;
        --su3)        GAUGE_GROUP="su3";  shift ;;
        --dry-run)    DRY_RUN=true;       shift ;;
        --skip-build) SKIP_BUILD=true;    shift ;;
        --beta-file)  BETA_FILE="$2";     shift 2 ;;
        --params-file) PARAMS_FILE="$2";  shift 2 ;;
        -h|--help)
            echo "Usage: $0 [--su2|--su3] [options]"
            echo "  --beta-file     Custom beta values file"
            echo "  --params-file   Custom base parameters file"
            echo "  --dry-run       Preview without executing"
            echo "  --skip-build    Skip build step"
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [ "$GAUGE_GROUP" = "su2" ]; then
    BETA_FILE="${BETA_FILE:-$PROJECT_DIR/input/beta_scan_su2.txt}"
    PARAMS_FILE="${PARAMS_FILE:-$PROJECT_DIR/input/base_params_su2.txt}"
    HEATBATH_BIN="mc_heatbath"
    PLAQ_FILE="plaquette.dat"
    SU3_FLAG=""
else
    BETA_FILE="${BETA_FILE:-$PROJECT_DIR/input/beta_scan_su3.txt}"
    PARAMS_FILE="${PARAMS_FILE:-$PROJECT_DIR/input/base_params_su3.txt}"
    HEATBATH_BIN="mc_heatbath_su3"
    PLAQ_FILE="plaquette_su3.dat"
    SU3_FLAG="--su3"
fi

RESULTS_DIR="$PROJECT_DIR/data/results"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m'

log_info()    { echo -e "${CYAN}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[OK]${NC} $1"; }
log_warn()    { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error()   { echo -e "${RED}[ERROR]${NC} $1"; }
log_step()    { echo -e "${MAGENTA}[STEP]${NC} $1"; }

has_working_conda_env() {
    command -v conda &> /dev/null && conda run -n "$CONDA_ENV" python3 --version &> /dev/null
}

generate_seed() {
    python3 -c "import time; print(int(time.time()*1e9))" 2>/dev/null | cut -c1-10 \
        || echo $(($(date +%s) * $$ % 2147483647))
}

read_param() {
    grep -E "^${1}\s+" "$2" | awk '{print $2}' | head -1
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
GAUGE_GROUP_UPPER=$(echo "$GAUGE_GROUP" | tr '[:lower:]' '[:upper:]')

echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}     ${GAUGE_GROUP_UPPER} Config Generation${NC}"
echo -e "${GREEN}======================================================${NC}"
echo ""

if [ ! -f "$BETA_FILE" ]; then
    log_error "Beta file not found: $BETA_FILE"
    exit 1
fi
if [ ! -f "$PARAMS_FILE" ]; then
    log_error "Params file not found: $PARAMS_FILE"
    exit 1
fi

SAVE_INTERVAL=$(read_param "save_interval" "$PARAMS_FILE")
T_SIZE=$(read_param "T" "$PARAMS_FILE")
L_SIZE=$(read_param "L" "$PARAMS_FILE")

log_info "Gauge group:     ${GAUGE_GROUP_UPPER}"
log_info "Beta file:       $BETA_FILE"
log_info "Params file:     $PARAMS_FILE"
log_info "Lattice:         ${T_SIZE}x${L_SIZE}^3"
log_info "Save interval:   $SAVE_INTERVAL"
log_info "Results dir:     $RESULTS_DIR"
echo ""

BETAS=($(grep -v '^#' "$BETA_FILE" | grep -v '^$'))
log_info "Beta values: ${BETAS[*]}"
log_info "Total runs: ${#BETAS[@]}"
echo ""

if [ "$DRY_RUN" = true ]; then
    log_warn "DRY RUN MODE"
    echo ""
fi

# ==============================================================================
# Build
# ==============================================================================
if [ "$SKIP_BUILD" = false ]; then
    log_step "Building $HEATBATH_BIN..."
    if [ "$DRY_RUN" = false ]; then
        mkdir -p "$PROJECT_DIR/build"
        cd "$PROJECT_DIR/build"
        cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
        make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4) "$HEATBATH_BIN" 2>&1 | tail -3
        cd "$PROJECT_DIR"
        log_success "Build complete"
    else
        echo "  Would run: cmake && make $HEATBATH_BIN"
    fi
    echo ""
fi

# ==============================================================================
# Main Loop
# ==============================================================================
mkdir -p "$RESULTS_DIR"

RUN_COUNT=0
TOTAL_RUNS=${#BETAS[@]}

for BETA in "${BETAS[@]}"; do
    RUN_COUNT=$((RUN_COUNT + 1))

    echo ""
    echo -e "${GREEN}------------------------------------------------------${NC}"
    echo -e "${GREEN}  Run $RUN_COUNT/$TOTAL_RUNS: beta = $BETA${NC}"
    echo -e "${GREEN}------------------------------------------------------${NC}"

    SEED=$(generate_seed)
    sleep 0.1

    BOUNDARY=$(read_param "boundary" "$PARAMS_FILE")
    BOUNDARY=${BOUNDARY:-periodic}

    RUN_NAME="T${T_SIZE}_L${L_SIZE}_b${BETA}_${BOUNDARY}_seed${SEED}"
    RUN_DIR="$RESULTS_DIR/$RUN_NAME"
    RUN_OUTPUT_DIR="$RUN_DIR/output"
    RUN_CONFIG_DIR="$RUN_DIR/configs"

    log_info "Run name: $RUN_NAME"
    log_info "Seed:     $SEED"

    if [ "$DRY_RUN" = true ]; then
        echo "  Would create: $RUN_DIR"
        echo "  Would run heatbath with beta=$BETA, seed=$SEED"
        echo "  Would detect thermalization"
        continue
    fi

    mkdir -p "$RUN_OUTPUT_DIR"
    mkdir -p "$RUN_CONFIG_DIR"

    TEMP_INPUT="$RUN_DIR/input.txt"
    create_input_file "$TEMP_INPUT" "$BETA" "$SEED" "10" "$RUN_OUTPUT_DIR" "$RUN_CONFIG_DIR"

    # Step 1: Heatbath
    log_step "Running heatbath..."
    cd "$PROJECT_DIR"
    ./build/bin/$HEATBATH_BIN -i "$TEMP_INPUT" 2>&1 | tee "$RUN_DIR/heatbath.log"
    log_success "Heatbath complete"

    # Step 2: Detect thermalization
    log_step "Detecting thermalization..."
    PLAQ_PATH="$RUN_OUTPUT_DIR/$PLAQ_FILE"
    if [ ! -f "$PLAQ_PATH" ]; then
        log_error "Plaquette file not found: $PLAQ_PATH"
        continue
    fi

    START_CONF=""
    if has_working_conda_env; then
        START_CONF=$(conda run -n "$CONDA_ENV" python3 "$SCRIPT_DIR/detect_thermalization.py" "$PLAQ_PATH" "$SAVE_INTERVAL" 2>"$RUN_DIR/thermalization.err") || START_CONF=""
    fi
    if [ -z "$START_CONF" ] && command -v python3 &> /dev/null; then
        START_CONF=$(python3 "$SCRIPT_DIR/detect_thermalization.py" "$PLAQ_PATH" "$SAVE_INTERVAL" 2>"$RUN_DIR/thermalization.err") || START_CONF=""
    fi
    if ! [[ "$START_CONF" =~ ^[0-9]+$ ]]; then
        log_warn "Thermalization detection failed, using default start_conf=50"
        START_CONF=50
    fi

    log_info "start_conf=$START_CONF"

    # Update input file with detected start_conf
    create_input_file "$TEMP_INPUT" "$BETA" "$SEED" "$START_CONF" "$RUN_OUTPUT_DIR" "$RUN_CONFIG_DIR"

    # Step 3: Thermalization plots
    log_step "Generating thermalization plots..."
    if has_working_conda_env; then
        conda run -n "$CONDA_ENV" python3 "$SCRIPT_DIR/plot_therm.py" "$RUN_DIR" $SU3_FLAG 2>/dev/null \
            && log_success "Thermalization plots saved" \
            || log_warn "Thermalization plot failed (non-fatal)"
    elif command -v python3 &> /dev/null; then
        python3 "$SCRIPT_DIR/plot_therm.py" "$RUN_DIR" $SU3_FLAG 2>/dev/null \
            && log_success "Thermalization plots saved" \
            || log_warn "Thermalization plot failed (non-fatal)"
    fi

    # Save run info
    cat > "$RUN_DIR/run_info.txt" << EOF
# Run Information — run_configs.sh
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

plaquette_file      $RUN_OUTPUT_DIR/$PLAQ_FILE
config_dir          $RUN_CONFIG_DIR
EOF

    log_success "Run $RUN_COUNT/$TOTAL_RUNS complete: $RUN_NAME"
done

# ==============================================================================
# Summary
# ==============================================================================
echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}     Config Generation Complete!${NC}"
echo -e "${GREEN}======================================================${NC}"
echo ""
log_info "Runs completed: $RUN_COUNT"
log_info "Results: $RESULTS_DIR"
echo ""
echo "Run directories:"
for BETA in "${BETAS[@]}"; do
    ls -d "$RESULTS_DIR"/T${T_SIZE}_L${L_SIZE}_b${BETA}_* 2>/dev/null || true
done
echo ""
log_info "Next step: measure topological charge with:"
echo "  ./scripts/run_topcharge_scan.sh --${GAUGE_GROUP}"
