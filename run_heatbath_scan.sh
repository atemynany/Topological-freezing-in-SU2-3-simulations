#!/bin/bash
# ==============================================================================
# run_heatbath_scan.sh  —  Phase 1: generate gauge configs for all beta values
# ==============================================================================
# Discovers all input/heatbath_input/setup_*_${GAUGE_GROUP}.txt files and runs heatbath for
# every (setup_file, beta) combination. Each setup file contains lattice params
# and one or more "beta X.X" lines.
#
# Usage:
#   ./run_heatbath_scan.sh [--su2|--su3] [options]
#
# Options:
#   --su2           SU(2) gauge group (default)
#   --su3           SU(3) gauge group
#   --dry-run       Preview without executing
#   --skip-build    Skip build step
#   --parallel      Run all beta values in parallel (uses all CPUs)
#   --jobs N        Run up to N beta values in parallel
#
# After this completes, run:
#   ./run_topcharge_scan.sh [--su2|--su3] [same options]
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"
BUILD_DIR="$PROJECT_DIR/build"
RESULTS_DIR="${HEATBATH_RESULTS_DIR:-$PROJECT_DIR/data/results}"
CONDA_ENV="master_thesis"

GAUGE_GROUP="su2"
DRY_RUN=false
SKIP_BUILD=false
PARALLEL_JOBS=0   # 0 = sequential (default)

while [[ $# -gt 0 ]]; do
    case $1 in
        --su2)           GAUGE_GROUP="su2"; shift ;;
        --su3)           GAUGE_GROUP="su3"; shift ;;
        --dry-run)       DRY_RUN=true;      shift ;;
        --skip-build)    SKIP_BUILD=true;   shift ;;
        --parallel)      PARALLEL_JOBS=$(nproc 2>/dev/null || echo 4); shift ;;
        --jobs)          PARALLEL_JOBS="$2"; shift 2 ;;
        -h|--help)       head -22 "$0" | tail -19; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [ "$GAUGE_GROUP" = "su3" ]; then
    HEATBATH_BIN="mc_heatbath_su3"
    PLAQ_FILE="plaquette_su3.dat"
    ACTION_FILE="action_density_su3.dat"
else
    HEATBATH_BIN="mc_heatbath"
    PLAQ_FILE="plaquette.dat"
    ACTION_FILE="action_density.dat"
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

generate_seed() {
    python3 -c "import time; print(int(time.time()*1e9) % 2147483646 + 1)" 2>/dev/null \
        || echo $(($(date +%s) * $$ % 2147483647))
}

create_input_file() {
    local output_file="$1" beta="$2" seed="$3" start_conf="$4"
    local run_output_dir="$5" run_config_dir="$6" setup_file="$7"

    local T=$(read_param "T" "$setup_file")
    local L=$(read_param "L" "$setup_file")
    local start_type=$(read_param "start_type" "$setup_file")
    local boundary=$(read_param "boundary" "$setup_file")
    local num_sweeps=$(read_param "num_sweeps" "$setup_file")
    local save_interval=$(read_param "save_interval" "$setup_file")
    local overrelax_steps=$(read_param "overrelax_steps" "$setup_file")

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
start_conf          ${start_conf:-50}
EOF

    if [ "$GAUGE_GROUP" = "su3" ] && [ -n "$overrelax_steps" ]; then
        echo "overrelax_steps     $overrelax_steps" >> "$output_file"
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
# Discover setup files
# ==============================================================================
SETUP_FILES=("$PROJECT_DIR"/input/heatbath_input/setup_*_${GAUGE_GROUP}.txt)

if [ ! -f "${SETUP_FILES[0]}" ]; then
    log_error "No setup files found matching input/heatbath_input/setup_*_${GAUGE_GROUP}.txt"
    exit 1
fi

# ==============================================================================
# Build flat job list: (setup_file, beta) pairs
# ==============================================================================
declare -a JOB_SETUP=()
declare -a JOB_BETA=()

for SETUP in "${SETUP_FILES[@]}"; do
    BETAS_IN_FILE=($(grep -E '^beta\s+' "$SETUP" | awk '{print $2}'))
    if [ ${#BETAS_IN_FILE[@]} -eq 0 ]; then
        log_warn "No beta values in $SETUP, skipping"
        continue
    fi
    for B in "${BETAS_IN_FILE[@]}"; do
        JOB_SETUP+=("$SETUP")
        JOB_BETA+=("$B")
    done
done

TOTAL=${#JOB_BETA[@]}
if [ "$TOTAL" -eq 0 ]; then
    log_error "No jobs found across setup files"
    exit 1
fi

# ==============================================================================
# Validation
# ==============================================================================
echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}     ${GAUGE_GROUP_UPPER} Heatbath Scan — Phase 1/2${NC}"
echo -e "${GREEN}======================================================${NC}"
echo ""

log_info "Gauge group:   ${GAUGE_GROUP_UPPER}"
log_info "Setup files:   ${#SETUP_FILES[@]}"
for SF in "${SETUP_FILES[@]}"; do
    log_info "  $(basename "$SF")"
done
log_info "Total jobs:    $TOTAL"
log_info "Results dir:   $RESULTS_DIR"
echo ""

# ==============================================================================
# Build (skipped in dry-run)
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
# Pre-generate seeds
# ==============================================================================
declare -a SEEDS=()
for i in "${!JOB_BETA[@]}"; do
    SEEDS+=("$(generate_seed)")
    sleep 0.1
done

# ==============================================================================
# Single-beta worker
# ==============================================================================
run_single_beta() {
    local SETUP_FILE="$1"
    local BETA="$2"
    local SEED="$3"
    local RUN_COUNT="$4"
    local TOTAL="$5"
    local TO_LOG="$6"

    local T_SIZE=$(read_param "T" "$SETUP_FILE")
    local L_SIZE=$(read_param "L" "$SETUP_FILE")
    local SAVE_INTERVAL=$(read_param "save_interval" "$SETUP_FILE")
    local BOUNDARY
    BOUNDARY=$(read_param "boundary" "$SETUP_FILE"); BOUNDARY=${BOUNDARY:-periodic}

    local RUN_NAME="T${T_SIZE}_L${L_SIZE}_b${BETA}_${BOUNDARY}_seed${SEED}_${GAUGE_GROUP}"
    local RUN_DIR="$RESULTS_DIR/$RUN_NAME"
    local RUN_OUTPUT_DIR="$RUN_DIR/output"
    local RUN_CONFIG_DIR="$RUN_DIR/configs"

    echo -e "${GREEN}------------------------------------------------------${NC}"
    echo -e "${GREEN}  [$RUN_COUNT/$TOTAL] beta = $BETA  (${T_SIZE}x${L_SIZE}^3, ${BOUNDARY}, seed=$SEED)${NC}"
    echo -e "${GREEN}------------------------------------------------------${NC}"

    mkdir -p "$RUN_OUTPUT_DIR" "$RUN_CONFIG_DIR"
    local TEMP_INPUT="$RUN_DIR/input.txt"
    create_input_file "$TEMP_INPUT" "$BETA" "$SEED" "10" "$RUN_OUTPUT_DIR" "$RUN_CONFIG_DIR" "$SETUP_FILE"

    # Heatbath
    log_step "[$RUN_COUNT/$TOTAL] Running heatbath for beta=$BETA..."
    cd "$PROJECT_DIR"
    if [ "$TO_LOG" = "tee" ]; then
        "$BUILD_DIR/bin/$HEATBATH_BIN" -i "$TEMP_INPUT" 2>&1 | tee "$RUN_DIR/heatbath.log"
    else
        "$BUILD_DIR/bin/$HEATBATH_BIN" -i "$TEMP_INPUT" > "$RUN_DIR/heatbath.log" 2>&1
    fi
    local heatbath_exit=${PIPESTATUS[0]}
    if [ "$heatbath_exit" -ne 0 ]; then
        log_error "[$RUN_COUNT/$TOTAL] Heatbath binary exited with code $heatbath_exit for beta=$BETA"
        return 1
    fi

    # Detect thermalization
    local PLAQ_PATH="$RUN_OUTPUT_DIR/$PLAQ_FILE"
    if [ ! -f "$PLAQ_PATH" ] && [ -f "$PROJECT_DIR/output/$PLAQ_FILE" ]; then
        mv "$PROJECT_DIR/output/$PLAQ_FILE" "$PLAQ_PATH"
    fi

    local START_CONF=50
    if [ -f "$PLAQ_PATH" ]; then
        START_CONF=$(conda run -n "$CONDA_ENV" python3 \
            "$SCRIPT_DIR/scripts/detect_thermalization.py" "$PLAQ_PATH" "${SAVE_INTERVAL:-10}" 2>/dev/null \
            || echo "50")
        [[ "$START_CONF" =~ ^[0-9]+$ ]] || START_CONF=50
    else
        log_warn "[$RUN_COUNT/$TOTAL] Plaquette file not found, using start_conf=50"
    fi

    # Update input.txt with start_conf (ready for topcharge scan)
    create_input_file "$TEMP_INPUT" "$BETA" "$SEED" "$START_CONF" "$RUN_OUTPUT_DIR" "$RUN_CONFIG_DIR" "$SETUP_FILE"

    cat > "$RUN_DIR/run_info.txt" << EOF
# Run Information — heatbath phase
# Generated: $(date)

gauge_group         ${GAUGE_GROUP_UPPER}
beta                $BETA
seed                $SEED
T                   $T_SIZE
L                   $L_SIZE
boundary            $BOUNDARY
save_interval       ${SAVE_INTERVAL:-10}
start_conf          $START_CONF
run_dir             $RUN_DIR

# Files
plaquette_file      $RUN_OUTPUT_DIR/$PLAQ_FILE
action_density_file $RUN_OUTPUT_DIR/$ACTION_FILE
config_dir          $RUN_CONFIG_DIR
EOF

    log_success "[$RUN_COUNT/$TOTAL] Beta=$BETA done (start_conf=$START_CONF) — configs in $RUN_CONFIG_DIR"
}

# ==============================================================================
# Dry-run preview
# ==============================================================================
if [ "$DRY_RUN" = true ]; then
    for i in "${!JOB_BETA[@]}"; do
        SETUP="${JOB_SETUP[$i]}"; BETA="${JOB_BETA[$i]}"; SEED="${SEEDS[$i]}"
        T_SIZE=$(read_param "T" "$SETUP"); L_SIZE=$(read_param "L" "$SETUP")
        BOUNDARY=$(read_param "boundary" "$SETUP"); BOUNDARY=${BOUNDARY:-periodic}
        RUN_NAME="T${T_SIZE}_L${L_SIZE}_b${BETA}_${BOUNDARY}_seed${SEED}_${GAUGE_GROUP}"
        echo "  Would create: $RESULTS_DIR/$RUN_NAME"
        echo "    Setup: $(basename "$SETUP"), beta=$BETA, seed=$SEED"
    done
    echo ""
    echo "Mode: $( [ "$PARALLEL_JOBS" -gt 0 ] && echo "parallel (max $PARALLEL_JOBS jobs)" || echo "sequential" )"
    exit 0
fi

# ==============================================================================
# Main loop
# ==============================================================================
if [ "$PARALLEL_JOBS" -gt 0 ]; then
    log_info "Parallel mode: running up to $PARALLEL_JOBS job(s) simultaneously"
    log_info "Progress visible in per-run heatbath.log files"
    echo ""

    declare -a PIDS=()
    declare -a PID_LABELS=()
    FAILED_COUNT=0

    for i in "${!JOB_BETA[@]}"; do
        SETUP="${JOB_SETUP[$i]}"; BETA="${JOB_BETA[$i]}"; SEED="${SEEDS[$i]}"
        COUNT=$((i+1))

        # Throttle: wait until a slot is free
        while [ "${#PIDS[@]}" -ge "$PARALLEL_JOBS" ]; do
            NEW_PIDS=(); NEW_LABELS=()
            for j in "${!PIDS[@]}"; do
                p="${PIDS[$j]}"
                if kill -0 "$p" 2>/dev/null; then
                    NEW_PIDS+=("$p"); NEW_LABELS+=("${PID_LABELS[$j]}")
                else
                    if wait "$p"; then
                        :
                    else
                        log_error "Job for ${PID_LABELS[$j]} failed — check heatbath.log"
                        FAILED_COUNT=$((FAILED_COUNT+1))
                    fi
                fi
            done
            PIDS=("${NEW_PIDS[@]}"); PID_LABELS=("${NEW_LABELS[@]}")
            [ "${#PIDS[@]}" -ge "$PARALLEL_JOBS" ] && sleep 0.5
        done

        log_info "Launching beta=$BETA [job $COUNT/$TOTAL]..."
        run_single_beta "$SETUP" "$BETA" "$SEED" "$COUNT" "$TOTAL" "log" &
        PIDS+=($!); PID_LABELS+=("beta=$BETA")
    done

    # Wait for remaining jobs
    for j in "${!PIDS[@]}"; do
        p="${PIDS[$j]}"
        if wait "$p"; then
            :
        else
            log_error "Job for ${PID_LABELS[$j]} failed — check heatbath.log"
            FAILED_COUNT=$((FAILED_COUNT+1))
        fi
    done

    echo ""
    if [ "$FAILED_COUNT" -gt 0 ]; then
        log_warn "$FAILED_COUNT job(s) failed. Check individual heatbath.log files."
    fi

else
    for i in "${!JOB_BETA[@]}"; do
        SETUP="${JOB_SETUP[$i]}"; BETA="${JOB_BETA[$i]}"; SEED="${SEEDS[$i]}"
        COUNT=$((i+1))
        echo ""
        run_single_beta "$SETUP" "$BETA" "$SEED" "$COUNT" "$TOTAL" "tee"
    done
fi

echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}  Heatbath scan complete. To measure topcharge run:${NC}"
echo -e "${GREEN}  ./run_topcharge_scan.sh --${GAUGE_GROUP}${NC}"
echo -e "${GREEN}======================================================${NC}"
