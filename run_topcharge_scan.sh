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
#   --dry-run       Preview without executing
#   --skip-build    Skip build step
#   --topcharge-file FILE
#                   Custom topcharge parameter file
#   --parallel      Run all beta values in parallel (uses all CPUs)
#   --jobs N        Run up to N beta values in parallel
#   --exclude PAT   Skip run directories whose path matches regex PAT
#
# Requires run dirs from:
#   ./run_heatbath_scan.sh [--su2|--su3] [same options]
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"
BUILD_DIR="$PROJECT_DIR/build"
RESULTS_DIR="${HEATBATH_RESULTS_DIR:-$PROJECT_DIR/data/results}"

GAUGE_GROUP="su2"
DRY_RUN=false
SKIP_BUILD=false
TOPCHARGE_FILE_PARAMS=""
PARALLEL_JOBS=0
EXCLUDE_PATTERN=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --su2)              GAUGE_GROUP="su2";          shift ;;
        --su3)              GAUGE_GROUP="su3";          shift ;;
        --dry-run)          DRY_RUN=true;               shift ;;
        --skip-build)       SKIP_BUILD=true;            shift ;;
        --topcharge-file)   TOPCHARGE_FILE_PARAMS="$2"; shift 2 ;;
        --parallel)         PARALLEL_JOBS=$(nproc 2>/dev/null || echo 4); shift ;;
        --jobs)             PARALLEL_JOBS="$2";         shift 2 ;;
        --exclude)          EXCLUDE_PATTERN="$2";       shift 2 ;;
        -h|--help)          head -25 "$0" | tail -22;   exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ==============================================================================
# File selection
# ==============================================================================
if [ -z "$TOPCHARGE_FILE_PARAMS" ]; then
    TOPCHARGE_FILE_PARAMS="input/heatbath_input/topcharge_params_${GAUGE_GROUP}.txt"
    if [ ! -f "$TOPCHARGE_FILE_PARAMS" ] && [ -f "${TOPCHARGE_FILE_PARAMS%.txt}_finished.txt" ]; then
        TOPCHARGE_FILE_PARAMS="${TOPCHARGE_FILE_PARAMS%.txt}_finished.txt"
    fi
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

[ ! -f "$TOPCHARGE_FILE_PARAMS" ] && { log_error "Topcharge params not found: $TOPCHARGE_FILE_PARAMS"; exit 1; }

# Discover all run dirs for this gauge group
RUN_DIRS=($(ls -d "$RESULTS_DIR"/T*_seed*_${GAUGE_GROUP} 2>/dev/null))
if [ -n "$EXCLUDE_PATTERN" ]; then
    RUN_DIRS=($(printf '%s\n' "${RUN_DIRS[@]}" | grep -v "$EXCLUDE_PATTERN"))
fi
if [ ${#RUN_DIRS[@]} -eq 0 ]; then
    log_error "No run directories found matching *_${GAUGE_GROUP} in $RESULTS_DIR"
    exit 1
fi

log_info "Gauge group:       ${GAUGE_GROUP_UPPER}"
log_info "Topcharge params:  $TOPCHARGE_FILE_PARAMS"
log_info "Run dirs found:    ${#RUN_DIRS[@]}"
log_info "Results dir:       $RESULTS_DIR"
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
# Single-beta worker
# ==============================================================================
run_topcharge_beta() {
    local RUN_DIR="$1"
    local COUNT="$2"
    local TOTAL="$3"
    local TO_LOG="$4"    # "log" = file only, "tee" = also stdout

    local RUN_OUTPUT_DIR="$RUN_DIR/output"
    local TEMP_INPUT="$RUN_DIR/input.txt"
    local RUN_NAME
    RUN_NAME=$(basename "$RUN_DIR")
    local BETA
    BETA=$(read_param "beta" "$TEMP_INPUT")

    echo -e "${GREEN}------------------------------------------------------${NC}"
    echo -e "${GREEN}  [$COUNT/$TOTAL] beta = $BETA  ($RUN_NAME)${NC}"
    echo -e "${GREEN}------------------------------------------------------${NC}"

    if [ ! -f "$TEMP_INPUT" ]; then
        log_error "[$COUNT/$TOTAL] input.txt missing in $RUN_DIR"
        return 1
    fi

    # Patch input.txt with topcharge params
    local START_CONF END_CONF CONF_STEP SMEAR_STEPS SMEAR_ALPHA SMEAR_INTERVAL EXCLUDE_BC EXCLUDE_BC_OPEN
    # Prefer the thermalization cut written by run_heatbath_scan.sh; fall back to
    # the shared params file for legacy run directories.
    START_CONF=$(read_param "start_conf"                   "$TEMP_INPUT")
    START_CONF=${START_CONF:-$(read_param "start_conf"     "$TOPCHARGE_FILE_PARAMS")}
    END_CONF=$(read_param "end_conf"                       "$TOPCHARGE_FILE_PARAMS")
    CONF_STEP=$(read_param "conf_step"                     "$TOPCHARGE_FILE_PARAMS")
    SMEAR_STEPS=$(read_param "smear_steps"                 "$TOPCHARGE_FILE_PARAMS")
    SMEAR_ALPHA=$(read_param "smear_alpha"                 "$TOPCHARGE_FILE_PARAMS")
    SMEAR_INTERVAL=$(read_param "smear_interval"           "$TOPCHARGE_FILE_PARAMS")
    EXCLUDE_BC=$(read_param "exclude_boundary_slices"      "$TEMP_INPUT")
    EXCLUDE_BC_OPEN=$(read_param "exclude_boundary_slices_open" "$TOPCHARGE_FILE_PARAMS")

    # Preserve an explicit per-run temporal subvolume. Otherwise use the shared
    # open-boundary default for open runs and the full lattice for periodic runs.
    if [ -n "$EXCLUDE_BC" ]; then
        log_info "[$COUNT/$TOTAL] Preserving exclude_boundary_slices=$EXCLUDE_BC from input.txt"
    elif [[ "$RUN_NAME" == *"_open_"* ]]; then
        EXCLUDE_BC=${EXCLUDE_BC_OPEN:-2}
        log_info "[$COUNT/$TOTAL] Detected open BC — exclude_boundary_slices=$EXCLUDE_BC"
    else
        EXCLUDE_BC=0
        log_info "[$COUNT/$TOTAL] Detected periodic BC — exclude_boundary_slices=0"
    fi

    for key in start_conf end_conf conf_step smear_steps smear_alpha smear_interval exclude_boundary_slices; do
        grep -v "^${key}[[:space:]]" "$TEMP_INPUT" > "${TEMP_INPUT}.tmp" && mv "${TEMP_INPUT}.tmp" "$TEMP_INPUT"
    done

    cat >> "$TEMP_INPUT" << EOF

# Topcharge params (from $TOPCHARGE_FILE_PARAMS)
start_conf          ${START_CONF:-100}
end_conf            ${END_CONF:-1000}
conf_step           ${CONF_STEP:-10}
smear_steps         ${SMEAR_STEPS:-20}
smear_alpha         ${SMEAR_ALPHA:-0.3}
exclude_boundary_slices ${EXCLUDE_BC}
EOF
    if [ -n "$SMEAR_INTERVAL" ]; then
        echo "smear_interval      $SMEAR_INTERVAL" >> "$TEMP_INPUT"
    fi

    # Run measurement
    log_step "[$COUNT/$TOTAL] Measuring topcharge for beta=$BETA..."
    cd "$PROJECT_DIR"
    if [ "$TO_LOG" = "tee" ]; then
        "$BUILD_DIR/bin/$TOPCHARGE_BIN" -i "$TEMP_INPUT" 2>&1 | tee "$RUN_DIR/topcharge.log"
    else
        "$BUILD_DIR/bin/$TOPCHARGE_BIN" -i "$TEMP_INPUT" > "$RUN_DIR/topcharge.log" 2>&1
    fi

    local TOPCHARGE_PATH="$RUN_OUTPUT_DIR/$TOPCHARGE_FILE"
    if [ ! -f "$TOPCHARGE_PATH" ] && [ -f "$PROJECT_DIR/output/$TOPCHARGE_FILE" ]; then
        mv "$PROJECT_DIR/output/$TOPCHARGE_FILE" "$TOPCHARGE_PATH"
    fi

    cat >> "$RUN_DIR/run_info.txt" << EOF

# Topcharge phase
topcharge_file      $RUN_OUTPUT_DIR/$TOPCHARGE_FILE
EOF

    log_success "[$COUNT/$TOTAL] Beta=$BETA done — results in $RUN_OUTPUT_DIR"
}

# ==============================================================================
# Main loop
# ==============================================================================
TOTAL=${#RUN_DIRS[@]}

if [ "$DRY_RUN" = true ]; then
    for RD in "${RUN_DIRS[@]}"; do
        echo "  Would run topcharge on: $RD"
    done
    echo ""
    echo "Mode: $( [ "$PARALLEL_JOBS" -gt 0 ] && echo "parallel (max $PARALLEL_JOBS jobs)" || echo "sequential" )"
    exit 0
fi

if [ "$PARALLEL_JOBS" -gt 0 ]; then
    log_info "Parallel mode: running up to $PARALLEL_JOBS job(s) simultaneously"
    log_info "Progress visible in per-run topcharge.log files"
    echo ""

    declare -a PIDS=()
    declare -a PID_LABELS=()
    FAILED_COUNT=0

    for i in "${!RUN_DIRS[@]}"; do
        RD="${RUN_DIRS[$i]}"; COUNT=$((i+1))

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
                        log_error "Job for ${PID_LABELS[$j]} failed — check topcharge.log"
                        FAILED_COUNT=$((FAILED_COUNT+1))
                    fi
                fi
            done
            PIDS=("${NEW_PIDS[@]}"); PID_LABELS=("${NEW_LABELS[@]}")
            [ "${#PIDS[@]}" -ge "$PARALLEL_JOBS" ] && sleep 0.5
        done

        log_info "Launching job $COUNT/$TOTAL ($(basename "$RD"))..."
        run_topcharge_beta "$RD" "$COUNT" "$TOTAL" "log" &
        PIDS+=($!); PID_LABELS+=("$(basename "$RD")")
    done

    for j in "${!PIDS[@]}"; do
        p="${PIDS[$j]}"
        if wait "$p"; then
            :
        else
            log_error "Job for ${PID_LABELS[$j]} failed — check topcharge.log"
            FAILED_COUNT=$((FAILED_COUNT+1))
        fi
    done

    echo ""
    if [ "$FAILED_COUNT" -gt 0 ]; then
        log_warn "$FAILED_COUNT job(s) failed. Check individual topcharge.log files."
    fi

else
    for i in "${!RUN_DIRS[@]}"; do
        RD="${RUN_DIRS[$i]}"; COUNT=$((i+1))
        echo ""
        run_topcharge_beta "$RD" "$COUNT" "$TOTAL" "tee"
    done
fi

echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}  Topcharge scan complete. To run analysis:${NC}"
echo -e "${GREEN}  python3 analysis/analysis.py --${GAUGE_GROUP}${NC}"
echo -e "${GREEN}======================================================${NC}"
