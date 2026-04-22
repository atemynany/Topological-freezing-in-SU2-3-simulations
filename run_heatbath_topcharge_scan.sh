#!/bin/bash
# ==============================================================================
# run_heatbath_topcharge_scan.sh  —  Fused SU(2) heatbath + in-memory topcharge
# ==============================================================================
# Replaces run_heatbath_scan.sh + run_topcharge_scan.sh for SU(2) when you don't
# want to persist gauge configs to disk. Discovers input/heatbath_topcharge_input/setup_*_su2.txt
# files, merges each with input/heatbath_input/topcharge_params_su2.txt, and runs
# heatbath_topcharge_su2 for every (setup_file, beta) combination.
#
# Usage:
#   ./run_heatbath_topcharge_scan.sh [options]
#
# Options:
#   --dry-run       Preview without executing
#   --skip-build    Skip build step
#   --jobs N        Run up to N beta values concurrently (default: 1)
#   --exclude PAT   Skip setup files whose path matches regex PAT
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"
BUILD_DIR="$PROJECT_DIR/build"
RESULTS_DIR="${HEATBATH_RESULTS_DIR:-$PROJECT_DIR/data/results}"

BINARY="heatbath_topcharge_su2"
TOPCHARGE_PARAMS="$PROJECT_DIR/input/heatbath_input/topcharge_params_su2.txt"

DRY_RUN=false
SKIP_BUILD=false
PARALLEL_JOBS=1
EXCLUDE_PAT=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)    DRY_RUN=true;    shift ;;
        --skip-build) SKIP_BUILD=true; shift ;;
        --jobs)       PARALLEL_JOBS="$2"; shift 2 ;;
        --exclude)    EXCLUDE_PAT="$2";   shift 2 ;;
        -h|--help)    head -18 "$0" | tail -15; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

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
    local output_file="$1" beta="$2" seed="$3" run_output_dir="$4" setup_file="$5"

    local T=$(read_param "T" "$setup_file")
    local L=$(read_param "L" "$setup_file")
    local start_type=$(read_param "start_type" "$setup_file")
    local boundary=$(read_param "boundary" "$setup_file")
    local num_sweeps=$(read_param "num_sweeps" "$setup_file")
    local save_interval=$(read_param "save_interval" "$setup_file")

    # Measurement params from shared topcharge_params file
    local smear_steps=$(read_param "smear_steps" "$TOPCHARGE_PARAMS")
    local smear_interval=$(read_param "smear_interval" "$TOPCHARGE_PARAMS")
    local smear_alpha=$(read_param "smear_alpha" "$TOPCHARGE_PARAMS")
    local exclude_boundary_slices_open=$(read_param "exclude_boundary_slices_open" "$TOPCHARGE_PARAMS")

    # Per-ensemble override in the setup file takes priority over the global default.
    # Convention: open setup files carry exclude_boundary_slices = (T_open - T_periodic)/2.
    local setup_exclude=$(read_param "exclude_boundary_slices" "$setup_file")

    local exclude_boundary_slices=0
    if [ -n "$setup_exclude" ]; then
        exclude_boundary_slices="$setup_exclude"
    elif [ "$boundary" = "open" ]; then
        exclude_boundary_slices="${exclude_boundary_slices_open:-0}"
    fi

    cat > "$output_file" << EOF
# Auto-generated input for heatbath_topcharge_su2
# setup: $(basename "$setup_file"), beta=$beta, seed=$seed
# generated: $(date)

T                        ${T:-16}
L                        ${L:-16}
beta                     $beta
seed                     $seed
start_type               ${start_type:-cold}
boundary                 ${boundary:-periodic}
num_sweeps               ${num_sweeps:-1000}
save_interval            ${save_interval:-10}

smear_steps              ${smear_steps:-40}
smear_interval           ${smear_interval:-10}
smear_alpha              ${smear_alpha:-0.5}
exclude_boundary_slices  $exclude_boundary_slices

output_dir               ${run_output_dir}/
EOF
}

# ==============================================================================
# Discover setup files
# ==============================================================================
if [ ! -f "$TOPCHARGE_PARAMS" ]; then
    log_error "Missing topcharge params file: $TOPCHARGE_PARAMS"
    exit 1
fi

SETUP_FILES=("$PROJECT_DIR"/input/heatbath_topcharge_input/setup_*_su2.txt)
if [ ! -f "${SETUP_FILES[0]}" ]; then
    log_error "No setup files found matching input/heatbath_topcharge_input/setup_*_su2.txt"
    exit 1
fi

declare -a JOB_SETUP=()
declare -a JOB_BETA=()

for SETUP in "${SETUP_FILES[@]}"; do
    if [ -n "$EXCLUDE_PAT" ] && [[ "$SETUP" =~ $EXCLUDE_PAT ]]; then
        log_info "Skipping $(basename "$SETUP") (excluded)"
        continue
    fi
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
    log_error "No jobs found"
    exit 1
fi

echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}    SU(2) Heatbath + In-Memory Topcharge Scan${NC}"
echo -e "${GREEN}======================================================${NC}"
log_info "Binary:        $BINARY"
log_info "Topcharge params: $(basename "$TOPCHARGE_PARAMS")"
log_info "Total jobs:    $TOTAL"
log_info "Results dir:   $RESULTS_DIR"
log_info "Parallel jobs: $PARALLEL_JOBS"
echo ""

# ==============================================================================
# Build
# ==============================================================================
if [ "$SKIP_BUILD" = false ] && [ "$DRY_RUN" = false ]; then
    log_step "Building $BINARY..."
    mkdir -p "$BUILD_DIR" && cd "$BUILD_DIR"
    cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
    make -j$(nproc 2>/dev/null || echo 4) "$BINARY" 2>&1 | tail -3
    cd "$PROJECT_DIR"
    log_success "Build complete"
    echo ""
fi

mkdir -p "$RESULTS_DIR"

# Pre-generate seeds
declare -a SEEDS=()
for i in "${!JOB_BETA[@]}"; do
    SEEDS+=("$(generate_seed)")
    sleep 0.1
done

# ==============================================================================
# Single-beta worker
# ==============================================================================
run_single_beta() {
    local SETUP_FILE="$1" BETA="$2" SEED="$3" RUN_COUNT="$4" TO_LOG="$5"

    local T_SIZE=$(read_param "T" "$SETUP_FILE")
    local L_SIZE=$(read_param "L" "$SETUP_FILE")
    local BOUNDARY=$(read_param "boundary" "$SETUP_FILE"); BOUNDARY=${BOUNDARY:-periodic}

    local RUN_NAME="T${T_SIZE}_L${L_SIZE}_b${BETA}_${BOUNDARY}_seed${SEED}_su2"
    local RUN_DIR="$RESULTS_DIR/$RUN_NAME"
    local RUN_OUTPUT_DIR="$RUN_DIR/output"

    echo -e "${GREEN}------------------------------------------------------${NC}"
    echo -e "${GREEN}  [$RUN_COUNT/$TOTAL] beta=$BETA  (${T_SIZE}x${L_SIZE}^3, ${BOUNDARY}, seed=$SEED)${NC}"
    echo -e "${GREEN}------------------------------------------------------${NC}"

    mkdir -p "$RUN_OUTPUT_DIR"
    local INPUT_FILE="$RUN_DIR/input.txt"
    create_input_file "$INPUT_FILE" "$BETA" "$SEED" "$RUN_OUTPUT_DIR" "$SETUP_FILE"

    cd "$PROJECT_DIR"
    if [ "$TO_LOG" = "tee" ]; then
        "$BUILD_DIR/bin/$BINARY" -i "$INPUT_FILE" 2>&1 | tee "$RUN_DIR/run.log"
    else
        "$BUILD_DIR/bin/$BINARY" -i "$INPUT_FILE" > "$RUN_DIR/run.log" 2>&1
    fi
    local rc=${PIPESTATUS[0]}
    if [ "$rc" -ne 0 ]; then
        log_error "[$RUN_COUNT/$TOTAL] $BINARY exited with code $rc for beta=$BETA"
        return 1
    fi

    cat > "$RUN_DIR/run_info.txt" << EOF
# Fused heatbath + topcharge run
# Generated: $(date)

gauge_group         SU2
beta                $BETA
seed                $SEED
T                   $T_SIZE
L                   $L_SIZE
boundary            $BOUNDARY
run_dir             $RUN_DIR

# Files
plaquette_file      $RUN_OUTPUT_DIR/plaquette.dat
topcharge_file      $RUN_OUTPUT_DIR/topcharge.dat
timeslice_file      $RUN_OUTPUT_DIR/topcharge_timeslice.dat
EOF

    log_success "[$RUN_COUNT/$TOTAL] beta=$BETA done — output in $RUN_OUTPUT_DIR"
}

# ==============================================================================
# Dry-run preview
# ==============================================================================
if [ "$DRY_RUN" = true ]; then
    for i in "${!JOB_BETA[@]}"; do
        SETUP="${JOB_SETUP[$i]}"; BETA="${JOB_BETA[$i]}"; SEED="${SEEDS[$i]}"
        T_SIZE=$(read_param "T" "$SETUP"); L_SIZE=$(read_param "L" "$SETUP")
        BOUNDARY=$(read_param "boundary" "$SETUP"); BOUNDARY=${BOUNDARY:-periodic}
        RUN_NAME="T${T_SIZE}_L${L_SIZE}_b${BETA}_${BOUNDARY}_seed${SEED}_su2"
        echo "  Would create: $RESULTS_DIR/$RUN_NAME"
        echo "    Setup: $(basename "$SETUP"), beta=$BETA, seed=$SEED"
    done
    echo ""
    echo "Mode: $( [ "$PARALLEL_JOBS" -gt 1 ] && echo "parallel (max $PARALLEL_JOBS jobs)" || echo "sequential" )"
    exit 0
fi

# ==============================================================================
# Main loop
# ==============================================================================
if [ "$PARALLEL_JOBS" -gt 1 ]; then
    log_info "Parallel mode: up to $PARALLEL_JOBS concurrent job(s)"
    log_info "Per-run output in <run_dir>/run.log"
    echo ""

    declare -a PIDS=()
    declare -a PID_LABELS=()
    FAILED_COUNT=0

    for i in "${!JOB_BETA[@]}"; do
        SETUP="${JOB_SETUP[$i]}"; BETA="${JOB_BETA[$i]}"; SEED="${SEEDS[$i]}"
        COUNT=$((i+1))

        while [ "${#PIDS[@]}" -ge "$PARALLEL_JOBS" ]; do
            NEW_PIDS=(); NEW_LABELS=()
            for j in "${!PIDS[@]}"; do
                p="${PIDS[$j]}"
                if kill -0 "$p" 2>/dev/null; then
                    NEW_PIDS+=("$p"); NEW_LABELS+=("${PID_LABELS[$j]}")
                else
                    if ! wait "$p"; then
                        log_error "Job for ${PID_LABELS[$j]} failed — check run.log"
                        FAILED_COUNT=$((FAILED_COUNT+1))
                    fi
                fi
            done
            PIDS=("${NEW_PIDS[@]}"); PID_LABELS=("${NEW_LABELS[@]}")
            [ "${#PIDS[@]}" -ge "$PARALLEL_JOBS" ] && sleep 0.5
        done

        log_info "Launching beta=$BETA [job $COUNT/$TOTAL]..."
        run_single_beta "$SETUP" "$BETA" "$SEED" "$COUNT" "log" &
        PIDS+=($!); PID_LABELS+=("beta=$BETA")
    done

    for j in "${!PIDS[@]}"; do
        p="${PIDS[$j]}"
        if ! wait "$p"; then
            log_error "Job for ${PID_LABELS[$j]} failed — check run.log"
            FAILED_COUNT=$((FAILED_COUNT+1))
        fi
    done

    echo ""
    if [ "$FAILED_COUNT" -gt 0 ]; then
        log_warn "$FAILED_COUNT job(s) failed. Check run.log files."
    fi
else
    for i in "${!JOB_BETA[@]}"; do
        SETUP="${JOB_SETUP[$i]}"; BETA="${JOB_BETA[$i]}"; SEED="${SEEDS[$i]}"
        COUNT=$((i+1))
        echo ""
        run_single_beta "$SETUP" "$BETA" "$SEED" "$COUNT" "tee"
    done
fi

echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}  Heatbath + topcharge scan complete.${NC}"
echo -e "${GREEN}  Results in: $RESULTS_DIR${NC}"
echo -e "${GREEN}======================================================${NC}"
