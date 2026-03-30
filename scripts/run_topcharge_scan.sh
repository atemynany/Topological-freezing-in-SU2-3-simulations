#!/bin/bash
# ==============================================================================
# run_topcharge_scan.sh — Topological Charge Measurement (Step 2)
# ==============================================================================
# Scans existing run directories in data/results/, runs meas_topcharge on each
# using the input.txt already present (with detected start_conf from Step 1).
# Optionally runs analysis afterwards.
#
# Usage:
#   ./scripts/run_topcharge_scan.sh [--su2|--su3] [options]
#
# Options:
#   --su2             SU(2) gauge group (default)
#   --su3             SU(3) gauge group
#   --dry-run         Preview without executing
#   --skip-build      Skip build step
#   --skip-analysis   Skip running analysis at the end
#   --results-dir     Custom results directory
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

GAUGE_GROUP="su2"
DRY_RUN=false
SKIP_BUILD=false
SKIP_ANALYSIS=false
RESULTS_DIR="$PROJECT_DIR/data/results"
CONDA_ENV="master_thesis"

while [[ $# -gt 0 ]]; do
    case $1 in
        --su2)           GAUGE_GROUP="su2";   shift ;;
        --su3)           GAUGE_GROUP="su3";   shift ;;
        --dry-run)       DRY_RUN=true;        shift ;;
        --skip-build)    SKIP_BUILD=true;     shift ;;
        --skip-analysis) SKIP_ANALYSIS=true;  shift ;;
        --results-dir)   RESULTS_DIR="$2";    shift 2 ;;
        -h|--help)
            echo "Usage: $0 [--su2|--su3] [options]"
            echo "  --results-dir     Custom results directory"
            echo "  --dry-run         Preview without executing"
            echo "  --skip-build      Skip build step"
            echo "  --skip-analysis   Skip analysis at the end"
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [ "$GAUGE_GROUP" = "su2" ]; then
    TOPCHARGE_BIN="meas_topcharge"
    TOPCHARGE_FILE="topcharge.dat"
else
    TOPCHARGE_BIN="meas_topcharge_su3"
    TOPCHARGE_FILE="topcharge_su3.dat"
fi

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

# ==============================================================================
# Validation
# ==============================================================================
GAUGE_GROUP_UPPER=$(echo "$GAUGE_GROUP" | tr '[:lower:]' '[:upper:]')

echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}     ${GAUGE_GROUP_UPPER} Topological Charge Measurement${NC}"
echo -e "${GREEN}======================================================${NC}"
echo ""

if [ ! -d "$RESULTS_DIR" ]; then
    log_error "Results directory not found: $RESULTS_DIR"
    log_info "Run ./scripts/run_configs.sh first to generate configurations."
    exit 1
fi

# Find run directories that have an input.txt and configs
RUN_DIRS=()
for d in "$RESULTS_DIR"/T*_L*_b*_seed*/; do
    [ -d "$d" ] || continue
    [ -f "$d/input.txt" ] || continue
    [ -d "$d/configs" ] || continue
    # Skip runs that already have topcharge output
    if [ -f "$d/output/$TOPCHARGE_FILE" ]; then
        log_info "Skipping (already measured): $(basename "$d")"
        continue
    fi
    RUN_DIRS+=("$d")
done

if [ ${#RUN_DIRS[@]} -eq 0 ]; then
    log_warn "No unmeasured run directories found in $RESULTS_DIR"
    log_info "Either all runs are already measured, or run ./scripts/run_configs.sh first."
    exit 0
fi

log_info "Results dir:  $RESULTS_DIR"
log_info "Runs to measure: ${#RUN_DIRS[@]}"
echo ""

if [ "$DRY_RUN" = true ]; then
    log_warn "DRY RUN MODE"
    echo ""
fi

# ==============================================================================
# Build
# ==============================================================================
if [ "$SKIP_BUILD" = false ]; then
    log_step "Building $TOPCHARGE_BIN..."
    if [ "$DRY_RUN" = false ]; then
        mkdir -p "$PROJECT_DIR/build"
        cd "$PROJECT_DIR/build"
        cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
        make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4) "$TOPCHARGE_BIN" 2>&1 | tail -3
        cd "$PROJECT_DIR"
        log_success "Build complete"
    else
        echo "  Would run: cmake && make $TOPCHARGE_BIN"
    fi
    echo ""
fi

# ==============================================================================
# Main Loop
# ==============================================================================
RUN_COUNT=0
TOTAL_RUNS=${#RUN_DIRS[@]}

for RUN_DIR in "${RUN_DIRS[@]}"; do
    RUN_COUNT=$((RUN_COUNT + 1))
    RUN_NAME=$(basename "$RUN_DIR")

    echo ""
    echo -e "${GREEN}------------------------------------------------------${NC}"
    echo -e "${GREEN}  Run $RUN_COUNT/$TOTAL_RUNS: $RUN_NAME${NC}"
    echo -e "${GREEN}------------------------------------------------------${NC}"

    INPUT_FILE="$RUN_DIR/input.txt"

    if [ "$DRY_RUN" = true ]; then
        echo "  Would run: $TOPCHARGE_BIN -i $INPUT_FILE"
        continue
    fi

    log_step "Measuring topological charge..."
    cd "$PROJECT_DIR"
    ./build/bin/$TOPCHARGE_BIN -i "$INPUT_FILE" 2>&1 | tee "$RUN_DIR/topcharge.log"

    # Verify output
    TOPCHARGE_PATH="$(grep -E '^output_dir' "$INPUT_FILE" | awk '{print $2}')$TOPCHARGE_FILE"
    if [ -f "$TOPCHARGE_PATH" ]; then
        log_success "Topcharge saved: $TOPCHARGE_PATH"
    else
        log_warn "Topcharge output file not found at expected location"
    fi

    # Update run_info.txt
    if [ -f "$RUN_DIR/run_info.txt" ]; then
        if ! grep -q "topcharge_file" "$RUN_DIR/run_info.txt"; then
            echo "topcharge_file      $TOPCHARGE_PATH" >> "$RUN_DIR/run_info.txt"
        fi
    fi

    log_success "Run $RUN_COUNT/$TOTAL_RUNS complete: $RUN_NAME"
done

# ==============================================================================
# Analysis
# ==============================================================================
if [ "$DRY_RUN" = false ] && [ "$SKIP_ANALYSIS" = false ]; then
    echo ""
    log_step "Running analysis..."
    if command -v conda &> /dev/null && conda run -n "$CONDA_ENV" python3 --version &> /dev/null 2>&1; then
        conda run -n "$CONDA_ENV" python3 "$PROJECT_DIR/analysis/analysis.py" --${GAUGE_GROUP} --results-dir "$RESULTS_DIR" \
            || log_warn "Analysis failed (non-fatal)"
    elif command -v python3 &> /dev/null; then
        python3 "$PROJECT_DIR/analysis/analysis.py" --${GAUGE_GROUP} --results-dir "$RESULTS_DIR" \
            || log_warn "Analysis failed (non-fatal)"
    else
        log_warn "Python not found, skipping analysis"
    fi
fi

# ==============================================================================
# Summary
# ==============================================================================
echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}     Topological Charge Measurement Complete!${NC}"
echo -e "${GREEN}======================================================${NC}"
echo ""
log_info "Runs measured: $RUN_COUNT"
log_info "Results: $RESULTS_DIR"
echo ""
log_info "To re-run analysis:"
echo "  python3 analysis/analysis.py --${GAUGE_GROUP} --results-dir $RESULTS_DIR"
