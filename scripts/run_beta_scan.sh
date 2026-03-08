#!/bin/bash
# ==============================================================================
# Beta Scan Job Script for Lattice QCD Simulations
# ==============================================================================
# Runs multiple simulations scanning over beta values.
# For each beta:
#   1. Auto-generates unique seed from Unix timestamp
#   2. Runs MC heatbath simulation
#   3. Detects thermalization point from plaquette data
#   4. Rounds up thermalization to ensure uncorrelated data
#   5. Runs topological charge measurement
#   6. Saves all results to organized output folder
#
# Usage:
#   ./run_beta_scan.sh [--su2|--su3] [options]
#
# Options:
#   --su2           Use SU(2) (default)
#   --su3           Use SU(3)
#   --dry-run       Show what would be run without executing
#   --skip-build    Skip the build step
#   --beta-file     Specify custom beta scan file
#   --params-file   Specify custom base parameters file
#
# Examples:
#   ./run_beta_scan.sh --su2
#   ./run_beta_scan.sh --su3 --beta-file input/my_betas.txt
#   ./run_beta_scan.sh --su2 --dry-run
#
# ==============================================================================

set -e  # Exit on error

# ==============================================================================
# Default Configuration
# ==============================================================================
GAUGE_GROUP="su2"
DRY_RUN=false
SKIP_BUILD=false
BETA_FILE=""
PARAMS_FILE=""
CONDA_ENV="master_thesis"  # Conda environment with numpy, matplotlib, etc.

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# ==============================================================================
# Parse Arguments
# ==============================================================================
while [[ $# -gt 0 ]]; do
    case $1 in
        --su2)
            GAUGE_GROUP="su2"
            shift
            ;;
        --su3)
            GAUGE_GROUP="su3"
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --skip-build)
            SKIP_BUILD=true
            shift
            ;;
        --beta-file)
            BETA_FILE="$2"
            shift 2
            ;;
        --params-file)
            PARAMS_FILE="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [--su2|--su3] [options]"
            echo ""
            echo "Options:"
            echo "  --su2           Use SU(2) gauge group (default)"
            echo "  --su3           Use SU(3) gauge group"
            echo "  --dry-run       Show planned actions without executing"
            echo "  --skip-build    Skip the build step"
            echo "  --beta-file     Specify custom beta scan file"
            echo "  --params-file   Specify custom base parameters file"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# ==============================================================================
# Set Files Based on Gauge Group
# ==============================================================================
if [ "$GAUGE_GROUP" = "su2" ]; then
    BETA_FILE="${BETA_FILE:-$PROJECT_DIR/input/beta_scan_su2.txt}"
    PARAMS_FILE="${PARAMS_FILE:-$PROJECT_DIR/input/base_params_su2.txt}"
    HEATBATH_BIN="mc_heatbath"
    TOPCHARGE_BIN="meas_topcharge"
    PLAQ_FILE="plaquette.dat"
    TOPCHARGE_FILE="topcharge.dat"
    CONFIG_DIR_NAME="configs"
else
    BETA_FILE="${BETA_FILE:-$PROJECT_DIR/input/beta_scan_su3.txt}"
    PARAMS_FILE="${PARAMS_FILE:-$PROJECT_DIR/input/base_params_su3.txt}"
    HEATBATH_BIN="mc_heatbath_su3"
    TOPCHARGE_BIN="meas_topcharge_su3"
    PLAQ_FILE="plaquette_su3.dat"
    TOPCHARGE_FILE="topcharge_su3.dat"
    CONFIG_DIR_NAME="configs_su3"
fi

# Results directory
RESULTS_DIR="$PROJECT_DIR/data/results"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m'

# ==============================================================================
# Utility Functions
# ==============================================================================

log_info() {
    echo -e "${CYAN}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

log_step() {
    echo -e "${MAGENTA}[STEP]${NC} $1"
}

# Generate unique seed from Unix timestamp (10 digits)
generate_seed() {
    # Use nanoseconds for uniqueness
    local seed=$(date +%s%N | cut -c1-10)
    # Ensure it's a valid integer
    echo "$seed"
}

# Read a parameter from the base params file
read_param() {
    local param_name="$1"
    local file="$2"
    grep -E "^${param_name}\s+" "$file" | awk '{print $2}' | head -1
}

# Create temporary input file with all parameters
create_input_file() {
    local output_file="$1"
    local beta="$2"
    local seed="$3"
    local start_conf="$4"
    local run_output_dir="$5"
    local run_config_dir="$6"
    
    # Read base parameters
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
    
    # SU(3) specific
    local overrelax_steps=$(read_param "overrelax_steps" "$PARAMS_FILE")
    local smear_interval=$(read_param "smear_interval" "$PARAMS_FILE")
    
    cat > "$output_file" << EOF
# Auto-generated input file
# Beta scan: beta=$beta, seed=$seed
# Generated: $(date)

# Lattice
T                   ${T:-16}
L                   ${L:-16}
beta                $beta

# MC parameters
seed                $seed
start_type          ${start_type:-cold}
boundary            ${boundary:-periodic}
num_sweeps          ${num_sweeps:-1000}
save_interval       ${save_interval:-10}
EOF

    # Add SU(3) specific parameter
    if [ "$GAUGE_GROUP" = "su3" ] && [ -n "$overrelax_steps" ]; then
        echo "overrelax_steps     $overrelax_steps" >> "$output_file"
    fi

    cat >> "$output_file" << EOF

# Topological charge
start_conf          ${start_conf:-50}
end_conf            ${end_conf:-1000}
conf_step           ${conf_step:-10}
smear_steps         ${smear_steps:-20}
smear_alpha         ${smear_alpha:-0.3}
exclude_boundary_slices ${exclude_boundary_slices:-0}
EOF

    # Add smear_interval for SU(2)
    if [ "$GAUGE_GROUP" = "su2" ] && [ -n "$smear_interval" ]; then
        echo "smear_interval      $smear_interval" >> "$output_file"
    fi

    cat >> "$output_file" << EOF

# Output directories
output_dir          ${run_output_dir}/
config_dir          ${run_config_dir}/
EOF

    # SU(3) topcharge uses output_file instead of output_dir
    if [ "$GAUGE_GROUP" = "su3" ]; then
        echo "output_file         ${run_output_dir}/topcharge_su3.dat" >> "$output_file"
    fi
}

# ==============================================================================
# Validation
# ==============================================================================

echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}     ${GAUGE_GROUP^^} Beta Scan - Lattice QCD Simulation${NC}"
echo -e "${GREEN}======================================================${NC}"
echo ""

# Check files exist
if [ ! -f "$BETA_FILE" ]; then
    log_error "Beta scan file not found: $BETA_FILE"
    exit 1
fi

if [ ! -f "$PARAMS_FILE" ]; then
    log_error "Base parameters file not found: $PARAMS_FILE"
    exit 1
fi

# Read configuration
SAVE_INTERVAL=$(read_param "save_interval" "$PARAMS_FILE")
T_SIZE=$(read_param "T" "$PARAMS_FILE")
L_SIZE=$(read_param "L" "$PARAMS_FILE")

log_info "Gauge group:     ${GAUGE_GROUP^^}"
log_info "Beta file:       $BETA_FILE"
log_info "Params file:     $PARAMS_FILE"
log_info "Lattice:         ${T_SIZE}x${L_SIZE}^3"
log_info "Save interval:   $SAVE_INTERVAL"
log_info "Results dir:     $RESULTS_DIR"
echo ""

# Read beta values (skip comments and empty lines)
BETAS=($(grep -v '^#' "$BETA_FILE" | grep -v '^$'))
log_info "Beta values to scan: ${BETAS[*]}"
log_info "Total runs: ${#BETAS[@]}"
echo ""

if [ "$DRY_RUN" = true ]; then
    log_warn "DRY RUN MODE - No commands will be executed"
    echo ""
fi

# ==============================================================================
# Build (if not skipped)
# ==============================================================================
if [ "$SKIP_BUILD" = false ]; then
    log_step "Building project..."
    if [ "$DRY_RUN" = false ]; then
        mkdir -p "$PROJECT_DIR/build"
        cd "$PROJECT_DIR/build"
        cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
        make -j$(nproc) $HEATBATH_BIN $TOPCHARGE_BIN 2>&1 | tail -3
        cd "$PROJECT_DIR"
        log_success "Build complete"
    else
        echo "  Would run: cmake && make $HEATBATH_BIN $TOPCHARGE_BIN"
    fi
    echo ""
fi

# ==============================================================================
# Create results directory
# ==============================================================================
mkdir -p "$RESULTS_DIR"

# ==============================================================================
# Main Loop: Process each beta value
# ==============================================================================
RUN_COUNT=0
TOTAL_RUNS=${#BETAS[@]}

for BETA in "${BETAS[@]}"; do
    RUN_COUNT=$((RUN_COUNT + 1))
    
    echo ""
    echo -e "${GREEN}------------------------------------------------------${NC}"
    echo -e "${GREEN}  Run $RUN_COUNT/$TOTAL_RUNS: beta = $BETA${NC}"
    echo -e "${GREEN}------------------------------------------------------${NC}"
    
    # Generate unique seed
    SEED=$(generate_seed)
    sleep 0.1  # Ensure unique seeds if running fast
    
    # Create run-specific output directory
    RUN_NAME="T${T_SIZE}_L${L_SIZE}_b${BETA}_seed${SEED}"
    RUN_DIR="$RESULTS_DIR/$RUN_NAME"
    RUN_OUTPUT_DIR="$RUN_DIR/output"
    RUN_CONFIG_DIR="$RUN_DIR/configs"
    
    log_info "Run name:     $RUN_NAME"
    log_info "Seed:         $SEED"
    log_info "Output dir:   $RUN_DIR"
    
    if [ "$DRY_RUN" = true ]; then
        echo "  Would create: $RUN_DIR"
        echo "  Would run MC heatbath with beta=$BETA, seed=$SEED"
        echo "  Would detect thermalization from plaquette data"
        echo "  Would run topological charge measurement"
        continue
    fi
    
    # Create directories
    mkdir -p "$RUN_OUTPUT_DIR"
    mkdir -p "$RUN_CONFIG_DIR"
    
    # Create temporary input file for heatbath (start_conf=10 placeholder)
    TEMP_INPUT="$RUN_DIR/input.txt"
    create_input_file "$TEMP_INPUT" "$BETA" "$SEED" "10" "$RUN_OUTPUT_DIR" "$RUN_CONFIG_DIR"
    
    # -------------------------------------------------------------------------
    # Step 1: Run MC Heatbath
    # -------------------------------------------------------------------------
    log_step "Running MC heatbath simulation..."
    cd "$PROJECT_DIR"
    
    # Run heatbath (show progress bar, also save to log)
    ./build/bin/$HEATBATH_BIN -i "$TEMP_INPUT" 2>&1 | tee "$RUN_DIR/heatbath.log"
    
    log_success "Heatbath complete"
    
    # -------------------------------------------------------------------------
    # Step 2: Detect thermalization
    # -------------------------------------------------------------------------
    log_step "Detecting thermalization point..."
    
    PLAQ_PATH="$RUN_OUTPUT_DIR/$PLAQ_FILE"
    if [ ! -f "$PLAQ_PATH" ]; then
        log_error "Plaquette file not found: $PLAQ_PATH"
        log_warn "Checking default output location..."
        # Fallback: check if it's in the default location
        if [ -f "$PROJECT_DIR/output/$PLAQ_FILE" ]; then
            mkdir -p "$RUN_OUTPUT_DIR"
            mv "$PROJECT_DIR/output/$PLAQ_FILE" "$PLAQ_PATH"
            log_info "Moved plaquette file from default location"
        else
            log_error "Cannot find plaquette file, skipping this run"
            continue
        fi
    fi
    
    # Run thermalization detection script (using conda environment)
    START_CONF=$(conda run -n "$CONDA_ENV" python3 "$SCRIPT_DIR/detect_thermalization.py" "$PLAQ_PATH" "$SAVE_INTERVAL" 2>&1 | tail -1)
    
    if ! [[ "$START_CONF" =~ ^[0-9]+$ ]]; then
        log_warn "Could not detect thermalization, using default start_conf=50"
        START_CONF=50
    fi
    
    log_info "Using start_conf=$START_CONF (conservative, rounded up)"
    
    # -------------------------------------------------------------------------
    # Step 3: Update input file with detected start_conf
    # -------------------------------------------------------------------------
    create_input_file "$TEMP_INPUT" "$BETA" "$SEED" "$START_CONF" "$RUN_OUTPUT_DIR" "$RUN_CONFIG_DIR"
    
    # -------------------------------------------------------------------------
    # Step 4: Run topological charge measurement
    # -------------------------------------------------------------------------
    log_step "Measuring topological charge..."
    
    # Run topcharge measurement (show progress bar, also save to log)
    ./build/bin/$TOPCHARGE_BIN -i "$TEMP_INPUT" 2>&1 | tee "$RUN_DIR/topcharge.log"
    
    # Verify output file exists (either in run dir or default location)
    TOPCHARGE_PATH="$RUN_OUTPUT_DIR/$TOPCHARGE_FILE"
    if [ ! -f "$TOPCHARGE_PATH" ]; then
        # Fallback: check default location
        if [ -f "$PROJECT_DIR/output/$TOPCHARGE_FILE" ]; then
            mv "$PROJECT_DIR/output/$TOPCHARGE_FILE" "$TOPCHARGE_PATH"
            log_info "Moved topcharge file from default location"
        else
            log_warn "Topcharge output file not found"
        fi
    fi
    
    log_success "Topological charge measurement complete"
    
    # -------------------------------------------------------------------------
    # Step 5: Save run metadata
    # -------------------------------------------------------------------------
    cat > "$RUN_DIR/run_info.txt" << EOF
# Run Information
# Generated: $(date)
# ==============================================================================

gauge_group         ${GAUGE_GROUP^^}
beta                $BETA
seed                $SEED
T                   $T_SIZE
L                   $L_SIZE
save_interval       $SAVE_INTERVAL
start_conf          $START_CONF
run_dir             $RUN_DIR

# Files
plaquette_file      $RUN_OUTPUT_DIR/$PLAQ_FILE
topcharge_file      $RUN_OUTPUT_DIR/$TOPCHARGE_FILE
config_dir          $RUN_CONFIG_DIR
EOF
    
    log_success "Run $RUN_COUNT/$TOTAL_RUNS complete: $RUN_NAME"
done

# ==============================================================================
# Summary
# ==============================================================================
echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}     Beta Scan Complete!${NC}"
echo -e "${GREEN}======================================================${NC}"
echo ""
log_info "Total runs completed: $RUN_COUNT"
log_info "Results saved to: $RESULTS_DIR"
echo ""
echo "Run directories:"
for BETA in "${BETAS[@]}"; do
    ls -d "$RESULTS_DIR"/T${T_SIZE}_L${L_SIZE}_b${BETA}_* 2>/dev/null || true
done
echo ""

# Run analysis automatically if not dry run
if [ "$DRY_RUN" = false ]; then
    log_step "Running analysis..."
    if command -v conda &> /dev/null && conda run -n "$CONDA_ENV" python3 --version &> /dev/null; then
        conda run -n "$CONDA_ENV" python3 "$PROJECT_DIR/analysis/analyze_beta_scan.py" --${GAUGE_GROUP} --results-dir "$RESULTS_DIR" || log_warn "Analysis failed (non-fatal)"
    elif command -v python3 &> /dev/null; then
        python3 "$PROJECT_DIR/analysis/analyze_beta_scan.py" --${GAUGE_GROUP} --results-dir "$RESULTS_DIR" || log_warn "Analysis failed (non-fatal)"
    else
        log_warn "Python not found, skipping analysis"
    fi
fi

log_info "To re-run analysis:"
echo "  python3 analysis/analyze_beta_scan.py --${GAUGE_GROUP} --results-dir $RESULTS_DIR"
