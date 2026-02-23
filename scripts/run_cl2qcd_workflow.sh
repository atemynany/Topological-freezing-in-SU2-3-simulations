#!/bin/bash
# run_cl2qcd_workflow.sh - Complete workflow for generating SU(3) configs with CL2QCD
#
# This script:
# 1. Reads the unified input file (cl2qcd_heatbath.txt format)
# 2. Generates CL2QCD native input format
# 3. Runs CL2QCD to generate gauge configurations (LIME format)
# 4. Converts LIME configs to our raw binary format
# 5. Optionally measures topological charge
#
# Usage:
#   ./scripts/run_cl2qcd_workflow.sh input/cl2qcd_heatbath.txt [--measure-topcharge]
#
# Prerequisites:
#   - Run ./scripts/build_cl2qcd.sh first to build executables

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="$PROJECT_DIR/build/bin"

# -----------------------------------------------------------------------------
# Parse arguments
# -----------------------------------------------------------------------------
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_file> [--measure-topcharge]"
    echo "Example: $0 input/cl2qcd_heatbath.txt --measure-topcharge"
    exit 1
fi

INPUT_FILE="$1"
MEASURE_TOPCHARGE=false

for arg in "$@"; do
    case $arg in
        --measure-topcharge)
            MEASURE_TOPCHARGE=true
            ;;
    esac
done

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found: $INPUT_FILE"
    exit 1
fi

# Check executables
if [ ! -f "$BUILD_DIR/cl2qcd_heatbath" ]; then
    echo "Error: cl2qcd_heatbath not found. Run ./scripts/build_cl2qcd.sh first"
    exit 1
fi

if [ ! -f "$BUILD_DIR/lime_to_raw" ]; then
    echo "Error: lime_to_raw not found. Run ./scripts/build_cl2qcd.sh first"
    exit 1
fi

echo "=============================================="
echo "CL2QCD Workflow"
echo "=============================================="
echo "Input file: $INPUT_FILE"
echo "Measure topological charge: $MEASURE_TOPCHARGE"
echo

# -----------------------------------------------------------------------------
# Parse input file
# -----------------------------------------------------------------------------
echo "[1/5] Parsing input file..."

parse_value() {
    grep "^$1" "$INPUT_FILE" | awk '{print $2}' | tr -d '\r'
}

NSPACE=$(parse_value "nSpace")
NTIME=$(parse_value "nTime")
BETA=$(parse_value "beta")
N_THERM=$(parse_value "nThermalizationSteps")
N_HEATBATH=$(parse_value "nHeatbathSteps")
SAVE_FREQ=$(parse_value "saveFrequency")
OUTPUT_DIR=$(parse_value "outputDir")
SMEAR_STEPS=$(parse_value "smearSteps")
SMEAR_RHO=$(parse_value "smearRho")

# Set defaults
NSPACE=${NSPACE:-16}
NTIME=${NTIME:-16}
BETA=${BETA:-6.0}
N_THERM=${N_THERM:-1000}
N_HEATBATH=${N_HEATBATH:-10000}
SAVE_FREQ=${SAVE_FREQ:-100}
OUTPUT_DIR=${OUTPUT_DIR:-"output/configs_cl2qcd"}
SMEAR_STEPS=${SMEAR_STEPS:-40}
SMEAR_RHO=${SMEAR_RHO:-0.12}

echo "  Lattice: ${NSPACE}^3 x $NTIME"
echo "  beta: $BETA"
echo "  Thermalization: $N_THERM"
echo "  MC steps: $N_HEATBATH"
echo "  Save frequency: $SAVE_FREQ"
echo "  Output dir: $OUTPUT_DIR"
echo

# Create directories
LIME_DIR="${OUTPUT_DIR}/lime"
RAW_DIR="${OUTPUT_DIR}/raw"
mkdir -p "$LIME_DIR" "$RAW_DIR"

# -----------------------------------------------------------------------------
# Generate CL2QCD input file
# -----------------------------------------------------------------------------
echo "[2/5] Generating CL2QCD input file..."

CL2QCD_INPUT="$OUTPUT_DIR/cl2qcd_input.txt"

cat > "$CL2QCD_INPUT" << EOF
# CL2QCD input file - auto-generated
# Generated from: $INPUT_FILE
# Date: $(date)

# Lattice dimensions
nSpace = $NSPACE
nTime = $NTIME

# Action parameters
beta = $BETA
fermionAction = none

# HMC/Heatbath parameters
startCondition = cold
nThermalizationSteps = $N_THERM
nConfigurations = $((N_HEATBATH / SAVE_FREQ))
writeFrequency = 1

# I/O
outputDirectory = $LIME_DIR
configurationPrefix = conf_
writeInIldgFormat = true

# Boundary conditions (periodic for now)
theta_fermion_temporal = 0.0
theta_fermion_spatial = 0.0

# Precision
precision = double
EOF

echo "  Generated: $CL2QCD_INPUT"
echo

# -----------------------------------------------------------------------------
# Run CL2QCD heatbath
# -----------------------------------------------------------------------------
echo "[3/5] Running CL2QCD heatbath simulation..."
echo "  This may take a while..."
echo

cd "$PROJECT_DIR"
"$BUILD_DIR/cl2qcd_heatbath" "$CL2QCD_INPUT" 2>&1 | tee "$OUTPUT_DIR/cl2qcd.log"

echo
echo "  Simulation complete. Check $OUTPUT_DIR/cl2qcd.log for details."
echo

# -----------------------------------------------------------------------------
# Convert LIME to raw format
# -----------------------------------------------------------------------------
echo "[4/5] Converting LIME configs to raw format..."

# Find all LIME configs
LIME_FILES=$(find "$LIME_DIR" -name "*.lime" -o -name "conf_*" 2>/dev/null | sort)
NUM_CONFIGS=0

for lime_file in $LIME_FILES; do
    if [ -f "$lime_file" ]; then
        # Extract config number from filename
        basename_file=$(basename "$lime_file")
        # Try to extract number
        config_num=$(echo "$basename_file" | grep -oE '[0-9]+' | tail -1)
        
        if [ -n "$config_num" ]; then
            raw_file="$RAW_DIR/config_${config_num}.bin"
        else
            raw_file="$RAW_DIR/${basename_file%.lime}.bin"
        fi
        
        echo "  Converting: $basename_file -> $(basename $raw_file)"
        "$BUILD_DIR/lime_to_raw" "$lime_file" "$raw_file"
        
        NUM_CONFIGS=$((NUM_CONFIGS + 1))
    fi
done

echo
echo "  Converted $NUM_CONFIGS configurations"
echo

# -----------------------------------------------------------------------------
# Measure topological charge (optional)
# -----------------------------------------------------------------------------
if [ "$MEASURE_TOPCHARGE" = true ]; then
    echo "[5/5] Measuring topological charge..."
    
    # Check if meas_topcharge_su3 exists
    if [ ! -f "$BUILD_DIR/meas_topcharge_su3" ]; then
        echo "  Error: meas_topcharge_su3 not found"
        echo "  Run 'cmake --build build' to build it"
    else
        # Create topcharge input file
        TOPCHARGE_INPUT="$OUTPUT_DIR/topcharge_input.txt"
        TOPCHARGE_OUTPUT="$OUTPUT_DIR/topcharge_results.dat"
        
        cat > "$TOPCHARGE_INPUT" << EOF
# Topological charge measurement input
# Auto-generated by run_cl2qcd_workflow.sh

nSpace $NSPACE
nTime $NTIME
beta $BETA
configDir $RAW_DIR
outputFile $TOPCHARGE_OUTPUT
smearSteps $SMEAR_STEPS
smearRho $SMEAR_RHO
EOF
        
        echo "  Running topological charge measurement..."
        
        # Measure for each config
        for raw_file in "$RAW_DIR"/*.bin; do
            if [ -f "$raw_file" ]; then
                config_num=$(basename "$raw_file" .bin | grep -oE '[0-9]+')
                echo "  Config $config_num..."
                "$BUILD_DIR/meas_topcharge_su3" "$TOPCHARGE_INPUT" "$raw_file" >> "$TOPCHARGE_OUTPUT"
            fi
        done
        
        echo "  Results written to: $TOPCHARGE_OUTPUT"
    fi
else
    echo "[5/5] Skipping topological charge measurement (use --measure-topcharge to enable)"
fi

echo
echo "=============================================="
echo "Workflow Complete!"
echo "=============================================="
echo
echo "Output:"
echo "  LIME configs: $LIME_DIR/"
echo "  Raw configs:  $RAW_DIR/"
if [ "$MEASURE_TOPCHARGE" = true ]; then
    echo "  Q results:    $OUTPUT_DIR/topcharge_results.dat"
fi
echo
echo "To analyze topological charge:"
echo "  python analysis/run_analysis_su3.py $OUTPUT_DIR/topcharge_results.dat"
echo
