#!/bin/bash
# ==============================================================================
# Lattice Gauge Theory Simulation Runner (SU(2) / SU(3))
# ==============================================================================
# This script runs the full simulation pipeline:
# 1. Build the project
# 2. Generate gauge configurations (MC heatbath)
# 3. Measure topological charge
# 4. Run analysis and plot results
#
# Usage:
#   ./run_simulation.sh [--su2|--su3] [input_file]
#
# Examples:
#   ./run_simulation.sh --su2 input/simulation.txt
#   ./run_simulation.sh --su3 input/simulation_su3.txt
#   ./run_simulation.sh                              # defaults to SU(2)
#
# ==============================================================================

set -e  # Exit on error

# Default values
GAUGE_GROUP="su2"
INPUT_FILE=""

# Parse arguments
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
        -h|--help)
            echo "Usage: $0 [--su2|--su3] [input_file]"
            echo ""
            echo "Options:"
            echo "  --su2    Run SU(2) simulation (default)"
            echo "  --su3    Run SU(3) simulation"
            echo ""
            echo "Default input files:"
            echo "  SU(2): input/simulation.txt"
            echo "  SU(3): input/simulation_su3.txt"
            exit 0
            ;;
        *)
            INPUT_FILE="$1"
            shift
            ;;
    esac
done

# Set defaults based on gauge group
if [ "$GAUGE_GROUP" = "su2" ]; then
    INPUT_FILE="${INPUT_FILE:-input/simulation.txt}"
    HEATBATH_BIN="mc_heatbath"
    TOPCHARGE_BIN="meas_topcharge"
    PLAQ_FILE="output/plaquette.dat"
    TOPCHARGE_FILE="output/topcharge.dat"
    ANALYSIS_SCRIPT="analysis/run_analysis.py"
else
    INPUT_FILE="${INPUT_FILE:-input/simulation_su3.txt}"
    HEATBATH_BIN="mc_heatbath_su3"
    TOPCHARGE_BIN="meas_topcharge_su3"
    PLAQ_FILE="output/plaquette_su3.dat"
    TOPCHARGE_FILE="output/topcharge_su3.dat"
    ANALYSIS_SCRIPT="analysis/run_analysis_su3.py"
fi

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

echo -e "${GREEN}================================================${NC}"
echo -e "${GREEN}${GAUGE_GROUP^^} Lattice Gauge Theory Simulation${NC}"
echo -e "${GREEN}================================================${NC}"
echo ""

# Check input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Error: Input file not found: $INPUT_FILE${NC}"
    exit 1
fi

echo -e "${CYAN}Gauge group:  ${GAUGE_GROUP^^}${NC}"
echo -e "${CYAN}Input file:   $INPUT_FILE${NC}"
echo ""

# ==============================================================================
# Step 1: Build
# ==============================================================================
echo -e "${GREEN}[Step 1/4] Building project...${NC}"
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null
make -j$(nproc) $HEATBATH_BIN $TOPCHARGE_BIN 2>&1 | tail -5
cd ..
echo -e "${GREEN}Build complete.${NC}"
echo ""

# ==============================================================================
# Step 2: Generate configurations
# ==============================================================================
echo -e "${GREEN}[Step 2/4] Generating gauge configurations...${NC}"
./build/bin/$HEATBATH_BIN -i "$INPUT_FILE"
echo ""

# ==============================================================================
# Step 3: Measure topological charge
# ==============================================================================
echo -e "${GREEN}[Step 3/4] Measuring topological charge...${NC}"
./build/bin/$TOPCHARGE_BIN -i "$INPUT_FILE"
echo ""

# ==============================================================================
# Step 4: Analysis (if Python available)
# ==============================================================================
echo -e "${GREEN}[Step 4/4] Running analysis...${NC}"
if command -v python3 &> /dev/null; then
    if [ -n "$ANALYSIS_SCRIPT" ] && [ -f "$ANALYSIS_SCRIPT" ]; then
        echo "  Running full analysis..."
        python3 "$ANALYSIS_SCRIPT" || echo -e "${YELLOW}Analysis script failed (non-fatal)${NC}"
    else
        echo -e "${YELLOW}Analysis script not found: $ANALYSIS_SCRIPT${NC}"
    fi
else
    echo -e "${YELLOW}Python3 not found, skipping analysis${NC}"
fi

echo ""
echo -e "${GREEN}================================================${NC}"
echo -e "${GREEN}${GAUGE_GROUP^^} Simulation complete!${NC}"
echo -e "${GREEN}================================================${NC}"
echo ""
echo "Output files:"
echo "  - Plaquette history:  $PLAQ_FILE"
if [ "$GAUGE_GROUP" = "su2" ]; then
    echo "  - Configurations:     output/configs/"
else
    echo "  - Configurations:     output/configs_su3/"
fi
echo "  - Topological charge: $TOPCHARGE_FILE"
echo "  - Figures:            output/figures/"
