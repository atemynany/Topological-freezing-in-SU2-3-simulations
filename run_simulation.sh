#!/bin/bash
# ==============================================================================
# SU(2) Lattice Gauge Theory Simulation Runner
# ==============================================================================
# This script runs the full simulation pipeline:
# 1. Build the project
# 2. Generate gauge configurations (MC heatbath)
# 3. Measure topological charge
# 4. Run analysis and plot results
#
# Usage:
#   ./run_simulation.sh [input_file]
#
# Default input file: input/simulation.txt
# ==============================================================================

set -e  # Exit on error

# Default input file
INPUT_FILE="${1:-input/simulation.txt}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}================================================${NC}"
echo -e "${GREEN}SU(2) Lattice Gauge Theory Simulation${NC}"
echo -e "${GREEN}================================================${NC}"
echo ""

# Check input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Error: Input file not found: $INPUT_FILE${NC}"
    exit 1
fi

echo -e "${YELLOW}Using input file: $INPUT_FILE${NC}"
echo ""

# ==============================================================================
# Step 1: Build
# ==============================================================================
echo -e "${GREEN}[Step 1/4] Building project...${NC}"
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null
make -j$(nproc) 2>&1 | tail -5
cd ..
echo -e "${GREEN}Build complete.${NC}"
echo ""

# ==============================================================================
# Step 2: Generate configurations
# ==============================================================================
echo -e "${GREEN}[Step 2/4] Generating gauge configurations...${NC}"
./build/bin/mc_heatbath -i "$INPUT_FILE"
echo ""

# ==============================================================================
# Step 3: Measure topological charge
# ==============================================================================
echo -e "${GREEN}[Step 3/4] Measuring topological charge...${NC}"
./build/bin/meas_topcharge -i "$INPUT_FILE"
echo ""

# ==============================================================================
# Step 4: Analysis (if Python available)
# ==============================================================================
echo -e "${GREEN}[Step 4/4] Running analysis...${NC}"
if command -v python3 &> /dev/null; then
    if [ -f "analysis/run_analysis.py" ]; then
        python3 analysis/run_analysis.py || echo -e "${YELLOW}Analysis script failed (non-fatal)${NC}"
    else
        echo -e "${YELLOW}Analysis script not found${NC}"
    fi
else
    echo -e "${YELLOW}Python3 not found, skipping analysis${NC}"
fi

echo ""
echo -e "${GREEN}================================================${NC}"
echo -e "${GREEN}Simulation complete!${NC}"
echo -e "${GREEN}================================================${NC}"
echo ""
echo "Output files:"
echo "  - Plaquette history: output/plaquette.dat"
echo "  - Configurations:    output/configs/"
echo "  - Topological charge: output/topcharge_*.dat"
echo "  - Figures:           output/results/figures/"
