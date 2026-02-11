#!/bin/bash
# ==============================================================================
# Step 2: Measure topological charge (with smearing)
# ==============================================================================
# Usage: ./scripts/step2_measure_topcharge.sh [input_file]
# ==============================================================================

set -e

INPUT_FILE="${1:-input/simulation.txt}"

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Error: Input file not found: $INPUT_FILE${NC}"
    exit 1
fi

echo -e "${GREEN}=== Step 2: Measure Topological Charge ===${NC}"
echo "Input file: $INPUT_FILE"
echo ""

# Build if needed
if [ ! -f "build/bin/meas_topcharge" ]; then
    echo "Building meas_topcharge..."
    mkdir -p build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null
    make -j$(nproc) meas_topcharge
    cd ..
fi

./build/bin/meas_topcharge -i "$INPUT_FILE"

echo ""
echo -e "${GREEN}Done. Results saved to output/topcharge.dat${NC}"
