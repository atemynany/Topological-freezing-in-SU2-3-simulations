#!/bin/bash
# ==============================================================================
# Step 1: Generate gauge configurations (MC heatbath)
# ==============================================================================
# Usage: ./scripts/step1_generate_configs.sh [input_file]
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

echo -e "${GREEN}=== Step 1: Generate Gauge Configurations ===${NC}"
echo "Input file: $INPUT_FILE"
echo ""

# Build if needed
if [ ! -f "build/bin/mc_heatbath" ]; then
    echo "Building mc_heatbath..."
    mkdir -p build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null
    make -j$(nproc) mc_heatbath
    cd ..
fi

./build/bin/mc_heatbath -i "$INPUT_FILE"

echo ""
echo -e "${GREEN}Done. Configurations saved to output/configs/${NC}"
