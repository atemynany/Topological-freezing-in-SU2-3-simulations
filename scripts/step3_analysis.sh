#!/bin/bash
# ==============================================================================
# Step 3: Run analysis and generate plots
# ==============================================================================
# Usage: ./scripts/step3_analysis.sh
# ==============================================================================

set -e

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${GREEN}=== Step 3: Analysis ===${NC}"
echo ""

# Check data files exist
if [ ! -f "output/plaquette.dat" ]; then
    echo -e "${RED}Error: output/plaquette.dat not found. Run step1 first.${NC}"
    exit 1
fi

if [ ! -f "output/topcharge.dat" ]; then
    echo -e "${RED}Error: output/topcharge.dat not found. Run step2 first.${NC}"
    exit 1
fi

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || true
conda activate master_thesis 2>/dev/null || true

if [ -f "analysis/run_analysis.py" ]; then
    cd analysis
    python3 run_analysis.py
    cd ..
else
    echo -e "${RED}Error: analysis/run_analysis.py not found${NC}"
    exit 1
fi

echo ""
echo -e "${GREEN}Done. Figures saved to output/figures/${NC}"
