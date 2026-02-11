#!/bin/bash
# ==============================================================================
# run_instanton_Q.sh
# 
# Script to build and run instanton topological charge computation,
# then plot the results.
#
# Usage:
#   ./scripts/run_instanton_Q.sh                           # Use defaults
#   ./scripts/run_instanton_Q.sh -i input/instanton.txt    # Use input file
#   L=20 RHO=4.0 ./scripts/run_instanton_Q.sh              # Override with env vars
# ==============================================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# ------------------------------------------------------------------------------
# Parse command line arguments
# ------------------------------------------------------------------------------
INPUT_FILE=""
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [-i input_file]"
            echo ""
            echo "Options:"
            echo "  -i, --input FILE    Read parameters from input file"
            echo "  -h, --help          Show this help message"
            echo ""
            echo "Environment variables (override input file):"
            echo "  L, T, RHO, A, NSMEAR, ALPHA, OUTPUT_DIR"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# ------------------------------------------------------------------------------
# Function to parse input file
# ------------------------------------------------------------------------------
parse_input_file() {
    local file="$1"
    if [ ! -f "$file" ]; then
        echo -e "${RED}Error: Input file not found: $file${NC}"
        exit 1
    fi
    
    while IFS= read -r line || [ -n "$line" ]; do
        # Skip comments and empty lines
        [[ "$line" =~ ^[[:space:]]*# ]] && continue
        [[ -z "${line// }" ]] && continue
        
        # Parse key-value pairs
        key=$(echo "$line" | awk '{print $1}')
        value=$(echo "$line" | awk '{print $2}')
        
        case "$key" in
            L)              FILE_L="$value" ;;
            T)              FILE_T="$value" ;;
            rho)            FILE_RHO="$value" ;;
            a)              FILE_A="$value" ;;
            smear_steps)    FILE_NSMEAR="$value" ;;
            smear_alpha)    FILE_ALPHA="$value" ;;
            output_dir)     FILE_OUTPUT_DIR="$value" ;;
            output_prefix)  FILE_OUTPUT_PREFIX="$value" ;;
        esac
    done < "$file"
}

# ------------------------------------------------------------------------------
# Load parameters: defaults -> input file -> environment variables
# ------------------------------------------------------------------------------

# Defaults
DEFAULT_L=16
DEFAULT_T=16
DEFAULT_RHO=3.0
DEFAULT_A=1.0
DEFAULT_NSMEAR=100
DEFAULT_ALPHA=0.5
DEFAULT_OUTPUT_DIR="output/instanton"
DEFAULT_OUTPUT_PREFIX="instanton_Q"

# Parse input file if provided
if [ -n "$INPUT_FILE" ]; then
    echo -e "${YELLOW}Reading parameters from: $INPUT_FILE${NC}"
    parse_input_file "$INPUT_FILE"
fi

# Final values: env var > input file > default
L=${L:-${FILE_L:-$DEFAULT_L}}
T=${T:-${FILE_T:-$DEFAULT_T}}
RHO=${RHO:-${FILE_RHO:-$DEFAULT_RHO}}
A=${A:-${FILE_A:-$DEFAULT_A}}
NSMEAR=${NSMEAR:-${FILE_NSMEAR:-$DEFAULT_NSMEAR}}
ALPHA=${ALPHA:-${FILE_ALPHA:-$DEFAULT_ALPHA}}
OUTPUT_DIR=${OUTPUT_DIR:-${FILE_OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}}
OUTPUT_PREFIX=${OUTPUT_PREFIX:-${FILE_OUTPUT_PREFIX:-$DEFAULT_OUTPUT_PREFIX}}

echo -e "${GREEN}=== Instanton Topological Charge Analysis ===${NC}"
echo ""

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="${PROJECT_DIR}/build"

# Step 1: Build the project
echo -e "${YELLOW}[1/4] Building project...${NC}"
cd "$PROJECT_DIR"

if [ ! -d "$BUILD_DIR" ]; then
    mkdir -p "$BUILD_DIR"
fi

cd "$BUILD_DIR"
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc) compute_instanton_Q

# Step 2: Create output directory
echo -e "${YELLOW}[2/4] Setting up output directory...${NC}"
mkdir -p "$PROJECT_DIR/$OUTPUT_DIR"
cd "$PROJECT_DIR/$OUTPUT_DIR"

# Step 3: Run the computation (just Q at step 0, no smearing loop)
echo -e "${YELLOW}[3/3] Computing instanton topological charge...${NC}"
echo "Parameters:"
echo "  Lattice: ${L}^3 x ${T}"
echo "  Instanton size rho: ${RHO}"
echo "  Lattice spacing a: ${A}"
echo ""

"${BUILD_DIR}/bin/compute_instanton_Q" \
    -L "$L" \
    -T "$T" \
    -rho "$RHO" \
    -a "$A" \
    -nsmear 0 \
    -alpha "$ALPHA" \
    -o "$OUTPUT_PREFIX"

echo ""
echo -e "${GREEN}=== Done! ===${NC}"
echo "Output files:"
echo "  Data: $OUTPUT_DIR/${OUTPUT_PREFIX}_instanton.dat"
echo "  Data: $OUTPUT_DIR/${OUTPUT_PREFIX}_anti_instanton.dat"

# Commented out: Q vs smearing and plotting
# To enable, set NSMEAR > 0 and uncomment below:
#
# "${BUILD_DIR}/bin/compute_instanton_Q" \
#     -L "$L" -T "$T" -rho "$RHO" -a "$A" \
#     -nsmear "$NSMEAR" -alpha "$ALPHA" -o "$OUTPUT_PREFIX"
#
# source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
# conda activate master_thesis 2>/dev/null || true
# PLOT_SCRIPT="${PROJECT_DIR}/../lattice_su2_Vilija_ansatz/Instanton/plotting_Q_instanton_smearsteps.py"
# python3 "$PLOT_SCRIPT" "$OUTPUT_DIR/${OUTPUT_PREFIX}_plot.png" \
#     --inst-file "$OUTPUT_DIR/${OUTPUT_PREFIX}_instanton.dat" \
#     --anti-file "$OUTPUT_DIR/${OUTPUT_PREFIX}_anti_instanton.dat"
