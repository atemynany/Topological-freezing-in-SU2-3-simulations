#!/bin/bash
# ==============================================================================
# FUCHS cluster job script — Topological charge measurement
# Goethe University Frankfurt, CSC FUCHS cluster
#
# Usage:
#   sbatch cluster/fuchs_topcharge.sh <input_file> [--su3]
#
# Example:
#   sbatch cluster/fuchs_topcharge.sh data/results/T16_L16_b2.5_periodic_seed12345/input.txt
#
# Build first on the login node:
#   ./scripts/build.sh release
# ==============================================================================

#SBATCH --job-name=lqcd_topcharge
#SBATCH --partition=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=512
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --mail-type=FAIL
#SBATCH --output=logs/topcharge_%j.out
#SBATCH --error=logs/topcharge_%j.err
# If you have a Goethe-HLR account uncomment and fill in:
##SBATCH --account=<your_FUCHS_group>

# ==============================================================================
# Configuration — edit these as needed
# ==============================================================================
GAUGE_GROUP="su2"   # su2 or su3
INPUT_FILE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --su3) GAUGE_GROUP="su3"; shift ;;
        --su2) GAUGE_GROUP="su2"; shift ;;
        *)     INPUT_FILE="$1";   shift ;;
    esac
done

if [ -z "$INPUT_FILE" ]; then
    echo "Error: no input file specified."
    echo "Usage: sbatch cluster/fuchs_topcharge.sh <input_file> [--su3]"
    exit 1
fi

# ==============================================================================
# Environment
# ==============================================================================
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

mkdir -p logs

# SU(3) measurement can use OpenMP for smearing loops
if [ "$GAUGE_GROUP" = "su3" ]; then
    export OMP_NUM_THREADS=20
    export OMP_PROC_BIND=close
fi

echo "======================================"
echo "FUCHS topcharge job"
echo "Job ID:       $SLURM_JOB_ID"
echo "Node:         $SLURMD_NODENAME"
echo "Gauge group:  $GAUGE_GROUP"
echo "Input file:   $INPUT_FILE"
echo "Started:      $(date)"
echo "======================================"

# ==============================================================================
# Run — delegate to the existing topcharge script
# ==============================================================================
bash scripts/run_topcharge.sh --${GAUGE_GROUP} "$INPUT_FILE"

echo "======================================"
echo "Finished: $(date)"
echo "======================================"
