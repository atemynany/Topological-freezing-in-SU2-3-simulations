#!/bin/bash
# ==============================================================================
# FUCHS cluster job script — Heatbath (gauge config generation)
# Goethe University Frankfurt, CSC FUCHS cluster
#
# Usage:
#   sbatch cluster/fuchs_heatbath.sh [--su2|--su3] [--beta <val>]
#
# Build first on the login node:
#   ./scripts/build.sh release
# ==============================================================================

#SBATCH --job-name=lqcd_heatbath
#SBATCH --partition=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=512
#SBATCH --time=48:00:00
#SBATCH --no-requeue
#SBATCH --mail-type=FAIL
#SBATCH --output=logs/heatbath_%j.out
#SBATCH --error=logs/heatbath_%j.err
# If you have a Goethe-HLR account uncomment and fill in:
##SBATCH --account=<your_FUCHS_group>

# ==============================================================================
# Configuration — edit these as needed
# ==============================================================================
GAUGE_GROUP="su2"       # su2 or su3
BETA=""                 # leave empty to use value from params file
PARAMS_FILE=""          # leave empty to use default (input/base_params_su2.txt)

# Parse any arguments passed via sbatch --export or directly
while [[ $# -gt 0 ]]; do
    case $1 in
        --su2)  GAUGE_GROUP="su2"; shift ;;
        --su3)  GAUGE_GROUP="su3"; shift ;;
        --beta) BETA="$2"; shift 2 ;;
        --params-file) PARAMS_FILE="$2"; shift 2 ;;
        *) shift ;;
    esac
done

# ==============================================================================
# Environment
# ==============================================================================
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

mkdir -p logs

# SU(3) uses OpenMP — pin threads to the cores on this node
if [ "$GAUGE_GROUP" = "su3" ]; then
    export OMP_NUM_THREADS=20
    export OMP_PROC_BIND=close
fi

echo "======================================"
echo "FUCHS heatbath job"
echo "Job ID:       $SLURM_JOB_ID"
echo "Node:         $SLURMD_NODENAME"
echo "Gauge group:  $GAUGE_GROUP"
echo "Started:      $(date)"
echo "======================================"

# ==============================================================================
# Run — delegate to the existing heatbath script with --skip-build
# ==============================================================================
ARGS="--${GAUGE_GROUP} --skip-build"
[ -n "$BETA" ]        && ARGS="$ARGS --beta $BETA"
[ -n "$PARAMS_FILE" ] && ARGS="$ARGS --params-file $PARAMS_FILE"

bash scripts/run_heatbath.sh $ARGS

echo "======================================"
echo "Finished: $(date)"
echo "======================================"
