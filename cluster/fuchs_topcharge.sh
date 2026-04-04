#!/bin/bash
# ==============================================================================
# FUCHS cluster job script — Topological charge scan (all beta values)
# Goethe University Frankfurt, CSC FUCHS cluster
#
# Measures topological charge on all configs from the heatbath scan.
# SU(2): runs all betas in parallel (one per CPU core — 20 simultaneous).
# SU(3): runs betas sequentially (each beta uses all 20 cores via OpenMP).
#
# Usage:
#   sbatch cluster/fuchs_topcharge.sh [--su2|--su3]
#
# Requires heatbath scan to have been run first:
#   sbatch cluster/fuchs_heatbath.sh [--su2|--su3]
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
# Configuration
# ==============================================================================
GAUGE_GROUP="su2"

while [[ $# -gt 0 ]]; do
    case $1 in
        --su2) GAUGE_GROUP="su2"; shift ;;
        --su3) GAUGE_GROUP="su3"; shift ;;
        *) shift ;;
    esac
done

# ==============================================================================
# Environment
# ==============================================================================
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

mkdir -p logs

echo "======================================"
echo "FUCHS topcharge scan"
echo "Job ID:       $SLURM_JOB_ID"
echo "Node:         $SLURMD_NODENAME"
echo "Gauge group:  $GAUGE_GROUP"
echo "CPUs:         ${SLURM_NTASKS:-20}"
echo "Started:      $(date)"
echo "======================================"

# ==============================================================================
# Run
# ==============================================================================
if [ "$GAUGE_GROUP" = "su3" ]; then
    export OMP_NUM_THREADS=${SLURM_NTASKS:-20}
    export OMP_PROC_BIND=close
    echo "SU(3) mode: sequential betas, OMP_NUM_THREADS=$OMP_NUM_THREADS"
    bash run_topcharge_scan.sh --su3 --skip-build
else
    JOBS=${SLURM_NTASKS:-20}
    echo "SU(2) mode: up to $JOBS betas in parallel"
    bash run_topcharge_scan.sh --su2 --skip-build --jobs "$JOBS"
fi

echo "======================================"
echo "Finished: $(date)"
echo "======================================"
