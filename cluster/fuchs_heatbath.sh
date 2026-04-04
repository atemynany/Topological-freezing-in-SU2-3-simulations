#!/bin/bash
# ==============================================================================
# FUCHS cluster job script — Heatbath scan (all beta values)
# Goethe University Frankfurt, CSC FUCHS cluster
#
# Runs the full beta scan from input/beta_scan_<gauge>.txt.
# SU(2): runs all betas in parallel (one per CPU core — 20 simultaneous).
# SU(3): runs betas sequentially (each beta uses all 20 cores via OpenMP).
#
# Usage:
#   sbatch cluster/fuchs_heatbath.sh [--su2|--su3]
#
# Build first on the login node:
#   ./scripts/build.sh release
# ==============================================================================

#SBATCH --account=agmisc
#SBATCH --job-name=lqcd_heatbath
#SBATCH --partition=fuchs
#SBATCH --qos=fuchs
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
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$PROJECT_DIR"

mkdir -p logs

echo "======================================"
echo "FUCHS heatbath scan"
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
    # SU(3): each beta uses all cores via OpenMP — run betas sequentially
    export OMP_NUM_THREADS=${SLURM_NTASKS:-20}
    export OMP_PROC_BIND=close
    echo "SU(3) mode: sequential betas, OMP_NUM_THREADS=$OMP_NUM_THREADS"
    bash run_heatbath_scan.sh --su3 --skip-build
else
    # SU(2): each beta is single-threaded — run all betas in parallel
    JOBS=${SLURM_NTASKS:-20}
    echo "SU(2) mode: up to $JOBS betas in parallel"
    bash run_heatbath_scan.sh --su2 --skip-build --jobs "$JOBS"
fi

echo "======================================"
echo "Finished: $(date)"
echo "======================================"
