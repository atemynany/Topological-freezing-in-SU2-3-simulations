#!/bin/bash
# ==============================================================================
# FUCHS cluster job script — Topcharge (single run, one beta)
# ==============================================================================
#SBATCH --account=agmisc
#SBATCH --job-name=lqcd_topcharge
#SBATCH --partition=fuchs
#SBATCH --qos=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=512
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --mail-type=FAIL
#SBATCH --output=logs/topcharge_%j.out
#SBATCH --error=logs/topcharge_%j.err

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
cd "$PROJECT_DIR"
mkdir -p logs

echo "======================================"
echo "Job ID:   $SLURM_JOB_ID"
echo "Node:     $SLURMD_NODENAME"
echo "Started:  $(date)"
echo "======================================"

# Build
bash scripts/build.sh release --march=avx

# Run
bash run_topcharge_scan.sh --su2 --skip-build

echo "======================================"
echo "Finished: $(date)"
echo "======================================"
