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
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --mail-type=FAIL
#SBATCH --output=logs/topcharge_%j.out
#SBATCH --error=logs/topcharge_%j.err

PROJECT_DIR="/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations"
cd "$PROJECT_DIR"
mkdir -p logs

echo "======================================"
echo "Job ID:   $SLURM_JOB_ID"
echo "Node:     $SLURMD_NODENAME"
echo "Started:  $(date)"
echo "======================================"

# Verify binary exists (build on login node first)
if [ ! -f "build/bin/meas_topcharge" ]; then
    echo "ERROR: binary not found — build on login node first:"
    echo "  cmake -B build -DCMAKE_BUILD_TYPE=Release -DMARCH=avx ."
    echo "  cmake --build build --parallel 2"
    exit 1
fi

# Run
bash run_topcharge_scan.sh --su2 --skip-build

echo "======================================"
echo "Finished: $(date)"
echo "======================================"
