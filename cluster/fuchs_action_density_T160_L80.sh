#!/bin/bash
#SBATCH --account=agmisc
#SBATCH --job-name=ad_T160_L80
#SBATCH --partition=fuchs
#SBATCH --qos=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=0
#SBATCH --time=48:00:00
#SBATCH --output=logs/action_density_T160_L80_%j.out
#SBATCH --error=logs/action_density_T160_L80_%j.err

set -euo pipefail

cd /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations
mkdir -p logs

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-20}"
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

INPUT_FILE="input/heatbath_input/action_density_T160_L80_open_su2.txt"
BINARY="./build/bin/meas_action_density_su2"

if [ ! -x "$BINARY" ]; then
    echo "Missing $BINARY; build it first with:"
    echo "  cmake --build build --target meas_action_density_su2 -j4"
    exit 1
fi

echo "Starting bulk action-density recalculation"
echo "Input:   $INPUT_FILE"
echo "Threads: $OMP_NUM_THREADS"
date

srun --cpu-bind=cores "$BINARY" -i "$INPUT_FILE"

date
echo "Done"
