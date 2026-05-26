#!/bin/bash
#SBATCH --account=mesonqcd
#SBATCH --job-name=b315_topcharge
#SBATCH --partition=general2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=02-00:00:00
#SBATCH --output=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs/b315_topcharge_%j.out
#SBATCH --error=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs/b315_topcharge_%j.err

set -euo pipefail

module purge

PROJECT_DIR=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations
cd "$PROJECT_DIR"
mkdir -p logs

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

if [ ! -x build/bin/meas_topcharge ]; then
    mkdir -p build
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --target meas_topcharge --parallel "$SLURM_CPUS_PER_TASK"
fi

bash run_topcharge_scan.sh --su2 --skip-build --jobs 1
