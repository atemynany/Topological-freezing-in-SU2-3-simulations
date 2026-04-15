#!/bin/bash
#SBATCH --account=agmisc
#SBATCH --job-name=lqcd_scan
#SBATCH --partition=fuchs
#SBATCH --qos=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=0
#SBATCH --time=240:00:00
#SBATCH --output=logs/scan_%j.out
#SBATCH --error=logs/scan_%j.err

cd /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations
mkdir -p logs

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

bash run_heatbath_scan.sh --su2 --skip-build --jobs 20
# Discovers all input/setup_*_su2.txt files automatically