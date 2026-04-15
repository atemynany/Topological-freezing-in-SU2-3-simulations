#!/bin/bash
#SBATCH --account=mesonqcd
#SBATCH --job-name=lqcd_scan
#SBATCH --partition=general2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=00-01:00:00
#SBATCH --output=/home/mesonqcd/barros/SU2_calc/Topological-freezing-in-SU2-3-simulations/logs/scan_%j.out
#SBATCH --error=/home/mesonqcd/barros/SU2_calc/Topological-freezing-in-SU2-3-simulations/logs/scan_%j.err

# Purge the environment to prevent dynamic linking conflicts
module purge

# Set working directory
cd /home/mesonqcd/barros/SU2_calc/Topological-freezing-in-SU2-3-simulations
mkdir -p logs

# Explicitly map OpenMP threads to the requested physical cores
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Execute the runtime script
bash run_heatbath_scan.sh --su2 --skip-build --jobs 1