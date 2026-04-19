#!/bin/bash
#SBATCH --account=mesonqcd
#SBATCH --job-name=lqcd_scan
#SBATCH --partition=general2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=11-00:00:00
#SBATCH --output=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs/scan_%j.out
#SBATCH --error=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs/scan_%j.err

# Purge the environment to prevent dynamic linking conflicts
module purge

# Output/results go to work; binary and script run from home
export HEATBATH_RESULTS_DIR=/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/data/results
mkdir -p /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/logs

# Set working directory to home (binary lives here)
cd /home/mesonqcd/barros/SU_Calc/Topological-freezing-in-SU2-3-simulations

# Explicitly map OpenMP threads to the requested physical cores
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Execute the runtime script
bash run_heatbath_scan.sh --su2 --skip-build --jobs 1