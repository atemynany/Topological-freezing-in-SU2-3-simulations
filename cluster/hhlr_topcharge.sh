#!/bin/bash
#SBATCH --account=mesonqcd
#SBATCH --job-name=lqcd_topcharge
#SBATCH --partition=general2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=00-12:00:00
#SBATCH --output=/home/mesonqcd/barros/SU2_calc/Topological-freezing-in-SU2-3-simulations/logs/topcharge_%j.out
#SBATCH --error=/home/mesonqcd/barros/SU2_calc/Topological-freezing-in-SU2-3-simulations/logs/topcharge_%j.err

# Purge the environment to prevent dynamic linking conflicts
module purge

# Set working directory
cd /home/mesonqcd/barros/SU2_calc/Topological-freezing-in-SU2-3-simulations
mkdir -p logs

# Topcharge measurement is single-threaded per run; parallelism comes from --jobs
export OMP_NUM_THREADS=1

bash run_topcharge_scan.sh --su2 --skip-build --jobs 40
