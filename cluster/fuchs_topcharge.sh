#!/bin/bash
#SBATCH --account=agmisc
#SBATCH --job-name=lqcd_topcharge
#SBATCH --partition=fuchs
#SBATCH --qos=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=12 
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=48:00:00
#SBATCH --output=logs/topcharge_%j.out
#SBATCH --error=logs/topcharge_%j.err

cd /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations
mkdir -p logs

export OMP_NUM_THREADS=1

# SU(2) — exclude large lattices still running heatbath (T65, T81)
bash run_topcharge_scan.sh --su2 --skip-build --jobs 12 --exclude "T65\|T81"

# Uncomment for SU(3):
# bash run_topcharge_scan.sh --su3 --skip-build --jobs 12
