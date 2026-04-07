#!/bin/bash
#SBATCH --account=agmisc
#SBATCH --job-name=lqcd_topcharge
#SBATCH --partition=fuchs
#SBATCH --qos=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --output=logs/topcharge_%j.out
#SBATCH --error=logs/topcharge_%j.err

cd /work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations
mkdir -p logs

bash run_topcharge_scan.sh --su2 --skip-build --jobs 20
# Scans all data/results/*_su2 directories automatically